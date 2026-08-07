# Code review — D13 (`d2a1c3a`), the `path_identity_key` / `collision_key` split

Scope: commit `d2a1c3a` only, on `fix/383-output-collision-preflight`. Files `src/io.rs`,
`src/cli.rs`. All behaviour below was executed against `target/release/trim_galore` (current
with the change). Scratch: `/private/tmp/claude-501/…/scratchpad/d13/`. No repo file other
than this report was touched.

## What the commit does

`io::collision_key` is factored into a private `lexical_normalise(&Path) -> PathBuf`
(absolutise, drop `.`, fold `..`) plus two public string keys over it:

- `path_identity_key` — case **preserved** → the three `Cli::validate` input-identity checks
  (`src/cli.rs:537`, `:549-550`, `:635`/`:638`).
- `collision_key` — `norm_path()`-folded → `io::preflight_output_collisions` (`src/io.rs:97`,
  `:102`) and the `--passthrough`-aliases-R1/R2 check (`src/cli.rs:939-941`).

## Verification setup

The machine's `/` is case-insensitive APFS. To test the Linux direction properly I built a
genuinely case-sensitive volume inside scratch and ran both sides:

```
hdiutil create -size 60m -fs "Case-sensitive APFS" -volname CSTEST cs.dmg
hdiutil attach cs.dmg -mountpoint "$PWD/csmnt" -nobrowse     # needed sandbox bypass
```
Confirmed case-sensitive (`Aa` and `aA` coexist as two files); the APFS side confirmed
case-insensitive (`r1.fq` and `R1.fq` share inode 194103634).

| # | FS | invocation | result |
|---|----|-----------|--------|
| B | APFS | `--paired ./a_R1.fq a_R1.fq` | refused, "appear to be the same file" ✓ |
| C | APFS | `a_R1.fq A_R1.fq` (SE) | refused, output-collision/APFS message ✓ |
| D | APFS | `--paired Sample_R1 Sample_R2 SAMPLE_R1 SAMPLE_R2` | output-collision/APFS message ✓ |
| E | APFS | `--paired S_R1 S_R2 S_R1 S_R2` | "duplicate of pair 1" ✓ (fires before pre-flight) |
| G | **case-sens** | the #216 guard, four genuinely distinct files | output-collision/APFS message ✓ |
| H | **case-sens** | `--paired r1.fq R1.fq`, two real files | runs; 2 outputs + 4 distinct reports ✓ |
| A | APFS | `--paired r1.fq R1.fq` (one physical file) | **runs** — see Q1 |
| I | **case-sens** | `--paired r1.fq r2.fq --passthrough R1.fq` | **refused** — see Q2 |

G is the load-bearing one: on a real case-sensitive filesystem the commit produces exactly the
message the CI guard asserts, and E shows the `Cli::validate` duplicate-pair check does pre-empt
the pre-flight when it fires — so the D8 failure mode is understood and closed.

Aside: the commit message's "not reproducible on macOS" is half right. The four-distinct-files
*premise* needs a case-sensitive FS, but which message fires is pure path math — case D above
reproduces the diagnostic swap on APFS. A local test was therefore available; see finding 3.

## Q1 — is case-sensitive input identity right in both directions?

**Yes, ship it.** Direction that matters is verified (G, H): on a case-sensitive filesystem two
files differing only in case are two files, they pair correctly, and the #216 guard reaches the
output pre-flight. Case-folding input identity rejects valid input, which is a worse failure
than the residual below.

The APFS false-negative is real and I reproduced it (A). `--paired r1.fq R1.fq` on APFS names
one physical file twice: `path_identity_key` sees two paths, and the outputs (`r1_val_1.fq` /
`R1_val_2.fq`) fold to different keys because the `_1`/`_2` suffix differs, so nothing catches
it. The run self-pairs the file and exits 0. Two consequences:

1. Nonsense-but-harmless output: R1 and R2 are the same reads.
2. **One trimming report is silently overwritten.** The run announces
   `r1.fq_trimming_report.txt` and `R1.fq_trimming_report.txt`; only one file exists afterwards
   and it contains `Input filename: R1.fq`. The R1 report is gone — the #383 harm class.

Two things keep this from being a blocker. First, it is **not a regression**: `dev` compared
inputs with raw `PathBuf ==`, also case-sensitive, so `dev` behaves identically. The parent
commit `080c1d0` covered this case only as a side effect of the overreach that broke CI.
Second, the sharp SE case *is* caught (C) — `a.fq A.fq` collides on the trimmed-output key with
the APFS message, and SE is where #383's data loss lived.

The reachable-harm set is also narrower than it looks: multi-pair case-only variants *are*
caught, because a `_val_1` path folds and collides (D). Only R1==R2-within-one-pair escapes,
precisely because `_1`/`_2` differ.

**Recommended follow-up (not this commit's fault):** the paired pre-flight's candidate list is
`vec![o1, o2]` plus unpaired/passthrough (`src/main.rs:740`) — it omits the report paths, while
the SE path includes them via `planned_secondary_outputs` (`src/main.rs:90-93`). Adding the two
report paths per pair converts the Q1 false-negative into a loud pre-flight error on APFS at
zero cost on Linux (H shows those four report names stay distinct on a case-sensitive FS, so no
false positive is possible). The principled fix for the whole class — filesystem identity via
`dev`+`ino` instead of string keys, which would also retire the symlink residual A8 — is bigger
than this PR should carry.

## Q2 — is every call site on the right key?

`grep` over `src/` finds no remaining raw path-equality comparison and no other consumer of
either key. Audit:

| site | question it asks | key | verdict |
|---|---|---|---|
| `cli.rs:537` R1≠R2 in pair | same file? | identity | correct |
| `cli.rs:549-550` duplicate pair | same file? | identity | correct |
| `cli.rs:635`/`:638` duplicate SE input | same file? | identity | correct |
| `io.rs:97` pre-flight input map | would an output land on an input? | collision | correct |
| `io.rs:102` pre-flight planned | same output file? | collision | correct |
| `cli.rs:939-941` `--passthrough` alias | **same file?** | collision | **finding 2** |

`guarded_inputs` (`src/main.rs:64`) feeds the pre-flight's input side and correctly inherits
the folded key — "would an output land on this input" should fold on APFS.

The `--passthrough` check is the one mismatch. Its own comment (`src/cli.rs:935-937`) says it
catches "case-only INPUT aliases" — an identity question wearing a collision key. Test I shows
the consequence: on a case-sensitive volume, with `R1.fq` genuinely distinct from `r1.fq`, the
run is refused with a message asserting `R1.fq` "aliases an input file", false on that
filesystem. That is verbatim the failure mode the commit fixes for the other three checks. It
is pre-existing and deliberate, it fails loudly, and the input shape is exotic, so I would not
block — but the commit's own rationale argues against it.

## The `..`-folding loop

Executed a verbatim copy of `lexical_normalise` over edge cases (`scratchpad/d13/norm.rs`):

```
"/.." -> "/"    "/../.." -> "/"    "/../x" -> "/x"    "//a" -> "/a"    "/a/" -> "/a"
"/a/b/../.." -> "/"    "/a/b/../../.." -> "/"    "/a/b/../c/./d/../e" -> "/a/c/e"
```

POSIX `/..` == `/` is right, and `..` cannot walk above root or leak a stray `..` into the key.
Two dead-code notes, both harmless: the `_ => stack.push(c)` unfoldable-`..` arm is unreachable
on POSIX, because `std::path::absolute` only fails on the empty path (verified:
`absolute("")` → `InvalidInput`, `absolute("..")` → `Ok`) and an empty path has no components;
and the `Component::CurDir => {}` arm is redundant, since `absolute()` plus `Components`'
own skipping of interior `.` already removes it. Windows `Prefix` handling looks right for
drive-rooted paths, and verbatim (`\\?\`) paths — where folding `..` would change meaning
because Windows treats it literally — are moot: the CI matrix is `ubuntu-latest` and
`macos-latest` only.

## Can the two new tests fail?

Mutation-checked in the standalone copy rather than by editing the repo.

- `identity_key_is_case_sensitive_but_collision_key_is_not` — **yes, load-bearing.** Re-folding
  case in `path_identity_key` (i.e. re-doing D8) flips the `assert_ne` to false. Both halves
  bite.
- `both_keys_normalise_spelling` — **yes, but narrowly.** Removing the `ParentDir` fold makes
  the `a/../x.fq` case fail. It does **not** catch removal of the `absolute()` call
  (`./x.fq` and `a/../x.fq` still normalise to `x.fq` without it), so the absolutisation that
  the original H1 finding was actually about is unpinned. Comparing `x.fq` against
  `current_dir().join("x.fq")` would close that.

## Findings

1. **APFS false-negative, low-medium, pre-existing, not blocking.** `--paired r1.fq R1.fq` on a
   case-insensitive FS self-pairs one file and silently overwrites one trimming report
   (reproduced). Same on `dev`. Cheapest fix: add the per-pair report paths to the paired
   pre-flight candidate list at `src/main.rs:740`.
2. **`--passthrough` alias check is an identity question on a collision key, low.**
   `src/cli.rs:939-941` false-positives on a case-sensitive filesystem (reproduced, test I).
   Pre-existing and loud. Either fold it into `path_identity_key` and accept the APFS gap
   symmetrically, or note in the comment why passthrough is treated differently from the other
   three identity checks.
3. **The regression that caused this commit has no unit-level guard, low, cheap to fix.** The
   two new tests pin the *keys*; what broke was the *wiring* in `cli.rs`. A test in the
   existing `Cli::parse_from(…).validate()` style — no files on disk needed, so it runs
   identically on macOS and Linux — would have failed on D8:

   ```rust
   let cli = Cli::parse_from(["trim_galore", "--paired",
       "Sample_R1.fastq.gz", "Sample_R2.fastq.gz", "SAMPLE_R1.fastq.gz", "SAMPLE_R2.fastq.gz"]);
   assert!(cli.validate().is_ok(), "case-only variants are distinct inputs");
   ```
   No existing test uses those four names against `validate()`; only the CI validation job
   covers it.
4. **Nits:** `both_keys_normalise_spelling` does not pin absolutisation (above); and
   `norm_path`'s doc block (`src/io.rs:33-42`) is still the only place describing the key
   family yet never mentions `path_identity_key` or the identity/collision distinction — one
   line would make the split discoverable from the function a reader hits first.

None of the five changes behaviour that this commit got wrong in the direction it set out to
fix, and none is a regression against `dev`.

MERGE
