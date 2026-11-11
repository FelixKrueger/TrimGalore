# CODE REVIEW — Reviewer B — `9a738a3` (`fix/440-library-refusals`)

**Target:** `git diff 872c4ae 9a738a3` — 9 files, +382/−25
**Plan:** `plans/08162026_pr441-library-presets/PLAN.md` v2 (incl. §12 implementation notes)
**Worktree:** `/private/tmp/claude-501/pr441-B` @ `9a738a333f1c2c75219b4f15ce2ebff418447495` (verified)

```
$ git rev-parse HEAD
9a738a333f1c2c75219b4f15ce2ebff418447495
```

---

# VERDICT: **ready with changes**

Both behavioural changes are correct and complete. I re-derived every precedence guarantee the plan
claims by running the binary, confirmed the "does not honour the clip flags" wording against
`specialty.rs` myself, and confirmed the banner and the text report cannot disagree on **any** of the
four paired shapes plus both single-end shapes. Six of the eleven mutation-matrix rows reproduced
exactly — including the honest version of V3, which §12 admits (D3) was only run as an equivalent.

Four changes I would want before push, none of them a defect in shipped behaviour:

| # | Fix | Size |
|---|---|---|
| **F-1** | Assert the **top-level** JSON clip keys in `test_json_report_records_all_four_overrides`. A permutation of those four key labels survives all 764 tests — **proved** below. This is the T-1 class the commit set out to close, in the one place it didn't look. | 4 lines, in a test that already has the right fixture |
| **F-3** | Add the refusal clause to `--library`'s `--help` text (`src/cli.rs:138-142`). `--fastqc` and `--rename` both document theirs; the plan only considered `flags.md`. | one clause |
| **F-4** | `src/cli.rs:1071-1072` now says something false — the `--fastqc` block is no longer "Last of the fallible checks". | one word |
| **F-2** | `tests/integration_library.rs:147-152`'s comment states an invariant the tool does not hold (**disproved** below). | reword + narrow the assertion |

Plus a decision to record for **F-5** (the refusal covers the derived flag but not the four primitives
it expands into — measured as a silent no-op) and **F-7** (docs list 4 of the 5 modes that refuse).

---

# 1. Build, suite, lints — reproduced

```
$ CARGO_PROFILE_DEV_DEBUG=0 cargo test --no-fail-fast
$ grep -c "^test result" B_baseline.log                       → 17
$ awk over the "test result" lines                            → TOTAL PASSED: 764
$ grep "^test result" B_baseline.log | grep -v " 0 failed"     → (no output)
```

**764 passing across 17 binaries, 0 failures** — matches §12 exactly. `--no-fail-fast` was used, so the
lib-target-truncation hazard does not apply.

```
$ cargo fmt --all -- --check
FMT: CLEAN
```

Clippy was run in **debug** with `--all-targets -- -D warnings`, not `--release`: the volume had
3.9 GB free and a release build already failed with `os error 28` in this feature's history. It was
**proved able to fail** before being believed — a canary (`if v.len() == 0`) appended to `src/library.rs`:

```
$ cargo clippy --all-targets -- -D warnings
CANARY EXIT=101
error: length comparison to zero
error: could not compile `trim-galore` (lib) due to 1 previous error
```

Canary removed, `git diff --stat src/` empty, re-run: exit 0, no warnings.

```
    Compiling trim-galore v2.3.0 (/private/tmp/claude-501/pr441-B)
    Finished `dev` profile [unoptimized] target(s) in 2.41s
```

**Caveat:** release-only lints are not covered by this run.

---

# 2. Correctness of the refusal — correct

**Siting.** `src/cli.rs:1112-1134`. `bail!` line numbers in `validate()` (which spans `605-1165`) run
`…, 1103, 1115, 1120, 1125, 1130` and then nothing before `Ok(())` at `1164`; everything from `1136`
to `1162` is `eprintln!`. The block **is** the last fallible statement in `validate()`, exactly as the
plan claims, so every precedence guarantee follows structurally rather than by luck.

**Every precedence claim verified by running the binary** (a function taking `"$@"`, so the zsh
word-splitting hazard does not apply):

```
### --hardtrim5 0 --library emseq          exit=1  Error: --hardtrim5 must be between 1 and 999, got 0
### --hardtrim5 30 --hardtrim3 30 --library exit=1  Error: --hardtrim5 and --hardtrim3 cannot be combined…
### --clock --library pbat (1 file)        exit=1  Error: --clock mode requires an even number of input files…
### --implicon --library pbat (1 file)     exit=1  Error: --implicon mode requires an even number of input files…
### --hardtrim5 30 --library emseq /nonexistent.fq.gz
                                          exit=1  Error: Input file not found: /nonexistent.fq.gz
### --clump_only --hardtrim5 30 --library  exit=1  Error: --clump_only does not clip; --library is not compatible
### --hardtrim5 30 --library NOSUCHKIT     exit=2  error: invalid value 'NOSUCHKIT' for '--library <PRESET>'
```

All four modes refuse, exit 1, and **write nothing** (V4 reproduced):

```
### --hardtrim5 30 --library emseq -o .   exit=1  Error: --hardtrim5 does not honour the clip flags; --library is not compatible   files:
### --hardtrim3 30 --library emseq -o .   exit=1  Error: --hardtrim3 does not honour the clip flags; --library is not compatible   files:
### --clock --paired --library pbat -o .  exit=1  Error: --clock does not honour the clip flags; --library is not compatible       files:
### --implicon --paired --library pbat    exit=1  Error: --implicon does not honour the clip flags; --library is not compatible    files:
```

**Can any input reach a specialty arm with a preset still set? No.** `main.rs:237` is `cli.validate()?`;
the earliest specialty dispatch is `main.rs:468`. Additional shapes checked beyond the plan's list:

- `--implicon=8 --library pbat` (equals form) → refused.
- Multi-pair (2 distinct pairs, 4 files) `--clock` / `--implicon` / `--hardtrim5` + `--library` → all
  three refused.
- `--fastqc --hardtrim5 --library` → reports the #421 message. Sensible, and it is the direct
  consequence of siting our block after #421's. **This is what falsifies the #421 comment (F-4).**

---

# 3. The wording — accurate, and I verified the claim independently

The plan's factual basis holds in full:

```
$ grep -n "clip_r1\|clip_r2\|three_prime_clip\|effective_clips\|clip_5prime\|clip_3prime" src/specialty.rs
107:            record.clip_5prime(start);
219:            record.clip_5prime(start);
```

Two hits, both `clip_5prime`; **zero** references to `clip_r1`/`clip_r2`/`three_prime_clip_*`/
`effective_clips`. `--clock` and `--implicon` clip by direct slicing rather than through
`clip_5prime`, which the grep above would have missed — so I checked those too:

```
314:                let r1_trimmed_seq = r1.seq[13..].to_string();
319:                let r2_trimmed_seq = r2.seq[15..].to_string();
413:                let r2_trimmed_seq = r2.seq[umi_length..].to_string();
```

and the doc comments say so out loud (`specialty.rs:241`, `:352`):

```
/// as `:R1:<umi1>:R2:<umi2>:F1:<fix1>:F2:<fix2>`, then clips:
/// - R1 at position 13 (removes UMI + CAGT + A)
/// - R2 at position 15 (removes UMI + CAGT + A + 2bp end-repair bias)
…
/// UMI barcode (`:BARCODE`), then clips R2 by `umi_length` bases.
```

**So all four modes clip, none reads `--clip_*`, and v1's "does not clip" would indeed have
contradicted the tool's own documentation inside a user-facing error.** The wording is correct.

**Consistency with the neighbouring families.** The shape `<subject> <incapability>; <flag> is not
compatible` matches the `--clump_only` family verbatim (`cli.rs:834` `--clump_only does not clip;
--clip_r1 is not compatible`, `:850` the same for `--library`), and the per-mode incapability idiom
matches the `-a2` block (`cli.rs:1053` `--clock performs no adapter trimming`). No inconsistency found.

O-2 (`--clump_only does not clip` vs `does not honour the clip flags` for the same flag) is a real
difference in wording between two families, but both sentences are true of their subject —
`--clump_only` genuinely performs no clipping, the specialty modes genuinely do. Agree with the plan:
no change.

---

# 4. The gate — correct at both call sites, and they cannot disagree

```
$ grep -rn "flag_summary" src/
src/library.rs:120:    pub fn flag_summary(&self, paired: bool) -> String {
src/library.rs:290:            r.flag_summary(true),          ← test
src/library.rs:299:            r.flag_summary(false),         ← test
src/report.rs:282:            library.flag_summary(config.paired)
src/main.rs:1121:            clips.flag_summary(cli.paired)
```

Exactly two non-test call sites, **both updated**. The plan's load-bearing detail holds.

**Is `paired` the right source of truth at each site?** Yes, and the two cannot diverge:

- Every `TrimConfig` that sets `library` also sets `paired` literally — `main.rs:1459/1471`
  (`run_single_file`, `false`), `:2249/2261` (`run_ubam_output_single`, `false`), `:2628/2640`
  (`write_paired_reports`, `true`). There is no fourth site.
- `write_paired_reports` is the **single funnel** for all four paired shapes — four call sites
  (`main.rs:1801`, `:2008`, `:2429`, `:2548`), one function, `paired: true`.
- Dispatch into all six report-writing paths branches on `cli.paired`, so `config.paired == cli.paired`
  at every one. `--paired` with one input is *required* to be a uBAM (`main.rs:340-346`), so there is
  no "effectively paired but `cli.paired == false`" shape for the banner to get wrong.

**Measured on every shape rather than argued.** Banner (stderr) vs text report:

| Shape | Banner | Report |
|---|---|---|
| SE FASTQ, `--library emseq --clip_R1 12` | `--clip_R1 12 --three_prime_clip_R1 10` | same |
| PE FASTQ two files, `--library accel` | `--clip_R1 10 --clip_R2 15 --three_prime_clip_R1 10 --three_prime_clip_R2 10` | same (both sides) |
| SE + `--output-format ubam`, `accel` | `--clip_R1 10 --three_prime_clip_R1 10` | same |
| PE two files + `--output-format ubam`, `accel` | all four | same (both sides) |
| **interleaved uBAM in**, PE, FASTQ out, `accel` | all four | same (both sides) |
| **interleaved uBAM in**, PE, uBAM out, `accel` | all four | same (both sides) |

`--clump_only` cannot reach the gate at all (`cli.rs:850` refuses `--library` first) — verified above.
Multi-pair runs go through the same `write_paired_reports`, so they are covered by row 2.

---

# 5. What the gate should NOT have touched — holds

**The JSON is untouched.** The diff contains **zero** changes to the JSON writer
(`src/report.rs:863-906`); the only non-test `report.rs` change is the text line at `:279-283`. No key
added, removed, renamed or reordered.

**A run with no `--library` cannot behave differently** — this is a reachability argument, not a
measured differential (a base binary would have needed another ~1.5 GB of `target/` on a volume with
3.9 GB free). The three source changes are: `flag_summary`, reached only under
`config.library.preset.is_some()` / `clips.preset.is_some()` (`report.rs:275-276`, `main.rs:1117`); the
`validate()` block, gated on `self.library.is_some()`; and tests. Confirmed by running `--rrbs --paired`
with no preset:

```
exit=0
occurrences of "Library preset|Clipping in force" in stderr:  0
occurrences of "Library preset|Clipping in force" in report:  0
library key present: True | value: None
top-level clip_r2: 2        ← the RRBS auto-clip still fires
```

**§8.3's reasoning holds.** `json["mode"]` exists and is derived from the same `config.paired` the gate
reads (`report.rs:808-813`, `"paired-end"` / `"single-end"`), so a machine consumer really can
distinguish the two cases; there is indeed no `paired` key. And the docs' claim to that effect is
accurate. The single-end run above emitted `mode= single-end` with all four values in the block.

**The override log lines are correctly left ungated** as a mechanism — but see **F-2** and **F-6**: the
*justification* the plan gives for leaving them ungated does not survive contact with the output.

---

# 6. Test quality — matrix spot-checks

Six rows re-run from scratch. **Every one reproduced exactly.** Each mutation was applied with a
`python3` assert on the match count (so a silent no-op mutation is impossible), the whole suite was run
with `--no-fail-fast`, and the file was restored and `cmp`-verified afterwards.

### V6a — only `main.rs:1121` reverted to `flag_summary(true)` ✅ exact

```
cargo exit=101
    explicit_clip_flag_overrides_the_preset_and_is_logged
test result: FAILED. 7 passed; 1 failed
```

**Exactly one test**, as §12 claims. The C-2 gap is genuinely closed — this is the single most
important row in the matrix and it is real.

### V8 — JSON block `clip_r2` / `three_prime_clip_r1` **key labels** swapped ✅ exact

```
cargo exit=101
    report::tests::test_json_report_records_all_four_overrides
test result: FAILED. 480 passed; 1 failed
```

T-1 closed for the block. (Nothing else catches it: the integration agreement loop runs on `pbat`
8/8/8/8, where a permutation is invisible — which is the seam in **F-1**.)

### V11 — `library` key omitted on no-preset runs ✅ exact

```
cargo exit=101
    a_run_without_a_preset_reports_no_library_block
    report::tests::test_json_report_library_is_null_without_a_preset
```

And the panic proves it is the **new** assertion doing the work, with the pre-existing `is_null()`
having passed on an absent key — so F-5a's "both stayed green before" is corroborated, not assumed:

```
panicked at src/report.rs:1553:9:            the library key must be present and null, not absent
panicked at tests/integration_library.rs:236:5: the library key must be present and null, not absent
```

### V13 — overrides comma guard removed ✅ exact

```
cargo exit=101
    preset_with_multiple_overrides_round_trips_through_json
    report::tests::test_json_report_records_all_four_overrides
    report::tests::test_json_report_records_library_block
```

3 red including the cardinality-1 test. **Keeping `test_json_report_records_library_block` was right**
(B's plan-review #5 vindicated; v1 would have deleted the only cardinality-1 coverage).

### V2 — refusal block disabled ✅ exact

```
cargo exit=101
    cli::tests::library_plus_clock_rejected
    cli::tests::library_plus_hardtrim3_rejected
    cli::tests::library_plus_hardtrim5_rejected
    cli::tests::library_plus_implicon_rejected
```

### V7 — R2 dropped unconditionally (`if read_2` without `&& !paired`) ✅ exact

```
cargo exit=101
    library::tests::flag_summary_lists_all_four_for_paired
    paired_report_lists_all_four_values
    report::tests::test_text_report_records_preset_name_and_expanded_values
```

### V3 — run **as the plan wrote it**, not as §12's D3 substitute ✅

§12's D3 records that V3 was run as an equivalent mutation (suppress the range `bail!`) rather than by
relocating the block. I ran the version the plan specifies: physically moved the whole refusal block
above the #386 check, so it now precedes the range check at what was `cli.rs:1016`.

```
$ grep -n "does not honour the clip flags" src/cli.rs   → 1011, 1016  (was 1116, 1121)
$ grep -n "#386 — dispatch runs" src/cli.rs             → 1031
cargo exit=101
    cli::tests::hardtrim5_range_error_precedes_the_library_refusal
test result: FAILED. 480 passed; 1 failed
```

**D3's equivalence claim holds.** One red rather than §12's "2" is correct and not a discrepancy: the
`#421` sibling only goes red under the *bail-suppression* variant, which changes the fastqc ordering
too; the relocation leaves that ordering alone.

### V15 — docs transcripts ✅ every quoted line verbatim

`trim_galore --library emseq --clip_R1 12 -o . sample.fq.gz`, with `grep -qF --` so the
leading-`--` hazard does not apply:

```
=== stderr (docs 53-54) ===
PRESENT: Library preset 'emseq' selected: --clip_R1 12 --three_prime_clip_R1 10
PRESENT: --clip_R1 12 was given on the command line and overrides the emseq preset value 10
=== text report (docs 68-70) ===
PRESENT: Library preset: emseq
PRESENT: Clipping in force: --clip_R1 12 --three_prime_clip_R1 10
PRESENT: --clip_R1 12 was given on the command line and overrides the emseq preset value 10
=== --clip_R2 count in stderr / report ===
0
0
```

The JSON block reproduced docs `:76-83` exactly — same keys, same order, same values — and `mode` was
`single-end`, so the sentence at `:86-88` is accurate.

---

# 7. Findings

## DEFECTS

### F-1 — **High** — a top-level JSON key permutation survives all 764 tests

`src/report.rs:907-921` (writer) · `src/report.rs:1489-1531` and `tests/integration_library.rs:203-214`
(the tests that should catch it)

The commit closes T-1 for the `library` **block** (V8, reproduced above) but not for the four
**top-level** clip keys. Proved by mutating the top-level labels at `report.rs:907` and `:912` — the
exact analogue of V8:

```
S8 mutation applied: TOP-LEVEL clip_r1 / three_prime_clip_r1 key labels swapped
cargo exit=0
=== totals (non-zero failures only) ===
(no output)
=== TOTAL PASSED ===
764
```

**Zero tests fail.** Two reasons, both inside code this commit wrote:

1. The new F-5b agreement loop (`tests/integration_library.rs:203-214`) runs on `pbat`, whose four
   values are **8/8/8/8** — so it cannot see a permutation on *either* side of the equality. It catches
   value drift (V12's `+100`) and nothing else. This is precisely the degenerate fixture the plan itself
   rejected for the unit test: §7 says "`pbat` (8/8/8/8) catches none". The T-1 lesson was applied in
   one place and not the other.
2. `test_json_report_records_all_four_overrides` **already builds** `extra` with distinct 1/2/3/4
   (`report.rs:1504-1510`) and then asserts only `json["parameters"]["library"][…]`, never
   `json["parameters"][…]`. The distinct values are constructed and thrown away.

This also means the commit message's "gating it would break the block/top-level agreement invariant
added here" overstates what the added test can detect. The *reasoning* still stands; the *guard* is
weaker than the sentence implies.

**Fix (4 lines, no new fixture needed)** — in `test_json_report_records_all_four_overrides`, after the
existing block assertions:

```rust
// The block and the top level describe the same resolved configuration.
for (k, v) in [("clip_r1", 1), ("clip_r2", 2), ("three_prime_clip_r1", 3), ("three_prime_clip_r2", 4)] {
    assert_eq!(json["parameters"][k], v, "top-level {k}");
    assert_eq!(json["parameters"][k], lib[k], "{k} disagrees with the library block");
}
```

That kills the S8 mutation, kills a top-level *value* swap, and makes the agreement invariant real
rather than fixture-limited.

### F-2 — **Medium** — a test comment asserts an invariant the tool does not hold

`tests/integration_library.rs:147-152`

```rust
// The banner is the other consumer of the same summary, and this is a
// single-end run: no Read 2 clipping may be announced anywhere on stderr.
assert!(
    !err.contains("--clip_R2") && !err.contains("--three_prime_clip_R2"),
```

"anywhere on stderr" is false. Measured, single-end, `--library pbat --clip_R2 5`:

```
Library preset 'pbat' selected: --clip_R1 8 --three_prime_clip_R1 8
--clip_R2 5 was given on the command line and overrides the pbat preset value 8
```

`--clip_R2` **is** announced on stderr — by the override line, which the plan deliberately leaves
ungated (§B2). The assertion passes only because this particular invocation happens to override
`--clip_R1`. Two consequences:

- The comment states a rule the code does not enforce, which is the CLAUDE.md failure mode ("state the
  fact") and the kind of sentence that gets copied into the next test.
- The assertion is coupled to the choice of override flag rather than to the thing it means to pin. It
  is the V6a killer, so it should be robust: target the banner line specifically, e.g. assert over the
  line that starts `Library preset '` rather than over all of stderr.

Reword to something the code actually guarantees — "the banner lists no Read 2 value on a single-end
run; the override line below is about what was typed, not what applies" — and narrow the assertion.

### F-3 — **Medium** — `--library`'s `--help` text was not updated with the refusal

`src/cli.rs:138-142`

The plan's step 23 covered `flags.md` and never considered the CLI help, which is the surface most
users read. Current text ends "...both numbers are named in the log and in the trimming reports.
Presets may change between releases; see CHANGELOG.md." — no mention of the refusal, and no mention
that the text report lists only the values in force.

This departs from the project's own established convention, set by the **very precedent this commit is
modelled on**:

```
--fastqc
    … Rejected with the specialty modes (--hardtrim5/3, --clock, --implicon) and with --paired
    given a single input file (the interleaved-uBAM path) writing FASTQ output, none of which
    write a report. …

--rename
    … Refused with --output-format ubam when any input is FASTQ, … also refused by --clump_only,
    which does not mutate IDs
```

Add one clause. Suggested: "Refused with `--hardtrim5`/`--hardtrim3`/`--clock`/`--implicon` and with
`--clump_only`, none of which honour the clip flags."

### F-4 — **Medium** — the neighbouring #421 comment is now false

`src/cli.rs:1071-1072`

```rust
// #421 — these arms reach no `fastqc::run` call site. Last of the fallible checks,
// so malformed values and invalid input lists report their own defect first.
```

After this commit the `--library` block at `1112-1134` is the last of the fallible checks, and the new
block's own comment does not say so. CLAUDE.md: "Current state only — never reference the state you are
fixing", and by the same logic a comment must not describe a state that has been superseded. Drop
"Last of the fallible checks," from `:1071` (the "so …" clause is still true of both blocks), or move
the phrase onto the new block.

### F-5 — **Medium**, scope decision needed — the refusal covers the derived flag, not the primitives

`src/cli.rs:1112-1134` vs the absence of any guard for `--clip_*` outside `if self.clump_only`

The error tells the user "`--hardtrim5` does not honour the clip flags". The workaround it implicitly
points at — writing those flags out by hand, which is exactly what `--library` exists to remove — is
still accepted, silently:

```
### hardtrim5 + raw clip flags (the literal expansion of --library emseq)
$ trim_galore --hardtrim5 30 --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 -o . R1.fq.gz
exit=0
(no warning of any kind on stderr)
files: BS-seq_10K_R1.30bp_5prime.fq.gz
len=30                       ← all four clip values discarded
```

```
$ grep -n "clip_r1\|clip_R1" src/cli.rs | grep -v "test\|library\|clump"
   → the only clip-flag refusals in the file are inside `if self.clump_only` (833, 839)
```

Three things make this worth recording rather than shrugging at:

- It is the same user-visible defect class as #440 — a clip request silently discarded — reachable by a
  different spelling of the same request.
- The adjacent `-a2` family already handles this situation, and handles it more cheaply: `cli.rs:1049-1068`
  emits `WARNING: -a2/--adapter2 was given but is not used in this mode (--clock performs no adapter
  trimming). Ignoring.` The clip flags get nothing.
- The `--phred64` comment the commit message cites as its precedent (`main.rs:273-277`) argues the
  opposite way to a partial application: *"Rejected uniformly rather than per-mode. `--hardtrim5/3` and
  `--clump_only` currently accept the flag as an inert no-op, so this does remove working invocations —
  but a mode-dependent rule … is worse to document, worse to test, and one refactor away from being
  wrong."*

Not introduced here and out of this commit's agreed scope — but this commit's error text is what makes
it conspicuous, and a user who hits the refusal is being steered straight into it. Per Felix's filing
bar this is a user-visible correctness gap rather than a review spin-off, so it is at least an issue;
a warning (the `-a2` shape) would close it without removing working invocations.

## STYLE / EDITORIAL OBSERVATIONS

### F-6 — **Low** — §B2's justification for leaving the override lines ungated does not fit them

The plan says the override line "is true about what the user typed and claims nothing about
application". But the verb is *overrides*, and the pair now reads:

```
Library preset 'pbat' selected: --clip_R1 8 --three_prime_clip_R1 8
--clip_R2 5 was given on the command line and overrides the pbat preset value 8
```

The second line claims a resolution outcome for a value the first line has just (correctly) declined to
list. It is literally defensible — `resolve()` really did pick 5 over 8, and `clips.clip_r2 == Some(5)` —
but a reader takes it as "my 5 is in effect", and on a single-end run nothing is. The gate made one of
the two lines honest and left the other making the stronger claim.

No code change recommended: behaviour is identical either way, and gating the override lines would lose
the "you typed this" record the plan wanted. But §B2's stated reason should be corrected to the real
one — *the line records the resolution, not the application* — so the next person does not inherit the
argument that it "claims nothing about application".

### F-7 — **Low** — docs and CHANGELOG list four of the five modes that refuse `--library`

`--clump_only` also refuses it (`src/cli.rs:850`), and neither page mentions `--clump_only` at all:

```
$ grep -n "clump" docs/src/content/docs/guide/library-presets.md docs/src/content/docs/guide/flags.md
(no output)
```

A reader of "It is refused with `--hardtrim5`, `--hardtrim3`, `--clock` and `--implicon`"
(`flags.md:13`) reasonably infers it is accepted with `--clump_only`. One word in each place.

### F-8 — **Low** — CHANGELOG's "For the same reason" is a loose connector

`CHANGELOG.md:18`. The preceding reason is "a preset is only a clip bundle", which is why the *mode*
flags stay separate flags. The refusal's reason is that those modes *do not honour the clip flags* —
related, but not the same proposition. Also, the paragraph now names `--clock` and `--implicon` twice
in two different roles four lines apart.

Separately: CLAUDE.md's "One sentence per entry, no exceptions" — this entry is a multi-paragraph block
inherited from #441 and the edit grows it by three lines. The plan consciously chose in-place editing
over a rewrite, and I agree that was the right call for the false sentence; but since these exact lines
were being touched anyway, it was the cheap moment to shrink them. Felix's call.

### F-9 — **Low** — `config.paired = true` in `test_json_report_records_all_four_overrides` is inert

`src/report.rs:1491`. The JSON is deliberately not gated on `paired`, and the test asserts nothing that
depends on it, so the line reads as though `paired` mattered to the assertions below — the opposite of
the point the test exists to make. It can be read as deliberate fixture realism; if kept, say so in one
clause, or drop it. (It would become load-bearing under F-1 only if a `mode` assertion were added.)

### F-10 — **Low** — two test doc comments justify the change rather than state the fact

`tests/integration_library.rs`: `/// … Guards against a gate that drops Read 2 unconditionally.` and
the two-line comment at `:147-148` ("The banner is the other consumer of the same summary…"). CLAUDE.md's
rule covers committed tests, and both sentences are arguments about the change rather than statements
about the current behaviour. The rationale is already in the commit message, where it belongs. Minor —
and F-2 requires rewriting one of them anyway.

---

# 8. Explicitly verified as correct

So that absence of a finding is distinguishable from absence of a check.

**Refusal**
- Sited as the last fallible statement in `validate()` — verified structurally by enumerating every
  `bail!` in the file, not by reading the plan.
- All five precedence guarantees the plan claims, each run against the binary (range, #386, `--clock`
  arity, `--implicon` arity, missing input) — plus two the plan does not mention (clap-level bad preset
  value; `--clump_only`'s short-circuit).
- Exit 1 and **zero files written** on all four modes.
- Refusal fires on the equals form `--implicon=8` and on multi-pair input for all three applicable modes.
- `validate()` precedes every specialty dispatch (`main.rs:237` vs `:468`).
- Message shape consistent with the `--clump_only` and `-a2` families.

**Wording**
- All four modes do clip — `clip_5prime` at `specialty.rs:107,219`, direct slicing at `:314,:319,:413`,
  documented at `:241,:352`.
- `specialty.rs` has **zero** references to `clip_r1`/`clip_r2`/`three_prime_clip_*`/`effective_clips`.
- Therefore "does not honour the clip flags" is accurate and "does not clip" would have been false.

**Gate**
- Exactly two non-test call sites; both updated.
- `config.paired` is set literally at all three `library`-bearing `TrimConfig` sites; there is no fourth.
- `write_paired_reports` is the single funnel for all four paired shapes (4 call sites, one function).
- Banner and text report agree, measured on **six** distinct shapes: SE FASTQ, PE FASTQ, SE uBAM-out,
  PE uBAM-out, interleaved-uBAM-in → FASTQ, interleaved-uBAM-in → uBAM.
- `--clump_only` cannot reach the gate.
- `flag_summary`'s new doc comment is two lines, states the fact, carries no rationale — the two
  justification sentences were removed as the plan required.

**JSON / non-regression**
- Zero changes to the JSON writer in the diff; no key added, removed, renamed or reordered.
- `library` block carries all four values on a single-end run (measured), with `mode: single-end`.
- `json["mode"]` exists and derives from the same `config.paired` — §8.3's consumer argument is sound.
- A no-preset run cannot reach any changed line (reachability argument stated in §5, plus a measured
  `--rrbs --paired` run showing 0 new lines and `"library": null` intact, RRBS auto-clip still firing).

**Tests**
- 764 pass, 17 binaries, `--no-fail-fast`, 0 failures.
- 6 matrix rows re-derived from scratch, all exact: V2, V3 (as written, not D3's substitute), V6a, V7,
  V8, V11, V13.
- V11's panic messages confirm the *new* presence assertion is what fires, so "both stayed green
  before" is corroborated rather than assumed.
- V13 confirms keeping the cardinality-1 test was correct.
- Every mutation restored and `cmp`-verified byte-identical; `git diff --stat src/` empty after each.

**Docs / CHANGELOG**
- Every line quoted in `library-presets.md` reproduced verbatim from a real run (`grep -qF --`).
- The JSON block at `:76-83` matches real output key-for-key and value-for-value.
- CHANGELOG: **70 unique `#N` refs before and after, identical sets** — `comm -23` and `comm -13` both
  empty. 39 link-form `[#N]` refs before and after. (The brief's figure of 52 does not match either
  count in this file; what matters is that the sets are identical.)
- The 5 `^-` lines are the reflowed prose of the entry being corrected — `--word-diff` shows no content
  dropped.
- The corrected sentence is now true: all four in JSON, values in force in the text report, single-end
  lists Read 1 only.

**Commit message** — checked claim by claim; every factual assertion holds.
- "the specialty arms return before `setup_trimming`, which is the only place a preset reaches a clip
  site" — `effective_clips()` has five consumers (`main.rs:393` the #408 guard, `:1116` setup_trimming,
  `:1471/2261/2640` report configs); only `:1116` feeds a clip site. ✓
- "Only `--clump_only` rejected it" ✓ (`cli.rs:850`).
- The #421 and #358 precedents are as described (`cli.rs:1073-1110`, `main.rs:273-277`) ✓.
- "`pbat` is 8/8/8/8" ✓; "12 new tests, 764 passing" ✓; "both call sites updated" ✓.
- One caveat, not an error: "A swapped JSON key label … now kills a test" is true of the block and says
  nothing about the top level, where it is still false (**F-1**).

---

# 9. Worktree state

```
$ git status --short
?? B_baseline.log
?? B_clippy.log
?? B_clippy_canary.log
?? B_mut_S8.log
?? B_mut_V11.log
?? B_mut_V13.log
?? B_mut_V2.log
?? B_mut_V3.log
?? B_mut_V6a.log
?? B_mut_V7.log
?? B_mut_V8.log
?? cli.rs.orig
?? library.rs.orig
?? main.rs.orig
?? mutM.log
?? mutcheck.sh
?? refs_new.txt
?? refs_old.txt
?? report.rs.orig
?? scratch/
?? testsummary.log
?? verify1.log
```

**No tracked file modified.** Every entry is untracked review scratch: my mutation logs, the four
`*.orig` snapshots used to restore mutated sources, the CHANGELOG ref lists, and `scratch/` (all binary
runs used `-o` into `scratch/`, never `test_files/`). `mutM.log`, `mutcheck.sh`, `testsummary.log` and
`verify1.log` are the pre-existing leftovers the brief said to ignore.

`cargo clean` run at the end of the review: `Removed 7623 files, 1.0GiB total`. Volume back to 7.1 GB
free. `git status --short` re-run afterwards: still 22 untracked entries, 0 tracked modifications.

**Review complete.** No GitHub writes were made; nothing was pushed, commented on, or opened.
