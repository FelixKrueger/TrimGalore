# Code Review — PR #374 (`fix(io): read gzipped input whose filename does not end in .gz`)

**Reviewer B** · base `dev` (`7201f03`) · head `7a74458` · +225/−58 across `CHANGELOG.md`, `src/fastq.rs`, `src/format.rs`, `src/main.rs`

Worktree used for all builds and runs:
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/db842fb4-f582-4e29-8459-cfbe1cf77ceb/scratchpad/wtB`

---

## Summary

**Merge recommendation: do not merge as-is. Request changes.**

The diagnosis is correct, the mechanism is sound, and the single-end fix is genuinely
right — verified byte-identical to the plain-input equivalent. The PR does not weaken
BAM discrimination and does not disturb byte-identity on any existing input.

But the fix is **incomplete in a way the PR title and CHANGELOG both overclaim**, and on
the most common paired invocation it is a **net regression**. `FastqReader::open` /
`open_threaded` / `sanity_check` have callers the PR did not update. The three that
matter:

| Path | Site | Status on PR head |
|---|---|---|
| `--paired` at `--cores 1` (**the default**) | `src/main.rs:1348-1349` | still broken, **and now destroys pre-existing output** |
| `--clump_only` | `src/clump_only.rs:290`, `430`, `432` | still broken, leaves a partial output file |
| `--passthrough` | `src/main.rs:743`, `1319`, `1358` | still broken (fails early, no files) |

The paired case is the serious one. Pre-fix, the run aborted in `sanity_check_any`
before any writer was constructed. Post-fix, `sanity_check_any` passes, so execution
reaches `FastqWriter::create` — which **truncates the output files** — and only then hits
the same `stream did not contain valid UTF-8`. I reproduced this destroying real prior
content (evidence in *Errors E1*). The PR's own description cites "before any output
directory was created" as a property of the bug; the PR removes that property on the path
it did not fix.

CI is green because all six new tests are unit-level on the two layers that were changed.
There is no integration test that runs the binary on a `.bgz` file, which is precisely
the level at which the bug was reported and at which the paired gap is visible.

Separately, one new test actively **pins the bug as a contract** — I proved it is the only
test that fails when you make `FastqReader::open` correct-by-default (*Logic L4*).

### Verified claims (all checked, not taken on trust)

| Claim | Verdict |
|---|---|
| Original bug exists on base `dev` | ✅ reproduced, SE and PE, `stream did not contain valid UTF-8` |
| PR fixes single-end | ✅ and output is byte-identical to plain/`.gz` equivalents, 1000 reads |
| PR fixes paired-end | ❌ only at `--cores >= 2`; broken at the `--cores 1` default |
| Deleted `main.rs` fns were byte-identical duplicates | ⚠️ behaviourally identical, *not* byte-identical — safe to delete |
| Verdict crosses the thread boundary correctly | ✅ `bool` is `Copy`, moved with the owned `PathBuf`; verified at `--cores 2` |
| BAM-vs-FASTQ discrimination not weakened | ✅ a `.bam` renamed `.fq.gz` is still read as BAM |
| Byte-identity preserved on existing inputs | ✅ 10 outputs across SE / PE-serial / PE-parallel / RRBS / hardtrim5 / clock, identical |
| Pre-existing non-seekable double-open bug | ✅ PR is exactly **neutral** — same open count, same messages |
| `cargo test` / `fmt` / `clippy -D warnings` | ✅ 375 unit + 87 integration pass; fmt and clippy clean |

---

## Issues by area

### Logic

**L1 (Critical) — The fix covers three *functions* but not their callers.**

The CHANGELOG states the verdict "was then discarded and re-derived from the *filename*
in three places (`FastqReader::open`, `FastqReader::open_threaded`, and
`FastqReader::sanity_check`)". That is an accurate description of the *mechanism* and a
misleading description of the *coverage*. Those three functions were reached from more
call sites than the two factories plus `sanity_check_any`. On PR head:

```
src/main.rs:743           FastqReader::sanity_check(pt_path)?          --passthrough
src/main.rs:1319          FastqReader::open_threaded(p)                --passthrough, cores>1
src/main.rs:1348-1349     FastqReader::open(input_r1 / input_r2)       --paired, cores==1  ← DEFAULT
src/main.rs:1358          FastqReader::open(p)                         --passthrough, cores==1
src/clump_only.rs:290     FastqReader::open(input)                     --clump_only SE
src/clump_only.rs:430/432 FastqReader::open(input_r1 / input_r2)       --clump_only PE
```

(`src/demux.rs:170` reads our own trimmed output and `src/specialty.rs:531`,
`src/clump_only.rs:1127`, `src/parallel.rs:*` are tests or self-produced files — those are
fine.)

Evidence, real `bgzip` input on the PR-head binary:

```
$ bgzip -c plain_R1.fastq > bg_R1.fq.bgz     # file(1): "Blocked GNU Zip Format (BGZF)"

$ trim_galore --paired -o out bg_R1.fq.bgz bg_R2.fq.bgz
Error: processing pair 1 of 1 (R1=bg_R1.fq.bgz, R2=bg_R2.fq.bgz)
Caused by: stream did not contain valid UTF-8          ← exit 1

$ trim_galore --clump_only --cores 2 -o out bg_R1.fq.bgz
clump-only: reordering 'bg_R1.fq.bgz' -> 'out/bg_R1.fq_clumped.fq' ...
Error: stream did not contain valid UTF-8              ← partial output file left behind

$ trim_galore --paired --cores 2 --passthrough index_I1.fq.bgz -o out bg_R1.fq.bgz bg_R2.fq.bgz
Error: stream did not contain valid UTF-8
```

`--cores 2` paired works, confirming the split is exactly serial-vs-threaded:

```
$ trim_galore --paired --cores 2 -o out bg_R1.fq.bgz bg_R2.fq.bgz
  → bg_R1.fq_val_1.fq  190535 bytes   ✅
```

**`src/clump_only.rs` is the sharpest instance.** `clump_only_single` already calls
`detect_input_format(input)` at line 265 (to reject uBAM), discards the verdict, then
re-derives it from the filename at line 290 — and derives the same fact a *third* time
from the filename at line 277 via `naming::is_gzipped`. The correct verdict is in scope,
25 lines above the broken call.

**L2 (Medium) — A second, undocumented behaviour change ships in the same diff.**

The change is bidirectional. A plain FASTQ *misnamed* `.fastq.gz` also changes behaviour:

```
                        pre-fix (dev)              post-fix (PR head)
plain content           Error: invalid gzip        SUCCESS →
named .fastq.gz         header (no output)         misnamed_trimmed.fq.gz
```

This is arguably the *right* behaviour (content wins), and the new test
`open_with_honours_a_plain_verdict_on_a_gz_name` covers the unit. But the end-to-end
consequence is not in the CHANGELOG, and the CHANGELOG's own reasoning cuts against it:
it defers the output-side change because doing it "would change behaviour for every run
(**including a plain file misnamed `.gz`**)". That is precisely the file whose *input*
behaviour this PR already changed. The PR half-applies content-authority and justifies
skipping the other half on a ground that already applies to the half it did.

Not a defect — flag it in the CHANGELOG so it isn't discovered in the field.

**L3 (Medium) — `io::is_gzipped`'s doc comment now asserts the opposite of the truth.**

`src/io.rs:12-18`, untouched by the PR:

> True iff `path` ends with a `.gz` extension. **The same heuristic that `FastqReader`
> uses to decide whether to wrap the input in a gzip decoder, so output naming + reader
> behaviour stay consistent.**

After this PR that is false — the factories override the filename heuristic, and output
naming and reader behaviour deliberately *diverge*. That divergence is the PR's own
"known limitation", and the one place a maintainer would look to understand it now
documents the pre-PR world. This is the highest-value two-line comment fix in the diff.

**L4 (High) — One new test pins the bug as a contract, and I proved it blocks the fix.**

`src/fastq.rs:1009`:

```rust
// And the filename-based entry point still gets it wrong, which is
// exactly why callers with a content-detected verdict must not use it.
assert!(FastqReader::sanity_check(&p).is_err());
```

Two problems.

First, `is_err()` asserts almost nothing — it would pass if the file were missing, if
permissions failed, or if the message changed entirely. It does not pin *why*.

Second and more important: it makes the correct fix look like a regression. I applied a
12-line experimental patch turning `is_gz_filename` into a 3-byte magic sniff, and ran
the suite. Result:

```
test fastq::tests::sanity_check_with_accepts_gzip_under_non_gz_extension ... FAILED
panicked at src/fastq.rs:1018:9
test result: FAILED. 374 passed; 1 failed
```

**One test failed, and it was this one.** Every other test in the crate is indifferent to
the improvement. That is the definition of a test that encodes a defect rather than a
requirement.

The same experiment fixed all three gaps from L1 in one change:

```
--paired --cores 1  .bgz  → bg_R1.fq_val_1.fq  190535 bytes   ✅
--clump_only        .bgz  → wrote 1000 records in 16 bins      ✅
--passthrough       .bgz  → val_1 / val_2 / reports written    ✅
```

Experiment reverted; worktree confirmed clean at `7a74458` (see *Fixes applied*).

**L5 (Low) — Answering the question posed about the pre-existing double-open bug: the PR
is exactly neutral.** It does not touch `detect_input_format`, changes no open counts, and
produces identical failure text. Confirmed on both binaries:

```
gzip -c reads | trim_galore -o out /dev/stdin
  pre:  Error: Failed to decompress first block of '/dev/stdin' for format detection
        Caused by: invalid gzip header
  post: identical

cat reads | trim_galore -o out /dev/stdin
  pre:  Error: Input file '/dev/stdin' doesn't seem to be in FastQ format ...
  post: identical
```

Note the gz-via-pipe shape is a **detection-stage** failure (`invalid gzip header` at
`format.rs:231`, because the 4-byte peek at `:207` is unrecoverable), distinct from the
plain-via-pipe shape (`first line doesn't start with '@'`). Two manifestations, one root
cause.

One forward-looking caveat: the PR threads the verdict as a `bool`. A proper seekability
fix threads a *buffered reader* (the `PeekReader` already noted as a TODO in the
`format.rs:207` comment), which would supersede the bool plumbing. The PR mildly entrenches
a shape that the eventual fix replaces. Not a reason to block.

### Efficiency

No concerns. `open_with` performs exactly one `File::open`, the same as the code it
replaces; `is_gz_filename` is called only in the three thin wrappers and is a string
compare. `open_reader` de-duplicates the decoder construction that was previously
copy-pasted between `open` and `open_direct`, which is a small genuine improvement. The
`bool` crossing the thread boundary is `Copy` — no allocation, no clone.

The `fmt == InputFormat::FastqGz` comparisons in `format.rs:265`/`283` could be
`matches!(fmt, InputFormat::FastqGz)` to avoid relying on `PartialEq`, but this is
cosmetic and costs nothing either way.

### Errors

**E1 (Critical) — The paired path now destroys pre-existing output before failing.**

This is the finding I would block on. Same command, same input, same output directory,
two binaries:

```
# seed the output dir with "good output from a previous run"
$ echo "PRECIOUS PREVIOUS OUTPUT - 12345" > d/bg_R1.fq_val_1.fq   # 33 bytes
$ echo "PRECIOUS PREVIOUS OUTPUT - 12345" > d/bg_R2.fq_val_2.fq   # 33 bytes

# PRE-FIX (base dev, c5c1834)
$ trim_galore --paired -o d bg_R1.fq.bgz bg_R2.fq.bgz
Error: stream did not contain valid UTF-8
$ ls -la d | grep val
  33 bg_R1.fq_val_1.fq        ← preserved
  33 bg_R2.fq_val_2.fq        ← preserved
$ cat d/bg_R1.fq_val_1.fq
PRECIOUS PREVIOUS OUTPUT - 12345          ✅ intact

# POST-FIX (PR head, 7a74458)
$ trim_galore --paired -o d2 bg_R1.fq.bgz bg_R2.fq.bgz
Error: processing pair 1 of 1 ... Caused by: stream did not contain valid UTF-8
$ ls -la d2 | grep val
   0 bg_R1.fq_val_1.fq        ← TRUNCATED
   0 bg_R2.fq_val_2.fq        ← TRUNCATED
$ cat d2/bg_R1.fq_val_1.fq
                                          ❌ destroyed
```

Mechanism: fixing `sanity_check_any` lets execution proceed past the entry guard into
`run_paired`, where `FastqWriter::create` (`src/main.rs:1350-1351`) truncates both output
paths, and the failure only surfaces on the first `next_record()`. The output-collision
pre-flight does not help — it guards against *collisions between* prospective outputs, not
against clobbering on a run that is going to fail.

Mitigating: the exit code is correctly `1`, and the failure is loud. There is **no
silently wrong output** — I checked. So this is data destruction on a failed run, not
corrupt results.

**E2 (Medium) — `--clump_only` leaves a partial output file on the same failure.**
`clump_only.rs` creates the output via `File::create` (~line 297) before the read fails,
so `bg_R1.fq_clumped.fq` is left behind. Same class as E1, smaller blast radius (opt-in
flag).

**E3 (Low) — No integration coverage at the level where the bug was reported.**
`tests/` has an established pattern (8 integration files, including
`integration_paired_format_guard.rs`, which is the natural neighbour). Nothing invokes the
binary on a `.bgz` input. A single integration test running SE + PE-serial + PE-parallel
against a `.fq.bgz` fixture would have caught L1/E1 and would have prevented the green CI
from being read as "fixed".

### Structure

**S1 (Medium) — `is_gz_filename` duplicates `io::is_gzipped` character-for-character.**

```rust
// src/io.rs:19  (pre-existing, public)
pub fn is_gzipped(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext == "gz")
}

// src/fastq.rs:37  (new, private)
fn is_gz_filename(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext == "gz")
}
```

Two names for one predicate, in a PR whose whole subject is "one fact was derived in
several places". The new doc comment is good and worth keeping — but it belongs on
`io::is_gzipped` (replacing the now-false claim in L3), with `fastq.rs` calling it.

**S2 (Low) — Naming: `open_with` reads ambiguously at call sites.**
`FastqReader::open_with(&p, true)` — `true` meaning what? The crate already has
`InputFormat`; a `Compression { Plain, Gzip }` two-variant enum local to `fastq.rs` (to
avoid a `fastq` → `format` module dependency) would make every call site self-documenting
and cost nothing. This is preference, not defect — but the boolean shows up in six new
test call sites where it reads worst.

**S3 (Low) — `open_direct` is now a two-line pass-through** over `open_reader` with one
caller, and the `String::with_capacity(512)` literal appears in both `open_with` and
`open_direct`. Pre-existing duplication of the `512`; the pass-through is new. Harmless,
and keeping it does minimise the diff.

**S4 (Low) — Test-helper duplication and asymmetry.** `gz_tmpdir` (`fastq.rs:945`) is an
exact copy of `fresh_tmpdir` (`format.rs`, pre-existing). Also
`threaded_reader_decompresses_gzip_under_non_gz_extension` omits the
`assert_eq!(detect_input_format(&p)?, InputFormat::FastqGz)` precondition that its sync
twin has — harmless, but the pair reads as if one is less complete.

### On the "byte-identical duplicates" claim

Verified by token-diffing both functions out of `git show dev:src/{main,format}.rs`, with
a sentinel to prove the comparison can fail:

```
open_threaded_reader: TOKEN-IDENTICAL  (87 vs 87 tokens)
open_sync_reader:     DIFFERS          (91 vs 90 tokens)
sentinel check (must say DIFFERS): DIFFERS      ← comparison is not blind
```

The single token is a trailing `,` — `main.rs` had the match arm as a bare expression,
`format.rs` had it brace-wrapped (rustfmt, because `crate::fastq::` pushes the line past
100 columns). Doc comments and path qualification also differ. So: **not byte-identical,
but behaviourally identical** — same match arms, same dispatch, same `preserve_tags`
handling. Deleting them is safe, and consolidating on the `format.rs` pair is the right
call. Only the PR's wording is imprecise.

### On the stated known limitation (output compression follows the filename)

Deferring it is **safe and correct**, and I would not block on it. A `.bgz` input yields
`reads.fq_trimmed.fq` — complete, correct, uncompressed data. Verified byte-identical to
the plain-input run:

```
plain input -> 41a9413c8850d6e8b68dc2bd9d329826
.bgz  input -> 41a9413c8850d6e8b68dc2bd9d329826
.gz   input -> 41a9413c8850d6e8b68dc2bd9d329826    (1000 reads each)
```

Uncompressed-but-correct output is strictly better than a failed run, the author flagged
it explicitly, and moving the output decision has real byte-identity risk against the
Perl validation matrix. Correct scoping. It is worth one sentence in the docs (not just
the CHANGELOG) since users will notice the missing `.gz`.

### Byte-identity against existing inputs

10 output files across SE, PE-serial, PE-parallel, RRBS, hardtrim5 and clock, comparing
base-`dev` and PR-head binaries on the repo fixtures (gzip-decompressed md5, reports
excluded because they embed version and cmdline):

```
=== files compared: 10 pre / 10 post ===
RESULT: IDENTICAL
sentinel OK: the comparison does detect differences
```

No existing input is read differently. The Perl v0.6.11 invariant is not at risk — which
matches the green `Validate vs Perl TrimGalore` job. BAM discrimination also holds: a
`.bam` renamed `.fq.gz` is still read as BAM (20 reads via `open_sync_reader`), so the
content-detection invariant is strengthened, not weakened.

---

## Fixes applied

**None.** This is an external contributor's branch; per instruction I recommended rather
than edited, including for the trivial items.

One experimental change was made and reverted: I temporarily rewrote `is_gz_filename`
(`src/fastq.rs:37`) as a 3-byte magic sniff to test recommendation **R1** below, built,
ran the affected commands and the unit suite, then reverted with
`git checkout -- src/fastq.rs` and **rebuilt from the reverted source**. Worktree confirmed
clean at `7a74458`; `grep EXPERIMENT src/fastq.rs` returns nothing; the binary's provenance
line reads `7a74458`. Nothing in the main checkout at `/Users/fkrueger/Github/TrimGalore`
was modified.

*(Method note, for whoever repeats this: I initially ran a smoke test after reverting the
source but before rebuilding, and got a false pass from the stale experimental binary. The
truncation evidence in E1 was re-run against a freshly built PR-head binary. This is the
`cargo` trap listed in the brief and it bites in the revert direction too.)*

---

## Recommendations

### Critical

**R1 — Fix the remaining call sites before merging.** Two options.

*Preferred — make the default correct, ~12 lines, fixes all three gaps at once.* Replace
the body of `is_gz_filename` (`src/fastq.rs:37`) with a magic-byte sniff, so
`FastqReader::open` / `open_threaded` / `sanity_check` are right for every caller,
including ones added later:

```rust
/// Content-based gzip check for constructors with no caller-supplied verdict.
/// Cheap enough (3 bytes) that filename guessing buys nothing.
fn sniff_gzip(path: &Path) -> bool {
    let mut m = [0u8; 3];
    File::open(path)
        .and_then(|mut f| std::io::Read::read(&mut f, &mut m))
        .is_ok_and(|n| n == 3 && m == [0x1F, 0x8B, 0x08])
}
```

I verified this fixes paired-`--cores 1`, `--clump_only` and `--passthrough`, and that no
test except the one in L4 objects. Deliberately no BAM logic, so no `fastq` → `format`
module dependency. Cost is one 3-byte read per file open.

*Alternative — thread the verdict explicitly at each site.* Update `main.rs:1348`, `1349`,
`1319`, `1358`, `743` and `clump_only.rs:290`, `430`, `432` to `*_with` forms.
`clump_only_single` already has the verdict in scope at line 265. More explicit, matches
the PR's existing style, but leaves the footgun live for future callers — so if you take
this route, also do **R3**.

**R2 — Do not let a failed run truncate output.** Whichever form R1 takes, the paired
serial path should not reach `FastqWriter::create` on input it cannot read. R1 fixes the
`.bgz` instance; the general shape (guard passes → writers truncate → read fails) is worth
a follow-up issue independent of this PR, since any post-guard read failure has the same
consequence.

### High

**R3 — Remove the footgun rather than testing it.** Delete
`assert!(FastqReader::sanity_check(&p).is_err())` (`src/fastq.rs:1009`). If R1's preferred
form is taken it becomes false anyway. If the explicit-threading alternative is taken,
make the filename entry points `#[deprecated]` or `pub(crate)`, or delete them so the
compiler finds every unconverted caller — the failure mode here is a *silent* wrong answer,
which is exactly what a type system should be made to catch.

**R4 — Add one integration test at the level the bug was reported.** A new
`tests/integration_bgz_input.rs` (or a case in `integration_paired_format_guard.rs`)
invoking the built binary on a `.fq.bgz` fixture for SE, PE `--cores 1`, and PE
`--cores 2`, asserting success and correct record counts. The six unit tests are good
tests of the new API; none of them can observe L1.

### Medium

**R5 — Fix the stale doc comment at `src/io.rs:12-18`.** It now claims the opposite of
the truth. Fold the new `is_gz_filename` docs into it and state the divergence plainly.

**R6 — Collapse `is_gz_filename` into `io::is_gzipped` (S1)** — or, if R1's sniff form
lands, keep the sniff in `fastq.rs` and update `io::is_gzipped`'s doc to say it is the
*output-naming* heuristic only, which is then its actual sole remaining job.

**R7 — Add the reverse-direction behaviour change to the CHANGELOG (L2).** One sentence:
a plain FASTQ misnamed `.gz` now reads successfully instead of erroring.

### Low

- **R8** — Reword the PR/CHANGELOG "byte-identical" claim to "behaviourally identical"; note the trailing-comma/brace difference so the next reader isn't misled into skipping the check.
- **R9** — Consider a `Compression { Plain, Gzip }` enum instead of the `bool` (S2).
- **R10** — Share one tmpdir helper between the `fastq.rs` and `format.rs` test modules, and add the missing `detect_input_format` assertion to the threaded factory test (S4).
- **R11** — Document the `.bgz` → plain-output limitation in the user docs, not only the CHANGELOG.
- **R12** — File the non-seekable-input double-open bug (`format.rs:207`/`:231`) as its own issue. The PR is neutral on it, but the `PeekReader` TODO already in that comment is the fix, and it would also resolve `<(zcat …)` and `/dev/stdin` input.

---

## Bottom line

Good diagnosis, correct mechanism, real fix for single-end, no collateral damage to
byte-identity or BAM detection, clean CI hygiene. But it ships under a title that promises
more than it delivers, and on the default paired invocation it converts a harmless early
abort into one that truncates the user's previous output. **R1 + R3 + R4** — a ~12-line
change, one deleted assertion, one integration test — turn it into a complete fix. Worth
sending back to the contributor with the evidence rather than merging and following up,
because the CHANGELOG as written would tell users the paired case is fixed.
