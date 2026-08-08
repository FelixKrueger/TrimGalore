# Code Review A — #388 paired report pre-flight

**Reviewer A** (independent; Reviewer B reviewed the same diff in parallel).
**Target:** `fix/388-paired-report-preflight` off `dev` @ `c4599f5`, uncommitted.
**Diff:** `src/main.rs` (two blocks), `tests/integration_output_collision.rs` (T1–T6), `CHANGELOG.md`.
**Nothing was fixed.** The tree is shared and uncommitted — every item below is a recommendation only.
**Repro binaries:** patched = `target/release/trim_galore` (verified to contain the fix);
baseline = scratch build of `HEAD` (= pre-change) at
`…/scratchpad/crev388a/btarget/release/trim_galore`, built from `git archive HEAD` (no repo writes).

## Verdict

The change is correct and does what it claims. Candidate expressions match the writers exactly at
both sites, the gate matches, no over-rejection was found in a five-shape sweep, and all four
rejection tests were confirmed discriminating against a real pre-change binary. Findings are one
scope/completeness item (a live sibling of the same defect on another paired path), and test/
diagnostic strengthening.

## Verified correct

- **Candidates == writers, FASTQ path.** New block `src/main.rs:765-772` pushes
  `report_name`/`json_report_name` for `chunk[0]`/`chunk[1]` with `output_dir`; the writer builds the
  identical four paths at `src/main.rs:1578-1585` from the same `output_dir` passed into `run_paired`,
  under the same `!cli.no_report_file` gate (`:1566`). Report names ignore `--basename`, so no
  discriminator is missed.
- **Candidates == writers, uBAM path.** New block `src/main.rs:1867-1873` vs writer
  `src/main.rs:2178-2185`, gate `:2164`. Same namers, same `output_dir`.
- **All four rejections fire on the report path, not a primary.** Each error message names
  `…_trimming_report.txt`, so the new candidates are what trip: T1/T4 `out/reads.fq_trimming_report.txt`,
  T2 `out/r1.fq… and out/R1.fq…`, T3 `out/r2.fq… and out/R2.fq…`.
- **Single-file interleaved paths correctly out of scope.** `run_paired_ubam_single_file`
  (`src/main.rs:1794-1803`) and `run_ubam_output_paired_single_file` write stem-derived
  `_R1`/`_R2` report names and take exactly one input, so neither a cross-input nor a self collision
  is constructible; they also have no pre-flight at all, consistently.
- **Specialty paired modes write no trimming reports** (`grep` of `src/specialty.rs`) — no hole there.
- **A4 (uBAM SE residual) is sound.** For SE, a report collision implies fold-equal filenames, which
  implies fold-equal stems, which collides `_trimmed.bam` first. No report-only hole exists there.
- **`guarded_inputs` / output-vs-input interaction is clean**, and strictly an improvement — see L3.
- **Efficiency:** four extra `PathBuf`s per pair into a list the pre-flight already hashes. Nothing.

## Negative controls (reproduced independently, all four)

Run against the baseline binary, i.e. the plan's §8 claim re-derived rather than trusted:

| shape | baseline | patched |
|---|---|---|
| T1 `--paired -o out A/reads.fq B/reads.fq` | exit 0, 2 `_val_` + **one** report pair | exit 1, 0 files |
| T2 `--paired -o out a/r1.fq b/R1.fq` | exit 0, surviving report says `Input filename: R1.fq` | exit 1 |
| T3 4-file cross-pair | exit 0, **3** reports for 4 inputs | exit 1 |
| T4 same as T1 + `--output-format ubam` | exit 0, `reads_val.bam` + one report pair | exit 1 |

T5/T6 pass on both, as guards should.

## High

### H1 — the identical defect is still live on `--clump_only --paired` (FASTQ)

Pre-existing, not introduced by this diff, but it is the same signature the fix exists to close:
two per-side reports whose names lack the positional discriminator the primaries carry. Verified
against the **patched** binary:

```
cd <scratch>/c1   # A/reads.fq and B/reads.fq, both plain FASTQ
trim_galore --clump_only --paired -o out A/reads.fq B/reads.fq   # exit 0
out/reads_clumped_1.fq
out/reads_clumped_2.fq
out/reads.fq_clumping_report.txt      # "Input: B/reads.fq" — R1's report is gone
```

`--clump_only --paired` routes through `run_specialty_paired` with `clumped_paired_output_names` as
the namer (`src/main.rs:521-544`), so only the two primaries are hashed, while
`clump_only::clump_only_paired` writes two reports (`src/clump_only.rs:535-536`).

Not affected, checked: the uBAM clump paired path writes **one** report per pair derived from R1
(`src/clump_only.rs:1082-1086`), a key finer than its BAM primary's, so the existing primary check
covers it.

Recommendation: either extend this branch to that path, or add it to §5 next to A4 as a named
residual with the reproducer. The CHANGELOG entry is literally accurate (clumping reports are a
distinct artifact), so no CHANGELOG change is needed if it stays out of scope — but leaving it
unnamed invites the same blind spot #388 came from.

## Medium

### M1 — the diagnostic is uninformative for the flagship shape

T1's own error:

```
Output path collision (…): out/reads.fq_trimming_report.txt and out/reads.fq_trimming_report.txt
would be written to the same file. Check that inputs produce distinct output paths
(e.g., different source directories or `--output_dir`).
```

The path is printed twice, neither colliding **input** is named, and the advice recommends the two
things that cannot help: the inputs already are in different source directories, and `--output_dir`
is what causes the collision. `preflight_output_collisions` takes a `hint` (used by the specialty
modes via `CWD_OUTPUT_HINT`); both new call sites pass `None`. Partly pre-existing for primary
collisions, but T1 and `--basename foo -o out A/reads.fq B/reads.fq` (also newly refused, verified)
now fail *only* for this reason, so it is the first thing a user will hit.
Recommendation: pass a hint on the two paired sites naming the actual remedy (rename one input, or
drop `--output_dir`), or have the message name the two inputs when the two planned paths are equal.

### M2 — T2/T3/T4 bypass the file's own rejection helper

`assert_rejected_cleanly` (`tests/integration_output_collision.rs:93`) asserts the `PREFIX`, the
expected wording, **and** that nothing was written — the module doc calls two reports beside one data
file "#383's most misleading artifact". T1 hand-rolls a weaker version of it; T2/T3/T4 assert only
exit status and `DUP_MSG`, so they would pass a refusal that had already written files. Recommendation:
`assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG)` in all four. (Empirically 0 leftovers today, so
this is future-proofing, not a live gap.)

### M3 — T6's acceptance cases don't assert what they claim to guard

The module doc states acceptance cases "assert content or filename rather than mere existence".
T6's second case asserts only `ok` — no filename at all; the first asserts only `foo_val_1.fq`.
Both comments claim a property about *report* names that neither test checks. Recommendation: assert
`r1.fq_trimming_report.txt` + `r2.fq_trimming_report.txt` in the basename case, and
`sample.fq_trimming_report.txt` + `sample.fastq_trimming_report.txt` in the same-stem case. That is
exactly the "report keys differ (full filename)" claim, and it would fail if `report_name` ever
became stem-based.

### M4 — the nearest true negative to T1 is untested

T1 minus `-o`: same filename as R1/R2 in different directories, no `--output_dir`. Verified accepted
(exit 0; reports land beside their own inputs in `A/` and `B/`, primaries both in `A/`). This is the
only shape that pins the new candidates honouring `output_dir = None`; an implementation that
resolved reports into one directory — as `run_paired_ubam_single_file` legitimately does with its own
`dir` — would over-reject a valid run and no test would notice. Recommendation: add it to T6.

### M5 — the premise #388 disproved is still stated unqualified, next to the fix

Three comments assert the SE argument in general terms:

- `src/io.rs:1211` "Assumption A2, in the direction that makes the pre-flight sufficient: where two
  inputs' PRIMARY output paths differ, every secondary output path must differ too … hashing
  primaries alone is enough."
- `src/io.rs:1285-1287` "the primary key is strictly *coarser* than the report key, which is why
  checking primaries covers reports".
- `src/main.rs:79` "assumption A2 argues that checking *primary* paths covers secondary paths".

All three are true only for single-end naming; `_val_1`/`_val_2` is the counterexample this branch
exists for. Recommendation: one line each, e.g. `(single-end naming only; paired `_val_N` primaries
break this — #388)`. Cheap, and it is the reasoning trap that produced the bug.

## Low

- **L1 — duplication.** The two new blocks are near-identical, and `planned_secondary_outputs`
  (`src/main.rs:83`) already computes the same pair for SE, so three sites must now stay in step with
  the writers. A `report_names(input, output_dir) -> [PathBuf; 2]` in `src/io.rs` beside `report_name`
  / `json_report_name` would give one. Acceptable as-is under the house pattern of explicit
  per-dispatch-path pre-flights; noted, not urged.
- **L2 — T5's `read_dir().count() == 2`** is correct and robust (verified: plain input ⇒ plain output,
  no FastQC, no extra artifacts), but on failure it can only say `2 != 3`. `assert_dir_holds_only`
  already exists in the file and names the intruder.
- **L3 — the CHANGELOG undersells the fix by one repair.** The paired report path was also missing
  from the output-vs-**input** check, so a report could overwrite a named input. Verified:
  `--paired r1.fq r1.fq_trimming_report.txt` (FASTQ content in the second file) on the baseline exits
  0, writes `r1.fq_trimming_report_val_2.fq`, then **overwrites its own R2 input** with R1's report —
  that input's reads are gone at exit 0. Patched refuses with the alias message. This makes the
  earlier unreleased entry's claim ("This affected `--paired` as well, and is now refused on every
  path") true, which it was not before this diff. Worth one sentence in the #388 entry.
- **L4 — placement.** The two newly-refused shapes are behaviour changes sitting under
  `#### Bug fixes` while the file has a `#### Changes` section. Consistent with the #383/#385 entry
  above it, so arguably fine; flagging for the author's call.

## Not findings

CHANGELOG claims checked line by line against the baseline and all hold, including "one report pair
for two inputs" and "re-using a filename across pairs lost a report the same way". Plan §8's
"no deviations" is accurate: the committed blocks are byte-identical to §4's listing.
