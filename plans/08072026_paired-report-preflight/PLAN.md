# Plan — Add report paths to the paired pre-flight (#388)

**Issue:** [#388](https://github.com/FelixKrueger/TrimGalore/issues/388)
**Branch point:** `dev` @ `c4599f5`
**Revision:** v2 (2026-08-07) — see §7

---

## 1. Goal

The `--paired` trim pre-flight checks `_val_`/`_unpaired_`/passthrough paths but not the
four report paths per pair, so a trimming report can silently overwrite another report at
exit 0. Add `report_name`/`json_report_name` for both sides of every pair to the candidate
list — on the paired **FASTQ** path and the paired **uBAM-output** path, which has the
identical live hole (verified) — gated on `!cli.no_report_file` to match the writers.

## 2. The defect, stated correctly

**Root cause (Reviewer B's formulation): paired primaries carry a positional discriminator
(`_val_1`/`_val_2`) that report names do not.** #385's SE argument — report keys are finer
than primary keys, so checking primaries suffices — held because SE has one primary per
input. In paired mode two inputs can have distinct primaries and colliding reports. v1
inherited the SE argument invalidly; every broken piece of v1 flowed from that.

Three verified reproducers on `dev` @ `c4599f5`, all exit 0 with one report pair lost:

1. **Within one pair, same filename, no case involved** (Reviewer B):
   `--paired -o out A/reads.fq B/reads.fq` → `reads_val_1.fq`, `reads_val_2.fq`, ONE
   report pair.
2. **Within one pair, fold-equal filenames** (Reviewer A): `--paired -o out a/r1.fq
   b/R1.fq` → surviving report says `Input filename: R1.fq`; R1's report gone.
3. **Cross-pair** (Reviewer A): `--paired a/r1.fq a/r2.fq b/R2.fq b/x.fq -o out` →
   3 report files for 4 inputs.

The APFS self-pair case from the issue (`r1.fq` + `R1.fq` naming one file) is **already
refused** by `Cli::validate`'s absolutising R1≠R2 check for same-spelling variants; the
fix additionally covers the genuinely-two-files variants above, which validate correctly
passes.

## 3. Behaviour changes (deliberate, need sign-off)

The fix converts all three reproducers from silent report loss into pre-flight refusals.
Two of them are *currently-working invocations* in the sense of exiting 0:

- **Same-filename R1/R2 into a shared output dir** (reproducer 1) — newly refused on
  every filesystem. This matches the SE precedent already committed
  (`se_trim_rejects_same_basename_across_dirs_with_output_dir`), and today's exit-0 run
  loses a report, so refusal is the correct outcome — but it is a visible change and
  gets its own CHANGELOG line.
- **Case-differing distinct inputs on a case-sensitive filesystem** (reproducer 2 run on
  ext4) — a true false positive there, since both reports could coexist. This is the
  documented #216 trade (`io.rs:45-47`: loud error over silent loss on APFS/NTFS), now
  extended to paired reports; stated, not denied.

## 4. Implementation outline

**Step 1 — paired FASTQ site.** `src/main.rs`, the `if cli.paired` FASTQ branch
(`:730-767`; note the search token `Pre-flight across pairs before any I/O` is ambiguous —
it also matches `run_specialty_paired` at `:2458`, which has no `candidates` binding and
would fail to compile). Inside the per-chunk loop after the passthrough block:

```rust
// #388 — report names carry no _val_ discriminator, so two inputs with
// distinct primaries can still collide on reports.
if !cli.no_report_file {
    for input in [&chunk[0], &chunk[1]] {
        candidates.push(naming::report_name(input, output_dir));
        candidates.push(naming::json_report_name(input, output_dir));
    }
}
```

**Step 2 — paired uBAM-output site.** `run_ubam_output`'s paired pre-flight loop
(`:1851-1859`), which pushes to `planned` (not `candidates`); writers at `:2163-2169` are
gated on the same flag. Same four paths per chunk, pushed to `planned`, same gate.

**Step 3 — tests**, `tests/integration_output_collision.rs`:

- **T1 (the fix's raison d'être):** reproducer 1 → non-zero exit, collision message,
  nothing written. Case-free, runs identically everywhere.
- **T2 (fold-equal):** reproducer 2 → refused. On case-insensitive FS this is a true
  positive, on case-sensitive a documented false positive; either way the *rejection* is
  asserted, so the test is filesystem-independent.
- **T3 (cross-pair):** reproducer 3 → refused.
- **T4 (uBAM twin):** reproducer 1 with `--output-format ubam` → refused.
- **T5 (gate):** reproducer 1 plus `--no_report_file` → **succeeds**, two `_val_` files,
  no reports (verified achievable today: validate passes these inputs). Pins that the
  candidates cannot drift ahead of the writers.
- **T6 (over-rejection guards, at the boundary):** single pair + `--basename` without
  `-o` → succeeds; same-stem-different-extension R1/R2 (`sample.fq` + `sample.fastq`)
  → succeeds. (Ordinary distinct pairs are already covered by
  `paired_accepts_two_distinct_pairs`.)

**Step 4 — CHANGELOG**, `#### Bug fixes`: the report-clobbering fix; plus the two §3
behaviour-change sentences.

## 5. Assumptions

- **A1 (corrected).** The newly-rejected set is exactly: *pairs (or pair-sets) whose
  report keys collide while primaries do not* — same-filename R1/R2 into one dir,
  fold-equal filenames into one dir, and cross-pair repeats. Reviewer-verified sweep of
  adjacent shapes found no other over-rejection: multi-pair repeated basenames + `-o`
  already collide on primaries; `--basename` >1 pair already rejected by validate;
  same-basename different dirs without `-o` keeps reports beside their inputs.
- **A2.** Candidate gating must equal writer gating (`!cli.no_report_file`, both sites).
  T5 pins it.
- **A3.** FastQC artifacts stay out of the candidate list (unchanged #385 rationale).
- **A4 (scope).** The uBAM **SE** path lacks `planned_secondary_outputs` too, but its
  primaries (`_trimmed.bam`) fold and collide first (reviewer-verified exit 1), so only
  a C2-class residual remains there; noted, not fixed here.
- **A4b (residual, found by code review).** `--clump_only --paired` has the identical
  defect — `_clumped_1`/`_clumped_2` primaries vs discriminator-free
  `clumping_report_name` — reproduced on this branch and filed as
  [#391](https://github.com/FelixKrueger/TrimGalore/issues/391). Not fixed here: the fix
  restructures `run_specialty_paired`'s signature, which is scope growth at review time.

## 6. Validation

T1–T6; gates (`fmt`, `clippy -D warnings`, full `cargo test`); negative control: comment
the Step-1 block → T1, T2, T3 must fail while T5, T6 still pass; comment the Step-2 block
→ T4 must fail. The §2 reproducers re-run post-build must all refuse.

## 7. Revision history

**v1 → v2** after dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`) — the two
reviews converged on near-identical findings, all verified by the orchestrator:

- All four v1 tests were broken: V1/V4 could not pass (validate's absolutising R1≠R2
  check fires first — a consequence of my own #385 D8/D13 change), V2 could not fail
  (the format detector rejects report-named inputs before the pre-flight; a genuine
  report is never valid FASTQ, so v1's "prior report fed back in" benefit was also
  near-unreachable), V3 could not detect over-rejection.
- v1's A1 ("report keys strictly finer, no new false positive possible") was false in
  both halves; B's root-cause formulation replaces it (§2), and the behaviour changes
  are now §3 items needing sign-off instead of being denied.
- The paired-uBAM twin came in scope: both reviewers independently verified the
  identical live bug there, and v1's exclusion reason ("mirrors the SE precedent") was
  not a reason.
- v1's insertion token was ambiguous (two hits); §4 pins the site by line range and
  notes the second site's compile-failure property.
- v1's §1 claim that "the SE path" gained secondary-output protection in #385 narrowed
  to "the SE **FASTQ** path" (A's I4).

---

## 8. Implementation notes

Implemented on `fix/388-paired-report-preflight` off `dev` @ `c4599f5`. All four §4 steps
as specified in v2, no deviations. Diff: `src/main.rs` (two blocks, 18 lines),
`tests/integration_output_collision.rs` (T1–T6), `CHANGELOG.md` (one entry naming both
§3 behaviour changes).

Gates: **546 tests** (was 540), fmt clean, `clippy -D warnings` clean. All three §2
reproducers refuse against the rebuilt binary with zero files written (the residue
counter first reported 2 — it was counting `ls`'s blank separator lines; `find -type f`
confirms 0).

Negative controls per §6, both discriminating: FASTQ block disabled → T1/T2/T3 FAIL
while T4/T5/T6 pass; uBAM block disabled → T4 FAILs alone. Reverted to green both times.

### Post-code-review round

Coverage: **COMPLETE, 20/20, first pass** — no gaps. Code review: A 0 Critical / 1 High,
B 0 Critical / 0 High; the High/Medium finding is the same on both sides
(`--clump_only --paired`, verified, filed as #391 + A4b above). Applied in-round:

- **A-M1:** both paired sites now pass `PAIRED_REPORT_HINT` — the generic advice
  ("different source directories or `--output_dir`") recommends the two things that
  cannot help for filename-keyed collisions. Same replace-don't-append mechanism as #384.
- **A-M2 / B-L1:** T2/T3/T4 route through `assert_rejected_cleanly` (nothing-written now
  asserted everywhere).
- **A-M3 / A-M4:** T6 asserts the report filenames it claims are distinct, and gains the
  T1-minus-`-o` acceptance case pinning `output_dir = None` handling.
- **A-M5:** the three comments still stating the single-end A2 premise unqualified are
  now scoped to SE naming with a #388 pointer — the reasoning trap that produced the bug.
- **B-M2 (declined):** the two candidate blocks stay near-duplicates; extracting a helper
  across `candidates`/`planned` buys little at two sites and the house precedent (eleven
  call sites before #385 consolidated them for a functional reason) tolerates it.

One iteration: the first test-patch script asserted against pre-rustfmt text and wrote
nothing (the assert fired before the save — caught because the test count didn't move);
re-applied against the on-disk shape with a count-verified regex.

Final gates: 546 tests, fmt clean, clippy clean.
