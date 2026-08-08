# Plan Coverage Report

**Mode:** B (code vs. plan §4 outline — no `IMPL.md` by design)
**Plan(s):** `plans/08072026_paired-report-preflight/PLAN.md` (v2)
**Date:** 2026-08-08
**Verdict:** COMPLETE

Audited: branch `fix/388-paired-report-preflight` off `dev` @ `c4599f5`, uncommitted working
diff (`CHANGELOG.md`, `src/main.rs`, `tests/integration_output_collision.rs`; +183 lines).

## Summary

- Total items: 20
- DONE: 20
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 0

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Report candidates added to paired **FASTQ** pre-flight, inside the per-chunk loop after the passthrough block, gated `!cli.no_report_file`, pushed to `candidates` | §4 Step 1 | DONE | `src/main.rs:765-772`, inside `if cli.paired` (`:730`), after the passthrough block (`:757-764`), before `planned.extend(candidates)` (`:773`). Code matches the plan's snippet verbatim including the comment. |
| 2 | Report candidates added to paired **uBAM-output** pre-flight, pushed to `planned`, same gate | §4 Step 2 | DONE | `src/main.rs:1867-1873`, in `run_ubam_output`'s `if cli.paired` loop (`:1852-1874`) after `paired_bam_output_name`. Pushes to `planned`, four paths per chunk, same gate. |
| 3 | T1 — reproducer 1 (same filename R1+R2, shared `-o`) refused, collision message, nothing written | §4 Step 3 | DONE | `paired_rejects_same_filename_r1_r2_into_shared_output_dir` (`tests/integration_output_collision.rs:860`). Asserts `!ok`, `DUP_MSG`, and empty output dir. Case-free. |
| 4 | T2 — reproducer 2 (fold-equal filenames as R1/R2) refused | §4 Step 3 | DONE | `paired_rejects_fold_equal_filenames_into_shared_output_dir` (`:883`). Asserts rejection only, so FS-independent as specified. |
| 5 | T3 — reproducer 3 (cross-pair filename reuse) refused | §4 Step 3 | DONE | `paired_rejects_cross_pair_report_collision` (`:899`). Exact 4-input invocation from §2. |
| 6 | T4 — uBAM twin: reproducer 1 with `--output-format ubam` refused | §4 Step 3 | DONE | `paired_ubam_rejects_same_filename_r1_r2_into_shared_output_dir` (`:922`). |
| 7 | T5 — gate: reproducer 1 + `--no_report_file` succeeds, two `_val_` files, no reports | §4 Step 3 | DONE | `paired_same_filename_accepted_when_reports_disabled` (`:950`). Asserts `ok`, both `reads_val_{1,2}.fq`, and `read_dir(out).count() == 2`. |
| 8 | T6 — over-rejection guards: `--basename` without `-o`; same-stem-different-extension R1/R2 | §4 Step 3 | DONE | `paired_report_candidates_do_not_over_reject` (`:979`). Both sub-cases present as specified (`r1.fq`+`r2.fq` with `--basename foo`; `sample.fq`+`sample.fastq`). |
| 9 | CHANGELOG entry for the report-clobbering fix under `#### Bug fixes` | §4 Step 4 | DONE | `CHANGELOG.md:112-128`, inside `#### Bug fixes` (heading at `:6`), immediately before `#### Changes` (`:129`). Names #388, the `_val_` discriminator root cause, both collision routes, and the `--no_report_file` gate. |
| 10 | CHANGELOG sentence for §3 behaviour change 1 (same-filename R1/R2 into shared `-o` newly refused on every filesystem) | §4 Step 4 / §3 | DONE | `CHANGELOG.md:119-121`: "A pair whose R1 and R2 share a filename written to a common output directory is now refused on every filesystem — that run previously lost a report." |
| 11 | CHANGELOG sentence for §3 behaviour change 2 (case-differing distinct inputs refused on case-sensitive FS; the #216 trade) | §4 Step 4 / §3 | DONE | `CHANGELOG.md:121-124`: names the case-sensitive false positive explicitly and attributes it to the #216 loud-error-over-silent-loss trade. |
| 12 | §6: T1–T6 all exist as named | §6 | DONE | Six new `#[test]` fns, one per T, at `:860/:883/:899/:922/:950/:979`. No T merged or dropped. |
| 13 | §6/§8: negative-control claims match what the controls can show | §6, §8 | DONE | Verified empirically without mutating the source — see "Negative-control verification" below. Both controls discriminate exactly as §8 states. |
| 14 | §6: the three §2 reproducers are covered by tests **and** refuse against the built binary | §6 | DONE | T1/T2/T3 are the three reproducers verbatim. Independently re-run against `target/debug/trim_galore`: all three exit 1, emit the collision message, and leave `find out -type f` = 0. Confirms §8's zero-residue claim. |
| 15 | §6 gates: `fmt`, `clippy -D warnings`, full `cargo test` | §6 | DONE | `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean; `cargo test` = **546 passed, 0 failed** across 14 binaries, matching §8's 546. |
| 16 | A1 — newly-rejected set is exactly the colliding-report shapes; no over-rejection of adjacent shapes | §5 A1 | DONE | Traced to tests: T6 (two boundary shapes) plus pre-existing `paired_accepts_two_distinct_pairs` (`:555`, two pairs into a shared `-o`) which still passes. A1's remaining adjacent shapes are carried as documented reviewer-verified rationale, as the plan states. |
| 17 | A2 — candidate gating equals writer gating at both sites | §5 A2 | DONE | FASTQ writers: `run_paired` `if !cli.no_report_file` (`:1566`) → `report_name`/`json_report_name` for `input_r1` and `input_r2` (`:1578-1584`). uBAM writers: `run_ubam_output_paired_two_files` `if !cli.no_report_file` (`:2164`) → same four calls (`:2178-2184`). Identical flag, identical four paths. T5 pins it. |
| 18 | A3 — FastQC artifacts stay out of the candidate list | §5 A3 | DONE | Diff adds no FastQC paths at either site; the #385 rationale is unchanged. Verified by inspection of both hunks. |
| 19 | A4 — uBAM **SE** path deliberately out of scope | §5 A4 | DONE | Scope honoured: the diff touches only the two paired sites. `run_ubam_output_single` (`:1934`) is untouched. Documented as a noted residual, not a gap. |
| 20 | §8 claims "all four §4 steps as specified, no deviations" — verify against the diff | §8 | DONE | No code deviation found: both blocks match §4's specified site, container (`candidates` vs `planned`), path set and gate. See the one non-blocking prose inaccuracy below. |

## Gaps (detail)

None. No PARTIAL, MISSING or DEVIATED items.

## Negative-control verification

§6 prescribes source mutation; §8 reports both controls discriminating. I verified the same
property **without editing the repo**, by using `--no_report_file` to switch off exactly the
candidates each block contributes — behaviourally equivalent to removing the block:

| Reproducer shape | With reports (as shipped) | Reports disabled | Conclusion |
|---|---|---|---|
| 1 — `A/reads.fq` + `B/reads.fq`, `-o out` | exit 1, collision msg, 0 files | exit 0 → `r1`/`reads_val_1.fq`, `reads_val_2.fq` | primaries do not collide; Step 1 is the sole cause of rejection |
| 2 — `a/r1.fq` + `b/R1.fq`, `-o out` | exit 1, collision msg, 0 files | exit 0 → `r1_val_1.fq`, `R1_val_2.fq` | same |
| 3 — cross-pair, `-o out` | exit 1, collision msg, 0 files | exit 0 → 4 distinct `_val_` files | same |
| 4 — reproducer 1 + `--output-format ubam` | exit 1 (T4) | exit 0 → single `reads_val.bam` | `paired_bam_output_name` yields ONE primary per pair, so Step 2 is the sole cause |

Corroborated by the naming code: `io::norm_path` (`src/io.rs:48-50`) lowercases, so keys are
case-folded; `paired_end_output_names` (`:240-276`) appends per-input `_val_1`/`_val_2`, which
keeps all primaries distinct in every shape above, while `report_name`/`json_report_name`
(`:485`, `:501`) key on the full input filename and collide. So §8's table — "FASTQ block
disabled → T1/T2/T3 FAIL while T4/T5/T6 pass; uBAM block disabled → T4 FAILs alone" — holds:
T4 is on a different function, T5 is gated off, T6 asserts success only.

## Test verification

| Test name | File | Status |
|-----------|------|--------|
| T1 `paired_rejects_same_filename_r1_r2_into_shared_output_dir` | tests/integration_output_collision.rs:860 | PASS |
| T2 `paired_rejects_fold_equal_filenames_into_shared_output_dir` | tests/integration_output_collision.rs:883 | PASS |
| T3 `paired_rejects_cross_pair_report_collision` | tests/integration_output_collision.rs:899 | PASS |
| T4 `paired_ubam_rejects_same_filename_r1_r2_into_shared_output_dir` | tests/integration_output_collision.rs:922 | PASS |
| T5 `paired_same_filename_accepted_when_reports_disabled` | tests/integration_output_collision.rs:950 | PASS |
| T6 `paired_report_candidates_do_not_over_reject` | tests/integration_output_collision.rs:979 | PASS |
| `paired_accepts_two_distinct_pairs` (A1 baseline, pre-existing) | tests/integration_output_collision.rs:555 | PASS |
| Full suite | `cargo test` | 546 passed, 0 failed, 0 ignored |

## Non-blocking observations

- §8 describes the `src/main.rs` change as "two blocks, **18 lines**". The diff adds **15**
  (8 in the FASTQ block, 7 in the uBAM block). Prose-only inaccuracy in the notes; the code
  itself matches §4 exactly. Not counted as a deviation.

## Verdict

**COMPLETE.** All four §4 steps are implemented at the sites the plan pins, T1–T6 all exist
under the plan's names and pass, both §3 behaviour changes are stated in the CHANGELOG under
`#### Bug fixes`, all four §5 assumptions trace to a test or to documented rationale, all
three §2 reproducers refuse with zero residue, both negative controls discriminate, and the
`fmt` / `clippy -D warnings` / 546-test gates are green. Nothing outstanding.
