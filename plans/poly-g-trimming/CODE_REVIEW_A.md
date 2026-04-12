# Code Review A: Poly-G Trimming Implementation

**Reviewer:** A  
**Date:** 2026-04-12  
**Scope:** 7 files under `optimus_prime/src/` implementing poly-G trimming with auto-detection  
**Plan:** `plans/poly-g-trimming/plan.md`  
**Build status:** Compiles cleanly, all 74 tests pass  

---

## Summary

The implementation is solid, well-structured, and faithfully follows the plan. The generalisation of `poly_a_trim_index()` into `homopolymer_trim_index()` is clean, the auto-detection piggybacking is efficient, the standalone `detect_poly_g()` fallback for user-specified adapters is correctly wired, and all code paths (sequential, parallel, single-end, paired-end) are covered. No correctness bugs were found. The issues below are mostly about minor missing plan items and code hygiene.

**Verdict: PASS** -- ready to merge after addressing the Medium items below.

---

## Issues by Area

### Logic

**No issues found.** The scoring algorithm, complement mapping, scan thresholds, strand logic (R1: poly-G 3' / R2: poly-C 5'), pipeline ordering (step 2.8), and CLI override semantics are all correct.

### Efficiency

**No issues found.** Poly-G detection piggybacks on the existing adapter scan with negligible overhead (`has_trailing_poly_g` is a simple reverse iteration that short-circuits on the first non-G). The standalone `detect_poly_g()` scan only runs when strictly necessary.

### Errors / Bugs

**(None are correctness bugs -- all items below are behavioral gaps.)**

### Structure

**[S1 - Low] Clippy warning: needless_range_loop in `homopolymer_trim_index`**  
File: `quality.rs` line 135  
The forward scan `for i in 0..n { if sequence[i] == target ... }` triggers a Clippy suggestion to use `.iter().enumerate()`. However, the index `i` is also used in the error-rate check (`errors * 5 <= i + 1`), so the iterator form would not be cleaner. I'd recommend suppressing with `#[allow(clippy::needless_range_loop)]` on the function, with a comment explaining why the index is needed.

**[S2 - Low] Clippy warning: type_complexity on `resolve_adapter` return type**  
File: `main.rs` line 192  
The return type `Result<(String, String, Option<String>, Option<(usize, usize)>)>` triggers the warning. Consider a named struct (e.g., `AdapterResolution`) to improve readability and silence Clippy. This is new in this PR since the 4th tuple element was added.

---

## Plan Coverage Gaps

**[P1 - Medium] Plan amendment A5 not implemented: poly-G report stats when explicitly enabled but count is 0**  
File: `report.rs` lines 207-213  
The plan amendment A5 states: "When `config.poly_g` is true, always print the poly-G stats line in the report (even if 0 reads were trimmed)." Currently, `write_run_stats()` only prints poly-G stats when `stats.poly_g_trimmed > 0`. The function doesn't receive the config, so it has no way to know whether poly-G was explicitly enabled.

**Recommendation:** Either pass a bool/config to `write_run_stats()` to conditionally always print the line, or accept this as a conscious simplification since poly-A follows the same `> 0` pattern. If the latter, note the deviation from the plan.

**[P2 - Low] Plan called for auto-detection stats in user-override messages**  
File: `main.rs` lines 142-146  
The plan specified that `--poly_g` / `--no_poly_g` messages should include the auto-detection percentage ("Auto-detection found X.X% poly-G reads in N.") so users can sanity-check. The implementation only prints the bare enable/disable message. This is because the auto-detection data isn't available when the user short-circuits with `--poly_g` / `--no_poly_g`.

**Recommendation:** This could be addressed by always running the scan (either piggybacked or standalone) before the override check, then printing stats in all branches. However, this adds an unnecessary file scan when `--poly_g` or `--no_poly_g` is explicitly set. The current behavior is a reasonable simplification -- document the deviation from the plan.

---

## Recommendations

**[R1 - Medium] Add a test for `detect_poly_g()` standalone function**  
The `has_trailing_poly_g()` helper has thorough unit tests, but the `detect_poly_g()` function itself (which opens a file and scans it) has no integration test. Since it's a critical fallback path (all users who specify `--illumina`, `--nextera`, `--adapter`, etc. rely on it for auto-detection), it deserves at least one test with a synthetic FASTQ file verifying the count and reads_scanned are returned correctly.

**[R2 - Low] Off-by-one in scan loop (pre-existing, consistently propagated)**  
Files: `adapter.rs` lines 96-99 and 236-239  
Both `autodetect_adapter()` and `detect_poly_g()` increment `reads_scanned` then check `> MAX_SCAN_READS`, meaning they actually scan 1,000,001 reads for a file with more than 1M reads. The `detect_poly_g()` function correctly copied this pattern from `autodetect_adapter()`, so the poly-G count from both paths is directly comparable. This is harmless (0.0001% extra work, no semantic impact) but worth noting for a future cleanup.

**[R3 - Low] Missing `--poly_g` in `Cli::validate()`**  
File: `cli.rs`  
There is no validation for `--poly_g` combined with `--no_poly_g` beyond clap's `conflicts_with`. This is fine since clap handles it, but consider adding a comment noting that clap enforces the mutual exclusion.

**[R4 - Low] Console summary in paired-end: poly-G stats only print when > 0**  
File: `main.rs` lines 445-452  
Like the report (P1 above), the paired-end console summary only shows poly-G stats when `> 0`. This is consistent with poly-A behavior, but per plan amendment A2, the intent was to "Print stats when poly-G was enabled, even if count is 0 (confirms the feature ran)." Consider always printing when `config.poly_g` is true. Same applies to `run_single()` (lines 294-298).

**[R5 - Low] Unit test coverage for `homopolymer_trim_index` edge cases**  
The plan listed a test for "Poly-G with errors: `GGGGGAGGGGG`" which is tested. However, a test with a very high error rate (just above the 20% threshold) that should NOT trim would strengthen confidence. For example: `ACGTGAGGG` has 2 errors in 5 G-region bases = 40% error -- verify no trimming occurs. Also, a sequence that is entirely non-G (`ACGTACGT` with `revcomp=true`) is tested, but a 2-base sequence (`CC` with `revcomp=true`) returning 0 (short-circuit) should be verified.

---

## Fixes Applied

None. All issues are recommendations; no unambiguous low-risk fixes were identified that warrant direct editing.

---

## What Went Well

1. **Generalisation done right.** `homopolymer_trim_index()` with a target-base parameter cleanly replaces hardcoded A/T logic, and the `poly_a_trim_index()` wrapper preserves backward compatibility.

2. **Standalone detect_poly_g() fallback.** This is a critical design fix from the plan review (amendment A1) and is correctly implemented -- users specifying `--illumina`, `--adapter`, etc. still get poly-G auto-detection.

3. **Parallel path parity.** Both `process_pairs()` and `process_reads()` in `parallel.rs` correctly collect poly-G stats, matching the sequential path exactly.

4. **Thorough unit tests.** 7 new poly-G/poly-C tests in `quality.rs` and 8 new `has_trailing_poly_g()` tests in `adapter.rs` cover the key scenarios well.

5. **Clear user messaging.** The auto-detection messages include raw counts, percentages, and actionable override flags (`--no_poly_g` / `--poly_g`).
