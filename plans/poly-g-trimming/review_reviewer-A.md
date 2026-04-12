# Plan Review: Poly-G Trimming with Auto-Detection

**VERDICT: APPROVE WITH SUGGESTIONS**

Overall this is an excellent, well-structured plan that follows the existing codebase patterns closely. The auto-detection algorithm is sound, the integration points are correct, and the design decisions are well-reasoned. The suggestions below range from one important correctness issue to several minor improvements.

---

## Critical

*None.*

## High

### H1. Auto-detection bypassed when adapter is user-specified or preset-selected

**Location:** Plan section 4 (`main.rs`), existing code in `main.rs` `resolve_adapter()`

**Issue:** The auto-detection scan (`autodetect_adapter()`) is only called when no adapter is explicitly specified (the `else` branch at the bottom of `resolve_adapter()`). When a user passes `--illumina`, `--nextera`, `--adapter ACGT`, or any preset, `autodetect_adapter()` is never called, which means `detection.poly_g_count` is never populated. The plan's `main.rs` wiring code assumes `detection` always exists with a valid `poly_g_count`.

This means poly-G auto-detection silently fails for any user who specifies an adapter explicitly -- which includes most production pipelines (e.g., `--illumina`). This is a significant gap because the users most likely to specify an adapter are also likely running on NovaSeq/NextSeq.

**Suggestion:** Either:
- (a) Always run the poly-G scan (even when adapter is user-specified) by extracting it into a separate function, or by refactoring `resolve_adapter()` to always scan and return the poly-G count. This is the cleanest approach.
- (b) Alternatively, add a separate lightweight scan function `detect_poly_g()` that only counts poly-G reads (no adapter counting), and call it whenever auto-detection is relevant (i.e., when `!cli.poly_g && !cli.no_poly_g`). This avoids redundant adapter scanning but adds a second pass over the file for the explicit-adapter case.

Option (a) is simpler and consistent -- the scan is fast (already bounded to 1M reads), and adapter counting is negligible overhead on top of it.

### H2. `--poly_g` and `--no_poly_g` need validation when combined with `--nextseq`

**Location:** Plan section 3 (`cli.rs`), section 4 (`main.rs`)

**Issue:** The plan says `--nextseq` and auto-detected poly-G can be active simultaneously and are "complementary, not redundant." This is true. However, the plan should explicitly document in the CLI help text that `--no_poly_g` does NOT disable `--nextseq`, and that `--nextseq` does NOT imply `--poly_g`. Users may be confused about the distinction between quality-based G-trimming (nextseq) and homopolymer-based G-trimming (poly_g). A note in the `--poly_g` help text saying "This is independent from --nextseq" would prevent user confusion.

**Suggestion:** Add a brief clarification to the CLI help strings. No code change needed beyond documentation.

---

## Medium

### M1. `homopolymer_trim_index()` signature: `revcomp` naming is misleading for poly-G

**Location:** Plan section 1 (`quality.rs`)

**Issue:** For poly-A, `revcomp=true` means "scan for poly-T on R2" because the poly-A tail on the sense strand appears as poly-T on the antisense strand. The `revcomp` parameter name makes biological sense for poly-A/poly-T.

For poly-G, the situation is different: poly-G artifacts are chemistry artifacts, not biology. Both R1 and R2 can have poly-G tails on the _same_ strand they're sequenced on (the sequencer produces G for no-signal regardless of strand). The plan says R2 gets poly-C trimming from the 5' end (reverse complement logic), which mirrors the poly-A treatment. This is correct for paired-end data where R2 is the reverse complement of R1 -- a poly-G tail on the fragment appears as poly-C at the start of R2.

However, this only holds when R2 reads the opposite strand. For some library preps (e.g., single-end data that was converted to "fake paired-end"), the assumption may not hold. The plan correctly handles R1/SE (poly-G from 3') and R2 (poly-C from 5'), which matches the standard Illumina paired-end chemistry, so this is fine in practice. But the `revcomp` naming could be clarified in the doc comment.

**Suggestion:** In the `homopolymer_trim_index()` doc comment, note that for poly-G the "revcomp" behavior trims poly-C from the 5' end (for R2 in standard paired-end). Consider naming the parameter `trim_5prime` instead of `revcomp`, since that's what it actually controls -- but this would be a broader refactor that also affects poly-A, so it may not be worth the churn.

### M2. Threshold calculation has an edge case with integer division for very small files

**Location:** Plan section 4 (`main.rs`)

```rust
let threshold = detection.reads_scanned / 10_000; // 0.01%
let threshold = threshold.max(10); // at least 10 reads
```

**Issue:** For files with fewer than 10,000 reads, `reads_scanned / 10_000` is 0, and then `.max(10)` makes the threshold 10. This means:
- 100 reads: need >10 poly-G reads to trigger (>10%) -- very conservative but OK
- 1,000 reads: need >10 poly-G reads to trigger (>1%) -- reasonable
- 10,000 reads: need >1 poly-G read... wait, 10_000/10_000 = 1, max(10) = 10. Need >10 reads. OK.
- 100,000 reads: need >10 poly-G reads to trigger (>0.01%) -- correct
- 1,000,000 reads: need >100 to trigger -- correct

The logic is actually fine. The plan acknowledges the "very short files" edge case already. No issue here on reflection.

### M3. Report output: poly-G stats should print even when count is 0 if poly-G was enabled

**Location:** Plan section 6 (`report.rs`)

**Issue:** The plan says `write_report_stats()` prints poly-G stats only `when stats.poly_g_trimmed > 0`. But if a user explicitly enables `--poly_g` and no reads are trimmed (e.g., data has no poly-G tails), the report should still mention that poly-G trimming was active. This matches the header which says "Poly-G trimming enabled" but the stats section would be silent.

**Suggestion:** Consider printing the poly-G stats line whenever `config.poly_g` is true (even if count is 0), to confirm that the feature ran. This is consistent with how adapter stats always print even when `reads_with_adapter == 0`.

### M4. Missing: console summary in `main.rs` for poly-G stats

**Location:** `main.rs`, `run_single()` and `run_paired()`

**Issue:** The plan covers adding poly-G stats to `report.rs` (the file-based report) but does not mention updating the on-screen summary in `main.rs`. Currently `main.rs` has `eprintln!` blocks for poly-A stats (lines 232-235 in `run_single()`, lines 373-379 in `run_paired()`). Analogous blocks are needed for poly-G.

**Suggestion:** Add the same pattern used for poly-A console output:
```rust
if stats.poly_g_trimmed > 0 {
    eprintln!("Reads with poly-G/C trimmed:     {:>10} ({:.1}%)", ...);
    eprintln!("  Poly-G/C bases removed:        {:>10}", ...);
}
```

---

## Low / Minor

### L1. `has_trailing_poly_g()` could be slightly more efficient with early exit

**Location:** Plan section 2 (`adapter.rs`)

**Issue:** The helper counts consecutive G's from the 3' end. It returns `true` as soon as `count >= min_len`, which is correct and efficient. However, `seq.iter().rev()` will visit every base until a non-G is found. For reads that are entirely poly-G (rare edge case), this scans the full read. This is fine performance-wise (bounded by read length, typically 150bp), just noting it for completeness.

### L2. Test suggestion: add a test for the combined `--poly_a --poly_g` case in `trimmer.rs`

**Location:** Plan section "Testing" item 6

**Issue:** The plan mentions testing `--poly_a --poly_g` together but only as an integration test. It would also be valuable to have a unit test in `trimmer.rs` (or `quality.rs`) for a sequence like `ACGTACGTAAAAGGGGG` to verify that poly-A removes the A's first, then poly-G removes the G's, and the final result is `ACGTACGT`. This confirms the ordering assumption from edge case analysis.

### L3. Detection message could include the raw count, not just percentage

**Location:** Plan section 4 (`main.rs`)

**Issue:** The auto-detection messages show the percentage but not the absolute count. For debugging and reproducibility, showing both would be helpful: "2,345 of 1,000,000 reads (0.23%) have poly-G tails (>=10bp)".

### L4. Plan does not mention updating `report::TrimConfig` for poly-G

**Location:** Plan section 6 (`report.rs`)

**Issue:** The plan mentions adding `pub poly_g: bool` to `ReportConfig` but the actual struct is named `report::TrimConfig` (line 91 of `report.rs`). The plan uses "ReportConfig" which doesn't exist. This is just a naming mismatch in the plan -- the implementation should use the actual struct name `report::TrimConfig`.

### L5. No mention of updating `Cargo.toml` or dependencies

The plan correctly identifies that no new dependencies are needed -- the implementation reuses existing functions and patterns. Good.

---

## Algorithm Assessment

### Auto-detection algorithm

The detection algorithm is sound:
- **Statistical basis:** P(10+ consecutive trailing G's by chance) = (1/4)^10 ~ 10^-6. In 1M reads, expected ~1. Threshold at 100 gives ~100x safety margin. On 2-colour instruments, typical rates are 0.5-5%, giving 5,000-50,000 hits. The signal-to-noise ratio is enormous (~10,000x above threshold).
- **Data-driven approach** is strictly better than read-name matching (fastp's approach), which is brittle against SRA renaming, platform evolution, and custom instrument names.
- **Detection strictness vs. trimming sensitivity** is a clever design: strict detection (10 consecutive, no mismatches) minimizes false positives; once detection triggers, permissive trimming (3+ bp with mismatches via Cutadapt scoring) catches even partial tails. This two-tier approach is well-reasoned.

### Trimming algorithm (Cutadapt scoring)

The +1/-2 scoring with 20% max error rate is appropriate for poly-G:
- Poly-G tails from 2-colour chemistry are typically very pure (the instrument reports no-signal as G with high confidence), so the 20% error tolerance is generous enough to handle rare basecall errors within the tail.
- The min 3bp threshold prevents false positives from random trailing G's while catching even short tails.
- Reusing the same algorithm as poly-A is a good choice for code consistency and maintainability.
- There is no biological reason poly-G would need different scoring parameters than poly-A (both are homopolymer tails; poly-G is actually *more* pure since it's a chemistry artifact, not biological with potential degradation).

### Pipeline position

Step 2.8 (after poly-A at 2.7, after adapter and RRBS) is correct. Adapter removal exposes poly-G tails that were hidden behind the adapter sequence; running poly-G after adapter trimming catches these.

---

## Summary

The plan is thorough, well-designed, and follows existing patterns closely. The one high-priority issue (H1: auto-detection bypass when adapter is explicitly specified) needs to be addressed before implementation, as it would silently disable auto-detection for most production users. The remaining suggestions are improvements, not blockers.
