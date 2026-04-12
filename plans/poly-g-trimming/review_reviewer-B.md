# Plan Review: Poly-G Trimming with Auto-Detection

**VERDICT: APPROVE WITH SUGGESTIONS**

Overall this is a well-designed plan that follows established patterns in the codebase closely. The auto-detection approach is scientifically sound and the implementation maps cleanly onto the existing poly-A infrastructure. The suggestions below range from one significant architectural concern to minor improvements.

---

## High Severity

### H1. Auto-detection is unavailable when adapter is user-specified or preset-selected

The `resolve_adapter()` function in `main.rs` (lines 139-185) has multiple early-return paths: when the user provides `--adapter`, `--nextera`, `--small_rna`, `--stranded_illumina`, `--bgiseq`, or `--illumina`, it returns immediately without calling `autodetect_adapter()`. This means the adapter scan loop -- which the plan proposes piggybacking on for poly-G detection -- **never runs** in these cases.

The plan's `main.rs` wiring (section 4) assumes `detection.poly_g_count` and `detection.reads_scanned` are always available, but they won't be when the user explicitly specifies an adapter.

This is a significant gap. A NovaSeq user who passes `--illumina` (or `--adapter AGATCGGAAGAGC`) would get no auto-detection and no poly-G trimming unless they also pass `--poly_g`.

**Recommendation:** Either:
- **(a)** Run a lightweight scan pass specifically for poly-G counting whenever auto-detection is skipped and `--no_poly_g` is not set. This could be a small standalone function that only counts trailing poly-G (fast, no adapter matching overhead).
- **(b)** Decouple poly-G counting from `autodetect_adapter()` entirely: add a separate `detect_poly_g(path) -> (usize, usize)` function that scans up to 1M reads. Call it from `main.rs` whenever poly-G auto-detection is needed, regardless of adapter path. This adds a second scan pass but keeps concerns separate.
- **(c)** Always run `autodetect_adapter()` even when the user specifies an adapter, just to get the poly-G count, then discard the adapter result. This would be the simplest code change but wastes ~10% of the scan doing adapter substring matching.

Option (a) seems best: fast, minimal code, no wasted work.

### H2. `poly_a_trim_index` generalisation loses information about which test exercises which path

The plan proposes keeping `poly_a_trim_index()` as a thin wrapper. This is good. But the plan should explicitly state that the existing poly-A unit tests remain unchanged and continue to call `poly_a_trim_index()` (not `homopolymer_trim_index` directly). New poly-G tests should call `homopolymer_trim_index(seq, b'G', ...)` directly to verify the generalisation. This avoids a subtle risk: if someone later removes the wrapper thinking it's redundant, the poly-A tests would break in a confusing way.

This is partly covered by the plan ("Keep `poly_a_trim_index()` as a thin wrapper for backwards compatibility with existing tests") but should be more explicit about the test strategy.

---

## Medium Severity

### M1. Off-by-one risk in the threshold calculation for small files

The plan proposes:
```rust
let threshold = detection.reads_scanned / 10_000; // 0.01%
let threshold = threshold.max(10); // at least 10 reads
```

For a file with exactly 10,000 reads: `10_000 / 10_000 = 1`, then `.max(10) = 10`. So you need >10 reads with poly-G out of 10,000, which is >0.1% -- 10x stricter than the stated 0.01% threshold. This is fine for small files (the plan acknowledges false negatives on tiny files are acceptable), but the discrepancy between the documented 0.01% and the actual threshold for files under 100K reads should be noted.

For a file with 1,000 reads: `1_000 / 10_000 = 0` (integer division), then `.max(10) = 10`. You'd need >10 poly-G reads out of 1,000 (>1%), which is very conservative. This is probably fine, but worth a comment in the code.

### M2. No mention of `--poly_g` + `--poly_a` interaction with reporting

The edge case section mentions `--poly_a --poly_g` ordering is correct. However, looking at the report output in `report.rs`, the poly-A stats say "Reads with poly-A/T trimmed". If both are active, a user might wonder whether poly-G bases are counted under poly-A stats. The plan should clarify that the report lines are completely separate and a single read can appear in both counts.

### M3. `has_trailing_poly_g` function should handle lowercase bases

The existing codebase appears to work with uppercase bases throughout (based on the patterns in `adapter.rs` and `quality.rs`). However, FASTQ files can technically contain lowercase bases. The `has_trailing_poly_g` function only checks `b'G'`. If this is intentional (matching the rest of the codebase), note it. If not, consider also matching `b'g'`.

Looking at the `contains_subsequence` function in `adapter.rs`, it also does exact byte matching without case folding, so this is consistent. No change needed, but worth a comment.

### M4. The `homopolymer_trim_index` function signature could be more future-proof

The plan derives the complement internally (`A<->T`, `G<->C`). If a future need arises (e.g., poly-N trimming or something exotic), the caller would need to modify the internals. Consider passing `target_5prime` as a second parameter instead of deriving it, or documenting that the complement is always the standard Watson-Crick complement.

This is minor -- the current design is clean and the A/T/G/C case covers all realistic homopolymers.

---

## Low Severity

### L1. The `--poly_g` CLI flag description mentions "auto-detected from the data" as the default

The help text says: "By default, poly-G trimming is auto-detected from the data." This is a behavior change from the perspective of existing users who are accustomed to no poly-G trimming. While this is the correct design, it should be prominently noted in release notes / changelog. The plan itself should flag this as a **behavior change** for the implementer to document.

### L2. Consider adding `--poly_g_min_length` for power users

fastp exposes `--poly_g_min_len` (default 10). The plan hardcodes the Cutadapt scoring parameters (min 3bp, +1/-2, 20% error). For the vast majority of use cases this is fine, but a single `--poly_g_min_length` parameter could allow power users to tune sensitivity. This is not essential for the initial implementation but worth noting as a future enhancement.

### L3. Report header line wording

The plan proposes:
```
"Poly-G trimming enabled: removing poly-G tails from 3' end of R1/SE reads, and poly-C heads from 5' end of R2 reads"
```

This parallels the poly-A line, which is good. But for single-end mode, mentioning "R2 reads" could be confusing. The poly-A report line has the same issue, so this is consistent, but a future improvement could conditionally adjust the wording.

### L4. Redundancy between `report.rs::TrimConfig` and `trimmer.rs::TrimConfig`

There are two separate `TrimConfig` structs (one in `report.rs` for report generation, one in `trimmer.rs` for the pipeline). Both need `poly_g: bool` added. The plan correctly identifies both, but this duplication is a pre-existing code smell worth noting -- not something this plan needs to fix.

### L5. Missing test: `--poly_g` with `--nextseq` combined behavior

The plan mentions the `--nextseq` + poly-G edge case in the design section but doesn't include a specific test for it. A test verifying that both features work in tandem (quality-based G-trimming + explicit poly-G removal) would catch any unexpected interaction. For example, after nextseq quality trimming removes some trailing G's, the poly-G step should handle any remaining G's that had quality scores above the cutoff.

### L6. The plan does not mention updating the `main.rs` on-screen paired-end summary

Looking at `main.rs` lines 373-379, there are explicit `eprintln!` statements for poly-A stats in the paired-end summary. Analogous lines for poly-G should be added. The plan covers `report.rs` reporting but does not explicitly mention the stderr summary in `main.rs::run_paired()` and `main.rs::run_single()`.

---

## Algorithm Assessment

The Cutadapt scoring approach (+1 match, -2 mismatch, 20% max error, min 3bp) is well-suited for poly-G trimming. The key question is whether poly-G tails from 2-colour chemistry have different error characteristics than biological poly-A tails.

**Analysis:** On 2-colour instruments, no-signal is called as G with high quality. The "errors" in a poly-G tail are typically real bases interspersed at the boundary between insert and artifact. These look similar to the boundary errors in poly-A tails. The +1/-2 scoring naturally handles this: it finds the best boundary between "mostly G" and "not G", tolerating up to 20% non-G bases within the tail.

fastp's approach (1 mismatch per 8bp, ~12.5% error rate) is slightly stricter but lacks the principled scoring that handles heterogeneous boundaries well. The Cutadapt algorithm is at least as good and likely better for impure tails.

The detection threshold (10 consecutive G's with no mismatches, >0.01% of reads) is mathematically sound. P(10 consecutive G's at 3' end by chance) = 0.25^10 = ~10^-6, giving ~1 expected hit per 1M reads on 4-colour data. The 100-read threshold (100x above noise) provides excellent specificity.

**One consideration:** Very short reads (e.g., 50bp) on 2-colour instruments may have poly-G tails that are proportionally longer (potentially 30-40bp of the 50bp read). The Cutadapt algorithm handles this correctly -- it would trim most of the read, which then gets caught by the length filter. This is the correct behavior.

---

## Summary

The plan is thorough, well-reasoned, and closely follows existing codebase patterns. The most important issue is **H1** -- the auto-detection being unavailable when adapter auto-detection is skipped. This affects a real-world use case (NovaSeq users who specify `--illumina`) and needs to be addressed in the plan before implementation. The remaining items are suggestions for robustness and completeness.
