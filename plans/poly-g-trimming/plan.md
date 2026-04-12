# Plan: Poly-G Trimming with Auto-Detection

## Context

On 2-colour chemistry instruments (NovaSeq, NextSeq, NovaSeqX), no-signal bases are called as G with high quality scores. This produces spurious poly-G tails at the 3' end of reads when the sequencer reads through to the end of short inserts. These are pure artifact and should be removed.

The Oxidized Edition already has `--poly_a` trimming (for RNA-seq poly-A tails) and `--nextseq` (quality-based G-trimming via BWA algorithm). Poly-G trimming is complementary to `--nextseq`:

- `--nextseq` overrides G quality to 0 during *quality* trimming — only works when the quality cutoff is high enough to trigger trimming.
- `--poly_g` explicitly removes poly-G homopolymer tails regardless of quality scores — catches long high-quality poly-G runs that `--nextseq` might miss.

## Design decisions

### 1. Data-driven auto-detection (not read-name matching)

fastp detects 2-colour instruments by matching read name prefixes (`@A`, `@NS`, `@VH`, etc.) against a hardcoded list. This is brittle:
- New instruments require code updates
- Read names may be rewritten by SRA, seqtk, or preprocessing tools
- Third-party platforms using 2-colour chemistry would be missed

Instead, we detect poly-G enrichment directly from the data — consistent with TrimGalore's adapter auto-detection philosophy. The first 1M reads are already scanned for adapter detection; we piggyback on that same pass.

**Detection algorithm:** During the adapter auto-detection loop (`adapter.rs`, the existing `while let Some(record) = reader.next_record()` loop that scans up to `MAX_SCAN_READS = 1_000_000`), count reads whose 3' end has ≥10 consecutive G bases.

**Why 10 consecutive G's (no mismatches) for detection:** The detection threshold is deliberately stricter than the trimming algorithm. For detection, we want zero false positives — 10 consecutive G's is unmistakable. The probability of 10+ trailing G's by chance on a 4-colour instrument is (0.25)^10 ≈ 10^-6, so in 1M reads we'd expect ~1. On a 2-colour instrument, the count will be in the thousands to tens of thousands.

**Threshold:** If >0.01% of scanned reads (>100 per 1M) have ≥10 trailing G's, enable poly-G trimming. This gives a ~100x margin above the random baseline (~1 expected) while easily catching even mild 2-colour artifacts.

**Is 1M reads enough?** Yes, more than enough. The signal-to-noise ratio is enormous:
- 4-colour instrument: ~1 read in 1M with 10+ trailing G's (random chance)
- 2-colour instrument: typically 5,000–50,000+ reads in 1M (0.5–5%+)
- Detection threshold: 100 reads (0.01%)
- Even a low-artifact 2-colour sample will have 10–100x more poly-G reads than the threshold

### 2. CLI: auto-detect by default, with overrides

```
--poly_g         Force-enable poly-G trimming (skip auto-detection)
--no_poly_g      Force-disable poly-G trimming (skip auto-detection)
(neither)        Auto-detect from first 1M reads (default)
```

This mirrors the adapter auto-detection model: most users never touch it, power users can override.

### 3. Separate from `--poly_a`

Different biology (poly-A = mRNA tails, poly-G = 2-colour no-signal artifact). A WGBS user on NovaSeq wants poly-G but not poly-A. An RNA-seq user on HiSeq wants the opposite. Can be combined freely.

### 4. Reuse the Cutadapt scoring algorithm for trimming

The existing +1 match, -2 mismatch, max 20% error rate, min 3bp algorithm is more principled than fastp's "1 mismatch per 8 bases" approach and is already implemented. It handles impure tails gracefully (e.g., `GGGGAGGGGG` with 1 error in 10 = 10% gets trimmed).

Note: the trimming threshold (min 3bp with mismatches allowed) is intentionally more sensitive than the detection threshold (10bp, no mismatches). Detection must be conservative to avoid false positives; once we know the data is from a 2-colour instrument, we want to catch even short tails.

### 5. Generalise `poly_a_trim_index()` rather than duplicate it

The current function hardcodes `b'A'` and `b'T'`. Refactor to accept the target base as a parameter, then call it for both poly-A and poly-G.

### 6. Paired-end strand logic

Same pattern as poly-A:
- R1/SE: trim poly-G from 3' end
- R2: trim poly-C from 5' end (reverse complement of the poly-G tail)

### 7. Pipeline position

Run at step 2.8 in `trim_read()`, immediately after poly-A (step 2.7). Both run after adapter trimming and RRBS, before N-trimming.

### 8. Report separately

Poly-G trimming stats are reported independently from poly-A.

## Files to change

| File | Change |
|------|--------|
| `quality.rs` | Generalise `poly_a_trim_index()` → `homopolymer_trim_index(seq, base, revcomp)` |
| `adapter.rs` | Add poly-G counting to the auto-detection loop; return poly-G count in `DetectionResult` |
| `cli.rs` | Add `--poly_g` and `--no_poly_g` flags |
| `trimmer.rs` | Add `poly_g: bool` to `TrimConfig`, add poly-G step at 2.8, extend `TrimResult` |
| `report.rs` | Add `poly_g_trimmed` / `poly_g_bases_trimmed` stats + report output |
| `main.rs` | Wire auto-detection result + CLI overrides → config |
| `parallel.rs` | Collect poly-G stats (same pattern as poly-A) |

## Step-by-step changes

### 1. `quality.rs` — Generalise the scoring function

Rename and refactor:

```rust
/// Find the start of a homopolymer tail (3' end) or head (5' end).
///
/// `target_3prime`: the base to match at the 3' end (e.g. b'A' or b'G').
/// The 5' complement is derived automatically (A↔T, G↔C).
///
/// Uses Cutadapt's scoring algorithm: +1 per matching base, -2 per mismatch.
/// Max 20% error rate. Tails shorter than 3 bases are ignored.
pub fn homopolymer_trim_index(sequence: &[u8], target_3prime: u8, revcomp: bool) -> usize
```

- When `revcomp == false`: scan from 3' end matching `target_3prime` (G for poly-G, A for poly-A)
- When `revcomp == true`: scan from 5' end matching the complement (C for poly-G, T for poly-A)
- The complement mapping: `A→T, T→A, G→C, C→G`

Keep `poly_a_trim_index()` as a thin wrapper for backwards compatibility with existing tests:
```rust
pub fn poly_a_trim_index(sequence: &[u8], revcomp: bool) -> usize {
    homopolymer_trim_index(sequence, b'A', revcomp)
}
```

Add unit tests for poly-G:
- Clear poly-G tail: `ACGTACGTGGGGGGGGG` → trim 9 G's
- No poly-G: `ACGTACGTACGT` → no trimming
- Short tail (2 G's): ignored
- Exactly 3 G's: trimmed
- Poly-G with errors: `GGGGGAGGGGG` → trimmed (1/11 = 9% error)
- Poly-C head (revcomp=true): `CCCCCCCCCCACGT` → clip 10 C's
- Mixed: sequence ending in `...AAAGGG` with poly-A already removed — poly-G should still find the G's

### 2. `adapter.rs` — Poly-G auto-detection

**Add `poly_g_count` to `DetectionResult`:**
```rust
pub struct DetectionResult {
    pub adapter: AdapterPreset,
    pub counts: Vec<(String, usize)>,
    pub reads_scanned: usize,
    pub message: String,
    pub suppressed: bool,
    /// Number of reads with ≥10 consecutive trailing G's (for poly-G auto-detection)
    pub poly_g_count: usize,
}
```

**Add counting inside the existing scan loop** (`autodetect_adapter()`):

```rust
let mut poly_g_count: usize = 0;

while let Some(record) = reader.next_record()? {
    reads_scanned += 1;
    if reads_scanned > MAX_SCAN_READS { break; }

    let seq = record.seq.as_bytes();

    // Existing adapter counting
    for (i, (_, adapter_seq)) in adapters.iter().enumerate() {
        if contains_subsequence(seq, adapter_seq.as_bytes()) {
            counts[i] += 1;
        }
    }

    // Poly-G detection: count reads with ≥10 consecutive trailing G's
    if has_trailing_poly_g(seq, 10) {
        poly_g_count += 1;
    }
}
```

**Helper function** (simple, fast — no scoring needed for detection):
```rust
/// Check if a sequence ends with `min_len` or more consecutive G bases.
fn has_trailing_poly_g(seq: &[u8], min_len: usize) -> bool {
    if seq.len() < min_len { return false; }
    let mut count = 0;
    for &base in seq.iter().rev() {
        if base == b'G' {
            count += 1;
            if count >= min_len { return true; }
        } else {
            break;
        }
    }
    false
}
```

**Add to the auto-detection report message** when poly-G is detected:
```
2-colour chemistry detected: X.X% of reads have poly-G tails (≥10bp). Poly-G trimming enabled.
```

### 3. `cli.rs` — Add flags

```rust
/// Trim poly-G tails from the 3' end of Read 1 (and single-end reads),
/// and poly-C heads from the 5' end of Read 2. Useful for data from
/// 2-colour instruments (NovaSeq, NextSeq) where no-signal bases are
/// called as high-quality G. By default, poly-G trimming is auto-detected
/// from the data. Use this flag to force-enable it.
#[clap(long = "poly_g", alias = "poly-g", alias = "polyG",
       conflicts_with = "no_poly_g")]
pub poly_g: bool,

/// Disable poly-G auto-detection and trimming.
#[clap(long = "no_poly_g", alias = "no-poly-g", alias = "no-polyG")]
pub no_poly_g: bool,
```

Place adjacent to the existing `--poly_a` flag.

### 4. `main.rs` — Wire auto-detection + CLI overrides

The logic after adapter auto-detection. **Every code path prints a message** so the user always knows what happened:

```rust
// Determine poly-G trimming: CLI overrides auto-detection
let poly_g_enabled = if cli.poly_g {
    // Force-enabled by user
    true
} else if cli.no_poly_g {
    // Force-disabled by user
    false
} else {
    // Auto-detect from adapter scan
    let threshold = detection.reads_scanned / 10_000; // 0.01%
    let threshold = threshold.max(10); // at least 10 reads (for small files)
    detection.poly_g_count > threshold
};

// --- On-screen message (always printed for record-keeping) ---
let poly_g_pct = if detection.reads_scanned > 0 {
    detection.poly_g_count as f64 / detection.reads_scanned as f64 * 100.0
} else { 0.0 };

if cli.poly_g {
    eprintln!(
        "Poly-G trimming: ENABLED (user-specified --poly_g). \
         Auto-detection found {:.2}% poly-G reads in {}.",
        poly_g_pct, detection.reads_scanned
    );
} else if cli.no_poly_g {
    eprintln!(
        "Poly-G trimming: DISABLED (user-specified --no_poly_g). \
         Auto-detection found {:.2}% poly-G reads in {}.",
        poly_g_pct, detection.reads_scanned
    );
} else if poly_g_enabled {
    eprintln!(
        "Poly-G trimming: ENABLED (auto-detected). \
         {:.2}% of reads have poly-G tails (≥10bp) — \
         consistent with 2-colour chemistry (NovaSeq/NextSeq). \
         To disable: --no_poly_g",
        poly_g_pct
    );
} else {
    eprintln!(
        "Poly-G trimming: not enabled (auto-detection found {:.2}% poly-G reads — \
         below 0.01% threshold). To force-enable: --poly_g",
        poly_g_pct
    );
}

// ... later, when building TrimConfig:
config.poly_g = poly_g_enabled;
```

This ensures:
- Users always see whether poly-G trimming is on or off
- The auto-detection percentage is always shown (even when overridden), so users can sanity-check
- The message explains *why* (user-specified vs auto-detected) and *how to change it*
- Pipeline logs capture the decision for reproducibility

### 5. `trimmer.rs` — Add to TrimConfig and trim_read()

**TrimConfig**: add `pub poly_g: bool` alongside `poly_a`.

**TrimResult**: add `pub poly_g_trimmed: usize` (bases trimmed, 0 = none).

**trim_read()**: After the existing poly-A block (step 2.7), add poly-G as step 2.8:

```rust
// 2.8. Poly-G / Poly-C trimming (2-colour chemistry artifact removal)
// R1/SE: trim poly-G from 3' end. R2: trim poly-C from 5' end.
let mut poly_g_trimmed: usize = 0;
if config.poly_g && !record.is_empty() {
    let revcomp = is_r2;
    let seq_len_before = record.seq.len();
    let idx = quality::homopolymer_trim_index(record.seq.as_bytes(), b'G', revcomp);
    if revcomp {
        if idx > 0 {
            record.clip_5prime(idx);
            poly_g_trimmed = idx;
        }
    } else {
        if idx < seq_len_before {
            record.truncate(idx);
            poly_g_trimmed = seq_len_before - idx;
        }
    }
}
```

Update stat collection in `run_single_end()` and `run_paired_end()`:
```rust
if result.poly_g_trimmed > 0 {
    stats.poly_g_trimmed += 1;
    stats.poly_g_bases_trimmed += result.poly_g_trimmed;
}
```

### 6. `report.rs` — Add stats and report lines

**TrimStats**: add two fields:
```rust
pub poly_g_trimmed: usize,
pub poly_g_bases_trimmed: usize,
```

**TrimStats::merge()**: add:
```rust
self.poly_g_trimmed += other.poly_g_trimmed;
self.poly_g_bases_trimmed += other.poly_g_bases_trimmed;
```

**ReportConfig**: add `pub poly_g: bool`.

**write_report_header()**: when `config.poly_g`:
```
"Poly-G trimming enabled: removing poly-G tails from 3' end of R1/SE reads, and poly-C heads from 5' end of R2 reads"
```

**write_report_stats()**: when `stats.poly_g_trimmed > 0`:
```
Reads with poly-G/C trimmed:    N (X.X%)
  Poly-G/C bases removed:       N
```

### 7. `parallel.rs` — Collect stats in parallel path

Same pattern as poly-A — add stat collection for `poly_g_trimmed` and `poly_g_bases_trimmed` in both `process_pairs()` and `process_reads()`.

## Testing

1. **Unit tests** (quality.rs): poly-G tail, poly-C head, edge cases (see step 1).

2. **Unit tests** (adapter.rs): `has_trailing_poly_g()` — sequences with and without trailing G's, edge cases.

3. **Auto-detection test**: Create a synthetic FASTQ where ~5% of reads end in 20 G's. Run without flags → verify poly-G trimming is auto-enabled. Run with `--no_poly_g` → verify it's disabled.

4. **Integration**: Run on NovaSeq test data and verify poly-G tails are removed. If no NovaSeq test data is available, create a synthetic FASTQ with known poly-G tails.

5. **Regression**: Run full test suite without `--poly_g` on existing 4-colour test data → verify auto-detection does NOT trigger, no behavior changes.

6. **Combined**: Test `--poly_a --poly_g` together to verify both run and report independently.

## Edge cases

- **Very short files** (< 1000 reads): The `threshold.max(10)` guard ensures at least 10 reads must have poly-G before auto-enabling. With 100 reads, even 10% poly-G wouldn't trigger — this is fine, false negatives on tiny files are harmless (user can always pass `--poly_g`).
- **Reads that are entirely poly-G**: Trimmed to zero length, then caught by the minimum length filter (default 20bp) and removed.
- **poly-A followed by poly-G** (`--poly_a --poly_g`): Poly-A runs first (step 2.7), removes the A tail, then poly-G runs (step 2.8) on the remainder. Order matters — a read ending `...AAAGGGGG` would have the A's trimmed first, then the G's. This is correct.
- **`--nextseq` + auto-detected poly-G**: Both can be active simultaneously. `--nextseq` handles the quality-trimming side, poly-G handles explicit homopolymer removal. They're complementary, not redundant.

## Review-driven amendments

The following changes address findings from the dual independent plan review (see `review_reviewer-A.md` and `review_reviewer-B.md`).

### A1. Standalone `detect_poly_g()` for when adapter auto-detection is skipped (H1 — both reviewers)

**Problem:** `resolve_adapter()` in `main.rs` has early-return paths for `--illumina`, `--adapter`, `--nextera`, etc. that skip `autodetect_adapter()`. The original plan piggybacked poly-G counting on that scan loop, so NovaSeq users who specify an adapter explicitly would get no auto-detection.

**Fix:** Add a standalone `detect_poly_g(path) -> (usize, usize)` function in `adapter.rs` that scans up to 1M reads counting only trailing poly-G. This runs whenever:
- Auto-detection was skipped (adapter was user-specified or preset-selected), AND
- Neither `--poly_g` nor `--no_poly_g` was passed

When `autodetect_adapter()` does run, the poly-G count from its piggybacked counting is used instead (no double scan).

```rust
/// Standalone poly-G detection scan (used when adapter auto-detection is skipped).
/// Returns (poly_g_count, reads_scanned).
pub fn detect_poly_g<P: AsRef<Path>>(path: P) -> Result<(usize, usize)>
```

### A2. Console summary for poly-G stats in `main.rs` (M4/L6 — both reviewers)

Add `eprintln!` blocks for poly-G stats in `run_single()` and `run_paired()`, matching the existing poly-A pattern. Print stats when poly-G was enabled, even if count is 0 (confirms the feature ran).

### A3. Fix `ReportConfig` → `report::TrimConfig` naming (L4 — Reviewer A)

The actual struct is `report::TrimConfig`, not `ReportConfig`. Plan text was inaccurate; implementation uses correct name.

### A4. Clarify `--nextseq` independence in CLI help (H2 — Reviewer A)

Add to `--poly_g` help text: "This is independent from --nextseq (quality-based G-trimming)."

### A5. Print poly-G report stats even when 0 if explicitly enabled (M3 — Reviewer A)

When `config.poly_g` is true, always print the poly-G stats line in the report (even if 0 reads were trimmed), matching adapter stats behavior. When auto-detected, only print if count > 0.

### A6. Add `--nextseq` + `--poly_g` combined test (L5 — Reviewer B)

Add a test verifying both features work in tandem.
