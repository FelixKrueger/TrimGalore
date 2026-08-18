//! Optional read reordering for tighter gzip compression (`--clumpify`).
//!
//! When `--clumpify` is enabled, the reader thread routes reads to N in-memory
//! bin buffers keyed by a canonical k-mer minimizer. When a bin's accumulated
//! raw bytes exceed the per-bin byte budget, the bin is sorted by minimizer
//! key (so reads sharing the same minimizer land adjacent inside the gzip
//! member) and shipped to a worker as one batch. Each flushed bin becomes
//! one gzip member; gzip's 32 KB sliding window then finds long redundant
//! runs of similar sequences within the sorted bin, shaving roughly 20–35%
//! off output size on typical Illumina FASTQ.
//!
//! Bin assignment is FNV-1a hash of the minimizer mod N (load-balancing); the
//! compression win comes from the **in-bin sort**, not from cross-bin
//! adjacency, so locality-preserving prefix bucketing is unnecessary (matches
//! stevekm/squish `dev2`'s clump bucket strategy).

use anyhow::{Result, bail};

use crate::fastq::FastqRecord;

/// k-mer length for canonical minimizer (matches stevekm/squish `dev2`).
/// Fixed at 16 so the encoded form fits exactly into a `u32` (16 bases × 2
/// bits each), letting canonicalisation and sliding-window updates run as
/// single-cycle integer ops.
const KMER_LEN: usize = 16;

/// Canonical minimizer key — the lexicographically smallest canonical
/// (forward-or-reverse-complement) k-mer over a read, encoded as 2-bit
/// packed `u32` (A=00, C=01, G=10, T=11). Wider than a byte so sort
/// comparisons reduce to a single `u32` compare instead of a 16-byte
/// memcmp; small enough to fit in a register.
pub type MinimizerKey = u32;

/// Minimum acceptable per-bin byte budget. Below this, gzip header overhead
/// starts to erode the compression win and we'd rather fail fast than ship a
/// degraded output. The dispatcher's `resolve_layout` bails when the derived
/// budget falls below this floor.
const MIN_BIN_BYTES: u64 = 1024 * 1024;

/// Map a sequence base to its 2-bit code. Anything that isn't ACGT
/// (including `N`, lowercase noise, ambiguous IUPAC codes) folds to A.
/// Real FASTQ data from any modern instrument has <0.1% N so the
/// resulting bin assignment skew is negligible.
#[inline(always)]
fn encode_2bit(b: u8) -> u32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Reverse-complement a 2-bit-packed 16-mer in O(1) bitwise ops.
///
/// 1. Bit-pair reverse: swap pair-within-byte, swap nibbles, swap bytes.
/// 2. Complement: each 2-bit pair flipped (A↔T, C↔G) — a single XOR with
///    `0xFFFF_FFFF` (each pair `00↔11`, `01↔10`).
#[inline(always)]
fn revcomp_2bit_u32(mut x: u32) -> u32 {
    x = ((x & 0xCCCC_CCCC) >> 2) | ((x & 0x3333_3333) << 2);
    x = ((x & 0xF0F0_F0F0) >> 4) | ((x & 0x0F0F_0F0F) << 4);
    x = x.swap_bytes();
    x ^ 0xFFFF_FFFF
}

#[inline(always)]
fn canonical_2bit(x: u32) -> u32 {
    let rc = revcomp_2bit_u32(x);
    if rc < x { rc } else { x }
}

/// Compute the canonical minimizer of `seq`: the lexicographically smallest
/// canonical (forward-or-revcomp) 16-mer over all positions.
///
/// Reads shorter than `KMER_LEN` are right-padded with `A` (the encoded
/// form of `N`) so every read still produces a deterministic key.
///
/// Implementation: maintain a 2-bit packed sliding window of the forward
/// k-mer; revcomp + canonical reduce to a handful of integer ops per
/// position. ~30× faster per base than the byte-array form this replaced.
pub fn canonical_minimizer(seq: &[u8]) -> MinimizerKey {
    if seq.len() < KMER_LEN {
        let mut x: u32 = 0;
        for i in 0..KMER_LEN {
            x = (x << 2)
                | if i < seq.len() {
                    encode_2bit(seq[i])
                } else {
                    0
                };
        }
        return canonical_2bit(x);
    }

    let mut fwd: u32 = 0;
    for &b in &seq[..KMER_LEN] {
        fwd = (fwd << 2) | encode_2bit(b);
    }
    let mut best = canonical_2bit(fwd);

    // Slide the window across the rest of the read. Shifting `fwd` left 2
    // bits naturally drops the oldest base (out the top of the u32) and
    // makes room for the new one — exactly what a 16-base sliding window
    // needs when the encoding is exactly 32 bits wide.
    for &b in &seq[KMER_LEN..] {
        fwd = (fwd << 2) | encode_2bit(b);
        let cand = canonical_2bit(fwd);
        if cand < best {
            best = cand;
        }
    }
    best
}

/// FNV-1a 32-bit hash of the minimizer key, mod `n_bins`.
///
/// We hash (rather than mod the key directly) so the bin distribution is
/// even even when minimizer keys are structurally biased. The dominant
/// compression win comes from the in-bin sort, not cross-bin locality —
/// even dispersion is the right objective for the bin index.
pub fn bin_for(key: MinimizerKey, n_bins: usize) -> usize {
    debug_assert!(n_bins > 0);
    let mut h: u32 = 0x811c_9dc5;
    for byte in key.to_be_bytes() {
        h ^= byte as u32;
        h = h.wrapping_mul(0x0100_0193);
    }
    (h as usize) % n_bins
}

// ─── Memory budget arithmetic ─────────────────────────────────────────────

/// Resolved layout for the bin dispatcher.
#[derive(Debug, Clone, Copy)]
pub struct ClumpLayout {
    pub n_bins: usize,
    pub bin_byte_budget: usize,
    /// Worker count this layout was sized for; the prediction is only valid for
    /// that count.
    cores: usize,
    /// Static reservation this layout was sized against.
    static_bytes: u64,
}

impl ClumpLayout {
    /// Layout with an explicit bin count and budget, for callers that bypass
    /// `resolve_layout`. The reservation is the trim-phase constant.
    ///
    /// Test-only: a budget below `MIN_BIN_BYTES` makes `predicted_peak_bytes`
    /// meaningless, so production layouts come from `resolve_layout` alone.
    #[cfg(test)]
    pub(crate) fn new(n_bins: usize, bin_byte_budget: usize, workers: usize) -> Self {
        debug_assert!(n_bins > 0 && bin_byte_budget as u64 >= MIN_BIN_BYTES);
        Self {
            n_bins,
            bin_byte_budget,
            cores: workers,
            static_bytes: STATIC_TRIM_BYTES,
        }
    }

    /// Predicted peak RSS in bytes. Inverse of the formula in `resolve_layout`,
    /// so the startup banner cannot disagree with the sizing.
    pub fn predicted_peak_bytes(&self) -> u64 {
        let dyn_bytes = self.bin_byte_budget as u64 * dyn_denominator(self.n_bins, self.cores) / 16;
        self.static_bytes + dyn_bytes
    }
}

/// `16 × (σ·n_bins + k·cores)` — the resident multiple of one bin budget, with
/// σ = 23/16 for the bin pool and k = 4 per worker.
fn dyn_denominator(n_bins: usize, cores: usize) -> u64 {
    23 * n_bins as u64 + 64 * cores as u64
}

/// Per-worker channel depths. `k` above charges for `work + 1 in hand +
/// result_per_core`, so their sum is the constant and neither may move alone.
#[derive(Debug, Clone, Copy)]
pub struct ChannelDepths {
    pub work: usize,
    pub result_per_core: usize,
}

/// A clumpy batch is a whole bin, so a shallower work queue caps resident bins.
pub fn channel_depths(clumpy: bool) -> ChannelDepths {
    if clumpy {
        ChannelDepths {
            work: 2, // VALIDATION 8 form (a): helper edited
            result_per_core: 2,
        }
    } else {
        ChannelDepths {
            work: 2,
            result_per_core: 2,
        }
    }
}

/// The three thread counts the memory model needs. They differ: `--clump_only`
/// is synchronous, and `--fastqc_args -t N` sets FastQC's threads independently.
#[derive(Debug, Clone, Copy)]
pub struct LayoutInputs {
    /// Drives `n_bins`. Keep it the user's `--cores` so bin grouping — and
    /// therefore output layout — does not shift with the other two.
    pub bin_cores: usize,
    /// Workers that can hold batches in flight; 1 on the synchronous paths.
    pub workers: usize,
    /// FastQC's thread count, or `None` when no report was requested.
    pub fastqc_threads: Option<usize>,
}

impl LayoutInputs {
    /// Every count equal, for the parallel trimming path.
    pub fn uniform(cores: usize, fastqc_threads: Option<usize>) -> Self {
        Self {
            bin_cores: cores,
            workers: cores,
            fastqc_threads,
        }
    }

    fn n_bins(&self) -> usize {
        (16_usize).max(4 * self.bin_cores)
    }
}

/// Static reservation for a run, in bytes.
fn static_bytes_for(inputs: &LayoutInputs) -> u64 {
    match inputs.fastqc_threads {
        Some(threads) => {
            let charged = threads.min(FASTQC_CHARGED_THREAD_CAP) as u64;
            STATIC_TRIM_BYTES + STATIC_FASTQC_PER_THREAD_BYTES.saturating_mul(charged)
        }
        None => STATIC_TRIM_BYTES,
    }
}

/// Minimum `--memory` (in bytes) these inputs need for the bin pool to clear
/// `MIN_BIN_BYTES`. Used by `main.rs` to decide whether to warn-and-fall-back to
/// plain mode rather than bail.
pub fn clumpify_min_memory_bytes(inputs: &LayoutInputs) -> u64 {
    let dyn_min = (MIN_BIN_BYTES * dyn_denominator(inputs.n_bins(), inputs.workers)).div_ceil(16);
    ((static_bytes_for(inputs) + dyn_min) * MARGIN_DEN).div_ceil(MARGIN_NUM)
}

/// Resident cost outside the bin pool: Rust runtime, gzip state, IO buffers and
/// allocator retention. Measured peak on 10 M-pair paired runs is 205 MiB on
/// macOS/arm64 and 173 MiB on Linux/x86_64.
const STATIC_TRIM_BYTES: u64 = 224 * 1024 * 1024;

/// FastQC's additional resident cost per thread. Its phase follows trimming but
/// the pool's pages stay resident, so it adds rather than maxes.
/// See `plans/08162026_clumpify-memory-accounting/phase2/CALIBRATION.md`.
const STATIC_FASTQC_PER_THREAD_BYTES: u64 = 24 * 1024 * 1024;

/// Threads beyond this are not charged, so the floor cannot run away on
/// high-core hosts. Scaling above it is unmeasured.
const FASTQC_CHARGED_THREAD_CAP: usize = 16;

/// The budget is sized to `MARGIN_NUM / MARGIN_DEN` of `--memory`, leaving ~9%
/// for run-to-run variation. Largest observed within-cell spread is 96 MiB.
const MARGIN_NUM: u64 = 10;
const MARGIN_DEN: u64 = 11;

/// Compute `(n_bins, bin_byte_budget)` from a memory budget, core count and
/// whether FastQC was requested.
///
/// The goal is **peak RSS ≤ memory_budget**. Resident memory is a fixed term
/// plus a multiple of one bin budget `B`:
///
/// 1. Reader's resident bins: `σ × n_bins × B`, text plus `Vec<FastqRecord>`
///    spine, σ = 23/16.
/// 2. Per worker: `k × B` for the queued batch, the batch in hand, and the
///    compressed output waiting on the result channel, k = 4.
///
///   n_bins          = max(16, 4 × bin_cores)
///   usable          = memory_budget × 10/11 − static reservation
///   bin_byte_budget = 16 × usable / (23 × n_bins + 64 × workers)
///
/// σ and k are fitted from a 10 M-pair paired-end matrix over `cores ∈ {2,3,4}`
/// with `n_bins` pinned at 16, taking the worst of two reps per cell; `cores ∈
/// {8,16}` were held out and predicted within 5%. Measured values are σ = 1.29
/// (Linux) / 1.41 (macOS) and k = 3.71 / 3.99; the constants take the worse of
/// each. The 10/11 factor is margin: the largest within-cell spread was 96 MiB.
/// See `plans/08162026_clumpify-memory-accounting/phase2/CALIBRATION.md`.
///
/// If the derived budget falls below `MIN_BIN_BYTES`, bails — better to fail
/// loudly than silently produce a degraded output.
pub fn resolve_layout(memory_budget_bytes: u64, inputs: &LayoutInputs) -> Result<ClumpLayout> {
    if inputs.bin_cores == 0 || inputs.workers == 0 {
        bail!("clumpify layout requires at least one worker core");
    }
    let n_bins = inputs.n_bins();
    let static_bytes = static_bytes_for(inputs);
    let budgeted = memory_budget_bytes.saturating_mul(MARGIN_NUM) / MARGIN_DEN;
    let usable = budgeted.saturating_sub(static_bytes);
    let denom = dyn_denominator(n_bins, inputs.workers);
    let bin_byte_budget = (usable.saturating_mul(16)) / denom;
    if bin_byte_budget < MIN_BIN_BYTES {
        bail!(
            "--memory budget too small for --cores {}: after reserving {} MiB \
             for static overhead (allocator, gzip state, IO buffers{}) and 9% margin, \
             the derived bin budget is {} bytes — below the {}-byte per-bin floor. \
             Increase --memory (try ≥ {} MiB) or decrease --cores.",
            inputs.bin_cores,
            static_bytes / (1024 * 1024),
            if inputs.fastqc_threads.is_some() {
                ", FastQC"
            } else {
                ""
            },
            bin_byte_budget,
            MIN_BIN_BYTES,
            clumpify_min_memory_bytes(inputs).div_ceil(1024 * 1024),
        );
    }
    Ok(ClumpLayout {
        n_bins,
        bin_byte_budget: bin_byte_budget as usize,
        cores: inputs.workers,
        static_bytes,
    })
}

/// Parse a memory-size string like `"512M"`, `"2G"`, `"1024K"`, `"1500"`.
///
/// Suffixes `K`/`M`/`G` (case-insensitive) are 1024-based (KiB / MiB / GiB).
/// No suffix means raw bytes.
pub fn parse_memory_size(s: &str) -> Result<u64> {
    let s = s.trim();
    if s.is_empty() {
        bail!("memory size cannot be empty");
    }
    let last = s.chars().last().unwrap();
    let (num_part, mult): (&str, u64) = match last {
        'K' | 'k' => (&s[..s.len() - 1], 1024),
        'M' | 'm' => (&s[..s.len() - 1], 1024 * 1024),
        'G' | 'g' => (&s[..s.len() - 1], 1024 * 1024 * 1024),
        c if c.is_ascii_digit() => (s, 1),
        _ => bail!("could not parse memory size '{s}': expected NUMBER[K|M|G]"),
    };
    let n: u64 = num_part.trim().parse().map_err(|_| {
        anyhow::anyhow!("could not parse memory size '{s}': expected NUMBER[K|M|G]")
    })?;
    n.checked_mul(mult)
        .ok_or_else(|| anyhow::anyhow!("memory size '{s}' overflows u64"))
}

// ─── Bin buffers + sort ───────────────────────────────────────────────────

/// Estimate the raw FASTQ-text bytes a record will occupy on disk.
/// Used by the dispatcher for per-bin byte-budget accounting.
pub fn estimated_record_bytes(rec: &FastqRecord) -> usize {
    // header + '\n' + seq + '\n' + "+\n" + qual + '\n'
    rec.id.len() + 1 + rec.seq.len() + 1 + 2 + rec.qual.len() + 1
}

/// Sort `records` in place by canonical minimizer first, then sequence,
/// quality, id, and finally input position for total determinism.
///
/// Why minimizer-primary instead of sequence-primary (alpha): empirically
/// (May 2026 benchmark) alpha sort beats minimizer sort by ~3 ppt on
/// amplicon-like data where most reads share a 5' prefix, but
/// underperforms by ~2 ppt on diverse short-read data (WGS/WES) because
/// (a) it only catches forward-strand prefix matches while minimizer
/// mode catches forward+revcomp anchors anywhere in the read, and (b)
/// streaming bin flushes fragment same-prefix runs across multiple gzip
/// members. Minimizer-primary is the robust default across data types.
///
/// `keys` is reordered to match — callers that no longer need them can
/// drop them after the sort. Sort is stable (`Vec::sort_by`), so equal
/// keys preserve insertion order — the deterministic input-position
/// tie-break is implicit.
pub fn sort_single_by_key(records: &mut Vec<FastqRecord>, keys: &mut Vec<MinimizerKey>) {
    debug_assert_eq!(records.len(), keys.len());
    let mut paired: Vec<(MinimizerKey, FastqRecord)> =
        keys.drain(..).zip(records.drain(..)).collect();
    paired.sort_by(|(ka, ra), (kb, rb)| {
        ka.cmp(kb)
            .then_with(|| ra.seq.cmp(&rb.seq))
            .then_with(|| ra.qual.cmp(&rb.qual))
            .then_with(|| ra.id.cmp(&rb.id))
    });
    for (k, r) in paired {
        keys.push(k);
        records.push(r);
    }
}

/// Sort paired-end batches in lockstep by R1's minimizer key. R1 and R2
/// are reordered identically so pair semantics are preserved.
pub fn sort_paired_by_key(
    r1: &mut Vec<FastqRecord>,
    r2: &mut Vec<FastqRecord>,
    keys: &mut Vec<MinimizerKey>,
) {
    debug_assert_eq!(r1.len(), r2.len());
    debug_assert_eq!(r1.len(), keys.len());
    let mut grouped: Vec<(MinimizerKey, FastqRecord, FastqRecord)> = keys
        .drain(..)
        .zip(r1.drain(..))
        .zip(r2.drain(..))
        .map(|((k, a), b)| (k, a, b))
        .collect();
    grouped.sort_by(|(ka, ra, _), (kb, rb, _)| {
        ka.cmp(kb)
            .then_with(|| ra.seq.cmp(&rb.seq))
            .then_with(|| ra.qual.cmp(&rb.qual))
            .then_with(|| ra.id.cmp(&rb.id))
    });
    for (k, a, b) in grouped {
        keys.push(k);
        r1.push(a);
        r2.push(b);
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, seq: &str, qual: &str) -> FastqRecord {
        FastqRecord {
            id: id.to_string(),
            seq: seq.to_string(),
            qual: qual.to_string(),
        }
    }

    #[test]
    fn minimizer_is_deterministic() {
        let s = b"ACGTACGTACGTACGTAACCGGTT";
        assert_eq!(canonical_minimizer(s), canonical_minimizer(s));
    }

    #[test]
    fn minimizer_canonicalises_revcomp() {
        // A read and its reverse complement must hash to the same minimizer
        // — this is what lets reads from opposite strands of the same
        // fragment share a bin.
        let fwd = b"ACGTACGTACGTACGTAAAAA";
        let rc: Vec<u8> = fwd
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect();
        assert_eq!(canonical_minimizer(fwd), canonical_minimizer(&rc));
    }

    #[test]
    fn minimizer_handles_short_reads() {
        // Sequences shorter than KMER_LEN are padded with N — must not panic
        // and must be deterministic.
        let s = b"ACGT";
        let k = canonical_minimizer(s);
        assert_eq!(canonical_minimizer(s), k);
    }

    #[test]
    fn minimizer_handles_empty_read() {
        // 0-length read pads to all-A (the encoded form of N) — exercise
        // it explicitly. canonical_2bit(0) is whichever of (0, revcomp(0))
        // is smaller; revcomp of 0 (all-A) is 0xFFFF_FFFF (all-T), so
        // canonical = 0.
        let k = canonical_minimizer(b"");
        assert_eq!(k, 0);
    }

    #[test]
    fn minimizer_handles_all_n_read() {
        // N folds to A (00) — same as empty read.
        let k = canonical_minimizer(&[b'N'; 50]);
        assert_eq!(k, 0);
    }

    #[test]
    fn minimizer_is_case_insensitive() {
        let lo = b"acgtacgtacgtacgta";
        let hi = b"ACGTACGTACGTACGTA";
        assert_eq!(canonical_minimizer(lo), canonical_minimizer(hi));
    }

    #[test]
    fn bin_for_distributes_reasonably() {
        // Realistic input: minimizers computed from 100K pseudo-random ACGT
        // sequences. Asserts no bin grows beyond ~1.7× the mean count — a
        // loose bound that just guards against catastrophic hash collapse,
        // not a quality claim.
        let n_bins = 64;
        let mut counts = vec![0u32; n_bins];
        let bases = *b"ACGT";
        let mut rng_state: u64 = 0x9E3779B97F4A7C15;
        let mut next_byte = || {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((rng_state >> 33) as usize) & 0x3]
        };
        for _ in 0..100_000_u32 {
            let mut seq = [0u8; 50];
            for b in seq.iter_mut() {
                *b = next_byte();
            }
            let key = canonical_minimizer(&seq);
            counts[bin_for(key, n_bins)] += 1;
        }
        let mean = 100_000.0 / n_bins as f64;
        let max = *counts.iter().max().unwrap() as f64;
        let min = *counts.iter().min().unwrap() as f64;
        assert!(
            max < mean * 1.7,
            "max bin count {max} exceeds 1.7× mean {mean}; counts={counts:?}"
        );
        assert!(
            min > mean * 0.4,
            "min bin count {min} below 0.4× mean {mean}; counts={counts:?}"
        );
    }

    #[test]
    fn parse_memory_size_handles_suffixes() {
        assert_eq!(parse_memory_size("0").unwrap(), 0);
        assert_eq!(parse_memory_size("1024").unwrap(), 1024);
        assert_eq!(parse_memory_size("1K").unwrap(), 1024);
        assert_eq!(parse_memory_size("1k").unwrap(), 1024);
        assert_eq!(parse_memory_size("512M").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory_size("2G").unwrap(), 2 * 1024 * 1024 * 1024);
        assert_eq!(parse_memory_size("  2g  ").unwrap(), 2 * 1024 * 1024 * 1024);
    }

    #[test]
    fn parse_memory_size_rejects_garbage() {
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("abc").is_err());
        assert!(parse_memory_size("12X").is_err());
        assert!(parse_memory_size("M").is_err());
    }

    /// Inputs where all three counts are the user's `--cores`, FastQC optional.
    fn ins(cores: usize, fastqc: bool) -> LayoutInputs {
        LayoutInputs::uniform(cores, fastqc.then_some(cores))
    }

    #[test]
    fn resolve_layout_default() {
        // 4 GiB + 2 cores: usable = 4 GiB×10/11 − 224 MiB, n_bins=16,
        // denom = 23×16 + 64×2 = 496, so B = 16 × usable / 496 ≈ 113 MiB.
        let layout = resolve_layout(4 * 1024 * 1024 * 1024, &ins(2, false)).unwrap();
        assert_eq!(layout.n_bins, 16);
        assert!(layout.bin_byte_budget >= 112 * 1024 * 1024);
        assert!(layout.bin_byte_budget < 114 * 1024 * 1024);
    }

    #[test]
    fn resolve_layout_scales_with_cores() {
        // 4 GiB + 8 cores: n_bins=32, denom = 23×32 + 64×8 = 1248, B ≈ 45 MiB.
        let layout = resolve_layout(4 * 1024 * 1024 * 1024, &ins(8, false)).unwrap();
        assert_eq!(layout.n_bins, 32);
        assert!(layout.bin_byte_budget >= 44 * 1024 * 1024);
        assert!(layout.bin_byte_budget < 46 * 1024 * 1024);
    }

    #[test]
    fn resolve_layout_pins_the_default_cell() {
        // 1 GiB + 4 cores is the cell #439 is about, so the sizing arithmetic may
        // not move it silently: 16 × (1 GiB×10/11 − 224 MiB) / 624.
        let layout = resolve_layout(1024 * 1024 * 1024, &ins(4, false)).unwrap();
        assert_eq!(layout.n_bins, 16);
        assert_eq!(layout.bin_byte_budget, 19_006_356);
    }

    #[test]
    fn sizing_constants_are_pinned_across_the_role_space() {
        // Hand-computed bin budgets. A coefficient edit — σ, k, the static
        // reservation, the FastQC per-thread charge, its cap, or the margin —
        // moves at least one of these, which the swept test below cannot see.
        // Do not trim this table: the 17-thread cell is the only one that
        // observes FASTQC_CHARGED_THREAD_CAP.
        let gib = 1024 * 1024 * 1024u64;
        let cases: [(u64, LayoutInputs, usize); 6] = [
            (gib, LayoutInputs::uniform(4, None), 19_006_356),
            (
                gib,
                LayoutInputs {
                    bin_cores: 4,
                    workers: 1,
                    fastqc_threads: None,
                },
                27_453_626,
            ),
            (
                gib,
                LayoutInputs {
                    bin_cores: 4,
                    workers: 4,
                    fastqc_threads: Some(16),
                },
                8_681_915,
            ),
            (
                gib,
                LayoutInputs {
                    bin_cores: 4,
                    workers: 4,
                    fastqc_threads: Some(17),
                },
                8_681_915,
            ),
            (gib, LayoutInputs::uniform(2, None), 23_911_222),
            (4 * gib, LayoutInputs::uniform(32, None), 11_761_649),
        ];
        for (budget, inputs, expected) in cases {
            let layout = resolve_layout(budget, &inputs).unwrap();
            assert_eq!(
                layout.bin_byte_budget, expected,
                "pin moved for budget {budget} inputs {inputs:?}"
            );
        }
    }

    #[test]
    fn clumpify_floor_at_the_default_cell_is_pinned_absolutely() {
        // Relative floor tests cannot see MIN_BIN_BYTES being lowered, because
        // both sides move together. This pins the constant and the figure the
        // docs publish for --cores 2 in one assertion.
        let floor = clumpify_min_memory_bytes(&ins(2, false));
        assert_eq!(floor.div_ceil(1024 * 1024), 281);
    }

    #[test]
    fn channel_depths_sum_to_the_per_worker_coefficient() {
        // k = 4 in dyn_denominator charges for the queued batch, the batch in
        // hand, and the result slots. Raising either depth invalidates the fit.
        let clumpy = channel_depths(true);
        assert_eq!(clumpy.work, 1);
        assert_eq!(clumpy.result_per_core, 2);
        assert_eq!(
            clumpy.work + 1 + clumpy.result_per_core,
            4,
            "channel depths must sum to k = 4 (64/16 in dyn_denominator); see \
             plans/08162026_clumpify-memory-accounting/phase2/CALIBRATION.md"
        );
        assert_eq!(channel_depths(false).work, 2);
    }

    #[test]
    fn resolve_layout_fastqc_reserves_more_and_shrinks_bins() {
        // FastQC's histograms are charged on top of the trim-phase reservation,
        // so the same budget yields a smaller pool.
        let budget = 1024 * 1024 * 1024;
        let plain = resolve_layout(budget, &ins(4, false)).unwrap();
        let with_qc = resolve_layout(budget, &ins(4, true)).unwrap();
        assert!(
            with_qc.bin_byte_budget < plain.bin_byte_budget,
            "FastQC layout {} should be smaller than {}",
            with_qc.bin_byte_budget,
            plain.bin_byte_budget
        );
        // The difference is the per-core reservation spread over the denominator.
        let expected = (24u64 * 4 * 1024 * 1024 * 16 / 624) as usize;
        let delta = plain.bin_byte_budget - with_qc.bin_byte_budget;
        assert!(
            delta.abs_diff(expected) <= 1024,
            "delta {delta}, expected ~{expected}"
        );
    }

    #[test]
    fn floor_and_sizing_round_consistently_without_overflow() {
        // Not a coefficient guard — `predicted ≤ budget` holds for any σ, k or
        // static, because the sizing and the prediction read the same helpers.
        // What this pins: the floor is the *least* resolvable budget across the
        // input space (a rounding direction that disagrees between the two makes
        // advertised budgets unusable), and the unchecked multiplies in
        // `predicted_peak_bytes` and `clumpify_min_memory_bytes` stay inside u64.
        for cores in [1usize, 2, 3, 4, 8, 16, 17, 32, 1024, 65536] {
            for fastqc in [false, true] {
                let inputs = ins(cores, fastqc);
                let floor = clumpify_min_memory_bytes(&inputs);

                let at = resolve_layout(floor, &inputs);
                assert!(at.is_ok(), "floor {floor} must resolve for {inputs:?}");
                assert!(at.unwrap().bin_byte_budget as u64 >= MIN_BIN_BYTES);
                assert!(
                    resolve_layout(floor - 1, &inputs).is_err(),
                    "one byte under {floor} must bail for {inputs:?}"
                );

                for budget in [floor, floor + 1, u64::MAX / 2, u64::MAX] {
                    if let Ok(layout) = resolve_layout(budget, &inputs) {
                        let predicted = layout.predicted_peak_bytes();
                        assert!(
                            predicted <= budget,
                            "predicted {predicted} > budget {budget} for {inputs:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn resolve_layout_bails_just_below_the_floor() {
        // One byte under the floor must bail, and the message must name the
        // reservation it actually charged.
        let inputs = ins(8, false);
        let floor = clumpify_min_memory_bytes(&inputs);
        let res = resolve_layout(floor - 1, &inputs);
        assert!(
            res.is_err(),
            "expected bail one byte under {floor}, got {res:?}"
        );
        let msg = format!("{}", res.unwrap_err());
        assert!(msg.contains("--memory"));
        assert!(msg.contains("static overhead"));
        assert!(
            msg.contains(&format!("{} MiB", STATIC_TRIM_BYTES / (1024 * 1024))),
            "bail should name the charged reservation: {msg}"
        );
    }

    #[test]
    fn resolve_layout_bails_below_the_reservation() {
        // Under the reservation itself, usable saturates to 0.
        let below = STATIC_TRIM_BYTES - 1024 * 1024;
        let res = resolve_layout(below, &ins(4, false));
        let msg = format!("{}", res.unwrap_err());
        assert!(msg.contains("--memory"));
    }

    #[test]
    fn clumpify_floor_matches_the_smallest_resolvable_budget() {
        // The advertised floor must itself resolve, one byte less must not, and
        // the MiB figure quoted to the user must resolve too — the caller rounds
        // bytes to MiB for the message, so a truncating round would advertise an
        // unusable budget.
        for cores in [1usize, 2, 3, 4, 8, 16, 32] {
            for fastqc in [false, true] {
                let inputs = ins(cores, fastqc);
                let floor = clumpify_min_memory_bytes(&inputs);
                assert!(
                    resolve_layout(floor, &inputs).is_ok(),
                    "advertised floor {floor} does not resolve at {cores} cores (fastqc {fastqc})"
                );
                assert!(
                    resolve_layout(floor - 1, &inputs).is_err(),
                    "a byte below the floor still resolves at {cores} cores (fastqc {fastqc})"
                );
                let advised_mib = floor.div_ceil(1024 * 1024);
                assert!(
                    resolve_layout(advised_mib * 1024 * 1024, &inputs).is_ok(),
                    "the advised {advised_mib} MiB does not resolve at {cores} cores \
                     (fastqc {fastqc})"
                );
            }
        }
    }

    #[test]
    fn fastqc_reservation_scales_with_threads_and_is_capped() {
        // Charged per FastQC thread, so two thread counts must differ...
        let budget = 4 * 1024 * 1024 * 1024;
        let four = LayoutInputs::uniform(4, Some(4));
        let eight = LayoutInputs::uniform(4, Some(8));
        assert!(
            static_bytes_for(&eight) > static_bytes_for(&four),
            "8 FastQC threads should reserve more than 4"
        );
        assert_eq!(
            static_bytes_for(&eight) - static_bytes_for(&four),
            STATIC_FASTQC_PER_THREAD_BYTES * 4
        );
        // ...and beyond the cap it stops growing, so the floor cannot run away.
        let capped = LayoutInputs::uniform(4, Some(FASTQC_CHARGED_THREAD_CAP));
        let over = LayoutInputs::uniform(4, Some(FASTQC_CHARGED_THREAD_CAP * 4));
        assert_eq!(static_bytes_for(&capped), static_bytes_for(&over));
        // Bin budget follows the reservation, not `--cores`.
        let a = resolve_layout(budget, &four).unwrap().bin_byte_budget;
        let b = resolve_layout(budget, &eight).unwrap().bin_byte_budget;
        assert!(a > b, "more FastQC threads must not yield a larger pool");
    }

    #[test]
    fn synchronous_paths_are_not_charged_for_absent_workers() {
        // `--clump_only` runs one worker while still binning by the user's cores,
        // so it must get a larger pool than the parallel path at the same budget.
        let budget = 1024 * 1024 * 1024;
        let parallel = resolve_layout(budget, &LayoutInputs::uniform(8, None)).unwrap();
        let sync = resolve_layout(
            budget,
            &LayoutInputs {
                bin_cores: 8,
                workers: 1,
                fastqc_threads: None,
            },
        )
        .unwrap();
        assert_eq!(sync.n_bins, parallel.n_bins, "bin grouping must not shift");
        assert!(
            sync.bin_byte_budget > parallel.bin_byte_budget,
            "one worker ({}) should buy a bigger pool than eight ({})",
            sync.bin_byte_budget,
            parallel.bin_byte_budget
        );
    }

    #[test]
    fn estimated_record_bytes_matches_disk_format() {
        // 4-line FASTQ block: header\nseq\n+\nqual\n
        let r = rec("@id", "ACGT", "IIII");
        // 3 + 1 + 4 + 1 + 2 + 4 + 1 = 16
        assert_eq!(estimated_record_bytes(&r), 16);
    }

    #[test]
    fn sort_single_groups_by_minimizer() {
        let mut records = vec![
            rec(
                "@a",
                "TTTTTTTTTTTTTTTTTTTTAAAAAA",
                "IIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            rec(
                "@b",
                "AAAAAAAAAAAAAAAAAAAAAAAAAA",
                "IIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            rec(
                "@c",
                "TTTTTTTTTTTTTTTTTTTTAAAAAA",
                "JJJJJJJJJJJJJJJJJJJJJJJJJJ",
            ),
        ];
        let mut keys: Vec<_> = records
            .iter()
            .map(|r| canonical_minimizer(r.seq.as_bytes()))
            .collect();
        sort_single_by_key(&mut records, &mut keys);
        // Records with the same minimizer key end up adjacent. Whichever
        // key compares smaller is first; verify by re-checking adjacency.
        let k0 = canonical_minimizer(records[0].seq.as_bytes());
        let k1 = canonical_minimizer(records[1].seq.as_bytes());
        let k2 = canonical_minimizer(records[2].seq.as_bytes());
        assert!(k0 <= k1 && k1 <= k2, "sorted keys must be non-decreasing");
    }

    #[test]
    fn sort_paired_preserves_lockstep() {
        let mut r1 = vec![
            rec("@a/1", "TTTTTTTTTTTTTTTTAAAAAA", "I".repeat(22).as_str()),
            rec("@b/1", "AAAAAAAAAAAAAAAAAAAAAA", "I".repeat(22).as_str()),
            rec("@c/1", "GGGGGGGGGGGGGGGGGGGGGG", "I".repeat(22).as_str()),
        ];
        let mut r2 = vec![
            rec("@a/2", "MATE_A", "IIIIII"),
            rec("@b/2", "MATE_B", "IIIIII"),
            rec("@c/2", "MATE_C", "IIIIII"),
        ];
        let mut keys: Vec<_> = r1
            .iter()
            .map(|r| canonical_minimizer(r.seq.as_bytes()))
            .collect();
        sort_paired_by_key(&mut r1, &mut r2, &mut keys);

        // For each i, R1[i].id and R2[i].id must still be mates ("@x/1"
        // ↔ "@x/2"). This is the invariant clumpify must never break.
        for (a, b) in r1.iter().zip(r2.iter()) {
            let a_stem = a.id.trim_end_matches("/1");
            let b_stem = b.id.trim_end_matches("/2");
            assert_eq!(a_stem, b_stem, "pair lockstep broken: {} ≠ {}", a.id, b.id);
        }
    }
}
