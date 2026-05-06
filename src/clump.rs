//! Optional read reordering for tighter gzip compression (`--clumpy`).
//!
//! When `--clumpy` is enabled, the reader thread routes reads to N in-memory
//! bin buffers keyed by a canonical k-mer minimizer. When a bin's accumulated
//! raw bytes exceed the per-bin byte budget, the bin is sorted by minimizer
//! key (so reads sharing the same minimizer land adjacent inside the gzip
//! member) and shipped to a worker as one batch. Each flushed bin becomes
//! one gzip member; gzip's 32 KB sliding window then finds long redundant
//! runs of similar sequences within the sorted bin, shaving roughly 20–35%
//! off output size on typical Illumina FASTQ.
//!
//! The routing model is the streaming N-bin dispatcher recommended in the
//! plan at `~/.claude/plans/i-m-curious-if-we-clever-codd.md`. Bin assignment
//! is FNV-1a hash of the minimizer mod N (load-balancing); the compression
//! win comes from the **in-bin sort**, not from cross-bin adjacency, so
//! locality-preserving prefix bucketing is unnecessary (matches stevekm/squish
//! `dev2`'s clump bucket strategy).

use anyhow::{Result, bail};

use crate::fastq::FastqRecord;

/// k-mer length for canonical minimizer (matches stevekm/squish `dev2`).
/// Fixed at 16 so the encoded form fits exactly into a `u32` (16 bases × 2
/// bits each), letting canonicalisation and sliding-window updates run as
/// single-cycle integer ops.
pub const KMER_LEN: usize = 16;

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
pub const MIN_BIN_BYTES: u64 = 1024 * 1024;

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
}

impl ClumpLayout {
    /// Predicted peak RSS in bytes for a run with this layout and `cores`
    /// workers. Inverse of the formula in `resolve_layout`; calling this from
    /// `main.rs` keeps the startup banner honest with what the layout was sized
    /// against.
    pub fn predicted_peak_bytes(&self, cores: usize) -> u64 {
        let dyn_bytes =
            self.bin_byte_budget as u64 * (5 * self.n_bins as u64 + 7 * cores as u64) / 4;
        STATIC_OVERHEAD_BYTES + dyn_bytes
    }
}

/// Fixed memory overhead (FastQC histograms, allocator retention, Rust runtime,
/// gzip decoder state, IO buffers) reserved out of `--memory` before sizing the
/// bin pool. Calibrated empirically against macOS `/usr/bin/time -l` "peak
/// memory footprint" on 31.5 M-record ATAC-seq runs (May 2026); Linux RSS
/// bookkeeping is leaner, so the constant is conservative there too.
const STATIC_OVERHEAD_BYTES: u64 = 512 * 1024 * 1024;

/// Compute `(n_bins, bin_byte_budget)` from a memory budget and core count.
///
/// The goal is **peak RSS ≤ memory_budget**. Three big chunks consume RAM
/// during a clumpy run, plus a roughly fixed overhead:
///
/// 1. Reader's resident bins: `n_bins × bin_byte_budget` of FASTQ text +
///    `Vec<FastqRecord>` spine entries (72 bytes per record on top of ~350 byte
///    text records ≈ 25% extra).
/// 2. Worker input batches in flight: up to `cores × bin_byte_budget` of text
///    + matching spine.
/// 3. Worker output Vecs growing during gzip-encode: roughly `0.5 ×` input
///    size at L1 (smaller at higher gzip levels, but use L1 as the upper bound).
///
/// Combine: peak ≈ STATIC + B × [(n_bins + cores) × 1.25 + cores × 0.5]
///                = STATIC + B × (5 × n_bins + 7 × cores) / 4
///
/// Solving for B given the user's budget gives the formula below. The 1.25 ×
/// spine factor is the doubling-free post-pre-size value; the 0.5 × output
/// factor is the L1 worst case.
///
///   n_bins          = max(16, 4 × cores)
///   usable          = memory_budget − STATIC_OVERHEAD_BYTES
///   bin_byte_budget = 4 × usable / (5 × n_bins + 7 × cores)
///
/// If the derived budget falls below `MIN_BIN_BYTES`, bails — better to fail
/// loudly than silently produce a degraded output.
pub fn resolve_layout(memory_budget_bytes: u64, cores: usize) -> Result<ClumpLayout> {
    if cores == 0 {
        bail!("clumpy layout requires at least one worker core");
    }
    let n_bins = (16_usize).max(4 * cores);
    let usable = memory_budget_bytes.saturating_sub(STATIC_OVERHEAD_BYTES);
    let denom = 5 * (n_bins as u64) + 7 * (cores as u64);
    let bin_byte_budget = (usable.saturating_mul(4)) / denom;
    if bin_byte_budget < MIN_BIN_BYTES {
        bail!(
            "--memory budget too small for --cores {cores}: after reserving {} MiB \
             for static overhead (FastQC, allocator, gzip state), the derived bin \
             budget is {} bytes — below the {}-byte per-bin floor. Increase --memory \
             (try ≥ {} MiB) or decrease --cores.",
            STATIC_OVERHEAD_BYTES / (1024 * 1024),
            bin_byte_budget,
            MIN_BIN_BYTES,
            (STATIC_OVERHEAD_BYTES + (MIN_BIN_BYTES * denom).div_ceil(4)) / (1024 * 1024),
        );
    }
    Ok(ClumpLayout {
        n_bins,
        bin_byte_budget: bin_byte_budget as usize,
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
        let bases = [b'A', b'C', b'G', b'T'];
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

    #[test]
    fn resolve_layout_default() {
        // 4 GiB + 2 cores: usable=3.5 GiB, n_bins=16, denom=5×16+7×2=94,
        // bin_byte_budget = 4×3584 MiB / 94 ≈ 152 MiB.
        let layout = resolve_layout(4 * 1024 * 1024 * 1024, 2).unwrap();
        assert_eq!(layout.n_bins, 16);
        assert!(layout.bin_byte_budget >= 150 * 1024 * 1024);
        assert!(layout.bin_byte_budget < 156 * 1024 * 1024);
    }

    #[test]
    fn resolve_layout_scales_with_cores() {
        // 4 GiB + 8 cores: usable=3.5 GiB, n_bins=32, denom=5×32+7×8=216,
        // bin_byte_budget = 4×3584 MiB / 216 ≈ 66 MiB.
        let layout = resolve_layout(4 * 1024 * 1024 * 1024, 8).unwrap();
        assert_eq!(layout.n_bins, 32);
        assert!(layout.bin_byte_budget >= 65 * 1024 * 1024);
        assert!(layout.bin_byte_budget < 70 * 1024 * 1024);
    }

    #[test]
    fn resolve_layout_predicted_peak_fits_budget() {
        // The whole point of the formula: predicted peak ≤ user-supplied --memory.
        for (mem_gib, cores) in [(2u64, 4), (4, 8), (8, 8), (16, 16)] {
            let budget = mem_gib * 1024 * 1024 * 1024;
            let layout = resolve_layout(budget, cores).unwrap();
            let predicted_peak = layout.predicted_peak_bytes(cores);
            assert!(
                predicted_peak <= budget,
                "predicted peak {} > budget {} for {} GiB / {} cores",
                predicted_peak,
                budget,
                mem_gib,
                cores
            );
        }
    }

    #[test]
    fn resolve_layout_bails_when_too_small() {
        // 560 MiB - 512 MiB STATIC = 48 MiB usable. With 8 cores: denom=216,
        // bin_byte_budget = 4 × 48 MiB / 216 ≈ 0.89 MiB < 1 MiB floor.
        let res = resolve_layout(560 * 1024 * 1024, 8);
        assert!(res.is_err(), "expected bail, got {res:?}");
        let msg = format!("{}", res.unwrap_err());
        assert!(msg.contains("--memory"));
        assert!(msg.contains("static overhead"));
    }

    #[test]
    fn resolve_layout_bails_below_static_overhead() {
        // Below 512 MiB STATIC, usable saturates to 0 → budget is 0 → fails.
        let res = resolve_layout(256 * 1024 * 1024, 4);
        let msg = format!("{}", res.unwrap_err());
        assert!(msg.contains("--memory"));
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
        // ↔ "@x/2"). This is the invariant clumpy must never break.
        for (a, b) in r1.iter().zip(r2.iter()) {
            let a_stem = a.id.trim_end_matches("/1");
            let b_stem = b.id.trim_end_matches("/2");
            assert_eq!(a_stem, b_stem, "pair lockstep broken: {} ≠ {}", a.id, b.id);
        }
    }
}
