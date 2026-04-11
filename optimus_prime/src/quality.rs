//! Quality trimming using the BWA algorithm.
//!
//! This implements the same quality trimming algorithm used by Cutadapt,
//! which is based on the algorithm from BWA (Li & Durbin, 2009).
//!
//! The algorithm finds the longest suffix of the read where the average
//! quality is below the cutoff, using a running sum approach.

/// Trim low-quality bases from the 3' end using the BWA algorithm.
///
/// Walks from the 3' end, computing a running sum of (cutoff - quality) for
/// each base. The trim point is where the running sum reaches its maximum.
/// This effectively finds the longest low-quality suffix to remove.
///
/// # Arguments
/// * `qualities` - ASCII quality string bytes
/// * `cutoff` - Phred quality score cutoff
/// * `phred_offset` - ASCII offset (33 for Phred33, 64 for Phred64)
///
/// # Returns
/// Number of bases to keep (trim position from 5' end).
pub fn quality_trim_3prime(qualities: &[u8], cutoff: u8, phred_offset: u8) -> usize {
    if qualities.is_empty() {
        return 0;
    }

    let mut max_sum: i64 = 0;
    let mut running_sum: i64 = 0;
    let mut trim_pos = qualities.len();

    for i in (0..qualities.len()).rev() {
        let q = qualities[i].saturating_sub(phred_offset) as i64;
        running_sum += (cutoff as i64) - q;
        // Early termination: if the sum drops below 0, the 3' end is high-quality
        // and there's nothing to trim. This matches BWA's bwa_trim_read and
        // Cutadapt's quality_trim_index behavior.
        if running_sum < 0 {
            break;
        }
        if running_sum > max_sum {
            max_sum = running_sum;
            trim_pos = i;
        }
    }

    if max_sum > 0 {
        trim_pos
    } else {
        qualities.len() // no trimming needed
    }
}

/// NextSeq/2-colour quality trimming variant.
///
/// Identical to `quality_trim_3prime` except that G bases have their quality
/// overridden to 0. On NextSeq/NovaSeq platforms, no-signal basecalls are
/// reported as high-quality G — this function ensures those trailing poly-G
/// artifacts are trimmed.
///
/// This matches Cutadapt's `--nextseq-trim` behavior.
pub fn quality_trim_3prime_nextseq(
    sequence: &[u8],
    qualities: &[u8],
    cutoff: u8,
    phred_offset: u8,
) -> usize {
    if qualities.is_empty() {
        return 0;
    }

    let mut max_sum: i64 = 0;
    let mut running_sum: i64 = 0;
    let mut trim_pos = qualities.len();

    for i in (0..qualities.len()).rev() {
        let q = if sequence[i] == b'G' {
            (cutoff as i64) - 1 // NextSeq: G bases contribute +1 to running sum (matches Cutadapt)
        } else {
            qualities[i].saturating_sub(phred_offset) as i64
        };
        running_sum += (cutoff as i64) - q;
        if running_sum < 0 {
            break;
        }
        if running_sum > max_sum {
            max_sum = running_sum;
            trim_pos = i;
        }
    }

    if max_sum > 0 {
        trim_pos
    } else {
        qualities.len()
    }
}

/// Find the start index of a poly-A tail (3' end) or poly-T head (5' end).
///
/// Uses Cutadapt's scoring algorithm: +1 per matching base, -2 per mismatch.
/// The best scoring position is returned, subject to a maximum 20% error rate.
/// Tails shorter than 3 bases are ignored.
///
/// When `revcomp` is false: scans from 3' end looking for poly-A tail,
/// returns the trim position (keep bases 0..index).
///
/// When `revcomp` is true: scans from 5' end looking for poly-T head,
/// returns the clip position (remove bases 0..index).
pub fn poly_a_trim_index(sequence: &[u8], revcomp: bool) -> usize {
    let n = sequence.len();
    if n < 3 {
        return if revcomp { 0 } else { n };
    }
    let mut best_score: i32 = 0;
    let mut score: i32 = 0;
    let mut errors: usize = 0;

    if revcomp {
        // Scan from 5' end for poly-T head
        let mut best_index: usize = 0;
        for i in 0..n {
            if sequence[i] == b'T' {
                score += 1;
            } else {
                score -= 2;
                errors += 1;
            }
            if score > best_score && errors * 5 <= i + 1 {
                best_score = score;
                best_index = i + 1;
            }
        }
        // Ignore heads shorter than 3bp
        if best_index < 3 { 0 } else { best_index }
    } else {
        // Scan from 3' end for poly-A tail
        let mut best_index: usize = n;
        for i in (0..n).rev() {
            if sequence[i] == b'A' {
                score += 1;
            } else {
                score -= 2;
                errors += 1;
            }
            if score > best_score && errors * 5 <= n - i {
                best_score = score;
                best_index = i;
            }
        }
        // Ignore tails shorter than 3bp
        if best_index > n - 3 { n } else { best_index }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_trimming_high_quality() {
        // All bases Q40 (ASCII 'I' = 73, Phred33 offset = 33, so Q = 40)
        // Cutoff 20: all bases are above, no trimming
        let qual = b"IIIIIIIIIII";
        let pos = quality_trim_3prime(qual, 20, 33);
        assert_eq!(pos, qual.len());
    }

    #[test]
    fn test_trim_all_low_quality() {
        // All bases Q2 (ASCII '#' = 35, Q = 2)
        // Cutoff 20: all bases below, trim everything
        let qual = b"#####";
        let pos = quality_trim_3prime(qual, 20, 33);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_trim_3prime_low_quality_tail() {
        // IIIII##### -> Q40,Q40,Q40,Q40,Q40,Q2,Q2,Q2,Q2,Q2
        // Should trim the last 5 bases
        let qual = b"IIIII#####";
        let pos = quality_trim_3prime(qual, 20, 33);
        assert_eq!(pos, 5);
    }

    #[test]
    fn test_bwa_running_sum_logic() {
        // The BWA algorithm considers the running sum, so a single low-quality
        // base in the middle of high-quality bases won't trigger trimming
        // if the sum doesn't become positive from the 3' end
        let qual = b"IIII#IIII";
        let pos = quality_trim_3prime(qual, 20, 33);
        // The running sum from 3' end:
        // pos 8: (20 - 40) = -20, sum = -20
        // pos 7: -20, sum = -40
        // pos 6: -20, sum = -60
        // pos 5: -20, sum = -80
        // pos 4: (20 - 2) = 18, sum = -62
        // ... sum never goes positive => no trimming
        assert_eq!(pos, qual.len());
    }

    #[test]
    fn test_empty_qualities() {
        let pos = quality_trim_3prime(b"", 20, 33);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_phred64() {
        // With Phred64 offset: ASCII 'h' = 104, Q = 104 - 64 = 40
        // ASCII 'B' = 66, Q = 66 - 64 = 2
        let qual = b"hhhhhhBBBB";
        let pos = quality_trim_3prime(qual, 20, 64);
        assert_eq!(pos, 6);
    }

    #[test]
    fn test_single_base_high() {
        let pos = quality_trim_3prime(b"I", 20, 33);
        assert_eq!(pos, 1);
    }

    #[test]
    fn test_single_base_low() {
        let pos = quality_trim_3prime(b"#", 20, 33);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_gradual_quality_drop() {
        // Simulate a gradual quality drop: Q30, Q25, Q20, Q15, Q10, Q5
        // Phred33: Q30='?'(63), Q25=':'(58), Q20='5'(53), Q15='0'(48), Q10='+'(43), Q5='&'(38)
        let qual = &[63u8, 58, 53, 48, 43, 38];
        let pos = quality_trim_3prime(qual, 20, 33);
        // From 3' end: running sums:
        // i=5: 20-5=15, sum=15, max=15, pos=5
        // i=4: 20-10=10, sum=25, max=25, pos=4
        // i=3: 20-15=5, sum=30, max=30, pos=3
        // i=2: 20-20=0, sum=30, 30 > 30? No -> pos stays at 3
        // i=1: 20-25=-5, sum=25
        // i=0: 20-30=-10, sum=15
        // => trim_pos = 3 (keep 3 bases: Q30, Q25, Q20)
        // Note: Cutadapt uses strict > (not >=), so ties keep the less
        // aggressive trim position, matching BWA's bwa_trim_read.
        assert_eq!(pos, 3);
    }

    // --- NextSeq / 2-colour tests ---

    #[test]
    fn test_nextseq_trailing_g_trimmed() {
        // High-quality trailing G bases should be trimmed in nextseq mode
        // Sequence: ACGTGGGG, all Q40 ('I')
        // Without nextseq: no trimming (all high quality)
        // With nextseq: G bases get Q=0, so the trailing GGGG are trimmed
        let seq = b"ACGTGGGG";
        let qual = b"IIIIIIII";
        let pos_normal = quality_trim_3prime(qual, 20, 33);
        let pos_nextseq = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos_normal, 8); // no trimming
        assert_eq!(pos_nextseq, 4); // trims 4 trailing Gs
    }

    #[test]
    fn test_nextseq_non_g_preserved() {
        // Non-G high-quality trailing bases are NOT trimmed in nextseq mode
        let seq = b"ACGTACGT";
        let qual = b"IIIIIIII";
        let pos = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos, 8); // no trimming — no trailing Gs
    }

    #[test]
    fn test_nextseq_mixed_tail() {
        // Mixed G and non-G at 3' end: only contiguous Gs from the end
        // ACGTACGG — the GG at the end are G, but the A before breaks the chain
        // With BWA running-sum, the G bases get Q=0 (cutoff-0=20 each),
        // but the A before has Q=40 (cutoff-40=-20), which resets the sum.
        let seq = b"ACGTACGG";
        let qual = b"IIIIIIII";
        let pos = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos, 6); // trims 2 trailing Gs
    }

    #[test]
    fn test_nextseq_low_quality_g_also_trimmed() {
        // G bases with low quality — should still be trimmed (Q=0 either way)
        let seq = b"ACGTGGGG";
        let qual = b"IIII####"; // trailing Gs already low quality
        let pos = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos, 4);
    }

    #[test]
    fn test_nextseq_all_g() {
        // Entire read is G — trims everything
        let seq = b"GGGGGGGG";
        let qual = b"IIIIIIII";
        let pos = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_nextseq_empty() {
        let pos = quality_trim_3prime_nextseq(b"", b"", 20, 33);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_nextseq_interior_g_not_overtrimmed() {
        // Interior G bases must NOT pull the trim point inward past
        // high-quality non-G bases. Cutadapt sets G quality to cutoff-1,
        // giving each G a contribution of +1 to the running sum — not
        // enough to overcome a high-quality non-G base (which contributes -20).
        // Sequence: AGGAGG (all Q40, cutoff 20)
        // Expected: trim trailing GG only → pos=4 (keep AGGA)
        // Bug (quality=0): would return 1 (interior GG overwhelms the A)
        let seq = b"AGGAGG";
        let qual = b"IIIIII"; // all Q40
        let pos = quality_trim_3prime_nextseq(seq, qual, 20, 33);
        assert_eq!(pos, 4); // trims only trailing GG
    }

    // --- Poly-A / Poly-T trimming tests ---

    #[test]
    fn test_poly_a_clear_tail() {
        // ACGTACGTAAAAAAAAAA — 10 trailing As
        let seq = b"ACGTACGTAAAAAAAAAA";
        let pos = poly_a_trim_index(seq, false);
        assert_eq!(pos, 8); // keep first 8 bases
    }

    #[test]
    fn test_poly_a_no_tail() {
        let seq = b"ACGTACGTACGT";
        let pos = poly_a_trim_index(seq, false);
        assert_eq!(pos, 12); // no trimming
    }

    #[test]
    fn test_poly_a_short_tail_ignored() {
        // Tails shorter than 3bp are ignored
        let seq = b"ACGTACGTAA";
        let pos = poly_a_trim_index(seq, false);
        assert_eq!(pos, 10); // 2 As — too short, no trimming
    }

    #[test]
    fn test_poly_a_exactly_3() {
        let seq = b"ACGTACGTAAA";
        let pos = poly_a_trim_index(seq, false);
        assert_eq!(pos, 8); // 3 As — minimum for trimming
    }

    #[test]
    fn test_poly_a_with_errors() {
        // Poly-A tail with one mismatch: AAAAAGAAAAA (G at position 5 in the tail)
        let seq = b"ACGTAAAAAGAAAAA";
        let pos = poly_a_trim_index(seq, false);
        // 10 As + 1 error = 1/11 = 9% error rate, under 20% threshold
        assert_eq!(pos, 4); // trims the entire A-rich tail including the G error
    }

    #[test]
    fn test_poly_a_all_a() {
        let seq = b"AAAAAAAAAA";
        let pos = poly_a_trim_index(seq, false);
        assert_eq!(pos, 0); // entire read is poly-A
    }

    #[test]
    fn test_poly_a_empty() {
        let pos = poly_a_trim_index(b"", false);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_poly_t_clear_head() {
        // TTTTTTTTTTACGTACGT — 10 leading Ts
        let seq = b"TTTTTTTTTTACGTACGT";
        let pos = poly_a_trim_index(seq, true);
        assert_eq!(pos, 10); // clip first 10 bases
    }

    #[test]
    fn test_poly_t_no_head() {
        let seq = b"ACGTACGTACGT";
        let pos = poly_a_trim_index(seq, true);
        assert_eq!(pos, 0); // no trimming
    }

    #[test]
    fn test_poly_t_short_head_ignored() {
        let seq = b"TTACGTACGT";
        let pos = poly_a_trim_index(seq, true);
        assert_eq!(pos, 0); // 2 Ts — too short
    }

    #[test]
    fn test_poly_t_with_errors() {
        // TTTTTCTTTTTACGT — 10 Ts with 1 error
        let seq = b"TTTTTCTTTTTACGT";
        let pos = poly_a_trim_index(seq, true);
        assert_eq!(pos, 11); // clips 11bp poly-T head (including the C error)
    }
}
