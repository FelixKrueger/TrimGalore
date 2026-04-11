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
}
