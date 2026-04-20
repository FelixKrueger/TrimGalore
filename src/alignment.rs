//! Semi-global adapter alignment engine.
//!
//! Reimplements Cutadapt's core adapter matching algorithm using
//! unit-cost semi-global dynamic programming alignment.

/// Result of an adapter alignment against a read.
#[derive(Debug, Clone, PartialEq)]
pub struct AdapterMatch {
    /// Start position in the read where the adapter match begins (0-based)
    pub read_start: usize,
    /// End position in the read (exclusive) where the adapter match ends
    pub read_end: usize,
    /// Number of adapter bases that matched (overlap length)
    pub overlap: usize,
    /// Number of errors (mismatches + insertions + deletions) in the match
    pub errors: usize,
}

/// Find the best 3' adapter match in a read using semi-global alignment.
///
/// The adapter is expected at the 3' end of the read. It may:
/// - Occur entirely within the read (full adapter match)
/// - Start within the read and extend beyond it (partial overlap at 3' end)
///
/// Both cases are handled by a single DP matrix:
/// - **Full adapter matches** from the last row `dp[m][j]`: scanned left-to-right,
///   the first (leftmost) valid match is selected — this maximizes trimming.
///   The error count is corrected via `dp[m][read_start + m]` to avoid inflation
///   from trailing deletion artifacts in the DP.
/// - **Partial 3' matches** from the last column `dp[i][n]`: scanned top-to-bottom
///   for the longest valid overlap.
///
/// Between the two, the match that trims more (smaller `read_start`) wins.
///
/// Uses unit-cost DP: mismatch = 1, insertion = 1, deletion = 1.
/// First row initialized to 0 (adapter can start at any read position).
/// Exact `read_start` is recovered via backtrace through the full DP matrix.
pub fn find_3prime_adapter(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f64,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    if read.is_empty() || adapter.is_empty() || min_overlap == 0 {
        return None;
    }

    let n = read.len();
    let m = adapter.len();

    // Build full DP matrix — for typical adapters (m=13) and reads (n≤150),
    // this is ~2KB. Storing the full matrix enables exact backtrace.
    //
    // dp[i][j] = min edit distance to align adapter[0..i] against some suffix of read[0..j]
    // First row: dp[0][j] = 0 for all j (adapter can start anywhere — free read prefix skip)
    // First column: dp[i][0] = i (aligning i adapter bases to nothing costs i deletions)
    let mut dp = vec![vec![0usize; n + 1]; m + 1];
    for (i, row) in dp.iter_mut().enumerate().skip(1) {
        row[0] = i;
    }

    for i in 1..=m {
        for j in 1..=n {
            let cost_sub = if adapter[i - 1] == read[j - 1] { 0 } else { 1 };

            dp[i][j] = (dp[i - 1][j - 1] + cost_sub) // match/mismatch
                .min(dp[i - 1][j] + 1) // deletion (gap in read)
                .min(dp[i][j - 1] + 1); // insertion (gap in adapter)
        }
    }

    // --- Find best full match: leftmost valid j, then correct error count ---
    // Cutadapt picks the leftmost valid match in the last row — this maximizes
    // trimming for 3' adapters. However, the leftmost j may have an inflated
    // error count due to a trailing deletion artifact (e.g., dp[m][j-1]=1 via
    // deletion when dp[m][j]=0 is exact). To get the true error count for the
    // alignment at that read_start, we also check dp[m][read_start + m].
    let max_errors_full = (max_error_rate * m as f64).floor() as usize;
    let mut best_full: Option<AdapterMatch> = None;
    for j in 1..=n {
        if dp[m][j] <= max_errors_full {
            let read_start = backtrace_start(&dp, adapter, read, m, j);
            // Correct the error count: check the "natural" alignment end
            // (adapter starting at read_start with no net indels)
            let natural_end = read_start + m;
            let errors = if natural_end <= n && dp[m][natural_end] <= max_errors_full {
                dp[m][natural_end].min(dp[m][j])
            } else {
                dp[m][j]
            };
            let read_end = if natural_end <= n && dp[m][natural_end] <= dp[m][j] {
                natural_end
            } else {
                j
            };
            best_full = Some(AdapterMatch {
                read_start,
                read_end,
                overlap: m,
                errors,
            });
            break; // Leftmost valid = most trimming
        }
    }

    // --- Find best partial match: longest valid overlap from last column ---
    let mut best_partial: Option<AdapterMatch> = None;
    let max_partial = m.min(n);
    for i in (min_overlap..=max_partial).rev() {
        let max_errors = (max_error_rate * i as f64).floor() as usize;
        if dp[i][n] <= max_errors {
            let read_start = backtrace_start(&dp, adapter, read, i, n);
            best_partial = Some(AdapterMatch {
                read_start,
                read_end: n,
                overlap: i,
                errors: dp[i][n],
            });
            break; // Longest valid = most trimming
        }
    }

    // Compare full match vs partial match: prefer the one that trims more
    // (smaller read_start). On ties, prefer the full match (longer overlap).
    match (best_full, best_partial) {
        (Some(f), Some(p)) => {
            if f.read_start <= p.read_start {
                Some(f)
            } else {
                Some(p)
            }
        }
        (Some(f), None) => Some(f),
        (None, p) => p,
    }
}

/// Backtrace through the DP matrix from cell (i, j) to row 0.
/// Returns the read start position (the column where the alignment begins).
///
/// Follows the optimal path: prefers diagonal (match/mismatch) over gaps,
/// and deletion (gap in read) over insertion (gap in adapter). This gives
/// the tightest alignment with the earliest start position.
fn backtrace_start(
    dp: &[Vec<usize>],
    adapter: &[u8],
    read: &[u8],
    mut i: usize,
    mut j: usize,
) -> usize {
    while i > 0 {
        // Prefer diagonal (match/mismatch) first
        if j > 0 {
            let cost = if adapter[i - 1] == read[j - 1] { 0 } else { 1 };
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                i -= 1;
                j -= 1;
                continue;
            }
        }
        // Then deletion (gap in read = adapter base not matched to any read base)
        if dp[i][j] == dp[i - 1][j] + 1 {
            i -= 1;
            continue;
        }
        // Then insertion (gap in adapter = extra read base consumed)
        if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            j -= 1;
            continue;
        }
        break; // shouldn't happen in a valid DP
    }
    j // column at row 0 = start position in the read
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match_at_3prime() {
        // Read ends with exact adapter sequence
        let read = b"ACGTACGTAAGAGATCGGAAGAGC";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 11);
        assert_eq!(m.overlap, 13);
        assert_eq!(m.errors, 0);
    }

    #[test]
    fn test_partial_adapter_at_3prime() {
        // Only first 5 bases of adapter present at end of read
        let read = b"ACGTACGTACGTAGATC";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.overlap, 5);
        assert_eq!(m.errors, 0);
        assert_eq!(m.read_start, 12);
    }

    #[test]
    fn test_no_match() {
        let read = b"ACGTACGTACGT";
        let adapter = b"TTTTTTTTTTTT";
        let m = find_3prime_adapter(read, adapter, 0.1, 5);
        assert!(m.is_none());
    }

    #[test]
    fn test_one_mismatch_allowed() {
        // Adapter with 1 mismatch in 13bp match = 7.7% error rate (< 10%)
        let read = b"ACGTACGTAGATCGGAAGAGX"; // X instead of C at end
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1);
        assert!(m.is_some());
        let m = m.unwrap();
        assert_eq!(m.errors, 1);
        assert_eq!(m.overlap, 13);
    }

    #[test]
    fn test_min_overlap_filter() {
        // Only 2 bases match but min_overlap is 5
        let read = b"ACGTACGTAG";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 5);
        assert!(m.is_none());
    }

    #[test]
    fn test_empty_inputs() {
        assert!(find_3prime_adapter(b"", b"AGATC", 0.1, 1).is_none());
        assert!(find_3prime_adapter(b"ACGT", b"", 0.1, 1).is_none());
    }

    /// Regression test: leftmost valid match wins over a better-quality match
    /// further right. This matches Cutadapt's 3' adapter behavior.
    #[test]
    fn test_leftmost_match_wins() {
        // Read has adapter at position 1 (1 mismatch) AND a perfect match later.
        // Cutadapt picks the leftmost valid match for 3' adapters.
        let read = b"AAGATCGGAAGAGCAACCAATGAGATCTCGTATGCCGTCTTCTGCTTG";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        // Must find the leftmost match (position ~1), not a later perfect match
        assert!(
            m.read_start <= 1,
            "Expected read_start <= 1, got {}",
            m.read_start
        );
        assert_eq!(m.overlap, 13);
    }

    /// Regression test: when adapter occurs twice in the read, the earlier
    /// occurrence is selected even if it has more errors (still within threshold).
    #[test]
    fn test_two_adapter_occurrences_picks_first() {
        // Adapter at pos 18 with 1 mismatch, adapter at pos 128 with 0 errors.
        // Cutadapt picks pos 18 (leftmost, trims more).
        let read = b"CTGTGCCTCTTCCGATTCAGATCGGAAGAGAGTAGAAGCTACACG";
        //                              ^AGATCGGAAGAGA (1 mismatch vs AGATCGGAAGAGC)
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 18);
        assert_eq!(m.overlap, 13);
        assert_eq!(m.errors, 1);
    }

    #[test]
    fn test_adapter_at_very_start_of_read() {
        // Entire read is adapter — should trim to 0
        let read = b"AGATCGGAAGAGC";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 0);
        assert_eq!(m.overlap, 13);
        assert_eq!(m.errors, 0);
    }

    #[test]
    fn test_single_base_read_with_adapter_start() {
        // 1-base read matching first base of adapter
        let read = b"A";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 0);
        assert_eq!(m.overlap, 1);
        assert_eq!(m.errors, 0);
    }

    #[test]
    fn test_full_match_beats_partial_when_earlier() {
        // Full adapter match at position 5 should beat a partial 3-base match
        // at the 3' end, because the full match trims more.
        let read = b"XXXXXAGATCGGAAGAGCXXXXXAGA";
        //                ^full match at 5      ^partial AGA at 3' end
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 5);
        assert_eq!(m.overlap, 13);
        assert_eq!(m.errors, 0);
    }

    #[test]
    fn test_partial_match_when_no_full_match() {
        // No full adapter match, but partial 3-base match at 3' end
        let read = b"XXXXXXXXXXXXXAGA";
        let adapter = b"AGATCGGAAGAGC";
        let m = find_3prime_adapter(read, adapter, 0.1, 1).unwrap();
        assert_eq!(m.read_start, 13);
        assert_eq!(m.overlap, 3);
        assert_eq!(m.errors, 0);
    }
}
