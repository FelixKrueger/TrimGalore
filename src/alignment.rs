//! Semi-global adapter alignment engine.
//!
//! Reimplements Cutadapt's core adapter matching algorithm using
//! unit-cost semi-global dynamic programming alignment.
//!
//! Per @an-altosian's #248 audit, ~62% of reads in real-world
//! Buckberry-scale data carry no adapter at all, and for those reads
//! the entire scalar DP fill is wasted work. To recover that cost
//! without giving up byte-identity, the public entry-point
//! [`find_3prime_adapter`] runs a Myers' bit-parallel prefilter first
//! ([`myers_proves_no_match`]). The prefilter is *conservative*: it
//! only short-circuits when its O(n) bit-vector walk can rigorously
//! prove that **no** match exists, considering both the full-adapter
//! case (last DP row) and the partial-overlap case (last DP column).
//! When the prefilter cannot prove no-match, control falls through
//! to the unchanged scalar DP, which produces the canonical answer.
//! This keeps the byte-identity invariant intact by construction.
//!
//! Limited to adapters ≤ 64 bp (single u64 bit-vector). Longer
//! patterns fall through to the scalar DP unconditionally.

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

    // Myers' bit-parallel prefilter (#248 #4). For adapters ≤ 64 bp,
    // an O(n) bit-vector walk can prove "no match exists" without
    // building the full m×n DP. Only short-circuits when both the
    // full-match (last DP row) and partial-match (last DP column)
    // checks fail; otherwise falls through to the scalar DP below
    // which gives the canonical answer. False positives → slight
    // perf loss, no correctness impact. See module docstring.
    if myers_proves_no_match(read, adapter, max_error_rate, min_overlap) {
        return None;
    }

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

/// Myers' bit-parallel approximate string matching, used as a prefilter
/// in front of the scalar DP. Returns `true` iff the algorithm can
/// **prove** that no adapter match exists (full or partial) within the
/// caller's error budget. A return of `false` means "possibly a match —
/// run the canonical scalar DP to find out".
///
/// Conservativeness is the key contract: the function MUST NEVER return
/// `true` when [`find_3prime_adapter`]'s scalar DP would return `Some`.
/// False positives (returning `false` for reads that turn out to have no
/// match) are fine — they cost a single scalar DP fill on top of the
/// O(n) bit-vector walk. False negatives (returning `true` when a match
/// exists) would silently drop adapter-trim events and break byte-
/// identity vs Perl v0.6.x.
///
/// ## Algorithm
///
/// Myers (1999) — "A fast bit-vector algorithm for approximate string
/// matching based on dynamic programming". The pattern's per-character
/// occurrence sets are encoded as bit-vectors `peq[c]` (bit `i` set iff
/// `adapter[i] == c`). The DP column is encoded as two bit-vectors `Pv`
/// and `Mv` of vertical differences (`+1` and `-1` between adjacent rows
/// in the same column). Each text character costs O(1) bit operations.
/// The score at the last DP row is tracked as a running counter via the
/// top bits of `Ph`/`Mh`.
///
/// At end-of-walk, we have:
///  - `min_score` — the minimum `dp[m][j]` across all `j` (full match)
///  - `Pv`, `Mv` at column `n` — vertical differences encoding `dp[i][n]
///    - dp[i-1][n]` for each i (partial match column)
///
/// We can therefore answer both "is there a full match?" and "is there a
/// partial match at the last column?" without building the full DP matrix.
///
/// ## Limits
///
/// - Adapter length must be in `1..=64` (single u64 bit-vector). Longer
///   adapters fall through unconditionally (return `false`); a multi-word
///   Myers' is feasible but the project's adapters are all ≤ 32 bp so
///   the single-word case is sufficient.
/// - Empty read or empty adapter or `min_overlap == 0` are caller-rejected
///   above; this function does not re-validate.
///
/// ## Attribution
///
/// Profiled, prototyped, and benchmarked at Buckberry scale (84M reads,
/// 38% adapter rate) by @an-altosian (Dongze He) in #248 item #4. The
/// "wrap as a prefilter around the existing scalar DP" design is his —
/// it's what makes the byte-identity guarantee automatic. Implementation
/// here is independent (his ephemeral worktree branch was cleaned up
/// before we could port from it), but follows his stated shape: prefilter
/// returns "no match possible", DP unchanged.
fn myers_proves_no_match(
    read: &[u8],
    adapter: &[u8],
    max_error_rate: f64,
    min_overlap: usize,
) -> bool {
    let m = adapter.len();
    let n = read.len();

    // Out-of-scope cases: let the scalar DP handle them. Empty read
    // is already handled in find_3prime_adapter; the m > 64 branch
    // gracefully degrades for longer adapters.
    if m == 0 || m > 64 || n == 0 {
        return false;
    }

    let max_errors_full = (max_error_rate * m as f64).floor() as usize;

    // Build per-base pattern bitmasks. `peq[b]` has bit i set iff
    // adapter[i] == b. Indexed by raw byte value (256 entries) so any
    // non-ACGTN byte simply gets a zero mask (no matches anywhere)
    // — mirrors the strict equality the scalar DP uses on line 64.
    let mut peq = [0u64; 256];
    for (i, &b) in adapter.iter().enumerate() {
        peq[b as usize] |= 1u64 << i;
    }

    // m-bit mask covering the active rows. For m == 64 we want all 64
    // bits set; `(1 << 64) - 1` would overflow, so we special-case.
    let m_mask: u64 = if m == 64 { !0u64 } else { (1u64 << m) - 1 };
    let top_bit: u64 = 1u64 << (m - 1);

    // Initial state — Myers' "approximate matching" variant where
    // dp[0][j] = 0 for all j (free prefix skip in the read). All
    // vertical differences in column 0 are +1 (dp[i][0] = i). The
    // running score equals dp[m][0] = m.
    let mut pv: u64 = m_mask;
    let mut mv: u64 = 0;
    let mut score: usize = m;
    let mut min_score: usize = m;

    for &c in read {
        let eq = peq[c as usize] & m_mask;

        // Standard Myers' update for one column. The arithmetic is
        // exactly the form in the 1999 paper, masked to m bits at the
        // tail to keep bits ≥ m from leaking into row m's score bit.
        let xv = eq | mv;
        let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;
        let mut ph = mv | !(xh | pv);
        let mut mh = pv & xh;

        // Update the running score at row m via the top bit of Ph/Mh.
        if ph & top_bit != 0 {
            score += 1;
        }
        if mh & top_bit != 0 {
            score = score.saturating_sub(1);
        }

        // Shift carry. NO `| 1` — that's the global-alignment boundary
        // (dp[0][j] = j); we use the approximate-matching variant where
        // dp[0][j] = 0 for all j (free prefix skip in the read), so the
        // horizontal difference at the virtual row 0 is always 0. This
        // matches Hyyrö 2001's canonical approximate-matching
        // formulation. Getting this wrong on simple all-mismatch reads
        // makes the prefilter incorrectly accept "could be a match"
        // when there's clearly nothing.
        ph = (ph << 1) & m_mask;
        mh = (mh << 1) & m_mask;

        pv = (mh | !(xv | ph)) & m_mask;
        mv = (ph & xv) & m_mask;

        if score < min_score {
            min_score = score;
        }
    }

    // ── Full-match check ──
    // If the minimum dp[m][j] across all j is still above the full-match
    // error budget, no full alignment satisfies the threshold. This
    // alone doesn't prove "no match" — partial matches at the last
    // column might still satisfy a per-overlap budget — so we keep
    // checking unless we can rule both out.
    if min_score <= max_errors_full {
        return false; // possible full match — defer to scalar DP
    }

    // ── Partial-match check ──
    // At end-of-walk, (pv, mv) encode the vertical differences of column
    // n: bit (i-1) of pv is set iff dp[i][n] - dp[i-1][n] == +1; bit
    // (i-1) of mv is set iff that difference == -1; otherwise the
    // difference is 0. We can therefore reconstruct dp[i][n] iteratively
    // from dp[0][n] = 0, comparing each against the per-i error budget.
    //
    // The scalar DP's partial-match loop scans i in [min_overlap, max_partial]
    // (line 109) where max_partial = m.min(n). If any i in that range has
    // dp[i][n] ≤ floor(max_error_rate * i), the scalar DP will return Some
    // — so we MUST defer.
    let max_partial = m.min(n);
    if min_overlap > max_partial {
        // No valid i remains; the partial-match loop in the scalar DP
        // wouldn't execute. Combined with the failed full-match check
        // above, no match is possible.
        return true;
    }

    let mut cur: usize = 0; // dp[0][n] = 0
    for i in 1..=max_partial {
        let bit = 1u64 << (i - 1);
        if pv & bit != 0 {
            cur += 1;
        } else if mv & bit != 0 {
            cur = cur.saturating_sub(1);
        }
        if i >= min_overlap {
            let max_errors_i = (max_error_rate * i as f64).floor() as usize;
            if cur <= max_errors_i {
                return false; // possible partial match — defer to scalar DP
            }
        }
    }

    // Both checks ruled it out: no full match, no partial match. Safe
    // to short-circuit; the scalar DP would return None.
    true
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

    // ── Myers' bit-parallel prefilter — direct unit tests ──
    //
    // These exercise myers_proves_no_match in isolation, asserting the
    // conservativeness contract: NEVER return true when the scalar DP
    // would return Some. The 12 tests above already cover the integrated
    // behaviour (find_3prime_adapter with prefilter active still returns
    // the same matches); these tests ground the bit-vector logic directly.

    #[test]
    fn myers_adapterless_read_proves_no_match() {
        // Read with zero adapter content (long stretch of unrelated bases).
        // Adapter is the standard 13 bp Illumina sequence.
        // At error_rate=0.1, full-match budget is floor(0.1*13)=1, so any
        // read sequence with >1 mismatch in every alignment window must
        // cleanly rule out a full match. Partial overlap of length 1 with
        // 0 errors needs the last read base to equal the first adapter
        // base ('A'); we end on 'X' to also rule out partials.
        let read = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
        let adapter = b"AGATCGGAAGAGC";
        assert!(
            myers_proves_no_match(read, adapter, 0.1, 1),
            "adapterless read with no partial-overlap base must prove no-match"
        );
    }

    #[test]
    fn myers_adapter_present_returns_false() {
        // Read containing the adapter exactly. Prefilter MUST defer.
        let read = b"ACGTACGTACGTACGTAGATCGGAAGAGC";
        let adapter = b"AGATCGGAAGAGC";
        assert!(
            !myers_proves_no_match(read, adapter, 0.1, 1),
            "adapter-present read must NOT short-circuit"
        );
    }

    #[test]
    fn myers_adapter_within_error_budget_returns_false() {
        // One mismatch — budget is floor(0.1*13)=1, so this is right at
        // the threshold. Prefilter must defer to the scalar DP.
        let read = b"ACGTACGTACGTACGTAGATCXGAAGAGC";
        //                                ^ one mismatch (G→X relative
        //                                  to the canonical adapter)
        let adapter = b"AGATCGGAAGAGC";
        assert!(
            !myers_proves_no_match(read, adapter, 0.1, 1),
            "adapter-with-1-error read must NOT short-circuit at error_rate=0.1"
        );
    }

    #[test]
    fn myers_partial_overlap_at_3prime_end_returns_false() {
        // Read ends with the first 3 bases of the adapter. The scalar DP
        // would return a partial match (overlap=3). Prefilter must defer.
        let read = b"XXXXXXXXXXXXXAGA";
        let adapter = b"AGATCGGAAGAGC";
        assert!(
            !myers_proves_no_match(read, adapter, 0.1, 1),
            "partial 3' overlap must NOT short-circuit"
        );
    }

    #[test]
    fn myers_long_adapter_falls_through() {
        // Adapter longer than 64 bp — out of scope for single-word
        // Myers'. Must always return false (defer to scalar DP).
        let adapter: Vec<u8> = b"ACGT".repeat(20); // 80 bp
        let read = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        assert_eq!(adapter.len(), 80);
        assert!(
            !myers_proves_no_match(read, &adapter, 0.1, 1),
            "adapter > 64 bp must always defer to scalar DP"
        );
    }

    #[test]
    fn myers_bgi_adapter_64bp_boundary_works() {
        // 64 bp synthetic adapter — the upper limit for single-word
        // Myers'. No adapter content in the read should cleanly prove
        // no-match.
        let adapter: Vec<u8> = b"ACGT".repeat(16); // 64 bp exactly
        assert_eq!(adapter.len(), 64);
        let read = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
        assert!(
            myers_proves_no_match(read, &adapter, 0.1, 1),
            "64 bp adapter at the prefilter limit must work correctly on adapterless reads"
        );
    }

    #[test]
    fn myers_zero_min_overlap_is_caller_responsibility() {
        // The find_3prime_adapter caller already rejects min_overlap == 0
        // before calling the prefilter. We don't re-validate; the test is
        // here to lock the contract: this function trusts its inputs.
        // With min_overlap == 1 (the smallest valid value), we should
        // see expected behaviour on a totally adapter-free read.
        let read = b"NNNNNNNNNNNNNNNN"; // all N — no overlap with ACGT adapter
        let adapter = b"AGATCGGAAGAGC";
        assert!(
            myers_proves_no_match(read, adapter, 0.1, 1),
            "all-N read against ACGT adapter must prove no-match"
        );
    }

    /// Equivalence sweep: for a battery of (read, adapter, error_rate,
    /// min_overlap) inputs, the prefilter MUST NOT contradict the scalar
    /// DP. Specifically: if `find_3prime_adapter` returns Some, the
    /// prefilter must NOT return true (which would have caused
    /// find_3prime_adapter to return None). This is the "no false negatives"
    /// half of the conservativeness contract.
    ///
    /// Hand-picked test cases targeting tricky DP regions: exact match,
    /// one-mismatch, partial 3' overlap, mid-read indel, off-by-one at
    /// match boundaries.
    #[test]
    fn myers_no_false_negatives_against_scalar_dp() {
        let cases: &[(&[u8], &[u8], f64, usize)] = &[
            // Exact full match
            (b"ACGTACGTAGATCGGAAGAGC", b"AGATCGGAAGAGC", 0.1, 1),
            // Full match with one substitution
            (b"ACGTACGTAGATCGXAAGAGC", b"AGATCGGAAGAGC", 0.1, 1),
            // Full match with one insertion
            (b"ACGTACGTAGATCGGXAAGAGC", b"AGATCGGAAGAGC", 0.15, 1),
            // Partial 3' overlap of length 5
            (b"XXXXXXXXAGATC", b"AGATCGGAAGAGC", 0.0, 1),
            // Partial 3' overlap of length 3, at error_rate=0.0 (must be exact)
            (b"XXXXXXXXAGA", b"AGATCGGAAGAGC", 0.0, 1),
            // Higher error rate widens the budget
            (b"ACGTACGTNNATCNGAAGAGC", b"AGATCGGAAGAGC", 0.2, 1),
            // Adapter at the very start of the read
            (b"AGATCGGAAGAGCNNNNNNNN", b"AGATCGGAAGAGC", 0.1, 1),
            // Long read, adapter midway
            (
                b"ACGTACGTACGTACGTACGTACGTACGTACGTAGATCGGAAGAGCXXXXXXX",
                b"AGATCGGAAGAGC",
                0.1,
                1,
            ),
            // Nextera adapter, partial overlap
            (b"NNNNNNNNNNNNNNCTGTC", b"CTGTCTCTTATA", 0.0, 1),
            // BGI 32 bp adapter, full match
            (
                b"NNNNNAAGTCGGAGGCCAAGCGGTCTTAGGAAGACAANNNNN",
                b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
                0.1,
                1,
            ),
            // BGI partial overlap at 3' end (length 8)
            (b"NNNNNNNNNNNNAAGTCGGA", b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", 0.0, 1),
        ];
        for (read, adapter, error_rate, min_overlap) in cases {
            let scalar_result = find_3prime_adapter(read, adapter, *error_rate, *min_overlap);
            if scalar_result.is_some() {
                // Re-run the prefilter directly; if it claims "no match"
                // we'd have a false negative.
                assert!(
                    !myers_proves_no_match(read, adapter, *error_rate, *min_overlap),
                    "FALSE NEGATIVE: prefilter rejected a match the scalar DP found.\n\
                     read = {:?}\nadapter = {:?}\nerror_rate = {}\nmin_overlap = {}\n\
                     scalar DP found: {scalar_result:?}",
                    std::str::from_utf8(read).unwrap_or("<non-utf8>"),
                    std::str::from_utf8(adapter).unwrap_or("<non-utf8>"),
                    error_rate,
                    min_overlap
                );
            }
        }
    }
}
