//! Read filtering by length, N-content, and max-length.

use crate::fastq::FastqRecord;

/// Filter result for a single read.
#[derive(Debug, PartialEq)]
pub enum FilterResult {
    /// Read passes all filters
    Pass,
    /// Read is too short (below length_cutoff)
    TooShort,
    /// Read is too long (above max_length)
    TooLong,
    /// Read has too many N bases
    TooManyN,
}

/// Check if a read passes length and N-content filters.
pub fn filter_single_end(
    record: &FastqRecord,
    length_cutoff: usize,
    max_length: Option<usize>,
    max_n: Option<MaxNFilter>,
) -> FilterResult {
    // N-content check first (matches TrimGalore order)
    if let Some(ref max_n_filter) = max_n
        && exceeds_n_threshold(record, max_n_filter)
    {
        return FilterResult::TooManyN;
    }

    // Length checks
    if record.len() < length_cutoff {
        return FilterResult::TooShort;
    }

    if let Some(max) = max_length
        && record.len() > max
    {
        return FilterResult::TooLong;
    }

    FilterResult::Pass
}

/// Filter result for a read pair.
#[derive(Debug, PartialEq)]
pub enum PairFilterResult {
    /// Both reads pass
    Pass,
    /// Pair discarded due to N-content (no unpaired rescue possible)
    TooManyN,
    /// One or both reads too short. Contains (r1_rescuable, r2_rescuable)
    TooShort { r1_ok: bool, r2_ok: bool },
    /// One or both reads too long (no unpaired rescue possible)
    TooLong,
}

/// Per-side minimum length thresholds for unpaired rescue under
/// `--retain_unpaired`. When a pair fails the length check, a read survives
/// on its own if its length meets the corresponding threshold.
#[derive(Debug, Clone, Copy)]
pub struct UnpairedLengths {
    pub r1: usize,
    pub r2: usize,
}

/// Check if a read pair passes filters.
///
/// Key asymmetry: N-content failure always discards the entire pair (no
/// unpaired rescue), but length failure can rescue individual reads via
/// --retain_unpaired.
pub fn filter_paired_end(
    r1: &FastqRecord,
    r2: &FastqRecord,
    length_cutoff: usize,
    max_length: Option<usize>,
    max_n: Option<MaxNFilter>,
    unpaired: UnpairedLengths,
) -> PairFilterResult {
    // N-content check — ALWAYS discards entire pair, no rescue
    if let Some(ref max_n_filter) = max_n
        && (exceeds_n_threshold(r1, max_n_filter) || exceeds_n_threshold(r2, max_n_filter))
    {
        return PairFilterResult::TooManyN;
    }

    // Length check — can rescue individual reads with --retain_unpaired.
    //
    // Perl-parity note (matches `0.6.11:trim_galore:2325-2343`):
    // when at least one mate fails the discard cutoff, each read is
    // independently routed to its unpaired file based solely on its own
    // length vs the per-side cutoff (`--length_1` / `--length_2`). There
    // is no requirement that the OTHER mate must have failed; both sides
    // can land in unpaired in the same record. Rust v2.1.0-beta.5 had an
    // extra `!r{1,2}_short` clause here that gated unpaired-rescue on
    // the read passing the discard cutoff, which silently dropped reads
    // when both mates failed `--length` but were individually long
    // enough for `--length_{1,2}`. See #245 (parity-hunt).
    let r1_short = r1.len() < length_cutoff;
    let r2_short = r2.len() < length_cutoff;

    if r1_short || r2_short {
        return PairFilterResult::TooShort {
            r1_ok: r1.len() >= unpaired.r1,
            r2_ok: r2.len() >= unpaired.r2,
        };
    }

    // Max-length check — discards entire pair
    if let Some(max) = max_length
        && (r1.len() > max || r2.len() > max)
    {
        return PairFilterResult::TooLong;
    }

    PairFilterResult::Pass
}

/// Configuration for N-content filtering.
#[derive(Debug, Clone)]
pub enum MaxNFilter {
    /// Absolute count: discard if N count > threshold
    Count(usize),
    /// Fraction of read length: discard if N fraction > threshold
    Fraction(f64),
}

/// Check if a record exceeds the N-content threshold.
fn exceeds_n_threshold(record: &FastqRecord, filter: &MaxNFilter) -> bool {
    if record.is_empty() {
        return true; // 0-length reads are filtered out
    }

    let n_count = record.n_count();

    match filter {
        MaxNFilter::Count(max) => n_count > *max,
        MaxNFilter::Fraction(max_frac) => {
            let frac = n_count as f64 / record.len() as f64;
            frac > *max_frac
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(seq: &str) -> FastqRecord {
        FastqRecord {
            id: "@test".to_string(),
            seq: seq.to_string(),
            qual: "I".repeat(seq.len()),
        }
    }

    #[test]
    fn test_filter_pass() {
        let rec = make_record("ACGTACGT");
        assert_eq!(filter_single_end(&rec, 5, None, None), FilterResult::Pass);
    }

    #[test]
    fn test_filter_too_short() {
        let rec = make_record("ACG");
        assert_eq!(
            filter_single_end(&rec, 5, None, None),
            FilterResult::TooShort
        );
    }

    #[test]
    fn test_filter_too_long() {
        let rec = make_record("ACGTACGTACGT");
        assert_eq!(
            filter_single_end(&rec, 5, Some(10), None),
            FilterResult::TooLong
        );
    }

    #[test]
    fn test_filter_too_many_n_count() {
        let rec = make_record("ACNNNGT");
        assert_eq!(
            filter_single_end(&rec, 1, None, Some(MaxNFilter::Count(2))),
            FilterResult::TooManyN
        );
    }

    #[test]
    fn test_filter_too_many_n_fraction() {
        let rec = make_record("NNNNACGT"); // 4/8 = 0.5
        assert_eq!(
            filter_single_end(&rec, 1, None, Some(MaxNFilter::Fraction(0.3))),
            FilterResult::TooManyN
        );
    }

    #[test]
    fn test_paired_n_filter_no_rescue() {
        let r1 = make_record("ACGTACGT");
        let r2 = make_record("NNNNNNNN"); // All N
        let result = filter_paired_end(
            &r1,
            &r2,
            5,
            None,
            Some(MaxNFilter::Count(2)),
            UnpairedLengths { r1: 35, r2: 35 },
        );
        assert_eq!(result, PairFilterResult::TooManyN);
    }

    #[test]
    fn test_paired_length_rescue() {
        let r1 = make_record("ACGTACGTACGT"); // 12bp, long enough
        let r2 = make_record("AC"); // 2bp, too short
        let result = filter_paired_end(&r1, &r2, 5, None, None, UnpairedLengths { r1: 10, r2: 10 });
        match result {
            PairFilterResult::TooShort { r1_ok, r2_ok } => {
                assert!(r1_ok); // R1 is rescuable (>= unpaired.r1)
                assert!(!r2_ok); // R2 is too short even for unpaired
            }
            _ => panic!("Expected TooShort"),
        }
    }

    /// Perl parity (0.6.11:trim_galore:2325-2343): when BOTH mates fail
    /// the discard `--length` cutoff but each is individually long enough
    /// for the per-side `--length_{1,2}` threshold, both reads route to
    /// their unpaired files. Regression for #245 — Rust previously had a
    /// `!r{1,2}_short` clause that gated unpaired rescue on the read
    /// itself passing the discard cutoff, which silently dropped reads
    /// when both mates failed `--length` together.
    #[test]
    fn test_paired_both_mates_fail_discard_both_routable_to_unpaired() {
        // 50 bp reads, --length 80 (both fail), --length_1 = --length_2 = 35
        // (default per-side); each individually >= 35 so both should rescue.
        let seq = "A".repeat(50);
        let r1 = make_record(&seq);
        let r2 = make_record(&seq);
        let result =
            filter_paired_end(&r1, &r2, 80, None, None, UnpairedLengths { r1: 35, r2: 35 });
        match result {
            PairFilterResult::TooShort { r1_ok, r2_ok } => {
                assert!(r1_ok, "R1 (50 bp >= 35) should rescue to unpaired");
                assert!(r2_ok, "R2 (50 bp >= 35) should rescue to unpaired");
            }
            other => panic!("Expected TooShort, got {:?}", other),
        }
    }
}
