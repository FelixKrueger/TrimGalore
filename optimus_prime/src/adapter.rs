//! Adapter definitions, presets, and auto-detection.
//!
//! Provides built-in adapter sequences for common sequencing platforms
//! and auto-detection by scanning the first 1M reads.

use anyhow::Result;
use crate::fastq::FastqReader;
use std::path::Path;

/// Built-in adapter presets.
#[derive(Debug, Clone)]
pub struct AdapterPreset {
    /// Human-readable name
    pub name: &'static str,
    /// Primary adapter sequence (used for R1, or both if no R2 adapter)
    pub seq: &'static str,
    /// Optional R2-specific adapter sequence
    pub seq_r2: Option<&'static str>,
}

/// All built-in adapter presets.
pub const ILLUMINA: AdapterPreset = AdapterPreset {
    name: "Illumina",
    seq: "AGATCGGAAGAGC",
    seq_r2: None,
};

pub const NEXTERA: AdapterPreset = AdapterPreset {
    name: "Nextera",
    seq: "CTGTCTCTTATA",
    seq_r2: None,
};

pub const SMALL_RNA: AdapterPreset = AdapterPreset {
    name: "smallRNA",
    seq: "TGGAATTCTCGG",
    seq_r2: Some("GATCGTCGGACT"),
};

pub const STRANDED_ILLUMINA: AdapterPreset = AdapterPreset {
    name: "Stranded Illumina",
    seq: "ACTGTCTCTTATA",
    seq_r2: None,
};

pub const BGISEQ: AdapterPreset = AdapterPreset {
    name: "BGI/DNBSEQ",
    seq: "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    seq_r2: Some("AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"),
};

/// Result of adapter auto-detection.
#[derive(Debug)]
pub struct DetectionResult {
    /// The detected adapter preset
    pub adapter: AdapterPreset,
    /// Number of reads containing each adapter type
    pub counts: Vec<(String, usize)>,
    /// Total reads scanned
    pub reads_scanned: usize,
    /// Human-readable report message
    pub message: String,
}

/// Maximum number of reads to scan for auto-detection.
const MAX_SCAN_READS: usize = 1_000_000;

/// Auto-detect the adapter type by scanning the first file.
///
/// Scans up to 1M reads from the first input file, counting exact substring
/// matches of the three canonical adapter sequences. Returns the most common
/// adapter, defaulting to Illumina on ties.
///
/// This matches TrimGalore's `autodetect_adapter_type()` behavior, which
/// only scans `$ARGV[0]` (the first file, i.e., R1 for paired-end).
pub fn autodetect_adapter<P: AsRef<Path>>(path: P) -> Result<DetectionResult> {
    let mut reader = FastqReader::open(path)?;

    let adapters = [
        ("Illumina", ILLUMINA.seq),
        ("Nextera", NEXTERA.seq),
        ("smallRNA", SMALL_RNA.seq),
    ];

    let mut counts = [0usize; 3];
    let mut reads_scanned = 0;

    while let Some(record) = reader.next_record()? {
        reads_scanned += 1;
        if reads_scanned > MAX_SCAN_READS {
            break;
        }

        let seq = record.seq.as_bytes();
        for (i, (_, adapter_seq)) in adapters.iter().enumerate() {
            if contains_subsequence(seq, adapter_seq.as_bytes()) {
                counts[i] += 1;
            }
        }
    }

    // Build counts vec for reporting
    let counts_vec: Vec<(String, usize)> = adapters
        .iter()
        .zip(counts.iter())
        .map(|((name, _), &count)| (name.to_string(), count))
        .collect();

    // Find the winner — default to Illumina on ties
    let (winner_idx, &winner_count) = counts
        .iter()
        .enumerate()
        .max_by(|(i_a, a), (i_b, b)| {
            a.cmp(b).then_with(|| {
                // On tie, prefer Illumina (index 0)
                i_b.cmp(i_a)
            })
        })
        .unwrap();

    let (winner_name, _) = adapters[winner_idx];

    // Determine which preset to return
    let preset = match winner_name {
        "Illumina" => ILLUMINA,
        "Nextera" => NEXTERA,
        "smallRNA" => SMALL_RNA,
        _ => ILLUMINA,
    };

    // Build second-best for the message
    let mut sorted_counts: Vec<(usize, &str, usize)> = counts
        .iter()
        .enumerate()
        .map(|(i, &c)| (i, adapters[i].0, c))
        .collect();
    sorted_counts.sort_by(|a, b| b.2.cmp(&a.2));

    let message = if winner_count == 0 {
        format!(
            "No adapters detected in the first {} reads. Defaulting to Illumina adapter.",
            reads_scanned
        )
    } else {
        format!(
            "Using {} adapter for trimming (count: {}). \
             Second best hit was {} (count: {}).",
            sorted_counts[0].1,
            sorted_counts[0].2,
            sorted_counts[1].1,
            sorted_counts[1].2,
        )
    };

    Ok(DetectionResult {
        adapter: preset,
        counts: counts_vec,
        reads_scanned,
        message,
    })
}

/// Check if `haystack` contains `needle` as a subsequence (exact substring match).
/// Uses simple byte-by-byte search. For auto-detection on 1M reads this is fast enough.
fn contains_subsequence(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.len() > haystack.len() {
        return false;
    }
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contains_subsequence() {
        assert!(contains_subsequence(b"ACGTAAGAGATCGGAAGAGCTTTT", b"AGATCGGAAGAGC"));
        assert!(!contains_subsequence(b"ACGTACGT", b"AGATCGGAAGAGC"));
        assert!(contains_subsequence(b"AGATCGGAAGAGC", b"AGATCGGAAGAGC"));
        assert!(!contains_subsequence(b"", b"AGATC"));
    }

    #[test]
    fn test_adapter_presets_r2() {
        // Small RNA and BGI should have R2 adapters
        assert!(SMALL_RNA.seq_r2.is_some());
        assert!(BGISEQ.seq_r2.is_some());
        // Illumina, Nextera, Stranded should not
        assert!(ILLUMINA.seq_r2.is_none());
        assert!(NEXTERA.seq_r2.is_none());
        assert!(STRANDED_ILLUMINA.seq_r2.is_none());
    }
}
