//! Adapter definitions, presets, and auto-detection.
//!
//! Provides built-in adapter sequences for common sequencing platforms
//! and auto-detection by scanning the first 1M reads.

use anyhow::{Context, Result};
use crate::fastq::FastqReader;
use std::io::BufRead;
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
    /// Whether adapter trimming was suppressed (--consider_already_trimmed)
    pub suppressed: bool,
    /// Number of reads with ≥10 consecutive trailing G's (for poly-G auto-detection)
    pub poly_g_count: usize,
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
pub fn autodetect_adapter<P: AsRef<Path>>(
    path: P,
    consider_already_trimmed: Option<usize>,
) -> Result<DetectionResult> {
    let mut reader = FastqReader::open(path)?;

    let adapters = [
        ("Illumina", ILLUMINA.seq),
        ("Nextera", NEXTERA.seq),
        ("smallRNA", SMALL_RNA.seq),
    ];

    let mut counts = [0usize; 3];
    let mut reads_scanned = 0;
    let mut poly_g_count: usize = 0;

    while let Some(record) = reader.next_record()? {
        reads_scanned += 1;
        if reads_scanned >= MAX_SCAN_READS {
            break;
        }

        let seq = record.seq.as_bytes();
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

    // Check if adapter trimming should be suppressed
    let suppressed = if let Some(threshold) = consider_already_trimmed {
        winner_count <= threshold
    } else {
        false
    };

    let message = if suppressed {
        format!(
            "No auto-detected adapter sequence exceeded the user-specified \
             'already adapter-trimmed' limit of {} counts. \
             Only quality trimming will be carried out.",
            consider_already_trimmed.unwrap()
        )
    } else if winner_count == 0 {
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

    // If suppressed, return an empty adapter to skip adapter trimming
    let final_preset = if suppressed {
        AdapterPreset {
            name: "already trimmed (adapter trimming suppressed)",
            seq: "",
            seq_r2: None,
        }
    } else {
        preset
    };

    Ok(DetectionResult {
        adapter: final_preset,
        counts: counts_vec,
        reads_scanned,
        message,
        suppressed,
        poly_g_count,
    })
}

/// Check if a sequence ends with `min_len` or more consecutive G bases.
/// Used for poly-G auto-detection (strict: no mismatches allowed).
fn has_trailing_poly_g(seq: &[u8], min_len: usize) -> bool {
    if seq.len() < min_len {
        return false;
    }
    let mut count = 0;
    for &base in seq.iter().rev() {
        if base == b'G' {
            count += 1;
            if count >= min_len {
                return true;
            }
        } else {
            break;
        }
    }
    false
}

/// Standalone poly-G detection scan (used when adapter auto-detection is skipped).
///
/// Scans up to 1M reads from the input file, counting reads whose 3' end
/// has ≥10 consecutive G bases. Returns (poly_g_count, reads_scanned).
///
/// This is called when the user specifies an adapter explicitly (--illumina,
/// --adapter, etc.), which causes `autodetect_adapter()` to be skipped. Without
/// this, poly-G auto-detection would silently fail for those users.
pub fn detect_poly_g<P: AsRef<Path>>(path: P) -> Result<(usize, usize)> {
    let mut reader = FastqReader::open(path)?;
    let mut reads_scanned = 0;
    let mut poly_g_count: usize = 0;

    while let Some(record) = reader.next_record()? {
        reads_scanned += 1;
        if has_trailing_poly_g(record.seq.as_bytes(), 10) {
            poly_g_count += 1;
        }
        if reads_scanned >= MAX_SCAN_READS {
            break;
        }
    }

    Ok((poly_g_count, reads_scanned))
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

/// Parse a raw adapter CLI string into a list of named adapter sequences.
///
/// Supported formats:
///   - Plain sequence: `"AGATCGGAAGAGC"` → `[("adapter_1", "AGATCGGAAGAGC")]`
///   - Embedded multi: `" SEQ1 -a SEQ2 -a SEQ3"` → 3 entries (leading space + `-a` separators)
///   - FASTA file: `"file:path/to/adapters.fa"` → entries from FASTA
pub fn parse_adapter_spec(raw: &str) -> Result<Vec<(String, String)>> {
    // Case 1: FASTA file reference
    if let Some(path) = raw.strip_prefix("file:") {
        return read_fasta_adapters(path.trim());
    }

    // Case 2: Embedded multi-adapter (contains " -a " separator)
    if raw.contains(" -a ") {
        let parts: Vec<&str> = raw.split(" -a ").collect();
        let mut adapters = Vec::with_capacity(parts.len());
        for (_i, part) in parts.iter().enumerate() {
            let seq = part.trim().to_uppercase();
            if seq.is_empty() {
                continue;
            }
            validate_adapter_sequence(&seq)?;
            adapters.push((format!("adapter_{}", adapters.len() + 1), seq));
        }
        if adapters.is_empty() {
            anyhow::bail!("No valid adapter sequences found in multi-adapter specification");
        }
        return Ok(adapters);
    }

    // Case 3: Single adapter sequence
    let seq = raw.trim().to_uppercase();
    if seq.is_empty() {
        anyhow::bail!("Empty adapter sequence");
    }
    validate_adapter_sequence(&seq)?;
    Ok(vec![("adapter_1".to_string(), seq)])
}

/// Read adapter sequences from a FASTA file.
///
/// Returns `Vec<(name, sequence)>` where name is the FASTA header (without `>`).
/// Sequences are uppercased and whitespace-stripped.
pub fn read_fasta_adapters<P: AsRef<Path>>(path: P) -> Result<Vec<(String, String)>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open adapter FASTA file: {}", path.display()))?;
    let reader = std::io::BufReader::new(file);

    let mut adapters: Vec<(String, String)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            // Flush previous entry
            if let Some(name) = current_name.take() {
                if current_seq.is_empty() {
                    anyhow::bail!("Empty sequence for adapter '{}' in {}", name, path.display());
                }
                validate_adapter_sequence(&current_seq)?;
                adapters.push((name, current_seq.clone()));
                current_seq.clear();
            }
            current_name = Some(header.trim().to_string());
        } else {
            current_seq.push_str(&line.to_uppercase());
        }
    }

    // Flush last entry
    if let Some(name) = current_name {
        if current_seq.is_empty() {
            anyhow::bail!("Empty sequence for adapter '{}' in {}", name, path.display());
        }
        validate_adapter_sequence(&current_seq)?;
        adapters.push((name, current_seq));
    }

    if adapters.is_empty() {
        anyhow::bail!("No adapter sequences found in FASTA file: {}", path.display());
    }

    Ok(adapters)
}

/// Validate that an adapter sequence contains only valid DNA characters.
fn validate_adapter_sequence(seq: &str) -> Result<()> {
    if !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'X')) {
        anyhow::bail!(
            "Adapter sequence must contain only DNA characters (A, C, G, T, N, X), got: '{}'",
            seq
        );
    }
    Ok(())
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
    fn test_has_trailing_poly_g() {
        // Clear poly-G tail
        assert!(has_trailing_poly_g(b"ACGTACGTGGGGGGGGGG", 10));
        // Exactly at threshold
        assert!(has_trailing_poly_g(b"ACGTGGGGGGGGGG", 10));
        // Below threshold
        assert!(!has_trailing_poly_g(b"ACGTGGGGGGGGG", 10)); // 9 Gs
        // No trailing Gs
        assert!(!has_trailing_poly_g(b"ACGTACGTACGT", 10));
        // Entirely poly-G
        assert!(has_trailing_poly_g(b"GGGGGGGGGGGGGGGG", 10));
        // Interrupted (non-G breaks the run)
        assert!(!has_trailing_poly_g(b"ACGTGGGGGAGGGGGGGG", 10)); // only 9 trailing
        // Short sequence
        assert!(!has_trailing_poly_g(b"GGG", 10));
        // Empty
        assert!(!has_trailing_poly_g(b"", 10));
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

    #[test]
    fn test_parse_single_adapter() {
        let result = parse_adapter_spec("AGATCGGAAGAGC").unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, "AGATCGGAAGAGC");
    }

    #[test]
    fn test_parse_single_adapter_lowercase() {
        let result = parse_adapter_spec("agatcggaagagc").unwrap();
        assert_eq!(result[0].1, "AGATCGGAAGAGC");
    }

    #[test]
    fn test_parse_multi_adapter_embedded() {
        let result = parse_adapter_spec(" AGCTAGCG -a TCTCTTATAT -a TTTATTCGGATTTAT").unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], ("adapter_1".to_string(), "AGCTAGCG".to_string()));
        assert_eq!(result[1], ("adapter_2".to_string(), "TCTCTTATAT".to_string()));
        assert_eq!(result[2], ("adapter_3".to_string(), "TTTATTCGGATTTAT".to_string()));
    }

    #[test]
    fn test_parse_fasta_adapter() {
        use std::io::Write;
        let dir = std::env::temp_dir().join("tg_test_fasta");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("adapters.fa");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">Illumina\nAGATCGGAAGAGC").unwrap();
        writeln!(f, ">Nextera\nCTGTCTCTTATA").unwrap();

        let result = read_fasta_adapters(&path).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], ("Illumina".to_string(), "AGATCGGAAGAGC".to_string()));
        assert_eq!(result[1], ("Nextera".to_string(), "CTGTCTCTTATA".to_string()));

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn test_parse_fasta_invalid_bases() {
        use std::io::Write;
        let dir = std::env::temp_dir().join("tg_test_fasta_invalid");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bad.fa");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">bad\nAGCTZZZ").unwrap();

        let result = read_fasta_adapters(&path);
        assert!(result.is_err());

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn test_parse_file_not_found() {
        let result = parse_adapter_spec("file:/nonexistent/adapters.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_adapter_with_file_prefix() {
        // Should dispatch to FASTA reader, which will fail on missing file
        let result = parse_adapter_spec("file:no_such_file.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_adapter_sequence_valid() {
        assert!(validate_adapter_sequence("ACGTNX").is_ok());
    }

    #[test]
    fn test_validate_adapter_sequence_invalid() {
        assert!(validate_adapter_sequence("ACGTZ").is_err());
    }
}
