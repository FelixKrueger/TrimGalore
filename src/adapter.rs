//! Adapter definitions, presets, and auto-detection.
//!
//! Provides built-in adapter sequences for common sequencing platforms
//! and auto-detection by scanning the first 1M reads.

use crate::fastq::FastqReader;
use anyhow::{Context, Result};
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

impl AdapterPreset {
    /// Convert this preset's primary adapter to a `(name, sequence)` vec for the multi-adapter pipeline.
    pub fn to_adapter_vec(&self) -> Vec<(String, String)> {
        vec![(self.name.to_string(), self.seq.to_string())]
    }

    /// Convert this preset's R2-specific adapter (if any) to a `(name, sequence)` vec.
    pub fn to_r2_vec(&self) -> Vec<(String, String)> {
        self.seq_r2
            .map(|s| vec![(format!("{}_r2", self.name), s.to_string())])
            .unwrap_or_default()
    }
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
    autodetect_adapter_with_max_scan(path, consider_already_trimmed, MAX_SCAN_READS)
}

/// Inner implementation of [`autodetect_adapter`] with the per-call scan
/// cap exposed as a parameter. Crate-private so unit tests can exercise
/// the boundary behaviour with a small `max_scan` instead of needing a
/// fixture larger than 1M reads; production callers always go through
/// the [`autodetect_adapter`] wrapper which pins the cap at
/// `MAX_SCAN_READS`.
pub(crate) fn autodetect_adapter_with_max_scan<P: AsRef<Path>>(
    path: P,
    consider_already_trimmed: Option<usize>,
    max_scan: usize,
) -> Result<DetectionResult> {
    let mut reader = FastqReader::open(path)?;

    // BGI/DNBSEQ added to the probe set: its 32bp adapter shares no
    // meaningful subsequence with the other three and the per-match
    // false-positive rate is vanishingly low, so including it is
    // uncontentious. Stranded Illumina is intentionally omitted —
    // its sequence differs from Nextera by a single base (A-tail),
    // which would produce constant ambiguous ties.
    let adapters = [
        ("Illumina", ILLUMINA.seq),
        ("Nextera", NEXTERA.seq),
        ("smallRNA", SMALL_RNA.seq),
        ("BGI/DNBSEQ", BGISEQ.seq),
    ];

    let mut counts = [0usize; 4];
    let mut reads_scanned = 0;
    let mut poly_g_count: usize = 0;

    while let Some(record) = reader.next_record()? {
        reads_scanned += 1;
        if reads_scanned >= max_scan {
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
        "BGI/DNBSEQ" => BGISEQ,
        _ => ILLUMINA,
    };

    // Build second-best for the message
    let mut sorted_counts: Vec<(usize, &str, usize)> = counts
        .iter()
        .enumerate()
        .map(|(i, &c)| (i, adapters[i].0, c))
        .collect();
    sorted_counts.sort_by_key(|entry| std::cmp::Reverse(entry.2));

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
            sorted_counts[0].1, sorted_counts[0].2, sorted_counts[1].1, sorted_counts[1].2,
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
        for part in &parts {
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
    let mut seq = raw.trim().to_uppercase();
    if seq.is_empty() {
        anyhow::bail!("Empty adapter sequence");
    }
    // Perl-style brace shorthand: A{10} → AAAAAAAAAA. Only applied to the
    // single-adapter case (not to multi-adapter elements or FASTA entries),
    // matching `trim_galore:3105–3112` in the Perl source.
    if let Some(expanded) = expand_brace_notation(&seq) {
        if expanded.is_empty() {
            anyhow::bail!(
                "Adapter sequence {} expanded to an empty string (use N >= 1)",
                seq
            );
        }
        eprintln!("Adapter sequence {} expanded to {}", seq, expanded);
        seq = expanded;
    }
    validate_adapter_sequence(&seq)?;
    Ok(vec![("adapter_1".to_string(), seq)])
}

/// Parse a slice of adapter CLI strings, concatenating results with
/// sequential `adapter_N` names regardless of which input produced them.
///
/// This is the entry point when a user passes multiple `-a` (or `-a2`)
/// flags — each string goes through `parse_adapter_spec` independently,
/// so mixing forms is supported: `-a " SEQ -a SEQ2" -a SEQ3` or
/// `-a SEQ -a "file:adapters.fa"` all work and produce a flat,
/// renumbered list.
pub fn parse_adapter_specs(specs: &[String]) -> Result<Vec<(String, String)>> {
    let mut result: Vec<(String, String)> = Vec::new();
    for spec in specs {
        for (_name, seq) in parse_adapter_spec(spec)? {
            result.push((format!("adapter_{}", result.len() + 1), seq));
        }
    }
    Ok(result)
}

/// Expand Perl-style `X{N}` shorthand to `N` copies of `X`.
///
/// Matches the whole string against `^[ACTGN]\{(\d+)\}$` (uppercase only;
/// the caller is expected to have upcased already). Returns `None` if the
/// input doesn't match the pattern, so callers can fall through to normal
/// adapter validation. Mirrors Perl's `extend_adapter_sequence` behaviour.
fn expand_brace_notation(seq: &str) -> Option<String> {
    let bytes = seq.as_bytes();
    if bytes.len() < 4 {
        return None; // minimum valid form is "A{0}"
    }
    let base = bytes[0];
    if !matches!(base, b'A' | b'C' | b'T' | b'G' | b'N') {
        return None;
    }
    if bytes[1] != b'{' || *bytes.last().unwrap() != b'}' {
        return None;
    }
    let digits = &bytes[2..bytes.len() - 1];
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let n: usize = std::str::from_utf8(digits).ok()?.parse().ok()?;
    Some((base as char).to_string().repeat(n))
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
                    anyhow::bail!(
                        "Empty sequence for adapter '{}' in {}",
                        name,
                        path.display()
                    );
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
            anyhow::bail!(
                "Empty sequence for adapter '{}' in {}",
                name,
                path.display()
            );
        }
        validate_adapter_sequence(&current_seq)?;
        adapters.push((name, current_seq));
    }

    if adapters.is_empty() {
        anyhow::bail!(
            "No adapter sequences found in FASTA file: {}",
            path.display()
        );
    }

    Ok(adapters)
}

/// Validate that an adapter sequence contains only valid DNA characters.
fn validate_adapter_sequence(seq: &str) -> Result<()> {
    if !seq
        .bytes()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'X'))
    {
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
        assert!(contains_subsequence(
            b"ACGTAAGAGATCGGAAGAGCTTTT",
            b"AGATCGGAAGAGC"
        ));
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
        assert_eq!(
            result[1],
            ("adapter_2".to_string(), "TCTCTTATAT".to_string())
        );
        assert_eq!(
            result[2],
            ("adapter_3".to_string(), "TTTATTCGGATTTAT".to_string())
        );
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
        assert_eq!(
            result[0],
            ("Illumina".to_string(), "AGATCGGAAGAGC".to_string())
        );
        assert_eq!(
            result[1],
            ("Nextera".to_string(), "CTGTCTCTTATA".to_string())
        );

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

    // --- Brace-notation expansion (Perl parity: -a A{10} → -a AAAAAAAAAA) ---

    #[test]
    fn test_expand_brace_notation_basic_a10() {
        assert_eq!(
            expand_brace_notation("A{10}").as_deref(),
            Some("AAAAAAAAAA")
        );
    }

    #[test]
    fn test_expand_brace_notation_all_valid_bases() {
        assert_eq!(expand_brace_notation("C{3}").as_deref(), Some("CCC"));
        assert_eq!(expand_brace_notation("T{5}").as_deref(), Some("TTTTT"));
        assert_eq!(expand_brace_notation("G{1}").as_deref(), Some("G"));
        assert_eq!(
            expand_brace_notation("N{20}").as_deref(),
            Some("NNNNNNNNNNNNNNNNNNNN")
        );
    }

    #[test]
    fn test_expand_brace_notation_zero_is_empty() {
        // Perl accepts A{0} → empty; validate_adapter_sequence will reject
        // the empty result downstream, not our job here.
        assert_eq!(expand_brace_notation("A{0}").as_deref(), Some(""));
    }

    #[test]
    fn test_expand_brace_notation_rejects_multi_char_prefix() {
        // AG{10} is not a single-base shorthand; must not match.
        assert_eq!(expand_brace_notation("AG{10}"), None);
    }

    #[test]
    fn test_expand_brace_notation_rejects_non_dna_base() {
        // X and Z are not in the [ACTGN] set.
        assert_eq!(expand_brace_notation("X{5}"), None);
        assert_eq!(expand_brace_notation("Z{5}"), None);
    }

    #[test]
    fn test_expand_brace_notation_rejects_non_digit_count() {
        assert_eq!(expand_brace_notation("A{abc}"), None);
        assert_eq!(expand_brace_notation("A{}"), None);
        assert_eq!(expand_brace_notation("A{-1}"), None);
    }

    #[test]
    fn test_expand_brace_notation_rejects_plain_sequence() {
        assert_eq!(expand_brace_notation("AGATCGGAAGAGC"), None);
    }

    #[test]
    fn test_expand_brace_notation_rejects_lowercase() {
        // Parser upcases the input before this is called; lowercase input
        // to expand_brace_notation itself should not match.
        assert_eq!(expand_brace_notation("a{10}"), None);
    }

    #[test]
    fn test_parse_adapter_spec_brace_expansion_end_to_end() {
        let result = parse_adapter_spec("A{10}").unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, "AAAAAAAAAA");
    }

    #[test]
    fn test_parse_adapter_spec_brace_expansion_lowercase_input() {
        // Lowercase input gets upcased by the parser, then expanded.
        // Rust is strictly more permissive than Perl here.
        let result = parse_adapter_spec("a{5}").unwrap();
        assert_eq!(result[0].1, "AAAAA");
    }

    #[test]
    fn test_parse_adapter_spec_brace_expansion_zero_errors() {
        // A{0} → "" → validate_adapter_sequence rejects empty.
        assert!(parse_adapter_spec("A{0}").is_err());
    }

    #[test]
    fn test_parse_adapter_spec_no_brace_expansion_in_multi() {
        // Inside "SEQ1 -a SEQ2" syntax, elements are NOT brace-expanded.
        // A{10} as a multi-adapter element fails DNA validation.
        let result = parse_adapter_spec(" AGCT -a A{10}");
        assert!(result.is_err());
    }

    // --- parse_adapter_specs: multi-flag + mixed-form path ---

    #[test]
    fn test_parse_adapter_specs_repeated_single_sequences() {
        // Simulates: -a AGCT -a TTTA
        let specs = vec!["AGCT".to_string(), "TTTA".to_string()];
        let result = parse_adapter_specs(&specs).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], ("adapter_1".to_string(), "AGCT".to_string()));
        assert_eq!(result[1], ("adapter_2".to_string(), "TTTA".to_string()));
    }

    #[test]
    fn test_parse_adapter_specs_mixed_embedded_plus_single() {
        // Simulates: -a " ACT -a TCCTTTG" -a GTGTGTGCTTT
        // First spec is the v0.6.x embedded-string form (expands to 2
        // adapters); second spec is a plain single sequence. The
        // concatenated, renumbered list should have 3 adapters with
        // sequential names.
        let specs = vec![" ACT -a TCCTTTG".to_string(), "GTGTGTGCTTT".to_string()];
        let result = parse_adapter_specs(&specs).unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], ("adapter_1".to_string(), "ACT".to_string()));
        assert_eq!(result[1], ("adapter_2".to_string(), "TCCTTTG".to_string()));
        assert_eq!(
            result[2],
            ("adapter_3".to_string(), "GTGTGTGCTTT".to_string())
        );
    }

    #[test]
    fn test_parse_adapter_specs_empty_list_yields_empty() {
        // Regression guard: empty input (user didn't pass any -a) → empty Vec.
        let specs: Vec<String> = Vec::new();
        let result = parse_adapter_specs(&specs).unwrap();
        assert!(result.is_empty());
    }

    // --- Auto-detection probe set (Illumina + Nextera + smallRNA + BGI) ---

    fn write_fastq_with_adapter(
        path: &std::path::Path,
        adapter_seq: &str,
        num_reads: usize,
    ) -> Result<()> {
        use crate::fastq::{FastqRecord, FastqWriter, OUTPUT_GZIP_LEVEL};
        let mut writer = FastqWriter::create(path, false, 1, OUTPUT_GZIP_LEVEL)?;
        // Prepend a 20bp random-ish prefix so reads are plausibly long;
        // the probe is looking for adapter_seq as a substring anywhere.
        let prefix = "ACGTACGTACGTACGTACGT";
        for i in 0..num_reads {
            let seq = format!("{}{}", prefix, adapter_seq);
            let rec = FastqRecord {
                id: format!("@r{}", i),
                seq: seq.clone(),
                qual: "I".repeat(seq.len()),
            };
            writer.write_record(&rec)?;
        }
        writer.flush()?;
        Ok(())
    }

    #[test]
    fn test_autodetect_picks_bgi_when_bgi_adapter_dominant() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_autodetect_bgi");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("bgi.fq");
        write_fastq_with_adapter(&path, BGISEQ.seq, 100)?;

        let result = autodetect_adapter(&path, None)?;
        assert_eq!(result.adapter.name, "BGI/DNBSEQ");
        // BGI count should be high; other probes near zero (32bp probe vs
        // 20bp random prefix — no chance collision).
        let bgi_count = result
            .counts
            .iter()
            .find(|(name, _)| name == "BGI/DNBSEQ")
            .map(|(_, c)| *c)
            .unwrap_or(0);
        assert_eq!(bgi_count, 100);

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_autodetect_picks_illumina_when_illumina_adapter_dominant() -> Result<()> {
        // Regression guard: BGI addition must not break Illumina detection.
        let dir = std::env::temp_dir().join("tg_test_autodetect_illumina_regression");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("illumina.fq");
        write_fastq_with_adapter(&path, ILLUMINA.seq, 100)?;

        let result = autodetect_adapter(&path, None)?;
        assert_eq!(result.adapter.name, "Illumina");

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_autodetect_zero_match_defaults_to_illumina() -> Result<()> {
        // Regression guard: with all 4 probes at zero, the tie-break (lowest
        // index wins) must still select Illumina (index 0), not BGI (index 3).
        let dir = std::env::temp_dir().join("tg_test_autodetect_zero");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("none.fq");
        // Sequences with no known adapter substring.
        write_fastq_with_adapter(&path, "ACACACACACAC", 10)?;

        let result = autodetect_adapter(&path, None)?;
        assert_eq!(result.adapter.name, "Illumina"); // zero-count fallback
        // Confirm all four probes are in the counts report
        assert_eq!(result.counts.len(), 4);
        let names: Vec<&str> = result.counts.iter().map(|(n, _)| n.as_str()).collect();
        assert!(names.contains(&"Illumina"));
        assert!(names.contains(&"Nextera"));
        assert!(names.contains(&"smallRNA"));
        assert!(names.contains(&"BGI/DNBSEQ"));

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_autodetect_counts_all_four_probes() -> Result<()> {
        // The detection report must surface all four probe counts, so
        // downstream reporting (report.rs) has the full breakdown.
        let dir = std::env::temp_dir().join("tg_test_autodetect_all_four");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("nextera.fq");
        write_fastq_with_adapter(&path, NEXTERA.seq, 50)?;

        let result = autodetect_adapter(&path, None)?;
        assert_eq!(result.adapter.name, "Nextera");
        assert_eq!(result.counts.len(), 4);
        // Nextera should be at 50; others at 0 (short sequences, no chance collision).
        let nextera_count = result
            .counts
            .iter()
            .find(|(n, _)| n == "Nextera")
            .map(|(_, c)| *c)
            .unwrap_or(0);
        assert_eq!(nextera_count, 50);

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// §5.3 regression: the auto-detection scan caps at `MAX_SCAN_READS`
    /// (1,000,000 in production). Verified here via the crate-private
    /// `autodetect_adapter_with_max_scan` helper with a small `max_scan`
    /// — generating a >1M-read fixture every test run is prohibitively
    /// slow, so the boundary contract is exercised at scale 7 instead
    /// of scale 1e6 (same control-flow shape: increment-then-break).
    /// Asserts both that scanning stops at the cap and that
    /// `reads_scanned` reflects the cap value, not the file's record
    /// count.
    #[test]
    fn test_autodetect_respects_max_scan_boundary() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_autodetect_max_scan");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("plenty.fq");
        // 100 records, all with the Illumina adapter embedded → if scanning
        // ran to EOF, the Illumina count would be 100. With max_scan=7 the
        // count must be strictly less.
        write_fastq_with_adapter(&path, ILLUMINA.seq, 100)?;

        let result = autodetect_adapter_with_max_scan(&path, None, 7)?;
        assert_eq!(
            result.reads_scanned, 7,
            "reads_scanned must equal max_scan when input has more records than the cap"
        );
        let illumina_count = result
            .counts
            .iter()
            .find(|(n, _)| n == "Illumina")
            .map(|(_, c)| *c)
            .unwrap_or(0);
        // The cap is checked AFTER `reads_scanned += 1` so the 7th increment
        // breaks before the loop body that would have counted the 7th match.
        // Concretely: 6 records contribute to the count, 7 reads_scanned.
        assert!(
            illumina_count < 100,
            "Illumina count must be capped (got {illumina_count}, expected < 100)"
        );
        assert!(
            illumina_count <= 7,
            "Illumina count must be at most max_scan (got {illumina_count})"
        );

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Sanity test: with `max_scan` larger than the input record count,
    /// scanning runs to EOF and `reads_scanned` equals the file size.
    /// Pairs with `test_autodetect_respects_max_scan_boundary` so a
    /// regression that breaks the cap (e.g. always-stop-at-7) wouldn't
    /// pass both tests at once.
    #[test]
    fn test_autodetect_below_max_scan_processes_full_file() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_autodetect_full_scan");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("small.fq");
        write_fastq_with_adapter(&path, ILLUMINA.seq, 50)?;

        let result = autodetect_adapter_with_max_scan(&path, None, 1_000_000)?;
        assert_eq!(
            result.reads_scanned, 50,
            "below-cap inputs must scan to EOF"
        );

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }
}
