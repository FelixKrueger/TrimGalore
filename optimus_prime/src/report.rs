//! Trimming report generation.
//!
//! Generates TrimGalore-compatible trimming reports that can be parsed by MultiQC.

use std::io::Write;

/// Statistics collected during trimming.
#[derive(Debug, Default)]
pub struct TrimStats {
    /// Total sequences processed
    pub total_reads: usize,
    /// Reads with adapter detected
    pub reads_with_adapter: usize,
    /// Total bases quality-trimmed
    pub bases_quality_trimmed: usize,
    /// Reads removed for being too short
    pub too_short: usize,
    /// Reads removed for being too long
    pub too_long: usize,
    /// Reads removed for too many Ns
    pub too_many_n: usize,
    /// Reads written to output
    pub reads_written: usize,
    /// Per-length adapter match counts (for the length distribution table)
    pub adapter_length_counts: Vec<usize>,
}

/// Statistics for paired-end validation.
#[derive(Debug, Default)]
pub struct PairValidationStats {
    /// Total read pairs analyzed
    pub pairs_analyzed: usize,
    /// Pairs removed (either too short, too long, or too many Ns)
    pub pairs_removed: usize,
    /// Pairs removed specifically for N-content
    pub pairs_removed_n: usize,
    /// Pairs removed for being too long
    pub pairs_removed_too_long: usize,
    /// R1 reads written to unpaired output
    pub r1_unpaired: usize,
    /// R2 reads written to unpaired output
    pub r2_unpaired: usize,
}

/// Configuration that was used for this trimming run (for the report header).
#[derive(Debug)]
pub struct TrimConfig {
    pub version: String,
    pub quality_cutoff: u8,
    pub adapter: String,
    pub adapter_r2: Option<String>,
    pub error_rate: f64,
    pub stringency: usize,
    pub length_cutoff: usize,
    pub max_length: Option<usize>,
    pub paired: bool,
    pub gzip: bool,
    pub trim_n: bool,
    pub phred_encoding: u8,
    pub command_line: String,
}

/// Write the trimming report header (parameter summary).
pub fn write_report_header<W: Write>(w: &mut W, config: &TrimConfig) -> std::io::Result<()> {
    writeln!(w, "SUMMARISING RUN PARAMETERS")?;
    writeln!(w, "=========================")?;
    writeln!(w, "Input filename: (from command line)")?;
    writeln!(w, "Trimming mode: {}", if config.paired { "paired-end" } else { "single-end" })?;
    writeln!(w, "Trim Galore version: {} (Optimus Prime)", config.version)?;
    writeln!(w, "Quality Phred score cutoff: {}", config.quality_cutoff)?;
    writeln!(w, "Quality encoding type selected: ASCII+{}", config.phred_encoding)?;
    writeln!(w, "Adapter sequence: '{}' (user-specified or auto-detected)", config.adapter)?;
    if let Some(ref a2) = config.adapter_r2 {
        writeln!(w, "Optional adapter 2 sequence: '{}'", a2)?;
    }
    writeln!(w, "Maximum trimming error rate: {}", config.error_rate)?;
    writeln!(w, "Minimum required adapter overlap (stringency): {} bp", config.stringency)?;
    writeln!(w, "Minimum required sequence length {}-end: {} bp",
        if config.paired { "for both reads before a sequence pair gets removed" } else { "single" },
        config.length_cutoff)?;
    if config.trim_n {
        writeln!(w, "Removing Ns from the end of reads")?;
    }
    if config.gzip {
        writeln!(w, "Output file will be GZIP compressed")?;
    }
    writeln!(w)?;

    Ok(())
}

/// Write the trimming run statistics.
pub fn write_run_stats<W: Write>(w: &mut W, stats: &TrimStats) -> std::io::Result<()> {
    writeln!(w, "=== Summary ===")?;
    writeln!(w)?;
    writeln!(w, "Total reads processed:              {:>10}", format_number(stats.total_reads))?;
    writeln!(w, "Reads with adapters:                {:>10} ({:.1}%)",
        format_number(stats.reads_with_adapter),
        percentage(stats.reads_with_adapter, stats.total_reads))?;
    writeln!(w)?;

    if stats.total_reads > 0 {
        writeln!(w, "Reads that were too short:          {:>10} ({:.1}%)",
            format_number(stats.too_short),
            percentage(stats.too_short, stats.total_reads))?;
        writeln!(w, "Reads that were too long:           {:>10} ({:.1}%)",
            format_number(stats.too_long),
            percentage(stats.too_long, stats.total_reads))?;
        writeln!(w, "Reads with too many N:              {:>10} ({:.1}%)",
            format_number(stats.too_many_n),
            percentage(stats.too_many_n, stats.total_reads))?;
    }

    writeln!(w, "Reads written (passing filters):    {:>10} ({:.1}%)",
        format_number(stats.reads_written),
        percentage(stats.reads_written, stats.total_reads))?;
    writeln!(w)?;

    Ok(())
}

/// Write paired-end validation statistics.
pub fn write_pair_validation_stats<W: Write>(
    w: &mut W,
    stats: &PairValidationStats,
) -> std::io::Result<()> {
    writeln!(w, "=== Paired-end validation ===")?;
    writeln!(w)?;
    writeln!(w, "Number of sequence pairs analysed:  {:>10}", format_number(stats.pairs_analyzed))?;
    writeln!(w)?;
    writeln!(w, "Number of sequence pairs removed because at least one read was shorter than the length cutoff ({:.2}%):",
        percentage(stats.pairs_removed, stats.pairs_analyzed))?;
    writeln!(w, "                             {:>10}", format_number(stats.pairs_removed))?;

    if stats.pairs_removed_n > 0 {
        writeln!(w, "Number of pairs removed for N content ({:.2}%):",
            percentage(stats.pairs_removed_n, stats.pairs_analyzed))?;
        writeln!(w, "                             {:>10}", format_number(stats.pairs_removed_n))?;
    }

    if stats.r1_unpaired > 0 || stats.r2_unpaired > 0 {
        writeln!(w)?;
        writeln!(w, "Unpaired read 1 kept:        {:>10}", format_number(stats.r1_unpaired))?;
        writeln!(w, "Unpaired read 2 kept:        {:>10}", format_number(stats.r2_unpaired))?;
    }

    writeln!(w)?;

    Ok(())
}

fn percentage(part: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 / total as f64 * 100.0
    }
}

fn format_number(n: usize) -> String {
    // Add comma separators for readability
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(ch);
    }
    result.chars().rev().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(0), "0");
        assert_eq!(format_number(100), "100");
        assert_eq!(format_number(1000), "1,000");
        assert_eq!(format_number(1000000), "1,000,000");
        assert_eq!(format_number(12345678), "12,345,678");
    }

    #[test]
    fn test_percentage() {
        assert_eq!(percentage(50, 100), 50.0);
        assert_eq!(percentage(0, 100), 0.0);
        assert_eq!(percentage(0, 0), 0.0);
    }
}
