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
    /// RRBS: reads trimmed 2bp from 3' end (adapter contamination at MspI site)
    pub rrbs_trimmed_3prime: usize,
    /// RRBS non-directional: reads trimmed 2bp from 5' end (CAA/CGA)
    pub rrbs_trimmed_5prime: usize,
    /// Reads with poly-A/poly-T tails trimmed
    pub poly_a_trimmed: usize,
    /// Total bases removed by poly-A/poly-T trimming
    pub poly_a_bases_trimmed: usize,
    /// Reads with poly-G/poly-C tails trimmed
    pub poly_g_trimmed: usize,
    /// Total bases removed by poly-G/poly-C trimming
    pub poly_g_bases_trimmed: usize,
    /// Reads discarded because no adapter was found (--discard-untrimmed)
    pub discarded_untrimmed: usize,
}

impl TrimStats {
    /// Merge another batch's stats into this accumulator.
    pub fn merge(&mut self, other: &TrimStats) {
        self.total_reads += other.total_reads;
        self.reads_with_adapter += other.reads_with_adapter;
        self.bases_quality_trimmed += other.bases_quality_trimmed;
        self.too_short += other.too_short;
        self.too_long += other.too_long;
        self.too_many_n += other.too_many_n;
        self.reads_written += other.reads_written;
        // Merge adapter length distribution (element-wise addition)
        if other.adapter_length_counts.len() > self.adapter_length_counts.len() {
            self.adapter_length_counts.resize(other.adapter_length_counts.len(), 0);
        }
        for (i, &count) in other.adapter_length_counts.iter().enumerate() {
            self.adapter_length_counts[i] += count;
        }
        self.rrbs_trimmed_3prime += other.rrbs_trimmed_3prime;
        self.rrbs_trimmed_5prime += other.rrbs_trimmed_5prime;
        self.poly_a_trimmed += other.poly_a_trimmed;
        self.poly_a_bases_trimmed += other.poly_a_bases_trimmed;
        self.poly_g_trimmed += other.poly_g_trimmed;
        self.poly_g_bases_trimmed += other.poly_g_bases_trimmed;
    }
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

impl PairValidationStats {
    /// Merge another batch's pair stats into this accumulator.
    pub fn merge(&mut self, other: &PairValidationStats) {
        self.pairs_analyzed += other.pairs_analyzed;
        self.pairs_removed += other.pairs_removed;
        self.pairs_removed_n += other.pairs_removed_n;
        self.pairs_removed_too_long += other.pairs_removed_too_long;
        self.r1_unpaired += other.r1_unpaired;
        self.r2_unpaired += other.r2_unpaired;
    }
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
    pub nextseq: bool,
    pub rrbs: bool,
    pub non_directional: bool,
    pub phred_encoding: u8,
    pub poly_a: bool,
    pub poly_g: bool,
    pub command_line: String,
}

/// Write the trimming report header (parameter summary).
pub fn write_report_header<W: Write>(w: &mut W, config: &TrimConfig) -> std::io::Result<()> {
    writeln!(w, "SUMMARISING RUN PARAMETERS")?;
    writeln!(w, "=========================")?;
    writeln!(w, "Input filename: (from command line)")?;
    writeln!(w, "Trimming mode: {}", if config.paired { "paired-end" } else { "single-end" })?;
    writeln!(w, "Trim Galore version: {} (Oxidized Edition)", config.version)?;
    if config.nextseq {
        writeln!(w, "2-colour high quality G-trimming enabled, with quality cutoff: --nextseq-trim={}", config.quality_cutoff)?;
    } else {
        writeln!(w, "Quality Phred score cutoff: {}", config.quality_cutoff)?;
    }
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
    if config.rrbs {
        writeln!(w, "File was specified to be an MspI-digested RRBS sample. Read 1 sequences with adapter contamination will be trimmed a further 2 bp from their 3' end, and Read 2 sequences will be trimmed by 2 bp from their 5' end to remove potential methylation-biased bases from the end-repair reaction")?;
    }
    if config.non_directional {
        writeln!(w, "File was specified to be a non-directional MspI-digested RRBS sample. Sequences starting with either 'CAA' or 'CGA' will have the first 2 bp trimmed off to remove potential methylation-biased bases from the end-repair reaction")?;
    }
    if config.poly_a {
        writeln!(w, "Poly-A trimming enabled: removing poly-A tails from 3' end of R1/SE reads, and poly-T heads from 5' end of R2 reads")?;
    }
    if config.poly_g {
        writeln!(w, "Poly-G trimming enabled: removing poly-G tails from 3' end of R1/SE reads, and poly-C heads from 5' end of R2 reads")?;
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
        if stats.discarded_untrimmed > 0 {
            writeln!(w, "Reads discarded as untrimmed:       {:>10} ({:.1}%)",
                format_number(stats.discarded_untrimmed),
                percentage(stats.discarded_untrimmed, stats.total_reads))?;
        }
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

    if stats.rrbs_trimmed_3prime > 0 {
        writeln!(w, "RRBS reads trimmed by additional 2 bp when adapter contamination was detected:\t{} ({:.1}%)",
            stats.rrbs_trimmed_3prime,
            percentage(stats.rrbs_trimmed_3prime, stats.total_reads))?;
    }
    if stats.rrbs_trimmed_5prime > 0 {
        writeln!(w, "RRBS reads trimmed by 2 bp at the start when read started with CAA or CGA in total:\t{} ({:.1}%)",
            stats.rrbs_trimmed_5prime,
            percentage(stats.rrbs_trimmed_5prime, stats.total_reads))?;
    }
    if stats.poly_a_trimmed > 0 {
        writeln!(w, "Reads with poly-A/T trimmed:    {:>10} ({:.1}%)",
            format_number(stats.poly_a_trimmed),
            percentage(stats.poly_a_trimmed, stats.total_reads))?;
        writeln!(w, "  Poly-A/T bases removed:       {:>10}",
            format_number(stats.poly_a_bases_trimmed))?;
    }
    if stats.poly_g_trimmed > 0 {
        writeln!(w, "Reads with poly-G/C trimmed:    {:>10} ({:.1}%)",
            format_number(stats.poly_g_trimmed),
            percentage(stats.poly_g_trimmed, stats.total_reads))?;
        writeln!(w, "  Poly-G/C bases removed:       {:>10}",
            format_number(stats.poly_g_bases_trimmed))?;
    }

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
