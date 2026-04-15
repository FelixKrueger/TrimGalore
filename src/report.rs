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
    /// Total basepairs in input reads (before any trimming)
    pub total_bp_processed: usize,
    /// Total basepairs written to output (after all trimming + filtering)
    pub total_bp_written: usize,
    /// Total basepairs after adapter/quality trimming but before length/pair filtering.
    /// This matches what Cutadapt would have reported as "Total written (filtered)".
    pub total_bp_after_trim: usize,
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
    /// RRBS directional PE: R2 reads that had 5' clip applied (auto-set clip_r2=2)
    pub rrbs_r2_clipped_5prime: usize,
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
        self.total_bp_processed += other.total_bp_processed;
        self.total_bp_written += other.total_bp_written;
        self.total_bp_after_trim += other.total_bp_after_trim;
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
        self.rrbs_r2_clipped_5prime += other.rrbs_r2_clipped_5prime;
        self.poly_a_trimmed += other.poly_a_trimmed;
        self.poly_a_bases_trimmed += other.poly_a_bases_trimmed;
        self.poly_g_trimmed += other.poly_g_trimmed;
        self.poly_g_bases_trimmed += other.poly_g_bases_trimmed;
        self.discarded_untrimmed += other.discarded_untrimmed;
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
    /// Input filename (basename only) — used in Cutadapt-compatible section for MultiQC
    pub input_filename: String,
}

/// Write the trimming report header (parameter summary).
pub fn write_report_header<W: Write>(w: &mut W, config: &TrimConfig) -> std::io::Result<()> {
    writeln!(w, "SUMMARISING RUN PARAMETERS")?;
    writeln!(w, "=========================")?;
    writeln!(w, "Input filename: {}", config.input_filename)?;
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

/// Write a Cutadapt-compatible section for MultiQC parsing.
///
/// MultiQC has no dedicated Trim Galore module — it uses its Cutadapt parser,
/// which requires "This is cutadapt" in the first 100 lines, a "Command line
/// parameters:" line for sample name extraction, bp statistics, and an adapter
/// length distribution table.
pub fn write_cutadapt_section<W: Write>(
    w: &mut W,
    config: &TrimConfig,
    stats: &TrimStats,
) -> std::io::Result<()> {
    writeln!(w)?;
    // Native identifier for MultiQC with Trim Galore v2.0 support
    writeln!(w, "Trim Galore {} (Oxidized Edition) — adapter trimming built in", config.version)?;
    // Backwards compatibility: older MultiQC discovers reports via "This is cutadapt"
    writeln!(w, "This is cutadapt 4.0 (compatible; for MultiQC backwards compatibility)")?;
    // Command line — MultiQC extracts the sample name from the input filename here
    writeln!(w, "Command line parameters: -j 1 -e {} -q {} -O {} -a {} {}",
        config.error_rate, config.quality_cutoff, config.stringency,
        config.adapter, config.input_filename)?;
    writeln!(w, "Processing reads on 1 core in single-end mode ...")?;
    writeln!(w)?;

    // Summary — regexes: "Total reads processed:\s*([\d,]+)" etc.
    writeln!(w, "=== Summary ===")?;
    writeln!(w)?;
    writeln!(w, "Total reads processed:           {:>10}", format_number(stats.total_reads))?;
    writeln!(w, "Reads with adapters:             {:>10} ({:.1}%)",
        format_number(stats.reads_with_adapter),
        percentage(stats.reads_with_adapter, stats.total_reads))?;
    if stats.discarded_untrimmed > 0 {
        writeln!(w, "Reads discarded as untrimmed:    {:>10} ({:.1}%)",
            format_number(stats.discarded_untrimmed),
            percentage(stats.discarded_untrimmed, stats.total_reads))?;
    }
    // In v0.6.x, Cutadapt wrote ALL reads — length/pair filtering was a separate
    // TrimGalore step. So Cutadapt's "Reads written" was total_reads minus only
    // --discard-untrimmed. We replicate that here for MultiQC backwards compatibility.
    let cutadapt_reads_written = stats.total_reads - stats.discarded_untrimmed;
    writeln!(w, "Reads written (passing filters): {:>10} ({:.1}%)",
        format_number(cutadapt_reads_written),
        percentage(cutadapt_reads_written, stats.total_reads))?;
    writeln!(w)?;

    // Basepair statistics — "Total basepairs processed:\s*([\d,]+) bp"
    writeln!(w, "Total basepairs processed:   {:>12} bp", format_number(stats.total_bp_processed))?;
    writeln!(w, "Quality-trimmed:             {:>12} bp ({:.1}%)",
        format_number(stats.bases_quality_trimmed),
        percentage(stats.bases_quality_trimmed, stats.total_bp_processed))?;
    // Same as reads: report bp after trimming but before length/pair filtering,
    // matching what Cutadapt would have written.
    writeln!(w, "Total written (filtered):    {:>12} bp ({:.1}%)",
        format_number(stats.total_bp_after_trim),
        percentage(stats.total_bp_after_trim, stats.total_bp_processed))?;
    writeln!(w)?;

    // Adapter section — "=== Adapter 1 ===" triggers section parsing
    writeln!(w, "=== Adapter 1 ===")?;
    writeln!(w)?;
    // "Sequence: ...; Type: regular 3'; Length: N; Trimmed: N times."
    writeln!(w, "Sequence: {}; Type: regular 3'; Length: {}; Trimmed: {} times.",
        config.adapter, config.adapter.len(), stats.reads_with_adapter)?;
    writeln!(w)?;

    // Allowed errors (cosmetic — MultiQC ignores this, but users may read it)
    write_allowed_errors(w, config.adapter.len(), config.error_rate)?;
    writeln!(w)?;

    // Adapter length distribution table
    // MultiQC parses: header with "length", "count", "expect"; then rows ^(\d+)\s+(\d+)\s+([\d\.]+)
    writeln!(w, "Overview of removed sequences")?;
    writeln!(w, "length\tcount\texpect\tmax.err\terror counts")?;

    for (length, &count) in stats.adapter_length_counts.iter().enumerate() {
        if length == 0 || count == 0 {
            continue;
        }
        // Expected: total_reads / 4^length (random match probability)
        let expect = stats.total_reads as f64 / 4f64.powi(length as i32);
        let max_err = (length as f64 * config.error_rate).floor() as usize;
        // Simplified error counts: assume all matches have 0 errors
        // (MultiQC only uses the first 3 columns)
        writeln!(w, "{}\t{}\t{:.1}\t{}\t{}", length, count, expect, max_err, count)?;
    }
    writeln!(w)?;
    writeln!(w)?;

    Ok(())
}

/// Write the "RUN STATISTICS" footer (TG-specific, not parsed by MultiQC).
pub fn write_run_footer<W: Write>(
    w: &mut W,
    config: &TrimConfig,
    stats: &TrimStats,
) -> std::io::Result<()> {
    writeln!(w, "RUN STATISTICS FOR INPUT FILE: {}", config.input_filename)?;
    writeln!(w, "=============================================")?;
    writeln!(w, "{} sequences processed in total", stats.total_reads)?;

    if stats.too_short > 0 {
        writeln!(w, "Sequences removed because they became shorter than the length cutoff of {} bp:\t{} ({:.1}%)",
            config.length_cutoff, stats.too_short, percentage(stats.too_short, stats.total_reads))?;
    }
    if stats.too_long > 0 {
        writeln!(w, "Sequences removed because they were longer than the upper length cutoff:\t{} ({:.1}%)",
            stats.too_long, percentage(stats.too_long, stats.total_reads))?;
    }
    if stats.too_many_n > 0 {
        writeln!(w, "Sequences removed because of too many N bases:\t{} ({:.1}%)",
            stats.too_many_n, percentage(stats.too_many_n, stats.total_reads))?;
    }

    if stats.rrbs_trimmed_3prime > 0 {
        writeln!(w, "RRBS reads trimmed by additional 2 bp when adapter contamination was detected:\t{} ({:.1}%)",
            stats.rrbs_trimmed_3prime, percentage(stats.rrbs_trimmed_3prime, stats.total_reads))?;
    }
    if stats.rrbs_trimmed_5prime > 0 {
        writeln!(w, "RRBS reads trimmed by 2 bp at the start when read started with CAA or CGA in total:\t{} ({:.1}%)",
            stats.rrbs_trimmed_5prime, percentage(stats.rrbs_trimmed_5prime, stats.total_reads))?;
    }
    if stats.poly_a_trimmed > 0 {
        writeln!(w, "Reads with poly-A/T trimmed:\t{} ({:.1}%); {} bp removed",
            stats.poly_a_trimmed, percentage(stats.poly_a_trimmed, stats.total_reads),
            format_number(stats.poly_a_bases_trimmed))?;
    }
    if stats.poly_g_trimmed > 0 {
        writeln!(w, "Reads with poly-G/C trimmed:\t{} ({:.1}%); {} bp removed",
            stats.poly_g_trimmed, percentage(stats.poly_g_trimmed, stats.total_reads),
            format_number(stats.poly_g_bases_trimmed))?;
    }

    writeln!(w)?;

    Ok(())
}

// ── JSON report ──────────────────────────────────────────────────────────────

/// Extra parameters for the JSON report not present on [`TrimConfig`].
///
/// These come from `trimmer::TrimConfig` and the CLI, and are passed separately
/// to avoid coupling report.rs to trimmer.rs.
#[derive(Debug)]
pub struct JsonReportParams {
    pub clip_r1: Option<usize>,
    pub clip_r2: Option<usize>,
    pub three_prime_clip_r1: Option<usize>,
    pub three_prime_clip_r2: Option<usize>,
    pub max_n: Option<f64>,
    pub discard_untrimmed: bool,
    pub consider_already_trimmed: Option<usize>,
}

/// Write a structured JSON trimming report.
///
/// Produces a machine-readable report for MultiQC's native TrimGalore module.
/// The schema is versioned (`schema_version: 1`) so MultiQC can handle future
/// changes without breaking.
pub fn write_json_report<W: Write>(
    w: &mut W,
    config: &TrimConfig,
    stats: &TrimStats,
    pair_stats: Option<&PairValidationStats>,
    read_number: u8,
    extra: &JsonReportParams,
) -> std::io::Result<()> {
    let i1 = "  ";
    let i2 = "    ";
    let i3 = "      ";

    writeln!(w, "{{")?;

    // ── Top-level fields ─────────────────────────────────────────
    json_str(w, i1, "tool", "Trim Galore", true)?;
    writeln!(w, "{}\"schema_version\": 1,", i1)?;
    json_str(w, i1, "trim_galore_version", &config.version, true)?;
    json_str(w, i1, "input_filename", &config.input_filename, true)?;
    json_str(w, i1, "mode", if config.paired { "paired-end" } else { "single-end" }, true)?;
    writeln!(w, "{}\"read_number\": {},", i1, read_number)?;
    json_str(w, i1, "command_line", &config.command_line, true)?;

    // ── Parameters ───────────────────────────────────────────────
    writeln!(w, "{}\"parameters\": {{", i1)?;
    writeln!(w, "{}\"quality_cutoff\": {},", i2, config.quality_cutoff)?;
    json_str(w, i2, "adapter", &config.adapter, true)?;
    json_opt_str(w, i2, "adapter_r2", config.adapter_r2.as_deref(), true)?;
    json_float(w, i2, "error_rate", config.error_rate, true)?;
    json_int(w, i2, "stringency", config.stringency, true)?;
    json_int(w, i2, "length_cutoff", config.length_cutoff, true)?;
    json_opt_int(w, i2, "max_length", config.max_length, true)?;
    writeln!(w, "{}\"phred_encoding\": {},", i2, config.phred_encoding)?;
    json_bool(w, i2, "trim_n", config.trim_n, true)?;
    json_bool(w, i2, "nextseq", config.nextseq, true)?;
    json_bool(w, i2, "rrbs", config.rrbs, true)?;
    json_bool(w, i2, "non_directional", config.non_directional, true)?;
    json_bool(w, i2, "poly_a", config.poly_a, true)?;
    json_bool(w, i2, "poly_g", config.poly_g, true)?;
    json_opt_int(w, i2, "clip_r1", extra.clip_r1, true)?;
    json_opt_int(w, i2, "clip_r2", extra.clip_r2, true)?;
    json_opt_int(w, i2, "three_prime_clip_r1", extra.three_prime_clip_r1, true)?;
    json_opt_int(w, i2, "three_prime_clip_r2", extra.three_prime_clip_r2, true)?;
    json_opt_float(w, i2, "max_n", extra.max_n, true)?;
    json_bool(w, i2, "discard_untrimmed", extra.discard_untrimmed, true)?;
    json_opt_int(w, i2, "consider_already_trimmed", extra.consider_already_trimmed, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── Read processing ──────────────────────────────────────────
    writeln!(w, "{}\"read_processing\": {{", i1)?;
    json_int(w, i2, "total_reads", stats.total_reads, true)?;
    json_int(w, i2, "reads_with_adapter", stats.reads_with_adapter, true)?;
    json_int(w, i2, "reads_written", stats.reads_written, true)?;
    json_int(w, i2, "reads_too_short", stats.too_short, true)?;
    json_int(w, i2, "reads_too_long", stats.too_long, true)?;
    json_int(w, i2, "reads_too_many_n", stats.too_many_n, true)?;
    json_int(w, i2, "reads_discarded_untrimmed", stats.discarded_untrimmed, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── Basepair processing ──────────────────────────────────────
    writeln!(w, "{}\"basepair_processing\": {{", i1)?;
    json_int(w, i2, "total_bp_processed", stats.total_bp_processed, true)?;
    json_int(w, i2, "quality_trimmed_bp", stats.bases_quality_trimmed, true)?;
    json_int(w, i2, "total_bp_written", stats.total_bp_written, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── Adapter trimming ─────────────────────────────────────────
    // For R2, use the R2-specific adapter if set (Small RNA, BGI presets)
    let adapter_seq = if read_number == 2 {
        config.adapter_r2.as_deref().unwrap_or(&config.adapter)
    } else {
        &config.adapter
    };
    writeln!(w, "{}\"adapter_trimming\": {{", i1)?;
    json_str(w, i2, "sequence", adapter_seq, true)?;
    json_str(w, i2, "type", "regular 3'", true)?;
    json_int(w, i2, "length", adapter_seq.len(), true)?;
    json_int(w, i2, "times_trimmed", stats.reads_with_adapter, true)?;

    // Sparse length distribution: omit index 0 and zero-count entries
    let entries: Vec<(usize, usize)> = stats
        .adapter_length_counts
        .iter()
        .enumerate()
        .filter(|&(i, &count)| i > 0 && count > 0)
        .map(|(i, &count)| (i, count))
        .collect();

    if entries.is_empty() {
        writeln!(w, "{}\"length_distribution\": {{}}", i2)?;
    } else {
        writeln!(w, "{}\"length_distribution\": {{", i2)?;
        for (idx, &(length, count)) in entries.iter().enumerate() {
            let comma = if idx < entries.len() - 1 { "," } else { "" };
            writeln!(w, "{}\"{}\": {}{}", i3, length, count, comma)?;
        }
        writeln!(w, "{}}}", i2)?;
    }
    writeln!(w, "{}}},", i1)?;

    // ── Poly-A trimming ──────────────────────────────────────────
    writeln!(w, "{}\"poly_a_trimming\": {{", i1)?;
    json_int(w, i2, "reads_trimmed", stats.poly_a_trimmed, true)?;
    json_int(w, i2, "bases_removed", stats.poly_a_bases_trimmed, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── Poly-G trimming ──────────────────────────────────────────
    writeln!(w, "{}\"poly_g_trimming\": {{", i1)?;
    json_int(w, i2, "reads_trimmed", stats.poly_g_trimmed, true)?;
    json_int(w, i2, "bases_removed", stats.poly_g_bases_trimmed, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── RRBS ─────────────────────────────────────────────────────
    writeln!(w, "{}\"rrbs\": {{", i1)?;
    json_int(w, i2, "trimmed_3prime", stats.rrbs_trimmed_3prime, true)?;
    json_int(w, i2, "trimmed_5prime", stats.rrbs_trimmed_5prime, true)?;
    json_int(w, i2, "r2_clipped_5prime", stats.rrbs_r2_clipped_5prime, false)?;
    writeln!(w, "{}}},", i1)?;

    // ── Pair validation ──────────────────────────────────────────
    // In PE mode: included in BOTH R1 and R2 reports (self-contained per file).
    // In SE mode: null.
    match pair_stats {
        Some(ps) => {
            writeln!(w, "{}\"pair_validation\": {{", i1)?;
            json_int(w, i2, "pairs_analyzed", ps.pairs_analyzed, true)?;
            json_int(w, i2, "pairs_removed", ps.pairs_removed, true)?;
            json_int(w, i2, "pairs_removed_n", ps.pairs_removed_n, true)?;
            json_int(w, i2, "pairs_removed_too_long", ps.pairs_removed_too_long, true)?;
            json_int(w, i2, "r1_unpaired", ps.r1_unpaired, true)?;
            json_int(w, i2, "r2_unpaired", ps.r2_unpaired, false)?;
            writeln!(w, "{}}}", i1)?;
        }
        None => {
            json_null(w, i1, "pair_validation", false)?;
        }
    }

    writeln!(w, "}}")?;
    Ok(())
}

// ── JSON writing helpers ─────────────────────────────────────────────────────

/// Escape a string for safe inclusion in JSON output.
fn json_escape_string(s: &str) -> String {
    let mut escaped = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                // Control characters → \u00xx
                escaped.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => escaped.push(c),
        }
    }
    escaped
}

fn json_str<W: Write>(w: &mut W, indent: &str, key: &str, value: &str, comma: bool) -> std::io::Result<()> {
    writeln!(w, "{}\"{}\": \"{}\"{}",
        indent, key, json_escape_string(value), if comma { "," } else { "" })
}

fn json_int<W: Write>(w: &mut W, indent: &str, key: &str, value: usize, comma: bool) -> std::io::Result<()> {
    writeln!(w, "{}\"{}\": {}{}",
        indent, key, value, if comma { "," } else { "" })
}

fn json_float<W: Write>(w: &mut W, indent: &str, key: &str, value: f64, comma: bool) -> std::io::Result<()> {
    if !value.is_finite() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("non-finite float value for JSON key \"{}\": {}", key, value),
        ));
    }
    writeln!(w, "{}\"{}\": {}{}",
        indent, key, value, if comma { "," } else { "" })
}

fn json_bool<W: Write>(w: &mut W, indent: &str, key: &str, value: bool, comma: bool) -> std::io::Result<()> {
    writeln!(w, "{}\"{}\": {}{}",
        indent, key, if value { "true" } else { "false" }, if comma { "," } else { "" })
}

fn json_null<W: Write>(w: &mut W, indent: &str, key: &str, comma: bool) -> std::io::Result<()> {
    writeln!(w, "{}\"{}\": null{}",
        indent, key, if comma { "," } else { "" })
}

fn json_opt_int<W: Write>(w: &mut W, indent: &str, key: &str, value: Option<usize>, comma: bool) -> std::io::Result<()> {
    match value {
        Some(v) => json_int(w, indent, key, v, comma),
        None => json_null(w, indent, key, comma),
    }
}

fn json_opt_str<W: Write>(w: &mut W, indent: &str, key: &str, value: Option<&str>, comma: bool) -> std::io::Result<()> {
    match value {
        Some(v) => json_str(w, indent, key, v, comma),
        None => json_null(w, indent, key, comma),
    }
}

fn json_opt_float<W: Write>(w: &mut W, indent: &str, key: &str, value: Option<f64>, comma: bool) -> std::io::Result<()> {
    match value {
        Some(v) => json_float(w, indent, key, v, comma),
        None => json_null(w, indent, key, comma),
    }
}

/// Write the "No. of allowed errors" section (cosmetic, not parsed by MultiQC).
fn write_allowed_errors<W: Write>(w: &mut W, adapter_len: usize, error_rate: f64) -> std::io::Result<()> {
    writeln!(w, "No. of allowed errors:")?;
    let mut parts = Vec::new();
    let mut range_start = 1;
    let mut current_max_err = 0usize;

    for len in 1..=adapter_len {
        let max_err = (len as f64 * error_rate).floor() as usize;
        if max_err != current_max_err {
            // Close previous range
            if range_start <= len - 1 {
                parts.push(format!("{}-{} bp: {}", range_start, len - 1, current_max_err));
            }
            range_start = len;
            current_max_err = max_err;
        }
    }
    // Close final range
    parts.push(format!("{}-{} bp: {}", range_start, adapter_len, current_max_err));
    writeln!(w, "{}", parts.join("; "))?;

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

    // ── JSON report tests ────────────────────────────────────────

    #[test]
    fn test_json_escape_string() {
        // Normal strings pass through unchanged
        assert_eq!(json_escape_string("hello"), "hello");
        assert_eq!(json_escape_string("AGATCGGAAGAGC"), "AGATCGGAAGAGC");

        // Backslash and quotes
        assert_eq!(json_escape_string(r#"say "hi""#), r#"say \"hi\""#);
        assert_eq!(json_escape_string(r"path\to\file"), r"path\\to\\file");

        // Control characters
        assert_eq!(json_escape_string("line1\nline2"), "line1\\nline2");
        assert_eq!(json_escape_string("col1\tcol2"), "col1\\tcol2");
        assert_eq!(json_escape_string("cr\r"), "cr\\r");

        // Other control chars get \u00xx encoding
        assert_eq!(json_escape_string("\x01"), "\\u0001");
        assert_eq!(json_escape_string("\x1f"), "\\u001f");

        // Empty string
        assert_eq!(json_escape_string(""), "");
    }

    /// Helper to build a minimal TrimConfig for testing.
    fn test_config() -> TrimConfig {
        TrimConfig {
            version: "2.0.0".to_string(),
            quality_cutoff: 20,
            adapter: "AGATCGGAAGAGC".to_string(),
            adapter_r2: None,
            error_rate: 0.1,
            stringency: 1,
            length_cutoff: 20,
            max_length: None,
            paired: false,
            gzip: true,
            trim_n: false,
            nextseq: false,
            rrbs: false,
            non_directional: false,
            phred_encoding: 33,
            poly_a: false,
            poly_g: false,
            command_line: "trim_galore sample.fq.gz".to_string(),
            input_filename: "sample.fq.gz".to_string(),
        }
    }

    /// Helper to build minimal JsonReportParams for testing.
    fn test_extra_params() -> JsonReportParams {
        JsonReportParams {
            clip_r1: None,
            clip_r2: None,
            three_prime_clip_r1: None,
            three_prime_clip_r2: None,
            max_n: None,
            discard_untrimmed: false,
            consider_already_trimmed: None,
        }
    }

    /// Helper to build TrimStats with some non-zero values.
    fn test_stats() -> TrimStats {
        let mut stats = TrimStats::default();
        stats.total_reads = 10000;
        stats.reads_with_adapter = 5234;
        stats.reads_written = 9800;
        stats.too_short = 150;
        stats.too_long = 0;
        stats.too_many_n = 50;
        stats.total_bp_processed = 1_000_000;
        stats.bases_quality_trimmed = 50_000;
        stats.total_bp_written = 900_000;
        stats.adapter_length_counts = vec![0, 1000, 500, 250, 0, 120];
        stats
    }

    #[test]
    fn test_write_json_report_se() {
        let config = test_config();
        let stats = test_stats();
        let extra = test_extra_params();

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();

        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("JSON output must be valid");

        assert_eq!(json["tool"], "Trim Galore");
        assert_eq!(json["schema_version"], 1);
        assert_eq!(json["trim_galore_version"], "2.0.0");
        assert_eq!(json["input_filename"], "sample.fq.gz");
        assert_eq!(json["mode"], "single-end");
        assert_eq!(json["read_number"], 1);

        // Parameters
        assert_eq!(json["parameters"]["quality_cutoff"], 20);
        assert_eq!(json["parameters"]["adapter"], "AGATCGGAAGAGC");
        assert!(json["parameters"]["adapter_r2"].is_null());
        assert_eq!(json["parameters"]["error_rate"], 0.1);
        assert_eq!(json["parameters"]["stringency"], 1);
        assert_eq!(json["parameters"]["length_cutoff"], 20);
        assert!(json["parameters"]["max_length"].is_null());
        assert_eq!(json["parameters"]["trim_n"], false);
        assert_eq!(json["parameters"]["discard_untrimmed"], false);
        assert!(json["parameters"]["consider_already_trimmed"].is_null());

        // Read processing
        assert_eq!(json["read_processing"]["total_reads"], 10000);
        assert_eq!(json["read_processing"]["reads_with_adapter"], 5234);
        assert_eq!(json["read_processing"]["reads_written"], 9800);
        assert_eq!(json["read_processing"]["reads_too_short"], 150);

        // Basepair processing
        assert_eq!(json["basepair_processing"]["total_bp_processed"], 1_000_000);
        assert_eq!(json["basepair_processing"]["quality_trimmed_bp"], 50_000);
        assert_eq!(json["basepair_processing"]["total_bp_written"], 900_000);

        // Adapter trimming
        assert_eq!(json["adapter_trimming"]["sequence"], "AGATCGGAAGAGC");
        assert_eq!(json["adapter_trimming"]["type"], "regular 3'");
        assert_eq!(json["adapter_trimming"]["length"], 13);
        assert_eq!(json["adapter_trimming"]["times_trimmed"], 5234);

        // Length distribution is sparse
        let ld = &json["adapter_trimming"]["length_distribution"];
        assert_eq!(ld["1"], 1000);
        assert_eq!(ld["2"], 500);
        assert_eq!(ld["3"], 250);
        assert_eq!(ld["5"], 120);
        assert!(ld.get("0").is_none()); // index 0 omitted
        assert!(ld.get("4").is_none()); // zero-count omitted

        // Sections with zeros are still present
        assert_eq!(json["poly_a_trimming"]["reads_trimmed"], 0);
        assert_eq!(json["poly_g_trimming"]["reads_trimmed"], 0);
        assert_eq!(json["rrbs"]["trimmed_3prime"], 0);

        // SE: pair_validation is null
        assert!(json["pair_validation"].is_null());
    }

    #[test]
    fn test_write_json_report_pe() {
        let mut config = test_config();
        config.paired = true;
        config.adapter_r2 = Some("TGGAATTCTCGG".to_string());
        config.input_filename = "sample_R1.fq.gz".to_string();

        let stats = test_stats();
        let extra = JsonReportParams {
            clip_r1: None,
            clip_r2: Some(2),
            ..test_extra_params()
        };

        let pair_stats = PairValidationStats {
            pairs_analyzed: 10000,
            pairs_removed: 200,
            pairs_removed_n: 50,
            pairs_removed_too_long: 0,
            r1_unpaired: 10,
            r2_unpaired: 15,
        };

        // R1 report
        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, Some(&pair_stats), 1, &extra).unwrap();
        let json: serde_json::Value = serde_json::from_slice(&buf).unwrap();

        assert_eq!(json["mode"], "paired-end");
        assert_eq!(json["read_number"], 1);
        // R1 uses the primary adapter
        assert_eq!(json["adapter_trimming"]["sequence"], "AGATCGGAAGAGC");
        // pair_validation is populated for R1 too
        assert_eq!(json["pair_validation"]["pairs_analyzed"], 10000);
        assert_eq!(json["pair_validation"]["pairs_removed"], 200);
        assert_eq!(json["pair_validation"]["r1_unpaired"], 10);
        assert_eq!(json["parameters"]["clip_r2"], 2);

        // R2 report
        config.input_filename = "sample_R2.fq.gz".to_string();
        let mut buf2 = Vec::new();
        write_json_report(&mut buf2, &config, &stats, Some(&pair_stats), 2, &extra).unwrap();
        let json2: serde_json::Value = serde_json::from_slice(&buf2).unwrap();

        assert_eq!(json2["read_number"], 2);
        // R2 uses the R2-specific adapter
        assert_eq!(json2["adapter_trimming"]["sequence"], "TGGAATTCTCGG");
        assert_eq!(json2["adapter_trimming"]["length"], 12);
        // pair_validation also present in R2
        assert_eq!(json2["pair_validation"]["pairs_analyzed"], 10000);
    }

    #[test]
    fn test_write_json_report_sparse_length_distribution() {
        let config = test_config();
        let extra = test_extra_params();

        // Empty adapter_length_counts → empty object
        let mut stats = TrimStats::default();
        stats.total_reads = 100;
        stats.reads_written = 100;

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();
        let json: serde_json::Value = serde_json::from_slice(&buf).unwrap();

        let ld = &json["adapter_trimming"]["length_distribution"];
        assert!(ld.is_object());
        assert_eq!(ld.as_object().unwrap().len(), 0);

        // Only non-zero entries appear
        stats.adapter_length_counts = vec![0, 0, 0, 42, 0, 0, 0, 7];
        let mut buf2 = Vec::new();
        write_json_report(&mut buf2, &config, &stats, None, 1, &extra).unwrap();
        let json2: serde_json::Value = serde_json::from_slice(&buf2).unwrap();

        let ld2 = json2["adapter_trimming"]["length_distribution"].as_object().unwrap();
        assert_eq!(ld2.len(), 2);
        assert_eq!(ld2["3"], 42);
        assert_eq!(ld2["7"], 7);
    }

    #[test]
    fn test_write_json_report_special_characters() {
        let mut config = test_config();
        config.command_line = r#"trim_galore --fastqc_args "--nogroup" "my file.fq.gz""#.to_string();
        config.input_filename = "my file.fq.gz".to_string();

        let stats = TrimStats::default();
        let extra = test_extra_params();

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();

        // Must produce valid JSON despite quotes and spaces in values
        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("JSON with special characters must be valid");

        assert_eq!(json["input_filename"], "my file.fq.gz");
        // Command line with embedded quotes is escaped correctly
        assert!(json["command_line"].as_str().unwrap().contains("--nogroup"));
    }

    #[test]
    fn test_write_json_report_all_zero_stats() {
        let config = test_config();
        let stats = TrimStats::default(); // all zeros
        let extra = test_extra_params();

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();

        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("All-zero JSON must be valid");

        assert_eq!(json["read_processing"]["total_reads"], 0);
        assert_eq!(json["read_processing"]["reads_written"], 0);
        assert_eq!(json["basepair_processing"]["total_bp_processed"], 0);
        assert_eq!(json["adapter_trimming"]["times_trimmed"], 0);
    }
}
