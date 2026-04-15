//! Trimming report generation.
//!
//! Generates TrimGalore-compatible trimming reports that can be parsed by MultiQC.

use std::io::Write;

/// Statistics collected during trimming.
#[derive(Debug, Default)]
pub struct TrimStats {
    /// Total sequences processed
    pub total_reads: usize,
    /// Reads where at least one adapter matched in any round
    pub total_reads_with_adapter: usize,
    /// Per-adapter: number of trimming events (one read may count for multiple adapters across rounds)
    pub reads_with_adapter: Vec<usize>,
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
    /// Per-adapter length distributions: adapter_length_counts[adapter_idx][match_len] = count
    pub adapter_length_counts: Vec<Vec<usize>>,
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
    /// Create a TrimStats pre-sized for `n` adapters.
    pub fn with_adapter_count(n: usize) -> Self {
        TrimStats {
            reads_with_adapter: vec![0; n],
            adapter_length_counts: vec![Vec::new(); n],
            ..Default::default()
        }
    }

    /// Merge another batch's stats into this accumulator.
    pub fn merge(&mut self, other: &TrimStats) {
        self.total_reads += other.total_reads;
        self.total_reads_with_adapter += other.total_reads_with_adapter;
        self.bases_quality_trimmed += other.bases_quality_trimmed;
        self.total_bp_processed += other.total_bp_processed;
        self.total_bp_written += other.total_bp_written;
        self.total_bp_after_trim += other.total_bp_after_trim;
        self.too_short += other.too_short;
        self.too_long += other.too_long;
        self.too_many_n += other.too_many_n;
        self.reads_written += other.reads_written;
        // Merge per-adapter stats (element-wise addition across 2D vecs)
        // Ensure self has enough adapter slots
        if other.reads_with_adapter.len() > self.reads_with_adapter.len() {
            self.reads_with_adapter.resize(other.reads_with_adapter.len(), 0);
            self.adapter_length_counts.resize(other.adapter_length_counts.len(), Vec::new());
        }
        for (i, &count) in other.reads_with_adapter.iter().enumerate() {
            self.reads_with_adapter[i] += count;
        }
        for (i, other_lengths) in other.adapter_length_counts.iter().enumerate() {
            if other_lengths.len() > self.adapter_length_counts[i].len() {
                self.adapter_length_counts[i].resize(other_lengths.len(), 0);
            }
            for (j, &count) in other_lengths.iter().enumerate() {
                self.adapter_length_counts[i][j] += count;
            }
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
    /// R1 adapters: (name, sequence). Single-element for normal usage.
    pub adapters: Vec<(String, String)>,
    /// R2 adapters: (name, sequence). Empty = use R1 adapters.
    pub adapters_r2: Vec<(String, String)>,
    /// Max rounds of adapter trimming per read (-n/--times).
    pub times: usize,
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
    /// Input filename (basename only) — used in text report and Cutadapt-compatible section
    pub input_filename: String,
    /// All input filenames — SE: one element, PE: both R1 and R2. Used in JSON report.
    pub input_filenames: Vec<String>,
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
    if config.adapters.len() == 1 {
        writeln!(w, "Adapter sequence: '{}' (user-specified or auto-detected)", config.adapters[0].1)?;
    } else {
        write!(w, "Adapter sequence(s):")?;
        for (name, seq) in &config.adapters {
            write!(w, " '{}' [{}]", seq, name)?;
        }
        writeln!(w, " (user-specified)")?;
    }
    if !config.adapters_r2.is_empty() {
        if config.adapters_r2.len() == 1 {
            writeln!(w, "Optional adapter 2 sequence (Read 2): '{}'", config.adapters_r2[0].1)?;
        } else {
            write!(w, "Adapter 2 sequence(s) (Read 2):")?;
            for (name, seq) in &config.adapters_r2 {
                write!(w, " '{}' [{}]", seq, name)?;
            }
            writeln!(w)?;
        }
    }
    if config.times > 1 {
        writeln!(w, "Trimming rounds per read (-n): {}", config.times)?;
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
        format_number(stats.total_reads_with_adapter),
        percentage(stats.total_reads_with_adapter, stats.total_reads))?;
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
pub fn write_cutadapt_compatible_section<W: Write>(
    w: &mut W,
    config: &TrimConfig,
    stats: &TrimStats,
    read_number: u8,
) -> std::io::Result<()> {
    // Select adapter list for this read
    let adapters = if read_number == 2 && !config.adapters_r2.is_empty() {
        &config.adapters_r2
    } else {
        &config.adapters
    };

    writeln!(w)?;
    // Native identifier for MultiQC with Trim Galore v2.0 support
    writeln!(w, "Trim Galore {} (Oxidized Edition) — adapter trimming built in", config.version)?;
    // Backwards compatibility: older MultiQC discovers reports via "This is cutadapt"
    writeln!(w, "This is cutadapt 4.0 (compatible; for MultiQC backwards compatibility)")?;
    // Command line — MultiQC extracts the sample name from the input filename here
    let first_adapter_seq = adapters.first().map(|(_, s)| s.as_str()).unwrap_or("");
    writeln!(w, "Command line parameters: -j 1 -e {} -q {} -O {} -a {} {}",
        config.error_rate, config.quality_cutoff, config.stringency,
        first_adapter_seq, config.input_filename)?;
    writeln!(w, "Processing reads on 1 core in single-end mode ...")?;
    writeln!(w)?;

    // Summary — regexes: "Total reads processed:\s*([\d,]+)" etc.
    writeln!(w, "=== Summary ===")?;
    writeln!(w)?;
    writeln!(w, "Total reads processed:           {:>10}", format_number(stats.total_reads))?;
    writeln!(w, "Reads with adapters:             {:>10} ({:.1}%)",
        format_number(stats.total_reads_with_adapter),
        percentage(stats.total_reads_with_adapter, stats.total_reads))?;
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

    // Per-adapter sections — "=== Adapter N ===" triggers MultiQC section parsing
    for (idx, (_name, seq)) in adapters.iter().enumerate() {
        let times_trimmed = stats.reads_with_adapter.get(idx).copied().unwrap_or(0);
        let length_counts = stats.adapter_length_counts.get(idx);

        writeln!(w, "=== Adapter {} ===", idx + 1)?;
        writeln!(w)?;
        writeln!(w, "Sequence: {}; Type: regular 3'; Length: {}; Trimmed: {} times.",
            seq, seq.len(), times_trimmed)?;
        writeln!(w)?;

        write_allowed_errors(w, seq.len(), config.error_rate)?;
        writeln!(w)?;

        writeln!(w, "Overview of removed sequences")?;
        writeln!(w, "length\tcount\texpect\tmax.err\terror counts")?;

        if let Some(counts) = length_counts {
            for (length, &count) in counts.iter().enumerate() {
                if length == 0 || count == 0 {
                    continue;
                }
                let expect = stats.total_reads as f64 / 4f64.powi(length as i32);
                let max_err = (length as f64 * config.error_rate).floor() as usize;
                writeln!(w, "{}\t{}\t{:.1}\t{}\t{}", length, count, expect, max_err, count)?;
            }
        }
        writeln!(w)?;
        writeln!(w)?;
    }

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
    // Input filenames as array: SE = ["file.fq.gz"], PE = ["R1.fq.gz", "R2.fq.gz"]
    write!(w, "{}\"input_filenames\": [", i1)?;
    for (i, fname) in config.input_filenames.iter().enumerate() {
        write!(w, "\"{}\"", json_escape_string(fname))?;
        if i + 1 < config.input_filenames.len() {
            write!(w, ", ")?;
        }
    }
    writeln!(w, "],")?;
    json_str(w, i1, "mode", if config.paired { "paired-end" } else { "single-end" }, true)?;
    writeln!(w, "{}\"read_number\": {},", i1, read_number)?;
    json_str(w, i1, "command_line", &config.command_line, true)?;

    // ── Parameters ───────────────────────────────────────────────
    writeln!(w, "{}\"parameters\": {{", i1)?;
    writeln!(w, "{}\"quality_cutoff\": {},", i2, config.quality_cutoff)?;
    // Adapters as array of objects
    write!(w, "{}\"adapters\": [", i2)?;
    for (i, (name, seq)) in config.adapters.iter().enumerate() {
        write!(w, "{{\"name\": \"{}\", \"sequence\": \"{}\"}}", json_escape_string(name), json_escape_string(seq))?;
        if i + 1 < config.adapters.len() { write!(w, ", ")?; }
    }
    writeln!(w, "],")?;
    write!(w, "{}\"adapters_r2\": [", i2)?;
    for (i, (name, seq)) in config.adapters_r2.iter().enumerate() {
        write!(w, "{{\"name\": \"{}\", \"sequence\": \"{}\"}}", json_escape_string(name), json_escape_string(seq))?;
        if i + 1 < config.adapters_r2.len() { write!(w, ", ")?; }
    }
    writeln!(w, "],")?;
    json_int(w, i2, "times", config.times, true)?;
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
    json_int(w, i2, "reads_with_adapter", stats.total_reads_with_adapter, true)?;
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
    // For R2, use the R2-specific adapter list if set
    let adapters = if read_number == 2 && !config.adapters_r2.is_empty() {
        &config.adapters_r2
    } else {
        &config.adapters
    };
    writeln!(w, "{}\"adapter_trimming\": [", i1)?;
    for (adapter_idx, (name, seq)) in adapters.iter().enumerate() {
        let times_trimmed = stats.reads_with_adapter.get(adapter_idx).copied().unwrap_or(0);
        writeln!(w, "{}{{", i2)?;
        json_str(w, i3, "name", name, true)?;
        json_str(w, i3, "sequence", seq, true)?;
        json_str(w, i3, "type", "regular 3'", true)?;
        json_int(w, i3, "length", seq.len(), true)?;
        json_int(w, i3, "times_trimmed", times_trimmed, true)?;

        // Sparse length distribution for this adapter
        let length_counts = stats.adapter_length_counts.get(adapter_idx);
        let entries: Vec<(usize, usize)> = length_counts
            .map(|counts| {
                counts.iter().enumerate()
                    .filter(|&(i, &count)| i > 0 && count > 0)
                    .map(|(i, &count)| (i, count))
                    .collect()
            })
            .unwrap_or_default();

        if entries.is_empty() {
            writeln!(w, "{}\"length_distribution\": {{}}", i3)?;
        } else {
            writeln!(w, "{}\"length_distribution\": {{", i3)?;
            let i4 = "        ";
            for (entry_idx, &(length, count)) in entries.iter().enumerate() {
                let comma = if entry_idx < entries.len() - 1 { "," } else { "" };
                writeln!(w, "{}\"{}\": {}{}", i4, length, count, comma)?;
            }
            writeln!(w, "{}}}", i3)?;
        }
        let comma = if adapter_idx < adapters.len() - 1 { "," } else { "" };
        writeln!(w, "{}}}{}", i2, comma)?;
    }
    writeln!(w, "{}],", i1)?;

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
            adapters: vec![("adapter_1".to_string(), "AGATCGGAAGAGC".to_string())],
            adapters_r2: Vec::new(),
            times: 1,
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
            input_filenames: vec!["sample.fq.gz".to_string()],
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

    /// Helper to build TrimStats with some non-zero values (single adapter).
    fn test_stats() -> TrimStats {
        let mut stats = TrimStats::with_adapter_count(1);
        stats.total_reads = 10000;
        stats.total_reads_with_adapter = 5234;
        stats.reads_with_adapter = vec![5234];
        stats.reads_written = 9800;
        stats.too_short = 150;
        stats.too_long = 0;
        stats.too_many_n = 50;
        stats.total_bp_processed = 1_000_000;
        stats.bases_quality_trimmed = 50_000;
        stats.total_bp_written = 900_000;
        stats.adapter_length_counts = vec![vec![0, 1000, 500, 250, 0, 120]];
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
        assert_eq!(json["input_filenames"][0], "sample.fq.gz");
        assert_eq!(json["mode"], "single-end");
        assert_eq!(json["read_number"], 1);

        // Parameters
        assert_eq!(json["parameters"]["quality_cutoff"], 20);
        assert_eq!(json["parameters"]["adapters"][0]["sequence"], "AGATCGGAAGAGC");
        assert_eq!(json["parameters"]["adapters_r2"].as_array().unwrap().len(), 0);
        assert_eq!(json["parameters"]["times"], 1);
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

        // Adapter trimming — now an array
        let at = &json["adapter_trimming"];
        assert!(at.is_array());
        assert_eq!(at.as_array().unwrap().len(), 1);
        assert_eq!(at[0]["sequence"], "AGATCGGAAGAGC");
        assert_eq!(at[0]["type"], "regular 3'");
        assert_eq!(at[0]["length"], 13);
        assert_eq!(at[0]["times_trimmed"], 5234);

        // Length distribution is sparse
        let ld = &at[0]["length_distribution"];
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
        config.adapters_r2 = vec![("smallrna_r2".to_string(), "TGGAATTCTCGG".to_string())];
        config.input_filename = "sample_R1.fq.gz".to_string();
        config.input_filenames = vec!["sample_R1.fq.gz".to_string(), "sample_R2.fq.gz".to_string()];

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
        // Both filenames in the array
        assert_eq!(json["input_filenames"][0], "sample_R1.fq.gz");
        assert_eq!(json["input_filenames"][1], "sample_R2.fq.gz");
        assert_eq!(json["input_filenames"].as_array().unwrap().len(), 2);
        // R1 uses the primary adapter
        assert_eq!(json["adapter_trimming"][0]["sequence"], "AGATCGGAAGAGC");
        // pair_validation is populated for R1 too
        assert_eq!(json["pair_validation"]["pairs_analyzed"], 10000);
        assert_eq!(json["pair_validation"]["pairs_removed"], 200);
        assert_eq!(json["pair_validation"]["r1_unpaired"], 10);
        assert_eq!(json["parameters"]["clip_r2"], 2);

        // R2 report — uses R2-specific adapter
        config.input_filename = "sample_R2.fq.gz".to_string();
        // R2 stats need their own adapter count matching adapters_r2
        let mut stats_r2 = TrimStats::with_adapter_count(1);
        stats_r2.total_reads = 10000;
        stats_r2.total_reads_with_adapter = 4000;
        stats_r2.reads_with_adapter = vec![4000];
        stats_r2.reads_written = 9800;
        stats_r2.total_bp_processed = 1_000_000;
        stats_r2.total_bp_written = 900_000;

        let mut buf2 = Vec::new();
        write_json_report(&mut buf2, &config, &stats_r2, Some(&pair_stats), 2, &extra).unwrap();
        let json2: serde_json::Value = serde_json::from_slice(&buf2).unwrap();

        assert_eq!(json2["read_number"], 2);
        // R2 uses the R2-specific adapter
        assert_eq!(json2["adapter_trimming"][0]["sequence"], "TGGAATTCTCGG");
        assert_eq!(json2["adapter_trimming"][0]["length"], 12);
        // pair_validation also present in R2
        assert_eq!(json2["pair_validation"]["pairs_analyzed"], 10000);
    }

    #[test]
    fn test_write_json_report_sparse_length_distribution() {
        let config = test_config();
        let extra = test_extra_params();

        // Empty adapter_length_counts → empty object
        let mut stats = TrimStats::with_adapter_count(1);
        stats.total_reads = 100;
        stats.reads_written = 100;

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();
        let json: serde_json::Value = serde_json::from_slice(&buf).unwrap();

        let ld = &json["adapter_trimming"][0]["length_distribution"];
        assert!(ld.is_object());
        assert_eq!(ld.as_object().unwrap().len(), 0);

        // Only non-zero entries appear
        stats.adapter_length_counts = vec![vec![0, 0, 0, 42, 0, 0, 0, 7]];
        let mut buf2 = Vec::new();
        write_json_report(&mut buf2, &config, &stats, None, 1, &extra).unwrap();
        let json2: serde_json::Value = serde_json::from_slice(&buf2).unwrap();

        let ld2 = json2["adapter_trimming"][0]["length_distribution"].as_object().unwrap();
        assert_eq!(ld2.len(), 2);
        assert_eq!(ld2["3"], 42);
        assert_eq!(ld2["7"], 7);
    }

    #[test]
    fn test_write_json_report_special_characters() {
        let mut config = test_config();
        config.command_line = r#"trim_galore --fastqc_args "--nogroup" "my file.fq.gz""#.to_string();
        config.input_filename = "my file.fq.gz".to_string();
        config.input_filenames = vec!["my file.fq.gz".to_string()];

        let stats = TrimStats::with_adapter_count(1);
        let extra = test_extra_params();

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();

        // Must produce valid JSON despite quotes and spaces in values
        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("JSON with special characters must be valid");

        assert_eq!(json["input_filenames"][0], "my file.fq.gz");
        // Command line with embedded quotes is escaped correctly
        assert!(json["command_line"].as_str().unwrap().contains("--nogroup"));
    }

    #[test]
    fn test_write_json_report_multi_adapter() {
        // Test JSON report with 3 adapters, each with distinct per-adapter stats.
        let mut config = test_config();
        config.adapters = vec![
            ("illumina".to_string(), "AGATCGGAAGAGC".to_string()),
            ("nextera".to_string(), "CTGTCTCTTATA".to_string()),
            ("smallrna".to_string(), "TGGAATTCTCGG".to_string()),
        ];
        config.times = 2;

        let mut stats = TrimStats::with_adapter_count(3);
        stats.total_reads = 10000;
        stats.total_reads_with_adapter = 6000; // 6000 reads had >=1 adapter
        stats.reads_with_adapter = vec![4000, 1500, 800]; // per-adapter trimming events
        stats.reads_written = 9500;
        stats.total_bp_processed = 1_000_000;
        stats.total_bp_written = 850_000;
        stats.adapter_length_counts = vec![
            vec![0, 1000, 500, 250],   // illumina length distribution
            vec![0, 0, 300, 200, 100], // nextera length distribution
            vec![0, 400, 200],         // smallrna length distribution
        ];

        let extra = test_extra_params();
        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();
        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("Multi-adapter JSON must be valid");

        // Parameters: adapters array with 3 elements
        let params_adapters = json["parameters"]["adapters"].as_array().unwrap();
        assert_eq!(params_adapters.len(), 3);
        assert_eq!(params_adapters[0]["name"], "illumina");
        assert_eq!(params_adapters[0]["sequence"], "AGATCGGAAGAGC");
        assert_eq!(params_adapters[1]["name"], "nextera");
        assert_eq!(params_adapters[2]["name"], "smallrna");
        assert_eq!(json["parameters"]["times"], 2);

        // read_processing uses the aggregate total
        assert_eq!(json["read_processing"]["reads_with_adapter"], 6000);

        // adapter_trimming: array with 3 elements, each with per-adapter stats
        let at = json["adapter_trimming"].as_array().unwrap();
        assert_eq!(at.len(), 3);

        // Adapter 0: illumina
        assert_eq!(at[0]["name"], "illumina");
        assert_eq!(at[0]["sequence"], "AGATCGGAAGAGC");
        assert_eq!(at[0]["length"], 13);
        assert_eq!(at[0]["times_trimmed"], 4000);
        assert_eq!(at[0]["length_distribution"]["1"], 1000);
        assert_eq!(at[0]["length_distribution"]["2"], 500);
        assert_eq!(at[0]["length_distribution"]["3"], 250);

        // Adapter 1: nextera
        assert_eq!(at[1]["name"], "nextera");
        assert_eq!(at[1]["sequence"], "CTGTCTCTTATA");
        assert_eq!(at[1]["length"], 12);
        assert_eq!(at[1]["times_trimmed"], 1500);
        assert_eq!(at[1]["length_distribution"]["2"], 300);
        assert_eq!(at[1]["length_distribution"]["3"], 200);
        assert_eq!(at[1]["length_distribution"]["4"], 100);

        // Adapter 2: smallrna
        assert_eq!(at[2]["name"], "smallrna");
        assert_eq!(at[2]["sequence"], "TGGAATTCTCGG");
        assert_eq!(at[2]["length"], 12);
        assert_eq!(at[2]["times_trimmed"], 800);
        assert_eq!(at[2]["length_distribution"]["1"], 400);
        assert_eq!(at[2]["length_distribution"]["2"], 200);
    }

    #[test]
    fn test_stats_merge_multi_adapter() {
        // Simulate two parallel workers producing per-adapter stats for 3 adapters,
        // then merging them into the main accumulator.
        let mut main_stats = TrimStats::with_adapter_count(3);

        // Worker A: saw 100 reads, adapter 0 matched 50 times, adapter 2 matched 10 times
        let mut worker_a = TrimStats::with_adapter_count(3);
        worker_a.total_reads = 100;
        worker_a.total_reads_with_adapter = 55;
        worker_a.reads_with_adapter = vec![50, 0, 10];
        worker_a.adapter_length_counts = vec![
            vec![0, 20, 30],     // adapter 0: 20 events at len=1, 30 at len=2
            vec![],              // adapter 1: no matches
            vec![0, 0, 0, 5, 5], // adapter 2: 5 events at len=3, 5 at len=4
        ];

        // Worker B: saw 80 reads, adapter 0 matched 30 times, adapter 1 matched 15 times
        let mut worker_b = TrimStats::with_adapter_count(3);
        worker_b.total_reads = 80;
        worker_b.total_reads_with_adapter = 40;
        worker_b.reads_with_adapter = vec![30, 15, 0];
        worker_b.adapter_length_counts = vec![
            vec![0, 10, 15, 5],  // adapter 0: longer length distribution than worker A
            vec![0, 0, 8, 7],    // adapter 1: 8 at len=2, 7 at len=3
            vec![],              // adapter 2: no matches
        ];

        main_stats.merge(&worker_a);
        main_stats.merge(&worker_b);

        // Totals
        assert_eq!(main_stats.total_reads, 180);
        assert_eq!(main_stats.total_reads_with_adapter, 95);

        // Per-adapter event counts
        assert_eq!(main_stats.reads_with_adapter, vec![80, 15, 10]);

        // Adapter 0 length distribution: worker_a has [0,20,30], worker_b has [0,10,15,5]
        // Merged: [0, 30, 45, 5]
        assert_eq!(main_stats.adapter_length_counts[0], vec![0, 30, 45, 5]);

        // Adapter 1: worker_a has [], worker_b has [0,0,8,7] → [0,0,8,7]
        assert_eq!(main_stats.adapter_length_counts[1], vec![0, 0, 8, 7]);

        // Adapter 2: worker_a has [0,0,0,5,5], worker_b has [] → [0,0,0,5,5]
        assert_eq!(main_stats.adapter_length_counts[2], vec![0, 0, 0, 5, 5]);
    }

    #[test]
    fn test_write_json_report_all_zero_stats() {
        let config = test_config();
        let stats = TrimStats::with_adapter_count(1); // all zeros
        let extra = test_extra_params();

        let mut buf = Vec::new();
        write_json_report(&mut buf, &config, &stats, None, 1, &extra).unwrap();

        let json: serde_json::Value = serde_json::from_slice(&buf)
            .expect("All-zero JSON must be valid");

        assert_eq!(json["read_processing"]["total_reads"], 0);
        assert_eq!(json["read_processing"]["reads_written"], 0);
        assert_eq!(json["basepair_processing"]["total_bp_processed"], 0);
        assert_eq!(json["adapter_trimming"][0]["times_trimmed"], 0);
    }
}
