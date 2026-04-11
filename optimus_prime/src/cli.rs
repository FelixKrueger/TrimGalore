//! Command-line argument parsing and validation.

use clap::Parser;
use std::path::PathBuf;

/// Optimus Prime — A fast, single-pass NGS adapter and quality trimmer.
///
/// Rust reimplementation of Trim Galore with single-pass paired-end processing.
/// Output is compatible with TrimGalore for use with MultiQC and existing pipelines.
#[derive(Parser, Debug)]
#[clap(name = "optimus_prime", version, about)]
pub struct Cli {
    /// Input FASTQ file(s). For paired-end, provide two files.
    #[clap(required = true)]
    pub input: Vec<PathBuf>,

    /// Quality trimming cutoff (Phred score). Bases below this are trimmed from 3' end.
    #[clap(short = 'q', long = "quality", default_value = "20")]
    pub quality: u8,

    /// Adapter sequence for trimming. Auto-detected if not specified.
    #[clap(short = 'a', long = "adapter")]
    pub adapter: Option<String>,

    /// Optional adapter sequence for Read 2 (paired-end only).
    /// Auto-set by --small_rna and --bgiseq presets.
    #[clap(long = "adapter2", alias = "a2")]
    pub adapter2: Option<String>,

    /// Use Illumina universal adapter (AGATCGGAAGAGC). This is the default.
    #[clap(long = "illumina", conflicts_with_all = &["nextera", "small_rna", "stranded_illumina", "bgiseq"])]
    pub illumina: bool,

    /// Use Nextera transposase adapter (CTGTCTCTTATA).
    #[clap(long = "nextera", conflicts_with_all = &["illumina", "small_rna", "stranded_illumina", "bgiseq"])]
    pub nextera: bool,

    /// Use Illumina Small RNA adapter (TGGAATTCTCGG). Sets --adapter2 for Read 2.
    #[clap(long = "small_rna", conflicts_with_all = &["illumina", "nextera", "stranded_illumina", "bgiseq"])]
    pub small_rna: bool,

    /// Use Illumina Stranded mRNA adapter (ACTGTCTCTTATA).
    #[clap(long = "stranded_illumina", conflicts_with_all = &["illumina", "nextera", "small_rna", "bgiseq"])]
    pub stranded_illumina: bool,

    /// Use BGI/DNBSEQ adapter. Sets --adapter2 for Read 2.
    #[clap(long = "bgiseq", conflicts_with_all = &["illumina", "nextera", "small_rna", "stranded_illumina"])]
    pub bgiseq: bool,

    /// Paired-end mode. Expects exactly 2 input files.
    #[clap(long = "paired")]
    pub paired: bool,

    /// Maximum allowed error rate for adapter matching (0-1).
    #[clap(short = 'e', long = "error_rate", default_value = "0.1")]
    pub error_rate: f64,

    /// Minimum overlap with adapter sequence required to trim (stringency).
    #[clap(long = "stringency", default_value = "1")]
    pub stringency: usize,

    /// Minimum required sequence length after trimming.
    /// Default: 20 (18 for smallRNA adapter).
    #[clap(long = "length")]
    pub length: Option<usize>,

    /// Maximum allowed sequence length (discard reads longer than this).
    #[clap(long = "max_length")]
    pub max_length: Option<usize>,

    /// Maximum number of N bases allowed in a read.
    /// Integer: absolute count. Decimal (0-1): fraction of read length.
    #[clap(long = "max_n")]
    pub max_n: Option<f64>,

    /// Trim N bases from both ends of reads.
    #[clap(long = "trim-n", alias = "trim_n")]
    pub trim_n: bool,

    /// Remove N bases from the 5' end of Read 1.
    #[clap(long = "clip_R1")]
    pub clip_r1: Option<usize>,

    /// Remove N bases from the 5' end of Read 2 (paired-end only).
    #[clap(long = "clip_R2")]
    pub clip_r2: Option<usize>,

    /// Remove N bases from the 3' end of Read 1.
    #[clap(long = "three_prime_clip_R1")]
    pub three_prime_clip_r1: Option<usize>,

    /// Remove N bases from the 3' end of Read 2 (paired-end only).
    #[clap(long = "three_prime_clip_R2")]
    pub three_prime_clip_r2: Option<usize>,

    /// NextSeq/NovaSeq 2-colour quality trimming. Trailing high-quality G bases
    /// are treated as no-signal artifacts and quality-trimmed. The value is the
    /// quality cutoff (replaces -q). Mutually exclusive with --quality.
    #[clap(long = "nextseq", alias = "2colour")]
    pub nextseq: Option<u8>,

    /// Use Phred+64 quality encoding (Illumina 1.5). Default is Phred+33.
    #[clap(long = "phred64")]
    pub phred64: bool,

    /// Use Phred+33 quality encoding (default, Sanger/Illumina 1.8+).
    #[clap(long = "phred33")]
    pub phred33: bool,

    /// Output directory for trimmed files.
    #[clap(short = 'o', long = "output_dir")]
    pub output_dir: Option<PathBuf>,

    /// Custom basename for output files (replaces input filename stem).
    #[clap(long = "basename")]
    pub basename: Option<String>,

    /// Do not gzip-compress output files.
    #[clap(long = "dont_gzip")]
    pub dont_gzip: bool,

    /// Suppress the trimming report.
    #[clap(long = "no_report_file")]
    pub no_report_file: bool,

    /// Retain unpaired reads when the mate is too short (paired-end only).
    #[clap(long = "retain_unpaired")]
    pub retain_unpaired: bool,

    /// Minimum length for unpaired Read 1 (with --retain_unpaired).
    #[clap(short = 'r', long = "length_1", default_value = "35", alias = "r1")]
    pub length_1: usize,

    /// Minimum length for unpaired Read 2 (with --retain_unpaired).
    #[clap(long = "length_2", default_value = "35", alias = "r2")]
    pub length_2: usize,

    /// Add clipped sequences to read IDs (format: :clip5:SEQ:clip3:SEQ).
    #[clap(long = "rename")]
    pub rename: bool,

    /// If auto-detected adapter count is at or below this threshold,
    /// skip adapter trimming (only quality trimming proceeds).
    /// Incompatible with explicit adapter presets.
    #[clap(long = "consider_already_trimmed",
           conflicts_with_all = &["illumina", "nextera", "small_rna", "stranded_illumina", "bgiseq"])]
    pub consider_already_trimmed: Option<usize>,

    /// RRBS mode for MspI-digested samples. Removes 2bp end-repair artifacts
    /// at MspI cut sites after adapter trimming. In paired-end directional mode,
    /// automatically sets --clip_R2 2 unless the user provides their own value.
    #[clap(long = "rrbs")]
    pub rrbs: bool,

    /// Non-directional RRBS libraries. Reads starting with CAA or CGA get 2bp
    /// trimmed from the 5' end. Requires --rrbs.
    #[clap(long = "non_directional", requires = "rrbs")]
    pub non_directional: bool,

    /// Run FastQC on the trimmed output files.
    #[clap(long = "fastqc")]
    pub fastqc: bool,

    /// Additional arguments to pass to FastQC. Implies --fastqc.
    #[clap(long = "fastqc_args")]
    pub fastqc_args: Option<String>,

    /// Number of compression threads for gzip output (default: 1).
    /// Values > 1 enable parallel gzip compression for faster I/O.
    #[clap(short = 'j', long = "cores", default_value = "1")]
    pub cores: usize,

    /// Trim poly-A tails from the 3' end of Read 1 (and single-end reads),
    /// and poly-T heads from the 5' end of Read 2. Runs after adapter trimming,
    /// so poly-A tails hidden behind adapters are also removed.
    /// Matches Cutadapt's --poly-a behavior.
    #[clap(long = "poly_a", alias = "poly-a", alias = "polyA")]
    pub poly_a: bool,
}

impl Cli {
    /// Validate CLI arguments after parsing.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.paired && self.input.len() != 2 {
            anyhow::bail!(
                "Paired-end mode requires exactly 2 input files, got {}",
                self.input.len()
            );
        }

        if !self.paired && self.input.len() > 1 {
            anyhow::bail!(
                "Single-end mode expects 1 input file, got {}. Use --paired for paired-end.",
                self.input.len()
            );
        }

        if self.error_rate < 0.0 || self.error_rate > 1.0 {
            anyhow::bail!("Error rate must be between 0 and 1, got {}", self.error_rate);
        }

        if self.stringency == 0 {
            anyhow::bail!("Stringency (minimum overlap) must be at least 1");
        }

        if self.nextseq.is_some() && self.quality != 20 {
            anyhow::bail!(
                "--nextseq/--2colour and -q/--quality are mutually exclusive. \
                 The nextseq value replaces the quality cutoff."
            );
        }

        if let Some(val) = self.nextseq {
            if val == 0 || val >= 200 {
                anyhow::bail!(
                    "NextSeq quality cutoff must be between 1 and 199, got {}",
                    val
                );
            }
        }

        if let Some(threshold) = self.consider_already_trimmed {
            if threshold > 10000 {
                anyhow::bail!(
                    "consider_already_trimmed value must be between 0 and 10000, got {}",
                    threshold
                );
            }
        }

        if self.cores == 0 {
            anyhow::bail!("--cores must be at least 1");
        }

        // Check input files exist
        for path in &self.input {
            if !path.exists() {
                anyhow::bail!("Input file not found: {}", path.display());
            }
        }

        Ok(())
    }

    /// Get the Phred encoding offset.
    pub fn phred_offset(&self) -> u8 {
        if self.phred64 { 64 } else { 33 }
    }

    /// Get the effective quality cutoff value.
    ///
    /// If `--nextseq` is set, uses that value as the cutoff.
    /// Otherwise uses the standard `--quality` value.
    pub fn effective_quality_cutoff(&self) -> u8 {
        self.nextseq.unwrap_or(self.quality)
    }
}
