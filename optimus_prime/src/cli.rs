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
    /// For Phred64, TrimGalore adds +31 to the cutoff value itself
    /// (because Cutadapt always operates in Phred33 space internally).
    /// Since we handle the offset directly in our quality trimmer,
    /// we just return the raw cutoff — the phred_offset parameter
    /// handles the encoding difference.
    pub fn effective_quality_cutoff(&self) -> u8 {
        self.quality
    }
}
