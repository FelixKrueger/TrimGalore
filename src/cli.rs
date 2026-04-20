//! Command-line argument parsing and validation.

use clap::Parser;
use std::path::PathBuf;

/// Trim Galore - Oxidized Edition: A fast, single-pass NGS adapter and quality trimmer.
///
/// Drop-in replacement for Trim Galore, rewritten in Rust. Produces byte-identical
/// output, compatible with MultiQC and existing pipelines.
#[derive(Parser, Debug)]
#[clap(
    name = "trim_galore",
    version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)"),
    long_version = concat!(
        env!("CARGO_PKG_VERSION"), " (Oxidized Edition)\n",
        env!("VERSION_BODY")
    ),
    about
)]
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
    #[clap(long = "fastqc_args", allow_hyphen_values = true)]
    pub fastqc_args: Option<String>,

    /// Number of worker threads for parallel processing (default: 1).
    /// Values > 1 run trimming and gzip compression across multiple threads.
    /// Near-linear speedup up to ~16 cores; diminishing returns beyond ~20
    /// (typically I/O-bound at that point, not algorithmic).
    #[clap(short = 'j', long = "cores", default_value = "1")]
    pub cores: usize,

    /// Trim poly-A tails from the 3' end of Read 1 (and single-end reads),
    /// and poly-T heads from the 5' end of Read 2. Runs after adapter trimming,
    /// so poly-A tails hidden behind adapters are also removed.
    #[clap(long = "poly_a", alias = "poly-a", alias = "polyA")]
    pub poly_a: bool,

    /// Trim poly-G tails from the 3' end of Read 1 (and single-end reads),
    /// and poly-C heads from the 5' end of Read 2. Useful for data from
    /// 2-colour instruments (NovaSeq, NextSeq) where no-signal bases are
    /// called as high-quality G. By default, poly-G trimming is auto-detected
    /// from the data. Use this flag to force-enable it.
    /// This is independent from --nextseq (quality-based G-trimming).
    #[clap(
        long = "poly_g",
        alias = "poly-g",
        alias = "polyG",
        conflicts_with = "no_poly_g"
    )]
    pub poly_g: bool,

    /// Disable poly-G auto-detection and trimming.
    #[clap(long = "no_poly_g", alias = "no-poly-g", alias = "no-polyG")]
    pub no_poly_g: bool,

    /// Number of adapter trimming rounds per read. With multiple adapters,
    /// this allows removing more than one adapter from the same read.
    /// Default: 1. Typical multi-adapter usage: -n 3.
    #[clap(short = 'n', long = "times", default_value = "1")]
    pub times: usize,

    /// Discard reads that did not contain an adapter sequence. Only reads
    /// where at least one adapter match was found are written to output.
    /// For paired-end, the pair is discarded if neither read had an adapter.
    #[clap(long = "discard_untrimmed", alias = "discard-untrimmed")]
    pub discard_untrimmed: bool,

    // --- Specialty modes (run-and-exit, bypass normal trimming) ---
    /// Hard-trim to keep only the first N bases from the 5' end.
    /// Bypasses adapter/quality trimming entirely.
    #[clap(long = "hardtrim5")]
    pub hardtrim5: Option<usize>,

    /// Hard-trim to keep only the last N bases from the 3' end.
    /// Bypasses adapter/quality trimming entirely.
    #[clap(long = "hardtrim3")]
    pub hardtrim3: Option<usize>,

    /// Epigenetic Clock mode (paired-end only). Extracts 8bp UMI + 4bp
    /// fixed sequence (CAGT) from both reads, appends to read IDs, and
    /// clips R1 at position 13, R2 at position 15. Bypasses normal trimming.
    #[clap(long = "clock", alias = "casio", alias = "breitling")]
    pub clock: bool,

    /// Transfer the first N bases from Read 2 as a UMI barcode to both
    /// read IDs, then clip R2 by N bases. Paired-end only.
    /// Bypasses normal trimming (IMPLICON preprocessing).
    /// Default UMI length: 8 (used when --implicon is given without a value).
    #[clap(long = "implicon", alias = "umi_from_r2",
           default_missing_value = "8", num_args = 0..=1, require_equals = true)]
    pub implicon: Option<usize>,

    /// Demultiplex reads after trimming based on 3' inline barcodes.
    /// Takes a barcode file (TSV: sample_name\tbarcode_sequence).
    /// Barcode is removed from the read and appended to the read ID.
    /// Single-end only.
    #[clap(long = "demux")]
    pub demux: Option<PathBuf>,

    // --- Deprecated flags (accepted for backwards compatibility, no-ops) ---
    /// [Deprecated] Output is gzipped by default in v2.0. Use --dont_gzip to disable.
    #[clap(long = "gzip", hide = true)]
    pub gzip: bool,

    /// [Deprecated] No longer needed — Cutadapt is built in.
    #[clap(long = "path_to_cutadapt", hide = true)]
    pub path_to_cutadapt: Option<String>,

    /// [Deprecated] No longer needed — Cutadapt is built in.
    #[clap(long = "cutadapt_args", hide = true, allow_hyphen_values = true)]
    pub cutadapt_args: Option<String>,

    /// [Deprecated] No longer needed — no Cutadapt subprocess.
    #[clap(long = "suppress_warn", hide = true)]
    pub suppress_warn: bool,

    /// [Deprecated] Reports are generated by default.
    #[clap(long = "report", hide = true)]
    pub report: bool,

    /// [Deprecated] Not yet supported in v2.0.
    #[clap(long = "keep", hide = true)]
    pub keep: bool,

    /// Easter egg (no-op).
    #[clap(long = "hulu", hide = true)]
    pub hulu: bool,
}

impl Cli {
    /// Validate CLI arguments after parsing.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.paired {
            // `#[clap(required = true)]` on `input` guarantees at least one file
            // reaches validate(), so no is_empty() check is needed.
            if self.input.len() % 2 != 0 {
                anyhow::bail!(
                    "Paired-end mode requires an even number of input files (R1/R2 pairs), got {}",
                    self.input.len()
                );
            }
            // Catch common user error: same file passed twice. Matches Perl's
            // byte-equality check at trim_galore:3208 (`$ARGV[$i] eq $ARGV[$i+1]`).
            // Does not follow symlinks or canonicalise.
            for chunk in self.input.chunks(2) {
                if chunk[0] == chunk[1] {
                    anyhow::bail!(
                        "Read 1 and Read 2 appear to be the same file: {}. \
                         Did you mean to pass distinct R1 and R2 files?",
                        chunk[0].display()
                    );
                }
            }
            // Catch accidental copy-paste of the same pair twice. Without this,
            // a duplicate pair trips the case-insensitive output-collision
            // pre-flight in main() with an APFS/NTFS message that misdirects
            // the user — here we emit a precise "duplicate of pair N" error.
            let pairs: Vec<(&std::path::PathBuf, &std::path::PathBuf)> =
                self.input.chunks(2).map(|c| (&c[0], &c[1])).collect();
            for (i, (r1, r2)) in pairs.iter().enumerate() {
                for (j, (pr1, pr2)) in pairs.iter().enumerate().take(i) {
                    if r1 == pr1 && r2 == pr2 {
                        anyhow::bail!(
                            "Pair {} ({}, {}) is a duplicate of pair {}. \
                             Did you mean to pass different files?",
                            i + 1,
                            r1.display(),
                            r2.display(),
                            j + 1
                        );
                    }
                }
            }
        }

        if !self.paired && self.input.len() > 1 && self.basename.is_some() {
            anyhow::bail!(
                "--basename cannot be used with multiple input files (ambiguous output naming)"
            );
        }
        if self.paired && self.input.len() > 2 && self.basename.is_some() {
            anyhow::bail!(
                "--basename cannot be used with multiple paired-end pairs (ambiguous output naming)"
            );
        }

        if self.error_rate < 0.0 || self.error_rate > 1.0 {
            anyhow::bail!(
                "Error rate must be between 0 and 1, got {}",
                self.error_rate
            );
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

        if self.times == 0 || self.times > 10 {
            anyhow::bail!("--times/-n must be between 1 and 10, got {}", self.times);
        }

        if self.cores == 0 {
            anyhow::bail!("--cores must be at least 1");
        }

        if let Some(n) = self.hardtrim5 {
            if n == 0 || n >= 1000 {
                anyhow::bail!("--hardtrim5 must be between 1 and 999, got {}", n);
            }
        }
        if let Some(n) = self.hardtrim3 {
            if n == 0 || n >= 1000 {
                anyhow::bail!("--hardtrim3 must be between 1 and 999, got {}", n);
            }
        }
        if self.clock && self.input.len() != 2 {
            anyhow::bail!("--clock requires exactly 2 paired-end input files");
        }
        if self.implicon.is_some() && self.input.len() != 2 {
            anyhow::bail!("--implicon requires exactly 2 paired-end input files");
        }
        if let Some(ref demux_file) = self.demux {
            if self.paired {
                anyhow::bail!("Demultiplexing is only allowed for single-end files");
            }
            if !demux_file.exists() {
                anyhow::bail!("Barcode file not found: {}", demux_file.display());
            }
        }

        // Check input files exist
        for path in &self.input {
            if !path.exists() {
                anyhow::bail!("Input file not found: {}", path.display());
            }
        }

        // Deprecation warnings for Perl-era flags
        if self.gzip {
            eprintln!(
                "WARNING: --gzip is deprecated in Trim Galore v2.0. Output is gzipped by default. Use --dont_gzip to disable. Ignoring."
            );
        }
        if self.path_to_cutadapt.is_some() {
            eprintln!(
                "WARNING: --path_to_cutadapt is deprecated in Trim Galore v2.0 (no external Cutadapt needed). Ignoring."
            );
        }
        if self.cutadapt_args.is_some() {
            eprintln!(
                "WARNING: --cutadapt_args is deprecated in Trim Galore v2.0 (no external Cutadapt needed). Ignoring."
            );
            eprintln!("         Note: --discard-untrimmed is now a native flag.");
        }
        if self.suppress_warn {
            eprintln!(
                "WARNING: --suppress_warn is deprecated in Trim Galore v2.0 (no Cutadapt subprocess). Ignoring."
            );
        }
        if self.keep {
            eprintln!(
                "WARNING: --keep is not yet supported in Trim Galore v2.0. RRBS reads below length cutoff will be removed. Ignoring."
            );
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    // All fixtures live in test_files/ and are guaranteed to exist under the
    // repo root during `cargo test` (cwd = crate root).
    const R1: &str = "test_files/BS-seq_10K_R1.fastq.gz";
    const R2: &str = "test_files/BS-seq_10K_R2.fastq.gz";
    const ALT_R1: &str = "test_files/SRR24766921_RRBS_R1.fastq.gz";
    const ALT_R2: &str = "test_files/SRR24766921_RRBS_R2.fastq.gz";

    #[test]
    fn test_validate_paired_odd_count_rejected() {
        for inputs in [
            vec![R1],
            vec![R1, R2, ALT_R1],
            vec![R1, R2, ALT_R1, ALT_R2, R1],
        ] {
            let mut argv = vec!["trim_galore", "--paired"];
            argv.extend(inputs.iter().copied());
            let cli = Cli::parse_from(argv);
            let err = cli.validate().unwrap_err().to_string();
            assert!(
                err.contains("even number of input files"),
                "expected even-number error, got: {err}"
            );
        }
    }

    #[test]
    fn test_validate_paired_r1_r2_equal_rejected_within_pair() {
        let cli = Cli::parse_from(["trim_galore", "--paired", R1, R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("appear to be the same file"),
            "expected within-pair duplicate error, got: {err}"
        );
    }

    #[test]
    fn test_validate_paired_duplicate_pair_rejected_across_pairs() {
        let cli = Cli::parse_from(["trim_galore", "--paired", R1, R2, R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("duplicate of pair"),
            "expected cross-pair duplicate error, got: {err}"
        );
    }

    #[test]
    fn test_validate_paired_basename_rejected_multi_pair() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--basename",
            "foo",
            R1,
            R2,
            ALT_R1,
            ALT_R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("basename cannot be used with multiple"),
            "expected multi-pair basename rejection, got: {err}"
        );
    }

    #[test]
    fn test_validate_paired_single_end_basename_still_allowed() {
        // Regression guard: SE `--basename` with a single input must still pass.
        let cli = Cli::parse_from(["trim_galore", "--basename", "foo", R1]);
        cli.validate()
            .expect("SE --basename with one input should validate");
    }

    #[test]
    fn test_validate_paired_two_files_accepted() {
        // Regression guard: the 2-file golden path must keep working.
        let cli = Cli::parse_from(["trim_galore", "--paired", R1, R2]);
        cli.validate().expect("two-file paired-end should validate");
    }

    #[test]
    fn test_validate_clock_still_requires_exactly_2() {
        // Regression guard: --clock must not have been widened by the
        // multi-pair --paired relaxation.
        let cli = Cli::parse_from(["trim_galore", "--clock", R1, R2, ALT_R1, ALT_R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("clock requires exactly 2"),
            "expected --clock strict-2 rejection, got: {err}"
        );
    }
}
