//! Command-line argument parsing and validation.

use clap::Parser;
use std::path::PathBuf;

/// Trim Galore - Oxidized Edition: A fast, single-pass NGS adapter and quality trimmer.
///
/// Drop-in replacement for Trim Galore, rewritten in Rust. Matches v0.6.x outputs
/// for the core feature set and extends it with poly-G / generic poly-A auto-trimming
/// and other additions. Compatible with MultiQC and existing pipelines.
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
    /// Supports A{N} shorthand for repeated single bases (e.g., -a A{10} → AAAAAAAAAA).
    /// For multiple adapters, repeat -a (e.g., -a SEQ1 -a SEQ2) or use "file:adapters.fa".
    #[clap(short = 'a', long = "adapter")]
    pub adapter: Vec<String>,

    /// Optional adapter sequence for Read 2 (paired-end only).
    /// Auto-set by --small_rna and --bgiseq presets.
    /// Supports A{N} shorthand for repeated single bases (e.g., -a2 T{150} → 150 T's).
    /// For multiple adapters, repeat -a2 (e.g., -a2 SEQ1 -a2 SEQ2) or use "file:adapters.fa".
    #[clap(long = "adapter2", alias = "a2")]
    pub adapter2: Vec<String>,

    /// Use Illumina universal adapter (AGATCGGAAGAGC). Also the auto-detect fallback.
    #[clap(long = "illumina", conflicts_with_all = &["nextera", "small_rna", "stranded_illumina", "bgiseq"])]
    pub illumina: bool,

    /// Use Nextera transposase adapter (CTGTCTCTTATA).
    #[clap(long = "nextera", conflicts_with_all = &["illumina", "small_rna", "stranded_illumina", "bgiseq"])]
    pub nextera: bool,

    /// Use Illumina Small RNA adapter (TGGAATTCTCGG).
    /// Also lowers --length default to 18 and sets --adapter2 (GATCGTCGGACT, Illumina small RNA 5').
    #[clap(long = "small_rna", conflicts_with_all = &["illumina", "nextera", "stranded_illumina", "bgiseq"])]
    pub small_rna: bool,

    /// Use Illumina Stranded mRNA adapter (ACTGTCTCTTATA).
    /// Not covered by auto-detection — must be set explicitly.
    #[clap(long = "stranded_illumina", conflicts_with_all = &["illumina", "nextera", "small_rna", "bgiseq"])]
    pub stranded_illumina: bool,

    /// Use BGI/DNBSEQ adapter. Sets --adapter2 for Read 2. Also probed by auto-detection.
    #[clap(long = "bgiseq", conflicts_with_all = &["illumina", "nextera", "small_rna", "stranded_illumina"])]
    pub bgiseq: bool,

    /// Paired-end mode. Accepts an even number of input files as consecutive R1/R2 pairs
    /// (e.g. --paired R1.fq.gz R2.fq.gz, or a glob matching multiple samples).
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
    /// In paired-end mode, both reads must pass; see --retain_unpaired to keep single survivors.
    #[clap(long = "length")]
    pub length: Option<usize>,

    /// Maximum allowed sequence length (discard reads longer than this).
    /// Typically only useful for smallRNA-seq to remove non-small-RNA reads.
    #[clap(long = "max_length")]
    pub max_length: Option<usize>,

    /// Maximum number of N bases allowed in a read.
    /// Integer: absolute count. Decimal (0-1): fraction of read length.
    /// In paired-end mode, either read over the limit removes the whole pair.
    #[clap(long = "max_n")]
    pub max_n: Option<f64>,

    /// Trim N bases from both ends of reads. Suppressed under --rrbs (matches Perl v0.6.x).
    #[clap(long = "trim-n", alias = "trim_n")]
    pub trim_n: bool,

    /// Remove N bases from the 5' end of Read 1. Useful for removing 5' quality-bias regions.
    #[clap(long = "clip_R1", alias = "clip_r1")]
    pub clip_r1: Option<usize>,

    /// Remove N bases from the 5' end of Read 2 (paired-end only).
    /// For paired-end bisulfite-seq, the end-repair step can introduce methylation bias; see Bismark User Guide.
    #[clap(long = "clip_R2", alias = "clip_r2")]
    pub clip_r2: Option<usize>,

    /// Remove N bases from the 3' end of Read 1, after adapter/quality trimming.
    #[clap(long = "three_prime_clip_R1", alias = "three_prime_clip_r1")]
    pub three_prime_clip_r1: Option<usize>,

    /// Remove N bases from the 3' end of Read 2 (paired-end only), after adapter/quality trimming.
    #[clap(long = "three_prime_clip_R2", alias = "three_prime_clip_r2")]
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

    /// Output directory for trimmed files. Created if it doesn't exist.
    #[clap(short = 'o', long = "output_dir")]
    pub output_dir: Option<PathBuf>,

    /// Custom basename for output files (replaces input filename stem).
    /// Only valid for a single file (single-end) or a single pair (paired-end).
    #[clap(long = "basename")]
    pub basename: Option<String>,

    /// Do not gzip-compress output files. Forces plain output regardless of
    /// input compression. By default, output compression mirrors the input
    /// (plain → plain, .gz → .gz; matches Perl v0.6.x behaviour).
    #[clap(long = "dont_gzip")]
    pub dont_gzip: bool,

    /// Reorder reads in the gzip output so reads sharing a canonical 16-mer
    /// minimizer land adjacent, letting gzip's dictionary find longer
    /// redundant runs. Typical saving: 16–55% on short-read FASTQ. Records
    /// are byte-identical — only their order on disk changes. Pair lockstep
    /// is preserved.
    ///
    /// Requires `--cores >= 2` and gzip output. Intended for short-read
    /// FASTQ; long-read inputs (ONT, PacBio) are unlikely to compress
    /// better. Combine with `--compression <N>` to trade speed against
    /// output size, and `--memory <N>` to grow the per-gzip-member sort
    /// run.
    #[clap(long = "clumpify")]
    pub clumpify: bool,

    /// Gzip compression level for output FASTQ (1–9). Default: 1 (fast,
    /// 75% larger files). Pass `--compression 6` for the gzip(1) default
    /// or `--compression 9` for archival use. Most useful in combination
    /// with `--clumpify`, where reordering plus a higher level
    /// compounds for substantially smaller output.
    #[clap(
        long = "compression",
        default_value_t = crate::fastq::DEFAULT_GZIP_LEVEL,
        value_parser = clap::value_parser!(u32).range(1..=9),
    )]
    pub compression: u32,

    /// Total memory budget for Trim Galore (e.g. `4G`, `512M`, `2G`).
    /// Currently used only by `--clumpify` for bin buffer sizing — bigger
    /// budget → bigger gzip members → better compression, up to a limit
    /// roughly equal to the uncompressed input size. Default: `4G`.
    /// Resolved bin layout and predicted peak RSS are printed at startup
    /// when `--clumpify` is set.
    #[clap(long = "memory", default_value = "4G")]
    pub memory: String,

    /// Suppress the trimming report.
    #[clap(long = "no_report_file")]
    pub no_report_file: bool,

    /// Retain unpaired reads when the mate is too short (paired-end only).
    /// Cutoff via --length_1 / --length_2 (default 35 each).
    #[clap(long = "retain_unpaired")]
    pub retain_unpaired: bool,

    /// Minimum length for unpaired Read 1 (with --retain_unpaired).
    #[clap(short = 'r', long = "length_1", default_value = "35", alias = "r1")]
    pub length_1: usize,

    /// Minimum length for unpaired Read 2 (with --retain_unpaired).
    #[clap(long = "length_2", default_value = "35", alias = "r2")]
    pub length_2: usize,

    /// Add clipped sequences to read IDs for --clip_R1/R2, --three_prime_clip_R1/R2, and --hardtrim5/3.
    /// Appends :clip5:SEQ and/or :clip3:SEQ to the read ID (each half only when that side was clipped). Useful for UMI handling.
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
    /// Do not use with Tecan Ovation RRBS kits — those use a diversity-trimming step instead.
    #[clap(long = "rrbs")]
    pub rrbs: bool,

    /// Non-directional RRBS libraries. Reads starting with CAA or CGA get 2bp
    /// trimmed from the 5' end. Requires --rrbs.
    /// Unlike directional --rrbs, does not auto-set --clip_R2 2 in paired-end mode.
    #[clap(long = "non_directional", requires = "rrbs")]
    pub non_directional: bool,

    /// Run FastQC on the trimmed output files (built in via the bundled
    /// fastqc-rust library; no external Java or FastQC binary needed).
    /// Produces FastQC 0.12.1-compatible *_fastqc.html / *_fastqc.zip artifacts.
    #[clap(long = "fastqc")]
    pub fastqc: bool,

    /// Additional arguments to pass to FastQC. Implies --fastqc.
    /// Common flags are translated to the bundled engine: --nogroup, --expgroup,
    /// --quiet, --svg, --nano, --nofilter, --casava, -t/--threads, -o/--outdir.
    /// Unrecognised flags emit a warning and are ignored.
    #[clap(long = "fastqc_args", allow_hyphen_values = true)]
    pub fastqc_args: Option<String>,

    /// Number of worker threads for parallel processing (default: 1).
    /// Values > 1 run trimming and gzip compression across multiple threads.
    /// Near-linear speedup through at least 16 cores; diminishing returns beyond ~20
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
    /// Bypasses adapter/quality trimming entirely. Output filenames end in .<N>bp_5prime.fq(.gz).
    #[clap(long = "hardtrim5")]
    pub hardtrim5: Option<usize>,

    /// Hard-trim to keep only the last N bases from the 3' end.
    /// Bypasses adapter/quality trimming entirely. Output filenames end in .<N>bp_3prime.fq(.gz).
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

    /// [Deprecated] v2.0 emits only essential progress output; use shell redirection if quieter output is needed.
    #[clap(long = "suppress_warn", hide = true)]
    pub suppress_warn: bool,

    /// [Deprecated] Reports are generated by default.
    #[clap(long = "report", hide = true)]
    pub report: bool,

    /// [Deprecated] The v2.0 single-pass architecture has no quality-trim intermediate file to keep.
    #[clap(long = "keep", hide = true)]
    pub keep: bool,

    /// Easter egg (no-op).
    #[clap(long = "hulu", hide = true)]
    pub hulu: bool,
}

/// Rewrite Perl-era multi-character short flags (`-r1`, `-r2`, `-a2`) as
/// their clap-compatible long-alias forms (`--r1`, `--r2`, `--a2`) before
/// parsing.
///
/// Clap derives single-character short flags only, so e.g. `-r1 40` would
/// parse as `-r=1` with `40` becoming a stray positional, producing a
/// confusing "odd count of input files" error. `-a2 SEQ` would similarly
/// parse as `-a=2` with `SEQ` becoming a positional input file. This
/// pre-parse hook transparently rewrites the exact tokens so Perl-era
/// invocations keep working.
///
/// Only exact-match tokens are rewritten — `-r10` (legitimate clap
/// `-r=10`) and any other value-suffixed form pass through unchanged.
pub fn rewrite_perl_short_flags<I>(args: I) -> Vec<String>
where
    I: IntoIterator<Item = String>,
{
    args.into_iter()
        .map(|a| {
            if a == "-r1" || a.starts_with("-r1=") {
                format!("--r1{}", &a[3..])
            } else if a == "-r2" || a.starts_with("-r2=") {
                format!("--r2{}", &a[3..])
            } else if a == "-a2" || a.starts_with("-a2=") {
                format!("--a2{}", &a[3..])
            } else {
                a
            }
        })
        .collect()
}

impl Cli {
    /// Shared validation for any paired-end mode (`--paired`, `--clock`,
    /// `--implicon`) that takes input files in pairwise (R1, R2, R1, R2, …)
    /// order. Checks:
    ///   1. Even count of input files.
    ///   2. Within each pair, R1 ≠ R2 byte-equal (matches Perl's
    ///      `$ARGV[$i] eq $ARGV[$i+1]` check at `trim_galore:3208`; does not
    ///      follow symlinks or canonicalise).
    ///   3. Across pairs, no duplicate pair (catches accidental copy-paste
    ///      and emits a precise error rather than the case-insensitive
    ///      output-collision pre-flight's APFS/NTFS message).
    ///
    /// `mode_label` is used in the user-facing error string, e.g.
    /// `"Paired-end"`, `"--clock"`, `"--implicon"`.
    fn validate_paired_input(&self, mode_label: &str) -> anyhow::Result<()> {
        if !self.input.len().is_multiple_of(2) {
            anyhow::bail!(
                "{} mode requires an even number of input files (R1/R2 pairs), got {}",
                mode_label,
                self.input.len()
            );
        }
        for chunk in self.input.chunks(2) {
            if chunk[0] == chunk[1] {
                anyhow::bail!(
                    "Read 1 and Read 2 appear to be the same file: {}. \
                     Did you mean to pass distinct R1 and R2 files?",
                    chunk[0].display()
                );
            }
        }
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
        Ok(())
    }

    /// Validate CLI arguments after parsing.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.paired {
            // `#[clap(required = true)]` on `input` guarantees at least one file
            // reaches validate(), so no is_empty() check is needed.
            self.validate_paired_input("Paired-end")?;
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

        if let Some(val) = self.nextseq
            && (val == 0 || val >= 200)
        {
            anyhow::bail!(
                "NextSeq quality cutoff must be between 1 and 199, got {}",
                val
            );
        }

        if let Some(threshold) = self.consider_already_trimmed
            && threshold > 10000
        {
            anyhow::bail!(
                "consider_already_trimmed value must be between 0 and 10000, got {}",
                threshold
            );
        }

        if self.times == 0 || self.times > 10 {
            anyhow::bail!("--times/-n must be between 1 and 10, got {}", self.times);
        }

        if self.cores == 0 {
            anyhow::bail!("--cores must be at least 1");
        }

        if self.clumpify {
            if self.cores < 2 {
                anyhow::bail!(
                    "--clumpify requires --cores >= 2 (the bin dispatcher feeds parallel workers)"
                );
            }
            if self.dont_gzip {
                anyhow::bail!(
                    "--clumpify and --dont_gzip are mutually exclusive (clumping plain text is pointless)"
                );
            }
            if self.clock {
                anyhow::bail!("--clumpify is not yet supported with --clock");
            }
            if self.implicon.is_some() {
                anyhow::bail!("--clumpify is not yet supported with --implicon");
            }
            if self.hardtrim5.is_some() {
                anyhow::bail!("--clumpify is not yet supported with --hardtrim5");
            }
            if self.hardtrim3.is_some() {
                anyhow::bail!("--clumpify is not yet supported with --hardtrim3");
            }
            if self.demux.is_some() {
                anyhow::bail!("--clumpify is not yet supported with --demux");
            }
            // Resolve & validate the layout up front so misconfigured runs
            // fail before any I/O. The same call is repeated at dispatch
            // time so the resolved values are passed through to workers.
            let memory_bytes = crate::clump::parse_memory_size(&self.memory)
                .map_err(|e| anyhow::anyhow!("--memory: {e}"))?;
            crate::clump::resolve_layout(memory_bytes, self.cores)?;
        }

        if let Some(n) = self.hardtrim5
            && (n == 0 || n >= 1000)
        {
            anyhow::bail!("--hardtrim5 must be between 1 and 999, got {}", n);
        }
        if let Some(n) = self.hardtrim3
            && (n == 0 || n >= 1000)
        {
            anyhow::bail!("--hardtrim3 must be between 1 and 999, got {}", n);
        }
        if self.clock {
            self.validate_paired_input("--clock")?;
        }
        if self.implicon.is_some() {
            self.validate_paired_input("--implicon")?;
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

    /// Resolve the gzip compression level for output FASTQ. Always
    /// `--compression`; defaults to `DEFAULT_GZIP_LEVEL` (1).
    pub fn gzip_level(&self) -> u32 {
        self.compression
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

    // ── Multi-pair widening for --clock and --implicon ──
    // (Replaces the earlier "strict-2" regression guard. Specialty
    // run-and-exit modes now share the same pairwise validation as
    // --paired itself.)

    #[test]
    fn test_validate_clock_two_pairs_accepted() {
        let cli = Cli::parse_from(["trim_galore", "--clock", R1, R2, ALT_R1, ALT_R2]);
        cli.validate()
            .expect("two distinct pairs should validate under --clock");
    }

    #[test]
    fn test_validate_clock_odd_count_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--clock", R1, R2, ALT_R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--clock") && err.contains("even number"),
            "expected --clock even-count rejection, got: {err}"
        );
    }

    #[test]
    fn test_validate_clock_r1_equal_r2_within_pair_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--clock", R1, R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("appear to be the same file"),
            "expected R1==R2 rejection under --clock, got: {err}"
        );
    }

    #[test]
    fn test_validate_clock_duplicate_pair_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--clock", R1, R2, R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("duplicate of pair"),
            "expected duplicate-pair rejection under --clock, got: {err}"
        );
    }

    #[test]
    fn test_validate_implicon_two_pairs_accepted() {
        let cli = Cli::parse_from(["trim_galore", "--implicon", R1, R2, ALT_R1, ALT_R2]);
        cli.validate()
            .expect("two distinct pairs should validate under --implicon");
    }

    #[test]
    fn test_validate_implicon_odd_count_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--implicon", R1, R2, ALT_R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--implicon") && err.contains("even number"),
            "expected --implicon even-count rejection, got: {err}"
        );
    }

    #[test]
    fn test_validate_implicon_duplicate_pair_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--implicon", R1, R2, R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("duplicate of pair"),
            "expected duplicate-pair rejection under --implicon, got: {err}"
        );
    }

    // ── Perl-migration short-flag rewrite (-r1 → --r1, -r2 → --r2) ──

    fn rewrite(args: &[&str]) -> Vec<String> {
        super::rewrite_perl_short_flags(args.iter().map(|s| s.to_string()))
    }

    #[test]
    fn test_rewrite_r1_bare() {
        assert_eq!(
            rewrite(&["trim_galore", "-r1", "40"]),
            vec!["trim_galore", "--r1", "40"]
        );
    }

    #[test]
    fn test_rewrite_r2_bare() {
        assert_eq!(
            rewrite(&["trim_galore", "-r2", "35"]),
            vec!["trim_galore", "--r2", "35"]
        );
    }

    #[test]
    fn test_rewrite_r1_equals_form() {
        assert_eq!(
            rewrite(&["trim_galore", "-r1=40"]),
            vec!["trim_galore", "--r1=40"]
        );
    }

    #[test]
    fn test_rewrite_r2_equals_form() {
        assert_eq!(
            rewrite(&["trim_galore", "-r2=35"]),
            vec!["trim_galore", "--r2=35"]
        );
    }

    #[test]
    fn test_rewrite_leaves_r_alone() {
        // -r 40 is valid clap short; must not be disturbed.
        assert_eq!(
            rewrite(&["trim_galore", "-r", "40"]),
            vec!["trim_galore", "-r", "40"]
        );
    }

    #[test]
    fn test_rewrite_leaves_r10_alone() {
        // -r10 is clap's short-with-value syntax (-r=10); must not be rewritten.
        assert_eq!(
            rewrite(&["trim_galore", "-r10"]),
            vec!["trim_galore", "-r10"]
        );
    }

    #[test]
    fn test_rewrite_leaves_r20_alone() {
        // -r20 is clap's short-with-value (-r=20); not a Perl `-r2` + value.
        assert_eq!(
            rewrite(&["trim_galore", "-r20"]),
            vec!["trim_galore", "-r20"]
        );
    }

    #[test]
    fn test_rewrite_leaves_unrelated_alone() {
        assert_eq!(
            rewrite(&["trim_galore", "--paired", "-a", "AGCT", "-o", "outdir"]),
            vec!["trim_galore", "--paired", "-a", "AGCT", "-o", "outdir"]
        );
    }

    #[test]
    fn test_rewrite_end_to_end_via_parse_from() {
        // Verify that after rewriting, Cli::parse_from successfully parses
        // -r1 / -r2 style invocations (this is the whole point of the rewrite).
        let args = rewrite(&[
            "trim_galore",
            "--paired",
            "--retain_unpaired",
            "-r1",
            "40",
            "-r2",
            "30",
            "test_files/BS-seq_10K_R1.fastq.gz",
            "test_files/BS-seq_10K_R2.fastq.gz",
        ]);
        let cli = Cli::parse_from(args);
        assert_eq!(cli.length_1, 40);
        assert_eq!(cli.length_2, 30);
    }

    #[test]
    fn test_rewrite_a2_bare() {
        assert_eq!(
            rewrite(&["trim_galore", "-a2", "GCAT"]),
            vec!["trim_galore", "--a2", "GCAT"]
        );
    }

    #[test]
    fn test_rewrite_a2_equals_form() {
        assert_eq!(
            rewrite(&["trim_galore", "-a2=GCAT"]),
            vec!["trim_galore", "--a2=GCAT"]
        );
    }

    #[test]
    fn test_rewrite_leaves_a10_alone() {
        // -a10 is clap's short-with-value (-a=10) — not a Perl `-a2` construct.
        // `10` isn't a valid DNA sequence but that's for validation to catch,
        // not for the rewrite to mangle.
        assert_eq!(
            rewrite(&["trim_galore", "-a10"]),
            vec!["trim_galore", "-a10"]
        );
    }

    #[test]
    fn test_rewrite_a2_end_to_end_via_parse_from() {
        let args = rewrite(&[
            "trim_galore",
            "--paired",
            "-a",
            "AGCT",
            "-a2",
            "GCAT",
            "-a2",
            "AAAA",
            "test_files/BS-seq_10K_R1.fastq.gz",
            "test_files/BS-seq_10K_R2.fastq.gz",
        ]);
        let cli = Cli::parse_from(args);
        assert_eq!(cli.adapter, vec!["AGCT"]);
        assert_eq!(cli.adapter2, vec!["GCAT", "AAAA"]);
    }

    /// Perl `trim_galore` accepts the lowercase clip-flag spellings
    /// (`--clip_r1` / `--clip_r2` / `--three_prime_clip_r1` /
    /// `--three_prime_clip_r2`) alongside the uppercase forms. The Rust port
    /// historically only matched the uppercase canonical, breaking every
    /// Perl-era pipeline using the lowercase spelling. Regression for #242.
    #[test]
    fn test_clip_flags_accept_lowercase_aliases() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--clip_r1",
            "5",
            "--clip_r2",
            "6",
            "--three_prime_clip_r1",
            "7",
            "--three_prime_clip_r2",
            "8",
            R1,
            R2,
        ]);
        assert_eq!(cli.clip_r1, Some(5));
        assert_eq!(cli.clip_r2, Some(6));
        assert_eq!(cli.three_prime_clip_r1, Some(7));
        assert_eq!(cli.three_prime_clip_r2, Some(8));

        // The canonical uppercase forms must of course still work. Mix a few
        // to confirm both aliases resolve to the same field.
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--clip_R1",
            "1",
            "--clip_r2",
            "2",
            "--three_prime_clip_R1",
            "3",
            "--three_prime_clip_r2",
            "4",
            R1,
            R2,
        ]);
        assert_eq!(cli.clip_r1, Some(1));
        assert_eq!(cli.clip_r2, Some(2));
        assert_eq!(cli.three_prime_clip_r1, Some(3));
        assert_eq!(cli.three_prime_clip_r2, Some(4));
    }

    // ── --clumpify / --compression validation ────────────────────────────

    #[test]
    fn test_clumpify_requires_cores_at_least_two() {
        let cli = Cli::parse_from(["trim_galore", "--clumpify", "--cores", "1", R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--clumpify requires --cores >= 2"),
            "got: {err}"
        );
    }

    #[test]
    fn test_clumpify_rejects_dont_gzip() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "2",
            "--dont_gzip",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--dont_gzip"), "got: {err}");
    }

    #[test]
    fn test_clumpify_rejects_clock() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "2",
            "--clock",
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--clock"), "got: {err}");
    }

    #[test]
    fn test_clumpify_rejects_implicon() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "2",
            "--implicon=8",
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--implicon"), "got: {err}");
    }

    #[test]
    fn test_clumpify_rejects_hardtrim5() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "2",
            "--hardtrim5",
            "30",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--hardtrim5"), "got: {err}");
    }

    #[test]
    fn test_clumpify_accepts_paired() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "4",
            "--paired",
            R1,
            R2,
        ]);
        cli.validate()
            .expect("clumpify + paired + cores=4 should validate");
    }

    #[test]
    fn test_clumpify_rejects_garbage_memory() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "2",
            "--memory",
            "garbage",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--memory"), "got: {err}");
    }

    #[test]
    fn test_clumpify_rejects_too_small_memory_for_cores() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "16",
            "--memory",
            "64M",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--memory"), "got: {err}");
    }

    #[test]
    fn test_compression_defaults_to_one() {
        let cli = Cli::parse_from(["trim_galore", R1]);
        assert_eq!(cli.compression, 1);
        assert_eq!(cli.gzip_level(), 1);
    }

    #[test]
    fn test_compression_explicit_level() {
        let cli = Cli::parse_from(["trim_galore", "--compression", "9", R1]);
        assert_eq!(cli.compression, 9);
        assert_eq!(cli.gzip_level(), 9);
    }

    #[test]
    fn test_compression_rejects_out_of_range() {
        let result = Cli::try_parse_from(["trim_galore", "--compression", "10", R1]);
        assert!(result.is_err(), "level 10 should be rejected by clap");
    }

    #[test]
    fn test_clumpify_with_compression_six() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--compression",
            "6",
            "--cores",
            "2",
            R1,
        ]);
        cli.validate()
            .expect("clumpify --compression 6 should validate");
        assert!(cli.clumpify);
        assert_eq!(cli.compression, 6);
    }
}
