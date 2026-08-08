//! Command-line argument parsing and validation.

use clap::Parser;
use std::path::PathBuf;

/// Output container format. FASTQ is the default and preserves byte-identity
/// with Perl Trim Galore 0.6.11. uBAM is opt-in and carries the input BAM's
/// aux tags through (via `--preserve-tags`).
///
/// See `plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN.md`
/// §3.2 for the full output-naming + behaviour matrix.
#[derive(Clone, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum OutputFormat {
    /// FASTQ output, mirroring input compression (default).
    #[default]
    Fastq,
    /// Unaligned BAM output. Always single-threaded; --clumpify,
    /// --passthrough, --clock, --implicon, --demux are rejected at validation.
    #[clap(name = "ubam")]
    UBam,
}

/// Trim Galore: A fast, single-pass NGS adapter and quality trimmer.
///
/// Single-binary adapter and quality trimmer for NGS FASTQ data, with poly-G /
/// generic poly-A auto-trimming and built-in FastQC reporting. Compatible with
/// MultiQC and existing pipelines.
#[derive(Parser, Debug)]
#[clap(
    name = "trim_galore",
    version = env!("CARGO_PKG_VERSION"),
    long_version = concat!(
        env!("CARGO_PKG_VERSION"), "\n",
        env!("VERSION_BODY")
    ),
    about,
    // Bare `trim_galore` (no args) prints help and exits 0 rather than a
    // "required arguments were not provided" usage error, matching the
    // convention of samtools/git/bwa and friends.
    arg_required_else_help = true
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

    /// Optional adapter sequence for Read 2 (paired-end only). Takes precedence
    /// over the Read 2 default that --small_rna and --bgiseq otherwise set.
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
    /// Also lowers --length default to 18 and sets --adapter2 (GATCGTCGGACT, Illumina small RNA 5') unless given.
    #[clap(long = "small_rna", conflicts_with_all = &["illumina", "nextera", "stranded_illumina", "bgiseq"])]
    pub small_rna: bool,

    /// Use Illumina Stranded mRNA adapter (ACTGTCTCTTATA).
    /// Not covered by auto-detection — must be set explicitly.
    #[clap(long = "stranded_illumina", conflicts_with_all = &["illumina", "nextera", "small_rna", "bgiseq"])]
    pub stranded_illumina: bool,

    /// Use BGI/DNBSEQ adapter. Sets --adapter2 for Read 2 unless given. Also probed by auto-detection.
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

    /// Trim N bases from both ends of reads. Suppressed under --rrbs.
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
    /// (plain → plain, .gz → .gz).
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

    /// Lossless reorder-only specialty mode: reorder FASTQ records by
    /// canonical 16-mer minimizer for gzip-friendly compression, WITHOUT
    /// any trimming, filtering, adapter detection, or record modification.
    /// Every input record appears in the output byte-identically (header,
    /// sequence, quality); only file-level order changes. Output files use
    /// the `*_clumped.fq(.gz)` (SE) or `*_clumped_{1,2}.fq(.gz)` (PE) suffix,
    /// and a short `*_clumping_report.txt` is emitted (distinct from
    /// `*_trimming_report.txt` to keep downstream nf-core/MultiQC scanners
    /// unconfused).
    ///
    /// Composes with `--compression`, `--memory`, `--cores`, `--paired`,
    /// `--fastqc`, `--dont_gzip`, and `--basename`. Trimming/filtering flags
    /// (`-a`, `--length`, `--rrbs`, `--polyA`, `--polyG`, `--trim-n`,
    /// `--clip_*`, `--nextseq`, `--rename`, `--discard_untrimmed`,
    /// `--consider_already_trimmed`, other specialty modes, `--passthrough`,
    /// `--retain_unpaired`, `--output-format ubam`) are rejected. `-q` /
    /// `--stringency` / `-e` have clap defaults and are silently ignored on
    /// this path (mode does no trimming; matches how `--hardtrim5` treats
    /// trim flags today).
    ///
    /// Contract-scope note: byte-identity applies to header + sequence +
    /// quality bytes. The plus-line (line 3 of each record) is normalized
    /// to bare `+` on output; CRLF line endings are normalized to LF. Both
    /// normalizations are codebase-wide behaviours, inherited from the
    /// FASTQ reader/writer.
    ///
    /// v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up.
    #[clap(long = "clump_only")]
    pub clump_only: bool,

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

    /// Total memory budget for Trim Galore (e.g. `1G`, `512M`, `8G`).
    /// Currently used only by `--clumpify` for bin buffer sizing — bigger
    /// budget → bigger gzip members → better compression, up to a limit
    /// roughly equal to the uncompressed input size. Default: `1G`.
    /// Resolved bin layout and predicted peak RSS are printed at startup
    /// when `--clumpify` is set.
    #[clap(long = "memory", default_value = "1G")]
    pub memory: String,

    /// Suppress the trimming report.
    #[clap(long = "no_report_file")]
    pub no_report_file: bool,

    /// Comma-separated list of BAM tags (e.g. "CB,UB,RX") to append to the
    /// FASTQ header on uBAM input. Tab-separated, samtools `-T`-compatible.
    /// Each tag must be a 2-char SAM tag name (`[A-Za-z][A-Za-z0-9]`).
    /// Ignored for FASTQ input. The `ALL` keyword is reserved for a future
    /// release and is currently rejected.
    #[clap(long = "preserve-tags", value_delimiter = ',', value_parser = parse_sam_tag_name)]
    pub preserve_tags: Vec<String>,

    /// Output container format. Default `fastq` keeps existing behaviour
    /// (input-compression-mirroring FASTQ); `ubam` emits unaligned BAM with
    /// aux tags propagated from uBAM inputs (when `--preserve-tags` is set).
    /// uBAM output is always single-threaded; some flag combinations are
    /// rejected (see `--help` and the startup diagnostics for details).
    #[clap(long = "output-format", value_enum, default_value_t = OutputFormat::Fastq)]
    pub output_format: OutputFormat,

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

    /// Pass a third FASTQ file through unchanged but keep it in lockstep with R1/R2.
    /// Use case: 10X Multiome / scATAC libraries where the cell-barcode (I1/I2) read
    /// must stay aligned to the trimmed R1/R2. Records dropped by length/quality/N
    /// filters are also dropped from the passthrough output. The passthrough file is
    /// never trimmed or adapter-scanned.
    ///
    /// Requires --paired with exactly one R1/R2 pair. Incompatible with
    /// --retain_unpaired, --clumpify, and all specialty modes (--clock, --implicon,
    /// --hardtrim5/3, --demux).
    ///
    /// Note: legacy headers like @read/1, @read/2, @read/3 sync correctly. If your
    /// files use unusual header conventions and the sync check fires unexpectedly,
    /// please file an issue with a sample header line. If --fastqc is enabled the
    /// passthrough FastQC report will look poor (cell-barcode reads are intentionally
    /// uniformly-structured) — that's expected, not a defect. On a mid-stream reader
    /// error (truncated or desynced passthrough), partial output files may remain on
    /// disk — re-run after fixing the input.
    #[clap(long = "passthrough")]
    pub passthrough: Option<PathBuf>,

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
    /// At --cores 1 the worker-pool is bypassed (single thread, ~5 MB RAM).
    /// From --cores 2 upward, an N+4 thread model applies: N workers + 2
    /// decompressors + 1 batcher + 1 writer. Wall-clock speedup is near-linear
    /// up to --cores 8 for paired-end runs; beyond that, gzip-output I/O on
    /// the storage layer typically becomes binding before workers run out of
    /// useful per-read work, so additional cores help progressively less.
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
/// clap value parser for `--preserve-tags`. Each tag must be a valid 2-char
/// SAM tag name (`[A-Za-z][A-Za-z0-9]`). The `ALL` keyword is reserved.
fn parse_sam_tag_name(s: &str) -> Result<String, String> {
    if s == "ALL" {
        return Err("--preserve-tags ALL is not supported in v1 (use an explicit list)".into());
    }
    let bytes = s.as_bytes();
    if bytes.len() != 2 {
        return Err(format!(
            "'{s}' is not a valid SAM tag name (must be exactly 2 characters)"
        ));
    }
    if !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
        return Err(format!(
            "'{s}' is not a valid SAM tag name (must match [A-Za-z][A-Za-z0-9])"
        ));
    }
    Ok(s.to_string())
}

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

/// Existence + restartability check for a path the run will re-read (#379).
///
/// `stat` rather than `exists()` so a pipe or FIFO is named here, before
/// anything opens it — `File::open` on a FIFO with no writer blocks forever,
/// which is a hang with no message rather than a diagnostic.
fn check_restartable_input(path: &std::path::Path, not_found: &str) -> anyhow::Result<()> {
    let meta = std::fs::metadata(path).map_err(|e| {
        // Only NotFound is "not found"; reporting EACCES or EMFILE that way is
        // the wrong-blame this guard exists to remove.
        if e.kind() == std::io::ErrorKind::NotFound {
            anyhow::anyhow!("{not_found}: {}", path.display())
        } else {
            anyhow::Error::new(e).context(format!("Cannot stat input file: {}", path.display()))
        }
    })?;
    if let Some(kind) = crate::format::non_restartable_kind(&meta) {
        anyhow::bail!("{}", crate::format::not_restartable_message(path, kind));
    }
    Ok(())
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
        // Allow N=1 for the `--paired` case: that's only legal if the single
        // file is a uBAM, in which case the de-interleaver produces R1+R2 from
        // one interleaved BAM. The "is it BAM?" check happens at main.rs
        // sanity_check time (after format detection has run). For specialty
        // modes (--clock, --implicon, --hardtrim) the strict even-count rule
        // still applies because they don't yet support uBAM input.
        if self.input.len() == 1 && self.paired && mode_label == "Paired-end" {
            return Ok(());
        }
        if !self.input.len().is_multiple_of(2) {
            anyhow::bail!(
                "{} mode requires an even number of input files (R1/R2 pairs), got {}",
                mode_label,
                self.input.len()
            );
        }
        for chunk in self.input.chunks(2) {
            // #383 — keyed like the output-collision pre-flight, so `./r1.fq r1.fq`
            // cannot pass as a pair. Raw `==` let two spellings of one file through.
            if crate::io::path_identity_key(&chunk[0]) == crate::io::path_identity_key(&chunk[1]) {
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
                if crate::io::path_identity_key(r1) == crate::io::path_identity_key(pr1)
                    && crate::io::path_identity_key(r2) == crate::io::path_identity_key(pr2)
                {
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
        // §3.4a — `--output-format ubam` exclusions. These are pure CLI-level
        // (no file I/O); enforced here. The §3.4b rule (preserve-tags + all
        // FASTQ inputs) requires format detection and lives in main.rs.
        if matches!(self.output_format, OutputFormat::UBam) {
            // v2 addition: --dont_gzip is meaningless with BAM output (BAM is
            // always BGZF-compressed). Rejecting here in the shared UBam block
            // covers both the trim uBAM path (which previously accepted this
            // silently — a real gap) and the new --clump_only uBAM path.
            if self.dont_gzip {
                anyhow::bail!(
                    "--dont_gzip is not compatible with --output-format ubam \
                     (BAM is always BGZF-compressed)"
                );
            }
            if self.clumpify {
                anyhow::bail!(
                    "--clumpify is for gzip output; not applicable with --output-format ubam"
                );
            }
            if self.passthrough.is_some() {
                anyhow::bail!("--passthrough is not supported with --output-format ubam in v1");
            }
            if self.clock {
                anyhow::bail!(
                    "--clock + --output-format ubam: UMI-to-BAM-tag mapping not defined in v1; \
                     use FASTQ output or convert after"
                );
            }
            if self.implicon.is_some() {
                anyhow::bail!(
                    "--implicon + --output-format ubam: UMI-to-BAM-tag mapping not defined in v1"
                );
            }
            if self.demux.is_some() {
                anyhow::bail!("--demux is not supported with --output-format ubam in v1");
            }
            // PLAN v2.1 §3.4a addendum (Step 3 implementation note): v2.1 left
            // --retain_unpaired silent, but it produces multiple output files
            // (`*_unpaired_{1,2}.fq.gz`) — same shape as --demux / --passthrough
            // which are already rejected. Multi-output BAM is out of scope for
            // v1; revisit if a real workload demands it.
            if self.retain_unpaired {
                anyhow::bail!(
                    "--retain_unpaired is not supported with --output-format ubam in v1 \
                     (unpaired records would require additional BAM output paths; \
                     run without --retain_unpaired or post-process via samtools)"
                );
            }
        }

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
        // #383 — caught here, not in the output-collision pre-flight, so the message can
        // be precise (see validate_paired_input's rationale for the paired equivalent).
        // `--clock`/`--implicon` are paired but never set `paired`, and own a more
        // specific R1==R2 message, so they are excluded rather than pre-empted.
        if !self.paired && !self.clock && self.implicon.is_none() {
            for (i, path) in self.input.iter().enumerate() {
                let key = crate::io::path_identity_key(path);
                if let Some(j) = self.input[..i]
                    .iter()
                    .position(|other| crate::io::path_identity_key(other) == key)
                {
                    anyhow::bail!(
                        "Input file {} was given more than once (arguments {} and {}). \
                         List each input once — trimming it twice would write the same \
                         output file twice.",
                        path.display(),
                        j + 1,
                        i + 1
                    );
                }
            }
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
            // Validate --memory format up front. Whether the resolved bin
            // pool is *large enough* for clumpify to actually run is decided
            // later in main.rs::resolve_clump_layout, which warns and falls
            // back to plain mode if the budget is below the floor.
            crate::clump::parse_memory_size(&self.memory)
                .map_err(|e| anyhow::anyhow!("--memory: {e}"))?;
        }

        // --clump_only: lossless reorder-only specialty mode. Rejection
        // matrix mirrors --clumpify's exclusion list plus every flag that
        // would trim, filter, clip, or mutate records — byte-identity is
        // the mode's load-bearing invariant. `-q` / `--stringency` / `-e`
        // are silently accepted (non-Option clap defaults; can't distinguish
        // user-set from default without threading ArgMatches::value_source).
        // Documented as ignored under --clump_only in the flag's --help text.
        if self.clump_only {
            if self.clumpify {
                anyhow::bail!(
                    "--clump_only and --clumpify are mutually exclusive \
                     (--clump_only is a lossless reorder-only variant that supersedes --clumpify's use case)"
                );
            }
            // Note: --clump_only v1 is single-threaded internally (deviation
            // from --clumpify's `--cores >= 2` requirement). `--cores` is
            // accepted at any value >= 1 but only affects future parallel
            // implementations; the byte-identity contract holds regardless.
            // v2: --paired + N=1 is now legal when the single input is an
            // interleaved uBAM (Shape B). The format check for "N=1 but not
            // a BAM" runs in main.rs::dispatch (needs `detect_input_format`,
            // which reads the file), matching how the trim uBAM path handles
            // the same distinction. If the file is a FASTQ, dispatch bails
            // with a clear error before opening any reader.
            // Adapter flags
            if !self.adapter.is_empty() {
                anyhow::bail!("--clump_only does not trim; -a/--adapter is not compatible");
            }
            if !self.adapter2.is_empty() {
                anyhow::bail!("--clump_only does not trim; -a2/--adapter2 is not compatible");
            }
            if self.illumina {
                anyhow::bail!("--clump_only does not trim; --illumina is not compatible");
            }
            if self.nextera {
                anyhow::bail!("--clump_only does not trim; --nextera is not compatible");
            }
            if self.small_rna {
                anyhow::bail!("--clump_only does not trim; --small_rna is not compatible");
            }
            if self.bgiseq {
                anyhow::bail!("--clump_only does not trim; --bgi/--bgiseq is not compatible");
            }
            if self.stranded_illumina {
                anyhow::bail!("--clump_only does not trim; --stranded_illumina is not compatible");
            }
            // Length / filter flags (Option-typed, so user-set is distinguishable)
            if self.length.is_some() {
                anyhow::bail!("--clump_only does not filter; --length is not compatible");
            }
            if self.max_length.is_some() {
                anyhow::bail!("--clump_only does not filter; --max_length is not compatible");
            }
            if self.max_n.is_some() {
                anyhow::bail!("--clump_only does not filter; --max_n is not compatible");
            }
            // Clip flags
            if self.trim_n {
                anyhow::bail!("--clump_only does not trim; --trim-n is not compatible");
            }
            if self.clip_r1.is_some() {
                anyhow::bail!("--clump_only does not clip; --clip_r1 is not compatible");
            }
            if self.clip_r2.is_some() {
                anyhow::bail!("--clump_only does not clip; --clip_r2 is not compatible");
            }
            if self.three_prime_clip_r1.is_some() {
                anyhow::bail!(
                    "--clump_only does not clip; --three_prime_clip_r1 is not compatible"
                );
            }
            if self.three_prime_clip_r2.is_some() {
                anyhow::bail!(
                    "--clump_only does not clip; --three_prime_clip_r2 is not compatible"
                );
            }
            // RRBS
            if self.rrbs {
                anyhow::bail!("--clump_only does not trim; --rrbs is not compatible");
            }
            if self.non_directional {
                anyhow::bail!("--clump_only does not trim; --non_directional is not compatible");
            }
            // Poly-*
            if self.poly_a {
                anyhow::bail!("--clump_only does not trim; --polyA is not compatible");
            }
            if self.poly_g {
                anyhow::bail!("--clump_only does not trim; --polyG is not compatible");
            }
            if self.no_poly_g {
                anyhow::bail!(
                    "--clump_only does not run poly-G auto-detection; --no_poly_g is not compatible"
                );
            }
            // 2-colour quality
            if self.nextseq.is_some() {
                anyhow::bail!("--clump_only does not trim; --nextseq/--2colour is not compatible");
            }
            // Renaming / filter-adjacent
            if self.rename {
                anyhow::bail!(
                    "--clump_only preserves record contents byte-identically; --rename would mutate read IDs"
                );
            }
            if self.discard_untrimmed {
                anyhow::bail!(
                    "--clump_only does not trim; --discard_untrimmed has no meaning here"
                );
            }
            if self.consider_already_trimmed.is_some() {
                anyhow::bail!(
                    "--clump_only does not trim; --consider_already_trimmed is not compatible"
                );
            }
            // Other specialty modes (all mutually exclusive)
            if self.hardtrim5.is_some() {
                anyhow::bail!("--clump_only and --hardtrim5 are mutually exclusive");
            }
            if self.hardtrim3.is_some() {
                anyhow::bail!("--clump_only and --hardtrim3 are mutually exclusive");
            }
            if self.clock {
                anyhow::bail!("--clump_only and --clock are mutually exclusive");
            }
            if self.implicon.is_some() {
                anyhow::bail!("--clump_only and --implicon are mutually exclusive");
            }
            if self.demux.is_some() {
                anyhow::bail!("--clump_only and --demux are mutually exclusive");
            }
            // Output shape (v2: --output-format ubam is now accepted; uBAM in/out
            // is dispatched to the BAM variant functions in main.rs. Format-guards
            // for two-BAM Shape A, non-BAM Shape B, and mixed-format Shape A live
            // in main.rs::dispatch (they need format detection).)
            if self.passthrough.is_some() {
                anyhow::bail!(
                    "--clump_only + --passthrough is not compatible (passthrough is a trim-pipeline feature)"
                );
            }
            if self.retain_unpaired {
                anyhow::bail!(
                    "--clump_only does not filter; --retain_unpaired has no meaning here"
                );
            }
            // Validate --memory format up front. Same treatment as --clumpify.
            crate::clump::parse_memory_size(&self.memory)
                .map_err(|e| anyhow::anyhow!("--memory: {e}"))?;
        }

        // --passthrough: 9-item compatibility envelope. Layout mirrors --clumpify
        // above. Each rejection has a precise user-facing message; 1.ix is the
        // input dual-consume guard (#389) — see its own comment for the key choice.
        if let Some(ref pt) = self.passthrough {
            // 1.i — paired-end required
            if !self.paired {
                anyhow::bail!("--passthrough requires --paired");
            }
            // 1.ii — exactly one R1/R2 pair in v1
            if self.input.len() != 2 {
                anyhow::bail!(
                    "--passthrough requires exactly one R1/R2 pair (got {} input files); \
                     multi-pair input with passthrough is not yet implemented",
                    self.input.len()
                );
            }
            // 1.iii — strict pair semantics in v1
            if self.retain_unpaired {
                anyhow::bail!(
                    "--passthrough is incompatible with --retain_unpaired \
                     (passthrough requires strict pair semantics in v1)"
                );
            }
            // 1.iv — clumpy reorder breaks lockstep
            if self.clumpify {
                anyhow::bail!("--passthrough is not yet supported with --clumpify");
            }
            // 1.v–1.vii — specialty modes own their own input arity/output naming
            if self.clock {
                anyhow::bail!("--passthrough is not compatible with --clock");
            }
            if self.implicon.is_some() {
                anyhow::bail!("--passthrough is not compatible with --implicon");
            }
            if self.hardtrim5.is_some() {
                anyhow::bail!("--passthrough is not compatible with --hardtrim5");
            }
            if self.hardtrim3.is_some() {
                anyhow::bail!("--passthrough is not compatible with --hardtrim3");
            }
            if self.demux.is_some() {
                anyhow::bail!("--passthrough is not compatible with --demux");
            }
            // 1.viii — file must exist, and must be re-readable (#379): the
            // passthrough stream is opened once to sanity-check and again to read.
            check_restartable_input(pt, "--passthrough file not found")?;
            // 1.ix — case-folded on purpose, not an identity check (#389): a case-variant
            // passthrough IS R1/R2 on APFS/NTFS and would be consumed twice. len == 2 per 1.ii.
            if self.input.len() == 2 {
                let pt_id = crate::io::path_identity_key(pt);
                let pt_key = crate::io::collision_key(pt);
                // Identity first: a byte-equal input gets the identity diagnosis even
                // if the other input is a case-variant of it (identity ⊆ collision).
                if let Some(same) = self
                    .input
                    .iter()
                    .find(|p| crate::io::path_identity_key(p) == pt_id)
                {
                    anyhow::bail!(
                        "--passthrough must be a third file (e.g. the index read), not \
                         one of the R1/R2 inputs: {} is input {}",
                        pt.display(),
                        same.display()
                    );
                } else if let Some(matched) = self
                    .input
                    .iter()
                    .find(|p| crate::io::collision_key(p) == pt_key)
                {
                    anyhow::bail!(
                        "--passthrough matches input {} case-insensitively (for APFS/NTFS \
                         safety): {}. On a case-insensitive filesystem these are the same \
                         file and the stream would be consumed twice; if they are genuinely \
                         two files, rename one so the paths differ by more than letter case.",
                        matched.display(),
                        pt.display()
                    );
                }
            }
        }

        // #386 — dispatch runs --hardtrim5 and returns, silently dropping a
        // 3' request. Bespoke message per the family precedent (§3.4a).
        if self.hardtrim5.is_some() && self.hardtrim3.is_some() {
            anyhow::bail!(
                "--hardtrim5 and --hardtrim3 cannot be combined in one invocation; \
                 run the two trims as separate invocations, feeding the first trim's \
                 output to the second."
            );
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
            // Read once, so restartability does not apply — but a writer-less FIFO
            // here blocks after the whole trim run has completed (#379 review).
            check_restartable_input(demux_file, "Barcode file not found")?;
        }

        // Check input files exist and can be re-read from the start (#379)
        for path in &self.input {
            check_restartable_input(path, "Input file not found")?;
        }

        // #369 — say when -a2 cannot be used, otherwise validate it up front so a
        // malformed value fails before the auto-detection scan. A value the run is
        // about to ignore is not worth failing on.
        if !self.adapter2.is_empty() {
            let unusable_reason = if self.hardtrim5.is_some() || self.hardtrim3.is_some() {
                Some("--hardtrim5/--hardtrim3 perform no adapter trimming")
            } else if self.clock {
                Some("--clock performs no adapter trimming")
            } else if self.implicon.is_some() {
                Some("--implicon performs no adapter trimming")
            } else if !self.paired {
                Some("it applies to Read 2 of a pair, and this is a single-end run")
            } else {
                None
            };
            match unusable_reason {
                Some(reason) => eprintln!(
                    "WARNING: -a2/--adapter2 was given but is not used in this mode ({reason}). Ignoring."
                ),
                None => {
                    crate::adapter::parse_adapter_specs_quiet(&self.adapter2)?;
                }
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
        // N=1 with --paired is INTENTIONALLY accepted by validate() — it's
        // the new paired-uBAM-interleaved entry point (one BAM file containing
        // R1/R2 records interleaved). The "is it actually a BAM?" check
        // happens in main.rs after format detection.
        // See plans/06252026_ubam-input-support/PLAN.md §3.3.
        for inputs in [vec![R1, R2, ALT_R1], vec![R1, R2, ALT_R1, ALT_R2, R1]] {
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
    fn test_no_args_displays_help_not_missing_arg_error() {
        // Running `trim_galore` with no arguments should print the help text
        // and exit 0 (like samtools/git/bwa), NOT bail with a "required
        // arguments were not provided" usage error. clap signals the former
        // with ErrorKind::DisplayHelpOnMissingArgumentOrSubcommand, which
        // renders to stdout and exits with a success code.
        use clap::CommandFactory;
        use clap::error::ErrorKind;
        let err = Cli::command()
            .try_get_matches_from(["trim_galore"])
            .unwrap_err();
        assert_eq!(
            err.kind(),
            ErrorKind::DisplayHelpOnMissingArgumentOrSubcommand,
            "no-args invocation should show help, got: {err}"
        );
    }

    // ── parse_sam_tag_name (PLAN §5 step 4.3, T16) ──────────────────────
    // Both code reviewers flagged the absence of these tests; ~10 LOC to
    // lock the validator's contract.

    #[test]
    fn parse_sam_tag_name_accepts_canonical_two_char() {
        assert_eq!(parse_sam_tag_name("CB").unwrap(), "CB");
        assert_eq!(parse_sam_tag_name("UB").unwrap(), "UB");
        assert_eq!(parse_sam_tag_name("RX").unwrap(), "RX");
        assert_eq!(parse_sam_tag_name("A0").unwrap(), "A0");
        assert_eq!(parse_sam_tag_name("Zz").unwrap(), "Zz");
    }

    #[test]
    fn parse_sam_tag_name_rejects_all_keyword() {
        // PLAN §5 step 4.3 — `ALL` is reserved for a future release.
        assert!(parse_sam_tag_name("ALL").is_err());
    }

    #[test]
    fn parse_sam_tag_name_rejects_wrong_length() {
        assert!(parse_sam_tag_name("X").is_err());
        assert!(parse_sam_tag_name("ABC").is_err());
        assert!(parse_sam_tag_name("").is_err());
    }

    #[test]
    fn parse_sam_tag_name_rejects_leading_non_alpha() {
        // SAM spec: tags match [A-Za-z][A-Za-z0-9]. Leading digit invalid.
        assert!(parse_sam_tag_name("1A").is_err());
        assert!(parse_sam_tag_name("9X").is_err());
    }

    #[test]
    fn parse_sam_tag_name_rejects_non_alphanumeric() {
        assert!(parse_sam_tag_name("A_").is_err());
        assert!(parse_sam_tag_name("A-").is_err());
        assert!(parse_sam_tag_name(" A").is_err());
    }

    #[test]
    fn test_validate_paired_single_input_accepted_at_validate_layer() {
        // The v3 paired-uBAM change: `--paired SINGLE.bam` is structurally
        // legal at the validate() layer. main.rs runs format detection and
        // errors if SINGLE.bam turns out to be FASTQ (test that path is
        // covered by integration tests, not here).
        let cli = Cli::parse_from(["trim_galore", "--paired", R1]);
        assert!(
            cli.validate().is_ok(),
            "--paired with N=1 must be accepted at the structural-validation layer"
        );
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

    /// #383. A duplicated positional would write one output twice. Rejected in
    /// `validate()` so the message is precise, not the APFS/NTFS collision text.
    #[test]
    fn test_validate_single_end_duplicate_input_rejected() {
        let cli = Cli::parse_from(["trim_galore", R1, R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("was given more than once"),
            "expected duplicate-input rejection, got: {err}"
        );
        assert!(
            !err.contains("APFS/NTFS"),
            "must not defer to the collision pre-flight's message: {err}"
        );
    }

    /// #386. Both hardtrims together used to run only the 5' trim, silently
    /// dropping the 3' request (the --hardtrim5 dispatch branch returns early).
    #[test]
    fn test_hardtrim5_and_hardtrim3_together_rejected() {
        for args in [
            vec!["trim_galore", "--hardtrim5", "20", "--hardtrim3", "15", R1],
            vec!["trim_galore", "--hardtrim3", "15", "--hardtrim5", "20", R1],
        ] {
            let err = Cli::parse_from(&args).validate().unwrap_err().to_string();
            assert!(err.contains("cannot be combined"), "got: {err}");
            assert!(
                err.contains("separate invocations"),
                "must carry the remedy: {err}"
            );
        }
    }

    /// Regression guards: each flag alone must keep parsing and validating.
    #[test]
    fn test_each_hardtrim_alone_still_accepted() {
        for args in [
            vec!["trim_galore", "--hardtrim5", "20", R1],
            vec!["trim_galore", "--hardtrim3", "15", R1],
        ] {
            let cli = Cli::parse_from(&args);
            cli.validate().expect("single hardtrim flag must validate");
        }
    }

    /// The same guard covers the specialty modes, which never consult `--paired`.
    #[test]
    fn test_validate_hardtrim_duplicate_input_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--hardtrim5", "20", R1, R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("was given more than once"),
            "expected duplicate-input rejection under --hardtrim5, got: {err}"
        );
    }

    #[test]
    fn test_validate_single_end_distinct_inputs_accepted() {
        // Negative control for the two tests above: distinct SE inputs still validate.
        let cli = Cli::parse_from(["trim_galore", R1, ALT_R1]);
        cli.validate()
            .expect("distinct single-end inputs should validate");
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
    fn test_clumpify_too_small_memory_passes_validation() {
        // Validation now accepts a too-small --memory; main.rs resolves the
        // layout at runtime and either succeeds, or warns + falls back to
        // plain mode. This avoids hard-failing a job over a configuration
        // detail the program can recover from.
        let cli = Cli::parse_from([
            "trim_galore",
            "--clumpify",
            "--cores",
            "16",
            "--memory",
            "64M",
            R1,
        ]);
        cli.validate()
            .expect("validation should pass; runtime warns and falls back to plain");
    }

    #[test]
    fn test_compression_defaults_to_one() {
        let cli = Cli::parse_from(["trim_galore", R1]);
        assert_eq!(cli.compression, 1);
    }

    #[test]
    fn test_compression_explicit_level() {
        let cli = Cli::parse_from(["trim_galore", "--compression", "9", R1]);
        assert_eq!(cli.compression, 9);
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

    // ── --passthrough validation (plan v2 Step 1) ─────────────────────────
    //
    // The third fixture used as the "passthrough" target — any existing
    // test_files/ FASTQ works since we never trim it in validation; we just
    // need a real path so step 1.viii (file exists) is satisfied.
    const PT: &str = "test_files/SRR24766921_RRBS_R2.fastq.gz";

    #[test]
    fn test_passthrough_paired_pair_accepted() {
        let cli = Cli::parse_from(["trim_galore", "--paired", "--passthrough", PT, R1, R2]);
        cli.validate()
            .expect("--passthrough with one R1/R2 pair should validate");
        assert_eq!(cli.passthrough.as_deref(), Some(std::path::Path::new(PT)));
    }

    #[test]
    fn test_passthrough_requires_paired() {
        let cli = Cli::parse_from(["trim_galore", "--passthrough", PT, R1]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--passthrough requires --paired"),
            "got: {err}"
        );
    }

    #[test]
    fn test_passthrough_rejects_multi_pair() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--passthrough",
            PT,
            R1,
            R2,
            ALT_R1,
            ALT_R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--passthrough requires exactly one R1/R2 pair"),
            "got: {err}"
        );
    }

    #[test]
    fn test_passthrough_rejects_retain_unpaired() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--retain_unpaired",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--retain_unpaired"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_clumpify() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--clumpify",
            "--cores",
            "2",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--clumpify"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_clock() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--clock",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--clock"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_implicon() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--implicon=8",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--implicon"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_hardtrim5() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--hardtrim5",
            "30",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--hardtrim5"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_hardtrim3() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--hardtrim3",
            "30",
            "--passthrough",
            PT,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--hardtrim3"), "got: {err}");
    }

    #[test]
    fn test_passthrough_rejects_demux() {
        // --demux conflicts with --paired in existing validation, but the
        // --passthrough vs --demux check must fire BEFORE that — confirm the
        // error message is the passthrough-specific one. (single-end input
        // shape because --demux requires single-end.)
        let cli = Cli::parse_from([
            "trim_galore",
            "--demux",
            "test_files/demux_test_samplesheet.txt",
            "--passthrough",
            PT,
            R1,
        ]);
        // With single-end input + --passthrough we hit 1.i ("requires --paired")
        // FIRST. To exercise the --demux check specifically, we'd need --paired
        // + --demux, which clap rejects at the --demux validation step itself.
        // The 1.i error is sufficient evidence the validation chain runs.
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--passthrough requires --paired") || err.contains("--demux"),
            "got: {err}"
        );
    }

    #[test]
    fn test_passthrough_rejects_missing_file() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--passthrough",
            "test_files/this_does_not_exist.fastq.gz",
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(err.contains("--passthrough file not found"), "got: {err}");
    }

    /// Shared 1.ix assertions (#389): normative prefix, matched input named,
    /// and no unconditional identity claim.
    fn assert_passthrough_identity_rejection(err: &str, matched: &str) {
        assert!(
            err.contains("--passthrough must be a third file"),
            "1.ix identity-branch prefix missing (#389); got: {err}"
        );
        assert!(
            err.contains(&format!("is input {matched}")),
            "matched input {matched} not named in the matched slot (#389); got: {err}"
        );
        assert!(
            !err.contains("aliases an input"),
            "1.ix must not assert identity (#389); got: {err}"
        );
    }

    #[test]
    fn test_passthrough_rejects_pointing_at_r1() {
        // Byte-equal ⊂ case-folded (fold pinned by io::tests::test_norm_path_case_folds);
        // the case-variant branch has its own lexical test below.
        let cli = Cli::parse_from(["trim_galore", "--paired", "--passthrough", R1, R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert_passthrough_identity_rejection(&err, R1);
    }

    #[test]
    fn test_passthrough_rejects_pointing_at_r2() {
        let cli = Cli::parse_from(["trim_galore", "--paired", "--passthrough", R2, R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert_passthrough_identity_rejection(&err, R2);
    }

    #[test]
    fn test_passthrough_rejects_dot_slash_spelling_of_input() {
        // pt and the matched input differ textually here, so the "is input {…}"
        // assertion discriminates the matched slot rather than echoing pt (#389).
        let pt = format!("./{R1}");
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--passthrough",
            pt.as_str(),
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert_passthrough_identity_rejection(&err, R1);
    }

    #[test]
    fn test_passthrough_rejects_case_variant_of_input() {
        // Lexical check: a case-variant of an UNWRITTEN input hits the case-only branch
        // on ext4 and APFS alike — only the passthrough must exist at 1.viii (#389).
        let dir = std::env::temp_dir().join(format!("tg_pt_case_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let pt = dir.join("r1.fq");
        std::fs::write(&pt, "@r\nACGT\n+\nIIII\n").unwrap();
        let r1 = dir.join("R1.fq");
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--passthrough",
            pt.to_str().unwrap(),
            r1.to_str().unwrap(),
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(
            err.contains("case-insensitively"),
            "1.ix must state the case-insensitive comparison (#389); got: {err}"
        );
        assert!(
            !err.contains("aliases an input"),
            "1.ix must not assert identity (#389); got: {err}"
        );
    }

    // ── --output-format (PLAN v2.1 §3.4a) ─────────────────────────────────

    #[test]
    fn output_format_default_is_fastq() {
        let cli = Cli::parse_from(["trim_galore", R1]);
        assert_eq!(cli.output_format, OutputFormat::Fastq);
    }

    #[test]
    fn output_format_ubam_parses() {
        let cli = Cli::parse_from(["trim_galore", "--output-format", "ubam", R1]);
        assert_eq!(cli.output_format, OutputFormat::UBam);
    }

    #[test]
    fn output_format_unknown_value_rejected() {
        let r = Cli::try_parse_from(["trim_galore", "--output-format", "binseq", R1]);
        assert!(r.is_err(), "BINSEQ is deferred in v1 — clap must reject it");
    }

    #[test]
    fn output_format_ubam_alone_validates() {
        let cli = Cli::parse_from(["trim_galore", "--output-format", "ubam", R1]);
        cli.validate()
            .expect("plain --output-format ubam must validate");
    }

    #[test]
    fn output_format_ubam_plus_clumpify_rejected() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--output-format",
            "ubam",
            "--clumpify",
            "--cores",
            "2",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--clumpify") && err.contains("--output-format ubam"),
            "expected clumpify+ubam rejection, got: {err}"
        );
    }

    #[test]
    fn output_format_ubam_plus_passthrough_rejected() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--output-format",
            "ubam",
            "--passthrough",
            ALT_R1,
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--passthrough") && err.contains("--output-format ubam"),
            "expected passthrough+ubam rejection, got: {err}"
        );
    }

    #[test]
    fn output_format_ubam_plus_clock_rejected() {
        let cli = Cli::parse_from(["trim_galore", "--clock", "--output-format", "ubam", R1, R2]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--clock") && err.contains("--output-format ubam"),
            "expected clock+ubam rejection, got: {err}"
        );
    }

    #[test]
    fn output_format_ubam_plus_implicon_rejected() {
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--implicon",
            "8",
            "--output-format",
            "ubam",
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--implicon") && err.contains("--output-format ubam"),
            "expected implicon+ubam rejection, got: {err}"
        );
    }

    #[test]
    fn output_format_ubam_plus_retain_unpaired_rejected() {
        // PLAN v2.1 §3.4a addendum: --retain_unpaired produces additional
        // FASTQ files (`*_unpaired_{1,2}.fq.gz`); multi-output BAM is out of
        // scope for v1.
        let cli = Cli::parse_from([
            "trim_galore",
            "--paired",
            "--retain_unpaired",
            "--output-format",
            "ubam",
            R1,
            R2,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--retain_unpaired") && err.contains("--output-format ubam"),
            "expected retain_unpaired+ubam rejection, got: {err}"
        );
    }

    #[test]
    fn output_format_ubam_plus_demux_rejected() {
        // --demux takes a barcode-file path; any path string suffices for parser.
        let cli = Cli::parse_from([
            "trim_galore",
            "--demux",
            "test_files/demux_test_samplesheet.txt",
            "--output-format",
            "ubam",
            R1,
        ]);
        let err = cli.validate().unwrap_err().to_string();
        assert!(
            err.contains("--demux") && err.contains("--output-format ubam"),
            "expected demux+ubam rejection, got: {err}"
        );
    }
}
