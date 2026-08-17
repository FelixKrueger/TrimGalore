use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use trim_galore::adapter;
use trim_galore::bam::BamReader;
use trim_galore::cli::{Cli, rewrite_perl_short_flags};
use trim_galore::clump;
use trim_galore::clump_only;
use trim_galore::demux;
use trim_galore::fastq::{FastqReader, FastqWriter, RecordSource};
use trim_galore::fastqc;
use trim_galore::filters::{MaxNFilter, UnpairedLengths};
use trim_galore::format::{
    InputFormat, PairedShape, detect_input_format, open_sync_reader, open_threaded_reader,
    reject_bam_format_mismatch_in_pair,
};
use trim_galore::io as naming;
use trim_galore::parallel;
use trim_galore::report;
use trim_galore::specialty;
use trim_galore::trimmer;

/// mimalloc replaces the platform allocator for the whole binary. Allocation
/// behaviour is the only thing that changes: no output byte differs, so the
/// Perl 0.6.11 validation matrix is unaffected.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Format-aware sanity check — dispatches `FastqReader::sanity_check` for
/// FASTQ input or peek-reads the first BAM record (asserting `is_unmapped()`)
/// for uBAM. The per-record aligned-BAM check in `BamReader::next_record`
/// catches mixed-aligned BAMs that slip past this fast-path.
fn sanity_check_any(path: &std::path::Path) -> Result<()> {
    match detect_input_format(path)? {
        // Hand the content-detected verdict down rather than letting
        // `sanity_check` re-derive it from the filename, which gets
        // `.fq.bgz` wrong.
        fmt @ (InputFormat::FastqPlain | InputFormat::FastqGz) => {
            FastqReader::sanity_check_with(path, fmt == InputFormat::FastqGz)
        }
        InputFormat::UnalignedBam => {
            let mut r = BamReader::open(path)?;
            match r.next_record()? {
                None => anyhow::bail!(
                    "Input file '{}' is a uBAM with no records (empty BAM).",
                    path.display()
                ),
                Some(_) => Ok(()),
            }
        }
    }
}

type AdapterList = Vec<(String, String)>;
type SetupResult = Result<(String, AdapterList, AdapterList, trimmer::TrimConfig)>;
type ResolvedAdapter = Result<(String, AdapterList, AdapterList, Option<(usize, usize)>)>;

/// Hint for the sites whose report names carry no positional discriminator:
/// paired trim (FASTQ and uBAM out) and `--clump_only --paired` (#391).
const PAIRED_REPORT_HINT: &str = "Every output of a pair — validated reads, reports, \
                                  and the --passthrough carrier — is written to \
                                  --output_dir, or to Read 1's directory when that is \
                                  not given. Names come from the input filename alone, \
                                  so inputs sharing a filename collide wherever they \
                                  live: rename one input, give the pairs distinct \
                                  --output_dir runs, or pass --no_report_file if only \
                                  the reports collide.";

/// CWD-output modes' remediation; see `PAIRED_REPORT_HINT` for the paired sites.
const CWD_OUTPUT_HINT: &str = "This mode writes output to the current working directory, \
                               so inputs sharing a basename collide whatever `--output_dir` \
                               is set to — run one invocation per input, or give the inputs \
                               distinct basenames.";

/// Every file a run reads and must therefore never write over: the positionals,
/// `--passthrough`, and `--demux`'s barcode file.
fn guarded_inputs(cli: &Cli) -> Vec<std::path::PathBuf> {
    let mut v = cli.input.clone();
    if let Some(ref pt) = cli.passthrough {
        v.push(pt.clone());
    }
    if let Some(ref bc) = cli.demux {
        v.push(bc.clone());
    }
    v
}

/// Secondary outputs an SE FASTQ trim run writes besides the trimmed file: the two
/// trimming reports, and — when `--demux` is set — every per-barcode file.
///
/// Needed because the pre-flight's output-vs-input check is only as complete as the
/// path list it is given: assumption A2 argues that checking *primary* paths covers
/// secondary paths for output-vs-**output** collisions, and says nothing about
/// output-vs-**input**. A secondary output can equal a named input while every
/// primary stays distinct.
fn planned_secondary_outputs(
    cli: &Cli,
    output_dir: Option<&Path>,
    gzip: bool,
) -> Result<Vec<naming::PlannedOutput>> {
    let mut v = Vec::new();
    for input in &cli.input {
        let src = || naming::OutputSource::Input(input.clone());
        if !cli.no_report_file {
            v.push((naming::report_name(input, output_dir), src()));
            v.push((naming::json_report_name(input, output_dir), src()));
        }
        if let Some(ref barcode_file) = cli.demux {
            let barcodes = demux::read_barcode_file(barcode_file)?;
            let trimmed =
                naming::single_end_output_name(input, output_dir, cli.basename.as_deref(), gzip);
            v.extend(
                demux::demux_output_paths(&trimmed, &barcodes, gzip, output_dir)
                    .into_iter()
                    .map(|p| (p, src())),
            );
        }
    }
    Ok(v)
}

/// Clumping-report paths a `--clump_only` run plans — one per report-keyed input,
/// nothing under `--no_report_file` (matching every writer's gate in clump_only.rs).
/// The pre-flight is only as complete as its candidate list: #391 was the paired
/// arm planning primaries alone while the writer also wrote per-mate reports.
fn clump_report_candidates(
    no_report_file: bool,
    report_inputs: &[impl AsRef<Path>],
    output_dir: Option<&Path>,
) -> Vec<naming::PlannedOutput> {
    if no_report_file {
        return Vec::new();
    }
    report_inputs
        .iter()
        .map(|input| {
            let input = input.as_ref();
            (
                naming::clumping_report_name(input, output_dir),
                naming::OutputSource::Input(input.to_path_buf()),
            )
        })
        .collect()
}

/// Prospective hardtrim output paths for every input, per resolved output format.
fn planned_hardtrim_outputs(
    cli: &Cli,
    keep: usize,
    end: specialty::HardtrimEnd,
    output_dir: Option<&Path>,
    gzip: bool,
) -> Vec<naming::PlannedOutput> {
    cli.input
        .iter()
        .map(|input| {
            let path = match cli.output_format {
                trim_galore::cli::OutputFormat::Fastq => {
                    specialty::hardtrim_output_name(input, keep, end, output_dir, gzip)
                }
                trim_galore::cli::OutputFormat::UBam => {
                    specialty::hardtrim_bam_output_name(input, keep, end, output_dir)
                }
            };
            (path, naming::OutputSource::Input(input.clone()))
        })
        .collect()
}

/// Resolve `--memory` into a `(n_bins, bin_byte_budget)` layout when
/// `--clumpify` is set, and emit a one-line startup notice with the
/// resolved values. Returns `None` if clumpify is off — or if the budget
/// is below the floor, in which case we print a loud warning and fall back
/// to plain mode rather than refusing to run.
fn resolve_clump_layout(cli: &Cli) -> Result<Option<clump::ClumpLayout>> {
    if !cli.clumpify {
        return Ok(None);
    }
    let memory_bytes =
        clump::parse_memory_size(&cli.memory).map_err(|e| anyhow::anyhow!("--memory: {e}"))?;
    let layout_inputs = clump::LayoutInputs::uniform(
        cli.cores,
        cli.fastqc_requested()
            .then(|| fastqc::resolved_threads(cli.fastqc_args.as_deref(), cli.cores)),
    );
    let min_required = clump::clumpify_min_memory_bytes(&layout_inputs);
    if memory_bytes < min_required {
        eprintln!();
        eprintln!(
            "WARNING: --memory {} is too small for --clumpify at --cores {} \
             (need ≥ {} MiB).",
            cli.memory,
            cli.cores,
            min_required.div_ceil(1024 * 1024),
        );
        eprintln!(
            "         Falling back to plain mode (no read reordering). \
             Increase --memory or drop --clumpify to silence this warning."
        );
        eprintln!();
        return Ok(None);
    }
    let layout = clump::resolve_layout(memory_bytes, &layout_inputs)?;
    eprintln!(
        "clumpify: {} bins × {} MiB; predicted peak ≈ {} MiB (gzip level {})",
        layout.n_bins,
        layout.bin_byte_budget / (1024 * 1024),
        layout.predicted_peak_bytes() / (1024 * 1024),
        cli.compression,
    );
    Ok(Some(layout))
}

fn main() -> Result<()> {
    env_logger::init();

    // Pre-parse rewrite for Perl-era `-r1`/`-r2` short flags before clap sees them.
    // A bare `trim_galore` (no args) prints the full help to stdout and exits 0,
    // like most modern CLIs, rather than a non-zero "missing arguments" error.
    // clap surfaces the help/version outcomes as `Err`, so intercept those three
    // kinds and route them to stdout with a success code; genuine usage errors
    // keep clap's default stderr + non-zero behaviour.
    let cli = match Cli::try_parse_from(rewrite_perl_short_flags(std::env::args())) {
        Ok(cli) => cli,
        Err(err) => {
            use clap::error::ErrorKind::{
                DisplayHelp, DisplayHelpOnMissingArgumentOrSubcommand, DisplayVersion,
            };
            if matches!(
                err.kind(),
                DisplayHelp | DisplayHelpOnMissingArgumentOrSubcommand | DisplayVersion
            ) {
                print!("{err}");
                std::process::exit(0);
            }
            err.exit();
        }
    };
    cli.validate()?;

    // Command-line as seen by the user (pre-clap-rewrite). Used by the
    // uBAM-output path to populate `@PG CL:`; consumed by the FASTQ path's
    // report writers via their own inline `std::env::args().collect(...)`.
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");

    // Input sanity check on first file
    eprintln!("\nTrim Galore v{}", env!("CARGO_PKG_VERSION"));
    eprintln!("{}", env!("VERSION_BODY"));
    eprintln!("==================================================\n");

    sanity_check_any(&cli.input[0])?;

    // ─── uBAM-specific format-dependent validation (T17) ───────────────────
    // Detect format of all inputs once; the call is cheap (24-byte peek +
    // at most one BGZF block decompress). Used for both validation and
    // dispatch below.
    let input_formats: Vec<InputFormat> = cli
        .input
        .iter()
        .map(|p| detect_input_format(p))
        .collect::<Result<_>>()?;
    let any_bam = input_formats
        .iter()
        .any(|f| matches!(f, InputFormat::UnalignedBam));

    // Issue #358 — `--phred64` declares the ASCII offset of the *input*. BAM
    // has no ASCII encoding to declare: it stores raw Phred (0–93) by spec,
    // and `BamReader` always yields Phred+33 ASCII internally. So the flag is
    // meaningless for BAM input in *every* mode, and actively harmful in two
    // ways: with quality trimming it subtracts 64 from Phred+33 data and
    // discards every read as low-quality; under `--clump_only` (which does no
    // quality arithmetic) it would silently zero the output quality once the
    // writer honours the offset.
    //
    // Rejected uniformly rather than per-mode. `--hardtrim5/3` and
    // `--clump_only` currently accept the flag as an inert no-op, so this does
    // remove working invocations — but a mode-dependent rule ("rejected unless
    // you're in hardtrim, or clump-only-to-FASTQ, …") is worse to document,
    // worse to test, and one refactor away from being wrong.
    //
    // Sited here, not in `Cli::validate()` §3.4a: validate() runs before input
    // format detection and cannot see `any_bam`. Same reason §3.4b below lives
    // in main.rs. Placed ahead of every dispatch branch — including the
    // early-returning specialty and clump-only paths — and ahead of
    // `ensure_output_dir`, so a rejected run creates nothing.
    if cli.phred64 && any_bam {
        anyhow::bail!(
            "--phred64 cannot be used with unaligned BAM input. BAM stores raw \
             Phred scores directly, so there is no ASCII encoding to declare — \
             the reader always yields Phred+33 internally. Passing --phred64 \
             here subtracts 64 from Phred+33 data, reducing every score by 31 \
             — which discards effectively the whole library as low-quality, or, \
             in modes that do no quality trimming such as --clump_only, silently \
             degrades the output quality instead. Drop --phred64 for BAM input."
        );
    }

    // INVARIANT ESTABLISHED HERE, RELIED ON FAR BELOW: past this point,
    // `cli.phred_offset()` is guaranteed to be 33 whenever any input is uBAM.
    // Every `cli.phred_offset()` call that feeds `BamWriter::create` depends on
    // it — the trim drivers (`run_ubam_output_*`), the specialty hardtrim
    // drivers, and the `--clump_only` BAM writers. Without this guard, a uBAM
    // input under `--phred64` would have 64 subtracted from Phred+33 data,
    // reducing every score by 31: read loss when quality trimming is active,
    // and a silent breach of `--clump_only`'s documented lossless-qual
    // invariant when it is not. Do not relax this to a warning.

    if !cli.preserve_tags.is_empty() && !any_bam {
        // PLAN v2.1 §3.4b — format-detection-time rule. With --output-format
        // ubam, --preserve-tags + all-FASTQ-inputs is a hard error because
        // there are no source tags to propagate AND the user explicitly
        // requested uBAM output (so the otherwise-silent loss has no
        // pressure-valve). Without --output-format ubam, it's just a
        // user-info warning (existing behaviour).
        if matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam) {
            anyhow::bail!(
                "--preserve-tags has no effect with all-FASTQ inputs; either remove \
                 the flag or convert at least one input to uBAM via 'samtools import'"
            );
        }
        eprintln!(
            "WARNING: --preserve-tags has no effect — no uBAM input was detected. \
             The flag is ignored for FASTQ input."
        );
    }
    if cli.passthrough.is_some() && any_bam {
        anyhow::bail!(
            "--passthrough is not supported with uBAM input in this release. \
             Either convert the uBAM to FASTQ via `samtools fastq` first, or drop --passthrough."
        );
    }
    // `any_bam` is computed from `cli.input`; the passthrough file is a separate
    // surface, so it needs its own format check — before any output is created.
    if let Some(ref pt) = cli.passthrough
        && matches!(detect_input_format(pt)?, InputFormat::UnalignedBam)
    {
        anyhow::bail!(
            "--passthrough is not supported with uBAM input in this release. \
             Either convert the uBAM to FASTQ via `samtools fastq` first, or drop --passthrough."
        );
    }
    if cli.paired && cli.input.len() == 1 && !matches!(input_formats[0], InputFormat::UnalignedBam)
    {
        anyhow::bail!(
            "--paired with a single input file is only legal if that file is a uBAM. \
             Provide R1 and R2 as separate files for FASTQ paired-end mode."
        );
    }

    // Issue #363 — both inputs of a `--paired` pair must agree on being BAM.
    //
    // Sited here for the same reason as the `--phred64` guard above: the
    // decision needs `input_formats`, which `Cli::validate()` cannot see. The
    // position matters beyond that, though — this is ahead of
    // `ensure_output_dir` (below), the three output-collision pre-flights, and
    // every dispatch branch. The three guards this replaces all sat *past* that
    // line: the trim-path one fired only after adapter auto-detection had
    // scanned up to 1M reads, the banner had printed, the output directory
    // existed, and on multi-pair input earlier pairs had already been written.
    //
    // Specialty modes are exempt: `--hardtrim5/3` never consult `cli.paired`
    // (each file is processed independently, so a mixed pair is genuinely
    // harmless there), and `--clock`/`--implicon` are left alone only because
    // changing them is out of scope for #363 — they pair for real and today
    // report a misleading read-count error on mixed input.
    //
    // N == 1 is excluded: that is single-file interleaved uBAM, validated
    // directly above.
    if cli.paired
        && cli.input.len() > 1
        && cli.hardtrim5.is_none()
        && cli.hardtrim3.is_none()
        && !cli.clock
        && cli.implicon.is_none()
    {
        let ubam_out = matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam);
        let shape = match (cli.clump_only, ubam_out) {
            (true, true) => PairedShape::ClumpOnlyUbamOut,
            (true, false) => PairedShape::ClumpOnlyFastqOut,
            (false, true) => PairedShape::TrimUbamOut,
            (false, false) => PairedShape::Trim,
        };
        reject_bam_format_mismatch_in_pair(&cli.input, &input_formats, shape)?;
    }

    // #408 — format-gated, so it cannot live in `Cli::validate()`; sited after the
    // two structural pair checks above so those report their own defect first. The
    // clip-flag list must match every `--rename`-driven `append_to_id` site
    // (`trimmer.rs` clip_5/clip_3, `specialty.rs` hardtrim) — without one of them
    // set, `--rename` appends nothing and there is nothing to lose.
    // The clip list reads through `effective_clips()` so a `--library` preset
    // (#440) counts here exactly as the four flags it stands in for would.
    if cli.rename
        && matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam)
        && (cli.effective_clips().any_set() || cli.hardtrim5.is_some() || cli.hardtrim3.is_some())
        && input_formats
            .iter()
            .any(|f| !matches!(f, InputFormat::UnalignedBam))
    {
        anyhow::bail!(
            "--rename is refused with --output-format ubam when any input is FASTQ. \
             Whether the annotation can be represented depends on the individual \
             header: the :clip5:/:clip3: suffix is appended to the end of the read ID, \
             so a header carrying text after the first space puts the annotation inside \
             that text, and BAM read names cannot contain whitespace, so none of that \
             tail reaches the output. The whole run is refused rather than decided per \
             record, because a per-record decision would keep the annotation for some \
             reads and silently drop it for others. Use \
             FASTQ output, or drop --rename. uBAM input is unaffected: BAM read names \
             carry no description, so there the annotation lands on the name itself."
        );
    }

    // Output gzip mode. Mirror Perl: by default the output's compression
    // matches the input's (plain `.fastq` → plain `.fq`, gzipped `.fastq.gz`
    // → gzipped `.fq.gz`). `--dont_gzip` overrides to always-plain. The
    // first input determines the mode for the whole run; mixing plain and
    // gzipped inputs in one invocation isn't a supported configuration.
    // See #245 for the parity rationale (Rust v2.1.0-beta.5 always gzipped
    // regardless of input, breaking pipelines that globbed `*.fq` no-gz).
    let gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0]);
    let output_dir = cli.output_dir.as_deref();

    // Auto-create --output_dir if it doesn't exist. See io::ensure_output_dir
    // for why this has to happen here and not lazily per-file.
    naming::ensure_output_dir(output_dir)?;

    // uBAM-output startup NOTEs — fired BEFORE dispatch so both specialty
    // (hardtrim) and normal-trim BAM paths see them. Code-review round-2
    // B-NIT-2 consolidation: previously these lived in run_ubam_output and
    // never fired on the hardtrim BAM paths (which return earlier).
    if matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam) {
        if cli.cores > 1 {
            eprintln!(
                "NOTE: --output-format ubam uses single-threaded compression in v1; \
                 --cores {} is ignored for BAM writing (FastQC, if requested, still \
                 uses it). For high-throughput uBAM output, run multiple invocations \
                 in parallel.",
                cli.cores
            );
        }
        if any_bam {
            // PLAN §3.3 step 4 — uBAM-in → uBAM-out cannot carry the BAM
            // missing-qual sentinel (0xFF) through the FASTQ intermediate.
            eprintln!(
                "NOTE: uBAM-in → uBAM-out: missing-qual sentinel (BAM 0xFF) does \
                 NOT round-trip — records with no original qual will emit Phred 0 \
                 ('!' × seq_len) on output instead of the 0xFF sentinel."
            );
        }
        if cli.phred64 {
            // Issue #358 O-3 — make the encoding assumption visible. Without
            // this, a Phred+33 file mistakenly run with --phred64 produces an
            // all-Phred-0 BAM that is indistinguishable downstream from
            // legitimate Q0 data. Also disambiguates the trimming report,
            // which prints "Quality encoding type selected: ASCII+64" (a
            // statement about the INPUT) next to a BAM holding raw scores.
            eprintln!(
                "NOTE: input is Phred+64 (ASCII+64); output BAM QUAL stores true \
                 Phred scores (0–93) per the SAM spec, not ASCII. Verify the input \
                 really is Phred+64 — running Phred+33 data with --phred64 reduces \
                 every score by 31 (floored at 0), turning typical Q2–Q41 into \
                 Q0–Q10. The result reads as plausible poor-quality data rather \
                 than obviously empty, so it is easy to miss downstream."
            );
        }
    }

    // Specialty modes — bypass normal trimming pipeline entirely
    if let Some(n) = cli.hardtrim5 {
        // #383 — collide across all inputs before the first write.
        naming::preflight_output_collisions(
            &planned_hardtrim_outputs(&cli, n, specialty::HardtrimEnd::Five, output_dir, gzip),
            &guarded_inputs(&cli),
            Some(CWD_OUTPUT_HINT),
        )?;
        for input in &cli.input {
            match cli.output_format {
                trim_galore::cli::OutputFormat::Fastq => specialty::hardtrim5(
                    input,
                    n,
                    gzip,
                    output_dir,
                    cli.rename,
                    cli.cores,
                    cli.compression,
                )?,
                trim_galore::cli::OutputFormat::UBam => specialty::hardtrim5_to_bam(
                    input,
                    n,
                    output_dir,
                    cli.rename,
                    &cli.preserve_tags,
                    &command_line,
                    cli.phred_offset(),
                )?,
            }
        }
        return Ok(());
    }
    if let Some(n) = cli.hardtrim3 {
        // #383 — collide across all inputs before the first write.
        naming::preflight_output_collisions(
            &planned_hardtrim_outputs(&cli, n, specialty::HardtrimEnd::Three, output_dir, gzip),
            &guarded_inputs(&cli),
            Some(CWD_OUTPUT_HINT),
        )?;
        for input in &cli.input {
            match cli.output_format {
                trim_galore::cli::OutputFormat::Fastq => specialty::hardtrim3(
                    input,
                    n,
                    gzip,
                    output_dir,
                    cli.rename,
                    cli.cores,
                    cli.compression,
                )?,
                trim_galore::cli::OutputFormat::UBam => specialty::hardtrim3_to_bam(
                    input,
                    n,
                    output_dir,
                    cli.rename,
                    &cli.preserve_tags,
                    &command_line,
                    cli.phred_offset(),
                )?,
            }
        }
        return Ok(());
    }
    if cli.clock {
        run_specialty_paired(
            &cli,
            "Clock",
            Some(CWD_OUTPUT_HINT),
            |r1, r2| {
                let src = naming::OutputSource::Pair(r1.to_path_buf(), r2.to_path_buf());
                vec![
                    (
                        specialty::clock_output_name(r1, "R1", output_dir, gzip),
                        src.clone(),
                    ),
                    (
                        specialty::clock_output_name(r2, "R2", output_dir, gzip),
                        src,
                    ),
                ]
            },
            |r1, r2| specialty::clock(r1, r2, gzip, output_dir, cli.cores, cli.compression),
        )?;
        return Ok(());
    }
    if let Some(umi_len) = cli.implicon {
        run_specialty_paired(
            &cli,
            "IMPLICON",
            Some(CWD_OUTPUT_HINT),
            |r1, r2| {
                let src = naming::OutputSource::Pair(r1.to_path_buf(), r2.to_path_buf());
                vec![
                    (
                        specialty::implicon_output_name(r1, umi_len, "R1", output_dir, gzip),
                        src.clone(),
                    ),
                    (
                        specialty::implicon_output_name(r2, umi_len, "R2", output_dir, gzip),
                        src,
                    ),
                ]
            },
            |r1, r2| {
                specialty::implicon(
                    r1,
                    r2,
                    umi_len,
                    gzip,
                    output_dir,
                    cli.cores,
                    cli.compression,
                )
            },
        )?;
        return Ok(());
    }
    if cli.clump_only {
        // --clump_only: lossless reorder-only specialty mode. Feature #353.
        // v1: FASTQ in/out. v2: uBAM in/out via --output-format ubam.
        let memory_bytes =
            clump::parse_memory_size(&cli.memory).map_err(|e| anyhow::anyhow!("--memory: {e}"))?;
        let basename = cli.basename.as_deref();
        match cli.output_format {
            trim_galore::cli::OutputFormat::Fastq => {
                // v1 FASTQ dispatch, unchanged except for one v2 guard.
                // Guard against B-C1 PANIC: --clump_only --paired with N=1
                // (single BAM) reaches `run_specialty_paired`, which does
                // `chunk[1]` on a length-1 chunk → index-out-of-bounds. The
                // general N=1 carve-out at the top of main.rs allows this
                // through because it exists for the trim uBAM path; here it
                // must be rejected because the FASTQ output arm can't handle
                // uBAM input (would drop aux tags) and no other N=1 shape
                // makes sense under --paired.
                if cli.paired && cli.input.len() == 1 {
                    anyhow::bail!(
                        "--clump_only --paired requires two FASTQ input files. \
                         Single-file interleaved uBAM input needs --output-format ubam \
                         (add `--output-format ubam` to the command line)."
                    );
                }
                if cli.paired {
                    run_specialty_paired(
                        &cli,
                        "--clump_only",
                        // #391 — reports carry no _clumped_N discriminator; the hint's
                        // --no_report_file remedy is the escape hatch when only they collide.
                        Some(PAIRED_REPORT_HINT),
                        |r1, r2| {
                            let (o1, o2) = naming::clumped_paired_output_names(
                                r1, r2, output_dir, basename, gzip,
                            );
                            let src =
                                naming::OutputSource::Pair(r1.to_path_buf(), r2.to_path_buf());
                            let mut v = vec![(o1, src.clone()), (o2, src)];
                            let pair_dir = naming::pair_output_dir(r1, output_dir);
                            v.extend(clump_report_candidates(
                                cli.no_report_file,
                                &[r1, r2],
                                Some(&pair_dir),
                            ));
                            v
                        },
                        |r1, r2| {
                            clump_only::clump_only_paired(
                                r1,
                                r2,
                                output_dir,
                                basename,
                                gzip,
                                cli.cores,
                                memory_bytes,
                                cli.compression,
                                cli.fastqc_requested(),
                                cli.fastqc_args.as_deref(),
                                cli.no_report_file,
                            )
                            .map(|_| ())
                        },
                    )?;
                } else {
                    // SE FASTQ pre-flight (issues #216, #383).
                    let mut planned: Vec<naming::PlannedOutput> = cli
                        .input
                        .iter()
                        .map(|input| {
                            (
                                naming::clumped_output_name(input, output_dir, basename, gzip),
                                naming::OutputSource::Input(input.clone()),
                            )
                        })
                        .collect();
                    // #391 — reports join so the input check covers them too.
                    planned.extend(clump_report_candidates(
                        cli.no_report_file,
                        &cli.input,
                        output_dir,
                    ));
                    naming::preflight_output_collisions(&planned, &guarded_inputs(&cli), None)?;
                    for input in &cli.input {
                        clump_only::clump_only_single(
                            input,
                            output_dir,
                            basename,
                            gzip,
                            cli.cores,
                            memory_bytes,
                            cli.compression,
                            cli.fastqc_requested(),
                            cli.fastqc_args.as_deref(),
                            cli.no_report_file,
                        )?;
                    }
                }
            }
            trim_galore::cli::OutputFormat::UBam => {
                // v2 uBAM dispatch. Format-guards for PE input shapes + explicit
                // multi-pair iteration + collision pre-flight before opening any
                // reader. Design details: PLAN_v2_ubam.md §Implementation outline
                // step 2 + §Behavior "PE steps".
                if cli.paired {
                    if cli.input.len() == 1 {
                        // Shape B: single interleaved uBAM. The "must be BAM"
                        // check that used to live here was already unreachable —
                        // its condition is a strict subset of the `--paired`
                        // N=1 non-BAM guard in main() — and was retired with
                        // #363. Retained note so the absence is intentional.
                        let mut planned = vec![(
                            naming::clumped_paired_bam_output_name(
                                &cli.input[0],
                                None,
                                output_dir,
                                basename,
                            ),
                            naming::OutputSource::Input(cli.input[0].clone()),
                        )];
                        // #391 — defensive symmetry: with one input the report (filename
                        // plus a suffix) can never alias it; the arm keeps its siblings' shape.
                        planned.extend(clump_report_candidates(
                            cli.no_report_file,
                            &cli.input,
                            output_dir,
                        ));
                        naming::preflight_output_collisions(&planned, &guarded_inputs(&cli), None)?;
                        clump_only::clump_only_paired_to_bam_one_pair(
                            &cli.input,
                            output_dir,
                            basename,
                            cli.cores,
                            memory_bytes,
                            &cli.preserve_tags,
                            &command_line,
                            cli.fastqc_requested(),
                            cli.fastqc_args.as_deref(),
                            cli.no_report_file,
                            cli.phred_offset(),
                        )?;
                    } else {
                        // Shape A: two-file paired (multi-pair supported for N=2, 4, 6, …).
                        // Two-BAM and mixed-format Shape A are rejected by
                        // `reject_bam_format_mismatch_in_pair` in main(), which
                        // runs before this dispatch (and before the pre-flight
                        // below). This branch previously carried its own copy of
                        // that check — the only one of the three that got the
                        // two-BAM-vs-mixed distinction right, and the model for
                        // the shared helper (#363).
                        // Multi-pair collision pre-flight (case-folded per issue #216).
                        let mut planned: Vec<naming::PlannedOutput> = Vec::new();
                        for chunk in cli.input.chunks(2) {
                            planned.push((
                                naming::clumped_paired_bam_output_name(
                                    &chunk[0],
                                    Some(&chunk[1]),
                                    output_dir,
                                    basename,
                                ),
                                naming::OutputSource::Pair(chunk[0].clone(), chunk[1].clone()),
                            ));
                            // #391 — ONE report per pair, keyed on the pair's first input
                            // (clump_only.rs writes clumping_report_name(inputs[0], …)).
                            // Keyed on R1 already, so `pair_output_dir` moves nothing here;
                            // routed through it so every paired site derives alike (#398).
                            let pair_dir = naming::pair_output_dir(&chunk[0], output_dir);
                            planned.extend(clump_report_candidates(
                                cli.no_report_file,
                                std::slice::from_ref(&chunk[0]),
                                Some(&pair_dir),
                            ));
                        }
                        naming::preflight_output_collisions(&planned, &guarded_inputs(&cli), None)?;
                        // Per-pair iteration. Mirrors the trim FASTQ multi-pair
                        // shape (main.rs:651-669) with pair-progress banner +
                        // per-pair sanity-checks + `.with_context()` error
                        // wrapping so a failure at pair 3/5 identifies the pair.
                        let total_pairs = cli.input.len() / 2;
                        for (pair_idx, chunk) in cli.input.chunks(2).enumerate() {
                            if total_pairs > 1 {
                                eprintln!("\n=== Pair {} of {} ===", pair_idx + 1, total_pairs);
                            }
                            // Sanity-check each input before opening any reader.
                            // R1 of pair 0 was already sanity-checked at main
                            // entry (main.rs:153); skip that. R2 of pair 0 and
                            // both files of every later pair need checking.
                            if pair_idx > 0 {
                                sanity_check_any(&chunk[0])?;
                            }
                            sanity_check_any(&chunk[1])?;
                            clump_only::clump_only_paired_to_bam_one_pair(
                                chunk,
                                output_dir,
                                basename,
                                cli.cores,
                                memory_bytes,
                                &cli.preserve_tags,
                                &command_line,
                                cli.fastqc_requested(),
                                cli.fastqc_args.as_deref(),
                                cli.no_report_file,
                                cli.phred_offset(),
                            )
                            .with_context(|| {
                                format!(
                                    "pair {} of {} ({} + {})",
                                    pair_idx + 1,
                                    total_pairs,
                                    chunk[0].display(),
                                    chunk[1].display()
                                )
                            })?;
                        }
                    }
                } else {
                    // SE BAM. Case-folded collision pre-flight across inputs.
                    let mut planned: Vec<naming::PlannedOutput> = Vec::new();
                    for input in &cli.input {
                        planned.push((
                            naming::clumped_bam_output_name(input, output_dir, basename),
                            naming::OutputSource::Input(input.clone()),
                        ));
                    }
                    // #391 — reports join so the input check covers them too.
                    planned.extend(clump_report_candidates(
                        cli.no_report_file,
                        &cli.input,
                        output_dir,
                    ));
                    naming::preflight_output_collisions(&planned, &guarded_inputs(&cli), None)?;
                    for input in &cli.input {
                        clump_only::clump_only_single_to_bam(
                            input,
                            output_dir,
                            basename,
                            cli.cores,
                            memory_bytes,
                            &cli.preserve_tags,
                            &command_line,
                            cli.fastqc_requested(),
                            cli.fastqc_args.as_deref(),
                            cli.no_report_file,
                            cli.phred_offset(),
                        )?;
                    }
                }
            }
        }
        return Ok(());
    }

    if cli.discard_untrimmed {
        eprintln!("Discarding reads without adapter match (--discard-untrimmed)");
    }

    if cli.cores > 1 {
        eprintln!(
            "Using {} worker threads (parallel trim + compress)",
            cli.cores
        );
    }

    // ─── Output format dispatch (PLAN v2.1 §5 step 4) ──────────────────────
    // uBAM output lives on a separate serial code path; the FASTQ path
    // below is the existing parallel-or-serial dispatch, untouched.
    if matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam) {
        return run_ubam_output(&cli, output_dir, &command_line);
    }

    // Paired-uBAM single-file de-interleaved path (T20-PE). The validation
    // above guarantees `cli.input[0]` is BAM when `--paired` and len == 1.
    // Routes through the single-file paired-uBAM helper which sets up its
    // own output paths and calls `BamReader::open_paired_interleaved`.
    if cli.paired && cli.input.len() == 1 {
        let outputs = naming::InterleavedFastqOutputs::new(
            &cli.input[0],
            output_dir,
            gzip,
            cli.retain_unpaired,
            !cli.no_report_file,
        );
        // `Input`, not `Pair`: one input means `Pair` could only be `Pair(x, x)`,
        // which renders `one.bam + one.bam`.
        let src = naming::OutputSource::Input(cli.input[0].clone());
        // Ahead of setup_trimming, which scans up to 1 M reads for adapters — a
        // refusal must not arrive after that wait.
        naming::preflight_output_collisions(
            &outputs.planned(&src),
            &guarded_inputs(&cli),
            Some(PAIRED_REPORT_HINT),
        )?;
        let (_label, adapters_r1, adapters_r2, config) = setup_trimming(&cli, &cli.input[0])?;
        run_paired_ubam_single_file(
            &cli,
            &cli.input[0],
            &outputs,
            &config,
            gzip,
            &adapters_r1,
            &adapters_r2,
        )?;
        return Ok(());
    }

    if cli.paired {
        // Pre-flight across pairs before any I/O; see io::collision_key for the key.
        let mut planned: Vec<naming::PlannedOutput> = Vec::new();
        for chunk in cli.input.chunks(2) {
            let (o1, o2) = naming::paired_end_output_names(
                &chunk[0],
                &chunk[1],
                output_dir,
                cli.basename.as_deref(),
                gzip,
            );
            // Uniform Pair for paired primaries (#397 decision C2): _val_1 takes its
            // stem from R1 but its directory from R1 too, and _val_2 mixes both — naming
            // one mate would encode a claim about which supplies the directory.
            let pair_src = naming::OutputSource::Pair(chunk[0].clone(), chunk[1].clone());
            let mut candidates = vec![(o1, pair_src.clone()), (o2, pair_src.clone())];
            if cli.retain_unpaired {
                let (u1, u2) = naming::unpaired_output_names(
                    &chunk[0],
                    &chunk[1],
                    output_dir,
                    cli.basename.as_deref(),
                    gzip,
                );
                candidates.push((u1, pair_src.clone()));
                candidates.push((u2, pair_src));
            }
            // --passthrough adds a third output path per pair. v1 only
            // supports a single pair (Cli::validate enforces input.len() == 2
            // when passthrough is set), so this either contributes zero or
            // one extra candidate to the collision set.
            if let Some(ref pt_input) = cli.passthrough {
                candidates.push((
                    naming::passthrough_output_name(
                        &chunk[0],
                        pt_input,
                        output_dir,
                        cli.basename.as_deref(),
                        gzip,
                    ),
                    naming::OutputSource::Input(pt_input.clone()),
                ));
            }
            // #388 — report names carry no _val_ discriminator, so two inputs
            // with distinct primaries can still collide on reports.
            //
            // `Input`, not `Pair`, even though the directory now comes from R1 (#398):
            // two distinct sources on one path select the message that says which two
            // inputs collided, where `Pair` would select the "list each input once"
            // text — wrong advice when renaming one input does fix it.
            if !cli.no_report_file {
                for input in [&chunk[0], &chunk[1]] {
                    let src = naming::OutputSource::Input(input.clone());
                    let (txt, json) = naming::paired_report_names(input, &chunk[0], output_dir);
                    candidates.push((txt, src.clone()));
                    candidates.push((json, src));
                }
            }
            planned.extend(candidates);
        }
        naming::preflight_output_collisions(
            &planned,
            &guarded_inputs(&cli),
            Some(PAIRED_REPORT_HINT),
        )?;

        // Adapter detection runs PER PAIR — intentional deviation from Perl
        // v0.6.x (which detected once on $ARGV[0] at trim_galore:2455). Shell-
        // glob invocations across mixed library types or 2-colour/4-colour
        // chemistries are common enough that per-pair detection is the safer
        // default. The header peek is microseconds per sample, so the cost is
        // negligible on I/O-bound workloads. Symmetrical with the single-end
        // loop below, which has always detected per file.
        let total_pairs = cli.input.len() / 2;
        for (pair_idx, chunk) in cli.input.chunks(2).enumerate() {
            if total_pairs > 1 {
                eprintln!("\n=== Pair {} of {} ===", pair_idx + 1, total_pairs);
            }
            // R1 of pair 0 was already sanity-checked at line 28; skip that.
            // R2 of pair 0 and both files of every later pair need checking
            // before setup_trimming reads the header.
            if pair_idx > 0 {
                sanity_check_any(&chunk[0])?;
            }
            sanity_check_any(&chunk[1])?;
            // --passthrough: sanity-check the third file once (Cli::validate
            // already enforces single-pair-only when passthrough is set).
            // `sanity_check_any`, not `FastqReader::sanity_check`, so the
            // passthrough path reaches the #379 restartability guard.
            if pair_idx == 0
                && let Some(ref pt_path) = cli.passthrough
            {
                sanity_check_any(pt_path)?;
            }

            let (_label, adapters_r1, adapters_r2, config) = setup_trimming(&cli, &chunk[0])?;

            run_paired(
                &cli,
                &chunk[0],
                &chunk[1],
                &config,
                gzip,
                output_dir,
                cli.basename.as_deref(),
                &adapters_r1,
                &adapters_r2,
            )
            .with_context(|| {
                format!(
                    "processing pair {} of {} (R1={}, R2={})",
                    pair_idx + 1,
                    total_pairs,
                    chunk[0].display(),
                    chunk[1].display()
                )
            })?;
        }
    } else {
        // Single-end: process each input file independently
        // (matches Perl TrimGalore behavior of looping over all positional args)
        // #383 — SE trim was the only trim path without the #216 pre-flight.
        let planned: Vec<naming::PlannedOutput> = cli
            .input
            .iter()
            .map(|input| {
                (
                    naming::single_end_output_name(
                        input,
                        output_dir,
                        cli.basename.as_deref(),
                        gzip,
                    ),
                    naming::OutputSource::Input(input.clone()),
                )
            })
            .collect();
        let mut planned = planned;
        planned.extend(planned_secondary_outputs(&cli, output_dir, gzip)?);
        naming::preflight_output_collisions(&planned, &guarded_inputs(&cli), None)?;
        for (i, input) in cli.input.iter().enumerate() {
            if i > 0 {
                eprintln!("\n--------------------------------------------------");
                sanity_check_any(input)?;
            }
            let (_label, adapters_r1, _, config) = setup_trimming(&cli, input)?;
            run_single_file(&cli, input, &config, gzip, output_dir, &adapters_r1)?;
        }
    }

    Ok(())
}

/// Set up adapter detection, poly-G scanning, and build TrimConfig for one input file.
///
/// When adapter is user-specified, auto-detection is skipped. When auto-detecting,
/// the poly-G scan piggybacks on the adapter scan. Returns (adapter_label,
/// adapters_r1, adapters_r2, TrimConfig).
fn setup_trimming(cli: &Cli, input_file: &Path) -> SetupResult {
    // Determine adapter
    let (adapter_label, adapters_r1, mut adapters_r2, autodetect_poly_g) =
        resolve_adapter(cli, input_file)?;

    // #369 — a user -a2 wins over the preset/auto-detected Read 2 candidate.
    let displaced_r2 = apply_adapter2_override(cli, &adapters_r1, &mut adapters_r2)?;

    // Display adapter info
    if adapters_r1.len() == 1 {
        eprintln!("Adapter: {} ({})", adapter_label, adapters_r1[0].1);
    } else {
        eprintln!(
            "Adapters ({}, {} sequences):",
            adapter_label,
            adapters_r1.len()
        );
        for (name, seq) in &adapters_r1 {
            eprintln!("  {}: {}", name, seq);
        }
    }
    // Read 2 adapters are only used in paired mode.
    if cli.paired && !adapters_r2.is_empty() {
        if adapters_r2.len() == 1 {
            eprintln!("Adapter 2 (Read 2): {}", adapters_r2[0].1);
        } else {
            eprintln!("Adapters R2 ({} sequences):", adapters_r2.len());
            for (name, seq) in &adapters_r2 {
                eprintln!("  {}: {}", name, seq);
            }
        }
    }
    if let Some(seq) = &displaced_r2 {
        eprintln!(
            "NOTE: Read 2 adapter taken from -a2; the {adapter_label} default ({seq}) is not used."
        );
    }
    if cli.times > 1 {
        eprintln!("Adapter trimming rounds per read (-n): {}", cli.times);
    }

    // Resolve length cutoff: smallRNA adapter auto-reduces to 18bp
    let first_adapter_seq = adapters_r1.first().map(|(_, s)| s.as_str()).unwrap_or("");
    let length_cutoff = cli.length.unwrap_or_else(|| {
        if first_adapter_seq == "TGGAATTCTCGG" {
            eprintln!("Reducing length cutoff to 18bp for small RNA-Seq reads because a cutoff of 20bp may remove some short species of small RNAs if they had been trimmed by 1,2 or 3bp");
            18
        } else {
            20
        }
    });
    eprintln!();

    // Build max_n filter. Values in (0.0, 1.0) are interpreted as a fraction
    // of the read length, matching Perl v0.6.8+ behaviour. The fraction case
    // is easy to enter accidentally (e.g. typing `--max_n 0.5` when meaning
    // "half a read"), so emit the same warning Perl does so users can see
    // which mode their invocation actually selected. See issue #243.
    let max_n = cli.max_n.map(|v| {
        if v >= 1.0 {
            MaxNFilter::Count(v as usize)
        } else {
            eprintln!("--max_n will be interpreted as a fraction of the read length ({v})");
            MaxNFilter::Fraction(v)
        }
    });

    // #440 — fold a `--library` preset into the four clip flags. Every clip
    // consumer below reads `clips`, never `cli.clip_*`, so a preset cannot reach
    // one path and miss another.
    let clips = cli.effective_clips();
    if let Some(preset) = clips.preset {
        eprintln!(
            "Library preset '{}' selected: {}",
            preset.canonical_name(),
            clips.flag_summary(cli.paired)
        );
        for o in &clips.overrides {
            eprintln!(
                "{} {} was given on the command line and overrides the {} preset value {}",
                o.flag,
                o.user_value,
                preset.canonical_name(),
                o.preset_value
            );
        }
    }

    // RRBS: auto-set --clip_r2 2 for directional paired-end mode. A preset that
    // supplies clip_R2 counts as "already set", same as an explicit flag.
    let clip_r2 = if cli.rrbs && !cli.non_directional && cli.paired && clips.clip_r2.is_none() {
        eprintln!(
            "Setting the option '--clip_r2 2' (to remove methylation bias from the start of Read 2)"
        );
        Some(2)
    } else {
        clips.clip_r2
    };

    if cli.rrbs {
        eprintln!(
            "File was specified to be an MspI-digested RRBS sample. Read 1 sequences with adapter contamination will be trimmed a further 2 bp from their 3' end, and Read 2 sequences will be trimmed by 2 bp from their 5' end to remove potential methylation-biased bases from the end-repair reaction"
        );
    }
    if cli.non_directional {
        eprintln!(
            "File was specified to be a non-directional MspI-digested RRBS sample. Sequences starting with either 'CAA' or 'CGA' will have the first 2 bp trimmed off to remove potential methylation-biased bases from the end-repair reaction"
        );
    }

    // Determine poly-G trimming: CLI overrides auto-detection
    let poly_g_enabled = if cli.poly_g {
        true
    } else if cli.no_poly_g {
        false
    } else {
        let (poly_g_count, reads_scanned) = if let Some((count, scanned)) = autodetect_poly_g {
            (count, scanned)
        } else {
            eprintln!("Scanning for poly-G content...");
            adapter::detect_poly_g(input_file)?
        };
        let threshold = (reads_scanned / 10_000).max(10);
        let enabled = poly_g_count > threshold;

        let poly_g_pct = if reads_scanned > 0 {
            poly_g_count as f64 / reads_scanned as f64 * 100.0
        } else {
            0.0
        };

        if enabled {
            eprintln!(
                "Poly-G trimming: ENABLED (auto-detected). \
                 {} of {} reads ({:.2}%) have poly-G tails (>=10bp) — \
                 consistent with 2-colour chemistry (NovaSeq/NextSeq). \
                 To disable: --no_poly_g",
                poly_g_count, reads_scanned, poly_g_pct
            );
        } else {
            eprintln!(
                "Poly-G trimming: not enabled (auto-detection found {} of {} reads ({:.2}%) \
                 with poly-G tails — below threshold). To force-enable: --poly_g",
                poly_g_count, reads_scanned, poly_g_pct
            );
        }
        enabled
    };

    if cli.poly_g {
        eprintln!("Poly-G trimming: ENABLED (user-specified --poly_g)");
    } else if cli.no_poly_g {
        eprintln!("Poly-G trimming: DISABLED (user-specified --no_poly_g)");
    }

    // Convert string adapters to bytes for the trimmer config
    let adapters_bytes: Vec<(String, Vec<u8>)> = adapters_r1
        .iter()
        .map(|(name, seq)| (name.clone(), seq.as_bytes().to_vec()))
        .collect();
    let adapters_r2_bytes: Vec<(String, Vec<u8>)> = adapters_r2
        .iter()
        .map(|(name, seq)| (name.clone(), seq.as_bytes().to_vec()))
        .collect();

    let config = trimmer::TrimConfig {
        adapters: adapters_bytes,
        adapters_r2: adapters_r2_bytes,
        times: cli.times,
        quality_cutoff: cli.effective_quality_cutoff(),
        phred_offset: cli.phred_offset(),
        error_rate: cli.error_rate,
        min_overlap: cli.stringency,
        length_cutoff,
        max_length: cli.max_length,
        max_n,
        trim_n: cli.trim_n,
        clip_r1: clips.clip_r1,
        clip_r2,
        three_prime_clip_r1: clips.three_prime_clip_r1,
        three_prime_clip_r2: clips.three_prime_clip_r2,
        rename: cli.rename,
        nextseq: cli.nextseq.is_some(),
        rrbs: cli.rrbs,
        non_directional: cli.non_directional,
        is_paired: cli.paired,
        poly_a: cli.poly_a,
        poly_g: poly_g_enabled,
        discard_untrimmed: cli.discard_untrimmed,
        gzip_level: cli.compression,
    };

    Ok((adapter_label, adapters_r1, adapters_r2, config))
}

/// Apply a user `-a2` over the preset/auto-detected Read 2 candidate (#369).
///
/// Returns the displaced default sequence, if the candidate was non-empty, so the
/// caller can report it. The override is skipped under `--consider_already_trimmed`
/// suppression: trimming R2 while R1 is left alone would be asymmetric, and the
/// mode announces that only quality trimming will happen.
fn apply_adapter2_override(
    cli: &Cli,
    adapters_r1: &AdapterList,
    adapters_r2: &mut AdapterList,
) -> Result<Option<String>> {
    // Not paired: `Cli::validate` has already warned, and Read 2 does not exist.
    if cli.adapter2.is_empty() || !cli.paired {
        return Ok(None);
    }

    // Suppression is the only way an R1 adapter reaches here with an empty sequence.
    let suppressed = adapters_r1.len() == 1 && adapters_r1[0].1.is_empty();
    if suppressed {
        eprintln!(
            "WARNING: -a2/--adapter2 not applied — adapter trimming is suppressed for this \
             library (--consider_already_trimmed). Ignoring."
        );
        return Ok(None);
    }

    let displaced = adapters_r2.first().map(|(_, seq)| seq.clone());
    *adapters_r2 = adapter::parse_adapter_specs(&cli.adapter2)?;
    Ok(displaced)
}

/// Returns (adapter_label, adapters_r1, adapters_r2, poly_g_from_autodetect).
/// `adapter_label` is for display purposes (e.g., "Illumina", "user-specified").
/// `adapters_r2` is a *candidate* — `apply_adapter2_override` may replace it with `-a2`.
/// `adapters_r1`/`adapters_r2` are `(name, sequence)` pairs; r2 is empty if not set.
/// The last element is Some((poly_g_count, reads_scanned)) when auto-detection ran,
/// None when the adapter was user-specified or preset-selected (poly-G must be
/// detected separately via `adapter::detect_poly_g()`).
fn resolve_adapter(cli: &Cli, input_file: &Path) -> ResolvedAdapter {
    if !cli.adapter.is_empty() {
        let adapters_r1 = adapter::parse_adapter_specs(&cli.adapter)?;
        let label = "user-specified".to_string();
        return Ok((label, adapters_r1, Vec::new(), None));
    }

    // Presets: single-adapter, use AdapterPreset methods
    if cli.nextera {
        return Ok((
            adapter::NEXTERA.name.to_string(),
            adapter::NEXTERA.to_adapter_vec(),
            Vec::new(),
            None,
        ));
    }
    if cli.small_rna {
        return Ok((
            adapter::SMALL_RNA.name.to_string(),
            adapter::SMALL_RNA.to_adapter_vec(),
            adapter::SMALL_RNA.to_r2_vec(),
            None,
        ));
    }
    if cli.stranded_illumina {
        return Ok((
            adapter::STRANDED_ILLUMINA.name.to_string(),
            adapter::STRANDED_ILLUMINA.to_adapter_vec(),
            Vec::new(),
            None,
        ));
    }
    if cli.bgiseq {
        return Ok((
            adapter::BGISEQ.name.to_string(),
            adapter::BGISEQ.to_adapter_vec(),
            adapter::BGISEQ.to_r2_vec(),
            None,
        ));
    }
    if cli.illumina {
        return Ok((
            adapter::ILLUMINA.name.to_string(),
            adapter::ILLUMINA.to_adapter_vec(),
            Vec::new(),
            None,
        ));
    }

    // Auto-detect (also piggybacks poly-G counting)
    eprintln!("Auto-detecting adapter type...");
    let detection = adapter::autodetect_adapter(input_file, cli.consider_already_trimmed)?;
    eprintln!("{}", detection.message);

    let poly_g_data = Some((detection.poly_g_count, detection.reads_scanned));
    let adapters_r1 = detection.adapter.to_adapter_vec();
    let adapters_r2 = detection.adapter.to_r2_vec();
    Ok((
        detection.adapter.name.to_string(),
        adapters_r1,
        adapters_r2,
        poly_g_data,
    ))
}

fn run_single_file(
    cli: &Cli,
    input: &Path,
    config: &trimmer::TrimConfig,
    gzip: bool,
    output_dir: Option<&Path>,
    adapters_r1: &[(String, String)],
) -> Result<()> {
    let output_path =
        naming::single_end_output_name(input, output_dir, cli.basename.as_deref(), gzip);
    let report_path = naming::report_name(input, output_dir);

    eprintln!("Trimming: {}", input.display());
    eprintln!("Output:   {}", output_path.display());

    let stats = if cli.cores > 1 || cli.clumpify {
        // Worker-pool parallel path: N workers each handle trim + compress.
        // `--clumpify` always routes here (validation enforces cores >= 2).
        let reader = open_threaded_reader(input, &cli.preserve_tags)?;
        let clump_layout = resolve_clump_layout(cli)?;
        parallel::run_single_end_parallel(
            reader,
            &output_path,
            config,
            cli.cores,
            gzip,
            clump_layout,
        )?
    } else {
        let mut reader = open_sync_reader(input, &cli.preserve_tags)?;
        let mut writer = FastqWriter::create(&output_path, gzip, 1, config.gzip_level)?;
        let stats = trimmer::run_single_end(reader.as_mut(), &mut writer, config)?;
        writer.finish()?;
        stats
    };

    // Print summary
    eprintln!("\n=== Summary ===\n");
    eprintln!("Total reads processed:           {:>10}", stats.total_reads);
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats.total_reads_with_adapter,
        pct(stats.total_reads_with_adapter, stats.total_reads)
    );
    if stats.discarded_untrimmed > 0 {
        eprintln!(
            "Reads discarded as untrimmed:    {:>10} ({:.1}%)",
            stats.discarded_untrimmed,
            pct(stats.discarded_untrimmed, stats.total_reads)
        );
    }
    eprintln!(
        "Reads too short:                 {:>10} ({:.1}%)",
        stats.too_short,
        pct(stats.too_short, stats.total_reads)
    );
    eprintln!(
        "Reads written (passing filters): {:>10} ({:.1}%)",
        stats.reads_written,
        pct(stats.reads_written, stats.total_reads)
    );
    if stats.rrbs_trimmed_3prime > 0 {
        eprintln!(
            "RRBS trimmed (3' end, adapter): {:>10} ({:.1}%)",
            stats.rrbs_trimmed_3prime,
            pct(stats.rrbs_trimmed_3prime, stats.total_reads)
        );
    }
    if stats.rrbs_trimmed_5prime > 0 {
        eprintln!(
            "RRBS trimmed (5' end, CAA/CGA): {:>10} ({:.1}%)",
            stats.rrbs_trimmed_5prime,
            pct(stats.rrbs_trimmed_5prime, stats.total_reads)
        );
    }
    if stats.poly_a_trimmed > 0 {
        eprintln!(
            "Reads with poly-A/T trimmed:     {:>10} ({:.1}%)",
            stats.poly_a_trimmed,
            pct(stats.poly_a_trimmed, stats.total_reads)
        );
        eprintln!(
            "  Poly-A/T bases removed:        {:>10}",
            stats.poly_a_bases_trimmed
        );
    }
    if stats.poly_g_trimmed > 0 {
        eprintln!(
            "Reads with poly-G/C trimmed:     {:>10} ({:.1}%)",
            stats.poly_g_trimmed,
            pct(stats.poly_g_trimmed, stats.total_reads)
        );
        eprintln!(
            "  Poly-G/C bases removed:        {:>10}",
            stats.poly_g_bases_trimmed
        );
    }

    // Write report
    if !cli.no_report_file {
        let input_filename = input
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let report_cfg = report::TrimConfig {
            version: env!("CARGO_PKG_VERSION").to_string(),
            quality_cutoff: cli.effective_quality_cutoff(),
            adapters: adapters_r1.to_vec(),
            adapters_r2: Vec::new(),
            times: cli.times,
            error_rate: cli.error_rate,
            stringency: cli.stringency,
            length_cutoff: config.length_cutoff,
            max_length: cli.max_length,
            paired: false,
            gzip,
            trim_n: cli.trim_n,
            nextseq: cli.nextseq.is_some(),
            rrbs: cli.rrbs,
            non_directional: cli.non_directional,
            phred_encoding: cli.phred_offset(),
            poly_a: cli.poly_a,
            poly_g: config.poly_g,
            command_line: std::env::args().collect::<Vec<_>>().join(" "),
            input_filename: input_filename.clone(),
            input_filenames: vec![input_filename.clone()],
            library: Some(cli.effective_clips()),
        };

        let file = File::create(&report_path)?;
        let mut w = BufWriter::new(file);
        report::write_report_header(&mut w, &report_cfg)?;
        report::write_cutadapt_compatible_section(&mut w, &report_cfg, &stats, 1)?;
        report::write_run_footer(&mut w, &report_cfg, &stats)?;
        eprintln!("\nReport: {}", report_path.display());

        // JSON report (use effective config values, not raw CLI values,
        // because setup_trimming() may override e.g. clip_r2 for RRBS)
        let json_path = naming::json_report_name(input, output_dir);
        let json_extra = report::JsonReportParams {
            clip_r1: config.clip_r1,
            clip_r2: config.clip_r2,
            three_prime_clip_r1: config.three_prime_clip_r1,
            three_prime_clip_r2: config.three_prime_clip_r2,
            max_n: cli.max_n,
            discard_untrimmed: config.discard_untrimmed,
            consider_already_trimmed: cli.consider_already_trimmed,
        };
        let json_file = File::create(&json_path)?;
        let mut jw = BufWriter::new(json_file);
        report::write_json_report(&mut jw, &report_cfg, &stats, None, 1, &json_extra)?;
        jw.flush()?;
        eprintln!("JSON report: {}", json_path.display());
    }

    // Run FastQC if requested (bundled fastqc-rust library — no shell-out)
    if cli.fastqc_requested() {
        fastqc::run(
            &output_path,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
    }

    // Demultiplex if requested
    if let Some(ref barcode_file) = cli.demux {
        eprintln!(
            "\nTrimming complete, starting demultiplexing procedure (based on 3' barcodes supplied as per file >{}<)",
            barcode_file.display()
        );
        let barcodes = demux::read_barcode_file(barcode_file)?;
        demux::demultiplex(
            &output_path,
            &barcodes,
            gzip,
            output_dir,
            cli.cores,
            cli.compression,
        )?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_paired(
    cli: &Cli,
    input_r1: &Path,
    input_r2: &Path,
    config: &trimmer::TrimConfig,
    gzip: bool,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    let (output_r1, output_r2) =
        naming::paired_end_output_names(input_r1, input_r2, output_dir, basename, gzip);

    eprintln!("Trimming (paired-end):");
    eprintln!("  R1: {}", input_r1.display());
    eprintln!("  R2: {}", input_r2.display());
    eprintln!("  Output R1: {}", output_r1.display());
    eprintln!("  Output R2: {}", output_r2.display());

    // --passthrough wiring (plan v2 Step 8). Compute the third output path
    // when active and eprintln! it for visibility, matching the R1/R2 idiom.
    let passthrough_input: Option<&Path> = cli.passthrough.as_deref();
    let passthrough_output: Option<std::path::PathBuf> = passthrough_input
        .map(|pt| naming::passthrough_output_name(input_r1, pt, output_dir, basename, gzip));
    if let (Some(pt_in), Some(pt_out)) = (passthrough_input, passthrough_output.as_deref()) {
        eprintln!("  Passthrough: {} → {}", pt_in.display(), pt_out.display());
    }

    // Compute unpaired output paths (needed for both parallel and sequential paths)
    let (unpaired_r1_path, unpaired_r2_path) = if cli.retain_unpaired {
        let (up1, up2) =
            naming::unpaired_output_names(input_r1, input_r2, output_dir, basename, gzip);
        eprintln!("  Unpaired R1: {}", up1.display());
        eprintln!("  Unpaired R2: {}", up2.display());
        (Some(up1), Some(up2))
    } else {
        (None, None)
    };

    // Internal-invariant backstop, NOT the user-facing rejection. Any BAM in a
    // two-file pair is rejected by `reject_bam_format_mismatch_in_pair` in
    // main(), which produces the shape-appropriate message (#363).
    //
    // This stays as an enforced check rather than a comment because the
    // sequential path below hard-codes `FastqReader::open` — and `--cores 1` is
    // the default. A BGZF BAM handed to `FastqReader` decompresses to binary and
    // parses as FASTQ, i.e. silent wrong output rather than an error. Deliberately
    // worded so it cannot be mistaken for the user-facing message.
    for p in [input_r1, input_r2] {
        if matches!(detect_input_format(p)?, InputFormat::UnalignedBam) {
            anyhow::bail!(
                "internal error: uBAM input reached run_paired ({}); the paired \
                 format guard in main() should have rejected it. Please report this at \
                 https://github.com/FelixKrueger/TrimGalore/issues",
                p.display()
            );
        }
    }

    let (stats_r1, stats_r2, pair_stats) = if cli.cores > 1 || cli.clumpify {
        // Worker-pool parallel path: N workers each handle trim + compress.
        // `--clumpify` always routes here (validation enforces cores >= 2).
        let reader_r1 = open_threaded_reader(input_r1, &cli.preserve_tags)?;
        let reader_r2 = open_threaded_reader(input_r2, &cli.preserve_tags)?;
        let reader_passthrough: Option<Box<dyn RecordSource>> = match passthrough_input {
            Some(p) => Some(Box::new(FastqReader::open_threaded(p)?)),
            None => None,
        };
        let clump_layout = resolve_clump_layout(cli)?;
        parallel::run_paired_end_parallel(
            reader_r1,
            reader_r2,
            reader_passthrough,
            &output_r1,
            &output_r2,
            passthrough_output.as_deref(),
            unpaired_r1_path.as_deref(),
            unpaired_r2_path.as_deref(),
            config,
            cli.cores,
            gzip,
            UnpairedLengths {
                r1: cli.length_1,
                r2: cli.length_2,
            },
            clump_layout,
        )?
    } else {
        // Sequential path (--cores 1, the default) — hard-codes FastqReader,
        // which is safe only because uBAM input is rejected by
        // `reject_bam_format_mismatch_in_pair` in main() and re-checked by the
        // internal-invariant backstop at the top of this function. Do not
        // remove that backstop: a BGZF BAM handed to FastqReader parses as
        // FASTQ rather than erroring.
        let mut reader_r1 = FastqReader::open(input_r1)?;
        let mut reader_r2 = FastqReader::open(input_r2)?;
        let level = config.gzip_level;
        let mut writer_r1 = FastqWriter::create(&output_r1, gzip, 1, level)?;
        let mut writer_r2 = FastqWriter::create(&output_r2, gzip, 1, level)?;

        // --passthrough: open the third reader + writer when active.
        // Plan v2 §Behavior Step 5/6: serial path opens via FastqReader::open
        // (not open_threaded — single-threaded mode).
        let mut reader_passthrough = match passthrough_input {
            Some(p) => Some(FastqReader::open(p)?),
            None => None,
        };
        let mut writer_passthrough = match passthrough_output.as_deref() {
            Some(p) => Some(FastqWriter::create(p, gzip, 1, level)?),
            None => None,
        };

        let (mut unpaired_w1, mut unpaired_w2) = match (&unpaired_r1_path, &unpaired_r2_path) {
            (Some(p1), Some(p2)) => (
                Some(FastqWriter::create(p1, gzip, 1, level)?),
                Some(FastqWriter::create(p2, gzip, 1, level)?),
            ),
            _ => (None, None),
        };

        let result = trimmer::run_paired_end(
            &mut reader_r1,
            &mut reader_r2,
            reader_passthrough
                .as_mut()
                .map(|r| r as &mut dyn RecordSource),
            &mut writer_r1,
            &mut writer_r2,
            writer_passthrough.as_mut(),
            unpaired_w1.as_mut(),
            unpaired_w2.as_mut(),
            config,
            UnpairedLengths {
                r1: cli.length_1,
                r2: cli.length_2,
            },
        )?;

        // Fixed order, and nothing fallible between the first and the last, so
        // the window in which only some of the set exists is the renames alone.
        writer_r1.finish()?;
        writer_r2.finish()?;
        if let Some(w) = writer_passthrough {
            w.finish()?;
        }
        if let Some(w) = unpaired_w1 {
            w.finish()?;
        }
        if let Some(w) = unpaired_w2 {
            w.finish()?;
        }

        result
    };

    // Print summary
    eprintln!("\n=== Summary (Read 1) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r1.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r1.total_reads_with_adapter,
        pct(stats_r1.total_reads_with_adapter, stats_r1.total_reads)
    );

    eprintln!("\n=== Summary (Read 2) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r2.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r2.total_reads_with_adapter,
        pct(stats_r2.total_reads_with_adapter, stats_r2.total_reads)
    );

    if stats_r1.poly_a_trimmed > 0 {
        eprintln!(
            "R1 reads with poly-A trimmed:    {:>10} ({:.1}%)",
            stats_r1.poly_a_trimmed,
            pct(stats_r1.poly_a_trimmed, stats_r1.total_reads)
        );
    }
    if stats_r2.poly_a_trimmed > 0 {
        eprintln!(
            "R2 reads with poly-T trimmed:    {:>10} ({:.1}%)",
            stats_r2.poly_a_trimmed,
            pct(stats_r2.poly_a_trimmed, stats_r2.total_reads)
        );
    }
    if stats_r1.poly_g_trimmed > 0 {
        eprintln!(
            "R1 reads with poly-G trimmed:    {:>10} ({:.1}%)",
            stats_r1.poly_g_trimmed,
            pct(stats_r1.poly_g_trimmed, stats_r1.total_reads)
        );
        eprintln!(
            "  R1 poly-G bases removed:       {:>10}",
            stats_r1.poly_g_bases_trimmed
        );
    }
    if stats_r2.poly_g_trimmed > 0 {
        eprintln!(
            "R2 reads with poly-G trimmed:    {:>10} ({:.1}%)",
            stats_r2.poly_g_trimmed,
            pct(stats_r2.poly_g_trimmed, stats_r2.total_reads)
        );
        eprintln!(
            "  R2 poly-G bases removed:       {:>10}",
            stats_r2.poly_g_bases_trimmed
        );
    }

    eprintln!("\n=== Paired-end validation ===\n");
    eprintln!(
        "Pairs analyzed:                  {:>10}",
        pair_stats.pairs_analyzed
    );
    eprintln!(
        "Pairs removed:                   {:>10} ({:.1}%)",
        pair_stats.pairs_removed,
        pct(pair_stats.pairs_removed, pair_stats.pairs_analyzed)
    );
    if pair_stats.r1_unpaired > 0 || pair_stats.r2_unpaired > 0 {
        eprintln!(
            "Unpaired R1 kept:                {:>10}",
            pair_stats.r1_unpaired
        );
        eprintln!(
            "Unpaired R2 kept:                {:>10}",
            pair_stats.r2_unpaired
        );
    }

    // Write reports — delegated to the shared `write_paired_reports` helper
    // so we don't drift against `run_paired_ubam_single_file`'s identical block.
    if !cli.no_report_file {
        let all_input_filenames: Vec<String> = [input_r1, input_r2]
            .iter()
            .map(|p| {
                p.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string()
            })
            .collect();

        // Both mates' reports land where the primaries do (#398). Each keeps its own
        // input-derived filename; only the directory is shared.
        let (r1_txt, r1_json) = naming::paired_report_names(input_r1, input_r1, output_dir);
        let (r2_txt, r2_json) = naming::paired_report_names(input_r2, input_r1, output_dir);
        let r1 = PairedReportFile {
            txt_path: r1_txt,
            json_path: r1_json,
            input_filename: all_input_filenames[0].clone(),
        };
        let r2 = PairedReportFile {
            txt_path: r2_txt,
            json_path: r2_json,
            input_filename: all_input_filenames[1].clone(),
        };

        let passthrough = match (cli.passthrough.as_deref(), passthrough_output.as_deref()) {
            (Some(pt_in), Some(pt_out)) => Some((pt_in, pt_out)),
            _ => None,
        };

        write_paired_reports(
            cli,
            config,
            gzip,
            r1,
            r2,
            &all_input_filenames,
            &stats_r1,
            &stats_r2,
            &pair_stats,
            adapters_r1,
            adapters_r2,
            passthrough,
        )?;
    }

    // Run FastQC if requested (bundled fastqc-rust library — no shell-out)
    if cli.fastqc_requested() {
        fastqc::run(
            &output_r1,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
        fastqc::run(
            &output_r2,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
        // --passthrough: also FastQC the carrier output. Note: cell-barcode
        // reads (16-28 bp of uniformly-structured sequence) will produce a
        // FastQC report with per-base content bias warnings — that's expected,
        // not a defect. The `--passthrough` help text documents this. (Plan
        // v2 Step 10.)
        if let Some(ref pt_out) = passthrough_output {
            fastqc::run(pt_out, cli.fastqc_args.as_deref(), output_dir, cli.cores)?;
        }
    }

    Ok(())
}

/// Paired-uBAM single-file pipeline (T20-PE). The input is one BAM file
/// containing interleaved R1/R2 records; `BamReader::open_paired_interleaved`
/// de-interleaves on the fly via the bounded `MAX_SLACK` buffer.
///
/// Output names come in from the caller, which built them for the pre-flight:
/// the input BAM's stem plus a suffix (`input.bam` → `input_val_1.fq[.gz]`).
#[allow(clippy::too_many_arguments)]
fn run_paired_ubam_single_file(
    cli: &Cli,
    input: &Path,
    outputs: &naming::InterleavedFastqOutputs,
    config: &trimmer::TrimConfig,
    gzip: bool,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    let output_r1 = outputs.val_1.clone();
    let output_r2 = outputs.val_2.clone();

    eprintln!("Trimming (paired-end, de-interleaved uBAM):");
    eprintln!("  Input:     {}", input.display());
    eprintln!("  Output R1: {}", output_r1.display());
    eprintln!("  Output R2: {}", output_r2.display());

    let (unpaired_r1_path, unpaired_r2_path) = match outputs.unpaired.clone() {
        Some((up1, up2)) => {
            eprintln!("  Unpaired R1: {}", up1.display());
            eprintln!("  Unpaired R2: {}", up2.display());
            (Some(up1), Some(up2))
        }
        None => (None, None),
    };

    let (stats_r1, stats_r2, pair_stats) = if cli.cores > 1 || cli.clumpify {
        let (r1, r2) = BamReader::open_paired_interleaved_with_tags(input, &cli.preserve_tags)?;
        let clump_layout = resolve_clump_layout(cli)?;
        parallel::run_paired_end_parallel(
            Box::new(r1),
            Box::new(r2),
            None, // --passthrough + BAM rejected earlier
            &output_r1,
            &output_r2,
            None,
            unpaired_r1_path.as_deref(),
            unpaired_r2_path.as_deref(),
            config,
            cli.cores,
            gzip,
            UnpairedLengths {
                r1: cli.length_1,
                r2: cli.length_2,
            },
            clump_layout,
        )?
    } else {
        // Serial path
        let (mut r1, mut r2) =
            BamReader::open_paired_interleaved_with_tags(input, &cli.preserve_tags)?;
        let level = config.gzip_level;
        let mut writer_r1 = FastqWriter::create(&output_r1, gzip, 1, level)?;
        let mut writer_r2 = FastqWriter::create(&output_r2, gzip, 1, level)?;
        let (mut unpaired_w1, mut unpaired_w2) = match (&unpaired_r1_path, &unpaired_r2_path) {
            (Some(p1), Some(p2)) => (
                Some(FastqWriter::create(p1, gzip, 1, level)?),
                Some(FastqWriter::create(p2, gzip, 1, level)?),
            ),
            _ => (None, None),
        };
        let result = trimmer::run_paired_end(
            &mut r1,
            &mut r2,
            None,
            &mut writer_r1,
            &mut writer_r2,
            None,
            unpaired_w1.as_mut(),
            unpaired_w2.as_mut(),
            config,
            UnpairedLengths {
                r1: cli.length_1,
                r2: cli.length_2,
            },
        )?;
        // Fixed order, nothing fallible in between — see run_paired.
        writer_r1.finish()?;
        writer_r2.finish()?;
        if let Some(w) = unpaired_w1 {
            w.finish()?;
        }
        if let Some(w) = unpaired_w2 {
            w.finish()?;
        }
        result
    };

    // Summary (mirrors run_paired's print shape).
    eprintln!("\n=== Summary (Read 1) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r1.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r1.total_reads_with_adapter,
        pct(stats_r1.total_reads_with_adapter, stats_r1.total_reads)
    );
    eprintln!("\n=== Summary (Read 2) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r2.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r2.total_reads_with_adapter,
        pct(stats_r2.total_reads_with_adapter, stats_r2.total_reads)
    );
    eprintln!("\n=== Paired-end validation ===\n");
    eprintln!(
        "Pairs analyzed:                  {:>10}",
        pair_stats.pairs_analyzed
    );
    eprintln!(
        "Pairs removed:                   {:>10} ({:.1}%)",
        pair_stats.pairs_removed,
        pct(pair_stats.pairs_removed, pair_stats.pairs_analyzed)
    );
    if pair_stats.r1_unpaired > 0 || pair_stats.r2_unpaired > 0 {
        eprintln!(
            "Unpaired R1 kept:                {:>10}",
            pair_stats.r1_unpaired
        );
        eprintln!(
            "Unpaired R2 kept:                {:>10}",
            pair_stats.r2_unpaired
        );
    }

    // Reports — text + JSON, via the shared write_paired_reports helper.
    // uBAM single-file paired-mode: only one input filename to record, but
    // both R1 and R2 report paths derive from a stem (_R1/_R2 suffix injected).
    if !cli.no_report_file {
        let input_filename = input
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let all_input_filenames = vec![input_filename.clone()];

        let reports = outputs
            .reports
            .as_ref()
            .expect("reports were requested, so the namer built them");
        let r1 = PairedReportFile {
            txt_path: reports.r1_txt.clone(),
            json_path: reports.r1_json.clone(),
            input_filename: input_filename.clone(),
        };
        let r2 = PairedReportFile {
            txt_path: reports.r2_txt.clone(),
            json_path: reports.r2_json.clone(),
            input_filename: input_filename.clone(),
        };

        // No --passthrough on the uBAM-paired path (rejected at Cli::validate).
        write_paired_reports(
            cli,
            config,
            gzip,
            r1,
            r2,
            &all_input_filenames,
            &stats_r1,
            &stats_r2,
            &pair_stats,
            adapters_r1,
            adapters_r2,
            None,
        )?;
    }

    Ok(())
}

// ─── uBAM-output dispatch (PLAN v2.1 §5 step 4) ─────────────────────────────
//
// Three internal branches matching the FASTQ dispatch shape: paired-uBAM
// single-file (one interleaved input), paired multi-pair (two-file pairs),
// and single-end. Each opens readers + a `BamWriter` and calls the
// `trimmer::run_*_to_bam` entry points added in Step 3.

/// uBAM-output dispatch entry point.
///
/// The cores-ignored + missing-qual NOTEs are emitted up-front in `main()`
/// (before specialty dispatch), so they cover both this path and the
/// hardtrim BAM paths uniformly. Code-review round-2 B-NIT-2.
fn run_ubam_output(cli: &Cli, output_dir: Option<&Path>, command_line: &str) -> Result<()> {
    // Paired-uBAM single-file de-interleaved input.
    if cli.paired && cli.input.len() == 1 {
        let outputs =
            naming::InterleavedBamOutputs::new(&cli.input[0], output_dir, !cli.no_report_file);
        // `Input`, not `Pair`: one input means `Pair` could only be `Pair(x, x)`,
        // which renders `one.bam + one.bam`.
        let src = naming::OutputSource::Input(cli.input[0].clone());
        // Ahead of setup_trimming, which scans up to 1 M reads for adapters — a
        // refusal must not arrive after that wait.
        naming::preflight_output_collisions(
            &outputs.planned(&src),
            &guarded_inputs(cli),
            Some(PAIRED_REPORT_HINT),
        )?;
        let (_label, adapters_r1, adapters_r2, config) = setup_trimming(cli, &cli.input[0])?;
        return run_ubam_output_paired_single_file(
            cli,
            &cli.input[0],
            &outputs,
            &config,
            output_dir,
            command_line,
            &adapters_r1,
            &adapters_r2,
        );
    }

    if cli.paired {
        // Two-file paired-BAM input (and mixed FASTQ+BAM pairs) are rejected by
        // `reject_bam_format_mismatch_in_pair` in main(), which runs before this
        // dispatch and before the pre-flight below, and which distinguishes the
        // two-BAM case from the mixed case (#363). The per-pair loop that used to
        // stand here emitted the two-BAM message for either.
        // Pre-flight: one BAM output per pair; collision on case-folded path.
        let mut planned: Vec<naming::PlannedOutput> = Vec::new();
        for chunk in cli.input.chunks(2) {
            planned.push((
                naming::paired_bam_output_name(
                    &chunk[0],
                    &chunk[1],
                    output_dir,
                    cli.basename.as_deref(),
                ),
                naming::OutputSource::Pair(chunk[0].clone(), chunk[1].clone()),
            ));
            // #388 — same report-collision hole as the FASTQ paired path.
            if !cli.no_report_file {
                for input in [&chunk[0], &chunk[1]] {
                    let src = naming::OutputSource::Input(input.clone());
                    let (txt, json) = naming::paired_report_names(input, &chunk[0], output_dir);
                    planned.push((txt, src.clone()));
                    planned.push((json, src));
                }
            }
        }
        naming::preflight_output_collisions(
            &planned,
            &guarded_inputs(cli),
            Some(PAIRED_REPORT_HINT),
        )?;

        let total_pairs = cli.input.len() / 2;
        for (pair_idx, chunk) in cli.input.chunks(2).enumerate() {
            if total_pairs > 1 {
                eprintln!("\n=== Pair {} of {} ===", pair_idx + 1, total_pairs);
            }
            // R1 of pair 0 was already sanity-checked at line ~128.
            if pair_idx > 0 {
                sanity_check_any(&chunk[0])?;
            }
            sanity_check_any(&chunk[1])?;

            let (_label, adapters_r1, adapters_r2, config) = setup_trimming(cli, &chunk[0])?;
            run_ubam_output_paired_two_files(
                cli,
                &chunk[0],
                &chunk[1],
                &config,
                output_dir,
                command_line,
                &adapters_r1,
                &adapters_r2,
            )
            .with_context(|| {
                format!(
                    "processing pair {} of {} (R1={}, R2={})",
                    pair_idx + 1,
                    total_pairs,
                    chunk[0].display(),
                    chunk[1].display()
                )
            })?;
        }
        return Ok(());
    }

    // Single-end loop.
    // #383 — same hole as the FASTQ SE loop.
    let mut planned: Vec<naming::PlannedOutput> = cli
        .input
        .iter()
        .map(|input| {
            (
                naming::single_end_bam_output_name(input, output_dir, cli.basename.as_deref()),
                naming::OutputSource::Input(input.clone()),
            )
        })
        .collect();
    // #409 — run_ubam_output_single writes both trimming reports too; without them
    // an input named like one is overwritten before it is read.
    if !cli.no_report_file {
        for input in &cli.input {
            let src = naming::OutputSource::Input(input.clone());
            planned.push((naming::report_name(input, output_dir), src.clone()));
            planned.push((naming::json_report_name(input, output_dir), src));
        }
    }
    naming::preflight_output_collisions(&planned, &guarded_inputs(cli), None)?;
    for (i, input) in cli.input.iter().enumerate() {
        if i > 0 {
            eprintln!("\n--------------------------------------------------");
            sanity_check_any(input)?;
        }
        let (_label, adapters_r1, _, config) = setup_trimming(cli, input)?;
        run_ubam_output_single(cli, input, &config, output_dir, command_line, &adapters_r1)?;
    }
    Ok(())
}

/// uBAM-output single-end driver. Mirrors `run_single_file` for the FASTQ
/// path; opens one reader + one `BamWriter`, calls
/// [`trimmer::run_single_end_to_bam`].
fn run_ubam_output_single(
    cli: &Cli,
    input: &Path,
    config: &trimmer::TrimConfig,
    output_dir: Option<&Path>,
    command_line: &str,
    adapters_r1: &[(String, String)],
) -> Result<()> {
    let output_path =
        naming::single_end_bam_output_name(input, output_dir, cli.basename.as_deref());
    let report_path = naming::report_name(input, output_dir);

    eprintln!("Trimming: {}", input.display());
    eprintln!("Output:   {}", output_path.display());

    let source_header = match detect_input_format(input)? {
        InputFormat::UnalignedBam => Some(trim_galore::bam::peek_header(input)?),
        InputFormat::FastqPlain | InputFormat::FastqGz => None,
    };
    if source_header.is_none() {
        eprintln!("NOTE: input is FASTQ; output uBAM records will have no aux fields.");
    }

    let mut reader = open_sync_reader(input, &cli.preserve_tags)?;
    let mut writer = trim_galore::bam::BamWriter::create(
        &output_path,
        source_header.as_ref(),
        &cli.preserve_tags,
        command_line,
        cli.phred_offset(),
    )?;
    let stats = trimmer::run_single_end_to_bam(reader.as_mut(), &mut writer, config)?;
    writer.finish()?;

    // Summary (shape mirrors `run_single_file`).
    eprintln!("\n=== Summary ===\n");
    eprintln!("Total reads processed:           {:>10}", stats.total_reads);
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats.total_reads_with_adapter,
        pct(stats.total_reads_with_adapter, stats.total_reads)
    );
    if stats.discarded_untrimmed > 0 {
        eprintln!(
            "Reads discarded as untrimmed:    {:>10} ({:.1}%)",
            stats.discarded_untrimmed,
            pct(stats.discarded_untrimmed, stats.total_reads)
        );
    }
    eprintln!(
        "Reads too short:                 {:>10} ({:.1}%)",
        stats.too_short,
        pct(stats.too_short, stats.total_reads)
    );
    eprintln!(
        "Reads written (passing filters): {:>10} ({:.1}%)",
        stats.reads_written,
        pct(stats.reads_written, stats.total_reads)
    );

    // Reports (text + JSON), same shape as `run_single_file`.
    if !cli.no_report_file {
        let input_filename = input
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let report_cfg = report::TrimConfig {
            version: env!("CARGO_PKG_VERSION").to_string(),
            quality_cutoff: cli.effective_quality_cutoff(),
            adapters: adapters_r1.to_vec(),
            adapters_r2: Vec::new(),
            times: cli.times,
            error_rate: cli.error_rate,
            stringency: cli.stringency,
            length_cutoff: config.length_cutoff,
            max_length: cli.max_length,
            paired: false,
            gzip: false, // uBAM is BGZF-framed; the gzip flag here is FASTQ-output-specific.
            trim_n: cli.trim_n,
            nextseq: cli.nextseq.is_some(),
            rrbs: cli.rrbs,
            non_directional: cli.non_directional,
            phred_encoding: cli.phred_offset(),
            poly_a: cli.poly_a,
            poly_g: config.poly_g,
            command_line: command_line.to_string(),
            input_filename: input_filename.clone(),
            input_filenames: vec![input_filename.clone()],
            library: Some(cli.effective_clips()),
        };

        let file = File::create(&report_path)?;
        let mut w = BufWriter::new(file);
        report::write_report_header(&mut w, &report_cfg)?;
        report::write_cutadapt_compatible_section(&mut w, &report_cfg, &stats, 1)?;
        report::write_run_footer(&mut w, &report_cfg, &stats)?;
        eprintln!("\nReport: {}", report_path.display());

        let json_path = naming::json_report_name(input, output_dir);
        let json_extra = report::JsonReportParams {
            clip_r1: config.clip_r1,
            clip_r2: config.clip_r2,
            three_prime_clip_r1: config.three_prime_clip_r1,
            three_prime_clip_r2: config.three_prime_clip_r2,
            max_n: cli.max_n,
            discard_untrimmed: config.discard_untrimmed,
            consider_already_trimmed: cli.consider_already_trimmed,
        };
        let json_file = File::create(&json_path)?;
        let mut jw = BufWriter::new(json_file);
        report::write_json_report(&mut jw, &report_cfg, &stats, None, 1, &json_extra)?;
        jw.flush()?;
        eprintln!("JSON report: {}", json_path.display());
    }

    // Run FastQC if requested. fastqc-rust dispatches on file EXTENSION
    // (not content), so this relies on output_path ending in `.bam`.
    if cli.fastqc_requested() {
        fastqc::run(
            &output_path,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
    }

    Ok(())
}

/// uBAM-output two-file paired-end driver (FASTQ+FASTQ or BAM+BAM input
/// pairs). Output is ONE interleaved BAM per PLAN §3.2; both stats sides
/// flow into [`write_paired_reports`].
#[allow(clippy::too_many_arguments)]
fn run_ubam_output_paired_two_files(
    cli: &Cli,
    input_r1: &Path,
    input_r2: &Path,
    config: &trimmer::TrimConfig,
    output_dir: Option<&Path>,
    command_line: &str,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    let output_path =
        naming::paired_bam_output_name(input_r1, input_r2, output_dir, cli.basename.as_deref());

    eprintln!("Trimming (paired-end, interleaved uBAM):");
    eprintln!("  R1:     {}", input_r1.display());
    eprintln!("  R2:     {}", input_r2.display());
    eprintln!("  Output: {}", output_path.display());

    // Source header from R1 if it's uBAM; R2's header is ignored (R1's @PG
    // chain is the canonical lineage — R2 is just the mate stream).
    //
    // Internal-invariant backstop (#363). Because the header comes from R1 alone
    // while each side is opened by its own per-file detection below, any BAM
    // reaching this two-file path is wrong output rather than an error: a mixed
    // pair would emit a BAM silently mixing FASTQ-derived records (no aux tags,
    // no source header) with BAM-derived ones, and a two-BAM pair would silently
    // discard R2's @HD/@PG chain and tag dictionary. So the enforced predicate is
    // "no BAM at all", matching the invariant that
    // `reject_bam_format_mismatch_in_pair` in main() actually establishes for
    // this shape — the guard is ~1700 lines away, hence enforcing rather than
    // documenting. (The code this replaced also rejected any BAM; narrowing it
    // to a mismatch check would have half-met the stated intent.)
    let fmt_r1 = detect_input_format(input_r1)?;
    let fmt_r2 = detect_input_format(input_r2)?;
    if matches!(fmt_r1, InputFormat::UnalignedBam) || matches!(fmt_r2, InputFormat::UnalignedBam) {
        anyhow::bail!(
            "internal error: uBAM input reached run_ubam_output_paired_two_files \
             ({} and {}); the paired format guard in main() should have rejected it. \
             Please report this at https://github.com/FelixKrueger/TrimGalore/issues",
            input_r1.display(),
            input_r2.display()
        );
    }
    let source_header = match fmt_r1 {
        InputFormat::UnalignedBam => Some(trim_galore::bam::peek_header(input_r1)?),
        InputFormat::FastqPlain | InputFormat::FastqGz => None,
    };

    let mut reader_r1 = open_sync_reader(input_r1, &cli.preserve_tags)?;
    let mut reader_r2 = open_sync_reader(input_r2, &cli.preserve_tags)?;
    let mut writer = trim_galore::bam::BamWriter::create(
        &output_path,
        source_header.as_ref(),
        &cli.preserve_tags,
        command_line,
        cli.phred_offset(),
    )?;
    let (stats_r1, stats_r2, pair_stats) = trimmer::run_paired_end_to_bam(
        reader_r1.as_mut(),
        reader_r2.as_mut(),
        &mut writer,
        config,
    )?;
    writer.finish()?;

    // Summaries.
    eprintln!("\n=== Summary (Read 1) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r1.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r1.total_reads_with_adapter,
        pct(stats_r1.total_reads_with_adapter, stats_r1.total_reads)
    );
    eprintln!("\n=== Summary (Read 2) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r2.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r2.total_reads_with_adapter,
        pct(stats_r2.total_reads_with_adapter, stats_r2.total_reads)
    );
    eprintln!("\n=== Paired-end validation ===\n");
    eprintln!(
        "Pairs analyzed:                  {:>10}",
        pair_stats.pairs_analyzed
    );
    eprintln!(
        "Pairs removed:                   {:>10} ({:.1}%)",
        pair_stats.pairs_removed,
        pct(pair_stats.pairs_removed, pair_stats.pairs_analyzed)
    );

    if !cli.no_report_file {
        let r1_filename = input_r1
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let r2_filename = input_r2
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let all_input_filenames = vec![r1_filename.clone(), r2_filename.clone()];

        let (r1_txt, r1_json) = naming::paired_report_names(input_r1, input_r1, output_dir);
        let (r2_txt, r2_json) = naming::paired_report_names(input_r2, input_r1, output_dir);
        let r1_desc = PairedReportFile {
            txt_path: r1_txt,
            json_path: r1_json,
            input_filename: r1_filename,
        };
        let r2_desc = PairedReportFile {
            txt_path: r2_txt,
            json_path: r2_json,
            input_filename: r2_filename,
        };

        write_paired_reports(
            cli,
            config,
            false, // uBAM output — gzip flag is FASTQ-output-specific.
            r1_desc,
            r2_desc,
            &all_input_filenames,
            &stats_r1,
            &stats_r2,
            &pair_stats,
            adapters_r1,
            adapters_r2,
            None, // --passthrough rejected upstream
        )?;
    }

    // Run FastQC if requested. fastqc-rust dispatches on file EXTENSION
    // (not content), so this relies on output_path ending in `.bam`.
    // PE uBAM output is a single interleaved BAM, so one call covers both mates.
    if cli.fastqc_requested() {
        fastqc::run(
            &output_path,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
    }

    Ok(())
}

/// uBAM-output single-file (one interleaved BAM in, one interleaved BAM out)
/// paired-end driver. Mirrors `run_paired_ubam_single_file`'s shape for the
/// FASTQ-output path.
#[allow(clippy::too_many_arguments)]
fn run_ubam_output_paired_single_file(
    cli: &Cli,
    input: &Path,
    outputs: &naming::InterleavedBamOutputs,
    config: &trimmer::TrimConfig,
    output_dir: Option<&Path>,
    command_line: &str,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    let output_path = outputs.val.clone();

    eprintln!("Trimming (paired-end, de-interleaved uBAM):");
    eprintln!("  Input:  {}", input.display());
    eprintln!("  Output: {}", output_path.display());

    let source_header = trim_galore::bam::peek_header(input)?;
    let (mut r1, mut r2) = BamReader::open_paired_interleaved_with_tags(input, &cli.preserve_tags)?;
    let mut writer = trim_galore::bam::BamWriter::create(
        &output_path,
        Some(&source_header),
        &cli.preserve_tags,
        command_line,
        cli.phred_offset(),
    )?;
    let (stats_r1, stats_r2, pair_stats) =
        trimmer::run_paired_end_to_bam(&mut r1, &mut r2, &mut writer, config)?;
    writer.finish()?;

    // Summaries (mirror `run_paired_ubam_single_file`).
    eprintln!("\n=== Summary (Read 1) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r1.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r1.total_reads_with_adapter,
        pct(stats_r1.total_reads_with_adapter, stats_r1.total_reads)
    );
    eprintln!("\n=== Summary (Read 2) ===\n");
    eprintln!(
        "Total reads processed:           {:>10}",
        stats_r2.total_reads
    );
    eprintln!(
        "Reads with adapters:             {:>10} ({:.1}%)",
        stats_r2.total_reads_with_adapter,
        pct(stats_r2.total_reads_with_adapter, stats_r2.total_reads)
    );
    eprintln!("\n=== Paired-end validation ===\n");
    eprintln!(
        "Pairs analyzed:                  {:>10}",
        pair_stats.pairs_analyzed
    );
    eprintln!(
        "Pairs removed:                   {:>10} ({:.1}%)",
        pair_stats.pairs_removed,
        pct(pair_stats.pairs_removed, pair_stats.pairs_analyzed)
    );

    if !cli.no_report_file {
        let input_filename = input
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let all_input_filenames = vec![input_filename.clone()];

        let reports = outputs
            .reports
            .as_ref()
            .expect("reports were requested, so the namer built them");
        let r1_desc = PairedReportFile {
            txt_path: reports.r1_txt.clone(),
            json_path: reports.r1_json.clone(),
            input_filename: input_filename.clone(),
        };
        let r2_desc = PairedReportFile {
            txt_path: reports.r2_txt.clone(),
            json_path: reports.r2_json.clone(),
            input_filename: input_filename.clone(),
        };

        write_paired_reports(
            cli,
            config,
            false, // uBAM output — gzip flag is FASTQ-output-specific.
            r1_desc,
            r2_desc,
            &all_input_filenames,
            &stats_r1,
            &stats_r2,
            &pair_stats,
            adapters_r1,
            adapters_r2,
            None,
        )?;
    }

    // Run FastQC if requested. fastqc-rust dispatches on file EXTENSION
    // (not content), so this relies on output_path ending in `.bam`.
    // PE uBAM output is a single interleaved BAM, so one call covers both mates.
    if cli.fastqc_requested() {
        fastqc::run(
            &output_path,
            cli.fastqc_args.as_deref(),
            output_dir,
            cli.cores,
        )?;
    }

    Ok(())
}

/// Per-report-side descriptor for the shared `write_paired_reports` helper.
///
/// The four callers compute their report paths differently — the two-file pairs
/// from the input file paths, the two single-file interleaved arms from a
/// synthesised stem — but the report-writing block downstream is identical.
/// This struct is the seam.
struct PairedReportFile {
    txt_path: std::path::PathBuf,
    json_path: std::path::PathBuf,
    /// Value placed in `report_cfg.input_filename`. For two-file FASTQ, this
    /// is the per-side input file's basename; for single-file uBAM-paired,
    /// the SAME BAM filename is used on both sides.
    input_filename: String,
}

/// Write the standard text + JSON paired-end reports for both R1 and R2.
///
/// One implementation for all four paired report sites. Pure I/O helper — no
/// policy or naming logic; callers supply paths via [`PairedReportFile`].
///
/// `passthrough` is `Some((pt_in, pt_out))` only when `--passthrough` is
/// active on the FASTQ paired path; uBAM paired path always passes `None`
/// (passthrough + uBAM is rejected at `Cli::validate`).
#[allow(clippy::too_many_arguments)]
fn write_paired_reports(
    cli: &Cli,
    config: &trimmer::TrimConfig,
    gzip: bool,
    r1: PairedReportFile,
    r2: PairedReportFile,
    all_input_filenames: &[String],
    stats_r1: &report::TrimStats,
    stats_r2: &report::TrimStats,
    pair_stats: &report::PairValidationStats,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
    passthrough: Option<(&Path, &Path)>,
) -> Result<()> {
    for (idx, (stats, descriptor)) in [(stats_r1, &r1), (stats_r2, &r2)].iter().enumerate() {
        let report_cfg = report::TrimConfig {
            version: env!("CARGO_PKG_VERSION").to_string(),
            quality_cutoff: cli.effective_quality_cutoff(),
            adapters: adapters_r1.to_vec(),
            adapters_r2: adapters_r2.to_vec(),
            times: cli.times,
            error_rate: cli.error_rate,
            stringency: cli.stringency,
            length_cutoff: config.length_cutoff,
            max_length: cli.max_length,
            paired: true,
            gzip,
            trim_n: cli.trim_n,
            nextseq: cli.nextseq.is_some(),
            rrbs: cli.rrbs,
            non_directional: cli.non_directional,
            phred_encoding: cli.phred_offset(),
            poly_a: cli.poly_a,
            poly_g: config.poly_g,
            command_line: std::env::args().collect::<Vec<_>>().join(" "),
            input_filename: descriptor.input_filename.clone(),
            input_filenames: all_input_filenames.to_vec(),
            library: Some(cli.effective_clips()),
        };

        let file = File::create(&descriptor.txt_path)?;
        let mut w = BufWriter::new(file);
        report::write_report_header(&mut w, &report_cfg)?;
        report::write_cutadapt_compatible_section(&mut w, &report_cfg, stats, (idx + 1) as u8)?;
        report::write_run_footer(&mut w, &report_cfg, stats)?;
        // Pair validation stats go in R2 report only (matches Perl behaviour).
        if idx == 1 {
            report::write_pair_validation_stats(&mut w, pair_stats)?;
            // --passthrough is FASTQ-paired-path only; uBAM-paired passes None.
            if let Some((pt_in, pt_out)) = passthrough {
                report::write_passthrough_stats(&mut w, pair_stats, pt_in, pt_out)?;
            }
        }
        eprintln!("Report: {}", descriptor.txt_path.display());

        // JSON report — pair_validation included in BOTH R1 and R2.
        let json_extra = report::JsonReportParams {
            clip_r1: config.clip_r1,
            clip_r2: config.clip_r2,
            three_prime_clip_r1: config.three_prime_clip_r1,
            three_prime_clip_r2: config.three_prime_clip_r2,
            max_n: cli.max_n,
            discard_untrimmed: config.discard_untrimmed,
            consider_already_trimmed: cli.consider_already_trimmed,
        };
        let json_file = File::create(&descriptor.json_path)?;
        let mut jw = BufWriter::new(json_file);
        report::write_json_report(
            &mut jw,
            &report_cfg,
            stats,
            Some(pair_stats),
            (idx + 1) as u8,
            &json_extra,
        )?;
        jw.flush()?;
        eprintln!("JSON report: {}", descriptor.json_path.display());
    }
    Ok(())
}

fn pct(part: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 / total as f64 * 100.0
    }
}

/// Multi-pair driver for the run-and-exit specialty modes (`--clock`,
/// `--implicon`, `--clump_only --paired`). Iterates over `cli.input.chunks(2)`,
/// prints a per-pair header for multi-pair invocations, and runs the supplied
/// per-pair function. Mirrors `--paired`'s output-collision pre-flight (case-
/// insensitive on full path); `pair_outputs` must return EVERY path the pair
/// will write — primaries plus gated secondaries (#391) — so two pairs that
/// would write to the same file fail loudly before any I/O.
fn run_specialty_paired<NameFn, RunFn>(
    cli: &Cli,
    mode_label: &str,
    hint: Option<&str>,
    mut pair_outputs: NameFn,
    mut run_pair: RunFn,
) -> Result<()>
where
    NameFn: FnMut(&Path, &Path) -> Vec<naming::PlannedOutput>,
    RunFn: FnMut(&Path, &Path) -> Result<()>,
{
    // Pre-flight across pairs before any I/O; see io::collision_key for the key.
    let mut planned: Vec<naming::PlannedOutput> = Vec::new();
    for chunk in cli.input.chunks(2) {
        planned.extend(pair_outputs(&chunk[0], &chunk[1]));
    }
    naming::preflight_output_collisions(&planned, &guarded_inputs(cli), hint)?;

    let total_pairs = cli.input.len() / 2;
    for (pair_idx, chunk) in cli.input.chunks(2).enumerate() {
        if total_pairs > 1 {
            eprintln!(
                "\n=== {} pair {} of {} ===",
                mode_label,
                pair_idx + 1,
                total_pairs
            );
        }
        run_pair(&chunk[0], &chunk[1]).with_context(|| {
            format!(
                "processing {} pair {} of {} (R1={}, R2={})",
                mode_label,
                pair_idx + 1,
                total_pairs,
                chunk[0].display(),
                chunk[1].display()
            )
        })?;
    }
    Ok(())
}
