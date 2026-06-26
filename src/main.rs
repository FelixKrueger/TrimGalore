use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use trim_galore::adapter;
use trim_galore::bam::BamReader;
use trim_galore::cli::{Cli, rewrite_perl_short_flags};
use trim_galore::clump;
use trim_galore::demux;
use trim_galore::fastq::{FastqReader, FastqWriter, RecordSource};
use trim_galore::fastqc;
use trim_galore::filters::{MaxNFilter, UnpairedLengths};
use trim_galore::format::{InputFormat, detect_input_format};
use trim_galore::io as naming;
use trim_galore::parallel;
use trim_galore::report;
use trim_galore::specialty;
use trim_galore::trimmer;

/// Format-aware sanity check — dispatches `FastqReader::sanity_check` for
/// FASTQ input or peek-reads the first BAM record (asserting `is_unmapped()`)
/// for uBAM. The per-record aligned-BAM check in `BamReader::next_record`
/// catches mixed-aligned BAMs that slip past this fast-path.
fn sanity_check_any(path: &std::path::Path) -> Result<()> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => FastqReader::sanity_check(path),
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

/// Open a threaded reader for `path`, dispatching by detected format.
/// `preserve_tags` is honoured for BAM input; silently ignored for FASTQ.
fn open_threaded_reader(
    path: &std::path::Path,
    preserve_tags: &[String],
) -> Result<Box<dyn RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => {
            Ok(Box::new(FastqReader::open_threaded(path)?))
        }
        InputFormat::UnalignedBam => Ok(Box::new(BamReader::open_threaded_with_tags(
            path,
            preserve_tags,
        )?)),
    }
}

/// Open a sync (single-threaded) reader for `path`, dispatching by detected
/// format. Used by the `--cores 1` (serial) path.
fn open_sync_reader(
    path: &std::path::Path,
    preserve_tags: &[String],
) -> Result<Box<dyn RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => Ok(Box::new(FastqReader::open(path)?)),
        InputFormat::UnalignedBam => Ok(Box::new(
            BamReader::open(path)?.with_preserved_tags(preserve_tags),
        )),
    }
}

type AdapterList = Vec<(String, String)>;
type SetupResult = Result<(String, AdapterList, AdapterList, trimmer::TrimConfig)>;
type ResolvedAdapter = Result<(String, AdapterList, AdapterList, Option<(usize, usize)>)>;

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
    let min_required = clump::clumpify_min_memory_bytes(cli.cores);
    if memory_bytes < min_required {
        eprintln!();
        eprintln!(
            "WARNING: --memory {} is too small for --clumpify at --cores {} \
             (need ≥ {} MiB).",
            cli.memory,
            cli.cores,
            min_required / (1024 * 1024),
        );
        eprintln!(
            "         Falling back to plain mode (no read reordering). \
             Increase --memory or drop --clumpify to silence this warning."
        );
        eprintln!();
        return Ok(None);
    }
    let layout = clump::resolve_layout(memory_bytes, cli.cores)?;
    eprintln!(
        "clumpify: {} bins × {} MB; predicted peak ≈ {} MB (gzip level {})",
        layout.n_bins,
        layout.bin_byte_budget / (1024 * 1024),
        layout.predicted_peak_bytes(cli.cores) / (1024 * 1024),
        cli.compression,
    );
    Ok(Some(layout))
}

fn main() -> Result<()> {
    env_logger::init();

    // Pre-parse rewrite for Perl-era `-r1`/`-r2` short flags before clap sees them.
    let cli = Cli::parse_from(rewrite_perl_short_flags(std::env::args()));
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
    if cli.paired && cli.input.len() == 1 && !matches!(input_formats[0], InputFormat::UnalignedBam)
    {
        anyhow::bail!(
            "--paired with a single input file is only legal if that file is a uBAM. \
             Provide R1 and R2 as separate files for FASTQ paired-end mode."
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
                 --cores {} is ignored. For high-throughput uBAM output, run \
                 multiple invocations in parallel.",
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
    }

    // Specialty modes — bypass normal trimming pipeline entirely
    if let Some(n) = cli.hardtrim5 {
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
                )?,
            }
        }
        return Ok(());
    }
    if let Some(n) = cli.hardtrim3 {
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
                )?,
            }
        }
        return Ok(());
    }
    if cli.clock {
        run_specialty_paired(
            &cli,
            "Clock",
            |r1, r2| {
                (
                    specialty::clock_output_name(r1, "R1", output_dir, gzip),
                    specialty::clock_output_name(r2, "R2", output_dir, gzip),
                )
            },
            |r1, r2| specialty::clock(r1, r2, gzip, output_dir, cli.cores, cli.compression),
        )?;
        return Ok(());
    }
    if let Some(umi_len) = cli.implicon {
        run_specialty_paired(
            &cli,
            "IMPLICON",
            |r1, r2| {
                (
                    specialty::implicon_output_name(r1, umi_len, "R1", output_dir, gzip),
                    specialty::implicon_output_name(r2, umi_len, "R2", output_dir, gzip),
                )
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
        let (_label, adapters_r1, adapters_r2, config) = setup_trimming(&cli, &cli.input[0])?;
        run_paired_ubam_single_file(
            &cli,
            &cli.input[0],
            &config,
            gzip,
            output_dir,
            &adapters_r1,
            &adapters_r2,
        )?;
        return Ok(());
    }

    if cli.paired {
        // Pre-flight: detect output-path collisions across pairs before any I/O.
        // Hash key is the full path, case-folded (ASCII lowercase) so that paths
        // differing only in letter-case — which alias the same file on APFS/NTFS
        // — collide here rather than silently overwriting on disk (issue #216).
        // Pragmatic trade-off: on opt-in case-sensitive APFS volumes this may
        // false-positive, but the penalty is a loud early error rather than
        // silent data loss.
        let mut out_paths: std::collections::HashMap<String, std::path::PathBuf> =
            std::collections::HashMap::new();
        let norm = |p: &std::path::Path| -> String { p.to_string_lossy().to_ascii_lowercase() };
        for chunk in cli.input.chunks(2) {
            let (o1, o2) = naming::paired_end_output_names(
                &chunk[0],
                &chunk[1],
                output_dir,
                cli.basename.as_deref(),
                gzip,
            );
            let mut candidates = vec![o1, o2];
            if cli.retain_unpaired {
                let (u1, u2) = naming::unpaired_output_names(
                    &chunk[0],
                    &chunk[1],
                    output_dir,
                    cli.basename.as_deref(),
                    gzip,
                );
                candidates.push(u1);
                candidates.push(u2);
            }
            // --passthrough adds a third output path per pair. v1 only
            // supports a single pair (Cli::validate enforces input.len() == 2
            // when passthrough is set), so this either contributes zero or
            // one extra candidate to the collision set.
            if let Some(ref pt_input) = cli.passthrough {
                candidates.push(naming::passthrough_output_name(
                    pt_input,
                    output_dir,
                    cli.basename.as_deref(),
                    gzip,
                ));
            }
            for p in candidates {
                if let Some(existing) = out_paths.insert(norm(&p), p.clone()) {
                    anyhow::bail!(
                        "Output path collision (case-insensitive, for APFS/NTFS safety): \
                         {} and {} would be written to the same file. \
                         Check that input pairs produce distinct output paths \
                         (e.g., different source directories or `--output-dir`).",
                        existing.display(),
                        p.display()
                    );
                }
            }
        }

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
            if pair_idx == 0
                && let Some(ref pt_path) = cli.passthrough
            {
                FastqReader::sanity_check(pt_path)?;
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
    let (adapter_label, adapters_r1, adapters_r2, autodetect_poly_g) =
        resolve_adapter(cli, input_file)?;

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
    if !adapters_r2.is_empty() {
        if adapters_r2.len() == 1 {
            eprintln!("Adapter 2 (Read 2): {}", adapters_r2[0].1);
        } else {
            eprintln!("Adapters R2 ({} sequences):", adapters_r2.len());
            for (name, seq) in &adapters_r2 {
                eprintln!("  {}: {}", name, seq);
            }
        }
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

    // RRBS: auto-set --clip_r2 2 for directional paired-end mode
    let clip_r2 = if cli.rrbs && !cli.non_directional && cli.paired && cli.clip_r2.is_none() {
        eprintln!(
            "Setting the option '--clip_r2 2' (to remove methylation bias from the start of Read 2)"
        );
        Some(2)
    } else {
        cli.clip_r2
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
        clip_r1: cli.clip_r1,
        clip_r2,
        three_prime_clip_r1: cli.three_prime_clip_r1,
        three_prime_clip_r2: cli.three_prime_clip_r2,
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

/// Returns (adapter_label, adapters_r1, adapters_r2, poly_g_from_autodetect).
/// `adapter_label` is for display purposes (e.g., "Illumina", "user-specified").
/// `adapters_r1`/`adapters_r2` are `(name, sequence)` pairs; r2 is empty if not set.
/// The last element is Some((poly_g_count, reads_scanned)) when auto-detection ran,
/// None when the adapter was user-specified or preset-selected (poly-G must be
/// detected separately via `adapter::detect_poly_g()`).
fn resolve_adapter(cli: &Cli, input_file: &Path) -> ResolvedAdapter {
    if !cli.adapter.is_empty() {
        let adapters_r1 = adapter::parse_adapter_specs(&cli.adapter)?;
        let adapters_r2 = adapter::parse_adapter_specs(&cli.adapter2)?;
        let label = "user-specified".to_string();
        return Ok((label, adapters_r1, adapters_r2, None));
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
        writer.flush()?;
        drop(writer);
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
    if cli.fastqc || cli.fastqc_args.is_some() {
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
    let passthrough_output: Option<std::path::PathBuf> =
        passthrough_input.map(|pt| naming::passthrough_output_name(pt, output_dir, basename, gzip));
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

    // Two-file paired-uBAM is rejected per PLAN §3.3 — paired uBAM expects
    // a single interleaved file (handled in main() before this fn is called).
    for p in [input_r1, input_r2] {
        if matches!(detect_input_format(p)?, InputFormat::UnalignedBam) {
            anyhow::bail!(
                "--paired with two BAM files is not supported. \
                 uBAM paired mode expects a single interleaved file: \
                 `trim_galore --paired interleaved.bam`. \
                 Got two BAM files; one of them is {}.",
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
        // Sequential path (--cores 1) — paired-BAM rejected above, so both
        // are FastqReader.
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

        writer_r1.flush()?;
        writer_r2.flush()?;
        if let Some(ref mut w) = writer_passthrough {
            w.flush()?;
        }
        if let Some(ref mut w) = unpaired_w1 {
            w.flush()?;
        }
        if let Some(ref mut w) = unpaired_w2 {
            w.flush()?;
        }
        drop(writer_r1);
        drop(writer_r2);
        drop(writer_passthrough);
        drop(unpaired_w1);
        drop(unpaired_w2);

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
            "R2 reads with poly-C trimmed:    {:>10} ({:.1}%)",
            stats_r2.poly_g_trimmed,
            pct(stats_r2.poly_g_trimmed, stats_r2.total_reads)
        );
        eprintln!(
            "  R2 poly-C bases removed:       {:>10}",
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

        let r1 = PairedReportFile {
            txt_path: naming::report_name(input_r1, output_dir),
            json_path: naming::json_report_name(input_r1, output_dir),
            input_filename: all_input_filenames[0].clone(),
        };
        let r2 = PairedReportFile {
            txt_path: naming::report_name(input_r2, output_dir),
            json_path: naming::json_report_name(input_r2, output_dir),
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
    if cli.fastqc || cli.fastqc_args.is_some() {
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
/// Output names: derived from the input BAM's stem (`input.bam` →
/// `input_val_1.fq[.gz]` + `input_val_2.fq[.gz]`).
#[allow(clippy::too_many_arguments)]
fn run_paired_ubam_single_file(
    cli: &Cli,
    input: &Path,
    config: &trimmer::TrimConfig,
    gzip: bool,
    output_dir: Option<&Path>,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    // Derive output paths from the single BAM input's stem.
    let stem = input
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let ext = if gzip { ".fq.gz" } else { ".fq" };
    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| input.parent().unwrap_or(Path::new(".")).to_path_buf());
    let output_r1 = dir.join(format!("{}_val_1{}", stem, ext));
    let output_r2 = dir.join(format!("{}_val_2{}", stem, ext));

    eprintln!("Trimming (paired-end, de-interleaved uBAM):");
    eprintln!("  Input:     {}", input.display());
    eprintln!("  Output R1: {}", output_r1.display());
    eprintln!("  Output R2: {}", output_r2.display());

    let (unpaired_r1_path, unpaired_r2_path) = if cli.retain_unpaired {
        let up1 = dir.join(format!("{}_unpaired_1{}", stem, ext));
        let up2 = dir.join(format!("{}_unpaired_2{}", stem, ext));
        eprintln!("  Unpaired R1: {}", up1.display());
        eprintln!("  Unpaired R2: {}", up2.display());
        (Some(up1), Some(up2))
    } else {
        (None, None)
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
        writer_r1.flush()?;
        writer_r2.flush()?;
        if let Some(ref mut w) = unpaired_w1 {
            w.flush()?;
        }
        if let Some(ref mut w) = unpaired_w2 {
            w.flush()?;
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

        let r1 = PairedReportFile {
            txt_path: dir.join(format!("{}_R1_trimming_report.txt", stem)),
            json_path: dir.join(format!("{}_R1_trimming_report.json", stem)),
            input_filename: input_filename.clone(),
        };
        let r2 = PairedReportFile {
            txt_path: dir.join(format!("{}_R2_trimming_report.txt", stem)),
            json_path: dir.join(format!("{}_R2_trimming_report.json", stem)),
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
        let (_label, adapters_r1, adapters_r2, config) = setup_trimming(cli, &cli.input[0])?;
        return run_ubam_output_paired_single_file(
            cli,
            &cli.input[0],
            &config,
            output_dir,
            command_line,
            &adapters_r1,
            &adapters_r2,
        );
    }

    if cli.paired {
        // Two-file paired-BAM input is rejected per PLAN §3.3 — uBAM paired
        // mode expects a single interleaved file (handled in the
        // input.len() == 1 branch above). Mirror the FASTQ-path rejection
        // and fail up-front BEFORE the collision pre-flight runs.
        // Code-review round-2 B-NIT-1 consolidation.
        for chunk in cli.input.chunks(2) {
            for p in [&chunk[0], &chunk[1]] {
                if matches!(detect_input_format(p)?, InputFormat::UnalignedBam) {
                    anyhow::bail!(
                        "--paired with two BAM files is not supported. \
                         uBAM paired mode expects a single interleaved file: \
                         `trim_galore --paired interleaved.bam`. \
                         Got two BAM files; one of them is {}.",
                        p.display()
                    );
                }
            }
        }
        // Pre-flight: one BAM output per pair; collision on case-folded path.
        let mut out_paths: std::collections::HashMap<String, std::path::PathBuf> =
            std::collections::HashMap::new();
        let norm = |p: &std::path::Path| -> String { p.to_string_lossy().to_ascii_lowercase() };
        for chunk in cli.input.chunks(2) {
            let out = naming::paired_bam_output_name(
                &chunk[0],
                &chunk[1],
                output_dir,
                cli.basename.as_deref(),
            );
            if let Some(existing) = out_paths.insert(norm(&out), out.clone()) {
                anyhow::bail!(
                    "Output path collision (case-insensitive, for APFS/NTFS safety): \
                     {} and {} would be written to the same file. \
                     Check that input pairs produce distinct output paths \
                     (e.g., different source directories or `--output-dir`).",
                    existing.display(),
                    out.display()
                );
            }
        }

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

    // Two-file paired-BAM rejection lives up-front in `run_ubam_output`
    // (code-review round-2 B-NIT-1) — both inputs are guaranteed FASTQ
    // here.

    eprintln!("Trimming (paired-end, interleaved uBAM):");
    eprintln!("  R1:     {}", input_r1.display());
    eprintln!("  R2:     {}", input_r2.display());
    eprintln!("  Output: {}", output_path.display());

    // Source header from R1 if it's uBAM; R2's header is ignored (R1's @PG
    // chain is the canonical lineage — R2 is just the mate stream).
    // (At this point both inputs are FASTQ; the BAM-input case is handled
    // separately by run_ubam_output_paired_single_file.)
    let source_header = match detect_input_format(input_r1)? {
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

        let r1_desc = PairedReportFile {
            txt_path: naming::report_name(input_r1, output_dir),
            json_path: naming::json_report_name(input_r1, output_dir),
            input_filename: r1_filename,
        };
        let r2_desc = PairedReportFile {
            txt_path: naming::report_name(input_r2, output_dir),
            json_path: naming::json_report_name(input_r2, output_dir),
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

    Ok(())
}

/// uBAM-output single-file (one interleaved BAM in, one interleaved BAM out)
/// paired-end driver. Mirrors `run_paired_ubam_single_file`'s shape for the
/// FASTQ-output path.
fn run_ubam_output_paired_single_file(
    cli: &Cli,
    input: &Path,
    config: &trimmer::TrimConfig,
    output_dir: Option<&Path>,
    command_line: &str,
    adapters_r1: &[(String, String)],
    adapters_r2: &[(String, String)],
) -> Result<()> {
    // Output stem derived from the BAM input's filename, same as the
    // FASTQ-output single-file paired path.
    let stem = input
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| input.parent().unwrap_or(Path::new(".")).to_path_buf());
    let output_path = dir.join(format!("{}_val.bam", stem));

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

        let r1_desc = PairedReportFile {
            txt_path: dir.join(format!("{}_R1_trimming_report.txt", stem)),
            json_path: dir.join(format!("{}_R1_trimming_report.json", stem)),
            input_filename: input_filename.clone(),
        };
        let r2_desc = PairedReportFile {
            txt_path: dir.join(format!("{}_R2_trimming_report.txt", stem)),
            json_path: dir.join(format!("{}_R2_trimming_report.json", stem)),
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

    Ok(())
}

/// Per-report-side descriptor for the shared `write_paired_reports` helper.
///
/// The two callers (`run_paired` for two-file FASTQ pairs, and
/// `run_paired_ubam_single_file` for one-file interleaved uBAM) compute their
/// report paths differently — one from the input file paths, the other from
/// a synthesised stem — but the report-writing block downstream is identical.
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
/// Centralises the ~50-line duplicated block previously inlined in both
/// `run_paired` and `run_paired_ubam_single_file`. Pure I/O helper — no
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
/// `--implicon`). Iterates over `cli.input.chunks(2)`, prints a per-pair
/// header for multi-pair invocations, and runs the supplied per-pair
/// function. Mirrors `--paired`'s output-collision pre-flight (case-
/// insensitive on full path) so two pairs that would write to the same
/// output file fail loudly before any I/O.
fn run_specialty_paired<NameFn, RunFn>(
    cli: &Cli,
    mode_label: &str,
    mut output_names: NameFn,
    mut run_pair: RunFn,
) -> Result<()>
where
    NameFn: FnMut(&Path, &Path) -> (std::path::PathBuf, std::path::PathBuf),
    RunFn: FnMut(&Path, &Path) -> Result<()>,
{
    // Pre-flight: detect output-path collisions across pairs (same shape
    // as the --paired pre-flight at lines 82–125).
    let mut out_paths: std::collections::HashMap<String, std::path::PathBuf> =
        std::collections::HashMap::new();
    let norm = |p: &std::path::Path| -> String { p.to_string_lossy().to_ascii_lowercase() };
    for chunk in cli.input.chunks(2) {
        let (o1, o2) = output_names(&chunk[0], &chunk[1]);
        for p in [o1, o2] {
            if let Some(existing) = out_paths.insert(norm(&p), p.clone()) {
                anyhow::bail!(
                    "Output path collision (case-insensitive, for APFS/NTFS safety): \
                     {} and {} would be written to the same file. \
                     Check that input pairs produce distinct output paths \
                     (e.g., different source directories or `--output_dir`).",
                    existing.display(),
                    p.display()
                );
            }
        }
    }

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
