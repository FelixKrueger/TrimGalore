use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use trim_galore::adapter;
use trim_galore::cli::{Cli, rewrite_perl_short_flags};
use trim_galore::clump;
use trim_galore::demux;
use trim_galore::fastq::{FastqReader, FastqWriter};
use trim_galore::fastqc;
use trim_galore::filters::{MaxNFilter, UnpairedLengths};
use trim_galore::io as naming;
use trim_galore::parallel;
use trim_galore::report;
use trim_galore::specialty;
use trim_galore::trimmer;

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

    // Input sanity check on first file
    eprintln!(
        "\nTrim Galore - Oxidized Edition v{}",
        env!("CARGO_PKG_VERSION")
    );
    eprintln!("{}", env!("VERSION_BODY"));
    eprintln!("==================================================\n");

    FastqReader::sanity_check(&cli.input[0])?;

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

    // Specialty modes — bypass normal trimming pipeline entirely
    if let Some(n) = cli.hardtrim5 {
        for input in &cli.input {
            specialty::hardtrim5(
                input,
                n,
                gzip,
                output_dir,
                cli.rename,
                cli.cores,
                cli.gzip_level(),
            )?;
        }
        return Ok(());
    }
    if let Some(n) = cli.hardtrim3 {
        for input in &cli.input {
            specialty::hardtrim3(
                input,
                n,
                gzip,
                output_dir,
                cli.rename,
                cli.cores,
                cli.gzip_level(),
            )?;
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
            |r1, r2| specialty::clock(r1, r2, gzip, output_dir, cli.cores, cli.gzip_level()),
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
                    cli.gzip_level(),
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
                FastqReader::sanity_check(&chunk[0])?;
            }
            FastqReader::sanity_check(&chunk[1])?;

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
                FastqReader::sanity_check(input)?;
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
        gzip_level: cli.gzip_level(),
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
        let clump_layout = resolve_clump_layout(cli)?;
        parallel::run_single_end_parallel(
            input,
            &output_path,
            config,
            cli.cores,
            gzip,
            clump_layout,
        )?
    } else {
        let mut reader = FastqReader::open(input)?;
        let mut writer = FastqWriter::create(&output_path, gzip, 1, config.gzip_level)?;
        let stats = trimmer::run_single_end(&mut reader, &mut writer, config)?;
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
            cli.gzip_level(),
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

    let (stats_r1, stats_r2, pair_stats) = if cli.cores > 1 || cli.clumpify {
        // Worker-pool parallel path: N workers each handle trim + compress.
        // `--clumpy` always routes here (validation enforces cores >= 2).
        let clump_layout = resolve_clump_layout(cli)?;
        parallel::run_paired_end_parallel(
            input_r1,
            input_r2,
            &output_r1,
            &output_r2,
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
        // Sequential path (--cores 1)
        let mut reader_r1 = FastqReader::open(input_r1)?;
        let mut reader_r2 = FastqReader::open(input_r2)?;
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
            &mut reader_r1,
            &mut reader_r2,
            &mut writer_r1,
            &mut writer_r2,
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
        drop(writer_r1);
        drop(writer_r2);
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

    // Write reports
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

        for (idx, (input, stats)) in [(input_r1, &stats_r1), (input_r2, &stats_r2)]
            .iter()
            .enumerate()
        {
            let report_path = naming::report_name(input, output_dir);
            let input_filename = input
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
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
                input_filename,
                input_filenames: all_input_filenames.clone(),
            };

            let file = File::create(&report_path)?;
            let mut w = BufWriter::new(file);
            report::write_report_header(&mut w, &report_cfg)?;
            report::write_cutadapt_compatible_section(&mut w, &report_cfg, stats, (idx + 1) as u8)?;
            report::write_run_footer(&mut w, &report_cfg, stats)?;
            // Pair validation stats go in R2 report only (matches Perl behavior)
            if idx == 1 {
                report::write_pair_validation_stats(&mut w, &pair_stats)?;
            }
            eprintln!("Report: {}", report_path.display());

            // JSON report — pair_validation included in BOTH R1 and R2
            // (use effective config values, not raw CLI, for clip/discard params)
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
            report::write_json_report(
                &mut jw,
                &report_cfg,
                stats,
                Some(&pair_stats),
                (idx + 1) as u8,
                &json_extra,
            )?;
            jw.flush()?;
            eprintln!("JSON report: {}", json_path.display());
        }
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
