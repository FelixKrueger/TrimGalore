use anyhow::Result;
use clap::Parser;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use optimus_prime::adapter;
use optimus_prime::cli::Cli;
use optimus_prime::fastq::{FastqReader, FastqWriter};
use optimus_prime::filters::MaxNFilter;
use optimus_prime::io as naming;
use optimus_prime::report;
use optimus_prime::specialty;
use optimus_prime::trimmer;

fn main() -> Result<()> {
    env_logger::init();

    let cli = Cli::parse();
    cli.validate()?;

    // Input sanity check on first file
    eprintln!("\nOptimus Prime v{}", env!("CARGO_PKG_VERSION"));
    eprintln!("=========================================\n");

    FastqReader::sanity_check(&cli.input[0])?;

    let gzip = !cli.dont_gzip;
    let output_dir = cli.output_dir.as_deref();

    // Specialty modes — bypass normal trimming pipeline entirely
    if let Some(n) = cli.hardtrim5 {
        for input in &cli.input {
            specialty::hardtrim5(input, n, gzip, output_dir, cli.rename, cli.cores)?;
        }
        return Ok(());
    }
    if let Some(n) = cli.hardtrim3 {
        for input in &cli.input {
            specialty::hardtrim3(input, n, gzip, output_dir, cli.rename, cli.cores)?;
        }
        return Ok(());
    }
    if cli.clock {
        specialty::clock(&cli.input[0], &cli.input[1], gzip, output_dir, cli.cores)?;
        return Ok(());
    }
    if let Some(umi_len) = cli.implicon {
        specialty::implicon(&cli.input[0], &cli.input[1], umi_len, gzip, output_dir, cli.cores)?;
        return Ok(());
    }

    // Determine adapter
    let (adapter_name, adapter_seq, adapter_r2_seq) = resolve_adapter(&cli)?;

    eprintln!("Adapter: {} ({})", adapter_name, adapter_seq);
    if let Some(ref a2) = adapter_r2_seq {
        eprintln!("Adapter 2 (Read 2): {}", a2);
    }

    // Resolve length cutoff: smallRNA adapter auto-reduces to 18bp
    // (matches TrimGalore behavior to preserve short miRNA species)
    let length_cutoff = cli.length.unwrap_or_else(|| {
        if adapter_seq == "TGGAATTCTCGG" {
            eprintln!("Reducing length cutoff to 18bp for small RNA-Seq reads because a cutoff of 20bp may remove some short species of small RNAs if they had been trimmed by 1,2 or 3bp");
            18
        } else {
            20
        }
    });
    eprintln!();

    // Build max_n filter
    let max_n = cli.max_n.map(|v| {
        if v >= 1.0 {
            MaxNFilter::Count(v as usize)
        } else {
            MaxNFilter::Fraction(v)
        }
    });

    // RRBS: auto-set --clip_r2 2 for directional paired-end mode
    // (removes filled-in C bases from Read 2 5' end at MspI sites)
    let clip_r2 = if cli.rrbs && !cli.non_directional && cli.paired && cli.clip_r2.is_none() {
        eprintln!("Setting the option '--clip_r2 2' (to remove methylation bias from the start of Read 2)");
        Some(2)
    } else {
        cli.clip_r2
    };

    if cli.rrbs {
        eprintln!("File was specified to be an MspI-digested RRBS sample. Read 1 sequences with adapter contamination will be trimmed a further 2 bp from their 3' end, and Read 2 sequences will be trimmed by 2 bp from their 5' end to remove potential methylation-biased bases from the end-repair reaction");
    }
    if cli.non_directional {
        eprintln!("File was specified to be a non-directional MspI-digested RRBS sample. Sequences starting with either 'CAA' or 'CGA' will have the first 2 bp trimmed off to remove potential methylation-biased bases from the end-repair reaction");
    }

    // Build trim config
    let config = trimmer::TrimConfig {
        adapter: adapter_seq.as_bytes().to_vec(),
        adapter_r2: adapter_r2_seq.as_ref().map(|s| s.as_bytes().to_vec()),
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
    };

    let basename = cli.basename.as_deref();

    if cli.cores > 1 && gzip {
        eprintln!("Using {} compression threads (parallel gzip)", cli.cores);
    }

    if cli.paired {
        run_paired(&cli, &config, gzip, output_dir, basename, &adapter_seq, &adapter_r2_seq)?;
    } else {
        run_single(&cli, &config, gzip, output_dir, basename, &adapter_seq)?;
    }

    Ok(())
}

fn resolve_adapter(cli: &Cli) -> Result<(String, String, Option<String>)> {
    if let Some(ref seq) = cli.adapter {
        let name = "user-specified".to_string();
        let r2 = cli.adapter2.clone();
        return Ok((name, seq.clone(), r2));
    }

    if cli.nextera {
        return Ok((adapter::NEXTERA.name.to_string(), adapter::NEXTERA.seq.to_string(), None));
    }
    if cli.small_rna {
        return Ok((
            adapter::SMALL_RNA.name.to_string(),
            adapter::SMALL_RNA.seq.to_string(),
            adapter::SMALL_RNA.seq_r2.map(|s| s.to_string()),
        ));
    }
    if cli.stranded_illumina {
        return Ok((
            adapter::STRANDED_ILLUMINA.name.to_string(),
            adapter::STRANDED_ILLUMINA.seq.to_string(),
            None,
        ));
    }
    if cli.bgiseq {
        return Ok((
            adapter::BGISEQ.name.to_string(),
            adapter::BGISEQ.seq.to_string(),
            adapter::BGISEQ.seq_r2.map(|s| s.to_string()),
        ));
    }
    if cli.illumina {
        return Ok((adapter::ILLUMINA.name.to_string(), adapter::ILLUMINA.seq.to_string(), None));
    }

    // Auto-detect
    eprintln!("Auto-detecting adapter type...");
    let detection = adapter::autodetect_adapter(&cli.input[0], cli.consider_already_trimmed)?;
    eprintln!("{}", detection.message);

    let r2 = detection.adapter.seq_r2.map(|s| s.to_string());
    Ok((
        detection.adapter.name.to_string(),
        detection.adapter.seq.to_string(),
        r2,
    ))
}

fn run_single(
    cli: &Cli,
    config: &trimmer::TrimConfig,
    gzip: bool,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    adapter_seq: &str,
) -> Result<()> {
    let input = &cli.input[0];
    let output_path = naming::single_end_output_name(input, output_dir, basename, gzip);
    let report_path = naming::report_name(input, output_dir);

    eprintln!("Trimming: {}", input.display());
    eprintln!("Output:   {}", output_path.display());

    let mut reader = if cli.cores > 1 {
        FastqReader::open_threaded(input)?
    } else {
        FastqReader::open(input)?
    };
    let mut writer = FastqWriter::create(&output_path, gzip, cli.cores)?;

    let stats = trimmer::run_single_end(&mut reader, &mut writer, config)?;
    writer.flush()?;
    drop(writer); // Finalize gzip stream (writes EOF trailer) before any downstream processing

    // Print summary
    eprintln!("\n=== Summary ===\n");
    eprintln!("Total reads processed:           {:>10}", stats.total_reads);
    eprintln!("Reads with adapters:             {:>10} ({:.1}%)",
        stats.reads_with_adapter,
        pct(stats.reads_with_adapter, stats.total_reads));
    eprintln!("Reads too short:                 {:>10} ({:.1}%)",
        stats.too_short, pct(stats.too_short, stats.total_reads));
    eprintln!("Reads written (passing filters): {:>10} ({:.1}%)",
        stats.reads_written, pct(stats.reads_written, stats.total_reads));
    if stats.rrbs_trimmed_3prime > 0 {
        eprintln!("RRBS trimmed (3' end, adapter): {:>10} ({:.1}%)",
            stats.rrbs_trimmed_3prime, pct(stats.rrbs_trimmed_3prime, stats.total_reads));
    }
    if stats.rrbs_trimmed_5prime > 0 {
        eprintln!("RRBS trimmed (5' end, CAA/CGA): {:>10} ({:.1}%)",
            stats.rrbs_trimmed_5prime, pct(stats.rrbs_trimmed_5prime, stats.total_reads));
    }
    if stats.poly_a_trimmed > 0 {
        eprintln!("Reads with poly-A/T trimmed:     {:>10} ({:.1}%)",
            stats.poly_a_trimmed, pct(stats.poly_a_trimmed, stats.total_reads));
        eprintln!("  Poly-A/T bases removed:        {:>10}", stats.poly_a_bases_trimmed);
    }

    // Write report
    if !cli.no_report_file {
        let report_cfg = report::TrimConfig {
            version: env!("CARGO_PKG_VERSION").to_string(),
            quality_cutoff: cli.effective_quality_cutoff(),
            adapter: adapter_seq.to_string(),
            adapter_r2: cli.adapter2.clone(),
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
            command_line: std::env::args().collect::<Vec<_>>().join(" "),
        };

        let file = File::create(&report_path)?;
        let mut w = BufWriter::new(file);
        report::write_report_header(&mut w, &report_cfg)?;
        report::write_run_stats(&mut w, &stats)?;
        eprintln!("\nReport: {}", report_path.display());
    }

    // Run FastQC if requested
    if cli.fastqc || cli.fastqc_args.is_some() {
        run_fastqc(&output_path, cli.fastqc_args.as_deref(), output_dir)?;
    }

    Ok(())
}

fn run_paired(
    cli: &Cli,
    config: &trimmer::TrimConfig,
    gzip: bool,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    adapter_seq: &str,
    adapter_r2_seq: &Option<String>,
) -> Result<()> {
    let input_r1 = &cli.input[0];
    let input_r2 = &cli.input[1];

    // Sanity check R2 as well
    FastqReader::sanity_check(input_r2)?;

    let (output_r1, output_r2) =
        naming::paired_end_output_names(input_r1, input_r2, output_dir, basename, gzip);

    eprintln!("Trimming (paired-end, single-pass):");
    eprintln!("  R1: {}", input_r1.display());
    eprintln!("  R2: {}", input_r2.display());
    eprintln!("  Output R1: {}", output_r1.display());
    eprintln!("  Output R2: {}", output_r2.display());

    let mut reader_r1 = if cli.cores > 1 {
        FastqReader::open_threaded(input_r1)?
    } else {
        FastqReader::open(input_r1)?
    };
    let mut reader_r2 = if cli.cores > 1 {
        FastqReader::open_threaded(input_r2)?
    } else {
        FastqReader::open(input_r2)?
    };
    let mut writer_r1 = FastqWriter::create(&output_r1, gzip, cli.cores)?;
    let mut writer_r2 = FastqWriter::create(&output_r2, gzip, cli.cores)?;

    // Optional unpaired writers
    let (mut unpaired_w1, mut unpaired_w2) = if cli.retain_unpaired {
        let (up1, up2) = naming::unpaired_output_names(input_r1, input_r2, output_dir, basename, gzip);
        eprintln!("  Unpaired R1: {}", up1.display());
        eprintln!("  Unpaired R2: {}", up2.display());
        (
            Some(FastqWriter::create(&up1, gzip, cli.cores)?),
            Some(FastqWriter::create(&up2, gzip, cli.cores)?),
        )
    } else {
        (None, None)
    };

    let (stats_r1, stats_r2, pair_stats) = trimmer::run_paired_end(
        &mut reader_r1,
        &mut reader_r2,
        &mut writer_r1,
        &mut writer_r2,
        unpaired_w1.as_mut(),
        unpaired_w2.as_mut(),
        config,
        cli.length_1,
        cli.length_2,
    )?;

    writer_r1.flush()?;
    writer_r2.flush()?;
    if let Some(ref mut w) = unpaired_w1 { w.flush()?; }
    if let Some(ref mut w) = unpaired_w2 { w.flush()?; }
    drop(writer_r1); // Finalize gzip streams before downstream processing
    drop(writer_r2);
    drop(unpaired_w1);
    drop(unpaired_w2);

    // Print summary
    eprintln!("\n=== Summary (Read 1) ===\n");
    eprintln!("Total reads processed:           {:>10}", stats_r1.total_reads);
    eprintln!("Reads with adapters:             {:>10} ({:.1}%)",
        stats_r1.reads_with_adapter,
        pct(stats_r1.reads_with_adapter, stats_r1.total_reads));

    eprintln!("\n=== Summary (Read 2) ===\n");
    eprintln!("Total reads processed:           {:>10}", stats_r2.total_reads);
    eprintln!("Reads with adapters:             {:>10} ({:.1}%)",
        stats_r2.reads_with_adapter,
        pct(stats_r2.reads_with_adapter, stats_r2.total_reads));

    if stats_r1.poly_a_trimmed > 0 {
        eprintln!("R1 reads with poly-A trimmed:    {:>10} ({:.1}%)",
            stats_r1.poly_a_trimmed, pct(stats_r1.poly_a_trimmed, stats_r1.total_reads));
    }
    if stats_r2.poly_a_trimmed > 0 {
        eprintln!("R2 reads with poly-T trimmed:    {:>10} ({:.1}%)",
            stats_r2.poly_a_trimmed, pct(stats_r2.poly_a_trimmed, stats_r2.total_reads));
    }

    eprintln!("\n=== Paired-end validation ===\n");
    eprintln!("Pairs analyzed:                  {:>10}", pair_stats.pairs_analyzed);
    eprintln!("Pairs removed:                   {:>10} ({:.1}%)",
        pair_stats.pairs_removed,
        pct(pair_stats.pairs_removed, pair_stats.pairs_analyzed));
    if pair_stats.r1_unpaired > 0 || pair_stats.r2_unpaired > 0 {
        eprintln!("Unpaired R1 kept:                {:>10}", pair_stats.r1_unpaired);
        eprintln!("Unpaired R2 kept:                {:>10}", pair_stats.r2_unpaired);
    }

    // Write reports
    if !cli.no_report_file {
        for (input, stats) in [(input_r1, &stats_r1), (input_r2, &stats_r2)] {
            let report_path = naming::report_name(input, output_dir);
            let report_cfg = report::TrimConfig {
                version: env!("CARGO_PKG_VERSION").to_string(),
                quality_cutoff: cli.effective_quality_cutoff(),
                adapter: adapter_seq.to_string(),
                adapter_r2: adapter_r2_seq.clone(),
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
                command_line: std::env::args().collect::<Vec<_>>().join(" "),
            };

            let file = File::create(&report_path)?;
            let mut w = BufWriter::new(file);
            report::write_report_header(&mut w, &report_cfg)?;
            report::write_run_stats(&mut w, stats)?;
            report::write_pair_validation_stats(&mut w, &pair_stats)?;
            eprintln!("Report: {}", report_path.display());
        }
    }

    // Run FastQC if requested
    if cli.fastqc || cli.fastqc_args.is_some() {
        run_fastqc(&output_r1, cli.fastqc_args.as_deref(), output_dir)?;
        run_fastqc(&output_r2, cli.fastqc_args.as_deref(), output_dir)?;
    }

    Ok(())
}

fn run_fastqc(output_path: &Path, fastqc_args: Option<&str>, output_dir: Option<&Path>) -> Result<()> {
    let mut cmd = std::process::Command::new("fastqc");
    if let Some(args) = fastqc_args {
        for arg in args.split_whitespace() {
            cmd.arg(arg);
        }
    }
    if let Some(dir) = output_dir {
        cmd.arg("-o").arg(dir);
    }
    cmd.arg(output_path);
    eprintln!("\nRunning FastQC on {}", output_path.display());
    let status = cmd.status();
    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(s) => {
            eprintln!("Warning: FastQC exited with status {}", s);
            Ok(())
        }
        Err(e) => {
            anyhow::bail!("Failed to run FastQC: {}. Is it installed and in PATH?", e);
        }
    }
}

fn pct(part: usize, total: usize) -> f64 {
    if total == 0 { 0.0 } else { part as f64 / total as f64 * 100.0 }
}
