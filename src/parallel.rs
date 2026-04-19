//! Worker-pool parallelism for trimming pipelines.
//!
//! When `--cores N` is specified (N > 1), reads are distributed across N
//! worker threads. Each worker handles trimming **and** gzip compression,
//! producing independently-compressed gzip blocks. The blocks are written
//! to output files in sequence order, producing valid gzip files (RFC 1952
//! permits concatenation of gzip members).
//!
//! This replaces the previous pipeline model (readers → single main → gzp
//! writers) with a worker-pool model that scales wall-time linearly with
//! core count, because the dominant cost (gzip compression, ~60% of runtime)
//! is distributed across workers instead of funneled through one thread.

use anyhow::{Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc;

use crate::fastq::{FastqReader, FastqRecord};
use crate::filters::{self, FilterResult, PairFilterResult};
use crate::report::{PairValidationStats, TrimStats};
use crate::trimmer::{self, TrimConfig, update_adapter_stats};

/// Number of read pairs per batch sent to each worker.
/// 4096 records × ~300 bytes ≈ 1.2 MB per batch — large enough to amortize
/// channel overhead and produce efficient gzip blocks, small enough for
/// even work distribution.
const BATCH_SIZE: usize = 4096;

// ─────────────────────────────── Paired-end ───────────────────────────────

/// Result of processing one batch of read pairs.
struct PairedBatchResult {
    seq: u64,
    compressed_r1: Vec<u8>,
    compressed_r2: Vec<u8>,
    compressed_unpaired_r1: Vec<u8>,
    compressed_unpaired_r2: Vec<u8>,
    stats_r1: TrimStats,
    stats_r2: TrimStats,
    pair_stats: PairValidationStats,
}

/// Run the paired-end trimming pipeline with N parallel workers.
///
/// Architecture:
/// - 2 background decompression threads (one per input file, via `open_threaded`)
/// - 1 reader/batcher thread (pulls decompressed records, creates numbered batches)
/// - N worker threads (each: trim both reads + compress output into gzip block)
/// - Main thread (collects compressed blocks in order, writes to output files)
///
/// Total threads: N + 4. Each worker independently handles trim + compress,
/// so the dominant bottleneck (gzip compression) scales linearly with N.
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end_parallel(
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
    unpaired_r1_path: Option<&Path>,
    unpaired_r2_path: Option<&Path>,
    config: &TrimConfig,
    cores: usize,
    gzip: bool,
    unpaired_length_r1: usize,
    unpaired_length_r2: usize,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    let retain_unpaired = unpaired_r1_path.is_some();

    // Per-worker channels (round-robin distribution — no MPMC dependency needed)
    type PairedWork = Option<(u64, Vec<FastqRecord>, Vec<FastqRecord>)>;
    let mut work_txs: Vec<mpsc::SyncSender<PairedWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<PairedWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<PairedWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }

    // Result channel: workers → main thread
    let (result_tx, result_rx) = mpsc::sync_channel::<PairedBatchResult>(cores * 2);

    std::thread::scope(|s| -> Result<(TrimStats, TrimStats, PairValidationStats)> {
        // ── Worker threads ──────────────────────────────────────────────
        // Each worker owns its receiver (mpsc::Receiver is !Sync, so we
        // move them rather than borrow).
        for rx in work_rxs.drain(..) {
            let rtx = result_tx.clone();
            s.spawn(move || {
                while let Ok(Some((seq, mut r1s, mut r2s))) = rx.recv() {
                    match process_paired_batch(
                        seq,
                        &mut r1s,
                        &mut r2s,
                        config,
                        gzip,
                        retain_unpaired,
                        unpaired_length_r1,
                        unpaired_length_r2,
                    ) {
                        Ok(result) => {
                            if rtx.send(result).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            eprintln!("Worker error: {}", e);
                            break;
                        }
                    }
                }
            });
        }
        // Drop main thread's sender so result_rx closes when all workers finish
        drop(result_tx);

        // ── Reader thread ───────────────────────────────────────────────
        // Uses open_threaded for both files: decompression runs on 2
        // background threads, this thread just batches the records.
        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader_r1 = FastqReader::open_threaded(input_r1)?;
            let mut reader_r2 = FastqReader::open_threaded(input_r2)?;
            let mut seq: u64 = 0;
            let mut batch_r1 = Vec::with_capacity(BATCH_SIZE);
            let mut batch_r2 = Vec::with_capacity(BATCH_SIZE);

            loop {
                let rec1 = reader_r1.next_record()?;
                let rec2 = reader_r2.next_record()?;
                match (rec1, rec2) {
                    (Some(r1), Some(r2)) => {
                        batch_r1.push(r1);
                        batch_r2.push(r2);
                        if batch_r1.len() >= BATCH_SIZE {
                            let br1 =
                                std::mem::replace(&mut batch_r1, Vec::with_capacity(BATCH_SIZE));
                            let br2 =
                                std::mem::replace(&mut batch_r2, Vec::with_capacity(BATCH_SIZE));
                            let idx = (seq as usize) % txs.len();
                            if txs[idx].send(Some((seq, br1, br2))).is_err() {
                                break;
                            }
                            seq += 1;
                        }
                    }
                    (None, None) => {
                        // Send final partial batch
                        if !batch_r1.is_empty() {
                            let idx = (seq as usize) % txs.len();
                            let _ = txs[idx].send(Some((seq, batch_r1, batch_r2)));
                        }
                        // Poison pill per worker
                        for tx in &txs {
                            let _ = tx.send(None);
                        }
                        break;
                    }
                    (Some(_), None) => {
                        bail!(
                            "Read 2 file is truncated — R1 has more reads than R2. \
                             Please check your paired-end input files!"
                        );
                    }
                    (None, Some(_)) => {
                        bail!(
                            "Read 1 file is truncated — R2 has more reads than R1. \
                             Please check your paired-end input files!"
                        );
                    }
                }
            }
            Ok(())
        });

        // ── Main thread: ordered collection + file writing ──────────────
        let mut out_r1 = File::create(output_r1)?;
        let mut out_r2 = File::create(output_r2)?;
        let mut out_up_r1 = unpaired_r1_path.map(File::create).transpose()?;
        let mut out_up_r2 = unpaired_r2_path.map(File::create).transpose()?;

        let mut expected: u64 = 0;
        let mut pending: BTreeMap<u64, PairedBatchResult> = BTreeMap::new();
        let mut total_r1 = TrimStats::with_adapter_count(config.adapters.len());
        let mut total_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
        let mut total_pair = PairValidationStats::default();

        while let Ok(result) = result_rx.recv() {
            pending.insert(result.seq, result);

            // Flush as many in-order blocks as possible
            while let Some(r) = pending.remove(&expected) {
                out_r1.write_all(&r.compressed_r1)?;
                out_r2.write_all(&r.compressed_r2)?;
                if let Some(ref mut f) = out_up_r1 {
                    if !r.compressed_unpaired_r1.is_empty() {
                        f.write_all(&r.compressed_unpaired_r1)?;
                    }
                }
                if let Some(ref mut f) = out_up_r2 {
                    if !r.compressed_unpaired_r2.is_empty() {
                        f.write_all(&r.compressed_unpaired_r2)?;
                    }
                }
                total_r1.merge(&r.stats_r1);
                total_r2.merge(&r.stats_r2);
                total_pair.merge(&r.pair_stats);
                expected += 1;
            }
        }

        // Propagate reader errors
        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => bail!("Reader thread panicked"),
        }

        Ok((total_r1, total_r2, total_pair))
    })
}

/// Process a batch of paired reads: trim, filter, compress into gzip blocks.
#[allow(clippy::too_many_arguments)]
fn process_paired_batch(
    batch_seq: u64,
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    config: &TrimConfig,
    gzip: bool,
    retain_unpaired: bool,
    unpaired_length_r1: usize,
    unpaired_length_r2: usize,
) -> Result<PairedBatchResult> {
    let mut stats_r1 = TrimStats::with_adapter_count(config.adapters.len());
    let mut stats_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
    let mut pair_stats = PairValidationStats::default();

    let cap = reads_r1.len() * 300;
    let mut buf_r1 = Vec::with_capacity(cap);
    let mut buf_r2 = Vec::with_capacity(cap);
    let mut buf_up_r1 = Vec::new();
    let mut buf_up_r2 = Vec::new();

    if gzip {
        // Scoped block: GzEncoders borrow the buffers; finish() before block
        // ends releases the borrows so we can return the buffers.
        {
            let mut gz_r1 = GzEncoder::new(&mut buf_r1, Compression::new(6));
            let mut gz_r2 = GzEncoder::new(&mut buf_r2, Compression::new(6));
            let mut gz_up_r1 = if retain_unpaired {
                Some(GzEncoder::new(&mut buf_up_r1, Compression::new(6)))
            } else {
                None
            };
            let mut gz_up_r2 = if retain_unpaired {
                Some(GzEncoder::new(&mut buf_up_r2, Compression::new(6)))
            } else {
                None
            };

            process_pairs(
                reads_r1,
                reads_r2,
                config,
                &mut stats_r1,
                &mut stats_r2,
                &mut pair_stats,
                &mut gz_r1,
                &mut gz_r2,
                gz_up_r1.as_mut(),
                gz_up_r2.as_mut(),
                unpaired_length_r1,
                unpaired_length_r2,
            )?;

            gz_r1.finish()?;
            gz_r2.finish()?;
            if let Some(gz) = gz_up_r1 {
                gz.finish()?;
            }
            if let Some(gz) = gz_up_r2 {
                gz.finish()?;
            }
        }
    } else {
        process_pairs(
            reads_r1,
            reads_r2,
            config,
            &mut stats_r1,
            &mut stats_r2,
            &mut pair_stats,
            &mut buf_r1,
            &mut buf_r2,
            if retain_unpaired {
                Some(&mut buf_up_r1)
            } else {
                None
            },
            if retain_unpaired {
                Some(&mut buf_up_r2)
            } else {
                None
            },
            unpaired_length_r1,
            unpaired_length_r2,
        )?;
    }

    Ok(PairedBatchResult {
        seq: batch_seq,
        compressed_r1: buf_r1,
        compressed_r2: buf_r2,
        compressed_unpaired_r1: buf_up_r1,
        compressed_unpaired_r2: buf_up_r2,
        stats_r1,
        stats_r2,
        pair_stats,
    })
}

/// Inner loop: trim + filter + write each read pair in a batch.
///
/// Generic over `W: Write` so it works with both `GzEncoder` (gzip output)
/// and `Vec<u8>` (plain text output) without runtime dispatch.
#[allow(clippy::too_many_arguments)]
fn process_pairs<W: Write>(
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    config: &TrimConfig,
    stats_r1: &mut TrimStats,
    stats_r2: &mut TrimStats,
    pair_stats: &mut PairValidationStats,
    writer_r1: &mut W,
    writer_r2: &mut W,
    mut writer_up_r1: Option<&mut W>,
    mut writer_up_r2: Option<&mut W>,
    unpaired_length_r1: usize,
    unpaired_length_r2: usize,
) -> Result<()> {
    for (r1, r2) in reads_r1.iter_mut().zip(reads_r2.iter_mut()) {
        stats_r1.total_reads += 1;
        stats_r2.total_reads += 1;
        stats_r1.total_bp_processed += r1.seq.len();
        stats_r2.total_bp_processed += r2.seq.len();
        pair_stats.pairs_analyzed += 1;

        let res_r1 = trimmer::trim_read(r1, config, false);
        let res_r2 = trimmer::trim_read(r2, config, true);

        stats_r1.bases_quality_trimmed += res_r1.quality_trimmed_bp;
        stats_r2.bases_quality_trimmed += res_r2.quality_trimmed_bp;
        update_adapter_stats(stats_r1, &res_r1);
        update_adapter_stats(stats_r2, &res_r2);
        if res_r1.rrbs_trimmed_3prime {
            stats_r1.rrbs_trimmed_3prime += 1;
        }
        if res_r1.rrbs_trimmed_5prime {
            stats_r1.rrbs_trimmed_5prime += 1;
        }
        if res_r2.rrbs_trimmed_3prime {
            stats_r2.rrbs_trimmed_3prime += 1;
        }
        if res_r2.rrbs_trimmed_5prime {
            stats_r2.rrbs_trimmed_5prime += 1;
        }
        if config.rrbs && res_r2.clip_5prime_applied {
            stats_r2.rrbs_r2_clipped_5prime += 1;
        }
        if res_r1.poly_a_trimmed > 0 {
            stats_r1.poly_a_trimmed += 1;
            stats_r1.poly_a_bases_trimmed += res_r1.poly_a_trimmed;
        }
        if res_r2.poly_a_trimmed > 0 {
            stats_r2.poly_a_trimmed += 1;
            stats_r2.poly_a_bases_trimmed += res_r2.poly_a_trimmed;
        }
        if res_r1.poly_g_trimmed > 0 {
            stats_r1.poly_g_trimmed += 1;
            stats_r1.poly_g_bases_trimmed += res_r1.poly_g_trimmed;
        }
        if res_r2.poly_g_trimmed > 0 {
            stats_r2.poly_g_trimmed += 1;
            stats_r2.poly_g_bases_trimmed += res_r2.poly_g_trimmed;
        }

        if config.discard_untrimmed && !res_r1.had_adapter && !res_r2.had_adapter {
            stats_r1.discarded_untrimmed += 1;
            stats_r2.discarded_untrimmed += 1;
            pair_stats.pairs_removed += 1;
            continue;
        }

        // Track bp after trimming but before pair/length filtering (Cutadapt-compatible stat)
        stats_r1.total_bp_after_trim += r1.seq.len();
        stats_r2.total_bp_after_trim += r2.seq.len();

        match filters::filter_paired_end(
            r1,
            r2,
            config.length_cutoff,
            config.max_length,
            config.max_n.clone(),
            unpaired_length_r1,
            unpaired_length_r2,
        ) {
            PairFilterResult::Pass => {
                stats_r1.total_bp_written += r1.seq.len();
                stats_r2.total_bp_written += r2.seq.len();
                r1.write_to(writer_r1)?;
                r2.write_to(writer_r2)?;
                stats_r1.reads_written += 1;
                stats_r2.reads_written += 1;
            }
            PairFilterResult::TooManyN => {
                pair_stats.pairs_removed += 1;
                pair_stats.pairs_removed_n += 1;
                stats_r1.too_many_n += 1;
                stats_r2.too_many_n += 1;
            }
            PairFilterResult::TooShort { r1_ok, r2_ok } => {
                pair_stats.pairs_removed += 1;
                stats_r1.too_short += 1;
                stats_r2.too_short += 1;
                if let Some(ref mut w) = writer_up_r1 {
                    if r1_ok {
                        r1.write_to(*w)?;
                        pair_stats.r1_unpaired += 1;
                    }
                }
                if let Some(ref mut w) = writer_up_r2 {
                    if r2_ok {
                        r2.write_to(*w)?;
                        pair_stats.r2_unpaired += 1;
                    }
                }
            }
            PairFilterResult::TooLong => {
                pair_stats.pairs_removed += 1;
                pair_stats.pairs_removed_too_long += 1;
                stats_r1.too_long += 1;
                stats_r2.too_long += 1;
            }
        }
    }
    Ok(())
}

// ─────────────────────────────── Single-end ───────────────────────────────

/// Result of processing one batch of single-end reads.
struct SingleBatchResult {
    seq: u64,
    compressed: Vec<u8>,
    stats: TrimStats,
}

/// Run the single-end trimming pipeline with N parallel workers.
pub fn run_single_end_parallel(
    input: &Path,
    output: &Path,
    config: &TrimConfig,
    cores: usize,
    gzip: bool,
) -> Result<TrimStats> {
    type SingleWork = Option<(u64, Vec<FastqRecord>)>;
    let mut work_txs: Vec<mpsc::SyncSender<SingleWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<SingleWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<SingleWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }

    let (result_tx, result_rx) = mpsc::sync_channel::<SingleBatchResult>(cores * 2);

    std::thread::scope(|s| -> Result<TrimStats> {
        // ── Worker threads ──────────────────────────────────────────────
        for rx in work_rxs.drain(..) {
            let rtx = result_tx.clone();
            s.spawn(move || {
                while let Ok(Some((seq, mut reads))) = rx.recv() {
                    match process_single_batch(seq, &mut reads, config, gzip) {
                        Ok(result) => {
                            if rtx.send(result).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            eprintln!("Worker error: {}", e);
                            break;
                        }
                    }
                }
            });
        }
        drop(result_tx);

        // ── Reader thread ───────────────────────────────────────────────
        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = FastqReader::open_threaded(input)?;
            let mut seq: u64 = 0;
            let mut batch = Vec::with_capacity(BATCH_SIZE);

            while let Some(record) = reader.next_record()? {
                batch.push(record);
                if batch.len() >= BATCH_SIZE {
                    let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
                    let idx = (seq as usize) % txs.len();
                    if txs[idx].send(Some((seq, full))).is_err() {
                        break;
                    }
                    seq += 1;
                }
            }

            if !batch.is_empty() {
                let idx = (seq as usize) % txs.len();
                let _ = txs[idx].send(Some((seq, batch)));
            }
            for tx in &txs {
                let _ = tx.send(None);
            }
            Ok(())
        });

        // ── Main thread: ordered collection + file writing ──────────────
        let mut out = File::create(output)?;
        let mut expected: u64 = 0;
        let mut pending: BTreeMap<u64, SingleBatchResult> = BTreeMap::new();
        let mut total = TrimStats::with_adapter_count(config.adapters.len());

        while let Ok(result) = result_rx.recv() {
            pending.insert(result.seq, result);
            while let Some(r) = pending.remove(&expected) {
                out.write_all(&r.compressed)?;
                total.merge(&r.stats);
                expected += 1;
            }
        }

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => bail!("Reader thread panicked"),
        }

        Ok(total)
    })
}

/// Process a batch of single-end reads: trim, filter, compress.
fn process_single_batch(
    batch_seq: u64,
    reads: &mut [FastqRecord],
    config: &TrimConfig,
    gzip: bool,
) -> Result<SingleBatchResult> {
    let mut stats = TrimStats::with_adapter_count(config.adapters.len());
    let cap = reads.len() * 300;
    let mut buf = Vec::with_capacity(cap);

    if gzip {
        {
            let mut gz = GzEncoder::new(&mut buf, Compression::new(6));
            process_reads(reads, config, &mut stats, &mut gz)?;
            gz.finish()?;
        }
    } else {
        process_reads(reads, config, &mut stats, &mut buf)?;
    }

    Ok(SingleBatchResult {
        seq: batch_seq,
        compressed: buf,
        stats,
    })
}

/// Inner loop: trim + filter + write each read in a batch.
fn process_reads<W: Write>(
    reads: &mut [FastqRecord],
    config: &TrimConfig,
    stats: &mut TrimStats,
    writer: &mut W,
) -> Result<()> {
    for record in reads.iter_mut() {
        stats.total_reads += 1;
        stats.total_bp_processed += record.seq.len();

        let result = trimmer::trim_read(record, config, false);
        stats.bases_quality_trimmed += result.quality_trimmed_bp;
        update_adapter_stats(stats, &result);
        if result.rrbs_trimmed_3prime {
            stats.rrbs_trimmed_3prime += 1;
        }
        if result.rrbs_trimmed_5prime {
            stats.rrbs_trimmed_5prime += 1;
        }
        if result.poly_a_trimmed > 0 {
            stats.poly_a_trimmed += 1;
            stats.poly_a_bases_trimmed += result.poly_a_trimmed;
        }
        if result.poly_g_trimmed > 0 {
            stats.poly_g_trimmed += 1;
            stats.poly_g_bases_trimmed += result.poly_g_trimmed;
        }

        if config.discard_untrimmed && !result.had_adapter {
            stats.discarded_untrimmed += 1;
            continue;
        }

        // Track bp after trimming but before length filtering (Cutadapt-compatible stat)
        stats.total_bp_after_trim += record.seq.len();

        match filters::filter_single_end(
            record,
            config.length_cutoff,
            config.max_length,
            config.max_n.clone(),
        ) {
            FilterResult::Pass => {
                stats.total_bp_written += record.seq.len();
                record.write_to(writer)?;
                stats.reads_written += 1;
            }
            FilterResult::TooShort => stats.too_short += 1,
            FilterResult::TooLong => stats.too_long += 1,
            FilterResult::TooManyN => stats.too_many_n += 1,
        }
    }
    Ok(())
}
