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
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use crate::clump::{self, ClumpLayout, MinimizerKey, canonical_minimizer, estimated_record_bytes};
use crate::fastq::{FastqReader, FastqRecord};
use crate::filters::{self, FilterResult, PairFilterResult, UnpairedLengths};
use crate::report::{PairValidationStats, TrimStats};
use crate::trimmer::{self, TrimConfig, update_adapter_stats};

/// Pipeline timing + memory counters, shared across reader/workers/main via
/// `Arc<...>`. Updated with `Relaxed` ordering — these are diagnostic
/// summaries, not control-flow — and printed at the end of the run when
/// `TG_PROFILE=1` is set in the environment, so default output is unchanged.
///
/// Per-record `Instant::now()` would add ~50 ns/record × 100M records ≈ 5 s
/// of pure overhead, so all timings are measured at coarse boundaries
/// (per-batch in the worker pool, per-`flush_bin` in the dispatcher) — except
/// the per-record minimizer scan in the reader, which is gated on
/// `profile_enabled()` to avoid the always-on cost when not profiling.
#[derive(Default, Debug)]
struct PipelineCounters {
    worker_active_ns: AtomicU64,
    worker_idle_ns: AtomicU64,
    reader_total_ns: AtomicU64,
    reader_minimizer_ns: AtomicU64,
    reader_dispatch_ns: AtomicU64,
    main_idle_ns: AtomicU64,
    batches_processed: AtomicU64,
    bin_flushes: AtomicU64,
    peak_pending_bytes: AtomicUsize,
}

/// Cached `TG_PROFILE=1` check. Read once per pipeline thread (not per record)
/// so the env-var lookup is cheap and the inner loop sees a plain `bool`.
fn profile_enabled() -> bool {
    std::env::var("TG_PROFILE").as_deref() == Ok("1")
}

impl PipelineCounters {
    fn print_summary(&self, label: &str, wall: Duration, n_workers: usize) {
        if !profile_enabled() {
            return;
        }
        let wall_ns = wall.as_nanos() as u64;
        let active = self.worker_active_ns.load(Ordering::Relaxed);
        let idle = self.worker_idle_ns.load(Ordering::Relaxed);
        let reader_total = self.reader_total_ns.load(Ordering::Relaxed);
        let reader_min = self.reader_minimizer_ns.load(Ordering::Relaxed);
        let reader_disp = self.reader_dispatch_ns.load(Ordering::Relaxed);
        let main_idle = self.main_idle_ns.load(Ordering::Relaxed);
        let theoretical_max = wall_ns.saturating_mul(n_workers as u64);
        eprintln!();
        eprintln!("=== {label} profile (TG_PROFILE=1) ===");
        eprintln!("Wall: {:.2}s   workers: {}", wall.as_secs_f64(), n_workers);
        eprintln!(
            "Workers active: {:.2}s ({:.0}% of {} × wall)   idle: {:.2}s",
            active as f64 / 1e9,
            100.0 * active as f64 / theoretical_max as f64,
            n_workers,
            idle as f64 / 1e9,
        );
        eprintln!(
            "Reader total:   {:.2}s   minimizer: {:.2}s   dispatch (channel send): {:.2}s",
            reader_total as f64 / 1e9,
            reader_min as f64 / 1e9,
            reader_disp as f64 / 1e9,
        );
        eprintln!(
            "Main idle (waiting for workers): {:.2}s ({:.0}% of wall)",
            main_idle as f64 / 1e9,
            100.0 * main_idle as f64 / wall_ns as f64,
        );
        eprintln!(
            "Batches: {}   bin flushes: {}   peak pending compressed bytes: {:.1} MB",
            self.batches_processed.load(Ordering::Relaxed),
            self.bin_flushes.load(Ordering::Relaxed),
            self.peak_pending_bytes.load(Ordering::Relaxed) as f64 / 1e6,
        );
    }
}

/// Number of read pairs per batch sent to each worker.
/// 4096 records × ~300 bytes ≈ 1.2 MB per batch — large enough to amortize
/// channel overhead and produce efficient gzip blocks, small enough for
/// even work distribution.
const BATCH_SIZE: usize = 4096;

/// Channel payload from the reader thread to a paired-end worker.
/// `None` is the poison-pill that tells a worker to exit cleanly.
type PairedWork = Option<(u64, Vec<FastqRecord>, Vec<FastqRecord>)>;
/// Channel payload from the reader thread to a single-end worker.
type SingleWork = Option<(u64, Vec<FastqRecord>)>;

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
    unpaired: UnpairedLengths,
    clump_layout: Option<ClumpLayout>,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    let retain_unpaired = unpaired_r1_path.is_some();
    let counters = Arc::new(PipelineCounters::default());
    let wall_start = Instant::now();

    // Per-worker channels (round-robin distribution — no MPMC dependency needed)
    let mut work_txs: Vec<mpsc::SyncSender<PairedWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<PairedWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<PairedWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }

    // Result channel: workers → main thread
    let (result_tx, result_rx) = mpsc::sync_channel::<PairedBatchResult>(cores * 2);

    let stats_paired =
        std::thread::scope(|s| -> Result<(TrimStats, TrimStats, PairValidationStats)> {
            // ── Worker threads ──────────────────────────────────────────────
            // Each worker owns its receiver (mpsc::Receiver is !Sync, so we
            // move them rather than borrow).
            for rx in work_rxs.drain(..) {
                let rtx = result_tx.clone();
                let counters = counters.clone();
                s.spawn(move || {
                    let mut active = Duration::ZERO;
                    let mut idle = Duration::ZERO;
                    let mut batches = 0u64;
                    let mut wait_start = Instant::now();
                    while let Ok(Some((seq, mut r1s, mut r2s))) = rx.recv() {
                        idle += wait_start.elapsed();
                        let work_start = Instant::now();
                        let result = process_paired_batch(
                            seq,
                            &mut r1s,
                            &mut r2s,
                            config,
                            gzip,
                            retain_unpaired,
                            unpaired,
                        );
                        active += work_start.elapsed();
                        batches += 1;
                        match result {
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
                        wait_start = Instant::now();
                    }
                    idle += wait_start.elapsed();
                    counters
                        .worker_active_ns
                        .fetch_add(active.as_nanos() as u64, Ordering::Relaxed);
                    counters
                        .worker_idle_ns
                        .fetch_add(idle.as_nanos() as u64, Ordering::Relaxed);
                    counters
                        .batches_processed
                        .fetch_add(batches, Ordering::Relaxed);
                });
            }
            // Drop main thread's sender so result_rx closes when all workers finish
            drop(result_tx);

            // ── Reader thread ───────────────────────────────────────────────
            // Uses open_threaded for both files: decompression runs on 2
            // background threads, this thread just batches the records.
            // Two routing modes:
            //   • Default: round-robin batches of BATCH_SIZE records to
            //     workers (input order).
            //   • Clumpy: bin records by canonical minimizer of R1; flush a
            //     bin when it fills to bin_byte_budget, sort it by minimizer
            //     key, then round-robin to a worker. Output ordering is
            //     bin-flush order, which makes each gzip member a sorted
            //     run of similar reads.
            let txs = std::mem::take(&mut work_txs);
            let reader_counters = counters.clone();
            let reader_handle = s.spawn(move || -> Result<()> {
                let reader_start = Instant::now();
                let mut reader_r1 = FastqReader::open_threaded(input_r1)?;
                let mut reader_r2 = FastqReader::open_threaded(input_r2)?;
                let result = if let Some(layout) = clump_layout {
                    read_pairs_clumpy(
                        &mut reader_r1,
                        &mut reader_r2,
                        &txs,
                        layout,
                        &reader_counters,
                    )
                } else {
                    read_pairs_round_robin(&mut reader_r1, &mut reader_r2, &txs, &reader_counters)
                };
                reader_counters
                    .reader_total_ns
                    .fetch_add(reader_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                result
            });

            // ── Main thread: ordered collection + file writing ──────────────
            let mut out_r1 = File::create(output_r1)?;
            let mut out_r2 = File::create(output_r2)?;
            let mut out_up_r1 = unpaired_r1_path.map(File::create).transpose()?;
            let mut out_up_r2 = unpaired_r2_path.map(File::create).transpose()?;

            let mut expected: u64 = 0;
            let mut pending: BTreeMap<u64, PairedBatchResult> = BTreeMap::new();
            let mut pending_bytes: usize = 0;
            let mut total_r1 = TrimStats::with_adapter_count(config.adapters.len());
            let mut total_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
            let mut total_pair = PairValidationStats::default();
            let mut main_idle = Duration::ZERO;

            loop {
                let recv_start = Instant::now();
                let Ok(result) = result_rx.recv() else { break };
                main_idle += recv_start.elapsed();
                pending_bytes += result.compressed_r1.len()
                    + result.compressed_r2.len()
                    + result.compressed_unpaired_r1.len()
                    + result.compressed_unpaired_r2.len();
                counters
                    .peak_pending_bytes
                    .fetch_max(pending_bytes, Ordering::Relaxed);
                pending.insert(result.seq, result);

                // Flush as many in-order blocks as possible
                while let Some(r) = pending.remove(&expected) {
                    pending_bytes -= r.compressed_r1.len()
                        + r.compressed_r2.len()
                        + r.compressed_unpaired_r1.len()
                        + r.compressed_unpaired_r2.len();
                    out_r1.write_all(&r.compressed_r1)?;
                    out_r2.write_all(&r.compressed_r2)?;
                    if let Some(ref mut f) = out_up_r1
                        && !r.compressed_unpaired_r1.is_empty()
                    {
                        f.write_all(&r.compressed_unpaired_r1)?;
                    }
                    if let Some(ref mut f) = out_up_r2
                        && !r.compressed_unpaired_r2.is_empty()
                    {
                        f.write_all(&r.compressed_unpaired_r2)?;
                    }
                    total_r1.merge(&r.stats_r1);
                    total_r2.merge(&r.stats_r2);
                    total_pair.merge(&r.pair_stats);
                    expected += 1;
                }
            }
            counters
                .main_idle_ns
                .fetch_add(main_idle.as_nanos() as u64, Ordering::Relaxed);

            // Propagate reader errors
            match reader_handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => return Err(e),
                Err(_) => bail!("Reader thread panicked"),
            }

            Ok((total_r1, total_r2, total_pair))
        })?;

    counters.print_summary("paired-end", wall_start.elapsed(), cores);
    Ok(stats_paired)
}
fn process_paired_batch(
    batch_seq: u64,
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    config: &TrimConfig,
    gzip: bool,
    retain_unpaired: bool,
    unpaired: UnpairedLengths,
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
        let level = config.gzip_level;
        // Scoped block: GzEncoders borrow the buffers; finish() before block
        // ends releases the borrows so we can return the buffers.
        {
            let mut gz_r1 = GzEncoder::new(&mut buf_r1, Compression::new(level));
            let mut gz_r2 = GzEncoder::new(&mut buf_r2, Compression::new(level));
            let mut gz_up_r1 = if retain_unpaired {
                Some(GzEncoder::new(&mut buf_up_r1, Compression::new(level)))
            } else {
                None
            };
            let mut gz_up_r2 = if retain_unpaired {
                Some(GzEncoder::new(&mut buf_up_r2, Compression::new(level)))
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
                unpaired,
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
            unpaired,
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

/// Outcome of trim+filter on a single read pair. Stats and the records
/// themselves are mutated in-place by `classify_paired`; this enum just
/// tells the caller what to *do* with the (now trimmed) records, so the
/// inner loop can stay sink-agnostic.
enum PairOutcome {
    /// Pair passed all filters — write both records.
    Pass,
    /// Pair dropped (discard_untrimmed, too-many-N, too-long, or
    /// too-short with no rescue side eligible).
    Discarded,
    /// Too-short pair with rescue eligibility per side. Caller writes
    /// only the eligible record(s) to its respective unpaired sink.
    Unpaired { r1_ok: bool, r2_ok: bool },
}

/// Trim + filter a single pair, mutating records and stats in place.
/// Returns what the caller should do with the (now trimmed) records.
/// Sink-agnostic: callers pick gzip writers, plain temp files, etc.
#[allow(clippy::too_many_arguments)]
fn classify_paired(
    r1: &mut FastqRecord,
    r2: &mut FastqRecord,
    config: &TrimConfig,
    stats_r1: &mut TrimStats,
    stats_r2: &mut TrimStats,
    pair_stats: &mut PairValidationStats,
    unpaired: UnpairedLengths,
) -> PairOutcome {
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
        return PairOutcome::Discarded;
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
        unpaired,
    ) {
        PairFilterResult::Pass => {
            stats_r1.total_bp_written += r1.seq.len();
            stats_r2.total_bp_written += r2.seq.len();
            stats_r1.reads_written += 1;
            stats_r2.reads_written += 1;
            PairOutcome::Pass
        }
        PairFilterResult::TooManyN => {
            pair_stats.pairs_removed += 1;
            pair_stats.pairs_removed_n += 1;
            stats_r1.too_many_n += 1;
            stats_r2.too_many_n += 1;
            PairOutcome::Discarded
        }
        PairFilterResult::TooShort { r1_ok, r2_ok } => {
            pair_stats.pairs_removed += 1;
            stats_r1.too_short += 1;
            stats_r2.too_short += 1;
            if r1_ok {
                pair_stats.r1_unpaired += 1;
            }
            if r2_ok {
                pair_stats.r2_unpaired += 1;
            }
            PairOutcome::Unpaired { r1_ok, r2_ok }
        }
        PairFilterResult::TooLong => {
            pair_stats.pairs_removed += 1;
            pair_stats.pairs_removed_too_long += 1;
            stats_r1.too_long += 1;
            stats_r2.too_long += 1;
            PairOutcome::Discarded
        }
    }
}

/// Inner loop: trim + filter + write each read pair in a batch via
/// `Write` sinks (gzip encoder or `Vec<u8>` for plain output).
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
    unpaired: UnpairedLengths,
) -> Result<()> {
    for (r1, r2) in reads_r1.iter_mut().zip(reads_r2.iter_mut()) {
        match classify_paired(r1, r2, config, stats_r1, stats_r2, pair_stats, unpaired) {
            PairOutcome::Pass => {
                r1.write_to(writer_r1)?;
                r2.write_to(writer_r2)?;
            }
            PairOutcome::Discarded => {}
            PairOutcome::Unpaired { r1_ok, r2_ok } => {
                if let Some(ref mut w) = writer_up_r1
                    && r1_ok
                {
                    r1.write_to(*w)?;
                }
                if let Some(ref mut w) = writer_up_r2
                    && r2_ok
                {
                    r2.write_to(*w)?;
                }
            }
        }
    }
    Ok(())
}

/// Round-robin paired reader: batches of `BATCH_SIZE` records dispatched in
/// input order. The default (non-clumpy) routing.
fn read_pairs_round_robin(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    txs: &[mpsc::SyncSender<PairedWork>],
    counters: &PipelineCounters,
) -> Result<()> {
    let mut seq: u64 = 0;
    let mut batch_r1 = Vec::with_capacity(BATCH_SIZE);
    let mut batch_r2 = Vec::with_capacity(BATCH_SIZE);
    let mut dispatch_dur = Duration::ZERO;

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;
        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                batch_r1.push(r1);
                batch_r2.push(r2);
                if batch_r1.len() >= BATCH_SIZE {
                    let br1 = std::mem::replace(&mut batch_r1, Vec::with_capacity(BATCH_SIZE));
                    let br2 = std::mem::replace(&mut batch_r2, Vec::with_capacity(BATCH_SIZE));
                    let idx = (seq as usize) % txs.len();
                    let send_start = Instant::now();
                    if txs[idx].send(Some((seq, br1, br2))).is_err() {
                        break;
                    }
                    dispatch_dur += send_start.elapsed();
                    seq += 1;
                }
            }
            (None, None) => {
                if !batch_r1.is_empty() {
                    let idx = (seq as usize) % txs.len();
                    let _ = txs[idx].send(Some((seq, batch_r1, batch_r2)));
                }
                for tx in txs {
                    let _ = tx.send(None);
                }
                break;
            }
            (Some(_), None) => bail!(
                "Read 2 file is truncated — R1 has more reads than R2. \
                 Please check your paired-end input files!"
            ),
            (None, Some(_)) => bail!(
                "Read 1 file is truncated — R2 has more reads than R1. \
                 Please check your paired-end input files!"
            ),
        }
    }
    counters
        .reader_dispatch_ns
        .fetch_add(dispatch_dur.as_nanos() as u64, Ordering::Relaxed);
    Ok(())
}

/// Clumpy paired reader: route each pair into a bin keyed by R1's canonical
/// minimizer; when a bin fills its byte budget, sort it by minimizer and
/// dispatch to a worker. Each flushed bin becomes one gzip member in the
/// output (re-uses the existing BTreeMap reassembly via a flush-order seq).
fn read_pairs_clumpy(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    txs: &[mpsc::SyncSender<PairedWork>],
    layout: ClumpLayout,
    counters: &PipelineCounters,
) -> Result<()> {
    let mut bins: Vec<PairedBin> = (0..layout.n_bins)
        .map(|_| PairedBin::with_budget(layout.bin_byte_budget))
        .collect();
    let mut seq: u64 = 0;
    let mut next_worker: usize = 0;
    let profile = profile_enabled();
    let mut minimizer_dur = Duration::ZERO;
    let mut dispatch_dur = Duration::ZERO;

    let flush_bin = |bin: &mut PairedBin,
                     seq: &mut u64,
                     next_worker: &mut usize,
                     dispatch_dur: &mut Duration|
     -> Result<bool> {
        let (mut r1s, mut r2s, mut keys) = bin.take();
        clump::sort_paired_by_key(&mut r1s, &mut r2s, &mut keys);
        let idx = *next_worker;
        let send_start = Instant::now();
        let send_result = txs[idx].send(Some((*seq, r1s, r2s)));
        *dispatch_dur += send_start.elapsed();
        if send_result.is_err() {
            return Ok(false);
        }
        counters.bin_flushes.fetch_add(1, Ordering::Relaxed);
        *seq += 1;
        *next_worker = (*next_worker + 1) % txs.len();
        Ok(true)
    };

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;
        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                let mz_start = profile.then(Instant::now);
                let key = canonical_minimizer(r1.seq.as_bytes());
                let bin_idx = clump::bin_for(key, layout.n_bins);
                if let Some(s) = mz_start {
                    minimizer_dur += s.elapsed();
                }
                bins[bin_idx].push(r1, r2, key);
                if bins[bin_idx].raw_bytes >= layout.bin_byte_budget
                    && !flush_bin(
                        &mut bins[bin_idx],
                        &mut seq,
                        &mut next_worker,
                        &mut dispatch_dur,
                    )?
                {
                    break;
                }
            }
            (None, None) => {
                for bin in bins.iter_mut() {
                    if !bin.is_empty()
                        && !flush_bin(bin, &mut seq, &mut next_worker, &mut dispatch_dur)?
                    {
                        break;
                    }
                }
                for tx in txs {
                    let _ = tx.send(None);
                }
                break;
            }
            (Some(_), None) => bail!(
                "Read 2 file is truncated — R1 has more reads than R2. \
                 Please check your paired-end input files!"
            ),
            (None, Some(_)) => bail!(
                "Read 1 file is truncated — R2 has more reads than R1. \
                 Please check your paired-end input files!"
            ),
        }
    }
    counters
        .reader_minimizer_ns
        .fetch_add(minimizer_dur.as_nanos() as u64, Ordering::Relaxed);
    counters
        .reader_dispatch_ns
        .fetch_add(dispatch_dur.as_nanos() as u64, Ordering::Relaxed);
    Ok(())
}

/// Per-bin buffer for the paired-end clumpy dispatcher.
#[derive(Default)]
struct PairedBin {
    r1: Vec<FastqRecord>,
    r2: Vec<FastqRecord>,
    keys: Vec<MinimizerKey>,
    raw_bytes: usize,
    budget: usize,
}

impl PairedBin {
    fn with_budget(budget: usize) -> Self {
        Self {
            budget,
            ..Self::default()
        }
    }

    /// Push one R1/R2 pair into the bin.
    ///
    /// On the first push since the bin was last flushed (`take`d), we extrapolate
    /// expected record count from this record's text size and the bin's budget,
    /// then `reserve_exact` to avoid the doubling cascade. With 32 bins all growing
    /// at once on a memory-tight host, the transient old+new spine overlap during
    /// `Vec` doublings was inflating peak memory by ~40% (measured: 8 GiB
    /// configured budget → 11.5 GiB peak footprint on 16 GiB Mac).
    fn push(&mut self, r1: FastqRecord, r2: FastqRecord, key: MinimizerKey) {
        let bytes_this = estimated_record_bytes(&r1) + estimated_record_bytes(&r2);
        if self.keys.capacity() == 0 {
            let predicted = self.budget.div_ceil(bytes_this).max(1);
            self.r1.reserve_exact(predicted);
            self.r2.reserve_exact(predicted);
            self.keys.reserve_exact(predicted);
        }
        self.raw_bytes += bytes_this;
        self.r1.push(r1);
        self.r2.push(r2);
        self.keys.push(key);
    }

    fn is_empty(&self) -> bool {
        self.r1.is_empty()
    }

    fn take(&mut self) -> (Vec<FastqRecord>, Vec<FastqRecord>, Vec<MinimizerKey>) {
        let r1 = std::mem::take(&mut self.r1);
        let r2 = std::mem::take(&mut self.r2);
        let keys = std::mem::take(&mut self.keys);
        self.raw_bytes = 0;
        (r1, r2, keys)
    }
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
    clump_layout: Option<ClumpLayout>,
) -> Result<TrimStats> {
    let counters = Arc::new(PipelineCounters::default());
    let wall_start = Instant::now();
    let mut work_txs: Vec<mpsc::SyncSender<SingleWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<SingleWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<SingleWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }

    let (result_tx, result_rx) = mpsc::sync_channel::<SingleBatchResult>(cores * 2);

    let total = std::thread::scope(|s| -> Result<TrimStats> {
        // ── Worker threads ──────────────────────────────────────────────
        for rx in work_rxs.drain(..) {
            let rtx = result_tx.clone();
            let counters = counters.clone();
            s.spawn(move || {
                let mut active = Duration::ZERO;
                let mut idle = Duration::ZERO;
                let mut batches = 0u64;
                let mut wait_start = Instant::now();
                while let Ok(Some((seq, mut reads))) = rx.recv() {
                    idle += wait_start.elapsed();
                    let work_start = Instant::now();
                    let result = process_single_batch(seq, &mut reads, config, gzip);
                    active += work_start.elapsed();
                    batches += 1;
                    match result {
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
                    wait_start = Instant::now();
                }
                idle += wait_start.elapsed();
                counters
                    .worker_active_ns
                    .fetch_add(active.as_nanos() as u64, Ordering::Relaxed);
                counters
                    .worker_idle_ns
                    .fetch_add(idle.as_nanos() as u64, Ordering::Relaxed);
                counters
                    .batches_processed
                    .fetch_add(batches, Ordering::Relaxed);
            });
        }
        drop(result_tx);

        // ── Reader thread ───────────────────────────────────────────────
        // Clumpy mode replaces round-robin with a bin dispatcher; see
        // `read_pairs_clumpy` for the paired-end shape of the same idea.
        let txs = std::mem::take(&mut work_txs);
        let reader_counters = counters.clone();
        let reader_handle = s.spawn(move || -> Result<()> {
            let reader_start = Instant::now();
            let mut reader = FastqReader::open_threaded(input)?;
            let result = if let Some(layout) = clump_layout {
                read_single_clumpy(&mut reader, &txs, layout, &reader_counters)
            } else {
                read_single_round_robin(&mut reader, &txs, &reader_counters)
            };
            reader_counters
                .reader_total_ns
                .fetch_add(reader_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
            result
        });

        // ── Main thread: ordered collection + file writing ──────────────
        let mut out = File::create(output)?;
        let mut expected: u64 = 0;
        let mut pending: BTreeMap<u64, SingleBatchResult> = BTreeMap::new();
        let mut pending_bytes: usize = 0;
        let mut total = TrimStats::with_adapter_count(config.adapters.len());
        let mut main_idle = Duration::ZERO;

        loop {
            let recv_start = Instant::now();
            let Ok(result) = result_rx.recv() else { break };
            main_idle += recv_start.elapsed();
            pending_bytes += result.compressed.len();
            counters
                .peak_pending_bytes
                .fetch_max(pending_bytes, Ordering::Relaxed);
            pending.insert(result.seq, result);
            while let Some(r) = pending.remove(&expected) {
                pending_bytes -= r.compressed.len();
                out.write_all(&r.compressed)?;
                total.merge(&r.stats);
                expected += 1;
            }
        }
        counters
            .main_idle_ns
            .fetch_add(main_idle.as_nanos() as u64, Ordering::Relaxed);

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => bail!("Reader thread panicked"),
        }

        Ok(total)
    })?;

    counters.print_summary("single-end", wall_start.elapsed(), cores);
    Ok(total)
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
            let mut gz = GzEncoder::new(&mut buf, Compression::new(config.gzip_level));
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

/// Outcome of trim+filter on a single single-end read. Like
/// `PairOutcome` but unary — `Pass` means write the (trimmed) record,
/// `Discarded` means drop it. Stats and the record are mutated in place
/// by `classify_single`.
enum SingleOutcome {
    Pass,
    Discarded,
}

/// Trim + filter a single record, mutating it and stats in place.
/// Returns whether the caller should write it. Sink-agnostic.
fn classify_single(
    record: &mut FastqRecord,
    config: &TrimConfig,
    stats: &mut TrimStats,
) -> SingleOutcome {
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
        return SingleOutcome::Discarded;
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
            stats.reads_written += 1;
            SingleOutcome::Pass
        }
        FilterResult::TooShort => {
            stats.too_short += 1;
            SingleOutcome::Discarded
        }
        FilterResult::TooLong => {
            stats.too_long += 1;
            SingleOutcome::Discarded
        }
        FilterResult::TooManyN => {
            stats.too_many_n += 1;
            SingleOutcome::Discarded
        }
    }
}

/// Inner loop: trim + filter + write each read in a batch via a
/// `Write` sink (gzip encoder or `Vec<u8>` for plain output).
fn process_reads<W: Write>(
    reads: &mut [FastqRecord],
    config: &TrimConfig,
    stats: &mut TrimStats,
    writer: &mut W,
) -> Result<()> {
    for record in reads.iter_mut() {
        if let SingleOutcome::Pass = classify_single(record, config, stats) {
            record.write_to(writer)?;
        }
    }
    Ok(())
}

/// Round-robin single-end reader: batches of `BATCH_SIZE` records dispatched
/// in input order. The default (non-clumpy) routing.
fn read_single_round_robin(
    reader: &mut FastqReader,
    txs: &[mpsc::SyncSender<SingleWork>],
    counters: &PipelineCounters,
) -> Result<()> {
    let mut seq: u64 = 0;
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut dispatch_dur = Duration::ZERO;

    while let Some(record) = reader.next_record()? {
        batch.push(record);
        if batch.len() >= BATCH_SIZE {
            let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
            let idx = (seq as usize) % txs.len();
            let send_start = Instant::now();
            if txs[idx].send(Some((seq, full))).is_err() {
                break;
            }
            dispatch_dur += send_start.elapsed();
            seq += 1;
        }
    }

    if !batch.is_empty() {
        let idx = (seq as usize) % txs.len();
        let _ = txs[idx].send(Some((seq, batch)));
    }
    for tx in txs {
        let _ = tx.send(None);
    }
    counters
        .reader_dispatch_ns
        .fetch_add(dispatch_dur.as_nanos() as u64, Ordering::Relaxed);
    Ok(())
}

/// Clumpy single-end reader: route each record into a bin keyed by its
/// canonical minimizer; when a bin fills its byte budget, sort it by
/// minimizer and dispatch to a worker.
fn read_single_clumpy(
    reader: &mut FastqReader,
    txs: &[mpsc::SyncSender<SingleWork>],
    layout: ClumpLayout,
    counters: &PipelineCounters,
) -> Result<()> {
    let mut bins: Vec<SingleBin> = (0..layout.n_bins)
        .map(|_| SingleBin::with_budget(layout.bin_byte_budget))
        .collect();
    let mut seq: u64 = 0;
    let mut next_worker: usize = 0;
    let profile = profile_enabled();
    let mut minimizer_dur = Duration::ZERO;
    let mut dispatch_dur = Duration::ZERO;

    let flush_bin = |bin: &mut SingleBin,
                     seq: &mut u64,
                     next_worker: &mut usize,
                     dispatch_dur: &mut Duration|
     -> Result<bool> {
        let (mut records, mut keys) = bin.take();
        clump::sort_single_by_key(&mut records, &mut keys);
        let idx = *next_worker;
        let send_start = Instant::now();
        let send_result = txs[idx].send(Some((*seq, records)));
        *dispatch_dur += send_start.elapsed();
        if send_result.is_err() {
            return Ok(false);
        }
        counters.bin_flushes.fetch_add(1, Ordering::Relaxed);
        *seq += 1;
        *next_worker = (*next_worker + 1) % txs.len();
        Ok(true)
    };

    while let Some(record) = reader.next_record()? {
        let mz_start = profile.then(Instant::now);
        let key = canonical_minimizer(record.seq.as_bytes());
        let bin_idx = clump::bin_for(key, layout.n_bins);
        if let Some(s) = mz_start {
            minimizer_dur += s.elapsed();
        }
        bins[bin_idx].push(record, key);
        if bins[bin_idx].raw_bytes >= layout.bin_byte_budget
            && !flush_bin(
                &mut bins[bin_idx],
                &mut seq,
                &mut next_worker,
                &mut dispatch_dur,
            )?
        {
            break;
        }
    }

    for bin in bins.iter_mut() {
        if !bin.is_empty() && !flush_bin(bin, &mut seq, &mut next_worker, &mut dispatch_dur)? {
            break;
        }
    }
    for tx in txs {
        let _ = tx.send(None);
    }
    counters
        .reader_minimizer_ns
        .fetch_add(minimizer_dur.as_nanos() as u64, Ordering::Relaxed);
    counters
        .reader_dispatch_ns
        .fetch_add(dispatch_dur.as_nanos() as u64, Ordering::Relaxed);
    Ok(())
}

/// Per-bin buffer for the single-end clumpy dispatcher.
#[derive(Default)]
struct SingleBin {
    records: Vec<FastqRecord>,
    keys: Vec<MinimizerKey>,
    raw_bytes: usize,
    budget: usize,
}

impl SingleBin {
    fn with_budget(budget: usize) -> Self {
        Self {
            budget,
            ..Self::default()
        }
    }

    /// Push one record into the bin. Pre-allocates capacity on the first push
    /// since (re)creation; see `PairedBin::push` for the rationale.
    fn push(&mut self, record: FastqRecord, key: MinimizerKey) {
        let bytes_this = estimated_record_bytes(&record);
        if self.keys.capacity() == 0 {
            let predicted = self.budget.div_ceil(bytes_this).max(1);
            self.records.reserve_exact(predicted);
            self.keys.reserve_exact(predicted);
        }
        self.raw_bytes += bytes_this;
        self.records.push(record);
        self.keys.push(key);
    }

    fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    fn take(&mut self) -> (Vec<FastqRecord>, Vec<MinimizerKey>) {
        let records = std::mem::take(&mut self.records);
        let keys = std::mem::take(&mut self.keys);
        self.raw_bytes = 0;
        (records, keys)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::{FastqReader, FastqWriter};
    use std::fs;
    use std::path::PathBuf;

    fn fresh_tmpdir(slug: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Minimal TrimConfig for the parity test: single adapter, no
    /// quality trimming or length filtering so the stats verifying
    /// equality see only the adapter-trim accounting we care about.
    fn parity_test_config() -> TrimConfig {
        TrimConfig {
            adapters: vec![("Illumina".to_string(), b"AGATCGGAAGAGC".to_vec())],
            adapters_r2: Vec::new(),
            times: 1,
            quality_cutoff: 0,
            phred_offset: 33,
            error_rate: 0.1,
            min_overlap: 1,
            length_cutoff: 0,
            max_length: None,
            max_n: None,
            trim_n: false,
            clip_r1: None,
            clip_r2: None,
            three_prime_clip_r1: None,
            three_prime_clip_r2: None,
            rename: false,
            nextseq: false,
            rrbs: false,
            non_directional: false,
            is_paired: false,
            poly_a: false,
            poly_g: false,
            discard_untrimmed: false,
            gzip_level: crate::fastq::DEFAULT_GZIP_LEVEL,
        }
    }

    /// #246 §5.2 regression: `trimmer::run_single_end` (serial path)
    /// and `parallel::run_single_end_parallel` (worker pool) must
    /// produce field-identical `TrimStats` on the same input. Beta.0/1
    /// had per-field stat drift between paths (commits 82d1e34,
    /// 3996fc5 fixed `total_bp_after_trim` / `rrbs_r2_clipped_5prime`);
    /// this test locks the invariant down. Closes the
    /// `parallel.rs` zero-tests module gap as a side effect.
    #[test]
    fn test_parallel_serial_trim_stats_parity() -> Result<()> {
        let dir = fresh_tmpdir("tg_parallel_serial_parity");
        let input_path = dir.join("input.fq");

        // 30 reads: alternating adapter-bearing and clean. Sized to
        // span more than one record but well under the 4096 batch
        // size — this exercises the merge path with a single non-full
        // batch, which is where stat drift bugs tend to live.
        let mut content = String::new();
        for i in 0..30 {
            let seq = if i % 2 == 0 {
                // 30 bp of random-looking bases + 13 bp Illumina adapter.
                "ACGTACGTACGTACGTACGTACGTACGTACAGATCGGAAGAGC"
            } else {
                // 34 bp clean (no adapter substring).
                "ACGTACGTACGTACGTACGTACGTACGTACGTAC"
            };
            let qual = "I".repeat(seq.len());
            content.push_str(&format!("@read_{i}\n{seq}\n+\n{qual}\n"));
        }
        fs::write(&input_path, content)?;

        let config = parity_test_config();

        // Serial path — open reader/writer directly, drive the
        // single-threaded `run_single_end` entry point.
        let serial_output = dir.join("serial_out.fq");
        let stats_serial = {
            let mut reader = FastqReader::open(&input_path)?;
            let mut writer =
                FastqWriter::create(&serial_output, false, 1, crate::fastq::DEFAULT_GZIP_LEVEL)?;
            let stats = crate::trimmer::run_single_end(&mut reader, &mut writer, &config)?;
            writer.flush()?;
            stats
        };

        // Parallel path — same input, 4 workers, plain output.
        let parallel_output = dir.join("parallel_out.fq");
        let stats_parallel =
            run_single_end_parallel(&input_path, &parallel_output, &config, 4, false, None)?;

        // Field-identical stats across paths. PartialEq derive on
        // TrimStats means a single equality check covers every field.
        assert_eq!(
            stats_serial, stats_parallel,
            "TrimStats must match between serial (cores=1) and parallel (cores=4) paths.\n\
             serial:   {stats_serial:#?}\n\
             parallel: {stats_parallel:#?}"
        );

        // Sanity-check that the test fixture actually exercised
        // adapter trimming — otherwise we'd be passing trivially.
        assert_eq!(stats_serial.total_reads, 30);
        assert!(
            stats_serial.total_reads_with_adapter > 0,
            "test fixture must have at least one adapter-bearing read"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// #246 follow-up: the worker pool emits each batch's output as
    /// its own independently-compressed gzip member, then concatenates
    /// them. Per RFC 1952 the result is a valid `.gz` file that
    /// decodes (via `MultiGzDecoder`) to the concatenation of every
    /// member's payload. `fastq::tests::test_multi_member_gzip_round_trip`
    /// already locks the **reader** side of this contract; this test
    /// locks the **producer** side — feeding parallel-trimmed gzip
    /// output back through `FastqReader::open` (which uses
    /// `MultiGzDecoder`) recovers every record. Sized to span more
    /// than one batch (>4096 records) so multiple gzip members are
    /// guaranteed in the output.
    #[test]
    fn test_parallel_gzip_output_decodes_multi_member() -> Result<()> {
        let dir = fresh_tmpdir("tg_parallel_gzip_multi_member");
        let input_path = dir.join("input.fq");

        let mut content = String::new();
        for i in 0..5000 {
            content.push_str(&format!(
                "@read_{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        fs::write(&input_path, content)?;

        let config = parity_test_config();
        let output_path = dir.join("out.fq.gz");
        let stats = run_single_end_parallel(&input_path, &output_path, &config, 4, true, None)?;

        assert_eq!(stats.total_reads, 5000);
        assert_eq!(stats.reads_written, 5000);

        let mut reader = FastqReader::open(&output_path)?;
        let mut count = 0;
        while reader.next_record()?.is_some() {
            count += 1;
        }
        assert_eq!(
            count, 5000,
            "all 5000 records must round-trip through the multi-member gzip output"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// #246 follow-up: single-record input is a boundary case for the
    /// worker pool — only one of the N workers does any work, the
    /// others see an empty channel and exit cleanly. Locks down that
    /// the merge path handles a single non-full batch correctly (no
    /// off-by-one on stats, no truncated output, no deadlock waiting
    /// on an extra batch from the reader).
    #[test]
    fn test_parallel_single_record_input() -> Result<()> {
        let dir = fresh_tmpdir("tg_parallel_single_record");
        let input_path = dir.join("input.fq");
        fs::write(
            &input_path,
            "@only_read\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
        )?;

        let config = parity_test_config();
        let output_path = dir.join("out.fq");
        let stats = run_single_end_parallel(&input_path, &output_path, &config, 4, false, None)?;

        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.reads_written, 1);
        assert_eq!(stats.total_reads_with_adapter, 0);

        let mut reader = FastqReader::open(&output_path)?;
        let rec = reader.next_record()?.expect("record must be present");
        assert_eq!(rec.id, "@only_read");
        assert_eq!(rec.seq, "ACGTACGTACGTACGTACGT");
        assert!(reader.next_record()?.is_none());

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// #246 follow-up: empty input must not deadlock the worker pool.
    /// The reader thread sees EOF on the first read, signals workers
    /// to exit, and the main thread returns clean zero-stats. A
    /// regression here would manifest as the test hanging — hence the
    /// pointed coverage.
    #[test]
    fn test_parallel_empty_input_no_deadlock() -> Result<()> {
        let dir = fresh_tmpdir("tg_parallel_empty");
        let input_path = dir.join("input.fq");
        fs::write(&input_path, "")?;

        let config = parity_test_config();
        let output_path = dir.join("out.fq");
        let stats = run_single_end_parallel(&input_path, &output_path, &config, 4, false, None)?;

        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.reads_written, 0);
        assert_eq!(stats.total_reads_with_adapter, 0);

        let mut reader = FastqReader::open(&output_path)?;
        assert!(
            reader.next_record()?.is_none(),
            "empty input must produce empty output (zero records)"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Higher gzip level smoke test: same input, level 1 vs level 9.
    /// The level-9 output must be smaller on disk, and decompressed
    /// output must be byte-identical between the two — gzip framing
    /// differs but payload is preserved (all gzip levels are lossless).
    #[test]
    fn test_higher_gzip_level_produces_smaller_output() -> Result<()> {
        let dir = fresh_tmpdir("tg_higher_gzip_level");
        let input_path = dir.join("input.fq");

        // 5_000 reads with a repeating-sequence motif so gzip can find
        // back-references at higher levels. Sized to span a single
        // batch so output is one gzip member regardless of cores. On
        // tiny inputs the level-9 dictionary/framing overhead can erase
        // the per-byte win, and flate2 patch bumps periodically nudge
        // that crossover; 5_000 reads gives level 9 ~325 KB to amortise
        // the overhead so the assertion is robust.
        let mut content = String::new();
        for i in 0..5_000 {
            content.push_str(&format!(
                "@read_{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        fs::write(&input_path, content)?;

        let config_default = parity_test_config();
        let mut config_high = parity_test_config();
        config_high.gzip_level = 9;

        let out_default = dir.join("out_default.fq.gz");
        let out_high = dir.join("out_high.fq.gz");

        run_single_end_parallel(&input_path, &out_default, &config_default, 2, true, None)?;
        run_single_end_parallel(&input_path, &out_high, &config_high, 2, true, None)?;

        assert!(
            fs::metadata(&out_high)?.len() < fs::metadata(&out_default)?.len(),
            "level-9 output should be smaller than level-1 default"
        );

        // Decompressed output must match: gzip levels are lossless;
        // only the framing differs.
        let read_records = |p: &std::path::Path| -> Result<Vec<(String, String, String)>> {
            let mut r = FastqReader::open(p)?;
            let mut v = Vec::new();
            while let Some(rec) = r.next_record()? {
                v.push((rec.id, rec.seq, rec.qual));
            }
            Ok(v)
        };
        assert_eq!(
            read_records(&out_default)?,
            read_records(&out_high)?,
            "decompressed records must be byte-identical between gzip levels 1 and 9"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    // ── --clumpy integration tests ───────────────────────────────────────

    use crate::clump::ClumpLayout;

    /// A modest layout for tests: enough bins to exercise multiple, but
    /// large enough per-bin that gzip can amortize its dictionary across
    /// the sorted runs (which is the point of clumpy). 1 MiB per bin
    /// matches the production minimum floor in `clump::MIN_BIN_BYTES`.
    fn small_clump_layout() -> ClumpLayout {
        ClumpLayout {
            n_bins: 16,
            bin_byte_budget: 1024 * 1024,
        }
    }

    fn read_all_records(path: &std::path::Path) -> Result<Vec<(String, String, String)>> {
        let mut r = FastqReader::open(path)?;
        let mut v = Vec::new();
        while let Some(rec) = r.next_record()? {
            v.push((rec.id, rec.seq, rec.qual));
        }
        Ok(v)
    }

    fn sorted(mut v: Vec<(String, String, String)>) -> Vec<(String, String, String)> {
        v.sort();
        v
    }

    /// Clumpy single-end: decompressed output is a permutation of the
    /// non-clumpy output (same multiset of records, possibly different
    /// order).
    #[test]
    fn test_clumpy_single_end_is_a_permutation() -> Result<()> {
        let dir = fresh_tmpdir("tg_clumpy_se_permutation");
        let input_path = dir.join("input.fq");

        // ~10K reads spanning multiple bins. Synthetic but realistic-ish:
        // mix of A-, C-, G-, T-prefix templates with random tails.
        let mut content = String::new();
        let prefixes = ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACGTACGT"];
        for i in 0..10_000_u32 {
            let prefix = prefixes[(i as usize) % prefixes.len()];
            let tail: String = (0..40)
                .map(|j| {
                    let r = i.wrapping_mul(2654435761).wrapping_add(j as u32);
                    [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                })
                .collect();
            content.push_str(&format!(
                "@read_{i}\n{prefix}{tail}\n+\n{}\n",
                "I".repeat(prefix.len() + tail.len())
            ));
        }
        std::fs::write(&input_path, &content)?;

        let config = parity_test_config();
        let plain_out = dir.join("plain.fq.gz");
        let clumpy_out = dir.join("clumpy.fq.gz");

        run_single_end_parallel(&input_path, &plain_out, &config, 4, true, None)?;
        run_single_end_parallel(
            &input_path,
            &clumpy_out,
            &config,
            4,
            true,
            Some(small_clump_layout()),
        )?;

        let plain_recs = read_all_records(&plain_out)?;
        let clumpy_recs = read_all_records(&clumpy_out)?;
        assert_eq!(plain_recs.len(), clumpy_recs.len());
        assert_eq!(
            sorted(plain_recs),
            sorted(clumpy_recs),
            "decompressed clumpy output must be a permutation of the plain output"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Clumpy paired-end: decompressed outputs are permutations of the
    /// non-clumpy outputs AND R1[i] / R2[i] are still mates after the
    /// reordering. This is the load-bearing invariant for paired clumpy.
    #[test]
    fn test_clumpy_paired_preserves_pair_lockstep() -> Result<()> {
        let dir = fresh_tmpdir("tg_clumpy_pe_lockstep");
        let r1_path = dir.join("r1.fq");
        let r2_path = dir.join("r2.fq");

        // Each pair shares an `@pair_N/1` / `@pair_N/2` ID so we can verify
        // mateship after the clumpy reorder.
        let prefixes = ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACGTACGT"];
        let mut r1_content = String::new();
        let mut r2_content = String::new();
        for i in 0..5_000_u32 {
            let prefix = prefixes[(i as usize) % prefixes.len()];
            let r1_tail: String = (0..40)
                .map(|j| {
                    let r = i.wrapping_mul(2654435761).wrapping_add(j as u32);
                    [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                })
                .collect();
            // R2 sequence different from R1 — mates aren't byte-identical.
            let r2_tail: String = (0..40)
                .map(|j| {
                    let r = i.wrapping_mul(0x9E3779B1).wrapping_add(j as u32);
                    [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                })
                .collect();
            r1_content.push_str(&format!(
                "@pair_{i}/1\n{prefix}{r1_tail}\n+\n{}\n",
                "I".repeat(prefix.len() + r1_tail.len())
            ));
            r2_content.push_str(&format!(
                "@pair_{i}/2\n{prefix}{r2_tail}\n+\n{}\n",
                "I".repeat(prefix.len() + r2_tail.len())
            ));
        }
        std::fs::write(&r1_path, r1_content)?;
        std::fs::write(&r2_path, r2_content)?;

        let mut config = parity_test_config();
        config.is_paired = true;

        let plain_r1_out = dir.join("plain_R1.fq.gz");
        let plain_r2_out = dir.join("plain_R2.fq.gz");
        let clumpy_r1_out = dir.join("clumpy_R1.fq.gz");
        let clumpy_r2_out = dir.join("clumpy_R2.fq.gz");

        run_paired_end_parallel(
            &r1_path,
            &r2_path,
            &plain_r1_out,
            &plain_r2_out,
            None,
            None,
            &config,
            4,
            true,
            UnpairedLengths { r1: 35, r2: 35 },
            None,
        )?;
        run_paired_end_parallel(
            &r1_path,
            &r2_path,
            &clumpy_r1_out,
            &clumpy_r2_out,
            None,
            None,
            &config,
            4,
            true,
            UnpairedLengths { r1: 35, r2: 35 },
            Some(small_clump_layout()),
        )?;

        let plain_r1 = read_all_records(&plain_r1_out)?;
        let plain_r2 = read_all_records(&plain_r2_out)?;
        let clumpy_r1 = read_all_records(&clumpy_r1_out)?;
        let clumpy_r2 = read_all_records(&clumpy_r2_out)?;

        // Permutation property on both files.
        assert_eq!(sorted(plain_r1.clone()), sorted(clumpy_r1.clone()));
        assert_eq!(sorted(plain_r2.clone()), sorted(clumpy_r2.clone()));

        // Pair lockstep: R1[i].id and R2[i].id must still share the same
        // `pair_N` stem. This is what would silently break if the bin
        // dispatcher routed R2 independently.
        assert_eq!(clumpy_r1.len(), clumpy_r2.len());
        for (a, b) in clumpy_r1.iter().zip(clumpy_r2.iter()) {
            let stem_a = a.0.trim_end_matches("/1");
            let stem_b = b.0.trim_end_matches("/2");
            assert_eq!(stem_a, stem_b, "pair lockstep broken: {} ≠ {}", a.0, b.0);
        }

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Clumpy must produce identical TrimStats to the non-clumpy run on
    /// the same input. Filtering and stat-merging are order-independent
    /// so this is the precise contract.
    #[test]
    fn test_clumpy_stats_parity_with_plain() -> Result<()> {
        let dir = fresh_tmpdir("tg_clumpy_stats_parity");
        let input_path = dir.join("input.fq");

        // Mix of adapter-bearing and clean reads, like the parallel parity
        // test but larger to span multiple bin flushes.
        let mut content = String::new();
        for i in 0..6_000_u32 {
            let seq = if i % 3 == 0 {
                "ACGTACGTACGTACGTACGTACGTACGTACAGATCGGAAGAGC"
            } else {
                "ACGTACGTACGTACGTACGTACGTACGTACGTAC"
            };
            content.push_str(&format!("@read_{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
        }
        std::fs::write(&input_path, content)?;

        let config = parity_test_config();
        let plain_out = dir.join("plain.fq.gz");
        let clumpy_out = dir.join("clumpy.fq.gz");

        let plain_stats = run_single_end_parallel(&input_path, &plain_out, &config, 4, true, None)?;
        let clumpy_stats = run_single_end_parallel(
            &input_path,
            &clumpy_out,
            &config,
            4,
            true,
            Some(small_clump_layout()),
        )?;

        assert_eq!(
            plain_stats, clumpy_stats,
            "stats must be field-identical between plain and clumpy paths"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// On highly clusterable input — many copies of a few templates with
    /// small mutations — clumpy output must be measurably smaller than
    /// plain. Loose bound (20%) so the test isn't fragile against future
    /// gzip / minimizer tweaks.
    #[test]
    fn test_clumpy_compresses_clusterable_input() -> Result<()> {
        let dir = fresh_tmpdir("tg_clumpy_compresses");
        let input_path = dir.join("input.fq");

        // 100 templates, 100 copies each = 10K reads. Each copy gets one
        // base mutated to keep gzip honest (trivially-identical reads
        // would over-state the win).
        let templates: Vec<String> = (0..100)
            .map(|t| {
                let seed = t as u64;
                (0..120)
                    .map(|i| {
                        let r = seed.wrapping_mul(2654435761).wrapping_add(i);
                        [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                    })
                    .collect()
            })
            .collect();

        let mut content = String::new();
        let mut counter: u32 = 0;
        for (t_idx, template) in templates.iter().enumerate() {
            for copy in 0..100 {
                let mut seq = template.clone();
                // Mutate one base per copy at a deterministic position.
                let pos = (copy * 7) % seq.len();
                let new_base = match seq.as_bytes()[pos] {
                    b'A' => 'C',
                    b'C' => 'G',
                    b'G' => 'T',
                    _ => 'A',
                };
                seq.replace_range(pos..pos + 1, &new_base.to_string());
                content.push_str(&format!(
                    "@read_{t_idx}_{copy}_{counter}\n{seq}\n+\n{}\n",
                    "I".repeat(seq.len())
                ));
                counter += 1;
            }
        }
        std::fs::write(&input_path, content)?;

        let config = parity_test_config();
        let plain_out = dir.join("plain.fq.gz");
        let clumpy_out = dir.join("clumpy.fq.gz");

        run_single_end_parallel(&input_path, &plain_out, &config, 4, true, None)?;
        run_single_end_parallel(
            &input_path,
            &clumpy_out,
            &config,
            4,
            true,
            Some(small_clump_layout()),
        )?;

        let plain_size = fs::metadata(&plain_out)?.len();
        let clumpy_size = fs::metadata(&clumpy_out)?.len();

        assert!(
            (clumpy_size as f64) < 0.80 * (plain_size as f64),
            "clumpy output ({clumpy_size}) must be < 80% of plain ({plain_size}) on \
             highly-clusterable input — got ratio {:.2}",
            clumpy_size as f64 / plain_size as f64
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Tiny input (fewer records than one bin's byte budget) must still
    /// produce correct, decompressible output. Exercises the EOF flush
    /// path where every bin is partial.
    #[test]
    fn test_clumpy_tiny_input() -> Result<()> {
        let dir = fresh_tmpdir("tg_clumpy_tiny");
        let input_path = dir.join("input.fq");

        // 5 records — well under any bin's byte budget.
        let mut content = String::new();
        for i in 0..5 {
            content.push_str(&format!(
                "@read_{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        std::fs::write(&input_path, content)?;

        let config = parity_test_config();
        let out = dir.join("out.fq.gz");
        let stats = run_single_end_parallel(
            &input_path,
            &out,
            &config,
            2,
            true,
            Some(small_clump_layout()),
        )?;

        assert_eq!(stats.total_reads, 5);
        assert_eq!(stats.reads_written, 5);

        // Round-trip must still decode (multi-member gzip is fine).
        let mut count = 0;
        let mut r = FastqReader::open(&out)?;
        while r.next_record()?.is_some() {
            count += 1;
        }
        assert_eq!(count, 5);

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }
}
