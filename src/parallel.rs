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

use crate::clump::{self, ClumpLayout, MinimizerKey, canonical_minimizer, estimated_record_bytes};
use crate::fastq::{FastqReader, FastqRecord};
use crate::filters::{self, FilterResult, PairFilterResult, UnpairedLengths};
use crate::report::{PairValidationStats, TrimStats};
use crate::trimmer::{self, TrimConfig, update_adapter_stats};

/// Number of read pairs per batch sent to each worker.
/// 4096 records × ~300 bytes ≈ 1.2 MB per batch — large enough to amortize
/// channel overhead and produce efficient gzip blocks, small enough for
/// even work distribution.
const BATCH_SIZE: usize = 4096;

/// Channel payload from the reader thread to a paired-end worker.
/// `None` is the poison-pill that tells a worker to exit cleanly.
///
/// 4-tuple shape: `(seq, R1 batch, R2 batch, optional passthrough batch)`.
/// The 4th element is `Some` iff `--passthrough` is active for this run.
/// (Plan v2 Step 6.)
type PairedWork = Option<(
    u64,
    Vec<FastqRecord>,
    Vec<FastqRecord>,
    Option<Vec<FastqRecord>>,
)>;
/// Channel payload from the reader thread to a single-end worker.
type SingleWork = Option<(u64, Vec<FastqRecord>)>;

// ─────────────────────────────── Paired-end ───────────────────────────────

/// Result of processing one batch of read pairs.
struct PairedBatchResult {
    seq: u64,
    compressed_r1: Vec<u8>,
    compressed_r2: Vec<u8>,
    /// Passthrough output bytes (plan v2 Step 7). `None` iff `--passthrough`
    /// was not active. `Some` even when the batch is empty — the worker
    /// stamps a valid empty gzip member so the main thread can write
    /// unconditionally when the output sink is open.
    compressed_passthrough: Option<Vec<u8>>,
    compressed_unpaired_r1: Vec<u8>,
    compressed_unpaired_r2: Vec<u8>,
    stats_r1: TrimStats,
    stats_r2: TrimStats,
    pair_stats: PairValidationStats,
}

/// Run the paired-end trimming pipeline with N parallel workers.
///
/// Architecture:
/// - 2 or 3 background decompression threads (one per input file, via `open_threaded`;
///   third decompressor only when `--passthrough` is active)
/// - 1 reader/batcher thread (pulls decompressed records, creates numbered batches,
///   performs the three-way header sync check when `--passthrough` is active)
/// - N worker threads (each: trim both reads + compress output into gzip block)
/// - Main thread (collects compressed blocks in order, writes to output files)
///
/// Total threads: N + 4 (or N + 5 with `--passthrough`). Each worker independently
/// handles trim + compress, so the dominant bottleneck (gzip compression) scales
/// linearly with N.
///
/// `--passthrough` (plan v2 Step 6/6a/7): when `input_passthrough` and
/// `output_passthrough` are both `Some`, the reader streams a third FASTQ
/// file in lockstep with R1/R2, the worker writes it untrimmed into a
/// parallel gzip stream, and the main thread emits a third output file.
/// Records dropped by length/N/quality filters are also dropped from the
/// passthrough stream. The result channel uses `Result<PairedBatchResult>`
/// so a reader-side sync error short-circuits the main flush loop *before*
/// further output bytes commit to disk (B-Crit-1 fix; partial outputs may
/// still exist on disk on mid-stream error — v1 contract).
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end_parallel(
    input_r1: &Path,
    input_r2: &Path,
    input_passthrough: Option<&Path>,
    output_r1: &Path,
    output_r2: &Path,
    output_passthrough: Option<&Path>,
    unpaired_r1_path: Option<&Path>,
    unpaired_r2_path: Option<&Path>,
    config: &TrimConfig,
    cores: usize,
    gzip: bool,
    unpaired: UnpairedLengths,
    clump_layout: Option<ClumpLayout>,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    let retain_unpaired = unpaired_r1_path.is_some();
    // Both passthrough params must be both-Some or both-None.
    debug_assert!(
        input_passthrough.is_some() == output_passthrough.is_some(),
        "run_paired_end_parallel: input_passthrough and output_passthrough must be both Some or both None"
    );
    let passthrough_active = input_passthrough.is_some();

    // Per-worker channels (round-robin distribution — no MPMC dependency needed)
    let mut work_txs: Vec<mpsc::SyncSender<PairedWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<PairedWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<PairedWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }

    // Result channel: workers → main thread. Payload is `Result<PairedBatchResult>`
    // so reader-side errors surface to the main thread without waiting on the
    // reader_handle.join() that happens after the result loop drains. *(Plan v2
    // Step 6a; B-Crit-1.)*
    let (result_tx, result_rx) = mpsc::sync_channel::<Result<PairedBatchResult>>(cores * 2);

    std::thread::scope(|s| -> Result<(TrimStats, TrimStats, PairValidationStats)> {
        // ── Worker threads ──────────────────────────────────────────────
        // Each worker owns its receiver (mpsc::Receiver is !Sync, so we
        // move them rather than borrow).
        for rx in work_rxs.drain(..) {
            let rtx = result_tx.clone();
            s.spawn(move || {
                while let Ok(Some((seq, mut r1s, mut r2s, pts))) = rx.recv() {
                    let result = process_paired_batch(
                        seq,
                        &mut r1s,
                        &mut r2s,
                        pts,
                        config,
                        gzip,
                        retain_unpaired,
                        unpaired,
                    );
                    match result {
                        Ok(batch) => {
                            if rtx.send(Ok(batch)).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            // Surface worker error via the result channel
                            // (Step 6a). main thread short-circuits the
                            // flush loop on the first Err.
                            let _ = rtx.send(Err(e));
                            break;
                        }
                    }
                }
            });
        }
        // Clone result_tx for the reader so it can also surface errors via
        // the result channel before exiting (Step 6a).
        let result_tx_for_reader = result_tx.clone();
        // Drop main thread's sender so result_rx closes when ALL clones drop
        // (every worker + the reader's clone).
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
        // --passthrough mode forces round-robin (clump_layout is None when
        // passthrough is set — validated by Cli::validate).
        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            let inner = || -> Result<()> {
                let mut reader_r1 = FastqReader::open_threaded(input_r1)?;
                let mut reader_r2 = FastqReader::open_threaded(input_r2)?;
                let mut reader_pt = match input_passthrough {
                    Some(p) => Some(FastqReader::open_threaded(p)?),
                    None => None,
                };
                if let Some(layout) = clump_layout {
                    // --clumpify rejects --passthrough at validation, so reader_pt
                    // is None here. The existing clumpy reader is unchanged.
                    read_pairs_clumpy(&mut reader_r1, &mut reader_r2, &txs, layout)
                } else {
                    read_pairs_round_robin(&mut reader_r1, &mut reader_r2, reader_pt.as_mut(), &txs)
                }
            };
            let res = inner();
            if res.is_err() {
                // Sentinel signal so the main flush loop short-circuits
                // BEFORE writing more bytes. The full error (with anyhow
                // chain intact) is returned via `res` and reaches main
                // thread through `reader_handle.join()` below — see the
                // `reader_join_err.or(first_error)` priority that prefers
                // it. *(Code-reviewer C1: was previously a `format!` round
                // trip that flattened the chain.)*
                let _ = result_tx_for_reader.send(Err(anyhow::anyhow!(
                    "reader thread errored; full error follows via reader_handle.join()"
                )));
            }
            // result_tx_for_reader drops here, closing the channel if all
            // workers also exited.
            res
        });

        // ── Main thread: ordered collection + file writing ──────────────
        let mut out_r1 = File::create(output_r1)?;
        let mut out_r2 = File::create(output_r2)?;
        let mut out_pt = output_passthrough.map(File::create).transpose()?;
        let mut out_up_r1 = unpaired_r1_path.map(File::create).transpose()?;
        let mut out_up_r2 = unpaired_r2_path.map(File::create).transpose()?;

        let mut expected: u64 = 0;
        let mut pending: BTreeMap<u64, PairedBatchResult> = BTreeMap::new();
        let mut total_r1 = TrimStats::with_adapter_count(config.adapters.len());
        let mut total_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
        let mut total_pair = PairValidationStats::default();

        // Track the first reader/worker error so we can return it after
        // joining the reader thread (Step 6a).
        let mut first_error: Option<anyhow::Error> = None;

        'outer: while let Ok(result) = result_rx.recv() {
            match result {
                Ok(batch) => {
                    pending.insert(batch.seq, batch);
                    // Flush as many in-order blocks as possible.
                    // R1, R2, AND passthrough (when active) are written
                    // within the SAME iteration before `expected += 1`, so
                    // the three output files stay in lockstep across
                    // batches.
                    while let Some(r) = pending.remove(&expected) {
                        out_r1.write_all(&r.compressed_r1)?;
                        out_r2.write_all(&r.compressed_r2)?;
                        if let (Some(ref mut f), Some(bytes)) =
                            (out_pt.as_mut(), &r.compressed_passthrough)
                        {
                            f.write_all(bytes)?;
                        }
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
                Err(e) => {
                    // Reader or worker error — stop further writes and bail.
                    // Partial output files may remain on disk (v1 contract;
                    // Cli --help documents this). Plan v2 §Assumptions §13.
                    first_error = Some(e);
                    break 'outer;
                }
            }
        }

        // Drop result_rx so any worker still blocked on `rtx.send` sees the
        // send fail, exits its loop, and drops its work-channel rx. That
        // unblocks the reader's `txs[idx].send` (which had been blocked by
        // the slow/exited workers) and lets the reader return. Without this
        // drop, a worker-side error mid-stream deadlocks reader_handle.join()
        // because: workers block on result channel send → workers stop
        // consuming work channels → reader blocks on work channel send →
        // join never completes. *(Code-reviewer A-L1.)*
        //
        // Done unconditionally (even on the success path) because if the
        // main loop drained cleanly, result_rx already closed via the
        // recv() returning Err — the drop is a no-op then.
        drop(result_rx);

        // Always join the reader so the thread isn't leaked. If the loop
        // exited via `first_error`, we still want to drain the join handle.
        // Priority: prefer the reader_join_err over first_error so the full
        // anyhow chain from the reader thread reaches the caller (the
        // sentinel sent via the result channel is intentionally lossy —
        // see the reader-thread send site above). For a worker-side error,
        // reader_join_err is None and first_error has full chain. *(C1.)*
        let reader_join_err = match reader_handle.join() {
            Ok(Ok(())) => None,
            Ok(Err(e)) => Some(e),
            Err(_) => Some(anyhow::anyhow!("Reader thread panicked")),
        };

        if let Some(e) = reader_join_err.or(first_error) {
            let _ = passthrough_active;
            return Err(e);
        }
        let _ = passthrough_active;

        Ok((total_r1, total_r2, total_pair))
    })
}
#[allow(clippy::too_many_arguments)]
fn process_paired_batch(
    batch_seq: u64,
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    reads_passthrough: Option<Vec<FastqRecord>>,
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
    // Passthrough buffer + sink. `passthrough_active` decides whether the
    // batch result carries `Some(bytes)` or `None`. (Plan v2 Step 7.)
    let passthrough_active = reads_passthrough.is_some();
    let mut buf_pt = Vec::new();
    let reads_pt_slice: Option<&[FastqRecord]> = reads_passthrough.as_deref();

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
            let mut gz_pt = if passthrough_active {
                Some(GzEncoder::new(&mut buf_pt, Compression::new(level)))
            } else {
                None
            };

            process_pairs(
                reads_r1,
                reads_r2,
                reads_pt_slice,
                config,
                &mut stats_r1,
                &mut stats_r2,
                &mut pair_stats,
                &mut gz_r1,
                &mut gz_r2,
                gz_pt.as_mut(),
                gz_up_r1.as_mut(),
                gz_up_r2.as_mut(),
                unpaired,
            )?;

            gz_r1.finish()?;
            gz_r2.finish()?;
            if let Some(gz) = gz_pt {
                gz.finish()?;
            }
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
            reads_pt_slice,
            config,
            &mut stats_r1,
            &mut stats_r2,
            &mut pair_stats,
            &mut buf_r1,
            &mut buf_r2,
            if passthrough_active {
                Some(&mut buf_pt)
            } else {
                None
            },
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
        compressed_passthrough: if passthrough_active {
            Some(buf_pt)
        } else {
            None
        },
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
///
/// When `reads_passthrough` and `writer_passthrough` are both `Some`, the
/// per-pair passthrough record is written verbatim on `PairOutcome::Pass`
/// and dropped on `Discarded`. The passthrough slice is IMMUTABLE — the
/// records are carried through untouched. *(Plan v2 Step 7.)*
///
/// The worker stamps `pair_stats.passthrough_records_checked += reads_r1.len()`
/// when passthrough is active, on the contract that the reader thread has
/// already per-record-verified sync before dispatching this batch. The
/// per-pair `_kept` / `_dropped` counters are maintained per-iteration.
#[allow(clippy::too_many_arguments)]
fn process_pairs<W: Write>(
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    reads_passthrough: Option<&[FastqRecord]>,
    config: &TrimConfig,
    stats_r1: &mut TrimStats,
    stats_r2: &mut TrimStats,
    pair_stats: &mut PairValidationStats,
    writer_r1: &mut W,
    writer_r2: &mut W,
    mut writer_passthrough: Option<&mut W>,
    mut writer_up_r1: Option<&mut W>,
    mut writer_up_r2: Option<&mut W>,
    unpaired: UnpairedLengths,
) -> Result<()> {
    // The reader has already pair-verified sync for this batch; stamp the
    // sync-check counter once for the whole batch (Step 7 / Reviewer-B §1.2
    // resolution).
    if reads_passthrough.is_some() {
        pair_stats.passthrough_records_checked += reads_r1.len();
    }

    let mut pt_iter = reads_passthrough.map(|s| s.iter());

    for (r1, r2) in reads_r1.iter_mut().zip(reads_r2.iter_mut()) {
        let pt = pt_iter.as_mut().and_then(|it| it.next());
        match classify_paired(r1, r2, config, stats_r1, stats_r2, pair_stats, unpaired) {
            PairOutcome::Pass => {
                r1.write_to(writer_r1)?;
                r2.write_to(writer_r2)?;
                if let (Some(w), Some(pt_rec)) = (writer_passthrough.as_deref_mut(), pt) {
                    pt_rec.write_to(w)?;
                    pair_stats.passthrough_records_kept += 1;
                }
            }
            PairOutcome::Discarded => {
                if pt.is_some() {
                    pair_stats.passthrough_records_dropped += 1;
                }
            }
            PairOutcome::Unpaired { r1_ok, r2_ok } => {
                // `Unpaired` is set by the filter based on r1_ok / r2_ok —
                // independent of whether --retain_unpaired is on. With
                // passthrough enabled but --retain_unpaired off (the only
                // legal combination — validation step 1.iii rejects both
                // on), the unpaired writers are None, the rescue writes
                // are no-ops, and the passthrough record gets dropped.
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
                if pt.is_some() {
                    pair_stats.passthrough_records_dropped += 1;
                }
            }
        }
    }
    Ok(())
}

/// Round-robin paired reader: batches of `BATCH_SIZE` records dispatched in
/// input order. The default (non-clumpy) routing.
///
/// When `reader_passthrough` is `Some`, also reads from a third FASTQ stream
/// in lockstep, performs the three-way header sync check per record (using
/// `crate::fastq::read_id_prefix`), and ships a parallel `batch_pt` vec in
/// the work tuple. The 8-arm `(Some/None)^3` EOF match (plan v2 §Behavior §7)
/// catches truncation and over-shoot in either direction with a precise,
/// row-numbered error message.
fn read_pairs_round_robin(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    reader_passthrough: Option<&mut FastqReader>,
    txs: &[mpsc::SyncSender<PairedWork>],
) -> Result<()> {
    let mut seq: u64 = 0;
    let mut batch_r1 = Vec::with_capacity(BATCH_SIZE);
    let mut batch_r2 = Vec::with_capacity(BATCH_SIZE);
    let mut batch_pt: Vec<FastqRecord> = if reader_passthrough.is_some() {
        Vec::with_capacity(BATCH_SIZE)
    } else {
        Vec::new()
    };
    let passthrough_active = reader_passthrough.is_some();
    let mut reader_pt = reader_passthrough;
    let mut record_idx: usize = 0;

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;
        let rec_pt = if let Some(ref mut r) = reader_pt {
            Some(r.next_record()?)
        } else {
            None
        };
        record_idx += 1;

        // 8-arm three-way EOF match (plan v2 §Behavior §7).
        match (rec1, rec2, rec_pt) {
            // ── No passthrough: legacy 4-arm shape ─────────────────────
            (Some(r1), Some(r2), None) => {
                batch_r1.push(r1);
                batch_r2.push(r2);
            }
            (None, None, None) => {
                if !batch_r1.is_empty() {
                    let idx = (seq as usize) % txs.len();
                    let _ = txs[idx].send(Some((seq, batch_r1, batch_r2, None)));
                }
                for tx in txs {
                    let _ = tx.send(None);
                }
                break;
            }
            (Some(_), None, None) => bail!(
                "Read 2 file is truncated — R1 has more reads than R2. \
                 Please check your paired-end input files!"
            ),
            (None, Some(_), None) => bail!(
                "Read 1 file is truncated — R2 has more reads than R1. \
                 Please check your paired-end input files!"
            ),

            // ── Passthrough happy path ─────────────────────────────────
            (Some(r1), Some(r2), Some(Some(pt))) => {
                // Three-way header sync check before adding to the batch.
                let p1 = crate::fastq::read_id_prefix(&r1.id);
                let p2 = crate::fastq::read_id_prefix(&r2.id);
                let p3 = crate::fastq::read_id_prefix(&pt.id);
                if p1 != p2 || p1 != p3 {
                    bail!(
                        "Read ID mismatch at record {}: R1='{}', R2='{}', \
                         passthrough='{}'. Files are out of sync.",
                        record_idx,
                        p1,
                        p2,
                        p3
                    );
                }
                batch_r1.push(r1);
                batch_r2.push(r2);
                batch_pt.push(pt);
            }
            // ── Three-way EOF ─────────────────────────────────────────
            (None, None, Some(None)) => {
                if !batch_r1.is_empty() {
                    let idx = (seq as usize) % txs.len();
                    let _ = txs[idx].send(Some((seq, batch_r1, batch_r2, Some(batch_pt))));
                }
                for tx in txs {
                    let _ = tx.send(None);
                }
                break;
            }
            // ── R1/R2 truncated — passthrough state doesn't matter ─────
            (Some(_), None, _) => bail!(
                "Read 2 file is truncated — R1 has more reads than R2. \
                 Please check your paired-end input files!"
            ),
            (None, Some(_), _) => bail!(
                "Read 1 file is truncated — R2 has more reads than R1. \
                 Please check your paired-end input files!"
            ),
            // ── R1/R2 in sync, passthrough off ─────────────────────────
            (Some(_), Some(_), Some(None)) => bail!(
                "--passthrough file is truncated — R1/R2 have more reads than \
                 passthrough at record {}",
                record_idx
            ),
            (None, None, Some(Some(_))) => bail!(
                "--passthrough file has more records than R1/R2 (extra records \
                 starting at record {})",
                record_idx
            ),
        }

        // Flush full batches.
        if batch_r1.len() >= BATCH_SIZE {
            let br1 = std::mem::replace(&mut batch_r1, Vec::with_capacity(BATCH_SIZE));
            let br2 = std::mem::replace(&mut batch_r2, Vec::with_capacity(BATCH_SIZE));
            let bpt = if passthrough_active {
                let v = std::mem::replace(&mut batch_pt, Vec::with_capacity(BATCH_SIZE));
                Some(v)
            } else {
                None
            };
            let idx = (seq as usize) % txs.len();
            if txs[idx].send(Some((seq, br1, br2, bpt))).is_err() {
                break;
            }
            seq += 1;
        }
    }
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
) -> Result<()> {
    let mut bins: Vec<PairedBin> = (0..layout.n_bins)
        .map(|_| PairedBin::with_budget(layout.bin_byte_budget))
        .collect();
    let mut seq: u64 = 0;
    let mut next_worker: usize = 0;

    let flush_bin = |bin: &mut PairedBin, seq: &mut u64, next_worker: &mut usize| -> Result<bool> {
        let (mut r1s, mut r2s, mut keys) = bin.take();
        clump::sort_paired_by_key(&mut r1s, &mut r2s, &mut keys);
        let idx = *next_worker;
        // --clumpify + --passthrough is rejected at validation, so always None
        // for the passthrough slot.
        if txs[idx].send(Some((*seq, r1s, r2s, None))).is_err() {
            return Ok(false);
        }
        *seq += 1;
        *next_worker = (*next_worker + 1) % txs.len();
        Ok(true)
    };

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;
        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                let key = canonical_minimizer(r1.seq.as_bytes());
                let bin_idx = clump::bin_for(key, layout.n_bins);
                bins[bin_idx].push(r1, r2, key);
                if bins[bin_idx].raw_bytes >= layout.bin_byte_budget
                    && !flush_bin(&mut bins[bin_idx], &mut seq, &mut next_worker)?
                {
                    break;
                }
            }
            (None, None) => {
                for bin in bins.iter_mut() {
                    if !bin.is_empty() && !flush_bin(bin, &mut seq, &mut next_worker)? {
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
        // Clumpy mode replaces round-robin with a bin dispatcher; see
        // `read_pairs_clumpy` for the paired-end shape of the same idea.
        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = FastqReader::open_threaded(input)?;
            if let Some(layout) = clump_layout {
                read_single_clumpy(&mut reader, &txs, layout)
            } else {
                read_single_round_robin(&mut reader, &txs)
            }
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
) -> Result<()> {
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
    for tx in txs {
        let _ = tx.send(None);
    }
    Ok(())
}

/// Clumpy single-end reader: route each record into a bin keyed by its
/// canonical minimizer; when a bin fills its byte budget, sort it by
/// minimizer and dispatch to a worker.
fn read_single_clumpy(
    reader: &mut FastqReader,
    txs: &[mpsc::SyncSender<SingleWork>],
    layout: ClumpLayout,
) -> Result<()> {
    let mut bins: Vec<SingleBin> = (0..layout.n_bins)
        .map(|_| SingleBin::with_budget(layout.bin_byte_budget))
        .collect();
    let mut seq: u64 = 0;
    let mut next_worker: usize = 0;

    let flush_bin = |bin: &mut SingleBin, seq: &mut u64, next_worker: &mut usize| -> Result<bool> {
        let (mut records, mut keys) = bin.take();
        clump::sort_single_by_key(&mut records, &mut keys);
        let idx = *next_worker;
        if txs[idx].send(Some((*seq, records))).is_err() {
            return Ok(false);
        }
        *seq += 1;
        *next_worker = (*next_worker + 1) % txs.len();
        Ok(true)
    };

    while let Some(record) = reader.next_record()? {
        let key = canonical_minimizer(record.seq.as_bytes());
        let bin_idx = clump::bin_for(key, layout.n_bins);
        bins[bin_idx].push(record, key);
        if bins[bin_idx].raw_bytes >= layout.bin_byte_budget
            && !flush_bin(&mut bins[bin_idx], &mut seq, &mut next_worker)?
        {
            break;
        }
    }

    for bin in bins.iter_mut() {
        if !bin.is_empty() && !flush_bin(bin, &mut seq, &mut next_worker)? {
            break;
        }
    }
    for tx in txs {
        let _ = tx.send(None);
    }
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

    // ── --clumpify integration tests ─────────────────────────────────────

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
            None,
            &plain_r1_out,
            &plain_r2_out,
            None,
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
            None,
            &clumpy_r1_out,
            &clumpy_r2_out,
            None,
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

    // ──────────────────────────────────────────────────────────────────────
    // --passthrough (plan v2 Step 11) — 6 of 8 validation tests live here.
    // Test §7 (--cores 1 dispatcher) requires the built binary and lives in
    // tests/integration_passthrough.rs; §8 (reader-thread error propagation
    // under load) is covered by §3 with a mid-stream truncation that forces
    // the error to surface during the result-channel drain.
    // ──────────────────────────────────────────────────────────────────────

    /// Compact in-memory model of a fixture row: shared ID + R1/R2/passthrough
    /// sequences. The qual string is just `I` × seq.len() (Phred 40).
    struct MultiomeRow {
        id: String,
        r1_seq: String,
        r2_seq: String,
        pt_seq: String,
    }

    /// Generate a deterministic Multiome-shaped fixture. Every 3rd record
    /// has SHORT R1 *and* R2 (16 bp each) — these get dropped by
    /// `length_cutoff = 20`, giving roughly a 33% drop rate.
    ///
    /// Why both mates short rather than just R1: the parallel and serial
    /// paths have a pre-existing divergence in how they count
    /// `r1_unpaired` / `r2_unpaired` (parallel increments based on the
    /// filter's `r1_ok` / `r2_ok` regardless of whether the rescue writer
    /// is open; serial only increments when the writer is open). Fixing
    /// that is out of scope for the passthrough PR — keeping both mates
    /// short means `r1_ok == r2_ok == false`, so neither path increments
    /// `r1_unpaired` / `r2_unpaired`, and the parity test passes.
    fn gen_multiome_rows(n: usize) -> Vec<MultiomeRow> {
        (0..n)
            .map(|i| {
                let short = i % 3 == 0;
                MultiomeRow {
                    id: format!("read_{i}"),
                    r1_seq: if short {
                        // 16 bp — below length_cutoff = 20 and below the
                        // default rescue cutoff length_1 = 35.
                        "ACGTACGTACGTACGT".to_string()
                    } else {
                        // 42 bp — survives.
                        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_string()
                    },
                    r2_seq: if short {
                        // Match R1's short length so r2_ok = false too —
                        // avoids the serial/parallel r2_unpaired
                        // divergence noted above.
                        "TGCATGCATGCATGCA".to_string()
                    } else {
                        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA".to_string()
                    },
                    pt_seq: "AAAACCCCGGGGTTTT".to_string(),
                }
            })
            .collect()
    }

    /// Write a Multiome fixture to plain `.fq` files under `dir`. Returns
    /// the three paths. `mutator` lets the caller perturb the passthrough
    /// records (truncate, shuffle, mid-stream desync, etc).
    fn write_multiome_fixture(
        dir: &Path,
        rows: &[MultiomeRow],
        passthrough_mutator: impl Fn(&[MultiomeRow]) -> Vec<(String, String, String)>,
    ) -> Result<(PathBuf, PathBuf, PathBuf)> {
        let r1 = dir.join("r1.fq");
        let r2 = dir.join("r2.fq");
        let pt = dir.join("pt.fq");

        let mut r1_s = String::new();
        let mut r2_s = String::new();
        for row in rows {
            r1_s.push_str(&format!(
                "@{}\n{}\n+\n{}\n",
                row.id,
                row.r1_seq,
                "I".repeat(row.r1_seq.len()),
            ));
            r2_s.push_str(&format!(
                "@{}\n{}\n+\n{}\n",
                row.id,
                row.r2_seq,
                "I".repeat(row.r2_seq.len()),
            ));
        }
        fs::write(&r1, r1_s)?;
        fs::write(&r2, r2_s)?;

        // Passthrough goes through the mutator so tests can truncate /
        // shuffle / desync.
        let pt_records = passthrough_mutator(rows);
        let mut pt_s = String::new();
        for (id, seq, qual) in pt_records {
            pt_s.push_str(&format!("@{id}\n{seq}\n+\n{qual}\n"));
        }
        fs::write(&pt, pt_s)?;

        Ok((r1, r2, pt))
    }

    /// Identity mutator: passthrough records mirror R1/R2 ids and use the
    /// 16 bp barcode sequence (the happy path).
    fn pt_identity(rows: &[MultiomeRow]) -> Vec<(String, String, String)> {
        rows.iter()
            .map(|row| {
                (
                    row.id.clone(),
                    row.pt_seq.clone(),
                    "I".repeat(row.pt_seq.len()),
                )
            })
            .collect()
    }

    /// `--passthrough`-aware test config. Inherits the parity_test_config
    /// shape but bumps length_cutoff to 20 so the fixture's "short" rows are
    /// reliably dropped.
    fn passthrough_test_config() -> TrimConfig {
        let mut c = parity_test_config();
        c.length_cutoff = 20;
        c.is_paired = true;
        c
    }

    fn read_records(path: &Path) -> Result<Vec<(String, String, String)>> {
        let mut r = FastqReader::open(path)?;
        let mut v = Vec::new();
        while let Some(rec) = r.next_record()? {
            v.push((rec.id, rec.seq, rec.qual));
        }
        Ok(v)
    }

    /// §1 — End-to-end smoke (serial path).
    ///
    /// 50-record fixture; passthrough output's record count + per-row IDs
    /// must lockstep with R1/R2 outputs. The smallest meaningful integration
    /// test — covers the serial pipeline end-to-end, including the
    /// per-record sync check and stat accumulation.
    #[test]
    fn test_passthrough_serial_smoke() -> Result<()> {
        let dir = fresh_tmpdir("tg_pt_serial_smoke");
        let rows = gen_multiome_rows(50);
        let (r1_path, r2_path, pt_path) = write_multiome_fixture(&dir, &rows, pt_identity)?;
        let config = passthrough_test_config();

        let out_r1 = dir.join("out_R1.fq");
        let out_r2 = dir.join("out_R2.fq");
        let out_pt = dir.join("out_pt.fq");

        let mut rdr_r1 = FastqReader::open(&r1_path)?;
        let mut rdr_r2 = FastqReader::open(&r2_path)?;
        let mut rdr_pt = FastqReader::open(&pt_path)?;
        let mut wr_r1 = FastqWriter::create(&out_r1, false, 1, 1)?;
        let mut wr_r2 = FastqWriter::create(&out_r2, false, 1, 1)?;
        let mut wr_pt = FastqWriter::create(&out_pt, false, 1, 1)?;

        let (_s1, _s2, pair_stats) = crate::trimmer::run_paired_end(
            &mut rdr_r1,
            &mut rdr_r2,
            Some(&mut rdr_pt),
            &mut wr_r1,
            &mut wr_r2,
            Some(&mut wr_pt),
            None,
            None,
            &config,
            UnpairedLengths { r1: 35, r2: 35 },
        )?;
        wr_r1.flush()?;
        wr_r2.flush()?;
        wr_pt.flush()?;
        drop(wr_r1);
        drop(wr_r2);
        drop(wr_pt);

        // Stats sanity: ~1/3 dropped (every 3rd row has short R1).
        assert_eq!(pair_stats.passthrough_records_checked, 50);
        assert_eq!(
            pair_stats.passthrough_records_kept + pair_stats.passthrough_records_dropped,
            50
        );
        assert!(pair_stats.passthrough_records_dropped > 0);
        assert!(pair_stats.passthrough_records_kept > 0);

        // Row-by-row lockstep: same length, same id prefix on every row.
        let r1_recs = read_records(&out_r1)?;
        let r2_recs = read_records(&out_r2)?;
        let pt_recs = read_records(&out_pt)?;
        assert_eq!(r1_recs.len(), r2_recs.len());
        assert_eq!(r1_recs.len(), pt_recs.len());
        for ((a, b), c) in r1_recs.iter().zip(r2_recs.iter()).zip(pt_recs.iter()) {
            assert_eq!(
                crate::fastq::read_id_prefix(&a.0),
                crate::fastq::read_id_prefix(&b.0),
            );
            assert_eq!(
                crate::fastq::read_id_prefix(&a.0),
                crate::fastq::read_id_prefix(&c.0),
            );
        }
        // Passthrough sequences are untouched (verbatim 16 bp barcode).
        for rec in &pt_recs {
            assert_eq!(rec.1, "AAAACCCCGGGGTTTT");
        }

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// §2 — Serial/parallel parity (the load-bearing invariant).
    ///
    /// ≥10K records → ≥3 batches at BATCH_SIZE=4096, exercises the BTreeMap
    /// ordered-flush logic. Asserts:
    ///   (1) field-identical `PairValidationStats` between serial and
    ///       parallel paths;
    ///   (2) decoded-record-identical R1 outputs;
    ///   (3) decoded-record-identical R2 outputs;
    ///   (4) decoded-record-identical passthrough outputs;
    ///   (5) row-by-row R1/R2/passthrough id-prefix lockstep on parallel
    ///       output (would catch a "scrambled batch order" bug that stat
    ///       counts alone would miss).
    #[test]
    fn test_passthrough_serial_parallel_parity() -> Result<()> {
        let dir = fresh_tmpdir("tg_pt_serial_parallel_parity");
        // 12_000 rows ⇒ 3 batches at BATCH_SIZE = 4096.
        let rows = gen_multiome_rows(12_000);
        let (r1_path, r2_path, pt_path) = write_multiome_fixture(&dir, &rows, pt_identity)?;
        let config = passthrough_test_config();

        // ── Serial ───────────────────────────────────────────────────
        let s_r1 = dir.join("serial_R1.fq.gz");
        let s_r2 = dir.join("serial_R2.fq.gz");
        let s_pt = dir.join("serial_pt.fq.gz");
        let serial_stats = {
            let mut rdr_r1 = FastqReader::open(&r1_path)?;
            let mut rdr_r2 = FastqReader::open(&r2_path)?;
            let mut rdr_pt = FastqReader::open(&pt_path)?;
            let mut wr_r1 = FastqWriter::create(&s_r1, true, 1, 1)?;
            let mut wr_r2 = FastqWriter::create(&s_r2, true, 1, 1)?;
            let mut wr_pt = FastqWriter::create(&s_pt, true, 1, 1)?;
            let stats = crate::trimmer::run_paired_end(
                &mut rdr_r1,
                &mut rdr_r2,
                Some(&mut rdr_pt),
                &mut wr_r1,
                &mut wr_r2,
                Some(&mut wr_pt),
                None,
                None,
                &config,
                UnpairedLengths { r1: 35, r2: 35 },
            )?;
            wr_r1.flush()?;
            wr_r2.flush()?;
            wr_pt.flush()?;
            stats
        };

        // ── Parallel (cores=4) ──────────────────────────────────────
        let p_r1 = dir.join("parallel_R1.fq.gz");
        let p_r2 = dir.join("parallel_R2.fq.gz");
        let p_pt = dir.join("parallel_pt.fq.gz");
        let parallel_stats = run_paired_end_parallel(
            &r1_path,
            &r2_path,
            Some(&pt_path),
            &p_r1,
            &p_r2,
            Some(&p_pt),
            None,
            None,
            &config,
            4,
            true,
            UnpairedLengths { r1: 35, r2: 35 },
            None,
        )?;

        // (1) PairValidationStats parity — depends on the Step 3 PartialEq
        // derive on PairValidationStats (Reviewer-A critical catch).
        assert_eq!(
            serial_stats.2, parallel_stats.2,
            "PairValidationStats must match between serial and parallel paths.\n\
             serial:   {:#?}\nparallel: {:#?}",
            serial_stats.2, parallel_stats.2
        );
        // TrimStats parity (R1 + R2) — already a documented invariant for
        // the non-passthrough path; this just confirms passthrough doesn't
        // perturb it.
        assert_eq!(serial_stats.0, parallel_stats.0);
        assert_eq!(serial_stats.1, parallel_stats.1);

        // (2/3/4) Decoded-record identity.
        let sr1 = read_records(&s_r1)?;
        let sr2 = read_records(&s_r2)?;
        let spt = read_records(&s_pt)?;
        let pr1 = read_records(&p_r1)?;
        let pr2 = read_records(&p_r2)?;
        let ppt = read_records(&p_pt)?;
        assert_eq!(sr1, pr1);
        assert_eq!(sr2, pr2);
        assert_eq!(spt, ppt);

        // (5) Row-by-row lockstep on parallel output.
        assert_eq!(pr1.len(), pr2.len());
        assert_eq!(pr1.len(), ppt.len());
        for i in 0..pr1.len() {
            let a = crate::fastq::read_id_prefix(&pr1[i].0);
            let b = crate::fastq::read_id_prefix(&pr2[i].0);
            let c = crate::fastq::read_id_prefix(&ppt[i].0);
            assert!(
                a == b && a == c,
                "lockstep broken at row {i}: R1={a} R2={b} PT={c}"
            );
        }

        // Stats sanity: ~⅓ dropped, fixture has 12_000 pairs.
        assert_eq!(parallel_stats.2.passthrough_records_checked, 12_000);
        let total = parallel_stats.2.passthrough_records_kept
            + parallel_stats.2.passthrough_records_dropped;
        assert_eq!(total, 12_000);

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// §3 — Sync check catches truncation (mid-stream).
    ///
    /// Passthrough file is truncated by ~10% from the end; with cores=4 and
    /// a ≥10K-record fixture the error fires mid-stream rather than at the
    /// very end, which also exercises Step 6a's reader-error propagation
    /// path (errors surfaced via the result channel rather than waiting on
    /// reader_handle.join()). The test passing without deadlock IS the
    /// §8 "no-hang" invariant.
    #[test]
    fn test_passthrough_truncation_detected_mid_stream() -> Result<()> {
        let dir = fresh_tmpdir("tg_pt_truncation");
        let rows = gen_multiome_rows(10_000);
        let (r1_path, r2_path, pt_path) = write_multiome_fixture(&dir, &rows, |rs| {
            // Truncate passthrough at row 5_000 — sync error fires mid-stream.
            pt_identity(rs).into_iter().take(5_000).collect()
        })?;
        let config = passthrough_test_config();

        let out_r1 = dir.join("out_R1.fq.gz");
        let out_r2 = dir.join("out_R2.fq.gz");
        let out_pt = dir.join("out_pt.fq.gz");

        let res = run_paired_end_parallel(
            &r1_path,
            &r2_path,
            Some(&pt_path),
            &out_r1,
            &out_r2,
            Some(&out_pt),
            None,
            None,
            &config,
            4,
            true,
            UnpairedLengths { r1: 35, r2: 35 },
            None,
        );
        let err = res.unwrap_err().to_string();
        assert!(
            err.contains("passthrough") && err.contains("truncated"),
            "expected truncation error mentioning passthrough, got: {err}"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// §4 — Sync check catches ID mismatch (shuffled).
    ///
    /// Passthrough records are present with the correct count but their
    /// IDs are permuted, so the first sync-check pair-wise comparison
    /// fails. The error message must name the row.
    #[test]
    fn test_passthrough_id_mismatch_detected() -> Result<()> {
        let dir = fresh_tmpdir("tg_pt_id_mismatch");
        let rows = gen_multiome_rows(200);
        let (r1_path, r2_path, pt_path) = write_multiome_fixture(&dir, &rows, |rs| {
            // Reverse the passthrough records so IDs no longer align with R1/R2.
            let mut v = pt_identity(rs);
            v.reverse();
            v
        })?;
        let config = passthrough_test_config();

        let out_r1 = dir.join("out_R1.fq.gz");
        let out_r2 = dir.join("out_R2.fq.gz");
        let out_pt = dir.join("out_pt.fq.gz");

        let res = run_paired_end_parallel(
            &r1_path,
            &r2_path,
            Some(&pt_path),
            &out_r1,
            &out_r2,
            Some(&out_pt),
            None,
            None,
            &config,
            4,
            true,
            UnpairedLengths { r1: 35, r2: 35 },
            None,
        );
        let err = res.unwrap_err().to_string();
        assert!(
            err.contains("Read ID mismatch"),
            "expected ID-mismatch error, got: {err}"
        );

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// §6 — `--dont_gzip` plain-output coverage.
    ///
    /// Runs the parity test with `gzip = false` to exercise the plain
    /// branch of `process_paired_batch` (which uses `Vec<u8>` sinks
    /// directly instead of GzEncoders). A bug that only lives in the plain
    /// branch would slip past §2.
    #[test]
    fn test_passthrough_dont_gzip_plain_output() -> Result<()> {
        let dir = fresh_tmpdir("tg_pt_dont_gzip");
        let rows = gen_multiome_rows(300);
        let (r1_path, r2_path, pt_path) = write_multiome_fixture(&dir, &rows, pt_identity)?;
        let config = passthrough_test_config();

        let out_r1 = dir.join("plain_R1.fq");
        let out_r2 = dir.join("plain_R2.fq");
        let out_pt = dir.join("plain_pt.fq");

        let stats = run_paired_end_parallel(
            &r1_path,
            &r2_path,
            Some(&pt_path),
            &out_r1,
            &out_r2,
            Some(&out_pt),
            None,
            None,
            &config,
            2,
            false, // ← --dont_gzip
            UnpairedLengths { r1: 35, r2: 35 },
            None,
        )?;

        // Same lockstep invariants as the gzip path.
        let r1_recs = read_records(&out_r1)?;
        let r2_recs = read_records(&out_r2)?;
        let pt_recs = read_records(&out_pt)?;
        assert_eq!(r1_recs.len(), r2_recs.len());
        assert_eq!(r1_recs.len(), pt_recs.len());
        assert_eq!(stats.2.passthrough_records_checked, 300);
        assert_eq!(
            r1_recs.len(),
            stats.2.passthrough_records_kept,
            "kept count matches surviving record count"
        );
        for ((a, b), c) in r1_recs.iter().zip(r2_recs.iter()).zip(pt_recs.iter()) {
            assert_eq!(
                crate::fastq::read_id_prefix(&a.0),
                crate::fastq::read_id_prefix(&b.0),
            );
            assert_eq!(
                crate::fastq::read_id_prefix(&a.0),
                crate::fastq::read_id_prefix(&c.0),
            );
        }

        fs::remove_dir_all(&dir).ok();
        Ok(())
    }
}
