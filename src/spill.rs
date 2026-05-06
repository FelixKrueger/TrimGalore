//! Disk-spill clumpy mode (`--clumpy_tmp [DIR]`).
//!
//! Two-phase external bucket sort, modelled on stevekm/squish's "external"
//! engine. Phase 1 streams trimmed records into per-bin temp files (plain
//! FASTQ text, no compression). Phase 2 loads each bin into RAM, sorts it
//! by full sequence (alpha sort), and emits one gzip member per bin.
//!
//! This mode eliminates the streaming-fragmentation limit of in-memory
//! `--clumpy`: same-prefix reads from anywhere in the input land in the
//! same gzip member, regardless of when they arrived. Compression
//! results match squish/bbmap-clumpify on diverse short-read data
//! (RNA-seq, WES) where the in-memory dispatcher leaves money on the
//! table because bins flush before all their records have arrived.
//!
//! ## Fusion safety
//!
//! Temp files default to `std::env::temp_dir()` (respects `$TMPDIR`),
//! which on every Nextflow + Seqera Fusion deployment resolves to a
//! local container scratch — NOT the Fusion mount. We never default to
//! the output directory because that IS the Fusion mount on those
//! deployments and writing intermediates there would trigger spurious
//! S3 uploads on every flush. Users who set `--clumpy_tmp /custom/path`
//! take responsibility for keeping that path off Fusion.

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::mpsc;

use crate::clump;
use crate::fastq::{FastqRecord, output_gzip_level};
use crate::filters::UnpairedLengths;
use crate::parallel::{PairOutcome, SingleOutcome, classify_paired, classify_single};
use crate::report::{PairValidationStats, TrimStats};
use crate::trimmer::TrimConfig;

/// Number of bins for the disk-spill bucket sort. 256 is squish's effective
/// "buckets used" range on real data (sequence-prefix bucketing with FNV
/// hashing converges to a few dozen non-empty buckets); going higher just
/// fragments members without compression benefit.
const N_BINS: usize = 256;

/// Per-worker buffer threshold. When a worker's per-bin local buffer
/// exceeds this size, it acquires the bin's file mutex and flushes the
/// buffer in one append. Bigger threshold = less lock contention but more
/// transient memory; 256 KB is a sweet spot for typical FASTQ workloads.
const FLUSH_THRESHOLD: usize = 256 * 1024;

/// Resolve `--clumpy_tmp` value into an actual directory path. The
/// `__SYSTEM_TMP__` sentinel (set by clap when the flag is passed without
/// a value) means "use `std::env::temp_dir()`".
pub fn resolve_tmp_root(arg: &str) -> PathBuf {
    if arg == "__SYSTEM_TMP__" {
        std::env::temp_dir()
    } else {
        PathBuf::from(arg)
    }
}

/// RAII guard around the per-run temp subdir. Created with a unique name
/// (`trim_galore_clumpy_<pid>_<random>/`); removed recursively on Drop,
/// even on panic. `tempfile::TempDir` handles the cleanup contract.
pub struct SpillTempDir {
    inner: tempfile::TempDir,
}

impl SpillTempDir {
    pub fn new(root: &Path) -> Result<Self> {
        std::fs::create_dir_all(root)
            .with_context(|| format!("Could not create --clumpy_tmp root: {}", root.display()))?;
        let inner = tempfile::Builder::new()
            .prefix("trim_galore_clumpy_")
            .tempdir_in(root)
            .with_context(|| {
                format!(
                    "Could not create per-run temp subdir under {}",
                    root.display()
                )
            })?;
        Ok(Self { inner })
    }

    pub fn path(&self) -> &Path {
        self.inner.path()
    }
}

// ─── Phase 1 plumbing ─────────────────────────────────────────────────────

/// Holds the per-bin temp file handles for Phase 1, protected by mutexes
/// so multiple worker threads can append to disjoint bins concurrently.
/// Optionally also holds two unpaired files (no binning) for retain_unpaired.
struct SpillRouter {
    bin_r1: Vec<Mutex<BufWriter<File>>>,
    bin_r2: Vec<Mutex<BufWriter<File>>>,
    unpaired_r1: Option<Mutex<BufWriter<File>>>,
    unpaired_r2: Option<Mutex<BufWriter<File>>>,
    paths_r1: Vec<PathBuf>,
    paths_r2: Vec<PathBuf>,
    unpaired_path_r1: Option<PathBuf>,
    unpaired_path_r2: Option<PathBuf>,
}

impl SpillRouter {
    fn new_paired(tmp: &Path, retain_unpaired: bool) -> Result<Self> {
        let mut bin_r1 = Vec::with_capacity(N_BINS);
        let mut bin_r2 = Vec::with_capacity(N_BINS);
        let mut paths_r1 = Vec::with_capacity(N_BINS);
        let mut paths_r2 = Vec::with_capacity(N_BINS);
        for i in 0..N_BINS {
            let p1 = tmp.join(format!("bin-{i:03}-r1.fq"));
            let p2 = tmp.join(format!("bin-{i:03}-r2.fq"));
            bin_r1.push(Mutex::new(BufWriter::with_capacity(
                64 * 1024,
                File::create(&p1)?,
            )));
            bin_r2.push(Mutex::new(BufWriter::with_capacity(
                64 * 1024,
                File::create(&p2)?,
            )));
            paths_r1.push(p1);
            paths_r2.push(p2);
        }
        let (unpaired_r1, unpaired_r2, up1_path, up2_path) = if retain_unpaired {
            let p1 = tmp.join("unpaired-r1.fq");
            let p2 = tmp.join("unpaired-r2.fq");
            (
                Some(Mutex::new(BufWriter::with_capacity(
                    64 * 1024,
                    File::create(&p1)?,
                ))),
                Some(Mutex::new(BufWriter::with_capacity(
                    64 * 1024,
                    File::create(&p2)?,
                ))),
                Some(p1),
                Some(p2),
            )
        } else {
            (None, None, None, None)
        };
        Ok(Self {
            bin_r1,
            bin_r2,
            unpaired_r1,
            unpaired_r2,
            paths_r1,
            paths_r2,
            unpaired_path_r1: up1_path,
            unpaired_path_r2: up2_path,
        })
    }

    fn new_single(tmp: &Path) -> Result<Self> {
        let mut bin_r1 = Vec::with_capacity(N_BINS);
        let mut paths_r1 = Vec::with_capacity(N_BINS);
        for i in 0..N_BINS {
            let p = tmp.join(format!("bin-{i:03}.fq"));
            bin_r1.push(Mutex::new(BufWriter::with_capacity(
                64 * 1024,
                File::create(&p)?,
            )));
            paths_r1.push(p);
        }
        Ok(Self {
            bin_r1,
            bin_r2: Vec::new(),
            unpaired_r1: None,
            unpaired_r2: None,
            paths_r1,
            paths_r2: Vec::new(),
            unpaired_path_r1: None,
            unpaired_path_r2: None,
        })
    }

    /// Flush + close all per-bin file handles so Phase 2 can re-open them
    /// for reading. Workers should call this collectively at the end of
    /// Phase 1.
    fn close_all(&self) -> Result<()> {
        for w in &self.bin_r1 {
            w.lock().unwrap().flush()?;
        }
        for w in &self.bin_r2 {
            w.lock().unwrap().flush()?;
        }
        if let Some(ref w) = self.unpaired_r1 {
            w.lock().unwrap().flush()?;
        }
        if let Some(ref w) = self.unpaired_r2 {
            w.lock().unwrap().flush()?;
        }
        Ok(())
    }
}

/// Per-worker local buffers — avoid lock contention on the shared
/// per-bin file mutexes by accumulating each bin's bytes locally and only
/// acquiring the lock when the buffer crosses the flush threshold.
struct LocalBuffers {
    bin_r1: Vec<Vec<u8>>,
    bin_r2: Vec<Vec<u8>>,
    unpaired_r1: Vec<u8>,
    unpaired_r2: Vec<u8>,
}

impl LocalBuffers {
    fn new_paired() -> Self {
        Self {
            bin_r1: (0..N_BINS)
                .map(|_| Vec::with_capacity(FLUSH_THRESHOLD))
                .collect(),
            bin_r2: (0..N_BINS)
                .map(|_| Vec::with_capacity(FLUSH_THRESHOLD))
                .collect(),
            unpaired_r1: Vec::new(),
            unpaired_r2: Vec::new(),
        }
    }

    fn new_single() -> Self {
        Self {
            bin_r1: (0..N_BINS)
                .map(|_| Vec::with_capacity(FLUSH_THRESHOLD))
                .collect(),
            bin_r2: Vec::new(),
            unpaired_r1: Vec::new(),
            unpaired_r2: Vec::new(),
        }
    }

    fn flush_full(&mut self, router: &SpillRouter) -> Result<()> {
        for (i, buf) in self.bin_r1.iter_mut().enumerate() {
            if buf.len() >= FLUSH_THRESHOLD {
                router.bin_r1[i].lock().unwrap().write_all(buf)?;
                buf.clear();
            }
        }
        for (i, buf) in self.bin_r2.iter_mut().enumerate() {
            if buf.len() >= FLUSH_THRESHOLD {
                router.bin_r2[i].lock().unwrap().write_all(buf)?;
                buf.clear();
            }
        }
        if self.unpaired_r1.len() >= FLUSH_THRESHOLD
            && let Some(ref w) = router.unpaired_r1
        {
            w.lock().unwrap().write_all(&self.unpaired_r1)?;
            self.unpaired_r1.clear();
        }
        if self.unpaired_r2.len() >= FLUSH_THRESHOLD
            && let Some(ref w) = router.unpaired_r2
        {
            w.lock().unwrap().write_all(&self.unpaired_r2)?;
            self.unpaired_r2.clear();
        }
        Ok(())
    }

    fn flush_all(&mut self, router: &SpillRouter) -> Result<()> {
        for (i, buf) in self.bin_r1.iter_mut().enumerate() {
            if !buf.is_empty() {
                router.bin_r1[i].lock().unwrap().write_all(buf)?;
                buf.clear();
            }
        }
        for (i, buf) in self.bin_r2.iter_mut().enumerate() {
            if !buf.is_empty() {
                router.bin_r2[i].lock().unwrap().write_all(buf)?;
                buf.clear();
            }
        }
        if !self.unpaired_r1.is_empty()
            && let Some(ref w) = router.unpaired_r1
        {
            w.lock().unwrap().write_all(&self.unpaired_r1)?;
            self.unpaired_r1.clear();
        }
        if !self.unpaired_r2.is_empty()
            && let Some(ref w) = router.unpaired_r2
        {
            w.lock().unwrap().write_all(&self.unpaired_r2)?;
            self.unpaired_r2.clear();
        }
        Ok(())
    }
}

/// Append a record's 4-line FASTQ form to a Vec<u8>. Same byte format
/// as `FastqRecord::write_to` so Phase 2 can re-parse with `FastqReader`.
fn append_record(buf: &mut Vec<u8>, rec: &FastqRecord) {
    buf.extend_from_slice(rec.id.as_bytes());
    buf.push(b'\n');
    buf.extend_from_slice(rec.seq.as_bytes());
    buf.push(b'\n');
    buf.extend_from_slice(b"+\n");
    buf.extend_from_slice(rec.qual.as_bytes());
    buf.push(b'\n');
}

// ─── Phase 2 plumbing ─────────────────────────────────────────────────────

/// Parse a plain-FASTQ temp file fully into memory. The file was written
/// by Phase 1 in `FastqRecord::write_to` byte form.
fn read_bin_file(path: &Path) -> Result<Vec<FastqRecord>> {
    let file = File::open(path)
        .with_context(|| format!("Could not open spill bin file {}", path.display()))?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let mut out = Vec::new();
    let mut id = String::new();
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    loop {
        id.clear();
        seq.clear();
        plus.clear();
        qual.clear();
        let n = reader.read_line(&mut id)?;
        if n == 0 {
            break;
        }
        reader.read_line(&mut seq)?;
        reader.read_line(&mut plus)?;
        reader.read_line(&mut qual)?;
        out.push(FastqRecord {
            id: id.trim_end_matches('\n').to_string(),
            seq: seq.trim_end_matches('\n').to_string(),
            qual: qual.trim_end_matches('\n').to_string(),
        });
    }
    Ok(out)
}

/// Encode a Vec of trimmed records as ONE gzip member, append to the
/// open output file. `gzip` controls whether output is gzip-encoded
/// (always true under --clumpy because --dont_gzip is rejected at
/// validation).
fn emit_gzip_member(out_file: &mut File, records: &[FastqRecord], gzip_level: u32) -> Result<()> {
    let mut buf = Vec::with_capacity(
        records
            .iter()
            .map(|r| r.id.len() + r.seq.len() + r.qual.len() + 5)
            .sum::<usize>(),
    );
    {
        let mut gz = GzEncoder::new(&mut buf, Compression::new(gzip_level));
        for r in records {
            r.write_to(&mut gz)?;
        }
        gz.finish()?;
    }
    out_file.write_all(&buf)?;
    Ok(())
}

// ─── Top-level paired-end disk-spill runner ─────────────────────────────

/// Run paired-end clumpy with disk-spill external sort. Two phases:
/// Phase 1 streams trimmed records to per-bin temp files; Phase 2 loads
/// each bin, alpha-sorts, and emits one gzip member per bin.
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end_clumpy_spill(
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
    unpaired_r1_path: Option<&Path>,
    unpaired_r2_path: Option<&Path>,
    config: &TrimConfig,
    cores: usize,
    unpaired: UnpairedLengths,
    tmp_dir: &Path,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    use crate::fastq::FastqReader;

    let retain_unpaired = unpaired_r1_path.is_some();
    let router = SpillRouter::new_paired(tmp_dir, retain_unpaired)?;
    let gzip_level = output_gzip_level(config.high_compression);

    eprintln!("clumpy disk-spill: tmp dir = {}", tmp_dir.display());
    eprintln!(
        "clumpy disk-spill: Phase 1 — trimming and routing to {} bins…",
        N_BINS
    );

    // ── Phase 1: parallel trim + route ────────────────────────────────
    type PairedWork = Option<(Vec<FastqRecord>, Vec<FastqRecord>)>;

    let mut work_txs: Vec<mpsc::SyncSender<PairedWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<PairedWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<PairedWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }
    let (stats_tx, stats_rx) = mpsc::channel::<(TrimStats, TrimStats, PairValidationStats)>();

    std::thread::scope(|s| -> Result<()> {
        let mut worker_handles = Vec::with_capacity(cores);
        for rx in work_rxs.drain(..) {
            let stx = stats_tx.clone();
            let router_ref = &router;
            let handle = s.spawn(move || -> Result<()> {
                let mut local = LocalBuffers::new_paired();
                let mut stats_r1 = TrimStats::with_adapter_count(config.adapters.len());
                let mut stats_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
                let mut pair_stats = PairValidationStats::default();
                while let Ok(Some((mut r1s, mut r2s))) = rx.recv() {
                    for (r1, r2) in r1s.iter_mut().zip(r2s.iter_mut()) {
                        match classify_paired(
                            r1,
                            r2,
                            config,
                            &mut stats_r1,
                            &mut stats_r2,
                            &mut pair_stats,
                            unpaired,
                        ) {
                            PairOutcome::Pass => {
                                let bin = clump::bin_for_seq(r1.seq.as_bytes(), N_BINS);
                                append_record(&mut local.bin_r1[bin], r1);
                                append_record(&mut local.bin_r2[bin], r2);
                            }
                            PairOutcome::Discarded => {}
                            PairOutcome::Unpaired { r1_ok, r2_ok } => {
                                if r1_ok && retain_unpaired {
                                    append_record(&mut local.unpaired_r1, r1);
                                }
                                if r2_ok && retain_unpaired {
                                    append_record(&mut local.unpaired_r2, r2);
                                }
                            }
                        }
                    }
                    local.flush_full(router_ref)?;
                }
                local.flush_all(router_ref)?;
                let _ = stx.send((stats_r1, stats_r2, pair_stats));
                Ok(())
            });
            worker_handles.push(handle);
        }
        drop(stats_tx);

        // Reader thread
        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            const BATCH_SIZE: usize = 4096;
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
                            if txs[idx].send(Some((br1, br2))).is_err() {
                                break;
                            }
                            seq += 1;
                        }
                    }
                    (None, None) => {
                        if !batch_r1.is_empty() {
                            let idx = (seq as usize) % txs.len();
                            let _ = txs[idx].send(Some((batch_r1, batch_r2)));
                        }
                        for tx in &txs {
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
        });

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => bail!("Reader thread panicked"),
        }
        // Propagate any worker error (e.g. disk-full during Phase 1).
        // Without this collection step, errors returned from worker
        // closures are silently lost — the original bug that produced
        // empty output files when /tmp filled up. A worker that errors
        // out also stops sending stats, so the run would otherwise
        // silently report partial counts.
        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => return Err(e),
                Err(_) => bail!("Worker thread panicked"),
            }
        }
        Ok(())
    })?;

    // Aggregate worker stats
    let mut total_r1 = TrimStats::with_adapter_count(config.adapters.len());
    let mut total_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
    let mut total_pair = PairValidationStats::default();
    while let Ok((s1, s2, sp)) = stats_rx.recv() {
        total_r1.merge(&s1);
        total_r2.merge(&s2);
        total_pair.merge(&sp);
    }

    router.close_all()?;

    // ── Phase 2: per-bin sort + emit ───────────────────────────────────
    // Parallelized across `cores` workers. Each worker reads one bin from
    // disk, alpha-sorts, gzip-encodes both R1 and R2 into `Vec<u8>`s, and
    // sends them back with the bin index. Main thread reassembles in
    // bin order via BTreeMap, writes to output, and deletes the temp
    // files once a bin is consumed (so disk peak stays below 1× input
    // during Phase 2 instead of growing).
    eprintln!(
        "clumpy disk-spill: Phase 2 — alpha-sorting and emitting {} bins (parallel, {} workers)…",
        N_BINS, cores
    );
    let mut out_r1 = File::create(output_r1)?;
    let mut out_r2 = File::create(output_r2)?;

    type PhaseTwoResult = (usize, Vec<u8>, Vec<u8>);
    // Snapshot which bins are non-empty BEFORE spawning workers — workers
    // delete the temp files after compressing, so a metadata check from
    // the main thread mid-flight would report "missing" and we'd
    // incorrectly skip the bin during reassembly. The snapshot avoids
    // racing the workers.
    let active_bins: Vec<usize> = (0..N_BINS)
        .filter(|&b| {
            std::fs::metadata(&router.paths_r1[b])
                .map(|m| m.len() > 0)
                .unwrap_or(false)
        })
        .collect();
    let (job_tx, job_rx) = mpsc::channel::<usize>();
    let job_rx = std::sync::Mutex::new(job_rx);
    let (result_tx, result_rx) = mpsc::sync_channel::<Result<PhaseTwoResult>>(cores * 2);

    for &bin in &active_bins {
        let _ = job_tx.send(bin);
    }
    drop(job_tx);

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..cores {
            let rtx = result_tx.clone();
            let rxr = &job_rx;
            let router_ref = &router;
            s.spawn(move || {
                loop {
                    let bin = {
                        let lock = rxr.lock().unwrap();
                        match lock.recv() {
                            Ok(b) => b,
                            Err(_) => return,
                        }
                    };
                    let res: Result<PhaseTwoResult> = (|| {
                        let p1 = &router_ref.paths_r1[bin];
                        let p2 = &router_ref.paths_r2[bin];
                        let mut r1 = read_bin_file(p1)?;
                        let mut r2 = read_bin_file(p2)?;
                        if r1.len() != r2.len() {
                            bail!(
                                "internal error: spill bin {} has mismatched R1/R2 lengths ({} vs {})",
                                bin, r1.len(), r2.len()
                            );
                        }
                        clump::sort_paired_alpha(&mut r1, &mut r2);
                        let approx_size: usize = r1
                            .iter()
                            .map(|r| r.id.len() + r.seq.len() + r.qual.len() + 5)
                            .sum();
                        let mut buf_r1 = Vec::with_capacity(approx_size);
                        let mut buf_r2 = Vec::with_capacity(approx_size);
                        {
                            let mut gz = GzEncoder::new(&mut buf_r1, Compression::new(gzip_level));
                            for r in &r1 {
                                r.write_to(&mut gz)?;
                            }
                            gz.finish()?;
                        }
                        {
                            let mut gz = GzEncoder::new(&mut buf_r2, Compression::new(gzip_level));
                            for r in &r2 {
                                r.write_to(&mut gz)?;
                            }
                            gz.finish()?;
                        }
                        // Free temp files eagerly so disk doesn't fill up during Phase 2.
                        let _ = std::fs::remove_file(p1);
                        let _ = std::fs::remove_file(p2);
                        Ok((bin, buf_r1, buf_r2))
                    })();
                    if rtx.send(res).is_err() {
                        return;
                    }
                }
            });
        }
        drop(result_tx);

        // Reassemble in bin order. Workers may finish out of order; we
        // hold compressed bytes in a BTreeMap until the next expected bin
        // (in `active_bins` order) arrives.
        use std::collections::BTreeMap;
        let mut pending: BTreeMap<usize, (Vec<u8>, Vec<u8>)> = BTreeMap::new();
        let mut next_idx = 0usize;
        while let Ok(item) = result_rx.recv() {
            let (bin, c1, c2) = item?;
            pending.insert(bin, (c1, c2));
            while next_idx < active_bins.len() {
                let expected = active_bins[next_idx];
                if let Some((c1, c2)) = pending.remove(&expected) {
                    out_r1.write_all(&c1)?;
                    out_r2.write_all(&c2)?;
                    next_idx += 1;
                } else {
                    break;
                }
            }
        }
        if !pending.is_empty() || next_idx != active_bins.len() {
            bail!(
                "internal error: Phase 2 reassembly mismatch (next_idx={}, active_bins.len()={}, pending={})",
                next_idx,
                active_bins.len(),
                pending.len()
            );
        }
        Ok(())
    })?;

    if retain_unpaired {
        if let (Some(in_path), Some(out_path)) = (&router.unpaired_path_r1, unpaired_r1_path) {
            emit_unpaired(in_path, out_path, gzip_level)?;
        }
        if let (Some(in_path), Some(out_path)) = (&router.unpaired_path_r2, unpaired_r2_path) {
            emit_unpaired(in_path, out_path, gzip_level)?;
        }
    }

    Ok((total_r1, total_r2, total_pair))
}

fn emit_unpaired(in_path: &Path, out_path: &Path, gzip_level: u32) -> Result<()> {
    if std::fs::metadata(in_path)
        .map(|m| m.len() == 0)
        .unwrap_or(true)
    {
        // Still create an empty output file for downstream tools.
        File::create(out_path)?;
        return Ok(());
    }
    let mut recs = read_bin_file(in_path)?;
    clump::sort_single_alpha(&mut recs);
    let mut out = File::create(out_path)?;
    emit_gzip_member(&mut out, &recs, gzip_level)?;
    let _ = std::fs::remove_file(in_path);
    Ok(())
}

// ─── Top-level single-end disk-spill runner ─────────────────────────────

pub fn run_single_end_clumpy_spill(
    input: &Path,
    output: &Path,
    config: &TrimConfig,
    cores: usize,
    tmp_dir: &Path,
) -> Result<TrimStats> {
    use crate::fastq::FastqReader;

    let router = SpillRouter::new_single(tmp_dir)?;
    let gzip_level = output_gzip_level(config.high_compression);

    eprintln!("clumpy disk-spill: tmp dir = {}", tmp_dir.display());
    eprintln!(
        "clumpy disk-spill: Phase 1 — trimming and routing to {} bins…",
        N_BINS
    );

    type SingleWork = Option<Vec<FastqRecord>>;
    let mut work_txs: Vec<mpsc::SyncSender<SingleWork>> = Vec::with_capacity(cores);
    let mut work_rxs: Vec<mpsc::Receiver<SingleWork>> = Vec::with_capacity(cores);
    for _ in 0..cores {
        let (tx, rx) = mpsc::sync_channel::<SingleWork>(2);
        work_txs.push(tx);
        work_rxs.push(rx);
    }
    let (stats_tx, stats_rx) = mpsc::channel::<TrimStats>();

    std::thread::scope(|s| -> Result<()> {
        let mut worker_handles = Vec::with_capacity(cores);
        for rx in work_rxs.drain(..) {
            let stx = stats_tx.clone();
            let router_ref = &router;
            let handle = s.spawn(move || -> Result<()> {
                let mut local = LocalBuffers::new_single();
                let mut stats = TrimStats::with_adapter_count(config.adapters.len());
                while let Ok(Some(mut reads)) = rx.recv() {
                    for record in reads.iter_mut() {
                        if let SingleOutcome::Pass = classify_single(record, config, &mut stats) {
                            let bin = clump::bin_for_seq(record.seq.as_bytes(), N_BINS);
                            append_record(&mut local.bin_r1[bin], record);
                        }
                    }
                    local.flush_full(router_ref)?;
                }
                local.flush_all(router_ref)?;
                let _ = stx.send(stats);
                Ok(())
            });
            worker_handles.push(handle);
        }
        drop(stats_tx);

        let txs = std::mem::take(&mut work_txs);
        let reader_handle = s.spawn(move || -> Result<()> {
            const BATCH_SIZE: usize = 4096;
            let mut reader = FastqReader::open_threaded(input)?;
            let mut seq: u64 = 0;
            let mut batch = Vec::with_capacity(BATCH_SIZE);
            while let Some(record) = reader.next_record()? {
                batch.push(record);
                if batch.len() >= BATCH_SIZE {
                    let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
                    let idx = (seq as usize) % txs.len();
                    if txs[idx].send(Some(full)).is_err() {
                        break;
                    }
                    seq += 1;
                }
            }
            if !batch.is_empty() {
                let idx = (seq as usize) % txs.len();
                let _ = txs[idx].send(Some(batch));
            }
            for tx in &txs {
                let _ = tx.send(None);
            }
            Ok(())
        });
        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => bail!("Reader thread panicked"),
        }
        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => return Err(e),
                Err(_) => bail!("Worker thread panicked"),
            }
        }
        Ok(())
    })?;

    let mut total = TrimStats::with_adapter_count(config.adapters.len());
    while let Ok(s) = stats_rx.recv() {
        total.merge(&s);
    }

    router.close_all()?;

    eprintln!(
        "clumpy disk-spill: Phase 2 — alpha-sorting and emitting {} bins (parallel, {} workers)…",
        N_BINS, cores
    );
    let mut out = File::create(output)?;
    type SinglePhaseTwo = (usize, Vec<u8>);
    let active_bins: Vec<usize> = (0..N_BINS)
        .filter(|&b| {
            std::fs::metadata(&router.paths_r1[b])
                .map(|m| m.len() > 0)
                .unwrap_or(false)
        })
        .collect();
    let (job_tx, job_rx) = mpsc::channel::<usize>();
    let job_rx = std::sync::Mutex::new(job_rx);
    let (result_tx, result_rx) = mpsc::sync_channel::<Result<SinglePhaseTwo>>(cores * 2);

    for &bin in &active_bins {
        let _ = job_tx.send(bin);
    }
    drop(job_tx);

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..cores {
            let rtx = result_tx.clone();
            let rxr = &job_rx;
            let router_ref = &router;
            s.spawn(move || {
                loop {
                    let bin = {
                        let lock = rxr.lock().unwrap();
                        match lock.recv() {
                            Ok(b) => b,
                            Err(_) => return,
                        }
                    };
                    let res: Result<SinglePhaseTwo> = (|| {
                        let p = &router_ref.paths_r1[bin];
                        let mut recs = read_bin_file(p)?;
                        clump::sort_single_alpha(&mut recs);
                        let approx_size: usize = recs
                            .iter()
                            .map(|r| r.id.len() + r.seq.len() + r.qual.len() + 5)
                            .sum();
                        let mut buf = Vec::with_capacity(approx_size);
                        {
                            let mut gz = GzEncoder::new(&mut buf, Compression::new(gzip_level));
                            for r in &recs {
                                r.write_to(&mut gz)?;
                            }
                            gz.finish()?;
                        }
                        let _ = std::fs::remove_file(p);
                        Ok((bin, buf))
                    })();
                    if rtx.send(res).is_err() {
                        return;
                    }
                }
            });
        }
        drop(result_tx);

        use std::collections::BTreeMap;
        let mut pending: BTreeMap<usize, Vec<u8>> = BTreeMap::new();
        let mut next_idx = 0usize;
        while let Ok(item) = result_rx.recv() {
            let (bin, c) = item?;
            pending.insert(bin, c);
            while next_idx < active_bins.len() {
                let expected = active_bins[next_idx];
                if let Some(c) = pending.remove(&expected) {
                    out.write_all(&c)?;
                    next_idx += 1;
                } else {
                    break;
                }
            }
        }
        if !pending.is_empty() || next_idx != active_bins.len() {
            bail!(
                "internal error: Phase 2 reassembly mismatch (next_idx={}, active_bins.len()={}, pending={})",
                next_idx,
                active_bins.len(),
                pending.len()
            );
        }
        Ok(())
    })?;

    Ok(total)
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_tmp_root_sentinel_uses_system_temp() {
        let p = resolve_tmp_root("__SYSTEM_TMP__");
        assert_eq!(p, std::env::temp_dir());
    }

    #[test]
    fn test_resolve_tmp_root_user_path() {
        let p = resolve_tmp_root("/some/custom/path");
        assert_eq!(p, PathBuf::from("/some/custom/path"));
    }

    #[test]
    fn test_spill_temp_dir_cleans_up_on_drop() {
        let root = std::env::temp_dir();
        let dir_path = {
            let guard = SpillTempDir::new(&root).expect("create temp dir");
            let p = guard.path().to_path_buf();
            assert!(p.exists());
            assert!(p.starts_with(&root));
            p
        };
        // After Drop, the dir is gone.
        assert!(
            !dir_path.exists(),
            "temp dir not cleaned up: {:?}",
            dir_path
        );
    }

    #[test]
    fn test_append_record_matches_write_to() {
        let r = FastqRecord {
            id: "@test/1".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        let mut via_append = Vec::new();
        append_record(&mut via_append, &r);
        let mut via_write = Vec::new();
        r.write_to(&mut via_write).unwrap();
        assert_eq!(via_append, via_write);
    }

    #[test]
    fn test_read_bin_file_round_trips() -> Result<()> {
        let tmp = SpillTempDir::new(&std::env::temp_dir())?;
        let path = tmp.path().join("test.fq");
        let recs = vec![
            FastqRecord {
                id: "@r1".to_string(),
                seq: "AAAA".to_string(),
                qual: "BBBB".to_string(),
            },
            FastqRecord {
                id: "@r2".to_string(),
                seq: "CCCC".to_string(),
                qual: "DDDD".to_string(),
            },
        ];
        let mut buf = Vec::new();
        for r in &recs {
            append_record(&mut buf, r);
        }
        std::fs::write(&path, &buf)?;
        let back = read_bin_file(&path)?;
        assert_eq!(back.len(), recs.len());
        for (a, b) in back.iter().zip(recs.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.seq, b.seq);
            assert_eq!(a.qual, b.qual);
        }
        Ok(())
    }

    // ── End-to-end pipeline tests ───────────────────────────────────────

    use crate::fastq::FastqReader;
    use crate::trimmer::TrimConfig;

    fn fresh_tmpdir(slug: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

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
            high_compression: false,
        }
    }

    fn read_records(path: &Path) -> Result<Vec<(String, String, String)>> {
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

    /// Disk-spill single-end output is a permutation of the plain output —
    /// same multiset of records, possibly different order on disk.
    #[test]
    fn test_disk_spill_single_end_is_a_permutation() -> Result<()> {
        let dir = fresh_tmpdir("tg_spill_se_permutation");
        let input_path = dir.join("input.fq");

        // ~6K reads spanning multiple bins.
        let mut content = String::new();
        let prefixes = ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACGTACGT"];
        for i in 0..6_000_u32 {
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
        let spill_out = dir.join("spill.fq.gz");
        let spill_tmp = SpillTempDir::new(&std::env::temp_dir())?;

        crate::parallel::run_single_end_parallel(&input_path, &plain_out, &config, 4, true, None)?;
        run_single_end_clumpy_spill(&input_path, &spill_out, &config, 4, spill_tmp.path())?;

        let plain = read_records(&plain_out)?;
        let spill = read_records(&spill_out)?;
        assert_eq!(plain.len(), spill.len());
        assert_eq!(
            sorted(plain),
            sorted(spill),
            "spill output must be a permutation of plain"
        );

        // Tmp dir was used during the run; after Drop it's cleaned up.
        let tmp_path = spill_tmp.path().to_path_buf();
        drop(spill_tmp);
        assert!(!tmp_path.exists(), "SpillTempDir must clean up on Drop");

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Disk-spill paired output preserves pair lockstep — R1[i]/R2[i] are
    /// still mates after the bin-route + alpha-sort.
    #[test]
    fn test_disk_spill_paired_preserves_pair_lockstep() -> Result<()> {
        let dir = fresh_tmpdir("tg_spill_pe_lockstep");
        let r1_path = dir.join("r1.fq");
        let r2_path = dir.join("r2.fq");

        let prefixes = ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACGTACGT"];
        let mut r1c = String::new();
        let mut r2c = String::new();
        for i in 0..3_000_u32 {
            let prefix = prefixes[(i as usize) % prefixes.len()];
            let r1_tail: String = (0..40)
                .map(|j| {
                    let r = i.wrapping_mul(2654435761).wrapping_add(j as u32);
                    [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                })
                .collect();
            let r2_tail: String = (0..40)
                .map(|j| {
                    let r = i.wrapping_mul(0x9E3779B1).wrapping_add(j as u32);
                    [b'A', b'C', b'G', b'T'][(r >> 28) as usize & 0x3] as char
                })
                .collect();
            r1c.push_str(&format!(
                "@pair_{i}/1\n{prefix}{r1_tail}\n+\n{}\n",
                "I".repeat(prefix.len() + r1_tail.len())
            ));
            r2c.push_str(&format!(
                "@pair_{i}/2\n{prefix}{r2_tail}\n+\n{}\n",
                "I".repeat(prefix.len() + r2_tail.len())
            ));
        }
        std::fs::write(&r1_path, r1c)?;
        std::fs::write(&r2_path, r2c)?;

        let mut config = parity_test_config();
        config.is_paired = true;
        let r1_out = dir.join("spill_r1.fq.gz");
        let r2_out = dir.join("spill_r2.fq.gz");
        let spill_tmp = SpillTempDir::new(&std::env::temp_dir())?;

        run_paired_end_clumpy_spill(
            &r1_path,
            &r2_path,
            &r1_out,
            &r2_out,
            None,
            None,
            &config,
            4,
            UnpairedLengths { r1: 35, r2: 35 },
            spill_tmp.path(),
        )?;

        let r1 = read_records(&r1_out)?;
        let r2 = read_records(&r2_out)?;
        assert_eq!(r1.len(), r2.len());
        for (a, b) in r1.iter().zip(r2.iter()) {
            let stem_a = a.0.trim_end_matches("/1");
            let stem_b = b.0.trim_end_matches("/2");
            assert_eq!(stem_a, stem_b, "pair lockstep broken: {} ≠ {}", a.0, b.0);
        }

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Disk-spill TrimStats are equal to the plain-pipeline TrimStats —
    /// the trim+filter logic is shared via `classify_single`, so this is
    /// the regression test that the refactor did not drift.
    #[test]
    fn test_disk_spill_stats_parity_with_plain() -> Result<()> {
        let dir = fresh_tmpdir("tg_spill_stats_parity");
        let input_path = dir.join("input.fq");
        let mut content = String::new();
        for i in 0..2_000_u32 {
            let seq = if i % 3 == 0 {
                "ACGTACGTACGTACGTACGTACGTACGTACAGATCGGAAGAGC"
            } else {
                "ACGTACGTACGTACGTACGTACGTACGTACGTAC"
            };
            content.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
        }
        std::fs::write(&input_path, content)?;

        let config = parity_test_config();
        let plain_out = dir.join("plain.fq.gz");
        let spill_out = dir.join("spill.fq.gz");
        let spill_tmp = SpillTempDir::new(&std::env::temp_dir())?;

        let plain_stats = crate::parallel::run_single_end_parallel(
            &input_path,
            &plain_out,
            &config,
            4,
            true,
            None,
        )?;
        let spill_stats =
            run_single_end_clumpy_spill(&input_path, &spill_out, &config, 4, spill_tmp.path())?;
        assert_eq!(plain_stats, spill_stats, "TrimStats must be identical");
        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    /// Tiny input (<< 1 bin's worth) still produces a valid round-trippable
    /// output. Exercises the empty-bin skip path in Phase 2.
    #[test]
    fn test_disk_spill_tiny_input() -> Result<()> {
        let dir = fresh_tmpdir("tg_spill_tiny");
        let input_path = dir.join("input.fq");
        let mut content = String::new();
        for i in 0..5 {
            content.push_str(&format!(
                "@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        std::fs::write(&input_path, content)?;

        let config = parity_test_config();
        let out = dir.join("out.fq.gz");
        let spill_tmp = SpillTempDir::new(&std::env::temp_dir())?;
        let stats = run_single_end_clumpy_spill(&input_path, &out, &config, 2, spill_tmp.path())?;
        assert_eq!(stats.total_reads, 5);
        assert_eq!(stats.reads_written, 5);

        // Multi-member gzip round-trip (no truncation).
        let mut count = 0;
        let mut r = FastqReader::open(&out)?;
        while r.next_record()?.is_some() {
            count += 1;
        }
        assert_eq!(count, 5);
        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }
}
