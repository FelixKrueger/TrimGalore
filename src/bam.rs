//! uBAM (unaligned BAM) input reader.
//!
//! Converts a `noodles::bam` record stream into the `FastqRecord` shape the
//! rest of the trimming pipeline expects. See `plans/06252026_ubam-input-support/PLAN.md`
//! §3.2 for the record-conversion contract.

use anyhow::{Context, Result, bail};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
use noodles::sam::header::record::value::{
    Map,
    map::{self, Program, header::Version, program::tag as program_tag},
};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::mpsc;
use std::thread::JoinHandle;

use crate::fastq::{FastqRecord, RecordSource};

impl RecordSource for BamReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        BamReader::next_record(self)
    }
}

/// Records-per-batch for the threaded reader channel. Matches the
/// `FastqReader` constant in `src/fastq.rs` (`READER_BATCH_SIZE`).
const BAM_READER_BATCH_SIZE: usize = 4096;
/// Lookahead batches buffered ahead in the channel. Matches `FastqReader`.
const BAM_READER_CHANNEL_BATCHES: usize = 4;

/// Per-side queue depth bound for the paired-interleaved de-interleaver.
/// Spike 2 surveyed all standard uBAM-emitting tools (samtools sort -n /
/// collate, Picard `FastqToSam`, fgbio `FastqToBam`) and found mate-adjacent
/// ordering universally. The slack is a guardrail for bespoke hand-spliced
/// pipelines that violate the contract — at ~1 KB / record it caps the
/// per-side buffer at ~1 MB, not unbounded.
///
/// See `plans/06252026_ubam-input-support/spikes/SPIKE_paired_ordering.md`.
const MAX_SLACK: usize = 1024;

/// Suggested remediation text when the de-interleaver detects grouped input.
const GROUPED_INPUT_ERR: &str = "Paired uBAM input has non-adjacent mates (one side queue exceeded MAX_SLACK = 1024 \
     records without seeing a matching mate). TrimGalore's paired-uBAM reader requires \
     mate-adjacent ordering to stream without unbounded buffering. \
     Re-interleave the input first: \
     `samtools collate -O input.bam tmp > interleaved.bam` \
     (or `samtools sort -n input.bam -o interleaved.bam`).";

/// Offset the **reader** adds when converting BAM raw Phred to Sanger ASCII.
///
/// Fixed at 33 and must stay fixed, but note precisely why: BAM `QUAL` is always
/// raw Phred, so there is no input-encoding variation to accommodate on the read
/// side; 33 is the right constant because TrimGalore's internal `FastqRecord` is
/// Sanger-encoded. (The SAM spec's ASCII+33 rule governs SAM *text*; BAM itself
/// stores raw bytes.)
///
/// Not to be confused with [`BamWriter::input_phred_offset`], which is 33 *or*
/// 64 and is the offset the **writer subtracts** from incoming FASTQ ASCII
/// (issue #358). Conflating the two is the defect that issue fixed.
const PHRED_OFFSET: u8 = 33;
/// Highest raw Phred score the SAM spec permits in `QUAL` (`93 + 33 = 126`,
/// i.e. `'~'`, the last printable ASCII character). Out-of-spec bytes are
/// clamped to this on read so the ASCII conversion stays single-byte — see
/// `bam_record_to_fastq`.
const MAX_PHRED_SCORE: u8 = 93;
const QUAL_MISSING_SENTINEL: u8 = 0xFF;
const QUAL_MISSING_REPLACEMENT: u8 = b'!';

/// Internal source for a `BamReader` — either direct (synchronous) reads or
/// a channel fed by a background decompression+parsing thread, mirroring
/// `crate::fastq::ReaderSource`.
enum BamReaderSource {
    Direct {
        inner: bam::io::Reader<bgzf::Reader<BufReader<File>>>,
        path: PathBuf,
        record_idx: usize,
    },
    Threaded {
        rx: mpsc::Receiver<Result<Vec<FastqRecord>>>,
        buffer: Vec<FastqRecord>,
        buf_pos: usize,
        _handle: JoinHandle<()>,
    },
}

/// uBAM reader. Yields `FastqRecord`s via `next_record()`, mirroring the
/// `FastqReader` shape so the worker pool can dispatch through the
/// `RecordSource` trait object.
///
/// All BAM records MUST have `BAM_FUNMAP` (`0x4`) set — aligned BAM input is
/// rejected per-record (see `bam_record_to_fastq` flag check). The first-record
/// fast-path check happens in `main.rs::sanity_check_any` (Task 19).
pub struct BamReader {
    source: BamReaderSource,
    preserved_tags: Vec<String>,
}

impl BamReader {
    /// Open a uBAM file for synchronous reading. Reads the BAM header
    /// eagerly; errors on malformed header.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let (inner, _header) = open_inner(&path)?;
        Ok(Self {
            source: BamReaderSource::Direct {
                inner,
                path,
                record_idx: 0,
            },
            preserved_tags: Vec::new(),
        })
    }

    /// Open a uBAM file with background decompression+parsing on a dedicated
    /// thread, mirroring `FastqReader::open_threaded`. Returns immediately;
    /// the thread batches `FastqRecord`s through a bounded channel.
    pub fn open_threaded<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::open_threaded_inner(path, Vec::new())
    }

    /// `open_threaded` plus tag preservation in a single call. The
    /// background thread snapshots `preserve_tags` at construction, so the
    /// builder-style `with_preserved_tags` would NOT be honoured after
    /// `open_threaded` — use this constructor instead when both are needed.
    pub fn open_threaded_with_tags<P: AsRef<Path>>(
        path: P,
        preserve_tags: &[String],
    ) -> Result<Self> {
        Self::open_threaded_inner(path, preserve_tags.to_vec())
    }

    /// Open a paired uBAM and de-interleave on the fly.
    ///
    /// Returns two `BamReader`s, both in `Threaded` mode, fed from a single
    /// background thread that:
    /// - reads BAM records sequentially from `path`,
    /// - routes by `BAM_FREAD1` (0x40) / `BAM_FREAD2` (0x80) flag bits,
    /// - emits matched pairs (R1, R2) to two parallel channels in lockstep,
    /// - bounds per-side buffering at `MAX_SLACK` records and errors with
    ///   `GROUPED_INPUT_ERR` if exceeded (catches `cat r1.bam r2.bam`-style
    ///   hand-spliced grouped input without exhausting memory).
    ///
    /// Per Spike 2: all standard tools (samtools sort -n / collate, Picard,
    /// fgbio) emit mate-adjacent paired BAM, so the happy path is `MAX_SLACK
    /// = 0` in practice. The bound is a guardrail.
    pub fn open_paired_interleaved<P: AsRef<Path>>(path: P) -> Result<(Self, Self)> {
        Self::open_paired_interleaved_with_tags(path, &[])
    }

    /// `open_paired_interleaved` plus tag preservation in a single call.
    pub fn open_paired_interleaved_with_tags<P: AsRef<Path>>(
        path: P,
        preserve_tags: &[String],
    ) -> Result<(Self, Self)> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            bail!("Input file not found: {}", path.display());
        }
        let (tx_r1, rx_r1) = mpsc::sync_channel(BAM_READER_CHANNEL_BATCHES);
        let (tx_r2, rx_r2) = mpsc::sync_channel(BAM_READER_CHANNEL_BATCHES);
        let tags_for_thread = preserve_tags.to_vec();
        let path_for_thread = path.clone();

        // Producer thread — the JoinHandle drops at end of this fn, which
        // detaches the thread (per Rust semantics; the thread runs to
        // completion regardless). Errors are surfaced to consumers via the
        // channels, so we don't need to keep the handle around for joining.
        let _producer_handle: JoinHandle<()> = std::thread::spawn(move || {
            let mut inner = match open_inner(&path_for_thread) {
                Ok((r, _header)) => r,
                Err(e) => {
                    let msg = format!("{:#}", e);
                    let _ = tx_r1.send(Err(anyhow::anyhow!("{msg}")));
                    let _ = tx_r2.send(Err(anyhow::anyhow!("{msg}")));
                    return;
                }
            };
            // Per-side queues for the de-interleaver. Records buffered here
            // are FastqRecord (already-converted) so the conversion cost is
            // accounted for at-read-time, not at-emit-time.
            let mut r1_queue: std::collections::VecDeque<FastqRecord> =
                std::collections::VecDeque::new();
            let mut r2_queue: std::collections::VecDeque<FastqRecord> =
                std::collections::VecDeque::new();
            // Output batches.
            let mut batch_r1: Vec<FastqRecord> = Vec::with_capacity(BAM_READER_BATCH_SIZE);
            let mut batch_r2: Vec<FastqRecord> = Vec::with_capacity(BAM_READER_BATCH_SIZE);
            let mut record_idx: usize = 0;

            const FREAD1: u16 = 0x40;
            const FREAD2: u16 = 0x80;

            // Helper to dual-emit a fatal error message.
            let send_pair_err = |tx_r1: &mpsc::SyncSender<Result<Vec<FastqRecord>>>,
                                 tx_r2: &mpsc::SyncSender<Result<Vec<FastqRecord>>>,
                                 msg: String| {
                let _ = tx_r1.send(Err(anyhow::anyhow!("{}", msg)));
                let _ = tx_r2.send(Err(anyhow::anyhow!("{}", msg)));
            };

            loop {
                let mut bam_rec = bam::Record::default();
                match inner.read_record(&mut bam_rec) {
                    Ok(0) => {
                        // EOF — flush remaining batches, then check de-interleaver.
                        if !batch_r1.is_empty() {
                            let _ = tx_r1.send(Ok(batch_r1));
                            let _ = tx_r2.send(Ok(batch_r2));
                        }
                        if !r1_queue.is_empty() || !r2_queue.is_empty() {
                            send_pair_err(
                                &tx_r1,
                                &tx_r2,
                                format!(
                                    "Paired uBAM '{}' ended with {} unmatched R1 and {} unmatched R2 records (orphans). Consider --retain_unpaired or re-interleave with `samtools collate`.",
                                    path_for_thread.display(),
                                    r1_queue.len(),
                                    r2_queue.len()
                                ),
                            );
                            return;
                        }
                        // EOF sentinels.
                        let _ = tx_r1.send(Ok(Vec::new()));
                        let _ = tx_r2.send(Ok(Vec::new()));
                        return;
                    }
                    Ok(_) => {
                        record_idx += 1;
                        let flags = bam_rec.flags().bits();
                        let is_r1 = flags & FREAD1 != 0;
                        let is_r2 = flags & FREAD2 != 0;
                        if is_r1 == is_r2 {
                            // Both or neither — invalid for paired uBAM.
                            send_pair_err(
                                &tx_r1,
                                &tx_r2,
                                format!(
                                    "Paired uBAM '{}': record {} has invalid flag bits (FREAD1={}, FREAD2={}, raw=0x{:04X}). Each record must have exactly one of READ1 (0x40) / READ2 (0x80) set.",
                                    path_for_thread.display(),
                                    record_idx,
                                    is_r1,
                                    is_r2,
                                    flags
                                ),
                            );
                            return;
                        }
                        // Convert (also runs per-record aligned-BAM check).
                        let fq = match bam_record_to_fastq(&bam_rec, &tags_for_thread) {
                            Ok(r) => r,
                            Err(e) => {
                                send_pair_err(
                                    &tx_r1,
                                    &tx_r2,
                                    format!(
                                        "{}: BAM record {}: {:#}",
                                        path_for_thread.display(),
                                        record_idx,
                                        e
                                    ),
                                );
                                return;
                            }
                        };
                        // Route to the correct side queue.
                        if is_r1 {
                            r1_queue.push_back(fq);
                        } else {
                            r2_queue.push_back(fq);
                        }
                        // Emit a matched pair if both queues are non-empty.
                        if !r1_queue.is_empty() && !r2_queue.is_empty() {
                            batch_r1.push(r1_queue.pop_front().unwrap());
                            batch_r2.push(r2_queue.pop_front().unwrap());
                            if batch_r1.len() >= BAM_READER_BATCH_SIZE {
                                let send_r1 = std::mem::replace(
                                    &mut batch_r1,
                                    Vec::with_capacity(BAM_READER_BATCH_SIZE),
                                );
                                let send_r2 = std::mem::replace(
                                    &mut batch_r2,
                                    Vec::with_capacity(BAM_READER_BATCH_SIZE),
                                );
                                if tx_r1.send(Ok(send_r1)).is_err()
                                    || tx_r2.send(Ok(send_r2)).is_err()
                                {
                                    return; // receiver(s) dropped
                                }
                            }
                        }
                        // Bound the buffer — if either side has exceeded MAX_SLACK
                        // pending records without a mate, the input is grouped.
                        if r1_queue.len() > MAX_SLACK || r2_queue.len() > MAX_SLACK {
                            send_pair_err(&tx_r1, &tx_r2, GROUPED_INPUT_ERR.to_string());
                            return;
                        }
                    }
                    Err(io_err) => {
                        send_pair_err(
                            &tx_r1,
                            &tx_r2,
                            format!(
                                "{}: failed to read BAM record {}: {}",
                                path_for_thread.display(),
                                record_idx + 1,
                                io_err
                            ),
                        );
                        return;
                    }
                }
            }
        });

        let r1 = Self {
            source: BamReaderSource::Threaded {
                rx: rx_r1,
                buffer: Vec::new(),
                buf_pos: 0,
                _handle: spawn_noop_joinhandle(),
            },
            preserved_tags: preserve_tags.to_vec(),
        };
        let r2 = Self {
            source: BamReaderSource::Threaded {
                rx: rx_r2,
                buffer: Vec::new(),
                buf_pos: 0,
                _handle: spawn_noop_joinhandle(),
            },
            preserved_tags: preserve_tags.to_vec(),
        };
        // `_producer_handle` drops here → thread detaches and continues.
        Ok((r1, r2))
    }

    fn open_threaded_inner<P: AsRef<Path>>(path: P, preserved_tags: Vec<String>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            bail!("Input file not found: {}", path.display());
        }
        let (tx, rx) = mpsc::sync_channel(BAM_READER_CHANNEL_BATCHES);
        let tags_for_thread = preserved_tags.clone();
        let path_for_thread = path.clone();

        let handle = std::thread::spawn(move || {
            let mut inner = match open_inner(&path_for_thread) {
                Ok((r, _header)) => r,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            };
            let mut batch = Vec::with_capacity(BAM_READER_BATCH_SIZE);
            let mut record_idx: usize = 0;
            loop {
                let mut record = bam::Record::default();
                match inner.read_record(&mut record) {
                    Ok(0) => {
                        if !batch.is_empty() {
                            let _ = tx.send(Ok(batch));
                        }
                        let _ = tx.send(Ok(Vec::new()));
                        break;
                    }
                    Ok(_) => {
                        record_idx += 1;
                        match bam_record_to_fastq(&record, &tags_for_thread) {
                            Ok(fq) => {
                                batch.push(fq);
                                if batch.len() >= BAM_READER_BATCH_SIZE {
                                    let full = std::mem::replace(
                                        &mut batch,
                                        Vec::with_capacity(BAM_READER_BATCH_SIZE),
                                    );
                                    if tx.send(Ok(full)).is_err() {
                                        break;
                                    }
                                }
                            }
                            Err(e) => {
                                let _ = tx.send(Err(e.context(format!(
                                    "{}: BAM record {}",
                                    path_for_thread.display(),
                                    record_idx
                                ))));
                                break;
                            }
                        }
                    }
                    Err(e) => {
                        let _ = tx.send(Err(anyhow::Error::from(e).context(format!(
                            "{}: failed to read BAM record {}",
                            path_for_thread.display(),
                            record_idx + 1
                        ))));
                        break;
                    }
                }
            }
        });

        Ok(Self {
            source: BamReaderSource::Threaded {
                rx,
                buffer: Vec::new(),
                buf_pos: 0,
                _handle: handle,
            },
            preserved_tags,
        })
    }

    /// Set the BAM tags to preserve in FASTQ headers. Default: no tags.
    /// Each tag must be a valid 2-char SAM tag name (validated at the CLI layer).
    ///
    /// Note: when transitioning to threaded mode, call `with_preserved_tags`
    /// BEFORE `open_threaded` — the threaded reader snapshots tags into the
    /// background thread at construction. For a `Direct` reader, this method
    /// works at any time.
    pub fn with_preserved_tags(mut self, tags: &[String]) -> Self {
        self.preserved_tags = tags.to_vec();
        self
    }

    /// Read the next record, converting to `FastqRecord`. Returns `None` at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        match &mut self.source {
            BamReaderSource::Direct {
                inner,
                path,
                record_idx,
            } => {
                let mut record = bam::Record::default();
                match inner.read_record(&mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => {
                        *record_idx += 1;
                        let rec = bam_record_to_fastq(&record, &self.preserved_tags).with_context(
                            || format!("{}: BAM record {}", path.display(), record_idx),
                        )?;
                        Ok(Some(rec))
                    }
                    Err(e) => Err(anyhow::Error::from(e).context(format!(
                        "{}: failed to read BAM record {}",
                        path.display(),
                        *record_idx + 1
                    ))),
                }
            }
            BamReaderSource::Threaded {
                rx,
                buffer,
                buf_pos,
                ..
            } => {
                if *buf_pos < buffer.len() {
                    let idx = *buf_pos;
                    *buf_pos += 1;
                    let record = std::mem::replace(
                        &mut buffer[idx],
                        FastqRecord {
                            id: String::new(),
                            seq: String::new(),
                            qual: String::new(),
                        },
                    );
                    return Ok(Some(record));
                }
                match rx.recv() {
                    Ok(Ok(batch)) => {
                        if batch.is_empty() {
                            Ok(None)
                        } else {
                            *buffer = batch;
                            *buf_pos = 1;
                            let record = std::mem::replace(
                                &mut buffer[0],
                                FastqRecord {
                                    id: String::new(),
                                    seq: String::new(),
                                    qual: String::new(),
                                },
                            );
                            Ok(Some(record))
                        }
                    }
                    Ok(Err(e)) => Err(e),
                    Err(_) => Ok(None),
                }
            }
        }
    }
}

// ─── BamWriter (uBAM output, Phase 1 v2.1) ─────────────────────────────────
//
// Writes one BAM record per `FastqRecord`. Aux tags are recovered from the
// textual tail of `FastqRecord.id` per PLAN §3.6 — the round-trip
// counterpart to `BamReader::bam_record_to_fastq`'s read-side encoding.
//
// This is a concrete struct, NOT a `RecordSink` trait impl: per PLAN v2 the
// FASTQ-output path stays on parallel.rs unchanged; uBAM output is a
// separate serial dispatch branch wired up in Step 4 (`main.rs`).

/// uBAM writer. Wraps `noodles::bam::io::Writer` and carries the SAM header
/// needed by the `write_alignment_record` trait method.
///
/// The writer does NOT track a tag-filter list. `BamReader::next_record`
/// already filtered tags into the textual `FastqRecord.id` tail in
/// user-specified order, so the tail IS the source of truth — the writer
/// just round-trips whatever it finds.
pub struct BamWriter {
    /// Declared before `pending` so the BGZF stream is closed before the
    /// temporary is removed on the abandon path.
    inner: bam::io::Writer<bgzf::Writer<BufWriter<File>>>,
    /// `None` once [`BamWriter::finish`] has taken it. While it is `Some`, a
    /// drop removes the temporary and no output file is created.
    pending: Option<crate::io::PendingOutput>,
    header: Header,
    /// ASCII offset of the incoming `FastqRecord.qual` — 33, or 64 under
    /// `--phred64`. Subtracted per byte on write, because BAM `QUAL` stores
    /// true Phred scores. See [`PHRED_OFFSET`] for why this is a separate,
    /// variable value rather than that constant.
    input_phred_offset: u8,
}

impl BamWriter {
    /// Open `path` for uBAM writing.
    ///
    /// - `source_header`: input uBAM's SAM header (propagated as-is plus a
    ///   trim_galore `@PG` line). Pass `None` for FASTQ input — a minimal
    ///   `@HD VN:1.6` + trim_galore `@PG` header is synthesised.
    /// - `_preserve_tags`: kept for API symmetry with `BamReader::open_*_with_tags`;
    ///   not consulted (the textual tail in `FastqRecord.id` is authoritative).
    /// - `command_line`: trim_galore invocation embedded in `@PG CL:`.
    /// - `input_phred_offset`: ASCII offset of the *incoming* `FastqRecord.qual`
    ///   (33, or 64 under `--phred64`); pass `cli.phred_offset()`. BAM `QUAL`
    ///   holds true Phred scores (0–93) per the SAM spec, so this is subtracted
    ///   per byte on write. Required rather than defaulted on purpose: issue
    ///   #358 was caused by a call site that never considered the offset, so
    ///   the compiler should force every site to state one.
    ///
    /// **PLAN §4 deviation:** the plan listed 3 params; we add a 4th
    /// (`command_line`) to keep header construction self-contained inside
    /// the writer instead of forcing the caller to clone the header just
    /// to mutate it. Documented in the v2.1 deviation log. (A 5th,
    /// `input_phred_offset`, was added later for issue #358 — see above.)
    pub fn create<P: AsRef<Path>>(
        path: P,
        source_header: Option<&Header>,
        _preserve_tags: &[String],
        command_line: &str,
        input_phred_offset: u8,
    ) -> Result<Self> {
        let path = path.as_ref();
        let header = build_output_header(source_header, command_line)?;
        let (pending, file) = crate::io::PendingOutput::create(path)?;
        let mut inner = bam::io::Writer::new(BufWriter::new(file));
        inner
            .write_header(&header)
            .with_context(|| format!("Failed to write BAM header to {}", path.display()))?;
        Ok(Self {
            inner,
            pending: Some(pending),
            header,
            input_phred_offset,
        })
    }

    /// Write one `FastqRecord` as a BAM record.
    ///
    /// `paired_side` selects flag bits per PLAN §3.3 step 2:
    /// - `None`     → single-end, flag `0x04` (FUNMAP).
    /// - `Some(1)`  → paired R1, flag `0x4D = 77` (PAIRED|UNMAPPED|MATE_UNMAPPED|READ1).
    /// - `Some(2)`  → paired R2, flag `0x8D = 141` (PAIRED|UNMAPPED|MATE_UNMAPPED|READ2).
    ///
    /// Aux fields are parsed from the textual tail of `record.id` per
    /// PLAN §3.6 — only A / Z / i / f scalar types are accepted; B / H
    /// arrays or malformed shape errors.
    pub fn write_record(&mut self, record: &FastqRecord, paired_side: Option<u8>) -> Result<()> {
        let (name, data) = parse_name_and_data(&record.id)
            .with_context(|| format!("malformed FastqRecord.id: {:?}", record.id))?;
        let flags = match paired_side {
            None => Flags::UNMAPPED,
            Some(1) => {
                Flags::SEGMENTED | Flags::UNMAPPED | Flags::MATE_UNMAPPED | Flags::FIRST_SEGMENT
            }
            Some(2) => {
                Flags::SEGMENTED | Flags::UNMAPPED | Flags::MATE_UNMAPPED | Flags::LAST_SEGMENT
            }
            Some(other) => bail!(
                "BamWriter::write_record: paired_side must be 1 or 2, got {}",
                other
            ),
        };
        // Sequence — validate per PLAN §3.3 step 3 (ACGTN-only, reject `=`,
        // coerce IUPAC degenerate to N). Mirrors the read-side
        // `bam_record_to_fastq` logic for symmetry. Code-review B-M1 fix.
        let normalized_seq = validate_and_normalize_seq_for_write(record.seq.as_bytes())
            .with_context(|| format!("while writing BAM record {:?}", record.id))?;
        // BAM stores raw Phred bytes; `FastqRecord.qual` is ASCII at the
        // input's declared offset — 33, or 64 under `--phred64`. Hardcoding 33
        // here was issue #358: Phred+64 input was under-subtracted by 31 and
        // stored as inflated-but-SAM-legal quality, so nothing errored.
        //
        // `saturating_sub` is load-bearing, not merely defensive: nothing in
        // `fastq.rs` validates the quality range (`FastqReader::sanity_check`
        // only checks the `@` prefix and rejects colorspace; the record path
        // only detects truncation), so a malformed sub-offset byte can reach
        // here and must floor to 0 rather than wrap.
        let raw_qual: Vec<u8> = record
            .qual
            .bytes()
            .map(|b| b.saturating_sub(self.input_phred_offset))
            .collect();
        debug_assert_eq!(
            raw_qual.len(),
            normalized_seq.len(),
            "FastqRecord seq/qual length mismatch reached BamWriter"
        );
        let rec = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_sequence(Sequence::from(normalized_seq))
            .set_quality_scores(QualityScores::from(raw_qual))
            .set_data(data)
            .build();
        // `record.id`, not `name` — that moved into `set_name` above. Shown under the
        // 60-char cap the tag notices use, with the byte count when it bites, because
        // a FASTQ id has no length limit and one rejection here is over-long. The
        // QNAME rule is appended only for `InvalidInput`, so a full disk or a closed
        // pipe is not relabelled as a read-name problem.
        self.inner
            .write_alignment_record(&self.header, &rec)
            .map_err(|e| {
                let shown = if record.id.chars().count() > 60 {
                    format!(
                        "{}… ({} bytes)",
                        record.id.chars().take(60).collect::<String>(),
                        record.id.len()
                    )
                } else {
                    record.id.clone()
                };
                let mut msg = format!("while writing BAM record \"{}\"", shown.escape_debug());
                if e.kind() == std::io::ErrorKind::InvalidInput {
                    msg.push_str(
                        " — BAM rejected it. A SAM read name must be 1-254 bytes of \
                         printable ASCII other than `@` and space. Rename the reads \
                         before trimming, or write FASTQ output, which accepts them.",
                    );
                }
                anyhow::Error::new(e).context(msg)
            })?;
        Ok(())
    }

    /// Flush and finalise the writer — writes the BGZF EOF marker, then
    /// publishes the output under its final name. Consumes `self` so the
    /// caller can't forget.
    ///
    /// `try_finish` before the rename, and `pending` out before the drop.
    pub fn finish(mut self) -> Result<()> {
        self.inner
            .try_finish()
            .context("failed to finalise uBAM output (BGZF EOF marker)")?;
        let pending = self.pending.take();
        drop(self);
        match pending {
            Some(p) => p.commit(),
            None => Ok(()),
        }
    }
}

/// Build the output header — either propagate `source` + append our `@PG`
/// via noodles' [`Programs::add`] (which auto-handles `PP:` chain linkage
/// AND ID disambiguation when the prefix collides), or synthesise a
/// minimal `@HD VN:1.6` + our `@PG` from scratch.
///
/// Re-processing a trim_galore-output BAM produces a second `@PG` with
/// ID `trim_galore-trim_galore` (PLAN §3.5 "provenance preserved by
/// adding to history, not preserving silence"). A 3-deep nesting (input
/// already has both `trim_galore` and `trim_galore-trim_galore`) returns
/// an `io::Result` error from `Programs::add`; surface it.
fn build_output_header(source: Option<&Header>, command_line: &str) -> Result<Header> {
    let version = env!("CARGO_PKG_VERSION");
    // Build @PG WITHOUT setting PP — `Programs::add` computes PP from
    // the chain's leaf program(s) automatically.
    let pg = build_pg_map(version, command_line);
    match source {
        Some(src) => {
            let mut h = src.clone();
            h.programs_mut().add("trim_galore", pg).with_context(|| {
                "failed to append trim_galore @PG line to the source header \
                 (likely a 3-deep re-trim chain collision); rename the source \
                 BAM or strip its @PG block first"
            })?;
            Ok(h)
        }
        None => Ok(Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_program("trim_galore", pg)
            .build()),
    }
}

/// Build a `Map<Program>` with PN / VN / CL. PP is NOT set here —
/// `Programs::add` computes it from the existing chain's leaves.
fn build_pg_map(version: &str, command_line: &str) -> Map<Program> {
    Map::<Program>::builder()
        .insert(program_tag::NAME, "trim_galore")
        .insert(program_tag::VERSION, version)
        .insert(program_tag::COMMAND_LINE, command_line)
        // `Program` is a unit struct with no required standard fields; .build()
        // is infallible in practice.
        .build()
        .expect("Map::<Program>::build is infallible for a unit Program")
}

/// Parse `FastqRecord.id` into `(name_bytes, aux_data)` per PLAN §3.3
/// step 1 + §3.6.
///
/// Shape: `@NAME[ DESC]?[\tTAG:TYPE:VALUE]*`. The leading `@` is stripped;
/// the first whitespace (TAB or SPACE) separates the record name from
/// either an aux tag tail (TAB) or a descriptive FASTQ-only annotation
/// (SPACE — discarded). The BAM spec rejects whitespace in QNAME, so
/// splitting on either is required for FASTQ-input → uBAM-output to
/// produce a valid BAM. Remaining tab-separated fields parse as
/// `TAG:TYPE:VALUE`; only A / Z / i / f scalar types are accepted.
fn parse_name_and_data(id: &str) -> Result<(Vec<u8>, Data)> {
    let stripped = id.strip_prefix('@').unwrap_or(id);
    let (name_str, tail) = if let Some(pos) = stripped.find(['\t', ' ']) {
        let sep = stripped.as_bytes()[pos];
        let name_str = &stripped[..pos];
        let rest = &stripped[pos + 1..];
        if sep == b'\t' {
            // Tab: rest is the aux tag tail (from BamReader::next_record).
            (name_str, Some(rest))
        } else {
            // Space: rest is descriptive FASTQ annotation — discard. Gated on
            // non-empty: a trailing space drops nothing (#406).
            if !rest.is_empty() {
                emit_description_dropped_once(rest);
            }
            (name_str, None)
        }
    } else {
        (stripped, None)
    };
    let name = name_str.as_bytes().to_vec();
    if name.is_empty() {
        bail!("FastqRecord.id has empty record name");
    }
    let mut data = Data::default();
    if let Some(tail) = tail {
        for field in tail.split('\t') {
            if field.is_empty() {
                continue; // tolerate stray trailing tabs
            }
            let mut parts = field.splitn(3, ':');
            let tag_name = parts
                .next()
                .ok_or_else(|| anyhow::anyhow!("malformed tag field: {:?}", field))?;
            let type_code = parts
                .next()
                .ok_or_else(|| anyhow::anyhow!("missing type code in tag field: {:?}", field))?;
            let value_str = parts
                .next()
                .ok_or_else(|| anyhow::anyhow!("missing value in tag field: {:?}", field))?;
            if tag_name.len() != 2 {
                bail!("tag name must be 2 chars, got {:?}", tag_name);
            }
            if type_code.len() != 1 {
                bail!("tag type must be 1 char, got {:?}", type_code);
            }
            let tag_bytes = tag_name.as_bytes();
            let tag = Tag::new(tag_bytes[0], tag_bytes[1]);
            let value = parse_tag_value(type_code.as_bytes()[0], value_str)
                .with_context(|| format!("while parsing tag {}:{}", tag_name, type_code))?;
            data.insert(tag, value);
        }
    }
    Ok((name, data))
}

/// Parse a single typed-tag value per the SAM aux encoding. Only A/Z/i/f
/// scalar types round-trip through TrimGalore's FASTQ intermediate; #317's
/// reader already rejects B/H on read, so they shouldn't appear in normal
/// flow — we still defend against hand-crafted tails.
fn parse_tag_value(type_code: u8, raw: &str) -> Result<BufValue> {
    match type_code {
        b'A' => {
            let mut bytes = raw.bytes();
            let c = bytes
                .next()
                .ok_or_else(|| anyhow::anyhow!("empty A: tag value"))?;
            if bytes.next().is_some() {
                bail!("A: tag value must be exactly one byte, got {:?}", raw);
            }
            Ok(BufValue::Character(c))
        }
        b'Z' => Ok(BufValue::from(raw)),
        b'i' => raw
            .parse::<i32>()
            .map(BufValue::Int32)
            .map_err(|e| anyhow::anyhow!("invalid i: tag value {:?}: {}", raw, e)),
        b'f' => raw
            .parse::<f32>()
            .map(BufValue::Float)
            .map_err(|e| anyhow::anyhow!("invalid f: tag value {:?}: {}", raw, e)),
        b'B' | b'H' => bail!(
            "BAM array (B:) and hex (H:) aux tag types do not round-trip through \
             --output-format ubam in v1; only A / Z / i / f scalars are supported"
        ),
        other => bail!(
            "unknown aux tag type code 0x{:02X} (only A / Z / i / f accepted)",
            other
        ),
    }
}

/// Validate + normalise a FASTQ-side sequence for BAM writing, per PLAN
/// §3.3 step 3. Coerces lowercase to uppercase, IUPAC degenerate bases
/// (R/Y/M/K/S/W/B/D/H/V) to N (write-side warning emitted once), and rejects `=`
/// plus any non-IUPAC byte. Symmetric with the read-side
/// `bam_record_to_fastq` validation.
fn validate_and_normalize_seq_for_write(seq: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(seq.len());
    let mut iupac_seen = false;
    for &b in seq {
        let upper = b.to_ascii_uppercase();
        match upper {
            b'A' | b'C' | b'G' | b'T' | b'N' => out.push(upper),
            b'=' => bail!(
                "FastqRecord sequence contains '=' (match-reference sentinel); \
                 invalid for uBAM output"
            ),
            b'R' | b'Y' | b'M' | b'K' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V' => {
                out.push(b'N');
                iupac_seen = true;
            }
            other => bail!(
                "FastqRecord sequence contains unrecognised base byte 0x{:02X} \
                 (expected ACGTN or IUPAC code)",
                other
            ),
        }
    }
    if iupac_seen {
        emit_iupac_warning_write_once();
    }
    Ok(out)
}

/// Construct a JoinHandle that wraps a no-op thread, used as a sentinel for
/// the paired-interleaved readers (only one of the two readers needs the
/// real handle; the other gets this). Cheap — the spawned thread exits
/// immediately. Workaround for the struct field needing a concrete value.
fn spawn_noop_joinhandle() -> JoinHandle<()> {
    std::thread::spawn(|| {})
}

/// Open a BAM file, read its header, return both the reader and the parsed
/// header. The header is needed by `BamWriter` to propagate `@HD`/`@PG`/`@CO`
/// lines through uBAM-in → uBAM-out (PLAN §3.5).
fn open_inner(path: &Path) -> Result<(bam::io::Reader<bgzf::Reader<BufReader<File>>>, Header)> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open uBAM input: {}", path.display()))?;
    let mut inner = bam::io::Reader::new(BufReader::new(file));
    let header = inner
        .read_header()
        .with_context(|| format!("Failed to read BAM header from {}", path.display()))?;
    Ok((inner, header))
}

/// Read just the SAM header from `path` without keeping the reader open.
///
/// Used by `main.rs::run_ubam_output` to extract the source header for
/// `BamWriter::create` before the streaming reader is constructed. The
/// file is opened and closed (one extra `open` syscall beyond the
/// streaming reader), at microseconds of overhead — far below the
/// per-record cost.
pub fn peek_header(path: &Path) -> Result<Header> {
    let (_reader, header) = open_inner(path)?;
    Ok(header)
}

/// Reject a BAM read name holding a byte that breaks the FASTQ id it becomes.
///
/// The chosen set is `u8::is_ascii_whitespace` — what `parse_name_and_data`,
/// `fastq::read_id_prefix` and the id's one-line framing split on. VT sits deliberately
/// outside it; `append_to_id`'s `trim_end` does strip a trailing one.
fn reject_whitespace_in_qname(name: &[u8]) -> Result<()> {
    if name.iter().any(u8::is_ascii_whitespace) {
        let shown = String::from_utf8_lossy(name);
        bail!(
            "BAM read name contains whitespace: \"{}\". Trim Galore reads a FASTQ header as \
             ending at the first whitespace, so any --preserve-tags aux tags after it are lost \
             on uBAM output, and a newline would break the FASTQ record structure. The SAM \
             specification forbids whitespace in QNAME. Rename the reads before trimming — the \
             uBAM section of the docs gives a samtools recipe for space-only names; a tab or \
             newline in a name cannot be expressed in SAM text at all.",
            shown.escape_debug()
        );
    }
    Ok(())
}

/// Reject bytes in an appended `TYPE:VALUE` that would break the FASTQ id holding it.
///
/// The QNAME boundary set minus space, which `Z`'s `[ !-~]*` grammar makes legal: a tab
/// would forge a tag field and a newline would break the single-line id.
fn reject_framing_in_tag_value(text: &str) -> Result<()> {
    if text.bytes().any(|b| b != b' ' && b.is_ascii_whitespace()) {
        // `Z` values are unbounded, so borrow the 60-char cap the description notice uses.
        let shown = if text.chars().count() > 60 {
            format!("{}…", text.chars().take(60).collect::<String>())
        } else {
            text.to_string()
        };
        bail!(
            "aux tag value cannot pass through the FASTQ header: \"{}\". A tab would forge a \
             tag field and a newline would break the single-line read ID; spaces are legal and \
             unaffected. Fix the tag at source.",
            shown.escape_debug()
        );
    }
    Ok(())
}

/// Convert one BAM record to a `FastqRecord`. Per-record validation:
/// - Must be unmapped (`is_unmapped()` true; per-record check resolves
///   PLAN-REVIEW B-Crit-4 first-record-only contradiction).
/// - Must NOT be reverse-complemented / secondary / supplementary.
/// - Sequence must be non-empty; IUPAC codes coerced to N on read, `=` rejected.
/// - Qual length must match seq length (or be all-`0xFF` / empty — missing).
fn bam_record_to_fastq(rec: &bam::Record, tags: &[String]) -> Result<FastqRecord> {
    let flags = rec.flags();
    if !flags.is_unmapped() {
        bail!(
            "aligned BAM record encountered (FUNMAP=0); Trim Galore only accepts \
             unaligned BAM (uBAM). Use 'samtools view -h -f 4' to filter or \
             'samtools fastq' to convert"
        );
    }
    if flags.is_reverse_complemented() {
        bail!("BAM record has reverse-complement flag (0x10) set; not valid in unaligned BAM");
    }
    if flags.is_secondary() {
        bail!("BAM record is secondary (flag 0x100); not valid in unaligned BAM");
    }
    if flags.is_supplementary() {
        bail!("BAM record is supplementary (flag 0x800); not valid in unaligned BAM");
    }

    // ID — `@{record_name}`. Empty name is an error.
    let name = rec
        .name()
        .ok_or_else(|| anyhow::anyhow!("BAM record has empty read name"))?;
    if name.is_empty() {
        bail!("BAM record has empty read name");
    }
    reject_whitespace_in_qname(name)?;
    let id_body = std::str::from_utf8(name).context("BAM record name is not valid UTF-8")?;
    let mut id = format!("@{}", id_body);

    // Tag preservation — append `\tTAG:TYPE:VALUE` in user-specified order.
    // Missing tags from a specific record are silently skipped (per-record).
    if !tags.is_empty() {
        let data = rec.data();
        for tag_name in tags {
            let bytes = tag_name.as_bytes();
            debug_assert_eq!(
                bytes.len(),
                2,
                "CLI parser must enforce 2-char SAM tag names"
            );
            let tag = Tag::new(bytes[0], bytes[1]);
            if let Some(value_result) = data.get(&tag) {
                let value = value_result.with_context(|| {
                    format!("failed to read tag '{}' from BAM record", tag_name)
                })?;
                id.push('\t');
                id.push_str(tag_name);
                id.push(':');
                append_tag_type_and_value(&mut id, &value)
                    .with_context(|| format!("aux tag '{}'", tag_name))?;
            }
        }
    }

    // Sequence — decode each byte (already ASCII per noodles' `Sequence::iter`).
    let seq_view = rec.sequence();
    if seq_view.is_empty() {
        bail!("BAM record has empty sequence");
    }
    let mut seq = String::with_capacity(seq_view.len());
    let mut iupac_seen = false;
    for base in seq_view.iter() {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' => seq.push(base as char),
            b'=' => bail!(
                "BAM record contains '=' (match-reference sentinel); invalid in unaligned BAM"
            ),
            b'R' | b'Y' | b'M' | b'K' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V' => {
                seq.push('N');
                iupac_seen = true;
            }
            other => bail!(
                "BAM record contains unrecognised base byte: 0x{:02X}",
                other
            ),
        }
    }
    if iupac_seen {
        emit_iupac_warning_read_once();
    }

    // Qual — raw Phred (0–93) + 33 → Sanger ASCII. Missing qual (BAM convention:
    // all-`0xFF`, or empty) → `'!'` × seq_len, matching `samtools fastq`.
    let qual_view = rec.quality_scores();
    let qual_raw: &[u8] = qual_view.as_ref();
    let qual = if qual_raw.is_empty() || qual_raw.iter().all(|&b| b == QUAL_MISSING_SENTINEL) {
        std::iter::repeat_n(QUAL_MISSING_REPLACEMENT as char, seq.len()).collect()
    } else {
        if qual_raw.len() != seq.len() {
            bail!(
                "BAM record qual length ({}) does not match seq length ({})",
                qual_raw.len(),
                seq.len()
            );
        }
        // The sentinel check above only fires when *every* byte is `0xFF`, so a
        // record with a MIXED qual (some real scores, some `0xFF`) reaches here.
        // Two hazards on that input, both handled per byte:
        //
        //  1. Plain `b + PHRED_OFFSET` overflows `u8` for any raw byte > 222 —
        //     panic under debug assertions, silent wrap in release.
        //  2. `saturating_add` avoids the overflow but is *worse*: `255u8 as
        //     char` is U+00FF, which occupies TWO bytes in UTF-8. That breaks
        //     the codebase-wide invariant that `qual` has one byte per base —
        //     `FastqRecord::truncate` / `clip_5prime` / `clip_3prime` index this
        //     `String` with byte offsets derived from `seq.len()`, so a
        //     multi-byte char makes them panic on a char boundary (in release
        //     too) or silently keep the wrong number of scores.
        //
        // So: map the sentinel to `QUAL_MISSING_REPLACEMENT` exactly as the
        // all-`0xFF` branch does — the mixed case now agrees with it rather
        // than diverging — and clamp anything else out of spec to the SAM
        // maximum (93 → `'~'`). Every output byte is single-byte ASCII, so the
        // one-byte-per-base invariant holds unconditionally.
        qual_raw
            .iter()
            .map(|&b| {
                if b == QUAL_MISSING_SENTINEL {
                    QUAL_MISSING_REPLACEMENT as char
                } else {
                    (b.min(MAX_PHRED_SCORE) + PHRED_OFFSET) as char
                }
            })
            .collect()
    };

    Ok(FastqRecord { id, seq, qual })
}

/// Read-side IUPAC coercion (uBAM in). Output-neutral on purpose: this fires for
/// uBAM→FASTQ and uBAM→uBAM alike, so it names where the coercion happens rather
/// than a destination that may not exist (#406).
fn emit_iupac_warning_read_once() {
    static SEEN: OnceLock<()> = OnceLock::new();
    SEEN.get_or_init(|| {
        eprintln!(
            "WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V); \
             coerced to N on read. This is legitimate in PacBio HiFi, ONT Dorado, \
             and 10x cellranger uBAMs."
        );
    });
}

/// Write-side IUPAC coercion (FASTQ in, uBAM out). Guarded separately from the
/// read side so a mixed FASTQ+uBAM run reports both, rather than whichever
/// direction happened to come first (#406).
fn emit_iupac_warning_write_once() {
    static SEEN: OnceLock<()> = OnceLock::new();
    SEEN.get_or_init(|| {
        eprintln!(
            "WARNING: input FASTQ contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V); \
             coerced to N in the uBAM output. This is legitimate in PacBio HiFi, \
             ONT Dorado, and 10x cellranger data."
        );
    });
}

/// One-time disclosure that a FASTQ header description was not carried into uBAM
/// output. Echoes the text actually dropped, because what is lost varies: an
/// instrument identifier, or an Illumina `1:N:0:INDEX` field.
///
/// "FASTQ" is unconditional: BAM input cannot reach this, because a whitespace QNAME is
/// refused and the tab introducing a tag tail always precedes any space in a tag value.
fn emit_description_dropped_once(dropped: &str) {
    static SEEN: OnceLock<()> = OnceLock::new();
    SEEN.get_or_init(|| {
        let shown: String = if dropped.chars().count() > 60 {
            format!("{}…", dropped.chars().take(60).collect::<String>())
        } else {
            dropped.to_string()
        };
        eprintln!(
            "NOTE: BAM read names cannot contain whitespace, so FASTQ header text after \
             the first space is not carried into uBAM output. First dropped: \"{}\"",
            shown
        );
    });
}

/// Format `{TYPE}:{VALUE}` for one BAM aux field, matching `samtools fastq -T`.
///
/// Every arm's output is scanned before returning, so no type code can route text into the
/// FASTQ id unchecked.
fn append_tag_type_and_value(out: &mut String, value: &Value<'_>) -> Result<()> {
    use std::fmt::Write;
    let start = out.len();
    match value {
        Value::Character(c) => write!(out, "A:{}", *c as char).unwrap(),
        Value::Int8(v) => write!(out, "i:{}", v).unwrap(),
        Value::UInt8(v) => write!(out, "i:{}", v).unwrap(),
        Value::Int16(v) => write!(out, "i:{}", v).unwrap(),
        Value::UInt16(v) => write!(out, "i:{}", v).unwrap(),
        Value::Int32(v) => write!(out, "i:{}", v).unwrap(),
        Value::UInt32(v) => write!(out, "i:{}", v).unwrap(),
        Value::Float(v) => write!(out, "f:{}", v).unwrap(),
        Value::String(bs) => {
            let s = std::str::from_utf8(bs).context("BAM string tag is not valid UTF-8")?;
            write!(out, "Z:{}", s).unwrap();
        }
        Value::Hex(_) => {
            // Hex tags (`H:…`) don't round-trip through the BamWriter's
            // §3.6 parser (which only accepts A/Z/i/f scalars). Reject
            // symmetrically with B: arrays — otherwise an H tag would be
            // emitted into the FASTQ tail on read and then hard-error
            // mid-stream on write. Code-review A-I2 fix.
            bail!(
                "BAM hex (H:) tag values are not supported in --preserve-tags v1 \
                 (does not round-trip through --output-format ubam); \
                 omit the hex tag from --preserve-tags or open an issue if needed"
            );
        }
        Value::Array(_) => {
            // Array tags (`B:…`) are uncommon for the typical use case
            // (CB / UB / RX are strings). Defer to a follow-up.
            bail!(
                "BAM array tag values are not supported in --preserve-tags v1; \
                 omit the array tag from --preserve-tags or open an issue if needed"
            );
        }
    }
    reject_framing_in_tag_value(&out[start..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_committed_se_fixture() -> Result<()> {
        let _r = BamReader::open("test_files/ubam_test.bam")?;
        Ok(())
    }

    #[test]
    fn open_committed_pe_fixture() -> Result<()> {
        let _r = BamReader::open("test_files/ubam_paired_test.bam")?;
        Ok(())
    }

    #[test]
    fn reads_all_se_records() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_test.bam")?;
        let mut count = 0;
        while r.next_record()?.is_some() {
            count += 1;
        }
        assert_eq!(count, 10);
        Ok(())
    }

    #[test]
    fn reads_all_pe_records() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_paired_test.bam")?;
        let mut count = 0;
        while r.next_record()?.is_some() {
            count += 1;
        }
        // 10 R1 + 10 R2 records (interleaved); both surface to next_record here
        // (the de-interleave logic lives in `open_paired_interleaved`, Task 15).
        assert_eq!(count, 20);
        Ok(())
    }

    #[test]
    fn first_record_id_starts_with_at_sign() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_test.bam")?;
        let rec = r.next_record()?.expect("fixture has ≥1 record");
        assert!(
            rec.id.starts_with('@'),
            "FASTQ ID must start with '@', got: {:?}",
            rec.id
        );
        // First record name from the BS-seq fixture
        assert_eq!(rec.id, "@SRR24827378.1");
        Ok(())
    }

    #[test]
    fn first_record_seq_and_qual_lengths_match() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_test.bam")?;
        let rec = r.next_record()?.expect("fixture has ≥1 record");
        assert!(!rec.seq.is_empty());
        assert_eq!(rec.seq.len(), rec.qual.len());
        // BS-seq reads are 65bp uniform in this fixture (SRR24827378 family)
        assert_eq!(rec.seq.len(), 65);
        Ok(())
    }

    #[test]
    fn qual_is_sanger_phred33() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_test.bam")?;
        let rec = r.next_record()?.expect("fixture has ≥1 record");
        // All qual chars must be in the printable Sanger range 33–126
        for c in rec.qual.chars() {
            let b = c as u8;
            assert!(
                (33..=126).contains(&b),
                "qual byte 0x{:02X} ('{}') outside printable Sanger range",
                b,
                c
            );
        }
        Ok(())
    }

    #[test]
    fn seq_contains_only_acgtn() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_test.bam")?;
        while let Some(rec) = r.next_record()? {
            for c in rec.seq.chars() {
                assert!(
                    matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'),
                    "non-ACGTN base in fixture output: '{}'",
                    c
                );
            }
        }
        Ok(())
    }

    #[test]
    fn paired_interleaved_de_interleaves_fixture() -> Result<()> {
        // The committed paired fixture is `samtools import`-emitted with
        // strict mate-adjacent ordering (R1 flag 77, R2 flag 141, repeating).
        // The de-interleaver should drain it cleanly with MAX_SLACK never
        // exceeded — 10 R1 / 10 R2 in lockstep → 10 matched pairs.
        let (mut r1_reader, mut r2_reader) =
            BamReader::open_paired_interleaved("test_files/ubam_paired_test.bam")?;

        let mut r1_records = Vec::new();
        let mut r2_records = Vec::new();
        // Pull from BOTH sides alternately — this is the access pattern
        // parallel.rs::read_pairs_round_robin uses.
        loop {
            let r1 = r1_reader.next_record()?;
            let r2 = r2_reader.next_record()?;
            match (r1, r2) {
                (Some(a), Some(b)) => {
                    // Same template name (mate-adjacent invariant).
                    assert_eq!(
                        a.id, b.id,
                        "paired-interleaved must emit matched-name pairs"
                    );
                    r1_records.push(a);
                    r2_records.push(b);
                }
                (None, None) => break,
                _ => panic!("paired readers desynchronised at EOF"),
            }
        }
        assert_eq!(r1_records.len(), 10, "expected 10 R1 records");
        assert_eq!(r2_records.len(), 10, "expected 10 R2 records");
        Ok(())
    }

    #[test]
    fn threaded_reads_match_direct() -> Result<()> {
        let mut direct = BamReader::open("test_files/ubam_test.bam")?;
        let mut threaded = BamReader::open_threaded("test_files/ubam_test.bam")?;
        loop {
            let d = direct.next_record()?;
            let t = threaded.next_record()?;
            match (d, t) {
                (None, None) => break,
                (Some(a), Some(b)) => {
                    assert_eq!(a.id, b.id);
                    assert_eq!(a.seq, b.seq);
                    assert_eq!(a.qual, b.qual);
                }
                _ => panic!("direct and threaded readers diverged at EOF"),
            }
        }
        Ok(())
    }

    #[test]
    fn paired_records_have_r1_then_r2_order() -> Result<()> {
        let mut r = BamReader::open("test_files/ubam_paired_test.bam")?;
        // First record should have READ1 flag → it'll be processed as-is via
        // next_record (no de-interleave happening here). Second should be the
        // matching R2 (same name). Verifies mate-adjacent ordering of the
        // fixture matches Spike 2's empirical finding for samtools import.
        let r1 = r.next_record()?.unwrap();
        let r2 = r.next_record()?.unwrap();
        assert_eq!(
            r1.id, r2.id,
            "fixture is mate-adjacent (same template name)"
        );
        Ok(())
    }

    // Aligned-BAM rejection cannot be tested against the committed fixture
    // (which is correctly all-unmapped). It is tested via integration test in
    // tests/integration_ubam.rs once we have an aligned-fixture pipeline.
    // Empty-seq / IUPAC / qual-mismatch / non-ACGTN tests would require
    // hand-built BAM bytes — deferred per IMPL.md "Notes for the
    // code-implementation agent" fallback (use committed fixture for happy
    // path; error paths covered at integration level).

    // ─── BamWriter + §3.6 parser tests (Phase 1 step 2.5) ──────────────────

    #[test]
    fn parse_name_only_no_tail() -> Result<()> {
        let (name, data) = parse_name_and_data("@SRR123.1")?;
        assert_eq!(name, b"SRR123.1");
        assert!(data.get(&Tag::new(b'C', b'B')).is_none());
        Ok(())
    }

    #[test]
    fn parse_name_strips_leading_at_sign_idempotently() -> Result<()> {
        let (n1, _) = parse_name_and_data("@read1")?;
        let (n2, _) = parse_name_and_data("read1")?;
        assert_eq!(n1, n2);
        assert_eq!(n1, b"read1");
        Ok(())
    }

    #[test]
    fn parse_name_splits_on_space_and_discards_annotation() -> Result<()> {
        // Real-world Illumina FASTQ headers: `@SRR.N <description>`.
        // The BAM spec rejects whitespace in QNAME, so the name must be
        // truncated at the first space. The annotation after is FASTQ-only
        // metadata that does NOT round-trip (it's not a tag-tail produced
        // by BamReader). PLAN §3.3 step 1.
        let (name, data) =
            parse_name_and_data("@SRR24827378.1 A00690:224:H7GC5DRX2:1:2101:1081:1016/1")?;
        assert_eq!(name, b"SRR24827378.1");
        assert!(
            data.get(&Tag::new(b'C', b'B')).is_none(),
            "no tags should be extracted from a space-separated annotation"
        );
        Ok(())
    }

    #[test]
    fn parse_tab_separator_still_treated_as_tag_tail() -> Result<()> {
        // Regression guard: the SPACE-or-TAB split must NOT downgrade the
        // \t-separated tag tail behaviour from #317 / BamReader output.
        let (name, data) = parse_name_and_data("@SRR1\tCB:Z:ATCG\tNM:i:5")?;
        assert_eq!(name, b"SRR1");
        assert!(data.get(&Tag::new(b'C', b'B')).is_some());
        assert!(data.get(&Tag::new(b'N', b'M')).is_some());
        Ok(())
    }

    #[test]
    fn parse_z_tag_recovered_as_string() -> Result<()> {
        let (name, data) = parse_name_and_data("@read1\tCB:Z:ATCGATCG-1")?;
        assert_eq!(name, b"read1");
        let val = data
            .get(&Tag::new(b'C', b'B'))
            .expect("CB tag missing from parsed data");
        let BufValue::String(s) = val else {
            panic!("expected Z (BufValue::String) variant");
        };
        let bytes: &[u8] = s.as_ref();
        assert_eq!(bytes, b"ATCGATCG-1");
        Ok(())
    }

    #[test]
    fn parse_i_tag_recovered_as_int32() -> Result<()> {
        let (_, data) = parse_name_and_data("@read1\tNM:i:42")?;
        let val = data.get(&Tag::new(b'N', b'M')).expect("NM tag missing");
        let BufValue::Int32(n) = val else {
            panic!("expected i (BufValue::Int32) variant");
        };
        assert_eq!(*n, 42);
        Ok(())
    }

    #[test]
    fn parse_f_tag_recovered_as_float() -> Result<()> {
        let (_, data) = parse_name_and_data("@read1\tXF:f:0.5")?;
        let val = data.get(&Tag::new(b'X', b'F')).expect("XF tag missing");
        let BufValue::Float(f) = val else {
            panic!("expected f (BufValue::Float) variant");
        };
        assert!((f - 0.5).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn parse_a_tag_recovered_as_character() -> Result<()> {
        let (_, data) = parse_name_and_data("@read1\tXA:A:+")?;
        let val = data.get(&Tag::new(b'X', b'A')).expect("XA tag missing");
        let BufValue::Character(c) = val else {
            panic!("expected A (BufValue::Character) variant");
        };
        assert_eq!(*c, b'+');
        Ok(())
    }

    #[test]
    fn parse_multiple_tags_all_recovered() -> Result<()> {
        let (_, data) = parse_name_and_data("@read1\tCB:Z:ATCG\tNM:i:5\tXF:f:1.5")?;
        assert!(data.get(&Tag::new(b'C', b'B')).is_some());
        assert!(data.get(&Tag::new(b'N', b'M')).is_some());
        assert!(data.get(&Tag::new(b'X', b'F')).is_some());
        Ok(())
    }

    #[test]
    fn parse_b_array_tag_rejected() {
        let err = parse_name_and_data("@read1\tXA:B:i,1,2,3").unwrap_err();
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("array") || msg.contains("B:"),
            "expected B-tag rejection, got: {}",
            msg
        );
    }

    #[test]
    fn parse_h_hex_tag_rejected() {
        let err = parse_name_and_data("@read1\tXH:H:DEADBEEF").unwrap_err();
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("hex") || msg.contains("H:"),
            "expected H-tag rejection, got: {}",
            msg
        );
    }

    #[test]
    fn read_side_emitter_rejects_hex_tag_value() {
        // Code-review A-I2 + round-2 R2-AGREE-1 regression guard: the
        // read-side `append_tag_type_and_value` must bail on `Value::Hex`
        // (symmetric with `Value::Array`). Reverting this fix would
        // silently emit `H:` into the FASTQ tag-tail, only for the
        // writer's `parse_tag_value` to abort the run mid-stream — the
        // asymmetric mid-stream-fatal failure mode round 1 caught.
        let hex_payload = bstr::BStr::new(b"DEADBEEF");
        let value = Value::Hex(hex_payload);
        let mut buf = String::new();
        let err = append_tag_type_and_value(&mut buf, &value)
            .expect_err("Value::Hex must be rejected by the read-side emitter");
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("hex") || msg.contains("H:"),
            "expected hex-tag rejection from read-side emitter, got: {}",
            msg
        );
    }

    #[test]
    fn qname_boundary_bytes_rejected() {
        for (label, name) in [
            ("space", b"read one".as_slice()),
            ("tab", b"read\tone"),
            ("LF", b"read\none"),
            ("CR", b"read\rone"),
            ("FF", b"read\x0Cone"),
        ] {
            assert!(
                reject_whitespace_in_qname(name).is_err(),
                "{} must be rejected in a QNAME",
                label
            );
        }
    }

    #[test]
    fn qname_non_boundary_bytes_accepted() {
        // VT is accepted by decision, not by derivation: widening the set here would
        // over-reject files that work today.
        for (label, name) in [
            ("plain", b"read_one".as_slice()),
            ("VT", b"read\x0Bone"),
            ("at-sign", b"read@one"),
            ("non-ASCII", "read\u{e9}one".as_bytes()),
        ] {
            assert!(
                reject_whitespace_in_qname(name).is_ok(),
                "{} must be accepted in a QNAME",
                label
            );
        }
    }

    #[test]
    fn qname_error_escapes_invisible_bytes() {
        let err = reject_whitespace_in_qname(b"read\tone").expect_err("tab must be rejected");
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("read\\tone"),
            "error must show the offending byte escaped, got: {}",
            msg
        );
    }

    #[test]
    fn tag_value_framing_bytes_rejected() {
        for (label, value) in [
            ("tab", "AAA\tCCC"),
            ("LF", "AAA\nCCC"),
            ("CR", "AAA\rCCC"),
            ("FF", "AAA\x0CCCC"),
        ] {
            assert!(
                reject_framing_in_tag_value(value).is_err(),
                "{} must be rejected in a Z tag value",
                label
            );
        }
    }

    #[test]
    fn tag_value_space_accepted() {
        // `Z` is `[ !-~]*`, so a space is legal — and the tab introducing the tag
        // precedes it, leaving `parse_name_and_data` able to split the tail.
        assert!(reject_framing_in_tag_value("has a space").is_ok());
    }

    #[test]
    fn read_side_emitter_rejects_newline_in_character_tag() {
        // `A` writes its raw byte into the same id as `Z`, so it needs the same guard.
        let value = Value::Character(b'\n');
        let mut buf = String::new();
        let err = append_tag_type_and_value(&mut buf, &value)
            .expect_err("a newline in an A value must be rejected by the read-side emitter");
        assert!(
            format!("{:#}", err).contains("cannot pass through the FASTQ header"),
            "expected the framing rejection, got: {:#}",
            err
        );
    }

    #[test]
    fn read_side_emitter_accepts_printable_character_tag() {
        let mut buf = String::new();
        append_tag_type_and_value(&mut buf, &Value::Character(b'+')).unwrap();
        assert_eq!(buf, "A:+");
    }

    #[test]
    fn read_side_emitter_rejects_newline_in_string_tag() {
        let value = Value::String(bstr::BStr::new(b"AAA\nCCC"));
        let mut buf = String::new();
        let err = append_tag_type_and_value(&mut buf, &value)
            .expect_err("a newline in a Z value must be rejected by the read-side emitter");
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("AAA\\nCCC"),
            "expected the escaped value in the message, got: {}",
            msg
        );
    }

    #[test]
    fn parse_unknown_type_code_rejected() {
        let err = parse_name_and_data("@read1\tXX:Q:foo").unwrap_err();
        let msg = format!("{:#}", err);
        assert!(msg.contains("unknown"), "got: {}", msg);
    }

    #[test]
    fn parse_empty_name_rejected() {
        assert!(parse_name_and_data("@").is_err());
        assert!(parse_name_and_data("").is_err());
    }

    #[test]
    fn parse_malformed_tag_missing_type_rejected() {
        assert!(parse_name_and_data("@read1\tCB").is_err());
        assert!(parse_name_and_data("@read1\tCB:").is_err());
    }

    #[test]
    fn parse_tolerates_trailing_empty_field() -> Result<()> {
        // A stray trailing `\t` should not break parsing — useful for
        // future tools that might emit a separator at the end.
        let (name, _) = parse_name_and_data("@read1\tCB:Z:ATCG\t")?;
        assert_eq!(name, b"read1");
        Ok(())
    }

    #[test]
    fn bam_writer_se_fastq_input_round_trip() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("se_out.bam");
        let records = [
            FastqRecord {
                id: "@read1".to_string(),
                seq: "ACGTACGTAC".to_string(),
                qual: "!!!!!!!!!!".to_string(),
            },
            FastqRecord {
                id: "@read2".to_string(),
                seq: "GGGCCCAATT".to_string(),
                qual: "IIIIIIIIII".to_string(),
            },
            FastqRecord {
                id: "@read3".to_string(),
                seq: "NNNNAAAATT".to_string(),
                qual: "5555AAAAII".to_string(),
            },
        ];
        let mut w = BamWriter::create(&out, None, &[], "trim_galore se test", 33)?;
        for r in &records {
            w.write_record(r, None)?;
        }
        w.finish()?;

        let mut r = BamReader::open(&out)?;
        for orig in &records {
            let got = r
                .next_record()?
                .expect("missing record on round-trip read-back");
            assert_eq!(got.id, orig.id, "id round-trip");
            assert_eq!(got.seq, orig.seq, "seq round-trip");
            assert_eq!(got.qual, orig.qual, "qual round-trip");
        }
        assert!(
            r.next_record()?.is_none(),
            "extra records on round-trip read-back"
        );
        Ok(())
    }

    #[test]
    fn bam_writer_pe_interleaved_round_trip() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("pe_interleaved.bam");
        let r1_records = [
            FastqRecord {
                id: "@pair1".to_string(),
                seq: "ACGTACGT".to_string(),
                qual: "IIIIIIII".to_string(),
            },
            FastqRecord {
                id: "@pair2".to_string(),
                seq: "GGGGCCCC".to_string(),
                qual: "AAAAAAAA".to_string(),
            },
        ];
        let r2_records = [
            FastqRecord {
                id: "@pair1".to_string(),
                seq: "TTTTGGGG".to_string(),
                qual: "55555555".to_string(),
            },
            FastqRecord {
                id: "@pair2".to_string(),
                seq: "AAAAATTT".to_string(),
                qual: "FFFFFFFF".to_string(),
            },
        ];
        let mut w = BamWriter::create(&out, None, &[], "trim_galore --paired test", 33)?;
        for (r1, r2) in r1_records.iter().zip(r2_records.iter()) {
            w.write_record(r1, Some(1))?;
            w.write_record(r2, Some(2))?;
        }
        w.finish()?;

        let (mut r1r, mut r2r) = BamReader::open_paired_interleaved(&out)?;
        for (orig_r1, orig_r2) in r1_records.iter().zip(r2_records.iter()) {
            let got1 = r1r.next_record()?.expect("missing R1 on read-back");
            let got2 = r2r.next_record()?.expect("missing R2 on read-back");
            assert_eq!(got1.id, orig_r1.id);
            assert_eq!(got1.seq, orig_r1.seq);
            assert_eq!(got1.qual, orig_r1.qual);
            assert_eq!(got2.id, orig_r2.id);
            assert_eq!(got2.seq, orig_r2.seq);
            assert_eq!(got2.qual, orig_r2.qual);
        }
        assert!(r1r.next_record()?.is_none());
        assert!(r2r.next_record()?.is_none());
        Ok(())
    }

    #[test]
    fn bam_writer_flag_bits_correct() -> Result<()> {
        // Verify the §3.3 step 2 flag values land in the file as raw bits.
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("flags.bam");
        let mut w = BamWriter::create(&out, None, &[], "trim_galore flag test", 33)?;
        w.write_record(
            &FastqRecord {
                id: "@se_read".to_string(),
                seq: "AC".to_string(),
                qual: "II".to_string(),
            },
            None,
        )?;
        w.write_record(
            &FastqRecord {
                id: "@pair".to_string(),
                seq: "AC".to_string(),
                qual: "II".to_string(),
            },
            Some(1),
        )?;
        w.write_record(
            &FastqRecord {
                id: "@pair".to_string(),
                seq: "GT".to_string(),
                qual: "II".to_string(),
            },
            Some(2),
        )?;
        w.finish()?;

        // Inspect raw flag bits via the bam crate (BamReader hides flags).
        let file = File::open(&out)?;
        let mut reader = bam::io::Reader::new(BufReader::new(file));
        let _h = reader.read_header()?;
        let mut rec = bam::Record::default();
        assert!(reader.read_record(&mut rec)? > 0);
        assert_eq!(rec.flags().bits(), 0x04, "SE flag must be 0x04 (FUNMAP)");
        assert!(reader.read_record(&mut rec)? > 0);
        assert_eq!(rec.flags().bits(), 0x4D, "R1 flag must be 0x4D = 77");
        assert!(reader.read_record(&mut rec)? > 0);
        assert_eq!(rec.flags().bits(), 0x8D, "R2 flag must be 0x8D = 141");
        assert_eq!(
            reader.read_record(&mut rec)?,
            0,
            "no extra records on flag-bits read-back"
        );
        Ok(())
    }

    #[test]
    fn bam_writer_aux_preservation_round_trip() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("aux.bam");
        let rec = FastqRecord {
            id: "@SRR1\tCB:Z:ATCGATCG-1\tUB:Z:AAAGGGCCC".to_string(),
            seq: "ACGTACGT".to_string(),
            qual: "IIIIIIII".to_string(),
        };
        let mut w = BamWriter::create(
            &out,
            None,
            &["CB".to_string(), "UB".to_string()],
            "trim_galore --preserve-tags CB,UB test",
            33,
        )?;
        w.write_record(&rec, None)?;
        w.finish()?;

        let mut r =
            BamReader::open(&out)?.with_preserved_tags(&["CB".to_string(), "UB".to_string()]);
        let got = r.next_record()?.expect("record missing on aux read-back");
        assert!(
            got.id.contains("CB:Z:ATCGATCG-1"),
            "CB tag not round-tripped: {:?}",
            got.id
        );
        assert!(
            got.id.contains("UB:Z:AAAGGGCCC"),
            "UB tag not round-tripped: {:?}",
            got.id
        );
        assert_eq!(got.seq, rec.seq);
        assert_eq!(got.qual, rec.qual);
        Ok(())
    }

    #[test]
    fn bam_writer_aux_typed_int_and_float_round_trip() -> Result<()> {
        // §3.6 explicitly calls out type fidelity for i / f tags. Verify
        // they land in the BAM as a real numeric type, NOT silently
        // downgraded to Z strings.
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("aux_typed.bam");
        let rec = FastqRecord {
            id: "@read1\tNM:i:42\tXF:f:0.5".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        let mut w = BamWriter::create(&out, None, &[], "trim_galore typed test", 33)?;
        w.write_record(&rec, None)?;
        w.finish()?;

        // Inspect typed Data at the BAM-record layer.
        let file = File::open(&out)?;
        let mut reader = bam::io::Reader::new(BufReader::new(file));
        let _h = reader.read_header()?;
        let mut bam_rec = bam::Record::default();
        assert!(reader.read_record(&mut bam_rec)? > 0);
        let data = bam_rec.data();

        let nm = data
            .get(&Tag::new(b'N', b'M'))
            .expect("NM missing")
            .expect("NM unreadable");
        let is_int = matches!(
            nm,
            Value::Int8(_)
                | Value::UInt8(_)
                | Value::Int16(_)
                | Value::UInt16(_)
                | Value::Int32(_)
                | Value::UInt32(_)
        );
        assert!(
            is_int,
            "NM should be a BAM integer type, not Z/string-encoded"
        );

        let xf = data
            .get(&Tag::new(b'X', b'F'))
            .expect("XF missing")
            .expect("XF unreadable");
        assert!(
            matches!(xf, Value::Float(_)),
            "XF should be BAM Float, not Z/string-encoded"
        );
        Ok(())
    }

    /// Issue #358 code-review follow-up — the read-side raw→ASCII conversion
    /// must emit exactly one byte per base.
    ///
    /// A mixed-sentinel `QUAL` (some real scores, some `0xFF`) bypasses the
    /// all-`0xFF` fast path. An earlier revision used `saturating_add`, which
    /// turned `0xFF` into `255u8 as char` = U+00FF — **two** UTF-8 bytes — so
    /// `qual.len() != seq.len()` and the byte-indexed helpers
    /// (`FastqRecord::truncate`, `clip_5prime`, `clip_3prime`) would panic on a
    /// char boundary or silently keep the wrong number of scores.
    ///
    /// Pins the invariant directly: byte length equals base count, and every
    /// byte is printable ASCII.
    #[test]
    fn qual_ascii_conversion_is_single_byte_per_base() {
        // Exercised through the same expression `bam_record_to_fastq` uses.
        let qual_raw: Vec<u8> = vec![40, QUAL_MISSING_SENTINEL, 93, 200, 0];
        let qual: String = qual_raw
            .iter()
            .map(|&b| {
                if b == QUAL_MISSING_SENTINEL {
                    QUAL_MISSING_REPLACEMENT as char
                } else {
                    (b.min(MAX_PHRED_SCORE) + PHRED_OFFSET) as char
                }
            })
            .collect();

        assert_eq!(
            qual.len(),
            qual_raw.len(),
            "qual must be one BYTE per base, got {} bytes for {} bases: {:?}",
            qual.len(),
            qual_raw.len(),
            qual.as_bytes()
        );
        assert!(
            qual.is_ascii(),
            "every quality byte must be single-byte ASCII, got {:?}",
            qual.as_bytes()
        );
        // 40→'I', sentinel→'!', 93→'~', 200 clamped to 93→'~', 0→'!'
        assert_eq!(qual, "I!~~!");
        // Byte-indexed slicing must not panic — the operation that regressed.
        assert_eq!(&qual[1..], "!~~!");
    }

    /// Issue #358 — the writer must subtract the *input's* ASCII offset, not a
    /// hardcoded 33. Both branches are pinned here: `'h'` (104) is Q40 under
    /// Phred+64 and `'I'` (73) is Q40 under Phred+33, so a correct writer
    /// stores raw 40 in both cases.
    ///
    /// Read back at the raw BAM-record layer rather than through `BamReader`,
    /// which re-adds `PHRED_OFFSET` and would mask the value under test. Same
    /// approach as `bam_writer_aux_typed_int_and_float_round_trip` above.
    #[test]
    fn bam_writer_subtracts_input_phred_offset() -> Result<()> {
        fn stored_qual(qual: &str, input_phred_offset: u8, name: &str) -> Result<Vec<u8>> {
            let dir = tempfile::tempdir()?;
            let out = dir.path().join(format!("{name}.bam"));
            let rec = FastqRecord {
                id: format!("@{name}"),
                seq: "ACGT".to_string(),
                qual: qual.to_string(),
            };
            let mut w = BamWriter::create(
                &out,
                None,
                &[],
                "trim_galore offset test",
                input_phred_offset,
            )?;
            w.write_record(&rec, None)?;
            w.finish()?;

            let file = File::open(&out)?;
            let mut reader = bam::io::Reader::new(BufReader::new(file));
            let _h = reader.read_header()?;
            let mut bam_rec = bam::Record::default();
            assert!(reader.read_record(&mut bam_rec)? > 0, "no record written");
            Ok(bam_rec.quality_scores().as_ref().to_vec())
        }

        // Phred+64: 'h' = 104 = Q40. Pre-fix this stored 104-33 = 71.
        assert_eq!(
            stored_qual("hhhh", 64, "p64")?,
            vec![40u8; 4],
            "Phred+64 input must store raw Phred 40, not 71"
        );

        // Phred+33 counterpart — an over-correction guard, not a bug guard:
        // this already passed before the fix and must keep passing.
        assert_eq!(
            stored_qual("IIII", 33, "p33")?,
            vec![40u8; 4],
            "Phred+33 input must still store raw Phred 40"
        );
        Ok(())
    }

    #[test]
    fn validate_seq_accepts_acgtn_upper() -> Result<()> {
        let out = validate_and_normalize_seq_for_write(b"ACGTN")?;
        assert_eq!(out, b"ACGTN");
        Ok(())
    }

    #[test]
    fn validate_seq_uppercases_lowercase() -> Result<()> {
        let out = validate_and_normalize_seq_for_write(b"acgtn")?;
        assert_eq!(out, b"ACGTN");
        Ok(())
    }

    #[test]
    fn validate_seq_coerces_iupac_to_n() -> Result<()> {
        // PLAN §3.3 step 3 — IUPAC degenerate bases coerce to N.
        // Code-review B-M1 regression guard.
        let out = validate_and_normalize_seq_for_write(b"ARYNKW")?;
        assert_eq!(out, b"ANNNNN");
        Ok(())
    }

    #[test]
    fn validate_seq_rejects_equals_sign() {
        let err = validate_and_normalize_seq_for_write(b"AC=T").unwrap_err();
        let msg = format!("{:#}", err);
        assert!(
            msg.contains("=") || msg.contains("match-reference"),
            "expected '=' rejection, got: {}",
            msg
        );
    }

    #[test]
    fn validate_seq_rejects_garbage_byte() {
        let err = validate_and_normalize_seq_for_write(b"AC?T").unwrap_err();
        let msg = format!("{:#}", err);
        assert!(msg.contains("unrecognised") || msg.contains("0x3F"));
    }

    #[test]
    fn build_output_header_disambiguates_existing_trim_galore_pg() -> Result<()> {
        // Re-processing a trim_galore-output BAM must NOT overwrite the
        // existing @PG; noodles' Programs::add() yields a disambiguated
        // ID `trim_galore-trim_galore` with PP pointing at the original.
        // Code-review AGREE-1 regression guard.
        let source = Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_program(
                "trim_galore",
                Map::<Program>::builder()
                    .insert(program_tag::NAME, "trim_galore")
                    .insert(program_tag::VERSION, "0.0.0-test")
                    .insert(program_tag::COMMAND_LINE, "trim_galore prior-run")
                    .build()
                    .unwrap(),
            )
            .build();

        let result = build_output_header(Some(&source), "trim_galore new-run")?;
        let ids: Vec<Vec<u8>> = result
            .programs()
            .as_ref()
            .keys()
            .map(|k| AsRef::<[u8]>::as_ref(k).to_vec())
            .collect();
        assert_eq!(ids.len(), 2, "expected two @PG entries: {:?}", ids);
        assert_eq!(ids[0], b"trim_galore", "original @PG ID preserved");
        assert_eq!(
            ids[1], b"trim_galore-trim_galore",
            "second @PG disambiguated"
        );
        Ok(())
    }

    #[test]
    fn bam_writer_synthesises_header_for_fastq_input() -> Result<()> {
        // No source header → BamWriter should still produce a valid BAM
        // file with a trim_galore @PG line in the synthesised header.
        let dir = tempfile::tempdir()?;
        let out = dir.path().join("synth_header.bam");
        let mut w = BamWriter::create(
            &out,
            None,
            &[],
            "trim_galore --output-format ubam input.fq",
            33,
        )?;
        w.write_record(
            &FastqRecord {
                id: "@r1".to_string(),
                seq: "AC".to_string(),
                qual: "II".to_string(),
            },
            None,
        )?;
        w.finish()?;

        // Re-open and read the header to verify trim_galore @PG is present.
        let file = File::open(&out)?;
        let mut reader = bam::io::Reader::new(BufReader::new(file));
        let header = reader.read_header()?;
        let programs = header.programs();
        let has_tg = programs
            .as_ref()
            .keys()
            .any(|k| AsRef::<[u8]>::as_ref(k) == b"trim_galore");
        assert!(
            has_tg,
            "@PG ID:trim_galore must appear in synthesised header"
        );
        Ok(())
    }
}
