//! FASTQ record type and streaming I/O.
//!
//! Provides gzip-aware reading and writing of FASTQ records with 64KB buffers
//! for efficient I/O throughput.

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use gzp::ZWriter;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Buffer size for gzip I/O — 64KB for throughput (flate2 default of 8KB is too small).
const BUF_SIZE: usize = 64 * 1024;

/// Default output gzip compression level. Set to 1 (fastest) — at
/// Buckberry-scale (84M reads, 38%
/// adapter rate, cores=8) the compression CPU dominated wall time on
/// saturated workers; lowering from level 6 to 1 measured ~−23% wall
/// and ~−43% user-CPU at byte-identity of the decompressed output
/// (gzip framing differs but `gzip -dc` yields the same bytes). Trade:
/// output `.fq.gz` files are ~75% larger. Users who care about storage
/// pass `--compression 6` (or `--compression 9`) — and stack
/// `--clumpify` on top for reordering plus a higher gzip level.
pub const DEFAULT_GZIP_LEVEL: u32 = 1;

/// Content-based gzip check for the constructors that have no caller-supplied
/// verdict.
///
/// Reads three bytes and looks for the gzip magic. That is cheap enough that
/// guessing from the filename buys nothing, and the filename is wrong in both
/// directions: `bgzip` output is commonly named `.bgz`, and a plain FASTQ is
/// sometimes misnamed `.fastq.gz`.
///
/// Deliberately contains no BAM logic, so `fastq` does not gain a dependency on
/// `format` and the `BAM\1` discrimination stays in `detect_input_format` where
/// it belongs. It cannot misfire on plain FASTQ, whose first byte is `@` and
/// never `0x1F`.
///
/// Also deliberately not `detect_input_format`, which *bails* on empty input.
/// `demux` reads back a trimmed output that is legitimately empty when every
/// read was filtered, and that must not become an error. Any failure here
/// (missing file, permissions, short file) answers "not gzip" and lets the
/// subsequent open report the real problem.
fn sniff_gzip(path: &Path) -> bool {
    let mut magic = [0u8; 3];
    File::open(path)
        .and_then(|mut f| std::io::Read::read(&mut f, &mut magic))
        .is_ok_and(|n| n == 3 && magic == [0x1F, 0x8B, 0x08])
}

/// A single FASTQ record with owned data.
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Header line including '@' prefix and any description (no trailing newline)
    pub id: String,
    /// Nucleotide sequence (no trailing newline)
    pub seq: String,
    /// Quality scores as ASCII (no trailing newline)
    pub qual: String,
}

impl FastqRecord {
    /// Returns the sequence length.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Returns true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Write this record to a FASTQ writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        // Single buffered write per record. Pre-format the entire 4-line
        // FASTQ block into a Vec<u8>, then issue one `write_all` instead
        // of four `writeln!` calls. Byte-identical output to the previous
        // form (same bytes, same order, same trailing `\n`s); the win is
        // amortising the per-call overhead — at Buckberry scale (84M
        // reads, 38% adapter rate) this is ~10% wall-clock at cores=8.
        // See #248 (item #2 in @an-altosian's perf audit).
        let mut buf = Vec::with_capacity(self.id.len() + self.seq.len() + self.qual.len() + 5);
        buf.extend_from_slice(self.id.as_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(self.seq.as_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(b"+\n");
        buf.extend_from_slice(self.qual.as_bytes());
        buf.push(b'\n');
        writer.write_all(&buf)?;
        Ok(())
    }

    /// Trim sequence and quality to `new_len` bases from the 5' end.
    /// Keeps the first `new_len` bases.
    pub fn truncate(&mut self, new_len: usize) {
        if new_len < self.seq.len() {
            self.seq.truncate(new_len);
            self.qual.truncate(new_len);
        }
    }

    /// Remove `n` bases from the 5' (left) end.
    /// Returns the clipped sequence for --rename support.
    pub fn clip_5prime(&mut self, n: usize) -> Option<String> {
        if n > 0 && self.seq.len() > n {
            let clipped = self.seq[..n].to_string();
            self.seq = self.seq[n..].to_string();
            self.qual = self.qual[n..].to_string();
            Some(clipped)
        } else {
            None
        }
    }

    /// Remove `n` bases from the 3' (right) end.
    /// Returns the clipped sequence for --rename support.
    pub fn clip_3prime(&mut self, n: usize) -> Option<String> {
        if n > 0 && self.seq.len() > n {
            let trim_to = self.seq.len() - n;
            let clipped = self.seq[trim_to..].to_string();
            self.seq.truncate(trim_to);
            self.qual.truncate(trim_to);
            Some(clipped)
        } else {
            None
        }
    }

    /// Count the number of 'N' bases in the sequence.
    pub fn n_count(&self) -> usize {
        self.seq.bytes().filter(|&b| b == b'N' || b == b'n').count()
    }

    /// Trim leading and trailing N bases from both ends (--trim-n).
    pub fn trim_ns(&mut self) {
        let seq_bytes = self.seq.as_bytes();
        let qual_bytes = self.qual.as_bytes();

        // Find first non-N from left
        let start = seq_bytes
            .iter()
            .position(|&b| b != b'N' && b != b'n')
            .unwrap_or(seq_bytes.len());

        // Find first non-N from right
        let end = seq_bytes
            .iter()
            .rposition(|&b| b != b'N' && b != b'n')
            .map(|p| p + 1)
            .unwrap_or(0);

        if start > 0 || end < seq_bytes.len() {
            let end = end.max(start); // handle all-N case
            self.seq = String::from_utf8_lossy(&seq_bytes[start..end]).to_string();
            self.qual = String::from_utf8_lossy(&qual_bytes[start..end]).to_string();
        }
    }

    /// Append clipping info to the read ID for --rename support.
    ///
    /// If the id carries a `\t`-separated tag tail (from a BamReader-fed
    /// uBAM input, where the format is `@NAME\tTAG:TYPE:VALUE\t...`), the
    /// suffix is spliced into the NAME portion BEFORE the first tab so
    /// the tail stays intact for the BamWriter round-trip. Otherwise the
    /// `:clip5:...` suffix would land inside the last tag's value and
    /// silently corrupt it. Caught by code-review finding C1 against the
    /// PLAN v2.1 §3.6 round-trip mechanism.
    pub fn append_to_id(&mut self, suffix: &str) {
        // Strip any trailing whitespace/newline from id before appending.
        let id = self.id.trim_end();
        if let Some(tab_pos) = id.find('\t') {
            self.id = format!("{}{}{}", &id[..tab_pos], suffix, &id[tab_pos..]);
        } else {
            self.id = format!("{}{}", id, suffix);
        }
    }
}

/// Extract the ID prefix used for paired-end / `--passthrough` sync checks.
///
/// Steps:
///   1. Strip a leading `@` if present (`FastqRecord::id` stores it).
///   2. Take everything before the first ASCII whitespace.
///   3. Strip a trailing `/1`, `/2`, or `/3` (legacy Illumina ID style — still
///      common in SRA / ENA archived data).
///
/// Examples:
///   `"@read1"`                → `"read1"`
///   `"@read1 1:N:0:CGATCG"`   → `"read1"`
///   `"@read1/1"`              → `"read1"`
///   `"@read1/2"`              → `"read1"`
///   `"@read1/3"`              → `"read1"`
///   `"@read1/1 1:N:0:CGATCG"` → `"read1"`
///   `"@read1/4"`              → `"read1/4"`  (only `/1` `/2` `/3` strip)
///
/// Used by both the serial paired-end pipeline (`trimmer::run_paired_end`)
/// and the parallel reader thread (`parallel::read_pairs_round_robin`) to
/// verify that R1, R2, and the optional `--passthrough` records share the
/// same template across all three streams. Per plan v2 §AB2 the trailing
/// `/[123]` strip is in v1 (not deferred) — covers the legacy
/// `@read/1` / `@read/2` / `@read/3` convention without false-positive risk
/// (legitimate read IDs never end in a literal `/1` `/2` `/3` token).
pub fn read_id_prefix(id: &str) -> &str {
    let stripped = id.strip_prefix('@').unwrap_or(id);
    let head = stripped.split_ascii_whitespace().next().unwrap_or("");
    head.rsplit_once('/')
        .filter(|(_, suf)| matches!(*suf, "1" | "2" | "3"))
        .map(|(p, _)| p)
        .unwrap_or(head)
}

/// Number of records per batch sent through the threaded reader channel.
/// Larger batches reduce channel overhead; 4096 records ≈ 1.5MB.
const READER_BATCH_SIZE: usize = 4096;

/// Number of batches buffered ahead in the channel.
/// 4 batches × 4096 records = ~16K records of lookahead.
const READER_CHANNEL_BATCHES: usize = 4;

/// Internal source for a FastqReader — either direct I/O or a channel
/// fed by a background decompression thread.
enum ReaderSource {
    Direct {
        reader: Box<dyn BufRead + Send>,
        line_buf: String,
    },
    Threaded {
        rx: std::sync::mpsc::Receiver<Result<Vec<FastqRecord>>>,
        buffer: Vec<FastqRecord>,
        buf_pos: usize,
        _handle: std::thread::JoinHandle<()>,
    },
}

/// Polymorphic input source — anything that can yield `FastqRecord`s.
///
/// Both `FastqReader` and `crate::bam::BamReader` implement this so the worker
/// pool in `parallel.rs` can dispatch through a single trait-object hand-off.
/// Spike 1 measured the dispatch overhead at +0.86–1.28% over a concrete-type
/// baseline on a 500K × 150bp parse loop — well under the 2% bar; see
/// `plans/06252026_ubam-input-support/spikes/SPIKE_recordsource.md`.
///
/// The `Send` bound is required because the reader is moved into the spawned
/// reader thread inside the parallel pipeline.
pub trait RecordSource: Send {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

impl RecordSource for FastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        FastqReader::next_record(self)
    }
}

/// A streaming FASTQ reader that handles both plain and gzipped files.
///
/// Supports two modes:
/// - **Direct** (`open`): reads synchronously on the calling thread.
/// - **Threaded** (`open_threaded`): spawns a background thread for
///   decompression, feeding records through a bounded channel. This
///   overlaps I/O decompression with trimming/compression on the main thread.
pub struct FastqReader {
    source: ReaderSource,
}

impl FastqReader {
    /// Open a FASTQ file for synchronous reading, deciding gzip from the
    /// filename.
    ///
    /// Prefer [`FastqReader::open_with`] when the caller already knows the
    /// format from file *content* (as `format::detect_input_format` does):
    /// the filename is the weaker signal and misses e.g. `.fq.bgz`.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        Self::open_with(path, sniff_gzip(path))
    }

    /// Open a FASTQ file for synchronous reading.
    ///
    /// `is_gzip` comes from the caller, normally `format::detect_input_format`,
    /// which inspects file content rather than the filename.
    pub fn open_with<P: AsRef<Path>>(path: P, is_gzip: bool) -> Result<Self> {
        let path = path.as_ref();
        let reader = Self::open_reader(path, is_gzip)?;

        Ok(FastqReader {
            source: ReaderSource::Direct {
                reader,
                line_buf: String::with_capacity(512),
            },
        })
    }

    /// Open a FASTQ file with background decompression on a dedicated thread,
    /// deciding gzip from the filename.
    ///
    /// Prefer [`FastqReader::open_threaded_with`] when the format is already
    /// known from file content.
    pub fn open_threaded<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        Self::open_threaded_with(path, sniff_gzip(path))
    }

    /// Open a FASTQ file with background decompression on a dedicated thread.
    ///
    /// Returns immediately. A background thread reads and decompresses
    /// the file, sending records through a bounded channel. The main
    /// thread calls `next_record()` which receives from the channel,
    /// overlapping decompression with processing.
    ///
    /// `is_gzip` comes from the caller; see [`FastqReader::open_with`].
    pub fn open_threaded_with<P: AsRef<Path>>(path: P, is_gzip: bool) -> Result<Self> {
        let path = path.as_ref().to_path_buf();

        // Validate the file exists before spawning the thread
        if !path.exists() {
            bail!("Input file not found: {}", path.display());
        }

        let (tx, rx) = std::sync::mpsc::sync_channel(READER_CHANNEL_BATCHES);

        let handle = std::thread::spawn(move || {
            let mut reader = match Self::open_direct(&path, is_gzip) {
                Ok(r) => r,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            };
            let mut batch = Vec::with_capacity(READER_BATCH_SIZE);
            loop {
                match Self::read_next_direct(&mut reader) {
                    Ok(Some(record)) => {
                        batch.push(record);
                        if batch.len() >= READER_BATCH_SIZE {
                            let full_batch = std::mem::replace(
                                &mut batch,
                                Vec::with_capacity(READER_BATCH_SIZE),
                            );
                            if tx.send(Ok(full_batch)).is_err() {
                                break; // receiver dropped
                            }
                        }
                    }
                    Ok(None) => {
                        // Send remaining records
                        if !batch.is_empty() {
                            let _ = tx.send(Ok(batch));
                        }
                        // Signal EOF with an empty batch
                        let _ = tx.send(Ok(Vec::new()));
                        break;
                    }
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        break;
                    }
                }
            }
        });

        Ok(FastqReader {
            source: ReaderSource::Threaded {
                rx,
                buffer: Vec::new(),
                buf_pos: 0,
                _handle: handle,
            },
        })
    }

    /// Internal: build the buffered byte reader for `path`.
    ///
    /// The single place that turns a gzip verdict into a decoder, so the
    /// verdict is applied identically on the sync and threaded paths.
    fn open_reader(path: &Path, is_gzip: bool) -> Result<Box<dyn BufRead + Send>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open input file: {}", path.display()))?;

        Ok(if is_gzip {
            Box::new(BufReader::with_capacity(
                BUF_SIZE,
                MultiGzDecoder::new(file),
            ))
        } else {
            Box::new(BufReader::with_capacity(BUF_SIZE, file))
        })
    }

    /// Internal: open a file and return the raw reader components (for use in threads).
    fn open_direct(path: &Path, is_gzip: bool) -> Result<(Box<dyn BufRead + Send>, String)> {
        let reader = Self::open_reader(path, is_gzip)?;
        Ok((reader, String::with_capacity(512)))
    }

    /// Internal: read one record from a direct reader (used by both Direct mode and threads).
    fn read_next_direct(
        state: &mut (Box<dyn BufRead + Send>, String),
    ) -> Result<Option<FastqRecord>> {
        let (reader, line_buf) = state;

        // Line 1: ID (starts with @)
        line_buf.clear();
        if reader.read_line(line_buf)? == 0 {
            return Ok(None); // EOF
        }
        let id = line_buf.trim_end_matches(['\n', '\r']).to_string();

        // Line 2: Sequence
        line_buf.clear();
        if reader.read_line(line_buf)? == 0 {
            bail!("Truncated FASTQ: missing sequence line after {}", id);
        }
        let seq = line_buf.trim_end_matches(['\n', '\r']).to_string();

        // Line 3: Plus line (discard)
        line_buf.clear();
        if reader.read_line(line_buf)? == 0 {
            bail!("Truncated FASTQ: missing '+' line after {}", id);
        }

        // Line 4: Quality
        line_buf.clear();
        if reader.read_line(line_buf)? == 0 {
            bail!("Truncated FASTQ: missing quality line after {}", id);
        }
        let qual = line_buf.trim_end_matches(['\n', '\r']).to_string();

        Ok(Some(FastqRecord { id, seq, qual }))
    }

    /// Read the next FASTQ record. Returns None at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        match &mut self.source {
            ReaderSource::Direct { reader, line_buf } => {
                // Line 1: ID (starts with @)
                line_buf.clear();
                if reader.read_line(line_buf)? == 0 {
                    return Ok(None); // EOF
                }
                let id = line_buf.trim_end_matches(['\n', '\r']).to_string();

                // Line 2: Sequence
                line_buf.clear();
                if reader.read_line(line_buf)? == 0 {
                    bail!("Truncated FASTQ: missing sequence line after {}", id);
                }
                let seq = line_buf.trim_end_matches(['\n', '\r']).to_string();

                // Line 3: Plus line (discard)
                line_buf.clear();
                if reader.read_line(line_buf)? == 0 {
                    bail!("Truncated FASTQ: missing '+' line after {}", id);
                }

                // Line 4: Quality
                line_buf.clear();
                if reader.read_line(line_buf)? == 0 {
                    bail!("Truncated FASTQ: missing quality line after {}", id);
                }
                let qual = line_buf.trim_end_matches(['\n', '\r']).to_string();

                Ok(Some(FastqRecord { id, seq, qual }))
            }
            ReaderSource::Threaded {
                rx,
                buffer,
                buf_pos,
                ..
            } => {
                // Return next record from current batch buffer
                if *buf_pos < buffer.len() {
                    let idx = *buf_pos;
                    *buf_pos += 1;
                    // Take the record out, replacing with a dummy to avoid clone
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
                // Buffer exhausted — receive next batch
                match rx.recv() {
                    Ok(Ok(batch)) => {
                        if batch.is_empty() {
                            Ok(None) // EOF signal
                        } else {
                            *buffer = batch;
                            *buf_pos = 1;
                            // Return first record from new batch
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
                    Err(_) => Ok(None), // channel closed
                }
            }
        }
    }

    /// Perform input sanity checks on the first record, deciding gzip from the
    /// filename.
    ///
    /// Prefer [`FastqReader::sanity_check_with`] when the format is already
    /// known from file content.
    pub fn sanity_check<P: AsRef<Path>>(path: P) -> Result<()> {
        let path = path.as_ref();
        Self::sanity_check_with(path, sniff_gzip(path))
    }

    /// Perform input sanity checks on the first record.
    /// Checks: FASTQ format validation, colorspace detection, empty file.
    ///
    /// This runs before anything else in `main()`, so it is the first place a
    /// wrong gzip verdict surfaces: on a `.fq.bgz` input it read the
    /// compressed bytes as text and failed with "stream did not contain valid
    /// UTF-8" before the reader factories were ever reached.
    pub fn sanity_check_with<P: AsRef<Path>>(path: P, is_gzip: bool) -> Result<()> {
        let path = path.as_ref();
        let mut reader = FastqReader::open_with(path, is_gzip)?;

        match reader.next_record()? {
            None => bail!(
                "Input file '{}' seems to be completely empty. Consider respecifying!",
                path.display()
            ),
            Some(record) => {
                // Check FASTQ format
                if !record.id.starts_with('@') {
                    bail!(
                        "Input file '{}' doesn't seem to be in FastQ format (first line doesn't start with '@')",
                        path.display()
                    );
                }

                // Check for colorspace (digits in sequence = SOLiD format)
                if record.seq.bytes().any(|b| b.is_ascii_digit()) {
                    bail!(
                        "File seems to be in SOLiD colorspace format which is not supported \
                         by Trim Galore (sequence is: '{}'). Colorspace data requires \
                         separate processing!",
                        record.seq
                    );
                }
            }
        }

        Ok(())
    }
}

/// The three concrete sinks a [`FastqWriter`] can write records into.
///
/// **An enum rather than `Box<dyn Write + Send>`, and the reason is teardown.**
/// Every sink here owes trailing bytes that `write`/`flush` never emit — the
/// gzip CRC/ISIZE trailer, the parallel compressor's final blocks — and the
/// call that emits them is inherent to the concrete type: neither
/// [`GzEncoder::try_finish`] nor [`gzp::ZWriter::finish`] is reachable through
/// a trait object. Behind a `Box<dyn Write>` the trailer was therefore written
/// by the sink's own `Drop`, which discards the error (`let _ =
/// self.try_finish()`) or panics — so an `ENOSPC` in the last few bytes of a
/// run produced a **truncated output file at its final name with exit 0**
/// (#434).
///
/// Naming the sinks makes the teardown callable, so [`FastqWriter::finish`]
/// can propagate its error *before* [`crate::io::PendingOutput::commit`]
/// publishes anything. Cost is a three-arm match per record in place of a
/// vtable dispatch.
enum Sink {
    Plain(BufWriter<File>),
    Gz(BufWriter<GzEncoder<File>>),
    ParGz(ParCompress<'static, Gzip, File>),
}

impl Sink {
    /// Close the sink, writing whatever trailer it owes, and return the error
    /// if that fails. Consumes `self`, so the sink's own `Drop` has nothing
    /// left to do — which is also how the success path sidesteps
    /// `ParCompress::drop`'s `.unwrap()`.
    fn finish(self) -> Result<()> {
        match self {
            // `into_inner` flushes the buffer and hands back the error if
            // that write fails, where `Drop` would have swallowed it.
            Sink::Plain(w) => {
                w.into_inner().map_err(|e| e.into_error())?;
            }
            Sink::Gz(w) => {
                let mut encoder = w.into_inner().map_err(|e| e.into_error())?;
                // The deflate end-of-stream plus the 8-byte CRC/ISIZE trailer.
                // `flush()` does NOT write these; only `try_finish` does.
                encoder.try_finish()?;
            }
            Sink::ParGz(mut w) => {
                w.finish()?;
            }
        }
        Ok(())
    }
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Sink::Plain(w) => w.write(buf),
            Sink::Gz(w) => w.write(buf),
            Sink::ParGz(w) => w.write(buf),
        }
    }

    /// Forwarded explicitly: `FastqRecord::write_to` issues exactly one
    /// `write_all` per record, and the default implementation would wrap it in
    /// a loop over `write` on the hot path.
    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match self {
            Sink::Plain(w) => w.write_all(buf),
            Sink::Gz(w) => w.write_all(buf),
            Sink::ParGz(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Sink::Plain(w) => w.flush(),
            Sink::Gz(w) => w.flush(),
            Sink::ParGz(w) => w.flush(),
        }
    }
}

/// A streaming FASTQ writer that handles both plain and gzipped output.
///
/// Records go to a temporary; [`FastqWriter::finish`] publishes them under the
/// final name. A writer that is dropped instead of finished leaves no output
/// file at all.
pub struct FastqWriter {
    /// Declared before `pending` so the compressor is closed — and its trailer
    /// written — before the temporary is removed on the abandon path.
    writer: Sink,
    pending: crate::io::PendingOutput,
}

impl FastqWriter {
    /// Create a new FASTQ writer. Gzip-compresses if `gzip` is true.
    /// When `cores` > 1 and gzip is enabled, uses parallel gzip compression.
    /// `gzip_level` is the deflate level (1 = fastest/largest, 6 = balanced,
    /// 9 = smallest/slowest). Default: `DEFAULT_GZIP_LEVEL` (1); the user
    /// overrides via `--compression`.
    pub fn create<P: AsRef<Path>>(
        path: P,
        gzip: bool,
        cores: usize,
        gzip_level: u32,
    ) -> Result<Self> {
        let path = path.as_ref();

        // Ensure parent directory exists
        if let Some(parent) = path.parent()
            && !parent.exists()
        {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("Failed to create output directory: {}", parent.display())
            })?;
        }

        let (pending, file) = crate::io::PendingOutput::create(path)?;

        let writer = if gzip {
            if cores > 1 {
                // Parallel gzip: split output into independently-compressed blocks
                Sink::ParGz(
                    ParCompressBuilder::<Gzip>::new()
                        .num_threads(cores)
                        .with_context(|| {
                            format!(
                                "Failed to create parallel compressor with {} threads",
                                cores
                            )
                        })?
                        .compression_level(Compression::new(gzip_level))
                        .from_writer(file),
                )
            } else {
                // Single-threaded gzip with zlib-rs SIMD backend
                Sink::Gz(BufWriter::with_capacity(
                    BUF_SIZE,
                    GzEncoder::new(file, Compression::new(gzip_level)),
                ))
            }
        } else {
            Sink::Plain(BufWriter::with_capacity(BUF_SIZE, file))
        };

        Ok(FastqWriter { writer, pending })
    }

    /// Write a FASTQ record.
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        record.write_to(&mut self.writer)
    }

    /// Close the writer and publish the output under its final name. Consumes
    /// `self` so a caller cannot forget.
    ///
    /// The sink is torn down *first* and its error propagated, so the final
    /// name appears only after every byte the output owes — including the
    /// gzip trailer — has been written successfully. A failed teardown returns
    /// here with `pending` still armed, so its `Drop` removes the temporary
    /// and nothing is published (#434). Before this, teardown ran in the
    /// sink's own `Drop` and its error was discarded.
    pub fn finish(self) -> Result<()> {
        let FastqWriter { writer, pending } = self;
        writer.finish()?;
        pending.commit()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A writer whose every write fails, which is what a full disk looks like
    /// from inside the sink. A `File` opened read-only is the shortest
    /// portable one — no size-limited filesystem required.
    fn unwritable_file(dir: &std::path::Path) -> File {
        let path = dir.join("read-only-sink");
        std::fs::write(&path, b"").unwrap();
        File::open(&path).unwrap()
    }

    #[test]
    fn gz_sink_finish_reports_a_failing_teardown() {
        // The bytes that fail here are the ones ONLY `try_finish` writes — the
        // deflate end-of-stream plus the 8-byte CRC/ISIZE trailer. Behind
        // `Box<dyn Write>` this error was written by `GzEncoder::drop` and
        // discarded (`let _ = self.try_finish()`), so a run that lost its
        // trailer still published the file and exited 0 (#434).
        let dir = tempfile::tempdir().unwrap();
        let sink = Sink::Gz(BufWriter::with_capacity(
            BUF_SIZE,
            GzEncoder::new(unwritable_file(dir.path()), Compression::new(1)),
        ));
        assert!(
            sink.finish().is_err(),
            "a gzip teardown that cannot write its trailer must not report success"
        );
    }

    #[test]
    fn plain_sink_finish_reports_a_failing_flush() {
        // The plain path was already covered by `BufWriter::flush` in the old
        // `finish`; this pins that `Sink::finish` did not lose it.
        let dir = tempfile::tempdir().unwrap();
        let mut sink = Sink::Plain(BufWriter::with_capacity(
            BUF_SIZE,
            unwritable_file(dir.path()),
        ));
        sink.write_all(b"@r\nACGT\n+\nIIII\n").unwrap(); // buffered, not yet written
        assert!(
            sink.finish().is_err(),
            "a buffered write that cannot reach the file must not report success"
        );
    }

    #[test]
    fn gz_sink_finish_writes_the_trailer_explicitly() {
        // Positive control for the above: the same call on a writable file
        // leaves a complete gzip member, so the teardown is doing the work
        // that `Drop` used to do rather than merely returning `Ok`.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("complete.gz");
        let sink = Sink::Gz(BufWriter::with_capacity(
            BUF_SIZE,
            GzEncoder::new(File::create(&path).unwrap(), Compression::new(1)),
        ));
        sink.finish().unwrap();

        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[..3], &[0x1f, 0x8b, 0x08], "gzip magic + deflate");
        // 8-byte CRC32 + ISIZE trailer, both zero for an empty member.
        assert_eq!(&bytes[bytes.len() - 8..], &[0u8; 8]);
    }

    #[test]
    fn test_append_to_id_no_tag_tail_appends_to_end() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        rec.append_to_id(":clip5:ACGT");
        assert_eq!(rec.id, "@read1:clip5:ACGT");
    }

    #[test]
    fn test_append_to_id_with_tag_tail_splices_before_first_tab() {
        // Code-review C1 regression guard: when the id carries a tag tail
        // (BamReader-fed uBAM input), the --rename suffix MUST splice into
        // the name portion, not append after the tags — otherwise the last
        // tag's value gets silently corrupted.
        let mut rec = FastqRecord {
            id: "@read1\tCB:Z:ATCGATCG-1\tUB:Z:GCTAGCTA".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        rec.append_to_id(":clip5:AATT");
        assert_eq!(rec.id, "@read1:clip5:AATT\tCB:Z:ATCGATCG-1\tUB:Z:GCTAGCTA");
    }

    #[test]
    fn test_append_to_id_trims_trailing_whitespace() {
        let mut rec = FastqRecord {
            id: "@read1  \n".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        rec.append_to_id(":clip5:A");
        assert_eq!(rec.id, "@read1:clip5:A");
    }

    #[test]
    fn test_fastq_record_clip_5prime() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACGTACGT".to_string(),
            qual: "IIIIIIII".to_string(),
        };
        let clipped = rec.clip_5prime(3);
        assert_eq!(clipped, Some("ACG".to_string()));
        assert_eq!(rec.seq, "TACGT");
        assert_eq!(rec.qual, "IIIII");
    }

    #[test]
    fn test_fastq_record_clip_3prime() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACGTACGT".to_string(),
            qual: "IIIIIIII".to_string(),
        };
        let clipped = rec.clip_3prime(3);
        assert_eq!(clipped, Some("CGT".to_string()));
        assert_eq!(rec.seq, "ACGTA");
        assert_eq!(rec.qual, "IIIII");
    }

    #[test]
    fn test_fastq_record_clip_too_short() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACG".to_string(),
            qual: "III".to_string(),
        };
        // Clip amount >= seq length: no clipping
        assert_eq!(rec.clip_5prime(3), None);
        assert_eq!(rec.clip_3prime(3), None);
        assert_eq!(rec.seq, "ACG");
    }

    #[test]
    fn test_fastq_record_n_count() {
        let rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACNGTnNAC".to_string(),
            qual: "IIIIIIIII".to_string(),
        };
        assert_eq!(rec.n_count(), 3);
    }

    #[test]
    fn test_fastq_record_trim_ns() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "NNACGTNN".to_string(),
            qual: "!!IIII!!".to_string(),
        };
        rec.trim_ns();
        assert_eq!(rec.seq, "ACGT");
        assert_eq!(rec.qual, "IIII");
    }

    #[test]
    fn test_fastq_record_trim_ns_all_n() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "NNNN".to_string(),
            qual: "!!!!".to_string(),
        };
        rec.trim_ns();
        assert_eq!(rec.seq, "");
        assert_eq!(rec.qual, "");
    }

    #[test]
    fn test_fastq_record_truncate() {
        let mut rec = FastqRecord {
            id: "@read1".to_string(),
            seq: "ACGTACGT".to_string(),
            qual: "IIIIIIII".to_string(),
        };
        rec.truncate(4);
        assert_eq!(rec.seq, "ACGT");
        assert_eq!(rec.qual, "IIII");
    }

    #[test]
    fn test_round_trip() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_fastq_round_trip");
        std::fs::create_dir_all(&dir)?;
        let out_path = dir.join("test_round_trip.fq");

        // Write
        let records = vec![
            FastqRecord {
                id: "@read1 description".to_string(),
                seq: "ACGTACGT".to_string(),
                qual: "IIIIIIII".to_string(),
            },
            FastqRecord {
                id: "@read2".to_string(),
                seq: "TGCA".to_string(),
                qual: "!!!!".to_string(),
            },
        ];

        {
            let mut writer = FastqWriter::create(&out_path, false, 1, DEFAULT_GZIP_LEVEL)?;
            for rec in &records {
                writer.write_record(rec)?;
            }
            writer.finish()?;
        }

        // Read back
        let mut reader = FastqReader::open(&out_path)?;
        let r1 = reader.next_record()?.expect("should have record 1");
        assert_eq!(r1.id, "@read1 description");
        assert_eq!(r1.seq, "ACGTACGT");
        assert_eq!(r1.qual, "IIIIIIII");

        let r2 = reader.next_record()?.expect("should have record 2");
        assert_eq!(r2.id, "@read2");
        assert_eq!(r2.seq, "TGCA");
        assert_eq!(r2.qual, "!!!!");

        assert!(reader.next_record()?.is_none());

        // Cleanup
        std::fs::remove_file(&out_path)?;

        Ok(())
    }

    /// §5.1 regression: the parallel writer (`--cores N`) emits each
    /// worker's chunk as its own independently-compressed gzip member,
    /// then concatenates them (RFC 1952 — multi-member gzip is a valid
    /// `.gz` file that decompresses to the concatenation of each
    /// member's payload). The reader has to use `MultiGzDecoder`, NOT
    /// the single-member `GzDecoder`, or it silently truncates output
    /// to whatever fits in the first gzip member when that file is fed
    /// back through Trim Galore as a follow-on input. Originally fixed
    /// in commit 9dcf519 (pre-beta.1) but never had a unit-level
    /// regression test. Locks down the contract: a manually-crafted
    /// 2-member gzip file round-trips through `FastqReader` with all
    /// records from both members.
    #[test]
    fn test_multi_member_gzip_round_trip() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_multi_member_gzip");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("two_member.fq.gz");

        // Member 1: 2 records.
        let mut member1: Vec<u8> = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut member1, Compression::default());
            enc.write_all(b"@m1_r1\nACGT\n+\nIIII\n")?;
            enc.write_all(b"@m1_r2\nGGCC\n+\n!!!!\n")?;
            enc.finish()?;
        }
        // Member 2: 2 more records, separately-compressed.
        let mut member2: Vec<u8> = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut member2, Compression::default());
            enc.write_all(b"@m2_r1\nTTAA\n+\nJJJJ\n")?;
            enc.write_all(b"@m2_r2\nNNNN\n+\n????\n")?;
            enc.finish()?;
        }
        // Concatenate the two gzip members and verify our crafted file
        // is genuinely multi-member (would-be a parsing-error point for
        // single-member decoders).
        std::fs::write(&path, [member1, member2].concat())?;

        let mut reader = FastqReader::open(&path)?;
        let mut ids: Vec<String> = Vec::new();
        while let Some(rec) = reader.next_record()? {
            ids.push(rec.id);
        }

        assert_eq!(
            ids,
            vec!["@m1_r1", "@m1_r2", "@m2_r1", "@m2_r2"],
            "MultiGzDecoder must yield records from both gzip members"
        );

        std::fs::remove_file(&path)?;
        Ok(())
    }

    // ── read_id_prefix (plan v2 Step 4) ──────────────────────────────────

    #[test]
    fn test_read_id_prefix_bare() {
        assert_eq!(read_id_prefix("@read1"), "read1");
    }

    #[test]
    fn test_read_id_prefix_with_description() {
        // Modern Illumina (post-CASAVA-1.8): "@HEADER 1:N:0:CGATCG"
        assert_eq!(read_id_prefix("@read1 1:N:0:CGATCG"), "read1");
    }

    #[test]
    fn test_read_id_prefix_strips_slash_one() {
        // Legacy SRA/ENA: "@read/1" for R1
        assert_eq!(read_id_prefix("@read1/1"), "read1");
    }

    #[test]
    fn test_read_id_prefix_strips_slash_two() {
        assert_eq!(read_id_prefix("@read1/2"), "read1");
    }

    #[test]
    fn test_read_id_prefix_strips_slash_three() {
        // 10X Multiome / scATAC: I1 read tagged "@read/3" in some pipelines.
        assert_eq!(read_id_prefix("@read1/3"), "read1");
    }

    #[test]
    fn test_read_id_prefix_strips_slash_one_with_description() {
        // Combination: legacy slash + modern description.
        assert_eq!(read_id_prefix("@read1/1 1:N:0:CGATCG"), "read1");
    }

    #[test]
    fn test_read_id_prefix_no_strip_slash_four() {
        // Only /1, /2, /3 are stripped. /4 (and any other suffix) is preserved.
        assert_eq!(read_id_prefix("@read1/4"), "read1/4");
    }

    #[test]
    fn test_read_id_prefix_empty() {
        assert_eq!(read_id_prefix(""), "");
    }

    #[test]
    fn test_read_id_prefix_at_only() {
        // Pathological but well-defined: bare '@' means an empty ID.
        assert_eq!(read_id_prefix("@"), "");
    }

    #[test]
    fn test_read_id_prefix_no_at_sign() {
        // Some callers may pass an already-stripped ID — should still work.
        assert_eq!(read_id_prefix("read1"), "read1");
        assert_eq!(read_id_prefix("read1/1"), "read1");
    }

    #[test]
    fn test_read_id_prefix_paired_three_way_sync() {
        // The load-bearing use case: R1, R2, passthrough share a prefix.
        // This is what the sync check across the three FASTQ streams hashes.
        assert_eq!(
            read_id_prefix("@SRR12345.1 1:N:0:ATCG"),
            read_id_prefix("@SRR12345.1 2:N:0:ATCG"),
        );
        assert_eq!(
            read_id_prefix("@SRR12345.1 1:N:0:ATCG"),
            read_id_prefix("@SRR12345.1 3:N:0:ATCG"),
        );
        // Legacy form, same template across three streams.
        assert_eq!(
            read_id_prefix("@SRR12345.1/1"),
            read_id_prefix("@SRR12345.1/2"),
        );
        assert_eq!(
            read_id_prefix("@SRR12345.1/1"),
            read_id_prefix("@SRR12345.1/3"),
        );
    }

    // ---- caller-supplied gzip verdict ----

    fn gz_tmpdir(slug: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_gz(path: &Path, body: &str) -> Result<()> {
        let mut enc = GzEncoder::new(File::create(path)?, Compression::default());
        write!(enc, "{body}")?;
        enc.finish()?;
        Ok(())
    }

    /// REGRESSION. `.fq.bgz` is gzip content under a non-`.gz` extension.
    /// `format::detect_input_format` classifies it `FastqGz` from the
    /// decompressed payload, but `FastqReader::open` derived the verdict from
    /// the filename and read the compressed bytes as text, failing with
    /// "stream did not contain valid UTF-8". `open_with` takes the caller's
    /// already-correct answer instead.
    #[test]
    fn open_with_reads_gzip_under_non_gz_extension() -> Result<()> {
        let dir = gz_tmpdir("tg_fastq_bgz");
        let p = dir.join("s.fq.bgz");
        write_gz(&p, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")?;

        let mut reader = FastqReader::open_with(&p, true)?;
        let first = reader.next_record()?.expect("first record");
        assert_eq!(first.id, "@r1");
        assert_eq!(first.seq, "ACGT");
        let second = reader.next_record()?.expect("second record");
        assert_eq!(second.id, "@r2");
        assert!(reader.next_record()?.is_none());
        Ok(())
    }

    /// The same file through the threaded constructor: the verdict has to
    /// cross the thread boundary, which is a separate code path.
    #[test]
    fn open_threaded_with_reads_gzip_under_non_gz_extension() -> Result<()> {
        let dir = gz_tmpdir("tg_fastq_bgz_threaded");
        let p = dir.join("s.fq.bgz");
        write_gz(&p, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")?;

        let mut reader = FastqReader::open_threaded_with(&p, true)?;
        assert_eq!(reader.next_record()?.expect("first record").id, "@r1");
        assert_eq!(reader.next_record()?.expect("second record").id, "@r2");
        assert!(reader.next_record()?.is_none());
        Ok(())
    }

    /// `sanity_check` runs before everything else in `main()`, so it is the
    /// first place the wrong verdict surfaced.
    #[test]
    fn sanity_check_with_accepts_gzip_under_non_gz_extension() -> Result<()> {
        let dir = gz_tmpdir("tg_fastq_bgz_sanity");
        let p = dir.join("s.fq.bgz");
        write_gz(&p, "@r1\nACGT\n+\nIIII\n")?;

        FastqReader::sanity_check_with(&p, true)?;
        // The un-suffixed entry point must agree: it sniffs the content too,
        // so there is no longer a "wrong" default to fall into.
        FastqReader::sanity_check(&p)?;
        Ok(())
    }

    /// A `.gz`-named plain file: the caller's `false` must win over the
    /// filename, the mirror image of the `.bgz` case.
    ///
    /// This direction is a user-visible behaviour change, not only an internal
    /// one: before, a plain FASTQ misnamed `.fastq.gz` failed with
    /// `invalid gzip header`, and now it reads. See the CHANGELOG entry.
    #[test]
    fn open_with_honours_a_plain_verdict_on_a_gz_name() -> Result<()> {
        let dir = gz_tmpdir("tg_fastq_plain_gz_name");
        let p = dir.join("plain.fq.gz");
        std::fs::write(&p, b"@r1\nACGT\n+\nIIII\n")?;

        let mut reader = FastqReader::open_with(&p, false)?;
        assert_eq!(reader.next_record()?.expect("record").id, "@r1");
        assert!(reader.next_record()?.is_none());
        Ok(())
    }
}
