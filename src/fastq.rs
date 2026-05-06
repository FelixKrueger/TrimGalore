//! FASTQ record type and streaming I/O.
//!
//! Provides gzip-aware reading and writing of FASTQ records with 64KB buffers
//! for efficient I/O throughput.

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use gzp::deflate::Gzip;
use gzp::par::compress::ParCompressBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Buffer size for gzip I/O — 64KB for throughput (flate2 default of 8KB is too small).
const BUF_SIZE: usize = 64 * 1024;

/// Default output gzip compression level. Set to 1 (fastest) — at
/// Buckberry-scale (84M reads, 38% adapter rate, cores=8) the
/// compression CPU dominated wall time on saturated workers; lowering
/// from level 6 to 1 measured ~−23% wall and ~−43% user-CPU at
/// byte-identity of the decompressed output (gzip framing differs
/// but `gzip -dc` yields the same bytes). Trade: output `.fq.gz`
/// files are ~75% larger. See #248 (item #1 in @an-altosian's perf
/// audit). Users who care about storage more than runtime can opt
/// into `--high_compression` to get level 6 instead.
pub const OUTPUT_GZIP_LEVEL: u32 = 1;

/// Gzip compression level used when `--high_compression` is set.
/// Matches Rust's `flate2::Compression::default()` (level 6) and
/// gzip(1)'s default — the "balanced" point on the size/speed curve.
/// Output bytes are ~75% smaller than level 1; per-block compression
/// CPU is roughly 2× higher.
pub const HIGH_COMPRESSION_GZIP_LEVEL: u32 = 9;

/// Returns the gzip output level to use given the `--high_compression`
/// flag state. Centralises the size-vs-speed trade so every output
/// callsite reads the same way: `Compression::new(output_gzip_level(...))`.
pub fn output_gzip_level(high_compression: bool) -> u32 {
    if high_compression {
        HIGH_COMPRESSION_GZIP_LEVEL
    } else {
        OUTPUT_GZIP_LEVEL
    }
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
    pub fn append_to_id(&mut self, suffix: &str) {
        // Strip any trailing whitespace/newline from id before appending
        let id = self.id.trim_end().to_string();
        self.id = format!("{}{}", id, suffix);
    }
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
        reader: Box<dyn BufRead>,
        line_buf: String,
    },
    Threaded {
        rx: std::sync::mpsc::Receiver<Result<Vec<FastqRecord>>>,
        buffer: Vec<FastqRecord>,
        buf_pos: usize,
        _handle: std::thread::JoinHandle<()>,
    },
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
    /// Open a FASTQ file for synchronous reading.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open input file: {}", path.display()))?;

        let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(BufReader::with_capacity(
                BUF_SIZE,
                MultiGzDecoder::new(file),
            ))
        } else {
            Box::new(BufReader::with_capacity(BUF_SIZE, file))
        };

        Ok(FastqReader {
            source: ReaderSource::Direct {
                reader,
                line_buf: String::with_capacity(512),
            },
        })
    }

    /// Open a FASTQ file with background decompression on a dedicated thread.
    ///
    /// Returns immediately. A background thread reads and decompresses
    /// the file, sending records through a bounded channel. The main
    /// thread calls `next_record()` which receives from the channel,
    /// overlapping decompression with processing.
    pub fn open_threaded<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();

        // Validate the file exists before spawning the thread
        if !path.exists() {
            bail!("Input file not found: {}", path.display());
        }

        let (tx, rx) = std::sync::mpsc::sync_channel(READER_CHANNEL_BATCHES);

        let handle = std::thread::spawn(move || {
            let mut reader = match Self::open_direct(&path) {
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

    /// Internal: open a file and return the raw reader components (for use in threads).
    fn open_direct(path: &Path) -> Result<(Box<dyn BufRead + Send>, String)> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open input file: {}", path.display()))?;

        let reader: Box<dyn BufRead + Send> = if path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(BufReader::with_capacity(
                BUF_SIZE,
                MultiGzDecoder::new(file),
            ))
        } else {
            Box::new(BufReader::with_capacity(BUF_SIZE, file))
        };

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

    /// Perform input sanity checks on the first record.
    /// Checks: FASTQ format validation, colorspace detection, empty file.
    pub fn sanity_check<P: AsRef<Path>>(path: P) -> Result<()> {
        let path = path.as_ref();
        let mut reader = FastqReader::open(path)?;

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

/// A streaming FASTQ writer that handles both plain and gzipped output.
pub struct FastqWriter {
    writer: Box<dyn Write + Send>,
}

impl FastqWriter {
    /// Create a new FASTQ writer. Gzip-compresses if `gzip` is true.
    /// When `cores` > 1 and gzip is enabled, uses parallel gzip compression.
    /// `gzip_level` is the deflate level (1 = fastest/largest, 6 = balanced,
    /// 9 = smallest/slowest). Use `output_gzip_level(high_compression)`
    /// to derive the right value from a `--high_compression` flag.
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

        let file = File::create(path)
            .with_context(|| format!("Failed to create output file: {}", path.display()))?;

        let writer: Box<dyn Write + Send> = if gzip {
            if cores > 1 {
                // Parallel gzip: split output into independently-compressed blocks
                Box::new(
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
                Box::new(BufWriter::with_capacity(
                    BUF_SIZE,
                    GzEncoder::new(file, Compression::new(gzip_level)),
                ))
            }
        } else {
            Box::new(BufWriter::with_capacity(BUF_SIZE, file))
        };

        Ok(FastqWriter { writer })
    }

    /// Write a FASTQ record.
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        record.write_to(&mut self.writer)
    }

    /// Flush and finalize the writer.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

impl Drop for FastqWriter {
    fn drop(&mut self) {
        let _ = self.writer.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
            let mut writer = FastqWriter::create(&out_path, false, 1, OUTPUT_GZIP_LEVEL)?;
            for rec in &records {
                writer.write_record(rec)?;
            }
            writer.flush()?;
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
}
