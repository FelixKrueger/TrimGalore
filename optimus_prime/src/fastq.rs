//! FASTQ record type and streaming I/O.
//!
//! Provides gzip-aware reading and writing of FASTQ records with 64KB buffers
//! for efficient I/O throughput.

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Buffer size for gzip I/O — 64KB for throughput (flate2 default of 8KB is too small).
const BUF_SIZE: usize = 64 * 1024;

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
        writeln!(writer, "{}", self.id)?;
        writeln!(writer, "{}", self.seq)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", self.qual)?;
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

/// A streaming FASTQ reader that handles both plain and gzipped files.
pub struct FastqReader {
    reader: Box<dyn BufRead>,
    line_buf: String,
}

impl FastqReader {
    /// Open a FASTQ file, auto-detecting gzip from the `.gz` extension.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open input file: {}", path.display()))?;

        let reader: Box<dyn BufRead> = if path
            .extension()
            .is_some_and(|ext| ext == "gz")
        {
            Box::new(BufReader::with_capacity(BUF_SIZE, GzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(BUF_SIZE, file))
        };

        Ok(FastqReader {
            reader,
            line_buf: String::with_capacity(512),
        })
    }

    /// Read the next FASTQ record. Returns None at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        // Line 1: ID (starts with @)
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            return Ok(None); // EOF
        }
        let id = self.line_buf.trim_end_matches(['\n', '\r']).to_string();

        // Line 2: Sequence
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Truncated FASTQ: missing sequence line after {}", id);
        }
        let seq = self.line_buf.trim_end_matches(['\n', '\r']).to_string();

        // Line 3: Plus line (discard)
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Truncated FASTQ: missing '+' line after {}", id);
        }

        // Line 4: Quality
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Truncated FASTQ: missing quality line after {}", id);
        }
        let qual = self.line_buf.trim_end_matches(['\n', '\r']).to_string();

        Ok(Some(FastqRecord { id, seq, qual }))
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
                         by Optimus Prime (sequence is: '{}'). Colorspace data requires \
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
    writer: Box<dyn Write>,
}

impl FastqWriter {
    /// Create a new FASTQ writer. Gzip-compresses if `gzip` is true.
    pub fn create<P: AsRef<Path>>(path: P, gzip: bool) -> Result<Self> {
        let path = path.as_ref();

        // Ensure parent directory exists
        if let Some(parent) = path.parent() {
            if !parent.exists() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("Failed to create output directory: {}", parent.display()))?;
            }
        }

        let file = File::create(path)
            .with_context(|| format!("Failed to create output file: {}", path.display()))?;

        let writer: Box<dyn Write> = if gzip {
            // Compression level 6 matches system gzip default
            Box::new(BufWriter::with_capacity(
                BUF_SIZE,
                GzEncoder::new(file, Compression::new(6)),
            ))
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
        let dir = std::env::temp_dir().join("optimus_prime_test");
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
            let mut writer = FastqWriter::create(&out_path, false)?;
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
}
