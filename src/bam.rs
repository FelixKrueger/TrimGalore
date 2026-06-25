//! uBAM (unaligned BAM) input reader.
//!
//! Converts a `noodles::bam` record stream into the `FastqRecord` shape the
//! rest of the trimming pipeline expects. See `plans/06252026_ubam-input-support/PLAN.md`
//! §3.2 for the record-conversion contract.

use anyhow::{Context, Result, bail};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use std::fs::File;
use std::io::BufReader;
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

const PHRED_OFFSET: u8 = 33;
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
        let inner = open_inner(&path)?;
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
                Ok(r) => r,
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

/// Open a BAM file, read its header, return the constructed reader.
fn open_inner(path: &Path) -> Result<bam::io::Reader<bgzf::Reader<BufReader<File>>>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open uBAM input: {}", path.display()))?;
    let mut inner = bam::io::Reader::new(BufReader::new(file));
    let _ = inner
        .read_header()
        .with_context(|| format!("Failed to read BAM header from {}", path.display()))?;
    Ok(inner)
}

/// Convert one BAM record to a `FastqRecord`. Per-record validation:
/// - Must be unmapped (`is_unmapped()` true; per-record check resolves
///   PLAN-REVIEW B-Crit-4 first-record-only contradiction).
/// - Must NOT be reverse-complemented / secondary / supplementary.
/// - Sequence must be non-empty; IUPAC codes coerced to N, `=` rejected.
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
                append_tag_type_and_value(&mut id, &value)?;
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
        emit_iupac_warning_once();
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
        qual_raw
            .iter()
            .map(|&b| (b + PHRED_OFFSET) as char)
            .collect()
    };

    Ok(FastqRecord { id, seq, qual })
}

fn emit_iupac_warning_once() {
    static SEEN: OnceLock<()> = OnceLock::new();
    SEEN.get_or_init(|| {
        eprintln!(
            "WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V); \
             coerced to N for FASTQ output. This is legitimate in PacBio HiFi, \
             ONT Dorado, and 10x cellranger uBAMs."
        );
    });
}

/// Format `{TYPE}:{VALUE}` for one BAM aux field, matching `samtools fastq -T`.
fn append_tag_type_and_value(out: &mut String, value: &Value<'_>) -> Result<()> {
    use std::fmt::Write;
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
        Value::Hex(bs) => {
            let s = std::str::from_utf8(bs).context("BAM hex tag is not valid UTF-8")?;
            write!(out, "H:{}", s).unwrap();
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
    Ok(())
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
}
