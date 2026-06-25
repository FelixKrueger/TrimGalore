//! Input-format detection for FASTQ vs uBAM.
//!
//! Two-stage check (PLAN §3.1 / §5 step 2):
//!
//! 1. **Cheap byte peek.** First byte `@` → plain FASTQ. First three bytes
//!    `1F 8B 08` → gzip family; fall through.
//! 2. **Decompress first block + payload check.** For the gzip family,
//!    decompress and read 4 bytes. If those equal `BAM\1` → unaligned BAM,
//!    otherwise gzipped FASTQ.
//!
//! The decompress step is **load-bearing**. Heuristic-only detection (matching
//! on the BGZF `BC` extra-field subfield) would misclassify `bgzip x.fq` —
//! BGZF-framed FASTQ — as BAM, because the framing is identical. The only
//! safe discriminator is the decompressed payload. Both plan reviewers
//! caught this (A-C2 + B-Crit-2).

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Classification of an input file's container format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    /// Plain (uncompressed) FASTQ.
    FastqPlain,
    /// Gzip-compressed FASTQ — includes both `gzip x.fq` (plain multi-member
    /// gzip) and `bgzip x.fq` (BGZF-framed gzip). Both decompress identically
    /// via `MultiGzDecoder` downstream.
    FastqGz,
    /// Unaligned BAM — BGZF-framed, first decompressed bytes are `BAM\1`.
    UnalignedBam,
}

/// Peek the first bytes of `path` and classify by content (NOT by filename).
///
/// See module-level docs for the algorithm and the rationale for the
/// decompress-and-check step.
pub fn detect_input_format(path: &Path) -> Result<InputFormat> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open input file: {}", path.display()))?;
    let mut peek = [0u8; 4];
    let n = file
        .read(&mut peek)
        .with_context(|| format!("Failed to read from {}", path.display()))?;

    if n == 0 {
        bail!("Input file '{}' is empty", path.display());
    }

    if peek[0] == b'@' {
        return Ok(InputFormat::FastqPlain);
    }

    // Gzip family — magic bytes `1F 8B 08`. The third byte is the
    // compression method; flate is what every modern gzip implementation
    // uses, and the BGZF spec also fixes it at 08.
    if n >= 3 && peek[0] == 0x1F && peek[1] == 0x8B && peek[2] == 0x08 {
        // Re-open from byte 0 (cheaper than seeking and sidesteps any
        // interaction between `MultiGzDecoder` and the already-consumed
        // prefix; the alternative `PeekReader` wrapper in PLAN §5 step 2.3
        // is a documented future optimization, not required for correctness).
        let file = File::open(path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut payload = [0u8; 4];
        let nread = decoder.read(&mut payload).with_context(|| {
            format!(
                "Failed to decompress first block of '{}' for format detection",
                path.display()
            )
        })?;
        if nread == 4 && &payload == b"BAM\x01" {
            return Ok(InputFormat::UnalignedBam);
        }
        return Ok(InputFormat::FastqGz);
    }

    bail!(
        "Input '{}' is not recognised as FASTQ (plain or gzipped) or unaligned BAM",
        path.display()
    );
}

/// Open a sync (single-threaded) reader for `path`, dispatching by detected
/// format. Returns the reader as a boxed trait object so callers can stay
/// reader-source-agnostic. `preserve_tags` is honoured for BAM input;
/// silently ignored for FASTQ.
pub fn open_sync_reader(
    path: &Path,
    preserve_tags: &[String],
) -> Result<Box<dyn crate::fastq::RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => {
            Ok(Box::new(crate::fastq::FastqReader::open(path)?))
        }
        InputFormat::UnalignedBam => Ok(Box::new(
            crate::bam::BamReader::open(path)?.with_preserved_tags(preserve_tags),
        )),
    }
}

/// Open a threaded (background-decompression) reader for `path`, dispatching
/// by detected format. `preserve_tags` is honoured for BAM input; silently
/// ignored for FASTQ.
pub fn open_threaded_reader(
    path: &Path,
    preserve_tags: &[String],
) -> Result<Box<dyn crate::fastq::RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => {
            Ok(Box::new(crate::fastq::FastqReader::open_threaded(path)?))
        }
        InputFormat::UnalignedBam => Ok(Box::new(crate::bam::BamReader::open_threaded_with_tags(
            path,
            preserve_tags,
        )?)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn fresh_tmpdir(slug: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn detect_fastq_plain_from_at_sign() -> Result<()> {
        let dir = fresh_tmpdir("tg_format_plain");
        let p = dir.join("a.fq");
        std::fs::write(&p, b"@read1\nACGT\n+\nIIII\n")?;
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqPlain);
        Ok(())
    }

    #[test]
    fn detect_fastq_gz_from_plain_gzip() -> Result<()> {
        let dir = fresh_tmpdir("tg_format_gz");
        let p = dir.join("a.fq.gz");
        let mut e = GzEncoder::new(std::fs::File::create(&p)?, Compression::default());
        e.write_all(b"@read1\nACGT\n+\nIIII\n")?;
        e.finish()?;
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqGz);
        Ok(())
    }

    #[test]
    fn detect_unaligned_bam_via_decompressed_magic() -> Result<()> {
        let p = std::path::PathBuf::from("test_files/ubam_test.bam");
        assert_eq!(detect_input_format(&p)?, InputFormat::UnalignedBam);
        Ok(())
    }

    #[test]
    fn detect_paired_ubam_via_decompressed_magic() -> Result<()> {
        let p = std::path::PathBuf::from("test_files/ubam_paired_test.bam");
        assert_eq!(detect_input_format(&p)?, InputFormat::UnalignedBam);
        Ok(())
    }

    /// LOAD-BEARING (PLAN_REVIEW_A C2 + PLAN_REVIEW_B Crit-2). `bgzip x.fq`
    /// produces BGZF-framed FASTQ — the same framing as BAM. The only safe
    /// discriminator is the decompressed payload (`@` for FASTQ vs `BAM\1`
    /// for BAM). This test uses noodles' BGZF writer to produce a genuine
    /// BGZF-framed file with FASTQ contents and confirms it classifies
    /// as `FastqGz`, NOT `UnalignedBam`.
    #[test]
    fn detect_bgzipped_fastq_is_fastq_not_bam() -> Result<()> {
        use noodles::bgzf;
        let dir = fresh_tmpdir("tg_format_bgzf_fq");
        let p = dir.join("a.fq.bgz");
        {
            let mut w = bgzf::Writer::new(std::fs::File::create(&p)?);
            w.write_all(b"@read1\nACGT\n+\nIIII\n")?;
            w.finish()?;
        }
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqGz);
        Ok(())
    }

    #[test]
    fn detect_rejects_random_binary() {
        let dir = fresh_tmpdir("tg_format_neg");
        let p = dir.join("noise.bin");
        std::fs::write(&p, [0u8, 1, 2, 3, 4, 5, 6, 7]).unwrap();
        let r = detect_input_format(&p);
        assert!(r.is_err(), "non-FASTQ non-BAM input must error");
    }

    #[test]
    fn detect_empty_file_errors() {
        let dir = fresh_tmpdir("tg_format_empty");
        let p = dir.join("empty.fq");
        std::fs::write(&p, b"").unwrap();
        let r = detect_input_format(&p);
        assert!(r.is_err());
    }
}
