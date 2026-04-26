//! Specialty trimming modes that bypass the normal trimming pipeline.
//!
//! These modes (--hardtrim5, --hardtrim3, --clock, --implicon) process
//! input files with simple fixed operations and exit immediately.

use anyhow::{Result, bail};
use std::path::{Path, PathBuf};

use crate::fastq::{FastqReader, FastqWriter};
use crate::io as naming;

/// Hard-trim every read to keep only the first `keep` bases from the 5' end.
///
/// Output filename: `*.{keep}bp_5prime.fq(.gz)`
pub fn hardtrim5(
    input: &Path,
    keep: usize,
    gzip: bool,
    output_dir: Option<&Path>,
    rename: bool,
    cores: usize,
) -> Result<()> {
    let output_path = hardtrim_output_name(input, keep, "5prime", output_dir, gzip);
    eprintln!(
        "Writing hard-trimmed (first {}bp) version of '{}' to '{}'",
        keep,
        input.display(),
        output_path.display()
    );

    let mut reader = FastqReader::open(input)?;
    let mut writer = FastqWriter::create(&output_path, gzip, cores)?;
    let mut count: usize = 0;

    while let Some(mut record) = reader.next_record()? {
        count += 1;
        if record.seq.len() > keep {
            if rename {
                let clipped = &record.seq[keep..];
                record.append_to_id(&format!(":clip5:{}", clipped));
            }
            record.truncate(keep);
        }
        writer.write_record(&record)?;
    }

    writer.flush()?;
    drop(writer);

    eprintln!("Finished writing {} sequences\n", count);
    Ok(())
}

/// Hard-trim every read to keep only the last `keep` bases from the 3' end.
///
/// Output filename: `*.{keep}bp_3prime.fq(.gz)`
pub fn hardtrim3(
    input: &Path,
    keep: usize,
    gzip: bool,
    output_dir: Option<&Path>,
    rename: bool,
    cores: usize,
) -> Result<()> {
    let output_path = hardtrim_output_name(input, keep, "3prime", output_dir, gzip);
    eprintln!(
        "Writing hard-trimmed (last {}bp) version of '{}' to '{}'",
        keep,
        input.display(),
        output_path.display()
    );

    let mut reader = FastqReader::open(input)?;
    let mut writer = FastqWriter::create(&output_path, gzip, cores)?;
    let mut count: usize = 0;

    while let Some(mut record) = reader.next_record()? {
        count += 1;
        if record.seq.len() > keep {
            let start = record.seq.len() - keep;
            if rename {
                let clipped = &record.seq[..start];
                record.append_to_id(&format!(":clip3:{}", clipped));
            }
            record.clip_5prime(start);
        }
        writer.write_record(&record)?;
    }

    writer.flush()?;
    drop(writer);

    eprintln!("Finished writing {} sequences\n", count);
    Ok(())
}

/// Epigenetic Clock preprocessing (paired-end only).
///
/// Read layout (both R1 and R2):
/// ```text
/// Position: 0        8    12 13                (R1)
///           UUUUUUUU CAGT A  FFFFFFF...
/// Position: 0        8    12   15              (R2)
///           UUUUUUUU CAGT A  BB FFFFF...       (BB = end-repair bias)
/// ```
///
/// Extracts 8bp UMI + 4bp fixed sequence from both reads, appends to read IDs
/// as `:R1:<umi1>:R2:<umi2>:F1:<fix1>:F2:<fix2>`, then clips:
/// - R1 at position 13 (removes UMI + CAGT + A)
/// - R2 at position 15 (removes UMI + CAGT + A + 2bp end-repair bias)
pub fn clock(
    input_r1: &Path,
    input_r2: &Path,
    gzip: bool,
    output_dir: Option<&Path>,
    cores: usize,
) -> Result<()> {
    let out1 = clock_output_name(input_r1, "R1", output_dir, gzip);
    let out2 = clock_output_name(input_r2, "R2", output_dir, gzip);

    eprintln!(
        "Writing Clock-processed version of '{}' to '{}'",
        input_r1.display(),
        out1.display()
    );
    eprintln!(
        "Writing Clock-processed version of '{}' to '{}'",
        input_r2.display(),
        out2.display()
    );

    let mut reader_r1 = FastqReader::open(input_r1)?;
    let mut reader_r2 = FastqReader::open(input_r2)?;
    let mut writer_r1 = FastqWriter::create(&out1, gzip, cores)?;
    let mut writer_r2 = FastqWriter::create(&out2, gzip, cores)?;

    let mut count: usize = 0;
    let mut filtered_count: usize = 0;

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;

        match (rec1, rec2) {
            (Some(mut r1), Some(mut r2)) => {
                count += 1;
                if count % 1_000_000 == 0 {
                    eprintln!("Processed {} sequences so far...", count);
                }

                // Extract UMI (8bp) and fixed sequence (4bp) from both reads
                let r1_seq = r1.seq.as_bytes();
                let r2_seq = r2.seq.as_bytes();

                if r1_seq.len() < 13 || r2_seq.len() < 15 {
                    bail!(
                        "Read too short for Clock processing at read {}: R1={}bp, R2={}bp (need >=13, >=15)",
                        count,
                        r1_seq.len(),
                        r2_seq.len()
                    );
                }

                let r1_umi = &r1.seq[..8];
                let r2_umi = &r2.seq[..8];
                let r1_fix = &r1.seq[8..12];
                let r2_fix = &r2.seq[8..12];

                // Count reads where fixed sequence is not CAGT
                if r1_fix != "CAGT" || r2_fix != "CAGT" {
                    filtered_count += 1;
                }

                // Append UMI info to both read IDs
                let suffix = format!(":R1:{}:R2:{}:F1:{}:F2:{}", r1_umi, r2_umi, r1_fix, r2_fix);
                r1.append_to_id(&suffix);
                r2.append_to_id(&suffix);

                // Clip: R1 at pos 13, R2 at pos 15
                let r1_trimmed_seq = r1.seq[13..].to_string();
                let r1_trimmed_qual = r1.qual[13..].to_string();
                r1.seq = r1_trimmed_seq;
                r1.qual = r1_trimmed_qual;

                let r2_trimmed_seq = r2.seq[15..].to_string();
                let r2_trimmed_qual = r2.qual[15..].to_string();
                r2.seq = r2_trimmed_seq;
                r2.qual = r2_trimmed_qual;

                writer_r1.write_record(&r1)?;
                writer_r2.write_record(&r2)?;
            }
            (None, None) => break,
            _ => bail!("Paired-end files have different numbers of reads!"),
        }
    }

    writer_r1.flush()?;
    writer_r2.flush()?;
    drop(writer_r1);
    drop(writer_r2);

    let perc = if count > 0 {
        format!("{:.2}", filtered_count as f64 / count as f64 * 100.0)
    } else {
        "N/A".to_string()
    };
    eprintln!("\nSequences processed in total: {}", count);
    eprintln!(
        "thereof had fixed sequence CAGT in both R1 and R2: {} ({}%)\n",
        filtered_count, perc
    );

    Ok(())
}

/// IMPLICON UMI transfer (paired-end only).
///
/// Transfers the first `umi_length` bases from R2 to both read IDs as a
/// UMI barcode (`:BARCODE`), then clips R2 by `umi_length` bases.
/// R1 sequence is unchanged (only the read ID gets the barcode appended).
pub fn implicon(
    input_r1: &Path,
    input_r2: &Path,
    umi_length: usize,
    gzip: bool,
    output_dir: Option<&Path>,
    cores: usize,
) -> Result<()> {
    let out1 = implicon_output_name(input_r1, umi_length, "R1", output_dir, gzip);
    let out2 = implicon_output_name(input_r2, umi_length, "R2", output_dir, gzip);

    eprintln!(
        "Writing UMI-trimmed version of '{}' to '{}'",
        input_r1.display(),
        out1.display()
    );
    eprintln!(
        "Writing UMI-trimmed version of '{}' to '{}'",
        input_r2.display(),
        out2.display()
    );

    let mut reader_r1 = FastqReader::open(input_r1)?;
    let mut reader_r2 = FastqReader::open(input_r2)?;
    let mut writer_r1 = FastqWriter::create(&out1, gzip, cores)?;
    let mut writer_r2 = FastqWriter::create(&out2, gzip, cores)?;

    let mut count: usize = 0;

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;

        match (rec1, rec2) {
            (Some(mut r1), Some(mut r2)) => {
                count += 1;
                if count % 1_000_000 == 0 {
                    eprintln!("Processed {} sequences so far...", count);
                }

                if r2.seq.len() < umi_length {
                    bail!(
                        "R2 read {} is shorter ({} bp) than UMI length ({})",
                        count,
                        r2.seq.len(),
                        umi_length
                    );
                }

                // Extract UMI from R2 5' end
                let barcode = r2.seq[..umi_length].to_string();

                // Append barcode to both read IDs
                let suffix = format!(":{}", barcode);
                r1.append_to_id(&suffix);
                r2.append_to_id(&suffix);

                // Clip R2 by umi_length bases from 5' end
                let r2_trimmed_seq = r2.seq[umi_length..].to_string();
                let r2_trimmed_qual = r2.qual[umi_length..].to_string();
                r2.seq = r2_trimmed_seq;
                r2.qual = r2_trimmed_qual;

                // Write (R1 sequence unchanged, R2 clipped)
                writer_r1.write_record(&r1)?;
                writer_r2.write_record(&r2)?;
            }
            (None, None) => break,
            _ => bail!("Paired-end files have different numbers of reads!"),
        }
    }

    writer_r1.flush()?;
    writer_r2.flush()?;
    drop(writer_r1);
    drop(writer_r2);

    eprintln!("\nSequences processed in total: {}\n", count);
    Ok(())
}

// --- Output naming helpers ---

fn hardtrim_output_name(
    input: &Path,
    keep: usize,
    end: &str,
    output_dir: Option<&Path>,
    gzip: bool,
) -> PathBuf {
    let stem = naming::strip_fastq_extensions(input);
    let ext = if gzip { ".fq.gz" } else { ".fq" };
    let filename = format!("{}.{}bp_{}{}", stem, keep, end, ext);
    match output_dir {
        Some(dir) => dir.join(filename),
        None => PathBuf::from(filename),
    }
}

pub fn clock_output_name(
    input: &Path,
    read: &str,
    output_dir: Option<&Path>,
    gzip: bool,
) -> PathBuf {
    let stem = naming::strip_fastq_extensions(input);
    let ext = if gzip { ".fq.gz" } else { ".fq" };
    let filename = format!("{}.clock_UMI.{}{}", stem, read, ext);
    match output_dir {
        Some(dir) => dir.join(filename),
        None => PathBuf::from(filename),
    }
}

pub fn implicon_output_name(
    input: &Path,
    umi_len: usize,
    read: &str,
    output_dir: Option<&Path>,
    gzip: bool,
) -> PathBuf {
    let mut stem = naming::strip_fastq_extensions(input);
    // TrimGalore strips _R1 from R1 stem and _R2/_R3/_R4 from R2 stem
    // to avoid double-R in the output filename (e.g., sample_R1_8bp_UMI_R1.fastq)
    if read == "R1" {
        if let Some(s) = stem.strip_suffix("_R1") {
            stem = s.to_string();
        } else if let Some(s) = stem.strip_suffix("R1") {
            stem = s.to_string();
        }
    } else {
        for suffix in &["_R2", "R2", "_R3", "R3", "_R4", "R4"] {
            if let Some(s) = stem.strip_suffix(suffix) {
                stem = s.to_string();
                break;
            }
        }
    }
    // TrimGalore uses .fastq(.gz) for implicon, not .fq(.gz)
    let ext = if gzip { ".fastq.gz" } else { ".fastq" };
    let filename = format!("{}_{}bp_UMI_{}{}", stem, umi_len, read, ext);
    match output_dir {
        Some(dir) => dir.join(filename),
        None => PathBuf::from(filename),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::FastqRecord;

    fn mk_rec(id: &str, seq: &str, qual: &str) -> FastqRecord {
        FastqRecord {
            id: id.to_string(),
            seq: seq.to_string(),
            qual: qual.to_string(),
        }
    }

    fn write_fastq(path: &Path, records: &[FastqRecord]) -> Result<()> {
        let mut w = FastqWriter::create(path, false, 1)?;
        for r in records {
            w.write_record(r)?;
        }
        w.flush()?;
        Ok(())
    }

    fn read_fastq(path: &Path) -> Result<Vec<FastqRecord>> {
        let mut reader = FastqReader::open(path)?;
        let mut out = Vec::new();
        while let Some(r) = reader.next_record()? {
            out.push(r);
        }
        Ok(out)
    }

    // --- --clock ---

    #[test]
    fn test_clock_happy_path_extracts_umi_and_clips() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_clock_happy");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("sample_R1.fq");
        let r2_path = dir.join("sample_R2.fq");

        // R1 layout: 8bp UMI + 4bp "CAGT" + 1bp + 7bp rest = 20bp, clipped at [13..] → 7bp rest
        // R2 layout: 8bp UMI + 4bp "CAGT" + 3bp + 5bp rest = 20bp, clipped at [15..] → 5bp rest
        write_fastq(
            &r1_path,
            &[
                mk_rec("@read1", "AAAAAAAACAGTNCCCCCCC", "IIIIIIIIIIIIIIIIIIII"),
                mk_rec("@read2", "GGGGGGGGCAGTNCCCCCCC", "IIIIIIIIIIIIIIIIIIII"),
            ],
        )?;
        write_fastq(
            &r2_path,
            &[
                mk_rec("@read1", "TTTTTTTTCAGTNNNGGGGG", "IIIIIIIIIIIIIIIIIIII"),
                mk_rec("@read2", "AAAAAAAACAGTNNNGGGGG", "IIIIIIIIIIIIIIIIIIII"),
            ],
        )?;

        clock(&r1_path, &r2_path, false, Some(&dir), 1)?;

        let out_r1 = dir.join("sample_R1.clock_UMI.R1.fq");
        let out_r2 = dir.join("sample_R2.clock_UMI.R2.fq");
        assert!(out_r1.exists(), "R1 output should exist at {:?}", out_r1);
        assert!(out_r2.exists(), "R2 output should exist at {:?}", out_r2);

        let r1_recs = read_fastq(&out_r1)?;
        let r2_recs = read_fastq(&out_r2)?;
        assert_eq!(r1_recs.len(), 2);
        assert_eq!(r2_recs.len(), 2);

        // Read 1: both IDs get the same suffix (from R1 and R2 UMIs/fixes)
        assert_eq!(
            r1_recs[0].id,
            "@read1:R1:AAAAAAAA:R2:TTTTTTTT:F1:CAGT:F2:CAGT"
        );
        assert_eq!(r2_recs[0].id, r1_recs[0].id);
        assert_eq!(r1_recs[0].seq, "CCCCCCC"); // [13..] of AAAAAAAACAGTNCCCCCCC
        assert_eq!(r2_recs[0].seq, "GGGGG"); // [15..] of TTTTTTTTCAGTNNNGGGGG
        // Quality clipped to same length
        assert_eq!(r1_recs[0].qual.len(), r1_recs[0].seq.len());
        assert_eq!(r2_recs[0].qual.len(), r2_recs[0].seq.len());

        // Read 2: different R1 UMI, same R2 UMI but different base
        assert_eq!(
            r1_recs[1].id,
            "@read2:R1:GGGGGGGG:R2:AAAAAAAA:F1:CAGT:F2:CAGT"
        );
        assert_eq!(r1_recs[1].seq, "CCCCCCC");
        assert_eq!(r2_recs[1].seq, "GGGGG");

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_clock_r1_too_short_errors() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_clock_r1_short");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("short_R1.fq");
        let r2_path = dir.join("short_R2.fq");

        // R1 of 12bp — below the 13bp minimum
        write_fastq(&r1_path, &[mk_rec("@r", "AAAAAAAACAGT", "IIIIIIIIIIII")])?;
        write_fastq(
            &r2_path,
            &[mk_rec("@r", "TTTTTTTTCAGTNNNN", "IIIIIIIIIIIIIIII")],
        )?;

        let result = clock(&r1_path, &r2_path, false, Some(&dir), 1);
        assert!(result.is_err(), "should bail on too-short R1");
        let err = format!("{}", result.unwrap_err());
        assert!(err.contains("too short"), "got: {}", err);

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_clock_r2_too_short_errors() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_clock_r2_short");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("short_R1.fq");
        let r2_path = dir.join("short_R2.fq");

        // R1 OK (≥13), R2 only 14bp (below 15bp minimum)
        write_fastq(
            &r1_path,
            &[mk_rec("@r", "AAAAAAAACAGTNC", "IIIIIIIIIIIIII")],
        )?;
        write_fastq(
            &r2_path,
            &[mk_rec("@r", "TTTTTTTTCAGTNN", "IIIIIIIIIIIIII")],
        )?;

        let result = clock(&r1_path, &r2_path, false, Some(&dir), 1);
        assert!(result.is_err(), "should bail on too-short R2");

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    // --- --implicon ---

    #[test]
    fn test_implicon_happy_path_umi_8() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_implicon_happy");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("sample_R1.fq");
        let r2_path = dir.join("sample_R2.fq");

        // R1 is untouched; R2 gets its first 8bp extracted as UMI, then clipped.
        write_fastq(
            &r1_path,
            &[
                mk_rec("@read1", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII"),
                mk_rec("@read2", "TTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIII"),
            ],
        )?;
        write_fastq(
            &r2_path,
            &[
                mk_rec("@read1", "AAAAAAAAGGGGGGGG", "IIIIIIIIIIIIIIII"),
                mk_rec("@read2", "CCCCCCCCGGGGGGGG", "IIIIIIIIIIIIIIII"),
            ],
        )?;

        implicon(&r1_path, &r2_path, 8, false, Some(&dir), 1)?;

        // Filename pattern: {stem_minus_R1_suffix}_{umi}bp_UMI_R[12].fastq
        let out_r1 = dir.join("sample_8bp_UMI_R1.fastq");
        let out_r2 = dir.join("sample_8bp_UMI_R2.fastq");
        assert!(out_r1.exists(), "R1 output should exist at {:?}", out_r1);
        assert!(out_r2.exists(), "R2 output should exist at {:?}", out_r2);

        let r1_recs = read_fastq(&out_r1)?;
        let r2_recs = read_fastq(&out_r2)?;
        assert_eq!(r1_recs.len(), 2);
        assert_eq!(r2_recs.len(), 2);

        // Read 1: R1 sequence unchanged; ID gets :AAAAAAAA barcode
        assert_eq!(r1_recs[0].id, "@read1:AAAAAAAA");
        assert_eq!(r1_recs[0].seq, "ACGTACGTACGTACGT"); // unchanged
        assert_eq!(r1_recs[0].qual, "IIIIIIIIIIIIIIII");
        // R2: same barcode on ID; sequence clipped by 8bp
        assert_eq!(r2_recs[0].id, "@read1:AAAAAAAA");
        assert_eq!(r2_recs[0].seq, "GGGGGGGG"); // R2[8..]
        assert_eq!(r2_recs[0].qual.len(), r2_recs[0].seq.len());

        // Read 2: different barcode
        assert_eq!(r1_recs[1].id, "@read2:CCCCCCCC");
        assert_eq!(r2_recs[1].id, "@read2:CCCCCCCC");
        assert_eq!(r2_recs[1].seq, "GGGGGGGG");

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_implicon_custom_umi_length_6() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_implicon_umi6");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("foo_R1.fq");
        let r2_path = dir.join("foo_R2.fq");

        write_fastq(&r1_path, &[mk_rec("@r", "ACGTACGT", "IIIIIIII")])?;
        write_fastq(&r2_path, &[mk_rec("@r", "TTTTTTGGG", "IIIIIIIII")])?;

        implicon(&r1_path, &r2_path, 6, false, Some(&dir), 1)?;

        let out_r1 = dir.join("foo_6bp_UMI_R1.fastq");
        let out_r2 = dir.join("foo_6bp_UMI_R2.fastq");
        let r1_recs = read_fastq(&out_r1)?;
        let r2_recs = read_fastq(&out_r2)?;

        assert_eq!(r1_recs[0].id, "@r:TTTTTT"); // 6-bp barcode from R2
        assert_eq!(r1_recs[0].seq, "ACGTACGT"); // R1 untouched
        assert_eq!(r2_recs[0].id, "@r:TTTTTT");
        assert_eq!(r2_recs[0].seq, "GGG"); // R2[6..]

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }

    #[test]
    fn test_implicon_r2_too_short_errors() -> Result<()> {
        let dir = std::env::temp_dir().join("tg_test_implicon_short");
        std::fs::create_dir_all(&dir)?;
        let r1_path = dir.join("short_R1.fq");
        let r2_path = dir.join("short_R2.fq");

        // R2 is 5bp, umi_length=8 → should bail
        write_fastq(&r1_path, &[mk_rec("@r", "ACGTACGT", "IIIIIIII")])?;
        write_fastq(&r2_path, &[mk_rec("@r", "ACGTA", "IIIII")])?;

        let result = implicon(&r1_path, &r2_path, 8, false, Some(&dir), 1);
        assert!(result.is_err(), "should bail when R2 < umi_length");
        let err = format!("{}", result.unwrap_err());
        assert!(
            err.contains("shorter") || err.contains("UMI"),
            "got: {}",
            err
        );

        std::fs::remove_dir_all(&dir).ok();
        Ok(())
    }
}
