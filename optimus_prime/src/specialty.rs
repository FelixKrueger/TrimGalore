//! Specialty trimming modes that bypass the normal trimming pipeline.
//!
//! These modes (--hardtrim5, --hardtrim3, --clock, --implicon) process
//! input files with simple fixed operations and exit immediately.

use anyhow::{bail, Result};
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
    eprintln!("Writing hard-trimmed (first {}bp) version of '{}' to '{}'",
        keep, input.display(), output_path.display());

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
    eprintln!("Writing hard-trimmed (last {}bp) version of '{}' to '{}'",
        keep, input.display(), output_path.display());

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

    eprintln!("Writing Clock-processed version of '{}' to '{}'", input_r1.display(), out1.display());
    eprintln!("Writing Clock-processed version of '{}' to '{}'", input_r2.display(), out2.display());

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
                        count, r1_seq.len(), r2_seq.len()
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
    eprintln!("thereof had fixed sequence CAGT in both R1 and R2: {} ({}%)\n", filtered_count, perc);

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

    eprintln!("Writing UMI-trimmed version of '{}' to '{}'", input_r1.display(), out1.display());
    eprintln!("Writing UMI-trimmed version of '{}' to '{}'", input_r2.display(), out2.display());

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
                        count, r2.seq.len(), umi_length
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

fn hardtrim_output_name(input: &Path, keep: usize, end: &str, output_dir: Option<&Path>, gzip: bool) -> PathBuf {
    let stem = naming::strip_fastq_extensions(input);
    let ext = if gzip { ".fq.gz" } else { ".fq" };
    let filename = format!("{}.{}bp_{}{}", stem, keep, end, ext);
    match output_dir {
        Some(dir) => dir.join(filename),
        None => PathBuf::from(filename),
    }
}

fn clock_output_name(input: &Path, read: &str, output_dir: Option<&Path>, gzip: bool) -> PathBuf {
    let stem = naming::strip_fastq_extensions(input);
    let ext = if gzip { ".fq.gz" } else { ".fq" };
    let filename = format!("{}.clock_UMI.{}{}", stem, read, ext);
    match output_dir {
        Some(dir) => dir.join(filename),
        None => PathBuf::from(filename),
    }
}

fn implicon_output_name(input: &Path, umi_len: usize, read: &str, output_dir: Option<&Path>, gzip: bool) -> PathBuf {
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
