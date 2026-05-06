//! Post-trimming demultiplexing based on 3' inline barcodes.
//!
//! After adapter/quality trimming, reads are split into per-sample files
//! based on barcode sequences at the 3' end of each read. The barcode is
//! removed from the sequence and appended to the read ID as `_BC:<barcode>`.
//! Reads that don't match any barcode go into a "NoCode" file.
//!
//! Single-end only (matches TrimGalore behavior).

use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::fastq::{FastqReader, FastqWriter};

/// A parsed barcode entry: barcode sequence → sample name.
#[derive(Debug)]
pub struct BarcodeEntry {
    pub barcode: String,
    pub sample_name: String,
}

/// Read and validate a barcode file (TSV: sample_name\tbarcode_sequence).
///
/// Returns the list of barcodes and validates that:
/// - Each line has exactly 2 tab-separated fields
/// - Barcode sequences contain only A, C, T, G, N
/// - All barcodes have the same length
pub fn read_barcode_file(path: &Path) -> Result<Vec<BarcodeEntry>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open barcode file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut entries = Vec::new();
    let mut expected_len: Option<usize> = None;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r').trim_end();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.splitn(2, '\t').collect();
        if parts.len() != 2 {
            bail!(
                "Barcode file {} line {}: expected tab-separated 'sample_name\\tbarcode_sequence', got: '{}'",
                path.display(),
                line_num + 1,
                line
            );
        }

        let sample_name = parts[0].to_string();
        let barcode = parts[1].to_string();

        // Validate barcode characters
        if !barcode
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'T' | b'G' | b'N'))
        {
            bail!(
                "Barcode file {} line {}: barcode must contain only A, C, T, G, N. Got: '{}'",
                path.display(),
                line_num + 1,
                barcode
            );
        }

        // Validate consistent barcode length
        match expected_len {
            None => expected_len = Some(barcode.len()),
            Some(len) if barcode.len() != len => {
                bail!(
                    "Barcode file {} line {}: barcode '{}' has length {} but expected {} (all barcodes must be the same length)",
                    path.display(),
                    line_num + 1,
                    barcode,
                    barcode.len(),
                    len
                );
            }
            _ => {}
        }

        eprintln!("Testing:\t{}\t{}", sample_name, barcode);
        entries.push(BarcodeEntry {
            barcode,
            sample_name,
        });
    }

    if entries.is_empty() {
        bail!("Barcode file {} is empty", path.display());
    }

    eprintln!(
        "Demultiplexing file {} seems to be in the correct format. Proceeding...",
        path.display()
    );
    Ok(entries)
}

/// Demultiplex a trimmed FASTQ file based on 3' inline barcodes.
///
/// For each read:
/// 1. Extract the last `barcode_length` bases from the sequence
/// 2. Look up the barcode in the sample map (exact match)
/// 3. Remove the barcode from sequence and quality
/// 4. Replace spaces in read ID with underscores, append `_BC:<barcode>`
/// 5. Write to the matching sample's output file (or NoCode)
///
/// Also writes a summary file with per-barcode counts.
pub fn demultiplex(
    trimmed_file: &Path,
    barcodes: &[BarcodeEntry],
    gzip: bool,
    output_dir: Option<&Path>,
    cores: usize,
    gzip_level: u32,
) -> Result<()> {
    let barcode_length = barcodes[0].barcode.len();
    eprintln!("Setting barcode length to {}", barcode_length);

    // Derive base name for output files:
    // trimmed file is e.g. "sample_trimmed.fq.gz" → strip .gz then .fq
    let trimmed_name = trimmed_file
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let mut base_name = trimmed_name.clone();
    if base_name.ends_with(".gz") {
        base_name = base_name[..base_name.len() - 3].to_string();
    }
    if base_name.ends_with(".fq") {
        base_name = base_name[..base_name.len() - 3].to_string();
    }

    let dir = output_dir.map(|d| d.to_path_buf()).unwrap_or_else(|| {
        trimmed_file
            .parent()
            .unwrap_or(Path::new("."))
            .to_path_buf()
    });

    // Open per-sample output writers + NoCode
    let mut writers: HashMap<String, FastqWriter> = HashMap::new();
    for entry in barcodes {
        let filename = output_filename(&base_name, &entry.sample_name, gzip);
        let path = dir.join(&filename);
        writers.insert(
            entry.barcode.clone(),
            FastqWriter::create(&path, gzip, cores, gzip_level)?,
        );
    }
    // NoCode writer for unmatched barcodes
    let nocode_filename = output_filename(&base_name, "NoCode", gzip);
    let nocode_path = dir.join(&nocode_filename);
    let mut nocode_writer = FastqWriter::create(&nocode_path, gzip, cores, gzip_level)?;

    // Per-barcode counts
    let mut counts: HashMap<String, usize> =
        barcodes.iter().map(|e| (e.barcode.clone(), 0)).collect();
    counts.insert("NoCode".to_string(), 0);

    // Process the trimmed file
    let mut reader = FastqReader::open(trimmed_file)?;
    let mut total: usize = 0;

    while let Some(mut record) = reader.next_record()? {
        total += 1;
        if total.is_multiple_of(10_000_000) {
            eprintln!("{} sequences processed", total);
        }

        if record.seq.len() < barcode_length {
            // Read too short for barcode — goes to NoCode
            let bc_tag = &record.seq;
            record.id = format!("{}_BC:{}", record.id.replace(' ', "_"), bc_tag);
            record.seq.clear();
            record.qual.clear();
            nocode_writer.write_record(&record)?;
            *counts.get_mut("NoCode").unwrap() += 1;
            continue;
        }

        // Extract barcode from 3' end
        let bc_start = record.seq.len() - barcode_length;
        let barcode_seq = record.seq[bc_start..].to_string();

        // Trim barcode from sequence and quality
        record.seq.truncate(bc_start);
        record.qual.truncate(bc_start);

        // Modify read ID: replace spaces with underscores, append barcode
        record.id = format!("{}_BC:{}", record.id.replace(' ', "_"), barcode_seq);

        // Write to matching sample or NoCode
        if let Some(writer) = writers.get_mut(barcode_seq.as_str()) {
            *counts.get_mut(barcode_seq.as_str()).unwrap() += 1;
            writer.write_record(&record)?;
        } else {
            *counts.get_mut("NoCode").unwrap() += 1;
            nocode_writer.write_record(&record)?;
        }
    }

    // Flush all writers
    for writer in writers.values_mut() {
        writer.flush()?;
    }
    nocode_writer.flush()?;
    drop(writers);
    drop(nocode_writer);

    eprintln!("Processed sequences from file >{trimmed_name}< in total: {total}");

    // Build sorted summary (descending by count)
    let mut sorted_counts: Vec<_> = counts.iter().collect();
    sorted_counts.sort_by(|a, b| b.1.cmp(a.1));

    // Build sample name lookup for display (barcode → sample_name)
    let display_map: HashMap<&str, &str> = barcodes
        .iter()
        .map(|e| (e.barcode.as_str(), e.sample_name.as_str()))
        .chain(std::iter::once(("NoCode", "NoCode")))
        .collect();

    // Print to stderr
    eprintln!("\nDemultiplexing summary\n======================");
    for (barcode, count) in &sorted_counts {
        let sample = display_map.get(barcode.as_str()).unwrap_or(&"unknown");
        eprintln!("{}\t{}\t{}", barcode, sample, count);
    }
    eprintln!("======================\n");

    // Write summary file
    let summary_path = dir.join(format!("{}_demultiplexing_summary.txt", base_name));
    let file = File::create(&summary_path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "Demultiplexing summary")?;
    writeln!(w, "======================")?;
    for (barcode, count) in &sorted_counts {
        let sample = display_map.get(barcode.as_str()).unwrap_or(&"unknown");
        writeln!(w, "{}\t{}\t{}", barcode, sample, count)?;
    }
    writeln!(w, "======================")?;
    writeln!(w, "      ¯\\_(ツ)_/¯")?;
    w.flush()?;

    Ok(())
}

/// Generate the output filename for a demux sample.
fn output_filename(base_name: &str, sample_name: &str, gzip: bool) -> String {
    if gzip {
        format!("{}_{}.fq.gz", base_name, sample_name)
    } else {
        format!("{}_{}.fq", base_name, sample_name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::Write;
    use std::path::PathBuf;

    /// Set up a unique tempdir for a test, removing any leftover from a
    /// prior aborted run. Mirrors the project's existing test convention
    /// (`std::env::temp_dir().join("tg_test_*")`).
    fn fresh_tmpdir(slug: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// §5.5 regression: Windows-authored barcode samplesheets use CRLF
    /// (`\r\n`) line endings. Perl strips the trailing `\r` before
    /// parsing; Rust's `read_barcode_file` does the same via
    /// `trim_end_matches('\r').trim_end()`. This test locks down that
    /// behaviour — without the trim, the barcode field would include a
    /// trailing `\r`, fail the ACGTN-only validator, and bail with a
    /// confusing "barcode must contain only A, C, T, G, N" error.
    #[test]
    fn test_read_barcode_file_handles_crlf_line_endings() {
        let dir = fresh_tmpdir("tg_demux_crlf_samplesheet");
        let path = dir.join("barcodes.tsv");
        // Three barcodes with explicit CRLF, plus a trailing empty line
        // to exercise the empty-line branch.
        let mut f = File::create(&path).unwrap();
        f.write_all(b"sample1\tACGTACGT\r\nsample2\tGTCAGTCA\r\nsample3\tTGCATGCA\r\n\r\n")
            .unwrap();
        f.flush().unwrap();

        let entries = read_barcode_file(&path).expect("CRLF samplesheet must parse");

        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0].sample_name, "sample1");
        assert_eq!(entries[0].barcode, "ACGTACGT");
        assert_eq!(entries[1].sample_name, "sample2");
        assert_eq!(entries[1].barcode, "GTCAGTCA");
        assert_eq!(entries[2].sample_name, "sample3");
        assert_eq!(entries[2].barcode, "TGCATGCA");
        // Critical: no stray '\r' must survive on the parsed barcode.
        for e in &entries {
            assert!(
                !e.barcode.contains('\r'),
                "barcode '{}' contains stray \\r — CRLF strip failed",
                e.barcode.escape_default()
            );
        }

        fs::remove_dir_all(&dir).unwrap();
    }

    /// §5.6 regression: when a read is shorter than the barcode length,
    /// `demultiplex` (src/demux.rs:178-187) routes it to the NoCode
    /// bucket rather than slicing past the read end (which would either
    /// panic on a non-UTF8 boundary or compare against an empty string).
    /// End-to-end via `demultiplex` because the routing decision is
    /// inline with file I/O — easier to verify by inspecting NoCode
    /// output than to refactor for in-isolation testability.
    #[test]
    fn test_demultiplex_short_reads_route_to_nocode() {
        let dir = fresh_tmpdir("tg_demux_short_nocode");
        let trimmed = dir.join("input_trimmed.fq");

        // Three reads: one too-short (5 bp < 8 bp barcode), one normal
        // length but non-matching barcode, one normal length matching
        // sample1's barcode at the 3' end.
        let mut f = File::create(&trimmed).unwrap();
        // Record 1: 5 bp — must route to NoCode (too short).
        // Record 2: 16 bp ending in ZZZZZZZZ — invalid bases, but the
        //           extracted 8 bp tail is "GGGGTTTT" → no match → NoCode.
        // Record 3: 16 bp ending in ACGTACGT — matches sample1's barcode.
        f.write_all(
            b"@short\nACGTA\n+\nIIIII\n\
              @nomatch\nNNNNNNNNGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n\
              @hit\nNNNNNNNNACGTACGT\n+\nIIIIIIIIIIIIIIII\n",
        )
        .unwrap();
        f.flush().unwrap();

        let barcodes = vec![BarcodeEntry {
            sample_name: "sample1".to_string(),
            barcode: "ACGTACGT".to_string(),
        }];

        demultiplex(
            &trimmed,
            &barcodes,
            false,
            Some(&dir),
            1,
            crate::fastq::DEFAULT_GZIP_LEVEL,
        )
        .unwrap();

        // Output base is the trimmed-file basename minus .gz/.fq suffixes:
        // "input_trimmed.fq" → "input_trimmed" + "_<sample>.fq".
        let nocode = dir.join("input_trimmed_NoCode.fq");
        let sample1 = dir.join("input_trimmed_sample1.fq");

        let nocode_text = fs::read_to_string(&nocode).unwrap();
        let sample1_text = fs::read_to_string(&sample1).unwrap();

        // Both the too-short read and the non-matching one land in NoCode.
        assert!(
            nocode_text.contains("@short_BC:"),
            "short read must route to NoCode. NoCode contents:\n{nocode_text}"
        );
        assert!(
            nocode_text.contains("@nomatch_BC:"),
            "non-matching read must route to NoCode. NoCode contents:\n{nocode_text}"
        );
        // The short read's seq + qual are cleared (no slice past EOR).
        let short_block = nocode_text
            .lines()
            .skip_while(|l| !l.starts_with("@short_BC:"))
            .take(4)
            .collect::<Vec<_>>();
        assert_eq!(
            short_block.len(),
            4,
            "short record must still emit a 4-line FASTQ block"
        );
        assert_eq!(short_block[1], "", "short read seq must be cleared");
        assert_eq!(short_block[3], "", "short read qual must be cleared");

        // The matching read lands in the sample bucket.
        assert!(
            sample1_text.contains("@hit_BC:ACGTACGT"),
            "matching read must route to sample1 bucket. sample1 contents:\n{sample1_text}"
        );

        fs::remove_dir_all(&dir).unwrap();
    }
}
