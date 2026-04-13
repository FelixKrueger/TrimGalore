//! Post-trimming demultiplexing based on 3' inline barcodes.
//!
//! After adapter/quality trimming, reads are split into per-sample files
//! based on barcode sequences at the 3' end of each read. The barcode is
//! removed from the sequence and appended to the read ID as `_BC:<barcode>`.
//! Reads that don't match any barcode go into a "NoCode" file.
//!
//! Single-end only (matches TrimGalore behavior).

use anyhow::{bail, Context, Result};
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
                path.display(), line_num + 1, line
            );
        }

        let sample_name = parts[0].to_string();
        let barcode = parts[1].to_string();

        // Validate barcode characters
        if !barcode.bytes().all(|b| matches!(b, b'A' | b'C' | b'T' | b'G' | b'N')) {
            bail!(
                "Barcode file {} line {}: barcode must contain only A, C, T, G, N. Got: '{}'",
                path.display(), line_num + 1, barcode
            );
        }

        // Validate consistent barcode length
        match expected_len {
            None => expected_len = Some(barcode.len()),
            Some(len) if barcode.len() != len => {
                bail!(
                    "Barcode file {} line {}: barcode '{}' has length {} but expected {} (all barcodes must be the same length)",
                    path.display(), line_num + 1, barcode, barcode.len(), len
                );
            }
            _ => {}
        }

        eprintln!("Testing:\t{}\t{}", sample_name, barcode);
        entries.push(BarcodeEntry { barcode, sample_name });
    }

    if entries.is_empty() {
        bail!("Barcode file {} is empty", path.display());
    }

    eprintln!("Demultiplexing file {} seems to be in the correct format. Proceeding...", path.display());
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

    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| trimmed_file.parent().unwrap_or(Path::new(".")).to_path_buf());

    // Open per-sample output writers + NoCode
    let mut writers: HashMap<String, FastqWriter> = HashMap::new();
    for entry in barcodes {
        let filename = output_filename(&base_name, &entry.sample_name, gzip);
        let path = dir.join(&filename);
        writers.insert(
            entry.barcode.clone(),
            FastqWriter::create(&path, gzip, cores)?,
        );
    }
    // NoCode writer for unmatched barcodes
    let nocode_filename = output_filename(&base_name, "NoCode", gzip);
    let nocode_path = dir.join(&nocode_filename);
    let mut nocode_writer = FastqWriter::create(&nocode_path, gzip, cores)?;

    // Per-barcode counts
    let mut counts: HashMap<String, usize> = barcodes
        .iter()
        .map(|e| (e.barcode.clone(), 0))
        .collect();
    counts.insert("NoCode".to_string(), 0);

    // Process the trimmed file
    let mut reader = FastqReader::open(trimmed_file)?;
    let mut total: usize = 0;

    while let Some(mut record) = reader.next_record()? {
        total += 1;
        if total % 10_000_000 == 0 {
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
    write!(w, "      ¯\\_(ツ)_/¯\n")?;
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
