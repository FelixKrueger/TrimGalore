//! `--clump_only`: lossless reorder-only specialty mode.
//!
//! Reorders FASTQ records by canonical 16-mer minimizer for gzip-friendly
//! compression, WITHOUT any trimming, filtering, adapter detection, or
//! record modification. Load-bearing invariant: every input record appears
//! in the output byte-identically (header, sequence, quality). Only
//! file-level order changes.
//!
//! Contract-scope note: the byte-identity claim covers the three
//! semantically-meaningful fields (header, sequence, quality). The
//! plus-line (line 3) is normalized to bare `+` on output and CRLF line
//! endings are normalized to LF — both codebase-wide `FastqReader` and
//! `FastqWriter` behaviours, not `--clump_only`-specific.
//!
//! Implementation shape (v1 — single-threaded):
//! 1. Stream records from the input FASTQ via `FastqReader`.
//! 2. For each record, compute `canonical_minimizer` and dispatch into a
//!    per-bin buffer keyed by `bin_for(minimizer, n_bins)`.
//! 3. When a bin's raw-bytes accumulator crosses the byte budget, sort
//!    the bin via `clump::sort_single_by_key` (stable, content-tiebreaker
//!    cascade) and write the sorted records as a single gzip member
//!    (or a plain-text chunk under `--dont_gzip`) appended to the
//!    output file.
//! 4. At EOF, flush remaining bins in bin-index order for determinism.
//!
//! Concatenated gzip members form a valid `.gz` file per RFC 1952. The
//! gzip-window reset at each member boundary is what makes clumping
//! effective: sorted records within a bin share long minimizer-anchored
//! substrings, so gzip's 32 KB sliding window finds intra-bin redundancy
//! that would be lost across bin boundaries.
//!
//! Parallelism is deferred to v1.1: single-threaded v1 delivers both
//! correctness and byte-identity; `--cores` is accepted for interface
//! parity but does not affect the current implementation. See issue
//! #353 for background.

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::clump::{
    self, MinimizerKey, canonical_minimizer, estimated_record_bytes, sort_paired_by_key,
    sort_single_by_key,
};
use crate::fastq::{FastqReader, FastqRecord};
use crate::fastqc;
use crate::format::{InputFormat, detect_input_format};
use crate::io as naming;

/// Reorder-only statistics — deliberately narrower than `TrimStats`.
///
/// Used by `write_clump_only_report` to render the reorder-only report at
/// `<stem>_clumping_report.txt`. No JSON companion (deliberately outside
/// nf-core's `*_trimming_report.*` scan glob; the report is short enough
/// to grep).
#[derive(Debug, Default)]
pub struct ClumpOnlyStats {
    pub total_records: u64,
    pub input_bytes: u64,
    pub output_bytes: u64,
    pub n_bins: usize,
    pub peak_bin_occupancy: usize,
    pub compression_level: u32,
    pub input_compressed: bool,
    pub output_compressed: bool,
}

/// Per-bin buffer for the single-end dispatcher.
///
/// Mirrors `parallel::PairedBin`'s shape: pre-`reserve_exact` on first
/// push to avoid the `Vec` doubling cascade that inflates peak memory
/// when many bins grow simultaneously.
#[derive(Default)]
struct SingleBin {
    records: Vec<FastqRecord>,
    keys: Vec<MinimizerKey>,
    raw_bytes: usize,
    budget: usize,
}

impl SingleBin {
    fn with_budget(budget: usize) -> Self {
        Self {
            budget,
            ..Self::default()
        }
    }

    fn push(&mut self, r: FastqRecord, key: MinimizerKey) {
        let bytes_this = estimated_record_bytes(&r);
        if self.keys.capacity() == 0 {
            let predicted = self.budget.div_ceil(bytes_this.max(1)).max(1);
            self.records.reserve_exact(predicted);
            self.keys.reserve_exact(predicted);
        }
        self.raw_bytes += bytes_this;
        self.records.push(r);
        self.keys.push(key);
    }

    fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    fn take(&mut self) -> (Vec<FastqRecord>, Vec<MinimizerKey>) {
        let records = std::mem::take(&mut self.records);
        let keys = std::mem::take(&mut self.keys);
        self.raw_bytes = 0;
        (records, keys)
    }
}

/// Per-bin buffer for the paired-end dispatcher. Mirrors `SingleBin`
/// but keeps R1/R2 in lockstep.
#[derive(Default)]
struct PairedBin {
    r1: Vec<FastqRecord>,
    r2: Vec<FastqRecord>,
    keys: Vec<MinimizerKey>,
    raw_bytes: usize,
    budget: usize,
}

impl PairedBin {
    fn with_budget(budget: usize) -> Self {
        Self {
            budget,
            ..Self::default()
        }
    }

    fn push(&mut self, r1: FastqRecord, r2: FastqRecord, key: MinimizerKey) {
        let bytes_this = estimated_record_bytes(&r1) + estimated_record_bytes(&r2);
        if self.keys.capacity() == 0 {
            let predicted = self.budget.div_ceil(bytes_this.max(1)).max(1);
            self.r1.reserve_exact(predicted);
            self.r2.reserve_exact(predicted);
            self.keys.reserve_exact(predicted);
        }
        self.raw_bytes += bytes_this;
        self.r1.push(r1);
        self.r2.push(r2);
        self.keys.push(key);
    }

    fn is_empty(&self) -> bool {
        self.r1.is_empty()
    }

    fn take(&mut self) -> (Vec<FastqRecord>, Vec<FastqRecord>, Vec<MinimizerKey>) {
        let r1 = std::mem::take(&mut self.r1);
        let r2 = std::mem::take(&mut self.r2);
        let keys = std::mem::take(&mut self.keys);
        self.raw_bytes = 0;
        (r1, r2, keys)
    }
}

/// Write a single bin's records to `out` as either one gzip member
/// (concatenatable per RFC 1952) or a plain-text chunk under `--dont_gzip`.
///
/// Byte-count return: number of bytes written to the output stream (not
/// the record count), used by the caller for report stats.
fn flush_bin_single<W: Write>(
    bin: &mut SingleBin,
    out: &mut W,
    gzip: bool,
    compression: u32,
) -> Result<u64> {
    let (mut records, mut keys) = bin.take();
    sort_single_by_key(&mut records, &mut keys);
    let bytes = write_records_member(&records, out, gzip, compression)?;
    Ok(bytes)
}

/// Paired-end variant. Sorts both mates in lockstep and writes each mate's
/// records as one gzip member to its own output. Returns `(r1_bytes, r2_bytes)`.
fn flush_bin_paired<W1: Write, W2: Write>(
    bin: &mut PairedBin,
    out_r1: &mut W1,
    out_r2: &mut W2,
    gzip: bool,
    compression: u32,
) -> Result<(u64, u64)> {
    let (mut r1, mut r2, mut keys) = bin.take();
    sort_paired_by_key(&mut r1, &mut r2, &mut keys);
    let b1 = write_records_member(&r1, out_r1, gzip, compression)?;
    let b2 = write_records_member(&r2, out_r2, gzip, compression)?;
    Ok((b1, b2))
}

/// Encode `records` as a single gzip member (or a plain-text chunk) and
/// write to `out`. Returns the number of output bytes written.
///
/// Under `gzip=true`, a fresh `GzEncoder` is created for these records
/// and finished before writing; the result is a self-contained gzip
/// member that concatenates cleanly with prior/subsequent members.
fn write_records_member<W: Write>(
    records: &[FastqRecord],
    out: &mut W,
    gzip: bool,
    compression: u32,
) -> Result<u64> {
    if gzip {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::new(compression));
        for rec in records {
            rec.write_to(&mut encoder)?;
        }
        let member = encoder.finish()?;
        let n = member.len() as u64;
        out.write_all(&member)?;
        Ok(n)
    } else {
        let mut n: u64 = 0;
        for rec in records {
            // Track bytes via a counting wrapper would be cleaner, but for
            // plain-text output the record's on-disk footprint is exactly
            // `estimated_record_bytes` — cheap to compute and matches what
            // write_to emits (id + '\n' + seq + '\n' + "+\n" + qual + '\n').
            let bytes_this = estimated_record_bytes(rec) as u64;
            rec.write_to(out)?;
            n += bytes_this;
        }
        Ok(n)
    }
}

/// Reject uBAM input for `--clump_only` (v1 is FASTQ in/out only; uBAM
/// in/out is deferred to v2 per PLAN §Resolved decisions).
fn reject_ubam(input: &Path) -> Result<()> {
    match detect_input_format(input)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => Ok(()),
        InputFormat::UnalignedBam => bail!(
            "uBAM input is not yet supported under --clump_only \
             (v1 is FASTQ in / FASTQ out only; see #353 for the v2 follow-up). \
             Input: {}",
            input.display()
        ),
    }
}

/// Run `--clump_only` on a single-end FASTQ input.
///
/// Byte-identity: every record R in the input file appears in the output
/// file with `R.id`, `R.seq`, `R.qual` byte-identical. Only file-level
/// order changes (records reordered by canonical minimizer, with a stable
/// content-tiebreaker cascade for cross-run determinism).
#[allow(clippy::too_many_arguments)]
pub fn clump_only_single(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip_output: bool,
    cores: usize,
    memory_budget_bytes: u64,
    compression: u32,
    fastqc: bool,
    fastqc_args: Option<&str>,
    no_report_file: bool,
) -> Result<ClumpOnlyStats> {
    reject_ubam(input)?;

    // Layout sizing uses a floor of 1 core (v1 is single-threaded).
    let layout = clump::resolve_layout(memory_budget_bytes, cores.max(1))?;

    let input_compressed = naming::is_gzipped(input);
    let output_path = naming::clumped_output_name(input, output_dir, basename, gzip_output);

    eprintln!(
        "clump-only: reordering '{}' -> '{}' ({} bins × {} MB budget, gzip level {}{})",
        input.display(),
        output_path.display(),
        layout.n_bins,
        layout.bin_byte_budget / (1024 * 1024),
        compression,
        if gzip_output { "" } else { ", --dont_gzip" },
    );

    let mut reader = FastqReader::open(input)
        .with_context(|| format!("Failed to open input: {}", input.display()))?;

    let mut bins: Vec<SingleBin> = (0..layout.n_bins)
        .map(|_| SingleBin::with_budget(layout.bin_byte_budget))
        .collect();

    let mut out = BufWriter::new(
        File::create(&output_path)
            .with_context(|| format!("Failed to create output: {}", output_path.display()))?,
    );

    let mut stats = ClumpOnlyStats {
        n_bins: layout.n_bins,
        compression_level: compression,
        input_compressed,
        output_compressed: gzip_output,
        ..Default::default()
    };
    let mut output_bytes: u64 = 0;

    // Streaming loop
    while let Some(record) = reader.next_record()? {
        let key = canonical_minimizer(record.seq.as_bytes());
        let bin_idx = clump::bin_for(key, layout.n_bins);
        bins[bin_idx].push(record, key);
        stats.total_records += 1;

        if bins[bin_idx].raw_bytes >= layout.bin_byte_budget {
            let occ = bins[bin_idx].records.len();
            if occ > stats.peak_bin_occupancy {
                stats.peak_bin_occupancy = occ;
            }
            output_bytes +=
                flush_bin_single(&mut bins[bin_idx], &mut out, gzip_output, compression)?;
        }
    }

    // Final flush: bin-index order for determinism
    for bin in bins.iter_mut() {
        if !bin.is_empty() {
            let occ = bin.records.len();
            if occ > stats.peak_bin_occupancy {
                stats.peak_bin_occupancy = occ;
            }
            output_bytes += flush_bin_single(bin, &mut out, gzip_output, compression)?;
        }
    }

    // Empty-input case: ensure the output file is a valid gzip stream
    // (empty file → invalid .gz that downstream `zcat`/gzip decoders reject).
    // For plain output an empty file is fine. Writing an empty gzip member
    // produces a valid 0-record .gz file.
    if output_bytes == 0 && gzip_output {
        output_bytes += write_records_member(&[], &mut out, gzip_output, compression)?;
    }

    out.flush()?;
    drop(out);

    stats.input_bytes = std::fs::metadata(input)
        .with_context(|| format!("Failed to stat input {}", input.display()))?
        .len();
    stats.output_bytes = output_bytes;

    eprintln!(
        "clump-only: wrote {} records in {} bins (peak {} records/bin), {} bytes -> {} bytes",
        stats.total_records,
        stats.n_bins,
        stats.peak_bin_occupancy,
        stats.input_bytes,
        stats.output_bytes,
    );

    // Reorder-only report (skipped when the caller passed --no_report_file).
    if !no_report_file {
        let report_path = naming::clumping_report_name(input, output_dir);
        write_clump_only_report(&stats, &report_path, input, &output_path)
            .with_context(|| format!("Failed to write report: {}", report_path.display()))?;
    }

    // Optional FastQC — new plumbing (specialty modes previously returned
    // before the fastqc::run call in main.rs).
    if fastqc {
        fastqc::run(&output_path, fastqc_args, output_dir, cores.max(1))?;
    }

    Ok(stats)
}

/// Run `--clump_only` on a paired-end FASTQ input pair.
///
/// Byte-identity is enforced per-record for each mate; pair lockstep is
/// preserved by sorting on R1's minimizer and reordering R2 in lockstep
/// (`sort_paired_by_key`).
#[allow(clippy::too_many_arguments)]
pub fn clump_only_paired(
    input_r1: &Path,
    input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip_output: bool,
    cores: usize,
    memory_budget_bytes: u64,
    compression: u32,
    fastqc: bool,
    fastqc_args: Option<&str>,
    no_report_file: bool,
) -> Result<ClumpOnlyStats> {
    reject_ubam(input_r1)?;
    reject_ubam(input_r2)?;

    let layout = clump::resolve_layout(memory_budget_bytes, cores.max(1))?;

    let input_compressed = naming::is_gzipped(input_r1) || naming::is_gzipped(input_r2);
    let (out_r1_path, out_r2_path) =
        naming::clumped_paired_output_names(input_r1, input_r2, output_dir, basename, gzip_output);

    eprintln!(
        "clump-only (paired): '{}' + '{}' -> '{}' + '{}' ({} bins × {} MB, gzip level {}{})",
        input_r1.display(),
        input_r2.display(),
        out_r1_path.display(),
        out_r2_path.display(),
        layout.n_bins,
        layout.bin_byte_budget / (1024 * 1024),
        compression,
        if gzip_output { "" } else { ", --dont_gzip" },
    );

    let mut reader_r1 = FastqReader::open(input_r1)
        .with_context(|| format!("Failed to open R1: {}", input_r1.display()))?;
    let mut reader_r2 = FastqReader::open(input_r2)
        .with_context(|| format!("Failed to open R2: {}", input_r2.display()))?;

    let mut bins: Vec<PairedBin> = (0..layout.n_bins)
        .map(|_| PairedBin::with_budget(layout.bin_byte_budget))
        .collect();

    let mut out_r1 = BufWriter::new(
        File::create(&out_r1_path)
            .with_context(|| format!("Failed to create R1 output: {}", out_r1_path.display()))?,
    );
    let mut out_r2 = BufWriter::new(
        File::create(&out_r2_path)
            .with_context(|| format!("Failed to create R2 output: {}", out_r2_path.display()))?,
    );

    let mut stats = ClumpOnlyStats {
        n_bins: layout.n_bins,
        compression_level: compression,
        input_compressed,
        output_compressed: gzip_output,
        ..Default::default()
    };
    let mut out_bytes_r1: u64 = 0;
    let mut out_bytes_r2: u64 = 0;

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;
        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                let key = canonical_minimizer(r1.seq.as_bytes());
                let bin_idx = clump::bin_for(key, layout.n_bins);
                bins[bin_idx].push(r1, r2, key);
                stats.total_records += 1;
                if bins[bin_idx].raw_bytes >= layout.bin_byte_budget {
                    let occ = bins[bin_idx].r1.len();
                    if occ > stats.peak_bin_occupancy {
                        stats.peak_bin_occupancy = occ;
                    }
                    let (b1, b2) = flush_bin_paired(
                        &mut bins[bin_idx],
                        &mut out_r1,
                        &mut out_r2,
                        gzip_output,
                        compression,
                    )?;
                    out_bytes_r1 += b1;
                    out_bytes_r2 += b2;
                }
            }
            (None, None) => break,
            (Some(_), None) => bail!(
                "Read 2 file is truncated — R1 has more reads than R2. \
                 Please check your paired-end input files!"
            ),
            (None, Some(_)) => bail!(
                "Read 1 file is truncated — R2 has more reads than R1. \
                 Please check your paired-end input files!"
            ),
        }
    }

    for bin in bins.iter_mut() {
        if !bin.is_empty() {
            let occ = bin.r1.len();
            if occ > stats.peak_bin_occupancy {
                stats.peak_bin_occupancy = occ;
            }
            let (b1, b2) =
                flush_bin_paired(bin, &mut out_r1, &mut out_r2, gzip_output, compression)?;
            out_bytes_r1 += b1;
            out_bytes_r2 += b2;
        }
    }

    // Empty-input case: emit an empty gzip member per mate so the output
    // files are valid `.gz` streams (see clump_only_single for rationale).
    if out_bytes_r1 == 0 && gzip_output {
        out_bytes_r1 += write_records_member(&[], &mut out_r1, gzip_output, compression)?;
    }
    if out_bytes_r2 == 0 && gzip_output {
        out_bytes_r2 += write_records_member(&[], &mut out_r2, gzip_output, compression)?;
    }

    out_r1.flush()?;
    out_r2.flush()?;
    drop(out_r1);
    drop(out_r2);

    let in_bytes_r1 = std::fs::metadata(input_r1)?.len();
    let in_bytes_r2 = std::fs::metadata(input_r2)?.len();
    stats.input_bytes = in_bytes_r1 + in_bytes_r2;
    stats.output_bytes = out_bytes_r1 + out_bytes_r2;

    eprintln!(
        "clump-only (paired): wrote {} pairs in {} bins (peak {} pairs/bin)",
        stats.total_records, stats.n_bins, stats.peak_bin_occupancy,
    );

    // Per-mate reports (mirrors --clumpify's per-input report convention).
    // Skipped when the caller passed --no_report_file.
    if !no_report_file {
        let r1_report = naming::clumping_report_name(input_r1, output_dir);
        let r2_report = naming::clumping_report_name(input_r2, output_dir);
        let r1_stats = ClumpOnlyStats {
            input_bytes: in_bytes_r1,
            output_bytes: out_bytes_r1,
            ..stats_shape_from(&stats)
        };
        let r2_stats = ClumpOnlyStats {
            input_bytes: in_bytes_r2,
            output_bytes: out_bytes_r2,
            ..stats_shape_from(&stats)
        };
        write_clump_only_report(&r1_stats, &r1_report, input_r1, &out_r1_path)?;
        write_clump_only_report(&r2_stats, &r2_report, input_r2, &out_r2_path)?;
    }

    if fastqc {
        fastqc::run(&out_r1_path, fastqc_args, output_dir, cores.max(1))?;
        fastqc::run(&out_r2_path, fastqc_args, output_dir, cores.max(1))?;
    }

    Ok(stats)
}

/// Copy the shared-across-mates fields (total_records, n_bins, peak,
/// compression, input/output compressed flags) so each per-mate report
/// carries the shape without hand-duplicating the assignment site.
fn stats_shape_from(src: &ClumpOnlyStats) -> ClumpOnlyStats {
    ClumpOnlyStats {
        total_records: src.total_records,
        input_bytes: 0,
        output_bytes: 0,
        n_bins: src.n_bins,
        peak_bin_occupancy: src.peak_bin_occupancy,
        compression_level: src.compression_level,
        input_compressed: src.input_compressed,
        output_compressed: src.output_compressed,
    }
}

/// Write the reorder-only text report. Format follows PLAN.md §Report shape:
/// header lines + a per-bin summary + a compression-ratio line that is
/// deliberately omitted when input/output compression states differ (gzip↔plain),
/// which would produce misleading numbers.
pub fn write_clump_only_report(
    stats: &ClumpOnlyStats,
    txt_path: &Path,
    input: &Path,
    output: &Path,
) -> Result<()> {
    let mut w = BufWriter::new(File::create(txt_path)?);
    writeln!(w, "Trim Galore version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "Mode: --clump_only (lossless reorder)")?;
    writeln!(
        w,
        "Input:  {} ({}, {} bytes)",
        input.display(),
        if stats.input_compressed {
            "gzip"
        } else {
            "plain"
        },
        stats.input_bytes,
    )?;
    writeln!(
        w,
        "Output: {} ({}, {} bytes)",
        output.display(),
        if stats.output_compressed {
            format!("gzip level {}", stats.compression_level)
        } else {
            "plain".to_string()
        },
        stats.output_bytes,
    )?;
    writeln!(w, "Records: {}", stats.total_records)?;
    writeln!(
        w,
        "Bins: {} (peak occupancy {} records)",
        stats.n_bins, stats.peak_bin_occupancy,
    )?;
    // Compression ratio only when both sides are gzip AND the on-disk
    // input bytes aren't a proxy for something else. When --dont_gzip
    // flipped output to plain, the ratio would be misleading; when input
    // was plain, likewise.
    if stats.input_compressed && stats.output_compressed && stats.output_bytes > 0 {
        let ratio = stats.input_bytes as f64 / stats.output_bytes as f64;
        writeln!(w, "Compression ratio: {:.2}x", ratio)?;
    }
    w.flush()?;
    Ok(())
}

/// Utility for callers that want to inspect the expected output path
/// (e.g. for output-collision pre-flight in `main.rs`).
pub fn expected_output_paths_single(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip_output: bool,
) -> PathBuf {
    naming::clumped_output_name(input, output_dir, basename, gzip_output)
}

/// Utility for callers that want to inspect the expected paired output
/// paths (both mates).
pub fn expected_output_paths_paired(
    input_r1: &Path,
    input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip_output: bool,
) -> (PathBuf, PathBuf) {
    naming::clumped_paired_output_names(input_r1, input_r2, output_dir, basename, gzip_output)
}

// ────────────────────────────── Tests ──────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use std::path::PathBuf;

    /// Build a synthetic FASTQ file at `path` (gzip if extension ends `.gz`)
    /// from the given records.
    fn write_synthetic_fastq(path: &Path, records: &[FastqRecord]) -> Result<()> {
        let is_gz = naming::is_gzipped(path);
        let file = File::create(path)?;
        if is_gz {
            let mut w = BufWriter::new(GzEncoder::new(file, Compression::new(6)));
            for r in records {
                r.write_to(&mut w)?;
            }
            w.flush()?;
        } else {
            let mut w = BufWriter::new(file);
            for r in records {
                r.write_to(&mut w)?;
            }
            w.flush()?;
        }
        Ok(())
    }

    /// Read every record from a plain-or-gzip FASTQ file back into a Vec.
    fn read_all(path: &Path) -> Result<Vec<FastqRecord>> {
        let mut reader = FastqReader::open(path)?;
        let mut out = Vec::new();
        while let Some(r) = reader.next_record()? {
            out.push(r);
        }
        Ok(out)
    }

    fn rec(id: &str, seq: &str, qual: &str) -> FastqRecord {
        FastqRecord {
            id: id.to_string(),
            seq: seq.to_string(),
            qual: qual.to_string(),
        }
    }

    /// Sort records into a canonical order (by id) for multiset comparison —
    /// the mode is a *permutation*, so sort-both-sides-then-compare is
    /// the load-bearing test. Returns tuples (id, seq, qual) since
    /// `FastqRecord` doesn't implement `PartialEq`.
    fn multiset_sorted(recs: Vec<FastqRecord>) -> Vec<(String, String, String)> {
        let mut out: Vec<(String, String, String)> =
            recs.into_iter().map(|r| (r.id, r.seq, r.qual)).collect();
        out.sort();
        out
    }

    /// Deterministic pseudo-random ACGT sequence generator that doesn't
    /// depend on `rand`. Same shape as the RRBS/clumpify test helpers.
    fn synth_records(n: usize, seq_len: usize) -> Vec<FastqRecord> {
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let seq: String = (0..seq_len)
                .map(|j: usize| {
                    let r = (i as u64).wrapping_mul(2654435761).wrapping_add(j as u64);
                    (*b"ACGT")[(r >> 28) as usize & 0x3] as char
                })
                .collect();
            let qual = "I".repeat(seq_len);
            out.push(rec(&format!("@read_{i}"), &seq, &qual));
        }
        out
    }

    fn tempdir(tag: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("tg_clump_only_test_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// Small memory budget large enough to pass the `resolve_layout`
    /// floor at cores=1.
    fn small_memory_budget() -> u64 {
        // resolve_layout requires STATIC_OVERHEAD (512 MiB) + enough for
        // MIN_BIN_BYTES per bin. Pick 1 GiB to comfortably exceed the floor.
        1024 * 1024 * 1024
    }

    // ── Load-bearing invariants ──────────────────────────────────────

    #[test]
    fn test_clump_only_single_permutation() -> Result<()> {
        let dir = tempdir("perm");
        let input = dir.join("in.fq.gz");
        let inputs = synth_records(500, 80);
        write_synthetic_fastq(&input, &inputs)?;

        let _ = clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        let output_path = naming::clumped_output_name(&input, Some(&dir), None, true);
        let outputs = read_all(&output_path)?;
        assert_eq!(outputs.len(), inputs.len(), "record count preserved");
        assert_eq!(
            multiset_sorted(inputs),
            multiset_sorted(outputs),
            "record multiset preserved byte-identically"
        );
        Ok(())
    }

    #[test]
    fn test_clump_only_paired_lockstep() -> Result<()> {
        let dir = tempdir("pair");
        let r1_path = dir.join("s_R1.fq.gz");
        let r2_path = dir.join("s_R2.fq.gz");
        let r1_inputs: Vec<FastqRecord> = (0..300)
            .map(|i| rec(&format!("@pair_{i}/1"), &"A".repeat(80), &"I".repeat(80)))
            .collect();
        let r2_inputs: Vec<FastqRecord> = (0..300)
            .map(|i| rec(&format!("@pair_{i}/2"), &"C".repeat(80), &"I".repeat(80)))
            .collect();
        write_synthetic_fastq(&r1_path, &r1_inputs)?;
        write_synthetic_fastq(&r2_path, &r2_inputs)?;

        let _ = clump_only_paired(
            &r1_path,
            &r2_path,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        let (out_r1, out_r2) =
            naming::clumped_paired_output_names(&r1_path, &r2_path, Some(&dir), None, true);
        let outputs_r1 = read_all(&out_r1)?;
        let outputs_r2 = read_all(&out_r2)?;

        // Pair lockstep: at every output index i, R1[i].id and R2[i].id
        // must share the same numeric suffix.
        for (i, (a, b)) in outputs_r1.iter().zip(outputs_r2.iter()).enumerate() {
            let a_num: usize =
                a.id.trim_start_matches("@pair_")
                    .trim_end_matches("/1")
                    .parse()
                    .unwrap();
            let b_num: usize =
                b.id.trim_start_matches("@pair_")
                    .trim_end_matches("/2")
                    .parse()
                    .unwrap();
            assert_eq!(a_num, b_num, "pair lockstep violated at output index {}", i);
        }
        // Per-mate multiset preserved
        assert_eq!(multiset_sorted(r1_inputs), multiset_sorted(outputs_r1));
        assert_eq!(multiset_sorted(r2_inputs), multiset_sorted(outputs_r2));
        Ok(())
    }

    #[test]
    fn test_clump_only_deterministic() -> Result<()> {
        // Two runs on the same input produce byte-identical output.
        let dir = tempdir("det");
        let input = dir.join("in.fq.gz");
        write_synthetic_fastq(&input, &synth_records(200, 60))?;

        let out1 = dir.join("run1");
        let out2 = dir.join("run2");
        std::fs::create_dir_all(&out1)?;
        std::fs::create_dir_all(&out2)?;

        clump_only_single(
            &input,
            Some(&out1),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        clump_only_single(
            &input,
            Some(&out2),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        let p1 = naming::clumped_output_name(&input, Some(&out1), None, true);
        let p2 = naming::clumped_output_name(&input, Some(&out2), None, true);
        let mut b1 = Vec::new();
        let mut b2 = Vec::new();
        File::open(&p1)?.read_to_end(&mut b1)?;
        File::open(&p2)?.read_to_end(&mut b2)?;
        assert_eq!(b1, b2, "cross-run byte-identity");
        Ok(())
    }

    #[test]
    fn test_clump_only_no_trimming_short_reads_kept() -> Result<()> {
        // Reads shorter than default --length filter (20) must be
        // preserved — --clump_only does not filter.
        let dir = tempdir("nolen");
        let input = dir.join("in.fq.gz");
        let inputs: Vec<FastqRecord> = (0usize..50)
            .map(|i| {
                let seq_len = if i.is_multiple_of(3) { 5 } else { 40 };
                rec(
                    &format!("@r_{i}"),
                    &"A".repeat(seq_len),
                    &"I".repeat(seq_len),
                )
            })
            .collect();
        write_synthetic_fastq(&input, &inputs)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        let outputs = read_all(&naming::clumped_output_name(&input, Some(&dir), None, true))?;
        assert_eq!(inputs.len(), outputs.len(), "no length filter applied");
        Ok(())
    }

    #[test]
    fn test_clump_only_no_adapter_detection() -> Result<()> {
        // Reads carrying obvious Illumina adapter sequences emerge with
        // the adapter bytes intact.
        let dir = tempdir("noadapt");
        let input = dir.join("in.fq.gz");
        let adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let inputs: Vec<FastqRecord> = (0..30)
            .map(|i| {
                let payload = "A".repeat(30);
                let seq = format!("{}{}", payload, adapter);
                rec(&format!("@r_{i}"), &seq, &"I".repeat(seq.len()))
            })
            .collect();
        write_synthetic_fastq(&input, &inputs)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        let outputs = read_all(&naming::clumped_output_name(&input, Some(&dir), None, true))?;
        // Every output record still carries the adapter.
        for r in &outputs {
            assert!(r.seq.ends_with(adapter), "adapter suffix lost: {}", r.seq);
        }
        Ok(())
    }

    #[test]
    fn test_clump_only_empty_input() -> Result<()> {
        let dir = tempdir("empty");
        let input = dir.join("in.fq.gz");
        write_synthetic_fastq(&input, &[])?;

        let stats = clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        assert_eq!(stats.total_records, 0);
        let outputs = read_all(&naming::clumped_output_name(&input, Some(&dir), None, true))?;
        assert_eq!(outputs.len(), 0);
        // Report exists.
        let report = naming::clumping_report_name(&input, Some(&dir));
        assert!(report.exists(), "report file created for empty input");
        Ok(())
    }

    #[test]
    fn test_clump_only_single_record() -> Result<()> {
        let dir = tempdir("one");
        let input = dir.join("in.fq.gz");
        let inputs = vec![rec("@only_read", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII")];
        write_synthetic_fastq(&input, &inputs)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        let outputs = read_all(&naming::clumped_output_name(&input, Some(&dir), None, true))?;
        assert_eq!(outputs.len(), 1);
        assert_eq!(outputs[0].id, inputs[0].id);
        assert_eq!(outputs[0].seq, inputs[0].seq);
        assert_eq!(outputs[0].qual, inputs[0].qual);
        Ok(())
    }

    #[test]
    fn test_clump_only_plain_output_no_gzip_extension() -> Result<()> {
        // --dont_gzip produces `*_clumped.fq` (no .gz), still byte-identical
        // records.
        let dir = tempdir("plain");
        let input = dir.join("in.fq.gz");
        let inputs = synth_records(100, 60);
        write_synthetic_fastq(&input, &inputs)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            false, // gzip_output = false
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;
        let out = naming::clumped_output_name(&input, Some(&dir), None, false);
        assert!(out.exists(), "plain output file exists: {}", out.display());
        assert_eq!(out.extension().and_then(|s| s.to_str()), Some("fq"));

        // Records preserved as a multiset.
        let outputs = read_all(&out)?;
        assert_eq!(inputs.len(), outputs.len());
        assert_eq!(multiset_sorted(inputs), multiset_sorted(outputs));
        Ok(())
    }

    #[test]
    fn test_report_filename_distinct_from_trimming_report() -> Result<()> {
        // Load-bearing: the clump-only report filename does NOT match
        // the `*_trimming_report.*` glob that nf-core / MultiQC scan.
        let dir = tempdir("report");
        let input = dir.join("in.fq.gz");
        write_synthetic_fastq(&input, &synth_records(20, 40))?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        let clumping_report = naming::clumping_report_name(&input, Some(&dir));
        let trimming_report = naming::report_name(&input, Some(&dir));
        assert!(clumping_report.exists(), "clumping report created");
        assert!(
            !trimming_report.exists(),
            "trimming report must NOT be created (nf-core scan protection)"
        );
        Ok(())
    }

    // ── Contract-scope normalizations (plus-line + CRLF) ──────────────
    // These guard the Resolved Decision 6 carve-outs: byte-identity is
    // scoped to header + sequence + quality; line 3 is normalized to
    // bare `+` on output, and CRLF endings are normalized to LF. Both
    // behaviours are inherited from `FastqReader`/`FastqWriter` — the
    // tests exist as regression guards so a future reader/writer change
    // doesn't silently drift the contract we documented in `--help`.

    /// Write a raw FASTQ file (bytes verbatim, bypassing `FastqRecord::write_to`)
    /// so tests can construct input with `+<header-repeat>` on line 3 or CRLF endings.
    fn write_raw_fastq_gz(path: &Path, contents: &[u8]) -> Result<()> {
        let file = File::create(path)?;
        let mut encoder = GzEncoder::new(file, Compression::new(6));
        encoder.write_all(contents)?;
        encoder.finish()?;
        Ok(())
    }

    #[test]
    fn test_clump_only_normalizes_plus_line() -> Result<()> {
        // Input with `+@read_1` on line 3 (an uncommon but valid FASTQ variant).
        // Output must have bare `+\n` on line 3 — behaviour inherited from
        // FastqWriter, which hard-codes the plus-line. Byte-identity claim in
        // Resolved Decision 6 explicitly excludes the plus-line.
        let dir = tempdir("plus_line");
        let input = dir.join("in.fq.gz");
        let raw = b"@read_1\nACGTACGTACGTACGT\n+@read_1\nIIIIIIIIIIIIIIII\n\
                    @read_2\nGGGGGGGGGGGGGGGG\n+@read_2\nJJJJJJJJJJJJJJJJ\n";
        write_raw_fastq_gz(&input, raw)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        // Read back the raw bytes of the decompressed output and check
        // that lines 3, 7, ... are all bare `+`.
        let out_path = naming::clumped_output_name(&input, Some(&dir), None, true);
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(File::open(&out_path)?).read_to_end(&mut decoded)?;
        let text = std::str::from_utf8(&decoded).unwrap();
        let lines: Vec<&str> = text.split('\n').collect();
        // Every 4th line (indexes 2, 6, 10, …) must be `+`, not `+@read_N`.
        for (idx, chunk) in lines.chunks(4).enumerate() {
            if chunk.len() < 3 || chunk[0].is_empty() {
                break; // trailing empty
            }
            assert_eq!(
                chunk[2], "+",
                "record {idx}: plus-line not normalized to bare `+`; got: {:?}",
                chunk[2]
            );
        }
        Ok(())
    }

    #[test]
    fn test_clump_only_normalizes_crlf() -> Result<()> {
        // Input with CRLF line endings — reader strips `\r`, writer emits LF only.
        // Byte-identity in Resolved Decision 6 explicitly covers only header +
        // sequence + quality bytes; CRLF → LF is codebase-wide normalization.
        let dir = tempdir("crlf");
        let input = dir.join("in.fq.gz");
        let raw = b"@read_1\r\nACGTACGTACGTACGT\r\n+\r\nIIIIIIIIIIIIIIII\r\n\
                    @read_2\r\nGGGGGGGGGGGGGGGG\r\n+\r\nJJJJJJJJJJJJJJJJ\r\n";
        write_raw_fastq_gz(&input, raw)?;

        clump_only_single(
            &input,
            Some(&dir),
            None,
            true,
            1,
            small_memory_budget(),
            6,
            false,
            None,
            false,
        )?;

        let out_path = naming::clumped_output_name(&input, Some(&dir), None, true);
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(File::open(&out_path)?).read_to_end(&mut decoded)?;
        // No `\r` bytes should survive the round-trip.
        assert!(
            !decoded.contains(&b'\r'),
            "CRLF was not normalized to LF; output still contains \\r bytes"
        );
        // Record contents (header/seq/qual) preserved.
        let text = std::str::from_utf8(&decoded).unwrap();
        assert!(text.contains("@read_1"));
        assert!(text.contains("ACGTACGTACGTACGT"));
        assert!(text.contains("GGGGGGGGGGGGGGGG"));
        Ok(())
    }
}
