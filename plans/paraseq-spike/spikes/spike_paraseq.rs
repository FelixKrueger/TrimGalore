//! plans/paraseq-spike — Throwaway perf benchmark.
//!
//! Question: does swapping paraseq in for our line-by-line `String`-allocating
//! parser net out faster, AFTER the moment we have to materialize owned records
//! to ship them across our worker-pool channel?
//!
//! Compares three loop bodies over the same fixture:
//!   (A) Baseline — mirrors src/fastq.rs:322-355 (BufReader::read_line into
//!       String, three Strings per record).
//!   (B) paraseq zero-copy — iterate RecordSet borrowed records, no alloc.
//!   (C) paraseq + owned — convert each borrowed record into an owned
//!       FastqRecord-shaped struct (the realistic integration cost).
//!
//! Each variant computes the same sentinel value (sum of seq lengths +
//! sum of last byte of each seq) so the optimizer can't elide parsing.
//!
//! Run (from this directory):
//!   cargo run --release -- build-fixture <num_reads> <out.fastq>
//!   cargo run --release -- bench <in.fastq>
//!
//! For the recorded results, fixture = 500,000 reads × 150bp synthesised
//! from the repo's 10K_150bp fixture by tiling.

use anyhow::{Context, Result, bail};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

const BUF_SIZE: usize = 64 * 1024;

// ─── Owned record shape — mirrors crate::fastq::FastqRecord ────────────────
#[derive(Debug)]
#[allow(dead_code)]
struct FastqRecord {
    id: String,
    seq: String,
    qual: String,
}

// ─── (A) Baseline parser — mirrors src/fastq.rs:322-355 exactly ────────────
fn baseline_parse(path: &Path) -> Result<(u64, u64, u64)> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(BUF_SIZE, file);
    let mut line_buf = String::with_capacity(512);
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;

    loop {
        line_buf.clear();
        if reader.read_line(&mut line_buf)? == 0 {
            break;
        }
        let _id = line_buf.trim_end_matches(['\n', '\r']).to_string();

        line_buf.clear();
        if reader.read_line(&mut line_buf)? == 0 {
            bail!("truncated: missing seq");
        }
        let seq = line_buf.trim_end_matches(['\n', '\r']).to_string();

        line_buf.clear();
        if reader.read_line(&mut line_buf)? == 0 {
            bail!("truncated: missing +");
        }

        line_buf.clear();
        if reader.read_line(&mut line_buf)? == 0 {
            bail!("truncated: missing qual");
        }
        let _qual = line_buf.trim_end_matches(['\n', '\r']).to_string();

        count += 1;
        seq_len_sum += seq.len() as u64;
        if let Some(&b) = seq.as_bytes().last() {
            last_byte_sum = last_byte_sum.wrapping_add(b as u64);
        }
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── (B) paraseq zero-copy ────────────────────────────────────────────────
fn paraseq_zero_copy(path: &Path) -> Result<(u64, u64, u64)> {
    use paraseq::Record;
    use paraseq::fastq;

    let mut reader = fastq::Reader::from_path(path)?;
    let mut record_set = reader.new_record_set();
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;

    while record_set.fill(&mut reader)? {
        for record in record_set.iter() {
            let record = record?;
            let seq = record.seq_raw();
            count += 1;
            seq_len_sum += seq.len() as u64;
            if let Some(&b) = seq.last() {
                last_byte_sum = last_byte_sum.wrapping_add(b as u64);
            }
        }
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── (C) paraseq → owned conversion (naïve: format! + str::from_utf8.to_string)
fn paraseq_owned_naive(path: &Path) -> Result<(u64, u64, u64)> {
    use paraseq::Record;
    use paraseq::fastq;

    let mut reader = fastq::Reader::from_path(path)?;
    let mut record_set = reader.new_record_set();
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;

    while record_set.fill(&mut reader)? {
        let mut batch: Vec<FastqRecord> = Vec::with_capacity(4096);
        for record in record_set.iter() {
            let record = record?;
            let id_bytes = record.id_str().as_bytes();
            let seq_bytes = record.seq_raw();
            let qual_bytes = record.qual().unwrap_or(b"");
            let owned = FastqRecord {
                id: format!("@{}", std::str::from_utf8(id_bytes).unwrap_or("")),
                seq: std::str::from_utf8(seq_bytes).unwrap_or("").to_string(),
                qual: std::str::from_utf8(qual_bytes).unwrap_or("").to_string(),
            };
            count += 1;
            seq_len_sum += owned.seq.len() as u64;
            if let Some(&b) = owned.seq.as_bytes().last() {
                last_byte_sum = last_byte_sum.wrapping_add(b as u64);
            }
            batch.push(owned);
        }
        // Simulate channel hand-off: drop the batch (caller would send across
        // a thread). The drop runs the String destructors so the work is
        // included in wall-clock.
        drop(batch);
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── (D) paraseq → owned conversion (lean: prepend '@' manually, From<Vec<u8>>)
//
// Avoids `format!` (which goes through Display/Write machinery) and the
// `&str → String::to_string()` round-trip. Each field is built from a
// pre-sized Vec<u8> with a single memcpy, then wrapped via
// `String::from_utf8(...)` (one UTF-8 validation pass, no extra alloc).
fn paraseq_owned_lean(path: &Path) -> Result<(u64, u64, u64)> {
    use paraseq::Record;
    use paraseq::fastq;

    let mut reader = fastq::Reader::from_path(path)?;
    let mut record_set = reader.new_record_set();
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;

    while record_set.fill(&mut reader)? {
        let mut batch: Vec<FastqRecord> = Vec::with_capacity(4096);
        for record in record_set.iter() {
            let record = record?;
            let id_bytes = record.id_str().as_bytes();
            let seq_bytes = record.seq_raw();
            let qual_bytes = record.qual().unwrap_or(b"");

            // id: prepend '@' manually so the call site emits a single
            // memcpy + a 1-byte unshift, not a Display-machinery round-trip.
            let mut id_v = Vec::with_capacity(id_bytes.len() + 1);
            id_v.push(b'@');
            id_v.extend_from_slice(id_bytes);
            let id = String::from_utf8(id_v).unwrap_or_default();

            let seq = String::from_utf8(seq_bytes.to_vec()).unwrap_or_default();
            let qual = String::from_utf8(qual_bytes.to_vec()).unwrap_or_default();

            let owned = FastqRecord { id, seq, qual };
            count += 1;
            seq_len_sum += owned.seq.len() as u64;
            if let Some(&b) = owned.seq.as_bytes().last() {
                last_byte_sum = last_byte_sum.wrapping_add(b as u64);
            }
            batch.push(owned);
        }
        drop(batch);
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── Timer helper — runs `f` `n_runs` times, returns ms per run ────────────
fn time_runs<F>(name: &str, n_runs: usize, mut f: F) -> Result<()>
where
    F: FnMut() -> Result<(u64, u64, u64)>,
{
    let mut timings_ms: Vec<f64> = Vec::with_capacity(n_runs);
    let mut last_result: Option<(u64, u64, u64)> = None;
    for _ in 0..n_runs {
        let t0 = Instant::now();
        let r = f()?;
        let elapsed = t0.elapsed().as_secs_f64() * 1000.0;
        timings_ms.push(elapsed);
        last_result = Some(r);
    }
    let median = {
        let mut s = timings_ms.clone();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        s[s.len() / 2]
    };
    let min = timings_ms.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = timings_ms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let (count, seq_sum, last_sum) = last_result.unwrap();
    println!(
        "{:30}  median={:7.1} ms   min={:7.1} ms   max={:7.1} ms   reads={}  seq_sum={}  sentinel={}",
        name, median, min, max, count, seq_sum, last_sum
    );
    Ok(())
}

// ─── Fixture builder — tile a small fastq up to N reads ────────────────────
fn build_fixture(num_reads: usize, out_path: &Path) -> Result<()> {
    // Source: repo root's 10K_150bp fixture, already-decompressed via gzip -dc
    // into a Vec<u8> in-memory. To keep this script dep-light we don't pull
    // in flate2 — instead we shell out the user. Simpler: synthesise from
    // scratch (deterministic content, no I/O wrangling).
    //
    // Each synthesised read:
    //   @read_<i>
    //   <150bp pseudo-random ACGT from a deterministic LFSR>
    //   +
    //   <150 chars of Phred '!'+33+i%30 → "BCDE..."-ish>
    let mut out = std::io::BufWriter::with_capacity(BUF_SIZE, File::create(out_path)?);
    let bases: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut lfsr: u32 = 0xACE1u32;
    let mut seq_buf = [0u8; 150];
    let mut qual_buf = [0u8; 150];
    for i in 0..num_reads {
        for j in 0..150 {
            // Galois LFSR step → 1 random bit per shift; consume 2 bits per base.
            for _ in 0..2 {
                let lsb = lfsr & 1;
                lfsr >>= 1;
                if lsb != 0 {
                    lfsr ^= 0xB400;
                }
            }
            seq_buf[j] = bases[(lfsr & 0x3) as usize];
            qual_buf[j] = b'!' + ((j as u8 + (i & 0x1F) as u8) % 40);
        }
        writeln!(out, "@read_{i} 1:N:0:CGATCG")?;
        out.write_all(&seq_buf)?;
        out.write_all(b"\n+\n")?;
        out.write_all(&qual_buf)?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        bail!("usage: spike_paraseq <build-fixture <n> <out.fastq> | bench <in.fastq>>");
    }
    match args[1].as_str() {
        "build-fixture" => {
            let n: usize = args.get(2).context("missing num_reads")?.parse()?;
            let out = Path::new(args.get(3).context("missing out path")?);
            eprintln!("Building fixture: {n} reads → {}", out.display());
            build_fixture(n, out)?;
            let bytes = std::fs::metadata(out)?.len();
            eprintln!("Done. {bytes} bytes.");
        }
        "bench" => {
            let path = Path::new(args.get(2).context("missing input path")?);
            let bytes = std::fs::metadata(path)?.len();
            eprintln!("Benching {} ({} bytes, plain FASTQ)", path.display(), bytes);
            eprintln!("Each variant runs 3× — reporting median / min / max.\n");

            time_runs("(A) baseline String × 3 ", 3, || baseline_parse(path))?;
            time_runs("(B) paraseq zero-copy   ", 3, || paraseq_zero_copy(path))?;
            time_runs("(C) paraseq → owned naïve", 3, || paraseq_owned_naive(path))?;
            time_runs("(D) paraseq → owned lean ", 3, || paraseq_owned_lean(path))?;
        }
        other => bail!("unknown subcommand: {other}"),
    }
    Ok(())
}
