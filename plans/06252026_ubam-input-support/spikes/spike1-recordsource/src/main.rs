//! Spike 1 — RecordSource dispatch overhead.
//!
//! Question: does `Box<dyn RecordSource>`-style dynamic dispatch in the inner
//! record-read loop measurably regress wall-clock vs a concrete-type baseline
//! or an `enum { Fastq(...), Bam(...) }` static-dispatch shim?
//!
//! Compares three loop bodies over the same fixture. All three mirror the
//! TrimGalore FastqReader parser shape (src/fastq.rs:322-355): BufReader,
//! four read_line calls per record, three Strings per record.
//!
//!   (A) Concrete `FastqLikeReader` struct — no abstraction (baseline).
//!   (B) Same parser called through `Box<dyn RecordSource>` trait object.
//!   (C) Same parser called via `enum Reader::Fastq(...)` + match.
//!
//! Each variant computes the same sentinel (sum of seq lengths + sum of last
//! seq byte) so the optimizer can't elide parsing. Drop pattern is identical
//! across variants so any difference is dispatch-attributable.
//!
//! NOTE: this spike doesn't need a real BAM reader — the question is dispatch
//! overhead, and the parser body is identical across variants. We only need
//! one concrete impl to exercise the trait/enum machinery; the "BAM-shaped"
//! second variant would only matter if monomorphization decisions hinged on
//! a second type. Since `Box<dyn>` is a single vtable slot regardless of how
//! many concrete impls exist, one variant is enough to measure overhead.
//!
//! Run (from this directory):
//!   cargo run --release -- build-fixture <num_reads> <out.fastq>
//!   cargo run --release -- bench <in.fastq>

use anyhow::{Context, Result, bail};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

const BUF_SIZE: usize = 64 * 1024;

// ─── Owned record shape — mirrors crate::fastq::FastqRecord ─────────────────
#[derive(Debug)]
#[allow(dead_code)]
struct FastqRecord {
    id: String,
    seq: String,
    qual: String,
}

// ─── The actual parser body, shared verbatim across all three variants ─────
//
// Splitting this into a free function ensures all three variants run literally
// the same instructions inside the loop — the only thing that varies is how
// they're *reached* (direct call / vtable / match).
#[inline(always)]
fn parse_one_record(
    reader: &mut dyn BufRead,
    line_buf: &mut String,
) -> Result<Option<FastqRecord>> {
    line_buf.clear();
    if reader.read_line(line_buf)? == 0 {
        return Ok(None);
    }
    let id = line_buf.trim_end_matches(['\n', '\r']).to_string();

    line_buf.clear();
    if reader.read_line(line_buf)? == 0 {
        bail!("truncated: missing seq");
    }
    let seq = line_buf.trim_end_matches(['\n', '\r']).to_string();

    line_buf.clear();
    if reader.read_line(line_buf)? == 0 {
        bail!("truncated: missing +");
    }

    line_buf.clear();
    if reader.read_line(line_buf)? == 0 {
        bail!("truncated: missing qual");
    }
    let qual = line_buf.trim_end_matches(['\n', '\r']).to_string();

    Ok(Some(FastqRecord { id, seq, qual }))
}

// ─── (A) Concrete baseline — no abstraction ────────────────────────────────
struct ConcreteReader {
    reader: BufReader<File>,
    line_buf: String,
}

impl ConcreteReader {
    fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(BUF_SIZE, File::open(path)?),
            line_buf: String::with_capacity(512),
        })
    }

    #[inline]
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        parse_one_record(&mut self.reader, &mut self.line_buf)
    }
}

fn run_concrete(path: &Path) -> Result<(u64, u64, u64)> {
    let mut r = ConcreteReader::open(path)?;
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;
    while let Some(rec) = r.next_record()? {
        count += 1;
        seq_len_sum += rec.seq.len() as u64;
        if let Some(&b) = rec.seq.as_bytes().last() {
            last_byte_sum = last_byte_sum.wrapping_add(b as u64);
        }
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── (B) Box<dyn RecordSource> trait object ────────────────────────────────
trait RecordSource {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

struct TraitFastqReader {
    reader: BufReader<File>,
    line_buf: String,
}

impl TraitFastqReader {
    fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(BUF_SIZE, File::open(path)?),
            line_buf: String::with_capacity(512),
        })
    }
}

impl RecordSource for TraitFastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        parse_one_record(&mut self.reader, &mut self.line_buf)
    }
}

fn run_trait(path: &Path) -> Result<(u64, u64, u64)> {
    // Construction through the trait object — exactly how main.rs would
    // dispatch between FastqReader and BamReader.
    let mut r: Box<dyn RecordSource> = Box::new(TraitFastqReader::open(path)?);
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;
    while let Some(rec) = r.next_record()? {
        count += 1;
        seq_len_sum += rec.seq.len() as u64;
        if let Some(&b) = rec.seq.as_bytes().last() {
            last_byte_sum = last_byte_sum.wrapping_add(b as u64);
        }
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── (C) enum Reader { Fastq(...), Bam(...) } + match ─────────────────────
struct EnumFastqReader {
    reader: BufReader<File>,
    line_buf: String,
}

impl EnumFastqReader {
    fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(BUF_SIZE, File::open(path)?),
            line_buf: String::with_capacity(512),
        })
    }

    #[inline]
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        parse_one_record(&mut self.reader, &mut self.line_buf)
    }
}

// Placeholder so the enum has the second arm it would in production. The
// arm is never constructed in this spike — the goal is to expose the match
// in the inner loop, which is what production code would have.
#[allow(dead_code)]
struct EnumBamReader;

impl EnumBamReader {
    #[inline]
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        // Unreachable in this spike, but the match arm must be present so
        // the compiler must emit the discriminant check on every call.
        Ok(None)
    }
}

enum Reader {
    Fastq(EnumFastqReader),
    #[allow(dead_code)]
    Bam(EnumBamReader),
}

impl Reader {
    #[inline]
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        match self {
            Reader::Fastq(r) => r.next_record(),
            Reader::Bam(r) => r.next_record(),
        }
    }
}

fn run_enum(path: &Path) -> Result<(u64, u64, u64)> {
    let mut r = Reader::Fastq(EnumFastqReader::open(path)?);
    let mut count: u64 = 0;
    let mut seq_len_sum: u64 = 0;
    let mut last_byte_sum: u64 = 0;
    while let Some(rec) = r.next_record()? {
        count += 1;
        seq_len_sum += rec.seq.len() as u64;
        if let Some(&b) = rec.seq.as_bytes().last() {
            last_byte_sum = last_byte_sum.wrapping_add(b as u64);
        }
    }
    Ok((count, seq_len_sum, last_byte_sum))
}

// ─── Timer helper ─────────────────────────────────────────────────────────
fn time_runs<F>(name: &str, n_runs: usize, mut f: F) -> Result<f64>
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
        "{:32}  median={:7.2} ms   min={:7.2} ms   max={:7.2} ms   reads={}  seq_sum={}  sentinel={}",
        name, median, min, max, count, seq_sum, last_sum
    );
    Ok(median)
}

// ─── Fixture builder — same LFSR pattern as paraseq spike ─────────────────
fn build_fixture(num_reads: usize, out_path: &Path) -> Result<()> {
    let mut out = std::io::BufWriter::with_capacity(BUF_SIZE, File::create(out_path)?);
    let bases: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut lfsr: u32 = 0xACE1u32;
    let mut seq_buf = [0u8; 150];
    let mut qual_buf = [0u8; 150];
    for i in 0..num_reads {
        for j in 0..150 {
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
        bail!("usage: spike1_recordsource <build-fixture <n> <out.fastq> | bench <in.fastq>>");
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

            // Warm-up: prime the OS page cache so the first variant doesn't pay
            // the cold-cache penalty by itself.
            let _ = run_concrete(path)?;

            let m_a = time_runs("(A) concrete baseline           ", 3, || run_concrete(path))?;
            let m_b = time_runs("(B) Box<dyn RecordSource>       ", 3, || run_trait(path))?;
            let m_c = time_runs("(C) enum Reader + match         ", 3, || run_enum(path))?;

            println!();
            println!("Δ vs. concrete baseline (median):");
            println!("  (B) trait : {:+.2}%   ({:+.2} ms)", (m_b - m_a) / m_a * 100.0, m_b - m_a);
            println!("  (C) enum  : {:+.2}%   ({:+.2} ms)", (m_c - m_a) / m_a * 100.0, m_c - m_a);
            println!();
            println!("Success criterion: trait wins if ≤2% slower than (A).");
        }
        other => bail!("unknown subcommand: {other}"),
    }
    Ok(())
}
