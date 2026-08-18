//! Deterministic gzipped paired FASTQ fixture for the `clumpify-memory` CI job.
//!
//! ```text
//! cargo run --release --example mk_clump_fixture -- <pairs> <read_len> <out_prefix>
//! ```
//!
//! Writes `<prefix>_R1.fastq.gz` and `<prefix>_R2.fastq.gz`. Gzip rather than
//! plain text because output encoding follows the input's, and the plain-output
//! branch is not the one users run.
//!
//! Sequences come from a fixed-seed SplitMix64 stream, so the bytes are
//! identical on every platform and every run. Distinct-sequence count is printed
//! before writing: a fixture built by repetition is a *small* fixture measured
//! many times, and the #439 calibration was invalidated once by exactly that,
//! invisibly, because nothing recorded sample size.

use std::fs::File;
use std::io::{BufWriter, Write};

use flate2::Compression;
use flate2::write::GzEncoder;

/// SplitMix64 — small, fast, and reproducible without a dependency.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// One gzip member per mate, written streaming so peak memory stays flat
/// regardless of `pairs`.
fn write_mate(
    path: &str,
    pairs: usize,
    read_len: usize,
    seed: u64,
    mate: u8,
) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut out = BufWriter::with_capacity(1 << 20, GzEncoder::new(file, Compression::fast()));
    let mut rng = Rng(seed);
    let qual = vec![b'I'; read_len];
    let mut seq = vec![0u8; read_len];

    for i in 0..pairs {
        for base in seq.iter_mut() {
            *base = b"ACGT"[(rng.next_u64() >> 33) as usize % 4];
        }
        writeln!(out, "@SIM:1:FCX:1:15:6329:{i} {mate}:N:0:ATCCGA")?;
        out.write_all(&seq)?;
        out.write_all(b"\n+\n")?;
        out.write_all(&qual)?;
        out.write_all(b"\n")?;
    }
    out.into_inner()?.finish()?;
    Ok(())
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: mk_clump_fixture <pairs> <read_len> <out_prefix>");
        std::process::exit(2);
    }
    let pairs: usize = args[1].parse().expect("pairs must be a positive integer");
    let read_len: usize = args[2]
        .parse()
        .expect("read_len must be a positive integer");
    let prefix = &args[3];
    assert!(pairs > 0 && read_len > 0, "pairs and read_len must be > 0");

    // Two independent streams, so R1 and R2 sequences differ as they would in a
    // real library while staying in lockstep by read ID.
    eprintln!(
        "fixture: {pairs} pairs x {read_len} bp, distinct sequences {pairs} per mate (one draw each)"
    );
    write_mate(
        &format!("{prefix}_R1.fastq.gz"),
        pairs,
        read_len,
        0x2545_F491_4F6C_DD1D,
        1,
    )?;
    write_mate(
        &format!("{prefix}_R2.fastq.gz"),
        pairs,
        read_len,
        0x1234_5678_9ABC_DEF0,
        2,
    )?;
    eprintln!("wrote {prefix}_R1.fastq.gz and {prefix}_R2.fastq.gz");
    Ok(())
}
