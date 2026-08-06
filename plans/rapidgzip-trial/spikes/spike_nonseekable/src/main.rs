//! Does rapidgzip-core v0.2.0 decode gzip from a non-seekable source, correctly?
//!
//! Decoded bytes go to **stdout** so the driver can `cmp` them against the
//! original; diagnostics go to **stderr**. Exit 0 on success, 1 on a decode
//! error, 2 on misuse.
//!
//! Modes:
//!   open <path>   `Decoder::open`, which routes by `supports_positional_reads`
//!   stream        `Decoder::stream_reader(io::stdin())`, forced sequential

use rapidgzip_core::Decoder;
use std::io::{self, Read, Write};

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_default();
    let threads: usize = std::env::var("RGZ_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(8);
    let decoder = match Decoder::builder().decoder_threads(threads).build() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("BUILD_ERROR {e}");
            std::process::exit(2);
        }
    };

    let opened = match mode.as_str() {
        "open" => match std::env::args().nth(2) {
            Some(p) => decoder.open(&p),
            None => {
                eprintln!("usage: spike open <path>");
                std::process::exit(2);
            }
        },
        "stream" => decoder.stream_reader(io::stdin()),
        _ => {
            eprintln!("usage: spike <open <path>|stream>");
            std::process::exit(2);
        }
    };

    let mut reader = match opened {
        Ok(r) => r,
        Err(e) => {
            eprintln!("OPEN_ERROR {e}");
            std::process::exit(1);
        }
    };
    let handle = reader.handle();

    // Stream to stdout rather than buffering, so a multi-MB fixture also
    // exercises the incremental path rather than one big read.
    let stdout = io::stdout();
    let mut out = stdout.lock();
    let mut buf = vec![0u8; 64 * 1024];
    let mut total: u64 = 0;
    // Sampled mid-stream: after the reader retires, both counters read 0, so a
    // post-EOF sample says nothing about which path ran.
    let mut mid: Option<(usize, usize)> = None;
    loop {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => {
                total += n as u64;
                if mid.is_none() && total > 8 * 1024 * 1024 {
                    let s = handle.stats();
                    mid = Some((s.active_workers, s.spawned_workers));
                }
                if out.write_all(&buf[..n]).is_err() {
                    eprintln!("STDOUT_ERROR");
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("READ_ERROR {e}");
                eprintln!("bytes_before_error={total}");
                std::process::exit(1);
            }
        }
    }
    let _ = out.flush();

    let stats = handle.stats();
    match mid {
        Some((a, s)) => eprintln!("mid_active={a} mid_spawned={s}"),
        None => eprintln!("mid_active=n/a mid_spawned=n/a (input under 8 MB)"),
    }
    eprintln!(
        "end_active={} end_spawned={} configured={}",
        stats.active_workers, stats.spawned_workers, threads
    );

    match reader.finish() {
        Ok(report) => {
            eprintln!("OK bytes={total} members={}", report.member_count);
        }
        Err(e) => {
            eprintln!("FINISH_ERROR {e}");
            std::process::exit(1);
        }
    }
}
