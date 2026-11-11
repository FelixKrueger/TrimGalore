//! Spike (#457): does the per-allocation overhead ride the text's 0..B
//! oscillation, or is it resident in full like the spine?
//!
//! usage: spike_occupancy retain <pairs> <read_len>
//!        spike_occupancy cycle  <pairs> <read_len> <n_bins> <bin_budget_bytes>
//!
//! `retain` allocates, drops, and re-allocates, printing a marker and waiting on
//! stdin at each step so the driver can sample *current* RSS from outside.
//! `ru_maxrss` is a high-water mark and cannot see a decrease, which is the
//! whole question here.
//!
//! `cycle` mirrors the dispatcher: a pool of bins, each `reserve_exact`ed on the
//! first push after a `mem::take`, flushed at `raw_bytes >= budget`. Bins are
//! chosen pseudo-randomly, because round-robin would make every bin reach its
//! budget on the same record and flush in lockstep — an artefact that would
//! report peak occupancy of 1.0 by construction.

use std::io::{BufRead, Write};

use trim_galore::clump::{MinimizerKey, estimated_record_bytes};
use trim_galore::fastq::FastqRecord;

// An example is its own binary and does not inherit main.rs's allocator.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

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

#[derive(Default)]
struct Bin {
    r1: Vec<FastqRecord>,
    r2: Vec<FastqRecord>,
    keys: Vec<MinimizerKey>,
    raw_bytes: usize,
}

impl Bin {
    /// Mirrors `PairedBin::push` at parallel.rs:901 — one `reserve_exact` per
    /// bin per flush cycle, sized from the first record after the take.
    fn push(&mut self, a: FastqRecord, b: FastqRecord, key: MinimizerKey, budget: usize) {
        let bytes_this = estimated_record_bytes(&a) + estimated_record_bytes(&b);
        if self.keys.capacity() == 0 {
            let predicted = budget.div_ceil(bytes_this).max(1);
            self.r1.reserve_exact(predicted);
            self.r2.reserve_exact(predicted);
            self.keys.reserve_exact(predicted);
        }
        self.r1.push(a);
        self.r2.push(b);
        self.keys.push(key);
        self.raw_bytes += bytes_this;
        if self.raw_bytes >= budget {
            let _ = std::mem::take(&mut self.r1);
            let _ = std::mem::take(&mut self.r2);
            let _ = std::mem::take(&mut self.keys);
            self.raw_bytes = 0;
        }
    }
}

fn record(i: usize, mate: u8, seq: &str, qual: &str) -> FastqRecord {
    FastqRecord {
        id: format!("@SIM:1:FCX:1:15:6329:{i} {mate}:N:0:ATCCGA"),
        seq: seq.to_string(),
        qual: qual.to_string(),
    }
}

/// Print a marker and block until the driver has sampled RSS.
fn checkpoint(name: &str, live: usize) {
    println!("MARK={name} LIVE={live}");
    std::io::stdout().flush().unwrap();
    let mut line = String::new();
    std::io::stdin().lock().read_line(&mut line).unwrap();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("");
    let pairs: usize = args[2].parse().expect("pairs");
    let read_len: usize = args[3].parse().expect("read_len");

    let base = "ACGT".repeat(read_len.div_ceil(4));
    let seq = base[..read_len].to_string();
    let qual = "I".repeat(read_len);

    match mode {
        "retain" => {
            let mut r1: Vec<FastqRecord> = Vec::new();
            let mut r2: Vec<FastqRecord> = Vec::new();
            let mut keys: Vec<MinimizerKey> = Vec::new();
            r1.reserve_exact(pairs);
            r2.reserve_exact(pairs);
            keys.reserve_exact(pairs);
            checkpoint("empty", 0);

            for i in 0..pairs {
                r1.push(record(i, 1, &seq, &qual));
                r2.push(record(i, 2, &seq, &qual));
                keys.push((i % 4096) as MinimizerKey);
            }
            checkpoint("filled", r1.len());

            // The production flush: mem::take, then drop.
            let _ = std::mem::take(&mut r1);
            let _ = std::mem::take(&mut r2);
            let _ = std::mem::take(&mut keys);
            checkpoint("dropped", r1.len());

            r1.reserve_exact(pairs);
            r2.reserve_exact(pairs);
            keys.reserve_exact(pairs);
            for i in 0..pairs {
                r1.push(record(i, 1, &seq, &qual));
                r2.push(record(i, 2, &seq, &qual));
                keys.push((i % 4096) as MinimizerKey);
            }
            checkpoint("refilled", r1.len());
        }
        "cycle" => {
            let n_bins: usize = args[4].parse().expect("n_bins");
            let budget: usize = args[5].parse().expect("bin_budget_bytes");
            let mut bins: Vec<Bin> = (0..n_bins).map(|_| Bin::default()).collect();
            let mut rng = Rng(0x5DEE_CE66_D1CE_2025);
            let mut flushes = 0usize;
            for i in 0..pairs {
                let b = (rng.next_u64() >> 33) as usize % n_bins;
                let before = bins[b].raw_bytes;
                bins[b].push(record(i, 1, &seq, &qual), record(i, 2, &seq, &qual),
                             (i % 4096) as MinimizerKey, budget);
                if bins[b].raw_bytes < before {
                    flushes += 1;
                }
            }
            let live: usize = bins.iter().map(|b| b.r1.len()).sum();
            let resident: usize = bins.iter().map(|b| b.raw_bytes).sum();
            println!("MODE=cycle PAIRS={pairs} READ_LEN={read_len} N_BINS={n_bins} \
                      BUDGET={budget} FLUSHES={flushes} LIVE={live} RESIDENT={resident}");
        }
        _ => {
            eprintln!("usage: spike_occupancy retain <pairs> <read_len>");
            eprintln!("       spike_occupancy cycle <pairs> <read_len> <n_bins> <budget>");
            std::process::exit(2);
        }
    }
}
