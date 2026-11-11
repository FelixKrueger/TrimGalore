//! Spike (#457): resident bytes per paired-bin slot.
//!
//! Builds exactly the shape `PairedBin` holds — two `Vec<FastqRecord>` and one
//! `Vec<MinimizerKey>`, all `reserve_exact`ed once — and holds it resident past
//! the process peak. The driver runs this at several record counts and takes the
//! *slope*, so the process baseline cancels; an absolute RSS reading would
//! confound the per-record cost with the static term, which is the ambiguity
//! this spike exists to resolve.
//!
//! usage: spike_slot_cost <pairs> <read_len> [full|spine]
//!
//! `spine` pushes records whose three `String`s are empty, so nothing is heap
//! text and the slope is the `size_of` claim on its own: 72 + 72 + 4 = 148 B per
//! pair. `full - spine` is then heap text plus whatever the allocator adds.

use trim_galore::clump::{MinimizerKey, estimated_record_bytes};
use trim_galore::fastq::FastqRecord;

// An example is its own binary and does not inherit main.rs's allocator, so the
// shipping allocator has to be declared here or the spike measures libmalloc.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args.len() > 4 {
        eprintln!("usage: spike_slot_cost <pairs> <read_len> [full|spine]");
        std::process::exit(2);
    }
    let pairs: usize = args[1].parse().expect("pairs");
    let read_len: usize = args[2].parse().expect("read_len");
    let spine_only = args.get(3).map(|s| s == "spine").unwrap_or(false);

    let base = "ACGT".repeat(read_len.div_ceil(4));
    let seq = base[..read_len].to_string();
    let qual = "I".repeat(read_len);

    let mut r1: Vec<FastqRecord> = Vec::new();
    let mut r2: Vec<FastqRecord> = Vec::new();
    let mut keys: Vec<MinimizerKey> = Vec::new();
    r1.reserve_exact(pairs);
    r2.reserve_exact(pairs);
    keys.reserve_exact(pairs);

    let mut text = 0usize;
    for i in 0..pairs {
        // The header matches mk_clump_fixture's, so text bytes per record match
        // the calibration matrix at the same read length and record count.
        let (a, b) = if spine_only {
            (FastqRecord { id: String::new(), seq: String::new(), qual: String::new() },
             FastqRecord { id: String::new(), seq: String::new(), qual: String::new() })
        } else {
            (FastqRecord {
                 id: format!("@SIM:1:FCX:1:15:6329:{i} 1:N:0:ATCCGA"),
                 seq: seq.clone(),
                 qual: qual.clone(),
             },
             FastqRecord {
                 id: format!("@SIM:1:FCX:1:15:6329:{i} 2:N:0:ATCCGA"),
                 seq: seq.clone(),
                 qual: qual.clone(),
             })
        };
        text += estimated_record_bytes(&a) + estimated_record_bytes(&b);
        r1.push(a);
        r2.push(b);
        keys.push((i % 4096) as MinimizerKey);
    }

    // Read every structure after the peak so nothing above can be elided.
    let guard = r1.len() + r2.len() + keys.len()
        + r1[pairs / 2].id.len() + r1[pairs / 2].seq.len()
        + r2[pairs / 2].qual.len() + keys[pairs / 2] as usize;
    let mode = if spine_only { "spine" } else { "full" };
    println!("PAIRS={pairs} READ_LEN={read_len} TEXT_BYTES={text} MODE={mode} GUARD={guard}");
}
