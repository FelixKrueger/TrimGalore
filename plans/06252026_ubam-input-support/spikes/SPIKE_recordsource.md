# SPIKE — RecordSource dispatch overhead (uBAM input support)

## Question

The uBAM-input-support plan must dispatch between two reader types
(`FastqReader`, `BamReader`) at the `main.rs` entry point. TrimGalore today
references the concrete `FastqReader` at ~10 call sites in `src/parallel.rs`.
Two options are on the table:

- **(B) Trait object:** `Box<dyn RecordSource> { fn next_record(&mut self) -> Result<Option<FastqRecord>>; }`
  — dynamic dispatch, smallest diff to `parallel.rs`.
- **(C) Enum + match:** `enum Reader { Fastq(FastqReader), Bam(BamReader) }` with a hand-rolled match
  in `next_record` — static dispatch, slightly larger diff but predictable codegen.

**Does dynamic dispatch in the inner record-read loop measurably regress wall-clock at production scale?**

## Success criteria

- **Trait wins:** ≤2% slower than the concrete baseline → recommend (B).
- **Enum wins:** trait is >2% slower than concrete → recommend (C).

## Scope boundary

- Out of scope: a real BAM parser. The question is dispatch overhead; the
  parser body is identical across variants. We mimic the TrimGalore
  FastqReader parser shape (`BufReader::read_line` into `String`, three
  `String`s per record — mirrors `src/fastq.rs:322-355`) and dispatch only
  through the three different shells.
- Out of scope: gzip-wrapped input, paired-end coordination, channel
  hand-off. The dispatch overhead is what we're measuring, and the
  baseline doesn't change with those layers wrapped around it.
- Out of scope: comparing against `paraseq` or other parsers — that's
  `plans/paraseq-spike`.

## Strategy

Build a standalone Cargo project under
`plans/06252026_ubam-input-support/spikes/spike1-recordsource/` (does NOT
modify the parent TrimGalore crate). Three variants over the same 500K×150bp
fixture, each computing the same sentinel checksum:

- **(A)** Concrete `FastqLikeReader` struct, no abstraction (baseline).
- **(B)** `Box<dyn RecordSource>` trait object (vtable dispatch).
- **(C)** `enum Reader { Fastq(...), Bam(...) }` + match (jump-table dispatch).

All three call the same `parse_one_record(&mut dyn BufRead, &mut String)`
function so the parser body is literally identical instructions — the only
thing that varies is how `next_record()` is reached. The enum variant has a
dead `Bam` arm wired in so the compiler must still emit the discriminant
check on each call (mirroring production code).

Release profile mirrors TrimGalore's `Cargo.toml`:
`opt-level = 3, lto = true, codegen-units = 1, strip = "debuginfo"`.

Each variant runs 3× plus a single pre-bench warm-up of (A) to prime the OS
page cache. Reports median / min / max wall-clock ms.

## Script

- Project: `plans/06252026_ubam-input-support/spikes/spike1-recordsource/`
- Entrypoint: `src/main.rs`
- Build: `cargo build --release`
- Fixture: `./target/release/spike1_recordsource build-fixture 500000 /tmp/claude/spike1_500k.fastq`
  (165 MB plain FASTQ; LFSR-generated, deterministic, same pattern as
  `plans/paraseq-spike/spikes/spike_paraseq.rs`)
- Bench: `./target/release/spike1_recordsource bench /tmp/claude/spike1_500k.fastq`

## Results

### Iteration 1

```
(A) concrete baseline             median=  86.93 ms   min=  86.90 ms   max=  88.07 ms
(B) Box<dyn RecordSource>         median=  88.04 ms   min=  87.99 ms   max=  89.74 ms
(C) enum Reader + match           median=  87.45 ms   min=  87.12 ms   max=  88.47 ms

Δ vs. concrete baseline (median):
  (B) trait : +1.28%   (+1.12 ms)
  (C) enum  : +0.60%   (+0.52 ms)
```

Sentinel match: all three variants → `reads=500000  seq_sum=75000000  sentinel=35787052`.
Parsers are doing identical work.

### Iteration 2 (stability confirmation, same binary)

```
(A) concrete baseline             median=  85.63 ms
(B) Box<dyn RecordSource>         median=  86.37 ms
(C) enum Reader + match           median=  86.65 ms

Δ vs. concrete baseline (median):
  (B) trait : +0.86%   (+0.74 ms)
  (C) enum  : +1.19%   (+1.02 ms)
```

In iteration 2 (B) is faster than (C). The trait/enum ordering flips between
runs — dispatch is below the noise floor.

### Meets criteria?

Yes — (B) trait dispatch is **+0.86% to +1.28%** vs. baseline, comfortably
under the 2% bar. Trait wins.

## Findings

1. **Dispatch overhead is below noise** at 500K reads (165 MB input). The
   max trait penalty across two runs was +1.28% (1.12 ms over 87 ms); the min
   was +0.86%. The trait-vs-enum ordering inverts between iterations, which
   means the difference between (B) and (C) is statistically indistinguishable
   — both add ≈1 ms over an 86–87 ms baseline.

2. **The inner-loop body dominates.** Each record costs four `BufRead::read_line`
   calls plus three `String` allocations (≈300 bytes alloc per record). One
   indirect call per record vs. one match arm vs. one direct call is invisible
   next to that workload. The vtable lookup hits L1 cache on the very first
   record and stays there — there's no realistic memory-bandwidth penalty.

3. **Extrapolation to Buckberry scale (≈84M reads):** linearly, at +0.86% the
   trait penalty is ~120 ms over ~14 s of pure-parse work. In practice the
   trim pipeline adds adapter alignment, quality trimming, gzip compression,
   and channel I/O — all of which dwarf 120 ms.

## Reference snippets

The dispatch boundary itself (the bit that an implementation agent would
otherwise re-derive):

```rust
// The trait — single method, owned-record return type. Send is needed
// because parallel.rs ships the reader handle into the worker pool.
trait RecordSource: Send {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

// Construction site at main.rs would look like:
let mut source: Box<dyn RecordSource + Send> = match input_kind(&path)? {
    InputKind::Fastq => Box::new(FastqReader::open(&path)?),
    InputKind::Bam   => Box::new(BamReader::open(&path)?),
};
```

The `parse_one_record` helper splitting trick is what made the variants
truly comparable — all three call into the same monomorphized function
body, so any measured delta is purely dispatch:

```rust
#[inline(always)]
fn parse_one_record(reader: &mut dyn BufRead, line_buf: &mut String)
    -> Result<Option<FastqRecord>> { ... }
```

(Note `inline(always)` here — the spike forces inlining so the trait
variant's dispatch happens once at the `next_record` call site, not also at
the parser body. In production code we don't need this; the natural inlining
is fine.)

## Recommendation

**Go with (B) — `Box<dyn RecordSource>`.** Trait-object dispatch costs
≈1 ms per 500K reads (under 1.3% of wall-clock, well within the 2%
success-criterion bar) and the trait gives `parallel.rs` the smallest
possible diff: the ~10 concrete `FastqReader` call sites only need their
type narrowed to `&mut dyn RecordSource`, not rewritten around a match.

The enum option (C) is a defensible fallback if the plan ever needs to
expose reader-specific methods (e.g., BAM tag inspection) without touching
the trait surface, but for the dispatch-overhead question that motivated
this spike, the two are statistically indistinguishable. Don't pay the
plan-diff cost without a separate reason.

## Limitations

1. **Plain FASTQ only.** Gzip-wrapped input adds ≈10–15× wall-clock per
   record (decoder time). The dispatch overhead would be an even smaller
   percentage of total time in the gzip path — i.e., this measurement is
   the *worst case* for dispatch overhead. Real BAM input adds htslib
   decompression on top of dispatch.

2. **No multi-threaded measurement.** The TrimGalore worker-pool path
   (`src/parallel.rs`) ships records across an `mpsc` channel. The channel
   hand-off cost (heap-allocated `FastqRecord` move + `Sender::send`) is
   already orders of magnitude above per-record dispatch, so threading
   doesn't change the conclusion — but I didn't measure it.

3. **macOS / Apple Silicon only.** The exact ms numbers won't repro on
   Linux/x86, but the dispatch-vs-noise ratio should — trait-object dispatch
   is a single indirect call regardless of platform, and 500K iterations is
   plenty for the dispatch cost to surface if it were going to.

4. **Single fixture.** 500K × 150bp. Real-world reads vary length-by-length;
   shorter reads would push the per-record dispatch fraction up, but only
   trivially — the four `read_line` calls and three `String` allocs survive
   regardless of seq length.
