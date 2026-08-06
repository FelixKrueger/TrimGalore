# SPIKE — rapidgzip-core v0.2.0: does non-seekable (streaming) gzip input work?

**Date:** 2026-08-05 · Darwin 25.5.0 / macOS 26.5.1, arm64 · rustc 1.95.0
**Prompted by:** Rob Patro on [TrimGalore#379](https://github.com/FelixKrueger/TrimGalore/issues/379#issuecomment-5181968638) — *"This should be fixed in the latest main. I will publish 0.2.0 today but it would be great if you could test main."*
**Under test:** `rapidgzip-core` at git tag **v0.2.0** = commit `f648e3bb`, which **is** `main`'s HEAD, so this tests both at once.
**Upstream issue:** [COMBINE-lab/rapidgzip-rust#6](https://github.com/COMBINE-lab/rapidgzip-rust/issues/6) — "Design: support non-seekable (streaming) compressed input", filed by @BenjaminDEMAILLE, closed COMPLETED.

## 1. Question, criteria, scope

**Question.** Does v0.2.0 decode gzip arriving on a non-seekable source — a pipe, a FIFO, `/dev/stdin`, a process substitution — correctly, for the container shapes Trim Galore cares about: plain single-member gzip, concatenated multi-member gzip, and BGZF?

**Success criteria.** Decoded bytes **byte-identical** to the original uncompressed data for every (format × source) pair; no hang; an explicit error rather than silent truncation when the input is damaged.

**Out of scope.** Throughput (settled at ~1.04x for Trim Galore by two earlier spikes — decompression is ~3.7% of realistic runtime and the serial stage is the output writer). Integration into Trim Galore. The write side, which a decoder-only library has nothing to do with.

## 2. Script

- `spikes/spike_nonseekable/` — a small cargo crate; `src/main.rs` decodes to **stdout** so the driver can `cmp`, and reports telemetry on **stderr**.
- `spikes/spike_nonseekable_driver.sh` — builds fixtures, runs the matrix, byte-compares.

```
cd plans/rapidgzip-trial/spikes/spike_nonseekable && cargo build --release
bash ../spike_nonseekable_driver.sh ./target/release/spike-nonseekable <scratch-dir>
```

Every invocation is bounded by `perl -e 'alarm 60; exec @ARGV'` — `alarm(2)` survives `execve`, and there is no `timeout(1)` on Darwin. A FIFO read that blocked would otherwise hang indefinitely.

## 3. Results

**11 passed, 0 failed.** Fixtures built from one 1 888 894-byte FASTQ; every row compares the decoder's output against that file with `cmp`.

| Case | Source | Format | Result |
|---|---|---|---|
| `regular_single` | regular file, `open()` | gzip, 1 member | byte-identical |
| `regular_bgzf` | regular file, `open()` | BGZF, 31 members | byte-identical |
| `pipe_single` | `cat … \|`, `stream_reader` | gzip | byte-identical |
| `pipe_multi` | `cat … \|`, `stream_reader` | **concatenated**, 2 members | byte-identical |
| `pipe_bgzf` | `cat … \|`, `stream_reader` | **BGZF** + 28-byte EOF member, 31 members | byte-identical |
| `fifo_open` | named FIFO, `open()` | gzip | byte-identical |
| `open_/dev/stdin_from_file` | `open /dev/stdin < f` | gzip | byte-identical |
| `open_/dev/stdin_from_pipe` | `cat f \| open /dev/stdin` | gzip | byte-identical |
| `open_process_substitution` | `open <(cat f)` | gzip | byte-identical |

Negative controls, both of which fired:

| Control | Result |
|---|---|
| Truncated gzip on a stream | **errors**: `invalid DEFLATE data at bit 320000: truncated DEFLATE stream` — not silent truncation |
| `cmp` can detect a difference | a 9-byte sentinel append is caught, so the nine PASSes above mean something |

**The answer to Rob's question is yes**, for every shape we care about, including the two that matter most to a bioinformatics consumer and are easy to get wrong: concatenated members and BGZF's trailing 28-byte EOF member.

## 4. Findings beyond the pass/fail

**`Decoder::open` routes by probing the opened handle.** `config.rs:774`:

```rust
pub fn open<P: AsRef<Path>>(&self, path: P) -> Result<DecoderReader, DecodeError> {
    let file = File::open(path).map_err(|error| DecodeError::input_io(0, error))?;
    if supports_positional_reads(&file) { self.reader(file) } else { self.stream_reader(file) }
}
```

So a consumer handing over a path needs no special case — which is the API shape that makes this usable. Note this is the same discriminator shape Trim Galore just shipped for #379: interrogate the *opened handle*, not the filename, and not the path.

**Non-seekable is sequential by design, not parallel.** The README is explicit: *"the four parallel decode paths all need positional reads, so a non-seekable source always uses the sequential path."* This is the load-bearing fact for anyone hoping streaming support brings parallel decode with it. It does not.

**Nothing is spooled.** Also from the README: input memory is one `input_page_size` window, nothing to memory or disk. That answers one of the open questions carried in our notes — memory cost per stream against a `--memory` budget — favourably.

**Dependency closure is 6 crates, all pure Rust:** `zlib-rs`, `libz-rs-sys`, `crossbeam-{utils,epoch,deque}`. No C zlib arrives by the back door, which is the property Trim Galore pins `flate2` to `zlib-rs` to preserve.

## 5. Iteration log

1. Built and ran the matrix first try. 11/11 byte-identical. But every row reported `spawned_workers=0` — including the regular-file rows that were supposed to be the parallel control.
2. Hypothesised the fixture was too small. Rebuilt at 57 MB uncompressed. Same zeros, and `active_workers=0` too — which *contradicts* the documented contract of exactly one active worker on a stream. That made the telemetry, not the library, the suspect.
3. Moved the sample from after-EOF to mid-stream (past 8 MB decoded). Now `mid_active=1` — matching the documented stream contract. But `mid_active=1, spawned=0` on the **regular file** as well, so the parallel path never engaged for either.

## 6. Limitations — read these before quoting anything above

- **The parallel path was never exercised, so the regular-file rows are not a control.** All three sources reported `mid_active=1, mid_spawned=0`, and 57 MB decoded in 0.03–0.05 s from both a file and a pipe with no measurable difference. The likeliest explanation is the fixture: 30 concatenated copies of the same text compresses to 1.7 MB, plausibly under whatever chunk threshold the parallel path needs. **This is a limitation of the experiment, not a defect claim about the library** — and confirming which it is was outside this spike's question.
- Consequently this spike says nothing about whether routing produces a *behavioural contrast*, only that both routes decode correctly.
- Darwin only. `/dev/fd/N` is a `dup` here and resolves through `/proc/self/fd/N` on Linux, so `open_/dev/stdin_from_file` in particular may take a different route there. Untested.
- One fixture shape (short fixed-length FASTQ, highly compressible). No real-world sequencing file, no gzip produced by other implementations, no zlib or raw-DEFLATE containers, no socket source.
- Timings are recorded only to show the file/pipe *comparison*; they are not a throughput measurement and should not be quoted as one.

## 7. Recommendation

**Report success to Rob**, with the parallel-path caveat stated as ours rather than his.

**Do not adopt on throughput grounds.** Unchanged from the two earlier spikes: decompression is ~3.7% of realistic runtime, so the ceiling is ~1.04x, and the serial stage is the output writer — which a decoder-only library cannot touch. Nothing here revisits that.

**The capability question is now open where it was closed before.** Streaming decode is a genuine prerequisite for ever supporting `<(zcat …)` in Trim Galore. It is *not* sufficient: #379's blocker is that we hand a **path** to five separate readers (format detection, sanity check, adapter auto-detection, the trim pass), so real pipe support needs a reader threaded through the call graph — options B and C in `plans/08032026_reject-non-restartable-input/PLAN.md` §10, which were declined on that cost, not on decompression. v0.2.0 removes one obstacle, not the binding one.

If pipe support is ever wanted, the sequencing is: thread a reader through the front end first (B/C), and only then consider whether the streaming decoder is the right thing to put behind it.

---

# Part 2 — the meaningful run, on oxy

The laptop run answered the correctness question but its fixture was worthless for anything else: 30 copies of one text compressed to 1.7 MB, so the parallel path never engaged and "0.04 s either way" said nothing. Re-run on `dockyard-oxy-0` with real data.

**Host:** `dockyard-oxy-0`, 64 physical / 128 logical cores, single socket, cargo/rustc 1.96.0.
**Scripts:** `spike_nonseekable_oxy.sh` (correctness matrix), `spike_nonseekable_oxy_timing.sh` (timings), `spike_nonseekable_oxy_multimember.sh` (real multi-member).

**Fixtures.** `~/benchmark_TG_oxy/real_files/` from the earlier spikes is **gone** and `/home` is 90% full, so the canonical file was re-downloaded from ENA — deliberately the same accession the two prior spikes used, so results are comparable:

| Fixture | Compressed | Decoded | Members |
|---|---|---|---|
| `SRR24827373_1.fastq.gz` (real ENA WGBS) | 2.13 GB | **16 577 221 356 B** | 1 |
| `synth_bclconvert_R1.fastq.gz` (different producer) | 0.22 GB | — | 1 |
| `synth_barcode_R1_val_1.fq.gz` (Trim Galore's own output) | 0.27 GB | — | **1** |
| `mm.gz` — 3 GB of real ENA reads, split and re-gzipped | 0.47 GB | 2.79 GB | **3** |

## Correctness: 17 of 17, all byte-identical

Every case compares an md5 of the decoder's output against `gzip -dc | md5sum` — an independent reference decoder, not a self-comparison. All three fixtures pass through **both** `open()` on a regular file and `stream_reader` on a pipe, as does the 3-member archive. A truncated 100 MB prefix of the real file errors (`invalid DEFLATE data at bit 800000000: truncated DEFLATE stream`), and the md5 control confirms the comparison can distinguish.

## The routed contrast, which Part 1 could not show

| Source | threads asked | `mid_active` | `mid_spawned` |
|---|---|---|---|
| regular file | 1 | 1 | 0 |
| regular file | 4 | 4 | **4** |
| regular file | 16 | 16 | **16** |
| regular file | 32 | 16 | 16 |
| regular file | 64 | 16 | 16 |
| **pipe** | 1 | 1 | 0 |
| **pipe** | 32 | **1** | **0** |

The parallel path engages on real data, and a non-seekable source stays sequential no matter how many threads are requested — exactly the documented contract, now measured rather than read. Worker count **caps at 16** however many are asked for; not investigated, and not obviously a problem.

## Timings — 2.13 GB compressed → 16.58 GB decoded, two reps each

| Path | rep1 | rep2 |
|---|---|---|
| file, 1 thread (sequential) | 10 208 ms | 10 311 ms |
| file, 4 threads | 8 659 ms | 8 187 ms |
| **file, 16 threads (parallel)** | **4 697 ms** | **4 701 ms** |
| file, 32 threads | 4 775 ms | 4 605 ms |
| **pipe, 32 threads (streaming)** | **10 973 ms** | **10 878 ms** |
| `gzip -dc` | 55 547 ms | 55 641 ms |
| `cat` (I/O floor) | 409 ms | 298 ms |

Four things fall out:

1. **Streaming costs ~2.3x against parallel** (10.9 s vs 4.7 s) on the same bytes. That is the price of non-seekable input, and it is a real cost rather than a rounding error.
2. **Even sequential, rapidgzip is 5.4x faster than system `gzip`** (10.2 s vs 55.6 s) — the zlib-rs backend, not the parallelism. Streaming is still 5.1x faster than `gzip -dc`.
3. **Parallelism buys 2.2x**, not 16x, at 16 workers on this file. Plateaus by 16.
4. `cat` at ~0.35 s means this is entirely CPU-bound; none of the above is I/O.

## What this does and does not mean for Trim Galore

**It does not change the ~1.04x conclusion, and there is a trap here worth naming.** A reader could see "parallel decode saves 5.5 s on this file" and conclude Trim Galore would get 5.5 s faster. It would not. Trim Galore decompresses on its **own thread**, overlapped with parallel trimming, so `T = max(T_reader, T_trim/N)`: the ~10 s of decode CPU is largely hidden behind compute, which is exactly why the earlier spike measured decompression at **3.7% of wall-clock** rather than the ~25% its CPU share would suggest. Those two numbers are consistent, not contradictory — and confusing them is the easiest way to misread this whole spike.

**Our current reader is already in rapidgzip's sequential class.** Trim Galore's `flate2` is pinned to `zlib-rs`, the same backend, single-threaded — so ~10 s for this file is roughly what we already do. Adopting the *streaming* path would be performance-neutral for us; adopting the *parallel* path would speed up a stage that is already overlapped.

## Limitations of Part 2

- One host, one accession, one read length. Linux/x86_64 only here; Part 1's Darwin rows cover the other platform, and `/dev/fd/N` semantics differ between them.
- The worker cap at 16 is unexplained. It may be a chunking consequence of this file, or an internal ceiling; either way it means "32 threads" and "64 threads" rows are really the 16-worker configuration.
- Timings are two reps on a shared host, not a controlled benchmark, and `taskset` pinning was not used. They are sound enough for the ~2x and ~5x conclusions and should not be quoted to two significant figures.
- BGZF was only tested at laptop scale (31 members, 1.9 MB). The oxy multi-member fixture is plain concatenated gzip, not BGZF.
- Trim Galore's own output turned out to be single-member (`--cores 1` produces one member), so "does rapidgzip decode our concatenated-member output" is answered by the constructed 3-member archive rather than by a genuine Trim Galore parallel-writer file. Close enough to call the question settled; not identical.

## Note on plan linkage

`plans/rapidgzip-trial/` has no `PLAN.md` — only `PROGRESS.md` and this `spikes/` directory — so there is no plan section to append to, per the spike skill's plan-linked branch. The two prior spikes in this line are `spikes/SPIKE_reader_bound.md` here and `plans/reader-serial-stages/spikes/SPIKE_serial_stages_oxy.md`.
