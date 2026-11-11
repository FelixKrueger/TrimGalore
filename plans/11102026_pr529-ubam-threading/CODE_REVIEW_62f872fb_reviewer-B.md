# Code review — PR #529 "feat(ubam): multithreaded uBAM input and output with `--cores`"

**Reviewer B** (independent; no coordination with Reviewer A)
**Head commit:** `62f872fb` · **Base:** `dev` @ `c0bbee4` · **Date:** 2026-11-10
**Source read:** pinned extract at `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/aecc6671-429c-462c-ad9c-38f8e4d35fbe/scratchpad/pr529_B`
**Base read:** `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/aecc6671-429c-462c-ad9c-38f8e4d35fbe/scratchpad/base_c0bbee4`

> No source file was edited. This report recommends only.

---

## Summary

The design is sound and the parts I was asked to scrutinise hardest hold up. I verified by reading
that:

- **Ordering is correct.** `spawn_conversion_pool` (reader side) and the two
  `run_*_parallel_to_bam` merge loops (writer side) both use a monotonic `seq` plus a `BTreeMap`
  reorder buffer keyed on `expected`. Crucially, ordering never depends on batch *contents*: a
  batch that yields zero records after filtering still occupies its `seq` and advances `expected`
  (`src/parallel.rs:419-427`, `:522-530`), and `read_raw_pairs` / `read_raw_records` never emit an
  empty batch (`src/bam.rs:159-165` and `:170-172`; `:210-216` and `:217-219`), which matters because `Ok(Vec::new())` is the
  reader channel's EOF sentinel — an empty batch would read as clean EOF.
- **Phred-offset parity across `--cores` holds.** `build_record_buf` is the single construction
  site; nothing hardcodes 33 on the write path any more. The serial writer passes
  `cli.phred_offset()`, the parallel encoder passes `config.phred_offset`, and
  `src/main.rs:1350` sets `phred_offset: cli.phred_offset()` on the one `TrimConfig` that reaches
  every `run_ubam_output_*` driver. Same value, both paths. The `batch_encoder_subtracts_input_phred_offset`
  test is genuinely discriminating (see "Test discrimination").
- **Header is written exactly once**, before any record block, in `ParallelBamSink::create`
  (`src/bam.rs:1098-1103`), and `BatchEncoder` never calls `write_header`, so worker batches are
  record bytes only.
- **Mate adjacency and FREAD1/FREAD2 survive.** Pairs are formed on the single reading thread
  (`read_raw_pairs`) and travel as a `RawPair`, so a worker can never split or reorder a pair;
  `process_paired_batch_bam` pushes `Some(1)` then `Some(2)` per surviving pair.
- **Backpressure is bounded and I found no deadlock.** Every channel is a `sync_channel`: 2 per
  conversion worker, `workers * 2` on the merge channel, `BAM_READER_CHANNEL_BATCHES = 4` per
  reader output, 2 per trim worker, `cores * 2` on the result channel, and noodles' own
  `bounded(worker_count)` pair inside the BGZF pool
  (`noodles-bgzf-0.35.0/src/multithreaded_writer/builder.rs:59-60`). The reorder `BTreeMap`s are
  therefore bounded by those capacities rather than unbounded. The paired `deliver`
  (`src/bam.rs:459`) sends R1 then R2 on separate 4-deep channels, which is safe because the
  consumer (`read_pairs_round_robin`) alternates one record per side — a lockstep consumer can
  never demand R2 while an R1 send is parked. Its `&&` short-circuit (skipping the R2 send when
  the R1 send fails) is also correct: a dropped receiver means teardown, and returning `false`
  stops the merger. Error propagation out of a dying worker I traced too, and it surfaces: the
  worker sends `Err` before breaking, the main loop records it as `first_error`, and the reader's
  subsequent send failure makes it `break` to `Ok(())` so `reader_join_err.or(first_error)` yields
  the real error rather than the placeholder.
- **De-interleaver behaviour is preserved.** `MAX_SLACK`, the grouped-input error, the
  invalid-flag-bits error and the orphan error all match the base (`src/bam.rs:110-182` vs base
  `src/bam.rs:212-322`), and the final partial batch is emitted *before* the orphan check bails, so
  nothing in front of the error is lost.

and I then **confirmed the central invariants by running the binary** — byte-identical record
streams across `--cores 1/2/3/4/8` on 4M single-end and 1M paired records, and correct `--phred64`
QUAL arithmetic end to end. See "Verified by execution" for the data, including the first
same-build `--cores 1` vs `--cores N` performance numbers for this PR (4–5.4× wall for ~11–20%
more CPU).

**One Critical finding**, and it is not a threading bug: `ParallelBamSink` writes straight to the
final output path with `File::create`, bypassing the `crate::io::PendingOutput` partial-output
guard that **every other record writer in the repo** uses. I reproduced both halves: a failed
`--cores 4` run publishes a *structurally valid, `samtools quickcheck`-clean, zero-record* BAM at
the final name, and it destroys a good output already sitting at that path. The serial path leaves
nothing behind and leaves the old file untouched. CI cannot see this because the uBAM tests only
exercise success paths.

Beyond that: a mid-stream write error aborts the process with `SIGABRT` and loses the real error
(reproduced); the thread count has tripled against a documented and nf-core-facing "N+4" contract;
plus a swallowed tail-flush error, an unwired `--compression`, and three stale documentation
statements with no CHANGELOG entry.

---

## Critical

### C-1 — `ParallelBamSink` bypasses the `PendingOutput` guard: a failed run publishes a valid-looking EMPTY BAM over the previous one

**`src/bam.rs:1091`** (`ParallelBamSink::create`), against **`src/bam.rs:912`** (`BamWriter::create`).

```rust
// ParallelBamSink::create — src/bam.rs:1091
let file = File::create(path)
    .with_context(|| format!("Failed to create uBAM output: {}", path.display()))?;
```

versus the serial writer, and **every other record writer in the repo**:

| writer | guard |
|---|---|
| `BamWriter::create` (serial uBAM) | `crate::io::PendingOutput::create` — `src/bam.rs:912` |
| `parallel.rs` FASTQ outputs | `PendingOutput::create` — `src/parallel.rs:220`, `:1300` |
| `FastqWriter` | `PendingOutput::create` — `src/fastq.rs:731` |
| `clump_only` SE + PE | `PendingOutput::create` — `src/clump_only.rs:343`, `:497-498` |
| **`ParallelBamSink::create`** | **none — `File::create` on the final path** |

`PendingOutput`'s contract (`src/io.rs:867-872`) is explicit: *"Dropped without committing, the
temporary is removed and the final path is never created — so a refused run leaves no output, and a
previous run's output at that path survives untouched."* The new sink breaks both halves.

It is worse than a plain truncated file. `noodles::bgzf::MultithreadedWriter` has a `Drop` impl
(`noodles-bgzf-0.35.0/src/multithreaded_writer.rs:149-158`) that calls `finish()` and therefore
**appends the BGZF EOF block**. So the abandoned file is a *structurally valid, complete, empty*
BAM — `samtools quickcheck`, which exists precisely to detect truncation by looking for that EOF
block, passes it.

**Verified by execution** (debug build of `62f872fb`, FASTQ of 6,000 reads where read 5,000 has a
301-byte QNAME that BAM must reject):

```
$ trim_galore --output-format ubam --cores 1 -o serial/ bad.fq
exit=1
serial/                       <- empty; no output, no .partial
Error: while writing BAM record "@XXXX…(301 bytes)" — BAM rejected it. …

$ trim_galore --output-format ubam --cores 4 -o par/ bad.fq
exit=1
-rw-r--r--  265  par/bad_trimmed.bam      <- published at the FINAL name
  file bytes: 265 · ends with BGZF EOF block: True · magic: b'BAM\x01' · 0 records
```

And the clobber half, seeding a known-good 10-record BAM at the output path first:

```
before: 507 bytes, md5 5ea4da8b67288e09feed5f2987d3406c   (10 records)
$ trim_galore --output-format ubam --cores 4 -o clobber/ bad.fq   # exit=1
after:  268 bytes, md5 9fce758d3e0c87acc1f2eae56ba9db5f   (0 records, valid EOF)
```

**Failure scenarios this covers** (any error after `create`, which is all of them):
over-long or illegal QNAME; a B/H aux tag in the input; an aligned record in a uBAM input; a
mate-sync error; grouped paired-uBAM input; `ENOSPC`; `SIGKILL`/OOM/`Ctrl-C` (here not even the
error text appears — just a valid empty-or-partial BAM at the final name).

**Why CI is green anyway:** `tests/integration_ubam_out.rs` exercises the `--cores 4` path only on
the success path, and the CI `validation-ubam` job likewise. No test asserts the *absence* of an
output file after a refused uBAM run on the parallel path — `plans/08132026_partial-output-on-refusal`
suggests that invariant was established deliberately, and nothing pins it here.

A second, quieter aggravation: the justification comment the reader reuses —
`src/parallel.rs:993`, `:1035` and `:1502`, *"Safe to drop: a dead receiver means that worker
already errored, and publication is gated on the worker join"* — is now attached to a path where
**publication is not gated on anything**. (The error itself does still surface correctly; I traced
that. It is only the output-file half of the premise that is false.)

**Recommendation:** give `ParallelBamSink` the same `PendingOutput` treatment as `BamWriter`: hold
the guard in the struct, `File::create` via `PendingOutput::create`, and `commit()` inside
`finish()` after `self.inner.finish()`. Note `MultithreadedWriter: Drop` means the temporary will
be a well-formed empty BAM when abandoned — that is fine, because `PendingOutput: Drop` deletes it.

**Confidence: certain** (reproduced; both halves).

---

## High

### H-1 — A mid-stream write error on the parallel sink aborts the process (SIGABRT) and the real error is never reported

**`src/bam.rs:1111-1113`** (`ParallelBamSink::write_encoded`), via
`noodles-bgzf-0.35.0/src/multithreaded_writer.rs:139`.

`MultithreadedWriter`'s internal `send()` does `write_tx.send(buffered_rx).unwrap()`. The writer
thread (`spawn_writer`) exits on the *first* `io::Error` from `write_frame`, dropping `write_rx`.
The next `write_all` on the main thread therefore panics on that `unwrap()`. Then, while
unwinding, `Drop for MultithreadedWriter` calls `finish()` → `flush()` → the same `send()` →
**the same `unwrap()` panics again**, which is a panic in a destructor during cleanup: Rust turns
that into a non-unwinding panic and aborts.

**Verified by execution** (debug build of `62f872fb`; the output path pre-created as a FIFO whose
reader exits after 1,000 bytes, giving a clean `EPIPE` once the pipe buffer fills — the same shape
as `ENOSPC` or `EIO`):

```
$ trim_galore --output-format ubam --cores 4 -o fifo/ big.fq     # 200,000 reads
thread 'main' panicked at …/noodles-bgzf-0.35.0/src/multithreaded_writer.rs:139:36:
called `Result::unwrap()` on an `Err` value: "SendError(..)"
…
panicked at core/src/panicking.rs:233:5:
panic in a destructor during cleanup
thread caused non-unwinding panic. aborting.
exit=134
```

The I/O error that caused all this — the thing the user needs to see, "No space left on device" —
is **discarded**: `spawn_writer`'s `io::Result<W>` is only recoverable via `finish()`, which never
gets to run cleanly. The user gets a Rust backtrace and `SIGABRT`.

By contrast the serial path (`bam::io::Writer<bgzf::Writer<BufWriter<File>>>`) has no `unwrap()`
anywhere in its write chain: every `write_record` returns `io::Result`, which `BamWriter` wraps
via `write_error_context` (`src/bam.rs:966`), and `finish()` adds
`.context("failed to finalise uBAM output …")`. So this is a regression in failure behaviour
introduced by the choice of writer, not pre-existing parity. (Verified by reading; I could not
construct an equivalent serial-path I/O failure to compare side by side, because `PendingOutput`
redirects the serial writer onto a regular `.partial` file.)

**Why CI misses it:** no test induces an output-side I/O error on any path.

**Recommendation:** wrap the writes so the panic cannot reach the user, or avoid
`MultithreadedWriter`. The cheapest containment that keeps the design:
`std::panic::catch_unwind` is not usable here (the second panic is in `Drop`), so the practical
options are (a) upstream a non-panicking `send` in noodles, or (b) use `bgzf::Writer` for framing
and keep parallelism in the encode stage only — which is where this PR already puts most of the
win, since `BatchEncoder` does the record building. Worth a note to the contributor either way:
this needs a decision, not a tweak. Also consider `ENOSPC` explicitly — it is the realistic
trigger on a cluster.

**Confidence: certain** on the behaviour (reproduced); **high** on the mechanism (noodles source
read).

### H-2 — Thread count is now 3N–4N, contradicting the documented and nf-core-facing "N+4" contract

**`src/parallel.rs:79`** (`/// Total threads: N + 4 (or N + 5 with --passthrough)`), and
`docs/src/content/docs/performance/threading.md`.

Counted from the code (`spawn_conversion_pool` at `src/bam.rs:241-262`; the scoped pools at
`src/parallel.rs:376` and `:481`; `ParallelBamSink::create` at `src/bam.rs:1093-1096`):

| configuration | threads at `--cores N` | at `--cores 8` |
|---|---|---|
| FASTQ in → FASTQ out (unchanged) | N + 4 | 12 |
| uBAM in → FASTQ out, SE | 2N + 4 | 20 |
| uBAM in → FASTQ out, PE two files | 3N + 6 | 30 |
| uBAM in → uBAM out, SE | 3N + 5 | 29 |
| uBAM in → uBAM out, PE two files | 4N + 7 | 39 |

(conversion pool = N workers + 1 reader + 1 merger, *per BAM input*; trim pool = N workers +
1 reader; BGZF sink = N deflaters + 1 writer; + main.)

`docs/…/performance/threading.md` is emphatic and quantitative about the opposite:

- page description: *"a fixed +4 thread infrastructure cost"*
- line 35: *"From `--cores 2` upward the **N+4 model** applies — each added core is exactly 1 worker thread"*
- line 46: *"At Perl `-j 8` vs v2.x `--cores 8`: up to ~27 vs **exactly 12 threads**"*
- line 89: *"Trim Galore v2.x can use the **full CPU allocation directly** … If a Nextflow process has 12 CPUs, just pass `--cores 12`."*
- line 39 has a comparison table whose *Perl* column is `~3N+3` — which is where the uBAM paths now sit.

So a BAM-input nf-core process following the documented advice (`--cores task.cpus`) will
oversubscribe its allocation 3–4×. That is not a cosmetic docs issue: it is the guidance the
threading page exists to give.

I am **not** arguing the pools are wrong — overlapping BAM decode with trimming is the right
idea, and the conversion pool is where the uBAM-input win comes from. The issue is that the
thread budget tripled silently. Options, in the order I'd consider them: cap the conversion pool
and the BGZF pool (e.g. `cores.div_ceil(2)` or a small constant) so the total stays near 2N;
or keep the pools and update `src/parallel.rs:79` plus the threading page with a per-format table
and revised nf-core advice. Either way the page cannot stay as it is.

**Confidence: high** on the arithmetic (read from code; I could not count threads at runtime —
`ps` is blocked in this sandbox — so these are calculated, not measured).

---

## Verified by execution

I built `62f872fb` in my own target dir and ran the binary. These are measurements, not estimates.

### The central ordering invariant holds

**uBAM in → uBAM out, single-end, 4,000,000 records (977 batches):**

| `--cores` | records | md5 of the record region (header excluded) |
|---|---|---|
| 1 | 4,000,000 | `b2642d4b39d1dc87044e67efdf693cfa` |
| 2 | 4,000,000 | `b2642d4b39d1dc87044e67efdf693cfa` |
| 4 | 4,000,000 | `b2642d4b39d1dc87044e67efdf693cfa` |
| 8 | 4,000,000 | `b2642d4b39d1dc87044e67efdf693cfa` |

**Paired interleaved uBAM in → uBAM out, 999,948 pairs in, 186 pairs dropped by filters** (so the
drop path is exercised, and `--cores 3` is deliberately a core count that does not divide the
batch count):

| `--cores` | records | md5 of the record region | first six flags |
|---|---|---|---|
| 1 | 1,999,524 | `03415e007afff3b15727c51f1c8c66bc` | `4d 8d 4d 8d 4d 8d` |
| 3 | 1,999,524 | `03415e007afff3b15727c51f1c8c66bc` | `4d 8d 4d 8d 4d 8d` |
| 8 | 1,999,524 | `03415e007afff3b15727c51f1c8c66bc` | `4d 8d 4d 8d 4d 8d` |

Byte-identical record streams, correct mate adjacency, correct FREAD1 (`0x4d`) / FREAD2 (`0x8d`)
alternation. The header differs only because `@PG CL:` embeds the `-o` path.

### Issue #358 is NOT reintroduced on the parallel path

The `config.phred_offset` → `BatchEncoder` wiring is untested in the suite (see M-5), so I tested
it. 20,000 reads with QUAL `'h'` (ASCII 104) under `--phred64 --output-format ubam`:

| `--cores` | distinct stored QUAL bytes |
|---|---|
| 1 | `{40: 680000}` |
| 4 | `{40: 680000}` |

`104 − 64 = 40` on both paths. Had `BatchEncoder` kept the hardcoded 33, `--cores 4` would have
stored 71. **Correct.**

### Performance (same build, same machine, M-series arm64, release)

The contributor supplied no numbers; these are mine.

**uBAM in → uBAM out, single-end, 4M records:**

| `--cores` | wall | speedup | CPU (user+sys) |
|---|---|---|---|
| 1 | 4.52 s | 1.00× | 4.03 s |
| 2 | 1.50 s | 3.02× | 4.54 s |
| 4 | 0.91 s | 4.96× | 4.83 s |
| 8 | 0.83 s | 5.41× | 4.81 s |

**Paired interleaved uBAM in → uBAM out, 1M pairs:**

| `--cores` | wall | speedup | CPU |
|---|---|---|---|
| 1 | 5.66 s | 1.00× | 6.43 s |
| 3 | 1.87 s | 3.03× | 6.98 s |
| 8 | 1.40 s | 4.04× | 7.11 s |

**The win is real: 4–5.4× wall for ~11–20% more CPU.** Two things worth saying about its shape,
because they affect how it should be described:

1. **The 1 → 2 step is superlinear (3.0×).** That is not scaling; it is the serial path having
   left BGZF *decode* and *encode* completely unoverlapped. Most of the headline number is
   pipelining, not parallelism — which is why it is largely spent by `--cores 2`.
2. **It plateaus at `--cores 4`** (4.96× → 5.41× for double the cores, SE). Combined with H-2's
   thread arithmetic, `--cores 4` buys ~92% of the available win at 17 threads where `--cores 8`
   spends 29 for the last 8%. Whatever is decided about H-2, the recommended setting for uBAM
   paths is lower than for FASTQ paths, and the threading page should say so.

I could **not** measure peak RSS: `ps` and `/usr/bin/time -l` are both blocked in this sandbox
(`sysctl kern.clockrate: Operation not permitted`). I am therefore reporting no memory number —
see M-2 for what the code implies, clearly marked as unmeasured.

---

## Medium

### M-1 — `ParallelBamSink::finish` throws away the writer it is handed, so a tail-flush error is silently swallowed

**`src/bam.rs:1115-1121`**

```rust
pub fn finish(mut self) -> Result<()> {
    self.inner
        .finish()                       // -> io::Result<BufWriter<File>>
        .context("failed to finalise uBAM output (BGZF EOF marker)")?;
    Ok(())                              // the BufWriter is dropped here, unflushed
}
```

`MultithreadedWriter::finish()` *returns the inner writer* (`noodles-bgzf-0.35.0/src/multithreaded_writer.rs:93`,
`writer_handle.join().unwrap()` yields `io::Result<W>`). The `?` produces that `BufWriter<File>`
and the `;` drops it immediately. `Drop for BufWriter` flushes and **discards the error**.

Concretely: the BGZF EOF block is 28 bytes and the preceding partial frame is under `BufWriter`'s
8 KiB default, so on `ENOSPC` at the very end the last frame and the EOF marker are lost, no error
is raised, and the run **reports success** with a BAM that `samtools quickcheck` will reject.

The serial `BamWriter::finish` (`src/bam.rs:948-958`) has the same shape (`try_finish()` then
`drop(self)` then `commit()`), so this is not a regression — but the fix on the new path is one
line and worth taking while the code is open:

```rust
let mut w = self.inner.finish().context("…")?;
std::io::Write::flush(&mut w).context("failed to flush uBAM output")?;
```

**Confidence: high** (noodles source + signature read; not reproduced — H-1's abort fires first
for mid-stream errors, and I could not provoke an error confined to the final flush).

### M-2 — Worker batches carry *uncompressed* BAM, so in-flight memory is materially higher than on the FASTQ path

**`src/parallel.rs:559`, `:587`** (`BatchEncoder::with_capacity(reads.len() * 400, …)`),
**`src/parallel.rs:341-347`, `:452-456`** (`BamBatchResult.bytes`, `SingleBamBatchResult.bytes`).

On the FASTQ path a worker gzips its batch, so what travels the result channel and sits in the
reorder `BTreeMap` is *compressed*. On the new uBAM path the worker emits uncompressed BAM and the
sink compresses later, so the same buffers hold roughly 4× the bytes — and the uBAM path adds a
second pool in front of it (the conversion pool, holding raw `bam::Record`s *and* converted
`FastqRecord`s).

Bounds read from the code at `--cores N`, paired, uBAM in → uBAM out:

- conversion pool: N×2 raw batches in the per-worker channels + N×2 on the merge channel + up to
  the same again in its reorder `BTreeMap` + 4 per side on the two output channels;
- trim pool: N×2 pair-batches queued + N in flight;
- result stage: N×2 `BamBatchResult` on the channel + up to N×2 more pending in the `BTreeMap`,
  each preallocating `4096 × 2 × 400 B ≈ 3.3 MB`.

That is tens of batches of multi-MB buffers rather than the FASTQ path's compressed equivalents.
**I did not measure this** — `ps` and `/usr/bin/time -l` are blocked in this sandbox — so I am
deliberately not quoting an RSS figure. What I can say from the code is that the structure
guarantees a step up, and that `docs/…/performance/threading.md:56-63` ("under ~100 MB across all
reasonable core counts", "~10 MB per additional worker") is a FASTQ-path statement that should not
be read as covering this path.

Two cheap mitigations if a measurement confirms it matters: bound the pending `BTreeMap` (it is
already bounded indirectly by the channel capacities, but only loosely), and size
`BatchEncoder::with_capacity` from the batch's actual bytes rather than a flat 400 B/record — see
L-4.

**Confidence: medium-high** on the mechanism, **unmeasured** on the magnitude.

### M-3 — `--compression` is silently ignored on the uBAM path, and that is the biggest remaining throughput lever

**`src/bam.rs:1093-1096`**

```rust
let mut inner = bgzf::multithreaded_writer::Builder::default()
    .set_worker_count(workers)
    .build_from_writer(BufWriter::new(file));
```

`Builder::default()` is `CompressionLevel::default()` = **6**
(`noodles-bgzf-0.35.0/src/writer/compression_level.rs:93-97`). TrimGalore's FASTQ default is
**1** (`src/fastq.rs:20-27`, `DEFAULT_GZIP_LEVEL`), chosen deliberately because compression CPU
dominated wall time, with `--compression` as the counter-lever. So uBAM output pays roughly 6×
the deflate CPU per byte of FASTQ output, and `--cores`-scaling it does not change that — it just
buys more threads to spend on it.

This is pre-existing (the serial `BamWriter` goes through `bgzf::Writer::new`, also level 6), so
the two paths at least agree. But the PR is specifically a throughput PR, and it introduces the
one `Builder` call site where honouring `cli.compression` is a single `.set_compression_level(…)`.
Given the measured plateau at `--cores 4`, level is probably the larger remaining win. Changing
the *default* would change output bytes, so the right move is to wire the flag and leave the
default alone — and say in the docs that `--compression` now applies to uBAM output.

**Confidence: certain** on the levels (both read from source); **suspected** on the size of the win
(not measured — I did not patch the level to test it, since I am not editing source).

### M-4 — The only end-to-end `--cores` + uBAM coverage runs a *single batch*

**`tests/integration_ubam_out.rs:311`** and the CI `validation-ubam` job.

`ubam_out_cores_gt_1_produces_output` runs `--cores 4` on `test_files/ubam_test.bam`, which I
counted: **10 records**. `BAM_READER_BATCH_SIZE` and `BATCH_SIZE` are both 4,096, so that is one
batch, handed to one worker. The reorder `BTreeMap`, the round-robin across the other three
workers, the `expected` advance, the final flush and the BGZF pool under load are all
**unexercised at the binary level**. The CI job uses the same 10-record fixture.

The multi-batch coverage that does exist is good but lives entirely in-crate
(`src/bam.rs:1789` (`pooled_readers_match_plain_readers`) at 10,000 records, `src/parallel.rs:2811` and `:2894` at 10,000 records). What is
missing is a test that drives the *binary*, through `main.rs`'s dispatch, across a batch boundary.

Also untested anywhere: **`ParallelBamSink::create(path, Some(&source_header), …)`**. Both parity
tests pass `None`, so the uBAM-in header/`@PG`-chain propagation through the *parallel* sink has
no coverage — exactly the surface the brief flags, since `assert_ubam_eq` ignores `@PG`.

The repo already contains `test_files/ubam_ws_qname_bulk.bam` — **5,002 records**, which crosses
one batch boundary. The concrete asks:

1. an integration test running `--cores 1` and `--cores 4` on `ubam_ws_qname_bulk.bam` with
   `--output-format ubam` and comparing the record streams (the md5-of-record-region comparison in
   "Verified by execution" is exactly this, and it passes today — it just is not pinned);
2. a parity test with `source_header = Some(…)` so the propagated header is compared between the
   serial writer and `ParallelBamSink`.

**Confidence: certain** (record counts measured; test bodies read).

### M-5 — The uBAM parity tests cannot detect a phred-offset wiring defect

**`src/parallel.rs:1708` (inside `parity_test_config`, `src/parallel.rs:1702`)** — `parity_test_config()` sets `phred_offset: 33`.

Both `test_parallel_ubam_output_matches_serial` and
`test_parallel_single_end_ubam_output_matches_serial` build their config from
`parity_test_config()`, so serial and parallel both pass 33, which is indistinguishable from a
hardcoded 33. `batch_encoder_subtracts_input_phred_offset` (`src/bam.rs:2552`) *is* genuinely
discriminating — `stored_qual("hhhh", 64)` expecting `[40; 4]` fails loudly if the encoder
hardcodes 33 — but it calls `BatchEncoder` directly, so it proves the encoder honours the offset
it is *given*, not that `process_*_batch_bam` passes the right one.

The wiring is in fact correct (I traced `src/main.rs:1350` → `config.phred_offset` →
`src/parallel.rs:559`/`:587`, **and confirmed it by execution** — see "Verified by execution").
This is a test-strength finding, not a bug: the one defect class this PR is most at risk of
reintroducing (#358) currently has no guard on the path it was reintroduced into. A
`--phred64 --output-format ubam --cores 4` integration test, or simply `phred_offset: 64` in one
of the two parity configs, closes it.

**Confidence: certain** (read; and the correctness half independently verified by running it).

---

## Low

### L-1 — Two comments the PR introduces break the project's comment rules

**`src/bam.rs:1421-1424`** — four lines where the rule is one (two with a reason), carrying both
evidence and historical narrative:

```rust
// One pass over the record's aux fields (looking each tag up separately
// re-scans the fields once per requested tag, which made this the
// dominant cost of reading a uBAM). The first occurrence of a tag wins,
// as with `Data::get`.
```

"looking each tag up separately re-scans the fields once per requested tag" describes the code
being replaced, and "which made this the dominant cost" is the measurement. The keeper is the
invariant — `// The first occurrence of a tag wins, as with `Data::get`.` — which constrains the
loop below it. The rest belongs in the commit message.

**`tests/integration_ubam_out.rs:312-313`** — textbook "reference the previous state you are
fixing":

```rust
// `--cores N>1` + uBAM-out runs on the worker pool (it used to warn that
// --cores was ignored) and produces the output BAM.
```

The rename from `ubam_out_cores_gt_1_warns_and_proceeds` to `ubam_out_cores_gt_1_produces_output`
already carries this; the parenthetical is noise to a reader who never saw the old assertion.

*Noted separately, not this PR's fault:* `build_record_buf` (`src/bam.rs:1011-1019`) inherits the
long #358 narrative block ("Hardcoding 33 here was issue #358: Phred+64 input was under-subtracted
by 31 …") verbatim from `BamWriter::write_record`. It is relocated, not introduced, so it is out of
scope here — but the move was the natural moment to shorten it, and now that two call sites share
the function the block is the only remaining place the rule is broken on this path.

### L-2 — Three user-facing statements are now false, and there is no CHANGELOG entry

The PR updates `CLAUDE.md` and `README.md` but leaves:

- **`docs/src/content/docs/guide/outputs.md:65`** — *"uBAM output is currently single-threaded —
  `--cores N` is ignored for BAM writing; FastQC (when `--fastqc` is requested) still uses the value."*
- **`wiki/overview.md:73`** — *"uBAM output is single-threaded in v1 and refuses seven flags …"*
- **`wiki/usage/ubam.md:80`** — *"Threading | Single-threaded compression in v1. `--cores` is ignored
  for BAM writing and an up-front note says so; FastQC still uses it."* (the "up-front note" no
  longer exists — the PR deletes it from `src/main.rs`)
- **`src/parallel.rs:79`** — *"Total threads: N + 4 (or N + 5 with `--passthrough`)"* — see H-2.
- **`CHANGELOG.md`** — no entry. Per the project's rule this is one sentence, e.g.
  `- [#529](url) - uBAM input and output now use `--cores` ([@Justype](url))`.

`CHANGELOG.md:888-889` and `docs/…/reference/changelog.md:894-895` also say "single-threaded in
v1", but those are the historical v2.1 record and are correct to leave alone.

### L-3 — The `cli.phred_offset()` invariant comment no longer enumerates all its consumers

**`src/main.rs:425-432`** declares *"INVARIANT ESTABLISHED HERE, RELIED ON FAR BELOW"* and then
lists who relies on it: *"Every `cli.phred_offset()` call that feeds `BamWriter::create` …"*.
`BatchEncoder` is now a second consumer, reached via `config.phred_offset` rather than a direct
`cli.phred_offset()` call — so a reader grepping for `BamWriter::create` will not find it. The
invariant still holds (verified); the comment that exists to protect it has gone stale. One clause
added to that list is the whole fix.

### L-4 — `BatchEncoder::with_capacity` sizes from a flat 400 bytes per record

**`src/parallel.rs:559`, `:587`** — `reads.len() * 400`, i.e. 1.6 MB (SE) / 3.3 MB (PE) per batch
regardless of read length. For 50–100 bp reads that over-allocates ~3–4×; for a long-read uBAM
(ONT/PacBio, which `samtools import` happily produces) it under-allocates by an order of magnitude
and the `Vec` reallocates repeatedly while multi-MB. Sizing from the batch's actual
`seq.len() + qual.len() + id.len()` — the crate already has `clump::estimated_record_bytes` for
exactly this — would fix both ends. Minor on its own; it compounds M-2.

### L-5 — Error precedence differs slightly from the serial de-interleaver

**`src/bam.rs:110-182`** vs base `src/bam.rs:212-322`. The base converts each record (which runs
the per-record `is_unmapped()` check) *before* queueing it, so an aligned record is rejected the
moment it is read. `read_raw_pairs` queues raw records and converts in the workers, so an aligned
record that never finds a mate is reported as the **orphan** error ("ended with N unmatched R1 …")
rather than the aligned-BAM error, and on grouped input the `MAX_SLACK` error can pre-empt an
aligned-record error that the base would have reported first. Both still fail; only the message
changes, and only on inputs that are broken in two ways at once. Cosmetic — noted so it is not
mistaken for a regression later.

### L-6 — `assert_eq!(n, 4_096, …)` silently hardcodes `BAM_READER_BATCH_SIZE`

**`src/bam.rs:1887`**. The assertion is one of the strongest in the PR (it is what proves the
merger delivers batch 0 before the error at record 4,500), but it breaks confusingly if
`BAM_READER_BATCH_SIZE` ever changes. `assert_eq!(n, BAM_READER_BATCH_SIZE, …)` reads the same and
survives.

### L-7 — Asymmetric handling of the final `emit` return value

**`src/bam.rs:217-219`** — `read_raw_records` calls `emit(batch)` for the final partial batch and
drops the `bool`, where `read_raw_pairs` (`src/bam.rs:170-172`) honours it. Harmless (nothing
follows either call but `Ok(())`), but the asymmetry invites a reader to wonder which is correct.

### L-8 — A panicking conversion worker degrades to clean EOF rather than an error

**`src/bam.rs:241-315`** (`spawn_conversion_pool`) with **`src/bam.rs:843`** (`Err(_) => Ok(None)`).

If a worker thread *panics*, the batch it held never reaches the merger. The merger blocks on
`expected`, never sees `PoolMsg::End`, and exits when `merge_rx` closes — **without calling
`deliver` at all**. The output channel's sender drops, and `BamReader::next_record`'s
`Err(_) => Ok(None)` reads a closed channel as clean EOF. Result: truncated output, exit 0.
`_handle: spawn_noop_joinhandle()` means the panic is never joined, so nothing surfaces it.

Narrow — only a panic triggers it, and `bam_record_to_fastq` returns `Result` for everything it
anticipates. The base has the same `Err(_) => Ok(None)` shape with a single producer thread; the
PR widens the surface from 1 thread to N+2 per BAM input. Worth a `catch_unwind` in the worker
body, or (cheaper) making the merger send `PoolMsg::Fail` if `merge_rx` closes before `End`
arrives.

**Confidence: high** on the mechanism (traced); **not reproduced** (I did not inject a panic).

---

## Test discrimination

The brief asked whether the new tests can actually fail. Assessed individually:

| test | discriminating? | note |
|---|---|---|
| `pooled_readers_match_plain_readers` (`src/bam.rs:1789`) | **Yes** | 10,000 records = 2.4 batches, unique `read_{i}` ids, compares `(id, seq, qual)` tuples — a reorder changes the tuple sequence and fails. |
| `pooled_readers_report_the_same_errors`, grouped half (`:1812`) | **Yes** | 2,000 R1 then 2,000 R2 genuinely exceeds `MAX_SLACK = 1024`. |
| same, aligned-record half | **Yes, strongest in the PR** | `assert_eq!(n, 4_096)` fails at `n == 0` if the merger emitted the error ahead of batch 0. Pins ordering *and* error placement in one assertion. |
| `batch_encoder_subtracts_input_phred_offset` (`:2552`) | **Yes, for the encoder** | `("hhhh", 64) → [40;4]` fails loudly on a hardcoded 33. Does not cover the caller — see M-5. |
| `test_parallel_ubam_output_matches_serial` (`src/parallel.rs:2811`) | **Mostly** | Byte-compares decompressed streams incl. header; fixture deliberately drops pairs. Blind to phred-offset wiring (M-5) and to header propagation (`source_header = None`, M-4). |
| `test_parallel_single_end_ubam_output_matches_serial` (`:2894`) | **Mostly** | Same strengths and same two blind spots. |
| `ubam_out_cores_gt_1_produces_output` (`tests/integration_ubam_out.rs:311`) | **No** | See below. |

**The one that cannot fail usefully** is `ubam_out_cores_gt_1_produces_output`. Its assertions are
(a) exit success, (b) `!stderr.contains("ignored")`, (c) the output file exists. (b) is a negative
assertion that passes for essentially any stderr — including one from a run that never touched the
worker pool — and it is the *only* assertion that gestures at the feature. Combined with M-4's
10-record fixture (one batch, one worker), this test would pass unchanged if the parallel path were
reverted to serial tomorrow, provided the NOTE stayed deleted. And per C-1, assertion (c) — "the
output BAM exists" — is now satisfied by a *failed* run too.

The fix is to assert something positive: run the same input at `--cores 1` and `--cores 4` and
compare the record streams (the comparison in "Verified by execution" is exactly that, and it
passes today).

I found **no case** where a new assertion points the same direction as the bug it guards — the
failure mode the brief warned about. The weakness here is a vacuous assertion, not an inverted one.

---

## Recommendations

| # | Finding | Priority | Action |
|---|---|---|---|
| C-1 | `ParallelBamSink` bypasses `PendingOutput`; a failed run publishes a valid empty BAM over the previous one | **Critical** | Blocking. Route `ParallelBamSink::create` through `PendingOutput::create` and `commit()` in `finish()`. Add a test asserting no output file after a refused `--cores N` uBAM run. |
| H-1 | Mid-stream write error → `SIGABRT`, real error lost | **High** | Blocking in my view: `ENOSPC` on a cluster is routine and this turns it into an abort with no diagnosis. Needs a design decision (upstream fix, or keep BGZF framing single-threaded and parallelise only the encode). |
| H-2 | Thread count 3N–4N vs documented N+4 | **High** | Either cap the conversion/BGZF pools or update `src/parallel.rs:79` and the threading page (incl. the nf-core `--cores task.cpus` advice). Given the measured plateau at `--cores 4`, capping is attractive. |
| M-1 | `finish()` drops the writer unflushed | Medium | One-line fix on the new path; mention the serial path has the same shape. |
| M-2 | In-flight memory step up (uncompressed batches + a second pool) | Medium | Measure peak RSS at `--cores 8` before deciding; I could not. Then size `with_capacity` from real bytes (L-4). |
| M-3 | `--compression` ignored on the uBAM path | Medium | Wire `cli.compression` into the new `Builder` call site; leave the default at 6 so output bytes do not move. |
| M-4 | Only end-to-end `--cores` uBAM test runs one batch; header propagation untested | Medium | Use the existing `ubam_ws_qname_bulk.bam` (5,002 records); add a parity test with `source_header = Some(…)`. |
| M-5 | Parity tests blind to phred-offset wiring | Medium | `phred_offset: 64` in one parity config, or a `--phred64 --cores 4` integration test. |
| L-1 | Two comment-rule violations | Low | Cut to the invariant line; drop the "it used to warn" parenthetical. |
| L-2 | Three false docs statements, no CHANGELOG entry | Low | Blocking for merge hygiene, trivial to do. |
| L-3 | Invariant comment does not name its new consumer | Low | One clause. |
| L-4 | Flat 400 B/record capacity | Low | Use `clump::estimated_record_bytes`. |
| L-5 | Error-precedence drift in the pooled de-interleaver | Low | Accept; note it. |
| L-6 | `4_096` hardcoded in an assertion | Low | Use the constant. |
| L-7 | Asymmetric `emit` return handling | Low | Make both honour it. |
| L-8 | Panicking conversion worker → clean EOF | Low | Make the merger emit `Fail` when `merge_rx` closes before `End`. |

### Overall

The threading design is correct and the win is real and measured (4–5.4× wall for ~11–20% more
CPU). Ordering, mate adjacency, flag bits and phred-offset parity all hold under execution at
several core counts, including a core count that does not divide the batch count — I could not
break them.

What stops this being a clean approval is not the concurrency: it is that the new output path
dropped two guarantees the rest of the codebase keeps. **C-1** (no `PendingOutput`) and **H-1**
(abort instead of an error on `ENOSPC`) both make a *failed* uBAM run worse than it is today, and
both are invisible to a test suite that only exercises success. C-1 in particular turns a refused
run into a published, `quickcheck`-clean, empty BAM that overwrites good data — the worst shape a
bioinformatics failure can take, because nothing downstream notices.

C-1 and L-2 are cheap. H-1 and H-2 need the maintainer to make a call, not the contributor.
