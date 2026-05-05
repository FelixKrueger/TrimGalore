---
title: Threading model
description: Why the Oxidized Edition is faster. Single-pass paired-end, worker-pool parallelism, and a fixed +4 thread infrastructure cost.
---

## Why the speedup is larger than "Rust is faster"

**The bottleneck is gzip, not trimming.** The actual trimming logic (adapter alignment, quality clipping) is ~5% of runtime; the other 95% is gzip compression (~60%) and decompression (~30%). Rust's speed advantage over Perl/Python only applies to that 5%.

The real wins come from **architectural differences**.

## 1. Single-pass vs three-pass

Trim Galore runs Cutadapt on R1, then R2, then pair-validates: reading and recompressing the data three separate times. The Oxidized Edition does everything in one pass.

## 2. Worker-pool parallelism

Each worker independently handles trimming **and** gzip compression for its batch of reads, producing independently-compressed gzip blocks concatenated in order (valid per RFC 1952). This distributes the dominant cost (compression) across N workers instead of funneling through one thread.

## 3. Fewer threads, more work per thread

The Oxidized Edition uses a single process with a fixed infrastructure cost of +4 threads:

```
Oxidized Edition --cores N thread breakdown:
  N worker threads (each: trim + gzip compress -> independent gzip block)
  2 decompression threads (one per input file)
  1 batcher thread (creates numbered batches of 4096 reads)
  1 main thread (collects blocks in order -> writes to output files)
                                                       = N+4 total
```

Legacy Perl Trim Galore (v0.6.x) instead orchestrated three subprocesses (Cutadapt + pigz compress + pigz/igzip decompress) that each spawned their own threads, peaking around `3N+3` and forcing the nf-core module to reserve `task.cpus - 4` cores for the `-j` flag. See the [migration guide](/reference/migration/#whats-new-in-v2) for the full v0.6.x → v2 contrast.

At `--cores 1`, the worker-pool is bypassed entirely: a single thread does everything with zero parallelism overhead (1 thread, 5 MB RAM). From `--cores 2` upward the **N+4 model** applies — each added core is exactly 1 worker thread plus ~10 MB of memory on top of the 4 fixed-cost threads (2 decompressors + 1 batcher + 1 writer). Wall-clock speedup is **near-linear up to `--cores 8` for paired-end runs**; beyond that, gzip-output I/O on the storage layer typically becomes binding before workers run out of useful per-read work, so each additional core helps progressively less.

For reference, the contrast against the legacy Perl 0.6.x model at the same `-j N` / `--cores N` setting:

| Cores | Perl 0.6.x threads (~3N+3) | Oxidized threads (N+4) |
|------:|---------------------------:|-----------------------:|
| 1 | up to ~6 | 1 |
| 4 | up to ~15 | 8 |
| 8 | up to ~27 | 12 |
| 16 | not measured | 20 |

At Perl `-j 8` vs Oxidized `--cores 8`: up to ~27 vs exactly 12 threads, yet **4.54× faster** wall and **5.93× less CPU** on the 84M-read Buckberry fixture (v2.1.0-beta.7).

Parallel efficiency at Buckberry scale: 100% (cores=1) → 72% (cores=8) → 34% (cores=16) → 22% (cores=24). Scaling is near-linear up to `--cores 8`; beyond that, gzip-output I/O on the storage layer typically becomes binding before workers run out of useful per-read work, so adding cores helps progressively less. **`--cores 8` is the sweet spot for nf-core / Snakemake / CWL workflows** — also the saturation point.

## Memory profile

Peak resident set size on the 84M-read Buckberry fixture (Trim Galore v2.1.0-beta.7, measured via `/usr/bin/time -v`):

| `--cores` | Peak RSS | Notes |
|----------:|---------:|-------|
| 1 | 5.1 MB | Worker-pool bypassed; single-threaded path. |
| 2 | 60.1 MB | Infrastructure threads come online (+4 fixed cost: 2 decompressors + 1 batcher + 1 writer). |
| 4 | 73.0 MB | |
| 8 | 91.6 MB | |

The c1 → c2 step is by far the largest (+55 MB) — that's the four infrastructure threads spinning up their I/O buffers. Past c2, each additional pair of workers adds ~7 MB on average; the bulk of which is the per-worker compression buffer. Memory growth is bounded and predictable, well-suited to cluster scheduling: even worst-case at the saturation point (`--cores 8`), the process never exceeds ~100 MB.

## Beyond speed

- **Zero external dependencies:** No Python, no Cutadapt, no pigz. Single static binary.
- **Simpler deployment:** `cargo install` or download a binary. No conda environment needed.
- **Single-pass paired-end:** Both reads processed together, with guaranteed synchronization, no temp files.
- **Lower memory:** 5 MB single-threaded, ~10 MB per additional worker. No Python interpreter, no subprocess pipes.
- **CPU-efficient:** Uses 5.9× to 13.5× less CPU time than Trim Galore (nf-core default to single-thread, on the 84M-read Buckberry fixture). Meaningful on shared HPC clusters where CPU-hours = money.
- **Reproducible:** Pure Rust with deterministic behaviour across platforms.
- **New features:** Poly-G trimming (auto-detected for 2-colour instruments like NovaSeq/NextSeq) and poly-A trimming, both built in without external tools.

## What this means for nf-core / Nextflow

The Oxidized Edition can use the **full CPU allocation directly**: no need to subtract cores for subprocess overhead, since everything runs in a single process. If a Nextflow process has 12 CPUs, just pass `--cores 12`.

For the historical nf-core pattern of `task.cpus - 4`, the equivalent Oxidized invocation is `--cores task.cpus`, with the fixed `+4` thread cost matching the existing CPU budget without manual subtraction.

## Reproducibility

`--cores N` produces byte-identical decompressed output for any N (verified via md5 across the benchmark range). The worker-pool emits independently-compressed gzip blocks in deterministic order, so the gzipped bytes themselves vary by core count, but decompressing them yields the same FASTQ content every time.
