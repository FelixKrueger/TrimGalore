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

When you pass `-j N` to Trim Galore, three separate programs each independently spawn threads. The theoretical maximum is approximately 3N+3, though not all threads are necessarily active simultaneously:

```
Trim Galore -j N thread breakdown (theoretical maximum):
  Cutadapt:                    N workers + 1 reader + 1 writer  =  N+2
  pigz (compress):             N threads                         =  N
  pigz/igzip (decompress):     up to N threads                   ≈  N
  Perl:                        1 main process                    =  1
                                                                 ≈ 3N+3 total
```

Note: The nf-core trimgalore module accounts for this by reserving `task.cpus - 4` cores for the `-j` flag (e.g., 12 allocated CPUs to `-j 8`). Thread counts above were observed via `ps` during benchmarking and represent approximate peak values.

The Oxidized Edition uses a single process with a fixed infrastructure cost of +4 threads:

```
Oxidized Edition --cores N thread breakdown:
  N worker threads (each: trim + gzip compress -> independent gzip block)
  2 decompression threads (one per input file)
  1 batcher thread (creates numbered batches of 4096 reads)
  1 main thread (collects blocks in order -> writes to output files)
                                                       = N+4 total
```

At `--cores 1`, the worker-pool is bypassed entirely: a single thread does everything with zero parallelism overhead (1 thread, 5 MB RAM). The infrastructure cost only applies from `--cores 2` upward, where each additional core adds exactly 1 thread and ~10 MB of memory.

| Cores | TG threads (up to ~3N+3) | Oxidized threads (N+4) |
|------:|-------------------------:|-----------------------:|
| 1 | up to ~6 | 1 |
| 4 | up to ~15 | 8 |
| 8 | up to ~27 | 12 |
| 16 | not measured | 20 |

At `-j 8` vs `--cores 8`: up to ~27 vs exactly 12 threads, yet 1.9x faster.

Parallel efficiency on the Xeon: 82% (2 cores) to 82% (4 cores) to 81% (8 cores) to 78% (16 cores) to 64% (24 cores). Scaling remains near-linear up to 16 cores, with diminishing returns beyond that. For most production use, `--cores 8` to `--cores 16` is the sweet spot. Beyond 16, additional cores still help but deliver progressively less benefit per core.

## Memory profile

| `--cores` | Memory | Notes |
|----------:|-------:|-------|
| 1 | 5 MB | Worker-pool bypassed; single-threaded path. |
| 2 | 43 MB | Infrastructure threads come online (+4 fixed cost). |
| 4 | 62 MB | |
| 8 | 100 MB | |
| 16 | 171 MB | |
| 24 | 157 MB | |

Each additional worker adds ~10 MB on average. The bulk of which is the per-worker compression buffer.

## Beyond speed

- **Zero external dependencies:** No Python, no Cutadapt, no pigz. Single static binary.
- **Simpler deployment:** `cargo install` or download a binary. No conda environment needed.
- **Single-pass paired-end:** Both reads processed together, with guaranteed synchronization, no temp files.
- **Lower memory:** 5 MB single-threaded, ~10 MB per additional worker. No Python interpreter, no subprocess pipes.
- **CPU-efficient:** Uses 2.6 to 5x less CPU time than Trim Galore. Meaningful on shared HPC clusters where CPU-hours = money.
- **Reproducible:** Pure Rust with deterministic behaviour across platforms.
- **New features:** Poly-G trimming (auto-detected for 2-colour instruments like NovaSeq/NextSeq) and poly-A trimming, both built in without external tools.

## What this means for nf-core / Nextflow

The Oxidized Edition can use the **full CPU allocation directly**: no need to subtract cores for subprocess overhead, since everything runs in a single process. If a Nextflow process has 12 CPUs, just pass `--cores 12`.

For the historical nf-core pattern of `task.cpus - 4`, the equivalent Oxidized invocation is `--cores task.cpus`, with the fixed `+4` thread cost matching the existing CPU budget without manual subtraction.

## Reproducibility

`--cores N` produces byte-identical decompressed output for any N (verified via md5 across the benchmark range). The worker-pool emits independently-compressed gzip blocks in deterministic order, so the gzipped bytes themselves vary by core count, but decompressing them yields the same FASTQ content every time.
