# Optimus Prime — Project Summary

## What was built

A complete **Rust rewrite of Trim Galore** that produces **byte-identical output** to the Perl original across every feature and test case. It's a true drop-in replacement — same CLI flags, same output filenames, same report format compatible with MultiQC.

**Architecture shift:** TrimGalore (Perl) is a wrapper that shells out to Cutadapt (Python/Cython) for adapter matching. Optimus Prime does everything in a single process — adapter detection, alignment, quality trimming, adapter removal, filtering — in one pass through the data. Paired-end reads are processed in a single pass rather than two sequential Cutadapt runs.

## Feature parity

| Feature | Status |
|---------|--------|
| Adapter auto-detection (Illumina/Nextera/smallRNA/BGI/Stranded) | Done |
| Quality trimming (BWA algorithm, Phred33/64) | Done |
| Adapter trimming (semi-global alignment with error rate) | Done |
| Paired-end single-pass processing | Done |
| `--rrbs` / `--non_directional` | Done |
| `--nextseq` / `--2colour` | Done |
| `--consider_already_trimmed` | Done |
| `--poly_a` trimming | Done |
| `--hardtrim5` / `--hardtrim3` | Done |
| `--clock` (Epigenetic Clock UMI) | Done |
| `--implicon` (UMI from R2) | Done |
| `--demux` (3' barcode demultiplexing) | Done |
| `--fastqc` / `--fastqc_args` | Done |
| `--rename`, `--trim-n`, `--max_n`, `--max_length` | Done |
| `--retain_unpaired`, `--clip_R1/R2`, `--three_prime_clip_R1/R2` | Done |
| `--cores N` (worker-pool parallelism) | Done (new) |
| Trimming reports (MultiQC-compatible) | Done |
| Colorspace rejection | Done |

## Performance: the benchmark story

**Test dataset:** 55.8M paired-end reads (SRR24827378), typical whole-genome bisulfite sequencing.
**Platform:** macOS, Apple M2 Pro (10 cores).

### Optimization progression (single-threaded TrimGalore as baseline)

| Configuration | Wall time | Speedup |
|---|---|---|
| TrimGalore `-j 1` (Perl + Cutadapt) | 27:04 | 1.0x |
| Optimus Prime, miniz_oxide, cores=1 | 20:06 | **1.35x** |
| + zlib-rs backend | 12:01 | **2.26x** |
| + pipeline parallelism (`--cores 2`) | 6:29 | **4.17x** |
| **Worker-pool (`--cores 2`)** | 7:16 | **3.73x** |
| **Worker-pool (`--cores 4`)** | **3:46** | **7.19x** |

### Fair comparison: multi-core TrimGalore vs Optimus Prime

| Configuration | Wall time | CPU time | Speedup (wall) |
|---|---|---|---|
| TrimGalore `-j 1` | 27:04 (1,624s) | 2,536s | 1.0x |
| TrimGalore `-j 2` | 13:50 (830s) | 2,741s | 2.0x |
| TrimGalore `-j 4` | 7:15 (435s) | 2,868s | 3.7x |
| OP `--cores 1` (sequential) | 11:46 (706s) | 695s | 2.3x |
| OP `--cores 2` (worker-pool) | 7:16 (436s) | 908s | 3.7x |
| **OP `--cores 4` (worker-pool)** | **3:46 (226s)** | **936s** | **7.2x** |

With 4 workers, Optimus Prime is **7.2x faster than TrimGalore `-j 1`** and **1.9x faster than TrimGalore `-j 4`**. CPU efficiency remains excellent: OP `--cores 4` uses 936 CPU-seconds (4.1 cores busy) while TrimGalore `-j 4` burns 2,868 CPU-seconds (6.6 cores busy) — a **3.1x difference in compute cost**. On shared clusters where CPU-hours = money, this matters.

## Why we initially expected 5-10x

The reasoning was sound in theory:
- **Rust vs Perl:** Rust is typically 50-100x faster than Perl for CPU-bound string processing
- **Single-pass vs two-pass:** TrimGalore runs Cutadapt separately for R1 and R2. Optimus Prime processes both in a single interleaved pass, halving file I/O
- **No subprocess overhead:** TrimGalore spawns Cutadapt as a child process, communicating via pipes and temp files. Optimus Prime does everything in-process
- **Cutadapt's Python overhead:** Even with Cython, there's interpreter overhead per read

A 5-10x improvement would be realistic if the trimming logic were the bottleneck.

## How we broke through the 4x ceiling

**The bottleneck is gzip, not trimming.** The actual trimming logic (adapter alignment, quality clipping) is ~5% of runtime — the other 95% is gzip compression (~60%) and decompression (~30%). Rust's speed advantage over Perl/Python only applies to that 5%.

**Phase 1: Pipeline model (capped at ~4x)**

The initial approach used a pipeline: reader threads → single main thread → writer threads. This hit a ceiling at `--cores 2` (389s wall-time) because the single main thread was funneling all data through one point. Adding more compression threads didn't help — the main thread couldn't feed them fast enough.

**Phase 2: Worker-pool model (linear scaling)**

The breakthrough was rethinking the architecture. Instead of a pipeline, each worker independently handles both trimming **and** gzip compression for its batch of reads. The workers produce independently-compressed gzip blocks, which are concatenated in order — valid per RFC 1952 (gzip members concatenate freely).

```
Architecture (--cores N):
  2 decompression threads (background, one per input file)
  1 batcher thread (creates numbered batches of 4096 reads)
  N worker threads (each: trim + gzip compress → independent gzip block)
  1 main thread (collects blocks in order → writes to output files)
```

This scales because the dominant cost (gzip compression) is now distributed across N workers instead of funneled through one thread:
- **cores=1**: 706s wall (sequential, no worker-pool overhead)
- **cores=2**: 436s wall (1.62x speedup, 81% efficiency)
- **cores=4**: 226s wall (3.13x speedup, 78% efficiency)

The ~78% parallel efficiency is limited by irreducible sequential work: input decompression, batch distribution, and output file writing all run on single threads (Amdahl's law). CPU time increases from 695s → 936s because independent gzip blocks are slightly less efficient than one continuous stream (smaller compression context per block).

## What Optimus Prime does better regardless of speed

- **Zero external dependencies:** No Python, no Cutadapt, no pigz. Single static binary.
- **Simpler deployment:** `cargo install` or download a binary. No conda environment needed.
- **Single-pass paired-end:** Both reads processed together — guaranteed synchronization, no temp files.
- **Lower memory:** No Python interpreter, no subprocess pipes.
- **Reproducible:** Pure Rust with deterministic behavior across platforms.
