---
title: Paired-end data
description: Single-pass paired-end processing, validation, and retaining singletons.
---

Trim Galore processes paired-end reads in a **single pass**. R1 and R2 are read, trimmed, and validated together. The per-read order is preserved, which most aligners require, and there are no temp files between Cutadapt invocations like the v0.6.x Perl wrapper used.

## Basic invocation

```bash
trim_galore --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

Outputs:

- `sample_R1_val_1.fq.gz` (validated R1)
- `sample_R2_val_2.fq.gz` (validated R2)
- `sample_R1.fastq.gz_trimming_report.txt` (and `.json`)
- `sample_R2.fastq.gz_trimming_report.txt` (and `.json`)

## Multiple pairs

`--paired` accepts any **even number** of input files and processes them as consecutive pairs in `R1, R2, R1, R2, …` order:

```bash
trim_galore --paired \
    s1_R1.fq.gz s1_R2.fq.gz \
    s2_R1.fq.gz s2_R2.fq.gz \
    s3_R1.fq.gz s3_R2.fq.gz
```

:::caution[Don't use `*_R1.fq.gz *_R2.fq.gz` globs]
A glob like `*_R1.fq.gz *_R2.fq.gz` produces **all R1s, then all R2s**, not interleaved pairs. Invoke once per sample, or interleave the pairs explicitly as shown above.
:::

## Pair validation

Quality and adapter trimming can shorten one read of a pair below the `--length` cutoff. The pair-validation step runs **after** per-read trimming and decides what happens to such pairs:

- If both reads pass `--length`, the pair is written to `*_val_1.fq.gz` and `*_val_2.fq.gz`.
- If at least one read is below `--length`, the pair is discarded by default.

The Read 2 trimming report carries the final pair counts at the end:

```
RUN STATISTICS FOR INPUT FILE: SLX_R2.fastq.gz
=============================================
1166076593 sequences processed in total

Total number of sequences analysed for the sequence pair length validation: 1166076593

Number of sequence pairs removed because at least one read was shorter than the length cutoff (20 bp): 3357967 (0.29%)
```

## Retaining singletons

When only one of the two reads dropped below the cutoff (e.g. one read had very poor qualities throughout), `--retain_unpaired` keeps the surviving read in a separate file:

```bash
trim_galore --paired --retain_unpaired sample_R1.fq.gz sample_R2.fq.gz
```

Additional outputs:

- `sample_R1_unpaired_1.fq.gz`. R1 reads whose mate was too short.
- `sample_R2_unpaired_2.fq.gz`. R2 reads whose mate was too short.

The per-side cutoff is governed by `--length_1` and `--length_2` (default 35 bp each).

## Per-side trimming

One side of a paired-end library can have biases at a specific end (e.g. RRBS Read 2 has fill-in artifacts at the 5' end). The clip flags remove fixed numbers of bases without touching the other read:

| Flag | Effect |
|------|--------|
| `--clip_R1 INT` | Remove `N` bases from the 5' end of Read 1. |
| `--clip_R2 INT` | Remove `N` bases from the 5' end of Read 2. |
| `--three_prime_clip_R1 INT` | Remove `N` bases from the 3' end of Read 1 (after quality / adapter trim). |
| `--three_prime_clip_R2 INT` | Remove `N` bases from the 3' end of Read 2. |

`--rrbs` in directional mode auto-sets `--clip_R2 2` to mask the 2 bp end-repair bias at the start of Read 2, unless you supply your own `--clip_R2`. `--non_directional` skips this auto-clip on purpose.

## `--discard_untrimmed`

`--discard_untrimmed` keeps only reads where at least one adapter match was found. For paired-end, the pair is discarded only when **neither** read had an adapter.
