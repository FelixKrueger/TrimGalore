---
title: Quick start
description: Run your first Trim Galore command. Single-end, paired-end, and RRBS examples.
---

Trim Galore reads FASTQ files (plain or gzip-compressed), trims adapters and low-quality bases, and writes trimmed FASTQ plus a per-file trimming report. The defaults work for most Illumina libraries.

## Single-end

```bash
trim_galore input.fastq.gz
```

Outputs:

- `input_trimmed.fq.gz` (trimmed reads)
- `input.fastq.gz_trimming_report.txt` (text report)
- `input.fastq.gz_trimming_report.json` (structured report for MultiQC)

## Paired-end

```bash
trim_galore --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

Outputs:

- `sample_R1_val_1.fq.gz` and `sample_R2_val_2.fq.gz` (validated paired reads)
- A trimming report per input file

## Parallel processing

Speedup is near-linear up to about 8 cores on v2.1.0-beta.7; beyond that, gzip-output I/O on the storage layer typically becomes binding and adding cores helps less. For nf-core / Snakemake / CWL workflows, `--cores 8` is also the saturation point.

```bash
trim_galore --cores 8 --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

## RRBS libraries

```bash
trim_galore --rrbs --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

For non-directional libraries, add `--non_directional`. See the [Bisulfite & RRBS guide](/rrbs/guide/) for the biology behind these modes.

## Run FastQC alongside

```bash
trim_galore --fastqc input.fastq.gz
```

FastQC is built in via the bundled `fastqc-rust` library: no Java or external `fastqc` install needed. Outputs are FastQC 0.12.1-compatible HTML + ZIP files.

## Common combinations

```bash
# Small RNA: auto-lowers --length to 18 bp
trim_galore --small_rna input.fastq.gz

# 2-colour aware quality trimming (NextSeq, NovaSeq): replaces -q
trim_galore --2colour 20 input.fastq.gz

# Force-disable poly-G trimming (it is auto-enabled on 2-colour data)
trim_galore --no_poly_g input.fastq.gz

# Trim a poly-A tail (mRNA-seq libraries)
trim_galore --poly_a input.fastq.gz

# Multiple adapters, up to 3 occurrences per read
trim_galore -a AGCTCCCG -a TTTCATTAT -a TTTATTCGGAT -n 3 input.fastq.gz

# Adapters from a FASTA file
trim_galore -a "file:./adapters.fa" input.fastq.gz
```

## Get the full flag list

```bash
trim_galore --help
```

For context on individual flags and how they interact, see the [user guide](/guide/overview/).
