---
title: Quick start
description: Run your first Trim Galore command. Single-end, paired-end, and RRBS examples.
---

Trim Galore reads FASTQ files (plain or gzip-compressed) or unaligned BAM (auto-detected by content), trims adapters and low-quality bases, and writes trimmed FASTQ — or uBAM with `--output-format ubam` — plus a per-file trimming report. The defaults work for most Illumina libraries.

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

FastQC is built in via the bundled `fastqc-rust` library: no Java or external `fastqc` install needed. Outputs are FastQC 0.12.1-compatible HTML + ZIP files. Works on both FASTQ and uBAM output paths, but is rejected on the specialty modes and on `--paired` with a single interleaved uBAM writing FASTQ — see [FastQC](/guide/outputs/#fastqc).

## Unaligned BAM (uBAM)

uBAM input is auto-detected — no flag needed. Paired reads must arrive as a single interleaved BAM with mates adjacent (samtools `sort -n` / `collate` / Picard / fgbio all produce this order). Two separate BAM files are not supported and are rejected before any output is written:

```bash
# uBAM in → FASTQ out (default)
trim_galore sample.bam
trim_galore --paired interleaved.bam

# uBAM in → uBAM out, preserving CB/UB aux tags (10X single-cell shape)
trim_galore --output-format ubam --preserve-tags CB,UB sample.bam
```

Paired uBAM output is a single interleaved BAM (`<stem>_val.bam`), matching samtools/Picard/fgbio convention. See [Output files](/guide/outputs/#ubam-output---output-format-ubam) for the full contract.

### Read names containing whitespace

A read name containing whitespace is refused, **on every output format** — including FASTQ output, where a space or tab loses nothing today but a newline silently desynchronises the whole file. Trim Galore reads a FASTQ header as ending at the first whitespace, so an aux-tag tail after it cannot be carried, and a newline splits one record across five lines. Whitespace other than a space inside a preserved aux-tag value is refused for the same reason; a space in a tag value is legal and unaffected.

The SAM specification forbids whitespace in QNAME, but samtools will build such a file. For names containing spaces only, samtools can rewrite them:

```bash
samtools view -h in.bam \
  | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} {gsub(/ /,"_",$1); print}' \
  | samtools view -b -o fixed.bam -
```

`FS="\t"` is essential. Without it `$1` stops at the first space, the substitution matches nothing, and the file comes out unchanged — `samtools view -b` accepts it and Trim Galore still refuses it. A name containing a tab or newline cannot be represented in SAM text at all, so it has to be corrected at source.

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
