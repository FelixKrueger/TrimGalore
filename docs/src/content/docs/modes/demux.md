---
title: Demultiplexing
description: 3' inline barcode demultiplexing with `--demux`.
---

Trim Galore can demultiplex FASTQ files by **3' inline barcode** in a post-trim pass. Pass a barcode file via `--demux` (single-end) or `--demux` plus `--barcode2` for paired-end libraries with R1 and R2 inline barcodes.

## How it works

After quality and adapter trimming, each read's 3' end is matched against the barcode list. Reads are written to per-barcode output files. Reads that do not match any barcode go to an "ambiguous" bucket.

For paired-end with `--barcode2`, both R1 and R2 must match their respective barcode lists. Pairs where one side fails to demultiplex go to ambiguous.

## Usage

```bash
trim_galore --demux barcodes.txt input.fq.gz
```

Per-barcode files are written alongside the trimmed FASTQ.

For paired-end:

```bash
trim_galore --demux r1_barcodes.txt --barcode2 r2_barcodes.txt --paired R1.fq.gz R2.fq.gz
```

## Notes

- Demultiplexing runs after the regular trimming pipeline. Quality, adapter, and length filters all apply first.
- Ambiguous reads are kept in a separate file, not discarded, so you can inspect what did not match.
- For high-throughput inline-demux workflows, dedicated tools (e.g. Cutadapt's demux mode) can be more flexible. `--demux` here is a built-in convenience for standard 3' barcode designs.
