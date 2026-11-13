---
title: How Trim Galore works
description: A walkthrough of the trimming pipeline. Quality, adapter, length, and specialty steps.
---

Trim Galore reads FASTQ (or unaligned BAM, auto-detected from content), applies quality and adapter trimming in a single pass, filters short reads, and writes trimmed FASTQ — or uBAM, with `--output-format ubam` — plus a per-input report. Specialty modes (RRBS, hard-trim, IMPLICON preprocessing, demultiplexing) layer on top. The pipeline is the same for any base-space high-throughput data; bisulfite/RRBS users should also read [Bisulfite & RRBS](/rrbs/guide/) for library-specific guidance.

## Adaptive quality and adapter trimming with Trim Galore

Trim Galore handles quality trimming, adapter detection, adapter removal, length filtering, and specialty modes in a single pass over the data, with optional post-trimming quality reporting via the bundled [fastqc-rust](https://crates.io/crates/fastqc-rust) library (FastQC 0.13.0-compatible HTML + ZIP, no Java required). Input format (FASTQ vs unaligned BAM) is auto-detected from file content, not extension; output is FASTQ by default and can be uBAM via `--output-format ubam`, with aux tags round-tripping through `--preserve-tags`. Originally a Perl wrapper around [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and an external [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) install, v2.x consolidates those passes into one binary, keeping the output filename conventions and all but a few of the v0.6.x flags.

:::note
Every step below happens in a single pass over the data rather than as sequential Cutadapt invocations. Outputs match v0.6.x for the core flag set, which CI checks byte-for-byte against Perl 0.6.11; newer capabilities such as `--poly_g` auto-detection and the generic `--poly_a` trimmer extend beyond it.
:::

The pipeline runs in this order:

1. [Quality trimming](/guide/quality/). Phred-based, from the 3' end (BWA algorithm).
2. [Adapter trimming](/guide/adapters/). Auto-detection or manual sequence; semi-global alignment with a configurable error rate and overlap.
3. RRBS-specific trimming (only with `--rrbs`). Removes 2 bp from adapter-trimmed reads at MspI sites. See the [RRBS mode page](/modes/rrbs/).
4. [Length filtering](/guide/length/). Drops reads (or read pairs) below the cutoff.
5. [Paired-end validation](/guide/paired-end/). Discards pairs where at least one read became too short.
6. Optional [demultiplexing](/modes/demux/) and FastQC.

[Specialty modes](/modes/rrbs/) (`--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`) run their own short-circuit pipeline and exit before quality and adapter trimming.

## Where to go next

- New to Trim Galore? Start with [Quick start](/quickstart/).
- Want to know what a flag does or how it interacts with another? See [Flag reference notes](/guide/flags/).
- Working with bisulfite or RRBS? Read [Bisulfite & RRBS](/rrbs/guide/) before tweaking flags.
- Looking up output files and the report format? See [Output files](/guide/outputs/) and [Trimming reports](/guide/reports/).
