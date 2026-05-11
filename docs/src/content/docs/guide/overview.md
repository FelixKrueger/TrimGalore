---
title: How Trim Galore works
description: A walkthrough of the trimming pipeline. Quality, adapter, length, and specialty steps.
---

Trim Galore reads FASTQ, applies quality and adapter trimming in a single pass, filters short reads, and writes trimmed FASTQ plus a per-input report. Specialty modes (RRBS, hard-trim, IMPLICON preprocessing, demultiplexing) layer on top. The pipeline is the same for any base-space high-throughput data; bisulfite/RRBS users should also read [Bisulfite & RRBS](/rrbs/guide/) for library-specific guidance.

## Adaptive quality and adapter trimming with Trim Galore

Trim Galore handles quality trimming, adapter detection, adapter removal, length filtering, and specialty modes in a single pass over the data, with optional post-trimming quality reporting via the bundled [fastqc-rust](https://crates.io/crates/fastqc-rust) library (FastQC 0.12.1-compatible HTML + ZIP, no Java required). Originally a Perl wrapper around [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and an external [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) install, v2.x is a faithful Rust rewrite that consolidates those passes while preserving the CLI and output filename conventions.

:::note
Trim Galore v2.x is a faithful Rust rewrite designed as a drop-in replacement for v0.6.x. All of the trimming steps below happen in a single pass over the data rather than as sequential Cutadapt invocations. Outputs match v0.6.x for the core feature set; newer capabilities (e.g. `--poly_g` auto-detection, a generic `--poly_a` trimmer) extend beyond the Perl version.
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
