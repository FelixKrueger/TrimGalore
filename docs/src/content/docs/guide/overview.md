---
title: How Trim Galore works
description: A walkthrough of the trimming pipeline. Quality, adapter, length, and specialty steps.
---

For all high throughput sequencing applications, we would recommend performing some quality control on the data, as it can often straight away point you towards the next steps that need to be taken (e.g. with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). Thorough quality control and taking appropriate steps to remove problems is vital for the analysis of almost all sequencing applications. This is even more critical for the proper analysis of RRBS libraries since they are susceptible to a variety of errors or biases that one could probably get away with in other sequencing applications. In our [brief guide to RRBS](/TrimGalore/rrbs/guide/) we discuss the following points:

* poor qualities. Affect mapping, may lead to incorrect methylation calls and/or mis-mapping.
* adapter contamination. May lead to low mapping efficiencies, or, if mapped, may result in incorrect methylation calls and/or mis-mapping.
* positions filled in during end-repair will infer the methylation state of the cytosine used for the fill-in reaction but not of the true genomic cytosine.
* paired-end RRBS libraries (especially with long read length) yield redundant methylation information if the read pairs overlap.
* RRBS libraries with long read lengths suffer more from all of the above due to the short size-selected fragment size.

Poor base call qualities or adapter contamination are however just as relevant for 'normal', i.e. non-RRBS, libraries.

## Adaptive quality and adapter trimming with Trim Galore

Trim Galore handles quality trimming, adapter detection, adapter removal, length filtering, and specialty modes in a single pass over the data, with optional post-trimming quality reporting via the bundled [fastqc-rust](https://crates.io/crates/fastqc-rust) library (FastQC 0.12.1-compatible HTML + ZIP, no Java required). Originally a Perl wrapper around [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and an external [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) install, v2.x is a faithful Rust rewrite that consolidates those passes while preserving the CLI and output filename conventions.

Even though Trim Galore works for any (base space) high throughput dataset (e.g. downloaded from the SRA) this section describes its use mainly with respect to RRBS libraries.

:::note
The Oxidized Edition (v2.x) is a faithful Rust rewrite designed as a drop-in replacement for v0.6.x. All of the trimming steps below happen in a single pass over the data rather than as sequential Cutadapt invocations. Outputs match v0.6.x for the core feature set; newer capabilities (e.g. `--poly_g` auto-detection, a generic `--poly_a` trimmer) extend beyond the Perl version.
:::

The pipeline runs in this order:

1. [Quality trimming](/TrimGalore/guide/quality/). Phred-based, from the 3' end (BWA algorithm).
2. [Adapter trimming](/TrimGalore/guide/adapters/). Auto-detection or manual sequence; semi-global alignment with a configurable error rate and overlap.
3. RRBS-specific trimming (only with `--rrbs`). Removes 2 bp from adapter-trimmed reads at MspI sites. See the [RRBS mode page](/TrimGalore/modes/rrbs/).
4. [Length filtering](/TrimGalore/guide/length/). Drops reads (or read pairs) below the cutoff.
5. [Paired-end validation](/TrimGalore/guide/paired-end/). Discards pairs where at least one read became too short.
6. Optional [demultiplexing](/TrimGalore/modes/demux/) and FastQC.

[Specialty modes](/TrimGalore/modes/rrbs/) (`--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`) run their own short-circuit pipeline and exit before quality and adapter trimming.

## Where to go next

- New to Trim Galore? Start with [Quick start](/TrimGalore/quickstart/).
- Want to know what a flag does or how it interacts with another? See [Flag reference notes](/TrimGalore/guide/flags/).
- Working with bisulfite or RRBS? Read [Bisulfite & RRBS](/TrimGalore/rrbs/guide/) before tweaking flags.
- Looking up output files and the report format? See [Output files](/TrimGalore/guide/outputs/) and [Trimming reports](/TrimGalore/guide/reports/).
