---
type: usage
title: Trimming single-end and paired-end reads
description: The ordinary single-end and paired-end runs, what happens without any flags, and the multi-pair case.
tags: [usage, single-end, paired-end, adapters]
page_id: trim-reads
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/main.rs
    commit: 26ec1cf
  - repo: TrimGalore
    path: src/cli.rs
    commit: bf28cbd
  - repo: TrimGalore
    path: src/adapter.rs
    commit: d2ea293
  - repo: TrimGalore
    path: src/trimmer.rs
    commit: bf0be75
---

## Summary

Reads are positional arguments. One file is a single-end run; `--paired` makes
consecutive files into R1/R2 pairs. With no other flags the adapter is
detected from the data, quality trimming uses the default cutoff, and a
trimming report is written next to the output.

## Single-end

```bash
trim_galore sample.fastq.gz
```

Several files in one command are processed independently, one after another:

```bash
trim_galore a.fastq.gz b.fastq.gz c.fastq.gz
```

## Paired-end

```bash
trim_galore --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

Consecutive files pair up, so more than one pair can go in one command. Each
pair is announced as "Pair N of M" as it starts:

```bash
trim_galore --paired s1_R1.fq.gz s1_R2.fq.gz s2_R1.fq.gz s2_R2.fq.gz
```

An odd number of files is refused, as is the same file twice, a pair whose two
files are identical, and a later pair that repeats an earlier one.

Before anything is written, a pre-flight checks that no two planned outputs
resolve to the same path — including case-only differences, which matter on
macOS and Windows filesystems. See `troubleshooting.md`.

## What happens without flags

| Behaviour | Detail |
|---|---|
| Adapter | Auto-detected by scanning the start of the file, up to one million reads. In paired mode this runs **per pair**, not once for the run. |
| Quality trimming | 3' trimming at the default Phred cutoff; `-q` changes it. |
| Poly-G | Switched on by itself when the same pre-scan finds enough reads with long trailing G runs, which is the 2-colour chemistry signature. `--no_poly_g` forces it off. |
| Compression | Output follows the input's encoding, not the compression level. A gzipped input gives gzipped output; `--dont_gzip` forces plain text. |
| Report | A trimming report is written per input, in text and JSON. `--no_report_file` suppresses it. |

Auto-detection probes the built-in adapter sets. `--stranded_illumina` stays
explicit-only, because its sequence is ambiguous with Nextera and guessing it
would be a coin flip.

## Choosing the adapter yourself

```bash
trim_galore -a AGATCGGAAGAGC sample.fastq.gz          # one sequence
trim_galore -a SEQ1 -a SEQ2 sample.fastq.gz           # several
trim_galore -a 'A{10}' sample.fastq.gz                # shorthand
trim_galore -a file:adapters.fa sample.fastq.gz       # from FASTA
```

`--adapter2` does the same for R2. When only `-a` is given, R2 uses the same
set.

## Going faster

```bash
trim_galore --cores 8 --paired sample_R1.fq.gz sample_R2.fq.gz
```

Each worker trims and compresses its own chunk, so the scaling holds until
output I/O becomes the limit rather than CPU. Note that output bytes can
differ between core counts: the records and their order do not change, but
gzip member boundaries do.

## Links

- `usage/rrbs.md` — RRBS and non-directional runs.
- `usage/clumpify.md` — compression and archival.
- `usage/ubam.md` — uBAM in and out.
- `reference/outputs.md` — what the output files are called.
- `reference/cli-flags.md` — every flag with its default.

## Open questions

- None.
