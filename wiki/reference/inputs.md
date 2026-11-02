---
type: reference
title: Accepted input files
description: How Trim Galore counts and pairs input files, how it decides what each file is, and the auxiliary input files the flags take.
tags: [input, fastq, ubam, formats]
page_id: inputs
subject_status: current
updated: 2026-11-02
sources:
  - repo: FelixKrueger/TrimGalore
    path: src/cli.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/main.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/format.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/fastq.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/bam.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/demux.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/adapter.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/io.rs
    commit: 6303152
---

## Summary

Read inputs are positional arguments, and at least one is required. Three
things are decided separately: how many files there are and how they pair up,
what each file contains, and whether the file can be read at all. The content
of the file decides its format; the filename never does.

## How many files, and how they pair

`Cli::validate` enforces these.

| Mode | Files | Rule |
|---|---|---|
| Single-end | 1 or more | Each file is processed independently, in a loop. |
| Paired-end (`--paired`) | Even number | Consecutive files form R1/R2 pairs. An odd count is refused. |
| Paired-end, one file | 1 | Legal only as the interleaved-uBAM path. FASTQ output on this shape writes no trimming report, so `--fastqc` is refused there; `--output-format ubam` does produce one. |
| Multi-pair | Even number above 2 | Each pair is announced as "Pair N of M". The specialty modes accept multi-pair too. |

`--basename` applies to one output stem, so `Cli::validate` restricts it once
more than one input (single-end) or more than one pair (paired-end) is given.

`Cli::validate` also refuses an input repeated on one command line, in three
separate forms: the two files of one pair being the same file, a later pair
duplicating an earlier one, and the same file given twice as single-end input.
Each would otherwise write one output twice.

## What each file is

`format::detect_input_format` reads the opening bytes and returns one of three
formats. Nothing here depends on the file extension.

| Format | How it is recognised | Notes |
|---|---|---|
| `FastqPlain` | First byte is `@` | Uncompressed FASTQ. |
| `FastqGz` | Gzip magic | Covers plain gzip and BGZF framing, so `bgzip`-compressed FASTQ lands here. |
| `UnalignedBam` | Gzip magic, then `BAM\1` in the decompressed payload | The payload probe is the discriminator: BGZF framing alone does not make a file a BAM. |

`format::open_sync_reader` and `format::open_threaded_reader` turn a detected
format into a reader, so callers share one input boundary.

uBAM input is read by `bam::BamReader`. Every record must be unmapped; an
aligned record is refused when it is reached, not only at the start, so a
mixed file cannot quietly produce trimmed output. Paired uBAM expects one
interleaved file with mates adjacent.

## Extensions

Extensions affect output naming only. `io::strip_fastq_extensions` recognises
`.gz`, `.bgz` and `.bgzf` alongside the FASTQ extensions, and matches them
case-insensitively, so an uppercase `.FASTQ.GZ` input names its outputs the
same way a lowercase one does.

## Auxiliary input files

These flags take a file that is not a read file.

| Flag | File | Content |
|---|---|---|
| `--passthrough` | Carrier read file | A third FASTQ read in lockstep with an R1/R2 pair, for 10X Multiome and scATAC cell-barcode reads. Requires `--paired` and exactly one pair. |
| `--demux` | Barcode file | TSV, two tab-separated fields per line: sample name, then barcode sequence. `demux::read_barcode_file` requires barcodes of equal length containing only A, C, T, G and N. |
| `-a`, `-a2` | Adapter FASTA | A value prefixed `file:` loads adapter sequences from a FASTA file, handled in `adapter.rs`. Both flags also accept `A{N}` shorthand for a repeated base, and are repeatable for multiple adapters. |

## Refused input

`main::sanity_check_any` checks the first input before dispatch, and routes by
detected format: FASTQ goes to `FastqReader::sanity_check_with`, uBAM to a
peek at the first record. On multi-pair input the remaining files are checked
as each pair is reached, so a bad file in a later pair is caught before that
pair is written.

The FASTQ check reads the first record and refuses an empty file, a first line
that does not start with `@`, and SOLiD colorspace data. Separately,
`check_restartable_input` refuses anything that cannot be re-read from the
start, because format detection, the sanity check and adapter auto-detection
each open the input independently.

`troubleshooting.md` lists each refusal with the remediation the tool prints.

## Links

- `troubleshooting.md` — the refusals above, and what to do about each.
- `reference/cli-flags.md` — the full flag surface, including every flag named here.
- `reference/outputs.md` — what the accepted input produces.
- `usage/ubam.md` — running uBAM input and uBAM output.

## Open questions

- None.
