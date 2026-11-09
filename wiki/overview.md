---
type: overview
title: What Trim Galore is for
description: The purpose, who runs it, and when it is and is not the right tool.
tags: [overview, ngs, trimming]
page_id: overview
subject_status: current
updated: 2026-11-09
sources:
  - repo: FelixKrueger/TrimGalore
    path: src/cli.rs
    commit: 6caf099
  - repo: FelixKrueger/TrimGalore
    path: src/format.rs
    commit: 6caf099
  - repo: FelixKrueger/TrimGalore
    path: src/fastq.rs
    commit: 6caf099
  - repo: FelixKrueger/TrimGalore
    path: src/bam.rs
    commit: 6caf099
  - repo: FelixKrueger/TrimGalore
    path: src/trimmer.rs
    commit: 6caf099
  - repo: FelixKrueger/TrimGalore
    path: src/main.rs
    commit: 6caf099
---

## Summary

Trim Galore is fast adapter and quality trimming for NGS data, with several
specialty modes on top, built-in QC, and a compression and archival mode. It
is one binary with no external runtime dependencies.

## Who runs it

Mostly pipeline users — nf-core workflows in particular — plus people running
it directly at the command line (developer's description).

The library target exists so that the binary and the integration tests can
reach internals. `src/lib.rs` states it carries no semver guarantee and is not
a published API, so it is not an interface to build against.

## When to use it

### Good fit

The developer's answers:

- Any standard Illumina adapter and quality trimming, with auto-detection of
  the common adapters.
- RRBS data, which has dedicated handling.
- Compression and archival, via clumpify.
- QC, which is built in rather than a separate step.

### Not a good fit

The developer's answer: extremely complex, non-standard sequencing use cases,
which are rare. Use Cutadapt for those.

### Limits in the code

| Limit | Enforced by |
|---|---|
| Input must be FASTQ, gzip or BGZF FASTQ, or unaligned BAM; nothing else is recognised | `format::detect_input_format` |
| SOLiD colorspace data is refused | `FastqReader::sanity_check_with` |
| Aligned BAM is refused, per record rather than only at the start | `bam::BamReader` |
| Input must be a regular file: pipes, FIFOs and process substitution are refused, because each input is opened more than once | `check_restartable_input` |
| Paired input needs an even file count; a single file is paired only as interleaved uBAM | `Cli::validate` |
| Paired uBAM must be one interleaved file with mates adjacent | `bam::BamReader::open_paired_interleaved` |
| Demultiplexing is single-end only | `demux` module |
| uBAM output is single-threaded in v1 and refuses seven flags, including `--clumpify` and `--demux` | `Cli::validate` |
| Aux tags round-trip as A, Z, i and f scalars only | `bam::parse_tag_value` |
| `--phred64` is refused with uBAM input, and verified against the input under `--fastqc` or with uBAM output | `lowest_quality_char` |
| `--trim-n` has no effect under `--rrbs` | `trimmer::trim_read` |
| Clumpify's memory reservation is calibrated above a floor read length; shorter reads are warned about, not covered | `main::warn_below_clumpify_floor` |

## Relationship to the Perl releases

Trim Galore v2.x is a Rust rewrite. The Perl v0.6.x version is deprecated
(developer's statement), and its source is reachable from the `0.6.11` tag
rather than from the working tree.

Byte-identity with Perl v0.6.11 is a hard invariant for the core flag paths,
not an aspiration: CI compares output against the Perl release for single-end,
paired-end, hardtrim5, clock and demux runs. Several deliberate behaviours
exist only to hold that line — the suppression of `--trim-n` under `--rrbs`,
and the point in the pipeline where the post-adapter read length is captured
for the report.

One deliberate divergence: adapter auto-detection runs per pair in paired
mode, where Perl v0.6.x detected once on the first file.

## Links

- `usage/install.md` — getting it, and which channel to use.
- `usage/trim-reads.md` — the ordinary run.
- `reference/cli-flags.md` — the full flag surface.
- `ENTRY.md` — routing by task.

## Open questions

- None.
