---
type: reference
title: Output files and reports
description: Where outputs land, how each one is named, whether they are gzipped, and what the trimming reports contain.
tags: [output, naming, reports, json]
page_id: outputs
subject_status: current
updated: 2026-10-13
sources:
  - repo: FelixKrueger/TrimGalore
    path: src/io.rs
    commit: c7e147d
  - repo: FelixKrueger/TrimGalore
    path: src/report.rs
    commit: 51fa55e
---

## Summary

Output names are derived from the input name, with the FASTQ extensions
stripped, plus a suffix per output kind. Reports are the exception: they keep
the input filename whole. Whether output is gzipped follows the input rather
than the compression level.

## Where outputs land

`io::pair_output_dir` gives the directory: `--output_dir` when given, else the
parent of R1. Every output of a paired run shares it — the two primaries, the
reports, and the `--passthrough` carrier — so a pair never scatters its
outputs across two directories.

## Read output names

The stem is `io::strip_fastq_extensions` of the input, or `--basename` when
given. `ext` is `.fq.gz` or `.fq`.

| Output | Name | Written by |
|---|---|---|
| Single-end | `<stem>_trimmed<ext>` | `io::single_end_output_name` |
| Paired-end | `<stem>_val_1<ext>`, `<stem>_val_2<ext>` | `io::paired_end_output_names` |
| Rescued singletons | `<stem>_unpaired_1<ext>`, `<stem>_unpaired_2<ext>` | `io::unpaired_output_names` |
| Passthrough carrier | `<stem>_passthrough<ext>` | `io::passthrough_output_name` |
| Clumped, single-end | `<stem>_clumped<ext>` | `io::clumped_output_name` |
| Clumped, paired | `<stem>_clumped_1<ext>`, `<stem>_clumped_2<ext>` | `io::clumped_paired_output_names` |
| uBAM, single-end | `<stem>_trimmed.bam` | `io::single_end_bam_output_name` |
| uBAM, paired | `<stem>_val.bam` | `io::paired_bam_output_name` |
| uBAM, clumped | `<stem>_clumped.bam` | `io::clumped_bam_output_name` |

Paired uBAM output is one interleaved file, not two, which matches the
samtools, Picard and fgbio convention.

## Whether output is gzipped

`io::output_is_gzipped` decides from the first input's encoding, unless
`--dont_gzip`. The configured compression level does not enter into it, so a
plain-text input produces plain-text output even when a level is set. Any
message that describes the output has to read this rather than the level, or
it can promise gzip over plain text.

## Report names

Reports keep the input filename intact, extensions included, because
`io::report_name` uses the file name rather than the stripped stem. An input
`sample.fastq.gz` therefore produces `sample.fastq.gz_trimming_report.txt`,
not `sample_trimming_report.txt`.

| Report | Name | Written by |
|---|---|---|
| Trimming report, text | `<input filename>_trimming_report.txt` | `io::report_name` |
| Trimming report, JSON | `<input filename>_trimming_report.json` | `io::json_report_name` |
| Clumping report | `<input filename>_clumping_report.txt` | `io::clumping_report_name` |

## What the text report contains

`report.rs` writes it in sections, each from its own function:

| Section | Function |
|---|---|
| Run parameters and adapter settings | `report::write_report_header` |
| Per-run counts and trimming totals | `report::write_run_stats` |
| Pair validation counts, paired runs only | `report::write_pair_validation_stats` |
| Carrier counts, `--passthrough` runs only | `report::write_passthrough_stats` |
| A section laid out like Cutadapt's | `report::write_cutadapt_compatible_section` |
| Closing summary | `report::write_run_footer` |

The Cutadapt-compatible section is what makes the report readable by tools
built for Cutadapt output, MultiQC included. Its "Total written (filtered)"
figure comes from the read length captured after adapter trimming, before the
later steps — see `internals/trimming-pipeline.md`.

`report::write_json_report` writes the JSON form from the same statistics.

## Before anything is written

Two guards run on the output side, both covered in `troubleshooting.md`:
`io::preflight_output_collisions` refuses a run whose planned outputs collide,
comparing paths case-folded through `io::collision_key`; and `io::PendingOutput`
publishes each file only once it is complete, so a refused run leaves no
half-written output.

## Links

- `reference/inputs.md` — what produces these outputs.
- `troubleshooting.md` — the collision refusal and the atomic-publish guard.
- `internals/trimming-pipeline.md` — where the report figures come from.

## Open questions

- None.
