---
title: Output files
description: Naming conventions for trimmed FASTQ files and trimming reports.
---

Trim Galore writes one trimmed FASTQ per input, plus a per-input text and JSON trimming report. File naming matches v0.6.x exactly, so existing pipelines continue to work without changes.

## Single-end

| File | Contents |
|------|----------|
| `INPUT_trimmed.fq.gz` | Trimmed reads. |
| `INPUT_trimming_report.txt` | Cutadapt-compatible text report. |
| `INPUT_trimming_report.json` | Structured JSON report (schema v1, for MultiQC). |

`INPUT` is the original filename minus the `.fastq` / `.fq` (and `.gz`) suffix.

## Paired-end

| File | Contents |
|------|----------|
| `R1_val_1.fq.gz` | Validated Read 1, post-trim and post-pair-validation. |
| `R2_val_2.fq.gz` | Validated Read 2. |
| `R1.fastq.gz_trimming_report.txt` / `.json` | Read 1 report. |
| `R2.fastq.gz_trimming_report.txt` / `.json` | Read 2 report. Carries the final pair counts. |

## Singleton (unpaired) reads

With `--retain_unpaired`:

| File | Contents |
|------|----------|
| `R1_unpaired_1.fq.gz` | Read 1 reads whose mate dropped below `--length_2`. |
| `R2_unpaired_2.fq.gz` | Read 2 reads whose mate dropped below `--length_1`. |

## Specialty modes

Specialty modes write to mode-specific filenames and exit before the main pipeline:

| Mode | Output |
|------|--------|
| `--hardtrim5 N` | `*.{N}bp_5prime.fq(.gz)` |
| `--hardtrim3 N` | `*.{N}bp_3prime.fq(.gz)` |
| `--clock` | `*.clock_UMI.R1.fq(.gz)` / `*.clock_UMI.R2.fq(.gz)` |
| `--implicon[=N]` | `*_{N}bp_UMI_R1.fastq(.gz)` / `*_{N}bp_UMI_R2.fastq(.gz)` |
| `--demux` | One file per barcode (and per side for paired-end). |

## Compression

Gzip-compressed input produces gzip-compressed output by default. Pass `--dont_gzip` to write plain FASTQ.

## Output directory

`--output_dir DIR` writes outputs to `DIR/` instead of the current working directory. The trimmed FASTQ filename stem is unchanged; only the parent directory differs.

## Renaming outputs

`--rename PREFIX` replaces the input filename stem in the output names. Useful for pipelines that thread sample IDs through trimming separately from input filenames.

## FastQC

`--fastqc` runs FastQC on the trimmed output files after pair validation. Reports land alongside the trimmed FASTQ. `--fastqc_args "..."` passes additional FastQC arguments through verbatim.
