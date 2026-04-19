# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TrimGalore is a wrapper around Cutadapt and FastQC for quality/adapter trimming of NGS data. It is a single Perl script (`trim_galore`, ~3,900 lines) with no build system.

## Running and Testing

### Dependencies
- **Perl 5.12+** (with strict/warnings)
- **Cutadapt** (required; v1.15+ for multi-core; needs Python 3 for multi-core)
- **FastQC** (optional, used with `--fastqc`)
- **pigz** (optional, parallel gzip for multi-core compression)

### Running
```bash
./trim_galore [options] <input.fastq.gz>
./trim_galore --paired file_R1.fastq.gz file_R2.fastq.gz
```

### Testing
There is no formal test suite. Tests are defined in `.travis.yml` and run manually:
```bash
# Basic single-end tests
trim_galore test_files/illumina_10K.fastq.gz
trim_galore --illumina --nextera --small_rna --trim-n test_files/4_seqs_with_Ns.fastq.gz

# Expected failure cases (should produce errors)
trim_galore test_files/colorspace_file.fastq       # colorspace rejection
trim_galore test_files/truncated.fq.gz              # truncated file
trim_galore test_files/empty_file.fastq             # empty file
```

Test FASTQ files are in `test_files/` (illumina, nextera, smallRNA, polyAT, clock, edge cases).

## Architecture

The entire tool is a single Perl script with this flow:

1. **`process_commandline()`** (line ~2715) — parses 50+ CLI options via Getopt::Long
2. **Speciality modes** (checked first, exit after processing):
   - `--hardtrim5` / `--hardtrim3` — hard-trim to fixed length
   - `--clock` — Epigenetic Clock UMI extraction (`clockwork()`, line ~261)
   - `--implicon` — IMPLICON UMI transfer (`implicon_umi()`, line ~109)
3. **Main trimming loop** — calls `trim()` (line ~605) per input file:
   - Quality trimming (Phred-based, from 3' end)
   - Adapter detection via `autodetect_adapter_type()` (line ~2453)
   - Adapter removal via Cutadapt subprocess
   - RRBS-specific 2bp removal at MspI sites
   - Length filtering
   - Generates `*_trimming_report.txt`
4. **Paired-end validation** — `validate_paired_end_files()` (line ~2015) ensures pair correspondence
5. **Optional post-processing** — demultiplexing (`demux()`, line ~1690), FastQC

### Key adapter sequences (auto-detected)
- **Illumina**: `AGATCGGAAGAGC` (default)
- **Nextera**: `CTGTCTCTTATA`
- **Small RNA**: `TGGAATTCTCGG`

### Output naming conventions
- Single-end: `*_trimmed.fq(.gz)`
- Paired-end: `*_val_1.fq(.gz)` / `*_val_2.fq(.gz)`
- Unpaired: `*_unpaired_1.fq(.gz)` / `*_unpaired_2.fq(.gz)`
