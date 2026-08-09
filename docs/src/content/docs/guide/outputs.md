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
| `--demux` | One file per barcode (single-end input only). |

## Compression

Gzip-compressed input produces gzip-compressed output by default. Pass `--dont_gzip` to write plain FASTQ.

## uBAM output (`--output-format ubam`)

Emit records as unaligned BAM (uBAM) instead of FASTQ. Useful for pipelines that carry BAM downstream (10X single-cell, methylation callers, anything that wants to preserve BAM aux tags through the trim step).

| Input | uBAM output |
|------|-------------|
| Single-end | `INPUT_trimmed.bam` |
| Paired-end | `INPUT_val.bam` — **one** interleaved BAM per pair |

Paired output follows the samtools/Picard/fgbio convention: R1 and R2 records interleave in a single BAM with `FREAD1` (`0x40`) / `FREAD2` (`0x80`) flag bits. No separate `_val_1.bam` / `_val_2.bam` files.

uBAM output is currently single-threaded — `--cores N` is ignored for BAM writing; FastQC (when `--fastqc` is requested) still uses the value.

### `@PG` chain preservation

The input `@HD` / `@PG` header chain is propagated verbatim, and a trim_galore `@PG` line is appended:

```
@PG	ID:trim_galore	VN:<version>	CL:<command-line>
```

uBAM-in → uBAM-out is therefore **not** byte-identical to the input — provenance is preserved by *adding* to history, same treatment as samtools / Picard. Record bodies round-trip losslessly on that path, apart from IUPAC degenerate bases, which are coerced to `N` on read (a warning is emitted).

Coming **from FASTQ**, uBAM output applies three further normalisations, because a BAM record cannot represent everything a FASTQ header can: header text after the first space is dropped (BAM read names cannot contain whitespace), lowercase bases are uppercased, and IUPAC codes are coerced to `N`. Aux tags are carried only when named in `--preserve-tags`. A one-time `NOTE:` is printed when a header description is dropped.

### Aux-tag round-trip with `--preserve-tags`

`--preserve-tags TAG1,TAG2,...` (samtools `-T`-compatible, comma-separated) carries listed aux tags from input records through trimming into output records. Supports `A` / `Z` / `i` / `f` scalar types; array (`B`) and hex (`H`) types are rejected because the FASTQ intermediate cannot textually encode them.

Typical uses:

- `--preserve-tags CB,UB` — 10X cell / UMI barcodes for single-cell
- `--preserve-tags RG,LB,BC` — read-group and library metadata

FASTQ input has no source tags to carry, so `--preserve-tags` does nothing there — and it is not silent: on its own it warns, and combined with `--output-format ubam` it is a hard error (the flag was asked for explicitly, so there is no pressure-valve for the otherwise-invisible no-op).

### Feature compatibility

Most FASTQ-shaped features work with uBAM output. Rejected at CLI-validate time (v1 scope):

- `--dont_gzip` — gzip flag is FASTQ-output-specific; uBAM is BGZF-framed unconditionally
- `--clock`, `--implicon` — specialty modes encode UMI in the FASTQ header, which can't round-trip through the BAM record
- `--demux` — one output file per barcode; uBAM demux is a planned follow-up
- `--passthrough` — carrier FASTQ shape is FASTQ-specific
- `--retain_unpaired` — would produce singleton BAMs; not in v1
- `--clumpify` — in-place reorder is FASTQ-shape-specific (note: `--clump_only --output-format ubam` **is** supported — see [Clump-only](/modes/clump-only/))

## Output directory

`--output_dir DIR` writes outputs to `DIR/` instead of the current working directory. The trimmed FASTQ filename stem is unchanged; only the parent directory differs.

## Renaming outputs

`--rename PREFIX` replaces the input filename stem in the output names. Useful for pipelines that thread sample IDs through trimming separately from input filenames.

## FastQC

`--fastqc` runs the bundled [`fastqc-rust`](https://crates.io/crates/fastqc-rust) library on the trimmed output files after pair validation, producing FastQC 0.12.1-compatible HTML + ZIP reports alongside the trimmed output. Works on both FASTQ output (default) and uBAM output (`--output-format ubam`) — `fastqc-rust` reads `.bam` natively, so no intermediate conversion. No Java or external `fastqc` install needed.

For paired-end runs, FastQC is invoked once per output file — so PE FASTQ produces two reports (one per mate), while PE uBAM produces a **single** report per pair (uBAM PE output is one interleaved BAM, and the report covers both mates pooled).

`--fastqc_args "..."` passes a subset of FastQC flags through. Currently supported: `--nogroup`, `--expgroup`, `--quiet`, `--svg`, `--nano`, `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`. Other flags emit a warning and are ignored.
