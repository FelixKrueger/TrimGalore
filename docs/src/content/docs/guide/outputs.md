---
title: Output files
description: Naming conventions for trimmed FASTQ files and trimming reports.
---

Trim Galore writes one trimmed FASTQ per input, plus a per-input text and JSON trimming report. File naming matches v0.6.x exactly. Output *locations* match it too, with one exception: see [Output directory](#output-directory) for where a paired run's files land when `--output_dir` is not given.

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

All four land in the same directory — `--output_dir` if given, otherwise Read 1's. Each report keeps its own input-derived name; only the directory is shared. See [Output directory](#output-directory).

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

Most FASTQ-shaped features work with uBAM output. Rejected at startup, before anything is written (v1 scope):

- `--dont_gzip` — gzip flag is FASTQ-output-specific; uBAM is BGZF-framed unconditionally
- `--clock`, `--implicon` — specialty modes encode UMI in the FASTQ header, which can't round-trip through the BAM record
- `--demux` — one output file per barcode; uBAM demux is a planned follow-up
- `--passthrough` — carrier FASTQ shape is FASTQ-specific
- `--retain_unpaired` — would produce singleton BAMs; not in v1
- `--clumpify` — in-place reorder is FASTQ-shape-specific (note: `--clump_only --output-format ubam` **is** supported — see [Clump-only](/modes/clump-only/))
- `--rename` — only when at least one input is FASTQ, and only when a clipping flag is set; uBAM input is accepted (see [Annotating read IDs](#annotating-read-ids) below)

## Output directory

`--output_dir DIR` writes every output to `DIR/`. The filename stems are unchanged; only the parent directory differs.

Without `--output_dir`, the trimming modes write beside the **input**: single-end output lands in that input's directory, and for a pair every output — both validated FASTQs, both trimming reports, and the `--passthrough` carrier — lands in **Read 1's** directory. The specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`) are the exception: they always write to the current working directory.

Because output names derive from the input filename alone, two inputs sharing a filename resolve to the same report path once they land in one directory, and the run is refused before anything is written. That affects the layout where reads are split by mate rather than by sample:

```
R1/sample1.fq  R1/sample2.fq      # refused without --output_dir:
R2/sample1.fq  R2/sample2.fq      # sample1.fq's two reports collide in R1/
```

Pass `--output_dir` to collect the outputs somewhere unambiguous. The per-sample layout needs no change, because each pair's Read 1 already sits in its own directory:

```
sampleA/reads_1.fq  sampleA/reads_2.fq      # unaffected
sampleB/reads_1.fq  sampleB/reads_2.fq
```

## Refused runs

A run that fails part-way through writes no output file. Records go to a hidden `.partial` sibling of the final name and are renamed into place only once the writer closes cleanly, so a refusal — a malformed read name a million records in, a truncated input — leaves nothing behind, and a previous run's output at the same path is untouched.

A leftover `.<name>.partial` means the process was killed mid-run, or that the final rename failed — in that case the error names it as holding the trimmed data. Otherwise it is safe to delete; the next run at the same output path overwrites it.

## Custom output basename

`--basename BASE` replaces the input filename stem in the **trimmed output** names:

| Mode | Output |
|------|--------|
| Single-end | `BASE_trimmed.fq.gz`, or `BASE_trimmed.bam` with `--output-format ubam` |
| Paired-end | `BASE_val_1.fq.gz` / `BASE_val_2.fq.gz`, or a single `BASE_val.bam` with `--output-format ubam` |

Useful for pipelines that thread sample IDs through trimming separately from input filenames. Only valid for one file (single-end) or one pair (paired-end); longer input lists are refused, because the output naming would be ambiguous.

Two limits worth knowing. Trimming reports keep their input-derived names (`INPUT_trimming_report.txt`) rather than following `BASE` — matching v0.6.x. And the specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`) ignore `--basename` entirely; their outputs are always named from the input's basename.

## Annotating read IDs

`--rename` does not affect filenames. It is a boolean that appends `:clip5:SEQ` and/or `:clip3:SEQ` to the **read IDs**, recording the bases removed by `--clip_R1/R2`, `--three_prime_clip_R1/R2` or `--hardtrim5/3` — each half only when that side was clipped. Without one of those flags it appends nothing. Commonly used to keep UMIs recoverable downstream.

Two combinations are refused. `--clump_only` rejects `--rename` at any output format, because it preserves record contents byte-identically. And `--output-format ubam` rejects it when **at least one input is FASTQ** and a clipping flag is set: the annotation is appended to the end of the read ID, so a header carrying text after the first space puts the annotation inside that text, and BAM read names cannot contain whitespace, so none of that tail reaches the output. Whether that happens depends on the individual header, so the whole run is refused rather than decided per record — a per-record decision would keep the annotation for some reads and silently drop it for others. uBAM input is accepted: BAM read names carry no description, so there the annotation lands on the name itself.

## FastQC

`--fastqc` runs the bundled [`fastqc-rust`](https://crates.io/crates/fastqc-rust) library on the trimmed output files after pair validation, producing FastQC 0.12.1-compatible HTML + ZIP reports alongside the trimmed output. Works on both FASTQ output (default) and uBAM output (`--output-format ubam`) — `fastqc-rust` reads `.bam` natively, so no intermediate conversion. No Java or external `fastqc` install needed.

Two situations write no report and reject the flag rather than ignoring it: the four specialty modes (`--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`), and `--paired` given a single interleaved uBAM with FASTQ output. For that last one, `--output-format ubam` does produce a report.

For paired-end runs, FastQC is invoked once per trimmed primary output — so PE FASTQ produces two reports (one per mate), while PE uBAM produces a **single** report per pair (uBAM PE output is one interleaved BAM, and the report covers both mates pooled). `--retain_unpaired`'s `*_unpaired_{1,2}` files get no report, and `--passthrough` adds a third report for its carrier output.

`--fastqc_args "..."` passes a subset of FastQC flags through. Currently supported: `--nogroup`, `--expgroup`, `--quiet`, `--svg`, `--nano`, `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`. Other flags emit a warning and are ignored.
