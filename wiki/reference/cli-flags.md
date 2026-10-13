---
type: reference
title: Command-line flags
description: Every flag Trim Galore accepts, grouped by what it affects, with defaults and the deprecated no-ops.
tags: [cli, flags, defaults]
page_id: cli-flags
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/cli.rs
    commit: bf28cbd
  - repo: TrimGalore
    path: src/library.rs
    commit: 51fa55e
  - repo: TrimGalore
    path: src/adapter.rs
    commit: d2ea293
---

## Summary

Flags are defined in one `Cli` struct in `cli.rs`, and cross-flag rules live
in `Cli::validate`. Reads are positional arguments, not a flag. A bare
`trim_galore` with no arguments prints the full help and exits zero rather
than erroring, matching the samtools and git convention.

Before parsing, `cli::rewrite_perl_short_flags` rewrites the Perl-era short
forms, so `-r1` and `-r2` still reach `--length_1` and `--length_2`.

## Trimming

| Flag | Default | Effect |
|---|---|---|
| `-q`, `--quality` | 20 | Phred cutoff for 3' quality trimming. |
| `-a`, `--adapter` | auto-detected | Adapter sequence. Repeatable; accepts `A{N}` shorthand and `file:adapters.fa`. |
| `--adapter2` | auto-detected | Same, for R2. |
| `-e`, `--error_rate` | 0.1 | Maximum error rate in an adapter match. |
| `--stringency` | 1 | Minimum overlap for an adapter match. |
| `-n`, `--times` | 1 | Adapter-removal rounds per read. |
| `--nextseq` | off | 2-colour quality trimming; G is treated as quality zero. |
| `--poly_a` | off | Trim poly-A (R1) and poly-T (R2). |
| `--poly_g` / `--no_poly_g` | auto | Force 3' poly-G trimming on or off, overriding auto-detection. |
| `--trim-n` | off | Trim Ns from both ends. Suppressed under `--rrbs`. |
| `--consider_already_trimmed` | unset | Takes a length; above it, adapter trimming is skipped and only quality trimming runs. Conflicts with the preset flags. |
| `--phred64` / `--phred33` | Phred+33 | Input quality encoding. |

## Adapter presets

Mutually exclusive with each other and with `--consider_already_trimmed`.
Each selects a built-in sequence instead of auto-detecting one.

| Flag | Preset |
|---|---|
| `--illumina` | Illumina |
| `--nextera` | Nextera |
| `--small_rna` | smallRNA |
| `--bgiseq` | BGI/DNBSEQ |
| `--stranded_illumina` | Stranded Illumina |

Auto-detection probes the first four only. Stranded Illumina is explicit-only:
its sequence is the Nextera sequence with one extra leading base, so guessing
between them from read data would be a coin flip.

## Filtering

| Flag | Default | Effect |
|---|---|---|
| `--length` | unset | Discard reads shorter than this after trimming. |
| `--max_length` | unset | Discard reads longer than this. |
| `--max_n` | unset | Discard reads with more than this many Ns. |
| `--discard_untrimmed` | off | Keep only reads where an adapter was found. |
| `--retain_unpaired` | off | Write rescued singletons to `_unpaired_1/2` files. |
| `-r`, `--length_1` | 35 | Minimum R1 length to be rescued as unpaired. |
| `--length_2` | 35 | Same for R2. |

## Clipping and presets

| Flag | Default | Effect |
|---|---|---|
| `--clip_R1`, `--clip_R2` | unset | Fixed 5' clip. |
| `--three_prime_clip_R1`, `--three_prime_clip_R2` | unset | Fixed 3' clip. |
| `--library` | unset | Names a library prep instead of setting four clip flags by hand; see `library.rs`. |
| `--rename` | off | Append the clipped bases to the read ID. |

## Output

| Flag | Default | Effect |
|---|---|---|
| `-o`, `--output_dir` | input's directory | Where outputs land. |
| `--basename` | derived from input | Output stem. Restricted once there is more than one input or pair. |
| `--dont_gzip` | off | Write plain text even when the input is gzipped. |
| `--compression` | crate default level | Gzip level for output. Does not decide *whether* output is gzipped. |
| `--output-format` | `fastq` | `ubam` writes unaligned BAM instead. |
| `--preserve-tags` | none | BAM aux tags to carry into the FASTQ header, samtools `-T`-compatible. |
| `--no_report_file` | off | Suppress the trimming report file. |

## Performance and clumping

| Flag | Default | Effect |
|---|---|---|
| `-j`, `--cores` | 1 | Worker count. Each worker trims and compresses its own chunk. |
| `--clumpify` | off | Reorder reads into bins for tighter compression while trimming. |
| `--clump_only` | off | Reorder only, with no trimming. |
| `--memory` | 1G | Budget the clumping bin pool is resolved against. |

## Modes

| Flag | Effect |
|---|---|
| `--paired` | Consecutive inputs form R1/R2 pairs. |
| `--rrbs` | MspI RRBS end-repair trimming; auto-sets `--clip_r2 2` on directional paired input. |
| `--non_directional` | Non-directional RRBS: reads starting `CAA` or `CGA` lose 2 bp at the 5' end. |
| `--hardtrim5`, `--hardtrim3` | Hard-trim to a fixed length from one end. Mutually exclusive. |
| `--clock` | Epigenetic-clock UMI handling. |
| `--implicon` | UMI taken from R2. |
| `--demux` | 3' inline barcode demultiplexing, single-end only. |
| `--passthrough` | Carry a third read file through in lockstep. Requires `--paired` and one pair. |

The four specialty modes — `--hardtrim5`, `--hardtrim3`, `--clock`,
`--implicon` — run and exit instead of entering the trimming pipeline, and
each owns its output naming.

## QC

| Flag | Effect |
|---|---|
| `--fastqc` | Run the bundled FastQC on the trimmed output. |
| `--fastqc_args` | A curated subset of FastQC flags; unknown flags warn and are ignored. |

`--fastqc` is refused on the four specialty modes and on one interleaved uBAM
with FASTQ output, because those paths write no report.

## Deprecated and no-op

Seven flags are accepted only so that v0.6.x command lines keep parsing. They
are listed in `reference/legacy-flags.md`.

## Links

- `reference/inputs.md` — the positional reads and the auxiliary input files.
- `reference/outputs.md` — what the output flags produce.
- `troubleshooting.md` — the flag combinations that are refused.

## Open questions

- `--library`'s preset names and the clip values each expands to are not
  listed here; they belong on a page for `library.rs`, deferred to a later
  update.
