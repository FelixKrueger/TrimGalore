---
title: Clump-only (lossless reorder)
description: Reorder FASTQ records for gzip-friendly compression without trimming (`--clump_only`).
---

`--clump_only` is a run-and-exit mode that reorders FASTQ records by canonical 16-mer minimizer for gzip-friendly compression, **without any trimming, filtering, or adapter detection**. Every input record appears in the output byte-identically; only file-level order changes.

Requested in [#353](https://github.com/FelixKrueger/TrimGalore/issues/353) as an archival/recompression path — same clumping win that [`--clumpify`](/performance/clumpy/) delivers as a compression stage inside the trim pipeline, but as a standalone lossless mode you can use to recompress raw reads for storage without ever mutating them.

## Usage

Single-end:

```bash
trim_galore --clump_only --compression 9 --cores 2 sample.fq.gz
```

Paired-end:

```bash
trim_galore --clump_only --paired --compression 9 --cores 2 sample_R1.fq.gz sample_R2.fq.gz
```

Output filenames:

- Single-end: `sample_clumped.fq.gz` (or `sample_clumped.fq` under `--dont_gzip`)
- Paired-end: `sample_R1_clumped_1.fq.gz` and `sample_R2_clumped_2.fq.gz`
- Report: `sample.fq.gz_clumping_report.txt` — a short text summary, deliberately distinct from the `_trimming_report.txt` glob that downstream nf-core/MultiQC pipelines scan (see the [Report shape](#report-shape) section).

## What's preserved and what isn't

**Byte-identical across the reorder:**

- Read ID (header line, including `@` and any description)
- Sequence
- Quality string

**Normalized codebase-wide (not `--clump_only`-specific):**

- **Plus-line (line 3 of each record)** is normalized to bare `+` on output. Inputs with `+<header-repeat>` on line 3 emerge with bare `+`. This applies to every Trim Galore code path, not just `--clump_only`, and is inherited from the `FastqReader`/`FastqWriter` behaviour.
- **CRLF line endings** are normalized to LF. Inputs with `\r\n` endings emerge with `\n` endings.

The byte-identity guarantee covers the three semantically-meaningful fields (header, sequence, quality). If you need byte-for-byte-exact recompression including line-3 headers or CRLF endings, `--clump_only` is not the right tool.

## Determinism

Two invocations with the same input, `--cores`, `--memory`, and `--compression` produce byte-identical output files. Downstream tools that md5 the output for archival integrity get stable hashes. This is enforced by a CI cross-run check.

## Compatibility

Composes with:

- `--paired`
- `--compression <1-9>` (gzip level; higher = smaller output, slower)
- `--memory <SIZE>` (bin-buffer sizing, e.g. `4G`)
- `--cores <N>` (accepted for interface parity — v1 is single-threaded internally; parallelism is planned for v1.1)
- `--fastqc` — the quality report is meaningful under reorder because record contents are unchanged
- `--dont_gzip` — produces `*_clumped.fq` (plain), diverging from `--clumpify`'s rejection of the same combination. `--clump_only`'s archival use case may want lossless plain output as an intermediate stage for a downstream compressor.
- `--basename BASE` — output becomes `BASE_clumped.fq(.gz)` (SE) or `BASE_clumped_{1,2}.fq(.gz)` (PE)

Rejected at CLI validation (would break byte-identity or has no meaning under this mode):

- `-a` / `-a2` / `--adapter` / adapter presets (`--illumina`, `--nextera`, `--small_rna`, `--bgi`, `--stranded_illumina`)
- `--length` / `--max_length` / `--max_n`
- `--trim-n` / `--clip_r1` / `--clip_r2` / `--three_prime_clip_r1` / `--three_prime_clip_r2`
- `--rrbs` / `--non_directional`
- `--polyA` / `--polyG` / `--no_poly_g`
- `--nextseq` / `--2colour`
- `--rename` (mutates read IDs)
- `--discard_untrimmed`, `--consider_already_trimmed`
- Other specialty modes: `--clumpify` (redundant), `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`, `--demux`
- `--output-format ubam`, `--passthrough`, `--retain_unpaired`

**Silently ignored** (kept in the CLI surface for parity with the trim path, but do nothing under this mode; documented so users don't wonder why they weren't rejected): `-q` / `--quality`, `--stringency`, `-e` / `--error`.

`uBAM input` is rejected at input-detection time: v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up.

## Report shape

The reorder report at `<stem>_clumping_report.txt` is deliberately narrower than the standard trimming report — no adapter counts, no quality stats, no filter counters, no JSON companion.

```
Trim Galore version: 2.x.x
Mode: --clump_only (lossless reorder)
Input:  <path> (<gzip|plain>, <bytes> bytes)
Output: <path> (<gzip level N | plain>, <bytes> bytes)
Records: <N>
Bins: <n_bins> (peak occupancy <M> records)
Compression ratio: <x.yy>x    [omitted if --dont_gzip or input was plain]
```

The `Compression ratio` line is emitted only when both input and output are gzipped; otherwise the ratio would be misleading (e.g. gzip input → plain output makes the "ratio" look like decompression). The report filename uses `_clumping_report` rather than `_trimming_report` specifically so downstream nf-core/MultiQC pipelines that scan `*_trimming_report.*` don't misclassify it as a trim result.

## When to use this vs. `--clumpify`

Use [`--clumpify`](/performance/clumpy/) when you're **already running a trim pass** and want to squeeze more out of gzip compression at the same time.

Use `--clump_only` when you want **only the clumping win**, without changing your reads at all. Typical scenarios:

- Archival storage of raw FASTQ (bit-for-bit-recoverable when read back, minus the plus-line normalization above)
- Recompression pipeline stages that require lossless reordering
- Testing the clumping mechanism against a fixed input (byte-identity + determinism guarantees make regression detection trivial)

## Performance

`--clump_only` skips adapter auto-detection (which normally scans the first 1 M reads to identify the adapter) and skips the trim/filter passes entirely. Per-record work reduces to `canonical_minimizer` + bin dispatch + a stable sort at flush time. In the ideal case, `--clump_only` at a given `--compression` level should be faster than `--clumpify` at the same level.

Compression ratio depends heavily on data type. See [Clumpy compression](/performance/clumpy/) for the per-data-type guidance — the same clumping mechanism is at work, so the same input types that clump well under `--clumpify` also clump well under `--clump_only`.
