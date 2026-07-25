---
title: Clump-only (lossless reorder)
description: Reorder FASTQ or uBAM records for gzip/BGZF-friendly compression without trimming (`--clump_only`).
---

`--clump_only` is a run-and-exit mode that reorders FASTQ or unaligned BAM (uBAM) records by canonical 16-mer minimizer for compression-friendly grouping, **without any trimming, filtering, or adapter detection**. Every input record appears in the output byte-identically (header/name + sequence + quality + preserved aux tags for uBAM); only file-level order changes.

Requested in [#353](https://github.com/FelixKrueger/TrimGalore/issues/353) as an archival/recompression path — same clumping win that [`--clumpify`](/performance/clumpy/) delivers as a compression stage inside the trim pipeline, but as a standalone lossless mode you can use to recompress raw reads for storage without ever mutating them.

## Usage

**FASTQ → FASTQ (default):**

```bash
# Single-end
trim_galore --clump_only --compression 9 --cores 2 sample.fq.gz

# Paired-end (two files)
trim_galore --clump_only --paired --compression 9 --cores 2 sample_R1.fq.gz sample_R2.fq.gz
```

**FASTQ → uBAM or uBAM → uBAM (v2):**

```bash
# Single-end, any input format
trim_galore --clump_only --output-format ubam --cores 2 sample.fq.gz
trim_galore --clump_only --output-format ubam --cores 2 sample.bam

# Paired-end, two FASTQ files → one interleaved BAM per pair
trim_galore --clump_only --paired --output-format ubam sample_R1.fq.gz sample_R2.fq.gz

# Paired-end, single interleaved uBAM → single interleaved uBAM
trim_galore --clump_only --paired --output-format ubam interleaved.bam

# Preserve aux tags through the reorder (uBAM only)
trim_galore --clump_only --output-format ubam --preserve-tags RG,BC,UB sample.bam
```

**Output filenames:**

- SE FASTQ: `sample_clumped.fq.gz` (or `.fq` under `--dont_gzip`)
- PE FASTQ: `sample_R1_clumped_1.fq.gz` and `sample_R2_clumped_2.fq.gz`
- SE uBAM: `sample_clumped.bam` (BGZF; `--dont_gzip` rejected)
- PE uBAM: **ONE** interleaved `sample_R1_clumped.bam` per pair (mate-adjacent, matching samtools/Picard/fgbio convention)
- Report: `<input-name>_clumping_report.txt` — a short text summary, deliberately distinct from the `_trimming_report.txt` glob that downstream nf-core/MultiQC pipelines scan (see the [Report shape](#report-shape) section).

## What's preserved and what isn't

**Byte-identical across the reorder:**

- Read ID (header line for FASTQ; name portion for BAM)
- Sequence
- Quality string
- **Aux tags** (uBAM only, via `--preserve-tags`) — A/Z/i/f scalars round-trip losslessly. B (array) and H (hex) tags are rejected at BAM-read time by the underlying reader; this is the same constraint as the trim uBAM path.

**Normalized codebase-wide (not `--clump_only`-specific), FASTQ only:**

- **Plus-line (line 3 of each record)** is normalized to bare `+` on output. Inputs with `+<header-repeat>` on line 3 emerge with bare `+`. This applies to every Trim Galore FASTQ code path, not just `--clump_only`, and is inherited from the `FastqReader`/`FastqWriter` behaviour.
- **CRLF line endings** are normalized to LF. Inputs with `\r\n` endings emerge with `\n` endings.

**Added on uBAM output** (not carried from input):

- **`@PG` line** — every uBAM run appends a new `@PG ID:trim_galore VN:<version> CL:<invocation>` record to the header, preserving the input `@PG` chain if any. Same behavior as the trim uBAM path. Provenance chain grows by one line per invocation.

The byte-identity guarantee covers the semantically-meaningful record fields. If you need byte-for-byte-exact recompression including FASTQ line-3 headers, CRLF endings, or BAM headers stripped of `@PG`, `--clump_only` is not the right tool.

## Determinism

Two invocations with the same input, `--cores`, `--memory`, and `--compression` produce **byte-identical FASTQ output files** and **record-body-identical BAM output files** (BAM whole-file hashes vary because the `@PG` line's `CL:` field varies with the invocation string; CI comparison uses a `@PG`-ignoring assertion, same treatment as the trim uBAM path). Downstream tools that md5 FASTQ output for archival integrity get stable hashes; tools that md5 BAM should compare via `samtools view -H | grep -v @PG` piped diff (or noodles-based tuple comparison). Both properties are enforced by CI.

## Compatibility

Composes with:

- `--paired` — two-file paired input (Shape A, multi-pair `N=4, 6, …` supported) OR single interleaved uBAM (Shape B, N=1 + `--output-format ubam`)
- `--compression <1-9>` (FASTQ gzip level; ignored for uBAM output — BAM is always BGZF)
- `--memory <SIZE>` (bin-buffer sizing, e.g. `4G`)
- `--cores <N>` (accepted for interface parity — v1 is single-threaded internally; parallelism is planned for a later release)
- `--fastqc` — runs on the reordered output; fastqc-rust reads both FASTQ and BAM natively
- `--dont_gzip` — FASTQ output only, produces `*_clumped.fq` (plain). Rejected with `--output-format ubam` (BAM is always BGZF-compressed).
- `--basename BASE` — output becomes `BASE_clumped.fq(.gz)` (SE FASTQ), `BASE_clumped_{1,2}.fq(.gz)` (PE FASTQ), or `BASE_clumped.bam` (uBAM)
- `--output-format ubam` — produces uBAM output. Composes with all of the above except `--dont_gzip`.
- `--preserve-tags TAG1,TAG2,…` — uBAM only. Aux tags round-trip through the reorder. A/Z/i/f scalars supported; B (array) and H (hex) rejected at BAM-read time.

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
- `--passthrough`, `--retain_unpaired`
- `--dont_gzip` + `--output-format ubam` (BAM is always BGZF; enforced at shared `Cli::validate` §3.4a — this also closes an existing gap on the trim uBAM path)

**PE input-shape rejections (uBAM output only):**

- Two BAM files under `--paired`: `--paired R1.bam R2.bam` is ambiguous. Use a single interleaved uBAM instead (matches samtools/Picard/fgbio convention).
- N=1 with a FASTQ file under `--paired --output-format ubam`: single-file paired mode requires a uBAM interleaved input.
- Mixed formats in the same pair: `--paired R1.fq.gz R2.bam` bails; both inputs of a pair must share format.

**FASTQ output path with uBAM input:** `--clump_only sample.bam` (without `--output-format ubam`) is rejected — routing uBAM through the FASTQ output would drop aux tags. Add `--output-format ubam` to keep tags intact.

**Silently ignored** (kept in the CLI surface for parity with the trim path, but do nothing under this mode; documented so users don't wonder why they weren't rejected): `-q` / `--quality`, `--stringency`, `-e` / `--error`.

## Report shape

The reorder report at `<stem>_clumping_report.txt` is deliberately narrower than the standard trimming report — no adapter counts, no quality stats, no filter counters, no JSON companion.

```
Trim Galore version: 2.x.x
Mode: --clump_only (lossless reorder)
Input:  <path> (<format>, <bytes> bytes)      # e.g. "FASTQ (gzip)" or "uBAM"
Output: <path> (<format>, <bytes> bytes)      # e.g. "gzip level 6" or "uBAM (BGZF)"
Records: <N>
Bins: <n_bins> (peak occupancy <M> records)
Compression ratio: <x.yy>x                    # omitted unless BOTH sides are compressed
Preserved tags: TAG1,TAG2,…                   # uBAM only; omitted if empty
```

The `Compression ratio` line is emitted only when both input and output are compressed (any of gzip / BGZF). Under `--dont_gzip` or when input was plain FASTQ, the ratio would be misleading and is omitted. The report filename uses `_clumping_report` rather than `_trimming_report` specifically so downstream nf-core/MultiQC pipelines that scan `*_trimming_report.*` don't misclassify it as a trim result. Multi-pair PE runs produce one report per pair, matching v1 SE's one-report-per-input convention.

## When to use this vs. `--clumpify`

Use [`--clumpify`](/performance/clumpy/) when you're **already running a trim pass** and want to squeeze more out of gzip compression at the same time.

Use `--clump_only` when you want **only the clumping win**, without changing your reads at all. Typical scenarios:

- Archival storage of raw FASTQ (bit-for-bit-recoverable when read back, minus the plus-line normalization above)
- Archival storage of aux-tag-carrying uBAM (e.g. 10X Chromium single-cell libraries where `CB`, `UB`, `BC` tags matter) — `--clump_only --output-format ubam --preserve-tags CB,UB,BC` reorders records for BGZF-friendly compression while round-tripping the tags losslessly
- Recompression pipeline stages that require lossless reordering
- Testing the clumping mechanism against a fixed input (byte-identity + determinism guarantees make regression detection trivial)

## Performance

`--clump_only` skips adapter auto-detection (which normally scans the first 1 M reads to identify the adapter) and skips the trim/filter passes entirely. Per-record work reduces to `canonical_minimizer` + bin dispatch + a stable sort at flush time. In the ideal case, `--clump_only` at a given `--compression` level should be faster than `--clumpify` at the same level.

Compression ratio depends heavily on data type. See [Clumpy compression](/performance/clumpy/) for the per-data-type guidance — the same clumping mechanism is at work, so the same input types that clump well under `--clumpify` also clump well under `--clump_only`.
