---
title: Clump-only (lossless reorder)
description: Reorder FASTQ or uBAM records for gzip/BGZF-friendly compression without trimming (`--clump_only`).
---

`--clump_only` is a run-and-exit mode that reorders FASTQ or unaligned BAM (uBAM) records by canonical 16-mer minimizer for compression-friendly grouping, **without any trimming, filtering, or adapter detection**. Every input record appears in the output byte-identically (header/name + sequence + quality + preserved aux tags for uBAM); only file-level order changes.

Requested in [#353](https://github.com/FelixKrueger/TrimGalore/issues/353) as an archival and recompression path. It applies the same minimizer-based reordering as [`--clumpify`](/performance/clumpy/), which runs as a compression stage inside the trim pipeline, but operates as a standalone mode that leaves record contents unchanged.

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
- Report: `<input-name>_clumping_report.txt` — a short text summary, named distinctly from the `_trimming_report.txt` glob that downstream nf-core/MultiQC pipelines scan (see [Reorder report](#reorder-report)).

## Record fidelity

The following fields are byte-identical between input and output:

- Read ID (header line for FASTQ; name portion for BAM)
- Sequence
- Quality string
- Aux tags (uBAM only, via `--preserve-tags`) — A/Z/i/f scalars round-trip losslessly. B (array) and H (hex) tags are rejected at BAM-read time by the underlying reader, the same constraint that applies to the trim uBAM path.

Two FASTQ fields are normalized by Trim Galore's writer on every code path, not only under `--clump_only`:

- The plus-line (line 3 of each record) is written as a bare `+`. Inputs carrying `+<header-repeat>` on line 3 emerge with bare `+`. This is inherited from `FastqReader`/`FastqWriter` behaviour.
- CRLF line endings are converted to LF.

One header record is added on uBAM output. Each run appends a `@PG ID:trim_galore VN:<version> CL:<invocation>` record, preserving any input `@PG` chain ahead of it, so the provenance chain grows by one line per invocation. This matches the trim uBAM path.

Byte-identity therefore covers the semantically meaningful record fields, not the whole file. Recompression that must reproduce FASTQ line-3 headers, CRLF endings, or a BAM header without the added `@PG` record is outside what `--clump_only` provides.

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

Accepted but inert: `-q` / `--quality`, `--stringency`, `-e` / `--error`. These are retained in the CLI surface for parity with the trim path and have no effect under this mode. They are listed here because silent acceptance is otherwise easy to mistake for an oversight.

## Reorder report

The reorder report at `<stem>_clumping_report.txt` is narrower than the standard trimming report by design — no adapter counts, no quality statistics, no filter counters, and no JSON companion.

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

## Relationship to `--clumpify`

[`--clumpify`](/performance/clumpy/) is the right choice when a trim pass is being run anyway and the trimmed output should also compress more tightly. `--clump_only` applies the reordering on its own, leaving record contents untouched. It suits:

- Archival storage of raw FASTQ, recoverable on read-back apart from the plus-line and line-ending normalization described above
- Archival storage of aux-tag-carrying uBAM, such as 10X Chromium single-cell libraries where the `CB`, `UB`, and `BC` tags are meaningful — `--clump_only --output-format ubam --preserve-tags CB,UB,BC` reorders records for BGZF-friendly grouping while round-tripping the tags losslessly
- Recompression pipeline stages that require lossless reordering
- Exercising the clumping mechanism against a fixed input, where the byte-identity and determinism guarantees make regressions straightforward to detect

## Performance

`--clump_only` skips adapter auto-detection, which otherwise scans the first 1 M reads, and skips the trim and filter passes entirely. Per-record work reduces to `canonical_minimizer`, bin dispatch, and a stable sort at flush time. At a given `--compression` level it therefore performs strictly less work per record than `--clumpify`.

Compression ratio depends heavily on data type. See [Clumpy compression](/performance/clumpy/) for the per-data-type guidance — the same clumping mechanism is at work, so the same input types that clump well under `--clumpify` also clump well under `--clump_only`.
