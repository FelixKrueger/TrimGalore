---
title: Clumpy compression
description: Reorder reads in the trimmed FASTQ output to make .fq.gz files significantly smaller — opt-in via --clumpify.
---

`--clumpify` is an opt-in mode that reorders reads inside each gzip member of the trimmed output so reads sharing a canonical 16-mer minimizer land adjacent on disk. gzip's 32 KB sliding window then finds long redundant runs of similar sequences, shrinking the `.fq.gz` by **15–55%** depending on the data type, with no information loss — only the on-disk order of records changes.

`--compression <N>` sets the gzip level independently (1–9, default 1). Combine the two for maximum effect: `--clumpify --compression 9` reorders reads **and** runs gzip at its slowest/smallest level.

## When to use it

| Data type | Typical saving (`--clumpify --compression 9`) | Recommendation |
|---|---|---|
| **ATAC-seq (paired)** | ~50% | ✅ **Strong yes** — Tn5 insertion bias creates very high fragment redundancy |
| **Ribo-seq (paired)** | ~45% | ✅ **Strong yes** — short ribosome-protected fragments are highly clustered |
| **MiSeq amplicon / CRISPR / sgRNA** | 30–37% | ✅ **Strong yes** — explicit amplification produces lots of duplicates |
| **Bisulfite-seq / RRBS (paired)** | ~25–35% | ✅ Yes — restriction-enzyme + bisulfite produces a redundant read fingerprint |
| **ChIP-seq (single-end)** | ~24% | ✅ Yes — peaks generate clustered reads |
| **RNA-seq (paired)** | 16–30% | ✅ Yes — highly-expressed transcripts create dense clusters; bigger savings at higher gzip levels |
| **WES / WGS (paired)** | 6–22% | 🟡 Modest — diverse coverage gives less clustering |
| **scRNA-seq (10x Chromium)** | varies | ✅ Yes for cell-barcode runs (R1 is short barcode + UMI; high R1 redundancy) |
| **Long-read (ONT, PacBio)** | ~0% | ❌ **No** — long reads are mostly unique fragments; clumpify doesn't help and adds wall time |
| **Variable-length / mixed amplicon** | ~0% | ❌ Skip — diversity defeats minimizer clustering |

## How to use it

```bash
# Reorder reads, default gzip level (1 — fastest)
trim_galore --clumpify <input>

# Maximum compression (slowest, smallest output)
trim_galore --clumpify --compression 9 <input>

# Compose for archival storage: max compression with full memory budget
trim_galore --clumpify --compression 9 --memory 16G <input>

# Higher gzip without reordering (gzip-only win, no clumping cost)
trim_galore --compression 6 <input>
```

`--clumpify` requires `--cores >= 2` (it feeds the existing parallel worker pool with binned batches) and gzip output (`--dont_gzip` is rejected). `--compression` is independent — it works with or without `--clumpify`, and applies to the regular trimming pipeline too.

## Performance and compression considerations

### Wall-time cost

Reordering itself is essentially free; the dominant cost at higher gzip levels is gzip's CPU time. Decoupling the two flags lets you pick the trade-off you actually want:

| Mode | Wall time vs plain |
|---|---|
| `--clumpify` (compression 1) | ~2× plain |
| `--clumpify --compression 6` | ~3–6× plain |
| `--clumpify --compression 9` | ~5–10× plain |
| `--compression 6` (no clumpify) | ~3–4× plain |
| `--compression 9` (no clumpify) | ~5–8× plain |

The minimizer computation uses 2-bit packed integer ops (one O(1) bitwise step per read position) and the per-bin sort is O(n log n) on small bins; both run in the parallel worker pool alongside trim+filter so they overlap with I/O.

### Memory

`--memory` (default `4G`) is a Trim Galore-wide memory budget. The clumpify dispatcher sizes the per-bin sort runs against it, with the formula picked so **predicted peak RSS ≤ `--memory`**:

```
n_bins          = max(16, 4 × cores)
usable          = memory − 512 MiB        (FastQC + allocator + runtime overhead)
bin_byte_budget = 4 × usable / (5 × n_bins + 7 × cores)
```

So with `--cores 8 --memory 4G` you get **32 bins × 66 MB**, predicted peak ≈ 4 GB — and the binary prints those resolved values at startup. On a 16 GiB ATAC-seq run this matches the actual `/usr/bin/time -l` "peak memory footprint" within 1%. Bigger budget → bigger per-gzip-member sort runs → better compression, asymptoting once the budget exceeds the uncompressed input size.

For the same input:

| `--memory` | Bin size | Saving (MiSeq amplicon `--compression 9`) |
|---|---|---|
| 1G | ~9 MB | 30% |
| 4G (default) | ~66 MB | 36% |
| 16G | ~290 MB | 37% |

The default `4G` is the sweet spot for typical ~1–10 GB inputs. For very large inputs (multi-GB compressed RNA-seq, 10x scRNA-seq), bumping to `8G` or `16G` adds a few percentage points if you have the headroom.

### What doesn't change

- Trimming algorithm and per-record output bytes — clumpify only changes the **order** of records.
- All `*_trimming_report.txt` / `.json` numbers are byte-identical between plain and clumpify runs (filter + stats code is order-independent).
- Multi-member gzip is RFC 1952 valid; `zcat`, `seqkit`, `samtools fastq`, and `MultiGzDecoder` all handle it transparently.
- Pair lockstep is preserved: R1[i] and R2[i] are still mates after clumpify reorders them.

### What can break

If a downstream tool walks reads in input order **and** depends on that order matching some external file (e.g. a separate barcode whitelist file with positional mapping), clumpify will break it. The trimmed FASTQ as a multiset of pairs is unchanged, but the line-by-line order is not.

For scRNA-seq with 10x Chromium reads, this is fine — the cell barcode + UMI are encoded in R1's sequence, and downstream tools (Cell Ranger, STARsolo, etc.) parse the barcode from each read independently. They don't depend on read order.

## Benchmark results

Real-world numbers from Phil's MacBook Pro (Apple Silicon, `--cores 8`, `--memory 4G`):

| Library | Reads | Plain size | `--clumpify --compression 6` saving | `--clumpify --compression 9` saving |
|---|---|---|---|---|
| MiSeq amplicon (CRISPR) | 4.4M SE | 500 MB | 16% | 35% |
| ChIP-seq (Illumina SE) | 28.6M SE | 1.5 GB | 24% | 27% |
| ATAC-seq (Illumina PE) | ~25M PE | 2.9 GB | 51% | 53% |
| Ribo-seq (Illumina PE) | ~30M PE | 4.0 GB | 44% | 46% |
| RNA-seq (Illumina PE) | 93M PE | 17.0 GB | 17% | 30% |
| scRNA-seq 10x Chromium (PE) | 392M PE | 39.6 GB | 13% | — |
| WES (Illumina SE) | 105M SE | 9.2 GB | 6% | 22% |
| Long-read (ONT/PacBio) | 100K SE | 558 MB | 0% | 0% |

The huge ATAC-seq and Ribo-seq savings come from the underlying biology: Tn5 transposase has insertion-site biases that produce densely clustered fragments; ribosome footprints are short and stack at start codons. Both produce data where many reads share long substrings, which is exactly what gzip's dictionary loves once those reads are co-located.

The modest WES saving reflects whole-genome diversity: ~4× coverage of 60 MB exome means each genomic position has ~4 forward + ~4 reverse reads, and minimizer clustering only reaches that small group.

The 0% long-read result holds the disclaimer line: ONT and PacBio reads are mostly unique multi-kb fragments, so there's nothing for clumpify to cluster.

## How it works (briefly)

For each (paired-end pair of) trimmed reads, the dispatcher:

1. Computes the **canonical 16-mer minimizer** of R1 — the lexicographically smallest of `min(forward, revcomp)` over all 16-mer windows in the read. R2 follows R1's bin assignment in lockstep so paired reads stay together.
2. Hashes the minimizer (FNV-1a) into one of `n_bins` per-bin in-memory buffers.
3. When a bin's accumulated raw bytes exceed `bin_byte_budget`, the bin is sorted by minimizer key (with sequence/quality/id tie-breaks for determinism) and shipped to a worker thread for trimming + gzip compression. Each worker emits one gzip member per batch.
4. The main thread reassembles compressed bytes in flush order via the existing BTreeMap-keyed queue, identical to the non-clumpify parallel pipeline.

Each gzip member therefore contains a sorted run of records sharing a minimizer hash bucket, which is what gzip's 32 KB sliding window exploits to find long back-references.

## Comparison with other tools

If you've used `bbmap clumpify` or [`stevekm/squish`](https://github.com/stevekm/squish), `--clumpify` produces compression results in the same ballpark on most data:

- On amplicon-type data, all three tools converge to ~37–38% saving at gzip level 9.
- On diverse data (WES, WGS), all three give ~20–30% saving — none of them works miracles on inherently diverse libraries.

The advantage of `--clumpify` over running a separate tool is that it's part of the trim pass, so there's no extra read-and-rewrite cycle. For an X GB input, a separate clumpify step would mean reading the trimmed output, sorting, and writing it back — typically +3–5× the trim wall time. `--clumpify` does it in one pass.

## Summary

- Default to `--clumpify --compression 6` for most short-read library types — modest extra wall time, meaningful compression.
- Use `--clumpify --compression 9` when storage cost dominates (archival, S3 transfer, tape).
- `--compression` works without `--clumpify` too — handy when you want a smaller output without paying for the reorder.
- Skip `--clumpify` for long-read data — it's a no-op there with non-trivial wall-time cost.
- Bump `--memory` (e.g. `--memory 16G`) on big inputs (multi-GB compressed) for the last few percentage points.
