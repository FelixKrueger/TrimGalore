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
| **scRNA-seq (10x Chromium)** | negative — output grows | ❌ **No** — R1 (cell barcode + UMI) reorders cleanly, but R2 (cDNA) follows R1's order to preserve pair lockstep and ends up scrambled vs the natural flowcell-cluster order. R2 disruption beats the R1 win. Use `--compression 6` without `--clumpify` for ~17% saving |
| **Long-read (ONT, PacBio)** | ~0% | ❌ **No** — long reads are mostly unique fragments; clumpify doesn't help and adds wall time |
| **Variable-length / mixed amplicon** | ~0% | ❌ Skip — diversity defeats minimizer clustering |

## How to use it

```bash
# Reorder reads, default gzip level (1 — fastest)
trim_galore --clumpify <input>

# Maximum compression (slowest, smallest output)
trim_galore --clumpify --compression 9 <input>

# Compose for archival storage: max compression with extra memory
trim_galore --clumpify --compression 9 --memory 4G <input>

# Higher gzip without reordering (gzip-only win, no clumping cost)
trim_galore --compression 6 <input>
```

`--clumpify` requires `--cores >= 2` (it feeds the existing parallel worker pool with binned batches) and gzip output (`--dont_gzip` is rejected).

`--compression` is independent: it works with or without `--clumpify`, and applies to the regular trimming pipeline too.

## Performance and compression considerations

### Wall-time cost

Reordering itself is essentially free; the dominant cost at higher gzip levels is gzip's CPU time. Decoupling the two flags lets you pick the trade-off you actually want:

| Mode | Wall time vs plain |
|---|---|
| `--clumpify` (compression 1, default) | ~1.0–1.4× plain (essentially free) |
| `--clumpify --compression 6` | ~1.5–6.4× plain |
| `--clumpify --compression 9` | ~5–10× plain |
| `--compression 6` (no clumpify) | ~1.6–6.4× plain |
| `--compression 9` (no clumpify) | ~5–8× plain |

The minimizer computation uses 2-bit packed integer ops (one O(1) bitwise step per read position) and the per-bin sort is O(n log n) on small bins; both run in the parallel worker pool alongside trim+filter so they overlap with I/O.

### Memory

`--memory` (default `1G`) is a Trim Galore-wide memory budget.

The clumpify dispatcher sizes the per-bin sort runs against it, with the formula picked in an attempt to get the predicted peak RSS ≤ `--memory`:

$$
\begin{aligned}
n_{\text{bins}}          &= \max(16,\; 4 \times \text{cores}) \\[0.8em]
\text{usable}            &= \text{memory} - 512\,\text{MiB}
  \quad\small\text{(FastQC + allocator + runtime overhead)} \\[0.8em]
\text{bin\_byte\_budget} &= \frac{4 \times \text{usable}}{5 \times n_{\text{bins}} + 7 \times \text{cores}}
\end{aligned}
$$

So with `--cores 8 --memory 1G` you get **32 bins × 12 MB**; with `--cores 8 --memory 4G` you get 32 bins × 66 MB. The binary prints those resolved values at startup.

In theory, bigger budget → bigger per-gzip-member sort runs → better compression. In practice, increasing memory doesn't seem to make very much difference in our tests.

#### Below-floor behaviour

`--clumpify` needs a minimum of ~535 MiB to run (mostly from a fixed 512 MiB reservation for FastQC, allocator, and runtime overhead). The exact floor varies slightly with `--cores` but stays in the 535–730 MiB range for any sensible core count.

If `--memory` is below the floor, Trim Galore prints a warning and falls back to plain mode:

```
WARNING: --memory 100M is too small for --clumpify at --cores 6 (need ≥ 552 MiB).
         Falling back to plain mode (no read reordering). Increase --memory or
         drop --clumpify to silence this warning.
```

The trim itself proceeds normally; only the read-reordering step is skipped.
Memory usage without `--clumpify` is typically significantly lower, around the 100MB mark.

### What doesn't change

- Trimming algorithm and per-record output bytes — clumpify only changes the **order** of records.
- All `*_trimming_report.txt` / `.json` numbers are byte-identical between plain and clumpify runs (filter + stats code is order-independent).
- Multi-member gzip is RFC 1952 valid; `zcat`, `seqkit`, `samtools fastq`, and `MultiGzDecoder` all handle it transparently.
- Pair lockstep is preserved: R1[i] and R2[i] are still mates after clumpify reorders them.

## Benchmark results

Real-world numbers from Phil's MacBook Pro (Apple Silicon, 16 GiB RAM, `--cores 6 --memory 1G`, all defaults). Each row shows three configurations vs the same dataset's plain trim:

| Library | Reads | Plain size | `--clumpify` <small>(L1, default)</small> | `--compression 6` <small>(no clumpify)</small> | `--clumpify --compression 6` |
|---|---|---|---|---|---|
| MiSeq amplicon (CRISPR) | 4.4M SE | 500 MB | −16.3% | −25.1% | **−31.7%** |
| ChIP-seq (Illumina SE) | 28.6M SE | 1.5 GB | ~0% | −14.1% | **−15.2%** |
| WES (Illumina SE) | 105M SE | 9.2 GB | −4.4% | −15.4% | **−20.5%** |
| Long-read (ONT/PacBio) | 100K SE | 558 MB | 0% | −8.6% | **−8.8%** |
| ATAC-seq (Illumina PE) | 31.5M PE | 2.9 GB | −34.8% | −21.7% | **−48.5%** |
| Ribo-seq (Illumina PE) | ~30M PE | 4.0 GB | −15.2% | −18.8% | **−32.0%** |
| RNA-seq (Illumina PE) | 93M PE | 17.0 GB | −15.1% | −14.0% | **−27.2%** |
| scRNA-seq 10x Chromium (PE) | 392M PE | 39.6 GB | **+6.2%** ⚠ | **−17.1%** | −12.1% |

For most short-read data, **`--clumpify --compression 6` is the headline configuration** — it saves the most, and on ATAC-seq and Ribo-seq it's also *faster* than `--compression 6` alone (the bin-sorted gzip members compress more efficiently, so deflate finds matches with less hash-chain work). The default `--clumpify` (compression 1) is essentially free in wall time and still gets you 15–35% saving on the most redundant data types.

The huge ATAC-seq and Ribo-seq savings come from the underlying biology: Tn5 transposase has insertion-site biases that produce densely clustered fragments; ribosome footprints are short and stack at start codons. Both produce data where many reads share long substrings, which is exactly what gzip's dictionary loves once those reads are co-located.

The modest WES saving reflects whole-genome diversity: ~4× coverage of 60 MB exome means each genomic position has ~4 forward + ~4 reverse reads, and minimizer clustering only reaches that small group.

**scRNA-seq 10x is the exception.** Cell-barcode R1 is highly redundant and reorders cleanly, but R2 (cDNA) is forced to follow R1's order to preserve pair lockstep, and ends up scrambled vs the natural flowcell-cluster order of the input. At L1 the R2 disruption flips the result to a 6% loss; at L6 the larger gzip dictionary recovers most of the R2 loss but `--clumpify` still costs you ~5 ppt vs `--compression 6` alone. **For 10x scRNA-seq, use `--compression 6` without `--clumpify`.**

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

## Summary — which configuration should I use?

The right configuration depends on **what happens to the trimmed FASTQ after Trim Galore finishes**.

### Pipeline intermediates (deleted after the run): use `--clumpify` alone

If trimmed FASTQ is an ephemeral artefact in a pipeline workdir — typical for nf-core / Snakemake / CWL workflows where the trim step's output feeds an aligner and gets garbage-collected after the pipeline completes — leave gzip compression at the default level 1 and just add `--clumpify`:

```bash
trim_galore --clumpify <input>
```

Wall-time penalty is essentially zero (1.0–1.4× plain on the datasets that benefit), and the smaller intermediate file pays you back in faster I/O for the *next* step (alignment reads less data from disk). Pure win for the dominant pipeline use case.

### Long-term storage / disk-constrained workdirs: add `--compression 6`

If trimmed FASTQ will be retained as a deliverable, archived to object storage, transferred over the network, or kept around because work-dir disk space is tight, the extra gzip CPU at level 6 is worth the run time:

```bash
trim_galore --clumpify --compression 6 <input>
```

This is roughly 4–6× plain wall on most data, but saves 15–50% on output size depending on data type. On the most redundant inputs (ATAC-seq, Ribo-seq) clumpify lets gzip find matches faster than it would on unsorted data, so `--clumpify --compression 6` can actually be *faster* than `--compression 6` alone — a strict win.

For maximum compression at archival cost, `--compression 9` is also available; it adds 2–5 ppt of saving for roughly double the L6 wall time.

### Special cases

- **10x scRNA-seq:** skip `--clumpify`. R2 fragmentation outweighs the R1 sort win at every gzip level on this data; use `--compression 6` without `--clumpify` for a clean ~17% saving.
- **Long-read (ONT / PacBio):** skip `--clumpify`. Long reads are unique fragments — there's nothing for the minimizer sort to cluster.
- **Memory:** the default `--memory 1G` is fine for almost all inputs. Bumping to `2G` adds ~2 ppt of saving on big datasets; going higher has sharply diminishing returns and risks memory pressure on 16 GiB hosts.
