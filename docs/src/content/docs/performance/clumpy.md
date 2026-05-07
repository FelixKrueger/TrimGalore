---
title: Clumpy compression
description: Reorder reads in the trimmed FASTQ output to make .fq.gz files significantly smaller — opt-in via --clumpify.
---

`--clumpify` is an opt-in mode that reorders reads inside each gzip member of the trimmed output so reads sharing a canonical 16-mer minimizer land adjacent on disk. gzip's 32 KB sliding window then finds long redundant runs of similar sequences, shrinking the `.fq.gz` by **15–55%** depending on the data type, with no information loss — only the on-disk order of records changes.

`--compression <N>` sets the gzip level independently (1–9, default 1). Combine the two for maximum effect: `--clumpify --compression 9` reorders reads **and** runs gzip at its slowest/smallest level.

## When to use it

### Clumpify

- ✅ Low complexity data: yes (ATAC-seq, ChIP-seq, Ribo-Seq, RNA-Seq, high sequencing depth WES)
- ❌ High complexity data: no (whole-genome sequencing) - has no effect
- ❌ Long reads: no (Oxford Nanopore) - has no effect
- ❌ Unusual paired-end formats: no - can have deleterious effect

### Compression

Whether to push up `--compression` level or not depends on what the trimmed FASTQ is used for:

- **Pipeline intermediates** (trimmed FASTQ is ephemeral, deleted after the pipeline finishes)
    - Leave as compression level 1, but can still use `--clumpify`.
    - The reorder is essentially free (1.0–1.4× slowdown on most data) and the smaller output makes the *next* step (typically an aligner) read less from disk — net I/O win for the whole pipeline.
- **Long-term storage or disk-constrained workdirs**
    - Add `--compression 6` (or `--compression 9` for archival)
    - `--clumpify --compression 6` can halve output file sizes (15–50% less) but makes the run time 4–6× slower.
    - Can be specified without `--clumpify`, but with redundant data types (ATAC-seq, Ribo-seq) it typically runs *faster* than with clumpify on, because deflate finds matches more cheaply on sorted runs.

### Data types

| Data type | Typical saving (`--clumpify --compression 9`) | Recommendation |
|---|---|---|
| **ATAC-seq (paired)** | ~50% | ✅ **Strong yes**: Tn5 insertion bias creates very high fragment redundancy |
| **Ribo-seq (paired)** | ~45% | ✅ **Strong yes**: short ribosome-protected fragments are highly clustered |
| **MiSeq amplicon / CRISPR / sgRNA** | 30–37% | ✅ **Strong yes**: explicit amplification produces lots of duplicates |
| **Bisulfite-seq / RRBS (paired)** | ~25–35% | ✅ **Yes**: restriction-enzyme + bisulfite produces a redundant read fingerprint |
| **ChIP-seq (single-end)** | ~24% | ✅ **Yes**: peaks generate clustered reads |
| **RNA-seq (paired)** | 16–30% | ✅ **Yes**: highly-expressed transcripts create dense clusters; bigger savings at higher gzip levels |
| **WES / WGS (paired)** | 6–22% | 🟡 **Modest**: diverse coverage gives less clustering |
| **scRNA-seq (10x Chromium)** | negative — output grows | ❌ **No**: R1 (cell barcode + UMI) reorders cleanly, but R2 (cDNA) follows R1's order to preserve pair lockstep and ends up scrambled vs the natural flowcell-cluster order. R2 disruption beats the R1 win. Use `--compression 6` without `--clumpify` for ~17% saving |
| **Long-read (ONT, PacBio)** | ~0% | ❌ **No**: long reads are mostly unique fragments; clumpify doesn't help and adds wall time |
| **Variable-length / mixed amplicon** | ~0% | ❌ **Skip**: diversity defeats minimizer clustering |

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

### Downstream BAM compression

The read-clustering effect carries through into downstream **unsorted** BAM files at essentially full strength.

Because aligners typically stream output reads out in the same order that they come in with, and BAM files useg gzip compression internally, the same clumping improvements hold true through alignment.

Here's an example using ATAC-seq data (31 M paired-end reads), aligning to a minimal index (chr22 only):

| BAM stage | Saving (clumpify vs plain) |
|---|---|
| Trimmed FASTQ (gzip level 1) | −34.8% |
| `samtools import` → uBAM (no alignment) | −36.4% |
| STAR 2.7.11b chr22 alignment → unsorted, aligned BAM | −34.2% |

Note that only unsorted BAMs benefit. If your pipeline coordinate-sorts the BAM immediately after alignment (e.g. `samtools sort` or `STAR --outSAMtype BAM SortedByCoordinate`), the read order is rearranged by genomic position and the input-order signal is erased. A coordinate-sorted BAM's size is determined by genomic distribution of reads, not by clumpify's clustering.

## Benchmark results

Real-world numbers from a benchmark using a MacBook Pro (Apple Silicon, 16 GiB RAM, `--cores 6 --memory 1G`, all defaults). Each dataset has 3 bars:

- `--clumpify`(L1, default)
- `--compression 6` (no clumpify)
- `--clumpify --compression 6`.

Plots show compression savings (how much smaller the resulting FastQ files are versus the regular run) and the Wall-time effect (how much slower the run was, 1x is original run).

<svg viewBox="0 0 1100 460" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system, BlinkMacSystemFont, sans-serif" fill="currentColor" role="img" aria-label="Compression saving vs plain, percentage. Bar chart with 8 datasets, 3 configurations each.">
  <text x="550" y="28" text-anchor="middle" font-weight="600" font-size="22">Compression saving vs plain (%)</text>
  <text x="550" y="50" text-anchor="middle" font-size="15" fill-opacity="0.65">positive = output smaller; negative = output grew</text>
  <line x1="70" y1="60" x2="70" y2="370" stroke="currentColor" stroke-opacity="0.5"/>
  <line x1="70" y1="370" x2="1040" y2="370" stroke="currentColor" stroke-opacity="0.5"/>
  <line x1="70" y1="318.33" x2="1040" y2="318.33" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="2 3"/>
  <text x="62" y="65" text-anchor="end" font-size="14" fill-opacity="0.75">50%</text>
  <text x="62" y="117" text-anchor="end" font-size="14" fill-opacity="0.75">40%</text>
  <text x="62" y="169" text-anchor="end" font-size="14" fill-opacity="0.75">30%</text>
  <text x="62" y="220" text-anchor="end" font-size="14" fill-opacity="0.75">20%</text>
  <text x="62" y="272" text-anchor="end" font-size="14" fill-opacity="0.75">10%</text>
  <text x="62" y="323" text-anchor="end" font-size="14" fill-opacity="0.75">0%</text>
  <text x="62" y="375" text-anchor="end" font-size="14" fill-opacity="0.75">−10%</text>
  <g font-size="19">
    <text x="130" y="385" text-anchor="end" transform="rotate(-25 130 385)">MiSeq</text>
    <text x="250" y="385" text-anchor="end" transform="rotate(-25 250 385)">ChIP-seq</text>
    <text x="370" y="385" text-anchor="end" transform="rotate(-25 370 385)">WES</text>
    <text x="490" y="385" text-anchor="end" transform="rotate(-25 490 385)">Long-read</text>
    <text x="610" y="385" text-anchor="end" transform="rotate(-25 610 385)">ATAC-seq</text>
    <text x="730" y="385" text-anchor="end" transform="rotate(-25 730 385)">Ribo-seq</text>
    <text x="850" y="385" text-anchor="end" transform="rotate(-25 850 385)">RNA-seq</text>
    <text x="970" y="385" text-anchor="end" transform="rotate(-25 970 385)">scRNA-seq 10x</text>
  </g>
  <rect x="87" y="234.10" width="26" height="84.23" fill="#2a4d8f"/>
  <rect x="117" y="188.65" width="26" height="129.68" fill="#c25400"/>
  <rect x="147" y="154.55" width="26" height="163.78" fill="#1a7f1a"/>
  <rect x="207" y="318.33" width="26" height="0.41" fill="#2a4d8f"/>
  <rect x="237" y="245.48" width="26" height="72.85" fill="#c25400"/>
  <rect x="267" y="239.80" width="26" height="78.53" fill="#1a7f1a"/>
  <rect x="327" y="295.60" width="26" height="22.73" fill="#2a4d8f"/>
  <rect x="357" y="238.76" width="26" height="79.57" fill="#c25400"/>
  <rect x="387" y="212.41" width="26" height="105.92" fill="#1a7f1a"/>
  <rect x="447" y="318.33" width="26" height="0.5" fill="#2a4d8f"/>
  <rect x="477" y="273.90" width="26" height="44.43" fill="#c25400"/>
  <rect x="507" y="272.86" width="26" height="45.47" fill="#1a7f1a"/>
  <rect x="567" y="138.53" width="26" height="179.80" fill="#2a4d8f"/>
  <rect x="597" y="206.21" width="26" height="112.12" fill="#c25400"/>
  <rect x="627" y="67.75" width="26" height="250.58" fill="#1a7f1a"/>
  <rect x="687" y="239.80" width="26" height="78.53" fill="#2a4d8f"/>
  <rect x="717" y="221.20" width="26" height="97.13" fill="#c25400"/>
  <rect x="747" y="152.99" width="26" height="165.34" fill="#1a7f1a"/>
  <rect x="807" y="240.31" width="26" height="78.02" fill="#2a4d8f"/>
  <rect x="837" y="245.99" width="26" height="72.34" fill="#c25400"/>
  <rect x="867" y="177.80" width="26" height="140.53" fill="#1a7f1a"/>
  <rect x="927" y="318.33" width="26" height="32.03" fill="#2a4d8f"/>
  <rect x="957" y="229.97" width="26" height="88.36" fill="#c25400"/>
  <rect x="987" y="255.83" width="26" height="62.50" fill="#1a7f1a"/>
  <g transform="translate(820, 80)" font-size="15">
    <rect x="0" y="0" width="18" height="14" fill="#2a4d8f"/>
    <text x="26" y="12">--clumpify (L1, default)</text>
    <rect x="0" y="20" width="18" height="14" fill="#c25400"/>
    <text x="26" y="32">--compression 6</text>
    <rect x="0" y="40" width="18" height="14" fill="#1a7f1a"/>
    <text x="26" y="52">--clumpify --compression 6</text>
  </g>
</svg>

<svg viewBox="0 0 1100 460" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system, BlinkMacSystemFont, sans-serif" fill="currentColor" role="img" aria-label="Wall time vs plain, multiplier. Bar chart with 8 datasets, 3 configurations each.">
  <text x="550" y="28" text-anchor="middle" font-weight="600" font-size="22">Wall time vs plain (×)</text>
  <text x="550" y="50" text-anchor="middle" font-size="15" fill-opacity="0.65">how much longer than the plain run took</text>
  <line x1="70" y1="60" x2="70" y2="370" stroke="currentColor" stroke-opacity="0.5"/>
  <line x1="70" y1="370" x2="1040" y2="370" stroke="currentColor" stroke-opacity="0.5"/>
  <line x1="70" y1="347.86" x2="1040" y2="347.86" stroke="currentColor" stroke-opacity="0.4" stroke-dasharray="3 3"/>
  <text x="1040" y="340" text-anchor="end" font-size="13" fill-opacity="0.65">1× (plain wall)</text>
  <text x="62" y="375" text-anchor="end" font-size="14" fill-opacity="0.75">0×</text>
  <text x="62" y="331" text-anchor="end" font-size="14" fill-opacity="0.75">2×</text>
  <text x="62" y="287" text-anchor="end" font-size="14" fill-opacity="0.75">4×</text>
  <text x="62" y="243" text-anchor="end" font-size="14" fill-opacity="0.75">6×</text>
  <text x="62" y="199" text-anchor="end" font-size="14" fill-opacity="0.75">8×</text>
  <text x="62" y="155" text-anchor="end" font-size="14" fill-opacity="0.75">10×</text>
  <text x="62" y="110" text-anchor="end" font-size="14" fill-opacity="0.75">12×</text>
  <text x="62" y="66" text-anchor="end" font-size="14" fill-opacity="0.75">14×</text>
  <g font-size="19">
    <text x="130" y="385" text-anchor="end" transform="rotate(-25 130 385)">MiSeq</text>
    <text x="250" y="385" text-anchor="end" transform="rotate(-25 250 385)">ChIP-seq</text>
    <text x="370" y="385" text-anchor="end" transform="rotate(-25 370 385)">WES</text>
    <text x="490" y="385" text-anchor="end" transform="rotate(-25 490 385)">Long-read</text>
    <text x="610" y="385" text-anchor="end" transform="rotate(-25 610 385)">ATAC-seq</text>
    <text x="730" y="385" text-anchor="end" transform="rotate(-25 730 385)">Ribo-seq</text>
    <text x="850" y="385" text-anchor="end" transform="rotate(-25 850 385)">RNA-seq</text>
    <text x="970" y="385" text-anchor="end" transform="rotate(-25 970 385)">scRNA-seq 10x</text>
  </g>
  <rect x="87" y="347.41" width="26" height="22.59" fill="#2a4d8f"/>
  <rect x="117" y="333.69" width="26" height="36.31" fill="#c25400"/>
  <rect x="147" y="336.57" width="26" height="33.43" fill="#1a7f1a"/>
  <rect x="207" y="339.00" width="26" height="31.00" fill="#2a4d8f"/>
  <rect x="237" y="261.29" width="26" height="108.71" fill="#c25400"/>
  <rect x="267" y="264.40" width="26" height="105.60" fill="#1a7f1a"/>
  <rect x="327" y="345.43" width="26" height="24.57" fill="#2a4d8f"/>
  <rect x="357" y="228.74" width="26" height="141.26" fill="#c25400"/>
  <rect x="387" y="232.28" width="26" height="137.72" fill="#1a7f1a"/>
  <rect x="447" y="346.31" width="26" height="23.69" fill="#2a4d8f"/>
  <rect x="477" y="296.93" width="26" height="73.07" fill="#c25400"/>
  <rect x="507" y="301.14" width="26" height="68.86" fill="#1a7f1a"/>
  <rect x="567" y="344.98" width="26" height="25.02" fill="#2a4d8f"/>
  <rect x="597" y="278.11" width="26" height="91.89" fill="#c25400"/>
  <rect x="627" y="317.52" width="26" height="52.48" fill="#1a7f1a"/>
  <rect x="687" y="346.31" width="26" height="23.69" fill="#2a4d8f"/>
  <rect x="717" y="253.32" width="26" height="116.68" fill="#c25400"/>
  <rect x="747" y="284.97" width="26" height="85.03" fill="#1a7f1a"/>
  <rect x="807" y="346.75" width="26" height="23.25" fill="#2a4d8f"/>
  <rect x="837" y="239.36" width="26" height="130.64" fill="#c25400"/>
  <rect x="867" y="258.84" width="26" height="111.16" fill="#1a7f1a"/>
  <rect x="927" y="342.77" width="26" height="27.23" fill="#2a4d8f"/>
  <rect x="957" y="236.26" width="26" height="133.74" fill="#c25400"/>
  <rect x="987" y="229.00" width="26" height="141.00" fill="#1a7f1a"/>
  <g transform="translate(820, 80)" font-size="15">
    <rect x="0" y="0" width="18" height="14" fill="#2a4d8f"/>
    <text x="26" y="12">--clumpify (L1, default)</text>
    <rect x="0" y="20" width="18" height="14" fill="#c25400"/>
    <text x="26" y="32">--compression 6</text>
    <rect x="0" y="40" width="18" height="14" fill="#1a7f1a"/>
    <text x="26" y="52">--clumpify --compression 6</text>
  </g>
</svg>

Datasets covered:

- **MiSeq amplicon (CRISPR)**: 4.4M SE, 500 MB plain output — [`ERR16944282`](https://www.ncbi.nlm.nih.gov/sra/?term=ERR16944282)
- **ChIP-seq (Illumina SE)**: 28.6M SE, 1.5 GB — [`SRX747791`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX747791)
- **WES (Illumina SE)**: 105M SE, 9.2 GB — [`SRR7890918`](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7890918)
- **Long-read (ONT)**: 100K SE, 558 MB — [`SRR37915503`](https://www.ncbi.nlm.nih.gov/sra/?term=SRR37915503)
- **ATAC-seq (Illumina PE)**: 31.5M PE, 2.9 GB — [`SRX2717909`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX2717909)
- **Ribo-seq (Illumina PE)**: ~30M PE, 4.0 GB — [`SRX11780879`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX11780879) ([`SRR15480782`](https://www.ncbi.nlm.nih.gov/sra/?term=SRR15480782))
- **RNA-seq (Illumina PE)**: 93M PE, 17.0 GB — [`SRX1603629`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX1603629)
- **scRNA-seq 10x Chromium (PE)**: 392M PE, 39.6 GB — [`pbmc8k_v2`](https://www.10xgenomics.com/datasets/8-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0) (10x Genomics public dataset)

## Comparison with other tools

If you've used [`bbmap clumpify`](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) or [`stevekm/squish`](https://github.com/stevekm/squish), `--clumpify` produces compression results in the same ballpark on most data:

- On amplicon-type data, all three tools converge to ~37–38% saving at gzip level 9.
- On diverse data (WES, WGS), all three give ~20–30% saving — none of them works miracles on inherently diverse libraries.

The advantage of `--clumpify` over running a separate tool is that it's part of the trim pass, so there's no extra read-and-rewrite cycle. For an X GB input, a separate clumpify step would mean reading the trimmed output, sorting, and writing it back; typically +3–5× the trim wall time as well as double the disk I/O. `--clumpify` does it in one pass.
