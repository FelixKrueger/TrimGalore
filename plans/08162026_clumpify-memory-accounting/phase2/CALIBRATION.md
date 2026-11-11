# Phase 2 calibration — `--clumpify` memory model (#439)

**Date:** 2026-08-16
**Binary:** dev tip `1595d59` (Phase 1 merged: clumpy work-channel depth 1, batch-sized output reserve)
**Harness:** `calibrate.py` (this directory)
**Raw data:** `results_cal_darwin.tsv` (macOS/arm64), `results_cal2_linux.tsv` (dockyard-oxy-0, glibc)
**Status:** calibration complete on both platforms; no code changed. Implementation is a separate step.

## Fixture

`synth_bclconvert` RRBS paired-end, 10,000,000 pairs, uniform 51 bp, 76-char headers,
**1,048,545 distinct sequences per 2M reads (~52%)**. Byte-identical on both hosts
(md5 `a922ad10…` / `15f6ef11…`, streamed from oxy and verified).

This replaces the fixture that invalidated the earlier attempt — 30,000 distinct
reads cycled 227×. The harness now prints reads, distinct sequences and mean
length before measuring anything, because the previous failure was invisible
precisely because nothing recorded sample size.

The ~48% duplication here is a property of RRBS (MspI sites concentrate reads),
not of fixture construction.

## Method

- Cells: `cores ∈ {2,3,4}` × `B ∈ {12.8, 19, 27.4, 37.9, 56.9}` MiB as **fitting**
  cells, `cores ∈ {8,16}` as **falsification** cells (predicted, never fitted).
- `n_bins = max(16, 4·cores)` pins `n_bins = 16` for `cores` 2–4, which is the only
  region where the bin-pool coefficient σ and the per-core coefficient `k` are
  separately identifiable. `--clumpify` requires `--cores >= 2`, so 2–4 is the
  whole identifying range.
- Model form: `peak ≈ f0 + B · (σ·n_bins + k·cores)`, mechanism-derived, with **no
  per-core fixed term** — the earlier "core-dependent intercept" was a
  slope↔intercept artefact.
- Each cell = **max over 2 reps** (the promise is one-sided; a mean fit cannot
  bound anything), and `f0` is taken as the **max** implied across fitting cells.
- Ratios computed in MiB on both sides; every cell asserts the `clumpify:` banner
  appeared, since a silent fallback to plain mode would report ~0.1.

### Harness failure found and fixed mid-run

The first Linux matrix reported a constant 1921.8 MiB across 12 of 18 cells and
fitted σ = −0.32. Cause: the fixture scan held 10M sequences in the measuring
process, and **on Linux a child's `ru_maxrss` is inherited from the forking
parent** — every cell reported the harness's own footprint. Confirmed directly
(child `ru_maxrss` = parent `ru_maxrss` = 71.2 MiB on Linux; 2.0 MiB regardless of
parent on macOS, which uses `posix_spawn`).

Fixed: the scan runs entirely in a shell pipeline, and the harness **aborts if its
own footprint exceeds 200 MiB**. macOS was never affected — verified by re-running
one cell under a 14.3 MiB parent (903.4 / 937.2 MiB against the matrix's 936.5,
inside run-to-run spread).

## Fitted coefficients

| platform | σ (per bin) | `k` (per core) | `f0` | largest within-cell spread |
|---|---|---|---|---|
| macOS/arm64 (mimalloc) | **1.41** | **3.99** | 204.7 MiB | 96.2 MiB |
| Linux/x86_64 (glibc + mimalloc) | **1.29** | **3.71** | 173.1 MiB | 36.9 MiB |
| shipped model assumes | 1.25 | 1.75 | 512 MiB (flat) | — |

Per-`cores` slopes, `n_bins` pinned at 16:

```
macOS   c2: 114.6 + 30.454·B    c3: 152.3 + 34.714·B    c4: 176.0 + 38.425·B
Linux   c2: 104.6 + 27.894·B    c3: 108.6 + 32.093·B    c4: 141.5 + 35.319·B
model      c2:       23.500·B      c3:       25.250·B      c4:       27.000·B
```

**The per-core term is the error, and it is platform-independent**: `k` measures
3.71–3.99 against an assumed 1.75, on two allocators and two operating systems.
σ is close to its assumed value. `f0` is 2.5–3× smaller than the flat reservation.

## Falsification (unfitted cells)

| cell | macOS: model vs measured | Linux: model vs measured |
|---|---|---|
| `cores 8`, B=19.0 | 1668 vs 1643 MiB (−1.5%) | 1521 vs 1560 MiB (+2.5%) |
| `cores 8`, B=37.9 | 3125 vs 3053 MiB (−2.3%) | 2863 vs 2792 MiB (−2.5%) |
| `cores 16`, B=19.0 | 3132 vs 2991 MiB (−4.5%) | 2870 vs 2236 MiB (**−22.1%**) |

Five of six predictions land within 5%. The `cores 16` Linux cell over-predicts by
22% — conservative (safe) but loose, and it points at reviewer A's alternative
hypothesis: with 64 bins, not all bins are simultaneously full, so the pool term is
sub-linear in `n_bins`. σ is probably a decreasing function of `n_bins` rather than
a constant. **This costs high-core users compression, not safety.**

## What honest sizing implies

Using the worst coefficients across both platforms (σ = 1.41, `k` = 3.99,
`f0` = 204.7 MiB) with a 10% margin — the margin exceeds the largest within-cell
spread on either host:

| cell | shipped `B` | honest `B` | change | modelled peak/budget (macOS / Linux) |
|---|---|---|---|---|
| 1 GiB / `c2` | 21.79 MiB | 23.78 MiB | **+9.1%** | 0.91 / 0.82 |
| 1 GiB / `c4` | 18.96 MiB | 18.85 MiB | **−0.6%** | 0.91 / 0.82 |
| 1 GiB / `c8` | 9.48 MiB | 9.43 MiB | −0.6% | 0.91 / 0.82 |
| 1 GiB / `c16` | 4.74 MiB | 4.71 MiB | −0.6% | 0.91 / 0.82 |
| 2 GiB / `c4` | 56.89 MiB | 43.02 MiB | **−24.4%** | 0.91 / 0.83 |
| 4 GiB / `c8` | 66.37 MiB | 45.68 MiB | **−31.2%** | 0.91 / 0.83 |

**The correction is free at the default budget and pays for itself at large ones.**
The flat 512 MiB reservation over-charges by ~307 MiB at every budget, which
happens to cancel the under-charged per-core term while the budget is small and
stops cancelling as the dynamic term grows. So:

- at 1 GiB, bin size moves by −0.6% (and *up* 9% at `--cores 2`);
- at 2 GiB and above, bins shrink 24–31% — exactly where the shipped formula
  breaches (measured worst ratios: macOS 1.193, Linux 1.091), and where the README
  already documents sharply diminishing returns;
- v3's constraint "do not move the default bin size" is satisfied **by**
  correctness rather than traded against it.

## Recommended implementation (not done — needs its own step)

1. `f0` becomes ~205 MiB with FastQC handled as a `max`, not an addend (FastQC runs
   after trimming, so the peak is `max(trim, fastqc)`); FastQC's cost under
   mimalloc is still unmeasured.
2. Denominator `4·(σ·n_bins + k·cores)` with σ = 1.41, `k` = 4.0 — i.e. roughly
   `(6·n_bins + 16·cores)` in the code's integer form.
3. An explicit named margin (10%) rather than sizing to the fitted mean.
4. Update the five affected `clump.rs` tests, the derivation docstring, and
   `resolve_layout_pins_the_default_cell` (which will move by −0.6%, deliberately).
5. Docs arithmetic: the formula, floor range, worked values and unit labels in
   `docs/performance/clumpy.md`.

## Limitations

- **One fixture, one read length, one data type.** 51 bp RRBS at ~52% distinct.
  Read length changes σ (the spine factor varies with record size), and clumpify
  is documented to *hurt* coverage-diverse data (WGBS, scRNA R2) — a
  low-duplication cell is still missing, and bin size moves that crossover.
- **No `--fastqc` series**, so the 512→205 MiB reduction is only justified for
  non-FastQC runs. Item 1 above must be measured, not assumed.
- **`cores 16` is over-predicted on Linux by 22%**, so the fixed-σ form is loose at
  high core counts.
- Single-end and `--clump_only` untested here; both share the layout.
- macOS spread is 96 MiB (2.6× Linux's), so the 10% margin is sized by the noisier
  platform.
