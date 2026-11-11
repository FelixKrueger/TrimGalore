# `--clumpify` predicted peak vs actual RSS — investigation

**Date:** 2026-08-15
**Trigger:** while measuring peak RSS for PR #437 (mimalloc global allocator), a
`--memory 1G` run whose banner promised `predicted peak ≈ 1023 MB` peaked at
1471 MB on macOS/arm64.
**Commit under test:** `115d84f` (dev tip) on both hosts.
**Status:** investigation complete, no code changed. Fix options at the end.

## Summary

The banner is honest on Linux/glibc and wrong on macOS/libmalloc:

| platform / allocator | `--cores 4` | `--cores 8` | `--cores 12` | `--cores 16` |
|---|---|---|---|---|
| macOS, system malloc | **1.48×** | **2.2–2.5×** | **2.77×** | **2.53×** |
| macOS, mimalloc | 1.05× | 0.99× | 0.73× | 0.69× |
| Linux, glibc | 0.91× | 0.91× | — | — |
| Linux, mimalloc | 0.96× | 1.00× | — | — |

Ratios are measured peak RSS ÷ the banner's predicted peak, all at
`--memory 1G --compression 6`, paired-end, on the diverse fixture. A user on
macOS asking for a 1 GB budget at `--cores 8` gets a ~2.5 GB process.

## Method

- Two binaries per host from the same commit: stock, and the same tree with
  PR #437's patch applied (mimalloc verified linked via `nm`; stock verified to
  have no `mi_malloc` symbols).
- Peak RSS from `os.wait4()` per-child `ru_maxrss` (bytes on macOS, KiB on
  Linux). Per-child rather than `RUSAGE_CHILDREN`, which would carry a running
  max across runs. Harness validated by separating a known-different pair
  (plain trim 76 MB vs clumpify 297 MB) before trusting any result.
- Two fixtures, byte-identical across hosts (md5-checked):
  - **duplicated** — `BS-seq_10K_R{1,2}` × 680 = 6.8 M pairs, 2.6 GB uncompressed.
  - **diverse** — `BS-seq_10K` + `clock_10K` + `SRR24766921_RRBS` cycled × 227 =
    6.81 M pairs, three libraries so minimizer diversity is 3× the above. On
    Linux a simulated RRBS pair (`synth_bclconvert`, 236 MB, 51 bp) was used as
    the diverse case.
- Linux host: `dockyard-oxy-0` (128 cores, 991 GB, AL2023, glibc).

Fixture diversity matters: the duplicated fixture compresses unusually well, so
worker output buffers are small and Linux ratios read 0.58–0.74×. On the diverse
fixture the same Linux cells read 0.91×. macOS was insensitive to this (1.48× on
both fixtures at `--cores 4`).

## Where the model goes wrong

`clump::resolve_layout` sizes bins so that

```
peak ≈ STATIC + B × [(n_bins + cores) × 1.25 + cores × 0.5]
     = STATIC + B × (5 × n_bins + 7 × cores) / 4
```

and `predicted_peak_bytes` inverts the same expression, so **the banner is
`≈ --memory` by construction**. Two independent errors hide inside that.

### 1. The per-core term counts one in-flight bin per worker; the real bound is three

`parallel.rs:119` gives every worker `mpsc::sync_channel::<PairedWork>(2)`, so a
worker can hold 2 queued bins plus 1 in hand. Results ride
`sync_channel(cores * 2)` (`parallel.rs:128`) and the main thread reassembles
through an unbounded `pending: BTreeMap` (`parallel.rs:227`). The model budgets
`1.25 × B` of input plus `0.5 × B` of output per core.

Measured with bin size pinned at 18 MB and `n_bins` held at 16 (so only `cores`
varies, `--cores 2` @ 958M vs `--cores 4` @ 1G):

| | per-core cost |
|---|---|
| model (`1.75 × B`) | 35 MB |
| measured, macOS mimalloc | 126 MB |
| measured, macOS system malloc | 150 MB |

The corrected arithmetic fits the Linux measurement directly: at `--cores 4`,
16 bins × 19 MB × 1.25 spine = 397 MB resident, plus 4 workers × 3 in-flight ×
19 MB × 1.25 = 298 MB, total 695 MB against 672 MB measured.

Checked and **rejected**: that a paired bin carries both mates against a
one-mate byte budget. `PairedBin::push` counts
`estimated_record_bytes(&r1) + estimated_record_bytes(&r2)`, so the budget is
per-pair already.

### 2. `STATIC_OVERHEAD_BYTES` is a flat 512 MiB that only FastQC spends

Measured on the diverse fixture, paired, `--cores 4`:

| run | peak RSS |
|---|---|
| plain trim, no `--clumpify` | 101 MB (Linux: 57 MB) |
| `--clumpify --memory 1G` | 1512 MB |
| `--clumpify --memory 1G --fastqc` | 1937 MB |

FastQC costs ~425 MB, which is roughly what the 512 MiB reservation was
calibrated against (see the constant's comment). Without `--fastqc` that
reservation is unspent, and it is precisely what keeps Linux under budget
despite error 1. It is not enough to absorb macOS allocator retention on top.

### Why the platforms diverge

Error 1 is a live-bytes error and exists on both platforms. Whether those bytes
become *resident* is the allocator's call: glibc reuses and returns freed bin
buffers promptly, so the excess mostly does not show up in RSS; macOS libmalloc
retains, so every extra live buffer becomes resident and stays. Hence the same
arithmetic error reads as 0.91× on Linux and 2.5× on macOS.

Consequence for PR #437: mimalloc **repairs** the macOS overshoot
(1.48 → 1.05 at `--cores 4`, 2.52 → 0.99 at `--cores 8`) and **consumes the
remaining Linux headroom** (0.91 → 1.00 on diverse data). Neither breaches the
budget, but on Linux the margin that currently absorbs error 1 goes to zero.

## Impact

- macOS users get 1.5–2.8× the memory they asked for; the overshoot grows with
  `--cores` and plateaus around 2.5–2.8× (retention saturates).
- Linux/glibc — where SLURM `--mem` actually bites — is within budget today, at
  0.91× on diverse data.
- Plain (non-clumpify) trimming makes no budget promise and peaks at ~100 MB
  (macOS) / ~57 MB (Linux) at `--cores 4`, so this is a `--clumpify`-only
  accounting defect.

## Fix options (not implemented)

1. **Count the real queue depth.** Replace the per-core coefficient with the
   channel-capacity truth (3 in-flight bins per worker, `cores × 2` results):
   `(5 × n_bins + 7 × cores) / 4` → roughly `1.25 × n_bins + 5 × cores`. Honest
   on both platforms, but it shrinks bins at a given `--memory`, which costs
   clumpify compression ratio unless the default budget rises with it.
2. **Bound the in-flight set instead.** `sync_channel(2)` → `(1)` per worker,
   and cap the reassembly `BTreeMap`. Keeps bins large, trades some throughput,
   and makes the existing formula closer to true.
3. **Reframe the banner.** Print the budget and a measured-typical range rather
   than a single "predicted peak", and say the estimate assumes glibc. Cheapest,
   least honest.
4. **Spend the static reservation deliberately.** Scale
   `STATIC_OVERHEAD_BYTES` by whether `--fastqc` is set, which frees ~425 MB of
   real bin budget on the common path — but removes the very cushion currently
   masking error 1, so it only makes sense together with option 1 or 2.

Recommendation: option 1 or 2, with a decision needed on whether to hold
compression ratio (raise default `--memory`) or hold memory (accept smaller
bins).

## Not established

- No real sequencing library was measured — both diverse fixtures are simulated
  or library-cycled. The oxy real-data fixtures from earlier benchmarking have
  been cleaned off the box.
- Only glibc and macOS libmalloc tested; musl is untested (and not a release
  target).
- The residual per-core gap (126 MB measured vs ~71 MB explained by 3 in-flight
  bins) is attributed to worker-side trimmed-text and compressed-output buffers
  by elimination, not by instrumentation.

## Raw data

`scratchpad/rss/` in session `720d3ced`: `results.tsv`, `results_huge.tsv`,
`results_investigate.tsv`, `results_fixedbin.tsv`, `results_reps.tsv`,
`results_div_diverse_darwin.tsv`, `oxy_results.tsv`, plus the harnesses
`bench.py`, `investigate.py`, `bench_host.py`, `bench_div.py`.
