# Drafts for review — not posted

Two artifacts: an issue for our `--memory` accounting defect, and a review reply
on PR #437. Post the issue first so the reply can cite its number.

---

## 1. Issue

**Title:** `--clumpify: --memory predicted peak understates actual RSS by up to 2.8× on macOS`

**Label:** `bug`

**Body:**

`--clumpify` prints a predicted peak at startup and sizes its bin pool so that
peak stays inside `--memory`. On macOS the actual peak runs 1.5–2.8× the printed
number; on Linux/glibc it stays inside the budget. Measured at commit `115d84f`,
paired-end, `--memory 1G --compression 6`, on a 6.8 M-pair fixture — ratios are
peak RSS ÷ the banner's own prediction:

| | `--cores 4` | `--cores 8` | `--cores 12` | `--cores 16` |
|---|---|---|---|---|
| macOS/arm64, system malloc | 1.48× | 2.2–2.5× | 2.77× | 2.53× |
| Linux/x86_64, glibc | 0.91× | 0.91× | — | — |

So `--memory 1G --cores 8` on macOS produces a ~2.5 GB process against a
promised ≈1023 MB.

Two causes, both in `clump.rs`:

The per-core term budgets one in-flight bin per worker (`1.25 × B` input plus
`0.5 × B` output). The real bound is three: `parallel.rs:119` gives each worker
`sync_channel::<PairedWork>(2)`, so two bins can be queued while a third is in
hand, and results then ride `sync_channel(cores * 2)` into an unbounded
reassembly `BTreeMap`. Holding bin size at 18 MB and `n_bins` at 16 so only
`cores` varies, the measured per-core cost is 126–150 MB against the model's
35 MB. The corrected arithmetic predicts the Linux measurement closely:
16 bins × 19 MB × 1.25 spine + 4 workers × 3 in-flight × 19 MB × 1.25 = 695 MB,
against 672 MB measured.

`STATIC_OVERHEAD_BYTES` is a flat 512 MiB, and only FastQC spends it — plain
trimming peaks at 101 MB, and `--fastqc` adds ~425 MB. On the common no-FastQC
path that reservation is unspent, and it is what keeps Linux inside budget
despite the first error. glibc also recycles freed bin buffers promptly, so the
excess live set never becomes resident; macOS libmalloc retains it, which is why
one arithmetic error reads as 0.91× on one platform and 2.5× on the other.

Repro:

```
trim_galore --paired --cores 8 --clumpify --compression 6 --memory 1G R1.fq.gz R2.fq.gz
# banner: clumpify: 32 bins × 9 MB; predicted peak ≈ 1023 MB
# actual:  2.2-2.5 GB peak RSS (macOS/arm64), 0.93 GB (Linux/glibc)
```

Peak RSS via per-child `getrusage`; input needs to exceed the budget or the bins
never fill and the overshoot hides.

Fix directions, with a trade-off that needs deciding:

1. Count the real queue depth — `(5 × n_bins + 7 × cores) / 4` becomes roughly
   `1.25 × n_bins + 5 × cores`. Honest on both platforms, but it shrinks bins at
   a given budget, which costs compression ratio unless the default `--memory`
   rises to compensate.
2. Bound the in-flight set instead — `sync_channel(2)` → `(1)` per worker and a
   cap on the reassembly map. Keeps bins large, trades some throughput.
3. Scale `STATIC_OVERHEAD_BYTES` by whether `--fastqc` is set, which frees
   ~425 MB of real bin budget — but only alongside (1) or (2), since that
   reservation is currently the cushion.

Surfaced while measuring #437; unrelated to that change, though mimalloc happens
to mask it on macOS (1.48× → 1.05×) and to spend the Linux cushion
(0.91× → 1.00×).

---

## 2. Review reply on #437

Thanks for retargeting — CI is readable now, and lint and audit are green as
expected.

On the measurements: I re-ran your comparison at 6.8 M read pairs, seventeen
times your fixture, on macOS/arm64 and on a 128-core Linux host with identical
input bytes. The gain survives scale on the `--clumpify` path — about 18% wall at
`--cores 4` on both platforms — but plain trimming moves only 1.5% on Linux
against 8% on macOS, which is the writer-bound picture our profiling showed. It
also depends heavily on how much of the run is gzip work: on a more
compression-heavy dataset the Linux clumpify gain fell to 4–6%. So the effect is
real but narrower than a flat 25%, and the commit message would be stronger for
saying which path and which platform the headline number comes from. Peak RSS
drops 27% at `--cores 4` on the memory-budgeted path on macOS and considerably
more at higher core counts, while on Linux it rises but stays inside budget —
worth recording, and it surfaced #NNN, which is ours rather than yours. Your musl
and bioconda question I can close: the recipe already declares
`{{ compiler('c') }}`, and no release target or container is musl, so the C build
costs us nothing on shipped channels.

Three small things and I am happy to merge. The manifest comment should be a
single line with the evidence left in the commit message; `default-features =
false` is a no-op that needs the claim dropped, since mimalloc's default feature
set is empty and `secure` is opt-in, which unfortunately means an amend now that
it is in the commit message; and I would rather pin `=0.1.52` to match how
`fastqc-rust` and `noodles` are handled, unless you think floating is right
here. The "fails on clean master" note in the description is also moot now.
