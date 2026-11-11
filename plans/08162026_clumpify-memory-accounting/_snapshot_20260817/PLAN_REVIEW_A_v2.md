# Plan review — `PLAN.md` v2 (honest `--clumpify` memory accounting, #439)

**Reviewer A**, independent, fresh context. Target: `plans/08162026_clumpify-memory-accounting/PLAN.md` (v2).
Repo `/Users/fkrueger/Github/TrimGalore`, branch `dev`, tip `2ab0918`. Read-only; nothing edited, no branch switched.
All line numbers below re-derived against `2ab0918`.

**Verdict: do not implement Phase 0 as written.** v2 is a large improvement on v1 — it stages the
measurement, it drops the unsupported coefficient, and it carries forward the compile and test
findings correctly (with one exception, below). But three things are wrong at a level that would waste
the whole Phase 0 measurement budget:

1. **The calibration fixture is 30,000 distinct reads cycled 227×.** Every compression number and the
   output-buffer contribution to every RSS number is an artifact. Measurement 6 — the sole evidence for
   the Goal's central constraint — does not survive this.
2. **The Goal is internally unsatisfiable on the plan's own fit**, and the plan does not say so.
3. **Phase 0's matrix cannot separate σ from k at `cores ≥ 4`**, so 10 of its 20 core cells buy one
   consistency check between them, while the four core counts that *could* identify the split
   (`cores ∈ {1,2,3,4}`, all of which pin `n_bins = 16`) are half unused.

Detail below, arithmetic shown.

---

## 1. Does the sweep support v2's conclusions?

### 1.1 The stated numbers are arithmetically correct

I refit `results_sweep.tsv` from scratch (OLS on all 7 points per series, `rss_mib ~ a + s·bin_mib_exact`):

```
cores=2  n=7:  rss = 151.14 (±29.55) + 32.709 (±0.809)·B    R²=0.99695
cores=4  n=7:  rss = 301.63 (±13.30) + 38.416 (±0.364)·B    R²=0.99955
```

matching the plan's `151.1 + 32.71·B` and `301.6 + 38.42·B`. Model charge
`(5·n_bins + 7·cores)/4` = `(80+14)/4 = 23.5` at c2 and `(80+28)/4 = 27.0` at c4; ratios
`32.71/23.5 = 1.392` and `38.42/27.0 = 1.423`. **Correct.**

`(1024 − 301.6)/38.42 = 18.80 MiB`; shipped `4·(1024−512)/108 = 18.963 MiB`. **Correct** — shipped is
0.87% *above* the fit's honest maximum.

Breach-grows-with-budget: `c4_B19` (budget 1025) → 1.025; `c4_B56.9` (budget 2048) → 1.215/1.211.
**Correct.**

`--compression` irrelevance: z1 reps 1033.1/1079.7 vs z6 1050.2 at 1G/c4. **Correct as a conclusion**,
though note the z1 rep-to-rep spread (4.5%) is larger than the z1-vs-z6 difference, so the data cannot
resolve a compression effect below ~5%. The plan's phrasing ("does not move RSS") is defensible.

The claimed ~1.37× for `B = 28 MiB` at a 1 GiB budget: I could not find that cell in
`results_sweep.tsv`; the nearest is `c4_B27.4` at *budget 1252*, ratio 1.102. If the 1.37× figure comes
from an unlogged run, it is not reproducible from the data I was pointed at — **cite it or drop it.**

### 1.2 The precision claimed for the central constraint is not there

Propagating the fit uncertainty (t₀.₀₂₅,₅ = 2.571):

| | point | 95% CI |
|---|---|---|
| slope c4 | 38.416 | [37.48, 39.35] |
| intercept c4 | 301.63 | [267.5, 335.8] |
| honest `B` at 1G/c4 | 18.80 | **[17.49, 20.18]** |

So `B_max = 18.8 ± 1.3 MiB`, i.e. ±7%. Measurement 3's "today's `B` is already within 1% of the honest
maximum" is a point-estimate coincidence, not a measurement. The *conclusion* (don't shrink `B` on the
strength of this fit) is right; the *precision* is false and it is load-bearing in the Goal, so fix
the wording. Any Phase 0 fit must report CIs or the same overclaim recurs.

### 1.3 The `--cores 2` series is too noisy to carry `k`

c2 residuals: `-21.3 -10.2 +41.3 +4.2 -15.8 -47.4 +49.2`. The two reps at `B=56.9` differ by
**96.6 MiB (4.9%)**. c4 reps differ by 6.7 and 7.3 MiB (0.3–0.9%). c2 is ~7× noisier.

Since `k = (s₄ − s₂)/2` and `σ = (2s₂ − s₄)/16`, c2's error dominates both:

```
se(k) = √(0.809² + 0.364²)/2 = 0.444   →  k = 2.855 ± 0.444,  95% CI ≈ [1.96, 3.75]
se(σ) = √(4·0.654 + 0.133)/16 = 0.104  →  σ = 1.687 ± 0.104,  95% CI ≈ [1.48, 1.90]
```

Read that carefully: **σ is convincingly under-charged** (CI [1.48, 1.90] vs model 1.25) — that finding
is solid. **k is not** — its 95% lower bound (1.96) barely clears the model's 1.75. The plan's headline
"under-charged ~1.4× consistently, stably across core counts" conflates the two by quoting the *total*
slope, which is 82% σ-driven at c2 and 70% at c4. The honest one-line summary is: *the bin-pool
coefficient is ~35% low; the per-core coefficient is not yet resolved.*

### 1.4 The linear form is not established, and the intercept is an extrapolation artifact

R² = 0.997/0.9996 over `B ∈ [12.8, 56.9]` says a line fits *in that window*. It says nothing about
`B → 0`, which is 4× below the smallest measured point and where the intercept lives. The intercept is
where all the physical interpretation is being drawn from, and it is the least-constrained quantity in
the fit.

The plan already knows this is shaky — its Phase 0 deliverable notes that ~75 MiB/core extrapolated
from c2/c4 contradicts a directly measured 0.94× at 1G/c8. **I confirmed that contradiction and it is
worse than stated.** Using the plan's own decomposition (`σ=1.687`, `k=2.855`, intercept `75.4·cores`)
against `results_fixedbin.tsv`'s `fixB_c8_mi` (budget 1535 MiB, n_bins=32, `B = 4·1023/216 = 18.94`):

```
predicted = 75.4·8 + (32·1.687 + 8·2.855)·18.94 = 603 + 76.8·18.94 = 2057 MiB
measured  = 1535 × (1.15 / 1.0486)              = 1682 MiB          [1.0486 = the 1e6-vs-MiB fix]
```

The extrapolated model over-predicts by **375 MiB (22%)** at c8. A constant ~150 MiB intercept instead
of `75.4·cores` gives 1604 MiB — much closer. So the c2/c4 fit does not extrapolate, and the plan's
own tabled "the intercept is core-dependent" (measurement 2) is the one finding of the six that the
wider data actively contradicts. Demote it to a hypothesis.

`results_reps.tsv` sharpens this: under mimalloc at a 1G budget, `div_c12_mi` = 0.73 and
`div_c16_mi` = 0.69 (corrected: 0.70 / 0.66). Ratios *fall* as cores rise, because `n_bins = 4·cores`
shrinks `B` faster than the real cost grows. A model with a genuine `+75 MiB/core` term would push
those the other way. Whatever the per-core cost is, it is much smaller than the c2→c4 intercept step
suggests.

**A mechanism worth naming, because it predicts exactly this shape:** the model assumes all `n_bins`
bins sit at full `bin_byte_budget` simultaneously. Bins fill at near-identical rates (FNV disperses
evenly) but flush at decorrelated times, so mean total occupancy is somewhere between `n_bins·B/2` and
`n_bins·B`, and the decorrelation grows with `n_bins`. That makes the *effective* σ decrease as
`n_bins` rises — which is a `cores`-dependent σ, not a `cores`-independent additive term. Phase 0
should test that form explicitly; it is a better candidate than a fixed per-core cost and it is
cheaper to make an upper bound out of (σ at its worst, i.e. small `n_bins`).

### 1.5 Fatal: the fixture is a 227× cycled library

This is the finding that governs everything else.

```
$ gzip -dc data/div_R1.fastq.gz | wc -l                       → 27,240,000   (6,810,000 records)
$ gzip -dc data/div_R1.fastq.gz | awk 'NR%4==1{print $1}' | sort -u | wc -l  →     30,000
$ gzip -dc data/div_R1.fastq.gz | awk 'NR%4==2' | head -2000000 | sort -u | wc -l → 28,659
```

6,810,000 / 30,000 = **227.0**. The head IDs are `SRR24827378.*`, the tail IDs are `SRR24766921.*`, so
it is two libraries of ~15 K reads each, concatenated and cycled 227 times. Consequences:

**(a) Measurement 6 is an artifact, and it is the only support for the Goal's central constraint.**
"44.5 MB at `B`=12.8 down to 20.2 MB at `B`=56.9" is not bin size buying compression on sequencing
data — it is bin size controlling how many of the 227 copies of each read land in the same sorted gzip
member. Sanity check on the absolute numbers: 34.6 MB of `_val_*.fq.gz` for 6.81 M pairs is **5.1
compressed bytes per pair** against ~370 raw bytes — 72× compression, where the *input's own* gzip
achieves 6.1×. No sequencing library behaves like that. On genuinely distinct reads the copies do not
exist, so the `dB/d(output)` sensitivity that motivates "compression must not regress" largely
evaporates. **The Goal's constraint is calibrated on a fixture engineered to maximise exactly the
quantity being protected.**

**(b) It biases the fitted per-core slope *downward*, which is the same class of error v1 made.**
The per-core term includes the worker's gzip output buffer, and the bytes actually *touched* in that
buffer are the compressed size. Per in-flight batch at `B = 19 MiB` (≈ 53,800 pairs at 370 B/pair):

```
fixture:   53,800 × 5.1  B ≈ 0.27 MB touched   (0.014·B)
real data: 53,800 × ~75  B ≈ 4.0  MB touched   (0.21·B)   [~5× compression, typical BS-seq]
```

That is ~0.19·B per in-flight batch missing from the fit, ≈ **+0.33·B per core** at the ~1.75
batches-in-flight the fit implies. Real-data `k` would be ≈ 3.19 rather than 2.855. Note where that
lands: v1's refuted `k = 3.0` may have been closer to right than v2's data suggests, for a reason
neither version identified. Calibrating on this fixture and shipping the result reproduces #439 in a
new place.

**(c) Validation row 4 ("no compression regression") is not runnable on it** in any meaningful sense,
and the plan's Remaining risk (b) ("only *relative* movement") understates the problem: on a 227×
cycled fixture the relative movement is also wrong, because it is dominated by cycle-period-vs-bin-size
interaction rather than by sequence similarity.

**Action (Critical): do not run Phase 0 on `div_R{1,2}.fastq.gz`.** Build a fixture with distinct
reads. The repo's own `test_files/` has real (small) libraries — BS-seq, RRBS, clock — and a 6.8 M-pair
diverse fixture can be synthesised with distinct sequences plus a realistic duplicate fraction, or
assembled by concatenating *different* SRA accessions rather than repeating two. Whatever is chosen,
record `sort -u` counts for IDs and sequences in the harness output so this class of error is visible
in the results file, not three layers down.

---

## 2. Can Phase 0's matrix identify the model?

### 2.1 The degeneracy is worse than the plan accounts for

`n_bins = max(16, 4·cores)`. Therefore:

| cores | 1 | 2 | 3 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| `n_bins` | 16 | 16 | 16 | 16 | 32 | 64 |

For the candidate form `RSS = c₀ + c₁·cores + (σ·n_bins + k·cores)·B`, at **any** `cores ≥ 4`
`n_bins = 4·cores`, so the B-slope collapses to `(4σ + k)·cores` — σ and k enter only through the sum
`4σ + k` and are **not separately identifiable from `cores ∈ {4, 8, 16}`**. Adding c8 and c16 to the
matrix therefore contributes:

- zero information about the σ/k split;
- two useful checks that the `{1,2,3,4}`-derived split *extrapolates* (which §1.4 already says it does
  not — so these are the most valuable cells in the matrix for falsification, and the least valuable
  for fitting).

All leverage on the σ/k split comes from `cores ∈ {1,2,3,4}`, where `n_bins` is pinned at 16 and
`cores` varies over a 4× range. **The plan uses 2 of those 4 points.** Adding `cores = 1` and
`cores = 3` triples the lever on `k` for four extra cells, and `cores = 1` is the single most
informative cell in the whole design because it nearly isolates σ (`slope = 16σ + k`, with `k`
minimally weighted).

### 2.2 The design is simultaneously insufficient and wasteful

- **Insufficient**: 2 core counts on the σ/k lever; no `cores = 1` (the σ-isolating cell); no
  uncertainty budget; `cores = 2` is the noisiest series and half the lever rests on it.
- **Wasteful**: 5 bin sizes × 4 core counts × 2 reps × 2 platforms = 80 runs. Five `B` levels test
  linearity, which is worth doing — but at one or two core counts, not four. Ten of the twenty core
  cells (c8, c16) are structurally degenerate for fitting.

**Recommended design (Important), ~21 cells per platform instead of ~20, far more informative:**

| block | cells | buys |
|---|---|---|
| `cores ∈ {1,2,3,4}` × `B ∈ {12.8, 19, 27.4, 37.9, 56.9}`, 2 reps at the two extremes | 20 | σ/k split with a 4× lever at fixed `n_bins`; linearity test at every core count; noise estimate |
| `cores ∈ {8, 16}` × `B ∈ {19, 37.9}`, 2 reps | 8 | falsifies extrapolation and the `n_bins`-dependent-σ hypothesis |
| `--fastqc` at `cores ∈ {1, 4}` × `B ∈ {19, 37.9}` | 4 | FastQC term, **under mimalloc** (see §5.4) |
| SE series at `cores ∈ {2, 4}` × `B ∈ {19, 37.9}` | 4 | one bin pool instead of two |

**Much better still (Important, and it removes the degeneracy outright):** add a hidden/env-gated
calibration override that sets `n_bins` directly, decoupling it from `cores`. Assumption 2 forbids
changing the *shipped* mapping, not the ability to vary it for measurement. With that, a 3×3 grid of
`(n_bins, cores)` identifies σ and k cleanly and independently at any core count, and the whole
`cores ≥ 4` degeneracy discussion disappears. This is a few lines of test-only plumbing that converts
a confounded observational sweep into a designed experiment. If Phase 0 does one thing differently,
it should be this.

### 2.3 A least-squares fit cannot deliver the plan's success criterion

The Goal is "an **upper bound** on actual peak RSS". Validation row 1 demands "ratio ≤ 1.0 in **every**
cell". Phase 0's deliverable is "a fitted model". OLS puts residuals on both sides of the line by
construction — the c4 fit already sits 18.7 MiB *below* the measured point at the default cell
(`301.63 + 38.416×19.0 = 1031.5` vs measured 1050.2). A mean fit can never satisfy an every-cell
upper-bound criterion.

Phase 0 must therefore deliver an **envelope**, not a fit: either max-residual padding, an explicit
margin term, or upper-quantile regression — and the margin must exceed the observed run-to-run spread,
which at c2/`B`=56.9 is 96.6 MiB (4.9%). **The plan never says this, and it is the difference between
Phase 1 succeeding and Phase 1 producing another honest-looking model that breaches.** Critical.

---

## 3. Is the central constraint right? No — as written the Goal is unsatisfiable

Stack the Goal's three requirements against the plan's own fit at the default cell (1G, c4):

```
(i)   actual RSS ≤ predicted peak        [the "upper bound" promise]
(ii)  predicted peak ≤ --memory         [the whole contract; also a tautology of resolve_layout]
(iii) resolved B unchanged at 18.96 MiB  [the Goal's constraint]

measured actual RSS at that cell = 1050.2 MiB
budget                           = 1024   MiB
```

(i) ∧ (ii) ⟹ actual ≤ 1024. Measurement says 1050. **Under (iii), the Goal cannot be met.** No choice
of coefficients fixes this — you are asking a formula to bound a number that is already larger than
the bound.

The only escapes are:

- **(A) reduce actual RSS structurally** by ≥ 26 MiB at that cell → §4;
- **(B) shrink `B`** — and once §2.3's margin is included this is not the 0.86% the plan's "within a
  few %" tolerates. Needing `301.6 + m + 38.42·B ≤ 1024` with a margin `m` of 25–50 MiB gives
  `B ≤ 17.5–18.2 MiB`, i.e. **4–8% below 18.96**, which violates Validation row 3;
- **(C) decouple the sizing formula from the printed prediction** — explicitly forbidden by
  Implementation step 3.

The plan picks none of them and does not acknowledge the conflict. **Critical: state the escape
explicitly.** My read is that (A) is the only one consistent with the rest of v2, which makes §4 not a
candidate but a prerequisite — and it should be said in the Goal, not in item 9.

Two further problems with fixing `B` at the default:

- **Over-constrained, on artifact evidence.** §1.5(a): the compression cost of shrinking `B` is
  measured on a 227×-cycled fixture where it is maximally steep. On distinct reads the penalty is far
  smaller, and a 3–5% `B` reduction may cost well under 1% of output size — which would make (B) a
  perfectly reasonable fix that the Goal has pre-emptively banned. **Re-derive the constraint on a real
  fixture before treating it as fixed.**
- **It is a constraint on one cell, presented as a constraint on the design.** Pinning `B` at
  `1G/c4` while recalibrating a multi-term model forces the recalibration into the terms that *don't*
  bind at that cell, which is how you get a model that is right at the default and wrong everywhere
  else. That is #437's failure mode — the plan's own measurement 4 ("#437 fixed the default cell, not
  the mechanism"). Fixing `B` at the default risks committing the same error one level up.

### 3.1 Interaction with FastQC — the plan's step 4 has the sign wrong for `clump_only`

Item 4 wants FastQC's ~425 MB in the fixed term rather than a flat reservation on every run, and item 5
notes the floor moves. Two things the plan misses:

- In `clump_only.rs`, FastQC is invoked with `cores.max(1)` threads (lines 375, 555–556, 848, 1096), so
  FastQC's own footprint **is** core-dependent there. A decision that `--clump_only` passes `cores = 1`
  to the layout (item 7) would then under-reserve for FastQC on exactly the path that has no fallback
  (§3.2). The two items interact and the plan treats them separately.
- The 425 MB figure itself is stock-allocator only — see §5.4.

### 3.2 `--clump_only` has no fallback, and a raised floor hard-fails it (unlisted)

`main.rs::resolve_clump_layout` (`src/main.rs:178`) warns and returns `None` when the budget is below
`clumpify_min_memory_bytes` (lines 184–199), so `--clumpify` degrades gracefully. The `--clump_only`
dispatch at `src/main.rs:587` **does not go through that function at all** — it calls
`clump_only::clump_only_*`, and all four of those call `clump::resolve_layout(memory_budget_bytes,
cores.max(1))?` (`src/clump_only.rs:275, 413, 748, 905`) with a bare `?`. A raised floor therefore
turns `--clump_only` from "works" into "**aborts**", at the default `--memory 1G`.

Blast radius, measured: `tests/integration_clump_only.rs` and `tests/integration_clump_only_ubam.rs`
pass `--cores 2` and **never pass `--memory`** (checked all 24 invocation sites), so they all run at
the 1 GiB default. Today `resolve_layout(1 GiB, 2)` gives `4·512/94 = 21.8 MiB` — comfortable. If the
fixed term rises past ~1000 MiB at `cores = 2` (entirely plausible if FastQC folds in, or if a per-core
term lands), every one of those tests aborts, and so does every real `--clump_only --memory 1G` run in
the field.

**Important:** the plan must (a) list this asymmetry, (b) decide whether `--clump_only` gets the same
warn-and-continue fallback, and (c) add it to Tests. `--clump_only` is documented as a lossless
reorder-only mode; silently aborting on a default-budget invocation is a worse regression than #439.

---

## 4. Should the structural fix be primary? Yes — and item 10 is mis-scoped

### 4.1 Item 9 (buffer recycling) is a prerequisite, not a candidate

§3 shows the Goal is unsatisfiable without reducing actual RSS. Item 9's mechanism is real and I
verified the code shape: `PairedBin::take()` (`src/parallel.rs:913–919`) `mem::take`s all three Vecs,
which are then `send`-ed to a worker (line 819–821) and dropped there; the next `push` into that bin
hits `self.keys.capacity() == 0` and `reserve_exact`s full capacity again
(`src/parallel.rs:897–902`). So every flush is an allocate-in-reader / free-in-worker cycle for
`~1.4·B` of spine + text handles. That is precisely the pattern mimalloc defers reclaiming
(cross-thread frees return to the owning heap lazily), and libmalloc retains outright.

Recycling makes total allocation `O(n_bins + cores)` for the run and — importantly — makes the
per-core term *structural* rather than *statistical*, which is the only way Validation row 1's
"ratio ≤ 1.0 in every cell" becomes provable rather than sampled. **Cost is modest**: a return channel
from workers to the reader carrying drained `Vec`s, plus a `PairedBin::refill(bufs)` that replaces
`reserve_exact`. The awkward part is that the buffers must survive `sort_paired_by_key`, which
currently `drain`s the three Vecs into a fresh `grouped: Vec<(MinimizerKey, FastqRecord, FastqRecord)>`
(`src/clump.rs:296–312`) — so recycling wants that sort reworked to sort in place over an index
permutation, or the scratch `grouped` buffer recycled too.

**Recommendation: promote item 9 to the primary fix and make Phase 0's *first* experiment a
throwaway prototype of it.** If it removes the per-core term, the model shrinks to
`c₀ + σ·n_bins·B + k'·cores·B` with a small, structurally-bounded `k'`, most of §2's identifiability
problem disappears, and the honest `B` at the default may *rise* — which satisfies the Goal with room
to spare instead of fighting it. That is a strictly better outcome than any recalibration, reached with
less measurement.

### 4.2 Item 8 (queue depth) is the same lever and should be re-rated

`mpsc::sync_channel::<PairedWork>(2)` at `src/parallel.rs:119` (and `:941` for SE). Each worker can
hold 2 queued batches plus 1 in hand = **3 × ~1.4·B per core** worst case, against the model's
`1.25·B` input charge. The fit's `k = 2.855` implies ~1.75 batches resident per core on average, i.e.
the queue is a first-order contributor to the very term the plan cannot pin down. With `capacity(0)`
in clumpy mode the worst case drops to `1 × 1.4·B` ≈ the model's existing charge.

The plan's throughput objection is weaker than it looks: wall times in `results_sweep.tsv` are flat in
`B` (c4: 10.06–10.62 s; c2: 17.9–19.8 s) and scale sublinearly in cores (c2 18 s → c4 10.5 s → c8
7.3 s → c16 6.7 s from the earlier TSVs), which says the **reader is the bottleneck at `cores ≥ 4`** —
it does the minimizer hash *and* the whole per-bin sort serially at `src/parallel.rs:819`. When the
reader is the bottleneck the queues are mostly empty, so shrinking their capacity costs approximately
nothing. Item 8's caveat about sharing the channels with plain round-robin mode is correct and the
gating is the right call; but "Phase 0 decides whether it is needed at all" understates it — it is one
of only two levers that can move actual RSS, and it is a one-line change to measure.

### 4.3 Item 10 is mis-classified as "cheap adjacent"

`let cap = reads_r1.len() * 300;` at `src/parallel.rs:349` (and `:1069` for SE) reserves two buffers of
`300 bytes × n_records` for **compressed** output. In clumpy mode the "batch" is an entire bin, so on
the 65 bp fixture (record ≈ 185 B, pair ≈ 370 B) this reserves `2 × 300/370 = 1.62·B` of address space
per in-flight batch against a model charge of `0.5·B`. Whether that reservation becomes *resident*
depends on allocator commit behaviour and I did not measure it — but the coincidence is worth noting:
the excess `1.62 − 0.5 = 1.12·B` per core is within 2% of the fit's per-core excess
`k − 1.75 = 1.105·B`. That may be a coincidence and I flag it as a hypothesis, not a finding.

Either way the classification is wrong. If the reservation *is* resident, item 10 alone accounts for
the entire per-core under-charge and fixing it drops RSS at 1G/c4 by `4 × 1.105 × 18.96 ≈ 84 MiB` —
three times the 26 MiB gap in §3, i.e. it closes the Goal's contradiction single-handed. If it is *not*
resident, item 10 is a genuine but minor over-reservation. **The two cases differ by 84 MiB and one
`git diff` line separates them.** Move item 10 into Phase 0 as experiment #2, immediately after the
recycling prototype: change `300` to a `raw_bytes`-derived cap, re-run 4 cells, read the answer. Do not
leave an 84 MiB question mark sitting in a "cheap adjacent fix" bullet while an 80-run matrix is
calibrated around it.

Also note the plan's stated motivation for item 10 (RRBS over-reservation at 51 bp) has the sensitivity
backwards: *shorter* reads make `300/bytes_per_record` **larger**, so 51 bp RRBS is the worst case and
150 bp is nearly fine — consistent with Assumption 4's spine factors (1.46 at 51 bp, 1.21 at 150 bp)
but the plan does not connect them.

---

## 5. Code-level claims

### 5.1 Verified correct

- `clump_only.rs` **four** `resolve_layout` call sites: `275, 413, 748, 905`. ✔
- All four have `fastqc_requested: bool` already in scope (params at `258, 396, 743, 893`). ✔
- `main.rs::resolve_clump_layout` is one function (`:178`) called from three loops (`:1344, :1580,
  :1859`), and takes `&Cli` so `cli.fastqc_requested()` is available. ✔
- `Cli::fastqc_requested()` at `src/cli.rs:592` is `self.fastqc || self.fastqc_args.is_some()`, and its
  doc comment states the `--fastqc_args`-alone case explicitly. ✔ Using it rather than `cli.fastqc` is
  right.
- `resolve_layout_predicted_peak_fits_budget` **is** a tautology: `predicted = STATIC +
  floor(4·usable/denom)·denom/4 ≤ STATIC + usable = budget` for all inputs, by the floor. ✔
- Banner prints MiB (`/(1024*1024)`) but labels `MB` — `src/main.rs:203–207`. ✔
- Docs `--cores 8 --memory 1G` → docs say "32 bins × 12 MB"; shipped gives
  `4·(1024−512)/216 = 9.481 MiB`. **Already wrong today.** ✔ (`--cores 8 --memory 4G` → 66.4 MiB, docs'
  "66 MB" is right.)
- Docs floor range: `cores=1` → 533 MiB (docs "~535"), `cores=32` → 728 MiB (docs "535–730"). The
  worked `--cores 6` warning example's "552 MiB" reproduces exactly under the code's truncating
  division. ✔ All three move if the fixed term changes.
- Docs LaTeX at `docs/src/content/docs/performance/clumpy.md:93–100` hard-codes `512 MiB` and the
  `(5·n_bins + 7·cores)` denominator. ✔ Line 108's "mostly from a fixed 512 MiB reservation for
  FastQC" also needs rewording if item 4 lands — the plan doesn't name that line.

### 5.2 The Tests section is still not complete — and "six" is a miscount

`src/clump.rs` has **five** affected tests, not six: `resolve_layout_default` (:445),
`resolve_layout_scales_with_cores` (:455), `resolve_layout_predicted_peak_fits_budget` (:465),
`resolve_layout_bails_when_too_small` (:483), `resolve_layout_bails_below_static_overhead` (:494). The
plan's sixth item — "a `"static overhead"` substring assertion" — is `src/clump.rs:490`, which lives
*inside* `resolve_layout_bails_when_too_small` and is therefore already counted. `grep -rn "static
overhead" src/ tests/` returns exactly that one assertion plus the message itself at `:207`. Fix the
count so nobody hunts a sixth test that isn't there.

The plan's catch that `resolve_layout_bails_below_static_overhead` **panics** rather than fails is
correct and worth keeping: it calls `res.unwrap_err()` (:496) with no prior `is_err()` guard, so if the
fixed term drops below ~283 MiB (`resolve_layout(256 MiB, 4)` returning `Ok` needs
`4·usable/108 ≥ 1 MiB`, i.e. `usable ≥ 27 MiB`) it panics on `unwrap_err`.

**Missing from Tests (Important):**

- `src/clump_only.rs:1185–1191` — `small_memory_budget()` returns 1 GiB with a comment that names the
  512 MiB constant. **13 call sites**: `:1208, 1247, 1298, 1310, 1351, 1384, 1410, 1438, 1467, 1498,
  1565, 1632, 1677`. All are `cores`-small, so they survive today, but they break the moment the fixed
  term at that core count exceeds ~1023 MiB — and the comment goes stale regardless.
- `tests/integration_clump_only.rs` (10 invocations) + `tests/integration_clump_only_ubam.rs` (13) —
  23 in total, all
  `--cores 2`, none passing `--memory` → all at the 1 GiB default, all abort if the floor rises
  (§3.2). The plan lists neither file.
- `tests/integration_preflight_tripwire.rs:592` asserts `stderr.contains("clumpify:")`. The plan cites
  this as the *guard that exists* for the new banner test — correct — but it is also itself sensitive
  to the fallback threshold moving (its comment at :590 says as much). It belongs in the affected list,
  not only in the "new tests" list.

### 5.3 Item 2's compile fix is the wrong remedy for the error it names

`ClumpLayout` (`src/clump.rs:136–140`) has exactly two fields, **both already `pub`**:

```rust
#[derive(Debug, Clone, Copy)]
pub struct ClumpLayout { pub n_bins: usize, pub bin_byte_budget: usize }
```

`src/parallel.rs:1611–1616` builds it with an exhaustive struct literal:

```rust
fn small_clump_layout() -> ClumpLayout {
    ClumpLayout { n_bins: 16, bin_byte_budget: 1024 * 1024 }
}
```

The plan says "any new field is `pub(crate)` or arrives via a constructor". Visibility is the fix for
**E0451** (private field in a literal). It is *not* the fix for **E0063** (missing field in a literal),
which fires on any new field regardless of visibility. So `pub(crate)` alone still does not compile —
`small_clump_layout()` must be updated, or the new state must arrive via a constructor and the literal
be replaced by a call to it. State that plainly; as written the plan repeats v1's mistake in a
different key.

Two more, both minor: `ClumpLayout` derives `Copy`, so a non-`Copy` field is a further break; and item
3 ("store what the layout was sized against and remove the `cores` parameter from
`predicted_peak_bytes`") *is itself* the new-field change, so items 2 and 3 are the same edit and
should be merged.

### 5.4 The FastQC 425 MB figure is stock-allocator, not mimalloc

Item 4's "measured cost is ~425 MB" traces to `results_investigate.tsv`: `pe_1g_c4_fastqc` 1936.9 −
`pe_1g_c4` 1512.0 = 424.9. Both rows are **`tg_stock`**. The shipped binary is mimalloc
(`pe_1g_c4_mi` = 1061.3, a 450 MB difference from stock on the *same* cell), and there is **no
`--fastqc` row under mimalloc anywhere in the scratchpad**. Whatever FastQC costs under the allocator
that ships is unmeasured, and quoting 425 MB as established repeats v1's pattern of a plausible number
carrying an unstated allocator assumption. Phase 0 covers it — but item 4 should say "unmeasured under
mimalloc", not "~425 MB". Also fold the divide-by-1e6 correction (1.0486) into any figure quoted from
these older TSVs; the plan applies that correction to ratios but not to the 425 MB difference.

---

## 6. Validation sufficiency, adversarially

**What passes all 9 rows and still ships a broken promise:**

1. **A mean fit passing "ratio ≤ 1.0 in every cell" by luck of the replay.** §2.3. Rows 1 and 2 replay
   *cells*, and a cell is one sample of a distribution whose observed spread reaches 4.9%. With ~29
   cells and a model tuned to their means, roughly half should breach; if they don't, it is because the
   replay reused the same seeds/ordering, not because the bound holds. **Fix: require each cell's
   ratio ≤ 1.0 at its *maximum* over ≥3 reps, and require a stated margin ≥ the largest observed
   within-cell spread.**
2. **Compression validated on the cycled fixture.** §1.5. Row 4 compares `_val_*.fq.gz` bytes on
   "duplicate-rich **and** low-duplication fixtures" — but if "duplicate-rich" is `div_R*`, the row
   passes on an artifact and says nothing about real libraries. **Fix: name the fixtures in the plan,
   and record ID/sequence `sort -u` counts alongside every compression number.**
3. **A slow writer filling the result channel.** `mpsc::sync_channel::<Result<PairedBatchResult>>(cores
   * 2)` at `src/parallel.rs:128`; each `PairedBatchResult` carries `compressed: Vec<u8>` whose
   *capacity* is the `n × 300` reservation from `:349`. The model has **no term at all** for results in
   flight. On a fast disk they drain immediately and never appear; on a slow or contended filesystem up
   to `3 × cores` results accumulate. Every measured cell here writes to fast local storage, so no row
   probes it. **Fix: add a cell with the output on a throttled/slow sink (or a deliberately stalled
   consumer), or bound the result channel's byte occupancy rather than its item count.**
4. **A read length outside the calibration.** Assumption 4 admits the spine factor spans 1.21–1.46
   across 150/65/51 bp and requires "more than one read length" — but no Validation row tests read
   length, and every measured cell to date is the same 65 bp fixture. **Fix: make read length an
   explicit axis in rows 1 and 2, with 51 bp as the worst case.**
5. **`--clump_only` aborting.** §3.2. No row exercises `--clump_only` at all.
6. **Row 3 passing while the promise fails.** "Default bin size unmoved: within a few % of 18.96 MiB"
   is satisfied by construction if the coefficients are chosen to hit it — it constrains the *formula*,
   not memory, exactly like the tautology the plan correctly calls out in
   `resolve_layout_predicted_peak_fits_budget`. Row 3 is only meaningful paired with a row asserting
   *measured* RSS ≤ budget at that same cell, which row 1 supplies only if the default cell is
   explicitly in the matrix. Make that explicit.
7. **Row 8's noise bar is unspecified.** "within measured noise" with "reps sufficient to beat observed
   variance" — but the wall-time spread in `results_sweep.tsv` at c4 is 10.06–10.62 s (5.5%) across
   *different* `B`, so the noise floor is already ~5%. A channel-capacity change costing 4% would pass.
   State the number of reps and the detectable effect size.

**On "no real sequencing library available on either host":** this is stated as a constraint in
Remaining risk (b), but `test_files/` holds real BS-seq, RRBS and clock fixtures — small, but real.
A 6.8 M-pair fixture built by *shuffling together* several distinct real accessions would be far more
defensible than one built by cycling two, and the accessions are named in the docs page's benchmark
list (`docs/src/content/docs/performance/clumpy.md:272–275`). The constraint is looser than the plan
assumes; downloading one library once is cheaper than an 80-run matrix on an artifact.

---

## 7. Where v2 is right

Recorded briefly, having probed each:

- Staging measurement ahead of code, with Phase 0 producing no diff — correct response to v1.
- Refusing to assume anything about the non-FastQC fixed cost (Assumption 5).
- `Cli::fastqc_requested()` over `cli.fastqc`, with the reason (§5.1).
- `resolve_layout_predicted_peak_fits_budget` is a tautology and defends nothing (§5.1).
- `resolve_layout_bails_below_static_overhead` panics rather than fails (§5.2).
- The docs' `32 bins × 12 MB` is already wrong; the banner mislabels MiB as MB (§5.1).
- Every cell must assert the `clumpify:` banner, because a silent fallback reads as a great ratio.
- Ratios in MiB on both sides; the 4.86% harness error named and corrected.
- Gating any channel-capacity change on clumpy mode rather than changing plain round-robin.
- Leaving the v2.2.0 CHANGELOG's falsified "within 1%" line alone as released history.
- The `BTreeMap` question is genuinely closed.

---

## Action items

### Critical

1. **Rebuild the calibration fixture.** `div_R{1,2}.fastq.gz` is 30,000 distinct reads cycled 227×
   (§1.5). Discard measurement 6, and with it the evidentiary basis for the Goal's "compression must
   not regress" constraint. Every compression number and the output-buffer share of every per-core
   number in v2 is an artifact. Record `sort -u` ID/sequence counts in the harness output from now on.
2. **Name the escape from the Goal's contradiction.** At 1G/c4, measured RSS 1050 MiB > budget
   1024 MiB, so "actual ≤ predicted ≤ budget with `B` unchanged" is unsatisfiable (§3). Choose (A)
   structural reduction, (B) a real `B` reduction with row 3 relaxed, or (C) explicit decoupling — and
   say which in the Goal.
3. **Phase 0 must deliver an upper envelope, not a least-squares fit** (§2.3), with a margin term
   exceeding the largest within-cell spread (currently 4.9%). Otherwise Validation row 1's "ratio ≤ 1.0
   in every cell" is unachievable by the artefact Phase 0 produces.
4. **Promote item 9 (buffer recycling) and item 10 (the `n × 300` reservation) to Phase 0's first two
   experiments, ahead of the matrix** (§4.1, §4.3). Item 10 is one line and settles an 84 MiB question;
   item 9 plausibly removes the per-core term and with it most of §2's identifiability problem.
   Calibrating an 80-run matrix around unresolved 26–84 MiB structural terms is the wrong order.
5. **Fix item 2's compile remedy.** Both `ClumpLayout` fields are already `pub`; the break is E0063
   (missing field), which visibility does not fix. `src/parallel.rs:1611–1616` must be updated or a
   constructor introduced (§5.3). Merge items 2 and 3 — they are one edit.

### Important

6. **Redesign the matrix** (§2.2). σ and k are not separately identifiable for `cores ≥ 4`; add
   `cores ∈ {1, 3}` (all of `{1,2,3,4}` pin `n_bins = 16`), thin `cores ∈ {8, 16}` to 2 budgets and
   treat them as falsification rather than fitting. Better: add a test-only override that sets `n_bins`
   independently of `cores` and turn the sweep into a designed experiment.
7. **Correct §1's statistics in the plan text** (§1.2, §1.3): `B_max` is 18.8 **±1.3** MiB, not "within
   1%"; σ is convincingly low but **k is not resolved** (95% CI [1.96, 3.75] against the model's 1.75).
   Stop quoting the total slope ratio as evidence for the per-core term.
8. **Demote measurement 2 to a hypothesis** (§1.4). A `75.4 MiB/core` intercept over-predicts
   `fixB_c8_mi` by 375 MiB (22%); `div_c12_mi`/`div_c16_mi` ratios *fall* with cores. Add the
   `n_bins`-dependent-σ (bin-flush decorrelation) alternative as an explicit candidate form.
9. **Add `--clump_only`'s missing fallback to scope** (§3.2). It bypasses `resolve_clump_layout`
   entirely and propagates `resolve_layout`'s error with `?`, so a raised floor aborts it — plus 23
   integration-test invocations at the default 1 GiB in `tests/integration_clump_only{,_ubam}.rs`.
10. **Complete the Tests section** (§5.2): five tests not six; add `clump_only.rs`'s
    `small_memory_budget()` and its 13 call sites, both `integration_clump_only*` files, and
    `integration_preflight_tripwire.rs:592` to the *affected* list.
11. **Re-rate item 8** (§4.2). Work channels are `sync_channel(2)`, so worst-case per-core input is
    `3 × 1.4·B` against a charged `1.25·B`. Wall-time scaling says the reader is the bottleneck at
    `cores ≥ 4`, so `capacity(0)` in clumpy mode likely costs nothing — measure it in Phase 0, don't
    defer the question.
12. **Restate the FastQC cost as unmeasured under mimalloc** (§5.4). The 425 MB figure is stock-
    allocator, uncorrected for the 1.0486 factor, and there is no mimalloc `--fastqc` row anywhere.
13. **Add the missing validation rows** (§6): max-over-reps rather than single-cell ratios; read length
    as an explicit axis with 51 bp; a slow-output-sink cell for result-channel occupancy; a
    `--clump_only` row; a quantified noise bar and rep count for row 8.

### Optional

14. Add a term (or an explicit bound) for results in flight — `sync_channel(cores*2)` of
    `PairedBatchResult`, each carrying an `n × 300` capacity (§6.3). Bounding by bytes rather than item
    count would make it disappear from the model.
15. Consider sorting in place over an index permutation instead of `sort_paired_by_key`'s
    `grouped: Vec<(MinimizerKey, FastqRecord, FastqRecord)>` (`src/clump.rs:296–312`), which duplicates
    the spine while the drained sources retain capacity. Needed anyway if item 9 lands.
16. Docs: `clumpy.md:108` ("mostly from a fixed 512 MiB reservation for FastQC") needs rewording under
    item 4, alongside the LaTeX block, the worked values, and the floor range the plan already lists.
17. Update the stale `512 MiB` reference in `src/clump_only.rs:1188`'s comment.
18. Item 10's stated motivation has the read-length sensitivity backwards in emphasis — 51 bp is the
    worst case because `300/bytes_per_record` grows as reads shorten; connect it to Assumption 4.

---

*Reviewer A. All arithmetic recomputed from `scratchpad/rss/results_sweep.tsv` and the earlier
`results*.tsv`; all line numbers re-derived against `2ab0918`.*

---

## Addendum (added after the report above, prompted by Reviewer B's fixture warning)

**My measurements are unaffected by B's truncation.** B truncated `scratchpad/t_R{1,2}.fastq.gz`; my
work used `scratchpad/rss/data/div_R{1,2}.fastq.gz`, which are a different inode pair with
`links=1` (never hardlinked) and size/mtime byte-identical to when I read them
(`206400204`/`Aug 16 06:31:59`, `205310831`/`Aug 16 06:32:04`). Nothing in §1 or §5 needs revisiting.
No new `results*.tsv` exists, so `results_sweep.tsv` remains the only fitted dataset and the §1 fit
review stands as written.

**Critical action item 1 has already been implemented, and the replacement fixture is sound.**
A new `scratchpad/rss/data/dist_R{1,2}.fastq.gz` appeared at 17 Aug 13:33, after the report above.
I verified it rather than assume it:

```
records                 10,000,000
unique IDs              10,000,000     (no cycling — the 227× defect is gone)
unique sequences         3,205,737     (32% distinct)
read length             51 bp, uniform
```

The 32% distinct-sequence rate is not a cycling artifact. The multiplicity histogram is smoothly
long-tailed — 2,124,217 singletons (66% of distinct sequences, 21% of reads), monotone decay through
mult=2 (388,462), mult=3 (181,026) … and an extreme tail of single sequences at 35,744 / 47,323 /
**178,424** copies. A replicated library produces a spike at one multiplicity and nothing else; this is
the shape of a real RRBS run (MspI-concentrated fragment starts, PCR duplicates, plus an
adapter-dimer/low-complexity sequence at 1.8% of the library). Input compression is 7.7× — normal —
so §1.5(b)'s concern about the output-buffer term being ~40× under-exercised is **resolved** by this
fixture. **Critical 1 needs no further work; do not rebuild again.**

Three caveats that survive, and one new one:

1. **The read length changed with the fixture (65 bp → 51 bp)**, so no before/after comparison across
   `div_*` and `dist_*` is interpretable — it conflates composition with read length. 51 bp is the
   *worst case* for the spine factor (Assumption 4: 1.46 at 51 bp vs 1.21 at 150 bp), which is the
   right conservative choice for an upper-bound model; but it means the fitted model will over-reserve
   at 150 bp by an unknown margin. **Important 13's read-length axis still stands** — one 150 bp cell
   is enough to bound the over-reservation.
2. **This is the duplicate-rich arm only.** The plan's own Phase 0 asks for "one low-duplication
   fixture cell per core count" because clumpify is documented to *hurt* WGBS, and Validation row 4
   needs both arms. A 32%-distinct RRBS library is not the low-duplication arm. Still missing.
3. **Headers are unusually long relative to the reads** — 75 characters (including a
   `:TTAGTTGT_1:N:0:NNNNNNNN` suffix) against a 51 bp sequence, so the ID is 41% of each record's
   182 raw bytes. `estimated_record_bytes` counts the id (`src/clump.rs:250–253`), so the spine factor
   fitted on this fixture is tuned to a high header-overhead ratio. Legitimate real data, but worth
   recording alongside the fit so the next reviewer does not have to rediscover it.
4. **New, from B — a direct test for §4.3's 84 MiB question, verified and refined.**
   `io::is_gzipped` (`src/io.rs:33`) matches on the filename extension only, while the reader sniffs
   bytes (`format::detect_input_format`); `src/io.rs:24–28` documents the divergence explicitly, with
   `.bgz` as the worked example ("decompressed correctly but produces plain output"). Output encoding
   is therefore selectable independently of input encoding, which lets a cell hold the input constant
   and turn gzip output off — separating the `cores × 0.5·B` output charge from the input and spine
   terms by subtraction, instead of inferring it from a slope. Two refinements after checking the call
   sites (complete list: `clump_only.rs:285, 431, 1145` and `main.rs:424`):

   - **The lever reaches `--clumpify` via `src/main.rs:424`**, not via the `clump_only.rs` sites:
     `let gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0]);` — the first input's filename
     sets output encoding for the whole trimming run.
   - **Do not measure this through `--clump_only`.** It is single-threaded (`--cores` is interface-only
     per the module comment at `clump_only.rs:33`), so it has no worker pool, no per-worker output
     buffers and no work channels. Isolating the output term there returns a clean null on a path that
     structurally lacks the term, and that null would read as evidence.
   - **`--dont_gzip` is the same lever on the same line and strictly cheaper**: it short-circuits
     before `is_gzipped`, needs no fixture aliasing and no new files (0 bytes, versus ~1.48 GB per mate
     for a decompressed input arm). It also changes output naming to `_val_1.fq`, which is irrelevant
     to an RSS measurement. Use it for the output-term cell; a genuinely decompressed input is only
     needed if the *input* arm must also be plain.

   This is a cheaper and more direct test of §4.3 than the `300`-byte source edit, and both should go
   into Phase 0.

Net effect on the report: Critical 1 is **done**; Critical 2–5 and all Important items stand unchanged.

**Provenance note on this addendum.** Two claims reached me this round that the filesystem and the
source contradicted, and both were traceable to truncated tool output — B's from a stale
`find -samefile`, mine from a `head -20`-truncated `grep` that cut `main.rs:424` and almost had me
report that the `is_gzipped` lever does not reach `--clumpify`. For the record, and because Phase 0
will run on shared fixtures: there are **no hardlinks anywhere** in this session's scratchpad
(`links=1` on all four fixtures; `find -samefile dist_R1` returns only itself), and
`scratchpad/{revA2,B2}/` and `scratchpad/t_R{1,2}.fastq.gz` do not exist. `dist_R{1,2}.fastq.gz` is
intact and passes `gzip -t`. I declined a request to `ln -sf` over paths attributed to me, since they
do not exist and creating them would have manufactured the hardlink hazard being warned about.
