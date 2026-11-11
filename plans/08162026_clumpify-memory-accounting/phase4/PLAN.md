# Plan — Phase 4: charge the per-record spine honestly (#457), v4

**Parent:** `../PLAN.md` (#439, closed) · **Phase 3 close-out:** `../phase3-close/PLAN.md` v3
**Issue:** [#457](https://github.com/FelixKrueger/TrimGalore/issues/457) — the last open issue of the arc
**Baseline:** `dev` at `c7dd106`
**Reviews:** `PLAN_REVIEW_reviewer-A_phase4.md`, `PLAN_REVIEW_reviewer-B_phase4.md` — both "do not treat §6 as decided"
**Status:** **decided — Option C, at `σ₁₆ = 33`, `k₁₆ = 92`.** The matrix ran (§10) and two spikes (§11,
§12) removed C's competitors. Dual review of v3 (`PLAN_REVIEW_v3_A.md`, `PLAN_REVIEW_v3_B.md`) reproduced
§6's arithmetic exactly and then **changed the integers**: v3's margin mixed two cells, and under the
correct statistic `(32, 89)` does not clear. Steps 1–3 are complete; Step 4c has not been started.

## Revision history

- **v1 — Option C chosen, on wrong arithmetic.** Its diagnosis and structural form were right; four of
  its numbers were not. Both reviewers, independently and without shared state, reached the same
  Critical: **the per-record cost is not an unknown to be fitted — it is `size_of`, and it is already
  written into `parallel.rs`.** v1's two-point fit gave `s = 28.87 B` against a structural **74 B**,
  2.6× low, and every cost figure downstream was ~3× too favourable.
- **v2 — Option C vs D reopened on the code-derived constant.** Replaced the fitted constants with the
  `size_of` decomposition, rebuilt both arithmetic tables with the static terms handled correctly,
  re-scoped the matrix from *discovering* `s` to *confirming* it, fixed a falsification cell that could
  not falsify, added Option D, and ran the free pre-screen.
- **v3 — Option C, and §6 re-derived from measurement.** The matrix falsified v2's
  central premise: `s = 108.74` fitted against `size_of`'s 74 (§10), and the spike decomposition put the
  true per-slot cost at 245.0–269.6 B/pair of which the spine is only 148 (§11). A second spike showed the
  remainder is not attenuated by occupancy but carries a read-length-dependent fragmentation term (§12).
  **Both findings kill a constant-per-record scheme and neither touches C's logic**, so §6 now derives
  C's denominator from the measured breach — a quantity that turns out to be insensitive to the one
  parameter the design cannot identify. §6's option comparison, `r_min`-by-fit, and Option D are
  withdrawn; §10–§12 are retained as the evidence trail.
- **v4 — this document. Same option, different integers, and a re-scoped Step 4c.** Two independent
  reviewers reproduced every figure in §6 and certified `f_required`; both then found the work downstream
  of the decision badly under-specified. Four changes of substance: **(1)** the spread margin mixed one
  cell's absolute spread with another's `f`; under the correct relative statistic `(32, 89)` fails and
  **`(33, 92)`** is adopted (§6b). **(2)** §6e's claim that the CI band is unaffected was **inverted** —
  lever A's magnitude is proportional to `B`, so at width 0.09 the margin against the only calibrated
  regression goes *negative*; the width narrows to 0.065 (§6e). **(3)** `calibrate.py` regressed on the
  banner's truncated MiB because v3's harness dropped phase 2's `b_exact_mib` column, so §10's `a`/`b`
  and §12's ratios need recomputing — §6's integers do not (§8 Step 4c.1). **(4)** Step 4c grew from "one
  line of arithmetic" to five failing tests, ten hard-coded expectations, a test that stays green while
  its invariant becomes false, three stale docstrings, and two CI gates that this refit exists to arm.

## 1. Goal

`--clumpify` charges the bin pool in **bytes** while part of its real cost is **per record**. At a fixed
byte budget a shorter read packs more records into the same bytes, so the charge does not move and the
reservation under-covers. Measured consequence at 25 bp / `--cores 2`: the printed prediction is
exceeded on both platforms, and the `--memory` budget on macOS/arm64.

Phase 4 charges the spine where it is incurred and re-establishes the promise across the supported range.

## 2. Context

### 2.1 The defect, in the code

`src/clump.rs:174`, used as `dyn_bytes = bin_byte_budget × dyn_denominator / 16` (`:167`):

```rust
fn dyn_denominator(n_bins: usize, cores: usize) -> u64 {
    23 * n_bins as u64 + 64 * cores as u64          // σ = 23/16 = 1.4375,  k = 64/16 = 4
}
```

`predicted_peak = static + dyn`, `static = 224 MiB` charged (measured 205 macOS / 173 Linux,
`:249-251`), `B = 16 × (memory × 10/11 − static) / dyn_denominator`.

### 2.2 Why Phase 2's fit could not see it — verified

`phase2/calibrate.py` varies `cores` and `budget`. It computes mean read length at `:65-73` and only
**prints** it at `:135`; every cell used one fixture. The fit had no length variation to attribute cost
to, so σ absorbed a length-dependent quantity as a constant. The defect follows from the calibration
design.

### 2.3 The spine is `size_of` — but the spine is not the whole per-record cost  [REVISED BY SPIKE]

**The correction that reshapes this plan**, reached independently by both reviewers.
`src/parallel.rs:901-911` (`PairedBin::push`; `SingleBin::push` at `:1286` is identical):

```rust
let bytes_this = estimated_record_bytes(&r1) + estimated_record_bytes(&r2);
if self.keys.capacity() == 0 {
    let predicted = self.budget.div_ceil(bytes_this).max(1);
    self.r1.reserve_exact(predicted);
    self.r2.reserve_exact(predicted);
    self.keys.reserve_exact(predicted);
}
```

`take()` (`:920`) uses `mem::take`, so capacity returns to zero and this **re-reserves once per bin per
flush cycle**. Three consequences:

- The spine is `reserve_exact`, so it is *exactly* `B/bytes_this` slots — not a doubling cascade, not an
  estimate. `size_of::<FastqRecord>()` = **72 B** (three `String`s; measured, not assumed),
  `MinimizerKey` = `u32` = **4 B**, so a paired slot is `72+72+4 = 148 B` and `bytes_this = 2r`. The
  spine is therefore **`74 × B/r`** bytes per bin: `s = 74` B per record-byte-of-`r`, **no fitting**.
- The spine is resident **in full from the first push after each flush**, while the bin's text
  oscillates over `0..B`. That is why the occupancy factor is below 1, not above it (§2.4).
- The same cost rides the **worker** term: an in-flight batch is a whole bin, text *and* spine. `k = 4`
  is documented as "queued batch + batch in hand + 2 result slots"; the first two carry records, the
  result slots carry compressed bytes. So `1/r` multiplies `n_bins` **and** `workers` — which v1's form
  did not.

**The decomposition, solved against the shipped constants** at the fitting fixture's `r` (147.9 B, §2.4a):

```
σ(r) = 0.937 + 74/r        →  σ×16 = 23.00 at r = 147.9   (shipped: 23)
k(r) = 2.999 + 148/r       →  k×16 = 64.00 at r = 147.9   (shipped: 64)
```

Two constants, one equation each, both landing on the shipped integers to two decimals. The `2.999`
resolves as 2 in-flight copies plus 1.00 compressed slot, matching `k`'s own docstring.

**Zero-parameter confirmation.** Both `1G_c4` cells share the same `B` (which does not depend on `r`),
so their difference isolates the per-record term:

```
74 × (n_bins + 2·workers) × B × (1/r₂₅ − 1/r₅₁)
  = 74 × 24 × 19,006,356 × (1/95.9 − 1/147.9) = 118.0 MiB
measured: 117.3 MiB Linux (0.6% off, 27 obs/cell) · 121.9 MiB macOS (3.2% off)
```

**Action:** pin `size_of::<FastqRecord>()` beside `estimated_record_bytes`, so the model's constant and
the struct cannot drift apart silently.

### 2.4 v1's two-point fit is withdrawn  [REVISED BY SPIKE]

v1 fitted `σ(r) = a + s/r` through two points, got `s ≈ 28.7`, `a ≈ 1.238`, and read `a` as "the pool
holds ~24% more than one copy of each bin's bytes". All three are withdrawn:

- **`s` was 2.6× low** because the fit mixed a *safety-rounded model constant* with a *raw measurement*.
  `1.4375` is not "σ measured at 51 bp" — it is `23/16`, the rounded-up envelope of Phase 2's measured
  1.29 (Linux) / 1.41 (macOS) (`clump.rs:281-287`). Anchoring one end of a slope on a deliberately
  conservative value biases the slope toward zero.
- **`a > 1` was not available.** A 72-byte spine per record is a hard floor on the per-record term, so
  any fit returning `s < 72` is absorbing something else. The mechanism points the other way: bins
  flush at `raw_bytes ≥ budget`, so text occupancy is ~uniform on `(0, B)` while the spine never is —
  hence `a = 0.937`, ~94% of one copy. Validation 2 should expect **`a < 1`, excluded from below**;
  v1's framing would have made a correct fit look like a failed one.
- **`1.55` at 25 bp cannot be reconciled** with the arc's own published ratios. Back-solving σ from
  macOS `25bp_1G_c4 = 1.07` of prediction gives **≈1.73**; from the c2 breach, **≈1.86–1.91**. Both land
  on the structural answer, not on 1.55.

### 2.4a `r` must be measured at the fixture size actually used

v1's `r` table came from 2000-pair fixtures. `examples/mk_clump_fixture` writes the record index into
the header (`@SIM:1:FCX:1:15:6329:{i} …`), so mean header length grows with record count — ~3.4 B
between 2 000 and 10 000 000 pairs. Corrected to the 10 M scale the job and the calibration use:

| read length | `r` (B/record) | v1 said |
|---:|---:|---:|
| 25 bp | **95.9** | 92.4 |
| 36 bp | **117.9** | 114.4 |
| 51 bp | **147.9** | 144.4 |
| 75 bp | **195.9** | 192.4 |
| 100 bp | **245.9** | 242.4 |
| 150 bp | **345.9** | 342.4 |

Header length is a real axis, not an artefact: it is per-record overhead that *raises* `r` and therefore
*lowers* σ, so real data with long headers is always safer than a synthetic fixture at the same length.

### 2.5 The lever's direction

`prediction = budget × 10/11` **by construction** — it does not move when σ moves. Raising the charge
enlarges the denominator, shrinking `B`, shrinking the pool and therefore the **actual** peak. The ratio
is corrected from above, and paid for in bin size, i.e. in compression.

### 2.6 Prerequisite: the calibration harness cannot run

`phase2/calibrate.py:35` expects `predicted peak ≈ (\d+) MiB`; the banner is now
`… ≈ 930 of 1024 MiB budget (gzip level 1)` after #451 and #468. It matches nothing and the harness
exits on every cell. Repair it so a future banner change fails **loudly**: named groups, and assert the
parse produced values.

## 3. Behaviour the refit must deliver

1. The spine is charged where it is incurred, so the promise holds across the supported range without
   depending on a constant that read length can invalidate.
2. Integer arithmetic and `saturating_mul` preserved, and `predicted_peak_bytes()` still the exact
   inverse of `resolve_layout` (currently within ~312 B) — the CI job's single-source parse and the
   banner's inability to disagree with the sizing both rest on it.
3. Output bytes change again: bin sizes move, so gzip member boundaries move. This is the **third** such
   change in one cycle (#442, #446, now), and the CHANGELOG should present it cumulatively.
4. `--clump_only` shares `dyn_denominator` and moves with it.

## 4. The calibration matrix — re-scoped

`s` is no longer being discovered, so the matrix's job changes from *fitting three parameters* to
**confirming a code-derived one and fitting the two occupancy terms** (`0.937`, `2.999`). That shrinks
it — and two design faults in v1's version must be fixed first.

**Fault 1 — v1 held pair count constant across the read-length axis.** Total bytes then vary with the
fitted variable, so short reads get far fewer flush cycles and sit further from steady state, which
**depresses their peak and biases the per-record term downward — the unsafe direction, on the axis being
fitted**:

| read length | `r` | bytes/mate at 10 M | flush cycles at `B` = 18.1 MiB |
|---:|---:|---:|---:|
| 25 bp | 95.9 | 0.96 GB | ~106 |
| 51 bp | 147.9 | 1.48 GB | ~163 |
| 150 bp | 345.9 | 3.46 GB | ~382 |

`SESSION_HANDOFF.md` §5 records the ramp this rides: 0.38 at 500 K pairs → ~0.95 at 10 M. **Equalise
total bytes, not pair count** (`pairs ∝ 1/r`, ~31 M pairs at 25 bp for parity with 8.7 M at 150 bp), and
state the disk cost — the arc's first trap is that disk exhaustion inflates peak RSS.

**Fault 2 — v1's falsification cell could not falsify.** Doubling a ~41-char header at 51 bp moves `r`
to ~189, giving Δσ ≈ 0.11 → **~32 MiB** of peak against within-cell spreads of 96 MiB (macOS) and 37 MiB
(Linux): 3× below the noise on the platform measured first. Either run it where it is detectable —
25 bp, `--cores 2`, 4–8× header contrast, ≥5 interleaved reps, Linux for its tighter spread, giving
Δσ ≈ 0.23 → ~84 MiB — or replace it with the exact test, which is asserting `size_of::<FastqRecord>()`
in code.

**Cells.** Read lengths {25, 36, 51, 100, 150} × `--cores` {2, 4} × 1G, byte-equalised, 3 reps, both
platforms. Falsification: `--cores` {8, 16} × 2G × {25, 100}, predicted from the fit and not fitted to.

**Discipline is not optional** — assert free disk before measuring, check every invocation's exit status,
hold records byte-identical across compared conditions, interleave reps, and check output *filenames*
rather than the banner.

## 5. Validation

Steps 1–3 are done; rows 1–4 record what they returned. Rows 5–16 are what the implementing PR must prove.
Rows marked **gate** are the ones whose absence would let the change ship unverified.

| # | Verify | How | Expected | State |
|---|---|---|---|---|
| 1 | The harness reports its own failure | broken pattern against real banner text | aborts | **done** (§10) |
| 2 | `s` matches `size_of` | 10-cell fit | — | **done — falsified**, 108.74 vs 74 (§10). Coefficients pending the Step 4c.1 recompute |
| 3 | Per-slot cost, measured directly | slope in record count | — | **done.** spine exactly 148 B/pair; total **245.0–269.6** B/pair (§11) |
| 4 | Whether the overhead rides occupancy | cycling pool, slope in `B` | — | **done — neither.** Mechanism attributed; **magnitude not transferable** (§12) |
| 5 | **gate — the promise holds at the calibrated floor** | re-run §4's 25 bp cells after the change, **≥10 reps on the binding cell** (25 bp c8 2G, 57 s each) | every rep `≤` printed prediction; predicted headroom 6.1% | to do |
| 6 | **gate — the CI promise gate is armed** | `25bp_1G_c2` and `25bp_1G_c4` move from `gating 0` to `1` in `ci.yml`, and the `#457` "until that refit lands" comment is deleted from the baseline file | both cells gate on `peak ≤ prediction` | to do |
| 7 | **gate — the band still detects lever A** | re-run the work-channel-depth revert once, as the 0.126–0.151 figures were obtained | measured delta ≥ width + a working margin at the new width 0.065 (§6e) | to do |
| 8 | Linearity in `B` held | compare re-measured ratios against §6c's predicted column | within the post-change **measured** within-cell spread from the same runs — **not** the fit's 19 MiB RMS, which is below the measurement noise | to do |
| 9 | Linux agrees, or is worse | §4's cells on Linux/x86_64, after Step 4c.1 | if Linux needs a larger `f`, take Linux's. Predicted band centres 0.607 / 0.678 double as a free check | to do — **blocks merge** |
| 10 | The budget axis is separated from cores | add **25 bp c8 at 1G**, and a second budget (512M) on one cell | `f` at c8 is a cores effect, not a budget effect; the second budget also identifies `static` | to do |
| 11 | `predicted_peak_bytes()` is still the exact inverse | the existing one-sided `predicted <= budget` test — **no test asserts a 312 B bound; `312` appears nowhere in `src/` or `tests/`** | passes; the residual bound rises to ~434 B at cores 32 | to do |
| 12 | The minimum viable budget moved, and the failure is a **silent capability loss** | `--clumpify` **warns and falls back to plain mode** (`main.rs:186-204`), it does not error. `--clump_only` **does** error | max cores keeping clumping engaged at 1G: **72 → 50**; with `--fastqc` **33 → 23**; `--clump_only` at 512M **41 → 28** (v4 first wrote 52 and 29, copied from a review rather than recomputed; the binary gives 50 and 28). Floors: c2 280.5 → 295.4 MiB, c16+fastqc 840.4 → 915.2, c32+fastqc 1012.0 → 1161.6 | to do |
| 13 | `--clump_only` moved with it | shares `dyn_denominator` (§3.4); it has no cell of its own and inherits `f` | banner and layout change consistently | to do |
| 14 | Output size — **the cost this change actually imposes** | two runs on a real RRBS/WGBS library from `docs/perf_data/`, `--memory 1G` vs `1363M`. **The Phase 4 fixture is 100.0% distinct and cannot measure a clumping gain** | quantify the savings lost; if material, say so in the docs | to do |
| 15 | The banner string survives | three integration tests assert `"930 of 1024 MiB budget"` verbatim | unchanged — `prediction` is `budget × 10/11` and does not move | to do |
| 16 | Perl byte-identity untouched | the `validation` job | green — assert, don't assume. Output *bytes* move (gzip member boundaries); no validated flag path runs `--clumpify` | to do |

## 6. The decision — resolved: **Option C** at `σ₁₆ = 33`, `k₁₆ = 92`

### 6a. Why C, and why the alternatives are gone

Each measurement removed an option, and none of them touched C:

| option | what it needed | what was measured |
|---|---|---|
| **D** — charge the spine in the flush threshold | a `size_of`-derived per-slot constant, exact at every read length | the spine *is* exactly 148 B/pair (§11), but that is ~55–60% of the **245.0–269.6** B/pair real cost, and the remainder varies ~10% across read lengths |
| **allocation-count variant** (proposed §11, withdrawn §12) | a per-record constant charged at the reserved slot count | a fragmentation term sits on top of that count and itself scales with `1/r` (§12) |
| **C** — a conservative floor | nothing about mechanism; only a bound that holds at the shortest supported input | unaffected. Every finding moves its numbers; none moves its logic |

The pattern is the same one twice: **every scheme that tries to be exact needs a constant that is right at
more than one read length, and no such constant exists.** C never claimed one.

### 6b. C derived from measured breach, and the margin rule stated

v3's §6b picked `r_min` and computed `σ(r_min)` from `σ(r) = a + s/r`. That route is closed: the model
omits fragmentation (§12) and its occupancy term is unidentifiable at one budget (§10 Finding 1). Derive
the denominator from what was measured instead.

`prediction = budget × 10/11` **by construction**, so it does not move. Listing which terms are constant
under a denominator scale `f` *before* computing the comparison — the guard v2's Self-Review said it
failed to apply:

| term | under a denominator scale `f` |
|---|---|
| `static` | **constant** — not part of the bin pool |
| `prediction` | **constant** — `budget × 10/11` by definition |
| `B` | scales `1/f` |
| `peak − static` | scales `1/f` |

```
f_required = (peak − static) / (prediction − static)
```

**The result is insensitive to `static`,** which is the point — `static` is the one quantity this design
cannot identify, and C's constant does not depend on it. `f = 1.435` tolerates `static` up to ≈434 MiB,
and the binding cell is stable across 173–450 MiB:

| `static` assumed | 173 | 205 | 224 | 260 |
|---|---:|---:|---:|---:|
| binding `f` | 1.3306 | 1.3370 | 1.3409 | 1.3486 |

Binding cell throughout: **25 bp, `--cores 8`, 2G, rep 2** (2419.1 MiB against a 1861 MiB prediction).
Per-cell `f_required` is in §6c's table.

**The margin rule — stated, because v3 applied one without stating it and got it wrong.** v3 took the
largest *absolute* within-cell spread (64.1 MiB, from the c16 cell), expressed it as a fraction of *that*
cell's dynamic part (3.1%), and applied it to the *c8* cell's `f`. Two cells, one number. The correct
statistic is the largest *relative* within-cell spread, and the rule is:

> `f_target = f_binding × (1 + worst relative within-cell spread observed anywhere)`

| construction | spread | `f_target` | does `(32, 89)` = 1.391 clear? |
|---|---:|---:|---|
| (a) the binding cell's **own** spread | 0.65% | 1.3457 | yes |
| (b) worst among 25 bp cells | 5.24% | 1.4071 | **no** |
| (c) **worst anywhere — adopted** | 5.55% | 1.4112 | **no** |

**(c) is adopted, and (a) is specifically rejected**, because the binding cell has the *smallest* spread of
all fourteen cells and only three reps. This arc already wrote the reason down: *"dispersion from n=2–5 is
biased **low**, not just imprecise — four cells looked tightest on few reps and turned out widest."*
Sizing the margin from the tightest 3-rep cell is that error exactly.

### 6c. The integers

`dyn_denominator(n_bins, cores) = σ₁₆·n_bins + k₁₆·cores`, today `23` and `64`.

| `σ₁₆` | `k₁₆` | `f_σ` | `f_k` | `B` cut | worst post-change ratio | clears (a)/(b)/(c) |
|---:|---:|---:|---:|---:|---|---|
| 32 | 89 | 1.391 | 1.391 | 28.1% | 0.965 | yes / **no** / **no** |
| **33** | **92** | **1.435** | **1.438** | **30.3%** | **0.939** | yes / yes / yes |
| 34 | 95 | 1.478 | 1.484 | 32.4% | 0.914 | yes / yes / yes |

**Take `σ₁₆ = 33`, `k₁₆ = 92`** — the cheapest pair clearing the adopted rule. Both terms scale together
(1.435, 1.438). The reason is **not** the `1/r` argument v3 gave, which rested on the withdrawn §2.3
decomposition: it is that the matrix contains only two `n_bins : cores` shapes, so `f_σ` and `f_k` are not
separately identifiable. A two-parameter solve is possible (`f_σ ≈ 1.12`, `f_k ≈ 1.65`, a 21% bin cut) but
strands the noisiest cell at ~3% headroom; **costed and rejected** in favour of uniform scaling.

All fourteen measured cells. `B` columns are **exact `bin_byte_budget`**, with the banner's floored integer
beside each — v3 divided the *banner* value by `f`, the exact thing Step 5 forbids:

| read len | cores | budget | `B` now | banner | `B` new | banner | `f` req | ratio now | ratio new |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 2 | 1G | 22.804M | 22 | 15.886M | 15 | 1.2574 | 1.201 | 0.903 |
| 25 | 4 | 1G | 18.126M | 18 | 12.623M | 12 | 1.2069 | 1.161 | 0.876 |
| 25 | 8 | 2G | 20.998M | 20 | 14.623M | 14 | **1.3370** | **1.300** | **0.939** |
| 25 | 16 | 2G | 10.499M | 10 | 7.312M | 7 | 1.2367 | 1.211 | 0.877 |
| 36 | 2 | 1G | 22.804M | 22 | 15.886M | 15 | 1.0964 | 1.075 | 0.816 |
| 36 | 4 | 1G | 18.126M | 18 | 12.623M | 12 | 1.1479 | 1.115 | 0.844 |
| 51 | 2 | 1G | 22.804M | 22 | 15.886M | 15 | 1.0030 | 1.002 | 0.765 |
| 51 | 4 | 1G | 18.126M | 18 | 12.623M | 12 | 1.0513 | 1.040 | 0.791 |
| 100 | 2 | 1G | 22.804M | 22 | 15.886M | 15 | 0.8257 | 0.864 | 0.669 |
| 100 | 4 | 1G | 18.126M | 18 | 12.623M | 12 | 0.8851 | 0.910 | 0.701 |
| 100 | 8 | 2G | 20.998M | 20 | 14.623M | 14 | 0.9802 | 0.982 | 0.718 |
| 100 | 16 | 2G | 10.499M | 10 | 7.312M | 7 | 0.9165 | 0.926 | 0.678 |
| 150 | 2 | 1G | 22.804M | 22 | 15.886M | 15 | 0.7113 | 0.775 | 0.607 |
| 150 | 4 | 1G | 18.126M | 18 | 12.623M | 12 | 0.7530 | 0.807 | 0.629 |

Headroom at the binding cell is **6.1%**, against that cell's own 0.65% spread and the worst spread
anywhere of 5.55%.

**Two promises, and only one is at risk.** Post-change the binding cell sits at 0.939 of the printed
prediction but **0.854 of `--memory` itself**, so the residual risk is a printed-prediction breach, not a
budget breach. That is the weaker of the two failures and worth stating rather than leaving implicit.

### 6d. The cost, stated plainly

**30.3% of bin size, for every input, at every read length.** v2's §6b claimed 19–20% for C; that came
from `s = 74` and is withdrawn. The cost rose because the real per-record cost is larger than `size_of`
suggested, not because the option changed.

**Long reads pay most in relative terms** — 150 bp ends at 0.607, so ~39% of the budget goes unused,
against 0.939 at the binding short-read cell. That is the flat-tax objection D existed to answer, and D is
not available.

**The compression cost is asserted nowhere in this plan, because Phase 4 cannot measure it.** The fixture
is **100.0% distinct sequences at every read length**, so it has no duplicates to cluster and its output
size cannot show a clumping gain or its loss. What the fixture *does* show is that an 18.2% bin cut moves
output by −0.11% to −0.24% — *downward* — which is evidence about all-distinct data only and must not be
quoted as "the cost is small". Validation 14 measures it on a real library. Until then the docs must not
claim a compression cost **or** its absence.

**The escape hatch is `--memory`, and the docs must give recomputed budgets, not a ratio.** `static` is a
fixed part of the budget, so scaling `--memory` by 1.435 over-shoots — the same error §9 names as its
sharpest lesson:

| to restore today's `B` at | needed | | |
|---|---:|---|---|
| 1G, `--cores 2` or `4` | **1363M** (+33.1%) | 2G, `--cores 8` | 2834M (+38.3%) |
| 4G, `--cores 2` | 5772M (+40.9%) | *a flat 1.435×* | *would say 1470M / 2938M — wrong* |

### 6e. Band coupling — v3's conclusion was inverted

v3 said: *"the band's own regression-detection power is unaffected, because it is set by width and by lever
A's magnitude, neither of which this change touches."* **Lever A's magnitude is exactly what this change
touches.** Lever A is the work-channel-depth revert; an extra queued batch is one whole bin per worker, so
its cost is `cores × σ_eff × B`, and **in ratio-of-prediction units it is proportional to `B`**.
`prediction` does not move and `B` falls by 1.435, so lever A's delta falls by 1.435 too:

| cell | half-range (27 obs) | → after | lever A delta | → after | margin at width 0.09 | at 0.065 |
|---|---:|---:|---:|---:|---:|---:|
| `1G_c4` | 0.0275 | 0.0192 | 0.126 | 0.0878 | **−0.0022** | +0.0228 |
| `25bp_1G_c4` | 0.0215 | 0.0150 | 0.128 | 0.0892 | **−0.0008** | +0.0242 |

At width 0.09 the margin goes **negative** — the width would exceed the regression it exists to catch, so
the band would read as armed while being blind to lever A. The noise side improves at the same time
(0.09 becomes 4.7× / 6.0× the half-range), which is exactly the room needed:

- **Narrow the width to 0.065 in the same PR.** 3× half-range is 0.0575 / 0.0450 ≤ 0.065 ✓, and the margin
  against a 0.0878 delta becomes 0.0228 — comparable to today's 0.036.
- **Narrow the down-side from −0.10 to ≈ −0.072** for the same reason.
- **Measure the post-change lever-A delta, don't infer it** (Validation 7). The 0.126–0.151 figures came
  from actually running the revert; the `B`-proportionality argument gives the direction, not the value.
- Re-baseline both cells from this PR's own runs; the baseline file's double-digit rule governs width, not
  centre. The `# BOOTSTRAP` marker already implements "band deliberately off, loudly".

**This compounds with Step 4c.4:** after the change nothing in the codebase ties `k` to `channel_depths`
either, so *both* of lever A's tripwires degrade in one commit.

### 6f. What is still not measured

v3 argued that `f`'s insensitivity to `static` licenses calibrating on one platform. **That is a non
sequitur** — insensitivity to a *constant* term says nothing about a platform-dependent *dynamic* term, and
Phase 2 measured σ = 1.29 Linux against 1.41 macOS, a 9% dynamic difference. Two valid arguments exist and
v3 made neither:

1. **The mechanism points the right way.** §12 attributes the excess to page-level fragmentation; darwin
   has 16 KiB pages against Linux's 4 KiB, so darwin should fragment *worse* and a darwin-derived floor
   should be conservative for Linux.
2. **The arc's own σ data agrees** — macOS 1.41 > Linux 1.29 at the same read length.

Neither is proof. **The rule stays Phase 2's: measure both, take the worse** (Validation 9).

| gap | consequence | blocks implementation? |
|---|---|---|
| **No Linux matrix** | integers calibrated on one platform | no — **blocks merge** |
| **Budget axis confounded with cores.** `f` varies 10.8% across the four 25 bp cells, non-monotonically in cores, and the highest `f` sits at the highest budget measured. No 1G/c8, 1G/c16 or 2G/c2 cell exists | the c8 result may be a budget effect | no — Validation 10, one 3-minute cell |
| **Linearity in `B` is questioned only downward.** Largest measured `B` is 22.8 MiB against 112.9 at 4G/c2, and the repo's own baseline file shows `1G_c4` 0.788 vs `2G_c4` 0.854 on Linux — **~9% super-linear upward**, flagged there as an open puzzle | large-budget cells may need more than `f` | no — Validation 8, and Validation 10's second budget |
| **Binding cell has 3 reps** | the margin rule uses the worst spread anywhere precisely because of this | no — Validation 5 raises it to ≥10 |
| **Shortest calibrated `r` is 96.642 B** | below it the bound is extrapolation | no — the warning, Step 4c.5 |
| **Cross-thread frees** (§12's open question) | part of the spike-to-pipeline gap | no — the constant comes from the pipeline |
| **No single-end cell** in 42 paired runs | constants are calibrated PE and applied to SE unchanged | no — Assumption 10 |

### 6g. Withdrawn

- **Option D** and the `text + spine = A` identity (§11); **the allocation-count variant** (§12).
- **`r_min` chosen by fitting `σ(r)`** — replaced by a measured `f`; `r_min` survives only as the warning
  threshold, at 97 B.
- **v2's §6b cost table** (19–20%) and **§6d band table**, both computed from `s = 74`.
- **v3's `(32, 89)`**, its 3.1% margin, its `B new` column, its "band unaffected" conclusion, and its
  one-platform justification.

## 7. Assumptions

1. **Fixed:** `prediction = budget × 10/11` stays the definition, so the change moves actual peak, not the
   printed prediction. Three integration tests assert the resulting banner verbatim (Validation 15).
2. **Fixed:** `--clumpify` requires `--cores ≥ 2` (`cli.rs:747`).
3. **Fixed:** integer arithmetic and `saturating_mul`. `B × dyn_denominator` stays inside `u64` at a
   `u64::MAX` budget under both old and new constants — no overflow is introduced.
4. **Fixed:** `size_of::<FastqRecord>()·2 + size_of::<MinimizerKey>() == 148` — **measured exactly** (§11)
   and pinned by Step 4c.6. What is withdrawn is the inference that 148 is the whole per-slot cost.
5. **Fixed:** `r` includes the header and must be measured at the record count actually used (§2.4a),
   confirmed to 0.011 B (§10).
6. **Fixed:** `f` is insensitive to `static`. `f = 1.435` tolerates `static` up to ≈434 MiB and the binding
   cell is stable across 173–450 MiB (§6b). This is what lets C be derived without identifying `static` —
   and **only that**; it says nothing about platform dependence of the dynamic term (§6f).
7. **Fixed:** the margin rule is `f_binding × (1 + worst relative within-cell spread anywhere)`, and the
   binding cell's own spread is explicitly **not** used (§6b).
8. **Open:** platform independence of the dynamic term. Only darwin/arm64 is measured. Resolution: take the
   worse of the two, as Phase 2 did (Validation 9).
9. **Open:** linearity in `B`, in **both** directions. Downward: `B` falls 30%, so flush-cycle count rises
   ~43% and §12 shows fragmentation tracks churn. Upward: the repo's baseline file shows ~9% super-linear
   behaviour between 1G and 2G (Validations 8, 10).
10. **Open, and unstated in v3:** the constants are calibrated **paired-end** and applied to **single-end**
    unchanged. The SE spine is 76 B/record against PE's 74 (one key per record rather than per pair), 2.7%
    apart, so the direction is safe — but no SE cell exists in 42 runs. `--clump_only` likewise inherits
    `f` with no cell of its own (Validation 13).
11. **Withdrawn:** `s = 74` as the per-record cost (§10, §11); `a = 0.937` (§12 — §2.4's own
    uniform-on-`(0,B)` argument gives 0.5, and 0.513 is measured); v1's assumption that no read-length
    floor need be chosen; and v3's claim that the `1/r` structure justifies uniform σ/k scaling (§6c gives
    the real reason — non-identifiability).

## 8. Implementation outline

### Steps 1–3 — done, see §10

1. ~~Repair the harness~~ — `phase4/calibrate.py`. **Carries a defect, fixed in Step 4c.1.**
2. ~~Run §4's matrix~~ — darwin only; Linux outstanding (Validation 9).
3. ~~Fit and confirm `s`~~ — **falsified**, which triggered §11 and §12 and settled §6.

### Step 4c — the change

**4c.1 — fix the harness before anything downstream depends on it again.** `calibrate.py`'s `design_row`
regresses on `b_banner`, the banner's **floored** MiB integer, because v3's harness dropped phase 2's
`b_exact_mib` column (`exact_b_mib`: 9 occurrences in `phase2/calibrate.py`, 0 in `phase4/`). At the c16
cell that is a ~5% error in the regressor. Restore an exact `bin_byte_budget` column and regress on it.
Consequences to record when it is re-run: §10's `a` shifts ≈−14%, `b` ≈+16%, one falsification cell flips
sign — so §10's "all in the unsafe direction" and §12's 1.22–1.75× ratios both need correcting.
**§6's integers are unaffected**: `f_required` is `(peak − static)/(prediction − static)` and never
references `B`. Do this **before** the Linux run so both platforms use the fixed harness.

**4c.2 — `src/clump.rs:174-176`, `dyn_denominator`: `23` → `33`, `64` → `92`.** The one line of arithmetic.

**4c.3 — the ten hard-coded expectations that move with it.** Five tests fail; recompute each against the
binary rather than from the table below, which is a cross-check only. Phase 3's PR 2 recorded five of six
pin values wrong when they were drafted from illustrative numbers.

| test | asserts | after `(33, 92)` |
|---|---|---|
| `resolve_layout_default` (`:574`) | `B ∈ [112, 114)` MiB at 4G/c2; comment names `denom = 496` | fails |
| `resolve_layout_scales_with_cores` (`:565`) | `B ∈ [44, 46)` MiB at 4G/c8; comment names `denom = 1248` | fails |
| `resolve_layout_pins_the_default_cell` (`:573`) | `B == 19_006_356` | fails |
| `sizing_constants_are_pinned_across_the_role_space` (`:583`) | six exact budgets | all six fail |
| `clumpify_floor_at_the_default_cell_is_pinned_absolutely` (`:632`) | `floor.div_ceil(MiB) == 281` | fails → 295 |
| `resolve_layout_fastqc_reserves_more_and_shrinks_bins` (`:670`) | **hard-codes the old denominator `624`** | fails |

**4c.4 — the test that stays green while its invariant becomes false.**
`channel_depths_sum_to_the_per_worker_coefficient` (`:641`) asserts `work + 1 + result_per_core == 4` with
the message *"must sum to k = 4 (64/16 in dyn_denominator)"*. `k` becomes `92/16 = 5.75`; the depths still
sum to 4 and the test still passes while guarding a relationship that no longer exists. **This is one of
lever A's two tripwires** (§6e), so it must not be left to rot. Either re-express the invariant against the
new `k₁₆`, or state in the test why the sum no longer equals `k` and what it now guards. The same falsehood
is in `channel_depths`'s docstring (`:177-178`, *"their sum is the constant and neither may move alone"*).

**4c.5 — the floor warning, at `r < 97 B`.** The calibrated floor is `r = 96.642 B`; v3's threshold of 96
left 0.64 B of calibrated range uncovered. Two corrections to v3's rationale:

- **`--small_rna` is the wrong justification.** The dispatcher bins **untrimmed** records
  (`parallel.rs:844-845`), so small-RNA input sits *above* the floor at input time. The real exposure is
  **re-processing already-trimmed FASTQ**, and terse headers: an SRA-style `@SRR1234567.1` puts ordinary
  25–36 bp data at `r = 68–90 B`, where the extrapolated requirement is ~1.56 against our 1.435 — a **~9%
  breach that is warned about but not prevented.** Say so in the warning text and in the docs.
- **`sanity_check_with` returns `()`**, so "reuse the record it already read" needs a signature change to
  hand the record back. Do that rather than opening the input again — v1's fresh sample added an I/O failure
  path to a function whose only current failure is `parse_memory_size`. One record is a weak estimator;
  note it in the warning's wording rather than implying precision.

**4c.6 — pin `size_of::<FastqRecord>()·2 + size_of::<MinimizerKey>() == 148` beside
`estimated_record_bytes`.** Called for in §2.3, §7 Assumption 4 and §11 rec 1; scheduled in no v3 step.

**4c.7 — the three stale docstrings and one stale comment.** `resolve_layout`'s docstring (`:273-286`)
names `σ = 23/16`, `k = 4`, the literal `23 × n_bins + 64 × workers`, and a provenance sentence
(*"fitted from a 10 M-pair paired-end matrix over cores ∈ {2,3,4} … take the worse of each"*) that becomes
false in every clause. `MARGIN_NUM`'s comment (`:263`) says *"largest observed within-cell spread is
96 MiB"*; the byte-equalised fixtures give 64.1 MiB. Keep each to one line, stating the current value —
not the one being replaced.

**4c.8 — arm the CI promise gates. This is the deliverable.** `ci.yml` passes `gating 0` for
`25bp_1G_c4` (`:444`) and `25bp_1G_c2` (`:449`), and the baseline file records twice that they are
un-gated *"until that refit lands"*. This is that refit. Move both to `gating 1` and delete the `#457`
comment, **in the same commit as the src change and the baseline re-record** — otherwise the fix ships with
no gate on the only cells it exists to fix.

**4c.9 — narrow the band to 0.065 and the down-side to ≈ −0.072** (§6e), re-baseline both cells from this
PR's runs, and **measure the post-change lever-A delta** rather than inferring it (Validation 7).

**4c.10 — re-measure and compare** (Validations 5, 8, 10): the 25 bp cells with ≥10 reps on the binding
cell, plus 25 bp c8 at 1G and one second budget. §6c's post-change column is a prediction until this runs.

### Step 5 — docs

**5.1** `docs/…/performance/clumpy.md`: worked values, the display maths (`:98`), examples (`:104`), the
floor table (`:114`) and the literal `need ≥ 311 MiB` (`:119`) all move. `modes/clump-only.md`'s `## Memory`
section too (`:85`), since it shares `dyn_denominator`. **Re-derive each against the binary.**

**5.2** State the long-read cost and give the `--memory` escape hatch as **recomputed budgets or a
formula** — 1G → 1363M (+33.1%), 2G/c8 → 2833M (+38.3%) — **never as a 1.435× multiplier**, because
`static` is a fixed part of the budget. Scaling a ratio with a fixed part is the error §9 names as its
sharpest lesson, and it would appear in shipped docs.

**5.3** Fix `clumpy.md`'s memory caution box, which currently exempts the 1 GiB default. Say nothing about
compression cost in either direction until Validation 14 measures it.

**5.4 CHANGELOG: one sentence**, cumulative about output bytes. This is the **third** change in one cycle
to move gzip member boundaries (#442, #446, now), so say it once for the cycle.

### Step 6 — independent of §6

**6.1 Cap the spine extrapolation** at **four** `reserve_exact` sites, not two: `parallel.rs:901`
(`PairedBin::push`) and `:1286` (`SingleBin::push`), plus `clump_only.rs:104` and `:147`, which duplicate
the same bin logic. `predicted` comes from one record's size, so a single short record after a flush
over-reserves the whole cycle. Derive it from a running mean of the previous cycle, or clamp it.

**6.2** §10 Finding 2: every Phase 4 fixture is fixed-length, so `bytes_this` is constant and **no
measurement here can see this defect** — it neither disturbs the calibration nor is validated by it. It
needs its own variable-length test.

## 9. Self-Review

**What v1 got wrong, and the pattern.** Four numbers, one cause: each was computed from the quantity I
had rather than the quantity that governed. `s` was fitted from a safety-rounded constant instead of
read from `reserve_exact`; `r` was measured at 2 000 pairs instead of 10 M; §7b scaled a ratio with a
fixed static part; §4 held pair count constant when total bytes was what mattered. Three of the four
erred toward making the plan look better.

**The sharpest lesson is about v1's §7b.** Its Self-Review *praised* that table for turning a caution
into a hard constraint "by computing it". Computing it is what stopped me checking whether the
computation's model was right — the static term does not scale, and holding it fixed inverts the
headline. The guard is cheap and I did not apply it: before trusting a computed comparison, list which
terms are constant under the change and which are not, and write that list down first.

> **v4 note — the rest of this subsection is v2's self-review, retained as the record and no longer
> current.** Its "Remaining risks" are five items about Option D, which §11 withdrew. Its "now verified"
> list below includes two claims that are now withdrawn: the `σ(r) = 0.937 + 74/r` decomposition and the
> `s = 74` per-record cost. The zero-parameter confirmation it cites is arithmetically correct but confirms
> a *component* (the 148 B spine), not the per-record cost. Read §10–§12 for what is actually established.

**What is now verified rather than assumed.** `size_of::<FastqRecord>()` = 72 (measured); the
`reserve_exact` / `mem::take` cycle (read); the decomposition reproducing 23 and 64 to two decimals; the
zero-parameter prediction at 0.6% / 3.2%; `r` at the 10 M scale; §6a's pre-screen from committed data;
and D's `text + spine = A` identity with its per-length table.

**Remaining risks.** (a) D's coefficient table is a first pass — only the identity is exact, and the
matrix must confirm the rest. (b) D raises peak RSS at long reads: defensible, but it needs saying in
the CHANGELOG rather than being discovered. (c) The occupancy terms may be platform-dependent
(Assumption 6), in which case the worse of the two is taken, as Phase 2 did for σ and `k`. (d) The
compression cost needs real libraries, deferred once already for the RRBS savings row — measure both in
one session. (e) D touches a hot path called once per record.

---

### v3 — what was checked, and the pattern behind two wrong constants

> **v4 note — historical.** This subsection reasons about `(32, 89)` and a 3.1% spread, both superseded
> in §6b/§6c. Its account of the `size_of`-as-a-floor error stands; its margin arithmetic does not.

**Applied the guard v2 said it had failed to apply.** §6b now lists which terms are constant under a
denominator scale *before* computing the comparison — `static` and `prediction` fixed, `B` and
`peak − static` scaling as `1/f`. That listing is what surfaced the result that makes C defensible on one
platform: `f` moves only 1.331 → 1.349 across a 173–260 MiB range of `static`, so C's constant does not
depend on the quantity the design cannot identify. Every figure in §6b/§6c is computed from
`results_phase4_darwin.tsv`, not derived by hand.

**Efficiency.** No hot-path change — that was D's cost, and D is gone. One integer pair in
`dyn_denominator`, and the floor warning reuses the record `sanity_check_with` has already read, so no
new I/O and no new failure path.

**Logic.** Both `σ₁₆` and `k₁₆` scale by the same factor because the breach rides `1/r`, which enters
`σ` and `k` alike (§2.3); scaling one alone would leave the other under-charging at exactly the cells
that breach. `(31, 86)` meets the bare requirement and was rejected for leaving 0.6% headroom against a
3.1% spread.

**Edge cases.** `min_required` carries `dyn_denominator` and so rises ~1.39× on its dynamic part — a
budget that worked before can now bail, which is user-visible and is Validation 9. Inputs below the
calibrated `r = 96 B` get a warning rather than a silent extrapolation. `--clump_only` shares the
denominator and moves with it. Long reads lose ~38% of their budget to a floor set by 25 bp, which is
real and belongs in the docs rather than in a user's bug report.

**Integration.** Band re-baseline in the same PR (§6e), docs re-derived against the binary rather than
scaled, CHANGELOG cumulative across #442/#446/this because it is the third gzip-boundary move in one
cycle, and Perl byte-identity asserted rather than assumed even though no validated flag path runs
`--clumpify`.

**The pattern behind both wrong constants, which is the lesson worth keeping.** v1 fitted `s = 28.87`
and was 2.6× low. v2 read `s = 74` off `size_of` and was 1.7× low. v2's own Self-Review named the
pattern — "each was computed from the quantity I had rather than the quantity that governed" — and then
repeated it, because `size_of` *is* an exact quantity and that made it feel like an answer. §2.4 even
wrote down that 72 B is a **floor** on the per-record term, one paragraph before §2.3 used it as the
value. A lower bound presented with certainty reads as a measurement.

What broke the cycle was not a better fit but a different instrument: measuring a **slope**, which has no
static term to trade against. Two fits on the same data gave `s = 108.74` or `s = 73.9` depending on
which parameters were held; the slope gave 148.0 for the spine with zero spread and then decomposed the
rest. **Where a design is degenerate, fit harder and you will get a plausible number for every prior you
bring; the fix is an experiment whose answer does not depend on the prior.**

**Remaining risks.** (a) One platform — Linux may need a larger `f`, and the rule is to take the worse.
(b) Linearity in `B` is assumed, and ~39% more flush cycles could move fragmentation (Validation 6).
(c) The long-read cost is the price of C and is not recoverable within C; `--memory` is the only lever.
(d) Step 6's spine cap is validated by nothing in Phase 4, because every fixture is fixed-length.

### v4 — what the dual review changed, and the one that changed the answer

**Both reviewers reproduced §6's arithmetic independently and it held** — all twelve `f` values, the
static-insensitivity table, the binding cell, every post-change ratio, and the `B` cut. One reviewer added
that `(32, 89)` was the *minimal* pair meeting v3's stated target. The derivation survived; everything
around it did not.

**The finding that changed the integers was a statistics error, not an arithmetic one.** v3's margin took
the largest *absolute* within-cell spread (64.1 MiB, from the c16 cell), divided it by *that* cell's dynamic
part, and applied the result to the *c8* cell's `f`. Two cells, one number. Under the correct relative
statistic the requirement is 1.411 and `(32, 89)` = 1.391 fails; `(33, 92)` clears at a 30.3% bin cut
instead of 28.1%. Worth noting *why* v3 got it wrong: 3.1% is a number that is easy to check and hard to
question, because recomputing it reproduces it exactly. Both reviewers did reproduce it; only one asked
whether it was the right quantity.

**The finding that reversed a conclusion was a missing proportionality.** §6e asserted the CI band's
detection power was untouched. Lever A is one extra whole bin per worker, so its magnitude in
ratio-of-prediction units scales with `B` — and `B` is the thing being cut. At width 0.09 the margin against
lever A goes *negative*. v3 reasoned about the band's width and noise floor and never asked whether the
signal moved.

**The finding that undermined the evidence was in my own instrument.** `calibrate.py` regressed on the
banner's floored MiB because the harness dropped a column phase 2 had. §10 exists because a fit could not be
trusted; the fit was also being fed a truncated regressor. The decision survives only because §11 measured
the per-slot cost directly and never depended on that fit — which is an argument for the slope method, not a
reason to relax about the fit.

**The pattern across all three: v3 was careful about the number it had chosen to defend and careless about
everything adjacent to it.** The margin, the band, and the harness were each treated as settled context
while the derivation got the scrutiny. A plan's weakest part is not its conclusion but whatever it stopped
questioning in order to reach one.

**Scope, honestly restated.** v3 called Step 4c "one line of arithmetic". It is one line plus five failing
tests, ten hard-coded expectations, a test that stays green while its invariant becomes false, three stale
docstrings, two CI gates to arm, a band width to narrow, a harness column to restore, four `reserve_exact`
sites, and a `size_of` pin that three sections asked for and no step scheduled.

**Remaining risks.** (a) Linux unmeasured — blocks merge, and §6f's mechanism argument gives direction only.
(b) Linearity in `B` unverified in **both** directions, with ~9% super-linearity already recorded upward in
the repo's own baseline file. (c) The binding cell has 3 reps and the smallest spread of fourteen; the margin
rule compensates, Validation 5 confirms. (d) The compression cost is unmeasured and unmeasurable on this
fixture. (e) Step 6's spine cap is validated by nothing in Phase 4.

## 10. Implementation notes

### Step 1 — harness repaired (done, 2026-08-19)

`phase4/calibrate.py`, verified on darwin/arm64 against `target/release/trim_galore` at
`3fa4de3`. Branch `perf/clumpify-spine-accounting`. **No source file is touched by this step.**

| §8 Step 1 requirement | how it is met | evidence |
|---|---|---|
| named-group parse | `BANNER` rebuilt against the shipped banner | banner captured from the binary *before* the pattern was written: `clumpify: 16 bins × 22 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 6)` |
| assert the parse produced values | `parse_banner` checks every group non-empty, `bins == max(16, 4·cores)`, banner budget == requested budget, `bin_mib`/`pred_mib` positive, and `pred == budget × 10/11` within 1 MiB | a pattern that matches but reads the wrong integers now fails; Phase 2 could not tell that from a healthy parse |
| prove it aborts on a broken pattern | `self_check` runs on the first cell's real output: the current pattern must parse it, and **both** Phase 2's rotted pattern and a mutated pattern must be rejected | smoke run printed `banner parses (16 bins, 22 MiB, pred 930 MiB, gzip level 6); rotted and mutated patterns both rejected` (§5 Verify 1) |
| `r` measured per fixture | one awk pass computes `len(id) + len(seq) + len(qual) + 5`, matching `clump::estimated_record_bytes`, and is compared against the analytic value for that pair count; a mismatch over 0.05 B aborts | measured == analytic at 25 bp and 51 bp; a one-byte header error moves it 1.000 B, so the gate has teeth |
| read length as an axis | 5 lengths × cores {2,4}, fixtures generated and deleted per length | `solve_pairs` gives 31,042,383 pairs at 25 bp and 8,673,732 at 150 bp for 3.0 GB/mate — §4's "~31 M" and "8.7 M" |

**The `r` model is confirmed against §2.4a independently.** `analytic_r(10 M, len)` reproduces all six
table values to 0.011 B, and the 2 000-pair deficit comes out at exactly 3.444 B against the plan's
"~3.4 B". Byte-equalisation is a fixed point, not a division, because `r` depends on the pair count
through the header — the §2.4a error class is now structurally impossible to repeat.

**The regression is verified against known answers, not against the smoke run.** The smoke fit returns
`a = −0.415`, `s = 1.35`: a 4 MB fixture never fills a 22 MiB bin, so RSS is static and the coefficients
are unidentified. Proving the solver therefore used synthetic cells generated from
`a = 0.937, b = 2.999, s = 74`: recovered exactly with static pinned, and with static fitted when the
design allows it. A mutated generator (`s = 40`) is recovered as 40, not 74 — the fit measures `s`
rather than asserting it.

**One bug found and fixed while writing it.** The distinct-sequence scan printed `$0` at the fourth line
of each block, which is the quality string, not the sequence. Every fixture would have reported one
distinct sequence — the same class of invisible fixture defect that invalidated the first #439
calibration. The sequence is now held from line 2. The `--distinct N` limiter also moved inside awk,
because `head` closing the pipe SIGPIPEs awk before the `END` block that prints the stats line.

### Deviation from §2.6: a new file, not an in-place repair

§2.6 says "repair it", meaning `phase2/calibrate.py`. Implemented as `phase4/calibrate.py` instead,
because `results_cal_darwin.tsv` and `results_cal2_linux.tsv` sit beside the Phase 2 script and were
produced by its `CELLS` table. Adding a read-length axis and byte-equalisation to that script would
leave two recorded TSVs that it can no longer reproduce. Phase 2 is left byte-identical; the new file
states its provenance in the docstring.

### Finding 1 — §4's matrix cannot identify the fixed term

At a single budget the design columns `1`, `B·n_bins` and `B·cores` all depend on `cores` alone
(`n_bins` is 16 throughout, and `B` is a function of cores and budget), so they span a two-dimensional
space and `static`, `a`, `b` are not separable. Only `s` is identifiable, because its column carries
`1/r`.

This is consistent with §6a — static is *held at its measured value*, not fitted — so the plan is
self-consistent as long as the fit pins it. The harness therefore pins `static` by default (205 MiB
darwin / 173 MiB linux, from `clump.rs:248-251`), prints the rank report either way, and **aborts rather
than returning junk coefficients** if asked to fit it on a degenerate design. Verified: the guard fires
on §4's matrix with `--fit-static`, and the four-parameter solve recovers static once a second budget is
added.

**Open for decision:** pinning static imports Phase 2/3's measurement as an assumption, and any error in
it biases `a` and `b`. Adding one more `--budgets` entry (e.g. `1024,2048`) makes it a measured quantity
instead, at roughly double the runtime. §4 as written says 1G only, so this is offered, not taken.

### Finding 2 — the matrix is blind to Step 8, in both directions

§8 Step 8 (cap the spine extrapolation) changes `predicted` from one record's size to a running mean.
On `mk_clump_fixture` every record is the same length, so `bytes_this` is constant and the cap changes
nothing. Two consequences the plan does not state:

1. **Ordering does not matter.** Step 8 can land before or after the matrix without invalidating the
   fit, so it is not a sequencing constraint on the decision.
2. **The fit cannot see the defect it is fitting around.** The occupancy terms are measured on
   fixed-length reads only, so they say nothing about variable-length input — which is every real
   library after trimming. Step 8's over-reservation is invisible to every cell in §4.

### Status

Step 1 done. Step 2 ran on darwin (below); Step 3's gate tripped, which is what sent the question to the
spikes in §11 and §12. §6 is now decided (C) and Step 4c is ready. Linux cells and the post-change
re-measurement are outstanding and both belong to the implementing PR.

### Steps 2–3 — matrix run on darwin/arm64, and the Step 3 gate is TRIPPED

42 runs, 10 fit cells (1G, cores {2,4}, read lengths {25,36,51,100,150}, byte-equalised at 3.0 GB/mate)
plus 4 falsification cells (2G, cores {8,16}), 3 reps each, exact distinct scans. ~115 min wall.
Raw: `results_phase4_darwin.tsv`. Static pinned at 205 MiB.

**§8 Step 3 says: stop and report if `s` does not confirm within a few percent. It does not.**

```
a (per bin)            =   0.786   plan 0.937
b (per worker)         =   2.777   plan 2.999
s (B per record-byte)  = 108.74    size_of-derived 74      (+47.0%)
```

`s = 108.74` means **217.5 B per paired slot against the 148 B** that `size_of` gives
(72 + 72 + 4). The excess is 69.5 B per pair, 34.7 B per record. §2.3's *functional form* survives;
its *value* does not.

**The measured defect, which is the result worth keeping.** 8 of 14 cells breach the promise. Seven of the
eight breach on **every** rep; **51 bp c2 breaches on 1 of 3** (0.9878, 0.9953, 1.0024) and is the single
marginal cell. Breach magnitudes run **2.2 to 558.1 MiB**, not the "150–250" v3 claimed, and 51 bp c2's
2.2 MiB breach is *below* that cell's own 13.6 MiB spread — the other seven all exceed theirs:

| read length | c2, 1G | c4, 1G | c8, 2G | c16, 2G |
|---:|---|---|---|---|
| 25 bp | **1.201** | **1.161** | **1.300** | **1.211** |
| 36 bp | **1.075** | **1.115** | — | — |
| 51 bp | 1.002 | **1.040** | — | — |
| 100 bp | 0.864 | 0.910 | 0.982 | 0.926 |
| 150 bp | 0.775 | 0.807 | — | — |

Ratio falls monotonically in `r` and crosses 1.0 at 51 bp — exactly where the shipped constants were
calibrated (`r ≈ 147.9`). That is the `1/r` signature §2.3 predicts, and it confirms the defect's
direction and size at scale.

**§4 Fault 1 is confirmed quantitatively.** Phase 3 measured 25 bp/c2 at 1.02–1.06 on 10 M pairs; the
byte-equalised cell (31 M pairs) reads 1.149–1.201. The old fixed-pair-count fixtures understated the
breach by ~0.12, in the unsafe direction, on the axis being fitted. (Consistent-with, not proof: fixture
scale and pair count moved together.)

**§5 Verify 4 also fails.** The falsification cells are under-predicted by the fit at +8.4%, +0.9%,
+16.0%, +9.3% — all in the unsafe direction. The occupancy terms do not transfer to 8/16 cores.

#### What is robust, and what is not

| quantity | robust? | evidence |
|---|---|---|
| `s ≈ 109`, not 74 | **yes** | identical to 2 dp at every static pin from 173 to 300 MiB, and `corr(B·n_bins, B·(n+2c)/r) = +0.022`, so `a` and `s` do not trade off |
| `a = 0.786`, `b = 2.777` | **no** | hostage to the pinned static: `a` runs 0.857 → 0.576 as static goes 173 → 300 MiB. Finding 1 biting — these need a second budget to be identified |
| curvature in `1/r` | **not established** | adding a `1/r²` column moves residual RMS 19.0 → 15.8 MiB, against a largest within-cell spread of 47.8 MiB |
| the linear-in-`1/r` form | survives | residual RMS is already below the measurement noise |

**The one alternative that would rescue `s = 74`.** Forcing `a` and `b` to their structural values leaves
only `s` free, and then `static = 260 MiB` returns `s = 73.9`. So "the excess is really an
under-pinned static term" is not yet excluded by this design — it is the same static/`a`/`b`
degeneracy as Finding 1, seen from the other side. A direct measurement has to settle it; the fit cannot.

#### Leading hypothesis for the 69.5 B/pair — unverified

`.to_string()` on a trimmed slice allocates exactly `len`, so `String` capacity is not the cause. But the
binary runs **mimalloc** (`main.rs:29`), which rounds every allocation to a size class, and there are
**three `String` allocations per record**. Neither `size_of` (inline struct bytes) nor
`estimated_record_bytes` (text length) can see that rounding. A 3-bit-mantissa bin model gives 24–62 B
per pair against the 69.5 B needed — the right order, but the bin table is a guess and the modelled waste
varies with read length (24 B at 25 bp, 62 B at 100 bp) where the data prefers a constant. **Do not
build on this until it is measured.**

Note the corollary if it holds: the *text* side is under-counted the same way, because `raw_bytes`
accumulates text lengths while resident heap holds rounded-up capacities.

#### What this does to §6e

**It undercuts D specifically, and leaves C's logic intact.** D's whole claim is the exact identity
`text + spine = A` with a `size_of`-derived 148 B per slot. If the real resident cost is ~217 B, then
charging 148 in the threshold leaves ~70 B/pair uncharged and the promise still breaches at short reads.
D can charge the measured value instead — but then the constant is *fitted*, which removes D's headline
advantage ("`r` leaves the model", no fitting), and if the excess is allocator rounding it varies with
read length, so no single constant is exact "at every read length" as §6c claims. C is a conservative
bound rather than a model, so a larger true per-record cost makes C *more expensive* without making it
wrong.

**Recommended re-plan, in order:**

1. **Measure the per-slot resident cost directly** — a small example that builds a `Vec<FastqRecord>` of
   known text size at each read length and reports the RSS delta. Deterministic, seconds to run, and it
   separates "148 B struct + text" from allocator overhead without a 2-hour campaign. This is the
   experiment that decides whether `s ≈ 109` or `static ≈ 260`.
2. **Measure static** rather than pinning it, by adding a second `--budgets` entry (Finding 1). Until
   then `a` and `b` are not quantities, and every cost figure that uses them inherits that.
3. **Only then re-open §6e.** Both options' arithmetic is computed from `s`, so both tables are stale.

**Not done, deliberately:** no source file is touched, and Step 4 is not started under either branch.

---

## 11. Spike Results — `spikes/SPIKE_slot_cost.md`

**Validates:** §8 Step 3's stop-and-report gate, and the premise under §2.3 that both §6b and §6c
compute their costs from.

**Outcome: NEEDS REVISION.** §2.3's `size_of` reasoning is confirmed exactly and is also incomplete.

| measurement | value | against |
|---|---|---|
| spine only (three empty `String`s per record) | **148.0 B/pair at every read length**, spread 0.0 | `size_of` says 148 — **confirmed exactly** |
| allocator cost, per `String` allocation | **17.8–21.9 B** (six per pair, so 107–132 B/pair) | invisible to both `size_of` and `estimated_record_bytes` |
| `s`, as `clump.rs` uses it | **122.5–134.8 B/record**, run-to-run spread 0.02 | §2.3's 74 **excluded**; the matrix's 108.74 **excluded** |
| allocator choice (libmalloc vs mimalloc) | 1–4 B/pair | not the driver — this was the leading hypothesis |
| variation across read length | 12.3 B/record (~10%) | no single constant is exact at every length |

**Method note.** The spike measures a *slope* in record count, so the process baseline cancels. That is
what makes it able to answer a question the matrix structurally cannot: with `a`/`b` free the matrix
gives `s = 108.74` at every static pin, and with `a`/`b` pinned to structural values `static = 260`
returns `s = 73.9`. A slope has no static term to trade against.

**Why the matrix's 108.74 and the spike's 122–135 are both right.** The spike holds everything resident.
In production the spine is resident in full from the first push after each flush, while the text — *and
the ~120 B/pair of allocator overhead attached to it* — oscillates over `0..B`. The matrix's single `s`
blends a fully-resident 148 B with a partially-resident ~120 B and lands between them.

**Root cause, restated.** The bin flushes when `raw_bytes` reaches `B`, and `raw_bytes` sums
`estimated_record_bytes` — text. Resident is that text, plus 148 B/pair of spine, plus ~120 B/pair the
allocator takes. At 25 bp the text term alone understates resident heap by 59%. Short reads breach
because they buy more *allocations* per byte, which is a sharper statement than "more records per byte".

### Recommended changes to this plan

1. **§2.3** — keep the derivation, correct the conclusion. `s = 74` is the spine, not the per-record
   cost. Assumption 4 ("`s = 74` is `size_of`-derived, not fitted — pinned in code") is **withdrawn**;
   the pinnable invariant is `size_of::<FastqRecord>()·2 + size_of::<MinimizerKey>() == 148`, which is
   worth pinning on its own terms and is now measured to hold.
2. **§6c — Option D as written is not viable.** "`text + spine = A`, exactly, at every read length"
   fails: charging 148 leaves ~45% of the per-slot cost uncharged, and the remainder varies 10% by read
   length so no constant is exact at more than one. D can charge ~268 B, but then the constant is
   measured rather than derived and D loses its headline advantage.
3. **§6b — Option C's logic survives; its numbers are stale.** C is a bound, not a model, so a larger
   true cost makes it more expensive without making it wrong. Recompute the floor against `s ≈ 135`
   (worst read length) instead of 74.
4. **A third option should be costed before §6 is re-opened:** charge the **allocation count**, not
   bytes. Resident is `text + 148·slots + ~20·3·records`, and the dispatcher already knows the record
   count it is about to reserve. This is D's shape — account rather than fit — aimed at the quantity
   measured to be constant rather than the one `size_of` happens to expose.
5. **Assumption 7 is resolved:** the spine cost is linear in record count, not stepped by page
   granularity — the slope is exact to 0.0 B/pair across a 4× range of record counts.
6. **New open question, from Finding 6:** does the allocator overhead ride the text's `0..B` oscillation,
   or does the allocator retain freed pages so it behaves as fully resident? That decides whether the
   charge is `120·slots` or `a·120·slots`, and it is the next measurement worth taking.

## 12. Spike Results 2 — `spikes/SPIKE_occupancy.md`

**Validates:** §11's open question — does the ~120 B/pair allocator overhead ride the text's `0..B`
oscillation, or is it resident in full?

**Outcome: REJECTED, both candidates.** The charge is `slots` **plus a read-length-dependent
fragmentation term**.

| measurement | value | reading |
|---|---|---|
| peak at 8 / 16 / 32 flush cycles | 1973.0 / 1973.0 / 1974.0 MiB | freed pages **are** reused; 4× the allocations, same peak. Steady state at 8 cycles |
| peak / one full pool | 1.33 (25 bp), 1.12 (51 bp), 1.04 (150 bp) | fragmentation **on top of** full residency, worse as `r` falls |
| text occupancy `a` | **0.513** | the mean of uniform on `(0, B)` — what random bin assignment plus flush-at-`B` gives |
| per-record charge, cycling | 266.9 B | against 129.1 B measured fully resident in §11 |
| `σ` vs the production 10-cell fit | 1.75× / 1.44× / 1.22× at 25 / 51 / 150 bp | shape and sign match, **level does not transfer** |

**Mechanism.** Allocations from different bins interleave in time, so one mimalloc page holds blocks
from many bins. When a bin flushes, its blocks are scattered across pages other bins still occupy, so
those pages neither release nor reuse cleanly. Shorter reads mean more allocations per byte, so more
scattering.

**§2.4 is marked `[REVISED BY SPIKE]`.** It argues text occupancy is "~uniform on `(0, B)`" and then
concludes `a = 0.937`. Uniform on `(0, B)` has mean **0.5**, and the isolated measurement gives 0.513.
The mechanism in §2.4 is right; the value 0.937 must be absorbing something else, and the two cannot both
be quoted as they stand.

**The magnitude does not transfer and must not be copied into `clump.rs`.** The isolated model is
single-threaded, drops batches at flush instead of handing them to a worker, and does not exercise
mimalloc's thread-local heaps — cross-thread frees being the classic fragmentation case. Attribute the
mechanism from the spike; take every constant from the pipeline.

### Recommended changes to this plan

1. **§6 should be re-opened with C as the presumptive answer.** The reason D died now arrives a second
   time by a different route: **no per-record constant is correct at more than one read length.** D
   failed because the spine is ~55% of the per-slot cost; the "third option" from §11 (charge the
   allocation count at a measured constant) fails because the fragmentation term sitting on top of that
   count is itself read-length-dependent. C is a conservative bound, so both findings only move its
   numbers.
2. **§11 recommendation 4 is withdrawn** — costing the allocation-count option is no longer worth doing.
3. **Recompute C's floor against measured `σ` at the shortest supported `r`,** taken from the pipeline
   (the 10-cell fit, or a re-run with a second budget per Finding 1), never from either spike.
4. **New open question:** how much of the 1.22–1.75× gap is cross-thread frees? It is the one mechanism
   present in production and absent from both spikes, and it is the remaining reason a pipeline-measured
   constant could still move under `--cores`.

---

## 13. Implementation notes — v4 Step 4c (2026-08-20)

Branch `perf/clumpify-spine-accounting`, off `origin/dev` @ `3fa4de3`.
**Verification: `cargo fmt` clean, `cargo clippy --all-targets --release -- -D warnings` clean,
805 tests pass.** Both gates were proved to fire before being trusted: a non-underscore unused
variable reds clippy, and the new spine-cap test fails with the fix reverted.

### What was implemented

| step | change | note |
|---|---|---|
| 4c.1 | `calibrate.py`: `exact_bin_mib()` added, `b_exact_mib` column restored, `design_row` regresses on it | banner-vs-exact error measured at **4.3–5.6%**, matching the review's ~5% |
| 4c.2 | `dyn_denominator`: `23` → `33`, `64` → `92` | one line |
| 4c.3 | six pinned expectations recomputed **off the binary** | see below — the plan's cross-check was wrong by 1 |
| 4c.4 | `channel_depths_sum_to_the_per_worker_coefficient` → `channel_depths_are_pinned_to_the_calibrated_values` | asserts the depths, not the false arithmetic |
| 4c.5 | floor warning at `r < 97 B`; `sanity_check_with` now returns the first record's size | verified live: a terse-header 25 bp record is 68 B and warns |
| 4c.6 | compile-time pin: `2·size_of::<FastqRecord>() + size_of::<MinimizerKey>() == 148` | a `const` assert, so it cannot drift silently |
| 4c.7 | 3 docstrings + `MARGIN_NUM`'s stale "96 MiB" | `resolve_layout`'s provenance paragraph rewritten |
| 4c.8 | `ci.yml`: `25bp_1G_c4` and `25bp_1G_c2` → `gating 1`; the `#457` comment deleted | **the deliverable** |
| 4c.9 | band width `0.09` → `0.065`, down-side `−0.10` → `−0.072`, in the enforced `awk` and the prose | baseline set to `# BOOTSTRAP` |
| Step 5 | `clumpy.md` + `clump-only.md`: formula, layout examples, floor table, both error messages, the caution box, and a new read-length section | every figure re-derived; the two error messages **captured from the binary** |
| Step 6 | spine cap at **all four** `reserve_exact` sites, plus two tests | one of them mutation-proved |

### Two things the plan got wrong, found by implementing it

1. **The absolute floor is 296 MiB, not the 294/295 the plan's cross-check table projected.** Step 4c.3's
   instruction to read pins off the binary rather than the table is what caught it. The other five pins
   matched the formula exactly, so the table was right about the mechanism and wrong about one rounding.
2. **`Ok(0)` after the exhaustive match in `sanity_check_with` is unreachable**, so the first attempt at
   4c.5 introduced a warning that `-D warnings` would have failed CI on. The match had to become the
   returning expression.

### Measured after the change — better than predicted

Re-running the two worst cells on darwin/arm64 with the same byte-equalised 31 M-pair 25 bp fixture:

| cell | before | §6c predicted | measured, 3 reps | worst | headroom |
|---|---:|---:|---|---:|---:|
| 25 bp c8 2G (**binding**) | 1.300 | 0.939 | 0.923, 0.909, 0.907 | **0.923** | 7.7% |
| 25 bp c2 1G | 1.201 | 0.903 | 0.821, 0.807, 0.815 | **0.821** | 17.9% |
| 25 bp c4 1G | 1.161 | 0.876 | 0.789, 0.816, 0.813 | **0.816** | 18.4% |

All nine reps clear the printed prediction; all nine breached it before. Relative within-cell spreads are
1.7%, 1.7% and 3.3%, against the 5.55% the margin rule in §6b assumed — so the rule remains conservative
on post-change data as well as pre-change.

**Validation 5 passes at the binding cell** — every rep under the printed prediction, where every rep was
over it before.

**Validation 8 (linearity in `B`) is answered, and the answer is "approximate and favourable".** All three
cells came in *below* prediction, so treating `peak − static` as scaling exactly `1/f` is conservative.
But they do not reconcile to a single `static` under a pure `1/f` scaling — back-solving gives ≈110 MiB
from the c8 cell and a negative value from the c2 cell — so the relationship is **not** linear. Note the
shape: the miss is smallest at 2G (−0.016) and largest at 1G (−0.082, −0.060), i.e. the smaller the bin
the more the prediction overshoots. That is consistent with §12's mechanism, where a smaller bin holds
fewer live allocations and fragments less — a second-order gain the `1/f` model does not carry. Recorded
as a conservative approximation, not a verified identity, and it means §6c's post-change column should be
read as an upper bound rather than an estimate.

### Not done, and why

- **Validation 7 (post-change lever-A delta) and 9 (Linux)** need CI. The band is deliberately in
  `# BOOTSTRAP`: the two committed centres are pre-refit and marked as untrustworthy, with the predicted
  values written in as a sanity check and explicitly *not* as baselines. First CI run records the real
  ratios; then replace both rows and delete the marker. Deleting it without re-recording reds the band
  immediately rather than passing silently, which is the safe failure.
- **Validation 10** (25 bp c8 at 1G, to separate budget from cores) and **Validation 14** (output size on a
  real RRBS/WGBS library) are measurements, not code, and neither gates the change.
- **§10 and §12's coefficients are not recomputed.** 4c.1 fixed the harness; re-running the fit is a
  measurement task. §6's integers never depended on it.

### §14 — Review recommendations applied (2026-08-20)

Dual code review (`CODE_REVIEW_v4_A.md`, `CODE_REVIEW_v4_B.md`) plus the coverage audit
(`COVERAGE_v4.md`, verdict INCOMPLETE) produced 2 MISSING, 5 PARTIAL and 24 review findings. **All are
applied.** Neither reviewer found a Critical.

**The floor warning was redesigned, not patched.** Five findings converged on it: it fired on the
calibration fixture itself (the first record's index is one digit, so `r` = 90 B against the fixture's
96.64 mean — it would have fired on both CI cells this PR arms), it skipped pairs 2..N, it printed before
the prediction it referenced, its `== 0` guard was unreachable, and it made `fastq.rs` reach into
`clump.rs`. The criterion is now **read length, not record bytes**: the calibrated floor is 25 bp, and
header bytes only *raise* the per-record size, so they make a run safer (§2.4a). That kills the false
positive at the floor, removes the header dependence, and let `sanity_check_with` revert to `Result<()>`.
It checks every input, is emitted after the banner, and covers `--clump_only` — which never reaches
`resolve_clump_layout` and so needed its own call. Four integration tests in
`tests/integration_clumpify_floor_warning.rs`, one of which is specifically the regression guard for the
old byte-based form.

**Residual `Vec`-doubling exposure removed** (B M6). `reserve_exact` landing within ±1 slot of the true
count means any downward drift in mean record size overflows it and doubles the allocation — ~half of
cycles on variable-length input, costed at +0.118 of prediction, wider than the 0.065 band this change
narrows to, and invisible to every fixture here because they are all fixed-length. Now
`predicted + predicted/16 + 1` at all four sites: 6.25% of spine against a 100% overshoot.

**Measured cost of that headroom**, 3 reps per cell:

| cell | without headroom | with headroom | worst-rep cost |
|---|---|---|---:|
| 25 bp c8 2G (binding) | 0.923 / 0.909 / 0.907 | **0.927 / 0.918 / 0.911** | +0.004 |
| 25 bp c2 1G | 0.821 / 0.807 / 0.815 | **0.808 / 0.838 / 0.838** | +0.017 |

The promise holds on every rep; the worst is 0.927, i.e. **7.3% headroom**. One thing to watch rather than
bury: the c2 cell's relative within-cell spread roughly doubled, 2.3% → **4.98%**, which is still inside
the 5.55% the §6b margin rule assumed but no longer comfortably. If the Linux run shows the same
widening, the margin rule's input needs re-deriving from post-change data rather than pre-change.

| finding | resolution |
|---|---|
| CHANGELOG (3× flagged) | new `#457` entry; the `#439` entry's claim that short reads are uncovered **corrected** |
| Validation 12's figures | 52 → **50**, 29 → **28**; the plan now records that v4 copied them from a review instead of recomputing |
| 4c.1 partial | the rank diagnostic, zero-parameter check and falsification block all moved to `b_exact` |
| 4c.9 partial | seven `+0.09` references reconciled; `1G_c2`'s verdict inverts (now clears the noise bar) and `25bp_1G_c2`'s tightens |
| bootstrap note (B M5) | now asks for **both** sides of the width rule — new centres *and* a measured lever-A delta |
| `1G_c2` lost its rationale | restored, on the real reason (tightest range of the six) |
| docs quoted a warning the binary does not print | the binary now emits two lines in the house style; the docs quote it verbatim |
| docs stated predictions as observations | 25 bp given as measured (0.81–0.92), 150 bp labelled an extrapolation |
| docs omitted the newly-falling-back config | `--fastqc` at `--cores 24`–`33` at the default budget named explicitly |
| `148` pin on 32-bit | `#[cfg(target_pointer_width = "64")]` — the figure is 64-bit |
| exact-capacity assertion | bounded rather than pinned; `reserve_exact` is documented as free to over-allocate |
| three of four sites untested | tests for `SingleBin` and both `clump_only` bins |
| test comment carried its own rationale | trimmed to the claim |
| §6d `2833M` | **2834M** — the docs had it right and the plan wrong |

**Validation 16 (Perl byte-identity) is a CI-only check** and was missing from §13's deferral list; it is
deferred with Validations 7, 9, 10 and 14.

**Verification after all of it:** `cargo fmt` clean, `cargo clippy --all-targets --release -- -D warnings`
clean, **810 tests pass**. Clippy caught one lint in the new test file (`ptr_arg`) that `cargo test`
compiled happily — a reminder that `cargo test` is not the gate CI runs.

### §15 — CI result, PR #474 (2026-08-20)

Run `32341837369`, commit `467865a`: **`completed/success`, all 10 jobs.** Verified on `conclusion`, not
on the watch command returning.

**Validation 9 (Linux) — passes, and Linux is not the worse platform.** Recorded ratios, one observation
each:

| cell | pre-refit Linux | post-refit | darwin, worst of 3 |
|---|---:|---:|---:|
| `1G_c4` | 0.788 | **0.636** | — |
| `1G_c4_fastqc` | 0.7069 | **0.541** | — |
| `25bp_1G_c4` (now gated) | 0.914 | **0.683** | 0.816 |
| `multipair_c4` | 0.7987 | **0.608** | — |
| `1G_c2` | 0.8789 | **0.648** | — |
| `25bp_1G_c2` (now gated) | **1.0301** | **0.748** | 0.838 |
| `2G_c4` | 0.854 | **0.608** | — |

`25bp_1G_c2` was the cell the baseline file named by issue number as #457's evidence — over the printed
prediction, inside the budget, on Linux. It now reads 0.748 and gates. Linux comes in 0.09–0.13 *below*
darwin at the same cells, so **the worse-of-two rule keeps darwin binding and `(33, 92)` stands
unchanged.** That is also the direction §6f's page-size argument predicted (darwin 16 KiB vs Linux 4 KiB,
so darwin fragments worse), which is weak corroboration of the mechanism rather than proof.

**Validation 16 (Perl byte-identity) — passes.** `Validate vs Perl TrimGalore` is green, so the md5
comparison against Perl 0.6.11 is unaffected. It was listed as deferred in §14; it is now done.

`Reproducibility (SOURCE_DATE_EPOCH)` and `Lint (fmt + clippy)` also pass, so the build stays
bit-identical and the local clippy result holds on the CI toolchain.

**The predicted centres were close but not committable.** §6f/§14 predicted ~0.62 and ~0.71 against
measured 0.636 and 0.683 — inside the 0.065 band, but the second is a third of the band's width out.
The baseline marked them "NOT to be committed as a baseline"; that was the right call.

**The band stays in `BOOTSTRAP`.** Real centres now exist, but the file's own arming rule needs two
things and only one is available: these are **n = 1**, and the 0.065 width was *inferred* by scaling
pre-change half-ranges by `1/f` rather than measured. The file's history is the argument against arming
on that basis — the width was widened to 0.09 in the first place because a tighter band false-fired on
another session's PR. The promise gate, which is the deliverable, is enforcing on all four gated cells
today.

**Still outstanding:** Validation 7 (post-refit lever-A delta, which needs a deliberate depth revert),
Validation 10 (25 bp c8 at 1G, to separate the budget axis from cores), Validation 14 (output size on a
library with duplicates), and a double-digit sample per banded cell.
