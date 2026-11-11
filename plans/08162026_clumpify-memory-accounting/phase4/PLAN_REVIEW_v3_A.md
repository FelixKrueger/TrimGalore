# Plan Review — Phase 4 v3 (#457), Reviewer A

**Target:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v3
**Repo state:** `perf/clumpify-spine-accounting` @ `3fa4de3`, no source modified (verified: `cargo test --lib clump::` → 25/25 pass)
**Raw data:** `results_phase4_darwin.tsv` (42 runs), `spikes/SPIKE_slot_cost.md`, `spikes/SPIKE_occupancy.md`
**Scope note:** Option C is the user's decision. This review does not re-argue C vs D. It audits whether §5–§8 are supported by §10–§12 and by the code.

---

## 0. Verdict

**The arithmetic in §6 is sound and reproduces from the TSV.** I recomputed every number in §6b and §6c
independently; the `f_required` derivation is algebraically correct, the static-insensitivity table is
exact, all 12 per-cell `f` values match to 4 dp, and every entry in §6c's `ratio new` column is exactly
reproducible. That is a real improvement on v1/v2 and the guard §9 says it applied did work.

**The gaps are all downstream of §6, in §5 and §8.** Four of them I would treat as blocking:

1. **The change silently un-clumps real configurations.** At the default `--memory 1G` the maximum
   `--cores` that keeps `--clumpify` engaged falls **72 → 52** without FastQC and **33 → 23** with it.
   `--cores 24 --fastqc --memory 1G` clumps today and falls back to plain mode after. §8 Step 7 tests
   `--cores 2` and `16` — both far from the boundary — and describes the failure mode as "bails", which
   is wrong: `main.rs:190` warns and falls back to plain mode.
2. **Three existing unit tests will fail and are not named in the plan**, plus a fourth pin the plan's
   count misses.
3. **The CI promise gate is never armed on the 25 bp cells** — the deliverable of #457. `ci.yml`
   passes `gating=0` for `25bp_1G_c4` and `25bp_1G_c2`, and the baseline file says in as many words that
   they are ungated "until that refit lands". §6e and §8 Step 8 only re-baseline the *drift band*.
4. **`channel_depths_sum_to_the_per_worker_coefficient` stays green while asserting a now-false
   invariant.** k goes 4 → 5.5625, which no longer decomposes into the channel slot counts the test and
   its docstring tie it to. A green test asserting a false invariant is the worst available outcome and
   is exactly the failure that test was written to prevent.

Below that: §6c's `B new` column violates the plan's own Step 5 rule (it scales the banner-rounded MiB
instead of recomputing — every row is 3.5–5% wrong), the 3.1% margin is the wrong statistic, §6f omits
the two largest un-measured axes (budget scaling and the `f`-varies-across-cells spread), the §6d cost
claim is contradicted by the plan's own `out_bytes` column and by the shipped docs, and §9's pre-`---`
block is stale v2 prose about Option D.

---

## 1. Reproduction of §6 — what checks out

Recomputed from `results_phase4_darwin.tsv`, worst rep per cell, `f = (peak − static)/(prediction − static)`.

### §6b static-insensitivity table — exact

| `static` | plan | recomputed |
|---:|---:|---:|
| 173 | 1.3306 | 1.33063 |
| 205 | 1.3370 | 1.33702 |
| 224 | 1.3409 | 1.34088 |
| 260 | 1.3486 | 1.34862 |

Binding cell is 25 bp / c8 / 2G, rep 2, RSS 2419.1 MiB — correct, and it is the max of that cell's three reps.

### §6b per-cell table — all 12 exact

1.3370 / 1.2574 / 1.2367 / 1.2069 / 1.1479 / 1.0964 / 1.0513 / 1.0030 / 0.9802 / 0.8257 / 0.7530 / 0.7113
all reproduce to 4 dp. **Two measured cells are missing from the table**: 100 bp c4 (0.8851) and 100 bp
c16 (0.9165). Both non-binding, so nothing changes, but a table headed "per-cell requirement" that
silently drops 2 of 14 cells invites the reader to wonder what else was dropped. Add them.

### §6c ratio-new column — exact

Every row reproduces. Denominator scaling is genuinely uniform across the design, which is worth stating
explicitly because it is not obvious: `(32·n_bins + 89·c)/(23·n_bins + 64·c)` = 1.3911 at c2,
1.39103 at c4, 1.39103 at c8, 1.39103 at c16. So `f` is a clean scalar on the whole denominator.

### The algebra is correct

`prediction_new = static_charged + (B/f)(f·denom)/16 = prediction_old`, so the prediction really is
invariant. Requiring `static_actual + dyn/f ≤ prediction` gives exactly
`f ≥ (peak − static_actual)/(prediction − static_actual)`. Note the relevant variable is *actual* static,
not *charged* — the plan gets this right, and it is why the charged 224 MiB never appears.

**The insensitivity claim is real but its cause should be stated**, because it is a bound, not a law:
`df/dstatic = (peak − prediction)/(prediction − static)²`. Insensitivity is proportional to the size of
the breach. At the binding cell that is 558/1656² = 2.0e-4 per MiB, hence 0.018 over an 87 MiB range.
Written that way a reader can tell when it would stop holding (a cell breaching by 2× would be ~4×
more sensitive), which the current "the result is insensitive to `static`" does not convey.

### One free corroboration the plan does not use

§12 measured isolated `σ(25 bp) = 3.342` and reports the isolated-to-pipeline gap as 1.75× at 25 bp,
implying a pipeline `σ(r=96.6) ≈ 1.91`. The chosen `σ₁₆ = 32` is `σ = 2.00`, i.e. 4.7% above it — an
independent route landing on the same ~3.5–5% margin as §6c's headroom. Weak (it leans on the 10-cell
fit the plan itself calls non-robust) but worth one line, because it is the only cross-check that σ and k
are individually sane rather than only jointly sufficient.

---

## 2. Logic review

### 2.1 Q1 — is the `f_required` derivation sound, and does it license one platform?

**The derivation: yes** (§1 above). **The one-platform licence: the stated justification is a non sequitur.**

§6f says: "`f`'s insensitivity to `static` (§6b) is what makes a one-platform derivation defensible in
the interim." Static-insensitivity says `f` does not depend on *which static value you assume*. It says
nothing about whether the *dynamic* term differs between platforms — and that is the whole question. The
platforms differ in both: Phase 2 measured σ = 1.29 Linux vs 1.41 macOS (`clump.rs:281-287`), a 9%
difference in the dynamic coefficient. Had that difference pointed the other way, Linux's `f` would be
~1.46 and the chosen 1.391 would be insufficient — and static-insensitivity would not have warned you.

**Two valid arguments exist and neither is made:**

- **Directional, from §12's own mechanism.** §12 attributes the excess to page-level fragmentation
  (blocks from many bins interleaved on one page). darwin/arm64 has 16 KiB pages against Linux/x86_64's
  4 KiB, so the same interleaving wastes ~4× more per page on darwin. The mechanism predicts darwin is
  the pessimistic platform.
- **Empirical, from the arc's own data.** Phase 2: σ 1.29 Linux / 1.41 macOS. §1: at 25 bp/c2 the
  prediction is exceeded on both platforms but the *budget* only on macOS. Both say darwin is worse.

Replace the static-insensitivity sentence with those. Static-insensitivity is still worth keeping — it
licenses not knowing `static` — just not this conclusion.

**Second-order, in the safe direction:** `static_actual` may itself rise with `--cores` (thread stacks,
mimalloc thread-local heaps). If it does, `f_required` at the binding c8 cell is *lower* than 1.337, so
using 205 (measured at c4) is conservative. Worth one clause so a future reader does not re-derive it.

### 2.2 Q2 — are σ₁₆ = 32, k₁₆ = 89 right?

**The level is right for the measured design points. Two things about how it is justified are not.**

**(a) The uniform-scaling justification is drawn from withdrawn reasoning.** §6c says both terms scale by
the same factor because "the breach rides `1/r`, which enters `σ` and `k` alike (§2.3)". But §6b
explicitly abandoned mechanism-based derivation, and §2.3's model is the one §12 rejected. Using
measurement to set the level and a withdrawn mechanism to set the split is inconsistent.

There is a sound reason and it is different: **the σ/k split is not identifiable from this design.** The
four cells contain only two distinct `n_bins : cores` shapes — 8:1 at c2 and 4:1 at c4/c8/c16 — so a
two-parameter solve is exactly determined with zero residual and no robustness. I ran it: solving the c2
and c8 constraints as equalities gives `f_σ ≈ 1.12`, `f_k ≈ 1.65` (≈ σ₁₆ 26, k₁₆ 106), which meets every
cell's bare requirement at a **21.0%** bin cut at c2 instead of 28.1%. It looks strictly cheaper — and it
is the wrong choice, because it spends the saving exactly where the data is noisiest: post-change headroom
at 25 bp c2 drops from 7.5% to 3.2%, against that cell's own 5.24% observed relative spread. **Uniform
scaling is right because it puts headroom where the spread is, not because of `1/r`.** Say that, and
record the non-uniform option as costed-and-rejected — otherwise the next reader will find the 7-point
saving and take it.

**(b) 3.1% is the wrong spread statistic, and the plan mixes two cells to get it.** §6b takes the largest
*absolute* within-cell spread (64.1 MiB, at 25 bp **c16**) and expresses it as a fraction of **that**
cell's dynamic part (3.13%), then applies it as a *relative* margin to a *different* cell's `f` (c8's
1.3370). Absolute MiB from one cell and a percentage from a second is not a quantity.

The defensible statistic is the largest *relative* within-cell spread. Recomputed, spread / (max RSS − static):

| cell | spread MiB | as % of dyn |
|---|---:|---:|
| 25 bp c2 | 47.8 | **5.24** |
| 25 bp c4 | 32.9 | 3.76 |
| 25 bp c8 (binding) | 14.3 | **0.65** |
| 25 bp c16 | 64.1 | 3.13 |
| worst of all 14 cells (100 bp c4) | 35.6 | **5.55** |

At `f = 1.391` the post-change headroom per cell, against that cell's own relative spread:

| cell | ratio new | headroom | own spread | margin |
|---|---:|---:|---:|---:|
| 25 bp c8 | 0.9655 | 3.45% | 0.65% | 5.3× |
| 25 bp c2 | 0.9251 | 7.49% | 5.24% | 1.43× |
| 25 bp c16 | 0.9011 | 9.89% | 3.13% | 3.16× |
| 25 bp c4 | 0.8967 | 10.33% | 3.76% | 2.75× |

So the choice is defensible **per cell**, and the plan's headline framing ("3.1% spread, 1.379 target")
is not how it is defensible. Two things make me uneasy about resting there:

- The binding constraint is the max of **3 reps of one cell**, and it is the cell with by far the
  *smallest* observed spread (0.65% against 3.1–5.2% elsewhere). An anomalously tight 3-draw sample
  setting the binding max is precisely where an under-sampled maximum bites. For n = 3 the expected
  range is ≈1.69σ; bounding a 99th percentile wants ≈2.33σ, i.e. ~1.4× the observed range.
- That cell runs in **57 s wall** (TSV `wall_s` 56.2/56.8/59.2). Ten reps is ~10 minutes. This is the
  cheapest risk reduction available anywhere in the plan.

**Recommendation:** run ≥10 reps of 25 bp / c8 / 2G before fixing the integers, restate the margin as
"max relative within-cell spread = 5.24% among 25 bp cells (5.55% across all)", and then either keep
(32, 89) with the per-cell table above as the explicit justification for accepting 3.45% at the binding
cell, or step to (33, 92) — `f` = 1.4359 at c8, headroom 6.1%, covering the 5.55% figure, at a 30.4%
bin cut. Note that the plan's own listed (33, 90) gives 5.4% headroom, which just misses 5.55%.

### 2.3 Q4 — does §6f's "does not block implementation" hold?

For the two gaps it lists, yes — both are honest and both are properly deferred to the PR. **§6f's
problem is what it omits.** Two axes carry more risk than either listed gap:

**(a) `f` is demonstrably not a constant across the design, and the design confounds cores with budget.**
Across the four 25 bp cells `f_required` runs 1.2069 → 1.3370, an 10.8% spread that is far larger than
the 0.65–5.24% within-cell noise, and it is **non-monotonic in cores** (c2 1.2574 > c4 1.2069 < c8
1.3370 > c16 1.2367). The two highest cells are the only two at 2 GiB. So `cores` and `budget` move
together throughout and neither effect is isolated. Taking the max over four points is conservative *at
those four points*; there is no cell at 1G/c8, 1G/c16, 2G/c2, or any budget above 2 GiB.

**(b) Linearity in `B` is only questioned downward, and the one upward data point contradicts it.**
Assumption 8 and Validation 6 ask only about shrinking `B` (more flush cycles). Nothing asks whether
`dyn ∝ B` holds *upward*, and users routinely pass `--memory 4G`/`8G`. Largest measured `B` is 22.80 MiB
(1G/c2); at 4G/c2 it is **112.89 MiB**, 5× outside the measured range.

The repo's own baseline file has the relevant data point and flags it as a puzzle. On Linux at 51 bp:
1G_c4 = 0.788 (B = 18.13 MiB, peak−static = 558 MiB) and 2G_c4 = 0.854 (B = 42.00 MiB, peak−static =
1416 MiB). B rises 2.317× while the dynamic part rises 2.538× — **~9% super-linear**, and the file itself
notes the 2G cell "reads *higher* than 1G_c4, not lower" despite non-saturation, which can only depress
it. If that 9% per budget doubling continues, the 3.45% headroom at the binding cell is gone by 4 GiB.

§10 Finding 1 already offered a second `--budgets` entry and recorded it as "offered, not taken" on
identifiability grounds. **The reason to take it has changed**: it is no longer about identifying
`static`, it is about whether the promise holds where users actually set `--memory`. One extra budget
(2G at c2/c4, or 4G at c4) on the 25 bp fixture is a few cells, and it converts the largest
extrapolation in the plan into a measurement. I would add it to §6f as a third gap and to §8 Step 9 as a
required cell; whether it blocks Step 4c is a judgement call, but it should not be invisible.

### 2.4 Q6 — internal contradictions left between v2 prose and what was measured

- **§9's entire pre-`---` block is stale v2 text about Option D**, which §6g withdrew. It lists under
  "Remaining risks": "(a) D's coefficient table is a first pass", "(b) D raises peak RSS at long reads",
  "(e) D touches a hot path called once per record" — risks for an option that no longer exists. Worse,
  under "**What is now verified rather than assumed**" it lists "the decomposition reproducing 23 and 64
  to two decimals; the zero-parameter prediction at 0.6% / 3.2%; ... and D's `text + spine = A` identity
  with its per-length table" — all of it computed from `s = 74`, which §6g withdrew and §11 excluded.
  Presenting withdrawn evidence under a "now verified" heading is the same failure mode the v3 addendum
  diagnoses two pages later ("a lower bound presented with certainty reads as a measurement"). Either
  delete the block or fence it as `### v2 — superseded`.
- **§2.3 is marked `[REVISED BY SPIKE]` but its body was not revised.** It still derives
  `σ(r) = 0.937 + 74/r → σ×16 = 23.00` and presents the zero-parameter confirmation as live. §12
  measured occupancy at 0.513, not 0.937, and Assumption 9 withdraws 0.937 explicitly. A heading marker
  is not a revision — as it stands a reader who starts at §2 is told 0.937 is measured and a reader who
  starts at §7 is told it is withdrawn. Same for §2.4, which states `a = 0.937` in its own argument.
- **§6a says "255–270 B/pair"; §11 measured 245.0–269.6** (mimalloc) / 246.0–270.7 (libmalloc). §5 row 3
  repeats "255–270". The low end is 245. This matters slightly because §6a's "only ~55%" is
  148/270; against 245 it is 60%.
- **§9's "Efficiency" claim is wrong on one point.** "The floor warning reuses the record
  `sanity_check_with` has already read, so no new I/O and no new failure path" — `fastq.rs:536` is
  `sanity_check_with(path, is_gzip) -> Result<()>`. It returns nothing. Reusing the record it read means
  changing that signature (or adding a sibling), and `main.rs:42` is a call site. That is not free, and
  the plan should say what it intends: change the return type to carry the first record, or add
  `sanity_check_with_record`. §8 Step 4c item 6 states the constraint correctly ("do not open the input
  again") but not that the plumbing does not exist yet.

### 2.5 §6d's cost claim is contradicted by the plan's own data and by the shipped docs

§2.5 says the cost is "paid for in bin size, i.e. in compression", and §6d states it as "28.1% of bin
size, for every input, at every read length" plus a docs sentence about 150 bp libraries losing ~38% of
budget. **Nothing measures the output-size effect, and the TSV has the column.**

At fixed read length, c2 vs c4 both have `n_bins = 16` and the same budget, so the pair isolates `B`:

| read length | `B` 22.80 → 18.13 MiB | `out_bytes` change |
|---:|---:|---:|
| 25 | −20.5% | **−0.228%** |
| 36 | −20.5% | −0.237% |
| 51 | −20.5% | −0.207% |
| 100 | −20.5% | −0.120% |
| 150 | −20.5% | −0.111% |

A 20.5% bin cut moves output by ≤0.24% — and *downward*. Varying `n_bins` instead (c2 16 bins vs c8 32
vs c16 64) moves it +1.3% / +1.7%, so `n_bins` is the driver and `n_bins` does not change under this
refit. The shipped docs already say the same thing in prose: *"In theory, bigger budget → bigger
per-gzip-member sort runs → better compression. In practice, increasing memory doesn't seem to make very
much difference in our tests."* (`docs/.../performance/clumpy.md`).

The honest caveat is that `mk_clump_fixture` produces all-distinct sequences (`distinct == pairs` in
every TSV row), so it is the wrong fixture for measuring a clumping *gain* — but that cuts both ways:
the plan is asserting a cost it has not measured, on a fixture that cannot measure it, in the direction
that makes the change look worse. `docs/perf_data/clumpify-rrbs-2026-05-07/` is the repo's real dataset
and there is precedent for running it. Either measure the compression cost there, or state the cost in
bin size only and stop calling it a compression cost. As written, §6d and the docs' existing "doesn't
make much difference" note will contradict each other in the same PR.

---

## 3. Assumptions

| # | plan's assumption | assessment |
|---|---|---|
| 1 | `prediction = budget × 10/11` invariant under `f` | **Verified algebraically and numerically.** Banner prints 930 of 1024 before and after at every configuration I checked (c2/c4/c8/c16, workers=1, fastqc=4). This matters more than the plan says — see §4.4. |
| 2 | `--clumpify` needs `--cores ≥ 2` | Holds; means the `cores = 1` floor row is unreachable in production. |
| 3 | integer arithmetic, `saturating_mul`, exact-inverse | Preserved by a coefficient-only edit. But see Validation 8 below: the "exact inverse ~312 B" property the plan promises to check **is not asserted anywhere in the tree**. |
| 4 | `2·size_of::<FastqRecord>() + size_of::<MinimizerKey>() == 148`, worth pinning | Agreed, and §11 measured it. **But no step in §8 pins it.** §2.3 says "Action:", §11 rec 1 says "worth pinning on its own terms", Assumption 4 says the same — and §8's twelve numbered items contain no such task, and §5 has no row for it. Dropped between analysis and outline. |
| 5 | `r` measured at the record count used | Confirmed; TSV `r_measured == r_analytic` to 3 dp in all 42 rows. |
| 6 | `f` insensitive to `static` | **Verified exactly.** Note the plan calls this Assumption 6 in §7 but the §9 v2 block calls the platform question "Assumption 6" — the numbering drifted between blocks. |
| 7 | platform independence — open, take the worse | Correct resolution; justification needs fixing (§2.1). |
| 8 | linearity in `B` — open | Correct to flag, but only half-flagged (§2.3(b)). |
| 9 | withdrawn: `s = 74`, `a = 0.937` | Correct — but §2.3/§2.4 still use both. |

**Unstated assumption worth surfacing: single-end behaves like paired.** All 42 runs are `--paired`. The
model's σ and k are shared with the SE path (`SingleBin::push`, `parallel.rs:1286`). The mechanics do look
equivalent — SE spine is 76 B/slot with `bytes_this = r`, so 76·B/r against paired's 74·B/r, and the
allocation rate per byte of text is identical (3 allocations per `r` bytes vs 6 per `2r`) — so I expect
it transfers. But SE at 25 bp *is* the `--small_rna`/miRNA case that §6f names as the sub-floor risk, and
it is the one combination the plan neither measured nor listed as an assumption. One SE 25 bp cell in
Validation 5 closes it.

**Unstated assumption: bins hold untrimmed records.** I checked — `parallel.rs:844-845` pushes records
straight from the reader and the worker trims, so the model's `r` is the *input* record size and Step
4c item 6's input-sampled warning measures the right quantity. This is load-bearing for the warning
being correct at all and is nowhere written down; if a future change moves trimming ahead of binning,
the warning silently starts measuring the wrong thing. Worth one line in §7.

---

## 4. Validation sufficiency (Q3 and Q5)

§5 rows 5–11 are the right *shape* — worst-rep, per-run gating, take-the-worse-platform. What follows are
the ways the change can produce a wrong or silently-degraded result that no row catches.

### 4.1 Critical — the CI promise gate is never armed on the 25 bp cells

`ci.yml`'s `measure` takes `$2` as a gating flag and, when it is 1, asserts `peak_kb -le pred_kb`.
Current flags:

```
measure 1G_c4        1   ...    # gated
measure 1G_c4_fastqc 1   ...    # gated
measure 25bp_1G_c4   0   ...    # NOT gated
measure 25bp_1G_c2   0   ...    # NOT gated
```

and `.github/clumpify-memory-baseline.tsv` says why, twice:

> "25 bp is above the prediction today: the per-bin coefficient is fitted per byte, not per record
> (#457). Recorded, not gated on the promise, **until that refit lands**."
> "25bp_1G_c2 is #457: over the printed prediction, inside the budget, on Linux. Its value will move when
> sigma is refitted, needing a re-baseline, **which is why it is not banded while that work is pending**."

§6e and §8 Step 8 address only the *drift band* centres. **Nothing in the plan flips those two gating
flags to 1**, and nothing updates the two baseline comments that will become false. The consequence is
concrete: the fix ships, the 25 bp promise breach is repaired, and CI still does not check it — so the
next change to touch the clump path can reintroduce #457 and go green. Add a Validation row and a Step:
arm `25bp_1G_c4` and `25bp_1G_c2` on the promise, band `25bp_1G_c2` (the baseline file already calls it
"eligible next pass"), and rewrite both comment blocks.

### 4.2 Critical — three unit tests fail, one pin is uncounted, and three comments go stale

§8 Step 5 says "the six-cell pin table and the absolute floor pin in `clump.rs` all move". That
undercounts. Recomputed against the code (all 25 tests pass today; values below are what the new
constants produce):

| test | today | after (32, 89) | outcome |
|---|---|---|---|
| `resolve_layout_default` | 112.89 MiB, window `[112, 114)` | **81.15 MiB** | **FAILS** — not mentioned |
| `resolve_layout_scales_with_cores` | 44.87 MiB, window `[44, 46)` | **32.25 MiB** | **FAILS** — not mentioned |
| `resolve_layout_fastqc_reserves_more_and_shrinks_bins` | hardcodes `/624`, delta 2 581 110 vs expected 2 581 110 | expected stays 2 581 110, actual becomes **1 855 544**, `abs_diff` 725 566 against a 1024 tolerance | **FAILS** — not mentioned; `624` is the *old* denominator written into the test body at `src/clump.rs:670` |
| `resolve_layout_pins_the_default_cell` | `19_006_356` | `13_663_555` | fails; a **seventh** pin, separate from the six-cell table, and its comment names `/624` |
| `sizing_constants_are_pinned_across_the_role_space` | 6 pins | `13_663_555`, `19_733_721`, `6_241_377`, `6_241_377`, `17_188_357`, `8_455_379` | covered by Step 5 |
| `clumpify_floor_at_the_default_cell_is_pinned_absolutely` | 281 | **294** | covered by Step 5 |
| `channel_depths_sum_to_the_per_worker_coefficient` | passes | **still passes, invariant now false** | see §4.3 |
| `floor_and_sizing_round_consistently_without_overflow` | passes | passes | coefficient-independent by construction (its own comment says so) |
| `synchronous_paths_are_not_charged_for_absent_workers` | passes | passes | relational |

Also stale after the edit: the arithmetic written into comments at `src/clump.rs:557` (`denom = 23×16 +
64×2 = 496`), `:566` (`= 1248`), `:576` (`/ 624`), and the `resolve_layout` docstring block at `:270-287`
which spells out σ = 23/16, k = 4, *and* the Phase 2 provenance sentence ("σ and k are fitted from a
10 M-pair paired-end matrix ... Measured values are σ = 1.29 (Linux) / 1.41 (macOS)"). That provenance
paragraph becomes wrong, not merely stale — the new constants are not fitted, they are a measured
conservative bound, and if only darwin runs they are not "both platforms" either.

None of these are silent — `cargo test` catches them. The cost of not listing them is that the
implementer discovers three red tests mid-change and has to decide on the spot whether each window is a
coefficient pin to update or a real regression. That is exactly the situation Step 5's warning about
Phase 3's five-of-six wrong pins exists to avoid.

### 4.3 Critical — a green test asserting a false invariant

`channel_depths_sum_to_the_per_worker_coefficient` asserts `work + 1 + result_per_core == 4` with the
message *"channel depths must sum to k = 4 (64/16 in dyn_denominator)"*, and `channel_depths`' docstring
says *"`k` above charges for `work + 1 in hand + result_per_core`, so their sum is the constant and
neither may move alone."* After the change k = 89/16 = 5.5625 and the sum is still 4. The test stays
green while the invariant it encodes is false.

This matters because that test is the guard against someone raising a channel depth without re-fitting —
the baseline file records a work-channel-depth revert as a measured +0.126 to +0.151 regression, so the
coupling is real. §8 Step 4c item 4 only mentions the `dyn_denominator` docstring.

The change needs a phrasing that keeps the coupling legible, e.g. k = 4 slots × 1.391 surcharge, with
the test asserting `(work + 1 + result_per_core) * 89 == 4 * 89` — or more usefully asserting the slot
count and the surcharge separately so a depth change still trips it. Whatever the form, this needs to be
a decision in the plan rather than a discovery in the diff, because the failure mode is silence.

### 4.4 Important — the "930 of 1024" banner survives, but by luck, and three tests assert it verbatim

`tests/integration_clump_only.rs:126`, `tests/integration_clump_only_ubam.rs:169` and
`tests/integration_preflight_tripwire.rs:605` all assert the literal string
`"predicted peak ≈ 930 of 1024 MiB budget"`. I verified it survives: predicted stays 976 128 919 B =
930 MiB at c2, c4, c8, c16, workers=1, and fastqc=4. That is a consequence of Assumption 1 (prediction is
`budget × 10/11` independent of the coefficients) and so is structural rather than lucky — but nothing in
§5 says "assert the printed prediction did not move", and if a future edit touches `MARGIN_NUM/DEN`
alongside the coefficients these three tests are the tripwire. Add the row; it costs a sentence.

### 4.5 Critical — the min-budget test is aimed at the wrong boundary and the failure mode is mis-stated

§8 Step 7 and Validation 9 say to check `--cores 2` and `16`, and describe the risk as "a `--memory`
value that worked before and now **bails**". Two corrections.

**It does not bail.** `main.rs:189-205` computes `clumpify_min_memory_bytes` and, below it, prints a
warning and `return Ok(None)` — falling back to plain mode. The output is valid; the user silently loses
the clumpify compression win, which is the whole feature. That is a quieter regression than an error and
should be described as such.

**Cores 2 and 16 are both far from the boundary.** Recomputed floors, MiB, current → new:

| `--cores` | no FastQC | with FastQC |
|---:|---|---|
| 2 | 281 → 294 (+4.8%) | 334 → 347 |
| 4 | 290 → 307 | 395 → 412 |
| 8 | 333 → 366 | 544 → 577 |
| 16 | 418 → **486** (+16.1%) | 841 → **908** |
| 32 | 590 → **724** (+22.8%) | 1012 → **1147** |
| 64 | 933 → **1202** (+28.8%) | 1356 → 1624 |

Inverted — the maximum `--cores` that keeps `--clumpify` engaged:

| budget | no FastQC | with FastQC |
|---|---|---|
| `--memory 512M` | 24 → **17** | 7 → **6** |
| `--memory 1G` (default) | 72 → **52** | 33 → **23** |
| `--memory 2G` | 167 → 120 | 128 → 92 |

**`--cores 24 --fastqc --memory 1G` clumps today and falls back to plain mode after this change.** On a
32-core workstation with defaults and FastQC on — an ordinary configuration — the feature switches
itself off. The boundary sits right where people are: `--cores 32 --fastqc --memory 1G` works today with
12 MiB to spare (1012 of 1024) and needs 1147 MiB afterwards.

**And for `--clump_only` the same shift is a hard error, not a fallback.** `docs/.../modes/clump-only.md:87`
records the deliberate asymmetry: *"Below the floor, `--clump_only` fails rather than degrading. This
differs from `--clumpify`, which warns and falls back to plain trimming"* — because reordering is the
entire job. So the plan's word "bails" is right for one mode and wrong for the other, and both need
testing. `--clump_only` floors (workers = 1), now → new:

| `--cores` | no FastQC | with FastQC |
|---:|---|---|
| 1–4 | 277 → 288 | 303 / 382 → 315 / 394 |
| 8 | 302 → 323 | 513 → 535 |
| 16 | 352 → **394** | 775 → 816 |
| 32 | 454 → **535** | 876 → 957 |
| 64 | 656 → **816** | 1078 → 1239 |

Max `--cores` before a hard error: at `--memory 512M`, **41 → 29**; at `--memory 1G`, 122 → 87 without
FastQC and **55 → 39** with it. So `--clump_only --cores 32 --memory 512M` succeeds today and errors
after — a working command starts failing, which is sharper than the `--clumpify` fallback and is the one
that will reach the issue tracker first.

Validation 9 should test the *boundary*, not two arbitrary core counts: for `--memory` ∈ {512M, 1G, 2G} ×
{FastQC, no FastQC}, find the largest `--cores` that still resolves, before and after. Then decide
deliberately: accept and document (the docs' "Below-floor behaviour" section already has the shape for
it), or clamp. Note `MIN_BIN_BYTES = 1 MiB` is what drives this and is itself unchanged — at 1G/c32 the
bin drops from 2.27 MiB to 1.63 MiB, so the margin over the floor narrows from 2.3× to 1.6×, which is
the mechanism behind the whole table.

### 4.6 Important — Validation 8 checks a property that is not asserted anywhere

Row 8 reads: "`predicted_peak_bytes()` is still the exact inverse of `resolve_layout` | existing unit
test across cores ≤ 32 | within the current ~312 B". I searched: **no test asserts a 312 B bound, and
`312` does not appear in `src/` or `tests/`.** The existing test
(`floor_and_sizing_round_consistently_without_overflow`) asserts `predicted <= budget`, one-sided, and its
own comment says this "holds for any σ, k or static, because the sizing and the prediction read the same
helpers". So row 8 as written cannot fail and cannot pass — there is nothing to run.

The 312 B figure is `max(denom)/16` at cores 32 with today's constants (4992/16 = 312.0), i.e. the
residual bound from the two integer divisions, not an assertion. **That bound grows with the
denominator**: 31.0 → 43.1 B at c2, 39.0 → 54.2 at c4, 78.0 → 108.5 at c8, 156.0 → 217.0 at c16,
**312.0 → 434.0 at c32**. Either restate row 8 as "the residual bound rises from 312 to 434 B at cores
32; `predicted <= budget` still holds" and leave it, or — better — add the two-sided assertion the row
implies (`budget - predicted <= denom/16 + 1`), which is the property being claimed and is currently
untested.

### 4.7 Important — Validation 6's tolerance and its attribution

Row 6 compares re-measured ratios against §6c's predicted column "within the 19 MiB residual RMS of
§10". Two issues. The 19.0 MiB is a *fit residual* on the pre-change data being reused as a tolerance for
a post-change ratio comparison — a category mix, and it is 2.1% of the pre-change dynamic part but 2.9%
of the post-change one, so it silently loosens. And the row attributes any systematic miss to
fragmentation ("~39% more flush cycles changed fragmentation"). There is a **second** candidate the row
does not name: a 28% smaller pool reaches steady state on fewer pairs, so a fixture that was marginally
saturating becomes more saturating and the ratio rises for that reason alone. The CI fixtures are the
exposed ones — `ci.yml`'s own comment says "10M pairs is where the bin pool reaches steady state; below
~6M every ratio is depressed". Both mechanisms push the measured ratio *above* prediction, so if the row
fires, distinguishing them changes the fix. Name both.

### 4.8 Important — no validation of the floor warning at all

Step 4c item 6 adds a warning below `r = 96 B` and §5 has no row for it: nothing checks it fires below,
nothing checks it stays silent above, and nothing checks it does not appear in the stderr of the three
integration tests that pattern-match on banner text.

The threshold also deserves a second look. `r = len(id) + len(seq) + len(qual) + 5` with the id including
`@` (`clump.rs:355-358`). The fixture's ~41.6 B header is what puts 25 bp at 96.6 B. Real data with a
terse header goes well below at ordinary read lengths:

| header | 25 bp | 36 bp | 50 bp |
|---|---:|---:|---:|
| `@SIM:1:FCX:1:15:6329:1234 1:N:0:ATCACG` (fixture, ~41 B) | 96.6 | 118.6 | 146.6 |
| `@SRR1234567.1 1 length=25` (fastq-dump, ~25 B) | 80 | 102 | 130 |
| `@SRR1234567.1` (`-F`, 13 B) | 68 | 90 | 118 |

So the warning fires on ordinary SRA-derived 25–36 bp data, and — the part that matters — the promise
does not hold there. Extrapolating the measured `f(r)` at c2 (`f ≈ 0.4997 + 73.19/r`, which reproduces
the five measured c2 cells to within 3.5%) gives `f_required ≈ 1.56` at r = 68 against the chosen 1.391:
a **~12% breach**, warned about but not prevented. §6f calls this "handled by the warning in Step 4c",
which is true only in the sense that the user is told. Given that the docs and CI both make the promise
absolute, consider whether a warning is the right instrument, or whether the denominator should carry a
`max(1, 96/r)` multiplier so short-`r` inputs are charged proportionally. That is not the per-record
exactness that §6a killed twice — it is the same conservative bound, extended below its calibration
floor by the one factor the mechanism agrees on. Worth a paragraph either way; silently promising
something that fails by 12% on SRA-style short reads is the kind of thing a user reports as a bug.

Also: `sanity_check` reads only R1's first record. For paired input with asymmetric mates, or for a
variable-length library whose first record is atypical, one record is a poor estimator — the same
one-record fragility Step 6 exists to fix on the reservation side. Sample a handful.

### 4.9 Optional — `--clump_only` and Perl byte-identity

Row 10 is right that `--clump_only` moves. I checked the surrounding CI: the `--clump_only` byte-identity
steps compare `paste - - - - | sort` multisets and the determinism steps compare two runs of the *same*
binary, so all of them survive a bin-size change. No committed golden clumped bytes exist. Row 11's
"no validated flag path runs `--clumpify`" is also correct — the `validation` matrix covers SE, PE,
hardtrim5, clock, demux. Both rows are fine as written; noting the multiset-vs-bytes distinction in row
10 would make it clear *why* they are safe rather than leaving it to the reader.

---

## 5. Efficiency

No concerns. The change is one integer pair in a function called once per run; §9's claim of no hot-path
cost is correct and was D's cost, not C's. The one efficiency claim that is wrong is the
"no new I/O and no new failure path" for the floor warning (§2.4 above) — it is still cheap, but it needs
a signature change, not a free read.

Step 6 / item 12 (cap the spine extrapolation) is correctly scoped as independent, correctly flagged as
invisible to every Phase 4 fixture, and correctly given its own variable-length test requirement. §10
Finding 2's observation that ordering does not matter is right.

---

## 6. Alternatives worth recording

1. **Non-uniform σ/k scaling** — costed in §2.2(a). `f_σ ≈ 1.12, f_k ≈ 1.65` meets every measured cell at
   a 21% rather than 28% bin cut, and should be recorded as rejected-because-unidentifiable rather than
   not considered. Two shapes in the design, zero residual, and it spends the saving at the noisiest cell.
2. **A short-`r` multiplier instead of a warning** — §4.8. Keeps C's logic, extends the bound below the
   calibration floor by the one factor `1/r` that survived both spikes.
3. **(33, 92) instead of (32, 89)** — §2.2(b). 6.1% headroom at the binding cell, covering the largest
   observed relative spread, at 30.4% bin cut instead of 28.1%. Given §4.5 shows the real user-visible
   cost is the fallback boundary rather than compression (§2.5), 2.3 more points of bin size may be
   cheaper than it looks — but it moves `--cores 24 --fastqc --memory 1G` further out of reach, so the
   two costs trade against each other and should be decided together.
4. **One extra budget cell** — §2.3(b). Already offered in §10 Finding 1 and declined for a reason that
   no longer applies.

---

## 7. Action items

### Critical

1. **Arm the CI promise gate on the 25 bp cells.** Flip `25bp_1G_c4` and `25bp_1G_c2` to `gating=1` in
   `.github/workflows/ci.yml`, band `25bp_1G_c2`, and rewrite the two baseline-file comments that say
   they are ungated "until that refit lands". Without this the deliverable of #457 ships unverified.
   Add as a §5 row and a §8 step.
2. **List every test and comment the coefficient edit breaks** (§4.2): `resolve_layout_default`,
   `resolve_layout_scales_with_cores`, `resolve_layout_fastqc_reserves_more_and_shrinks_bins` (hardcoded
   `624` at `src/clump.rs:670`), `resolve_layout_pins_the_default_cell` (a seventh pin, not in the
   six-cell table), plus the comment arithmetic at `:557`, `:566`, `:576` and the fitted-provenance
   paragraph in the `resolve_layout` docstring at `:270-287`. New values are computed in §4.2 above.
3. **Decide what happens to `channel_depths`' k = 4 invariant** (§4.3). It stays green while becoming
   false. Needs a phrasing and a test form in the plan, not in the diff.
4. **Re-aim Validation 9 at the boundary, and treat the two modes separately** (§4.5). `--clumpify`
   warns and falls back to plain mode (it does not bail); `--clump_only` errors, deliberately. Max
   `--cores` at the default `--memory 1G` drops 72 → 52 without FastQC and **33 → 23** with it, so
   `--cores 24 --fastqc --memory 1G` stops clumping; for `--clump_only`, `--cores 32 --memory 512M` goes
   from working to erroring (41 → 29 max cores). Test {512M, 1G, 2G} × {FastQC, no} × {`--clumpify`,
   `--clump_only`}, then accept-and-document or clamp.
5. **Recompute §6c's `B new` column against the binary** (§4.2 / §2). Every row scales the
   banner-*rounded* MiB, which is the exact practice Step 5 forbids. Correct banner values are 16 / 13 /
   15 / 7 MiB where the plan says 15.8 / 12.9 / 14.4 / 7.2 — the 25 bp c8 row is 15 MiB, not 14.4, because
   its exact `B` today is 21.00 MiB and the banner floors it to 20.

### Important

6. **Replace §6f's one-platform justification** (§2.1). Static-insensitivity does not address a
   platform-dependent *dynamic* term. Use §12's page-size mechanism (darwin 16 KiB vs Linux 4 KiB) and
   Phase 2's σ 1.29/1.41, both of which say darwin is the pessimistic platform.
7. **Restate the margin with the right statistic and run more reps at the binding cell** (§2.2b). 3.1%
   mixes two cells; the largest relative within-cell spread is 5.24% among 25 bp cells, 5.55% overall. The
   binding cell has only 3 reps and the *smallest* spread of any cell, and it runs in 57 s — take 10.
   Then justify (32, 89) with the per-cell headroom table, or step to (33, 92).
8. **Add two gaps to §6f** (§2.3): (a) `f` varies 10.8% across the four 25 bp cells, non-monotonically in
   cores, with budget confounded — no 1G/c8, 1G/c16 or 2G/c2 cell exists; (b) linearity in `B` is
   unverified *upward*, largest measured `B` is 22.80 MiB against 112.89 at 4G/c2, and the baseline
   file's own 1G_c4 = 0.788 vs 2G_c4 = 0.854 pair is ~9% super-linear on Linux at 51 bp. Take §10
   Finding 1's offered second budget.
9. **Fix the σ/k split justification** (§2.2a): drop the withdrawn `1/r` argument, state
   non-identifiability (two shapes), and record the non-uniform solve as costed-and-rejected.
10. **Either measure the compression cost on real data or stop calling it one** (§2.5). The TSV shows
    ≤0.24% output change for a 20.5% bin cut, in the *downward* direction, and the shipped docs already
    say memory barely affects compression. `docs/perf_data/clumpify-rrbs-2026-05-07/` is the right fixture.
11. **Add a §8 step and §5 row for the 148 B `size_of` pin** (§3, Assumption 4). Three sections call for
    it; no step schedules it.
12. **Fix Validation 8** (§4.6). The ~312 B exact-inverse bound is asserted nowhere; the existing test is
    one-sided `predicted <= budget` and is coefficient-independent. The residual bound rises 312 → 434 B
    at cores 32. Either restate the row or add the two-sided assertion it implies.
13. **Add validation rows for the floor warning** (§4.8) — fires below 96 B, silent above, absent from the
    three integration tests that match banner text — and reconsider the threshold: a terse SRA header puts
    ordinary 25–36 bp data at r = 68–102 B, where the extrapolated requirement is ~1.56 against the chosen
    1.391, a ~12% breach that is warned about but not prevented. Also sample more than one record, and note
    that `sanity_check_with` currently returns `()`, so "reuse the record it read" needs a signature change.
14. **Add one single-end 25 bp cell to Validation 5** (§3). All 42 runs are paired; SE 25 bp is the
    `--small_rna` case the plan names as the sub-floor risk.
15. **Fence or delete §9's pre-`---` block** (§2.4). It lists Option D risks and cites `s = 74`-derived
    results under "now verified rather than assumed", both withdrawn by §6g.
16. **Actually revise §2.3 and §2.4, not just their headings** (§2.4). Both still assert `a = 0.937` and
    the `σ(r) = 0.937 + 74/r` decomposition that Assumption 9 withdraws.
17. **Add a §5 row asserting the printed prediction is unchanged** (§4.4). It is — I verified 930 of 1024
    survives everywhere — and three integration tests depend on it verbatim.
18. **Name the second mechanism in Validation 6** (§4.7): a smaller pool saturates on fewer pairs, so a
    marginally-saturating fixture reads higher. Also reconsider reusing a pre-change fit residual as the
    post-change tolerance.

### Optional

19. Add the two missing cells (100 bp c4 = 0.8851, c16 = 0.9165) to §6b's table.
20. Correct "255–270 B/pair" to 245–270 in §6a and §5 row 3.
21. State §6b's insensitivity as the derivative `(peak − prediction)/(prediction − static)²` so a reader
    can see when it would stop holding.
22. Add §12's independent corroboration of `σ₁₆ = 32` (implied pipeline `σ(96.6) ≈ 1.91` vs charged 2.00).
23. Record in §7 that bins hold *untrimmed* records (`parallel.rs:844-845`), which is what makes the
    input-sampled `r` warning measure the right quantity.
24. Note in Validation 10 *why* the `--clump_only` CI steps survive (multiset and same-binary comparisons,
    no committed golden bytes).
25. `modes/clump-only.md:85` needs its floors moved too (277 → 288 at `--cores` 1–4, 454 → **535** at
    `--cores 32`), and `:87`'s fail-rather-than-degrade paragraph is why Validation 9 must cover both modes.
26. Docs inventory for Step 5: `performance/clumpy.md` display maths at `:100-103` carries `23 × n_bins +
    64 × workers` (and the page is the KaTeX positive control in CI), `:104` the 32 bins × 9/44 MiB
    examples (→ 6 / 32 MiB), `:109`'s "reduced bin sizes by 24–31%" caution, `:114`'s full floor table
    (281/590/334/1012/933/1356 → 294/724/347/1147/1202/1624), `:119`'s literal `need ≥ 311 MiB` example
    at `--cores 6` (→ 336), and the "fitted ... on both macOS/arm64 and Linux/x86_64" sentence at `:92`,
    which is only true if the Linux cells run.
