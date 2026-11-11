# Plan Review — Phase 4 v3 (#457), Reviewer B

**Target:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v3
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `perf/clumpify-spine-accounting`
**Reviewer:** B (independent; no contact with Reviewer A)
**Date:** 2026-08-19

Every number in §6 was recomputed from
`/private/tmp/.../scratchpad/p4matrix/results_phase4_darwin.tsv` without reference to the plan's own
arithmetic. §10's fit was re-derived from scratch until it reproduced. The two spikes' method arguments
were checked against their code. `src/clump.rs`, `src/parallel.rs`, `src/clump_only.rs`,
`.github/clumpify-memory-baseline.tsv` and `.github/workflows/ci.yml` were read. No source file was
modified and nothing in §8 was implemented.

**Verdict: the decision is sound and the integers are right. The plan around them is not ready to
implement.** Option C's derivation reproduces exactly, `(32, 89)` is the minimal integer pair meeting the
stated target, and the choice is insulated from every instrument defect I found. But §6e states the exact
opposite of what the CI band arithmetic gives, §8 Step 4c under-scopes the change by five failing tests,
and §5's acceptance tolerance is set below the measurement noise — the same error §4 Fault 2 exists to
name.

---

## 1. What reproduced, and what did not

### 1.1 Reproduced exactly

| §6 claim | independent recompute | verdict |
|---|---|---|
| `f_required = (peak−static)/(prediction−static)`, worst rep per cell, static 205 | all **12** listed cells match to 4 dp | ✓ |
| binding `f = 1.3370` at 25 bp, c8, 2G, **rep 2** | 2419.1 → (2419.1−205)/1656 = 1.33702 | ✓ |
| static-insensitivity table 1.3306 / 1.3370 / 1.3409 / 1.3486 | 1.33063 / 1.33702 / 1.34093 / 1.34860 | ✓ |
| binding cell is the same at every static | stable at (25, c8, 2G, rep 2) for static 173→450 | ✓ |
| largest 25 bp within-cell spread 64.1 MiB, 3.1% of that cell's dynamic part | 2253.0−2188.9 = 64.1; 64.1/2048.0 = 3.130% | ✓ |
| `f_σ`, `f_k` for all four candidate pairs | 1.348/1.344, 1.391/1.391, 1.435/1.406, 1.478/1.469 | ✓ |
| `B` cut 25.7 / **28.1** / 29.6 / 32.1 % | 25.71 / 28.11 / 29.73 / 32.17 % | ✓ |
| worst post-change ratio 0.994 / **0.965** / 0.946 / 0.917 | 0.9940 / 0.9655 / 0.9462 / 0.9171, all at (25, c8, 2G, rep 2) | ✓ |
| §6c's twelve `ratio new` values | all match to 3 dp | ✓ |
| `B now` column (22 / 18 / 20 / 10 MiB) | matches the banner and the integer formula | ✓ |
| §10's `a = 0.786`, `b = 2.777`, `s = 108.74`, RMS 19.0, corr +0.022 | reproduced **only** with `B` taken as the truncated banner integer, cell maxima, static 205 — see §2.7 | ✓ (with a caveat) |
| §10's ratio table | all fourteen are the cell **maxima** | ✓ |
| §12's σ ratios 1.75 / 1.44 / 1.22 | 1.748 / 1.437 / 1.224 against `σ(r) = 0.786 + 108.74/r` | ✓ |
| §2.3's zero-parameter 118.0 MiB | 74×24×19,006,356×(1/95.9−1/147.9) = 118.0 MiB | ✓ (but see §2.11) |
| `r` model | `analytic_r(10M, 25) = 95.889`; measured 96.642 at 31.04 M pairs | ✓ |
| `size_of::<FastqRecord>() = 72` | three `String`s, no other fields (`fastq.rs:58-65`) | ✓ |
| `estimated_record_bytes = id+seq+qual+5` | `clump.rs:355-358`; the harness's awk matches | ✓ |
| `min_required` rises on its **dynamic part** only | c2 280.5→293.8 MiB, c16 418.0→485.1, c32 589.6→723.8 | ✓ — §8/§9's phrasing is correct, §5 row 9's "grows 1.39×" is loose |
| `B × dyn_denominator` stays inside u64 at `u64::MAX` budget | peaks at u64::MAX − 15, both old and new constants | ✓ no overflow introduced |

**`(32, 89)` is the minimal pair meeting `f ≥ 1.379`**: `σ₁₆ ≥ 31.72`, `k₁₆ ≥ 88.26`. The choice is not
just adequate, it is the cheapest adequate one. Good.

### 1.2 Did not reproduce

| §6 claim | recompute | consequence |
|---|---|---|
| `B new` = 15.8 / 12.9 / 14.4 / 7.2 MiB | **16.392 / 13.031 / 15.095 / 7.548 MiB** (banner 16/13/15/7) | §2.5 — wrong, and by the mechanism §8 Step 5 forbids |
| "Every measured cell after the change" (12 rows) | 14 cells measured; **100 bp c4 1G** and **100 bp c16 2G** are absent | §2.6 |
| "8 of 14 cells breach the promise, **on every rep**" | 51 bp c2 breaches on **1 of 3** reps (0.9878, 0.9953, 1.0024) | §2.11 |
| "within-cell spreads (5–64 MiB) far below the breach (150–250 MiB)" | breaches run 2.2 → 558.1 MiB; 51 bp c2's breach (2.2) is *below* its spread (13.6) | §2.11 |
| §6e "every banded cell moves by 0.2–0.3" | 1G_c4 0.788→0.619 (−0.169); 25bp_1G_c4 0.914→0.709 (−0.205) | direction right, high end overstated; conclusion unaffected |
| §6f "shortest supported `r` is 96.6 B" vs §8 Step 6's threshold "96 B" | 0.64 B of calibrated range left uncovered by the warning | §3.4 |

---

## 2. Logic review

### 2.1 CRITICAL — §6e's central claim about the CI band is the opposite of the truth

§6e: *"**The band's own regression-detection power is unaffected**, because it is set by width and by
lever A's magnitude, neither of which this change touches."* Step 8 follows it: *"width stays 0.09; both
cells stay armed."*

**Lever A's magnitude is exactly what this change touches.** Lever A is the work-channel-depth revert
(`channel_depths(true).work` 1 → 2). An extra queued batch is *one whole bin* per worker, so its cost is
`cores × σ_eff × B`, and its magnitude **in ratio-of-prediction units is proportional to `B`**.
`prediction = budget × 10/11` does not move; `B` falls by `f = 1.391`. Therefore lever A's delta falls by
1.391 too.

From `.github/clumpify-memory-baseline.tsv`'s own recorded figures:

| | half-range (27 obs) | lever A delta | margin at width 0.09 |
|---|---:|---:|---:|
| `1G_c4` today | 0.0275 | 0.126 | **0.036** |
| `1G_c4` after | 0.0198 | **0.0906** | **0.0006** |
| `25bp_1G_c4` today | 0.0215 | 0.128 | **0.038** |
| `25bp_1G_c4` after | 0.0155 | **0.0920** | **0.0020** |

The baseline file states the rule the band must satisfy: `width ≥ 3× observed half-range` **and**
`width ≤ measured regression delta − a working margin`. After the change, width 0.09 fails the second
condition on both cells: the band would be nominally armed while sitting within 0.002 of the only
regression it was ever calibrated to catch. One bad draw and lever A passes silently.

The noise side improves (0.09 becomes 4.5× / 5.8× the half-range), which is precisely the room needed:

> **Narrow the width to ≈ 0.065 in the same PR.** 3× half-range is 0.059 / 0.047 ≤ 0.065 ✓, and the
> margin against a 0.0906 delta becomes 0.026 — comparable to today's 0.036. Consider narrowing the
> −0.10 down-side to ≈ −0.072 for the same reason.

Two caveats, both worth writing down rather than assuming: (i) the post-change lever-A delta is inferred
from `B`-proportionality, not measured — the PR can measure it directly by re-running the revert once,
which is how 0.126–0.151 was obtained in the first place; (ii) the direction is not in doubt either way,
so §6e cannot stand as written even if the exact width lands elsewhere.

**This compounds with §2.2 below**: after the change nothing in the codebase ties `k` to
`channel_depths` any more either. Both of lever A's tripwires degrade in the same commit, and the plan
notices neither.

### 2.2 CRITICAL — Step 4c is not "one line of arithmetic": five tests fail and one silently rots

§8 Step 4c reads as a two-integer edit plus a docstring, and Step 5 accounts for *"the six-cell pin table
and the absolute floor pin"*. The real footprint in `src/clump.rs`:

| test | what it asserts | after `(32, 89)` |
|---|---|---|
| `resolve_layout_default` (`:574`) | `B ∈ [112, 114) MiB` at 4 GiB/c2; comment names `denom = 496` | **81.15 MiB — FAILS** |
| `resolve_layout_scales_with_cores` (`:565`) | `B ∈ [44, 46) MiB` at 4 GiB/c8; comment names `denom = 1248` | **32.25 MiB — FAILS** |
| `resolve_layout_pins_the_default_cell` (`:573`) | `B == 19_006_356`; comment names `/ 624` | **13,663,555 — FAILS** |
| `sizing_constants_are_pinned_across_the_role_space` (`:583`) | six exact budgets | **all six — FAILS** |
| `clumpify_floor_at_the_default_cell_is_pinned_absolutely` (`:632`) | `floor.div_ceil(MiB) == 281` | **294 — FAILS** |
| `channel_depths_sum_to_the_per_worker_coefficient` (`:641`) | `work + 1 + result == 4`, message *"must sum to k = 4 (64/16 in dyn_denominator)"* | **still passes, rationale now false** |

That is **five failing tests and ten hard-coded expectations**, not one. The last row is the dangerous
one: `k` becomes `89/16 = 5.5625`, the depths still sum to 4, and the assertion keeps passing while
guarding a relationship that no longer exists. The same falsehood is in `channel_depths`'s own docstring
(`:177-178`, *"their sum is the constant and neither may move alone"*) and in `resolve_layout`'s
docstring (`:273-286`: `σ = 23/16`, `k = 4`, the literal `23 × n_bins + 64 × workers`, and the provenance
sentence *"σ and k are fitted from a 10 M-pair paired-end matrix over cores ∈ {2,3,4} … the constants
take the worse of each"*, which will be false in every clause). `MARGIN_NUM`'s comment (`:263`, *"Largest
observed within-cell spread is 96 MiB"*) is also stale — the byte-equalised fixtures give 64.1 MiB.

New pin values, computed from the integer formula as a **cross-check only** — Step 5 is right that they
must be read off the binary:

| site | old | new |
|---|---:|---:|
| 1 GiB, `uniform(4)` (two sites) | 19,006,356 | 13,663,555 |
| 1 GiB, `bin_cores 4 / workers 1` | 27,453,626 | 19,733,721 |
| 1 GiB, `bin_cores 4 / workers 4 / fastqc 16` and `17` | 8,681,915 | 6,241,377 |
| 1 GiB, `uniform(2)` | 23,911,222 | 17,188,357 |
| 4 GiB, `uniform(32)` | 11,761,649 | 8,455,379 |
| absolute floor at c2 (MiB, `div_ceil`) | 281 | 294 |

### 2.3 IMPORTANT — the budget axis is confounded with cores, and the highest `f` sits at the highest budget measured

§6f lists four things "still not measured". The one it omits is the one that could invalidate the chosen
integers: **`f` was measured at 1 GiB for `cores ∈ {2,4}` and at 2 GiB for `cores ∈ {8,16}`, and nowhere
else.** No read length has two budgets at one core count. So the observed pattern

```
25 bp @ 1G:  c2 1.2574   c4 1.2069        (mean 1.206)
25 bp @ 2G:  c8 1.3370   c16 1.2367       (mean 1.271)
```

is equally consistent with "`f` peaks near c8" and with "`f` rises with pool size". §6b argues `f` is
budget-invariant *under the model* (every term in `peak − static` is ∝ `B`, so the ratio is not), and
that argument is correct — but it is the same proportionality that §6f/Assumption 8 concede is
unverified, and §12 gives a mechanism (fragmentation from cross-bin page interleaving) that has no
obligation to be scale-free. If half of the 0.065 gap is the budget axis, an 8 GiB run — routine on WGBS
hosts — would need `f ≈ 1.45` and `(32, 89)` would breach.

**The experiment is ~3 minutes.** One cell — **25 bp, `--cores 8`, `--memory 1G`** — shares its core
count with the existing 2 GiB cell and isolates the budget axis outright (`min_required(c8) = 332 MiB`,
so 1 GiB resolves; `B` = 9.06 MiB; measured wall at c8 was 56 s/rep). A second point, 25 bp/c4 at 512 MiB
(`B` = 6.19 MiB, a 2.9× lever arm), simultaneously **identifies `static`** — §10 Finding 1's declined
offer — and **verifies linearity in `B`**, which is Assumption 8 and Validation 6's whole premise.

Note the sequencing this exposes: the plan declines a second budget *before* the change (Finding 1,
"offered, not taken") and then mandates a full post-change re-measurement (Validation 6) that costs more.
Spending 15 minutes now de-risks the integers; spending two hours later only tells you whether you were
right.

### 2.4 IMPORTANT — Validation 6's tolerance is below the measurement noise, which is §4 Fault 2's own error

§5 row 6: compare re-measured ratios against §6c's predicted column, *"within the 19 MiB residual RMS of
§10"*.

19 MiB is the **fit's residual RMS**, not the measurement's dispersion. The within-cell spreads in the
same TSV run **3.9 to 64.1 MiB** (25 bp c16: 64.1; 100 bp c8: 56.6; 25 bp c2: 47.8), and they will still
be ~46 MiB after scaling by 1/1.391. A ±19 MiB acceptance band therefore fails on noise alone in most
cells, and a red result would carry no information.

§4 Fault 2 makes exactly this criticism of v1's falsification cell — *"~32 MiB of peak against
within-cell spreads of 96 MiB (macOS) and 37 MiB (Linux): 3× below the noise"* — and §5 row 6 reproduces
it. Set the tolerance from the post-change within-cell spread measured in the same runs (≥ ~50 MiB, or
state it as "no cell's worst rep above its §6c prediction by more than its own observed spread").

### 2.5 IMPORTANT — §6c's `B new` column violates §8 Step 5's own rule

The four `B new` values (15.8 / 12.9 / 14.4 / 7.2 MiB) are the **banner-truncated** old values divided by
`f`. Recomputed through `resolve_layout`'s integer arithmetic they are **16.392 / 13.031 / 15.095 /
7.548 MiB**, and the banner will print **16 / 13 / 15 / 7**. The errors are 0.09–0.72 MiB (up to 4.9%),
largest exactly where truncation was largest (c8 at 20.998 → 20, c16 at 10.499 → 10).

`§8 Step 5`: *"Recompute every pinned cell against the binary — do not scale the old ones … Phase 3's PR
2 recorded five of six pin values wrong when they were drafted from illustrative numbers."* §6c is the
table an implementer will reach for. The `B cut` (28.1%) and every `ratio new` are unaffected because
they use exact ratios — only the MiB column is wrong — but it is wrong in the precise way the plan warns
about, one section later.

### 2.6 IMPORTANT — "Every measured cell after the change" is 12 of 14

§6c omits **100 bp c4 1G** (0.910 → 0.716) and **100 bp c16 2G** (0.926 → 0.696); §6b's per-cell `f`
table omits the same two (0.8851 and 0.9165). Both are non-binding, so nothing in the decision changes —
but §10's ratio table has all fourteen, so §6 and §10 disagree on the size of the evidence base, and a
claim of completeness that is 86% complete is the kind of thing the next reviewer has to re-derive to
trust anything else. Add both rows.

### 2.7 IMPORTANT — `calibrate.py` regresses on the truncated banner MiB, not `bin_byte_budget`

`calibrate.py:465` (`design_row`): `b, n, c = r["b_banner"], r["bins"], r["cores"]`. `b_banner` is
`layout.bin_byte_budget / (1024*1024)` **as printed and integer-truncated by the banner** (`main.rs:217`).
The true `B` values are 22.804, 18.126, 20.998, 10.499 MiB, so the design matrix carries a **−0.7% to
−4.8% error that is correlated with the cores axis** (c4 barely truncated, c8/c16 ~4.8%).

Refitting on exact `bin_byte_budget`, everything else identical:

| | plan (banner `B`) | exact `B` | shift |
|---|---:|---:|---:|
| `a` | 0.786 | **0.676** | −14% |
| `b` | 2.777 | **3.234** | +16% |
| `s` | 108.74 | **106.48** | −2.1% |
| residual RMS | 19.0 MiB | 18.3 MiB | — |
| falsification errors | +8.4 / +0.9 / +16.0 / +9.3 % | **+4.6 / −2.6 / +11.6 / +5.1 %** | one flips sign |

Consequences, in order of importance:

1. **§6's chosen integers are unaffected.** `f` comes from raw `(peak − static)/(prediction − static)`
   and never touches the fit. Say so explicitly — it is the strongest thing in §6b's favour.
2. **§10's "all in the unsafe direction" is an artefact.** With exact `B`, 25 bp c16 is *over*-predicted
   by 2.6%. §5 row 4's and §10's framing of Verify 4 needs the correction.
3. **§12's "1.22–1.75×" becomes ≈1.37–1.88×**, because it divides the spike σ by `σ(r) = a + s/r` from
   the biased fit. The figure appears in §6f, §7 and §12. Nothing depends on it.
4. **Fix it before the Linux run.** Otherwise the two platforms' `a`/`b` are compared through a
   platform-dependent truncation error (different `B`, different rounding loss).

The irony is worth one line in the plan: a harness whose stated purpose is that *"a pattern that matches
but yields the wrong numbers is the failure Phase 2 could not distinguish from a healthy parse"* then
treats a deliberately lossy printed integer as the model's exact input. `parse_banner` already knows
`budget`, `bins` and `cores` — it can compute `B` exactly and assert the banner truncates to it.

### 2.8 IMPORTANT — no validation row for the thing the change actually costs users

§6d states the cost as "28.1% of bin size". The user-visible consequence is **output size**, and §5's
eleven rows do not mention it. §2.5 says the cost is "paid for in bin size, i.e. in compression" and then
nobody measures compression.

What the TSV does say, holding `n_bins` at 16 and cutting `B` from 22.804 → 18.126 MiB (−20.5%):

| read length | Δ output |
|---|---:|
| 25 bp | −0.23% |
| 36 bp | −0.24% |
| 51 bp | −0.21% |
| 100 bp | −0.12% |
| 150 bp | −0.11% |

Reassuring — except **`distinct == pairs` in every row of the TSV**: `mk_clump_fixture` emits 100%
distinct sequences, which is the exact opposite of the duplicate-rich case `--clumpify` exists for and
where its 20–35% savings come from. Smaller bins mean fewer near-duplicate neighbours per sorted run and
a smaller effective window, so this fixture cannot see the loss.

v2's §9 carried this as a live risk (*"(d) The compression cost needs real libraries, deferred once
already for the RRBS savings row"*). **v3's addendum drops it.** Restore it as a validation row, and make
it cheap: Step 9 is already re-running cells, so add one RRBS/WGBS fixture at `--memory 1G` (new
`B` = 16.39 MiB) against `--memory 1328M` (which restores today's 22.80 MiB exactly), and diff
`out_bytes`. That isolates `B`'s effect on savings with two runs.

Related docs consequence: `docs/…/performance/clumpy.md`'s caution box says the v2.4 refit *"reduced bin
sizes by 24–31% at budgets of 2 GiB and above. Figures at the default `--memory 1G` are unaffected."*
After this change the default 1 GiB **is** affected (−28%), so that exemption becomes wrong and every
savings row goes stale. Step 10 doesn't name it.

### 2.9 IMPORTANT — Step 6 cannot reuse what it says it reuses

§8 Step 6: *"**Reuse the record `sanity_check_with` already reads** — do not open the input again … v1's
Step 4.7 sampled afresh, which adds an I/O failure path."*

`FastqReader::sanity_check_with` (`fastq.rs:536`) returns `Result<()>`. It reads the record, inspects it,
and drops it. So does `main.rs::sanity_check_any` (`:37`). There is nothing to reuse without changing
both signatures — plus every caller and the `sanity_check_with_accepts_gzip_under_non_gz_extension`
test — and the uBAM branch returns a `noodles` record, not a `FastqRecord`, so `estimated_record_bytes`
is not directly available there.

The design instinct is right and the ordering works (`sanity_check_any` runs well before the layout at
`main.rs:184`). But the step needs to say what it changes: return `Option<usize>` (or the record) from
both functions, and decide what the uBAM path contributes. Also state that `r` will be sampled from
**R1 only** while the pool charges `bytes_this = r1 + r2`, which is exact only for equal-length mates.

### 2.10 IMPORTANT — the CI gate this refit exists to arm is never armed

`.github/workflows/ci.yml:444`:

```
# 25 bp is above the prediction today: the per-bin coefficient is
# fitted per byte, not per record (#457). Recorded, not gated on the
# promise, until that refit lands.
measure 25bp_1G_c4   0 --clumpify --paired --cores 4 --memory 1G "$F25R1" "$F25R2"
```

**This is that refit.** Nothing in §3, §5 or §8 flips `25bp_1G_c4` from `gating 0` to `gating 1` or
deletes that comment. The same applies to `25bp_1G_c2` (`gating 0`, the baseline file calls it *"#457:
over the printed prediction … not banded while that work is pending"* at 1.0301 → ~0.79 after the change)
and arguably to `1G_c2` (*"too tight to gate"* at 0.879 → ~0.63, no longer tight).

Arming those gates is the visible closure of #457 in CI, and it is the difference between "we believe the
promise holds at 25 bp" and "CI fails if it stops holding".

One clause for Step 8 while you are there: the src change and the baseline re-record must be in the **same
commit**, not merely the same PR, or the intervening commit's CI is red on the −0.10 down-side.

### 2.11 IMPORTANT — §9's first half is stale v2 text whose risks are all about a withdrawn option

§9's opening block is v2's Self-Review, and it survives unedited into a document where Option D is
withdrawn (§6g). Its "What is now verified rather than assumed" still lists *"the zero-parameter
prediction at 0.6% / 3.2%"* and *"D's `text + spine = A` identity with its per-length table"* as
achievements. Its "Remaining risks" (a)–(e) are **five items about Option D**: *"D's coefficient table is
a first pass"*, *"D raises peak RSS at long reads"*, *"D touches a hot path called once per record"*. An
implementer reading §9 top-to-bottom meets five live risks for code that will never be written, and then
a second "Remaining risks" list four paragraphs later that contradicts them.

Same class of problem in §2, at a lower grade:

- **§2.3** is marked `[REVISED BY SPIKE]` and then presents, unqualified, *"Two constants, one equation
  each, both landing on the shipped integers to two decimals"* and a **"Zero-parameter confirmation"**
  whose zero parameters are `s = 74` and `a = 0.937` — both withdrawn in §7 Assumption 9. The arithmetic
  is right (I reproduce 118.0 MiB); what it confirms is that `74 × (n+2c) × B × Δ(1/r)` happens to land
  near a measured difference, which after §11 is a coincidence needing an explanation, not a
  confirmation. One sentence saying so would do.
- **§10's evidence sentence overstates**: *"8 of 14 cells breach the promise, on every rep, with
  within-cell spreads (5–64 MiB) far below the breach (150–250 MiB)"*. 51 bp c2 breaches on 1 of 3 reps
  (0.9878 / 0.9953 / 1.0024) and its breach (2.2 MiB) is *below* its spread (13.6 MiB). Actual breaches
  run 2.2 → 558.1 MiB. The conclusion is safe — the binding cell breaches by 558 MiB against a 14.3 MiB
  spread, which is a 39:1 signal — so say that instead.
- **§5 row 3's "total 255–270 B/pair"** is a hybrid: non-text resident per pair from §11's table is
  254.9 / 275.2 / 279.6 / 276.4 / 270.1, i.e. **255–280**; `2s` is 245–270. Pick one.

### 2.12 IMPORTANT — the one action §2.3 names is in no implementation step

§2.3: *"**Action:** pin `size_of::<FastqRecord>()` beside `estimated_record_bytes`, so the model's
constant and the struct cannot drift apart silently."* §7 Assumption 4 keeps 148 as an invariant "worth
pinning as an invariant on its own terms". §11 recommendation 1 says the same. **No step in §8 does it.**

It is a two-line `const _: () = assert!(...)` or a unit test, and it is the only thing standing between a
future `FastqRecord` field addition and a silent 24-byte-per-record under-charge. Per CLAUDE.md's
"put an invariant on the thing it constrains", it belongs above `estimated_record_bytes` in `clump.rs`,
not in the plan.

### 2.13 IMPORTANT — Step 12 names two of four `reserve_exact` sites

§8 Step 12 (cap the spine extrapolation) cites `parallel.rs:901-911`, and §2.3 adds `SingleBin::push` at
`:1286`. But `src/clump_only.rs` carries its **own duplicated** `SingleBin` (`:89`, `push` at `:104`) and
`PairedBin` (`:131`, `push` at `:147`), each with the identical `reserve_exact`-once-per-cycle pattern —
the comment at `clump_only.rs:85` even says *"Mirrors `parallel::PairedBin`'s shape"*. Fixing two of four
leaves `--clump_only` over-reserving on variable-length input, and the divergence is invisible because
the two files never reference each other's constants. `clump_only.rs` also has four `resolve_layout` call
sites, so Validation 10's "banner and layout change consistently" has four banners to check.

### 2.14 OPTIONAL — the `f ≥ 1.379` step is the one number in §6 not derivable from a stated rule

The margin construction is: take the **largest absolute** within-cell spread among the 25 bp cells
(64.1 MiB, from **c16**), express it as a fraction of **that** cell's dynamic part (3.13%), then apply
that fraction to a **different** cell's `f` (c8's 1.3370) → 1.379. Three reasonable rules give three
answers:

| rule | target `f` | `(32, 89)` at 1.391 |
|---|---:|---|
| each cell margined by its own relative spread, take the max | 1.346 | passes comfortably |
| the plan's hybrid | 1.379 | passes |
| largest **relative** spread among 25 bp cells (5.24%, at c2) applied to the binding `f` | **1.407** | **fails** → `(33, 91)` |

`(32, 89)` satisfies two of three, so I am not disputing the choice — but the plan states 1.379 as though
it followed from something, and the reader cannot tell which rule produced it or that a neighbouring rule
rejects the answer. Two things fix it: name the rule, and add the observation that makes the residual risk
small — **the promise the user actually holds is `peak ≤ --memory`, and at the binding cell the
post-change peak is 1796.7 MiB against a 2048 MiB budget, i.e. 12.3% headroom.** A 5% run-to-run excess
would breach the *printed prediction* (0.965 → ~1.008) and stay far inside the *budget*. §6 never
separates the two promises even though §1 names both.

### 2.15 OPTIONAL — §5 row 8 asks for a check that does not exist

Row 8: *"`predicted_peak_bytes()` is still the exact inverse of `resolve_layout` | existing unit test
across cores ≤ 32 | within the current ~312 B"*.

No test asserts 312 B. `floor_and_sizing_round_consistently_without_overflow` (`:679`) asserts
`predicted ≤ budget` and its comment says explicitly *"Not a coefficient guard — `predicted ≤ budget`
holds for any σ, k or static"*. The 312 B is `max(denom)/16` at cores ≤ 32 = 4992/16; under `(32, 89)`
it becomes **6944/16 = 434 B** (and 868 KiB at the test's cores = 65536). So row 8 will pass whatever you
do, and its stated tolerance is wrong. Either restate it as `denom/16` (312 → 434 B) or point at the
assertion that exists.

Verified separately: `bin_byte_budget * dyn_denominator` peaks at `u64::MAX − 15` under both constant
sets — structurally safe, no overflow introduced, because `B = floor(sat/denom)`.

---

## 3. Assumptions

### 3.1 Assumption 6 (`f` insensitive to `static`) — sound, and **under-claimed**

The mechanism: `df/dS = (peak − prediction)/(prediction − S)²`. At the binding cell that is
558.1/1656² = 2.0×10⁻⁴ per MiB, hence 0.018 over an 87 MiB range. Verified.

The plan shows a 173–260 MiB table. **The actual statement available is far stronger: `f = 1.391`
tolerates any true `static` up to ≈ 434 MiB** (solving `(2419.1 − S)/(1861 − S) = 1.391`), and the
binding cell does not change anywhere in 173–450. That is worth saying, because it converts "we sampled
a plausible range" into "the answer is right unless `static` is more than double its measurement".

**What would break it**, in the order I'd rank them:

1. **`static` not being additive and `B`-independent.** `STATIC_TRIM_BYTES`'s own docstring includes
   *"allocator retention"* — and allocator retention is exactly the pool-proportional fragmentation §12
   measured. If part of the 205 MiB is really pool-proportional, it shrinks under the change too and `f`
   is conservative (safe direction). If the reverse — a `cores`-scaling-but-`B`-independent term hiding
   inside `peak − static` — `f` under-delivers. §2.3's second-budget cell settles which.
2. **Small budgets.** The insensitivity is a property of the *estimate*, and it is strong because
   `prediction − static` is large (1656 MiB) at the cells measured. At a budget near `min_required`
   (`prediction − static` ≈ 90 MiB) the same derivation would be ~340× more sensitive to `static`. That
   does not invalidate the chosen `f` (which is budget-invariant under the model), but it does mean the
   insensitivity argument gives no cover at all near the floor — and the floor is where Validation 9
   lives.
3. **`--fastqc`.** `static_bytes_for` adds up to 384 MiB, so a `--fastqc` run's `static` is 224–608 MiB —
   above the range the table sampled and near the 434 MiB tolerance. The charge scales with it correctly
   by construction, and the CI's `1G_c4_fastqc` cell is gated on the promise at 0.707, so the risk is
   low. But Validation 9 should include a `--fastqc` cell, since that is where `min_required` is largest
   (c16 + `--fastqc`: 840 → 908 MiB, so `--memory 900M --cores 16 --fastqc` loses clumpify).

### 3.2 Assumption 8 (linearity in `B`) — partly measurable from data already collected

§12 varied flush cycles 8 → 16 → 32 at fixed `B` and got 1973.0 / 1973.0 / 1974.0 MiB. That is direct
evidence that the fragmentation factor **saturates in cycle count**, which is most of what Validation 6
worries about (~39% more cycles). The plan doesn't credit it.

What §12 did *not* report is the diagnostic it already computed: `spike_occupancy.py:124` does
`slope, _ = lstsq(bs, peaks)` over three budgets (12, 22, 40 MiB) and **discards the intercept**. A
near-zero intercept and small 3-point residuals would be first-order evidence that peak is proportional
to `B`; a large intercept would say it is not. Print it — it costs nothing and it is the only
`B`-linearity evidence in the whole plan.

### 3.3 Unstated: the constants are calibrated paired-end and applied to single-end unchanged

Every cell in §4/§10 is `--paired`. `SingleBin`'s slot is `72 + 4 = 76 B` against `bytes_this = r`, i.e.
`76/r` versus the paired `148/(2r) = 74/r` — 2.7% apart, and the per-byte allocation count is identical
(3 per `r` bytes either way). So SE almost certainly transfers, and I'd expect no problem. But §7 should
say it rather than leaving a reader to derive it, and §5 has no SE row.

### 3.4 The floor warning's threshold and its stated justification

Two small things. (i) §8 Step 6 warns below **96 B** while §6f names the calibrated floor as **96.6 B**
(measured 96.642) — 0.64 B of calibrated range is left uncovered by a plan that is otherwise fussy about
exactly this. Use 97, or the measured value. (ii) The justification *"`--small_rna` reaches there"` is
weak: the clump dispatcher bins **untrimmed** input records (`parallel.rs:844-845`, reader-side), so a
small-RNA library's `r` is its raw 50 bp+ read, ~105 B and above the floor. Real Illumina headers run
45–60 B, so at 25 bp real data sits at `r ≈ 105` and never warns. The honest case for the warning is
already-trimmed FASTQ being re-processed, and synthetic/minimal-header input. Say that — it is still
worth having, but not for the reason given.

---

## 4. Efficiency and second-order effects

`B` falls 28.1%, so flushes per bin rise 1.391× (+39.1%) — the plan's figure is right. What follows:

| effect | assessment |
|---|---|
| **gzip member count** | +39%. Members are ~18 B of framing each, negligible in bytes; the real cost is the 32 KB window resetting 39% more often. |
| **output size** | see §2.8 — measured harmless on a 100%-distinct fixture, unmeasured on the duplicate-rich data clumpify targets. **The only unquantified user-visible cost.** |
| **CI band** | §2.1 — the plan has this exactly backwards. |
| **throughput** | unmeasured and probably fine. Sort work *falls* slightly (`N·log(n/1.391) < N·log n`); 39% more deflate stream inits is the cost. The TSV has no within-cores `B` contrast so nothing can be inferred, but it already records `wall_s`/`cpu_s`, so Step 9 gets this for free — add it. |
| **`min_required`** | +4.8% at c2 (280.5→293.8 MiB), +16% at c16, +23% at c32. Well under the 1 GiB default, so only explicit small budgets regress, and they *warn and fall back to plain mode* rather than bail. Sharpest case is `--fastqc` (§3.1). |
| **exact-inverse tolerance** | 312 → 434 B at cores ≤ 32 (§2.15). |
| **`--clump_only`** | shares `dyn_denominator`, `f = 1.3912` identically. Zero Phase 4 cells measured it; it is single-threaded with `workers = 1`, so fewer in-flight copies makes a uniform `f` conservative there — but it inherits a 28% bin cut, hence an unmeasured compression cost, on no evidence of its own. Validation 10 checks consistency, not correctness. Worth one sentence of acknowledgement. |

**A trap in Step 10 worth naming explicitly.** §6d's escape hatch is *"anyone on long reads who wants the
old bin size raises the budget"*. If the docs say "raise `--memory` by ~1.39×" they will be **wrong, by
the exact error §9 names as its sharpest lesson** — `static` does not scale. To restore today's `B` at
1 GiB/c2 you need `--memory 1328M`, i.e. **+29.7%**, not +39%; at 1424 MiB (1.391×) you overshoot to
24.8 MiB. Give the formula or the worked value, and note it is budget-dependent.

**One free win the plan misses.** The §6e re-baseline runs the two banded cells on **Linux** in CI
anyway. Their predicted post-change centres are **0.619** (`1G_c4`) and **0.709** (`25bp_1G_c4`) at
`static = 173`. If the PR's runs land near those, that is a Linux, independent, zero-extra-cost partial
confirmation of both Validation 6 (linearity in `B`) and Validation 7 (platform agreement). Write the
predictions down before the run so the check is a prediction rather than a rationalisation.

---

## 5. On §9's lesson about degenerate designs

> *"Where a design is degenerate, fit harder and you will get a plausible number for every prior you
> bring; the fix is an experiment whose answer does not depend on the prior."*

The diagnosis is right and the instinct — change instrument rather than estimator — is the correct
takeaway. Two ways it is drawn too broadly:

1. **The fit was not uniformly untrustworthy; it was honest about which directions were identified.**
   `s` came out at 105–109 at *every* static pin from 173 to 300 MiB (I confirm it is exactly invariant,
   because the two `cores` columns span the constant vector at a single budget), and the harness *aborts*
   rather than returning junk for `a`/`b`/`static`. The identified direction is what falsified `s = 74`.
   The failure was reading the *unidentified* coefficients as quantities. So the lesson is narrower and
   more useful: **know which directions of your design are identified, and quote only those** — which
   `calibrate.py`'s rank report already implements.
2. **"An experiment whose answer does not depend on the prior" is not what the spikes were.** They traded
   a statistical prior for a **configuration** prior — single-threaded, batches dropped at flush, no
   cross-thread frees — and that prior was wrong by 1.22–1.75× (1.37–1.88× after §2.7). §12's own
   limitations say so. The accurate version: *a slope removes the confounded intercept but not the
   configuration gap; it can settle a structural question (is the spine exactly 148?) and cannot settle a
   magnitude.* As written, §9's lesson would license the next reader to take a spike number into
   `clump.rs`, which is the one thing §12 spends a section forbidding.

Also: §9's supporting sentence — *"the slope gave 148.0 for the spine with zero spread and then
decomposed the rest"* — overstates. The spine was zero-spread; "the rest" (~120 B/pair) varies 10% across
read lengths and its cycling counterpart was ~2× different. The clean result and the messy one are being
quoted with one voice.

**On the boundary §12 draws (does a spike number leak into §6?) — it does not.** I checked every figure
in §6b/§6c against the TSV: all of it is raw `peak`, `prediction` and `static`. The spikes appear in §6a
and §5 rows 3–4 only as qualitative grounds for killing D and the allocation-count variant, which is the
correct use. The 96 B threshold comes from the matrix fixture, not a spike. §6 respects the boundary.

**On §5 rows 1–4 being "already satisfied"** — rows 1 and 2 hold as stated (I re-ran the parse logic and
reproduced the falsification). Row 3 holds for the spine and mis-states the total (§2.11). **Row 4 is
marked "done" without its qualifier**: §12 answered it in an isolated model whose magnitude explicitly
does not transfer, and the additive form `slots + fragmentation(r)` was *imposed* on three data points
rather than distinguished from the equally consistent multiplicative form `F(r)·(a + s_full/r)` — §12's
own `peak/pool` column (1.33 / 1.12 / 1.04) is the multiplicative reading. Nothing downstream depends on
it, so this is a labelling fix: "done — mechanism attributed, magnitude not transferable".

---

## 6. Action items

### Critical

1. **§6e and Step 8 — reverse the CI-band conclusion and narrow the width.** Lever A's magnitude is
   proportional to `B` and falls 28% with it, so at width 0.09 the margin against the only calibrated
   regression collapses from 0.036 to ~0.001–0.002. Narrow to ≈0.065 (3× half-range is 0.047–0.059, so
   both conditions in the baseline file hold), and re-measure the depth revert once to confirm rather
   than infer. §2.1.
2. **§8 Step 4c/5 — re-scope the change.** Five tests fail (`resolve_layout_default`,
   `resolve_layout_scales_with_cores`, `resolve_layout_pins_the_default_cell`,
   `sizing_constants_are_pinned_across_the_role_space`,
   `clumpify_floor_at_the_default_cell_is_pinned_absolutely`), ten hard-coded expectations move, and
   `channel_depths_sum_to_the_per_worker_coefficient` keeps passing while asserting a relationship
   (`k = 64/16 = 4`) that becomes false. Update that test, `channel_depths`'s docstring,
   `resolve_layout`'s docstring including its provenance paragraph, and `MARGIN_NUM`'s "96 MiB". §2.2.

### Important

3. **Run one cell before committing the integers: 25 bp, `--cores 8`, `--memory 1G`** (~3 min). It
   separates the budget axis from the cores axis — currently confounded, with the highest `f` at the
   highest budget measured — and a second point at 512 MiB additionally identifies `static` and tests
   linearity in `B`. §2.3.
4. **Fix Validation 6's tolerance.** 19 MiB is the fit's residual RMS; the measurement's within-cell
   spreads are 3.9–64.1 MiB. Set the band from the post-change spread measured in the same runs. §2.4.
5. **Fix `calibrate.py`'s `design_row` to use exact `bin_byte_budget`, not `b_banner`** — before the
   Linux run. `a` shifts −14%, `b` +16%, and one falsification cell flips sign, so §10's "all in the
   unsafe direction" and §12's 1.22–1.75× both need correcting. State explicitly that §6's integers are
   unaffected. §2.7.
6. **Add an output-size validation row** and restore the compression risk v3 dropped. Two runs on a real
   RRBS/WGBS fixture (`--memory 1G` vs `--memory 1328M`) isolate `B`'s effect on savings; the Phase 4
   fixture is 100% distinct sequences and cannot see it. Fix the clumpy.md caution box, which currently
   exempts the 1 GiB default. §2.8.
7. **Arm the CI gates this refit exists for**: `25bp_1G_c4` and `25bp_1G_c2` from `gating 0` to `1`, and
   delete the `#457` "until that refit lands" comment. Land the src change and the baseline re-record in
   one commit. §2.10.
8. **Fix §6c's `B new` column** (16.392 / 13.031 / 15.095 / 7.548 MiB, banner 16/13/15/7) — it currently
   scales the truncated banner value, the exact thing Step 5 forbids. §2.5.
9. **Add the two missing 100 bp cells** to §6b and §6c, or stop claiming "every measured cell". §2.6.
10. **Rewrite §9's stale first block.** Its "Remaining risks" are five items about withdrawn Option D and
    its "now verified" list includes two withdrawn claims. Qualify §2.3's "zero-parameter confirmation",
    and correct §10's "8 of 14 … on every rep" (51 bp c2 is 1 of 3) and "breach 150–250 MiB" (actual
    2.2–558.1). §2.11.
11. **Add the `size_of` pin as a numbered step.** §2.3, §7 Assumption 4 and §11 rec 1 all call for it;
    §8 contains no step. §2.12.
12. **Step 12 must cover four `reserve_exact` sites**, not two — `clump_only.rs:104` and `:147` duplicate
    `parallel.rs`'s bins. §2.13.
13. **Extend Validation 9 to a `--fastqc` cell**, where `min_required` is largest (c16 + `--fastqc`:
    840 → 908 MiB). §3.1.

### Optional

14. State the margin rule that produces `f ≥ 1.379`, and add the two-promise observation: the binding
    cell lands at 0.877 of `--memory` post-change, so the residual risk is a printed-prediction breach,
    not a budget breach. §2.14.
15. Strengthen Assumption 6 to its real content: `f = 1.391` tolerates `static` up to ≈434 MiB, and the
    binding cell is stable across 173–450 MiB. §3.1.
16. Print the discarded intercept and 3-point residuals from `spike_occupancy.py:124` — free evidence for
    Assumption 8. Credit §12's 8/16/32-cycle invariance as partial support. §3.2.
17. Record the predicted Linux band centres (0.619, 0.709) before the re-baseline run, so it doubles as a
    zero-cost Linux check on Validations 6 and 7. §4.
18. Fix §5 row 8 (no test asserts 312 B; the tolerance becomes 434 B), row 3's 255–270 → 255–280, and
    row 4's "done" → "mechanism attributed, magnitude not transferable". §2.15, §2.11, §5.
19. Move the floor warning to 97 B (calibrated floor is 96.642) and replace the `--small_rna`
    justification — the dispatcher bins untrimmed reads, so small-RNA input sits above the floor; the
    real case is re-processed already-trimmed FASTQ. §3.4.
20. State in §7 that the constants are calibrated paired-end and applied to single-end unchanged (76/r
    vs 74/r, 2.7% apart), and add a `--clump_only` acknowledgement: it inherits `f` with no cell of its
    own. §3.3, §4.
21. Docs (Step 10): give the `--memory` escape hatch as +29.7% at 1 GiB or as a formula, not 1.39× —
    scaling a ratio with a fixed `static` part is the error §9 names as its sharpest lesson. §4.
22. Narrow §9's lesson: name the identified-vs-unidentified distinction, and say that a slope removes the
    confounded intercept but not the configuration gap. §5.
