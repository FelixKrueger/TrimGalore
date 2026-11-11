# Plan Coverage Report — Phase 4 (#457) v4

**Mode:** B — code vs. implementation plan (§8 "Implementation outline" + §5 "Validation" of `PLAN.md` v4)
**Plan:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` (v4)
**Code:** branch `perf/clumpify-spine-accounting`, uncommitted working tree, diffed against `origin/dev` @ `3fa4de3`
**Date:** 2026-08-20
**Verdict:** **INCOMPLETE — 2 MISSING, 5 PARTIAL, 4 documented deferrals, 1 pending CI**

> **The audited tree was being modified during this audit.** `src/parallel.rs` and
> `src/clump.rs` changed under me at 08:09–08:10 (local), while I was reading them, and
> `src/clump_only.rs` at 08:10:18. Two defects I had found were fixed by that writer before I
> finished; both are recorded under "Defects found, fixed mid-audit" below rather than dropped,
> because they were real in the tree I read. Everything in this report was then re-verified
> against the tree as of **08:12**, and that state is snapshotted at
> `/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/f5588c8e-a7ba-45ff-a13c-b4cc5aa3a1c0/scratchpad/audited_diff.patch` with per-file SHA-256 in `/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/f5588c8e-a7ba-45ff-a13c-b4cc5aa3a1c0/scratchpad/audited_hashes.txt`. If the
> tree has moved again since, re-check against those hashes before acting on this report.

## Summary

- Total ledger items: **32** (§8: 4c.1–4c.10, 5.1–5.4, 6.1–6.2 = 16; §5: rows 1–16)
- DONE: **20**
- PARTIAL: **5** — 4c.1, 4c.9, 4c.10, Validation 5, Validation 8
- MISSING: **2** — Step 5.4 (CHANGELOG), Validation 12
- DEFERRED, documented in §13: **4** — Validations 7, 9, 10, 14
- PENDING CI, *not* documented in §13: **1** — Validation 16

Build state verified independently, not taken from §13:

- `cargo test` — **805 pass, 0 fail** (506 lib + 299 across 15 integration binaries). §13's figure reproduces exactly.
- `cargo clippy --all-targets --release -- -D warnings` — **exit 0, zero diagnostics**, re-run after `touch`ing the five changed sources so it was a real check and not a cache hit.
- `cargo build --release` — clean; the binary was then used to confirm every published figure.

## Coverage ledger — §8 Implementation outline

| # | Item | Status | Notes |
|---|------|--------|-------|
| 4c.1 | Fix `calibrate.py` to regress on exact `bin_byte_budget` | **PARTIAL** | `design_row` fixed; three other consumers still read the floored banner value. Detail below. |
| 4c.2 | `dyn_denominator`: `23` → `33`, `64` → `92` | DONE | `src/clump.rs:174-176`. |
| 4c.3 | Ten hard-coded expectations recomputed off the binary | DONE | All six tests verified; all pinned literals independently recomputed. Sub-table below. |
| 4c.4 | The test that stays green while its invariant becomes false | DONE | Renamed and re-scoped; see note on residual tripwire strength. |
| 4c.5 | Floor warning at `r < 97 B` | DONE | Verified live at 68 B. Three notes below. |
| 4c.6 | Pin `2·size_of::<FastqRecord>() + size_of::<MinimizerKey>() == 148` | DONE | Compile-time `const _`. Introduces one collateral defect, below. |
| 4c.7 | Three stale docstrings + one stale comment | DONE | All four verified; no stale coefficient survives in `src/` or `docs/`. Sub-table below. |
| 4c.8 | Arm the CI promise gates, delete the `#457` comment | DONE | **The deliverable.** Both cells `gating 1`; the "until that refit lands" wording is gone from `ci.yml` and the baseline; same working tree as the src change. |
| 4c.9 | Narrow the band to 0.065 / −0.072, re-baseline, measure lever A | **PARTIAL** | Enforced `awk` and the two banded rows updated; the rest of the file's analysis still judges at the retired 0.09. Re-baseline deferred via `# BOOTSTRAP`. |
| 4c.10 | Re-measure the 25 bp cells | **PARTIAL** | 3 reps on the binding cell where Validation 5 asks for ≥10; Validation 10's two extra cells not run (documented). |
| 5.1 | Re-derive every docs figure against the binary | DONE | Every figure independently reproduced from the binary — list below. |
| 5.2 | Long-read cost + `--memory` escape hatch as budgets, never a multiplier | DONE | And the docs **corrected a second plan-table error**; see below. |
| 5.3 | Fix the caution box; say nothing about compression cost | DONE | 1 GiB exemption removed; silent in both directions. |
| 5.4 | CHANGELOG, one sentence, cumulative about output bytes | **MISSING** | `CHANGELOG.md` is untouched, and its existing `#439` entry now carries a false claim. |
| 6.1 | Cap the spine extrapolation at four `reserve_exact` sites | DONE | All four verified. |
| 6.2 | Variable-length test for the cap | DONE | Two tests; coverage is one of the four sites. |

### 4c.3 — the six named tests, audited individually

Every pinned literal was recomputed from the formula independently of the plan's cross-check table and independently of the test suite, then confirmed against the passing binary.

| test | asserts after `(33, 92)` | recomputed | Status |
|---|---|---|---|
| `resolve_layout_default` (`clump.rs:563`) | `B ∈ [78, 80)` MiB, comment names `denom = 712` | 78.64 MiB, denom 712 | DONE |
| `resolve_layout_scales_with_cores` (`:573`) | `B ∈ [31, 33)` MiB, comment names `denom = 1792` | 31.25 MiB, denom 1792 | DONE |
| `resolve_layout_pins_the_default_cell` (`:582`) | `B == 13_236_569`, comment names `896` | 13,236,569 | DONE |
| `sizing_constants_are_pinned_across_the_role_space` (`:591`) | six exact budgets | 13,236,569 / 19,128,978 / 6,046,334 / 6,046,334 / 16,657,256 / 8,191,148 — **all six exact** | DONE |
| `clumpify_floor_at_the_default_cell_is_pinned_absolutely` (`:640`) | `floor.div_ceil(MiB) == 296` | exact floor 295.35 MiB → `div_ceil` 296 | DONE |
| `resolve_layout_fastqc_reserves_more_and_shrinks_bins` (`:664`) | denominator `896` | 896 | DONE |

**§13's account of the 1 MiB discrepancy is accurate.** §8 4c.3's cross-check table projected `→ 295`; the exact floor is 295.35 MiB and `div_ceil` gives **296**. §5 Validation 12's own figure (295.4) was right; 4c.3's table was the wrong one. The committed value follows the binary.

### 4c.7 — the three docstrings and one comment, audited individually

| location | was | is | Status |
|---|---|---|---|
| `dyn_denominator` docstring (`clump.rs:172-173`, fn at `:174`) | `σ = 23/16 … k = 4` | `σ = 33/16 … k = 92/16` | DONE |
| `channel_depths` docstring (`:178-179`) | *"their sum is the constant and neither may move alone"* | *"`k` charges for these slots plus the allocator's per-record cost, so a deeper queue raises peak beyond `k`"* | DONE |
| `resolve_layout` docstring (`:267-291`, fn at `:292`) | `23/16`, `k = 4`, the literal `23 × n_bins + 64 × workers`, and the four-clause provenance paragraph | all four clauses rewritten; provenance now names the byte-equalised matrix, the `r = 96.642 B` floor, and the phase 4 plan | DONE |
| `MARGIN_NUM` comment (`:262-263`) | *"largest observed within-cell spread is 96 MiB"* | 64 MiB | DONE |

Swept for survivors: `23 ×`, `64 ×`, `23/16`, `64/16`, `k = 4`, `denom = 624`, `= 496`, `= 1248` appear nowhere in `src/`; the retired figures `281/277/590/454/933/1356/334/311/775 MiB` and `309806` appear nowhere in `docs/`. No other docs page carries an affected figure (checked all seven pages mentioning clumpify; `benchmarks.md`, `quality.md`, `outputs.md`, `passthrough.md`, `library-presets.md` publish no memory numbers, and `README.md` is qualitative only).

## Coverage ledger — §5 Validation rows

| # | Row | Status | Notes |
|---|-----|--------|-------|
| 1 | Harness reports its own failure | DONE | §10, pre-existing. |
| 2 | `s` matches `size_of` | DONE (falsified) | Coefficient recompute deferred, documented. |
| 3 | Per-slot cost measured directly | DONE | And now pinned in code by 4c.6 — the spike's 148 B is a compile-time invariant. |
| 4 | Whether the overhead rides occupancy | DONE | §12. |
| 5 | **gate** — promise holds at the calibrated floor, ≥10 reps on the binding cell | **PARTIAL** | §13 reports 3 reps (0.923 / 0.909 / 0.907) where the row requires **≥10**; 25 bp c2 has 1 rep (0.821). All under prediction, all where every rep was over it before — the promise is demonstrated, the rep count is not met, and no run log is committed. |
| 6 | **gate** — CI promise gate armed | DONE | Verified in `ci.yml:441,446` and the baseline. |
| 7 | **gate** — band still detects lever A at width 0.065 | DEFERRED | Documented; needs CI. Reason holds — but see the follow-up gap below. |
| 8 | Linearity in `B` held | **PARTIAL** | §13 answers it in prose ("approximate and favourable", explicitly *not* a verified identity) with no committed data. Honest, and short of the row's "compare against the post-change measured within-cell spread from the same runs". |
| 9 | Linux agrees, or is worse | DEFERRED | Documented. The plan marks this **blocks merge**; it still does. |
| 10 | Budget axis separated from cores | DEFERRED | Documented. |
| 11 | `predicted_peak_bytes()` still the exact inverse | DONE | The one-sided test at `clump.rs:706-713` exists and passes; `312` appears nowhere in `src/` or `tests/`, as the row states. Aside: the row's "~434 B at cores 32" is wrong — the residual is bounded by `denom/16`, which at cores 32 uniform is 7168/16 = **448 B**. Nothing asserts it, so no code consequence. |
| 12 | Minimum viable budget moved; failure is a silent capability loss | **MISSING** | Not mentioned anywhere in §13. Audited here — see the gap detail. Two of its three max-cores figures do not reproduce. |
| 13 | `--clump_only` moved with it | DONE | `clump_only.rs` calls `clump::resolve_layout` at four sites, so it inherits the coefficient; docs updated; its three banner tests still pass. |
| 14 | Output size on a real library | DEFERRED | Documented, and Step 5.3 correctly leaves the docs silent in both directions meanwhile. Cleanest of the four deferrals. |
| 15 | The banner string survives | DONE | Three tests assert `"predicted peak ≈ 930 of 1024 MiB budget"` verbatim (`integration_clump_only.rs:126`, `integration_clump_only_ubam.rs:169`, `integration_preflight_tripwire.rs:605`) and pass. Confirmed live in the banner. |
| 16 | Perl byte-identity untouched | PENDING CI | Not in §13's "not done" list, so it is an **undocumented** deferral. Its premise is now verified locally: the `validation` job contains **zero** `--clumpify` invocations, so no validated flag path touches the changed code. |

## Gaps (detail)

### Step 5.4 — CHANGELOG entry (MISSING)

**Expected:** one sentence, cumulative about output bytes moving, naming this as the third such change in the cycle (#442, #446, now).

**Found:** nothing. `git diff origin/dev -- CHANGELOG.md` is empty and `git status` does not list the file. §13 does not list 5.4 under "Not done, and why", so this is a silent drop, not a deferral.

**Compounding:** the existing `#439` entry (`CHANGELOG.md:69-79`) still ends:

> Very short reads are not yet covered: the per-bin coefficient is fitted per byte rather than per record, so a 25 bp library packs more records into the same bin bytes and can exceed the budget at `--cores 2` ([#457]).

This change makes that false. The shipped changelog would tell a reader the defect is still open while the binary has fixed it. **Both the new entry and the correction to the old one are outstanding.**

### Validation 12 — the moved floor and the silent capability loss (MISSING)

**Expected:** confirm the max cores keeping clumping engaged at 1G moves 72 → 52, with `--fastqc` 33 → 23, `--clump_only` at 512M 41 → 29; floors c2 280.5 → 295.4, c16+fastqc 840.4 → 915.2, c32+fastqc 1012.0 → 1161.6; and that `--clumpify` warns-and-falls-back while `--clump_only` errors.

**Found:** no mention in §13. Audited here by reimplementing `clumpify_min_memory_bytes` / `resolve_layout` exactly and sweeping the core axis:

| quantity | plan says | reproduces as |
|---|---|---|
| max cores, 1G, no `--fastqc` | 72 → **52** | 72 → **50** |
| max cores, 1G, `--fastqc` | 33 → 23 | 33 → 23 ✓ |
| max cores, `--clump_only` 512M | 41 → **29** | 41 → **28** |
| floor c2 | 280.5 → 295.4 MiB | ✓ (`div_ceil` 296) |
| floor c16 + fastqc | 840.4 → 915.2 MiB | ✓ (916) |
| floor c32 + fastqc | 1012.0 → 1161.6 MiB | ✓ (1162) |

The behavioural half of the row is verified live: `--clumpify --memory 100M --cores 6` prints `WARNING: … (need ≥ 339 MiB). Falling back to plain mode` and exits 0; `--clump_only --memory 700M --cores 16 --fastqc` errors out. **None of the max-cores figures is published in code or docs, so nothing shipped is wrong** — the gap is that the row was never worked and two of its numbers are off by 1–2 cores.

### 4c.1 — the harness fix is partial

**Expected:** restore an exact `bin_byte_budget` column and regress on it, because regressing on the banner's floored MiB is a ~5% error in the regressor at the c16 cell.

**Found:** the core of it is done and correct. `exact_bin_mib()` added at `calibrate.py:110`, `SIGMA16, K16 = 33, 92` (`:88`) so it tracks the new binary, `b_exact` computed per row (`:316`), `b_exact_mib` written to the TSV (`:427`, `:468`), and `design_row` regresses on `r["b_exact"]` (`:493`) with a docstring that says *"`b_exact`, never `b_banner`: the banner floors to whole MiB."*

**Gap:** three consumers immediately below that docstring still read `b_banner`:

- `:510` — the rank diagnostic's `distinct_cores_b` set. Harmless-to-conservative.
- `:546` — the **zero-parameter check**, which multiplies `b` by `A_EXPECTED`/`B_EXPECTED`/`S_EXPECTED` to predict RSS.
- `:556` — the **falsification block**, which predicts the held-out cells from the fit.

The last two are model predictions scaled by `B`, so they carry the same ~5% error 4c.1 exists to remove — and 4c.1 specifically named "one falsification cell flips sign" as a consequence of that error. Leaving the falsification comparison on the floored value means the sign-flip cannot be checked from this harness. 4c.1 also says *"Do this before the Linux run so both platforms use the fixed harness"*, so this interacts directly with the Validation 9 deferral.

### 4c.9 — the band narrowing is only half-applied to the prose

**Expected:** narrow to 0.065 / −0.072 "in the enforced `awk` and the prose", re-baseline both cells, and measure the post-change lever-A delta.

**Found:** the enforced `awk` is correct (`ci.yml:428-431`), the file header explains the narrowing and its `B`-proportionality reason, and the two banded rows are annotated `pre-refit; at +0.065 post-refit: …`.

**Gap:** the rest of the file still judges at the retired width, unmarked:

- `:45` the column header still reads `at +0.09`, now contradicting the two rows beneath it.
- `:51` `multipair_c4` — "1.2x at +0.09; needs +0.121 for 3x".
- `:58` `1G_c4_fastqc` — "1.2x at +0.09".
- `:73` `2G_c4` — "Lever A moves it ~+0.09 to +0.13".
- `:77` `1G_c2` — "2.9x at +0.09 — just under the bar". At 0.065 against a half-range that also scales by 1.435 this becomes ~3.8x, i.e. it would now **clear** the bar it is recorded as missing.
- `:78` `25bp_1G_c2` — "4.2x at +0.09 — eligible next pass". Post-refit this is ~3.0x, i.e. **right at** the bar, not comfortably eligible.
- `:85` the diagnostic rule of thumb — "a depth revert moves every cell together by 12 points or more (measured)" — becomes ~8.4 points under the same proportionality.

So the two c2 verdicts are now wrong **in opposite directions**, and the diagnosing-a-red heuristic overstates the signal by 1.4x. These are the judgements the next person will use to decide what to band, and §13's claim that the narrowing landed "in the enforced `awk` and the prose" overstates what was done.

### 4c.10 / Validation 5 — rep count and evidence

Validation 5 is written as a **gate** and specifies "≥10 reps on the binding cell (25 bp c8 2G, 57 s each)". §13 reports three. The three are unanimous and comfortably under prediction, and §6b's margin rule was deliberately sized from the worst spread anywhere precisely because the binding cell has few reps — so the *direction* is not in doubt. But the row's own reason for demanding ten is the arc's recorded lesson that dispersion from n=2–5 is biased low, and that lesson applies to these three readings too.

There is also **no committed artifact** for the re-measurement: `phase4/` holds no TSV or log, so §13's four numbers are prose only. I can confirm the code, every published figure, and every test; I cannot independently confirm those four measurements. Given that this arc has already had a fabricated peak-RSS figure survive four review rounds, that is worth stating rather than assuming.

## Judgement on the four documented deferrals

The task asked whether each deferral is genuinely documented *and* whether the reason holds.

| deferral | documented? | reason holds? | verdict |
|---|---|---|---|
| Validation 7 (lever-A delta) | yes, §13 | yes — it needs a Linux CI run | **Accept, with a follow-up gap.** The `# BOOTSTRAP` marker means nothing is falsely armed. But 0.065 was chosen from the `B`-proportionality *argument*, which §6e itself says "gives the direction, not the value", and the BOOTSTRAP note in the baseline tells the next person only to "record the ratios from this PR's own Linux run, replace both rows, and delete this line" — **it does not tell them to run the depth revert first.** Delete the marker on ratios alone and the band arms at a width validated by inference, not measurement. |
| Validation 9 (Linux) | yes, §13 | yes | **Accept — and it still blocks merge.** The plan says so twice (§5 row 9, §6f). Note the interaction with the 4c.1 gap: the Linux run would use a harness whose falsification and zero-parameter blocks still consume the floored `B`. |
| Validation 10 (budget vs cores) | yes, §13 | yes, narrowly | **Accept, but "neither gates the change" understates it.** §6f lists the budget/cores confound as a live gap, and the adopted `f` comes from the 25 bp c8 **2G** cell — the one cell the confound touches, and the highest budget measured. The promise at that cell was re-measured directly, so the change is not unverifiable; what stays unresolved is the provenance of the chosen integers. |
| Validation 14 (output size) | yes, §13 | yes | **Accept, cleanly.** The deferral is honoured in the artifact: Step 5.3 leaves the docs silent about compression cost in both directions, which is exactly what §6d demanded. |
| `# BOOTSTRAP` band state | yes, §13 and §6e | yes | **Accept.** The mechanism pre-exists, fails safe (band off, loud warning, and deleting the marker without re-recording reds immediately), and §6e sanctioned it by name. State the trade honestly though: this PR ships with **no drift band at all** where `dev` had two cells armed at 0.09 — traded for the stronger `peak ≤ prediction` gate now armed on the two 25 bp cells. Net coverage is better, but the drift half is temporarily zero. |
| §10/§12 coefficient recompute | yes, §13 | yes | **Accept.** The claim "§6's integers never depended on it" is structurally checkable and true: `f_required = (peak − static)/(prediction − static)` references no `B`. |

## Defects found, fixed mid-audit

Both were real in the tree I read and both were repaired by the concurrent writer before I
finished. Recorded so the record shows they were caught, and so a re-run of this audit does not
report them as new.

**1. `estimated_record_bytes` had lost its doc comment to the 4c.6 pin.** The `const _` assert was
inserted *between* the existing doc comment and the function it documented, so all four `///`
lines documented the anonymous const and the public function was undocumented. **Now fixed** —
`src/clump.rs:354-355` gives the const its own two-line comment and `:361-362` restores
`estimated_record_bytes`'s. The pin itself was sound throughout and is genuinely mutation-proof:
`FastqRecord` is three `String`s (72 B) and `MinimizerKey` is `u32`, so 2×72+4 = 148 and any added
field breaks the build.

**2. `PairedBin::push`'s docstring was stale after Step 6.1.** It still read *"we extrapolate
expected record count from this record's text size and the bin's budget"*, which Step 6.1 made
false. **Now fixed** — `src/parallel.rs:898-900` reads *"expected pair count from the previous
cycle's mean pair size and the bin's budget"*, which also corrects the record-vs-pair unit that the
original got wrong. `clump_only.rs`'s two field comments were polished in the same pass to
"Mean bytes per reserved slot", likewise the correct unit.

## Notes on items marked DONE

Recorded so a reader knows what was and was not checked, not as gaps.

**5.1 — every published figure reproduced from the binary.** I reimplemented the sizing and floor arithmetic independently and, where possible, ran the binary:

| figure | source | verified |
|---|---|---|
| `clumpy.md` floors, cores 2–32, no fastqc: 296 → 740 MiB | recomputed | ✓ |
| with `--fastqc`: 349 → 1162 MiB | recomputed | ✓ |
| cores 64: 1232 MiB / 1655 MiB | recomputed | ✓ |
| `32 bins × 6 MiB` at c8/1G; `× 31 MiB` at c8/4G | recomputed 6.31 / 31.25 MiB | ✓ |
| `need ≥ 339 MiB` at cores 6 | **run live** | ✓ verbatim |
| floor warning at 68 B | **run live** on a terse-header 25 bp fixture | ✓ fires |
| `clump-only.md` floors 290 MiB (c1–4) → 544 MiB (c32) | recomputed | ✓ |
| `clump-only.md` error: `215908 bytes`, `try ≥ 821 MiB`, `608 MiB static` | **run live** at `--memory 700M --cores 16 --fastqc` | ✓ verbatim |
| `~30%` bin cut at every budget | recomputed: 30.35–30.36% across all six pinned cells | ✓ |

**5.2 — the docs caught a second plan-table error that §13 does not record.** §6d and Step 5.2 both give the 2G/c8 escape hatch as **2833M**; the docs publish **2834M**. The docs are right: at 2833 MiB the new formula yields `B` = 22,014,590 against the old 22,021,658, i.e. it does *not* restore the bin size, and 2834 MiB does (22,023,101). So the "re-derive against the binary" discipline caught two plan errors, not the one §13 reports.

**5.2 — two published ratios are predictions, not measurements.** *"at 150 bp a run peaks around 0.61 of its printed prediction, against 0.94 at 25 bp"* quotes §6c's predicted column. The 150 bp cell was not re-measured, and the 25 bp binding cell actually measured **0.923**, not 0.94. Both are hedged with "around", and the direction is conservative (measurement came in below prediction), so this is a wording observation rather than a wrong claim.

**4c.4 — the tripwire is weaker than it was, as the plan permitted.** The renamed `channel_depths_are_pinned_to_the_calibrated_values` asserts `(work, result_per_core) == (1, 2)` and explains what it now guards, which satisfies 4c.9's second option. What it no longer does is tie the depths to `k₁₆` by any arithmetic — so if `k₁₆` moves again, this test stays green. That is inherent to the option chosen and is not a deviation.

**4c.5 — three observations.** (a) No test covers `warn_below_clumpify_floor`; a regression that silenced it would pass CI. The plan did not ask for one. (b) The ~9% extrapolated breach that 4c.5 quantifies is published in neither the warning nor the docs — the substance ("peak memory may exceed the printed prediction" / "warns rather than silently extrapolating") is there, the magnitude is not. (c) The warning fires before the below-floor fallback decision, so `--clumpify --memory 100M` prints a calibrated-floor warning about a reservation that then does not happen. Verified live; cosmetic.

**6.2 — the cap is tested at one of its four sites.** `bin_reservation_follows_the_previous_cycle_not_one_short_record` and `bin_reservation_falls_back_to_the_first_record_before_any_cycle` both exercise `parallel::PairedBin`. `parallel::SingleBin` and both `clump_only` bins are untested. The plan asked for "its own variable-length test" (singular) and two exist, so the step is met; the coverage asymmetry is worth knowing.

## Verdict

**INCOMPLETE.** The change itself is sound and the deliverable is armed: the coefficients moved, all ten pinned expectations are correct against the binary (independently recomputed, not taken from the plan's table), the two 25 bp CI cells now gate on `peak ≤ prediction`, the spine cap is at all four sites, and every published docs figure reproduces from the binary — including two the plan's own tables got wrong. 805 tests pass and clippy is clean under `-D warnings`, both re-verified here rather than accepted from §13.

Outstanding, in the order they should be addressed:

1. **Step 5.4 — write the CHANGELOG entry, and fix the old one.** One sentence, cumulative about gzip member boundaries (#442, #446, now). Separately, the `#439` entry's closing sentence still tells readers very short reads are uncovered and can exceed the budget at `--cores 2`; that is now false and ships as-is.
2. **4c.9 — finish the band prose.** The `at +0.09` column header and six further judgements still evaluate at the retired width; the two `1G_c2` / `25bp_1G_c2` eligibility verdicts are now wrong in opposite directions, and the "12 points or more" diagnosing-a-red heuristic overstates the signal by 1.4x.
3. **4c.1 — finish the harness fix.** `calibrate.py:546` and `:556` still multiply the floored `b_banner` by fitted coefficients, so the zero-parameter check and the falsification block retain the ~5% error the step exists to remove — including for the cell 4c.1 says flips sign. Do this before the Linux run, as 4c.1 instructs.
4. **Validation 12 — work the row.** Its behavioural half is verified (warn-and-fall-back vs error) and its three floors reproduce, but two of the three max-cores figures do not (72 → 50 not 52; 41 → 28 not 29). Nothing shipped depends on them.
5. **Validation 5 — raise the binding cell to ≥10 reps** as the row specifies, and commit the run log. Three unanimous reps establish the direction; they do not meet the gate as written, and there is currently no artifact behind the four numbers in §13.
6. **Before deleting `# BOOTSTRAP`, run Validation 7.** The marker's own instructions say to record ratios and delete — they omit re-measuring the lever-A delta, so following them arms a width justified by inference alone.
7. **Validation 9 (Linux) still blocks merge**, per the plan.

Nothing in §13's "Not done, and why" was found to be misrepresented, and its two "what the plan got wrong" claims both check out. The omissions this audit adds are Step 5.4, Validation 12, and Validation 16 — none of which §13 mentions at all.
