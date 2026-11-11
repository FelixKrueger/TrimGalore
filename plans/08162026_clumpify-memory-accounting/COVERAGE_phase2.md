# Plan Coverage Report — Phase 2

**Mode:** B (design plan standing in for an implementation plan; no `IMPL.md` exists)
**Plan(s):** `plans/08162026_clumpify-memory-accounting/PLAN.md` (v3) + `phase2/CALIBRATION.md`
**Scope:** Phase 2 only. Phase 1 is merged (`1595d59`, identical in content to branch HEAD `70e75e3`); Phase 3 (CI guard, docs arithmetic) is future work and is not counted as a gap.
**Date:** 2026-08-17
**Verdict:** INCOMPLETE — 12 unresolved (8 PARTIAL, 4 MISSING), plus 4 undocumented deviations

## Provenance of what was audited

`git diff 1595d59 HEAD -- src/ tests/ CHANGELOG.md` is empty, so the uncommitted working-tree
diff is exactly Phase 2: `CHANGELOG.md`, `src/clump.rs`, `src/clump_only.rs`, `src/main.rs`,
`src/parallel.rs` — 152 insertions, 81 deletions. No file under `tests/` and no file under
`docs/` is touched.

Verified independently by me:

- `cargo test --release` → **732 passed, 0 failed** across 16 test binaries (matches the recorded number).
- `cargo fmt --all -- --check` → clean.
- `cargo clippy --all-targets --release -- -D warnings` → 0 warnings.

All arithmetic below was recomputed from the shipped constants (`STATIC_TRIM_BYTES` = 224 MiB,
`STATIC_FASTQC_PER_CORE_BYTES` = 24 MiB, `MARGIN_NUM/DEN` = 10/11, `dyn_denominator` =
`23·n_bins + 64·cores`, `/16`) and cross-checked against `phase2/results_cal_darwin.tsv` and
`phase2/results_cal2_linux.tsv` (36 data rows each, 2 reps × 18 cells).

## Summary

- Total items: 37 counted (+1 out of scope)
- DONE: 14
- PARTIAL: 8
- MISSING: 4
- DEVIATED: 11 (7 documented, 4 undocumented)

## Coverage ledger

### A. Model and code obligations

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| A1 | Denominator replaced with the calibrated form (σ, k) | CALIB item 2 | DONE | Code ships σ = 23/16 = 1.4375 and k = 64/16 = 4.0, both at or above the worse fitted values (1.41 / 3.99). Refines CALIB's coarse "roughly `6·n_bins + 16·cores`" (σ = 1.5) to 16ths; still conservative against the fit. |
| A2 | `f0` becomes ~205 MiB | CALIB item 1 | DONE | 224 MiB against a measured max of 204.7 MiB (macOS) / 173.1 (Linux). Rounded up. |
| A3 | FastQC handled as a **`max`**, not an addend | CALIB item 1; PLAN "Model form" §2 | DEVIATED (undocumented) | Shipped as an **addend**: `static = 224 MiB + 24 MiB × cores`. Conservative for safety, but it is the direct cause of the only floor regression (A/C7) and neither the notes nor CALIBRATION record that the `max` form was abandoned or why. |
| A4 | FastQC's cost "must be measured, not assumed" | CALIB item 1 + Limitations | PARTIAL | Measurement was performed (notes: 76 / 139 / 283 MiB non-pool at c4 / c8 / c16; 191.9 MiB whole-process at c2) and 24 MiB/core is conservative against all three (19.0 / 17.4 / 17.7). **No raw artifact exists**: `phase2/` holds only the two pre-implementation TSVs, and `calibrate.py` has no `--fastqc` path, so the constant's basis is prose-only and not reproducible from the committed harness. |
| A5 | Explicit named margin that **exceeds the largest observed spread** | PLAN "Staging" Phase 2; CALIB item 3 | PARTIAL | Margin is `×10/11` (9.09%), named in code and in the bail text. But at the shipped default budget the margin is **93.1 MiB against a largest observed within-cell spread of 96.2 MiB** — the stated criterion fails for every budget below ~1.06 GiB. CALIBRATION asserts the margin "exceeds the largest within-cell spread on either host" with no budget qualifier; the `MARGIN_NUM` comment places "~9%" and "96 MiB" adjacent without reconciling them. |
| A6 | "State which quantile the banner means" | PLAN "Staging" Phase 2 | PARTIAL | The `resolve_layout` docstring records max-over-2-reps fitting and the 10/11 margin. Nothing user-facing states it. Note the banner's predicted peak is now a constant **90.9% of `--memory`** in every resolvable cell (verified: 0.909 at all six shipped cells), because the prediction is the exact inverse of the sizing. |
| A7 | `predicted_peak_bytes()` is the exact inverse of the sizing | Notes §1 | DONE | Same integer ops, same constants; truncation loss ≤ `denom/16` bytes. Banner cannot disagree with sizing. |
| A8 | `resolve_layout` / `clumpify_min_memory_bytes` take `fastqc`, threaded from `Cli::fastqc_requested()` and the four `clump_only.rs` sites | Notes §3 | DONE | `main.rs:184,201`; `clump_only.rs:275,413,748,905`. `fastqc_requested()` = `fastqc \|\| fastqc_args.is_some()` (`cli.rs:592`), so the `--fastqc_args`-alone edge case is covered. |
| A9 | Confirm the falsifiable prediction that lever A takes `k` from 3σ to 2σ (expect a fall of ~1.3–1.4) | PLAN "Model form" §1 — "a falsifiable prediction the calibration must confirm" | **MISSING** | CALIBRATION measures `k` = 3.71 / 3.99 on the depth-1 (post-lever-A) binary, i.e. ≈ 2.6–2.8σ, and never compares it against the 2σ prediction; no depth-2 `k` measurement exists to test the "falls by roughly σ" claim. The nearest prior figure (reviewer B's `k` = 4.61) implies a fall of 0.62, not 1.3–1.4. The code docstring instead re-derives `k` = 4 from three resident batch-equivalents (queued batch, batch in hand, compressed output), replacing the prediction without stating that it failed. |
| A10 | Update `resolve_layout_pins_the_default_cell`, "which will move by −0.6%, deliberately" | CALIB item 4 | DEVIATED (undocumented magnitude) | Pin moves 19,884,107 → 19,006,356 B = **−4.4%**, not −0.6%. Every cell lands ~3.8% below CALIBRATION's "honest B" because `f0` (224 vs 204.7) and σ (1.4375 vs 1.41) were rounded up: 1 GiB/c2 **+4.7%** (projected +9.1%), 1 GiB/c4 −4.4% (−0.6%), 1 GiB/c8 −4.4% (−0.6%), 2 GiB/c4 **−26.2%** (−24.4%), 4 GiB/c8 **−32.4%** (−31.2%). Direction is conservative; the magnitude difference is recorded nowhere. |
| A11 | `--clump_only`'s missing warn-and-fall-back | PLAN outline item 5 | DEVIATED (documented) | Deferred again, with an explicit rationale and the concrete now-reachable case (`--clump_only --fastqc --cores 16 --memory 800M` aborts; floor 840.4 MiB). Named as a behaviour decision for the user. |
| A12 | Docs arithmetic (formula, floor range, worked values, unit labels) | CALIB item 5 | OUT OF SCOPE (Phase 3) | Not counted. Confirmed still stale: `docs/src/content/docs/performance/clumpy.md:108` says "~535 MiB … fixed 512 MiB reservation … 535–730 MiB range" against actual floors of 281–418 MiB (no FastQC) and up to 841 MiB (FastQC). The plan lists this under Phase 3 and under "Open — docs savings table". |

### B. Test obligations

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| B1 | The five named existing `clump.rs` tests updated | PLAN "Tests" | DONE | `resolve_layout_default`, `_scales_with_cores`, `_bails_when_too_small`, `_bails_below_static_overhead`, `_predicted_peak_fits_budget` — all updated, plus `_pins_the_default_cell` (six rewritten). |
| B2 | The `"static overhead"` substring assertion survives the bail-message rewrite | PLAN "Tests" | DONE | Message retains the phrase; assertion at `src/clump.rs:567`. |
| B3 | New tests for the FastQC term and the advertised floor | Notes §4 | DONE | `resolve_layout_fastqc_reserves_more_and_shrinks_bins`, `clumpify_floor_matches_the_smallest_resolvable_budget` (loops cores {2,4,8,16} × fastqc {false,true}, asserting the floor resolves and floor − 1 MiB does not). |
| B4 | Notes claim "**Nine** `clump.rs` tests rewritten or added" | Notes §4 | DEVIATED (undocumented) | Eight: six rewritten + two added. Suite moved 730 → 732, consistent with two new tests. |
| B5 | `clump_only.rs`'s `small_memory_budget()` helper and its 13 call sites | PLAN "Tests" (named as affected) | PARTIAL | Functionally fine — 1 GiB clears every new floor — but the helper was left untouched and its comment still reads "resolve_layout requires STATIC_OVERHEAD (512 MiB)" (`src/clump_only.rs:1188`), naming a constant this diff deletes. Only in-tree reference to the removed name. |
| B6 | Both `tests/integration_clump_only*.rs` default-budget invocations still resolve | PLAN "Tests" | DONE | Unchanged and passing. Floor at `--cores 2` = 280.5 MiB; the `--clump_only --fastqc_args` test at default cores clears its 302 MiB floor. Both under the 1 GiB default. |
| B7 | Banner guard in `tests/integration_preflight_tripwire.rs` | PLAN "Tests" | DONE | Banner text unchanged by Phase 2; the Phase-1 guard pinning both MiB labels still applies. |
| B8 | Retire or relocate the tautological `resolve_layout_predicted_peak_fits_budget` | PLAN "Tests" | DEVIATED (documented) | Explicitly assigned to Phase 3 in the Phase 1 post-review notes ("deleting a passing test before its replacement exists would reduce coverage"). Phase 2 extended it with two FastQC rows; it remains formula-vs-its-own-inverse, now passing with a fixed 9% slack. |
| B9 | Harness: banner-present assertion per cell, distinct-read logging, parent-footprint guard | PLAN "Tests"; CALIB "Method" | DONE | `calibrate.py`: `BANNER` regex with hard exit per cell, `n_bins` cross-check, prints reads / distinct sequences / mean length, aborts if its own RSS exceeds 200 MiB, cell summary = max over reps. |
| B10 | Suite, fmt, clippy state | Notes | DONE | Re-run by me: 732 passed / 0 failed; fmt clean; clippy 0 warnings. |

### C. Validation table (rows 1–11)

| # | Row | Status | Notes |
|---|---|---|---|
| C1 | 1 — promise holds worst case, **whole matrix, both platforms**, ratio ≤ 1.0 − margin | PARTIAL | Post-implementation evidence is **6 macOS/arm64 cells**: worst peak/prediction 0.976, worst peak/budget 0.887 ≤ 0.909, so the margin criterion holds where measured. But no Linux cell was re-measured against the shipped constants (the Linux TSV is the pre-implementation fitting matrix), and 6 cells is not the matrix. No raw artifact for the 6 cells — the numbers exist only in PLAN.md prose. The notes head the table "macOS/arm64" but nowhere record the cross-platform gap as outstanding. |
| C2 | 2 — held-out cells ratio ≤ 1.0 | DONE | CALIBRATION's falsification section: 3 unfitted cells × 2 platforms, 5 of 6 within 5%; the loose one (Linux c16, −22.1%) **over**-predicts, i.e. safe, and the looseness is stated. Post-implementation, 1 GiB/c8 and 2 GiB/c16 also lie outside the fitting range. |
| C3 | 3 — lever A on glibc (`1G/c4`, `2G/c4`) | DEVIATED (documented) | Never measured. PLAN "Remaining risks" (b) names it explicitly: "Lever A is macOS-only so far; glibc confirmation is validation row 3". Phase 1 lever; the Linux calibration ran on the Phase-1 binary but never as a with/without comparison. |
| C4 | 4 — no compression regression: distinct-read fixture **and** a low-duplication cell, with and without `--fastqc` | PARTIAL | (a) Low-duplication cell **missing**, acknowledged in CALIBRATION Limitations. (b) FastQC comparison: the only figure given is +0.9% output vs the non-FastQC run (improved from +2.9%) — an increase. (c) "Output no larger than baseline" is contradicted by the harness's own `out_bytes` column, which is monotone in `B` on this fixture (c4: 56.9 MiB → 507.29 MB, 37.9 → 518.84, 27.4 → 528.97, 19.0 → 540.10, 12.8 → 552.07). Interpolating, the 2 GiB/c4 shrink 56.89 → 42.00 MiB costs roughly **+1.5–2%** output, and the −4.4% default-cell shrink roughly +0.3%. The evidence is in the committed TSVs; no document states it. |
| C5 | 5 — output byte-identical to baseline | **MISSING** | No Phase 2 evidence. Phase 1 could claim byte-identity because it changed no sizing; Phase 2 changes `B` in every cell, so bin grouping — and therefore compressed bytes — must change, as the TSV's `out_bytes` column demonstrates. Nothing in the notes records that row 5 no longer applies, and no substitute check (same records, reordered) is reported. |
| C6 | 6 — throughput, plain **and** clumpy paths | **MISSING** | No wall/CPU numbers for the shipped constants, and none acknowledged as outstanding. The plain path is untouched by this diff (the `parallel.rs` change is a test-helper signature only), so only the clumpy half is unevidenced — bin sizes move by up to −32%, which changes flush frequency. |
| C7 | 7 — `--clumpify --fastqc`, 2 budgets × 2 core counts; ratio ≤ 1.0; **floor no higher than today** | PARTIAL | Ratios: three FastQC cells (1G/c4, 1G/c8, 2G/c16) at 0.813 / 0.679 / 0.750 of prediction, 0.739 / 0.617 / 0.682 of budget — all pass; the grid is 3 cells over 2 budgets × 3 core counts, not a filled 2×2 (no 2G/c4 or 2G/c8 FastQC cell). **The floor clause fails on the FastQC path at high cores**: c16 goes 620 → **840.4 MiB** (+35.6%), with the crossover between 8 and 12 cores. Non-FastQC floors all drop (c2 535.5 → 280.5, c4 539.0 → 289.3, c8 566.0 → 332.2, c16 620.0 → 418.0 MiB). The rise **is** acknowledged in the Deviation paragraph, which quotes the prior floor as 536 MiB — that is the `--cores 2` value; the correct c16 comparison is 620 MiB. |
| C8 | 8 — single-end series | DEVIATED (documented) | Not measured. CALIBRATION Limitations: "Single-end and `--clump_only` untested here; both share the layout." `calibrate.py` hardcodes `--paired`. |
| C9 | 9 — `--clump_only` at the default budget with `--cores 2` resolves, does not abort | DONE | Floor at c2 = 280.5 MiB against the 1 GiB default; both `integration_clump_only*.rs` suites pass. Verified by arithmetic and by the suite, not by an RSS measurement (acknowledged in Limitations). |
| C10 | 10 — multi-pair run, ratio ≤ 1.0 | **MISSING** | Every measured cell is a single pair; `calibrate.py` accepts one R1/R2. The Self-Review lists "multi-pair per-process peak" as an itemised edge case, but neither the Phase 2 notes nor CALIBRATION's Limitations record row 10 as unverified. |
| C11 | 11 — read length 51 bp **and** 150 bp; σ within the assumed range | DEVIATED (documented) | Only 51 bp measured; acknowledged in CALIBRATION Limitations ("One fixture, one read length"). Note the shipped σ = 1.4375 sits **below** the plan's own prior 51 bp figure (~1.46) and above its 150 bp figure (~1.21); nothing reconciles the new fit against those numbers. |

### D. Phase 2 fixture and matrix obligations

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| D1 | Rebuild the fixture on distinct reads; log distinct counts | PLAN "Phase 2 fixture and matrix" | DONE | `synth_bclconvert` RRBS PE, 10 M pairs, 51 bp, 1,048,545 distinct sequences per 2 M reads (~52%), byte-identical on both hosts (md5s recorded). Harness prints reads / distinct / mean length before measuring. |
| D2 | Identifying cells `cores ∈ {1,2,3,4}` | PLAN | DEVIATED (documented) | Ran 2–4; `--clumpify` requires `--cores ≥ 2`, and CALIBRATION states this ("2–4 is the whole identifying range"). |
| D3 | Falsification cells `cores ∈ {8,16}`, **2 budgets each** | PLAN | PARTIAL | Three cells, not four: c8 at B = 19.0 and 37.9, c16 at B = 19.0 only. The shortfall is not acknowledged. |
| D4 | Test-only override setting `n_bins` independently of `cores` | PLAN ("Better still") | DEVIATED (optional) | Not implemented. Phrased as a recommendation, not an obligation; not blocking. |
| D5 | Read length as an explicit axis; `--fastqc` series in the calibration | PLAN | DEVIATED (documented) | Neither is in the matrix; both are named in CALIBRATION Limitations. The FastQC series was instead measured at implementation time — see A4 for the missing artifact. |
| D6 | CHANGELOG rewritten to claim the budget is now kept | working-tree diff | DEVIATED (undocumented) | The diff reverses the Phase-1 post-review wording ("not yet a guaranteed bound") to "keeps peak memory inside the `--memory` budget it prints, at every core count and budget measured". The claim is scoped by "measured", but the measurement behind it is macOS-only (C1), and the Phase 2 notes' "Built" list does not mention the CHANGELOG edit at all. |

## Gaps (detail)

### A9: the falsifiable prediction was never evaluated

**Expected:** PLAN "Model form" §1 states that lever A takes `k` from 3σ to 2σ and calls it "a
falsifiable prediction the calibration must confirm (expect `k` to fall by roughly σ, i.e.
~1.3–1.4)".
**Found:** CALIBRATION fits `k` = 3.99 (macOS) / 3.71 (Linux) on the post-lever binary and
compares it only against the *shipped* assumption of 1.75. Nothing compares it against 2σ
(≈ 2.8), and no depth-2 measurement exists on the rebuilt fixture, so the "falls by ~1.3–1.4"
claim has no test. The code's docstring substitutes a different derivation (three resident
batch-equivalents per worker).
**Gap:** either evaluate the prediction (a depth-2 `k` on the same fixture) or record in the plan
that it was superseded, and by what.

### A5: the margin does not exceed the largest observed spread at the default budget

**Expected:** "Margin must exceed the largest observed spread" (PLAN "Staging" Phase 2).
CALIBRATION: "a 10% margin — the margin exceeds the largest within-cell spread on either host".
**Found:** the margin is proportional (`budget × 1/11`) while the spread is absolute (96.2 MiB
macOS). At 1 GiB the margin is 93.1 MiB; the criterion holds only for budgets ≥ ~1.06 GiB. At
2 GiB it is 186.2 MiB and comfortable.
**Gap:** either an absolute floor under the margin, or an explicit statement that the criterion is
met only above ~1.06 GiB. Measured headroom at 1 GiB/c4 (peak/budget 0.862) is larger than the
margin, so this is a stated-criterion gap rather than an observed breach.

### A3: FastQC as an addend rather than a `max`

**Expected:** both PLAN "Model form" §2 and CALIB item 1 require FastQC in a `max`, on the stated
mechanism that it runs after trimming.
**Found:** `static_bytes_for()` adds 24 MiB × cores on top of the 224 MiB trim reservation.
**Gap:** the departure is not recorded. Its consequence is the C7 floor regression (620 → 840 MiB
at c16 + FastQC), which in turn is what makes A11's deferred `--clump_only` fallback reachable —
under a `max` form the floor would not rise. The notes document the second-order effect (the
abort case) but not the first-order decision that produced it.

### C5 / C6 / C10: three validation rows with neither evidence nor acknowledgement

**Expected:** byte-identical output (row 5), throughput on both paths (row 6), and a multi-pair
cell (row 10).
**Found:** nothing for any of the three in the Phase 2 notes, and none of them appear in
CALIBRATION's Limitations or in the plan's Remaining risks. Row 5 is additionally *unsatisfiable*
as written, because Phase 2 moves `B` in every cell.
**Gap:** rows 6 and 10 need measurement or an explicit deferral; row 5 needs restating (e.g. as a
record-level equivalence check) or retiring for Phase 2 with the reason recorded.

### C4: the compression evidence in the repo points the other way

**Expected:** "output no larger than baseline", on the distinct-read fixture **and** a
low-duplication cell, with and without `--fastqc`.
**Found:** the low-duplication cell is missing (acknowledged). The `out_bytes` column of both
committed TSVs is monotone decreasing in `B`, so the Phase 2 shrink at large budgets increases
output size by roughly 1.5–2% at 2 GiB/c4 and ~0.3% at the default cell. The only figure quoted
anywhere is the +0.9% FastQC-vs-non-FastQC comparison.
**Gap:** state the output-size cost of the new bin sizes at 2 GiB and above, and run the
low-duplication cell — CALIBRATION already flags that bin size moves the clumpify crossover.

### A4: the FastQC constant has no raw data behind it

**Expected:** CALIB item 1 — "FastQC's cost under mimalloc is still unmeasured … must be measured,
not assumed."
**Found:** it was measured, and the notes give the numbers (76 / 139 / 283 MiB at c4 / c8 / c16;
191.9 MiB at c2), but `phase2/` contains no TSV for them and `calibrate.py` cannot produce them —
it has no `--fastqc` path and no core-count-only sweep.
**Gap:** commit the measurement (a TSV, or the harness flag that produced it), or the 24 MiB/core
constant cannot be re-derived by anyone else.

### B5: stale comment naming a deleted constant

**Expected:** PLAN "Tests" names `small_memory_budget()` and its 13 call sites as affected.
**Found:** untouched; the comment still reads "resolve_layout requires STATIC_OVERHEAD (512 MiB)"
at `src/clump_only.rs:1188`. The value (1 GiB) is still correct against every new floor.
**Gap:** the comment; no behavioural gap.

### D3: one falsification cell short

**Expected:** `cores ∈ {8, 16}`, 2 budgets each — 4 cells.
**Found:** 3 cells (c16 has only B = 19.0), on both platforms.
**Gap:** the second c16 budget, or a note that it was dropped.

## Test verification

| Test | File | Status |
|---|---|---|
| `resolve_layout_default` | `src/clump.rs:484` | PASS (rewritten: 4 GiB/c2 → ~113 MiB) |
| `resolve_layout_scales_with_cores` | `src/clump.rs:494` | PASS (rewritten: 4 GiB/c8 → ~45 MiB) |
| `resolve_layout_pins_the_default_cell` | `src/clump.rs:503` | PASS (pin updated to 19,006,356 B; recomputed independently and it matches the shipped constants) |
| `resolve_layout_fastqc_reserves_more_and_shrinks_bins` | `src/clump.rs:512` | PASS (new) |
| `resolve_layout_predicted_peak_fits_budget` | `src/clump.rs:534` | PASS — still tautological (formula vs its own inverse, now with fixed 9% slack); retirement deferred to Phase 3 |
| `resolve_layout_bails_when_too_small` | `src/clump.rs:560` | PASS (rewritten: 256 MiB/c8 against a 332 MiB floor); retains the `"static overhead"` assertion |
| `resolve_layout_bails_below_static_overhead` | `src/clump.rs:571` | PASS (rewritten: 128 MiB/c4) |
| `clumpify_floor_matches_the_smallest_resolvable_budget` | `src/clump.rs:579` | PASS (new; 8 cells × both directions) |
| `small_clump_layout()` helper (E0063 site) | `src/parallel.rs:1615` | PASS — updated for the 4th field |
| `tests/integration_clump_only.rs` (16 tests) | — | PASS unchanged |
| `tests/integration_clump_only_ubam.rs` | — | PASS unchanged |
| `tests/integration_preflight_tripwire.rs` banner guard | — | PASS unchanged |
| Whole suite | — | **732 passed, 0 failed** |

**Not present:** any integration-level test for the new FastQC floor or for the fallback path at a
budget between the old and new floors (e.g. `--clumpify --fastqc --cores 16 --memory 800M`, which
now warns and falls back). The plan does not require one; noting it because it is the only
user-visible behaviour change Phase 2 makes that no test pins.

## Adjacent observation (not a plan obligation)

`main.rs:192` prints the fallback floor with floor division (`min_required / (1024*1024)`) while
`clump.rs:250` bails with `.div_ceil(1024*1024)`. At c16 + FastQC the true floor is 840.4 MiB, so
the warning advises "need ≥ 840 MiB" — a value that itself falls back. `clumpify_floor_matches_the_
smallest_resolvable_budget` guarantees the floor in bytes, not the MiB figure the user is shown.

## Verdict

**INCOMPLETE — 12 unresolved (8 PARTIAL, 4 MISSING), plus 4 undocumented deviations.**

The calibrated model itself is implemented as specified and is internally consistent: the
constants are at or above the fitted worst-case values on both platforms, the prediction is the
exact inverse of the sizing, the FastQC term is threaded through every call site including the
`--fastqc_args`-alone path, the eight affected `clump.rs` tests are updated, and the suite, fmt
and clippy are clean at 732 passing. What is unresolved is evidence and disclosure, not arithmetic.

Must be addressed:

1. **A9** — evaluate the plan's own falsifiable prediction (`k` = 2σ) or record that it was
   superseded and by what derivation.
2. **A5** — the margin (93.1 MiB at 1 GiB) does not exceed the largest observed spread (96.2 MiB)
   at the default budget; either add an absolute floor or qualify the claim in both CALIBRATION
   and the `MARGIN_NUM` comment.
3. **A3** — record the decision to charge FastQC as an addend rather than a `max`, since it is the
   cause of the C7 floor regression and of A11's now-reachable abort.
4. **C1** — no Linux cell was measured against the shipped constants; the CHANGELOG's new promise
   (D6) rests on six macOS cells, and no raw artifact exists for them.
5. **C5, C6, C10** — byte-identity, throughput and multi-pair have neither evidence nor an
   acknowledged deferral. Row 5 needs restating, as Phase 2 makes it unsatisfiable as written.
6. **C4** — the committed TSVs show output size rising as `B` falls; quantify the cost at 2 GiB
   and run the still-missing low-duplication cell.
7. **A4** — commit the FastQC measurement or the harness flag that produced it.
8. **A10, B4, D3, B5** — small factual corrections: the default cell moved −4.4% (not −0.6%),
   eight tests were touched (not nine), one falsification cell is missing, and
   `small_memory_budget()`'s comment still names the deleted 512 MiB constant.

Gaps the plan **does** own, and which are therefore not counted against it: `--clump_only`'s
fallback (A11), lever A on glibc (C3), single-end (C8), read length (C11), the `cores ∈ {2,3,4}`
identifying range (D2), the tautological test's Phase 3 retirement (B8), and all Phase 3 docs
arithmetic (A12).
