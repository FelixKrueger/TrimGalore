# Plan — honest `--clumpify` memory accounting (#439), v3

**Issue:** [#439](https://github.com/FelixKrueger/TrimGalore/issues/439) (ratios corrected 2026-08-16)
**Investigation:** `plans/08152026_clumpify-memory-overshoot/FINDINGS.md`
**Reviews:** v1 → `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`; v2 → `PLAN_REVIEW_A_v2.md`, `PLAN_REVIEW_B_v2.md` (all four: "do not implement as written")
**Spike:** `spikes/SPIKE_memory-levers.md` — the measurement this revision is built on
**Baseline:** dev tip `2ab0918`; PR #437 (mimalloc) merged, no prerequisites outstanding

## Revision history

- **v1 — refuted.** Split the static reservation by FastQC and set the per-core
  coefficient to 3.0, claiming 48% larger bins. The gain came from spending the
  cushion that was absorbing the under-count; also would not have compiled.
- **v2 — refuted.** Staged the measurement and pinned bin size at the default.
  Both reviewers independently proved the **Goal was unsatisfiable**: at `1G/c4`
  actual RSS (1050 MiB) already exceeded the budget (1024 MiB) with `B` held, so
  no formula could bound it. v2 also rested on a fitting artefact (a "missing
  per-core fixed term") and on a fixture that is 30,000 distinct reads cycled
  227×.
- **v3 — this document.** The spike found the escape: cutting the work-channel
  depth removes 12.6–13.9% of actual RSS at zero throughput cost, which makes the
  promise satisfiable at the default cell. Structural reduction comes first;
  calibration follows, against a much smaller residual, on a rebuilt fixture.

## Goal

Make the printed predicted peak a true upper bound on actual peak RSS across the
`(cores, --memory)` space on both platforms, by **reducing actual memory first**
and calibrating the model against what remains.

The escape from v2's contradiction is named and measured: lever A below. This is
explicitly *not* a pure-accounting change — v2 proved accounting alone cannot
satisfy the promise at the shipped default.

## Evidence

### Established by the spike (macOS/arm64, `2ab0918`, 2 reps, max-over-reps ratio)

| cell | base | A: work channel 2→1 | B: reserve from record bytes | A+B |
|---|---|---|---|---|
| `1G/c4` | **1.045** | 0.913 | 0.999 | **0.857** |
| `2G/c4` | 1.236 | 1.064 | 1.235 | **1.056** |
| `1G/c8` | 0.941 | 0.850 | 0.881 | 0.796 |

- Output is **byte-identical** to baseline under both levers — these are pure
  memory/scheduling changes.
- Lever A costs **no wall clock** (−3.96% at `2G/c4`, i.e. faster). v2 deferred it
  on a throughput fear that measured as zero.
- Lever B's benefit **does not scale with bin size** (−47 MiB at 1G, −3.6 MiB at
  2G), so it is a correctness fix, not a memory lever, until explained.
- **Large budgets still breach**: `2G/c4` = 1.056 with both levers. The residual
  is ~6%, where it was ~24%.

### Invalidated, and why it matters

- **Every coefficient measured so far is provisional twice over**: measured
  *without* the levers, and on a fixture of 30,000 distinct reads cycled 227×.
  That fixture inflates compression-dependent quantities, which include the
  output-buffer share of RSS. It invalidates v2's measurement 6 and equally
  reviewer B's pooled fit (`f0` = 244 MiB, σ = 1.33, `k` = 4.61, rmse 4%).
- **There is no "missing per-core fixed term."** The 151 → 302 MiB intercept shift
  is slope↔intercept trade-off within each series; a 3-parameter model with no
  such term fits all 19 macOS cells including the `c8` cell v2 called
  contradictory. v2's measurement 2 is deleted.
- **Precision was overstated**: `B_max` is 18.8 **± 1.3** MiB, not "within 1%",
  and `k`'s 95% CI is [1.96, 3.75] — not resolved by the two-series difference.
- **`1G/c4` is not "the default cell"** — the default is `--cores 1`-ish usage
  across a range of budgets. Pinning a constraint at one cell protects the wrong
  thing, and is how #437 came to fix one cell and not the mechanism.

## Model form

Adopt the mechanism-derived form rather than a curve fit (reviewer B):

```
peak ≈ f0(fastqc) + B · (σ · n_bins + k · cores)
   σ ≈ spine factor (text + Vec spine), read-length dependent
   k ≈ (work-channel depth + 1) · σ
```

`FINDINGS.md` already derives and arithmetically confirms this on Linux (695
predicted vs 672 measured). Two consequences:

1. **Lever A changes `k` by construction**: depth 2 → 1 takes `k` from `3σ` to
   `2σ`, a falsifiable prediction the calibration must confirm (expect `k` to fall
   by roughly σ, i.e. ~1.3–1.4).
2. **FastQC belongs in a `max`, not an addend** — it runs sequentially, after
   trimming, so the peak is `max(trim_phase, fastqc_phase)`. Its cost is also
   currently unmeasured under mimalloc (the ~425 MB figure is stock-allocator and
   uncorrected for the 4.86% unit error).

## Staging

**Phase 1 — structural reduction (code, small, measurable).** Levers A and B,
gated on clumpy mode. This is the part that makes the promise reachable.

**Phase 2 — calibrate the residual.** Fit `f0`, σ, `k` as an **upper envelope with
an explicit named margin**, not a least-squares mean — "ratio ≤ 1.0 in every cell"
is unachievable from a mean fit whose within-cell spread is already 4.9%. Margin
must exceed the largest observed spread; state which quantile the banner means.

**Phase 3 — defend and correct.** CI check, docs arithmetic, banner unit label.

### Phase 2 fixture and matrix

- **Rebuild the fixture on distinct reads.** Use oxy's `synth_bclconvert` pair
  (236 MB each, simulated but genuinely diverse), streamed to macOS so both hosts
  share bytes. Record `sort -u` distinct-ID and distinct-sequence counts in the
  harness output from now on — the artefact was invisible because nothing logged
  sample size.
- **Identifying cells:** `n_bins = max(16, 4·cores)`, so σ and `k` are separately
  identifiable only where `n_bins` is pinned — `cores ∈ {1, 2, 3, 4}`. Sweep `B`
  across all four.
- **Falsification cells, not fitting cells:** `cores ∈ {8, 16}`, 2 budgets each.
- **Better still:** add a test-only override that sets `n_bins` independently of
  `cores`, turning the sweep into a designed experiment. Cheaper than the runs it
  saves.
- **Bounded by the host:** some corners (high cores × largest `B`) exceed 16 GiB
  on the macOS box; cap the matrix accordingly rather than discovering it midway.
- Include read length as an explicit axis (σ measures ~1.21 at 150 bp, ~1.39 at
  65 bp, ~1.46 at 51 bp), and a `--fastqc` series under mimalloc.

## Implementation outline (Phase 1)

1. **`src/parallel.rs` — lever A.** Work-channel depth 2 → 1, **gated on clumpy
   mode** so plain round-robin keeps its jitter absorption: the two construction
   sites are shared by both paths. Re-derive line numbers at implementation time
   (`#434`/`#436` have already shifted this file once).
2. **`src/parallel.rs` — lever B.** Size the per-batch output reserve from actual
   record bytes via `estimated_record_bytes` (already imported) rather than
   `len * 300`. Note the paired site applies `cap` to **both** mate buffers, so
   size R2 from R2 rather than reusing R1's sum.
3. **`ClumpLayout` constructibility.** Both fields are already `pub`; the hazard
   is **E0063 (missing field)**, which visibility does not fix — any new field
   requires updating the struct literal in `parallel.rs`'s test helper or adding a
   constructor. v1 and v2 both got this remedy wrong.
4. **Drop the drift surface.** Store what the layout was sized against and remove
   the `cores` parameter from `predicted_peak_bytes`, so no caller can evaluate
   the prediction against a different core count than the sizing used.
5. **`--clump_only` needs a fallback decision.** It bypasses
   `main.rs::resolve_clump_layout` entirely and propagates `resolve_layout`'s
   error with a bare `?`, so any raised floor turns it from "works" into
   "aborts" — including ~23 integration invocations that run at the default
   1 GiB with `--cores 2`. Decide explicitly: same warn-and-continue fallback, or
   documented hard failure. It is also single-threaded while sizing bins from the
   user's `--cores`, so decide what `cores` it passes.
6. **Banner unit label.** It prints MiB labelled `MB`. Fix it — this is the slip
   that produced a 4.86% error in my own harness, and it sits in the exact line
   #439 is about.

Phase 2 then replaces the denominator with the calibrated form; no coefficient is
written into this plan, deliberately.

## Tests

Affected existing tests — **five** in `src/clump.rs`, not six:
`resolve_layout_default`, `resolve_layout_scales_with_cores`,
`resolve_layout_bails_when_too_small`, `resolve_layout_bails_below_static_overhead`
(**panics** rather than fails if the reservation drops), and
`resolve_layout_predicted_peak_fits_budget` (needs any new argument). Plus a
`"static overhead"` substring assertion that a bail-message rewrite breaks.

Also affected, and missed by v2: `clump_only.rs`'s `small_memory_budget()` helper
and its 13 call sites, both `tests/integration_clump_only*.rs` files (default-budget
invocations — see item 5), and the banner guard in
`tests/integration_preflight_tripwire.rs`.

`resolve_layout_predicted_peak_fits_budget` **passes for any coefficient** — it
compares the formula against its own inverse. Retire it or move the assertion
somewhere it bites; Phase 3's CI check is the real replacement.

New: a banner-present assertion in every measured cell (a silent fallback to plain
mode reads as a ratio of ~0.1), and a regression test pinning resolved `B` at a
named cell so a future coefficient edit cannot move bin size unnoticed.

## Validation

| # | Verify | How | Expected |
|---|---|---|---|
| 1 | Promise holds, worst case | **max over reps** per cell, whole matrix, both platforms | ratio ≤ 1.0 − margin in every cell |
| 2 | Held-out cells | 3 cells excluded from the fit | ratio ≤ 1.0 |
| 3 | Lever A on glibc | oxy, `1G/c4` and `2G/c4` | RSS reduction confirmed; glibc retains less, so the benefit may differ |
| 4 | No compression regression | distinct-read fixture **and** a low-duplication cell, with and without `--fastqc` | output no larger than baseline |
| 5 | Output unchanged | byte-compare `_val_*.fq.gz` against baseline | identical |
| 6 | Throughput | wall/CPU, plain **and** clumpy paths, reps beating the observed noise bar | plain path untouched; clumpy within noise |
| 7 | FastQC path | `--clumpify --fastqc`, mimalloc, 2 budgets × 2 core counts | ratio ≤ 1.0; floor no higher than today |
| 8 | Single-end | SE series (lever A patches the SE channel too) | ratio ≤ 1.0 |
| 9 | `--clump_only` | default-budget invocation at `--cores 2` | resolves, does not abort |
| 10 | Multi-pair | the promise is per-process; every measurement so far is per-pair | ratio ≤ 1.0 across a multi-pair run |
| 11 | Read length | 51 bp and 150 bp cells | σ within the assumed range |

Baselines must be captured **before** Phase 1 lands, or rows 4–6 are not runnable.

## Assumptions

1. **Fixed:** `--memory` bounds peak RSS of the process, not virtual size.
2. **Fixed:** `n_bins = max(16, 4·cores)` is unchanged — altering it would change
   output grouping and therefore compression.
3. **Configurable:** `--memory` default stays `1G`; `--compression` is irrelevant
   to RSS (measured, both levers and baseline).
4. **Open:** σ is read-length dependent and the calibration must span more than
   one length, or the model is tuned to one fixture.
5. **Not assumed:** any coefficient. v1 and v2 both failed by naming one.

## Questions or ambiguities

- **Open — does the residual need calibration or another structural cut?** `2G/c4`
  sits at 1.056 with both levers. Charging more per bin (smaller bins at large
  budgets) and buffer recycling are the two candidates; compare them on data
  before choosing.
- **Open — `--clump_only`'s fallback and its `cores`** (item 5).
- **Resolved — why lever B does not scale** with bin size. Reviewer A measured it
  on a 20,000 × 10 kb-read fixture at `--cores 4`: the plain-text branch peaked
  802 MiB and the gzip branch 508 MiB from an identical 82 MB reserve, the
  difference being exactly the bytes-written difference. Large over-reserves are
  mmapped and never faulted in; small ones come out of an already-committed segment
  and are resident. Hence −47 MiB at 1G and −3.6 MiB at 2G. B is tidiness.
- **Open — docs savings table.** The formula, floor range, worked values and five
  unit labels are pure arithmetic and can be fixed now. The published savings
  figures need real libraries and will go stale if Phase 2 moves bin size.
- **Resolved by spike:** lever A is primary and free; recycling is demoted;
  byte-identity holds.
- **Resolved by review:** no per-core fixed term; the reassembly `BTreeMap` needs
  no cap (round-robin plus channel depth bounds the reader's lead).

## Implementation notes — Phase 1 (2026-08-16)

Implemented on `dev` at `2ab0918`. `git diff --stat`: `CHANGELOG.md`, `src/clump.rs`,
`src/clump_only.rs`, `src/main.rs`, `src/parallel.rs`,
`tests/integration_preflight_tripwire.rs` — 49 insertions, 24 deletions.

**Built:**

1. **Lever A, gated on clumpy mode.** `queue_depth = if clump_layout.is_some() { 1 }
   else { 2 }`, applied at both work-channel construction sites (paired and
   single-end). Plain round-robin keeps depth 2 and its jitter absorption.
2. **Lever B, ungated — a correctness fix, not a memory lever.** Both output
   buffers are reserved from their own mate's `estimated_record_bytes` sum instead
   of a shared `len * 300`. Sizing R2 from R2 fixes the review's note that reusing
   R1's count under-reserves unequal mates. Ungated because sizing from the batch's
   own bytes is right in every mode.

   Both code reviewers established that its *memory* effect is data-dependent and
   inverts: mean `estimated_record_bytes` is 337.9 B on a 150 bp fixture against
   the old 300 (a **larger** reserve), and 192.6 B on 75 bp BS-seq (36% smaller),
   with break-even near 125 bp. Since every clumpy run is a gzip run and gzip
   buffers hold compressed bytes, lever B is close to memory-neutral in clumpy
   mode. Claim it as tidiness — exact bound in the plain branch, no doubling
   cascade on long reads — and not as part of the memory fix.
3. **Drift surface removed.** `ClumpLayout` now stores the `cores` it was sized
   for and `predicted_peak_bytes()` takes no argument, so the banner cannot be
   evaluated against a different core count than the sizing used. Added
   `ClumpLayout::new()`; the private field makes the struct literal in
   `parallel.rs`'s test helper impossible, which is the E0063 hazard both earlier
   reviews flagged — resolved by the constructor rather than by visibility.
4. **Five unit labels corrected** MB → MiB (one in `main.rs`, four in
   `clump_only.rs`); all five divide by 1024². Guarded by a new assertion in
   `se_trim_clumpify` that would have failed before the change.

**Measured (shipped gated binary vs `2ab0918`, diverse fixture, max-over-reps):**

| cell | base | Phase 1 | ΔRSS |
|---|---|---|---|
| clumpy `1G/c4` | 1033.6 MiB, ratio 1.008 | 884.1 MiB, **ratio 0.863** | −14.5% |
| clumpy `2G/c4` | 2403.2 MiB, ratio 1.173 | 2077.2 MiB, **ratio 1.014** | −13.6% |
| clumpy `1G/c8` | 982.9 MiB, ratio 0.959 | 817.6 MiB, **ratio 0.798** | −16.8% |
| plain `c4` | 96.8 MiB | 91.2 MiB | −5.8% |

Output is **byte-identical** to `2ab0918` on both the clumpy and plain paths
(34,629,640 and 407,440,179 bytes respectively). `cargo test --release`: 729
passed, 0 failed. `cargo fmt --check` clean, `cargo clippy --all-targets
--release -- -D warnings` clean.

**Iteration log:**

- `#1` Gated lever A on `clump_layout.is_some()` rather than patching the depth
  unconditionally as the spike did, so plain mode keeps its queue. Clumpy results
  reproduced the spike (0.863 vs 0.857 at `1G/c4`); plain mode kept depth 2.
- `#2` First plain-mode timing showed wall 54.2 → 74.9 s, a 38% apparent
  regression. Re-measured interleaved with 3 reps: CPU +0.47%, wall mean +5.94%
  with a 10.1 s within-binary spread against base's 1.9 s. The outlier rep's CPU
  was normal, so the extra wall time was waiting on a 407 MB output flush, not
  compute. No regression; CPU is the measure to quote on this box.
- `#3` Verification harness mis-parsed the baseline's banner (it prints `MB`,
  Phase 1 prints `MiB`). Harness bug, not an implementation one; regex now accepts
  both.

**Deviation from the plan:** item 5 (`--clump_only`'s missing fallback) is
**deferred to Phase 2**, not implemented. The abort path only becomes reachable
when the floor moves, and Phase 1 changes no sizing arithmetic —
`clumpify_min_memory_bytes` is untouched, so every budget that resolved before
still resolves. The decision it needs (warn-and-continue vs documented hard
failure) is a behaviour choice for Felix, and making it now would be guessing
ahead of the floor change that motivates it.

**Not attempted in Phase 1** (unchanged from the plan's staging): the calibrated
model, the fixture rebuild, the docs arithmetic, and the CI check. `2G/c4` still
sits at 1.014 — inside the budget but with ~1% of margin, which is not a margin.
Phase 2 remains necessary.

### Post-review fixes (2026-08-16)

Dual code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) and the coverage audit
(`COVERAGE.md`, verdict INCOMPLETE — 3 items) produced no correctness findings;
both reviewers independently confirmed the gate is sound, part 3 behaviour-neutral,
and fmt/clippy/tests clean. Applied:

1. **CHANGELOG rewritten.** It claimed "keeps its peak memory inside the `--memory`
   budget it prints" — a bound `2G/c4` = 1.014 contradicts. Both reviewers rated
   this High; A's failure scenario was a user reading the entry, running
   `--memory 2G --cores 4` under a 2 GiB allocation and being OOM-killed. Now
   claims the measured 13–17% reduction and the MiB labels, and says plainly that
   the budget is not yet a guaranteed bound at large values.
2. **The `resolve_layout` model docstring corrected.** It still described in-flight
   work as `cores × B`; lever A makes it `2 × cores × B`, and the docstring now
   also states that the formula charges one batch per worker where two can be
   resident — i.e. the under-estimate is explicit. This matters because Phase 2
   calibrates `k` from that derivation.
3. **Lever B re-described** as a correctness fix (above).
4. **Banner guard now pins both labels.** The original assertion's two clauses both
   keyed on the *bins* label, so reverting the peak label alone would have passed.
   Now `contains("MiB; predicted peak ≈") && contains("MiB (gzip level")`.
5. **Comment style.** Dropped "rather than passed in again" from the `cores` field
   doc (it referenced the removed parameter), deleted the single-end reserve comment
   (it narrated its own line; the paired twin earns its keep because *per-mate* is
   the non-obvious part), and restated the test comment as fact rather than
   evidence.
6. **`resolve_layout_pins_the_default_cell` added** — asserts 1 GiB/4 cores resolves
   to exactly **19,006,356** bytes, so no future coefficient edit can move the cell
   #439 is about without breaking a test. 730 tests pass. (This note read 19,884,107
   until 2026-08-17; that was the pre-Phase-2 value and the shipped test never
   asserted it.)

**Coverage gaps resolved:** the `B`-pinning test is item 6 above. The tautological
`resolve_layout_predicted_peak_fits_budget` was **deliberately left in place and
assigned to Phase 3**, since the plan pairs its retirement with "Phase 3's CI check
is the real replacement" — deleting a passing test before its replacement exists
would reduce coverage, not improve it. **Discharged 2026-08-17** in Phase 3's PR 3,
which is the PR that adds the replacement, so the precondition held.

**Logged, not fixed** (all pre-existing or out of Phase 1 scope):

- `docs/performance/clumpy.md` quotes the banner as "32 bins × 12 MB", now doubly
  wrong: the figure was already incorrect against the shipped formula, and the unit
  label has changed. Phase 3.
- `docs/performance/threading.md` claims memory "never exceeds ~100 MB"; reviewer A
  measured **802 MiB** in plain mode at `--cores 4` on 10 kb reads, because
  `BATCH_SIZE` is 4096 *records* regardless of length. Pre-existing and unbounded
  in read length — worth its own issue, and sharper for ONT/PacBio users than #439.
- The reserve is exact in the plain branch **except under `--rename`**, which
  appends ~14 B/record; consequence is one realloc.
- `ClumpLayout::new` is public with no production caller and validates neither
  invariant `resolve_layout` enforces.
- Reviewer B found a pre-existing bug: the EOF bin flush swallows a failed send.
- **Next structural lever for the `2G/c4` residual** (reviewer A, measured): the
  result channel `sync_channel(cores * 2)` plus the `pending` map hold *written*,
  therefore resident, compressed buffers and are untouched by this diff. Cheaper to
  try than buffer recycling.

## Implementation notes — Phase 2 (2026-08-17)

Implemented on top of `1595d59`. Calibration inputs: `phase2/CALIBRATION.md`.

**Built:**

1. **Calibrated sizing.** `usable = budget × 10/11 − static`, and
   `bin_byte_budget = 16 × usable / (23 × n_bins + 64 × cores)` — σ = 23/16 and
   k = 4, the worse of the two platforms' fits. `predicted_peak_bytes()` is the
   exact inverse, so the banner cannot disagree with the sizing.
2. **`STATIC_OVERHEAD_BYTES` (flat 512 MiB) replaced** by a 224 MiB trim-phase
   reservation (measured max 205 MiB across both platforms) plus, when FastQC is
   requested, **24 MiB per core**.
3. **`resolve_layout` and `clumpify_min_memory_bytes` take `fastqc`**, threaded
   from `Cli::fastqc_requested()` in `main.rs` and from each of the four
   `clump_only.rs` call sites' own `fastqc_requested` parameter.
4. **Nine `clump.rs` tests** rewritten or added, including
   `clumpify_floor_matches_the_smallest_resolvable_budget`, which asserts the
   advertised floor resolves and a MiB below it does not — that test caught a real
   bug on first run (see iteration log).

**Measured, macOS/arm64, 10 M-pair distinct-read fixture, worst of 2 reps:**

| cell | `B` | peak / prediction | peak / budget |
|---|---|---|---|
| 1 GiB / c4 | 18 MiB | 0.949 | 0.862 |
| **2 GiB / c4** | 42 MiB | 0.976 | **0.887** (was 1.173 pre-Phase-1, 1.014 after) |
| 1 GiB / c8 | 9 MiB | 0.912 | 0.828 |
| 1 GiB / c4 `--fastqc` | 15 MiB | 0.813 | 0.739 |
| 1 GiB / c8 `--fastqc` | 6 MiB | 0.679 | 0.617 |
| 2 GiB / c16 `--fastqc` | 8 MiB | 0.750 | 0.682 |

Every cell is inside both the printed prediction and the budget. `cargo test
--release`: 732 passed. fmt and clippy clean.

**The FastQC term was measured, not assumed — and the first attempt was wrong.**
I shipped a flat 288 MiB addend, reasoning from a pre-#437 measurement where
FastQC appeared to add ~425 MB. Measuring the shipped binary showed that was an
allocator-retention artifact: at `--cores 2` with a minimal pool the whole process
peaks at 191.9 MiB, *below* the 224 MiB trim reservation. Probing across core
counts showed FastQC's non-pool cost is ~76 / 139 / 283 MiB at c4 / c8 / c16 — it
scales with `--cores`, because `fastqc::run` is handed `cli.cores`. Replacing the
flat addend with 24 MiB/core took the `1G/c4 --fastqc` bin budget from 10 → 15 MiB
and its output-size penalty from +2.9% to +0.9% against the non-FastQC run.

**Iteration log:**

- `#1` Test literal `288 * 1024 * 1024 * 16` overflowed `i32` at compile time;
  typed it `u64`.
- `#2` `clumpify_floor_matches_the_smallest_resolvable_budget` failed: the
  advertised floor did not resolve. Cause was my own margin arithmetic —
  `resolve_layout` computed `budget / 11 * 10` while the floor inverted
  `× 11 / 10`, so the two rounded apart. Multiplying first (`budget * 10 / 11`)
  makes them exact inverses.
- `#3` `800M/c16 --fastqc` fell back to plain mode during verification. Correct
  behaviour — the floor at that cell is ~841 MiB — but it proved the harness
  assertion works: a fallback would otherwise have been recorded as a beautiful
  ratio of 0.2.

### Post-review fixes — Phase 2 (2026-08-17)

Dual code review (`CODE_REVIEW_{A,B}_phase2.md`) and the coverage audit
(`COVERAGE_phase2.md`, INCOMPLETE — 12 unresolved, 4 undocumented deviations).
Both reviewers independently brute-forced the two load-bearing properties and
confirmed them: the predicted peak never exceeds the budget (swept to `u64::MAX`
× cores 1…2²⁰ × fastqc, 0 violations, 0 overflows), and the advertised floor is
the *least* resolvable budget. Applied:

1. **`main.rs`'s below-floor warning was advertising an unusable `--memory`.** It
   truncated to MiB while the sibling bail message — which this work had already
   switched to `div_ceil` — rounded up. B measured 124 of 126 cells giving advice
   that fails on retry. Now `div_ceil`, and
   `clumpify_floor_matches_the_smallest_resolvable_budget` asserts the **advised
   MiB figure** resolves, not just the byte-exact floor, so the two messages
   cannot drift apart again.
2. **`cores` split into three roles** (`LayoutInputs`): `bin_cores` drives
   `n_bins` (pinned, so bin grouping does not shift), `workers` counts batches in
   flight, and `fastqc_threads` is FastQC's real thread count from
   `fastqc::resolved_threads`. This fixes two findings at once — `--fastqc_args
   "-t 16"` used to reserve for `--cores` while FastQC ran 16 threads (the
   prediction stopped being a bound), and `--clump_only` was charged 64 units per
   worker for a synchronous mode, squeezing bins to 2.07 MiB at `--cores 16
   --fastqc`, a hair over `MIN_BIN_BYTES`.
3. **The FastQC term is capped at 16 charged threads**, so the floor cannot run
   away on high-core hosts: at `--cores 32 --fastqc` it is 1012 MiB rather than
   1434 MiB. Scaling above 16 threads is unmeasured and the constant says so.
4. **Four new/strengthened tests.** The two bail tests were inert — they passed
   under all 17 of B's constant perturbations, including reverting the static to
   512 MiB — and now sit one byte under the floor and assert the message names the
   charged reservation. Added `fastqc_reservation_scales_with_threads_and_is_capped`
   and `synchronous_paths_are_not_charged_for_absent_workers`, the per-core scaling
   and worker-split behaviours that previously had no test that could see them.

5. **Byte-identity measured, not assumed** (coverage row C5). Against the Phase 1
   binary on the same input, at `1G/--cores 4`: `.fq.gz` bytes **differ**, and the
   record multiset is **identical** on both mates. Output grew +0.27% / +0.24% —
   the compression cost of the default cell's bins moving 18.96 → 18.13 MiB. So
   Phase 1's "byte-identical output" claim correctly stops applying here, and the
   CHANGELOG now says so.

6. **`--fastqc_args "-t N"` charged, verified end to end.** At `--cores 2`,
   `--memory 1025M`: bins are 22 MiB plain, 21 MiB under `--fastqc` (2 threads
   charged), and 10 MiB under `--fastqc_args "-t 16"` (16 charged). Before the
   split the third case reserved for 2 threads while FastQC ran 16.

**Recorded deviation from the calibration spec.** `CALIBRATION.md` item 1
recommended treating FastQC as a `max`; this ships an **addend**, and the plan
failed to say so. The reviewers split on it: A held the spec was right, B held the
spec was wrong because the pool's pages stay resident across the
trim→FastQC boundary under mimalloc, so `max` would under-reserve. B's reading is
what the probes show — at `--cores 16` with only 64 MiB of pool the process peaked
at 438.8 MiB, i.e. the trim reservation *plus* ~215 MiB of FastQC. The addend
stands; the spec is the document that needs amending.

**Corrections to the Phase 2 notes above.**

- The default-cell pin moved 19,884,107 → 19,006,356 B = **−4.4%**, not the −0.6%
  `CALIBRATION.md` forecast. The forecast used the exact fitted values; the shipped
  constants round σ up to 23/16 and `f0` up to 224 MiB, and those roundings cost
  the difference.
- "Nine tests" was **eight** (six rewritten, two added); the sizing-related count
  is now twelve.
- The margin rule as written in this plan ("margin must exceed the largest observed
  spread") is **not met at the default budget**: the margin is proportional (1/11)
  while the spread is absolute (96.2 MiB), so at `--memory 1G` it is 93.09 MiB and
  only clears 96 MiB above ~1.03 GiB. No measured cell breaches, but the rule
  should read "≈9% of budget, which exceeds the observed spread above ~1 GiB".
- The floor comparison in the deviation note below compared c16's new floor against
  **c2's** old one. At `--cores 16` the old floor was 620 MiB, not 536 MiB.

**Deviation:** `--clump_only`'s missing warn-and-fall-back (Phase 1's deferred
item) is **still deferred**, and it now matters more: the FastQC floor rises with
cores (~841 MiB at `--cores 16 --fastqc` against 536 MiB before), and
`--clump_only` propagates a layout error with a bare `?` rather than degrading.
Its own integration tests run at the default 1 GiB with `--cores 2`, where the
floor is 281 MiB, so nothing breaks today — but a `--clump_only --fastqc
--cores 16 --memory 800M` invocation that used to work will now abort. That is a
behaviour decision for Felix, not a guess for me.

## Self-Review

**Logic.** Three revisions, each killed by evidence: v1 by review, v2 by its own
unsatisfiable Goal, and both by measurement. v3 inverts the order that failed
twice — reduce actual memory first, then fit what remains, and name no
coefficient. The one structural claim it does make (lever A takes `k` from `3σ` to
`2σ`) is written as a falsifiable prediction for Phase 2 to confirm.

**Efficiency.** Lever A measured free; lever B measured free. Phase 2's matrix is
sized by identifiability rather than by factorial completeness, which is what made
v2's 80 runs mostly redundant.

**Edge cases.** Carried and now itemised: `--fastqc_args` alone (use
`fastqc_requested()`), the FastQC-path floor, `--clump_only`'s abort path and
single-threadedness, single-end's one bin pool, multi-pair per-process peak,
read-length dependence, and silent fallback masquerading as a good ratio.

**Integration.** Docs arithmetic can land independently. Perl byte-identity is
unaffected — `--clumpify` has no Perl counterpart. Output byte-identity under both
levers is measured, not assumed.

**Remaining risks.** (a) Every coefficient measured to date is provisional until
the fixture is rebuilt on distinct reads — including the ones a reviewer fitted.
(b) Lever A is macOS-only so far; glibc confirmation is validation row 3.
(c) The promise is verified per-cell, never proven; cells outside the matrix can
still breach, which is why rows 2 and 10 exist.
