# Plan review — Phase 4 (σ refit, #457) — **Reviewer B**

**Target:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v1
**Repo state:** `dev` at `c7dd106` (not switched; verified with `git rev-parse`)
**Method:** code reading + independent arithmetic. **No `trim_galore` was run**; no fixture was
created; no tracked file was modified. Every claim below is tagged **code-verified**,
**arithmetic-verified** or **unverified**.

My arithmetic replication of `resolve_layout` / `clumpify_min_memory_bytes` reproduces **all six
committed pin values** (19,006,356 / 27,453,626 / 8,681,915 ×2 / 23,911,222 / 11,761,649) and **all
six published docs floors** (281 / 334 / 590 / 1012 / 933 / 1356 MiB) exactly, so the numbers I quote
below come from a formula that matches the shipped one.

---

## Verdict

The plan's diagnosis is right, its structural form is *almost* right, and its two headline arithmetic
tables (§7a, §7b) are both wrong in the same direction: **they scale quantities that contain the
224 MiB static term as if those quantities were pure bin pool.** That single error class produces
three separate mis-statements (§6's floor, §7b's band table, Step 10's docs delta) and, more
seriously, it combines with a 2.5×-too-small `s` to make the chosen option look ~3× cheaper than it
is.

The decisive number behind choosing Option C — *"a ~5–6% bin shrink is a much smaller price than §2.5
led me to expect, and it is the number that makes this option viable rather than merely simple"* — does
not survive contact with the arc's own published measurements. Phase 3 §0.4's macOS 25 bp runs
back-solve to **σ ≈ 1.72–1.91**, not 1.55, and the code says why: `PairedBin::push` `reserve_exact`s
**148 bytes per pair** of spine up front. The honest price of Option B/C at the CI fixture's own record
size is **13–20% of bin size**, not 5–6%.

At that price the option comparison inverts, so my primary recommendation is: **do not treat §6 as
decided.** Keep it as a recommendation pending Step 3, add the alternative in **A1** to the option set,
and fix §7a/§7b before either is quoted again.

What the plan gets right and should keep: the mechanism → form → fit chain (§2.1–2.3); catching that
the harness is broken (§2.6); insisting `a` be fitted rather than assumed 1; the counter-intuitive
lever direction (§2.5); the observation that the refit and the drift band are coupled (§7b's
*conclusion*, if not its table); that output bytes change a second time (§3.4); and choosing a warning
site downstream of `sanity_check_any` (§8/Step 7 — that ordering is correct, see A5).

---

## Critical

### C1 — `s ≈ 28.7 B/record` is ~2.5× too small, and the code states the right value

**code-verified + arithmetic-verified.**

`FastqRecord` is three `String`s (`fastq.rs:58-65`) = 72 B inline. `PairedBin::push`
(`parallel.rs:900-914`) and `SingleBin::push` (`parallel.rs:1286-1296`), and their twins in
`clump_only.rs:103-113` / `:147-159`, do this on the first push after every flush:

```rust
let predicted = self.budget.div_ceil(bytes_this).max(1);
self.r1.reserve_exact(predicted);  self.r2.reserve_exact(predicted);  self.keys.reserve_exact(predicted);
```

So each bin immediately reserves `B/(2r)` slots × (72 + 72 + 4) = **148 B per pair = 74 B per
record** — paired; **76 B per record** single-end (`Vec<FastqRecord>` + `Vec<u32>`). This is
allocated in full right after each flush, so unlike the text it does **not** scale with fill. The
mechanistic form is therefore

```
σ(r) = φ·ρ + 74/r          (φ = mean fill fraction, ρ = allocator rounding on the text)
```

which *is* `a + s/r` — the form is right — but with **s ≈ 74 plus per-record size-class rounding**,
not 28.7.

Two independent checks, both from data already in the arc:

1. **Back-solve Phase 3 §0.4's macOS 25 bp measurements** through the model's own form
   `peak = f0 + B(16σ + 4w)` with `f0 = 204.7` (CALIBRATION.md) and B from `resolve_layout`:

   | cell | measured peak | implied σ |
   |---|---:|---:|
   | 1G/c2 (min of 4 reps) | 1048.1 MiB | **1.81** |
   | 1G/c2 (max of 4 reps) | 1082.6 MiB | **1.91** |
   | 1G/c4 | 994.3–1000.0 MiB | **1.72** |

   With `s = 74` and macOS's fitted σ(51 bp) = 1.41, the prediction is σ(25 bp) = 1.41 − 74/147.9 +
   74/95.9 = **1.70** — inside the measured range. With the plan's `s = 28.7` it is 1.55, which
   **under-predicts those same runs by 50 MiB at c4 and 96–130 MiB at c2.**

2. `a` then comes out **below** 1.0 (≈0.86–0.93), not 1.238. §2.4's interpretation — *"the pool
   appears to hold ~24% more than one copy of each bin's bytes"* — is very likely sign-inverted: a
   bin flushes *at* B and is then empty, so mean occupancy is below one copy. §2.3's derivation
   (`pool per bin = B + s·(B/r)`) silently assumes every bin holds a full `B` of text at the peak, and
   the fill fraction it drops is exactly what `a` measures (**unverified** as to value, but the
   omission is structural).

**Action:** state the code-derived prior for `s` in §2.3, replace §2.4's two-point derivation with
it plus the Phase 3 back-solve, and drop the "more than one copy" reading of `a` until the fit
speaks.

### C2 — the price of Option B/C is 13–20%, not 5–6%, so §6's decision rests on the wrong number

**arithmetic-verified** (from published measurements; assumes the measured peak's non-`f0` part is
∝ B, which is true by construction of the fit).

"What σ×16 makes the promise hold at the CI fixture's 25 bp?", solved from each published
measurement (`allowed B_new/B_old = (930 − f0)/(peak − f0)`):

| platform / cell | peak | needs σ×16 | B shrink at that σ (c2 / c4) |
|---|---:|---:|---|
| Linux `25bp_1G_c2` mean 1.0301 | 958.0 | 25 | 6.1% / 4.9% |
| Linux `25bp_1G_c2` max 1.050 | 976.5 | 25 | 6.1% / 4.9% |
| macOS 1G/c4 | 1000.0 | 27 | 9.3% |
| macOS 1G/c2 min | 1048.1 | 29 | 16.2% / 13.3% |
| macOS 1G/c2 max | 1082.6 | **30** | **18.4% / 15.2%** |

The plan's 4.9–6.1% is the **Linux-only** answer. §4 says to take the worse of the two platforms, and
the handoff records the platform offset as ~17 points; the binding cell is macOS `1G/c2`, at
σ×16 ≈ 29–30. `--clump_only` (workers = 1) shrinks harder still — 18.2% at σ×16 = 29, 20.6% at 30 —
because its denominator is smaller (B = 26.2 MiB at 1G against 18.1 at c4).

And that is at `r_min` = the fixture's own 95.9 B. §7a's other rows (86 B for `--small_rna`, 61 B) land
at σ×16 ≈ 32 and ≈ 38 once `s` is right — 22% and 33% at c2.

**Action:** §6's "For/Against" needs re-costing before the option is called decided. A flat 13–20% bin
shrink for every user, at every read length, to protect the short tail, is a materially different
trade from 5%. Note this is not an argument that C is wrong — it is an argument that the plan chose
before it had the number it says the choice depends on, and the estimate it used is off in the
flattering direction.

### C3 — §7b's table is wrong at the row it calls decisive; the conclusion survives, the reasoning does not

**arithmetic-verified.**

The table scales the *whole* ratio by the bin shrink. Peak = static + pool, and the static term does
not shrink. It also leaves lever A at a flat +0.126 when lever A is itself B-proportional (one extra
queued batch per worker), so it should shrink too. Redoing both, with the committed baseline 0.7883,
trigger 0.8783, and lever A's weakest measured form (+0.126 at `1G_c4`):

| σ×16 | new `1G_c4` mean (static 173 / 205 / 224 MiB) | + lever A unscaled | + lever A scaled | verdict |
|---:|---|---:|---:|---|
| 25 | 0.759 / 0.761 / 0.762 | 0.885 / 0.887 / 0.888 | 0.879 / 0.880 / 0.881 | **FIRES** (by 0.000–0.009) |
| 26 | 0.745 / 0.748 / 0.749 | 0.871 / 0.874 / 0.875 | 0.862 / 0.865 / 0.866 | misses by 0.012–0.016 |
| 28 | 0.720 / 0.724 / 0.726 | 0.846 / 0.850 / 0.852 | 0.832 / 0.835 / 0.838 | misses by 0.041–0.047 |

The plan's row 1 says "0.877 — misses by 0.001". In all six treatments it is **≥ 0.878, i.e. it still
fires.** So the claim *"at every candidate the band would go green on the very regression it exists
for"* is not established for σ×16 = 25, which is the row §9 singles out as the one an estimate would
have got wrong. §9's self-assessment (*"Computing it turned that into a hard constraint"*) is
therefore resting on the arithmetic error it congratulates itself for avoiding.

Two refinements that go the plan's way and should replace the current argument:

- For σ×16 ≥ 26 the miss is 0.012–0.047 — robust to every static assumption. Since C2 puts the likely
  landing zone at 28–30, the plan's **conclusion** (remove the rows, re-band later from ≥10
  observations) is right; only its justification for the σ=25 row needs deleting.
- **The failure is not always silent.** The baseline file also fails a cell that drops more than 0.10
  below baseline. At σ×16 ≥ ~31 the `1G_c4` drop exceeds 0.10 and the refit PR goes **red on the
  lower rail** — a self-announcing failure, not a silent one. Worth saying, because it changes what a
  reader should expect to see in the refit PR's CI.

### C4 — §4 holds pair count constant across the read-length axis, which biases σ(short r) downward

**arithmetic-verified + handoff-sourced.**

`SESSION_HANDOFF.md` §5: *"Fixture size is the design. Ratios climb monotonically: 0.38 at 500 K
pairs … ~0.95 at 10 M."* The ramp is driven by flush cycles (≈ total bytes / B) and allocator
retention across them. At a constant 10 M pairs the total bytes vary with the fitted variable:

| read length | r (10 M fixture) | bytes/mate | flush cycles at B = 18.1 MiB |
|---:|---:|---:|---:|
| 25 bp | 95.9 | 0.96 GB | ~106 |
| 51 bp | 147.9 | 1.48 GB | ~163 |
| 150 bp | 345.9 | 3.46 GB | ~382 |

So the short-read cells — the ones the whole refit is about — get **3.6× fewer flush cycles** and sit
further from steady state, which depresses their measured peak and therefore biases the fitted `s`
**down**: the unsafe direction, on the axis being fitted. This may be part of why the plan's `s` is
low.

**Action:** equalise total bytes rather than pair count (`pairs ∝ 1/r`, e.g. ~31 M pairs at 25 bp for
parity with 8.7 M at 150 bp), or demonstrate a plateau with a pair-count ramp at the shortest length
before trusting any short-r cell. Either way, state the disk cost — the arc's own first trap is that
disk exhaustion inflates peak RSS.

### C5 — Step 1's harness repair is under-scoped, and its stated rationale is false

**code-verified.**

`calibrate.py` still encodes the **pre-Phase-3 model**:

| line | content | shipped model now |
|---|---|---|
| `:33` | `STATIC_MIB = 512` | 224 MiB (+24/FastQC thread) |
| `:53` | `denom_for = 5*n_bins + 7*cores` | `23*n_bins + 64*workers` |
| `:56` | `budget_mib_for = round(target_b*denom/4) + 512` | plus a 10/11 margin that does not exist here |
| `:61` | `exact_b_mib = 4*(budget − 512)/denom` | `16*(budget*10/11 − static)/denom` |

The fit regresses RSS on **`b_exact`** (this stale arithmetic), not on **`b_banner`** (what the binary
actually used), and nothing asserts the two agree — while `assert bins == n_bins_for(cores)` (`:107`)
*does* guard the bin count. Fix only the regex, as Step 1 says, and the harness runs, prints, and
fits **silently against the wrong x-axis**. That is the exact failure mode §4 warns about.

The highest-value single line in the repair is not the named-group regex — it is
`assert abs(b_banner - b_exact) < 1` per cell.

Also in this file, unmentioned by Step 1: `REPS = 2` (§4 wants 3); `--compression 6` hard-coded at
`:96` while the shipped default is level 1 and every CI cell measures at level 1 (writer speed changes
backpressure, hence queue occupancy, hence peak — so Phase 2's constants were fitted at a compression
level the CI cells they now gate do not use).

Finally, Step 1's justification — *"A silent no-match is what made this harness unusable without
anyone noticing"* — is **false**. `run_cell` (`:103-106`) raises `SystemExit` naming the log when
`m is None`. It failed loudly; §2.6's own text says so ("the harness exits on every cell"). Correct the
rationale so it points at the real silent path above.

---

## Important

**I1 — §3 contradicts the chosen option.** §3.1 requires `dyn_denominator` to "become a function of
`r`" and §3.2 requires "a rational `s/r` term" to keep the overflow and exact-inverse properties.
Option C does the opposite: `r_min` is folded into one compile-time integer and no `s/r` term ever
enters the code. Rewrite §3, or a reader implementing from §3 builds Option A. **code-verified**
against §6/§8.

**I2 — Validation 6's constant is wrong by construction.** The "~312 B" is exactly
`dyn_denominator/16` at cores = 32 (`23·128 + 64·32 = 4992`; `4992/16 = 312`). Raising σ raises the
bound: 328 B at σ×16 = 25, 352 B at 28. Assert `≤ dyn_denominator/16 + 1`, not a literal 312.
(Observed gaps in my sweep: 9–219 B, non-monotone because both sides truncate.) **arithmetic-verified.**

**I3 — the warning site misses `--clump_only`, and cannot see which input it is resolving.**
`clump_only.rs` calls `clump::resolve_layout` directly at `:286`, `:434`, `:779`, `:947` with its own
banners and never goes through `main.rs::resolve_clump_layout`. So a warning sited there covers 3 of
7 layout sites and skips the mode with the **largest B** (26.2 MiB at 1G against 18.1 at c4) — the
mode where a short-`r` breach is largest in absolute terms, and a *lossless* mode users may point at
arbitrary input. Separately, `resolve_clump_layout(cli)` takes only `cli`: it is called **per
input/pair** (`main.rs:1386`, `:1623`, `:1902`) but already consults `cli.input[0]` for the encoding
decision, so sampling `cli.input[0]` would warn about the wrong file for pairs 2..N. Passing the real
path is a 3-site signature change — which removes Step 7's stated reason for preferring this site
("avoids changing `sanity_check_with`'s signature at its five call sites"). The sampler must also go
through `format::open_sync_reader`: uBAM **input** with `--clumpify` is legal (only uBAM *output* is
rejected, `cli.rs:623-625`) and `FastqReader` would fail on it. **code-verified.**

**I4 — nothing validates the warning.** §5's nine rows have no cell for it. Needed: fires below
`r_min`; silent above; silent on the below-floor fallback path (`resolve_clump_layout` returns
`Ok(None)` before any banner); fires per pair in multi-pair runs. Note it *will* fire on the tiny
fixtures in `integration_preflight_tripwire.rs::se_trim_clumpify` and the `clump_only` integration
tests — harmless, because every banner assertion in the suite is `contains` (`tests/…:126`, `:169`,
`:592-605`), and no test asserts the absence of `WARNING`. Say so in the plan rather than discovering
it. **code-verified.**

**I5 — §4's fitted matrix drops the axis that identifies the slope.** At a single `--memory 1G` with
cores ∈ {2, 4}, `B` takes exactly two values (22.80 and 18.13 MiB) and each is tied to its cores
value, so per read length `f0` and the slope are determined by **two points** — the very weakness the
plan criticises in §2.4 ("two points cannot test linearity"). Phase 2 avoided it with five `B` values
× three cores, and its own notes warn about a "slope↔intercept artefact". Reinstate ≥3 budgets (or
≥3 cores) for at least the shortest and longest `r`. **code-verified** (`n_bins = max(16, 4·cores)`
gives 16 for both c2 and c4).

**I6 — the per-record cost lives in the worker term too, so the form is mis-specified.**
`channel_depths(true) = { work: 1, result_per_core: 2 }`, and `k = 4` charges *queued + in hand + 2
result slots*. The queued batch and the in-hand batch are `Vec<FastqRecord>` — spine-bearing; only
the result slots hold `Vec<u8>`. So roughly **two of k's four units carry the same 74 B/record**, and
the correct form is

```
dyn = B · [ (a + s/r)·n_bins + (a_k + s_k/r)·workers ]
```

Fitting the per-record term onto `n_bins` alone makes it absorb the worker part *at the fitted core
counts* and mis-extrapolate: the mis-attribution ratio is `(n_bins + 2w)/n_bins` = 1.25 at c2 but 1.5
at every c ≥ 4 (since `n_bins = 4c` there). Consistent symptom: the σ back-out from Phase 3's 25 bp
data is **core-dependent** (1.72 at c4, 1.81–1.91 at c2), which one σ(r) cannot fit. Sanity check the
other way: `a_k ≈ 2.6`, `s_k ≈ 148` reproduces the fitted `k` = 3.6–4.0 at 51 bp and predicts c4's
25 bp peak within 0.5%. Mechanism **code-verified**; the c2-vs-c4 gap is suggestive only
(**unverified** — macOS within-cell spread is 96 MiB, comparable to the gap).

**I7 — §7a's `--small_rna` rationale is a category error.** The bins hold **untrimmed input**
records: `read_pairs_clumpy` reads → keys → bins → flushes to a worker, and the worker trims
(`parallel.rs:836-870`). So `r` is the *input* mean record size; what a preset trims *to* never
reaches the pool. The reachability argument has to be about short **inputs** — 25–36 bp runs, 26 bp
10x/Multiome R1s, re-trimming already-trimmed data — and that changes the `r_min` number. (The
conclusion "short records are user-reachable" survives; the stated reason and the 86 B figure do
not.) **code-verified.**

**I8 — §2.3's `r` table was measured on 2000-pair fixtures, not the 10 M ones §4 will run.** This
generator puts the record index in the header (`@SIM:1:FCX:1:15:6329:{i} {mate}:N:0:ATCCGA`,
`examples/mk_clump_fixture.rs:57`), so mean header length is **37.44** chars at 2000 pairs and
**40.89** at 10 M. Every `r` is 3.44 B low; the fixtures that matter are 95.89 / 117.89 / 147.89 /
195.89 / 245.89 / 345.89 B. The plan's table reproduces the 2000-pair values to 2 dp, so this is what
happened. It matters directly: "92 B = 25 bp, the CI fixture" is really 95.9 B, and §2.3 itself sets
the rule that was broken ("any fixture used for fitting must have its `r` measured rather than
assumed"). **arithmetic-verified.**

**I9 — §2.4's two points are not two measurements of the same quantity.** (a) `s = 28.7` comes from
`r = 92/144`; the table's own 92.4/144.4 give 28.87. (b) σ = 1.4375 is the **shipped constant** — the
worse of two platforms, rounded *up* to the 1/16 grid — not a fitted value (fitted: 1.41 macOS, 1.29
Linux, CALIBRATION.md). (c) That constant was fitted on a fixture with **76-char headers**, i.e.
`r ≈ 183`, not 144.4. Pairing a rounded design constant taken at one `r` with a rough measurement at
another `r` gives a slope with no defined meaning; using 183 instead drops `s` to 21. Say which
quantity each point is, or drop the derivation and let the code-derived prior (C1) justify the form.
**arithmetic-verified.**

**I10 — the outline never flips the CI 25 bp cells to gating.** `ci.yml` carries
`measure 25bp_1G_c4 0` and `measure 25bp_1G_c2 0` with the comment *"25 bp is above the prediction
today: the per-bin coefficient is fitted per byte, not per record (#457). Recorded, not gated on the
promise, **until that refit lands**."* Step 4 removes baseline rows and adds a warning but never turns
the promise gate on for the cells #457 is about, and never deletes that comment. Without it the fix
ships unprotected — the arc's central lesson. **code-verified.**

**I11 — §7b reinvents `BOOTSTRAP`.** The baseline file and job already implement exactly the state
§7b describes: a `# BOOTSTRAP` line makes the job emit `::warning::baseline is in BOOTSTRAP mode -
drift band not enforced`, skip the armed-check loop, and skip the per-cell drift comparison, while
ratios still print and the promise assertion still gates. That is "deliberately off, loudly" with no
new bash in a job that has already been burned twice by silent-parse bugs — and it leaves the old
centres in the file for the re-band. Related mechanics worth stating correctly in §7b: `BANDED` is
**not** the enforcement mechanism; the `want=$(awk -v c="$label" …)` lookup inside `measure` gates on
the **rows**. Removing rows disarms the band; emptying `BANDED` only disables the arming check (and
leaving `BANDED` populated with the rows removed makes the job **error**, which is why both must move
together). **code-verified.**

**I12 — the falsification cells sit where the fixture cannot saturate.** §4 puts them at
`--memory 2G`, and the baseline file says of `2G_c4`: *"a 2 GiB pool needs ~2.3x this fixture to
saturate, so the cell cannot near its printed peak and gating it would prove little."* A
non-saturated cell will confirm the fit for the wrong reason — and §4 elsewhere acknowledges exactly
this ("a non-saturated cell understates the pool"). Either grow the fixture for those cells or
falsify at 1G on held-out (cores, length) combinations. **code-verified** (baseline file) +
**arithmetic-verified** (B at 2G/c8 = 21.0 MiB, pool ≈ 672 MiB of text).

**I13 — §4's run count is double.** 6 lengths × 2 cores × 1 budget × 3 reps = **36 per platform** (72
is both platforms), plus 12 falsification and 3 control. Worth fixing because the runtime/disk budget
should be explicit: 7 fixture pairs at 10 M pairs is roughly 6–9 GB once 150 bp is included, and
`fixture_stats` runs `sort -u` over 10 M sequences **per fixture**, which needs its own temp space.
**arithmetic-verified.**

---

## Optional

**O1 — the 1/16 quantisation is itself a design choice worth revisiting.** `σ×16 = ceil(16σ)` throws
away up to 6.25% of σ; §7a's own 86 B row is σ = 1.5718 → 25.15 → 26, and that rounding *doubles* its
cost (7.1% vs 4.9%). Scaling by 64 or 256 keeps integer arithmetic, cuts quantisation waste to 1.6%
or 0.4%, and makes `r_min` a continuous decision. The existing overflow test
(`floor_and_sizing_round_consistently_without_overflow`, budgets to `u64::MAX`, cores to 65536, run
under debug overflow checks) still covers it. **arithmetic-verified.**

**O2 — §6's floor claim is wrong: 281 → 283 MiB, not "~290".** *"Floors rise by the same ~5%"* is
false; the floor is `(224 MiB + denom/16 MiB) × 11/10`, so a 6.1% pool rise moves it **0.7%**. Values:
c2 281 → 283 (σ25) / 286 (σ28); c4 290 → 292 / 295. Note 290 *is* the current c4 floor, which may be
where the number came from. **arithmetic-verified.**

**O3 — Step 10's docs delta is the same error.** *"the worked values, the floor range and the
bin-size examples all move ~5–14%"*: the published floors move **+0.6% to +3.8%** at σ×16 = 25
(281→283, 590→608, 933→968, 1356→1391) and +1.5% to +9.4% at 28. Only the bin sizes move 5–14%.
Step 10's instruction (re-derive against the binary) is right either way. **arithmetic-verified.**

**O4 — §7a's c2 shrink at σ×16 = 28 is 13.9%, not 14.0%.** **arithmetic-verified.**

**O5 — no single-end and no `--clump_only` cell anywhere in §4/§5,** though both share
`dyn_denominator`, SE's per-record spine is 76 B rather than 74, and `--clump_only`'s B is 45% larger.
Phase 2's limitations flagged the same gap. One cell each at the shortest `r` would close it.

**O6 — Step 6 is feasible, and the plan should say so rather than leave it open.** Both the noise and
lever A scale with B, so the window `3× half-range ≤ width ≤ lever A − margin` survives the refit: at
σ×16 = 28, `1G_c4`'s range 0.055 → ~0.049 (3× half-range 0.074) against lever A 0.126 → 0.112. Tight
but open. **arithmetic-verified** under the assumption that dispersion is pool-proportional
(**unverified**).

---

## Alternatives the plan does not consider

### A1 — Option D: make the flush accounting memory-aware, and σ becomes a genuine constant

This is the fix I would put ahead of C.

The defect is not that σ is a constant. It is that **the bin flushes on `estimated_record_bytes`
(on-disk text) while the model charges memory.** Fix the accounting and the mismatch disappears:

```rust
/// Resident cost of holding one record in a bin: its text plus the Vec spine and key slot.
pub fn resident_record_bytes(rec: &FastqRecord) -> usize {
    estimated_record_bytes(rec) + PER_RECORD_RESIDENT   // 76 SE, 74 per record paired
}
```

used for `raw_bytes` and for `predicted` at the four push sites (`parallel.rs:903`, `:1287`,
`clump_only.rs:105`, `:148`). It must be a **new** helper, not a change to
`estimated_record_bytes`, which is also the on-disk byte counter (`clump_only.rs:235`) and the output
buffer sizer (`parallel.rs:353`, `:1080`). The `raw_bytes >= budget` comparisons need no edit.

What that buys:

- A bin's resident bytes are bounded by `B` **regardless of `r`**. σ loses its `r` term, so the
  promise holds at every read length: no `r_min`, no sampling, no warning, no bounded promise, no
  documented range, no §7a decision, no §2.3 header caveat (headers are text and already counted).
- It fixes **I6 for free**: the queued and in-hand batches are also bounded by `B` in memory terms,
  so `k`'s units become true memory bounds too. Under Option C, one σ cannot cover c2 and c4 at short
  `r` (C2's table needs 29–30 at c2 but 27 at c4) precisely because the worker term is mis-charged;
  Option D removes the cause instead of fitting around it.
- Determinism, reproducibility, tripwire and CI surfaces are untouched for the same reason Option C
  claims: `prediction = budget × 10/11` never involves σ. **code-verified** — the three pinned
  assertions (`tests/integration_clump_only.rs:126`, `…_ubam.rs:169`,
  `integration_preflight_tripwire.rs:602-605`) all pin the *prediction*, never `B`.

Honest costs: the σ refit is still needed (it is needed under every option); output bytes change
once, as §3.4 already accepts; the banner's `× N MiB` becomes a memory figure rather than a text
figure, which the docs must say; and the compression profile changes shape rather than vanishing.
With σ refitted to ~1.0 (σ×16 = 16, so `B` grows by 624/512 = 1.22× at c4), text per member relative
to today is ≈ 0.69 at 25 bp, 0.81 at 51 bp, 0.94 at 100 bp, 1.00 at 150 bp — i.e. **the cost is paid
by short-read data, in proportion to how much extra memory short reads actually consume**, and
long-read users pay nothing. Option C at the honest σ×16 = 29–30 charges *every* length a flat
13–20%. Against a safe Option C, Option D holds ~21% less text per member at 25 bp, ~6% less at 51 bp,
and **8–16% more** at 100–150 bp. Given that clumpify's documented gains are on ≥50 bp libraries,
that trade looks strictly better. (σ_D ≈ 1.0 is **unverified** — it is `φ·ρ`, which only the fit can
supply; the *shape* of the argument does not depend on its exact value.)

### A2 — pin `s` from the code instead of fitting it

A unit test in `parallel.rs::tests` (which already has `use super::*` access to `PairedBin` /
`SingleBin`) can push synthetic records with a known `r` into a bin and assert the reserved capacity
is `budget/(2r)` slots — turning the per-record constant into a **pinned property with zero
measurement noise**, and answering §5's Validation 2 (linearity in `1/r`) analytically. That leaves
§4 fitting only `(a, k, f0)`, which is what RSS measurement is actually good for. The motivating
arithmetic: at 1G/c4 one unit of σ is `16 × 18.13 MiB` = 290 MiB, so the arc's documented within-cell
spreads (96 MiB macOS, 37 MiB Linux) are **±0.33 and ±0.13 in σ** — the instrument is coarse relative
to the quantity, and 36 runs will not fix that.

### A3 — the header-length control cannot falsify as specified

Doubling a ~41-char header at 51 bp moves `r` from 147.9 to ~189, so with `s = 74` the effect is
Δσ ≈ 0.11 → `16 × 18.13 MiB × 0.11` = **32 MiB** of peak, against within-cell spreads of 96 MiB
(macOS) and 37 MiB (Linux). The cell §4 calls *"the cell that distinguishes the hypothesis from a
read-length correlation"* is 3× below the noise on the platform the plan measures on first. Make it
detectable: run it at the **shortest** length and the **largest B** (25 bp, `--cores 2`, where the
same doubling gives Δσ ≈ 0.23 → ~84 MiB), with a 4–8× header contrast, ≥5 interleaved reps, and
prefer Linux for its tighter spread — or replace it with A2, which is exact.
**arithmetic-verified.**

### A4 — cap the spine extrapolation (a defect independent of the refit)

`predicted = budget.div_ceil(bytes_this)` extrapolates a whole bin's reservation from the **first
record after each flush**. On variable-length input — already-trimmed data, which the plan mentions
only as an argument *against* Option A — a short first record inflates the reservation without bound:
a 20 bp first record in 150 bp data reserves ~5× the slots needed, worth roughly +2.3 in σ-equivalent,
i.e. a breach **no `r_min` can bound and no warning would predict**. The bin already knows its
previous fill count after the first flush; using that (or clamping `predicted`) is a small, local fix.
Worth its own issue regardless of which option Phase 4 takes. **code-verified** mechanism,
**unverified** magnitude.

### A5 — one thing about the site that is right and should be kept

`sanity_check_any` runs at `main.rs:270`, long before the three `resolve_clump_layout` calls
(`:1386`, `:1623`, `:1902`), and those calls are **per input/pair** — so a sampler there sees an
already-validated input and can warn per pair. If Option C survives, keep that ordering and say why;
it is the strongest part of Step 7. **code-verified.**

---

## What I did not check

- **No measurement of any kind.** Every peak-RSS number I use is quoted from `phase3/PLAN.md` §0.4,
  `phase3-close/PLAN.md`, `phase2/CALIBRATION.md` or `.github/clumpify-memory-baseline.tsv`.
- The σ back-solves in C1/C2 assume the measured peak's non-`f0` part is exactly ∝ `B` and take
  `f0` = 204.7 MiB (macOS) / 173.1 MiB (Linux) from CALIBRATION.md. Both are the fit's own
  definitions, but neither is re-verified here.
- `s ≈ 74` is the *reserved spine* only. I have **not** measured the allocator-rounding component
  (mimalloc size classes on `id`/`seq`/`qual`), which adds a further per-record term of order 10–25 B
  and is a **sawtooth in read length**, not smooth in `1/r` — two files with the same `r` can differ
  if their individual lengths fall either side of a size class. §4's chosen lengths straddle those
  boundaries unevenly (25→32, 100→112, 150→160), so residual "curvature" in Validation 2 may be
  quantisation rather than model error. **unverified**, but it is the most likely reason a clean
  `a + s/r` fit will show structure.
- Option D's compression figures depend on the unfitted σ_D; I have not built or run anything.
- Whether the 27-observation baselines are Linux-only (the job is `runs-on: ubuntu-latest`, and the
  file's `2G_c4` note says "First Linux reading", so I treated them as Linux).

---

## Action items, prioritised

| # | Action | Where | Priority |
|---|---|---|---|
| 1 | Replace §2.4's `s` with the code-derived 74 B/record prior + the Phase 3 back-solve; drop the "more than one copy" reading of `a`; add the fill fraction to §2.3 | §2.3, §2.4 | **Critical** |
| 2 | Re-cost Option B/C at σ×16 ≈ 29–30 (13–20%, 18–21% for `--clump_only`); reopen §6 as a recommendation pending Step 3 | §6, §7a | **Critical** |
| 3 | Recompute §7b holding static fixed and scaling lever A; delete the σ×16 = 25 "misses by 0.001" claim and §9's paragraph built on it; note the lower-rail red at large σ | §7b, §9 | **Critical** |
| 4 | Equalise total bytes (not pair count) across the read-length axis, or prove a plateau at the shortest length | §4 | **Critical** |
| 5 | Extend Step 1 to `STATIC_MIB`, `denom_for`, `budget_mib_for`, `exact_b_mib`, the 10/11 margin, `REPS`, `--compression`; add `assert b_banner ≈ b_exact`; fix the false "silent no-match" rationale | Step 1, §2.6 | **Critical** |
| 6 | Add Option D (A1) to §6 and decide against it explicitly if it is rejected | §6 | **Critical** |
| 7 | Rewrite §3.1/§3.2 to describe Option C | §3 | Important |
| 8 | Express Validation 6 as `≤ dyn_denominator/16 + 1` | §5 | Important |
| 9 | Cover `--clump_only`'s four layout sites; pass the real input path; use `format::open_sync_reader` | Step 7 | Important |
| 10 | Add validation rows for the warning (fires / silent / fallback / per pair) | §5 | Important |
| 11 | Reinstate a `--memory` or third-cores axis in the fitted cells | §4 | Important |
| 12 | Fit the per-record term on `n_bins + 2·workers`, or justify not doing so | §2.3, §4 | Important |
| 13 | Re-argue `r_min` from short **inputs**, not from what `--small_rna` trims to | §7a, Assumption 7 | Important |
| 14 | Re-measure `r` on the 10 M fixtures (95.89 … 345.89 B) | §2.3, §7a | Important |
| 15 | Label §2.4's anchor as the shipped constant, at `r ≈ 183`, not a 144.4 B fit | §2.4 | Important |
| 16 | Flip `25bp_1G_c4` (and `25bp_1G_c2` if it clears) to `gating=1` and delete the "until that refit lands" comment | Step 4/8 | Important |
| 17 | Use the existing `# BOOTSTRAP` marker instead of new empty-`BANDED` bash; state that the rows, not `BANDED`, arm the band | §7b, Step 8 | Important |
| 18 | Move falsification off the non-saturating 2G cells, or grow the fixture there | §4 | Important |
| 19 | Fix the run count (36/platform) and state the disk and `sort -u` budget | §4 | Important |
| 20 | Strengthen the header control (shortest length, largest B, 4–8× contrast, ≥5 reps) or replace with A2's unit test | §4, §5 | Important |
| 21 | Consider a 1/64 or 1/256 σ grid; fix the floor figures (283, not 290) and the docs delta (+0.6–3.8%); add SE and `--clump_only` cells; note Step 6 stays feasible; file A4 as its own issue | §6, §7a, Steps 5/6/10 | Optional |
