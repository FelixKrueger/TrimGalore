# Code review B — Phase 2, `--clumpify` memory accounting (#439)

**Reviewer:** B (independent; did not read `CODE_REVIEW_A_phase2.md`)
**Target:** uncommitted diff on `fix/439-clumpify-peak-memory`, over `70e75e3` (Phase 1)
**Touched:** `src/clump.rs`, `src/clump_only.rs`, `src/main.rs`, `src/parallel.rs`, `CHANGELOG.md`
**Recommend only — no edits applied** (shared working tree).

## Verification performed

| check | result |
|---|---|
| `cargo test --release` | **732 passed, 0 failed** (16 suites; matches the recorded state) |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean — first run returned in 0.26 s off cache, so I `touch`ed `src/clump.rs` and re-ran; it recompiled (`Checking trim-galore`) and stayed clean |
| model replica + brute force | Python replica of `resolve_layout` / `clumpify_min_memory_bytes` / `predicted_peak_bytes`, exact u64 semantics incl. `saturating_*`; swept `cores` 1–128 × `fastqc` × budgets from the per-cell floor to `u64::MAX` |
| mutation test | 17 perturbations of σ, k, `STATIC_TRIM_BYTES`, `STATIC_FASTQC_PER_CORE_BYTES`, the margin and `MIN_BIN_BYTES`, scored against all 8 sizing assertions |

## Verdict on the two headline correctness claims

**Probe 1 — `predicted_peak_bytes()` is genuinely the inverse, and the banner cannot exceed the budget. PASS.**

Proof: `B = ⌊usable·16/D⌋` ⟹ `B·D ≤ usable·16` ⟹ `⌊B·D/16⌋ ≤ usable`, so
`peak = static + ⌊B·D/16⌋ ≤ static + usable = budgeted = ⌊budget·10/11⌋ ≤ budget`. The two
roundings can only both round *down*, so they cannot conspire against the inequality.

Saturation is safe rather than accidentally safe. When `budget·10` saturates, `budgeted`
becomes `u64::MAX/11`; when `usable·16` saturates, `B = ⌊u64::MAX/D⌋`, and `B·D ≤ u64::MAX`
holds by construction, so the plain multiply in `predicted_peak_bytes` (`clump.rs:162`) cannot
overflow. Brute force over `{floor(c,f)} ∪ {256 MiB … 64 GiB} ∪ {2^40, 2^50, 2^60, 2^62, 2^63,
u64::MAX/10, u64::MAX/10+1, u64::MAX/2, u64::MAX}` × `cores ∈ {1,2,3,4,8,16,64,256,1024,2^20}` ×
`fastqc`: **0 violations, 0 overflows**. `peak/budget` is a flat 0.909 everywhere.

**Probe 2 — the advertised floor is exactly the smallest resolvable budget. PASS, and it is exact,
not merely conservative.**

`B ≥ M ⟺ usable·16 ≥ M·D ⟺ usable ≥ ⌈M·D/16⌉ = dyn_min`, and `⌊budget·10/11⌋ ≥ S ⟺ budget ≥
⌈S·11/10⌉`. `clumpify_min_memory_bytes` computes exactly that composition, so it is the least
resolvable budget, not an upper bound on it. Brute-forced `cores` 1–128 × `fastqc`: floor always
resolves, `floor − 1` never does (**0 violations**). The `× 11 / 10`-before-divide form flagged in
the iteration log is the fix that makes this exact — see Low 3 for the mutation evidence.

The `resolve_layout` bail message's suggested minimum (`div_ceil(1024*1024)`, `clump.rs:250`) is
achievable at every cell I checked. **The `main.rs` copy of that advice is not — see High 1.**

One consequence worth recording: because `main.rs:184` uses the exact floor as its fallback
predicate, `resolve_layout` at `main.rs:201` can no longer fail for the too-small reason at all,
and `cores == 0` is already rejected by `Cli::validate` (`cli.rs:727`). Its `?` is unreachable in
practice; `--clump_only` is the only caller that can see the bail.

---

## High

### H1. The below-floor warning advises a `--memory` that does not work

`src/main.rs:186-193` — `min_required / (1024 * 1024)` truncates, so the warning names a figure
below the true floor whenever the floor is not a whole number of MiB. The diff *fixed exactly this
class of bug* in the sibling message (`clump.rs:250` now uses `div_ceil`) and left the reachable
copy truncating.

Measured over `cores` 2–64 × `fastqc`: **124 of 126 cells give unusable advice.** Only
`--cores 16` without `--fastqc` lands on an exact MiB (418.00).

Concrete failure:

```
$ trim_galore --clumpify --cores 2 --memory 200M R1.fq.gz
WARNING: --memory 200M is too small for --clumpify at --cores 2 (need ≥ 280 MiB).
         Falling back to plain mode (no read reordering). …

$ trim_galore --clumpify --cores 2 --memory 280M R1.fq.gz
WARNING: --memory 280M is too small for --clumpify at --cores 2 (need ≥ 280 MiB).
         Falling back to plain mode (no read reordering). …
```

The true floor at `cores 2` is 294,125,568 B = **280.5 MiB**; `280M` = 293,601,280 B. The user
follows the instruction and gets the identical instruction back. Same at `cores 4` (289.30 → says
289), `cores 6` (310.75 → says 310), `cores 8` (332.20 → says 332), `cores 12`, `cores 32`, and
every `--fastqc` variant.

**Fix:** `min_required.div_ceil(1024 * 1024)` at `main.rs:192`.
**Test gap:** `clumpify_floor_matches_the_smallest_resolvable_budget` works in bytes, so it cannot
see this. Add, inside the same loop, an assertion that the *advised MiB figure* resolves:
`resolve_layout(floor.div_ceil(MIB) * MIB, cores, fastqc).is_ok()` — that is the number the user
actually types.

### H2. `--fastqc_args "-t N"` silently bypasses the per-core FastQC reservation

`src/clump.rs:174-180` derives the FastQC term from the `cores` it is handed, i.e. `cli.cores`.
But `src/fastqc.rs:85-90` lets `--fastqc_args` overwrite the thread count outright:

```rust
"-t" | "--threads" => {
    if let Some(v) = tokens.get(i + 1).and_then(|t| t.parse().ok()) {
        config.threads = v;   // overrides threads: cores.max(1) set at fastqc.rs:45
```

So the reservation and the phase it is reserving for can be sized from different numbers, in the
unsafe direction.

Concrete failure: `trim_galore --clumpify --cores 2 --memory 1G --fastqc_args "-t 16" R1.fq.gz`.
Sizing reserves `224 + 2×24 = 272 MiB`, resolves `B = 21.25 MiB`, and the banner promises a peak of
≈ 930 MiB. FastQC then runs with 16 threads, whose non-pool cost the Phase 2 notes measured at
**~283 MiB** — against 48 MiB reserved. The printed prediction stops being a bound, which is the
one property Phase 2 exists to establish.

Note `fastqc_requested()` (`cli.rs:592`) is `self.fastqc || self.fastqc_args.is_some()`, so
`--fastqc_args "-t 16"` alone reaches this path with no `--fastqc` present.

**Fix:** compute the effective thread count once and use it for both — either hoist the `-t` parse
so `resolve_clump_layout` can see it, or clamp `config.threads` to `cli.cores` in `fastqc::run` and
say so in `--fastqc_args`' help. The second is smaller and also makes `--cores` mean what the
threading docs say it means.

---

## Medium

### M1. The margin is proportional; the variation it is sized against is absolute — and the default budget falls on the wrong side

`src/clump.rs:201-204`. The doc comment justifies 10/11 with "Largest observed within-cell spread
is 96 MiB", and `PLAN.md:100` states the rule as "Margin must exceed the largest observed spread".
A 1/11 slice only exceeds a *fixed* 96.2 MiB once the budget passes `96.2 × 11 ≈ 1058 MiB`.

At the shipped default `--memory 1G` the margin is `1073741824 − 976128930 = 97,612,894 B =
93.09 MiB` — **3.1 MiB short of its own justification**, and it degrades from there: at the 281 MiB
floor the margin is 25.5 MiB against a 96 MiB spread.

No breach follows today (the default cell measured 0.862 of budget, ~141 MiB spare), because the
constants are conservative elsewhere. But the stated design rule is not satisfied over
[floor, ~1.03 GiB], which is the whole band a memory-constrained user operates in — and it is
precisely the band where a 96 MiB excursion is fatal.

**Recommend** one of: (a) `budgeted = budget − max(budget/11, 96 MiB)`, which makes the constant's
own sentence true at every budget; or (b) restate the doc comment as "9% of budget; exceeds the
96 MiB observed spread above ~1 GiB" so the limitation is on the record rather than implied.

### M2. `--clump_only` is charged a worker term it does not have, and the FastQC addend makes it bite

`--clump_only` is fully single-threaded — `clump_only.rs:33` says so, and every reader is sync
(`FastqReader::open` at 290/431, `open_sync_reader` at 767/969); `cores` reaches only
`fastqc::run`. Yet all four sites pass it to `resolve_layout` (`clump_only.rs:275, 413, 748, 905`),
which charges `64 × cores` for in-flight worker batches *and* sizes `n_bins = 4 × cores`.

Without `--fastqc` this costs about 4.4% of bin size — noise. With `--fastqc` the per-core
reservation compounds it:

| `--clump_only --memory 1G` | bin budget before | after | Δ |
|---|---|---|---|
| `--cores 2 --fastqc` | 21.79 MiB | 21.26 MiB | −2.4% |
| `--cores 4 --fastqc` | 18.96 MiB | 15.66 MiB | −17.4% |
| `--cores 8 --fastqc` | 9.48 MiB | 6.60 MiB | −30.4% |
| `--cores 16 --fastqc` | 4.74 MiB | **2.07 MiB** | −56.3% |

(without `--fastqc` the same cells move +4.7% / −4.4% / −4.4% / −4.4%.)

2.07 MiB is a hair above `MIN_BIN_BYTES`, whose own doc comment (`clump.rs:34-37`) says that below
1 MiB "gzip header overhead starts to erode the compression win". A `--clump_only --cores 16
--fastqc` user gets materially larger output than before, from a reservation for workers that do
not exist. The same invocation at `--memory 800M` now aborts outright, because these four sites
propagate the bail with a bare `?`.

This is the deferred Phase 1 item arriving. **The deferral is defensible** — the choice between
warn-and-fall-back and documented hard failure is Felix's, the plan says so plainly at
`PLAN.md:416-423`, and no test or documented invocation breaks today (the integration suite runs at
1 GiB / `--cores 2`, floor 280.5 MiB). But two things belong in the plan alongside it: that the
regression is a *compression* loss as well as an abort, and that `--clump_only` could sidestep both
by passing `1` for the worker term while keeping the user's `cores` for the FastQC term — the split
the code already has the parameters for.

### M3. FastQC as an addend rather than a `max` — the right call, undocumented as a deviation

`CALIBRATION.md:114-118` item 1 specifies FastQC "handled as a `max`, not an addend (FastQC runs
after trimming, so the peak is `max(trim, fastqc)`)". The implementation ships an addend
(`clump.rs:174-180`).

I think the implementation is **right and the calibration recommendation is wrong**: `fastqc::run`
is called after `run_single_end`/`run_paired_end` returns (`main.rs:1481-1488`), so the bin pool is
dropped — but under mimalloc the pages are retained rather than returned, so the FastQC working set
stacks on a still-resident trim peak. `max()` would under-reserve. The Phase 2 measurements
(191.9 MiB total at `--cores 2` with a minimal pool, i.e. *below* the 224 MiB reservation; then
76/139/283 MiB of non-pool cost at c4/c8/c16) are consistent with an addend and were the right
thing to measure.

The problem is bookkeeping: `PLAN.md:392-400` records the measurement and the 288 MiB→24 MiB/core
correction but never records that it overrode the calibration's stated form. Everything downstream
of that choice — the raised floor, M2's compression loss, the 0.682 ratio at `2 GiB/c16 --fastqc` —
traces to it, and a future reader comparing CALIBRATION item 1 to the code will read a silent
contradiction. One paragraph in the Phase 2 notes closes it.

Also worth flagging while the term is fresh: 24 MiB/core is a straight-line fit through
76/139/283, but the measured per-core rates are 19.0 / 17.4 / 17.7 MiB — i.e. flat, and the
constant is 35% above the worst of them. At `--cores 16` that reserves 384 MiB against 283 MiB
measured, which is where the `2 GiB/c16 --fastqc` cell's 32% unused budget comes from. The output
penalty was measured at c4 (+0.9%); at c16 the bins are 2 MiB and the penalty is unmeasured.

### M4. CHANGELOG — narrower than the entry it replaces, still ahead of the evidence

`CHANGELOG.md:8-11` now claims the budget is kept "at every core count and budget measured".
Against the measured cells that is **literally true**: `PLAN.md:378-388` lists six cells, all inside
budget (0.862 / 0.887 / 0.828 / 0.739 / 0.617 / 0.682), and no cell contradicts it. The trailing
qualifier is real hedging, which is a clear improvement on the High-rated over-claim it replaces.

What I would still change is the bolded clause, because the bold ends *after* the promise and
before the hedge, and the hedge is doing all the load-bearing work. The shipped model was verified
on one platform (macOS/arm64), one fixture (51 bp RRBS PE, ~52% distinct), paired-end only, at
`--cores` 4/8/16 and 1–2 GiB. Not verified with the new constants:

- **`--cores 2`** — the mode's own documented minimum (`cli.rs:732`), and the *only* cell where the
  new model makes bins **larger** (22.80 MiB against 21.79 MiB at 1 GiB, `CALIBRATION.md:95`). The
  one direction in which this change could reduce safety is the one core count not measured. It was
  a *fitting* cell, so the fit covers it (`114.6 + 30.454·B` ⟹ ~809 MiB at 1 GiB, ratio 0.79), but
  no shipped-binary run exists.
- **Linux with the new constants** — calibrated, not re-verified.
- Single-end, `--clump_only`, multi-pair (`PLAN.md` validation rows 8–10), low-duplication data.
- **Read length.** σ = 1 + spine/text. At 51 bp with 76-char headers that is ≈ 1.41, under the
  shipped 23/16 = 1.4375; at 150 bp it falls to ~1.21, so long reads are safe. Short reads are not:
  a 25 bp `--small_rna` record gives σ ≈ 1.55, **8% over** the pool coefficient. The pool is 59% of
  the denominator at `cores 4`, so ~5% of under-charge against a 9% margin that M1 already shows is
  thin at the default budget. `--small_rna` is a first-class mode and is absent from the matrix.

**Recommend** naming the matrix rather than the universe — e.g. "…inside the `--memory` budget it
prints across the measured matrix (macOS, `--cores` 4–16, 1–2 GiB)" — or dropping the bound framing
for a tracking claim ("the printed prediction is now within 5–10% of measured peak"). Length is
fine: at four lines this is the shortest entry in the file, and the file's local convention is far
more verbose than `CLAUDE.md`'s one-sentence rule (pre-existing, not this diff's problem).

---

## Low

### L1. `resolve_layout_predicted_peak_fits_budget` is not just weak — it is fully redundant

`clump.rs:534-559`. Mutation-tested against 17 perturbations; it passes **16 of 17**:

| mutation | T5 verdict |
|---|---|
| σ 23 → 20 (the old 5/4), 22, 24 | PASS |
| k 64 → 28 (the old 7/4), 63, 65 | PASS |
| `STATIC_TRIM_BYTES` 224 → 205, 225, 512 (the old value) | PASS |
| `STATIC_FASTQC_PER_CORE_BYTES` 24 → 0, 23, 25, 288 | PASS |
| margin removed (10/10), margin 9/10 | PASS |
| `MIN_BIN_BYTES` 1 MiB → 2 MiB | PASS |
| **margin inverted (11/10)** | **FAIL** |

So it catches exactly one thing — a margin sign error — and `resolve_layout_default` and
`resolve_layout_pins_the_default_cell` both catch that too. It is therefore redundant with tests
already in the file, not merely tautological. Keeping it until Phase 3's CI check exists is
defensible (`PLAN.md:332-336`), but it should not be counted toward coverage of the new constants.

### L2. Both bail tests are inert

`resolve_layout_bails_when_too_small` (`clump.rs:561-571`, 256 MiB / `cores 8`, floor 332.20 MiB —
76 MiB of slack) and `resolve_layout_bails_below_static_overhead` (`clump.rs:573-579`, 128 MiB /
`cores 4`, floor 289.30 MiB) **pass under all 17 mutations**, including reverting
`STATIC_TRIM_BYTES` to 512 MiB and doubling `MIN_BIN_BYTES`. They pin only "a very small budget
errors". Recommend `clumpify_min_memory_bytes(8, false) - 1` as the input, which makes them assert
the boundary they are named for.

### L3. The floor test uses `floor − 1 MiB` where `floor − 1` is available

`clump.rs:581-596`. The floor is provably exact (probe 2), so the weaker form passes — but it would
also pass if the advertised floor over-stated by up to 1 MiB − 1 B, which is the same order as the
error H1 is about.

Credit where due: I confirmed this test **does** bite on the two arithmetic-form bugs it was written
for. Replaying iteration-log #2 (`resolve_layout` doing `budget/11*10` while the floor inverts
`× 11/10`) fails it at `cores` 2/4/8/16; so does making the floor truncate instead of `div_ceil`.
It is the only test in the file that covers the floor↔sizing inverse relationship.

### L4. Stale comment naming a constant that no longer exists

`src/clump_only.rs:1188` — "resolve_layout requires STATIC_OVERHEAD (512 MiB)". The constant is
`STATIC_TRIM_BYTES` and it is 224 MiB. Test-only, but it is the "state the current fact" rule.

### L5. `ClumpLayout::new`'s fourth parameter feeds the safety claim and is unvalidated

`clump.rs:150`. The constructor is public, has no production caller (already logged in Phase 1's
"logged, not fixed"), validates neither invariant `resolve_layout` enforces, and now takes
`static_bytes` — the one input to `predicted_peak_bytes()` that is not derived. `parallel.rs:1615`
hard-codes `224 * 1024 * 1024`, a silent duplicate of `STATIC_TRIM_BYTES` with no compile-time link;
if the constant moves, the helper disagrees and nothing notices (the helper never calls
`predicted_peak_bytes`, so today it is harmless). Suggest `static_bytes_for(2, false)` at the call
site, or have `new` derive the reservation from `(cores, fastqc)` so it cannot be handed a number
the sizing never used.

### L6. Overflow hardening on an unbounded `--cores`

`clump.rs:176` — `STATIC_FASTQC_PER_CORE_BYTES * cores as u64`, a plain multiply. `--cores` is an
unbounded `usize` (`cli.rs:351-352`); `validate()` rejects only 0, and `< 2` for `--clumpify`.
`--clumpify --fastqc --cores 800000000000` overflows: panic in debug, silent wrap in release
(which would produce a *small* reservation and defeat the bound). `4 * cores` in `n_bins` has the
same shape and is pre-existing, but the new term overflows ~6 million times sooner. `saturating_mul`
costs nothing on this path.

### L7. Comment style vs `CLAUDE.md`

The three new constants (`clump.rs:191-204`) and the `resolve_layout` docstring
(`clump.rs:206-229`) carry measured numbers in source: "205 MiB on macOS/arm64 and 173 MiB on
Linux/x86_64", "76 MiB at 4 cores, 139 at 8 and 283 at 16", "σ = 1.29 (Linux) / 1.41 (macOS) and
k = 3.71 / 3.99", "largest within-cell spread was 96 MiB". `CLAUDE.md` asks for one line, two max,
stating the fact, with measurements in the commit message rather than the source.

**My judgement: the pointer is right, the numbers are a second copy.** For a calibrated model the
`See plans/…/phase2/CALIBRATION.md` line at `clump.rs:226` is exactly correct and should stay —
someone editing σ needs to know a document governs it. But with the pointer present, every measured
figure inline is a duplicate that can drift from the doc it cites, which is the failure mode the
house rule exists to prevent. Recommend one line per constant stating what it is, the pointer for
the derivation, and the measurements only in the commit message.

Two mitigating notes: the constant this replaced (`STATIC_OVERHEAD_BYTES`) had the same four-line
measured-value style, so this is local precedent rather than a new habit; and the σ/k derivation
block is genuinely load-bearing for anyone touching the denominator. If any of it stays, that block
is the part to keep.

### L8. The default cell moved 7× further than the calibration predicted

`CALIBRATION.md:96` and item 4 (`:123`) both say the default cell moves **−0.6%**, "deliberately".
The shipped constants move it **−4.4%** (19,884,107 → 19,006,356 B), because they round σ up from
the measured 1.41 to 23/16 = 1.4375 and k up from 3.99 to 4.0. Conservative and intended, but the
"the correction is free at the default budget" framing (`CALIBRATION.md:102`) is now off by an order
of magnitude, and the same 4.4% applies at `cores` 4/8/16 without `--fastqc`. One line in the
Phase 2 notes recording the actual figure keeps the two documents in step.

### L9. Banner cosmetics

`main.rs:203-208` truncates both figures, so the printed line cannot be reconstructed: at
`1G/c4` it reads "16 bins × 18 MiB; predicted peak ≈ 930 MiB", where `B` is really 18.13 MiB and
`224 + 18 × 39 = 926 ≠ 930`. Harmless, but the prediction now sits a flat 9.1% below *every*
budget, which will read as "why didn't it use the 1G I gave it?". One clause — "(9% held back for
run-to-run variation)" — turns a support question into a documented choice.

---

## Docs that this diff makes wrong (Phase 3 scope — stated plainly, as asked)

`docs/src/content/docs/performance/clumpy.md`:

| line | published | now |
|---|---|---|
| 93-99 | `usable = memory − 512 MiB` | `memory × 10/11 − (224 MiB, + 24 MiB × cores with --fastqc)` |
| 93-99 | `bin_byte_budget = 4 × usable / (5 n_bins + 7 cores)` | `16 × usable / (23 n_bins + 64 cores)` |
| 99 | "(FastQC + allocator + runtime overhead)" | FastQC is no longer in the base reservation |
| 102 | `--cores 8 --memory 1G` → "**32 bins × 12 MB**" | 32 bins × **9 MiB** (was already wrong: the old formula gave 4 MiB) |
| 102 | `--cores 8 --memory 4G` → "32 bins × 66 MB" | 32 bins × **44 MiB** |
| 108 | "minimum of ~535 MiB … stays in the 535–730 MiB range" | **276–590 MiB** without `--fastqc`, **302–1434 MiB** with it, over `--cores` 1–32 — and the floor now depends on `--fastqc`, which the sentence does not mention |
| 108 | "mostly from a fixed 512 MiB reservation for FastQC, allocator, and runtime overhead" | wrong on both the figure and the FastQC attribution |
| 112-116 | sample warning "at `--cores 6` (need ≥ 552 MiB)" | **310 MiB** as the code prints it today, **311 MiB** once H1 is fixed |
| 102, 108 | units written `MB` | banner now prints `MiB` |

`docs/src/content/docs/reference/changelog.md:199` also cites the 512 MiB reservation — that is a
historical release entry on the synced changelog page and should be left alone.

---

## Summary

The two properties Phase 2 set out to establish both hold, and hold exactly rather than by luck:
`predicted_peak_bytes()` is the algebraic inverse of the sizing with both roundings pointing the
safe way (verified to `u64::MAX`), and `clumpify_min_memory_bytes` is the least resolvable budget
rather than an upper bound on it (verified `cores` 1–128 × `fastqc`). The constants are the
worse-of-two-platforms choice the calibration recommended, the FastQC term is measured rather than
assumed and — despite contradicting `CALIBRATION.md` item 1 — is the *more* defensible form under
mimalloc retention. 732 tests, fmt and clippy all clean as recorded.

Two things I would not ship as-is: the below-floor warning tells the user a `--memory` value that
does not work (H1, one-line fix, and the diff already fixed the twin), and `--fastqc_args "-t N"`
lets FastQC's thread count diverge from the number the reservation was sized on, which breaks the
bound (H2). Below that, the margin's own justification is false at and below the default budget
(M1), and the FastQC addend charges a single-threaded `--clump_only` for workers it does not have,
halving its bins at high core counts (M2).

On test quality: `resolve_layout_pins_the_default_cell` and
`resolve_layout_fastqc_reserves_more_and_shrinks_bins` are the two that carry the new constants —
between them they fail under every perturbation I tried, including ±1 MiB on the FastQC term.
`clumpify_floor_matches_the_smallest_resolvable_budget` catches the exact arithmetic-form bug it was
written for. `resolve_layout_predicted_peak_fits_budget` survives 16 of 17 mutations and the one it
catches is caught twice elsewhere — it is redundant, not merely vacuous. The two bail tests survive
all 17 and sit 76 MiB from the boundary they name.
