# Plan — Phase 3: defend and correct the `--clumpify` memory promise (#439), v2

**Parent plan:** `../PLAN.md` v3 (Phase 3 = its "defend and correct" stage)
**Issue:** [#439](https://github.com/FelixKrueger/TrimGalore/issues/439) — deliberately still open
**Reviews:** v1 → `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` (both: "implement with the listed changes")
**Baseline:** `dev` at `fdade7e`, confirmed equal to `origin/dev`
**Prerequisites:** none outstanding. Nothing in flight.

## Revision history

- **v1 — reviewed, three Criticals, two of them agreed by both reviewers.** Its docs
  arithmetic, its `file:line` references (all twelve verified by both reviewers) and its
  three corrections to `SESSION_HANDOFF.md` were right. Both of its *guards* were not:
  the swept property test was blind to every coefficient, and the CI fixture was undersized
  by 13–17 points of ratio. v1 also drew a structural conclusion ("the CI job cannot be a
  drift tripwire") from the undersized measurement.
- **v2 — this document.** Fixes both guards, adds the two decisions taken after review
  (D3, D4), and records a defect found while resolving the one place the reviewers
  contradicted each other: **every measurement in this arc — Phase 2's calibration, both
  reviewers', and v1's — was taken on the plain-text output branch or compared against
  numbers that were.** See §0.

## 0. What resolving the reviewers' contradiction turned up

Reviewer A reported that `peak ≤ prediction` **does not hold** today: one rep in six at
10 M pairs / `1G/c4` measured 990.9 MiB against a printed 930 MiB. Reviewer B never ran
10 M and assumed the ratio was ~0.95, i.e. that it holds. Since the choice of CI assertion
turns on this, it was measured again with both agents finished and the machine idle.

**A's breach did not reproduce at c4** — 10 reps, max 908.3 MiB (0.977 of prediction).
**But `1G/c8` breached reproducibly**: 5 reps at 971.3 / 1042.8 / 1064.3 / 1073.0 /
1081.8 MiB — 4 of 5 over the 930 MiB prediction and **3 of 5 over the 1024 MiB budget**.

Chasing that found the cause, and it is not the memory model:

| cell (10 M pairs, 51 bp, 100% distinct, macOS/arm64) | plain output | **gz output** | of prediction | Phase 2 calibration |
|---|---|---|---|---|
| `1G/c4` | max 908.3 MiB (10 reps) | max **873.2** MiB (3 reps) | 0.939 | 0.949 ✓ |
| `1G/c8` | max **1081.8** MiB (5 reps) | max **814.4** MiB (3 reps) | 0.876 | 0.912 ✓ |

On the gz path the promise holds with 6–12% headroom and reproduces Phase 2's calibration
to within ~1–4 points. Every breach lives on the plain-text output branch, which is the
branch the parent plan already documented as pathological ("the plain-text branch peaked
802 MiB and the gzip branch 508 MiB from an identical 82 MB reserve") and which nobody
noticed they were measuring.

**Why it went unnoticed for three rounds.** Every harness in this arc asserts that the
`clumpify:` banner appeared, because a silent fallback to plain *mode* reads as a ratio of
~0.1. Nobody guarded the adjacent axis — clumpify engaged, on the plain-*output* branch —
and the banner prints `(gzip level 1)` while writing uncompressed bytes, so the one artifact
everyone checked actively asserted the thing that was not true.

### The defect this exposes (D3: its own issue)

```
$ trim_galore --clumpify --cores 2 in.fastq          → in_trimmed.fq   (plain; clumpify runs)
$ trim_galore --clumpify --gzip --cores 2 in.fastq   → in_trimmed.fq   (--gzip ignored)
$ trim_galore --clumpify --dont_gzip --cores 2 in.fastq
  Error: --clumpify and --dont_gzip are mutually exclusive (clumping plain text is pointless)
```

`--gzip` warns *"deprecated in Trim Galore v2.0. Output is gzipped by default. Use
--dont_gzip to disable. Ignoring."* — and that sentence is false for plain input: output is
plain, and there is no other way to request gzip. So `--clumpify` refuses the explicit
request to clump plain text as pointless, then does exactly that whenever the input happens
to be plain, on the one branch where the `--memory` promise measurably breaks. Note also
that `modes/clump-only.md:78` documents `--dont_gzip` as **supported** under `--clump_only`
(`*_clumped.fq`), so the two clumping modes take opposite positions on plain output.

Three coupled parts, all verified on `fdade7e`:

1. `--gzip`'s deprecation message is wrong whenever input is plain.
2. There is no way to obtain gzipped output from plain input.
3. `--clumpify` runs the plain-output branch it declares pointless, and there peak RSS
   reaches **1.057 × `--memory`** at `--cores 8` (10 M pairs, measured, 3 of 5 reps).

**Routing (D3):** its own issue. Phase 3 stays a defend-and-correct phase and does not grow
a behaviour change. Phase 3 gains exactly one requirement from it: **the CI fixture must be
gzipped**, which both reviewers wanted for disk and sensitivity and which is now known to be
required for the measurement to be valid at all.

## Decisions

| # | Decision | Chosen | Rationale |
|---|---|---|---|
| D1 | `--clump_only` below the layout floor | **Documented hard failure** | Reordering *is* the whole job; a fallback would write a byte-copy of the input at full I/O cost and exit 0. The message is already actionable. Behaviour is unchanged, so the change carries zero regression risk. |
| D2 | Banner margin visibility | **Name the budget**: `predicted peak ≈ 930 of 1024 MiB budget` | Makes the margin visible arithmetic without the banner asserting anything about it, and hands a cluster user the `--mem` number directly. |
| D3 | The plain-output defect (§0) | **Its own issue**; Phase 3 unchanged apart from a gzipped CI fixture | Keeps #439 closable on merged work and keeps Phase 3 free of a behaviour change with its own calibration question. |
| D4 | The docs corrections | **Ship first, as a docs-only PR** | Pure arithmetic, verified against the binary, blocked by nothing, and the repo rule requires docs PRs to contain only docs changes. Stops eleven wrong figures being published while the guards are reworked. |

## Context

### Corrections to the handoff (both reviewers verified all three)

| Handoff | Actual | Why it differs |
|---|---|---|
| banner prints `≈ 931 MiB` | **930 MiB** | integer truncation in `predicted_peak_bytes() / (1024*1024)` |
| floor range `276–590` / `302–1012` MiB | **281–590** / **334–1012** MiB, for `--cores` 2–32 | 276/302 are the `--cores 1` values, reachable only through `--clump_only` |
| `--clump_only --fastqc --cores 16 --memory 800M` now aborts | **it resolves** — that cell's floor is **775 MiB** | Phase 2's post-review `workers: 1` split lowered it after the note was written |

Reviewer B adds a consequence worth carrying: the third row means the *only* concrete case
Phase 2's deviation note offered for D1's urgency does not exist. D1 is still right; its
urgency was overstated upstream.

### The `--clump_only` regression surface (two cells corrected after review)

Both reviewers independently recomputed this table and found the same two errors. All
values `workers: 1`, `fastqc_threads = cores`, `div_ceil` to MiB as the binary advises:

| `--cores` | old floor (v2.3.0) | new floor | previously-working budgets that now abort |
|---:|---:|---:|---|
| 2 | 535 MiB | **329 MiB** (v1 said 277 — the *plain* floor, wrong series) | none — floor dropped |
| 8 | 566 MiB | 513 MiB | none |
| 10 | 580 MiB | **579 MiB** (v1 said 578 — truncated where the binary uses `div_ceil`) | none |
| 11 | 586 MiB | 611 MiB | [587, 611) MiB |
| 16 | 620 MiB | 775 MiB | [620, 775) MiB |
| 32 | 728 MiB | 876 MiB | [728, 876) MiB |

Widest window is **154.4 MiB at `--cores 16`** — v1's "≤ 150 MiB" was wrong. Without
`--fastqc` the floor dropped at every core count (620 → 352 MiB at `--cores 16`), so nothing
regresses there. Both reviewers confirmed the crossover is exactly `--cores ≥ 11` (swept
1–256 and 1–64 respectively). The narrowness of this surface is what makes D1 cheap.

### Files this phase touches

| File | Role | PR |
|---|---|---|
| `Docs/…/performance/clumpy.md` | the arithmetic | docs (D4) |
| `Docs/…/modes/clump-only.md` | D1's documented behaviour; the `--cores` parity contradiction | docs (D4) |
| `src/main.rs:207-213` | the `clumpify:` banner (D2) | code |
| `src/clump.rs` | `channel_depths()`; pinned-cell table; renamed sweep | code |
| `src/parallel.rs:117`, `:943`, `:130`, `:959` | four depth literals → the new helper | code |
| `src/clump_only.rs:275/421/764/929` | D1 bail context; banner-wording tidy (Q4) | code |
| `tests/integration_preflight_tripwire.rs:599` | arithmetic assertion, not a label assertion | code |
| `.github/workflows/ci.yml` | new `clumpify-memory` job | CI |
| `examples/mk_clump_fixture.rs` | new, deterministic **gzipped** fixture generator | CI |
| `phase2/CALIBRATION.md` | amend item 1 | either |
| `CHANGELOG.md` | extend #439's existing bullet | code |

`Docs/…/reference/changelog.md:198-199` cites the 512 MiB reservation and "within 1%" and
is **deliberately left alone** — it is the synced changelog page. It and `clumpy.md:96`/`:108`
are the only 512 references in `Docs/`; a grep sweep confirms it.

## Behavior

### 1. Banner (D2)

```
clumpify: 32 bins × 9 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
```

Second figure is `memory_bytes / (1024*1024)` — the parsed `--memory`, so `1500M` prints
1500 (verified). Sub-MiB budgets truncate to 0 but are unreachable: the smallest resolvable
budget is 281 MiB and everything below falls back to plain mode before the banner. Reviewer
A notes the echo is a truncation of the parsed value rather than of the user's string
(`1073741825` prints 1024); harmless for the same reason.

**Guard consequence found by reviewer B.** After D2 the predicted peak has **no unit token of
its own** — the `MiB` belongs arithmetically to the budget. v1's replacement clause
`"MiB budget (gzip level"` pins the *budget's* divisor, not the peak's, so changing
`predicted_peak_bytes() / (1024*1024)` to `/ (1000*1000)` yields `≈ 975 of 1024 MiB budget`
and **both clauses still pass**. That re-opens the divisor bug Phase 2's post-review fix #4
closed. Fix in Step 1.

### 2. `--clump_only` below the floor (D1)

Behaviour unchanged — non-zero exit, no output file. What changes: the bail names the mode
(today it says only "--memory budget too small for --cores 16"), and the behaviour is
documented. Reviewer B verified the context wrap cannot break a test: nothing asserts on the
bail text and `small_memory_budget()` is 1 GiB, above every `workers: 1` floor.

Edge case: at `--cores 1` the floor is 277 MiB and lowering `--cores` is not available, yet
the message still ends "or decrease --cores". The docs paragraph must say raising `--memory`
is the only remedy there.

### 3. Docs arithmetic (D4)

Eleven wrong statements in `clumpy.md`, in file order. Both reviewers verified every line
number and every replacement value.

| Line | Wrong | Correct |
|---|---|---|
| 91 | "the formula picked in an attempt to get the predicted peak RSS ≤ `--memory`" | drop the hedge: the coefficients are fitted to measured peak RSS on both platforms |
| 96 | `usable = memory − 512 MiB` | `usable = memory × 10/11 − static` |
| 96–97 | "(FastQC + allocator + runtime overhead)" as one flat term | static = 224 MiB, **plus** 24 MiB per FastQC thread only when a report is requested, capped at 16 charged threads |
| 98 | `4 × usable / (5·n_bins + 7·cores)` | `16 × usable / (23·n_bins + 64·workers)` |
| 98 | `cores` in the denominator | `workers` — `--cores` on the trim path, 1 under `--clump_only` |
| 102 | "32 bins × 12 MB" at `1G/c8` | **32 bins × 9 MiB** |
| 102 | "32 bins × 66 MB" at `4G/c8` | **32 bins × 44 MiB** |
| 108 | "needs a minimum of ~535 MiB" | **281 MiB** at `--cores 2` |
| 108 | "stays in the 535–730 MiB range for any sensible core count" | **281–590 MiB** without `--fastqc` and **334–1012 MiB** with, **for `--cores` 2–32** |
| 108 | "mostly from a fixed 512 MiB reservation for FastQC, …" | 224 MiB trim reservation; FastQC conditional and capped |
| 113 | "need ≥ 552 MiB" at `--cores 6` | **need ≥ 311 MiB** |

Six additions the reviewers found:

- **Name the core span.** At `--cores 64` the floor is **933 MiB** plain / **1356 MiB** with
  `--fastqc` — above the 1 GiB default. "Any sensible core count" is carrying real weight,
  so the corrected sentence must say 2–32.
- **Define the two new terms.** Step 3 introduces `× 10/11` and `workers` into a page that
  defines neither. One clause each.
- **Give the 9% margin a documented home.** Phase 2's reviewer B predicted D2's support
  question verbatim ("why didn't it use the 1G I gave it?"). One sentence naming the ~9% as
  deliberate run-to-run headroom discharges it without touching D2.
- **Date the stale savings figures** (open question 1, both reviewers). The RRBS row
  advertises "~24% at default `1G`; up to ~31% at `4G+`" while the `4G/c8` bin moved
  66 → 44 MiB, a 33% cut in exactly the input to that claim — and line 104 ("increasing
  memory doesn't seem to make very much difference") already contradicts the RRBS row, and
  did so before Phase 2. Do not fabricate replacements: add one clause dating them to
  pre-refit sizing, resolve the internal contradiction, and re-measure the RRBS row before
  the next release. Dating a claim is the difference between "known stale" and "published
  wrong".
- **`clumpy.md:110`** ("prints a warning and falls back to plain mode") needs a one-clause
  `--clump_only` exception, since D1 makes that path a hard failure and the pages cross-link
  at line 11.
- **`clumpy.md:119`** ("Memory usage without `--clumpify` … around the 100MB mark") loses
  its figure and points at the threading page. The number is known wrong and the correct one
  needs Step 7's measurement.

And in `clump-only.md`: **line 76 claims `--cores` is "accepted for interface parity"**,
which Step 2's new paragraph directly contradicts — `--cores` sets `n_bins` and therefore
moves the floor (277 MiB at c2 → 454 MiB at c32). Amend it.

### 4. Deterministic drift guards — rebuilt

v1's two guards are replaced. Both reviewers proved v1's sweep blind to **every** constant;
B ran the matrix:

```
sigma 23→20  PASS   k 64→32  PASS   static 224→512 MiB  PASS
fastqc 24→0  PASS   cap 16→1024  PASS   margin REMOVED  PASS
margin inverted  FAIL (650 violations)   floor/resolve rounding disagree  FAIL (28)
```

The reason is algebraic, so no extra cells fix it:
`b = floor(usable×16/d)` ⇒ `b·d ≤ usable×16` ⇒ `dyn ≤ usable` ⇒ `predicted ≤ budgeted ≤ budget`
for **any** σ, k, static. v1's headline assertion is a restatement of `MARGIN_NUM ≤ MARGIN_DEN`.
`Err`-below-floor/`Ok`-above is likewise a restatement of `clumpify_min_memory_bytes` being
the algebraic inverse of `resolve_layout`; both read the same `dyn_denominator` and
`static_bytes_for`, so they move together.

So Phase 3 ships **three** guards, each named for what it actually proves:

- **4a — pinned-cell table.** `resolve_layout_pins_the_default_cell` (`clump.rs:546`) is
  today the only test in the repo that fails under a coefficient perturbation (both
  reviewers confirmed; Phase 2's reviewer B recorded it failing under ±1 MiB on the FastQC
  term). Extend that *style* to a small table of hand-computed
  `(budget, LayoutInputs) → exact bin_byte_budget`, covering: the default cell,
  `workers ≠ bin_cores` (the `--clump_only` shape), the FastQC cap boundary at 16 **and** 17
  threads, and `--cores 2` and `--cores 32`. Six cells. Values are one Python line each.
- **4b — the sweep, renamed.** Keep it for the two properties it does prove: **rounding
  direction** across the input space (it catches the real Phase 2 iteration-#2 bug —
  `budget / 11 * 10` vs `budget * 10 / 11` — at 28 violations) and **saturation** at
  `u64::MAX`. The latter is not idle: `predicted_peak_bytes()` (`clump.rs:162`) and
  `clumpify_min_memory_bytes` (`:216-217`) use unchecked `*`, safe today only because
  `resolve_layout`'s `saturating_mul(16)` bounds `B`. In a debug `cargo test` an overflow
  there panics, so the sweep pins three unchecked multiplies. Name it for that.
- **4c — channel-depth pin, both channels.** v1 pinned only the work channel. `clump.rs:247-248`
  states what `k = 4` pays for: *"the queued batch, the batch in hand, and the compressed
  output waiting on the result channel"* — and the third is `sync_channel(cores * 2)` at
  `parallel.rs:130` and `:959`, which v1 left unguarded. The parent plan names it as the
  next lever someone will pull for the `2G/c4` residual. One helper returning both depths,
  both arms pinned.

**Stated honestly, per reviewer A:** 4c catches an edit to the helper, not a re-inlined
literal at the call sites — and the literal is what stands there today, so re-inlining is
the path of least resistance. v1 called this "the guard for the exact regression the handoff
fears"; that claim is dropped. The extraction is still worth doing: it removes the
duplication and puts the `k`-coupling note beside `dyn_denominator`, which is where the repo's
comment rule wants an invariant.

### 5. CI guard — resized, gzipped, and hardened

v1's §5 concluded the job "structurally cannot be a drift tripwire". Both reviewers refuted
that as a fixture artefact. B measured the ramp at `1G/c4`, monotone throughout:

```
500 K → 0.378    1.5 M → 0.754    3 M → 0.804    6 M → 0.875    10 M → 0.949 (Phase 2)
```

No cell is at steady state at 3 M, and the gap to the 10 M calibration is 13–17 points, not
the "~3 points" v1 claimed. The consequence inverts v1's design conclusion: at 3 M a
lever-A revert takes `1G/c4` from ~776 to ~880 MiB, under the 930 MiB prediction and
invisible; at 10 M the base is ~873–883 MiB and the same +104 MiB **breaches the prediction**,
so the single assertion v1 already specified *is* the drift tripwire — with no
runner-specific baseline and no threshold to tune. **Open question 2 dissolves.**

v1's supporting argument ("buying sensitivity means buying fixture size, without bound") was
about scaling to larger budgets. Sensitivity is needed at *one* cell, not every budget, and
the smallest cell is the cheapest to saturate — the opposite of v1's conclusion.

**Shape of the job:**

- **Fixture: 10 M pairs, 51 bp, gzipped.** Gzip is not an optimisation here — §0 shows the
  plain-output branch is a different memory regime, and CI must measure the branch users
  run. It also takes the fixture from 2.96 GB to ~490 MB.
- **Cells: `1G/c4` (gating) and `1G/c4 --fastqc`.** `2G/c4` is **dropped** on B's argument:
  its pool is 672 MiB, needing ~2.3× the data to saturate, so at any affordable size it is
  the least informative of the three, costs the most RSS, and is the sole reason
  Assumption 4's "≥ 4 GiB RAM" is load-bearing.
- **Assertions:** `peak ≤ prediction` (gating); `prediction == budget × 10/11 ± 1 MiB`;
  `peak ≥ 0.5 × prediction` (lower bound); the `clumpify:` banner is present; **the output
  is `.fq.gz`**. Log `peak/prediction` and `peak/budget` every run.
- **Expected ratios:** ~0.94 at `1G/c4` and lower under `--fastqc`, on macOS. Linux measured
  ~9 points below macOS in Phase 2, so expect ~0.85 / ~0.72 on the runner.

**Two self-referential holes reviewer B found, both closed by the assertions above.** First,
Step 6.4 sells parsing the budget off the banner as "the job never restates a number the
binary owns" — but then the binary owns *both* sides of the inequality, and since
prediction ≡ `budget × 10/11` identically, deleting the margin inflates the prediction to the
budget, real RSS grows ~10%, and `peak ≤ prediction` still passes. Hence the margin
assertion. Second, nothing guarded the *parse*: an empty capture without `set -euo pipefail`
compares as 0 and **passes**. Hence non-empty integer checks — the same discipline v1's own
Validation-10 preamble demands of everything else.

**The lower bound replaces v1's Validation 11 eyeball.** It catches a silent fallback *and*
a fixture that stopped filling the pool, and unlike a ratio ceiling it is runner-independent,
because it detects a 3× miss rather than a 5% drift.

### 6. The `threading.md` issue — re-measure before filing

v1's headline figure is input-limited, which v1 flagged for its 50 kb cell and then failed to
apply to the 10 kb row it was quoting. A 400 MB fixture of 10 kb reads is ~4.7 batches of
4096 records against ~20 slots of resident capacity. Both reviewers measured past it:

| records / input | reviewer A | reviewer B |
|---|---|---|
| 10 kb, 400 MB | 768.5 MiB | 813.2 MiB (c4), 817.4 (c8) |
| 10 kb, 764 MB–1.2 GB | 1366.0 MiB | 1808.8 (c4), 2052.3 (c8) |
| 10 kb, 1.5 GB | 2121.2 MiB | — |
| 10 kb, 3.0 GB | **2282.6 MiB** ← plateau | — |

So the honest headline is **~20–23× the documented ~100 MB**, not 8×, and v1's "flat in cores
past 4" is a fixture artefact — at 1.2 GB, c8 exceeds c4 by 244 MiB. The `--cores 1 →
~10 MiB` figure reproduces (10.1 / 10.2 / 10.6 across three harnesses) and is the sharpest
fact in the issue: the pool is bypassed entirely.

Also for the issue:

- **`threading.md:63` carries a second wrong claim in the same sentence** that v1 missed:
  *"Peak RSS scales with read length, **not input file size**"*. The tables above refute it
  at fixed read length. Correcting the ceiling and leaving this standing is a half-fix of one
  sentence, and it is the sentence that tells a scheduler user they can ignore input size.
- **Three batch constants, not two.** v1 named `fastq.rs:220` and `bam.rs:38`;
  `parallel.rs:33 BATCH_SIZE = 4096` is the one sizing the round-robin batches in the path
  being measured, and its own docstring says "4096 records × ~300 bytes ≈ 1.2 MB per batch"
  — the same read-length-blind assumption.
- **Lead with the gzip-output figure**, quoting plain text as the worst case. §0 is the
  reason: `threading.md`'s claim is about the normal path, and all three harnesses measured
  plain output, so a first responder gzipping their input will fail to reproduce.

Scope: **file only.** The fix (a byte-budgeted batch, as clumpy bins already use) touches the
path every mode shares and needs its own plan.

### 7. Issues to file

| Issue | Content | Bar |
|---|---|---|
| **plain-output `--clumpify`** (D3, §0) | the three coupled parts, the measured breach table, and the `--clump_only` asymmetry at `clump-only.md:78` | data-integrity: breaks the `--memory` promise at 1.057× budget |
| **`threading.md` memory claims** (§6) | the re-measured table, the plateau, both wrong sentences, three batch constants | user-visible wrong documentation |
| **EOF bin flush swallows a failed send** | `parallel.rs:714`, `:754`, `:1202` use `let _ = …send(…)` where `:794` and `:825` check `.is_err()`; a dead worker silently drops the final partial batch and the process can exit 0 | **both reviewers**: clears the repo's issue-filing bar more easily than the docs sentence does |
| **`--clump_only` predicted peak** (Q3) | it prints bins × MiB but no prediction, so #439's visibility never reaches the mode whose below-floor path D1 just made a hard failure | reviewer B: a recommendation living only in a plan document is the shape of thing that rots |

Carried from Phase 1's "logged, not fixed" and now routed: `ClumpLayout::new` (`clump.rs:150`)
is public, has no production caller, and validates neither invariant `resolve_layout`
enforces — so any future caller can build a layout whose `predicted_peak_bytes()` is
arbitrary, and the CI guard would faithfully report the lie and pass. Steps 4a–4c already
open the file; make it `pub(crate)` or have it assert `bin_byte_budget >= MIN_BIN_BYTES`.

### 8. Unverified axes — disposition (both reviewers flagged the silence)

`SESSION_HANDOFF.md` §2's last row was dropped from v1 entirely. Two of the four are
load-bearing for a phase about defending the promise:

- **Short reads.** σ ≈ 1.55 measured at 25 bp against the shipped 1.4375 means there is a
  *known* read length at which the model under-reserves — the promise-breaking direction.
  The CI fixture is 51 bp, the one length σ was fitted on, so the guard is structurally
  unable to see it. **Add a second read length to the generator** (25 bp) as a logged,
  non-gating cell; one extra generated fixture in the same job.
- **Multi-pair.** This follows from Phase 2's own central discovery: FastQC is an addend
  rather than a `max` *because mimalloc keeps the pool's pages resident across a phase
  boundary*. Bins are re-allocated per pair, so by the same mechanism pair 2's pool may sit
  on pair 1's retained pages — which would breach a per-process promise, and nothing tests
  it. **Add one multi-pair cell** (the same fixture passed twice).
- **Single-end and low-duplication/WGBS:** explicitly deferred, named here so they stay in
  the paper trail rather than falling out of it.

## Signature

```rust
// src/clump.rs, beside dyn_denominator so the coupling is visible where k lives.

/// Work-channel and result-channel depths. A clumpy batch is a whole bin, so shallower
/// queues cap resident bins; `k` in `dyn_denominator` was fitted at these depths.
pub fn channel_depths(clumpy: bool) -> ChannelDepths {
    if clumpy {
        ChannelDepths { work: 1, result_per_core: 2 }
    } else {
        ChannelDepths { work: 2, result_per_core: 2 }
    }
}
```

```rust
// src/main.rs — D2
eprintln!(
    "clumpify: {} bins × {} MiB; predicted peak ≈ {} of {} MiB budget (gzip level {})",
    layout.n_bins,
    layout.bin_byte_budget / (1024 * 1024),
    layout.predicted_peak_bytes() / (1024 * 1024),
    memory_bytes / (1024 * 1024),
    cli.compression,
);
```

```
examples/mk_clump_fixture.rs
    cargo run --release --example mk_clump_fixture -- <pairs> <read_len> <out_prefix>

Writes <prefix>_R1.fastq.gz and <prefix>_R2.fastq.gz — gzipped, because the plain-output
branch is a different memory regime (§0). Sequences from a fixed-seed inline LCG; no `rand`
dependency; byte-identical across platforms and runs. Prints pairs, read length and
distinct-sequence count to stderr before writing.

Must clear `cargo clippy --all-targets --release -- -D warnings`, which covers `examples/` —
a hand-rolled LCG is `needless_range_loop` / `cast_possible_truncation` territory.
```

## Implementation outline

Three PRs. D4 puts docs first; reviewer B's ordering finding forbids the CI job from riding
in the docs PR.

### PR 1 — docs only (D4)

1. `clumpy.md`: the eleven corrections plus the six additions in §3.
2. `clump-only.md`: D1's below-floor subsection, and amend line 76's `--cores` parity claim.
   **Placement:** line 75 is a bullet inside the "Composes with:" list under
   `## Compatibility` — inserting a subsection there splits the list. Put it after the list
   or give it its own `## Memory` heading.
3. Leave `reference/changelog.md:198-199` alone (synced page).
4. `npm run build` in `Docs/`, including the `docs-build` job's no-silent-render-corruption
   check. Verify the LaTeX block renders.

### PR 2 — code and guards

1. **D2 banner** (`main.rs:207-213`) per Signature.
2. **Tripwire, arithmetic not labels** (`integration_preflight_tripwire.rs:599`). The case
   runs `["--clumpify", "--cores", "2", …]` at the default `--memory 1G`, which prints
   `predicted peak ≈ 930 of 1024 MiB budget` deterministically (verified). Assert that exact
   substring — it pins both divisors, the margin and the wording at once. Keep the
   `"MiB; predicted peak ≈"` clause for the bins label.
3. **D1 bail context** at the four `clump_only.rs` sites (275/421/764/929):
   `.with_context(|| "--clump_only cannot run within the given --memory")?`. `anyhow::Context`
   is already in use. No fallback, no exit-path change.
4. **Q4 banner tidy**, in this PR rather than later (reviewer A: shipping D2 alone leaves a
   released binary where "budget" means `--memory` on one banner and per-bin bytes on three
   others). `clump_only.rs:289`, `:784`, `:1031` → `N bins × M MiB each`; `:436` already omits
   "budget". Five stderr lines change in total, so Integration's "one user-visible byte
   sequence" claim is corrected accordingly.
5. **4c `channel_depths`** in `clump.rs`; four call sites in `parallel.rs` (`:117`, `:943`,
   `:130`, `:959`); pin both arms with a failure message stating that raising either
   invalidates `k`.
6. **4a pinned-cell table** — six cells, hand-computed.
7. **4b sweep**, renamed for rounding and saturation.
8. **`ClumpLayout::new`** → `pub(crate)`, or assert `bin_byte_budget >= MIN_BIN_BYTES`.
9. **CHANGELOG**: extend #439's existing bullet in `### Unreleased` with one clause naming
   the banner, rather than adding a second bullet for the same issue.

### PR 3 — CI guard

1. `examples/mk_clump_fixture.rs` per Signature.
2. New job `clumpify-memory`, `runs-on: ubuntu-latest`,
   `if: github.event_name != 'schedule'`, **`timeout-minutes: 30`**, `shell: bash` with
   `set -euo pipefail`. Reuse the checkout / toolchain / cargo-cache block but **omit
   `needs:`** — there is no existing job to attach to and the job is meant to run in
   parallel.
3. Generate 10 M pairs × 51 bp gzipped into `$RUNNER_TEMP`; plus the 25 bp and multi-pair
   cells from §8 as logged, non-gating additions.
4. Run each cell under `/usr/bin/time -v` (`command -v` check, or
   `sudo apt-get install -y time`), **`--paired`** explicit, stderr to its own log.
5. Parse `Maximum resident set size (kbytes)` and
   `predicted peak ≈ (\d+) of (\d+) MiB budget`. Assert all captures are non-empty positive
   integers **before** comparing, and compare in **KiB on both sides** to avoid a second
   truncation.
6. Assert, per cell: banner present; output is `.fq.gz`; `peak ≤ prediction`;
   `peak ≥ 0.5 × prediction`; `prediction == budget × 10/11 ± 1 MiB`. Echo both ratios
   whether or not the cell passes.
7. Comment that the measuring parent must stay trivial: on Linux a child's `ru_maxrss` is
   inherited from the forking parent, which produced a full matrix of fictional numbers
   during Phase 2 calibration. `/usr/bin/time` is ~2 MB and safe; a Python or Node wrapper
   holding the fixture is not.
8. Cleanup step. Note the rationale is *this job's own disk*, not later jobs — each GitHub
   job gets a fresh runner VM (reviewer B).

### Issue filing (any time)

The four issues in §7, and the `CALIBRATION.md` item-1 amendment (§9 below).

### 9. `CALIBRATION.md` amendment

Item 1 of "Recommended implementation" recommends FastQC as a `max`. The shipped addend is
correct and the spec is wrong: at `--cores 16` with a 64 MiB pool the process peaks at
438.8 MiB — the trim reservation *plus* ~215 MiB — because mimalloc keeps the pool's pages
resident across the trim→FastQC boundary. Amend in place with a dated note; do not rewrite,
the document is the record of what was believed when the fit was made. Add a second dated
note recording §0: the fitted ratios are gz-path numbers and the plain-output branch is a
different regime.

## Efficiency

- PRs 1 and 2 are text or compile-time. `channel_depths` is a branch on a `bool` evaluated
  four times per run.
- 4a is six integer comparisons; 4b's ~1400 integer cells are trivial.
- **PR 3's cost, corrected after review.** v1 said "~25 s of measurement and ~7 s of
  generation to a job that must build release anyway" — there is no such job;
  `clumpify-memory` is new and builds release itself. On a cold cargo cache (noodles,
  fastqc-rust, zip, mimalloc) that is minutes, hence `timeout-minutes: 30` against v1's 20
  (`rust-tests` allows 45). The cells were timed on macOS/arm64; a 4-vCPU runner will be
  several times slower, especially the FastQC cell. Disk with a gzipped 10 M fixture is
  ~490 MB in and ~600 MB out per cell — well inside the runner's ~14 GB.
- The fixture is 100% distinct. **v1's justification for that is withdrawn**: it claimed
  random sequence maximises the compressed-output term, but reviewer A's 10 M/100%-distinct
  ratios (c2 0.952, c8 0.900) sit essentially on Phase 2's 10 M/52%-distinct ones (0.949,
  0.912), and my own gz-path c4/c8 figures agree. Distinctness is roughly neutral; **pair
  count and output encoding are the axes**. The choice stands as deterministic and cheap, and
  the warning that this fixture must never be reused for a compression claim stays.

## Integration

- **Reads:** nothing new at runtime.
- **Writes:** five changed stderr lines (one `--clumpify` banner, four `--clump_only`
  banners). No output-file bytes change; Perl byte-identity is untouched (`--clumpify` has no
  Perl counterpart).
- **Downstream:** exactly one in-repo assertion depends on the banner (verified by both
  reviewers via grep; the tripwire is a file-set diff test, not an stderr golden).
- **Ordering, corrected.** v1's Integration said "Step 6 can land in either" PR while its
  Self-Review said the CI job must follow the banner change. Both cannot hold: if the job
  lands before D2, its regex matches nothing, and with no non-empty guard that is a **silent
  pass**. Also `ci.yml` is not a docs file. So: **PR 3 lands after PR 2 merges**, never
  inside PR 1.

## Assumptions

1. **Fixed:** no coefficient moves in Phase 3 — `resolve_layout_pins_the_default_cell` must
   pass unmodified. (Narrowed from v1's broader "the model is correct as shipped", which §0
   shows is true only on the gz path.)
2. **Fixed:** `predicted_peak_bytes()` remains the exact inverse of `resolve_layout`. Both
   reviewers verified it symbolically and numerically (within 5–312 B of `budget × 10/11`
   across every cell tried). Note this same exactness is what made v1's sweep vacuous — one
   property, load-bearing in two places with opposite signs.
3. **Fixed:** `--clumpify` requires `--cores ≥ 2`, so `clumpy.md`'s floor range starts there.
4. **Fixed:** `ubuntu-latest` has ≥ 4 usable cores and ≥ 4 GiB RAM — and dropping `2G/c4`
   removes the RAM half. Positive: the job is **insensitive to runner core count**, because
   `--cores` is explicit and FastQC's threads derive from it via `resolved_threads`. The only
   runner-dependent inputs are allocator page behaviour and I/O timing.
5. **Reclassified:** fixture sizing is not a free dial — it determines whether the job
   guards anything (§5). The sizing rationale is the ramp, not a cost preference.
6. **Open:** all ratios here are macOS/arm64. Linux measured ~9 points lower in Phase 2.
7. **Not assumed:** that the plain-output branch is in scope. It is D3's issue.
8. **New:** that `--clumpify` output is gzipped. False for plain input (§0) — which is why
   PR 3 asserts the output extension rather than trusting the banner's `(gzip level 1)`.

## Validation

| # | Verify | How | Expected |
|---|---|---|---|
| 1 | Banner arithmetic, not wording | `--memory 1G` and `--memory 1500M` | `≈ 930 of 1024`; `≈ 1363 of 1500` |
| 2 | No coefficient moved | `cargo test resolve_layout_pins_the_default_cell`, unmodified | passes |
| 3 | Tripwire pins the divisors | change `predicted_peak_bytes()/(1024*1024)` to `/(1000*1000)` | fails (this is what v1's label clause missed) |
| 4 | Channel-depth pin bites | set the clumpy work arm to 2, then the result arm to 4 | the pin fails in both cases |
| 5 | **4a bites where 4b cannot** | perturb σ 23→20, `k` 64→32, static 224→512 MiB, FastQC 24→0, cap 16→1024, one at a time | **the pinned-cell table fails for every one**; 4b is *expected to pass* — that is the point of splitting them |
| 6 | 4b bites on its own property | invert the margin, and make the floor/resolve rounding disagree | 4b fails both |
| 7 | `--clump_only` below the floor | a below-floor cell | non-zero exit, no output file, message names `--clump_only` |
| 8 | `--clump_only` above the floor, incl. high cores and the two uBAM sites | existing suites plus `--fastqc --cores 16 --memory 800M` | all resolve |
| 9 | Docs figures match the binary | run each cell; for the floor-range endpoints, which no invocation prints, reimplement `clumpify_min_memory_bytes` **or** read the top end from `--cores 32 --fastqc --memory 1011M`'s warning | identical |
| 10 | Docs build | `npm run build` in `Docs/` | clean; LaTeX renders |
| 11 | CI guard passes and reports | first run | banner present, output `.fq.gz`, `peak ≤ prediction`, both ratios logged |
| 12 | **CI guard catches the regression it is sold on** | revert lever A (work depth → 2) and run the job | **red**. If it is not, the fixture is still too small — learn it now, not in a year. This replaces v1's `STATIC_TRIM_BYTES → 32 MiB` mutation, which proves only that the plumbing works |
| 13 | CI guard fails loudly on a parse miss | feed it a truncated log | non-zero exit, not a silent pass |
| 14 | Linux baseline | read the ratios off the first runs | ~0.85 / ~0.72 expected; **if any cell exceeds 1.0 on unmodified code, stop and treat it as a finding, not a threshold to tune** |
| 15 | Issue figures before they are public | re-measure the 10 kb cell at ≥ 8× the record count; confirm the plateau; measure gz-output as the headline | plateau ≥ 2 GiB, still climbing |
| 16 | The §0 defect reproduces as written | the three commands in §0 | as shown, before filing |

Validations 12, 13 and 16 are not optional. Four false negatives in one recent session came
from harnesses that could not fail — prove each of these can.

## Questions or ambiguities

- **Resolved (D1–D4):** `--clump_only` hard failure; the banner names the budget; the
  plain-output defect gets its own issue; docs ship first.
- **Resolved (Q2 — CI ratio ceiling):** dissolved. At a steady-state fixture the existing
  `peak ≤ prediction` assertion is itself the drift tripwire, and reviewer A's data kills the
  ceiling independently — the ratio at a stressful fixture spans 0.849–1.065, so 0.92 would be
  red most runs. The runner-independent half of the idea — a **lower** bound — is adopted
  instead.
- **Resolved (Q3, Q4):** both adopted; Q3 as an issue filed in §7 rather than a plan-only
  recommendation, Q4 folded into PR 2 so no release ships the mixed terminology.
- **Resolved (Q1 — stale savings figures):** date them, resolve the pre-existing line-104
  contradiction, re-measure the RRBS row before release. Both reviewers converged on this
  independently.
- **Open — is plain output intended at all in v2?** The `--gzip` deprecation message says
  output is "gzipped by default", which is false for plain input; `--clumpify` rejects
  `--dont_gzip` as pointless while `clump-only.md:78` documents it as supported. D3's issue
  should ask, not assert — the answer decides whether the fix is "always gzip" or "honour
  `--gzip`".
- **Open — does the §0 breach appear on Linux and above `--cores 8`?** All of §0 is
  macOS/arm64 at `--cores` 2/4/8. The issue should say so rather than generalise.

## Self-Review

**What changed from v1, and why.** Both reviewers independently proved v1's two guards
useless, by different methods (A by worked arithmetic on the assertion set, B by replicating
the sizing in Python and running six perturbations). Neither finding is a matter of degree:
the swept test could not fail, and the CI job could not see the regression it was written
for. v1's error in both cases was the same — it reasoned from one measurement instead of a
ramp, and it then promoted the result to a structural claim. v2 replaces both claims with the
ramp and the perturbation matrix.

**The measurement lesson is bigger than either guard.** §0's plain-output confound survived a
calibration, an implementation, two code reviews, two coverage audits and two plan reviews.
It survived because every harness checked the one thing that would catch a fallback to plain
*mode* and nobody questioned the output encoding — while the banner asserted `(gzip level 1)`
over uncompressed bytes. Hence Assumption 8 and the `.fq.gz` assertion in PR 3: the guard
now checks the axis that fooled everyone, not just the axis we already knew about.

**Efficiency.** PR 3 got cheaper and stronger at once: gzipping the fixture cuts it from
2.96 GB to ~490 MB *and* puts it on the branch users run, and dropping `2G/c4` removes both
the least informative cell and the runner-RAM assumption.

**Edge cases.** Carried and now itemised: the banner's unit-token loss under D2 (a guard
regression, closed by the arithmetic assertion), the `--cores 1` remedy being unavailable
under D1, suffix-less `--memory` (unreachable, verified), the FastQC cap boundary at 16/17
threads (now pinned), the parse-miss silent pass, short reads as the known under-charging
axis, multi-pair via the same allocator retention that forced FastQC to be an addend, and
`ClumpLayout::new` as an unguarded route to a meaningless prediction.

**Remaining risks.** (a) All of §0 and §5 is macOS/arm64; the Linux baseline is unmeasured,
and Validation 14 exists because v1's equivalent row had no branch for an intermittent red.
(b) 4a pins six cells, so a coefficient edit that happens to preserve all six passes — the
pins are a tripwire, not a proof. (c) 4c catches a helper edit, not a re-inlined literal.
(d) The stale savings figures stay stale until re-measured; dating them is honest, not a fix.
(e) The §0 defect is filed, not fixed, so `--clumpify` on plain input keeps breaching the
budget until that issue is picked up — which is a deliberate scope decision (D3), not an
oversight.
