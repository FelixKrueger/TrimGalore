# Plan Review B — Phase 3: defend and correct the `--clumpify` memory promise (#439)

**Reviewer:** B (independent; no shared state with Reviewer A)
**Target:** `plans/08162026_clumpify-memory-accounting/phase3/PLAN.md`
**Date:** 2026-08-17
**Baseline verified:** `dev` at `fdade7e`; `target/release/trim_galore` reports
`fdade7e — macos/aarch64 — built 2026-08-17T14:38:36Z`, so the binary matches the tree.

## Verdict

**Implement with the listed changes.**

This is the strongest plan in the #439 arc. Every one of its three corrections to
`SESSION_HANDOFF.md` is right, all twelve of its `file:line` references are exact, and its
central premise about the CI guard (that the prediction is identically `budget × 10/11`) is
algebraically sound. I re-derived all of that independently and it holds.

But the phase exists to *defend* the promise, and **both of its defensive mechanisms fail to
do that as specified**:

- **Step 4's swept property test is blind to every coefficient.** I replicated the sizing
  arithmetic exactly and ran the plan's own assertion set against the perturbations
  Validation #5 names. σ 23→20, `k` 64→32, static 224→512 MiB, and *removing the margin
  entirely* all **pass**. Validation #5 is false on all three of the cases it claims.
- **The CI fixture at 3 M pairs is not at steady state.** I measured all three cells at
  3 M and 6 M pairs: every ratio rises 7–8 points. So §5's "3 M reaches steady state and
  lands within ~3 points of the 10 M-pair calibration matrix" is refuted, and the
  conclusion drawn from it — "structurally cannot be a drift tripwire" — is a property of
  the chosen fixture size, not of the design.

Taken together, Phase 3 as written would ship **no** working guard against the exact
regression it was written to prevent, while stating in three places that it does. The good
news is that the fix is cheap and the plan already contains it: raising the fixture to
~10 M pairs makes the single assertion the plan already specifies catch a lever-A revert
with no threshold tuning, and dissolves open Question 2 along the way.

Third Critical: the public issue in Step 7 would carry an understated headline figure and a
plateau claim that does not exist.

---

## Method

Where I state a number as verified, this is how:

| Method | Used for |
|---|---|
| Ran `target/release/trim_galore` at `fdade7e` into `$TMPDIR` | every banner figure, every floor (read out of the bail/warning message), the `--cores 1` refusal, suffix-less `--memory` |
| Exact Python replication of `clump.rs`'s integer arithmetic, validated against the binary (reproduces `bin_byte_budget = 19_006_356` at 1G/c4 and all 14 measured floors) | the swept-test perturbation matrix, the old-vs-new floor comparison, the algebraic proof below |
| Peak RSS via `os.posix_spawn` + `os.wait4` (`ru_maxrss`, bytes on Darwin), measuring parent kept trivial | fixture-size sweeps, long-read cells |
| `git show 2ab0918:src/clump.rs` | the v2.3.0 formula for the old floors |

Fixture sizes were matched to the plan's: my 10 kb fixture is 400.2 MB against the plan's
400.3 MB, and my 3 M-pair fixture is 493 MB/mate against the plan's 489 MB/mate. I
reproduced the plan's own figures to within 1–4% before changing anything, which is what
makes the divergences below attributable to the variable I changed.

---

## 1. Logic review

### 1.1 The corrections to the handoff — all three verified

| Plan's claim | Verified |
|---|---|
| banner prints `930`, not `931` | ✅ `clumpify: 32 bins × 9 MiB; predicted peak ≈ 930 MiB (gzip level 1)` at 1G/c8; 930 at every `--cores` for a 1 GiB budget |
| floor range is `281–590` / `334–1012`, not `276–590` / `302–1012` | ✅ measured `--clumpify` floors: c2 = 281 / 334, c32 = 590 / 1012 |
| `--clump_only --fastqc --cores 16 --memory 800M` **resolves**; that cell's floor is 775 MiB | ✅ 774M bails advising "≥ 775 MiB", 775M and 800M both resolve |
| unchanged: 9 MiB at 1G/c8, 44 MiB at 4G/c8, `need ≥ 311 MiB` at c6, 1012 MiB at c32 `--fastqc` | ✅ all four reproduced exactly |

The plan is right and the handoff was stale. Note the third row matters more than the plan
says: it means the *only* concrete case Phase 2's deviation note offered for why D1 was
urgent does not exist. D1 is still the right call, but its urgency was overstated upstream
and the plan is correct to say so.

### 1.2 CRITICAL — Step 4's swept test cannot see any coefficient

Validation #5 promises: *"perturb σ (23→20), `k` (64→32), the static (224→512 MiB) one at a
time → at least one assertion fails for each perturbation."*

I built the plan's assertion set — `Err` below the floor / `Ok` at or above it,
`predicted_peak_bytes() ≤ budget`, `bin_byte_budget ≥ MIN_BIN_BYTES`, no overflow — over the
plan's stated input space (`cores ∈ {1,2,3,4,8,16,17,32,1024,65536}` × powers of two to
`u64::MAX/2` × floor±1 × `fastqc` both ways) and ran it against each perturbation:

```
BASELINE (as shipped)                    -> PASS
sigma 23 -> 20                           -> PASS   <-- test is blind
k 64 -> 32                               -> PASS   <-- test is blind
static 224 -> 512 MiB                    -> PASS   <-- test is blind
fastqc per-thread 24 -> 0 MiB            -> PASS   <-- test is blind
fastqc cap 16 -> 1024                    -> PASS   <-- test is blind
margin 10/11 -> 1/1 (margin REMOVED)     -> PASS   <-- test is blind
margin 10/11 -> 11/10 (inverted)         -> FAIL (650 violations)
floor/resolve rounding disagree          -> FAIL (28 violations)
```

The reason is algebraic, not a sweep-coverage gap, so no amount of extra cells fixes it:

```
b         = floor(usable × 16 / d)          ⇒ b·d ≤ usable × 16
dyn       = floor(b × d / 16)               ⇒ dyn ≤ usable
predicted = static + dyn                    ⇒ predicted ≤ static + usable = budgeted
budgeted  = floor(budget × 10/11)           ⇒ predicted ≤ budget   for ANY σ, k, static
```

`predicted_peak_bytes() ≤ budget` is therefore a restatement of `MARGIN_NUM ≤ MARGIN_DEN`.
And `Err below the floor / Ok at it` is a restatement of `clumpify_min_memory_bytes` being
the exact inverse of `resolve_layout` — both sides move together when a constant moves,
because both call the same `dyn_denominator` and `static_bytes_for`.

**Failure scenario.** An implementer writes Step 4 exactly as specified, reaches Validation
#5, perturbs σ, and the suite stays green. Best case they lose half a day and come back
asking what the test was for. Worst case — and this is the realistic one, because the plan
asserts three times that this test "can fail" — they conclude their perturbation was
somehow wrong, mark Validation #5 done, and the repo carries a 1,400-cell test that proves
`10 ≤ 11`. The parent plan's assignment ("retire it or move the assertion somewhere it
bites") would then be recorded as discharged when it is not.

**What the swept test *is* worth keeping for** — and the plan should say this instead: it
pins the *rounding direction* across the whole input space, which is exactly the Phase 2
iteration-#2 bug (`budget / 11 * 10` vs `budget * 10 / 11`). I confirmed it catches that
shape: 28 violations. It also pins saturation behaviour at `u64::MAX`. Both are real. Neither
is coefficient drift.

**What would actually bite.** The only in-repo test that fails under a coefficient
perturbation today is `resolve_layout_pins_the_default_cell` (Phase 2 post-review item 6),
and Phase 2's reviewer B recorded that it plus
`resolve_layout_fastqc_reserves_more_and_shrinks_bins` "fail under every perturbation I
tried, including ±1 MiB on the FastQC term". The correct Step 4 is therefore to **pin more
cells the way `resolve_layout_pins_the_default_cell` pins one** — a small table of
`(budget, LayoutInputs) → exact bin_byte_budget` covering `workers ≠ bin_cores`, the FastQC
cap boundary at 16 and 17 threads, and `--cores 2` and `--cores 32` — *and* keep the sweep
for the rounding/overflow property under an honest name. Two tests, both cheap, and between
them every constant in the file is load-bearing for at least one assertion.

### 1.3 CRITICAL — 3 M pairs is not steady state; the drift-tripwire conclusion is a fixture artefact

Same fixture geometry as the plan's, `--paired --clumpify --cores 4`, macOS/arm64,
single rep, peak RSS via `wait4`:

| cell | 3 M (mine) | 3 M (plan §5) | **6 M (mine)** | 10 M (Phase 2 calibration) |
|---|---|---|---|---|
| `1G/c4` | 0.804 | 0.835 | **0.875** | 0.949 |
| `2G/c4` | 0.810 | 0.855 | **0.888** | 0.976 |
| `1G/c4 --fastqc` | 0.679 | 0.724 | **0.750** | 0.813 |

(peak / printed prediction. My 3 M column reproduces the plan's to within 3–5%, which is
inside the documented 96 MiB macOS within-cell spread — so the 3 M→6 M→10 M trend is the
signal.)

And the ramp below it, `1G/c4`:

```
500 K pairs  -> 351.2 MiB  ratio 0.378
1.5 M pairs  -> 700.8 MiB  ratio 0.754
3 M pairs    -> 748.1 MiB  ratio 0.804
6 M pairs    -> 813.6 MiB  ratio 0.875
10 M pairs   -> (calibration) ratio 0.949
```

Monotone throughout. **No cell is at steady state at 3 M pairs**, and the gap to the 10 M
calibration is 14.5 / 16.6 / 13.4 points — not "~3 points". This also fully explains the
9–12 point discrepancy between §5's table and the parent plan's Phase 2 table, which the
plan does not remark on: it is pair count, not reps and not diversity.

Three consequences, in increasing importance:

1. **§5's headroom is wrong.** At 3 M the guard has 930 − 748 = 182 MiB (19.6%) of slack,
   not the ~154 MiB the plan implies. It is even less sensitive than the plan's own
   pessimistic bullet.

2. **The plan's sensitivity argument is right at 3 M and wrong at 10 M.** "Reverting lever A
   at `1G/c4` adds ~104 MiB: 776 → ~880 MiB, still under the 930 MiB prediction" — true at
   3 M. At 10 M the base is ~883 MiB, so a lever-A revert lands at ~987 MiB and **breaches
   the 930 MiB prediction**. The single assertion the plan already specifies catches it.

3. **So "it is a promise guard, not a drift tripwire" is false as stated.** The plan calls
   this structural: *"buying sensitivity means buying fixture size, without bound."* It is
   not without bound — you need sensitivity at *one* cell, not at every budget, and
   `1G/c4` reaches ~95% of prediction at 10 M pairs. That is ~3.3 GB of transient fixture
   and well under a minute of generation, inside the plan's own 20-minute budget. The
   supporting claim that "sensitivity to lever A scales with `B`, `B` scales with the
   budget" is exactly the reason to size the *smallest* cell properly and drop the largest,
   which is the opposite of the plan's conclusion.

**Recommended reshape of Step 6, which also answers Question 2:**

- Raise `1G/c4` to ~10 M pairs. Read the ratio off two runs. At ~0.95 the existing
  `peak ≤ prediction` assertion *is* the drift tripwire, with no runner-specific baseline
  and no ratio ceiling to tune. Question 2 dissolves.
- **Drop `2G/c4`.** Its pool is 672 MiB (16 × 42), so it needs ~2.3× the data of `1G/c4`
  to saturate; at any affordable fixture size it is the least informative of the three, it
  costs the most RSS, and it is the sole reason Assumption 4's "≥ 4 GiB RAM" is
  load-bearing.
- Keep `1G/c4 --fastqc` (its pool is 240 MiB, so it saturates soonest of the three) and
  expect its ratio to sit lowest — see §4.3, because Validation 11 will fire on it.
- If disk is the binding constraint at 10 M pairs, have the generator emit `.fastq.gz`
  (the reader handles it) and state which form the local sizing run used, so the two agree.

I would not soften §5 — I would delete its third bullet and replace it with the measured
ramp above, because the ramp is the whole argument for the fixture size and it is currently
asserted in the one direction the data does not support.

### 1.4 CRITICAL — Step 7's public issue would carry a wrong headline and a false plateau

The plan re-measured 812 MiB rather than citing 802 "since it is going into a public issue".
I reproduced it, on the same fixture size (400.2 MB vs the plan's 400.3 MB), and then
tripled the input at the same read length:

| read length | input | `--cores 1` | `--cores 4` | `--cores 8` |
|---|---|---|---|---|
| 10 kb | 400 MB | 10.1 MiB | **813.2 MiB** | 817.4 MiB |
| 10 kb | 1.2 GB | 10.1 MiB | **1808.8 MiB** | **2052.3 MiB** |

So:

- **812 MiB is input-limited**, exactly the defect the plan correctly flags for its 50 kb
  cell and then does not apply to the row it is quoting. A 400 MB fixture of 10 kb reads
  holds ~4.7 batches of 4096 records; the plan's own §6 text says "roughly 20
  batch-equivalents are resident". The fixture cannot fill the pipeline it is measuring.
  The honest headline is *at least* 2 GiB, still climbing — 20× the published ~100 MB
  claim, not 8×.
- **"`--cores 8` is 812.2 MiB, i.e. flat in cores past 4" is a fixture artefact.** At 1.2 GB
  c8 exceeds c4 by 244 MiB. The plan hedges with "on this fixture", but the issue text
  would carry a plateau in cores that does not exist, and a plateau claim is the single
  most likely thing a maintainer triaging the issue would act on.
- `--cores 1 → 10.1 MiB` reproduces the plan's 10.6 MiB. That part of the issue is solid,
  and it is the sharpest fact in it: the pool is bypassed entirely.

One more thing to put in the issue, from the parent plan's own resolved question: the
plain-*text* output branch peaked 802 MiB where the gzip branch peaked 508 MiB from an
identical reserve. My runs (and, from the reproduction agreement, the plan's) wrote
uncompressed output. `threading.md`'s claim is about the normal gzip path, so the issue
should lead with the gzip-output figure and quote plain-text as the worst case, or the first
responder will fail to reproduce it.

### 1.5 IMPORTANT — D2 removes the predicted peak's own unit label, undoing Phase 2 post-review fix #4

Phase 2's post-review deliberately strengthened `integration_preflight_tripwire.rs:599`
because *"the original assertion's two clauses both keyed on the bins label, so reverting the
peak label alone would have passed."*

After D2 the banner is:

```
clumpify: 16 bins × 22 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
```

The predicted peak now has **no unit token of its own** — the `MiB` belongs arithmetically to
`memory_bytes / 1024²`, a different quantity. Step 1.2's replacement clause
`"MiB budget (gzip level"` pins the *budget's* divisor, not the peak's. Change
`predicted_peak_bytes() / (1024*1024)` to `/ (1000*1000)` and the banner reads
`predicted peak ≈ 975 of 1024 MiB budget` — **both clauses still pass**. That is the exact
class of slip that produced the 4.86% harness error and the original `MB`/`MiB` bug, and
Phase 3 would re-open it while calling the change guard-neutral.

This is not an argument against D2. The fix is one line and strictly stronger than either
clause: the tripwire case runs `["--clumpify", "--cores", "2", "a.fastq"]` with the default
`--memory 1G`, and I verified 1G/c2 prints `predicted peak ≈ 930 MiB` deterministically. So
assert the arithmetic:

```rust
out.stderr.contains("predicted peak ≈ 930 of 1024 MiB budget")
```

That pins both divisors, the margin, and the new wording in one substring. Keep the
`"MiB; predicted peak ≈"` clause for the bins label.

### 1.6 IMPORTANT — the `--clump_only` regression table has two wrong cells

Every other floor in the plan checks out. These two do not:

| `--cores` | plan "new floor" | measured `--fastqc` floor | measured plain floor |
|---:|---:|---:|---:|
| 2 | 277 | **329** | 277 |
| 10 | 578 | **579** | 315 |

The `--cores 2` row is the **plain** floor dropped into a table whose other five rows are
`--fastqc` floors. The `--cores 10` row is off by one because the plan truncated where the
binary uses `div_ceil` (and the same row's "old floor 580" is a `div_ceil` of an old
arithmetic that truncated to 579 — the two columns use different rounding).

Neither changes a verdict: at c2 the FastQC floor still drops 535 → 329, and at c10 the
exact floors are 579.50 → 578.05, so still no regression. But the plan opens the section
with "Every number below was re-derived by running `target/release/trim_galore`", and these
two were not — which matters because this table is the evidence that makes D1 cheap.

The crossover claim is exactly right, and I confirmed it across `--cores` 1–64: c9 is
572.75 → 545.33, c10 is 579.50 → 578.05, c11 is 586.25 → 610.78. **`≥ 11` is the precise
boundary.** Two small fixes to the surrounding prose: the widest window is **154.4 MiB at
`--cores 16`**, not "≤ 150 MiB", and the c11 window in whole MiB is `[587, 611)` — 586 MiB
did not resolve under the old formula either.

### 1.7 IMPORTANT — Step 3 guards one half of `k` and leaves the other unguarded

`clump.rs:247-248` states what `k = 4` pays for: *"the queued batch, the batch in hand, and
the compressed output waiting on the result channel"*. Step 3 pins the first of those
(`queue_depth`, `parallel.rs:117`/`:943`). The third is `sync_channel(cores * 2)` at
`parallel.rs:130` and `:959`, and it has no guard at all — and the parent plan's Phase 1
post-review names it explicitly as the *next* thing someone will touch:

> **Next structural lever for the `2G/c4` residual** (reviewer A, measured): the result
> channel `sync_channel(cores * 2)` plus the `pending` map hold *written*, therefore
> resident, compressed buffers … Cheaper to try than buffer recycling.

**Failure scenario.** Someone picks up the `2G/c4` residual, raises the result channel to
`cores * 4` because it is the cheapest lever on the list, gains throughput, and every
deterministic guard in the repo stays green — `work_queue_depth` is untouched, the swept
test is blind (§1.2), and the CI job at 3 M pairs has 182 MiB of slack (§1.3). `k` is now
wrong and nothing says so.

Fix: fold the result-channel depth into the same helper (`clump::result_queue_depth`, or a
single `work_channel_depths(clumpy)` returning both), pin both arms, and put the
"raising either invalidates `k`" note on the helper next to `dyn_denominator` as the plan
already intends. Same cost, twice the coverage, and it satisfies the plan's own stated
principle of putting the invariant beside the constant it constrains.

### 1.8 IMPORTANT — the ordering rule in Integration contradicts the Self-Review

Self-Review: *"Step 6 must follow Step 1 — stated in the outline rather than left implicit."*
Integration: *"Steps 5 and 7's docs edit can land as a separate docs-only PR — and should…
**Step 6 can land in either.**"*

Both cannot hold. If Step 6 lands in the docs PR before the code PR merges, the job's regex
`predicted peak ≈ (\d+) of (\d+) MiB budget` matches nothing on that branch — and with no
non-empty guard on the parse (§2.2) that is a *silent pass*, not a failure. Separately,
`.github/workflows/ci.yml` is not a docs file, so putting Step 6 in the "docs-only PR"
violates the same repo rule the plan invokes one clause earlier.

State it once: **Step 6 lands with Steps 1–4, or in its own third PR after they merge.**

### 1.9 IMPORTANT — dropped scope

Against `SESSION_HANDOFF.md` §2, six of seven rows are discharged (D1, D2, CI guard, docs
arithmetic, threading issue, `CALIBRATION.md`). The seventh is **not addressed anywhere** —
not in Goal, not in the explicit non-goals, not in Questions, not in Remaining Risks:

> **Unverified axes** — single-end, multi-pair (promise is per-process, all cells were
> single-pair), short reads (σ ≈ 1.55 on 25 bp vs coefficient 1.4375), low-duplication/WGBS
> data.

Two of those are load-bearing for a phase whose subject is defending the promise:

- **Short reads.** σ ≈ 1.55 at 25 bp against a shipped 1.4375 means there is a *known* read
  length at which the model under-reserves. Declining to change σ is a fine non-goal;
  silently dropping the *verification* item is different. And the CI fixture is 51 bp —
  the one length the model was fitted on — so the guard is structurally unable to see it.
  One extra generated read length is nearly free in the same job.
- **Multi-pair.** This is the sharpest one, and it follows from Phase 2's own central
  discovery. FastQC is an addend rather than a `max` *because mimalloc keeps the bin pool's
  pages resident across a phase boundary*. Bins are re-allocated per pair, so by the same
  mechanism pair 2's pool may sit on top of pair 1's retained pages — which would breach the
  promise on any multi-pair run, and nothing in the repo or the plan tests it. One extra CI
  cell (the same fixture passed twice as two pairs) settles it.

Also dropped from the parent plan's Phase 1 "Logged, not fixed" list, without mention:

- **`ClumpLayout::new` is public, has no production caller, and validates neither invariant
  `resolve_layout` enforces.** This belongs in a "defend the promise" phase: any future
  caller can construct a layout whose `predicted_peak_bytes()` is arbitrary, and the CI
  guard — which parses the banner — would faithfully report the lie and pass. Steps 3 and 4
  already open `clump.rs`. Make it `pub(crate)`/`#[cfg(test)]`, or have it assert
  `bin_byte_budget >= MIN_BIN_BYTES`.
- **Reviewer B's pre-existing bug: the EOF bin flush swallows a failed send.** That is
  silently dropped reads — a data-integrity defect. Phase 3 files an issue for a wrong
  sentence in the docs and not for this. Under the repo's own issue-filing bar
  ("data-integrity/user-visible defects"), the flush bug clears it more easily than
  `threading.md` does. File both, or say why not.

---

## 2. Assumptions

### 2.1 Stated assumptions — checked

| # | Verdict |
|---|---|
| 1 (model correct as shipped; `resolve_layout_pins_the_default_cell` must pass unmodified) | ✅ Sound, and it is currently the *only* test that fails under a coefficient perturbation — see §1.2. Worth saying so in the plan, since it makes the assumption load-bearing rather than incidental. |
| 2 (`predicted_peak_bytes()` remains the exact inverse) | ✅ Verified algebraically and numerically: prediction = `floor(budget × 10/11)` to within 5–154 bytes across 1G/2G/4G/1500M × c4/c8/c32 × `fastqc`. This premise is the plan's strongest and I could not break it. |
| 3 (`--clumpify` requires `--cores ≥ 2`) | ✅ `--clumpify --cores 1` errors: *"--clumpify requires --cores >= 2 (the bin dispatcher feeds parallel workers)"*. |
| 4 (runner ≥ 4 cores, ≥ 4 GiB) | ✅ Holds for public-repo `ubuntu-latest` (4 vCPU / 16 GiB). **And it is the only reason `2G/c4` matters** — dropping that cell (§1.3) removes the assumption. Positive note: the job is *insensitive* to runner core count, because `--cores` is explicit and FastQC's thread count derives from it via `resolved_threads(None, cores)`. Worth stating; it is a genuine strength of the design. |
| 5 (cells and fixture size are dials, not correctness) | ⚠️ True in principle, but §1.3 shows the fixture size *is* what determines whether the job guards anything. Reclassify as "sizing determines coverage; the sizing rationale is in §5". |
| 6 (macOS ratios; Linux will differ) | ✅ Correct and correctly flagged. Linux measured ~9 points lower in Phase 2, so at 3 M pairs expect ~0.71 / 0.72 / **0.59** — see §4.3. |
| 7 (not assumed that CI covers lever A) | ⚠️ The *claim* is honest but the *reason* given is wrong (§1.3). Keep the caveat, replace the reasoning. |

### 2.2 Unstated assumptions

- **That the log parse succeeds.** Nothing in Step 6 guards against an empty capture. If
  `Maximum resident set size` or the banner regex misses — a `/usr/bin/time` version change,
  a banner reword, Step 6 landing before Step 1 (§1.8) — a shell comparison on an empty
  string either errors or, without `set -euo pipefail`, compares as 0 and **passes**. The
  plan's own Validation-10 preamble names this failure family ("four false negatives in one
  recent session came from harnesses that could not fail"); the job needs the same
  discipline applied to itself. Assert both captures are non-empty positive integers before
  comparing, and set `shell: bash` + `set -euo pipefail`.
- **That the binary cannot move its own goalpost.** Step 6.4 sells parsing the budget off
  the banner as a virtue: *"the job never restates a number the binary owns."* But then the
  binary owns *both* sides of the inequality. Since prediction ≡ `budget × 10/11`
  identically, removing `MARGIN_NUM/MARGIN_DEN` inflates the prediction to `budget`, real
  RSS grows ~10%, and `peak ≤ prediction` still passes. Add the one-line closure:
  `prediction == budget × 10/11` (± 1 MiB). Two integers already in hand; no new
  measurement.
- **That `/usr/bin/time` exists on `ubuntu-latest`.** It does today. Failure mode is loud,
  so this is cheap risk — but note "Maximum resident set size (kbytes)" is KiB, and doing
  the comparison in KiB on both sides (`pred_mib * 1024`) avoids a second truncation, in a
  codebase where a unit slip in this exact arithmetic has already cost a calibration run.
- **That each job's disk is shared.** Step 6.7 justifies cleanup "so the runner's disk is
  not left at +1.2 GB for later jobs in the same workflow". Each GitHub job gets a fresh
  runner VM; nothing is shared. Harmless step, wrong reason — and the reason is the kind of
  thing that grows into `if: always()` scaffolding for nothing.
- **That the `--memory` default is 1G in the tripwire case.** True (`cli.rs:249`), and it is
  what makes the exact-substring assertion in §1.5 available. Worth stating, since it is the
  premise of the strongest available guard.
- **That the fixture generator is held only to "it compiles".** `cargo clippy --all-targets
  --release -- -D warnings` covers `examples/`. A hand-rolled LCG is prime
  `needless_range_loop`/`cast_possible_truncation` territory.
- **That `--paired` is passed.** Step 6.3 names the cells (`1G/c4` etc.) but never the flag;
  §5's numbers are paired, and a single-end run has one bin pool and half the data.

Two assumption-level things the plan gets right and deserves credit for: the suffix-less
`--memory` edge case is correctly dispositioned (I verified `--memory 1500` warns and falls
back to plain mode *before* the banner, so no guard is needed and the plan is right to say
so explicitly), and Validation 10 genuinely works — perturbing `STATIC_TRIM_BYTES` to 32 MiB
leaves the prediction at 930 MiB while raising real RSS by ~192 MiB, so `peak > prediction`
by construction. That is the right choice of perturbation and the plan is right that it is
not optional.

---

## 3. Efficiency

Steps 1–4 and 7–9 are text or compile-time; `work_queue_depth` is a branch on a `bool`
executed twice per run. All fine and correctly assessed.

Two overstatements in the Efficiency section, both harmless but worth fixing so nobody
budgets from them:

- *"adds ~25 s of measurement and ~7 s of generation to a job that must build release
  anyway"* — there is no existing job it attaches to; `clumpify-memory` is new and must
  build release itself. On a cold cargo cache (noodles, fastqc-rust, zip, mimalloc) that is
  minutes, and `timeout-minutes: 20` is tight — `rust-tests` allows 45 for a comparable
  shape. And the three cells were timed on macOS/arm64 (7–9.5 s each); a 4-vCPU
  `ubuntu-latest` will be several times that, especially the FastQC cell.
- *"~1.2 GB of transient disk"* counts the fixture but not the three cells' outputs
  (~250 MB of `.fq.gz` per cell, plus FastQC artefacts). Real figure is ~2 GB at 3 M pairs,
  ~5–6 GB if the fixture goes to 10 M — still fine against ~14 GB free alongside a release
  `target`, but gzipping the fixture is the cheap lever if it gets tight.

The Efficiency note on fixture *diversity* is the best-reasoned paragraph in the section and
I agree with it: 100% distinct is conservative for a memory guard, and the warning not to
reuse it for a compression claim is exactly right. My measurements confirm diversity is not
what drives the ratio — pair count is (§1.3).

The swept test's ~1,400 integer cells are trivially cheap. Once it is renamed to what it
actually proves, the sweep to `u64::MAX / 2` is worth keeping: saturation is the one place
the arithmetic could stop being conservative, and it is the cheapest possible check on it.

---

## 4. Validation sufficiency

### 4.1 Where an implementation could follow this plan exactly and still ship something wrong

1. **Coefficient drift ships unguarded** (§1.2 + §1.3). This is the phase's entire purpose.
   Rows 5 and 9–11 are the only ones aimed at it; row 5 does not fire, and rows 9–11 pass at
   3 M pairs with 182 MiB of slack.
2. **A peak-figure unit regression ships unguarded** (§1.5). Validation 3 ("revert the banner
   by hand, run `integration_preflight_tripwire` → fails") is satisfied by reverting the
   *wording*; it does not exercise the divisor.
3. **The CI job silently passes on a parse miss** (§2.2). Validation 10 proves the job can
   fail *when the numbers are there*. It does not prove it fails when they are not — and
   "Step 6 can land in either PR" (§1.8) makes that reachable on the first run.
4. **A multi-pair breach ships unguarded** (§1.9), through the same allocator-retention
   mechanism Phase 2 discovered.
5. **The result-channel depth moves and `k` silently stops applying** (§1.7).

### 4.2 Rows that are well-designed

Row 10 (prove the guard can fail, with `STATIC_TRIM_BYTES` → 32 MiB) is the right
perturbation for the right reason, and I verified it works. Row 2 is the correct
tripwire-on-scope-creep. Row 7 is worth having and cheap — and I confirmed it is safe:
no test in `src/clump_only.rs` or either `tests/integration_clump_only*.rs` asserts on the
bail text, and `small_memory_budget()` is 1 GiB, comfortably above every `workers: 1` floor,
so Step 2's `.with_context()` cannot break an existing test. Row 12 (docs build) is the sort
of row plans forget.

### 4.3 Row 11 will fire, and the plan should expect it

Row 11: *"within a few points of the macOS ratios in §5; if ≪ 0.7, raise the pair count."*
Linux measured ~9 points below macOS in Phase 2. At 3 M pairs on the runner, expect
~0.71 / 0.72 / **0.59**. The FastQC cell trips the stated threshold. That is Row 11 working
correctly — but the plan reads as though a trip would be a surprise, and Row 11 is a *human
eyeball on a log line*, not an assertion.

Make it an assertion. A lower bound (`peak ≥ 0.5 × prediction`) is strictly better than the
ratio ceiling of Question 2 for a first pass: it catches a silent fallback to plain mode
*and* a fixture that stopped filling the pool, and unlike a ceiling it is runner-independent,
because the failure it detects is a 3× miss rather than a 5% drift. Land the floor now;
defer the ceiling, or let §1.3's 10 M fixture make it unnecessary.

### 4.4 Missing rows

- Banner arithmetic, not wording: `--memory 1G` ⇒ `≈ 930 of 1024`; `--memory 1500M` ⇒
  `≈ 1363 of 1500` (both verified against the binary today).
- Job fails loudly on an unparseable log (feed it a truncated log; expect non-zero).
- `--clump_only` above the floor at high cores — Row 7 covers the default `--cores 2`
  path, but D1's whole surface is `--fastqc --cores ≥ 11`, and Step 2 touches all four
  `resolve_layout` sites including the two uBAM ones (`clump_only.rs:764`, `:929`), which
  no listed validation reaches.

---

## 5. Alternatives worth considering

- **For Step 4: pin cells, don't sweep properties.** Detailed in §1.2. `resolve_layout_pins_
  the_default_cell` is proof that the pinning style bites; the sweep style provably does not.
  Keeping both, honestly named, is strictly better than replacing one tautology with another.
- **For Step 6: one steady-state cell beats three under-filled ones.** Detailed in §1.3.
  Same wall time, same disk order of magnitude, and it converts the job from a promise guard
  into the drift guard the plan says is impossible.
- **For the margin: assert it in CI for free.** `prediction == budget × 10/11` closes the
  self-referential loop (§2.2) with two integers already parsed.
- **For §1.7: one helper returning both channel depths**, rather than one helper for the one
  lever that already moved.
- **Rejected alternative — a `--memory`-aware batch size for the plain path.** Tempting
  given §1.4's numbers, and correctly out of scope: it is the path every mode shares, and
  the plan is right that it needs its own plan. Filing the issue is the correct move; my
  objection is only to the figures in it.

---

## 6. The four open questions — my recommendations

**Q1 — stale savings figures in `clumpy.md`.** *Plan: leave the numbers; re-measure before
the next release.* **Partly disagree.** Two facts make "leave it" the inconsistent choice in
a PR whose subject is memory honesty. First, the RRBS row is the only row advertising a
*memory-dependent* benefit ("~24% at default `1G`; up to ~31% at `4G+`"), and the 4G/c8 bin
moved 66 → 44 MiB — a 33% cut in sort-run size, i.e. precisely the input to that claim.
Second, line 104 ("increasing memory doesn't seem to make very much difference") already
**contradicts** the RRBS row in the same file, and did so before Phase 2 — so there is a
pre-existing internal inconsistency the arithmetic pass is walking past.

Recommendation: don't fabricate numbers, but *date* them. One clause on the RRBS row —
"measured before the v2.4 memory refit; re-measurement pending" — and resolve the line-104
contradiction while you are in the file. Dating a claim is not measuring it, costs one
clause, and is the difference between "known stale" and "published wrong". Then re-measure
before release as the plan says.

**Q2 — a ratio ceiling in CI.** *Plan: land without it, observe 2–3 runs, then add with ≥2×
the observed spread.* **Agree on sequencing, but the question mostly dissolves.** With a
steady-state `1G/c4` cell (§1.3) the headroom is ~5% and a lever-A revert breaches the
existing assertion outright — no threshold, no runner-specific baseline, nothing to tune.
What I *would* add now is the **lower** bound from §4.3, which is the runner-independent half
of the same idea and turns Row 11 into a guard instead of a log line. If the fixture stays
at 3 M, then the plan's staged approach is the right one and 2× the observed spread is the
right margin — but be explicit that a ceiling on an under-filled fixture is calibrating
against a moving baseline.

**Q3 — should `--clump_only` print a predicted peak?** *Plan: yes, but a separate small PR.*
**Agree, with one addition: file it as an issue in this phase.** The scope argument is
sound. But note it is the natural companion to D1 rather than an unrelated nicety: D1 makes
below-floor a hard failure whose remedy is "raise `--memory` / lower `--cores`", and a user
has no way to see how close they are, which is exactly what makes the abort feel abrupt.
Step 7 already files one issue; file this as a second one-liner. A recommendation living
only in a plan document is the shape of thing that rots.

**Q4 — "budget" terminology collision.** *Plan: drop "budget" from the three `--clump_only`
banners so all four read `N bins × M MiB`.* **Agree.** Verified the asymmetry: `clump_only.rs`
`:289`, `:784`, `:1031` say "budget"; `:436` does not. Two refinements. (a) If you are
touching four lines anyway, `N bins × M MiB each` costs five characters and removes the
per-bin/total ambiguity that both spellings share today. (b) It makes the Integration claim
"one user-visible byte sequence changes" into five — either fold that into the claim or
defer Q4, but don't leave the claim standing beside a change that contradicts it.

**Consequence of D2 the plan does not carry (not a re-litigation).** Phase 2's reviewer B
predicted this exact banner: *"the user will see a predicted peak below the budget, which
will read as 'why didn't it use the 1G I gave it?'. One clause — '(9% held back for
run-to-run variation)' — turns a support question into a documented choice."* D2 names the
budget but not why the two numbers differ, and Step 5 does not give the margin a documented
home either — it corrects the formula to `× 10/11` and leaves the reader to infer what that
factor is. One sentence in `clumpy.md`'s memory section naming the ~9% as deliberate
run-to-run headroom discharges the consequence without touching D2.

---

## 7. Action items

### Critical

- **C1 — Redesign Step 4.** The swept property test cannot see σ, `k`, the static
  reservation, the FastQC term, the cap, or margin removal (verified by exact replication;
  proof in §1.2). Replace with (a) a pinned-cell table in the style of
  `resolve_layout_pins_the_default_cell`, covering `workers ≠ bin_cores`, the FastQC cap at
  16/17 threads, and `--cores 2`/`32`; plus (b) the sweep retained under an honest name for
  the rounding/saturation property it *does* prove. **Rewrite Validation #5**: as written it
  cannot be satisfied.
- **C2 — Resize the CI fixture and rewrite §5's third bullet.** 3 M pairs is not steady
  state; measured ratios rise 7–8 points at 6 M and 13–17 points to the 10 M calibration
  (§1.3). Raise `1G/c4` to ~10 M pairs, **drop `2G/c4`**, keep `1G/c4 --fastqc`. Replace
  "structurally cannot be a drift tripwire" with the measured ramp — at ~0.95 of prediction
  the assertion the plan already specifies catches a lever-A revert, and Question 2
  dissolves.
- **C3 — Re-measure before filing Step 7's issue.** 812 MiB at 10 kb is input-limited: same
  400 MB fixture reproduces 813 MiB, but 1.2 GB gives 1809 MiB at c4 and 2052 MiB at c8
  (§1.4). Drop the "flat in cores past 4" claim — it is a fixture artefact. State the
  headline as ≥ 2 GiB and still climbing, quote the gzip-output figure as the headline with
  plain-text as the worst case, and apply the 50 kb caveat to the whole table.

### Important

- **I1** — Replace Step 1.2's label clause with the arithmetic assertion
  `contains("predicted peak ≈ 930 of 1024 MiB budget")`. D2 removes the peak figure's own
  unit token, re-opening the divisor bug Phase 2 post-review fix #4 closed (§1.5).
- **I2** — Harden Step 6's parse: `shell: bash` + `set -euo pipefail`, assert both captures
  are non-empty positive integers, compare in KiB on both sides, and add
  `prediction == budget × 10/11 (± 1 MiB)` so the binary cannot move its own goalpost
  (§2.2). Add a validation row that feeds the job a truncated log and expects non-zero.
- **I3** — Fix two cells in the `--clump_only` regression table: `--cores 2` new
  `--fastqc` floor is **329 MiB** (277 is the plain floor), `--cores 10` is **579** not 578
  (§1.6). Verdicts unchanged; the "every number re-derived" claim is not.
- **I4** — Extend Step 3 to the result-channel depth (`parallel.rs:130`, `:959`). `k = 4`
  charges for it per `clump.rs:247-248`, and the parent plan names it as the next lever
  someone will pull (§1.7).
- **I5** — Resolve the ordering contradiction: Step 6 lands with Steps 1–4 or in its own PR
  after them, never in the docs-only PR (§1.8).
- **I6** — Give the docs floor range a stated core span. Measured `--clumpify` floors at
  `--cores 64` are **933 MiB** plain / **1356 MiB** with `--fastqc` — above the 1 GiB
  default — while the sentence being corrected says "any sensible core count".
- **I7** — Amend `clump-only.md:76` ("`--cores <N>` accepted for interface parity"), which
  Step 2.2's new paragraph directly contradicts; `--cores` sets `n_bins` and therefore moves
  the floor (measured 277 MiB at c2 → 454 MiB at c32).
- **I8** — Step 7's issue must name **three** batch constants. `src/parallel.rs:33
  BATCH_SIZE = 4096` is the one that sizes the round-robin batches in the plain path being
  measured (`parallel.rs:683`); the plan lists only `fastq.rs:220` and `bam.rs:38`.
- **I9** — Address or explicitly defer `SESSION_HANDOFF` §2's "Unverified axes" row. Two
  cheap CI cells cover the sharpest halves: a second read length (the fixture is 51 bp, the
  one length σ was fitted on, and σ ≈ 1.55 at 25 bp is a known under-reserve) and a
  multi-pair run (same mimalloc retention that forced FastQC to be an addend) (§1.9).
- **I10** — Carry forward two dropped items from the parent plan's "Logged, not fixed":
  `ClumpLayout::new` is public and validates neither invariant (squarely in "defend the
  promise", and Steps 3–4 already open the file), and the swallowed send on the EOF bin
  flush is a data-integrity defect that clears the repo's issue-filing bar more easily than
  the `threading.md` sentence does (§1.9).
- **I11** — Step 5 should say explicitly that
  `Docs/src/content/docs/reference/changelog.md:194-199` — which cites the 512 MiB
  reservation and "within 1%" — is left alone, per both Phase 2 reviewers and the
  synced-changelog rule. It is the only other 512 MiB reference in `Docs/`, so a `grep`
  sweep will find it.
- **I12** — Give the 9% margin a documented home in Step 5. Step 5 corrects the formula to
  `× 10/11` without ever saying what the factor is, while D2 puts two numbers 9% apart on
  stderr — the support question Phase 2's reviewer B predicted verbatim (§6).

### Optional

- **O1** — Widest regression window is **154.4 MiB at `--cores 16`**, not "≤ 150 MiB"; the
  c11 window in whole MiB is `[587, 611)`, since 586 MiB did not resolve under the old
  formula either.
- **O2** — Step 6.7's cleanup rationale is wrong: each job gets a fresh runner VM, so
  nothing is shared with later jobs. Keep the step, drop the reason.
- **O3** — Efficiency understates the job: new job, cold-cache release build, macOS-timed
  cells on a 4-vCPU runner, and ~2 GB of disk rather than 1.2 GB. Consider
  `timeout-minutes: 30`.
- **O4** — `examples/mk_clump_fixture.rs` is inside `cargo clippy --all-targets --release
  -- -D warnings`; a hand-rolled LCG needs to clear that bar.
- **O5** — Step 6.3 does not say `--paired`. All of §5's cells are paired.
- **O6** — On Q4, `N bins × M MiB each` costs five characters and removes the
  per-bin/total ambiguity; and note the Integration claim becomes five changed stderr lines,
  not one.
- **O7** — Worth stating as a strength: the job is insensitive to runner core count, because
  `--cores` is explicit and FastQC's threads derive from it via `resolved_threads`. The only
  runner-dependent inputs are allocator page behaviour and I/O timing.
- **O8** — `clumpy.md:110` ("prints a warning and falls back to plain mode") now needs a
  one-clause `--clump_only` exception, since D1 makes that path a hard failure and the two
  pages cross-link at line 11.
- **O9** — If the fixture grows per C2, have the generator emit `.fastq.gz` and state which
  form the local sizing run used, so CI and the sizing measurement agree.

---

## 8. What I could not fault

Recorded so the plan gets credit where the evidence supports it, rather than a manufactured
finding per section:

- All three handoff corrections are right, and I reproduced each against the binary.
- All twelve `file:line` references check out exactly: `main.rs:207-213`, `parallel.rs:117`
  and `:943`, `clump_only.rs:275`/`421`/`764`/`929`,
  `integration_preflight_tripwire.rs:599`, `clumpy.md` lines 91/96/97/98/102/108/113/119,
  `threading.md:61` and `:63`, `clump-only.md` ≈75.
- "Exactly one in-repo assertion depends on the banner" — verified by grep across `src/` and
  `tests/`; the tripwire is a file-set diff test, not an stderr golden, so nothing else
  moves.
- `predicted_peak_bytes()` is the exact inverse and prediction ≡ `floor(budget × 10/11)` to
  within 154 bytes across every cell I tried. The single-assertion premise is sound.
- The suffix-less `--memory` edge case is correctly dispositioned as unreachable — verified,
  `--memory 1500` warns and falls back before the banner.
- Validation 10's perturbation genuinely works, for the structural reason that the
  prediction is static-independent while real cost is not.
- `examples/` autodiscovery needs no `[[example]]` stanza (`mk_bam_fixture.rs` proves it).
- Step 2's `.with_context()` cannot break an existing test: no assertion anywhere touches
  the `clump_only` bail text, and `small_memory_budget()` is 1 GiB.
- The `≥ 11` crossover for the `--clump_only --fastqc` regression is exactly right across
  `--cores` 1–64, and "without `--fastqc` the floor dropped at every core count" is right
  too (620 → 352 MiB at c16 confirmed; the drop widens with cores).
- The Efficiency note on why 100%-distinct sequence is the conservative fixture for a
  *memory* guard, and must not be reused for a compression claim, is correct and is the kind
  of caveat that usually goes unsaid.
