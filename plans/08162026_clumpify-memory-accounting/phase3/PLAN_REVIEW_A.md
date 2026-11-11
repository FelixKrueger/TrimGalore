# Plan review A — Phase 3: defend and correct the `--clumpify` memory promise (#439)

**Reviewer:** A (independent; no shared state with Reviewer B)
**Target:** `plans/08162026_clumpify-memory-accounting/phase3/PLAN.md`
**Baseline reviewed:** `dev` at `fdade7e`, release binary `target/release/trim_galore` (`fdade7e — macos/aarch64 — built 2026-08-17T14:38:36Z`)
**Date:** 2026-08-17

## Verdict

**Implement with the listed changes.**

Three of the four things this plan does are ready: D1 and D2 are clean, the docs
arithmetic is correct down to the line number, and every correction the plan makes to
`SESSION_HANDOFF.md` is right — I re-derived each one independently and each held.

The two guards are not ready, and both fail in the same direction: they look like
protection and are not.

- **The CI fixture is undersized by enough to invert §5's own conclusion.** I reproduced
  the plan's 3 M-pair numbers to within 0.4 MiB, then ran the identical generator at
  10 M pairs: `1G/c4` peak went from 776.6 MiB (ratio 0.835) to a 6-rep max of
  **990.9 MiB — ratio 1.065, i.e. above the printed prediction on unmodified code.**
  So (a) the gate as specified has ~17 points of unearned slack at 3 M, and (b) at a
  properly-sized fixture the gate is red on `fdade7e` today. Both halves of §5's
  design rest on the smaller number.
- **Step 4's replacement test cannot fail under the perturbations Validation #5 names.**
  All three of its assertions are invariant under σ, `k` and static changes, because the
  sizing, the prediction and the floor all read the same two helpers. It is exactly as
  tautological as the test it retires.

Neither is a reason to re-plan — the fixes are small and local, and everything else in
the phase can proceed unchanged.

## What I verified, and how

Method stated per item, because several of these figures are load-bearing.

| Claim | Method | Result |
|---|---|---|
| Banner prints `≈ 930 MiB`, not 931, at `--memory 1G` | ran the binary at `1G × c2/c4/c8` | **930** ✓ plan right, handoff stale |
| `32 bins × 9 MiB` at `1G/c8`; `32 bins × 44 MiB` at `4G/c8` | ran the binary | both ✓ |
| Prediction is `budget × 10/11` for *every* 1 GiB cell | ran `c2/c4/c8` ± `--fastqc`; re-derived symbolically | ✓ 930 in all six; prediction ∈ (budget×10/11 − denom/16, budget×10/11], i.e. within 312 B at ≤ c32 |
| Floor range `281–590` / `334–1012` MiB over cores 2–32 | reimplemented `clumpify_min_memory_bytes` in Python from `clump.rs` | ✓ exactly (280.50→281, 589.60→590, 333.30→334, 1012.00) |
| `276`/`302` are the `--cores 1` values | same | ✓ (276.10 / 302.50) |
| `need ≥ 311 MiB` at `--cores 6` | ran `--clumpify --cores 6 --memory 100M` | ✓ verbatim |
| `--clump_only --fastqc --cores 16 --memory 800M` resolves; floor 775 MiB | ran both 800M and 774M | ✓ 800M runs (`64 bins × 1 MiB`); 774M exits 1 advising `≥ 775 MiB` |
| Handoff's `~841 MiB` was the pre-`workers: 1` figure | recomputed with `workers: 16` | ✓ 840.4 MiB — the plan's diagnosis is right |
| All 11 `clumpy.md` line references | read the file | ✓ all 11 exact (91, 96, 96–97, 98×2, 102×2, 108×3, 113), plus 119 |
| `threading.md:61` and `:63` | read the file | ✓ both exact |
| Four `resolve_layout` sites at 275/421/764/929; tripwire at 599; `main.rs:207-213` | grep | ✓ all exact |
| `tests/integration_clump_only.rs` exists; `examples/` autodiscovery | `ls`, `grep Cargo.toml` | ✓ file exists; no `[[example]]` stanzas, `mk_bam_fixture.rs` present |
| §5's three cells at 3 M pairs | own generator (51 bp, 56-char headers, 100% distinct), `os.wait4` | `1G/c4` **776.6** MiB vs plan's 776.2 — reproduces exactly; `2G/c4` 1452.5 vs 1590.9 (0.781 vs 0.855); `1G/c4 --fastqc` 633.8 vs 673.0 |
| §6 read-length table | own fixtures at the plan's byte sizes, `os.wait4` | 150 bp **77.8** (plan 71.1); 1 kb **248.8** (240.0); 10 kb **768.5/728.3/777.8** (812.0); c1 **10.2** (10.6); c8 **783.0** (812.2) — all reproduce within rep spread |
| `--clump_only` regression window: first at c11, plain floor dropped everywhere | swept cores 1–256, both formulas | ✓ first window at **c11**; plain floor never exceeds the old floor |

Two harness notes so the numbers can be re-derived: `/usr/bin/time -l` is blocked under
the sandbox, so all RSS figures come from a trivial Python parent doing
`os.fork` + `execv` + `os.wait4` (macOS `ru_maxrss` is in bytes); the parent holds
nothing, and the fact that my 3 M `1G/c4` figure landed 0.4 MiB from the plan's is the
harness's own calibration.

---

## 1. Logic review

### 1.1 The CI guard's fixture size is wrong, and it inverts §5's conclusion — **Critical**

§5 asserts: *"3 M pairs reaches steady state and lands within ~3 points of the 10 M-pair
calibration matrix."* Both halves fail.

**Against the plan's own documents.** Phase 2 measured, on 10 M pairs at the same 51 bp
(`phase2/CALIBRATION.md` confirms: 10,000,000 pairs, uniform 51 bp):

| cell | Phase 2, 10 M | Phase 3 §5, 3 M | gap |
|---|---|---|---|
| `1G/c4` | 0.949 | 0.835 | 11.4 pts |
| `2G/c4` | 0.976 | 0.855 | 12.1 pts |
| `1G/c4 --fastqc` | 0.813 | 0.724 | 8.9 pts |

**By measurement.** Same generator, same cell, only the pair count changed:

```
3 M pairs, 1G/c4:  776.6, 731.6 MiB                        → max ratio 0.835
10 M pairs, 1G/c4: 789.3, 789.3, 814.5, 819.2, 820.1, 990.9 → max ratio 1.065
10 M pairs, 1G/c2: 851.8, 880.4, 885.8                      → max ratio 0.952
10 M pairs, 1G/c8: 791.3, 836.7                             → max ratio 0.900
10 M pairs, 2G/c4: 1715.7                                   → ratio 0.922
```

The 10 M ratios sit right on Phase 2's (0.952 / 0.900 against its 0.949 / 0.912), so the
driver is **fixture size, not distinctness**, and 3 M is 11 points (median) to 23 points
(max) short of steady state. The plan's `500 K → 0.32` data point established that too
small is bad; it did not establish that 3 M is enough, and 3 M is not.

**Why this inverts the design.** §5's third bullet argues the job *structurally* cannot
catch a lever-A revert: "Reverting lever A at `1G/c4` adds ~104 MiB: 776 → ~880 MiB,
still under the 930 MiB prediction." The arithmetic is right (4 workers × 1 slot × 1.44 ×
18 MiB = 103.7 MiB) but the baseline is the undersized one. At Phase 2's baseline for the
same cell — 0.949 × 930 = 883 MiB — the same +104 MiB lands at **987 MiB > 930**, and the
plain `peak ≤ prediction` assertion **catches it**. So "it is a promise guard, not a drift
tripwire" is not a structural property of the design; it is a consequence of choosing
3 M pairs. The accompanying argument ("buying sensitivity means buying fixture size,
without bound") is about scaling to *larger budgets*; at a fixed 1 GiB budget the data
needed to fill the pool is fixed, and 10 M supplies it.

**Failure scenario.** An implementer follows §5, lands the job at 3 M, and it goes green.
Six months later someone reverts lever A while refactoring `parallel.rs` (the two
`queue_depth` literals are precisely the shape that invites it). Peak at `1G/c4` moves
776 → 880 MiB, under the 930 MiB gate, job still green, and the plan has told the team in
writing that this job is "the only thing in the repo that catches a change making peak
RSS exceed the printed prediction." The pins in item 4 are supposed to cover this, and
per §1.3 below they do not either.

**Fix.** Size the fixture to reproduce the calibration matrix's ratio at `1G/c4` — 10 M
pairs does it, and gzipping the input (see §3.2) buys ~5 points more for a sixth of the
disk. Then re-read §1.2 before choosing the assertion, because at that size the assertion
as specified is red today.

### 1.2 `peak ≤ prediction` is not currently satisfied at `--memory 1G` — **Critical**

The plan's Goal opens: *"Phases 1 and 2 made the printed prediction a true upper bound on
measured peak RSS."* At 10 M pairs, `1G/c4`, one rep in six measured **990.9 MiB against a
printed 930 MiB** — ratio 1.065. Still inside the 1024 MiB *budget* (so the CHANGELOG's
user-facing claim survives), but above the *prediction*, which is the exact quantity
Step 6.5 proposes to gate on.

This was foreseeable from the parent plan without measuring. `PLAN.md:481-483` records:

> The margin rule as written in this plan ("margin must exceed the largest observed
> spread") is **not met at the default budget**: the margin is proportional (1/11) while
> the spread is absolute (96.2 MiB), so at `--memory 1G` it is 93.09 MiB and only clears
> 96 MiB above ~1.03 GiB.

Phase 3 then proposes to make `peak ≤ prediction` at `--memory 1G` a merge gate. My
observed spread at that cell is **201.6 MiB** (789.3 → 990.9) — 2.1× the 96 MiB the margin
was sized against, on a 100%-distinct fixture the calibration never used. The fitted line
sharpens it: `CALIBRATION.md`'s macOS c4 slope is `176.0 + 38.425·B`, so at B = 18.13 MiB
the model over-predicts a *typical* run by only ~58 MiB. 58 MiB of model slack plus 93 MiB
of margin against a 202 MiB spread is not a bound.

**Failure scenario.** The job lands with a fixture big enough to be sensitive. It goes
red on `dev` for no reason anyone can reproduce on demand, roughly one run in six. Someone
raises the threshold or drops the assertion, and the guard is gone — the standard fate of
a flaky gate. Alternatively it lands at 3 M, never fires, and §1.1 applies.

**Fix.** Split the two roles the single assertion is carrying:

1. **Gate on `peak ≤ budget`.** That is the promise the CHANGELOG makes, it is what a
   cluster user cares about, and it held in every one of my ~20 runs including the 990.9
   outlier. §5 calls it "redundant"; it is redundant only if the stronger claim is true,
   and today it is not.
2. **Log `peak/prediction` and `peak/budget` every run** (the plan already does this) but
   do not gate on the ratio.
3. **Run the 990.9 MiB outlier down before Phase 3 claims the prediction is an upper
   bound** — 6 reps at `1G/c4` on a 10 M-pair 100%-distinct fixture reproduces it. It may
   be a `--fastqc`-free tail-flush transient, it may be real allocator retention; either
   way the Goal paragraph and the §5 framing need to say what is actually true.

This also settles open question 2: a `≤ 0.92 × prediction` ceiling is unworkable. The
ratio's own variation at a sensitive fixture spans 0.849–1.065 at that cell, and the only
fixture where 0.92 is quiet is one too weak to test anything.

### 1.3 Step 4's swept test cannot fail under the perturbations Validation #5 names — **Critical**

§4 says the retired test "compares the formula against its own inverse, so it passes for
any coefficients", and that the sweep is "strictly more coverage … and unlike them it can
fail". The sweep has the same defect, for the same reason. `resolve_layout`,
`predicted_peak_bytes()` and `clumpify_min_memory_bytes` all read the *same*
`dyn_denominator` (`clump.rs:169`) and `static_bytes_for` (`clump.rs:202`). Perturb a
coefficient and all three move together. Taking the four assertions in Step 4.2 against
the three perturbations in Validation #5:

| assertion | σ 23→20 | `k` 64→32 | static 224→512 MiB |
|---|---|---|---|
| `predicted_peak_bytes() ≤ budget` on `Ok` | passes | passes | passes |
| `Err` below floor, `Ok` at/above | passes | passes | passes |
| `bin_byte_budget ≥ MIN_BIN_BYTES` on `Ok` | passes | passes | passes |
| no overflow | passes | passes | passes |

Worked: at `1G/c4` with σ = 20, denom = 576, so `B` = 20,590,220 and predicted =
234,881,024 + 20,590,220 × 576/16 = 976,128,944 B = 930.9 MiB ≤ 1 GiB. Passes. The
identity `predicted ≈ budget × 10/11` is coefficient-independent by construction, and the
floor is the algebraic inverse of the same expression, so it tracks too.

What *does* catch all three is the test the plan promises not to touch —
`resolve_layout_pins_the_default_cell` (`clump.rs:546`). I checked each: σ→20 gives
20,590,220; `k`→32 gives 23,911,222; static→512 MiB gives 11,263,026; the pin is
19,006,356. All three fail it. So Validation #5's "at least one assertion fails for each
perturbation" is technically satisfiable only by the *pre-existing* test, while Assumption
1 simultaneously requires that test to pass unmodified. The two rows contradict each
other, and the row's parenthetical — "the property the retired test could not see" — is
false.

**Failure scenario.** An implementer writes the sweep, runs Validation #5, sees
`pins_the_default_cell` fail and the new sweep pass, and reasonably concludes the
perturbation "worked". The plan's claim to have replaced a tautology with a real guard
goes into the PR description and the phase closes with the coverage gap intact.

**Fix.** Keep the sweep — it has real, narrow value (see §3.1) — but stop selling it as
the coefficient guard, and add the guard it is missing: extend
`pins_the_default_cell` into a small table of hand-computed
`(budget, cores, fastqc, workers) → exact bin_byte_budget` pins. Four to six cells is
enough, and it must include one `--fastqc` cell and one `workers: 1` cell so the two
Phase-2 post-review behaviours are pinned too. Values are one Python line each; I
computed the whole family above. Then rewrite Validation #5 to name which test fails.

### 1.4 The queue-depth pin does not observe the call sites — **Important**

`work_queue_depth(true) == 1` fails only if someone edits the helper. The regression the
handoff fears is a revert at the *call site* — and the literal `2` is what stands at
`parallel.rs:117` and `:943` today, so re-inlining it is the path of least resistance for
anyone touching that code. That leaves the helper untouched and the test green.

The plan calls this "the guard for the exact regression the handoff fears" (§4) and
Validation #4 tests it by editing the helper, which is the one mutation it does catch.
Extracting the constant is still worth doing — it removes the duplication and puts the
`k`-coupling comment next to `dyn_denominator`, which is right per the repo's
comment-placement rule — but the plan should say what it covers. If site coverage is
actually wanted, a `#[test]` asserting `!include_str!("parallel.rs").contains("sync_channel::<PairedWork>(2)")`
does bite; it is ugly, and I would rather see the honest limitation written down.

### 1.5 §6's 10 kb figure is input-limited and understates the defect ~2.8× — **Important**

The plan applies the input-limited caveat to the 50 kb cell only ("that fixture held only
4000 reads, i.e. fewer than one full 4096-record batch, so it is an underestimate"). The
same reasoning applies to the 10 kb cell and the plan does not say so. At 10 kb a
4096-record batch is ~82 MB, so a 382 MB fixture is **4.7 batches** against a queue
capacity of ~20 slots (1 reader batch + 4 workers × depth 2 + 4 in hand + `cores × 2`
result slots). It cannot reach the residency the plan's own mechanism predicts.

Measured, `--cores 4`, gzip output, 10 kb reads, my generator:

| records | input | peak RSS |
|---:|---:|---:|
| 19,100 | 382 MB | 768.5 MiB (the plan's cell) |
| 38,200 | 764 MB | 1366.0 MiB |
| 76,400 | 1.5 GB | 2121.2 MiB |
| 152,800 | 3.0 GB | **2282.6 MiB** ← plateau |

The plateau at ~2.3 GiB *is* "roughly 20 batch-equivalents" — the plan's mechanism is
right and its measurement never got there. Filing the issue with 812 MiB and "8× out"
understates the defect: the steady state is **~23×** the documented ~100 MB. My 150 bp
(77.8 MiB) and 1 kb (248.8 MiB) rows reproduce and *are* genuine steady states — 190,000
× 1 kb is 46 batches, comfortably past capacity — so only the 10 kb and 50 kb rows need
redoing.

**Fix.** Re-measure the 10 kb row at ≥ 152,800 records before filing, and state the
plateau as the ceiling. This is a public issue about a user-facing memory claim; the
number in it should be the one an ONT user will hit, and 812 MiB is not.

### 1.6 The issue misses the second wrong claim in the same sentence — **Important**

`threading.md:63` reads, in full: *"Trim Galore v2 stays under 100 MB across all
reasonable core counts … **Peak RSS scales with read length, not input file size:** at
`--cores 8`, a 1M-read 50-bp fixture sits at ~28 MB and an 84M-read 150-bp fixture at
~92 MB — bounded by the per-worker batch + channel buffers, which don't grow with total
input length."*

My §1.5 table refutes the bolded claim head-on: read length held at 10 kb, peak grows
768 → 2283 MiB purely with input size, until the queue saturates. The plan lists the
"~100 MB" halves of lines 61 and 63 and not this one. An issue that corrects the ceiling
and leaves "peak RSS does not grow with input file size" standing leaves the page wrong in
a way that matters more for scheduling, because it is the sentence that tells a user they
can ignore input size.

### 1.7 Two wrong rows in the `--clump_only` regression table, and prose that contradicts it — **Important**

The table says "Both floors computed from their own formulas, `workers: 1`". Recomputing
from `clump.rs` (workers = 1, `fastqc_threads` = cores, `div_ceil` to MiB as the binary
advises):

| cores | plan's new floor | recomputed | note |
|---:|---:|---:|---|
| 2 | 277 | **329** | 277 is the *no-FastQC* floor — wrong series |
| 8 | 513 | 513 | ✓ |
| 10 | 578 | **579** | truncated where the binary uses `div_ceil` |
| 11 | 611 | 611 | ✓ |
| 16 | 775 | 775 | ✓ (confirmed live: 774M aborts advising 775) |
| 32 | 876 | 876 | ✓ |

I confirmed the c2 slip directly: `--clump_only --memory 200M` (no `--fastqc`) exits
non-zero advising "try ≥ **277** MiB". So row c2 mixes the two series in a table whose
whole point is the FastQC comparison. Every old-floor value checks out, and both
conclusions survive — first regression at c11 ✓ (I swept 1–256), plain floor dropped at
every core count ✓ — so this is accuracy, not direction.

The prose does not survive: "inside a ≤ 150 MiB window" against a widest window of
**154.4 MiB at c16**, which the plan's own row shows as [620, 775). Small, but it is a
bound stated as a fact, in the paragraph that justifies D1 being cheap.

### 1.8 Handoff §2's "Unverified axes" row is dropped without disposition — **Important**

The handoff lists as pending: single-end, multi-pair (the promise is per-process; every
cell measured so far is single-pair), **short reads (σ ≈ 1.55 measured at 25 bp against
the shipped 1.4375)**, and low-duplication/WGBS data. None appears in Phase 3 — not in
scope, not in Remaining risks, not deferred to a named place.

The short-read one is not a nicety. σ 8% above the constant means the model *under*-charges
the bin pool on 25–35 bp libraries, which is the direction that breaks the promise, and
short reads are not exotic in this tool's user base. Phase 3 need not measure it — but a
phase whose stated goal is "defend the promise" should not let the one known
under-charging axis fall out of the paper trail. One line in Remaining risks, or fold it
into the Step 7 issue, or file a second issue.

### 1.9 A data-integrity item from Phase 1's "logged, not fixed" list is unrouted — **Important**

Phase 1's post-review recorded: *"Reviewer B found a pre-existing bug: the EOF bin flush
swallows a failed send."* Still live, and it is not in Phase 3 or anywhere else. The final
EOF batch is sent as `let _ = txs[idx].send(...)` at `src/parallel.rs:714` (paired,
no-passthrough), `:754` (passthrough) and `:1202` (single-end), while the *mid-stream*
sends at `:794` and `:825` correctly check `.is_err()`. If a worker has died, the last
partial batch is dropped and the process can still exit 0 — silent output truncation.

Phase 3 files exactly one issue. This one clears the same bar the threading issue does
(user-visible data integrity, arguably sharper) and costs one more `gh issue create`.
Either add it to Step 7 or state why it stays unfiled.

### 1.10 Consequences of D1 and D2 that were not visible when they were decided

Not re-litigating either; flagging what carrying them out exposes.

- **D2 raises the bar on the prediction being right.** Today the banner prints one number
  a user cannot interpret. After D2 it prints `930 of 1024 MiB budget` side by side, which
  is the point — but it also means a user who measures 990 MiB can see, unaided, that the
  tool overshot the figure it printed. §1.2's outlier is currently invisible; D2 makes it
  legible. That is an argument for running the outlier down *before* shipping D2, not
  against D2.
- **D2's `--memory` echo is exact but only by luck of the suffixes.** `1500M` → 1500 ✓
  (I ran it). The plan is right that sub-MiB budgets are unreachable and needs no guard.
  But `--memory 1500K` parses to 1,536,000 B → 1 MiB, and `--memory 1073741825` (1 GiB +
  1 B) prints 1024 — so the echoed figure is a truncation of the user's string, not the
  string. Harmless here because both fall below the floor and never reach the banner;
  worth one clause in the plan since the plan already reasons about this edge and stops
  one step short.
- **D1's remedy is unavailable at `--cores 1`.** The plan's docs paragraph says "`--cores`
  only sets bin count here, so lowering it is free and lowers the floor". At `--cores 1`
  the floor is 277 MiB and there is nowhere to go, yet the bail still ends "or decrease
  --cores" (I reproduced it). The paragraph should say that raising `--memory` is the only
  remedy at `--cores 1`.

### 1.11 Ordering and structure

The one ordering hazard the plan names (Step 6 must follow Step 1, because D2 is what puts
the budget on the banner line for a single-source parse) is correct and handled. Two
smaller ones it does not name:

- **`clump-only.md` line 75 is a bullet inside the "Composes with:" list** under
  `## Compatibility` (I read the file). Inserting a subsection "after the `--memory`
  bullet" splits that list in half. Put it after the list, or give it its own
  `## Memory` section.
- **The `needs:` key.** Step 6.2 says to reuse `validation-ubam`'s step block; that job
  carries `needs: [rust-tests, lint]`. Copying it contradicts Integration's "parallel with
  the rest, not on the critical path". Omit `needs` deliberately, or drop the claim.

---

## 2. Assumptions

Taking the plan's seven in turn.

1. **"The calibrated model is correct as shipped."** Fixed, and the right stance — but
   §1.2 shows the *promise* built on it (prediction as upper bound) is not currently
   satisfied at the default budget on a 100%-distinct fixture. The assumption should be
   narrowed to "no coefficient moves in Phase 3", which is what Phase 3 actually needs,
   rather than a correctness claim the phase then gates CI on.
2. **"`predicted_peak_bytes()` remains the exact inverse of `resolve_layout`."** Verified
   symbolically: predicted ∈ (budget×10/11 − denom/16, budget×10/11], so ≤ 312 B below
   `budget × 10/11` at cores ≤ 32. Sound. Note this exactness is *also* why Step 4's sweep
   cannot fail (§1.3) — the same property that makes the CI parse single-source makes the
   unit test vacuous. Worth stating in one place rather than being load-bearing in two
   sections with opposite signs.
3. **"`--clumpify` requires `--cores ≥ 2`."** Verified: the `--cores 1` floor is reachable
   only via `--clump_only`, and the plan's page assignment is right.
4. **"ubuntu-latest has ≥ 4 usable cores and ≥ 4 GiB RAM."** Holds on the current
   4 vCPU / 16 GiB public runner, and `2G/c4` needs ~1.7 GiB by my measurement. The
   under-stated constraint is **disk**: the plan quotes the fixture as 886 MB in §5 and
   ~1.2 GB in Efficiency (mine came to 978 MB for 3 M pairs with 56-char headers), and
   ubuntu-latest ships ~14 GB free against a restored cargo cache plus `target/`. §3.2's
   gzip recommendation removes the question entirely.
5. **"The CI cells and the 3 M fixture size are cost/coverage dials, not correctness."**
   This is the assumption §1.1 refutes. At 3 M the job cannot see a lever-A revert; at
   10 M it can. That is not a dial, it is the difference between the job doing its
   advertised work and not.
6. **"The measured ratios are macOS; Linux will differ."** Correct and well handled by
   Validation #11 — but the row's advice is one-directional ("if ≪ 0.7, raise the pair
   count") with nothing for the case §1.2 predicts, an intermittent red. Add the other
   branch.
7. **"Not assumed: that the CI guard covers lever A."** Honestly stated — and per §1.1
   the underlying claim ("cannot", structurally) is wrong in the plan's favour: with the
   right fixture it does cover lever A.

**One unstated assumption worth surfacing:** that plain-text fixture input is
representative. It is not, and it is the conservative direction (§3.2) — every real
`--clumpify` run reads `.gz`, which adds decompressor threads and buffers that the static
reservation is charged for.

---

## 3. Efficiency

The plan's efficiency reasoning is sound where it matters — Steps 1–4 and 7–9 are
text-or-compile-time, `work_queue_depth` is a branch executed twice per run, and the
sweep's ~1400 integer cells are nothing. Two refinements.

### 3.1 The sweep's real value is narrower than claimed, and worth keeping for it

Step 4.2's overflow clause earns its place, just not for the reason given.
`predicted_peak_bytes()` at `clump.rs:162` uses an unchecked `*`, as do
`clumpify_min_memory_bytes`'s two multiplies at `:216-217`. Today they are safe by a tight
argument: `resolve_layout`'s `saturating_mul(16)` bounds `B ≤ u64::MAX / denom`, so
`B * denom` cannot overflow — and I confirmed the bound holds up to `cores` 65536, the top
of the plan's proposed range. In a debug `cargo test` an overflow there is a panic, so a
sweep to `u64::MAX / 2` genuinely pins three unchecked multiplies against a future edit
that breaks the saturation chain. Say *that* is what the sweep buys, and add the
coefficient pins from §1.3 for the rest.

### 3.2 Gzip the CI fixture: more sensitivity, a sixth of the disk

Measured, 3 M pairs, `1G/c4`, same records, only the input encoding differs:

| input | reps | max ratio | on disk |
|---|---|---|---|
| plain `.fastq` | 776.6, 731.6 MiB | 0.835 | 978 MB |
| `.fastq.gz` (level 1) | 825.4, 818.1 MiB | **0.888** | **160 MB** |

Gz input buys ~5 points of sensitivity, halves the spread in my two reps, matches how the
tool is actually invoked, and takes the fixture from 978 MB to 160 MB — which also makes
Step 6.7's cleanup step and Assumption 4's disk question moot. `examples/mk_clump_fixture.rs`
should write `.fastq.gz` (or the workflow should gzip after generating). Combined with
10 M pairs this lands the job in the 0.95 region — sensitive enough to see lever A, which
is the whole argument of §1.1.

Note the plan's stated rationale for 100% distinct sequences — *"the conservative choice …
the per-worker term holds compressed output, which random sequence maximises"* — is not
supported by the data. My 10 M/100%-distinct ratios (c2 0.952, c8 0.900) sit essentially
on Phase 2's 10 M/52%-distinct ones (0.949, 0.912). Distinctness is roughly neutral here;
size is the axis. The choice is fine, the justification should go or be corrected, and the
warning that this fixture must never be reused for a compression claim should stay.

---

## 4. Validation sufficiency

Where an implementer could follow this plan exactly and still ship something wrong.

| # | Gap | What slips through |
|---|---|---|
| a | **Validation #5's expected outcome is wrong** (§1.3) | The phase closes believing the tautology was replaced. This is the single largest hole, because it is the validation whose job is to prove a guard bites. |
| b | **Nothing validates the banner's *numbers*, only its substrings** | Validation #1 is a manual eyeball. Step 1.2's tripwire clause after D2 pins the *budget*'s unit, not the peak's — the peak has no unit of its own anymore. Print `predicted_peak_bytes()` without dividing by 1024², or echo `memory_bytes` raw, and both clauses still pass. One line closes it: extract both integers with `predicted peak ≈ (\d+) of (\d+) MiB budget`, assert `peak < budget` and `budget == 1024` for the `--memory 1G` cell. |
| c | **Validation #10 proves the CI script can fail, not that it can fail *usefully*** | Setting `STATIC_TRIM_BYTES` to 32 MiB is a ~200 MiB perturbation. A guard with ~17 points of slack (§1.1) passes #10 and still cannot see a 104 MiB regression. The right form of #10 is: *revert lever A and confirm the job goes red* — the actual regression it is sold as protecting against. If that fails, the fixture is too small, and you learn it before merging rather than in a year. |
| d | **Validation #11 has no branch for an intermittent red** (§1.2) | "if ≪ 0.7, raise the pair count" walks the job straight into flakiness with no stated response. Add: if any cell exceeds 1.0 on unmodified code, stop and treat it as a finding, not a threshold to tune. |
| e | **Nothing checks `/usr/bin/time -v` exists on the runner** | `/usr/bin/time` comes from the `time` package and the bash builtin does not accept `-v`. This fails loudly on the first run rather than silently, so it is cheap — but it is a wasted CI round trip that `sudo apt-get install -y time` (or a `command -v` check) avoids. |
| f | **Nothing validates the issue's headline number before it is public** (§1.5) | Step 7 files 812 MiB. The steady state is ~2283 MiB. There is no validation row covering Step 7's content at all. Add one: re-measure the 10 kb cell at ≥ 8× the record count and confirm the plateau. |
| g | **Validation #8 ("docs figures match the binary") does not cover the two figures the plan reasons about but does not print** | The floor *range* endpoints (281/590/334/1012) are not printed by any invocation — they come from `clumpify_min_memory_bytes`, and only the single advised floor for the current cell surfaces. #8 as written ("run the cell and diff the banner/warning against the doc") cannot check them. I verified them by reimplementing the formula; say that is the method, or a `--cores 32 --fastqc --memory 1011M` invocation pins the top end from the warning text. |

The parts of Validation that *are* sufficient: #2 (coefficient pin unmodified), #3
(tripwire reverts), #6/#7 (`--clump_only` behaviour either side of the floor), #12 (docs
build). And the framing note under the table — "Validation 10 is not optional … prove this
one can [fail]" — is exactly the right instinct; it just needs the sharper mutation in (c).

---

## 5. Alternatives worth considering

1. **Gate on `peak ≤ budget`, log `peak/prediction`.** (§1.2) The budget is the promise the
   CHANGELOG makes and the number a cluster user puts in an sbatch script; it held in every
   run I did. The prediction is an internal artefact that happens to be printed. Gating the
   weaker, user-meaningful claim and logging the stronger one gives a job that is both
   non-flaky and legible.
2. **Make the drift guard a pin table rather than a sweep.** (§1.3) Four to six
   hand-computed `bin_byte_budget` values catch every coefficient edit, run in
   microseconds, and read as documentation of what the model resolves to. The sweep stays
   for the overflow bound.
3. **Skip the 3 M/10 M question by choosing a cheaper stressful cell.** Sensitivity scales
   with the baseline ratio, and the ratio is highest where the model is tightest —
   `1G/c2` measured 0.952 for me, versus 0.879 median at c4, and `CALIBRATION.md` records
   the model over-predicting by 22% at c16 on Linux. A `1G/c2` cell is the most sensitive
   per byte of fixture. Worth considering as the gating cell with c4 and `--fastqc` as
   logged-only.
4. **Do not file the threading issue from this plan's numbers.** (§1.5, §1.6) Re-measure
   first. An issue is a public artefact; 812 MiB when the plateau is 2283 MiB, plus a
   missed refutation of the "not input file size" claim in the same sentence, is a
   correction someone will have to make in public later.
5. **Consider folding D1's docs paragraph into `clumpy.md` as well.** The `--clumpify`
   page documents warn-and-fall-back (line 110) and, after Step 5, the corrected floor
   range. A reader who lands there and then uses `--clump_only` has no signal that the same
   budget behaves differently. One cross-reference sentence, not a section.

---

## 6. The four open questions — my recommendations

**Q1 — stale savings figures in `clumpy.md`.** Agree with leaving the numbers; disagree
with leaving them unmarked. After Step 5 the page will be self-contradictory: line 104
says "increasing memory doesn't seem to make very much difference", while the corrected
worked value two lines up shows the `4G/c8` bin moving 66 → 44 MiB — a 33% reduction in
exactly the quantity that claim is about. Add one sentence marking the savings table and
line 104 as measured against pre-refit sizing, with a pointer to #439. That is honest,
costs nothing, and does not fabricate a replacement. Re-measure the RRBS memory-scaling
row before the next release, as the plan says — it is the only row advertising a
memory-dependent benefit.

**Q2 — a ratio ceiling in CI.** Do not add it, and do not defer it either — close it.
My data kills it outright: at a fixture stressful enough to be worth gating, the ratio at
`1G/c4` spans 0.849–1.065 across 6 reps. A 0.92 ceiling would be red most runs, and the
only fixture where 0.92 is quiet is one too weak to detect anything (§1.1). The
sensitivity the ceiling is meant to buy belongs in the deterministic pins (§1.3), which
cost nothing and cannot flake. Measuring the cells on oxy would not change this — the
spread is the problem, not the offset.

**Q3 — should `--clump_only` print a predicted peak?** Yes, separate PR — agree, and note
it is cheaper than the plan implies. `predicted_peak_bytes()` already exists and
`--clump_only`'s layout already carries `workers: 1` and its own `static_bytes`, so the
printed figure would be correct with no model work at all: it is one added format
argument at each of the four `eprintln!` sites. Bundle it with Q4, since both edit the same
four lines.

**Q4 — the "budget" terminology collision.** Agree: drop "budget" from
`clump_only.rs:289`, `:784` and `:1031`, leaving all four banners as `N bins × M MiB`
(`:436` already reads that way — I confirmed all four). No test asserts on those
substrings. One addition: do it **in the same PR as D2**, not later. If D2 ships alone,
there is a released binary in which "budget" means `--memory` on one banner and per-bin
bytes on three others, which is the exact confusion the change exists to prevent.

---

## 7. Action items

### Critical

1. **Re-size the CI fixture, or stop claiming the job guards drift.** (§1.1) 3 M pairs
   measures 0.835 at `1G/c4` against 0.949 for a properly-sized fixture; the gate has ~17
   points of unearned slack and would pass a lever-A revert. 10 M pairs plus gzipped input
   (§3.2) lands in the 0.95 region. Then re-derive §5's third bullet — with that fixture
   the plain `peak ≤ prediction` assertion *does* catch lever A, which is a better story
   than the one the plan tells.
2. **Change the gating assertion to `peak ≤ budget`, and run the 990.9 MiB outlier down
   before Phase 3 asserts the prediction is an upper bound.** (§1.2) One rep in six at
   10 M/`1G/c4` measured 990.9 MiB against a printed 930 MiB. `PLAN.md:481-483` already
   recorded that the 1/11 margin does not cover the observed spread at `--memory 1G`; this
   is that, measured. Gate the budget, log the prediction ratio, and fix the Goal paragraph
   so it states what is true.
3. **Fix Step 4 and Validation #5.** (§1.3) None of the four assertions can fail under
   σ 23→20, `k` 64→32 or static 224→512 MiB — the sizing, the prediction and the floor all
   read the same two helpers, so they move together. Keep the sweep for the overflow bound
   (`clump.rs:162`, `:216-217` all use unchecked `*`), add a 4–6 cell table of
   hand-computed `bin_byte_budget` pins including one `--fastqc` and one `workers: 1` cell,
   and rewrite Validation #5 to name the test that fails.

### Important

4. **Re-measure the 10 kb row before filing the issue.** (§1.5) 382 MB is 4.7 batches
   against ~20 slots of capacity; the plateau at 152,800 records is 2282.6 MiB, so the
   issue's "8× out" is really ~23×. Apply the plan's own 50 kb caveat to the 10 kb cell.
5. **Add `threading.md:63`'s "peak RSS scales with read length, not input file size" to the
   issue.** (§1.6) The §1.5 table refutes it at fixed read length. Correcting the ceiling
   and leaving the scaling claim is a half-fix of one sentence.
6. **Correct two rows and one prose bound in the `--clump_only` regression table.** (§1.7)
   c2 = 329 MiB not 277 (277 is the no-FastQC floor — I reproduced it live); c10 = 579 not
   578 (the binary uses `div_ceil`); widest window is 154.4 MiB at c16, so "≤ 150 MiB"
   contradicts the plan's own row.
7. **Say what the queue-depth pin covers.** (§1.4) It catches an edit to the helper, not a
   re-inlined literal at `parallel.rs:117`/`:943` — which is the shape the regression will
   actually take. Keep the extraction; drop "the guard for the exact regression the handoff
   fears".
8. **Give handoff §2's "Unverified axes" a disposition.** (§1.8) Short reads at σ ≈ 1.55
   against the shipped 1.4375 under-charge the pool — the promise-breaking direction. One
   line in Remaining risks or a second issue; not silence.
9. **Route the EOF-flush swallowed send.** (§1.9) `src/parallel.rs:714`, `:754`, `:1202`
   use `let _ = …send(…)` where `:794` and `:825` check `.is_err()`; a dead worker means a
   silently dropped final batch and exit 0. Add it to Step 7 or say why not.
10. **Two CI mechanics.** (§1.11, §4e) Omit `needs:` when copying `validation-ubam`'s block
    or drop the "not on the critical path" claim; and install/verify `/usr/bin/time`. Also
    reconcile the fixture size quoted twice (§5 "886 MB" vs Efficiency "~1.2 GB").
11. **Sharpen Validation #10 to "revert lever A and confirm the job goes red."** (§4c) A
    192 MiB static perturbation proves the plumbing works; the 104 MiB regression is the
    one that matters.

### Optional

12. **Assert the banner's numbers, not just its substrings.** (§4b) `predicted peak ≈ (\d+)
    of (\d+) MiB budget`, then `peak < budget` and `budget == 1024` at `--memory 1G`.
13. **Move the `clump-only.md` below-floor subsection out of the "Composes with:" list.**
    (§1.11) Line 75 is a bullet; a heading there splits the list.
14. **Say that raising `--memory` is the only remedy at `--cores 1`.** (§1.10) The bail's
    "or decrease --cores" is degenerate there.
15. **Name the core-count domain of the docs floor range.** 281–590 MiB holds for
    `--cores` 2–32; at `--cores 64` the floor is 933 MiB. "Any sensible core count" is
    carrying real weight.
16. **Define `10/11` and `workers` on the docs page.** (§1.3 of the plan's table) Step 5
    introduces both terms into `clumpy.md` and the page defines neither. One clause each.
17. **Extend #439's existing CHANGELOG bullet rather than adding a second.** It is still in
    `### Unreleased`, so two bullets for one issue in one section reads oddly. One clause
    naming the banner.
18. **Drop the "100% distinct is conservative" justification.** (§3.2) My 10 M/100% ratios
    sit on Phase 2's 10 M/52% ones; distinctness is neutral here, size is the axis. Keep
    the "never reuse this for a compression claim" warning.
19. **`ClumpLayout::new` (`clump.rs:150`) is still public with no production caller and
    validates neither invariant `resolve_layout` enforces** — carried from Phase 1's
    "logged, not fixed". Out of Phase 3's scope, but it is the API through which a future
    caller builds a layout whose prediction means nothing.

---

## 8. What is sound, briefly

Worth saying plainly, because most of this plan is right and the findings above are
concentrated in two places.

- **Every correction to the handoff is correct**, and I verified each against the binary
  rather than the arithmetic: 930 not 931; the floor ranges being `--cores 1` values on a
  page that requires `--cores ≥ 2`; the `800M` cell resolving because of Phase 2's
  `workers: 1` split (I reproduced the stale 840.4 MiB figure by putting `workers: 16`
  back). Re-deriving instead of citing was the right call and it caught real staleness.
- **All 11 `clumpy.md` line references are exact**, as are `threading.md:61`/`:63`, the
  four `resolve_layout` sites, `main.rs:207-213` and the tripwire at `:599`. No stale
  `file:line` refs — notable given the handoff's own warning that two earlier reviews
  caught them.
- **§5's 3 M measurements reproduce to within 0.4 MiB**, which is what let me trust the
  same harness at 10 M. The *measurements* are good; it is the inference from them that
  needs work.
- **§6's mechanism is right in every respect I could test**: the `--cores 1` bypass
  (10.2 MiB against the plan's 10.6), flatness past c4 (783 MiB at c8), the batch-count ×
  record-length story, and "roughly 20 batch-equivalents" — which turns out to be exactly
  where the peak plateaus once the fixture is big enough to get there.
- **D1's rationale is well-argued and its risk assessment is right**: a `--clump_only`
  fallback would emit a byte-copy of the input at full I/O cost and exit 0, and the error
  already names the budget to retry with. The bail genuinely does not name the mode today
  (I reproduced: `--memory budget too small for --cores 16: …`), so Step 2's context wrap
  fixes something real at zero behaviour risk.
- **The self-review is unusually honest** — naming the guard's blind spot, the macOS-only
  ratios, and the stale savings figures as choices rather than oversights. Findings 1–3
  are places where that honesty stopped one measurement short, not places where it was
  absent.

**Scope note:** read-only on source, as instructed. `git status` on `src`, `tests`,
`Docs`, `.github` and `Cargo.toml` is clean; the tree is still at `fdade7e`; no branch
switch, stash or worktree. All measurement ran in the session scratch directory. My only
write is this file. (A cleanup `rm` of the scratch fixtures was declined by the permission
gate; the scratch directory is session-scoped, so they go with it.)
