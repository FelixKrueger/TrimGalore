# Plan review — `phase3/PLAN.md` **v3** (defend and correct the `--clumpify` memory promise, #439)

**Reviewer A**, round 3. Target: `plans/08162026_clumpify-memory-accounting/phase3/PLAN.md` at v3.
Repo `/Users/fkrueger/Github/TrimGalore`, branch `dev`, HEAD **`fdade7e`** (verified equal to the plan's
baseline). Read-only: nothing in `src/`, `tests/`, `docs/`, `.github/` or `plans/` edited except this
file; no branch switched; no `chmod` run anywhere.

> **Provenance.** I ran **zero `trim_galore` invocations** — no fixture, reps, free disk or exit codes
> to report, because I generated no measurements. Findings are tagged:
> **[V-code]** source read at `fdade7e`, command quoted · **[V-arith]** Python reimplementation from
> the constants as read · **[NV]** not verified by me. Every RSS figure in the plan is **[NV]**; §0 I
> read as a retraction and did not re-measure, per instruction.

## Verdict

**Implement.** v3 is materially correct: I re-derived every arithmetic claim I could reach and found
**one stale number, and it is in the parent plan rather than here** (§5.1). The guard rework is sound,
§5.3's band is the right answer to the platform hole, and §5.4 is a real defect class that v2 missed.

Four findings, none of them a redesign: **one Critical** (4d's deletion is ordered before its
replacement, which the parent plan explicitly forbade in the same sentence that created the
assignment), and three Important — two coverage gaps in §4a where the guard boundary sits one layer
below the behaviour Phase 2 actually fixed, and an unreconciled factor of 1.44 in §5.3's delta.

---

## 1. §4a's six pinned cells — are they the right six?

**For `resolve_layout` itself: yes, and I would not change them.** The set covers every coefficient
with at least one cell, `workers ≠ bin_cores` separates the two `cores` roles, and the 16/17-thread
pair is the only way to catch a `FASTQC_CHARGED_THREAD_CAP` change — the plan is right to say the
table must not be trimmed, and right to put that on the table itself.

I also checked the `n_bins` mapping, which the perturbation matrix does not list **[V-arith]**:
`max(16, 4·bin_cores)` → changing the `16` breaks the `--cores 2` cell (n_bins 16→8), and changing the
`4` multiplier breaks the default cell (`max(16,16)=16` → `max(16,20)=20`). Both covered. Good.

**But the guard boundary is drawn one layer too low, and the two behaviours Phase 2 is credited with
fixing sit above it.** 4a pins `resolve_layout` *given* a `LayoutInputs`. Nothing pins the
*construction* of that struct — which is where both Phase 2 fixes live.

### 1a. `LayoutInputs` construction is completely unguarded  **[Important]**

**[V-code]** `resolved_threads` has five call sites — `src/main.rs:187` and
`src/clump_only.rs:281, 427, 770, 935` — and **zero tests**. `src/fastqc.rs`'s test module
(`mod tests` at `:126`) contains eight `test_apply_fastqc_args_*` tests and nothing for
`resolved_threads`. So:

- A revert to `fastqc_threads: Some(cores)` when `--fastqc_args -t 2` was given — the exact bug
  PROGRESS records Phase 2 as fixing — **passes all six 4a pins, passes 4b, and passes the CI job**.
- A revert to `workers: cores` in `clump_only.rs` — the other Phase 2 fix — likewise passes everything.

Worse, the CI job **structurally cannot see either**: its cells are `--cores 4` with default FastQC
threading, where `bin_cores == workers == fastqc_threads`, so all three roles are numerically equal and
no substitution among them is observable.

**Fix (cheap, pure, no fixture):** two unit tests — `resolved_threads(Some("-t 2"), 8) == 2`, and one
asserting `clump_only`'s inputs are built with `workers: 1` while `bin_cores` tracks `--cores`. Add a
sentence to §4a saying the table guards the *resolution* and these guard the *inputs*; otherwise a
reader of §4a reasonably concludes the six cells are the whole guard.

### 1b. Lowering `MIN_BIN_BYTES` is undetectable  **[Important]**

All three floor-related tests are **relative to the function under test** — the same algebraic-inverse
vacuity §4 spends a section proving about v1's sweep **[V-code]**:

```rust
// src/clump.rs:603-610
let floor = clumpify_min_memory_bytes(&inputs);
let res = resolve_layout(floor - 1, &inputs);
```

Both sides read `MIN_BIN_BYTES`, so they move together. Direction matters:

- **Raising** it 1 MiB → 4 MiB is caught, because four of §4a's six pins sit below 4 MiB
  (18.13 / 3.36 / 2.07 / **1.95** / 22.80 / 2.27 MiB) and would bail instead of resolving. ✔
- **Lowering** it 1 MiB → 256 KiB is caught by **nothing**. Every pin still resolves to the same value,
  both bail tests recompute the floor and still pass, and 4b's `Err`-below/`Ok`-above property is
  preserved by construction. Yet every documented floor moves — the 281 / 311 / 334 / 590 / 1012 MiB
  figures that D4 is publishing in the same phase.

**Fix, and it buys two things at once:** one absolute assertion, ideally
`clumpify_min_memory_bytes(&ins(2, false)) == 281 MiB` expressed in bytes. That pins the constant *and*
pins the headline docs figure PR 1 publishes, so the docs and the binary cannot drift apart silently.
A bare `assert_eq!(MIN_BIN_BYTES, 1024*1024)` would also work but guards less.

---

## 2. §4d — does it discharge the assignment?  **[Critical, ordering only]**

**In substance, yes.** The parent plan's assignment (`../PLAN.md:170`, `:332-336`) offers "retire it or
move the assertion somewhere it bites"; 4d retires it and folds the intent into 4a, and Validation 6a's
`grep -c … == 0` plus striking the parent note is a proper closure check. The diagnosis is also right,
and I reached it independently on the parent plan: the test is a restatement of
`MARGIN_NUM ≤ MARGIN_DEN`, and #446 made it *more* vacuous — the chain is now
`predicted ≤ budget × 10/11 ≤ budget`, clearing with 9% to spare **[V-code]** (`clump.rs:577`).

**But as scheduled it does the one thing the parent plan explicitly forbade.** The parent states the
condition and its rationale in the same breath (`../PLAN.md:332-336`):

> "…**assigned to Phase 3**, since the plan pairs its retirement with 'Phase 3's CI check is the real
> replacement' — **deleting a passing test before its replacement exists would reduce coverage, not
> improve it**."

v3 schedules the deletion in **PR 2 step 8** (`:518`) and the replacement in **PR 3**, with
**"Ordering: PR 3 after PR 2"** (`:568`). So between those two merges the repo has the deletion and not
the replacement — precisely the state the parent ruled out. It is harmless *in substance*, because what
is deleted is vacuous and therefore no coverage is really lost; but the plan is currently asserting
compliance with a precondition it inverts, and Validation 6a will report "discharged" during the window.

**Fix — pick one and say which:** (a) move the deletion to PR 3, which already touches nothing in
`clump.rs` but could; or (b) keep it in PR 2 and record the deviation explicitly — "the parent's
precondition protects against losing real coverage; this test has none, so the ordering is accepted".
(b) is defensible and cheaper. What is not fine is leaving it implicit.

---

## 3. §5.3 — the band is right; its delta is derived two ways and only one is stated  **[Important]**

**The design is correct and the platform argument is the sharpest thing in v3.** Gating on
`peak ≤ prediction` alone cannot see a lever-A revert at the Linux baseline, the band catches drift in
both directions, and the lower bound catching a fallback or an unsaturated fixture is a genuinely better
instrument than a tuned ceiling. I have no objection to D5.

**The +104.2 MiB figure needs its provenance stated.** The plan derives it as
`4 × 1.4375 × 18.13 MiB` — i.e. `cores × σ × B`. But the model's own per-worker unit is `1 × B`, not
`σ × B`: `dyn = B × (23·n_bins + 64·workers)/16 = B × (1.4375·n_bins + 4·workers)` **[V-arith]**, so
`k = 4` charges `4·B` per worker and one extra queued batch is `+1·B` per worker
= `4 × 18.13 = 72.5 MiB`. Two defensible readings, differing by 44%:

| | multiplier | delta at `1G/c4` | Linux 790 + delta | vs 930 gate | band move |
|---|---|---:|---:|---|---:|
| physical bytes held (plan's) | σ·B | 104.2 MiB | 894 | green, misses by 35 | 11.2 pts |
| the model's own charge | 1·B | 72.5 MiB | 862 | green, misses by 67 | 7.8 pts |

**§5.3's conclusion is robust either way** — both stay under the 930 gate (so the absolute assertion
alone is insufficient, as claimed) and both exceed the `+0.05` upper bound (so the band fires). Say
that, rather than resting on one number: the argument is *stronger* stated as a range, because it holds
across the ambiguity. And have **Validation 12 record the measured delta**, so the band's headroom
becomes a known quantity instead of an assumed one — it is the one number in §5.3 that the job itself
can supply.

Aside, not an action: if a queued batch really holds `σ·B` while `k` charges `1·B` for it, that is a
small standing under-charge in `k`'s queued component. It is inside the 9% margin and not worth
reopening; noting it so it is not rediscovered as novel.

---

## 4. §5.4 — correct, and one implementation trap worth naming

The mechanism is right and I have nothing to add to the diagnosis: a stalled writer backs up
`sync_channel(cores * 2)` and the `pending` map, both holding written-but-unflushed buffers, so disk
pressure presents as a memory-model breach. Making a *parseable log from a failed run* the explicit
target of Validation 14 — distinct from Validation 13's unparseable log — is exactly the right
distinction, and it is the one that would have caught v2's §0.

The three clauses are sufficient as written (assert free space, truncate between cells, check every
exit code). One trap for the implementer: **`/usr/bin/time -v` exits with the child's status**, so a
single check on the wrapper covers the binary — but only if the pipeline does not swallow it. With
`set -euo pipefail` and `/usr/bin/time -v … 2> log`, the status propagates; piping into `grep` or `tee`
would replace it with the filter's status, and `grep` succeeding on a failed run's parseable banner is
the exact hole being closed. Worth one clause in PR 3 step 4, since the plan already mandates
`set -euo pipefail`.

---

## 5. What the phase inherits from the parent — systematic sweep

The lead asked for this; here is the complete set. **[V-code]** `grep -n -i "phase 3" ../PLAN.md`
returns four hits, and I checked each against v3:

| parent | assignment | v3 |
|---|---|---|
| `:102` | "CI check, docs arithmetic, banner unit label" | all three present — CI §5, docs §3, and the unit label is already shipped (`main.rs:208` prints `MiB`, not `MB`) ✔ |
| `:170` | retire `resolve_layout_predicted_peak_fits_budget` | **4d** ✔ (ordering caveat §2) |
| `:332-336` | same, with the precondition | **4d** ✔ (this is the sentence §2 quotes) |
| `:340-342` | `clumpy.md`'s "32 bins × 12 MB" is "doubly wrong: the figure … and the unit label" | §3's line-102 row corrects **both** — 12 → 9 *and* MB → MiB ✔ |

**No further dropped assignments.** 4d was the only one, and it is now in.

### 5.1 The one stale number, and it is in the parent plan

**[V-arith]** `../PLAN.md:328-330` describes `resolve_layout_pins_the_default_cell` as asserting
**19,884,107 bytes**. The shipped test asserts **19,006,356** (`src/clump.rs:551`), and recomputing from
the constants gives exactly that:

```
budgeted = 1073741824 × 10 / 11 = 976,128,930
usable   = 976,128,930 − 234,881,024 = 741,247,906
denom    = 23×16 + 64×4 = 624
B        = 16 × 741,247,906 / 624 = 19,006,356      (18.1259 MiB)
```

v3's §4a baseline pins lead with `19_006_356` — **correct**. So this is parent-plan doc rot, not a
phase defect, but it is the number a reader would quote when describing the guard, and PR 2 is already
striking a parent-plan note two paragraphs above it. **Fix it in the same pass.** [Optional]

---

## 6. §0 — read, not re-measured

The retraction is well made and I am not going to relitigate it. Two notes where v3 quotes my
round-2 report, both minor and both in my favour to correct rather than leave:

1. **The budget at that cell is 1025 MiB, not 1024.** The sweep harness computes
   `round(19 × 108/4) + 512 = 1025`, and `results_sweep.tsv` records `budget_mib = 1025`,
   `predicted_mib = 1025` **[V-arith]**. So the gap was 25.2 MiB, not 26. Reviewer B caught this
   against my report and I corrected it there; §0.3 has inherited the 1024.
2. **"Unsatisfiable *because of* the pre-refit coefficients" understates it slightly.** My claim was
   conditional on `B` being *pinned* at 18.96 MiB: with behaviour pinned, `actual ≈ 1050` is a physical
   fact, so no coefficient vector — margin included — can make `predicted` both ≥ 1050 and ≤ 1025, since
   a margin only *shrinks* `predicted` and tightens the bind. B worked this through independently and
   withdrew the objection. The resolution came from **`B` moving 4.40%** (18.96 → 18.126 MiB), not from
   the margin alone. That distinction is the actionable half: the parent Goal's "do not move the bin
   size at the default" constraint is what had to give, and #446 gave it. v3's framing reads as though
   the margin absorbed it.

Neither changes a decision. I flag them only because §0.3 is the paper trail for a retraction and
should be exact about what it retracts.

---

## 7. Where v3 is right about things I probed

Recorded briefly; all **[V-code]** or **[V-arith]** at `fdade7e`.

- The `k = 4` decomposition closes exactly: work 1 + 1 in hand + `result_per_core` 2 = 4
  (`clump.rs:246-248`; `parallel.rs:117`, `:130`, `:943`, `:959`). §Signature's `channel_depths`
  matches shipped behaviour, so 4c is a pure refactor and Integration's "no output bytes change" holds.
- Every docs figure in §3 reproduces under `div_ceil`: 281 / 311 / 334 / 590 / 1012 / 933 / 1356 MiB,
  the `clump_only` series 329 / 513 / 579 / 611 / 775 / 876, the plain 352 / 454, crossover exactly
  `--cores ≥ 11`, the 154.4 MiB window at `c16`, bins 9 MiB at `1G/c8` and 44 MiB at `4G/c8`, and the
  banner's 930 and 1363.
- §Context's correction of v2's own rounding error (the `--cores 1` floor is 277 / 303 by `div_ceil`,
  not the 276 / 302 truncations) is right — and catching your own prior table two rows from where you
  fixed the same class of error is the kind of thing that usually goes unnoticed.
- `ClumpLayout::new` is `pub` with exactly one caller, in test code (`parallel.rs:1615`), so
  `pub(crate)` is safe.
- The EOF-send issue's line numbers and the two additions I suggested (the `:1193`/`:1202` proximity,
  the #434 linkage) are carried correctly into §7.
- §5.1's inversion of v2's gzip rationale — that gz is the *more sensitive* fixture at the gating cell,
  which is what matters once §5.3 makes sensitivity binding — is a better argument than the one it
  replaces, and it is honest about that.

---

## Action items

### Critical
1. **Resolve 4d's ordering** (§2). Either move the deletion to PR 3, or record the deviation from
   `../PLAN.md:332-336`'s explicit "deleting a passing test before its replacement exists" condition.
   Do not leave it implicit while Validation 6a reports discharge.

### Important
2. **Add two unit tests for `LayoutInputs` construction** (§1a) — `resolved_threads` has five call
   sites and zero tests, so both Phase 2 behaviour fixes are unguarded at every level, and the CI cells
   at `--cores 4` cannot distinguish the three `cores` roles. Note in §4a that the table guards
   resolution, not inputs.
3. **Pin `MIN_BIN_BYTES` absolutely** (§1b) — lowering it is invisible to all six pins and both bail
   tests, and it moves every floor figure D4 publishes. `clumpify_min_memory_bytes(&ins(2,false)) ==
   281 MiB` guards the constant and the docs number together.
4. **State §5.3's delta as a range** (§3) — `72.5–104.2 MiB` depending on whether the queued batch is
   charged at `1·B` (the model's unit) or `σ·B` (physical); the conclusion holds across both, which is
   a stronger argument. Have Validation 12 record the measured delta.

### Optional
5. Add the exit-status clause to PR 3 step 4 (§4): don't pipe `/usr/bin/time` into `grep`/`tee`, or the
   child's status is replaced by the filter's — the precise hole §5.4 exists to close.
6. Fix `../PLAN.md:329`'s stale `19,884,107` → `19,006,356` (§5.1) in the same pass that strikes the
   `:332-336` note.
7. Correct §0.3's two inherited details (§6): budget 1025 not 1024, and that `B` moving 4.40% is what
   resolved the constraint.

---

## Operational note — `chmod 444` does not protect an artifact from this toolchain

Answering the lead's factual question, since it bears on snapshot policy. **I did not chmod anything.**
I tested the mechanism in my own scratchpad rather than speculate:

| action on a mode-444 file | result |
|---|---|
| shell in-place append (`>>`) | **blocked** — "permission denied" |
| shell replace-by-rename (`> tmp; mv -f`) | **succeeds**, and mode becomes `-rw-r--r--` — a visible tell |
| **`Edit` tool** | **succeeds** (18 → 40 bytes) and **mode stays `-r--r--r--`** |

That reproduces the parent report's state exactly — 50,656 bytes, still `-r--r--r--`, no chmod by me.
So the 20:59:37 write was my status-block `Edit` going straight through the lock at 20:57:58.

**The consequence for you is worse than "the lock is weak": a bypassed lock leaves no trace.** The mode
still reads 444 afterwards, so `stat` cannot distinguish a respected lock from a bypassed one — only a
size or content comparison against the snapshot can. Snapshots are your only protection, and they are
also your only *detector*. If you want a real lock, the directory needs to be read-only (the rename in
row 2 succeeds because the *directory* is writable), or the artifact needs to live somewhere tracked.

I have not touched any locked file since, and will ask before amending one.

---

*Reviewer A, round 3. Code claims verified at `fdade7e` with commands quoted; all floor, bin, banner and
pin arithmetic reimplemented from the constants. No `trim_galore` invocation was executed by me, so no
RSS figure here is mine; every measurement in the plan is **[NV]**.*
