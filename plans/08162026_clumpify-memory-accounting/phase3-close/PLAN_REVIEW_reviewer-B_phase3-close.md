# Plan review — Reviewer B — `phase3-close/PLAN.md`

**Target:** `plans/08162026_clumpify-memory-accounting/phase3-close/PLAN.md`
**Reviewed against:** `dev` = `242546a`, branch `ci/clumpify-memory-guard-439` tip `3f1a930`, PR #452
(`mergeable: false`, `mergeable_state: dirty` — confirmed via REST).
**Measurement performed:** none. No fixture created, no `trim_galore` invocation. Every number
below is either read from committed source, recomputed from figures already in the plan corpus, or
taken from the GitHub Actions API.

**Verdict: sound plan, approve after fixing two silent-inertness holes.** The rebase argument is
correct and I can strengthen it to a proof. The two-condition band claim is correct. A3's
peak/budget conversions are all correct. But the plan's own central deliverable — a drift band that
actually gates — has **two paths to shipping inert while appearing live**, and one of them is the
very validation the plan wrote to prevent that. Both fixes are a few lines.

---

## What I verified as correct — do not re-litigate these

| Claim | Status |
|---|---|
| `1b2579e` is an ancestor of the branch; `origin/feat/clumpify-guards-439` is **not** | code-verified |
| `1b2579e..branch` = exactly the four `ci(clump):` commits | code-verified |
| `44dc732` vs `1b2579e` differs in `CHANGELOG.md` only | code-verified |
| `git merge-tree` → `CONFLICT (content): CHANGELOG.md`, `Auto-merging src/clump.rs` | code-verified, reproduced verbatim |
| The band needs **both** a data row and `bootstrap=0`; deleting `# BOOTSTRAP` alone leaves it inert | code-verified, and I ran the extracted band logic |
| TSV has comment rows only; `# BOOTSTRAP` occurs exactly once (line 14); no comment row can match `$1==label` | code-verified |
| `--clumpify` requires `--cores >= 2` (`src/cli.rs:742-746`) | code-verified |
| Workflow triggers are `push: [master, dev]` + `pull_request: [master, dev]`, so B1's rep 4 from the merge push is real | code-verified |
| A3's four peak/budget conversions (0.76/0.80/0.89/0.94) | arithmetic-verified against the shipped model |
| The job's `3m59s` | verified from the Actions API for run `32065624366` |
| Assumption 3 (appending cannot perturb earlier cells) | code-verified — sequential `measure` calls, per-process `ru_maxrss`, no shared mutable state |
| The margin pin does **not** misfire at `--cores 2` | arithmetic-verified — see I5 |

Two of these I could only confirm by going further than the plan does, and the extra step matters:

**The rebase is provably clean, and the plan's stated reason does not prove it.** The plan argues
from "the only file differing between `1b2579e` and #451's tip is `CHANGELOG.md`". That is the wrong
comparison — you are rebasing onto `origin/dev`, which also carries #448. The proof is disjointness
of two sets:

```
git diff --name-only 1b2579e origin/dev
  CHANGELOG.md
  docs/src/content/docs/modes/clump-only.md
  docs/src/content/docs/performance/clumpy.md

git diff --name-only 1b2579e origin/ci/clumpify-memory-guard-439
  .github/clumpify-memory-baseline.tsv
  .github/workflows/ci.yml
  examples/mk_clump_fixture.rs
  src/clump.rs
```

The intersection is empty, so for every file the four replayed patches touch, `origin/dev`'s blob is
**identical** to `1b2579e`'s. Each patch applies against byte-identical context. Conflict-free is
guaranteed, not merely likely. Put this two-command check in the plan in place of the current
argument — it is shorter *and* it is a proof.

**The same check also closes a hole in B1 that the plan does not know it has.** See I2.

---

## Critical

### C1. Validation 5's perturbation fires the wrong clause — it cannot detect a broken upper bound

Marked: **code-verified** (I extracted the band logic from the job and ran it).

Validation 5 says: *"after B, perturb one row by +0.20 and confirm the cell goes red."* The band is

```bash
awk -v r="$ratio" -v w="$want" 'BEGIN{exit !(r > w + 0.05)}'   # upper — catches regressions
awk -v r="$ratio" -v w="$want" 'BEGIN{exit !(r < w - 0.10)}'   # lower — catches fixture decay
```

Perturbing the row **upward** by 0.20 makes `want = 0.976` against an observed `r ≈ 0.776`, so
`r < want − 0.10` holds and the **lower** clause fires. I ran it:

```
want=0.976  observed=0.776  -> LOWER (fixture clause)
want=0.899  observed=0.699  -> LOWER (fixture clause)     # 1G_c4_fastqc
want=0.576  observed=0.776  -> UPPER (regression clause)  # i.e. -0.20
```

The upper clause is the half the plan itself calls load-bearing ("the upper bound is the load-bearing
half", B2). **A typo'd, inverted, or deleted upper clause passes Validation 5 unchanged.** This is
the same shape as the harness failure already recorded in memory — a mutation whose direction agreed
with its assertion, so the test could only ever pass — and the same shape as risk (b)'s "asserted to
work, not proven to work". The plan is spending its most important validation on the clause it cares
least about.

**Fix (three lines of plan text).** Perturb each direction and assert the *message*, since the two
clauses emit distinguishable text:

| Perturbation | Expected `::error::` substring |
|---|---|
| row **− 0.20** | `is >0.05 above baseline` |
| row **+ 0.20** | `is >0.10 below baseline` |

Boundary behaviour, verified: `r = w + 0.050` passes, `r = w + 0.051` fails — the comparison is
strict `>`, so a perturbation of exactly the band width does not fire. Use ±0.20, not ±0.05.

### C2. Nothing asserts that a baseline row was actually found — a label typo silently re-inerts the band

Marked: **code-verified** (ran the lookup against a synthetic TSV; and against a missing file).

The plan correctly identifies a two-condition trap and then walks into a third form of it. After B,
`want` empty is indistinguishable from "this cell is deliberately unbanded":

```
D: missing label
  want=''
  reached end of check for 25bp_1G_c2 (fail=0)     # silently skipped, job green
```

Three ways Item B ships inert with the `BOOTSTRAP` warning **removed**, i.e. with the log now
positively implying the band is live:

1. A label typo in the TSV (`1G_c4_fastq` for `1G_c4_fastqc`) — one cell silently unbanded.
2. Someone renames a cell in `ci.yml` later without touching the TSV — same, and this is *likely*,
   because A1 is itself a precedent for editing the cell list.
3. **The baseline file is missing or renamed.** Verified: `grep -q '^# BOOTSTRAP' "$base"` on a
   missing file is inside an `if`, so `set -e` is exempt; `bootstrap` stays `0`; every `want` is
   empty; the job exits green. The only trace is one `grep: No such file` line. Before Item B this
   cost nothing. After Item B it silently discards the entire drift half of D5.

The plan's mitigation is Validation 5's *"confirm `want` is non-empty for all four rows"* — a one-off
manual `awk` at activation time. It does not survive the next edit to either file, which is exactly
the lifetime over which the file is supposed to gate.

**Fix — make it an in-job assertion, ~4 lines.** Two options; the second is stronger:

```bash
# in `measure`, replacing the bare `if [ -n "$want" ] …`
if [ "$bootstrap" = 0 ] && [ -z "$want" ] && [ "$gating" = 1 ]; then
  echo "::error::$label: no baseline row in $base and BOOTSTRAP is gone"; fail=1
fi
```

plus, in the step preamble, `test -f "$base" || { echo "::error::baseline file missing"; exit 1; }`.
The `gating=1` qualifier keeps the deliberately-unbanded c2 cells legal while making an accidental
omission on a gating cell loud. If you would rather not couple the two flags, hold an explicit
`banded` list and diff it against the labels actually measured.

This is worth doing inside Item B rather than deferring: Item B's entire purpose is to convert a
loudly-inert gate into a live one, and it currently swaps a printed warning for a silent failure
mode.

---

## Important

### I1. The escalation threshold is 1.13; the budget-breach line is 1.101

Marked: **arithmetic-verified** against `src/clump.rs` and `src/main.rs:209`.

A3's edge case: *"if `25bp_1G_c2` returns ≥ 1.13 … then the budget is breached on Linux too."*

The printed prediction at `--memory 1G` is 930 MiB and the budget 1024 MiB for **both** c2 and c4
(I recomputed the layout: c2 → `n_bins=16`, `denom=496`, `B=23,911,222`, pred `976,128,906 B` =
930 MiB; c4 → `denom=624`, `B=19,006,356`, pred 930 MiB). So

```
peak/budget = peak/pred × 930/1024 = peak/pred × 0.908203
peak/budget = 1.000  ⟺  peak/pred = 1.1011
```

The budget is breached at **peak/pred ≥ 1.101**, not 1.13. The plan has conflated two different
questions in one sentence: *1.13 is "macOS's value reproduced"* (it is macOS's second-lowest rep),
whereas *1.101 is "the budget is breached"*. There is a live ~3-point window, `[1.101, 1.130)`, in
which Linux breaches the budget and the plan's stated escalation criterion says nothing — and
Item A exists specifically to decide #457's urgency, so mis-stating the breach line defeats the cell.

Two thresholds, both worth naming explicitly:

| Ratio | Meaning |
|---|---|
| ≥ 1.101 | **budget breached on Linux** → #457 escalates to "reproduces on both platforms" |
| ≥ 1.13 | macOS's magnitude reproduced → the platform offset does not apply at c2 at all |

### I2. B1's "four reps" needs a code-identity precondition it cannot verify in advance

Marked: **code-verified for today, unverifiable for the future.**

Assumption 2 verifies that reps 1–2 measured identical code. B1 then adds reps 3 and 4 — from the
force-push and the merge — and never applies the same bar to them, although that is the assumption
B1 actually rests on. It holds today, and I checked:

```
git diff --name-only 1b2579e origin/dev -- src Cargo.toml Cargo.lock   →  (empty)
```

So the rebase changes nothing the binary is built from, and reps 1–4 are genuinely comparable. But
this is a property of `dev` *at this moment*. `wt-450` is sitting at `4f7b907` on
`fix/450-inert-clip-refusal`; if #450 (or anything else) merges into `dev` before #452, rep 3/4's
binary differs from rep 1/2's, and the "four reps" mean silently averages two binaries. Nothing in
the plan would notice.

**Fix.** Record the commit each rep came from in the TSV's range comment (the plan already commits to
a per-row comment, so this is free), and make the harvest step in B assert
`git diff --quiet <rep1-sha> <rep-n-sha> -- src Cargo.lock` for each rep it uses. If a rep fails that,
drop it and take another push — cheap, since reps keep arriving (see O5).

### I3. Validation 2 checks the wrong ref, so the assertion the plan calls essential is vacuous

Marked: **code-verified.**

Validation 2 is billed as *"the regression the rebase exists to prevent — assert it, don't assume it"*
and its mechanism is `after A: git show origin/dev:CHANGELOG.md | grep -c 457`. Run after A but
before #452 merges, `origin/dev` is untouched by anything A did, so it passes no matter how badly the
rebase went. The check that catches a botched rebase is on the **branch** side:

```bash
git diff --name-only origin/dev...HEAD    # must NOT list CHANGELOG.md
git diff --quiet origin/dev HEAD -- CHANGELOG.md && echo "CHANGELOG untouched by the branch"
```

The outline's step 3 already lists the four expected files, which is the real assertion — so the
coverage exists but Validation 2 points at the wrong thing. Retarget it, and keep the
`grep -c 457` form as a *post-merge* check on `dev`.

### I4. "+0.05 is ≥ 3× the widest observed spread" is 2.94×, and the two claims that actually matter are stronger

Marked: **arithmetic-verified.**

The widest spread in the plan's own table is 1.7 points (`25bp_1G_c4`). `0.05 / 0.017 = 2.94`, so
"≥ 3×" is false on the plan's own definition, and it is asserted twice — in Context §1 and again in
B2 as the criterion a future editor should re-check. A future editor applying a false criterion is
worse than one applying none.

Replace it with the two sensitivity numbers, which are the ones a reader needs and which are more
reassuring than "3×". Using mean-as-centre and half-range as the excursion unit:

| cell | mean | half-range | false positive needs | missing a minimum lever-A (+0.078) needs |
|---|---:|---:|---:|---:|
| `1G_c4` | 0.776 | 0.0060 | 8.3× excursion | 4.7× downward |
| `1G_c4_fastqc` | 0.699 | 0.0020 | 25× | 14× |
| `25bp_1G_c4` | 0.932 | 0.0085 | 5.9× | 3.3× |
| `multipair_c4` | 0.790 | 0.0025 | 20× | 11× |

Both failure directions need a 3.3×-or-larger excursion beyond anything observed, on the worst cell,
against Phase 3 §5.3's *minimum* lever-A delta of 7.8 points. That answers the lead's question 4
affirmatively and it answers risk (a) with a number instead of a hedge.

### I5. "Both non-gating (`0`)" and Validation 4 overstate the isolation — five other paths can fail the job

Marked: **code-verified** (and the risk turns out to be nil today, which is the interesting part).

`$gating` guards exactly one line: `test "$peak_kb" -le "$pred_kb"`. A cell with `gating=0` still
sets `fail=1` on: a non-zero exit code; a missing `^clumpify:` banner; output that is not `*.fq.gz`;
an empty parse field; and the `pred == budget × 10/11` margin pin. So Validation 4's *"the new cells
cannot fail the job"* is true only of the ratio, and risk (c)'s *"they detect nothing"* is wrong in
the plan's own disfavour — the cells **do** newly guard the plain-mode fallback, the gzip promise and
the margin at `--cores 2`, which is real added coverage the plan is declining to claim.

I checked the two paths that could plausibly bite at `1G/c2` and both are safe:

- **Margin pin.** `expect_kb = 953,250`, `pred_kb = 952,320`, diff `−930 KiB`, inside the ±2048 KiB
  slack. It passes at c2 because `denom = 496` is divisible by 16, exactly as at c4 (`624`).
- **Plain-mode fallback** (`main.rs:195`, which would kill the `^clumpify:` check). The floor at c2
  is `(224 MiB + 31 MiB) × 11/10 ≈ 282 MiB`, far below 1024 MiB. Read length does not enter the
  floor, so the 25 bp cell is equally safe.

So this is a wording fix, not a defect: say "the ratio check is off" rather than "non-gating", and
move the five surviving checks from risk (c) into A1 as a benefit.

### I6. The plan does not say *where* to run the rebase, and the obvious place is a tree three agents are reading

Marked: **unverified as a hazard; verified as a shared tree.** `git worktree list` shows this
checkout at `3f1a930` plus `wt-450` and `wt-450-base`, and the brief states another process shares
this tree. The outline's step 1 is a bare `git rebase --onto …`, which rewrites the working tree and
`HEAD` under whatever else is reading it. Memory already records that a branch switch mid-verification
silently invalidated agent results with identical-looking test counts; a rebase is a strict
superset of that.

**Fix:** one line in the outline — do the rebase in a dedicated worktree, or after confirming no
agent is mid-verification. Note that `git worktree remove` needs a sandbox bypass here.

---

## Optional

**O1. `--force-with-lease` with no argument is weaker than it looks in this specific sequence.**
The bare form uses the remote-tracking ref as the expected value, and step 1's `origin/dev` reference
requires a `git fetch` — which refreshes that very ref. So a fetch immediately before the push
re-baselines the lease against whatever is on the remote *now*, defeating the protection it exists
for. You know the expected tip; pin it: `git push --force-with-lease=ci/clumpify-memory-guard-439:3f1a930`.
(Marked: mechanism-verified from git's documented semantics; I did not push anything.) The branch has
`origin/ci/clumpify-memory-guard-439` configured as upstream and `push.default` is unset (`simple`),
so a bare `git push` does resolve correctly — that part is fine.

**O2. `51bp_1G_c2` breaks the existing label convention, and labels become committed baseline keys.**
The four current 51 bp cells carry no length prefix (`1G_c4`, `multipair_c4`); only the 25 bp cell is
prefixed. Phase 3 §5.3 calls this cell `1G/c2`. Choosing `1G_c2` keeps the convention and matches the
parent plan's prose; choosing `51bp_1G_c2` is defensible for explicitness but then the two existing
51 bp cells are misnamed. Either way, settle it in A1 rather than in B — a rename after B has written
baseline rows silently unbands the cell (see C2).

**O3. Two of B2's four placeholders sit exactly on a rounding tie.** `(0.941+0.924)/2 = 0.9325` and
`(0.792+0.787)/2 = 0.7895`; the plan writes `0.933` and `0.790`, i.e. half-up. Harmless for
placeholders, but B2's recomputation instruction should name the convention, because a ±0.001 choice
on a committed threshold is exactly the class of detail this arc keeps losing. (Marked:
arithmetic-verified.)

**O4. A4's "if the job lands over ~12 min, raise the timeout" has no referent.** `timeout-minutes`
is 30, so 12 minutes is not a pressure point. And the estimate is conservative: from the Actions API
for run `32065624366`, `Measure and gate` took **2m19s for four cells** (build 57s, fixtures 26s,
total 3m59s). Two c2 cells at roughly double a c4 cell's wall time is **+1.5 to 2 min**, landing the
job near 5m30s. Say "the 30-minute timeout has ~5× headroom; revisit only if the job passes ~20 min"
and drop the 12. (Marked: verified from the Actions API.)

**O5. More than four reps will exist for free, and Assumption 5 is the plan's own stated weakness.**
The job runs on every PR targeting `master`/`dev` and every push to them, so by the time B is written
there will be additional reps from unrelated PRs at no cost. Use them — with I2's code-identity
filter, and preferring runs on `push: dev` (which measure a merged binary) over PR runs that touch
`src/parallel.rs` or `src/clump.rs`. n=6–8 turns Assumption 5 from "bounds the spread" toward
"estimates it" for zero effort.

**O6. Validation 7 needs a throwaway *PR*, not just a throwaway branch.** Triggers are
`push: [master, dev]` and `pull_request: [master, dev]`, so a branch pushed without a PR does not run
the job at all. A draft PR targeting `dev` is the cheapest form. One clause in Validation 7 saves the
next person a confused ten minutes. (Marked: code-verified from `ci.yml` lines 3–15.)

**O7. Stale header.** "**Baseline:** `dev` at `4f7b907`; PRs #448, #451, #452 open" — `dev` is
`242546a` and #448/#451 are merged. The body is written correctly for the post-merge world ("after
#451 merges, not before"), so this is only a header refresh, but a plan whose stated baseline is two
merges out of date invites a reader to re-verify everything downstream of it.

**O8. "`n_bins = 16` at c2" is true but not distinctive.** `n_bins = max(16, 4·cores)` gives 16 at
**both** c2 and c4; it first exceeds 16 at `--cores 5`. What makes `B` largest at c2 is the smaller
`denom` (496 vs 624) from fewer workers, with the pool's share rising from 59% to 74%. The mechanism
in A1 is right; the stated cause is the clamp, which c4 shares. (Marked: arithmetic-verified.)

---

## Answers to the eight questions

**1. The rebase.** Correct, and I upgraded the plan's "guaranteed clean" from an incomplete argument
to a proof — see the header section. `--force-with-lease` is *nearly* sufficient; pin the OID (O1).
No silent-drop or duplication failure mode: `git rebase` drops upstream-duplicate commits by
patch-id, and none of the four `ci(clump):` patches exists on `dev`. The tracked tree is clean, so
the rebase will not refuse, and no replayed path collides with an untracked file. The real unstated
hazard is *where* it runs (I6).

**2. The two-condition band claim.** Correct, code-verified against the job's actual shell, and I
also confirmed the awk lines survive `set -e` — worth knowing, because with `bootstrap=1` in both
measured runs **this code has never executed in CI**. There is a **third** inertness condition the
plan does not identify (a missing baseline file), and a fourth (a label typo) — C2.

**3. The arithmetic in A3.** The four peak/budget conversions are all correct (×0.908203). The
escalation threshold is wrong: 1.101, not 1.13 (I1). Two smaller issues with the offsets:

- The macOS input ranges are narrower than the corpus supports. Phase 3 §0.4 lists six 25 bp reps at
  c2 spanning **1.119–1.164**; the plan quotes 1.127–1.164, dropping reviewer A's lowest. Likewise
  25 bp c4 is **1.054–1.075** across both harnesses, not 1.069–1.075. The offset is therefore
  **11–15 points**, not 13–15.
- The 51 bp prediction blends the 25 bp offset with the 51 bp one (0.84 uses ~17, 0.88 uses ~13),
  which widens the interval with a length-specific quantity applied to the wrong length. And "~17"
  is the top of the 15.7–16.9 range, not its centre.

The **conclusion survives all of this**, which is the important part. Widening the inputs to the full
evidence gives 25 bp/c2 ≈ 0.97–1.05 → peak/budget 0.88–0.96: still no breach. And an independent
route agrees. macOS decomposes almost perfectly additively — the c4→c2 core increment is +0.070 at
51 bp and +0.074 at 25 bp, i.e. length-independent — so:

```
Linux 25bp/c2 ≈ 0.932 (measured c4) + core increment
  core increment = macOS's +0.070            → 1.002  (peak/budget 0.910)
  core increment scaled by Linux's 1.225×
  larger read-length effect (+0.086)          → 1.018  (peak/budget 0.925)
```

Both land inside the plan's stated interval, and the additive and multiplicative platform models
agree to within 0.01 — so A3's answer is robust to the modelling choice it never states.

That yields the sharpest falsifiable form of A3, which I recommend substituting for the current
prose: **a Linux budget breach at 25 bp/c2 requires a core increment of +0.169, i.e. 2.4× macOS's
+0.070.** That is a claim a single number can refute.

**4. Four reps, and mean-not-max.** Defensible on both counts. Four is the top of the TSV's own
"two or three more runs" bar. Mean is right, and I4 gives the numbers that justify it better than the
plan does: both false-positive and false-negative directions need ≥3.3× the observed excursion. The
real weakness is not *n* but the reps' clustering in time and the code-identity gap (I2); O5 fixes
the first two for free.

**5. Comparability under appending.** Yes — code-verified. `measure` calls are sequential in one
shell, peak RSS is per-process via `ru_maxrss`, and cells 1–4 complete before the new cells start, so
causality alone protects them. No shared disk or page-cache state flows backwards. Assumption 3 is
sound *for cells 1–4*, which is what B2 needs. The comparability risk that does exist is in the
binary, not the cell order — I2.

**6. Validation sufficiency.** Two gaps rise to Critical (C1, C2). Validation 7 is correctly
identified as the one that decides whether D5 bought what it cost, and correctly flagged as never
run; O6 fixes its mechanics. Lever A's minimum delta (+7.8 points from §5.3's range, not the
headline +11) does clear the +0.05 band with ~2 points to spare even from the bottom of the observed
spread, so Validation 7 should pass — but "should" is the word risk (b) exists to reject, so it must
not be dropped. Beyond C1/C2 I found nothing high-risk unguarded: the exit-code, disk, engagement,
gzip and parse clauses are all present and fire regardless of `gating`.

**7. Two c2 cells or one — justified, and the plan under-argues it.** Keep both. But the stated
reason is weaker than the true one. "A single 25 bp/c2 number cannot separate read length from core
count" understates what is already known: the read-length effect at fixed cores is *already measured
on both platforms* at c4 — 0.941 vs 0.770 = +17 points on Linux, +13 on macOS. What is genuinely
unmeasured on Linux is the **c4→c2 core increment**, against a macOS reference of +0.070 that is
length-independent. `51bp_1G_c2` measures precisely that missing term, at a read length where macOS
shows no budget breach, so it is a clean control. Without it, a single 25 bp/c2 number confounds
"Linux's core increment differs" with "Linux's read-length effect differs at c2" and cannot support
the Phase 4 scoping decision the plan wants it for. Secondary payoff: `51bp/c2` is Phase 3's
gating-candidate-that-was-too-tight (macOS 1.009, 1 of 12 reps over), so the Linux number also tells
you whether it could ever be gated. Rewrite A1's third bullet around the core increment; the
decision itself is right.

**8. Overconfident or unverifiable claims.** Collected: "≥ 3× the spread" (I4, false at 2.94×);
"≥ 1.13 → budget breached" (I1, wrong line); "the expected finding is that #457 breaches the
prediction but not the budget" — the plan's own interval is 0.98–1.03 and straddles 1.0, so
"breaches the prediction" is true of only the upper half of its own range, and the additive route
puts the centre at 1.00–1.02; "both non-gating, so they detect nothing" (I5, understates); "the fix
is guaranteed clean" (true, but not by the reason given); "~17 points" at 51 bp (16.3 by the two-rep
mean); the macOS ranges narrowed by dropping reviewer A's reps (question 3); and the stale header
(O7). None of these changes a decision except I1.

---

## What I did not check

- **No measurement.** I ran no `trim_galore` invocation, created no fixture, and produced no new
  peak-RSS figure. Every ratio here is arithmetic on figures already in the corpus.
- I did not verify the **macOS peak figures** in Phase 3 §0.4 (1048.1 etc.) — I took them as given and
  only recomputed ratios from them.
- I did not verify Phase 3's **lever-A delta** of 72–104 MiB / 7.8–11.2 points; I used it as stated.
- I did not run the rebase, push anything, or create a worktree. C1/C2's fixes are untested code
  sketches — the band logic I ran was the *existing* clauses extracted verbatim, not my replacements.
- I did not check whether **#450** or any other open PR touches the clump memory path; I only
  observed that its worktree exists and that a merge ahead of #452 would break reps 1–2's
  comparability with reps 3–4 (I2).
- I did not confirm the **c2 cells' wall time** empirically; the +1.5–2 min in O4 is a 2× scaling of
  the measured c4 cells, not a measurement.
