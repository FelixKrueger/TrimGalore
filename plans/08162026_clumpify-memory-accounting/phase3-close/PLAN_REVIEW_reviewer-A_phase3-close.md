# Plan review — Reviewer A — `phase3-close/PLAN.md`

**Target:** `plans/08162026_clumpify-memory-accounting/phase3-close/PLAN.md`
**Reviewed at:** repo `/Users/fkrueger/Github/TrimGalore`, branch `ci/clumpify-memory-guard-439` @ `3f1a930`, `origin/dev` @ `242546a`
**Date:** 2026-08-17
**Measurement performed:** none. No `trim_galore` invocation, no fixture generated, no tracked file touched. All numeric findings are arithmetic on figures already in the plan / parent plan, or on constants read out of `src/clump.rs`.

Every finding is tagged **[code-verified]**, **[arithmetic-verified]** or **[unverified]**. Section 0 answers the eight questions put to me; sections 1–5 are the review proper; section 6 is the prioritised action list.

---

## 0. Direct answers to the eight questions

**0.1 — Is the rebase correct and genuinely conflict-free? Is `--force-with-lease` sufficient?**

The rebase is correct, and it is conflict-free by *proof* rather than by simulation. **[code-verified]**

```
git diff --name-only 1b2579e 44dc732            → CHANGELOG.md          (plan's claim holds)
git diff --name-only 1b2579e origin/dev         → CHANGELOG.md
                                                  docs/…/modes/clump-only.md
                                                  docs/…/performance/clumpy.md
git diff --name-only 1b2579e origin/ci/…-439    → .github/clumpify-memory-baseline.tsv
                                                  .github/workflows/ci.yml
                                                  examples/mk_clump_fixture.rs
                                                  src/clump.rs
```

The set of files the new base changed and the set the four replayed commits touch are **disjoint**. A rebase 3-way-merges each commit against its own parent, so every touched file has an identical base and target — there is nothing for git to reconcile. This is stronger than the plan's `merge-tree` simulation, which tests a merge rather than a replay. Worth putting in the plan, because it converts "verified clean by simulation" into "cannot conflict".

No silent-drop or duplicate risk: `git rebase` drops commits whose patch-id already exists upstream, and none of the four patches touch any file `1b2579e..origin/dev` modified, so none can match. `1b2579e..HEAD` is exactly the four `ci(clump):` commits **[code-verified]**. Local and remote branch tips agree (`3f1a930` both), the tree has no modified tracked files, and the branch is checked out in the main worktree — so no branch-switch or autostash hazard. **[code-verified]**

`--force-with-lease` is **not** sufficient as written — see **I3**. Two gaps: the bare form takes its expected value from the remote-tracking ref, which a `fetch` by the other process sharing this tree would silently update; and step 1 never fetches, so `origin/dev` could be stale and the rebase would land on the wrong base with no assertion catching it (see **I2**).

**0.2 — The two-condition band claim.**

Correct, and the sharper half of it is correct too. **[code-verified]** against the job's shell on the branch and the TSV on the branch:

- The guard is `if [ -n "$want" ] && [ "$bootstrap" = 0 ]`.
- `want=$(awk -v c="$label" '$1==c {print $2}' "$base")`. The TSV on the branch is 31 lines, **all** beginning `#`, so `$1` is `#` or `#word` for every line and matches no label. `want` is empty for all four cells.
- `bootstrap=1` because `grep -q '^# BOOTSTRAP'` matches TSV line 14.

Deleting `# BOOTSTRAP` alone sets `bootstrap=0` **and removes the `::warning::`**, so the job stops announcing that the band is off while `[ -n "$want" ]` still short-circuits it. The plan's phrase "inert *and looking activated*" is exactly right, and the disappearance of the warning is what makes it dangerous. B3's line numbers (TSV 14–17) are also correct.

One residual the plan does not name: the sentinel is `^# BOOTSTRAP`, so if B3's "retitled provenance block" happens to begin `# BOOTSTRAP history …` the sentinel survives and the band stays inert. B3 should assert `! grep -q '^# BOOTSTRAP'` after editing.

**0.3 — The arithmetic in A3.**

Mostly sound; the peak/budget conversion is exact; the escalation threshold is wrong. **[arithmetic-verified]**

Verified correct:

| claim | check | verdict |
|---|---|---|
| spreads 1.2 / 0.4 / 1.7 / 0.5 pts | differences of the tabulated reps | ✔ all four exact |
| 25 bp offset 13–15 pts | 1.069−0.941 = 0.128, 1.075−0.924 = 0.151 | ✔ (12.8 rounds to 13) |
| `25bp_1G_c2` → 0.98–1.03 | 1.127−0.15 = 0.977; 1.164−0.13 = 1.034 | ✔ exact |
| peak/budget = ratio × 10/11 | prediction is `budget×10/11` by construction, and the job pins it to ±2048 KiB | ✔ code-verified |
| `25bp_1G_c2` → 0.89–0.94 | 0.98×10/11 = 0.891; 1.03×10/11 = 0.936 | ✔ |
| `51bp_1G_c2` → 0.76–0.80 | 0.84×10/11 = 0.764; 0.88×10/11 = 0.800 | ✔ |
| macOS 25 bp/c2 = 1.127–1.164 | handoff §2 gives 1.02–1.06× *budget*, 7 reps; ×11/10 → 1.122–1.166 | ✔ consistent |
| macOS 51 bp/c2 = 1.009, 1-of-12 | phase3 `PLAN.md:440` verbatim | ✔ faithful |

Wrong or under-derived:

- **The escalation threshold should be 1.10, not 1.13.** I computed `1024/930.91 = 1.1000` and `1024/930 = 1.1011`. Peak exceeds the *budget* the moment peak/pred exceeds ~1.100. The plan's ≥1.13 trigger is derived from "macOS's number reproduces unchanged", which is a different question, and it leaves the window **[1.100, 1.13)** as a real Linux budget breach that the plan's own decision rule would read as "below threshold, no escalation". See **I1**.
- **The 0.84–0.88 band for `51bp_1G_c2` is not derived from the offset the plan states.** The text says the 51 bp offset "is ~17"; applied consistently that gives 1.009 − 0.169 = **0.840**, a point estimate, not 0.84–0.88. The 0.88 end requires the *25 bp* offset (1.009 − 0.129). The plan is silently taking the union of both offset ranges. Harmless in direction (it predicts "no breach", and a lower value is more so) but it hides which offset is being trusted. See **I10**.
- **The two macOS anchors are different statistics.** 1.009 is the *max of 12 reps*; 1.127–1.164 is the *full range of 7*. The plan then subtracts them to get "a 12–15 point read-length effect at fixed cores" — max-vs-range understates the effect, and that difference is the quantity the 51 bp control exists to measure on Linux. See **I10**.
- Minor: "at 51 bp it is ~17" uses rep 1 alone (0.939 − 0.770 = 16.9); the reps 1–2 mean the sentence claims to use gives 0.939 − 0.776 = 16.3. **[arithmetic-verified]**, cosmetic.

**0.4 — Is four reps defensible, and is mean the right statistic?**

**Yes to both, but the plan's stated reason for mean is the weaker of the two available, and "full sensitivity" is overstated.** **[arithmetic-verified]**

Treating each 2-rep range as an estimator of σ (E[range₂] = 1.128σ) and comparing a fresh run against a mean of four reps (sd = σ√1.25):

| cell | range | σ̂ | sd(new − mean₄) | upper trigger at +5.0 | lever A weakest form (7.8) |
|---|---:|---:|---:|---|---|
| `1G_c4` | 1.2 | 1.06 | 1.19 | z = 4.21, p ≈ 1e-5 | z = 2.35, miss ≈ 0.9 % |
| `1G_c4_fastqc` | 0.4 | 0.35 | 0.40 | z = 12.6 | z = 7.1 |
| `25bp_1G_c4` | 1.7 | 1.51 | 1.68 | z = 2.97, **p ≈ 0.15 %** | z = 1.66, **miss ≈ 4.8 %** |
| `multipair_c4` | 0.5 | 0.44 | 0.50 | z = 10.1 | z = 5.7 |

Two things follow that the plan does not say:

1. **Mean is not merely preferable, it is necessary.** The parent plan puts lever A's weakest form at **7.8 points** (`phase3/PLAN.md:432`, from the 72.5 MiB variant), not the plan's "8–11". Against a +5.0 trigger that leaves **2.8 points** of detection margin. Max-centring would move the centre up by roughly the observed spread (~1.7 pts on the worst cell) and spend 1.7 of those 2.8 points, taking the weakest-form miss rate from ~5 % to ~23 % on that cell. So the plan's conclusion is right and its margin is thinner than "full sensitivity" implies.
2. **A false red is likely once over the job's life, and it will land on the 25 bp cell.** ~0.15 % per run → 14 % over 100 runs, 26 % over 200. §5.4 calls this the most expensive failure mode. The per-row range comment is the right mitigation and should be joined by a one-line note in the TSV saying *what to do* when a single cell reds by 5–7 points on unmodified code (re-run once; if it clears, it was spread).

The redundancy across four cells is what makes 2.8 points adequate — lever A moves every cell together, so the job reds if any one fires. That is a further argument for keeping all four calibrated, and it is worth stating because it is currently accidental rather than designed. Caveat: the four cells share one VM and one allocator history, so their deviations are correlated to an unknown degree; the independence that redundancy relies on is **[unverified]**.

**0.5 — Does appending really preserve comparability?**

On the axis the plan names, yes — and correctly reasoned: each `measure` is a separate process, peak RSS is per-process, and both new cells run after all four existing ones, so no earlier cell can be affected. **[code-verified]** against the job's shell.

But the plan checks only the cell-order axis, and reps 3–4 differ from reps 1–2 on **two more axes**:

- **The base changes.** I verified this one holds: `git diff --name-only 1b2579e origin/dev -- src/ examples/ Cargo.toml Cargo.lock build.rs tests/` is **empty**, so the post-rebase release binary is behaviourally identical to the one reps 1–2 measured (`44dc732` is CHANGELOG-only; #448 is docs-only). **[code-verified]** — the conclusion is safe, but a reader cannot get there from the plan, and had `dev` moved on `src/` the four-rep mean would have pooled two binaries silently. Make it a stated precondition (**I8**).
- **The runner image.** `ubuntu-latest` is a rolling label. Assumption 6 treats image drift as a multi-week risk, but reps 3 and 4 are days after reps 1–2, so it applies *inside* the four reps used to set the centre. Cheap fix in **I8**.

**0.6 — Validation sufficiency, especially Validation 7.**

The highest-risk thing in Item B is unguarded: **nothing permanent checks that the band is actually armed** (see **C1**). Validation 5 catches it once, at activation, by hand — and B3 explicitly defers the c2 rows to "a later pass", which guarantees a later hand-edit of this file by someone who will not re-run Validation 5.

Validation 5 itself tests the wrong half: perturbing a row by **+0.20** fires only the *lower* bound. I ran the job's two awk lines directly:

```
ratio=0.776  want=0.976 (row +0.20)   upper=no                lower=FIRES
ratio=0.776  want=0.576 (row −0.20)   upper=FIRES             lower=no
ratio=0.854  want=0.776 (lever A min) upper=FIRES             lower=no
```

So Validation 5 as written proves the mechanism is live but never exercises the upper comparison — the half the plan itself calls load-bearing. **[verified by direct execution of the job's awk]** Fix: perturb −0.20 (or both directions). Mitigating note: I confirmed above that the upper line fires correctly for lever A's magnitude, so the residual risk is small and the fix is cheap insurance rather than an open hole.

Validation 7 is correctly identified as the one that decides whether D5 bought what it cost, and the plan is right to refuse to drop it. Two things stand in its way that the plan does not mention:

- **It cannot be run on a throwaway branch by pushing.** The workflow's triggers are `push: [master, dev]`, `pull_request: [master, dev]`, `schedule` — there is **no `workflow_dispatch`** **[code-verified]**. A throwaway branch must be opened as a (draft) PR against `dev` to make the job run at all. This friction is a plausible reason Phase 3's Validation 12 was never run, so naming it is worth a line (**I5**). Same constraint applies to Validation 5.
- **Its arithmetic does work out.** Baseline 0.776, trigger 0.826; lever A's forms land at 0.854 / 0.888. Red both times, by 2.8 / 6.2 points. **[arithmetic-verified]** So Validation 7 should pass — it is worth predicting that in the plan so a red-for-the-wrong-reason is distinguishable.

Also unguarded: **Item A produces no durable artefact** (**I6**). The two c2 cells print to a job log that expires, B gives them no TSV row, and no step records the answer on #457 — yet "decide how urgent #457 is" is Item A's stated purpose.

Validation 4's *claim* ("The new cells cannot fail the job") is false as written — five paths in `measure` set `fail=1` irrespective of `gating` (**I7**). Its *expectation* is correctly narrow.

**0.7 — Two c2 cells or one?**

**Keep both, but fix the justification — it is weaker and more specific than the plan says.**

The plan asserts "a single 25 bp/c2 number cannot separate the two candidate causes". That overstates it. The runner already has both c4 cells, so:

- read length at c4 on Linux: 0.776 vs 0.933 = **15.7 points** — already measured;
- cores at 25 bp: `25bp_1G_c2` − `25bp_1G_c4` — the one new cell gives this.

Those two contrasts already decompose the **main effects**. What only `51bp_1G_c2` supplies is the **interaction** — whether the read-length effect at c2 matches the one at c4 — which is precisely the "if Linux reproduces that *difference*" test the plan wants, but is not the same as "cannot separate the two causes". Completing a 2×2 for the last cell is a reasonable ~2 min/run, and the plan already surfaces the recurring cost honestly, so I would keep it. Not scope creep; just mis-argued. Note too that on the numbers already in hand both platforms show ~13–16 points of read-length effect at c4, so the mechanism arguably *already* transfers at c4 and the marginal value of the fourth cell is narrower than "#457's mechanism transfers" suggests.

**0.8 — Claims that are wrong, unverifiable, or over-confident.**

`C1`, `I1`, `I2`, `I4`, `I7`, `I9`, `I10`, `O1`, `O2`, `O3`, `O4` below. The one flatly wrong mechanism claim is **O1** (`n_bins` is 16 at c4 as well, so it is not what makes c2 different).

---

## 1. Logic

The two-item split at a merge boundary is right and its two dependencies are both real and both stated: reps (B needs runs A produces) and cell order (B's centres must be calibrated against the final sequence). The A-before-B ordering argument is sound.

The rebase rationale is the strongest part of the plan. It correctly narrows "the diff re-applies its content" to the actual hazard — a `CHANGELOG.md` conflict whose resolution toward #452 reverts `44dc732` and guts D7's deliverable — and I confirmed each supporting fact. The one place the logic does not close is that the *verification* of that hazard (Validation 2) does not test the branch (**I2**).

The A3 prediction table is good practice — recording a falsifiable prediction before measuring is exactly right, and Assumption 7 states it explicitly. The flaw is in the decision rule attached to it (**I1**), not in the predictions.

The Item B ordering has one unstated consequence: B's own PR run is rep 5, judged against a centre built from reps 1–4. That is correct as a design (it is the first genuinely independent test of the band) and B4 says so, but the plan should say the quantity being tested is "does a fifth draw fall within +5/−10 of the mean of four" — with the ~0.15 %/run false-red rate from §0.4, a red on B4 is not automatically evidence that "the band width or the rep count is wrong" as B4 asserts; it could be one draw in the tail. B4's instruction should be "re-run once; if it clears, record it as spread and revisit the width".

## 2. Assumptions

| # | Plan's status | My verdict |
|---|---|---|
| 1 | Fixed — `--clumpify` needs `--cores ≥ 2` | ✔ **[code-verified]** `src/cli.rs:747-751` bails. c2 is genuinely the lowest measurable cell. |
| 2 | Fixed — reps 1–2 measured identical code | ✔ **[code-verified]** only diff `69d8f67`→`3f1a930` is the TSV comment block. |
| 3 | Fixed — appending does not perturb earlier cells | ✔ on the axis named; **incomplete** — see §0.5 and **I8**. |
| 4 | Fixed — `f51`/`f25` already generated | ✔ **[code-verified]** both at 10 M pairs in the "Generate fixtures" step. |
| 5 | Open — four reps enough for a centre | ✔ honest, and §0.4 quantifies it. Add the derived false-red rate. |
| 6 | Open — runner allocator stable over weeks | ✔ but under-scoped: it bites *within* the four reps too (**I8**). The "all four cells move together" signature is a genuinely good discriminator and worth keeping. |
| 7 | Not assumed — #457 reproduces on Linux | ✔ correctly framed. |

**Two assumptions the plan makes without listing:**

- **A8 (implicit): the post-rebase binary is the same binary reps 1–2 measured.** True — I verified it — but load-bearing for B1's four-rep pooling and never checked in the plan.
- **A9 (implicit): `--clumpify --cores 2 --memory 1G` clears the layout floor and engages.** If it did not, `main.rs:189` warns and falls back to plain mode, no `clumpify:` banner is printed, and `measure` sets `fail=1` **regardless of `gating=0`** — a "non-gating" cell would turn the PR-blocking check red. It does clear: `n_bins = max(16, 4×2) = 16`, `denom = 23×16 + 64×2 = 496`, `dyn_min = (1 MiB × 496)/16 = 31 MiB`, floor `= (224 + 31) × 11/10 ≈ **281 MiB**` against a 1024 MiB budget. **[arithmetic-verified]** from `src/clump.rs:38,174,244,290`. And c2's floor is *lower* than c4's (289 MiB), so a cell that already passes at c4 cannot fail this way at c2. Worth one line in the plan, because it is the only way an appended non-gating cell could break the job.

## 3. Efficiency

Proportional and honestly costed. +~4 min/run is a safe upper bound: the c2 cells roughly double per-cell wall time against c4, and the four existing cells are a minority of the current 3m59s (which also carries a cargo build and two 10 M-pair fixture generations). No new fixture, no new build, no new dependency. **[unverified]** — I did not time anything.

Disk is unchanged by appending: `measure` does `rm -f "$out"/*.fq.gz` per cell, so peak disk stays `fixtures + one cell's output`. **[code-verified]** Two notes: the 8 GiB assertion runs *before* fixture generation, so it does not assert free disk at measurement time (**O6**), and the residual per-cell `*_trimming_report.{txt,json}` files are not cleaned but are negligible.

`timeout-minutes: 30` has ample headroom; A4's advice to raise it at ~12 min is incoherent (**O2**).

## 4. Alternatives worth weighing

1. **Assert the band is armed, in the job, permanently** (see **C1**) instead of relying on a one-time manual check. This is the highest value-per-line change available and it directly closes the trap the plan's own Context section identifies.
2. **Asymmetric centring** — mean for the lower bound, max-of-n for the upper — was presumably considered and rejected. §0.4 shows rejecting it is right: it would spend 1.7 of the 2.8 available detection points. Say so with the number, since "taking the max would spend sensitivity" is currently an assertion.
3. **Widen the upper bound to +0.06 and drop the false-red rate by ~5×.** Detection margin against lever A's weakest form falls 2.8 → 1.8 points (miss ≈ 14 % on the 25 bp cell, from 5 %). Given four correlated cells fire together, this trade may well be worth taking — but it needs the correlation measured first, so it is a post-activation decision, not one for this plan. Worth recording as a lever with its price attached.
4. **Give the 25 bp cell its own wider band** rather than one global ±. It is the cell with 3–4× the spread of the others *and* the one #457 is about. A per-row margin column would let the tight cells stay tight. Costs a format change; probably not now, but it is the natural next move if the 25 bp cell flakes.
5. **Record the c2 numbers on #457 as a comment** rather than only in the TSV — the issue is where the urgency decision gets made, and comments do not expire (**I6**).

## 5. What I did not check

- Any wall-time or memory measurement. No `trim_galore` run, no fixture.
- Whether the actual CI run logs contain what the plan says they do (I took the reps 1–2 table and run IDs on trust; the *commits* they name I did verify).
- #457's and #452's GitHub state (no `gh` calls); I used the pinned state in the brief.
- The correlation between the four cells' run-to-run deviations — needed to know how much the four-cell redundancy really buys.
- Whether reps 3–4 will in fact land on the same `ubuntu-latest` image as reps 1–2.

---

## 6. Action items

### Critical

**C1 — Item B ships a gate with no self-check that it is armed. [code-verified]**
`want` empty → the band silently skips that cell, with no error and (post-`BOOTSTRAP`) no warning. A typo'd label, a renamed cell, or a provenance line that lost its `#` all produce a green job with the gate off. Validation 5 catches this at activation only, and B3 *deliberately defers* the c2 rows to "a later pass" — so a later hand-edit by someone who will not re-run Validation 5 is planned, not hypothetical. Add a permanent assertion after the `measure` calls, e.g.:

```bash
# Every cell with a committed baseline must have been compared against it.
if [ "$bootstrap" = 0 ]; then
  for label in 1G_c4 1G_c4_fastqc 25bp_1G_c4 multipair_c4; do
    awk -v c="$label" '$1==c {n++} END{exit !(n==1)}' "$base" \
      || { echo "::error::$label has no single baseline row - the drift band is not armed for it"; fail=1; }
  done
fi
```

This also catches the duplicate-row case (**O5**) and makes B3's sentinel-survival risk visible. Without it, Item B's entire deliverable can be inert while every run is green — which is the exact failure the plan opens by describing.

### Important

**I1 — Escalation threshold: use 1.10, not 1.13. [arithmetic-verified]**
`1024/930.91 = 1.100`. Anything above ~1.100 breaches the budget on Linux. The window **[1.100, 1.13)** is a real breach that the plan's stated rule reads as "no escalation". Restate the edge case with all three bars named: **1.000** = the promise breaks (Phase 3's Validation 15 already calls that a finding); **1.100** = `--memory` is exceeded → #457 escalates to both platforms; **1.13** = macOS's number reproduced unchanged, i.e. the platform offset did not apply at c2. The plan's A3 table row "budget breach? **no**" is right for 0.98–1.03 but should carry the 1.100 line so a reader can place an out-of-band result.

**I2 — Validation 2 asserts on the wrong ref. [code-verified reasoning]**
`git show origin/dev:CHANGELOG.md | grep -c 457` tests whether #451 merged into `dev` — not whether the rebase preserved it *on the branch*. Rebase onto a stale `origin/dev` and this passes while the branch's CHANGELOG lacks the mechanism sentence. Assert on the rebased branch: `git show HEAD:CHANGELOG.md | grep -c 457` plus the "~50 bp and longer" qualifier, and add `git diff --quiet origin/dev HEAD -- CHANGELOG.md` (must be quiet: #452 must not touch the file at all).

**I3 — Pin the lease and the base. [code-verified]**
Bare `--force-with-lease` takes its expected value from `refs/remotes/origin/…`; the brief says another process shares this tree, and any `fetch` between the rebase and the push silently re-leases onto whatever is there now. Use `git push --force-with-lease=ci/clumpify-memory-guard-439:3f1a930`. Also add to step 1: `git fetch origin` immediately before, and assert `git rev-parse origin/dev` = `242546a` — nothing downstream catches a rebase onto a stale base.

**I4 — Validation 5 tests the lower bound only. [verified by executing the job's awk]**
`+0.20` fires the "fixture may have stopped filling the pool" branch, not the regression branch. Perturb **−0.20** to exercise the upper bound, or do both and state which message each is expected to produce. I confirmed the upper line fires correctly at lever A's magnitude, so this is cheap insurance, not an open hole.

**I5 — Validations 5 and 7 need a PR, not a branch push. [code-verified]**
Triggers are `push: [master, dev]`, `pull_request: [master, dev]`, `schedule`; there is no `workflow_dispatch`. "On a throwaway branch, on the runner" requires opening a draft PR against `dev`. Say so — this friction is a plausible reason Phase 3's Validation 12 never ran, and Validation 7 is the plan's own designated must-not-drop item. Optionally add `workflow_dispatch` to make it a one-command check in future.

**I6 — Item A produces no durable record of its own deliverable. [code-verified]**
The c2 ratios print only to a job log that expires; B3 gives the c2 cells no TSV row; no step reports them on #457. Add a step after reps 3–4 land: record both c2 cells in the TSV comment block (same treatment A2 gives rep 2) **and** post the numbers plus the A3 prediction outcome as a comment on #457. Otherwise "decide how urgent #457 is" ends as an unrecorded reading of a transient log.

**I7 — "The new cells cannot fail the job" is false. [code-verified]**
Five paths in `measure` set `fail=1` regardless of `gating`: non-zero exit, missing `clumpify:` banner, absent `*.fq.gz`, an empty parse field, and the `budget × 10/11` margin pin. Restate Validation 4 as "a c2 ratio above 1.0 must not fail the job" and add the reassurance that makes the addition safe: the c2/1G layout floor is **~281 MiB** (below c4's ~289 MiB), so the warn-and-fall-back path at `main.rs:189` — the one that would trip "clumpify never engaged" — is unreachable here. **[arithmetic-verified]** from `src/clump.rs`.

**I8 — Make rep comparability across the rebase an explicit precondition, and record the runner image. [code-verified]**
Assumption 3 covers cell order only. Add: (a) `git diff --name-only 1b2579e origin/dev -- src/ examples/ Cargo.toml Cargo.lock build.rs tests/` must be **empty** before pooling reps 1–2 with 3–4 — it is today, so this is a check to record rather than a problem to fix, but pooling two binaries would be silent; (b) echo the runner image version (`/imagegeneration/imagedata.json`, or `$ImageVersion`) in the job and note it per rep in the TSV, so Assumption 6's image-bump signature is checkable *within* the four reps rather than only across weeks.

**I9 — Replace "full sensitivity" with the number. [arithmetic-verified]**
Lever A's weakest form is **7.8** points (`phase3/PLAN.md:432`), not 8. Against a +5.0 trigger that is 2.8 points of margin ≈ 1.7σ on the 25 bp cell → ~5 % single-cell miss, and a ~0.15 %/run false-red rate on that cell (≈26 % over 200 runs). Both belong in B2's rationale and in remaining-risk (a), because they turn "mean is preferable" into "mean is necessary — max-centring spends 1.7 of the 2.8 points". Add a line to the TSV telling a future reader what to do on a lone 5–7 point red (re-run once; if it clears, it was spread).

**I10 — State which statistic each macOS anchor is, and derive 51 bp/c2 from the 51 bp offset. [arithmetic-verified]**
1.009 is a max of 12 reps; 1.127–1.164 is the range of 7. Subtracting them understates the read-length effect the 51 bp control exists to measure. And 0.84–0.88 is the union of the 25 bp and 51 bp offsets, not the "~17" the text states — the consistent figure is 0.840. Label each anchor (max / range / single rep) and use like for like on both ends of the difference.

### Optional

**O1 — `n_bins` is not what distinguishes c2 from c4. [code-verified]** `n_bins = max(16, 4 × bin_cores)` = **16 at both**. `B` is larger at c2 because `denom = 23·n_bins + 64·workers` shrinks from 624 to 496 → `B` = 22.80 MiB vs 18.13 MiB, and the pool term's share rises from 59 % to 74 %. Conclusion right, stated mechanism wrong; fix the parenthetical in A1.

**O2 — A4's timeout advice is incoherent.** `timeout-minutes: 30`; "if the job lands over ~12 min, raise the timeout" would have an implementer raising a limit with 18 minutes of slack. Either drop it or say what 12 min would actually mean (the job becoming a cost worth revisiting).

**O3 — The header is stale.** "`dev` at `4f7b907`; PRs #448, #451, #452 open" — #448 and #451 are merged and `dev` is `242546a`. A reader could act on it. **[code-verified]**

**O4 — "~17 points" at 51 bp uses rep 1 alone** (0.939 − 0.770 = 16.9); the reps 1–2 mean the sentence claims gives 16.3.

**O5 — The TSV lookup has no comment-skip and no duplicate guard. [code-verified]** `awk '$1==c'` scans every line; safe today only because every line starts `#`. A provenance line written without `#` becomes a silent data row, and two matching rows make `want` two lines, whereupon awk quietly numeric-coerces the first. C1's assertion covers both.

**O6 — The 8 GiB disk assertion precedes fixture generation**, so it does not assert free disk *at measurement time*, which is what §5's rule is about. Observed 88 GiB makes this low-risk; a second `df` check immediately before the first `measure` would close it for one line.

**O7 — `git merge-tree --write-tree` is not read-only** — it writes tree objects (no refs, no worktree). The plan calls the simulation read-only. Immaterial, but the §0.1 disjointness proof is a better argument anyway and does not need the simulation.

**O8 — B2's placeholder means round exact .0005 ties upward** (0.9325 → 0.933, 0.7895 → 0.790). Tie handling is unspecified for the real computation; they are placeholders, so this only matters as a note on the recomputation step.

---

## 7. Overall

Land it, with **C1** folded into Item B before it is written, and **I1**–**I4** applied to Item A's text before the rebase. The plan is unusually careful where this arc has been bitten before — the rebase hazard is diagnosed exactly right, the placeholder warning on B2 is the correct reflex, and Assumption 7 makes A3 falsifiable rather than decorative. Its weak spots are all of one kind: **decision rules and validations stated one notch more confidently than the arithmetic behind them supports** — the 1.13 threshold, "full sensitivity", "cannot fail the job", "cannot separate the two causes", and a Validation 5 that fires for the wrong reason. None of those change the shape of the work; all of them change what a future reader would conclude from a number.

The one thing I would not let slide is **C1**. Item B's whole purpose is to stop the band being inert-while-looking-armed, and as drafted it fixes today's instance of that by hand while leaving the mechanism that produced it in place — with a later hand-edit of the same file already scheduled.
