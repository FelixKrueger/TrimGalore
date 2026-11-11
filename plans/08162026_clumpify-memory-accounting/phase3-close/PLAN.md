# Plan — close out Phase 3 of #439: land #452 correctly, then make the drift band real — v3

**Parent:** `../PLAN.md` (#439) · **Phase 3 plan:** `../phase3/PLAN.md` v3
**Issue:** [#439](https://github.com/FelixKrueger/TrimGalore/issues/439), still open · **Live defect:** [#457](https://github.com/FelixKrueger/TrimGalore/issues/457)
**Baseline:** `dev` at **`bc48bf0`** — all three #439 PRs merged (#448 `efc1b11`, #451 `242546a`, #452 `bc48bf0`). Item A complete; **Item B is what remains**
**Reviews:** `PLAN_REVIEW_reviewer-A_phase3-close.md`, `PLAN_REVIEW_reviewer-B_phase3-close.md` — both "land it after fixes"; neither ran a measurement
**Not in scope:** #457's σ refit (Phase 4-sized, deliberately unstarted)

## Revision history

- **v1 — reviewed, two Criticals, both agreed on substance.** Its rebase diagnosis was right but
  its *proof* used the wrong comparison; its Validation 5 fired the wrong clause; and it left the
  band with no permanent self-check that it is armed. Four numeric claims were wrong or narrower
  than the corpus supports: the escalation threshold (1.13 → **1.101**), "≥ 3× the widest spread"
  (**2.94×**, false), the macOS anchor ranges (narrowed by dropping a reviewer's lowest rep in both
  cells), and the `n_bins` mechanism.
- **v2 — reviewed and implemented for Item A.** Folded both reviews in, and added one finding
  neither reviewer reached: **gating ≠ banded**, so B's proposed fix would have left two of the four
  calibrated cells unprotected.
- **v3 — this document.** Item A shipped (#452 merged as `bc48bf0`); reps 3 and 4 landed and changed
  Item B's premise. `1G_c4_fastqc` turned out to be the **least** stable cell, not the most, so it
  cannot be banded at ±0.05; B2's placeholders are replaced by real four-rep means; and the banded
  set is now chosen by a stated rule rather than "every cell with a row". B2 also retracts the
  interaction claim reported after rep 3. **Item B's mechanism is unchanged from v2** — the `BANDED`
  list, the `test -f` preamble, the sentinel assertion and Validations 5–6 all survive review; only
  the list's membership narrows, in the conservative direction, on data that arrived after the
  review. No second dual review.

## Goal

Two items, two PRs, separated by a merge boundary because the second needs measurements the first
produces.

- **Item A** — land #452 on a correct base, and give the `clumpify-memory` job the `--cores 2`
  cells that decide how urgent #457 is. Rides #452's mandatory rebase.
- **Item B** — retire the baseline file's `BOOTSTRAP` mode so the drift half of D5 actually gates,
  **and give it a self-check that it stays armed**.

Neither touches `src/`. No new dependency, no new fixture, no coefficient moves.

## Context

### Placement

`.github/workflows/ci.yml` — the `clumpify-memory` job, lines 271–434 on
`ci/clumpify-memory-guard-439`; the `measure` calls are lines 425–431, `exit $fail` at 434.
`.github/clumpify-memory-baseline.tsv` — 31 lines, all comments, no data rows.

### The rebase is mandatory — and the proof is disjointness, not the file-difference argument

#452's branch carries `1b2579e` — #451's **first** commit — but not `44dc732`, #451's tip.
Squash-merging #451 therefore left #452 `CONFLICTING`: predicted by `git merge-tree`, then
confirmed by GitHub (`mergeable: false`, `mergeable_state: dirty`). The conflict is `CHANGELOG.md`
alone; `src/clump.rs` auto-merges. **Resolving it toward #452 reverts `44dc732`**, the commit
pointing the short-read caveat at #457 — D7's deliverable in weakened form.

v1 argued the rebase is clean because "the only file differing between `1b2579e` and #451's tip is
`CHANGELOG.md`". **Both reviewers independently flagged that as the wrong comparison** — the rebase
lands on `origin/dev`, which also carries #448. The correct argument is set disjointness:

```
git diff --name-only 1b2579e origin/dev        →  CHANGELOG.md
                                                  docs/…/modes/clump-only.md
                                                  docs/…/performance/clumpy.md

git diff --name-only 1b2579e origin/ci/…-439   →  .github/clumpify-memory-baseline.tsv
                                                  .github/workflows/ci.yml
                                                  examples/mk_clump_fixture.rs
                                                  src/clump.rs
```

The intersection is **empty**, so for every file the four replayed patches touch, `origin/dev`'s
blob is byte-identical to `1b2579e`'s. Each patch applies against identical context: conflict-free
is *guaranteed*, not merely likely. Shorter than v1's argument and actually a proof.

No silent-drop risk either: `git rebase` drops commits whose patch-id already exists upstream, and
none of the four patches touches a file `1b2579e..origin/dev` modified, so none can match.

### What is already measured — reps 1 and 2

`rep1 = 69d8f67` (run 32065004600), `rep2 = 3f1a930` (run 32065624366). The **only** difference
between those commits is the TSV's comment block, so binary, workflow steps and fixture generator
are byte-identical: a genuine repeat measure ~7 minutes apart. Both runs: free disk 88 GiB, probe
OK, every cell `exit=0`.

| cell | rep 1 | rep 2 | mean | half-range |
|---|---:|---:|---:|---:|
| `1G_c4` | 0.770 | 0.782 | 0.776 | 0.0060 |
| `1G_c4_fastqc` | 0.697 | 0.701 | 0.699 | 0.0020 |
| `25bp_1G_c4` | 0.941 | 0.924 | 0.9325 | 0.0085 |
| `multipair_c4` | 0.792 | 0.787 | 0.790 | 0.0025 |

**v1's claim that `+0.05` is "≥ 3× the widest observed spread" is false — it is 2.94×**
(`0.05 / 0.017`), and v1 asserted it twice, including as the criterion a future editor should
re-check. A future editor applying a false criterion is worse than one applying none. Replaced with
the two sensitivity figures, which are both stronger and actually the ones a reader needs — excursions
expressed in multiples of each cell's observed half-range, against Phase 3 §5.3's **minimum**
lever-A delta of **+0.078**:

| cell | a false red needs | missing a minimum lever A needs |
|---|---:|---:|
| `1G_c4` | 8.3× | 4.7× |
| `1G_c4_fastqc` | 25× | 14× |
| `25bp_1G_c4` | **5.9×** | **3.3×** |
| `multipair_c4` | 20× | 11× |

Both failure directions need ≥ 3.3× the largest excursion yet observed, on the worst cell. Reviewer
A reaches the same conclusion by a different route — treating each 2-rep range as a σ estimator
(`E[range₂] = 1.128σ`) gives a **~0.15 %/run false-red rate on the 25 bp cell**, ~26 % over 200
runs. The two models are complementary, not contradictory: B's answers "how far beyond anything
seen?", A's answers "how often will this bite?". Both belong in the file.

### The band is inert for **four** independent reasons, not two

v1 named two. Reviewer B found two more, and I reproduced both locally:

1. `# BOOTSTRAP` present → `bootstrap=1`, band skipped (and a `::warning::` says so).
2. No data row → `want` empty → `[ -n "$want" ]` short-circuits, **silently**.
3. **The baseline file missing or renamed.** `grep -q '^# BOOTSTRAP' "$base"` sits inside an `if`,
   so `set -e` is exempt; `bootstrap` stays `0`; `awk` has `2>/dev/null || true`; every `want` is
   empty; **the job exits 0**. Verified locally. The only trace is one `grep: No such file` line.
4. **A label typo, or a cell renamed in `ci.yml` without touching the TSV.** Same silent skip — and
   A1 is itself a precedent for editing the cell list.

Deleting `# BOOTSTRAP` removes the warning as well as the flag, so after Item B every one of paths
2–4 is **silent and green while the log positively implies the band is live**. That is the exact
failure this plan opens by describing, and v1 fixed only today's instance of it, by hand, while
scheduling a further hand-edit of the same file. Hence **C1** in Item B.

## Behaviour

### Item A — the `--cores 2` cells (rides the rebase into #452)

**A1. Append two non-gating cells** after `multipair_c4` (line 431), before `exit $fail`:

```bash
measure 1G_c2        0 --clumpify --paired --cores 2 --memory 1G "$F51R1" "$F51R2"
measure 25bp_1G_c2   0 --clumpify --paired --cores 2 --memory 1G "$F25R1" "$F25R2"
```

**Label `1G_c2`, not `51bp_1G_c2`** (reviewer B). The four existing cells treat unprefixed as 51 bp
and mark only the short-read case (`1G_c4`, `1G_c4_fastqc`, `multipair_c4` vs `25bp_1G_c4`), and
Phase 3 §5.3 calls this cell `1G/c2`. Labels become committed baseline keys, so this must be settled
now: a rename after B writes rows silently unbands the cell via inertness path 4.

Three properties, each with a reason:

- **Appended last, not inserted beside `25bp_1G_c4`.** Cell order sets the disk and page-cache state
  each later cell starts from. Both reviewers code-verified that appending cannot perturb cells 1–4
  (sequential `measure` calls, per-process `ru_maxrss`, no shared mutable state), so reps 1–2 stay
  poolable with reps 3–4 — which is what B2 needs.
- **Both have the ratio check off (`gating=0`).** §5.3 already concluded `--cores 2` is too tight to
  gate on — macOS 51 bp/c2 measured **1.009× prediction**, the max of 12 reps, 1 of which was over.
  Note "non-gating" is a narrow statement: `$gating` guards exactly one line. See A5.
- **`1G_c2` is the control, and what it supplies is the *core increment*.** v1 claimed a single
  25 bp/c2 number "cannot separate the two candidate causes". Both reviewers called that an
  overstatement, and B pinned down why: the read-length effect at fixed cores is **already measured
  on both platforms** at c4 (Linux 0.776 vs 0.9325 = +15.7 points; macOS ~+13). What is genuinely
  unmeasured on Linux is the **c4 → c2 core increment**, against a macOS reference of **+0.070** that
  is length-independent (+0.070 at 51 bp, +0.077 at 25 bp). `1G_c2` measures exactly that term, at a
  read length where macOS shows no budget breach, so it is a clean control. Without it, one 25 bp/c2
  number confounds "Linux's core increment differs" with "Linux's read-length effect differs at c2",
  and cannot support the Phase 4 scoping decision it is being run for. Secondary payoff: `1G_c2` is
  Phase 3's gating-candidate-that-was-too-tight, so the Linux number also says whether it could ever
  be gated.

No new fixture: `f51` and `f25` are both generated already.

**A2. Record rep 2 in the TSV comment block**, with the commit each rep came from (see A6).

**A3. The falsifiable prediction.** v1 predicted an interval and called the expected finding
"breaches the prediction but not the budget" — but its own interval straddled 1.0, so that was true
of only half its range. Reviewer B's sharper form, which one number can refute, replaces it.

Correcting v1's inputs first: it quoted macOS 25 bp/c2 as 1.127–1.164 and 25 bp/c4 as 1.069–1.075,
**both narrowed by dropping reviewer A's lowest rep**. The corpus (`../phase3/PLAN.md` §0.4) gives
**1.119–1.164** at c2 and **1.054–1.075** at c4, so the platform offset is **11–15 points, not
13–15**. Each anchor is also a different statistic and must be labelled: macOS 51 bp/c2 `1.009` is a
**max of 12**; 25 bp/c2 is a **range of 6**.

Two independent routes, agreeing to within 0.01 — which is worth stating because it makes the answer
robust to the modelling choice v1 never named:

| route | `25bp_1G_c2` prediction | peak/budget |
|---|---:|---:|
| multiplicative (offset applied to macOS) | 0.97–1.05 | 0.88–0.96 |
| additive (Linux c4 + macOS core increment) | 1.002–1.018 | 0.910–0.925 |

**The claim to test: a Linux budget breach at 25 bp/c2 requires a core increment of +0.169 — 2.4×
macOS's +0.070.** `1G_c2` measures the increment directly, so the two cells together refute or
confirm this in one run.

**A4. Three thresholds, all named** — v1 gave one and attached the wrong consequence to it. The
budget-breach line is `1024/930 = **1.1011**`, a property of the margin; v1's 1.13 answers the
different question "did the platform offset apply at all?", leaving `[1.101, 1.130)` as a real Linux
budget breach its own rule would wave through.

| peak/pred | meaning |
|---:|---|
| ≥ **1.000** | the printed promise breaks. Phase 3's Validation 15 already calls this a finding, not a threshold to tune |
| ≥ **1.101** | **`--memory` is exceeded** → #457 escalates to "reproduces on both platforms" |
| ≥ **1.130** | macOS's magnitude reproduced → the platform offset does not apply at c2 at all |

**A5. What `gating=0` does and does not buy.** It guards exactly one line
(`test "$peak_kb" -le "$pred_kb"`). Five paths still set `fail=1`: non-zero exit, missing
`^clumpify:` banner, output not `*.fq.gz`, an empty parse field, and the `pred == budget × 10/11`
margin pin. Both reviewers verified the two that could plausibly bite at c2, and both are safe:

- **Margin pin:** `expect_kb = 953,250` vs `pred_kb = 952,320`, diff **−930 KiB**, inside the
  ±2048 KiB slack. It passes at c2 because `denom = 496` is divisible by 16, exactly as `624` is.
- **Plain-mode fallback**, which would trip the banner check: the c2 layout floor is **~281 MiB**,
  below c4's ~289 and far below 1024. Read length does not enter the floor, so the 25 bp cell is
  equally safe.

So v1's risk (c) — "they detect nothing" — was wrong in the plan's own disfavour: these cells newly
guard the plain-mode fallback, the gzip promise and the margin pin **at `--cores 2`**. Claim it.

**A6. Rep comparability is a precondition, not an assumption.** v1 verified reps 1–2 measured
identical code and then pooled reps 3–4 without applying the same bar. It holds today —
`git diff --name-only 1b2579e origin/dev -- src/ examples/ Cargo.toml Cargo.lock build.rs tests/` is
**empty** — but that is a property of `dev` *right now*. `wt-450` is sitting on
`fix/450-inert-clip-refusal`; if #450 or anything else touching `src/` merges into `dev` before
#452, reps 3–4 measure a different binary and B2's mean silently pools two of them.

So: record each rep's commit in the TSV comment, and make B's harvest step assert
`git diff --quiet <rep1-sha> <rep-n-sha> -- src Cargo.lock` per rep. A rep that fails is dropped and
another taken — cheap, because reps keep arriving (A7).

**A7. Extra reps are free.** The job runs on every PR targeting `master`/`dev` and every push to
them, so by the time B is written there will be reps beyond four, from unrelated PRs, at no cost.
Prefer `push: dev` runs (a merged binary) and apply A6's filter. n = 6–8 moves Assumption 5 from
"bounds the spread" toward "estimates it" for zero effort.

**A8. Where the rebase runs.** This checkout is shared — `git worktree list` shows it plus `wt-450`
and `wt-450-base`, and review agents were reading it minutes ago. A rebase rewrites `HEAD` and the
working tree, a strict superset of the branch-switch that has already silently invalidated agent
results here. Do it in a dedicated worktree, or after confirming nothing is mid-verification.
(`git worktree remove` needs a sandbox bypass.)

**A9. Timeout.** `timeout-minutes: 30`, unchanged, with ~5× headroom. From the Actions API for run
32065624366: build 57 s, fixtures 26 s, `Measure and gate` **2m19s for four cells**, total 3m59s.
Two c2 cells at roughly double a c4 cell's wall time is **+1.5–2 min**, landing near 5m30s. v1's
"raise the timeout past ~12 min" had no referent against a 30-minute limit; revisit only past ~20 min.

**A10. Item A must leave a durable record.** Its stated purpose is deciding #457's urgency, but the
c2 ratios print only to a job log that expires, and B3 gives the cells no TSV row. Once reps 3–4
land: record both c2 cells in the TSV comment block, **and post the two numbers plus A3's
prediction outcome as a comment on #457**. Otherwise the deliverable ends as an unrecorded reading
of a transient log.

### Item B — exit `BOOTSTRAP`, and keep the band armed (follow-up PR, after #452 merges)

**B1. Prerequisite — satisfied.** Four reps of the original four cells are in hand: reps 1–2 from
#452's PR runs, rep 3 from Item A's force-push, rep 4 from the merge-to-`dev` push. All four pass
A6's code-identity check against `69d8f67`. The two c2 cells have reps 3–4 only.

**B2. Which cells get a data row — a stated rule, not "every cell measured".**

> **Band a cell only if it has ≥ 3 reps AND a false-red multiple ≥ 3×** — i.e. a spurious red would
> need an excursion at least three times the largest yet observed for that cell.

| cell | reps | false red needs | verdict |
|---|---:|---:|---|
| `1G_c4` | 4 | 6.7× | **band** at `0.777` |
| `25bp_1G_c4` | 4 | 3.4× | **band** at `0.923` |
| `multipair_c4` | 4 | 20× | **band** at `0.790` |
| `1G_c4_fastqc` | 4 | **2.1×** | **exclude — spread** |
| `1G_c2` | 2 | 3.7× | exclude — rep count |
| `25bp_1G_c2` | 2 | 20× | exclude — rep count |

```
1G_c4	0.777
25bp_1G_c4	0.923
multipair_c4	0.790
```

**These are real four-rep means, not placeholders.** Half-up to 3 dp, stated because a ±0.001 choice
on a committed threshold is the class of detail this arc keeps losing: `1G_c4` 3.108/4 = 0.7770;
`25bp_1G_c4` 3.691/4 = 0.92275 → 0.923; `multipair_c4` 3.160/4 = 0.7900.

- **Excluding `1G_c4_fastqc` does not leave it unguarded**, which is what makes exclusion cheaper
  than a per-row margin column. It is `gating=1`, so `peak ≤ prediction` still asserts on it every
  run, with 0.745 against 1.000 of headroom. What exclusion forgoes is *drift* detection on that one
  cell — and that signal is redundant, because lever A adds ~72–104 MiB to **every** clumpify cell:
  `multipair_c4` and `25bp_1G_c2` would each need 11.2× their observed half-range to miss it.
- **A per-row margin column was considered and rejected for now.** It costs a TSV format change plus
  job changes, and — the deciding objection — it invites sizing a width from n = 4, the same
  biased-low estimate that just failed at n = 2. When a cell's spread is unknown the answer is more
  reps, not a guessed width. Hold the per-row column for a cell that needs its own width for a
  structural reason rather than a noise one.
- **`25bp_1G_c4` is banded despite a 1.9× miss margin.** The redundancy above covers the miss
  direction, its 3.4× false-red rate is acceptable, and it is the cell #457 concerns, so drift
  coverage there has independent value.
- **The c2 cells stay out on rep count** even though `25bp_1G_c2` looks tightest of all six. Its
  0.5-point range on two draws is precisely the trap `1G_c4_fastqc` sprang: two draws cannot land far
  apart as often as they can land close, so dispersion from n = 2 is biased low, not merely
  imprecise. Two more pushes and it qualifies.
- **Value = mean, not max.** The upper bound is the load-bearing half, and lever A's weakest form is
  **+0.078** against a +5.0 trigger. Max-centring would move each centre up by roughly its observed
  spread and spend much of that margin.
- **Each row gets a comment giving the observed range, the rep count and each rep's commit** (A6),
  plus the instruction for a lone 5–7 point red: **re-run once; if it clears, it was spread.**
- **Record why the three excluded cells are excluded, in the file.** An absent row is otherwise
  indistinguishable from a typo — which is exactly what B4 exists to catch, so the two must agree.

**B3. Delete the `# BOOTSTRAP` block** (TSV lines 14–17), then **assert `! grep -q '^# BOOTSTRAP'`**
— the sentinel is anchored, so a retitled provenance block beginning `# BOOTSTRAP history …` would
keep the band inert. Keep the observations block as provenance, retitled as history.

**B4. C1 — a permanent self-check that the band is armed.** The highest-value change in either
review, and the one that stops Item B swapping a printed warning for a silent failure.

**One finding neither reviewer reached: `gating` ≠ `banded`.** Reviewer B's fix couples the
assertion to `gating=1`, but the two axes are independent, and after B2's rule they barely overlap:
`1G_c4_fastqc` is gating and **not** banded, while `25bp_1G_c4` and `multipair_c4` are banded and
not gating. Coupling them would leave both banded non-gating cells unprotected — including the cell
#457 concerns. So the banded set is listed explicitly:

```bash
# Cells with a committed baseline. A label here without exactly one row means the
# drift band is not armed for it.
BANDED="1G_c4 25bp_1G_c4 multipair_c4"
```

- In the step preamble: `test -f "$base" || { echo "::error::baseline file missing"; exit 1; }` —
  closes inertness path 3.
- After the `measure` calls, when `bootstrap = 0`: for each label in `$BANDED`, assert the TSV holds
  **exactly one** matching row (`awk '$1==c {n++} END{exit !(n==1)}'`), else `::error::` and
  `fail=1`. Closes paths 2 and 4, and the duplicate-row case where `awk` would quietly coerce the
  first of two.

**B5. B's own run is rep 5, judged against a centre built from four.** It must be green on
unmodified code — but v1's "if it is not, the band width or the rep count is wrong" overstates it:
with a ~0.15 %/run false-red rate, one red could be a tail draw. Correct instruction: **re-run once;
if it clears, record it as spread and revisit the width only if it recurs.**

## Implementation outline

**PR A — into #452, on the rebase** (after #451 merges — done — and per A8, not in a tree something
else is reading):

1. `git fetch origin`, then **assert `git rev-parse origin/dev` = `242546a`** (or whatever the
   intended base is). Nothing downstream catches a rebase onto a stale base.
2. `git rebase --onto origin/dev 1b2579e ci/clumpify-memory-guard-439`.
3. Confirm the replay: `git log --oneline origin/dev..HEAD` = exactly the four `ci(clump):` commits,
   no `feat(clump):`.
4. Confirm the branch no longer carries #451's content: `git diff --name-only origin/dev...HEAD`
   lists only the four expected files, and **`git diff --quiet origin/dev HEAD -- CHANGELOG.md`**
   (the branch must not touch it at all).
5. Append A1's two `measure` lines plus a one-line comment saying the ratio check is off because c2
   is too tight to gate. Repo comment rule: one line, current state, no history, no justification.
6. Update the TSV comment block per A2 and A6.
7. Extract the step's `run` script and `bash -n` it; validate the YAML with
   `Docs/node_modules/js-yaml`. Both have caught real breakage in this job.
8. `git push --force-with-lease=ci/clumpify-memory-guard-439:3f1a930` — **pin the OID**. The bare
   form takes its expected value from the remote-tracking ref, which step 1's fetch has just
   refreshed, defeating the protection it exists for.

**PR B — after #452 merges:**

9. Harvest reps 3–4 (and any extras per A7) via
   `gh api /repos/{owner}/{repo}/actions/jobs/{id}/logs` — `gh run view --log` refuses while any job
   in the run is in progress. Apply A6's code-identity filter to each rep.
10. Recompute the means (B2's warning, half-up), write the rows with their range/commit comments,
    delete the `BOOTSTRAP` block and assert the sentinel is gone.
11. Add B4's `BANDED` list, the `test -f` preamble check and the exactly-one-row assertion.
12. Open the PR; require its own `Clumpify peak RSS` green, per B5.
13. Post A10's comment on #457.

## Efficiency

Item A: **+1.5–2 min** per run (measured basis in A9), no new fixture, no new build, no new
dependency. Item B: text plus ~6 lines of shell, zero runtime cost. Neither touches `src/`, so no
test-suite and no Perl byte-identity surface.

Disk is unchanged by appending — `measure` does `rm -f "$out"/*.fq.gz` per cell, so peak disk stays
`fixtures + one cell's output`. Two residuals, both minor: the 8 GiB assertion runs *before* fixture
generation, so it does not assert free disk **at measurement time**, which is what §5's rule is
about — a second `df` immediately before the first `measure` closes it for one line; and per-cell
`*_trimming_report.{txt,json}` are not cleaned but are negligible.

## Integration

- **Writes:** `.github/workflows/ci.yml` (two `measure` lines + a comment in A; `BANDED`, the
  `test -f` and the row assertion in B); `.github/clumpify-memory-baseline.tsv` (comments in A, rows
  in B).
- **Reads:** the job reads the TSV once per cell. Nothing else in the repo reads either file.
- **Downstream:** once B lands, every PR to `master`/`dev` and every push to them is band-gated on
  four cells. That is the first time a `--clumpify` memory regression fails CI for someone who does
  not know this job exists, so the failure text carries the load — the two existing messages name
  the likely cause, and B2 adds the re-run instruction.
- **Ordering:** A strictly before B. B first would calibrate centres against a cell sequence about
  to change.
- **Note:** with `bootstrap=1` in both measured runs, **the band code has never executed in CI**.
  Item B is its first execution, which is why B4 and Validation 5 matter more than their line count
  suggests.

## Assumptions

1. **Fixed:** `--clumpify` requires `--cores >= 2` (`src/cli.rs:742-751`), so c2 is the lowest
   measurable cell.
2. **Fixed:** reps 1–2 measured identical code — only diff `69d8f67`→`3f1a930` is the TSV comment.
3. **Fixed:** appending cannot perturb cells 1–4 — sequential `measure`, per-process `ru_maxrss`, no
   shared mutable state. Both reviewers code-verified.
4. **Fixed:** `f51`/`f25` already generated at 10 M pairs, so c2 needs no new generation.
5. **Fixed:** `--clumpify --cores 2 --memory 1G` clears the layout floor (~281 MiB) and engages, so
   the plain-mode fallback that would trip the banner check is unreachable. c2's floor is *below*
   c4's, so a cell passing at c4 cannot fail this way at c2.
6. **Open:** that four reps is enough for a centre. It is the top of the file's own bar, and the
   sensitivity table quantifies it, but n = 4 **bounds** the spread rather than estimating it. A7
   improves this for free.
7. **Open:** that the runner image is stable across reps. `ubuntu-latest` is a rolling label and reps
   3–4 are days after 1–2, so this bites *within* the four reps, not only across weeks. Echo the
   image version (`$ImageVersion` / `/imagegeneration/imagedata.json`) and note it per rep, so the
   "all four cells move together" signature is checkable.
8. **Open:** the four cells' run-to-run deviations are correlated to an unknown degree, so the
   redundancy that makes 2.8 points adequate — lever A moves every cell, and the job reds if any one
   fires — is currently accidental rather than measured.
9. **Not assumed:** that #457 reproduces on Linux at c2. A3 records a falsifiable prediction, and A4
   names the line that would settle it.

## Validation

| # | Verify | How | Expected |
|---|---|---|---|
| 1 | The rebase dropped the duplicate | `git log --oneline origin/dev..HEAD` | 4 commits, none `feat(clump):` |
| 2 | **The rebase did not touch `CHANGELOG.md`** | on the rebased branch: `git diff --quiet origin/dev HEAD -- CHANGELOG.md`, and `git diff --name-only origin/dev...HEAD` must not list it | quiet; absent. v1 asserted on `origin/dev`, which passes no matter how badly the rebase went — keep the `grep -c 457` form as a *post-merge* check |
| 3 | Both new cells engage and measure | job log | a `clumpify:` banner and a `peak/pred` line for `1G_c2` and `25bp_1G_c2`; no `::error::` |
| 4 | A c2 ratio above 1.0 does not fail the job | the ratio check is off; the other five paths must still pass | conclusion `success`. Note this is narrower than "cannot fail" — see A5 |
| 5 | **The band fires in both directions, with the right message** | after B, perturb one row ±0.20 and match the error text. The comparison is strict `>`, so ±0.05 would not fire — verified | `−0.20` → `is >0.05 above baseline`; `+0.20` → `is >0.10 below baseline`. v1 tested only `+0.20`, which fires the *lower* clause — the half it called load-bearing was never exercised |
| 6 | **The band is armed and stays armed** | delete a row; typo a label; rename the TSV; add a duplicate row | each fails loudly with a message naming the cell. This is the one that survives the next edit |
| 6a | **The excluded cells stay silent, not accidental** | confirm `1G_c4_fastqc`, `1G_c2`, `25bp_1G_c2` are absent from `BANDED` **and** from the data rows, and that the file says why | no error for them, and the omission is documented. B4 makes an undocumented omission loud; this checks the documented one is deliberate |
| 7 | The activating run is green unmodified | B's own PR run | green; on a lone red, re-run once before concluding (B5) |
| 8 | **The band still catches lever A** | Phase 3's Validation 12, both forms — edit the helper, and re-inline the literal at `parallel.rs:117`/`:943` — via a **draft PR against `dev`**: there is no `workflow_dispatch`, so a pushed branch alone does not run the job | red both times. Arithmetic says it should be: baseline 0.776, trigger 0.826, lever A lands at 0.854/0.888 → red by 2.8/6.2 points. Predicting that is what makes a red-for-the-wrong-reason visible. **Outstanding since Phase 3** |
| 9 | Disk and exit codes | the job's existing clauses, plus a `df` immediately before the first `measure` | 8 GiB assertion passed; every cell `exit=0` |
| 10 | The sentinel is really gone | `! grep -q '^# BOOTSTRAP' "$base"` | true |

Validation 8 decides whether D5 bought what it cost and cannot run until B lands. Validation 6 is
the one that makes Item B durable rather than a one-time hand-check.

## Questions or ambiguities

- **Resolved — two c2 cells, keep both.** Both reviewers agree; the justification is rewritten
  around the **core increment** (A1), which is the term genuinely unmeasured on Linux.
- **Resolved — label `1G_c2`**, matching the existing convention and Phase 3's prose, settled before
  any baseline row exists.
- **Resolved — `25bp_1G_c2` and `1G_c2` wait for a later pass** on rep count, per B2's rule. Two
  more pushes and both qualify.
- **Resolved — `1G_c4_fastqc` is excluded rather than given a per-row margin**, because exclusion
  leaves it promise-gated and the drift signal is redundant across cells (B2). The per-row margin
  column stays in reserve.
- **Resolved — the RRBS `--memory 4G+` savings row** stays dated, and its re-measurement is batched
  with #457's refit: both need real libraries on oxy, so it is one session rather than two.
- **Resolved — Validation 8 runs immediately after Item B**, as a throwaway draft PR against `dev`
  (no `workflow_dispatch`), and is **#439's close condition**. It is the one item that decides
  whether D5's committed baseline file was worth breaking the "never restate a number the binary
  owns" principle for, and its Phase 3 ancestor was never run.
- **Not filed — `1G_c4_fastqc`'s 4.8-point variance.** ~48 MiB of run-to-run spread, well inside the
  budget, so no promise is at risk and nothing is user-visible: an instrument property, not a defect.
  It belongs as a note in the TSV explaining why the cell is unbanded.
- **Open, deferred to post-activation — widen the upper bound to +0.06?** Would cut the false-red
  rate ~5× at ~1 point of detection margin. Needs Assumption 8's correlation measured first.
- **Excluded — #457's σ refit.** Phase 4-sized, deliberately unstarted. Item A measures the axis; it
  touches no coefficient.

## Self-Review

**What v1 got wrong, and the pattern.** Both reviewers independently reached the same two structural
findings — Validation 5 firing the wrong clause, and no permanent check that the band is armed — and
the same four numeric corrections. The pattern in every one is the same: **a decision rule or
validation stated one notch more confidently than the arithmetic behind it supported.** 1.13 came
from measurement history rather than from the arithmetic of the thing being asserted; "≥ 3×" was
never divided out; "cannot fail the job" and "they detect nothing" were both asserted without
reading the five `fail=1` paths; the macOS ranges were narrowed by silently dropping a reviewer's
reps. None changed the shape of the work; all changed what a reader would conclude from a number.

**Logic.** A-before-B has two real dependencies, both stated. The rebase argument is now a proof
rather than an incomplete comparison.

**Efficiency.** Costs are measured rather than estimated (A9).

**Edge cases.** Named: four inertness paths; `gating` ≠ `banded`; the strict-`>` band boundary; a c2
ratio over 1.0; the margin pin and plain-mode floor at c2; a rounding tie on a committed threshold;
a rebase in a shared tree; a binary change between reps; and the runner image rolling mid-series.

**Remaining risks.**

- **(a)** Four reps bounds the spread rather than estimating it, and the derived ~0.15 %/run
  false-red rate on the 25 bp cell means one red over the job's life is likely (~26 % over 200 runs).
  The per-row range comment plus B2's re-run instruction are the mitigation; A7 shrinks the risk for
  free.
- **(b)** Validation 8 has never been run, so D5's band remains asserted-to-work. Its arithmetic now
  predicts a pass, which is a weaker claim than a measurement and is labelled as such.
- **(c)** Assumption 8's correlation is unmeasured, so the four-cell redundancy the margin leans on
  is not yet quantified.
- **(d)** Item A rides a force-push to a reviewed PR, so #452's merged history differs from what was
  reviewed: two `measure` lines and comments. Anything larger belongs in its own PR.
- **(e)** Both reviewers ran zero `trim_galore` invocations, so every figure here rests on the
  corpus and on reps 1–2 from CI. No measurement in this plan is independently reproduced.

## Implementation notes — Item A (2026-08-18)

Branch `ci/clumpify-memory-guard-439` rebased onto `dev` at `242546a`; **PR #452 is now
`MERGEABLE`** (was `CONFLICTING`). Commit `7d676b7` sits on the four replayed `ci(clump):` commits
(`060526a`, `7261bef`, `43ec2d2`, `17afa7a`). 2 files, 28 insertions, 9 deletions. Item B is
deliberately not in this PR.

**Built:** A1 (two cells, labelled `1G_c2` / `25bp_1G_c2` per reviewer B's O2), A2 and A6 (rep 2
plus the commit behind each rep in the TSV comment block).

**A8 was discharged via its second branch rather than a worktree.** A peer session was busy in the
repo, so I established *where*: the main tree was clean with no `target/.cargo-lock`, and the peer
is on `fix/450-inert-clip-refusal` in `wt-450` — a different branch in a different worktree. The
rebase in the main tree was then safe, and the ref-mismatch awkwardness of a detached-HEAD worktree
was avoided.

**The disjointness proof held exactly.** `comm -12` over the two file lists is empty, and the rebase
replayed 4/4 with no conflict — no `CHANGELOG.md` prompt at any point.

**Validations run:** 1 ✓ (four commits, zero `feat(clump):`); 2 ✓ in all three forms
(`git diff --quiet origin/dev HEAD -- CHANGELOG.md` quiet, the branch changes only the four expected
files, and the `#457` caveat plus the "~50 bp and longer" qualifier are present on the branch);
4 ✓ empirically; YAML parses and all four extracted `run` steps pass `bash -n`.

**Validation 4 was measured, not asserted.** I extracted the job's own clauses and drove them with
`peak = 1050 MiB` at c2 — ratio 1.129, `peak/budget` 1.025, a genuine budget breach. The margin pin
passed at exactly the predicted −930 KiB diff, the ratio check was off, the band skipped, and the
job exited **0**. So an appended c2 cell cannot turn the PR red on the ratio axis, which is what
made appending safe rather than merely intended.

**Iteration log:**

- `#1` `git fetch` failed with `nc: connection failed, SOCKS error 2`. The remote is ssh, which the
  sandbox blocks, while `gh` works because it is https to an allowlisted host. The first
  `origin/dev` assertion had therefore been made against a cached ref that could have been stale;
  re-fetched with the sandbox bypass and re-asserted `242546a` before rebasing. Step 1 of the
  outline exists for precisely this, and it fired on its first real use.

### Item A's measurement — rep 3, run 32109372844 (job passed)

`7d676b7`. Free disk 88 GiB, all six cells `exit=0`, `BOOTSTRAP` warning printed as designed.
A6 checked: `git diff --name-only 69d8f67 7d676b7 -- src/ examples/ Cargo.toml Cargo.lock build.rs`
is **empty**, so rep 3 pools with reps 1–2.

| cell | rep 1 | rep 2 | rep 3 | peak/budget (rep 3) |
|---|---:|---:|---:|---:|
| `1G_c4` | 0.770 | 0.782 | 0.785 | 0.713 |
| `1G_c4_fastqc` | 0.697 | 0.701 | 0.707 | 0.642 |
| `25bp_1G_c4` | 0.941 | 0.924 | 0.914 | 0.831 |
| `multipair_c4` | 0.792 | 0.787 | 0.789 | 0.717 |
| **`1G_c2`** | — | — | **0.867** | 0.788 |
| **`25bp_1G_c2`** | — | — | **1.031** | **0.937** |

**A3's prediction holds, on both cells and both routes.** `1G_c2` 0.867 sits inside the predicted
0.84–0.88 (budget 0.788 inside 0.76–0.80); `25bp_1G_c2` 1.031 sits inside 0.97–1.05 (budget 0.937
inside 0.88–0.96). The additive route's narrower 1.002–1.018 is the one that under-shot.

**The answer Item A was run for: #457 breaks the printed promise on Linux but does not exceed
`--memory`.** `25bp_1G_c2` = 1.031 of prediction (> 1.000, so A4's first bar is crossed) but 0.937
of budget, against a breach line of 1.101. So #457 is a **macOS budget breach plus a cross-platform
prediction breach** — the user-facing `--memory` guarantee holds on Linux; what fails there is the
number the banner prints.

**The interaction the control cell was added for is real, and it is larger on Linux.** All four
figures from the same run:

| | c4 | c2 | increment |
|---|---:|---:|---:|
| 51 bp | 0.785 | 0.867 | **+0.082** (macOS +0.070) |
| 25 bp | 0.914 | 1.031 | **+0.117** (macOS +0.077) |

macOS's core increment is length-independent (+0.070 vs +0.077, a 0.007 spread); Linux's is not
(+0.082 vs +0.117, a **+0.035** spread — 5× macOS's). So the refit's calibration matrix needs
**both** platforms: a fit calibrated on macOS would carry over the wrong assumption that the
per-worker term is length-neutral. That is the Phase 4 scoping answer `1G_c2` was added to produce,
and one cell could not have produced it.

**Non-gating was load-bearing, not cautious.** `25bp_1G_c2` at 1.031 would have failed the promise
assertion. Had §5.3's "c2 is too tight to gate" been overridden, this PR would be red now on
correct code.

**Rep 3 erodes Item B's sensitivity margins, and this is the honest update.** Three of four spreads
widened once a third draw landed — exactly what "n = 2 bounds rather than estimates" predicts.
Recomputed on three reps:

| cell | half-range (2 reps → 3 reps) | false red needs | missing a minimum lever A needs |
|---|---|---:|---:|
| `1G_c4` | 0.0060 → 0.0075 | 6.7× (was 8.3×) | 3.7× (was 4.7×) |
| `1G_c4_fastqc` | 0.0020 → 0.0050 | 10× (was 25×) | 5.6× (was 14×) |
| `25bp_1G_c4` | 0.0085 → **0.0135** | **3.7×** (was 5.9×) | **2.1×** (was 3.3×) |
| `multipair_c4` | 0.0025 → 0.0025 | 20× | 11× |

The band still works everywhere, but the 25 bp cell's miss margin is now ~2.1× its observed
half-range rather than 3.3×. That sharpens two questions Item B had deferred: the per-row margin
column, and whether `25bp_1G_c4` should be banded at all in the first pass. **Do not treat B2's
placeholders as converging** — the widest cell moved 2.7 points across three reps.

**Deferred as planned:** A10's #457 comment waits for rep 4 (the merge-to-`dev` push), so the
public note carries more than one rep of the new cells. Rep 3 is not written into the TSV either —
that belongs with Item B's data rows, keeping #452 at the two `measure` lines plus comments that
risk (d) bounds it to.

### Rep 4 — run 32110691552, the merge-to-`dev` push (`bc48bf0`, job passed)

A6 checked: `git diff --name-only 69d8f67 origin/dev -- src/ examples/ Cargo.toml Cargo.lock
build.rs` is **empty**, so rep 4 pools. (The 4d test deletion was in the *first* of the four CI
commits, so reps 1–3 already carried it — no `cfg(test)` reasoning needed.)

| cell | rep1 | rep2 | rep3 | rep4 | mean | range | false red needs | missing min lever A needs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1G_c4` | 0.770 | 0.782 | 0.785 | 0.771 | 0.7770 | 0.015 | 6.7× | 3.7× |
| `1G_c4_fastqc` | 0.697 | 0.701 | 0.707 | **0.745** | 0.7125 | **0.048** | **2.1×** | **1.2×** |
| `25bp_1G_c4` | 0.941 | 0.924 | 0.914 | 0.912 | 0.9227 | 0.029 | 3.4× | 1.9× |
| `multipair_c4` | 0.792 | 0.787 | 0.789 | 0.792 | 0.7900 | 0.005 | 20× | 11.2× |
| `1G_c2` | — | — | 0.867 | 0.894 | 0.8805 | 0.027 | 3.7× | 2.1× |
| `25bp_1G_c2` | — | — | 1.031 | 1.026 | 1.0285 | 0.005 | 20× | 11.2× |

**Two conclusions reported after rep 3 are corrected here.**

**1. The read-length interaction does not survive.** Rep 3 gave a `+0.035` interaction against
macOS's `+0.007` and I reported Linux's core increment as read-length dependent. Rep 4 gives
**−0.009** — opposite sign:

| | 51 bp increment | 25 bp increment | interaction |
|---|---:|---:|---:|
| rep 3 | +0.082 | +0.117 | **+0.035** |
| rep 4 | +0.123 | +0.114 | **−0.009** |
| mean | +0.1025 | +0.1155 | +0.0130 |

Against a per-rep spread of 0.044 in the interaction itself, `+0.013` is not distinguishable from
macOS's `+0.007`. **Withdrawn.** The driver was rep 3's low `1G_c2` reading (0.867 vs rep 4's
0.894, a 2.7-point swing on the control cell) — a single-run decomposition resting on the noisiest
term in it.

What may survive is weaker and not yet assertable: the increment's *magnitude* looks larger on
Linux (mean +0.10 vs macOS +0.07) at **both** lengths, but with 51 bp readings of +0.082 and +0.123
that is a 4-point spread on a 10-point quantity. The refit still wants both platforms — for the
increment's size, not for a length interaction — and the c2 cells need more reps before either
claim is made.

**2. `1G_c4_fastqc` is the least stable cell, not the most — and Item B cannot band it at ±0.05.**
It read 0.697 / 0.701 / 0.707 and then **0.745**, 3.8 points above its previous max. Its range is
now 4.8 points, the **widest of all six**, having looked like the tightest on two reps
(half-range 0.002, and v2's most reassuring row at "25× / 14×").

Concretely: mean 0.7125, so a +0.05 band triggers at 0.7625 against an observed max of 0.745 —
**1.75 points of headroom against a 4.8-point observed range.** Banding it as specified would red on
correct code, which is §5.4's most expensive failure.

Mechanism is plausible rather than mysterious: FastQC's cost scales with `--cores` and mimalloc
retains the pool's pages across the trim→FastQC boundary — the reason Phase 2 made FastQC an addend
rather than a `max`. That boundary is timing-sensitive in a way the pure-trim cells are not.

**So Item B needs a design change before it is written**, and this is the shape both reviewers
anticipated while naming the wrong cell (A's "per-row margin column", B's "give the 25 bp cell its
own wider band" — the data says the cell that needs it is the **FastQC** one):

| cell | verdict on banding at ±0.05 |
|---|---|
| `multipair_c4`, `25bp_1G_c2` | safe — 20× / 11× |
| `1G_c4` | safe — 6.7× / 3.7× |
| `25bp_1G_c4` | marginal — 3.4× / 1.9× |
| `1G_c4_fastqc` | **do not band** at ±0.05: needs a per-row margin or first-pass exclusion |

**What is solid, and it is the decision-relevant part.** `25bp_1G_c2` is the *most* stable cell of
the six (1.031 / 1.026, half-range 0.0025), so #457's answer is two-rep confirmed: the printed
promise is broken on Linux (both reps > 1.000) and `--memory` is **not** exceeded (0.932–0.937
against a 1.101 breach line). A10's #457 comment can now be written on two reps as planned.

## Implementation notes — Item B (2026-08-18)

Branch `ci/clumpify-baseline-band-439` off `dev` at `bc48bf0`. **PR
[#462](https://github.com/FelixKrueger/TrimGalore/pull/462)**, commit `59fe05b`. 2 files, 68
insertions, 27 deletions.

**Built:** B2's three banded cells at the four-rep means, B3's `BOOTSTRAP` deletion, and B4's three
guards — `BANDED="1G_c4 25bp_1G_c4 multipair_c4"`, the `test -f` file-exists check, and the
exactly-one-row assertion.

**Every validation was run against the committed clauses, not a paraphrase.** The band logic had
never executed in CI before this PR, so this is its first exercise:

| # | result |
|---|---|
| 5 | band fires **both** directions with distinguishable messages. `+0.050` passes and `+0.051` fires (strict `>`, as reviewer B predicted); lever A's weakest form (+0.078) reds all three banded cells; a large drop fires the "below baseline" clause |
| 6 | all five armed-check failure modes caught: deleted row, typo'd label, renamed file, duplicate row, and a cell renamed in `ci.yml` only |
| 6a | the three excluded cells are absent from `BANDED`, absent from the data rows, and documented in the file |
| 10 | `! grep -q '^# BOOTSTRAP'` — sentinel gone |
| — | YAML parses; all four extracted `run` steps pass `bash -n` |

**Deviation from the outline, with its reason.** B4 placed the armed-check *after* the `measure`
calls (reviewer A's sketch). It is implemented **before** them, so a broken baseline appears at the
top of the log rather than after six minutes of measurement — but with `fail=1` rather than `exit 1`,
so the ratios still print for whoever has to repair the file. A hard exit there would break the
legitimate workflow of adding a cell and reading its ratio to seed a baseline.

**Rep 5 — the activating run (32113210868) — passed with the band live**, and no `BOOTSTRAP`
warning, so enforcement is real rather than announced:

| cell | rep 5 | band | margin to upper |
|---|---:|---|---:|
| `1G_c4` | 0.785 | [0.677, 0.827] | 0.042 |
| `25bp_1G_c4` | 0.929 | [0.823, 0.973] | 0.044 |
| `multipair_c4` | 0.808 | [0.690, 0.840] | 0.032 |

**Rep 5 moved a number the first draft had already published, and it was caught before merge.**
`multipair_c4` came in at 0.808, 1.6 points above its previous max, widening its range 0.005 → 0.021
and its miss multiple from **11.2× to 2.7×**. The first commit message, PR body and TSV all named
`multipair_c4` as the cell carrying the miss direction "at 11.2×". All three were corrected: the
carrier is `1G_c4` at 3.7×, and the file now names *the widest margin currently measured* rather than
a fixed cell, so the argument cannot go stale the same way again. This is the fourth cell in this arc
whose spread widened as reps arrived.

**Two consequences recorded rather than acted on:**

- **The rule needed a clause rep 5 exposed.** `1G_c2` (3.7×) and `25bp_1G_c2` (5.9×) now clear the
  3× multiple on three reps — but their third rep *is* rep 5, this PR's own run, so banding them
  would let a cell certify itself. The rule now reads "≥ 3 reps **that predate the PR adding its
  row**". `25bp_1G_c2` also tripled its range on that third draw, which is the signal to wait.
- **Baselines stay at the four-rep means.** Folding rep 5 in would move `1G_c4` 0.777 → 0.779 and
  cost B5's independent test for 0.002. Recorded in the file as excluded-by-design, to be folded in
  at the next deliberate re-baseline.

**Iteration log:**

- `#1` None. The band, the guards and all five armed-check cases worked first time, which is
  attributable to having run the extracted clauses locally before pushing rather than to the code
  being obvious.

### A6 judged, not proxied — `dev` moved under #462

While #462 was in CI, `dev` gained `9370305` and `51fa55e` from the `--library`/`--rename` session
(#460, #461). Those touch `src/cli.rs`, `src/library.rs`, `src/main.rs`, `src/report.rs`, so **A6's
proxy (`git diff --quiet … -- src Cargo.lock`) would flag reps 6+ as non-poolable.**

Judged at source instead: the only clump-related hunk in that diff is a **doc comment** mentioning
`--clump_only` — no `clump`, `channel_depth`, `sync_channel`, `resolved_threads`, `bin_byte`,
`LayoutInputs` or `predicted_peak` code changes. The clump memory path is untouched, so reps 6+
remain comparable to 1–5 in substance.

This is the case the TSV's closing line now anticipates: *"when `src/` moves, judge whether the
change can reach the clump path before pooling."* The proxy is deliberately conservative — it exists
to force the judgement, not to make it. Left unexamined it would have retired a usable rep series on
a docs-comment edit.

#462 needs no rebase: it is `.github/`-only and `dev`'s new commits are `src/`-only, so the two are
disjoint and GitHub reports `MERGEABLE`.

### Note for #457's σ refit — `src/clump_only.rs`'s banners move in #464

PR #464 (issue #456, another session) adds the predicted peak to all four `--clump_only` banners, so
each of those four format strings gains two arguments. It uses `memory_budget_bytes` and
`predicted_peak_bytes()`, both already in scope, so **no model work and no peak-RSS change** — and
`--clumpify`'s banner is byte-unchanged, with
`integration_preflight_tripwire.rs`'s `predicted peak ≈ 930 of 1024 MiB budget` assertion untouched.

Two consequences for the refit, if #464 merges first:

- The four `clump_only.rs` banner sites the refit will touch have moved and now carry the prediction.
- The CI job's parse (`predicted peak ≈ (\d+) of (\d+) MiB budget`, `head -1`) gains a second
  possible source of banners. No current cell runs `--clump_only`, so nothing changes today — but the
  first-match convention is now load-bearing against two modes rather than one, and a future
  `--clump_only` cell would become parseable rather than needing new plumbing.

No impact on Item B's bands: every banded cell is a `--clumpify` run.

## Validation 8 — discharged, both forms (2026-08-18)

Outstanding since Phase 3, where `BOOTSTRAP` made it untestable. Run on a throwaway draft PR
(#467, closed; branch deleted) against `dev` at `45876c5`, since the workflow has no
`workflow_dispatch`.

| | `1G_c4` | `25bp_1G_c4` | `multipair_c4` | unit suite | band |
|---|---:|---:|---:|---|---|
| baseline | 0.777 | 0.923 | 0.790 | — | — |
| **form (b)** literal re-inlined at `parallel.rs:117`/`:945` | 0.903 | 1.051 | 0.941 | **503/503 pass** | **red ×3** |
| **form (a)** `channel_depths` edited | 0.922 | 1.114 | 0.917 | 4c sum pin fails | **red ×3** |

**The result that justifies D5 retrospectively.** In *both* runs the two `gating=1` cells stayed
**under** the printed prediction — 0.903 / 0.814 and 0.922 / 0.809 against 1.000 — so the absolute
promise assertion **passed on every cell it guards while the regression was live**. Combined with
form (b)'s 503/503, a repo gating only on `peak <= prediction` plus its unit suite would have shipped
this regression silently. §5.3 argued exactly that four rounds ago and could never test it; it is now
measured, and D5's committed baseline file has earned the principle it broke.

**The plan under-estimated the delta, in the safe direction.** §5.3 put lever A at 7.8–11.2 points
and v2 used +0.078 as the conservative bound for every sensitivity figure. Measured: **+12.6 to
+19.1** points on the banded cells. So detection margin against the +0.05 trigger is ~7.6 points
rather than the 2.8 v2 worried about, and `25bp_1G_c4`'s "1.9× miss margin" is really ~4.4×. Every
sensitivity number in v2 and in the TSV is conservative by ~1.6–2.4×; they are left as committed,
since a tighter estimate would only relax a guard that is working.

**Two incidental observations, neither shipped behaviour:**

- Under the regression, `25bp_1G_c2` reaches 1.115–1.117 of prediction and **1.013–1.015 of budget**
  — i.e. lever A would push Linux over `--memory` at 25 bp / `--cores 2`. Relevant as context for
  #457's severity: the shipped code has ~8 points of headroom there, not none.
- `1G_c4_fastqc` moved +11 points too (0.701 → 0.809/0.814), so had it been banded it would also have
  fired. Excluding it cost nothing in detection on the regression it exists to catch, which is the
  redundancy argument B2 made, confirmed.

**#439's close condition is met.** Every deliverable is merged (#448, #451, #452, #462) and the guard
is proven rather than asserted. Not closed here — that is Felix's call, and #453 is still open with a
behaviour decision outstanding.

## The band is mis-tuned — 27 observations, not 4 (2026-08-18)

**A false red landed on #466** (another session's PR, run 32122024334): `multipair_c4` at 0.852
against baseline 0.790, +6.2 points, firing alone while `1G_c4` sat at +1.2 and `25bp_1G_c4` at
−1.9. The re-run passed 10/10 with no code change. The TSV's diagnostic worked exactly as written —
one cell moving alone is spread, lever A moves every cell — but it cost a real re-run on someone
else's PR.

That prompted collecting **every** `clumpify-memory` observation to date (31 job logs, Validation 8's
two deliberately-regressed runs excluded), rather than reasoning from the 4–5 reps the baselines were
set from:

| cell | n | mean | min | max | range | false-red multiple | banded |
|---|---:|---:|---:|---:|---:|---:|---|
| `1G_c4` | 27 | 0.7883 | 0.770 | 0.825 | 0.055 | **1.8×** | yes (0.777) |
| `25bp_1G_c4` | 27 | 0.9138 | 0.898 | 0.941 | 0.043 | **2.3×** | yes (0.923) |
| `multipair_c4` | 27 | 0.7987 | 0.771 | 0.852 | 0.081 | **1.2×** | yes (0.790) |
| `1G_c4_fastqc` | 27 | 0.7069 | 0.688 | 0.768 | 0.080 | 1.2× | no |
| `1G_c2` | 25 | 0.8789 | 0.860 | 0.894 | 0.034 | 2.9× | no |
| `25bp_1G_c2` | 25 | 1.0301 | 1.007 | 1.050 | 0.043 | 2.3× | no |

**Not one cell clears the 3× bar B2's own rule requires**, including all three that were banded on
it. `1G_c4` came within 0.002 of its trigger (0.825 against 0.827). And the baselines are visibly
off-centre — 1.1, 0.9 and 0.9 points from their 27-observation means — because they were set from
four draws that happened to cluster.

**This is the same error as every other in the arc, at the largest scale yet.** A dispersion estimate
from n = 4 was used to pick a threshold, and n = 4 estimates dispersion low by construction. It was
caught three times in miniature (`25bp_1G_c4` at rep 3, `1G_c4_fastqc` at rep 4, `multipair_c4` at
rep 5) and each time the response was to note the widening rather than to conclude the *method* was
wrong. The rule "band only if the false-red multiple is ≥ 3×" was itself computed on the biased
estimate, so it certified cells it should have excluded.

**What makes a fix tractable is Validation 8.** It measured lever A at **+0.126 / +0.128 / +0.151**
per cell, where v2's sensitivity analysis had assumed the predicted **+0.078**. There is far more
room between "wider than the noise" and "narrower than the regression" than the plan believed:

| cell | 3× half-range needs | at a +0.09 upper bound |
|---|---:|---|
| `1G_c4` | +0.083 | observed max clears by 5.3 pts; lever A fires by 3.6 — **OK** |
| `25bp_1G_c4` | +0.065 | clears by 6.3; fires by 3.8 — **OK** |
| `multipair_c4` | **+0.122** | clears by 3.7; fires by 6.1 — only 2.2×, **below the bar** |

**Recommended:** re-centre the baselines on the 27-observation means, widen the upper bound
**+0.05 → +0.09**, and **unband `multipair_c4`** — its 8.1-point range is the worst of the six and a
single global width cannot cover it at 3× without approaching the regression it must catch. Two
banded cells remain, both clearing 3× and both measured firing on lever A. Keep the lower bound at
−0.10. Unbanding costs little for the same reason it did for `1G_c4_fastqc`: a depth revert moves
every cell, so either survivor catches it.

**Restate the rule two-sidedly.** As written it constrains only the noise side. It should read: the
band width must be ≥ 3× the observed half-range **and** ≤ the measured regression delta minus a
working margin — and the observed half-range must come from a double-digit sample, not from four
reps.

Not implemented. Needs a decision on the widening, since it relaxes a shipped gate.
