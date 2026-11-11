# Progress: Narrow the #408 `--rename` + uBAM refusal for inert Read 2 clip flags (#450)

**Last updated:** 2026-08-18

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **v2**, revised after dual review |
| Plan Review | ✅ Complete | `PLAN_REVIEW_reviewer-A_g1.md` (358 lines), `PLAN_REVIEW_reviewer-B_g1.md` (409 lines) |
| Impl Plan | ❌ Excluded | Implemented straight from `PLAN.md`'s 18-step outline |
| Implementation | ✅ Complete | **Merged to `dev` as `51fa55e`** via PR #461. 8 files, +449/−88, 796 tests on the rebased base. Notes in `PLAN.md` §12 |
| Code Review | ✅ Complete | `CODE_REVIEW_impl_reviewer-{A,B}_g1.md` — 0 Critical, 0 High; 9 recommendations applied |
| Coverage | ✅ Complete | `COVERAGE.md` — **Verdict: COMPLETE**, 52 items, 0 MISSING, 0 PARTIAL |

## Review outcome

Both reviewers verified all **8** load-bearing claims TRUE and found no defect in the proposed
behaviour. A: 2 Critical, 6 Important, 10 Optional, 1 out-of-scope defect. B: 1 Critical, 5 Important,
4 Optional.

| Finding | Both? | Disposition in v2 |
|---|---|---|
| V4's mutation pointed the wrong way; V4 ≡ V5 | **A only** — B blessed V4-V6 | Rewritten; V6b added. The one contradiction between reviewers, and A was right |
| Docs check was a false-negative grep; 3 surfaces state the gate | Both (Critical) | Steps 16-17, three named edits; concept-level grep pattern recorded |
| Read ownership at 6 sites, not 3 | Both | **Scope widened to the full fold-in** (Q-4) |
| Ladder's emission order pinned by no test | Both | Step 15, committed assertions + an all-four case |
| `--three_prime_clip_R1` tested nowhere | A only | Row 4 + step 13 |
| Mixed-ownership row missing | B only (A's I4 overlaps) | Row 7 — the only row whose failure mode is the original #408 data defect |
| V7's matrix has impossible cells; identical bails compare equal | A only | V7 cell rules |
| `any_effective` docstring overstates | B only | Narrowed to read shape only |
| `user_clips()` is a public bypass | B only | `pub(crate)` + doc-scoped |
| `any_set` is published API | B only | Removal stated as deliberate (step 3) |
| Assumption 6 stale; 8 citations off | Both | Corrected; provenance note now lists which files differ |

## Decisions — all closed

| Ref | Decision |
|-----|----------|
| Approach | Narrow the guard, over suppressing the warning or documenting only |
| Q-1 | `--clip_R1 0` over-refuses too — **leave unfiled, out of scope.** Not to be folded in |
| Q-2 | CHANGELOG: **separate bullet** citing #450 |
| Q-3 | **Filed as #450** |
| Q-4 | **Full six-site fold-in**, after both reviewers rejected the three-site middle |
| Q-5 | X1 (`--rrbs` under-refusal) **filed as #459**, out of scope |
| Q-6 | Fix the `--help` overstatement at `cli.rs:312-313` in this commit |

## Issues filed this session

- **#449** — `-a2`'s warn-and-ignore behaviour undocumented (plan Q-1 of the #447 work)
- **#450** — this plan's target
- **#459** — the `--rrbs` under-refusal, found by Reviewer A

## History

- 2026-08-18: **#461 merged to `dev` (`51fa55e`)**, #450 closed. Rebased onto `bc48bf0` first to pick up
  the new `Clumpify peak RSS` CI job; 10/10 checks pass. Three worktrees cleaned up (3.7 GB reclaimed)
- 2026-08-18: #449 closed — fixed by contributor PR #460 (`9370305`), with one commit pushed to their
  branch correcting a single-cause claim that was false for three of the four inert shapes

- 2026-08-17: Implementation → ✅ Complete, uncommitted in `wt-450`. All 18 steps, 5 deviations (D1–D5),
  mutation matrix V4/V5/V6/V6b each red at its target, V7 differential 27 cells / 1 differing / 0 vacuous
- 2026-08-17: Plan → ✅ v2 written; Plan Review → ✅ Complete; Q-4/Q-5/Q-6 resolved; #459 filed
- 2026-08-17: Plan Review → 🚧 dual reviewers launched, after correcting five `src/main.rs` line
  citations read off the working tree rather than `origin/dev` @ `4f7b907`. Reviewer A died mid-task and
  was recovered via `SendMessage`
- 2026-08-17: Q-1/Q-2/Q-3 resolved by Felix; #450 filed
- 2026-08-17: Plan → ✅ Complete (v1; ordering defect found and fixed during self-review)
- 2026-08-17: Felix chose "narrow the guard" over suppressing the warning or documenting only
