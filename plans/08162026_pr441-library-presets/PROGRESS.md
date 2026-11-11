# Progress: `--library` presets (#440 / PR #441) and the follow-up fixes

**Last updated:** 2026-08-17

> **This directory runs the pipeline out of order, deliberately.** The code review came *first*,
> because the implementation arrived as an outside contribution (PR #441 by `BenjaminDEMAILLE`) rather
> than from this pipeline. `PLAN.md` therefore plans the **follow-up fixes to #441**, not #441 itself.
>
> Note for a future `/progress update`: the review files carry `_reviewer-{A,B}_<generation>` suffixes,
> not the skill's `CODE_REVIEW_A.md` / `PLAN_REVIEW_A.md` names, so automatic detection will not match
> them. There are **two generations** of code review here: `CODE_REVIEW_reviewer-{A,B}_145432e.md`
> review the contributor's commit, and `CODE_REVIEW_impl_reviewer-{A,B}_9a738a3.md` review ours. The
> Code Review row below refers to the second pair.

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **v2**, revised after dual plan review |
| Plan Review | ✅ Complete | `PLAN_REVIEW_reviewer-A_g1.md` (583 lines), `PLAN_REVIEW_reviewer-B_g1.md` (565 lines) |
| Impl Plan | ❌ Excluded | PLAN.md v2 carries the numbered outline; no separate IMPL.md |
| Implementation | ✅ Complete | `2142577` on `fix/440-library-refusals`; 9 files, +403/−28, 12 new tests, 764 passing. **Unpushed** |
| Code Review | ✅ Complete | `CODE_REVIEW_impl_reviewer-{A,B}_9a738a3.md` — A: READY TO PUSH; B: ready with changes (1 High, fixed) |
| Coverage | ✅ Complete | `COVERAGE.md` — **Verdict: COMPLETE**, 100 items, 0 MISSING, 0 PARTIAL |

## Upstream review of PR #441 (separate from the rows above)

| Item | State |
|---|---|
| `CODE_REVIEW_reviewer-A_145432e.md` | ✅ Complete, 519 lines, 6 mutations |
| `CODE_REVIEW_reviewer-B_145432e.md` | ✅ Complete, 12 mutations. Session died on the final revert; annotated with an editor's note |
| Verdict (both, independent) | Mergeable with changes. No spec violations |
| PR #441 merge | ✅ Merged 2026-08-16 as `872c4ae` (squash, author Benjamin Demaille preserved) |

## Plan review outcome

Both reviewers: **READY WITH CHANGES**, no rework. All 14 line anchors verified. Reviewer A applied the
plan's own B2 patch and measured the blast radius: 751 tests, exactly one failure, the predicted one.

Three Critical items (A, all confirmed by B) and Felix's two decisions are folded into v2 — see its
Revision history table. Felix chose: refusal wording `"does not honour the clip flags"`, and editing
#441's CHANGELOG entry in place rather than adding a second entry.

## History

- 2026-08-16: Plan → ✅ Complete **v2** (revised after dual plan review; 3 Critical + 7 Important folded in)
- 2026-08-16: Plan Review → ✅ Complete (dual, independent; one contradiction resolved on priority)
- 2026-08-16: Plan → ✅ Complete v1 (four scope decisions taken first)
- 2026-08-16: #441 dual code review complete — 2 agreed findings, 5 test seams, 1 contradiction
  resolved against Reviewer B (docs register)

- 2026-08-17: Coverage → ✅ COMPLETE; Code Review → ✅ Complete (dual). B's High (top-level JSON key
  permutation invisible to all 764 tests) reproduced, fixed and re-proved; commit amended to `2142577`
- 2026-08-17: Implementation → ✅ Complete (`9a738a3`, amended `2142577`); 11-row mutation matrix
- 2026-08-16: PR #441 merged as `872c4ae` after rebase + CHANGELOG conflict resolution
