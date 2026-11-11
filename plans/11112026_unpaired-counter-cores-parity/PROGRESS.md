# Progress: `r1_unpaired` / `r2_unpaired` parity between the serial and parallel paths

**Last updated:** 2026-11-11

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` revision 2 — issue #544 |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — 1 Critical + 2 blocking, all addressed in rev 2 |
| Impl Plan | ❌ Excluded | Outline carried in `PLAN.md` §5, as for #541 |
| Implementation | ✅ Complete | 3 files, +225/−32, uncommitted on `dev`; notes in `PLAN.md` §12 |
| Code Review | ✅ Complete | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` — 0 Critical, 0 High; 4 changes applied |
| Coverage | ✅ Complete | `COVERAGE.md` — **Verdict: COMPLETE**, all 10 steps, all 11 validations |

## History

- 2026-11-11: Code Review → ✅ Complete (both recommend as-is; post-review pass applied, gates re-run)
- 2026-11-11: Coverage → ✅ Complete (`COVERAGE.md`, verdict COMPLETE; all 5 mutations reproduced independently)
- 2026-11-11: Implementation → ✅ Complete (845 tests, clippy and fmt clean, 5 mutation checks)

- 2026-11-11: Plan → ✅ Complete (revision 2: four fixture shapes, host moved to a non-passthrough sibling, JSON gating claim corrected, Perl grounding and #245 shape added)
- 2026-11-11: Plan Review → ✅ Complete (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`)
- 2026-11-11: Plan → ✅ Complete (`PLAN.md` revision 1 created)
