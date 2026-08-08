# Progress: `--clump_only --paired` clumping reports join the collision pre-flight (#391)

**Last updated:** 2026-08-08

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md r2 — review findings + Felix's two decisions folded in (shared-constant hint fix, clump_report_candidates helper, 4 prose corrections, 4 new tests incl. Shape A discriminator + A3 property test); awaiting his sign-off |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md ("approve with changes"), PLAN_REVIEW_B.md ("no criticals") — fix confirmed sound; gaps in validation + several plan-text errors to correct; 2 forks for Felix (hint text, helper-vs-inline) |
| Impl Plan | ✅ Complete | Folded into PLAN.md r2 implementation outline (no separate IMPL.md, per pipeline convention) |
| Implementation | ✅ Complete | 2704b52 + review batch d99f23e; expected-fail control 6/6 |
| Code Review | ✅ Complete | CODE_REVIEW_A.md (approve, 3 Medium additive), CODE_REVIEW_B.md (approve, 1 Medium docs) — batch applied in d99f23e |
| Coverage | ✅ Complete | COVERAGE.md: INCOMPLETE (1 assertion) → resolved in d99f23e, addendum recorded |

## History

- 2026-08-08: Implementation + Code Review + Coverage → ✅ Complete (review batch d99f23e; coverage gap closed)
- 2026-08-08: Plan revised to r2 (dual-review findings + Felix's hint/helper decisions incorporated; re-presented)
- 2026-08-08: Plan Review → ✅ Complete (dual reviewers; fix sound, validation gaps + plan-text corrections needed; revision pending)
- 2026-08-08: Plan → ✅ Complete (PLAN.md written; presented for manual review)
