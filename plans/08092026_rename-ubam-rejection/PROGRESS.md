# Progress: reject `--rename` with `--output-format ubam` (#408)

**Last updated:** 2026-08-09

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md **r2** — dual review folded in; scope now format-gated (Felix), which relocates the check to main.rs §3.4b and resolves the --clump_only pre-emption for free |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md + _B.md — both REVISE; Perl-parity premise CONFIRMED by both; A2 found false (C1 guard asserts success); B applied r1's patch and ran the suite |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | Only on exact trigger |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## History

- 2026-08-09: Plan revised to r2 (dual-review findings + format-gated scope decision)
- 2026-08-09: Plan Review → ✅ Complete (dual reviewers)
- 2026-08-09: Plan → ✅ Complete (Perl-parity evidence gathered from the 0.6.11 tag; reject-vs-repair decided; PLAN.md written)
