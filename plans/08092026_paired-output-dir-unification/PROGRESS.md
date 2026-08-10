# Progress: one output directory per pair (#398)

**Last updated:** 2026-08-10

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md **r2** — dual review folded in; behavioural footprint now bidirectional (a refusal disappears), second inverted guard named, output-vs-input class added, step 7's SE contract fixed |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md + PLAN_REVIEW_B.md — both REVISE, six Criticals, no contradictions. `PLAN_REVIEW_A2.md` is a stopped redundant respawn; ignore it |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | Only on exact trigger |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Notes for the implementer

Three test changes are **predicted** in r2 §Context. Each is a red or newly-green test whose obvious "fix" is the wrong one:

1. **Two inverted guards**, not one — `tests/integration_output_collision.rs:1029-1041` (trim, block 3 only) and `:1206-1240` (clump twin, whole body). Both must be **rewritten to assert the refusal, never narrowed**; narrowing restores the split layout on the candidate side while the writers move.
2. **A refusal disappears** — `clump_paired_rejects_shared_mate_report_without_output_dir` (`:1247-1275`) flips to asserting success, and the untested trim twin of that shape needs a new test.
3. **A new output-vs-input hazard** in #409's exact shape — a mate named like the other mate's prospective report. Needs a test per paired report arm plus `--no_report_file` siblings.

Two mechanical traps: use `assert_dir_holds_only`, **not** `assert_rejected_cleanly` (which asserts an *empty* directory — "fixing" that failure damages eight other tests); and `clump_report_candidates` has **five** call sites, two of which pass the whole single-end input list and must keep per-input parents.

## History

- 2026-08-10: Plan → r2 (dual review folded in; six Criticals. Felix reconfirmed the refusal decision after being shown the mate-side-directory layout, with both layouts to be named in plan, CHANGELOG and docs)
- 2026-08-09: Plan Review → ✅ Complete (A 565 lines, B 680, both REVISE; A wrote non-incrementally so it looked stalled and a redundant A2 respawn was started then stopped)
- 2026-08-09: Plan → ✅ Complete (r1, after two critical scope questions were resolved by Felix; the pre-existing `--passthrough` third-directory rule and the new collision were both found during investigation, neither was in the issue text)
