# Progress: curated GA release notes, and a workflow that cannot ship without them

**Last updated:** 2026-11-18

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **revision 2** — Step 3 and Step 5 rewritten after review (§11) |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — 3 Critical found by both, 0 contradictions |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | — |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Notes

- Policy decision taken 2026-11-11; manual review 2026-11-13; agent review 2026-11-18.
- **The review's headline: revision 1's own remedy would have reintroduced the v2.3.0 partial ship.**
  Step 3's snippet used `$VERSION`/`$SHA`/`$TITLE`, which the real step never defines. `bash -e` (not
  `-u`) expands them to empty, the expansion exits 0, and `gh release create ""` fails *after* the tag
  push. A plan about not stranding tags proposed a remedy that strands tags.
- Revision 1's dry-run proof was wrong three ways, each a different kind of mistake: an assertion that
  could not fail (no dry run pushes a tag, wherever the check sits), an impossible expected outcome
  (job 4 is *skipped* on a dry run, because `smoke-test-docker` is dry-run-gated and the workflow has
  no status functions), and a proof that could not be run from where the work lands (a dispatch runs
  the workflow on the dispatched ref, and `master` has no notes step until the FF).
- The placement never needed a dispatch: all eight other jobs transitively `needs: check-release`, and
  with no status functions a failed job 0 skips them all. That is now 5a, asserted in CI.
- Every line number revision 1 cited was verified exact; the GA/prerelease gate is tight across all
  five reachable dispatch combinations.
- Three new open questions in §12: 5b versus the single-commit convention, adding `actionlint`, and
  `cancel-in-progress` cancelling a run inside the tag-push window.

## History

- 2026-11-18: Plan → ✅ **Revision 2** (Step 3 rewritten, Step 5 replaced by 5a/5b, Step 2a added,
  3 wrong claims corrected, 3 open questions raised)
- 2026-11-18: Plan Review → ✅ Complete (both reports; 3 Critical found independently by both)
- 2026-11-13: Plan → ✅ Manual review complete (§9 resolved, §10 folded in)
- 2026-11-11: Plan → ✅ Complete (`PLAN.md` written)
