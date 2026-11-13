# Progress: Replace `bgzf::MultithreadedWriter` on the parallel uBAM output path (#542)

**Last updated:** 2026-11-13

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` — revised after review (§13) |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` |
| Impl Plan | 📋 Planned | — |
| Implementation | ✅ Complete | 882 tests, fmt + clippy clean |
| Code Review | ✅ Complete | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` — agreed findings fixed (PLAN §15) |
| Coverage | 🚧 Implementing | `COVERAGE.md` — INCOMPLETE; no code gap, 2 edge-case tests left out of scope |

## Notes

- Scope decisions taken before the plan was finalised: thread-model documentation is **out of
  scope** and filed as #568; #546 was closed with a summary comment.
- The plan corrects the design sketch in #542 itself — a scoped `bgzf::Writer` appends an EOF
  block on drop, so the sketch would have put one in every batch (PLAN §3.1).
- Both reviewers upheld all four §3 facts and both found the same two blocking defects, in the
  plan's own validation rather than its design: §9.1 could not fail on the mistake it existed to
  catch, and §9.3 was not constructible. Both rewritten; see PLAN §13.
- Two further decisions taken 2026-11-13: split out `BgzfSink<W: Write>` so the write-error path is
  testable, and accept the throughput regression rather than lowering the compression level.

## History

- 2026-11-13: Coverage → 🚧 (COVERAGE.md INCOMPLETE: no code gap; 2 edge-case tests deliberately out of scope)
- 2026-11-13: Code Review → ✅ Complete (both reports; 6 agreed findings fixed, 1 rejected as false)
- 2026-11-13: Implementation → ✅ Complete (882 tests; 6 mutations verified to kill their guards)
- 2026-11-13: Plan Review → ✅ Complete (both reports written; PLAN.md revised, §13)
- 2026-11-13: Plan → ✅ Complete (`PLAN.md` written)
