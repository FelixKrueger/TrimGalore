# Progress: uBAM (unaligned BAM) input support

**Last updated:** 2026-06-25

**Tracking issue:** [#316](https://github.com/FelixKrueger/TrimGalore/issues/316)

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md v3 (cheap fixes + Spike 1 + Spike 2 all consolidated) |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md, PLAN_REVIEW_B.md (dual independent reviewers) |
| Impl Plan | ✅ Complete | IMPL.md (26 tasks, TDD, full plan coverage checklist verified) |
| Implementation | 🚧 Implementing | **18 of 26 tasks done.** SE uBAM works end-to-end through production binary. Architectural refactor (RecordSource trait + parallel.rs polymorphic) complete. 281 tests pass, fmt+clippy clean. PE uBAM, integration tests, CI, docs pending. See IMPL.md §Implementation log Session 1 + 2. |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## History

- 2026-06-25 (session 2): Tasks 16, 18, 19, 20-SE landed. RecordSource trait through entire pipeline (parallel.rs + trimmer.rs + main.rs). --preserve-tags CLI flag. Format-aware dispatch in main.rs. PE-BAM rejected for v1 (defers to T15). End-to-end SE uBAM works through production binary.
- 2026-06-25 (session 1): Tasks 1–14 landed. Deps, format module, BamReader with sync + threaded readers, +17 unit tests.
- 2026-06-25: Impl Plan → ✅ Complete (IMPL.md v1, 26 TDD tasks, 40 plan items covered)
- 2026-06-25: Plan → v3 (review fixes + spike results consolidated; supersedes intermediate v1.1 / v1.2)
- 2026-06-25: Spike 2 (paired-uBAM ordering) → ✅ Complete (`spikes/SPIKE_paired_ordering.md`). A-C4's OOM premise disproved; bounded-buffer design adopted defensively.
- 2026-06-25: Spike 1 (RecordSource dispatch) → ✅ Complete (`spikes/SPIKE_recordsource.md`). Trait wins (+0.86–1.28% vs concrete, under 2% bar).
- 2026-06-25: Plan Review → ✅ Complete (Reviewer A + Reviewer B, independent, 5 agreed-critical findings)
- 2026-06-25: Plan → ✅ Complete (PLAN.md written, v1)
