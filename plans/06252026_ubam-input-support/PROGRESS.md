# Progress: uBAM (unaligned BAM) input support

**Last updated:** 2026-06-25

**Tracking issue:** [#316](https://github.com/FelixKrueger/TrimGalore/issues/316)

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md v3 (cheap fixes + Spike 1 + Spike 2 all consolidated) |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md, PLAN_REVIEW_B.md (dual independent reviewers) |
| Impl Plan | ✅ Complete | IMPL.md (26 tasks, TDD, full plan coverage checklist verified) |
| Implementation | 🚧 Implementing | **21 of 26 tasks done.** SE + PE uBAM both work end-to-end through production binary. 283 tests pass, fmt+clippy clean. Remaining: integration tests (T21-T23), CI job (T24), reproducibility (T25), docs (T26). See IMPL.md §Implementation log Sessions 1–3. |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## History

- 2026-06-25 (session 3): Tasks 15, 17, 20-PE landed. Paired-interleaved BAM reader with MAX_SLACK=1024 de-interleaver. CLI validation rules for --preserve-tags / --passthrough+BAM / --paired+1-file. Single-file paired-uBAM dispatch in main.rs. PE-uBAM works end-to-end.
- 2026-06-25 (session 2): Tasks 16, 18, 19, 20-SE landed. RecordSource trait through entire pipeline (parallel.rs + trimmer.rs + main.rs + adapter.rs). --preserve-tags CLI flag. Format-aware dispatch. SE uBAM works end-to-end. Milestone committed (eed17d3).
- 2026-06-25 (session 1): Tasks 1–14 landed. Deps, format module, BamReader with sync + threaded readers, +17 unit tests.
- 2026-06-25: Impl Plan → ✅ Complete (IMPL.md v1, 26 TDD tasks, 40 plan items covered)
- 2026-06-25: Plan → v3 (review fixes + spike results consolidated; supersedes intermediate v1.1 / v1.2)
- 2026-06-25: Spike 2 (paired-uBAM ordering) → ✅ Complete (`spikes/SPIKE_paired_ordering.md`). A-C4's OOM premise disproved; bounded-buffer design adopted defensively.
- 2026-06-25: Spike 1 (RecordSource dispatch) → ✅ Complete (`spikes/SPIKE_recordsource.md`). Trait wins (+0.86–1.28% vs concrete, under 2% bar).
- 2026-06-25: Plan Review → ✅ Complete (Reviewer A + Reviewer B, independent, 5 agreed-critical findings)
- 2026-06-25: Plan → ✅ Complete (PLAN.md written, v1)
