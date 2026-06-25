# Progress: uBAM (unaligned BAM) input support

**Last updated:** 2026-06-25

**Tracking issue:** [#316](https://github.com/FelixKrueger/TrimGalore/issues/316)

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md v3 (cheap fixes + Spike 1 + Spike 2 all consolidated) |
| Plan Review | ✅ Complete | PLAN_REVIEW_A.md, PLAN_REVIEW_B.md (dual independent reviewers) |
| Impl Plan | ✅ Complete | IMPL.md (26 tasks, TDD, full plan coverage checklist verified) |
| Implementation | ✅ Complete | **26 of 26 tasks done.** 296 tests pass, fmt + clippy clean, reproducibility intact (md5 verified), end-to-end SE / PE / --preserve-tags all working. See IMPL.md §Implementation log Sessions 1–4. |
| Code Review | ✅ Complete | CODE_REVIEW_A.md, CODE_REVIEW_B.md (dual independent reviewers). Cheap fixes (parse_sam_tag_name unit tests, specialty.rs format-aware) folded in. Remaining medium/low findings documented as follow-up items. |
| Coverage | ✅ Complete | COVERAGE.md (plan-manager Mode B audit). Session-3 verdict was INCOMPLETE (T21–T26 expected gaps + 2 bugs caught); session-4 closed all gaps. Re-run plan-manager for fresh COMPLETE verdict. |

## History

- 2026-06-25 (session 4): Tasks 21–26 + reviewer-flagged fixes landed. Integration tests (8 new tests in tests/integration_ubam.rs), CI validation-ubam job, reproducibility check verified, full docs pass (CLAUDE.md + README.md + CHANGELOG.md). Specialty modes now format-aware (PLAN §3.4 invariant holds). parse_sam_tag_name unit tests added. 296 tests pass total.
- 2026-06-25 (session 3): Tasks 15, 17, 20-PE landed. Paired-interleaved BAM reader with MAX_SLACK=1024 de-interleaver. CLI validation rules for --preserve-tags / --passthrough+BAM / --paired+1-file. Single-file paired-uBAM dispatch in main.rs. PE-uBAM works end-to-end. Committed as e842614.
- 2026-06-25 (session 2): Tasks 16, 18, 19, 20-SE landed. RecordSource trait through entire pipeline. --preserve-tags CLI flag. Format-aware dispatch. SE uBAM works end-to-end. Milestone committed (eed17d3).
- 2026-06-25 (session 1): Tasks 1–14 landed. Deps, format module, BamReader with sync + threaded readers, +17 unit tests.
- 2026-06-25: Impl Plan → ✅ Complete (IMPL.md v1, 26 TDD tasks, 40 plan items covered)
- 2026-06-25: Plan → v3 (review fixes + spike results consolidated; supersedes intermediate v1.1 / v1.2)
- 2026-06-25: Spike 2 (paired-uBAM ordering) → ✅ Complete (`spikes/SPIKE_paired_ordering.md`). A-C4's OOM premise disproved; bounded-buffer design adopted defensively.
- 2026-06-25: Spike 1 (RecordSource dispatch) → ✅ Complete (`spikes/SPIKE_recordsource.md`). Trait wins (+0.86–1.28% vs concrete, under 2% bar).
- 2026-06-25: Plan Review → ✅ Complete (Reviewer A + Reviewer B, independent, 5 agreed-critical findings)
- 2026-06-25: Plan → ✅ Complete (PLAN.md written, v1)
