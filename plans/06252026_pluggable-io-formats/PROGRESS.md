# Progress: Pluggable I/O formats for TrimGalore + Bismark

**Last updated:** 2026-06-26

**Tracking issues:** [TrimGalore#315](https://github.com/FelixKrueger/TrimGalore/issues/315), [Bismark#1025](https://github.com/FelixKrueger/Bismark/issues/1025)

## Phase Status

| # | Phase | Status | Directory | Notes |
|---|-------|--------|-----------|-------|
| 1 | TrimGalore uBAM output (re-scoped, v2.1) | 📝 Planning | `phase1-trimgalore-formats/` | PLAN.md v2.1 written. uBAM input ✅ (TrimGalore#317); uBAM output (opt-in, interleaved paired) pending. Awaiting one more lightweight reviewer pass on v2.1 deltas, then impl trigger. |
| 2 | Bismark uBAM input | 📋 Planned | `phase2-bismark-formats/` | Depends on `RecordSource` from #317 (already stable), NOT on Phase 1's new code. |
| 3 | Cross-tool integration | 📋 Planned | `phase3-cross-tool-integration/` | Depends on #1 and #2. |
| 4 | BINSEQ + mim (deferred from #1) | ⏸️ Deferred | — | 4 objective gating conditions per EPIC §8. Re-evaluate every 6 months from Phase 1 ship. Sunset at 12 months if no progress. |

## History

- 2026-06-26: Phase 1 PLAN → v2.1. Addressed all dual-review (C+D) findings: tag-source plumbing decided (parse from FastqRecord.id tail; §3.6 documents round-trip); §3.4 split into CLI-validate-time and format-detection-time rules; `--demux` rejection table row added; missing-qual semantics decided (raw Phred always); `write_paired_reports` helper extraction added as §5 step 0; `assert_ubam_eq` helper specified for fixture comparison. EPIC.md → v3 (Phase 2 wording corrected; Phase 4 gates made objective with 12-month sunset).
- 2026-06-26: Dual plan-review on Phase 1 PLAN v2 → PARTIAL verdict (2 critical + 5 important; no architectural blockers). Reports at `phase1-trimgalore-formats/PLAN_REVIEW_{C,D}.md`.
- 2026-06-25: Phase 1 PLAN → v2 (re-scope). Dropped BINSEQ + mim + `RecordSink` trait per dual plan-review findings (PLAN_REVIEW_A.md + PLAN_REVIEW_B.md). New scope: opt-in uBAM output via `--output-format ubam`, interleaved paired BAM, no trait abstraction. Epic Phase 4 added for deferred BINSEQ + mim work.
- 2026-06-25: Dual plan-review on Phase 1 PLAN v1 → 5 critical findings, plan v2 required. Reports at `phase1-trimgalore-formats/PLAN_REVIEW_{A,B}.md`.
- 2026-06-25: Epic created. Spawned from cross-repo discussion in TrimGalore#315 + Bismark#1025. uBAM input for TrimGalore (#317) already landed and is the first deliverable of Phase 1.
