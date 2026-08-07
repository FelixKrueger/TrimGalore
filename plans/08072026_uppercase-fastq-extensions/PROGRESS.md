# Progress: Match FASTQ extensions case-insensitively (#384)

**Last updated:** 2026-08-07

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2 — 3 Criticals + 8 Importants folded in |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`; Open 1 answered fold-both by both, independently |
| Impl Plan | ❌ Excluded | `PLAN.md` §5 is the outline; too small for a separate `IMPL.md` |
| Implementation | ✅ Complete | 538 tests, fmt + clippy clean; V6 controls both discriminate; notes in `PLAN.md` §13 |
| Code Review | ✅ Complete | A + B: 0 Critical, 0 High; all Mediums applied |
| Coverage | ✅ Complete | INCOMPLETE→resolved: the one PARTIAL (missing `.bgz` assertion) applied same session |

## History

- 2026-08-07: Plan → v2. A: caller enumeration missed all three `clump_only.rs` sites; v1's V4 collision was already rejected today (couldn't fail); naive helper panics on multi-byte names. B: Perl parity *holds* on compression and the fold breaks it knowingly; the real new refusals are mixed-spelling pairs; V5 was self-referential. Both: fold both on Open 1.

- 2026-08-07: Plan → ✅ Complete. Scope grew during planning: `is_gzipped` is case-sensitive too and decides output *compression*, so folding only `strip_fastq_extensions` would trade a naming mismatch for silently uncompressed output (`PLAN.md` §2.3, A3).
- 2026-08-07: Decision recorded on #384 — fold case rather than pin, because Perl v0.6.11 produced *agreeing* names here and the current mismatch is the rewrite's own.
