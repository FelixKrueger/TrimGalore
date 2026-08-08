# Progress: Add report paths to the paired pre-flight (#388)

**Last updated:** 2026-08-07

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2 — validation section rebuilt; uBAM twin in scope |
| Plan Review | ✅ Complete | A + B converged: 4 Criticals each (all four v1 tests broken, A1 false); 3 reproducers verified |
| Impl Plan | ❌ Excluded | §3 is the outline; ~6 lines of source |
| Implementation | ✅ Complete | 546 tests; both §6 controls discriminate; §8 notes |
| Code Review | ✅ Complete | A + B: 0 Critical; shared High → #391 filed; M1–M5 applied |
| Coverage | ✅ Complete | **COMPLETE 20/20, first pass** |

## History

- 2026-08-07: Plan → v2. Both reviewers: every v1 test unable to pass or fail (validate's absolutising R1≠R2 fires first; format detector precedes pre-flight); A1 false — paired primaries carry a positional discriminator reports lack, so the SE finer-than argument doesn't transfer. uBAM twin verified live, brought in scope. Two behaviour changes named for sign-off.

- 2026-08-07: Plan → ✅ Complete
