# Progress: Run the `--phred64` pre-flight whenever output is uBAM

**Last updated:** 2026-11-08

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` — plan-review findings in §12, implementation notes §13, fix round §14 |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` |
| Impl Plan | ❌ Excluded | Went straight from design plan to implementation; `PLAN.md` §5 served as the spec |
| Implementation | ✅ Complete | §5.1–5.7 plus the post-review fix round. Uncommitted on `dev` |
| Code Review | ✅ Complete | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` — no functional defect from either |
| Coverage | ✅ Complete | `COVERAGE.md` — **verdict COMPLETE**, 7/7 on §5 and 7/7 on §9 |

## Verification

- `cargo test --release` — **842 passed**, exit 0 (re-run after the fix round).
- `cargo fmt --all -- --check` — exit 0.
- Clippy — cold gate established by the auditor's run (45 compiled / 148 checked, zero warnings);
  re-checked clean for `trim-galore` after the fix round. See §14's caveat: my own "cold" re-run
  was warm, because subagents share the parent scratchpad and the target dir was already populated.
- All three prescribed mutations killed exactly their predicted targets (§13).

## Outstanding, deliberately

- **Prose Lows not fixed** (B4–B7, A's six) — listed in §14. B7 (`CHANGELOG.md:35-40` still says
  the refusal happens "at the FastQC step") predates this change and is stale regardless.
- **Docs site not built** — `npm install` rewrites `docs/package-lock.json`. Link target verified
  by inspection only.

## Spin-off

- **FASTQ-path `-q` band** (§2.4, §10): a mislabelled Phred+33 input with `-q` inside the band
  silently deletes good bases on the FASTQ path — measured 39 bp vs 59 bp at `-q 5`. Out of scope
  here, no good gate available. **Not yet filed.**

## History

- 2026-11-08: Implementation + Code Review + Coverage → ✅ Complete (fix round applied, 842 green)
- 2026-11-08: Plan → ✅ Complete (revised; 2 Critical + 5 Important findings folded in)
- 2026-11-07: Plan Review → ✅ Complete (PLAN_REVIEW_A.md, PLAN_REVIEW_B.md)
- 2026-11-06: Plan → ✅ Complete (PLAN.md created)
