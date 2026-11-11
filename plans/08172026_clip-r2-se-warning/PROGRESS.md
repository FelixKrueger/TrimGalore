# Progress: Read-2 clip warning on single-end runs, and the refusal's precedence tests

**Last updated:** 2026-08-17

Follow-ups to `plans/08162026_pr441-library-presets/` §13 — Reviewer A's F-1 and F-2 from the
post-implementation review of `2142577`. Base `dev` @ `fdade7e`.

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **v2** — full `-a2` treatment, revised after dual plan review |
| Plan Review | ✅ Complete | `PLAN_REVIEW_reviewer-A_g1.md` (379 lines), `..._B_g1.md` (462 lines) — A: 1 Critical; B: 2 Critical |
| Impl Plan | ❌ Excluded | `PLAN.md` carries the numbered outline; no separate `IMPL.md` |
| Implementation | ✅ Complete | `ddc8230` (amended from `3161b8a`); 8 files, +495/−22, 13 new tests, 781 passing. **Unpushed** |
| Code Review | ✅ Complete | `CODE_REVIEW_impl_reviewer-{A,B}_3161b8a.md` — both READY WITH CHANGES, no behavioural defect; A 1 High + B 1 unique Medium, all fixed bar M-3 |
| Coverage | ✅ Complete | `COVERAGE.md` — **Verdict: COMPLETE**, 65 items, 0 MISSING, 0 PARTIAL |

## Scope

| Item | In scope |
|---|---|
| A's F-1 — warn when `--clip_R2` / `--three_prime_clip_R2` is given on a single-end run | ✅ |
| A's F-2 — three ordering tests pinning the refusal's precedence | ✅ |
| B's F-5 (user-visible half) — the four primitive clip flags are silent no-ops under the specialty modes | ✅ **now in scope** — Felix chose the full `-a2` treatment, so they warn with a mode-accurate reason. Refusing them stays out (plan R-3) |

## History

- 2026-08-17: Plan → ✅ Complete. Key design decision recorded: the warning reads the user's raw
  `clip_r2` fields, never `effective_clips()`, or every single-end preset run would warn

- 2026-08-17: Plan → ✅ **v2**. Felix chose the full `-a2` treatment. v1 shipped the exact state
  `integration_adapter2.rs:414-422` documents as a bug (warn, then announce the displaced value), and
  told a two-file `--clock` run "this is a single-end run". v2 gates the override lines on both
  surfaces and adds the reason ladder; scope now absorbs F-5's user-visible half
- 2026-08-17: Plan Review → ✅ Complete (dual). Both found the juxtaposition unfixed, from different
  surfaces: A the text report, B the stderr stream. Baseline corrected 764 → 768 by both

- 2026-08-17: Implementation → ✅ `3161b8a`. 18/18 V-rows verified, incl. each gate reverted separately
  (V8a/V8b). V18 differential: the only behavioural output change in 24 files is one report line

- 2026-08-17: Code Review + Coverage → ✅. Two test-coverage gaps in the code's own branches (A's H-1,
  M-1) shipped green and are now closed; two prose claims false of the JSON report corrected. B's M-3
  (the #408 guard is a third read-awareness site) left open as a behaviour decision — see §13
