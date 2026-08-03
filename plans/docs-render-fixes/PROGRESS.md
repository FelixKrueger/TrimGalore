# Progress: Two docs-render bugs — KaTeX eating prose, and a broken image path

**Last updated:** 2026-08-01

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` — **revised**, all review findings folded in |
| Plan Review | ✅ Complete | `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` |
| Impl Plan | 📋 Planned | Not used — implemented directly from `PLAN.md` |
| Implementation | ✅ Complete | +31/−3 across 6 files, **uncommitted**. See `PLAN.md` §10 |
| Code Review | ✅ Complete | `CODE_review_reviewer-A.md`, `CODE_review_reviewer-B.md`. Verdict: **content fixes ship, CI step needs amending first** |
| Coverage | ✅ Complete | `COVERAGE.md` — **COMPLETE**, 0 MISSING, 0 PARTIAL, 36 items |

## Stage-5 outcome: the plan was wrong, not the implementation

Coverage says COMPLETE and both code reviewers say don't ship the CI step. Those do not conflict — coverage audits *conformance to the plan*, and the step conforms exactly to §4 step 5. The defect is in the plan.

**§8 split V1 (`class="katex"`) from V2 (`katex-error`) precisely because neither sentinel detects the other's defect. §4 step 5 then carried only V2's pattern into CI.** So the assertion added to prevent recurrence cannot detect the changelog corruption this change exists to fix — that page renders 3 `class="katex"` and **0** `katex-error`. Both reviewers proved it exits 0 on a tree whose only defect is that corruption. It passes today only because the unrelated image bug sits on the same page.

Second, agreed independently: the guard runs in `ci.yml` (gates PRs) but not `docs.yml` (deploys on every push to `dev`, no assertions at all). Both original defects reached production by that path.

## Review verdict: diagnosis right, one instructed step harmful

Both reviewers confirmed the two bugs are real, reproduced them in built output, and confirmed `singleDollarTextMath: false` is the correct fix — verified against the installed `remark-math@6.0.0` and `micromark-extension-math`, with `$$` display math **structurally** immune (options reach `mathText` only, never `mathFlow`).

But both independently found the same three critical defects:

1. **§3.3 / §4 step 3 introduces a new bug.** The `\~` escapes at `benchmarks.md:129` guard **GFM single-tilde strikethrough**, not `$` math. Removing them renders the sentence inside `<del>` with a `**` bold pair collapsed to literal asterisks — and **all eight validation checks pass green on it**. The plan's premise ("now redundant… verify the rendered text is unchanged") is wrong.
2. **§2.2's headline evidence is false.** `grep -c '1000-sample cohort'` returns **1** pre-fix, not 0 — the math node ends before that phrase. V4 therefore cannot fail.
3. **A second math node on the benchmarks page renders as a red `katex-error` span.** `katex-error` carries no `application/x-tex`, so V1/V3/V7 are structurally blind to failed math.

Also agreed: `grep -c` counts lines not occurrences; A2 is exhaustively verified (nothing relies on single-dollar math); A3/A4 verified by inspection and curl; and **§7 Q2's `dev` pin is wrong** — both in-repo precedents (`README.md:3-4`, `astro.config.mjs:38`) pin `master`, and the asset serves 200 on both branches.

## What the revision changed

| Area | Change |
|---|---|
| §3.3 | **Reversed.** The `\~` escapes stay, with a source comment naming what they guard. §4 step 3 no longer removes them |
| §2.2 | Corrected the false grep claim; documented the second math node and its red `katex-error` render |
| §2.5, A6 | New — GFM single-tilde strikethrough as a fourth divergence class, with the flanking mechanism |
| §3.2 | Image now **self-hosted** from `docs/public/images/`, dissolving the branch-pin question rather than answering it. Plus an `alt` attribute |
| §3.4, A7 | New — why `singleTilde: false` is rejected at the config layer, and the residual `$$` hazard |
| §4 step 4 | One line, not two — the math defect on the tracked page is fixed by the config change alone |
| §4 step 5 | New — three permanent `ci.yml` assertions (`katex-error`, `src="docs/`, `<del>`), the only change that prevents recurrence |
| §8 | Rewritten: 8 → 12 checks. Sentinel changed to `class="katex`; V4/V6 replaced (both were no-ops); `grep -F` and occurrence counting; scratch-directory baselines; GitHub-side image check added |
| §6 | A1–A5 promoted to verified with evidence; A6 and A7 added |

## Open decision

**§7 Q1 — ordering.** This plan should land before the sibling automation, at the cost of one throwaway line in the tracked changelog page. That cost halved on review (the math defect there is fixed by the config change alone), which strengthens the recommendation.

## History

- 2026-08-01: `PLAN.md` revised — all review findings folded. 198 → 283 lines
- 2026-08-01: Plan Review → ✅ Complete (dual agents). Verdict: **do not implement as written** — one instructed step introduced a regression invisible to every listed check
- 2026-08-01: Plan → ✅ Complete (PLAN.md created). Split out of `plans/changelog-mirror-automation/` per user decision
