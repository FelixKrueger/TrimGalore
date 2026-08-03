# Progress: Generate the published changelog page from `CHANGELOG.md`

**Last updated:** 2026-08-01

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md`, revised twice — manual review then dual agent review |
| Plan Review | ✅ Complete | `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | **Blocked on implementation trigger.** Also sequenced behind `plans/docs-render-fixes/` — see §3.4 |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## The problem, measured

| File | Lines | Last touched |
|---|---|---|
| `CHANGELOG.md` | 1331 | 2026-07-31, `7201f03` |
| `docs/src/content/docs/reference/changelog.md` | 1141 | **2026-07-25**, `4f86c0f` |

Body diff, root H1 stripped from both sides: **196 lines root-only, 0 lines page-only, no missing `###` sections.** Both reviewers went further — every diff hunk is a pure append, so the page body is a **strict subsequence** of the root body. "Purely behind, never divergent" is exact, not an approximation, which is what makes a one-way overwrite lossless and deleting the tracked copy safe. That fact is A1, and the plan says to re-verify it — with a falsifiable harness — immediately before the deletion.

*(The figure was 199 before manual review; the extra 3 were the root's H1 and its blank lines, which the generator now strips. See §2.1.)*

Currently unpublished: `--clump_only` uBAM in/out, both `--phred64` fixes, the mixed-pair diagnosis (#363), the **`-a2` fix (#369)**, and the whole `#### Infrastructure` section. #369 is the one that matters — a public issue whose fix is invisible on the public changelog.

## Approach

**Generate the page at build time and delete the checked-in copy**, so drift is unrepresentable rather than merely detectable. A ~30-line Node generator in `docs/scripts/`, wired into `build` and `dev` with explicit `&&`.

Both reviewers endorsed the rejection of a custom Astro content loader, and both independently proposed the nearest rival — keep the page tracked plus a `git diff --exit-code` CI drift check — which is now recorded and rejected in §5.1.

## What review changed

Three of the plan's *"Verified"* claims were wrong. The design was not.

| Found by | Change |
|---|---|
| Manual review | The leading H1 would have shipped **two H1s** (Starlight renders `title` as `<h1>`); V1 as written was unsatisfiable alongside correct rendering. New §3.0 |
| Manual review | Banner linked to `master` while the page is generated from `dev`. §3.1, V1b, A6 |
| **Both agents** | **A4 wrong** — `docs/README.md:39` also references the page, in a section that becomes false and names a stale *second transformation*. §4 step 6 now mandatory |
| **Both agents** | **§3.2 wrong on both counts** — braces are not all in code spans, and `CHANGELOG.md:1164` is a raw `<img>` in prose. Restated with a divergence table |
| **Both agents** | **V9 could not fail** — `^### ` is already 49/49 on the stale page. Replaced with `^#### ` and body line counts |
| **Both agents** | Move the one-time checks into the generator as permanent assertions. Called the highest-value change by both. §4 step 2a |
| Reviewer B | **The `editUrl` breakage** — the largest single miss; GitHub's `/edit/` on a deleted path opens the *new-file* editor, inviting a contributor to recreate the copy the plan deleted. §4 step 4a |
| Reviewer B | **KaTeX `$` bug** — a live production defect. Split to `plans/docs-render-fixes/` |
| Reviewer B | `slug:`-vs-`link:` soft dependency **overturned §9 Q1**; one `test -f` line in `ci.yml` now §4 step 7 |
| Reviewer A | Line-1-anchored strip hole (§9 Q5); V5's mis-specified expectation; two leftover `199`s the manual fold had missed |

Reviewers contradicted each other on **how to fix the broken `<img>`**: A left generator-side rewriting on the table, B forecloses it to protect §3.0's one-transformation invariant. **B's position adopted.**

## Sequencing

`plans/docs-render-fixes/PLAN.md` should land **first**. Automation makes the page an unattended artefact, so every future `$` in a changelog entry would mangle with nobody in the loop. Do not automate over a known silent-corruption path. No file overlap between the two.

## History

- 2026-08-01: Plan Review → ✅ Complete (dual agents). `PLAN.md` revised: §3.0/§3.2/§3.4 rewritten, §4 steps 2a/4a/7 added, §5.1 added, A4/A5/A7 corrected, A8 added, V3–V6/V9 revised, V11 added, §9 Q1 overturned, Q4/Q5 recorded
- 2026-08-01: Manual review → `PLAN.md` revised (H1 strip, banner branch retarget)
- 2026-07-31: Plan → ✅ Complete (PLAN.md created)
