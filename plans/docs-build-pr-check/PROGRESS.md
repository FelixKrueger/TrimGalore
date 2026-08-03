# Progress: Build the docs site in PR CI

**Last updated:** 2026-07-30

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **v2** — review findings incorporated, revision history in §11 |
| Plan Review | ✅ Complete | `PLAN_review_reviewer-{A,B}.md` — 2 Criticals, 6 Importants, read-only |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | **Awaiting implementation trigger** |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Notes

- Found while merging dependabot PR #368 (`astro` 7.1.1 → 7.1.3, `@astrojs/starlight` 0.41.3 → 0.41.4, **`satori` 0.28 → 0.29`**): all eight PR checks green, **none of them built the docs site**. `docs.yml` triggers only on `push` to master/dev, and pushes to `dev` auto-deploy to www.trimgalore.com — so the first validation of any docs change is the live deploy.
- #368 was validated by hand before merging: `npm ci && npm run build` in a scratch worktree gave 31 pages, 30 OG PNGs, none zero-byte, exit 0. This plan automates that step.
- **Hosting decision: `ci.yml`, not a `pull_request` trigger on `docs.yml`.** The latter looks smaller but grants `pages: write` to PR runs, queues PR builds against the `pages` concurrency group (serialising them behind live deploys), and needs an event guard on the deploy job. `ci.yml` already runs on `pull_request` for both trunks with no pages permissions.
- **The OG-image assertion is the point.** satori can fail per-image while Astro still exits 0, so "build succeeded" is insufficient evidence — the same silent-success shape that let #369 survive three releases behind a green parity gate. V2/V3 require the new gate to be *shown* to fail on a broken build and on an empty PNG.

## Related

- Same class of gap as the `-a2` CI case added in #370: a check suite that looks green without exercising the thing at risk.
- Third of three items parked after #364 closed. The other two: the v2.4.0 release (pending clump-only feedback; `CITATION.cff` two releases stale) and `docs/reference/changelog.md` mirror automation.

## History

- 2026-07-30: Plan → ✅ Complete (PLAN.md created)

## v1 → v2: what review changed

Both Criticals came from **Reviewer B**, and both landed on the checks meant to prove the gate works:

1. **V2's deliberate break did not break the build.** v1 proposed a malformed frontmatter key; Starlight is not strict about *unknown* keys, so `npm run build` exits **0** — verified. V2 would have gone green and proved nothing. Replaced with removing a required `title:` (verified exit 1) or a bad sidebar slug.
2. **The assertion could be a silent no-op that the validation table could not detect.** `working-directory: docs` means paths must be `dist/…`, but v1 wrote `docs/dist/…` throughout — and a wrong path under `bash -e` without `pipefail` prints to stderr and exits **0**. Worse, V3/V4 were to be rehearsed from the repo root, a different cwd from the job, so both could pass while the step was inert.

**Reviewer A's Importants**, all verified independently: the OG path form covered **3 of 30** images (27 are nested, and v1's V3 targeted a top-level file so it would pass against the bug); the emptiness check is vacuous when `dist/og` is absent; `-size -1k` **misses a 500-byte file** through block rounding, so `-size -10240c` is needed; and the decisive hosting argument was missing — a path-filtered `pull_request` workflow reports *nothing* on non-docs PRs, which is worse than skipped and settles §2.2.

**Six v1 claims corrected**, tabulated in `PLAN.md` §11 — including two false statements about `ci.yml` (`needs:` does exist; `audit` has no schedule guard) and a CHANGELOG precedent that was **backwards** (#247 *did* get a contributor-facing Infrastructure section).

**Added:** `permissions: {contents: read}`; V3b/V3c/V8; and a fix for `justfile:86,90,94`'s `cd Docs` — three sites, 0 tracked files under `Docs/` versus 154 under `docs/`, broken on Linux. That is the same case mismatch already recorded as an environment gotcha, here a real committed bug in a file that exists for local CI parity.

## Session notes worth carrying

Two memories were written/corrected off the back of this work: `feedback_mv_alias_hazard` (**`cp` is aliased `-i` too** — the prior version recommended `cp` as the safe alternative, which cost a 10-minute timeout) and `feedback_harness_false_negatives` (four false negatives in one session; prove a check can fail before trusting that it passed).
