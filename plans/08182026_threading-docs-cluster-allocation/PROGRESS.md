# Progress: paired-end memory figures and cluster-allocation guidance in the threading docs

**Last updated:** 2026-08-19

## Live couplings to honour when this resumes (recorded 2026-08-19)

Shipped in #472 (`951deab` on `dev`, published to trimgalore.com):

1. **A forward-looking claim is now live in two places** — `threading.md`'s *"The next release uses
   roughly twice this on the plain path"* and the matching CHANGELOG entry. Both are true only while
   `mimalloc` stays the global allocator. **If the allocator decision is revisited before the next
   release, both must be revised**, or the docs ship a false claim. The decision was: keep the speed,
   document the cost.
2. **The sentence both plan reviewers flagged is still live**: *"well below the noise floor of cluster
   scheduler memory allocation"* (`threading.md`, in the cross-tool comparison paragraph). #472 fixed
   the *data* (the 65 bp mislabel) but not the *framing*, and reviewers identified this framing as what
   made a 1 GB budget look generous. It is a one-line change, deliberately left out of #472 to keep that
   PR scoped to factual corrections.
3. **Two known figure inconsistencies remain**, both needing the measurement work this plan gates on:
   `threading.md`'s ~75 MiB at 150 bp/`--cores 4` (a small 400 MB fixture, possibly input-limited) now
   sits beside 73.0 MB at 65 bp/`--cores 4`; and `clumpy.md`'s 104 MiB for 10M x 51 bp paired at
   `--cores 4` (a random-ACGT fixture, so incompressible and an upper bound).

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` — **premise falsified by review; needs rewrite** |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — both reviewers independently rejected the central premise |
| Impl Plan | ⏸️ Deferred | Blocked: see Blockers below |
| Implementation | 📋 Planned | — |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Blockers

1. **The plan's central premise is false.** Both reviewers found, and it was independently verified,
   that the published `threading.md` memory table is **paired-end** (`docs/perf_data/buckberry-2026-04-29/memory_summary.txt`:
   `84M PE reads`; `scripts/benchmark.sh:6,11`: `--paired`). The plan's "the figure was single-end, so
   paired is ~2× it" reasoning has no basis and must be re-derived.
2. **A possible memory regression outranks this plan.** `threading.md:57` publishes 91.6 MB at
   `--cores 8` paired (beta.7, Linux); `src/clump.rs:249-250` records 173 MiB for paired runs on
   Linux/x86_64. Both paired, both Linux, ~1.9× apart, unreconciled. If a reconciliation run against
   the Buckberry pair does not reproduce ~92 MB, that is a code finding and this docs plan should stop.
3. **Version/allocator confound.** `mimalloc` became the global allocator in `2ab0918` (#437,
   2026-08-16), *after* the `v2.3.0` tag, and #437 reports a −24% peak-RSS move on the clumpify path.
   Figures must be attributed to a specific binary and allocator, and the nf-core container pins
   `trim-galore:2.3.0` (stock allocator).

## Notes

- Docs-only follow-up to nf-core/rnaseq#1903 (TRIMGALORE OOM-killed at a 1 GB cap); the config fix
  shipped as nf-core/rnaseq#1905 and is unaffected by the findings above.
- Two public statements need correcting: the #1903 comment and the #1905 PR body both say the cited
  ~100 MB figure was single-end. Pending maintainer sign-off.
- Builds on #454 / #465 (commit `0e86e2a`), which re-qualified the long-read ceiling the same day.

## History

- 2026-08-18: Impl Plan → ⏸️ Deferred (premise falsified; regression question must resolve first)
- 2026-08-18: Plan Review → ✅ Complete (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`)
- 2026-08-18: Plan → ✅ Complete (`PLAN.md` created)
