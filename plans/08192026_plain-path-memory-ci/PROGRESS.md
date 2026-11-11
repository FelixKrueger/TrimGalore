# Progress: a peak-RSS baseline for the plain trimming path in CI

**Last updated:** 2026-08-19

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — both verdicts: right shape, revise before implementing |
| Impl Plan | ⏸️ Deferred | Plan needs revision first — see Review outcomes |
| Implementation | 📋 Planned | Two passes: bootstrap (record only), then band |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Review outcomes (2026-08-19)

Both reviewers: the gap is real and the design is the right shape, but revise before implementing.

**Agreed, verified in source:**

1. `plain_c1` fails its own engagement assertion — the banner is guarded by `if cli.cores > 1`
   (`src/main.rs:858`), so a `--cores 1` run prints none.
2. The banner is a flag echo (`Using {cli.cores} worker threads`, printed before any pool exists), so it
   is **not** the `^clumpify:` analogue the plan claims and cannot see a pool that collapsed to serial.
3. The banding rule's regression side is mis-specified — wrong tree (A) and wrong unit (B).

**Scope correction (Reviewer A, verified):** `1595d59`'s own commit message records
`plain -c4 97 -> 91 MiB`, and the 2x2 confirms −6.4 MiB under mimalloc. So on the shipping tree #442 is
a ~6% *improvement*, not an +18 MiB regression. A CI band cannot catch a 6% move against runner noise,
so this job's honest target class is **allocator-scale (~40%+) moves only**. It would have caught #437;
it would **not** have caught #442 on this tree.

**Design addition (Reviewer B):** a between-cell ordering gate — `peak(c1) < peak(c2) < peak(c4)` — is
host- and allocator-independent because all cells run in one job on one runner. Armable on run #1 with
no baseline and no sample, and computed from the measured quantity, so unlike the banner it catches a
collapsed pool. Structural, not drift: complements the absolute band rather than replacing it.

**Open question closed (Reviewer B):** separate job, decisively — the repo is public so runner minutes
are free, and the model job's whole run is 5 min (73 s build, 30 s fixtures, 214 s measure).

**Bootstrap method (Reviewer A):** replace passive accumulation across `dev` with n>=15 draws on one
pinned commit via `workflow_dispatch` + `gh run rerun`. Pooling across `dev` mixes trees with different
true means on a path most commits can reach, inflating the width with the signal it exists to catch.

## Notes

- Motivated by two plain-path memory events that shipped to `dev` unseen this week: `2ab0918` (#437,
  mimalloc, ~1.9x at `--cores 8`) and `1595d59` (#442, a shared-path buffer reservation whose sign
  flips with the allocator). Both surfaced only via nf-core/rnaseq#1903.
- Mirrors the `clumpify-memory` job in `.github/workflows/ci.yml` and reuses its `BOOTSTRAP` mechanism,
  its GNU-time capability probe, and its free-disk assertion.
- **Deliberately a two-PR shape.** The first pass records and does not enforce, because the existing
  baseline file's own rule requires a double-digit sample per cell before a width can be chosen. The
  3-rep figures from `dockyard-oxy-0` cannot set a band for `ubuntu-latest`.
- Measurement evidence for the two known deltas lives in
  `../08182026_threading-docs-cluster-allocation/MEASUREMENT.md`.

## Related

- `plans/08182026_threading-docs-cluster-allocation/` — the docs follow-up from the same investigation;
  its own plan is ⏸️ Deferred pending a rewrite, but the allocator question that blocked it is settled.
- FelixKrueger/TrimGalore#472 — the docs correction PR (65 bp mislabel + allocator cost), open on `dev`.

## History

- 2026-08-19: Impl Plan → ⏸️ Deferred (plan revision required before implementation)
- 2026-08-19: Plan Review → ✅ Complete (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`)
- 2026-08-19: Plan → ✅ Complete (`PLAN.md` created)
