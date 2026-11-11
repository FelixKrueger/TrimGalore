# Progress: Honest `--clumpify` memory accounting (#439)

**Last updated:** 2026-08-17

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v3 |
| Plan Review | ✅ Complete | v1 ×2, v2 ×2 (`PLAN_REVIEW_*`) |
| Spike | ✅ Complete | `spikes/SPIKE_memory-levers.md` |
| Phase 1 (structural) | ✅ Merged | PR #442, squashed to `1595d59` |
| Phase 2 calibration | ✅ Complete | `phase2/CALIBRATION.md` + both TSVs |
| Phase 2 implementation | ✅ Complete | Uncommitted on `dev`; 734 tests pass |
| Code Review (Phase 2) | ✅ Complete | `CODE_REVIEW_{A,B}_phase2.md`; both High findings + all Mediums applied |
| Coverage (Phase 2) | ✅ Complete | `COVERAGE_phase2.md` INCOMPLETE→addressed: 4 undocumented deviations recorded, C5 byte-identity measured |
| Phase 3 | 📝 Planning | `phase3/PLAN.md` **v2** — reviewed ×2, both "implement with changes"; v2 rebuilds both guards |
| Plan Review (Phase 3) | ✅ Complete | `phase3/PLAN_REVIEW_{A,B}.md`; 3 Criticals each, 2 agreed |

## Where #439 stands

The `2 GiB / --cores 4` cell — the worst measured breach — has gone
**1.173 → 1.014 → 0.887** of budget across Phase 1 and Phase 2. Every measured
cell is now inside both the printed prediction and the `--memory` budget, with
and without `--fastqc`.

## Notes

- Calibrated coefficients: σ = 23/16, k = 4 per core (worse of macOS 1.41/3.99 and
  Linux 1.29/3.71), static 224 MiB + 24 MiB/core when FastQC is requested, and a
  10/11 margin sized by the largest observed within-cell spread (96 MiB).
- The FastQC term was **measured after an initial wrong guess**: a flat 288 MiB
  addend derived from a pre-#437 number that turned out to be allocator retention.
  FastQC's real cost scales with `--cores` (~76/139/283 MiB at c4/c8/c16).
- Two harness failures were found and fixed during this phase, both instrument
  rather than code: a fixture of 30,000 distinct reads cycled 227× (replaced by a
  10 M-pair fixture with ~1.05M distinct per 2M), and a Linux measurement where the
  child's `ru_maxrss` was inherited from a bloated forking parent. The harness now
  aborts if its own footprint exceeds 200 MiB.
- **Decided 2026-08-17 (D1):** `--clump_only` below the floor stays a **hard
  failure**, documented. Reordering is the whole job there, so a fallback would
  write a byte-copy of the input and exit 0. Measured regression surface is
  narrower than the earlier note implies: only `--clump_only --fastqc --cores ≥ 11`
  regresses, inside a ≤ 150 MiB band ([620, 775) MiB at `--cores 16`). Without
  `--fastqc` the floor **dropped** at every core count.
- **Decided 2026-08-17 (D2):** the banner names the budget —
  `predicted peak ≈ 930 of 1024 MiB budget`. (It prints 930, not 931.)

## Post-review state

- `cores` is now three values (`LayoutInputs`): bins / workers / FastQC threads.
  Fixes `--fastqc_args -t N` reserving for the wrong count, and `--clump_only`
  paying for workers a synchronous mode never starts.
- FastQC charge capped at 16 threads (floor at `--cores 32 --fastqc`: 1012 MiB
  rather than 1434).
- Output is **no longer byte-identical** to earlier releases under `--clumpify`:
  measured +0.27% / +0.24% size, record multiset unchanged. Stated in the CHANGELOG.
- 734 tests pass; the two previously-inert bail tests now sit one byte under the
  floor.

## Found while reviewing Phase 3 (2026-08-17)

**Every RSS measurement in this arc was taken on the plain-text output branch**, or compared
against numbers that were — Phase 2's calibration, both plan reviewers', and Phase 3 v1's.
On the real gz path the promise holds (`1G/c4` 0.939, `1G/c8` 0.876 of prediction, matching
the calibration). On the plain-output branch `1G/c8` reaches **1.057 × budget**, reproducibly
(3 of 5 reps at 10 M pairs).

It hid because every harness asserts the `clumpify:` banner appeared — which catches a
fallback to plain *mode* — while the banner prints `(gzip level 1)` over uncompressed bytes.

Root cause is a CLI defect, not the memory model: `--gzip` is deprecated-and-ignored with the
message "Output is gzipped by default", which is false for plain input; there is no other way
to request gzip; and `--clumpify` rejects `--dont_gzip` as "pointless" while doing plain
output whenever the input is plain. **Routed to its own issue (D3), not Phase 3.**

## History

- 2026-08-17: Phase 3 plan review ×2 → ✅ Complete; PLAN.md → v2 (both guards rebuilt; D3, D4 decided)
- 2026-08-17: Phase 3 → 📝 Planning (`phase3/PLAN.md`; D1 and D2 decided)
- 2026-08-17: Phase 2 implementation → ✅ Complete (calibrated model, FastQC per-core term)
- 2026-08-17: Phase 2 calibration → ✅ Complete (both platforms, distinct-read fixture)
- 2026-08-16: Phase 1 → ✅ Merged as `1595d59` (PR #442, all 9 checks green)
- 2026-08-16: Code Review + Coverage (Phase 1) → ✅ Complete
- 2026-08-16: Plan v3, spike, plan reviews (v1 ×2, v2 ×2) → ✅ Complete
