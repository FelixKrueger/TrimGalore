# Progress: Honest `--clumpify` memory accounting (#439)

**Last updated:** 2026-08-18

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
| Phase 3 | ✅ **Complete** | `phase3/PLAN.md` v3 — **#448 docs and #451 guards merged** (`efc1b11`, `242546a`); **#452** CI open and `MERGEABLE` after Item A's rebase (`7d676b7`). Six issues filed (#453–#458) |
| Plan Review (Phase 3) r1 | ✅ Complete | `phase3/PLAN_REVIEW_{A,B}.md`; 3 Criticals each, 2 agreed |
| Plan Review (Phase 3) r2 | ✅ Complete | Two reviewers of phase3 v2 — `PLAN_REVIEW_B_v2.md` (40 runs) and `PLAN_REVIEW_A_v2_measured.md` (**62 verified runs**) — plus `PLAN_REVIEW_A_v2.md`, which reviews the **parent** plan at `2ab0918` and ran zero invocations. Two A-lineages were in flight; the `1050.2 MiB` figure belongs to the parent lineage. |
| Plan Review (Phase 3) r3 | ✅ Complete | `phase3/PLAN_REVIEW_A_phase3_v3.md` — verdict **implement**. One Critical (4d inverted the parent's stated precondition → deletion moved to PR 3) plus two unguarded layers above §4a's pins. |
| Phase 3 close-out plan | ✅ Complete | `phase3-close/PLAN.md` **v3** — Item A: the #452 rebase + two `--cores 2` cells; Item B: the `BOOTSTRAP` exit + a self-check that the band stays armed |
| Item A (rebase + c2 cells) | ✅ **Merged** | #452 squashed to `bc48bf0`; all three #439 PRs now in. Two recorded `--cores 2` cells. Reps 3–4: **#457 breaks the printed promise on Linux (1.026–1.031) but not `--memory` (0.932–0.937 of budget)** |
| Item B (BOOTSTRAP exit) | ✅ **Merged** (#462 → `45876c5`) | Bands **three** cells by a stated rule (≥3 reps AND false-red ≥3×): `1G_c4` 0.777, `25bp_1G_c4` 0.923, `multipair_c4` 0.790. `1G_c4_fastqc` excluded on spread, both c2 cells on rep count. Mechanism unchanged from v2, so no second dual review |
| Validation 8 (band catches lever A) | ✅ **Discharged, both forms** | Draft PR #467 (closed, branch deleted). Band red ×3 in both forms; **both `gating=1` cells stayed under the prediction and form (b) passed 503/503**, so promise+suite alone would have shipped it. Measured delta +12.6 to +19.1 points, above the 7.8–11.2 assumed |
| Phase 4 plan (σ refit) | 📝 **v2**, decision reopened | `phase4/PLAN.md` v2 — dual-reviewed. Both reviewers independently found the same Critical: **the per-record cost is `size_of`, already in `parallel.rs`'s `reserve_exact`** (74 B per record-byte-of-r), not something to fit. v1's two-point fit was 2.6× low, so its cost figures were ~3× too favourable and Option C's price is ~19–20%, not 5–6%. A fourth option (D — account for the spine in the flush threshold, removing `r` from the model) is cheaper at every read length. Awaiting a C-vs-D decision |
| Plan Review (Phase 4) | ✅ Complete | `phase4/PLAN_REVIEW_reviewer-{A,B}_phase4.md` — both "do not treat §6 as decided". Four of v1's numbers corrected; two matrix design faults found (constant pair count biases short-r `s` **downward**; the falsification cell was 3× below noise) |
| Plan Review (close-out) | ✅ Complete | `PLAN_REVIEW_reviewer-{A,B}_phase3-close.md` — both "land it after fixes". Two Criticals, **both agreed**: Validation 5 fired the wrong clause, and the band had no permanent armed-check. Four numeric corrections agreed. Neither reviewer measured. |
| ⚠️ Artifact lost | — | A confirms it destroyed the parent's **16 Aug** `PLAN_REVIEW_A_v2.md` — its brief named that filename, the lease guard pointed at `PLAN_REVIEW_A.md` instead, and it did not `ls` first. Unrecoverable (`plans/` untracked; `_snapshot_20260817/` captured A's own file 2 min after the write). Parent v3's header cites all four reviews, so findings were likely absorbed. **Fix: fresh filename per generation, snapshot before writing.** |

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

## Live defect found 2026-08-17 — short reads breach the budget

**The `--memory` promise is broken at 25 bp / `--cores 2`**: 1048–1083 MiB against a 1024 MiB
budget, 7 reps across two independent harnesses, gz path, ample disk, all `exit=0` with 10 M
records both mates. Narrows with cores (c4 0.97× budget, c8/c16 0.88×) because `B` shrinks and
the bin pool stops dominating. User-reachable via `--small_rna` (miRNA trims to ~22 nt).

Mechanism: σ is treated as a constant but is a **spine** factor — 25 bp packs 1.57× more records
per bin byte than 51 bp, so the under-charge scales with record count, not byte count. Phase 3
v3 estimated this at ~33 MiB and concluded a gated cell "would not fire"; that was ~3× too small
in the reassuring direction.

**Falsifies `CHANGELOG.md:30` as written** (still `### Unreleased`, so nothing shipped).
**D7: scope the claim, ship Phase 3, file the σ refit separately.**

## Retracted (2026-08-17)

An earlier entry here claimed the plain-text output branch was a distinct memory regime
reaching **1.057 × budget**. **That is withdrawn.** A controlled 2×2 (records byte-identical,
reps interleaved, 40 runs) plus independent re-measurement put `1G/c8` plain output at
756–826 MiB against a 930 MiB prediction — no breach in any encoding combination, and output
encoding measures neutral (0.2% at c8, gz *higher* at c4).

Cause was the harness: those readings were taken with ~1.6 GiB free disk while each rep wrote
~2.9 GB, and `io::PendingOutput` writes a temp file beside the previous output. A stalled
writer backs up the result channel and `pending` map, so **disk exhaustion inflates peak RSS**
— and the harness grepped `peak=` without ever checking the exit code, so a failed run's
parseable banner read as data. See [[feedback_harness_false_negatives]] item 7.

What survives is a smaller CLI defect with no memory claim: `--gzip`'s deprecation message
says "Output is gzipped by default", which is false for plain input, and `--clumpify` rejects
`--dont_gzip` as "pointless" while `clump-only.md:78` documents it as supported. Two of the
original three claims were also wrong — naming the input `*.gz` *does* yield gzipped output,
and "always gzip" shipped in v2.1.0-beta.5 and was reverted under #245. Filed narrowed (D6).

## #439 CLOSED 2026-08-18

All deliverables merged (#442, #446, #448, #451, #452, #462) and the guard verified against the
regression it exists to catch, in both forms. Closing comment:
[#issuecomment-5326448016](https://github.com/FelixKrueger/TrimGalore/issues/439#issuecomment-5326448016).

**#453 CLOSED.** Flag decision ([#issuecomment-5326443441](https://github.com/FelixKrueger/TrimGalore/issues/453#issuecomment-5326443441)):
`--clumpify --dont_gzip` stays refused — `--clump_only --dont_gzip` prints `1926073 bytes ->
1926073 bytes`, its own proof that clumping without a compressor buys nothing. The `--clump_only`
asymmetry is deliberate and is the escape hatch for clumped-but-plain output. Both false banner strings **fixed and merged** in #468
(`c7e147d`): `--clumpify` on plain input now prints `(plain output)` plus a warning, and
`--clump_only --dont_gzip` no longer names a level beside the flag disabling it.
`io::output_is_gzipped` and `clump_only::banner_encoding` are the single sources, so message and
behaviour cannot drift apart again. Gzipped banners byte-unchanged, so the CI parse and the tripwire
are untouched; both guards proven able to fail.

**#453 was closed with one string still live**: the `--gzip` deprecation notice ("Output is gzipped
by default", false on plain input) is fixed in **#466 — open, green, mergeable, another session's
PR**. The closing comment says to reopen if #466 is abandoned. An earlier comment of mine on that
issue claimed #466 "has already fixed it"; corrected there.

Still open: #454, #456 (#464 merged, issue not yet closed), #457, #458, and #466 as the last piece
of #453. #455 retired as non-reproducing.

## Correction to this arc's own output — #455 does not reproduce (2026-08-18)

**#455 was filed by this arc and its premise is wrong.** The `--library`/`--rename` session tested it
rather than reading it: a panic injected into `process_single_batch`, `--cores 4` over 20,000 records
→ exit 1, no output published, no `.partial` left, while the control writes all 20,000.

#434 joins the workers **before** `pending_output.commit()`, so every route to a dropped
work-receiver already surfaces an error — a worker `Err` via the result channel, a panic via
`join().is_err()`, or the main thread having already dropped `result_rx` because it is bailing. So
the three `let _ = …send(…)` EOF sites cannot silently truncate output.

The `let _ =` is still real, and is safe **only** through that non-local invariant, which nothing in
the code states. That is a readability concern, not the data-integrity defect we filed.

Worth recording against the filing bar: both round-1 reviewers flagged #455 as clearing it *"more
easily than the docs sentence"*, and it was filed on that basis. Two independent readers agreeing
that a `let _ =` on a send looks unsafe was not a substitute for injecting a fault and watching what
happens — the same lesson as this arc's measurement incidents, one layer over: **the reviewers'
agreement raised confidence without adding evidence.** Its framing was also wrong in a second way:
the mid-stream sites do not surface the error either, they `break` and let the worker's own report
carry it.

## #458 closed by measuring it — #470 (2026-08-19)

`dev` `c7dd106`. `2G_c4` is now a recorded cell on every PR and push, non-gating and unbanded.

The cell was dropped from CI on a **gating** argument — a 672 MiB pool needs ~2.3x the fixture to
saturate, so it cannot fail the absolute promise informatively. That does not transfer to a
*recorded* cell, which only needs reproducibility, and it reuses the existing 51 bp fixture.

First Linux readings: **0.854 and 0.836** of the printed 1861 MiB (1589 / 1556 MiB peak), i.e.
0.776 / 0.760 of a 2048 MiB budget against macOS's 0.887 — an ~11 point offset, in family with the
11–17 seen elsewhere.

**Two of my own claims were corrected by the measurement**, both before merge:

- I documented that non-saturation would make the cell "read low". It reads **higher** than `1G_c4`
  (0.836–0.854 against 0.787–0.809 in the same runs), because absolute peak more than doubles with
  the pool. "Depressed relative to a larger fixture" is not "depressed relative to the other cells",
  and the comment claimed the second. Same shape as the 1.13-vs-1.101 threshold error: a comparative
  with an unstated referent.
- I had guessed it would be a weak drift detector because its prediction is larger. `B` is **42.0 MiB**
  here against 18.1 at `1G/c4`, so lever A moves it ~+0.09 to +0.13 — the same proportional
  sensitivity as the banded cell.

Also the first exercise of the margin assertion at a **second budget**: 2048 x 10/11 prints 1861,
837 KiB from exact, inside the ±2048 KiB tolerance. It had only ever run at 1 GiB.

**#455 was closed by another session** (#471, `a107833`) as a docs change stating why the EOF flush may
drop a failed send — the non-local invariant that makes it safe, which nothing had recorded.

**#457 is the only issue left open** in the arc.

## Band re-tuned and shipped — #469 (2026-08-18)

`dev` `d773131`. The band from #462 was calibrated from four reps per cell and was too tight: it
false-fired on an unrelated PR (`multipair_c4` 0.852 against a 0.840 trigger, moving alone while the
other cells went +0.012 and −0.019; re-run clean). Recomputed from **27 observations per cell** across
31 job logs, Validation 8's two regressed runs excluded:

| cell | n | mean | range | 3× width | banded |
|---|---:|---:|---:|---:|---|
| `1G_c4` | 27 | 0.7883 | 0.055 | 0.083 | **yes, 0.788** |
| `25bp_1G_c4` | 27 | 0.9138 | 0.043 | 0.065 | **yes, 0.914** |
| `multipair_c4` | 27 | 0.7987 | 0.081 | 0.121 | no |
| `1G_c4_fastqc` | 27 | 0.7069 | 0.080 | 0.120 | no |
| `1G_c2` | 25 | 0.8789 | 0.034 | 0.051 | no |
| `25bp_1G_c2` | 25 | 1.0301 | 0.043 | 0.065 | no |

Upper bound **+0.05 → +0.09**, lower bound unchanged. Affordable only because Validation 8 *measured*
lever A at +0.126–0.151 where the sizing had assumed the predicted +0.078.

**The rule is now two-sided**, which was the actual defect: `width ≥ 3× observed half-range` **and**
`≤ the measured regression delta − a margin`, with the half-range from a double-digit sample. As first
written it constrained only the noise side — and was evaluated using the very biased estimate it
existed to guard against, so it certified the cells it should have excluded. Verified against the
committed clauses: all 54 observations of the banded cells pass, both red at the measured lever-A
deltas, unbanded cells are not checked.

## #453 fully closed out — all three strings shipped

#466 merged as `bf28cbd`, discharging the condition #453 was closed under. Verified on `dev` by
running each case:

- `--gzip` on plain input → "Output compression follows the input, so this plain-text input produces
  plain-text output and no flag changes that." (#466)
- `--clumpify` on plain input → `(plain output)` plus the no-gain warning (#468)
- `--clump_only --dont_gzip` → `…, plain output)` (#468)

## Decided 2026-08-18

- **Item B bands three cells, not four.** Rule: ≥ 3 reps AND a false-red multiple ≥ 3×.
  `1G_c4_fastqc` is excluded on spread — and that costs little, because it stays promise-gated
  (`gating=1`, 0.745 against 1.000) and the drift signal is redundant: lever A moves every cell, and
  `multipair_c4` needs 11.2× its observed half-range to miss it.
- **A per-row margin column was rejected for now.** It invites sizing a width from n = 4, the same
  biased-low estimate that failed at n = 2. Held in reserve for a cell needing its own width for a
  structural rather than a noise reason.
- **Validation 8 (lever A with the band live) runs right after Item B** as a throwaway draft PR
  against `dev`, and is **#439's close condition**.
- **The RRBS `--memory 4G+` savings row** stays dated; its re-measurement is batched with #457's
  refit, since both need real libraries on oxy.
- **`1G_c4_fastqc`'s 4.8-point variance is not filed** — inside the budget, no promise at risk,
  nothing user-visible. A TSV note, not an issue.
- **#457 commented** ([#issuecomment-5325138291](https://github.com/FelixKrueger/TrimGalore/issues/457#issuecomment-5325138291)):
  the budget holds on Linux at `--cores 2` (0.932–0.937), the printed prediction does not
  (1.026–1.031). The prediction defect is general; the budget breach stays macOS-only.

## Open work after Phase 3

- **#457 is the live defect, and it is now measured on both platforms.** On **macOS** it exceeds
  the `--memory` budget at 25 bp / `--cores 2` (1.02–1.06×, 7 reps, two harnesses). On **Linux** it
  does **not**: measured 1.031 of the printed prediction but **0.937 of budget** (breach line
  1.101). So it breaks the number the banner prints on both platforms and the `--memory` guarantee
  on macOS only. σ is fitted per byte but is a spine factor, so short reads pack more records into
  the same bin bytes. Needs a calibration matrix across read lengths **on both platforms** — Linux's
  core increment is length-dependent where macOS's is not. Deliberately not bundled into Phase 3 (D7).
- **The CI baseline is in BOOTSTRAP mode** (`.github/clumpify-memory-baseline.tsv`): the drift
  band is inert until the values land and the marker is deleted. Only the absolute promise gates
  until then. **Two reps now exist**, and they measured identical code — the only diff between
  `69d8f67` (run 32065004600) and `3f1a930` (run 32065624366) is the TSV's comment block:

  | cell | rep 1 | rep 2 | spread |
  |---|---:|---:|---:|
  | `1G_c4` | 0.770 | 0.782 | 1.2 pts |
  | `1G_c4_fastqc` | 0.697 | 0.701 | 0.4 pts |
  | `25bp_1G_c4` | 0.941 | 0.924 | 1.7 pts |
  | `multipair_c4` | 0.792 | 0.787 | 0.5 pts |

  **The runner's spread is 0.4–1.7 points, not macOS's ~10**, so D5's `+0.05` bound is ≥ 3× the
  widest observed spread and the band is workable as specified.
- **The band is inert for two independent reasons, not one.** `want` comes from
  `awk '$1==c {print $2}'` over a file with only comment rows, and the guard is
  `[ -n "$want" ] && [ "$bootstrap" = 0 ]`. Deleting `# BOOTSTRAP` alone leaves the band inert
  *while looking activated*.
- **The platform offset is ~17 points, not ~9** (runner 0.770 vs macOS 0.939 at the gating cell),
  so the absolute promise has more slack on Linux than the plan assumed — which is the argument
  for the band, not against it.
- **#457 does not reproduce on Linux at `--cores 4`** (0.941, inside both bounds). The worst
  macOS cell, `--cores 2`, is still unmeasured there; #457 is commented with this.
- **RETRACTED (rep 4): the read-length interaction.** Rep 3 gave a +0.035 interaction against
  macOS's +0.007 and this file reported Linux's core increment as read-length dependent. Rep 4 gives
  **−0.009**, the opposite sign; the two-rep mean is +0.013 against a per-rep spread of 0.044, so it
  is not distinguishable from macOS. The driver was rep 3's low `1G_c2` reading (0.867 vs 0.894). The
  refit still wants both platforms, for the increment's **magnitude** (Linux mean +0.10 vs macOS
  +0.07 at both lengths), but that too needs more reps before it is asserted.
- **`1G_c4_fastqc` cannot be banded at ±0.05, and this blocks Item B as designed.** It read
  0.697/0.701/0.707 then **0.745** — range 4.8 points, the widest of all six cells, having looked
  like the tightest on two reps. Mean 0.7125 puts the trigger at 0.7625 against an observed max of
  0.745: **1.75 points of headroom on a 4.8-point range.** Banding it would red on correct code.
  Plausible mechanism: FastQC's cost scales with `--cores` and mimalloc retains the pool's pages
  across the trim→FastQC boundary, which is why Phase 2 made FastQC an addend rather than a `max`.
  Item B needs a per-row margin or a first-pass exclusion — the shape both reviewers anticipated
  while naming the wrong cell.
- **Banding verdict on four reps:** `multipair_c4` and `25bp_1G_c2` safe (20×/11×), `1G_c4` safe
  (6.7×/3.7×), `25bp_1G_c4` marginal (3.4×/1.9×), `1G_c4_fastqc` unsafe (2.1×/1.2×). Spreads widened
  for three of four cells as reps landed — n=2 bounds dispersion, it does not estimate it.
- **#452 was stacked on #451's *first* commit, not its tip**, so squash-merging #451 left it
  `CONFLICTING` on `CHANGELOG.md` (`src/clump.rs` auto-merges clean). Predicted by
  `git merge-tree`, then confirmed by GitHub. Resolving that conflict toward #452 would revert
  `44dc732` — the commit pointing the short-read caveat at #457. Fix is
  `git rebase --onto origin/dev 1b2579e ci/clumpify-memory-guard-439`, guaranteed clean because
  `CHANGELOG.md` is the only file differing between `1b2579e` and #451's tip and #452 never
  touches it. Planned as Item A in `phase3-close/PLAN.md`.
- **Phase 3's Validation 12 — "the drift band catches lever A" — has never been run**, because
  BOOTSTRAP made it untestable. It carries forward as `phase3-close` Validation 7. Until it runs,
  D5's band is asserted-to-work rather than proven-to-work.

## History

- 2026-08-18: **#462 merged** (`45876c5`) — drift band live; **Validation 8 discharged in both forms**, #439's close condition met
- 2026-08-18: `phase3-close/PLAN.md` → **v3** (three-cell banding rule; interaction claim retracted); #457 commented with the Linux `--cores 2` result
- 2026-08-18: **#452 merged** (`bc48bf0`) after Item A's rebase + two `--cores 2` cells; reps 3–4 measured; Item B blocked on the FastQC cell's spread
- 2026-08-17: **#448 merged** (`efc1b11`), **#451 merged** (`242546a`); #452 left `CONFLICTING` by the stale stack; `phase3-close/PLAN.md` written, dual review launched
- 2026-08-17: Phase 3 → v3 (§0 retracted, CI gate → baseline band per D5, issue narrowed per D6); round-2 review partial
- 2026-08-17: Phase 3 plan review ×2 → ✅ Complete; PLAN.md → v2 (both guards rebuilt; D3, D4 decided)
- 2026-08-17: Phase 3 → 📝 Planning (`phase3/PLAN.md`; D1 and D2 decided)
- 2026-08-17: Phase 2 implementation → ✅ Complete (calibrated model, FastQC per-core term)
- 2026-08-17: Phase 2 calibration → ✅ Complete (both platforms, distinct-read fixture)
- 2026-08-16: Phase 1 → ✅ Merged as `1595d59` (PR #442, all 9 checks green)
- 2026-08-16: Code Review + Coverage (Phase 1) → ✅ Complete
- 2026-08-16: Plan v3, spike, plan reviews (v1 ×2, v2 ×2) → ✅ Complete
