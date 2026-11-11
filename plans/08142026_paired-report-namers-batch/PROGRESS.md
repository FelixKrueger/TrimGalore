# PROGRESS — #423 → #424 → #429, one PR

**Plan:** `PLAN.md` v3 (dual review folded in; O1 decided)
**Issues:** [#423](https://github.com/FelixKrueger/TrimGalore/issues/423),
[#424](https://github.com/FelixKrueger/TrimGalore/issues/424),
[#429](https://github.com/FelixKrueger/TrimGalore/issues/429)
**Branch base:** `dev` @ `d2ea293`
**Test baseline:** 707 macOS / 706 Linux

| Step | State | Notes |
|---|---|---|
| Sites surveyed | **done** | Eight `pair_output_dir` callers outside `io.rs` (review corrected my six); four in scope, four a different scheme |
| Premises verified | **done** | Neither interleaved arm reads `cli.basename`; the three `Arm::Unplanned` cases located |
| Plan written | **done** | v1, superseded by v2 |
| Felix's manual review | **done** | 2026-08-14: "approve" |
| Dual plan review | **done** | `PLAN_REVIEW_A.md` + `PLAN_REVIEW_B.md`; both "not safe as written", 4 Criticals each, 9 in common |
| Plan revised to v2 | **done** | §12 records what changed |
| O1 decided | **done** | Felix, 2026-08-14: **the struct**. Plan committed to it as v3 |
| Implementation — #423 | **done** | Commit `6569652` on `fix/423-paired-report-namers`. 707 → **714** tests, fmt + clippy clean. V1a + V1b both identical |
| Implementation — #424 | not started | Awaiting confirmation to continue in the same branch |
| Implementation — #429 | not started | Independent of the other two |
| Dual code review + coverage audit | not started | Automatic after implementation |
| PR / merge / close ×3 | not started | Squash-merge into `dev`; manual close of all three (master is default) |

## Open items carried in the plan

- **O1** — **answered: the struct.** `InterleavedReports` + `InterleavedFastqOutputs` +
  `InterleavedBamOutputs`, with `--retain_unpaired` and `--no_report_file` as `Option` fields and a
  `planned()` that opens with an exhaustive destructure, so a field added later cannot escape the
  candidate list without a build error.
- **O2** — **closed** in favour of wrapping: two distinct noodles failure modes, and the name is
  reachable at the wrap site.
- **O3** — which `io::ErrorKind` does each #429 failure carry? Resolve by probe as the first
  implementation action; affects only how much the message explains.
- **O4** — one PR, three commits. Assumption taken.

## Resolved by the review, no longer open

- The old O3 (one-input collision message) is **resolved**: it reads correctly under
  `OutputSource::Input`. It would have read `one.bam + one.bam` under `Pair`, which v1's step 5
  actively prescribed.

## Scope corrections found while planning

- **#423's "four places" is four in-scope sites out of eight `pair_output_dir` callers** at
  `d2ea293` (my v1 said six; the review found `clump_only.rs:537` and `:1089` too). `main.rs:619`
  plans **two** per-mate clumping reports and `:746` plans **one** keyed on R1 — v1 had those two
  descriptions swapped. All four clump sites are a different scheme, so folding them in would be a
  behaviour change, not a refactor.
- **#429 has narrowed by half.** Its body pairs the opaque message with a partial output file;
  #428 (`d2ea293`) fixed the residue, so only the message remains.

## The thing most likely to go wrong

**Composing the interleaved namer from `io.rs::report_name`.** `report_name` keys on
`file_name()`; both interleaved arms key on `file_stem()`. Divergent by design, so the natural
composition renames every interleaved report `ip.bam_R1_trimming_report.txt` and breaks MultiQC
globs. Nothing in `tests/`, `docs/` or `.github/` pins those four filenames — zero grep hits — and
after step 5 the writer and the candidate list both derive from the same namer, so a joint rename
keeps `created ⊆ planned` true and every promoted case green. **V1b is the only guard.**

Second most likely: a candidate list that hardcodes `.fq`. All three promoted cases pass while
every `.gz`-named-input run is mis-planned.

## Implementation note (#423 commit)

`planned()` is **not** in this commit. It is only consumed by #424's pre-flight, and an unused
method fails CI under `-D warnings` — the same trap the review found for the tripwire's
`fn unplanned`. It lands in the #424 commit, where it has a caller.

Two things worth carrying forward:

- **The baselines had to be recaptured.** `target/release` still held the PR-436 build, so the
  first attempt would have compared #423 against someone else's patch. Rebuilt from `dev` first.
- **My two failing unit tests were wrong, not the code.** `Path::new("ip.bam").parent()` is
  `Some("")`, not `None`, so `pair_output_dir`'s `unwrap_or(".")` never fires and a bare input
  name yields bare output names. The created-file baseline already recorded that; the test
  expectations were written from intent instead of from the artifact.
