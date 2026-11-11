# Plan Coverage Report — `total_bp_after_trim` parity between the serial and parallel paths

**Plan(s):** PLAN.md (revision 2)
**Verdict:** COMPLETE

**Mode:** B — code vs. plan (post-implementation; the plan carries the implementation outline, there
is no separate `IMPL.md`)
**Date:** 2026-11-10
**Audited:** working tree at `/Users/fkrueger/Github/TrimGalore`, branch `dev`, base `c0bbee4`,
change uncommitted (`git diff --numstat`: `src/parallel.rs` +202/−5, `CHANGELOG.md` +10,
`docs/src/content/docs/reference/changelog.md` +10)

## Summary

- Total items: 18 (9 outline steps + 9 §9 validations)
- DONE: 17
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 1 (Step 5 — declared in §12, deviation 1)

No source file was edited by this audit; the tree remains on `dev` with the change in place.

## Coverage ledger — §5 implementation outline

| # | Item | Status | Evidence |
|---|------|--------|----------|
| 1 | Three accumulation sites switched to `bp_after_cutadapt` | DONE | `src/parallel.rs:542`, `:543`, `:1157`. All nine sites across `parallel.rs` + `trimmer.rs` now read `bp_after_cutadapt` (serial FASTQ `trimmer.rs:350`, `:536`, `:537`; serial uBAM `:672`, `:780`, `:781`; parallel `:542`, `:543`, `:1157`). The three sites are exhaustive: `trim_read` is called at `parallel.rs:495`, `:496`, `:1133` only, and `classify_single`/`classify_paired` have one caller each (`:1194` in `process_reads`, `:628` in `process_pairs`). |
| 2 | One-line invariant comment at each site | DONE | `src/parallel.rs:541` and `:1156`, identical text, one line, current-state only. Serial comments untouched. |
| 3 | `parity_test_config_post_adapter()`, **all four** clip values | DONE | `src/parallel.rs:1456-1465`; `clip_r1: Some(5)` `:1458`, `clip_r2: Some(5)` `:1459`, `three_prime_clip_r1: Some(3)` `:1460`, `three_prime_clip_r2: Some(3)` `:1461`. The review's Critical finding (R2 undiscriminated) is closed — see item 5. `parity_test_config()` (`:1424`) is byte-for-byte unchanged and still has 12 non-helper callers (assumption 5 holds). |
| 4 | Single-end parity test strengthened | DONE | New sibling `test_parallel_serial_bp_after_trim_parity_with_clipping` at `src/parallel.rs:1547`. The existing clean-config run survives: `test_parallel_serial_trim_stats_parity` (`:1474`, config at `:1496`) is unmodified and passes. |
| 5 | Paired parity test covering **both** R1 and R2 | DEVIATED (declared) | New sibling `test_parallel_serial_bp_after_trim_parity_paired` at `src/parallel.rs:1622` instead of extending `test_passthrough_serial_parallel_parity`. R1 asserted at `:1693-1697`, R2 at `:1699-1703`, with separate per-side anchors at `:1709-1712`. §10 delegated extend-vs-sibling to the implementer; deviation declared in §12. The passthrough test (`:2468`) is untouched and passes. |
| 6 | Absolute-value anchor, hand-derived | DONE | `src/parallel.rs:1597` (`EXPECTED_BP_AFTER_TRIM = 15*30 + 15*34`), `:1707-1708` (`EXPECTED_R1 = 15*30 + 15*34`, `EXPECTED_R2 = 30*33`). **Independently re-derived from the fixture literals, not read back from the code:** the adapter-bearing R1 string is 43 bp and ends in exactly the 13 bp `AGATCGGAAGAGC`, leaving a 30 bp prefix; the clean R1 string is 34 bp; the R2 string is 33 bp. Hence 960 (SE and PE R1) and 990 (PE R2). The derived `total_bp_processed` (15×43 + 15×34 = 1155) and `total_bp_written` (15×22 + 15×26 = 720; R2 30×25 = 750) match §12's mutation output exactly. Anchors are asserted on **both** paths in both tests. |
| 7 | Fixture discrimination guard in both tests | DONE (form substituted) | `src/parallel.rs:1609-1612` (SE) and `:1716-1720` (PE, both reads conjoined). The implementation asserts `total_bp_after_trim > total_bp_written` where the plan's snippet showed `total_bp_after_trim < total_bp_processed`. The substituted form matches Step 7's stated intent (that the config still makes *the two quantities* differ) rather than only proving adapter trimming occurred; it is sound here because `length_cutoff: 0` means no read is filtered, so `total_bp_written` is the sum of final clipped lengths. §12 records the actual form; it is not listed among the declared deviations. |
| 8 | CHANGELOG entry | DONE | `CHANGELOG.md:6-14` under the previously empty `### Unreleased`, as a new `#### Bug fixes` subsection — the same subsection name the file already uses (`:233`, and `#### Changes` at `:18`). One sentence after the bold title, in the file's `- **Title** ([#NNN](url)). Prose.` format, with no migration prose. It also carries the §10-sanctioned mention of #232 rather than rewriting that entry. |
| 9 | Mutation checks confirming the tests discriminate | DONE (recorded, not re-run here) | §12 records four; re-running requires editing a source file, which this audit is barred from doing. Each is logically possible as specified and would fail for the stated reason — see "Step 9" below. |

## Coverage ledger — §9 validation table

| # | Validation | Status | Evidence |
|---|---|---|---|
| 1 | SE figure no longer varies with `--cores` | DONE — independently verified | Post-fix reports on disk: `scratchpad/bpexp2/clip_c{1,4}/…_trimming_report.txt:28` both `644,525 bp (99.2%)`; `bpexp2/rrbs_c{1,4}/…:29` both `490,527 bp (96.2%)`. Pre-fix `bpexp/clip_c4` shows `594,545 bp (91.5%)`, reproducing §2's measured delta. |
| 2 | **Paired** figure no longer varies, both report files | DONE — independently verified | `bpexp2/pair_c{1,4}/BS-seq_10K_R1…:28` both `644,525`; `…_R2…:28` both `645,073`. R1 ≠ R2, so `:542` and `:543` are separately exercised end to end. |
| 3 | Records untouched | DONE — independently verified | `gzip -dc | shasum -a 256` over all four clip runs (pre-fix ×2, post-fix ×2): `53d619cd94e5…` identical, 9,996 records. All four RRBS runs: `f49dea91ec11…`, 9,729 records. Matches §12's quoted prefixes. |
| 4 | `total_bp_written` unchanged | DONE — independently verified | `total_bp_written: 594545` in the JSON report for `bpexp/clip_c1`, `bpexp/clip_c4`, `bpexp2/clip_c1`, `bpexp2/clip_c4`. |
| 5 | SE site genuinely covered (revert `:1157`) | DONE — recorded in §12 | Not re-run (no source edits permitted). Logically necessary: with `:1157` reverted, parallel `total_bp_after_trim` becomes 720 against serial 960, failing both the whole-struct `assert_eq!` (`:1587`) and the parallel anchor (`:1602`). |
| 6 | **Each** paired site genuinely covered (`:542`, then `:543`) | DONE — recorded in §12 | Not re-run. Reverting `:542` alone breaks only the R1 assertions (960 vs 720, `:1693`/`:1709`); reverting `:543` alone breaks only R2 (990 vs 750, `:1699`/`:1710`). Both numbers follow from the fixture independently of §12. This is the finding revision 1 could not have caught. |
| 7 | The anchor catches a serial regression | DONE — variant | §9 specifies flipping `trimmer.rs:350` alone; §12 flipped `trimmer.rs:350` **and** `parallel.rs:1157` together. The variant is strictly more informative — flipping the serial line alone would also break the `serial == parallel` assertion, so it would not isolate the anchor's unique value, whereas the pair makes both paths agree at 720 and leaves only the anchor to fire. Recorded in §12's mutation table, not in its deviations list. |
| 8 | No collateral damage, `cargo test` | DONE — partially re-verified | Implementer's log `scratchpad/full_suite_final.log`: 21 `test result: ok` lines summing to 844 passed, 0 failed, 0 ignored, no `FAILED`/`panicked` occurrences; both new tests appear as `... ok` (`:419`, `:421`), as do the two pre-existing parity tests (`:423`, `:523`). Independently re-ran the two new tests in an isolated `CARGO_TARGET_DIR`: `2 passed; 0 failed; 516 filtered out`, exit 0. |
| 9 | CI gates | DONE — independently verified | `cargo fmt --all -- --check` exit 0, empty output. `cargo clippy --all-targets -- -D warnings` exit 0, zero `warning`/`error` lines (dev profile, isolated target dir); the implementer's `scratchpad/clippy.log` shows the release-profile run finishing clean. Local rustc 1.95 vs CI 1.99 remains the stated caveat. |

## Scope boundaries

Respected. The diff touches three files and, inside `src/parallel.rs`, only: the two comments, the
three accumulation statements, one test helper and two tests. Nothing from **#542**
(`MultithreadedWriter`) and nothing from **#529** (uBAM threading) appears — no writer, BAM or
threading code is modified, consistent with §7's finding that #529 inherits the fix through
`classify_single`/`classify_paired` on rebase.

## Beyond the plan / undeclared extras

No undeclared source additions. Both declared deviations are real and accurately described:

1. **Focused paired test instead of extending the passthrough test** — confirmed. The new test
   drives `crate::trimmer::run_paired_end` and `run_paired_end_parallel` with `None` in the three
   passthrough slots, and `test_passthrough_serial_parallel_parity` is unmodified and still passes.
2. **Docs changelog mirrored** — confirmed and correctly motivated. At `c0bbee4` the two changelog
   bodies were already in sync (`diff` of root tail-from-line-2 against the docs file from line 9:
   one blank line only); the root-only edit created the desync, and the mirror restored byte parity
   on the same single-blank-line offset. The docs page does declare itself a copy in its own preamble
   (`docs/src/content/docs/reference/changelog.md:6-8`).

Three smaller divergences from the outline's literal text, none of which affects coverage and all of
which are visible in §12's step table rather than its deviations list:

- Step 7's guard expression (`> total_bp_written` rather than `< total_bp_processed`) — see item 7.
- Step 5 said to set `length_cutoff` as `passthrough_test_config()` does (20); the new paired test
  leaves it at `parity_test_config()`'s 0 and sets only `is_paired: true` (`:1645-1647`). With these
  read lengths nothing is filtered at either value, so the anchors and the guard are unaffected.
- Validation 7 was run as a two-line mutation rather than the specified one-line flip — see item 7.

## Does §12 match the code?

Yes, with one numeric slip:

- **Wrong:** §12 line 346 gives the diff as `src/parallel.rs` **(+207/−5)**. `git diff --numstat`
  reports **+202/−5**; 207 is `--stat`'s total of changed lines (202 + 5). The other two file
  counts (+10, +10) are correct.
- Every other §12 claim checked out: the three sites and the nine-site consistency, the one-line
  comments, all four clip values in the helper, both new test names, the hand-derived anchors
  (960 / 960 / 990, independently re-derived here from the fixture literals), the discrimination
  guards, the changelog entry and its docs mirror, and all six verification rows — the SE and PE
  figures at both core counts, the two output hashes and record counts, `total_bp_written` 594,545,
  the 844/0 suite result, and both CI gates.
- §12's mutation table could not be re-executed under this audit's no-edit constraint; the two code
  reviewers may re-run it. Its four reported figures (960/720 and 990/750) are the ones the fixture
  arithmetic predicts, so they are at least self-consistent with the code as it stands.

## Verdict

**COMPLETE.** All nine outline steps are present in the code; all nine §9 validations are satisfied,
six of them re-verified independently here (1–4, 8, 9) and three resting on §12's record plus a
logical check of the assertions they target (5–7). The one declared behavioural deviation is within
the latitude §10 granted, and the only inaccuracy found in §12 is the `+207/−5` line count.
