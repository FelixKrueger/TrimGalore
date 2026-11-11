# Plan Coverage Report

**Mode:** B (code vs. plan — `PLAN.md` v2 carries the numbered implementation outline, so the design
plan *is* the spec; there is no separate `IMPL.md`)
**Plan:** `plans/08162026_pr441-library-presets/PLAN.md` (v2, 2026-08-16, incl. §11 dispositions and
§12 implementation notes)
**Code:** `9a738a333f1c2c75219b4f15ce2ebff418447495` on `fix/440-library-refusals`, parent `872c4ae`
(`feat(cli): --library presets expand into end clipping (#440) (#441)`, author `Benjamin Demaille`)
**Worktree audited:** `/private/tmp/claude-501/impl440` (read-only; nothing modified)
**Diff:** `git diff 872c4ae 9a738a3` — 9 files, +382 / −25
**Date:** 2026-08-17
**Verdict:** COMPLETE

## Summary

- Total items: 100 — 34 outline steps (A), 13 behaviour clauses and edge cases (B), 17 validation rows
  (C), 7 target findings (D), 8 assumptions (E, tallied as DONE where UPHELD), 8 §11 dispositions (F),
  4 deviations (G), 9 further plan claims (H)
- DONE: 95
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 5 — V3 plus D1–D4, all documented in §12; each description checked against the code

Measured in this audit, not quoted from §12:

| Check | Command | Result |
|---|---|---|
| Suite | `cargo test --no-fail-fast` | **764 passed, 0 failed**, 17 result lines (15 with tests + `main.rs` and doc-tests at 0) |
| New tests | `#[test]` count per file, `872c4ae` → `9a738a3` | **+12** (cli +5, library +2, report +2, integration_library +2, integration_ubam_out +1) |
| Formatting | `cargo fmt --all -- --check` | clean |
| V1 | `cargo test library_plus_` | 4 passed |
| V5 | the two tests run separately (see note) | 3 passed; both plan-named test names resolve |
| V4 | release binary, all four modes, empty dir | exit 1, correct message, **0 files** ×4 |
| V15 | release binary vs. the docs transcripts | every quoted line verbatim; 0 `--clip_R2` in SE stderr and text report |
| Clippy | not run | **disk at 3.6 GB on a 100%-full volume**; `--all-targets --release` would rebuild 16 test binaries under the release profile. §12's "clean" is author-reported and unverified here. `cargo fmt --all -- --check` and the full test suite were both run. |

The 764 total decomposes as 481 (lib) + 15 + 13 + 16 + 11 + 8 + 1 + 8 + 57 + 14 + 21 + 3 + 57 + 19 + 40
= 752 base + 12 new, matching §12 exactly.

**Mutation rows V2, V6, V6a, V7, V8–V13 are author-reported.** Re-running them means editing source in
`impl440`, which this audit is forbidden to do, and no spare worktree is affordable at 3.6 GB free. For
each, the *target test* was confirmed to exist and pass, and the *property* was independently verified
where a read-only route existed (V3, V4, V15, V16 below). This is recorded as a limit on the evidence,
not as a doubt about the claims.

---

## A. Implementation outline — steps 1–23

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| A1 | `flag_summary` gated on `paired`; R2 entries skipped when `!paired`; array becomes a 3-tuple | Step 1 | DONE | `src/library.rs:118-137`. Body matches the plan's code verbatim |
| A2 | The two justification sentences deleted from `flag_summary`'s doc comment | Step 1 | DONE | Old "Printed in the log and in both trimming reports, so a run stays reproducible…" replaced by the plan's exact two-line comment. Satisfies part of §11's C-1 |
| A3 | `src/main.rs:1121` → `clips.flag_summary(cli.paired)` | Step 2 | DONE | `src/main.rs:1121` |
| A4 | `src/report.rs:279` → `library.flag_summary(config.paired)` | Step 2 | DONE | `src/report.rs:276-283` |
| A5 | `cargo fmt` run after step 2's reformat | Step 2 | DONE | `report.rs`'s `writeln!` is split across four lines; `cargo fmt --all -- --check` clean (verified) |
| A6 | Refusal block inserted after the `fastqc_requested` block, before the deprecation warnings | Step 3 | DONE | `src/cli.rs:1111-1133`; the `// The specialty arms return before setup_trimming…` comment matches the plan verbatim |
| A7 | Four `bail!`s naming `--hardtrim5` / `--hardtrim3` / `--clock` / `--implicon`, wording "does not honour the clip flags; --library is not compatible" | Step 3 | DONE | `cli.rs:1115/1120/1125/1130`. Wording confirmed against the running binary, all four modes |
| A8 | `library_plus_hardtrim5_rejected` | Step 4.1 | DONE | `cli.rs:2418`; asserts the error contains both `--library` and `--hardtrim5` |
| A9 | `library_plus_hardtrim3_rejected` | Step 4.2 | DONE | `cli.rs:2428` |
| A10 | `library_plus_clock_rejected` (`--clock --paired … R1 R2`) | Step 4.3 | DONE | `cli.rs:2438`; argv shape as specified |
| A11 | `library_plus_implicon_rejected` — **bare `--implicon`**, no separate value | Step 4.4 | DONE | `cli.rs:2458`; bare flag as specified, with a doc comment stating why (`require_equals`) |
| A12 | Section header comment `// ── #440: --library refused where no clip site is reachable ──` | Step 4 | DONE | `cli.rs:2415`, placed immediately after the four `fastqc_plus_<mode>_rejected` tests the plan anchored on (`cli.rs:2313-2351`) |
| A13 | `hardtrim5_range_error_precedes_the_library_refusal` | Step 4.5 | DONE | `cli.rs:2478`; asserts `between 1 and 999` **and** `!err.contains("--library")` — stronger than the plan's wording |
| A14 | `test_json_report_records_library_block` kept unchanged (cardinality-1 coverage) | Step 6 | DONE | Zero diff hits for the name; still at `report.rs:1456` with `emseq` / `clip_r1: 12` |
| A15 | `test_json_report_records_all_four_overrides` — `accel` + 1/2/3/4, four block keys asserted | Step 7 (T-1) | DONE | `report.rs:1491`; asserts `clip_r1..three_prime_clip_r2` = 1,2,3,4 |
| A16 | Same test asserts the four `ClipOverride` records' flag names and displaced preset values 10/15/10/10 | Step 7 (T-2) | DONE | `report.rs:1521-1533`; ordered `expected` array, `overrides.len() == 4` |
| A17 | `test_text_report_omits_r2_clipping_on_single_end` | Step 8 | DONE | `report.rs:1404`; asserts the `pbat` SE line and **absence of both** `--clip_R2` and `--three_prime_clip_R2` |
| A18 | `test_text_report_records_preset_name_and_expanded_values` stays green unchanged as the paired half | Step 9 | DONE | `report.rs:1380`, `config.paired = true`, untouched by the diff (0 diff hits), asserts all four |
| A19 | `test_json_report_library_is_null_without_a_preset` gains a `contains_key` presence assertion | Step 10 (F-5a) | DONE | `report.rs:1552-1559` |
| A20 | Two direct unit tests for `flag_summary` — one paired, one single-end | Step 11 | DONE | `flag_summary_lists_all_four_for_paired` and `flag_summary_omits_read_2_for_single_end`, `src/library.rs:287/296`. Landed beside the function in `library.rs` rather than under step 11's `report.rs` heading; substance is exactly as specified |
| A21 | `both_reports_record_the_preset_and_its_expanded_values` expectation changed to `--clip_R1 8 --three_prime_clip_R1 8` | Step 12 | DONE | `tests/integration_library.rs:194-197`, now anchored on the full `Clipping in force: ` prefix |
| A22 | F-5b agreement loop added to the same test, all four keys | Step 13 | DONE | `integration_library.rs:207-220`; asserts block == top-level **and** == 8 per key |
| A23 | `explicit_clip_flag_overrides_the_preset_and_is_logged` gains the stderr negative assertion | Step 14 (C-2) | DONE | `integration_library.rs:147-152`; `!err.contains("--clip_R2") && !err.contains("--three_prime_clip_R2")` — stronger than the plan's single clause |
| A24 | `a_run_without_a_preset_reports_no_library_block` gains the presence assertion | Step 15 | DONE | `integration_library.rs:233-241` |
| A25 | New `paired_report_lists_all_four_values` — text report **and** stderr banner both list all four | Step 16 | DONE | `integration_library.rs:258`; real `--paired --library accel` run, both asserted |
| A26 | New `preset_with_multiple_overrides_round_trips_through_json` | Step 17 | DONE | `integration_library.rs:287`; real binary, four overrides, asserts `len == 4` plus flag / preset_value / user_value per entry |
| A27 | T-3 test: `--rename --output-format ubam --library emseq` on FASTQ input refused, **error text** asserted, whole output listing empty | Step 18 | DONE | `tests/integration_ubam_out.rs:1044`; asserts `--rename is refused` + `--output-format ubam`, and `dir_listing(&dir) == vec!["sp.fastq"]` — the exact pattern at `integration_ubam_out.rs:1032-1038` |
| A28 | Docs: override transcript at `library-presets.md:48-54` replaced with the real two-flag output | Step 19 | DONE | Now `--clip_R1 12 --three_prime_clip_R1 10`; reproduced verbatim from the binary (V15) |
| A29 | Docs: text-report transcript at `:62-66` same fix | Step 20 | DONE | `Clipping in force: --clip_R1 12 --three_prime_clip_R1 10`; verbatim (V15) |
| A30 | Docs: the falsified `:60` prose rewritten, and the text/JSON explanation sited in the 60-79 passage | Step 21 | DONE | `library-presets.md:62-64` now says the text report "lists the values in force for the run, which for this single-end example is Read 1 only" — inside the passage. A second paragraph after the JSON block (`:86-88`) explains the asymmetry and points at the top-level `mode` key |
| A31 | Docs: "What a preset does *not* set" gains the specialty-mode bullet | Step 22 | DONE | `library-presets.md:42` |
| A32 | `flags.md:14` refusal clause — decision taken | Step 23 / O-4 | DONE | Clause added per the plan's own recommendation. §12 logs this as D2; implementing the plan's stated recommendation is conformance, not deviation |
| A33 | CHANGELOG: the falsified "all four values in force" sentence rewritten | Step 9 | DONE | 0 occurrences of the old sentence; the replacement states all four in JSON and "in the text report the values in force for the run, so a single-end report lists Read 1 only" |
| A34 | CHANGELOG: refusal folded into the "A preset is only a clip bundle…" paragraph, **no second entry** | Step 9 | DONE | Single hunk, 0 lines matching `^\+- \[` (no new bullet). Unique issue refs 70 → 70, so R-4's dropped-entry risk did not materialise |

## B. Behaviour clauses B1 and B2, with every edge case

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| B1 | The new block is the **last fallible check** in `validate()` | B1.1 | DONE | Verified structurally: the last `bail!` in `cli.rs` is now `1130` (ours); everything from `1134` to `Ok(())` is `eprintln!` deprecation warnings. Every precedence guarantee follows |
| B2 | Bails when `library.is_some()` and any of the four modes is set, naming the mode and `--library` | B1.2 | DONE | Confirmed against the running binary for all four modes |
| B3 | Exit non-zero and **no output file created** | B1.3 | DONE | Independently re-verified (V4): all four modes exit 1 with 0 files in the output dir |
| B4 | `--hardtrim5 0` reports the `1 and 999` range error first | B1 edge | DONE | Binary: `Error: --hardtrim5 must be between 1 and 999, got 0`; pinned by `hardtrim5_range_error_precedes_the_library_refusal` |
| B5 | `--hardtrim5` + `--hardtrim3` reports #386's bail first | B1 edge | DONE | Structural: `cli.rs:1009` precedes `:1111`. No test required by the plan |
| B6 | `--clock` / `--implicon` arity errors report first | B1 edge | DONE | Structural: `cli.rs:1027`/`:1030` precede `:1111`. No test required by the plan |
| B7 | Missing input files report first | B1 edge | DONE | Structural: `cli.rs:1042-1044` precede `:1111`. No test required by the plan |
| B8 | `--clump_only` short-circuits before the new block | B1 edge / Assumption 8 | DONE | Binary: `--clump_only --hardtrim5 30 --library emseq` → `Error: --clump_only does not clip; --library is not compatible` (`cli.rs:850`). The new block is unreachable there, exactly as stated |
| B9 | `flag_summary` takes `paired: bool`; R1 entries unaffected in both modes | B2.1/B2.2 | DONE | `flag_summary_lists_all_four_for_paired` + `flag_summary_omits_read_2_for_single_end` pin both directions |
| B10 | Paired output byte-identical to today | B2.3 | DONE | `test_text_report_records_preset_name_and_expanded_values` (unchanged, paired, all four) stays green; `paired_report_lists_all_four_values` asserts the exact 4-flag string in report **and** banner. Only report/stderr text is touched — read bytes are not on this path |
| B11 | Single-end becomes `--clip_R1 8 --three_prime_clip_R1 8` for `pbat` | B2.3 | DONE | Asserted in `test_text_report_omits_r2_clipping_on_single_end` and in `both_reports_record_…` |
| B12 | **Override log lines deliberately NOT gated** | B2 edge | DONE | `main.rs`'s `for o in &clips.overrides` loop is untouched. Verified live: an SE `--library emseq --clip_R1 12` run still prints `--clip_R1 12 was given on the command line and overrides the emseq preset value 10`, and the docs transcript keeps that line under the gated banner |
| B13 | **The JSON keeps all four values on single-end** | B2 edge / §8.3 | DONE | Verified live on an SE run: `mode: single-end`, library block and top level both `12/10/10/10`, and they agree — which is the F-5b invariant the gate would have broken. Pinned by A22 |

## C. Validation rows V1–V16 (plus V6a)

Legend for the evidence column: **re-run** = executed in this audit; **§12** = evidenced by the plan's
own 11-row mutation matrix, not reproducible read-only.

| # | Row | Status | Evidence | Notes |
|---|---|---|---|---|
| C1 | V1 — refusal fires on all four modes | DONE | re-run | `cargo test library_plus_` → 4 passed, 0 failed |
| C2 | V2 — those tests can fail | DONE | §12 | Refusal block disabled → all 4 `library_plus_*_rejected` red |
| C3 | V3 — ordering is pinned | DEVIATED (documented, D3) | §12 + re-run | §12 suppressed the range `bail!` rather than relocating the refusal block. The *property* was independently re-verified: `hardtrim5_range_error_precedes_the_library_refusal` passes, and the binary reports `--hardtrim5 must be between 1 and 999, got 0` with no `--library` mention. Structurally airtight too — `cli.rs:1019` precedes `:1111` |
| C4 | V4 — no output on refusal | DONE | re-run | All four modes in an empty dir: exit 1, correct message, **0 files** each time. Also re-confirmed with the `--clump_only` short-circuit case |
| C5 | V5 — the SE gate works, and the names match real tests | DONE | re-run | Both plan-named tests exist and pass. **The plan's command as written cannot work** — `cargo test A B` passes two filters and matches nothing (§12 iteration 2 hit exactly this). Run separately: `test_text_report_omits_r2_clipping_on_single_end` 1 passed, `flag_summary_` 2 passed. Plan-command defect, not an implementation gap |
| C6 | V6 — the SE gate can fail | DONE | §12 | `flag_summary` ignoring `paired` → 4 red (both new unit tests + `both_reports_record_…` + `explicit_clip_flag_overrides_…`) |
| C7 | **V6a** — the stderr half can fail | DONE | §12 | Reverting **only** `main.rs:1121` → exactly 1 red: `explicit_clip_flag_overrides_the_preset_and_is_logged`. This is the C-2 gap and it is closed; the assertion exists at `integration_library.rs:147-152` |
| C8 | V7 — the paired half can fail | DONE | §12 | R2 dropped unconditionally → 3 red incl. `paired_report_lists_all_four_values` |
| C9 | V8 — T-1 closed | DONE | §12 | JSON `clip_r2`/`three_prime_clip_r1` **key** labels swapped → `test_json_report_records_all_four_overrides` red; §12 records that this mutation survived all 751 pre-change tests |
| C10 | V9 — T-2 closed | DONE | §12 | `--clip_R2`/`--three_prime_clip_R1` **flag** labels swapped in `pick` → 2 red |
| C11 | V10 — T-3 closed | DONE | §12 | #408 guard reverted to the explicit clip list → `rename_into_ubam_refused_for_fastq_input_via_library_preset` red |
| C12 | V11 — F-5a closed | DONE | §12 | JSON `library` key deleted → both presence tests red (both stayed green before) |
| C13 | V12 — F-5b closed | DONE | §12 | Top-level `three_prime_clip_r2` +100 → the agreement loop red |
| C14 | V13 — cardinality 1 still covered | DONE | §12 | Overrides comma guard removed → 3 red incl. `test_json_report_records_library_block`, vindicating step 6's decision to keep it |
| C15 | V14 — suite and lints | DONE | re-run (tests, fmt); §12 (clippy) | **764 passed / 0 failed / 17 result lines** and `cargo fmt --all -- --check` clean, both measured here. Clippy `--all-targets --release` was **not** run: 3.6 GB free on a 100%-full volume, and the task forbids new builds. §12's "clippy clean" stands unverified by this audit |
| C16 | V15 — docs transcripts are real | DONE | re-run | Ran `--library emseq --clip_R1 12` single-end against the release binary. All four quoted lines reproduce verbatim (`grep -qF --`): the banner, the override line, `Library preset: emseq`, `Clipping in force: --clip_R1 12 --three_prime_clip_R1 10`. The docs' JSON block also matches the emitted JSON (12/10/10/10, one override). 0 `--clip_R2` in stderr or the text report |
| C17 | V16 — `--rrbs` unaffected | DONE | §12 (D4) + re-run of the mechanism | §12: 5 shapes, 24 files, all FASTQ and `.txt` byte-identical; JSON identical once `command_line` is excluded (D4). Independently re-verified the *mechanism* the plan itself calls the real argument: an `--rrbs --paired` run writes **0** `Library preset` / `Clipping in force` lines in either report, so the new line cannot enter those bytes |

**No V-row is without evidence.** V5 is the one row whose plan-side command is broken; its substance is
re-verified here.

## D. The seven findings this work was meant to close

| # | Finding | Status | Test(s) that close it |
|---|---|---|---|
| D-a | Specialty-mode refusal (both code reviewers) | DONE | `library_plus_hardtrim5_rejected`, `library_plus_hardtrim3_rejected`, `library_plus_clock_rejected`, `library_plus_implicon_rejected`, `hardtrim5_range_error_precedes_the_library_refusal` |
| D-b | Single-end reporting line (both code reviewers) | DONE | `test_text_report_omits_r2_clipping_on_single_end`, `flag_summary_omits_read_2_for_single_end`, `flag_summary_lists_all_four_for_paired`, `both_reports_record_the_preset_and_its_expanded_values`, `explicit_clip_flag_overrides_the_preset_and_is_logged` (stderr half), `paired_report_lists_all_four_values` |
| D-c | T-1 — JSON key permutation invisible | DONE | `test_json_report_records_all_four_overrides` (four distinct values 1/2/3/4) |
| D-d | T-2 — override flag labels unasserted | DONE | `test_json_report_records_all_four_overrides` (flag + preset_value per record) and `preset_with_multiple_overrides_round_trips_through_json` |
| D-e | T-3 — `--rename` + uBAM guard not exercised via a preset | DONE | `rename_into_ubam_refused_for_fastq_input_via_library_preset` |
| D-f | F-5a — `is_null()` cannot distinguish absent from present-and-null | DONE | `test_json_report_library_is_null_without_a_preset` and `a_run_without_a_preset_reports_no_library_block`, both with `contains_key` |
| D-g | F-5b — JSON block vs. top level could silently disagree | DONE | the agreement loop in `both_reports_record_the_preset_and_its_expanded_values` |

## E. Assumptions 1–8 — does the code contradict any?

| # | Assumption | Status | Verification |
|---|---|---|---|
| E1 | #441 merges first and unmodified | UPHELD | Parent `872c4ae` = `feat(cli): --library presets expand into end clipping (#440) (#441)`, author `Benjamin Demaille <41732816+BenjaminDEMAILLE@…>` preserved |
| E2 | The four specialty modes do not honour the clip flags | UPHELD | `grep -cE 'clip_r1|clip_r2|three_prime_clip|effective_clips' src/specialty.rs` → **0**. The refusal wording is therefore accurate, and "does not clip" would indeed have been false |
| E3 | Text report and JSON may differ on single-end; a consumer uses `json["mode"]`, and there is no `paired` key | UPHELD | `"mode"` emitted at `report.rs:810`; live SE run yields `mode: single-end`. Grep finds no `paired` JSON key. The F-5b agreement holds on the same run, which is the constraint the assumption rests on |
| E4 | Override log lines stay ungated | UPHELD | `main.rs` override loop untouched; live SE run still prints the override sentence |
| E5 | Four distinct override values are acceptable test data, named as such in the doc comment | UPHELD | `report.rs:1489-1490`: "Four distinct values so that any permutation of the four keys is visible; no real preset has four distinct values, so they come from overrides" |
| E6 | `PresetClips`'s fields stay plain `usize` | UPHELD | All four fields are `pub …: usize`. The step-13 invariant's precondition holds |
| E7 | Docs changes limited to what the change falsifies plus one factual bullet | UPHELD | Three doc hunks only: the two falsified transcripts + the falsified prose, the one new bullet, and the `flags.md` clause (step 23). Nothing else touched |
| E8 | `--clump_only` short-circuits before the new block | UPHELD | Live: `--clump_only --hardtrim5 30 --library emseq` → `--clump_only does not clip; --library is not compatible` (`cli.rs:850`), never reaching `:1111` |

## F. §11 — findings deliberately NOT fixed: was each decision honoured?

Out of scope by decision, so **absence here is conformance, not a gap**. Confirmed genuinely absent
from the diff — none was silently implemented anyway.

| # | Finding | Disposition | Honoured? |
|---|---|---|---|
| F1 | F-5c — no `--library` differential test on the uBAM / interleaved paths | Out | YES — 0 diff hits for `md5` / `assert_ubam_eq` / `byte` in `integration_ubam_out.rs`. The one test added there is a refusal test, not a differential |
| F2 | C-1 — source comments justify the change | Partly done | YES, as scoped — step 1 removed the two `flag_summary` doc-comment sentences (A2); nothing else was rewritten. Spot-checked the new comments: all state current facts, none narrates the prior state or pre-empts an objection |
| F3 | L-1 — text report records clip values only when a preset was used | Out | YES — `write_report_header`'s clip output is still the single preset-gated `Clipping in force:` line; the pre-existing asymmetry is unchanged |
| F4 | F-6 — preset names are case-sensitive | Out | YES — `grep -n ignore_case src/cli.rs` → no hits |
| F5 | F-9 — preset log line repeats per input file | Out | YES — the `main.rs` diff touches only the `flag_summary` argument; the log line's placement is unchanged |
| F6 | F-7 — `json_null()` exists but `"library": null` is hand-written | Out | YES — `json_null` still at `report.rs:1208`; `writeln!(w, "{}\"library\": null,", i2)` still hand-written at `report.rs:905` |
| F7 | L-2 — `match … and_then` vs. the file's let-chain idiom | Out (style) | YES — 0 `and_then` hits anywhere in the diff |
| F8 | L-4 — a failing `Vec<u8>` comparison dumps ~39.7 MB | Out, noted | YES — 0 byte-level comparisons among the new integration assertions; they compare report text as the plan said they would |

## G. Deviations D1–D4 — is each description accurate against the code?

| # | Deviation | Status | Accuracy check |
|---|---|---|---|
| G1 | D1 — rustfmt reflowed three places the plan did not anticipate: the four `bail!` strings, one test array in `cli.rs`, one `assert_eq!` in `library.rs` | DEVIATED (documented) | Accurate in class and consequence (no behaviour change); `cargo fmt --all -- --check` is clean. One imprecision: **two** `cli.rs` test arrays are multi-line, not one — `library_plus_clock_rejected` and `library_plus_implicon_rejected` (the 7-element argv lists), while the two 6-element `hardtrim` ones stayed single-line. Cosmetic; nothing follows from it |
| G2 | D2 — step 23 was a "decide" item, implemented per the plan's own recommendation | DEVIATED (documented) | Accurate. `flags.md:14` carries the half-sentence the plan recommended. Substantively this is conformance rather than deviation (see A32) |
| G3 | D3 — V3 run as an equivalent mutation (suppress the range `bail!`) rather than relocating the refusal block | DEVIATED (documented) | Accurate, and the equivalence claim holds: both mutations test "the range error must be reported first". Independently re-verified by test and by the binary (C3) |
| G4 | D4 — V16 required excluding `command_line` from the JSON comparison | DEVIATED (documented) | Plausible and self-consistent: the two differential binaries live at different paths, and `command_line` records the invocation. Not independently reproduced (a base-commit build is unaffordable at 3.6 GB free). The stronger, cheaper argument was re-verified instead — an `--rrbs` run emits no preset lines at all (C17) |

## H. Plan claims outside the numbered steps

| # | Claim | Source | Status | Notes |
|---|---|---|---|---|
| H1 | Signature is `pub fn flag_summary(&self, paired: bool) -> String` with the specified doc comment | §Signature | DONE | `src/library.rs:118-120`, both doc lines verbatim |
| H2 | All three report sites carry `paired` | §Integration | DONE | `main.rs:1459` (`paired: false`), `:2249` (`paired: false`), `:2628` (`paired: true`), each with `library: Some(cli.effective_clips())` at `:1471` / `:2261` / `:2640`. The plan's `1449/2239/2618` anchors were quoted pre-rebase and have shifted +10 uniformly; the substance holds |
| H3 | JSON untouched | §Integration | DONE | The `report.rs` diff changes one `writeln!` in the **text** header plus tests. No JSON key added, removed, renamed, or re-typed |
| H4 | The refusal precedes `ensure_output_dir` and every dispatch arm | §Integration | DONE | Confirmed behaviourally by C4 (0 files written on all four refusals) |
| H5 | No hot path; four `Option::is_some()` checks short-circuited on `library.is_some()` | §Efficiency | DONE | Implementation is exactly that shape |
| H6 | O-1 — JSON does **not** omit R2 on single-end | §Questions | DONE | Live SE run keeps all four (B13) |
| H7 | O-2 — `--clump_only` keeps "does not clip", specialty modes say "does not honour the clip flags"; no change | §Questions | DONE | `cli.rs:850` vs. `cli.rs:1115-1131`; both messages confirmed live, wording as decided |
| H8 | O-3 — excluded editorial items stay plan notes, not issues | §Questions | DONE | §11 carries all eight; no new issue is referenced in the diff |
| H9 | R-4 — CHANGELOG three-way convergence resolved without dropping entries | §Self-Review | DONE | Single hunk inside #441's entry, no new bullet, unique issue refs 70 → 70 |

## Test verification

All 12 new tests, plus the 4 amended ones the plan required. Every row was executed in this audit.

| # | Test name | File | Status |
|---|---|---|---|
| 1 | `library_plus_hardtrim5_rejected` | `src/cli.rs:2418` | PASS |
| 2 | `library_plus_hardtrim3_rejected` | `src/cli.rs:2428` | PASS |
| 3 | `library_plus_clock_rejected` | `src/cli.rs:2438` | PASS |
| 4 | `library_plus_implicon_rejected` | `src/cli.rs:2458` | PASS |
| 5 | `hardtrim5_range_error_precedes_the_library_refusal` | `src/cli.rs:2478` | PASS |
| 6 | `flag_summary_lists_all_four_for_paired` | `src/library.rs:287` | PASS |
| 7 | `flag_summary_omits_read_2_for_single_end` | `src/library.rs:296` | PASS |
| 8 | `test_text_report_omits_r2_clipping_on_single_end` | `src/report.rs:1404` | PASS |
| 9 | `test_json_report_records_all_four_overrides` | `src/report.rs:1491` | PASS |
| 10 | `paired_report_lists_all_four_values` | `tests/integration_library.rs:258` | PASS |
| 11 | `preset_with_multiple_overrides_round_trips_through_json` | `tests/integration_library.rs:287` | PASS |
| 12 | `rename_into_ubam_refused_for_fastq_input_via_library_preset` | `tests/integration_ubam_out.rs:1044` | PASS |

Amended tests (not counted in the 12):

| Test name | File | Amendment required by the plan | Status |
|---|---|---|---|
| `test_json_report_library_is_null_without_a_preset` | `src/report.rs:1538` | Step 10 — `contains_key` presence assertion | PASS |
| `both_reports_record_the_preset_and_its_expanded_values` | `tests/integration_library.rs:180` | Steps 12 + 13 — SE expectation and the F-5b agreement loop | PASS |
| `explicit_clip_flag_overrides_the_preset_and_is_logged` | `tests/integration_library.rs:130` | Step 14 — stderr negative assertion (C-2 fix) | PASS |
| `a_run_without_a_preset_reports_no_library_block` | `tests/integration_library.rs:227` | Step 15 — `contains_key` presence assertion | PASS |
| `test_text_report_records_preset_name_and_expanded_values` | `src/report.rs:1380` | Step 9 — must stay green **unchanged** | PASS (0 diff hits) |
| `test_json_report_records_library_block` | `src/report.rs:1456` | Step 6 — must be kept **as it is** | PASS (0 diff hits) |

## Gaps (detail)

None. No item is PARTIAL or MISSING.

## Verdict

**Verdict:** COMPLETE

All 23 implementation-outline steps, both behaviour clauses with all 13 of their stated edge cases, all
17 validation rows, all seven target findings, and all eight assumptions are satisfied. The four
deviations are documented in §12 and each description checks out against the code (one cosmetic
imprecision in D1's count, noted at G1). All eight §11 dispositions were honoured — none was silently
implemented.

Three things a reader should carry forward, none of them a gap:

1. **Clippy was not re-run** (3.6 GB free on a 100%-full volume; `--all-targets --release` would rebuild
   16 test binaries). Tests and `cargo fmt --all -- --check` were both run and are clean. Worth one
   `cargo clippy --all-targets --release -- -D warnings` before the PR, once disk allows.
2. **V5's command in the plan cannot work** — `cargo test A B` passes two filters and matches nothing,
   the same trap §12's iteration log records at item 2. Both named tests exist and pass when run
   separately. If the plan text is reused, split that command.
3. **The mutation rows (V2, V6, V6a, V7, V8–V13) remain author-reported.** Reproducing them requires
   editing `impl440`, which a read-only audit may not do. Their target tests were each confirmed to
   exist and pass, and every property with a read-only route (V3, V4, V15, V16) was independently
   re-verified here.

