# Plan Coverage Report

**Mode:** B (code vs. plan — `PLAN.md` v2 carries the numbered outline and the validation matrix, so it is the spec)
**Plan(s):** `plans/08172026_clip-r2-se-warning/PLAN.md` (v2, 2026-08-17)
**Code:** `3161b8a541c09d76d77e4e140f7b2638b8333d5c`, base `fdade7efb8b579fd9b25899d04b1609dba860b20`, worktree `/private/tmp/claude-501/fu440`
**Date:** 2026-08-17
**Verdict:** COMPLETE

## Summary

- Total items: 65
- DONE: 60
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 5 (all four documented deviations D1–D4, plus the one test D3 narrowed)

**Measured, not assumed:**

| Measurement | Result |
|---|---|
| `cargo test --no-fail-fast` | **780 passed, 0 failed, 0 ignored, 18 test binaries** — matches §12 exactly |
| New test binary | `tests/integration_inert_clip_flags.rs` present and run |
| `cargo fmt --all -- --check` | clean (no diff) |
| `cargo clippy --all-targets -- -D warnings` (debug) | exit 0, zero `warning`/`error` lines |
| Diff shape | 8 files, +437 / −18 — matches §12 |

Baseline arithmetic checks out: 768 (base) + 12 new = 780.

---

## Coverage ledger

### A. Implementation outline (§1–§11, 13 numbered items + 5 section-level changes)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `ClipOverride` gains `read_2: bool` with the specified doc comment | §1 | DONE | `library.rs:92-93`, comment verbatim from the plan |
| 2 | `pick` gains a `read_2` param; four call sites pass `false, true, false, true` | §1 | DONE | `library.rs:169-200`. Call order `clip_r1(false)`, `clip_r2(true)`, `three_prime_clip_r1(false)`, `three_prime_clip_r2(true)` — same order and values as `flag_summary`'s array |
| 3 | `main.rs` override loop gated on `cli.paired \|\| !o.read_2` | §2 | DONE | `main.rs:1128`, exactly the plan's filter expression |
| 4 | `report.rs` override loop gated on `config.paired \|\| !o.read_2` | §2 | DONE | `report.rs:284-288`. Nothing else in either loop changed |
| 5 | The reason ladder inserted after the last `bail!` in `validate()` | §3 | DONE | `cli.rs:1138-1174`, before the Perl-era deprecation block. Ladder order is specialty-before-`!paired` as specified; one-line comment states the fact |
| 6 | Test 1 `hardtrim_conflict_precedes_the_library_refusal` | §4.1 | DONE | `cli.rs:2531`. Asserts `"cannot be combined in one invocation"` **and** `!err.contains("--library")` |
| 7 | Test 2 `clock_arity_error_precedes_the_library_refusal` | §4.2 | DONE | `cli.rs:2551`. Asserts `"requires an even number of input files"` + the `--library` absence. Args reflowed by rustfmt (D2) |
| 8 | Test 3 `missing_input_precedes_the_library_refusal` | §4.3 | DONE | `cli.rs:2569`. Asserts `"Input file not found"` + the `--library` absence |
| 9 | All three ordering tests sit beside `hardtrim5_range_error_precedes_the_library_refusal` | §4 | DONE | Sibling is at `cli.rs:2517`; the three new tests follow at 2531 / 2551 / 2569 |
| 10 | Test 4 `inert_clip_flag_is_a_warning_not_an_error` | §5 | DONE | `cli.rs:2589`. Both specified invocations (`--clip_R2 5 R1`, `--hardtrim5 20 --clip_R1 3 R1`) asserted `Ok` |
| 11 | Test 5 `every_override_knows_which_read_it_belongs_to` | §6 | DONE | `library.rs:297`. Resolves `Accel` with all four overridden, asserts `len() == 4`, then `o.read_2 == o.flag.ends_with("_R2")` per override. Flag name used as oracle only in the test; production never inspects it |
| 12 | New file `tests/integration_inert_clip_flags.rs`, modelled on `integration_adapter2.rs` | §7 | DONE | 186 lines, 4 tests, module doc explains the two inert sets and why a bare `!paired` predicate is wrong |
| 13 | Every invocation in the new file passes `-o <tempdir>` | §7 | DONE | All four tests build a `tempdir(tag)` and pass `-o`; nothing writes into tracked `test_files/` |
| 14 | Test 6 `inert_clip_flag_warns_with_a_mode_accurate_reason` — 5-row mode table | §7.6 | DONE | Rows `hardtrim5`, `hardtrim3`, `clock`, `implicon`, `single_end`. Each asserts the run succeeds, stderr contains `"is not used in this mode"`, **and** that row's reason substring (failure message *"warning must give a reason true for this mode"*, mirroring `-a2`). Adds a fourth assertion that the flag is named — stronger than specified |
| 15 | Specialty rows pass two files and no `--paired`, and use `--clip_R1` | §7.6 | DONE | `hardtrim5`/`hardtrim3` use `BS-seq_10K_R{1,2}`; `clock`/`implicon` use `clock_10K_R{1,2}` as the plan directed. All four specialty rows use `--clip_R1`, which proves the all-four candidate set |
| 16 | Test 7 `paired_clip_flags_do_not_warn` | §7.7 | DONE | `--paired --clip_R2 5 --three_prime_clip_R2 7`: asserts `!err.contains(INERT)` |
| 17 | Test 8 `preset_alone_on_single_end_does_not_warn` | §7.8 | DONE | `--library emseq R1`: no `"is not used in this mode"`. The V4/V6 regression guard |
| 18 | Test 9 `each_inert_flag_gets_its_own_line` | §7.9 | DONE | Asserts the warning-line count is **exactly 2** and that each flag+value (`--clip_R2 5`, `--three_prime_clip_R2 7`) is named |
| 19 | Test 10 `single_end_announces_no_read_2_clip` | §8.10 | DEVIATED | Documented as **D3**. stderr assertion targets the override *sentence* (`"--clip_R2 5 was given on the command line"`) rather than the bare flag name, plus a positive assertion that the warning **is** present. The report half is **not** narrowed — still `!text.contains("--clip_R2")`. See §F |
| 20 | Test 11 `paired_still_announces_r2_overrides` | §8.11 | DONE | Asserts the override line is present in **both** stderr and the text report — guards an over-eager gate |
| 21 | Test 12 `test_text_report_omits_r2_override_lines_on_single_end` | §9.12 | DONE | `report.rs:1430`. `paired = false`, `emseq` with `clip_r2: Some(5)`; also sets `clip_r1: Some(3)` and asserts the R1 override **is** still named — stronger than specified |
| 22 | Docs: new bullet in the cross-flag-interactions list | §10 (item 13) | DONE | `docs/src/content/docs/guide/flags.md:15`, immediately after the `--library` bullet. Same shape as the `--fastqc` bullet at `:19` (bolded flag names, the behaviour, then the reason). Names all four flags, both inert sets, and the report consequence |
| 23 | CHANGELOG: one clause folded into #440's existing entry, **no second bullet** | §11 | DONE | Single entry still; clause added to the paragraph describing the refusal. Wording is a slight **superset** of the plan's draft — it also carries the B2 half ("with nothing on either surface announcing a Read 2 clip the run did not perform"), which the plan's draft clause omitted |

### B. Behaviour clauses

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 24 | Ladder computed after every `bail!`, specialty modes before `!paired` | B1.1 | DONE | Verified structurally: the last `bail!`/`return Err` in `validate()` is at `cli.rs:1132`; `awk` over 1139→end-of-`validate()` finds **zero** further `bail!`/`return Err`. Ladder order matches `-a2`'s |
| 25 | The reason determines the candidate flag set (all four vs. the two R2) | B1.2 | DONE | `all_four` bool gates `clip_r1`/`three_prime_clip_r1` to `None` on the `!paired` arm |
| 26 | One warning per candidate flag the user actually set, then continue; never an error | B1.3 | DONE | Live: `--hardtrim5 30 --clip_R1 3` exits **0** with 1 warning. Test 4 pins `Ok` at the `validate()` layer |
| 27 | Paired runs, and runs with no candidate flag set, emit nothing | B1.4 | DONE | Live-verified both halves: single-end with only `--clip_R1 3` → **0** warnings; `--paired` with no clip flags → **0** warnings. Test 7 pins the paired-with-R2-flags case |
| 28 | Reason wording reuses *does not honour the clip flags*, not "performs no clipping" | B1 | DONE | All three specialty arms use the shipped `--library` vocabulary. Justified: `specialty.rs:107` and `:219` do call `clip_5prime`, so "performs no clipping" would be false — both cited lines verified present |
| 29 | Edge case: `--clip_R2 0` warns | B1 | DONE | **Live-verified.** `--clip_R2 0` on single-end emits `WARNING: --clip_R2 0 was given but is not used in this mode (…single-end run). Ignoring.` The `if let Some(n)` form admits `Some(0)`. No dedicated test — the plan specified none |
| 30 | Edge case: multi-file single-end warns once per flag per invocation, not per file | B1 | DONE | **Live-verified.** Two input files × two R2 flags → exactly **2** warning lines, not 4. Structural: the ladder lives in `validate()`, which runs once |
| 31 | Edge case: `--paired` with one interleaved uBAM stays silent | B1 | DONE | By construction: `self.paired` is true, so the ladder falls through every arm to `None`. No candidate set is computed, so nothing can warn |
| 32 | Edge case: a preset's values are invisible here (raw fields) | B1 | DONE | **Live-verified independently of the test.** `--library emseq --clip_R2 5` on single-end warns about `--clip_R2` **only** — even though `emseq` also sets `three_prime_clip_r2: 10`. That silent second value is direct proof the raw fields are read, not `effective_clips()` |
| 33 | `ClipOverride` gains `read_2: bool` | B2.1 | DONE | = item 1 |
| 34 | Both override-printing loops skip `read_2 && !paired` | B2.2 | DONE | = items 3, 4. Both surfaces named by the plan (`main.rs:1128`, `report.rs:284`) are the ones gated |
| 35 | Paired output unchanged; single-end loses only R2 override lines | B2.3 | DONE | Item 20 pins the paired half; V18's differential (§12) measured it — one line differs, on one shape |
| 36 | The invariant: no surface announces a Read 2 clip on a single-end run | B2 | DONE | **Live-verified.** Single-end `--library emseq --clip_R2 5`: stderr banner reads `Library preset 'emseq' selected: --clip_R1 10 --three_prime_clip_R1 10` with no R2 override line; the text report contains **zero** occurrences of `--clip_R2` and `Clipping in force: --clip_R1 10 --three_prime_clip_R1 10` |
| 37 | The report-bytes change is acknowledged, not hidden | B2.4 | DONE | §Integration states it plainly; §12's V18 quantifies it as the single differing line |
| 38 | Three `Cli::validate()` ordering tests, each also asserting `--library` is unmentioned | B3 | DONE | = items 6, 7, 8. All three carry `!err.contains("--library")` |

### C. Validation matrix V1–V18

Every row has evidence. The 8 mutation rows are recorded in §12's matrix; V16 and V18 have their own §12 subsections; the remaining rows are pass/fail assertions I re-ran directly.

| # | V | Claim | Evidence | Status |
|---|---|---|---|---|
| 39 | V1 | Reason true in every mode — 5 rows pass | Re-ran: `inert_clip_flag_warns_with_a_mode_accurate_reason ... ok` (all 5 rows in one test) | DONE |
| 40 | V2 | Can fail: ladder loses its three specialty arms | §12 matrix — the table test goes red | DONE |
| 41 | V3 | All four flags named under the specialty modes | Re-verified in source: all four specialty rows pass `--clip_R1`; test passes | DONE |
| 42 | V4 | Can fail: `all_four` hard-coded `false` | §12 matrix — same table red | DONE |
| 43 | V5 | No noise on a preset-only single-end run | Re-ran: `preset_alone_on_single_end_does_not_warn ... ok` | DONE |
| 44 | V6 | Can fail: read `effective_clips()` instead of raw | §12 matrix — **2** red (`preset_alone_on_single_end_does_not_warn` + `explicit_clip_flag_overrides_the_preset_and_is_logged`). Independently corroborated by item 32's live check | DONE |
| 45 | V7 | No surface announces an R2 clip on single-end | Re-ran both: `single_end_announces_no_read_2_clip ... ok`, `test_text_report_omits_r2_override_lines_on_single_end ... ok`. Plus item 36's live check | DONE |
| 46 | V8 | Can fail: revert either gate **separately** | §12 matrix splits it: **V8a** (`main.rs` only) → `single_end_announces_no_read_2_clip`; **V8b** (`report.rs` only) → the report unit test + the integration test. Each surface independently covered, so neither gate can be reverted silently | DONE |
| 47 | V9 | The paired half survives | Re-ran: `paired_still_announces_r2_overrides ... ok` | DONE |
| 48 | V10 | Can fail: gate drops R2 unconditionally | §12 matrix — `paired_still_announces_r2_overrides` red | DONE |
| 49 | V11 | `read_2` matches the flag name | Re-ran: `every_override_knows_which_read_it_belongs_to ... ok` | DONE |
| 50 | V12 | Can fail: swap two `read_2` args at `pick`'s call sites | §12 matrix — **6** red, including the oracle test | DONE |
| 51 | V13 | One line per flag — exactly 2 | Re-ran: `each_inert_flag_gets_its_own_line ... ok`; the test asserts `count() == 2`. Corroborated live in item 30 | DONE |
| 52 | V14 | The three precedence relations hold | Re-ran all three ordering tests — ok | DONE |
| 53 | V15 | Each can fail | §12 matrix — **9** red (3 new + 6 pre-existing sharing those checks). Executed as an equivalent (D4) | DONE |
| 54 | V16 | The stderr/text/JSON triple | **Independently re-measured.** stderr warns and names `--clip_R2 5`; text report has 0 mentions of `--clip_R2`; JSON has `mode = single-end` and `parameters.clip_r2 = 5`. All four hold, exactly as §12 records | DONE |
| 55 | V17 | Suite and lints | **Re-measured:** 780 passed / 0 failed / 18 binaries; `cargo fmt --all -- --check` clean; `cargo clippy --all-targets -- -D warnings` exit 0 with zero warning lines. Release-profile clippy is evidenced by §12's explicit claim; I did not re-run it, because the brief forbade release builds on the 100 %-full volume (5.7 GiB free at audit end) | DONE |
| 56 | V18 | Paired output byte-identical | §12: 5 shapes, 24 files, none empty; every FASTQ and `.txt` identical on all paired shapes and on single-end without an R2 override; all 8 JSON identical once `command_line` (which records the two binary paths) is excluded. The sole behavioural change is one line, quoted in §12. Not independently re-runnable here — it needs a base-commit build, which the disk constraint rules out | DONE |

### D. Assumptions 1–6

| # | Assumption | Status | Notes |
|---|---|---|---|
| 57 | 1 — warning reads raw `self.clip_r*`, never `effective_clips()` | DONE | Source confirms; independently proven live by item 32 (the preset's own `three_prime_clip_r2: 10` produces no warning) |
| 58 | 2 — `--clip_R2 0` warrants a warning | DONE | Live-verified, item 29 |
| 59 | 3 — one line per flag | DONE | Item 18 / V13 |
| 60 | 4 — placement after every `bail!` keeps a `validate()`-level refusal quiet | DONE | Zero `bail!`/`return Err` between the insertion point and the end of `validate()`; last one is at `cli.rs:1132` |
| 61 | 5 — the specialty modes honour **none** of the four clip flags | DONE | `grep` for `clip_r1\|clip_r2\|three_prime_clip\|effective_clips` over `src/specialty.rs` returns **nothing**. The mode's own fixed clipping is the two `clip_5prime` calls at `:107` and `:219`, exactly as cited |
| 62 | 6 — `read_2` set at `pick`'s call sites, not derived from the flag string; item 11's test keeps the two lists in agreement | DONE | **Checked against the actual test, as requested.** Production code contains no `ends_with`/string inspection of `flag`; the only `flag.ends_with("_R2")` in the crate is inside the test, where it is the oracle. The test first asserts `overrides.len() == 4`, so a dropped override fails too, not just a mis-set bool. V12's 6-red result confirms the guard bites |

### E. The three findings this work closes

| # | Finding | Mapped tests | Status |
|---|---|---|---|
| 63 | **A's F-1** — a single-end run announced a Read 2 clip it did not perform | `single_end_announces_no_read_2_clip` (item 19) + `test_text_report_omits_r2_override_lines_on_single_end` (item 21), one per surface. V8a/V8b prove each independently | DONE |
| 64 | **A's F-2** — the refusal's precedence was unpinned | The three ordering tests (items 6–8), each also asserting `--library` is unmentioned. V15 red on all three | DONE |
| 65 | **B's F-5, user-visible half** — all four clip flags are inert under the specialty modes, not just the R2 pair | The table's four specialty rows use `--clip_R1` (item 15), which can only warn because of the mode. V4 (`all_four = false`) reddens exactly those rows | DONE |

### F. Deviations D1–D4 — descriptions checked against the code

All four are DEVIATED, documented. Each description is accurate:

| D | Plan's description | Accurate? |
|---|---|---|
| **D1** | `cargo build` failed with `error: cached cgu … should have an object file, but doesn't`; a corrupted incremental cache from earlier disk exhaustion, fixed by `cargo clean` | **Accurate as a process note.** Build-environment only, no code consequence. Nothing in the tree reflects it, and nothing should |
| **D2** | rustfmt reflowed one test's argument array (`clock_arity_error_precedes_the_library_refusal`) | **Accurate.** That test's `parse_from` array is the one broken one-arg-per-line in the diff. `cargo fmt --all -- --check` is clean, so the reflow was applied, not left pending. No behaviour change |
| **D3** | Step 10's assertion narrowed: the plan said stderr contains no `--clip_R2`, but the warning names the flag on purpose, so the test targets the override *sentence* and separately asserts the warning is present | **Accurate, and the narrowing is exactly what §12 says.** Code asserts `!err.contains("--clip_R2 5 was given on the command line")` (the override sentence, not the flag name) and `err.contains("is not used in this mode")`. Critically, the narrowing is **stderr-only** — the report half is still the full `!text.contains("--clip_R2")`, so the invariant is unweakened where it can be held. The plan's original assertion would indeed have contradicted the feature: I confirmed live that stderr must contain `--clip_R2 5` for the warning to do its job |
| **D4** | V15 executed as an equivalent (suppress the three earlier checks) rather than by relocating the refusal block | **Accurate**, and the property proved is the same one V15 names: each earlier check must report first. 9 tests red, which is a broader blast radius than relocation would have produced |

### G. Explicitly deferred decisions — honoured as written, not silently implemented

| Item | Status |
|---|---|
| **Q-1** (documented, decided) — document the clip warning even though `-a2`'s stays undocumented | Honoured. The new `flags.md` bullet documents the clip case; `-a2`'s warning remains undocumented. No `-a2` docs were added, so the recorded asymmetry is what shipped |
| **Q-2** (open, deferred) — the JSON still records an inapplicable `clip_r2` | **Honoured, verified live.** Single-end `--library emseq --clip_R2 5` still writes `parameters.clip_r2 = 5` with `mode = single-end`. No JSON change was slipped in: `report.rs`'s diff touches only the text header loop and adds a test. The new `read_2` field does **not** appear in the serialized `overrides` array (`{"flag": "--clip_R2", "preset_value": 10, "user_value": 5}`), so the schema is untouched for MultiQC consumers |
| **R-3** — the specialty modes warn but are **not** refused | **Honoured, verified live.** `--hardtrim5 30 --clip_R1 3` exits **0** and warns. Test 4 pins `validate()` returning `Ok` for both the single-end and the `--hardtrim5` shape. No refusal was added |

---

## Gaps (detail)

None. No PARTIAL and no MISSING items.

## Test verification

All 12 new tests exist and pass. Names, files and assertions match what the outline specified; where they differ they assert **more**, not less. Verified from a full `cargo test --no-fail-fast` run (780 passed / 0 failed / 18 binaries).

| # | Test name | File | Status |
|---|-----------|------|--------|
| 1 | `hardtrim_conflict_precedes_the_library_refusal` | `src/cli.rs:2531` | PASS |
| 2 | `clock_arity_error_precedes_the_library_refusal` | `src/cli.rs:2551` | PASS |
| 3 | `missing_input_precedes_the_library_refusal` | `src/cli.rs:2569` | PASS |
| 4 | `inert_clip_flag_is_a_warning_not_an_error` | `src/cli.rs:2589` | PASS |
| 5 | `every_override_knows_which_read_it_belongs_to` | `src/library.rs:297` | PASS |
| 6 | `inert_clip_flag_warns_with_a_mode_accurate_reason` | `tests/integration_inert_clip_flags.rs` | PASS (5 rows) |
| 7 | `paired_clip_flags_do_not_warn` | `tests/integration_inert_clip_flags.rs` | PASS |
| 8 | `preset_alone_on_single_end_does_not_warn` | `tests/integration_inert_clip_flags.rs` | PASS |
| 9 | `each_inert_flag_gets_its_own_line` | `tests/integration_inert_clip_flags.rs` | PASS |
| 10 | `single_end_announces_no_read_2_clip` | `tests/integration_library.rs` | PASS (assertion narrowed per D3) |
| 11 | `paired_still_announces_r2_overrides` | `tests/integration_library.rs` | PASS |
| 12 | `test_text_report_omits_r2_override_lines_on_single_end` | `src/report.rs:1430` | PASS |

Tests asserting **more** than the outline required: #6 adds a "the warning must name the flag" assertion; #11 covers both surfaces where the plan named the override line generally; #12 adds a positive assertion that the Read 1 override is still printed; #5 asserts `overrides.len() == 4` before the per-override oracle check.

## Verdict

**Verdict:** COMPLETE

Every numbered outline item (1–13 across §1–§11), including the docs step (§10) and the CHANGELOG step (§11), landed. All three behaviour clauses and all eight of their edge cases are implemented, and the four I could exercise from outside the test suite (`--clip_R2 0`, multi-file single-end, the preset-invisibility property, and the no-surface-announces-R2 invariant) behave as specified against the release binary. All 18 validation rows carry evidence — the 8 mutation rows in §12's matrix, V16 and V18 in their own §12 subsections, and the 8 pass/fail rows re-run here. No assumption is contradicted by the code; Assumption 6's two-list guard is real (production contains no flag-string inspection, and the oracle lives only in the test). The three findings each map to implemented tests. D1–D4 are documented and, checked line by line against the code, accurately described — including D3, whose narrowing is stderr-only and leaves the report-side assertion at full strength. Q-2 and R-3 were honoured as deferred, and neither was silently implemented.

### Non-blocking observations (not gaps, no action implied)

1. **§12's "23 outline steps" is not derivable from the outline as written.** The outline has 11 sections carrying 13 numbered items. 11 + 12 new tests = 23, which is presumably the arithmetic. Every step landed either way; only the label is ambiguous.
2. **The CHANGELOG edit leaves one line longer than its neighbours** (`… and both numbers are named in the log; a preset that supplies --clip_R2 suppresses the --rrbs directional`, ~110 chars against a ~95-char surrounding wrap). Cosmetic reflow artefact; no content problem.
3. **The CHANGELOG clause is a superset of the plan's draft**, adding the B2 gating half. It stays one entry with no second bullet, which is what §11 required.
4. **Release-profile clippy was not re-run** (only debug). §12 claims it clean; the brief forbade release builds on the full volume. This is the one V-row assertion resting solely on the implementer's record rather than on a measurement I repeated.
