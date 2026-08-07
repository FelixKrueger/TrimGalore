# Plan Coverage Report

**Mode:** B (code vs. implementation plan — §5's 13 steps are the implementation plan; no separate `IMPL.md` exists by design)
**Plan(s):** `plans/08062026_se-output-collision-preflight/PLAN.md` (v2)
**Date:** 2026-08-07
**Verdict:** COMPLETE

Audited state: branch `fix/383-output-collision-preflight` off `dev` @ `72624c7`, change uncommitted
(7 modified files + new untracked `tests/integration_output_collision.rs`).

`cargo test` from the crate root: **521 passed, 0 failed** (399 lib + 122 integration across 12
integration binaries; `tests/integration_output_collision.rs` contributes 20).

## Summary

- Total items: 25
- DONE: 25
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 0 (D1–D6 are all documented in §13 and all verified against the code)

Breakdown: 13 §5 steps · 7 §10 validation groups (V1 9/9, V2 11/11, V3 7/7, V4, V5, V6 3/3, V7) ·
2 call-site items · 2 items for §3.5 and §4 · 1 item for §13's deviations.

## Coverage ledger

### §5 Implementation outline — steps 1–13

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `collision_key` + `preflight_output_collisions` in `io.rs`; extend `norm_path` doc-comment | Step 1 | DONE | `io.rs:48` (private) and `:57` (pub), placed directly after `norm_path` as specified. `norm_path`'s doc-comment now names `collision_key` (`:37-38`); the stale "`main::run` paired-end" locator is gone. Input keys held in a `HashMap<String,&PathBuf>` rather than §4.1's `HashSet<String>` — forced by §3.2's requirement that the alias message name the input; superset, not a shortfall. |
| 2 | Delete `preflight_collision_bam`; re-point 3 calls; doc-comment rejoins `resolve_clump_layout` | Step 2 | DONE | Function gone (no grep hit anywhere). Baseline `72624c7` had calls at `:533`/`:566`/`:614`; all three now call the shared helper (`:548`, `:585`, `:633`). The clumpify-budget doc-comment is attached to `resolve_clump_layout` at `main.rs:91-95`. |
| 3 | Convert the 4 hand-rolled copies to the shared helper | Step 3 | DONE | All four baseline `Output path collision` string literals (base `:490`, `:718`, `:1812`, `:2420`) are deleted. Converted: `:513` clump_only SE FASTQ (one candidate/input); `:734` paired FASTQ (accumulates across **all** `chunks(2)` into one `Vec` before a single call, keeps `--retain_unpaired` `:709-719` and `--passthrough` `:724-731` candidates); `:1825` paired uBAM (one candidate/pair); `:2430` `run_specialty_paired` (both names per pair, all pairs, one call). Every site passes `guarded_inputs(...)`. §2.5 side effect confirmed: `grep -rn "output-dir" src/ ci.yml` (excluding `output_dir`) returns nothing. |
| 4 | `HardtrimEnd` enum + `pub` on both hardtrim namers; update 4 internal sites + 1 test | Step 4 | DONE | `specialty.rs:19` enum, `:25-29` `as_str()` → `"5prime"`/`"3prime"`. Both namers `pub fn` (`:442`, `:465`). Four internal sites updated: `:45`, `:88`, `:139`, `:191`. `bgz_stem_reaches_specialty_output_names` updated (`:844`, plan anchor `:822` shifted). |
| 5 | `planned_hardtrim_outputs` in `main.rs` | Step 5 | DONE | `main.rs:71-89`; matches §4.3 signature exactly, `match cli.output_format` dispatches FASTQ vs uBAM namer. |
| 6 | Guard `--hardtrim5` with `HardtrimEnd::Five` + `HARDTRIM_HINT` | Step 6 | DONE | `main.rs:359-364`, first statement in the `if let Some(n) = cli.hardtrim5` block, before the `for input in &cli.input` loop at `:365`. `HARDTRIM_HINT` at `:56-59` names both real remedies ("one invocation per input", "distinct basenames"). |
| 7 | Guard `--hardtrim3` with `HardtrimEnd::Three` | Step 7 | DONE | `main.rs:390-395`, identical shape, `HardtrimEnd::Three`. |
| 8 | Guard SE trim FASTQ | Step 8 | DONE | `main.rs:791-799`, in the `} else {` SE arm before the loop at `:800`. A9 verified: the candidate expression at `:796` is argument-identical to the writer's at `:1103`. Uses `guarded_inputs(&cli)` in place of the plan's `&cli.input` (D2). |
| 9 | Guard SE trim uBAM | Step 9 | DONE | `main.rs:1863-1869`. A9 verified: candidate `:1867` argument-identical to writer `:1893`. |
| 10 | `cli.rs` rejects duplicate SE inputs with its own message | Step 10 | DONE | `cli.rs:625-642`, inside `Cli::validate()` (hence before `ensure_output_dir` at `main.rs:314`). Message "Input file … was given more than once (arguments N and M)." plus the why-clause — matches §11 Open 1's proposed shape. Predicate `!self.paired && !self.clock && self.implicon.is_none()` per D1. |
| 11 | Tests: unit (io.rs, cli.rs) + `tests/integration_output_collision.rs` | Step 11 | DONE | 12 unit tests in `io.rs:936-1164`, 3 in `cli.rs:1201-1233`, 20 in the new integration file. Follows `integration_paired_format_guard.rs` conventions (`CARGO_BIN_EXE_trim_galore`, `tempdir(tag)` keyed on `process::id()`, `(bool, String)` from `Command::output()`). Step 11's absolutised-fixture-path clause is satisfied by a stronger route: every fixture is written in-test (file header `:16-18`), so no `current_dir` case ever references a crate-root-relative path. |
| 12 | CI: three validation steps | Step 12 | DONE | `ci.yml:642`, `:658`, `:677` — inserted directly after the `#216` guard block ending at `:640`. Idiom matches the existing guards exactly (`set +e` / `rc=${PIPESTATUS[0]}` / `set -e` / `test $rc -ne 0` / `grep -q "Output path collision"` / `test -z "$(ls -A …)"`). |
| 13 | CHANGELOG: bug fixes + changes | Step 13 | DONE | `#### Bug fixes`: #383 plus both v2 defects (raw-string key, output-vs-output-only), each with its reachability. `#### Changes`: the A6 duplicate-input refusal (flagged as a behaviour change, noting the `--clock`/`--implicon` carve-out) and the `--output-dir` → `--output_dir` message correction. |

### §10 Validation

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| V1 | 9 unit sub-cases in `src/io.rs` | §10 V1 | DONE | **9/9 present**, plus 1 bonus. 1→`preflight_accepts_distinct_outputs`; 2→`preflight_rejects_identical_outputs` (asserts the path is named); 3→`preflight_rejects_case_only_variants`; 4→`preflight_rejects_dot_slash_alias`; 5→`preflight_rejects_absolute_versus_relative_alias` (builds the key off `current_dir()`); 6→`preflight_rejects_output_that_aliases_an_input` (asserts alias wording **and** the absence of duplicate wording); 7→`preflight_accepts_dotdot_alias_known_limitation`; 8→`preflight_appends_hint_only_when_given` (both directions); 9→`preflight_accepts_empty_and_single`. Bonus: `preflight_rejects_input_alias_across_spellings`. |
| V2 | 11 integration rejection sub-cases | §10 V2 | DONE | **11/11 present**, plus 2 bonus, in `tests/integration_output_collision.rs`. 1→`se_trim_rejects_shared_stem`; 2→`se_trim_rejects_dot_slash_alias`; 3→`se_trim_rejects_same_basename_across_dirs_with_output_dir`; 4→`se_trim_ubam_rejects_shared_stem`; 5→`hardtrim5_rejects_same_basename_across_dirs` (also greps the hint); 6→`hardtrim3_rejects_same_basename_across_dirs`; 7→`hardtrim5_ubam_rejects_shared_stem`; 8→`hardtrim3_ubam_rejects_shared_stem`; 9→`se_trim_rejects_output_that_aliases_an_input`; 10→`paired_rejects_output_that_aliases_an_input`; 11→`se_trim_rejects_gz_and_bgz_sharing_a_stem`. Bonus: `se_trim_rejects_absolute_versus_relative_alias`, `duplicate_input_gets_a_precise_message`. The empty-output-directory assertion is centralised in `assert_rejected_cleanly` (`:95-116`) and used by 9 of the 11; the two no-`-o` cases assert exact directory contents instead, which is the stronger form where output lands beside the inputs. Shape note: V2.5/V2.6 add `-o` on top of `current_dir(tempdir)` — the plan named only `current_dir` — for the same reason V6.2 states explicitly (the residue check needs a dedicated directory, and `-o` does not rescue a hardtrim collision). Consequence recorded under §3.5 row 11 below. |
| V3 | 7 integration acceptance sub-cases | §10 V3 | DONE | **7/7 present**, every one asserting content or filename. 1→`se_trim_accepts_same_basename_across_dirs` (four assertions: each output holds 40 of its own reads and 0 of the other's); 2→`se_trim_accepts_single_input` (read count); 3→`hardtrim5_accepts_distinct_stems_and_names_output` (`x.20bp_5prime.fq` + `y.20bp_5prime.fq` by name); 4→`hardtrim3_accepts_distinct_stems_and_names_output` (`3prime` by name); 5→`hardtrim3_ubam_accepts_distinct_stems_and_names_output` (`*.20bp_3prime.bam`); 6→`se_trim_ubam_accepts_distinct_stems`; 7→`paired_accepts_two_distinct_pairs`. The `3prime`/`5prime` filename assertions are what make the §4.2 enum's purpose testable. |
| V4 | `src/cli.rs` duplicate-SE-input unit tests | §10 V4 | DONE | `cli.rs:1201-1233`: `test_validate_single_end_duplicate_input_rejected` (asserts the dedicated wording **and** `!contains("APFS/NTFS")`), `test_validate_hardtrim_duplicate_input_rejected` (specialty-mode coverage), `test_validate_single_end_distinct_inputs_accepted` (negative control). The pre-existing duplicate-pair and R1≠R2 tests still pass — they are what caught D1, and all 399 lib tests are green. |
| V5 | A2-as-a-table unit test | §10 V5 | DONE | `distinct_primary_outputs_imply_distinct_secondary_outputs` (`io.rs:1071-1123`) asserts the plan's stated direction (distinct primaries ⇒ distinct secondaries) over the full 2×2×2 `--basename`/`--dont_gzip`/`-o` matrix across 3 inputs, for `report_name`, `json_report_name`, `clumping_report_name`, and the demux stem via `demux::demux_base_name` (D6). The `--basename`-collapses-everything case is skipped with a comment pointing at the two guards that reject it. FastQC's zip name is documented as not separately assertable (`:1065-1070`) with the `--fastqc_args "-o DIR"` exception named, so the row cannot be "fixed" away. `primary_output_key_is_coarser_than_secondary_keys` retained as the complement, and covers the uBAM primary too. |
| V6 | 3 CI validation steps | §10 V6 | DONE | `ci.yml:642` (SE collision, `illumina_10K.fastq.gz` copied to `sample.fastq.gz` + `sample.fq.gz`), `:658` (`--hardtrim5 30`, two subdirectories, with `-o`, and it greps the mode-specific hint as well as the collision prefix), `:677` (output-aliases-input, which additionally asserts the named input is md5-byte-identical afterwards and that the directory still holds exactly 2 files). All three verified to run and pass locally against the release binary as reproductions (a)/(c)/(f1) below. |
| V7 | Gates + manual reproduction re-runs | §10 V7 | DONE | Verified independently this audit, not taken on report. `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` exit 0, no warnings; `cargo test` from the crate root **521 passed / 0 failed**. `target/release/trim_galore` rebuilt, then all **eight** §2.2 invocations re-run from a scratch directory: each exits 1 with the correct message and writes nothing. The V3.1 positive control accepts with zero crosstalk. Full transcript in the table below. |

### The 11 call sites

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| C | 3 re-pointed + 4 converted + 4 new = 11; no `cli.input` loop left unguarded | §1 / §5 | DONE | Exactly **11** `naming::preflight_output_collisions` calls in `src/main.rs`: `:360`, `:391`, `:513`, `:548`, `:585`, `:633`, `:734`, `:799`, `:1825`, `:1869`, `:2430`. Cross-checked against `72624c7`, which had 3 `preflight_collision_bam` calls (`:533`/`:566`/`:614`) + 4 hand-rolled `Output path collision` literals (`:490`/`:718`/`:1812`/`:2420`) = 7 guarded sites. 7 + 4 new = 11, and the arithmetic 3 re-pointed + 4 converted + 4 new holds exactly. |
| C2 | No fifth hole — every `cli.input` traversal accounted for | §2.1 | DONE | `grep -nE "for .* in .*cli\.input\|cli\.input\.iter\(\)\|cli\.input\.chunks\("` on `src/main.rs` returns 15 loops. **5** are candidate-accumulation loops feeding a guard (`:577`→`:585`, `:630`→`:633`, `:700`→`:734`, `:1817`→`:1825`, `:2425`→`:2430`). The other **10** are output-producing, and each has an immediately preceding guard: `:365`←`:360`, `:396`←`:391`, `:514`←`:513`, `:591`←`:585`, `:634`←`:633`, `:744`←`:734`, `:800`←`:799`, `:1828`←`:1825`, `:1870`←`:1869`, `:2433`←`:2430`. This matches §2.1's "ten `cli.input` loops" independently. The 11th call (`:548`) guards the loop-free single-interleaved-uBAM `--clump_only` branch. The four remaining `cli.input.len() == 1` branches (`:475`, `:536`, `:677`, `:1796`) each produce outputs from a single input, so no duplicate and no self-alias is constructible — §2.1 marks them n/a and that holds. `--clump_only` PE FASTQ is not a separate hole: it delegates to `run_specialty_paired` (`:483`), guarded at `:2430`. |

### §3.5 edge cases and §4 signatures

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| E | 15-row edge-case table traces to test / assumption / out-of-scope | §3.5 | DONE | All 15 rows trace. 1 one input→`preflight_accepts_empty_and_single` + `se_trim_accepts_single_input`. 2 `.fastq.gz`+`.fq.gz`→`se_trim_rejects_shared_stem`, V6.1, repro (a). 3 `.bgz`→`se_trim_rejects_gz_and_bgz_sharing_a_stem`, repro (b). 4 `./x`→2 tests + repro (e2). 5 `/abs/x`→2 tests + repro (e1). 6 `a/../x` accepted→`preflight_accepts_dotdot_alias_known_limitation` + repro, A8 held. 7 output==input→3 tests, V6.3, repros (f1)(f2). 8 same file twice→2 unit + 1 integration + repro, A6. 9 case-only→unit only, correctly, since two case variants cannot coexist on APFS. 10 same basename/diff dirs/no `-o`/SE accepted→`se_trim_accepts_same_basename_across_dirs` (content-attributing) + V3.1 repro. 11 same basename/diff dirs/no `-o`/hardtrim rejected→**not pinned by `cargo test`** (both hardtrim FASTQ rejection tests add `-o`); pinned only by V7 repro (c), which I ran and which rejects with the hint. Low risk, since the writer and the candidate builder share one namer (A9 + the §4.2 enum), so the two arms cannot disagree — but it is the mode's documented `--hardtrim5 N */*.fastq.gz` pattern, and the automated pin is weaker than the plan's wording implies. 12 with `-o`→`se_trim_rejects_same_basename_across_dirs_with_output_dir` (V2.3), the highest-value case, present. 13 `--basename` >1 SE input→`cli.rs:620` untouched, lib tests green. 14 zero inputs→clap `required` (A5). 15 `--hardtrim5`+`--hardtrim3`→unchanged, §11 Open 3, explicitly left undone in §13. |
| S | 4 signature groups match the code | §4.1–§4.4 | DONE | §4.1 both signatures match character-for-character, including `hint: Option<&str>` and the private `collision_key`; the doc-comment enumerates the call sites and adds "A new dispatch path needs a call here too" as §4.1 intended. §4.2 `HardtrimEnd { Five, Three }` with `as_str()` → `"5prime"`/`"3prime"`, both namers `pub fn` taking `end: HardtrimEnd`. §4.3 `planned_hardtrim_outputs` matches, `std::path::PathBuf` spelled in full as §8 requires (`main.rs:5` still imports `Path` only). §4.4 present in `Cli::validate()` with its own message; placed at `:625` (after the `--basename` check) rather than literally adjacent to the duplicate-pair check at `:534-558`, but all four properties §4.4 argued for hold — precise message, runs before `ensure_output_dir` (`main.rs:314`), one place for all non-paired modes, unit-testable. §4.4's "covers every mode in one place" is narrowed by D1's `--clock`/`--implicon` carve-out, which is documented. |

### §13 deviations D1–D6

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| D | Each of D1–D6 is real, matches the code, and is the only deviation | §13 | DONE | **D1** predicate is `!self.paired && !self.clock && self.implicon.is_none()` (`cli.rs:629`) — verified, and the three `--clock`/`--implicon` tests it protects are green. **D2** `guarded_inputs()` at `main.rs:62`; `grep -c "guarded_inputs("` = 12 = 1 definition + 11 call sites, and all 11 pre-flight calls pass it, so §3.2's input set (positionals + `--passthrough`) reaches every guarded path. **D3** both defect-1 integration tests now run **without** `-o` (`:144`, `:162`, each with a comment saying why) and close with `assert_dir_holds_only(&dir, &["sample.fastq", "sample.fq"])`. **D4** `tempdir()` ends in `std::fs::canonicalize` (`:35`) with the `/tmp` → `/private/tmp` rationale spelled out. **D5** confirmed and confirmed *consistent*: the pre-existing `#216` step at `ci.yml:617` has the identical unquoted-`#` quirk, so all four collision steps truncate the same way. **D6** `demux::demux_base_name` at `demux.rs:122`, called from both `demultiplex` (`:153`) and the V5 test; the separate `trimmed_name` binding that feeds the Perl-byte-identity-relevant "Processed sequences from file" message survives at `:148-152`. |

## Gaps (detail)

No coverage gaps. Every §5 step, every §10 validation sub-case, all 11 call sites, all 15 §3.5
edge-case rows and all four §4 signature groups are implemented as specified, and every
deviation is recorded in §13.

Two observations that are **not** coverage gaps, recorded so they are not lost:

### Observation 1 — `demultiplex` lost its doc-comment to D6's extraction (documentation only)

**What:** `demux.rs:106-115` is `demultiplex`'s original doc-comment (the 5-step algorithm plus
"Also writes a summary file with per-barcode counts"). D6 inserted `demux_base_name`'s own
doc-comment immediately after it (`:116-121`), so the whole block now documents
`demux_base_name`, and `pub fn demultiplex` at `:137` has no doc-comment at all.

**Why it is not a coverage gap:** the plan never specifies `demux.rs`; D6 is itself a documented
deviation and its stated claim — "pure refactor, behaviour-preserving" — holds. Nothing
executable is affected, and `cargo doc` warnings are not gated.

**Why it is still worth one line:** this is precisely the defect class the plan's own Step 2 went
out of its way to repair ("removing it rejoins the comment to `resolve_clump_layout`"), so
leaving the mirror image of it in `demux.rs` is inconsistent with the change's own standard.
Fix is a one-line move. Flagged for the code reviewers, who own that call.

### Observation 2 — §3.5 row 11 is pinned by a manual run, not by `cargo test`

**What:** the row "same basename, different dirs, **no** `-o` — hardtrim → reject" is the
documented `--hardtrim5 N */*.fastq.gz` pattern and hardtrim's widest collision class. Both
hardtrim FASTQ rejection tests (`hardtrim5_rejects_same_basename_across_dirs`,
`hardtrim3_rejects_same_basename_across_dirs`) add `-o`, so `cargo test` exercises only the
`output_dir`-join arm of `hardtrim_output_name`, never the bare-filename/CWD arm.

**Why it is not a coverage gap:** V2.5/V2.6 as written require `current_dir(tempdir)`, which the
tests do set; `-o` was added for the reason V6.2 states in the plan itself (the "nothing written"
residue check needs a dedicated directory). The uncovered arm cannot silently regress, because
the candidate builder and the writer share one namer through `HardtrimEnd` (A9 + §4.2), so a
namer change moves both together. And I re-ran the no-`-o` case directly this audit — repro (c)
below — where it rejects with the mode-specific hint and writes nothing to the CWD.

## Test verification

`cargo test` from the crate root, single run, no filters: **521 passed, 0 failed, 0 ignored**
(399 lib + 122 integration). New tests: 12 unit in `src/io.rs`, 3 unit in `src/cli.rs`,
20 integration in `tests/integration_output_collision.rs`.

| Test | File | Status |
|---|---|---|
| `preflight_accepts_distinct_outputs` (V1.1) | `src/io.rs` | PASS |
| `preflight_rejects_identical_outputs` (V1.2) | `src/io.rs` | PASS |
| `preflight_rejects_case_only_variants` (V1.3) | `src/io.rs` | PASS |
| `preflight_rejects_dot_slash_alias` (V1.4) | `src/io.rs` | PASS |
| `preflight_rejects_absolute_versus_relative_alias` (V1.5) | `src/io.rs` | PASS |
| `preflight_rejects_output_that_aliases_an_input` (V1.6) | `src/io.rs` | PASS |
| `preflight_accepts_dotdot_alias_known_limitation` (V1.7, pins A8) | `src/io.rs` | PASS |
| `preflight_appends_hint_only_when_given` (V1.8) | `src/io.rs` | PASS |
| `preflight_accepts_empty_and_single` (V1.9) | `src/io.rs` | PASS |
| `preflight_rejects_input_alias_across_spellings` (bonus) | `src/io.rs` | PASS |
| `distinct_primary_outputs_imply_distinct_secondary_outputs` (V5) | `src/io.rs` | PASS |
| `primary_output_key_is_coarser_than_secondary_keys` (V5 complement) | `src/io.rs` | PASS |
| `test_validate_single_end_duplicate_input_rejected` (V4) | `src/cli.rs` | PASS |
| `test_validate_hardtrim_duplicate_input_rejected` (V4) | `src/cli.rs` | PASS |
| `test_validate_single_end_distinct_inputs_accepted` (V4 control) | `src/cli.rs` | PASS |
| `se_trim_rejects_shared_stem` (V2.1) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_rejects_dot_slash_alias` (V2.2) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_rejects_same_basename_across_dirs_with_output_dir` (V2.3) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_ubam_rejects_shared_stem` (V2.4) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim5_rejects_same_basename_across_dirs` (V2.5, +hint grep) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim3_rejects_same_basename_across_dirs` (V2.6) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim5_ubam_rejects_shared_stem` (V2.7) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim3_ubam_rejects_shared_stem` (V2.8) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_rejects_output_that_aliases_an_input` (V2.9) | `tests/integration_output_collision.rs` | PASS |
| `paired_rejects_output_that_aliases_an_input` (V2.10) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_rejects_gz_and_bgz_sharing_a_stem` (V2.11) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_rejects_absolute_versus_relative_alias` (bonus) | `tests/integration_output_collision.rs` | PASS |
| `duplicate_input_gets_a_precise_message` (bonus, A6) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_accepts_same_basename_across_dirs` (V3.1) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_accepts_single_input` (V3.2) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim5_accepts_distinct_stems_and_names_output` (V3.3) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim3_accepts_distinct_stems_and_names_output` (V3.4) | `tests/integration_output_collision.rs` | PASS |
| `hardtrim3_ubam_accepts_distinct_stems_and_names_output` (V3.5) | `tests/integration_output_collision.rs` | PASS |
| `se_trim_ubam_accepts_distinct_stems` (V3.6) | `tests/integration_output_collision.rs` | PASS |
| `paired_accepts_two_distinct_pairs` (V3.7) | `tests/integration_output_collision.rs` | PASS |

### V7 gates and reproductions, re-run this audit

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0, no warnings |
| `cargo test` (crate root) | 521 passed, 0 failed |
| `cargo build --release` then §2.2 (a) shared stem | exit 1, "Output path collision", only the 2 inputs left |
| §2.2 (b) `.gz` + `.bgz` | exit 1, "Output path collision" |
| §2.2 (c) hardtrim5 same basename, 2 dirs, **no `-o`** | exit 1, "one invocation per input", nothing in CWD |
| §2.2 (d) hardtrim5 uBAM arm | exit 1, "Output path collision", nothing written |
| §2.2 (e1) absolute vs relative | exit 1, "Output path collision", only the 2 inputs left |
| §2.2 (e2) `./` prefix | exit 1, "Output path collision", only the 2 inputs left |
| §2.2 (f1) SE output aliases an input | exit 1, "which is also one of its inputs", the named input's 40 PRIOR reads intact |
| §2.2 (f2) `--paired` output aliases an input | exit 1, "which is also one of its inputs", `a_R1_val_1.fq`'s 40 reads intact |
| §10 V3.1 positive control (no `-o`) | exit 0; `dirA` 40 DIRA / 0 DIRB, `dirB` 40 DIRB / 0 DIRA |
| A6 same file twice | exit 1, "was given more than once", 0 occurrences of "APFS/NTFS" |
| A8 residual `a/../x` + `y` | exit 0 — accepted, as A8 states |

Transcript and driver script:
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d0ee0171-fc72-476b-bdba-c46c90a5bacd/scratchpad/reaudit/v7.sh`

## Verdict

**Verdict:** COMPLETE

All 25 audited items are DONE: §5 steps 1–13, §10 V1–V7 (with 9/9, 11/11 and 7/7 named sub-cases
present and individually located, not merely "a test file exists"), all 11 call sites with the
3-re-pointed / 4-converted / 4-new arithmetic reconciled against `72624c7`, the "no fifth hole"
claim re-derived independently from the 15 `cli.input` traversals, all 15 §3.5 edge-case rows,
all four §4 signature groups, and D1–D6 each verified against the code.

The two items the earlier audit raised are closed. `norm_path`'s doc-comment names
`collision_key` and the stale `main::run` locator is gone (`io.rs:33-45`). V5 is now
`distinct_primary_outputs_imply_distinct_secondary_outputs`, which asserts A2 in the direction
the plan states across the full `--basename`/`--dont_gzip`/`-o` matrix and reaches the demux stem
through D6's extracted function rather than a copy of its formula.

Nothing is blocked and nothing needs implementing to satisfy the plan. Two non-blocking
observations are recorded above for the code reviewers: `demultiplex`'s doc-comment is now
attached to `demux_base_name` (documentation only, one-line fix), and §3.5 row 11 — hardtrim
colliding with no `-o` — is pinned by a manual reproduction rather than by `cargo test`.

