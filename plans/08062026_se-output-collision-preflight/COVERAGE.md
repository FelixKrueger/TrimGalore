# Plan Coverage Report

**Mode:** B (code vs. implementation plan — the plan's own §5 is the ledger source; no separate `IMPL.md`, by design)
**Plan(s):** `plans/08062026_se-output-collision-preflight/PLAN.md` (revision v2, §13 notes D1–D6 + post-audit closures)
**Date:** 2026-08-07
**Verdict:** COMPLETE

> **Provenance.** This is a third, independent audit, written without reading `COVERAGE_round1.md`,
> `SNAPSHOT_coverage_reaudit.md`, or any `CODE_REVIEW_*.md`. It replaces the previous `COVERAGE.md`,
> which is preserved byte-identically as `SNAPSHOT_coverage_reaudit.md`
> (both md5 `0e709a371a890fa482660441a9f7272c`) — nothing was lost.
>
> `PLAN.md` was edited while this audit was in progress: it grew from 703 to 726 lines, gaining **D6**
> and the "Post-audit gap closures" section. The ledger below is verified against the **726-line**
> version. This audit's one deviation finding is documented there, which is why the verdict is
> COMPLETE rather than INCOMPLETE.

## Summary

- Total items: 65
- DONE: 64
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 1 — item 47 (V5), documented at §13 "Post-audit gap closures" ¶2, therefore acceptable

Audited state: `fix/383-output-collision-preflight` off `dev` @ `72624c7`, uncommitted.
Actual diff: **7 files changed, +595/−111**, plus the new `tests/integration_output_collision.rs`
(604 lines, 20 tests).

Gates re-run by this audit:

| Gate | Result |
|---|---|
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean (exit 0) |
| `cargo test` (crate root) | **521 passed, 0 failed** (399 lib + 122 integration) |
| `cargo build --release` | succeeds |
| `.github/workflows/ci.yml` parses | valid YAML, **41** steps in `validation` |

## Coverage ledger

### §5 — the 13 implementation steps

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `collision_key` + `preflight_output_collisions` in `io.rs`, after `norm_path`; extend `norm_path`'s doc-comment to name the new caller | Step 1 | DONE | `src/io.rs:47-51` (key), `:53-90` (helper). `norm_path` doc-comment names `collision_key` at `io.rs:37-38`, and the stale "`main::run` paired-end" locator is gone. |
| 2 | Delete `preflight_collision_bam`; re-point its three calls | Step 2 | DONE | Symbol absent from the whole tree. Three calls now at `main.rs:548`, `:585`, `:633`. The incidental fix landed: `resolve_clump_layout`'s doc-comment rejoined its function. |
| 3 | Convert the four hand-rolled copies, preserving each site's semantics | Step 3 | DONE | `main.rs:513` (clump_only SE FASTQ), `:734` (paired FASTQ — candidates accumulated across **all** chunks into one `Vec` before a single call, so cross-pair collisions are still caught; `--retain_unpaired` and `--passthrough` extras kept), `:1825` (paired uBAM), `:2430` (`run_specialty_paired`). Every site passes `guarded_inputs()`. Both `--output-dir` message sites are gone (§2.5). |
| 4 | `HardtrimEnd` enum + `pub` on both namers; update the 4 internal call sites and the bgz test | Step 4 | DONE | `specialty.rs:19-31` (enum, `as_str`, derives `Debug, Clone, Copy, PartialEq, Eq`); `pub fn` at `:442` and `:465`; call sites `:45`, `:88`, `:139`, `:191`; test updated at `:838-846`. |
| 5 | `planned_hardtrim_outputs` matching on `cli.output_format` | Step 5 | DONE | `main.rs:69-88`. `std::path::PathBuf` spelled in full per §4.3/§8. |
| 6 | Guard `--hardtrim5`, first statement in the block, with `HARDTRIM_HINT` | Step 6 | DONE | `main.rs:358-364`, immediately before the write loop at `:365`. `HARDTRIM_HINT` at `:55-59`, naming the real remedies (one invocation per input, or distinct basenames). |
| 7 | Guard `--hardtrim3` identically with `HardtrimEnd::Three` | Step 7 | DONE | `main.rs:389-395`, before the write loop at `:396`. |
| 8 | Guard SE trim FASTQ in the `} else {` arm before the loop | Step 8 | DONE | `main.rs:791-799`. The namer call is character-identical to the writer's, closing A9 by construction. |
| 9 | Guard SE trim uBAM, same shape with `single_end_bam_output_name` | Step 9 | DONE | `main.rs:1863-1869`. |
| 10 | Reject duplicate SE inputs in `cli.rs` with its own message | Step 10 | DONE | `cli.rs:625-642`, with D1's corrected predicate `!self.paired && !self.clock && self.implicon.is_none()`. |
| 11 | Tests: unit in `io.rs` (V1) and `cli.rs` (V4); new `tests/integration_output_collision.rs` (V2, V3) | Step 11 | DONE | 10 new `io.rs` preflight unit tests + 2 A2-table tests; 3 new `cli.rs` tests; 20 integration tests. Conventions match `tests/integration_paired_format_guard.rs`: `env!("CARGO_BIN_EXE_trim_galore")`, `tempdir(tag)` keyed on `process::id()`, `(bool, String)` from `Command::output()`. The `current_dir` hazard is handled more strongly than the plan asked — fixtures are built in-test, so no case references `test_files/` at all. |
| 12 | CI: three validation steps | Step 12 | DONE | `ci.yml:642-696`, inserted after the existing `/tmp/rust_case/out` guard as §10 V6 specified. |
| 13 | CHANGELOG: `#### Bug fixes` for #383 + the two defects; `#### Changes` for A6 and the `--output-dir` correction | Step 13 | DONE | `CHANGELOG.md:63-117`. All four required entries present; `#### Changes` reused, not duplicated. |

### §4 — specified signatures

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 14 | `fn collision_key(p: &Path) -> String` = `norm_path(absolute(p).unwrap_or_else(…))` | §4.1 | DONE | Matches the plan exactly, including the `Err` fallback that degrades to the raw-string key rather than skipping the check. |
| 15 | `pub fn preflight_output_collisions(planned: &[PathBuf], inputs: &[PathBuf], hint: Option<&str>) -> Result<()>` | §4.1 | DONE | Signature exact. Doc-comment enumerates the call sites as §4.1 requires and adds "A new dispatch path needs a call here too." Inputs are held as `HashMap<String, &PathBuf>` rather than the planned `HashSet<String>` — that is what lets the alias message name the offending input, which §3.2 requires. Planned paths use `HashMap::with_capacity` per §6. |
| 16 | `pub enum HardtrimEnd { Five, Three }`, `.as_str()` → `"5prime"`/`"3prime"`; both namers `pub` | §4.2 | DONE | As specified. `fn` → `pub fn` trips no lint, confirmed by the clean clippy run (§8's prediction holds). |
| 17 | `fn planned_hardtrim_outputs(cli, keep, end, output_dir, gzip) -> Vec<std::path::PathBuf>` | §4.3 | DONE | Five arguments, under clippy's `too_many_arguments` threshold. |
| 18 | A6 duplicate-SE-input check in `Cli::validate`, own message in the style of `cli.rs:535` | §4.4 | DONE | "Input file {} was given more than once (arguments {} and {}). List each input once — trimming it twice would write the same output file twice." Placed at `:629`, after the `--basename` check rather than literally beside the duplicate-pair check, but inside `Cli::validate` and before `ensure_output_dir` — both properties the plan required (A7). |

### §10 V1 — unit tests in `src/io.rs` (9 specified, 9 found)

| # | Item | Source | Status | Test |
|---|------|--------|--------|------|
| 19 | distinct outputs → `Ok` | V1.1 | DONE | `preflight_accepts_distinct_outputs` |
| 20 | byte-identical outputs → `Err`, message names both | V1.2 | DONE | `preflight_rejects_identical_outputs` — asserts the duplicate wording **and** that the path is named |
| 21 | case-only variants → `Err` | V1.3 | DONE | `preflight_rejects_case_only_variants`; the plan's "unreachable from integration — cannot coexist on APFS" rationale is carried into the doc-comment |
| 22 | `./x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Err` | V1.4 | DONE | `preflight_rejects_dot_slash_alias` |
| 23 | `<cwd>/x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Err` | V1.5 | DONE | `preflight_rejects_absolute_versus_relative_alias` |
| 24 | planned == an input → `Err` with the **alias** wording, not the duplicate wording | V1.6 | DONE | `preflight_rejects_output_that_aliases_an_input` — asserts the alias wording present *and* the duplicate wording absent, exactly as specified |
| 25 | `a/../x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Ok`, pinning A8 | V1.7 | DONE | `preflight_accepts_dotdot_alias_known_limitation`, with a comment instructing that A8 be retired deliberately if it ever starts failing |
| 26 | `hint: Some(…)` appears in the duplicate message; `None` does not add it | V1.8 | DONE | `preflight_appends_hint_only_when_given` — both directions asserted |
| 27 | empty `planned` → `Ok` | V1.9 | DONE | `preflight_accepts_empty_and_single` (also covers the single-input no-op) |

Beyond plan: `preflight_rejects_input_alias_across_spellings` — the alias check is keyed the same way, so a `./` spelling still catches it.

### §10 V2 — integration rejections (11 specified, 11 found)

All assert non-zero exit, the `Output path collision (case-insensitive, for APFS/NTFS safety)` prefix,
the case-appropriate wording, **and** that nothing was written — via `assert_rejected_cleanly`
(output directory empty) or `assert_dir_holds_only` (exact directory contents, for the no-`-o` cases
that D3 introduced). This is the "stronger than the primary is absent" contract §10 V2 asked for.

| # | Item | Source | Status | Test |
|---|------|--------|--------|------|
| 28 | SE trim FASTQ, `sample.fastq.gz` + `sample.fq.gz` | V2.1 | DONE | `se_trim_rejects_shared_stem` |
| 29 | SE trim FASTQ, `./sample.fastq.gz` + `sample.fq.gz` (defect 1) | V2.2 | DONE | `se_trim_rejects_dot_slash_alias` — runs without `-o` per D3, so the output inherits the input's spelling and the test cannot pass for the wrong reason |
| 30 | SE trim FASTQ, `dirA/same` + `dirB/same` with `-o <third>` | V2.3 | DONE | `se_trim_rejects_same_basename_across_dirs_with_output_dir` — the plan's highest-value single case, pinning `output_dir` (not `None`) in Step 8 |
| 31 | SE trim uBAM, colliding stems | V2.4 | DONE | `se_trim_ubam_rejects_shared_stem` |
| 32 | `--hardtrim5 20` FASTQ, same basename in two dirs, `current_dir` | V2.5 | DONE | `hardtrim5_rejects_same_basename_across_dirs` |
| 33 | `--hardtrim3 20` FASTQ, same shape | V2.6 | DONE | `hardtrim3_rejects_same_basename_across_dirs` |
| 34 | `--hardtrim5 20 --output-format ubam` | V2.7 | DONE | `hardtrim5_ubam_rejects_shared_stem` |
| 35 | `--hardtrim3 20 --output-format ubam` | V2.8 | DONE | `hardtrim3_ubam_rejects_shared_stem` |
| 36 | SE trim, output aliases an input — alias message (defect 2) | V2.9 | DONE | `se_trim_rejects_output_that_aliases_an_input` |
| 37 | `--paired`, output aliases an input (the §2.2(f) four-file invocation) | V2.10 | DONE | `paired_rejects_output_that_aliases_an_input` — also asserts the pre-existing `a_R1_val_1.fq` still holds its 40 `OLD1` reads: byte-level proof the converted paired site gained the protection |
| 38 | `.gz` + `.bgz` sharing a stem (#382 regression pin) | V2.11 | DONE | `se_trim_rejects_gz_and_bgz_sharing_a_stem` |

Beyond plan: `se_trim_rejects_absolute_versus_relative_alias` — the second of D3's two defect-1 integration tests.

### §10 V3 — integration acceptances (7 specified, 7 found)

Every case asserts content or filename, never mere existence.

| # | Item | Source | Status | Test |
|---|------|--------|--------|------|
| 39 | SE trim, `dirA/same` + `dirB/same`, **no** `-o` → exit 0, two outputs, content-attributed | V3.1 | DONE | `se_trim_accepts_same_basename_across_dirs` — asserts all four directions (dirA holds 40 `DIRA_` and 0 `DIRB_`; dirB the converse). Zero crosstalk, and the case that would have passed on existence alone while #383 was live |
| 40 | SE trim, single input → exit 0, one output | V3.2 | DONE | `se_trim_accepts_single_input` — asserts 40 reads present |
| 41 | `--hardtrim5 20 -o <tmp>`, two distinct stems → both `*.20bp_5prime.fq` | V3.3 | DONE | `hardtrim5_accepts_distinct_stems_and_names_output` |
| 42 | `--hardtrim3 20 -o <tmp>`, two distinct stems → both `*.20bp_3prime.fq` | V3.4 | DONE | `hardtrim3_accepts_distinct_stems_and_names_output` — the first output-producing coverage `--hardtrim3` has ever had, as §10 V3.4 noted |
| 43 | `--hardtrim3 20 --output-format ubam`, two distinct stems → both `*.20bp_3prime.bam` | V3.5 | DONE | `hardtrim3_ubam_accepts_distinct_stems_and_names_output` |
| 44 | SE trim uBAM, two distinct stems → both `<stem>_trimmed.bam` | V3.6 | DONE | `se_trim_ubam_accepts_distinct_stems` — catches an over-rejecting candidate list |
| 45 | `--paired`, two distinct pairs → four `_val_` outputs | V3.7 | DONE | `paired_accepts_two_distinct_pairs` |

The filename assertions on the hardtrim acceptances are what actually close A9 for the enum: a
candidate list built with the wrong end discriminator would pass every V2 rejection test and fail these.

### §10 V4–V7

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 46 | Duplicate SE input rejected with the dedicated message, not the APFS/NTFS wording; distinct SE inputs still validate; existing duplicate-pair and R1≠R2 tests still pass | V4 | DONE | `test_validate_single_end_duplicate_input_rejected` (asserts the precise message **and** the absence of "APFS/NTFS"), `test_validate_single_end_distinct_inputs_accepted` (negative control), plus beyond-plan `test_validate_hardtrim_duplicate_input_rejected` and integration-level `duplicate_input_gets_a_precise_message`. The three `--clock`/`--implicon` tests D1 broke are green in the 399-test lib run. |
| 47 | A2 as a table: `report_name`, `json_report_name`, demux base name, and default-config `<stem>_fastqc.zip` all differ where primaries differ, across `--basename`/`--dont_gzip`/`-o` on and off; comment naming the `--fastqc_args --outdir` exception | V5 | **DEVIATED** (documented) | See Deviations below. |
| 48 | CI step 1 — SE trim collision: non-zero exit, collision grep, empty output dir | V6.1 | DONE | `ci.yml:642-656`. `set +e` / `rc=${PIPESTATUS[0]}` / `set -e` / `test $rc -ne 0` / `grep -q "Output path collision"` / `test -z "$(ls -A …)"` — the existing guards' idiom exactly, on `illumina_10K.fastq.gz` copied to `sample.fastq.gz` + `sample.fq.gz` as specified. |
| 49 | CI step 2 — `--hardtrim5 30` collision, same basename in two subdirs, with `-o` (not `cd`) | V6.2 | DONE | `ci.yml:658-675`. All three assertions present, plus a beyond-plan `grep -q "one invocation per input"` that pins the mode-specific hint actually reaching the user. |
| 50 | CI step 3 — SE output-aliases-input | V6.3 | DONE | `ci.yml:677-696`. Non-zero exit, `grep -q "which is also one of its inputs"`, and a **stronger** no-side-effects check than specified: md5 of the aliased input must match before/after, and `ls -A \| wc -l` must still be 2. |
| 51 | Gates; rebuild; re-run the §2.2 reproductions; the V3.1 positive control | V7 | DONE | `fmt` clean, `clippy -D warnings` clean, `cargo test` 521/521, `cargo build --release` succeeds — all four re-run by this audit. **Scope note:** the six §2.2 shell reproductions were *not* re-run by this audit (the command was denied by the permission gate). They are covered equivalently by V2.1/V2.2/V2.11/V2.5/V2.7/V2.9/V2.10 and the V3.1 positive control, all of which I ran and confirmed passing; §13's "all eight reproductions exit non-zero" remains the implementer's claim, corroborated but not independently reproduced here. |

### §3.5 — edge-case table (15 rows, all traced)

| # | Case | Expected | Status | Trace |
|---|------|----------|--------|-------|
| 52 | One input, no alias | no-op | DONE | `preflight_accepts_empty_and_single`, `se_trim_accepts_single_input` |
| 53 | `sample.fastq.gz` + `sample.fq.gz` | reject | DONE | `se_trim_rejects_shared_stem`, CI step 1 |
| 54 | `sample.fastq.gz` + `sample.fastq.bgz` | reject | DONE | `se_trim_rejects_gz_and_bgz_sharing_a_stem` |
| 55 | `./x.fastq.gz` + `x.fq.gz` | reject (new in v2) | DONE | `preflight_rejects_dot_slash_alias` + `se_trim_rejects_dot_slash_alias` |
| 56 | `/abs/x.fastq.gz` + `x.fq.gz`, cwd `/abs` | reject (new in v2) | DONE | `preflight_rejects_absolute_versus_relative_alias` + `se_trim_rejects_absolute_versus_relative_alias` |
| 57 | `a/../x.fastq.gz` + `x.fq.gz` | accept — known residual | DONE | `preflight_accepts_dotdot_alias_known_limitation`, pinning A8 |
| 58 | Output path equals an input path | reject, own message | DONE | `preflight_rejects_output_that_aliases_an_input`, `se_trim_rejects_output_that_aliases_an_input`, `paired_rejects_output_that_aliases_an_input`, CI step 3 |
| 59 | Same file listed twice | reject in `Cli::validate`, own message, before `ensure_output_dir` | DONE | `cli.rs:629-642`; `test_validate_single_end_duplicate_input_rejected`, and `duplicate_input_gets_a_precise_message`, which also asserts `dup_trimmed.fq` does not exist — i.e. pins the "before `ensure_output_dir`" half |
| 60 | Case-only path variants | reject | DONE | `preflight_rejects_case_only_variants` (unit only — two case variants cannot coexist on APFS, which is the plan's stated reason the helper lives in the library) |
| 61 | Same basename, different dirs, no `-o` — SE trim | accept | DONE | `se_trim_accepts_same_basename_across_dirs` |
| 62 | Same basename, different dirs, no `-o` — hardtrim | reject | DONE | `hardtrim5_rejects_same_basename_across_dirs`, `hardtrim3_rejects_same_basename_across_dirs` |
| 63 | Same basename, different dirs, **with** `-o` | reject; V2.3 pins this | DONE | `se_trim_rejects_same_basename_across_dirs_with_output_dir` |
| 64 | `--basename` with >1 SE input | unchanged, already rejected | DONE | `cli.rs:620-624` verified present and untouched. No dedicated SE unit test exists and the plan required none ("unchanged"); the behaviour is also relied on by the V5 test's `--basename` carve-out comment. |
| 65 | Zero inputs / `--hardtrim5 N --hardtrim3 M` | unreachable / unchanged | DONE | Zero inputs: clap `required`, assumption A5 — no test warranted. Dual hardtrim: unchanged, the `hardtrim5` branch still returns early at `main.rs:387`; explicitly deferred as §11 Open 3 and confirmed in §13 "Left undone, deliberately". |

## Deviations (detail)

### Item 47 (V5): the `<stem>_fastqc.zip` row of the A2 table — DEVIATED, documented

**Expected.** §10 V5 specifies a table asserting that where two inputs' primary output paths differ,
`report_name`, `json_report_name`, the demux base name, **and the default-configuration
`<stem>_fastqc.zip`** all differ too — across `--basename` / `--dont_gzip` / `-o`, each on and off —
with a comment naming the `--fastqc_args --outdir` exception so a future reader does not delete a row.

**Found.** `distinct_primary_outputs_imply_distinct_secondary_outputs` (`src/io.rs`) asserts
`report_name`, `json_report_name`, `clumping_report_name` (beyond plan) and
`crate::demux::demux_base_name`. The flag matrix is complete: `basename in [None, Some("fixed")]` ×
`gzip in [true, false]` × `output_dir in [None, Some]`. The `--fastqc_args "-o DIR"` exception comment
is present as required. The demux assertion goes through the same function `demultiplex` now calls —
which is what D6's `demux.rs` extraction exists for, so the test cannot pass against a stale copy of
the formula. A second beyond-plan test, `primary_output_key_is_coarser_than_secondary_keys`,
establishes the direction that makes the argument non-circular (three spellings collapsing to one
primary while keeping three distinct report names), repairing the flaw the plan itself identified in
v1's V6.

The `<stem>_fastqc.zip` assertion is **absent**, with a nine-line rationale in the test's doc-comment:
the bundled crate derives the artifact name from the primary path handed to it, so there is no second
key of ours to compare, and writing the formula out would only test the formula.

**Why this is acceptable.** §13 "Post-audit gap closures" ¶2 records the substitution explicitly,
including "FastQC's zip name is documented as not-separately-assertable rather than faked." That
satisfies the plan's own convention that deviations be registered in §13, so this is
DEVIATED-with-reason, not a gap. No action required.

*(Audit history: this item was the sole finding that would have made the verdict INCOMPLETE. It was
already closed in §13 before this report was finalised — the §13 text post-dates the 703-line PLAN.md
this audit started from.)*

## Test verification

### Suites (`cargo test`, crate root)

| Suite | Tests | Status |
|---|---|---|
| `src/lib.rs` unit | 399 | PASS |
| `tests/integration_output_collision.rs` (new) | 20 | PASS |
| `tests/integration_adapter2.rs` | 15 | PASS |
| `tests/integration_clump_only.rs` | 12 | PASS |
| `tests/integration_clump_only_ubam.rs` | 15 | PASS |
| `tests/integration_gzip_non_gz_extension.rs` | 7 | PASS |
| `tests/integration_no_args_help.rs` | 1 | PASS |
| `tests/integration_non_restartable_input.rs` | 8 | PASS |
| `tests/integration_paired_format_guard.rs` | 11 | PASS |
| `tests/integration_passthrough.rs` | 2 | PASS |
| `tests/integration_ubam.rs` | 8 | PASS |
| `tests/integration_ubam_out.rs` | 23 | PASS |
| **Total** | **521** | **0 failed** |

### Named tests behind the plan's V-items

| Test name | File | Covers | Status |
|-----------|------|--------|--------|
| `preflight_accepts_distinct_outputs` | `src/io.rs` | V1.1 | PASS |
| `preflight_rejects_identical_outputs` | `src/io.rs` | V1.2 | PASS |
| `preflight_rejects_case_only_variants` | `src/io.rs` | V1.3 | PASS |
| `preflight_rejects_dot_slash_alias` | `src/io.rs` | V1.4 | PASS |
| `preflight_rejects_absolute_versus_relative_alias` | `src/io.rs` | V1.5 | PASS |
| `preflight_rejects_output_that_aliases_an_input` | `src/io.rs` | V1.6 | PASS |
| `preflight_accepts_dotdot_alias_known_limitation` | `src/io.rs` | V1.7 (A8) | PASS |
| `preflight_appends_hint_only_when_given` | `src/io.rs` | V1.8 | PASS |
| `preflight_accepts_empty_and_single` | `src/io.rs` | V1.9 | PASS |
| `preflight_rejects_input_alias_across_spellings` | `src/io.rs` | beyond plan | PASS |
| `distinct_primary_outputs_imply_distinct_secondary_outputs` | `src/io.rs` | V5 | PASS |
| `primary_output_key_is_coarser_than_secondary_keys` | `src/io.rs` | beyond plan | PASS |
| `se_trim_rejects_shared_stem` | `tests/integration_output_collision.rs` | V2.1 | PASS |
| `se_trim_rejects_dot_slash_alias` | same | V2.2 | PASS |
| `se_trim_rejects_absolute_versus_relative_alias` | same | beyond plan (D3) | PASS |
| `se_trim_rejects_same_basename_across_dirs_with_output_dir` | same | V2.3 | PASS |
| `se_trim_ubam_rejects_shared_stem` | same | V2.4 | PASS |
| `hardtrim5_rejects_same_basename_across_dirs` | same | V2.5 | PASS |
| `hardtrim3_rejects_same_basename_across_dirs` | same | V2.6 | PASS |
| `hardtrim5_ubam_rejects_shared_stem` | same | V2.7 | PASS |
| `hardtrim3_ubam_rejects_shared_stem` | same | V2.8 | PASS |
| `se_trim_rejects_output_that_aliases_an_input` | same | V2.9 | PASS |
| `paired_rejects_output_that_aliases_an_input` | same | V2.10 | PASS |
| `se_trim_rejects_gz_and_bgz_sharing_a_stem` | same | V2.11 | PASS |
| `se_trim_accepts_same_basename_across_dirs` | same | V3.1 | PASS |
| `se_trim_accepts_single_input` | same | V3.2 | PASS |
| `hardtrim5_accepts_distinct_stems_and_names_output` | same | V3.3 | PASS |
| `hardtrim3_accepts_distinct_stems_and_names_output` | same | V3.4 | PASS |
| `hardtrim3_ubam_accepts_distinct_stems_and_names_output` | same | V3.5 | PASS |
| `se_trim_ubam_accepts_distinct_stems` | same | V3.6 | PASS |
| `paired_accepts_two_distinct_pairs` | same | V3.7 | PASS |
| `duplicate_input_gets_a_precise_message` | same | V4 (integration) | PASS |
| `test_validate_single_end_duplicate_input_rejected` | `src/cli.rs` | V4 | PASS |
| `test_validate_hardtrim_duplicate_input_rejected` | `src/cli.rs` | beyond plan | PASS |
| `test_validate_single_end_distinct_inputs_accepted` | `src/cli.rs` | V4 negative control | PASS |
| `test_validate_clock_r1_equal_r2_within_pair_rejected` | `src/cli.rs` | D1 regression | PASS |
| `test_validate_clock_duplicate_pair_rejected` | `src/cli.rs` | D1 regression | PASS |
| `test_validate_implicon_duplicate_pair_rejected` | `src/cli.rs` | D1 regression | PASS |
| `bgz_stem_reaches_specialty_output_names` | `src/specialty.rs` | Step 4 (updated for the enum) | PASS |

### Call-site count (§1, §5, §4.1's "eleven call sites")

`grep -c 'naming::preflight_output_collisions' src/main.rs` → **11**.

| Line | Dispatch path | Class |
|---|---|---|
| 360 | `--hardtrim5` | new |
| 391 | `--hardtrim3` | new |
| 513 | `--clump_only` SE FASTQ | converted |
| 548 | `--clump_only` PE uBAM, 1-file interleaved | re-pointed |
| 585 | `--clump_only` PE uBAM, multi-pair | re-pointed |
| 633 | `--clump_only` SE uBAM | re-pointed |
| 734 | `--paired` trim FASTQ | converted |
| 799 | SE trim FASTQ | new |
| 1825 | `--paired` trim uBAM | converted |
| 1869 | SE trim uBAM | new |
| 2430 | `run_specialty_paired` — serves `--clock`, `--implicon`, **and** `--clump_only --paired` FASTQ (routed at `main.rs:483`) | converted |

**4 new + 4 converted + 3 re-pointed = 11**, exactly as §1 claims. All 16 `cli.input` /
`cli.input.chunks(2)` loop sites in `main.rs` were enumerated; the 10 that write output are each
preceded by a guard, and the rest are candidate-building loops feeding those guards.
`preflight_collision_bam` appears nowhere in `src/`.

Guard-before-reader ordering (§3.3 contract item 1) verified on all four trim paths: SE FASTQ guard
`:799` precedes `setup_trimming` at `:805`; paired FASTQ `:734` precedes `:765`; SE uBAM `:1869`
precedes `:1875`; paired uBAM `:1825` precedes `:1838`. So no guarded path opens a reader — and
therefore runs adapter auto-detection — before the collision check.

## Observations (not ledger items)

Recorded for the orchestrator; none is a plan-coverage gap.

1. **§13's headline statistics are stale.** Line 634 says "**520 tests pass** (398 lib + 122
   integration)" — actual **521** (399 lib + 122). Line 638 says "6 files changed, +503/−99" —
   actual **7 files, +595/−111** (D6's `src/demux.rs` is the seventh, and is documented as a
   deviation but not folded into the headline count). Line 720 says "40 steps in `validation`" —
   actual **41**. Worth a one-line refresh before commit; no bearing on correctness.

2. **`demux.rs`'s doc-comment was re-parented by D6.** `demux_base_name` was inserted between
   `demultiplex`'s doc-comment and `demultiplex` itself, so the "1. … 5. Write to the matching
   sample's output file (or NoCode)" list and "Also writes a summary file with per-barcode counts."
   now document `demux_base_name`. Cosmetic; already noted in `PROGRESS.md`.

3. **The helper's doc-comment is marginally broader than the code.** It says the helper is "Called by
   every `main.rs` dispatch path that writes more than one file." The `--paired` +
   single-interleaved-BAM-input path with FASTQ output (`main.rs:677-689`) writes two `_val_` files
   and has no call. §2.1 classified that row as "n/a — one output" and both plan reviewers signed off
   on "no fifth hole exists", so it is outside this ledger — and no collision is constructible there,
   because a single input yields two names differing by a fixed `_val_1`/`_val_2` suffix, neither of
   which can equal the input. One clause of the doc-comment would make it exact.

4. **D5 confirmed as described.** Parsing `ci.yml` shows the three new step names truncating at
   `(issue`, as does the pre-existing `(issue #216)` step — consistent, as D5 says.

5. **Out of scope for coverage, flagged for routing.** `PROGRESS.md` records a unanimous code-review
   High — that `--demux` can still overwrite a named input at exit 0, contradicting §3.2's argument
   that the barcode file is excluded "by argument". This audit neither reproduced nor evaluated it;
   it is a defect in the *plan's* reasoning rather than a shortfall of the code against the plan, so
   it does not affect the verdict. It does mean §3.2 and the corresponding CHANGELOG sentence may
   need revising, which is a plan/CHANGELOG action rather than a coverage one.

## Verdict

**COMPLETE.**

Everything the plan specified is implemented and passing: all 13 of §5's steps; all four §4 signature
groups as written; 9/9 V1 cases, 11/11 V2 cases, 7/7 V3 cases; V4; all three V6 CI steps with the
assertions they were specified to make; V7's gates; and all 15 rows of §3.5's edge-case table traced
to a test, an assumption, or an explicit out-of-scope entry. The eleven call sites are verified by
count and by checking every `cli.input` loop for a preceding guard, with the guard proven to precede
the reader on all four trim paths.

The one deviation — V5's `<stem>_fastqc.zip` row, replaced by a `clumping_report_name` assertion and
a documented not-separately-assertable note — is registered in §13's "Post-audit gap closures", which
is the plan's own register for deviations. Nothing is unresolved.

The implementation exceeds the plan in several places worth keeping: the CI hardtrim step asserts the
mode-specific hint text reaches the user; CI step 3 asserts md5 byte-identity of the aliased input
rather than just its presence; six tests exist beyond those specified; and the integration fixtures
are built in-test, which is a stronger fix for the `current_dir` hazard than the absolutised paths
Step 11 called for.

Two follow-ups, neither blocking: refresh §13's stale statistics (Observation 1), and route the
`--demux` review finding (Observation 5) as a plan/CHANGELOG amendment rather than a coverage item.
