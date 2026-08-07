# Plan Coverage Report

**Mode:** B (code vs. plan — §5 is the implementation plan; no `IMPL.md` by design)
**Plan(s):** `plans/08062026_se-output-collision-preflight/PLAN.md` (v2)
**Date:** 2026-08-07
**Verdict:** INCOMPLETE — 2 items unresolved

## Summary

- Total items: 22
- DONE: 19
- PARTIAL: 2
- MISSING: 0
- DEVIATED: 1 (documented in §13 as D1)

Scope per instruction: §5 steps 1–13, §10 V1–V6 with named sub-cases counted, and the
eleven call sites cross-checked against every loop over `cli.input`. §3.5 and §4 were
spot-confirmed only, and are marked as such.

## Coverage ledger

### Part 1 — §5 implementation outline (steps 1–13)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `collision_key` + `preflight_output_collisions` in `io.rs` | Step 1 | PARTIAL | Both fns present (`io.rs:48-92`), signature matches §4.1 exactly. Sub-clause "Extend `norm_path`'s doc-comment to name the new caller" **not done** — `io.rs:33-43` still lists only `Cli::validate` and the paired pre-flight. |
| 2 | Delete `preflight_collision_bam`, re-point its 3 calls | Step 2 | DONE | Fn gone; calls now at `main.rs:548`, `:585`, `:633`. Incidental doc-comment reunion with `resolve_clump_layout` confirmed. |
| 3 | Convert the four hand-rolled copies | Step 3 | DONE | `main.rs:513` (clump_only SE FASTQ), `:734` (paired FASTQ, candidates accumulated across all chunks into one `Vec` before one call), `:1825` (paired uBAM), `:2430` (`run_specialty_paired`, both names per pair across all pairs). All pass `guarded_inputs()`. Both `--output-dir` message sites gone (§2.5 fix landed). |
| 4 | `HardtrimEnd` enum + `pub` on both hardtrim namers | Step 4 | DONE | See Part 3 note. |
| 5 | `planned_hardtrim_outputs` in `main.rs` | Step 5 | DONE | `main.rs:71-89`, matches on `cli.output_format`, `std::path::PathBuf` spelled in full per §4.3. |
| 6 | Guard `--hardtrim5` | Step 6 | DONE | `main.rs:358-364`, first statement in the block, `HardtrimEnd::Five`, `Some(HARDTRIM_HINT)`. `HARDTRIM_HINT` const at `:56-60`. |
| 7 | Guard `--hardtrim3` | Step 7 | DONE | `main.rs:389-395`, identical with `HardtrimEnd::Three`. |
| 8 | Guard SE trim FASTQ | Step 8 | DONE | `main.rs:791-799`, in the `} else {` arm before the loop, `single_end_output_name`. |
| 9 | Guard SE trim uBAM | Step 9 | DONE | `main.rs:1863-1869`, `single_end_bam_output_name`. |
| 10 | `cli.rs`: reject duplicate SE inputs | Step 10 | DEVIATED (documented) | Predicate widened to `!self.paired && !self.clock && self.implicon.is_none()` — recorded as D1 in §13. |
| 11 | Tests: V1 unit in `io.rs`, V4 unit in `cli.rs`, new `tests/integration_output_collision.rs` (V2, V3) | Step 11 | DONE | 20 new tests in the integration file; detail in Part 2 items 14–17. The "`current_dir` cases must pass absolutised fixture paths" caveat is honoured — the file writes its own fixtures via `write_fastq` into the tempdir instead of referencing `test_files/`, and `tempdir()` is canonicalised (D4). |
| 12 | CI: three validation steps | Step 12 | DONE | Detail in Part 2 item 19. |
| 13 | CHANGELOG entries | Step 13 | DONE | `#### Bug fixes` (`CHANGELOG.md:63-107`) covers #383 plus both defects and the hardtrim CWD widening; `#### Changes` (`:106+`) covers the A6 duplicate-input refusal and the `--output-dir` → `--output_dir` message correction. Both required sections and both required `#### Changes` entries present. |

### Part 2 — §10 validation

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 14 | V1 unit, 9 sub-cases in `io.rs` | §10 V1 | DONE | **9/9 named cases present**, all in `src/io.rs` mod tests: V1.1 `preflight_accepts_distinct_outputs`; V1.2 `preflight_rejects_identical_outputs`; V1.3 `preflight_rejects_case_only_variants`; V1.4 `preflight_rejects_dot_slash_alias`; V1.5 `preflight_rejects_absolute_versus_relative_alias`; V1.6 `preflight_rejects_output_that_aliases_an_input` (asserts alias wording **and** absence of duplicate wording); V1.7 `preflight_accepts_dotdot_alias_known_limitation`; V1.8 `preflight_appends_hint_only_when_given`; V1.9 `preflight_accepts_empty_and_single`. +1 bonus `preflight_rejects_input_alias_across_spellings`. |
| 15 | V2 integration rejections, 11 sub-cases | §10 V2 | DONE | **11/11 present** in `tests/integration_output_collision.rs`: 1 `se_trim_rejects_shared_stem`; 2 `se_trim_rejects_dot_slash_alias`; 3 `se_trim_rejects_same_basename_across_dirs_with_output_dir`; 4 `se_trim_ubam_rejects_shared_stem`; 5 `hardtrim5_rejects_same_basename_across_dirs`; 6 `hardtrim3_rejects_same_basename_across_dirs`; 7 `hardtrim5_ubam_rejects_shared_stem`; 8 `hardtrim3_ubam_rejects_shared_stem`; 9 `se_trim_rejects_output_that_aliases_an_input`; 10 `paired_rejects_output_that_aliases_an_input`; 11 `se_trim_rejects_gz_and_bgz_sharing_a_stem`. +1 bonus `se_trim_rejects_absolute_versus_relative_alias`. Empty-output-dir assertion is real: `assert_rejected_cleanly` (`:95-116`) checks exit≠0, the `Output path collision (case-insensitive, for APFS/NTFS safety)` prefix, the case-specific wording, and an empty dir; the two no-`-o` defect-1 cases use `assert_dir_holds_only` per D3. |
| 16 | V3 integration acceptances, 7 sub-cases | §10 V3 | DONE | **7/7 present**, and every one asserts content or filename, never bare existence: 1 `se_trim_accepts_same_basename_across_dirs` (40 own reads / 0 crosstalk each way); 2 `se_trim_accepts_single_input` (read count); 3 `hardtrim5_accepts_distinct_stems_and_names_output` (`*.20bp_5prime.fq`); 4 `hardtrim3_accepts_distinct_stems_and_names_output` (`*.20bp_3prime.fq`); 5 `hardtrim3_ubam_accepts_distinct_stems_and_names_output` (`*.20bp_3prime.bam`); 6 `se_trim_ubam_accepts_distinct_stems` (`<stem>_trimmed.bam`); 7 `paired_accepts_two_distinct_pairs` (four `_val_` filenames). |
| 17 | V4 unit, duplicate SE input in `cli.rs` | §10 V4 | DONE | `test_validate_single_end_duplicate_input_rejected` (dedicated message **and** `!contains("APFS/NTFS")`), `test_validate_single_end_distinct_inputs_accepted` (negative control), plus `test_validate_hardtrim_duplicate_input_rejected`. Existing duplicate-pair / R1≠R2 tests all still pass. Integration analogue: `duplicate_input_gets_a_precise_message`. |
| 18 | V5 unit, A2 as a table | §10 V5 | PARTIAL | See Gap 2. `primary_output_key_is_coarser_than_secondary_keys` exists and can fail, but covers only `report_name`, `json_report_name`, `clumping_report_name`, `single_end_bam_output_name` — the plan also named the **demux base name** and the **default-config `<stem>_fastqc.zip`**, and required the matrix "across `--basename` / `--dont_gzip` / `-o` on and off". Premise is also inverted relative to §10 V5. The `--fastqc_args -o` exception comment **is** present. |
| 19 | V6 CI, three validation steps | §10 V6 | DONE | `.github/workflows/ci.yml`, all three after the existing `rust_case` guard: `Validate single-end collision pre-flight`, `Validate --hardtrim5 collision pre-flight and its own remediation` (also greps the hint string `one invocation per input`), `Validate an output may not overwrite an input`. Idiom matches (`set +e` / `rc=${PIPESTATUS[0]}` / `set -e` / `test $rc -ne 0` / `grep -q "Output path collision"` / `ls -A` residue). Step 3 substitutes a stronger check — md5 of the aliased input unchanged + exactly 2 files — for the bare `ls -A`. |

### Part 3 — the eleven call sites, and spot-confirmations

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 20 | Eleven `preflight_output_collisions` call sites; no loop over `cli.input` missed | §1, §4.1, §5 | DONE | See the enumeration below. `grep -c naming::preflight_output_collisions src/main.rs` = **11**, matching §1's 3 re-pointed + 4 converted + 4 new. Every writing loop over `cli.input` has a guard immediately in front of it; the two unguarded `cli.input.len() == 1` interleaved-BAM branches produce a single output (§2.1 marks them n/a). No fifth hole. |
| 21 | §3.5 edge-case table | §3.5 | DONE (spot-confirmed) | Spot-confirmed, not deep-audited, per instruction. 12 of the 15 rows are directly pinned by a V1/V2/V3/V4 case above; the three unpinned rows are the pre-existing/unreachable ones (`--basename` with >1 SE input, zero inputs, `--hardtrim5 N --hardtrim3 M`), all marked "unchanged" in the plan. |
| 22 | §4 signatures | §4 | DONE (spot-confirmed) | `preflight_output_collisions(&[PathBuf], &[PathBuf], Option<&str>) -> Result<()>` and `collision_key` match §4.1 verbatim; `HardtrimEnd` + `pub` namers match §4.2; `planned_hardtrim_outputs` matches §4.3 including the fully-spelled `std::path::PathBuf`; the `cli.rs` check sits beside the duplicate-pair check per §4.4. `cargo fmt --check` clean; 520 tests pass. |

**Call-site enumeration** (`src/main.rs`, current working tree):

| Guard | Dispatch path | Writing loop it protects |
|---|---|---|
| `:360` | `--hardtrim5` (both formats) | `:365` | 
| `:391` | `--hardtrim3` (both formats) | `:396` |
| `:513` | `--clump_only` SE FASTQ | `:514` |
| `:548` | `--clump_only` PE uBAM, 1-file interleaved | single output |
| `:585` | `--clump_only` PE uBAM, multi-pair | `:591` |
| `:633` | `--clump_only` SE uBAM | `:634` |
| `:734` | `--paired` trim FASTQ (candidates accumulated across all chunks) | `:744` |
| `:799` | SE trim FASTQ (**new**) | `:800` |
| `:1825` | `--paired` trim uBAM | `:1828` |
| `:1869` | SE trim uBAM (**new**) | `:1870` |
| `:2430` | `run_specialty_paired` — `--clock`, `--implicon`, `--clump_only --paired` FASTQ | `:2433` |

Unguarded-by-design: `:677`/`:1796` (`--paired` single interleaved BAM → one output) and `:475-481` (`--clump_only --paired` with one FASTQ input → hard `bail!`).

## Gaps (detail)

### Gap 1 — item 1: `norm_path`'s doc-comment was not extended (Step 1)

**Expected:** §5 Step 1 — "add `collision_key` and `preflight_output_collisions` …
**Extend `norm_path`'s doc-comment to name the new caller.**"
**Found:** both functions added exactly as specified (`src/io.rs:48-92`). `norm_path`'s
doc-comment (`src/io.rs:33-43`) is byte-unchanged and still enumerates only two callers:
`Cli::validate()`'s `--passthrough` check and "`main::run` paired-end output-collision
pre-flight … (issue #216)". It does not mention `collision_key`, `preflight_output_collisions`,
or #383 — and the paired pre-flight it names no longer lives in `main::run`, so the
existing text is now stale as well as incomplete.
**Gap:** one doc-comment edit. Documentation-only — no behavioural effect, nothing a test
could catch. Not recorded in §13.

### Gap 2 — item 18: V5's A2 table is narrower than specified, and tests the converse

**Expected:** §10 V5 — "For path pairs whose primaries **differ**, assert that
`report_name`, `json_report_name`, **the demux base name**, and the *default-configuration*
**`<stem>_fastqc.zip`** all differ too, **across `--basename` / `--dont_gzip` / `-o` on and
off**."
**Found:** `primary_output_key_is_coarser_than_secondary_keys` in `src/io.rs` mod tests. It
takes three variants whose primaries are **equal** (`d/sample.{fastq.gz,fq.gz,fastq.bgz}`),
asserts the primaries collapse to one, then asserts `report_name`, `json_report_name` and
`clumping_report_name` stay pairwise distinct, and that the uBAM primary also collapses. The
`--fastqc_args -o` exception is named in a comment, as required.
**Gap:** three things. (a) Two of the four named namers are absent — the demux base name
(`demux.rs:142-147`, which A2 explicitly claims to have verified) and the default-config
`<stem>_fastqc.zip`; `clumping_report_name` was substituted for them, which the plan did not
ask for. (b) The `--basename` / `--dont_gzip` / `-o` on-and-off matrix is absent — every
call passes `None, None, true` or `None`. (c) The premise runs the other way: the test
establishes "primaries equal → secondaries distinct", whereas A2's claim, and V5's stated
form, is "primaries distinct → secondaries distinct". The implemented assertion is failable
(so it is not v1's V6 defect returning), but it does not exercise the implication A2 rests
on. Not recorded in §13.

Neither gap affects the fix's behaviour: the 11 guards, both defect fixes, the 27 named
V1–V4 sub-cases and all three CI steps are in place and green.

## Test verification

`cargo test` from the crate root: **520 passed, 0 failed, 0 ignored** across 14 binaries —
exactly the count §13 claims. `cargo fmt --all -- --check` clean.

| Plan case | Test name | File | Status |
|---|---|---|---|
| V1.1 | `preflight_accepts_distinct_outputs` | `src/io.rs` | PASS |
| V1.2 | `preflight_rejects_identical_outputs` | `src/io.rs` | PASS |
| V1.3 | `preflight_rejects_case_only_variants` | `src/io.rs` | PASS |
| V1.4 | `preflight_rejects_dot_slash_alias` | `src/io.rs` | PASS |
| V1.5 | `preflight_rejects_absolute_versus_relative_alias` | `src/io.rs` | PASS |
| V1.6 | `preflight_rejects_output_that_aliases_an_input` | `src/io.rs` | PASS |
| V1.7 | `preflight_accepts_dotdot_alias_known_limitation` | `src/io.rs` | PASS |
| V1.8 | `preflight_appends_hint_only_when_given` | `src/io.rs` | PASS |
| V1.9 | `preflight_accepts_empty_and_single` | `src/io.rs` | PASS |
| V1 bonus | `preflight_rejects_input_alias_across_spellings` | `src/io.rs` | PASS |
| V2.1 | `se_trim_rejects_shared_stem` | `tests/integration_output_collision.rs` | PASS |
| V2.2 | `se_trim_rejects_dot_slash_alias` | same | PASS |
| V2.3 | `se_trim_rejects_same_basename_across_dirs_with_output_dir` | same | PASS |
| V2.4 | `se_trim_ubam_rejects_shared_stem` | same | PASS |
| V2.5 | `hardtrim5_rejects_same_basename_across_dirs` | same | PASS |
| V2.6 | `hardtrim3_rejects_same_basename_across_dirs` | same | PASS |
| V2.7 | `hardtrim5_ubam_rejects_shared_stem` | same | PASS |
| V2.8 | `hardtrim3_ubam_rejects_shared_stem` | same | PASS |
| V2.9 | `se_trim_rejects_output_that_aliases_an_input` | same | PASS |
| V2.10 | `paired_rejects_output_that_aliases_an_input` | same | PASS |
| V2.11 | `se_trim_rejects_gz_and_bgz_sharing_a_stem` | same | PASS |
| V2 bonus | `se_trim_rejects_absolute_versus_relative_alias` | same | PASS |
| V3.1 | `se_trim_accepts_same_basename_across_dirs` | same | PASS |
| V3.2 | `se_trim_accepts_single_input` | same | PASS |
| V3.3 | `hardtrim5_accepts_distinct_stems_and_names_output` | same | PASS |
| V3.4 | `hardtrim3_accepts_distinct_stems_and_names_output` | same | PASS |
| V3.5 | `hardtrim3_ubam_accepts_distinct_stems_and_names_output` | same | PASS |
| V3.6 | `se_trim_ubam_accepts_distinct_stems` | same | PASS |
| V3.7 | `paired_accepts_two_distinct_pairs` | same | PASS |
| V4 | `test_validate_single_end_duplicate_input_rejected` | `src/cli.rs` | PASS |
| V4 | `test_validate_single_end_distinct_inputs_accepted` | `src/cli.rs` | PASS |
| V4 bonus | `test_validate_hardtrim_duplicate_input_rejected` | `src/cli.rs` | PASS |
| V4 (integration) | `duplicate_input_gets_a_precise_message` | `tests/integration_output_collision.rs` | PASS |
| V5 | `primary_output_key_is_coarser_than_secondary_keys` | `src/io.rs` | PASS, but PARTIAL coverage — see Gap 2 |
| V6.1 | `Validate single-end collision pre-flight` | `.github/workflows/ci.yml` | PRESENT (not run here) |
| V6.2 | `Validate --hardtrim5 collision pre-flight and its own remediation` | same | PRESENT (not run here) |
| V6.3 | `Validate an output may not overwrite an input` | same | PRESENT (not run here) |

## Verdict

**Verdict:** INCOMPLETE — 2 items unresolved

Both unresolved items are shortfalls against the plan text that are undocumented in §13; neither
changes the behaviour of the fix.

1. **Item 1 / Gap 1** — Step 1's sub-clause "Extend `norm_path`'s doc-comment to name the new
   caller" was not done. `src/io.rs:33-43` is unchanged and now both incomplete (no mention of
   `collision_key` / `preflight_output_collisions` / #383) and stale (it locates the paired
   pre-flight in `main::run`, which no longer holds it). One comment edit.
2. **Item 18 / Gap 2** — §10 V5's A2 table is narrower than specified and tests the converse
   implication. Missing: the demux base name, the default-config `<stem>_fastqc.zip`, and the
   `--basename` / `--dont_gzip` / `-o` on-and-off matrix. `clumping_report_name` was substituted
   for the two absent namers. The existing assertion is failable, so v1's V6 defect has not
   returned, but A2's actual claim ("primaries distinct → secondaries distinct") is not exercised.

Everything else is in place: all 11 `preflight_output_collisions` call sites (3 re-pointed,
4 converted, 4 new), no unguarded writing loop over `cli.input`, both defect fixes, the
`HardtrimEnd` enum, `planned_hardtrim_outputs`, the `Cli::validate` duplicate-input check (with
D1's documented predicate widening), all 27 named V1–V4 sub-cases, all three V6 CI steps, and the
CHANGELOG entries in both required sections. 520 tests pass; `cargo fmt --check` clean.

Also noted, and correctly scoped out by the plan itself: §11 Open 3 (`--hardtrim5` + `--hardtrim3`
silently ignoring the 3′ trim) and Open 4 (the `modes/hardtrim.md` CWD-output sentence) remain
unaddressed, as §13 "Left undone, deliberately" states.

---

## Addendum — resolution of both gaps

**Everything above this line is the auditing agent's own work, and its verdict was correct when
rendered. This section was written by the orchestrating session and is not an independent
finding.** Both unresolved items were closed after the audit ran; the verdict line above is
therefore a snapshot, not the current state.

**Gap 1 closed.** `src/io.rs:33-43` now names `collision_key` as the absolutising key that every
`main.rs` pre-flight hashes, and drops the stale "`main::run` paired-end" locator the audit
flagged. Written as plain backticks rather than a rustdoc intra-doc link: `norm_path` is `pub`
and `collision_key` is private, so the link form emits `public documentation links to private
item` (confirmed with `cargo doc`; no CI gate runs it, but the repo already carries one such
warning and there was no reason to add a second).

**Gap 2 closed, addressing all three sub-points (a), (b) and (c).**

- **(c) the inverted premise.** New `distinct_primary_outputs_imply_distinct_secondary_outputs`
  asserts A2 in its stated direction — primaries distinct ⇒ every secondary distinct. The
  original coarser-than test is kept as its complement, since coarseness is *why* checking
  primaries covers reports rather than merely coinciding with them.
- **(b) the matrix.** That test runs the full 2 × 2 × 2 over `--basename` / `--dont_gzip` / `-o`,
  skipping the `--basename` combinations where primaries legitimately collapse (the pre-flight
  rejects those, and `cli.rs:620` rejects them earlier still for multi-input SE).
- **(a) the missing namers.** The demux base name is now genuinely covered. Its derivation was
  inline in `demultiplex`, so asserting it would have meant copying the formula into the test —
  a test that cannot detect drift in the code it claims to check. Extracted as
  `demux::demux_base_name`, called from both; recorded as deviation **D6** in `PLAN.md` §13, and
  verified behaviour-preserving by re-running `--demux` end-to-end (output stems and the
  full-filename summary line unchanged). The default-config `<stem>_fastqc.zip` is deliberately
  *not* asserted, and the test says why: the bundled crate derives it from the primary path we
  hand it, so there is no second key of ours to compare, and writing the formula out would only
  test the formula. Its one real exception, `--fastqc_args "-o DIR"`, remains the known residual
  in §9 A2 and §11.

Writing the matrix surfaced a subtlety the plan had not stated, and the new test caught it on
first run. Under `--basename` with no `-o`, `d/fixed_trimmed.fq.gz` and `e/fixed_trimmed.fq.gz`
are distinct primaries whose demux *stems* are identical, because `demux_base_name` drops the
directory. The demux *paths* still differ, because `demultiplex` resolves its directory to `-o`
else the primary's parent. So the stem alone is not the key, and the assertion is now conditioned
on the two primaries sharing a parent — exactly when a stem collision would be a path collision.

**Post-fix gates:** 521 tests pass (was 520), `cargo fmt --all -- --check` clean,
`cargo clippy --all-targets --release -- -D warnings` clean, `cargo doc` adds no new warning.
Negative control on the new test: breaking `report_name` makes it fail, reverting makes it pass.

**Outstanding action:** re-run `plan-manager` in a healthy session. Six consecutive subagent
failures meant this audit's own gaps section, test table and verdict came within one attempt of
being written by the session that wrote the code — which would not have been an independent
check. The ledger above is complete enough to make a re-audit short.
