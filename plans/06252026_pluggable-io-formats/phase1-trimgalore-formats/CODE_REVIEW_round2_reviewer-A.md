# Code Review (Round 2 — re-review of fixes) — Phase 1 uBAM output — Reviewer A

**Target:** working-tree diff (`git diff HEAD` + untracked `tests/integration_ubam_out.rs`, `test_files/ubam_out_{se,pe}_REFERENCE.bam`)
**Spec:** PLAN v2.1 (`.../phase1-trimgalore-formats/PLAN.md`)
**Round-1 reports:** `CODE_REVIEW_reviewer-A.md` (NEEDS_REWORK, gated on C1) + `CODE_REVIEW_reviewer-B.md` (APPROVE_WITH_FIXES).
**Reviewer:** A (independent re-review; ran in parallel with B, no coordination)
**Date:** 2026-06-26

## Verification gate (re-run locally)

| Check | Result |
|---|---|
| `cargo build --release` | clean (only sandbox `xcrun`/MacOSX.sdk noise — NOT a code warning; verified the "1 warning" is the SDK-path message, not lint) |
| `cargo test --release` | **349 pass** (326 lib + 0 + 2 + 8 + 13 + 0 doc); 0 failed. Matches the implementer's claim. |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --release --all-targets -- -D warnings` | clean (two new `#[allow(clippy::too_many_arguments)]` on report/dispatch helpers — justified; consistent with the codebase's existing config-threading pattern) |

## Fixes verified

| Round-1 finding | Verdict | Note |
|---|---|---|
| **C1 / A-C1** — `--rename` corrupts last preserved aux tag (CRITICAL, silent data corruption) | **PASS** | `src/fastq.rs:148-156` `append_to_id` now splices `:clip5:…` BEFORE the first `\t` when a tag tail is present, else appends to end. The rename annotation correctly lands in the BAM QNAME (`read1:clip5:AATT` — valid, no whitespace), tail intact. **Empirically confirmed both guard tests FAIL on revert**: the lib test `test_append_to_id_with_tag_tail_splices_before_first_tab` and the integration test `ubam_out_rename_with_preserve_tags_keeps_tags_intact` (which literally surfaced `UB = "GCTAGCTA:clip5:AATTA"` corruption when I reverted the fix). |
| **I1 / AGREE-1** — re-processing collapses `@PG` chain (raw `insert` overwrite) | **PASS** | `src/bam.rs:622-642` now uses `programs_mut().add("trim_galore", pg)` with `.with_context(...)`; `build_pg_map` no longer sets PP (noodles computes it). `last_pg_id` helper deleted (grep confirms no definition). **Empirically confirmed** the guard test `build_output_header_disambiguates_existing_trim_galore_pg` FAILS on revert to raw `insert` (yields one `@PG` instead of `[trim_galore, trim_galore-trim_galore]`). The doc comment correctly documents the 3-deep collision → `io::Result` error surfacing. |
| **I2 / B-I2** — two BAM inputs silently accepted under `--paired --output-format ubam` | **PASS** | `src/main.rs:1625-1635` rejects when either input is detected as `UnalignedBam`, with the exact same message + condition as the FASTQ path (`main.rs:911-921`) — byte-for-byte consistent. Guarded by `ubam_out_two_bam_pair_rejected` (asserts non-zero exit + "two BAM files is not supported"/"interleaved file" in stderr). |
| **I2 / A-I2** — `H:` hex tags emitted on read, rejected on write (asymmetric → mid-stream abort) | **PASS (with a minor coverage note)** | `src/bam.rs:957-968` read-side `append_tag_type_and_value` now bails on `Value::Hex`, matching the `Value::Array` treatment — symmetric fail-fast at input. Correct by inspection. *Coverage note:* the lib test `parse_h_hex_tag_rejected` (bam.rs:1273) exercises the **write-side** parser (`parse_name_and_data`), which already rejected `H` before the fix — so it does not specifically guard the **read-side** `Value::Hex` arm that was changed. Constructing a noodles `Value::Hex` in a unit test is awkward, so this is understandable, but the A-I2 read-side change has no direct regression guard. See NIT-1. |
| **M1 / B-M1** — write path doesn't enforce ACGTN-only / reject `=` / IUPAC→N | **PASS** | New `validate_and_normalize_seq_for_write` (bam.rs:760-785): uppercases, ACGTN pass-through, IUPAC (R/Y/M/K/S/W/B/D/H/V)→N with a `OnceLock`-gated single warning, `=` → hard error, any other byte → error. Wired into `write_record` (bam.rs:577) with context. Five guard tests (bam.rs:1545-1584) cover accept/uppercase/IUPAC/`=`/garbage. The `debug_assert_eq!` seq/qual-length invariant still holds (the helper is a 1:1 byte map, so `normalized_seq.len() == record.seq.len()`). |
| **M3 / B-M3** — missing-qual stderr NOTE (§3.3 step 4) not emitted | **PASS** | `run_ubam_output` (main.rs:1387-1397) emits the NOTE when `any_bam` is true; `any_bam` is computed once in `main()` (main.rs:144) and threaded through the call (main.rs:297). The FASTQ-input note ("no aux fields") is also retained (main.rs:1507). |
| **AGREE-2 + §9 gaps** — two missing §9 test rows | **PASS** | `ubam_out_cores_gt_1_warns_and_proceeds` (asserts exit 0 + "ignored" in stderr + BAM produced) and `ubam_out_preserve_tags_mixed_batch_allowed` (asserts uBAM-side record carries CB, FASTQ-side record does not) both present and passing. Plus the two new lock-in tests `ubam_out_rename_with_preserve_tags_keeps_tags_intact` and `ubam_out_two_bam_pair_rejected`. All 13 integration tests pass. |

**Empirical revert-testing note:** I temporarily reverted C1 (`append_to_id`) and I1 (`build_output_header`) in the working tree to confirm their guard tests actually fail, then restored both files byte-identically (verified `diff` against backup is empty; `programs_mut().add` line present; full suite green again). The working tree is back to the as-submitted state.

---

## NEW findings

### IMPORTANT

#### A-R2-I1 — Golden-reference fixtures `ubam_out_{se,pe}_REFERENCE.bam` are untracked; two reference tests will hard-fail in CI

**Where:** `tests/integration_ubam_out.rs:154` (`ubam_out_se_ubam_input_matches_reference`) and `:171` (`ubam_out_pe_one_ubam_interleaved_matches_reference`) call `assert_ubam_eq(&out_bam, Path::new("test_files/ubam_out_{se,pe}_REFERENCE.bam"))`.

**Problem:** `git status` shows both `test_files/ubam_out_se_REFERENCE.bam` and `test_files/ubam_out_pe_REFERENCE.bam` as **untracked** (`??`), and `git check-ignore` returns exit 1 (NOT gitignored — simply never `git add`ed). Unlike the tag-bearing tests (which `eprintln!("SKIP")` + `return` when their fixture is absent), these two reference tests are **not** gated on file existence — `assert_ubam_eq` → `header_minus_pg(expected)` opens the path unconditionally and will panic on a missing file. They pass locally only because the files are present in the working tree.

**Effect:** on a clean CI checkout (which sees only committed files), both tests panic at file-open. This is a merge-blocker for the test suite.

**Scope caveat:** this is **pre-existing** from the original implementation round (the fix patch didn't introduce or touch these fixtures), but it surfaces under the brief's "test coverage matches / would survive CI" check, and the consequence (red CI) is real. **Fix:** `git add test_files/ubam_out_se_REFERENCE.bam test_files/ubam_out_pe_REFERENCE.bam` before merge. (No code change needed — purely a missing `git add`.) The other three BAM fixtures the new tests depend on — `ubam_test.bam`, `ubam_paired_test.bam`, `ubam_test_with_tags.bam` — ARE tracked, so the C1/B-I2 lock-in tests will run in CI.

### NIT

#### A-R2-N1 — A-I2 read-side `Value::Hex` rejection has no direct regression guard

The lib test named for the H-tag fix (`parse_h_hex_tag_rejected`, bam.rs:1273) actually exercises the **write-side** parser, which rejected `H` before the fix too. The arm that the A-I2 fix changed — read-side `append_tag_type_and_value`'s `Value::Hex => bail!` (bam.rs:957) — is unguarded. Reverting it would not fail any test. Constructing a noodles `Value::Hex<'_>` in a unit test is awkward (it wraps a borrowed hex slice), so this is a low-priority gap; the fix is correct by inspection and symmetric with the already-tested `Value::Array` arm. Optional: a small unit test feeding a `Value::Hex` to `append_tag_type_and_value` directly.

#### A-R2-N2 — `_preserve_tags` dead field (carried over from round 1, N1)

`src/bam.rs:76` `_preserve_tags` is still stored-but-never-read (the tail is the source of truth). Round-1 NIT, untouched, still cosmetic. No action required.

#### A-R2-N3 — Report filename `.bam` leak (B-N2) intentionally NOT fixed

Per the re-review brief, the pre-existing `ubam_test.bam_trimming_report.txt` naming (shared with the FASTQ-output BAM-input path) was flagged out-of-scope in round 1 and intentionally left. Confirmed present in the gate output (line 15: `ubam_test_with_tags.bam_trimming_report.txt`); NOT re-raised as new, per instruction.

---

## PLAN §3 / §9 coverage check

All §3.4a CLI rejection rules present in `src/cli.rs:513-545` (clumpify / passthrough / clock / implicon / demux / retain_unpaired), each with a unit test. §3.4b format-detection-time rule (`--preserve-tags` + all-FASTQ) wired at main.rs:148 keyed on `any_bam`, with `ubam_out_preserve_tags_all_fastq_rejected` integration test. §3.6 tag round-trip, §3.3 flag bits / qual / seq validation, §3.5 header propagation + cores warning + missing-qual NOTE all covered. No unmet §3/§9 promise found other than the untracked-fixture issue above (A-R2-I1), which affects whether two §9 rows actually execute in CI.

---

## Verdict

**APPROVE_WITH_FIXES.**

All seven round-1 findings (C1, I1, I2/two-BAM, I2/H-tag, M1, M3, AGREE-2 §9 rows) are correctly and completely fixed. The two highest-severity fixes (C1 silent corruption, I1 provenance collapse) were empirically revert-tested and their guard tests genuinely catch regressions. The fix patch introduces no new defects in the code itself.

The single remaining blocker is **A-R2-I1**: `git add` the two untracked golden-reference BAM fixtures so `ubam_out_{se,pe}_*_matches_reference` don't hard-fail on a clean CI checkout. This is a one-command fix with no code change. With that done, this is a clean APPROVE.
