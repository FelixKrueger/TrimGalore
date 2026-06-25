# Plan Coverage Report

**Mode:** B (code vs. implementation plan)
**Plan(s):** `plans/06252026_ubam-input-support/IMPL.md` (26 tasks, 40 plan items mapped)
**Date:** 2026-06-25
**Verdict:** INCOMPLETE — 6 tasks pending (T21–T26), 2 unexpected test failures in partially-landed T21/T22.

## Summary

- Total tasks: 26
- DONE: 18 (T1–T20)
- DEVIATED (documented): 2 (T6 PeekReader optional, T9–T13 in-memory builder fallback)
- PARTIAL: 1 (T21/T22 — files exist but contain failing tests; reference files not git-tracked)
- MISSING: 5 (T23, T24, T25, T26; T21 reference-snapshot also untracked)
- Test count check: `cargo test --release` reports **283 passed, 0 failed** (281 lib + 2 integration_passthrough). The implementation log's "283 tests" claim is corroborated, but the **6 tests in `tests/integration_ubam.rs` are not in that count** — that file is not yet `git add`ed (still untracked) AND 2 of its 6 tests fail at runtime when invoked.

## Coverage ledger

| # | Task | Source | Status | Notes |
|---|------|--------|--------|-------|
| T1 | Add `noodles = "=0.88.0"` dep + binary-size baseline | §5 step 1 | DONE | `Cargo.toml:48` has the exact-pinned line; binary size unchanged (0-byte delta per session log; dep was already transitive via fastqc-rust). |
| T2 | Register `bam` + `format` modules in `src/lib.rs` | §5 step 2/3 | DONE | `src/lib.rs:3` `pub mod bam;`, line 10 `pub mod format;`. |
| T3 | Test fixtures + `test_files/README.md` recipe | §5 step 7 | DONE | `test_files/ubam_test.bam` (507 B), `test_files/ubam_paired_test.bam` (746 B). `test_files/README.md` has the samtools recipe. |
| T4 | `InputFormat` enum + `detect_input_format` signature | §4, §5 step 2 | DONE | `src/format.rs:25-34` + `:40`. |
| T5 | Two-stage detection (peek + decompress + `BAM\1` magic) | §3.1, §5 step 2 | DONE | `src/format.rs:40-83`; unit tests `detect_fastq_plain_from_at_sign`, `detect_fastq_gz_from_plain_gzip`, `detect_unaligned_bam_via_decompressed_magic`, `detect_bgzipped_fastq_is_fastq_not_bam` (LOAD-BEARING). |
| T6 | `PeekReader` to avoid double-decompression | §5 step 2.3 | DEVIATED | Documented at `src/format.rs:60-63` as optional future micro-optimization. The IMPL.md plan itself marks this refactor as optional. Re-open-from-byte-0 used instead — correct, slightly slower. |
| T7 | Format-detection negative tests | §5 step 2.4 | DONE | `detect_rejects_random_binary`, `detect_empty_file_errors`, `detect_paired_ubam_via_decompressed_magic` all present in `src/format.rs::tests`. |
| T8 | `BamReader::open` + header parse smoke test | §5 step 3.2 | DONE | `src/bam.rs:86-97`. Test `open_committed_se_fixture` + `open_committed_pe_fixture`. |
| T9 | `bam_record_to_fastq` ID conversion + empty-name guard | §3.2.1 | DEVIATED | Logic at `src/bam.rs:526-534`. Happy-path test `first_record_id_starts_with_at_sign`. Empty-name error test deferred to integration per documented fallback at `src/bam.rs:821-827`. |
| T10 | Seq decode + IUPAC coerce + reject `=` | §3.2.2 | DEVIATED | Logic at `src/bam.rs:561-585`. Happy-path test `seq_contains_only_acgtn`. Error-path tests for `=`, IUPAC, empty-seq deferred to integration per fallback note. |
| T11 | Qual offset + missing-qual + length check | §3.2.3 | DEVIATED | Logic at `src/bam.rs:589-605`. Happy-path test `qual_is_sanger_phred33`. Mismatch-error test deferred. |
| T12 | Flag rejection (rev-comp / secondary / supplementary / per-record unmapped) | §3.2.4, B-Crit-4 | DEVIATED | Logic at `src/bam.rs:508-524`. All four guards present. Negative tests deferred to integration per fallback note. |
| T13 | Tag preservation (user-order, missing-skip-silent) | §3.2.5 | DEVIATED | Logic at `src/bam.rs:536-558` + `append_tag_type_and_value:622-651`. No happy-path unit test against a tagged fixture — deferred to T23 (which is itself MISSING). |
| T14 | `BamReader::open_threaded` + bounded channel | §5 step 3.4 | DONE | `src/bam.rs:102-104` + `open_threaded_inner:323-398`. Test `threaded_reads_match_direct`. |
| T15 | `open_paired_interleaved` w/ MAX_SLACK=1024 | §3.3, §5 step 3.5 | DEVIATED | Logic at `src/bam.rs:131-321`. `MAX_SLACK = 1024` at line 40, `GROUPED_INPUT_ERR` at 43-48, bound-check at 279-281. Happy-path test `paired_interleaved_de_interleaves_fixture`. MAX_SLACK trip-test deferred (would need an in-memory grouped BAM builder). |
| T16 | `--preserve-tags` CLI flag + tag-name validation | §5 step 4 | PARTIAL | Field at `src/cli.rs:189-190`, parser at `:386-402`, "ALL" rejection at `:387`. **No unit tests** for the parser (`preserve_tags_parses`, `preserve_tags_rejects_bad_name`, `preserve_tags_rejects_all_keyword`) — IMPL.md prescribed three; none found. Validation happens at clap-derive time, so it works behaviourally — just lacks the unit-test layer. |
| T17 | `Cli::validate()` rules (`--passthrough` + BAM, `--paired` + 1 file) | §5 step 4.2 | DONE | Rules at `src/main.rs:143-161` (note: implemented in main.rs after `sanity_check_any`, not in `Cli::validate()` — this is a documented deviation in Session-3 log). CLI test `test_validate_paired_single_input_accepted_at_validate_layer` at `src/cli.rs:759`. |
| T18 | `RecordSource` trait + polymorphic `parallel.rs` | §5 step 5.3 | DONE | Trait at `src/fastq.rs:213-215`, `FastqReader` impl at `:217-221`, `BamReader` impl at `src/bam.rs:20-24`. `parallel.rs` takes `Box<dyn RecordSource>` (lines 94-96, 646-648, 904, etc.). Test helper `open_fq` at `src/parallel.rs:1227`. |
| T19 | `sanity_check_any` format dispatch | §5 step 5.4 | DONE | `src/main.rs:26-40`. Called at `:128`. BAM path peek-reads first record. |
| T20 | `main.rs` dispatch + reader factory + paired-uBAM helper | §5 step 5.1/5.2/5.3 | DONE | Factories `open_threaded_reader:44-57` + `open_sync_reader:61-71`. SE dispatch at lines 656/667. PE-uBAM helper `run_paired_ubam_single_file` at `:1183-1268`. Two-file paired-BAM rejection at `:870`. |
| T21 | Integration test: `single_end_ubam` golden | §6.2 | PARTIAL | `tests/integration_ubam.rs:46-63` exists. Reference at `test_files/ubam_test_trimmed_REFERENCE.fq`. **Files are git-untracked** (status `??`). Test passes when run. |
| T22 | Integration test: `paired_end_interleaved_ubam` + clumpify + specialty | §6.2 | PARTIAL | Tests at `integration_ubam.rs:66-141`. `paired_end_interleaved_ubam_matches_reference` + `clumpify_plus_ubam_runs_clean` PASS. **`hardtrim5_plus_ubam_runs_clean` FAILS** — specialty hardtrim path does not go through the format-aware reader: errors with `"stream did not contain valid UTF-8"`. **`passthrough_plus_ubam_rejected` FAILS** — gets a different error message ("multi-pair input with passthrough not yet implemented") than the expected `--passthrough is not supported with uBAM input`. |
| T23 | Integration test: `preserve_tags_roundtrip` golden | §6.2 | MISSING | No `preserve_tags_roundtrip_matches_golden` test. No `test_files/ubam_test_with_tags.bam` fixture. No `test_files/ubam_test_with_tags_REFERENCE.fq.gz`. |
| T24 | CI job `validation-ubam` | §6.3 | MISSING | `grep -n validation-ubam .github/workflows/ci.yml` → 0 hits. `scripts/compare_fastq_tuples.sh` not present. |
| T25 | Reproducibility check | §9 | MISSING | No verification step run/recorded. (No code change required; verification only.) |
| T26 | Documentation (CLAUDE.md / README.md / CHANGELOG.md) | §5 step 8 | MISSING | `grep -n uBAM CLAUDE.md README.md CHANGELOG.md` → 0 hits. None mention uBAM input. |

## Plan-item coverage (40 items → 26 tasks)

All 18 done tasks cover plan items 1–4, 7–14, 16–17, 19–25, 27–31 (the items mapped to T1–T20). Items mapped to T21–T26 (34–40) are not addressed.

## Gaps (detail)

### T16 — missing parser unit tests (PARTIAL)
**Expected:** Three unit tests (`preserve_tags_parses`, `preserve_tags_rejects_bad_name`, `preserve_tags_rejects_all_keyword`) per IMPL.md.
**Found:** Field + parser fn exist and work via clap; no unit tests for the parser fn.
**Gap:** Add the three `#[test]` cases in `src/cli.rs::tests`.

### T21 — golden reference untracked (PARTIAL)
**Expected:** Committed golden reference for SE uBAM.
**Found:** `test_files/ubam_test_trimmed_REFERENCE.fq` exists but is untracked (`git status` shows `??`). Same for `tests/integration_ubam.rs`.
**Gap:** `git add tests/integration_ubam.rs test_files/ubam_*_REFERENCE.fq`.

### T22 — two failing integration tests (PARTIAL)
**Expected:** `hardtrim5_plus_ubam_runs_clean` and `passthrough_plus_ubam_rejected` pass.
**Found:**
- `hardtrim5_plus_ubam_runs_clean` FAILS: `Error: stream did not contain valid UTF-8`. Specialty modes (`hardtrim5/3 / clock / implicon`) at `src/main.rs:178-...` loop over `cli.input` calling `specialty::hardtrim5(input, …)` — that path uses a FASTQ-direct reader, not the format-aware factory. Per IMPL.md row #33 plan item, specialty + BAM should "just work" but does not.
- `passthrough_plus_ubam_rejected` FAILS: the BAM-passthrough check at `main.rs:149-154` is reached only AFTER `Cli::validate`'s multi-pair check, which fires first and emits a different error message. The test's assertion on stderr text mismatches.
**Gap:** Either (a) route specialty modes through `open_sync_reader` / `open_threaded_reader` (fix the wiring) or (b) explicitly reject specialty+BAM (update PLAN §3.4). Adjust the passthrough test to either assert on the actually-fired error or move the BAM-passthrough guard before the multi-pair check.

### T23 — preserve-tags golden snapshot (MISSING)
**Expected:** Tagged fixture + golden reference + integration test.
**Found:** Nothing — no fixture, no reference, no test.
**Gap:** Build `test_files/ubam_test_with_tags.bam`, run `trim_galore --preserve-tags CB,UB ...`, commit `test_files/ubam_test_with_tags_REFERENCE.fq[.gz]`, add `preserve_tags_roundtrip_matches_golden` to `tests/integration_ubam.rs`.

### T24 — CI validation-ubam job (MISSING)
**Expected:** Job in `.github/workflows/ci.yml` that runs `samtools fastq` vs trim_galore-direct, tier-1 content-tuple GATING + tier-2 md5 INFORMATIONAL.
**Found:** No job, no helper script.
**Gap:** Add the YAML block per IMPL.md §Task 24, plus `scripts/compare_fastq_tuples.sh`.

### T25 — reproducibility check (MISSING)
**Expected:** Run `SOURCE_DATE_EPOCH=1700000000 cargo build --release` twice, diff the binaries.
**Found:** Not run/documented.
**Gap:** Execute the verification and document the result in the implementation log.

### T26 — documentation (MISSING)
**Expected:** CLAUDE.md module-map row for `src/bam.rs` and `src/format.rs`; README.md uBAM usage example; CHANGELOG.md entry.
**Found:** None of these.
**Gap:** Add per IMPL.md §Task 26.

## Test verification (Mode B)

| Test name | File | Status |
|-----------|------|--------|
| `cargo test --release` (281 lib + 2 integration_passthrough) | crate root | **PASS** — 283 passed, 0 failed (matches the IMPL.md claimed count). |
| `cargo test --release --test integration_ubam` (6 tests) | tests/integration_ubam.rs | **PARTIAL** — 4 pass, 2 FAIL: `hardtrim5_plus_ubam_runs_clean`, `passthrough_plus_ubam_rejected`. File is untracked. |
| All `src/format.rs::tests` (7 tests) | src/format.rs | PASS (via cargo test --release). |
| All `src/bam.rs::tests` (10 tests including paired_interleaved_de_interleaves_fixture) | src/bam.rs | PASS. |
| Clippy `cargo clippy --release --all-targets -- -D warnings` | — | Not re-verified now; session log claims clean. |

## Verdict

**INCOMPLETE — 6 unresolved items.**

The user-known pending tasks T21–T26 are confirmed as the gap surface, with two surprises:

1. **T21/T22 are partially landed** (untracked file `tests/integration_ubam.rs` plus three untracked reference `.fq` files in `test_files/`). Two of the six tests in that file fail: `hardtrim5_plus_ubam_runs_clean` (real wiring bug — specialty modes still use FASTQ-direct readers) and `passthrough_plus_ubam_rejected` (validation order swallows the BAM-specific error).
2. **T16 unit tests for `parse_sam_tag_name` are missing** — the parser itself works (verified via clap-derive behaviour), but the three prescribed test cases were not added.

Specific items the implementer should address next:
- Resolve the wiring bug for specialty modes + BAM (or change the PLAN to reject specialty+BAM with a clear error).
- Reorder the validation rules so BAM-passthrough check fires before the multi-pair check (or update the test assertion).
- `git add` the untracked integration test + reference files.
- Add the three `parse_sam_tag_name` unit tests (T16 prescription).
- Build the tagged fixture + reference (T23), wire the CI job (T24), run the reproducibility check (T25), update CLAUDE.md / README.md / CHANGELOG.md (T26).
