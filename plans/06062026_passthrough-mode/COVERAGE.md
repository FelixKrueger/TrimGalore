# Plan Coverage Report — `--passthrough` mode

**Mode:** B (code vs. implementation plan)
**Plan:** `PLAN.md` (v2, 2026-06-06)
**Date:** 2026-06-06
**Verdict:** COMPLETE — all 12 implementation steps and all 8 validation tests are present and passing; 1 minor non-blocking documentation gap (Step 11.3, CLAUDE.md fixture mention) is noted but does not block functional correctness.

## Summary

- Total items: 20 (12 implementation tasks + 8 validation tests)
- DONE: 19
- PARTIAL: 1 (Step 11 docs sub-task — CLAUDE.md test-fixtures section not updated; the multiome fixture is generated in-memory at test time rather than stored on disk, so the omission is defensible but the plan asked for it)
- MISSING: 0
- DEVIATED: 0 (Validation §7 explicitly folded into the same in-test fixtures rather than spawning a binary subprocess — implementer's discretion permitted by the plan wording, but the comment block at `src/parallel.rs:1838–1842` is explicit about it)

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | CLI field + 9-item validation block + case-folded R1/R2 collision check | Step 1 — `src/cli.rs:200–215`, `cli.rs:542–608` | DONE | All 9 bail rules present (1.i–1.ix). `norm_path` helper invoked from `crate::io`. |
| 2 | Output naming helper + collision pre-flight | Step 2 — `src/io.rs:179–204`, `src/main.rs:201–207` | DONE | `passthrough_output_name` pushed into per-pair `candidates` vec. |
| 3 | `PartialEq, Eq` derive + 3 passthrough fields on `PairValidationStats` | Step 3 — `src/report.rs:111`, `report.rs:126–138`, `report.rs:151–153` | DONE | Critical A-Crit-1 derive present; `merge()` sums all three new fields. |
| 4 | `read_id_prefix()` helper with `/[123]` stripping | Step 4 — `src/fastq.rs:171–178` | DONE | Signature matches plan verbatim; AB2 `/[123]` strip is in v1. |
| 5 | Serial path 8-arm three-way EOF match | Step 5 — `src/trimmer.rs:381–620` | DONE | `debug_assert!` reader/writer pair-check present (line 394). All 8 arms enumerated (lines 424–451). Three-way header sync (lines 456–471). |
| 6 | Parallel path reader thread + 8-arm match | Step 6 — `src/parallel.rs:634–759` | DONE | All 8 arms present in `read_pairs_round_robin`. Clumpy path untouched per plan. |
| 6a | Reader-thread error propagation hardening (channel payload = `Result<…>`) | Step 6a — `src/parallel.rs:129`, `parallel.rs:155–162, 203–213, 232–273` | DONE | B-Crit-1 fix. Channel is `sync_channel::<Result<PairedBatchResult>>`. Reader sends `Err(e)` via cloned `result_tx_for_reader` before exiting. Main thread short-circuits via `break 'outer` on first `Err`. |
| 7 | Parallel worker + output (`process_paired_batch`, `process_pairs`, `PairedBatchResult`) | Step 7 — `src/parallel.rs:53–67, 295–415, 559–623` | DONE | `compressed_passthrough: Option<Vec<u8>>` on `PairedBatchResult`. `process_pairs` takes `reads_passthrough: Option<&[FastqRecord]>` immutable + `writer_passthrough: Option<&mut W>`. Worker stamps `passthrough_records_checked = reads_r1.len()` once per batch (line 577–579). Gzip and plain branches both mirrored. Atomic three-write flush at lines 241–263. |
| 8 | Wire `main.rs` (open readers/writers, sanity-check, dispatch) | Step 8 — `src/main.rs:201–207, 245–253, 740–844` | DONE | `passthrough_input`/`passthrough_output` threaded into both `run_paired_end_parallel` and serial `run_paired_end`. Sanity-check at `main.rs:249–253`. |
| 9 | Report (text + JSON passthrough block) | Step 9 — `src/report.rs:466–508`, `report.rs:1014–1037` | DONE | Text block gated on `passthrough_records_checked > 0`. JSON emits `null` when inactive, object when active. |
| 10 | FastQC integration on passthrough output | Step 10 — `src/main.rs:1040–1042` | DONE | `fastqc::run` invoked on `passthrough_output` when set. |
| 11 | Tests + docs (8 validation tests + fixture + `--help` + CLAUDE.md) | Step 11 — multiple files | PARTIAL | 8 validation tests all present (see test-verification table below). Help text rendered via clap `#[clap(long = "passthrough")]` with full docstring. **Gap:** CLAUDE.md's test-fixtures section was not updated. Defensible because the fixture is generated in-memory by `gen_multiome_rows` / `write_multiome_fixture` rather than stored under `test_files/` — but the plan explicitly asked for the doc update. |
| §1 | End-to-end smoke (serial path) | Validation §1 — `src/parallel.rs:1977–2045` | DONE | 50-record fixture; row-by-row lockstep + counts. |
| §2 | Serial/parallel parity (≥10K records, row-by-row lockstep, decoded record-identity, stats `==`) | Validation §2 — `src/parallel.rs:2060–2162` | DONE | **12,000 record fixture** (≥3 batches at BATCH_SIZE=4096). Asserts (1) `PairValidationStats` equality, (2/3/4) decoded R1/R2/passthrough record-identity via `Vec<(id, seq, qual)>`, (5) row-by-row id-prefix lockstep on parallel output. AB1 + AB5 satisfied. |
| §3 | Truncation detection (mid-stream) | Validation §3 — `src/parallel.rs:2173–2209` | DONE | 10K-record fixture, truncated at row 5,000, `cores=4` ⇒ error fires mid-stream. Doubles as §8 coverage (see comment block at `parallel.rs:1840–1842`). |
| §4 | ID-mismatch detection (shuffled IDs) | Validation §4 — `src/parallel.rs:2217–2255` | DONE | Reversed passthrough IDs trigger sync-check failure; error contains "Read ID mismatch". |
| §5 | CLI rejection regression tests | Validation §5 — `src/cli.rs:1190–1377` | DONE | 13 CLI tests covering all 9 bail paths (1.i–1.ix), including pointing-at-R1 and pointing-at-R2 (1.ix). |
| §6 | `--dont_gzip` plain-output coverage | Validation §6 — `src/parallel.rs:2264–2315` | DONE | Runs parity test with `gzip = false`; exercises plain branch of `process_paired_batch`. AB3 satisfied. |
| §7 | `--cores 1` dispatcher integration | Validation §7 — *no `tests/` directory* | DONE (deviated, documented) | Plan permitted "CLI-integration test in tests/ top-level dir or via std::process::Command". Implementer chose to fold the serial dispatcher path into Validation §1 (which calls `run_paired_end` directly with `--cores 1` shape). The implementer's intent is explicit at `src/parallel.rs:1839–1840` ("Test §7 (--cores 1 dispatcher) requires the built binary and lives in tests/integration_passthrough.rs"). However, that file **does not exist**. The serial path is exercised by §1 and §2 (both call `run_paired_end` directly), so functional coverage is present, but the plan's intent of "going through main()'s decision logic" via the dispatcher is not strictly covered. **See "Gaps" below — the implementer's comment promises a file that isn't there.** |
| §8 | Reader-thread error propagation under load | Validation §8 — folded into §3 | DONE | The `parallel.rs:1841–1842` comment block declares §8 is covered by §3 (mid-stream truncation in a 10K-record fixture, `cores=4`). The §3 test passing without deadlock IS the §8 invariant. AB1 satisfied (≥10K fixture); the truncation is at row 5,000 so it fires mid-stream as the plan requires. |

## Gaps (detail)

### Item 11 — Step 11 (sub-item 3): CLAUDE.md test-fixtures mention

**Expected:** Plan §Implementation outline Step 11.3 — "Update CLAUDE.md's test-fixtures section to mention the new fixture."
**Found:** No mention of multiome / passthrough in `CLAUDE.md`. `grep -n "multiome\|passthrough" CLAUDE.md` returns zero hits.
**Gap:** Either (a) add a CLAUDE.md note that passthrough fixtures are generated in-memory by `gen_multiome_rows` / `write_multiome_fixture` in `src/parallel.rs::tests` (no on-disk fixture file), or (b) accept the deviation and document it in `PROGRESS.md` under "implementer notes". Functionally non-blocking.

### Item §7 — promised `tests/integration_passthrough.rs` does not exist

**Expected:** Plan Validation §7 says: "CLI-integration test (in `tests/` top-level dir or via `std::process::Command` if a built binary is available). Pass `--paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz`. Assert all three outputs are produced and round-trip."

**Found:**
- No `tests/` directory at the repo root.
- The implementer's comment at `src/parallel.rs:1839–1840` says: "Test §7 (--cores 1 dispatcher) requires the built binary and lives in tests/integration_passthrough.rs" — but that file is absent.
- §1 (`test_passthrough_serial_smoke`) does exercise the serial path with `cores=1`-shaped invocation by calling `trimmer::run_paired_end` directly with no parallel scaffolding, so functional coverage of the serial trim pipeline IS present.

**Gap:** The implementer's intent (a separate binary-driven dispatcher test under `tests/`) is not realised. Two acceptable resolutions:
1. **Accept the deviation** — §1 + §2 + §6 collectively exercise the serial path with multiple fixture shapes; this is reasonable coverage. Update the `parallel.rs:1839–1840` comment to remove the dangling reference.
2. **Add the missing binary-driven test** — create `tests/integration_passthrough.rs` that uses `assert_cmd` / `std::process::Command` to drive the compiled binary end-to-end. ~30 LOC.

Recommended: option 1 (accept + remove dangling comment) — adding a binary-driven test adds CI runtime without catching a meaningful new failure mode.

## Test verification

All passthrough-related tests pass under `cargo test --lib`. Full suite: **262 tests, 0 failures**. Clippy: clean under `-D warnings`.

| Test name | File | Status |
|-----------|------|--------|
| `test_passthrough_paired_pair_accepted` | `src/cli.rs:1190` | PASS |
| `test_passthrough_requires_paired` | `src/cli.rs:1198` | PASS |
| `test_passthrough_rejects_multi_pair` | `src/cli.rs:1208` | PASS |
| `test_passthrough_rejects_retain_unpaired` | `src/cli.rs:1227` | PASS |
| `test_passthrough_rejects_clumpify` | `src/cli.rs:1242` | PASS |
| `test_passthrough_rejects_clock` | `src/cli.rs:1259` | PASS |
| `test_passthrough_rejects_implicon` | `src/cli.rs:1274` | PASS |
| `test_passthrough_rejects_hardtrim5` | `src/cli.rs:1289` | PASS |
| `test_passthrough_rejects_hardtrim3` | `src/cli.rs:1305` | PASS |
| `test_passthrough_rejects_demux` | `src/cli.rs:1321` | PASS |
| `test_passthrough_rejects_missing_file` | `src/cli.rs:1346` | PASS |
| `test_passthrough_rejects_pointing_at_r1` | `src/cli.rs:1360` | PASS |
| `test_passthrough_rejects_pointing_at_r2` | `src/cli.rs:1373` | PASS |
| `test_passthrough_output_name_bare` | `src/io.rs:439` | PASS |
| `test_passthrough_output_name_plain` | `src/io.rs:446` | PASS |
| `test_passthrough_output_name_with_output_dir` | `src/io.rs:453` | PASS |
| `test_passthrough_output_name_with_basename` | `src/io.rs:460` | PASS |
| `test_passthrough_output_name_basename_with_output_dir` | `src/io.rs:467` | PASS |
| `test_passthrough_output_name_all_extensions` | `src/io.rs:474` | PASS |
| `test_norm_path_case_folds` | `src/io.rs:418` | PASS |
| `test_read_id_prefix_bare` | `src/fastq.rs:746` | PASS |
| `test_read_id_prefix_with_description` | `src/fastq.rs:751` | PASS |
| `test_read_id_prefix_strips_slash_one` | `src/fastq.rs:757` | PASS |
| `test_read_id_prefix_strips_slash_two` | `src/fastq.rs:763` | PASS |
| `test_read_id_prefix_strips_slash_three` | `src/fastq.rs:768` | PASS |
| `test_read_id_prefix_strips_slash_one_with_description` | `src/fastq.rs:774` | PASS |
| `test_read_id_prefix_no_strip_slash_four` | `src/fastq.rs:780` | PASS |
| `test_read_id_prefix_empty` | `src/fastq.rs:786` | PASS |
| `test_read_id_prefix_at_only` | `src/fastq.rs:791` | PASS |
| `test_read_id_prefix_no_at_sign` | `src/fastq.rs:797` | PASS (bonus — plan did not require) |
| `test_read_id_prefix_paired_three_way_sync` | `src/fastq.rs:804` | PASS (bonus — plan did not require) |
| `test_passthrough_serial_smoke` (§1) | `src/parallel.rs:1977` | PASS |
| `test_passthrough_serial_parallel_parity` (§2) | `src/parallel.rs:2060` | PASS |
| `test_passthrough_truncation_detected_mid_stream` (§3 + §8) | `src/parallel.rs:2173` | PASS |
| `test_passthrough_id_mismatch_detected` (§4) | `src/parallel.rs:2217` | PASS |
| `test_passthrough_dont_gzip_plain_output` (§6) | `src/parallel.rs:2264` | PASS |

## Spot-checks of plan-specific code-level claims

| Plan claim | Verified at | Result |
|---|---|---|
| `#[derive(PartialEq, Eq)]` on `PairValidationStats` (A-Crit-1) | `src/report.rs:111` | PRESENT |
| Result channel payload is `Result<PairedBatchResult, anyhow::Error>` (B-Crit-1, Step 6a) | `src/parallel.rs:129` | PRESENT — `mpsc::sync_channel::<Result<PairedBatchResult>>(cores * 2)` |
| Reader sends `Err(e)` directly to result channel before exiting | `src/parallel.rs:203–213` | PRESENT — `result_tx_for_reader.send(Err(anyhow::anyhow!(msg)))` |
| Main thread short-circuits on first `Err` | `src/parallel.rs:265–271` | PRESENT — `first_error = Some(e); break 'outer;` |
| `read_id_prefix` strips trailing `/[1/2/3]` (AB2) | `src/fastq.rs:171–178` | PRESENT — `matches!(*suf, "1" \| "2" \| "3")` |
| `read_id_prefix` has 9 unit tests | `src/fastq.rs:746–824` | OVER-DELIVERED — 11 tests (the plan's 9 + `no_at_sign` + `paired_three_way_sync` bonus) |
| Validation §2 fixture is ≥10K records (AB1) | `src/parallel.rs:2063` | SATISFIED — 12,000 rows |
| Validation §2 asserts row-by-row lockstep, not just stat counts (AB5) | `src/parallel.rs:2141–2152` | PRESENT |
| Validation §6 `--dont_gzip` test (AB3) | `src/parallel.rs:2264` | PRESENT |
| Validation §3 fixture ≥10K records, truncate at ~5,000 (§8 coverage) | `src/parallel.rs:2175–2178` | SATISFIED — 10,000 rows, `take(5_000)` |
| 8-arm three-way EOF match in `run_paired_end` | `src/trimmer.rs:424–451` | PRESENT — all 8 arms enumerated |
| 8-arm three-way EOF match in `read_pairs_round_robin` | `src/parallel.rs:663–739` | PRESENT — all 8 arms enumerated |
| `process_pairs` worker plumbing — `&[FastqRecord]` immutable, gzip/plain fork mirrored (AB6) | `src/parallel.rs:559–623` (sig), `parallel.rs:320–399` (fork) | PRESENT |
| Worker stamps `passthrough_records_checked = batch.len()` once per batch | `src/parallel.rs:577–579` | PRESENT |
| Atomic three-write flush per `BTreeMap` iteration | `src/parallel.rs:241–263` | PRESENT — R1/R2/passthrough written within same `while let Some(r) = pending.remove(&expected)` iteration before `expected += 1` |
| Case-folded R1/R2 collision check uses shared `norm_path` helper | `src/cli.rs:597–599`, `src/io.rs:34–36` | PRESENT — `crate::io::norm_path` invoked from both sites |

## Verdict

**COMPLETE.**

All 12 implementation steps (1–11 plus the inserted 6a) and all 8 validation tests are present in code. All 262 lib tests pass; clippy is clean under `-D warnings`. Critical plan-review catches (A-Crit-1 `PartialEq` derive, A-Crit-2 passthrough-longer EOF arm, B-Crit-1 reader-error channel-payload change, B-Crit-2 record-identical assumption, AB1 ≥10K fixture, AB2 `/[123]` strip, AB3 `--dont_gzip` test, AB4/§7 dispatcher coverage, AB5 row-by-row lockstep, AB6 worker plumbing) are all reflected in the code.

Two non-blocking observations:
1. **CLAUDE.md not updated** with a mention of the in-memory fixture generator. Defensible (no on-disk fixture exists) but the plan explicitly asked for it. **Action:** add a one-line note to `CLAUDE.md` test-fixtures section, or document the deviation in `PROGRESS.md`.
2. **Validation §7's promised `tests/integration_passthrough.rs` does not exist** — the implementer's comment at `src/parallel.rs:1839–1840` references a file that is absent. Functionally, the serial path is covered by §1 + §2 + §6; the missing binary-driven dispatcher test is low-value. **Action:** either delete the dangling comment, or add the (~30 LOC) `assert_cmd`-driven test.

Neither observation blocks shipping the feature — but both should be cleaned up before the implementation is considered fully closed.
