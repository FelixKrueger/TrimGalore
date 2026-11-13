# Plan Coverage Report

**Mode:** B (code vs. plan)
**Plan(s):** `plans/11132026_issue-542-bgzf-writer/PLAN.md`
**Date:** 2026-11-13
**Verdict:** INCOMPLETE
4 items unresolved — all of them validation/test-coverage items. Every code change in §2, §4, §5 and §6 landed.

**Input deviation:** no `IMPL.md` exists for this work (`implementation-planner` was never run), so the
audit ledger is taken from the design plan itself: §2 (file-change table + seven doc-comment sites),
§5 (behaviour and edge cases), §6 (seven-step outline), §9 (validations 9.1–9.7) and §10 (assumptions).
§14's claims were verified against the code rather than accepted.

**What was audited:** uncommitted working-tree changes on `dev`, base `d862f7d` —
`src/bam.rs` (+303), `src/parallel.rs` (+27/-), `src/main.rs` (-3 args, 1 comment),
`.github/workflows/ci.yml` (+28), `CHANGELOG.md` and `docs/src/content/docs/reference/changelog.md`
(+6 each).

## Summary

- Total items: 44
- DONE: 38
- PARTIAL: 1
- MISSING: 2
- DEVIATED: 3

## Coverage ledger

### §2 — file changes

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `src/bam.rs`: `BgzfSink<W>`, `ParallelBamSink` reduced to a wrapper, `workers` dropped, `BgzfBlocks` + `into_bgzf` | §2 table | DONE | `src/bam.rs:1105-1250` |
| 2 | `src/parallel.rs`: batches compressed in the workers, two call sites take the newtype | §2 table | DONE | `:345`, `:455`, `:568`, `:609`, `:421`, `:524` |
| 3 | `src/main.rs`: three `create` calls lose `cli.cores`; stale thread comment rewritten | §2 table | DONE | `:2381`, `:2555`, `:2701`; comment `:2698` |

### §2 — the seven doc-comment sites (and the two to leave alone)

All seven rewritten, each to current state; none narrates the change. Verified against the diff text,
not against §14's claim.

| # | Site (pre-change line) | Was | Now | Status |
|---|---|---|---|---|
| 4 | `bam.rs:1109-1110` struct doc | "compressed to BGZF by a pool of threads" | "batches arrive already encoded and compressed, and are appended in order" | DONE |
| 5 | `bam.rs:1112-1113` field-order rationale | "so the BGZF stream is closed before the temporary is removed" | "so the writer is flushed before the temporary is removed on the abandon path" | DONE |
| 6 | `bam.rs:1120` `create` doc | "`workers` is the number of compression threads" | "Same header rules as [`BamWriter::create`]." | DONE |
| 7 | `bam.rs:1156` `finish` doc | "Flush, stop the compression threads, write the BGZF EOF block…" | "Terminate the BGZF stream, then publish the file under its final name." | DONE |
| 8 | `parallel.rs:354-355` PE doc | "the sink compresses them to BGZF on its own pool of `cores` threads" | "trim, filter, encode and BGZF-compress … the main thread appends" | DONE |
| 9 | `parallel.rs:460-462` SE twin | same claim | same rewrite | DONE |
| 10 | `main.rs:2700` | "N workers trim and encode, N threads compress" | "N workers trim, encode and compress" | DONE |
| 11 | `parallel.rs:74` (FASTQ architecture block) — **must stay** | — | untouched (no diff hunk; `src/parallel.rs:70-80` still describes the FASTQ worker/gzip design) | DONE |
| 12 | `bam.rs:1853` (reader pools) — **must stay** | — | untouched; now at `src/bam.rs:1936`, still the `pooled_readers_match_plain_readers` doc | DONE |

Residual-reference check: `grep -rn "write_encoded\|MultithreadedWriter\|multithreaded_writer" src/ tests/ .github/`
returns nothing outside the changelog prose. No stale name survives in code.

### §4 — signatures

| # | Item | Status | Notes |
|---|------|--------|-------|
| 13 | `BgzfBlocks` newtype with `len`/`is_empty` | DONE | `bam.rs:1120-1135`; `as_bytes` added beyond the plan (§14 deviation 5, declared). `len()` has no caller — declared in §4, so DONE as specified. |
| 14 | `BatchEncoder::into_bgzf` | DONE | `bam.rs:1109` |
| 15 | `BatchEncoder::into_bytes` retained | DONE | `bam.rs:1104`; used by 9.1 |
| 16 | `BgzfSink<W: Write>` — `new`/`header`/`write_compressed`/`finish` | DONE | `bam.rs:1159-1195`, all four present with the planned semantics |
| 17 | `ParallelBamSink::create` without `workers`; `header`; `write_compressed`; `finish` | DONE | `bam.rs:1205-1250`. `inner` is `Option<BgzfSink<…>>` with two `expect`s (§14 deviation 4, declared and justified — `finish` consumes by value). |

### §5 — behaviour

| # | Item | Status | Notes |
|---|------|--------|-------|
| 18 | Create: `PendingOutput` → `BufWriter<File>` → header as its own BGZF member | DONE | `BgzfSink::new` encodes through a throwaway `bam::io::Writer<Vec<u8>>` then `bgzf_members` |
| 19 | Per batch in the worker: `into_bgzf()` (write → flush → `into_inner`) | DONE | `bgzf_members`, `bam.rs:1138-1146` |
| 20 | Per batch on the main thread: verbatim append, reorder buffer unchanged | DONE | `parallel.rs:421`, `:524`; no ordering code added |
| 21 | Finish order: terminator → `BufWriter::flush` → `commit` | DONE | `BgzfSink::finish` `bam.rs:1189-1194` then `ParallelBamSink::finish` `:1240-1250`. The ordering §12 flagged is correct. |
| 22 | Error: `Err` up the existing `first_error` path; `PendingOutput` still `Some`, so `Drop` removes the temporary | DONE | unchanged loop structure; `pending` only taken in `finish` |

### §5 — edge cases

| # | Case | Status | Notes |
|---|------|--------|-------|
| 23 | **Empty batch** (every read filtered) | MISSING (test) | The behaviour is implemented, but nothing asserts that `into_bgzf` on an empty encoder yields an empty `Vec` and no stray member. The plan's shape choice (`flush` + `into_inner` over `try_finish`) was decided *by* this case (§12), and it is the one edge case with no coverage and no declared waiver. |
| 24 | **Zero records overall** | MISSING (test) | No test drives the parallel uBAM path with empty input. `test_parallel_empty_input_no_deadlock` (`parallel.rs:1905`) is the FASTQ path; `tests/integration_ubam_out.rs` has no empty-input case. The 102 tests there all assert ≥1 record. |
| 25 | **Batch > 64 KiB** (multi-member batches) | DONE | Covered by the existing parity tests (10–12 k records, `parallel.rs:2860`/`:2936`) and `ubam_out_cores_gt_1_matches_single_core_across_batches` (12 000 records, `assert_ubam_eq` at `--cores 1` vs `4`), plus the new write-error test's 2 000-record batch. |
| 26 | **Worker fails mid-batch** | DONE | Unchanged by design, and already guarded end-to-end by `ubam_out_cores_gt_1_failure_leaves_existing_output_untouched` (`tests/integration_ubam_out.rs:369`), which fails in the second batch and asserts the pre-existing file survives. |
| 27 | **Write error (`ENOSPC`/`EPIPE`)** | DONE | `write_error_is_returned_not_panicked` |
| 28 | **`--cores 1`** | DONE | All three dispatch sites still gate on `cli.cores > 1` (`main.rs:2379`, `:2552`, `:2697`); the new CI step runs `--cores 1` as the comparison arm. |

### §6 — implementation outline

| # | Step | Status | Notes |
|---|------|--------|-------|
| 29 | 1. Two helpers, `bgzf_members` and `bgzf_terminator`, one postcondition each | DONE | `bam.rs:1137-1153`; doc comments state the postconditions ("Never emits the end-of-file block" / "and nothing else"). One extra inline comment inside `bgzf_members` explains the `flush`/`into_inner` mechanics that step 1 said belong in the commit message — noted, not scored. |
| 30 | 2. `BatchEncoder::into_bgzf` | DONE | — |
| 31 | 3. `BgzfSink` + wrapper; `NonZeroUsize` dance gone; field-order comment rewritten | DONE | the `NonZeroUsize`/`set_worker_count` block is removed entirely |
| 32 | 4. Workers compress; both result structs carry `BgzfBlocks` | DONE | — |
| 33 | 5. Three `main.rs` call sites + the comment | DONE | — |
| 34 | 6. Existing tests updated for the dropped argument | DONE | `bam.rs:2816` (was `:2733`), `parallel.rs:2860` (was `:2861`), `:2936` (was `:2937`) |
| 35 | 7. Changelog: one sentence under `### Unreleased`, naming the panic fix and the N+1 thread reduction | DONE | `CHANGELOG.md:8-13` and the docs copy at `docs/src/content/docs/reference/changelog.md:14-19`, identical text, both under `### Unreleased` → `#### Changes`. Names the write-error fix, the double-panic abort, and 29 → 20 threads. |

### §9 — validations

| # | Validation | Status | Notes |
|---|------------|--------|-------|
| 36 | **9.1** batch round-trips (the flush-omission mode) | DONE | `into_bgzf_round_trips_a_small_batch` (`bam.rs:3015`) — 3 records, positive round-trip via `bgzf::Reader`, plus an explicit non-empty assertion. Built as specified, not as the defective first draft. **PASSES.** |
| 37 | **9.1a** no stray terminator mid-stream | DONE | `into_bgzf_emits_no_terminator` (`bam.rs:3047`) — `windows()` scan for the 28-byte block **and** a non-empty assertion, so it cannot pass vacuously. **PASSES.** |
| 38 | **9.2** output readable and record-identical to serial | DONE | The two inflate-and-compare parity tests needed only the dropped argument; no "exactly one EOF block" assertion was added to the serial side, as §9.2 required (§3.5). The new terminator-count test is pointed at `ParallelBamSink` only and says so in its doc comment. **PASSES.** |
| 39 | **9.2a** an htslib reads the new framing (CI) | DEVIATED | Step added in the `validation-ubam` job (`ci.yml:1280-1307`), which installs samtools and builds `--release`. Two deviations: (a) `samtools quickcheck` dropped — **declared** (§14 deviation 1) and correct, it exits 8 on any unaligned BAM; (b) **undeclared** — the fixture is `BS-seq_10K_R1.fastq.gz` decompressed, not the `test_files/ubam_test.bam` the plan named. That substitution makes the step *stronger*: `ubam_test.bam` is 507 bytes (~10 records), one sub-64-KiB batch, so it could not have exercised multi-member concatenation; 10 000 reads gives ~3 batches. The step asserts `samtools view -c` > 0, equal record counts and equal record-body md5 across `--cores 1`/`4`, and exactly one terminator at EOF. The step runs FASTQ-in → uBAM-out rather than uBAM-in → uBAM-out; the sink path is identical either way. |
| 40 | **9.3** write failure reports an error and publishes nothing | DEVIATED | First half exactly as planned: `write_error_is_returned_not_panicked` (`bam.rs:3112`) injects `FailingWriter` (budget 4096, then `io::Error::other`) into `BgzfSink<W>` and asserts the error propagates with the underlying message. **PASSES.** Second half **undeclared deviation**: the plan said to drive `ParallelBamSink` against a path under a read-only directory; the code instead has `a_failed_run_publishes_nothing` (`bam.rs:3135`), which drops the sink without `finish`. That exercises the abandon path of §5.5 rather than a write failure, so it is a different mechanism — but the end-to-end case the plan was reaching for is already covered by the pre-existing `ubam_out_cores_gt_1_failure_leaves_existing_output_untouched`. Substance covered; method changed without a note. |
| 41 | **9.4** file-size cost ≤ 2 % | DEVIATED (declared) | §14 reports **+0.045 %** (44 869 908 → 44 890 204 bytes, 500 k reads), far inside the 2 % threshold the plan set up front. Measured serial-vs-parallel on one build rather than before-vs-after binaries — **declared** (§14 deviation 2) with a sound rationale: `--cores 1` vs `--cores 8` isolates the per-batch block-boundary cost directly. Not independently re-measured in this audit. |
| 42 | **9.5** throughput cost | PARTIAL → **not run** | **Declared** (§14 deviation 3): the volume is at 100 % capacity with 5.1 GiB free, and this repo has had a benchmark corrupted by disk exhaustion before. The reasoning for not running it is sound, but the plan made this a validation and §7 took a deliberate decision ("a regression up to ~25 % at `--cores 8` is accepted") that nothing has yet checked. Outstanding. |
| 43 | **9.6** thread count fell by exactly N+1 | PARTIAL | The *after* figure is measured and matches the prediction exactly: **20 peak threads** at `--cores 8` SE, 62 samples, = 2N+4. The *before* figure (29) is taken from #542 rather than re-measured, which §9.6 specified as "measured, not derived" — though the same row also calls the before-figure a cross-check against #542's number, so this is a partial reading rather than a contradiction. Disclosed in §14's measurement table but not listed among the deviations. |
| 44 | **9.7** whole suite, fmt, clippy | DONE | **Independently re-verified in this audit** (not taken from §14): `cargo test` → **881 passed, 0 failed, 0 ignored** across 20 binaries, exit 0 — matching §14's 876 + 5 exactly. `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release -- -D warnings` exit 0 with an empty log (fully cached, so the cache matches the current sources — cargo would have recompiled otherwise). |

### §10 — assumptions

| # | Assumption | Status | Notes |
|---|---|---|---|
| — | Compression level stays the noodles default (6) on both sides | DONE | No `set_compression_level` anywhere on the BGZF path; `grep` finds levels only in `fastq.rs` (gzip) and `clump_only.rs`. |
| — | `--cores 1` never reaches this path | DONE | verified at all three dispatch sites (item 28) |
| — | `clump_only` and the serial `BamWriter` untouched | DONE | no diff hunk in `src/clump_only.rs`; `BamWriter` unchanged |

### §14 — the implementer's claims, checked

| Claim | Verdict |
|---|---|
| 881 tests pass; fmt clean; clippy clean | **Confirmed independently** (item 44) |
| Four mutation tests, each killing one new test | Not re-run (the mutations were reverted; `src/bam.rs` is byte-identical to the state audited). The four tests exist, assert what the table says they assert, and each one's assertion points *against* the mutation it names — the failure mode that this repo's memory flags (a mutation pointing the same way as its assertion) is not present here. |
| Six deliberate deviations | 4 and 5 verified in code (items 13, 17); 1 and 2 verified as declared (items 39, 41); 3 verified as declared and outstanding (item 42); 6 verified site-by-site (items 4–12). |
| 9.5 not run | Confirmed, and the reason stands. |
| Two further deviations **not** on that list | 9.2a's fixture substitution (item 39) and 9.3's second-half method change (item 40). Both are benign — one strengthens the test, one is already covered elsewhere — but neither was declared. |

## Gaps (detail)

### Item 23: empty-batch round-trip has no test

**Expected:** §5's first edge case — `into_bgzf` on an empty encoder yields an empty `Vec`, `write_compressed` writes 0 bytes, no stray member. §12 records that this case is what decided the `flush()`+`into_inner()` shape over `try_finish`.
**Found:** The behaviour is implemented and is correct by inspection, but no test pins it. A future change to `bgzf_members` that emitted a member or a terminator on empty input would pass the whole suite.
**Gap:** One assertion — `BatchEncoder::with_capacity(…).into_bgzf()?.is_empty()`.

### Item 24: zero-records-overall has no test on this path

**Expected:** §5's second edge case — header member + EOF block, read back as a valid 0-record BAM.
**Found:** No empty-input test for `--output-format ubam --cores > 1`. `test_parallel_empty_input_no_deadlock` covers the FASTQ path only; none of the 102 tests in `tests/integration_ubam_out.rs` uses empty input.
**Gap:** Either a unit test over `ParallelBamSink` with no batches, or an integration run on an empty FASTQ.

### Item 42: 9.5 (throughput) not run

**Expected:** `--cores 4` and `8`, three runs each, wall-clock, on a machine with ≥2N free cores; a regression beyond ~25 % to be brought back before merging.
**Found:** Not run. Declared, with a defensible reason (volume at 100 %, 5.1 GiB free, and a prior benchmark in this repo corrupted by disk exhaustion).
**Gap:** Free space, rebuild both binaries, md5 them, then measure. §7 accepted a real regression as the price of honouring `--cores N`; the size of it is still unknown.

### Item 43: 9.6's before-figure not re-measured

**Expected:** "measured, not derived" — 29 → 20, with the before-figure cross-checking the method against #542.
**Found:** 20 measured over 62 samples. 29 taken from #542.
**Gap:** One sampling run on a `d862f7d` binary, if the cross-check is still wanted. Low value — the after-figure matching 2N+4 to the thread already confirms the model.

## Test verification (Mode B)

| Test | File | Status |
|------|------|--------|
| `into_bgzf_round_trips_a_small_batch` | `src/bam.rs:3015` | **PASS** |
| `into_bgzf_emits_no_terminator` | `src/bam.rs:3047` | **PASS** |
| `parallel_sink_terminates_the_file_exactly_once` | `src/bam.rs:3067` | **PASS** |
| `write_error_is_returned_not_panicked` | `src/bam.rs:3112` | **PASS** |
| `a_failed_run_publishes_nothing` | `src/bam.rs:3135` | **PASS** |
| `batch_encoder_subtracts_input_phred_offset` (updated for the dropped arg) | `src/bam.rs:2816` | **PASS** |
| PE serial-vs-parallel uBAM parity (updated) | `src/parallel.rs:2860` | **PASS** |
| SE serial-vs-parallel uBAM parity (updated) | `src/parallel.rs:2936` | **PASS** |
| `ubam_out_cores_gt_1_matches_single_core_across_batches` | `tests/integration_ubam_out.rs:328` | **PASS** |
| `ubam_out_cores_gt_1_failure_leaves_existing_output_untouched` | `tests/integration_ubam_out.rs:369` | **PASS** |
| Whole suite | — | **881 passed, 0 failed, 0 ignored** |
| `cargo fmt --all -- --check` | — | clean |
| `cargo clippy --all-targets --release -- -D warnings` | — | exit 0, no warnings |
| CI step "Validate parallel uBAM-out — htslib reads the worker-compressed framing" | `.github/workflows/ci.yml:1280` | present in `validation-ubam`; not executed in this audit (needs samtools + CI runner) |
| Empty-batch assertion | — | **MISSING** |
| Zero-record parallel uBAM output | — | **MISSING** |

## Verdict

**Verdict:** INCOMPLETE — 4 items unresolved.

No code gap. Every file change, signature, behaviour step and doc-comment site the plan specified is
present, and the two sites §2 said to leave alone are untouched. The whole suite, fmt and clippy were
re-verified independently rather than taken from §14.

What remains:

1. **§9.5 (throughput) has not been run** — declared, with a sound reason (no disk). §7 accepted a
   real regression by decision; its magnitude is still unmeasured.
2. **§5's empty-batch case has no test** — the case that decided the `flush`+`into_inner` shape.
3. **§5's zero-records case has no test** on the parallel uBAM path.
4. **§9.6's before-figure (29 threads) was not re-measured** — taken from #542. Low value now that the
   after-figure matches 2N+4 exactly.

Two further deviations were found that §14 does not declare, both benign and neither a gap: §9.2a runs
on a 10 000-read FASTQ instead of the 507-byte `ubam_test.bam` the plan named (stronger — the named
fixture is a single sub-64-KiB batch and could not have tested multi-member concatenation), and
§9.3's "publishes nothing" half is tested by dropping the sink without `finish` rather than through a
read-only directory (the end-to-end case is already covered by a pre-existing integration test).
