# Code Review A — uBAM Input Support

**Reviewer:** A (independent)
**Scope:** `eed17d3` (SE milestone) + uncommitted working tree (PE additions)
**Files audited:** `src/bam.rs`, `src/format.rs`, `src/fastq.rs`, `src/parallel.rs`, `src/trimmer.rs`, `src/cli.rs`, `src/main.rs`, `src/adapter.rs`, `Cargo.toml`.
**Build state:** `cargo check` clean. `cargo test --lib bam format` 14 tests pass. T21–T26 still pending per IMPL.md (expected, not penalised).

## Summary

The feature is well-structured: format detection is content-based and load-bearing-tested against bgzipped-FASTQ confusion; the `RecordSource` trait gives a clean dispatch point; per-record aligned-BAM rejection inside `bam_record_to_fastq` correctly resolves the first-record-only contradiction raised in plan review B. The paired-interleaved de-interleaver design is sound — `MAX_SLACK = 1024` is conservative but enough to catch hand-spliced grouped input cleanly. The committed fixtures are real `samtools import` BAMs; the happy path is genuinely exercised. Test count went from 267 baseline + 14 BAM/format additions, all green.

The main concerns are (i) one race-condition resource leak in the paired de-interleaver, (ii) at least three error paths whose behaviour is asserted in source comments but never exercised by any test (and that will be load-bearing once aligned BAMs land in the wild), and (iii) ~150 lines of duplicated report generation in `run_paired_ubam_single_file` that mirrors `run_paired`'s shape closely enough to invite drift.

## Issues by area

### Logic

**A-H1 (HIGH) — Aligned-BAM detection in paired-interleaved relies on the per-record check fire-and-forget.** In `open_paired_interleaved_with_tags` (`src/bam.rs:153–299`), the producer thread calls `bam_record_to_fastq` per record and propagates errors via `send_pair_err` on both channels. But the producer **discards both errors** with `let _ = tx.send(...)` and returns. The consumer side (`next_record` in `Threaded` mode) at `src/bam.rs:474` does `Ok(Err(e)) => Err(e)` — only the *first* received error is propagated; the second channel still has a queued error that is silently dropped when the reader is dropped. This is acceptable behaviour (one error is enough), but worth documenting because the same record's error is duplicated in both channels and the consumer that wins depends on scheduling. Not a correctness bug; a clarity issue. *(Recommendation: add a comment clarifying this is intentional.)*

**A-H2 (HIGH) — `sanity_check_any` uses a fresh `BamReader::open` (`Direct`) that consumes one record, then is dropped.** `src/main.rs:30–37`. The intent is "peek the first record and validate." But the very next line in main.rs spawns the real reader (`open_threaded_reader` or `open_paired_interleaved_with_tags`) which re-opens the file and re-reads byte 0. This is correct but wastes one BGZF block decompression. More importantly: **for paired uBAM single-file mode, the sanity check sees only R1 of record 0.** If the file is grouped (all R1s followed by all R2s), the sanity check would pass because record 0 is a valid unmapped R1, and the grouped-input error would only fire MAX_SLACK records later inside the producer thread, AFTER the worker pool is spinning. This is consistent with the design (the de-interleaver IS the place that catches this) but the implementation log says the first-record check is the fast path; for paired BAM it gives a false sense of early-exit and the worker pool startup cost is paid before the failure surfaces. *(Recommendation: in main.rs around line 128, when the detected format is `UnalignedBam` and `cli.paired && cli.input.len() == 1`, peek the first TWO records and verify they are mate-adjacent — saves the worker-pool startup cost on the common grouped-input failure mode.)*

**A-M1 (MEDIUM) — `bam_record_to_fastq` does not validate `FPAIRED` (0x1) consistency in paired mode.** `src/bam.rs:507–525`. The de-interleaver routes by FREAD1/FREAD2 only. A BAM record with FREAD1=1 but FPAIRED=0 is malformed (BAM spec requires FPAIRED when FREAD1 or FREAD2 is set) and would be accepted. Real-world impact: low — no standard tool produces this — but a hand-built BAM could slip through. *(Recommendation: add `if !flags.is_segmented() { /* paired */ }` check, or document that we rely on FREAD1/FREAD2 alone. Defer is fine for v1.)*

**A-M2 (MEDIUM) — Empty-seq vs missing-qual interaction.** `src/bam.rs:561–605`. `seq_view.is_empty()` errors. But: BAM records can legitimately have `*` for sequence (encoded as empty in noodles), e.g. when seq was lost in processing but qual is preserved. samtools rejects these for FASTQ conversion the same way. This is consistent — just noting it's correctly rejecting.

**A-M3 (MEDIUM) — `name.is_empty()` is dead code after `name()?`.** `src/bam.rs:527–532`. `noodles::bam::Record::name()` returns `Option<&BStr>` where `None` represents the `*` sentinel; an empty-byte name shouldn't be reachable through the public API. The `is_empty()` check at line 530 cannot fire in practice. Either remove it or document defensively.

**A-L1 (LOW) — `record_idx` mismatch in error vs success path.** `src/bam.rs:421–434`. In the `Direct` mode `next_record`, on success we increment `*record_idx` THEN format the error context using `record_idx` (the field, not the post-increment). On the Err(e) path at line 433 we use `*record_idx + 1`. This is correct because on read-failure we haven't incremented yet. But a casual reader could miss this — recommend renaming `record_idx` → `records_read_ok` or adding a comment.

**A-L2 (LOW) — `paired_records_have_r1_then_r2_order` test name is misleading.** `src/bam.rs:806–819`. The test asserts that R1 and R2 share the *same template name*, NOT that R1 precedes R2 in flag order. The first two records of `ubam_paired_test.bam` are presumably (R1, R2) by `samtools import` convention, but the assertion as written would also pass if both records were R1 of two different templates with the same name. *(Recommendation: assert also that one is FREAD1 and the other FREAD2.)*

### Errors / coverage

**A-H3 (HIGH) — IUPAC, `=`-base, qual-mismatch, and aligned-BAM rejection have no unit-level coverage at all.** Per IMPL.md these were explicitly deferred to "integration level," but `tests/` has no `integration_ubam.rs` (T21 is pending). The result is that four documented error paths (`src/bam.rs:570`, `:573`, `:594`, `:509`) are **completely untested**. The fixtures contain only happy-path records. Risk: a future refactor of `bam_record_to_fastq` that breaks any of these checks lands green. Mitigations would be (a) hand-craft small BAM bytes in a test (~50 lines of noodles writer code, feasible), or (b) commit a tiny "tainted" fixture per case. With only the SE path landing in a beta and PE following soon, this is a near-term hole. *(Recommendation: at minimum add one hand-built-BAM unit test for the aligned-BAM rejection, since it's the most likely real-world failure mode and the validation is per-record not per-file.)*

**A-M4 (MEDIUM) — `Box<dyn RecordSource>` dispatch via BAM is exercised by ZERO tests.** `src/parallel.rs` has `open_fq()` (line 1227) which boxes a `FastqReader`, but no equivalent helper for `BamReader`. Every parallel test goes through the FASTQ path. The trait-object dispatch overhead measurement from Spike 1 is well documented, but **the BAM path through `run_single_end_parallel` and `run_paired_end_parallel` is currently exercised only via the binary, not via any unit/integration test.** The unit-level BAM tests all use `BamReader::next_record` directly (concrete type). A regression that broke `<BamReader as RecordSource>::next_record` would land green at `cargo test --lib`. *(Recommendation: add a `parallel.rs` test that opens a uBAM via `Box::new(BamReader::open_threaded("test_files/ubam_test.bam")?)` and runs `run_single_end_parallel`; takes ~30 lines and locks down the integration point.)*

**A-M5 (MEDIUM) — Channel-error-path edge case in threaded BAM reader.** `src/bam.rs:332–386`. If the producer thread's `tx.send(Err(...))` fails because the receiver was already dropped (e.g., the consumer panicked), the error is `let _`-discarded and the thread silently exits. The receiver side then sees the channel closed and returns `Ok(None)` (line 475) — **effectively swallowing both the read error and the EOF distinction**. In practice the consumer rarely panics, but if it does, the user gets "EOF" instead of the real BAM error. Documented as acceptable in `FastqReader` too (same pattern), so this is parity, not a new bug. Worth a comment.

### Efficiency

**A-L3 (LOW) — `format::open_sync_reader` and `format::open_threaded_reader` duplicate `main.rs::open_sync_reader` / `open_threaded_reader`.** Two copies of the same factory. `src/format.rs:89–119` and `src/main.rs:44–71` are character-for-character equivalent. main.rs imports `InputFormat` and re-implements rather than calling `format::open_*_reader`. *(Fix: delete the main.rs copies and import `format::open_threaded_reader` / `format::open_sync_reader`. Easy DRY win; I have not applied it because it slightly changes main.rs's import surface and the diff is mid-implementation.)*

**A-L4 (LOW) — `detect_input_format` is called THREE times per single-end input on the worker-pool path.** Once in `sanity_check_any` (main.rs:128), once in the input_formats vec (line 137), once inside `open_threaded_reader` (line 48). Each call opens the file, peeks 4 bytes, opens it again to MultiGzDecoder, decompresses one block. Cost is microseconds, but the redundancy is ugly. *(Recommendation for follow-up: detect once at the top of main, pass the format through.)*

**A-L5 (LOW) — `spawn_noop_joinhandle()` spawns two real OS threads per paired uBAM open.** `src/bam.rs:486–488`. Each `BamReader` from `open_paired_interleaved_with_tags` carries a `_handle: JoinHandle<()>` whose only purpose is to satisfy the struct shape; the actual producer is detached. The cost is OS thread setup/teardown (~10–100µs on Linux/macOS, negligible). No reaper/zombie concern — these are pthreads, not processes; they exit and are joined when the BamReader drops. Slightly cleaner: refactor `BamReaderSource::Threaded._handle` to `Option<JoinHandle<()>>` so the paired-interleaved readers can pass `None` and skip the noop spawn entirely. *(Cosmetic; no fix needed.)*

### Structure / duplication

**A-H4 (HIGH) — `run_paired_ubam_single_file` duplicates `run_paired`'s report block at ~150 lines.** `src/main.rs:1182–1405`. The setup (output-path derivation, eprintln-trim banner, parallel-vs-serial branch) is genuinely different (single BAM stem vs two FASTQ stems), but the report generation from line 1324 down is nearly identical to `run_paired`'s lines 1056–1147 — same `report::TrimConfig` construction, same `for (idx, ...)` loop, same `write_pair_validation_stats` placement, same JSON shape. Drift risk is real: a future change to `run_paired`'s report (e.g., adding a new stat) will almost certainly miss this function. *(Recommendation: extract a `write_paired_reports(input_filenames: &[Path], stems: &[Path], stats_r1, stats_r2, pair_stats, config_for_reports) -> Result<()>` helper that both call sites use. ~100 lines of net reduction. Defer to a post-feature DRY pass per the review brief, but flag it now so it doesn't slip.)*

**A-M6 (MEDIUM) — Tag preservation tests are absent.** `src/bam.rs::tests` has no test for `with_preserved_tags`, `open_threaded_with_tags`, or the tag-formatting in `append_tag_type_and_value`. Per IMPL.md the golden test is T23. But the `parse_sam_tag_name` validator at the CLI layer is also untested in `src/cli.rs::tests` (or at least I see no test exercising it via clap's `try_parse_from`). At minimum a unit test that exercises `parse_sam_tag_name("CB")` (pass), `parse_sam_tag_name("X")` (fail length), `parse_sam_tag_name("9X")` (fail leading alpha), and `parse_sam_tag_name("ALL")` (reserved) would lock the validator. ~10 lines.

**A-L6 (LOW) — Error message inconsistency.** `src/bam.rs:142` uses `bail!("Input file not found: {}", ...)` matching `FastqReader::open_threaded`. But `src/bam.rs:326` (`open_threaded_inner`) uses the same string. Three identical bails. Could be a helper.

## Fixes applied

None — all findings are either documentation/coverage gaps (better surfaced in a report than committed silently) or architectural recommendations that touch multiple call sites. The build is clean and the existing tests pass; nothing here blocks landing.

## Recommendations with priority

1. **(CRITICAL — before T21 lands)** Add a hand-built-BAM unit test for the aligned-BAM rejection. Single most likely real-world failure mode; the per-record check is the only line standing between users and silent corruption-by-soft-clip. ~30 lines using noodles' `bam::io::Writer`.

2. **(HIGH)** Add a `parallel.rs` test that boxes a `BamReader` and runs it through `run_single_end_parallel`. Locks down the trait-object dispatch path for BAM, which is currently zero-tested at unit level (A-M4).

3. **(HIGH)** Refactor `run_paired_ubam_single_file`'s report block into a shared helper with `run_paired`. Drift risk is the main concern (A-H4); the duplication is acceptable for one release but should not survive into v2.3.

4. **(MEDIUM)** Add fast-path two-record peek in main.rs for paired-uBAM single-file input (A-H2). Catches the common "grouped input" mistake before the worker pool spins up.

5. **(MEDIUM)** Add unit tests for `parse_sam_tag_name` validator (A-M6). ~10 lines, locks the CLI surface.

6. **(MEDIUM)** Pull the duplicated `open_*_reader` factories in main.rs and call `format::open_*_reader` instead (A-L3). One-line replacement per call site.

7. **(LOW)** Document the "first error wins" semantics in `open_paired_interleaved_with_tags` (A-H1), and clarify the `record_idx` accounting comment in `BamReader::next_record` Direct path (A-L1).

## Verdict

**The implementation is sound and matches the plan.** Format detection, the trait-object refactor, and per-record validation are all well-executed. The hard contract — bgzipped FASTQ vs unaligned BAM disambiguation — is correctly enforced via payload-byte check.

**The coverage gaps are real and worth closing before the integration-test pass.** Specifically: zero error-path unit tests for `bam_record_to_fastq`, zero trait-object tests for the BAM dispatch path, and zero tests for `parse_sam_tag_name`. Each is a ~10–30 line fix; together they would catch ~four real regression classes.

**Drift risk in `run_paired_ubam_single_file` should be tracked as a follow-up** but is not blocking.

Report file: `/Users/fkrueger/Github/TrimGalore/plans/06252026_ubam-input-support/CODE_REVIEW_A.md`
