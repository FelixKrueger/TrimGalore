# CODE_REVIEW_B — uBAM input support (#316)

Reviewer B, independent review (no coordination with Reviewer A).
Scope: `git diff master...HEAD` (commit `eed17d3` SE milestone) + uncommitted PE additions on top.
Lens: concurrency correctness, build/dep surface, test-coverage gaps, byte-identity invariant, CLI ergonomics, reproducibility.

## Summary

The implementation is in good shape against the v3 plan. `cargo build --release`, `cargo test --lib` (281 passing), and `cargo clippy --all-targets --release -- -D warnings` all succeed. The dependency surface is minimal and clean: `noodles 0.88.0 (features=["bam"], default-features=false)` shares the exact transitive tree that fastqc-rust 1.0.1 already pulls; `cargo tree` confirms NO `tokio` is pulled, and `find ~/.cargo/registry/src` confirms **no `build.rs`** in any of `noodles-bam-0.73.0`, `noodles-bgzf-0.35.0`, `noodles-sam-0.69.0`, `noodles-core-0.16.0`, or `noodles-csi-0.42.0` — so no new build-time non-determinism risk. Release binary is 7.5 MB.

The concurrency model in `BamReader::open_paired_interleaved` is sound: lockstep matched-pair emission + capacity-equal mpsc channels + `is_err() || is_err()` short-circuit prevents the classic asymmetric-channel deadlock. Byte-identity to v0.6.11 (CI `validation` job) is preserved — format detection is a side branch; FASTQ flow is unchanged through trait-object dispatch (Spike-1 cost ~1%, no semantic change).

The biggest concrete gap is **test coverage for error/edge paths in `bam_record_to_fastq` and the de-interleaver overflow** — explicitly noted as deferred per IMPL.md "code-implementation agent notes" (use committed fixture for happy paths; cover error paths at integration). T21–T26 (integration tests, CI job, repro check, docs) are explicitly pending and not blockers for this review.

## Issues

### Logic / Concurrency

**B1 — [Medium] Detached paired-uBAM producer thread loses panics silently.**
`src/bam.rs:153,319` — `_producer_handle: JoinHandle<()>` is dropped at function end, detaching the thread. The doc-comment at line 149-152 says this is intentional. Confirmed: send-Err paths surface to consumers correctly via the `is_err() || is_err()` short-circuit (line 270-273). However, a `panic!` inside `bam_record_to_fastq` or noodles parsing would terminate the producer without surfacing to either consumer — they would see `rx.recv()` close (returning `None`) and treat it as EOF, producing a silently-truncated output. This is unlikely with noodles' lazy decode (panics would be bugs, not user-data issues), but it deviates from the rest of the codebase's `reader_handle.join()` discipline in `src/parallel.rs:292-296`. Recommendation: store the producer `JoinHandle` in a `OnceLock<JoinHandle>` accessible to both readers, and join it on the second-to-drop. Or accept the limitation and add a `std::panic::catch_unwind` wrapper inside the producer closure that surfaces panic-as-error via `send_pair_err`. Not blocking for v1; flag for follow-up.

**B2 — [Low] `spawn_noop_joinhandle()` leaks two kernel threads per paired-uBAM open.**
`src/bam.rs:486-488` — Each call to `open_paired_interleaved_with_tags` spawns 2 throwaway threads as `_handle` field placeholders. Trivial code smell; cleaner fix is `_handle: Option<JoinHandle<()>>` in `BamReaderSource::Threaded` (set to `None` for the paired-side readers). Not user-visible.

**B3 — [Low] MAX_SLACK threshold message off-by-one.**
`src/bam.rs:279, 43-48` — The check triggers when `r1_queue.len() > MAX_SLACK` (= 1025), but the error message says "exceeded MAX_SLACK = 1024 records." Either change to `>=` or update the message to "queue depth 1025". Cosmetic.

**B4 — [Medium] Adapter auto-detect on paired-uBAM scans interleaved R1+R2 records.**
`src/main.rs:261` (paired-uBAM single-file path) → `setup_trimming(&cli, &cli.input[0])` → `adapter::autodetect_adapter(input)` → `crate::format::open_sync_reader(path, &[])` which opens BAM as plain SE (not de-interleaved). So adapter detection samples both R1 and R2 sequences mixed together. Illumina R2 carries the SAME adapter as R1 in standard chemistry (`AGATCGGAAGAGC` at 3'), so signal is preserved — but for `--small_rna` / asymmetric chemistries this could dilute or bias the auto-detect. Not a correctness bug for the common case; worth a code-comment noting the implication or restricting auto-detect on paired uBAM to a R1-only filter.

### Build / Dependency Surface

**B5 — [Positive finding, no fix needed] Dep tree is clean.**
- `cargo tree -e normal -f "{p} {f}"` confirms `noodles 0.88.0` is **already pulled transitively** via fastqc-rust 1.0.1 at the exact same version (and noodles-bam 0.73.0, noodles-bgzf 0.35.0, noodles-sam 0.69.0, noodles-core 0.16.0, noodles-csi 0.42.0 — all reused, no parallel-version compile).
- `grep -c tokio Cargo.lock` = **0** — no async runtime pulled.
- `find ... -name build.rs` in all 5 noodles crates = **empty** — no proc-macro / codegen / wall-clock risk to `SOURCE_DATE_EPOCH` reproducibility. The existing CI `reproducibility` job's invariants are intact.

### Test Coverage Gaps

**B6 — [Medium] `bam_record_to_fastq` error paths are untested in `cargo test`.**
The committed fixture `test_files/ubam_test.bam` is correctly all-unmapped, well-formed, single-end ACGTN — so the unit tests cover happy paths only. Per `src/bam.rs:821-827` comment this is explicit deferral to integration tests (T21–T23, pending). The error branches with NO test coverage today:
1. Aligned-BAM rejection per-record (line 509-515) — the v3 plan B-Crit-4 fix.
2. Reverse-complemented / secondary / supplementary rejection (517-524).
3. Empty seq (line 562).
4. `=` rejection (line 570).
5. IUPAC coerce-to-N + `emit_iupac_warning_once` (line 573-585).
6. Qual length mismatch (line 594).
7. Unknown base byte (line 577).
8. Empty read name (line 529, 531).
9. Tag value: array-tag bail (line 644-647), hex-tag round-trip (637-640), all integer types.

Recommendation: add hand-crafted BAM bytes (or vendor a tiny `samtools view -h | samtools view -bS` round-trip into the build script) for these. Won't block the v1 merge if T21/T22 integration tests land before release.

**B7 — [Medium] `parse_sam_tag_name` value parser has zero unit tests.**
`src/cli.rs:386-402` — the function rejects "ALL", non-2-char, and non-`[A-Za-z][A-Za-z0-9]` patterns. Reviewed by reading; logic correct. But unit tests for the four rejection paths AND the happy path (e.g. `"CB"`, `"UB"`, `"A1"`) are entirely missing. Should live alongside the existing tag-related tests in `mod tests`. ~10 LOC.

**B8 — [Medium] `MAX_SLACK` overflow path has no test.**
The `GROUPED_INPUT_ERR` branch (line 279-282) is dead code in the test suite — the committed paired fixture is strictly mate-adjacent (line 751-756 doc-comment confirms it). A synthesised grouped BAM (R1×N then R2×N, N > 1024) would exercise this. Without it, regression in the bounded-buffer guardrail is invisible.

**B9 — [Medium] No test for `Cli::validate` rule `--paired SINGLE.bam`.**
`src/cli.rs:759-769` validates the structural-only case (any single file). The actual format-detection-gated check is in `src/main.rs:155-161`. No integration test fires that error path with a single FASTQ + `--paired`. Could be a 5-line test if `main.rs` were factored to expose the format check; otherwise an integration test covers it.

**B10 — [Low] `--preserve-tags` clap-level rejection of `ALL` has no end-to-end test.**
`Cli::try_parse_from(["...", "--preserve-tags", "ALL"])` should error. Currently inferred via reading the value-parser; would be 3 LOC to test.

### CLI Ergonomics

**B11 — [Low] `--preserve-tags` error messages are slightly inconsistent.**
- `parse_sam_tag_name` uses single quotes: `'{s}' is not a valid SAM tag name`.
- The other `Cli::validate` messages in the file use varied quoting: backticks (`--clumpify`), no-quotes (`--paired`).
Not a functional bug; cosmetic consistency would be nice.

**B12 — [Low] Help-text says `--preserve-tags` is "Ignored for FASTQ input" but main.rs:143-148 emits a WARNING.**
Mild contradiction. The warning is appropriate (user-error visibility), but the help should say "Emits a warning if no uBAM input is present" rather than "Ignored." Tiny doc nit.

### Compatibility with FASTQ paths (byte-identity)

**B13 — [Positive finding, no fix needed] FASTQ-in → FASTQ-out byte path is preserved.**
The architectural change moves all readers (FASTQ and BAM) through `Box<dyn RecordSource>` (Spike-1: ~+1% overhead — accepted). The FASTQ reader internals (`src/fastq.rs:230+`) are unchanged. The worker pool body in `parallel.rs` is reader-agnostic; the trait-object hop happens BEFORE the channel pipeline, not inside the hot loop. Output bytes for a FASTQ input MUST be identical to pre-PR. Spot-check: `./target/release/trim_galore test_files/BS-seq_10K_R1.fastq.gz` produces `BS-seq_10K_R1_trimmed.fq.gz` with md5 `f8d879e31d9b797667796cec0f208dd5` — matches the historic byte-output expectation (one run; CI's `validation` job will confirm).

### Reproducibility

**B14 — [Positive finding, no fix needed] No new build-time non-determinism.**
Verified the noodles tree has no `build.rs` files (5 crates checked). `SOURCE_DATE_EPOCH`-driven reproducibility is intact. The existing CI `reproducibility` job should continue to pass with no changes. T25 (the explicit repro check) is correctly listed as pending in IMPL.md.

## Fixes applied

None. All findings are either positive-confirmations (B5, B13, B14), low-severity polish (B2, B3, B10, B11, B12), or test-coverage gaps explicitly deferred to T21–T26 (B6, B8, B9). The two medium-severity gaps that ARE actionable in this PR (B7 — `parse_sam_tag_name` tests; B11 — help text) are recommendations, not requirements for v1 ship.

## Recommendations with priority

| # | Finding | Priority | Type | Where |
|---|---|---|---|---|
| B1 | Detached producer panic loss | Medium | Concurrency | `src/bam.rs:153,319` |
| B2 | `spawn_noop_joinhandle` smell | Low | Refactor | `src/bam.rs:486-488` |
| B3 | MAX_SLACK off-by-one in error | Low | Polish | `src/bam.rs:43-48,279` |
| B4 | Auto-detect scans both R1+R2 on paired uBAM | Medium | Logic | `src/main.rs:261` + `adapter.rs` |
| B6 | `bam_record_to_fastq` error-path tests | Medium | Coverage | `src/bam.rs:tests` |
| B7 | `parse_sam_tag_name` unit tests | Medium | Coverage | `src/cli.rs:tests` |
| B8 | MAX_SLACK overflow test | Medium | Coverage | `src/bam.rs:tests` |
| B9 | `--paired SINGLE.bam` integration test | Medium | Coverage | `tests/integration_ubam.rs` (pending T22) |
| B10 | `--preserve-tags ALL` clap-rejection test | Low | Coverage | `src/cli.rs:tests` |
| B11 | Error-message quoting consistency | Low | Polish | `src/cli.rs:386-402` |
| B12 | "Ignored for FASTQ" doc-vs-runtime mismatch | Low | Doc | `src/cli.rs:184-188` |

## Verdict

**Ship-ready for the v1 uBAM milestone, modulo T21–T26 landing.** No critical findings. Concurrency reasoning holds up. Byte-identity invariants are preserved. Dep tree is clean. Recommend B6/B7/B8 as the smallest unit-test additions that would significantly tighten the safety net before release.
