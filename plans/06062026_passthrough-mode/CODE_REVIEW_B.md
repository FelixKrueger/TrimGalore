# Code review — `--passthrough` mode (Reviewer B)

**Reviewer:** B (independent of A)
**Date:** 2026-06-06
**Scope:** Working-tree diff against `master` for the `--passthrough` feature.
**PLAN:** `/Users/fkrueger/Github/TrimGalore/plans/06062026_passthrough-mode/PLAN.md` (v2)
**Build state:** `cargo test --lib` → 262 passed / 0 failed. `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release -- -D warnings` clean. The 5 new parallel-path tests (`parallel::tests::test_passthrough_*`) all pass.

## Summary

The implementation is solid and faithful to the v2 plan. The B-Crit-1 reader-error-propagation harness is in place and exercised by `test_passthrough_truncation_detected_mid_stream`; the 8-arm three-way EOF match is genuinely exhaustive in both serial and parallel paths (Rust's exhaustiveness check enforces it, and I walked all 12 inhabitants by hand); `read_id_prefix` strips `/[123]` per AB2; `PartialEq, Eq` on `PairValidationStats` is in place; the JSON-emitter `pair_validation` trailing-comma edit is correct in every reachable mode. Diff is large (~1.6k LOC added) but mostly mechanical, and the parity test now uses a 12 000-record fixture that genuinely spans three batches at `BATCH_SIZE = 4096`.

The most significant gap is **not a bug in the new code** — it's an acknowledged pre-existing parallel/serial divergence on `r1_unpaired`/`r2_unpaired` accounting that the new parity test sidesteps via fixture design (both mates short, so neither side is rescue-eligible). This is documented in a comment on `gen_multiome_rows` (`src/parallel.rs:1860–1865`) but the underlying bug is left for a follow-up. Several smaller issues below: a deferred CLI-binary integration test that was named in plan §7 but never landed, two `let _ = passthrough_active;` shadow-suppression lines that look like dead code from a refactor, and a string-roundtrip in the reader-error propagation path that loses the `anyhow` error chain.

## Issues by area

### Logic

**L1 — `r1_unpaired` / `r2_unpaired` parallel/serial divergence (Medium; pre-existing, not introduced by this PR).**
Location: `src/parallel.rs:524–534` (`classify_paired`) vs `src/trimmer.rs:580–591` (`run_paired_end`).

The parallel `classify_paired` unconditionally bumps `pair_stats.r1_unpaired += 1` / `r2_unpaired += 1` whenever the filter returns `TooShort { r1_ok: true, r2_ok: _ }` — independent of whether `--retain_unpaired` opened any unpaired writer. The serial path only increments these counters inside the `if let Some(w) = unpaired_r1 && r1_ok` blocks — so when `--retain_unpaired` is off, the serial counters stay at 0 and the parallel counters can be non-zero. This is the bug `gen_multiome_rows` works around by making **both** mates short (so `r1_ok == r2_ok == false`).

Status: acknowledged in the implementation (`src/parallel.rs:1858–1865` comment); explicitly deferred. I'd call this **acceptable scope-cut for the passthrough PR** because (a) the bug already lives on `master`, (b) the parity test is honest about the workaround, and (c) fixing it is a one-line gate in `classify_paired` but would need fresh test coverage to land safely. **Recommendation:** file a follow-up issue immediately so it doesn't get forgotten — and consider making `gen_multiome_rows` accept a config knob so a future fix can flip a single boolean to exercise the case.

**L2 — Missing `tests/integration_passthrough.rs` (Medium).**
Plan §Validation §7 (AB4) required a `--cores 1` dispatcher-level CLI-integration test that runs `main.rs::run_paired` end-to-end via `std::process::Command` on the built binary, to confirm the routing decision and the wiring at `src/main.rs:738–787` (output path computation, sanity check, FastQC sweep) works. The comment at `src/parallel.rs:1841–1843` claims this test "lives in `tests/integration_passthrough.rs`" — but the file does not exist (I checked the working tree).

Consequence: no end-to-end CLI exercise of:
- the `cli.passthrough.is_some()` branch in the output-collision pre-flight (`src/main.rs:201–211`)
- `passthrough_output_name` integration with `output_dir` + `basename` at the call site
- the `pair_idx == 0` gating of the passthrough sanity-check (line 248)
- the `eprintln!("  Passthrough: ...")` user-visible logging
- the report generation `write_passthrough_stats` call (`src/main.rs:980–987`)
- the FastQC sweep extension (`src/main.rs:1035–1041`)

The unit-level tests cover the serial and parallel libraries but bypass `main.rs::run_paired` entirely. **Recommendation:** before merging, either (a) add the integration test the plan promised, or (b) rewrite the misleading comment to accurately state what's missing and link the follow-up issue.

**L3 — `pair_idx == 0` guard is dead-code-by-validation (Low).**
`src/main.rs:248`:
```rust
if pair_idx == 0 && let Some(ref pt_path) = cli.passthrough {
    FastqReader::sanity_check(pt_path)?;
}
```
`Cli::validate()` step 1.ii rejects `--passthrough` with `input.len() != 2`, so the paired-end dispatcher's outer loop only iterates once and `pair_idx` is always 0. The check is harmless defence-in-depth but rather than a `pair_idx == 0` filter, a `debug_assert!(pair_idx == 0)` would express the invariant correctly and crash loudly if the validation envelope ever relaxes to multi-pair without a corresponding wiring update.

**L4 — Serial-path passthrough drop accounting in the `(None, None, _)`-without-passthrough EOF (Low).**
`src/trimmer.rs:430` correctly maps `(None, None, Some(None))` to `(None, None, None)` and then breaks at line 607. But `(None, None, None)` (no-passthrough EOF case) also maps to that same tuple via the wildcard match at line 426. Both paths break cleanly. Walked through and the 12-arm enumeration is sound.

### Efficiency

**E1 — `read_id_prefix` UTF-8 boundary work on every record (Low — performance acceptable, but I want to flag the cost in the parallel reader).**
`src/fastq.rs:175`:
```rust
let head = stripped.split_ascii_whitespace().next().unwrap_or("");
head.rsplit_once('/')
    .filter(|(_, suf)| matches!(*suf, "1" | "2" | "3"))
    .map(|(p, _)| p)
    .unwrap_or(head)
```
`split_ascii_whitespace` walks the string char-by-char on Unicode boundary checks, even though ASCII whitespace is what we want. On 1M paired records this is 3 × 1M ≈ 3M scans. Plan §Efficiency says ~90 ms / 1M reads which is fine, but a `head.bytes().position(|b| b == b' ' || b == b'\t')` would shave ~30% off the inner cost. Not a blocker, but consider it for v2 if the reader thread becomes the bottleneck at high `--cores`. Treat as a deferred optimisation.

**E2 — Two `let _ = passthrough_active;` shadow-suppression lines (Low / cleanup).**
`src/parallel.rs:286, 289`:
```rust
if let Some(e) = first_error.or(reader_join_err) {
    let _ = passthrough_active;
    return Err(e);
}
let _ = passthrough_active;
```
`passthrough_active` at line 114 is bound but never read — these `let _ = ...` lines look like leftovers from a refactor where the variable was used for output-naming or pre-flight checks before that logic moved into `read_pairs_round_robin`. Either remove the line-114 binding and the two suppressions, or actually use the variable (e.g., to gate a debug log on entry/exit). Compiler doesn't complain because the suppressions silence the unused-binding warning.

### Errors

**E-Err-1 — Reader-thread error round-trips through `String`, losing the anyhow chain (Medium).**
`src/parallel.rs:206–209`:
```rust
let res = inner();
if let Err(ref e) = res {
    let msg = format!("{e:#}");
    let _ = result_tx_for_reader.send(Err(anyhow::anyhow!(msg)));
}
```
`format!("{e:#}")` produces a multi-line string with the full error chain, but wrapping it back into `anyhow::anyhow!(msg)` discards the original error type — any downstream `e.downcast_ref::<SomeErr>()` introspection breaks, and the chain becomes a single flat string.

Mitigating factor: the **same** `res` is returned by the reader closure (line 213) and reaches the main thread via `reader_handle.join()` at line 278–282. So the main thread still has access to the unflattened error on the `reader_join_err` path — `first_error` only catches the flattened one when the result-channel arm fires first. In practice, for any real reader error, the main loop drains a result-channel `Err` first, sets `first_error`, then `first_error.or(reader_join_err)` picks the flattened one. The unflattened reader error is silently dropped.

For human-readable error messages this is fine (the flattened form contains everything). For programmatic error handling (which Trim Galore doesn't really need today), the chain is lost. **Recommendation:** clone the `anyhow::Error` via `Err(anyhow::anyhow!("{:#}", e))` is equivalent — but to preserve the chain, use a wrapper that supports `Clone` (e.g. `Arc<anyhow::Error>`) or rely entirely on the join handle's return value. As-is, the priority order `first_error.or(reader_join_err)` should probably be inverted: prefer the **unflattened** reader error from the join handle over the flattened one from the channel, since they encode the same root cause.

**E-Err-2 — Partial output file cleanup is explicitly v1-deferred — documented (Acceptable).**
Per plan §Assumption §13, mid-stream reader errors leave `*_val_1.fq.gz` / `*_val_2.fq.gz` / `*_passthrough.fq.gz` partials on disk. The `--passthrough` help text mentions this. Acceptable per plan, but worth a runtime warning to the user something like:

> `--passthrough` failed mid-stream — partial output files at `<path>` may exist; re-run after fixing inputs.

Right now the user just gets the bail message without the cleanup hint. **Recommendation:** in `run_paired` (`src/main.rs:738+`), wrap the `parallel::run_paired_end_parallel` call's `Err` arm with an `eprintln!` that lists the partial output paths.

### Structure

**S1 — `usize` vs `u64` convention follow-through (Praise).**
The plan §Signature called for `u64` on the three new counters; the implementer correctly identified that the existing `PairValidationStats` fields are `usize` and matched the existing convention (`src/report.rs:127–138` with an explicit note explaining the choice). Good call — the plan's `u64` was an oversight.

**S2 — `passthrough_records_checked` accounting consistency (Confirmed).**
- Serial (`src/trimmer.rs:470`): per-record increment, gated on `pt_record.is_some()`.
- Parallel (`src/parallel.rs:578`): once-per-batch via `pair_stats.passthrough_records_checked += reads_r1.len();`, gated on `reads_passthrough.is_some()`.

The parity test (`test_passthrough_serial_parallel_parity`) asserts these stats match — I confirmed the test passes with the 12 000-record fixture. The contract that the reader thread has already per-record-verified sync before dispatching each batch holds: `read_pairs_round_robin` only pushes to `batch_pt` after a successful sync check (`src/parallel.rs:691–706`), and on EOF (line 670, line 712) the dispatched `batch_r1.len()` equals the actual record count. **Accounting is consistent.**

**S3 — JSON emitter's `pair_validation` trailing-comma adjustment (Praise).**
`src/report.rs:1002` changes `}` to `},` and line 1005 changes `json_null` `comma=false` to `comma=true`. The new `passthrough` block at line 1014 is now the last field. I traced all four reachable cases — SE/no-pt, SE/pt (impossible per validation), PE/no-pt, PE/pt — and all produce valid JSON. The existing `serde_json::from_slice` deserialization test would have caught a malformed comma, and `cargo test report::tests::test_write_json_report_paired_end_includes_pair_validation` passes.

**S4 — `passthrough_active` flag duplication (Low / cleanup).**
The flag is computed three times:
- `src/parallel.rs:114` (`run_paired_end_parallel`)
- `src/parallel.rs:316` (`process_paired_batch`, from `reads_passthrough.is_some()`)
- `src/parallel.rs:648` (`read_pairs_round_robin`, from `reader_passthrough.is_some()`)

Each is locally derived from a different source-of-truth (input_passthrough, reads_passthrough, reader_passthrough). Could be hoisted into a single computed value passed by argument, but the current shape is defensible — each function owns its own logic. Cosmetic.

**S5 — Specialty `test_passthrough_rejects_demux` is testing the wrong thing (Low).**
`src/cli.rs` test (line ≈1320): comment admits "we hit 1.i first; the 1.i error is sufficient evidence the validation chain runs." That's not actually exercising the `--demux` rejection — it's exercising the `--paired` precondition. A real test for the `--demux` rejection would need `--demux <samplesheet>` + `--paired` + `--passthrough`, but per the test comment, `--demux` already conflicts with `--paired` at validation. **Recommendation:** delete the test (it's misleading) or restructure to demonstrate that `--demux` + `--passthrough` would fail (it does, via the existing demux/paired check). The current test name is a lie.

**S6 — Hardcoded test fixture path in CLI tests (Low).**
`src/cli.rs:1184`: `const PT: &str = "test_files/SRR24766921_RRBS_R2.fastq.gz";`. Reuses an RRBS test fixture as a stand-in passthrough file. Fine for unit tests where we don't actually consume the file, but if `test_files/SRR24766921_RRBS_R2.fastq.gz` is ever renamed (it has been before; see the git history mentioning `Phase-1B Perl-parity hunt`), these 13 tests all break with a confusing "file not found" error. Consider a generic `test_files/passthrough_test_input.fq.gz` symlinked to whatever fixture is convenient.

## Fixes applied

None applied directly — all findings are either documented decisions (S1, S2, S3 = praise; E-Err-2 = plan-accepted), pre-existing bugs (L1), deferred work (L2), or low-priority cleanups requiring judgement (E1, E2, S4, S5, S6) where I prefer to flag for the author rather than touch the diff.

## Recommendations with priority

### Critical
None. The implementation is correct on the load-bearing paths (sync check, error propagation, ordered flush, stats parity), all 262 tests pass, and the diff is `cargo fmt` / `cargo clippy` clean.

### High
- **L2** — File the `tests/integration_passthrough.rs` test that plan §7 promised, **or** edit the misleading comment at `src/parallel.rs:1841–1843` so it doesn't claim coverage that doesn't exist. Picking (b) over (a) is fine for this PR if a follow-up issue is filed.

### Medium
- **L1** — File a follow-up issue for the pre-existing parallel/serial `r1_unpaired`/`r2_unpaired` divergence. The fixture workaround in `gen_multiome_rows` is honest but it's a footgun for anyone extending the parity test later.
- **E-Err-1** — Invert the priority in `first_error.or(reader_join_err)` to prefer the unflattened reader-handle error, **or** stop round-tripping through `String` in the reader-error path. Both are short edits.
- **E-Err-2** — Add a user-facing partial-output warning to the failure path so users know to clean up `*_val_1` / `*_passthrough` etc. before re-running.

### Low
- **L3** — Make the `pair_idx == 0` gate a `debug_assert!`.
- **E1** — Replace `split_ascii_whitespace().next()` with a byte-level scan in `read_id_prefix` if the reader thread becomes the parallel-path bottleneck.
- **E2** — Remove the two `let _ = passthrough_active;` shadow-suppressions and the now-unused binding at `src/parallel.rs:114`.
- **S4** — Consolidate the three `passthrough_active` computations into a single arg-passed flag (cosmetic).
- **S5** — Fix the `test_passthrough_rejects_demux` test to either really test demux+passthrough or be removed.
- **S6** — Decouple the hardcoded fixture path in the CLI tests from a specific RRBS fixture.

## Verdict

**APPROVE with high-priority follow-up.** The new code is correct, tested, idiomatic, and faithful to the v2 plan. The two issues I'd want addressed before tag/release (rather than as a follow-up issue) are L2 (the missing integration test, or at least the false comment) and E-Err-1 (the error-chain flattening). Everything else can ship as low/medium follow-up work.
