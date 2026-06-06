# Code Review A — `--passthrough` mode

**Reviewer:** Reviewer A (independent, fresh context)
**Date:** 2026-06-06
**Scope:** all `src/` changes vs `master`, plan v2 (`PLAN.md`) is the spec.
**Tests:** `cargo test --lib` → 262 pass, 0 fail. `cargo clippy --all-targets --release -- -D warnings` → clean. `cargo fmt --all -- --check` → clean. The 5 new `parallel::tests::test_passthrough_*` tests all pass.

---

## Summary

The implementation is materially faithful to plan v2. The biggest delta — the `Result<PairedBatchResult>` channel for fail-fast reader-error propagation (Step 6a) — works correctly for the reader-error path that the tests exercise. The 8-arm three-way EOF match is exhaustive and behaves sensibly. The serial/parallel parity test (§2) is the load-bearing invariant and it passes on a 12K-record fixture that genuinely forces multi-batch ordering.

That said, the new error-handling architecture has **one untested edge case that can deadlock** (worker-side error while the reader is still producing), one piece of **dead code** that should be removed, and two **divergence/coverage issues** that aren't blocking but should be acknowledged in follow-up work.

The plan's `u64` → implementation's `usize` deviation for the new stats counters is documented inline and is the right call — it matches the surrounding fields in `PairValidationStats`. The case-folded R1/R2 collision test is weak (only exercises byte-equal aliasing) but the case-fold helper itself is well-covered by `io::tests::test_norm_path_case_folds`.

---

## Issues by area

### Logic

#### L1 (High) — Potential deadlock when a **worker** errors mid-stream

**File:** `src/parallel.rs:128–292`

The plan's Step 6a fix correctly addresses the *reader*-error case: the reader bails, returns from `read_pairs_round_robin`, sends `Err` on `result_tx_for_reader`, and the closure scope drops `txs`, which unblocks every worker. Test §3 (`test_passthrough_truncation_detected_mid_stream`) exercises this path and passes without hanging.

The **worker-error case is not handled symmetrically.** When a worker errors:

1. Worker sends `Err` on `rtx`, breaks out of `while let Ok(...)`, scope drops its own `rtx` and `rx`.
2. Main thread receives `Err`, sets `first_error`, hits `break 'outer` (line 270).
3. Main thread is no longer draining `result_rx`.
4. Reader thread is **still running** — still producing batches via `txs[idx].send(...)`.
5. Other workers continue processing and start blocking on `rtx.send(...)` once `result_rx` capacity (`cores * 2`) fills.
6. With worker channels each capped at 2 (`mpsc::sync_channel::<PairedWork>(2)` at line 120), the reader's next sends to live workers will eventually block too.
7. `reader_handle.join()` at line 278 then waits forever.

This is a real deadlock — not just a slow exit. It is not covered by any test in the new suite. The reader-error case the tests *do* cover is structurally different because the reader returns by itself before `join()` is called.

**Severity rationale:** worker errors in `process_paired_batch` are rare (gzip-encoder write failures on a `Vec<u8>` sink essentially never happen; the only realistic failure is OOM during compression, which is exceptional). But "exceptional" ≠ "impossible," and the failure mode is "process hangs" rather than "process exits with error" — which makes diagnosis on a user machine miserable.

**Recommendation:** after `break 'outer`, drain `result_rx` to completion (with `try_recv` in a loop, plus `recv_timeout` shutdown polling) before calling `reader_handle.join()`. Or have the main thread close the work channels first by signalling the reader to stop. The cleanest fix is probably to wrap the `'outer` loop body so that after `first_error` is set, the loop *continues* draining `result_rx` (discarding any further `Ok`/`Err`) until the channel closes — that lets the reader and workers all exit cleanly, then join.

#### L2 (Low) — `passthrough_records_checked` not stamped on empty batches

**File:** `src/parallel.rs:577`

```rust
if reads_passthrough.is_some() {
    pair_stats.passthrough_records_checked += reads_r1.len();
}
```

If a passthrough-active batch has `reads_r1.len() == 0` (theoretically possible at exact-multiple-of-`BATCH_SIZE` EOF), this is a no-op. The reader-side flushing logic at line 670 and 709 only sends partial batches when `!batch_r1.is_empty()`, so this can't actually happen in practice. Safe, but worth a `debug_assert!(reads_r1.len() == reads_passthrough.as_ref().map_or(0, |s| s.len()))` to catch future drift.

#### L3 (Medium) — Serial/parallel `r1_unpaired`/`r2_unpaired` divergence is masked, not fixed

**File:** `src/parallel.rs:524–535` (parallel `classify_paired`) vs `src/trimmer.rs:579–595` (serial `run_paired_end`)

The parallel path increments `r1_unpaired`/`r2_unpaired` whenever the filter says `r{1,2}_ok = true`, regardless of whether the unpaired writer is open. The serial path only increments when the writer is also open. With `--retain_unpaired` rejected under `--passthrough`, neither writer is open in the new use case — yet the parallel path still increments the counter and the serial path doesn't.

The fixture works around this by making both R1 and R2 mates equally short (`gen_multiome_rows` comment at `parallel.rs::tests::gen_multiome_rows`), so `r{1,2}_ok` is always `false` for dropped pairs. That's a deliberate fixture choice, and the comment is honest about it.

**This is pre-existing — not introduced by this PR.** But it does mean the §2 parity test is narrower than it claims: it asserts `PairValidationStats` equality, but only for fixtures where the divergence happens to be invisible. A different `--length`/`--length_1`/`--length_2` configuration would still trip it.

**Recommendation (non-blocking):** open a follow-up issue to fix the underlying divergence — `classify_paired` should not increment `r{1,2}_unpaired` when no unpaired writer is open. The fix is local to `classify_paired` and would let the parity test be widened beyond the current narrow fixture. Acceptable to defer to a follow-up PR.

#### L4 (Low) — 8-arm three-way match is exhaustive but order-sensitive

**Files:** `src/trimmer.rs:424–451`, `src/parallel.rs:654–739`

The serial-path version writes:

```rust
match (rec1, rec2, rec_pt) {
    (r1, r2, None) => (r1, r2, None),   // passthrough disabled fallthrough
    (Some(r1), Some(r2), Some(Some(pt))) => ...
    (None, None, Some(None)) => ...
    (Some(_), None, _) => bail!(...),    // R2 truncation — wildcard, fires for any pt state
    (None, Some(_), _) => bail!(...),
    (Some(_), Some(_), Some(None)) => bail!(...),
    (None, None, Some(Some(_))) => bail!(...),
}
```

The wildcards `_` on the third position in the R1/R2-truncation arms work because the more-specific arms (the `(r1, r2, None)` fallthrough and the happy/EOF arms) are listed earlier. Match-arm ordering matters here in a non-obvious way — flipping the order would change behavior (e.g., putting `(Some(_), None, _)` first would swallow the `(Some(_), None, None)` no-passthrough path through the same arm, which is actually still correct because the error message handles both cases generically). Still, a comment block at the top of the match enumerating the (passthrough-disabled, happy, EOF, error) categories would help future maintenance.

**Verdict:** correct as written. Cosmetic only.

#### L5 (Info) — `(Some, Some, None)` is treated as passthrough-disabled, not as an error

The serial-path catch-all `(r1, r2, None) => (r1, r2, None)` fires when `rec_pt` is the bare `None` — which only happens when `reader_pt.is_none()` (passthrough disabled). If passthrough is enabled, `rec_pt = Some(reader.next_record()?)` and is either `Some(Some(rec))` or `Some(None)`. So this arm is unambiguous. The invariant could be expressed as a `debug_assert!(reader_pt.is_some() == matches!(rec_pt, Some(_)))` at the top of the match, but it's small enough that the current shape is fine.

### Efficiency

#### E1 (Low, cleanup) — Dead `passthrough_active` local in `run_paired_end_parallel`

**File:** `src/parallel.rs:114`, `286`, `289`

```rust
let passthrough_active = input_passthrough.is_some();
...
let _ = passthrough_active;  // line 286
...
let _ = passthrough_active;  // line 289
```

Lines 286 and 289 are `let _ = ...` suppressions for an unused variable — clippy will warn without them. But `passthrough_active` is never read anywhere in the function body (the inner `process_paired_batch` has its own local of the same name, line 313). This is leftover scaffolding from an earlier draft.

**Fix (apply directly):** delete line 114 (`let passthrough_active = ...`) and the two `let _ = passthrough_active;` suppressions at 286/289.

#### E2 (Info) — Per-batch `pair_stats.passthrough_records_checked += reads_r1.len()` is correct

**File:** `src/parallel.rs:577`

The plan questioned whether per-batch accounting matches per-record. It does: the reader's happy path arm pushes one record to each batch in lockstep, so `reads_r1.len() == reads_passthrough.unwrap().len()` is invariant. Verified by tracing `read_pairs_round_robin` — all three batch vecs are pushed together inside the `(Some, Some, Some(Some))` arm, and the EOF/error arms send before reaching the BATCH_SIZE flush. The §2 parity test (12K records, 3 batches) is the operational confirmation.

Could be hardened with a `debug_assert_eq!(reads_r1.len(), s.len())` after `let s = reads_passthrough.unwrap()`, but the invariant is structural.

### Errors / safety

#### Err1 (Medium) — Reader-error path uses string-round-trip via `anyhow::anyhow!(format!("{e:#}"))`

**File:** `src/parallel.rs:207–208`

```rust
let msg = format!("{e:#}");
let _ = result_tx_for_reader.send(Err(anyhow::anyhow!(msg)));
```

The reader formats its error to a string, then constructs a new `anyhow::Error` from the string. This loses the error chain (source/cause) — downstream `e.chain()` calls won't recover the original. For a simple `bail!("Read ID mismatch...")` that's a single-link chain anyway, so the user-visible message is preserved. But if the reader bubbles an I/O error from `FastqReader::next_record()?`, the original `io::Error` source is gone.

The reason it's structured this way is that `anyhow::Error` is `!Sync` — sending it through `mpsc::SyncSender<Result<...>>` requires a `Send + Sync` payload, and `anyhow::Error` is `Send` but not `Sync` (because the inner trait object is `dyn Error + Send + 'static`, which is not `Sync`).

Actually wait — `mpsc::SyncSender<T>` requires `T: Send`, not `T: Send + Sync`. Let me double-check.

After verification: `mpsc::SyncSender<T>` requires `T: Send`. `anyhow::Error` IS `Send`. So the channel can send `anyhow::Error` directly. The string round-trip is unnecessary.

**Fix recommendation:** change

```rust
let msg = format!("{e:#}");
let _ = result_tx_for_reader.send(Err(anyhow::anyhow!(msg)));
```

to (consume `res` by moving):

```rust
if res.is_err() {
    // Take the error out of `res`, send it on the result channel, and
    // recompute `res` for the return path.
    let e = res.unwrap_err();
    let msg = format!("{e:#}");
    let _ = result_tx_for_reader.send(Err(e));
    res = Err(anyhow::anyhow!(msg));
}
```

…or restructure so the error is consumed and a `bail!`-equivalent message is returned. The current shape is *correct* in behavior (the user-visible message survives), it just unnecessarily flattens the error chain.

**Severity:** Medium because it's a quality issue, not a correctness issue. The test §3 passes because the assertion is `err.contains("passthrough") && err.contains("truncated")` — message text only. If a downstream caller ever wanted to `.downcast_ref::<io::Error>()` on a reader-side IO failure, they wouldn't recover it. Today no caller does that.

#### Err2 (Low) — `eprintln!("Worker error: ...")` removed; worker errors now silent except via channel

**File:** `src/parallel.rs:155–161`

Pre-PR, worker errors printed to stderr. Now they only surface via the result channel. If the result-channel `send` fails (very unlikely, but possible if main has already broken out), `let _ = rtx.send(Err(e))` silently swallows the error. The pre-PR `eprintln!` was a belt-and-suspenders backup; removing it makes worker errors invisible if both the channel send fails AND the main thread is in a state where it can't surface them.

In practice this is a non-issue (the channel send won't fail because main thread is blocked in `recv` until something is sent), but if L1 above is fixed by draining `result_rx` after `break 'outer`, the worker errors get drained and re-checked. A defensive `eprintln!` before the `rtx.send` would cost nothing.

#### Err3 (Info) — Partial-output cleanup deferred per plan §13 — acceptable

Documented in the help text and in `PLAN.md` §Assumptions §13. Fine.

### Structure / style

#### S1 (Low) — `test_passthrough_rejects_demux` weakens its own assertion

**File:** `src/cli.rs:1325–1345`

```rust
let err = cli.validate().unwrap_err().to_string();
assert!(
    err.contains("--passthrough requires --paired") || err.contains("--demux"),
    "got: {err}"
);
```

The test acknowledges via comments that it can't drive the `--demux` rejection arm cleanly (clap-level constraints) and so the assertion accepts either message. This is honest, but it means the `if self.demux.is_some()` rejection at `cli.rs:594` is **not actually exercised** by any test. If someone deletes that arm tomorrow, no test fails.

**Recommendation:** convert this into a test that builds a `Cli` struct directly (bypassing `Cli::parse_from`) with `paired = true`, `demux = Some(path)`, `passthrough = Some(path)` — then call `validate()` and assert on the `--demux`-specific message. Or accept the gap (the arm is one line; low risk).

#### S2 (Low) — Stats field type deviation from plan documented inline

**File:** `src/report.rs:124–137`

Plan said `u64`, implementation used `usize` to match existing siblings (`pairs_analyzed`, `pairs_removed`, etc., all `usize`). The deviation is documented in a comment block. This is the right call — codebase convention beats plan literalism here. Accept.

#### S3 (Low) — Plan said `format!("{}{}", stem, ext)` but implementation matches; no issue

Verified: `passthrough_output_name` body in `io.rs:191–204` matches the plan's signature spec verbatim. No drift.

#### S4 (Low) — Variable naming `reader_pt` vs `reader_passthrough`

The serial path uses `reader_pt` (line 405) and the parallel reader thread uses both `reader_pt` (line 191) and `reader_passthrough` (the parameter, line 637). Minor inconsistency, but the local-vs-parameter rename is sensible because the parameter is `Option<&mut FastqReader>` while the local is owned by the closure. Acceptable.

#### S5 (Low) — Test fixture cleanup elides on failure

All parallel tests do `fs::remove_dir_all(&dir).ok();` *before* `Ok(())`. If an assertion fails, the dir isn't cleaned. This matches the convention of surrounding tests in `parallel.rs` (see `test_parallel_serial_trim_stats_parity`), so no change needed — but worth noting that CI tmp dirs accumulate on failure. The `fresh_tmpdir` helper presumably handles this on next run.

---

## Fixes applied

None — all findings are recommendations rather than direct fixes. The user's harness explicitly opted for review-then-implement. Two changes that are safe to apply directly (E1 dead code, Err1 error-chain) would be one-line edits but are left for the user to authorize per the global "implementation requires explicit trigger" rule.

---

## Recommendations with priority

### Critical
*(none)*

### High
- **L1** — Drain `result_rx` after `break 'outer` to prevent deadlock on worker-side errors. Add a regression test that simulates a worker error mid-stream (e.g., by injecting a bad gzip level or an OOM-like failure into one worker's batch — feasible via a test-only feature flag).

### Medium
- **L3** — Track the pre-existing serial/parallel divergence on `r{1,2}_unpaired` accounting in a follow-up issue. The current fixture-side workaround is honest but narrow. Fix in `classify_paired` to mirror the serial path (only count when the writer is open).
- **Err1** — Drop the `format!("{e:#}")` string round-trip in the reader-error propagation; send the `anyhow::Error` directly (the channel supports `Send`-only payloads). Preserves the error chain for downstream `.chain()` consumers.

### Low
- **E1** — Delete dead `passthrough_active` local + suppression in `run_paired_end_parallel` (lines 114, 286, 289).
- **L2** — Add `debug_assert_eq!(reads_r1.len(), reads_passthrough.as_ref().map_or(reads_r1.len(), |s| s.len()))` at the top of `process_pairs` to catch future drift between R1 batch length and passthrough batch length.
- **S1** — Strengthen `test_passthrough_rejects_demux` so the `--demux` rejection arm in `Cli::validate()` is actually exercised (build `Cli` directly, bypass clap).
- **Err2** — Restore the `eprintln!("Worker error: ...")` belt-and-suspenders before the channel `send`, so worker errors stay visible even if the channel send is dropped.
- **L4** — Add a category-comment block above the 8-arm match (both serial and parallel) — listing (disabled / happy / EOF / R1-truncation / R2-truncation / pt-short / pt-long) — to make the order-sensitivity self-documenting.

### Info / non-blocking
- **§7 integration test deferred:** acceptable scope cut. The §1 smoke test + §2 parity test cover the underlying library paths; `main.rs::run_paired` is a thin wrapper around them. A 5-line `tests/integration_passthrough.rs` that invokes the built binary would still be valuable as a smoke check (and would catch wiring regressions in `main.rs`), but its absence isn't a blocker.
- **S2 `u64` → `usize` deviation:** documented inline, matches codebase convention. Accept.

---

## Verdict

**Approve with a High-severity follow-up** (L1 deadlock window). The PR is shippable for the documented v1 use case (reader-side sync errors, which the tests cover comprehensively), but the worker-error deadlock should be fixed before this is exercised in a busy production workflow. Everything else is cleanup or coverage-expansion work that can land in a follow-up PR.

The plan-v2 → implementation correspondence is high: every plan-step has an implementation site, the dual plan-review fixes (AB1–AB6, A-Crit-1, A-Crit-2, B-Crit-1, B-Crit-2) are all visible in the diff, and the 5 new parallel-path tests are well-shaped. The implementation team picked up the right conventions (e.g., `usize` over plan-specified `u64`) and documented deviations inline.
