# Code Review — `--clump_only` (Reviewer B)

**Target:** `feature/clump-only` (off `dev`)
**Scope:** `src/clump_only.rs`, `src/cli.rs`, `src/main.rs`, `src/io.rs`, `src/lib.rs`, `tests/integration_clump_only.rs`, `.github/workflows/ci.yml`, `docs/src/content/docs/modes/clump-only.md`, `CHANGELOG.md`
**Date:** 2026-07-25

---

## 1. Verdict

**APPROVE WITH REVISIONS.**

The core byte-identity path is correct: streaming reader → per-bin dispatch → stable, tiebreaker-cascaded sort → RFC-1952 concatenated gzip members. Empty-input edge case is handled (`src/clump_only.rs:334-336, 488-493`). Rejection matrix in `Cli::validate()` covers every load-bearing flag from PLAN §Rejection matrix. No cross-cutting impact on `--clumpify`, the trim pipeline, or Perl-parity CI (verified by inspection: no shared mutable state, no changes to sort primitives, no touch on `trimmer.rs` / `parallel.rs`).

The revisions below are test-coverage and behavior-parity gaps, not correctness holes in the shipping code path.

---

## 2. Critical findings

### 2.1 [Medium] `--no_report_file` is silently ignored on the `--clump_only` path
- **Evidence:** `src/main.rs:326-388` never inspects `cli.no_report_file`; `clump_only_single`/`_paired` unconditionally call `write_clump_only_report` (`src/clump_only.rs:356-358, 511-524`). Every other trim/specialty path guards report writing with `if !cli.no_report_file` (`src/main.rs:890, 1209, 1429, 1647, 1787, 1902`).
- **Impact:** User passes `--clump_only --no_report_file` → clumping report file is still written. Inconsistent with the rest of the tool; users piping into archival workflows will see an unexpected sidecar file.
- **Fix:** Thread `no_report: bool` into `clump_only_single`/`_paired`, or gate the `write_clump_only_report` calls in `main.rs`. Add a regression test.

### 2.2 [Medium] Two PLAN-mandated unit tests missing; deviation undocumented
- **Evidence:** PLAN §8 lists `test_clump_only_normalizes_plus_line` and `test_clump_only_normalizes_crlf` — neither exists (`src/clump_only.rs:717-1027`; grep `test_clump_only_` returns 8 tests only). The contract-scope note is prominent in docs and `--help`, but the codebase-wide normalizations it depends on (plus-line → bare `+`, CRLF → LF) are unpinned for this mode.
- **Impact:** A future refactor of `FastqReader`/`FastqWriter` that changed either behaviour would break `--clump_only`'s documented contract without failing any test. The PROGRESS.md deviations log does not mention the omission.
- **Fix:** Add both tests (fixtures can be synthesized in-memory via `write_synthetic_fastq` variant that emits `+HEADER` or `\r\n`); or document in PROGRESS.md why they were dropped.

### 2.3 [Medium] uBAM-input rejection has no integration test despite PLAN §8 calling it out
- **Evidence:** `src/clump_only.rs:1029-1033` explicitly comments "covered by integration tests" but `tests/integration_clump_only.rs` contains no such test. PLAN §8 line 295: "Assert exit != 0 with a clear 'uBAM input is not yet supported under --clump_only' message." An existing uBAM fixture exists (`test_files/` has `.bam` fixtures per repo layout).
- **Impact:** The `reject_ubam` path (`src/clump_only.rs:233-243`) is untested end-to-end. Silent regression risk if `detect_input_format` semantics change.
- **Fix:** Add one integration test spawning the binary against a uBAM fixture, asserting non-zero exit and error-message substring.

### 2.4 [Low] Silent-accept regression guard covers only `-q`; not `--stringency` or `-e`
- **Evidence:** `tests/integration_clump_only.rs:331-360` tests `-q 30` only. PLAN §Silently ignored (line 122-124) lists all three flags with identical rationale.
- **Impact:** Regression against `--stringency` or `-e` accidentally being wired into a code path on the `--clump_only` dispatch would slip past CI.
- **Fix:** Extend `silently_accepts_quality_flag` to a parameterized helper or add two more tests.

### 2.5 [Low] Peak memory doubles briefly per bin flush under `gzip=true`
- **Evidence:** `src/clump_only.rs:207-215` allocates a full `Vec::new()` for the gzip encoder output, then `out.write_all(&member)`. For a bin at the `--memory`-scaled budget (multi-hundred-MB), the peak instantaneously carries `records + member` before the member drops.
- **Impact:** Not a correctness bug; only inefficient. Users on tight RSS budgets may see peak overshoot vs. `--clumpify`'s worker pool which streams into a shared writer.
- **Fix (optional):** Wrap `out` in a byte-counting writer and construct `GzEncoder::new(&mut counter, ...)`; call `.finish()?` and read `counter.written`. Slightly more plumbing; not required for v1.

### 2.6 [Low] Rejection-matrix has no per-flag unit tests in `cli.rs`
- **Evidence:** PLAN §Validation item 5: "unit tests per rejected flag in `src/cli.rs::tests`". `grep test.*clump_only src/cli.rs` returns zero. Only 3 rejections (`--length`, `-a`, `--rename`) are covered indirectly by the integration suite.
- **Impact:** Adding a new trim/filter flag to `Cli` in a future PR could silently miss the `--clump_only` rejection block. A per-flag test matrix would catch that.
- **Fix:** Add table-driven test in `src/cli.rs::tests` iterating over `(flag_argv, expected_substring)` pairs.

---

## 3. Notable but non-blocking

- **`main.rs:359, 2074`** duplicate the `norm` closure inline instead of using `io::norm_path` (already available and unit-tested at `src/io.rs:34`). Consistency nit.
- **`clump_only.rs:363`** — `fastqc::run` is called sequentially for R1 then R2 (PE); the trim path parallelizes elsewhere. Minor perf gap under `--fastqc`; not user-visible for record correctness.
- **`clump_only.rs:334-336, 488-493`** — the empty-input gzip-member emission is the right call; a plain zero-byte `.gz` file breaks downstream `zcat`/`gunzip`. Deviation is documented.
- **`clump_only.rs:395-405`** — startup line contains "clump-only" as required by PLAN §V.9 startup-sentinel. Good.
- **Deviation-log accuracy:** The "cores >= 2 relaxed to >= 1" deviation is defensible; `clump.rs::resolve_layout` still enforces `cores >= 1`, so no undefined behaviour. Reasonable.
- **Cross-run determinism:** Verified — sort primitives at `clump.rs:271, 289` use `sort_by` (stable) with a full content cascade (`key → seq → qual → id`). Same-input → same-output holds absolutely.

---

## 4. What the implementation does well

The byte-identity contract is enforced at three complementary layers — sort determinism at the primitive layer (already-stable + content-tiebreaker cascade at `clump.rs:271,289`), gzip-member concatenation at the writer layer (`write_records_member` per-bin fresh encoder + `.finish()`), and rejection at the CLI layer (~20 explicit `bail!` arms with per-flag messaging). The mode is a genuinely isolated code path — zero shared mutable state with the trim pipeline, zero touch on `parallel.rs::read_single_clumpy` / `read_pairs_clumpy`, and no risk of regressing Perl-parity byte-identity. Contract-scope framing (plus-line and CRLF normalization scoped out honestly) is well-handled in docs, `--help`, and the module doc-comment. The integration test's `paste - - - - | sort | diff` pattern correctly catches the class of bug where per-line sort would miss cross-record line-mixups (e.g. R1 quality on R2 sequence).
