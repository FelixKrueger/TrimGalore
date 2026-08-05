# Plan Coverage Report

**Mode:** B (code vs. plan)
**Plan:** `plans/08032026_reject-non-restartable-input/PLAN.md` (v3)
**Implementation spec:** §5 (5.1–5.4) + §9; requirements audited from §3, §3.4, §8, §12
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`, uncommitted working tree
**Date:** 2026-08-04
**Verdict:** **COMPLETE** — 0 items MISSING, 0 PARTIAL. 4 DEVIATED (all documented, all justifications verified). 1 validation row (V3) is CI-gated and cannot be run locally.

## Summary

| Status | Count |
|---|---|
| Total items | 62 |
| DONE | 55 |
| PARTIAL | 0 |
| MISSING | 0 |
| DEVIATED | 4 |
| DEFERRED (CI-gated) | 1 |
| N/A (retired / out of scope by the plan) | 2 |

Test run: `cargo test --release` from the crate root, **exit 0, 375 unit + 92 integration, 0 failed**
(`integration_non_restartable_input`: 5 passed in 0.06 s). `cargo fmt --all -- --check` clean;
`cargo clippy --all-targets --release -- -D warnings` clean. All 14 tests named below were observed
running and passing by name — none is reported as verified on inspection alone.

Files audited: `src/format.rs`, `src/cli.rs`, `src/main.rs`, `src/demux.rs`,
`tests/integration_non_restartable_input.rs`, `.config/nextest.toml`,
`.github/workflows/ci.yml`, `CHANGELOG.md`.

---

## A. §5.1 Source changes (10 items)

| # | Item | Source | Status | Notes |
|---|---|---|---|---|
| S1a | Layer 1: replace the `cli.input` `path.exists()` loop with an `fs::metadata()` match | §5.1.1 | DONE | `cli.rs:947-950` → `check_restartable_input(path, "Input file not found")`; helper at `cli.rs:482-495`. Loop is **unconditional** inside `Cli::validate`, so every mode is covered |
| S1b | Layer 1 at `--passthrough` | §5.1.1 | DONE | `cli.rs:900-902`, replaces `if !pt.exists()` |
| S1c | Layer 1 at the `--demux` barcode file | §5.1.1 | DEVIATED | Not applied; `cli.rs:942` still `if !demux_file.exists()`. Documented as **D1** — justification independently verified, see §H |
| S2 | `main.rs:768` `FastqReader::sanity_check` → `sanity_check_any` | §5.1.2 | DONE | `main.rs:763-771`. `sanity_check_any` (`main.rs:29-43`) dispatches through `detect_input_format` then still calls `FastqReader::sanity_check` on the FASTQ arms, exactly as §3.2 specified. Grep confirms **no** `FastqReader::sanity_check` call remains on any passthrough path |
| S3 | Layer-2 type check after `File::open`, **before** the read; `#[cfg(unix)]`-gated with a `#[cfg(not(unix))]` no-op | §5.1.3 | DONE | `format.rs`: `is_non_restartable_type` has both cfg arms; called in `detect_input_format` between the open and `read_filled`. Ordering verified in the diff. Minor form difference: the `use std::os::unix::fs::FileTypeExt` sits inside the cfg-gated fn rather than at module scope — equivalent, and no bare unix `use` escapes a cfg gate anywhere in the change |
| S4 | Loop the first read until 4 bytes or EOF, retrying `Interrupted` | §5.1.4 | DONE | `read_filled()`; used for **both** reads |
| S5 | Add `reopen_restarts`, private, above `detect_input_format`, per §4 | §5.1.5 / §4 | DONE | Signature `fn reopen_restarts(path: &Path, first: &[u8]) -> Result<bool>` matches §4 exactly; private; placed above `detect_input_format`; §4's doc comment reproduced; `.with_context()` on the second open carries the §4 wording verbatim |
| S6 | Call it after the `n == 0` check and before the `peek[0] == b'@'` branch | §5.1.6 | DONE | Verified in the diff; empirically confirmed by the empty-file and short-file probes below |
| S7 | Reword the `format.rs:227-230` comment ("not required for correctness" is now false) | §5.1.7 | DONE | Now: "Re-open rather than seek: sidesteps a `MultiGzDecoder` interaction with the consumed prefix, and non-restartable input is rejected above." Two lines, per CLAUDE.md |
| S8 | CHANGELOG entry under `#### Fixes` in `Unreleased` | §5.1.8 | DONE | `CHANGELOG.md:86`, inside `### Unreleased` (line 4) → `#### Fixes` (line 84). Covers the three reproductions, the `--passthrough` case, and the tty residual §11 asked to be mentioned |

**Undocumented but benign form deviation.** §4 said "the type checks stay **inline** … each is two lines".
The implementation factors them into `format::is_non_restartable_type` + `format::not_restartable_message`
(both `pub(crate)`) plus `cli::check_restartable_input`. This is not in §12's deviation list, but it
*serves* §3.3's stated requirement that the two layers' wording "cannot drift", and has no behavioural
consequence. Recorded, not counted as a gap.

## B. §5.2 CI backstop (2 items)

| # | Item | Status | Notes |
|---|---|---|---|
| C1 | `timeout-minutes: 30` on the `rust-tests` job | DONE | `ci.yml:32`, on the matrix job that runs both `cargo nextest run` (`:70`) and `cargo test --release` (`:77`) — the four job-steps §5.2 identified |
| C2 | `.config/nextest.toml` with `slow-timeout = { period = "60s", terminate-after = 4 }` | DONE | New file, byte-for-byte the value §5.2 specified, under `[profile.default]` (the profile CI actually uses) |

## C. §5.3 FIFO test construction — the six rules (6 items)

| # | Rule | Status | Notes |
|---|---|---|---|
| R1 | A writer must open the FIFO, spawned before the call | DONE | `spawn_fifo_writer` (`format.rs` test mod) opens write-only **on a spawned thread**, precisely the trap §5.3 warns about. Used by `detect_rejects_a_fifo_with_the_pipe_message`. Correctly *not* used by the layer-1 test or the integration tests, which need writer-less FIFOs |
| R2 | Bound the whole test body, setup included; thread + `mpsc::recv_timeout`; no `timeout(1)` | DONE | `bounded()` uses a thread + `recv_timeout(10s)` and panics; `run_bounded()` uses `try_wait` + a 20 s deadline. No `timeout(1)` anywhere. **Form deviation:** setup (`fresh_tmpdir`, `make_fifo`, `spawn_fifo_writer`) sits *outside* the bound. Substance holds — the wedge shape R2 was written against (Shape B, `child.wait()` on the writer) does not exist here: `mkfifo` exits immediately, the writer's blocking `open` runs on a detached thread, and the `JoinHandle` is never joined. No setup operation can block on the calling thread |
| R3 | The writer must tolerate `EPIPE` | DONE | `if let Ok(mut f) = …open(&p)` and `let _ = f.write_all(…)` — neither the open nor the write is unwrapped. No spawned `cat` has `.success()` asserted (the two `sh -c` pipe tests assert `!ok` on the *binary*, and `sh` reports the last pipeline member's status) |
| R4 | Assert the fixture is a FIFO (`mkfifo` exit status **and** `is_fifo()`) | DONE | Both `make_fifo` implementations (unit + integration) check `st.success()` *and* `fs::metadata(path).file_type().is_fifo()`. This is also what satisfies V9 |
| R5 | One `fresh_tmpdir` slug per FIFO test | DONE | Repo-wide sweep of all 99 `temp_dir().join(…)` / `fresh_tmpdir(…)` / `fresh_dir(…)` sites: **zero duplicate slugs**. The new slugs are `tg_format_fifo_detect`, `tg_format_fifo_stat`, `tg_format_reopen`, `tg_format_devfd`, `tg_format_devfd_repeat`, `tg_379_fifo_in`, `tg_379_stdin_plain`, `tg_379_stdin_gz`, `tg_379_passthrough`, `tg_379_control` |
| R6 | Integration recipe: `spawn()`, keep the `Child`, `kill()` then `wait()` on timeout — never `output()` | DONE | `run_bounded()` at `tests/integration_non_restartable_input.rs:60-91`: `spawn()`, stderr drained on a separate thread (so a chatty child cannot deadlock the wait loop), `try_wait` poll to a deadline, then `kill()` + `wait()` + panic. No `output()` in the file. The "kill the writer child" clause is moot — the integration tests use writer-less FIFOs and spawn no writer child |

**Harness proven capable of failing (independently, without touching source).** §12 records a mutation
check; because this audit must not modify source, I verified the two bounding mechanisms directly by
compiling standalone `rustc` programs that replicate `bounded()` and `run_bounded()` verbatim in shape
and point them at a writer-less FIFO:

```
bounded()      → panicked after 3.004 s  ("blocked; it must fail, not hang")
run_bounded()  → killed the child + panicked after 3.042 s
```

Both convert a genuine indefinite `open(2)` block into a failure. V5's structural claim holds.

## D. §5.4 Tests (10 items)

All ten confirmed present **and observed passing by name** in the `cargo test --release` run.

| # | Required test | Status | Test name / location |
|---|---|---|---|
| T1 | A regular plain FASTQ is accepted (happy-path regression guard) | DONE | `detect_fastq_plain_from_at_sign` (pre-existing; now also guards the new code) |
| T2 | A regular gzipped FASTQ is accepted | DONE | `detect_fastq_gz_from_plain_gzip` |
| T3 | A FIFO is rejected with the pipe message — **one** test | DONE | `detect_rejects_a_fifo_with_the_pipe_message`. Exactly one, per §5.3; asserts both `cannot be re-read from the start` **and** `is a pipe or FIFO` |
| T4 | An empty regular file still produces the *empty* error, not either new message | DONE | `detect_empty_file_errors`, strengthened per **D3** |
| T5 | `reopen_restarts(&regular_fq, b"@rea")? == true` | DONE | `reopen_restarts_answers_both_directions` |
| T6 | `reopen_restarts(&regular_fq, b"XXXX")? == false` | DONE | same test |
| T7 | `reopen_restarts("/dev/urandom", …)? == false` | DONE | `reopen_restarts_false_for_a_streaming_char_device` |
| T8 | `/dev/fd/N` over a 4-byte file rejected, `#[cfg(target_os = "macos")]` | DONE | `detect_rejects_dev_fd_over_a_four_byte_file` |
| T9 | Integration: FIFO, `/dev/stdin` from a pipe, `--passthrough <fifo>`; non-zero exit + `cannot be re-read from the start` | DONE | `fifo_input_is_rejected_and_creates_no_output`, `dev_stdin_from_a_plain_pipe_is_rejected`, `dev_stdin_from_a_gzipped_pipe_is_rejected`, `passthrough_fifo_is_rejected_and_creates_no_output`. Two of the three surfaces additionally assert `!out.exists()` |
| T10 | Do **not** assert on `/dev/stdin` redirected from a *regular* file (row B, platform-divergent) | DONE | No such assertion exists. Verified by reading the whole file |

**Beyond the plan (both verified running and passing):**
- `path_stat_names_a_fifo_without_opening_it` — isolates layer 1's decisive property (`fs::metadata`
  does not block where `File::open` would) and asserts `is_non_restartable_type` in **both** directions.
- `detect_rejects_dev_fd_when_the_leading_bytes_repeat` — the test §12 says retires A2. Confirmed it
  isolates the start-offset condition: with `@a@a@a@a`, the byte comparison alone would return
  restartable, so only the offset check can reject it.
- `regular_file_still_runs` (integration) — the positive control without which every rejection
  assertion would pass on a build that rejected everything.

## E. §3.4 Edge cases (12 items)

Each row checked for the **guard or branch**, not merely the function. Rows marked "probed" were
confirmed against the built `target/release/trim_galore`.

| # | Case | Required handling | Status | Evidence |
|---|---|---|---|---|
| E1 | Empty regular file | Existing `is empty`, unchanged | DONE | Probed: `Error: Input file '…/zero.fq' is empty`. `n == 0` precedes the second-open call in `detect_input_format`. Test T4 asserts the restartability message is **absent** |
| E2 | File shorter than 4 bytes | Length-aware compare; regular short file falls through; the same file via `/dev/fd/N` now **rejected** | DONE | Probed both directions: `@` (1 B) → `Truncated FASTQ: missing sequence line after @`; `@abc` (4 B) → same class; `X` (1 B) → `is not recognised as FASTQ…`. The **same 4-byte file via `/dev/fd/9`** → `cannot be re-read from the start: opening it a second time did not return to the beginning`. That is the row-B-shaped input that fstats as *regular*, so it proves the offset condition is load-bearing in production, not only in the macOS unit tests |
| E3 | Empty FIFO (writer opened, wrote nothing, exited) | **Behaviour change:** pipe message, previously `is empty` | DONE | Realised by construction: the type check precedes `read_filled`, verified in the diff. No test isolates this specific row — noted, but §5.4 asked for one FIFO test and the code path is identical to T3's |
| E4 | FIFO with **no writer at all** | Rejected by layer 1 without opening (v2's Medium residual) | DONE | Probed: rejected promptly, no output dir. Covered by `path_stat_names_a_fifo_without_opening_it` (unit) and `fifo_input_is_rejected_and_creates_no_output` (integration, writer-less by design) |
| E5 | Socket | Rejected by layer 2 step 2, or layer 1 | DONE (inspection) | `ft.is_socket()` present in `is_non_restartable_type`, which both layers call. Empirical check not possible here: the sandbox denies `AF_UNIX` `bind` (`PermissionError: Operation not permitted`). Flagged as a harness limitation, not a finding |
| E6 | `/dev/null` | Char device, `n == 0` → existing `is empty` | DONE | Probed: `Error: Input file '/dev/null' is empty` |
| E7 | `/dev/zero` | Offset 0, bytes equal → falls through → existing not-recognised error | DONE | Probed: `Error: Input '/dev/zero' is not recognised as FASTQ (plain or gzipped) or unaligned BAM` |
| E8 | `/dev/urandom` | Offset 0, bytes differ → the **new** re-open message | DONE | Probed: exactly the §3.3 re-open message. Also T7 |
| E9 | `/dev/stdin` at an interactive terminal | Residual, out of scope; §11 asked for a CHANGELOG sentence | DONE | Not fixed (correct — §11 declared it out of scope for option A) and named in the CHANGELOG entry, with the `isatty` reason |
| E10 | Directory | `read` fails with `EISDIR` — unchanged | DONE | Probed: `Failed to read from …` / `Is a directory (os error 21)`. Layer 1's `fs::metadata` succeeds for a directory exactly as `exists()` did, so nothing shifted |
| E11 | Second open fails with an unexpected `io::Error` | Propagated as itself **with context**, not as non-restartable (§10 Q3) | DONE (inspection) | `File::open(path).with_context(ctx)?` in `reopen_restarts`, and `reopen_restarts(…)?` at the call site propagates `Err` rather than treating it as `false`. Context string matches §4 verbatim |
| E12 | Regular file being concurrently appended | Restartable; out of scope | N/A | Explicitly out of scope in the plan |

## F. §9 Validation table (10 items)

The brief flagged V5, V6, V7 and V10 as the rows most likely to "quietly not hold". All four were
checked against something capable of failing.

| # | Row | Status | Evidence |
|---|---|---|---|
| V1 | Three #379 reproductions say the right thing, and each **stops saying its own wrong thing** | DONE | All three reproduced against the release binary: `<(cat plain.fastq)` → `Input '/dev/fd/11' is a pipe or FIFO…`; `gzip -c … \| trim_galore /dev/stdin` → `Input '/dev/stdin' is a pipe or FIFO…`; `mkfifo f; cat small.fq > f & trim_galore f` → pipe message, promptly. Per-row negative controls are in the integration tests (`doesn't seem to be in FastQ format` for the plain row, `Failed to decompress first block` for the gz row) — i.e. the vacuous-pass defect v2 had is fixed. **Note:** the automated suite substitutes a plain pipe for `<(…)`, because `<(…)` is not POSIX `sh`; I verified the process-substitution form manually and it is caught by layer 1 |
| V2 | Regular files unaffected; outputs byte-identical to `dev` | DONE | Full suite green (467 tests, including the output-content integration suites). Diff inspection confirms **no** transformation, writer or report code is touched — every added branch is a pre-reader rejection. §12 records the 6-output + 6-report byte-identity comparison against a `HEAD` worktree build across SE-gz / SE-plain / PE `--cores 1` / PE `--cores 4`, with a sentinel negative control that demonstrably fired (iteration log item 3). **I did not re-run that md5 comparison**; it is recorded, not re-derived |
| V3 | Byte-identity against Perl 0.6.11 | **DEFERRED — CI-gated** | Requires the CI `validation` / `validation-ubam` jobs; cannot run locally. §12 already records this as the one outstanding pre-merge check. Not an implementation gap |
| V4 | The empty-file error still wins | DONE | Probed (`Input file '…' is empty`, no restartability message) and asserted in both directions by T4 (D3's strengthening) |
| V5 | **No FIFO test hangs**, and nothing can hang for long | DONE | New integration file completes in **0.06 s**, whole suite in seconds. Both CI bounds in place (C1, C2). And the bounding mechanisms were **proven capable of firing** — standalone replicas of `bounded()` and `run_bounded()` pointed at a writer-less FIFO panicked/killed at 3.00 s and 3.04 s rather than hanging. This is the structural claim §5.2 exists to protect, and it holds |
| V6 | The writer-exited FIFO — the case that hangs on `dev` today | DONE | Probed both variants: live-writer FIFO (`cat small.fq > f &`) and **writer-less** FIFO — both rejected promptly with the pipe message and no output. The writer-less case is the strictly harder one (it is what used to hang forever), and is what the integration test uses. §12's "one thing the plan did not anticipate" correctly explains why the literal `cat x > f &` shape cannot be used inside a harness now: layer 1 never opens the FIFO, so the writer blocks in `open` forever |
| V7 | Every guard can fail, in **both** directions | DONE | Layer 1: writer-less FIFO rejected without opening (`path_stat_names_a_fifo_without_opening_it`, observed passing) / regular file proceeds (same test asserts `!is_non_restartable_type`, plus `regular_file_still_runs` end-to-end). Layer 2 type check: FIFO → pipe message (T3) / regular → passes (T1, T2). Second-open check: `true` for `b"@rea"`, `false` for `b"XXXX"`, `false` for `/dev/urandom`, `false` for macOS `/dev/fd/N` over 4 bytes, `false` for the repeating-prefix case. Additionally the `/dev/fd/9` probe shows the offset condition rejecting real production input, so this predicate is no longer "only ever observed passing" |
| V8 | uBAM input still detected | DONE | `detect_unaligned_bam_via_decompressed_magic`, `detect_paired_ubam_via_decompressed_magic`, `detect_bgzipped_fastq_is_fastq_not_bam` all pass, plus `integration_ubam` (8) and `integration_ubam_out` (23) and `integration_clump_only_ubam` (15). Every added layer is format-agnostic |
| V9 | The Linux mechanism is measured, not assumed | DONE (code); Linux leg CI-gated | The required assertion exists in **both** `make_fifo` helpers: `fs::metadata(path).file_type().is_fifo()` on the fixture immediately before it is handed to the binary, plus the `mkfifo` exit-status check. Green on Darwin here. The `ubuntu-latest` leg is the first CI run, as the plan intended — that is by design, not a gap. Row B correctly left unasserted (T10) |
| V10 | **Every** input-bearing CLI surface reaches a guard | DEVIATED | `cli.input` and `cli.passthrough` are each asserted rejected for a FIFO by their own integration test. The third surface, the `--demux` barcode file, is argued out rather than covered — **D1**, justification verified (§H). Two observations: (a) the enumeration is *by test*, not *by mechanism*, so V10's forward-looking intent — "the next flag that takes a path should be caught by a test rather than by a reviewer" — is not structurally met; nothing fails if a new path-bearing flag is added unguarded. (b) A writer-less FIFO passed to `--demux` still blocks in `File::open` with no message. That is not a restartability defect and not a regression (it predates this change), but it is the one remaining "hang with no message" on an input-bearing flag |

## G. §8 Assumptions (9 items)

| # | Assumption | Status | Evidence |
|---|---|---|---|
| A1 | Every currently-working input is restartable; nothing in `tests/`, `.github/`, `justfile` uses `/dev/stdin`, `mkfifo` or process substitution | DONE | Re-verified by grep over `tests/ .github/ justfile`: **zero** hits outside the new `integration_non_restartable_input.rs`. Positive control run (the pattern does match the new file) so the negative is trustworthy |
| A2 | *Retired in v3* | N/A | Retired by the offset condition. §12's extra test `detect_rejects_dev_fd_when_the_leading_bytes_repeat` is what makes the retirement testable, and it passes |
| A3 | `File::open` on a writer-less FIFO blocks | DONE | Re-confirmed a fourth time here: the standalone `bounded()` replica blocked until the 3 s bound fired |
| A4 | `detect_input_format` is on no hot path | DONE | Unchanged by this work; call-site inventory intact |
| A5 | For a regular file the `n == 0` check still precedes the second-open check | DONE | Verified in code order and probed (E1). For a FIFO the type check wins instead, as §3.4 intends (E3) |
| A6 | `fstat` names a FIFO as a FIFO on both platforms | DONE (Darwin) | Verified on Darwin by test and probe. Linux leg is V9's CI run, which is how the plan chose to discharge it |
| A7 | Unix-only — "handled rather than flagged" | **DONE — claim verified** | `is_non_restartable_type` has a `#[cfg(unix)]` arm and a `#[cfg(not(unix))]` no-op returning `false`. Every unix-only `use` is inside a cfg-gated item: `FileTypeExt` inside the gated fn and inside `#[cfg(unix)] fn make_fifo`; `AsRawFd` inside the two `#[cfg(target_os = "macos")]` tests. The integration file carries a file-level `#![cfg(unix)]`. `cli::check_restartable_input` is portable. No bare unix `use` at module scope, so the `src/lib.rs` compile-break §5.1.3 warned about does not exist |
| A8 | Both reads must return everything asked for — "handled by the looped reads, not assumed" | **DONE — claim verified** | `read_filled()` loops to buffer-full-or-EOF and retries `ErrorKind::Interrupted`; it is used for **both** the first read in `detect_input_format` and the second read in `reopen_restarts`. The comparison is length-aware (`n == first.len() && &again[..n] == first`), closing the ≤4-byte prefix hole — probed live via `/dev/fd/9` over a 4-byte file. The claim that the same loop hardens the pre-existing `n >= 3` gzip test also holds: that branch now sees a filled `peek` |
| A9 | No input-bearing CLI surface bypasses the guard — "asserted by V10 rather than by inspection" | **DONE for the two live surfaces; see V10/D1** | Layer 1's loop is unconditional inside `Cli::validate`, and `cli.validate()` has exactly **one** call site (`main.rs:166`) ahead of all dispatch — so every mode (paired, SE, specialty, clump-only, uBAM-out) inherits it for `cli.input`. `cli.passthrough` covered at `cli.rs:900`. The `--demux` barcode file is the one surface *not* covered, deliberately (D1) |

## H. §12 Deviations — are they deviations, or gaps relabelled? (3 items)

| # | Deviation | Verdict | Why |
|---|---|---|---|
| D1 | The `--demux` barcode file is **not** checked; V10's enumeration drops from three surfaces to two | **DEVIATED — justification holds** | Verified independently against `src/demux.rs`. `read_barcode_file` (`demux.rs:31-34`) performs exactly **one** `File::open` and streams the file with `BufReader::lines()` in a single pass. Grep confirms exactly one production call site (`main.rs:1265`), and `demux::demultiplex` receives already-parsed `BarcodeEntry` values — it never re-opens the path. So the invariant this guard enforces ("re-opening yields the same bytes from the start") genuinely does not apply, and D1's conclusion that rejecting a FIFO barcode file would be an unjustified regression is correct. This narrows §3.0/§5.1's stated scope, and the narrowing is sound. **One residual worth recording** (not a gap against this plan): a *writer-less* FIFO barcode file still blocks in `File::open` with no message — a hang, but not a restartability failure, and identical to `dev` today |
| D2 | A `stream_position()` error falls through to the byte comparison instead of propagating | **DEVIATED — justification holds** | `if file.stream_position().is_ok_and(\|pos\| pos != 0) { return Ok(false); }`. Strictly less aggressive than §3.1 step 5, as claimed: only a *successfully read* non-zero offset is decisive. This is what preserves §3.4's `/dev/null` and `/dev/zero` rows, both of which I probed and both of which land where §3.4 predicted. It does not weaken §10 Q3, which concerns the second `open` failing — that still propagates with context (E11) |
| D3 | The pre-existing `detect_empty_file_errors` test was **strengthened**, not merely left alone | **DONE — a genuine improvement, correctly recorded** | It previously asserted only `is_err()`, which would have passed if an empty file had started being reported as non-restartable. It now asserts `is empty` is present **and** `cannot be re-read` is absent — which is exactly what V4 demands. Observed passing |

---

## Gaps (detail)

None. No item is MISSING and none is PARTIAL.

The four DEVIATED rows are two underlying decisions, each recorded in §12 with a justification that
this audit verified against the code:

- **S1c / D1 / V10** — one decision, seen from three angles: the `--demux` barcode file is excluded
  because it is read once, not re-read. Verified in `src/demux.rs`.
- **D2** — a `stream_position()` error defers to the byte comparison rather than becoming a new hard
  error for exotic character devices. Verified by probe against `/dev/null`, `/dev/zero`, `/dev/urandom`.

## Outstanding before merge (not implementation gaps)

1. **V3 — Perl 0.6.11 byte-identity.** Needs the CI `validation` and `validation-ubam` jobs; not
   runnable locally. §12 flags it as the invariant that must not move. Evidence that it will hold is
   strong (no writer or transformation path is touched, and V2's local byte-identity comparison
   passed with a firing sentinel control) — but it is unverified until CI runs.
2. **V9 / A6 — the Linux leg.** The `ubuntu-latest` matrix entry is the measurement the plan chose
   for A6. Deliberate, and the fixture assertion that makes a platform surprise name itself is in
   place.
3. **#374 merge order.** §7.1 recommends merging Benjamin Demaille's PR first and re-checking the
   `main.rs` line references. This work was implemented against `dev` @ `40f77ec`; the one changed
   `main.rs` line (`:770`) is not among the lines #374 edits, and layer 1 lives in `cli.rs`, which
   #374 does not touch.

## Verification commands run for this audit

```
cargo test --release                                  # exit 0 — 375 unit + 92 integration, 0 failed
cargo fmt --all -- --check                            # clean
cargo clippy --all-targets --release -- -D warnings   # clean
```

Plus, against `target/release/trim_galore`: the three #379 reproductions (process substitution,
gzipped pipe → `/dev/stdin`, live-writer FIFO), a writer-less FIFO, `/dev/null`, `/dev/zero`,
`/dev/urandom`, an empty regular file, a directory, 1- and 4-byte regular files, and the same 4-byte
file via `/dev/fd/9`. Plus two standalone `rustc` programs replicating `bounded()` and `run_bounded()`
to prove each converts a real indefinite block into a failure.

Not verifiable in this environment: a Unix-socket input (sandbox denies `AF_UNIX bind`), and
`/dev/stdin` at an interactive terminal (no tty).

## Verdict

**COMPLETE.**

Every item in §5.1–§5.4, §5.2's CI backstop, §5.3's six FIFO rules, §3.4's twelve edge cases and
§9's validation table is implemented as specified or deviates for a reason recorded in §12 and
confirmed here against the code. Nothing is missing; nothing is partial.

The four rows the brief singled out as most likely to "quietly not hold" all hold, and each was
checked against something capable of failing rather than merely observed passing:

- **V5** — the two bounding harnesses were proven to fire on a real `open(2)` block (3.00 s, 3.04 s).
- **V6** — both FIFO variants rejected promptly, including the writer-less case that used to hang forever.
- **V7** — all six directions exercised, and the offset predicate additionally caught real production
  input (`/dev/fd/9` over a 4-byte file), so it is no longer a check only ever seen to pass.
- **V10** — two of three surfaces asserted; the third excluded on a justification verified against
  `src/demux.rs`. The residual is that the enumeration is by test, not by mechanism, so a *future*
  path-bearing flag added unguarded would not fail anything.

One verification remains genuinely open and cannot be closed locally: **V3** (Perl 0.6.11
byte-identity), which needs the CI `validation` / `validation-ubam` jobs before merge.
