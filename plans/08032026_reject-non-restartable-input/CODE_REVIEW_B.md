# CODE_REVIEW_B — Reject non-restartable input (#379)

**Reviewer:** B (independent, fresh context)
**Target:** uncommitted working-tree diff on `dev` @ `40f77ec` — `src/format.rs`, `src/cli.rs`, `src/main.rs`, `.github/workflows/ci.yml`, `CHANGELOG.md`, `tests/integration_non_restartable_input.rs` (new), `.config/nextest.toml` (new)
**Spec:** `plans/08032026_reject-non-restartable-input/PLAN.md` v3
**Constraint honoured:** no file was modified. Every fix below is given as an exact edit for the author to apply.

---

## 1. Summary

The three-layer guard is **behaviourally correct**. I probed all twelve §3.4 edge-case rows that are reachable on Darwin and every one lands where the plan says it should — including the four that the plan's own reviewers argued about (`/dev/null` → `is empty`, `/dev/zero` → not-recognised, `/dev/urandom` → Reopen message, `/dev/stdin` from a regular-file redirect → Reopen message). I found no input that is wrongly accepted, and no input that terminates in neither direction other than the tty residual the plan already books as a known case.

I also confirmed the guard's testing is not self-satisfied: **both** layers have a test that fails if the layer is deleted, in the fast direction, not the wedge direction.

One finding is material and I would not merge without it:

> **H1 — the layer-0 swap at `src/main.rs:770` makes `--passthrough <uBAM>` strictly worse.** `sanity_check_any` *accepts* a BAM passthrough file (it routes to `BamReader`), where `FastqReader::sanity_check` rejected it. The run then proceeds, creates the R1/R2 and passthrough output files, and dies with `stream did not contain valid UTF-8`. Reproduced. The plan's justification for this being safe (§3.2: "`--passthrough` + uBAM is rejected at `main.rs:253`") is false, and §2.4 of the same plan says why.

Everything else is Medium or below: a test helper that reports assertion failures as hangs, an errno thrown away at layer 1, and a latent precondition in `reopen_restarts`.

### Claims I re-ran rather than accepted

| Claim | Verdict |
|---|---|
| `cargo fmt --check` clean | Confirmed. |
| `cargo clippy --all-targets --release -D warnings` clean | Confirmed (cache replay, no warnings). |
| Full suite, 0 failures | Confirmed. `cargo nextest run`: **467 tests, 467 passed, 8.727 s**. `cargo test --release`: exit 0. |
| "a mutation makes the FIFO test fail in 10 s rather than hang" | Sound, and stronger than stated — it fails on **either** branch of the race. If the second open blocks, `bounded` times out; if the writer is still alive so the second open succeeds, `reopen_restarts` reads mid-stream, returns `false`, and the `is a pipe or FIFO` assertion fails outright. |
| Byte-identity vs `HEAD` for 6 outputs + 6 reports | Not re-run; out of practical scope for this review. No writer path is touched by the diff, which is consistent. V3 (Perl 0.6.11) is correctly still booked as CI-only. |
| §5.3 rule 5 (one `fresh_tmpdir` slug per FIFO test) | Confirmed. All five new unit slugs and all five new integration slugs are unique, and no other test file uses a `tg_379_*` or `tg_format_fifo*` / `tg_format_devfd*` / `tg_format_reopen` slug. |

### What I verified end to end

Run against `./target/release/trim_galore` built from the working tree, each bounded with `perl -e 'alarm N; exec @ARGV'`:

| Input | Result | Plan row |
|---|---|---|
| writer-less FIFO | `is a pipe or FIFO`, exit 1, prompt, no output dir | §3.4 "no writer at all" ✓ |
| `cat f \| … /dev/stdin` | `is a pipe or FIFO`, exit 1 | §1 row 2 ✓ |
| `<(cat f)` → `/dev/fd/63` | `is a pipe or FIFO`, exit 1 | §1 row 1 ✓ |
| `/dev/stdin < f` (macOS row B) | `…opening it a second time did not return to the beginning of the data` | §2.2.1 row B ✓ |
| `/dev/null` | `Input file '/dev/null' is empty` | §3.4 ✓ |
| `/dev/zero` | `is not recognised as FASTQ …` | §3.4 ✓ |
| `/dev/urandom` | Reopen message | §3.4 ✓ |
| directory | `Failed to read from <dir>` / `Is a directory (os error 21)` | §3.4 ✓ (unchanged) |
| `--passthrough <pipe>` | `is a pipe or FIFO`, **no output dir created** | §2.4 fixed ✓ |
| `--passthrough /dev/stdin < f` | Reopen message | the case layer 0 uniquely buys (see L3) |

The rendered message is correct, including the `\x20   ` indent trick surviving the `\`-continuations:

```
Error: Input '/tmp/.../stream.fq' is a pipe or FIFO, not a regular file, so it cannot be re-read from the start.

Trim Galore reads each input more than once — format detection, the initial
sanity check, and adapter auto-detection each open it independently — so it
requires a regular file.

Common causes are pipes, FIFOs, process substitution such as
`<(zcat reads.fq.gz)`, and /dev/stdin. Write the stream to a file first:

    zcat reads.fq.gz > reads.fq && trim_galore [options] reads.fq
```

---

## 2. Issues by area

### 2.1 Logic

#### H1 — `--passthrough <uBAM>` now passes the sanity check and fails deep, with partial output on disk

`src/main.rs:763-771`:

```rust
// `sanity_check_any`, not `FastqReader::sanity_check`, so the
// passthrough path reaches the #379 restartability guard.
if pair_idx == 0
    && let Some(ref pt_path) = cli.passthrough
{
    sanity_check_any(pt_path)?;
}
```

`sanity_check_any` (`src/main.rs:29-42`) dispatches on `detect_input_format`, and its `UnalignedBam` arm opens a `BamReader` and returns `Ok(())`. The passthrough file is then read, unconditionally and on both paths, by `FastqReader` — `src/main.rs:1346` (`open_threaded`) and `src/main.rs:1385` (`open`).

The plan asserts this arm is unreachable:

> §3.2: "`--passthrough` + uBAM is rejected at `main.rs:253`, so `sanity_check_any`'s BAM arm is unreachable here and harmless."

`src/main.rs:253` is `if cli.passthrough.is_some() && any_bam`, and `any_bam` is built at `src/main.rs:184-191` by mapping over **`cli.input`**. The passthrough file is not in `cli.input` — which is the whole premise of §2.4, and which §2.4 states in terms:

> "the `--passthrough` + uBAM rejection at `main.rs:253` keys off `any_bam`, computed from `cli.input`, so the passthrough file's container format is never inspected anywhere either."

So §3.2 contradicts §2.4, and the implementation followed §3.2. Reproduced on the working-tree binary:

```
$ trim_galore --paired -o $D/out --passthrough test_files/ubam_test.bam $D/r1.fastq $D/r2.fastq
  Passthrough: test_files/ubam_test.bam → …/out/ubam_test_passthrough.fq
Error: processing pair 1 of 1 (R1=…/r1.fastq, R2=…/r2.fastq)

Caused by:
    stream did not contain valid UTF-8
exit 1

$ ls $D/out
r1_val_1.fq  r2_val_2.fq  ubam_test_passthrough.fq      # three empty files
```

Before the change, `FastqReader::sanity_check` read the BGZF bytes as text at that same line and errored **before** `run_paired`, so nothing was created. The diff therefore moves one `--passthrough` case from "errors before any output" to "errors after three output files exist, with an opaque message" — the exact failure shape the CHANGELOG entry claims this change removes for `--passthrough`.

This is not a reason to revert the layer-0 swap: it is the only thing that catches the *Reopen* class on the passthrough path (`--passthrough /dev/stdin` — verified rejected, no output dir). The fix is to give the passthrough file the BAM rejection that `cli.input` gets.

**Exact fix** — `src/main.rs:767-771`, replace:

```rust
            if pair_idx == 0
                && let Some(ref pt_path) = cli.passthrough
            {
                sanity_check_any(pt_path)?;
            }
```

with:

```rust
            if pair_idx == 0
                && let Some(ref pt_path) = cli.passthrough
            {
                // main.rs:253 rejects passthrough + uBAM from `any_bam`, which is
                // built from cli.input only — the passthrough file needs its own check.
                if matches!(detect_input_format(pt_path)?, InputFormat::UnalignedBam) {
                    anyhow::bail!(
                        "--passthrough is not supported with uBAM input in this release. \
                         Either convert '{}' to FASTQ via `samtools fastq` first, or drop \
                         --passthrough.",
                        pt_path.display()
                    );
                }
                sanity_check_any(pt_path)?;
            }
```

`detect_input_format` is idempotent and cheap here (§6/A4), and running it twice keeps the restartability guard as the first thing the passthrough file meets. Both symbols are already imported at `src/main.rs:16-18`.

Add a test alongside the others in `tests/integration_non_restartable_input.rs`, or better in `tests/integration_passthrough.rs` next to its siblings:

```rust
#[test]
fn passthrough_ubam_is_rejected_before_output() {
    let dir = fresh_dir("tg_379_passthrough_ubam");
    let (r1, r2) = (dir.join("r1.fastq"), dir.join("r2.fastq"));
    write_fastq(&r1);
    write_fastq(&r2);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("--paired").arg("-o").arg(&out)
        .arg("--passthrough").arg("test_files/ubam_test.bam")
        .arg(&r1).arg(&r2);
    let (ok, stderr) = run_bounded(cmd, "trim_galore --passthrough on a uBAM");

    assert!(!ok, "a uBAM passthrough must be rejected; stderr was: {stderr}");
    assert!(stderr.contains("not supported with uBAM input"), "unexpected message: {stderr}");
    assert!(
        !stderr.contains("valid UTF-8"),
        "must not fail deep in the reader: {stderr}"
    );
    assert!(!out.join("r1_val_1.fq").exists(), "must reject before writing output");
}
```

Priority: **High**. Exotic invocation, but it is a diagnostic regression inside the change whose sole purpose is diagnostics, and it leaves partial output where none was left before.

#### M3 — `reopen_restarts`'s `PEEK_LEN` clamp answers silently and wrongly outside its unstated precondition

`src/format.rs:303-306`:

```rust
    let mut again = [0u8; PEEK_LEN];
    let want = first.len().min(again.len());
    let n = read_filled(&mut file, &mut again[..want]).with_context(ctx)?;
    Ok(n == first.len() && &again[..n] == first)
```

Two silent answers outside `1 <= first.len() <= PEEK_LEN`:

- `first.len() > PEEK_LEN` → `want == 4`, so `n <= 4 < first.len()` and the function returns `Ok(false)` **always**, i.e. "not restartable" for a perfectly good regular file. It cannot return `true`.
- `first.len() == 0` → `want == 0`, `read_filled` returns `Ok(0)` without reading, and `0 == 0 && [] == []` returns `Ok(true)` vacuously.

Neither is reachable today: the sole caller is `src/format.rs:334`, which passes `&peek[..n]` with `1 <= n <= PEEK_LEN` (the `n == 0` bail at `:330` guarantees the lower bound). But the doc comment at `:278-287` states only the FIFO precondition, so the next caller has nothing to read. This is the "silently produce wrong results" class that `CLAUDE.md` asks code to fail loudly on.

The `.min()` is what hides it: it converts an out-of-contract argument into a plausible-looking answer instead of a compile or runtime error.

**Exact fix** — `src/format.rs:288`, extend the doc block and assert. Add after the existing `Callers must reject FIFOs first:` paragraph:

```rust
/// `first` must be the 1..=`PEEK_LEN` bytes the first handle actually returned.
```

and replace `:303-304`:

```rust
    let mut again = [0u8; PEEK_LEN];
    let want = first.len().min(again.len());
```

with:

```rust
    debug_assert!(
        (1..=PEEK_LEN).contains(&first.len()),
        "reopen_restarts: `first` must be 1..={PEEK_LEN} bytes, got {}",
        first.len()
    );
    let mut again = [0u8; PEEK_LEN];
    let want = first.len().min(PEEK_LEN);
```

Note the tension the author will hit: this file already argues at `:119-128` that a `debug_assert` "would compile out of every shipped binary". That argument was about a `pub` function reachable from a release build; this is a private helper with one caller, so a `debug_assert` plus the doc line is proportionate. If you disagree, a plain `assert!` costs one comparison per `detect_input_format` call, which is nothing against the extra `open` it already does.

Priority: **Medium** (latent, not live).

#### Not a hole — the `is_ok_and` fall-through (D2)

`src/format.rs:299`:

```rust
    if file.stream_position().is_ok_and(|pos| pos != 0) {
        return Ok(false);
    }
```

I probed for a gap here and did not find one. The three ways `stream_position` can fail to be decisive all land safely:

- `Err(ESPIPE)` on a pipe → falls through to the byte comparison, which would read from the stream. Unreachable: `is_non_restartable_type` at `:322` has already bailed, on both the path (`cli.rs:490`) and the handle. The doc comment at `:286-287` states the precondition and both callers honour it.
- `Err(_)` on an exotic character device → falls through, byte comparison decides. `/dev/zero` (accept) and `/dev/urandom` (reject) both verified end to end, and both go through the byte half regardless because their offset is `Ok(0)`.
- `Ok(0)` on a `dup`-shared handle whose offset happens to be 0 → the byte comparison is the backstop, and `detect_rejects_dev_fd_when_the_leading_bytes_repeat` is the test that keeps the offset half honest in the other direction.

D2 is strictly more conservative than the plan (`Err` → defer rather than propagate), and I agree with the deviation for the reason given: propagating would turn a class that falls through today into a hard error.

`read_filled` (`:265-276`) is also correct. `Ok(0)` breaks on EOF, `Interrupted` retries, `n` is monotone so it terminates, and `Ok(n)` with `n < buf.len()` is a genuine short read the caller can distinguish. The only unbounded case is an unending `Interrupted` storm, which is what `std::io::Read::read_exact` does too.

The `n == 0` empty check at `:330` correctly precedes the reopen check at `:334`, so an empty regular file keeps its own message — asserted in both directions by the strengthened `detect_empty_file_errors` (D3), which is the right strengthening.

### 2.2 Errors and diagnostics

#### M2 — layer 1 throws away the errno it now has in hand

`src/cli.rs:487-489`:

```rust
fn check_restartable_input(path: &std::path::Path, not_found: &str) -> anyhow::Result<()> {
    let meta =
        std::fs::metadata(path).map_err(|_| anyhow::anyhow!("{not_found}: {}", path.display()))?;
```

`map_err(|_| …)` collapses every `io::Error` into "not found" and drops the source, so `EACCES` on a parent directory, `ELOOP` on a symlink cycle and `EMFILE` under fd pressure all print `Input file not found: reads.fq`.

This is byte-for-byte the behaviour of the `Path::exists()` it replaces (`exists()` is `fs::metadata(path).is_ok()`), so it is not a regression — but it is the one site in this diff that acquired the errno and discarded it, and the plan rejects exactly this pattern one layer down:

> §10 Q3: "`Ok(false)` would report `EMFILE` or `EACCES` as 'your input is a pipe' — the same wrong-blame this plan exists to remove, one layer down."

The same reasoning applies to reporting `EACCES` as "not found" one layer *up*. `src/cli.rs:1764` asserts on `--passthrough file not found`, so the `NotFound` wording must stay.

**Exact fix** — replace `src/cli.rs:488-489` with:

```rust
    let meta = std::fs::metadata(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            anyhow::anyhow!("{not_found}: {}", path.display())
        } else {
            anyhow::Error::new(e).context(format!("Failed to stat input '{}'", path.display()))
        }
    })?;
```

This keeps `test_passthrough_rejects_missing_file` (`src/cli.rs:1753-1765`) green and matches the `Failed to stat input file` wording already used at `src/format.rs:319-321`.

Priority: **Medium**.

#### L2 — the pipe message never names a socket

`src/format.rs:210-212` documents `NotRestartable::Pipe` as "A pipe, FIFO **or socket**, named by its file type", and `is_non_restartable_type` (`:249-253`) matches `is_fifo() || is_socket()`. The message at `:222-226` says only:

> `Input 'X' is a pipe or FIFO, not a regular file, so it cannot be re-read from the start.`

The plan spends a paragraph (§3.3) on v2's noun being wrong for the majority case; a socket gets the same treatment here, in the same change. It is still an improvement — before this, a socket path reached `File::open` and returned a bare `ENXIO` — so this is cosmetic.

**Exact fix** (optional): add a `Socket` variant to `NotRestartable`, have `is_non_restartable_type` return which one, and give it `"Input '{}' is a socket, not a regular file, so it cannot be re-read from the start."`. Cheaper alternative, which keeps the `is a pipe or FIFO` substring the two tests assert on: leave it, and drop "or socket" from the `:210` doc so the code and the message agree.

Priority: **Low**.

#### L1 — the Reopen class leaves an empty output directory

`naming::ensure_output_dir(output_dir)` runs at `src/main.rs:315`, ahead of every layer-2 site (`sanity_check_any` at `:760`/`:762`/`:770`). So the pipe/FIFO class, caught by layer 1 inside `Cli::validate()` at `src/main.rs:166`, creates nothing — verified, and asserted by `fifo_input_is_rejected_and_creates_no_output` and `passthrough_fifo_is_rejected_and_creates_no_output`. The Reopen class does not:

```
$ trim_galore --paired -o $D/o1 --passthrough /dev/stdin r1 r2 < idx.fastq
Error: Input '/dev/stdin' cannot be re-read from the start: …
$ ls -la $D/o1
total 0            # created, empty
```

The CHANGELOG says of `--passthrough`, "it is now rejected before anything is written". True for the class the entry is about (a FIFO), not for `/dev/fd`-shaped input. An empty directory is harmless; the claim is the thing to soften.

**Exact fix**: no code change. In `CHANGELOG.md`, change "it is now rejected before anything is written" to "it is now rejected before any output file is written".

Priority: **Low**.

### 2.3 Tests — can any of them hang CI, and do the bounds bound what they claim?

Short answer: **no new test can hang CI**, and the bounds are real. Longer answer, per mechanism:

**`bounded()` (`src/format.rs:505-514)`.** The body runs on a detached thread, `recv_timeout(10 s)` bounds the wait, and a panic on the harness thread is a test failure. The leaked body thread cannot hold the process open — the Rust test harness returns from `main`, which exits the process regardless of live threads. Under nextest each test is its own process, so it is moot there too. Correct.

**`spawn_fifo_writer()` (`src/format.rs:537-545)`.** Correct, and correct for the non-obvious reason: the write-only `open` is off the main thread, so it rendezvous with the reader's `open` instead of self-deadlocking (§5.3 rule 1 + the trap the plan calls out). `if let Ok(mut f) = …` and `let _ = f.write_all(…)` give rule 3 (EPIPE tolerance) without an `unwrap`. The handle is bound to `_writer`, not `_`, so it is not dropped at the end of the statement — that distinction matters here and is easy to break in a later edit, though nothing depends on the handle actually being joined.

**`make_fifo()` (`src/format.rs:517-530` and `tests/…:46-57`).** Implements rule 4 exactly: `mkfifo` exit status *and* a positive `is_fifo()` assertion on the created node. This is the check that stops a `$PATH` failure from being reported as a guard bug.

**`run_bounded()` (`tests/integration_non_restartable_input.rs:60-91`).** Implements rule 6 correctly and for the right reason — `spawn()` not `output()`, stderr drained on its own thread so a chatty child cannot deadlock the `try_wait` loop, `kill()` then `wait()` on the 20 s deadline. Two rough edges, both reachable only in the regression case:

- **L5a**: the `panic!` at `:84` fires *before* `rx.recv_timeout` at `:89`, so the failure message carries no stderr at all — on the one failure that matters most, you get `trim_galore on a FIFO never exited; it must fail, not hang` and nothing about what the binary said. Fix: drain first.
- **L5b**: for the two `sh -c "… | trim_galore …"` tests, `child.kill()` signals only `sh`. A wedged `trim_galore` and its `cat` survive as orphans still holding the inherited stderr pipe, so nextest's `leak-timeout` (default 100 ms) would add a LEAK on top of the failure. Still red, just noisier.

**Exact fix for L5a** — `tests/integration_non_restartable_input.rs:80-85`, replace:

```rust
            None if Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                panic!("{what} never exited; it must fail, not hang");
            }
```

with:

```rust
            None if Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                let stderr = rx.recv_timeout(Duration::from_secs(5)).unwrap_or_default();
                panic!("{what} never exited; it must fail, not hang. stderr was: {stderr}");
            }
```

Priority: **Low**.

#### M1 — `bounded()` reports an assertion failure as a hang

`src/format.rs:510-513`:

```rust
        match rx.recv_timeout(std::time::Duration::from_secs(10)) {
            Ok(v) => v,
            Err(_) => panic!("{what} blocked; it must fail, not hang"),
        }
```

`recv_timeout` has two error variants and `Err(_)` swallows the distinction. If the body panics — `expect_err("a FIFO must be rejected")` firing because a regression made `detect_input_format` *succeed* on a FIFO, or `expect("stat must not block")` firing because the fixture vanished — the sender drops without sending, `recv_timeout` returns `Disconnected` immediately, and the test reports `detect_input_format on a FIFO blocked; it must fail, not hang`. It did not block. It failed in 200 µs, for a different reason.

The body's own panic message does reach the captured output via the default hook, so the information is recoverable — but the headline is wrong, and the wrong headline is "the design decision this whole plan is about has regressed into a hang", which is the most alarming possible misdiagnosis. In a change whose thesis is *do not describe the failure inaccurately*, the harness describes the failure inaccurately.

**Exact fix** — `src/format.rs:506-514`, replace the body of `bounded` with:

```rust
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(body());
    });
    use std::sync::mpsc::RecvTimeoutError;
    match rx.recv_timeout(std::time::Duration::from_secs(10)) {
        Ok(v) => v,
        Err(RecvTimeoutError::Timeout) => panic!("{what} blocked; it must fail, not hang"),
        Err(RecvTimeoutError::Disconnected) => {
            panic!("{what} panicked; see the panic message above")
        }
    }
```

Priority: **Medium**.

#### L3 — the layer-0 change is the only edit with no test that fails if it is reverted

Delete-one-thing check, run by inspection over all four edits:

| Reverted | What fails |
|---|---|
| layer 1 (`cli.rs:949`/`:902`) | `fifo_input_is_rejected_and_creates_no_output` and `passthrough_fifo_is_rejected_and_creates_no_output` — the binary blocks in `File::open` on the writer-less FIFO, `run_bounded` kills at 20 s and panics. Fails in 20 s, not 6 hours. ✓ |
| layer 2 type check (`format.rs:322`) | `detect_rejects_a_fifo_with_the_pipe_message`, in ≤10 s, on either branch of the writer race. ✓ |
| layer 2 offset condition (`format.rs:299`) | `detect_rejects_dev_fd_when_the_leading_bytes_repeat` (macOS leg). ✓ |
| layer 2 byte condition (`format.rs:306`) | `reopen_restarts_false_for_a_streaming_char_device`. ✓ |
| **layer 0 (`main.rs:770`)** | **Nothing.** `passthrough_fifo_is_rejected_and_creates_no_output` uses a writer-less FIFO, which layer 1 already rejects inside `Cli::validate` — long before line 770 runs. |

So the edit that H1 shows to be actively harmful is also the only one with no test defending it. What layer 0 uniquely buys is the *Reopen* class on the passthrough path, which I verified by hand (`--passthrough /dev/stdin < f` → Reopen message, no output dir) and which nothing asserts.

**Exact fix**: fold into H1's new test, and add a macOS-only companion so the Reopen class on the passthrough surface is covered:

```rust
/// The class layer 1 cannot see: a `/dev/fd`-shaped passthrough is a regular
/// file to `stat`, so only `sanity_check_any` at main.rs:770 rejects it.
#[cfg(target_os = "macos")]
#[test]
fn passthrough_dev_stdin_is_rejected() {
    let dir = fresh_dir("tg_379_passthrough_devfd");
    let (r1, r2, idx) = (dir.join("r1.fastq"), dir.join("r2.fastq"), dir.join("idx.fastq"));
    write_fastq(&r1);
    write_fastq(&r2);
    write_fastq(&idx);
    let out = dir.join("out");

    let mut cmd = Command::new("sh");
    cmd.arg("-c").arg(format!(
        "'{}' --paired -o '{}' --passthrough /dev/stdin '{}' '{}' < '{}'",
        binary().display(), out.display(), r1.display(), r2.display(), idx.display()
    ));
    let (ok, stderr) = run_bounded(cmd, "trim_galore --passthrough /dev/stdin");
    assert!(!ok, "stderr was: {stderr}");
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
}
```

Priority: **Low** as a coverage gap in its own right; the substance is H1.

#### L4 — the per-test bound covers only one of the two CI test legs

`.config/nextest.toml` sets `[profile.default] slow-timeout = { period = "60s", terminate-after = 4 }`, which applies to `cargo nextest run` (`ci.yml:71`). The **release** leg, `cargo test --release` (`ci.yml:74`), has no per-test bound — libtest has none to give. A wedge there is caught only by `timeout-minutes: 30`, and reports as a job timeout rather than a named test failure.

The plan's §5.2 says the job bound is "the one backstop that survives an implementer forgetting every rule in §5.3", and that is exactly right — I raise this only so the asymmetry is on the record, not as a defect. `timeout-minutes: 30` is the load-bearing bound; `terminate-after` is a nicety on one leg.

On the blast radius of a repo-wide `terminate-after = 4` (240 s kill): measured margin is ~52x. Slowest test on this machine is `integration_adapter2::a2_wins_on_every_resolution_path` at **4.579 s**; whole suite 8.727 s for 467 tests. Even a 10x-slower cold runner has room. No existing test is at risk.

Priority: **Low** (informational).

#### L6 — `/dev/urandom` negative control is probabilistic

`reopen_restarts_false_for_a_streaming_char_device` (`src/format.rs:599-608`) asserts that four random bytes differ from four other random bytes. False-pass probability 2⁻³² ≈ 2.3 × 10⁻¹⁰ per run. Recording it rather than recommending a change — there is no better portable non-restartable source that cannot block, which the plan established (V7), and reading more bytes only changes the exponent.

Priority: **Low** (informational, no action).

### 2.4 Does the CLI-layer change preserve existing behaviour?

Yes. Checked all three `exists()` sites named in the prompt:

| Site | Message | Status |
|---|---|---|
| `cli.rs:902` (`--passthrough`) | `--passthrough file not found: {}` | Preserved verbatim. Asserted by `src/cli.rs:1764` (`test_passthrough_rejects_missing_file`) — still green. |
| `cli.rs:949` (`cli.input`) | `Input file not found: {}` | Preserved verbatim. No test asserts it (the identical string at `src/fastq.rs:281`, `src/bam.rs:168` and `:352` is unrelated and untouched). |
| `cli.rs:942` (`--demux` barcode file) | `Barcode file not found: {}` | **Unchanged** — deliberately, per D1. |

**D1 verified independently.** `demux::read_barcode_file` (`src/demux.rs:31-34`) opens once and streams with a `BufReader`, and `grep` shows exactly one call site (`src/main.rs:1265`). So a FIFO barcode file does work today and rejecting it would have been an unjustified regression. D1 is the right call and its stated reason is accurate.

Ordering is unchanged too: `check_restartable_input` sits where the `exists()` it replaced sat, so the `--passthrough` case-fold collision check (1.ix, `cli.rs:903-919`) still runs after the existence check, and the `cli.input` loop still runs after the `--demux` block.

#### L7 — a third multi-read input surface exists that neither the plan nor D1 mentions

V10 claims "**Every** input-bearing CLI surface reaches a guard", narrowed by D1 to `cli.input` + `cli.passthrough`. There is one more: an adapter FASTA given as `-a file:adapters.fa` / `-a2 file:…`.

`adapter::parse_adapter_spec_inner` (`src/adapter.rs:313-316`) calls `read_fasta_adapters(path)` on every invocation, and the file is read once per call — but the call happens more than once: `Cli::validate` runs `parse_adapter_specs_quiet` for the `-a2` #369 pre-validation, and `setup_trimming` runs `parse_adapter_specs` again, **per pair**. So a FIFO adapter FASTA is read 2+ times and would block or come back empty on the second.

Out of #379's scope (the issue is about read input), and I would not add a guard for it — an adapter FASTA from a stream is not a thing anyone does. The point is the documentation claim: V10's "every" should be narrowed the way D1 narrowed it, so the next reader does not treat V10 as proof that no surface is left.

**Exact fix**: append to D1 in `PLAN.md` §12, or add a D4:

> **D4 — an adapter FASTA (`-a file:x.fa`) is also read more than once** (`adapter.rs:314`, once in `Cli::validate` for the `-a2` pre-check and again per pair in `setup_trimming`) and is **not** guarded. Deliberate: out of scope for #379, and a streamed adapter FASTA is not a real invocation. V10's enumeration is `cli.input` + `cli.passthrough`, not "every path-valued flag".

Priority: **Low** (documentation accuracy).

### 2.5 Efficiency

No concerns. `detect_input_format` gains one `fstat` on an already-open handle plus, for a regular file, one `open` + one `lseek` + a ≤4-byte `read`. §6's accounting is accurate and A4 (not a hot path) holds — the call count is per input file, not per record. The `[u8; 4]` probe buffer allocates nothing. `read_filled` replacing a single `read` is free in the common case (one syscall either way).

The declined optimisation in §6 (inspect the gzip branch's existing second handle) is the right call for the reason given: it would split the predicate across two sites and make it conditional on the branch taken.

### 2.6 Structure and style

**Comments, against `CLAUDE.md`'s "one line, state the fact":** mostly compliant, and the ones that are two lines earn it by pinning an ordering constraint that the code cannot show. Specifically good:

- `src/format.rs:317-318` — "Before the first `read`: a FIFO with a live writer but no data yet blocks there, and a second open with no writer blocks forever." Two facts, each of which is why the statement is *there* rather than three lines down. Keep.
- `src/main.rs:765-766` — "`sanity_check_any`, not `FastqReader::sanity_check`, so the passthrough path reaches the #379 restartability guard." Names a non-obvious call choice. Keep. (Its claim is what H1 disputes, but the comment itself is the right shape.)
- `src/format.rs:346-347` — the reworded gzip comment is 2 lines where the original was 4, and drops the now-false "not required for correctness". Exactly what §5.1 step 7 asked for.

Two I would trim:

- `src/format.rs:297-298` — "An inherited offset is the exact signal. A device that cannot report a position falls through to the byte comparison below." The second sentence narrates the next four lines, which is the thing `CLAUDE.md` names. One line does it: `// An inherited offset is the exact signal; anything else defers to the bytes.`
- `.github/workflows/ci.yml:30-31` — "without a bound that is a 6-hour job, not a red X" is the evidence, not the fact. `# A blocked FIFO open would otherwise burn GitHub's 360-minute ceiling.` One line. That said, this file's neighbouring comments run to 20 lines (`ci.yml:56-70`), so by "read like the surrounding code" the current form is in keeping — take or leave it.

Doc comments (`///`) are longer, and that is fine: they carry preconditions (`reopen_restarts`' FIFO rule, `is_non_restartable_type`'s reason for taking `Metadata`) that a caller cannot derive. The one gap is the one M3 names — `first`'s length contract is missing from the doc that documents everything else about it.

**Naming.** `is_non_restartable_type`, `not_restartable_message`, `reopen_restarts`, `read_filled`, `check_restartable_input` all read correctly and none of them lies about what it does. `reopen_restarts` returning `bool` (not `Result<()>`) keeps the "unexpected error is `Err`, not `false`" distinction that Q3 turns on. Good.

**Duplication.** `not_restartable_message` shared between `cli.rs` and `format.rs` is the right factoring — the plan wanted the wording unable to drift, and `pub(crate)` is the minimum visibility that achieves it. The `#[cfg(unix)]` / `#[cfg(not(unix))]` pair on `is_non_restartable_type` correctly avoids the bare `use FileTypeExt` that A7 identified as a hard break for a non-Unix library consumer.

One structural note, no action: on non-Unix, `is_non_restartable_type` returns `false` for everything, so a Windows named pipe would reach `reopen_restarts`, whose `File::open` can block. There is no `cfg(unix)` elsewhere in `src/`, CI is Ubuntu + macOS only, and the crate does not claim Windows support — so this is theoretical. Worth knowing if Windows is ever targeted.

---

## 3. Fixes recommended

None applied (per the override in my brief). In priority order:

| # | Priority | File:line | Fix |
|---|---|---|---|
| H1 | **High** | `src/main.rs:767-771` | Reject a uBAM `--passthrough` file before `sanity_check_any`, plus an integration test. Exact code in §2.1. Without it, `--passthrough <ubam>` regresses from a pre-output error to `stream did not contain valid UTF-8` with three empty output files on disk. |
| M1 | Medium | `src/format.rs:510-513` | Split `RecvTimeoutError::Timeout` from `Disconnected` in `bounded()` so a panicking body is not reported as a hang. |
| M2 | Medium | `src/cli.rs:488-489` | Keep the `NotFound` message for `NotFound` only; propagate other `io::Error`s with `Failed to stat input '{}'`. |
| M3 | Medium | `src/format.rs:288, 303-304` | Document `first`'s 1..=`PEEK_LEN` contract and `debug_assert!` it, so an out-of-contract call cannot silently answer `false` (or vacuously `true`). |
| L1 | Low | `CHANGELOG.md` | "before anything is written" → "before any output file is written" (`ensure_output_dir` at `main.rs:315` precedes the Reopen-class rejection). |
| L2 | Low | `src/format.rs:210` or `:222-226` | Either add a `Socket` variant with its own noun, or drop "or socket" from the `NotRestartable::Pipe` doc so code and message agree. |
| L3 | Low | `tests/integration_non_restartable_input.rs` | Add the macOS `--passthrough /dev/stdin` test so the layer-0 edit has a test that fails when it is reverted. |
| L5 | Low | `tests/integration_non_restartable_input.rs:80-85` | Drain stderr before the timeout `panic!` so the hang message carries what the binary said. |
| L7 | Low | `PLAN.md` §12 | Add D4: the `-a file:x.fa` adapter FASTA is a third multi-read surface, deliberately unguarded; narrow V10's "every input-bearing CLI surface" accordingly. |
| L8 | Low | `src/format.rs:297-298`, `.github/workflows/ci.yml:30-31` | Trim two 2-line comments to 1 line each. Optional. |

No Critical findings. H1 is the only one I would block a merge on.

---

## 4. What is right, and worth not losing in a refactor

Recording these so a later change does not undo them by accident:

1. **The layer ordering inside `detect_input_format` is load-bearing** — `fstat` before the first `read` (`:319-324` before `:327`), and `n == 0` before the second open (`:330` before `:334`). Moving either breaks a documented case, and the second one silently: an empty file would start reporting as non-restartable. D3's strengthened `detect_empty_file_errors` is what catches that, and it catches it because it asserts the *absence* of the new message, not just `is_err()`.
2. **`reopen_restarts` interrogates the second handle, not the first.** A future reader will be tempted to "simplify" it to a seekability probe on the handle already open. That is #379's own proposal and it accepts macOS row B. The doc comment at `:282-284` says so; keep it attached to the code.
3. **`bounded` + `run_bounded` are the reason a dead guard fails instead of wedging.** The wedge is the failure mode this design is most exposed to, and it is the one that costs 6 hours of CI rather than showing a red X. Do not replace `run_bounded` with `Command::output()`, and do not move the write-only FIFO `open` onto the main thread.
4. **Layer 1 is what makes a writer-less FIFO diagnosable at all.** `fs::metadata` never blocks where `File::open` blocks forever — that asymmetry is the entire mechanism, and it is why `check_restartable_input` must stay a `stat` and never become an open.
