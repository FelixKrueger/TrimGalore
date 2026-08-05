# CODE REVIEW A — Reject non-restartable input (#379)

**Reviewer:** A (fresh context, no shared state with Reviewer B)
**Target:** uncommitted working-tree diff on `dev` @ `40f77ec` — `src/format.rs`, `src/cli.rs`, `src/main.rs`, `.github/workflows/ci.yml`, `CHANGELOG.md`, new `tests/integration_non_restartable_input.rs`, new `.config/nextest.toml`
**Spec:** `PLAN.md` v3 (§12 deviations D1–D3)
**Mode:** report only — no files modified. Every recommendation below carries the exact edit.

---

## 1. Summary

The three-layer guard is sound and, on everything I could construct, correct. Layer 1 covers **more** surfaces than the plan's own inventory claims: I verified rejection through `--paired` (R2 position), `--hardtrim5`, `--clump_only` and `--passthrough`, not just plain single-end. The §3.4 edge-case table reproduces exactly, row for row. The message renders with the intended wrapping and 4-space command indent. `read_filled` is correct. The empty-file message still wins for a regular file. `fmt` and `clippy --all-targets --release -D warnings` are clean (re-run, not taken on trust), and the suite is 467/467 under nextest with a 4.9 s slowest test.

One **High** finding: the layer-0 edit at `src/main.rs:770` (`FastqReader::sanity_check` → `sanity_check_any`) makes a **uBAM `--passthrough` file pass a check it previously failed**, so that run now fails *after* three output files exist instead of before any. I verified this against a purpose-built `HEAD` worktree binary. `PLAN.md` §3.2 asserts this arm is unreachable; §2.4 of the same document already records why it is not.

Two **Medium** findings, both about diagnostics in a change whose entire purpose is diagnostics: the unit-test bound reports an assertion failure as a hang (proven with a standalone repro), and layer 1 discards the `stat` errno so `EACCES`/`EMFILE` are reported as "Input file not found".

The claim I probed hardest and could **not** break: *no input that previously succeeded now fails*. I could not construct a counterexample.

### Verified claims (challenged, not accepted)

| Claim | Result |
|---|---|
| fmt + clippy clean | Re-ran both. Clean. |
| 375 unit + 92 integration, 0 failures | `cargo nextest run`: **467 tests, 467 passed**, 9.1 s total, slowest 4.87 s. Matches. |
| No test hangs | New integration file: 5 passed in 0.38 s. `format::tests` #379 rows: 0.01 s. No nextest LEAK reports despite three detached threads. |
| Mutation test fails in 10 s rather than hanging | Mechanism confirmed by inspection: with the type check dead, the writer thread has already returned, so `reopen_restarts`' second `File::open` blocks and `bounded` fires. Structurally sound. |
| `terminate-after = 4` (240 s) is a safe backstop | ~50× margin over the slowest observed test. Safe. |
| A4 — `detect_input_format` is on no hot path | Confirmed: every call site in `adapter.rs`, `specialty.rs`, `clump_only.rs`, `main.rs` is per-file setup. None per-record or per-chunk. |
| §3.4 edge-case table | Reproduced exactly for empty / 1-byte / `/dev/null` / `/dev/zero` / `/dev/urandom` / directory / symlink-to-FIFO / missing file. |
| Not-found strings preserved | Byte-identical to the `exists()` versions; `cli.rs:1764` (`--passthrough file not found`) still passes. |

---

## 2. Issues by area

### 2.1 Logic

**H1 — `sanity_check_any` on `--passthrough` accepts a uBAM the old check rejected, and defers the failure past output creation.** `src/main.rs:770`

`PLAN.md` §3.2 states: *"`--passthrough` + uBAM is rejected at `main.rs:253`, so `sanity_check_any`'s BAM arm is unreachable here and harmless."* That is false, and §2.4 of the same plan says why:

> the `--passthrough` + uBAM rejection at `main.rs:253` keys off `any_bam`, computed from `cli.input`, so the passthrough file's container format is never inspected anywhere either.

FASTQ R1/R2 with a uBAM passthrough leaves `any_bam == false`, so `main.rs:253` does not fire and `sanity_check_any(pt_path)` takes the `UnalignedBam` arm, opens a `BamReader`, reads one record and returns `Ok(())`. Measured, same fixture (`test_files/ubam_test.bam`) both sides:

```
HEAD (FastqReader::sanity_check)        working tree (sanity_check_any)
─────────────────────────────────      ──────────────────────────────────────────
Error: stream did not contain          … Passthrough: pt.bam → outw/pt_passthrough.fq
valid UTF-8                            Error: processing pair 1 of 1 (R1=r1.fq, R2=r2.fq)
                                       Caused by: stream did not contain valid UTF-8
files in output dir: (none)            files in output dir: pt_passthrough.fq
                                                            r1_val_1.fq
                                                            r2_val_2.fq
```

So the change moves this failure from "before any output file" to "after three output files", which is the exact shape §2.4 identifies as the thing to eliminate (*"it hangs with partial output on disk. Exactly the failure mode §1 says this plan removes"*). The configuration is odd — FASTQ mates with a uBAM index read — but it is reachable, it is a regression against `HEAD`, and the fix is small. **High**, not Critical: it still errors, and no wrong output is produced.

**L1 — `reopen_restarts` silently answers `false` for `first.len() > PEEK_LEN`.** `src/format.rs:304-306`

```rust
let want = first.len().min(again.len());
let n = read_filled(&mut file, &mut again[..want]).with_context(ctx)?;
Ok(n == first.len() && &again[..n] == first)
```

With `first.len() > PEEK_LEN`, `want` clamps to 4, so `n ≤ 4 < first.len()` and the function reports **not restartable** for a perfectly good regular file. Only the sole caller keeps this unreachable (`&peek[..n]`, `n ≤ 4`); the `.min()` makes the function *look* as though it handles a longer slice, which is the trap. The doc comment states the FIFO precondition but not the length one. Reject-direction only — it cannot wrongly accept. **Low**, but the mis-answer is silent, which is what makes it worth a line.

**L2 — `--demux <writer-less FIFO>` still hangs with no message, after the full output is on disk.** `src/cli.rs:942-944` (D1)

Measured:

```
$ mkfifo f5; trim_galore -o outd --demux f5 r1.fq      # bounded by alarm(8)
…
JSON report: outd/r1.fq_trimming_report.json
Trimming complete, starting demultiplexing procedure (based on 3' barcodes supplied as per file >f5<)
*** blocked in File::open until killed at 8 s ***
```

D1's rationale — *"`read_barcode_file` opens the barcode file once … a FIFO barcode file works today, and rejecting it would be an unjustified regression"* — is true **only with a live writer**. Without one, `demux::read_barcode_file` (`main.rs:1265`) blocks in `File::open` with no message, after trimming, FastQC and the JSON report have all completed. That is not a restartability violation, so D1 is right that layer 2's invariant does not apply; it is a *blocking-open* hang, which layer 1 would have closed for free. V10's "every input-bearing CLI surface reaches a guard" is falsified here. Pre-existing, not a regression. **Low** — but either take the guard or amend D1 to record the writer-less case.

**L3 — `-a2 file:<stream>` is a second-open surface and is unguarded.** `src/cli.rs:972` + `src/main.rs:1017`

`Cli::validate` quiet-parses `--a2` (which reads a `file:` FASTA), and `main.rs:1017` parses it again. That is a genuine two-open path on a user-supplied file, outside every layer. Confirmed empirically: with a three-shot FIFO writer, `--a2 file:<fifo>` completes; with a single-shot writer the second parse would see EOF. Nobody streams an adapter FASTA, so this is **Low** — but it is the second counterexample to V10, and V10's stated purpose was *"the next flag that takes a path should be caught by a test rather than by a reviewer"*.

**L4 — D2's `is_ok_and` fall-through re-opens A2's residual, narrowly.** `src/format.rs:299`

D2 is the right call and is strictly less aggressive than the plan, as advertised. Worth recording precisely what it costs: on any handle whose `stream_position()` **errors**, the byte comparison becomes the sole condition, and the byte comparison is exactly what A2 was retired for (a repeating 4-byte prefix defeats it). No real input reaches that state — `lseek` fails with `ESPIPE` for pipes, FIFOs and sockets, and all three are rejected by the type check at `format.rs:322` before `reopen_restarts` is called, leaving only regular files and character devices, for which `stream_position()` succeeds. So this is a documentation point, not a defect. Consider adding "A2's residual returns on this branch" to D2.

**No hole found in `read_filled`.** `Ok(0)` → EOF; `Interrupted` retried; every other error propagated with `.with_context`. The `Interrupted` arm can in principle spin forever, but that is `std::io::Read::read_exact`'s own contract, not a new risk. It does harden the pre-existing `n >= 3` gzip test as A8 claims.

**No hole found in the `n == 0` / new-check interaction.** For a regular file the empty check still precedes the second open, verified end to end (`Input file 'empty.fq' is empty`, no restartability text). D3's strengthened `detect_empty_file_errors` pins both halves.

### 2.2 Errors and error-path quality

**M2 — layer 1 discards the `stat` errno, so `EACCES`/`EMFILE`/`ELOOP` are reported as "Input file not found".** `src/cli.rs:488-489`

```rust
let meta =
    std::fs::metadata(path).map_err(|_| anyhow::anyhow!("{not_found}: {}", path.display()))?;
```

Measured with a file inside a `chmod 000` directory (`stat` → `EACCES`):

```
Error: Input file not found: locked/hidden.fq
```

This is *identical* to the old `exists()` behaviour, so it is not a regression — but it is the same wrong-blame the plan rules out one layer down. `PLAN.md` §10 Q3, on `reopen_restarts`:

> **`Err`** … `Ok(false)` would report `EMFILE` or `EACCES` as "your input is a pipe" — the same wrong-blame this plan exists to remove, one layer down.

Layer 1 now does precisely that, one layer up, at a site the change touched and where the `io::Error` was in hand. `reopen_restarts` honours Q3 correctly (`with_context`, path named). **Medium** on the strength of the plan's own principle, and because it is three lines.

**No mis-blame found in layer 2.** `detect_input_format`'s three fallible steps each carry distinct context (`Failed to open input file`, `Failed to stat input file`, `Failed to read from`), and `reopen_restarts` carries `Failed to re-open input file '…' for the restartability check` on both the open and the read, exactly as §4 specifies. A transient `EMFILE` on the second open surfaces as itself, not as "your input is a pipe".

### 2.3 Tests — can any of them hang CI, and do the bounds bound what they claim?

**M1 — `bounded()` reports an assertion failure inside the body as "blocked; it must fail, not hang".** `src/format.rs:510-513`

```rust
match rx.recv_timeout(std::time::Duration::from_secs(10)) {
    Ok(v) => v,
    Err(_) => panic!("{what} blocked; it must fail, not hang"),
}
```

`Err(_)` conflates `Timeout` with `Disconnected`. If the body panics — e.g. `detect_input_format` returns `Ok`, so `expect_err("a FIFO must be rejected")` fires — the sender drops, `recv_timeout` returns `Disconnected` **immediately**, and the reported failure is a hang. Standalone repro (`rustc -O --edition 2021`):

```
thread '<unnamed>' panicked at probe_bounded.rs:18:13:
a FIFO must be rejected: Err value not returned
thread 'main' panicked at probe_bounded.rs:10:19:
detect_input_format on a FIFO blocked; it must fail, not hang  [real recv error: Disconnected]
--- outcome after 203.334µs
```

So the one failure mode this harness exists to distinguish — hang versus reject — is reported backwards for the *other* regression, the one where the FIFO is silently accepted. The inner panic is still visible in libtest's captured output, so a careful reader recovers the truth; the headline does not. In §5.3's own terms this is the mirror image of "a blocked test looks like a slow job". **Medium.**

**Nothing here can hang CI.** Checked directly:

- `bounded` panics on the test thread; the blocked body thread is detached, never joined, and dies with the process. libtest exits via `process::exit` and does not wait on non-main threads, so neither `cargo test --release` nor nextest can wedge on it.
- `spawn_fifo_writer` returns a `JoinHandle` bound to `_writer` (a real binding, not `_`), which is dropped-and-detached at scope end, not joined. No self-deadlock: the writer's `O_WRONLY` open rendezvous with the reader's open, and it cannot exit before the reader's open has returned, so the ordering the test relies on is guaranteed rather than raced. `write_all` is `let _ =`, so the expected `EPIPE` is swallowed (rule 3).
- `run_bounded` (integration, `tests/…:60-91`) drains stderr on its own thread, so a chatty child cannot deadlock against the wait loop; the 20 s deadline `kill()`s then `wait()`s; the post-exit `recv_timeout(5 s)` degrades to an empty string rather than blocking. For the `sh -c "cat … | trim_galore …"` rows, `sh` reaps both pipeline members before exiting, so stderr reaches EOF.
- No nextest LEAK reports across the whole suite, despite three detached threads.
- All 10 temp-dir slugs are unique (`tg_format_fifo_detect`, `tg_format_fifo_stat`, `tg_format_reopen`, `tg_format_devfd`, `tg_format_devfd_repeat`, `tg_379_fifo_in`, `tg_379_stdin_plain`, `tg_379_stdin_gz`, `tg_379_passthrough`, `tg_379_control`) — rule 5 honoured, and the 19-slug audit recorded at `ci.yml:52-56` still holds.

**The bounds do not cover setup — which is safe here, but the comment says otherwise.** `src/format.rs:501-502` claims *"each bounds its body"*. `bounded()` wraps only the call under test; `make_fifo` and `spawn_fifo_writer` run unbounded on the test thread, and the integration tests bound only the child process, not `fresh_dir`/`make_fifo`. `PLAN.md` §5.3 rule 2 asked for the whole body *"setup included"*, because the measured wedge (Shape B, `child.wait()` on a writer) was in setup. It is fine as written — the writer's blocking `open` is on its own thread, and `mkfifo(1)`, `remove_dir_all`, `create_dir_all` and `fs::write` cannot block — but that is a property worth stating rather than a bound. **Low**, comment accuracy only.

**L8 — `spawn_fifo_writer` swallows an open failure, then the reader blocks and the panic blames the guard.** `src/format.rs:540`. `if let Ok(mut f) = …open(&p)` discards the error; if the write-only open ever fails, no writer arrives, the reader's `File::open` blocks and `bounded` panics with "blocked; it must fail, not hang" — pointing at the guard when the fault is in the harness. That is the class rule 4 guards against for `mkfifo`, unguarded one line later. **Low.**

**L5 — the `/dev/urandom` negative is nondeterministic at 2⁻³².** `src/format.rs:601`. Four random bytes could match; the plan chose this knowingly. Noting only so it is not mistaken for a deterministic control.

**Coverage note (informational).** `detect_rejects_dev_fd_when_the_leading_bytes_repeat` is the only test isolating the offset condition — v3's central design change — and it is `#[cfg(target_os = "macos")]`. On `ubuntu-latest` reverting `stream_position()` fails nothing. The matrix does include `macos-latest`, so CI covers it; there is just no redundancy, and §12 already records this as the Darwin-only leg.

### 2.4 CI

**L6 — `timeout-minutes: 30` may be tight on a cold cache.** `.github/workflows/ci.yml:32`. `lto = true` + `codegen-units = 1` means every one of the now-**ten** release test binaries gets its own full LTO link. Measured locally (M-series, warm registry, cold `target`): `cargo build --release` 1 m 12 s; `cargo test --release --no-run` a further **1 m 00 s wall / 175 s CPU**. Scaling 175 CPU-seconds to a 2-vCPU `ubuntu-latest` and adding a cold `cargo install cargo-nextest --locked`, the debug build, and both test runs puts a cold-cache job in the 15–22 minute band — a 1.4–2× margin, not the 10× the number implies. The bound's stated job is to turn a 360-minute wedge into a fast red, and the debug leg is already bounded at 240 s per test by nextest, so only a release-profile-only hang depends on this number; 45 buys the same protection with room. **Low.**

`.config/nextest.toml` is correctly placed (workspace-root `.config/`), not gitignored, and `[profile.default]` is the profile `cargo nextest run` uses with no `--profile` flag. Note for the record that `terminate-after` does **not** cover the `cargo test --release` step (`ci.yml:74`) — the job bound is the only backstop there, which is what §5.2 says.

### 2.5 Structure and style

Against `CLAUDE.md`'s "one line, state the fact" rule, the new comments sit inside a file whose existing register is far more discursive, and they are consistent with it. One exception:

`.github/workflows/ci.yml:30-31` — *"A test that opens a FIFO with no writer blocks instead of failing; without a bound that is a 6-hour job, not a red X."* The second clause is evidence and colour; `CLAUDE.md` puts that in the commit message. One line would do.

`.config/nextest.toml`'s two-line comment states the fact and earns its second line. The `read_filled` and `reopen_restarts` doc comments carry rationale, which is what doc comments in this crate do. The rewritten `format.rs:346-347` is a genuine improvement — four lines to two, and it now says the true thing.

`make_fifo` and `fresh_tmpdir`/`fresh_dir` are duplicated between `src/format.rs`'s test module and the integration file. Unavoidable across crate boundaries without a test-support crate, and `fresh_tmpdir` duplication is already the repo idiom in five places. No action.

The `REJECTION` constant in the integration file is the right call — it is the one substring both message variants share, so a wording drift fails in one place.

### 2.6 Efficiency

One `fstat` on an already-open handle plus, for regular files, one extra `open` + `stream_position()` + ≤4-byte `read` per `detect_input_format` call. Peak concurrent handles rise from 1 to 2 (3 on the gzip branch). Against the five-plus opens and full-file read per input this is unmeasurable, and A4 holds — no call site is per-record. Layer 1 is genuinely zero extra syscalls: `Path::exists()` is `fs::metadata(path).is_ok()`, so the replacement is the same call. Nothing to raise.

### 2.7 Byte-identity invariant

Nothing in the diff touches quality trimming, adapter alignment, filtering, naming or writing. Every new code path terminates in `bail!` or falls through unchanged. The plan's V2 result (6 outputs + 6 reports byte-identical to a `HEAD` build with a firing sentinel control) is the right test and I see no path by which it could be wrong. V3 (Perl 0.6.11) still needs the CI `validation` / `validation-ubam` run before merge, as §12 says.

---

## 3. Fixes recommended

Ordered by priority. None applied — reporting only, per the review brief.

### High

**F1 — reject a uBAM `--passthrough` file on its own detected format (H1).** `src/main.rs`, immediately after the existing block at `:253-258` and therefore before `ensure_output_dir` at `:315`:

```rust
    // `any_bam` above is computed from `cli.input`; the passthrough file is a
    // separate surface and needs its own format check.
    if let Some(ref pt) = cli.passthrough
        && matches!(detect_input_format(pt)?, InputFormat::UnalignedBam)
    {
        anyhow::bail!(
            "--passthrough is not supported with uBAM input in this release. \
             Either convert the uBAM to FASTQ via `samtools fastq` first, or drop --passthrough."
        );
    }
```

This restores the fail-before-any-output-file property, reuses the existing wording verbatim, and makes `main.rs:770`'s BAM arm genuinely unreachable — which is what §3.2 assumed. It also runs the #379 guard on the passthrough file earlier than `main.rs:770` does. Worth a regression test asserting no output file exists after a uBAM passthrough is rejected; `test_files/ubam_test.bam` is a ready fixture.

### Medium

**F2 — distinguish `Timeout` from `Disconnected` in `bounded` (M1).** `src/format.rs:510-513`:

```rust
        match rx.recv_timeout(std::time::Duration::from_secs(10)) {
            Ok(v) => v,
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                panic!("{what} blocked; it must fail, not hang")
            }
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                panic!("{what} panicked; see the panic above")
            }
        }
```

**F3 — keep the errno at layer 1 (M2).** `src/cli.rs:488-489`:

```rust
    let meta = std::fs::metadata(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            anyhow::anyhow!("{not_found}: {}", path.display())
        } else {
            anyhow::Error::new(e).context(format!("Cannot stat input file: {}", path.display()))
        }
    })?;
```

`cli.rs:1764` and the "Input file not found" wording are preserved for the `NotFound` case, which is the only case any test asserts on.

### Low

**F4 — pin `reopen_restarts`' length precondition (L1).** `src/format.rs`, before line 303, and add the precondition to the doc comment alongside the FIFO one:

```rust
    debug_assert!(
        first.len() <= PEEK_LEN,
        "reopen_restarts compares at most PEEK_LEN bytes"
    );
```

**F5 — either guard the `--demux` barcode file or amend D1 (L2).** The guard is one line, `src/cli.rs:942-944`:

```rust
            check_restartable_input(demux_file, "Barcode file not found")?;
```

My recommendation is to take it. D1 correctly observes that the restartability invariant does not apply, but the blocking-open hang does, and the working configuration it would break — a live-writer FIFO as a two-column barcode TSV — has no plausible user, whereas the hang leaves a completed trim run wedged with no message. If you decline, add the writer-less case to D1 so the deviation records what it costs, and drop V10's "every input-bearing CLI surface" to the two surfaces it actually covers.

**F6 — record `-a2 file:` as an unguarded second-open surface (L3).** No code change needed for this release; V10's enumeration should name it, or it will be rediscovered as a bug report.

**F7 — surface a failed writer open in `spawn_fifo_writer` (L8).** `src/format.rs:540`:

```rust
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .open(&p)
                .expect("FIFO writer must open");
```

A panic on this thread disconnects the channel; with F2 applied, `bounded` then says "panicked", not "blocked" — the two fixes compose.

**F8 — raise the CI job bound (L6).** `.github/workflows/ci.yml:32`: `timeout-minutes: 45`.

**F9 — trim the ci.yml comment to the fact (style).** `.github/workflows/ci.yml:30-31`:

```yaml
    # A FIFO test with no writer blocks rather than fails, so bound the job.
```

**F10 — correct `src/format.rs:501-502`** from "each bounds its body" to what is true: each bounds the call under test, and no setup step can block.

---

## 4. Plan-text corrections worth folding into §12

Not code, but they are claims a future reader will rely on:

1. **§3.2 is wrong** that `sanity_check_any`'s BAM arm is unreachable for `--passthrough`. §2.4 already contains the refutation.
2. **V10 is falsified** by two surfaces: `--demux` (D1, deliberate) and `-a2 file:` (not considered).
3. **D1's rationale is incomplete** — it addresses restartability but not the writer-less blocking open, which is the failure the plan is named after.
4. **D2 should note** that A2's repeating-prefix residual returns on the `stream_position()`-error branch, even though nothing real reaches it.
