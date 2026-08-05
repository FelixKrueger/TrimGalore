# PLAN — Reject non-restartable input with an accurate diagnostic

**Issue:** [#379](https://github.com/FelixKrueger/TrimGalore/issues/379)
**Branch:** `fix/reject-non-restartable-input` off `dev`
**Scope decision:** option **A** of three. This does **not** add pipe support — it replaces a misleading error with a correct one. See §10 for what B and C would have cost.

**Revision history**

| v | What changed |
|---|---|
| v1 | Byte-comparison probe alone. |
| v2 | Measurement showed that probe **hangs** on the commonest FIFO shape (`MEASUREMENTS_reopen_probe.md` §1 row E). Added a non-blocking `fstat` pre-filter in front of it. |
| **v3** | Dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). Both reviewers independently found the guard **misses `--passthrough` entirely** (§2.4) and that the FIFO test specification was self-contradictory and unbounded (§5.3). Adds a path-level check in `Cli::validate` that also closes v2's Medium residual; swaps the byte comparison for an **offset-or-bytes** predicate; folds in the now-merged shape of [#374](https://github.com/FelixKrueger/TrimGalore/pull/374) (§7). |

---

## 1. Goal

When the input cannot be re-read from the start, fail immediately with a message that says so and explains what to do. Today such input fails with one of three wrong outcomes:

| Input | Today |
|---|---|
| `<(cat plain.fastq)` | `doesn't seem to be in FastQ format (first line doesn't start with '@')` — `fastq.rs:485` |
| `cat reads.fq.gz \| … /dev/stdin` | `Failed to decompress first block of '…' for format detection` — `format.rs:236` |
| `mkfifo f; cat reads.fq > f &` | **hangs**, no message |

All three blame the user's data, or say nothing, for a limitation of ours.

Non-goal: making pipes, FIFOs, process substitution or `/dev/stdin` work. That is a capability change (§10).

---

## 2. Context

### 2.1 Why the input is read more than once

Trim Galore opens each input at least **five times** for a gzipped file, and reads the first million records twice:

| Pass | Site | Opens | Consumes |
|---|---|---|---|
| Format detection | `format.rs:208`, `:231` | 2 | 4 bytes, then one decompressed block |
| Sanity check | `fastq.rs:474` (`FastqReader::sanity_check`) | 1 | 1 record, then **drops the reader** |
| Adapter auto-detection | `adapter.rs:115` / `:268` via `open_sync_reader` | 1 | up to `MAX_SCAN_READS = 1_000_000` (`adapter.rs:84`) |
| Trim pass | `format.rs:256`/`:273` factories | 1 | everything |

The true count is higher and rising: `detect_input_format` runs at *every* one of those sites; `cli.input[0]` passes through it twice before any reader is built (`main.rs:178` → `:30`, then the all-inputs map at `:187`); and after #374 the bare `FastqReader::open` / `open_threaded` / `sanity_check` entry points each perform **two** opens rather than one, because `sniff_gzip` opens the path to read three magic bytes and then the stream is opened again (§7).

The exact number does not matter. What matters is that it is greater than one, at every layer, and that each pass opens the path afresh. That is fine for a regular file and wrong for anything else.

### 2.2 The invariant that breaks, and the guards that do not test it

`format.rs:227-230` states the assumption explicitly:

> Re-open from byte 0 (cheaper than seeking and sidesteps any interaction between `MultiGzDecoder` and the already-consumed prefix; the alternative `PeekReader` wrapper in PLAN §5 step 2.3 is a documented future optimization, **not required for correctness**).

The real invariant is **re-opening the path yields the same bytes from the start**. Measured on Darwin (kernel 25.5.0, macOS 26.5.1, arm64 — `MEASUREMENTS_reopen_probe.md` §1), with the first handle still open, as `format.rs:231` leaves it:

| # | input | open1 | open2 | restarts |
|---|---|---|---|---|
| A | regular file | `@MS2` | `@MS2` | **yes** |
| B | `/dev/stdin` redirected from a regular file | `@MS2` | `_rea` | no |
| C | `/dev/stdin` from a pipe | `@MS2` | `_rea` | no |
| D | `<(cat reg.fq)` process substitution | `@MS2` | `_rea` | no |
| E | named FIFO, **writer exited** after writing | `@MS2` | *blocked >3 s* | **hangs** |
| F | named FIFO, writer still open | `@MS2` | `_rea` | no |

Row A is the positive control: the comparison is capable of returning true. Both reviewers reproduced row E independently.

**A seekability probe** (`lseek(fd, 0, SEEK_CUR)` on the *first* handle) accepts row B — that fd *is* seekable — which then fails with the original misleading message. On Darwin `/dev/stdin` is `fd/0`, and opening `/dev/fd/N` is `dup(N)`, so the second open shares the file offset. This is what killed the approach #379 proposes. Note the distinction from §3's actual predicate, which interrogates the **second** handle (§2.2.2).

**A file-type check cannot replace it** either: row B fstats as a *regular file* (§2.2.1).

**And the byte-level check cannot stand alone**, because row E shows its own second open blocking indefinitely. POSIX `open(O_RDONLY)` on a FIFO waits for a *writer*; the reader we are already holding does not satisfy it. `mkfifo f; cat reads.fq > f &` is the commonest FIFO shape, and its writer exits as soon as the payload fits the 64 KiB pipe buffer. A guard whose purpose is to print a message would print nothing and wedge.

Rows E and F differ only in whether the writer is still alive when the second open happens — i.e. in payload size against the pipe buffer. That is a race, not a property of the input, and §5.3 is written around it.

### 2.2.1 The file-type checks

`File::metadata()` on the handle `detect_input_format` already holds (`MEASUREMENTS_reopen_probe.md` §2):

| # | input | `fstat` on the open handle |
|---|---|---|
| A | regular file | regular |
| B | `/dev/stdin` from a regular file | **regular** |
| C | `/dev/stdin` from a pipe | FIFO/pipe |
| D | process substitution | FIFO/pipe |
| E, F | named FIFO | FIFO/pipe |

`is_fifo()` names every case where a second open would block (E) or consume stream bytes (C, D, F), and reports *regular* for B — which is why it guards the second-open check rather than replacing it. `std::os::unix::fs::FileTypeExt` is std; no new dependency. Reviewer B measured that `File::metadata` follows a symlink to its target's type, so a symlink to a FIFO is caught and there is no TOCTOU window on a handle-level check.

`is_socket()` **is** included, on Reviewer B's recommendation. v2 excluded it as unreachable (`File::open` on a socket path failing first), but that claim was never measured — the measurement file has no socket row — and on Linux a socket reached through `/proc/self/fd/N` returns `ENXIO`, which §10 Q3's `Err` decision would surface as a bare OS error. The check costs one `||` and removes a dependency on an unmeasured claim.

`is_char_device()` is **not** included: `/dev/null` is a character device and *is* restartable, so rejecting the class would be wrong. Character devices fall through, and §3.2 records where each one lands.

The argument against a type check in v1 ("it is platform-dependent") was aimed at the wrong target. The *behaviour* is platform-dependent: Linux resolves `/dev/stdin` through `/proc/self/fd/0` and may genuinely restart, Darwin `dup`s and cannot. Only a behavioural check accepts what works on each platform and rejects what does not. A pipe, by contrast, is a pipe on both — so the type filter is safe to layer underneath.

### 2.2.2 The second-layer predicate: offset, then bytes

v2 compared the leading bytes of the two opens. Reviewer B showed that predicate has a hole and a sharper replacement.

The hole: with a single `read` and a prefix comparison — the natural reading of v2's wording — a `/dev/fd/N` over a **≤4-byte** file is *wrongly accepted*, because the second read returns 0 bytes and an empty prefix compares equal. Measured:

| file size | n1 | n2 | prefix compare | length-aware compare |
|---|---|---|---|---|
| 4 B | 4 | 0 | **true — wrongly accepted** | false |
| 1 B | 1 | 0 | **true — wrongly accepted** | false |
| 5 B | 4 | 1 | false | false |

The replacement: `Seek::stream_position()` on the **second** handle. `0` means a genuine fresh open; non-zero means the open inherited an offset, which is exactly what rows B, D and F are. Measured by Reviewer B:

| input | second open `stream_position()` | offset verdict |
|---|---|---|
| regular 4 B / 6 B / 40 B | `Ok(0)` | restartable |
| `/dev/fd/N` over 4 B | `Ok(4)` | **not restartable** (the byte compare missed this) |
| `/dev/fd/N` over 6 B / 40 B | `Ok(4)` | not restartable |
| `/dev/null`, `/dev/zero` | `Ok(0)` | restartable |
| `/dev/urandom` | `Ok(0)` | restartable |

This is std, needs no `libc`, and needs no second `read`. It is **not** the seekability probe §2.2 rejects: that one asks the *first* handle "are you seekable" and accepts row B; this asks the *second* "did you start at zero", which is the property that actually differs.

**Decision: take both conditions, OR'd** — not restartable if `stream_position() != 0` **or** the length-aware byte comparison fails. Rationale:

- The offset condition is exact, so it retires A2's pathological-repeating-prefix residual and closes the ≤4-byte hole.
- Keeping the byte condition preserves v2's treatment of a *streaming* character device: `/dev/urandom` is `Ok(0)` on offset but differs on bytes, so it stays rejected. Dropping the byte half would silently change that to the not-FASTQ error.
- Neither condition can fire on a genuine regular file (measured, both directions).

### 2.3 Placement, and what the v2 inventory got wrong

`format::detect_input_format` is the function every **`cli.input`** path passes through: the two factories (`format.rs:256`, `:273`), `main::sanity_check_any` (`main.rs:30`) and its siblings at `:51`/`:68`, the all-inputs map at `main.rs:187`, `clump_only.rs:265`/`:402`/`:747`/`:910`/`:949`/`:950`, `specialty.rs:128`/`:180`, and `main.rs:1328`/`:1900`/`:2050`/`:2051`.

v2 claimed this covered "every input path" and stamped it **Verified**. That was verified against `cli.input`, not against every input-bearing CLI surface. Corrections, all confirmed against source:

| v2 claim | Correction |
|---|---|
| `specialty.rs:531` is a call site | It is inside `#[cfg(test)] mod tests` (module opens at `:508`) — a `read_fastq` helper, not a production path |
| `demux.rs:170` is downstream of a detect | It opens Trim Galore's **own trimming output** (name derived at `demux.rs:129`). Safe because we wrote it, not because a detect preceded it. This is also why #374's `sniff_gzip` deliberately avoids `detect_input_format` (§7) |
| `main.rs:1373`/`:1374`/`:1383` are all downstream of a detect | `:1373`/`:1374` are R1/R2 and are. **`:1383` is the `--passthrough` file and is not** (§2.4) |
| "the guard fires first in every flow" | False. See §2.4 |

The three `clump_only.rs` pairs do check out — `:290` is preceded by `:265`, and `:430`/`:432` by `:402`, on the same paths — and none of `clump_only`'s detect sites is on a self-created intermediate.

### 2.4 The hole both reviewers found: `--passthrough`

`cli.passthrough` is a separate `Option<PathBuf>` (`cli.rs:282`), not an element of `cli.input`, so it is in none of the inventories above. Its whole lifecycle:

| Step | Site | What it does |
|---|---|---|
| 1 | `cli.rs:884` | `if !pt.exists()` — a `stat`; a FIFO satisfies it |
| 2 | `main.rs:768` | `FastqReader::sanity_check(pt_path)` — **not** `sanity_check_any`, so no `detect_input_format`. One record, then drops the reader |
| 3 | `main.rs:1344` (parallel) / `main.rs:1383` (serial) | `FastqReader::open_threaded(p)` / `open(p)` |

The internal-invariant loop at `main.rs:1327` iterates `[input_r1, input_r2]` only. So the passthrough file is opened twice — four times after #374, which adds a sniff open to each bare entry point — with **no guard at any layer**. Consequences if this plan shipped without covering it:

- **`--passthrough <writer-exited FIFO>` still hangs with no message**, and by then `ensure_output_dir` has run and the R1/R2 writers exist — so it hangs *with partial output on disk*. Exactly the failure mode §1 says this plan removes.
- **`--passthrough <live-writer FIFO or /dev/fd/N>` desyncs**: step 2 consumes a record, step 3 resumes mid-stream, and the three-way header sync check at `parallel.rs:640` fires with a message about a truncated or desynced passthrough — blaming the user's data for a limitation of ours, one flag over.

Reviewer B found the same gap is structural rather than incidental: the `--passthrough` + uBAM rejection at `main.rs:253` keys off `any_bam`, computed from `cli.input`, so the passthrough file's container format is never inspected anywhere either.

---

## 3. Behavior

Two placements, deliberately. Both fail closed; each covers what the other cannot.

### 3.0 Layer 1 — `Cli::validate` (path level, never opens)

`cli.rs:933` already stats every input unconditionally (`if !path.exists()`), reached from `main.rs:166` before any dispatch branch. Reviewer A measured the decisive property: **`fs::metadata` on a writer-less FIFO returns promptly while `File::open` blocks.**

```
FIFO, no writer ever:   Path::exists = true
                        fs::metadata = Ok(fifo=true)
                        File::open   = *** BLOCKED >3s ***
regular file (control): fs::metadata = Ok(fifo=false);  File::open returned
```

So replacing that `exists()` with an `fs::metadata()` match costs **zero extra syscalls** and buys three things v2 could not:

1. It closes v2's Medium residual. `trim_galore fifo` with no writer blocked at the first open with no message, and v2 called that unfixable without opening every input `O_NONBLOCK`. `stat(2)` never opens, so the FIFO is named and rejected before anything can block.
2. It covers `--passthrough` (`cli.rs:884`) and the `--demux` barcode file (`cli.rs:926`) at the same site.
3. It fires earliest, so nothing is created before the message.

It does **not** replace layer 2: `Cli::validate` is bypassed by a library consumer calling e.g. `clump_only::clump_only_single` without a `Cli`, and a path stat is a TOCTOU proxy where layer 2's `fstat` is not. Two cheap checks that fail closed.

### 3.1 Layer 2 — `detect_input_format` (handle level)

The order is load-bearing:

1. Open the path (unchanged).
2. **New:** `fstat` the handle. If `is_fifo() || is_socket()`, bail with the §3.3 pipe message. This precedes the `read`, so a FIFO with a live writer and no data yet is rejected instead of blocking in `read`.
3. Read up to 4 bytes into `peek`, **looping** until the buffer is full or EOF and retrying `ErrorKind::Interrupted` (§8 A8).
4. `n == 0` → the existing empty-input error (unchanged).
5. **New:** open the path a second time. Not restartable if `stream_position() != 0`, or if a looped read of `n` bytes does not return exactly `n` bytes equal to `peek[..n]`. On not-restartable, bail with the §3.3 re-open message.
6. Everything from here is untouched: `@` → `FastqPlain`; gzip magic → decompress a block and probe for `BAM\1`; anything else → the existing "not recognised as FASTQ" error.

### 3.2 Layer 0 — route `--passthrough` through the guard

`main.rs:768`: `FastqReader::sanity_check(pt_path)` → `sanity_check_any(pt_path)`.

`sanity_check_any` dispatches through `detect_input_format` (`main.rs:30`) and then calls `FastqReader::sanity_check` for the FASTQ arms, so the existing check still runs. Two notes:

- **This does not undo #374.** Both `detect_input_format` and #374's `sniff_gzip` are content-based, so a `.bgz` passthrough file keeps working; the run simply gains one redundant content check.
- The empty-input message changes from `seems to be completely empty` (`fastq.rs:478`) to `is empty` (`format.rs:216`) for an empty passthrough file. Both are errors, both are accurate, and an empty passthrough would desync R1/R2 anyway.

~~`--passthrough` + uBAM is rejected at `main.rs:253`, so `sanity_check_any`'s BAM arm is unreachable here and harmless.~~

**Corrected in the fix round (§13, H1) — this was false, and §2.4 above already said why.** `main.rs:253` tests `any_bam`, computed from `cli.input`, so a uBAM *passthrough* file never reached it. `sanity_check_any`'s BAM arm therefore **accepted** the file (it opens a `BamReader` and returns `Ok`) where `FastqReader::sanity_check` had rejected it, and the run went on to die in the FASTQ reader with three empty output files on disk. Both code reviewers found this independently. The arm is genuinely unreachable now, because §13's fix adds a passthrough-specific format check beside `main.rs:253` — i.e. the claim was made true rather than merely asserted.

### 3.3 The messages

Two openings, one shared tail. Both contain the substring **`cannot be re-read from the start`**, so a single assertion covers either path.

Pipe or FIFO (layer 1, and layer 2 step 2):

```
Input 'X' is a pipe or FIFO, not a regular file, so it cannot be re-read from
the start.
```

v2 said "named pipe (FIFO)". Reviewer B pointed out that two of the three #379 reproductions — `<(zcat reads.fq.gz)` and `cat x | trim_galore /dev/stdin` — are **anonymous** pipes, so v2's noun was wrong for the majority case, in a plan whose premise is that we should not describe the user's input inaccurately. Nothing caught it, because V1 asserts only the shared tail.

Re-open mismatch (layer 2 step 5):

```
Input 'X' cannot be re-read from the start: opening it a second time did not
return to the beginning of the data.
```

Shared tail:

```
Trim Galore reads each input more than once — format detection, the initial
sanity check, and adapter auto-detection each open it independently — so it
requires a regular file.

Common causes are pipes, FIFOs, process substitution such as
`<(zcat reads.fq.gz)`, and /dev/stdin. Write the stream to a file first:

    zcat reads.fq.gz > reads.fq && trim_galore [options] reads.fq
```

"Common causes are" rather than "This affects", because the list is not exhaustive — a streaming character device reaches the same message (§3.4). Per `CLAUDE.md`, the reasoning lives here and in the commit message rather than in a source comment.

### 3.4 Edge cases

| Case | Handling |
|---|---|
| Empty regular file | Existing `Input file 'X' is empty` (`format.rs:216`), unchanged — a regular file passes both type checks and `n == 0` precedes the second-open check |
| File shorter than 4 bytes | Length-aware comparison (§2.2.2); a 1-byte regular file is restartable and falls through to the existing format errors. The same file via `/dev/fd/N` is now correctly **rejected** — v2 would have accepted it |
| Empty FIFO (writer opened, wrote nothing, exited) | **Behaviour change:** the pipe message, previously `is empty`. Both true; naming the pipe is more useful, and step 2 must precede the read to avoid blocking |
| FIFO with **no writer at all** | Rejected by layer 1 without opening. This was v2's Medium residual |
| Socket | Rejected by layer 2 step 2, or by layer 1 |
| `/dev/null` | Character device, `n == 0` → existing `is empty` ✓ |
| `/dev/zero` | Character device, offset 0 and bytes equal → falls through → existing not-recognised error ✓ |
| `/dev/urandom` | Character device, offset 0 but bytes differ → the **new** re-open message. Accurate (it genuinely cannot be re-read), and the tail no longer claims to be exhaustive |
| `/dev/stdin` **at an interactive terminal** | Character device; neither type check fires, and the first `read` blocks with no message. Residual — see §11 |
| Directory | `read` fails with `EISDIR` — unchanged from today |
| Second open fails with an unexpected `io::Error` | Propagated as itself with context, **not** reported as non-restartable — §10 Q3 |
| Regular file being concurrently appended | Restartable; out of scope, and no worse than today |

---

## 4. Signature

```rust
/// True iff opening `path` a second time returns to the beginning of the data.
///
/// Trim Galore opens each input several times (format detection, sanity check,
/// adapter auto-detection, the trim pass), so a path whose re-open does not
/// restart cannot be processed. `/dev/fd/N` shares the file offset on some
/// platforms even though it is seekable, which is why this interrogates the
/// second handle rather than probing the first for seekability.
///
/// Callers must reject FIFOs before calling this: a second open on a FIFO with
/// no live writer blocks indefinitely.
fn reopen_restarts(path: &Path, first: &[u8]) -> Result<bool>
```

Private to `format.rs`. Takes the bytes already read so the caller does not read twice from the first handle.

- `Ok(false)` means *opened successfully and did not restart* — a non-zero start offset, or a length-aware byte mismatch.
- An unexpected `io::Error` is `Err`, carrying `.with_context(|| format!("Failed to re-open input file '{}' for the restartability check", path.display()))` in the style of `format.rs:209`/`:213`. Without that context the improvement Q3 bought is only half-delivered: a bare `Too many open files (os error 24)` with no path is still a bad diagnostic.

The type checks stay **inline** — layer 1 in `Cli::validate`, layer 2 in `detect_input_format`. Each is two lines and layer 2 needs the open handle, which no helper has. §9 V7 exercises them through their callers.

---

## 5. Implementation outline

### 5.1 Source

1. **`Cli::validate` (layer 1).** Replace the `path.exists()` check at `cli.rs:933` with an `fs::metadata()` match: `Err` → the existing not-found error; `Ok(m)` where `m.file_type().is_fifo() || is_socket()` → the §3.3 pipe message; otherwise proceed. Apply the same at `cli.rs:884` (`--passthrough`) and `cli.rs:926` (`--demux` barcode file). No new syscalls — these sites already stat.
2. **`main.rs:768` (layer 0).** `FastqReader::sanity_check(pt_path)` → `sanity_check_any(pt_path)`.
3. **`detect_input_format` type check (layer 2 step 2).** After the `File::open` at `format.rs:208` and **before** the `read` at `:211`. `use std::os::unix::fs::FileTypeExt;` at module scope, `#[cfg(unix)]`-gated with a `#[cfg(not(unix))]` no-op (§8 A7) — the crate ships a `src/lib.rs`, so a bare `use` is a hard compile break for any non-Unix consumer.
4. **Loop the first read** at `format.rs:211` until 4 bytes or EOF, retrying `Interrupted`.
5. **Add `reopen_restarts`**, private, above `detect_input_format`, per §4.
6. **Call it** after the `n == 0` check at `format.rs:215-217` and before the `peek[0] == b'@'` branch at `:219`.
7. **Reword the comment at `format.rs:227-230`** — "not required for correctness" is now false in the general case and true only because the guard rejects the inputs where it would matter. One line.
8. **CHANGELOG** entry under `#### Fixes` in `Unreleased`. Reviewer B checked `docs/`, `README.md` and `--help` for `/dev/stdin`, `mkfifo` and `… | trim_galore` and found nothing, so the CHANGELOG is the entire documentation surface.

### 5.2 CI backstop — do this first, not last

Neither exists today. Both reviewers filed it, and Reviewer B gave the argument that makes it structural rather than hygiene: **if the layer-2 type check is ever accidentally dead** — placed after the `read`, or the predicate inverted — the FIFO test does not fail, it *hangs*, because the second open blocks (row E). The plan's central design decision fails, if it fails, in the one mode the plan itself calls hardest to notice.

```yaml
# .github/workflows/ci.yml — rust-tests job
timeout-minutes: 30
```

```toml
# .config/nextest.toml — no such file today
[profile.default]
slow-timeout = { period = "60s", terminate-after = 4 }
```

Verified: `grep -rn "timeout-minutes" .github/workflows/` returns nothing, there is no `.config/`, and the job runs `cargo nextest run` (`ci.yml:67`) **and** `cargo test --release` (`:74`) on both matrix legs. An escaped wedge therefore burns GitHub's default 360-minute ceiling on four job-steps. This is the one backstop that survives an implementer forgetting every rule in §5.3.

### 5.3 FIFO test construction — the part v2 got wrong

v2 asked for "both writer states — writer exited, and writer still open". **Drop that.** Reviewer A measured why it is not a specification anyone can implement:

```
Shape A  writer held open in-process for the whole call  -> message in 144 µs   deterministic
Shape B  child.wait() for the writer BEFORE the call     -> *** BLOCKED >3s ***  in SETUP
Shape C  spawn writer, then call immediately             -> message in 28 ms     deterministic
Shape D  v1's second open, writer gone                   -> *** BLOCKED >3s ***  (row E)
```

"Writer exited" is **unreachable as a pre-condition**: a writer's `open(O_WRONLY)` blocks until a reader opens, and the only reader is the call under test, so Shape B blocks in `child.wait()`. What row E observed was a writer exiting *between* the first and second open — which v2's pre-filter deletes by construction. And Reviewer B showed v2's own rule 2 ("payload under 64 KiB so the writer exits") *contradicted* the writer-still-open case, making it the scheduler race the measurement file identifies. With the pre-filter, both states traverse identical code: there is one test.

**The construction.** Hold the write end open on a spawned thread inside the test process:

```rust
let w = { let p = fifo.clone(); thread::spawn(move || {
    let f = std::fs::OpenOptions::new().write(true).open(&p).unwrap();
    std::thread::park();            // hold the writer open
    drop(f);
}) };
```

The write-only `open` rendezvous with our reader's `open`, so both complete. **It must not be on the main thread** — a write-only `open` on a FIFO blocks until a reader appears, so doing it before `detect_input_format` on the same thread self-deadlocks. That is the obvious way to write this test and it is a trap.

**Rules, each from a measured failure mode.** Every one of these replaced or corrected a v2 rule:

1. **A writer must open the FIFO**, or the *first* open blocks (§3.4). Spawn it before the call; its `open` blocks until ours succeeds, so there is no race.
2. **Bound the whole test body, setup included** — not "the call". v2 bounded the call; the measured wedge (Shape B) is in setup, outside that bound. Run the body on a thread and `mpsc::recv_timeout`. Do **not** rely on `timeout(1)`: absent on Darwin and on GitHub's `macos-latest` runner.
3. **The writer must tolerate `EPIPE`.** The guard bails and drops the reader before the writer's `write_all` lands, so an in-process writer must not `unwrap()` its write, and a spawned `cat` must not have `.success()` asserted — it is killed by `SIGPIPE`, payload-size-dependently. (This is the real reason v2's payload rule was needed; the reason v2 gave was already stale.)
4. **Assert the fixture is a FIFO.** Check `mkfifo`'s exit status and `assert!(fs::metadata(&f)?.file_type().is_fifo())`. If `mkfifo` fails, `cat x > f` creates a regular file, `detect_input_format` returns `FastqPlain`, and the test fails pointing at the guard when the fault is in the harness.
5. **One `fresh_tmpdir` slug per FIFO test.** `fresh_tmpdir` opens with `remove_dir_all`, and the debug leg is process-per-test under nextest, so a shared slug lets one process unlink another's FIFO while a writer is blocked on it. `ci.yml:52-56` records an audit that all 19 existing `temp_dir().join(…)` sites use unique slugs — do not be the first exception.
6. **Integration tests get their own recipe.** The repo idiom is `Command::new(binary()).output()` (`tests/integration_no_args_help.rs:21-24` and every `integration_clump_only*` case). `output()` reads both pipes to EOF and cannot be interrupted, so a child blocked in `open` holds them open and `output()` never returns. Use `spawn()`, keep the `Child`, and on timeout `kill()` then `wait()`. Give the writer child `Stdio::null()` and kill it in cleanup — a leaked child holding output handles trips nextest's `leak-timeout` (default 100 ms) and is reported LEAK.

### 5.4 Tests

**Unit, in `src/format.rs`'s existing `mod tests`** (which already has `fresh_tmpdir`):

- a regular plain FASTQ is accepted; a regular gzipped FASTQ is accepted (happy-path regression guards);
- a FIFO is rejected with the pipe message — **one** test, per §5.3;
- an empty regular file still produces the *empty* error, not either new message;
- `reopen_restarts(&regular_fq, b"@rea")? == true` — positive direction;
- `reopen_restarts(&regular_fq, b"XXXX")? == false` — negative direction, portable (see V7);
- `reopen_restarts(Path::new("/dev/urandom"), …)? == false` — negative direction against a real non-restartable source that cannot block;
- a `/dev/fd/N` over a 4-byte file is rejected, `#[cfg(target_os = "macos")]` — the case v2 would have accepted (§2.2.2).

**FIFO creation:** `std::process::Command::new("mkfifo")` (confirmed at `/usr/bin/mkfifo` on Darwin, present on both runners). Correcting v2's stated rationale, which both reviewers flagged: dev-dependencies are `serde_json`, `tempfile` **and `bstr`** (`Cargo.toml:50-56`), not "tempfile only", and `libc 0.2.184` is **already in `Cargo.lock`** transitively via `rustix`/`tempfile` — so `libc::mkfifo` would add one dev-dep line and *zero* new compiled crates. The honest comparison is "a `$PATH` dependency with a failure mode that masquerades as a guard bug (rule 4)" versus "one dev-dep line on a crate already in the graph, and an immediate errno". `Command` is chosen for needing no `unsafe`; the trade-off is close.

**Integration** `tests/integration_non_restartable_input.rs`: invoke the built binary with a FIFO, with `/dev/stdin` from a pipe, and with `--passthrough <fifo>`; assert non-zero exit and stderr containing `cannot be re-read from the start`. Rules 1–6 apply. Do **not** assert on `/dev/stdin` redirected from a *regular file* — that is row B, whose correct outcome is platform-divergent (§2.2.1).

---

## 6. Efficiency

Per input:

- **Layer 1:** zero extra syscalls. `cli.rs:884`/`:926`/`:933` already stat; `fs::metadata` replaces `exists()`, which is the same call.
- **Layer 2, FIFO or socket:** one `fstat` on an already-open handle, then bail. No extra open.
- **Layer 2, regular file:** one `fstat`, one extra `open`, one `stream_position()`, and a ≤4-byte `read`.

`detect_input_format` runs several times per input (§2.1), so the regular-file cost is one extra open per *call*, not per input; on the gzip branch a call now performs **three** opens (`:208`, the probe, `:231`). Against five-plus opens and a full-file read this is unmeasurable. No allocation — the probe buffer is a `[u8; 4]` on the stack.

An optimisation exists and is **not** taken: the gzip branch already re-opens at `:231`, so its handle's offset could be inspected for free, making the check zero-cost there and paid only on the plain branch. That would split the check across two sites and make it conditional on the branch taken. One `open` is not worth that. (Reviewer A's ALT-4 — memoising detection per path — would remove the redundant opens altogether and make this moot, but it collides with #374 and is a larger change than this plan wants.)

---

## 7. Integration

**Reads:** the input path. **Writes:** nothing.

**Order:** layer 1 in `Cli::validate` (`main.rs:166`), before any dispatch; layer 2 inside `detect_input_format`, before any reader is constructed. Only one call site changes (`main.rs:768`).

**Downstream:** an input that previously failed later with a wrong message — or hung with no message — now fails immediately with a right one. No input that previously *succeeded* can newly fail: restartability is exactly the property every existing successful run already had (§8 A1).

### 7.1 #374 — no longer "not at all"

v2 asserted the two changes do not interact. That was true of v2's shape and is now wrong in three specific ways. Benjamin Demaille (`@BenjaminDEMAILLE`) pushed `16d47ca` and `a293a1b` on 2026-08-04, and the PR now takes the content-sniff route:

1. **`FastqReader`'s bare entry points gain an open each.** `is_gz_filename` becomes `sniff_gzip`, and `open` / `open_threaded` / `sanity_check` each call it with a path, so each performs its own `File::open` for three magic bytes before the stream is opened. That is a **new instance of the very invariant this plan guards** — on a pipe the sniff consumes three bytes and the stream open loses them; on a writer-exited FIFO the sniff's open blocks. `sniff_gzip` swallows every error by design (`.is_ok_and(…)`, "any failure here answers 'not gzip'"), so it cannot report either. This makes §2.4's unguarded passthrough path strictly worse — four opens, no guard — and is an additional argument for layer 0 and layer 1, which stop the input before `FastqReader` ever sees it.
2. **`detect_input_format` itself is untouched**, so layer 2's edit site is conflict-free. But #374 does edit `src/format.rs` (`:258`, `:278`, `:296`) and deletes 37 lines from `src/main.rs` — the duplicate `open_sync_reader`/`open_threaded_reader`, i.e. two of §2.3's inventory entries (`main.rs:51`, `:68`) cease to exist. **Every `main.rs` line number in this plan shifts.** Benjamin's own description already cites the post-merge numbers for the passthrough sites: `main.rs:743`/`1319`/`1358`, against this plan's `768`/`1344`/`1383`.
3. **`sniff_gzip` deliberately avoids `detect_input_format` because it bails on empty input**, since `demux.rs:170` reads back a trimmed output that is legitimately empty when every read was filtered. Record this: it is a standing reason the guard cannot later be pushed down into `FastqReader`, and it independently confirms §2.3's correction about `demux.rs:170`.

**Merge order: #374 first.** It is complete, green (468 tests) and an external contributor's; this plan should absorb the churn rather than impose it. Then rebase, correct the `main.rs` references, and re-verify §2.4's three sites. Note that layer 1 lives in `cli.rs`, which #374 does not touch at all — so the part of this change that closes the most cases cannot conflict with it either way.

---

## 8. Assumptions

- **A1.** Every currently-working input is restartable. Follows from all the passes in §2.1 already depending on it. Nothing in `tests/`, `.github/` or `justfile` uses `/dev/stdin`, `mkfifo` or process substitution; **re-verified** by grep with a positive control, and with a tighter pattern than v2 used — a bare `<(` matches Rust tuple syntax (`Vec<(…)>`) throughout `tests/`. Reviewer A adds a sharpening: on Linux, row B genuinely restarts, so `trim_galore /dev/stdin < reads.fq` works on `ubuntu-latest` today and continues to, because the offset check returns 0 there. No regression on either platform. Exactly true only once A8 lands.
- **A2.** *Retired.* v2 assumed comparing ≤4 bytes was sufficient and accepted a pathological repeating-prefix residual. The offset condition (§2.2.2) is exact, so the residual is gone.
- **A3.** `File::open` on a FIFO with no writer blocks. **Verified** three times — `MEASUREMENTS_reopen_probe.md` §1 row E, and independently by both reviewers.
- **A4.** `detect_input_format` is on no hot path — a handful of times per input file, never per record. Confirmed by the call-site inventory (§2.3) and independently by Reviewer B.
- **A5.** For a **regular** file the `n == 0` empty check still precedes the second-open check, so empty-file behaviour is unchanged. Fixed by construction in §3.1 steps 4–5. For a FIFO the type check now wins instead (§3.4, deliberate).
- **A6.** `fstat` names a FIFO as a FIFO on both target platforms. **Verified on Darwin; reasoned on Linux.** The POSIX `S_IFIFO` guarantee covers *named* FIFOs (rows E, F); rows C and D on Linux go through `/proc/self/fd/N`, a different mechanism. Reviewer A's kernel-level reading: re-opening a pipefs inode yields `S_IFIFO`, and `fifo_open`'s partner-wait is gated on `!is_pipe`, so an anonymous pipe does not block on re-open while a named FIFO does — matching row E on both platforms. V9 makes CI prove this rather than assume it.
- **A7.** Unix-only. `FileTypeExt` does not exist on Windows, and there is no `cfg(unix)` anywhere in `src/`, `tests/` or `build.rs` today. Handled rather than flagged: §5.1 step 3 gates it, because the crate ships a `src/lib.rs` and a bare `use` would be the first hard compile break for a non-Unix consumer.
- **A8.** *New.* Both reads must return everything asked for, or the comparison is wrong in both directions: a short second read wrongly rejects a good regular file (NFS, FUSE, `EINTR` surfacing as `ErrorKind::Interrupted` rather than being retried), and a prefix reading of the comparison wrongly *accepts* `/dev/fd/N` over a ≤4-byte file. Handled by the looped reads in §3.1 steps 3 and 5, not assumed. The same loop hardens the pre-existing `n >= 3` gzip test at `format.rs:226`, which today misclassifies a gzip file on a short first read.
- **A9.** *New, and the load-bearing one.* No input-bearing CLI surface bypasses the guard. v2 assumed this implicitly and was wrong (§2.4). Now covered by layer 1, which sees every path `Cli::validate` sees, and asserted by V10 rather than by inspection.

---

## 9. Validation

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | The three #379 reproductions say the right thing, and **each stops saying its own wrong thing** | Per row: `<(cat plain.fastq)`, `cat reads.fq.gz \| … /dev/stdin`, `mkfifo f; cat small.fq > f &` | All: non-zero exit, stderr contains `cannot be re-read from the start`. Per-row negative control — row 1 must not contain `doesn't seem to be in FastQ format` (`fastq.rs:485`); row 2 must not contain `Failed to decompress first block` (`format.rs:236`); row 3 must return promptly, since today it produces no stderr at all. v2 asserted only row 1's string for all three, which passes vacuously on rows 2 and 3 — it would have passed on unpatched `dev` |
| V2 | Regular files are unaffected | Full `cargo test`; plus a trim of `10K_150bp.fastq.gz` and its uncompressed twin | All pass; outputs byte-identical to `dev`. **Prove the comparison can fail** by diffing against a sentinel-appended copy |
| V3 | Byte-identity against Perl 0.6.11 holds | CI `validation` and `validation-ubam` jobs | Green. The guard adds no transformation, but this is the invariant that must not move |
| V4 | The empty-file error still wins | `trim_galore` on a 0-byte regular file | Stderr contains `Input file` … `is empty` (`format.rs:216`), **not** either §3.3 message. Not `seems to be completely empty` — that is `fastq.rs:478`, pre-empted at `main.rs:30` |
| V5 | **No FIFO test hangs**, and nothing can hang for long | New tests bounded per §5.3 rule 2; plus §5.2's `timeout-minutes` and nextest `terminate-after` in place | Complete in seconds. Structural, not lucky: for a FIFO the second open never happens, and if that regresses the CI bound turns a 6-hour wedge into a fast failure |
| V6 | The writer-exited FIFO — the case that hangs on `dev` today | `mkfifo f; cat small.fq > f & trim_galore f` | The pipe message, promptly. On `dev` this hangs at `fastq.rs:474`; v1 of this plan would have hung inside the guard |
| V7 | Every guard can fail, in both directions | Layer 1: `fs::metadata` on a writer-less FIFO → rejected without opening; on a regular file → proceeds. Layer 2 type check: FIFO → pipe message; regular → passes through. Second-open check: `reopen_restarts(&regular, b"@rea")` → `true`, `(&regular, b"XXXX")` → `false`, `/dev/urandom` → `false`, and (macOS only) `/dev/fd/N` over 4 bytes → `false` | All exercised. v2 specified `/dev/stdin` from a pipe as the negative case; both reviewers showed that is not constructible — a unit test cannot control its own fd 0, it blocks on a terminal, it never reaches the check anyway because the type filter catches it first, and on `ubuntu-latest` the negative branch would run **zero** times. Record in §9 that `/dev/fd/N` is not a portable substitute (dup on Darwin, fresh open on Linux) and that this predicate has no portable real-world negative reproduction |
| V8 | uBAM input still detected | Existing uBAM tests, plus a real `.bam` and a `bgzip`-ed FASTQ | Formats unchanged; every layer is format-agnostic |
| V9 | The Linux mechanism is measured, not assumed | Integration test asserts `fs::metadata(…).file_type().is_fifo()` on the descriptor it is about to hand the binary, before invoking it | Green on `ubuntu-latest` **and** `macos-latest`. The first CI run becomes the Linux measurement for A6, and a platform surprise names itself instead of appearing as a message mismatch. Row B stays unasserted (§5.4) |
| V10 | **Every input-bearing CLI surface reaches a guard** | Enumerate them — `cli.input`, `cli.passthrough`, the `--demux` barcode file — and assert each is rejected for a FIFO. Not "exercise one of them" | All rejected. This is the row whose absence hid §2.4: nine rows of input-shape testing could not see a whole flag being unguarded, and the next flag that takes a path should be caught by a test rather than by a reviewer |

V5, V6, V7 and V10 are the ones that would quietly not hold: V5/V6 because a hanging test looks like a slow job, V7 because a check only ever observed passing is not a check, V10 because it tests the plan's placement claim rather than its mechanism.

---

## 10. Questions or ambiguities

**Resolved before writing (scope):** three targets were costed and **A** chosen.

- **A — reject accurately** (this plan). ~70 lines plus tests, across `cli.rs`, `main.rs` (one line) and `format.rs`. Fixes the diagnostic, not the capability.
- **B — pipes when `-a` is given.** Requires threading a reader rather than a path through `format.rs`, `fastq.rs`, `bam.rs` and `main.rs`, folding `sanity_check` into the trim pass, and forbidding adapter auto-detection on streams. Benjamin's #374 review independently reached the same conclusion — "whoever takes it will want to hand a reader downstream instead of a path" — and confirmed a `bool` verdict is neither a step toward it nor an obstacle.
- **C — full single-pass front end.** All of B plus buffering the `MAX_SCAN_READS` window: ~200 MB at 65 bp, to be reconciled with `--memory` and `--clumpify`.

**Note for whoever attempts B:** this guard is in the way by construction. All three layers reject the exact inputs B would enable, and the second-open check consumes ≤4 bytes from a stream before rejecting it. They must come out, not be worked around.

**Open (assumption taken):**

1. **Should the message name workarounds per input type?** Taken: one generic remedy in a shared tail, with only the first sentence varying (§3.3). Enumerating FIFO vs `/dev/fd` vs `/dev/stdin` fixes is longer and each reduces to "write it to a file".
2. **How much of the stream should the second-open check compare?** Taken: none of it, primarily — the start offset is exact and length-insensitive (§2.2.2). The ≤4-byte comparison is retained only as a second reject condition, to keep streaming character devices rejected.
3. **`Ok(false)` vs `Err` when the second open fails.** **`Err`**, changed in v2 and kept. `Ok(false)` would report `EMFILE` or `EACCES` as "your input is a pipe" — the same wrong-blame this plan exists to remove, one layer down. With `.with_context` (§4) so the path is named.
4. **Whitelist or blacklist the file type?** Taken: blacklist (`is_fifo() || is_socket()`). Reviewer A costed the whitelist (`is_file() || is_char_device()`), which would also give sockets and directories an accurate message instead of a raw `errno`. Declined: it is a behaviour change for input classes nobody has complained about (block devices), and it converts "unknown file type" from *fall through and let the format check decide* into a hard error. Conservative is right for a diagnostic-only fix.
5. **Compare `(st_dev, st_ino)` across the two handles?** Declined, and recorded because it is the obvious third idea. On Darwin `open("/dev/fd/N")` is a `dup`, so the second handle reports the same device and inode. The offset is the only thing that differs — i.e. Q2.

**Critical:** none.

---

## 11. Self-Review

**Logic.** Three layers test one invariant, at the three points where an input can enter: `Cli::validate` (every path the binary is given), `detect_input_format` (every path a library consumer detects), and the one call site that had neither. Each covers what the others cannot — the path stat cannot be TOCTOU-proof but never blocks; the handle `fstat` is exact but requires an open; the offset check is the only thing that catches a `dup`-ed regular file. The ordering constraints are pinned in §3.1: type check before the `read`, empty check before the second open.

**What review changed, three times over.**

*Before v1:* I intended a seekability probe, which is what #379 proposes. Measurement killed it — row B is seekable and does not restart.

*Manual review of v1:* the byte comparison alone **hangs** on row E, the commonest FIFO shape and one of the three input classes the message names. v1 recorded this only as a test-hygiene risk. It was a design defect, and the fix was the type check v1 had dismissed on an axis that did not matter.

*Dual agent review of v2:* two independent reviewers both found, as their Critical, that the guard **misses `--passthrough` entirely** — and that v2 stamped the claim "Verified" when it had been verified against `cli.input` only. Both also found the FIFO test specification unbuildable: v2 asked for a writer state that cannot exist, bounded the wrong scope, and contradicted itself between two adjacent rules. Neither finding touched the mechanism; both touched claims I had marked as settled. The pattern is the same each time — a check verified against a narrower set than it names.

**Edge cases.** §3.4, expanded from six rows to twelve: the two new type-check classes, three character devices with three different outcomes, the terminal case, the ≤4-byte `/dev/fd/N` case v2 would have accepted, and the FIFO-with-no-writer case v2 called unfixable.

**Efficiency.** Zero extra syscalls at layer 1; one `fstat` and (regular files only) one open per `detect_input_format` call at layer 2. The available optimisation is declined in §6 with a reason.

**Integration.** One call site changes. #374's interaction is now specified rather than denied (§7.1), with a merge order and the reason.

**Remaining risks.**

- *Medium:* `trim_galore /dev/stdin` **at an interactive terminal** is a character device, so no type check fires and the first `read` blocks with no message. Distinguishing a tty needs `isatty`, i.e. `libc`. It is a plausible thing for someone who has just read #379 to try, and `/dev/stdin` is named in our own remedy text — so this is the residual most likely to generate a follow-up issue. Out of scope for A; worth a sentence in the CHANGELOG entry.
- *Medium:* the FIFO tests remain the most delicate part of the change. §5.3's six rules and §5.2's CI bound are what keep them honest; the bound is the only one that survives an implementer skipping the rules.
- *Low:* three layers is more surface than one. Mitigated by all of them emitting the same two messages and by V7 exercising each in both directions.
- *Low:* layer 1 is bypassed by a library consumer that never builds a `Cli`. That is what layer 2 is for.
- *Very low:* an exotic filesystem where a regular file does not restart. Would already be broken today.

---

## 12. Implementation notes

Implemented 2026-08-04 on `dev` (uncommitted at time of writing). `.github/workflows/ci.yml`, `.config/nextest.toml` (new), `src/cli.rs`, `src/format.rs`, `src/main.rs`, `CHANGELOG.md`, `tests/integration_non_restartable_input.rs` (new) — 340 insertions, 19 deletions.

### Deviations from the plan

**D1 — the `--demux` barcode file is *not* checked.** ~~§3.0 and §5.1 step 1 said to apply layer 1 at `cli.rs:926` as well.~~ **Reversed in the fix round (§13).** The reasoning below is correct about *restartability* and was silent about *blocking*: Reviewer A measured that a writer-less FIFO barcode file blocks in `File::open` with no message **after** trimming, FastQC and the JSON report have all completed. Layer 1 closes that for one line, and the working configuration it removes — a live-writer FIFO used as a barcode TSV — has no plausible user. The check is now applied, so layer 1 covers all three surfaces.

The original reasoning, retained because it is why the *restartability* argument alone did not justify the guard: `read_barcode_file` (`demux.rs:31-34`) opens the barcode file **once** and streams it with a `BufReader`; there is no second open. So the invariant this guard enforces genuinely does not hold there — the guard is justified by the blocking open, not by restartability.

**D2 — a `stream_position()` error falls through instead of propagating.** §3.1 step 5 implies the offset is always available. A character device that cannot report a position would have become a new hard error for a class that today falls through to the format checks, so the implementation only treats a *successfully read* non-zero offset as decisive (`is_ok_and`) and otherwise defers to the byte comparison. This is strictly less aggressive than the plan and preserves existing behaviour for exotic devices.

*Cost, per Reviewer A:* on any handle whose `stream_position()` **errors**, the byte comparison becomes the sole condition — and the byte comparison is exactly what A2 was retired for, so A2's repeating-prefix residual returns on that branch. Nothing real reaches it: `lseek` fails with `ESPIPE` for pipes, FIFOs and sockets, and all three are rejected by the type check before `reopen_restarts` is called, leaving only regular files and character devices, for which `stream_position()` succeeds.

**D3 — the pre-existing `detect_empty_file_errors` test was strengthened**, not just added to. It asserted only `is_err()`, which would have passed if the empty file had started being reported as non-restartable. It now asserts the empty message is present *and* the restartability message is absent (V4).

### One thing the plan did not anticipate

With layer 1 in place, the binary **never opens** a FIFO given on the command line. A test harness that starts a writer (`cat x > fifo &`) therefore leaves that writer blocked in `open` forever, because no reader ever arrives. §5.3 rule 6 already said to kill the writer child rather than wait on it — but for the EPIPE reason, not this one. The first smoke-test script hit exactly this and hung on its own `wait`, not on the binary. The integration tests avoid it entirely by using writer-less FIFOs, which is also the harder case: it is the one that used to hang forever. Layer 2's unit test still needs a live writer, because `detect_input_format` must get past `File::open` to reach the `fstat`.

### Verification

| Check | Result |
|---|---|
| V1 — three #379 reproductions | All three now emit `cannot be re-read from the start`. Per-row negative controls in the integration test: the plain-pipe row asserts absence of `doesn't seem to be in FastQ format`, the gz row absence of `Failed to decompress first block` |
| V2 — regular files unaffected | **6 output files and 6 trimming reports byte-identical** to a `HEAD` worktree build across SE-gz, SE-plain, PE `--cores 1`, PE `--cores 4`. Both comparisons carry a sentinel negative control that fired correctly |
| V4 — empty file | `Input file 'X' is empty`, restartability message absent |
| V5 — no test hangs | Full suite 375 unit + 92 integration, 0 failed. New tests complete in 0.01 s (unit) and 0.92 s (integration) |
| V5 — **mutation check** | With `is_non_restartable_type` forced to `false`, `detect_rejects_a_fifo_with_the_pipe_message` **fails in 10.01 s** with `blocked; it must fail, not hang`. This is the check the plan cared most about: a dead guard fails loudly instead of wedging for 6 hours |
| V6 — writer-exited FIFO | Rejected promptly. Also the writer-*less* FIFO, which v2 recorded as an unfixable residual |
| V7 — every guard, both directions | Layer 1 (`stat` on a writer-less FIFO vs a regular file), layer 2 type check, `reopen_restarts` positive + negative, `/dev/urandom`, and two macOS `/dev/fd` cases |
| V9 — no accidental platform conditionality | Verified on Darwin only; the `#[cfg(target_os = "macos")]` rows are the two `/dev/fd` cases. CI provides the Linux leg |
| V10 — every input surface | `cli.input` and `cli.passthrough` both covered by integration tests; see D1 for the third |
| `cargo fmt --all -- --check` | Clean (two of my own hunks reformatted; `git status` confirms fmt did not stray) |
| `cargo clippy --all-targets --release -D warnings` | Clean |
| V3 — Perl 0.6.11 byte-identity | **Not run locally** — needs the CI `validation` / `validation-ubam` jobs. V2's result is strong evidence it holds (no writer path touched), but this is the invariant that must not move, so it needs the CI run before merge |

### Test added beyond the plan

`detect_rejects_dev_fd_when_the_leading_bytes_repeat` — a file whose bytes 4..8 repeat bytes 0..4. The byte comparison alone calls it restartable, so only the start-offset condition rejects it. Without this, nothing isolated the offset check: the ≤4-byte case in §5.4 is caught by the length-aware comparison too, so §2.2.2's central design change had no test that would fail if it were reverted. This is the test that retires A2.

### Iteration log

1. Source + tests written; `cargo build --release` clean first time.
2. Smoke script hung on its own `wait` for blocked FIFO writers (see above). Not a defect in the binary — all six cases had already produced the right message. Rewritten to use writer-less FIFOs.
3. `--quiet` used in the V2 harness is not a Trim Galore flag; all eight baseline runs failed silently because the loop redirected stderr. Caught by the sentinel negative control, which reported `CONTROL FAILED`. Re-run without it: 8/8 exit 0.
4. zsh `nomatch` aborted the md5 loop on `*.fq` with no matches; `NULL_GLOB`.
5. Report comparison reported three files as `DIFFERS` — `diff: /dev/fd/11: Operation not permitted`, the sandbox blocking process substitution. Re-run with temp files: all six identical. A tooling failure that reads as a real difference.

### Follow-ups

- **#374 merge order.** §7.1 recommends merging Benjamin's PR first; this work was implemented against `dev` @ `40f77ec` and will need the `main.rs` line references re-checked after that lands. Nothing in this change touches the lines #374 edits.
- **The tty residual** (§11) is in the CHANGELOG entry as a known case rather than fixed.

---

## 13. Fix round — dual code review + coverage audit

`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` (independent, no shared state) and `COVERAGE.md`. Coverage verdict was **COMPLETE** (62 items, 0 MISSING, 0 PARTIAL). Both reviewers filed the same High finding and the same three Mediums. All findings taken.

### H1 (High) — `--passthrough <uBAM>` regressed. Both reviewers, independently.

Measured before/after on purpose-built binaries:

| | pre-change | after the layer-0 swap | after this fix |
|---|---|---|---|
| message | `stream did not contain valid UTF-8` | same, wrapped in `processing pair 1 of 1` | `--passthrough is not supported with uBAM input…` |
| on disk | output dir, no files | **three empty files** | **no output dir at all** |

Root cause: §3.2's justification was false and §2.4 contained the refutation. Fixed beside `main.rs:253` — Reviewer A's placement, chosen over B's because it precedes `ensure_output_dir` (`main.rs:315`) and makes `main.rs:770`'s BAM arm genuinely unreachable. Defended by `passthrough_ubam_is_rejected_before_output`, which asserts the message, the *absence* of `valid UTF-8`, and that no output file exists.

### Mediums, all three from both reviewers

- **M1** — `bounded()`'s `Err(_)` conflated `Timeout` with `Disconnected`, so a *panicking* body was reported as `blocked; it must fail, not hang`. Reviewer A reproduced it firing in 203 µs with the wrong headline. That is the mirror image of the failure the harness exists to identify, in a change about not misdescribing failures. Now split, with a distinct message per variant.
- **M2** — layer 1's `map_err(|_| …)` reported `EACCES`/`EMFILE`/`ELOOP` as `Input file not found`. Not a regression (identical to the `exists()` it replaced) but the same wrong-blame §10 Q3 rules out one layer down, at the one site in the diff that had the errno in hand. `NotFound` keeps its wording; everything else propagates with `Cannot stat input file`.
- **M3** — `reopen_restarts`' `.min()` clamp answered `false` forever for `first.len() > PEEK_LEN` and `true` vacuously for `0`. Unreachable from the sole caller, but the `.min()` made the function look as though it handled a longer slice. Contract documented and `debug_assert`ed.

### Also taken

| Source | Change |
|---|---|
| A | `--demux` barcode file now guarded — see the revised **D1** |
| A | `timeout-minutes` 30 → **45**. A measured `lto = true` + `codegen-units = 1` across ten release test binaries at 175 CPU-seconds locally; scaled to a 2-vCPU runner with a cold cache that is a 15–22 minute job, so 30 was a 1.4–2× margin, not the 10× it looked like |
| A | `spawn_fifo_writer` now `expect`s its open. A silent failure there left the reader blocking, and the panic blamed the guard. Composes with M1: the failure now reads "panicked", not "blocked" |
| A | Corrected the test-module comment: each test bounds *the call under test*, not its setup — and no setup step can block |
| B | `NotRestartable::Socket` added with its own noun. The doc said "pipe, FIFO or socket" while the message said only "pipe or FIFO" — §3.3's own standard, applied to us |
| B | CHANGELOG: "before anything is written" → "before any output **file** is written". `ensure_output_dir` (`main.rs:315`) precedes the layer-2 sites, so the Reopen class does leave an empty directory |
| B | `run_bounded` drains stderr **before** the timeout `panic!`, so the failure that matters most carries the binary's own output |
| B | `passthrough_dev_stdin_is_rejected` (macOS) — B's delete-one-thing check found layer 0 was the only edit with no test that fails when reverted, and it was the edit H1 shows to be harmful |

### D4 — a third multi-read surface, deliberately unguarded

Both reviewers found it independently: an adapter FASTA given as `-a file:x.fa` / `-a2 file:…` is read by `adapter::parse_adapter_spec_inner` once in `Cli::validate` (the `-a2` #369 pre-check) and again per pair in `setup_trimming`. That is a genuine two-open path on a user-supplied file, outside every layer. **Not guarded:** out of scope for #379, and a streamed adapter FASTA is not a real invocation.

**V10's wording is narrowed accordingly.** It claimed "every input-bearing CLI surface"; it covers `cli.input`, `cli.passthrough` and (now) the `--demux` barcode file, and it does *not* cover `-a file:`. Both reviewers also noted the enumeration is by test rather than by mechanism, so a future path-valued flag added unguarded still fails nothing — recorded rather than solved.

### Verification after the fix round

| Check | Result |
|---|---|
| Full suite | 375 unit + **95** integration (8 in the #379 file), 0 failed |
| `cargo fmt --all -- --check`, `clippy --all-targets --release -D warnings` | Clean |
| H1 | Rejected up front with the right message and **no output directory** |
| Byte-identity | Re-checked against the retained `HEAD` baseline after the `main.rs`/`cli.rs` edits: SE-gz and PE `--cores 4` outputs identical, sentinel control firing (1 914 872 vs 1 914 882 bytes) |
| V3 | **Still CI-gated.** Perl 0.6.11 byte-identity needs the `validation` / `validation-ubam` jobs |

### Harness note for the next session

Two more false negatives, both caught by controls rather than by inspection: zsh does **not** word-split unquoted variables, so `set -- $pair` gave one word and a byte-identity loop compared zero files while reporting success; and a negative control comparing two *unrelated* outputs reported "equal" because both md5 inputs were empty. Compare a sentinel-modified copy of a real file, never two different files.
