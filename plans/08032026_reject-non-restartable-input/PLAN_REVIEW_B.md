# PLAN_REVIEW_B — Reject non-restartable input with an accurate diagnostic

**Reviewer:** B (independent, fresh context)
**Target:** `plans/08032026_reject-non-restartable-input/PLAN.md` revision **v2**
**Supporting evidence reviewed:** `MEASUREMENTS_reopen_probe.md`
**Source read:** `src/format.rs`, `src/main.rs`, `src/fastq.rs`, `src/clump_only.rs`, `src/specialty.rs`, `src/demux.rs`, `src/cli.rs`, `Cargo.toml`, `.github/workflows/ci.yml`, `tests/`
**Own measurements:** taken on this machine, kernel Darwin 25.5.0 / macOS 26.5.1, arm64, rustc from `~/.cargo/bin`. Probe sources and outputs quoted inline. Nothing that opens a FIFO was run.

**Verdict:** the design is right and the v2 fix is the correct fix. One completeness claim in §2.3 is false and leaves the bug reachable on one flag; the FIFO test rule set contradicts one of its own test cases; and one §9 row is not implementable as written. All three are fixable in the plan without changing the approach.

---

## 1. Logic review

### 1.1 The two-layer guard: is it sound?

Yes, and the ordering argument in §2.2/§3 is correct. I re-derived it independently:

- The invariant the code depends on is *re-opening the path yields the same bytes from the start* (`format.rs:208` + `format.rs:231`, plus `fastq.rs:474` and `adapter.rs:115`/`:268` each opening afresh). §2.1's inventory of that is accurate — I confirmed all four rows against source.
- A type check alone cannot work, because row B (`/dev/stdin` from a regular file on Darwin) fstats as regular.
- A byte comparison alone cannot work, because the second `open` on a writer-less FIFO blocks. This is the v2 discovery and it is a genuine design defect in v1, not a test-hygiene issue. §11's account of that is honest.
- Therefore fstat-then-bytes, with the fstat before the `read`. Correct, and §3's step ordering pins it.

I also confirmed the two subtler claims that hold the design up:

- **fstat on the handle, not a path stat** (A6) — right, and it also disposes of symlinks: `File::metadata` follows to the target's type. Measured: `symlink → /dev/null` reports `chr=true`, so a symlink to a FIFO reports FIFO. No TOCTOU window on the type check.
- **The `n == 0` empty check must precede the probe** (A5) — verified against source: `format.rs:215-217` is the `Input file '{}' is empty` bail, and §3 puts the probe after it.

### 1.2 The hole: `--passthrough` never reaches `detect_input_format`

§2.3 states:

> The direct `FastqReader::open` calls in `specialty.rs:531`, `clump_only.rs:290`/`:430`/`:432`, `demux.rs:170` and `main.rs:1373`/`:1374`/`:1383` are all downstream of a `detect_input_format` on the same path, so the guard fires first in every flow. **Verified** by call-site grep.

`main.rs:1383` is not R1 or R2. It is the `--passthrough` file:

```rust
// src/main.rs:1382-1385
let mut reader_passthrough = match passthrough_input {
    Some(p) => Some(FastqReader::open(p)?),
    None => None,
};
```

and its `--cores > 1` / `--clumpify` twin at `src/main.rs:1343-1346` (`FastqReader::open_threaded(p)`) is not in the inventory at all. Nothing detects that path:

- `src/main.rs:187` maps `detect_input_format` over **`cli.input`** only; `cli.passthrough` is a separate `Option<PathBuf>` (`cli.rs:282`).
- `src/main.rs:768` sanity-checks it with `FastqReader::sanity_check(pt_path)?` — **not** `sanity_check_any`, which is the wrapper that calls `detect_input_format` (`main.rs:29-30`).
- The internal-invariant loop at `src/main.rs:1328` iterates `[input_r1, input_r2]` only.
- `cli.rs:884` checks only `pt.exists()`, which a FIFO satisfies.

So `trim_galore --paired --passthrough <(zcat I1.fq.gz) R1 R2` is unchanged by this plan: it reaches `FastqReader::sanity_check` → `FastqReader::open` (`fastq.rs:474`) and produces the old `doesn't seem to be in FastQ format` message; a writer-exited FIFO hangs there, and hangs again at 1383/1344. That is the exact bug #379 reports, left reachable on one flag, in a plan whose stated placement rationale is "the one place every input passes through".

Corroborating observation, out of scope but it shows the gap is structural rather than a one-off: the `--passthrough` + uBAM rejection at `main.rs:253` keys off `any_bam`, which is computed from `cli.input`, so the passthrough file's own container format is never inspected anywhere.

Fix is one line — `sanity_check_any(pt_path)?` at `main.rs:768`, or an explicit `detect_input_format(pt_path)?` beside it — plus a §9 row. See C1.

### 1.3 The rest of the §2.3 inventory

Two more entries are wrong, harmlessly, but they are two of five items in a list stamped "Verified":

- `specialty.rs:531` is inside `#[cfg(test)] mod tests` (module opens at `specialty.rs:508`) — a `read_fastq` test helper, not a production path.
- `demux.rs:170` opens `trimmed_file`, which is Trim Galore's **own trimming output** (name derived at `demux.rs:129`), not a user input. It is safe because it is a regular file we just wrote, not because a `detect_input_format` preceded it.

The three `clump_only.rs` entries do check out: `:290` is preceded by `:265`, and `:430`/`:432` by `:402`, on the same paths. I also checked every `detect_input_format` call site in `clump_only.rs` (`:265`, `:402`, `:747`, `:910`, `:949`, `:950`) — all are on user inputs, none on a self-created intermediate, so the guard cannot fire on a temp file Trim Galore made itself.

### 1.4 Comparison semantics: the probe can silently pass

§4 says "return whether the slices are equal"; §3.2 says "Probe compares `peek[..n]` against the same length from the second open". Neither pins what happens when the second read is **shorter**. Read as a common-prefix compare — a natural reading of "the same length" — the probe accepts the case it exists to reject. Measured, single `read` of `n1` bytes on the second handle, `/dev/fd/N` over a regular file (the row-B shape at small sizes):

| file size | n1 | n2 | prefix compare | length-aware compare |
|---|---|---|---|---|
| 4 B | 4 | 0 | **true — wrongly accepted** | false |
| 1 B | 1 | 0 | **true — wrongly accepted** | false |
| 5 B | 4 | 1 | false | false |
| 8 B | 4 | 4 | false | false |

The required predicate is `n2 == first.len() && buf[..n2] == first`. Blast radius is small — an input of ≤4 bytes fails downstream anyway — but this is the only place the guard can pass unnoticed, and nothing in §5 step 5 exercises `n2 < n1`. §3.2's row "a 1-byte file is restartable and falls through to the existing format errors" is right for a *regular* 1-byte file (measured: n1=1, n2=1, equal) and wrong for the same file reached through `/dev/fd/N`.

### 1.5 Short reads make A1 false in the corner

A single `read` also means a genuine short read on a restartable file is reported as non-restartable. `Read::read` is permitted to return fewer bytes than asked, and `File::read` surfaces `EINTR` as `ErrorKind::Interrupted` rather than retrying — which, under §10 Q3's v2 decision, becomes a hard `Err`. Rare (4 bytes, local FS), realer on NFS/FUSE. Looping the read until `first.len()` or EOF costs three lines and makes A1 ("no input that previously succeeded can newly fail") exactly true instead of nearly true.

### 1.6 Message accuracy

§3.1's FIFO opening is:

```
Input 'X' is a named pipe (FIFO), so it cannot be re-read from the start.
```

Two of the three #379 reproductions — `<(zcat reads.fq.gz)` and `cat x | trim_galore /dev/stdin` — are **anonymous** pipes. Both fstat as `S_IFIFO` (MEASUREMENTS §2 rows C and D), so this message fires for them and asserts something untrue about the user's input. A plan whose opening premise is that the current message "blames the user's data for a limitation of ours" should not ship a wrong noun, and nothing catches it: V1 asserts only the shared tail `cannot be re-read from the start`, which both messages contain by design (§3.1). "is a pipe or FIFO, not a regular file" covers all three shapes.

### 1.7 Termination, in both directions

I looked for inputs that terminate in neither direction and found one the plan does not list:

- **FIFO, no writer ever** — blocks at the first `open`. Listed (§3.2, §11 Medium). Correct and correctly scoped out.
- **Character device that blocks on `read`** — `trim_galore /dev/stdin` at an interactive terminal fstats as a character device, so the pre-filter does not fire, and the **first** `read` at `format.rs:211` blocks with no message. This is a plausible thing for a user who has just read #379 to try, and `/dev/stdin` is one of the forms §3.1's own remedy text names. Same class as the FIFO-with-no-writer residual, and §3.2's character-device row currently claims char devices are handled "both correct". Distinguishing a tty needs `isatty`, i.e. `libc`, so listing it as a residual is the right answer — but it should be listed.
- Everything else terminates: directories fail at the `read` with `EISDIR` (unchanged), block devices fail at `open` with `EACCES` for a normal user, `/dev/zero` and `/dev/urandom` return promptly.

### 1.8 Character devices: §3.2's row is incomplete

Measured on this machine:

| input | fstat | n1 | byte probe | resulting error |
|---|---|---|---|---|
| `/dev/null` | char | 0 | not reached | existing `is empty` ✓ |
| `/dev/zero` | char | 4 | **true** | existing not-recognised ✓ |
| `/dev/urandom` | char | 4 | **false** | **new** `cannot be re-read from the start` |
| symlink → `/dev/null` | char (followed) | 0 | not reached | existing `is empty` ✓ |

§3.2 says character devices produce "existing empty / not-FASTQ errors, both correct" and names only `/dev/null` and `/dev/zero`. A streaming character device gets the new message, whose tail names "pipes, FIFOs, process substitution … and /dev/stdin" — none of which it is. Harmless, but the row claims a completeness it does not have. The silver lining is in §3.4 below: `/dev/urandom` is the portable negative control the plan needs.

---

## 2. Assumptions

| # | Verdict | Note |
|---|---|---|
| A1 | **Qualified** | Holds for `cli.input`. The promise it underwrites — "the guard fires first in every flow" — is false for `cli.passthrough` (§1.2). Also only approximately true until the short-read loop lands (§1.5). Its grep half checks out: no `/dev/stdin` or `mkfifo` in `tests/`, `.github/`, `justfile`. The `<(` half of that grep matches Rust generics (`Vec<(…)>`) throughout `tests/`, so re-run it with a tighter pattern before restating it as verified. |
| A2 | **Superseded** | The pathological repeating-prefix case A2 discusses is not the realistic hole; §1.4's short-second-read case is, and A2 does not mention it. Both vanish under the Alt-1 offset check. |
| A3 | **Verified** | MEASUREMENTS §1 row E, plus POSIX. Correctly promoted from unverified assumption to the reason the type check exists. |
| A4 | **Verified independently** | 18 production `detect_input_format` call sites, all per-file. Nothing per-record. |
| A5 | **Verified against source** | `format.rs:216` is the `is empty` bail and precedes everything. V4's string is now correct; `fastq.rs:478` is indeed the other one (`seems to be completely empty`), and `fastq.rs:485` is `doesn't seem to be in FastQ format` — so V1's negative assertion targets a string that really exists. |
| A6 | **Verified on Darwin, asserted on Linux** | Acceptable: the ubuntu integration test is self-verifying. The important half — fstat on the held handle, not a path stat — is right, and my symlink measurement confirms `File::metadata` resolves to the target type. |
| A7 | **Verified** | `grep -rn "cfg(unix)\|cfg(windows)\|cfg(target_os" src/ tests/ build.rs` returns nothing. Unix-only is consistent with the tree as it stands. |

### Unstated assumptions worth surfacing

1. **That a `read` returns everything asked for.** §1.5. It is also the shape the existing code at `format.rs:211` already has, so the plan inherits it rather than introducing it — but the new code is where it produces a *new wrong message*.
2. **That `is_socket()` is unreachable.** §2.2.1 asserts `File::open` on a socket path fails with `ENXIO`, cited to MEASUREMENTS §2 — whose table has rows A–E and no socket row, even though `probe2.rs` prints `sock=`. I tried to measure it; the sandbox blocks `bind(2)`, so I could not settle it either. See I6: including the check costs one `||` and removes the dependency on an unmeasured claim.
3. **That the character-device set is `{/dev/null, /dev/zero}`.** §1.8.
4. **That a unit test can observe its own `/dev/stdin`.** §3.4 / C3.
5. **That "reap it, do not `wait` on it while holding the reader"** (§5 step 7 rule 2) is actionable. Reaping *is* waiting; the sentence needs to say what to do instead (kill, then wait).

---

## 3. Validation sufficiency (§9)

The nine rows are well chosen and §9's closing paragraph correctly identifies V5/V6/V7 as the fragile ones. Three ways the change could still ship broken with everything green:

### 3.1 A hang is not bounded by anything mechanical

`.github/workflows/ci.yml` has **no `timeout-minutes` on any job** — I grepped all workflows, zero hits. GitHub's default is 360 minutes, on both matrix legs. The debug leg is `cargo nextest run` (ci.yml:66-67), and nextest's default `slow-timeout` only *warns*; it does not terminate. So V5 ("Complete in seconds") is a human observation, and the only thing standing between a wedged test and a 6-hour job is every future test author remembering §5 step 7 rule 3.

This matters more than it looks, because of the chain: **if the fstat pre-filter is accidentally dead** — placed after the `read`, or the predicate inverted — the writer-exited FIFO test does not fail, it *hangs*, because the byte probe's second open blocks (MEASUREMENTS §1 row E). The plan's central v2 design decision therefore fails, if it fails, in the one mode the plan itself calls hardest to notice. Two one-line backstops that do not depend on discipline:

```yaml
# .github/workflows/ci.yml, rust-tests job
timeout-minutes: 30
```

```toml
# .config/nextest.toml  (no such file today)
[profile.default]
slow-timeout = { period = "60s", terminate-after = 4 }
```

### 3.2 Rule 3 does not work for the integration test

The repo's idiom is `Command::new(binary()).output()` on the test's main thread, with no bound — `tests/integration_no_args_help.rs:21-24` and every case in `tests/integration_clump_only*.rs`. `Command::output()` reads both pipes to EOF and cannot be interrupted; a child blocked in `open` on a writer-less FIFO holds those pipes open, so `output()` never returns. §5 step 8 says only "Same three rules from step 7 apply", and rule 3's recipe (`recv_timeout` on a thread) leaves the child alive after the timeout — the test fails, and a wedged `trim_galore` is stranded on the runner. The missing step is: `spawn()` rather than `output()`, keep the `Child`, and on timeout `kill()` then `wait()`.

The same applies to the FIFO **writer** child: give it `Stdio::null()` and kill it in cleanup. Under nextest a leaked child holding output handles trips `leak-timeout` (default 100 ms) and the test is reported LEAK; under the `cargo test --release` leg it can stall output collection.

### 3.3 §5 step 7's rule set contradicts §5 step 5

Rule 2: "Keep the payload well under 64 KiB so the writer completes into the pipe buffer and **exits** rather than lingering."
Step 5: "a FIFO is rejected with the FIFO message, **in both writer states** — writer exited, and **writer still open**."

Under rule 2, whether the writer is still open when we fstat is a scheduler race — which is precisely the property MEASUREMENTS §1 identifies as making E vs F "a scheduler coin-flip". So as specified, the writer-still-open case is nondeterministic, and it is nondeterministic in the direction that produces the *other* test's conditions rather than failing loudly.

The deterministic construction, and it needs no child process: hold the write end **in the test process on a spawned thread** —

```rust
let w = { let p = fifo.clone(); thread::spawn(move || {
    let f = std::fs::OpenOptions::new().write(true).open(&p).unwrap();
    std::thread::park();            // hold the writer open
    drop(f);
}) };
```

— because the write-only `open` rendezvous with our reader's `open`, so both complete. It must not be on the main thread: a write-only `open` on a FIFO blocks until a reader appears, so doing it before `detect_input_format` on the same thread self-deadlocks. That trap is not in the rule set, and it is the obvious way to write the test.

### 3.4 V7's fourth cell is not implementable as a unit test

> Byte probe: `reopen_restarts` on a regular file → `true`; on `/dev/stdin` from a pipe → `false`

A unit test cannot control its own process's fd 0. Under `cargo test --release` run from a developer's terminal, `File::open("/dev/stdin")` + `read` **blocks waiting for typed input**; under nextest fd 0 is whatever the harness supplies. So this row is either a hang or a coin-flip depending on how the suite was invoked.

Two portable substitutes, both measured here:

1. `reopen_restarts(&regular_file, b"ZZZZ")` → `Ok(false)`. Pure comparison logic, no I/O tricks, no platform dependence.
2. `/dev/urandom` → `Ok(false)` on both platforms; present on both runners; cannot block; and it is the §1.8 case the edge-case table forgot.

Do **not** use `/dev/fd/{fd}` of a file you opened (the shape I used to reproduce row B): it is a `dup` on Darwin but a fresh open on Linux, so it yields `false` on macOS and `true` on ubuntu. The byte probe's negative direction has no portable *real-world* reproduction, which is worth saying in §9 rather than leaving V7 to be discovered at implementation time.

### 3.5 One more missing row

§9 tests the guard on inputs, not the property that broke in §1.2 — *that every input-bearing CLI surface reaches the guard*. Add a row that enumerates them (`cli.input`, `cli.passthrough`) rather than exercising one of them, so the next flag that takes a path is caught by the test rather than by a reviewer.

Rows that are sound as written: V2's sentinel-appended negative control is exactly the right instinct; V3 is the right invariant to name even though the guard adds no transformation; V4's string is now correct against `format.rs:216`; V8 is adequate. I found no docs debt to add to §5 — `grep` over `docs/`, `README.md` and `--help` finds no `/dev/stdin`, no `mkfifo`, and no `… | trim_galore`, so the CHANGELOG entry really is the whole documentation surface.

---

## 4. Efficiency

§6 is sound and the declined optimisation is declined for the right reason. Two small corrections:

- On the **gzip** branch a `detect_input_format` call now performs **three** opens (`:208`, the probe, `:231`), not two. Immaterial against a full-file read, but the declined-optimisation paragraph reads as if the count were two.
- The offset variant in Alt-1 is strictly cheaper still: no second `read` at all.

No allocation, no per-record cost, `fstat` on an already-open handle is free. A4 holds.

---

## 5. Alternatives

### Alt-1 — discriminate on the second handle's **file offset**, not on bytes (recommended)

After the second `open`, `Seek::stream_position()` (std, no `libc`, no `read`) is `0` for a genuine fresh open and non-zero when the open shared an offset — which is exactly what rows B, D and F are. Measured on this machine:

| input | second open `stream_position()` | offset verdict | byte verdict |
|---|---|---|---|
| regular 4 B / 6 B / 40 B | `Ok(0)` | restartable | restartable |
| `/dev/fd/N` over 4 B | `Ok(4)` | **not restartable** | **accepted by a prefix compare** (§1.4) |
| `/dev/fd/N` over 6 B | `Ok(4)` | not restartable | not restartable |
| `/dev/fd/N` over 40 B | `Ok(4)` | not restartable | not restartable |
| `/dev/zero`, `/dev/null` | `Ok(0)` | restartable | restartable |
| `/dev/urandom` | `Ok(0)` | restartable | not restartable |

Why it is better: it is **exact rather than heuristic**, so it retires A2 entirely (no pathological repeating-prefix residual); it is length-insensitive, so §1.4 and §1.5 both disappear; and it needs no read on the second handle. Why it does not replace the type check: the second `open` on a FIFO still blocks, and `seek` on a FIFO returns `ESPIPE`. So the two-layer structure and its ordering survive unchanged — this only swaps the second layer's predicate.

Note this is *not* the seekability probe §2.2 rightly rejects. That one asks "is this fd seekable" of the **first** handle and accepts row B. This asks "did the **second** open start at zero", which is the property that actually differs.

Cheapest form if you would rather not re-argue the design: keep the byte compare and add the offset as an extra reject condition — *not restartable if `pos != 0` **or** the bytes differ*. One line, closes §1.4's hole exactly, and preserves the `/dev/urandom`-class behaviour you already have.

### Alt-2 — type check only, drop the byte probe (correctly rejected)

Row B fstats as regular, so it would be accepted and then fail with the original misleading message. §2.2 already rebuts this. Keep the rebuttal in the plan: it is the first thing a reviewer of the diff will propose, and the rebuttal is not obvious.

### Alt-3 — compare `(st_dev, st_ino)` across the two handles (does not work)

On Darwin `open("/dev/fd/N")` is a `dup`, so the second handle reports the same device and inode as the first. The offset is the only thing that differs — i.e. Alt-1. Worth one line in §10 so nobody re-invents it.

---

## 6. Action items

### Critical

**C1. Close the `--passthrough` bypass, and re-do the §2.3 inventory.**
`src/main.rs:768` uses `FastqReader::sanity_check`, not `sanity_check_any`, so `cli.passthrough` never reaches `detect_input_format`; its readers are built at `src/main.rs:1344` and `src/main.rs:1383`. §2.3's "the guard fires first in every flow. **Verified** by call-site grep" is false, and 1344 is not even listed. Route the passthrough path through the guard (`sanity_check_any(pt_path)?`, or an explicit `detect_input_format(pt_path)?` beside it), and add the §3.5 validation row. While re-doing the inventory, mark each entry as production / test-only / self-produced-file — `specialty.rs:531` is test code and `demux.rs:170` reads our own output (§1.3).

**C2. Fix the FIFO test rule set: rule 2 contradicts the writer-still-open case, and the deterministic construction is not stated.**
Specify the in-process writer thread (§3.3), and state explicitly that the write-only `open` must not happen on the test's main thread. Add the fourth rule the set is missing: after a timeout, `kill()` + `wait()` the child; `Stdio::null()` for the writer's stdio. Replace "reap it, do not `wait` on it while holding the reader" with the action to take.

**C3. Replace V7's `/dev/stdin`-from-a-pipe cell with something a unit test can actually do.**
It is unobservable from inside the harness and blocks on a terminal (§3.4). Use `reopen_restarts(&regular_file, b"ZZZZ") == Ok(false)` plus `/dev/urandom`. Record in §9 that `/dev/fd/N` is *not* a portable substitute (dup on Darwin, fresh open on Linux) and that the byte probe's negative direction has no portable real-world reproduction.

### Important

**I1. Pin the comparison as length-aware.** §4 and §3.2 leave it ambiguous, and the natural prefix reading silently accepts every `/dev/fd/N` input of ≤4 bytes (§1.4, measured). Spec it as `n2 == first.len() && buf[..n2] == first`.

**I2. Loop the second read** until `first.len()` bytes or EOF, retrying `Interrupted`, so a short read on a restartable file is not reported as non-restartable (§1.5). This is what makes A1 exactly true.

**I3. Add mechanical bounds to CI**: `timeout-minutes` on `rust-tests` and a nextest `terminate-after`. Today nothing bounds a hang, and a dead type check manifests as a hang rather than a failure (§3.1).

**I4. Give §5 step 8 its own recipe.** `Command::output()` on a blocked child never returns; "same three rules" does not cover a subprocess (§3.2).

**I5. Reword the FIFO message.** "named pipe (FIFO)" is false for two of the three #379 reproductions, and no assertion catches it because V1 targets the shared tail (§1.6). "is a pipe or FIFO, not a regular file".

**I6. Include `is_socket()` in the type check.** One `||`. It removes the plan's reliance on an unmeasured claim (§2.2.1 cites MEASUREMENTS §2, which has no socket row), and it converts a possible raw-errno path — Linux `open` of a socket through `/proc/self/fd/N` returns `ENXIO`, which §10 Q3 turns into a bare OS error — into the intended diagnostic.

**I7. Fix §3.2's character-device row and add the blocking-char-device residual.** `/dev/urandom` gets the *new* message, not an existing one (§1.8, measured), and `trim_galore /dev/stdin` at a terminal still blocks in the first `read` with no message (§1.7) — a residual of the same class as FIFO-with-no-writer, involving an input form §3.1's own remedy text names.

### Optional

**O1.** §5 step 6: dev-dependencies are `serde_json`, `tempfile`, **`bstr`** (`Cargo.toml:50-56`), not "tempfile only". The argument (no `libc`, no `nix`) is unaffected. `mkfifo` confirmed at `/usr/bin/mkfifo` on Darwin.

**O2.** §5 step 4: the comment is `format.rs:227-230`; `:226` is the `if`.

**O3.** State "one `fresh_tmpdir` slug per FIFO test" explicitly. It is load-bearing rather than hygienic here: the debug leg is process-per-test under nextest, and `fresh_tmpdir` opens with `remove_dir_all`, so a shared slug lets one process unlink another's FIFO while a writer is blocked on it, stranding the writer. ci.yml:52-56 records a completed audit that all 19 existing `temp_dir().join(...)` sites use unique slugs — worth not being the first exception.

**O4.** §6: say three opens on the gzip branch (§4).

**O5.** Add Alt-3 to §10 so the dev/ino idea is pre-rebutted.

**O6.** `MEASUREMENTS_reopen_probe.md` header says "Darwin 26.5.1"; 26.5.1 is the macOS product version and the kernel is Darwin 25.5.0. Trivial, but a measurement doc should name its platform unambiguously so a future reader can tell which runner it does or does not speak for. The larger point stands and the plan already concedes it: every row in that file is Darwin, and the plan's Linux claims (A6, row B "may genuinely restart") are reasoned from POSIX rather than measured. That is acceptable here only because the ubuntu integration test verifies the one that matters.

---

## 7. What v2 got right, and what v2 introduced

**Right:** the fstat pre-filter is the correct fix for the v1 hang, it needs no dependency, it costs nothing, and §2.2.1's explanation of why it *guards* rather than *replaces* the byte probe (row B fstats regular) is exactly the argument that makes the two-layer design necessary rather than belt-and-braces. §11's account of why v1's dismissal of the type check was aimed at the wrong target is honest and correct. The Q3 change from `Ok(false)` to `Err` is right for the reason given — reporting `EMFILE` as "your input is a pipe" would be the same wrong-blame one layer down — and it is now safe precisely because the type check handles the cases v1 needed `Ok(false)` for.

**Introduced by v2:** nothing in the design. Two in the specification. The Q3 change to `Err` widened the surface on which an unexpected `io::Error` becomes a bare OS message, which is what makes I6 worth doing (`ENXIO` on a socket now lands there). And v2's §5 step 7, written to make FIFO tests deterministic, added rule 2 without reconciling it against step 5's writer-still-open case — so the rule introduced to remove a race specifies one (§3.3).

Everything the plan asserts about its own line numbers and error strings checks out this time: `format.rs:208/211/215-217/219/227-230/231/256/273`, `fastq.rs:474/478/485`, `adapter.rs:115/268` with `MAX_SCAN_READS = 1_000_000`, `main.rs:30/51/68/178/187/1328/1900/2050/2051`, and the quoted `format.rs` comment. V4's empty-file string is correct.
