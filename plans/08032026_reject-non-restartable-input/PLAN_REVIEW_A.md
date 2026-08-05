# PLAN_REVIEW_A — Reject non-restartable input with an accurate diagnostic

**Reviewer:** A (independent, fresh context)
**Target:** `plans/08032026_reject-non-restartable-input/PLAN.md` revision **v2**
**Supporting evidence reviewed:** `MEASUREMENTS_reopen_probe.md`
**Date:** 2026-08-04 · Darwin 26.5.1 (arm64) · rustc 1.95.0

**Verdict: the design is sound and the v2 fix is the right one. Two Critical items must
be resolved before implementation — neither invalidates the two-layer guard.** One is a
whole input path that bypasses the guard (`--passthrough`), contradicting a claim §2.3
marks "Verified". The other is that the FIFO test specification cannot construct one of
the two states it requires, and the wedge it produces lands in test *setup*, which §5
step 7's bound does not cover — with no CI backstop, that is a 6-hour job timeout on
four job-steps.

---

## 0. What I verified independently

Every line reference and quoted string in the plan checks out except where noted below.
I reproduced the plan's measurements rather than accepting them, and added the cases §3.2
omits.

| Plan claim | Result |
|---|---|
| `format.rs:208` open, `:211` read, `:215-217` `n==0`, `:219` `@`, `:226-230` comment, `:231` second open | ✅ all exact |
| `fastq.rs:474` `FastqReader::open` in `sanity_check`, `:478` `seems to be completely empty` | ✅ exact (V4's string is right this time) |
| `fastq.rs:485` `doesn't seem to be in FastQ format (first line doesn't start with '@')` | ✅ exact — the only occurrence in `src/` |
| `format.rs:216` `Input file '{}' is empty` | ✅ exact |
| §2.3 `detect_input_format` inventory (17 sites) | ✅ complete — my grep found no site the plan missed |
| `adapter.rs:115`/`:268` via `open_sync_reader`; `MAX_SCAN_READS = 1_000_000` (`adapter.rs:84`) | ✅ |
| A1: nothing in `tests/`, `.github/`, `justfile` uses `/dev/stdin`, `mkfifo`, or `<(…)` | ✅ re-verified with positive control (the `<(` hits are Rust tuple syntax) |
| A3 / MEASUREMENTS row E: second open on a writer-exited FIFO blocks | ✅ **reproduced** — blocked >3 s on this box |
| CI matrix `[ubuntu-latest, macos-latest]` (`ci.yml:28`) | ✅ |
| §5 step 6 "dev-dependencies are `tempfile` only" | ❌ they are `serde_json`, `tempfile`, `bstr`; and `libc 0.2.184` is **already in `Cargo.lock`** (see O1) |
| §2.3 "`specialty.rs:531`" as a production call site | ❌ it is inside `#[cfg(test)] mod tests` (a `read_fastq` helper) |
| §2.3 / §7 "every input path passes through `detect_input_format`" | ❌ **`--passthrough` does not** (see C1) |

**I ran the plan's §3 guard verbatim** (open → fstat → read → `n==0` → reopen-and-compare)
against every case in §3.2 plus the ones it omits. Darwin:

```
reg.fq            FALLS-THROUGH peek=[64,114,49,10]   type=[reg]
tiny.fq (1 byte)  FALLS-THROUGH peek=[64]             type=[reg]      ← §3.2 row confirmed
/dev/stdin <pipe  step2 FIFO-MESSAGE                  type=[fifo]     ← caught by pre-filter
/dev/stdin <file  step5 REOPEN-MESSAGE                type=[reg]      ← byte probe is load-bearing
<(cat reg.fq)     step2 FIFO-MESSAGE  (/dev/fd/63)    type=[fifo]
symlink→fifo      step2 FIFO-MESSAGE                  type=[fifo]     ← follows the link, fstats the target
symlink→reg       FALLS-THROUGH                                       ← positive control
/dev/null         step4 EMPTY-MESSAGE                 type=[chardev]  ← as §2.2.1 claims
/dev/zero         FALLS-THROUGH peek=[0,0,0,0]        type=[chardev]  ← as §2.2.1 claims
/dev/urandom      step5 REOPEN-MESSAGE                type=[chardev]  ← not in the plan's list
adir/ (directory) step3 READ-FAILED "Is a directory"  type=[dir]      ← unchanged from today
/dev/tty          step1 OPEN-FAILED ENXIO             type=[chardev]
```

So: the two-layer guard behaves exactly as §3 specifies on Darwin, including the
awkward row B, and it does not disturb the sub-4-byte, `/dev/null`, `/dev/zero`,
directory or symlink cases. **The core design is correct.** My findings are about
coverage, tests, and validation — not about the mechanism.

Probe sources are in `$TMPDIR` (`rev_a_probe3.rs`, `rev_a_probe5.rs`, `rev_a_probe6.rs`,
`rev_a_probe7.rs`); each bounds every wait with `mpsc::recv_timeout(3s)` internally, per
the reviewer brief.

---

## 1. Logic review

### C1 — `--passthrough` is an input path that never reaches `detect_input_format`

§2.3 states, and marks **Verified**:

> `format::detect_input_format` is the function every input path passes through […] so
> the guard fires first in every flow.

and §7 repeats it: *"inside `detect_input_format`, which every input path already calls
before any reader is constructed"*.

That is false for `--passthrough`. The passthrough file lives in `cli.passthrough`, not
`cli.input`, so it is not in the `main.rs:187` all-inputs map. Its full lifecycle is:

| Step | Site | What it does |
|---|---|---|
| 1 | `cli.rs:884` | `if !pt.exists()` — a `stat`; a FIFO passes |
| 2 | `main.rs:768` | `FastqReader::sanity_check(pt_path)` — **open #1**, one record, drop |
| 3 | `main.rs:1344` (parallel) / `main.rs:1383` (serial) | `FastqReader::open{_threaded}(p)` — **open #2** |

Two opens, no `detect_input_format`, no guard. Consequences after this change ships:

- **`--passthrough <writer-exited FIFO>` hangs with no message** at step 3 — and by then
  `ensure_output_dir` has run and the R1/R2 writers exist, so it hangs *with partial
  output on disk*. This is the exact failure mode §1 says the plan removes.
- **`--passthrough <live-writer FIFO or /dev/fd/N>`** desyncs: step 2 consumes a record,
  step 3 resumes mid-stream, and the three-way header sync check (`parallel.rs:640`)
  fires with a message about a truncated or desynced passthrough — blaming the user's
  data for a limitation of ours, one flag over.

The plan's own words for this class: *"a diagnostic that blames the user's data for a
limitation of ours"*.

This is not a reason to move the guard, but the claim must be corrected and the path
covered. Cheapest in-place fix: route the passthrough through the same entry point —
replace `main.rs:768`'s `FastqReader::sanity_check(pt_path)` with `sanity_check_any(pt_path)`
(it already dispatches through `detect_input_format` at `main.rs:30`; passthrough + uBAM is
rejected at `main.rs:253`, so the BAM arm is unreachable and harmless). Better still, see
**ALT-1**, which covers this *and* the §11 Medium residual in one place.

**Evidence that would settle it:** already gathered — `grep -n "passthrough" src/main.rs`
shows `:768`, `:1344`, `:1383` and no `detect_input_format` between them.

### C2 — the FIFO test specification cannot build one of its two states, and its bound does not cover the wedge

§5 step 5 requires:

> a FIFO is rejected with the FIFO message, **in both writer states** — writer exited, and
> writer still open. With the type check these are deterministic; without it the first one
> hangs

and §5 step 7 rule 3 is the backstop:

> **Bound every FIFO test in-process** — run **the call** on a thread and `mpsc::recv_timeout`

I implemented all three readings of "writer exited" and measured them. Results:

```
Shape A  writer held open in-process for the whole call  -> FIFO-MESSAGE in 144 µs   deterministic
Shape B  child.wait() for the writer BEFORE opening      -> *** BLOCKED >3s ***      in SETUP
Shape C  spawn writer, then call immediately (rule 1)    -> FIFO-MESSAGE in 28 ms    deterministic
Shape D  v1's second open, writer gone                   -> *** BLOCKED >3s ***      (row E reproduced)
```

Three separate problems:

1. **"Writer exited" is unreachable as a pre-condition.** A writer's `open(O_WRONLY)` on a
   FIFO blocks until a reader opens, and the only reader is `detect_input_format`'s own
   `File::open`. So the writer cannot have exited before the call under test begins —
   Shape B blocks in `child.wait()`. What MEASUREMENTS row E actually observed is
   "writer exited between the *first* and the *second* open", which v2 deletes by
   construction: with the fstat pre-filter there is no second open, so Shapes B and C
   traverse identical code and the "two states" collapse into one test.

2. **The wedge is in setup, and rule 3 bounds only "the call".** An implementer who reads
   "writer exited" literally writes `child.wait()` (or an unbounded `OpenOptions::write`)
   in setup, outside the `recv_timeout`. Shape B is that code, and it blocks. Rule 3 as
   worded does not catch it. The rule must be *bound the whole test body*, not the call.

3. **There is no CI backstop at all**, so an escaped wedge is maximally expensive. I
   checked: no `.config/nextest.toml` exists (so nextest's `slow-timeout` is warn-only,
   `terminate-after` unset), and `grep -rn "timeout-minutes" .github/workflows/` returns
   **nothing**. The `rust-tests` job runs `cargo nextest run` *and* `cargo test --release`
   on both matrix legs (`ci.yml:66`, `:74`), so one wedged FIFO test burns the GitHub
   default **360-minute** job timeout on **four** job-steps. §9's "a hanging test looks
   like a slow job in CI" understates this by two orders of magnitude.

A fourth rule is also missing, and I hit it: **the writer must tolerate `EPIPE`.** In
Shape A the guard bails and drops the reader before the writer's `write_all` lands, so my
writer thread panicked with `BrokenPipe`. Joining that thread fails the test; not joining
it hides a real setup failure. Symmetrically, a `sh -c 'cat … > fifo'` writer is killed by
SIGPIPE, so `assert!(child.wait()?.success())` is payload-size-dependent — flaky. (Rule 2's
"keep the payload well under 64 KiB" turns out to still be necessary, but for *this*
reason, not the stale second-open reason it gives.)

And a fifth, because `mkfifo` is shelled out to: **if `mkfifo` fails, `cat x > f` silently
creates a regular file**, `detect_input_format` returns `FastqPlain`, and the test fails
with "expected the FIFO message" — a diagnosis pointing at the guard when the fault is in
the harness. Assert the FIFO exists and *is* a FIFO in setup.

---

## 2. Assumptions

**A1 (every currently-working input is restartable) — holds.** Re-verified independently
with a positive control. Note the reason it holds is stronger than the plan says: on
Linux, row B (`/dev/stdin` from a regular file) genuinely restarts, so
`trim_galore /dev/stdin < reads.fq` works on `ubuntu-latest` today and continues to,
because the byte probe returns `true` there. On Darwin it fails today and fails after.
No regression either way. ✅

**A2 (≤4 bytes is enough) — holds, but is missing a sibling.** The plan reasons about
*content* (a stream repeating its first 4 bytes) and never about *length*. See I1.

**A3 (`File::open` on a writerless FIFO blocks) — verified, twice over.** I reproduced
row E (Shape D) and also confirmed the writerless case directly: `File::open` blocked
>3 s while `Path::exists()` and `fs::metadata()` both returned promptly. ✅

**A6 (`fstat` names a FIFO on both platforms) — overstated.** The POSIX `S_IFIFO`
guarantee the plan invokes covers *named* FIFOs, which is rows E/F. Rows C and D are a
different mechanism, and the plan itself says so in §2.2.1 (Darwin `dup`s `/dev/fd/N`;
Linux resolves through `/proc/self/fd/N`). So "`is_fifo()` fires for `/dev/stdin` from a
pipe **on Linux**" is an inference, not a measurement, and the whole measurement table is
Darwin-only. See I5 — the inference is almost certainly right, and it is cheap to make CI
prove it rather than assume it.

**A7 (Unix-only) — correct and adequately flagged.** Worth noting the crate ships a
`src/lib.rs`, so a bare `use std::os::unix::fs::FileTypeExt;` is the first hard *compile*
break for a Windows consumer. Three lines of `#[cfg(unix)]` now is cheaper than finding
out later (O5).

**Missing — A8: both `read` calls return the same number of bytes.** The probe's
correctness depends on it and nothing states it. See I1.

**Missing — A9: no input path bypasses `detect_input_format`.** §2.3 asserts this as
verified; it is false (C1). It deserves to be an explicit, testable assumption rather
than a parenthetical, because it is the single load-bearing claim behind the placement
decision.

---

## 3. Efficiency

§6 is accurate and I have nothing to add to the numbers: one `fstat` on an already-open
handle, plus one extra `open` and a ≤4-byte `read` per `detect_input_format` call on the
regular-file path; no allocation; a FIFO costs the `fstat` only. Against ~5 opens and a
full-file read per input this is unmeasurable, and the declined gzip-branch optimisation
is correctly declined — splitting the check across two sites and making it
branch-conditional is not worth one `open`.

One observation the plan could make in passing: `detect_input_format` is called
*redundantly* on the same path (`main.rs:178`→`:30`, then `:187`, then `:1328`, then twice
more via the reader factories), so the guard also runs redundantly. That is a small
efficiency wart the plan correctly leaves alone, and incidentally a benefit — the guard
gets several chances to fire. It does mean the same message can be emitted from several
call frames; since the text is identical, no test is affected.

---

## 4. Validation sufficiency (§9)

Could the change ship broken with all nine rows green? **Yes, in four ways.**

### I4 — V1's negative assertion is vacuous on two of its three rows

V1 asserts, for all three #379 reproductions: *"**no** occurrence of `doesn't seem to be
in FastQ format`"*. But that string is only what two of the three produce today:

- `<(cat plain.fastq)` → `detect_input_format` returns `FastqPlain` at `format.rs:220`
  without a second open, then `fastq.rs:485` fires. ✅ the assertion has teeth.
- `cat reads.fq.gz | … /dev/stdin` → the gzip branch re-opens at `format.rs:231` and
  `MultiGzDecoder` chokes on the offset-4 prefix → **`Failed to decompress first block of
  '…' for format detection`** (`format.rs:236`). The V1 negative assertion passes
  vacuously — it would pass on unpatched `dev` too.
- `mkfifo f; cat small.fq > f &` → today this **hangs**; there is no stderr to assert on.

So V1 as written cannot distinguish "fixed" from "failed differently". Fix: per-row
negative controls — assert absence of `doesn't seem to be in FastQ format` on the plain
row, absence of `Failed to decompress first block` on the gz row, and for the FIFO row
assert *promptness* (which V6 does, informally). §1's summary sentence should likewise
name the gz message; right now it generalises one row's diagnostic to all three.

### I3 — V7's byte-probe negative direction is not constructible

V7 requires: *"Byte probe: `reopen_restarts` on a regular file → `true`; on `/dev/stdin`
from a pipe → `false`"*. The second half cannot be a unit test:

- `/dev/stdin` from a pipe is a FIFO, so in `detect_input_format` the type check bails
  first and `reopen_restarts` is never called on it;
- calling `reopen_restarts("/dev/stdin", …)` directly from a test means depending on
  whatever `cargo test` / `cargo nextest` hands the test process as fd 0 — inherited under
  one runner, `/dev/null` under another, and nextest's per-test process isolation differs
  again. Non-deterministic and platform-divergent, i.e. exactly what §5 step 8 rightly
  refuses to do for row B.

This matters because §9 itself names V7 as one of the three rows that "would quietly not
hold". A probe only ever observed returning `true` is not a probe — and on
`ubuntu-latest` the byte probe will in practice *never* return `false`, since row B
restarts there and every pipe is caught by the type check. Linux CI would exercise the
`false` branch **zero times**.

Deterministic, portable replacement that keeps the teeth: call the helper with a
deliberately wrong `first`.

```rust
// negative control: the comparison itself must be able to say "no"
assert_eq!(reopen_restarts(&regular_fq, b"XXXX")?, false);
assert_eq!(reopen_restarts(&regular_fq, b"@rea")?, true);   // positive, same file
```

Pair it with a Darwin-only `#[cfg(target_os = "macos")]` row-B case if you want the real
shared-offset mechanism covered — but the mismatched-`first` test is the one that must run
on both legs.

### I5 — the Linux leg is asserted, not measured

§5 step 8 commits the integration test to asserting `cannot be re-read from the start` for
`/dev/stdin` from a pipe on **both** legs, and V9 asserts both legs green. All six
measurement rows are Darwin. I could not run a Linux probe here (the sandbox denies the
Docker socket), so I reasoned it instead:

- `/proc/self/fd/N` for an anonymous pipe re-opens the pipefs inode; `fstat` gives
  `S_IFIFO`, so `is_fifo()` fires and the guard produces the FIFO message. Consistent with
  Darwin's `dup` route reaching the same verdict by a different mechanism.
- The re-open does **not** block even with the writer closed: the kernel's `fifo_open`
  blocking wait for a partner is gated on `!is_pipe`, and `is_pipe` is true for a pipefs
  inode. Named FIFOs (`is_pipe == false`) do block — matching row E on both platforms.

I believe the plan is right. But "I believe" is the wrong footing for a CI assertion that
gates a merge, and the fix costs nothing: **have the test assert the mechanism instead of
trusting it.** In the integration test, before invoking the binary, assert
`fs::metadata("/proc/self/fd/0" or "/dev/stdin")?.file_type().is_fifo()` — or simply have
the unit test assert `is_fifo()` on a *named* FIFO (POSIX-guaranteed, both legs) and let
the integration test assert only the shared substring, which §5 step 8 already sensibly
does. Then the first CI run *is* the Linux measurement, and a platform surprise reports
itself as a specific failed assertion rather than as an unexplained message mismatch.

### Also missing from §9

- **No row for `--passthrough`** (C1) — the gap is invisible to all nine rows.
- **No negative control on the FIFO test harness** (C2, point 5): nothing proves the
  fixture is a FIFO rather than a regular file `>`-created by a failed `mkfifo`.
- **No row asserting the guard's cost is not paid on the hot path.** A4 is asserted from
  the call-site inventory, which is fine — but if anyone later moves the check into
  `FastqReader::open`, nothing fails. Optional: assert `detect_input_format` is not
  reachable from `next_record` by inspection only; not worth a test.

V2, V3, V4, V6 and V8 are well-specified. **V2's "prove the comparison can fail by
diffing against a sentinel-appended copy" and V4's corrected string are both exactly
right** — that discipline is what caught v1's error, and it is applied here.

---

## 5. Alternatives

### ALT-1 (recommended, additive): a non-blocking **path** stat in `Cli::validate`

The measurement file already contains the evidence for this and the plan does not use it:
`MEASUREMENTS_reopen_probe.md` §2's second column is `fs::metadata(path)`, and it names
FIFO/pipe for rows C, D and E — identically to the fstat-on-handle column.

I confirmed the decisive property the plan never tests: **a path stat does not block on a
FIFO with no writer.**

```
FIFO, no writer ever:   Path::exists = true
                        fs::metadata = Ok((fifo=true, file=false, chardev=false, sock=false))
                        File::open   = *** BLOCKED >3s ***
regular file (control): fs::metadata = Ok((fifo=false, file=true, …));  File::open returned
```

`cli.rs:933` already stats every input (`if !path.exists()`), unconditionally, for every
mode, at `main.rs:166` — before any `File::open`. Replacing that `exists()` with a
`fs::metadata()` match costs **zero extra syscalls** and buys three things:

1. **It kills §11's Medium residual.** `trim_galore fifo` with no writer ever currently
   blocks at the first open with no message, and §11 calls this "unfixable without opening
   every input `O_NONBLOCK`". That is not so — `stat(2)` never opens, so the FIFO is named
   and rejected before anything can block. No `libc`, no `unsafe`, no change to the open
   path for any input.
2. **It covers `--passthrough`** (C1) if applied at `cli.rs:884` too, and the `--demux`
   barcode file at `cli.rs:926` for free.
3. It fires earliest, so nothing is created before the message.

§2.3's stated objection — *"A guard in CLI validation would miss `clump_only` and
`specialty`, which call `detect_input_format` directly"* — does not hold for the binary:
`main()` calls `cli.validate()` at `main.rs:166`, before every dispatch branch. It holds
only for a third-party library consumer calling `clump_only::clump_only_single` without a
`Cli`, which is a real but narrow case.

So the right answer is **both layers, not either**: keep the handle-level `is_fifo()` in
`detect_input_format` (defence in depth, covers library callers, TOCTOU-proof as A6
notes) *and* add the path-level check in `Cli::validate`. The TOCTOU concern A6 raises
against a path stat is benign for a *rejection* guard in a layered design: if the path
changed type between validate and open, the handle-level check is the backstop, and vice
versa. Two cheap checks that fail closed.

Cost: ~10 lines and one more test. Payoff: one Critical finding and one Medium residual
risk, both closed.

### ALT-2 (considered, don't take): whitelist regular + char devices instead of blacklisting FIFOs

`if !(ft.is_file() || ft.is_char_device()) { bail }` would also give sockets and
directories an accurate message instead of today's `ENXIO` / `Is a directory (os error 21)`.
I verified both of those outcomes are unchanged by the plan. But a whitelist is a
behaviour change for input classes nobody has complained about (block devices), and it
converts "unknown file type" from *fall through and let the format check decide* into a
hard error. The plan's blacklist is the conservative choice and the right one for a
diagnostic-only fix. Worth one line in §10 as a considered-and-declined alternative,
since a reviewer will ask.

### ALT-3 (declined, but the cost was mis-stated): `O_NONBLOCK`

§11 rejects this as needing `libc` and changing the open path for all inputs. The second
half is a good objection. The first half is not: `libc 0.2.184` is **already in
`Cargo.lock`** (transitively, via `rustix`/`tempfile`). ALT-1 achieves the same outcome
with neither cost, so `O_NONBLOCK` stays declined — on the better reason.

### ALT-4 (out of scope, worth a sentence): memoise detection per path

Caching `detect_input_format`'s result would remove the redundant opens §2.1 documents
*and* make the probe genuinely free. It collides with #374 and is a bigger change than
this plan wants. Noting it because it is the change that would make §6's declined
optimisation moot.

---

## 6. Action items

### Critical

**C1. Correct §2.3/§7 and cover `--passthrough`.**
The claim that every input path passes through `detect_input_format` is false; `--passthrough`
is opened twice (`main.rs:768`, then `main.rs:1344`/`:1383`) with no guard, so
`--passthrough <fifo>` still hangs with no message *after* output files exist. Either
change `main.rs:768` to `sanity_check_any(pt_path)`, or take ALT-1 (which covers it at
`cli.rs:884` and closes the §11 residual too). Add a §9 row. Downgrade "Verified by
call-site grep" to what the grep actually verified — `cli.input`, not every input.

**C2. Rewrite §5 steps 5 and 7 for the FIFO tests, and add a CI backstop.**
- Drop "in both writer states". It is unreachable as a pre-condition (measured: Shape B
  blocks in `child.wait()`), and v2's pre-filter makes the two states traverse identical
  code. Specify **one** deterministic test, and say how: hold the write end open on an
  in-process thread for the duration of the call (Shape A, 144 µs), or spawn the writer
  and call immediately (Shape C, 28 ms). Add one sentence recording *why* "writer exited"
  is not a testable state, so the next reader does not re-add it.
- Reword rule 3 from "bound **the call**" to **bound the whole test body, setup
  included** — the measured wedge is in setup, outside the current bound.
- Add rule 4: **the writer must tolerate `EPIPE`.** The guard drops the reader before the
  writer finishes, so an in-process writer must not `unwrap()` its `write_all`, and a
  spawned `cat` must not have `.success()` asserted unconditionally. (Rule 2's payload
  advice survives, but for this reason — update its stale rationale.)
- Add rule 5: **assert the fixture is a FIFO** — check `mkfifo`'s exit status and
  `assert!(fs::metadata(&f)?.file_type().is_fifo())`. If `mkfifo` fails, `cat x > f`
  creates a regular file and the test fails pointing at the guard.
- Add to §5 (new step) and §9: **`timeout-minutes` on the `rust-tests` job** and a
  `.config/nextest.toml` with `slow-timeout = { period = "30s", terminate-after = 2 }`.
  Neither exists today; the job runs `cargo nextest run` and `cargo test --release` on
  both legs, so an escaped wedge costs the GitHub default 360 minutes × 4 job-steps. This
  is the structural backstop that survives an implementer forgetting rules 3–5, and it is
  the one thing here that turns the highest-consequence failure mode into a fast, legible
  failure.

### Important

**I1. Defend against a short `read` on either open, or state the assumption.**
§4/§5 step 2 say "read up to `first.len()` bytes […] return whether the slices are
equal". If open #1 returns 4 bytes and open #2 returns 3 (NFS, FUSE), the comparison is
`false` and a perfectly good regular file is **newly rejected** with "cannot be re-read
from the start" — the plan's own wrong-blame failure, introduced by the fix. Low
probability, but the fix is a 3-line fill loop on both reads, and the same loop hardens
the *existing* `n >= 3` gzip test at `format.rs:226`, which today misclassifies a gzip
file on a short first read. Add it, or add an explicit A8 saying you accept the risk.

**I2. Give `reopen_restarts`'s `Err` path a context line.**
§5 step 2 says "Propagate open/read failure as `Err`" with no `with_context`. v2 changed
Q3 to `Err` precisely to avoid mis-blaming `EMFILE`/`EACCES` — but a bare
`Too many open files (os error 24)` with no path is only half the improvement. Add
`.with_context(|| format!("Failed to re-open input file '{}' for the restartability check", path.display()))`,
matching the style of `format.rs:209`/`:213`.

**I3. Replace V7's byte-probe negative direction.** `/dev/stdin` from a pipe is not
constructible in a unit test and never reaches the probe anyway. Use
`reopen_restarts(&regular, b"XXXX") == Ok(false)` as the portable negative control, plus
the existing positive on the same file. Optionally add a `#[cfg(target_os = "macos")]`
row-B case for the real shared-offset mechanism. Without this, the `false` branch is
exercised zero times on `ubuntu-latest`.

**I4. Give V1 per-row negative controls.** The single assertion "no occurrence of
`doesn't seem to be in FastQ format`" is vacuous on the gz row (which fails with
`Failed to decompress first block of '…' for format detection`, `format.rs:236`) and on
the FIFO row (which produces no stderr today because it hangs). Assert the row-specific
wrong message is absent. Same correction to §1's opening sentence, which generalises the
plain row's diagnostic to all three.

**I5. Stop asserting the Linux mechanism and have CI prove it.** A6's POSIX `S_IFIFO`
guarantee covers named FIFOs only; rows C/D on Linux go through `/proc/self/fd/N`, a
different mechanism, and no row of the measurement table is Linux. My kernel-level reading
says the plan is right (pipefs re-open yields `S_IFIFO` and does not block even with the
writer closed, because `fifo_open`'s partner-wait is gated on `!is_pipe`) — but make the
integration test assert `is_fifo()` on the descriptor it is about to hand over, so the
first CI run *is* the measurement and a surprise names itself.

### Optional

**O1. Fix the dev-dependency facts in §5 step 6.** They are `serde_json`, `tempfile`,
`bstr` — not "tempfile only" — and `libc 0.2.184` is **already in `Cargo.lock`**. So
`libc::mkfifo` adds a line to `[dev-dependencies]` and *zero* new compiled crates, which
makes the stated trade-off ("adding `libc` […] is the heavier option", "no new
supply-chain surface") point the wrong way. `Command::new("mkfifo")` is still defensible
— it needs no `unsafe` — but the honest comparison is "`$PATH` dependency and a failure
mode that masquerades as a guard bug (C2 rule 5)" versus "one dev-dep line on a crate
already in the graph, and an immediate errno". Either choice is fine; state it correctly.

**O2. `specialty.rs:531` is a `#[cfg(test)]` helper**, not a production call site. Drop
it from §2.3's inventory — the remaining entries are all genuine and I verified each.

**O3. Add `/dev/urandom` to §2.2.1's char-device list.** The plan names `/dev/null`
(→ empty) and `/dev/zero` (→ not-FASTQ), both confirmed. `/dev/urandom` reaches the
**new** re-open message, which is defensible but is a third outcome the list implies does
not exist.

**O4. Record ALT-2 (whitelist) in §10** as considered-and-declined. A reviewer will ask
why sockets and directories keep their raw `errno` messages when the plan is about
accurate diagnostics; one line pre-empts it. (The socket claim itself is right in
substance — `open()` on a socket fails — though `ENXIO` is the Linux errno; Darwin
differs. The branch is unreachable either way.)

**O5. Consider `#[cfg(unix)]` now rather than flagging it.** A7 is honest, but the crate
ships a `src/lib.rs`, and a bare `use std::os::unix::fs::FileTypeExt;` is the first hard
compile break for any non-Unix target. Three lines (`#[cfg(unix)]` on the check plus a
`#[cfg(not(unix))]` no-op) keeps the door open at no cost to the Unix path.

---

## 7. What v2 got right, and did v2 introduce anything new

The v2 fix is correct and I could not break it. The fstat pre-filter is the right shape,
in the right order, and it does exactly what §2.2.1 claims — I reproduced the whole §3
flow and every row lands where the plan says, including row B needing the byte probe and
row E never reaching it. Q3's change from `Ok(false)` to `Err` is a genuine improvement
and correctly reasoned. V4's corrected string is right. The scope discipline (A over B/C,
disjoint from #374, with a note for whoever attempts B) is exemplary.

**Did v2 introduce anything new?** Two things, both small and both listed above:

- **I1** is new to v2 in *consequence*, not in mechanism: v1 also compared bytes, but v2
  made that comparison the sole gate for non-FIFO input, so a short read now turns a
  working regular file into a hard rejection. Nothing in the plan mentions read length.
- **The stale rationale for rule 2.** v2 removed the second open that rule 2 existed to
  survive, but kept the rule with its old justification. It is still needed (EPIPE /
  writer exit status), so the rule is right and its reason is wrong — the kind of drift
  that gets a rule deleted by the next reader.

Neither is a design defect. C1 and C2 are the two items I would hold implementation on:
C1 because the plan's central placement argument is provably incomplete and an implementer
will not re-check a claim marked "Verified", and C2 because an unbounded wedge in test
setup with no CI timeout is, as the plan itself says, the failure that looks like a slow
job — for six hours, four times over.
