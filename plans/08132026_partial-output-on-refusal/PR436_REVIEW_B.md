# Code review B — PR #436 (`fix(io): publish an output only if its teardown succeeded`, #434)

Reviewer B, fresh context. HEAD `b8aa7d8`, parent `63f5d44` (= `dev`). Built and run in a
sha-pinned private worktree (`$TMPDIR/pr436-b`, `$TMPDIR/pr436-b-base` for the parent), private
`CARGO_TARGET_DIR` for each. No file in the repository was edited; the PR was not touched.

**Verdict: merge with changes.** The code is sound and the three #428 invariants survive. What
is wrong is what the PR *says* about the code: the blast radius of the byte change is understated
in both the CHANGELOG and the PR body, in the one direction that costs a downstream consumer
real time. Two claim corrections and one comment correction; no code change is required for
correctness.

---

## Gates — verified, not assumed

| Gate | Result |
|---|---|
| `cargo test --release` on `b8aa7d8` | 16 suites, **710 passed, 0 failed** |
| `cargo test --release` on `63f5d44` (parent) | 16 suites, **707 passed, 0 failed** |
| `cargo fmt --all -- --check` | clean (exit 0) |
| `cargo clippy --all-targets --release -- -D warnings` | clean, no diagnostics |

**The test-count disagreement resolves; nobody was wrong.** Parent 707 → PR 710 is exactly the
three new tests. The author's 708 is a Linux number and reconciles through platform gating: three
macOS-only tests (`src/format.rs:689`, `src/format.rs:711`,
`tests/integration_non_restartable_input.rs:252`) against one Linux-only test (`src/cli.rs:1950`),
so macOS = Linux + 2. Linux parent 705 + 3 = 708. The "709 Linux baseline" figure in the review
brief was the only wrong number.

## Hot path — measured, because §10 O5 made evidence the precondition

The plan deferred this design partly because "it touches `write_record` on the hot path so it
needs `scripts/benchmark.sh` numbers rather than an assertion". The PR honestly declines to
supply them. I supplied them on different hardware and the precondition clears.

1M reads (100 × `test_files/illumina_10K.fastq.gz`), `--cores 1`, gzip output, alternating runs,
macOS, release builds differing only by this patch:

| | median | min | max | n |
|---|---|---|---|---|
| `63f5d44` | 12.170 s | 10.912 s | 12.263 s | 5 |
| `b8aa7d8` | 11.769 s | 9.771 s | 12.418 s | 5 |

Patched build is nominally **faster** at median and min — the opposite sign from the author's VPS
run, which is the point: the within-build spread swamps the between-build gap on both machines.
Two independent boxes now agree that the vtable→match swap is not measurable. Mechanically that
is expected: the match is inlinable and monomorphic per writer, the indirect call was not.

**No performance objection to this design stands.** That collapses the strongest argument for the
alternative shape below, so I am not raising it as a required change.

---

## Findings

### H1 — Critical claim error: "`--cores N` output is byte-identical either way" is false

`CHANGELOG.md:26`, and the same sentence in the PR body ("`--cores N` output is byte-identical
either way").

The `ParGz` sink changes its bytes too. Measured, base vs PR:

| Command | compressed md5 | decompressed md5 | Δ bytes |
|---|---|---|---|
| `--hardtrim5 30 --cores 2` (ParGz) | **differs** (`e81e727c…` → `a5243d8d…`) | identical | 11 |
| `--demux … --cores 2`, every per-barcode file (ParGz) | **differs** | identical | — |
| SE `--cores 2` trimming | identical | identical | 0 |
| PE `--cores 1` trimming, both mates | differs | identical | 12 each |
| plain (non-gzip) output | identical | identical | 0 |

The mechanism is the same one the PR identified, applied to the other sink: the old teardown
reached `ParCompress::flush` (→ `flush_last(false)`) *twice* — once in `FastqWriter::finish`, once
in the deleted `FastqWriter::drop` — before `ParCompress::drop` ran `finish()` (→
`flush_last(true)`). Each extra `flush_last` inserts a block boundary. The new code calls
`flush_last(true)` once. Output stays a single valid gzip member (`gzip -t` clean, one `1f 8b 08`
occurrence, 10000/10000 records) and decompresses identically.

**The correct generalisation:** bytes change for **every gzipped output written through
`FastqWriter`** — that is serial-gzip *and* `ParGz`, so all four specialty modes and `--demux` at
any `--cores` value, plus SE/PE trimming (which always constructs `FastqWriter` with `cores = 1`;
see `main.rs:1332`, `:1586`). Byte-identity holds only for plain output and for the SE/PE
`--cores N` trimming path, which bypasses `FastqWriter` entirely — its workers gzip into a `Vec`
and the main thread writes raw bytes (`parallel.rs:996`, `:1010`).

Why it matters: the CHANGELOG is what a downstream maintainer reads to decide whether a pinned
compressed hash breaks. nf-core/modules pins compressed md5s in its nf-test snapshots. As written,
the entry tells someone their `--cores N` snapshots are safe when for specialty and demux modes
they are not.

**Fix — `CHANGELOG.md:22-26`, replace the second paragraph:**

```markdown
  **One visible consequence:** the old teardown reached `GzEncoder::flush` / a second
  `ParCompress::flush`, which emit an extra deflate block boundary before the final block; the
  explicit `try_finish` / `finish` do not. Every gzipped output written through `FastqWriter`
  therefore differs from v2.3.0 by a few bytes — measured on a 1M-read output, 101,067,014 bytes
  before and 101,067,003 after, with an identical decompressed md5. Plain output and the SE/PE
  `--cores N` trimming path are byte-identical, so anyone holding a compressed md5 of a
  specialty-mode, `--demux` or serial-gzip output needs to re-baseline it.
```

The PR body's equivalent sentence needs the same correction.

### H2 — "every byte it owes was written" overstates what is delivered: `close(2)` is still unchecked

`CHANGELOG.md:8` (the entry's headline) and `src/fastq.rs:722-724`.

In all three arms the last thing that happens after the trailer is that the `File` drops, and
`File::drop` discards the `close(2)` result. On NFS, CIFS and some quota implementations a
deferred `ENOSPC`/`EIO` is reported *only* at close. So a full disk can still publish a short
output — the window is narrowed from "the whole trailer" to "the close", not closed. There is no
`sync_all` or `sync_data` anywhere in `src/`.

I would not add `sync_all`: an fsync per output file is a real cost on a trimming run and nobody
asked for durability. Soften the claim instead.

**Fix — `CHANGELOG.md:8`:**

```markdown
- **An output file is now published only if its writer closed cleanly, including the gzip
  trailer** ([#434](https://github.com/FelixKrueger/TrimGalore/issues/434)).
```

**Fix — `src/fastq.rs:722-724`,** replace "the final name appears only after every byte the output
owes — including the gzip trailer — has been written successfully" with "the final name appears
only after the sink's trailer has been written successfully (`close(2)` is still unchecked, so a
network filesystem can defer an error past this point)".

### M1 — the `ParGz` arm turns a teardown error into a panic, so `finish`'s doc is wrong for one of three sinks

`src/fastq.rs:611-613`.

```rust
Sink::ParGz(mut w) => {
    w.finish()?;
}
```

`gzp-2.0.2/src/par/compress.rs:377-388`: `ZWriter::finish` runs `self.flush_last(true)?` **before**
it `take()`s `tx_compressor`, `tx_writer` and `handle`. `flush_last`'s only failure is a channel
`send` to a compressor/writer thread that has already died — i.e. precisely the ENOSPC case. So
that failure returns `Err` with all three fields still `Some`; `w` then drops at the end of the
match arm, `ParCompress::drop` (`:396-401`) sees three `Some`s and calls `self.finish().unwrap()`,
whose `flush_last` sends to the same dead receiver and **panics**.

So `finish`'s doc claim — "A failed teardown returns here with `pending` still armed" — does not
hold for `ParGz`: control leaves by unwinding, not by return.

**Data safety survives.** It is a first panic, so unwinding runs `PendingOutput::drop`, `released`
is still `false`, the `.partial` is removed and nothing is published. The symptom is exit 101 with
an `unwrap` panic message instead of a clean error.

Two mitigations against raising this further: it is **pre-existing** (the old `Box<dyn Write>`
dropped the same `ParCompress` through the same `Drop`, so `--cores N` specialty/demux runs
panicked identically before this PR), and the plan's §4 edge-case table already names
`ParCompress::drop`'s `.unwrap()`. Also note the sibling case *does* work: if `flush_last` succeeds
and the writer thread's error arrives via `handle.join()`, all three fields are `None`, `Drop`
no-ops, and the error propagates cleanly. On a real full disk the two outcomes race.

If you want it clean — optional, and it leaks one fd plus one dead thread handle on a fatal path:

```rust
Sink::ParGz(mut w) => {
    if let Err(e) = w.finish() {
        // `finish` leaves its channels armed when it fails early; `Drop` would
        // `unwrap()` the same failure.
        std::mem::forget(w);
        return Err(e.into());
    }
}
```

Otherwise, correct the doc at `src/fastq.rs:725-726` to say the guarantee holds for `Plain` and
`Gz` and that `ParGz` may panic out through `gzp`'s `Drop`.

### M2 — `gz_sink_finish_reports_a_failing_teardown` fails on the gzip *header*, not the trailer its comment names

`src/fastq.rs:749-761`.

Nothing is written to the sink before `finish()`, so `BufWriter`'s buffer is empty, `into_inner()`
flushes nothing, and `GzEncoder::try_finish` → `write_header()` is the **first** write to reach the
read-only `File`. In flate2 1.1.9 `src/gz/write.rs`, `header` is drained only by `write_header`
(`:130-136`), and `try_finish` calls it first (`:89-91`); no `write` preceded it to drain it.

Consequences:

- The comment's central claim — "The bytes that fail here are the ones ONLY `try_finish` writes —
  the deflate end-of-stream plus the 8-byte CRC/ISIZE trailer" — is false. The failing bytes are
  the ten header bytes, which `flush()` also writes (`:148-152`).
- The test therefore does **not** distinguish fixed from broken behaviour. The old
  `finish()` → `self.writer.flush()?` → `BufWriter::flush` → `GzEncoder::flush` → `write_header()`
  would have returned `Err` for this same scenario. The regression #434 describes — record writes
  and header succeed, the *trailer* write fails — is not reached by any of the three tests.
- The same applies to `gz_sink_finish_writes_the_trailer_explicitly`: the old teardown also left
  eight zero bytes at the tail, so the positive control passes either way. It is still a useful
  assertion that the sink emits a complete member; it just is not discriminating.

This is not fixable at the `Sink` level, and that is worth saying out loud: the variants are
concrete over `File`, so there is no seam for a writer that succeeds for N bytes and then fails.
Position-dependent failure injection would need `Sink` to be generic or to gain a test-only
variant, which is more machinery than the fault deserves. The PR body is already explicit that the
failure is not reproduced by a test — the defect is only that this comment claims more than the
test does.

**Fix — replace `src/fastq.rs:750-754` with:**

```rust
        // `try_finish` is the only call that writes the header and trailer, and
        // `Sink::finish` must report its failure rather than discard it (#434).
```

### M3 — the comments are long, and most of their length is historical narrative

A house-rule call for the maintainer rather than a defect, but it is the pattern CLAUDE.md singles
out ("current state only — never reference the state you are fixing", "never justify the change in
the source"):

| Location | Issue |
|---|---|
| `src/fastq.rs:570-586` | 17-line `Sink` doc; a bold thesis sentence, an ENOSPC anecdote, `let _ = self.try_finish()` quoted from code that no longer exists, an issue number |
| `src/fastq.rs:725-726` | "Before this, teardown ran in the sink's own `Drop` and its error was discarded." |
| `src/parallel.rs:948-955` | 8 lines, 6 of them about what the code used to do ("printed 'Worker error' and broke") |
| `src/parallel.rs:310-319`, `:1039-1047` | both worker-join comments spend most of their length on "used to leave a truncated output … exiting 101" |

In fairness: this matches the density of the surrounding code. `io.rs` already carries
multi-paragraph rationale, and the 12-line `drop(result_rx)` comment at `parallel.rs:277-292` that
the single-end path copies is pre-existing. So the PR reads like its neighbours; the neighbours
just do not follow the rule either.

Suggested one-liners if you want them tightened:

- `src/fastq.rs:570-586` → `/// The three concrete sinks a [`FastqWriter`] can write into; named rather than boxed so each one's teardown is callable and its error propagable.`
- `src/fastq.rs:725-726` → delete the last sentence; keep "A failed teardown returns here with `pending` still armed, so its `Drop` removes the temporary and nothing is published."
- `src/parallel.rs:948-955` → `// `Result<_>` in the channel: the main thread decides whether to publish, so it needs the worker's error. Mirrors the paired path.`
- `src/parallel.rs:310-319` and `:1039-1047` → `// Join the workers before publishing: `thread::scope` joins on the way out of this closure, i.e. after the commit.`

### L1 — the byte delta is not a constant

Measured 12 bytes on two different serial-gzip outputs (10K fixture, and the 1M-read input), 7
bytes on `--hardtrim5 30 --cores 1`, 11 bytes on `--hardtrim5 30 --cores 2`; the author measured
11. The CHANGELOG quotes two absolute sizes for one output, which is accurate as written — the note
is only that the figure must not be generalised to "11 bytes" if the wording is ever tightened.
The variability is expected: the removed sync marker is 5 bytes plus 0–7 bits of alignment, and
the final block repacks differently.

### L2, L3, L4 — notes, no user-visible symptom, not worth issues

- **L2** `bail!("Worker thread panicked")` (`parallel.rs:320`, `:1050`) discards the panic payload,
  so the worker's own message is lost. Consistent with the existing "Reader thread panicked".
- **L3** The join loop `bail!`s on the first panicked worker and leaves the rest unjoined, so if
  two workers panic, `thread::scope` re-raises the second on the way out and the run exits 101
  instead of returning the clean `bail!`. Likewise, when the reader errors *and* a worker panicked,
  the early `return Err(reader_err)` skips the worker joins entirely and `thread::scope` panics.
  In both cases the closure's locals — including the `PendingOutput` guards — drop before `scope`
  re-raises, so nothing is published either way.
- **L4** In the single-end path a worker error is now dropped entirely when a reader error also
  exists: `reader_join_err.or(first_error)` discards `first_error`, and the `eprintln!` that used
  to print it is gone. The paired path has always behaved this way, so this is consistency rather
  than regression — but the diagnostic in that one combination is strictly less than before.

Per the issue-filing bar, none of L1–L4 names a symptom worth a GitHub issue. If they ever
accumulate with other small items, L3 and L4 would batch together as "parallel-path error
reporting loses information in multi-fault cases".

---

## What I verified and found correct

**The three #428 invariants hold.**

1. *Commit exactly where the writer is closed, never later.* `FastqWriter::finish`
   (`src/fastq.rs:729-733`) still calls `pending.commit()` in the same function that closes the
   sink, so the four in-run consumers still see the final path. The sharpest of them —
   `clump_only.rs:823-828` and `:1069-1072`, `writer.finish()?` immediately followed by a
   `?`-propagating `std::fs::metadata(&output_path)` — still works, and now fails *earlier* and
   more informatively if teardown fails, since `finish()` returns the teardown error before the
   `metadata` call can produce a confusing `ENOENT`. Confirmed by the clump_only suite passing.
2. *`commit` sets `released` before the rename.* `src/io.rs` is untouched by this PR
   (`git diff --stat` = CHANGELOG.md, src/fastq.rs, src/parallel.rs only).
3. *The abandon path removes the temporary and never publishes,* for all three variants. Field
   order still puts `writer` before `pending`, so drop order is sink-then-guard. `Plain`:
   `BufWriter::drop` flushes (error ignored) — identical in effect to the deleted
   `FastqWriter::drop`'s `let _ = self.writer.flush()`. `Gz`: `BufWriter::drop` flushes into the
   encoder, then `GzEncoder::drop` runs `let _ = try_finish()`. `ParGz`: `ParCompress::drop` runs
   `finish().unwrap()` — same as before, see M1. Then `PendingOutput::drop` sees `released == false`
   and unlinks. **Deleting `impl Drop for FastqWriter` loses nothing**; its only statement was a
   `flush` that `BufWriter::drop` already performs, and removing it is what permits `finish` to
   destructure `self` (a manual `Drop` impl forbids partial moves).

**`Sink`'s partial `Write` forwarding is safe.** The only method that ever reaches the sink is
`write_all`: `FastqWriter::write_record` → `FastqRecord::write_to` (`src/fastq.rs:79-97`), which
issues exactly one `write_all` per record. No call site reaches `write_vectored` or `write_fmt` on
a `FastqWriter` — the only `write!`/`writeln!` in the neighbourhood are `demux.rs:285-292` (a
report `File`) and `fastq.rs:1120` (a test, straight onto a `GzEncoder`). So the two unforwarded
methods cannot reach a default implementation where `Box<dyn Write>` would have dispatched to
`BufWriter`'s override. `Sink` is private, so no external implementor can change that.

**Neither new worker join can deadlock.** In both functions the order is: main loop → unconditional
`drop(result_rx)` → reader join → worker joins. After `drop(result_rx)` a worker parked in
`rtx.send` sees the send fail immediately and exits; after the reader has returned, `txs` is
dropped, so a worker parked in `rx.recv()` sees a closed channel and exits. By the time the worker
joins run, every worker is finished or unblocked on both the clean and the error path. The
single-end path's new `drop(result_rx)` is what makes this true there — without it the early
`break 'outer` would deadlock the reader join, exactly as its comment says.

**`first_error` precedence matches the paired path** — both use `reader_join_err.or(first_error)`,
so a worker error cannot be reported instead of a more informative reader error. The single-end
reader does not send the paired path's lossy sentinel, so its main loop keeps writing already-queued
batches after a reader error; that is harmless, because the reader's error is caught by the join
and returned before `commit`, and the temporary is discarded.

**Output correctness across the byte change.** Every affected output decompresses to an identical
md5, passes `gzip -t`, and holds the same record count. Serial gzip on the 10K fixture and on 1M
reads, PE both mates, `--hardtrim5` at `--cores 1` and `--cores 2`, and all seven `--demux`
per-barcode files.

**No consumer of raw `.gz` bytes in this repo is broken.** The two raw-md5 comparisons in
`.github/workflows/ci.yml` — the alias check at `:685-694` and the `--clump_only` cross-run
determinism check at `:836-837` — both compare one build against itself, so they are unaffected.
Every parity comparison against Perl 0.6.11 goes through `gzip -dc` (`:385-386`, `:397-398`,
`:416-417`, `:766-767`, `:776-777`, `:788-789`). The `--clump_only` path writes through its own
`BufWriter<File>` (`clump_only.rs:297-298`), not `FastqWriter`, so it is untouched. The two
`std::fs::read` sites in `tests/integration_gzip_non_gz_extension.rs` (`:388`, `:492`) assert only
the `1f 8b` magic. `--clump_only`'s and `process_single_batch`'s own encoders already call
`finish()` explicitly (`clump_only.rs:224`, `parallel.rs:1076`), so no other teardown error is
being discarded. The only external exposure is nf-core/modules' nf-test snapshots, already known.

**Is the enum the right design, and did #434 over-specify?** The enum is serviceable and, per the
measurement above, free. A cheaper shape exists — a local object-safe trait, which keeps today's
dispatch untouched:

```rust
trait Sink: Write + Send {
    fn finish_sink(self: Box<Self>) -> Result<()>;
}
```

`self: Box<Self>` receivers are object-safe, so `Box<dyn Sink>` works; `dyn Sink` implements its
`Write` supertrait, so `write_record` compiles unchanged and still dispatches through the vtable to
`BufWriter`'s overrides. That removes the need to re-implement `Write` at all (and with it M2's
forwarding question and the whole hot-path debate), at the cost of a trait plus three small impls.
So yes, #434 over-specified mildly by naming the enum — but since the hot-path cost measures as
nothing, this is a preference, not a fix, and I am not asking for it.

**Is the problem worth a type change at all?** Yes. The failure class is broader than a full disk
in the last few bytes — `EDQUOT` on a quota, `EIO`, an `--output_dir` on a filling volume — and the
pre-fix symptom is exit 0 with a trailerless `.gz`, which is the one failure mode a trimmer must
not have. It is not silent *corruption*, because a downstream `gzip -dc` does reject a missing
trailer; it is a downstream failure with the blame pointed at the wrong tool. Worth fixing, and the
fix is small.

---

## Recommendations, prioritised

| # | Priority | Action |
|---|---|---|
| H1 | **High** | Correct the byte-identity claim in `CHANGELOG.md:22-26` and the PR body — `ParGz` (specialty modes, `--demux`) changes bytes too; only plain output and the SE/PE `--cores N` trimming path are identical. Replacement text above. |
| H2 | Medium | Soften "every byte it owes was written" in `CHANGELOG.md:8` and `src/fastq.rs:722-724`: `close(2)` is still unchecked. |
| M2 | Medium | Correct the comment at `src/fastq.rs:750-754` — the test fails on the header, not the trailer, and does not discriminate fixed from broken. |
| M1 | Medium | Either add the `mem::forget` guard at `src/fastq.rs:611-613` or note in `finish`'s doc that `ParGz` may panic out through `gzp`'s `Drop`. Pre-existing; data safety is unaffected. |
| M3 | Low | Optional: tighten the five comment blocks listed above to current-state one-liners. Maintainer's call, since the surrounding code reads the same way. |
| L1–L4 | Low | Notes only; no issues. |

No change is required to the logic. H1 is the one I would block a merge on, because the CHANGELOG
is the artefact a downstream maintainer acts on.
