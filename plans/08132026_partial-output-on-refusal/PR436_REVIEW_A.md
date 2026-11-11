# PR #436 — Code Review A

**PR:** #436 `fix(io): publish an output only if its teardown succeeded (#434)`
**Author:** `sujeito-operator` (fork, non-collaborator, disclosed as agent-written)
**Reviewed at:** `b8aa7d8`, parent `63f5d44` (= `dev`)
**Reviewer:** A (independent; Reviewer B in parallel, no shared state)
**Method:** own sha-pinned worktree at `/private/tmp/claude-501/pr436-a` with a private
`CARGO_TARGET_DIR`; baseline `dev` binary built separately at
`/private/tmp/claude-501/pr436-a-base` for byte comparison. Every number below was measured on
this machine, not inferred.

---

## Verdict

**Merge with changes.** No Critical. One High, which is a one-clause fix to a user-facing
statement that is factually wrong. Three Mediums, all small. The change itself is correct: I
could not break the invariant it claims, and the two parallel-path changes are sound including
the deadlock question.

**Top three findings**

1. **High-1 — the output-byte change is not serial-gzip-only.** The `ParGz` arm shifts too, and
   both the PR body and the CHANGELOG say it does not. Measured: `--hardtrim5 30 --cores 2`
   gives 209,305 bytes on `dev` and 209,294 on the patch, decompressed md5 identical.
2. **Medium-1 — `Sink::ParGz`'s error path panics rather than returning the error**, so the
   PR's central claim ("`finish()` runs the concrete teardown and returns its error") holds for
   two of three arms. Safety is unaffected; the promise is not delivered.
3. **Medium-3 — the comments violate the house rule, and repeat a defect #428's own review
   round already corrected** (PLAN §12: "Three comments exceeded the house rule"). Five sites,
   mostly narrating the state being fixed.

---

## Gates — verified, not assumed

| Gate | Result |
|---|---|
| `cargo test --release` | **710 passed, 0 failed**, 16 suites (441 lib + 269 integration) |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean (exit 0; a real 15 s check, not a cache hit — it recompiled the dep graph) |
| `cargo build --release` | clean |

### The test-count discrepancy is fully resolved, and nobody was wrong

The author reports 708 (440 lib + 268 integration); the requester measured 710; the brief
guessed Linux would give 709. The platform delta is **+2 for macOS**, not +1:

| Gate | Platform-only tests |
|---|---|
| macOS-only | `src/format.rs:689`, `src/format.rs:711` (lib), `tests/integration_non_restartable_input.rs:252` (integration) |
| Linux-only | `src/cli.rs:1950` (lib) |

So macOS lib = Linux lib + 1 → Linux lib 440 ✓ (author), and macOS integration = Linux
integration + 1 → Linux integration 268 ✓ (author). **710 (macOS) = 708 (Linux) + 2.** The
author's arithmetic is correct for Linux, the requester's for macOS. The brief's 709 was the
only wrong figure — it counted the lib-side delta and missed the macOS-only integration test.

---

## Attack-by-attack findings

### 1. Is the enum the right design, or did #434 over-specify?

**The problem is worth a type change.** It is the `close(2)`-can-fail problem: the symptom is
exit 0 plus an unreadable `.gz`, the window is narrow but the failure is silent, and the
alternative — documenting it — was already tried in #428 and withdrawn during review (PLAN §12:
"The docs claimed a full disk leaves nothing behind… Clause deleted"). Silent beats narrow.

**#434 did over-specify by naming the enum.** A cheaper shape exists and nobody costed it:

```rust
trait Sink: Write + Send {
    fn finish(self: Box<Self>) -> Result<()>;
}
// writer: Box<dyn Sink>
```

A consuming method on a trait object is legal via `self: Box<Self>`, so teardown becomes
callable **without touching the write path at all** — the vtable dispatch stays exactly as it
is, and the hot-loop question the author spent a benchmark section failing to resolve never
arises. Diff size is comparable.

**That said, the enum is not worse and should not be sent back.** Static dispatch lets
`write_all` inline per variant and drops a heap allocation, and `write_to` issues one
`write_all` per record over a ~250-byte buffer, so the dispatch cost is noise either way. The
author's benchmark section is honest about being inconclusive and is *right* that it does not
matter. Worth the maintainer knowing the trait alternative existed; not worth blocking on.

**No `Send` regression.** The old field was `Box<dyn Write + Send>`; the enum carries no bound.
I compile-probed it: `Sink` and `FastqWriter` are both still `Send` (probe appended to
`src/fastq.rs`, run, then reverted with `git checkout --`; worktree confirmed clean).

### 2. The output bytes — HIGH-1

**Serial gzip: claim verified, figure is fixture-dependent.**

| | bytes | decompressed md5 |
|---|---|---|
| `dev` @ `63f5d44` | 341,686 | `12fde5c8…` |
| `b8aa7d8` | 341,674 | `12fde5c8…` |

12 bytes on `test_files/BS-seq_10K_R1.fastq.gz`, where the author measured 11 on theirs. The
delta depends on the bit alignment of the final deflate block, so **"11 bytes" is not a
constant** and the CHANGELOG should not read as if it were. Decompressed output is genuinely
unchanged.

**HIGH-1: the `ParGz` arm changes too, and the PR says it does not.** Both the PR body and the
CHANGELOG assert "`--cores N` output is byte-identical either way".

| `--hardtrim5 30 --cores 2` | bytes | decompressed md5 |
|---|---|---|
| `dev` @ `63f5d44` | 209,305 | `3c9ca888…` |
| `b8aa7d8` | 209,294 | `3c9ca888…` |

Mechanism, read from gzp 2.0.2 and confirmed by the numbers above. `ParCompress`'s `Write::flush`
is `flush_last(false)` (`src/par/compress.rs:465-467`), and `is_last` selects the deflate flush
mode (`src/deflate.rs:96-100`: `FlushCompress::Finish` vs `FlushCompress::Sync`). The old
`FastqWriter::finish` called `self.writer.flush()?` first, so the tail chunk went out with
`Sync` and `Drop`'s `finish()` then added an empty `Finish` chunk. The new code encodes the tail
chunk with `Finish` directly — one sync marker and one empty block shorter.

**Which arms.** `Sink::ParGz` is reachable wherever `cores` is threaded into
`FastqWriter::create`: `specialty.rs:54, 96, 268, 269, 379, 380` (`--hardtrim5`, `--hardtrim3`,
`--clock`, `--implicon`) and `demux.rs:197, 203`. Every `main.rs` site hardcodes `1`
(`:1332, 1586, 1587, 1597, 1603, 1604, 1870, 1871, 1874, 1875`), and the `--cores N` *trimming*
path never builds a `Sink` at all — `parallel.rs` compresses per batch into a `Vec` and writes
through a raw `File`. **That is why the author's `--cores 2` measurement came back
byte-identical: they measured the one gzip path this change cannot reach.** The claim is a
sound observation generalised one step too far.

**No consumer of raw `.gz` bytes in this repo — checked, including the ones that would hide it.**
Every Perl-parity comparison goes through `gzip -dc` (`ci.yml:385-398, 416-417, 766-789`). The
two md5-stability steps compare a build against itself (`ci.yml:685-694` re-run;
`:836-837` clump-only determinism). The only tests that read raw `.gz` bytes assert the magic
number and nothing else (`tests/integration_gzip_non_gz_extension.rs:388-393, 492-495`). And
**no CI step runs a specialty mode or `--demux` with `--cores >= 2`** — the only `--cores`
usages are paired multi-pair and `--clump_only` (`ci.yml:697-704, 804-913`), and `clump_only`
owns its own writer rather than `FastqWriter`. So CI could not have caught the `ParGz` shift
even if it compared bytes. Nothing in-repo breaks; the gap is entirely in what users are told.

**Fix — `CHANGELOG.md:25-30`.** Replace:

> **One visible consequence, on the serial-gzip path only:** … `--cores N` output is byte-identical either way.

with:

> **One visible consequence, on the gzip output paths:** the old teardown reached
> `GzEncoder::flush`, which emits a deflate sync marker before the final block; the explicit
> `try_finish` does not. The `.gz` bytes therefore differ slightly from v2.3.0 — a dozen bytes
> shorter on a 1M-read output, with an identical decompressed md5. The specialty modes and
> `--demux` shift the same way under `--cores >= 2`; the `--cores N` trimming output does not
> change.

### 3. Deleting `impl Drop for FastqWriter` — claim verified

Its only body was `let _ = self.writer.flush()`. Per variant on the abandon path:

| Variant | What replaces the deleted flush |
|---|---|
| `Plain(BufWriter<File>)` | `BufWriter::drop` flushes the buffer |
| `Gz(BufWriter<GzEncoder<File>>)` | `BufWriter::drop` flushes into the encoder, then `GzEncoder::drop` does `let _ = try_finish()` — strictly more than the old flush |
| `ParGz(ParCompress)` | no `BufWriter`; `ParCompress::drop` → `finish()` → `flush_last(true)`, which flushes. Previously `flush()` ran first and then the same `Drop`; both write everything |

Field order still puts `writer` before `pending`, so declaration-order field drop closes the
sink before the temporary is removed. Removing the explicit `Drop` does not change that —
without a `Drop` impl the fields drop in the same order.

**Abandon path verified end-to-end on the `ParGz` arm**, which is the one with no unit
coverage. `--hardtrim5 10 --cores 2` on `test_files/ubam_ws_qname_late.bam`: exit **1**, the
whitespace-QNAME error with the record ordinal intact, and the output directory contains
**neither `*.10bp_5prime.fq.gz` nor any `.partial`** (checked with `ls -la`, dotfiles included).
The 21 tests in `tests/integration_partial_output.rs` cover the rest and pass.

### 4. The `ParGz` arm — MEDIUM-1

The comment at `fastq.rs:594-597` claims `finish()` "Consumes `self`, so the sink's own `Drop`
has nothing left to do". **True on the success path, false on the error path.** gzp 2.0.2,
`src/par/compress.rs:377-388`:

```rust
fn finish(&mut self) -> Result<W, GzpError> {
    self.flush_last(true)?;              // early return, BEFORE the takes
    drop(self.tx_compressor.take());
    drop(self.tx_writer.take());
    match self.handle.take().unwrap().join() { … }
}
```

`flush_last` fails only by a channel `send` failing (`:347-356`), which happens precisely when
the writer thread has already died — the ENOSPC case. So `finish()` returns `Err` with all
three fields still `Some`, and `Drop` (`:396-401`) tests exactly those three fields and calls
`self.finish().unwrap()`. The second `flush_last` fails the same way → **panic**.

Traced through: `w.finish()?` propagates, `w` drops at the early return, `ParCompress::drop`
panics, the panic unwinds through `FastqWriter::finish` where `pending` is a live local, so
`PendingOutput::drop` removes the temporary. **Nothing is published, which is the property that
matters** — and it is identical to `dev`, where `flush()?` then field-drop reached the same
panic. Not a regression.

But the PR body and CHANGELOG both say `finish()` "runs the concrete teardown and returns its
error", and on `--cores >= 2` with a full disk the user gets a panic and a backtrace note
instead of the clean `No space left on device` anyhow chain. There is no good way to engineer
around it: neutralising `Drop` needs `mem::forget` (leaks the compressor threads and their
channels) or `catch_unwind`, both worse than the panic. **Recommend documenting it, not fixing
it** — `fastq.rs:594-597`:

```rust
    /// Close the sink, writing whatever trailer it owes, and return the error.
    /// A failed `ParGz` teardown panics instead: `ParCompress::finish` errors
    /// before it disarms its own `Drop`.
```

The panic *cannot* escalate to an abort here, because the `?` early return is not itself an
unwind. The pre-existing abort case — `ParCompress::drop` panicking while already unwinding,
which leaves a `.partial` — is unchanged and is PLAN §4's "Abort during unwinding" row.

### 5. `Sink`'s `Write` impl — no issue

`write`, `write_all` and `flush` are forwarded; `write_vectored`, `write_fmt`,
`is_write_vectored` are not. Nothing reaches a defaulted method:

- `writer` is private and touched at exactly one place, `fastq.rs:716`
  (`record.write_to(&mut self.writer)`).
- `FastqRecord::write_to` (`fastq.rs:79-97`) formats the whole 4-line block into a `Vec` and
  issues **one** `write_all`. Nothing else.
- `write_vectored` / `write_fmt` / `is_write_vectored` are never called on any writer in the
  crate.

The forwarded `write_all` matters and the author is right to have written it: the default would
loop over `write`, which for `Sink` would re-match per iteration. `BufWriter` overrides
`write_all`; `ParCompress` does not, so it got the default through the old vtable too — no
change there.

### 6. The two parallel changes — no deadlock, no error masking

**Deadlock: cannot happen at either join.** In both functions the order is
`drop(result_rx)` → reader join → worker joins. By the time the worker joins run, (a) the
result channel is disconnected, so a worker parked in `rtx.send` returns `Err` immediately and
breaks, and (b) the reader has returned, and `txs` was moved into the reader closure, so every
work-channel sender is dropped and `rx.recv()` returns `Err`. Every worker is therefore
terminating or terminated. SE: `parallel.rs:1029, 1031, 1048`. Paired: `:291, 300, 318`.

An early `?` from `out.write_all` bypasses the explicit `drop(result_rx)` — but `result_rx` is
captured by move into the `thread::scope` closure (it is consumed by `drop(result_rx)`), so it
is dropped when the closure body exits, which is before `scope`'s own join loop. No deadlock
there either.

**Error precedence: matches the paired path, and no masking.** Both use
`reader_join_err.or(first_error)` (SE `:1038`, paired `:306`). The question is whether a
worker error can be masked by a reader error manufactured by the worker's own exit. It cannot:
`read_single_round_robin` treats a failed `send` as `break` and returns `Ok(())`
(`:1190-1192, 1204`), and `read_single_clumpy` does the same via `return Ok(false)`
(`:1226, 1252`). So `reader_join_err` is `None` in that scenario and the real worker error is
returned. The reader returns `Err` only from `reader.next_record()?`, i.e. a genuine input
fault, which correctly outranks a downstream worker error.

**No worker error can be lost.** `mpsc` delivers queued items before reporting disconnection,
and a worker drops its sender only after its `send` attempt, so a sent `Err` is always received
before the loop sees a closed channel.

**The worker-join fix is real, and slightly better than the PR claims.** I expected
`thread::scope` to re-raise and override the `bail!`; it does not. std takes the result out of
the `Packet` in `join`, so an explicitly joined panic is "handled" and never sets
`a_thread_panicked`. Measured with two standalone programs:

| Case | Outcome |
|---|---|
| one panicking worker, joined | closure's `Err` propagates cleanly — exit 1, `Worker thread panicked` |
| two panicking workers, `bail!` on the first leaves the second unjoined | `scope` panics `a scoped thread panicked` — exit 101 |

Commit is skipped in both. So the common case converts "commit a truncated output, then exit
101" into "publish nothing, exit 1 with a clean message" — a genuine improvement, and the
`bail!` message is reachable rather than dead.

**Note, below the filing bar: the SE path has no reader sentinel where the paired path does.**
On a reader error the SE main loop drains and writes every already-queued batch into the
temporary before it notices, where paired short-circuits on a sentinel `Err`
(`parallel.rs:195-206`). Pre-existing, and harmless — the wasted bytes go to a file that is then
unlinked. Recording it so it is not rediscovered as a defect.

### 7. #428's invariants — both hold

- **Commit where the writer is closed, never later.** `FastqWriter::finish` still does teardown
  then `pending.commit()` in one call, and every driver still calls `finish()` at the same
  place. The four in-run consumers of the final path (FastqQC at `main.rs:1459`/`:2241` and
  `clump_only.rs:374`/`:847`, `demux::demultiplex` at `main.rs:1474`, the `?`-propagating
  `std::fs::metadata` at `clump_only.rs:828`/`:1072`) all still run after `finish()` returns.
  The 710-test suite includes the T7/T10 cases that would fail on a commit-ordering slip.
- **`released` before the rename.** `io.rs` is not in the diff; `commit` at `io.rs:734-743`
  still sets `self.released = true` before `std::fs::rename`, so a failed rename preserves the
  data and `Drop` (`:751-757`) does not delete it. Re-read to confirm rather than assumed.

One ordering detail worth stating because it changed shape: the file descriptor is still closed
before the rename in all three arms — `Plain` discards the `File` returned by `into_inner`,
`Gz` drops the `GzEncoder` at the end of its match arm, `ParGz` discards the `W` returned by
`finish` — all inside `Sink::finish`, which returns before `commit()`.

### 8. The measured claims and the three new tests

Counts and gates: resolved above, all correct for their platform. The timing table is honestly
presented — the author states plainly that the spread within each build exceeds the gap between
them and that the box cannot resolve it, and declines to claim neutrality. That is the right
way to report an inconclusive measurement; I have no criticism of it, and per Question 1 the
question it fails to answer does not matter.

**Do the three tests exercise the failure they claim?** Mostly yes.

- `plain_sink_finish_reports_a_failing_flush` — genuine. 15 bytes into a 64 KB `BufWriter` stay
  buffered, so `into_inner`'s `flush_buf` is the first write to the read-only fd and it fails.
- `gz_sink_finish_reports_a_failing_teardown` — genuine, but the comment overstates it
  (**LOW-2**). With no records written, flate2 has not yet emitted the gzip header, so the
  failing write inside `try_finish` covers header + final block + trailer, not "the bytes that
  fail here are the ones ONLY `try_finish` writes". It does still isolate `try_finish` — nothing
  else wrote — so the test is sound and the sentence is not.
- `gz_sink_finish_writes_the_trailer_explicitly` — a weak but non-vacuous control
  (**LOW-3**). `assert_eq!(&bytes[len-8..], &[0u8; 8])` holds for an empty member's zero
  CRC/ISIZE but does not prove well-formedness. It does discriminate against the
  no-`try_finish` case: without it the file is 0 bytes and the slice panics.

**MEDIUM-2: `Sink::ParGz` has no test at all.** It is the arm whose bytes changed (High-1),
whose error path panics (Medium-1), and whose teardown is the only one that joins threads. Two
arms got a failure test and one got a control; `ParGz` got neither. A positive control is cheap
and would have caught High-1:

```rust
#[test]
fn pargz_sink_finish_writes_a_complete_member() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("par.gz");
    let mut sink = Sink::ParGz(
        ParCompressBuilder::<Gzip>::new()
            .num_threads(2)
            .unwrap()
            .compression_level(Compression::new(1))
            .from_writer(File::create(&path).unwrap()),
    );
    sink.write_all(b"@r\nACGT\n+\nIIII\n").unwrap();
    sink.finish().unwrap();

    let mut got = String::new();
    MultiGzDecoder::new(File::open(&path).unwrap())
        .read_to_string(&mut got)
        .unwrap();
    assert_eq!(got, "@r\nACGT\n+\nIIII\n");
}
```

A `MultiGzDecoder` round-trip is also the stronger form Low-3 wants; consider using it for the
`Gz` control too.

The author's three "what this does not claim" disclosures are accurate and I found nothing
hidden behind them.

---

## MEDIUM-3 — comment length and historical narrative

The house rule: one line, two maximum with a reason; state the fact, not the evidence; **current
state only, never the state being fixed**; never justify the change in the source. PLAN §12
records that #428's own review round trimmed three comments for exactly this — so this is a
repeat, and the same instinct (preserving hard-won reasoning in the wrong venue) is producing
it. The reasoning belongs in the commit message.

Five sites, with proposed replacements.

**`fastq.rs:570-586` — 17 lines, most of them the history of the boxed writer.**

```rust
/// The three concrete sinks a [`FastqWriter`] can write records into.
///
/// Named rather than boxed so each sink's teardown — the gzip trailer, the
/// parallel compressor's final blocks — is callable and its error propagable
/// before the output is published.
```

**`fastq.rs:722-727` — ends on pure narrative.** Delete the last sentence ("Before this,
teardown ran in the sink's own `Drop` and its error was discarded.") and the `(#434)` reference:

```rust
    /// The sink is torn down first and its error propagated, so the final name
    /// appears only once every byte the output owes — including the gzip
    /// trailer — has been written. A failed teardown returns with `pending`
    /// still armed, so its `Drop` removes the temporary and nothing is
    /// published.
```

**`parallel.rs:948-954` — 7 lines, 4 of them the old single-end behaviour.**

```rust
    // `Result<_>` in the channel so a worker error reaches the main thread,
    // which is the thread that decides whether to publish. Mirrors the paired path.
```

**`parallel.rs:310-317` and `:1042-1047` — the two join comments**, 8 and 6 lines, each
narrating the old bug and the old exit code. Both become:

```rust
        // Join the workers before publishing: `thread::scope` would otherwise
        // join them after the commit, letting a panicking worker publish a
        // short output.
```

**`fastq.rs:750-754` — the test comment** restates #434's mechanism in five lines. One line,
and it fixes Low-2 at the same time:

```rust
        // Only `try_finish` writes the deflate end-of-stream and the CRC/ISIZE trailer.
```

---

## LOW-1 — the CHANGELOG entry is 20 lines and two paragraphs

The standing rule is one sentence per entry, no rationale, no mechanism, no measurements. This
entry has all three. **But the repo's own file has multi-paragraph precedent immediately below
it** — the #398 entry runs ~20 lines with a second `**This refuses some invocations…**`
paragraph, and the #428 entry above it is three lines. So the local convention is already
looser than the rule, and cutting this one to a sentence while its neighbour stays an essay
would look arbitrary. **Maintainer's call.** My recommendation: cut the first paragraph to its
first sentence plus the `--cores N` clause, keep the second (byte-change) paragraph because
users act on it, and let the mechanism live in the commit message. Either way High-1's factual
correction has to land.

---

## Efficiency

Nothing in the hot loop regresses in a way I can substantiate, and nothing new allocates per
record. The per-record cost moves from one vtable indirection to one three-arm branch, on a
path that issues one `write_all` per ~250-byte record — noise. `worker_handles` adds one `Vec`
of `cores` elements per run. No extra syscalls. The teardown does slightly *less* work than
before on both gzip arms (one fewer sync-flush round), which is where the 11–12 bytes come
from.

---

## Nothing filed

Per the filing bar, none of the above is proposed as a new issue. High-1 and Medium-1/2/3 are
changes to this PR. The two notes (the SE reader sentinel; the pre-existing
abort-during-unwinding `.partial`) have no user-visible symptom I can name beyond what PLAN §4
already records, so they stay in this report.

---

## Summary table

| # | Priority | Finding | Fix |
|---|---|---|---|
| High-1 | High | `--cores N` byte-identity claim is false for `ParGz` (specialty + `--demux`, `--cores >= 2`); measured 209,305 → 209,294 | one clause in `CHANGELOG.md:25-30` and the matching PR-body line |
| Medium-1 | Medium | `Sink::ParGz`'s failed teardown panics rather than returning `Err`; `ParCompress::finish` errors before disarming its `Drop` | two-sentence comment at `fastq.rs:594-597`; do not engineer around it |
| Medium-2 | Medium | `Sink::ParGz` has no test — the arm with the byte change and the panicking error path | add the `MultiGzDecoder` round-trip control above |
| Medium-3 | Medium | Five comments exceed the house rule and narrate the state being fixed; a repeat of a #428 review finding | replacements above |
| Low-1 | Low | 20-line CHANGELOG entry vs the one-sentence rule; local precedent is looser | maintainer's call |
| Low-2 | Low | `gz_sink_finish_reports_a_failing_teardown`'s comment overstates which bytes fail | folded into Medium-3 |
| Low-3 | Low | Trailer control asserts 8 zero bytes rather than decoder integrity | use the round-trip form |
| — | note | SE parallel path has no reader sentinel where paired does; drains before noticing | none; nothing is published |

**No CI has run on this PR.** The four gates above were run by me on `b8aa7d8` in an isolated
worktree; the `validation` job (Perl 0.6.11 md5 parity) I did **not** run — it needs conda and
the pinned upstream script. My byte analysis says it cannot be affected, because every
comparison in it goes through `gzip -dc` and no step reaches a `ParGz` sink, but that is an
argument rather than a green run and should be treated as one.
