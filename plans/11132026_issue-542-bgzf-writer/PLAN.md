# PLAN — Replace `bgzf::MultithreadedWriter` on the parallel uBAM output path (#542)

**Status:** revised after dual plan review; awaiting an implementation trigger (no code written)
**Issue:** [#542](https://github.com/FelixKrueger/TrimGalore/issues/542)
**Base:** `dev` @ `d862f7d`
**Related:** #529 (landed, was the blocker), #541, #568 (docs follow-up, filed from this plan's thread accounting)
**Reviews:** `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — see §13 for what changed in response.

---

## 1. Goal

Remove `bgzf::MultithreadedWriter` from the parallel uBAM output path. Compression moves into
the existing trim workers, which hand the sink opaque, already-compressed BGZF bytes; the sink
becomes an ordered append plus a terminator.

This retires three defects at once:

| Finding | Why it goes away |
| --- | --- |
| **C-1** — a write error surfaces as a panic, then a double panic aborts | No pool, no `send().unwrap()`. Every write is an `io::Result` on the main thread. |
| **H-2** — `Drop` discards the finish error, leaving a truncated file that looks structurally valid | `finish()` is explicit and fallible; the EOF block is written only on the success path, before `PendingOutput::commit`. |
| Thread-count overrun | The pool's N deflaters + 1 writer disappear — exactly N+1 threads. |

Records, their order and the header are unchanged. The `.bam` container bytes change (see §8).

---

## 2. Context

### Placement

Three files, plus tests:

| File | Change |
| --- | --- |
| `src/bam.rs` | New generic `BgzfSink<W: Write>` holding the compression and terminator logic; `ParallelBamSink` becomes a thin wrapper pairing it with `PendingOutput`. `workers` param dropped. New `BgzfBlocks` newtype + `BatchEncoder::into_bgzf`. |
| `src/parallel.rs` | `process_{single,paired}_batch_bam` compress before returning; the two `sink.write_encoded` call sites take the newtype. |
| `src/main.rs` | three `ParallelBamSink::create` calls lose `cli.cores`; the stale thread comment at `:2700`. |

### Doc comments that become false

Seven sites, all of which currently describe a compression pool that will not exist. Every one is
rewritten to the new current state — never to narrate the change:

`src/bam.rs:1109-1110` (struct doc, "compressed to BGZF by a pool of threads"), `:1112-1113`
(field-order rationale — the order still matters, but not because "the BGZF stream is closed"),
`:1120` ("`workers` is the number of compression threads"), `:1156` ("stop the compression
threads"); `src/parallel.rs:354-355` ("the sink compresses them to BGZF on its own pool of `cores`
threads") and `:460-462` (the single-end twin); `src/main.rs:2700` ("N threads convert the input
records, N workers trim and encode, N threads compress" — the third clause goes).

`src/parallel.rs:74` describes the FASTQ path and stays. `src/bam.rs:1853` is about the reader
pools and stays.

### Not in scope

- **`clump_only.rs` is unaffected.** `--clump_only --output-format ubam` uses the serial
  `BamWriter` (`src/clump_only.rs:698`, `:943`), not `ParallelBamSink`. The CI jobs at
  `.github/workflows/ci.yml:1075-1100` exercise that path and are untouched by this change.
- **The serial `BamWriter` is unaffected.** It wraps `bgzf::Writer` directly and already
  finalises through a fallible `try_finish`.
- **Thread-model documentation.** Decided: out of scope, filed as #568.

### Existing pattern this mirrors

The parallel FASTQ path already does exactly this — a scoped `GzEncoder` per chunk in the worker
(`src/parallel.rs:645-661`), the main thread appending opaque bytes (`:250-251`, `:1036`). BGZF is
a gzip-member stream, so concatenating per-batch members is valid, the same way the FASTQ path
concatenates gzip members.

---

## 3. Facts established before planning

All four were checked rather than assumed. **The first contradicts the design sketch in #542.**

### 3.1 The issue's code sketch appends an EOF block to every batch

#542 proposes a scoped writer and comments `// buf holds complete BGZF blocks, no EOF marker`.
That is false. `bgzf::Writer::drop` is:

```rust
fn drop(&mut self) {
    if self.inner.is_some() { let _ = self.try_finish(); }
}
```

`try_finish` appends `BGZF_EOF`. Leaving the scope therefore terminates every batch. Measured
(probe run against `noodles-bgzf-0.35.0`, then reverted):

| Shape | Bytes | Ends with EOF block |
| --- | --- | --- |
| scoped `write_all` + `flush`, dropped at scope end | 66 | **yes** |
| `write_all` + `flush` + `into_inner()` | 38 | no |
| fresh empty writer, `try_finish()` alone | 28 | is exactly the EOF block |

The 66 and 38 are for one short payload and scale with it; only the 28 is a constant.

**`into_inner()` is load-bearing, not stylistic**: it takes `inner`, leaving `None`, which disarms
`Drop`. The plan uses `flush()` then `into_inner()`.

**Severity, corrected.** The first draft called a stray per-batch EOF the headline hazard. Both
reviewers measured it at ~28 bytes per ~433 KB batch — about 0.01% — and confirmed samtools, bgzip,
gzip and noodles all read such a file correctly. It is a cosmetic defect. The dangerous mistake is
the **opposite** one: omitting the `flush()`, which on a sub-64 KiB batch makes `into_inner()`
return an *empty* `Vec` — a structurally valid, `quickcheck`-clean BAM with every record of that
batch missing, exit 0. §9 is ranked on that basis.

### 3.2 The terminator needs no private constant

#542 says `BGZF_EOF` being `pub(crate)` means we need "either a final empty writer's `try_finish`
or our own copy of the 28 bytes". The first option works: a fresh writer over the output, with an
empty staging buffer, emits exactly those 28 bytes and nothing else (row 3 above). **No copied
constant.** A hardcoded copy would be a second source of truth for a format detail that noodles
already owns.

### 3.3 Compression level is already equal on both sides — nothing to measure

#542 asks for the level to be compared rather than assumed. Reading settles it:
`bgzf::writer::Builder` and `bgzf::multithreaded_writer::Builder` both default to
`CompressionLevel::default()`, which is `6` (`writer/compression_level.rs:93-97`), and #529 sets
only `worker_count`. The serial `BamWriter`'s `bgzf::Writer` is the same default. All three paths
are level 6 before and after.

### 3.4 No byte-identity invariant exists on the `.bam` container

This is what makes the change safe to land at all:

- `assert_ubam_eq` (`tests/integration_ubam_out.rs:125`) compares the header minus `@PG` and
  decoded record tuples — not file bytes.
- The CI determinism checks use `samtools view | sort | md5sum` — record bodies.
- The Perl-parity validation matrix covers FASTQ output only and runs `--cores 1`.
- **The existing serial-vs-parallel parity tests already compare inflated streams**
  (`src/parallel.rs:2881`), under the comment *"the BGZF framing may differ; the decompressed BAM
  stream must not"*. Framing divergence is therefore tolerated **by design**, and this is the real
  data-loss guard — it would catch an omitted `flush()`, because its fixtures are 10–12k records.

So the golden fixtures and the `--cores 1` vs `--cores 4` comparison at
`tests/integration_ubam_out.rs:357` stay valid unchanged. Nothing in the repo asserts a `.bam`
byte layout.

### 3.5 Two things already true of the shipping code

Both found in review, both verified here, and both constrain the tests rather than the design:

- **The serial `BamWriter` emits two EOF blocks.** `bam::io::Writer::try_finish` delegates to
  `bgzf::Writer::try_finish` without taking `inner`, and `bam::io::Writer` has no `Drop`, so the
  field's drop fires `try_finish` again. Confirmed on real output — the last 56 bytes of
  `ubam_test_trimmed.bam` are two identical 28-byte terminators. Harmless, pre-existing, and **not
  in scope**; it means an "exactly one EOF block" assertion must never be pointed at serial output.
- **`ParallelBamSink` has no CI coverage at all.** Every `--cores` + `--output-format ubam` line in
  `ci.yml` is `--clump_only`, which uses the serial writer. No htslib ever reads this sink's
  framing. §9.2a adds that.

### 3.6 The error path is not masked

Checked because it decides whether a guard is needed: on a sink write failure the main loop sets
`first_error` and breaks, dropping `result_rx`; the reader's sends then fail, but
`read_{single,pairs}_round_robin` `break` and return `Ok(())` (`src/parallel.rs:1226-1242`). So
`reader_join_err.or(first_error)` surfaces the real write error rather than a generic reader error.
No change needed.

---

## 4. Signatures

```rust
/// BGZF-compressed bytes: whole members, no EOF block.
///
/// A newtype because `ParallelBamSink::write_compressed` appends its argument verbatim.
/// Uncompressed records passed there would produce a corrupt file with no error, so the
/// boundary is typed rather than documented.
pub struct BgzfBlocks(Vec<u8>);

impl BgzfBlocks {
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
}

impl BatchEncoder {
    /// Encode, then compress to BGZF members. Replaces `into_bytes` on the parallel path.
    pub fn into_bgzf(self) -> Result<BgzfBlocks>;
}

/// The BGZF layer over any writer: header member, ordered appends, terminator.
///
/// Generic so a test can inject a writer that fails. `ParallelBamSink::create` builds its
/// own `File` from a `Path`, which leaves the write-error path — the defect #542 exists to
/// fix — with no way to be exercised.
pub struct BgzfSink<W: Write> {
    inner: W,
    header: Header,
}

impl<W: Write> BgzfSink<W> {
    /// Writes the header as its own BGZF member.
    pub fn new(inner: W, header: Header) -> Result<Self>;
    pub fn header(&self) -> &Header;
    /// Append one batch's BGZF members. Callers must append in input order.
    pub fn write_compressed(&mut self, blocks: &BgzfBlocks) -> Result<()>;
    /// Append the EOF block and flush. Does NOT publish — that is the wrapper's job.
    pub fn finish(self) -> Result<W>;
}

impl ParallelBamSink {
    /// `workers` is gone — compression happens in the caller's workers.
    pub fn create<P: AsRef<Path>>(
        path: P,
        source_header: Option<&Header>,
        command_line: &str,
    ) -> Result<Self>;

    pub fn header(&self) -> &Header;
    pub fn write_compressed(&mut self, blocks: &BgzfBlocks) -> Result<()>;

    /// Terminate the BGZF stream, then publish under the final name.
    pub fn finish(self) -> Result<()>;
}
```

`ParallelBamSink` keeps its existing public shape — the three `main.rs` call sites change only by
dropping an argument — and holds `BgzfSink<BufWriter<File>>` plus `Option<PendingOutput>`. The
split exists so §9.3 is constructible; it is not refactoring for its own sake.

`BatchEncoder::into_bytes` stays — `§9.1`'s test needs the uncompressed bytes, and removing a
public method is a wider change than this issue needs.

---

## 5. Behaviour

1. **Create.** `PendingOutput::create` as now; wrap the handle in `BufWriter<File>`. Build the
   header, encode it through a throwaway `bam::io::Writer<Vec<u8>>`, then compress those bytes
   through a scoped `bgzf::Writer` and append — same `flush` + `into_inner` discipline as a batch,
   so the header is its own BGZF member.
2. **Per batch, in the worker.** After the last record is pushed, `into_bgzf()`: write the encoded
   bytes into a `bgzf::Writer` over a fresh `Vec`, `flush()`, `into_inner()`. Returns whole members
   with no EOF block.
3. **Per batch, on the main thread.** `write_compressed` appends verbatim. The existing reorder
   buffer already delivers batches in sequence, so the sink is only ever an append — no ordering
   work is added or needed.
4. **Finish.** A fresh `bgzf::Writer` over the inner writer, `try_finish()`, `into_inner()` →
   appends exactly the EOF block. Then `BufWriter::flush()`, then `PendingOutput::commit()`.
5. **Error.** Any step returns `Err` up the existing `first_error` path. `PendingOutput` is still
   `Some`, so its `Drop` removes the temporary and no file is published.

### Edge cases

| Case | Handling |
| --- | --- |
| **Empty batch** (every read filtered out) | `into_bgzf` on an empty encoder: `flush()` is a no-op on an empty staging buffer, `into_inner()` returns an empty `Vec`. `write_compressed` writes 0 bytes. No stray member. |
| **Zero records overall** | Header member + EOF block. Same as the serial writer on empty input; `assert_ubam_eq` reads it as a valid 0-record BAM. |
| **Batch > 64 KiB uncompressed** (always, at 4096 reads) | `bgzf::Writer::write` flushes a block whenever staging fills, so a batch becomes several members. Nothing special. |
| **A worker fails mid-batch** | `enc.push` returns `Err` before `into_bgzf` is reached; the partial buffer is dropped. Unchanged from today. |
| **Write error (`ENOSPC`/`EPIPE`)** | Returned as `Err` from `write_compressed` on the main thread — the defect this issue exists to fix. Today it panics in a pool thread and double-panics into an abort. |
| **`--cores 1`** | Does not reach this path at all; `cli.cores > 1` gates it in all three dispatch functions. |

---

## 6. Implementation outline

1. **`src/bam.rs` — add `BgzfBlocks` and *two* helpers, because there are two distinct operations.**
   - `fn bgzf_members(raw: &[u8]) -> io::Result<Vec<u8>>` — `write_all` → `flush` → `into_inner`.
     Emits members, never a terminator. Both the header and every batch go through it.
   - `fn bgzf_terminator<W: Write>(w: &mut W) -> io::Result<()>` — a fresh empty `bgzf::Writer`,
     `try_finish()`, `into_inner()`. Emits the terminator, nothing else.

   The first draft claimed one helper meant "exactly one place that can omit the flush". That was
   wrong: the terminator never went through it, so there were two paths and only one was guarded.
   Splitting them makes each one's postcondition stateable and testable in isolation — and getting
   the terminator wrong is invisible, since `try_finish()` followed by a drop emits *two* EOF
   blocks, which is already what the serial writer does (§3.5) and which no reader complains about.

   One-line comment on each, stating the postcondition — not the reasoning, and not the `Drop`
   mechanics, which belong in the commit message.
2. **`src/bam.rs` — `BatchEncoder::into_bgzf`.** `self.into_bytes()` → `bgzf_members` → `BgzfBlocks`.
3. **`src/bam.rs` — add `BgzfSink<W: Write>` and reduce `ParallelBamSink` to a wrapper.**
   `BgzfSink::new` writes the header through `bgzf_members`; `write_compressed` appends;
   `finish` calls `bgzf_terminator`, flushes and returns `W`. `ParallelBamSink` holds
   `BgzfSink<BufWriter<File>>` + `Option<PendingOutput>`, drops the `workers` parameter and the
   `NonZeroUsize` dance, and its `finish` calls through then commits. Field order still matters —
   rewrite the `:1112` comment to say the writer is flushed before the temporary is removed, since
   "the BGZF stream is closed" stops being accurate.
4. **`src/parallel.rs` — workers compress.** `BamBatchResult.bytes` and
   `SingleBamBatchResult.bytes` become `BgzfBlocks`; `process_*_batch_bam` call `into_bgzf()?`.
   The two main-thread sites call `sink.write_compressed(&r.bytes)`.
5. **`src/main.rs` — three call sites.** Drop the `cli.cores` argument at `:2381`, `:2556`, `:2703`.
   Rewrite the now-wrong comment at `:2702` (`"N threads convert the input records, N workers trim
   and encode, N threads compress"` — the third clause is gone).
6. **Tests.** Update `src/bam.rs:2733` and `src/parallel.rs:2861`, `:2937` for the dropped argument,
   then add §9.
7. **Changelog.** One sentence, `### Unreleased`, naming the panic-on-write-error fix and the N+1
   thread reduction. Per the register rule the thread figure is the entry's substance, so it stays;
   nothing else does.

---

## 7. Efficiency

- **CPU: this is a real reduction in parallel capacity, not a redistribution.** The first draft had
  this backwards. Today N trim workers *and* N deflaters both do work — 2N threads. Afterwards N
  threads do both jobs. Level-6 BGZF measures ~14.1 ms per 1.17 MB batch (~85 MB/s/core) against
  3.1 ms at level 1, so the compression being absorbed is not cheap.
- **The FASTQ-path precedent is weaker than §2 claims.** That path compresses at
  `DEFAULT_GZIP_LEVEL = 1` (`src/fastq.rs:29`); uBAM output compresses at 6. The *structure*
  transfers — per-chunk member, ordered append — but the per-byte cost does not.
- **Decision taken:** keep level 6 and accept a throughput regression on machines with spare cores,
  because `--cores N` promises a thread budget the current code exceeds by about 2x, and the
  threading docs tell nf-core users to pass `--cores task.cpus`. A run that honours its allocation
  is the point of the change. §9.5 measures the cost and the changelog records it; dropping the
  level to 1 was considered and rejected, as it would change output size for existing users.
- **Memory:** lower. The pool's bounded channels of `BytesMut` blocks go away; each worker gains one
  `Vec` roughly half the encoded batch size.
- **Allocation:** `into_bgzf` allocates one `Vec` per batch. A reusable per-worker buffer would
  avoid that, but it means threading a scratch buffer through `process_*_batch_bam`, and the
  existing FASTQ path allocates per batch too. Not worth the coupling; revisit only if §9.5 shows
  a regression.
- **File size:** a small increase. Each batch boundary forces a block boundary — the same trade the
  FASTQ path already makes. Quantified in §9.4.

---

## 8. Integration

- **Reads/writes:** same output path, same `PendingOutput` publication, same header.
- **Downstream:** records and order unchanged, so MultiQC, `samtools` and the golden fixtures are
  unaffected. The container's block layout changes; §3.4 establishes nothing depends on it.
- **Public API:** `ParallelBamSink::create` loses a parameter and `write_encoded` is renamed. Both
  are `pub` on a library crate, so this is a breaking library change — consistent with the crate's
  own framing of `--output-format ubam` as v1 of that surface, and it has no external callers in
  this repo.

---

## 9. Validation

Reviewer note: the first draft's 9.1 and 9.3 were both defective — 9.1 could not fail on the
mistake it existed to catch, and 9.3 could not be built. Both are rewritten.

| # | What to verify | How | Expected |
| --- | --- | --- | --- |
| **9.1** | **A batch round-trips** — the flush-omission mode of §3.1, which silently drops a whole batch. | Unit test in `src/bam.rs`: `BatchEncoder::into_bgzf` over a **small** batch (a few records, well under 64 KiB — the size at which an omitted flush yields an empty buffer), then inflate the bytes with `bgzf::Reader` and assert the decoded records equal what went in. **A positive round-trip, not a negated `ends_with`:** `ends_with` is false for an over-long needle, so "does not end with the EOF block" is also true of an *empty* buffer, which is exactly what the bug produces. | Records recovered intact. An omitted `flush()` fails here on a count of 0. |
| **9.1a** | No stray terminator mid-stream. The cosmetic half of §3.1. | Same test: assert the batch bytes do not contain the 28-byte EOF block **and** are non-empty. | Absent. Secondary to 9.1 — a stray block costs ~0.01% and every reader tolerates it. |
| **9.2** | Output is readable and record-identical to serial. | The existing comparisons at `src/parallel.rs:2861`/`:2937` already inflate both streams and compare (§3.4); they need only the dropped-argument edit. **Do not add an "exactly one EOF block" assertion to the serial side** — it has two (§3.5). | Decompressed streams equal; stats equal. |
| **9.2a** | **An htslib actually reads the new framing.** The sink has no CI coverage today (§3.5). | Add a CI step running `--output-format ubam --cores 4` on `test_files/ubam_test.bam`, then `samtools quickcheck` + `samtools view \| sort \| md5sum` against a `--cores 1` run. | Both md5s equal; `quickcheck` clean. |
| **9.3** | **A write failure reports an error instead of aborting**, and publishes nothing. | Now constructible via the §4 split: build a `BgzfSink<W>` over a `W` whose `write` returns `ErrorKind::Other` after N bytes, and assert `Err` rather than a panic. Separately, drive `ParallelBamSink` against a path under a read-only directory and assert the final path does not exist. (`unwritable_file` at `src/fastq.rs:814` is private to that module's tests and cannot be reused.) | `Err` surfaced, no panic, no abort; nothing published. **This is the guard for the defect the issue exists for** — without it the fix ships unverified. |
| **9.4** | File-size cost of per-batch block boundaries. | Trim a 1M×150bp input at `--cores 4` before and after; compare `.bam` sizes. | **Threshold stated up front: ≤2%.** Above that, batches are too small and the result goes in the plan before proceeding. |
| **9.5** | Throughput cost, which §7 now expects to be real rather than nil. | Same input, `--cores 4` and `8`, wall-clock, three runs each, on a machine with ≥2N free cores — the worst case for this change. **md5 the binaries first**; a stale binary has produced a fake benchmark in this repo before. | A regression up to ~25% at `--cores 8` is accepted per §7 and recorded in the changelog. Beyond that, bring the number back before merging. |
| **9.6** | **Thread count fell by exactly N+1** — measured, not derived. | `--cores 8`, uBAM in → uBAM out SE, sample thread count during the run. | 29 → 20 (3N+5 → 2N+4). The before-figure cross-checks the method against #542's measurement. |
| **9.7** | Whole suite, formatting, lints. | `cargo test`, `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`. | 876+ pass; clean. Local rustc 1.95 vs CI 1.99 — necessary, not sufficient. |

---

## 10. Assumptions

**Fixed (verified, §3):** level 6 both sides; `into_inner()` suppresses the EOF `Drop`; a bare
`try_finish()` emits exactly the EOF block; no byte-identity invariant on the container;
`clump_only` does not use this sink; the reorder buffer already guarantees order.

**Fixed (by construction):** BGZF members concatenate (RFC 1952, and the FASTQ path relies on it);
`--cores 1` never reaches this path.

**Assumed, not verified:** that no external consumer depends on `ParallelBamSink`'s current
signature. Checked within this repo only.

**Configurable, deliberately left alone:** the compression level stays the noodles default. #542
raises no case for exposing it, and `--high_compression` was removed in v2.2.0.

---

## 11. Questions and ambiguities

- **Critical — none outstanding.** Four decisions were put to Felix and answered: thread-model docs
  out of scope → #568; #546 closed with a summary comment; the sink is split into
  `BgzfSink<W: Write>` so the write-error path is testable; and the throughput regression is
  accepted and recorded rather than bought back by lowering the compression level.
- **Open (recorded, assumption taken):** whether `BatchEncoder::into_bytes` should be removed now
  that only tests use it on this path. Taken: keep it — §9.1 needs it and removing a public method
  widens the diff past the issue.
- **Open:** whether 4096 reads/batch is still the right size once a batch is also a compression
  unit. Taken: leave it. §9.4 will show whether the block-boundary cost argues otherwise, and
  changing it would move a figure the FASTQ path shares.
- **Open — #542's fourth ask, dropped from the first draft.** The issue closes by noting that
  `--cores 1` and `--cores N` would still carry two compression implementations obliged to produce
  equivalent output, and calls that "a standing liability". The §4 split narrows it — both sides now
  go through `bgzf_members` — but the serial `BamWriter` still finalises via `bam::io::Writer`, and
  that is the asymmetry which gives it two EOF blocks (§3.5). Taken: out of scope here; worth its
  own issue once this lands, which would retire the duplicate terminator at the same time.

---

## 12. Self-review

**Checked and adjusted:**

- **Logic.** Traced create → per-batch → finish → commit against the error path. Caught that
  `finish` must flush the `BufWriter` *after* the terminator and *before* `commit`, since the EOF
  block is written through a short-lived `bgzf::Writer` that does not own the `BufWriter`.
- **The sketch in the issue is wrong** (§3.1) — but it is the *cheap* error. Copying it puts a
  terminator in every batch, costing ~0.01% and nothing else; every reader tested handles it. The
  expensive error is the mirror image, omitting the `flush()`, which drops a whole small batch
  while leaving a file that passes `quickcheck`. The first draft ranked these the wrong way round
  and aimed its headline test at the harmless one; §9.1 now tests the round-trip instead.
- **Silent-corruption boundary.** `write_encoded` taking `&[u8]` means a caller that passes
  uncompressed bytes gets a corrupt file and no error — the sink can no longer tell. Added the
  `BgzfBlocks` newtype rather than relying on the rename. This is the same class as §3.1: wrong
  bytes that still parse.
- **Edge cases.** Empty batch, zero records, oversized batch, mid-batch failure, write failure,
  `--cores 1`. The empty-batch case decided the `flush()`-then-`into_inner()` shape over
  `try_finish`: on an empty encoder the latter would emit a stray EOF block where the former emits
  nothing.
- **Integration.** Confirmed `clump_only` and the serial writer are out of the blast radius, and
  that no test or CI job asserts container bytes (§3.4) — without which this change could not land
  as-is. The corollary, found in review, is that no CI job reads this sink's output at all (§3.5);
  hence 9.2a.
- **Testability as a design constraint.** The reason §4 grew `BgzfSink<W: Write>` is that the
  path-taking API left the issue's headline defect unobservable. A fix for a write-error bug that
  cannot be shown to handle a write error is not finished work.

**Remaining risks:**

- **9.4/9.5 are measurements, and could argue for a larger batch size.** That would be a follow-up,
  not a reason to hold this.
- **The validation matrix cannot see this path** — it is FASTQ-only and `--cores 1`, and CI does not
  exercise this sink at all (§3.5). The §9 tests plus the new 9.2a step are the entire guard.
- **9.6's before-figure depends on a current binary.** The tree's `target/release/trim_galore` is
  from `6315e96`, three refactors behind `HEAD`; it must be rebuilt before being measured or
  benchmarked.

---

## 13. Revision — what the dual review changed

Reviewers A and B agreed on six points; both re-derived §3's four facts independently and upheld
all four, including that #542's own sketch is defective.

| Finding | Response |
| --- | --- |
| **§9.1 could not fail on the mistake it existed to catch** — a negated `ends_with` is satisfied by the empty buffer a missing `flush()` produces | Rewritten as a positive round-trip on a deliberately small batch; the stray-terminator check demoted to 9.1a |
| **§9.3 was not constructible** — `unwritable_file` is private to `fastq.rs`'s tests (at `:814`, not `:836`), and `create` takes a `Path` | §4 splits out `BgzfSink<W: Write>`; 9.3 now injects a failing writer |
| **Severity was inverted** — a stray per-batch EOF costs ~0.01% and reads correctly everywhere; the flush omission is the data-loss mode | §3.1 and §9 re-ranked |
| **§7 was wrong in direction** — 2N working threads become N, and the FASTQ precedent runs at level 1 against uBAM's 6 | §7 rewritten with the measured figures and an explicit accepted-regression decision |
| **The sink has zero CI coverage**; no htslib reads its framing | New step 9.2a |
| **Four doc comments unlisted, plus `bam.rs:1112` the plan said to keep** | §2 now lists all seven sites |
| **The terminator sat outside the shared helper**, so §6's "exactly one place that can omit the flush" was false | §6 step 1 now specifies two helpers, `bgzf_members` and `bgzf_terminator`, one postcondition each |

Also folded in: the existing parity tests already compare inflated streams and tolerate framing
divergence by design (§3.4) — which is what re-ranked §9; the serial writer's pre-existing double
EOF block (§3.5), which constrains what 9.2 may assert; the error path is not masked (§3.6); the
66/38 byte figures are payload-dependent; and #542's fourth ask is recorded in §11.

Unresolved between the reviewers: only severity wording — A called the plan "not safe to implement
as written", B "safe, but not as written", on the same two blocking items. B ranked §9.3 Critical
against A's Important, on the grounds that it shapes §4; that reading was taken.

---

## 14. Implementation notes

Implemented 2026-11-13 on `dev` @ `d862f7d`. Tracked diff: `src/bam.rs`, `src/parallel.rs`,
`src/main.rs`, `.github/workflows/ci.yml`, `CHANGELOG.md` + the docs copy. **881 tests pass**
(876 + 5 new); `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -D warnings`
clean.

### Mutation testing — the part that matters

Five new tests all passed on first run, which in this repo means nothing until each is shown to fail
on the bug it guards. Each mutation was applied to a pristine snapshot and reverted; `src/bam.rs`
was confirmed byte-identical (md5 `54b01e67…`) afterwards.

| Mutation | Test that must die | Result |
| --- | --- | --- |
| `flush()` removed from `bgzf_members` | `into_bgzf_round_trips_a_small_batch` | **killed** — "compressed to 0 bytes", the exact mode the first draft's §9.1 would have passed |
| `into_inner()` → scope drop (the #542 sketch) | `into_bgzf_emits_no_terminator` | **killed** — "a batch carries a BGZF end-of-file block". The round-trip test *passed* under this mutation, confirming the two cover distinct modes |
| `bgzf_terminator` call omitted from `finish` | `parallel_sink_terminates_the_file_exactly_once` | **killed** — left 0, right 1 |
| write error swallowed in `write_compressed` | `write_error_is_returned_not_panicked` | **killed** |

### Measurements

| # | Result |
| --- | --- |
| **9.6** | **29 → 20 peak threads** at `--cores 8`, uBAM in → uBAM out SE, both measured on a baseline build of `d862f7d` and the current build (63 and 57 samples). Matches the predicted 3N+5 → 2N+4 exactly; reduction N+1 = 9. |
| **9.4** | **+0.045%** (44,869,908 → 44,890,204 bytes on 500 k reads), against a 2% threshold. |
| **9.2a** | Two separate runs, conflated in the first draft of this table and corrected here. **CI fixture** (`BS-seq_10K_R1`, 9,996 records, 218,514 bytes): record-body md5 identical at `--cores 1` and `--cores 4` (`24fb8d5f…`). **500 k-read input** (43 MB): `samtools view -c` decodes **496,150** records. First time real htslib has read this sink's output. |
| **9.5** | **Measured** after reclaiming disk. `--cores 4`: median 1.373 s → 1.577 s, **+14.9%** (best-of-5 +19.5%). `--cores 8`: 1.25–1.29 s → 1.26–1.28 s, **no measurable change**. Both within §7's accepted ≤25%. |

Measurements were re-taken against the final binary (md5 `b81af3fa…`) after the fmt and lint fixes,
because the first set used a build from before them.

### Deviations from the plan

1. **`samtools quickcheck` is useless here and was dropped from 9.2a.** It exits 8 with "had no
   targets in header" on *any* unaligned BAM — the very property #542 notes it cannot see past. The
   step instead asserts `samtools view -c` > 0, equal record-body md5s across `--cores`, and exactly
   one terminator at end-of-file.
2. **9.4 was measured as serial-vs-parallel on one build**, not before-vs-after binaries. The
   quantity of interest is the per-batch block-boundary cost, which `--cores 1` against `--cores 8`
   isolates directly, and it avoids a second multi-GB target directory.
3. **9.5 (throughput) was deferred, then measured** once `target/debug` and the leftover
   `target_repro_b/` were reclaimed (5.1 → 6.5 GiB free). No second target directory was needed:
   the baseline binary was built by restoring `d862f7d`'s three source files in place and rebuilding,
   keeping only the two 8 MB binaries. **Caveat on the figures:** this machine has 10 physical cores,
   so `--cores 8`'s worst case (≥2N free cores) is not representable — the old code's 29 threads were
   already oversubscribed there, which is why the regression shows at `--cores 4` (where the old
   code's 8 working threads fitted) and disappears at `--cores 8`. On a ≥16-core box the `--cores 8`
   regression would be real.
4. **`ParallelBamSink::inner` is `Option<BgzfSink<BufWriter<File>>>`.** `BgzfSink::finish` consumes
   itself to return `W`, so the wrapper must move it out; the plan's §4 sketch implied a plain
   field. `header()` and `write_compressed` carry an `expect` documented as unreachable through the
   public API.
5. **`BgzfBlocks::as_bytes` added** beyond the planned `len`/`is_empty` — the sink and the framing
   tests both need the bytes.
6. **Seven doc-comment sites rewritten, as §2 listed.** One of my own new comments was caught in a
   late audit for narrating the state being fixed (it named `MultithreadedWriter` and its double
   panic); six more were trimmed for carrying plan-section references and rationale. Current state
   only.

### Iteration log

- **#1** Core change built first time; `cargo build --release` clean after the three-file edit.
- **#2** `cargo fmt` reflowed `BgzfSink::finish`; clippy then failed on `io_other_error` in the new
  `FailingWriter` test. Switched to `std::io::Error::other`. Clean.
- **#3** First thread measurement returned a peak of **0** — `ps` is sandbox-blocked, so the loop
  was counting denials, not threads. Re-run with the bypass: 20. A plain "0" would have read as a
  plausible result.
- **#4** The CI step's nested `python3 - <<'EOF'` inside a YAML block would have shipped
  ten-space-indented Python and died on `IndentationError`. Replaced with `python3 -c`, then ran the
  whole step's logic locally rather than trusting it.

---

## 15. Post-implementation review round

Two code reviewers plus a coverage audit, all independent. **All three: safe to commit, no Critical
or High defect.** Reports: `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`, `COVERAGE.md`
(`**Verdict:** INCOMPLETE`, 38 of 44 items DONE, no code gap).

Both reviewers independently re-verified the suite rather than taking §14's word: 881 passed, fmt
clean, clippy exit 0. Reviewer B walked the BGZF member chain byte by byte (parallel: 27 members,
exactly one terminator at end of file; serial: 25 members, two terminators) and exercised the paired
parallel path and a zero-record input by hand — both correct.

### Fixed in this round

| Finding | Raised by | Fix |
| --- | --- | --- |
| `a_failed_run_publishes_nothing` **passed with `PendingOutput::drop` deleted** — the temporary is `.unpublished.bam.partial`, so the final path never existed either way | A (L9) and B (M-2), independently | Asserts the temporary exists before the drop and is gone after. Mutation 5 (drop removed) now kills it; it did not before. |
| **`finish()`'s error path had no test** — a mutation discarding it passed all five new tests, the integration suite and the new CI step, because every one finalises onto a writer that succeeds. This is defect H-2, which the change claims to retire | A (M1) | New `finish_error_is_returned` with a flush-failing writer. Mutation 6 kills it. |
| The `Option<BgzfSink>` and its three `expect`s are removable; **§14 deviation 4's reasoning was wrong** — `ParallelBamSink` has no `Drop` impl, so `finish(self)` can destructure | A (L7) and B (M-1), independently | `inner`/`pending` are now plain fields; `finish(self)` destructures. Three `expect`s gone. |
| `// flush stages the tail block` **inverts noodles' terminology** — `write` stages, `flush` emits | A (L4) and B (L-1) | Reworded. |
| `BgzfBlocks::len()` is dead code | A (L3) and B | Removed. |
| §14's 9.2a row **welded two runs together** — md5 `24fb8d5f…` is the 9,996-record/218 KB run, not the 496,150-record/43 MB one | A (M2) | Both runs now listed separately. Confirmed by re-measuring both. |
| §9.5's acceptance criterion said the regression would be "recorded in the changelog", and the entry recorded only the thread win; "29 threads" was also stated as fact though never re-measured | B (H-1), A (L8b) | 9.5 measured (above); the entry now carries both figures and is 4 lines, matching its longest siblings. |

### One reviewer finding rejected

**A's L8(c), "it links an issue where every sibling links a PR", is false.** Inside the Unreleased
section issue links outnumber PR links **76 to 24**, and the sibling entries for #541 and #544 — from
this same review cycle — both link issues. The entry matches its neighbours.

On the changelog the reviewers contradicted each other: A called it the longest in its section and
in the wrong register, B checked the house style and judged that it fits. Measured: the entry was 5
lines against siblings' 3–4, so A was right on length and B right on shape. It is now 4 lines.

### Deliberately not done (outside the approved scope)

The write error still does not name the output path (A L1) though `PendingOutput::tmp_path()` is
available; the library API break is unrecorded in the changelog (A L11, premise verified — `pub mod
bam` is exported and the crate publishes to crates.io); the empty-batch and zero-record cases have
no committed test (COVERAGE 2 and 3, though B verified both behave correctly); CI covers the
single-end parallel path only; and the thread table in #568 will need its numbers updated once this
lands. A's L8(d) — that a panic-to-error fix belongs under `#### Bug fixes` — is left as written.
