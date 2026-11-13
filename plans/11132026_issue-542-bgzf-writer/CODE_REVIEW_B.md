# CODE REVIEW B — #542, per-batch BGZF compression on the parallel uBAM output path

**Reviewer:** B (independent; A and a coverage auditor ran concurrently on the same tree)
**Target:** uncommitted working tree, `dev` @ `d862f7d`
**Files:** `src/bam.rs`, `src/parallel.rs`, `src/main.rs`, `.github/workflows/ci.yml`, `CHANGELOG.md` + docs copy
**Date:** 2026-11-13

---

## Verdict

**Safe to commit.** No Critical or High *defect*. The BGZF framing is correct and I verified it
against real htslib rather than inferring it; the error and abandon paths publish nothing and
surface the underlying `io::Error`; the five new tests are all capable of failing.

One **High** item is not a defect but a decision for you: plan §9.5 (throughput) was not measured,
and the CHANGELOG records only the thread-count win. Everything else is Medium and below.

---

## What I verified independently

Not taken from PLAN §14 — re-run or re-derived here.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | **exit 0, zero warnings** (my first attempt captured no status — zsh has no `PIPESTATUS`; re-run and captured properly) |
| `cargo test` | exit 0, **881 passed / 0 failed** summed across 25 test binaries — matches §14's claim exactly |
| Release binary | the on-disk binary was **stale** (11:43 vs `src/bam.rs` 11:45); rebuilt (md5 `7394b6b1…`) before any probe |
| SE, htslib 1.21, `--cores 4` vs `--cores 1` | 9,996 records both; record-body md5 `24fb8d5f03ae08e2abbe8b007356b02d` **identical**. Independently reproduces §14's `24fb8d5f…` |
| BGZF member chain, parallel file | walks cleanly to EOF: **27 members, exactly 1 EOF block, at end of file**. Header is its own 235-byte member |
| Same walk, serial file | 25 members, **2 EOF blocks** — pre-existing, confirms the asymmetry the tests encode |
| **PE** parallel vs serial (not covered by CI) | 19,992 records, 9,996/9,996 FREAD1/FREAD2, md5 **identical**, 1 terminator |
| **Zero-record** run (`--length 10000 --cores 4`) | 276-byte header-only BAM, 1 terminator, htslib reads 0 records, exit 0 |
| The new CI step's python assertion | **live and discriminating**: exits 1 with "found 2" on the serial file, exits 0 on the parallel file |
| File-size cost | **+0.047%** (218,416 → 218,518 B). §9.4 reported +0.045% |

Probe artefacts: `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/a2f653ef-8bfc-4b3c-9a4c-4b2b32c6aa15/scratchpad/revB_probe/`

### Framing correctness, from the dependency's source

Read `noodles-bgzf-0.35.0/src/writer.rs` rather than inferring:

- `bgzf_members` — `write_all` → `Writer::write` only calls `flush()` when `staging_buf` reaches
  `MAX_BUF_SIZE` (~65505 B); the explicit `flush()` emits the tail block; `into_inner()` takes
  `inner`, leaving `None`, so `Drop`'s `try_finish` is a no-op. **No terminator, no lost tail.**
- `bgzf_terminator` — fresh `Writer` over `&mut W`, empty staging buffer, so `try_finish()` writes
  `BGZF_EOF` and nothing else; `into_inner()` disarms `Drop`. **Exactly one terminator.**
- Compression level is unchanged: `writer::Builder::default()` and
  `multithreaded_writer::Builder::default()` both take `CompressionLevel::default()` (level 6).

### Drop order and the abandon path

`ParallelBamSink` declares `inner` before `pending` and has **no `Drop` impl**, so fields drop in
declaration order: `BufWriter<File>` flushes and closes, *then* `PendingOutput::drop` removes
`.<name>.partial`. Correct, and correct on Windows. On the `finish` error path the returned
`BufWriter` is dropped as a statement temporary before `commit()` renames, so the file is closed
before the rename. `BgzfSink::finish` calls `self.inner.flush()` explicitly, so the discarded
`BufWriter`'s error-swallowing `Drop` has an empty buffer and nothing is lost silently.

### Error masking (plan §3.6), re-derived

I checked this before reading §3.6 and reached the same conclusion. On a sink write failure the
main loop sets `first_error` and breaks; `drop(result_rx)` kills the workers; the reader's
`txs[idx].send(...)` then fails, but `read_single_round_robin` (`src/parallel.rs:1225-1241`) and
`read_pairs_round_robin` (`:947`) **`break` and return `Ok(())`** on a dead receiver. So
`reader_join_err` is `None` and `reader_join_err.or(first_error)` surfaces the real `ENOSPC`. The
one residual: a *genuinely* truncated input concurrent with a disk-full would report the reader's
error. Edge-of-edge, pre-existing, not worth changing.

### Can the five tests fail?

Audited §14's mutation table against the dependency source rather than trusting it.

- `into_bgzf_round_trips_a_small_batch` — **genuinely kills** the flush-omission bug. Proven from
  source: three ~50-byte records never reach `MAX_BUF_SIZE`, so with `flush()` removed `into_inner()`
  returns an **empty** `Vec` and `assert!(!blocks.is_empty())` fires. This is the mode that produces
  a structurally valid, `quickcheck`-clean BAM with a whole batch missing — the right test to get right.
- `into_bgzf_emits_no_terminator` — kills the `into_inner()` → scope-drop mutation. Distinct from
  the above, as §14 claims.
- `parallel_sink_terminates_the_file_exactly_once` — kills both an omitted and a doubled terminator
  (0 and 2 both fail the `assert_eq!(count, 1)`). Also the guard that stops anyone "tidying away"
  `let _ = term.into_inner();`.
- `write_error_is_returned_not_panicked` — fails loudly if it ever stops erroring (`expect_err`
  panics), so it cannot pass vacuously. Its margin rests on 2,000 records compressing to more than
  the ~3.9 KB of budget the header leaves; comfortable (~12 KB), but see L-6.
- `a_failed_run_publishes_nothing` — **the one weak test.** See M-2.

No vacuous passes found.

---

## Findings

### H-1 — Throughput was not measured, and the CHANGELOG records only the upside

**High (not a defect — your call).** `src/parallel.rs`, `CHANGELOG.md`.

Compression moves from a dedicated `cores`-wide pool onto the `cores` trim workers. On a machine
with more cores than `--cores`, the output stage no longer runs in parallel with trimming, and
level-6 BGZF over ~1.2 MB of raw BAM per 4,096-read batch is not a rounding error. Plan §7 predicted
a real regression and §9.5 set the bar: *"A regression up to ~25% at `--cores 8` is accepted per §7
**and recorded in the changelog**."*

It was neither measured (§14 deviation 3, disk full at 5.1 GiB free — correctly refused rather than
faked) nor recorded. The shipped entry advertises "29 threads to 20" and says nothing about the
trade. This is the one §9 acceptance criterion that is neither met nor waived in the artifact users read.

I did not measure it either — an honest benchmark needs an isolated target directory and a
before-binary, and this tree is shared with two other agents. Analytically the main thread does
*less* work than before (it writes ~45 MB compressed instead of memcpy-ing 67 MB uncompressed into
the pool's staging buffer), and on a box with exactly `--cores` cores this should be *faster*; the
exposure is a `--cores 8` run on a ≥16-core box.

**Pick one:** (a) free space and run §9.5 before committing; (b) commit now and add half a clause to
the entry — "…at a throughput cost not yet measured on hosts with spare cores"; (c) drop the thread
figures from the entry and let #568's docs follow-up carry the thread story. (b) is cheap and honest.

---

### M-1 — `Option<BgzfSink<..>>` and its three `expect`s are unnecessary

**Medium.** `src/bam.rs:1196-1246`. Not a live bug — "unreachable through the public API" is **true
today**: nothing calls `header()` or `write_compressed()` with `inner == None`, because only
`finish(self)` takes it and the value is gone afterwards. No `catch_unwind`, no `Clone`.

But §14 deviation 4's justification is wrong: *"`BgzfSink::finish` consumes itself to return `W`, so
the wrapper must move it out; the plan's §4 sketch implied a plain field."* `ParallelBamSink` has
**no `Drop` impl** (verified: no `impl Drop` anywhere in `src/bam.rs`), so a partial move out of
`self` is allowed and the plain field works:

```rust
pub struct ParallelBamSink {
    /// Declared before `pending` so the writer is flushed before the temporary is removed on the
    /// abandon path.
    inner: BgzfSink<BufWriter<File>>,
    pending: crate::io::PendingOutput,
}

pub fn finish(self) -> Result<()> {
    let Self { inner, pending } = self;
    inner.finish()?;   // on Err, `pending` drops here and removes the temporary
    pending.commit()
}
```

Both `Option`s go, all three `expect`s go, the `match` goes, and `header()`/`write_compressed()`
become one-liners. Drop order on the abandon path is unchanged — it is still the *field* order that
governs it, so keep `inner` first and keep the comment that says why.

Worth doing because a panic whose unreachability has to be argued in a comment is a standing
liability; this one can simply not exist.

---

### M-2 — `a_failed_run_publishes_nothing` passes with `PendingOutput::drop` deleted

**Medium.** `src/bam.rs:3139-3152`. The test creates a sink, drops it, and asserts
`!out.exists()`. But `PendingOutput::create` writes to `partial_output_name()` — `.unpublished.bam.partial`
in the same directory — so the final path was **never going to exist**. Remove the whole of
`impl Drop for PendingOutput` and this test still passes, leaving the temporary on disk. The
assertion it makes is true by construction, not by behaviour.

The repo's dominant failure mode is checks that pass because they never tested the thing, so this
one is worth a line:

```rust
assert!(
    !dir.path().join(".unpublished.bam.partial").exists(),
    "the abandoned temporary survived"
);
```

Mitigating context, which stops this being High: the pre-existing integration test
`ubam_out_cores_gt_1_failure_leaves_existing_output_untouched`
(`tests/integration_ubam_out.rs:369`) already drives the **real** `ParallelBamSink` to a mid-stream
failure — an over-long read name in batch 1, after batch 0 has gone through `write_compressed` —
and asserts a pre-existing file at the final name survives. So "nothing published" *is* covered
end-to-end; it is this unit test specifically that doesn't earn its keep. The narrow remaining gap
is a *write* failure (vs an encode failure) through the real sink, which M-2's sibling halves
(`write_error_is_returned_not_panicked` + the integration test) cover between them. I would not add
a test for the intersection.

---

### L-1 — The comment on the load-bearing line says the opposite of what happens

`src/bam.rs:1142`:

```rust
// `flush` stages the tail block; `into_inner` then takes the buffer, terminator unwritten.
```

`flush()` **writes** the tail block — `flush_block()` compresses the staging buffer and emits the
frame. `write()` is what stages. The second clause is right and is the important half. Suggest:

```rust
// `flush` emits the tail block; `into_inner` then takes the buffer, terminator unwritten.
```

One word, but this is the single line explaining the invariant the whole change rests on.

### L-2 — `bgzf_terminator`'s `into_inner()` reads as dead code

`src/bam.rs:1151` — `let _ = term.into_inner();` is what prevents `Drop` writing a *second*
terminator, and the plan calls it "load-bearing, not stylistic". The code gives no hint; a future
cleanup removes it. `parallel_sink_terminates_the_file_exactly_once` catches that, so this is
comfort not safety. Half a line, current-state only:

```rust
// `into_inner` disarms the `Drop` that would append a second terminator.
```

### L-3 — `BgzfBlocks::len()` is dead

`src/bam.rs:1121-1123` — unused in `src/` and `tests/`. It exists to pair with `is_empty()`, but
clippy's lint is `len_without_is_empty`; a lone `is_empty()` is fine. Delete `len()`.

### L-4 — "Failed to write BAM header" appears three times in one error chain

`BgzfSink::new` attaches it to both the compress and the write step, then `ParallelBamSink::create`
adds `"Failed to write BAM header to {path}"`. A real failure prints the same sentence three times
and the user cannot tell compression from I/O. Suggest `"Failed to compress the BAM header"` on the
`bgzf_members` call and dropping the inner `write_all` context (the outer one already names the path).

### L-5 — The new CI step covers single-end only

`--cores` defaults to 1, so the existing `Validate PE uBAM-out` step never touches
`ParallelBamSink`; the new step is SE. PE parallel output therefore has no htslib coverage in CI,
only the inflate-and-compare unit test at `src/parallel.rs:2857`. I verified PE correctness locally
(19,992 records, md5-identical, one terminator), so this is coverage not correctness. Cheapest fix:
add `--cores 4` to the existing PE step — it costs nothing and the flag-bit assertions already there
become parallel-path assertions.

### L-6 — `write_error_is_returned_not_panicked`'s comment overstates the budget

`src/bam.rs:3121` — *"Enough for the header member, nothing for the records."* There is ~3.9 KB
left after the ~200-byte header member; the test works because 2,000 records compress to more than
that, not because the budget is exhausted. Either set `budget` to the header size or say
"not enough for the records".

### L-7 — The diff invalidates open issue #568's thread table

#568 ("The documented N+4 thread model does not hold for uBAM input") tabulates
`uBAM in → uBAM out, SE or paired interleaved` at **3N+5**; this change makes it 2N+4, and
`src/cli.rs:357-362` still promises "an N+4 thread model applies". So in one release the CHANGELOG
says "29 threads to 20" while `--help` says 12 at `--cores 8`. Tracked and deliberate (plan §2
cites #568 as the docs follow-up), but plan artifacts rot silently — drop a one-line comment on #568
when this lands so its table isn't stale when someone picks it up.

### L-8 — Register of the CHANGELOG entry

I checked the house style before raising this, and the entry **fits** it: every neighbouring entry
is a bold headline plus one explanatory sentence over 3-4 lines, and #529 and #566 both carry
measurements. So this is not a violation of the one-sentence rule.

The one difference: neighbours name user-visible symptoms, this one names the dependency's
internals (`bgzf::MultithreadedWriter`, `Drop`, "inside its pool"). A reader deciding whether an
upgrade affects them doesn't need the mechanism — and "signalled failure by panicking inside its
pool" is slightly off: the pool's writer thread returns `Err` and exits, and the panic is the
*caller's* `send().unwrap()` on the next flush. Optional tighten:

> `bgzf::MultithreadedWriter` aborted the process without naming the `ENOSPC` or `EPIPE`; the
> workers now compress their own batches, which also drops nine threads at `--cores 8`.

### L-9 — Three-line comment on the new CI step

`.github/workflows/ci.yml:1280-1282`. Sentences two and three justify the step's existence, which
belongs in the commit message under the project's comment rule. The first sentence alone carries the
fact. Minor; the content is accurate and current-state.

---

## Not findings — checked and clear

- **Framing.** Header member, per-batch members, one terminator — walked the whole chain, clean.
- **Empty batch.** All-filtered batch → zero-length `BgzfBlocks` → zero-length `write_all`, no empty
  BGZF block emitted (`Writer::flush` early-returns on an empty staging buffer). Verified end-to-end.
- **Stale references.** `write_encoded` is gone from the entire repo; `MultithreadedWriter` is gone
  from `src/`. The three surviving "pool" comments are readers/demux and correctly out of scope.
- **`src/main.rs:2698`.** The third clause of the thread comment was dropped as §2 required.
- **Memory.** `pending: BTreeMap<_, BgzfBlocks>` now holds compressed batches — ~3× smaller than
  before. Per-worker transient peak rises slightly (raw + compressed alive together inside
  `into_bgzf`); net win.
- **Ordering.** Unchanged; the reorder buffer still delivers batches in sequence and the sink is a
  pure append. The 12,000-record integration test and both parity unit tests cover it.
- **CI step mechanics.** Lands in `validation-ubam` on `ubuntu-latest` with samtools installed and
  the release binary built; `set -e` plus the explicit `|| { exit 1; }` guards mean every comparison
  is live; a wrong output filename would fail loudly, not pass silently. 10,000 reads at
  `BATCH_SIZE = 4096` is genuinely 3 batches, so multi-member framing is exercised, not just
  header + 1. `samtools view` without `-h` excludes `@PG`, so the md5 won't shift on a version bump.
- **The 29-thread baseline.** Confirmed in #542's body ("reviewing #529 measured 29 threads at
  `--cores 8` single-end"). §14 labels it as not re-measured; the changelog doesn't carry
  provenance, which is normal for a changelog.

## Fixes applied

**None.** Two other reviewers were working this tree concurrently and the caller requires sight of
findings before any edit. No file under `src/`, `tests/`, `.github/` or the repo root was touched;
`git status` is unchanged apart from this report. All probes went to the scratchpad.

## Priority order

1. **H-1** — decide: measure §9.5, or add the caveat to the CHANGELOG. Half a clause either way.
2. **M-1** — drop the `Option` and the three `expect`s.
3. **M-2** — one assertion line so the test tests something.
4. **L-1** — one word, on the line that matters most.
5. L-2 … L-9 — take or leave.

None of 2-5 blocks a commit.
