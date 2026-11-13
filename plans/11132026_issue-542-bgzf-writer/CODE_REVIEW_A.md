# CODE REVIEW A — #542, per-batch BGZF compression on the parallel uBAM sink

**Reviewer:** A (independent; B and a coverage auditor ran concurrently)
**Target:** uncommitted working tree, `dev` @ `d862f7d`
**Scope reviewed:** `src/bam.rs`, `src/parallel.rs`, `src/main.rs`, `.github/workflows/ci.yml`, `CHANGELOG.md` + docs copy
**Constraint honoured:** no file under `src/`, `tests/`, `.github/` or the repo root was edited. `git status --porcelain` at the end of the review is identical to the session-start snapshot. All probes live in `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/a2f653ef-8bfc-4b3c-9a4c-4b2b32c6aa15/scratchpad/`.

---

## Verdict

**Safe to commit.** No Critical or High finding. The BGZF framing is correct, and I verified it
against `noodles-bgzf-0.35.0`'s source rather than its documentation, and then against real htslib.
One Medium is a missing regression guard (not a defect in shipped behaviour); one Medium is a
misattributed measurement in the plan's own record. Everything else is Low and cosmetic.

The change is the rare refactor that is smaller *and* more correct than what it replaces: a panic
path became a `Result` path, a thread pool disappeared, and the object graph lost a layer.

---

## 1. Framing correctness — verified, sound

I read `noodles-bgzf-0.35.0/src/writer.rs` for `write`, `flush`, `flush_block`, `try_finish`,
`into_inner` and `Drop`, and every claim the implementation rests on holds:

| Claim the code depends on | Source evidence | Verdict |
| --- | --- | --- |
| `into_inner()` disarms the terminator-on-drop | `into_inner` is `self.inner.take()`; `Drop` is `if self.inner.is_some() { let _ = self.try_finish(); }` | ✅ |
| `flush()` on an empty staging buffer writes **nothing** | `fn flush` → `if self.staging_buf.is_empty() { Ok(()) }` | ✅ — this is what makes the empty-batch edge case safe |
| A fresh writer's `try_finish()` emits exactly the 28-byte EOF and nothing else | `try_finish` = `self.flush()?` (no-op) then `inner.write_all(BGZF_EOF)` | ✅ |
| `bgzf::Writer::flush` never flushes the *inner* writer | `flush` only calls `flush_block`; nothing touches `inner.flush()` | ✅ — hence the explicit `self.inner.flush()` in `BgzfSink::finish`, which is load-bearing |
| A batch > 64 KiB splits into several members automatically | `write` calls `self.flush()` once `!has_remaining()` | ✅ |

**Empty batch is the one I went looking for and it is clean.** Had `flush()` emitted an empty block
on an empty staging buffer, `bgzf_members(&[])` would have produced 28 bytes byte-identical to the
EOF marker, planted mid-file. It does not: an all-filtered batch yields a zero-length `Vec` and
`write_compressed` writes nothing. `src/bam.rs:1138-1144`.

**End-to-end, against real htslib** (local samtools 1.x, the shipped
`target/release/trim_galore`, md5 `b81af3fa9966fb5f17b4c73bddb395b9` — i.e. the binary the plan's
§14 measurements were taken on, so no stale-binary risk):

```
records: --cores 1 = 9996, --cores 4 = 9996
record-body md5:   24fb8d5f03ae08e2abbe8b007356b02d (identical)
--cores 4: terminator count 1, at end of file
--cores 1: terminator count 2 (pre-existing serial behaviour, out of scope)
```

10 000 reads at 4096/batch is 3 batches, so the ordered multi-member append is genuinely exercised,
not just a single-member degenerate case.

**Header member.** `BgzfSink::new` writes the header as its own member. This is strictly *more*
conservative than the old behaviour (where the header shared a block with the first records), and
`bam::io::Reader::read_header` plus `samtools view` both accept it. No finding.

---

## 2. Error and abandon paths — verified, sound

Traced create → write → finish → commit, and every early return, against `crate::io::PendingOutput`
(`src/io.rs:878-932`).

- `ParallelBamSink::create`: if `BgzfSink::new` fails, `?` drops `pending` → `Drop` removes the
  temporary. The final path was never created. ✅
- `write_compressed` failure → `parallel.rs` sets `first_error`, `break 'outer` → propagates out of
  `run_*_parallel_to_bam` → `main.rs`'s `?` drops the sink → `inner` (the `BufWriter<File>`) drops
  first, `pending` second, temporary removed. ✅ Field order is correct and the comment on it
  states the current fact.
- `finish`: `inner.finish()?` consumes the `BgzfSink`, so the `File` is closed **before**
  `pending.commit()` renames. On the error path `self.pending` is still `Some` and its `Drop`
  removes the temporary. ✅ This is the H-2 fix and it is structurally right.
- A `try_finish` that fails part-way through writing the EOF block leaves `term.inner` as `Some`,
  so `Drop` retries the write and corrupts the tail — but the error has already propagated and the
  temporary is deleted, so nothing is published. Checked; **not** a finding.

`run_ubam_output_{single,paired_two_files,paired_single_file}` all use `sink.finish()?`, so no
finalisation error is swallowed at any call site.

---

## 3. The `Option<BgzfSink<..>>` and its three `expect`s — "unreachable" is true

`inner` becomes `None` only inside `finish(mut self)`, which consumes the receiver. After the
`take()` the only remaining use of `self` is `self.pending.take()`; no `&self` or `&mut self`
escapes, and the type is not `Clone`. **The three `expect`s are genuinely unreachable through the
public API.** See L7 for a way to delete them anyway.

---

## 4. Can the five new tests fail? — audited

I re-derived each mutation from the noodles source rather than taking §14's table on trust. All five
pass locally (`cargo test --lib`, 548 filtered siblings), and four of them die on a real mutation:

| Test | Mutation it kills | My independent check |
| --- | --- | --- |
| `into_bgzf_round_trips_a_small_batch` | `flush()` removed | ✅ 3 records never fill `MAX_BUF_SIZE` (65 498 B), so `into_inner()` returns an empty `Vec` and `assert!(!blocks.is_empty())` fires. The positive round-trip is the right shape; the plan was right that a negated `ends_with` would have passed here. |
| `into_bgzf_emits_no_terminator` | `into_inner()` → scope drop (#542's own sketch) | ✅ `Drop` → `try_finish` → `BGZF_EOF` appended; the window search finds it. And the round-trip test *would* survive this mutation (an EOF block inflates to nothing), so §14's claim that the two cover distinct modes is correct. |
| `parallel_sink_terminates_the_file_exactly_once` | `bgzf_terminator` call omitted (count 0) **and**, unlisted in §14, deletion of `let _ = term.into_inner()` (count 2) | ✅ kills both directions |
| `write_error_is_returned_not_panicked` | write error swallowed | ✅ but see L2 on its margin |
| `a_failed_run_publishes_nothing` | only a `create` that writes to the final path directly | ⚠️ weakest of the five — see L9 |

No test passes vacuously. The one gap is a path no test reaches at all: **M1**.

---

## 5. The new CI step — it runs, and it fails on a regression

The step lands in the `validation-ubam` job (`.github/workflows/ci.yml:1138`), which installs
`samtools` (`:1150`) and builds the release binary before it. I ran the step's whole body locally
against the shipped binary and it behaves as written:

- Output naming is right: `in.fq` → `in_trimmed.bam` via
  `naming::single_end_bam_output_name` → `strip_fastq_extensions` (`src/io.rs:312-320`). Confirmed
  on disk, both `--cores` values.
- `--cores 1` takes the serial `BamWriter` (`cli.cores > 1` gates the parallel branch), so the
  comparison is genuinely serial-vs-parallel, and `samtools view` (no `-h`) correctly excludes the
  `@PG CL:` line, which differs between the two runs.
- Portability: `md5sum | cut -d' ' -f1`, `gzip -dc`, `python3` are all present on `ubuntu-latest`;
  the default `bash -e` shell means the `python3 -c` assertion failing fails the step.
- The `python3 -c` one-liner is correct — in `-c` mode `sys.argv[1]` is the file. §14's iteration
  log #4 (a heredoc would have shipped ten-space-indented Python) is a good catch and the
  replacement is sound.
- Dropping `samtools quickcheck` is right: it exits 8 with "had no targets in header" on any
  unaligned BAM.

One honest limitation worth recording somewhere: the step proves htslib **accepts** the framing and
that records survive, but like the unit tests it cannot prove a mid-stream terminator is harmful,
because htslib skips empty blocks rather than stopping at them. That is fine — `into_bgzf_emits_no_terminator`
is the guard for that, and the plan already ranks it as cosmetic.

---

## Findings

### Medium

**M1 — `finish()`'s error path has no test, so half the fix is unguarded.**
`src/bam.rs:1185-1192`. The change claims to retire three defects; the write-error path (C-1) now
has `write_error_is_returned_not_panicked`, but the finalisation path (H-2 — "`Drop` discards the
finish error, leaving a truncated file that looks structurally valid") has none. A mutation that
swallowed it — `let _ = bgzf_terminator(&mut self.inner);` or `let _ = self.inner.flush();` — passes
all five new tests, all the integration tests and the new CI step, because every one of them
finalises onto a writer that succeeds. The `BgzfSink<W>` split exists precisely to make this
constructible, and the test is four lines:

```rust
let mut sink = BgzfSink::new(FailingWriter { budget: 300 }, header.clone())?;   // header fits
let err = sink.finish().expect_err("a failing writer must fail the terminator");
assert!(format!("{err:#}").contains("simulated device failure"));
```

Same reasoning covers `BgzfSink::new`'s own header-write error (`budget: 0` → `new` must return
`Err`), which is likewise untested. Neither is a defect in the code as written; both are the
regression guard the plan's own §12 standard ("a fix for a write-error bug that cannot be shown to
handle a write error is not finished work") asks for, applied to the second error path.

**M2 — PLAN §14's §9.2a row mixes two different runs' figures.**
It reads: "`samtools view -c` decodes **496,150** records from a 43 MB parallel-written uBAM, and
the record-body md5 is identical at `--cores 1` and `--cores 4` (`24fb8d5f…`)". The md5 of the
record bodies of the **9,996**-record, 218 KB `BS-seq_10K_R1` output — the fixture the new CI step
uses — is `24fb8d5f03ae08e2abbe8b007356b02d`, which I reproduced independently above. Two different
record sets cannot share an md5, so the count/size and the md5 in that row come from different runs
(the 496,150 / 43 MB figures match §9.4's 500 k-read, 44,869,908-byte input). The *code* is fine and
both underlying runs evidently happened; the written record is wrong, and given this repo's history
of figures that survive review rounds, it should be corrected before the plan is archived. Either
split the row in two or re-state it against the CI fixture.

### Low

**L1 — a write failure still does not name the output file.** `write_compressed` →
`"failed to write BAM records"`, `finish` → `"failed to finalise uBAM output (BGZF EOF marker)"` /
`"failed to flush uBAM output"`. None names a path, and no enclosing frame in `main.rs` adds one
(`create` does — `"Failed to write BAM header to {path}"`). #542's complaint is literally that the
abort happened "without naming the `ENOSPC` or `EPIPE`"; the errno now surfaces but the file does
not. `ParallelBamSink` holds the `PendingOutput`, which exposes `tmp_path()` for exactly this.
Mitigating: `eprintln!("Output:   {}", output_path.display())` runs a few lines earlier, so the user
can infer it.

**L2 — `write_error_is_returned_not_panicked` rests on an unstated size inequality.**
`src/bam.rs:3113-3134`. `FailingWriter { budget: 4096 }` must cover the header member and *not*
cover the records. Reconstructing the 2 000 records' BAM encoding and deflating at level 6
(`scratchpad/revA_bgzf_margin.py`) gives ≈5.4 KB of members against ≈3.9 KB of remaining budget — a
margin of about 1.4x, on highly compressible synthetic data whose ratio depends on the flate2
backend. If that inequality ever flips the test **fails loudly** rather than passing vacuously, so
this is fragility, not a false negative. A `budget` just above the header member, or a
`fail_after_header: bool`, removes the coupling. The dual-purpose `budget` is also what makes the
comment `// Enough for the header member, nothing for the records` inaccurate — there are ~3.9 KB
left for records.

**L3 — `BgzfBlocks::len()` is dead.** `src/bam.rs:1122-1124`. No caller in `src/`, `tests/` or
`examples/`; only `is_empty()` (tests) and `as_bytes()` (the sink and the tests) are used. `pub` on
a library crate so neither `dead_code` nor clippy flags it. Delete it, or keep it only if it is
wanted as API surface.

**L4 — the comment on `bgzf_members`' `flush()` inverts noodles' terminology.**
`src/bam.rs:1141`: "`flush` stages the tail block". In `noodles-bgzf`, `write` *stages* (into
`staging_buf`) and `flush` *emits* the staged bytes as a block. Since the whole correctness of this
helper turns on that one call, the comment should say what it does: "`flush` emits the staged tail
block; `into_inner` takes the buffer with no terminator."

**L5 — `let _ = term.into_inner();` is load-bearing and unannotated.** `src/bam.rs:1151`. It is the
line that disarms `Drop`'s `try_finish`; delete it and every file gets two terminators. A reader
cleaning up "a discarded return value" has no way to know. A one-liner stating the current fact
("`into_inner` leaves `Drop` nothing to finish, so the block is written once") is within the comment
rules — it describes present behaviour, not the state being fixed. The behaviour is test-covered,
so this is about the next maintainer, not about correctness.

**L6 — `BgzfSink::new` labels a compression failure as a write failure, and the caller repeats the
sentence.** `src/bam.rs:1170-1174` uses `"Failed to write BAM header"` for *both* `bgzf_members`
(compression) and `write_all`, and `ParallelBamSink::create` wraps the whole thing in
`"Failed to write BAM header to {path}"`. A real failure prints
`Failed to write BAM header to X: Failed to write BAM header: <io>`. Give the compression arm its
own wording (`"Failed to compress the BAM header"`) and drop the inner duplicate.

**L7 — the `Option` wrapper is avoidable.** `ParallelBamSink` has no `Drop` impl (`grep "impl Drop"
src/bam.rs` → nothing), so `finish(mut self)` can move the field out directly —
`let inner = self.inner; inner.finish()?; match self.pending.take() { … }` — which deletes the
`Option`, all three `expect`s and the comment explaining them, with the same drop ordering on both
the success and the error path. Worth one compile to confirm. The counter-argument for keeping the
consuming `BgzfSink::finish(self) -> Result<W>` is that it makes a double-finish (two terminators) a
compile error, which is a real guarantee; a `finish(&mut self) -> Result<()>` would also remove the
`Option` but give that up. Either way the current shape is correct, just heavier than it needs to be.

**L8 — the CHANGELOG entry is the longest in its section, narrates the replaced mechanism, carries
an un-remeasured figure, and is arguably in the wrong section.** Four small things, one entry:
(a) five lines against a section whose entries run three to four, and the clause
"`bgzf::MultithreadedWriter` signalled failure by panicking inside its pool, whose `Drop` then
panicked again" is the previous state's internals, which the repo's own rule says belongs in the
commit message — note that §14 deviation 6 records catching exactly this in a source comment, so the
corrected artifact reproduced the defect it corrected; (b) "29 threads" is #542's figure, never
re-measured (§9.6 says so plainly), while the changelog states it as fact — either re-measure
against a base build or write "about a third fewer threads"; (c) it links an issue where every
sibling links a PR; (d) a panic-to-error fix reads as a bug fix, and the closest analogue in the
file — #565, "A reader thread that dies is no longer read as end of input" — sits under
`#### Bug fixes` two sections down. The thread-count half is a Change, so splitting the entry would
resolve both. I checked the house style before raising this: the bold-headline-plus-one-sentence
shape *is* the local convention, so only the length, the figure and the placement are at issue, not
the form.

**L9 — `a_failed_run_publishes_nothing` is near-redundant and checks less than its name.**
`src/bam.rs:3137-3151`. It asserts the final path is absent, which is also true if `Drop` never
removed the temporary; it does not assert the `.unpublished.bam.partial` is gone. `PendingOutput`'s
own tests (`src/io.rs:1147`) cover the temporary's lifecycle, and
`tests/integration_ubam_out.rs::ubam_out_cores_gt_1_failure_leaves_existing_output_untouched` covers
the same property end-to-end on a real failure. Keep it (it is nearly free) but assert the temporary
is gone too, which is the part nothing else checks at this layer.

**L10 — the CI step's comments run to seven lines across three blocks.** Against the repo's
one-line-default rule for committed config. The three-line block above `run:` and the three-line
block above `samtools view -c` both carry reasoning that belongs in the commit message; the
`quickcheck` explanation is the one worth keeping, as it pre-empts a reviewer adding it back.

**L11 — the library API break is unrecorded.** `ParallelBamSink::create` loses a parameter and
`write_encoded` is renamed to `write_compressed`; both are `pub` on a crate published to crates.io.
§8 of the plan acknowledges it and there are no in-repo external callers, but the changelog says
nothing. One clause on the existing entry would cover it.

---

## Checked and found sound (no action)

- Batch ordering: the reorder buffer is untouched; a batch that errors never lets its successors
  reach the sink, because `expected` cannot advance past it.
- `BgzfBlocks: Send` (newtype over `Vec<u8>`) — required by the `mpsc` hand-off inside
  `thread::scope`; auto-derived.
- Memory: `pending` now holds compressed rather than raw batches, and the pool's bounded channels of
  `BytesMut` are gone. Strictly lower.
- `bam::io::Writer<Vec<u8>>` has no `Drop`, so `head.into_inner()` yields exactly the raw header
  bytes — unchanged from the previous implementation.
- `into_bgzf`'s `Result` is effectively infallible (the sink is a `Vec`), which is fine; the
  signature keeps the boundary honest.
- `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean
  (ran after waiting on another agent's build lock — a warm cache, but cargo re-lints changed files,
  and the tree is unchanged since).
- The seven doc-comment rewrites listed in §2 are all present and all state current behaviour. The
  one I would still question is `BgzfSink`'s "Generic over the writer so a test can inject one that
  fails" — it justifies a design choice in source — but it documents an invariant a maintainer needs
  (don't collapse the generic), so I would leave it.

---

## Recommended order

1. **M2** — correct the §9.2a row (documentation integrity; costs nothing).
2. **M1** — add the two four-line error-path tests for `BgzfSink::finish` and `BgzfSink::new`.
3. **L8** — split or trim the changelog entry and settle the section.
4. **L1, L4, L5, L6, L3** — a single tidy-up commit.
5. **L2, L7, L9, L10, L11** — optional; none blocks the commit.

None of these changes the framing, the error semantics or the output, so the commit need not wait on
them if the plan prefers to land the fix and follow up.
