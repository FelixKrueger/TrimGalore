# Code review — PR #529 `feat(ubam): multithreaded uBAM input and output with --cores`

**Reviewer A** · head `62f872fb` · base `dev` @ `c0bbee4` · 2026-11-10
Pinned tree: `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/aecc6671-429c-462c-ad9c-38f8e4d35fbe/scratchpad/pr529_A`
Built and executed at `target_revA` (real compile, 1m33s).

> No fixes applied. Recommendations only.

---

## Summary

The core parallelisation is **sound, and I confirmed it by running it.** Ordering,
phred-offset parity, header emission and the de-interleaver contract are all correct. On a
400,000-read uBAM-in → uBAM-out run — ~98 batches, so the reorder buffer is genuinely
exercised — `--cores 1` and `--cores 4` produce **byte-identical records in identical order**
and identical trimming reports. Throughput on that input scales 1.56s → 0.79s → 0.50s → 0.42s
across `--cores 1/2/4/8`, so the feature delivers what it claims.

What CI does not catch sits almost entirely in the **failure paths**, and there the parallel
path is materially weaker than the serial path it now replaces above `--cores 1`:

1. **H1 — `ParallelBamSink` bypasses `PendingOutput`, so a failed run publishes a truncated
   BAM at the final output name.** Verified: 400,000 reads with one bad record at 390,000 →
   `--cores 1` exits 1 and writes **no file**; `--cores 4` exits 1 and writes
   **`sample_trimmed.bam` containing 376,832 of 400,000 records**, with a byte-exact BGZF EOF
   block, which `samtools view -c` reads cleanly (exit 0). `samtools quickcheck` returns the
   same exit 8 for this file and for a known-good 10,000-record output, so it cannot tell them
   apart. This is the one finding I would gate the merge on.
2. **H2 — a write error mid-stream panics rather than returning.**
   `MultithreadedWriter::send` calls `write_tx.send(..).unwrap()`; once the writer thread has
   returned `Err` (ENOSPC, EDQUOT, EIO) that channel is disconnected and the next buffer flush
   panics. `--cores 1` returns a contextualised `anyhow` error.
3. **M1 — `total_bp_after_trim` varies with `--cores` on the uBAM path.** Verified:
   `--rrbs --output-format ubam` reports `400,000 bp (75.5%)` at `--cores 1` and
   `380,000 bp (71.7%)` at `--cores 4`. With `--clip_R1 5` the gap is wider:
   `400,000 (75.5%)` vs `350,000 (66.0%)`. The BAM records are bit-identical, so this is
   report-only — but `--clip_R1` is a routine flag, not a corner case.
4. **M2 — thread count is now up to ~4 × `--cores`** while `--help` still documents "N+4".
5. **M3 — three user-facing docs still say uBAM output is single-threaded,** including the
   published site under `docs/`.

No ordering, parity, determinism or data-corruption defect found. H1 is the one that matters.

---

## Critical

None.

---

## High

### H1 — `ParallelBamSink` bypasses `PendingOutput`: a failed run publishes a truncated BAM at the final output name

**`src/bam.rs:1083-1105`** (`ParallelBamSink::create`) vs **`src/bam.rs:903-923`**
(`BamWriter::create`):

```rust
// BamWriter::create — serial
let (pending, file) = crate::io::PendingOutput::create(path)?;   // writes to <path>.partial

// ParallelBamSink::create — parallel
let file = File::create(path)                                     // writes to <path>
```

`PendingOutput` (`src/io.rs:887-932`) writes to `partial_output_name(final_path)` and only
`rename`s into place on `commit()`; its `Drop` removes the temp when `commit` was never
reached. `ParallelBamSink` has neither — no temp, no `Drop` — and `finish()` is reached only
on success (`main.rs:2348-2350`, `2524-2526`, `2670-2677`: `run_*_to_bam(..)?` precedes
`sink.finish()?`, so an `Err` skips it).

What makes this worse than an ordinary leftover partial file: `bgzf::MultithreadedWriter`
implements `Drop` as `if !matches!(self.state, State::Done) { let _ = self.finish(); }`
(noodles-bgzf 0.35.0, `src/multithreaded_writer.rs`), and `finish()` flushes the buffer, joins
the deflaters and lets the writer thread append `BGZF_EOF`. The abandoned file is therefore
not merely truncated — it is a **structurally complete BAM**.

**Verified.** 400,000 reads, one 301-byte read name at record 390,000 (BAM's QNAME limit is
254 bytes — exactly the case `write_error_context` exists to annotate):

| | exit | published output |
|---|---|---|
| `--cores 1 --output-format ubam` | 1 | **no file** |
| `--cores 4 --output-format ubam` | 1 | **`sample_trimmed.bam`, 376,832 / 400,000 records** |

Reproduced 4/4 runs at 376,832 records (= 92 × 4,096; deterministic here because batch cost is
uniform, not because anything guarantees it). The published file:

- ends with the byte-exact 28-byte BGZF EOF block;
- `samtools view -c` → `376832`, **exit 0**;
- `samtools quickcheck -v` → exit 8 "had no targets in header" — **but a known-good
  10,000-record uBAM from the same binary returns the identical exit 8**, because unaligned
  BAM has no `@SQ` lines. So quickcheck cannot distinguish complete from truncated output
  here.

On a smaller input the leftover can instead be a **valid, empty BAM**: with the bad record at
9,500 of 10,000 the published file was 262 bytes, 0 records, exit 0 from `samtools view -c`,
6/6 runs. (The error short-circuits the merge loop at `src/parallel.rs:430-433` (paired) / `531-534` (single-end) — the
`Err` arm `break 'outer`s without draining already-completed earlier batches out of the
result channel, so a small run loses everything while a large one keeps a prefix.)

**Why it matters.** A pipeline that checks the exit code is fine. One that checks for the
output file's existence, or a user who re-runs and does not notice stderr scrolled past,
silently proceeds on 94% of their data with no signal from any standard BAM validity check.
The repo already treats this class as a defect: `plans/08132026_partial-output-on-refusal/`.

**Confidence: high — verified by execution.**

**Recommendation.** Give `ParallelBamSink` the same temp-then-commit discipline as
`BamWriter`: take `(pending, file)` in `create`, hold `pending: Option<PendingOutput>`
declared *after* `inner` so the BGZF stream closes before the temp is removed (the field-order
comment on `BamWriter` at `src/bam.rs:868` already records why that ordering matters), and
`commit()` in `finish()`. Secondarily, consider draining the already-completed batches before
honouring the error, so the failure mode does not depend on batch size.

### H2 — a mid-stream write error panics instead of returning an error

**`src/bam.rs:1111-1113`** (`ParallelBamSink::write_encoded`) → `MultithreadedWriter::write`
→ `MultithreadedWriter::send` (noodles-bgzf 0.35.0):

```rust
fn send(&mut self) -> io::Result<()> {
    ...
    write_tx.send(buffered_rx).unwrap();     // panics on a disconnected channel
    ...
    deflate_tx.send(message).unwrap();
}
```

`spawn_writer`'s body is `while let Ok(buffered_rx) = write_rx.recv() { ... write_frame(..)? }`
— the `?` ends the thread on the first I/O error, dropping `write_rx`. The next time
`write_encoded` fills a 64 KiB buffer, `write_tx.send(..).unwrap()` panics.

**Failure scenario.** `--cores 4 --output-format ubam` onto a filesystem that fills up, or an
NFS/quota'd path returning `EDQUOT`, part-way through a large run: a Rust panic and a non-zero
exit instead of `Failed to write BAM records: No space left on device`. At `--cores 1` the
same condition returns a clean contextualised error from `BamWriter::write_record`.

`ParallelBamSink::finish()` *does* handle this correctly — `inner.finish()` surfaces the
writer thread's `io::Result` and it is `.context(..)`-wrapped. The gap is only for errors
raised before `finish`, which for any output larger than one BGZF block is the common case.

**Confidence: high** that the panic path exists as described (read from the noodles source);
**medium** that a real ENOSPC reaches it — it needs at least one further buffer flush after
the writer thread dies, near-certain mid-run but not guaranteed on the final block. I did not
fill a filesystem to confirm.

**Recommendation.** This is a noodles API wart rather than a TrimGalore bug, but it is now on
a user-reachable path where `--cores 1` behaves better. Minimum: document the difference. The
proper fix is upstream — `MultithreadedWriter::send` returning `io::ErrorKind::BrokenPipe`
instead of unwrapping — and the project's standing preference is to send that patch rather
than file an issue.

---

## Medium

### M1 — `total_bp_after_trim` differs between `--cores 1` and `--cores N` on the uBAM path

**`src/parallel.rs:1435`** (`classify_single`) and **`src/parallel.rs:820-821`**
(`classify_paired`):

```rust
stats.total_bp_after_trim += record.seq.len();
```

**`src/trimmer.rs:672`** (`run_single_end_to_bam`) and **`src/trimmer.rs:780-781`**
(`run_paired_end_to_bam`):

```rust
stats.total_bp_after_trim += result.bp_after_cutadapt;
```

`bp_after_cutadapt` is captured at `src/trimmer.rs:165`, explicitly *before* RRBS truncation,
poly-G, poly-A, `trim_n` and the fixed `--clip_*` steps — the comment there and
`test_bp_after_cutadapt_skips_rrbs_truncation` (`src/trimmer.rs:1076`) pin that as deliberate
(#232: it has to match Perl v0.6.x on RRBS samples). `record.seq.len()` read after `trim_read`
returns is the *final* length. They coincide only when none of those later steps trims
anything.

**Verified**, 10,000 adapter-bearing reads, `--output-format ubam`, cutadapt-section
"Total written (filtered)" (`src/report.rs:635-636`):

| flags | `--cores 1` | `--cores 4` |
|---|---|---|
| `--rrbs` | 400,000 bp (75.5%) | 380,000 bp (71.7%) |
| `--clip_R1 5` | 400,000 bp (75.5%) | 350,000 bp (66.0%) |

The emitted BAM records are **bit-identical** between the two (`samtools view` output
compared with `cmp`), so this is report-only, not data corruption. But `--clip_R1` is a
routine flag — far more widely used than `--rrbs` — which makes the blast radius larger than
the RRBS framing suggests.

**Scope.** The divergence predates this PR: `classify_single`/`classify_paired` are untouched
by the diff, and I reproduced the identical 400,000 vs 380,000 split on the **FASTQ** path,
which has had `--cores N` all along. What this PR changes is that the uBAM path — previously
*always* serial, hence always `bp_after_cutadapt` — now takes the `parallel.rs` branch above
`--cores 1`. So it is newly reachable on a path it could not reach before, and it should be
fixed on both paths rather than only the new one.

**Why the new tests cannot see it.** `parity_test_config()` (`src/parallel.rs:1702-1729`) sets
`rrbs: false`, `poly_a: false`, `poly_g: false`, `trim_n: false`, and all four `clip_*` to
`None`. Under that config `bp_after_cutadapt == record.seq.len()` by construction
(`test_bp_after_cutadapt_equals_seq_len_when_no_rrbs`, `src/trimmer.rs:1102`), so
`assert_eq!(s1, p1)` in both new parity tests is satisfied for the wrong reason. This is the
"assertion points the same direction as the bug" shape the brief asked about — the assertions
are fine; the *fixture config* is what cannot discriminate. The pre-existing
`test_parallel_serial_trim_stats_parity` (`src/parallel.rs:1739`) shares the config and the
blind spot, which is presumably how the FASTQ-path divergence survived commit `3996fc5`, whose
message (quoted at `src/parallel.rs:1735`) claims to have fixed `total_bp_after_trim`.

**Confidence: high — verified by execution.**

**Recommendation.** Change `classify_single`/`classify_paired` to add
`result.bp_after_cutadapt`, and give at least one parity test a config with `rrbs: true` or a
non-`None` `clip_r1`. Fixing it changes FASTQ-path `--cores N` report output for RRBS and
clip users, so it wants its own CHANGELOG line and arguably its own PR — but it should not
ship *silently widened* to uBAM.

### M2 — thread count is now up to ~4 × `--cores`; `--help` still documents "N+4"

`src/cli.rs:357-363`, unchanged by this PR, still reads:

> From `--cores 2` upward, an N+4 thread model applies: N workers + 2 decompressors + 1
> batcher + 1 writer.

What `--cores N` now spawns for **uBAM in → uBAM out, paired, two files**
(`main.rs:2515-2527`):

| source | threads |
|---|---|
| `open_threaded_reader(input_r1, .., cli.cores)` → `open_threaded_pooled` | 1 read + N convert + 1 merge |
| `open_threaded_reader(input_r2, .., cli.cores)` | 1 read + N convert + 1 merge |
| `run_paired_end_parallel_to_bam(.., cores)` | N trim/encode workers + 1 reader |
| `ParallelBamSink::create(.., cli.cores)` | N deflaters + 1 bgzf writer |

≈ **4N + 6**. At `--cores 8` that is ~38 threads where the documented model promises 12. The
comment the PR adds at `main.rs:2658` ("N threads convert the input records, N workers trim
and encode, N threads compress") states the new fact correctly, but `--help` — what a user
sizing an HPC job actually reads — does not.

This matters for Slurm/LSF users who set `--cores` to `cpus-per-task`: the job oversubscribes
its allocation ~3×, which under a cgroup means contention rather than an error, so it presents
as a performance regression rather than a misconfiguration. The repo has
`plans/08182026_threading-docs-cluster-allocation` on exactly this.

Memory scales similarly: each conversion pool holds up to 2 raw batches per worker
(`sync_channel(2)`), plus `workers * 2` converted batches in the merge channel, plus the
reorder `BTreeMap`, plus `BAM_READER_CHANNEL_BATCHES = 4` on the consumer channel — at 4,096
records/batch and two readers for paired input, roughly linear in `--cores`.

My timing run (400k reads, uBAM→uBAM, single run, synthetic input — indicative, not a
benchmark) showed 1.56s / 0.79s / 0.50s / 0.42s at `--cores 1/2/4/8`, i.e. returns are already
flat past 4. That is weak evidence that N conversion threads *per reader* buys anything over
2 or 3.

**Confidence: high** on the thread counts (read from the call sites and
`spawn_conversion_pool`); the oversubscription *impact* is inference.

**Recommendation.** Update the `--cores` help text and threading docs to the real model; and
consider not multiplying the budget three ways — cap the conversion pool, or split `cores`
across the three stages rather than giving each the full count.

### M3 — three user-facing docs still say uBAM output is single-threaded

The PR updates `CLAUDE.md` and `README.md` but not:

- `docs/src/content/docs/guide/outputs.md:65` — "uBAM output is currently single-threaded —
  `--cores N` is ignored for BAM writing; FastQC (when `--fastqc` is requested) still uses the
  value."
- `wiki/usage/ubam.md:80` — "Threading | Single-threaded compression in v1. `--cores` is
  ignored for BAM writing **and an up-front note says so**; FastQC still uses it."
- `wiki/overview.md:73` — "uBAM output is single-threaded in v1 and refuses seven flags…"

`docs/` is the published site (auto-deploys from `dev`), so this ships a direct contradiction
of the binary — including a promise that users will see a startup note that this PR removes
(`base_c0bbee4/src/main.rs:575-583`).

**Confidence: high** (grep + read). **Recommendation.** Update all three in this PR.

### M4 — no CHANGELOG entry

The diff touches 8 files; `CHANGELOG.md` is not among them. `--cores` going from "silently
ignored, with a startup note" to "fully applied" on the uBAM path is precisely what a user
deciding whether to upgrade needs to see. One sentence, in the file's existing
`- [#529](url) - … ([@Justype](url))` format.

---

## Low

### L1 — `spawn_conversion_pool` can strand the merger if a worker panics

**`src/bam.rs:241-336` (reading thread at 272-280)**. The reading thread increments `seq` *before* it knows the send
succeeded:

```rust
let res = produce(&mut |raw: RB| {
    let idx = (seq as usize) % work_txs.len();
    seq += 1;
    work_txs[idx].send((seq - 1, raw)).is_ok()
});
...
let _ = merge_tx.send((seq, last));
```

If `convert` ever panics, that worker's `rx` drops, the batch it held never produces a merge
message, and the merger's `while let Some(m) = pending.remove(&expected)` never advances past
the missing sequence number. The merger then drains to channel-close and returns **without
ever calling `deliver`**, so the consumer's channel drops and `BamReader::next_record` hits
`Err(_) => Ok(None)` (`src/bam.rs:843`) — a panicked conversion surfaces as **clean EOF**,
silently truncating the input.

I could not find a panicking operation in `bam_record_to_fastq`: the `text[start..end]` slices
are constructed at `String::len()` boundaries, the `debug_assert_eq!` on tag length compiles
out in release, and every fallible step returns `Result`. So: **latent, not live.** Noting it
because it is the kind of structure a future edit to the conversion function turns into silent
data loss, with no signal when it happens.

**Recommendation.** Increment `seq` only on a successful send, and have the merger treat
channel-close with `expected` still pending as an error delivered to the consumer rather than
as EOF. The second half also hardens `BamReaderSource::Threaded`'s `Err(_) => Ok(None)`, the
shared mechanism behind both this and L2.

### L2 — `&&` short-circuit drops the R2 batch and ends the R2 side as clean EOF

**`src/bam.rs:459`** (`open_paired_interleaved_pooled`):

```rust
Ok(Some((r1s, r2s))) => tx_r1.send(Ok(r1s)).is_ok() && tx_r2.send(Ok(r2s)).is_ok(),
```

If the R1 consumer hangs up while the R2 consumer is still reading, `tx_r1.send` fails, `&&`
short-circuits, `r2s` is dropped un-sent, `deliver` returns `false`, the merger returns, and
`tx_r2` drops — so the R2 reader sees `Err(_) => Ok(None)`, a clean EOF mid-file.

Not reachable from the current callers: `read_pairs_round_robin` (`src/parallel.rs:955`) owns
both readers, reads them in lockstep, and drops them together on any error. The base's
equivalent (`base_c0bbee4/src/bam.rs:296-299`) uses
`if tx_r1.send(..).is_err() || tx_r2.send(..).is_err() { return; }`, which has the same
property — so **not a regression.** Noting it only because the new code preserves a shape whose
safety rests on a caller invariant that nothing states.

### L3 — all aux fields are now decoded, not just the requested ones

**`src/bam.rs:1432-1447`** (`bam_record_to_fastq`). The old code called `data.get(&tag)` per
requested tag; the new single-pass loop is:

```rust
for field in rec.data().iter() {
    let (tag, value) = field.context("failed to read an aux field from BAM record")?;
```

— so a record carrying a malformed or unrepresentable aux field the user did **not** request
now fails the whole run, where before it was skipped. Arguably more correct, but it is a
behaviour change on `--preserve-tags` that nothing documents or tests.

**Confidence: medium.** It turns on how much of a field's value noodles decodes eagerly in
`Data::iter` versus how much `Data::get` had to decode anyway to skip preceding fields; I did
not read `noodles-sam`'s `Data` far enough to settle it. Worth a line in the PR description.

I did check the two things here that *would* have been bugs, and both are fine: tag emission
order still follows the user's `--preserve-tags` order (the second loop iterates `tags`, not
the aux fields), and a duplicated request (`--preserve-tags CB,CB`) still emits the tag twice,
matching the old per-tag lookup.

### L4 — comment-style violations the PR introduces

Project rule: one line by default, current state only, never reference the state being fixed,
never justify the change in the source.

- **`tests/integration_ubam_out.rs:311-312`** — "`--cores N>1` + uBAM-out runs on the worker
  pool **(it used to warn that `--cores` was ignored)** and produces the output BAM." The
  parenthetical is textbook historical narrative. Cut to "…runs on the worker pool and
  produces the output BAM."
- **`src/bam.rs:1421-1425`** — the new four-line note in `bam_record_to_fastq`: "One pass over
  the record's aux fields (looking each tag up separately re-scans the fields once per
  requested tag, **which made this the dominant cost of reading a uBAM**). The first
  occurrence of a tag wins, as with `Data::get`." The bolded clause justifies the change in
  the source; it belongs in the commit message. The "first occurrence wins" sentence is the
  keeper — it is an invariant on the loop below it.

Separately (**not this PR's fault**, flagged per the brief): the 5-line doc comment on
`write_error_context` (`src/bam.rs:961-966`) is pre-existing text relocated by the refactor. It
reads as justification now that it sits above a function rather than inline, and it is a
reasonable thing to trim while those lines are already being touched.

### L5 — the integration-test assertion was inverted rather than retargeted

**`tests/integration_ubam_out.rs:328-332`**: `assert!(stderr.contains("ignored"))` became
`assert!(!stderr.contains("ignored"))`. A negative substring assertion on a common English
word is weak: it passes if the note is removed, but it also passes if the note is merely
reworded, and it fails spuriously if any future diagnostic uses the word "ignored" (several
already do, for unknown `--fastqc_args`). Asserting something positive about the parallel path
— e.g. that the record count matches a `--cores 1` run — would carry more.

---

## What I verified vs. what I suspect

**Verified by execution** (binary at `target_revA`, fixtures under
`.../scratchpad/exp_h1`, `exp_m1`, `exp_big`, `exp_parity`):

- H1's end state on disk, 4/4 runs on the large fixture and 6/6 on the small one, including
  the byte-exact BGZF EOF block, `samtools view -c` exit 0, and `samtools quickcheck`
  returning the same code for good and truncated output.
- M1's two report divergences, and that the BAM records are bit-identical across them, and
  that the same divergence reproduces on the FASTQ path.
- **Ordering and content parity end-to-end**: 400,000-record uBAM-in → uBAM-out at `--cores 1`
  vs `--cores 4` gave identical records in identical order (`samtools view` + `cmp`) across
  ~98 batches, with identical trimming reports.
- Small-fixture paired interleaved uBAM-in → uBAM-out parity (`ubam_paired_test.bam`, 20
  records): records and reports identical.
- Throughput scaling (one run each, synthetic input): 1.56 / 0.79 / 0.50 / 0.42 s at
  `--cores 1/2/4/8`.

**Verified by reading** (TrimGalore + noodles-bgzf 0.35.0 + noodles-bam 0.89.0):

- The ordered merge in `run_{single,paired}_end_parallel_to_bam` is correct: `BTreeMap` plus a
  monotone `expected`, stats merged *inside* the reorder loop so accumulation order is input
  order; zero-record batches write `&[]` harmlessly; the final partial batch is dispatched
  before the `None` sentinels.
- Phred parity holds. `config.phred_offset` and `BamWriter::create`'s `input_phred_offset`
  both come from `cli.phred_offset()` at one construction site (`main.rs:1351`).
  `build_record_buf` is the only place the subtraction happens, and both callers pass the
  input offset — nothing hardcodes 33 on the write side. The `const PHRED_OFFSET: u8 = 33` at
  `src/bam.rs:71` is the *read* side and correctly stays fixed. The new
  `batch_encoder_subtracts_input_phred_offset` test is genuinely discriminating: `hhhh` at
  offset 64 and `IIII` at offset 33 must both give 40, which a hardcoded 33 would fail.
- Header is written exactly once, in `ParallelBamSink::create`, before any record block, via
  the same `build_output_header` the serial writer uses — so the `@HD`/`@PG` chain and `ID`
  disambiguation are identical across paths. `BamWriter::create`'s `_preserve_tags` argument
  is unused, so `ParallelBamSink` not taking it is correct, not an omission.
- Compression level matches: `bam::io::Writer::new` → `bgzf::io::Writer::new` →
  `CompressionLevel::default()`; `multithreaded_writer::Builder::default()` → the same. BGZF
  block boundaries are driven by `MAX_BUF_SIZE` in `MultithreadedWriter::write`, i.e. by the
  byte stream alone and not by worker count, so output is deterministic across `--cores`.
- `read_raw_pairs` preserves the base de-interleaver's contract exactly — same
  exactly-one-of-FREAD1/FREAD2 check with the same message, same `> MAX_SLACK` bound applied
  after the pairing step, same orphan check after the final batch flush, same
  `GROUPED_INPUT_ERR` text. Compared line by line against `base_c0bbee4/src/bam.rs:162-345`.
- `--cores` defaults to 1 (`src/cli.rs:364`), so the default uBAM path stays on the serial
  writer and the CI golden fixtures are unaffected.
- Worker errors propagate, and `read_{single,pairs}_round_robin` return `Ok(())` when a send
  fails, so `reader_join_err.or(first_error)` does not mask the real worker error with a
  generic reader-side one.
- No deadlock on the happy path: in-flight batches are bounded by the per-worker
  `sync_channel(2)`, the `cores * 2` result channel and `BAM_READER_CHANNEL_BATCHES`; the
  reorder `BTreeMap` is bounded by those.
- Mate adjacency and FREAD1/FREAD2 survive by construction — `process_paired_batch_bam`
  pushes `Some(1)` then `Some(2)` into one batch's byte buffer, and batches concatenate in
  sequence order.
- Empty input (0 records) terminates cleanly on both pooled readers.
- `pooled_readers_report_the_same_errors` is discriminating: `assert_eq!(n, 4_096)` pins that
  the first batch is delivered *before* the error from the second, which is the ordering
  property that matters.

**Suspected, not executed:** H2's behaviour under a real ENOSPC (mechanism read, filesystem
not filled); M2's performance impact on a constrained node; L3's dependence on noodles' lazy
aux-decode semantics.

---

## Recommendations, in priority order

1. **H1** — give `ParallelBamSink` the `PendingOutput` temp-then-commit discipline.
   *Gate the merge on this.*
2. **H2** — handle, or at minimum document, the mid-stream write panic; prefer an upstream
   patch to `MultithreadedWriter::send`.
3. **M1** — switch `classify_single`/`classify_paired` to `result.bp_after_cutadapt`, and give
   one parity test a config with `rrbs: true` or a non-`None` `clip_r1`.
4. **M3** — update `docs/src/content/docs/guide/outputs.md`, `wiki/usage/ubam.md`,
   `wiki/overview.md`.
5. **M2** — fix the `--cores` help text; consider not multiplying the thread budget three ways.
6. **M4** — add the CHANGELOG line.
7. **L1 / L2** — make `BamReaderSource::Threaded` distinguish channel-close-while-pending from
   EOF; increment `seq` only on a successful send.
8. **L4 / L5** — trim the two comments; strengthen the inverted assertion.
