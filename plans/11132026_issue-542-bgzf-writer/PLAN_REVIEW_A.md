# PLAN REVIEW A — #542, replace `bgzf::MultithreadedWriter` on the parallel uBAM output path

**Reviewer:** A (independent; no contact with Reviewer B)
**Plan:** `/Users/fkrueger/Github/TrimGalore/plans/11132026_issue-542-bgzf-writer/PLAN.md`
**Tree:** `dev` @ `d862f7d`, working tree clean (`git diff --stat` empty before and after this review)
**Probes:** `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/a2f653ef-8bfc-4b3c-9a4c-4b2b32c6aa15/scratchpad/bgzfprobe/`
(standalone crate against `noodles-bgzf =0.35.0` / `noodles-bam =0.73.0` / `noodles-sam =0.69.0`, the
exact versions in `Cargo.lock`; nothing in the repo was modified)

---

## Verdict

**Safe to implement as written? No — not until §9 is strengthened. The design is right; the validation is not.**

The replacement design is correct and the three defects it retires are real. Every one of §3's four
facts holds under re-derivation. But the plan's own risk ranking inside §3.1 is **inverted**, and as a
direct consequence the §9.1 test — the one the plan calls "the single easiest thing to get wrong in the
whole change, and nothing else would catch it" — is, as specified, **satisfied by the more dangerous of
the two failure modes it exists to catch**. Separately, §9.3 (the regression guard for C-1, the defect
this issue exists to fix) is **not constructible as specified**, so C-1's fix would ship unguarded.

Fix Critical-1 and Important-1/2 and this is a clean, well-scoped change.

---

## 1. §3 re-derived, fact by fact

I read the pinned source rather than accepting the plan, and wrote probes where a probe would settle it.

### §3.1 — "the issue's sketch appends an EOF block to every batch" → **HOLDS**

Source confirms it exactly. `noodles-bgzf-0.35.0/src/writer.rs`:

- `Drop for Writer`: `if self.inner.is_some() { let _ = self.try_finish(); }`
- `try_finish`: `self.flush()?` then `inner.write_all(BGZF_EOF)`
- `into_inner`: `self.inner.take().unwrap()` — leaves `None`, which disarms that `Drop`
- `Write::flush`: `if self.staging_buf.is_empty() { Ok(()) } else { self.flush_block() }`

Probe (`src/main.rs` in the probe crate), 21-byte payload:

| Shape | Bytes | EOF blocks |
| --- | --- | --- |
| scoped `write_all` + `flush`, dropped at scope end | 77 | 1 (at end) |
| `write_all` + `flush` + `into_inner()` | 49 | 0 |
| `write_all` + `into_inner()`, **no flush** | **0** | 0 |
| empty writer, `try_finish()` + `into_inner()` | 28 | 1 (is exactly the EOF block) |
| empty writer, `try_finish()` then **dropped** | **56** | **2** |
| empty writer, `flush()` + `into_inner()` | 0 | — |

Two notes on the plan's table. Its figures are 66/38/28; mine are 77/49/28. The difference is payload
length only — the structural conclusions are identical, and **only the 28 is a constant**. Don't let the
66/38 figures reach a test assertion.

And `into_inner()` is load-bearing in *two* places, not one — see Important-1.

### §3.2 — "the terminator needs no private constant" → **HOLDS**

`BGZF_EOF` is `pub(crate) static` at `writer.rs:28`. A fresh writer over the output with an empty
staging buffer and a bare `try_finish()` emits exactly those 28 bytes (probe row 4) because `flush()`
short-circuits on an empty staging buffer. The plan is right to refuse a copied constant.

### §3.3 — "compression level is already equal on both sides" → **HOLDS**

Settled by reading, as the plan says:

- `writer/compression_level.rs`: `impl Default for CompressionLevel { fn default() -> Self { Self(6) } }`
- `writer/builder.rs`: `Builder { compression_level: CompressionLevel }` with `#[derive(Default)]`
- `multithreaded_writer/builder.rs`: `impl Default` sets `compression_level: CompressionLevel::default()`
- Both paths funnel into the *same* `deflate::encode(src, compression_level, dst)` with the same
  `CompressionLevelImpl` type alias — the `libdeflate` feature is crate-wide, so it cannot differ
  between the two writers even if it were enabled.
- `src/bam.rs:1128-1130` sets only `set_worker_count`.

Level 6 before and after, on all three paths. Nothing to measure. The issue's request to "compare
against a serial run" is correctly dismissed.

### §3.4 — "no byte-identity invariant exists on the `.bam` container" → **HOLDS, and is understated**

Confirmed and then some:

- `assert_ubam_eq` (`tests/integration_ubam_out.rs:125`) compares header-minus-`@PG` plus decoded record
  tuples. ✓
- CI's clump-only determinism check is `samtools view … | sort | md5sum` (`ci.yml:1082-1083`). ✓
- The Perl matrix is FASTQ-only at `--cores 1`. ✓
- The two serial-vs-parallel unit tests (`src/parallel.rs:2881-2890`, `:2947-2955`) compare the
  **decompressed** stream byte-for-byte, with the comment "the BGZF framing may differ; the
  decompressed BAM stream must not". That is exactly the right invariant and it survives this change
  untouched — block boundaries move, the uncompressed byte stream does not.

What the plan misses, in its favour on safety but against it on coverage: **no CI job exercises
`ParallelBamSink` at all.** The uBAM-out validation block (`ci.yml:1259-1380`) runs
`--output-format ubam` with no `--cores`, i.e. the serial `BamWriter`. The only `--cores 2` uBAM CI
steps are the four `--clump_only` ones (`ci.yml:1038, 1058, 1078, 1095`), and §2 is correct that
`clump_only` uses the serial writer (`src/clump_only.rs:698, :943`). So the guard is not "the §9 tests
plus a FASTQ-blind matrix" — it is three `cargo test` tests and nothing else. See Important-3.

### Scope claims in §2 — all check out

`ParallelBamSink` has exactly six references repo-wide (3 × `main.rs`, 2 × `parallel.rs` tests,
1 × `bam.rs` test) — the plan's call-site list is complete. All three dispatch sites are gated on
`if cli.cores > 1` (`main.rs:2379, :2553, :2699`), so §10's "`--cores 1` never reaches this path" holds.
There is no `impl Drop for ParallelBamSink` to remove. The `main.rs:2700` comment the plan wants
rewritten reads exactly as quoted.

---

## 2. Critical

### C-1 — §9.1 as specified is **passed** by the worse of the two failure modes

This is the one thing that must change before implementation.

§3.1's hazard has two sides, and they are not equally dangerous. The plan treats the EOF-appended side
as the catastrophe and never names the other side as a testable failure at all. Measured, with
`samtools` 1.x as an independent reader:

| Variant | samtools `view -c` | Exit | `quickcheck` | Size cost |
| --- | --- | --- | --- | --- |
| Correct (members + one EOF) | 10 / 10 records | 0 | pass | — |
| **Per-batch EOF** (the #542 sketch) | **10 / 10 records** | 0 | pass | +28 B per batch |
| **Omitted flush, small batch** | **0 records** | **0** | **pass** | file is header + EOF |
| Omitted flush, 4096-record batch | `error reading file` | 1 | — | tail block lost |

So:

- **Per-batch EOF is nearly harmless.** samtools reads every record; noodles' reader reads straight
  through a mid-stream EOF block and returns an identical payload. The cost is 28 bytes per batch —
  against ~433 KB of encoded records per 4096-read batch, that is 0.006%. The plan's "would have shipped
  as a silently bloated file" overstates it by about four orders of magnitude.
- **Omitted flush is the catastrophe**, and only on small inputs, which is what makes it vicious. Below
  ~64 KiB per batch nothing has flushed yet, so `into_inner()` returns an **empty Vec**: a structurally
  valid, `quickcheck`-clean BAM with *every record silently gone*, exit 0. Above 64 KiB the full blocks
  have already gone out and only the tail is lost, which corrupts the record stream and fails loudly.

Now read §9.1's assertion against that: *"assert the bytes do **not** end with the 28-byte EOF block."*
An empty `Vec` does not end with the EOF block. **The broken implementation passes.** And the end-to-end
check in the same row ("exactly one appears in a finished file") also passes, because a file of
header + EOF has exactly one EOF block at the end.

The plan says of §9.1: "this is the single easiest thing to get wrong in the whole change, and nothing
else would catch it." Both halves of that sentence are wrong in the same direction. The flush omission
*is* the easiest thing to get wrong — the issue itself warns about it twice — and the existing
decompressed-stream comparisons at `src/parallel.rs:2881`/`:2947` *would* catch it (their fixtures are
10,000 and ~12,000 records, so the loud variant), whereas they would **not** catch per-batch EOF (an EOF
block decompresses to zero bytes, so the decompressed streams stay identical). The plan's reasoning for
why §9.1 is needed is sound; its reasoning about which bug it catches is backwards.

**Required:** §9.1 must assert, over a *small* batch (a few records, so the staging buffer never fills
on its own):

1. `!blocks.is_empty()` — kills the flush omission, which no other assertion in the plan or the repo does;
2. inflating `blocks` yields **exactly** `BatchEncoder::into_bytes()` for the same records — the positive
   statement, which subsumes (1) and also catches a short write;
3. **zero** EOF blocks anywhere in the batch bytes (count occurrences, not just the suffix);
4. exactly one EOF block in a finished file, at the end — keep this, it is the assertion that earns its
   keep (it is the only guard on the terminator path; see Important-1).

For (3)/(4), derive the 28 expected bytes in the test by running `try_finish()` on an empty
`bgzf::Writer` over a `Vec` — the §3.2 trick — rather than hardcoding them. §3.2 argues against a second
source of truth for the constant in production code; the same argument applies to the test, and the
trick is two lines.

---

## 3. Important

### I-1 — the terminator is a **second** place the `into_inner` discipline can be wrong, and it sits outside the shared helper

§6 step 1's stated justification for `bgzf_members` is "there is exactly one place that can omit the
flush". That is not true as designed: `finish()` uses `try_finish()` + `into_inner()`, which does **not**
go through `bgzf_members` (the helper's shape, `write_all` → `flush` → `into_inner`, produces nothing for
an empty input). So the discipline lives in two places, and the helper's whole rationale is weakened.

Getting the terminator wrong is **invisible**: `try_finish()` followed by a plain drop emits two EOF
blocks (probe row 5, 56 bytes), and samtools reads such a file without complaint. Worse — and this is
worth knowing before anyone writes the test — **the shipping serial `BamWriter` already does exactly
this.** `BamWriter::finish` calls `inner.try_finish()` then `drop(self)`, and the `bgzf::Writer` inside
still holds `Some(inner)`, so its `Drop` runs `try_finish` a second time. Confirmed on real output from
the tree's binary:

```
$ xxd -s -56 ubam_test_trimmed.bam | tail -4
0000028c: 1f8b 0804 0000 0000 00ff 0600 4243 0200  ............BC..
0000029c: 1b00 0300 0000 0000 0000 0000 1f8b 0804  ................
000002ac: 0000 0000 00ff 0600 4243 0200 1b00 0300  ........BC......
000002bc: 0000 0000 0000 0000                      ........
```

Two back-to-back 28-byte EOF blocks. Three consequences for this plan:

- §2's "the serial `BamWriter` is unaffected … already finalises through a fallible `try_finish`" is true
  about *fallibility* but positions the serial writer as a reference to copy. On this detail it is not one.
- If the implementer factors §9.1's "exactly one EOF block" assertion into a helper and applies it to both
  sides of the §9.2 serial-vs-parallel comparison — the natural thing to do — **it will fail on the serial
  side** and send them chasing a pre-existing, harmless quirk. Say so in the plan.
- It independently reinforces §3.4: the serial path (2 EOF blocks) and the current parallel path (1, because
  `MultithreadedWriter::finish` sets `State::Done` so its `Drop` is a no-op) *already* differ in container
  bytes, and nothing has ever noticed.

**Suggested:** make both paths pure functions returning bytes, so the sink only ever does `write_all`:

```rust
fn bgzf_members(raw: &[u8]) -> io::Result<Vec<u8>>   // write_all -> flush -> into_inner
fn bgzf_eof() -> io::Result<Vec<u8>>                 // try_finish -> into_inner, 28 bytes
```

That restores §6 step 1's "exactly one place" property, makes the terminator unit-testable on its own
(`assert_eq!(bgzf_eof()?.len(), 28)`), and leaves `ParallelBamSink` with no BGZF knowledge at all.

### I-2 — §9.3 is not constructible as specified, so C-1's fix ships with no regression guard

Two separate blockers:

1. **The cited precedent is unreachable.** `unwritable_file` is at `src/fastq.rs:814` (not `:836` — that
   line is inside a test that *uses* it), inside the `#[cfg(test)] mod tests` that starts at `:807`. An
   integration test in `tests/integration_ubam_out.rs` cannot see it.
2. **More fundamentally, there is no injection point.** `unwritable_file` returns a read-only `File`, and
   `ParallelBamSink::create` takes a **path**, not a writer — `PendingOutput::create` makes its own handle.
   A read-only *directory* makes `PendingOutput::create` fail, which tests `create`, not `write_compressed`.
   Short of a full filesystem or Linux's `/dev/full`, there is no portable way to make a write to the
   `.partial` file fail.

The structural lesson is right there in the precedent: `src/fastq.rs`'s `Sink` is testable **because it is
constructed from a writer**. `ParallelBamSink` is not.

This matters more than a missing test usually would. The existing
`ubam_out_cores_gt_1_failure_leaves_existing_output_untouched`
(`tests/integration_ubam_out.rs:368`) induces an **encoder** error (a 300-byte read name), not an I/O
error — so it passes today *with* C-1 present. There is currently no test anywhere that exercises a write
failure on this path, and §9.3 as written cannot become one. C-1 is the defect the issue exists to fix;
it should not land unverified.

**Options, cheapest first:**

- Hold `inner: Box<dyn Write + Send>` and add a `#[cfg(test)]` constructor taking a writer plus an optional
  `PendingOutput`. The sink and the test both live in `src/bam.rs`, so this is a few lines. Dynamic
  dispatch costs one vtable call per ~1 MB write — unmeasurable.
- Make the sink generic over `W: Write`. Cleaner, but it pushes a type parameter onto
  `run_{single,paired}_end_parallel_to_bam`, which widens the diff past what this issue needs.
- Declare the write-error path untestable and say so explicitly in §12's remaining risks. Honest, but it
  leaves the issue's headline fix unguarded.

I'd take the first. A stub writer that fails on the Nth `write` also gives something today's panic-based
design cannot have at all: a *deterministic* assertion that the error text reaches the user. (I did check
the surrounding error plumbing and it is sound — a `write_compressed` error sets `first_error`, breaks
`'outer`, and is not masked by `reader_join_err`, because `read_{single,pairs}_round_robin` return `Ok(())`
on a failed send rather than `Err` (`src/parallel.rs:1226-1228, :867-872`). So the real `ENOSPC` will
surface. That is worth a test precisely because it is non-obvious.)

### I-3 — nothing outside `cargo test` sees this path, and every in-repo reader is the library that wrote the file

Per §3.4 above: no CI job runs `--output-format ubam` with `--cores > 1` through `ParallelBamSink`. And
every in-repo check — `assert_ubam_eq`, `bam_tuples`, both `inflate` closures — reads the output with
**noodles**, the same library that produced it. A container-level defect that noodles tolerates and htslib
does not would be invisible to the entire suite.

That is not hypothetical here: this change is specifically about container framing, and the plan's whole
safety argument rests on "the block layout changes; nothing depends on it".

**Suggested, one line of YAML:** add `--cores 4` to the existing "Validate SE uBAM-out — content-tuple
parity vs FASTQ-out path" step (`ci.yml:1259-1278`). It already pipes through `samtools fastq` and diffs
against the FASTQ path, so it becomes an independent-reader check on the parallel sink for free.

---

## 4. Optional

- **O-1 — §7's memory claim is directionally right for the wrong reason, and unmeasured.** "Each worker
  gains one `Vec` roughly half the encoded batch size" names the only *increase* and attributes the
  decrease to the pool's channels. The dominant term is neither: the `result_tx` channel is bounded at
  `cores * 2` and today carries **raw** batches (`with_capacity(reads * 2 * 400)` = 3.2 MB for paired),
  where after the change it carries **compressed** ones. At `--cores 8` that is ~16 slots dropping from
  ~3.2 MB to ~1 MB — tens of MB, swamping both the ~0.5 MB of pool channels removed and the ~8 MB of
  transient per-worker compression buffers added. Also, "half" is a loose ceiling; BAM seq+qual at level 6
  lands nearer a quarter. Either restate the reason or drop the claim — §9 measures neither.
- **O-2 — §9.4's threshold is too loose to detect anything.** Each batch boundary adds one extra partial
  block: ~26 bytes of framing plus slightly worse compression on a sub-64 KiB block, against ~433 KB of
  encoded records per batch. The honest expectation is **well under 1%**; "single-digit percent at most"
  would wave through a 5% regression that means something is actually wrong. Also, no 1M×150bp fixture
  exists in `test_files/` — say where it comes from (the untracked `10K_150bp.fastq.gz` in the tree is
  10K reads, not 1M).
- **O-3 — the `BgzfBlocks` doc comment in §4 breaks the project's comment rule.** Four lines, two of them
  justifying the change ("A newtype because … so the boundary is typed rather than documented"). Per
  CLAUDE.md that argument belongs in the commit message. One line: `/// BGZF members, no EOF block;
  [`ParallelBamSink::write_compressed`] appends them verbatim.` (§6 step 1's one-line invariant comment
  and §6 step 3's kept field-order comment are both fine as described — they state current invariants on
  the thing they constrain.)
- **O-4 — §6 step 6 undersells the test edits.** `src/bam.rs:2733` is not just a dropped argument:
  `sink.write_encoded(&enc.into_bytes())?` becomes `sink.write_compressed(&enc.into_bgzf()?)?`. The
  compiler catches it, but the step reads as mechanical when one line is not.
- **O-5 — the issue's fourth ask is silently dropped.** #542's closing paragraph asks for `--cores 1` and
  `--cores N` to share one compression implementation ("Two paths obliged to produce identical bytes are a
  standing liability"), citing the shared `build_record_buf` precedent. After this change they are still
  two implementations producing *different* container bytes by design (§7, §8). §11 says "Critical — none
  outstanding" and lists only the two scope questions put to Felix. Either state that unification is out of
  scope and why, or file it — the serial path's double EOF block (I-1) is a small live example of what two
  implementations drift into.
- **O-6 — §9.6 measures the smaller of the two overruns.** The issue reports 29 threads SE and **39**
  paired. The arithmetic checks out (SE = 1 main + (N+2) reader pool + (N+1) trim + (N+1) deflate pool =
  3N+5 = 29 at N=8; minus N+1 = 2N+4 = 20), but §1's table claims the overrun "goes away" for both shapes
  and only SE is measured. Also pin *how* the count is sampled, so the before-figure is comparable to the
  issue's. And note that even after the change the SE path is 2N+4 against a documented N+4 (`src/cli.rs:358-361`)
  — reduced, not reconciled; #568 should say so.
- **O-7 — minor.** `BufWriter<File>` earns little once every `write_compressed` is hundreds of KB (writes
  larger than the buffer bypass it); it is still right for the 49-byte header member and the 28-byte
  terminator, so keep it, but a reader will wonder. And `write_compressed`'s context string (inherited
  from `write_encoded`: `"failed to write BAM records"`) carries no path — worth adding, since reporting
  the real `ENOSPC` usefully is half of what C-1 is about.
- **O-8 — §8's public-API break.** `ParallelBamSink` must stay `pub` (the bin crate reaches it through
  `trim_galore::bam::`), so the break is unavoidable — but the crate ships on crates.io, so it is
  semver-relevant. Decide whether the changelog entry mentions it. §7's changelog guidance ("one sentence,
  the thread figure is the substance") is otherwise exactly right.

---

## 5. Things the plan got right that are worth not losing in review

- **§3.1 and §12's self-review caught a real defect in the issue's own sketch.** That is the single most
  valuable thing in the document, and it is correct. My quarrel is only with the severity ranking.
- **`flush()` + `into_inner()` over `try_finish()` for batches**, decided on the empty-batch case. Probe
  rows 4 and 6 confirm the reasoning precisely: `try_finish` on an empty encoder emits a stray 28-byte EOF
  block, `flush` + `into_inner` emits nothing.
- **The `BgzfBlocks` newtype.** Right call, and for the right reason — the sink genuinely cannot tell
  uncompressed bytes from compressed ones, and the rename alone would not have stopped a caller.
- **The order of `finish`:** terminator → `BufWriter::flush()` → `commit()`. Note this also fixes something
  the plan does not claim: today `ParallelBamSink::finish` takes the `BufWriter` back out of
  `MultithreadedWriter::finish()` and **drops it**, so a failure on that final flush is swallowed by
  `BufWriter::drop`. The explicit flush closes that hole too. Worth adding to the §1 table under H-2.
- **§12's "the validation matrix cannot see this path … the tests in §9 are the entire guard"** is correct
  and, per I-3, truer than stated.
- **§9.5's "md5 the binaries first"** and §12's note that `target/release/trim_galore` is stale (I confirmed:
  dated 12 Nov, three refactors behind `HEAD`). Exactly the right instinct.

---

## 6. Action items

**Critical — do before implementing**

1. Rewrite §9.1's assertions per C-1: non-empty blocks, inflate-equals-`into_bytes()`, **zero** EOF blocks
   in a batch, exactly one at the end of a finished file; derive the 28 bytes from an empty writer. Fix the
   accompanying prose in §3.1 and §12 — the omitted flush is the silent catastrophe (0 records, exit 0,
   `quickcheck` clean, on any input under ~64 KiB per batch), and per-batch EOF costs 28 bytes per batch
   and is read correctly by both samtools and noodles.

**Important — decide before implementing**

2. Give the terminator a pure helper (`bgzf_eof() -> io::Result<Vec<u8>>`) so the `into_inner` discipline
   lives in one shape, and note in the plan that §9.1's "exactly one EOF block" must not be applied to the
   serial writer, which emits two (verified on shipping output).
3. Resolve §9.3: either give `ParallelBamSink` a writer-shaped test constructor (`Box<dyn Write + Send>`
   inner + `#[cfg(test)]` constructor) and assert the error surfaces with its context, or state in §12 that
   C-1's fix lands unguarded. Correct the `unwritable_file` reference to `src/fastq.rs:814` and note it is
   `#[cfg(test)]`-scoped, hence unit-test-only.
4. Add `--cores 4` to one samtools-based uBAM-out CI step, so something other than noodles reads a
   `ParallelBamSink` file. Record in §3.4 that no CI job currently exercises this sink at all.

**Optional**

5. O-1 (restate or drop the memory claim) · O-2 (tighten §9.4 to <1%, name the fixture) ·
   O-3 (one-line `BgzfBlocks` comment) · O-4 (spell out the `bam.rs:2733` edit) ·
   O-5 (state the serial/parallel unification decision, or file it) · O-6 (sample paired threads too;
   pin the method) · O-7 (path in the error context) · O-8 (changelog and the crates.io API break) ·
   plus: add the swallowed final-flush fix to §1's H-2 row, and drop the payload-dependent 66/38 figures
   from §3.1 so they never reach an assertion.
