# PLAN_REVIEW_B — `--passthrough` mode

**Reviewer:** B (independent, fresh context).
**Plan under review:** `plans/06062026_passthrough-mode/PLAN.md`.
**Grounding:** read `src/cli.rs`, `src/trimmer.rs`, `src/parallel.rs`, `src/io.rs`, `src/main.rs`, `src/report.rs`, `src/fastq.rs`, `src/fastqc.rs` end-to-end before reviewing.

Overall: this is a solid, careful plan. The architecture (Option-threaded third reader/writer; reuse `classify_paired` / `process_pairs` as the sink-agnostic seam; ride the new stats inside `PairValidationStats`; piggy-back on the existing `BTreeMap<u64, PairedBatchResult>` ordered flush) is the right shape, and the choice of incompatibility envelope is defensible. The plan is also right that serial/parallel parity is the load-bearing invariant here, not Perl byte-identity.

There are, however, several real gaps that will bite the implementer if not surfaced now. I rank them at the bottom.

---

## 1. Logic review

### 1.1 Step ordering inside `classify_paired` / `process_pairs` — the plan is slightly misaligned with the actual code

The plan (Behavior §7) says: "Read three records → sync-check → `classify_paired` → on Pass write passthrough; on Discarded drop". That's right at the *write* layer (`process_pairs`, `parallel.rs:435`), but it implies the sync-check happens *inside* `process_pairs`. Look at the actual codebase: in the parallel path, batches of `Vec<FastqRecord>` arrive at the worker already-decoded; the worker has no source-side reader to drive a "three readers in lockstep" loop. The plan does later acknowledge this (Implementation §6.3: "sync check happens in the **reader thread**"), but Behavior §7 is sloppy enough that the implementer could get confused and add a worker-side check too.

**Concrete suggestion:** rewrite Behavior §7.iii to say: "(a) in the serial path, sync-check inline at the point of `classify_paired`; (b) in the parallel path, sync-check in the reader thread before the per-pair `batch_passthrough[i]` index is locked in; the worker side just consumes the parallel `Vec` by index." This also clarifies the stat-tracking site (see §1.2 below).

### 1.2 `passthrough_records_checked` accounting in the parallel path is under-specified

Under-current plan: the per-batch `pair_stats: PairValidationStats` is initialised in `process_paired_batch` (`parallel.rs:206`) — workers increment it as they classify. But the plan's reader-side sync check is the producer of `passthrough_records_checked`. Two questions the plan does not answer:

- Does the reader thread increment `passthrough_records_checked` into a shared atomic, or does it stamp the count onto each `PairedWork` payload? Either works, but neither is mentioned.
- More importantly, the reader currently does not produce any `PairValidationStats` at all — it just sends `PairedWork = (u64, Vec<>, Vec<>)`. To merge a reader-side counter cleanly the plan needs *either* (a) an `assumed_checked` field on `PairedBatchResult` populated by the worker from the batch length, or (b) a final counter merged into `total_pair` in the main thread after the reader joins.

**Concrete suggestion:** decide. The cleanest answer is "the worker stamps `pair_stats.passthrough_records_checked = batch.len()` for the worker's own batch when the third stream is present, because at that point the reader has already verified per-record sync." That keeps the merge logic completely local to `PairedBatchResult::pair_stats` and avoids a new cross-thread counter.

### 1.3 Edge case: passthrough EOF in `(None, None)` arm is asymmetric with R1/R2 EOF semantics

Implementation §5.5 says "The `(None, None)` arm: also consume from the passthrough reader and bail if it returns `Some`." Good. But the symmetric case — passthrough returns `None` while R1 returns `Some` — has to be checked at the *top* of the iteration, **before** matching on `(rec1, rec2)`. The current plan suggests doing the three-way next + match-on-mismatch at step 7.ii, but the existing serial code (`trimmer.rs:392`) does `let rec1 = ... ; let rec2 = ...; match (rec1, rec2)`. The implementer needs to thread the passthrough `next_record()` into that pattern, e.g. read all three first, then match on the 3-tuple. The plan doesn't enumerate the 8 cases of the 3-tuple, only the 4 R1/R2 cases. The implementer will have to invent the cross-product. Recommend the plan explicitly list all eight `(Some/None)^3` arms and what each does (six are "bail: file X is shorter than file Y", one is `(Some,Some,Some)` proceed, one is `(None,None,None)` exit).

### 1.4 `_R1.fq.gz` passthrough collision is plausible, not "won't actually collide"

Behavior §4 parenthetically says: "e.g. `I1.fq.gz` and `R1.fq.gz` both stripping to `I1_passthrough.fq.gz` and `R1_passthrough.fq.gz` — won't actually collide". True in *that* example. But in the realistic Multiome layout, the user has `R1.fq.gz`, `R2.fq.gz`, `R3.fq.gz` (10X labels the cell-barcode file `R3`, not `I1`). With `passthrough_output_name` using the passthrough input's basename via `strip_fastq_extensions`, `R3.fq.gz` → `R3_passthrough.fq.gz`, fine. But if the user passes the wrong file ordering — e.g. `--passthrough R1.fq.gz` and positional `R1.fq.gz R2.fq.gz` — the pre-flight does NOT currently catch the **input** duplication. Validation §1 step 1.ix only catches it if the implementer is precise about "byte-equal `PathBuf`"; on case-insensitive filesystems "byte-equal" is unsafe. Recommend §1.ix says: "passthrough cannot equal R1 or R2 — compare via the same case-folded normalisation the output-collision pre-flight uses, since pointing `--passthrough` at R1 on APFS would silently double-consume one input." This is a logic gap, not a typo.

### 1.5 The "byte-identical passthrough output" claim is wrong; canonicalised-record-identical is the truth

Assumption #9 ("The passthrough record is written byte-identically to its input — no trimming, no re-encoding"). This is *not* true given the codebase's I/O layer. `FastqReader::next_record` (`fastq.rs:325–355`) reads four lines, trims `\n`/`\r`, and constructs a `FastqRecord { id, seq, qual }` — it *discards* the third (`+...`) line entirely. `FastqRecord::write_to` (`fastq.rs:53–70`) re-emits the canonical four-line block with `\n` separators and a bare `b"+\n"` on line 3. Consequences:

- **Line 3 is canonicalised.** If the passthrough input has `+SRR123.1` on line 3 (some sequencers and downstream tools do), it becomes `+`. Probably fine for 10X Multiome (the cell-barcode files are Illumina-standard with bare `+`), but the plan should not promise byte-identity it can't deliver.
- **Line endings are canonicalised.** `\r\n` (rare but real on Windows-touched FASTQ) becomes `\n`. Same caveat.
- **No-op rewrite of the FASTQ structure.** Records ARE preserved record-identically (id, seq, qual all faithful), but the file bytes are not.

**Concrete suggestion:** soften Assumption #9 to: "The passthrough record's id/sequence/quality are preserved verbatim. Line endings are canonicalised to `\n` and the `+` line is rewritten to bare `+` (consequence of going through the standard FASTQ reader/writer; matches the canonicalisation R1/R2 already get)." This also affects Validation §1's "every passthrough record's ID matches R1/R2 IDs of the same row" — that test will work, but a test checking byte-identity to input would silently fail.

### 1.6 `--basename` interaction with multi-pair is already-rejected, but the plan double-asserts it

Validation step 1.ii says "input.len() != 2 → bail". Existing `Cli::validate` (`cli.rs:438`) already bails when `paired && input.len() > 2 && basename.is_some()`. Together this means there is exactly one code path where `--passthrough` runs (single pair) and `--basename` rejects the multi-pair case redundantly. Fine, just note the double-guard for clarity in the implementation outline. (Not a bug, observation.)

### 1.7 `--cores 1 + --clumpify` interaction is fine but not explicit

The plan asserts (Validation §6) "—clumpy + passthrough rejected at validation". Good. Note that the existing dispatcher (`main.rs:732`) goes through `run_paired_end_parallel` whenever `cli.cores > 1 || cli.clumpify`. So with `--passthrough` you have:
- `--cores 1` → serial `trimmer::run_paired_end` (plan §5 path).
- `--cores >= 2` → parallel `run_paired_end_parallel` with `clump_layout = None` (plan §6/§7 path).

The plan covers both, but it doesn't say "verify both branches are exercised by the test matrix". With the existing test fixture being a 50-record file, a single-batch `--cores 4` run will only emit one `PairedBatchResult` and thus not stress the BTreeMap ordering. **The parity test §2 should use a fixture large enough to span multiple batches** — at minimum 5,000 records (>`BATCH_SIZE` = 4096), ideally 10,000–20,000 to get 3–5 batches and exercise the out-of-order arrival path. Plan §11.2's "50 records" fixture is insufficient for §2.

### 1.8 `--retain_unpaired` rejection is over-tight

The plan rejects `--retain_unpaired` because "a half-rescued mate has no clean interpretation for the index file". Defensible, but the actual semantic is well-defined: if R1 is rescued and R2 is dropped, write the passthrough record to a third *unpaired-passthrough* sink — or just write it to the regular passthrough sink (since by definition the cell-barcode is for that single fragment regardless of which mate survived). Recommend the plan either (a) say explicitly "this was considered and deferred to v2 because no user has asked for it yet", or (b) write the passthrough record on **any** `r1_ok || r2_ok` rescue. I'd lean (a) — keep v1 conservative — but the plan should not pretend there's no clean interpretation.

---

## 2. Assumptions

### Explicit assumptions that hold up

- §1 strict 3-file arity, §2 paired required, §3 FASTQ-only, §4 gzip autodetect, §6 strict pair semantics, §7 sync mandatory — all reasonable, all consistent with the existing codebase patterns. The model on `--clumpify` (also: cores >= 2, paired required, specialty modes rejected) is correctly mirrored.

### Hidden / implicit assumptions the plan does not surface

1. **Output gzip mode is decided from R1 only.** `main.rs:84` derives the global `gzip` from `is_gzipped(&cli.input[0])` (which is R1 in paired-end mode). If R1 is plain `.fq` and `--passthrough I1.fq.gz` is gzipped, the passthrough **output** will be plain `.fq` (following the R1-driven decision), not `.fq.gz`. Document this. Better: in the validation block, warn if R1's compression and the passthrough input's compression disagree. Best: make passthrough output compression *follow R1's* per design, with a note.

2. **The reader thread allocates `Vec<FastqRecord>` for passthrough records but workers may not touch them.** Memory footprint of the parallel batch grows from "2 vecs of 4096 records" to "3 vecs of 4096 records" regardless of whether any pair survives filtering — this matters only at very high `--cores` (24+) and small RAM. Plan §"Efficiency" handwaves "negligible" (1.2 MB extra per batch); fine, but worth flagging that the passthrough records are held in memory through the whole worker-pool latency tail (decompress → channel → worker → result channel → flush), so the worst-case in-flight memory is `cores * 2 * 1.2 MB` ≈ 38 MB at `--cores 16`. Not "negligible" if a user is also running `--clumpify` (rejected here — fine).

3. **The header-prefix sync check stops at first whitespace AND first `/`.** The plan says "split on whitespace only" (§7.iii) but then defers stripping `/1`/`/2`/`/3` to a future change. With strict whitespace-only matching, the legacy `@read1/1` (R1), `@read1/2` (R2), `@read1/3` (passthrough) case **will fail the sync check**. The plan acknowledges this but defers. I'd argue this is a *day-1* bug for any data older than 2014. The fix is genuinely one line:

   ```rust
   pub fn read_id_prefix(id: &str) -> &str {
       let s = id.strip_prefix('@').unwrap_or(id).split_ascii_whitespace().next().unwrap_or(id);
       s.rsplit_once('/').filter(|(_, suf)| matches!(*suf, "1" | "2" | "3")).map(|(p, _)| p).unwrap_or(s)
   }
   ```

   Strongly recommend doing this in v1.

4. **`--fastqc` on the passthrough file may be misleading.** Cell-barcode reads are 16–28 bp of uniformly-structured sequence — FastQC will scream about per-base sequence content bias, adapter content (none, but the report won't say "expected"), and N-content (often 0% by design). The plan says "verify the bundled FastQC engine handles the file without complaint" — it will *run*, but the output report is going to look catastrophically bad. Consider either suppressing FastQC for the passthrough file by default with a flag to opt in, or documenting prominently that the report is expected to look like a horror show.

5. **Empty passthrough records are written verbatim.** §"Edge cases" mentions "Passthrough record with empty sequence: No special handling". Worth adding to validation: a fixture where a passthrough record has `seq == ""` should still run. (`FastqRecord::write_to` will emit a 4-line block with an empty seq + empty qual line, which is valid FASTQ but some downstream tools choke on it.)

6. **The `plans/` directory is gitignored.** The plan notes this at the bottom — good catch. But the plan also creates `multiome_*.fq.gz` fixtures under `test_files/`, which **is** tracked. The fixture-creation step should be its own commit to keep the diff reviewable.

---

## 3. Efficiency

### What the plan gets right

- The classification of cost (`split_ascii_whitespace` ~10 ns; ~30 ms per 1M reads) is accurate. Sub-percentage of total runtime.
- The memory math (~1.2 MB extra per batch) is correct.
- Reuse of the existing `BTreeMap<u64, PairedBatchResult>` ordered flush is the right call — no new ordering primitive needed.

### What the plan handwaves or misses

1. **Channel pressure.** The current `mpsc::sync_channel::<PairedWork>(2)` has capacity 2 per worker (`parallel.rs:85`). Adding a third `Vec<FastqRecord>` to each payload grows the per-channel-slot memory by 50%. At `--cores 16`, in-flight worker queue = 16 × 2 = 32 slots × 3 × 1.2 MB = ~115 MB just sitting in channels. Pre-existing channel was ~77 MB. Not catastrophic but worth a one-line note in §Efficiency, and consider whether `sync_channel` capacity should be reduced to 1 when the passthrough payload is active.

2. **`Vec<FastqRecord>` allocations for passthrough records that are dropped.** In `process_paired_batch`, the per-batch passthrough `Vec` is consumed by the worker; on `PairOutcome::Discarded` the worker just lets the record drop. That's an allocation lifecycle that never gets pooled. With a 50% drop rate (plausible on aggressive `--length 50` filtering of 36 bp cell-barcode reads — although they wouldn't be at risk; only R1/R2 lengths drive `PairOutcome`) you'd allocate-then-free `Vec<FastqRecord>` worth `BATCH_SIZE * 0.5` per batch. Fine on `mimalloc`/`jemalloc`, not free on glibc. Not a blocker, but a measurable allocator pressure when scaled to 10X-scale 1B-read runs.

3. **The third gzip encoder per worker.** Each worker now spins up an additional `GzEncoder` per batch. The plan says "no contention" — true for *threads*, but each `GzEncoder::new` allocates ~256 KB of internal state (deflate window + dictionary). At `--cores 16`, you're spinning up 16 more gz encoders × ~256 KB = ~4 MB additional steady-state RAM. Reasonable. But if any pair is dropped (Discarded), the corresponding passthrough `Vec<u8>` allocation in `PairedBatchResult` is still made (and gzip-finalised) with zero records. The "empty gzip member" gz overhead is ~20 bytes per empty member — negligible. Fine.

4. **Reader thread is now the throughput bottleneck.** Currently the reader thread does `next_record()` × 2 per pair. With passthrough, it does `next_record()` × 3 + a `split_ascii_whitespace` × 3 per record. At >100 MB/s decompressed FASTQ on `open_threaded`, this adds maybe 5% reader-thread overhead. Not blocking, but the reader is already the implicit serial bottleneck in `run_paired_end_parallel`; adding work there matters more than adding work in workers. Worth measuring on the benchmark fixture once implementation lands.

---

## 4. Validation sufficiency

The plan proposes 4 + 1 optional tests. They are individually well-designed, but the **coverage of the parallel-path-specific failure modes is incomplete**.

### What the proposed tests catch

- §1 end-to-end serial smoke — basic lockstep correctness.
- §2 serial/parallel parity — stat parity AND decoded-bytes parity. Strong invariant.
- §3 truncation — passthrough shorter than R1/R2.
- §4 shuffled IDs — sync-check active.
- §5 (opt) — validation rejections (the 9 `bail!` cases).

### What the proposed tests miss

1. **Cross-batch ordering hazard.** With BATCH_SIZE=4096 and a 50-record fixture, the parallel path emits exactly one batch and never exercises the BTreeMap-reassembly logic. **A passthrough-mode-specific cross-batch parity test is needed:** ≥10,000 paired records (3 batches at BATCH_SIZE=4096), drop 30% via `--length`, verify the passthrough output records emerge in the original input order AND in lockstep with R1/R2 — not just by stat counts but by zipping the decoded outputs and comparing row-by-row to a known-good filter trace.

2. **`--cores` jitter test.** Run the same input at `--cores 2`, 4, 8, 16 and assert all four produce identical passthrough output bytes (or at least identical-decoded records). The current parity test compares serial vs `--cores 4` only — but inter-cores parity catches batch-distribution bugs the cores-1-vs-N comparison can mask if the reader's distribution logic has an off-by-one.

3. **Reader-thread error propagation.** When the reader thread detects a truncated passthrough file, it must (a) abort *before* the workers fill the result channel with already-processed batches, and (b) ensure the main thread receives the error rather than hanging on `result_rx.recv()`. The current code's `reader_handle.join()` pattern (parallel.rs:186–190) joins **after** `result_rx.recv()` loop drains — so a reader error mid-stream will produce the wrong behaviour: the reader returns `Err`, but the result channel keeps draining (workers haven't been signalled), and main thread proceeds to writes happily with truncated output before checking `reader_handle.join()`. The plan §6.4 says "fail fast" but the existing harness does not actually fail fast on reader errors — it fails at *join time*, after all in-flight batches finish writing. **Add a test that asserts no output files are created (or are explicitly removed) when the reader bails mid-stream.** This is the most likely silent-failure mode of the whole feature.

4. **Header sync false-positive coverage.** Test that `@read 1:N:0:ATCG` and `@read 2:N:0:ATCG` and `@read 3:N:0:ATCG` (modern Illumina, post-2014, three reads of one cluster) successfully pass the sync check. This is the "golden path" of the Multiome use case and there's no test for it as currently planned.

5. **`--dont_gzip` interaction.** No test covers plain-text output. Add one.

6. **`--output_dir` interaction with collision pre-flight.** Two pairs (`--paired` rejected at v1, so this is implicit only-with-multi-pair, but the multi-pair rejection itself needs a test — covered by §5 optional).

7. **Header lines with unicode.** Real-world FASTQ rarely has unicode, but `split_ascii_whitespace` is well-defined on byte sequences. The plan uses `id: &str` (UTF-8 valid) for the prefix helper — `FastqRecord::id` is `String`, so this is fine. No test needed, but worth confirming the implementer doesn't try to use `split_whitespace` (which counts Unicode whitespace and would diverge subtly between R1 IDs containing a `U+00A0` vs not — extremely unlikely but should be ruled out by the helper's name).

### The most important test the plan does not have

**Reader-thread `bail!` propagation under load.** This is the test most likely to catch a real bug: sync error at record 5,000 of 10,000, `--cores 4`, ensure (a) the function returns `Err`, (b) no partial passthrough output remains, (c) all worker threads cleanly exit, (d) no deadlock. The current scope's `result_tx.clone()` + drop pattern survives most reader failures but the *unpaired* output writes happen unconditionally in the flush loop — there is a window where partial output can land on disk before the reader's error reaches the main thread.

---

## 5. Alternatives worth considering

### 5.1 CLI shape: positional vs flag

The plan uses `--passthrough <FILE> R1 R2`. Reasonable, mirrors `--demux <FILE>`. An alternative is positional triple: `--paired --passthrough R1 R2 I1` (third positional consumed only when `--passthrough` is a *bool* flag). The bool-flag form is harder to typo (no value-without-flag risk) but breaks the symmetric "`--paired R1 R2`" pattern. The `--passthrough <FILE>` form is the right call. **No change recommended.**

### 5.2 Header sync check: strip `/[123]`

Already addressed in §2.3 — strongly recommend doing this in v1.

### 5.3 Stat tracking placement: new struct vs extend `PairValidationStats`

The plan extends `PairValidationStats`. Alternative: introduce a `PassthroughStats` struct with `kept`, `dropped`, `checked` fields, threaded as a fourth return tuple slot. Trade-offs:

- **Extending `PairValidationStats`:** 3 fields added; existing serialisation needs an `Option<PassthroughBlock>` JSON wrapper anyway to keep MultiQC happy; the stats are conceptually paired-validation-stats (one per pair, gated by passthrough mode).
- **New `PassthroughStats` struct:** cleaner separation; touches `run_paired_end`/`run_paired_end_parallel` return signature (extra slot or wrapped enum); requires a separate `merge()` call.

I narrowly favour the plan's choice. The new fields are paired-mode-specific (`passthrough_records_checked == pair_stats.pairs_analyzed` when enabled, 0 otherwise), so their lifecycle matches `PairValidationStats` exactly. **No change recommended.**

### 5.4 Worker payload: `Option<Vec<FastqRecord>>` vs always-present empty Vec

The plan picks `Option<Vec<FastqRecord>>` for `PairedWork`. Alternative: always include an empty Vec when passthrough is disabled. The `Option` form is marginally cleaner (the worker can skip the per-record write loop entirely when `None`), but adds a layer of `if let Some(p) = passthrough_records { ... }` boilerplate to the worker. The plan also picks `Vec<u8>` (empty when disabled) for `PairedBatchResult::compressed_passthrough` — which is the opposite choice. The asymmetry is slightly inelegant but defensible (`PairedBatchResult` has no awareness of whether passthrough is active; the main thread does). **No change recommended; just call out the asymmetry in the implementation outline so the reviewer doesn't flag it.**

### 5.5 An alternative architecture: passthrough as a side-channel post-filter

A radically different approach: do *not* thread the passthrough through `run_paired_end_parallel`. Instead, after the trimming pass writes R1/R2 outputs, do a second pass that streams the passthrough input + decoded R1 output, matching read IDs and emitting the passthrough records whose R1 ID survived. Trade-offs:

- **Pro:** zero changes to `run_paired_end_parallel`. The feature is fully orthogonal.
- **Pro:** trivial to add `--retain_unpaired` support later (just stream against R1+R2_unpaired union).
- **Con:** second pass = 2x read of the passthrough input. ~10–15% wall-clock overhead.
- **Con:** the sync check moves from "early fail" to "post-trim fail" — by the time you detect a desync, R1/R2 outputs are already written.
- **Con:** memory cost is O(N) for the R1-output-id set (or O(1) with streaming if both files are already in input order; but you have to *trust* that they are, which is the whole problem).

The plan's in-line approach is the right call given the constraints (early-fail, single-pass, no temp-set memory). **Mention this alternative in the plan's "Things considered but not chosen" subsection** so the implementer understands why the in-line threading is the chosen shape.

### 5.6 Defer the parallel path entirely for v1

Genuinely worth asking: is the v1 use case (10X Multiome ATAC pre-processing for nf-core/atacseq) bottlenecked by single-thread throughput? A typical Multiome ATAC library is 200–500 M paired reads. At a serial throughput of ~3M reads/min (cores=1 trimming, no compression bottleneck) that's 60–170 min — uncomfortably long. **Parallel path is necessary in v1.** OK, ignore this alternative.

---

## 6. Action items, prioritised

### Critical (must address before implementation)

1. **Plan the reader-thread error propagation explicitly.** The plan handwaves this. The implementer must (a) signal workers to stop on reader error, (b) drain the result channel without writing partial outputs to disk, (c) ensure `reader_handle.join()` is checked *before* output files are committed. Add a test specifically for this. (§4.3, §4.7)

2. **Enumerate the 8-arm three-way EOF/Some match.** Behavior §7 leaves it implicit; the implementer will reinvent it (probably correctly, but maybe not). Spell out all 8 cases. (§1.3)

3. **Decide who owns `passthrough_records_checked` accumulation in the parallel path.** Reader thread, worker, or main thread? The plan is silent. Recommend "worker stamps from batch length on the assumption the reader has already verified sync." (§1.2)

4. **Fix the passthrough-equals-R1/R2 collision check to be case-folded, not byte-equal.** Otherwise APFS users will silently dual-consume one input. Validation §1.ix needs the `normalize_path` helper from the existing pre-flight. (§1.4)

5. **Soften Assumption #9 from "byte-identical" to "record-identical with line-3 and line-ending canonicalisation".** Otherwise the plan asserts an invariant the codebase cannot deliver. (§1.5)

### Important (should address before implementation)

6. **Bigger parity test fixture (≥10K records, multi-batch).** A 50-record fixture cannot exercise the BTreeMap reassembly in `--cores >= 2`. (§1.7)

7. **Decide on `/1` `/2` `/3` suffix stripping in the header sync check.** One-line change; recommend v1. (§2.3)

8. **Add a multi-cores jitter parity test (`--cores 2/4/8` produce identical decoded output).** Catches batch-distribution off-by-ones. (§4.2)

9. **Add a passthrough-mode large-fixture cross-batch test that asserts row-by-row passthrough/R1/R2 ID parity in the output.** Stat counts alone don't catch a "passthrough records emerge in scrambled order vs R1/R2" bug. (§4.1)

10. **Document the gzip-mode-driven-by-R1 behaviour.** Add to Assumptions; consider warning in the validation block when R1 and passthrough disagree on compression. (§2.1)

### Optional (improves quality but not blocking)

11. Mention the FastQC output will look poor on cell-barcode reads — either suppress by default or document loudly. (§2.4)

12. Add a `--dont_gzip` end-to-end smoke test. (§4.5)

13. Consider lowering the `sync_channel` capacity from 2 to 1 when the passthrough payload is active, to bound in-flight memory at high `--cores`. Measure first. (§3.1)

14. In the plan's "Things considered but not chosen" subsection, mention the side-channel post-filter alternative (§5.5) so reviewers understand why in-line threading was chosen.

15. The "Things explicitly OUT of scope for v1" section omits `--demux` + `--passthrough` and `--retain_unpaired` + `--passthrough` rationale. Add the "what would a clean rescue/demux semantic look like" sentence so the deferral is reasoned, not just "rejected". (§1.8)

16. The `passthrough` field name on `Cli` collides linguistically with the existing `passthrough` semantics ("specialty modes pass through to their own dispatchers"). Not a real bug, but consider `passthrough_file` or `carrier` or `barcode_passthrough` for the field name to keep `cli.passthrough` from being grep-confusing. (Tiny nit.)

---

## Summary

The plan is well-structured and most of the architectural choices are correct. The most serious gaps are around the parallel-path's reader-thread error propagation, the under-specified `passthrough_records_checked` accounting in the parallel path, the incomplete 3-way EOF match enumeration, and an over-confident "byte-identical" assumption that the codebase's I/O layer does not actually support. The test plan needs a multi-batch fixture and a reader-error test to cover the parallel-path-specific failure modes. Header sync should strip `/[123]` in v1, not v2.

None of these are blocking; all are correctable in the plan itself before implementation.
