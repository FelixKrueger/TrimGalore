# PLAN_REVIEW_A — Phase 1: TrimGalore input + output formats

**Reviewer:** A (independent, fresh context window)
**Plan file:** `plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN.md` v1
**Date:** 2026-06-25

---

## Critical findings

### C1 — §5 step 1's "pure refactor" claim is architecturally false: `parallel.rs` workers do not have a writer-trait seam today

**Severity:** CRITICAL.

The plan asserts (§5 step 1.3, §6, and §7):

> Update `parallel.rs::run_*_parallel` to take `Box<dyn RecordSink>` for output (mirror of how it took `Box<dyn RecordSource>` in #317). Same move-by-value pattern…
> Step 6 verification gate: `cargo test --release` — 296 tests still pass with no functional change (this is a pure refactor; FASTQ-out byte-identity must be exact).

**This is not symmetric to the input side, and the symmetry argument breaks under inspection.**

Reading `src/parallel.rs::process_paired_batch` (lines 305–427) and `process_single_batch` (lines 981–1006), each worker:

1. Allocates a per-batch `Vec<u8>` buffer.
2. Wraps it in a `GzEncoder<&mut Vec<u8>>` (when `gzip` is true) — the compressor lives **inside the worker**, on a per-batch basis, owned by the worker thread.
3. Calls `record.write_to(&mut gz)` per record through a generic `W: Write` parameter (`process_pairs`, `process_reads`).
4. Calls `gz.finish()` to seal the gzip member.
5. Returns the resulting `Vec<u8>` (already-compressed gzip bytes) through `mpsc::SyncSender<PairedBatchResult>` to the main thread.
6. The **main thread** does the actual `File::create` + `write_all` of the raw compressed bytes, **with no writer object at all** — the file is just `std::fs::File`.

The `RecordSink` trait as proposed in §4 (single `fn write_record(&mut self, record: &FastqRecord) -> Result<()>` + `fn flush(&mut self) -> Result<()>`) does NOT fit this model. Concretely:

- **`FastqWriter::create` wraps the underlying `File` in either `ParCompressBuilder<Gzip>` (multi-thread parallel gzip) or `BufWriter<GzEncoder>` (single-thread).** Constructing a `FastqWriter` per worker would *double-gzip* the output (each worker emits a gzip member, then the writer would re-gzip it on top), AND would dispatch parallel compression on top of parallel workers. Both are wrong.
- **Constructing a single `FastqWriter` on the main thread and sharing it across workers** requires `&mut dyn RecordSink` to be shared across threads with ordering — that's exactly the BTreeMap ordered-flush pattern that `mpsc::Receiver<PairedBatchResult>` already implements, and it requires the worker→main payload to be **bytes**, not records. The proposed trait method signature operates on records.
- **The current worker keeps records and the compressor in the same thread for cache locality** — moving the compressor out of the worker would funnel compression back through the main thread, undoing the architecture commit message that explicitly says "the dominant cost (gzip compression, ~60% of runtime) is distributed across workers instead of funneled through one thread" (parallel.rs lines 10–12).

The plan needs to choose:
- **Option A:** `RecordSink` is a *main-thread-side* abstraction that owns an `Fn(&Path) -> Result<Box<dyn RecordSink>>` factory; workers continue to use a separate, lower-level **chunked-bytes contract** (e.g., a trait that produces a `Vec<u8>` chunk per record-batch). Then `RecordSink` is for the serial path and specialty paths; the parallel path uses a parallel-specific seam.
- **Option B:** Pull the writer out of the worker entirely. Workers produce trimmed `Vec<FastqRecord>` batches and hand them to the main thread, which dispatches through `&mut dyn RecordSink`. This is correct but undoes the architecture in `parallel.rs` and would be a major perf regression — gzip would re-funnel through the main thread.
- **Option C:** Make `RecordSink` two-layered — `write_chunk(records: &[FastqRecord]) -> Result<Vec<u8>>` for the worker (returns the encoded bytes for that chunk) + `write_chunk_bytes(bytes: &[u8]) -> Result<()>` on the main-thread sink. This is workable but it's not the single-method trait the plan promises.

Until this is resolved, "296 tests pass with no functional change" cannot be the verification gate, because the refactor proposed isn't actually feasible as written. **This is the load-bearing architectural decision for Phase 1 and must be made explicit before any implementation begins.**

### C2 — `parallel.rs` `--cores > 1` FASTQ output uses `ParCompressBuilder` for *serial* path but per-worker `GzEncoder` for *parallel* path; the `RecordSink` factory must distinguish

Closely related to C1, but separable. The serial path (`run_single_file` line 668 in main.rs) goes:
```
FastqWriter::create(&path, gzip=true, cores=cli.cores, level)
  → ParCompressBuilder if cores>1, else GzEncoder
```
The parallel path (`run_single_end_parallel`) does NOT use `FastqWriter` at all — it bypasses it and builds a `GzEncoder` per worker inside `process_single_batch`, then `File::create(output)` + `write_all` on the main thread.

So `FastqWriter` impls `RecordSink` only services the serial path. The parallel path needs an entirely different sink shape. The plan does not distinguish these, and silently conflates them in §5 step 1.

The `MimIndexWriter<W: RecordSink>` design in §4 / §5 step 4 inherits this confusion. If `MimIndexWriter` wraps a `FastqWriter` to track byte offsets, that only works in the **serial path**; in the parallel path, byte offsets cannot be known until after the main thread concatenates the per-worker `Vec<u8>` chunks (and even then, the offsets are *compressed* offsets, not record offsets). This is precisely what Spike A2 is supposed to discover, but Spike A2 needs to know it's looking at **two different write paths**, not one.

### C3 — Spike A2 success criterion is undefined

The plan (§5 step 0, Spike A2):

> Check whether that's tractable inside `flate2`'s GzEncoder (which manages compression internally) or whether we need BGZF-framed output for mim to work at all… **This spike's outcome may push mim sidecar production into a separate post-pass over the produced FASTQ.gz, in which case adjust §3.2 / §5 step 4 accordingly.**

There's no measurable gate. What does "tractable" mean — is "we can call `GzEncoder::total_in()` after every record" sufficient, or do we need byte-precise compressed offsets per record? The mim format itself isn't characterised in this plan — Phase 1 cannot scope until we know:

- Are mim offsets *uncompressed* byte offsets in the FASTQ stream (then `total_in()` works) or *compressed* byte offsets in the gzip stream (then they require BGZF framing because mid-block decode is impossible)?
- Does mim sidecar work record-level or block-level?
- Does the parallel-output path (independently-compressed gzip members) defeat mim entirely, even with BGZF?

Without a documented mim-format spec inside the plan, this spike is open-ended. A reasonable success criterion: "produce a `.mim` sidecar against `parallel.rs`-generated multi-member gzip output and verify a downstream mim reader can seek to record N within ±0 bytes." If that's infeasible, *neither* Path A nor Path B in §5 step 4 helps, because both assume the offsets are knowable.

### C4 — Spike A3's gate threshold and consequence are stated but not load-bearing

§9 row "Binary size delta" says "≤ 2 MB (per Spike A3)", and §10 Q "Default-on vs feature-gated for BINSEQ / mim" says "Plan assumes default-on for prebuilt; flip if delta > 2 MB."

This is fine as far as it goes, but the consequence of *flunking* the gate is under-specified. If `binseq` adds 8 MB, what happens?

- Cargo features stay `default = ["binseq", "mim"]`?
- Prebuilt release binary drops them (the most likely sensible answer)?
- Plan is rewritten?

The CI matrix gain in §5 step 7 ("two new jobs: `rust-tests-slim` and `rust-tests-full`") and §9 cover the *build* side, but not the prebuilt-release-binary decision. Recommend a clear branch: if delta > 2 MB, prebuilt enables only `noodles+bam` (already present) and source builds opt-in to `binseq` / `mim`.

---

## Open findings (worth discussing, not blocking)

### O1 — §3.5 / §11 `--preserve-tags X,Y` + FASTQ-input + `--output-format ubam` hard-error may misfire on legitimate mixed-input pipelines

The §11 self-review explicitly tightens this:

> Originally allowed `--preserve-tags X,Y` with FASTQ input + `--output-format ubam` (just emit empty tags). Tightened to hard-error: if user requests tag preservation with no source tags, that's a bug in their pipeline, not something to silently accept.

This is reasonable for single-input runs. But TrimGalore supports **multiple input files in one invocation** (single-end loop over each input; `--paired` consumes pairs of files). What about:

```
trim_galore --paired --output-format ubam --preserve-tags CB,UB \
  fastq_pair_R1.fq.gz fastq_pair_R2.fq.gz \
  ubam_pair_R1.bam   ubam_pair_R2.bam
```

The user is asserting "carry CB/UB if the input has them; otherwise produce empty tags" — that's a legitimate batch-job policy, not a bug. If we hard-error on the FASTQ pair, we break the batch.

Recommend: hard-error only when *every* input is FASTQ AND `--preserve-tags` is set AND `--output-format ubam`. If at least one input pair is uBAM, accept and emit empty tags for the FASTQ pairs (with a stderr note per pair).

Alternatively: just emit empty tags + stderr warn, like the original design. The "bug in their pipeline" framing isn't airtight when batch-job policy is the use case.

### O2 — §3.3 paired-uBAM output as two separate BAMs deviates from BAM ecosystem convention

The §11 self-review notes:

> Originally had paired-uBAM output as one interleaved BAM. Changed to two separate BAMs (mirrors FASTQ paired-end output shape) — symmetric naming, simpler.

But the ecosystem (samtools, picard, fgbio, ChromBPNet, Cell Ranger) overwhelmingly treats paired BAM as **mate-adjacent in one file**. The TrimGalore uBAM-input pathway in #317 already reads single-interleaved BAMs via `bam::open_paired_interleaved` (see `src/bam.rs:131`). Emitting two BAMs from one input BAM would require downstream consumers to re-interleave them before any aligner can consume them — which is a weird shape.

Recommend: revisit. Either emit one interleaved BAM (matches ecosystem; symmetric with input shape; mate-adjacency preserved naturally) or **emit one interleaved BAM by default with `--separate-bam-outputs` as an opt-in for FASTQ-symmetric naming**. The "simpler" framing isn't actually simpler for users downstream.

### O3 — §3.4 `--output-format binseq` + `--demux` not in the rejection list

§3.4 enumerates blocked combinations: `--clumpify`, `--passthrough`, `--emit-mim` + non-FASTQ. But `--demux` is also a per-format output path (writes one FASTQ per barcode to a directory). What's the behaviour of `--demux` + `--output-format binseq`? The plan doesn't say.

Two reasonable answers — either is fine, but pick one:
- Reject with a clear error (matches the `--clumpify` pattern).
- Allow it: write one `.binseq` per barcode.

Same for `--demux` + `--output-format ubam` (one `.bam` per barcode). Currently silent.

### O4 — §5 step 6 underestimates specialty-mode integration cost

§5 step 6.2 says:

> Thread the factory through `run_single_file` / `run_paired` / `run_paired_ubam_single_file` / specialty modes.

Reading `src/specialty.rs`, all four specialty modes construct `FastqWriter::create` directly (lines 33, 76, 137–138, 250–251) and assume FASTQ output. Threading a writer factory through means:

1. **Output naming** is currently per-specialty (`hardtrim_output_name`, `clock_output_name`, `implicon_output_name` in `src/specialty.rs:309-340`). Each name encodes `_hardtrim5`, `_clock_UMI`, `_implicon`. These naming functions assume a single output suffix. Extending to `.bam` / `.binseq` / `.fq.gz` per format means either (a) bake format-aware suffixes into each helper, or (b) lift naming into the factory. The plan doesn't say which.
2. **Specialty modes have semantic invariants that may conflict with BAM/BINSEQ output.** `--clock` writes a UMI into the FASTQ header via the `_UMI_<seq>` suffix — that's a FASTQ-header convention. In BAM, a UMI goes into a BAM aux tag (`RX:Z:<seq>`). The plan doesn't address whether `--clock` + `--output-format ubam` translates the UMI into an aux tag, drops it, or rejects.
3. **`--implicon` has the same UMI-in-header pattern.**
4. **uBAM output for specialty modes also need `source_header: Option<&Header>`** — which file's header gets propagated when the specialty mode consumes multi-pair input?

This isn't "low-cost" — it's a 4-mode × 3-format × naming + semantics matrix that needs design. Either flag specialty-mode format support as out-of-scope for Phase 1 (FASTQ-only) or add an explicit sub-step for it.

### O5 — §6 efficiency claim about `RecordSink` parallels `RecordSource` overhead measurement, but only for the serial path

§6 says:

> `RecordSink` trait-object dispatch mirrors the `RecordSource` overhead measured in #317's Spike 1 (+1.28%), well within budget.

This is true *only if* the trait-object dispatch happens at record granularity (which is what the proposed signature does). But the parallel path currently dispatches at *chunk* granularity through `mpsc` (one send per batch of 4096 records). Switching to per-record dispatch on the main thread would be a 4096× channel-call increase. The plan should either keep the parallel path's chunk-level dispatch (which means `RecordSink` doesn't help there, see C1) or measure the per-record-dispatch overhead under load.

### O6 — `--output-format ubam` + `--paired` with two FASTQ files: BAM header source ambiguous

Per §4's `BamWriter::create` signature:

```rust
pub fn create<P>(path: P, source_header: Option<&noodles::sam::Header>, …) -> Result<Self>;
```

In `--paired` mode with two FASTQ inputs, `source_header` is `None` for both outputs (synthesised header). That's fine. But what about `--paired` with mixed input (one FASTQ pair + one uBAM pair, in the same invocation)? The plan implicitly assumes a single input pair per process. Recommend: explicit confirmation that `--paired` BAM-output uses **per-pair header synthesis** — each output pair gets its own synthesised header based on its own input's source. State it.

### O7 — `--output-format ubam` + paired-end interleaved uBAM input (#317's path) — round-trip story

#317 added `open_paired_interleaved` to de-interleave a single uBAM into two `BamReader` streams. If a user runs:

```
trim_galore --paired --output-format ubam input_interleaved.bam
```

…they get two separate output BAMs per §3.3. The downstream consumer who originally provided the interleaved input now has two files. This is the symmetry-break I flagged in O2, but it has a specific consequence here: the round-trip `interleaved-ubam → trim_galore → interleaved-ubam` is impossible without an explicit re-interleave step.

This may be fine (we don't owe users a 1:1 round-trip), but document it clearly. If we decide to emit interleaved BAM (O2), this concern dissolves.

### O8 — `--emit-mim` + serial path (`--cores 1`)

§10 question item:

> What happens if `--emit-mim` is specified with `--cores 1` (serial path)? — Plan assumes it works the same as parallel path. Mim writer is sink-side; serial vs parallel doesn't matter.

Given the analysis in C1 / C2, this question becomes more pointed: in the serial path `MimIndexWriter` wraps a `FastqWriter` (which itself wraps `BufWriter<GzEncoder<File>>`) — record-to-compressed-byte-offset tracking is feasible because there's one continuous gzip stream. In the parallel path, the output is N gzip members concatenated by the main thread, and the worker that compressed batch K doesn't know the starting byte offset of batch K in the final file. The plan needs to acknowledge that `--emit-mim` in the parallel path requires either (a) post-pass indexing or (b) the main thread tracks per-batch compressed sizes as it writes (Phase 1 should add this — it's cheap).

---

## Validation gaps

The §9 validation table is reasonable but missing:

- **`process_paired_batch` parity test for each output format.** The PR #246 / commit 82d1e34 saga (referenced in `parallel.rs` tests around line 1278) shows how easy it is to drift stats between serial and parallel paths. Each new output format adds a new variant of this risk. Recommend: extend `test_parallel_serial_trim_stats_parity` to cover `--output-format binseq` and `--output-format ubam` (these tests use the in-process API, not the binary).
- **`output_format_binseq_plus_clumpify_rejected`** is in §5 step 8 but the same test pattern is missing for `passthrough`, `emit-mim` × non-FASTQ, and the `--preserve-tags` rejections from §3.5. Each rejection rule needs its own test (mirroring the `parse_sam_tag_name` discipline from #317).
- **Reproducibility delta.** CI has a `reproducibility` job that builds twice with `SOURCE_DATE_EPOCH` set and asserts bit-identity. The plan says it "runs unchanged" — but new optional crates (binseq, mim-index) MUST also be reproducible. Each crate's build-graph determinism (no wall-clock timestamps in compile-time codegen, no `OnceCell::new()` with random seeds) needs verification. Add this to Spike A3's measurement loop.
- **Validation matrix for `--output-format fastq` explicit-default.** The plan should assert byte-identity between `--output-format fastq` (explicit) and no flag (implicit default). It's the kind of off-by-one regression CI is built to catch.

---

## Alternatives worth considering

### A1 — Two writer traits, not one

Given C1, consider splitting the abstraction:

- `RecordSink` (record-level, used by serial path + specialty + `MimIndexWriter` wrapping): `write_record(&FastqRecord)`.
- `ChunkSink` (chunk-level, used by `parallel.rs` workers): `encode_chunk(&[FastqRecord]) -> Result<Vec<u8>>` + a main-thread `write_chunk_bytes(&[u8]) -> Result<()>`.

The factory in §6 step 6 returns *both* — they're constructed from the same underlying writer state. This is less elegant on paper but matches the actual runtime architecture.

### A2 — Keep mim out of Phase 1

Given C3, the mim format spec gap, and the spike's lack of a success criterion, mim-sidecar is the riskiest piece of Phase 1. A reasonable response: ship Phase 1 without mim (uBAM-in + uBAM-out + BINSEQ-in + BINSEQ-out), and make mim a follow-up minor once Spike A2 has a clear outcome.

This is more conservative than the plan, but it lets BINSEQ + uBAM-out land cleanly without being gated on the mim unknowns. The plan would then be: Phase 1a (uBAM-out + BINSEQ in/out), Phase 1b (mim, gated on spike).

### A3 — `--output-format ubam` + interleaved output by default (O2 above)

Already covered above. Worth restating: emit interleaved BAM by default, FASTQ-symmetric two-file output as opt-in. Matches ecosystem; preserves round-trip; simpler for downstream.

---

## Action items (prioritized)

### Critical (block implementation until resolved)

1. **Resolve C1: `RecordSink`-vs-`parallel.rs`-architecture mismatch.** Pick between Options A / B / C (or A1's two-trait split). Spike this *before* committing to the trait shape in §4.
2. **Resolve C3: define mim format spec + Spike A2 success criterion in concrete terms.** Document what a `.mim` sidecar contains and what offset semantics it requires. Without this, Spike A2 cannot exit.

### Important (resolve in plan v2)

3. **Re-litigate O2 (paired uBAM = two files vs one interleaved).** Plan v2 should justify the choice with reference to ecosystem norms.
4. **Loosen the §3.5 `--preserve-tags` + FASTQ-input + uBAM-out rejection (O1)** to handle batch jobs with mixed-input.
5. **Add `--demux` × output-format rejection rules (O3) to §3.4.**
6. **Specialty-mode integration (O4): either flag as out-of-scope for Phase 1 or add explicit naming + semantics design for each of the 4 modes × 3 formats.**
7. **Clarify Spike A3 consequence (C4): what specifically happens if binary delta > 2 MB?**
8. **Add the validation gaps listed above** — stat parity per output format, `--output-format fastq` explicit-default round-trip, reproducibility on new crates.

### Optional (worth thinking about)

9. **Consider A1's two-trait split** if C1's resolution is messy.
10. **Consider A2 — defer mim to a follow-up minor** if Spike A2 turns out tricky.
11. **Resolve O5 (per-record vs per-chunk trait dispatch overhead)** with a measurement.
12. **Document O6 (per-pair BAM header source) and O7 (interleaved round-trip break).**

---

## Bottom line

The plan is conceptually sound — the symmetric pluggable I/O architecture is the right shape for this problem space. But two critical architectural questions need resolution before implementation:

- **C1** (`RecordSink`-vs-`parallel.rs` architecture mismatch) is the load-bearing one. The plan's "pure refactor, 296 tests pass unchanged" framing is incorrect. Workers do their own gzip compression and ship raw bytes back; there is no writer-trait seam to extend on the parallel side.
- **C3** (mim spec gap + Spike A2 success criterion) makes a third of Phase 1 unmeasurable until clarified.

The other open findings (O1–O8) are tractable scope/semantics decisions but should be made explicit in the plan v2.

Recommend the plan be revised to either (a) resolve C1 with one of the architectural options above, or (b) split Phase 1 into two: a `RecordSink` foundation phase that ships uBAM-out + BINSEQ in/out with the resolved seam, and a follow-up mim phase once the format spec is concrete.
