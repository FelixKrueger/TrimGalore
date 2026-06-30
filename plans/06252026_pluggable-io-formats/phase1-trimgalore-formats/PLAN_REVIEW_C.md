# PLAN_REVIEW_C — Phase 1 (TrimGalore uBAM output, v2)

Reviewer C. Lens: **architectural soundness of v2's asymmetric design** — stress-testing the claim "FASTQ and uBAM outputs live on different code paths, no abstraction binds them."

Reviewed plan: `plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN.md` (v2, 2026-06-25).

---

## Summary verdict

The v2 architectural pivot (drop `RecordSink`, keep `parallel.rs` untouched, route uBAM to a serial branch with a concrete `BamWriter`) is sound and addresses the v1 reviewer findings. But the v2 plan has **one critical contradiction with the existing #317 codebase** (the `BamWriter::write_record(source_aux: &Data)` signature can't actually be threaded through the existing `RecordSource`/`FastqRecord` pipeline) and **three important integrity gaps** (where validation actually fires, fixture/`@PG`-bump fragility detail, and §3.4 completeness). Recommend addressing before implementation.

---

## Logic / architecture

### CRITICAL — `BamWriter::write_record(source_aux: &Data)` is a layer violation that the v2 pipeline cannot satisfy

§4 specifies:

```rust
pub fn write_record(
    &mut self,
    record: &FastqRecord,
    paired_side: Option<u8>,
    source_aux: Option<&noodles::sam::record::Data>,
) -> Result<()>;
```

The `source_aux: Option<&noodles::sam::record::Data>` parameter requires the caller (the trimmer entry point) to have access to the **original `noodles::sam::record::Data`** of each input record.

But in the #317 pipeline, that data is **already lost** by the time the trimmer sees a record:

- `BamReader::next_record` (`src/bam.rs:507`) converts the BAM record into a `FastqRecord` via `bam_record_to_fastq`. The aux data is rendered into a textual `\tTAG:TYPE:VALUE` tail on `FastqRecord.id` (`src/bam.rs:534-558`). The `noodles::sam::record::Data` value is dropped after this point.
- `RecordSource::next_record` (in both `BamReader` and `FastqReader`) returns `Result<Option<FastqRecord>>` — there is no surface to carry the original `Data` through.
- `format::open_sync_reader` returns `Box<dyn RecordSource>` — even the static type erases BAM-ness.

The trimmer therefore has only the `FastqRecord.id` string. To make §4 work, the `BamWriter` would have to **re-parse** the aux tail out of the `FastqRecord.id` string (the inverse of `bam.rs:534-558`'s textual encoding) — which is a brittle, lossy round-trip that the plan does not describe.

**This is more than a comment-style point.** Section 5 step 3 says "Bodies are mostly copies of the FASTQ-writer entry points; the only difference is the writer type and the per-record `paired_side` argument" — but if the only difference is the writer type, then the trimmer has no way to thread `source_aux` through, because the FASTQ-writer code path never had it.

**Fix options** (pick one and write it into the plan):

1. **Parse tags out of the FASTQ-side id tail.** `BamWriter::write_record` takes only `&FastqRecord` and `paired_side`. Internally it splits `record.id` on `\t`, parses the `TAG:TYPE:VALUE` triples, and emits them as aux fields. The textual round-trip is by-definition lossless for the types `append_tag_type_and_value` already produces — but the round-trip must be explicitly designed and tested, especially around `Z:` values that contain whitespace (unlikely but possible).
2. **Extend the `RecordSource` surface** so BAM input can return aux alongside the FastqRecord — e.g. an associated `extra: Option<Data>`. This is heavier (touches `parallel.rs`-adjacent types) but cleaner.
3. **Stash `Data` in a side channel** keyed by record name — fragile and reorder-sensitive; do not do this.

Option 1 is the lowest-touch fix and aligns with §3.3 step 5's note ("Tag values are propagated unchanged (not re-encoded)") — but the plan needs to spell out the parse function and add a unit test for the textual-tag round-trip (e.g. CB/UB/RX through writer→reader→writer should be bytewise stable).

**Severity:** Critical. Without resolving this, the §4 signature is not implementable, and step 3 ("bodies are mostly copies") understates the actual work.

### CRITICAL — §3.4 preserve-tags rule cannot fire at `Cli::validate()`

§3.4 row 5: rejects `--preserve-tags X,Y + ALL inputs FASTQ + --output-format ubam` "at `Cli::validate()`" (§3.4's framing line: *"Excluded combinations (rejected with clear errors at `Cli::validate()`)"*).

But determining "ALL inputs are FASTQ" requires **content-based format detection** (`detect_input_format`), which is per-file file I/O. `Cli::validate()` today does not do format detection — see `src/cli.rs:483` — that work happens in `main.rs` after `Cli::validate()` returns (`src/main.rs:130-148`).

So this rejection cannot live in `Cli::validate()`. It has to be raised in `main.rs`, after the existing `input_formats: Vec<InputFormat>` block at `src/main.rs:134-138`. The other four rejections (clumpify, passthrough, clock, implicon + ubam) ARE pure-CLI and DO belong in `Cli::validate()`.

**Fix:** Plan §3.4 should split: four rejections at `Cli::validate()` + one rejection at the `main.rs` format-detection block (next to the existing `--passthrough + any uBAM` check at `src/main.rs:149-154`). The implementation outline (§5 step 1) should make the same split.

**Severity:** Critical for correctness; small effort to fix.

### IMPORTANT — `run_ubam_output` will duplicate `run_paired_ubam_single_file` and `run_single_file` more than the plan acknowledges

§5 step 4 says "Single-end loop: …; Paired-end (2 FASTQ files): …; Paired-uBAM-single-file (1 BAM file + --paired): open via BamReader::open_paired_interleaved, then same shape." But:

- `run_single_file` (`src/main.rs:638-822`) is **184 lines** including stats printing, report writing (text + JSON), FastQC, and demux dispatch.
- `run_paired_ubam_single_file` (`src/main.rs:1183-1405`) is **222 lines**, with its own variant of stats printing + R1/R2 report writing + paired stats footer.
- `run_paired` (referenced at `src/main.rs:825+`) is a third variant of the same shape.

Reviewer A in v1 flagged drift risk on `run_paired_ubam_single_file` vs `run_paired`. v2 adds a fourth and fifth variant on top of these (the SE-uBAM-output and PE-uBAM-output bodies inside `run_ubam_output`). The "FastQC, demux dispatch, --retain_unpaired naming, paired stats footer, R1/R2 dual-report" patterns will need to be re-implemented in each.

The plan claims "bodies are mostly copies of the FASTQ-writer entry points" (§5 step 3), but in practice it's 4 nearly-identical 200-line functions instead of 2 — drift is realistic and silent (no fixture catches "wrong report header" by design unless we assert on it).

**Recommendation:** Before implementation, factor out the report-writing + summary-printing tails into shared helpers (one for SE, one for PE). This is a single small refactor that pays for itself across 4 call sites and makes the new uBAM branches genuinely thin. The plan should explicitly call this out under §5 step 4 (currently §5 step 3 lists it as "Refactor opportunity (NOT required for v1): factor the trim+filter loop body into a shared inner function. Defer if it complicates review." — but the loop body is the easy part; the reports + summary are the duplication tax).

**Severity:** Important — affects maintainability, not correctness.

### IMPORTANT — Missing-qual detection in §3.3 step 4 conflates two cases

§3.3 step 4: "subtract 33 from each FastqRecord.qual byte (reverse of `BamReader::next_record`'s +33 conversion). Missing-qual handling: if all bytes would be 0xFF (Phred 0 placeholder marker via the `BAM_FUNMAP` + missing convention), emit the BAM sentinel. Otherwise emit raw Phred bytes."

Two issues:

1. **The detection is wrong.** `FastqRecord.qual` is always ASCII-Phred+33. After subtracting 33, "all bytes would be 0xFF" would mean every input qual byte is `33 + 0xFF = 0x120` — impossible for a u8. The intended check is presumably "all subtracted bytes are 0" (i.e. all input quals are `!` = ASCII 33), or "the qual is the synthetic `!` × seq_len that the input-side reader emits for the BAM missing-qual sentinel" (`src/bam.rs:589-593`). The plan's wording doesn't pick a side.
2. **The round-trip is one-way-lossy by design and the plan doesn't acknowledge it.** Input-side: `qual_raw == [0xFF; n]` OR `qual_raw.is_empty()` → emit `!` × seq_len (`src/bam.rs:591-593`). On the output side, `!` × seq_len is indistinguishable from "input read genuinely had Phred 0 across the board" (low-quality but legitimate). The plan needs to decide: do we round-trip the missing-qual sentinel (i.e. detect `!`-only and emit `[0xFF; n]` on output)? Or do we always emit raw subtracted bytes (so missing-qual becomes Phred-0)?

The "right" choice depends on what samtools/Picard do on `samtools fastq input.bam | samtools import - > output.bam`. Worth confirming empirically — both choices are defensible but they must be picked and tested. The current text is ambiguous and would produce bugs.

**Recommendation:** Replace §3.3 step 4 with concrete code. Pseudocode:

```text
phred_bytes = record.qual.bytes().map(|b| b - 33).collect::<Vec<u8>>()
if phred_bytes.iter().all(|&b| b == 0) {
    // All-`!` input — could be missing-qual sentinel OR legitimate Phred 0
    // PLAN MUST DECIDE: emit [0xFF; n] (preserve missing-qual semantics)
    //                   or emit [0u8; n] (lose missing/zero distinction)
}
```

And add a fixture-based round-trip test that pins the chosen semantics.

**Severity:** Important — silent encoding bug if left ambiguous.

---

## Assumptions

### IMPORTANT — §3.4 completeness: mixed-format input batches with `--output-format ubam`

§3.4 doesn't enumerate this combination:

- `--output-format ubam` with multiple inputs where some are FASTQ and some are uBAM.

Currently #317 doesn't allow mixing (the input format is per-file via `detect_input_format`, but the trimming loop is single-input-at-a-time). With uBAM-out, the question is: do all outputs share one `@PG` chain, or does each output BAM get its own header derived from its own input's source header? §3.5 ("uBAM input → uBAM output, header propagation") describes single-input semantics; it doesn't say what happens for two inputs where one is FASTQ and one is uBAM (each gets a separate output BAM — and the FASTQ-input one synthesizes a header per §5 step 2, the uBAM-input one propagates).

**Recommendation:** Add to §3.5: "**Multi-input + uBAM-output:** Each input file's output BAM is independent — uBAM-sourced outputs propagate that input's header + append `@PG`; FASTQ-sourced outputs synthesize a minimal header. No cross-input header merging."

Tangentially: `--rrbs` + `--output-format ubam` is not enumerated either. The plan likely intends this to "just work" (different trim behaviour, but output format is orthogonal) — confirm explicitly and add a test, or reject. RRBS is downstream-sensitive enough that silent "trim differs, output BAM looks the same" could surprise users.

**Severity:** Important — documentation gap; would surface as user confusion or surprise output.

### IMPORTANT — `--demux` rejection should be in §3.4, not buried in §7

§7 says "`--demux`: v1 rejection — demux outputs multiple FASTQ files by 3' barcode; mixing with uBAM output is out of scope." But §3.4's rejection table doesn't list `--demux + --output-format ubam`. Same for §5 step 1 ("Cli::validate() adds the 5 exclusions from §3.4"). The plan needs to be consistent: either lift demux into §3.4 (making it 6 exclusions) or remove the §7 sentence.

**Severity:** Important — internal inconsistency would lead to a missed validation check at implementation time.

### OPTIONAL — Specialty-mode rejection timing (verified: it fires early)

The brief asks whether `--clock + --output-format ubam` rejects at `Cli::validate()` (early — before any I/O) or runtime. Per §3.4 framing it's at `Cli::validate()`. Verified by reading `src/cli.rs:483-661` — both `--clock` and `--implicon` are checked unconditionally as flags in `validate()`, and clap parses `--output-format` to a typed enum field before `validate()` runs. So the rejection genuinely is early. **No issue here** — but worth a one-line note in the plan: "`Cli::validate()` runs before any input file is touched, so this rejection fires within milliseconds of CLI parse, before sanity-checks, format detection, header emit, etc."

(The exception is the preserve-tags rule, see Critical #2.)

**Severity:** Optional — confirmation, not action.

---

## Efficiency

### OPTIONAL — Tag pass-through cost

§6 claims "Aux tag propagation: O(#tags) per record. With typical N≤3 preserved tags, negligible." But under fix-option-1 (parse tags out of `FastqRecord.id` tail), the per-record cost is "scan id for `\t`, split, parse each TAG:TYPE:VALUE". For typical N=2-3 tags this is still microseconds — but the plan should re-state §6 once the parse path is decided. If the round-trip is implemented as a `split('\t').skip(1).map(...)` it's fine; if it's a regex it might be a measurable hot-path cost.

**Severity:** Optional — re-validate once §4 signature is fixed.

---

## Validation sufficiency

### IMPORTANT — §9 validation gate: the `@PG`-on-version-bump fixture issue lacks concrete CI shell

§11 "Remaining risks" says: "Mitigation: assert on (header.minus_@PG_lines, records) tuples in CI, not on the full byte sequence; document fixture-regen recipe."

This is fine as an architectural decision but the CI step is non-trivial to implement cleanly with stock samtools. The most likely shell incantation:

```bash
# Compare bodies (records-only)
diff <(samtools view actual.bam) <(samtools view expected.bam)

# Compare header minus @PG lines
diff <(samtools view -H actual.bam | grep -v '^@PG') \
     <(samtools view -H expected.bam | grep -v '^@PG')

# Assert that exactly one new @PG line exists
test "$(samtools view -H actual.bam | grep -c '^@PG.*ID:trim_galore')" -eq 1
```

But this is fragile: if a future release adds a second `@PG` line, the third assertion breaks; if `samtools view` introduces new header lines in normalisation, the second assertion breaks. The plan should pin the samtools version in CI (matrix entry) AND document the exact shell, OR write a small Rust integration test that uses noodles to do the comparison structurally.

**Recommendation:** Add a concrete CI block to §5 step 6 (or §9), and prefer noodles-based structural comparison in `tests/integration_ubam_out.rs` over samtools-pipe-diff in the CI YAML — the latter is too easy to get wrong silently.

**Severity:** Important — without concrete CI commands, the §9 "validation gate" is aspirational, not landed.

### IMPORTANT — Reproducibility job risk

§6 says "FASTQ-output path: unchanged. The v2 architectural pivot was specifically chosen to avoid touching `parallel.rs`. Performance characteristics on the FASTQ path are bit-identical to today." — true for FASTQ output. But the CI `reproducibility` job builds the **binary** and checks it's byte-identical with the same `SOURCE_DATE_EPOCH`. Adding `OutputFormat` enum + `BamWriter` struct + new `run_ubam_output` function changes the binary. The plan asserts (§9: "Reproducibility intact — Existing CI `reproducibility` job runs unchanged (no new deps) — Bit-identical with same SOURCE_DATE_EPOCH") — that should still hold because no new wall-clock or absolute-path inputs are introduced, but the plan should add a one-line note that the `@PG` line's `VN:<version>` is `env!("CARGO_PKG_VERSION")` (compile-time constant, reproducible) and the `CL:<command-line>` is provided per-invocation (not baked into the binary). Confirm explicitly so the reviewer doesn't have to re-derive this.

**Severity:** Important — small documentation addition prevents a future "did we just break reproducibility?" investigation.

---

## Alternatives

### OPTIONAL — Consider an intermediate `BamRecordLite` type

The Critical #1 fix-option-1 (parse the `id` tail) is fine for v1, but exposes a subtlety: BAM input → BAM output has to **textually round-trip aux tags through `FastqRecord.id`**. If a downstream feature ever needs binary BAM-data types (e.g. `B:` array tags, currently rejected by `append_tag_type_and_value` at `src/bam.rs:641-647`), this design will need to be extended.

An alternative: introduce a small `BamRecordExtras { aux: Option<Data> }` and add a `next_record_with_extras` method to `BamReader` (FASTQ reader returns `None`). Trimmer code can pass extras alongside the FastqRecord. This is more code now but a cleaner foundation for the deferred Phase 4 (BINSEQ-out / mim-out, both of which will face the same "carry source metadata to output" question).

Not a Critical for v1 — option 1 works — but worth mentioning so the team knows option 2 exists.

**Severity:** Optional — architectural foresight.

---

## Action items

### Critical (block implementation)

1. **C1.** Fix the §4 `BamWriter::write_record(source_aux: Option<&Data>)` signature — either redesign to take `&FastqRecord` only and parse tags from the id tail (with explicit round-trip semantics + unit test), or extend the `RecordSource` trait to carry aux through. Pick one in writing, update §4 and §5 step 2, and add a corresponding test under §5 step 5.

2. **C2.** Split §3.4's rejection rules into two groups: **CLI-level** (clumpify, passthrough, clock, implicon — fire in `Cli::validate()`) and **format-detection-level** (preserve-tags-with-no-uBAM-source — fires in `main.rs` after `detect_input_format`, alongside the existing `passthrough + any_bam` check at `src/main.rs:149-154`). Update §5 step 1 to match.

### Important (resolve before implementation lands)

3. **I1.** Plan the dispatch-function refactor in §5 step 4: factor SE summary+report+JSON+FastQC into a shared helper, and the PE equivalent. Don't add a 4th and 5th 200-line variant of `run_single_file`/`run_paired_ubam_single_file`.

4. **I2.** Rewrite §3.3 step 4 with concrete pseudocode, decide the missing-qual round-trip semantics (preserve `!`-only-as-missing? or always emit raw Phred?), and add a fixture-based round-trip test.

5. **I3.** Add to §3.4: `--demux + --output-format ubam` rejection (currently only in §7). Add to §3.5: multi-input behaviour with mixed FASTQ + uBAM sources. Confirm `--rrbs + --output-format ubam` is allowed (with a one-line test).

6. **I4.** Add concrete CI shell (or noodles-based integration test) for the "header minus @PG + records" comparison. Pin the samtools version if using `samtools view` in CI.

7. **I5.** Add a one-line note to §6 confirming reproducibility (`VN:` is compile-time, `CL:` is per-invocation, no wall-clock).

### Optional (nice to have)

8. **O1.** Add a sentence to §3.4 (or §5 step 1) explicitly noting that `Cli::validate()` runs pre-I/O, so the CLI-level rejections fire within milliseconds of parse.

9. **O2.** Mention the `BamRecordExtras { aux: Option<Data> }` alternative under §10 "Open" so the deferred Phase 4 design has prior art to consider.

---

## What I'm NOT flagging (per scope)

- **Completeness vs v1 reviews:** Reviewer D's lane.
- **BINSEQ / mim re-inclusion:** Locked decision.
- **Trait-drop decision itself:** Load-bearing v2 change; both prior reviewers agreed.
