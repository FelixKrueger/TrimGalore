# PLAN — Phase 1: TrimGalore uBAM output (re-scoped, v2.1)

**Epic:** `06252026_pluggable-io-formats/EPIC.md`, Phase 1 — TrimGalore uBAM output.
**Depends on:** — (foundational). uBAM input from [TrimGalore#317](https://github.com/FelixKrueger/TrimGalore/pull/317) is already landed and counts as the first deliverable of this phase.

## 0 · Revision history at a glance

- **v1** — proposed `RecordSink` trait + BINSEQ + mim. Dual review (`PLAN_REVIEW_A.md` + `PLAN_REVIEW_B.md`) caught 5 critical issues that invalidated the architecture.
- **v2** — re-scoped to drop the trait, defer BINSEQ + mim, settle paired-uBAM as one interleaved file. Dual review (`PLAN_REVIEW_C.md` + `PLAN_REVIEW_D.md`) verified all v1 criticals were fixed, but flagged 2 new criticals (tag-source plumbing + validate-time vs format-detection-time rule split) plus 5 important specification gaps.
- **v2.1 (this revision)** — concrete decisions on the v2-review findings. The architectural pivot from v1→v2 (no trait, FASTQ-path untouched, serial uBAM branch) is preserved. v2.1 fills in specification gaps.

### v2.1 changes from v2

| # | v2 issue | Source | v2.1 resolution |
|---|---|---|---|
| 1 | `BamWriter::write_record(source_aux: &Data)` was unsatisfiable — the post-#317 `RecordSource` returns `FastqRecord` with tags folded textually into `id`, not structured `Data` | C-C1 + D-#3 (agreed CRITICAL) | **Pick (a): parse tags from `FastqRecord.id` tail.** §3.6 added documenting the tag round-trip. BamWriter API simplified — no `source_aux` parameter; tag values are recovered from the textual `\tCB:Z:...\tUB:Z:...` form #317 emits. Limitation: only Z/i/f/A scalar tag types supported (B/H arrays already rejected in #317). |
| 2 | `--preserve-tags` "ALL inputs FASTQ" check cannot fire at `Cli::validate()` (needs file I/O) | C-C2 CRITICAL | **§3.4 split into 3.4a (CLI-validate-time, 4 rules) and 3.4b (format-detection-time, 1 rule).** Rejection layer is now explicit. |
| 3 | Code duplication understated — `run_single_file` (184 lines) + `run_paired_ubam_single_file` (222 lines) + new `run_ubam_output` is 4th/5th near-copy | C-I1 | **§5 step 0 added** — extract `write_paired_reports` helper from `run_paired` and `run_paired_ubam_single_file` FIRST, before adding any new dispatch variant. Addresses A-O3 carry-over from v1 review. |
| 4 | Missing-qual wording wrong — "all bytes 0xFF after subtracting 33" is impossible | C-I2 | **§3.3 step 4 rewritten.** Decision: preserve the missing-qual sentinel (uBAM-in with missing qual → uBAM-out with missing qual). FASTQ-in with `!`-only qual → uBAM-out with `!` qual as actual Phred 0 (since FASTQ has no missing-qual sentinel). |
| 5 | §3.4 incomplete — missing `--demux`, `--rrbs`, mixed-format multi-input | C-I3 + D-#5 | **3 new rejection rows added** to §3.4a/b. |
| 6 | `@PG`-fixture-regen mitigation didn't name the helper or show CI shell | C-I4 + D-#4 | **§5 step 5 specifies `assert_ubam_eq` helper** (Rust integration-test fn using noodles, NOT a shell pipe-diff). §9 row updated with the helper reference. |
| 7 | EPIC.md still says "trait stability" for Phase 2 dependency; Phase 4 gates subjective | D-#1 + D-#2 | **EPIC.md updated** (separate change — see EPIC revision history). |

## 1 · Goal

Add **opt-in uBAM output** (`--output-format ubam`) to TrimGalore. After this phase lands:

- **Inputs accepted (unchanged from #317):** FASTQ (plain / gzip / BGZF-framed), uBAM.
- **Outputs available:**
  - FASTQ (default — byte-identical to today, no change to parallel.rs).
  - **uBAM (opt-in via `--output-format ubam`)** — new in this phase. Carries aux tags from input uBAM into output records (via the existing `--preserve-tags` flag from #317), or empty aux when input is FASTQ.
- **Architecture:** NO trait abstraction. FASTQ-output stays on the existing parallel.rs path unchanged. uBAM-output is a **serial** dispatch branch with a concrete `BamWriter`. Future formats (BINSEQ, mim, etc.) get their own branches when added.

Out of scope (deferred to a future phase): BINSEQ in/out, mim sidecar output, opportunistic mim consumption.

## 2 · Context

### Where this sits in the codebase

- `src/bam.rs` — already has `BamReader` from #317. This phase adds `BamWriter` as a sibling.
- `src/main.rs` — dispatch. Adds branching on `cli.output_format`: FASTQ → existing path (parallel or serial); uBAM → new serial-only branch.
- `src/cli.rs` — adds `--output-format <fastq|ubam>` enum flag + validation rules.
- `src/parallel.rs` — **unchanged.** This is the load-bearing v2 simplification — FASTQ output keeps the per-worker-compress + concatenate model that exists today.
- `src/trimmer.rs` — small extension: `run_single_end_to_bam` and `run_paired_end_to_bam` entry points that take a `BamWriter` instead of a `FastqWriter`. The trim logic itself is unchanged.
- `Cargo.toml` — no new deps (noodles 0.88 already pulled).

### Files this plan touches

- Modified: `src/bam.rs` (add `BamWriter`), `src/cli.rs` (flag + validation), `src/main.rs` (dispatch), `src/trimmer.rs` (BamWriter-typed entry points), `.github/workflows/ci.yml` (CI extension).
- New: golden-reference uBAM-output fixtures under `test_files/`, integration tests in `tests/integration_ubam_out.rs`.

## 3 · Behavior

### 3.1 — `--output-format` CLI flag

| CLI | Behavior |
|---|---|
| (no `--output-format`) | FASTQ output mirroring input compression (existing behaviour preserved — `--dont_gzip` overrides) |
| `--output-format fastq` | Same as default; provided as an explicit-default option for scripts that want to be explicit |
| `--output-format ubam` | uBAM output (`*.bam`). Always single-threaded (see §3.4). `--compression`, `--clumpify`, `--dont_gzip` ignored with stderr note |

### 3.2 — Output naming

| Output | SE example | PE example |
|---|---|---|
| FASTQ (default) | `sample_trimmed.fq.gz` | `sample_val_1.fq.gz` + `sample_val_2.fq.gz` |
| **uBAM (new)** | `sample_trimmed.bam` | **`sample_val.bam`** — ONE interleaved BAM with FREAD1/FREAD2 flag bits per record (mate-adjacent ordering — matches samtools/Picard/fgbio/CellRanger convention and #317's `open_paired_interleaved` reader expectation) |

Specialty modes:
- `--hardtrim5 N` + uBAM-out: `sample_<N>bp_trimmed.bam` (single-end shape from the existing `--hardtrim5` naming).
- `--hardtrim3 N` + uBAM-out: `sample_<N>bp_trimmed.bam`.

### 3.3 — Per-record conversion (FastqRecord → bam::Record)

For each trimmed FastqRecord on the uBAM-output path:

1. **Name:** split `FastqRecord.id` on the first TAB or SPACE — everything before is the BAM record `name` (after stripping the leading `@`). Everything after is either preserved-tag tail (§3.6) or descriptive metadata that's discarded.
2. **Flag:** `0x4` (`FUNMAP`) for single-end. For paired-end:
   - R1: `0x1 | 0x4 | 0x8 | 0x40` (PAIRED + UNMAPPED + MATE_UNMAPPED + READ1) = 0x4D = 77.
   - R2: `0x1 | 0x4 | 0x8 | 0x80` = 0x8D = 141.
   - These match what `samtools import -1 r1.fq -2 r2.fq` produces (verified empirically in #317's fixture build).
3. **Sequence:** packed via noodles' 4-bit encoder. Validation mirrors the input direction — ACGTN only; reject IUPAC and `=`.
4. **Quality:** convert each FastqRecord.qual byte by subtracting 33 (reverse of `BamReader::next_record`'s +33). **Missing-qual handling:**
   - If the FastqRecord originated from a uBAM input record with missing qual, `#317`'s `BamReader` emits ASCII `'!'` × seq_len for the FASTQ qual (Phred 0). When writing BACK to uBAM, the BamWriter cannot distinguish "input was missing-qual sentinel" from "input was real Phred 0" — the textual encoding is lossy in this exact case.
   - **Decision: emit raw Phred bytes always.** A read that was missing-qual on input becomes a real Phred-0-everywhere read on output. This is a one-way information loss documented at the §6.3 level — `BamReader` adds `'!'` for missing qual on read; `BamWriter` interprets `!` as Phred 0 on write. Equivalent to how `samtools fastq | samtools import` would behave.
   - Stderr note at startup when uBAM-in → uBAM-out is detected: `"NOTE: missing-qual sentinel in input BAM is converted to Phred 0 through the FASTQ intermediate. Records with no original qual will emit '!' × seq_len rather than the BAM 0xFF sentinel."`
5. **Aux data:** populated via the §3.6 tag round-trip mechanism. Non-empty only when ALL of:
   - `--preserve-tags X,Y,...` was specified
   - The FastqRecord.id carries the textual tag tail produced by `BamReader::next_record` (i.e. input was uBAM)
   - At least one requested tag is found in the textual tail
   For specialty modes that encode metadata in FASTQ headers (`--clock` UMI, `--implicon` UMI), v1 rejects the combination — see §3.4a.

### 3.6 — Aux tag round-trip mechanism (NEW)

Tags flow through TrimGalore via a textual encoding in `FastqRecord.id`. This is asymmetric on the type axis but symmetric on the content axis:

| Direction | Mechanism |
|---|---|
| **Read (#317):** BAM aux → FastqRecord.id | `BamReader::bam_record_to_fastq` appends `\t{TAG}:{TYPE}:{VALUE}` to id for each tag in `--preserve-tags` order (per #317 `append_tag_type_and_value`). Only A/Z/i/f scalar types emitted; B/H rejected. |
| **Write (v2.1):** FastqRecord.id → BAM aux | `BamWriter::write_record` parses `\t{TAG}:{TYPE}:{VALUE}` from the id tail and reconstructs typed `bam::record::Data` fields. Only A/Z/i/f scalar types parsed; anything else (B/H arrays, malformed) errors. |

Parser shape:

1. Strip leading `@` from `FastqRecord.id`.
2. Split on first `\t` — everything before is the record name; everything after is the tag tail (may be empty).
3. Split the tag tail on `\t` to get individual fields.
4. For each field, match the regex `^([A-Za-z][A-Za-z0-9]):([AZif]):(.*)$`. Capture is (tag-name, type-code, value).
5. Construct a `noodles::sam::record::data::field::Value`:
   - `A`: `Value::Character(value.bytes().next().unwrap())`
   - `Z`: `Value::String(value.as_bytes())`
   - `i`: `Value::Int32(value.parse::<i32>()?)` (BAM aux integer types are pre-decided by noodles' encoder)
   - `f`: `Value::Float(value.parse::<f32>()?)`
6. Insert into `bam::record::Data` keyed by the 2-byte tag.

**This round-trip is lossy in two cases — both documented:**
- BAM `B:i,1,2,3` array tags are stripped at read time (#317 already rejects them on the read path). They don't round-trip.
- BAM `i:` tags that were originally `int8`/`int16`/`int64` are re-emitted as `int32` on write. Bit width may change. Decompressed value is byte-identical only for values within `i32` range (which covers all biologically meaningful tag values).

**This shape is preserved unchanged across `--preserve-tags` orderings:** the user-specified order in the textual tail is the order written back to BAM. If the input BAM had tags in a different order, the textual tail (which #317 emits in user-specified order) is the source of truth on output.

### 3.4a — Excluded combinations rejected at `Cli::validate()` (CLI-level, no file I/O)

These rules check only the parsed CLI structure — no input format detection needed. They fire BEFORE `sanity_check_any` runs.

| Combination | Reason | Error message |
|---|---|---|
| `--output-format ubam` + `--clumpify` | clumpify reorders reads to improve gzip compression; meaningless for BAM | `"--clumpify is for gzip output; not applicable with --output-format ubam"` |
| `--output-format ubam` + `--passthrough` | passthrough writes a parallel FASTQ stream; mixed-format pipelines aren't supported in v1 | `"--passthrough is not supported with --output-format ubam in v1"` |
| `--output-format ubam` + `--clock` | --clock encodes UMI metadata in FASTQ headers; translation to BAM aux tags is undefined in v1 | `"--clock + --output-format ubam: UMI-to-BAM-tag mapping not defined in v1; use FASTQ output or convert after"` |
| `--output-format ubam` + `--implicon` | Same as --clock — UMI encoding undefined for BAM | `"--implicon + --output-format ubam: UMI-to-BAM-tag mapping not defined in v1"` |
| `--output-format ubam` + `--demux` | demux emits multiple FASTQ files by 3' barcode; multi-output BAM is out of scope for v1 | `"--demux is not supported with --output-format ubam in v1"` |
| `--output-format ubam` + `--hardtrim5` or `--hardtrim3` | Allowed — see §3.2 naming. No rule here. | — |
| `--output-format ubam` + `--rrbs` (or `--non_directional`) | Allowed — RRBS is a trim-behaviour flag, not an output-format flag. Trim happens at FastqRecord layer; BamWriter just writes whatever the trimmer emits. | — |

### 3.4b — Excluded combinations rejected after format detection (in `main.rs`)

These require knowing whether inputs are FASTQ or uBAM, which only happens once `detect_input_format` has run on each `cli.input[i]`. Fire in `main.rs` between the existing block at line 130-148 (post-#317) and the dispatch branch.

| Combination | Reason | Error message |
|---|---|---|
| `--preserve-tags X,Y,...` + ALL inputs detected as FASTQ + `--output-format ubam` | No source tags to preserve | `"--preserve-tags has no effect with all-FASTQ inputs; either remove the flag or convert at least one input to uBAM via 'samtools import'"` |

Note on the v1→v2 loosening (A-O1): mixed-input batches (some FASTQ, some uBAM) are LEGITIMATE — `--preserve-tags` works per-record, so uBAM-source records get their tags preserved and FASTQ-source records get empty aux. This rule only fires when NO input could conceivably produce a tag. Tested by `tests/integration_ubam_out.rs::preserve_tags_mixed_batch_allowed`.

### 3.5 — Edge cases

- **uBAM output from a FASTQ input** (no source tags): proceeds normally; aux data is empty per record. Stderr note at startup: `"input is FASTQ; output uBAM records will have no aux fields"`.
- **`--output-format ubam` with `--cores N > 1`:** stderr warning then proceed serially. `"--output-format ubam uses single-threaded compression in v1; --cores N is ignored. For high-throughput uBAM output, run multiple invocations in parallel."`
- **uBAM input → uBAM output, header propagation:** the input's SAM header (the `@HD`, `@PG`, `@CO` lines — never `@SQ` since uBAM is unaligned) is propagated to the output, with an additional `@PG ID:trim_galore VN:<version> PP:<previous-PG-id> CL:<command-line>` line appended. Output is NOT byte-identical to the input — provenance is preserved by adding to history, not preserving silence.
- **Paired-uBAM input → paired-uBAM output:** input's `open_paired_interleaved` produces ordered (R1, R2) pairs; output emits them as a single interleaved BAM in the same order. Mate-adjacent.
- **Empty input → uBAM output:** existing input sanity check fires first ("input is empty"); never reaches the BamWriter.

## 4 · Signatures

```rust
// src/bam.rs — NEW BamWriter (sibling to existing BamReader)

pub struct BamWriter {
    inner: noodles::bam::io::Writer<bgzf::Writer<BufWriter<File>>>,
    preserve_tags: Vec<String>, // user-specified tag list; empty = no aux propagation
}

impl BamWriter {
    /// Open `path` for BAM writing. `source_header` is the input uBAM's SAM
    /// header to propagate (with a trim_galore @PG line appended); pass None
    /// for FASTQ input (a minimal synthesised header is used).
    pub fn create<P: AsRef<Path>>(
        path: P,
        source_header: Option<&noodles::sam::Header>,
        preserve_tags: &[String],
    ) -> Result<Self>;

    /// Write one FastqRecord as a BAM record. Flag is `paired_side`-dependent
    /// (None = single-end FUNMAP; Some(1) = paired R1; Some(2) = paired R2).
    ///
    /// Aux tags are recovered from the textual tail of `record.id` per §3.6
    /// (the encoding `BamReader::bam_record_to_fastq` produces on the read
    /// side). The caller does NOT thread aux data separately — `record.id`
    /// is the single source of truth.
    pub fn write_record(
        &mut self,
        record: &FastqRecord,
        paired_side: Option<u8>,
    ) -> Result<()>;

    /// Flush and finalize the writer. Must be called exactly once.
    pub fn finish(self) -> Result<()>;
}
```

Note that this is a **concrete struct**, NOT a trait. Per the v2 architecture decision, FASTQ and uBAM outputs live on different code paths; no abstraction binds them.

**API simplification from v2 → v2.1:** dropped the `source_aux: Option<&Data>` parameter. The post-#317 `RecordSource` returns `FastqRecord` with tags folded textually into `id`, so the caller never has a typed `Data` to pass anyway. Recovering tags from `id` (per §3.6) is the correct layer.

## 5 · Implementation outline

### Step 0 — Refactor `write_paired_reports` helper (pre-work)

Before adding `run_ubam_output` as a 3rd near-copy of `run_paired` / `run_paired_ubam_single_file`, extract the report-generation block into a shared helper. This addresses Reviewer A's v1-finding A3 (which v2 didn't action) and Reviewer C's I1.

1. New function in `src/main.rs`:
   ```rust
   fn write_paired_reports(
       cli: &Cli,
       config: &TrimConfig,
       gzip: bool,
       inputs: PairedInputDescriptors, // helper enum: TwoFastqFiles { r1, r2 } | OneInterleaved { bam }
       output_dir: Option<&Path>,
       stats_r1: &TrimStats,
       stats_r2: &TrimStats,
       pair_stats: &PairValidationStats,
       adapters_r1: &[(String, String)],
       adapters_r2: &[(String, String)],
   ) -> Result<()>;
   ```
2. Both `run_paired` and `run_paired_ubam_single_file` switch to calling this helper (their existing inline ~50-line report blocks become a 3-line call site).
3. `run_ubam_output` (added in Step 4) calls the same helper.
4. **Verification gate:** all 296 existing tests pass unchanged. This step is a pure refactor; no behaviour change.

### Step 1 — `--output-format` CLI flag + validation

1. Add `OutputFormat` enum in `src/cli.rs`: `Fastq`, `UBam`.
2. New `Cli::output_format: OutputFormat` field with clap `value_enum`, default `Fastq`.
3. **`Cli::validate()` adds the 5 rules from §3.4a only** (the 6th — preserve-tags + all-FASTQ — fires at format-detection time in `main.rs` per §3.4b).
4. **`main.rs` between sanity_check_any + dispatch (lines 130-148 area, post-#317)** adds the §3.4b rule. Format-aware check using the already-computed `input_formats: Vec<InputFormat>`.
5. Unit tests in `cli.rs::tests`: each §3.4a rule has a rejection test (mirror the `parse_sam_tag_name` test set from #317).
6. Integration test in `tests/integration_ubam_out.rs::preserve_tags_all_fastq_rejected`: invokes the binary with `--preserve-tags CB --output-format ubam` + a FASTQ-only input, asserts non-zero exit + specific error message.

### Step 2 — `BamWriter` struct

1. New `BamWriter` type in `src/bam.rs`. Wraps `noodles::bam::io::Writer<bgzf::Writer<BufWriter<File>>>`.
2. `BamWriter::create(path, source_header, preserve_tags)`:
   - When `source_header` is `Some`, propagate the input header + append a `trim_galore` `@PG` line.
   - When `None`, synthesise a minimal header with just the `@HD` line + the trim_galore `@PG` line.
3. `BamWriter::write_record(record, paired_side, source_aux)`:
   - Build a `bam::Record` from the FastqRecord (per §3.3 spec).
   - Set flag bits per `paired_side`.
   - Copy aux fields from `source_aux` if `preserve_tags` is non-empty and `source_aux` is `Some`.
4. `BamWriter::finish(self)`: flush + write BGZF EOF marker. Consumes self (Rust idiom for forced finalisation).
5. Unit tests:
   - `bam_writer_se_fastq_input`: synthesise 3 FastqRecords, write SE BAM, read back via `BamReader`, verify (id, seq, qual) tuples match.
   - `bam_writer_pe_interleaved`: synthesise 3 pairs, write interleaved BAM, read back via `BamReader::open_paired_interleaved` (existing from #317), verify pair-alignment.
   - `bam_writer_aux_preservation`: uBAM-in fixture with CB/UB tags → write back with `--preserve-tags CB,UB` → read back → tags match exactly.
   - `bam_writer_flag_bits`: SE record has flag 0x4; paired R1 has 0x4D=77; paired R2 has 0x8D=141.

### Step 3 — Trimmer entry points for uBAM output

1. Add `run_single_end_to_bam(reader: &mut dyn RecordSource, writer: &mut BamWriter, config: &TrimConfig) -> Result<TrimStats>`.
2. Add `run_paired_end_to_bam(reader_r1: ..., reader_r2: ..., writer: &mut BamWriter, ...)`.
3. Bodies are mostly copies of the FASTQ-writer entry points; the only difference is the writer type and the per-record `paired_side` argument.
4. Refactor opportunity (NOT required for v1): factor the trim+filter loop body into a shared inner function. Defer if it complicates review.

### Step 4 — `main.rs` dispatch

1. After format detection + CLI validation, branch on `cli.output_format`:
   ```rust
   match cli.output_format {
       OutputFormat::Fastq => {
           // existing dispatch — parallel for cores > 1 OR clumpify,
           // serial otherwise. Zero changes to this branch.
       }
       OutputFormat::UBam => {
           // NEW — always serial. Warn if cli.cores > 1 ("ignored").
           run_ubam_output(...)?;
       }
   }
   ```
2. New `run_ubam_output(cli, gzip, output_dir)` function:
   - Single-end loop: for each input, open reader (FASTQ or BAM), open BamWriter, call `trimmer::run_single_end_to_bam`.
   - Paired-end (2 FASTQ files): open two readers, open BamWriter, call `trimmer::run_paired_end_to_bam`.
   - Paired-uBAM-single-file (1 BAM file + `--paired`): open via `BamReader::open_paired_interleaved`, then same shape.
3. Output naming via new helpers in `src/io.rs`: `single_end_bam_output_name(path, ...)`, `paired_bam_output_name(path, ...)`.

### Step 5 — Tests + golden fixtures

1. New file `tests/integration_ubam_out.rs`:
   - `ubam_out_se_fastq_input` — FASTQ in → uBAM out, verify samtools view round-trip gives same (id, seq, qual) tuples as the FASTQ-out path.
   - `ubam_out_se_ubam_input` — uBAM in → uBAM out, verify aux tags propagated when `--preserve-tags` set.
   - `ubam_out_pe_two_fastq` — two FASTQ in → interleaved uBAM out, mate-adjacent verified by alternating FREAD1/FREAD2 flag.
   - `ubam_out_pe_one_ubam_interleaved` — one uBAM in (`--paired`) → interleaved uBAM out, round-trip verifies mate adjacency preserved.
   - `ubam_out_hardtrim5` — `--hardtrim5 20` + uBAM out works (specialty mode that's allowed).
   - `ubam_out_clock_rejected` — `--clock` + uBAM out errors clearly.
   - `ubam_out_clumpify_rejected` — `--clumpify` + uBAM out errors clearly.
   - `ubam_out_preserve_tags_no_source_rejected` — FASTQ-only input + `--preserve-tags` + uBAM out errors clearly.
2. Golden fixtures: `test_files/ubam_out_se_REFERENCE.bam` and `test_files/ubam_out_pe_REFERENCE.bam`. Committed.
3. **`assert_ubam_eq` test helper** (NEW — addresses C-I4 + D-#4):

   ```rust
   // tests/integration_ubam_out.rs
   /// Compare two BAM files for content equivalence, ignoring `@PG` lines
   /// (trim_galore's @PG carries the package version, which would otherwise
   /// break the fixture on every routine version bump).
   ///
   /// Uses noodles to parse both files; compares:
   ///   - header lines EXCEPT @PG (full equality)
   ///   - record stream (name, flag, seq, qual, sorted aux fields) tuple-equality
   ///
   /// Tuple-equality avoids transient field-ordering noise in aux Data; sorting
   /// by tag-name on both sides normalises that.
   fn assert_ubam_eq(actual: &Path, expected: &Path);
   ```

   Implementation uses `noodles::bam::io::Reader` for both files. Tuple-compares records via `(name, flags, seq.iter().collect::<Vec<u8>>(), qual.iter().collect::<Vec<u8>>(), sorted_aux_tuples)`. NOT a `samtools view -H | grep -v @PG` shell pipe — that's fragile to whitespace + ordering. The Rust helper is precise and CI-stable.

4. Document the fixture-regen recipe in `test_files/README.md` for the case where committed-byte fixtures DO need refresh (e.g. major header schema change).

### Step 6 — CI validation extension

1. Extend `.github/workflows/ci.yml::validation` job:
   - New step: SE uBAM-out content-tuple parity (samtools view round-trip should produce identical (id, seq, qual) tuples to the existing FASTQ-out path).
   - New step: PE interleaved uBAM-out — `samtools view -c -f 0x40` and `-f 0x80` counts equal; `samtools flagstat` confirms `paired in sequencing` count.
   - New step: aux-tag preservation — `samtools view` greps for the committed `CB:Z:ATCGATCG-1` tag.
2. NO Perl byte-identity check needed (no Perl equivalent for uBAM output).

### Step 7 — Documentation

1. `CLAUDE.md` §Conventions: add a line on uBAM output's serial-only stance and the `@PG`-bumps-on-version-change behaviour.
2. `README.md` Usage section: one new example `# uBAM output (preserve aux tags through trimming):\ntrim_galore --preserve-tags CB,UB --output-format ubam sample.bam`.
3. `CHANGELOG.md`: extend the existing uBAM entry under "Unreleased" with the `--output-format ubam` flag.

## 6 · Efficiency

- **FASTQ-output path: unchanged.** The v2 architectural pivot was specifically chosen to avoid touching `parallel.rs`. Performance characteristics on the FASTQ path are bit-identical to today.
- **uBAM-output path: serial-only.** This is a deliberate v1 trade-off. uBAM-in → uBAM-out is bottlenecked by input BGZF decompression (single-threaded in noodles) regardless of output parallelism, so the loss is small for the dominant use case. For FASTQ-in → uBAM-out workloads with large inputs and many cores, this is slower than it could be — flagged as a known limitation; users who need throughput can post-process the FASTQ output to BAM via `samtools import`.
- **Per-record conversion cost:** noodles' 4-bit-packed seq encoding + Phred subtraction is well under the gzip compression cost in the FASTQ path. BamWriter throughput should be similar to or better than FastqWriter at gzip level 1.
- **Aux tag propagation:** O(#tags) per record. With typical N≤3 preserved tags, negligible.

## 7 · Integration

- **Input path:** unchanged from #317. Format detection + reader factory route through `format::open_*_reader`.
- **Output path:** new `--output-format ubam` branch in main.rs dispatch. FASTQ branch unchanged.
- **Stats / reports:** unchanged. Trimming reports (text + JSON) are emitted the same way regardless of output format.
- **Validation matrix:** existing Perl 0.6.11 byte-identity assertion stays for FASTQ output. New CI steps cover uBAM output (no Perl baseline since the legacy tool doesn't support BAM output).
- **Specialty modes:** `--hardtrim5/3` work with uBAM-out; `--clock` / `--implicon` rejected (UMI-to-aux-tag mapping deferred).
- **`--demux`:** v1 rejection — demux outputs multiple FASTQ files by 3' barcode; mixing with uBAM output is out of scope.

## 8 · Assumptions

### Carried forward from the epic

- Single static binary invariant — no new external runtime deps. uBAM via the already-pulled `noodles` 0.88 umbrella; no new crate needed.
- Byte-identity baseline for FASTQ output to Perl 0.6.11 preserved.
- BINSEQ + mim deferred until crate ecosystem matures (per the §0 v2 revision summary).

### Phase-specific (v2)

- **NO `RecordSink` trait.** FASTQ and uBAM outputs live on different code paths. Future formats add new branches, not new trait impls.
- **uBAM output is always serial in v1.** Parallel uBAM output would require a different writer architecture (likely a single-writer + ordered-channel-from-workers design); deferred until a real workload justifies the complexity.
- **Paired uBAM output is one interleaved file.** Matches samtools/Picard/fgbio/CellRanger convention + #317's reader-side expectation.
- **`@PG` line is appended on uBAM output**, so uBAM-in → uBAM-out is provenance-preserving but NOT byte-identical to input. Documented.

## 9 · Validation

| What | How | Expected |
|---|---|---|
| FASTQ output unchanged | `cargo test --release` + CI Perl 0.6.11 byte-identity assertion | All 296 tests pass; FASTQ-out md5 matches Perl |
| BamWriter round-trip (SE) | Integration test: FASTQ → BamWriter → BamReader → assert tuple parity | Pass on 10-record fixture |
| BamWriter round-trip (PE interleaved) | Integration test: 2 FASTQ → BamWriter → `BamReader::open_paired_interleaved` → assert mate adjacency | Pass on 10-pair fixture |
| Aux tag preservation | uBAM-in-with-tags → uBAM-out with `--preserve-tags CB,UB` → samtools view confirms aux fields | Pass |
| Empty aux on FASTQ-in | FASTQ-in → uBAM-out → samtools view confirms no aux fields | Pass + stderr note |
| Flag bits correct | samtools view confirms flag 0x4 / 0x4D / 0x8D for SE / R1 / R2 | Pass |
| Five rejection rules | Per-rule unit tests in `cli.rs::tests` | Pass |
| `--cores N>1` + uBAM-out warns + still works | Integration test: invoke with --cores 4 + --output-format ubam, assert exit 0 + stderr contains "ignored" | Pass |
| Reproducibility intact | Existing CI `reproducibility` job runs unchanged (no new deps) | Bit-identical with same SOURCE_DATE_EPOCH |
| samtools content parity (the v1 validation gate) | CI: `samtools view sample.bam` from FASTQ-out path + uBAM-out path produce identical (name, seq, qual) triples | Pass — this catches encoding bugs the BamReader-only round-trip couldn't (Reviewer B's Important #5) |
| Committed-golden parity (version-bump-resilient) | Integration test calls `assert_ubam_eq(actual, expected_golden)` — compares header-minus-@PG + record tuples via noodles | Pass; resilient to `@PG VN:` bumps on routine version increments |
| `--preserve-tags` mixed-batch | Integration test: one FASTQ input + one uBAM-with-tags input + `--preserve-tags CB,UB --output-format ubam`; assert FASTQ-side records have empty aux, uBAM-side records have CB+UB tags | Pass — confirms A-O1 loosening works |
| Aux tag type-correctness | uBAM-in with `i` (integer) + `f` (float) tags → round-trip → samtools view confirms types preserved (NOT silently downgraded to `Z` strings) | Pass — covers §3.6 round-trip type fidelity |

## 10 · Questions or ambiguities

| Q | Priority | Disposition |
|---|---|---|
| uBAM output trigger (opt-in vs mirror vs sidecar vs deferred) | Critical | **Resolved 2026-06-25 (v1)** — opt-in via `--output-format ubam`. |
| Paired uBAM output: two BAMs vs one interleaved | Critical | **Resolved 2026-06-25 (v2)** — one interleaved BAM (per Reviewer A-O2 + Reviewer B-Crit2, and matches #317's reader). |
| BINSEQ / mim inclusion in Phase 1 | Critical | **Resolved 2026-06-25 (v2)** — both deferred to future phase pending crate-ecosystem maturity. |
| RecordSink trait | Critical | **Resolved 2026-06-25 (v2)** — dropped. FASTQ and uBAM live on different code paths. |
| Parallel uBAM output | Open | Deferred to follow-up. v1 is serial-only with a stderr note. |
| `--clock` / `--implicon` + uBAM-out UMI translation | Open | v1 rejects the combination; revisit when there's a real workload demanding it. |
| `--demux` + uBAM-out | Open | v1 rejects (multi-output demux + BAM is a different scope). |
| Mate-info fields on paired uBAM output | Open | v1 emits with `RNEXT=*, PNEXT=0, TLEN=0` (per BAM spec for unaligned mate pairs). Could populate richer mate info but uBAM by definition has no positions to point at. |

## 11 · Self-Review

### v2 changes from v1, mapped to reviewer findings

| v1 design | v2 change | Source |
|---|---|---|
| `RecordSink` per-record-write trait | **Dropped.** FASTQ output unchanged. uBAM output is a separate serial code path with concrete `BamWriter` | A-C1 + B-Crit1 (agreed) |
| `parallel.rs::run_*_parallel` accepts `Box<dyn RecordSink>` | **Dropped.** parallel.rs unchanged | A-C1/C2 + B-Crit1 |
| Paired uBAM-out = `sample_val_1.bam` + `sample_val_2.bam` (two BAMs) | **One interleaved `sample_val.bam`** matching samtools/Picard/fgbio/CellRanger convention | A-O2 + B-Crit2 (agreed) |
| BINSEQ input/output (BinseqReader + BinseqWriter + Cargo feature) | **Dropped.** Defer to future phase. binseq crate is unstable + pulls zstd-sys C build.rs | B-Crit4 |
| mim sidecar output (`--emit-mim`) + mim opportunistic input | **Dropped.** Defer to future phase. mim-index requires Rust 1.91; project floor is 1.88 | B-Crit3 + A-C3 |
| Three pre-work spikes (A1 BINSEQ, A2 mim, A3 binary-size) | **Dropped.** All three were measuring deferred features | n/a (consequence of scope cut) |
| `--preserve-tags X,Y` + FASTQ input + `--output-format ubam` = hard error | **Loosened.** Hard-error only if ALL inputs are FASTQ; mixed batches with at least one uBAM are allowed | A-O1 |
| Specialty mode interaction underspecified | **Sharpened.** `--hardtrim5/3` work with uBAM-out; `--clock` / `--implicon` rejected with a clear "UMI-to-aux-tag mapping deferred" error | A-O4 |
| §9 validation = "round-trip through BinseqReader" (self-referential) | **Replaced.** v2 validates via `samtools view` (external reference) — catches encoder bugs that the BamReader round-trip alone wouldn't | B-Important #5 |
| §6 BINSEQ encoding cost claims | **Removed.** BINSEQ deferred | B-Important #6 |
| §5 step 7 feature-gated CLI message | **Removed.** No Cargo features in v2 (no new deps) | B-Important #7 |

### Self-review of v2 itself

- **Efficiency** — FASTQ path is unchanged (zero perf impact). uBAM path is serial in v1; the dominant uBAM workload (uBAM-in → uBAM-out) is bottlenecked by input decompression anyway, so the loss is small. Flagged as known limitation.
- **Logic** — main.rs dispatch branches cleanly between FASTQ (existing) and uBAM (new); no shared abstraction means no risk of cross-format bugs. Specialty mode interactions explicitly enumerated.
- **Edge cases** — §3.5 covers FASTQ-in → uBAM-out (empty aux), uBAM-in → uBAM-out (header propagation + @PG addition), paired interleaved input → interleaved output, empty input.
- **Integration** — Perl byte-identity gate stays intact (no parallel.rs change). New CI assertions cover uBAM-out without competing with the Perl oracle (no Perl uBAM-out exists).

### Remaining risks

- **`@PG`-on-version-bump breaks fixture byte-identity.** Mitigation: assert on (header.minus_@PG_lines, records) tuples in CI, not on the full byte sequence; document fixture-regen recipe.
- **Serial uBAM output may be slow for FASTQ-in → uBAM-out with large inputs.** Documented; users who need throughput can post-process via `samtools import`. Parallel uBAM output is in §10 as a future-work item.
- **Mate-info fields for paired uBAM-out are empty (RNEXT=`*`, PNEXT=0).** Per BAM spec for unaligned mate pairs; documented but worth validating that downstream tools accept this shape (samtools view should; verify in CI).

---

## Revision History

- **v1, 2026-06-25** — initial plan. RecordSink trait, BINSEQ + mim included, three pre-work spikes. Reviewed (PLAN_REVIEW_A.md + PLAN_REVIEW_B.md) — 5 critical findings invalidated the architecture.
- **v2, 2026-06-25** — re-scoped per reviewer findings. RecordSink trait dropped (FASTQ + uBAM on different code paths). BINSEQ + mim deferred to future phase. Paired uBAM-out is now one interleaved BAM (ecosystem convention). Five rejection rules sharpened. Validation moved from self-referential reader-round-trip to samtools-anchored content parity. Reviewed (PLAN_REVIEW_C.md + PLAN_REVIEW_D.md) — verdict PARTIAL with 2 critical + 5 important specification gaps.
- **v2.1, 2026-06-26** — addresses all v2-review findings. Tag-source plumbing resolved by parsing the textual tag tail from `FastqRecord.id` (option (a) — keeps `RecordSource` trait + `FastqRecord` shape stable). New §3.6 documents the round-trip. §3.4 split into CLI-validate-time (3.4a, 5 rules incl. new --demux row) vs format-detection-time (3.4b, 1 rule). Missing-qual semantics decided (raw Phred always; sentinel preservation impossible through FASTQ intermediate). §5 step 0 added: refactor `write_paired_reports` helper before adding any new dispatch variant. §5 step 5 specifies `assert_ubam_eq` helper (noodles-based, NOT samtools pipe-diff). §9 gains four new validation rows. BamWriter API simplified — `source_aux: Option<&Data>` parameter dropped (was unsatisfiable per C-C1).
