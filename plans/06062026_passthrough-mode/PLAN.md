# PLAN — `--passthrough` mode (Multiome/10X cell-barcode carrier)

**Status:** v2, revised 2026-06-06 after dual independent plan-reviewer agent feedback. Awaiting manual re-review → implementation trigger.

---

## Revision history

### v2 (2026-06-06) — incorporates dual plan-reviewer feedback

Source reviews: `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`.

**Consensus items (both reviewers, highest confidence):**
- AB1: Parity-test fixture bumped from 50 records to ≥10,000 (forces ≥3 batches at `BATCH_SIZE = 4096`; exercises BTreeMap ordered-flush).
- AB2: `read_id_prefix()` now strips trailing `/[123]` in v1 (was Q&A #2, deferred to v2). One-line change; covers legacy SRA/ENA headers (`@read/1` / `@read/2` / `@read/3`).
- AB3: Added Validation §6 — `--dont_gzip` plain-output test (plain branch of `process_paired_batch` diverges from gzip branch).
- AB4: Added Validation §7 — `--cores 1` dispatcher-level CLI integration test (direct `run_paired_end` invocation bypasses `main.rs::run_paired` routing).
- AB5: Validation §2 (parity) strengthened — asserts row-by-row R1/R2/passthrough ID-lockstep, not just stat counts and decoded passthrough equality.
- AB6: Worker-thread plumbing for `process_pairs` spelled out: `passthrough: &[FastqRecord]` (immutable; written verbatim) vs R1/R2's `&mut [FastqRecord]` (mutated in-place during trim); gzip/plain branch fork mirrored explicitly.

**Reviewer A — critical bugs caught:**
- A-Crit-1: `PairValidationStats` (`src/report.rs:106`) currently derives only `Debug, Default`. Step 3 now adds `#[derive(PartialEq, Eq)]` as an explicit sub-step. Without this, Validation §2's `assert_eq!(serial_stats, parallel_stats)` does not compile.
- A-Crit-2: EOF arm for passthrough-longer-than-R1/R2 enumerated explicitly in `read_pairs_round_robin` (Step 6.5). Plan v1 only covered this for the serial path.
- A-Other: §Edge cases corrected — `FastqReader::sanity_check` rejects empty files at startup (`fastq.rs:407–438`); there is no "clean exit, zero-records report" path for empty passthrough input.
- A-Other: §Efficiency sync-check cost corrected from "~30 ms / 1M reads" to "~90 ms / 1M reads" (still well under 0.1% overhead; corrected for honesty).

**Reviewer B — critical bugs caught:**
- B-Crit-1: Reader-thread error propagation hardened. Plan v1 silently inherited the existing parallel.rs pattern where `reader_handle.join()` runs **after** the `result_rx.recv()` loop drains — so a mid-stream reader error lands partial output on disk before reaching the caller. New Implementation Step 6a specifies the required fix shape (reader-error visibility to the main thread before write commits).
- B-Crit-2: Assumption #9 softened from "byte-identical" to "record-identical". `FastqReader::next_record` (`fastq.rs:325–355`) discards the `+...` description on line 3; `FastqRecord::write_to` (`fastq.rs:53–70`) re-emits bare `b"+\n"` and canonicalises `\r\n` → `\n`. Records (id/seq/qual) are preserved verbatim; file bytes are canonicalised.
- B-Other: Output gzip mode is driven from R1's compression alone (`main.rs:84`: `is_gzipped(&cli.input[0])`), not the passthrough input's extension. Documented in §Assumptions.
- B-Other: Behavior §7 now enumerates all 8 arms of the three-way `(Some/None)^3` EOF match (1 proceed, 1 exit, 6 bail).
- B-Other: `passthrough_records_checked` accounting site decided: the worker stamps it from `batch.len()` when passthrough is present, on the contract that the reader has already verified per-record sync before dispatching the batch. No new cross-thread counter needed.

**Contradictions resolved (per main-loop judgement):**
- *Stats placement* (A: new `PassthroughStats` struct ↔ B: extend `PairValidationStats`) — kept B's framing. Lifecycle matches (per-pair, gated by passthrough mode); avoids rippling a 4th tuple slot through every caller.
- *Severity of `--passthrough == R1/R2` collision check* (A: out-of-scope, doc-only ↔ B: real logic gap) — adopted B's framing. Validation step 1.ix now requires the same case-folded `norm()` helper the output-collision pre-flight uses (`main.rs:178`), consistent with issue #216 protection.

### v1 (2026-06-06)
Initial draft.

---

## Goal

Add a new CLI flag `--passthrough <FILE>` to Trim Galore. In paired-end mode, it accepts a third input file (typically an Illumina I1/I2 index or a 10X Multiome cell-barcode + UMI read) that is **not trimmed** — instead, its records are carried through 1:1 to a parallel output file, kept in lockstep with R1/R2 such that any pair dropped by quality/adapter/length/N filtering also drops the corresponding passthrough record. The headline driver is the 10X Multiome ATAC-seq use case where Cellranger Arc cannot detect the modified ATAC adapter; users want to trim R1/R2 against their own `-a/-a2` sequence while keeping the I1 barcode file aligned to the survivors.

The outcome is a single invocation like:

```
trim_galore --paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz
```

producing `R1_val_1.fq.gz`, `R2_val_2.fq.gz`, and `I1_passthrough.fq.gz`, all with identical record counts and matching read-ID order.

---

## Context

### Placement

This feature is purely additive — no behavioural change in the absence of the flag. Touch points:

- `src/cli.rs` — new `pub passthrough: Option<PathBuf>` field on `Cli`; new validation block in `Cli::validate()`.
- `src/io.rs` — new `passthrough_output_name()` function next to `paired_end_output_names()`.
- `src/main.rs` — extend the paired-end output-collision pre-flight (lines ~170–220) to include the passthrough output path; thread the passthrough input through to `run_paired()`.
- `src/trimmer.rs` — extend `run_paired_end()` (serial path) to accept an optional third reader/writer + sync-check stats. Add `read_id_prefix()` helper.
- `src/parallel.rs` — extend `run_paired_end_parallel()`, `PairedWork`, `PairedBatchResult`, `process_paired_batch`, `process_pairs`, and `read_pairs_round_robin` to carry the optional third stream. Harden reader-thread error propagation.
- `src/report.rs` — add `#[derive(PartialEq, Eq)]` to `PairValidationStats`; add three new passthrough counter fields + corresponding merge logic; emit passthrough block in text + JSON reports.
- `src/fastqc.rs` — extend the FastQC sweep entry point to include the passthrough output when `--fastqc` is set.

### Existing patterns this aligns with

- **`Cli::validate()` incompatibility blocks**: the `--clumpify` rejection chain (`cli.rs:488–520`) is the model — each incompatible flag becomes one `anyhow::bail!`. We use exactly the same shape.
- **`classify_paired` is sink-agnostic** (`parallel.rs:309`) — it returns `PairOutcome::Pass / Discarded / Unpaired` and lets the caller decide what to write. The passthrough sink slots in as a third immutable slice argument to `process_pairs`, with the existing `Pass / Discarded` branches gated on whether the third sink is `Some`.
- **Output-collision pre-flight** (`main.rs:166–215`): single `HashMap<case_folded_path, PathBuf>` that bails on duplicate. Both the passthrough output candidate and the passthrough-equals-input collision check reuse the same `norm()` closure (`main.rs:179`).
- **Output naming** (`io.rs:80–119`): `paired_end_output_names()` is the template; the new `passthrough_output_name()` follows the same `(input, output_dir, basename, gzip) → PathBuf` shape and reuses `strip_fastq_extensions()`.
- **Report stats accumulation** (`parallel.rs:155–183`): `pending: BTreeMap<u64, PairedBatchResult>` reassembles in order, `merge()` per stats struct. The passthrough stat fields ride along inside `PairValidationStats`.
- **Parity contract test** (`parallel.rs:test_parallel_serial_trim_stats_parity` at line 1008): the new parity test for passthrough mode follows this exact shape — serial vs. parallel produce field-identical stats on the same input — with the added strengthenings called out in Validation §2.

### Compatibility envelope (v1)

**Allowed:**
- `--paired` + any adapter/quality flags (`-a`, `-a2`, `--quality`, `--nextseq`, `--illumina`, `--nextera`, `--small_rna`, `--stranded_illumina`, `--bgiseq`, `--rrbs`, `--non_directional`, `--max_n`, `--max_length`, `--length`, `--clip_R1/R2`, `--three_prime_clip_R1/R2`, `--trim-n`, `--poly_a`, `--poly_g`, `--discard_untrimmed`, `--times`, `--stringency`, `--error_rate`, `--consider_already_trimmed`, `--rename`).
- `--cores N` (single-end-style parallelism via worker pool).
- `--fastqc` and `--fastqc_args` (FastQC runs on all three outputs).
- `--output_dir`, `--basename`, `--dont_gzip`, `--compression <N>`.

**Rejected in v1 (one-line `anyhow::bail!` each in `Cli::validate()`):**
- `--retain_unpaired` — passthrough requires strict pair semantics; a half-rescued mate has no clean interpretation for the index file in v1. (Note: a clean v2 semantic exists — emit the passthrough record on any `r1_ok || r2_ok` survivor — deferred until a user asks.)
- `--clumpify` — the clumpy reorder shuffles records by minimizer, which breaks the lockstep contract. The clumpy dispatcher would need a parallel third stream and matching reordering; deferred.
- Multi-pair input (input.len() > 2 with `--paired --passthrough`) — single triple only for v1.
- `--clock`, `--implicon`, `--hardtrim5`, `--hardtrim3`, `--demux` — specialty modes own their own input arity and output naming; combining them with `--passthrough` is murky and not in the Multiome use case.
- `--paired` is **required** when `--passthrough` is set (no single-end + passthrough — there's no Multiome shape that needs this).

### Why now / why this matters

User-reported: 10X Multiome ATAC-seq with a modified adapter that Cellranger Arc fails to detect. Users currently work around this with awk scripts that filter the I1 file post-hoc; that's error-prone and slow. Trim Galore already does the pair-aware filtering work; carrying a third file is a small extension that turns a multi-tool workflow into a single invocation.

---

## Behavior

Numbered steps for the paired-end pipeline when `--passthrough <PATH>` is set:

1. **Parse and validate.** `clap` parses `--passthrough` into `Cli::passthrough: Option<PathBuf>`. `Cli::validate()` runs the new incompatibility block:
   1. If `passthrough.is_some()` and `!paired` → bail with `"--passthrough requires --paired"`.
   2. If `passthrough.is_some()` and `input.len() != 2` → bail with `"--passthrough requires exactly one R1/R2 pair (got N input files); multi-pair input with passthrough is not yet implemented"`.
   3. If `passthrough.is_some()` and `retain_unpaired` → bail with `"--passthrough is incompatible with --retain_unpaired (passthrough requires strict pair semantics in v1)"`.
   4. If `passthrough.is_some()` and `clumpify` → bail.
   5. If `passthrough.is_some()` and `clock` → bail.
   6. If `passthrough.is_some()` and `implicon.is_some()` → bail.
   7. If `passthrough.is_some()` and `hardtrim5.is_some()` / `hardtrim3.is_some()` / `demux.is_some()` → bail.
   8. If the passthrough file does not exist → bail with `"--passthrough file not found: {path}"`.
   9. **Case-folded collision check** — if `case_fold(passthrough_path) == case_fold(input[0])` or `== case_fold(input[1])` → bail with `"--passthrough cannot point at R1 or R2 (case-insensitive match on APFS/NTFS)"`. Use the same `norm()` closure as the output-collision pre-flight at `main.rs:179`. *(Resolves contradiction #2 from review.)*
2. **Sanity-check the passthrough file.** Call `FastqReader::sanity_check(passthrough_path)` once at startup, alongside the existing R1/R2 sanity check in `main.rs` (line ~265). An empty passthrough file is rejected here, before any output is written.
3. **Compute the passthrough output path.** Use the new `io::passthrough_output_name(input_passthrough, output_dir, basename, gzip)` → e.g. `<stem>_passthrough.fq.gz`. The `--basename` semantics match `_val_{1,2}`: when `--basename foo` is set, all three outputs use the `foo_*` stem. Without `--basename`, the passthrough output stem is derived from the passthrough input filename via `strip_fastq_extensions()`. **Note:** the gzip extension is determined by the global `gzip` boolean (which itself derives from R1's compression at `main.rs:84`), not from the passthrough input's extension. See Assumptions §11.
4. **Extend the output-collision pre-flight.** In `main.rs` lines ~177–214, push the passthrough output path into the per-pair `candidates` vec before the `HashMap::insert` loop.
5. **Open the passthrough reader.** In the parallel path, open a third `FastqReader::open_threaded(passthrough_path)` in the reader thread. In the serial path, open a third `FastqReader::open(passthrough_path)` inline.
6. **Open the passthrough writer.** In both paths, create a third `FastqWriter::create(passthrough_output_path, gzip, compression_level, ...)` (serial) or a third `File::create()` + per-batch `GzEncoder` (parallel).
7. **Lockstep loop with explicit 8-arm three-way EOF match.** For each iteration, read one record from R1, one from R2, and one from the passthrough file, then `match (rec1, rec2, rec3)`:
   1. `(Some(r1), Some(r2), Some(idx))` → **proceed.**
      1. **Header-prefix sync check.** Apply `read_id_prefix()` (definition in §Signature) to all three IDs. If any pair-wise comparison disagrees, bail with `"Read ID mismatch at record N: R1='{id1}', R2='{id2}', passthrough='{idx_id}'. Files are out of sync."` (Names the offending row so the user can grep their input.)
      2. **Record-counted.** Increment `pair_stats.passthrough_records_checked` by 1 in the serial path; in the parallel path the worker stamps `pair_stats.passthrough_records_checked = batch.len()` after the batch is dispatched (the reader has already per-record-verified sync).
      3. Run `classify_paired(r1, r2, …)` (parallel) or its serial equivalent. Return is `PairOutcome::{Pass, Discarded, Unpaired}`.
      4. On `Pass`: write the passthrough record to its writer; `pair_stats.passthrough_records_kept += 1`.
      5. On `Discarded`: drop the passthrough record; `pair_stats.passthrough_records_dropped += 1`.
      6. On `Unpaired`: cannot occur (rejected at validation step 1.iii). `debug_assert!(false, "Unpaired outcome impossible under --passthrough")`.
   2. `(None, None, None)` → **exit cleanly.** All three EOF simultaneously, expected end.
   3. `(Some(_), None, _)` → bail: `"Read 2 file is truncated — R1 has more reads than R2"`.
   4. `(None, Some(_), _)` → bail: `"Read 1 file is truncated — R2 has more reads than R1"`.
   5. `(Some(_), Some(_), None)` → bail: `"--passthrough file is truncated — R1/R2 have more reads than passthrough at record N"`.
   6. `(None, None, Some(_))` → bail: `"--passthrough file has more records than R1/R2 (extra records starting at record N)"`.
   7. `(Some(_), None, Some(_))` → bail: `"Read 2 file is truncated — R1 has more reads than R2"` (R2 truncation takes precedence over passthrough mismatch in the error message; either way both are off).
   8. `(None, Some(_), Some(_))` → bail: `"Read 1 file is truncated — R2 has more reads than R1"` (same precedence).

   *(Resolves Reviewer B §1.3: enumerate all 8 arms; A-Crit-2: passthrough-longer detection.)*
8. **Final stats.** `PairValidationStats::merge()` sums the three new fields across batches in the parallel path; the serial path increments them directly.
9. **Report.** The trimming report (`*_trimming_report.txt` for R1) gains a passthrough block:

       === Passthrough file ===
       Input:                          I1.fq.gz
       Output:                         I1_passthrough.fq.gz
       Records carried through:        9,872 (98.7%)
       Records dropped (pair filtered): 128 (1.3%)
       Sync checks performed:          10,000 (all matched)

   The JSON report (`*_trimming_report.json`) gains a `"passthrough"` object with the same fields. MultiQC compatibility: the new lines live in a dedicated block at the end of the existing report.
10. **FastQC.** When `--fastqc` is set, the post-trim FastQC sweep additionally runs over the passthrough output file. Note: cell-barcode reads (16–28 bp of uniformly-structured sequence) will produce a "bad-looking" FastQC report — per-base content bias, zero adapter content, near-zero N content. This is expected; the `--help` text for `--passthrough` warns about it. (Future: optional `--fastqc-skip-passthrough` flag if users complain.)

### Edge cases

- **Passthrough file shorter than R1/R2:** Caught at step 7.v (named file, names the row).
- **Passthrough file longer than R1/R2:** Caught at step 7.vi.
- **First record IDs do not match:** Caught at step 7.i.1 on the first iteration (loud, before any output is written in the serial path; before workers receive any batch in the parallel path).
- **Empty input files:** `FastqReader::sanity_check` at step 2 rejects empty input with a hard error. There is no zero-records path through the trimming pipeline. *(Corrects A-Other empty-input edge case.)*
- **Single record input:** One iteration, one partial batch in the parallel path (`BATCH_SIZE = 4096`).
- **All pairs dropped (`--discard_untrimmed` + zero adapters):** Output files are valid empty gzip members. Passthrough output is also a valid empty gzip member. Stats show 100% dropped.
- **Mid-file passthrough parse error:** Propagates via `?` through the reader thread. In the parallel path, the new error-propagation harness (Implementation Step 6a) signals workers and the main thread, which stops further output writes. Partial output files may remain on disk before the error reaches the caller — see assumption §13 for the v1 contract.
- **`@read/1` vs `@read/2` vs `@read/3` legacy headers:** Handled. `read_id_prefix()` strips trailing `/[123]` so the three streams sync up correctly. *(Resolves AB2 and Reviewer B §2.3.)*
- **`@read 1:N:0:CGATCG` modern Illumina headers:** Handled. `read_id_prefix()` splits on whitespace; all three streams share the prefix before the space.
- **Passthrough record with empty sequence:** No special handling — the passthrough writer emits a valid 4-line block with empty seq + empty qual. Valid FASTQ per spec.
- **Mixed gzip/plain inputs:** Each `FastqReader` autodetects per-file. **Output compression is uniform across all three outputs and is driven by R1's input compression** (see §Assumptions §11). If R1 is plain but the passthrough input is gzipped, the passthrough *output* will still be plain; the user should expect this.

---

## Signature

### New CLI field (`src/cli.rs`)

```rust
/// Pass a third FASTQ file through unchanged but keep it in lockstep with R1/R2.
/// Use case: 10X Multiome / scATAC libraries where the cell-barcode (I1/I2) read
/// must stay aligned to the trimmed R1/R2. Records dropped by length/quality/N
/// filters are also dropped from the passthrough output. The passthrough file is
/// never trimmed or adapter-scanned.
///
/// Requires --paired with exactly one R1/R2 pair. Incompatible with
/// --retain_unpaired, --clumpify, and all specialty modes (--clock, --implicon,
/// --hardtrim5/3, --demux).
///
/// Note: legacy headers like @read/1, @read/2, @read/3 sync correctly. If your
/// files use unusual header conventions and the sync check fires unexpectedly,
/// please file an issue with a sample header line. If --fastqc is enabled the
/// passthrough FastQC report will look poor (cell-barcode reads are intentionally
/// uniformly-structured) — that's expected, not a defect.
#[clap(long = "passthrough")]
pub passthrough: Option<PathBuf>,
```

### New output-naming helper (`src/io.rs`)

```rust
/// Generate the passthrough output filename for paired-end mode with --passthrough.
///
/// The naming follows the same `--basename` semantic as `paired_end_output_names`:
/// when `--basename foo` is set the output is `foo_passthrough.{fq,fq.gz}`; otherwise
/// the stem derives from the passthrough input filename via `strip_fastq_extensions`.
///
/// Note: `gzip` is the resolved global flag (driven by R1's input compression at
/// `main.rs:84`), not the passthrough input's extension.
pub fn passthrough_output_name(
    input_passthrough: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input_passthrough));
    let ext = if gzip { "_passthrough.fq.gz" } else { "_passthrough.fq" };
    let filename = format!("{}{}", stem, ext);

    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| {
            input_passthrough
                .parent()
                .unwrap_or(Path::new("."))
                .to_path_buf()
        });
    dir.join(&filename)
}
```

### Header-prefix sync-check helper (`src/trimmer.rs` or `src/fastq.rs`)

```rust
/// Extract the ID prefix used for paired/passthrough sync checks.
///
/// Steps:
///   1. Strip leading `@` if present (FastqRecord::id stores it).
///   2. Take everything before the first ASCII whitespace.
///   3. Strip a trailing `/1`, `/2`, or `/3` (legacy Illumina ID style).
///
/// Examples:
///   "@read1"                  → "read1"
///   "@read1 1:N:0:CGATCG"     → "read1"
///   "@read1/1"                → "read1"
///   "@read1/2"                → "read1"
///   "@read1/3"                → "read1"
///   "@read1/1 1:N:0:CGATCG"   → "read1"
///
/// *(Strips `/[123]` in v1 per dual-review consensus AB2; covers legacy SRA/ENA.)*
pub fn read_id_prefix(id: &str) -> &str {
    let stripped = id.strip_prefix('@').unwrap_or(id);
    let head = stripped.split_ascii_whitespace().next().unwrap_or("");
    head.rsplit_once('/')
        .filter(|(_, suf)| matches!(*suf, "1" | "2" | "3"))
        .map(|(p, _)| p)
        .unwrap_or(head)
}
```

**Unit tests required** (Step 4):
- `read_id_prefix("@read1") == "read1"`
- `read_id_prefix("@read1 1:N:0:CGATCG") == "read1"`
- `read_id_prefix("@read1/1") == "read1"`
- `read_id_prefix("@read1/2") == "read1"`
- `read_id_prefix("@read1/3") == "read1"`
- `read_id_prefix("@read1/1 1:N:0:CGATCG") == "read1"`
- `read_id_prefix("@read1/4") == "read1/4"` (only 1/2/3 strip)
- `read_id_prefix("") == ""`
- `read_id_prefix("@") == ""`

### Serial path signature (`src/trimmer.rs::run_paired_end`)

```rust
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    reader_passthrough: Option<&mut FastqReader>,  // NEW
    writer_r1: &mut FastqWriter,
    writer_r2: &mut FastqWriter,
    writer_passthrough: Option<&mut FastqWriter>,  // NEW
    mut unpaired_r1: Option<&mut FastqWriter>,
    mut unpaired_r2: Option<&mut FastqWriter>,
    config: &TrimConfig,
    unpaired: crate::filters::UnpairedLengths,
) -> Result<(TrimStats, TrimStats, PairValidationStats)>
```

`reader_passthrough` and `writer_passthrough` must be both `Some` or both `None`; mismatched `Some/None` is a programming error caught by `debug_assert!` at function entry.

### Parallel path signature (`src/parallel.rs::run_paired_end_parallel`)

```rust
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end_parallel(
    input_r1: &Path,
    input_r2: &Path,
    input_passthrough: Option<&Path>,    // NEW
    output_r1: &Path,
    output_r2: &Path,
    output_passthrough: Option<&Path>,   // NEW
    unpaired_r1_path: Option<&Path>,
    unpaired_r2_path: Option<&Path>,
    config: &TrimConfig,
    cores: usize,
    gzip: bool,
    unpaired: UnpairedLengths,
    clump_layout: Option<ClumpLayout>,
) -> Result<(TrimStats, TrimStats, PairValidationStats)>
```

### Worker payload + result (`src/parallel.rs`)

```rust
// Worker payload: third Vec is None when passthrough disabled.
type PairedWork = Option<(u64, Vec<FastqRecord>, Vec<FastqRecord>, Option<Vec<FastqRecord>>)>;

// Batch result: optional compressed third stream + a flag for whether passthrough
// was active so the main thread can decide whether to write `compressed_passthrough`.
// (Asymmetric to `compressed_unpaired_*` empty-Vec idiom — see Reviewer-A C3
// resolution: `Option<Vec<u8>>` makes the disabled case unambiguous.)
struct PairedBatchResult {
    seq: u64,
    compressed_r1: Vec<u8>,
    compressed_r2: Vec<u8>,
    compressed_passthrough: Option<Vec<u8>>,  // NEW; None ↔ disabled
    compressed_unpaired_r1: Vec<u8>,
    compressed_unpaired_r2: Vec<u8>,
    stats_r1: TrimStats,
    stats_r2: TrimStats,
    pair_stats: PairValidationStats,
}
```

### Worker function signature (`src/parallel.rs::process_pairs`)

```rust
#[allow(clippy::too_many_arguments)]
fn process_pairs<W: Write>(
    reads_r1: &mut [FastqRecord],
    reads_r2: &mut [FastqRecord],
    reads_passthrough: Option<&[FastqRecord]>,   // NEW — IMMUTABLE; written verbatim
    config: &TrimConfig,
    stats_r1: &mut TrimStats,
    stats_r2: &mut TrimStats,
    pair_stats: &mut PairValidationStats,
    writer_r1: &mut W,
    writer_r2: &mut W,
    writer_passthrough: Option<&mut W>,          // NEW — same gzip/plain fork as r1/r2
    mut writer_up_r1: Option<&mut W>,
    mut writer_up_r2: Option<&mut W>,
    unpaired: UnpairedLengths,
) -> Result<()>
```

**Inner iteration shape** (Step 7 detail):

```rust
let pt_iter = reads_passthrough.map(|s| s.iter());
let mut pt_iter = pt_iter; // dancing with `mut` for the `.next()`

for (r1, r2) in reads_r1.iter_mut().zip(reads_r2.iter_mut()) {
    let pt = pt_iter.as_mut().and_then(|it| it.next());
    match classify_paired(r1, r2, config, stats_r1, stats_r2, pair_stats, unpaired) {
        PairOutcome::Pass => {
            r1.write_to(writer_r1)?;
            r2.write_to(writer_r2)?;
            if let (Some(w), Some(rec)) = (writer_passthrough.as_deref_mut(), pt) {
                rec.write_to(w)?;
                pair_stats.passthrough_records_kept += 1;
            }
        }
        PairOutcome::Discarded => {
            if pt.is_some() {
                pair_stats.passthrough_records_dropped += 1;
            }
        }
        PairOutcome::Unpaired { .. } => debug_assert!(false, "impossible under --passthrough"),
    }
}
```

Same call shape is used in **both** the gzip and plain branches of `process_paired_batch` — the gzip branch wraps writers in scoped `GzEncoder`s, the plain branch passes `&mut Vec<u8>` directly. Implementer **must** mirror the existing two-branch fork in `parallel.rs:213–277`; the third writer is wrapped exactly the same way as `writer_r1` / `writer_r2`. *(Resolves AB6 worker-plumbing detail.)*

### Stats fields (`src/report.rs::PairValidationStats`)

```rust
#[derive(Debug, Default, PartialEq, Eq)]  // NEW: PartialEq, Eq for Validation §2 parity assert.
pub struct PairValidationStats {
    // Existing fields:
    pub pairs_analyzed: u64,
    pub pairs_removed: u64,
    pub pairs_removed_n: u64,
    pub pairs_removed_too_long: u64,
    pub r1_unpaired: u64,
    pub r2_unpaired: u64,

    // New fields (default 0; non-zero only when --passthrough is set):
    pub passthrough_records_kept: u64,
    pub passthrough_records_dropped: u64,
    pub passthrough_records_checked: u64,
}
```

`PairValidationStats::merge()` sums all three new fields. *(A-Crit-1: explicit `PartialEq, Eq` derive — Validation §2 won't compile without it.)*

---

## Implementation outline

Ordered for incremental testing. Each step has a checkpoint; commit (or at least `cargo test`) between steps.

### Step 1 — CLI + validation

1. Add `passthrough: Option<PathBuf>` to `Cli` in `src/cli.rs`.
2. Add the 9-item validation block (steps 1.i–1.ix above) to `Cli::validate()`, immediately after the existing `--clumpify` block.
3. **Step 1.ix uses the case-folded normalisation** consistent with `main.rs:179`'s output-collision pre-flight — declare a small helper, e.g. `fn norm_path(p: &Path) -> String { p.to_string_lossy().to_ascii_lowercase() }`, and use it in both places (or expose `main.rs`'s closure as a `pub(crate)` function).
4. Add unit tests in `src/cli.rs::tests` covering each rejection path. Add positive test `test_passthrough_paired_accepted` mirroring `test_validate_paired_two_files_accepted`. Add test for case-folded R1/R2 collision (`--passthrough r1.fq.gz` with positional `R1.fq.gz`).
5. **Checkpoint:** `cargo test cli::tests` and `cargo clippy --all-targets --release -- -D warnings` pass.

### Step 2 — Output naming + collision pre-flight

1. Add `passthrough_output_name()` to `src/io.rs`.
2. Add unit tests covering: bare input, `--output_dir`, `--basename`, `--dont_gzip`, all four input extensions.
3. In `src/main.rs` paired-end branch (line ~177), conditionally push the passthrough output path into the per-pair `candidates` vec when `cli.passthrough.is_some()`.
4. **Checkpoint:** `cargo test io::` passes.

### Step 3 — Stats fields + PartialEq

1. **Add `#[derive(PartialEq, Eq)]` to `PairValidationStats` at `src/report.rs:106`.** *(A-Crit-1: required for Validation §2.)*
2. Add the three new `u64` fields to `PairValidationStats`.
3. Extend `PairValidationStats::merge()` to sum the new fields.
4. Add a JSON-serialization entry for the new fields (whatever the existing pattern is).
5. **Checkpoint:** `cargo build` passes; existing tests still pass (new fields default to 0; `PartialEq` derive doesn't break anything because `TrimStats` already derives it).

### Step 4 — `read_id_prefix()` helper

Add the function (signature above) to `src/trimmer.rs` (or `src/fastq.rs` — implementer's choice). Add the 9 unit tests listed in §Signature.

**Checkpoint:** `cargo test` passes.

### Step 5 — Serial path (`run_paired_end`)

1. Extend the function signature with the two new `Option<&mut FastqReader>` / `Option<&mut FastqWriter>` parameters.
2. Add the `debug_assert!` at entry that both are `Some` or both `None`.
3. Inside the loop, replace the existing `match (rec1, rec2)` with the **8-arm three-way match** from Behavior §7. When `reader_passthrough.is_none()`, fall back to the 4-arm 2-tuple match (i.e., conditionally select the match shape).
4. On `Pass`, write the passthrough record + increment `passthrough_records_kept`. On `Discarded`, increment `passthrough_records_dropped`.
5. Increment `passthrough_records_checked` per iteration when passthrough is active.
6. **Checkpoint:** Run Validation §1 (end-to-end serial smoke). Stats accounting verified.

### Step 6 — Parallel path: reader thread

1. Extend `PairedWork` to include `Option<Vec<FastqRecord>>` for passthrough records.
2. Extend `read_pairs_round_robin` to accept an optional third `&mut FastqReader`. When `Some`, push one passthrough record per pair into a parallel `batch_passthrough` vec, and include it in the channel send.
3. **Sync check happens here** (reader thread, before dispatching to a worker). Use `read_id_prefix()` across all three records of each iteration.
4. Replace the existing `match (rec1, rec2)` with the **8-arm three-way match** from Behavior §7. The reader thread is the producer of every truncation/mismatch error.
5. **A-Crit-2: `(None, None, _)` arm must drain passthrough.** Even on a clean R1/R2 EOF, call `reader_passthrough.next_record()?` and bail if it returns `Some` (passthrough-longer case).
6. The passthrough records ride the channel verbatim — workers don't sync-check.
7. **Do not touch `read_pairs_clumpy`** — `--clumpify` + `--passthrough` is rejected at validation.

### Step 6a — Reader-thread error propagation hardening *(B-Crit-1: required)*

**Problem:** the existing `run_paired_end_parallel` pattern (`parallel.rs:185–193`) checks `reader_handle.join()` **after** the `result_rx.recv()` loop drains. A mid-stream reader error (truncation, sync mismatch, parse error) lets workers continue processing already-dispatched batches and the main thread writes them to disk before observing the reader error.

**Fix shape (v1, minimum-viable):** change the result channel payload from `PairedBatchResult` to `Result<PairedBatchResult, anyhow::Error>`. The reader thread, on any error, sends `Err(e)` directly to `result_tx` (capture a clone before moving `result_tx` into the worker scope) and then drops its worker channels so workers exit. The main thread loop checks `match result_rx.recv() { Ok(Ok(batch)) => ..., Ok(Err(e)) => return Err(e), Err(_) => break }` — propagating the reader error stops further writes immediately.

**Cleanup contract:** v1 does **not** delete partial output files on error. The user gets a clear error message and may need to clean up `*_val_1.fq.gz` / `*_val_2.fq.gz` / `*_passthrough.fq.gz` themselves. Document this in the `--help` text and the trimming report's failure path. *(Future v2: optional auto-cleanup with `--cleanup-on-error`.)*

**Test:** Validation §8 (new) — reader-thread error mid-stream returns `Err` and stops the main thread before all batches are written.

### Step 7 — Parallel path: worker + output

1. Extend `process_paired_batch` to accept `passthrough_records: Option<Vec<FastqRecord>>` and `retain_unpaired` is **already** in scope. Mirror the existing 2–4 GzEncoder pattern — when `passthrough_records.is_some()`, instantiate a 5th encoder (or 3rd in the no-retain-unpaired case).
2. Extend `process_pairs` per the §Signature subsection — `passthrough: Option<&[FastqRecord]>` (immutable), `writer_passthrough: Option<&mut W>`.
3. The inner zip-iteration shape is in §Signature.
4. Extend `PairedBatchResult` with `compressed_passthrough: Option<Vec<u8>>`.
5. In the main thread's flush loop, when `output_passthrough.is_some()` write `r.compressed_passthrough.as_ref().expect("active passthrough ⇒ Some")` to the passthrough output file. **Important:** the three writes (R1, R2, passthrough) must happen *within the same `while let Some(r) = pending.remove(&expected)` iteration*, before `expected += 1`. This guarantees the three output files stay in lockstep across batches.
6. **Worker stamps `pair_stats.passthrough_records_checked = batch.len()` once per batch** when `passthrough_records.is_some()`. *(Resolves Reviewer B §1.2.)*
7. **Checkpoint:** Validation §2 parity test passes (serial vs. parallel produce field-identical stats AND decoded-record-identical R1+R2+passthrough output).

### Step 8 — Wire main.rs

1. In `src/main.rs::run_paired` (or wherever the paired-end dispatcher lives), thread `cli.passthrough.as_deref()` and the computed passthrough output path through to `run_paired_end_parallel` / `run_paired_end`.
2. Open the passthrough reader/writer at the same call site as R1/R2.
3. The sanity-check on the passthrough file happens just before the call.

### Step 9 — Report

1. In `src/report.rs::generate_trimming_report`, conditionally emit the `=== Passthrough file ===` block at the bottom of the text report when `pair_stats.passthrough_records_checked > 0`. Include input filename, output filename, kept/dropped counts, percentages, sync-checks-performed count.
2. In the JSON report, conditionally add a `"passthrough": { ... }` object with the same fields.

### Step 10 — FastQC integration

In `src/main.rs` post-trim FastQC sweep, when `cli.fastqc` and `cli.passthrough.is_some()`, add the passthrough output path to the FastQC file list.

### Step 11 — Tests + docs

1. Add the 8 validation tests (§Validation §1–§8).
2. Add a new fixture `test_files/multiome_*` (see Validation §1 spec).
3. Update CLAUDE.md's test-fixtures section to mention the new fixture.
4. Verify CLI `--help` renders the new flag and the legacy-header / FastQC-output-bias notes.
5. **Do not** add this combination to the CI `validation` matrix — Perl 0.6.11 has no counterpart.

---

## Efficiency

- **Memory:** Each parallel batch grows from 2 vecs of records to 3. At `BATCH_SIZE = 4096` records × ~80 bytes/passthrough (16–28 bp barcodes are typical) ≈ +330 KB per batch. With `--cores 16` × 2 channel capacity = 32 in-flight batches × 330 KB ≈ +10 MB peak. Negligible. *(Note: passthrough records carried by 16 workers × `Vec<FastqRecord>` ≈ ~10 MB additional RSS — bounded and well below the user's `--memory` budget.)*
- **CPU:** Sync check is `split_ascii_whitespace().next()` + `rsplit_once('/')` filter ≈ ~30 ns × 3 calls per record. **At 1M reads: ~90 ms total** (corrected from v1's 30 ms — Reviewer A §3.4). Well under 0.1% of total runtime.
- **Channel pressure:** `sync_channel(2)` per worker; payload grows from 2 vecs to 3. No regression beyond the +330 KB per slot noted above.
- **Throughput:** Parallel-path saturation point at `--cores ≈ 8` is gzip-output-I/O-bound; adding the passthrough write is independent and does not contend with R1/R2 writes (separate files, separate encoders per worker). No expected throughput regression.
- **GzEncoder count:** 2–4 → 3–5 encoders per worker per batch. ~256 KB internal buffer each. +~4 MB steady state at `--cores 16`. Negligible.
- **Reader-thread bottleneck:** the reader does `next_record()` × 3 + `read_id_prefix()` × 3 per iteration instead of × 2. Reader was already the implicit serial bottleneck in `run_paired_end_parallel`; this adds ~5% reader-thread work. Worth measuring post-implementation on the benchmark fixture but not blocking.

---

## Integration

### What is read

- A new third FASTQ file (`cli.passthrough.as_deref()`), opened via `FastqReader::open_threaded` (parallel) or `FastqReader::open` (serial).

### What is written

- One new output file per invocation: `<stem>_passthrough.fq(.gz)`.
- Optionally one FastQC zip + HTML for the passthrough file when `--fastqc` is set.
- Modified trimming reports (text + JSON) gain a `=== Passthrough file ===` block conditionally.

### Order relative to other steps

- The passthrough record is read **at the same time** as R1/R2 (lockstep) in the reader thread.
- The passthrough write happens **after** `filter_paired_end` returns `Pass` (i.e. only for surviving pairs) — **in the same flush iteration** as R1/R2 (parallel) or the same loop iteration (serial).
- The passthrough writer is flushed at the same point as R1/R2 writers (at function exit / batch-end gzip `finish()`).

### Downstream impact

- **MultiQC:** purely additive content at the bottom of the report; the existing TrimGalore module won't break. Worth a manual MultiQC test post-implementation.
- **nf-core pipelines:** nf-core/atacseq or future nf-core/scatacseq could consume `--passthrough` output. Out of scope here but flagged for the relevant maintainers.
- **Reproducible builds:** No new wall-clock, hostname, or absolute-path uses.
- **CI `validation` job:** Not affected (no Perl 0.6.11 counterpart). `cargo test` gains 8 new tests.

---

## Assumptions

### Fixed rules

1. **Three-way input arity is strict.** `--passthrough` + `--paired` requires exactly 2 positionals (R1, R2) + 1 `--passthrough` value.
2. **Paired-end is required.** `--passthrough` without `--paired` is rejected.
3. **The passthrough file is FASTQ.** Enforced by `FastqReader::sanity_check` at startup.
4. **The passthrough file may be gzipped or plain**, independently of R1/R2.
5. **Output compression is uniform** across all three outputs and follows R1's input compression (resolved via `main.rs:84`'s `is_gzipped(&cli.input[0])`), not the passthrough input's extension. *(B-Other: documented.)*
6. **Strict pair semantics: no rescue.** `--retain_unpaired` rejected.
7. **Header-prefix sync check is mandatory and per-record.** No opt-out flag.
8. **The header-prefix definition:** strip `@`, take everything before first whitespace, strip trailing `/[123]`.
9. **The passthrough record is written record-identically** — id/seq/qual preserved verbatim. **NOT byte-identical** to the input file: line 3 (`+...`) is canonicalised to bare `+`, line endings canonicalise to `\n`. This matches the canonicalisation R1/R2 already get. *(B-Crit-2: softened.)*
10. **`--rename` does not affect the passthrough record's ID** (documented as R1/R2-only).
11. **Output gzip is driven by R1's input compression.** A user with plain R1 + gzipped passthrough input gets plain `_passthrough.fq` output. To force gzip output across mixed input compression, use `--compression <N>` (which implies gzip). *(B-Other.)*
12. **`PairValidationStats` derives `PartialEq, Eq`** for the load-bearing Validation §2 parity assert. *(A-Crit-1.)*
13. **Reader-error mid-stream may leave partial output files on disk.** v1 does not auto-clean. `--help` documents this; user can re-run after fixing the input. *(B-Crit-1, v1 contract.)*

### Configurable behavior

- `--passthrough <FILE>` — the flag itself.
- `--basename foo` — affects all three output stems uniformly (`foo_val_1`, `foo_val_2`, `foo_passthrough`).
- `--output_dir <DIR>` — all three outputs land in the same directory.
- `--dont_gzip` / `--compression <N>` — apply globally.
- `--fastqc` — runs FastQC on all three outputs when set.

### Things explicitly OUT of scope for v1

- Multi-pair input.
- Per-pair passthrough flags.
- Clumpy reordering with passthrough.
- Trimming the passthrough file.
- Custom passthrough output suffix (e.g. `--passthrough-suffix _index`).
- `--retain_unpaired` + passthrough rescue semantics. *(A clean v2 semantic exists — emit on any surviving mate — but not in v1 until a user asks.)*
- Auto-cleanup of partial output on reader error.

---

## Validation

Eight targeted tests covering the key failure points.

### §1 — End-to-end smoke (serial path)

**What:** Small 3-file fixture (R1, R2, I1 with 50 records each) runs to completion. The passthrough output has the expected record count and every passthrough record's ID matches the R1/R2 IDs of the same row.

**How:** Test fixture in `test_files/multiome_*.fq.gz`. Drive `run_paired_end` directly. Decode all three outputs. Assert (a) all three lengths equal, (b) row-by-row ID-prefix match.

**Expected:** All pass.

### §2 — Serial/parallel parity *(strengthened per AB1 + AB5)*

**What:** Serial vs. parallel produce field-identical `PairValidationStats` AND decoded-record-identical R1, R2, AND passthrough outputs on the same large input.

**How:**
- **Fixture: ≥10,000 paired records** (not 50). At `BATCH_SIZE = 4096` this forces ≥3 batches in the parallel path, exercising the `BTreeMap` ordered-flush logic.
- ≥30% of records must be filtered (e.g., `--length` cutoff at the right threshold) so the test actually exercises the drop path.
- Decode the passthrough output **and** R1+R2 outputs from both paths. Assert:
  - `serial_stats == parallel_stats` (depends on Step 3's `PartialEq` derive).
  - `decoded_serial_r1 == decoded_parallel_r1` (record-identical via `Vec<(String, String, String)>` of (id, seq, qual)).
  - Same for R2.
  - Same for passthrough.
  - **Row-by-row lockstep:** for every `i`, `read_id_prefix(passthrough[i].id) == read_id_prefix(r1_out[i].id) == read_id_prefix(r2_out[i].id)`.

**Expected:** All equalities hold. This is the load-bearing invariant.

### §3 — Sync check catches truncation (passthrough shorter)

**What:** Truncated passthrough (10 fewer records) → loud error naming the file.

**How:** Fixture `multiome_I1_truncated.fq.gz` with 40 records, R1/R2 with 50. Call `run_paired_end_parallel`. Assert `Err` containing "passthrough" and "truncated".

### §4 — Sync check catches ID mismatch (shuffled IDs)

**What:** Same record count, permuted IDs → hard error.

**How:** Fixture `multiome_I1_shuffled.fq.gz`. Call `run_paired_end_parallel`. Assert `Err` containing "Read ID mismatch" and the row number.

### §5 — CLI rejection regression tests

One test per `anyhow::bail!` added in `Cli::validate()` (steps 1.i–1.ix). Follow the `test_clumpify_rejects_*` pattern. **Includes the case-folded R1/R2 collision check (step 1.ix).**

### §6 — `--dont_gzip` plain-output coverage *(new per AB3)*

**What:** Run the parity test (§2 shape) with `gzip = false`. Verify the plain branch of `process_paired_batch` produces identical decoded output to the gzip branch.

**How:** Same fixture as §2, but configure `gzip = false`. Decode the plain output (no `MultiGzDecoder` wrapping). Assert record-identity against the gzip-decoded version.

### §7 — `--cores 1` dispatcher integration *(new per AB4)*

**What:** Exercise the `main.rs::run_paired` dispatcher path with `cores = 1` so the routing to `run_paired_end` (serial) goes through `main()`'s decision logic, not direct invocation.

**How:** CLI-integration test (in `tests/` top-level dir or via `std::process::Command` if a built binary is available). Pass `--paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz`. Assert all three outputs are produced and round-trip.

### §8 — Reader-thread error propagation *(new per B-Crit-1)*

**What:** A mid-stream reader error (sync mismatch at record ~5,000 of 10,000, `--cores 4`) returns `Err` to the caller, and the function does not hang. Partial outputs may exist on disk (v1 contract — see Assumption #13) — assert only on the error return, not on file state.

**How:** Use the §4 shuffled fixture but with the mismatch placed mid-stream (e.g., at row 5,000 in a 10,000-record file, so it falls mid-batch in the parallel path). Call `run_paired_end_parallel` with `--cores 4`. Assert (a) `Err` is returned, (b) the error message names the file and the row, (c) the call returns within a reasonable wall-time (no deadlock).

---

## Questions or ambiguities

### Resolved in v2

1. ~~Header-prefix definition — strip `/[123]`?~~ → **Yes, in v1.** Strips trailing `/1` `/2` `/3` per dual-review consensus AB2.
2. ~~Trimming report block visual style.~~ → Match existing `=== ... ===` blocks (implementer's choice).
3. ~~Passthrough output stem — passthrough input's basename or pair's basename?~~ → **Passthrough input's basename** (consistent with `--basename foo` overriding all three uniformly).
4. ~~Stats placement (`Option<PassthroughStats>` vs extending `PairValidationStats`)?~~ → **Extend `PairValidationStats`.** *(Contradiction resolved per main-loop judgement.)*
5. ~~Severity of `--passthrough == R1/R2` collision check?~~ → **Case-folded check, required.** *(Contradiction resolved per main-loop judgement.)*

### Still open (non-blocking)

- **FastQC for cell-barcode reads produces a "bad-looking" report** (uniform sequence content → per-base bias warning). v1 accepts this; `--help` documents it. If users complain, add `--fastqc-skip-passthrough` opt-out in a follow-up.
- **Auto-cleanup of partial outputs on reader error.** Deferred. v1 leaves partial files on disk + clear error message.

### Critical — NONE

All previously-flagged critical items resolved.

---

## Self-Review

### Efficiency

- Memory growth quantified per-batch and at `--cores 16` scale. Bounded and well under typical `--memory` budgets.
- Sync-check cost corrected to ~90 ms / 1M reads (was 30 ms in v1; A-Other catch).
- Reader-thread overhead flagged for post-implementation measurement, but small enough not to block.

### Logic

- 8-arm three-way EOF match enumerated explicitly (B §1.3 catch).
- Reader-thread error propagation has a concrete v1 fix shape (B-Crit-1 catch).
- `passthrough_records_checked` accounting site is decided (B §1.2 catch): worker stamps from batch length.
- Worker-thread plumbing for `process_pairs` spelled out: `&[FastqRecord]` immutable for passthrough vs `&mut [FastqRecord]` for R1/R2; gzip/plain fork mirrored (AB6 catch).

### Edge cases

- Empty input correctly handled by `sanity_check` rejection at startup (A-Other catch — v1 was wrong).
- Passthrough-longer-than-R1/R2 has an explicit arm in the 8-arm match (A-Crit-2 catch).
- Legacy `/1` `/2` `/3` headers handled via `read_id_prefix()` stripping (AB2 catch).
- Mid-stream reader error: v1 documents partial-output contract; tests for fail-fast behaviour.

### Integration

- No regression on existing `--paired` runs.
- `PartialEq, Eq` derive on `PairValidationStats` is purely additive — no existing test depends on the absence of these traits (A-Crit-1).
- Validation §2 fixture sized to genuinely exercise multi-batch ordering, not just stat-counting (AB1).

### Remaining risks

- **MultiQC parsing of the new report block** — purely additive content, low risk, but worth a 2-minute manual check post-implementation.
- **Reader-thread error propagation fix (Step 6a)** is the highest-leverage change in this revision. The implementer should be careful to test the fail-fast path (Validation §8) and confirm the existing parallel R1/R2-only tests still pass after the channel-payload change from `PairedBatchResult` to `Result<PairedBatchResult>`.
- **FastQC output for cell-barcode reads will look poor.** Documented; if users complain post-release, opt-out flag is a 5-line follow-up.
- **`--retain_unpaired` rejection** is explicit about the v2 path being open if a user asks (no longer pretending there's "no clean interpretation"; B §1.8 catch).
