# Design plan — `--clump_only` Phase 2 (uBAM in/out)

**Issue:** [FelixKrueger/TrimGalore#353](https://github.com/FelixKrueger/TrimGalore/issues/353) (contributor: @wkgardner)
**Prior work:** `plans/clump-only/PLAN.md` (v1 FASTQ, 438 lines) + PR #355 (in-flight)
**Branch:** `feature/clump-only-ubam` (off `feature/clump-only`)
**Target PR:** stacks on #355; opens after #355 merges

---

## Goal

Extend the v1 FASTQ-only `--clump_only` specialty mode with unaligned BAM (uBAM) input and output, matching the trim path's existing `--output-format ubam` support. After v2, `--clump_only` accepts uBAM input (via content-based `format::detect_input_format`), accepts `--output-format ubam` for uBAM output, and honors `--preserve-tags` for aux-tag round-trip through the FASTQ-record intermediate. Byte-identity of the semantically-meaningful fields (id + sequence + quality + preserved aux tags) is preserved; the `@PG` line in the BAM header adds provenance and is the sole reason cross-run byte-identity of the output file requires `@PG`-ignoring comparison (same treatment as the trim uBAM path).

## Context

### Placement

`--clump_only` gets a **uBAM-output dispatch branch** alongside the existing FASTQ-output branch. On the input side, `format::open_sync_reader` already returns a `Box<dyn RecordSource>` that dispatches by content type, so uBAM input works transparently once the current v1 `reject_ubam` guard is removed. Sibling to how `hardtrim5` grew a `hardtrim5_to_bam` variant.

### Verified references (verified against source at start of this plan)

- `src/specialty.rs::hardtrim5_to_bam` (line 106) — the load-bearing template for "specialty mode + uBAM output" functions. Signature: `(input: &Path, keep: usize, output_dir: Option<&Path>, rename: bool, preserve_tags: &[String], command_line: &str) -> Result<()>`.
- `src/bam.rs::BamWriter::create(path, source_header: Option<&Header>, preserve_tags: &[String], command_line: &str) -> Result<Self>` (line 528).
- `src/bam.rs::BamWriter::write_record(&mut self, record: &FastqRecord, paired_side: Option<u8>) -> Result<()>` (line 555) — `paired_side` is `None` for SE, `Some(1)` / `Some(2)` for PE. Mate-adjacent PE output happens by alternating `Some(1)` and `Some(2)` calls.
- `src/bam.rs::BamWriter::finish(self) -> Result<()>` (line 601) — consumes self; produces a valid BAM even for zero records.
- `src/bam.rs::BamReader::open` / `open_threaded` / `open_paired_interleaved` / `open_threaded_with_tags` (lines 96–146).
- `src/bam.rs::peek_header(path) -> Result<Header>` (line 813) — used to grab the input `@HD`/`@PG` chain before opening a reader.
- `src/bam.rs::build_output_header(source: Option<&Header>, command_line: &str) -> Result<Header>` (line 619) — appends our `@PG ID:trim_galore VN:… CL:…` to the source chain.
- `src/format.rs::open_sync_reader(path, preserve_tags) -> Result<Box<dyn RecordSource>>` — polymorphic reader factory.
- `src/format.rs::detect_input_format(path) -> Result<InputFormat>` — returns `FastqPlain | FastqGz | UnalignedBam` from content.
- `src/fastq.rs::RecordSource` — trait implemented by both `FastqReader` and `BamReader`.
- `src/main.rs::run_ubam_output` (line 1481) + `run_ubam_output_single` (1588) + `run_ubam_output_paired_two_files` (1709) — the trim-path uBAM dispatch, informative but NOT reused directly (trim-specific).

### v1 primitives to reuse unchanged

- `clump::canonical_minimizer`, `bin_for`, `sort_single_by_key`, `sort_paired_by_key` — sort primitives already stable with content-tiebreaker cascade (confirmed by both plan-reviewers).
- `clump::resolve_layout` and `clump::estimated_record_bytes` — memory sizing and per-record byte accounting (see Open Question 3 for the BAM-aux-tag-byte question).
- `clump_only::SingleBin` and `PairedBin` — bin buffer structs from v1; reused verbatim.
- `clump_only::write_clump_only_report` — text report writer; needs a small extension for BAM-output byte accounting (see Behavior §Report shape).

## Behavior

### Contract (extended from v1)

For every input record R:

- `R.id` (name portion + preserved aux tags) round-trips byte-identically through the FASTQ-record intermediate. **A/Z/i/f scalar tags** round-trip losslessly; **B (array)** and **H (hex)** tags are rejected at BAM-read time per the existing `bam.rs` constraint. This is codebase-wide behaviour inherited from the trim uBAM path, not `--clump_only`-specific.
- `R.seq`, `R.qual` byte-identical.
- Input aux tags NOT listed in `--preserve-tags` are dropped (same as trim uBAM path).

For the **output file as a whole**:

- Header `@HD` + input `@PG` chain preserved verbatim.
- New `@PG ID:trim_galore VN:<version> CL:<cmdline>` line appended (same as trim uBAM path).
- Cross-run byte-identity of the whole file **does NOT hold** because `CL:` in `@PG` varies with the exact invocation string. CI byte-identity check uses `assert_ubam_eq`-shape comparison (ignores `@PG`; matches trim uBAM CI pattern).
- Record contents **do** cross-run-byte-identically (identical input + identical `--cores`/`--memory`/`--compression` → identical body bytes; only the header's `@PG.CL` varies).

### SE steps (`--clump_only --output-format ubam` on SE input)

1. Parse CLI → early-dispatch to `clump_only_single_to_bam` (bypass trim pipeline).
2. **Do not** call `reject_ubam` (v1's guard is gone on the BAM path); do accept FASTQ input too (uBAM output from FASTQ input is a legitimate v2 combo — same as trim's `--output-format ubam` with FASTQ input).
3. `peek_header(input)` if input is uBAM; `None` if FASTQ.
4. Open reader via `format::open_sync_reader(input, &cli.preserve_tags)`.
5. Create output writer: `BamWriter::create(&output_path, source_header.as_ref(), &cli.preserve_tags, &command_line)`.
6. Stream records into per-bin `SingleBin` buffers, computing `canonical_minimizer` per record — same loop as v1's `clump_only_single`.
7. On bin overflow: sort via `sort_single_by_key`; **for each sorted record**, `writer.write_record(&record, None)?`. Records land in BAM natively via `paired_side: None`.
8. At EOF: flush remaining bins in bin-index order.
9. `writer.finish()?` — produces a valid BAM (including for zero records — verified by BamWriter's design; see Open Question 4).
10. Write reorder-only report (see Report shape).

### PE steps (`--clump_only --paired --output-format ubam`)

Two input shapes to support, plus one to reject:

**Shape A: two files, both FASTQ.** `--clump_only --paired --output-format ubam R1.fq(.gz) R2.fq(.gz)`. Mirrors `--paired` FASTQ-input convention. Extends to N=4, 6, … files (multi-pair, matching v1's `run_specialty_paired` iteration over `.chunks(2)`).

**Shape B: single interleaved uBAM.** `--clump_only --paired --output-format ubam interleaved.bam`. Matches how the trim uBAM path handles `--paired` + N=1 uBAM (`Cli::validate`'s N=1-with-uBAM carve-out at `cli.rs:466`). Uses `BamReader::open_paired_interleaved_with_tags` to de-interleave into two internal R1/R2 streams.

**Rejected: two BAM files under `--paired`.** `--clump_only --paired --output-format ubam R1.bam R2.bam` bails with a clear "uBAM paired mode expects a single interleaved file" message, mirroring the trim path's rejection at `main.rs:1502-1514`. Rationale: BAM's pair semantics are file-internal (mate flag bits + query-name grouping); two-file paired BAM is semantically ambiguous.

**Rejected: N=1 with non-BAM input under `--paired`.** `--clump_only --paired --output-format ubam one.fq.gz` (or any non-BAM N=1) bails with a clear message. Rationale: `BamReader::open_paired_interleaved` requires BAM input; opening a FASTQ file with it would produce a cryptic BAM-magic error deeper in noodles.

**Rejected: mixed input formats in Shape A.** `--paired R1.fq.gz R2.bam` (or vice versa) bails with a clear "input pair must be same format" message. Rationale: cross-format pairing is nonsense; catching it explicitly beats letting one reader open successfully then the other fail with a cryptic error.

In accepted shapes (A and B), the output is a **single interleaved BAM per pair** at `<stem>_clumped.bam` — mate-adjacent, matching samtools/Picard/fgbio convention and the existing `paired_bam_output_name` convention (`*_val.bam` on the trim path). Multi-pair Shape A produces one output BAM per input pair.

Steps (per pair):

1. **Pre-loop input-shape guard** — validate all inputs at once before opening any reader (see rejection matrix above). Emit clear error and bail if any guard fires.
2. **Pre-loop collision pre-flight** — for each pair (Shape A) or single input (Shape B), compute the prospective output-BAM path, hash case-folded via `naming::norm_path`, insert into a `HashMap`. On duplicate → bail with the standard "output path collision (case-insensitive, for APFS/NTFS safety)" message. Mirrors v1's SE pre-flight pattern.
3. `peek_header(input)` on the pair's R1 (Shape A) or the interleaved file (Shape B). None for FASTQ-input.
4. Open readers: `format::open_sync_reader(r1, preserve_tags)?` + `format::open_sync_reader(r2, preserve_tags)?` (Shape A) OR `BamReader::open_paired_interleaved_with_tags(input, preserve_tags)?` (Shape B).
5. Create ONE output writer: `BamWriter::create(&pair_output_path, source_header.as_ref(), preserve_tags, command_line)`.
6. Stream (R1, R2) pairs into `PairedBin` buffers — same loop as v1's `clump_only_paired`.
7. On bin overflow: `sort_paired_by_key`; for each sorted pair index i, write `writer.write_record(&r1[i], Some(1))?` then `writer.write_record(&r2[i], Some(2))?`. Mate-adjacent per samtools convention.
8. At EOF: flush remaining bins in bin-index order.
9. `writer.finish()?`.
10. Write reorder-only report — one per pair (multi-pair Shape A produces one report per pair, matching v1 SE's one-report-per-input convention).
11. If `--fastqc`: `fastqc::run(&pair_output_path, cli.fastqc_args, output_dir, cores)?` — runs on the reordered BAM (fastqc-rust reads BAM natively, verified against crate source).

### Rejection matrix (updated from v1)

**Now accepted** under `--clump_only`:
- `--output-format ubam` — new dispatch branch.
- uBAM input (SE or PE-interleaved N=1) — no longer bailed by `clump_only::reject_ubam`.
- `--preserve-tags TAG,…` — composes with uBAM in/out; ignored (with existing warning) on FASTQ input.
- `--fastqc` + `--output-format ubam` — runs FastQC on the reordered BAM (fastqc-rust reads BAM natively; verified against crate source at `~/.cargo/registry/src/…/fastqc-rust-1.0.1/src/runner.rs`). This diverges from the trim uBAM path today, which silently skips FastQC on BAM output — the trim skip is an oversight worth fixing in a separate follow-up issue, not v2 scope.

**Newly rejected — hoisted into `Cli::validate`'s shared §3.4a `if OutputFormat::UBam { … }` block** so it applies to ALL uBAM output paths (not just `--clump_only`):
- `--dont_gzip` + `--output-format ubam` — BAM is always BGZF-compressed; `--dont_gzip` has no meaning. **Note:** this closes an existing gap — the trim uBAM path currently does NOT reject this combination (only `--clumpify` does, at `cli.rs:650-652`). Hoisting into §3.4a is a real correctness improvement for the trim path as a side benefit.

**Newly rejected — under `--clump_only` specifically** (see PE-steps above for user-facing messages):
- Two-BAM input under `--paired` (Shape A with both files BAM).
- N=1 non-BAM input under `--paired` (Shape B with FASTQ).
- Mixed input formats in Shape A (FASTQ + uBAM in the same pair).

**Unchanged from v1** — still rejected under `--clump_only`:
- All trim / filter / clip / adapter / poly-* / nextseq / rename / discard_untrimmed / consider_already_trimmed / other-specialty-modes / passthrough / retain_unpaired flags.
- `-q` / `--stringency` / `-e` still silently ignored (clap defaults on non-Option types).

### Report shape (extended from v1)

`ClumpOnlyStats` struct extended with:

```rust
pub struct ClumpOnlyStats {
    // v1 fields:
    pub total_records: u64,
    pub input_bytes: u64,
    pub output_bytes: u64,
    pub n_bins: usize,
    pub peak_bin_occupancy: usize,
    pub compression_level: u32,
    pub input_compressed: bool,
    pub output_compressed: bool,
    // v2 additions:
    pub input_format: InputFormatLabel,      // "FASTQ (plain)" | "FASTQ (gzip)" | "uBAM"
    pub output_format: OutputFormatLabel,    // "FASTQ (plain)" | "FASTQ (gzip level N)" | "uBAM (BGZF)"
}
```

Report content:

```
Trim Galore version: 2.x.x
Mode: --clump_only (lossless reorder)
Input:  <path> (<format>, <bytes>)
Output: <path> (<format>, <bytes>)                   [PE-BAM: one line, one file]
Records: <N>                                          [PE-BAM: N pairs]
Bins: <n_bins> (peak occupancy <M>)
Compression ratio: <x.yy>x   [omitted if input/output compression states aren't comparable]
Preserved tags: TAG1,TAG2,…                          [omitted if empty]
```

- PE-BAM output emits ONE report (one interleaved output file → one output_bytes value → one report). Diverges from v1's per-mate PE report.
- Compression ratio is emitted only when the comparison is meaningful (both sides gzip or both sides BGZF; input-uBAM → output-uBAM is a meaningful ratio; FASTQ→BAM is not — see Open Question 2 for exact rule).

## Signature

```rust
// New functions in src/clump_only.rs

#[allow(clippy::too_many_arguments)]
pub fn clump_only_single_to_bam(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    cores: usize,
    memory_budget_bytes: u64,
    preserve_tags: &[String],
    command_line: &str,
    fastqc: bool,
    fastqc_args: Option<&str>,
    no_report_file: bool,
) -> Result<ClumpOnlyStats>;

/// Runs on one pair of inputs. Caller (main.rs) drives multi-pair iteration
/// over `cli.input.chunks(2)` (or a single-file Shape B special-case for N=1).
#[allow(clippy::too_many_arguments)]
pub fn clump_only_paired_to_bam_one_pair(
    inputs: &[PathBuf],           // 1 element = interleaved uBAM (Shape B); 2 elements = R1 + R2 (Shape A)
    output_dir: Option<&Path>,
    basename: Option<&str>,
    cores: usize,
    memory_budget_bytes: u64,
    preserve_tags: &[String],
    command_line: &str,
    fastqc: bool,
    fastqc_args: Option<&str>,
    no_report_file: bool,
) -> Result<ClumpOnlyStats>;

// New filename helpers in src/io.rs (mirroring paired_bam_output_name convention)
pub fn clumped_bam_output_name(input: &Path, output_dir: Option<&Path>, basename: Option<&str>) -> PathBuf;
pub fn clumped_paired_bam_output_name(input_r1: &Path, _input_r2: Option<&Path>, output_dir: Option<&Path>, basename: Option<&str>) -> PathBuf;
```

Notes on signatures:

- No `gzip_output` parameter — BAM is always BGZF; the concept doesn't apply.
- No `compression` parameter — BAM's compression level is fixed at BGZF's default (verify against `BamWriter::create` internals during implementation; if it accepts a level, thread it through).
- **`fastqc: bool` + `fastqc_args: Option<&str>` restored.** Per §Resolved decision 1, `--fastqc` runs on the reordered BAM (fastqc-rust supports BAM natively).
- Added `preserve_tags: &[String]` and `command_line: &str` — mirrors `hardtrim5_to_bam`.
- Kept `no_report_file: bool` from v1's post-review remediation (fix #2).
- **Signature name change**: `clump_only_paired_to_bam` → `clump_only_paired_to_bam_one_pair` to signal that it processes ONE pair; caller drives multi-pair iteration. Restores v1's N=4, 6, … multi-pair support that a "flat inputs slice" signature would have regressed.

## Implementation outline

Ordered, each step independently reviewable.

### 1. CLI validation edits (four precise sites)

**Site 1 — `src/clump_only.rs::reject_ubam` (currently line 233):** delete this function AND remove its two call sites (`clump_only_single`, `clump_only_paired` — currently lines 264 and 391-392). v1's per-function uBAM guard is gone; format dispatch now happens transparently via `format::open_sync_reader`.

**Site 2 — `src/cli.rs`'s `if self.clump_only { … }` block (currently around line 685):** delete the block that reads `if matches!(self.output_format, OutputFormat::UBam) { bail!("--clump_only + --output-format ubam is not yet supported in v1 …") }`.

**Site 3 — `src/cli.rs::validate_paired_input` guard for N=1 (around line 701-707 in the `--clump_only` block):** delete the v1 rejection that bails on `--clump_only --paired` with N=1. Shape B (N=1 interleaved BAM) is now legal. **BUT** replace it with a narrower guard that only bails when N=1 AND the single input is NOT a BAM (see PE-steps §Rejected: N=1 with non-BAM input). This narrower guard is a format-detect at CLI-validate time, matching how the trim uBAM path handles the same distinction at `cli.rs:466`.

**Site 4 — `src/cli.rs`'s shared §3.4a `if OutputFormat::UBam { … }` block (currently around `cli.rs:516-580`):** ADD a new rejection: `if self.dont_gzip { anyhow::bail!("--dont_gzip is not compatible with --output-format ubam (BAM is always BGZF-compressed)"); }`. Placement in §3.4a (not inside the `if self.clump_only` block) means this rejection also applies to the existing trim uBAM path — closing an existing gap (the trim path currently accepts `--dont_gzip --output-format ubam` silently, which is a real inconsistency).

**Keep intact:** all other v1 rejections in the `--clump_only` block (adapter flags, length/max_length/max_n, clip, RRBS, poly-*, nextseq, rename, discard_untrimmed, consider_already_trimmed, other specialty modes, passthrough, retain_unpaired).

**Add:** Shape A two-BAM rejection and mixed-format Shape A rejection. These are format-dependent so live at the input-detection layer (early in the dispatch in `main.rs`), not in `Cli::validate` — see Step 2.

### 2. Main dispatch update (src/main.rs)

Extend the existing `if cli.clump_only { … }` block to dispatch on `cli.output_format`. The BAM branch has explicit multi-pair iteration + Shape B special-case + collision pre-flight + format-guard checks.

```rust
if cli.clump_only {
    let memory_bytes = clump::parse_memory_size(&cli.memory)?;
    let basename = cli.basename.as_deref();
    match cli.output_format {
        OutputFormat::Fastq => { /* existing v1 dispatch, unchanged */ }
        OutputFormat::UBam => {
            if cli.paired {
                // Shape B: N=1 interleaved BAM (Cli::validate already guaranteed
                // N=1 → BAM via the narrowed guard at site 3 above).
                if cli.input.len() == 1 {
                    // Collision pre-flight (single output path, but keeps the
                    // pattern consistent).
                    let out = naming::clumped_paired_bam_output_name(
                        &cli.input[0], None, output_dir, basename);
                    // (SE-style HashMap pre-flight — for consistency; degenerate
                    // when N=1 but harmless).
                    preflight_collision(&[out.clone()])?;
                    clump_only::clump_only_paired_to_bam_one_pair(
                        &cli.input, output_dir, basename, cli.cores,
                        memory_bytes, &cli.preserve_tags, &command_line,
                        cli.fastqc, cli.fastqc_args.as_deref(),
                        cli.no_report_file,
                    )?;
                } else {
                    // Shape A: two-file paired, N=2/4/6/… (multi-pair).
                    // Format-guard: reject two-BAM Shape A + mixed-format Shape A.
                    for chunk in cli.input.chunks(2) {
                        let fmt_r1 = detect_input_format(&chunk[0])?;
                        let fmt_r2 = detect_input_format(&chunk[1])?;
                        if matches!(fmt_r1, InputFormat::UnalignedBam)
                            || matches!(fmt_r2, InputFormat::UnalignedBam)
                        {
                            let is_bam = |f| matches!(f, InputFormat::UnalignedBam);
                            if is_bam(fmt_r1) && is_bam(fmt_r2) {
                                bail!(
                                    "--paired with two BAM files is not supported. \
                                     uBAM paired mode expects a single interleaved file: \
                                     `trim_galore --clump_only --paired interleaved.bam`. \
                                     Got two BAM files; one of them is {}.",
                                    chunk[0].display(),
                                );
                            } else {
                                bail!(
                                    "--paired requires both input files to be the same format. \
                                     Got mixed: {} and {}.",
                                    chunk[0].display(), chunk[1].display(),
                                );
                            }
                        }
                    }
                    // Multi-pair collision pre-flight (case-folded per issue #216).
                    let mut planned: Vec<PathBuf> = Vec::new();
                    for chunk in cli.input.chunks(2) {
                        planned.push(naming::clumped_paired_bam_output_name(
                            &chunk[0], Some(&chunk[1]), output_dir, basename));
                    }
                    preflight_collision(&planned)?;
                    // Per-pair iteration.
                    for chunk in cli.input.chunks(2) {
                        clump_only::clump_only_paired_to_bam_one_pair(
                            chunk, output_dir, basename, cli.cores,
                            memory_bytes, &cli.preserve_tags, &command_line,
                            cli.fastqc, cli.fastqc_args.as_deref(),
                            cli.no_report_file,
                        )?;
                    }
                }
            } else {
                // SE BAM (possibly multi-input). Reuse v1's SE collision pre-flight
                // shape but with the BAM naming helper.
                let mut planned: Vec<PathBuf> = Vec::new();
                for input in &cli.input {
                    planned.push(naming::clumped_bam_output_name(input, output_dir, basename));
                }
                preflight_collision(&planned)?;
                for input in &cli.input {
                    clump_only::clump_only_single_to_bam(
                        input, output_dir, basename, cli.cores,
                        memory_bytes, &cli.preserve_tags, &command_line,
                        cli.fastqc, cli.fastqc_args.as_deref(),
                        cli.no_report_file,
                    )?;
                }
            }
        }
    }
    return Ok(());
}
```

**`preflight_collision(paths: &[PathBuf]) -> Result<()>`** — a small helper (either inline in main.rs or hoisted to `io::preflight_collision`) that hashes case-folded paths via `naming::norm_path` and bails on duplicate with the standard "output path collision (case-insensitive, APFS/NTFS safety)" message. Same logic as v1's SE pre-flight; extract-and-share to avoid duplication.

### 3. Core BAM implementations (src/clump_only.rs)

Add two new functions mirroring the v1 pair but with BAM I/O:

**`clump_only_single_to_bam`**:

1. `peek_header(input)?` if uBAM; `None` if FASTQ.
2. Resolve `ClumpLayout` from `memory_budget_bytes` (same as v1).
3. Output path = `naming::clumped_bam_output_name(input, output_dir, basename)`.
4. Open reader = `format::open_sync_reader(input, preserve_tags)?` (dispatches by content type).
5. Create writer = `BamWriter::create(&output_path, source_header.as_ref(), preserve_tags, command_line)?`.
6. Same bin-dispatch loop as `clump_only_single` (streaming reader → per-bin `SingleBin` → flush on overflow).
7. On bin flush: sort via `sort_single_by_key`; **for each record**, call `writer.write_record(&rec, None)?`. Do NOT reuse `write_records_member` (that's the gzip-member helper for FASTQ path).
8. At EOF: flush remaining bins in bin-index order.
9. `writer.finish()?`.
10. Populate `ClumpOnlyStats` with input_bytes / output_bytes / total_records / n_bins / peak / input_format / output_format labels.
11. Write report if `!no_report_file`.

**`clump_only_paired_to_bam`**:

Same shape but with the input-shape dispatch:

```rust
let (mut reader_r1, mut reader_r2, source_header): (
    Box<dyn RecordSource>, Box<dyn RecordSource>, Option<Header>,
) = match inputs.len() {
    1 => {
        // Shape B: single interleaved BAM.
        let header = peek_header(&inputs[0])?;
        let (r1, r2) = BamReader::open_paired_interleaved_with_tags(&inputs[0], preserve_tags)?;
        (Box::new(r1), Box::new(r2), Some(header))
    }
    2 => {
        // Shape A: two files (uBAM or FASTQ each).
        let header = match detect_input_format(&inputs[0])? {
            InputFormat::UnalignedBam => Some(peek_header(&inputs[0])?),
            _ => None,
        };
        let r1 = format::open_sync_reader(&inputs[0], preserve_tags)?;
        let r2 = format::open_sync_reader(&inputs[1], preserve_tags)?;
        (r1, r2, header)
    }
    _ => bail!("--clump_only --paired --output-format ubam requires 1 (interleaved) or 2 (paired) input files, got {}", inputs.len()),
};
```

Then the standard pair-loop from v1 `clump_only_paired`, replacing FASTQ writers with a single `BamWriter` and writing mate-adjacent with `Some(1)` / `Some(2)`.

### 4. Filename helpers (src/io.rs)

Add:

- `clumped_bam_output_name(input, output_dir, basename) -> PathBuf` → `<stem>_clumped.bam`.
- `clumped_paired_bam_output_name(input_r1, _input_r2, output_dir, basename) -> PathBuf` → `<stem>_clumped.bam` (ONE file for PE, matching `paired_bam_output_name`'s `_val.bam` convention).

Both honor `basename` when set (`<basename>_clumped.bam`).

### 5. Report writer extension (src/clump_only.rs)

Add `input_format: InputFormatLabel` and `output_format: OutputFormatLabel` enum fields to `ClumpOnlyStats`. Update `write_clump_only_report` to render them in the text output. Adjust the "Compression ratio" line's emit-condition: only emit when input and output are both gzip-family (gzip↔gzip, BGZF↔BGZF, or gzip↔BGZF as an approximation — decide at implementation time; erring toward "omit unless meaningful" is safer).

### 6. `estimated_record_bytes` audit (src/clump.rs)

**Verified during plan-review** — both plan-reviewers (A and B) confirmed that `estimated_record_bytes` at `src/clump.rs:250-252` counts the full `rec.id.len()` in its byte total, which correctly accounts for the tab-separated aux-tag suffix when `--preserve-tags` is set (aux tags are appended to the id via `bam_record_to_fastq` in `bam.rs`). **No change needed.** Step retained for the audit trail.

### 7. Tests (unit + integration)

Unit tests in `src/clump_only.rs`:
- `test_clump_only_single_to_bam_permutation` — build synthetic FASTQ input, run to BAM, decode via `noodles::bam::Reader`, assert record multiset preserved.
- `test_clump_only_paired_to_bam_lockstep` — synthetic FASTQ pair → interleaved BAM output, decode, assert R1[i]/R2[i] pair lockstep at every output position (mate-adjacent).
- `test_clump_only_bam_deterministic_records` — two runs, decode both, assert record-body byte-identity (ignoring `@PG` in header).
- `test_clump_only_ubam_in_ubam_out_tag_roundtrip` — synthetic uBAM input with `--preserve-tags CB,UB`, verify output uBAM has the same tags on the same records.

Integration tests in `tests/integration_clump_only_ubam.rs` (new file):
- `se_ubam_out_from_fastq_in` — FASTQ input, `--clump_only --output-format ubam`, decode output via `samtools view` in CI (or via noodles in the test), assert record multiset.
- `se_ubam_in_ubam_out_pg_chain` — uBAM input with pre-existing `@PG`, verify output preserves it and appends TrimGalore's `@PG`.
- `pe_ubam_out_interleaved_from_fastq_pair` — two FASTQ inputs, `--clump_only --paired --output-format ubam`, assert output is ONE interleaved BAM with mate-adjacent records.
- `pe_ubam_in_interleaved_output` — single interleaved uBAM input via N=1-under-`--paired`, output ONE interleaved BAM, tag round-trip verified.
- `rejects_dont_gzip_with_ubam_output` — `--clump_only --output-format ubam --dont_gzip` → exit non-zero with clear message (rejection at cli.rs §3.4a).
- `rejects_two_bam_paired` — `--clump_only --paired --output-format ubam R1.bam R2.bam` → exit non-zero; error mentions "single interleaved file" (mirrors trim uBAM rejection at main.rs:1502-1514).
- `rejects_non_bam_n1_paired` — `--clump_only --paired --output-format ubam one.fq.gz` → exit non-zero; clear error identifying the format mismatch.
- `rejects_mixed_format_paired` — `--clump_only --paired --output-format ubam R1.fq.gz R2.bam` → exit non-zero; clear error identifying the format mismatch.
- `multi_pair_pe_bam_produces_one_output_per_pair` — `--clump_only --paired --output-format ubam A_R1.fq A_R2.fq B_R1.fq B_R2.fq` produces `A_R1_clumped.bam` and `B_R1_clumped.bam` (two files, one per pair). Regression guard for the v1 multi-pair capability.
- `fastqc_produces_report_on_ubam_out` — `--clump_only --output-format ubam --fastqc` produces a `*_fastqc.html` + `*_fastqc.zip` alongside the BAM output. Confirms the §Resolved decision 1 direct-run behavior.
- `pe_bam_collision_preflight_case_folded` — two input pairs whose case-folded output paths collide (e.g. `A_R1.fq`/`a_R1.fq` on APFS) → exit non-zero at pre-flight (before any reader opens).

### 8. CI validation steps (.github/workflows/ci.yml)

New steps in the `validation` job:
- **SE uBAM byte-identity**: run `--clump_only --output-format ubam` on a uBAM fixture; decode both input and output via `samtools view` (already available on the CI runner); sort by read name; diff. Assert record multiset preserved.
- **PE uBAM interleaving**: run `--clump_only --paired --output-format ubam` on a paired uBAM fixture; assert mate-adjacency (samtools flagstat + custom check).
- **`@PG` chain preservation**: assert output BAM's header contains both the source `@PG` and our TrimGalore `@PG`.
- **Cross-run record determinism** (body-only, ignoring `@PG`): two runs, `samtools view <out>` piped to `sort` on both, md5 compare.

### 9. Documentation (docs/src/content/docs/modes/clump-only.md)

Update the existing v1 page:
- New "uBAM in / uBAM out" section describing the composition with `--output-format ubam` and `--preserve-tags`.
- Document the `@PG` line addition and the "byte-identity ignores `@PG`" carve-out.
- Document the `--dont_gzip` incompatibility with `--output-format ubam` (now enforced at cli.rs §3.4a for all uBAM output paths).
- Note that `--fastqc` runs on BAM output natively (fastqc-rust supports BAM).
- Document the PE input-shape acceptance / rejection matrix: Shape A (two-file FASTQ paired, multi-pair supported), Shape B (single interleaved uBAM), rejected combinations (two-BAM Shape A, non-BAM N=1, mixed-format Shape A).
- Remove the "v1 is FASTQ in / FASTQ out only" caveat.
- Add PE-uBAM output naming note: **ONE interleaved BAM per pair**, matching samtools/Picard/fgbio. Multi-pair produces one output BAM per input pair.

### 10. CHANGELOG.md

Draft entry under `### Unreleased` (or a new versioned header if v2 gets its own release):

> - **`--clump_only` now supports uBAM in/out** — extending the v1 lossless reorder mode with unaligned BAM support via `--output-format ubam` (input auto-detected by content). Aux-tag round-trip via `--preserve-tags TAG1,TAG2,…` (A/Z/i/f scalars; B/H tags rejected at BAM-read time, same as trim path). PE-uBAM output is ONE interleaved BAM per pair (mate-adjacent, matching samtools/Picard/fgbio). Multi-pair input under `--paired` is preserved from v1 (`N=4, 6, …` FASTQ inputs → one output per pair). Output BAM's `@PG` line records the invocation; input `@PG` chain preserved. `--dont_gzip + --output-format ubam` now rejected (previously silently accepted on trim uBAM path — this closes an existing gap). `--fastqc` runs on BAM output (fastqc-rust reads BAM natively). Byte-identity preserved for record contents (id + seq + qual + preserved tags); cross-run byte-identity of the whole file uses `@PG`-ignoring comparison (identical to the trim uBAM path's CI treatment).

## Efficiency

- Same O(N log N / n_bins) per-bin sort as v1. BAM I/O throughput is dominated by BGZF encode/decode, which uses `noodles`'s exact-pinned `=0.88.0` — inherits its perf.
- **No shared-buffer temporary** for BAM output — `BamWriter::write_record` streams directly to the output BGZF encoder (unlike v1's per-bin gzip-member buffered-Vec pattern), so the peak-memory doubling issue Reviewer B flagged in v1 (findings #7) does not apply to the BAM path. Bonus.
- Aux-tag round-trip adds proportional cost to the `id` string per record (typically 40–150 additional bytes per record for common tags). `estimated_record_bytes` should account for this via the `id.len()` term; verify at implementation time (see step 6).

## Integration

### Reads

- SE uBAM: `BamReader::open` via `format::open_sync_reader`.
- PE uBAM (Shape A, two files): two independent `format::open_sync_reader` calls.
- PE uBAM (Shape B, one interleaved file): `BamReader::open_paired_interleaved_with_tags` — bounded-buffer per-side de-interleave, MAX_SLACK=1024, existing.

### Writes

- SE uBAM: ONE `<stem>_clumped.bam` file per input.
- PE uBAM: ONE interleaved `<stem>_clumped.bam` file per pair (both shapes A and B produce the same output shape; multi-pair Shape A produces multiple output files). Shape B's report filename derives from the single interleaved input's stem (`<stem>_clumped.bam` where `<stem>` = `strip_bam_extension(input)`).
- Report: `<stem>_clumping_report.txt` — SE and multi-pair PE emit one per input/pair; Shape B (N=1 interleaved) emits one report for the pair (same stem as the BAM output).

### Order relative to other steps

- Runs before any trim-pipeline code (early dispatch from `main.rs::main`, returns immediately).
- FastQC hook: on FASTQ output, same as v1. On BAM output, `fastqc::run` invoked on the reordered BAM (fastqc-rust supports BAM natively).
- `run_ubam_output` (trim path): completely separate code path; no interaction.

### Downstream impact

- Existing `--clumpify` mode: zero effect.
- Existing trim uBAM path: zero effect.
- Perl-parity CI validation: unaffected (this is a new code path with its own invariants).
- nf-core / MultiQC pipelines: continue to see `*_clumping_report.txt` (not scanned by trim-report globs).

## Assumptions

Fixed rules (inherited from v1 + trim uBAM path):

- `noodles = "=0.88.0"` unchanged (exact-pinned for byte-identity of uBAM I/O per CLAUDE.md).
- BAM output is always BGZF-compressed. No plain-BAM option.
- Aux-tag preservation is A/Z/i/f scalars only. B (array) and H (hex) rejected at BAM-read time by existing `bam.rs` code.
- Aligned BAM input is rejected per-record by existing `BamReader::next_record` (via `is_unmapped()` check in `bam.rs`); mixed-aligned uBAM inputs cannot silently produce wrong output.
- `@PG` line format matches the trim uBAM path exactly (`ID:trim_galore VN:<CARGO_PKG_VERSION> CL:<command_line>`).
- `--paired` + N=1 with uBAM input is legal (existing carve-out in `Cli::validate::validate_paired_input`); non-BAM N=1 under `--paired` is rejected explicitly by the new format-detect guard.
- `BamWriter::finish` produces a valid header-only BAM for zero records (verified during plan-review at `bam.rs:1641-1647`).

Configurable (via CLI):

- `--preserve-tags TAG1,TAG2,…` — which aux tags to round-trip.
- `--cores` — accepted at any value ≥ 1 (v1 is single-threaded internally; **carrying over v1's deviation from `--clumpify`'s `>= 2` rule**. `--cores` currently has no effect on the clump-only code path but is retained for interface parity; parallelism is a v1.1/v2.1 follow-up).
- `--memory` — bin buffer sizing (same as v1).
- `--basename BASE` — overrides the output stem: `<BASE>_clumped.bam`.
- `--no_report_file` — suppresses the reorder-only report.
- `--fastqc` / `--fastqc_args` — runs fastqc-rust on the output BAM (fastqc-rust supports BAM natively).

## Validation

Failure-point-targeted checks. Each should be a test or CI step.

1. **Record byte-identity through BAM round-trip** — decode input and output BAMs to `(id, seq, qual, tags)` tuples; sort by id; assert exact per-record match. Kills: id-string mangling, seq/qual mishandling, tag drop or reorder in the FastqRecord intermediate. Location: unit test + CI step.
2. **Pair lockstep in interleaved output** — for PE, decode output BAM; assert record[i] and record[i+1] share the same name (mate pair) at every even i. Kills: unpaired R1 or R2 landing next to unrelated mate; interleaving error. Location: unit test + CI step.
3. **`@PG` chain preservation** — decode output BAM header; assert both the input's `@PG` records AND our new TrimGalore `@PG` record are present, in order. Kills: header replacement (destroying provenance) or missing appended `@PG`. Location: integration test.
4. **Cross-run record determinism (ignoring `@PG`)** — two runs; `samtools view` piped to `sort` on both; md5 diff on the body. Kills: unstable sort, non-deterministic BAM writer ordering, header-line reordering in the body region. Location: CI step.
5. **Rejection matrix on `--dont_gzip` + `--output-format ubam`** — assert exit non-zero with clear message. Kills: silent acceptance of an incompatible combination. Location: unit test in `cli.rs::tests` + integration test.
6. **`--preserve-tags` round-trip** — synthetic uBAM with `RG:Z:sample`, `BC:Z:AAAA`, `NM:i:0`, `AS:f:1.5`; run with `--preserve-tags RG,BC,NM,AS`; assert all four tags present on every output record with same values. Kills: tag ordering swap, type-code mistranslation, value-string corruption. Location: unit test.
7. **PE input shape dispatch** — one-file interleaved input AND two-file input both produce the same interleaved output shape. Kills: dispatch confusion between the two shapes. Location: two integration tests.
8. **`--fastqc` produces a report on BAM output** — assert `*_fastqc.html` and `*_fastqc.zip` exist alongside the BAM. Confirms fastqc-rust reads the BAM directly (per §Resolved decision 1). Kills: regression to the old "warn-and-skip" behavior we ruled out. Location: integration test.
9. **Empty-input BAM handling** — empty uBAM input (header only, zero records); output is a valid BAM decodable by `samtools view` with zero records. Kills: BamWriter panicking on zero records, or producing a truncated/invalid BAM. Location: unit test.

## Resolved decisions

All open questions resolved by Felix before dual plan-review. Locked below.

**Critical ambiguities: none.** All resolutions internally consistent with v1 and with the trim uBAM path.

1. **`--fastqc` on uBAM output → run FastQC on the BAM directly.** fastqc-rust natively supports `.bam`/`.ubam`/`.sam` (verified against `~/.cargo/registry/src/.../fastqc-rust-1.0.1/src/runner.rs` — accepts `--format bam|sam|bam_mapped|sam_mapped|fastq` and recognizes `.bam` and `.ubam` extensions). No warn-and-skip, no temp-FASTQ intermediate. Note: the existing trim uBAM path in `main.rs::run_ubam_output_*` does NOT call `fastqc::run` today; that's an oversight to fix in a separate follow-up issue, not scope for this v2 PR.

2. **Compression-ratio line: emit when both sides are compressed** (any of gzip / BGZF), omit when either side is plain. BGZF and gzip both use deflate under the hood, so cross-family ratio comparisons are meaningful for measuring the clumping win. Report-writer code should conditionally emit based on `input_compressed && output_compressed`.

3. **`estimated_record_bytes` behavior**: **verified during plan-review** — both plan-reviewers confirmed the current implementation at `src/clump.rs:250-252` counts the full `rec.id.len()` in the byte total, which correctly accounts for tab-separated aux-tag suffixes. No change needed at implementation.

4. **Empty-input BAM handling**: **verified during plan-review** — Reviewer B cited `src/bam.rs:1641-1647` as the code that produces a valid header-only BAM for zero-record inputs. `BamWriter::finish` handles the empty case correctly. Unit test #9 kept as a regression guard.

5. **`@PG` chain unbounded growth under repeated invocations**: accepted. Matches the trim uBAM path's behavior exactly. Provenance is valuable and deduping by `ID` would destroy it; deduping by exact-CL would create a divergence from the trim path that's worse than the "verbose chain" it fixes.

## Self-Review

**Efficiency**: No unnecessary passes. Reuses v1's `SingleBin` / `PairedBin` / sort primitives / `resolve_layout`. `BamWriter::write_record` streams directly to BGZF (no per-bin Vec buffering — sidesteps v1's memory-doubling concern from Reviewer B). Aux-tag round-trip adds a proportional per-record cost to the FastqRecord's id string, handled by the existing `format::open_sync_reader` machinery.

**Logic**: Steps in Behavior are ordered (parse → format-guard PE inputs → collision pre-flight → per-pair loop → peek header → open reader(s) → open writer → stream/bin → sort → write mate-adjacent → finish → report → optional FastQC). Multi-pair Shape A iterates over `.chunks(2)` (restoring v1's N=4, 6, … support that a "flat inputs slice" signature would have regressed). Shape B (N=1 interleaved BAM) is a top-level special case in the PE branch.

**Edge cases**:
- Empty input → validation test #9 (verified during plan-review at `bam.rs:1641-1647`).
- Single-record input → covered by generic permutation test with N=1.
- All-N sequences → same as v1 (canonical_minimizer folds N to A; no BAM-specific concern).
- Malformed BAM → BamReader's existing error handling propagates.
- BGZF decode error mid-stream → propagates through `BamReader::next_record` as `Result::Err`.
- Disk full mid-write → `BamWriter::write_record` propagates the io error; `writer.finish()` may not run (via `?`), which could leave a truncated BAM. Same behavior as trim uBAM path; not `--clump_only`-specific.
- **Two-BAM Shape A / non-BAM Shape B / mixed-format Shape A** — all rejected at pre-loop format-guard with clear messages.

**Integration**:
- No shared state with trim pipeline.
- No changes to v1 FASTQ path (all v1 unit + integration tests should stay green).
- No shared writer state with `--clumpify` (each mode owns its output).
- Perl-parity CI unaffected.
- **Side benefit**: hoisting the `--dont_gzip + --output-format ubam` rejection into `cli.rs` §3.4a closes a pre-existing gap where the trim uBAM path silently accepted the combination.

**Adjustments made during post-plan-review revision (dual reviewer feedback incorporated)**:
- **Fixed §Rejection matrix vs §Resolved decisions plan drift** — §Rejection matrix, §Signature, §Implementation, §Validation, §Integration all now agree that `--fastqc` runs on BAM output (matching §Resolved decision 1). The pre-decision "warn-and-skip" position has been purged from all sections.
- **Explicit multi-pair PE support restored** — signature renamed `clump_only_paired_to_bam_one_pair`; main.rs dispatch iterates `.chunks(2)`. Regression from v1's N=4, 6, … support that a flat-slice signature would have introduced.
- **Explicit format-guard rejections** for the two PE shape holes both reviewers surfaced: two-BAM Shape A (matching trim uBAM path's rejection at `main.rs:1502-1514`), non-BAM Shape B (format-detect at CLI-validate time, mirroring `cli.rs:466`), and mixed-format Shape A (new).
- **Explicit collision pre-flight** for PE-BAM and multi-SE-BAM paths (Reviewer A finding C4). Shared `preflight_collision(paths: &[PathBuf])` helper.
- **`--dont_gzip + --output-format ubam` rejection hoisted** from the `--clump_only`-specific block to the shared §3.4a block in `Cli::validate`. Closes a pre-existing gap on the trim uBAM path (Reviewer B finding #5).
- **Precise CLI edit sites** — Reviewer B finding #4: the reject-uBAM guard lives in `clump_only.rs::reject_ubam` (line 233), NOT in `cli.rs`. §Implementation outline step 1 now enumerates four precise edit sites with line references.
- **Q3 and Q4 both marked as verified during plan-review** — no code changes at implementation for either.

**Remaining risks**:
- **Shape B report-filename convention**: single-file interleaved input uses input's stripped-BAM-extension stem for `<stem>_clumped_report.txt`. If the input's basename overlaps with a Shape A pair name in the same directory, the report file could clash with a pair-run report. Not a common workflow (users don't typically mix Shape A + Shape B invocations in the same directory), but noted.
- `noodles = "=0.88.0"` pin is unchanged; if a future security patch requires a bump, that's a broader byte-identity re-verification concern that this v2 inherits from v1 but doesn't originate.
- **Trim uBAM path FastQC-skip is a pre-existing bug** highlighted by §Resolved decision 1's rationale. Fix is out of scope for v2 (a separate small PR) but worth filing as a follow-up.
