# Design plan — `--clump_only` (lossless reorder-only mode)

**Issue:** [FelixKrueger/TrimGalore#353](https://github.com/FelixKrueger/TrimGalore/issues/353) (contributor: @wkgardner, 2026-07-24)
**Branch:** `feature/clump-only` (off `dev`)
**Target release:** unreleased on `dev`; ships in the next minor/patch after merge

---

## Goal

Add a run-and-exit specialty mode `--clump_only` that reorders FASTQ records by canonical 16-mer minimizer (reusing the existing `--clumpify` machinery) **without any trimming, filtering, adapter detection, or record modification**. Every input record appears in the output byte-identically; only its position in the file changes. Purpose: enable Trim Galore to be used as a standalone lossless-recompression stage in archival pipelines that need `--clumpify`'s gzip-window benefit but require record-level fidelity.

## Context

### Placement

`--clump_only` is a **run-and-exit specialty mode**, siblings-wise with `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`, `--demux`. Dispatch pattern already established in `src/main.rs` (early-`if` on the flag, call the handler, `return`). The trim pipeline is never entered.

### Existing primitives to reuse (all in `src/clump.rs`)

- `pub fn canonical_minimizer(seq: &[u8]) -> MinimizerKey` — line 83
- `pub fn bin_for(key: MinimizerKey, n_bins: usize) -> usize` — line 123
- `pub fn sort_single_by_key(records: &mut Vec<FastqRecord>, keys: &mut Vec<MinimizerKey>)` — line 271
- `pub fn sort_paired_by_key(...)` — line 289
- `pub fn resolve_layout(memory_budget_bytes: u64, cores: usize) -> Result<ClumpLayout>` — line 196
- `pub fn clumpify_min_memory_bytes(cores: usize) -> u64` — line 157
- `pub fn parse_memory_size(s: &str) -> Result<u64>` — line 226

### Existing clumpy workers to *not* reuse directly

`src/parallel.rs::read_single_clumpy` (line 1129) and `read_pairs_clumpy` (line 776) bake in the full trim → clip → filter → gzip pipeline. We do **not** thread a "skip trim" boolean through these; we own a parallel path in `src/specialty.rs` (or a new sibling module) that reuses only the sort primitives.

### Existing CLI validation to reference

`src/cli.rs::validate()` at line 512 already rejects `--clumpify` in combination with: `--output-format ubam` (517), `--clock`, `--implicon`, `--hardtrim5`, `--hardtrim3`, `--demux`, `--dont_gzip`, `--passthrough` (614-641, 674), and enforces `--cores >= 2` (617). `--clump_only`'s rejection matrix mirrors this shape.

### Existing output naming (`src/io.rs`)

- SE trim: `<stem>_trimmed.fq(.gz)`
- PE trim: `<stem>_val_1.fq(.gz)` / `<stem>_val_2.fq(.gz)`
- SE uBAM: `<stem>_trimmed.bam`

`--clump_only` proposes: SE `<stem>_clumped.fq(.gz)`, PE `<stem>_clumped_1.fq(.gz)` / `<stem>_clumped_2.fq(.gz)`.

### Existing report (`src/report.rs`)

`TrimStats` struct with per-adapter / per-quality / per-filter counters, emitted as `<stem>_trimming_report.txt` and `<stem>_trimming_report.json`. **Not reused here** — see report filename in Behavior below.

## Behavior

### Contract (load-bearing invariant)

**For every input record R in the input FASTQ, exactly one record R' appears in the output such that**:
- `R'.header == R.header` (byte-exact)
- `R'.sequence == R.sequence` (byte-exact)
- `R'.quality == R.quality` (byte-exact)

**And the multiset of output records equals the multiset of input records** (permutation, no adds, no drops, no dedup).

**File-level byte-identity across runs**: on the same input, two invocations with the same `--cores` / `--memory` / `--compression` produce byte-identical output files. This requires deterministic ordering (see below).

### Contract-scope note: plus-line and line-endings are normalized

The byte-identity claim above is scoped to the three semantically-meaningful fields (header, sequence, quality). Two normalizations happen codebase-wide (not `--clump_only`-specific), inherited from `src/fastq.rs`:

- **Plus-line**: `FastqRecord` (`src/fastq.rs:32-39`) has no plus-line field; the reader discards line 3 and the writer hard-codes `+\n`. Inputs with `+<header-repeat>` on line 3 emerge with bare `+`. This is existing TrimGalore behavior for all read paths.
- **Line-endings**: reader's `trim_end_matches(['\n', '\r'])` accepts CRLF and LF input; writer emits LF-only output. CRLF inputs are normalized to LF.

Both normalizations are documented in the `--clump_only` help text and docs page so the user-facing byte-identity claim is honest about what it covers.

### SE steps

1. Parse CLI; if `--clump_only` set → dispatch to `run_clump_only_single` (bypass all trim-pipeline setup).
2. Skip adapter auto-detection (never called on this path).
3. Resolve `ClumpLayout` from `--memory`, `--cores` (same call as `--clumpify` today).
4. **Read** input FASTQ streamingly, computing `canonical_minimizer` per record; buffer records into per-bin `Vec<FastqRecord>` up to memory budget.
5. **Sort** each bin by `MinimizerKey`, with a **secondary tiebreaker** that guarantees determinism (see Determinism below).
6. **Write** bins in bin-index order to `<stem>_clumped.fq(.gz)` using the resolved compression level. Each bin becomes one gzip member (as today), so output is a valid RFC-1952 concatenation.
7. Emit a **reorder report** to `<stem>_clumping_report.txt` (and `.json`) — see Report shape below.
8. If `--fastqc` was set, run bundled FastQC on the reordered output (same integration point as the trim path, minus the trimming-report dep).

### PE steps

Same as SE with per-pair records. `sort_paired_by_key` already exists in `clump.rs` and preserves R1/R2 lockstep. Pair lockstep is guaranteed by construction: sort is by the **pair's** minimizer (typically R1's; verify by reading the current `sort_paired_by_key` header) and R1[i] and R2[i] move together.

### Determinism (within-bin ordering)

The current `sort_single_by_key` / `sort_paired_by_key` — must be verified to use **stable sort** (Rust's `sort_by_key`, not `sort_unstable_by_key`) so that records sharing a minimizer preserve input order. If they currently use unstable sort, we either:
- (a) switch them to stable sort (may regress `--clumpify` perf a hair — acceptable, verify via `cargo bench` if a bench exists),
- (b) add a secondary tiebreaker (e.g. `(minimizer_key, input_index)`) applied only in the `--clump_only` path.

**Recommendation:** stable sort everywhere (option a). Simpler, guarantees same-input-same-output for `--clumpify` too, and any perf delta is measurable.

### Rejection matrix (all → error in `Cli::validate`)

`--clump_only` is incompatible with:

**Trimming / clipping / filtering:**
- `-a`, `-a2`, `--adapter`, `--adapter2`
- `--illumina`, `--nextera`, `--small_rna`, `--bgi`, `--stranded_illumina`
- `--length`, `--max_length`, `--max_n`
- `--trim-n`, `--clip_r1`, `--clip_r2`, `--three_prime_clip_r1`, `--three_prime_clip_r2`
- `--rrbs`, `--non_directional`
- `--polyA`, `--polyG`, `--no_poly_g`
- `--nextseq` / `--2colour`
- `--rename` (mutates read IDs — would break byte-identity)
- `--discard_untrimmed` (drops reads that weren't adapter-trimmed — meaningless here)
- `--consider_already_trimmed`

**Other specialty modes** (all mutually exclusive):
- `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`, `--demux`, `--clumpify` (redundant with `--clump_only`)

**Output-shape combinations:**
- `--output-format ubam` — deferred to v2
- `--passthrough` — nothing to trim, passthrough concept doesn't apply
- `--retain_unpaired` — no unpaired discards happen here

**Input-format:**
- **uBAM input** — deferred to v2. Reject at `format::detect_input_format` dispatch time (not just CLI validation) with a clear "uBAM input is not yet supported under `--clump_only`" message. FASTQ-only in v1.

**Silently ignored** (see "Default-value flags" note below):
- `-q` / `--quality` — has clap default (20); can't distinguish user-set from default at validate time. Docs and `--help` explicitly note this is ignored under `--clump_only`.
- `--stringency` — has clap default; same story.
- `-e` / `--error` — has clap default; same story.

**Compatible** (composes normally):
- `--paired`
- `--fastqc` (opt-in)
- `--compression` (1–9)
- `--memory` (any parseable size)
- `--cores` (>= 1; accepted for interface parity but v1 is single-threaded internally — see PROGRESS.md §"Documented deviations")
- `--dont_gzip` — **allowed** (unlike `--clumpify`, which rejects it). Clump-only's archival use case may legitimately want lossless plain-FASTQ output as an intermediate stage for a downstream compressor. Note that this diverges from `--clumpify`'s rejection ("clumping plain text is pointless"): here the reorder itself is the primary output, not a gzip-window optimisation.
- `--basename` — must respect user-provided output stem; naming becomes `<basename>_clumped.fq(.gz)` (SE) or `<basename>_clumped_{1,2}.fq(.gz)` (PE). Verify at implementation time that filename helpers honor `cli.basename`.

**Default-value flags note**: `-q`, `--stringency`, `-e` are clap-typed with default values (not `Option<T>`), so `Cli::validate()` cannot distinguish "user passed" from "default used" without threading `ArgMatches::value_source()`. Rather than refactoring the derive-parser pattern, we silent-accept and document ignore-behavior in `--help` + docs. This matches how existing specialty modes (`--hardtrim5`, `--clock`, `--implicon`) treat trim-related flags today — the flags exist in the CLI surface but are not exercised on the specialty code path.

### Report shape

**Filename**: `<stem>_clumping_report.txt` (SE, per-input). PE report layout mirrors `--clumpify`'s current behavior — determine at implementation time by reading how the existing clumpify path names its per-input reports.

**No JSON**: the trimming report's `.json` exists specifically for nf-core/rnaseq's MultiQC parser. This report is deliberately **outside** the `*_trimming_report.*` scan glob (see below), so machine-parseable JSON is ceremonial — the text report is short enough to grep. Fewer serialization formats to keep in sync as the report shape evolves.

**Content (text)**:
```
Trim Galore version: 2.3.x
Mode: --clump_only (lossless reorder)
Input:  <path> (<gzip|plain>, <bytes>, <records> reads)
Output: <path> (<gzip level N | plain>, <bytes>, <records> reads)
Bins:   <n_bins> (peak occupancy <M> records)
Compression ratio: <input_bytes>/<output_bytes> = <x.yx>   [omitted if --dont_gzip or input was plain]
```

Omit the "Compression ratio" line when input and output compression states differ (gzip→plain or plain→gzip): the ratio becomes misleading or nonsensical. Emit it only when both input and output are gzip.

**Why not `_trimming_report.txt`**: nf-core/rnaseq and other pipelines scan `*_trimming_report.txt` by convention for MultiQC. Producing an empty-of-trim-stats file under that name would either confuse those parsers or silently produce misleading stats. New filename → new file-glob → no false positive for existing pipelines.

## Signature

```rust
// src/specialty.rs (or new src/clump_only.rs — see Implementation outline)

/// Reorder-only mode: read FASTQ, sort by canonical 16-mer minimizer, write
/// records byte-identically to a new file. No trimming, no filtering.
///
/// Byte-identity invariant: every input record appears in the output with
/// header, sequence, and quality string byte-exact. Only file-level order
/// changes. Deterministic for identical (input, cores, memory, compression).
pub fn clump_only_single(
    input: &Path,
    output_dir: Option<&Path>,
    gzip_output: bool,
    cores: usize,
    memory_budget_bytes: u64,
    compression: u32,
) -> Result<ClumpOnlyStats>;

pub fn clump_only_paired(
    r1: &Path,
    r2: &Path,
    output_dir: Option<&Path>,
    gzip_output: bool,
    cores: usize,
    memory_budget_bytes: u64,
    compression: u32,
) -> Result<ClumpOnlyStats>;

/// Filename helpers, mirroring specialty.rs conventions.
/// `gzip = false` (via `--dont_gzip`) produces `*_clumped.fq` instead of `*_clumped.fq.gz`.
pub fn clump_only_output_name(input: &Path, output_dir: Option<&Path>, gzip: bool) -> PathBuf;
pub fn clump_only_output_name_pe(input: &Path, mate: &str /* "1" | "2" */, output_dir: Option<&Path>, gzip: bool) -> PathBuf;

/// Reorder-report statistics — deliberately narrower than TrimStats.
#[derive(Debug, Default)]
pub struct ClumpOnlyStats {
    pub input_bytes: u64,
    pub output_bytes: u64,
    pub total_records: usize,
    pub n_bins: usize,
    pub peak_bin_occupancy: usize,
    pub compression_level: u32,
    pub input_compressed: bool,
    // ... whatever else the report needs
}
```

## Implementation outline

Ordered, each step independently reviewable.

### 1. CLI (src/cli.rs)

- Add `pub clump_only: bool` field to `Cli` with `#[clap(long = "clump_only")]`.
- Add doc comment mirroring the tone of other specialty modes (samples: `hardtrim5`, `clock`).
- Extend `Cli::validate()`:
  - Reject with clear message when `--clump_only` combined with any of the "Trimming / clipping / filtering" flags (loop over the group; one bail per group category with the specific flag named for good UX).
  - Reject when combined with other specialty modes and with `--output-format ubam` / `--passthrough` / `--retain_unpaired`.
  - Enforce `--cores >= 1` (v1 single-threaded; see PROGRESS.md §"Documented deviations" for the rationale for relaxing from `--clumpify`'s `>= 2` rule).
  - Enforce that `--memory` parses and is `>=` `clumpify_min_memory_bytes(cores)`.

### 2. Main dispatch (src/main.rs)

Insert an early-dispatch block **before** the existing `--hardtrim5` block (or wherever fits the ordering, alphabetical within specialty modes is fine):

```rust
if cli.clump_only {
    // No adapter auto-detection; no trim pipeline setup.
    let memory_bytes = clump::parse_memory_size(&cli.memory)?;
    if cli.paired {
        run_specialty_paired(...);  // wraps clump_only_paired
    } else {
        for input in &cli.inputs {
            let stats = specialty::clump_only_single(...)?;
            report::write_clump_only_report(...)?;
        }
    }
    return Ok(());
}
```

### 3. Core implementation (src/specialty.rs, OR new src/clump_only.rs)

**Decision:** create **`src/clump_only.rs`** as a sibling module — the reorder logic is substantial enough (streaming reader, bin dispatcher, sort, writer) that colocating with hardtrim/clock would clutter `specialty.rs`. `mod.rs`-level exposure via `pub mod clump_only;` in `src/lib.rs`.

Contents:

- `clump_only_single(...)`: streaming loop over `FastqReader`, per-record: compute `canonical_minimizer`, dispatch into per-bin `Vec<FastqRecord>` (up to `layout.records_per_bin_target`), flush bins that hit budget through the same pool pattern `read_single_clumpy` uses today — but with the trim call replaced by a **no-op write** (just `writer.write_record(&r)?;`).
- Same for `clump_only_paired(...)`, using `read_pairs_clumpy` as the structural reference.
- Filename helpers: `clump_only_output_name*` produce `*_clumped.fq(.gz)` / `*_clumped_{1,2}.fq(.gz)`. Honor `cli.basename` when set (produces `<basename>_clumped.fq(.gz)` / `<basename>_clumped_{1,2}.fq(.gz)`).
- **Do not use** `trimmer::run_single_end` / `run_paired_end` — those are the trim pipeline. This module is standalone.
- **uBAM input rejection**: after `format::detect_input_format` classifies the input, if it's uBAM, bail with a clear "uBAM input is not yet supported under --clump_only; see #353 for the v2 follow-up" error. Do this before opening any reader, so no partial work happens.
- **`--fastqc` wiring is NEW work**, not a copy-paste: `fastqc::run` is called only from the trim paths (`main.rs:882, 1188`); existing specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`, `--demux`) return before that call. `--clump_only` needs an explicit post-write call to `fastqc::run` on the output files.

### 4. Sort determinism (src/clump.rs) — **already in place, verify only**

Both plan reviewers verified: `sort_single_by_key` (`src/clump.rs:271`) and `sort_paired_by_key` (`src/clump.rs:302`) are **already stable** and **already carry a content-based tiebreaker cascade** (`key → seq → qual → id`). No conversion or new tiebreaker needed.

Action: read the two sort functions during implementation to confirm the cascade still holds; add one unit test in `clump_only.rs` that pins the property (same input twice → same output bytes) as a regression guard.

### 5. Report (src/report.rs OR src/clump_only.rs)

- New `ClumpOnlyStats` struct (see Signature).
- `write_clump_only_report(stats: &ClumpOnlyStats, txt_path: &Path) -> Result<()>`.
- Text format as in Behavior. **No JSON** — deliberately outside nf-core's `*_trimming_report.*` scan; anyone needing machine-parseable stats can grep the text report.

### 6. Filename generation (src/io.rs)

Add:
- `pub fn clumped_output_name(input: &Path, output_dir: Option<&Path>, gzip: bool) -> PathBuf` — mirrors `trimmed_output_name` at line 65-95 of `io.rs`.
- `pub fn clumped_val_output_names(r1, r2, output_dir, gzip) -> (PathBuf, PathBuf)` — mirrors the existing PE naming helper.

### 7. Output-collision pre-flight (src/main.rs)

Existing PE pre-flight (case-folded hash of prospective output paths per issue #216) applies to `--clump_only` too — add the new output names to the pre-flight set. Same code path.

### 8. Tests (unit + integration)

- **Unit tests in `src/clump_only.rs`**:
  - `test_clump_only_single_permutation`: build a small FASTQ, run `clump_only_single`, reconstitute 4-line records, assert `sort output.records == sort input.records` (byte-exact per record).
  - `test_clump_only_paired_lockstep`: build synthetic PE, run, assert R1[i]/R2[i] pair lockstep + permutation per-mate.
  - `test_clump_only_deterministic`: run twice on same input, assert md5(output1) == md5(output2).
  - `test_clump_only_no_trimming`: use a FASTQ with reads shorter than default `--length` filter's floor (20), assert none are dropped.
  - `test_clump_only_no_adapter_detection`: use a FASTQ with obvious Illumina adapter, assert output records preserve the adapter bytes (adapter stays in the sequence).
  - `test_clump_only_empty_input`: empty FASTQ → empty output + zero-stats report.
  - `test_clump_only_single_record`: 1-record input → 1-record output identical.
  - `test_clump_only_ignores_quality_flag`: `--clump_only -q 30` produces output byte-identical to `--clump_only` alone (regression guard for the silent-accept decision on `-q`/`--stringency`/`-e`).
  - `test_clump_only_normalizes_plus_line`: input with `+<header-repeat>` on line 3 emerges with bare `+` on output; document as expected behavior (Contract-scope note).
  - `test_clump_only_normalizes_crlf`: input with `\r\n` line endings emerges with `\n` endings on output.

- **Integration test in `tests/integration_clump_only.rs`** (new file):
  - Spawn built binary with `--clump_only <fixture>`, assert exit 0, decompress output, reconstitute 4-line records via `paste - - - -` equivalent, sort both sides, byte-diff.
  - PE variant.
  - Assert `<stem>_clumping_report.txt` exists and contains the expected sentinel string "Mode: --clump_only".
  - Assert `<stem>_trimming_report.txt` does **not** exist (guarantees no downstream nf-core confusion).
  - `--clump_only --dont_gzip` variant: assert output ends in `.fq` (not `.fq.gz`), still byte-identical records.
  - `--clump_only` on uBAM input: assert exit != 0 with a clear "uBAM input is not yet supported under --clump_only" message.

### 9. CI validate job (.github/workflows/ci.yml)

Add a new step to the existing `validate` job (or a new small job "Validate clump-only byte-identity"):

**Important**: use `paste - - - -` to reconstitute 4-line FASTQ records into single tab-separated lines **before** sorting. Sorting individual lines would catch dropped/added records but NOT catch records with mixed-up lines (e.g. R1's quality attached to R2's sequence — a real byte-identity bug that a per-line sort would miss).

```yaml
- name: --clump_only byte-identity (SE)
  run: |
    ./target/release/trim_galore --clump_only --cores 2 test_files/BS-seq_10K_R1.fastq.gz -o /tmp/co-se
    zcat test_files/BS-seq_10K_R1.fastq.gz | paste - - - - | sort > /tmp/co-in.sorted
    zcat /tmp/co-se/BS-seq_10K_R1_clumped.fq.gz | paste - - - - | sort > /tmp/co-out.sorted
    diff -q /tmp/co-in.sorted /tmp/co-out.sorted

- name: --clump_only byte-identity (PE)
  run: |
    ./target/release/trim_galore --clump_only --paired --cores 2 \
      test_files/BS-seq_10K_R{1,2}.fastq.gz -o /tmp/co-pe
    for mate in 1 2; do
      zcat test_files/BS-seq_10K_R${mate}.fastq.gz | paste - - - - | sort > /tmp/co-r${mate}.sorted
      zcat /tmp/co-pe/BS-seq_10K_R${mate}_clumped_${mate}.fq.gz | paste - - - - | sort > /tmp/co-r${mate}-out.sorted
      diff -q /tmp/co-r${mate}.sorted /tmp/co-r${mate}-out.sorted
    done

- name: --clump_only determinism
  run: |
    ./target/release/trim_galore --clump_only --cores 2 test_files/BS-seq_10K_R1.fastq.gz -o /tmp/co-det1
    ./target/release/trim_galore --clump_only --cores 2 test_files/BS-seq_10K_R1.fastq.gz -o /tmp/co-det2
    md5sum /tmp/co-det1/*.fq.gz /tmp/co-det2/*.fq.gz
    [ "$(md5sum /tmp/co-det1/*.fq.gz | awk '{print $1}')" = "$(md5sum /tmp/co-det2/*.fq.gz | awk '{print $1}')" ]
```

### 10. Documentation

- `docs/src/content/docs/modes/clump-only.md` — new Astro Starlight page. Sidebar entry in `docs/astro.config.mjs` under "Specialty modes".
- Cross-link from `docs/src/content/docs/performance/clumpy.md` (the existing `--clumpify` page) — brief mention that `--clump_only` exists for the lossless use case.

### 11. CHANGELOG.md

Add under `### Unreleased`:

> - **New `--clump_only` specialty mode** — reorders FASTQ records by canonical 16-mer minimizer for gzip-friendly compression, without any trimming, filtering, or adapter detection. Output records are byte-identical to input records; only file-level order changes. Composes with `--compression`, `--memory`, `--cores`, `--fastqc`, `--paired`, and `--dont_gzip`. Produces `*_clumped.fq(.gz)` output and a short `*_clumping_report.txt` (text-only, deliberately distinct from `*_trimming_report.*` so downstream pipelines don't scan an empty-of-trim-stats file). Byte-identity + determinism enforced by CI (`zcat | sort | diff` + `md5sum` cross-run). Requested in #353. FASTQ in/out only in v1; uBAM in/out is a natural follow-up.

## Efficiency

- **Time complexity**: Same as `--clumpify` — O(N log N / n_bins) per bin sort, dominated by I/O + gzip encoding.
- **Memory**: Peak = `layout.records_per_bin_target × avg_record_bytes × n_bins`. Bounded by `--memory` budget same as `--clumpify`.
- **Skipping trim work is a net perf win** vs `--clumpify` today — no adapter alignment, no quality-trim pass, no filter checks. In the ideal case, `--clump_only` should be faster than `--clumpify` on the same input.
- **No parallelism regressions**: reuses the `--cores >= 2` bin-worker pool pattern.
- **Stable sort** (if switching): Rust's `slice::sort_by_key` is O(N log N), same as unstable. Constant-factor slower but rarely measurable on the bin sizes we operate at.

## Integration

### Reads

- Input FASTQ (plain or gzip), same detection path as trim mode via `format::detect_input_format` / `fastq::FastqReader`.
- No adapter database read, no upstream sample sheet.

### Writes

- New output files: `*_clumped.fq(.gz)` (SE) or `*_clumped_{1,2}.fq(.gz)` (PE). `.gz` extension present unless `--dont_gzip` is set.
- New report file: `*_clumping_report.txt` per input file (SE). PE layout mirrors `--clumpify`'s current per-input report shape (verify at implementation time). No JSON.
- Optional FastQC output if `--fastqc` set: `*_clumped_fastqc.html` + `*_clumped_fastqc.zip` (via bundled fastqc-rust, same integration point as trim mode)

### Order relative to other steps

- Runs **before** any trim-pipeline code — dispatched from `main.rs` and returns immediately.
- FastQC hook runs **after** `clump_only_*` completes (same order as the trim path's post-trim FastQC call).

### Downstream impact

- Existing trim pipeline: **zero effect** (no shared code, no shared state).
- Existing `--clumpify` mode: **zero effect if we don't touch the sort primitives; small no-op behavior change if we switch `sort_single_by_key` to stable sort** (output would be deterministic across runs even for `--clumpify`, which is arguably a bug fix). Verify existing `--clumpify` test snapshots stay green either way.
- nf-core/rnaseq and other MultiQC parsers: **safe** — new report filename (`*_clumping_report.txt`) doesn't collide with the scanned-glob `*_trimming_report.txt`.

## Assumptions

Attributed to context Felix already established (issue #353 discussion + this planning session):

- The `sort_single_by_key` / `sort_paired_by_key` primitives in `src/clump.rs` correctly reorder records; if they're not currently stable, we make them stable (or add a tiebreaker).
- `canonical_minimizer` handles all edge cases (short seqs, all-N, ambiguous bases) as it does today for `--clumpify` — no new correctness requirements introduced.
- The bin-writer pool in `read_pairs_clumpy` / `read_single_clumpy` produces gzip-member-concatenated output that's a valid `.gz` file per RFC 1952. Reusing this same output shape for `--clump_only` inherits that guarantee.
- Users understand that "byte-identical records" refers to header + sequence + quality bytes, **not** to whitespace variations in the FASTQ delimiter lines. If input uses `+<repeat-header>` on the plus line, output does too; if input uses bare `+`, output does too. FastQ parser preserves this if the reader emits records that carry the plus-line bytes verbatim — verify.
- `--memory` values that resolve below `clumpify_min_memory_bytes(cores)` fail validation at CLI parse, not at runtime (existing behavior).
- Test fixture `test_files/BS-seq_10K_R{1,2}.fastq.gz` has enough diversity (10K reads, real BS-seq) to exercise a meaningful bin distribution.

## Validation

Failure-point-targeted checks. Each should be a test or CI step.

1. **Record byte-identity** — Sort input and output records by header; assert byte-exact match per record. Kills: sequence/quality/header corruption in the reader or writer. Location: unit test `test_clump_only_single_permutation` + CI step.
2. **Multiset preservation (no dedup / no drops / no adds)** — Assert `count(input records) == count(output records)`. Kills: silent record loss in binning, missed final-bin flush, incorrect end-of-stream handling. Location: unit test + integration test.
3. **Pair lockstep** — For PE: assert R1[i]'s header matches R2[i]'s header (modulo /1 vs /2 or space-separated flag) at every output position. Kills: independent reordering of R1 and R2 (would break downstream mapping). Location: unit test `test_clump_only_paired_lockstep`.
4. **Cross-run determinism** — Two runs on the same input produce byte-identical `.gz` output (md5 match). Kills: unstable sort, hash-order-dependent bin dispatch, non-deterministic gzip flush ordering. Location: CI step `--clump_only determinism`.
5. **Rejection matrix** — Assert every listed incompatible flag combination fails `Cli::validate` with a clear error message. Kills: silent acceptance of a contradictory combination (e.g. `--clump_only -q 20` silently doing no quality trim). Location: unit tests per rejected flag in `src/cli.rs::tests`.
6. **Short-read preservation** — Input contains reads shorter than the default length filter (20 bp); output contains all of them. Kills: any accidental invocation of `filters::length_filter` in the clump-only path. Location: unit test `test_clump_only_no_trimming`.
7. **Adapter preservation** — Input contains reads with obvious Illumina adapter sequences; output records still have the adapter bytes. Kills: accidental adapter-detect + trim in the clump-only path. Location: unit test `test_clump_only_no_adapter_detection`.
8. **Report filename** — Assert `*_clumping_report.txt` is created and `*_trimming_report.txt` is **not**. Kills: report reuse from the trim path silently producing a nf-core-poison file. Location: integration test.
9. **Startup notice sanity** — Startup diagnostic prints "Mode: --clump_only" (or similar sentinel). Kills: silent fallthrough into a non-clump-only path when the flag is set. Location: `assert_cmd`-style integration test on stderr.

## Resolved decisions

All open design decisions have been resolved. Locked below.

### Pre-review decisions (Felix)

1. **`--dont_gzip` allowed with `--clump_only`** — diverges from `--clumpify`'s rejection. Rationale: the archival use case explicitly asks for lossless output and shouldn't be refused just because it "seems pointless" as with `--clumpify`.
2. **Within-bin ordering: stable sort only (v1)**. Cross-input-order determinism (records-in-any-input-order → identical output) is a follow-up if requested; changes the semantic guarantee.
3. **Report filename: `<stem>_clumping_report.txt`** — text-only, no JSON. Deliberately outside nf-core's `*_trimming_report.*` scan; report is short enough to grep. Verify no existing pipeline scans `*_clumping*` at implementation time.
4. **`--fastqc` unchanged**: opt-in, runs on reordered output. Report is meaningful since record contents are unchanged.
5. **PE report layout mirrors `--clumpify`'s current per-input behavior**. Determine specifics at implementation time by reading the existing clumpify path.

### Post-review decisions (dual plan-reviewer feedback incorporated)

6. **Byte-identity contract scoped to header + sequence + quality only**. Plus-line normalized to bare `+`; CRLF normalized to LF. These normalizations are codebase-wide (not `--clump_only`-specific), inherited from `src/fastq.rs`, and documented explicitly in the `--help` text and docs.
7. **CI byte-identity check uses `paste - - - -`** to reconstitute 4-line FASTQ records before sort — otherwise per-line sort would miss records where lines are mixed up.
8. **Rejection matrix additions**: `--rename`, `--discard_untrimmed`, `--consider_already_trimmed`, and uBAM input all rejected explicitly. uBAM rejection happens at `format::detect_input_format` dispatch time, not just CLI validation.
9. **`-q` / `--stringency` / `-e` silently accepted** (option ii — matches existing specialty-mode conventions). These have clap default values on non-`Option` types so `Cli::validate()` can't distinguish user-set from default without refactoring the parser pattern. Doc + `--help` explicitly note they are ignored under `--clump_only`.
10. **Sort determinism is already in place** — `sort_single_by_key` (`src/clump.rs:271`) and `sort_paired_by_key` (`src/clump.rs:302`) are already stable with a content-based tiebreaker cascade (`key → seq → qual → id`). Implementation Step 4 is now verification-only.
11. **`--basename` compatibility**: honor when set; produces `<basename>_clumped.fq(.gz)` / `<basename>_clumped_{1,2}.fq(.gz)`. Filename helpers must read `cli.basename`.
12. **Report "Compression ratio" line omitted** when input/output compression states differ (gzip↔plain), avoiding misleading numbers.

**Remaining critical ambiguities: none.** Plan is implementation-ready.

## Self-Review

**Efficiency**: No unnecessary passes. Reuses existing streaming reader, existing per-bin dispatcher, existing sort primitives. Skipping trim work is a strict win.

**Logic**: Steps in Behavior are ordered (parse → dispatch → skip auto-detect → resolve layout → stream/bin → sort → write → report → optional FastQC). No missing prerequisites — `resolve_layout` is called before the streaming loop, and the layout drives buffer sizes.

**Edge cases**: Empty FASTQ → validated in unit test 6. Single record → validated. Short reads (below length filter) → validated. Reads with adapter → validated. PE lockstep → validated. Cross-run determinism → validated by CI. All eight validation targets in the "Validation" section map directly to identified failure modes.

**Integration**: No shared mutable state with trim pipeline; no shared writer state with `--clumpify` (each mode owns its own output). Zero impact on existing byte-identity invariant (Perl 0.6.11 parity CI job). New CI byte-identity job for clump-only is additive.

**Adjustments made during review**:
- Moved `clump_only` code out of `specialty.rs` into its own `src/clump_only.rs` module (originally proposed as `specialty.rs`, then noted that the reorder logic is substantial enough to warrant separation).
- Added the **rejection matrix startup notice sanity** validation (item 9) — catches a specific failure mode where a flag combination sneaks through validation but doesn't actually take the clump-only code path.
- Called out `_trimming_report.txt` naming as a load-bearing decision (nf-core scanning), separate from the aesthetic "which name reads better" question.

**Remaining risks**:
- If `sort_single_by_key` is currently unstable and switching it to stable causes a measurable perf regression on `--clumpify`, we need Option B (secondary tiebreaker in the clump-only path only). Not a blocker; discoverable at implementation time.
- The `plus-line preservation` assumption (bare `+` vs `+<header-repeat>`) needs verification — depends on how `FastqReader::next()` stores the plus-line content. If it doesn't preserve verbatim, we may need a small reader/writer tweak. Not a blocker; verify early in implementation.
