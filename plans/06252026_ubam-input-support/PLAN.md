# PLAN — uBAM (unaligned BAM) input support

**Tracking issue:** [#316](https://github.com/FelixKrueger/TrimGalore/issues/316). Spun out of [#315](https://github.com/FelixKrueger/TrimGalore/issues/315) where uBAM was tagged as a capability ask, not a perf one.

## 1 · Goal

Accept **uBAM** (unaligned BAM — every read has `BAM_FUNMAP` / `0x4` set) as an additional input format alongside FASTQ, so users of 10x Genomics, Element AVITI, Fulcrum/fgbio, PacBio HiFi, and ONT pipelines can feed Trim Galore directly without an upstream `samtools fastq` step.

**Hard invariants preserved:**
- v2.x "single static binary, no external runtime deps" — implemented via `noodles-bam` (pure Rust). `rust-htslib` is ruled out (would link libhts).
- Output remains FASTQ-flavoured (`*_trimmed.fq.gz`, `*_val_{1,2}.fq.gz`, …). We do **not** write BAM.
- Existing FASTQ paths are byte-identical to today (this is purely additive).

## 2 · Context

### Current input architecture

- `src/fastq.rs::FastqReader::open` and `::open_threaded` are the only input entry points. Both yield owned `FastqRecord { id: String, seq: String, qual: String }` via `next_record()`.
- `src/main.rs` calls `FastqReader::sanity_check` once per input as the entry-side guard against truncated/empty/colorspace input.
- The worker pool in `src/parallel.rs` consumes `Vec<FastqRecord>` batches over an `mpsc::SyncSender` channel. The pool is reader-agnostic.
- Paired-end mode (`--paired`) takes consecutive R1/R2 file pairs; `--passthrough` adds a third optional stream with three-way `read_id_prefix` sync.

### Library choice (locked)

**`noodles` 0.88 umbrella crate with `features = ["bam"]`** (pure Rust, MIT, well-maintained by @zaeleus). Sync API used **without** the `async` feature so we don't pull `tokio`/`futures`. Crucially, `noodles 0.88`, `noodles-bam 0.73`, `noodles-bgzf 0.35`, and `noodles-sam 0.69` are **already in `Cargo.lock`** transitively via `fastqc-rust 1.0.1` — using the umbrella crate with the `bam` feature reuses that tree rather than forcing a parallel-version compile (verified during plan-review v2 incorporation). Rust toolchain floor stays at 1.88.

### Files this plan touches

- New: `src/bam.rs`, `src/format.rs`, `test_files/ubam_test.bam`
- Modified: `src/lib.rs` (module decl), `src/main.rs` (dispatch), `src/cli.rs` (`--preserve-tags`, paired-uBAM validation), `src/parallel.rs` (call-site type-parameter change — `Box<dyn RecordSource>` per Spike 1), `Cargo.toml` (dep), `.github/workflows/ci.yml` (new `validation-ubam` job)
- Untouched: `src/trimmer.rs`, `src/filters.rs`, all adapter/alignment/quality logic. The worker-pool body stays reader-agnostic; only the input-type parameter at the entry sites changes.

## 3 · Behavior

### 3.1 — Auto-detection (no flag required)

Detection is by **magic bytes**, not filename. The reasoning: users routinely rename files; filename should be a hint, not the source of truth.

| Detection step | Decision |
|---|---|
| First byte is `@` (printable ASCII) | `FastqPlain` |
| First 3 bytes are `1F 8B 08` (gzip family — includes `bgzip x.fq.gz` plain-FASTQ-in-BGZF) | **Decompress the first BGZF/gzip block via `noodles_bgzf::Reader` (or `MultiGzDecoder`) and peek 4 bytes of the decompressed payload**: if those equal `42 41 4D 01` (ASCII `BAM\1`) → `UnalignedBam`. Otherwise → `FastqGz`. |
| anything else | error, existing colorspace/empty messages take precedence |

**Why the decompress step (v3 revision per A-C2 + B-Crit-2):** `bgzip x.fq` produces a BGZF-framed FASTQ that carries the same `BC` extra-field subfield as BAM (it *is* BGZF, after all — BCL Convert and modern Illumina runs are increasingly BGZF-FASTQ rather than plain-gzip-FASTQ). The only safe discriminator is the decompressed payload — BAM starts with `BAM\1`, FASTQ starts with `@`. The first BGZF block is bounded (max ~64 KB compressed), so the cost is small and one-time. Both plan reviewers caught this independently — high-confidence the heuristic-only approach is wrong.

If the file detects as BAM but **any read** has `BAM_FUNMAP` (`0x4`) cleared — i.e. it has been aligned — error out with: `"input '{path}' is an aligned BAM; Trim Galore only accepts unaligned BAM (uBAM). Use 'samtools view -h -f 4' to filter or 'samtools fastq' to convert"`.

### 3.2 — Record conversion (`bam::Record` → `FastqRecord`)

For each BAM record:

1. **ID** — `format!("@{}", record.name())`. If `name()` is empty, error with `"BAM record at offset N has empty read name"`.
2. **Seq** — decode the 4-bit-packed sequence to ASCII via noodles' `record.sequence()` iterator. Each nibble maps `=ACMGRSVTWYHKDBN` per BAM spec. Emit ACGTN verbatim; **coerce IUPAC degenerate codes** (`R Y M K S W B D H V`) **to `N`** with a single warning to stderr at first encounter (don't spam per-record). Reject `=` (BAM "match-reference" sentinel — never appears in uBAM). Rationale (v3 revision per A-C3): PacBio HiFi, ONT Dorado, and 10x cellranger uBAMs legitimately carry IUPAC codes in some workflows; the original "error on anything non-ACGTN" policy would reject real-world data.
3. **Qual** — read raw Phred bytes from `record.quality_scores()`, add `33` to each (Sanger / Phred+33). If qual is missing (encoded as `0xFF` per byte), emit `'!'` (Phred 0) repeated to match seq length — matches `samtools fastq` behaviour. If `record.quality_scores().len() != record.sequence().len()`, error.
4. **Flags check** — error if `flag.is_reverse_complemented()` is set (un-aligned should never have this). Error if `flag.is_secondary()` or `is_supplementary()` (also impossible for uBAM).
5. **Tag preservation** — if `--preserve-tags TAG1,TAG2,…` is set, walk the BAM aux data, and for each tag in the **user-specified order**, append `\t{TAG}:{TYPE}:{VALUE}` to the header (literal tab between fields, matching `samtools fastq -T` output). Missing tags from a specific record are silently skipped (per-record). Array tags (`B:…`) and binary tags (`H:…`) stringify per SAM spec.

### 3.3 — Single-end vs paired

| `--paired`? | Input files | Behaviour |
|---|---|---|
| No | 1 BAM file | All records → single-end pipeline. Records carrying R1/R2 flags ARE warned ("uBAM with R1/R2 flags but `--paired` was not specified — treating as single-end"). |
| No | N BAM files | Each handled independently as a single-end input, same as multiple `.fastq.gz` files today. |
| Yes | **1 BAM file** | De-interleave: records with `BAM_FREAD1` (0x40) → R1 stream, `BAM_FREAD2` (0x80) → R2 stream. **Requires mate-adjacent ordering** (records of the same template arrive adjacent — `samtools sort -n`, `samtools collate`, Picard `FastqToSam`/`IlluminaBasecallsToSam`, fgbio `FastqToBam`/`ZipperBams` all comply per Spike 2 survey). Both flags simultaneously set or neither set → error. |
| Yes | 2 BAM files | Error: `"--paired with two BAM files is not supported; uBAM paired mode expects a single interleaved file"`. |
| Yes | 2 FASTQ files | Unchanged — current behaviour. |

### 3.4 — Combinations excluded in v1

- `--passthrough` with BAM input → reject at validation. (The three-way ID sync would need re-thinking against BAM record offsets.)
- `--clumpify` with BAM input → allowed; the dispatcher already parses to compute its key, and `FastqRecord` is the same shape regardless of source.
- `--hardtrim5/3 / --clock / --implicon` specialty modes with BAM input → allowed (these all loop over `next_record()` and don't care about format).

### 3.5 — Edge cases

- Empty BAM (no records): same error as empty FASTQ (`"Input file '{path}' seems to be completely empty"`).
- Truncated BGZF (EOF mid-block): noodles surfaces an `io::Error`; we propagate with `.with_context()`.
- BAM with `record.sequence().is_empty()`: error.
- BAM with `record.flags().is_unmapped() == false`: aligned BAM, error per §3.1.
- Paired-uBAM with unequal R1/R2 counts (orphans): with `--retain_unpaired`, route to unpaired sinks as today. Without it, error at the first orphan: `"R1/R2 record count mismatch at BAM offset N; orphaned reads — consider --retain_unpaired"`.
- Paired-uBAM with mates not adjacent (grouped — e.g. hand-spliced pipelines): error when the de-interleaver's per-side slack exceeds `MAX_SLACK = 1024` records. Resolved per Spike 2: standard tools (samtools/Picard/fgbio) all emit mate-adjacent already, so this is a guardrail against bespoke pipelines, not a common-case failure. Error suggests `samtools collate -O input.bam tmp > interleaved.bam` or `samtools sort -n input.bam -o interleaved.bam`.
- Tag in `--preserve-tags` doesn't exist in any record: emit one warning to stderr at first encounter, never repeat.

## 4 · Signatures

```rust
// src/format.rs — new module

pub enum InputFormat {
    FastqPlain,
    FastqGz,
    UnalignedBam,
}

/// Peek the first 18 bytes of `path` and classify the input format.
/// Filename extension is consulted as a fallback when magic bytes are
/// ambiguous (e.g. uncompressed FASTQ with non-@ first byte → would normally
/// fail the colorspace check at sanity_check; that path is unchanged).
pub fn detect_input_format(path: &Path) -> Result<InputFormat>;
```

```rust
// src/bam.rs — new module

pub struct BamReader { /* private */ }

impl BamReader {
    /// Open a uBAM file. Errors if the BAM is aligned (any record has
    /// `BAM_FUNMAP` cleared) — checked lazily on the first record read.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self>;

    /// Open with a background decompression+parsing thread, mirroring
    /// `FastqReader::open_threaded`. Returns immediately; the thread
    /// batches FastqRecords through a bounded channel.
    pub fn open_threaded<P: AsRef<Path>>(path: P) -> Result<Self>;

    /// For paired-uBAM: split a single interleaved BAM into two
    /// FastqRecord-producing readers (R1 and R2).
    pub fn open_paired_interleaved<P: AsRef<Path>>(path: P) -> Result<(Self, Self)>;

    /// Set tag-preservation policy. Default: no tags emitted.
    pub fn with_preserved_tags(self, tags: &[String]) -> Self;

    /// Read the next record. Returns `None` at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}
```

The shape mirrors `FastqReader` exactly so `main.rs` can dispatch through a small `enum Reader { Fastq(FastqReader), Bam(BamReader) }` rather than touching the worker pool.

## 5 · Implementation outline

### Step 1 — Dependency + binary-size baseline

1. Add to `Cargo.toml`:
   ```toml
   noodles = { version = "=0.88.0", default-features = false, features = ["bam"] }
   ```
   - Exact pin (`=0.88.0`) matches the version already resolved transitively via `fastqc-rust 1.0.1` (`grep -A 2 'name = "noodles"' Cargo.lock` confirms `noodles 0.88.0`, `noodles-bam 0.73.0`, `noodles-bgzf 0.35.0`, `noodles-sam 0.69.0`). Using the umbrella crate with the `bam` feature pulls the same tree, not a parallel-version one.
   - `default-features = false` keeps the `async` feature off. Verify with `cargo tree -f "{p} {f}"` that `tokio` is NOT pulled (it would be a regression on the single-binary invariant).
   - Treat any future bump as a deliberate-test event (same convention as `fastqc-rust` pin).
2. Build `--release` twice (before and after the dep is added) and record the binary-size delta. **Predicted: < 500 KB** because the bgzf/sam paths are already linked via fastqc-rust; only the bam-specific code paths and record builder add new symbols. Acceptance: **≤ 1 MB increase** (tightened from v1's ≤ 2 MB on the basis that the major transitive cost is already paid). If exceeded, gate behind a Cargo feature `ubam` (default-on but optional for slim builds).

### Step 2 — Format detection (`src/format.rs`)

1. New file `src/format.rs`. Pub-export from `src/lib.rs`.
2. `detect_input_format(path)` opens the file and classifies via a **two-stage** check (v3 revision per A-C2 / B-Crit-2):
   - **Stage A (cheap byte peek):** read first byte.
     - `@` → `FastqPlain`. Return.
     - `1F` (start of gzip family) → fall through to Stage B.
     - anything else → `Err` (caller surfaces existing colorspace/empty error).
   - **Stage B (decompress one block):** wrap the file in `noodles_bgzf::Reader` (which also handles plain-gzip framing). Read up to 4 bytes of decompressed payload.
     - If those 4 bytes are `[0x42, 0x41, 0x4D, 0x01]` (`BAM\1`) → `UnalignedBam`.
     - Otherwise → `FastqGz`. Both `bgzip x.fq` (BGZF-framed FASTQ) and ordinary `gzip x.fq` (plain gzip) decompress with the same `MultiGzDecoder` downstream, so the distinction doesn't matter further.
3. **Avoid double-decompression** — pass the already-decompressed first block forward to the reader factory rather than re-opening the file from byte 0. (A small `PeekReader` wrapper that owns the already-consumed bytes is the cleanest pattern.)
4. Unit tests covering: plain FASTQ, gzipped FASTQ, **BGZF-framed FASTQ** (i.e. `bgzip x.fq` — the load-bearing case both reviewers flagged), uBAM, random binary file (negative).

### Step 3 — `src/bam.rs` reader

1. New file `src/bam.rs`. Pub-export from `src/lib.rs`.
2. `BamReader::open` wraps `noodles_bam::io::Reader::new(BufReader::new(File::open(path)?))`. Reads the header eagerly (errors on malformed BAM header → `bail!`).
3. Implement record conversion (§3.2). Pull out a private `bam_record_to_fastq(rec: &bam::Record, tags: &[String]) -> Result<FastqRecord>` so it's directly unit-testable.
4. `open_threaded` mirrors `FastqReader::open_threaded`: background thread reads + converts, sends batches of 4096 `FastqRecord` over an `mpsc::sync_channel(4)`.
5. `open_paired_interleaved` returns `(BamReader, BamReader)` where each presents the appropriate FREAD1/FREAD2 stream from the same physical file. Implementation: one background reader thread, two `mpsc::SyncSender` (one per side), routes by flag bits. **Bounded-buffer de-interleaver per Spike 2:** each side-queue has `MAX_SLACK = 1024` records; if depth exceeds slack (i.e., a long run of records arrive without their mates), error with: `"Paired uBAM input has non-adjacent mates (>1024 records arrived on one side without a mate). Re-interleave first: 'samtools collate -O input.bam tmp > interleaved.bam' or 'samtools sort -n input.bam -o interleaved.bam'."` Use `std::sync::Arc<AtomicBool>` for cooperative shutdown.

### Step 4 — CLI surface (`src/cli.rs`)

1. New flag:
   ```
   --preserve-tags <TAGS>   Comma-separated list of BAM tags (e.g.
                            "CB,UB,RX") to append to the FASTQ header
                            on uBAM input. Tab-separated, samtools
                            -T-compatible. Ignored for FASTQ input.
   ```
2. `Cli::validate()` adds:
   - `--preserve-tags` set without any BAM input → warning ("`--preserve-tags` ignored because no BAM input is present"), not error.
   - `--passthrough` + BAM input → hard error.
   - `--paired` + exactly 1 input file → only legal if that file is BAM.
3. Tag-name validation: each tag must match `^[A-Za-z][A-Za-z0-9]$` (SAM-spec tag shape). Reject `ALL` keyword in v1 (defer to a follow-up if requested).

### Step 5 — Dispatch (`src/main.rs`)

1. After arg-parse, before sanity_check: call `detect_input_format` on each input file.
2. Build a reader enum per input:
   ```rust
   enum Reader {
       Fastq(FastqReader),
       Bam(BamReader),
   }
   ```
3. `Reader::next_record(&mut self) -> Result<Option<FastqRecord>>` delegates. The worker pool (`run_single_end_parallel`, `run_paired_end_parallel`) accepts `Box<dyn RecordSource> { fn next_record(&mut self) -> Result<Option<FastqRecord>>; }`. **Resolved per Spike 1:** trait-object dispatch lands at +0.86–1.28% vs. concrete baseline on a 500K × 150bp fixture, well under the 2% bar — the parser body (four `read_line` + three `String` allocs per record) dominates by orders of magnitude over one indirect call. The enum-of-readers alternative was statistically indistinguishable. Trait wins on ergonomics and `parallel.rs` diff size. Report at `plans/06252026_ubam-input-support/spikes/SPIKE_recordsource.md`.
4. Update `sanity_check` to be format-aware: for BAM, peek-read the first record and assert `is_unmapped()`. **In addition** (resolving the §3.1-vs-§5.5.4 contradiction surfaced by B-Crit-4), the per-record `BamReader::next_record` body MUST also check `flags().is_unmapped()` on EVERY record and error on the first violation. Mixed-aligned-BAMs (first record unmapped, later records aligned — possible from some `samtools merge` outputs) would otherwise slip through the sanity_check fast-path and silently produce wrong-trimmed output. The per-record check is one branch in an already-iterating hot loop — effectively free.

### Step 6 — Tests

#### 6.1 — Unit tests (in-tree)

- `src/format.rs::tests::detects_fastq_plain | fastq_gz | ubam | rejects_random_binary`
- `src/bam.rs::tests::record_to_fastq_id_seq_qual` — synthesise a minimal BAM record in-memory via noodles' builder, convert, assert all three fields byte-equal to expected.
- `src/bam.rs::tests::qual_offset_added` — raw Phred 30 byte → ASCII `'?'` (63).
- `src/bam.rs::tests::missing_qual_emits_bangs` — `0xFF` qual sentinel → `!` for full seq length.
- `src/bam.rs::tests::preserve_tags_tab_separated` — record with CB and UB tags, preserve order honoured.
- `src/bam.rs::tests::rejects_aligned_bam` — record with FUNMAP cleared → error.
- `src/bam.rs::tests::paired_interleaved_routes_by_flag` — 4 records (2 R1, 2 R2) → two streams, each gets its 2 records in BAM order.

#### 6.2 — Integration tests (`tests/integration_ubam.rs`)

- `single_end_ubam` — pre-recorded `test_files/ubam_test.bam` (10 reads, mix of clean + adapter-bearing), run through CLI, assert output FASTQ matches a committed reference.
- `paired_end_interleaved_ubam` — `test_files/ubam_paired_test.bam` (10 pairs, interleaved), run with `--paired`, assert R1/R2 outputs match.
- `preserve_tags_roundtrip` — `samtools view test_files/ubam_test.bam` shows tags; running with `--preserve-tags CB,UB` produces FASTQ with those tags in the header in samtools-T format.

#### 6.3 — CI validation job (`.github/workflows/ci.yml`)

- New job `validation-ubam`. Steps (v3 revision per B-Crit-3):
  1. Install `samtools` on the CI runner via apt or conda. **Build-time only** — never a runtime dep of the produced binary.
  2. `samtools fastq -n test_files/ubam_test.bam | gzip > /tmp/ref.fq.gz` — the `-n` flag suppresses the default `/1` `/2` ID-suffix-addition. Without `-n` samtools appends `/1` or `/2` to records (legacy SRA-style), which would break header parity against our TrimGalore-direct path that uses bare record names. **Do NOT use `trim_galore -` or any stdin pipe** — TrimGalore does not read from stdin; both runs must consume real files.
  3. Run with `--preserve-tags` OFF (default):
     - `trim_galore /tmp/ref.fq.gz -o /tmp/via-samtools/`
     - `trim_galore test_files/ubam_test.bam -o /tmp/via-trimgalore/`
  4. **Tier-1 assertion: content-tuple parity** (resilient to gzip framing differences). Decompress each output and assert: same record count, same per-record `(id, seq, qual)` tuples in same order. This is what semantically must hold.
  5. **Tier-2 assertion (optional, may break legitimately): md5 parity.** When `--preserve-tags` is off, the two paths should be byte-identical. If md5 ever diverges, investigate — but the failure mode is often benign gzip-framing differences (worker chunk boundaries shift between runs) rather than a real semantic divergence. Document the tier-2 assertion as informational, not gating, until we confirm it's stable.
  6. Tag preservation is validated separately (§6.2 `preserve_tags_roundtrip` integration test against a committed golden snapshot, NOT against samtools — tag order is policy-defined by the `--preserve-tags` argument, not by BAM file order, so md5-parity-with-samtools is impossible when tags are on).

### Step 7 — Test fixtures

1. `test_files/ubam_test.bam` — 10 single-end reads, synthesised from the existing `BS-seq_10K_R1.fastq.gz` top 10 records by running:
   ```
   head -40 BS-seq_10K_R1.fastq.gz | samtools import -0 - -o test_files/ubam_test.bam
   samtools index test_files/ubam_test.bam   # optional
   ```
   Commit the BAM (small, ~5 KB).
2. `test_files/ubam_paired_test.bam` — interleaved version using R1+R2 from the same fixture.
3. Document the recreate-from-source recipe in `test_files/README.md` (create if missing).

### Step 8 — Documentation

1. Update `CLAUDE.md` (project) §Architecture and §Conventions for new module.
2. Update `README.md` with a uBAM usage example.
3. CHANGELOG entry referencing #316.

## 6 · Efficiency

- BAM read path uses noodles' lazy record parsing (records are decoded on iteration, not eagerly loaded). Memory footprint is the same as the FASTQ path: one record at a time, 4096-record batches.
- Tag preservation walks the aux data once per record. Length cost: O(#tags-requested × #aux-fields). Typical BAM has < 20 aux fields; user-requested tags ≤ 10. Negligible.
- Decompression: BAM is BGZF, which is concatenated gzip blocks. We get the same single-threaded decompress cost as `.fastq.gz` — slightly different framing, same magnitude. (BGZF could in principle be parallel-decompressed; out of scope here per the paraseq-spike findings — compression dominates output, not input.)
- Binary size: noodles-bam crate is small (~80 KB), but transitively pulls noodles-sam + noodles-bgzf + noodles-core. Predicted release binary delta: 500 KB – 1.5 MB. Acceptance gate: ≤ 2 MB.

## 7 · Integration

- **Input path:** new branch in `main.rs` before the trim pipeline. All worker-pool code untouched (the channel still ships `Vec<FastqRecord>`).
- **Output path:** completely unchanged. FASTQ writers, gzip multi-member output, parallel compression — all the same.
- **Stats/reports:** counts are over `FastqRecord`s after conversion, so reports look identical whether input was FASTQ or BAM.
- **Validation CI matrix:** new `validation-ubam` job runs alongside the existing Perl-0.6.11 `validation` job. The existing job is unchanged — uBAM has no Perl baseline by definition.
- **Specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`):** work transparently because they all iterate via `next_record()`.

## 8 · Assumptions

### Decisions (locked)
- **Tag format:** opt-in via `--preserve-tags`, tab-separated, samtools `-T`-compatible. Chosen 2026-06-25 by user (Q1 in plan-writer pause).
- **Library:** `noodles` 0.88 umbrella crate with `features = ["bam"]`, pure-Rust, sync API, `default-features = false`. Exact-pinned (`=0.88.0`) to match the version already pulled transitively via `fastqc-rust 1.0.1`. `rust-htslib` is excluded by the single-binary invariant.
- **Output format:** unchanged — FASTQ in, FASTQ out, BAM in, FASTQ out. No BAM writer.
- **Paired-uBAM convention:** interleaved single file (BAM standard), de-interleaved at read time. Two-BAM-files paired mode is rejected.

### Defaults (flippable in review)
- **Default-on (no Cargo feature gate):** assumed unless binary size grows by > 2 MB. If exceeded, gate behind `ubam` Cargo feature (default-on for prebuilt releases, opt-out for source builds wanting minimal size).
- **Aligned-BAM input:** hard error, not warn-and-coerce. Rationale: silently consuming aligned data with reverse-complemented R2 sequences would produce wrong trimming results; loud is correct.
- **`--passthrough` + BAM:** rejected in v1. Could be added in a follow-up if there's demand.
- **`--preserve-tags ALL`** keyword: NOT supported in v1 (per-tag list only). Defer to follow-up.

### Constraints
- BAM spec is well-defined; we trust noodles' parsing.
- BAM quality scores are 0–93 raw Phred; +33 gives printable ASCII 33–126 (Sanger range). No quality scores above 93 in legitimate data.
- BAM seq decode normalises to `A/C/G/T/N` (v3 revision per A-C3). IUPAC degenerate bases (`R Y M K S W B D H V`) are coerced to `N` with a single warning on first encounter — they appear legitimately in PacBio HiFi, ONT Dorado, and 10x cellranger uBAMs. The BAM `=` (match-reference) nibble is rejected — it never appears in unaligned data.

## 9 · Validation

| What | How | Expected |
|---|---|---|
| Format detection — happy paths | Unit test: feed each of `*.fq`, `*.fq.gz`, `*.bam` byte-prefixes; assert the right `InputFormat` is returned | All three classified correctly |
| Aligned BAM rejected | Unit test: synthesise a single-record aligned BAM (FUNMAP cleared), call `BamReader::next_record` | Errors with the §3.1 message |
| Qual offset correctness | Unit test: raw Phred byte `30u8` in BAM → FastqRecord qual byte `'?'` (63) | Pass |
| Missing qual handled | Unit test: BAM record with all-`0xFF` qual → FastqRecord qual is `!` × seq_len | Pass |
| Preserve-tags ordering | Unit test: BAM record with `[CB:Z:X, UB:Z:Y, RX:Z:Z]` aux; preserve-tags `UB,CB` → header has `UB:Z:Y\tCB:Z:X`, RX dropped | Pass |
| samtools round-trip parity | CI `validation-ubam` job: md5(`trim_galore ubam_test.bam`) == md5(`samtools fastq ubam_test.bam \| trim_galore -`) | Equal |
| Paired interleaved | Integration test: 5 pairs in one BAM, `--paired`, assert R1/R2 outputs | Both files contain 5 records, ID-paired |
| Binary size delta | CI step that builds before/after `--release` and prints/asserts delta | ≤ 2 MB |
| Reproducibility intact | Existing CI `reproducibility` job runs unchanged with new dep | Bit-identical builds with same SOURCE_DATE_EPOCH |

## 10 · Questions or ambiguities

| Q | Priority | Disposition |
|---|---|---|
| Tag format (tab-sep / space-sep / sidecar) | Critical | **Resolved** 2026-06-25 — tab-separated, opt-in via `--preserve-tags`. |
| Trait dispatch vs. enum-of-readers in `parallel.rs` | Critical | **Resolved** via Spike 1 — `Box<dyn RecordSource>` (+0.86–1.28% vs concrete, under the 2% bar; statistically indistinguishable from enum). Report: `spikes/SPIKE_recordsource.md`. |
| Paired-uBAM ordering convention | Critical | **Resolved** via Spike 2 — all standard tools (samtools sort -n / collate, Picard, fgbio) emit mate-adjacent; A-C4's OOM-on-grouped-input premise was empirically disproved. Reader uses `MAX_SLACK = 1024` bounded buffer + loud error pointing at `samtools collate` for the hand-spliced-pipeline edge case. Report: `spikes/SPIKE_paired_ordering.md`. |
| Default-on vs feature-gated | Open | Plan assumes default-on; binary-size delta tightened to ≤ 1 MB (major transitive cost already paid via fastqc-rust). |
| Aligned BAM: hard error vs warn-and-coerce | Open | Plan assumes hard error. Per-record check (v3 revision per B-Crit-4). |
| `--preserve-tags ALL` keyword | Open | Plan defers to follow-up. If you want it in v1, ~10 LOC addition. |
| `--passthrough` + BAM | Open | Plan rejects in v1. If you want it later, three-way sync becomes BAM-offset-aware instead of header-prefix-aware. |
| BAM with R1/R2 flags but no `--paired` | Open | Plan warns and treats as SE. Alternative: error. Felix's call. |
| `samtools` requirement in CI | Open | The `validation-ubam` job installs `samtools` on the runner. This is **build-time only**, not runtime. Acceptable per invariant? |

## 11 · Self-Review

Checked before delivery:

- **Efficiency** — BAM decompression and parsing are O(records). No worse than FASTQ. Tag walk is O(tags × aux-fields), bounded small. Memory footprint unchanged.
- **Logic** — End-to-end trace: detect → reader-instantiation → background thread → channel → worker pool (unchanged) → output. No step duplicated, no step missing. Paired-interleaved de-interleave is the one new control-flow seam; covered by §5.3 step 5 + §6.1 unit test.
- **Edge cases** — empty BAM, truncated BGZF, missing qual, aligned BAM, orphans in paired, mixed paired/single, tag in flag but not in record, tag name validation — all documented in §3.5 / §3 / §8.
- **Integration** — worker pool, output writers, reports, specialty modes, reproducibility, CI all considered. The one integration risk is the `Reader` enum vs trait choice; I picked a small trait (`RecordSource`) over an enum to avoid touching `parallel.rs` signatures. Documented in §5.5 step 3.

**Adjusted during self-review:**
- Originally had `samtools` as a runtime fixture-prep dep. Moved to CI-only (test fixtures are committed as `.bam` blobs, not generated at test time).
- Originally was going to support `--preserve-tags ALL`. Realised that's a separate feature with its own tag-formatting edge cases (B-arrays, H-binary); deferred.
- Originally had `--passthrough` + BAM as "should work". On closer read of `read_pairs_round_robin` the three-way sync is built around header prefixes that are FASTQ-specific; rejecting in v1 is cleaner.

**Remaining risks:**
- `noodles` 0.88 umbrella API surface may shift before 1.0; 0.x semver means breaking changes possible on minor bumps. Exact-pinned (`=0.88.0`) à la `fastqc-rust`, treat upstream bumps as deliberate-test events.
- Binary size could surprise — measurement step in Step 1 catches that early; if it blows past 1 MB the feature-gate is the contingency.
- Bioconda recipe will need a refresh once shipped (the prebuilt-binary recipe doesn't care about the dep, but anyone building from source will need the same Cargo.lock).
- **Hand-spliced pipelines that violate mate-adjacent ordering** would trip the `MAX_SLACK = 1024` guardrail. Per Spike 2 this is not a common-tool issue, but custom shell-fu (e.g., `cat r1.bam r2.bam`) would hit it. The error message is explicit about the remediation.

### v3 changes (this revision)

Incorporated from dual plan-review (`PLAN_REVIEW_A.md` / `PLAN_REVIEW_B.md`) AND Spike 1 / Spike 2:

| Change | Section | Source finding |
|---|---|---|
| BGZF detection now decompresses first block and checks `BAM\1` magic at payload offset 0 (not heuristic-only) | §3.1, §5 step 2 | A-C2 + B-Crit-2 (agreed) |
| `noodles` umbrella crate at `=0.88.0` (shares fastqc-rust's tree); binary-size budget tightened to ≤ 1 MB; Rust 1.88 floor preserved | §2, §5 step 1, §8 | A-C1 + B-Crit-1 (agreed; verified in Cargo.lock) |
| `parallel.rs` moved from "untouched" to "modified" (single type-parameter change); `Box<dyn RecordSource>` chosen for dispatch | §2, §5 step 5.3 | A-C5 + B-Crit-5 (agreed); **resolved by Spike 1** (+0.86–1.28% vs concrete, under 2% bar) |
| Sequence decoding coerces IUPAC to N with one warning (previously errored) | §3.2 step 2, §8 Constraints | A-C3 (unique to A; real-world data) |
| Aligned-BAM check explicitly per-record (not first-record-only) | §5 step 5.4 | B-Crit-4 (unique to B; internal contradiction) |
| Validation rewritten: `samtools fastq -n` (pin flags), no stdin, content-tuple primary assertion, md5 demoted to informational | §6.3 | B-Crit-3 (unique to B; trim_galore has no stdin) |
| Paired-uBAM accepts mate-adjacent ordering (all standard tools comply); guards against hand-spliced grouped input via `MAX_SLACK = 1024` bounded buffer | §3.3, §3.5, §5 step 5.3 | A-C4 (unique to A) **disproved-but-still-guarded** by Spike 2 — A-C4's OOM premise was empirically false; the design is more defensive than v1 even though the threat was overstated |

---

## Spike Results

### Spike 1 — RecordSource dispatch overhead (2026-06-25)

**Plan step validated:** §5 Step 5 (line 183) — the trait-vs-enum choice for
the `Reader` boundary feeding `parallel.rs`. Both plan reviewers flagged this
as the one architectural question worth spiking before commit.

**Outcome:** **confirmed.** Trait-object dispatch (`Box<dyn RecordSource>`) is
≤2% slower than a concrete-type baseline at production scale.

**Key data points** (500K × 150bp plain FASTQ, release profile mirroring TrimGalore's):
- (A) concrete baseline: median 85.6–86.9 ms
- (B) `Box<dyn RecordSource>`: +0.86% to +1.28% vs. (A)
- (C) `enum Reader` + match: +0.60% to +1.19% vs. (A)
- Trait-vs-enum ordering inverted between iteration 1 and iteration 2 — the
  difference between (B) and (C) is below the noise floor.

**Plan changes:** none. §5 Step 5's pick of `RecordSource` trait over a
`Reader` enum is endorsed — the dispatch cost is invisible next to the
parser body (four `read_line` calls + three `String` allocs per record),
and the trait keeps the `parallel.rs` diff minimal.

**Report:** `plans/06252026_ubam-input-support/spikes/SPIKE_recordsource.md`

---

## Spike Results — Spike 2: Paired-uBAM record ordering (C4)

**Plan step / decision validated:** §3.3 (paired BAM behaviour) + §5.3 step 5 (`open_paired_interleaved`). Triggered by PLAN_REVIEW_A finding **C4** (multi-GB OOM hazard claim from `samtools sort -n` grouping).

**Outcome:** PLAN §5.3 step 5 is **confirmed with one revision**. PLAN_REVIEW_A C4's stated premise was **disproved**: `samtools sort -n` produces mate-adjacent output (queryname sort with R1/R2 flag tiebreak), not grouped output. Survey of all common uBAM producers (samtools sort -n, samtools collate, Picard `IlluminaBasecallsToSam`, Picard `FastqToSam`, fgbio `FastqToBam`, fgbio `ZipperBams`) found mate-adjacent ordering universally. BCL Convert / 10x / PacBio HiFi / ONT Dorado either don't emit paired uBAM at all or emit aligned BAM (out of v1 scope).

**Key data points:**
- Empirically verified with samtools 1.21: `samtools sort -n` on adversarial grouped input (`R1×10` then `R2×10`) produces strict interleaving `R1_01, R2_01, R1_02, R2_02, …`.
- samtools-sort docs: *"Records with the same name will be ordered according to the values of the READ1 and READ2 flags"* — R1 before R2 within a queryname tie.
- samtools-fasta docs list both `samtools collate` and `samtools sort -n` as valid pre-interleavers.

**Plan changes:**
1. §3.3 paired-BAM row: contract is "mate-adjacent ordering required" (not "strict interleaving" — `samtools collate` produces mate-adjacent but inter-pair-scrambled output, which is also fine).
2. §5.3 step 5: keep the "one reader thread + two `mpsc::SyncSender`" topology; add a bounded-buffer de-interleaver with `MAX_SLACK = 1024` (per-side queue depth) inside the reader thread. Algorithm + Rust sketch in the spike report.
3. §3.5 edge cases: add a "grouped paired uBAM" bullet with the exact error message (suggests both `samtools collate` and `samtools sort -n` as remediations).
4. PLAN_REVIEW_A C4 should be downgraded from **CRITICAL** to **IMPORTANT** in the next review reconciliation — the hazard is real but its trigger is much narrower than the review claimed (bespoke hand-spliced pipelines, not standard tooling).

**Report:** `plans/06252026_ubam-input-support/spikes/SPIKE_paired_ordering.md`

---

## Revision History

- **v1, 2026-06-25** — initial plan. Tag format resolved to tab-separated opt-in. Defaults documented as flippable.
- **v3, 2026-06-25** — dual plan-review (A + B) + Spike 1 (RecordSource dispatch) + Spike 2 (paired-uBAM ordering) all consolidated. 5 cheap fixes applied: BGZF decompress detection, `noodles` umbrella dep, IUPAC coerce-to-N, per-record aligned-BAM check, validation rewrite. 2 architectural choices resolved by spikes: trait dispatch chosen (Spike 1), bounded-buffer de-interleaver with `MAX_SLACK = 1024` (Spike 2 — A-C4's OOM premise was empirically disproved but the guardrail kept defensively). See §11 v3-changes table for the full mapping. Supersedes the intermediate v1.1 / v1.2 (which only documented spike outcomes without folding them into the prose).
- **v1.1, 2026-06-25** — Spike 1 (RecordSource dispatch) results appended. §5 Step 5 confirmed; no plan changes.
- **v1.2, 2026-06-25** — Spike 2 (paired-uBAM ordering / C4) results appended. C4 premise disproved empirically; §5.3 step 5 confirmed with bounded-buffer (`MAX_SLACK=1024`) addition; §3.5 edge-case bullet added.
