# PLAN_REVIEW_B — Phase 1 Pluggable I/O formats

**Reviewer:** B (external crate risk + back-compat + CLI ergonomics + binary/repro lens)
**Plan under review:** `phase1-trimgalore-formats/PLAN.md` v1, 2026-06-25
**Mandate:** Bring a different lens from Reviewer A — focus on external-crate maturity, FASTQ byte-identity stress-testing of the §5 step 1 refactor, CLI discoverability, and the "negligible / free I/O" claims.

---

## TL;DR

The plan has good bones — symmetric `RecordSource`/`RecordSink`, spike-gated, feature-gated. But it depends on **two brand-new 0.x crates with unfavourable risk profiles**, claims a "pure refactor" in §5 step 1 that crosses a non-obvious seam in `parallel.rs`, and quietly proposes a paired-uBAM output shape that disagrees with every industry tool listed in TrimGalore's own §bam.rs prior art. Several "negligible" / "≤ 2 MB" / "free I/O" assertions are unsupported. Spike A1/A2/A3 — fine. But the plan needs three additional pre-commit gates: a **transitive-dep audit on the binseq/mim crates**, a **paired-uBAM-shape decision**, and a **back-compat assertion that goes deeper than `gzip -dc | md5`**.

---

## 1 · External crate risk surface (the lens this reviewer was asked for)

### 1.1 BINSEQ crate (`binseq` v0.9.0, Arc Institute / Noam Teyssier)

Findings from `crates.io` and the upstream repo (verified 2026-06-25):

- **Created 2025-04-08; 15 versions in ~9 months.** Highest-stable = 0.9.0 (published to crates.io 2026-01-24; tagged on GitHub 2026-05-05 — note the 4-month tag-vs-publish lag, which signals release-process flakiness). Total downloads: 15,979; recent: 6,106.
- **Format-spec stability is poor.** v0.7.5 release notes (Oct 2025) literally titled the release "Binseq v2" — meaning the **on-disk format** versioned. v0.8.0 release notes (Dec 2025) explicitly say "Breaking `sheader` and `xheader` API". That is **two breaking format/API changes in ~6 months**. The plan's "exact-pin per fastqc-rust convention" mitigates pinning churn but does not mitigate "we picked a format that's still settling".
- **Bus factor 1**, sole maintainer is `noamteyssier`. (1 contribution from `drbh` in v0.6.5.)
- **Transitive dep surface is significant** — the crate pulls (latest 0.9.0):
  - `zstd = "0.13.3"` with `zstdmt` — `zstd-sys` has a **`build.rs` that compiles C/asm**. This directly contradicts the plan's §9 claim "Reproducibility intact — bit-identical with same `SOURCE_DATE_EPOCH`". The TrimGalore `build.rs` is currently the only `build.rs` in tree (with explicit `SOURCE_DATE_EPOCH` handling). Adding `zstd-sys` introduces a second build.rs with platform-dependent C compilation — empirical reproducibility check required before claiming "intact".
  - `memmap2`, `bytemuck`, `sucds`, `bitnuc`, `auto_impl`, `byteorder`, `itoa`, `parking_lot`, `rand`, `paraseq`. `memmap2` adds platform-specific mmap surface (not necessarily bad but raises the "no platform-specific I/O" bar.
- **`default = ["paraseq", "anyhow"]`** — the epic explicitly parks paraseq (epic §2 out-of-scope). Phase 1 MUST set `default-features = false` and pull only the non-paraseq subset. The plan does not mention this. **§5 step 7 needs to specify `default-features = false`** explicitly, or paraseq sneaks in.

### 1.2 mim crate (`mim-index` v0.1.1, COMBINE-lab / Ragnar Groot Koerkamp)

- **Created 2026-03-05; one usable version (0.1.0 was yanked).** **45 total downloads, 15 recent.** 4 months old.
- **MSRV conflict:** the crate declares `rust-version = "1.91"`. TrimGalore's Cargo.toml declares `rust-version = "1.88"` and `CLAUDE.md` calls 1.88 the floor. Adding `mim-index` as a `[dependencies]` entry **breaks the project's stated MSRV** — must either bump TrimGalore's `rust-version` to ≥1.91 (and re-check CI's `dtolnay/rust-toolchain@stable` posture) or vendor the relevant bits. The plan does not mention this.
- **Bus factor 1**, single author. License BSD-3-Clause is fine.
- **Default features pull `needletail` + `paraseq`** — must also `default-features = false`. Plan doesn't say.
- **Transitive deps in the *non-optional* set**: `clap` 4, `ctrlc`, `tracing`, `tracing-subscriber`, `serde_cbor` 0.11 (note: `serde_cbor` is **unmaintained** per RUSTSEC — its last release was 2021 and the author archived it), `bincode` 2.0.1, `blake3`, `serde_json`, `libz-rs-sys` + `flate2 zlib-rs`. Several of these are **binary-only deps** (clap, ctrlc, tracing-subscriber) that should arguably be `optional`/dev-only in a library — but they are not.
  - **`serde_cbor` brings a RUSTSEC advisory surface.** TrimGalore already maintains a `cargo audit` workflow (recent commit `d22b4b3` suppressed RUSTSEC-2024-0436 for paste). Adding mim-index will likely introduce a new advisory to suppress or work around. Worth surfacing in the spike.
  - **`clap` 4 transitive collision risk.** TrimGalore already depends on clap 4 directly. mim-index's clap pull is likely compatible, but the version range must be vetted — a minor disagreement causes two-version build.
- **Bus factor + maturity:** A 4-month-old, 45-download crate with one usable version, on a co-maintainer who is the spec author, is the load-bearing element of the entire "free indexing" story. If it goes unmaintained, TrimGalore's `--emit-mim` flag has no upstream fix path.

### 1.3 Recommendation

Add a **pre-spike** in front of the existing three spikes:

> **Spike A0 — transitive-dep audit.** For both `binseq` v0.9.0 (`default-features = false`) and `mim-index` v0.1.1 (`default-features = false`), produce `cargo tree --no-default-features --features binseq,mim` against a throwaway TrimGalore branch. Report (a) every new C/asm `build.rs` pulled, (b) any deps with RUSTSEC advisories, (c) any deps that conflict with existing TrimGalore versions (clap 4, flate2, parking_lot, serde), (d) MSRV implied by the union of all crates, (e) total transitive crate count delta. The plan as written hand-waves "exact-pinned per fastqc-rust convention" — `fastqc-rust` brings noodles transitively and is the only existing example; adding two more crates with their own dep graphs is materially different.

---

## 2 · Back-compat (§5 step 1 "pure refactor" claim)

**Claim under test:** §5 step 1 — "this is a pure refactor; FASTQ-out byte-identity must be exact." Verified by reading `src/parallel.rs`.

### 2.1 The seam is non-trivial

`parallel.rs::run_paired_end_parallel` (cores>1 path) **does not use `FastqWriter` at all** for the main output. Each worker:
1. Builds `GzEncoder::new(&mut buf_r1, Compression::new(level))` per batch (line 336).
2. Writes records into the encoder.
3. Calls `gz.finish()` per batch, which seals one gzip member.
4. Returns the `Vec<u8>` compressed bytes to the main thread.
5. Main thread `out_r1.write_all(&r.compressed_r1)` directly to a raw `File`.

The serial path (`--cores 1`) uses `FastqWriter::create(..., gzip, cores=1, ...)` which is `BufWriter<GzEncoder<File>>` — fundamentally different.

The plan §5 step 1 says: *"Update `parallel.rs::run_*_parallel` to take `Box<dyn RecordSink>` for output (mirror of how it took `Box<dyn RecordSource>` in #317)."*

This is **architecturally larger than the input-side mirror**. The input side (#317) had one polymorphic read primitive (`next_record() -> Option<FastqRecord>`). The output side has:
- **Sequence-aware ordering invariant** (the `BTreeMap<u64, batch>` keyed by batch seq number).
- **Per-batch gzip member framing** (each worker `finish()`es its own member — this is load-bearing for RFC 1952 multi-member concat).
- **Compression-buffer ownership boundary** between worker (compresses) and main (writes).

A `Box<dyn RecordSink>::write_record` called directly in workers would either need (a) one Sink per worker (then how does main merge?), or (b) one shared Sink with locking (kills parallelism). The natural answer is (c) keep the per-worker `Vec<u8>` compression model and have the `RecordSink` trait apply only at the **per-batch finalisation seam**, with a separate `Sink::open_per_batch_buffer() -> Box<dyn Write>` style — but that's not what the plan signature in §4 shows. **The signature in §4 (`write_record` + `flush`) is the wrong shape for the parallel.rs world.**

### 2.2 Specific byte-identity risks

Even granting a correct refactor, several places could shift bytes:

- **Gzip member count.** Today, parallel writes one gzip member per worker batch (4096 records). If the new `RecordSink` model accidentally moves compression to the main thread, member count changes → md5 of *raw gzip bytes* changes. (Saved by §9 validation's `gzip -dc | md5sum` — see §3 below.)
- **Gzip-level constant.** `parallel.rs` reads `config.gzip_level` per batch. `FastqWriter::create` reads it once. If the refactor consolidates them and the level differs by one anywhere (e.g. `--cores 1` vs `--cores 4` accidentally diverging), no `gzip -dc | md5` would catch it, but downstream tooling depending on file size would notice.
- **`BufWriter::with_capacity(BUF_SIZE, ...)`** in `FastqWriter` — `BUF_SIZE = 64*1024`. The worker pool path uses `Vec::with_capacity(reads.len() * 300)`. If the refactor forces serial-path FastqWriter through the new sink shape, buffering changes could shift flush boundaries → not a byte-identity bug, but a memory-footprint regression. Worth measuring.

### 2.3 Validation: §9 is too weak

The CI assertion is `gzip -dc /tmp/tg/out.fq.gz | md5sum`. That assertion is **deliberately frame-insensitive**: it compares decompressed content, not gzip bytes. **Good** for tolerating multi-member reordering. **But:**
- It will not catch any **record-ordering** drift across the refactor (e.g. if a future worker-pool change reorders records). Today's parallel path preserves input order via the BTreeMap seq counter; the refactor must preserve that. Recommend an explicit test that asserts record order is preserved across a `cores=4` run with a 10K-record fixture (compare to `cores=1` output).
- It will not catch any **single-vs-parallel divergence** — currently the validation matrix uses default `--cores`, so a refactor that subtly diverges `--cores 1` vs `--cores 4` slips through. **Add `--cores 1` and `--cores 4` to the validation matrix as separate jobs.**

### 2.4 Recommendation

Before claiming "pure refactor" in §5 step 1:

1. **Re-design §4's RecordSink trait shape** for the parallel-worker model. A `write_record(&FastqRecord) -> Result<()>` shape doesn't fit. The trait probably needs a `make_batch_buffer() -> Box<dyn Write + Send>` + `commit_batch(seq: u64, buf: Vec<u8>) -> Result<()>` shape, or equivalent. Without this, the implementer will end up shoehorning trait-object dispatch into a hot loop that wasn't designed for it.
2. **Add `--cores 1` vs `--cores 4` cross-check to validation matrix** — assert same `gzip -dc | md5` for both, on at least one SE and one PE fixture.
3. **Add raw-gzip-bytes md5 check** for `--cores 1` to lock the serial-path output frame. Multi-core path's bytes are intentionally non-deterministic (worker scheduling) so don't lock those, but DO lock serial.

---

## 3 · CLI ergonomics & discoverability

### 3.1 Feature-gated rejection error message

§5 step 7.4: *"CLI rejects `--output-format binseq` cleanly when the `binseq` feature is disabled ('this binary was built without BINSEQ support; rebuild with --features binseq')."*

Clap's `value_enum` derive produces variants at compile time. If `OutputFormat::Binseq` is `#[cfg(feature = "binseq")]`-gated out, the variant doesn't exist → clap rejects `--output-format binseq` with `"invalid value 'binseq' for '--output-format <FORMAT>'; possible values: fastq, ubam"`. The user gets **no hint about rebuilding**. The plan's promised error message requires the `Binseq` variant to exist unconditionally and gate only the backend impl behind `cfg`. **The plan signature in §5 step 1 contradicts §5 step 7.4** unless implementer notices and keeps the variant always-defined.

**Fix:** Specify explicitly in §5 step 2 that `OutputFormat::Binseq` and `OutputFormat::UBam` variants are unconditionally compiled, and the runtime check happens after parsing.

### 3.2 `--emit-mim` rejection specificity

§3.4: `--emit-mim + --output-format <not fastq>` rejected. Good. But §3.5/§5 step 2.4 doesn't show the exact error message. Recommend: *"`--emit-mim` produces a sidecar for FASTQ output; not compatible with `--output-format ubam`. Drop one or the other."* (vs. just "incompatible flags").

### 3.3 **Paired-uBAM output naming is wrong**

§3.3: paired uBAM emits `sample_val_1.bam` + `sample_val_2.bam`. The §3.3 note says "NOT interleaved — emits two separate BAMs to keep the file shape symmetric with FASTQ output".

This conflicts with TrimGalore's own `src/bam.rs` prior art and with industry convention:
- `src/bam.rs` line 32-48: *"all standard tools (samtools sort -n / collate, Picard `FastqToSam`, fgbio `FastqToBam`) emit mate-adjacent paired BAM"* — the **input** side already commits to interleaved.
- 10x CellRanger: emits interleaved uBAM with `RG:CB:UB` tags.
- fgbio `FastqToBam`: outputs **single interleaved** BAM.
- Picard `FastqToSam`: same.
- Bismark (Phase 2 target): reads interleaved uBAM per the prior-art note.

**The asymmetry guarantees** that any user who runs `trim_galore --paired R1.fq R2.fq --output-format ubam` and pipes it back into anything BAM-aware **has to re-interleave with `samtools collate`** before it works. This defeats the "no FASTQ round-trip" goal stated in EPIC §1 (*"TrimGalore → Bismark pipeline can run with any format pair without a FASTQ round-trip in between"*).

The §11 self-review says *"Originally had paired-uBAM output as one interleaved BAM. Changed to two separate BAMs (mirrors FASTQ paired-end output shape) — symmetric naming, simpler."* The simpler-naming argument is weaker than the cross-tool-compatibility argument. **Reconsider this decision.**

**Recommendation:** Emit one interleaved BAM by default (`sample_val.bam` or `sample_paired.bam`), or accept a `--bam-output-shape {interleaved,split}` flag with `interleaved` as the default. Two separate BAMs as the only option is the wrong default for the Phase 3 integration story.

### 3.4 Validation matrix extension is round-trip-only — "tests test the test"

§5 step 9 / §9 validation:
- FASTQ: existing Perl md5 baseline (good).
- BINSEQ: *"round-trip through `BinseqReader` and assert record-count parity"*.
- uBAM: *"`samtools view -c` count parity"* (uses external samtools — fine for CI).
- mim: *"verify the sidecar exists and is parseable by the mim crate's verify entry-point"*.

The BINSEQ and mim assertions are **self-referential**: the writer + the reader are both from the same crate. If both have the same encoding bug, the round-trip passes. **Record-count parity** is a weak invariant — it doesn't catch sequence-byte corruption, qual-byte corruption, off-by-one trims, or order shuffles.

**Recommendation:**
1. For BINSEQ: after round-trip, assert the **decoded (id, seq, qual) tuples** match the original FASTQ tuples, not just the count. This is what the existing #317 integration tests do for uBAM and the plan's own §5 step 5 promises — extend that pattern to the CI matrix, not just the cargo tests.
2. For uBAM: use `samtools view` + diff on the SEQ/QUAL/QNAME fields, not just `view -c`. (Or even better: `samtools fastq` the uBAM and md5 it against the FASTQ baseline.)
3. For mim: build a tiny verifier test in the cargo suite that opens the sidecar with `mim-index`'s reader and asserts it can seek to record N and recover the same `(id, seq, qual)` as the FASTQ reader sees at offset N.

### 3.5 Specialty modes — coverage gap

§7 says: *"Specialty modes (`--hardtrim5/3` / `--clock` / `--implicon`)... will also need to honor `--output-format` — Step 6 extends their writer construction."*

§5 step 6 is one sentence on this. **The four specialty modes own their output naming** (per `CLAUDE.md`: "Each mode owns its output naming"). Each currently has its own hand-rolled output path. Extending `--output-format binseq` to `--hardtrim5` means each specialty mode needs its own per-format naming logic *and* its own per-format dispatch.

**Recommendation:** Either (a) explicitly enumerate the specialty-mode-by-format matrix as out-of-scope for this phase (then `--output-format ubam --hardtrim5 30` → reject with "specialty modes emit FASTQ only in v1"), or (b) add a separate Step 6.5 that walks through each of the four specialty modes' output construction. The current "Step 6 extends their writer construction" is hand-wave-y for code that touches `src/specialty.rs` (which `CLAUDE.md` warns is non-trivial).

---

## 4 · Binary size & reproducibility

### 4.1 §9 "Reproducibility intact" — verify don't assume

The reproducibility CI job builds twice with the same `SOURCE_DATE_EPOCH` and asserts bit-identity. Adding `binseq` introduces `zstd-sys` (C compile via `build.rs`) and `bitnuc`. Adding `mim-index` introduces `blake3` (with SIMD asm), `libz-rs-sys`, `libc`.

**zstd-sys's build.rs invokes `cc` to compile the bundled zstd C source.** Reproducibility under `SOURCE_DATE_EPOCH` depends on (a) the bundled C source not having `__DATE__` / `__TIME__` / `__FILE__` macros that leak wall-clock, and (b) the system `cc` producing deterministic output. zstd upstream is generally well-behaved, but I have not personally verified it under TrimGalore's reproducibility job. **§9 cannot claim "intact" without running the doubled-build assertion with these crates pulled.** Move this from claim → spike output.

**blake3** is more concerning. Its `*-asm` features (auto-detected on x86) produce different code paths depending on the build host's CPU features. This may or may not interact with `SOURCE_DATE_EPOCH` reproducibility — depends on whether feature detection happens at compile time (deterministic) or runtime (deterministic in binary, just slower on older CPUs). Worth verifying.

**Recommendation:** Spike A3 (binary size) must also assert **reproducibility under `SOURCE_DATE_EPOCH`** with the new crates pulled. If it fails, document which crate broke it and either replace or vendor.

### 4.2 §6 "BINSEQ encoding cost should be similar to FASTQ"

§6 says: *"BINSEQ encoding — Phil's headline number is '90× parallel data-access', measured 2× alignment-time win. Per-record encoding cost should be similar to FASTQ; we trade gzip CPU for the BINSEQ packing."*

The 90× number is **read-side** (mmap + bitnuc decode). The encoding side does 4-bit packing **plus** zstd compression (per BINSEQ format). zstd at any level is slower than gzip level 1 (TrimGalore's default). Especially with `zstdmt` enabled in the crate's dep tree, threading overhead exists. **"Per-record encoding cost should be similar to FASTQ" is unsupported by the cited 90×-decode number.** Spike A1 should measure encode wall-time on a real fixture, not assume.

### 4.3 §6 "mim sidecar — Same-pass: O(1) per record"

§6 says: *"mim sidecar writing — cost depends on Spike A2 outcome. Same-pass: O(1) per record (just track an offset counter)."*

The "just track an offset counter" claim is **only true if** the underlying writer is a streaming `GzEncoder` whose decompressed-byte-position equals the byte-offset we'd want to index. Inside `flate2`'s `GzEncoder`, `inner.write(buf)` does not expose the uncompressed-byte-offset directly — you must track it externally. The bigger question is what mim's sidecar **actually indexes**: BGZF-block offsets, or compressed-stream offsets, or decompressed-record offsets. Each implies a different writer architecture.

Looking at `mim-index`'s description (*"Small index enabling multithreaded fastq.gz decompression"*), it almost certainly indexes **gzip-member starts in the compressed stream** so that workers can each decompress a member independently. **That requires the writer to know its current compressed-byte offset.** Today's parallel.rs path emits multiple gzip members per batch but the main thread concatenates them as opaque `Vec<u8>` blobs — the per-batch byte size IS known, so this is tractable, but it requires plumbing those sizes through the new `RecordSink`. The serial path uses a single-member `GzEncoder<File>` — the entire output is one member and mim's "multi-threaded decompression" payoff is zero. **`--emit-mim --cores 1` would produce a sidecar that doesn't unlock parallelism.** Document this asymmetry, or add a `--cores 1 + --emit-mim` warning.

**Spike A2 must** read mim-index's spec and identify *exactly* what it indexes. The current plan handwaves "byte offsets within the wrapped sink's output file" (§4) without committing to which kind.

---

## 5 · Other findings

### 5.1 §3.1 detection table magic-bytes ambiguity

The table says: *"First 4 bytes `BINSEQ\0` magic (or whatever the BINSEQ spec defines — verified by Spike A1 below)"*.

"BINSEQ\0" is 7 bytes, not 4. The actual magic per `binseq` 0.9.0 is a 4-byte format-version header (see `bitnuc`'s file-header definitions). **Plan should not commit to "BINSEQ\0"** even as a placeholder — Spike A1 verifies. Replace with `(per BINSEQ format spec — Spike A1)`.

Also: BINSEQ has **two file types in the latest crate** — `.bq` (block-encoded) and `.vbq` (variable-block). The plan says `.binseq` extension. **Neither of these matches the actual binseq crate's file extensions.** §10 question 5 acknowledges this is open. Spike A1 must answer it.

### 5.2 §3.5 edge cases — uBAM-out behaviour conflict

§3.5 says:
- *"uBAM output requested but input was FASTQ: produces records with no aux tags + a stderr note"*.
- *"`--preserve-tags X,Y` + `--output-format ubam` + FASTQ input: error (no source tags to preserve)"*.

§11 self-review says the second one was hardened from "silently accept" to "hard error". Good. But the first one — FASTQ-in + uBAM-out **without** `--preserve-tags` — still silently produces empty-tag BAMs. That seems inconsistent. If the user explicitly asked for uBAM output knowing the input was FASTQ, the empty-tag BAM is presumably what they want (a wrapper for downstream tools that demand BAM). OK to keep, but add the stderr note text explicitly to §5 step 3 so it's not lost.

### 5.3 §5 step 3 — uBAM-out byte-identity with input uBAM

§11 self-review acknowledges *"`@PG ID:trim_galore VN:<version>` ... means uBAM-in → uBAM-out won't md5-match the input"*. Good. But §5 step 3.5 says: *"round-trip an in-memory uBAM record through BamReader → FastqRecord → BamWriter and assert the second BAM is bit-identical to the first (in unaligned-only fields; @PG line will differ)"*.

**"Bit-identical except for some fields"** is not bit-identical. The test must define the exact field-set being compared. Suggest a helper that strips @HD/@PG/@CO lines and then bytewise-compares the remaining record stream. Otherwise the test becomes a moving target as noodles' BAM writer evolves.

### 5.4 §5 step 7 default-features

§5 step 7.1 Cargo.toml stanza:

```toml
default = ["binseq", "mim"]
binseq = ["dep:binseq"]
mim = ["dep:mim-index"]
```

Two problems:
1. `dep:binseq` enables `binseq` with **its default features** → pulls paraseq + anyhow (§1.1 above). Must be `binseq = ["dep:binseq"]` paired with `binseq = { version = "...", default-features = false, features = ["..."] }` in `[dependencies]`. The plan doesn't say.
2. `mim = ["dep:mim-index"]` — same problem with `needletail` + `paraseq` defaults.

**Fix:** Update §5 step 7.1 to show the full `[dependencies]` block with `default-features = false` on both, and to enumerate which sub-features of each are actually needed (e.g. `binseq = { ..., features = [] }` if anyhow conversion isn't needed).

### 5.5 §5 step 7.5 CI matrix expansion impact

*"CI matrix gains two new jobs: `rust-tests-slim` (default-features-off) and `rust-tests-full` (--all-features)"*.

Current CI already has: build, test, fmt, clippy, audit, reproducibility, validation. Adding two more jobs (slim + full) doubles the test-job count. Validation alone is the slowest job. **Worth quantifying CI wall-time impact** — if CI goes from 8min to 14min, that's friction on every PR. Recommend the slim job as a smoke (build + a 5-test subset) rather than full test run.

### 5.6 §8 assumption — "BINSEQ encoding is per-record stateless"

*"BINSEQ encoding is per-record stateless (no cross-record compression dictionary like gzip). If this proves false for the chosen crate, BinseqWriter may need buffering across records — Spike A1 catches this."*

**This is almost certainly false** — binseq 0.9.0 uses zstd (per Cargo.toml). zstd absolutely uses cross-record dictionary state at the block level. The block-level batching is the whole point of the format. So `BinseqWriter` will need batching internally. **This is not "if proves false" — it's "verify the block size and buffering shape in Spike A1"**. The §11 risk should be re-framed as "BinseqWriter must respect the block-level encoding boundary; design `flush()` to seal blocks correctly".

This also interacts with parallel.rs (§2 above): if each worker emits its own BINSEQ block, you get N separately-zstd-compressed blocks per batch. Whether they concatenate into a valid `.binseq` file is a format-spec question — gzip permits concat (RFC 1952); zstd permits concat (zstd magic-prefix repeats); but **whether BINSEQ-the-container permits it must be verified.**

### 5.7 #1 question is "resolved" but Q7 still open

§10 question 7: *"`--output-format ubam` with FASTQ input | Plan emits records with empty aux + a stderr note. Alternative: hard error. Felix's call."*

This was *partially* resolved in §11 self-review (the `--preserve-tags` sub-case was hardened to error), but the bare `--output-format ubam` + FASTQ-input case remains "open" in §10. Either close it (and update §10) or leave it open (and remove the §11 self-review's "tightened to hard-error" claim). Currently the two sections disagree.

---

## 6 · Action items (prioritised)

### Critical (block implementation)

1. **C-1.** Re-design `RecordSink` trait shape in §4 to fit the parallel.rs worker-pool model (per §2.1 of this review). The current `write_record(&FastqRecord) -> Result<()>` shape forces serial dispatch and conflicts with the per-batch gzip-buffering architecture. Without this, §5 step 1 cannot be a "pure refactor".
2. **C-2.** Decide paired-uBAM output shape (interleaved vs split). The current plan's "two separate BAMs" disagrees with every industry tool (fgbio, Picard, samtools, CellRanger) and with TrimGalore's own input-side commitment to interleaved. This is locked-in behaviour after Phase 1 lands — Bismark Phase 2 will inherit. **Resolve before implementation.**
3. **C-3.** Add **Spike A0 (transitive-dep audit)** in front of A1/A2/A3 (per §1.3 of this review). Surface RUSTSEC advisories, build.rs files, MSRV deltas, dep version conflicts. Without this, the "negligible" surface-area claim is unverified.
4. **C-4.** Resolve mim-index's `rust-version = "1.91"` vs TrimGalore's MSRV-1.88 conflict. Bump TrimGalore MSRV, vendor, or drop mim from Phase 1.

### Important (must address before sign-off)

5. **I-1.** Update §5 step 7 Cargo.toml stanza to specify `default-features = false` on both `binseq` and `mim-index` (per §5.4 of this review). Otherwise paraseq sneaks in, contradicting the epic's out-of-scope decision.
6. **I-2.** Strengthen §9 validation: round-trip BINSEQ and uBAM checks must compare `(id, seq, qual)` tuples, not just counts (per §3.4). Add `--cores 1` vs `--cores 4` cross-check on the existing FASTQ matrix (per §2.3).
7. **I-3.** Spike A3 must also assert **reproducibility under `SOURCE_DATE_EPOCH`** with the new crates pulled (per §4.1). Not just binary size.
8. **I-4.** Spike A1 must identify the **exact file magic, extension(s) `.bq`/`.vbq`/`.binseq`, and block-buffering shape** for the binseq crate (per §5.1 and §5.6). The plan currently asserts statelessness — almost certainly wrong given zstd.
9. **I-5.** Specify the specialty-mode × output-format matrix explicitly in §5 step 6 — either out-of-scope-with-clear-reject, or per-mode dispatch design (per §3.5).
10. **I-6.** Reconcile §10 question 7 with §11 self-review (per §5.7).

### Optional (worth considering but not blocking)

11. **O-1.** Add an explicit `--cores 1 + --emit-mim` warning (the sidecar's parallel-decompression payoff is zero in this case; per §4.3).
12. **O-2.** Quantify CI-matrix wall-time delta from slim+full jobs (per §5.5); consider slim as a smoke instead of full re-run.
13. **O-3.** Soften `OutputFormat` variant gating to keep the variants compile-unconditionally (per §3.1), so the "rebuild with --features" error message is actually reachable.
14. **O-4.** §5 step 3.5's "bit-identical except some fields" test needs an explicit field-comparator helper (per §5.3).
15. **O-5.** Replace "BINSEQ\0" placeholder in §3.1 with `(per BINSEQ format spec — Spike A1)` to avoid implementer cargo-culting the placeholder string (per §5.1).

---

## 7 · Out-of-scope notes for Reviewer A (this reviewer's lens, just so we don't double-bill)

Things this reviewer did NOT closely examine (deferred to Reviewer A's architectural lens):
- Detailed trait-object Drop / borrow-check shape for `Box<dyn RecordSink>`.
- Whether the three pre-work spikes (A1/A2/A3) are scoped correctly.
- Whether the symmetric RecordSource/RecordSink shape will survive the Phase 3 shared-crate decision.
- The §5 step 6 writer-factory closure idiom (`Box<dyn Fn(&Path) -> Result<Box<dyn RecordSink>>>`) — this has known lifetime gotchas with generic closures but I leave that to Reviewer A.

---

## Reviewer B summary

The plan is well-intentioned, properly spike-gated, and identifies most of the right risk axes — but it under-weights the **maturity of the two new external crates** (one barely 4 months old with 45 downloads, the other on its 2nd format-breaking spec rev in 6 months), under-specifies the **§5 step 1 refactor's trait shape** for the parallel-worker model (the largest hidden complexity), and makes a **paired-uBAM output choice** that disagrees with every cited industry tool (and with Phase 2/3 of the same epic). Resolve C-1/C-2/C-3/C-4 before implementation.
