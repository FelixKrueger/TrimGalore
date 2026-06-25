# Plan Review B — uBAM input support

**Reviewer:** B (independent, no shared state with Reviewer A)
**Plan reviewed:** `plans/06252026_ubam-input-support/PLAN.md` v1, 2026-06-25
**Verdict:** **MAJOR — five Critical findings; not ready to implement as written.** The locked decisions (noodles-bam, FASTQ-out, opt-in tag preservation, interleaved paired BAM) are sound. The detection layer, dispatch refactor, and md5 round-trip validation all have load-bearing flaws that need to be re-thought before code lands.

---

## 1 · Critical findings

### B-Crit-1 — `noodles-bam = "0.90"` collides with the version pulled in by `fastqc-rust`; binary-size assumptions are wrong

PLAN §5 Step 1 / §6 / §11 — predicted delta "500 KB – 1.5 MB, acceptance ≤ 2 MB".

`Cargo.lock` already contains the full noodles tree (lines 1046–1053+): `noodles 0.88.0`, `noodles-bam 0.73.0`, `noodles-bgzf 0.35.0`, `noodles-core 0.16.0`, `noodles-sam 0.69.0`, `noodles-csi 0.42.0`. They were dragged in by `fastqc-rust 1.0.1` (a direct dependency that ships the bundled FastQC, see `Cargo.lock` `name = "fastqc-rust"` block).

Adding `noodles-bam = "0.90"` as a top-level dep does **not** reuse those crates because the requested version range is incompatible with the transitively-required `0.73.x`. Cargo will compile *both* noodles trees (`0.73`/`0.88` for fastqc-rust + `0.90` for our direct dep), doubling the relevant code-gen footprint. The binary-size measurement step in §5.1 will catch that the delta is bigger than predicted, but the plan currently has no plan for *resolving* it. Two realistic options:

1. **Pin to the version fastqc-rust uses.** Add `noodles-bam = "=0.73.0"` (matching `noodles-bgzf 0.35`, `noodles-sam 0.69`). De-duplication then makes the delta near zero. Trade: locked to whatever fastqc-rust upgrades to next; if 1.0.2 bumps to `noodles 0.89`, we're forced along.
2. **Use the umbrella `noodles` crate already in the tree.** Enable `noodles`'s `bam` feature (which `fastqc-rust` may or may not enable — needs verification). Then no new direct dep at all, and the user-facing API is `noodles::bam::*`.

Either way, the plan's claim that we get a ≤ 2 MB delta + the §11 "feature-gate is the contingency" mitigation are both based on a wrong premise. Re-spec §5 Step 1 with the version-alignment strategy spelled out.

**Secondary surprise:** noodles-bam 0.90 has `rust_version = "1.89.0"`. The project's CLAUDE.md says the crate's floor is 1.88. Either bump the project floor (and verify CI/MSRV story) or downgrade to noodles-bam 0.88, which the index lists with `rust_version = "1.88.0"`. Plan doesn't mention this at all.

### B-Crit-2 — Magic-byte detection scheme false-positives FASTQ.gz and misses real-world uBAM variants

PLAN §3.1 / §5 Step 2:

> Bytes `[0,1] == [0x1F, 0x8B]` and `[2] == 0x08` (gzip): inspect bytes 12–17 (FEXTRA + BGZF `BC` subfield). If present, return `UnalignedBam`. Otherwise `FastqGz`.

Problems:

- **The detection rule confuses BGZF with BAM.** Any BGZF-framed file (`fastq.gz` written by `bgzip`, VCF.gz, GFF.gz, gzipped FASTQs from current Illumina BCL Convert which use BGZF blocks by default since 2024) has the same `BC` subfield as a BAM. A `.fastq.gz` produced by `bgzip` would be misclassified as `UnalignedBam`, then noodles would error on the missing `BAM\1` magic at byte 0 of the decompressed stream — a confusing error for users with valid FASTQ input.
  - The real BAM magic is the four bytes `BAM\1` (`0x42 0x41 0x4D 0x01`) at the start of the *decompressed* BGZF payload, not anything in the gzip header. The plan needs to: (a) decompress the first block, (b) check for `BAM\1` at offset 0, then classify.
- **Detection rule false-negatives some real BAM files.** Old BAM files produced by tools that don't write the BGZF `BC` subfield (rare but they exist — some old picard versions, some early ONT tools) won't be classified as BAM under §3.1.
- **`detect_input_format` signature is `Result<InputFormat>` but the table caller wants four behaviors** (FastqPlain / FastqGz / UnalignedBam / fall-through-to-existing-sanity-check). The "anything else" row of §3.1 says "error, existing colorspace/empty messages take precedence" — but if `detect_input_format` returns Err, the caller never gets to `sanity_check`. Either return a fourth `Unknown` variant or document that the caller must catch the error.

Fix: peek the BGZF magic (`1F 8B 08 04` + a sanity BSIZE), decompress just the first block (BGZF blocks are typically ≤ 64 KiB; first one is much smaller), then look for `BAM\1` in the payload. The BGZF format spec guarantees the first block decompresses standalone.

### B-Crit-3 — `samtools fastq` md5-parity test (§6.3) cannot succeed without explicit flag matching

PLAN §6.3 / §9 (samtools round-trip parity row):

> md5(`trim_galore ubam_test.bam`) == md5(`samtools fastq ubam_test.bam \| trim_galore -`)

There are at least four issues that will make this fail until the plan addresses them:

1. **`samtools fastq` defaults to appending `/1` and `/2` to read names** when input has BAM_FREAD1/FREAD2 set (and the user passes `-N`, or by default depending on samtools version). To match samtools' default-default (no flags), TrimGalore must NOT add `/1`/`/2`. To match the common workflow that uses `samtools fastq -1 r1.fq -2 r2.fq`, TrimGalore would need to. The plan does not specify which behavior is chosen.
2. **`samtools fastq` writes `+\n` (no name on the `+` line)**; the plan implicitly assumes the same via `record.write_to`, which writes `+\n`. ✓ This one is fine, but worth pinning down in tests.
3. **`samtools fastq` reads stdin via `-`**, but TrimGalore has no stdin support — `grep -n "stdin\|input \"-\"" src/cli.rs` returns nothing. The §6.3 command `samtools fastq … | trim_galore -` would never work in the current binary. Either add stdin support as part of this plan or remove that pipe and use a temp file (cosmetic).
4. **`samtools fastq` does not preserve aux tags by default; with `-T`, the field order is samtools' internal aux iteration order, not user-specified order.** The plan's `--preserve-tags TAG1,TAG2,…` honours **user-specified order** (§3.2 step 5). When `--preserve-tags` is on, the md5 round-trip parity test will fail unless the plan documents that the parity claim only holds with `--preserve-tags` off, and only with very specific `samtools fastq` flags. §9's row is silent on this constraint.

The fix is to either:
- explicitly document the flag set required on the samtools side (e.g. `samtools fastq -n` to suppress `/1`/`/2`), and limit the parity claim to the no-tag case;
- or drop the md5 parity claim and replace it with: "decoded FASTQ contains identical (name, seq, qual) tuples in the same order", which is what we actually care about. md5 parity is a stronger claim than is needed.

### B-Crit-4 — Aligned-BAM check is described as "any read has FUNMAP cleared" but is checked only once at `sanity_check` (§5 Step 5)

PLAN §3.1: *"If the file detects as BAM but **any read** has `BAM_FUNMAP` (`0x4`) cleared … error out."*

PLAN §5 Step 5.4: *"`sanity_check` … for BAM, peek-read 1 record and assert `is_unmapped()`."*

These two are contradictory. Either:
- We enforce "every record must be unmapped" by checking each record as we stream — at the cost of one extra branch per record (cheap, no measurable cost). Fail on the first FUNMAP-cleared record encountered.
- Or we accept the one-record sanity-check semantics and rename §3.1 to "if the **first** read has FUNMAP cleared".

The first reading is the correct one (mixed BAMs with some aligned, some unaligned reads are a real failure mode — a user could feed a coordinate-sorted BAM that has a handful of unmapped reads at the start). The plan should commit to per-record checking and add unit-test coverage for "BAM with first read unmapped, second read aligned → error on second record". Currently §6.1 only tests `rejects_aligned_bam` (single-record case).

Scanning all records is fast: it's one flag bit comparison after parsing, and noodles already parses the flags for every record on its iteration path.

### B-Crit-5 — `RecordSource` trait at the parallel-pool dispatch site leaks BAM-specific concerns into a load-bearing perf path

PLAN §5 Step 5.3: *"Use a small inner trait `RecordSource { fn next(&mut self) -> Result<Option<FastqRecord>>; }` so we keep the existing parallel module surface and just swap the concrete type at the call site."*

Three problems:

1. **The current FASTQ reader path is hot.** `parallel.rs::read_pairs_round_robin` calls `reader_r1.next_record()` once per record at line 673, then writes batches. Wrapping that in a trait object (`Box<dyn RecordSource>`) inserts a vtable lookup per call → 84M dynamic dispatches per Buckberry-scale run. The bench-measured wins from #248 (10% wall-clock at cores=8 from amortising per-call overhead via `record.write_to`) suggest the dispatch site IS sensitive to this overhead.
   - If the trait is monomorphic (generic `<R: RecordSource>`), then `parallel.rs`'s `read_pairs_round_robin` and `read_single_round_robin` need to grow generic parameters. That's a bigger refactor than the plan suggests — and `read_pairs_round_robin` is already called from the reader-thread closure which captures `input_r1: &Path`. The generic path requires re-typing reader-thread-spawn closures throughout.
2. **`open_threaded` interaction is unspecified.** Today, the worker pool calls `FastqReader::open_threaded` which spawns a background decompression thread (`fastq.rs` lines 244–302). The plan's `BamReader::open_threaded` (§4) says it mirrors that. But the `RecordSource` trait abstraction can't carry the threading distinction without leaking it; either the trait requires the caller to manage thread spawning, or the trait wraps the threaded model behind it. Plan doesn't say which.
3. **`open_paired_interleaved` returns `(BamReader, BamReader)` (§4)** — two `BamReader`s reading from the same file. Implementation §5.3 step 5 says "one background reader thread, two `mpsc::SyncSender` … routes by flag bits." That's a fundamentally different topology than the current pool, which assumes two **independent** readers each running their own background thread. The `read_pairs_round_robin` reader thread (`parallel.rs:188-189`) opens R1 and R2 separately and pulls from each. If the paired-interleaved path collapses to one reader+two-mpsc, then either:
   - the worker pool needs to know it's reading from a fused source (so it doesn't try to spawn two reader threads), or
   - the BamReader pair internally manages a shared reader thread, which means the second `BamReader` is just a channel receiver — and then mixing it with the `RecordSource` trait creates an odd asymmetry (the "two readers" model has the *first* spawning the thread but the *second* depending on the first being polled).

Recommendation: defer the trait/refactor decision to a spike (§5 §6.1 / §11 already calls out this as the one new control-flow seam). A 1-day spike that builds the smallest end-to-end path (single-end uBAM through the existing parallel pool) would surface the right abstraction.

---

## 2 · Important findings (worth fixing but not blocking)

### B-Imp-1 — Tag-walk efficiency claim "negligible" is true per-record but the constant matters at scale

PLAN §6: *"Length cost: O(#tags-requested × #aux-fields). Typical BAM has < 20 aux fields; user-requested tags ≤ 10. Negligible."*

At 84M reads (Buckberry-scale), 10 tags × 20 aux fields × 84M = 16.8 billion comparisons. Even at 1 ns each, that's ~17 s of wall time. Worth either:
- Building a `HashMap<TagId, AuxValue>` per record (allocate-per-record cost may be worse than the linear scan)
- Or: build a `&[u8; 2] → Option<index>` lookup once from `--preserve-tags`, then for each aux field, lookup in O(1) (small linear scan over <=10 entries is fine and cache-friendly). Inner-loop becomes O(#aux-fields × 1).

The plan's "negligible" deserves a bench note before commit. Not blocking — this is an optimisation that can land in a follow-up if the bench shows it matters.

### B-Imp-2 — PacBio HiFi qualities can exceed 60; plan says "0–93" but should verify noodles behavior on out-of-range bytes

PLAN §8 Constraints: *"BAM quality scores are 0–93 raw Phred; +33 gives printable ASCII 33–126 (Sanger range). No quality scores above 93 in legitimate data."*

PacBio HiFi (CCS consensus) routinely produces Q60+ scores. Q93 = ASCII 126 = `~` which is the boundary; per the SAM spec, anything above 93 (raw byte > 93) would push past `~` into non-printable territory. The plan should explicitly validate-and-bail on bytes > 93 in `bam_record_to_fastq` rather than silently emit unprintable ASCII into the FASTQ. Add a unit test: raw Phred 94 → error.

Per `samtools fastq`'s behavior, it just adds 33 and emits whatever falls out, including non-printable characters. So strict TrimGalore behavior diverges; the md5 parity test (§6.3) would catch that — but only if the fixture includes a high-Q record, which the proposed top-10-records-from-BS-seq fixture (§5 §7) does not. Add an explicit fixture for high-Q to §6.1 unit tests.

### B-Imp-3 — All-`!` quality records can downstream-confuse quality trimming

PLAN §3.2 step 3: *"If qual is missing … emit `'!'` (Phred 0) repeated to match seq length"*

Quality trimming (`src/quality.rs`, BWA-style 3' trim) on an all-Q0 sequence will trim everything. With default `-q 20`, the read becomes length 0. Then the length filter kicks it. The output for a uBAM-with-no-qual file is effectively "everything dropped to too_short".

This is technically the right behavior (no qual = no signal = trim everything), but it's surprising. Mitigation:
- Emit a one-time warning to stderr at first missing-qual record: "`{n}/{total}` records have missing quality scores; quality trimming will discard them. Consider `-q 0` to disable quality trimming for this input."

### B-Imp-4 — Empty BAM error message inconsistency

PLAN §3.5: *"Empty BAM (no records): same error as empty FASTQ"*

For uBAM, "empty" can mean (a) zero records but a valid header (a legitimate state for `samtools view -b` with no matching reads), or (b) truncated file (BGZF EOF block present but missing). The current `sanity_check` cannot distinguish these; for FASTQ they're indistinguishable because there's no header. For BAM, the empty-but-valid-header case is real. Worth emitting different messages, even if the action is the same (bail).

### B-Imp-5 — Concurrent BAM read by `sanity_check` and the actual read path opens the file twice

PLAN §5 Step 5.4: *"`sanity_check` … for BAM, peek-read 1 record."*

For BAM, the `sanity_check` opens the file, parses the header (potentially expensive — BAM headers in real data can be MB-sized for high-read-group projects), then reads one record. Then the actual `BamReader::open` re-opens and re-parses. That's wasteful for small files (test fixtures it's fine; for big production BAMs it's seconds wasted).

Mitigation: have `sanity_check` return enough state that the subsequent read path can reuse it. Or just skip the peek check in BAM mode since `BamReader::open` already validates the header eagerly and the first-record FUNMAP check happens on the first `next_record()` anyway.

### B-Imp-6 — Detection runs *before* sanity_check, but PLAN §5 Step 5 says "After arg-parse, before sanity_check"

PLAN §5 Step 5.1: *"After arg-parse, before sanity_check: call `detect_input_format` on each input file."*

This is fine for the first input but currently `sanity_check` is called multiple times (per-pair at `main.rs:243-245`, per-file at `main.rs:283`). Each of those sites needs `detect_input_format` plumbed through. The plan only mentions the first one. Also: `is_gzipped()` in `naming::is_gzipped(&cli.input[0])` at `main.rs:84` is used to decide output gzip mode — for BAM input, this returns true (BAM file extension isn't `.gz`, but neither is it the wanted source-of-truth). Plan needs to specify: when input is BAM, what's the output gzip default? Probably `true` (BAM is binary-compressed input, so gzipped output matches the user's storage profile) but the plan doesn't say.

### B-Imp-7 — `--preserve-tags` tag name regex `^[A-Za-z][A-Za-z0-9]$` is wrong per SAM spec

PLAN §5 Step 4.3.

Per SAM spec §1.5: "Each TAG can only appear once in one alignment line. A TAG containing lowercase letters are reserved for end users." Tag tags are 2 characters of `[A-Za-z][A-Za-z0-9]` — actually the spec says `[A-Za-z][A-Za-z0-9]`, so the regex is correct. But the plan should be aware that **lowercase first letters are user-defined tags** (perfectly valid; e.g. `xC`, `nM` etc. used by some tools), so the regex is fine but ought to be tested with a user-tag like `xC` to make sure validation doesn't choke.

Also: SAM tags `CB`, `CR`, `UB`, `UR`, `RX`, `BX` etc. — verify the plan's example list `CB,UB,RX` are all the right tags. CR (cell barcode raw) and UR (UMI raw) are also commonly used; mention them in the help text.

---

## 3 · Open / discussion points

### B-Open-1 — Multi-stream gzip FASTQ files

The detection rule in §3.1 looks at bytes 12–17 to spot BGZF's `BC` subfield. The current TrimGalore output is multi-member gzip (multiple `gzip` members concatenated; see `fastq.rs:701-742` test). If a user feeds a TrimGalore output BACK through TrimGalore, the first member's gzip header is at byte 0, the `XLEN`/`SI1`/`SI2` slot at bytes 12–17 may or may not contain BGZF's `BC` — depends on whether `gzp` writes BGZF-flavoured headers or plain gzip headers. Looking at `gzp::deflate::Gzip`, it's plain gzip, no BGZF subfield, so this is fine. But: worth a regression test "TrimGalore output → TrimGalore input" classified as `FastqGz`, not `UnalignedBam`.

### B-Open-2 — Specialty modes (`--clock`, `--implicon`) with BAM input

PLAN §3.4 says specialty modes "work transparently because they all iterate via `next_record()`." But specialty modes today directly call `FastqReader::open(input)?` (see `src/specialty.rs:32`, `:75`, `:135`, `:136`, `:248`, `:249`). They don't go through the dispatch in §5 Step 5. So either:
- The dispatch refactor needs to push the `Reader` enum into specialty.rs as well, or
- Specialty + BAM input remains unsupported in v1 (which §3.4 says is supported — inconsistent).

The plan §3.4 implicitly assumes the specialty modes can be made format-aware, but doesn't list `src/specialty.rs` in §2 "Files this plan touches" as modified. Either add specialty.rs to the change list, or reject BAM input in specialty modes for v1.

### B-Open-3 — `--paired` with `1` BAM file: pre-flight collision check needs adjustment

PLAN §3.3: `--paired` with 1 BAM file → de-interleave to R1/R2.

The pre-flight at `src/main.rs:178-225` walks `cli.input.chunks(2)` to build output paths. For paired-uBAM, `cli.input` has length 1, so `chunks(2)` yields one chunk of length 1, and `chunk[1]` would panic on index. The dispatch needs a special pre-flight path for single-file paired BAM that computes the R1/R2 output names from the single input. Plan §5 Step 5 doesn't mention this.

### B-Open-4 — `noodles` lazy parsing claim

PLAN §6: *"BAM read path uses noodles' lazy record parsing (records are decoded on iteration, not eagerly loaded)."*

The current noodles-bam 0.90 `io::Reader::records()` returns an iterator of `io::Result<Record>`. The `Record` is decoded eagerly when iterating. The plan says "lazy" — that's roughly true (one record decoded at a time, not the whole file), but the word "lazy" is misleading. Re-word to "streaming" or "incremental".

### B-Open-5 — Reproducibility of the BAM fixture (§5 §7)

The plan says to commit `test_files/ubam_test.bam` as a binary blob. CI also runs the `reproducibility` job which builds the binary twice and asserts bit-identity. The BAM fixture itself is not affected, but the validation md5 comparison (§6.3) compares an output gzipped FASTQ. The `samtools fastq | gzip` output's gzip framing depends on the version of `gzip` on the runner. Use `--no-name --no-timestamp` flags or `gzip -n` to suppress timestamp/filename, OR decompress before md5. The current `validation` job (`.github/workflows/ci.yml:291-292`) uses `gzip -dc | md5sum` — which compares decompressed bytes. §6.3 step 5 should also decompress before md5; the plan says "md5-compare the two outputs" without specifying. Be explicit.

### B-Open-6 — Validation strategy gap: no test asserts that adapter trimming actually fires on BAM input

The plan §9 validates:
- Format detection
- Quality offset
- Missing qual
- Tag ordering
- Round-trip parity
- Paired interleaved

But not: "uBAM input with embedded Illumina adapters → trimming actually runs, output has adapters removed." This is the load-bearing functional path — the whole point of the feature. The §6.2 integration test mentions "mix of clean + adapter-bearing", but doesn't pin down assertions. Add an assertion: count of reads with adapter ≥ N.

---

## 4 · Validation sufficiency

The proposed validation matrix is reasonable in shape but has the following gaps:

| Gap | Severity | Suggested test |
|---|---|---|
| samtools-flag mismatch (§6.3) | Critical | Pin samtools flags or replace md5 with content-equivalence |
| Aligned-BAM detected late, not on first aligned record | Critical | New unit test: first record unmapped, second aligned → error |
| Q > 93 silent passthrough | Important | Unit test: raw Phred 94 → error |
| All-missing-qual run kills all reads | Important | Integration test: BAM with all-`0xFF` qual, default `-q 20` → all reads dropped, warning emitted |
| BAM-as-passthrough rejection | Important | CLI validation test: `--paired R1.bam R2.bam --passthrough X.bam` → error message |
| Output gzip default for BAM input | Important | Integration test: BAM-in → `.fq.gz` output by default |
| TrimGalore output round-trip (multi-member gzip) misclassified as BAM | Open | Unit test on the detection rule |
| Specialty mode (`--clock`) with BAM input | Open | If supported per §3.4: integration test; if not: validation error |
| `bgzip`'d FASTQ misclassified as BAM | Critical | Unit test: bgzip-compressed FASTQ → `FastqGz`, not `UnalignedBam` |
| Pre-flight collision check for single-file paired BAM | Open | Integration test: BAM with output-dir set, verify collision logic doesn't panic |

---

## 5 · Alternatives worth considering

### A1 — Skip the trait, use a runtime enum and pay the match cost

`enum Reader { Fastq(FastqReader), Bam(BamReader) }` with a `next_record()` method that does a `match self`. Modern Rust compilers will collapse this to a tagged dispatch that's faster than a `dyn` vtable, and the parallel pool sees the same `Reader` shape regardless. This is conceptually simpler than the `RecordSource` trait and avoids the generic-parameter blow-up in `parallel.rs`. The cost is one branch per record instead of a vtable lookup — likely lost in noise.

### A2 — Ship without `--paired`-on-single-BAM in v1

Interleaved BAM de-interleaving is the trickiest single piece of this work (single source, two outputs, flag-based routing, shared reader thread with two SyncSender outputs). Phase 1: single-end uBAM only. Phase 2: paired interleaved. Cuts ~30% of the implementation surface and removes the most novel control-flow seam. The user feedback in #316 will tell us whether paired support is urgent — probably yes for 10x Genomics workflows, probably no for ONT / HiFi which are SE by nature.

### A3 — Use stdin instead of file paths for BAM

Add stdin support (`-`) to TrimGalore in this same plan. Then `samtools fastq` → `trim_galore -` is a real workflow, and the validation comparison in §6.3 is meaningful. The stdin path also handles "I already have BAM but I want to subset first" (e.g. `samtools view -b -f 4 in.bam | trim_galore -`). Larger scope but more user value than "BAM-in directly" in isolation. Cost: stream-based sanity_check needs reworking.

### A4 — Defer tag preservation entirely to a follow-up

Tag preservation adds non-trivial complexity (regex validation, B:/H: stringification, ordering, per-record absence handling). Land BAM-in *without* tags first; ship tag preservation as `--preserve-tags` in a follow-up. The locked decision in §10 ("Resolved 2026-06-25 — tab-separated, opt-in") is for the *format* of the feature when it ships, not necessarily for v1. The plan would be smaller and ship faster.

---

## 6 · Action items, prioritized

### Critical (resolve before implementation)

1. **B-Crit-1** — Align noodles-bam version with what fastqc-rust already pulls (`=0.73` or umbrella `noodles` with `bam` feature); update §5 Step 1 binary-size acceptance gate; resolve `rust_version 1.89` vs project floor 1.88.
2. **B-Crit-2** — Replace magic-byte detection with: BGZF first-block decompression + `BAM\1` magic check. Add a fourth `InputFormat::Unknown` variant or contract that detection failure falls through cleanly to `sanity_check`.
3. **B-Crit-3** — Either pin `samtools fastq` flags explicitly for the §6.3 parity test, or downgrade the assertion from md5 parity to content-tuple parity. Document the constraint that parity holds with `--preserve-tags` off only.
4. **B-Crit-4** — Reconcile §3.1 ("any read") vs §5 Step 5.4 ("peek-read 1 record"). Commit to per-record FUNMAP checking; add a "first unmapped, second mapped" unit test.
5. **B-Crit-5** — Spike the trait/enum dispatch decision before committing to `RecordSource`. Verify the threaded reader topology for paired-interleaved is workable through the existing pool, or document that the pool grows a fused-source path.

### Important (fix before merge)

6. **B-Imp-1** — Hashmap or O(1) tag-lookup; add a perf note even if "negligible".
7. **B-Imp-2** — Validate Phred byte ≤ 93; reject Q > 93 with a clear error.
8. **B-Imp-3** — Emit a one-time stderr warning on first missing-qual record encountered.
9. **B-Imp-4** — Distinguish "empty BAM, valid header" from "truncated BAM" in error messages.
10. **B-Imp-5** — Either reuse parse state from sanity_check or skip the BAM sanity_check entirely (let BamReader's eager header parse + first-record FUNMAP check do the work).
11. **B-Imp-6** — Specify output gzip default for BAM input; route detect_input_format through ALL sanity_check call sites in main.rs.
12. **B-Imp-7** — Test tag-name validation with a user-tag (e.g. `xC`); broaden `--preserve-tags` help text to mention `CR`/`UR` as common companion tags.

### Open (discuss, may defer)

13. **B-Open-1** — Regression test for multi-member gzip TrimGalore output → classified as FastqGz, not BAM.
14. **B-Open-2** — Decide specialty modes + BAM (support or reject); update §2 file-change list accordingly.
15. **B-Open-3** — Specify single-file paired-BAM pre-flight collision logic; the `chunks(2)` walk panics.
16. **B-Open-4** — Re-word "lazy parsing" → "streaming" / "incremental".
17. **B-Open-5** — Be explicit about gzip-decompress-before-md5 step in §6.3.
18. **B-Open-6** — Add an integration test that asserts adapter trimming actually fires on BAM input.

---

## 7 · Self-review

- **Stayed independent of Reviewer A.** I have not read PLAN_REVIEW_A.md (per the dual-review protocol).
- **Verified claims against the codebase** where possible (`src/fastq.rs`, `src/parallel.rs`, `src/specialty.rs`, `Cargo.lock`).
- **Crates.io verification** — confirmed noodles-bam 0.90 exists, has rust_version 1.89, and is brought in by a different version range than fastqc-rust's existing `noodles 0.88`.
- **Skipped trying to read noodles' source** for the lazy-parsing claim; took the documented behavior at face value.

**My biggest concern is B-Crit-1 + B-Crit-2** — the version conflict and the wrong detection rule. Either alone would cause a bad first day of implementation. Both together suggest the §5 Step 1 binary-size-baseline step won't even run cleanly without picking a side on the noodles version conflict.

— Reviewer B, 2026-06-25
