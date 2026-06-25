# PLAN_REVIEW_A — uBAM input support (Reviewer A)

**Plan reviewed:** `/Users/fkrueger/Github/TrimGalore/plans/06252026_ubam-input-support/PLAN.md` (v1, 2026-06-25)
**Reviewer:** Reviewer A (fresh context)
**Date:** 2026-06-25

---

## Top-level verdict

The plan is structurally sound and the integration seam (a thin `RecordSource` trait at the dispatch site, keeping `parallel.rs` unchanged) is the right call. The locked decisions (opt-in `--preserve-tags`, `noodles-bam` pure-Rust, FASTQ-only output, single-file interleaved paired) are all defensible.

However, there are several **CRITICAL** problems in §3.1 (format detection magic bytes), §3.2 (sequence/quality conversion), §3.3 (paired-interleaved validation), §5.3 step 5 (the buffer-unbounded de-interleave design), §5.1 (the `noodles-bam` version pin doesn't exist), and §6.3 (the validation-ubam CI parity claim has a samtools-fastq ordering hazard). These need resolution before implementation starts; most are tractable.

Findings are organized by severity:

- **CRITICAL** — design or implementation will produce wrong code, deadlock, or fail to compile.
- **IMPORTANT** — likely-correct intent but the plan as written hides a real risk; needs an explicit decision.
- **OPTIONAL** — discussion-worthy, won't block.

---

## CRITICAL findings

### C1. `noodles-bam = "0.90"` does not exist; pin is fictional (§5.1, §8 “Library”)

`noodles-bam` 0.x semver hasn't reached anywhere near 0.90 as of 2026-06-25 (latest published is in the 0.60-ish family; check `cargo search noodles-bam` before committing). The plan asserts `noodles-bam = "0.90"` as if it's a known-good pin. Two consequences:

1. The dep line in §5.1 step 1 will fail `cargo add`.
2. The crate authority claim (“@zaeleus, MIT”) is right, the version is not. The plan should resolve the actual current latest version, pin it exactly (`=X.Y.Z`, mirroring the `fastqc-rust` posture per CLAUDE.md), and confirm the surface used (`io::Reader`, `bam::Record::name() / sequence() / quality_scores() / data()`) matches that version. Some noodles minor versions have renamed accessors (e.g. `name()` ↔ `read_name()`, `sequence()` returning a borrow vs an owned `Sequence`).

**Action:** swap “0.90” for the actual current crates.io version, exact-pinned, and re-validate every API call site in §3.2 / §5.3.

---

### C2. Magic-byte BGZF discrimination is wrong (§3.1, §5.2 step 2)

§3.1 says: classify as `UnalignedBam` when bytes 0-3 are `1F 8B 08 04` **and** the FEXTRA section contains a `BC` subfield; otherwise `FastqGz`. Two correctness problems:

1. **A plain `gzip --rsyncable` file also sets `FLG.FEXTRA` (bit 0x04)**. So `1F 8B 08 04 …` is necessary but not sufficient — and rsyncable FASTQ files are common from `pigz --rsyncable` pipelines. The plan needs to inspect the FEXTRA subfields, not just check that FEXTRA is present.

2. **`samtools view -h` / `bgzip` emit BGZF that *also* often shows up as the first frame of a BAM**, but the plan needs to verify the BGZF "BC" subfield is *plus* the BAM file magic in the **decompressed** first block (`BAM\1` = bytes `42 41 4D 01`). The 18-byte peek won't reach the BAM file magic (the first BGZF block is up to 64 KB compressed). Without decoding the first block, you cannot reliably distinguish:
   - Plain-`bgzip`ped FASTQ (BGZF wrapper, decompressed payload starts with `@`)
   - Actual BAM (BGZF wrapper, decompressed payload starts with `BAM\1`)

The plan in §3.1 currently classifies *both* as `UnalignedBam`, then would crash inside noodles when it tries to parse a non-BAM payload.

**Action:** detection needs a two-step probe: (a) is this BGZF? (b) decompress just the first BGZF block and check for `BAM\1`. `noodles_bgzf::Reader` already exposes block-at-a-time reads, so this is a few lines, not architecturally hard — but the plan as written will misclassify bgzipped FASTQ and emit a confusing noodles error instead of a clear "this looks like bgzipped FASTQ, not BAM".

---

### C3. Sequence-decode error policy will reject legitimate uBAM (§3.2 step 2)

The plan says: "decode the 4-bit-packed sequence … we only emit the ACGTN subset … any other letter → error". The BAM spec defines the 4-bit nibble alphabet as `=ACMGRSVTWYHKDBN`. The IUPAC degenerate codes (`MRWSYKVHDB`) **do appear in legitimately unaligned BAMs** — specifically:

- PacBio HiFi CCS sometimes encodes consensus ambiguity in the uBAM produced by `ccs --hifi-kinetics`.
- ONT basecallers (Dorado, Guppy) can emit `N` and occasionally other ambiguity codes for low-confidence positions.
- 10x Genomics' `cellranger` produces uBAM with `N` in the cell-barcode and UMI regions.

The plan claims degenerate bases "imply consensus / alignment, not raw reads" — that's not correct. Ambiguity in the read sequence (not the reference) is permitted in uBAM and is independent of alignment status.

If we error out, we'll reject several of the exact pipelines §1 names as motivations (10x, PacBio HiFi, ONT). The right policy is one of:
- Pass IUPAC codes through verbatim (the trimming pipeline currently treats anything not ACGT as a non-match; check `quality.rs` and `alignment.rs` — they should already handle this gracefully).
- Coerce non-ACGTN to `N` with a one-line stderr summary at end-of-run.

**Action:** pick one of the two coercion strategies; do not error.

---

### C4. Paired-interleaved de-interleave design has a memory blow-up hazard (§5.3 step 5)

The plan says `open_paired_interleaved` returns `(BamReader, BamReader)` backed by **one** background reader thread, with two `mpsc::SyncSender`s — one per side, routed by FREAD1/FREAD2 flag. This works only if the consumer drains both sides at roughly the same rate.

In practice the trim worker pool drains R1 and R2 **strictly in lockstep** (one R1, one R2, per outer-loop iteration in `read_pairs_round_robin`, parallel.rs:673-674). That's exactly the right interleaving in BAM-standard interleaved files (`R1, R2, R1, R2, …`), so a small bounded channel works.

**But the plan doesn't verify the BAM interleaving order**. If the input is grouped (`R1, R1, …, R1, R2, R2, …, R2` — which `samtools sort -n` *will* produce for "queryname-sorted") then the reader has to buffer all R1 records before any R2 record arrives. With a 10M-read uBAM, that's potentially 1-2 GB resident before the first R2 is delivered → OOM on small CI runners, and silently slow even on big iron.

**Action:** at minimum, the plan must:
- Define the contract (strict pairwise interleaving required; anything else is rejected).
- Detect the violation: if record N+1 has the same FREAD1/2 flag as record N for K consecutive records (K = small buffer, say 64), bail with a precise error pointing at `samtools collate` as the remediation.
- Document this in §3.5 edge cases and `--paired` help text.

Optionally: support `samtools sort -n` grouped order, but only with explicit `--paired-grouped` opt-in and a documented memory cost.

---

### C5. `RecordSource` trait collides with the actual reader-thread architecture (§5.3 step 4, §5.5 step 3)

The plan claims `parallel.rs` "stays reader-agnostic" via a small `RecordSource` trait. Looking at parallel.rs:

- `read_pairs_round_robin` (line 654) takes **`reader_r1: &mut FastqReader`** by concrete type.
- `read_pairs_clumpy` (line 785) same.
- The reader threads call `FastqReader::open_threaded` (lines 188, 957) — concrete-typed.

Making this trait-based requires either:
1. Genericizing `read_pairs_*` over `R: RecordSource` (two functions × 1 type param each — modest churn but it changes the public-or-not signature inside the crate).
2. Wrapping at the call site with `enum Reader { Fastq(FastqReader), Bam(BamReader) }` and having `parallel.rs` accept `enum Reader` directly (this is the path the plan's §5.3 step 2 actually sketches, but step 4 in §5.5 says "small inner trait" — those are different designs).

The plan equivocates between the two. Pick one. The enum is simpler and gives slightly worse codegen; the generic is cleaner. The current `parallel.rs` is 2336 lines and is gnarly enough that I'd vote enum-wrap; but the **plan must commit to one design** because changes to `read_pairs_round_robin` etc. for the trait path are non-trivial and currently invisible in the "files this plan touches" list (§2 still says `src/parallel.rs` is "Untouched").

**Action:** decide enum vs. trait. Add `src/parallel.rs` to the touched-files list either way. Spell out the exact signature change.

---

### C6. CI `validation-ubam` byte-identity claim is unlikely to hold (§6.3, §9)

The plan asserts `md5(trim_galore ubam.bam) == md5(trim_galore via samtools fastq | trim_galore ...)`. Two ways this breaks:

1. **samtools fastq emits paired reads as `/1` / `/2` suffixes by default** (controlled by `-N`), TrimGalore's path will preserve the BAM's `record.name()` without those suffixes (per §3.2). So the **read IDs differ** between the two routes — the headers won't md5-match.

2. **samtools fastq is non-deterministic on multi-threaded sort** and even single-threaded its output ordering for unpaired/orphan reads depends on internal hash ordering in some versions (the `@HD VN:` line metadata also leaks into stdout in certain modes).

To get byte-identity you need an exact matching `samtools fastq` invocation (likely `samtools fastq -n -0 /dev/stdout ubam.bam | gzip -n`) AND a fixture committed with known-stable ordering. The plan glosses over this. As written, the CI job will fail on the first run and the team will start adding `-n` here, `-N` there until it passes — at which point we've drifted from the documented invariant.

**Action:** either (a) make the validation an in-tree golden file comparison instead of a round-trip through samtools — commit a small known-good `.fq.gz` reference and md5 against it; or (b) write out the exact samtools invocation in the plan and prove it byte-matches in a spike before committing the CI step.

---

## IMPORTANT findings

### I1. Aligned-BAM detection: "any read has FUNMAP cleared" is impractical (§3.1, §3.2 step 4)

The plan errors out on the first record with FUNMAP cleared. But `BamReader::sanity_check` in §5.5 step 4 says "peek-read 1 record and assert is_unmapped()". This is a partial check: a uBAM that has **mostly** unmapped reads but a few aligned ones (e.g. user did `samtools view -b -f 4 in.bam > broken.bam` but the filter was applied wrong) will pass sanity_check and then fail mid-stream after partial output is written.

This is bad UX. Either:
- Walk the entire file in sanity_check (slow but loud and early).
- Document that the check is best-effort and accept that mid-stream failures with partial output can happen.
- Require the uBAM header to contain no `@SQ` lines (uBAM convention: no reference sequences ⇒ no alignment) — this is a O(1) check.

**Action:** prefer the `@SQ`-line absence check + lazy per-record FUNMAP. Document the choice in §3.5.

---

### I2. Quality offset assumption is wrong-direction (§3.2 step 3, §8 constraints)

The plan says: "BAM quality scores are 0–93 raw Phred; +33 gives printable ASCII 33–126 (Sanger range)". Two issues:

1. Some BAM files (notably old Illumina pipelines, PacBio HiFi `qs` re-encoded data) have Phred scores stored **already-offset** at +33 in BAM — this is a spec violation but real. The plan's blanket +33 produces nonsense ASCII (out of Sanger range, possibly 156-219). Trim Galore's quality trimmer (BWA-style) would silently produce garbage.
2. PacBio HiFi reads can legitimately have Phred values above 93 (the new "QV99" CCS reads). The plan's claim "no quality scores above 93 in legitimate data" is dated.

**Action:** add a guard: after `+33`, validate every qual byte is in `[33, 126]`. If any byte exceeds 126, emit a clear error: "BAM qual byte > 93 (likely double-offset or QV99 reads) — file `samtools view -h | head` to inspect."

---

### I3. `--preserve-tags` array (`B:`) and binary (`H:`) tag stringification under-specified (§3.2 step 5)

The plan says "Array tags (`B:…`) and binary tags (`H:…`) stringify per SAM spec." But:

- `B:i,1,2,3,4,...` array can be **arbitrarily long** (e.g. ML/MM methylation arrays in modified-base uBAMs are O(read-length) bytes per record). At 10M reads × 150 bp this adds GB to the FASTQ header. Should there be a max length / warn / refuse?
- `H:` hex-encoded byte strings have variable length too.
- The plan's example tags (`CB,UB,RX`) are all `Z:` (string) — small. The realistic 10x / methylation tags will trip §6.2 `preserve_tags_roundtrip` if the reference file uses the wrong stringification.

**Action:** define the per-type stringification explicitly (cite samtools' rules), and document a max emitted-tag-length or commit to "we emit whatever noodles renders, period".

---

### I4. Background reader-thread error-surfacing for `open_threaded` (§5.3 step 4)

`FastqReader::open_threaded` (fastq.rs:244) has a careful error path that sends `Err` through the channel on read failure. The plan's `BamReader::open_threaded` says "mirrors `FastqReader::open_threaded`" — but doesn't address:

- noodles errors are `io::Error`, not `anyhow::Error`. Wrap with `.with_context()`.
- The reader-handle join error-priority pattern in parallel.rs:301-309 has specific semantics that the new BAM reader must respect (full anyhow chain wins over the in-channel sentinel). For `open_paired_interleaved` with ONE background thread fanning to two channels, you need to be careful that BOTH channels get an error-sentinel send + the thread-join error path is wired the same way.

**Action:** write out the BAM open_threaded body as a code block in the plan, not "mirrors FastqReader". The "mirror" claim has nontrivial details.

---

### I5. Empty-BAM detection is harder than empty-FASTQ (§3.5)

The plan says "Empty BAM (no records): same error as empty FASTQ". An empty (zero-record) BAM is a valid file with a header and a BGZF EOF marker (~28 bytes minimum). FastqReader's sanity_check (fastq.rs:444-468) fails because `next_record` returns None at byte zero.

For BAM, the reader will successfully parse the header, then `next_record` returns None. The error message "seems to be completely empty" is misleading — the file isn't empty, it just has no records. The user might be looking at a header-only BAM produced by a partial workflow.

**Action:** distinguish "BAM with 0 records" from "BAM with no parseable header" in error messages. The current plan conflates them.

---

### I6. Cargo.lock churn from noodles' transitive deps unlocks reproducibility risk (§7 last bullet, §11 "remaining risks")

`noodles-bam` 0.x pulls in `noodles-sam`, `noodles-core`, `noodles-bgzf`, `noodles-util`. Each on its own 0.x cadence. Even with `=0.90.0` pin on the top-level, the transitives are loose by default. Cargo.lock will pin them, but:

- Bioconda recipe builds from source and the lockfile is what produces the bit-identity guaranteed by the reproducibility CI job.
- If the bioconda recipe doesn't carry Cargo.lock (some recipes regenerate it), the build can drift.

**Action:** verify our bioconda recipe vendors Cargo.lock; if not, this dep adds risk. Either pin all transitives explicitly (lots of typing), or document the requirement in the bioconda recipe.

---

### I7. Defer-to-follow-up safety check for `--passthrough` + BAM (§3.4, §8, §11)

The plan rejects `--passthrough` + BAM in v1 because "the three-way ID sync would need re-thinking against BAM record offsets". On read, the three-way sync uses `read_id_prefix` (fastq.rs:171) — which works on `record.id`. For BAM records converted to `FastqRecord`, the `id` is just `@{record.name()}`. The three-way sync would work identically.

So **the deferral rationale is weak**. The real questions are:
- What does the third file (passthrough) look like? A second BAM with cell barcodes? A separate FASTQ with cell barcodes? The 10x Multiome case is currently a separate I1 FASTQ — combining with BAM input is an odd mix.
- Strict pair semantics in v1 (cli.rs:561-566) — same as plain `--passthrough`.

**Action:** either reject for a more honest reason (e.g. "no demonstrated user need"), or admit it'll just work and add a test fixture.

---

### I8. Multi-input single-end + BAM detection ordering (§3.3 row 2, §5.5 step 1)

Plan §3.3: "N BAM files → each handled independently as a single-end input, same as multiple .fastq.gz files today." Plan §5.5 step 1: "call detect_input_format on each input file" — implying per-file detection.

What if the user passes a mix? `trim_galore a.fastq.gz b.bam c.fq` is currently legal under TrimGalore's "loop over single-end inputs" model. The plan doesn't say:
- Is mixed-format input supported?
- Does the output-gzip default (set from input[0]'s extension at main.rs:84) need to be re-derived per input?
- What about the per-file paired-uBAM case where one BAM gets paired and the next is single-end?

**Action:** explicitly allow or reject mixed-format input, write it into §3.3, and update the `gzip` derivation logic in main.rs if needed.

---

## OPTIONAL findings

### O1. Tag-preservation R2 path (§3.2 step 5)

For paired-uBAM with `--preserve-tags`, are tags emitted on R1 only, R2 only, or both? 10x has different tag sets per mate (`CR`/`CY` are on R1, sometimes `UB` on R2 only). The plan doesn't say. samtools fastq -T emits tags on both mates. Probably mirror that.

---

### O2. `--preserve-tags` interaction with `--rename` and ID rewriters (§3.2)

`FastqRecord::append_to_id` (fastq.rs:140) appends a suffix for `--rename`. If `--preserve-tags` runs first and adds `\tCB:Z:…` to the header, then `--rename` appends after, the result is `@id\tCB:Z:…<rename-suffix>`. That breaks `samtools fastq -T` round-trip semantics. Decide an order.

---

### O3. binary size baseline (§5.1 step 2, §6 efficiency)

The plan says "≤ 2 MB increase, contingency = feature gate". The CI `reproducibility` job is the canonical check that bit-identity holds. There's no CI step that asserts the binary-size delta — only "build twice and record". Make this an actual `cargo build --release && stat target/release/trim_galore` step in CI with a hard threshold.

---

### O4. Cargo feature-gate naming (§5.1 step 2)

If gated behind a Cargo feature, "ubam" is a fine name but consider "bam-input" for clarity vs. "noodles" (matches the upstream package family) or just "bam". The plan should pick a name and commit.

---

### O5. The `is_gzipped(path)` driver of output-gzip default (main.rs:84)

`gzip = !cli.dont_gzip && is_gzipped(cli.input[0])` (main.rs:84). `is_gzipped` checks the filename extension `.gz`. A `.bam` file isn't `.gz`, so the default output would be **plain `.fq`**, not `.fq.gz` — even though BAM is BGZF-compressed and users almost certainly expect gzipped output. The plan doesn't address this.

**Action:** treat BAM as gzip-equivalent for the output-default heuristic.

---

### O6. Test fixture recipe in §7 is incorrect

`head -40 BS-seq_10K_R1.fastq.gz | samtools import -0 - …` — this `head`s a binary file, doesn't decompress it. The recipe should be `zcat BS-seq_10K_R1.fastq.gz | head -40 | samtools import -0 -`. Will catch the implementer.

---

### O7. `--preserve-tags ALL` keyword deferral (§8)

Plan defers to follow-up. Worth checking: does samtools fastq -T support `*` or `ALL`? If so, anticipate compatibility. The samtools manpage uses `-T *` (asterisk). Choose the same convention to avoid users learning two systems.

---

## Coverage / test-gap analysis (§6)

Tests that **are missing**:

- **BGZF mid-record truncation** — only mentioned in §3.5 prose, no test. Add a unit test that truncates a known-good BAM at byte N and asserts a clean error message.
- **Orphan handling under `--retain_unpaired`** — §3.5 calls out the behavior but §6.1/§6.2 has no test for the orphan-rescue path on BAM input.
- **Mixed paired/single-end BAM input** (§3.3 row 2, multiple BAMs single-end) — no test.
- **Output-gzip default with BAM input** — see O5; if we change the default, lock it down with a test.
- **`--clumpify` + BAM** — §3.4 says "allowed" but no integration test.
- **Specialty modes + BAM** (`--clock`, `--implicon`, `--hardtrim5/3`) — §7 (integration) says they "work transparently"; verify with at least one smoke test per mode.
- **`@SQ`-line-bearing BAM rejection** — if we adopt I1's recommendation.
- **Degenerate-base BAM (IUPAC `M/R/W/S/Y/K`)** — if we adopt C3's coerce-to-N path, test it.
- **Qual > 93 / double-offset BAM** — if we adopt I2's guard.

---

## Action items (prioritized)

### Must resolve before implementation (Critical)
1. **C1** Fix the `noodles-bam` version pin — verify it exists on crates.io, then exact-pin.
2. **C2** Replace magic-byte BGZF discrimination with a "decompress first BGZF block, look for `BAM\1`" check.
3. **C3** Pick a policy for IUPAC degenerate bases (pass-through or coerce-to-N), not error.
4. **C4** Define the strict-interleaving contract for paired-uBAM, add a buffer-bound + clear error on violation.
5. **C5** Pick enum-wrap vs. trait-generic for the reader abstraction; update §2 file-touched list to include `parallel.rs`.
6. **C6** Either commit golden FASTQ references for `validation-ubam`, or write the exact `samtools fastq` invocation and prove byte-identity in a pre-merge spike.

### Should resolve before implementation (Important)
7. **I1** Use `@SQ`-line absence as the cheap aligned-BAM gate; document.
8. **I2** Guard against qual bytes > 126 after +33; clear error message.
9. **I3** Pin down B/H tag stringification rules + max emitted-tag length.
10. **I4** Write out `BamReader::open_threaded` body (don't hand-wave with "mirrors FastqReader").
11. **I5** Distinguish "no header" vs. "zero records" in error messages.
12. **I6** Confirm bioconda recipe carries Cargo.lock; if not, plan for transitive-dep risk.
13. **I7** Either fix the `--passthrough` + BAM deferral rationale or admit it's just untested.
14. **I8** Spell out mixed-format input semantics in §3.3.

### Worth discussing (Optional)
15. **O1–O7** as listed above.

---

## Summary

The seam choice (thin layer at dispatch, parallel.rs unchanged in spirit) is right, locked decisions are defensible, but six CRITICAL items need resolution. The largest one is **C4 (paired-interleaved memory hazard)** — a real correctness/robustness bug waiting to happen. **C2 (BGZF discrimination)** is the next-most-likely silent failure. **C1 (fictional version pin)** is trivial to fix but blocks implementation start.

Recommend: revise plan to v2 addressing C1-C6 and I1-I8, then re-review.

---

**Report written to:** `/Users/fkrueger/Github/TrimGalore/plans/06252026_ubam-input-support/PLAN_REVIEW_A.md`
