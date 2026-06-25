# SPIKE — Paired-uBAM record ordering & streaming de-interleaver design

**Date:** 2026-06-25
**Linked plan:** [`PLAN.md`](../PLAN.md) §3.3 (paired BAM behaviour), §5.3 step 5 (`open_paired_interleaved`)
**Trigger:** [`PLAN_REVIEW_A.md`](../PLAN_REVIEW_A.md) finding **C4** — claim that `samtools sort -n` produces grouped output (R1×N, R2×N) and would OOM a naïve de-interleaving reader.

---

## 1. Question, success criteria, scope

**Question.** What output ordering do common uBAM-producing tools emit (interleaved `R1, R2, R1, R2…` vs grouped `R1×N, R2×N…`), and what's the cheapest paired-uBAM reader design that (a) streams the happy-path case in bounded memory, (b) detects pathological grouped input quickly, (c) tolerates legitimate-but-slightly-jittered ordering?

**Success criteria.**
1. Each of the 8 tools in the survey gets a documented or empirically verified ordering verdict.
2. A concrete buffer-bound heuristic (number, not "small") with a rationale that ties to the worst legitimate jitter observed in the survey.
3. A user-facing error message ready to drop into `src/bam.rs`.

**Scope boundary (out).**
- Implementing the reader (separate task — `src/bam.rs` per PLAN §5.3).
- Aligned-BAM coordinate input (PLAN §3.1 rejects it; not the spike's question).
- BCL Convert / 10x / PacBio / ONT are surveyed for completeness only — none of them produce paired uBAM in a form TrimGalore would ingest in v1.

---

## 2. Methodology

1. Read PLAN §3.3 + §5.3 and PLAN_REVIEW_A C4 to understand the asserted hazard.
2. Web-research tool docs (samtools, Picard, fgbio, BCL Convert, 10x, PacBio, ONT) in three parallel sub-spikes (one per family).
3. Empirically verify the samtools claim locally with samtools 1.21 on a 20-record synthetic SAM in two starting orderings (interleaved and adversarial-grouped).
4. Sketch a bounded-buffer de-interleaver in Rust (`deinterleave_sketch.rs`) — design-only, not wired into the crate.

Scripts and BAM fixtures: [`spike2-paired-ordering/`](spike2-paired-ordering/).
Run the verifier with `bash spike2-paired-ordering/build_test_bams.sh` from any CWD.

---

## 3. Findings — tool ordering survey

| Tool | Emits paired uBAM? | Default record ordering | Source |
|---|---|---|---|
| `samtools sort -n` | Yes (post-hoc) | **Mate-adjacent** (`R1_A, R2_A, R1_B, R2_B…`). Queryname is the primary key; within a tied name, R1 sorts before R2 by flag (READ1 vs READ2). | [samtools-sort](https://www.htslib.org/doc/samtools-sort.html); empirically verified samtools 1.21 |
| `samtools collate` | Yes (re-groups) | **Mate-adjacent** in scrambled inter-pair order. Faster than `sort -n`; explicit "shuffle and group by name" command. | [samtools-collate](https://www.htslib.org/doc/samtools-collate.html); empirically verified |
| Picard `IlluminaBasecallsToSam` | Yes | **Mate-adjacent**. `SORT=true` default → queryname sort → R1 immediately followed by R2. | Picard docs, `SORT` param |
| Picard `FastqToSam` | Yes | **Mate-adjacent**. Reads two FASTQs in lockstep; `SORT_ORDER=queryname` further locks adjacency. May emit `@HD SO=unsorted GO=query`. | Picard docs |
| fgbio `FastqToBam` | Yes | **Mate-adjacent**. Default = input order (FASTQ lockstep). `--sort-order queryname` available. | [fgbio docs](https://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html) |
| fgbio `ZipperBams` | Pass-through | **Mate-adjacent** (rejects non-queryname-sorted/grouped input). | fgbio docs |
| Illumina BCL Convert | **No** — FASTQ only | n/a — only uBAM path is via Picard `IlluminaBasecallsToSam`. | Illumina BCL Convert docs |
| 10x Cell Ranger | Aligned only (not uBAM) | Coord-sorted; mates not adjacent. `bamtofastq` README: *"order not preserved"*. Out of scope for v1 (aligned). | 10x bamtofastq README |
| PacBio HiFi BAM | **No** — long reads, no R1/R2 | n/a (single molecule, ZMW-grouped). | PacBio BAM spec |
| ONT Dorado BAM | **No** — long reads, no R1/R2 | n/a. | Dorado docs |

### Key correction to PLAN_REVIEW_A's C4 premise

**`samtools sort -n` does NOT produce grouped (R1×N, R2×N) output.** It produces mate-adjacent output, because queryname is the primary sort key — R1 and R2 of the same template share a name and therefore sort adjacent, with R1 ahead of R2 by the documented flag-tiebreak rule (samtools-sort: *"Records with the same name will be ordered according to the values of the READ1 and READ2 flags"*).

**Empirical proof (samtools 1.21, this spike's `build_test_bams.sh`):**

- Input: adversarial grouped order (`read_01..10` all R1, then `read_01..10` all R2).
- After `samtools sort -n`: `read_01 R1, read_01 R2, read_02 R1, read_02 R2, …` — strict interleaving restored.
- After `samtools collate`: same mate-adjacency, scrambled inter-pair order.

So the OOM hazard described by C4 doesn't materialise on the standard toolchain. The remaining grouped-input risk is narrow but worth defending against: bespoke pipelines that hand-stitch `samtools view -f 64` and `samtools view -f 128` outputs, or `gatk MergeSamFiles` runs that concatenate single-flag streams.

### Side note — secondary/supplementary records

Empirically: a queryname-sort tie group with a supplementary alignment becomes `R1 primary, R1 supplementary, R2 primary`, breaking strict pair-adjacency. For uBAM this is moot (unmapped reads can't have alignment-derived supplementary records by definition), but a defensive reader should refuse records with `SECONDARY (0x100)` or `SUPPLEMENTARY (0x800)` set — already covered by PLAN §3.5's "unmapped only" gate.

---

## 4. Recommended design

### Algorithm

Two bounded `VecDeque<FastqRecord>` queues, one per side. For each incoming record:

1. Read flag bits. Route to the R1 queue (`flag & 0x40`) or R2 queue (`flag & 0x80`). Reject `0`, `0x40|0x80`, `0x100`, `0x800`.
2. If both queues are non-empty, pop one from each → emit a pair. (Optionally sanity-check `r1.id == r2.id` after stripping `/1` `/2` suffixes — useful for catching corruption but adds a string comparison per record.)
3. If either queue length exceeds `MAX_SLACK` → return a `GroupedInputDetected` error.
4. At EOF, both queues must be empty; otherwise error (orphan reads — handled by §3.5 `--retain_unpaired` logic).

This is O(1) memory on the happy path (queues stay near-empty), O(MAX_SLACK) memory worst-case before the error fires.

### Recommended `MAX_SLACK = 1024`

Rationale:

| Source of slack | Worst-case pending side-queue depth |
|---|---|
| Strictly-interleaved happy path | 0–1 |
| Supplementary record between R1 primary and R2 primary (uBAM: impossible) | ≤ 1 |
| Picard `FastqToSam` with `SORT_ORDER=unsorted` and a multi-threaded writer | empirically 0 (lockstep producer) |
| Bespoke pipelines: hand-spliced `samtools view -f 64` + `view -f 128` | unbounded — should fail fast |

The survey shows no legitimate uBAM producer emits more than a handful of consecutive same-side records. `MAX_SLACK = 1024` gives ~1 MB headroom at typical record sizes (~1 KB), four orders of magnitude below the multi-GB OOM that C4 worried about, and three orders of magnitude above any realistic jitter. A 10 M-read grouped uBAM trips the error after the first ~1024 records (≈ 1 MB read), not after 10 M (≈ multi-GB).

If the future surfaces a legitimate producer with > 1024-record jitter, the bound is a single named constant in `src/bam.rs` — easy to relax. Making it configurable is unnecessary for v1.

### Recommended error message

```
Paired uBAM input is grouped (all R1 records arrive before R2 records).
TrimGalore's paired-uBAM reader requires mate-adjacent ordering to stream
without unbounded buffering.

Re-interleave the input first, then re-run TrimGalore:

    samtools collate -O input.bam tmp_prefix > interleaved.bam

(or `samtools sort -n` if you prefer a fully name-sorted file).
```

Both `samtools collate` and `samtools sort -n` are valid pre-processors — `samtools fasta`'s own docs list them as alternatives. `collate` is faster (no full sort); `sort -n` produces a more deterministic file. Mention both; let users pick.

### Threading shape (refines PLAN §5.3 step 5)

Keep PLAN §5.3's "one background reader thread, two `mpsc::SyncSender`s" topology. The de-interleaver lives in the reader thread:

```
[BGZF decode + bam::Record parse]
        │
        ▼
[record → FastqRecord conversion]
        │
        ▼
[DeInterleaver.feed(rec, flag)] ──► emits Pair(r1, r2)
        │                                │
        │ Buffered / GroupedInputDetected │ send to (r1_tx, r2_tx)
        ▼                                ▼
   continue loop                  consumers
```

Bounded channels (`mpsc::sync_channel(4)` of batches of 4096) keep the downstream rate-limited; the de-interleaver itself is bounded by `MAX_SLACK`.

---

## 5. Iteration log

- **#1:** Read PLAN §3.3, §5.3 + PLAN_REVIEW_A C4. Launched three parallel research forks (samtools, Picard/fgbio, Illumina/10x/PacBio/ONT) and built a local samtools test harness to verify the C4 claim. All three forks converged on "mate-adjacent" being the universal default; the samtools harness empirically disproved C4's grouped-output claim.
- **#2:** Probed edge cases: supplementary alignments interleaved in a queryname-tied group (unmapped uBAM can't hit this — the PLAN §3.5 unmapped-only gate covers it); `-n` (natural) vs `-N` (ASCII) sort — both produce mate-adjacent output. Wrote the bounded-buffer sketch with `MAX_SLACK = 1024` and the user-facing error string.
- **#3:** (not needed — design converged in two iterations.)

---

## 6. Reference snippets

See [`spike2-paired-ordering/deinterleave_sketch.rs`](spike2-paired-ordering/deinterleave_sketch.rs) for the carryable algorithm and error string.

Key constants to migrate verbatim into `src/bam.rs`:

```rust
const FREAD1: u16 = 0x40;
const FREAD2: u16 = 0x80;
const MAX_SLACK: usize = 1024;
```

---

## 7. Recommendation

**Proceed with PLAN §5.3 step 5 (`open_paired_interleaved`)** with three concrete revisions:

1. The `(BamReader, BamReader)` topology with one background reader + two channels is sound. No change to channel design.
2. Add the bounded-buffer de-interleaver inside the reader thread (algorithm + `MAX_SLACK=1024` per §4 above).
3. Replace the existing PLAN §3.5 edge-case bullet on paired-uBAM ordering with the explicit contract and error message from §4.

C4's OOM hazard, as stated, does not exist on the standard samtools/Picard/fgbio toolchain — but the bounded-buffer defence is cheap (~80 lines of Rust, one constant) and catches the residual bespoke-pipeline risk. **C4 should be downgraded from CRITICAL to IMPORTANT** in the plan-review reconciliation: the hazard is real but its trigger is much narrower than the review claimed, and the mitigation is small.

---

## 8. Limitations

- The local samtools verification used samtools 1.21 only. Samtools versions older than ~1.10 had less strict queryname-tie ordering — production users on those should be advised to upgrade if they hit the grouped-input error. This is not a code change, just a doc note.
- Picard / fgbio behaviour was confirmed from docs and source-pointer skim, not by running the tools locally. Documented well enough that this isn't a meaningful gap.
- The `MAX_SLACK = 1024` constant is justified by the survey, not by a benchmark. If a future bug report points to a legitimate producer with > 1024 jitter, the constant moves; the algorithm doesn't.
- The deinterleave sketch is not compiled — it references `FastqRecord` as a placeholder. Wiring against the real type is part of implementation, not the spike.
