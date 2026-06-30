# PLAN_REVIEW_D — Phase 1 Pluggable I/O formats (v2)

**Reviewer:** D (independent, fresh context window — focus: completeness vs. v1 findings + integrity of re-scoping)
**Plan under review:** `phase1-trimgalore-formats/PLAN.md` v2, 2026-06-25
**Date:** 2026-06-26

---

## TL;DR

**Verdict: PARTIAL — close to ready but two genuine gaps remain.**

v2 cleanly resolves all five v1 CRITICAL findings — drop the trait, defer the unstable crates, switch to interleaved paired-uBAM, sharpen the rejection rules, switch validation to an external oracle. The re-scoping is principled and the §11 self-review table is honest about what changed and why.

However, the re-scoping introduces **two new issues** v1 didn't have, plus one moderate carry-over:

1. **Phase 2 inheritance is now under-defined.** EPIC v1 said Phase 2 depends on Phase 1's "dispatch pattern" / "trait stability" (EPIC §4 sub-plan table; §7 sequencing notes). With v2 explicitly dropping the trait, what shape does Bismark inherit? The EPIC needs an aligned update.
2. **Phase 4 has no objective gate.** "When binseq stabilises" / "when mim-index MSRV ≤ 1.88" are subjective. A real risk: Phase 4 lives in indefinite limbo and the deferred work is never re-litigated.
3. **A-O1 carry-over: §3.4's "all FASTQ" rule is correct but the §3.3 source-aux lookup state is "TBD at impl time"** (§4 doc comment) — that ambiguity blocks impl.

Recommend addressing #1 and #2 in v2.1 or in the EPIC before implementation begins. #3 is a §4 cleanup.

---

## 1 · Verification of v1 findings against v2 (per the mandate)

### A-C1 / B-Crit1 — RecordSink trait wrong shape — **✅ FIXED**

v2 drops the trait entirely (§0 row 1; §4 "concrete struct, NOT a trait"; §11 row 1). This is **architecturally sufficient and correct**. I verified `src/parallel.rs` — the per-worker `GzEncoder<&mut Vec<u8>>` + `mpsc::SyncSender` model that v1 reviewers correctly identified would have been violated by a `write_record` trait is **untouched** in v2 (§1 "load-bearing v2 simplification — FASTQ output keeps the per-worker-compress + concatenate model"; §5 step 4 "Zero changes to this branch"). uBAM gets its own serial branch with a concrete `BamWriter`. No abstraction = no risk of cross-format bugs (§11 logic note).

**Does dropping the trait "punt the architectural question by inviting code duplication"?** §5 step 3 ("Bodies are mostly copies … only difference is the writer type and per-record `paired_side` argument") flags this, with step 3.4 explicitly deferring a shared inner function. That's a reasonable v1 trade-off — duplication is shallow (the trim loop body is ~50 LOC; the divergence is just which writer is called), and the alternative (force a trait now) is what v1 was killed for. Acceptable.

**Risk that does survive:** future formats (Phase 4 BINSEQ; speculative SAM-text) each get a new branch, not a new trait impl (§1 "Future formats … get their own branches when added"). At N=3 formats this is fine; at N=5 it gets messy. But that's a Phase-4-onward problem, not a Phase 1 problem.

### A-C3 — Spike A2 had no success criterion — **✅ FIXED (with caveat)**

Mim is deferred (§0 row 5), spike dropped. Underlying mim-format-spec question is real, and v2 §0 explicitly acknowledges: *"The work isn't lost — BINSEQ + mim become a future phase (likely Phase 4) gated on crate-ecosystem maturity. The deferred phase can land independently without re-litigating Phase 1's architecture."* That's an honest punt rather than a sleight-of-hand resolution. **Caveat:** the Phase 4 plan, when written, will need to re-answer Spike A2's question. EPIC §4 sub-plan table row 4 should be augmented with a one-line "Phase 4 must include a mim-format-spec spike before committing to any sidecar architecture."

### A-O1 — `--preserve-tags` hard-error misfire — **✅ FIXED**

§3.4 row 5 wording: *"`--preserve-tags X,Y,...` set + **no** input is uBAM (i.e. ALL inputs are FASTQ) + `--output-format ubam`"*. The "ALL inputs are FASTQ" condition exactly matches v1's recommendation. §3.4 note at the bottom calls out the v1→v2 loosening explicitly. Good.

### B-Crit3 — mim-index MSRV conflict — **✅ FIXED**

Mim deferred (§0 row 4; §11 row 5). v2's EPIC §3 Phase-4 paragraph correctly enumerates *both* gating conditions ("mim-index 0.1.1 requires Rust 1.91; project MSRV floor is 1.88. Hard showstopper until either side moves"). **However:** "either side moves" is not an objective gate (see §2 below).

### B-Crit4 — binseq instability — **✅ FIXED (with caveat)**

Binseq deferred (§0 row 4; §11 row 4). EPIC §3 enumerates the specific concerns (15 versions in 9 months, 2 breaking format revs, zstd-sys C build.rs breaks reproducibility). Same caveat as B-Crit3 — "when binseq stabilises" is not objective.

### B-Important #5 — validation self-referential — **✅ FIXED**

v2 §9 row "samtools content parity (the v1 validation gate)" explicitly states: *"`samtools view sample.bam` from FASTQ-out path + uBAM-out path produce identical (name, seq, qual) triples — catches encoding bugs the BamReader-only round-trip couldn't (Reviewer B's Important #5)"*. §11 row 9 echoes the change. **Is samtools external enough?** Yes — samtools is libhts-based (C, not noodles); they're separate codebases with separate decoders. The risk would be if both were Rust noodles wrappers, which they aren't. ✅

A secondary BamReader round-trip remains as a unit test (§5 step 2 unit tests + step 5 integration tests), which is also fine — it's testing within-codebase symmetry, which is a useful separate property. Defense in depth.

### Other v1 findings — quick check

- A-C2 (ParCompress vs GzEncoder confusion): **✅ MOOT** — trait dropped, refactor not happening.
- A-C4 (Spike A3 binary-size threshold gate): **✅ MOOT** — no new deps in v2 (noodles already pulled).
- A-O2 (interleaved paired uBAM): **✅ FIXED** — §3.2 "ONE interleaved BAM with FREAD1/FREAD2 flag bits".
- A-O3 (`--demux` × output-format): **⚠️ PARTIAL** — §7 last bullet rejects `--demux` + uBAM-out "out of scope" but this isn't in §3.4's enumerated rejection table. The §3.4 table is sold as the canonical rejection list; missing demux is inconsistent. Add it.
- A-O4 (specialty modes): **✅ FIXED** — §3.4 rejects `--clock` and `--implicon`; §3.2 explicitly supports `--hardtrim5/3` with naming spec. Note: `--implicon` UMI is on R2, not the FASTQ header — minor wording quibble in §3.4 row 4 ("UMI encoding undefined for BAM") is fine.
- A-O5 (per-record trait dispatch cost): **✅ MOOT** — no trait.
- A-O6 (per-pair header source): **⚠️ PARTIAL** — §3.5 covers uBAM-in → uBAM-out header propagation but doesn't explicitly state per-pair-header logic for batches with multiple input pairs (a single TrimGalore invocation can take many pairs; each output BAM needs its own header). §4's `BamWriter::create(source_header: Option<&Header>)` signature implies this works per-writer, but it's not spelled out. Worth one line in §3.5.
- A-O7 (interleaved round-trip): **✅ FIXED** by the O2 interleaved decision.
- A-O8 (`--emit-mim` + serial): **✅ MOOT** — mim deferred.
- B-Important #6/#7: **✅ FIXED** (BINSEQ removed; Cargo features section gone).
- B-§5.1 (BINSEQ magic), B-§5.4 (default-features), B-§5.6 (zstd statelessness): **✅ MOOT** — deferred.

### Summary table

| v1 finding | Severity | v2 disposition |
|---|---|---|
| A-C1 / B-Crit1 RecordSink trait shape | CRIT | ✅ FIXED |
| A-C2 ParCompress/GzEncoder | CRIT | ✅ MOOT |
| A-C3 Spike A2 success criterion | CRIT | ✅ FIXED (caveat: re-surface in Phase 4) |
| A-C4 Spike A3 gate consequence | CRIT | ✅ MOOT |
| B-Crit2 paired uBAM = two files | CRIT | ✅ FIXED |
| B-Crit3 mim MSRV | CRIT | ✅ FIXED (but gate subjective — see §2) |
| B-Crit4 binseq instability | CRIT | ✅ FIXED (but gate subjective — see §2) |
| A-O1 preserve-tags hard error | OPEN | ✅ FIXED |
| A-O2 paired uBAM shape | OPEN | ✅ FIXED |
| A-O3 demux × format | OPEN | ⚠️ PARTIAL (§7 mentions it, §3.4 table doesn't) |
| A-O4 specialty modes | OPEN | ✅ FIXED |
| A-O5 trait dispatch cost | OPEN | ✅ MOOT |
| A-O6 per-pair header source | OPEN | ⚠️ PARTIAL (signature implies, prose doesn't say) |
| A-O7 interleaved round-trip | OPEN | ✅ FIXED |
| A-O8 emit-mim + serial | OPEN | ✅ MOOT |
| B-Important #5 self-referential validation | IMP | ✅ FIXED (samtools is genuinely external) |
| B-Important #6 BINSEQ encoding claim | IMP | ✅ MOOT |
| B-Important #7 feature-gated CLI | IMP | ✅ MOOT |

**Net: every CRITICAL is fixed. Two OPEN items remain partial.**

---

## 2 · NEW issues introduced by re-scoping

### N1 — Phase 4 has no objective gate — **IMPORTANT**

EPIC §4 sub-plan table row 4: *"slots in when binseq stabilises and mim-index MSRV is compatible"*.

What does "stabilises" mean concretely? `binseq` could ship 1.0 next month, six months from now, or never. With 15 versions in 9 months and a sole maintainer (bus factor 1, per B-§1.1), the realistic outcome is **deferred indefinitely**. The "re-evaluate when ready" framing is fine for new work in flight, but for a deferred-and-justified-defer it creates a stale-issue problem: nobody owns "is binseq stable yet?" so nobody asks.

Two objective gates that would work:

- **Time-bounded re-evaluation.** "Re-evaluate Phase 4 disposition 6 months after Phase 3 ships. If binseq has reached 1.0 with at least 3 patch versions and no on-disk format changes since 1.0, light Phase 4; otherwise either close as won't-do or set a new 6-month checkpoint." This avoids zombie phases.
- **Specific version targets.** "Phase 4 lights when `binseq ≥ 1.0.0` (no breaking format revs since) AND `mim-index ≥ 0.2.0` with `rust-version ≤ 1.88`." Either condition unmet → stays deferred. Specific, checkable, no judgment call.

Right now Phase 4 reads like an indefinite parking lot. The bus-factor-on-binseq concern from v1 still applies — if the maintainer disappears, Phase 4 must be re-litigated from scratch (BINSEQ-the-format vs BINSEQ-the-crate). Worth a sentence in EPIC §3 Phase 4 paragraph or §8 Open questions making this explicit.

### N2 — Phase 2 Bismark inheritance is now under-defined — **IMPORTANT**

EPIC §3 Phase 2: *"Same `RecordSource` shape Bismark consumes — open question whether we extract the trait to a small `recordsource-rs` crate published from one repo and pulled by the other, or duplicate it minimally with a CI cross-check."*

EPIC §4 sub-plan table row 2 (Phase 2): *"Depends on #1 (dispatch pattern)."*

EPIC §7: *"Phase 2 is gated on Phase 1's trait stability."*

Three things to reconcile:

1. **`RecordSource` already exists** (verified: `src/fastq.rs:213`, `src/bam.rs:20`, used in `parallel.rs`, `trimmer.rs`, `format.rs`). It pre-dates this epic (came with #317). So Phase 2 inheriting `RecordSource` is well-defined and unchanged by v2.
2. **`RecordSink` was the new abstraction v1 proposed and v2 dropped.** Phase 2 doesn't need a sink trait — Bismark's output is its existing aligned SAM/BAM, which isn't being touched (EPIC §3 Phase 2: "alignment-output behaviour is unchanged").
3. **But EPIC §3/§7 still uses the singular phrase "trait shape" / "trait stability"** as if Phase 2's inheritance is undefined-until-Phase-1-decides. With v2, the decision IS made: trait = RecordSource only (already stable since #317), no sink trait.

So EPIC §3 Phase 2 and §7 should be updated to read something like: *"Phase 2 inherits the `RecordSource` trait shape, already stable since #317. The Phase-3 shared-crate-vs-mirrored question still applies."* Without this update, future readers of EPIC.md will be confused about what "dispatch pattern" Phase 2 inherits — they'll look for the RecordSink trait that was deferred.

**This is an EPIC-doc-drift issue, not a Phase-1 plan bug.** But it matters for "epic coherence" (per the review mandate).

### N3 — `@PG` on version bump → fixture regen on every TrimGalore release — **IMPORTANT**

§9 last row + §11 "Remaining risks" first bullet acknowledges: *"`@PG`-on-version-bump breaks fixture byte-identity. Mitigation: assert on (header.minus_@PG_lines, records) tuples in CI, not on the full byte sequence; document fixture-regen recipe."*

**Is "assert on tuples not bytes" implementable cleanly?** Yes — noodles exposes the SAM header as a structured object with separate sections (HD, PG, CO, SQ — though no SQ here since uBAM). A helper like:

```rust
fn strip_pg_lines(header: &sam::Header) -> sam::Header { … }
fn assert_ubam_eq(actual: &Path, expected: &Path) -> Result<()> {
    let (a_h, a_recs) = read_bam_header_and_records(actual)?;
    let (e_h, e_recs) = read_bam_header_and_records(expected)?;
    assert_eq!(strip_pg_lines(&a_h), strip_pg_lines(&e_h));
    assert_eq!(a_recs, e_recs);
}
```

…is straightforward. Worth being explicit in §5 step 5 about adding this helper — currently §5 step 5 just lists test names without the helper.

**CHANGELOG entry per release:** if the helper strips `@PG`, no — fixture stays valid across version bumps. If implementer accidentally locks bytes instead, every release breaks CI. **Recommend: §5 step 5 explicitly add a sub-bullet "Add `assert_ubam_eq` helper that strips @PG lines before comparison; commit the helper alongside the golden fixtures."** This pre-empts the most likely bug.

### N4 — uBAM input parallel / uBAM output serial asymmetry — **OPTIONAL**

This was raised in the mandate explicitly. Let me think it through:

- uBAM **input** today (post-#317): goes through `parallel.rs` worker pool via `Box<dyn RecordSource>` (verified: `src/parallel.rs:94-96, 646-648, 904, 1098`). So `trim_galore --cores 8 input.bam` parallelises trimming + gzip output.
- uBAM **output** in v2: serial-only with `--cores N>1` warning (§3.4 row 6, §3.5 second bullet).

User-facing surprise: *"I gave you `--cores 8`, parallel-read-parallel-trim works, but switching `--output-format ubam` quietly goes serial."*

The warning in §3.5 covers this: *"`--output-format ubam uses single-threaded compression in v1; --cores N is ignored. For high-throughput uBAM output, run multiple invocations in parallel.`"* That's a good UX warning. **But:** for pipelines that auto-tune `--cores` from the host (nf-core does this), the warning will fire on every run and just get logged-and-ignored. Not a bug, just noise.

**Marginal recommendation:** consider making the warning a one-shot stderr line at startup (per-invocation, not per-batch) phrased as *"uBAM output is single-threaded in v1; --cores 8 is ignored for compression. Trimming throughput is unaffected."* — that explicitly tells the user the parallel work that DOES happen (trim) vs the serial work (compress). Helps nf-core users grep CI logs without alarm.

Not a blocker.

### N5 — Is Phase 4 even worth keeping alive? — **OPTIONAL**

The mandate raises this directly. Let me reason about it:

EPIC §1 goal: *"TrimGalore → Bismark pipeline can run with any format pair without a FASTQ round-trip in between."*

Post-v2: TrimGalore can read uBAM (#317) and emit uBAM (Phase 1). Bismark Phase 2 reads uBAM. So the **uBAM-out → uBAM-in handshake satisfies the EPIC's goal entirely.** BINSEQ and mim are nice-to-haves layered on top.

**If Phase 4 never lights**, the cross-tool TrimGalore→Bismark pipeline still works. The EPIC's headline outcome is delivered by Phases 1+2+3 alone.

**This is a real shift from EPIC v1's framing**, where all four formats were the deliverable. v2's EPIC §3 Phase 4 paragraph treats it as a "deferred but coming" item, but functionally it's now **optional / aspirational**. The EPIC should explicitly acknowledge this — either:

- *"Phase 4 is now optional. Phases 1+2+3 deliver the headline cross-tool handshake. Phase 4 adds BINSEQ + mim as additional format pairs if/when their crate ecosystems mature."*
- Or close Phase 4 entirely and re-light as a new epic later.

Per the review mandate ("should the epic explicitly note that Phase 4 is now optional / aspirational?"): yes, this would be clarifying. Combined with N1 (objective gate), this is a clean disposition: time-box the re-evaluation, frame as optional, no zombie phase.

### N6 — `--preserve-tags` mixed-batch source-aux lookup state is "TBD at impl time" — **IMPORTANT**

§4 `BamWriter` doc comment:

> // For mixed-input batches: tag_source[i] = Some(BamReader) when input i was uBAM,
> // None when FASTQ. Used to look up aux fields per record.
> // (Detail TBD at impl time — may collapse to per-call lookup instead of stored state.)

This is the v2 loosening of A-O1 (per-pair tag preservation in mixed batches). The mechanism is genuinely unclear — `BamWriter` is per-output-file (one per pair), so a single `BamWriter` doesn't see "the batch"; it sees one pair's records. The mixed-batch case is handled at the **main.rs dispatch loop** level (one BamWriter per output, opened with the appropriate source_aux strategy).

But: how does `write_record` know whether the current record came from a uBAM input (carry aux) or FASTQ input (empty aux)? The §4 signature passes `source_aux: Option<&sam::record::Data>` per call, which means **the caller** (the trimmer entry point in `run_*_to_bam`) must know whether the current record's input was uBAM. That requires the trimmer to either:

- Take a second `Option<&BamReader>` per input position to source aux from, OR
- Receive aux as a side-channel on the FastqRecord (FastqRecord doesn't carry aux today).

Neither is in v2's §4 signatures. The "TBD at impl time" comment punts a real design question — not a showstopper, but **§5 step 3 should be explicit about which path is taken**. Otherwise the implementer has to invent it during step-3 coding.

Recommend: extend §5 step 3 with a sub-bullet specifying *"per-input tag-source plumbing: trimmer entry points take an additional `Option<&BamRecord>` per record OR FastqRecord gains an `aux: Option<sam::record::Data>` field. Recommend the latter — keeps the trimmer signature stable."* This is the same "make the design decision in the plan, not at the keyboard" rule the v1 reviewers applied to A-C1.

---

## 3 · Other observations

### `--cores` warning vs. clumpify rejection — UX inconsistency

`--clumpify` + uBAM → hard error (§3.4 row 1).
`--cores N>1` + uBAM → warn-and-continue (§3.5 second bullet).

Both are about parallelism-flag-meaningless-for-uBAM, but one rejects and the other warns. The asymmetry is fine — `--clumpify` is opt-in (user explicitly asked for it), while `--cores` is often set globally by pipelines (nf-core sets it from host CPU count, not by user intent). Different ergonomic profiles. Just worth noting that the **rationale** for the asymmetry should be in §3.4 / §3.5 in one sentence. Currently it reads as inconsistency.

### §3.3 step 1 "name strip" — encoding gotcha

§3.3 step 1: *"strip the leading `@` from `FastqRecord.id`"*. **Check that FastqRecord.id already excludes the `@`** — many parsers store the line minus the prefix. If so, the strip-`@` step is wrong (would chop the first character of the actual ID). Verifiable by reading `src/fastq.rs::FastqReader`. Not a plan bug per se, but a §5 step 2.4 unit test for `bam_writer_se_fastq_input` should explicitly assert the ID round-trips byte-for-byte from input FASTQ to output BAM `name`. (And v1 review's "no off-by-one on @" discipline applies.)

### §5 step 6 CI extension — `samtools` availability

Step 6 step 1 (CI new step): uses `samtools view -c -f 0x40` and `-f 0x80`. The validation job already installs Perl trim_galore via conda; samtools is presumably in the same conda env or PATH. **Worth verifying samtools is on the validation runner** before the plan lands — if it's not, this requires a conda channel addition, which is a small but separate CI change. One-line check / step in §6.

---

## 4 · Validation gap check (per skill mandate §4)

§9 covers the high-risk failure modes well:
- FASTQ byte-identity (Perl oracle) — unchanged from today, locks the no-regression invariant.
- BamWriter SE + PE round-trip — within-codebase sanity.
- Aux tag preservation, empty-aux on FASTQ-in, flag bits — explicit assertions.
- Five rejection rules — per-rule unit tests.
- `--cores N>1` + uBAM-out warning — integration test.
- Reproducibility — existing CI job runs unchanged.
- samtools content parity — external oracle.

**Missing:** `(header.minus_@PG, records)` helper assertion explicitly named (per N3). Recommend §9 add a row:

| `@PG`-tolerant fixture comparison | Helper `assert_ubam_eq` strips @PG lines before bytewise comparison | Pass across version bumps without fixture regen |

---

## 5 · Action items (prioritized)

### Important (address before implementation)

1. **N6** — §5 step 3: specify the per-input tag-source plumbing mechanism (FastqRecord gets an `aux` field, OR trimmer takes a parallel `Option<&BamRecord>` stream). Don't leave "TBD at impl time."
2. **N2** — EPIC §3 Phase 2 + §7 wording update: Phase 2 inherits `RecordSource` (already stable since #317), NOT a future `RecordSink`. Reword "trait stability" to clarify which trait. This is a doc-coherence fix in EPIC.md, not in PLAN.md.
3. **N1** — EPIC §3 Phase 4 paragraph: add an objective re-evaluation gate ("re-evaluate 6 months post-Phase-3 ship" OR specific version targets like `binseq ≥ 1.0` + `mim-index ≥ 0.2 with MSRV ≤ 1.88`). Avoid the zombie-phase outcome.
4. **N3** — §5 step 5 + §9: explicitly add the `assert_ubam_eq` helper that strips `@PG` lines, so fixture byte-identity survives TrimGalore version bumps without regen-on-every-release.
5. **A-O3 carry-over** — §3.4: add `--demux` + `--output-format ubam` to the enumerated rejection table (currently only in §7 prose). The table is sold as canonical.

### Optional

6. **N5** — EPIC: explicitly acknowledge Phase 4 is now optional (Phases 1+2+3 deliver the headline goal). Either re-frame Phase 4 as "additional format pairs, optional" or close it and re-light later as a new epic.
7. **A-O6 carry-over** — §3.5: one explicit sentence on per-pair header synthesis ("each output BAM gets its own header derived from its own input pair's source").
8. **N4** — soften `--cores N>1` + uBAM-out warning wording so it's pipeline-friendly (one-shot at startup, explicit about which work IS parallel).
9. **§3.3 step 1** — verify `FastqRecord.id` storage convention (with or without leading `@`); add explicit ID byte-identity assertion in `bam_writer_se_fastq_input`.
10. **§5 step 6** — verify samtools is on the validation CI runner; add to env if missing.

### Nothing in "Critical" bucket — the re-scoping itself is sound

---

## 6 · Bottom line

**v2 is a substantial improvement over v1.** Every CRITICAL is genuinely fixed (not just claimed-fixed), and the re-scoping rationale (§0 + §11) is honest about the trade-offs. Dropping the RecordSink trait is the right call given parallel.rs's actual architecture, and the BINSEQ/mim deferral is justified by concrete crate-maturity findings.

The remaining items are:
- One real design ambiguity (N6 — per-input tag-source plumbing).
- Two EPIC-coherence issues (N1 Phase-4 gate, N2 Phase-2 inheritance wording) — neither blocks Phase 1 impl but both should be addressed in the EPIC concurrently with Phase 1 landing.
- One CI-resilience hardening (N3 — `@PG`-stripping helper named explicitly).
- A few minor table/prose consistency cleanups.

None of these warrant another v3 rev. They can land as a v2.1 cleanup or, for the EPIC items, in EPIC.md directly.

**Verdict: PARTIAL — ready for the implementation trigger after N6 + N3 are addressed in PLAN.md, and N1 + N2 in EPIC.md.**

---

## 7 · Reviewer D scope notes (just so we don't double-bill with Reviewer C)

Things this reviewer did NOT closely examine (deferred to Reviewer C's different lens):
- Detailed implementation correctness of §3.3 per-record conversion (Phred subtraction, IUPAC rejection, sequence-encoding edge cases). Spot-checked but not exhaustive.
- noodles 0.88 API surface for `bam::io::Writer` + header propagation specifics.
- Detailed test enumeration in §5 step 5 vs. the existing `integration_passthrough.rs` / `integration_ubam_in.rs` patterns.
- The `--passthrough` rejection's downstream consequences for nf-core/methylseq pipelines.
- Specific samtools-version requirements for `view -f 0x40` flag-bit semantics.
- File extension `.bam` collision risk with existing fixtures (`test_files/ubam_*.bam`).

Reviewer D focused on: (1) v2 → v1 finding traceability, and (2) NEW issues introduced by re-scoping.
