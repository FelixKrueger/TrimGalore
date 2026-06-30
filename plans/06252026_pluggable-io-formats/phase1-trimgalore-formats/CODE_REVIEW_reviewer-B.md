# Code Review — Phase 1 uBAM Output (`--output-format ubam`) — Reviewer B

**Scope:** working-tree diff (`git diff HEAD` + untracked new files) implementing Steps 2–7 of PLAN v2.1.
**Spec:** `plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN.md` (v2.1).
**Method:** read every changed hunk; cross-referenced against PLAN §3.3/§3.4/§3.6/§9; re-ran build, clippy, full test suite; wrote and ran throwaway probes for sequence-encoding, qual-mismatch, and `@PG`-chain behaviour (all probes deleted).

## Summary

The implementation is well-structured, faithful to the v2 architectural decision (no trait; FASTQ path untouched; serial uBAM branch), and the `write_paired_reports` refactor cleanly removes the duplication PLAN §5 step 0 called out. Verification gate reproduced locally: `cargo build --release` clean, `cargo clippy --release --all-targets -- -D warnings` clean, full suite **336 tests pass** (317 lib + 0 + 2 + 8 + 9 + 0). The new `BamWriter`, the per-record conversion, the flag bits (0x04 / 0x4D / 0x8D), the tag round-trip parser, and the rejection rules are all correct on the happy path and exercised by good tests.

However, there are **two real correctness gaps** (both IMPORTANT, both verified empirically) and several smaller spec/coverage gaps. None are CRITICAL — read data is never corrupted on the common path — but the `@PG`-collapse and the missing two-BAM-pair guard each produce silently-wrong-shaped output for specific (legitimate) inputs.

**Verdict: APPROVE_WITH_FIXES.**

---

## IMPORTANT

### I1 — Re-processing a Trim-Galore uBAM silently collapses the `@PG` provenance chain
`src/bam.rs:174-188` (`build_output_header`) + `src/bam.rs:211-218` (`last_pg_id`).

`build_output_header` propagates the source header by **cloning it and calling the raw `IndexMap` insert**:
```rust
h.programs_mut().as_mut().insert("trim_galore".into(), pg);
```
It bypasses noodles' purpose-built `Programs::add()` (noodles-sam 0.69.0 `header/programs.rs:71`), which exists precisely to (a) chain `PP:` to the current leaf and (b) disambiguate a duplicate ID by appending `-<prev-id>`.

**Verified:** I wrote a probe that runs a FASTQ→uBAM pass (producing `@PG ID:trim_galore`), then feeds that BAM's header as the source for a second pass. Result:
```
second.bam @PG count=1   keys=["trim_galore"]
```
The first pass's `@PG` is **overwritten**, not chained. noodles' `add()` would have yielded `["trim_galore", "trim_galore-trim_galore"]` with a proper `PP:` link. So uBAM-in → uBAM-out where the input was *itself* a Trim Galore product loses prior provenance — directly contradicting PLAN §3.5 ("provenance is preserved by adding to history, not preserving silence").

Separately, `last_pg_id` uses `programs().as_ref().keys().last()` (raw insertion-order last key) rather than the chain **leaf**. For the common single-chain uBAM header (one `@PG` from `samtools import`/Picard/fgbio) these coincide, so `PP:` is correct in practice; but for a branched/multi-chain source header the `PP:` target may be wrong.

**Fix:** use `header.programs_mut().add("trim_galore", pg_without_pp)` and let noodles compute `PP:` + ID-disambiguation, instead of hand-rolling `last_pg_id` + raw `insert`. (`add` returns `io::Result`; surface the error.) This also makes the `last_pg_id` helper unnecessary.

### I2 — Two BAM files in `--paired --output-format ubam` are silently accepted; the FASTQ path rejects the same input
`src/main.rs:1587` (`run_ubam_output_paired_two_files`) vs `src/main.rs:911-921` (`run_paired`).

The FASTQ-output paired path explicitly rejects two-BAM-file pairs:
```rust
"--paired with two BAM files is not supported. uBAM paired mode expects a single interleaved file…"
```
The uBAM-output paired path has **no equivalent guard**. `run_ubam_output` (main.rs:1392) only special-cases `cli.input.len() == 1`; with two BAM inputs it falls through to `run_ubam_output_paired_two_files`, which opens each BAM as an independent R1/R2 stream.

**Verified:**
```
trim_galore --paired --output-format ubam ubam_test.bam ubam_test_with_tags.bam   →  EXIT=0
```
It "works", but it treats two distinct single-end uBAMs as a synthetic pair — exactly the input the architecture says is unsupported (uBAM paired is one interleaved file, matching #317's reader). Behaviour now diverges by output format for identical input, which is surprising and undocumented. (The same run also exposes the cosmetic report-name `.bam` leak — see N2.)

**Fix:** add the same two-BAM rejection at the top of `run_ubam_output_paired_two_files` (or, better, in `run_ubam_output` before dispatch, mirroring main.rs:911-921).

---

## MINOR / spec gaps (between IMPORTANT and NIT)

### M1 — Write path does not enforce the ACGTN-only / reject-`=`-and-IUPAC rule that PLAN §3.3 step 3 specifies
`src/bam.rs:149-155` (`write_record`). The sequence is handed to `Sequence::from(record.seq.as_bytes().to_vec())` with no validation. noodles' 4-bit encoder accepts all 16 IUPAC slots.

**Verified** (probe, deleted):
- input `acgtRYN` → output BAM seq `ACGTRYN` (lowercase silently upper-cased — benign/good; **IUPAC `R`/`Y` pass through unchanged** — not coerced to `N`).
- input `AC=T` → output BAM seq `AC=T` (the `=` match-reference sentinel survives).

The **read** path (`bam_record_to_fastq`, src/bam.rs:852-867) is strict: IUPAC→N, `=`→hard error. So the two directions are asymmetric. This only bites for **FASTQ input containing non-ACGTN bases** (legal in FASTQ, e.g. soft-masked lowercase or IUPAC ambiguity codes) — `FastqReader::sanity_check` only screens for colorspace, not IUPAC. For uBAM input the read side already normalised, so this is a no-op there.

Impact is bounded (uncommon input), but it is a documented spec requirement that isn't met, and a `=` in unaligned BAM is meaningless. Recommend either enforcing §3.3 step 3 in `write_record` (cheap: map non-ACGTN→N, bail on `=`) or amending the PLAN to record the relaxation. Note `Sequence::from`'s implicit upper-casing is actually a desirable normalisation worth keeping.

### M2 — Two PLAN §9 validation rows have no test
- **`--cores N>1` + uBAM-out warns + still works** (§9, §3.5). The warning path (main.rs:1369) is correct by inspection and runs serially, but no test asserts exit 0 + "ignored" on stderr for `--cores 4 --output-format ubam`. The only `--cores 2` invocation (`ubam_out_clumpify_rejected_at_cli`) errors out *before* reaching the warning.
- **`--preserve-tags` mixed-batch** (§9 + §3.4b note, which names `tests/integration_ubam_out.rs::preserve_tags_mixed_batch_allowed` explicitly). No such test exists. The A-O1 loosening (one FASTQ + one uBAM input → FASTQ-side empty aux, uBAM-side tags preserved) is therefore unverified end-to-end.

Both are coverage gaps, not code defects.

### M3 — Missing-qual stderr NOTE (PLAN §3.3 step 4) not emitted
`grep` for "missing-qual"/"Phred 0"/"0xFF sentinel" finds nothing in main.rs/bam.rs/trimmer.rs. PLAN §3.3 step 4 specified a startup stderr note when uBAM-in → uBAM-out is detected, warning that the missing-qual sentinel becomes Phred 0. The **behaviour** (Phred-0 emission via `saturating_sub(33)`) is correct and documented in CHANGELOG/CLAUDE.md; only the runtime heads-up note is absent. Low impact; add the note or drop the requirement from the PLAN.

---

## NIT

### N1 — `_preserve_tags` field is stored but never read
`src/bam.rs:76`. Intentional (documented as a symmetry hint), prefixed `_` to silence the warning. Fine as-is, but it's dead state; could be dropped. Non-blocking.

### N2 — Report filename carries the `.bam` extension for BAM input
`ubam_test.bam_trimming_report.txt`. Caused by `io::report_name` using the full `file_name()` (src/io.rs:255-261). **Pre-existing** behaviour shared with the FASTQ-output BAM-input path (not introduced by this diff), so consistent — but cosmetically odd. Out of scope for this review; flagging only because it surfaced in the I2 repro.

### N3 — `debug_assert_eq!` on qual/seq length (src/bam.rs:144-148) is compiled out in release
Harmless: I verified that a deliberate length mismatch makes noodles' `write_alignment_record` return an `Err` (not panic, not silent) in `--release`, so the invariant is still enforced at the boundary. The comment is accurate. No action.

---

## Things checked and found correct

- **Flag bits** 0x04 / 0x4D / 0x8D — correct, asserted at the raw-byte layer (`bam_writer_flag_bits_correct`).
- **Qual conversion** `saturating_sub(33)` — correct inverse of read-side `+33`; defensive against underflow.
- **Tag round-trip parser** (`parse_name_and_data` / `parse_tag_value`) — TAB-vs-SPACE split correct (the bug-fix in deviation #4 is sound and regression-guarded); A/Z/i/f accepted, B/H/unknown rejected; empty-name, missing-type, trailing-tab edge cases all covered by lib tests. Leading-space and bare-`@` ids correctly bail.
- **`open_inner` → `(Reader, Header)` refactor** — all three call sites (open, open_threaded closure, open_paired_interleaved closure) destructure as `Ok((r, _header))`; no information lost; threaded closures unchanged in behaviour.
- **Header synthesis for FASTQ input** — `@HD VN:1.6` + trim_galore `@PG`, verified present (`bam_writer_synthesises_header_for_fastq_input`).
- **Typed aux fidelity** — `i`→Int32, `f`→Float land as real BAM numeric types, not Z-strings (`bam_writer_aux_typed_int_and_float_round_trip`).
- **Trimmer entry points** (`run_single_end_to_bam` / `run_paired_end_to_bam`) — stats accumulation, filter handling, discard-untrimmed, truncation-error arms all match the FASTQ counterparts line-for-line, minus the correctly-rejected passthrough/unpaired branches.
- **CLI rejection rules** (§3.4a: clumpify/passthrough/clock/implicon/demux/retain_unpaired) — all present, all with rejection tests; the §3.4b format-detection-time rule wired in main.rs:148-165 with an integration test.
- **`--retain_unpaired` rejection** (deviation #2) — sound; matches the multi-output-BAM-out-of-scope pattern.
- **hardtrim5/3 → BAM naming** (deviation #3) — `<stem>.<N>bp_{5,3}prime.bam` preserves the 5'/3' discriminator; CWD-when-no-output_dir behaviour matches the FASTQ `hardtrim_output_name` exactly (both `PathBuf::from(filename)`), so no divergence.
- **CI** (`validation-ubam` new steps) — SE content-parity (`samtools fastq -n` + `awk 'NR%4!=3'`), PE FREAD1==FREAD2 count + flagstat, aux-tag grep — all correct; `-n` correctly strips `/1` suffixes for name parity.
- **`BamWriter::create`'s 4th param `command_line`** (deviation #1) — reasonable; keeps header construction self-contained.

---

## Recommended action order
1. **I1** — switch `build_output_header` to `Programs::add()` (fixes provenance-collapse + leaf-PP correctness; removes `last_pg_id`).
2. **I2** — add the two-BAM-pair guard to the uBAM-output path.
3. **M1** — enforce §3.3 step 3 in `write_record` (or amend PLAN).
4. **M2 / M3** — add the two missing §9 tests; add or drop the missing-qual NOTE.

Verdict: **APPROVE_WITH_FIXES** — ship after I1 and I2; M1–M3 can follow but should not be left silently open against the PLAN.
