# Plan Coverage Report — v2 (uBAM in/out) audit

**Mode:** B (code vs. design plan `PLAN_v2_ubam.md`)
**Plan:** `plans/clump-only/PLAN_v2_ubam.md` (551 lines)
**Progress log:** `plans/clump-only/PROGRESS.md` §"v2 Implementation log"
**Branch:** `feature/clump-only-ubam` (v2 diff = uncommitted working-tree edits vs. `feature/clump-only`)
**Date:** 2026-07-25
**Verdict:** **INCOMPLETE WITH MINOR GAPS** — 5 test-coverage gaps against the plan; core implementation is complete.

## Summary

- Total items audited: 34 (5 resolved decisions + 10 impl-outline steps + 9 validation items + 4 rejection-matrix items + 3 contract invariants + 6 self-review adjustments — with three §Resolved-decisions duplicating checks elsewhere collapsed to a single row each).
- DONE: 29
- PARTIAL: 4 (§Impl-outline step 7, §Impl-outline step 8, §Validation #6 preserve-tags round-trip, §Validation #7 shape dispatch)
- MISSING: 1 (§Validation #4 cross-run record determinism)
- DEVIATED (documented): 0 material (only the `PairedInputSetup` struct — noted in PROGRESS.md, same intent as plan)

## Coverage ledger

### §Resolved decisions (5)

| # | Item | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | `--fastqc` on uBAM → run FastQC directly on BAM | `src/clump_only.rs` `clump_only_single_to_bam` + `clump_only_paired_to_bam_one_pair` call `fastqc::run(&output_path, …)`; integration test `fastqc_produces_report_on_ubam_out` | DONE |
| 2 | Compression-ratio: emit only when both sides compressed | `src/clump_only.rs::write_clump_only_report` gates on `stats.input_compressed && stats.output_compressed && stats.output_bytes > 0`; `input_is_compressed` treats gzip + BGZF as compressed | DONE |
| 3 | `estimated_record_bytes` audit (verified during plan-review) | No code change (per plan) | DONE |
| 4 | Empty-input BAM handling (verified during plan-review) | No code change; `BamWriter::finish` handles zero records | DONE |
| 5 | `@PG` chain unbounded growth accepted | `BamWriter::create(&output, source_header.as_ref(), …)` propagates `source` chain then appends our `@PG` via `build_output_header` (existing) | DONE |

### §Implementation outline (10 steps)

| # | Step | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | CLI edits — 4 precise sites | (a) `reject_ubam` and its two call sites removed from `src/clump_only.rs` (replaced with narrower "add --output-format ubam" hint on FASTQ path); (b) `--clump_only + --output-format ubam` v1-defer bail deleted from `src/cli.rs`; (c) N=1+`--paired` guard removed from `src/cli.rs` clump-only block (dispatch-time format check in main.rs replaces it); (d) `--dont_gzip + --output-format ubam` bail added to shared §3.4a `OutputFormat::UBam` block in `src/cli.rs:544-554` | DONE |
| 2 | Main dispatch update | `src/main.rs` `if cli.clump_only { … match cli.output_format { UBam => … } }` — Shape B / Shape A / SE branches with format-guards; `preflight_collision_bam` helper hashes case-folded paths | DONE |
| 3 | Core BAM implementations | `src/clump_only.rs::clump_only_single_to_bam` + `clump_only_paired_to_bam_one_pair` + `flush_bin_single_to_bam` + `flush_bin_paired_to_bam` + `PairedInputSetup` — mirror v1 shape with `BamWriter::write_record(&rec, None / Some(1) / Some(2))` | DONE |
| 4 | Filename helpers | `src/io.rs::clumped_bam_output_name` + `clumped_paired_bam_output_name` + 7 unit tests | DONE |
| 5 | Report writer extension | `ClumpOnlyStats` extended with `input_format_label` / `output_format_label` / `preserved_tags`; `write_clump_only_report` renders new fields; compression-ratio gated | DONE |
| 6 | `estimated_record_bytes` audit | Verified during plan-review — no code change (per plan §Resolved 3) | DONE |
| 7 | Tests (unit + integration) | Integration: 9 tests present in `tests/integration_clump_only_ubam.rs`. Plan asked for 11 integration + 4 dedicated BAM unit tests | **PARTIAL** — see Gaps |
| 8 | CI validation steps | `.github/workflows/ci.yml` — SE uBAM record-parity + PE uBAM mate-adjacent record-parity + `@PG` chain check. Plan asked for 4 CI steps | **PARTIAL** — see Gaps |
| 9 | Documentation | `docs/src/content/docs/modes/clump-only.md` rewritten with dual FASTQ/uBAM usage, aux-tag round-trip, `@PG` line note, PE shape acceptance/rejection matrix, updated report shape | DONE |
| 10 | CHANGELOG.md | New `### Unreleased` bullet describing v2 + `--dont_gzip` rejection hoist | DONE |

### §Validation (9 items)

| # | Item | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | Record byte-identity through BAM round-trip | `se_ubam_in_ubam_out_pg_chain` (integration) + `pe_ubam_in_interleaved_output` (integration) both `assert_eq!(sorted_bam_tuples(&input), sorted_bam_tuples(&out))`; CI SE + PE record-parity | DONE |
| 2 | Pair lockstep in interleaved output | `pe_ubam_out_interleaved_from_fastq_pair` asserts flag 0x4D/0x8D pattern and equal names at every even i | DONE |
| 3 | `@PG` chain preservation | `se_ubam_in_ubam_out_pg_chain` (asserts `output_pg_count == input_pg_count + 1` and `has_trim_galore_pg`); CI `@PG` step | DONE |
| 4 | Cross-run record determinism (ignoring `@PG`) | Not implemented as a v2-specific CI step or unit test | **MISSING** |
| 5 | Rejection: `--dont_gzip + --output-format ubam` | `src/cli.rs:544-554` bail; integration test `rejects_dont_gzip_with_ubam_output` | DONE |
| 6 | `--preserve-tags` round-trip (RG:Z + BC:Z + NM:i + AS:f) | No dedicated multi-type test with synthetic tag values as plan specified. `se_ubam_in_ubam_out_pg_chain` + `pe_ubam_in_interleaved_output` incidentally round-trip whatever tags exist on fixtures via `sorted_bam_tuples`, but no explicit A/Z/i/f cross-type check | **PARTIAL** |
| 7 | PE input-shape dispatch (Shape A + Shape B both produce interleaved output) | `pe_ubam_out_interleaved_from_fastq_pair` (Shape A) + `pe_ubam_in_interleaved_output` (Shape B) both verify output shape | DONE |
| 8 | `--fastqc` produces `*_fastqc.html` + `.zip` on BAM | `fastqc_produces_report_on_ubam_out` | DONE |
| 9 | Empty-input BAM handling | Verified during plan-review — no test kept for the clump-only BAM path (per user's audit direction, NOT a gap) | DONE |

### §Rejection matrix updates (4 items)

| # | Item | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | `--dont_gzip + --output-format ubam` (hoisted to shared §3.4a) | `src/cli.rs:544-554` + integration test | DONE |
| 2 | Shape A two-BAM | `src/main.rs` dispatch loop (`r1_is_bam && r2_is_bam` bail); integration `rejects_two_bam_paired` (note: fixture uses same file for R1+R2 which is caught first by identical-file guard — test still asserts exit != 0 per plan) | DONE |
| 3 | Shape B non-BAM (N=1 FASTQ under `--paired --output-format ubam`) | `src/main.rs` dispatch (`if cli.input.len() == 1 && !UnalignedBam { bail! }`); integration `rejects_non_bam_n1_paired` | DONE |
| 4 | Mixed-format Shape A (FASTQ + BAM in same pair) | `src/main.rs` dispatch (`r1_is_bam != r2_is_bam { bail! }`); NO integration test | **PARTIAL** — code present, test missing |

### §Behavior/Contract invariants (3)

| # | Item | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | Record byte-identity (id + seq + qual + preserved tags) | Integration `sorted_bam_tuples` equality checks; CI record-parity | DONE |
| 2 | `@PG` chain preservation (input chain + new line) | `BamWriter::create` propagates via `build_output_header`; test asserts `output_pg_count == input_pg_count + 1` | DONE |
| 3 | Mate-adjacent PE-BAM output | `flush_bin_paired_to_bam` alternates `Some(1)`/`Some(2)`; test asserts flag/name pattern | DONE |

### §Self-Review Adjustments made during post-plan-review revision (6)

| # | Adjustment | Where implemented | Status |
|---|------|-------------------|--------|
| 1 | Purge "warn-and-skip" FastQC position everywhere | §Rejection matrix + §Signature + §Impl + §Validation all agree; `fastqc::run` invoked on BAM path | DONE |
| 2 | Multi-pair PE support via `_one_pair` signature | Function renamed; `main.rs` iterates `cli.input.chunks(2)`; integration `multi_pair_pe_bam_produces_one_output_per_pair` | DONE |
| 3 | Explicit format-guard rejections (two-BAM Shape A / non-BAM Shape B / mixed Shape A) | All 3 code guards in `main.rs`; 2 of 3 integration tests present | PARTIAL (Rejection #4 test missing) |
| 4 | Explicit collision pre-flight for PE-BAM + multi-SE-BAM | `preflight_collision_bam` helper called on all 3 BAM branches; no integration test proves case-folded matching | PARTIAL (test missing) |
| 5 | `--dont_gzip + --output-format ubam` hoisted to shared §3.4a | `src/cli.rs:544-554`; applies to trim uBAM path too | DONE |
| 6 | Precise CLI edit sites | All 4 sites reflected in diff (see §Impl outline step 1) | DONE |

## Gaps (detail)

### Gap 1 — §Impl outline step 7: 4 dedicated BAM unit tests missing (PARTIAL)

**Expected:** `test_clump_only_single_to_bam_permutation`, `test_clump_only_paired_to_bam_lockstep`, `test_clump_only_bam_deterministic_records`, `test_clump_only_ubam_in_ubam_out_tag_roundtrip` (unit tests inside `src/clump_only.rs`).
**Found:** No new unit tests were added to `src/clump_only.rs` for the BAM path; the 11 v1 unit tests remain. Coverage is delegated entirely to integration tests.
**Gap:** Unit-level coverage of `flush_bin_single_to_bam` / `flush_bin_paired_to_bam` / determinism at the function boundary is absent; if the binary shell breaks (e.g. arg-parsing regresses) the integration tests fail non-diagnostically. The tag-roundtrip case in particular has no direct verification.
**Remediation:** Add the four unit tests as specified in the plan, using in-memory `Cursor` writers and `noodles::bam::io::Reader` for the decode side. Lowest-value-add is `deterministic_records` (implicitly covered by integration path); highest is `ubam_in_ubam_out_tag_roundtrip` (see Gap 3).

### Gap 2 — §Impl outline step 7 + §Rejection matrix #4: two integration tests missing (PARTIAL)

**Expected:** `rejects_mixed_format_paired` (Shape A with `R1.fq.gz` + `R2.bam`) and `pe_bam_collision_preflight_case_folded` (two input pairs whose case-folded output paths collide).
**Found:** Neither test exists in `tests/integration_clump_only_ubam.rs`. The underlying code guards (`r1_is_bam != r2_is_bam` bail in `main.rs`; `preflight_collision_bam` helper) are present.
**Remediation:** Add both integration tests. `rejects_mixed_format_paired` should assert exit != 0 with a message containing "same format" or "mixed". `pe_bam_collision_preflight_case_folded` should provide two input-pair combinations whose case-folded output stems collide (e.g. on non-APFS runners simulate by passing `A_R1.fq` twice with `--basename A`) and assert exit != 0 with the "output path collision" message.

### Gap 3 — §Validation item 6: preserve-tags round-trip has no dedicated multi-type test (PARTIAL)

**Expected:** Synthetic uBAM with `RG:Z:sample`, `BC:Z:AAAA`, `NM:i:0`, `AS:f:1.5`; run with `--preserve-tags RG,BC,NM,AS`; assert all four tags present on every output record with same values. Kills tag-ordering swap, type-code mistranslation, value-string corruption.
**Found:** No test explicitly constructs a uBAM with A/Z/i/f tags of the values specified. Existing integration tests (`se_ubam_in_ubam_out_pg_chain`, `pe_ubam_in_interleaved_output`) use whatever tags happen to be on `test_files/ubam_test.bam` / `ubam_paired_test.bam` fixtures.
**Remediation:** Add a unit test in `src/clump_only.rs` that constructs an in-memory uBAM with the four tag types, runs `clump_only_single_to_bam` with `--preserve-tags`, decodes output, and asserts type + value preservation on every record. This is Validation Item 6's stated purpose: "Kills: tag ordering swap, type-code mistranslation, value-string corruption."

### Gap 4 — §Impl outline step 8 + §Validation item 4: cross-run BAM determinism CI step missing (MISSING)

**Expected:** CI step running the same `--clump_only --output-format ubam` invocation twice, `samtools view <out>` piped to `sort` on both, `md5sum` compare on the record-body portion. Kills: unstable sort, non-deterministic BAM writer ordering.
**Found:** No such CI step. The three added CI steps (SE record-parity, PE record-parity, `@PG` chain) do not verify record-body determinism across two runs.
**Remediation:** Add a fourth CI step under `validation` job:

```yaml
- name: Validate --clump_only uBAM cross-run record determinism (ignoring @PG)
  run: |
    rm -rf /tmp/co-bam-det; mkdir -p /tmp/co-bam-det/a /tmp/co-bam-det/b
    ./target/release/trim_galore --clump_only --output-format ubam --cores 2 \
      -o /tmp/co-bam-det/a test_files/ubam_test.bam
    ./target/release/trim_galore --clump_only --output-format ubam --cores 2 \
      -o /tmp/co-bam-det/b test_files/ubam_test.bam
    samtools view /tmp/co-bam-det/a/ubam_test_clumped.bam | sort | md5sum > /tmp/co-bam-det/a.md5
    samtools view /tmp/co-bam-det/b/ubam_test_clumped.bam | sort | md5sum > /tmp/co-bam-det/b.md5
    diff -q /tmp/co-bam-det/a.md5 /tmp/co-bam-det/b.md5 \
      || { echo "MISMATCH: --clump_only uBAM cross-run determinism failed"; exit 1; }
    echo "--clump_only uBAM cross-run determinism (record body): OK"
```

Or attach it as a unit-level determinism test in `src/clump_only.rs` (`test_clump_only_bam_deterministic_records` from the plan).

## Verified deviations (documented in PROGRESS.md)

- **`PairedInputSetup` struct extraction** — same intent as the plan's `Signature` block; introduced to sidestep `clippy::type_complexity` on the underlying 6-tuple. No contradiction with the plan.

## Verdict rationale

Core implementation (all 10 impl-outline steps' code) is present and matches the plan's shape precisely — the 4 CLI edit sites, the multi-pair dispatch, the `PairedInputSetup`-wrapped BAM function pair, the filename helpers, the extended stats + report writer, the docs, and the CHANGELOG all landed. All 3 §Behavior/Contract invariants (record byte-identity, `@PG` chain, mate-adjacent PE) are exercised by tests. All 5 §Resolved decisions and all 6 §Self-Review adjustments are honored in code.

The gaps are all test-coverage gaps against explicit plan items:
- 4 dedicated BAM unit tests (Gap 1) — none present; delegated to integration
- 2 integration tests (Gap 2) — `rejects_mixed_format_paired` + `pe_bam_collision_preflight_case_folded`
- 1 preserve-tags cross-type round-trip test (Gap 3) — no dedicated multi-type assertion
- 1 CI cross-run determinism step (Gap 4) — no such step

None of these gaps affect correctness of the shipped code; they affect regression-guard breadth. Merging as-is is defensible if the caller accepts the reduced regression surface; otherwise, adding the 4 gaps (roughly 200 lines of test code + one CI step) closes them without touching production code.
