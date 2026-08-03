# Plan Coverage Report — `--clump_only`

**Mode:** B (code vs. design plan) — no separate impl-plan file; PLAN.md is both
**Plan:** `plans/clump-only/PLAN.md`
**Branch:** `feature/clump-only` (uncommitted working tree off `dev`)
**Verdict:** **INCOMPLETE WITH MINOR GAPS** — 2 missing tests (⚠️), 1 minor doc/impl inconsistency (⚠️). No load-bearing invariant broken.

## Summary

- Total items audited: **34** (12 Resolved decisions + 11 Impl steps + 9 Validation items + 2 auxiliary — rejection matrix rolled up, contract invariants rolled up)
- ✅ Fully covered: **31**
- ⚠️ Partial / minor gap: **3**
- ❌ Missing: **0**

## Coverage matrix

### §Resolved decisions (12)

| # | Decision | Where implemented | Status |
|---|---|---|---|
| 1 | `--dont_gzip` allowed under `--clump_only` | `cli.rs:685-822` never rejects `--dont_gzip`; `main.rs:326-387` threads `gzip`; `io.rs:265-289`; integration `se_byte_identity_dont_gzip` | ✅ |
| 2 | Stable sort only in v1 | `clump_only.rs:44-46` reuses `sort_single_by_key` / `sort_paired_by_key` (already-stable per PROGRESS) | ✅ |
| 3 | Report filename `*_clumping_report.txt`, no JSON | `io.rs:326-333`; `clump_only.rs:554-601` writes text only | ✅ |
| 4 | `--fastqc` opt-in, runs on reordered output | `clump_only.rs:362-364, 526-529` | ✅ |
| 5 | PE report mirrors `--clumpify`'s per-input layout | `clump_only.rs:511-524` writes one report per mate | ✅ |
| 6 | Byte-identity scoped to header+seq+qual; plus-line + CRLF normalization documented | `cli.rs:200-204` (--help), `docs/.../clump-only.md:30-43` | ⚠️ *see Gap 3* |
| 7 | CI uses `paste - - - -` before sort | `.github/workflows/ci.yml:638-668` | ✅ |
| 8 | Rejection matrix: `--rename`, `--discard_untrimmed`, `--consider_already_trimmed`, uBAM input rejected explicitly | `cli.rs:772-786`; `clump_only.rs:233-243` (uBAM at `detect_input_format` time, not just CLI) | ✅ |
| 9 | `-q` / `--stringency` / `-e` silently accepted | `cli.rs:685-822` (never rejected); integration `silently_accepts_quality_flag` | ✅ |
| 10 | Sort determinism already in place — verification-only | Documented in PROGRESS.md; unit `test_clump_only_deterministic` + CI `Validate --clump_only cross-run determinism` | ✅ |
| 11 | `--basename` compatibility (`<basename>_clumped.fq(.gz)`) | `io.rs:265-311` accepts `basename`; `main.rs:331, 342, 377` threads it; io.rs unit tests `test_clumped_output_name_with_basename` + `test_clumped_paired_output_names_with_basename` | ✅ |
| 12 | Report "Compression ratio" line omitted when compression states differ | `clump_only.rs:595-598` (only emitted when both `input_compressed && output_compressed`) | ✅ |

### §Implementation outline (11 steps)

| # | Step | Where | Status |
|---|---|---|---|
| 1 | CLI: `clump_only: bool` field + `Cli::validate` rejection matrix; `--cores >= 2`; `--memory` parse | `cli.rs:207-208`, `cli.rs:685-822` | ⚠️ *see Gap 1 — `--cores >= 2` deliberately relaxed* |
| 2 | Main dispatch: early-return block, no adapter auto-detect | `main.rs:326-388` | ✅ |
| 3 | Core impl in new `src/clump_only.rs` (single + paired + uBAM reject + fastqc wiring) | `src/clump_only.rs` (1034 lines) + `lib.rs` `pub mod clump_only` | ✅ |
| 4 | Sort determinism verify-only | Verification confirmed in PROGRESS.md §Documented deviations | ✅ |
| 5 | Report: `ClumpOnlyStats` + text writer | `clump_only.rs:59-69, 554-601` | ✅ |
| 6 | Filename helpers in `src/io.rs` | `io.rs:265-333` (`clumped_output_name`, `clumped_paired_output_names`, `clumping_report_name`) | ✅ |
| 7 | Output-collision pre-flight for `--clump_only` | `main.rs:354-372` (SE case-folded hash); PE via `run_specialty_paired` | ✅ |
| 8 | Unit + integration tests (10 unit tests listed, 6 integration scenarios) | 9 unit tests in `clump_only.rs:719-1027` + 9 integration tests in `tests/integration_clump_only.rs` | ⚠️ *see Gap 2* |
| 9 | CI validate job: SE byte-identity + PE byte-identity + determinism | `.github/workflows/ci.yml:638-691` (3 steps + report-filename check inline) | ✅ |
| 10 | Docs: new Astro page + sidebar + cross-link | `docs/src/content/docs/modes/clump-only.md`; `docs/astro.config.mjs:78`; `docs/.../performance/clumpy.md:10-11` | ✅ |
| 11 | CHANGELOG entry under `### Unreleased` | `CHANGELOG.md:8-19` | ✅ |

### §Validation (9 items)

| # | Validation | Where | Status |
|---|---|---|---|
| 1 | Record byte-identity | Unit `test_clump_only_single_permutation` + integration `se_byte_identity_gzip` + CI SE step | ✅ |
| 2 | Multiset preservation (no dedup/drops/adds) | Same test bank; asserts on `outputs.len() == inputs.len()` and sorted-multiset diff | ✅ |
| 3 | Pair lockstep | Unit `test_clump_only_paired_lockstep` + integration `pe_byte_identity` | ✅ |
| 4 | Cross-run determinism | Unit `test_clump_only_deterministic` + integration `se_deterministic_across_runs` + CI determinism step | ✅ |
| 5 | Rejection matrix (per-flag) | Integration `rejects_length_flag`, `rejects_adapter_flag`, `rejects_rename_flag`; CLI-side asserts at `cli.rs:685-822` | ⚠️ *see Gap 2* |
| 6 | Short-read preservation | Unit `test_clump_only_no_trimming_short_reads_kept` | ✅ |
| 7 | Adapter preservation | Unit `test_clump_only_no_adapter_detection` | ✅ |
| 8 | Report filename discipline | Unit `test_report_filename_distinct_from_trimming_report` + integration `report_uses_clumping_not_trimming_filename` + inline CI check `ci.yml:687-690` | ✅ |
| 9 | Startup sentinel `Mode: --clump_only` | Written by `clump_only.rs:562`; asserted by integration `report_uses_clumping_not_trimming_filename` (line 242) and CI-implicit via report grep | ✅ (interpreted: sentinel lives in report, not stderr — no other stderr sentinel implemented; matches PROGRESS's design) |

### §Behavior contract invariants

| Invariant | Where | Status |
|---|---|---|
| Header/seq/qual byte-exact per record | `write_records_member` writes via `FastqRecord::write_to`, no mutation; permutation tests confirm | ✅ |
| Multiset equality (permutation, no adds/drops/dedup) | Streaming loop pushes every record; final flush covers un-flushed bins; permutation tests confirm | ✅ |
| Same-input-same-output determinism | Unit + integration + CI determinism checks | ✅ |
| Plus-line and CRLF normalization documented | `cli.rs:200-204` (--help), `docs/.../clump-only.md:38-43` | ⚠️ *see Gap 3* |

### §Rejection matrix (rolled up)

Every entry in the PLAN.md rejection matrix (lines 94-135) is enforced at `cli.rs:685-822`. Spot-checked: adapter flags (`-a`, `-a2`, presets), length/max_length/max_n, trim-n/clip flags, RRBS/non_directional, poly-A/G/no_poly_g, nextseq/2colour, rename, discard_untrimmed, consider_already_trimmed, hardtrim5/3, clock, implicon, demux, clumpify (mutually-exclusive), output-format ubam, passthrough, retain_unpaired. All rejected with per-flag messages. uBAM input rejected inside `clump_only::reject_ubam` (line 233) at `detect_input_format` time — matches PLAN.md line 119.

## Gaps (detail)

### Gap 1 — `--cores >= 2` requirement relaxed to `>= 1` (⚠️ documented deviation)

**Plan line:** PLAN.md 131, 218 ("Enforce `--cores >= 2` (identical to existing `--clumpify` rule)").
**Found:** `cli.rs:692-695` — comment states the deviation ("v1 is single-threaded internally … `--cores` is accepted at any value >= 1 but only affects future parallel implementations").
**Assessment:** Documented deviation in PROGRESS.md §"Documented deviations" (item 1). Rationale is defensible (single-threaded v1 has no use for `>= 2` cores; requiring it would be user-hostile). No technical inconsistency in the deviation itself — byte-identity contract holds regardless of `cores` value.
**Remediation:** No code change required. Consider a one-line PLAN.md footnote pinning the deviation to `feedback_bioconda_pep440_version.md`-style provenance so a future auditor doesn't re-open it.

### Gap 2 — Two unit tests from PLAN.md §Implementation outline step 8 not created (⚠️)

**Plan lines:** PLAN.md 286 (`test_clump_only_normalizes_plus_line`) and 287 (`test_clump_only_normalizes_crlf`).
**Found:** Neither unit test exists. `grep -n "plus_line\|normalizes_plus\|normalizes_crlf\|CRLF\|crlf"` returns only the doc-comment mention at `clump_only.rs:11`.
**Impact:** Low. The normalization behavior is inherited from `FastqReader`/`FastqWriter` (codebase-wide, not `--clump_only`-specific) and is covered by existing FASTQ-layer tests. But the PLAN specifically listed these as regression guards for the documented contract-scope carve-out in Resolved decision 6.
**Remediation:** Add two short unit tests to `clump_only.rs::tests`: (a) write a FASTQ with `+<header-repeat>` on line 3 → assert output has bare `+`; (b) write a FASTQ with `\r\n` line endings → assert output has `\n` endings. Both should take ~15 lines each and reuse `write_synthetic_fastq` + `read_all`. Alternatively, if these belong in `fastq.rs` tests (since the behavior is inherited), reword PLAN.md step 8 to reflect that.

Also, `test_clump_only_ignores_quality_flag` (PLAN.md 285) lives as an integration test (`silently_accepts_quality_flag`, line 331) rather than a unit test. Functionally equivalent, not a real gap — flagging for completeness.

Note: PLAN.md 295 also lists `--clump_only` on uBAM input as an integration test. The unit-test file (`clump_only.rs:1029-1033`) explicitly punts this to `tests/integration_ubam.rs` pattern, but no such uBAM integration test for `--clump_only` was actually added. `reject_ubam` (clump_only.rs:233) is exercised only by internal path; no end-to-end binary test asserts exit≠0 on a uBAM file passed to `--clump_only`. **Recommend adding one integration test** using a small uBAM fixture (or ubam-shaped bytes) to close this loop.

### Gap 3 — `--help` mentions plus-line and CRLF normalization; unit tests don't back it (⚠️, overlap with Gap 2)

**Plan line:** PLAN.md 69 ("Both normalizations are documented in the `--clump_only` help text and docs page").
**Found:** Help text lives at `cli.rs:200-204`; docs page at `docs/.../clump-only.md:38-43`. Both present.
**Gap:** No test enforces that the docs claim is honored. Combined with Gap 2, the "docs promise / no test / no assertion" chain is a soft spot.
**Remediation:** Same as Gap 2 — the plus-line / CRLF unit tests would close both this and Gap 2.

## Verified deviations (from PROGRESS.md §"Documented deviations")

1. **`--cores >= 2` → `>= 1`**: Consistent with §Signature (`cores: usize`, no floor stated) and §Compatible list ("`--cores` (>= 2, same as `--clumpify` today)" contradicts, but the softer stance in §Behavior/§Signature is preserved). No internal PLAN contradiction beyond that single §Compatible line — flag for a one-line PLAN.md update. Not a technical inconsistency; both interfaces produce byte-identical output.
2. **Empty-input emits empty gzip member**: Confirmed at `clump_only.rs:334-336, 488-493` (only when `gzip_output`). Consistent with §Contract's "valid RFC-1952 gzip file" claim. No PLAN contradiction — PLAN was silent on empty-input edge case; this deviation is a strict improvement.
3. **Sort determinism verify-only**: Matches PLAN.md 254-257 (Resolved decision 10). No contradiction.

## Recommendation

Verdict is **INCOMPLETE WITH MINOR GAPS**. The load-bearing invariants (byte-identity, permutation, determinism, rejection matrix) are all enforced by code, unit tests, integration tests, AND CI. Missing items are:

1. **(non-blocking)** Two normalization regression tests (plus-line + CRLF) from PLAN.md step 8 lines 286-287.
2. **(non-blocking)** An end-to-end uBAM-rejection integration test from PLAN.md step 8 line 295.
3. **(optional)** Reword PLAN.md line 131 ("--cores (>= 2 …)") to reflect the documented `>= 1` deviation, so future audits don't re-flag it.

None of these block the PR from opening. All three could be a follow-up commit on the same branch before merge.
