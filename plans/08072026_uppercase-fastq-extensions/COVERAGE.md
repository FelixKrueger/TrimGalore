# Plan Coverage Report

**Mode:** B (code vs. plan — §5 steps, §9 validation, §3.1 edge cases, §13, CHANGELOG)
**Plan(s):** `plans/08072026_uppercase-fastq-extensions/PLAN.md` (v2 + §13)
**Date:** 2026-08-07
**Verdict:** INCOMPLETE — 1 item unresolved

## Summary

- Total items: 29
- DONE: 28
- PARTIAL: 1
- MISSING: 0
- DEVIATED: 0

Branch `fix/384-uppercase-extensions`, uncommitted working tree on `dev` @ `533582a`.
Diff is exactly the 3 files §13 claims. `cargo test`: **538 passed, 0 failed** (exit 0),
matching §13. `cargo fmt --check` clean, `cargo clippy --all-targets --release -D warnings`
clean. The single PARTIAL is one absent named assertion inside an otherwise-complete test;
no behavioural gap.

## Coverage ledger

### §5 implementation outline

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `strip_suffix_ignore_ascii_case` added beside `strip_fastq_extensions`, byte-wise form | §5 Step 1 / §4 / §12 C3 | DONE | `src/io.rs:516-526`, body identical to §4's listing (`checked_sub` + `eq_ignore_ascii_case` on bytes + `.then()`); 3-line doc-comment states the multi-byte-panic reason |
| 2 | `strip_fastq_extensions` uses the helper for both suffix lists; `file_stem` fallback and two-sequential-strips structure intact | §5 Step 2 | DONE | `src/io.rs:548-557`; only the two `find_map`/`strip_suffix` calls changed, list contents and order preserved |
| 3 | `is_gzipped` → `eq_ignore_ascii_case("gz")`, doc-comment extended with the must-agree rationale | §5 Step 3 | DONE | `src/io.rs:29-36`; comment names #384 and the fold-one-not-the-other failure mode |
| 4 | Unit tests V1–V3 in `io.rs` | §5 Step 4 | DONE | 3 tests added; see V1/V2/V3 rows for per-case verification |
| 5 | Integration test in `tests/integration_gzip_non_gz_extension.rs` | §5 Step 5 | DONE | 3 tests added (V5, V4, V8) in the file #374/#382 use |
| 6 | CHANGELOG entry under `#### Changes` | §5 Step 6 | DONE | `CHANGELOG.md:125-148`, under the `#### Changes` heading at line 112 (not `#### Bug fixes` at line 6) |

### §9 validation

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 7 | V1 — 7 named stems fold, retained case preserved | §9 V1 | DONE | `test_strip_fastq_extensions_folds_case`: all 7 named inputs present, plus `SAMPLE.FASTQ.BGZF` (superset). `Sample.FastQ.Gz`→`Sample` and `sample_R1.FQ.GZ`→`sample_R1` pin the no-lowercasing property |
| 8 | V2 — `is_gzipped` true for `.gz`/`.GZ`/`.Gz`; false for `.fastq`, `.bgz`, `.BGZ`, `sample.gz.fastq` | §9 V2 | **PARTIAL** | 6 of 7 named cases asserted. **Lowercase `.bgz` is asserted nowhere** — not in `test_is_gzipped_folds_case` (which has `.BGZ`) nor in pre-existing `test_is_gzipped`. See Gaps |
| 9 | V3 — existing strip tests unchanged; `.bam`/`.txt` fallback un-folded; non-ASCII + shorter-than-suffix names | §9 V3 | DONE | `test_strip_fastq_extensions` / `_bgz` absent from the diff ⇒ unchanged, and pass. `test_strip_fastq_extensions_non_ascii_and_short_names`: `😀` (the §4 panic input, idx 1 mid-character), `é.fq`, `a` (shorter than any suffix), `sample.bam`, `sample.txt` |
| 10 | V4 — the two mixed-spelling refusals, not pure case-variants | §9 V4 | DONE | `uppercase_and_lowercase_spellings_of_one_stem_now_collide` loops `sample.fq.gz` then `SAMPLE.FQ.GZ` against `SAMPLE.FASTQ.GZ`; asserts non-zero exit, `Output path collision` on stderr, and an empty output dir. Both discriminate: pre-fold stems differ (`SAMPLE.FASTQ` vs `SAMPLE`/`sample`), post-fold they collide. Filesystem-independent as required |
| 11 | V5 — literal filenames, no `--dont_gzip`, halves fail independently | §9 V5 | DONE | `uppercase_extension_names_and_compresses_like_lowercase`: no `--dont_gzip`; literal `SAMPLE_trimmed.fq.gz` + `SAMPLE.FASTQ.GZ_trimming_report.txt`, no function-under-test in the expectations; `is_file()` resolution (with dir listing in the failure message) precedes the `1F 8B` byte check |
| 12 | V6 — manual one-offs; §13 records both controls discriminating | §9 V6 / §13 | DONE | §13 table: stem-fold reverted → V1+V5 both FAILED; `is_gzipped` reverted → V2+V5+V8 all three FAILED; both reverted-to-green. V5 failing at path resolution under control 2 recorded as predicted. Committed guard for the compression half correctly identified as V2 |
| 13 | V7 — fmt, clippy, cargo test, §2.1 re-run | §9 V7 | DONE | Re-run in this audit: `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` exit 0, zero warnings (fingerprints post-date `src/io.rs`, so the cached pass covers this source); `cargo test` 538 passed / 0 failed. §2.1 reproduction without `--dont_gzip` recorded in §13 and equivalently covered by V5 |
| 14 | V8 — `--clump_only` report content and stem for `.GZ` input | §9 V8 | DONE | `clump_only_uppercase_gz_reports_gzip_and_ratio`: `SAMPLE_clumped.fq.gz` exists, report contains `gzip`, does **not** contain `plain`, contains `Compression ratio:`. The `!contains("plain")` assertion is what discriminates the label flip; §13 confirms it failed under control 2 |

### §3.1 edge-case table

| # | Input row | Traces to | Status | Notes |
|---|-----------|-----------|--------|-------|
| 15 | `sample.fastq.gz` → unchanged | pre-existing `test_strip_fastq_extensions` | DONE | unchanged and passing |
| 16 | `SAMPLE.FASTQ.GZ` → `SAMPLE_trimmed.fq.gz` | V1 unit + V5 integration | DONE | both halves (stem, compression) asserted |
| 17 | `Sample.FastQ.Gz` → `Sample` | V1 unit | DONE | the retained-case case |
| 18 | `SAMPLE.FQ` → `SAMPLE`, plain in / plain out | V1 unit | DONE | stem asserted; the plain-out half rests on `!is_gzipped("SAMPLE.FASTQ")` — the same branch, `.FQ` itself not asserted |
| 19 | `SAMPLE.FASTQ.BGZ` → `SAMPLE`, still plain | V1 unit + V2 unit | DONE | both directions of A4's asymmetry pinned for the uppercase spelling |
| 20 | `sample.BAM` → `file_stem` fallback, unchanged | §3.1 rationale + V3 test | DONE | test uses lowercase `sample.bam`; `file_stem` is case-agnostic and `.bam` is in neither suffix list, so the uppercase variant cannot diverge |
| 21 | `sample.txt` → fallback, unchanged | V3 test | DONE | asserted |
| 22 | `SAMPLE.FASTQ.GZ` + `sample.fastq.gz` → **already** rejected, fold changes nothing | §3.1 rationale + §12 C2 | DONE | load-bearing row; §12 C2 records the pre-change verification (exit 1 on `dev`, same inode on APFS). V4 explicitly forbids using it as a test — correctly absent |
| 23 | `SAMPLE.FASTQ.GZ` + `sample.fq.gz` → **newly** rejected | V4 test, case 1 | DONE | load-bearing row |
| 24 | `SAMPLE.FASTQ.GZ` + `SAMPLE.FQ.GZ` → **newly** rejected | V4 test, case 2 | DONE | load-bearing row; not case-dependent, so it holds on case-sensitive filesystems too |

### §13 notes and CHANGELOG (§7)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 25 | §13 "all six §5 steps done as specified — no deviations" holds against the diff | §13 | DONE | `git status`/`git diff` touch exactly `src/io.rs`, `tests/integration_gzip_non_gz_extension.rs`, `CHANGELOG.md`. No undocumented deviation found: helper body, both folds, doc-comment, test placement and test count all as §5/§4 specify. Bookkeeping nit only — §13 says io.rs `+70/−3`; `git diff --numstat` gives `67/3` (70 is the `--stat` total). CHANGELOG `+24` and tests `+111` exact |
| 26 | §7's instruction to re-confirm by grep that the validation matrix is unaffected | §7 | DONE | Verified in this audit rather than in §13: no uppercase FASTQ/FQ/GZ extension in `.github/workflows/ci.yml` or `test_files/`. Requirement holds; §13 does not record having run the grep |
| 27 | CHANGELOG has four bullets: rename, compression flip, new refusals, `--clump_only` | §7 / §12 | DONE | all four present as sub-bullets under the #384 entry |
| 28 | Perl-parity departure stated (B's I2) | §7 / §12 | DONE | "Perl v0.6.11's compression test was case-sensitive, so this knowingly departs from v0.6.11 for uppercase names" |
| 29 | No MultiQC-grouping claim (B's I4) | §7 / §12 | DONE | Entry claims only internal consistency and explicitly notes the report filename is unchanged. No mention of MultiQC or downstream grouping anywhere in the diff |

## Gaps (detail)

### Item 8: V2's lowercase `.bgz` false-assertion

**Expected:** V2 names seven cases for `is_gzipped` — true for `.gz`, `.GZ`, `.Gz`; false for
`.fastq`, `.bgz`, `.BGZ`, `sample.gz.fastq`.

**Found:** Six. Every `is_gzipped` assertion in `src/io.rs` (lines 967-976 pre-existing,
983-989 new) was enumerated: no assertion takes a lowercase `.bgz` path. The pre-existing
`test_is_gzipped` covers `.gz` true and `.fastq`/`sample.gz.fastq` false but never `.bgz`;
the new `test_is_gzipped_folds_case` covers `.BGZ` only.

**Gap:** one assertion, e.g. `assert!(!is_gzipped(Path::new("sample.fq.bgz")));` in
`test_is_gzipped_folds_case`.

**Severity:** documentation-level. `.BGZ` exercises the identical branch, so A4's
"`.bgz` does not imply gzipped output" asymmetry is already guarded against regression —
just for the uppercase spelling rather than the lowercase anchor the plan named. No
behaviour is unverified.

## Test verification

| Test | File | Status |
|------|------|--------|
| `io::tests::test_strip_fastq_extensions_folds_case` (V1) | `src/io.rs:982` | PASS |
| `io::tests::test_is_gzipped_folds_case` (V2) | `src/io.rs:982` | PASS (1 named case absent) |
| `io::tests::test_strip_fastq_extensions_non_ascii_and_short_names` (V3) | `src/io.rs:1015` | PASS |
| `io::tests::test_strip_fastq_extensions` (V3 regression) | `src/io.rs` | PASS, unchanged |
| `io::tests::test_strip_fastq_extensions_bgz` (V3 regression) | `src/io.rs` | PASS, unchanged |
| `io::tests::test_is_gzipped` (V3 regression) | `src/io.rs:965` | PASS, unchanged |
| `uppercase_extension_names_and_compresses_like_lowercase` (V5) | `tests/integration_gzip_non_gz_extension.rs:352` | PASS |
| `uppercase_and_lowercase_spellings_of_one_stem_now_collide` (V4) | `tests/integration_gzip_non_gz_extension.rs:388` | PASS |
| `clump_only_uppercase_gz_reports_gzip_and_ratio` (V8) | `tests/integration_gzip_non_gz_extension.rs:421` | PASS |
| **Suite total** | — | **538 passed, 0 failed, 0 ignored** (exit 0) |

Gates: `cargo fmt --all -- --check` clean · `cargo clippy --all-targets --release -- -D warnings`
exit 0, zero warnings · `cargo test` exit 0.

## Verdict

**INCOMPLETE — 1 item unresolved.**

One item remains, and it is a single line of test code:

1. **Item 8 (V2)** — add the lowercase `.bgz` false-assertion to
   `test_is_gzipped_folds_case` in `src/io.rs`, e.g.
   `assert!(!is_gzipped(Path::new("sample.fq.bgz")));`

Everything else in §5, §9, §3.1, §13 and the §7 CHANGELOG requirements is DONE. In
particular the load-bearing §3.1 rows (items 22-24) are all accounted for — the
already-rejected pair by verified rationale, both newly-rejected pairs by V4 — §13's
"no deviations" claim holds against the diff, and the CHANGELOG carries four bullets with
the Perl compression-parity departure stated and no MultiQC-grouping claim.
