# Plan Coverage Report — `--phred64` uBAM fixes (issue #358)

**Mode:** B (code vs. plan)
**Plan:** `plans/phred64-ubam/PLAN.md` (v2)
**Implementation:** uncommitted working tree on `fix/phred64-ubam`, base `dev` @ `2e3da11`
**Date:** 2026-07-25
**Verdict:** **COMPLETE** — all 12 implementation steps landed, all 13 validation checks accounted for, both declared deviations are real and correctly scoped. 4 minor observations, 0 functional gaps.

---

## Summary

| | Count |
|---|---|
| Total items audited | 46 |
| DONE | 41 |
| PARTIAL | 1 |
| DEVIATED (declared) | 2 |
| DEVIATED (undeclared, non-functional) | 2 |
| MISSING | 0 |

Independently reproduced, not taken from §13: `cargo fmt` clean · `cargo clippy --all-targets --release -- -D warnings` clean · `cargo test` **418 passed / 0 failed** (359 + 12 + 15 + 1 + 2 + 8 + 21) · baseline at `2e3da11` = **409** `#[test]` · **zero** test functions removed · diffstat **9 files, +500/−14** plus the untracked fixture — exactly as §13 claims.

---

## 1. §5.2 — the 12 implementation steps

| Step | Item | Status | Evidence |
|---|---|---|---|
| 1 | Guard `cli.phred64 && any_bam` in the T17 block | **DONE** | `src/main.rs:211-220`, sited at `:211` — after `any_bam` (`:187`), before `naming::ensure_output_dir` (`:266`) and before every dispatch branch (specialty `:308`). Error text is the plan's verbatim. |
| 2 | `input_phred_offset: u8` on `struct BamWriter` | **DONE** | `src/bam.rs:521` with the §4 doc comment cross-referencing `PHRED_OFFSET`. |
| 3 | `create` gains the parameter + doc | **DONE** | `src/bam.rs:547-551` doc block; param at `:566`. **The pre-existing "PLAN §4 deviation" paragraph was preserved verbatim and extended** ("A 5th, `input_phred_offset`, was added later for issue #358") — this was an explicit plan requirement that v1's block silently dropped. |
| 4 | Fix + comment rewrite correcting **both** falsehoods | **DONE** | `src/bam.rs:600-614`. Neither "FastqRecord.qual is Sanger ASCII (+33)" nor "FastqReader rejects sub-33 input" survives. The replacement's factual claim is **independently verified**: `FastqReader::sanity_check` (`src/fastq.rs:477-495`) checks only the `@` prefix and colorspace; the only `qual`-related `bail!`s in the file are truncation (`:382`, `:416`). No quality-range validation exists. |
| 5 | One-line comment marking `PHRED_OFFSET` read-side-only | **DONE** | `src/bam.rs:60-66` — 7 lines rather than one, naming the conflation as the defect being fixed. Exceeds spec. |
| 6 | `main.rs` 3 direct `create` sites pass `cli.phred_offset()` | **DONE** | `:1852`, `:1994`, `:2126` — enclosing fns `run_ubam_output_single`, `…_paired_two_files`, `…_paired_single_file`, all of which already took `cli: &Cli`. |
| 7 | `specialty.rs` hardtrim5/3 gain the param; callers updated | **DONE** (naming deviation, see §7 below) | `src/specialty.rs:113`, `:163` signatures; `:134`, `:184` create args; callers `src/main.rs:327`, `:352`. Params are **7**, and no `#[allow(clippy::too_many_arguments)]` was added — the plan's prediction that none is needed holds (clippy passes clean). |
| 8 | `clump_only.rs` two fns + 3 main callers + 1 test caller | **DONE** | `src/clump_only.rs:745`, `:895` signatures; `:772`, `:989` create args; callers `src/main.rs:518`, `:585`, `:617`; test caller `src/clump_only.rs:1546`. |
| 9 | `bam.rs:920` `+` → `saturating_add` | **DONE** (declared residual) | `src/bam.rs:965` now `b.saturating_add(PHRED_OFFSET)`; the all-`0xFF` sentinel branch above is untouched. Residual documented in place at `:957-961`. |
| 10 | O-3 startup `NOTE:` in the existing uBAM block | **DONE** | `src/main.rs:291-304`, inside the `OutputFormat::UBam` block, after the `--cores` and sentinel NOTEs. **Verified firing** by direct invocation. |
| 11 | CHANGELOG under `### Unreleased` → append to existing `#### Fixes`, plus a distinct behaviour-change line | **DONE** | `#### Fixes` at `CHANGELOG.md:55` (pre-existing, appended at `:64`); new `#### Behaviour changes` at `:84`, still inside `### Unreleased` (`:4`) and above `### Version 2.3.0` (`:111`). Names both regressing invocations explicitly, and carries the §11 "archived output cannot be repaired by re-running" paragraph. |
| 12 | Docs → `guide/quality.md` "Phred encoding" | **DONE** | 4 paragraphs at `docs/src/content/docs/guide/quality.md:33-42`. `outputs.md` is **unmodified** and contains no `phred` mention — the plan's prohibition on filing the guard under its "Rejected at CLI-validate time" list is honoured. The optional cross-reference was skipped, which the plan permits. |

---

## 2. §5.1 — the ordering constraint (guard before writer fix)

**Status: DONE by construction; ⚠️ §13's standalone-verification claim is not independently verifiable.**

The branch has **zero commits** — `git rev-parse HEAD` is `2e3da11`, the base — and `git reflog` shows only the branch checkout. There is therefore no git artifact recording the sequence in which the two halves were written, and I cannot corroborate §13's "the guard was implemented and verified in isolation before the writer was touched."

What I *can* establish, and what the constraint actually requires:

1. **The dangerous intermediate state does not exist in the repository.** The entire change is a single uncommitted diff and will become one commit. §5.1's operative requirement — "no intermediate commit may silently corrupt a working path" — is satisfied structurally, independently of authoring order.
2. **§5.1's "no compile dependency either way" claim is confirmed.** The guard reads only `cli.phred64` and `any_bam`, both pre-existing at `:187`; it touches no `bam.rs` symbol. It was implementable and testable standalone.
3. **The silent all-Q0 path is closed in the shipped state.** Verified live on all four BAM-input routes:

   | Invocation | Exit | Output dir created |
   |---|---|---|
   | `--phred64 -q 20 <bam>` | 1 | **no** |
   | `--phred64 --output-format ubam <bam>` | 1 | **no** |
   | `--clump_only --phred64 --output-format ubam <bam>` | 1 | **no** |
   | `--hardtrim5 20 --phred64 <bam>` | 1 | **no** |

Recommendation for the eventual commit: keep it as **one** commit, or if split, the guard must be the first. Nothing in the current state prevents this.

---

## 3. §9 — the 13 validation checks

| # | Check | Plan location | Status | Evidence / note |
|---|---|---|---|---|
| 1 | Bug 1, trim SE | `tests/integration_ubam_out.rs` | **DONE** | `phred64_ubam_out_se_stores_true_phred`. Expectation is `[40×24, 2×10]` rather than `[40;34]` — consequence of the declared fixture deviation. Passes `-q 0 --length 0` so the Q2 tail reaches the writer. |
| 2 | Phred+33 unchanged (raw 40 from `'I'`) | same file | **DONE**, location deviated | Implemented as the second assertion of the `src/bam.rs` unit test (`stored_qual("IIII", 33, …)` → `vec![40;4]`), not as an integration test. The plan's own check 10 row already specified this same assertion at unit level, so nothing is lost; integration-level Phred+33 qual-byte equality is separately covered by the two `assert_ubam_eq` golden tests. **The classification is honoured explicitly in-code**: "an over-correction guard, not a bug guard: this already passed before the fix and must keep passing." |
| 3 | Bug 1, specialty | `integration_ubam_out.rs` | **DONE** | `phred64_ubam_out_specialty_stores_true_phred` — asserts 20 bytes, all raw 40. |
| 4 | Bug 1, clump-only | `integration_clump_only_ubam.rs` | **DONE** | `phred64_clump_only_ubam_out_stores_true_phred`. |
| 5 | Guard: BAM in, FASTQ out | `integration_ubam_out.rs` | **PARTIAL** | `phred64_bam_input_rejected_fastq_out` asserts non-zero exit and that stderr contains both `--phred64` and `BAM`. **The plan's third clause — "no output dir created" — is not asserted.** The `fresh_tmpdir` helper (`tests/integration_ubam_out.rs:27-32`) pre-creates the directory with `create_dir_all`, so the property is untestable through it. I verified the property holds in reality (table in §2 above); it is simply unguarded against regression. |
| 6 | Guard: BAM in, uBAM out | same | **DONE** | `phred64_bam_input_rejected_ubam_out`. |
| 7 | Guard: mixed input | same | **DONE** | `phred64_mixed_input_rejected` (`--paired phred64_test.fastq ubam_test.bam`). Asserts exit only — which is all the plan's row required. |
| 8 | Guard: clump-only + BAM + uBAM out | `integration_clump_only_ubam.rs` | **DONE** | `phred64_clump_only_bam_input_rejected`, with a comment explaining why this is the sharpest case. |
| 9 | Guard: early-returning mode | `integration_ubam_out.rs` | **DONE** | `phred64_bam_input_rejected_early_returning_mode`, and the test comment cites the B-NIT-2 precedent the plan's check-9 rationale invoked. |
| 10 | Unit: offset arithmetic | `src/bam.rs` tests | **DONE** | `bam_writer_subtracts_input_phred_offset`. Reads back via `bam::io::Reader` + `Record::quality_scores()` — the exact layer §9 prescribed, with the in-file precedent named in the doc comment. Pins **both** branches (`'h'`/64 → 40, `'I'`/33 → 40). The "cannot compile pre-fix" property holds literally: the call passes 5 arguments. |
| 11 | FastQC on uBAM | `integration_ubam_out.rs` **or manual** | **⚠️ MANUAL only** | No automated test. I reproduced it: `--phred64 --output-format ubam --fastqc -q 0 --length 0` → `Encoding  Sanger / Illumina 1.9`, base 1 mean `40.0`, bases 25-34 mean `2.0`. Matches §13 exactly. The row's "or manual" permits this; note only that §9's header sentence counts check 11 among the integration tests, so the header's "7 integration tests" is 8 automated integration tests + this one manual. |
| 12 | No regression | `cargo test` / clippy / fmt | **DONE** | fmt clean; clippy `-D warnings` clean; **418 passed, 0 failed**; baseline 409 confirmed by `git grep '#[test]' 2e3da11`; a name-level diff of every baseline test fn against the current tree shows **0 removed**, 9 added (1 unit + 6 uBAM-out + 2 clump-only) — matching §13's breakdown. |
| 13 | Golden fixtures unchanged | `assert_ubam_eq` references | **DONE** | `git status test_files/` shows only `README.md` modified and `phred64_test.fastq` untracked; `ubam_out_se_REFERENCE.bam` and `ubam_out_pe_REFERENCE.bam` are byte-untouched, and both consuming tests pass. §7's no-op claim holds. Correctly *not* expected to fail pre-fix. |

**Pre-fix-failure classification:** matches the plan. Checks 2 and 13 are the two non-failing guards, and both are labelled as such in the implementation (check 2 in the test comment, check 13 by leaving the fixtures alone). I did not re-run the suite against pre-fix code — deliberately, because two code reviewers are working the same working tree and stashing would disturb them. The classification is verifiable by inspection: pre-fix, checks 1/3/4 would read 71/33, checks 5-9 would exit 0, and check 10 would not compile.

---

## 4. §8 — assumptions

| # | Assumption | Status |
|---|---|---|
| 1 | BAM `QUAL` holds true Phred 0..93 | **DONE** — the fix's premise, encoded in `write_record`. |
| 2 | `BamReader` normalises to Phred+33 and stays that way | **DONE** — read side unchanged apart from step 9's `saturating_add`; `PHRED_OFFSET` remains a fixed 33 and is now commented as such. |
| 3 | FASTQ output stays encoding pass-through | **DONE** — `src/fastq.rs` is **unmodified**; no offset arithmetic added. Now also stated in user docs. |
| 4 | `--phred64` describes input encoding only | **DONE** — asserted in the docs, the O-3 NOTE and the guard's error text. |
| 5 | Hard error, not warn-and-ignore | **DONE** — `anyhow::bail!`; exit 1 on all four routes. |
| 6 | Required `create` parameter, not a defaulted setter | **DONE** — positional 5th param; all 14 call sites compiler-enumerated. §13 records two sites missed by a first regex pass and caught by the compiler, which is precisely the property §4 chose this design for. |
| 7 | Test sites pass literal `33`, not the const | **DONE** — all 7 pre-existing test `create` calls pass `33`; `PHRED_OFFSET` is not reused on the write side. |
| 8 | **RESOLVED**: four sites need one added parameter, zero intermediate threading | **DONE — reality matched exactly.** Verified by enclosing-function analysis: all 5 wrapper callers (`main.rs:327, 352, 518, 585, 617`) sit directly inside `fn main()` with `cli` in scope, and the 3 `run_ubam_output_*` fns already took `cli: &Cli`. **No intermediate function required a new parameter.** |

---

## 5. §10 — the three open questions

| | Disposition in plan | Status |
|---|---|---|
| **O-1** `--phred33 --phred64` together silently yields 64 | deferred, leave as-is | **DONE — nothing changed.** `src/cli.rs` is entirely unmodified: no `conflicts_with` added, `phred_offset()` still `if self.phred64 { 64 } else { 33 }` (`:970-972`). §3.1's benign-interaction edge case verified live: `--phred33 --phred64 <bam>` → exit 1 (guard fires on `self.phred64`). |
| **O-2** no extra rejection for `--clump_only` + FASTQ | deferred, leave as-is | **DONE — nothing changed.** `--clump_only --phred64 <fastq>` still exits 0 and writes `*_clumped.fq`. The §3.1 governing principle (reject on input *format*, not on mode inertness) is reproduced verbatim in the guard's code comment, so the reconciliation is recorded where a future reader will find it. |
| **O-3** adopted — startup `NOTE:` | must land with the plan's SAM-spec wording | **DONE.** Landed at `main.rs:298-303` as a superset of the plan's suggested sentence: *"NOTE: input is Phred+64 (ASCII+64); output BAM QUAL stores true Phred scores (0–93) per the SAM spec, not ASCII. Verify the input really is Phred+64 — running Phred+33 data with `--phred64` yields an all-zero-quality BAM."* The plan's wording was prefixed "e.g.", so the addition is within scope; the added second sentence is what actually delivers O-3's stated *rationale* (mitigating the silent all-Q0 floor). The reporting-ambiguity resolution ("Quality encoding type selected: ASCII+64" describes the input) is captured in the code comment and, user-visibly, in `quality.md:37`. |

Also verified: all four §3.1 edge cases behave as specified — mixed input rejected, FASTQ-only unaffected (exit 0), BAM-without-`--phred64` unaffected (exit 0), `--phred33 --phred64` + BAM rejected.

---

## 6. §7 — integration claims

| Claim | Where the plan said it must land | Status |
|---|---|---|
| Guard preserves `clump_only`'s lossless-qual invariant | code + **PR body** | **DONE in code**, `PENDING` for the PR body. The rationale is recorded in the guard's comment (`main.rs:196-198`) and, at length, in `phred64_clump_only_bam_input_rejected`'s test comment; the CHANGELOG behaviour-change entry states it would "silently write an all-zero-quality BAM while reporting success, breaking that mode's documented lossless guarantee." |
| FastQC on uBAM is an unlisted beneficiary | validation row + **CHANGELOG clause** | **DONE.** CHANGELOG: "`--fastqc` on uBAM output was affected too — per-base-quality plots were centred on Q71." Validation row is check 11 (manual — see §3). |
| Byte-identity invariants untouched; v0.6.11 matrix out of scope | **PR body** | **DONE in substance**, `PENDING` for the PR body. Confirmed mechanically: `src/fastq.rs` and `src/quality.rs` are both unmodified, and the only ASCII↔raw conversions remain the two in `bam.rs`. |
| No golden-fixture churn; grep argument self-sufficient | verified claim + **PR-body pre-emption of the regen-trigger question** | **DONE in substance**, `PENDING` for the PR body. Re-verified: `.github/` still has **zero** `phred64` hits, so no CI path passes the flag; the only new hits are the two test files. Reference BAMs unmodified and their tests pass. |
| Withdrawn fixture-encoding heuristic; §9's consequence — do not assert a blanket Phred+33 claim, list `truncated.fq.gz` as the exception | `test_files/README.md` | **DONE, and fully.** A new "Quality encoding" section lists both non-Phred+33 fixtures with provenance, and closes with an explicit prohibition: "Do not assert 'all fixtures are Phred+33' anywhere — the two above are counter-examples." This is exactly the §7 caveat, discharged rather than merely avoided. |

⚠️ **Four PR-body items are unclosed** simply because the branch is unpushed and no PR exists: the one-PR-is-a-correctness-constraint statement (§5.1 item 2), the `clump_only` invariant note, the byte-identity statement, and the fixture-regen-trigger pre-emption. There is no PR-notes artifact in `plans/phred64-ubam/`. These are a checklist for whoever opens the PR, not implementation gaps.

---

## 7. §13 — declared deviations, and the search for undeclared drift

### The two declared deviations are real and correctly scoped

**D-1 — fixture composition (`'h'` uniform → 24×`'h'` + 10×`'B'`).** Confirmed: `test_files/phred64_test.fastq` is 4 reads × 34 bp, `hhhhhhhhhhhhhhhhhhhhhhhhBBBBBBBBBB`. The stated cause is reproducible — I ran check 11 on the shipped fixture and got `Sanger / Illumina 1.9` with means 40.0/2.0, so `fastqc-rust`'s minimum-byte heuristic is satisfied. The deviation is documented in three places (§13, `test_files/README.md`, and the test comments), including the non-obvious `-q 0 --length 0` requirement without which only one offset is exercised. This deviation **strengthens** the plan's intent rather than diluting it: check 1 now pins two distinct quality values instead of one.

**D-2 — `bam.rs` mixed-sentinel residual documented, not fixed.** Confirmed: `saturating_add` landed (the overflow is gone); the residual — a mixed-sentinel byte saturating to 255 rather than mapping to `QUAL_MISSING_REPLACEMENT` — is commented in place at `src/bam.rs:957-961` with the reason for deferring. Step 9's stated goal was "remove the overflow," which is met; the residual is a scope boundary correctly drawn and disclosed.

### Undeclared drift found: 2 items, both non-functional

**U-1 — parameter naming.** §5.2 steps 7/8 specified `phred_offset: u8` for the four `specialty.rs`/`clump_only.rs` functions; they were implemented as `input_phred_offset: u8`. Consistent with the struct field §4.1 deliberately named `input_phred_offset` to defuse the `PHRED_OFFSET` collision, so this is the plan's own §4.1 reasoning applied one layer out. Trivial, and an improvement — but not declared.

**U-2 — check 2's layer.** Moved from `tests/integration_ubam_out.rs` to the `src/bam.rs` unit test. No assertion lost (check 10's row already called for the same `'I'` → 40 assertion at that layer), but the relocation is undeclared.

### Explicitly permitted latitude, correctly used — not drift

- `outputs.md` cross-reference skipped — step 12 marked it "Optionally", and critically the guard was **not** filed under the "Rejected at CLI-validate time" list, which the plan warned would state something false.
- No CI step added — §9 marked CI "optional" and argued the integration tests cover the wiring more precisely.
- §3.2's Solexa/Illumina-1.0 flooring note appears only in the plan, not in docs or CHANGELOG. The plan said it "belongs in the record" without naming a target, and required no change. The docs caution covers the adjacent Phred+33-mistake case. **Informational only** — if the author wants it user-visible, `quality.md`'s caution block is the natural home.

### Nothing else changed

I read the complete diff of all 9 modified files. Every hunk maps to a numbered §5.2 step, a §9 validation check, or the §4/§4.1 documentation requirements. No incidental refactoring, no unrelated files, no drive-by edits.

---

## 8. Counts

| Claim | Source | Verified |
|---|---|---|
| 7 production `BamWriter::create` sites | §2 | ✅ `main.rs:1852/1994/2126`, `specialty.rs:134/184`, `clump_only.rs:772/989` |
| 7 test call sites | §2 | ✅ 6 pre-existing in `bam.rs` (`:1395, 1445, 1473, 1528, 1568, 1744`) + 1 in `clump_only.rs:1546`; a 7th `bam.rs` site (`:1627`) is the *new* unit test, so the "6 updated + 1" arithmetic is right |
| 12 production edit points, not 7 | §13 | ✅ 8 new `cli.phred_offset()` sites in `main.rs` (5 wrapper callers + 3 direct) + 2 `specialty.rs` + 2 `clump_only.rs` = 12. `grep -c 'cli.phred_offset()' src/main.rs` = 12 total, of which 4 (`:929, 1148, 1906, 2271`) are pre-existing report-path uses |
| 19 total edit points | §13 | ✅ 12 + 7 |
| 409 baseline tests → 418 | §9 / §13 | ✅ both confirmed at the source |
| 9 files, +500/−14 | §13 | ✅ exact |

§13's correction of §2's headline figure is itself accurate, and the accounting is transparent about what it counts (call sites, not the 4 signature edits inside the same 4 functions). **Nothing was missed.**

---

## Test verification

| Test | File | Status |
|---|---|---|
| `bam_writer_subtracts_input_phred_offset` | `src/bam.rs:1616` | PASS |
| `phred64_ubam_out_se_stores_true_phred` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_ubam_out_specialty_stores_true_phred` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_bam_input_rejected_fastq_out` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_bam_input_rejected_ubam_out` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_bam_input_rejected_early_returning_mode` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_mixed_input_rejected` | `tests/integration_ubam_out.rs` | PASS |
| `phred64_clump_only_ubam_out_stores_true_phred` | `tests/integration_clump_only_ubam.rs` | PASS |
| `phred64_clump_only_bam_input_rejected` | `tests/integration_clump_only_ubam.rs` | PASS |
| `ubam_out_se_ubam_input_matches_reference` (check 13) | `tests/integration_ubam_out.rs` | PASS, fixture unmodified |
| `ubam_out_pe_one_ubam_interleaved_matches_reference` (check 13) | `tests/integration_ubam_out.rs` | PASS, fixture unmodified |
| Full suite | — | 418 passed / 0 failed |
| FastQC-on-uBAM (check 11) | — | **no automated test**; reproduced manually |

---

## Verdict

**COMPLETE.** Every §5.2 step landed. Every §9 check is accounted for, and the pre-fix-failure classification (checks 2 and 13 as non-failing guards) is honoured explicitly in the implementation. §8's Assumption 8 was resolved correctly and reality matched it exactly. Both §10 deferrals are genuinely untouched, and O-3 landed with the required SAM-spec wording. §13's self-report is accurate throughout — every claim I could check independently checked out, including the corrected counts and the FastQC numbers.

Nothing requires fixing before this can proceed. The following are optional hardening / housekeeping, in descending order of value:

1. **Check 5's "no output dir created" clause is unasserted** (the only PARTIAL). The property holds today but nothing guards it, and it is the clause that proves the guard sits ahead of `ensure_output_dir` — a placement the plan itself flagged as historically error-prone in this file (B-NIT-2). A test would need a path the helper has not pre-created, e.g. `fresh_tmpdir("…").join("nested")` passed to `-o`, asserting `!path.exists()` after the failed run.
2. **Check 11 has no automated test.** Permitted by its row, but it is the one user-visible artefact (`fastqc_data.txt` means) that no unit test can reach, and the fixture deviation exists *because* of it.
3. **Four §7/§5.1 statements are assigned to a PR body that does not exist yet** — the one-PR correctness constraint, the `clump_only` lossless-invariant note, the byte-identity/v0.6.11-matrix statement, and the fixture-regen-trigger pre-emption. Worth capturing before the PR is opened; the plan predicts each will be asked about.
4. **Two undeclared micro-deviations** (parameter name `input_phred_offset`; check 2 moved to the unit layer) are worth one line each in §13 for completeness. Neither changes behaviour.
5. **Keep this as a single commit** (or guard-first if split) so §5.1's ordering constraint remains satisfied in history as well as in the working tree.
