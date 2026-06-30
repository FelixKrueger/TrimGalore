# Code Review (Round 2 — fix re-review) — Phase 1 uBAM output (`--output-format ubam`) — Reviewer B

**Target:** working-tree diff (`git diff HEAD` + untracked `tests/integration_ubam_out.rs`, `test_files/ubam_out_{se,pe}_REFERENCE.bam`)
**Spec:** PLAN v2.1 (`.../phase1-trimgalore-formats/PLAN.md`)
**Original reports:** `CODE_REVIEW_reviewer-A.md`, `CODE_REVIEW_reviewer-B.md`
**Reviewer:** B (independent re-review; ran in parallel with Reviewer A, no coordination)
**Date:** 2026-06-26

## Verification gate (re-run locally)

| Check | Result |
|---|---|
| `cargo build --release` | **clean** (only sandbox `xcrun` FS-event noise, not code; exit 0) |
| `cargo test --release` | **349 pass / 0 fail** — `326 lib + 2 passthrough + 8 ubam-in + 13 ubam-out`. Matches the implementer's claim exactly. |
| `cargo clippy --release --all-targets -- -D warnings` | **clean** (exit 0) |
| `cargo fmt --all -- --check` | **clean** (exit 0) |

Test count delta from round 1: 336 → 349 (+13). The +9 lib tests (C1 splice ×3, @PG disambiguation ×1, seq-validation ×5) and +4 integration tests (cores>1, mixed-batch, rename+tags, two-bam-pair) account for the increase.

## Fixes verified

| Original finding | Fix location | Verdict | Notes |
|---|---|---|---|
| **C1 / A-C1** — `--rename` corrupts last preserved tag | `src/fastq.rs:148-156` (`append_to_id`) | **PASS** | Splices suffix *before* the first `\t` when a tag tail exists; appends to end otherwise. 3 lib tests (`fastq.rs:589/600/615`) cover both branches + whitespace trim. Integration guard `ubam_out_rename_with_preserve_tags_keeps_tags_intact` asserts UB value has no `:clip5:`. Clipped seq is ACGTN so it can never introduce a stray tab — splice is safe. Covers the hardtrim `--rename` path too (`specialty.rs:140,189` use the same fn). |
| **I1 / AGREE-1** — duplicate/overwriting `@PG` on re-processing | `src/bam.rs:622-642` (`build_output_header`) | **PASS** | Switched to `h.programs_mut().add("trim_galore", pg)` with `.with_context`; `last_pg_id` helper dropped; PG built without PP so noodles computes the chain. Test `build_output_header_disambiguates_existing_trim_galore_pg` (`bam.rs:1586`) asserts 2 entries `["trim_galore", "trim_galore-trim_galore"]`. Note: the disambiguation suffix is noodles' `-<prev>` scheme, not samtools' `.1`/`.2` — Reviewer A suggested `.1` but the PLAN/doc-comment now explicitly documents `trim_galore-trim_galore`, so this is an intentional, internally-consistent choice. |
| **I2 / B-I2** — two BAM files silently accepted under `--paired --output-format ubam` | `src/main.rs:1625-1635` (in `run_ubam_output_paired_two_files`) | **PASS** | Per-input `detect_input_format` check bails with the same message as the FASTQ path. Integration guard `ubam_out_two_bam_pair_rejected` asserts non-zero exit + message. (Minor: the guard fires *after* the caller's collision pre-flight loop at `main.rs:1418`, so it's per-pair-late rather than a single up-front gate — functionally correct, see NIT-1.) |
| **I2 / A-I2** — `H:` hex tags emitted on read, rejected on write (asymmetric, mid-stream fatal) | `src/bam.rs:957-968` (`append_tag_type_and_value`) | **PARTIAL** | Read side now `bail!`s on `Value::Hex(_)`, symmetric with the existing `Value::Array(_)` rejection — correct fail-fast-at-input behaviour. **However the read-side bail has no regression test** (see NEW-1). The only H-tag test (`parse_h_hex_tag_rejected`, `bam.rs:1273`) exercises the *write-side parser* `parse_name_and_data`, not the read-side emitter the fix touched. Fix is correct; reverting it would not break any test. |
| **M1 / B-M1** — write path doesn't enforce ACGTN-only / reject `=` / IUPAC→N | `src/bam.rs:760-786` (`validate_and_normalize_seq_for_write`), called from `write_record:577` | **PASS** | New helper mirrors the read-side `bam_record_to_fastq` validation: uppercases, ACGTN pass, IUPAC→N (warn-once), `=` and garbage bytes bail. 5 lib tests (`bam.rs:1546-1584`) cover all five arms. Covers the specialty hardtrim BAM paths too (they route through `write_record`). |
| **M3 / B-M3** — missing-qual stderr NOTE not emitted | `src/main.rs:1387-1397` (`run_ubam_output`), `any_bam` threaded from `main:144`/`297` | **PASS** | One-shot NOTE printed when `any_bam` is true. Wording rewritten vs the PLAN draft but conveys the same Phred-0/0xFF semantics. (Minor over/under-fire edge: see NIT-2.) |
| **AGREE-2 + §9 gaps** — missing `--cores>1` warn-and-proceed + mixed-batch tests | `tests/integration_ubam_out.rs:311,340,394,443` | **PASS** | All four promised tests present, each a genuine guard: `cores_gt_1` asserts exit 0 + "ignored" + output exists; `mixed_batch` asserts uBAM-side has CB and FASTQ-side does not; `rename_with_preserve_tags` is the C1 binary-level guard; `two_bam_pair_rejected` is the B-I2 binary-level guard. All referenced fixtures (`ubam_test_with_tags.bam`, `ubam_test.bam`, `ubam_paired_test.bam`) are git-tracked, so the `SKIP`-if-absent branches will not silently no-op in CI. |

**All P0/P1/P2 findings from round 1 are resolved.** No fix introduced a build/clippy/fmt/test regression.

---

## NEW findings

### IMPORTANT — none.

The fix patch does not introduce any new correctness defect, regression, or missed edge case at the IMPORTANT level. The headline C1 corruption is genuinely closed and guarded at both the unit and binary level.

### NIT

**NEW-1 — A-I2 read-side `H:` rejection has no regression test.**
`src/bam.rs:957` (`append_tag_type_and_value`, `Value::Hex(_)` arm). The fix is correct (fail-fast at read instead of mid-stream-fatal at write), but no test feeds a Hex-typed BAM aux through the read path to assert the bail. `parse_h_hex_tag_rejected` (`bam.rs:1273`) only covers the *write-side parser*. If the read-side `bail!` were reverted to the old "emit `H:` into the id tail" behaviour, the suite would stay green. Low-impact (hex tags are rare; CB/UB/RX are all `Z`), but it is the one fix whose guard test is absent. Recommend a small lib test constructing a `Value::Hex` and asserting `append_tag_type_and_value` returns `Err`.

**NIT-1 — two-BAM rejection (B-I2) fires per-pair, after the caller's collision pre-flight.**
`src/main.rs:1413-1435` runs `paired_bam_output_name` + a case-folded collision HashMap over *all* input pairs first; only then does each `run_ubam_output_paired_two_files` call (line 1449) hit the two-BAM bail at `:1625`. For two BAM inputs this means the collision loop and an extra `detect_input_format` run before the rejection. Functionally correct (the run still errors before any record is written), and it costs only microseconds, but a single up-front guard in `run_ubam_output` before the dispatch (mirroring the FASTQ path's `main.rs:911-921`) would be cleaner and fail earlier. Not blocking.

**NIT-2 — M3 NOTE fires on `any_bam`, slightly over- and under-eager.**
`src/main.rs:1387`. The NOTE keys on `any_bam` (true if *any* input is uBAM):
- **Over-eager:** in a mixed FASTQ+uBAM batch, or when the uBAM input has complete quality scores, the "missing-qual sentinel does not round-trip" NOTE still prints even though no record actually had a missing-qual sentinel. It is an informational heads-up, not a per-record claim, so this is benign — but it will print on every uBAM-input run regardless of whether any qual was actually missing.
- **Under-eager:** the specialty BAM paths (`hardtrim5_to_bam`/`hardtrim3_to_bam`) `return` at `main.rs:217`/`241` *before* reaching `run_ubam_output`, so the NOTE never prints for `--hardtrim5/3 --output-format ubam` on a uBAM input — even though those paths have the identical missing-qual round-trip behaviour. The PLAN §3.3 step 4 wording says "when uBAM-in → uBAM-out is detected", which the hardtrim path satisfies but doesn't honour.

Both are cosmetic-NOTE-placement issues, not data-correctness issues (the actual Phred-0 conversion is correct everywhere via `saturating_sub(33)`). Flagging for completeness.

**NIT-3 (carry-over, unchanged) — `_preserve_tags` field on `BamWriter` is dead.**
`src/bam.rs:544`. Both round-1 reviewers flagged this; it was intentionally left as a documented symmetry marker. No action needed.

### Explicitly NOT re-flagged

- **B-N2** (report filename carries `.bam` extension for BAM input) — confirmed pre-existing/out-of-scope per the review brief; not introduced by this patch.

---

## PLAN §3 / §9 coverage after the patch

- **§3.3 step 3** (ACGTN-only / reject `=` / IUPAC→N on write) — now MET via `validate_and_normalize_seq_for_write` (was the M1 gap).
- **§3.3 step 4** (missing-qual NOTE) — MET for the main `run_ubam_output` path; not emitted on the specialty hardtrim BAM path (NIT-2). PLAN wording arguably wants it there too, but the behaviour it warns about is unchanged, so this is a NOTE-placement nicety.
- **§3.5** (provenance preserved by adding to history) — now MET via `Programs::add()` (was the I1 gap).
- **§9** — both previously-missing rows (`--cores N>1` warn-and-proceed; `--preserve-tags` mixed-batch) now have tests. No remaining §9 row is untested.

No §3 / §9 promise is left silently open.

---

## Verdict

**APPROVE**

Every P0/P1/P2 finding from round 1 is correctly fixed, each headline fix (C1, I1, I2/two-BAM, M1) carries a regression guard that would fail if reverted, and the verification gate is fully green (349 tests, clippy/fmt clean). The fix patch introduces no new correctness defect. The three remaining items are all NITs: the A-I2 read-side bail lacks a direct test (NEW-1), and two cosmetic NOTE-/guard-placement observations (NIT-1, NIT-2). None block the merge; NEW-1 is the only one I'd suggest closing before this path sees heavy hex-tag traffic, and it's a one-test addition.
