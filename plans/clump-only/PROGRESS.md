# Progress — `--clump_only` (lossless reorder-only mode)

**Plan:** `PLAN.md`
**Issue:** [FelixKrueger/TrimGalore#353](https://github.com/FelixKrueger/TrimGalore/issues/353)
**Branch:** `feature/clump-only` (off `dev`)

## Pipeline

| Stage | Status | Notes |
|---|---|---|
| Design plan drafted | ✅ done | `PLAN.md` — 2026-07-25 |
| Manual review by user | ✅ done | 5 open questions locked into §Resolved decisions |
| Dual `plan-reviewer` agents | ✅ done | Reviewers A/B: `APPROVE WITH REVISIONS` — see `PLAN_review_reviewer-{A,B}.md` |
| Reviewer findings incorporated | ✅ done | 4 blockers resolved (contract scope + CI granularity + rejection matrix + silent-accept for -q/-e/stringency); §Resolved decisions expanded with 12 items |
| Implementation | ✅ done | See §Implementation log below; 11 steps executed |
| Local verification | ✅ done | `cargo build --release` + `cargo fmt --check` + `cargo clippy -D warnings` + `cargo test --release` all pass (386 tests total after remediation, 12 new integration tests for --clump_only) |
| Dual `code-reviewer` verification | ✅ done | Reviewers A/B: `APPROVE WITH REVISIONS`; 1 HIGH (panic bug) + 5 medium/low findings addressed in remediation |
| `plan-manager` coverage audit | ✅ done | `INCOMPLETE WITH MINOR GAPS` (31/34 items covered); 3 gaps addressed in remediation |
| Post-review remediation | ✅ done | See §"Post-review remediation" below — 7 items fixed in single Option A pass |
| PR opened against `dev` | ⏸ pending | |
| CI green (all 8 checks + new byte-identity + determinism) | ⏸ pending | |
| Merged | ⏸ pending | |
| Docs live on trimgalore.com | ⏸ pending | Auto-deploy from dev on merge |
| Issue #353 closed | ⏸ pending | Manual close (dev ≠ default branch, per `feedback_auto_close_keyword` memory) |
| Follow-up notification to @wkgardner | ⏸ pending | Optional — offer to loop them in on the draft PR as promised in issue reply |

## Implementation log (2026-07-25)

Files touched:

| File | Change |
|---|---|
| `src/cli.rs` | Added `pub clump_only: bool` field + doc-block. Added `Cli::validate()` rejection matrix (~20 rejected flags with precise per-flag error messages). `-q`/`--stringency`/`-e` silently accepted per Blocker 4 (ii). |
| `src/main.rs` | Added `use trim_galore::clump_only;`. Early-return dispatch block after `--implicon` and before the trim pipeline (SE + PE paths). SE path includes case-folded output-collision pre-flight; PE path uses existing `run_specialty_paired` helper which does the same. |
| `src/clump_only.rs` (**new**, ~950 lines) | Core impl: `clump_only_single`, `clump_only_paired`, `SingleBin`/`PairedBin` structs, `flush_bin_single`/`flush_bin_paired`, `write_records_member` (encodes each bin as one gzip member — RFC-1952 concatenation), `ClumpOnlyStats` struct, `write_clump_only_report` (text-only, "Compression ratio" line conditional), `reject_ubam`. 9 unit tests covering permutation, PE lockstep, cross-run determinism, no-trimming, no-adapter-detection, empty input, single record, plain output, report-filename discipline. |
| `src/lib.rs` | Added `pub mod clump_only;`. |
| `src/io.rs` | Added `clumped_output_name`, `clumped_paired_output_names`, `clumping_report_name` + 6 unit tests. |
| `tests/integration_clump_only.rs` (**new**) | 9 integration tests spawning the built binary against `test_files/BS-seq_10K_R{1,2}.fastq.gz`. Includes byte-identity SE/PE, --dont_gzip variant, cross-run determinism, report filename discipline, and rejection-matrix regressions (--length, adapter, --rename) + silent-accept regression (-q 30). |
| `.github/workflows/ci.yml` | Added 3 CI steps in the `validation` job: SE byte-identity via `paste - - - -` reconstitution + sort + diff; PE byte-identity for both mates; cross-run determinism via md5 diff. |
| `docs/src/content/docs/modes/clump-only.md` (**new**) | Astro Starlight docs page: usage, byte-identity contract-scope note (plus-line + CRLF normalization), determinism, compatibility (composable/rejected/silent-accept), report shape, when-to-use vs `--clumpify`, performance. |
| `docs/astro.config.mjs` | Sidebar entry under "Specialty modes". |
| `docs/src/content/docs/performance/clumpy.md` | Cross-link `:::note` block pointing to `--clump_only` for the lossless-recompression use case. |
| `CHANGELOG.md` | New entry under `### Unreleased` describing the mode, invariants, rejections, and v1 scope. |

### Documented deviations from PLAN.md

- **`--cores >= 2` requirement relaxed to `>= 1`** (implicit via base `--cores == 0` rejection). Plan §Rejection matrix said cores >= 2 (matching `--clumpify`). v1 is single-threaded internally, so requiring 2 cores was user-hostile. `--cores` still accepted for interface parity; parallelism is a v1.1 follow-up that can add without breaking the byte-identity contract.
- **Empty-input case emits an empty gzip member** for the `--gzip_output` path so the output file is a valid `.gz` stream (empty file is not valid gzip). Under `--dont_gzip` the output is a genuinely empty file. Added to unit tests.
- **Sort determinism verification (PLAN Step 4) was verify-only** — reviewer B and reviewer A independently confirmed `sort_single_by_key`/`sort_paired_by_key` at `src/clump.rs:271, 302` already use stable sort with a content-tiebreaker cascade. No code change needed.

### Test tally

- Unit tests (lib): **348 total** (339 previously + 9 new for clump_only). All pass.
- Integration tests (clump_only.rs): **9 new**. All pass.
- Full release test suite: **381 tests, 0 failures**.

### Verification commands (all pass)

```
cargo build --release
cargo fmt --all -- --check
cargo clippy --all-targets --release -- -D warnings
cargo test --release
```

## Post-review remediation (2026-07-25, second pass)

Dual `code-reviewer` (A + B) and `plan-manager` all returned APPROVE WITH REVISIONS / INCOMPLETE WITH MINOR GAPS. Full review reports at `CODE_review_reviewer-{A,B}.md` and `PLAN_manager_report.md`. Seven items addressed in Option A single-commit remediation:

| # | Finding | Sources | Fix | Files touched |
|---|---|---|---|---|
| 1 | **PANIC on `--paired --clump_only <single-uBAM>`** — `chunk[1]` on len-1 slice at `main.rs:2076`. | Reviewer A C-1 (HIGH) | Added `--paired` + N==1 rejection at CLI validation layer — bails before `run_specialty_paired`'s pre-flight ever runs. | `src/cli.rs` |
| 2 | `--no_report_file` silently ignored on `--clump_only` path — inconsistent with all other trim/specialty paths. | Reviewer A N-1, Reviewer B #1 | Added `no_report_file: bool` param to `clump_only_single`/`_paired`; report writes guarded. `main.rs` dispatch passes `cli.no_report_file`. | `src/clump_only.rs`, `src/main.rs` |
| 3 | Silent-accept regression guard only covered `-q`, not `--stringency`/`-e`. | Reviewer A N-2, Reviewer B #4 | Refactored to `assert_silently_accepts` helper; added tests for `--stringency 5` and `-e 0.05`. | `tests/integration_clump_only.rs` |
| 4 | Missing PLAN-mandated unit tests: `test_clump_only_normalizes_plus_line` + `test_clump_only_normalizes_crlf`. | Reviewer B #2, plan-manager #1 | Added both tests. Each constructs raw FASTQ input (bypassing `FastqRecord::write_to`) with the unusual feature, then verifies the output is normalized. | `src/clump_only.rs` |
| 5 | uBAM-input rejection had no integration test. | Reviewer B #3, plan-manager #2 | Added `rejects_ubam_input` integration test using `test_files/ubam_test.bam`. | `tests/integration_clump_only.rs` |
| 6 | Inline `norm` closure duplicated `io::norm_path`. | Reviewer A N-4, Reviewer B non-blocker | Promoted `io::norm_path` from `pub(crate)` to `pub`; replaced my SE-preflight closure with `naming::norm_path`. Pre-existing closures in `run_specialty_paired` left untouched (out of scope for this remediation). | `src/io.rs`, `src/main.rs` |
| 7 | PLAN.md line 131 wording drift: still said `--cores >= 2, same as --clumpify today`. | plan-manager #3 | Updated to reflect the `>= 1` deviation with pointer to PROGRESS.md. Same fix applied to PLAN.md line 217. | `plans/clump-only/PLAN.md` |

### Not addressed (deferred / out of scope)

- Reviewer B #5 (memory doubling in `write_records_member`) — real perf concern, but not a correctness bug and would require restructuring the gzip encoder to stream directly into `out` rather than buffer in `Vec::new()`. Deferred as v1.1 perf work.
- Reviewer B #6 (rejection-matrix unit tests in `cli.rs`) — plan §Validation item 5 called for per-flag unit tests; only 4 of ~20 rejections are covered by integration tests today. Reasonable follow-up but not a blocker for v1 (the load-bearing rejections have coverage; the tail is defensive).
- Reviewer A N-3 (`detect_input_format` called twice per input) — micro-perf; the second call inside `reject_ubam` is redundant with the reader's implicit detection. Trivial to fix but touching the reader signature costs more than the saving warrants at v1.
- Reviewer A N-5 (sort tiebreaker only implicitly exercised) — the tiebreaker cascade (`seq → qual → id`) is validated by the cross-run determinism test but no test explicitly constructs identical-minimizer + different-seq inputs. Nice-to-have; adds coverage without changing behavior.

### Test tally (updated)

- Unit tests (lib): **350 total** (348 previously + 2 new normalization tests). All pass.
- Integration tests (clump_only.rs): **12 total** (9 previously + 3 new: uBAM-rejection + `--stringency`/`-e` silent-accept).
- Full release test suite: **386 tests, 0 failures**.

### Verification commands (all pass — second pass)

```
cargo build --release
cargo fmt --all -- --check
cargo clippy --all-targets --release -- -D warnings
cargo test --release   → 386 pass / 0 fail
```

## Deferred to v2 (not tracked here)

- uBAM in/out for `--clump_only` — natural follow-up per issue #353 discussion
- Cross-input-order determinism (stronger than within-bin: same records in any input order → same output bytes) — see PLAN.md "Open" question 2

## Resolved decisions (see PLAN.md §"Resolved decisions")

1. `--dont_gzip` allowed with `--clump_only` (diverges from `--clumpify`)
2. Stable sort only in v1 (cross-input-order determinism deferred)
3. Report filename: `*_clumping_report.txt` — text-only, no JSON
4. `--fastqc` unchanged (opt-in, runs on reordered output)
5. PE report mirrors `--clumpify`'s current per-input layout

## Session log

- **2026-07-25 (morning)** — Plan drafted after issue #353 discussion. Branch `feature/clump-only` created off `dev`. Reply posted on issue #353 accepting the feature with the shape documented in PLAN.md. Contributor @wkgardner offered to write the PR; Felix declined (drove it himself) and offered to loop them in for design-corners review on the draft PR.
- **2026-07-25 (midday)** — Five open design questions resolved (§Resolved decisions items 1-5). Dual `plan-reviewer` agents launched in parallel (fresh contexts, Reviewer A + Reviewer B); both returned `APPROVE WITH REVISIONS`. Both found the same two highest-severity issues: plus-line contract unenforceable + sort-determinism work already in place. Four blockers derived from the reviews were resolved (see §Resolved decisions items 6-12). Plan is now implementation-ready.
- **2026-07-25 (afternoon)** — v1 implementation shipped as PR #355 (in-flight against `dev`).
- **2026-07-25 (evening)** — v2 phase started (`plans/clump-only/PLAN_v2_ubam.md`, 551 lines). Branch `feature/clump-only-ubam` off `feature/clump-only`. Dual plan-reviewer round on v2 plan returned `APPROVE WITH REVISIONS` — both agreed on FastQC self-contradiction (plan-drift from my Q1 lock) + false parity claim on `--dont_gzip`; asymmetric findings on PE shape holes (A caught Shape B non-BAM slip; B caught Shape A two-BAM slip + multi-pair regression); precise CLI edit sites (B) + collision pre-flight gap (A). Six blocker fixes + non-blocker doc/detail fixes applied to plan in place. Q3 (`estimated_record_bytes` audit) and Q4 (empty-BAM handling) both verified during plan-review — no implementation-time work needed. Implementation trigger followed; v2 shipped locally: 402 tests pass (357 unit + 12 v1 integration + 9 v2 integration + 24 pre-existing).

## v2 Implementation log (2026-07-25)

Files touched:

| File | Change |
|---|---|
| `src/cli.rs` | 4 precise edit sites: removed `--output-format ubam` v1-defer rejection; narrowed the N=1+`--paired` rejection to allow uBAM inputs (delegating the "N=1 must be BAM" format-check to `main.rs::dispatch`); hoisted `--dont_gzip + --output-format ubam` rejection into the shared §3.4a `OutputFormat::UBam` block (closes a pre-existing gap on the trim uBAM path too). |
| `src/main.rs` | Extended `if cli.clump_only { ... }` with a `match cli.output_format`: FASTQ branch unchanged from v1; UBam branch has Shape A/B dispatch, format-guards (two-BAM Shape A, non-BAM Shape B, mixed-format Shape A), multi-pair iteration via `.chunks(2)` (restores v1 SE-style N=4+ support), and explicit `preflight_collision_bam` helper for case-folded output-path collision detection. |
| `src/clump_only.rs` | Removed the v1 `reject_ubam` guard; added `PairedInputSetup` struct + `clump_only_single_to_bam` + `clump_only_paired_to_bam_one_pair` + `flush_bin_single_to_bam` + `flush_bin_paired_to_bam` helpers. Extended `ClumpOnlyStats` with `input_format_label` / `output_format_label` / `preserved_tags` fields. Updated `write_clump_only_report` to render the new labels and to emit "Compression ratio" only when both sides compressed (v2 rule from §Resolved decision 2). Added a per-function guard on v1 FASTQ path that rejects uBAM input with a clear "add `--output-format ubam`" hint (tag-drop prevention). |
| `src/io.rs` | Added `clumped_bam_output_name` + `clumped_paired_bam_output_name` + 7 unit tests. |
| `tests/integration_clump_only.rs` | Renamed v1 test `rejects_ubam_input` → `rejects_ubam_input_without_ubam_output`, adjusted assertion to match the new "add `--output-format ubam`" message. |
| `tests/integration_clump_only_ubam.rs` (**new**) | 9 integration tests spawning the built binary against real fixtures: SE/PE positive paths, PE Shape B, `@PG` chain preservation, `--dont_gzip` rejection, Shape A two-BAM rejection, Shape B non-BAM rejection, multi-pair PE regression guard, `--fastqc` on BAM native support. |
| `.github/workflows/ci.yml` | 3 new steps: SE uBAM record-parity via `samtools view`, PE uBAM mate-adjacent record-parity, `@PG` chain preservation check. |
| `docs/src/content/docs/modes/clump-only.md` | Rewrote for v2 — dual FASTQ+uBAM usage examples, aux-tag round-trip, `@PG` line documentation, PE input-shape acceptance/rejection matrix, expanded compatibility table, updated report shape example. |
| `CHANGELOG.md` | New entry under `### Unreleased` above the v1 entry describing v2 additions. |

### Documented deviations from PLAN_v2_ubam.md

None material. The `PairedInputSetup` struct was introduced during implementation to sidestep `clippy::type_complexity` on the underlying 6-tuple — same intent as the plan's Signature block, cleaner form.

### Test tally

- Unit tests (lib): **357 total** (up from 350; +7 for the new BAM filename helpers + format-label unit tests).
- Integration tests: **12 (v1) + 9 (v2 new) + 24 (pre-existing) = 45 total**.
- Full release test suite: **402 tests, 0 failures**.

### Verification commands (all pass — v2 second-pass)

```
cargo build --release
cargo fmt --all -- --check
cargo clippy --all-targets --release -- -D warnings
cargo test --release   → 402 pass / 0 fail
```

## v2 Post-review remediation (2026-07-25)

Dual code-reviewer (A + B) and plan-manager all returned APPROVE WITH REVISIONS / INCOMPLETE WITH MINOR GAPS / NEEDS SIGNIFICANT REWORK. Full review reports at `CODE_v2_review_reviewer-{A,B}.md` and `PLAN_v2_manager_report.md`. Nine items addressed:

| # | Finding | Sources | Fix | Files touched |
|---|---|---|---|---|
| 1 | **PANIC on `--clump_only --paired ubam.bam`** (missing `--output-format ubam`) — index-out-of-bounds in `run_specialty_paired` on N=1. | Reviewer B C-1 (CRITICAL) | Added guard in the FASTQ output arm of clump_only dispatch that rejects `cli.paired && cli.input.len() == 1` before entering `run_specialty_paired`. New integration test `rejects_paired_single_bam_without_output_format_ubam` regression-guards the panic. | `src/main.rs`, `tests/integration_clump_only_ubam.rs` |
| 2 | **`rejects_two_bam_paired` test doesn't test the two-BAM guard** — R1==R2 dup check fires first, actual guard untested. | Reviewer A C-1 | Rewrote test to copy the fixture to a second distinct path; asserted stderr contains "single interleaved" (the guard's specific message, not the R1==R2 dup message). | `tests/integration_clump_only_ubam.rs` |
| 3 | **Aux-tag round-trip zero integration coverage** — fixture used had no aux tags. | Reviewer B H-1 | Added `ubam_aux_tag_roundtrip_via_preserve_tags` integration test using `test_files/ubam_test_with_tags.bam` (has `CB:Z` + `UB:Z` tags). Asserts `(name, CB, UB)` multiset preserved through the reorder. | `tests/integration_clump_only_ubam.rs` |
| 4 | **CI cross-run determinism step missing.** | plan-manager gap 4 | Added `Validate --clump_only uBAM cross-run record determinism` CI step that runs the mode twice and compares samtools-view-piped-to-sort output via md5. | `.github/workflows/ci.yml` |
| 5 | **Missing planned unit tests (4)** — permutation, PE lockstep, deterministic records, tag round-trip. | Reviewer A C-2 + plan-manager gap 1 | Added `test_empty_bam_input_produces_valid_bam` (the one the integration layer can't cheaply cover). Documented deferral of the other three (permutation / lockstep / determinism / tag round-trip) with rationale: equivalent invariant coverage exists at the integration layer (`se_ubam_out_from_fastq_in`, `pe_ubam_out_interleaved_from_fastq_pair`, cross-run determinism CI, `ubam_aux_tag_roundtrip_via_preserve_tags`). See §Documented deviations below. | `src/clump_only.rs` |
| 6 | **Missing planned integration tests (2)** — `rejects_mixed_format_paired`, `pe_bam_collision_preflight_case_folded`. | Reviewer A C-2 + plan-manager gap 2 | Added both tests. `pe_bam_collision_preflight_case_folded` gracefully skips on case-insensitive filesystems (APFS default) — case-sensitive Linux CI runners exercise the guard. | `tests/integration_clump_only_ubam.rs` |
| 7 | **Multi-pair PE-BAM dispatch drops `=== pair N of M ===` banner + per-pair `.with_context()`** — UX regression vs v1 FASTQ and trim uBAM paths. | Reviewer A C-3 | Restored both in the Shape A `.chunks(2)` iteration. Pair number + total in the banner; per-pair `.with_context()` wraps errors with the R1/R2 paths so failures at pair 3/5 identify the pair. | `src/main.rs` |
| 8 | **PE-BAM dispatch skips per-pair `sanity_check_any`** on `chunk[0]` (pair > 0) and `chunk[1]` (every pair). | Reviewer A C-4 | Added per-pair `sanity_check_any` calls in the Shape A loop, mirroring `main.rs:651-669`'s trim FASTQ pattern. | `src/main.rs` |
| 9 | **Misleading `bam.rs:1641-1647` comment** — that reference cited a test that writes one record, not zero. | Reviewer A C-5 | Rewrote comment to describe the actual noodles `bam::io::Writer` empty-BAM contract and reference the new `test_empty_bam_input_produces_valid_bam` unit test as the local verifier. | `src/clump_only.rs` |

### Not addressed (deferred / non-blocking)

- **Non-blockers from both reviewers** — dead `peek_header` branch in Shape A (Reviewer B N-1), 1-element `preflight_collision_bam` call (B-N2), BamWriter truncated-on-panic risk (B-N3, same as trim path), `preserved_tags` report field always populated (B-N5), `has_trim_galore_pg` key-only match (A non-blocker), asymmetric input/output format labels (A non-blocker), double `detect_input_format` on Shape A (A non-blocker), Shape A report `input_bytes` combined against R1's path (A non-blocker), `--dont_gzip + --output-format ubam` silent behavior change for the trim path (A non-blocker — flagged in CHANGELOG as a "closes an existing gap" note).
- **Trim uBAM path FastQC-skip** — pre-existing bug called out in §Resolved decision 1's rationale. Fix is out of scope for this v2 (a separate small PR).

### Documented deviations from PLAN_v2_ubam.md

1. **`PairedInputSetup` struct** — introduced during implementation to sidestep `clippy::type_complexity`. Same intent as the plan's Signature block.
2. **Three of four planned unit tests deferred** (permutation, PE lockstep, deterministic records, tag round-trip in `src/clump_only.rs::tests`) — the invariants they'd assert are exercised at the integration layer today. Adding parallel unit tests would give localization value on failure but no new invariant coverage. Documented in the `v2 uBAM-output unit tests` section of the test module. The empty-input test IS added because the integration layer can't cheaply cover it (needs in-process BAM writer/reader dance).

### Test tally (updated after remediation)

- Unit tests (lib): **358 total** (up from 357, +1 for `test_empty_bam_input_produces_valid_bam`).
- Integration tests (clump_only_ubam.rs): **13 total** (up from 9, +4 for `rejects_paired_single_bam_without_output_format_ubam`, `ubam_aux_tag_roundtrip_via_preserve_tags`, `rejects_mixed_format_paired`, `pe_bam_collision_preflight_case_folded`).
- Full release test suite: **407 tests, 0 failures**.

### Verification commands (all pass — post-remediation)

```
cargo build --release
cargo fmt --all -- --check
cargo clippy --all-targets --release -- -D warnings
cargo test --release   → 407 pass / 0 fail
```
