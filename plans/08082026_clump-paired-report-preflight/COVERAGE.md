# Plan Coverage Report

**Mode:** B (code vs plan — no separate IMPL.md; ledger built from the plan's Implementation outline, Behavior, and Validation sections)
**Plan(s):** `plans/08082026_clump-paired-report-preflight/PLAN.md` (r2, incl. Implementation notes)
**Date:** 2026-08-08
**Verdict:** INCOMPLETE — 1 item unresolved

Branch `fix/391-clump-report-preflight` @ `2704b52`; base `dev` @ `7ea2741`; PR [#396](https://github.com/FelixKrueger/TrimGalore/pull/396) (OPEN).
Diff scope: `CHANGELOG.md`, `src/io.rs`, `src/main.rs`, `tests/integration_output_collision.rs` — exactly the four files the plan names, and no others.

## Summary

- Total items: 29 (11 Implementation-outline steps + 7 Behavior items + 11 Validation rows)
- DONE: 28
- PARTIAL: 1
- MISSING: 0
- DEVIATED: 0 undocumented (2 deviations documented in the plan's Implementation notes; both verified present and both non-behavioural)

The single unresolved item is one omitted **test assertion**, not a behavioural gap: `clump_se_rejects_report_that_aliases_an_input` does not assert the alias filename in stderr, which outline step 8 specifies. The pre-flight does emit it (verified by hand, below), so the assertion would pass if added.

## Coverage ledger

### Implementation outline (11 steps)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| O1 | `run_specialty_paired`: widen `NameFn` to `-> Vec<PathBuf>`, `planned.extend(pair_outputs(…))`, update doc comment (mode list + "every path the pair writes") | Outline 1 | DONE | `main.rs:2524-2556`. Doc comment now reads "(`--clock`, `--implicon`, `--clump_only --paired`)" and "`pair_outputs` must return EVERY path the pair will write — primaries plus gated secondaries (#391)". Loop body is a single `extend`. |
| O2 | Add `clump_report_candidates` near `planned_secondary_outputs`, matching its doc-comment style | Outline 2 | DONE | `main.rs:114-130`, immediately after `planned_secondary_outputs` (ends `:112`). Signature matches Behavior 5 verbatim: `(no_report_file: bool, report_inputs: &[impl AsRef<Path>], output_dir: Option<&Path>) -> Vec<PathBuf>`. Doc comment carries the same "the pre-flight is only as complete as its candidate list" rationale as its neighbour. |
| O3 | `--clock`/`--implicon` call sites wrap their two names in `vec![…]` | Outline 3 | DONE | `main.rs:485-491` (clock), `:500-506` (implicon). These are the only other `run_specialty_paired` callers (`grep` finds exactly three: `:483`, `:498`, `:547`). |
| O4 | Clump paired-FASTQ arm: namer returns primaries + `clump_report_candidates(cli.no_report_file, &[r1, r2], output_dir)`; hint `None` → `Some(PAIRED_REPORT_HINT)`; replace the comment that justified `None` | Outline 4 | DONE | `main.rs:547-564`. Namer builds `vec![o1, o2]` then extends with the helper over `&[r1, r2]`. Old comment ("clumped_paired_output_names uses input.parent(), not the CWD") replaced by the #391 rationale. |
| O5 | Task 2: four sibling arms extend via the helper with the per-arm inputs; Shape B framed as defensive symmetry in comment and CHANGELOG | Outline 5 | DONE | SE FASTQ `main.rs:589-594` (`&cli.input`); Shape B `:630-637` (`&cli.input`, len==1 guarded at `:618`) with the "defensive symmetry" comment; Shape A `:670-676` (`std::slice::from_ref(&chunk[0])`); SE BAM `:726-731` (`&cli.input`). Each matches its writer's key — see the writer audit below. |
| O6 | Reword `PAIRED_REPORT_HINT` per Behavior 3; rewrite its doc comment (first constant only), naming the user sites | Outline 6 | DONE | `main.rs:55-62`. Copy is byte-identical to Behavior 3's specified text (confirmed against live stderr, below). Doc comment names all three users: "paired trim (FASTQ and uBAM out) and `--clump_only --paired` (#391)" — and `grep` confirms exactly three use sites (`:552`, `:838`, `:1942`). `CWD_OUTPUT_HINT`'s comment (`:64`) left untouched, as specified. No stale "regardless of source directory" text survives anywhere in `src/`, `tests/`, or `docs/`. |
| O7 | io.rs property test(s) extended over the clump namers via `collision_key`; `_clumped_N`/#391 note added to the doc comment | Outline 7 | DONE | `distinct_primary_outputs_imply_distinct_secondary_outputs` (`io.rs:1229+`) now loops three primary namers (`single_end_output_name`, `clumped_output_name`, `clumped_bam_output_name`) with **all** comparisons on `collision_key`, plus a new case-variant input `d/SAME.fastq.gz`. `primary_output_key_is_coarser_than_secondary_keys` (`:1317`) adds all three clump primary namers incl. `clumped_paired_bam_output_name`. Doc comment `:1310-1315` gained "and `_clumped_N` inverts it the same way, #391". See the coverage-split note under Notes. |
| O8 | Nine integration tests in `tests/integration_output_collision.rs`, unique `tempdir` tag each | Outline 8 | **PARTIAL** | All nine exist and pass; tags all unique (44 tags in the file, no tag is a prefix-path of another, so `remove_dir_all` cannot cross-wipe). **Gap:** test 7 (`clump_se_rejects_report_that_aliases_an_input`) omits the specified stderr assertion on the alias filename. Detail in Gaps. |
| O9 | Expected-fail control on every rejection test vs the unpatched build; results recorded in the PR body | Outline 9 | DONE | Plan Implementation notes claim it (all six rejection tests exited 0 unpatched, the two alias cases destroying an input; three acceptance tests passed unpatched) — treated as the evidence per audit instructions. Independently confirmed the record reached PR #396's body verbatim. The six/three split matches the actual test set (6 rejection + 3 acceptance = 9). |
| O10 | CHANGELOG under Unreleased → `#### Bug fixes`: reports join the clump pre-flights; credit the free report-vs-input closure; Shape B as defensive symmetry | Outline 10 | DONE | `CHANGELOG.md:8-18`, under the existing `#### Bug fixes` heading. All three required content points present, including "The single-interleaved-BAM shape gets the same candidate line as defensive symmetry only; with one input its report can never alias it." |
| O11 | `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test` | Outline 11 | DONE | Run live this audit: fmt clean; clippy clean (release, `-D warnings`, no output); `cargo test` 557 passed / 0 failed across 14 binaries. Step 11's note re-verified: no test pins `GENERIC_ADVICE` or any hint text; the only hint fragment pinned by the new tests is `--no_report_file`. |

### Behavior (7 items)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| B1 | Pre-flight surface: per pair, two primaries + (gated) two per-mate reports; list spans all pairs; sibling arms add per-input / `chunk[0]` / single-input candidates | Behavior 1 | DONE | One `planned` vec accumulated across `cli.input.chunks(2)`, pre-flight called once after the loop (`main.rs:2547-2552`). Cross-pair coverage empirically confirmed by test 4. |
| B2 | Newly rejected shapes (a) issue repro, (b) cross-pair fold-equal under `-o`, (c) no-`-o` shared-mate, (d) report-vs-input | Behavior 2 | DONE | (a) test 1, (b) test 4, (c) test 5, (d) tests 7 (SE FASTQ) + 9 (Shape A). All five rejection paths verified passing; (c) reaching the pre-flight also confirms the plan's claim that `validate_paired_input` permits one file as the mate of two different pairs. |
| B3 | Hint fixed at the shared constant, with the specified copy; clump arm switches `None` → `Some`; first doc comment rewritten (not swapped) | Behavior 3 | DONE | Verified in live stderr: the reworded hint renders exactly as the plan specifies, em-dash intact, with correct word gaps across the `\` line continuations. |
| B4 | `--no_report_file` parity — report candidates exist iff the writer writes them, via one shared gate; repro succeeds under the flag | Behavior 4 | DONE | Single early-return gate in the helper (`main.rs:123-125`). All four writers gate on `!no_report_file` (`clump_only.rs:365`, `:534`, `:836`, `:1082`). Tests 2 and 8 are the acceptance controls. |
| B5 | Shared helper with the specified signature, feeding five call sites with the per-arm report inputs | Behavior 5 | DONE | Five call sites verified against the four writers: paired FASTQ `&[r1, r2]` ↔ `clumping_report_name(input_r1/input_r2)`; SE FASTQ `&cli.input` ↔ `(input)`; SE BAM `&cli.input` ↔ `(input)`; Shape A `from_ref(&chunk[0])` ↔ `(inputs[0])`; Shape B `&cli.input` (len 1) ↔ `(inputs[0])`. `grep clumping_report_name src/` finds no fifth production writer (the other two hits are in `clump_only.rs`'s own `mod tests`). |
| B6 | No naming changes — output paths, report paths, and report contents byte-identical on accepted runs | Behavior 6 | DONE | `clump_only.rs` (all writers) untouched by the diff; every `io.rs` hunk falls inside `mod tests` (starts `:569`, hunks at `:1236`/`:1308`/`:1348`), so no naming function changed. `main.rs` changes are confined to candidate lists, the hint constant, and the namer's return type. Empirically pinned by test 3's exact-set assertion and tests 2/3/8's content verification. |
| B7 | Edge cases: odd counts and N=1 `--paired` still rejected earlier; `--basename` reports still key on input filenames; duplicate-mate reuse tested | Behavior 7 | DONE | The N=1 `--clump_only --paired` guard (`main.rs:539-545`) is untouched context in the diff. `--basename` variant is test 6; duplicate-mate reuse is test 5. |

### Validation (11 rows)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| V1 | Issue repro rejected, attributed to the report | Validation 1 | DONE | `clump_paired_rejects_shared_report_name` passes: `assert_rejected_cleanly` (non-zero exit + collision prefix + `DUP_MSG` + empty out dir) plus separate asserts on `reads.fq_clumping_report.txt` and `--no_report_file`. |
| V2 | Gate parity | Validation 2 | DONE | `clump_paired_accepts_shared_report_name_with_no_report_file`: exit 0, both primaries content-verified at 40 reads with source-attributable IDs, `assert_dir_holds_only` proves zero report files. |
| V3 | Candidate list == written set | Validation 3 | DONE | `clump_paired_accepts_distinct_filenames` asserts exactly the four expected filenames (two primaries + two reports) and content-verifies both primaries. |
| V4 | Cross-pair surface under `-o` | Validation 4 | DONE | `clump_paired_rejects_cross_pair_report_collision`: rejection only, no filename pin — matching the plan's reasoning about candidate order. |
| V5 | Case-free no-`-o` class, on ext4 AND APFS | Validation 5 | DONE | `clump_paired_rejects_shared_mate_report_without_output_dir` passes locally (APFS). Fixtures are all-lowercase and the collision is exact-string, so the case dimension is absent by construction and the result cannot differ by filesystem. The ext4 confirmation is CI's `Rust Tests (ubuntu-latest)`, **pending** at audit time (run 31254241237) — an external confirmation in flight, not an implementation gap. |
| V6 | Report-vs-input, SE: rejection + `ALIAS_MSG` + inputs intact; acceptance sibling | Validation 6 | DONE | Both tests pass with everything this row's expectation names. (The assertion the plan asked for beyond this row lives in O8 — see Gaps.) |
| V7 | `chunk[0]` keying, Shape A | Validation 7 | DONE | `clump_paired_bam_rejects_report_that_aliases_an_input` passes: rejection, `ALIAS_MSG`, all four inputs present with nothing extra written (`assert_dir_holds_only`), content-verified for the two a pre-patch run would have clobbered. The discriminator holds: keyed on `chunk[1]`, pair 1 would plan `y.fq_clumping_report.txt`, collide with nothing, and the run would proceed — so a mis-keyed helper call fails this test. |
| V8 | Report path isolated from stems | Validation 8 | DONE | `clump_paired_rejects_shared_report_name_with_basename`: `--basename foo` makes the primaries `foo_clumped_1/2` — non-colliding by construction — and the run is still rejected, naming the report. |
| V9 | A3 held by CI | Validation 9 | DONE | Both io.rs property tests pass (`cargo test --lib io::tests` — 54 passed; second test confirmed individually). Non-vacuity witness: the new `d/SAME.fastq.gz` input fold-collides with `d/same.fastq.gz`, so it is skipped only because the guard is now `collision_key`-based — under the previous `PathBuf` guard that pair would reach the assertion and fail it. The metric change is therefore load-bearing, not cosmetic. |
| V10 | Fix actually fires (expected-fail control) | Validation 10 | DONE | Per audit instructions, the plan's Implementation notes are the evidence; they claim all six rejection tests exited 0 against the unpatched tree. Independently confirmed the same record appears in PR #396's body. |
| V11 | No collateral | Validation 11 | DONE | Live: 557 tests green (0 failed), fmt clean, clippy `-D warnings` clean. The three existing tests the plan singled out all exist and are in the green set: `pe_byte_identity` (`tests/integration_clump_only.rs:142`), `multi_pair_pe_bam_produces_one_output_per_pair` (`tests/integration_clump_only_ubam.rs:477`), `pe_bam_collision_preflight_case_folded` (`:630`). |

### Documented deviations (Implementation notes)

Both are recorded in the plan and verified present; neither changes behaviour, so per the audit rules neither is a gap.

| Deviation | Verified |
|---|---|
| `assert_dir_holds_only`'s failure message generalized to "directory must hold exactly the expected files", since it now also serves acceptance-side exact-set assertions | Yes — `tests/integration_output_collision.rs:86-89`, with the doc comment updated to describe both uses |
| `distinct_primary_outputs_imply_distinct_secondary_outputs` restructured to loop three primary namers; the demux-stem sub-check stays scoped to the SE trim namer and is now guarded by `collision_key` inequality | Yes — `src/io.rs:1247-1301`; sub-check reads `let se_primaries = &primary_sets[0];` and guards on `collision_key(a) != collision_key(b) && a.parent() == b.parent()` |

## Gaps (detail)

### Item O8: `clump_se_rejects_report_that_aliases_an_input` omits the specified alias-filename assertion

**Expected** (outline step 8, seventh bullet): "assert `ALIAS_MSG` **+ the alias filename** (the input branch names one path only), `count_reads_from`-verified intact inputs, `assert_dir_holds_only(&dir, &[both inputs])`."

**Found** (`tests/integration_output_collision.rs:1206-1221`): the test asserts non-zero exit, `stderr.contains(ALIAS_MSG)`, `count_reads_from` on the alias file (40 reads, `REPORTY`), and `assert_dir_holds_only(&dir, &["s.fq", "s.fq_clumping_report.txt"])`. There is no assertion that stderr contains `s.fq_clumping_report.txt`.

**Gap:** one assertion. The pre-flight's output-vs-input branch does print the path — verified by hand against the built binary:

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety): this run would
write output to …/se/s.fq_clumping_report.txt, which is also one of its inputs. …
```

So adding `assert!(stderr.contains("s.fq_clumping_report.txt"), …)` would pass as written. Without it, the test would still pass if the rejection came from some other alias in the run rather than the clumping-report candidate — which is precisely the attribution the plan wanted pinned, and which its sibling tests (1, 5, 6) do pin on their own messages.

Two related observations, recorded for completeness rather than as gaps, since the plan's own Validation rows do not ask for them:

- Content verification in tests 7 and 9 covers the inputs a pre-patch run destroyed (the alias files) rather than literally all inputs; the remaining inputs are covered for existence and no-extra-writes by `assert_dir_holds_only`. Outline step 8's wording ("intact inputs" / "all four inputs intact") is satisfied in substance.
- Outline step 7 permits "and/or", and all three clump primary namers are covered across the two property tests — but the split is asymmetric: `clumped_output_name` and `clumped_bam_output_name` are checked in **both** tests, while `clumped_paired_bam_output_name` appears only in `primary_output_key_is_coarser_than_secondary_keys`. The paired-BAM namer is therefore checked in the "spellings collapse to one primary" direction but not in `distinct_primary_outputs_imply_distinct_secondary_outputs`' "distinct primaries ⇒ distinct reports, on `collision_key`" direction. Adding it there would be a one-line `map` (`|p| clumped_paired_bam_output_name(p, None, output_dir, basename)`), matching how test 2 already invokes it.

## Test verification (Mode B)

`cargo test --test integration_output_collision clump` — 9 passed, 0 failed, 34 filtered out.

| Test name | File | Status |
|-----------|------|--------|
| `clump_paired_rejects_shared_report_name` | `tests/integration_output_collision.rs:1016` | PASS |
| `clump_paired_accepts_shared_report_name_with_no_report_file` | `:1048` | PASS |
| `clump_paired_accepts_distinct_filenames` | `:1077` | PASS |
| `clump_paired_rejects_cross_pair_report_collision` | `:1112` | PASS |
| `clump_paired_rejects_shared_mate_report_without_output_dir` | `:1141` | PASS |
| `clump_paired_rejects_shared_report_name_with_basename` | `:1175` | PASS |
| `clump_se_rejects_report_that_aliases_an_input` | `:1206` | PASS (assertion set short of spec — see Gaps) |
| `clump_se_accepts_report_alias_with_no_report_file` | `:1225` | PASS |
| `clump_paired_bam_rejects_report_that_aliases_an_input` | `:1260` | PASS |

`cargo test --lib io::tests` — 54 passed, 0 failed.

| Test name | File | Status |
|-----------|------|--------|
| `distinct_primary_outputs_imply_distinct_secondary_outputs` | `src/io.rs:1229` | PASS |
| `primary_output_key_is_coarser_than_secondary_keys` | `src/io.rs:1317` | PASS (confirmed individually as well) |

Whole-suite and gate checks, run live during this audit:

| Check | Result |
|---|---|
| `cargo test` (all 14 test binaries) | 557 passed, 0 failed, 0 ignored |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean |
| PR #396 CI | Lint, Coverage, Docs, Reproducibility, Security audit pass; `Rust Tests (ubuntu-latest)` and `(macos-latest)` pending at audit time |

## Verdict

**INCOMPLETE — 1 item unresolved.**

Everything the plan specifies is implemented and verified except one test assertion:

1. **`tests/integration_output_collision.rs`, `clump_se_rejects_report_that_aliases_an_input`** — add the stderr assertion on the alias filename that outline step 8 specifies:

   ```rust
   assert!(
       stderr.contains("s.fq_clumping_report.txt"),
       "the aliased path must be the clumping report:\n{stderr}"
   );
   ```

   The message already contains that path (verified against the built binary), so this is an added assertion with no production change. Either add it, or record the omission in the plan's Implementation notes as a deliberate deviation — the plan's other five rejection tests all pin their message this way, so the asymmetry is currently undocumented.

Nothing else remains. No plan item is MISSING; the two deviations in the Implementation notes are both present in the code and both non-behavioural; the ext4 leg of Validation 5 is an in-flight CI job on a test whose construction is case-free, not an implementation gap.

---

## Post-audit resolution (2026-08-08, caller)

The single unresolved item (O8: missing alias-filename assertion in
`clump_se_rejects_report_that_aliases_an_input`) was implemented in commit
`d99f23e` exactly as this report specified, and passes. The same commit applied
the review batch Felix approved (fold-equal + over-rejection tests, cross-pair
class assertion, `clumped_paired_bam_output_name` in the distinct-primaries
property test — closing this report's second recorded observation — plus
CHANGELOG/docs/comment precision items). Full suite after the batch: 559
passed, 0 failed; fmt and clippy `-D warnings` clean. With the gap closed, the
plan's ledger stands fully covered.
