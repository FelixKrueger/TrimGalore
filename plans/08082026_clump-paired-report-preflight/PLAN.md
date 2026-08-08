# Plan: `--clump_only --paired` clumping reports join the collision pre-flight (#391)

**Issue:** [#391](https://github.com/FelixKrueger/TrimGalore/issues/391) — filed from #388's code review (both reviewers, independently; reproduced on the #388 branch).
**Base:** `dev` @ `7ea2741`.

## Revision history

- **r2 (2026-08-08):** incorporates dual plan-review (PLAN_REVIEW_A.md "approve with changes", PLAN_REVIEW_B.md "no criticals") and two maintainer decisions: **fix the shared `PAIRED_REPORT_HINT` constant** (not a new clump-specific one) and **use a shared `clump_report_candidates` helper** (not five inline blocks). Corrects four r1 prose errors (Shape B exposure over-claim A§1.5; duplicate-inputs over-claim A§1.7; "swap the comments" instruction A§1.6; `--basename` mechanism B-A4). Adds: the A3 property test (A§4.2/B-V1, both reviews' highest-value item), the Shape A representative test (A-Critical-1/B-V3), the no-`-o` cross-pair same-mate test (A§4.3), the test-5 respec (A§4.4/B-V2), the exact-set upgrade of test 3 (B-V4), the `--basename` variant of test 1 (A§4.4), expected-fail controls on every rejection test (A§4.4), unique tempdir tags (B-V5b), and the corrected A4/A5 assumptions.
- r1 (2026-08-08): initial plan; Task 2 scope confirmed IN by Felix.

## Goal

`--clump_only --paired -o out A/reads.fq B/reads.fq` exits 0 with **one** clumping report: primaries carry positional discriminators (`reads_clumped_1.fq` / `reads_clumped_2.fq`) so they never collide, but `clumping_report_name` keys on the bare input filename, so both mates plan `out/reads.fq_clumping_report.txt` and R2's report silently overwrites R1's. Make the per-mate clumping-report names part of the collision pre-flight — gated on `!no_report_file` to match the writer — so this fails loudly before any I/O. Task 1 thereby also closes report-vs-**input** on the paired arm for free (the pre-flight runs its input check on every candidate — A§1.8). Task 2 extends the same candidates to the sibling clump arms, whose real exposure is the report-vs-input direction (#385 defect-2 class).

## Context

- **Dispatch:** `main.rs:526-550` — the `--clump_only --paired` FASTQ arm calls `run_specialty_paired` with a namer returning only the two primaries.
- **Driver:** `main.rs:2476-2518` — the namer feeds only the pre-flight (called once, at `main.rs:2490`); written paths are re-derived by the writers (`clump_only.rs:415-416`, `535-536`). Verified three ways by review A§1.2: widening the namer's return type cannot change a written byte. Doc comment's mode list (`main.rs:2470-2471`) still reads "(`--clock`, `--implicon`)" and must gain `--clump_only --paired`.
- **Writers and gates (all four, verified A§1.3):** SE FASTQ `clump_only.rs:365` (keyed per input); paired FASTQ `:534` (both mates in one gate, keyed on `input_r1`/`input_r2`); SE BAM `:836` (per input); paired BAM Shapes A+B `:1082` (ONE report per pair, keyed on `inputs[0]` — `clump_only.rs:1083`). `grep clumping_report_name src/` confirms no fifth writer; `--clumpify` on the trim path writes no clumping report.
- **Namer:** `io.rs:469-482` — full input filename + `_clumping_report.txt`, into `output_dir` or **that input's own** `parent()`. Ignores `--basename`.
- **Directory asymmetry (review finding, load-bearing):** without `-o`, the paired arm's *primaries* both go to `input_r1.parent()` (`io.rs:407-409`) while each *report* follows its own mate's parent (`io.rs:478-481`). Consequences: the hint must not claim directory-independence (see Behavior 3), and a collision class exists with no `-o` at all (see Behavior 2c). Same asymmetry exists at #388's trim site — follow-up material, not fixed here.
- **Pre-flight:** `io.rs:96-130`. Output-vs-output errors name both paths + hint; the output-vs-**input** branch names only the input path (`io.rs:109-115`) — constrains test assertions (r1's A5 was half-wrong).
- **`--clock`/`--implicon`:** verified unaffected — `specialty.rs` has zero report/fastqc writes.
- **uBAM arms and report-vs-report (A3):** re-derived independently by both reviewers and airtight — report and primary key off the *same* input with the *same* directory rule, and `strip_fastq_extensions` is fold-equivariant (`strip_suffix_ignore_ascii_case` + positional `file_stem`), so a report collision implies a primary collision the pre-flight already rejects. Premises now stated explicitly: `gzip` is global, derived from `cli.input[0]` (`main.rs:351`); multi-input `--basename` never reaches these arms — rejected at `cli.rs:624` (SE) / `cli.rs:651` (paired), which is the *actual* guard (r1 cited a nonexistent primary-collision mechanism).
- **Shape B is not exposed at all** (r1 over-claimed): it requires `input.len() == 1`, and `--passthrough`/`--demux` are rejected under `--output-format ubam` (`cli.rs:588/602`), so `guarded_inputs` is exactly the one input and the report path is that filename plus a suffix — never fold-equal to it. Task 2's Shape B line is **defensive symmetry**, and the CHANGELOG/PR must say so, not claim a closed hole.
- **Duplicate inputs (r1 over-claimed):** `validate_paired_input` rejects only within-pair R1==R2 and *exact* duplicate pairs; the general duplicate scan (`cli.rs:631`) is SE-only. **A file reused as a mate of two different pairs passes validate** — which enables the case-free collision in Behavior 2c.
- **Test home:** `tests/integration_output_collision.rs` (question closed in r1's favour — B§5.5: its module doc narrates this defect family and its helpers/constants are exactly what these tests need). Shared constants `DUP_MSG`/`ALIAS_MSG`/`PREFIX` at `:89-91`; models: `assert_rejected_cleanly` `:95`, `assert_dir_holds_only` `:77`, `se_trim_rejects_output_that_aliases_a_report_input` `:674`.

## Behavior

1. **Pre-flight surface.** Candidates per paired-FASTQ pair: the two primaries plus — when `!cli.no_report_file` — the two per-mate clumping-report paths; the list spans all pairs. The four sibling arms add their report-keyed candidates the same way (per input for SE FASTQ/BAM; `chunk[0]` per pair for Shape A; the single input for Shape B).
2. **Newly rejected shapes** (all exit 0 with silent loss today):
   a. the issue's repro — same filename mates, `-o`, reports collide;
   b. cross-pair fold-equal report names under `-o`;
   c. **no `-o`, same file as the mate of two pairs** (A§4.3): reports collide in the shared mate's own parent while all four primaries stay distinct — case-free and filesystem-independent;
   d. report-vs-input: an input named like a clumping report is no longer silently overwritten (paired arm via Task 1; SE FASTQ/BAM and Shape A via Task 2).
3. **Hint (maintainer decision: fix the shared constant).** `PAIRED_REPORT_HINT`'s "collide **regardless of source directory**" is false wherever reports follow their inputs' parents (this arm and #388's site alike), and "named from the input filename alone" points at the primaries, which are fine. Reword the constant once, accurately for all its sites (trim paired FASTQ `main.rs:784`, trim→uBAM `main.rs:1888`, and now this arm):
   > "Outputs and reports are named from the input filename alone, so inputs sharing a filename can collide when `--output_dir` (or a shared input directory) sends them to one place — rename one input, or pass `--no_report_file` if only the reports collide."
   The clump arm's hint argument changes `None` → `Some(PAIRED_REPORT_HINT)`. **Doc-comment fix is a rewrite of the FIRST comment only** (`main.rs:55-56`) — r1's "swap them" would make it self-referential (A§1.6); the second comment is already correct. The rewritten comment names all users.
4. **`--no_report_file` parity.** Report candidates exist iff the writer will write them — one shared gate in the helper (Behavior 5). The issue's repro succeeds under `--no_report_file`.
5. **Shared helper (maintainer decision).** One seam both auditable against the writers and unit-testable:
   ```rust
   /// Clumping-report paths a clump-only run plans — one per report-keyed input,
   /// nothing when --no_report_file (matching every writer's gate).
   fn clump_report_candidates(
       no_report_file: bool,
       report_inputs: &[impl AsRef<Path>],
       output_dir: Option<&Path>,
   ) -> Vec<PathBuf>
   ```
   Call sites: paired FASTQ namer closure (`&[r1, r2]`), SE FASTQ (`&cli.input` — note this arm currently builds `planned` as an iterator chain, `main.rs:553-557`; extend after collecting), SE BAM (`&cli.input`), Shape A (`std::slice::from_ref(&chunk[0])` per chunk), Shape B (`&cli.input`). The "which inputs" parameter is where a `chunk[1]` mistake could hide — that is exactly what the Shape A test discriminates (Validation 7).
6. **No naming changes.** Output paths, report paths, and report contents are byte-identical on every accepted run.
7. **Edge cases.** Odd input counts and N=1 `--paired` are rejected earlier (existing checks, unchanged). `--basename`: reports still key on input filenames (A1), so the repro shape rejects identically — and with `--basename` the primaries *cannot* collide (`foo_clumped_1/2`), which makes the `--basename` variant the sharpest isolation of the report path (Validation 8). Duplicate-mate reuse is *not* rejected by validate (see Context) — Behavior 2c depends on that and tests it.

## Signature

`run_specialty_paired`'s namer widens from `(PathBuf, PathBuf)` to `Vec<PathBuf>` — "every path this pair writes". The tuple encoded the bug: it cannot express a third planned output (B§5.1). Parameter renamed `output_names` → `pair_outputs` (not r1's `planned_outputs`, one character from the local `planned` — A§1.8). Callers: clock/implicon wrap their two names in `vec![…]`; the clump arm appends `clump_report_candidates(…)`. Rejected alternatives unchanged from r1 (optional secondary closure; duplicate pre-flight in the arm), plus one added for the record (review request): **discriminating the report name instead of rejecting** — rejected because report filenames are a published contract (`*_clumping_report.txt` deliberately dodges MultiQC's `*_trimming_report.*` glob, `io.rs:463-468`) and #388 chose rejection for the identical class; two answers to one defect class is worse than either.

## Implementation outline

1. `main.rs::run_specialty_paired` — widen `NameFn` to `-> Vec<PathBuf>`; `planned.extend(pair_outputs(&chunk[0], &chunk[1]))`; update the doc comment (mode list + "every path the pair writes" contract).
2. `main.rs` — add `clump_report_candidates` (Behavior 5) near `planned_secondary_outputs` (`main.rs:89`), whose doc-comment style it follows.
3. `--clock`/`--implicon` call sites — wrap existing names in `vec![…]`.
4. `--clump_only --paired` FASTQ arm — namer returns primaries + `clump_report_candidates(cli.no_report_file, &[r1, r2], output_dir)`; hint `None` → `Some(PAIRED_REPORT_HINT)`; update the line-530 comment (it justified `None` against `CWD_OUTPUT_HINT`).
5. Task 2 — the four sibling arms extend their planned lists via the helper (inputs per Behavior 5). Shape B framed as defensive symmetry in comment and CHANGELOG.
6. Hint constant — reword `PAIRED_REPORT_HINT` per Behavior 3; rewrite its doc comment (first constant only), naming the three-plus-one user sites.
7. **io.rs property test (both reviews' top item)** — extend `distinct_primary_outputs_imply_distinct_secondary_outputs` (`io.rs:1229`) and/or `primary_output_key_is_coarser_than_secondary_keys` (`io.rs:1293`) so the clump namers (`clumped_output_name`, `clumped_bam_output_name`, `clumped_paired_bam_output_name` vs `clumping_report_name`) are covered, comparing via **`collision_key`** (the existing `PathBuf` inequality misses the fold dimension — A§4.2). Add `_clumped_N`/#391 alongside the "`_val_N` inverts this, #388" note in the doc comment at `io.rs:1287-1291`. This machine-checks A3, which licenses the uBAM report-vs-report non-fix AND proves Task 2 rejection-neutral (its entire behavioural delta is report-vs-input — B§1.4).
8. Integration tests in `tests/integration_output_collision.rs` (unique `tempdir` tag per test — tags share `process::id()` and `remove_dir_all` first, so a shared tag lets one test wipe another mid-run and make test 2 pass vacuously — B-V5b):
   - `clump_paired_rejects_shared_report_name` — the repro; `assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG)` **plus** separate asserts: stderr contains `reads.fq_clumping_report.txt` (pins the collision to the report) and contains `--no_report_file` (the durable hint fragment only — B-A2).
   - `clump_paired_accepts_shared_report_name_with_no_report_file` — same fixtures + flag; success; both primaries content-verified; zero report files. The matched negative control.
   - `clump_paired_accepts_distinct_filenames` — upgraded to `assert_dir_holds_only(&out, &[the four exact filenames])` (B-V4): pins "candidate list == written set", the invariant whose absence is this whole family's root cause.
   - `clump_paired_rejects_cross_pair_report_collision` — pairs (`a/x.fq`,`a/y.fq`)+(`b/y.fq`,`b/x.fq`), `-o out`; assert rejection ONLY — candidate order makes the flagged filename pair 2's report, so no filename pin (A§4.4).
   - **New** `clump_paired_rejects_shared_mate_report_without_output_dir` — pairs (`p/a.fq`,`d/x.fq`)+(`q/b.fq`,`d/x.fq`), no `-o` (Behavior 2c). Case-free, filesystem-independent; also pins the per-mate report-directory rule.
   - **New** `--basename` variant of test 1 — `--basename foo` forces distinct primaries by construction; reports still collide; rejected.
   - `clump_se_rejects_report_that_aliases_an_input` — **respec'd** (r1's shape couldn't run: `assert_rejected_cleanly` demands an empty dir and the alias must sit in the report's directory): no `-o`, one dir holding `s.fq` + `s.fq_clumping_report.txt` (both valid FASTQ — content-based detection accepts a `.txt` name, and validity keeps the rejection attributable to the pre-flight, not a parse error); assert `ALIAS_MSG` + the alias filename (the input branch names one path only), `count_reads_from`-verified intact inputs, `assert_dir_holds_only(&dir, &[both inputs])`. Model: `:674`. Plus a `--no_report_file` acceptance sibling.
   - **New** `clump_paired_bam_rejects_report_that_aliases_an_input` (Shape A — A-Critical-1): inputs `a/x.fq a/y.fq a/x.fq_clumping_report.txt a/z.fq` (all valid FASTQ), `--clump_only --paired --output-format ubam`, no `-o`; pair 1's report == input 3 → rejected, `ALIAS_MSG`, all four inputs intact. **Discriminates `chunk[0]` vs `chunk[1]`**: keyed wrong, the planned report would be `a/y.fq_clumping_report.txt`, nothing collides, and the run proceeds.
9. **Expected-fail controls on every rejection test** (A§4.4): run each against the unpatched build first; each must fail there (A verified all exit 0 pre-patch; the SE alias case *destroys an input* pre-patch). Record the pre-patch results in the PR body.
10. `CHANGELOG.md` — Unreleased → `#### Bug fixes` (the heading's actual name): reports join the clump pre-flights (#391); credit Task 1's free report-vs-input closure; Shape B named as defensive symmetry.
11. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`. (r1's Validation-7 grep is pre-answered by both reviewers: nothing pins `GENERIC_ADVICE`; the only hint pins are `CWD_OUTPUT_HINT`'s at `:625`/`:661`, untouched. The reworded `PAIRED_REPORT_HINT` has no test pins today.)

## Efficiency

Candidate list grows from 2P to ≤4P (paired FASTQ) and by ≤1 per input/pair elsewhere; the pre-flight stays two hash passes, O(paths). One `Vec` allocation per pair replaces a tuple, on a path about to read gigabytes. No new I/O; the pre-flight still precedes every reader open.

## Integration

- Writers untouched; report names/contents unchanged; byte-identity on accepted runs is the safety claim, which is why the "one shared derivation for namer and writer" refactor (B§5.2) is deliberately NOT taken — the exact-set assertion (step 8, test 3) buys most of that protection for one line.
- `run_specialty_paired` is `main.rs`-internal; no test references its signature.
- Perl-parity CI matrix does not cover `--clump_only`. Existing tests verified unaffected by both reviewers (incl. `pe_byte_identity`, `multi_pair_pe_bam_produces_one_output_per_pair`, `pe_bam_collision_preflight_case_folded` — the last still trips on a *primary*, which precedes its report in candidate order).

## Assumptions

- **A1:** `clumping_report_name` ignores `--basename` — verified (`io.rs:469-482`).
- **A2:** clock/implicon plan nothing beyond two primaries — verified (zero report/fastqc writes in `specialty.rs`).
- **A3:** uBAM/SE report-vs-report collisions imply primary collisions — re-derived by both reviewers; premises now explicit (fold-equivariant stripping; same-input same-dir keying; global `gzip` from `cli.input[0]`, `main.rs:351`; multi-input `--basename` rejected at `cli.rs:624/651`). Machine-checked by step 7 from this change on.
- **A4 (sharpened):** FastQC side-outputs on the clump paths derive from the discriminated primaries (`clump_only.rs:551-554`) and so cannot collide where primaries don't — benign here, not merely "out of scope".
- **A5 (corrected):** the pre-flight's output-vs-output error names both paths; the output-vs-**input** branch names only the input (`io.rs:109-115`). For report-vs-report between two mates the two displayed paths are the *same string* — pre-existing, unfixable without input provenance in the pre-flight; follow-up material (B-A3).
- **A6 (new):** content-based format detection accepts FASTQ content under any filename (first-byte `@` probe) — the premise of both alias tests.
- **A7 (new):** Shape A's report keys on `inputs[0]` (`clump_only.rs:1083`), never `inputs[1]` — pinned by the Shape A test.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Issue repro rejected, attributed to the report | test 1 (DUP_MSG + report filename + hint fragment) | exit ≠ 0, out dir empty |
| 2 | Gate parity | test 2 (same fixtures, `--no_report_file`) | exit 0, primaries content-verified, zero reports |
| 3 | Candidate list == written set | test 3 exact-set assertion | exit 0, exactly the four expected files |
| 4 | Cross-pair surface (`-o`) | test 4 | exit ≠ 0 (no filename pin) |
| 5 | Case-free no-`-o` class (Behavior 2c) | shared-mate test | exit ≠ 0 on ext4 AND APFS |
| 6 | Report-vs-input, SE | respec'd alias test + acceptance sibling | exit ≠ 0, `ALIAS_MSG`, inputs intact / sibling exit 0 |
| 7 | `chunk[0]` keying, Shape A | Shape A alias test | exit ≠ 0, `ALIAS_MSG`, four inputs intact |
| 8 | Report path isolated from stems | `--basename` variant of test 1 | exit ≠ 0 with primaries provably non-colliding |
| 9 | A3 held by CI | io.rs property test over clump namers via `collision_key` | passes; fails if any clump namer's report key ever becomes finer than its primary key |
| 10 | Fix actually fires | expected-fail control on EVERY rejection test vs the unpatched build | each fails pre-patch; results recorded in PR body |
| 11 | No collateral | full `cargo test` + clippy + fmt | green |

## Questions or ambiguities

- **[Resolved r2 — Felix]** Hint: fix the shared constant. Helper vs inline: shared helper. Test placement: `integration_output_collision.rs` (closed by both reviewers in r1's favour).
- **[Open]** None blocking. Copy of the reworded hint is reviewable prose; its durable fragment (`--no_report_file`) is the only pinned part.

## Follow-ups (out of scope, to file separately)

- Report-vs-report collisions display the same path twice and cannot name which *input* to rename — input provenance in `preflight_output_collisions` (B-A3).
- Report/primary directory asymmetry on the paired arms (reports follow each mate's parent; primaries follow R1's) — shared with #388's site (B§1.3).

## Implementation notes (2026-08-08, branch `fix/391-clump-report-preflight` @ `2704b52`)

Implemented as planned: `pair_outputs: … -> Vec<PathBuf>` widening, `clump_report_candidates(no_report_file, &[impl AsRef<Path>], output_dir)` helper feeding all five arms, `PAIRED_REPORT_HINT` reworded at the constant with its doc comment rewritten (first comment only), clump-paired hint `None` → `Some(PAIRED_REPORT_HINT)`, io.rs property tests extended over the clump namers via `collision_key` (+ case-variant input, + `_clumped_N`/#391 doc note), nine integration tests, CHANGELOG under `#### Bug fixes`. Full suite 557 green; fmt + clippy `-D warnings` clean.

**Expected-fail control (validation 10):** all six rejection tests run against the unpatched tree first — all six exited 0 there (the SE and Shape A alias cases destroying an input in the process); the three acceptance tests passed unpatched. Post-fix 9/9 green.

**Deviations (documented, none behavioural):**
- `assert_dir_holds_only`'s failure message generalized ("directory must hold exactly the expected files") since it now also serves acceptance-side exact-set assertions (B-V4).
- The `distinct_primary_outputs_imply_distinct_secondary_outputs` restructure loops over three primary namers; the demux-stem sub-check stays scoped to the SE trim namer (its semantics are SE-trim-specific) and is now guarded by `collision_key` inequality like the rest.

**Iteration log:**
- #1: initial implementation (`2704b52`); expected-fail control 6/6 pre-patch, 9/9 post.
- #2 (`d99f23e`): review batch per dual code review + coverage audit — coverage's one gap (alias-filename assertion) closed; A-M1/M2/M3, B-R1/R3/R4/R5/R6/R7, A-L1/L3/L4 applied; suite 559 green. Review items deliberately NOT taken: A-L5/B-I9 (helper stays in main.rs beside its neighbours), A-L6/B-duplication (trim-path unification deferred to when those blocks next change), A-L8 (CHANGELOG bullet order — Felix's call), B-I7 (helper doc history sentence kept for local consistency), B-I10 (half-line on why CWD_OUTPUT_HINT is wrong here — covered by the constants' own docs).

- **Logic:** every reviewer-cited anchor re-verified against the tree before adoption (io.rs test names/lines, integration helpers/constants, `cli.rs:588/602/624/651` guards). The helper's one hazard (wrong `report_inputs` at a call site) is pinned by Validation 7, which was designed to discriminate exactly that mistake.
- **Corrections owned:** r1's four prose errors are corrected in place and called out in the revision history rather than silently rewritten.
- **Edge cases:** Behavior 2c added the shape r1's `-o`-only analysis missed; `--basename` and duplicate-mate interactions now stated with their true mechanisms.
- **Remaining risks:** the reworded hint is shared prose across three pre-existing sites — its accuracy there was reviewed (the old clause was equally false at #388's site), but implementation must re-read those two sites' surrounding comments when swapping the text.
