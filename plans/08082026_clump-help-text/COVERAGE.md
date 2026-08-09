# Plan Coverage Report

**Mode:** B (code vs plan). No separate `IMPL.md` — ledger built from the plan's **Implementation outline** (7 steps), **Behavior** (3 claims), **Assumptions** (A1–A4) and **Validation** table (6 rows), plus the explicitly "deliberately not added" help-text regression test, recorded as a decision rather than a gap.
**Plan:** `/Users/fkrueger/Github/TrimGalore/plans/08082026_clump-help-text/PLAN.md` (r2)
**Audited commit:** `0676fdf` "docs(cli): make --clump_only's help text format-aware (#400)", diffed against `ac08a97`
**Date:** 2026-08-08
**Verdict:** **INCOMPLETE — 2 items unresolved** (both PARTIAL, both trace to one clause in the final sentence of paragraph 4; no functional impact, prose completeness only)

## Summary

- Total audited items: 20
- DONE: 18
- PARTIAL: 2
- MISSING: 0
- DEVIATED: 0
- Recorded decisions (not counted): 1

Diff footprint is `CHANGELOG.md` (+15) and `src/cli.rs` (+28/−19) only — doc-comments and changelog prose, no executable lines. `git diff --stat ac08a97..0676fdf` confirms no other file was touched, so A2's "no companion docs edit" and the off-limits changelog/synced-docs records hold by construction.

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Paragraph 1 (`cli.rs:181-191`) — `-h` surface made format-aware: opener names uBAM, byte-identity enumeration names aux tags as a `--preserve-tags` opt-in, filenames given per format incl. the PE-uBAM single-interleaved case, `*_clumping_report.txt` sentence retained | Outline 1 | DONE | All four sub-requirements present. Verified rendered on the `-h` path (line 62 of `-h`, single paragraph as predicted). `_clumping_report.txt` sentence kept; emission confirmed empirically. |
| 2 | Paragraph 2 (`cli.rs:193-203`) — `--compression`/`--dont_gzip` qualified FASTQ-output-only; `--output-format ubam` moved from rejected to composes; `--preserve-tags` (uBAM only) added; `--dont_gzip` + ubam kept visible as an explicit rejection; trim/filter rejection list and the `-q`/`--stringency`/`-e` sentence unchanged | Outline 2 | DONE | All six sub-requirements present verbatim. `--output-format ubam` no longer in the rejected list; the `--dont_gzip` + ubam exclusion is stated, not dropped. |
| 3 | Paragraph 3 (`cli.rs:205-210`) — one added sentence: `@PG` appended on uBAM output, whole-file identity not preserved, record bodies are | Outline 3 | DONE | Exactly one sentence added, matching `clump_only.rs:720-723` and docs §Record fidelity (61-63). |
| 4 | Paragraph 4 (`cli.rs:212-215`) — stale sentence replaced; must say uBAM input **requires** `--output-format ubam`; **state the two paired shapes inline** — two files (Shape A, **multi-pair supported**) or one interleaved uBAM (Shape B); no docs-site pointer | Outline 4 | **PARTIAL** | Stale sentence gone ✓; "requires" wording ✓; both shapes stated inline ✓; no pointer ✓ (grep for `http`/`/modes/`/`docs.` inside the rendered block returns 0). **Missing:** the "multi-pair supported" detail the step names. Multi-pair *is* supported (`main.rs:546` routes `--clump_only --paired` through `run_specialty_paired`, which iterates `cli.input.chunks(2)`) and docs:73 states it, so this is an omitted true detail, not a correction. See Gap 1. |
| 5 | `OutputFormat::UBam` variant doc (`cli.rs:17-19`) — add the two missing rejections (`--retain_unpaired`, `--dont_gzip`) | Outline 5 | DONE | Variant doc now lists all 7 rejections; matches the shipped §3.4a set exactly (`cli.rs:585` dont_gzip, `:591` clumpify, `:596` passthrough, `:599` clock, `:605` implicon, `:610` demux, `:618` retain_unpaired). Renders under "Possible values" in `--help` (line 128). |
| 6 | CHANGELOG — one bullet under `### Unreleased` → `#### Changes`, naming the filename correction explicitly; avoid the `#### Bug fixes` / `#### Fixes` mis-landing hazard | Outline 6 | DONE | Bullet at `CHANGELOG.md:145`, inside `#### Changes` (`:143`) inside `### Unreleased` (`:4`) — the hazard sections are `:6` and `:310`, both avoided. Names the filename correction and the `_clumped_{1,2}.fq.gz` scripting risk explicitly; also records the variant-doc addition; closes "No behaviour change." |
| 7 | `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test` | Outline 7 | DONE | All three green — see Test verification. |
| 8 | `--help` and `-h` output change | Behavior | DONE | Both rebuilt-and-rendered. `-h` shows the rewritten paragraph 1; `--help` shows all four paragraphs. |
| 9 | `--output-format` "Possible values" text gains two flag names | Behavior | DONE | `--retain_unpaired` and `--dont_gzip` appear in the rendered `ubam` value description. |
| 10 | No runtime behaviour change, no flag semantics change | Behavior | DONE | Diff touches only `///` doc-comments and `CHANGELOG.md`. No executable statement, no `#[clap]` attribute, no validation logic changed. |
| 11 | A1 — nothing pins the help text | Assumption | DONE (holds) | `tests/integration_no_args_help.rs:36` still pins only `"Usage: trim_galore"` and `"--adapter"`; neither string is in the edited block. No new help-content assertion anywhere in `tests/` or `src/`. The full suite passes unmodified. |
| 12 | A2 — the docs page is the specification; no companion edit | Assumption | DONE (holds) | `docs/src/content/docs/modes/clump-only.md` is untouched by the diff. Read end to end for row 4; it remains accurate on every point at issue. |
| 13 | A3 — help prose is exempt from the terse-comment convention; the block deliberately grows | Assumption | DONE (holds) | Block grew from 26 to 38 doc-comment lines; growth is intentional per the plan. |
| 14 | A4 — leave `--cores` unqualified, out of scope | Assumption | DONE (holds) | `--cores` appears in the composes list with no format or threading qualifier, as directed. |
| 15 | V1 — stale phrases gone from live text (`grep -rn "natural follow-up\|FASTQ in / FASTQ out" src/ tests/`) | Validation 1 | DONE | **Zero hits** (exit 1). The two off-limits historical records survive untouched: `CHANGELOG.md:294` and `docs/src/content/docs/reference/changelog.md:29`. `main.rs:524`'s accurate code comment left in place as directed. |
| 16 | V2 — both help surfaces correct after a rebuild | Validation 2 | DONE | Rebuilt (`cargo build --release`, exit 0), then rendered both. Paragraph 1 is no longer FASTQ-only: on the `-h` path a user now sees the uBAM opener, the aux-tag scoping, and the per-format filenames incl. `ONE interleaved *_clumped.bam per pair`. `--help` additionally carries paragraphs 2-4 and the variant doc. |
| 17 | V3 — build integrity (explicitly *not* evidence about the copy) | Validation 3 | DONE | fmt clean, clippy `-D warnings` clean, 561 tests pass across 13 binaries. |
| 18 | V4 — whole-block accuracy vs the specification; pass criterion: **every claim is either true for both output formats or explicitly scoped to one** | Validation 4 | **PARTIAL** | 22 of 23 claims pass. One clause fails the criterion: "Paired mode takes two files" carries no format scoping, and under `--output-format ubam` two BAM files are a hard error. See Gap 2 (and the three observations below it, which pass). |
| 19 | V5 — no newly-claimed composing flag is actually rejected | Validation 5 | DONE | Checked all 9 composing flags against the `if self.clump_only` block (`cli.rs:751-892`) and the shared §3.4a `UBam` block (`cli.rs:580-625`). No matching `anyhow::bail!` for `--compression`, `--memory`, `--paired`, `--fastqc`, `--basename`, `--output-format ubam`, `--preserve-tags`. `--dont_gzip`'s only bail is the ubam one the copy explicitly documents. `--cores`' only bail is the global `cores == 0` floor (`:706`), not a mode rejection. `--memory` is parsed for format only (`:890`). |
| 20 | V6 — filename claims confirmed against reality | Validation 6 | DONE | All three invocations run from a temp dir with copied fixtures. See Test verification. |
| — | Deliberately not added: a help-text regression test | Plan decision | N/A (decision) | Recorded as a documented decision, not a gap. Confirmed no such test was added: `tests/` gained no file and no help-content assertion. |

## Gaps (detail)

### Gap 1 — Outline step 4: the "multi-pair supported" detail is absent

**Expected:** step 4 specifies the two paired shapes be stated inline as "two files (Shape A, **multi-pair supported**) or one interleaved uBAM (Shape B)". The docs page, which the plan designates as the specification, states it at line 73: "`--paired` — two-file paired input (Shape A, multi-pair `N=4, 6, …` supported) OR single interleaved uBAM (Shape B, N=1 + `--output-format ubam`)".

**Found:** `cli.rs:214-215` — "Paired mode takes two files, or — with `--output-format ubam` — a single interleaved uBAM." Both shapes are named; the multi-pair capability is not.

**Gap:** multi-pair input under `--clump_only --paired` is genuinely supported (`main.rs:546` → `run_specialty_paired`, which iterates `cli.input.chunks(2)`), so this is an omitted true capability rather than a deliberate correction of the plan. The plan's Revision history frames the four-paragraph rewrite as the deliverable and names this parenthetical explicitly, so its absence is a shortfall against the spec. Nothing in the shipped copy contradicts multi-pair support.

### Gap 2 — Validation row 4: "Paired mode takes two files" is not format-scoped

**Expected:** row 4's pass criterion — "every claim is either true for both output formats or explicitly scoped to one".

**Found:** `cli.rs:214-215` — "Paired mode takes two files, or — with `--output-format ubam` — a single interleaved uBAM." The "takes two files" clause is unqualified. Under `--output-format ubam` it holds only when both files are FASTQ:

- Two BAM files: hard error, `format.rs:182-191` — "`--clump_only --paired` with two BAM files is not supported. uBAM paired mode expects a single interleaved file" (docs:99).
- N=1 FASTQ under `--paired --output-format ubam`: rejected (docs:100).
- Mixed FASTQ/BAM in one pair: rejected, `format.rs:196-203` (docs:101).

**Gap:** the sentence would let a user with two uBAM files read "uBAM input requires `--output-format ubam`" plus "Paired mode takes two files" and attempt `--clump_only --paired --output-format ubam R1.bam R2.bam`, which bails. This is the same failure shape the plan's own Self-Review flags against r1 ("implied a format combination the code refuses"), one step removed: the copy does not assert the refused shape works, but it does not scope the clause that suggests it. The runtime error is clear and offers a `samtools merge -n` remediation, so the consequence is a wasted invocation, not wrong output. Same source sentence as Gap 1.

### Observations from row 4 that pass the criterion (recorded, not gaps)

- **Paragraph 3's contract-scope enumeration** still reads "byte-identity applies to header + sequence + quality bytes" — FASTQ vocabulary, and it omits the aux tags paragraph 1 now names. Not falsified on either format (those three fields *are* byte-identical on both paths), and the plan scoped step 3 to adding one sentence, so this is per-spec.
- **The A/Z/i/f aux-tag constraint** (docs:54) is not in the help. No help claim contradicts it: B (array) and H (hex) tags are rejected at BAM-read time, so "appears byte-identically" cannot be silently breached — the run fails instead.
- **`--cores`'s "single-threaded internally"** (docs:76) is not in the help. Explicitly out of scope per A4; `--cores` is accepted on both format paths, so this is a depth qualification, not a format-scoping defect.
- **`--basename`** is listed as composing without noting it bails on multiple inputs/pairs (`cli.rs:633`, `:660`). Pre-existing and format-independent; the docs page (79) is equally unqualified. Outside this plan's scope.

## Test verification

| Check | Command / file | Status |
|---|---|---|
| V1 — stale phrases in live text | `grep -rn "natural follow-up\|FASTQ in / FASTQ out" src/ tests/` | PASS (zero hits) |
| V1 — off-limits records preserved | `CHANGELOG.md:294`, `docs/.../reference/changelog.md:29` | PASS (both present, untouched) |
| V2 — `-h` renders paragraph 1 only, format-aware | `./target/release/trim_galore -h` | PASS |
| V2 — `--help` renders all four paragraphs + variant doc | `./target/release/trim_galore --help` | PASS |
| V3 — formatting | `cargo fmt --all -- --check` | PASS (exit 0) |
| V3 — lints | `cargo clippy --all-targets --release -- -D warnings` | PASS (exit 0) |
| V3 — test suite | `cargo test` | PASS — 561 tests, 0 failed (lib 410; `integration_clump_only` 12; `integration_clump_only_ubam` 15; `integration_ubam_out` 23; `integration_output_collision` 45; `integration_adapter2` 15; `integration_gzip_non_gz_extension` 11; `integration_paired_format_guard` 11; `integration_non_restartable_input` 8; `integration_ubam` 8; `integration_passthrough` 2; `integration_no_args_help` 1; doc-tests 0) |
| V5 — no composing flag is rejected | read-through of `cli.rs:751-892` + `cli.rs:580-625` | PASS (9/9 flags clear) |
| V6a — SE uBAM filename | `--clump_only --output-format ubam BS-seq_10K_R1.fastq.gz` | PASS — wrote `BS-seq_10K_R1_clumped.bam` (10 000 records, 16 bins, BGZF) exactly as the copy claims |
| V6b — PE uBAM writes exactly one BAM | `--clump_only --paired --output-format ubam BS-seq_10K_R{1,2}.fastq.gz` | PASS — **one** `BS-seq_10K_R1_clumped.bam` (10 000 pairs); `ls \| grep -c "_clumped.bam$"` = 1, confirming the "ONE interleaved per pair" claim |
| V6c — `--dont_gzip` + ubam rejected as described | `--clump_only --dont_gzip --output-format ubam …` | PASS — exit 1, "`--dont_gzip` is not compatible with `--output-format ubam` (BAM is always BGZF-compressed)", matching the copy's "(BAM is always BGZF)" |
| Filename claims vs source | `io.rs:356` / `:384` / `:417` / `:442` | PASS — `_clumped.fq(.gz)` SE, `_clumped_{1,2}` PE, `<stem>_clumped.bam` SE-BAM, single `<stem>_clumped.bam` PE-BAM |
| Report emission on the BAM arm | observed `*_clumping_report.txt` in V6a/V6b; `clump_only.rs:837`, `:1084` | PASS |
| A1 — help text unpinned | `tests/integration_no_args_help.rs:36` | PASS (pins only `Usage: trim_galore` and `--adapter`; test green) |

## Verdict

**INCOMPLETE — 2 items unresolved.** Both are PARTIAL, both editorial, and both originate in the same clause of the final sentence at `cli.rs:214-215`. Every mechanical and behavioural requirement of the plan is satisfied: the five falsified claims the plan enumerates are all corrected, the variant doc matches the shipped 7-rejection set, the changelog bullet landed in the right section, and rows 1, 2, 3, 5 and 6 pass — including the empirical filename check that confirms a paired uBAM run writes exactly one interleaved BAM.

To close:

1. **Outline step 4 / Gap 1** — add the multi-pair detail the step names, e.g. that the two-file shape accepts multiple pairs (`N = 4, 6, …`), matching docs:73.
2. **Validation row 4 / Gap 2** — scope the "takes two files" clause by output format, so a user with two uBAM files is not led toward `--paired --output-format ubam R1.bam R2.bam` (rejected at `format.rs:182-191`). One qualifier — that the two-file shape is FASTQ input — satisfies row 4's criterion and closes both gaps in the same edit.

Neither item can affect runtime; the whole diff is doc-comments plus changelog prose.

## Notes

- The plan's header records **Base: `dev` @ `ce13761`**, while this audit was run against base `ac08a97` (per the audit assignment). The diff reviewed is `ac08a97..0676fdf`. The discrepancy does not affect any finding — the audited hunks are confined to the `--clump_only` doc-comment, the `OutputFormat::UBam` variant doc, and one `CHANGELOG.md` bullet.
- Row 6 was run from a temp directory with copied fixtures, as the plan requires (specialty modes write beside the input).
