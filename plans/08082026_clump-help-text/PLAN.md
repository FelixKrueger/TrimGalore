# Plan: `--clump_only` help text — make the whole block format-aware (#400)

**Issue:** [#400](https://github.com/FelixKrueger/TrimGalore/issues/400) — filed from #391's code review.
**Base:** `dev` @ `ce13761`.

## Revision history

- **r2 (2026-08-08):** incorporates dual plan review (PLAN_REVIEW_A.md "REVISE", PLAN_REVIEW_B.md "APPROVE WITH CHANGES") and two maintainer decisions: **no docs-site pointer — inline the shapes** and **fold in the `OutputFormat::UBam` variant doc**. Both reviewers rejected r1's two-sentence scope and both caught that r1's own proposed copy would have introduced **two new false statements**. Scope grows from 2 sentences to the whole four-paragraph block plus the variant doc (~12 lines of prose). Line references corrected (block is `cli.rs:180-206`; r1 cited 183-207, and 207 is the `#[clap]` attribute).
- r1 (2026-08-08): initial plan — two sentences. Superseded.

## Goal

The clap rustdoc on `clump_only` (`cli.rs:180-206`) renders as `--help` and predates the shipped uBAM arms. **Five** claims in it are falsified by those arms, not the two #400 named. Make the whole block format-aware so no claim is FASTQ-only-but-unqualified, and fix the same defect class in the `OutputFormat::UBam` variant doc. No runtime behaviour changes.

## Context

- **The block spans four paragraphs, and they are two different surfaces** (B §1.1, verified empirically from the built binary): `--help` renders all four; **`-h` renders only paragraph 1**, because clap's short help for an argument is the first doc paragraph. r1 touched paragraphs 2 and 4 only — so on the `-h` path the fix would have been invisible while paragraph 1 kept naming the wrong files.
- **The five falsified claims:**
  1. `cli.rs:184-188` — output filenames given unqualified as `*_clumped.fq(.gz)` (SE) / `*_clumped_{1,2}.fq(.gz)` (PE). Under `--output-format ubam`: `<stem>_clumped.bam` (SE, `io.rs:417-431`) and **ONE interleaved** `<stem>_clumped.bam` per pair (PE, `io.rs:442-458` — a different output *count*, not just extension). Both reviewers rate this the highest-impact item in the block; users script against filenames.
  2. `cli.rs:195` — `--output-format ubam` listed among *rejected* flags. `Cli::validate` has no such rejection; `cli.rs:866-869` carries the opposite comment ("v2: `--output-format ubam` is now accepted").
  3. `cli.rs:183-184` — byte-identity enumerated as "(header, sequence, quality)", which omits uBAM aux tags (they round-trip **only** when named in `--preserve-tags`, A/Z/i/f scalars only, `clump_only.rs:714-718`). A 10X user reordering a `CB`/`UB`-carrying BAM without `--preserve-tags` loses them silently under a heading that says lossless.
  4. `cli.rs:200-204` — the contract-scope paragraph closes on "**Both** normalizations", reading as exhaustive, but the uBAM path appends a per-run `@PG` record so whole-file identity does not hold (`clump_only.rs:720-723`).
  5. `cli.rs:206` — "v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up."
- **Two flags in the "Composes with" list are format-dependent**, which r1's step 1 would have contradicted by appending `--output-format ubam` to the same list: `--dont_gzip` + ubam is a hard error (`cli.rs:576-581`), and `--compression` is inert on the BAM path (never passed to the three BAM entry points; recorded as level `0` in `clump_only.rs:780`/`:1018`).
- **Only 3 of 4 format combinations are legal.** uBAM-in → FASTQ-out is rejected in both arms (`clump_only.rs:265-272` SE; `format.rs:165-172` PE via `PairedShape::ClumpOnlyFastqOut`) because FASTQ output would drop aux tags. r1's proposed "FASTQ and uBAM are both supported on input and output" implied 2×2 and would have invited exactly the refused invocation.
- **`OutputFormat::UBam`'s variant doc** (`cli.rs:17-18`) lists 5 of the 7 ubam rejections — missing `--retain_unpaired` (`cli.rs:609`) and `--dont_gzip` (`cli.rs:576`) — and **does** render in `--help` under "Possible values". In scope per maintainer decision.
- **`docs/src/content/docs/modes/clump-only.md` is the specification.** Both reviewers independently read it end to end and verified it accurate on every point at issue: per-format filenames, the one-interleaved-BAM PE shape, `--compression` ignored under ubam, `--dont_gzip` rejected with ubam, aux tags via `--preserve-tags` with the A/Z/i/f constraint, the added `@PG`, the three PE input-shape rejections, and the FASTQ-out-with-uBAM-in rejection. No companion edit needed; use it as the source.
- **Zero pinning risk, resolved now rather than at implementation** (both reviewers): `tests/integration_no_args_help.rs:36` pins only `"Usage: trim_galore"` and `"--adapter"`; nothing else in `tests/` or `src/` asserts on help output; no snapshot tests; `release.yml:226`/`:356` run `--help` as an exit-status smoke test with no content grep; no `cargo doc` CI job and no `[lints]` in `Cargo.toml`, so rustdoc lints are not a gate.
- **Two hits that must NOT be edited:** `CHANGELOG.md:278-279` and its synced mirror `docs/src/content/docs/reference/changelog.md:28-29` carry the same stale-sounding wording as *historical release notes* (true at the release they describe; the synced page is off-limits by convention). Keep the verification grep scoped to `src/ tests/`. Also `main.rs:524`'s code comment ("v1: FASTQ in/out. v2: uBAM in/out via `--output-format ubam`") is accurate and not user-facing — leave it.

## Behavior

`--help` and `-h` output change; the `--output-format` "Possible values" text gains two flag names. No runtime behaviour changes, no flag semantics change.

## Implementation outline

Re-read the block from disk before editing (post-`cargo fmt` anchor drift is a recorded session trap), then rewrite all four paragraphs against the docs page. B's §6 offers copy verified line-by-line against the code; adopt it with the maintainer's no-pointer decision applied.

1. **Paragraph 1** (`cli.rs:180-188`) — the `-h` surface. Make the opener format-aware ("FASTQ or unaligned BAM (uBAM) records"), extend the byte-identity enumeration to name aux tags as a `--preserve-tags` opt-in, and give the output filenames per format including the PE-uBAM single-interleaved case. Keep the `*_clumping_report.txt` sentence (verified true on every arm, `clump_only.rs:836`/`:1082`).
2. **Paragraph 2** (`cli.rs:190-198`) — qualify `--compression` and `--dont_gzip` as FASTQ-output-only, move `--output-format ubam` out of the rejected list into the composes list, add `--preserve-tags` (uBAM only), and keep the `--dont_gzip` + ubam exclusion **visible** as an explicit rejection rather than silently dropping it. Leave the trimming/filtering rejection list and the `-q`/`--stringency`/`-e` sentence as-is (both verified accurate).
3. **Paragraph 3** (`cli.rs:200-204`) — add one sentence: on uBAM output a `@PG` record is appended, so whole-file identity is not preserved; record bodies are.
4. **Paragraph 4** (`cli.rs:206`) — replace the stale sentence. Must say uBAM input **requires** `--output-format ubam` (not that it is one of two options), and **state the two paired shapes inline** — two files (Shape A, multi-pair supported) or one interleaved uBAM (Shape B) — with **no docs-site pointer** (maintainer decision; would be the first such reference in `cli.rs`, and `CHANGELOG.md:291-293` records deliberately removing unresolvable pointers from flag descriptions).
5. **`OutputFormat::UBam` variant doc** (`cli.rs:17-18`) — add the two missing rejections (`--retain_unpaired`, `--dont_gzip`) so the "Possible values" text matches the shipped rejection set.
6. **CHANGELOG** — one bullet under `### Unreleased` → `#### Changes` (`CHANGELOG.md:143`). Precedent: `CHANGELOG.md:291-293` ("**Tidied the `--help` text** … No behaviour change") — same section, same kind of change, a stronger match than r1's `--output-dir` citation (which is a runtime stderr message). A separate bullet is justified: this corrects a **false capability claim**, and notably that tidy pass itself missed line 206. Name the filename correction explicitly — a user who scripted against `_clumped_{1,2}.fq.gz` for uBAM output was misled. **Hazard:** `### Unreleased` contains both `#### Bug fixes` and `#### Fixes` (duplication tracked in #401); land in `#### Changes`.
7. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`.

## Efficiency / Integration

Nil / doc-string only. No call sites, no runtime, no output-size effect.

## Assumptions

- **A1 (verified, no longer deferred):** nothing pins the help text — see Context.
- **A2 (verified by both reviewers):** the docs page is accurate post-#396 and serves as the specification; no companion edit.
- **A3:** the terse-comment convention (one line default, two max) governs **code comments**, not user-facing `--help` prose, where completeness is the governing virtue — the project demonstrates this with 138 lines of docs prose for this one flag. This plan deliberately grows the block.
- **A4 (resolved contradiction):** the reviewers disagreed on whether `--cores` in the composes list is accurate. Checked directly: `cli.cores` **is** passed to all three BAM entry points (`main.rs:642`, `:699`, `:737`) and drives FastQC threads, while the reorder itself is single-threaded per `cli.rs:749-752`. Both were partly right; either way `--cores` is not falsified by uBAM and is **out of scope** — leave it unqualified.

## Validation

Rows 1–3 from r1 verify mechanics only; both reviewers flagged that the table had no row asserting the new copy is **true** — exactly where the risk lay. Rows 4–6 close that.

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Stale phrases gone from live text | `grep -rn "natural follow-up\|FASTQ in / FASTQ out" src/ tests/` (scope deliberate — the two `CHANGELOG`/synced-docs hits are historical and off-limits) | zero hits |
| 2 | **Both** help surfaces correct | rebuild, then render `-h` **and** `--help` | the uBAM facts a user needs appear on the surface that user hits; paragraph 1 no longer FASTQ-only |
| 3 | Build integrity (**not** a validation of the copy — no test asserts this text) | fmt + clippy `-D warnings` + `cargo test` | green |
| 4 | Whole-block accuracy vs the specification | read all four paragraphs + the variant doc side by side with `clump-only.md` §Output filenames (39-45), §Record fidelity (47-63), §Compatibility (69-105). **Pass criterion: every claim is either true for both output formats or explicitly scoped to one** | no claim contradicts the docs page |
| 5 | No newly-claimed flag is actually rejected | for each flag named as composing, confirm no matching `anyhow::bail!` in the `if self.clump_only` block (`cli.rs:742-883`) or the shared §3.4a `UBam` block (`cli.rs:571-616`) | none found (this is the mechanical check that would have caught r1's `--dont_gzip` defect at plan time) |
| 6 | Filename claims confirmed against reality, not just against docs | `--clump_only --output-format ubam test_files/BS-seq_10K_R1.fastq.gz` → expect `BS-seq_10K_R1_clumped.bam`; `--clump_only --paired --output-format ubam test_files/BS-seq_10K_R{1,2}.fastq.gz` → expect **exactly one** `_clumped.bam`; `--clump_only --dont_gzip --output-format ubam …` → expect the rejection the new copy describes | as stated (run outside the repo root — specialty modes write to the CWD) |

**Deliberately not added** (both reviewers, independently): a help-text regression test. Prose assertions rot faster than the prose they guard, the repo has zero help-content assertions by design, and the root cause is doc-drift-after-feature, which no unit test catches.

## Questions or ambiguities

- **[Resolved r2 — Felix]** No docs-site pointer; inline the shapes. `OutputFormat::UBam` variant doc folded in.
- **[Open, non-critical]** Exact copy is reviewable prose; the fixed requirements are the six accuracy points in Context and the "true for both formats or scoped to one" criterion.

## Follow-ups (out of scope, noted)

- `--phred64` + FASTQ-in + uBAM-out threads `cli.phred_offset()` into the BAM writers while the FASTQ→FASTQ arm ignores phred entirely (byte copy). Quality *values* are preserved; the ASCII string necessarily is not, since BAM stores raw Phred. Semantically correct, asymmetric with the FASTQ arm, immaterial to help accuracy (A §7).

## Implementation notes (2026-08-08, branch `fix/400-clump-help-text`)

Base correction: the plan header says `dev @ ce13761`, but the branch forked from `ac08a97` (dev moved when the #385 re-audit artifacts landed). No finding is affected — the audited hunks are the `--clump_only` doc-comment, the `OutputFormat::UBam` variant doc, and one CHANGELOG bullet.

- `0676fdf` — the four-paragraph rewrite + variant doc + CHANGELOG, as planned. Validations 1/2/3/5/6 all run, incl. the three smoke invocations (SE → `R1_clumped.bam`; PE → exactly one `_clumped.bam`; `--dont_gzip` + ubam rejected).
- `e8e79ca` — review batch. **Both reviewers APPROVED WITH CHANGES, and between them found two distinct ways the rewrite's own byte-identity claim was false** — a defect this PR introduced by making a FASTQ-only claim format-spanning: A found that uBAM output uppercases bases and coerces IUPAC to `N` (`bam.rs:800-820`); B found that it discards the space-separated header description (`bam.rs:717`), so a standard Illumina `1:N:0:INDEX` field does not survive. Confirmed on a fixture rather than by reading — `@READ_001 1:N:0:INDEX` / `acgtRYKM` emerges as QNAME `READ_001` / SEQ `ACGTNNNN`. Copy now states FASTQ→FASTQ as byte-identical and uBAM as lossless-per-record-body-but-not-byte-identical, with paragraph 3 scoped to FASTQ output.
  Also in this commit: the uBAM-input-requires-ubam-output rule moved into paragraph 1 (both reviewers — paragraph 1 is the only one `-h` renders, so short-help readers were told BAM input works and then refused; verified present in rendered `-h`); `--preserve-tags` reworded from "(uBAM only)" to "requires at least one uBAM input; an error otherwise" (it is a hard error, not inert — `main.rs:286`); and the paired clause gained both multi-pair support and the two-separate-uBAMs rejection, **which also closes the coverage audit's two PARTIAL items** (they traced to that same clause).
- `5407a37` + `d19aa75` — the docs page carried the identical false claim in both its Record-fidelity bullets and its intro paragraph. The plan designated that page as the specification (A2), so leaving it would have had the docs site contradicting the corrected help text. A2 is therefore **falsified in part** and recorded as such.
- **Not taken → filed:** #406 (the IUPAC warning states the wrong direction on FASTQ→uBAM runs, and neither the uppercase nor the description-drop normalization is warned about at all) and #407 (`--preserve-tags`' own help says "Ignored for FASTQ input" while the code hard-errors under ubam output).
- **Known, not addressed:** B measured that clap is built without `wrap_help` (`Cargo.toml:29`), so nothing wraps — paragraphs 1 and 2 are now the two longest lines in the whole `--help`. Accuracy was preferred over brevity here; if that becomes a readability complaint the fix is enabling `wrap_help`, not shortening the facts.

## Self-Review (r2)

- **What r1 got wrong, owned:** it inherited the issue body's "two-line fix" as a *scope conclusion* rather than a size estimate, and its Context recorded the `--dont_gzip` + ubam exclusion without any implementation step consuming it — so its own proposed copy would have asserted two mutually exclusive flags compose, and implied a format combination the code refuses. Both reviewers caught both independently.
- **Load-bearing facts re-verified against the tree before adoption:** the `-h`/`--help` paragraph split, `io.rs:417-458`'s single-interleaved-BAM naming, `cli.rs:576-581`, the two uBAM-in→FASTQ-out rejection sites, the variant doc's missing two, and the `--cores` disagreement (A4).
- **Traps checked:** grep scope keeps the two historical changelog records off-limits; validation 2 rebuilds before rendering (stale-binary trap); validation 3 is explicitly labelled a build check, not evidence about the copy (the project's "prove the check can fail" rule); line references corrected so an implementer editing "207" doesn't hit the `#[clap]` attribute.
- **Remaining risks:** none structural — the change cannot affect runtime. The residual is editorial: prose accuracy, which validations 4–6 target directly.
