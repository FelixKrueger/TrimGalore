# Code Review B — #400: make `--clump_only`'s help text format-aware

**Reviewer:** B (independent, fresh context)
**Target:** commit `0676fdf` (`docs(cli): make --clump_only's help text format-aware (#400)`), base `ac08a97`
**Files in diff:** `src/cli.rs` (clap rustdoc on `clump_only`, 4 paragraphs; `OutputFormat::UBam` variant doc), `CHANGELOG.md` (one bullet)
**Plan:** `plans/08082026_clump-help-text/PLAN.md` (r2)
**Governing virtue for this review:** completeness + accuracy of user-facing `--help` prose (per plan A3, the terse-code-comment rule does not apply here).

## Verdict

**APPROVE WITH CHANGES.** The rewrite fixes all five claims the plan set out to fix, and every new claim I checked against the code holds — the filename shapes, the `--compression`/`--dont_gzip` FASTQ-only scoping, the `@PG` sentence, and the "uBAM input requires `--output-format ubam`" rule are all verified true, several empirically. The variant doc's seven rejections are exactly the seven the `§3.4a` block enforces. Build gates are green and the plan's own validation rows 1 and 6 pass.

Two things keep this from a clean APPROVE, and both are residual instances of the exact defect class #400 was filed to eliminate:

1. **Paragraph 3 still carries an unqualified FASTQ-only claim** — "byte-identity applies to header + sequence + quality bytes". Under `--output-format ubam` the space-separated FASTQ description is **discarded**; I demonstrated it. This contradicts paragraph 1 of the same block (which correctly says "name") and contradicts the docs page that the plan designated as the specification.
2. **The one cross-format legality rule is invisible on `-h`.** Paragraph 1 — the only paragraph `-h` renders — now advertises uBAM input with no hint that it requires `--output-format ubam`. That constraint lives in paragraph 4. The plan's own Context named this `-h`/`--help` asymmetry as r1's fatal flaw; the shipped copy reproduces it one level down.

Both are one-clause edits. Nothing here is a correctness or runtime risk.

---

## Verification performed

| Check | Method | Result |
|---|---|---|
| Build | `cargo build --release` | ok |
| Format | `cargo fmt --all -- --check` | clean |
| Lint | `cargo clippy --all-targets --release -- -D warnings` | clean |
| Tests | `cargo test` | **561 passed, 0 failed** (410+15+12+15+11+1+8+45+11+2+8+23) |
| Plan validation row 1 | `grep -rn "natural follow-up\|FASTQ in / FASTQ out" src/ tests/` | zero hits |
| Plan validation row 2 | rebuilt, rendered `-h` and `--help` | `-h` = paragraph 1 only (confirmed); `--help` = all four |
| Plan validation row 6a | `--clump_only --output-format ubam BS-seq_10K_R1.fastq.gz` | `BS-seq_10K_R1_clumped.bam` ✓ |
| Plan validation row 6b | `--clump_only --paired --output-format ubam R1 R2` | **exactly one** `BS-seq_10K_R1_clumped.bam` ✓ |
| Plan validation row 6c | `--clump_only --dont_gzip --output-format ubam …` | rejected with the documented message ✓ |
| Header-fidelity probe (mine, not in the plan) | FASTQ with `@READ_001 some:description here` → uBAM; BAM parsed directly | QNAME = `READ_001` — **description dropped** |

Code read for the accuracy audit: the five `--clump_only` dispatch arms (`src/main.rs:522-757`), the `§3.4a` `UBam` rejection block (`src/cli.rs:579-624`), the `if self.clump_only` block (`src/cli.rs:751-892`), `src/main.rs:286-295` (`--preserve-tags` §3.4b), the four clumped namers (`src/io.rs:356-464`), `src/clump_only.rs:250-272` and `:690-800`, `src/bam.rs:588-640` + `:707-723`, `src/format.rs:165-192`, and `docs/src/content/docs/modes/clump-only.md` end to end.

---

## Area 1 — Accuracy, claim by claim

### Claims that hold (verified against code, not the plan)

| New claim | Verified against | Verdict |
|---|---|---|
| "reorder FASTQ or unaligned BAM (uBAM) records" | 5 dispatch arms in `main.rs:522-757` | true |
| aux tags round-trip only when named in `--preserve-tags` | `clump_only.rs:715-717` rustdoc; `preserve_tags` threaded to the 3 BAM entry points only (`main.rs:644`, `:701`, `:739`) | true |
| SE FASTQ `*_clumped.fq(.gz)` | `io::clumped_output_name` (`io.rs:356-377`) | true |
| PE FASTQ `*_clumped_{1,2}.fq(.gz)` | `io::clumped_paired_output_names` (`io.rs:384-409`) | true |
| SE uBAM `*_clumped.bam` | `io::clumped_bam_output_name` (`io.rs:417-435`); empirically confirmed | true |
| **ONE interleaved** `*_clumped.bam` per pair (PE uBAM) | `io::clumped_paired_bam_output_name` (`io.rs:442-458`, `_input_r2` unused) + `flush_bin_paired_to_bam` writing `R1, Some(1)` then `R2, Some(2)` (`clump_only.rs:702-710`); empirically confirmed | true — the highest-value correction in the diff |
| `*_clumping_report.txt` emitted | `io::clumping_report_name`; both FASTQ and BAM arms | true (see NIT-9 for the two qualifiers it omits) |
| `--compression` **(FASTQ output only)** | never passed to any BAM entry point; recorded as level `0` in `clump_only.rs:780` | true |
| `--dont_gzip` **(FASTQ output only)** | `cli.rs:583-588` bail | true |
| `--memory`, `--cores`, `--paired`, `--fastqc`, `--basename` compose | all five threaded to both FASTQ and BAM arms | true |
| `--output-format ubam` composes (moved out of the rejected list) | no `output_format` rejection in the `clump_only` block; `cli.rs:872-875` comment confirms | true — this was the headline falsehood |
| rejected list (`-a`, `--length`, `--rrbs`, `--polyA`, `--polyG`, `--trim-n`, `--clip_*`, `--nextseq`, `--rename`, `--discard_untrimmed`, `--consider_already_trimmed`, other specialty modes, `--passthrough`, `--retain_unpaired`) | each has a matching `anyhow::bail!` in `cli.rs:768-885` | true |
| `--dont_gzip` + `--output-format ubam` rejected (BAM always BGZF) | `cli.rs:583-588`; empirically confirmed | true |
| `-q` / `--stringency` / `-e` silently ignored | `cli.rs:746-750` comment; no bail | true (unchanged) |
| `@PG` appended → whole-file identity not preserved, record bodies are | `clump_only.rs:720-723`; `bam::build_output_header`; empirically confirmed (`@HD VN:1.6` + one `@PG` with `CL:`) | true |
| uBAM input **requires** `--output-format ubam` | `clump_only.rs:265-272` (SE) and `format.rs:165-177` via `PairedShape::ClumpOnlyFastqOut` (PE) — the help's parenthetical "(the FASTQ output path would drop aux tags)" matches both bail messages in substance | true |
| Variant doc's seven rejections | `cli.rs:583-624`: `--dont_gzip`, `--clumpify`, `--passthrough`, `--clock`, `--implicon`, `--demux`, `--retain_unpaired` — exactly seven | complete and correct |

### MEDIUM-1 — Paragraph 3's byte-identity sentence is still FASTQ-only and unqualified (and is now false on the uBAM path)

Paragraph 3 was left almost intact:

> Contract-scope note: byte-identity applies to **header** + sequence + quality bytes.

On `--output-format ubam` this is false whenever a FASTQ header carries a space-separated description. `src/bam.rs:707-722`:

```rust
if sep == b'\t' {
    // Tab: rest is the aux tag tail (from BamReader::next_record).
    (name_str, Some(rest))
} else {
    // Space: rest is descriptive FASTQ annotation — discard.
    (name_str, None)
}
```

Demonstrated end to end. Input:

```
@READ_001 some:description here
```

Output BAM, header + records parsed directly from the BGZF stream:

```
@HD	VN:1.6
@PG	ID:trim_galore	PN:trim_galore	VN:2.3.0	CL:… --clump_only --output-format ubam desc.fastq
QNAME: 'READ_002'
QNAME: 'READ_001'
```

`some:description here` is gone. This is not an exotic shape — it is the standard Illumina header, `@<instrument>:…:<pos> 1:N:0:<INDEX>`, so the read-number / filter / index field is what gets dropped. A user who reads paragraph 3 as written and archives FASTQ into uBAM under a heading that says byte-identical loses it silently.

Three aggravating factors:

- **It contradicts paragraph 1 of the same block**, which the rewrite changed to "name, sequence, quality". Paragraph 1 is right for BAM; paragraph 3 is right for FASTQ; the block now says both without reconciling them.
- **It contradicts the designated specification.** `docs/src/content/docs/modes/clump-only.md:52` is precise: "Read ID (header line for FASTQ; name portion for BAM)". The plan's Context (and both plan reviewers) treated that page as accurate on every point at issue and as the source to write from. This sentence did not get written from it.
- **It fails the plan's own pass criterion** for validation row 4: "every claim is either true for both output formats or explicitly scoped to one".

**Trivial fix**, replacing the first sentence of paragraph 3:

> Contract-scope note: byte-identity covers the read ID, sequence and quality bytes — on FASTQ output the whole header line; on uBAM output the name only, since a space-separated FASTQ description has no BAM field to land in.

and, if paragraph 1 is to stay terse, "name" there can then read "read ID" and defer.

### MINOR-3 — "Paired mode takes two files" invites a refused invocation

Paragraph 4:

> FASTQ and uBAM input are both supported, but uBAM input requires `--output-format ubam` … Paired mode takes two files, or — with `--output-format ubam` — a single interleaved uBAM.

The preceding sentence has just told the reader uBAM input is supported; "takes two files" then reads format-agnostic. Two uBAM files under `--paired --output-format ubam` is rejected (`src/format.rs:181-192`, "`{mode} with two BAM files is not supported`"), and so is a mixed-format pair. The refusal is loud and its message is good, so the cost is a failed invocation rather than wrong output — but this is the same trap the plan flagged in r1's copy ("would have invited exactly the refused invocation"), just narrower.

**Trivial fix:** "Paired mode takes two FASTQ files, or — with `--output-format ubam` — either two FASTQ files or one interleaved uBAM; two separate uBAMs are refused as ambiguous."

### MINOR-4 — "`--preserve-tags` (uBAM only)" names the wrong side

The precondition is uBAM **input**, not uBAM output. With `--output-format ubam` and all-FASTQ inputs, `--preserve-tags` is a **hard error**, not inert (`src/main.rs:286-295`):

```rust
anyhow::bail!(
    "--preserve-tags has no effect with all-FASTQ inputs; either remove \
     the flag or convert at least one input to uBAM via 'samtools import'"
);
```

So `--clump_only --output-format ubam --preserve-tags CB sample.fq.gz` — a plausible reading of "composes with `--preserve-tags` (uBAM only)" — fails. **Trivial fix:** "(uBAM input only)".

*Adjacent, out of scope:* `--preserve-tags`'s own help (`cli.rs:246-248`) says "Ignored for FASTQ input", which is false under `--output-format ubam` for the same reason. Same defect class as #400, different flag — worth its own issue rather than scope creep here.

### MINOR-5 — Paragraph 1's "name" under-claims for FASTQ output

Old text said "(header, sequence, quality)"; new says "name, sequence, quality". On the FASTQ→FASTQ path the entire header line, description included, is byte-identical, so "header" was the accurate word there and "name" quietly gives up a true guarantee. "name" is the accurate word for BAM. Neither covers both. Folds into MEDIUM-1's fix — name both formats once and let the other paragraph defer.

### MINOR-6 — Paragraph 3's normalization sentences are unscoped FASTQ-writer behaviour

"The plus-line (line 3 of each record) is normalized to bare `+` on output; CRLF line endings are normalized to LF." There is no line 3 and no line ending on the BAM path. Not misleading (a BAM user finds it irrelevant, not wrong), but unscoped against the plan's stated criterion. Optional: prefix "On FASTQ output, ".

### Out-of-scope observation — `--cores`

The plan's A4 resolved this and put it out of scope; I agree with the outcome and record only that the help is now less precise than its own docs page. `--cores` is threaded to all three BAM entry points and to `clump::resolve_layout`, so it affects bin layout and FastQC threads, but the reorder itself is single-threaded on every arm. `clump-only.md:83` says so explicitly ("accepted for interface parity — v1 is single-threaded internally"); the help lists `--cores` unqualified among things it "composes with". Equally (im)precise for both formats, so #400's remit does not reach it. Leaving it is correct for this PR.

### Variant doc completeness

Complete and correct for flag-vs-flag rejections. One further uBAM-output-only hard error is not represented: `--preserve-tags` with all-FASTQ inputs (`main.rs:286`). Omitting it is defensible — it is conditional rather than flag-vs-flag, and it fires after format detection rather than "at validation" as the sentence says — so I would **not** ask for a change. Recorded for completeness only.

One cosmetic note: the seven are listed with `--dont_gzip` last, though it is the first to bail. Irrelevant to a reader.

---

## Area 2 — Old vs new, side by side

**New inaccuracies introduced:** none that I could find. Every claim in the new copy is true; the defects above are a *retained* false claim (MEDIUM-1), a scoping gap (MINOR-3, MINOR-4, MINOR-6) and an under-claim (MINOR-5).

**True statements dropped:**

- "header" → "name" (MINOR-5) — the one real loss.
- "gzip-friendly compression" → "compression-friendly grouping": a gain, not a loss; BGZF is not gzip in the user's mental model and the docs page uses the same neutral phrasing.
- "Output files use the … suffix" → "Output is …": neutral.
- `--output-format ubam` removed from the rejected list: correct, it is accepted.
- "v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up.": correctly deleted; grep confirms zero hits in `src/ tests/` and the two historical `CHANGELOG` / synced-docs records are untouched, as the plan required.

**Plan-coverage gap (NIT-10):** plan step 4 asked the paired-shapes sentence to state "two files (Shape A, **multi-pair supported**) or one interleaved uBAM (Shape B)". The shipped sentence drops the multi-pair clause. Multi-pair *is* supported on the uBAM PE path (`main.rs:658-756` iterates `chunks(2)` with per-pair banners) and the docs page says so (`clump-only.md:80`). Whether this was deliberate compression or an oversight is worth a one-line confirmation from the implementer.

---

## Area 3 — CHANGELOG

**Section: correct.** Lands under `### Unreleased` → `#### Changes` (line 143), not the `#### Bug fixes`/`#### Fixes` pair the plan flagged as a hazard (lines 6 and 310). Ordering within the section is newest-first (#400, #389, #383), consistent with its neighbours.

**Register: correct.** Bold lead clause, issue link in the first sentence, concrete before/after, `No behaviour change.` closer — matches the #389 and #383 bullets directly below it. No banned AI-slop phrasing.

**Accuracy: every claim verified true.** The `--output-format ubam`-in-the-rejected-list history, the "v1 is FASTQ in / FASTQ out only" quote, the one-interleaved-`*_clumped.bam` consequence, the `_clumped_{1,2}.fq.gz` scripting hazard, the FASTQ-output-only scoping of `--compression`/`--dont_gzip`, the `@PG` rationale, the requires-`--output-format ubam` rule, and the two added variant-doc rejections all check out against the diff and the code.

Two observations, neither blocking:

- **Length.** At 15 lines this is the longest bullet in `#### Changes` for the smallest change in it. The precedent the plan cited (`CHANGELOG.md:291-293`, "**Tidied the `--help` text** … No behaviour change") is three lines. The repo's recent entries do run long, so this is within house style — but the middle third ("The whole block is now format-aware: …" enumerating six items) duplicates what a reader gets from `--help` itself and could lose 4 lines without losing information.
- **Omission.** The entry does not mention the byte-identity-enumeration correction (header → name + `--preserve-tags` aux tags). That is the subtlest change in the diff, it is the one MEDIUM-1 shows is still incomplete, and it is the one a user archiving into uBAM most needs to know changed. If MEDIUM-1 is fixed, this sentence should be added at the same time.

---

## Area 4 — Structure and `--help` readability

**Paragraph 1 does *not* carry everything a short-help user needs — MEDIUM-2.**

Confirmed empirically that `-h` renders paragraph 1 only; the rendered text stops at "…scanners unconfused)" and the next line is `--compression`. Paragraph 1 now opens:

> Lossless reorder-only specialty mode: reorder **FASTQ or unaligned BAM (uBAM)** records …

and mentions `--output-format ubam` only in the output-filename clause. The rule that uBAM *input* **requires** `--output-format ubam` is in paragraph 4, which `-h` never shows. So the `-h` reader is told BAM input works, is shown that `--output-format ubam` changes output filenames, and is given nothing that connects the two — `trim_galore --clump_only sample.bam` looks sanctioned and is refused.

This is the same failure mode the plan diagnosed in r1 ("on the `-h` path the fix would have been invisible while paragraph 1 kept naming the wrong files"). The rewrite fixed the filenames on that surface but introduced a new `-h`-only gap in their place: the one legality rule that governs the input side.

**Trivial fix**, appending to paragraph 1's existing `--output-format ubam` clause: "… ONE interleaved `*_clumped.bam` per pair (PE); uBAM input requires `--output-format ubam`." Nine words, on the surface that needs them.

**Wrapping — MINOR-7, no change required.** `clap` is built without the `wrap_help` feature (`Cargo.toml:29`, `features = ["derive"]`), so clap performs **no** wrapping at all: each paragraph is emitted as one physical line and the terminal soft-wraps it mid-word. Measured on the built binary:

| Line | Chars |
|---|---|
| paragraph 1 (`--clump_only`) | **728** (was ~545) |
| paragraph 2 (`--clump_only`) | **709** (was ~620) |
| next longest line in all of `--help` | 513 |

So this flag now owns the two longest lines in the entire help output by ~40%, and `-h` shows a 728-character wall for a boolean flag. This is a pre-existing property of the CLI rather than a regression, and completeness was the agreed governing virtue (plan A3), so I am not asking for cuts. If the maintainer wants relief without losing content, moving the `*_clumping_report.txt` sentence from paragraph 1 to paragraph 2 removes ~180 characters from the `-h` surface — the report filename is not something a short-help reader needs before running the mode, whereas the input-legality rule from MEDIUM-2 is.

**NIT-8 — `--dont_gzip` appears twice in paragraph 2**, once as composing "(FASTQ output only)" and once in the rejection clause. Deliberate per plan step 2, and both statements are accurate, but on a single unwrapped 709-character line the two verdicts land close together and read as a contradiction on first scan. Optional reshape: drop the parenthetical from the composes list and let the rejection clause carry the whole story ("…`--dont_gzip` is rejected together with `--output-format ubam` (BAM is always BGZF), so it applies to FASTQ output only").

**NIT-9 — "A short `*_clumping_report.txt` is emitted"** is unconditional in prose but conditional in fact: `--no_report_file` suppresses it, PE FASTQ writes **two** (one per mate) while the uBAM shapes write **one per pair** (`clump-only.md:141`). Pre-existing wording, unchanged by this diff, and the glob form is not wrong. Mentioning `--no_report_file` in the composes list would close the loop cheaply.

**Paragraph split otherwise reads well.** Four paragraphs, one job each: what it does and what it produces / what composes and what does not / what byte-identity does not cover / input-format legality. That is the right decomposition, and paragraph 4's promotion from a one-line "future work" note to a real input-format paragraph is a clear improvement.

---

## Area 5 — Efficiency

Nil. Doc strings and a changelog bullet; no call sites, no runtime path, no output bytes. `cargo build --release` succeeds and the binary's behaviour is unchanged by construction.

---

## Recommendations by priority

### Should fix before merge

1. **MEDIUM-1** — scope paragraph 3's byte-identity sentence per format. It is currently false for `--output-format ubam` (FASTQ description dropped, empirically shown), it contradicts paragraph 1 of the same block, and it contradicts the docs page the plan designated as the specification. This is the last surviving instance of the defect class #400 was filed to remove. *Trivial fix.*
2. **MEDIUM-2** — add the "uBAM input requires `--output-format ubam`" rule to **paragraph 1**, the only paragraph `-h` renders. Nine words. *Trivial fix.*

### Should fix, low cost

3. **MINOR-4** — "`--preserve-tags` (uBAM only)" → "(uBAM input only)". With `--output-format ubam` and all-FASTQ inputs it is a hard error, not inert. *Trivial fix.*
4. **MINOR-3** — qualify "Paired mode takes two files" so two separate uBAMs are not implied to work. *Trivial fix.*
5. **CHANGELOG** — if MEDIUM-1 is fixed, add the byte-identity-enumeration correction to the bullet; it is currently the one substantive change the entry does not mention.

### Optional / maintainer's call

6. **MINOR-5 + MINOR-6** — say "read ID" rather than "name" in paragraph 1, and prefix paragraph 3's plus-line/CRLF sentences with "On FASTQ output". Both fall out of MEDIUM-1's rewrite for free.
7. **NIT-10** — confirm whether dropping plan step 4's "multi-pair supported" clause was deliberate.
8. **NIT-8** — reshape paragraph 2 so `--dont_gzip` states its verdict once.
9. **MINOR-7** — if `-h` bulk is a concern, move the `*_clumping_report.txt` sentence out of paragraph 1 (~180 chars off the short-help surface, nothing lost).
10. **Separate issue** — `--preserve-tags`'s own help says "Ignored for FASTQ input", false under `--output-format ubam`. Same defect class, different flag; out of scope here.

### Explicitly not recommended

- A help-text regression test. Both plan reviewers argued against it independently, the repo has zero help-content assertions by design, and MEDIUM-1 is proof that the failure mode is doc-drift-after-feature, which no prose assertion catches. Agreed.
