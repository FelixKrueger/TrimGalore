# Code Review A — #400 `--clump_only` help text

**Reviewer:** A (independent; Reviewer B works in parallel with no shared state)
**Target commit:** `0676fdf` — `docs(cli): make --clump_only's help text format-aware (#400)`
**Base:** `dev` @ `ac08a97`
**Files in diff:** `src/cli.rs` (variant doc + four-paragraph `clump_only` rustdoc), `CHANGELOG.md` (one entry)
**Method:** every claim in the new copy re-derived from the code (`main.rs` five `--clump_only` dispatch arms, `cli.rs::Cli::validate` §3.4a + clump block, `io.rs` namers, `clump_only.rs` writers, `bam.rs` reader/writer), cross-read against `docs/src/content/docs/modes/clump-only.md`, then **empirically exercised** with the built binary (10 invocations, listed in §Evidence).

---

## Verdict

**APPROVE WITH CHANGES.** The rewrite is a real improvement and lands all five falsified claims the plan set out to fix; I verified each one independently against the code and, for the filename claims, against actual output files. Nothing in the diff can affect runtime, and `fmt` / `clippy -D warnings` / `test` are green.

Two accuracy problems survive into the shipped copy, both of the same class #400 exists to correct — an unqualified claim that is only true on one output format:

1. **A1 (High)** — paragraph 1's byte-identity claim now covers uBAM output, where the BAM writer **silently uppercases lowercase bases** and **coerces IUPAC degenerate bases to N**. Demonstrated on a 3-record fixture: FASTQ→FASTQ preserves both, FASTQ→uBAM does not. The old copy escaped this only by (falsely) claiming FASTQ-only scope; the rewrite extends the claim into territory the code does not honour.
2. **A2 (Medium)** — paragraph 1 read alone is the `-h` surface, and alone it implies uBAM-in → FASTQ-out is legal. That is the illegal 2×2 the plan explicitly caught in r1's draft; the constraint that rules it out now sits in paragraph 4, which `-h` does not render.

Plus one flag mis-scoped (**A3**, `--preserve-tags` needs uBAM *input*, not just uBAM output — hard error otherwise, verified), one plan step under-delivered (**A4**, multi-pair), and small editorial items. The CHANGELOG entry is accurate, correctly sectioned, and in house register.

---

## Area 1 — Accuracy, claim by claim

### Verified TRUE (no action)

| Claim (new copy) | Verified against | Result |
|---|---|---|
| "reorder FASTQ or unaligned BAM (uBAM) records" | `format::open_sync_reader` at `clump_only.rs:764`; five dispatch arms `main.rs:529-748` | OK |
| "compression-friendly grouping" (was "gzip-friendly compression") | BGZF on the BAM arms, gzip on FASTQ | OK, better than old |
| SE FASTQ `*_clumped.fq(.gz)` | `io.rs:356-379` | OK |
| PE FASTQ `*_clumped_{1,2}.fq(.gz)` | `io.rs:384-409` | OK |
| SE uBAM `*_clumped.bam` | `io.rs:417-431`; **ran it** → `BS-seq_10K_R1_clumped.bam` | OK |
| PE uBAM **ONE** interleaved `*_clumped.bam` per pair | `io.rs:442-458`; **ran it** → exactly one file; multi-pair N=4 → 2 files | OK (highest-value fix in the diff) |
| `*_clumping_report.txt` emitted, distinct from `*_trimming_report.txt` | `io.rs:463+`; **ran it** on every arm | OK |
| `--compression` **(FASTQ output only)** | never passed to the three BAM entry points (`main.rs:638-650`, `:695-707`, `:733-745`); recorded as level `0` at `clump_only.rs:780` | OK |
| `--dont_gzip` **(FASTQ output only)** + rejected with `--output-format ubam` | `cli.rs:585-590`; **ran it** → exact message the copy describes | OK |
| `--memory`, `--cores`, `--paired`, `--fastqc`, `--basename` compose | all threaded to both format arms; **ran `--fastqc` on the BAM arm** → `_clumped_fastqc.{html,zip}` produced | OK |
| `--output-format ubam` moved out of the rejected list | no such bail in `cli.rs:751-892`; `cli.rs:875-878` states the opposite | OK (falsified claim #2 fixed) |
| Rejected trimming/filtering list (`-a`, `--length`, `--rrbs`, `--polyA`, `--polyG`, `--trim-n`, `--clip_*`, `--nextseq`, `--rename`, `--discard_untrimmed`, `--consider_already_trimmed`, specialty modes, `--passthrough`, `--retain_unpaired`) | one bail each, `cli.rs:769-888` | OK, all 14 present |
| `-q` / `--stringency` / `-e` silently ignored | `cli.rs:744-750` rationale; no bail | OK |
| `@PG` appended on uBAM output → whole-file identity not preserved | `bam.rs:677-681` `add_program`; **ran `samtools view -H`** → `@PG ID:trim_galore VN:2.3.0 CL:…` present | OK (falsified claim #4 fixed) |
| "uBAM input requires `--output-format ubam` (the FASTQ output path would drop aux tags)" | SE `clump_only.rs:265-272`, PE `format.rs:165-172`; **ran both** → message matches the copy near-verbatim | OK (falsified claim #5 fixed) |
| Paired "a single interleaved uBAM" with `--output-format ubam` | **ran Shape B** on the interleaved BAM produced above → 10 000 pairs, one output | OK |
| Variant doc's seven rejections | §3.4a is the **only** `OutputFormat::UBam` block in `cli.rs` (`cli.rs:580` is the single `matches!` site) and holds exactly seven bails: `--dont_gzip` `:585`, `--clumpify` `:591`, `--passthrough` `:596`, `--clock` `:599`, `--implicon` `:605`, `--demux` `:610`, `--retain_unpaired` `:618` | OK — **complete and correct**, matches the code 1:1 |

Validation row 1 also passes: `grep -rn "natural follow-up\|FASTQ in / FASTQ out" src/ tests/` → zero hits.

### A1 (High) — "byte-identically … sequence" is false on the uBAM-output path

New paragraph 1:

> Every input record appears in the output byte-identically — name, sequence, quality, plus any aux tags named in `--preserve-tags` on uBAM; only file-level order changes.

`bam.rs::validate_and_normalize_seq_for_write` (`:800-826`) runs on the write side of **every** uBAM output arm:

- `b.to_ascii_uppercase()` — lowercase input bases are **silently** uppercased, no warning;
- `R Y M K S W B D H V` → `N`, with a once-per-run warning;
- `=` and any non-IUPAC byte → hard error.

The read side (`bam.rs:927-946`) applies the same IUPAC→N coercion, so **uBAM → uBAM is affected too** — which is the docs page's own headline archival use case (10X `CB`/`UB` BAMs; the warning text itself names "PacBio HiFi, ONT Dorado, and 10x cellranger uBAMs" as places IUPAC codes appear).

Proved, not inferred. A 3-record fixture with one lowercase read and one `RYKM` read:

```
FASTQ -> FASTQ (control):  acgtACGT...  preserved    ACGTRYKM...  preserved
FASTQ -> uBAM:             ACGTACGT...  UPPERCASED   ACGTNNNN...  COERCED
```

The old copy was not exposed to this, because it declared "v1 is FASTQ in / FASTQ out only" — and on the FASTQ arm the record is a byte copy, so the claim held exactly. By making the block format-aware **without** scoping the fidelity sentence, the rewrite converted a true statement into one that is false on the new path. That is precisely the defect class the issue was filed against, and it fails the plan's own pass criterion ("every claim is either true for both output formats or explicitly scoped to one").

`docs/…/modes/clump-only.md` §Record fidelity is silent on this too, so the specification the plan leaned on does not cover it — a genuine gap, not a transcription slip.

**Recommendation.** Not a trivial fix: it needs a maintainer call on wording and probably a companion docs edit. Minimal in-help form, appended to the fidelity sentence:

> … only file-level order changes. On uBAM output the sequence is uppercased and IUPAC degenerate bases are coerced to `N` (BAM's 4-bit alphabet); FASTQ output is a byte copy.

Two follow-ups worth filing regardless of the help wording:
- the docs page's §Record fidelity list should carry the same caveat;
- **the lowercase coercion is silent** while IUPAC warns. For a mode whose selling point is archival losslessness, silent case loss on the BAM path deserves at least the same one-time warning.

### A2 (Medium) — paragraph 1 alone implies the illegal uBAM-in → FASTQ-out combination

Confirmed from the built binary that `-h` renders **only** paragraph 1, as the plan predicted. Read on its own, paragraph 1 says the mode reorders "FASTQ or unaligned BAM (uBAM) records" and then gives the FASTQ filenames as the unmarked default with the BAM ones behind `--output-format ubam`. The natural inference is the full 2×2: uBAM in, FASTQ out, giving `*_clumped.fq.gz`.

That combination is rejected on both arms. The plan's Context flagged exactly this inference as the reason r1's draft was unacceptable ("implied 2×2 and would have invited exactly the refused invocation") — the shipped copy fixes it in paragraph 4, the surface `-h` users never see. So the defect was not removed, it was relocated to the short-help surface.

Mitigation in the field is decent: the rejection message is precise and self-remediating. But the point of separating `-h` from `--help` is that many users never run the long form.

**Recommendation (trivial fix).** A few words in paragraph 1, e.g. in the filename sentence:

> … or, under `--output-format ubam` (**required for uBAM input**), `*_clumped.bam` (SE) or ONE interleaved `*_clumped.bam` per pair (PE).

### A3 (Medium) — `--preserve-tags` is scoped to uBAM *output*, but it needs uBAM *input*

New paragraph 2 lists `` `--preserve-tags` (uBAM only) `` immediately after two entries scoped `(FASTQ output only)`. In that parallel construction "(uBAM only)" reads as "uBAM output only" — and that reading is wrong:

```
$ trim_galore --clump_only --output-format ubam --preserve-tags CB sample.fq.gz
Error: --preserve-tags has no effect with all-FASTQ inputs; either remove the flag
       or convert at least one input to uBAM via 'samtools import'
```

`main.rs:286-298` — with `--output-format ubam` set, `--preserve-tags` on all-FASTQ inputs is a **hard error**, and that guard sits at `main.rs:286`, ahead of the `--clump_only` dispatch at `:522`. So `--preserve-tags` requires at least one uBAM **input**; uBAM output alone is not enough.

Paragraph 1's phrasing ("aux tags named in `--preserve-tags` on uBAM") is fine — it is the composes-list entry that mis-scopes.

The docs page carries the same imprecision ("`--preserve-tags TAG1,TAG2,…` — uBAM only"), so this was inherited rather than invented — but this diff is the first time the flag appears in `--help` at all, so it ships the imprecision to a new audience.

Compounding it, `--preserve-tags`' own clap doc (unchanged, rendered in `-h`) says "**Ignored for FASTQ input**" — flatly contradicting the hard error above. Out of this diff's scope, but the same defect class as #400 and worth a follow-up issue.

**Recommendation (trivial fix).** `` `--preserve-tags` (uBAM input only) ``.

### A4 (Low–Medium) — plan step 4 asked for multi-pair; the copy dropped it

Paragraph 4 ships as "Paired mode takes two files, or — with `--output-format ubam` — a single interleaved uBAM." Plan step 4 specified "two files (Shape A, **multi-pair supported**) or one interleaved uBAM (Shape B)".

Multi-pair is real and I ran it: `--clump_only --paired --output-format ubam R1 R2 R1' R2'` prints `=== Pair 1 of 2 ===` / `=== Pair 2 of 2 ===` and writes two `_clumped.bam` files (`main.rs:683-717`); the FASTQ arm gets the same shape from `run_specialty_paired`. "Takes two files" is true but reads as an upper bound, and the docs page states the multi-pair support explicitly.

Second, smaller gap in the same sentence: it never says two separate uBAM files are refused. A user with `R1.bam R2.bam` gets a good error, but the help does not pre-empt it.

**Recommendation (trivial fix).** "Paired mode takes two FASTQ files (or several pairs — an even count), or — with `--output-format ubam` — a single interleaved uBAM; two separate uBAM files are not accepted."

### A5 (Low) — paragraph 3 still says "header", paragraph 1 now says "name"

Paragraph 1 changed `(header, sequence, quality)` → `name, sequence, quality`; paragraph 3 kept "byte-identity applies to **header** + sequence + quality bytes".

Both words matter, in opposite directions, and `bam.rs::parse_name_and_data` (`:707-722`) is why: on a FASTQ header `@NAME DESC`, the space-separated description is **discarded** on BAM output (`:717-718`, "Space: rest is descriptive FASTQ annotation — discard"). Confirmed in the dump — `@read1 some description here` came out as QNAME `read1`, description gone.

So:
- paragraph 1's "name" is the **more accurate** word for uBAM output, but under-promises on FASTQ output, where the entire header line including the description *is* preserved byte-identically (verified in the control run);
- paragraph 3's retained "header" **over-promises** on uBAM output.

The docs page already has the resolution, per format: "Read ID (header line for FASTQ; name portion for BAM)". As shipped, the two paragraphs disagree with each other and no reader can resolve that from the help alone.

**Recommendation (trivial fix).** Make paragraph 3 match the docs formulation — "byte-identity applies to the read ID (the full header line on FASTQ output, the name portion on uBAM output), sequence, and quality bytes" — or at minimum use one word in both paragraphs.

### A6 (Low) — paragraph 3's plus-line / CRLF sentences are unscoped FASTQ claims

"The plus-line (line 3 of each record) is normalized to bare `+` on output; CRLF line endings are normalized to LF. Both normalizations are codebase-wide behaviours, inherited from the FASTQ reader/writer."

Neither concept exists on the BAM path — BAM has no plus-line and no line endings. The statements are inapplicable rather than false, so the harm is nil, but under the plan's own pass criterion they are the last unscoped FASTQ-only claims in the block, one sentence away from the `@PG` sentence that *is* scoped. A five-word lead-in ("On FASTQ output, the plus-line…") closes it. Optional.

### Adjacent bug found while verifying (out of scope — recommend a follow-up issue)

The IUPAC warning fires with its **direction reversed** on the write side. Running FASTQ input → uBAM output:

```
WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V);
         coerced to N for FASTQ output.
```

The input was FASTQ and the output was BAM. `bam.rs` shares one `emit_iupac_warning_once()` between the read path (`:944`, where the wording is right) and the write path (`:823`, where it is exactly backwards). Not caused by this diff; same user-facing-accuracy family as #389/#400.

---

## Area 2 — Old vs new, side by side

Nothing true was deleted. Paragraph by paragraph:

| Old | New | Assessment |
|---|---|---|
| P1 "reorder FASTQ records … gzip-friendly compression" | "FASTQ or unaligned BAM (uBAM) … compression-friendly grouping" | improvement |
| P1 "byte-identically (header, sequence, quality)" | "byte-identically — name, sequence, quality, plus any aux tags named in `--preserve-tags` on uBAM" | the aux-tag half is the fix for falsified claim #3 and is right; "name" trades one imprecision for another (**A5**); the sentence's scope is now wrong on sequence (**A1**) |
| P1 "Output files use the … suffix" (FASTQ only) | per-format, incl. the ONE-interleaved-BAM PE shape | the diff's most valuable change |
| P2 `--compression`, `--dont_gzip` unqualified | both "(FASTQ output only)" | correct |
| P2 `--output-format ubam` in the **rejected** list | in the composes list; `--dont_gzip`+ubam kept visible as an explicit rejection | correct, and keeping the exclusion visible rather than silently dropping it is the right call |
| P2 rejected trimming list | unchanged apart from removing `--output-format ubam` | all 14 entries still accurate |
| P3 | + one `@PG` sentence | correct |
| P4 "v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up." | two sentences on input formats and paired shapes | false statement replaced by true ones; **A4** notes what the plan asked for and the copy omitted |
| Variant doc: 5 rejections | 7 | complete, exactly matches §3.4a |

**One newly introduced inaccuracy: A1.** **One newly introduced internal inconsistency: A5.** No other regressions.

Variant-doc completeness footnote (Low, optional): two other uBAM-output rejections live in `main.rs`, not `validate()` — `--phred64` + BAM input (`:264`) and `--preserve-tags` + all-FASTQ input (`:293`). The doc's "rejected at validation" arguably scopes them out, and the plan deliberately limited the edit to §3.4a, so I would not block on it; the `--preserve-tags` one is the trap most worth a mention if the line is ever revisited.

---

## Area 3 — CHANGELOG

**Correct as shipped.** No changes required.

- **Section:** `#### Changes` under `### Unreleased` (`CHANGELOG.md:143`), which is right — no behaviour change, and it avoids the duplicated `#### Bug fixes` / `#### Fixes` hazard the plan flagged (tracked in #401). Placed at the head of the section, matching the newest-first ordering of its neighbours (#400 → #389 → #383).
- **Register:** bold lead clause, parenthesised issue link, what-was-wrong-then-what-changed, closing "No behaviour change." Reads like its neighbours; no AI-slop phrasing; the `issues/` link form matches every sibling entry.
- **Accuracy:** all eight factual assertions in the bullet check out against the diff and the code — the old text's rejected-list placement, the "v1 is FASTQ in / FASTQ out only" close, the one-interleaved-BAM PE shape, the two FASTQ-output-only qualifiers, `--preserve-tags` as the aux-tag opt-in, the `@PG` reason, uBAM-input-requires-the-flag, and the two added variant-doc rejections.
- Leading with the filename correction as "most consequential" is the right emphasis and matches what I independently judged the highest-value fix.
- Only observation: the entry does not mention the P1 fidelity-enumeration change (`header` → `name` + aux tags). If **A1** is addressed that sentence needs a clause anyway, so this resolves itself.

---

## Area 4 — Structure and render

Rendered from the built binary at 80 and 100 columns.

- **`-h` renders paragraph 1 only** — confirmed, exactly as the plan assumed. Paragraph 1 is format-aware, so the fix is visible on the surface r1 would have missed. Right call.
- **`--help` renders all four** paragraphs, blank-line separated, correctly indented. Em-dashes render as literal em-dashes in all three places, matching house convention.
- **Length.** Paragraph 1 grew from ~540 to 833 characters (≈10 wrapped lines at 80 cols). That is a lot for a "short" help entry — but this binary's `-h` already prints every flag's full first paragraph (`--compression` and `--memory` are comparable), so it is consistent with its surroundings rather than an outlier. Not a blocker.
- **Does paragraph 1 carry what a short-help user needs?** Mostly yes — what the mode does, both input formats, the fidelity contract, per-format output filenames, the report name. It does **not** carry the uBAM-input constraint, which is **A2** and the one thing I would add. Everything else missing from `-h` (the `--dont_gzip`+ubam rejection, the paired shapes) produces a precise error at runtime, so the omission is tolerable.
- **Ordering within paragraph 1.** The filename list — the highest-value item for anyone scripting — now sits at roughly character 470 of 833, behind the fidelity clause. Editorial only; worth a thought if the paragraph is reopened for **A1**/**A2**.
- **Grammar nit (trivial fix).** "… or under `--output-format ubam` a single `*_clumped.bam` (SE) **and** ONE interleaved `*_clumped.bam` per pair (PE)" — the SE/PE alternation wants "or", not "and"; as written it can read as though both are produced. Also "a single" adds nothing for SE and slightly conflicts with SE multi-input runs, which write one BAM per input.
- **Variant doc line.** `- ubam:  Unaligned BAM output. Always single-threaded; …` renders as one unwrapped ~180-character line when stdout is not a TTY (clap only wraps against a real terminal width); in a terminal it wraps normally. Adding two flag names made a long line longer, not a new problem. Fine.

---

## Area 5 — Efficiency

Nil, as expected: the diff is doc comments and a changelog entry. No call sites, no runtime paths, no allocation or I/O effect. Nothing to review.

---

## Evidence

Build and checks at `0676fdf`: `cargo build --release` clean; `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean; `cargo test` exit 0.

Invocations run (outputs to a scratch directory via `-o`):

1. `--clump_only --output-format ubam` (SE FASTQ) → `BS-seq_10K_R1_clumped.bam`
2. `--clump_only --paired --output-format ubam R1 R2` → **exactly one** `BS-seq_10K_R1_clumped.bam` + one report
3. `--clump_only --paired --output-format ubam R1 R2 R1' R2'` → `Pair 1 of 2` / `Pair 2 of 2`, two BAMs (**A4**)
4. `--clump_only --paired --output-format ubam interleaved.bam` (Shape B) → 10 000 pairs, one output
5. `--clump_only --dont_gzip --output-format ubam` → rejected, message as documented
6. `--clump_only --output-format ubam --preserve-tags CB` + FASTQ input → **hard error** (**A3**)
7. `--clump_only sample.bam` (no `--output-format ubam`) → rejected, message as documented
8. `--clump_only --paired --output-format ubam a.bam b.bam` → two-BAM rejection
9. `--clump_only --output-format ubam --fastqc` → `_clumped_fastqc.{html,zip}`
10. lowercase + IUPAC fixture, FASTQ→FASTQ vs FASTQ→uBAM, then `samtools view -h` → **A1** and **A5** proven

---

## Recommendations by priority

**High**
1. **A1** — scope the byte-identity sentence, or state the uBAM-output sequence normalisation (uppercasing + IUPAC→N). Needs a wording decision; also file the docs-page gap and consider warning on the silent lowercase coercion.

**Medium**
2. **A2** *(trivial fix)* — put the "uBAM input requires `--output-format ubam`" constraint into paragraph 1, so the `-h` surface stops implying the illegal 2×2.
3. **A3** *(trivial fix)* — `--preserve-tags` → "(uBAM **input** only)".

**Low**
4. **A4** *(trivial fix)* — say multi-pair is supported, and that two separate uBAM files are not, in paragraph 4.
5. **A5** *(trivial fix)* — reconcile paragraph 3's "header" with paragraph 1's "name"; the docs page's per-format formulation is the model.
6. **A6** *(trivial fix, optional)* — scope paragraph 3's plus-line / CRLF sentences to FASTQ output.
7. Grammar nit in paragraph 1's filename sentence: "and" → "or"; drop "a single" (Area 4).
8. Follow-up issues, all outside this diff: the reversed IUPAC warning direction (`bam.rs:823`); `--preserve-tags`' own "Ignored for FASTQ input" line, which contradicts the hard error; `--output-format`'s "see `--help`" self-reference inside `--help`.

**No changes required**
- CHANGELOG entry.
- Variant doc's seven rejections — complete and correct.
- All output-filename claims, confirmed against real output files.

No files in the worktree were edited, per instruction.
