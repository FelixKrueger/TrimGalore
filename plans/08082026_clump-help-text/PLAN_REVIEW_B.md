# Plan Review B — `--clump_only` help text (#400)

**Plan:** `plans/08082026_clump-help-text/PLAN.md`
**Reviewer:** B (independent; Reviewer A ran separately — no shared state)
**Base:** `dev` @ `ce13761`
**Verdict:** APPROVE WITH CHANGES — both named defects are real, but the named set is
**incomplete** and the **proposed replacement sentence is itself inaccurate**.

---

## Summary

The plan correctly identifies two stale statements and correctly scopes the change as
doc-string-only with no runtime effect. Two problems keep it from being ready:

1. The defect set is incomplete in the way that matters most. The **first** paragraph of
   the same doc block promises `*_clumped.fq(.gz)` / `*_clumped_{1,2}.fq(.gz)` output — and
   that paragraph is the **only** one `-h` renders (verified empirically below). Under
   `--output-format ubam` the real names are `<stem>_clumped.bam` (SE) and **one
   interleaved** `<stem>_clumped.bam` per pair. Land the plan as written and `--help`
   advertises uBAM output in paragraph 2 while paragraph 1 still tells the user to look for
   two `_clumped_{1,2}.fq.gz` files.
2. The proposed replacement copy says "FASTQ and uBAM are both supported on input and
   output". uBAM **input** on the FASTQ **output** path is rejected — in both the SE and PE
   arms, by design, because it would drop aux tags. The plan would swap one falsehood for
   another in the one place users check.

Everything else the plan asserts checks out: the two line-level defects exist, no test or
CI job pins the phrases, the CHANGELOG heading exists, and the docs page needs no companion
edit.

---

## What I verified as correct in the plan

| Plan claim | Verdict | Evidence |
|---|---|---|
| `--output-format ubam` sits in the *rejected* list | **True** | `src/cli.rs:195` |
| A "v1 is FASTQ in / FASTQ out only" closing sentence exists | **True** | `src/cli.rs:206` (not 207) |
| uBAM in/out is shipped for SE, PE Shape A, PE Shape B | **True** | five arms at `src/main.rs:522-750` |
| `--preserve-tags` composes | **True** | threaded into all three BAM entry points |
| `--dont_gzip` is rejected *with* uBAM output | **True** | `src/cli.rs:576-581` (shared §3.4a block) |
| **A1** no test pins the stale phrases | **True — and knowable now** | see §3 |
| **A2** the docs page is already accurate | **True** | see §2 |
| `#### Changes` exists under `### Unreleased` | **True** | `CHANGELOG.md:4` / `:143` |
| Help text is not covered by the Perl-parity matrix | **True** | CI only runs `--help` for exit status |

The plan also does **not** conflate the project's terse-code-comment convention with
user-facing help prose — step 4 explicitly reasons "help text is user-facing". That is the
right call, and it matters, because my recommendation below *grows* the block. The one-line
comment rule must not be applied to `--help` copy.

---

## 1. Logic review

### 1.1 Critical — the fix is invisible on the `-h` path, and paragraph 1 stays false

Rendered from the built binary at `target/release/trim_galore` (source text matches
`src/cli.rs` verbatim, so it is representative):

- `--help` renders **all four** paragraphs of the `clump_only` rustdoc.
- `-h` renders **only the first paragraph** — clap's short help for an argument is the
  first doc paragraph; `long_help` is `--help`-only.

The plan touches only paragraph 2 (line 195) and paragraph 4 (line 206). So:

- A user who types `-h` learns nothing about uBAM support either before or after the fix,
  and is still told the output is `*_clumped.fq(.gz)` / `*_clumped_{1,2}.fq(.gz)`.
- A user who types `--help` gets a paragraph 2/4 that advertise uBAM output and a paragraph
  1 that contradicts them on filenames.

The naming claim is not a stylistic nit — it is the sentence that determines which files a
user goes looking for, and the PE case is the sharp one: **one** interleaved BAM per pair,
not two mate files.

Ground truth (`src/io.rs`):

- `clumped_bam_output_name` (`io.rs:417-431`) → `<stem>_clumped.bam`
- `clumped_paired_bam_output_name` (`io.rs:442-458`) → a **single** `<stem>_clumped.bam`,
  stem from R1, no `_1`/`_2` suffix ("Returns a SINGLE interleaved BAM path")

Paragraph 1's opening clause ("reorder **FASTQ** records") carries the same FASTQ-only
framing; the docs page title is already "Reorder **FASTQ or uBAM** records".

**Action:** paragraph 1 must be in scope. This is the sibling falsehood the fix would
otherwise leave behind.

### 1.2 Critical — the proposed replacement sentence introduces a new inaccuracy

Plan step 2 proposes: *"FASTQ and uBAM are both supported on input and output; select uBAM
output with `--output-format ubam` …"*

Read naturally, that offers four combinations. Only three exist. uBAM input on the FASTQ
output path is rejected in both arms:

- SE: `src/clump_only.rs:265-272` — `"uBAM input under --clump_only requires
  --output-format ubam (using the FASTQ output path with uBAM input would drop aux tags)"`
- PE: `src/format.rs:165-172`, via `PairedShape::ClumpOnlyFastqOut` — same message, fired
  for **any** BAM in the pair regardless of the other side

The docs page states this as its own bullet ("**FASTQ output path with uBAM input:** …
is rejected"). For a change whose entire purpose is accuracy, shipping this sentence would
be the worst outcome of the three (stale text, corrected text, or newly-wrong text).

**Action:** the sentence must say uBAM input *requires* `--output-format ubam`, not that it
is merely one of two options.

### 1.3 Important — adding `--output-format ubam` to "Composes with" contradicts two entries already in that list

Plan step 1 extends the "Composes with" sentence with `--output-format ubam`. The resulting
list would read `--compression, --memory, --cores, --paired, --fastqc, --dont_gzip,
--basename, --output-format ubam` — asserting that `--dont_gzip` and `--output-format ubam`
both compose, when they are **mutually exclusive**:

- `src/cli.rs:576-581` rejects `--dont_gzip` under `--output-format ubam` outright.

And `--compression` is **inert** on the uBAM path, not merely less useful:

- `main.rs:574` / `:604` pass `cli.compression` to the FASTQ entry points; the three BAM
  calls (`main.rs:638-650`, `695-707`, `733-745`) **never pass it**.
- `clump_only.rs:780` / `:1018` hardcode `compression_level: 0` with the comment
  "BGZF; label carries the story".

The plan's Context section already knows the `--dont_gzip` fact but the implementation
outline does not carry it into the copy. The docs page qualifies both inline
("`--compression <1-9>` (FASTQ gzip level; ignored for uBAM output …)"; "`--dont_gzip` —
FASTQ output only … Rejected with `--output-format ubam`"); the help text should too.

### 1.4 Important — the losslessness paragraphs are FASTQ-scoped and become incomplete once uBAM is advertised

Paragraph 1's "byte-identically (header, sequence, quality)" and paragraph 3's
contract-scope note are written for a FASTQ-only mode. On the uBAM path:

- Aux tags round-trip **only** for tags named in `--preserve-tags`, and only A/Z/i/f
  scalars (`clump_only.rs:714-718`). A 10X user who runs `--clump_only --output-format
  ubam` on a `CB`/`UB`-carrying BAM **without** `--preserve-tags` silently loses them —
  while the help text calls the mode lossless.
- A `@PG` record is appended, so whole-file identity does not hold
  (`clump_only.rs:720-723`: "Cross-run byte-identity of the whole file is NOT guaranteed …
  but record-body byte-identity IS").
- "The plus-line (line 3 of each record) is normalized … inherited from the FASTQ
  reader/writer" is meaningless for BAM output.

"Lossless" is this mode's entire selling point, so the carve-out belongs in the help text,
not only on the docs page (which does cover it: "Byte-identity therefore covers the
semantically meaningful record fields, not the whole file"). Two clauses suffice — see §6.

### 1.5 Non-issues I checked so the implementer does not chase them

- **`main.rs:524`** — the code comment `// v1: FASTQ in/out. v2: uBAM in/out via
  --output-format ubam.` is *accurate* (it describes the version history of the dispatch it
  sits on) and is not user-facing. Leave it. Note the plan's grep pattern
  `"FASTQ in / FASTQ out"` (spaced slash) will not match it anyway; a looser pattern would.
- **`CHANGELOG.md:278-279`** and **`docs/src/content/docs/reference/changelog.md:28-29`**
  carry the same "FASTQ in/out only in v1; uBAM in/out is a natural follow-up" wording.
  These are **historical release notes** — true at the release they describe — and the docs
  changelog page is a synced artifact. They must **not** be edited. The plan's grep is
  scoped to `src/ tests/` so it will not surface them; keep it that way.
- **`cli.rs:10-11`** — the `OutputFormat` enum-level doc contains a
  `plans/…/PLAN.md §3.2` pointer. clap renders only *variant* docs under "Possible values",
  so this does not leak into `--help`. Verified in the rendered output. Not a defect.
- The `*_clumping_report.txt` claim holds on every arm, including the BAM ones
  (`clump_only.rs:836`, `:1082`).
- `--fastqc` genuinely composes on the BAM path (`clump_only.rs:845`, `:1090`).
- `--cores` genuinely affects the BAM path the same way it affects FASTQ (bin-layout sizing
  and FastQC threads), so its presence in "Composes with" is not falsified by uBAM.

---

## 2. Assumptions

**A1 — "no test pins the two stale phrases."** Holds, and the answer was knowable before
implementation, so deferring the grep costs nothing but gains nothing either:

- `tests/integration_no_args_help.rs:36` is the only help assertion in the tree, and it
  pins exactly two substrings: `"Usage: trim_galore"` and `"--adapter"`.
- No other file under `tests/` or `src/` asserts on `--help` output.
- `.github/workflows/release.yml:226` and `:356` run `trim_galore --help` as **exit-status**
  smoke tests with no content grep.

Conclusion: the doc block is unpinned prose. The plan's step-3 grep is safe to keep as a
belt-and-braces check, but should not be presented as the thing that de-risks the change.

**A2 — "the docs site's clump-only page is already accurate and needs no companion edit."**
Holds. `docs/src/content/docs/modes/clump-only.md` is comprehensive and correct on every
point above: per-format output filenames, the one-interleaved-BAM PE shape, `--compression`
ignored for uBAM, `--dont_gzip` rejected with uBAM, aux tags via `--preserve-tags` with the
A/Z/i/f constraint, the added `@PG` record, the three PE input-shape rejections, and the
FASTQ-output-with-uBAM-input rejection. **Use it as the source for the replacement copy** —
that is the cheapest path to an accurate help block and it keeps the two surfaces aligned.

**Unstated assumption worth surfacing:** the plan treats "the `--help` text" as one
surface. It is two (`-h` = paragraph 1, `--help` = all four), and which paragraph a fact
lives in decides who sees it. §1.1 is a direct consequence.

**Second unstated assumption:** that a "two-line fix" framing from the issue body should
bound the diff. The issue says "two-line fix" as a size estimate, not a scope contract; the
issue text itself came from an out-of-scope observation in another review and was never a
full audit of the block. Growing the scope to the whole doc block is the correct reading of
the issue's intent (make `--help` stop contradicting the shipped code), and the PR body
should say so in one sentence.

---

## 3. Efficiency

Nil, as the plan states. A doc-comment change has no runtime, allocation, or output-size
effect; the only cost is a recompile of the crate. Nothing further to review here.

---

## 4. Validation sufficiency

The three-row table is adequate for the change **as scoped**, but it cannot catch the
defects in §1.1–§1.4, and one of its rows has near-zero discriminating power.

- **Row 1 (grep for the stale phrases)** — will pass. Verifies removal, not correctness of
  the replacement. A grep cannot detect that the new sentence is also wrong (§1.2).
- **Row 2 ("new sentences present, old absent")** — the only row with real power, but it is
  stated as an eyeball check with no pass criteria for the rest of the block. As written, an
  implementer who reads only the two sentences they changed would sign it off with
  paragraph 1 still contradicting paragraph 2.
- **Row 3 (`cargo test` + fmt + clippy)** — **will pass no matter what the copy says.** No
  test asserts on this text, so this row proves the crate still compiles, nothing more. Per
  the project's own "prove the check can fail before trusting that it passed" rule, this
  should be labelled a build check, not a validation of the change.

**Strengthen to:**

1. Render **both** `-h` and `--help` from a freshly rebuilt binary (the plan's freshness
   note already covers the rebuild) and confirm the uBAM facts a user needs are present on
   the surface that user will hit.
2. Read the **entire four-paragraph block** side by side with the docs page's *Output
   filenames*, *Record fidelity*, and *Compatibility* sections. Pass criterion, stated
   explicitly: **every claim in the block is either true for both output formats or
   explicitly scoped to one.** That is the check that catches §1.1, §1.3 and §1.4.
3. Two smoke runs to confirm the copy against reality rather than against the docs — the
   fixtures already exist:
   - `--clump_only --output-format ubam test_files/BS-seq_10K_R1.fastq.gz` → confirm
     `BS-seq_10K_R1_clumped.bam` appears (paragraph 1's naming claim, SE).
   - `--clump_only --paired --output-format ubam test_files/BS-seq_10K_R{1,2}.fastq.gz` →
     confirm exactly **one** `_clumped.bam` (paragraph 1's naming claim, PE — the sharp case).
   - Optionally `--clump_only --dont_gzip --output-format ubam …` → confirm the rejection
     the new copy will describe.

**On pinning help text in a test:** I recommend **against** it. Prose assertions rot faster
than the prose they guard, and the mode already has an authoritative, maintained docs page.
The durable guard here is that the help text points at that page (the plan's copy already
does), not a substring assertion.

---

## 5. Alternatives

1. **Minimal, as planned (two sentences).** Rejected — leaves paragraph 1 false, leaves
   `-h` users uninformed, and ships the §1.2 inaccuracy.
2. **Make the whole four-paragraph block format-aware.** *Recommended.* Still one file, one
   doc comment, zero behaviour change, and it ends with a block that is internally
   consistent for both output formats. Cost: ~10 lines of prose instead of two, against an
   issue that estimated "two-line fix" — worth one sentence of explanation in the PR body.
3. **Shrink the help block and delegate the matrix to the docs page** — keep paragraph 1
   plus "Full compatibility matrix, output filenames and uBAM shapes: <docs URL>". Tempting:
   the block is already ~230 words of `--help`, and the docs page is authoritative. But
   `--help` should stand alone offline (cluster users without a browser), and this would
   *reduce* what `-h` users learn. Recommend the middle path in option 2: format-aware short
   lists in help, docs pointer retained for the full matrix.
4. **Fix only the docs, leave help alone.** Rejected — the false "rejected" listing actively
   deters users from a shipped feature, which is precisely the harm #400 was filed about.

---

## 6. Suggested copy (verified line by line against the code)

Offered because the plan flags exact wording as an open, reviewable question. Every claim
below was checked against the sources cited in §1.

**Paragraph 1** — add format-awareness and the BAM naming:

> Lossless reorder-only specialty mode: reorder FASTQ or unaligned BAM (uBAM) records by
> canonical 16-mer minimizer for compression-friendly grouping, WITHOUT any trimming,
> filtering, adapter detection, or record modification. Every input record appears in the
> output byte-identically — name, sequence, quality, plus aux tags named in
> `--preserve-tags` on uBAM; only file-level order changes. Output is `*_clumped.fq(.gz)`
> (SE) or `*_clumped_{1,2}.fq(.gz)` (PE), or under `--output-format ubam` a single
> `*_clumped.bam` (SE) and ONE interleaved `*_clumped.bam` per pair (PE). A short
> `*_clumping_report.txt` is emitted (distinct from `*_trimming_report.txt` to keep
> downstream nf-core/MultiQC scanners unconfused).

**Paragraph 2** — qualify the two format-dependent entries, move `--output-format ubam`
across, and keep the exclusion visible:

> Composes with `--compression` (FASTQ output only), `--memory`, `--cores`, `--paired`,
> `--fastqc`, `--dont_gzip` (FASTQ output only), `--basename`, `--output-format ubam`, and
> `--preserve-tags` (uBAM only). Trimming/filtering flags (`-a`, `--length`, `--rrbs`,
> `--polyA`, `--polyG`, `--trim-n`, `--clip_*`, `--nextseq`, `--rename`,
> `--discard_untrimmed`, `--consider_already_trimmed`, other specialty modes,
> `--passthrough`, `--retain_unpaired`) are rejected, as is `--dont_gzip` together with
> `--output-format ubam` (BAM is always BGZF). `-q` / `--stringency` / `-e` have clap
> defaults and are silently ignored on this path (mode does no trimming; matches how
> `--hardtrim5` treats trim flags today).

**Paragraph 3** — one added sentence for the BAM side of the contract:

> … Both normalizations are codebase-wide behaviours, inherited from the FASTQ
> reader/writer. On uBAM output a `@PG` record is appended to the header, so whole-file
> identity is not preserved — record bodies are.

**Paragraph 4** — replaces the stale line 206:

> FASTQ and uBAM input are both supported, but uBAM input requires `--output-format ubam`
> (the FASTQ output path would drop aux tags). Paired uBAM output takes either two FASTQ
> files or one interleaved uBAM; see the clump-only docs page for the full shape matrix.

---

## 7. CHANGELOG placement

`#### Changes` under `### Unreleased` is the right home — the heading exists
(`CHANGELOG.md:143`, inside the Unreleased span that runs to `:501`).

Two refinements to the plan's step 4:

- **Cite a better precedent.** `CHANGELOG.md:291-293` — *"**Tidied the `--help` text**:
  removed developer-internal references … No behaviour change."* — is in the **same
  section** and is the **same kind of change** (help text only, no behaviour). That is a
  stronger precedent than the `--output-dir` entry (`:168-169`), which is a runtime *stderr
  message*, not help text.
- **A separate bullet is still justified,** rather than extending that one. The tidy pass
  was about removing internal noise; #400 corrects a **false capability claim** that
  changes what a user believes they can run. Worth noting in passing that the tidy pass
  demonstrably *missed* line 206 — the same class of leftover version-roadmap note it
  claimed to have swept — which is mild support for auditing the whole block now rather
  than patching two lines.
- **Minor hazard:** `### Unreleased` contains both `#### Bug fixes` (`:6`) and `#### Fixes`
  (`:295`). That duplication is already tracked separately (#401) and is not this plan's
  problem, but the implementer should be careful not to land the entry in either of those
  by mistake. `#### Changes` is correct.

---

## 8. Action items

### Critical

1. **Bring paragraph 1 into scope and correct the output-naming sentence.** Under
   `--output-format ubam` the names are `<stem>_clumped.bam` (SE) and **one interleaved**
   `<stem>_clumped.bam` per pair (`io.rs:417-458`). This is the only paragraph `-h` renders,
   so leaving it FASTQ-only means the fix is invisible on the short-help path while
   actively misdirecting `--help` readers. (§1.1)
2. **Do not ship the proposed sentence "FASTQ and uBAM are both supported on input and
   output."** uBAM input on the FASTQ output path is rejected (`clump_only.rs:265-272`,
   `format.rs:165-172`). Say that uBAM input **requires** `--output-format ubam`. (§1.2)

### Important

3. **Qualify `--compression` and `--dont_gzip` when `--output-format ubam` joins the
   "Composes with" list.** `--dont_gzip` + ubam is a hard error (`cli.rs:576-581`);
   `--compression` is never passed to the BAM entry points and is recorded as level `0`
   (`clump_only.rs:780`, `:1018`). As drafted, the list would assert two mutually exclusive
   flags both compose. (§1.3)
4. **Scope the losslessness claim for uBAM:** aux tags survive only when named in
   `--preserve-tags` (A/Z/i/f only), and a `@PG` record is appended so whole-file identity
   does not hold. Silent `CB`/`UB` loss under a heading that says "lossless" is the
   user-visible risk. (§1.4)
5. **Replace validation row 3's framing and add the paragraph-level read-through.**
   `cargo test`/fmt/clippy pass regardless of what the copy says — no test asserts on this
   text. The load-bearing check is: render `-h` *and* `--help`, then read the whole block
   against the docs page with the explicit pass criterion "true for both formats, or
   scoped to one". Add the two `test_files/BS-seq_10K_R{1,2}.fastq.gz` uBAM smoke runs to
   confirm the naming claim rather than assuming it. (§4)

### Optional

6. **Fix the line references** in the plan: the rustdoc is `cli.rs:180-206` (not 183-207)
   and the "natural follow-up" sentence is line **206**. `cli.rs:195` is correct. (§1.5)
7. **Adjacent, same bug class — decide in or out explicitly.** The `OutputFormat::UBam`
   variant doc (`cli.rs:17-18`) lists five of the seven ubam rejections; `--retain_unpaired`
   (`cli.rs:609`) and `--dont_gzip` (`cli.rs:576`) are missing, and this text *does* render
   in `--help` under "Possible values". One line to fix in the same PR, or file separately —
   but make it a decision, not an omission.
8. **Warn the implementer off the historical hits.** `CHANGELOG.md:278-279` and
   `docs/src/content/docs/reference/changelog.md:28-29` contain the same stale-sounding
   wording as release history and must not be edited; `main.rs:524`'s comment is accurate.
   Keep the grep scoped to `src/ tests/`. (§1.5)
9. **Cite `CHANGELOG.md:291-293` ("Tidied the `--help` text") as the precedent** instead of
   the `--output-dir` message entry — same section, same kind of change. (§7)
10. **Add one PR-body sentence** explaining why a "two-line fix" issue produced a
    ~10-line doc diff: the sibling claims live in the same doc block, and fixing two of four
    paragraphs would leave the block self-contradictory. (§2)
