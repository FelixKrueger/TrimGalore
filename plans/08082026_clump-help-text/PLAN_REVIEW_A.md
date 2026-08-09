# Plan Review A — `--clump_only` help text (#400)

**Plan:** `plans/08082026_clump-help-text/PLAN.md`
**Reviewer:** A (independent)
**Base:** `dev` @ `ce13761`
**Verdict:** REVISE before implementing. Premise correct, every mechanical claim verified true — but the fix is incomplete in the same paragraph, and both proposed sentences introduce new inaccuracies.

---

## 1. Logic review

### 1.1 What I verified as correct

Everything the plan asserts about the *world* checks out, with one off-by-one:

| Plan claim | Verified | Evidence |
|---|---|---|
| `--output-format ubam` wrongly in the rejected list | **True** | `src/cli.rs:195`. `Cli::validate` has no ubam rejection under `clump_only`; instead `src/cli.rs:866-869` carries an explicit comment: *"v2: `--output-format ubam` is now accepted; uBAM in/out is dispatched to the BAM variant functions in main.rs."* |
| "v1 is FASTQ in / FASTQ out only. uBAM in/out is a natural follow-up." is stale | **True** | `src/cli.rs:206` (not 207 — see 1.3). Five shipped arms in `src/main.rs:522-750`. |
| A2: docs page already accurate, no companion edit | **True** | I read all 138 lines of `docs/src/content/docs/modes/clump-only.md`. It is accurate against the code on every point I checked, including the PE shapes, the `--dont_gzip` split, `--compression` inertness under uBAM, and the `@PG` caveat. **It is the correct reference truth for the rewrite** — see 5.2. |
| A1: no test pins the phrases | **True, and knowable now** | See 3. |
| CHANGELOG `#### Changes` under `### Unreleased` | **True** | `CHANGELOG.md:143`. Precedent is stronger than the plan claims — see 4. |

So the plan is not wrong about anything it looked at. The problem is what it did not look at.

### 1.2 The fix as written ships two *new* false statements

**(a) `--dont_gzip` + `--output-format ubam` — the plan's own Context knows this and no step acts on it.**

The sentence being edited (`cli.rs:190-191`) currently reads:

> Composes with `--compression`, `--memory`, `--cores`, `--paired`, `--fastqc`, `--dont_gzip`, and `--basename`.

Step 1 says to *extend this same sentence* with `--output-format ubam` and `--preserve-tags`. The result is one sentence asserting that `--dont_gzip` and `--output-format ubam` both compose — while `src/cli.rs:576-581` bails on exactly that pair:

```
--dont_gzip is not compatible with --output-format ubam (BAM is always BGZF-compressed)
```

Plan Context line 14 records this fact verbatim ("`--dont_gzip` is rejected *with* uBAM output") but neither implementation step consumes it. For a PR whose entire purpose is deleting a false claim from this paragraph, adding one back into the same sentence is the worst available outcome. The docs page already has the wording to copy (`clump-only.md:78`, `:80`, `:95`).

**(b) The replacement copy overstates the format matrix.**

Step 2's proposed sentence — "FASTQ and uBAM are both supported on input and output" — reads as 2×2. Only 3 of the 4 combinations are legal:

| | FASTQ out | uBAM out |
|---|---|---|
| FASTQ in | ✅ | ✅ |
| uBAM in | ❌ **rejected** | ✅ |

uBAM-in → FASTQ-out is refused (`src/main.rs:537-544` for the `--paired` N=1 shape, plus the SE path pinned by `tests/integration_clump_only.rs:304` `rejects_ubam_input_without_ubam_output`, which asserts the error routes the user to `--output-format ubam`). Rationale in both places: FASTQ output would silently drop aux tags. `clump-only.md:103` states this explicitly. The replacement sentence needs to say uBAM input *requires* `--output-format ubam`, otherwise it invites the exact invocation the code rejects.

### 1.3 Incomplete defect set — the same paragraph has more uBAM-falsified claims

The prompt's worry was justified. #400's own criterion is "predates the shipped uBAM arms and now contradicts the docs site." Applying that criterion to the rest of the doc block turns up three more hits, one of them more consequential than either the plan names.

**(i) The output-filename sentence — `cli.rs:184-188`. Highest user impact of anything in this review.**

> Output files use the `*_clumped.fq(.gz)` (SE) or `*_clumped_{1,2}.fq(.gz)` (PE) suffix, and a short `*_clumping_report.txt` is emitted …

Unqualified, and false under `--output-format ubam`:

- SE uBAM → `*_clumped.bam` (`src/io.rs:417-431`)
- PE uBAM → **ONE** interleaved `*_clumped.bam`, no `_1`/`_2` suffix at all (`src/io.rs:442-458`, whose own doc comment says *"Returns a SINGLE interleaved BAM path … no `_1`/`_2` suffix"*)

The PE case differs in output **count**, not just extension. Users script against filenames; this is worse than a stale forward-looking sentence, and it is two lines above one of the two the plan is fixing. `clump-only.md:41-45` has the full correct list.

**(ii) The contract-scope paragraph is FASTQ-only but reads as exhaustive — `cli.rs:200-204`.**

> … The plus-line (line 3 of each record) is normalized to bare `+` on output; CRLF line endings are normalized to LF. **Both** normalizations are codebase-wide behaviours, inherited from the FASTQ reader/writer.

On the uBAM path there is no plus-line, and there is a *third* deviation from whole-file identity that this paragraph does not mention: every run appends a `@PG ID:trim_galore VN:<version> CL:<invocation>` record (`&command_line` threaded into `clump_only_single_to_bam` / `clump_only_paired_to_bam_one_pair` at `src/main.rs:645`, `:702`, `:740`; documented at `clump-only.md:61-63` and `CHANGELOG.md:226-231`). Whole-file BAM hashes therefore vary per invocation — which is why CI uses a `@PG`-ignoring assertion.

This matters *because* of step 2: planting a "uBAM is supported" sentence immediately below a paragraph whose closing word is "**Both**" tells the reader the exhaustive caveat list covers uBAM. It doesn't.

**(iii) The byte-identity enumeration omits aux tags — `cli.rs:183-184.`**

> Every input record appears in the output byte-identically (header, sequence, quality)

For a uBAM record those are not all the fields. Aux tags round-trip **only** with `--preserve-tags` and are otherwise dropped. Compare `clump-only.md:6`: "(header/name + sequence + quality + preserved aux tags for uBAM)". Severity is lower than (i)/(ii) because the enumeration doesn't claim tags are kept — and step 1's `--preserve-tags` addition largely mitigates it *if* worded as an opt-in rather than a bare list item.

### 1.4 Claims I checked that are fine — no edit needed

Recording these so the implementer doesn't over-reach:

- **`--fastqc`** — genuinely composes on every arm including BAM (`main.rs:646`, `:703`, `:741`); fastqc-rust reads BAM natively (`clump-only.md:77`). Accurate.
- **`--basename`** — threaded through all five arms. Accurate.
- **`-q` / `--stringency` / `-e` silently ignored** (`cli.rs:195-198`) — matches `validate`'s comment at `cli.rs:738-741` and `clump-only.md:105`. Accurate.
- **`--output-format`'s own help** (`cli.rs:243-249`) — says "some flag combinations are rejected (see `--help` and the startup diagnostics for details)". Still true; no collateral edit required. Good news for scope.
- **"other specialty modes" in the rejected list** — loosely covers `--clumpify`/`--hardtrim*`/`--clock`/`--implicon`/`--demux`, all of which do bail (`cli.rs:743`, `851-865`). Acceptable as-is.

### 1.5 Line-number drift in the plan

The plan says the block is `cli.rs:183-207` and cites `cli.rs:207` twice for the follow-up sentence. Actual: the doc block spans **180-206**; the stale sentence is at **206**; **207 is the `#[clap(long = "clump_only")]` attribute**. Trivial, but an implementer editing "line 207" touches the attribute. Line 195 is correct as cited.

---

## 2. Assumptions

- **A1 (no test pins the phrases) — verified TRUE. The plan defers this to implementation; it is already answerable, so the deferral buys nothing.** See 3 for the evidence.
- **A2 (docs page accurate post-#396) — verified TRUE**, and stronger than the plan uses it for. The plan treats the docs page as "needs no companion edit"; it should also treat it as the *specification* for the rewrite (5.2).
- **Unstated assumption: "two statements" is the complete defect set.** This is the plan's load-bearing implicit assumption and it is **false** (1.3). It comes from taking the issue body's "two-line fix" framing as a scope finding rather than as the reporter's first-glance estimate.
- **Unstated assumption: the terseness convention applies here.** The plan never cites it, but its "Fix the two statements; change nothing else" instinct pattern-matches to the project's one-line-comment rule. That rule governs **code comments**; clap rustdoc is user-facing `--help` prose where completeness is the governing virtue. The project demonstrates this itself — 138 lines of docs prose for this one flag. The terseness rule must not be used to justify leaving 1.3(i) in place.

---

## 3. Test / CI pinning — safe, and the answer is knowable now

The plan defers the grep to implementation step 3. **Safe to defer, but unnecessary — here is the answer.** I checked wider than the plan's `tests/ src/` scope:

- **The two phrases:** `src/cli.rs:206` is the only source hit. The other two hits are **historical release records** — `CHANGELOG.md:279` and its synced mirror `docs/src/content/docs/reference/changelog.md:29`, both inside the shipped `--clump_only` feature entry. **Both must be left alone** (a changelog records what was true at release; the synced page is additionally off-limits by project convention). The plan's grep scope (`tests/ src/`) would never surface them — fine for the "does a test pin this" question, but the plan should say so, or a diligent implementer greps wider, finds two hits, and either panics or wrongly "fixes" the changelog.
- **`tests/integration_no_args_help.rs`** pins only `"Usage: trim_galore"`, `"--adapter"`, and the absence of `"error:"` (lines 35-42). Nothing from the clump block.
- **No test anywhere renders or asserts long help** — no `render_long_help` / `long_about` / `verbatim_doc_comment` hits in `src/`, no snapshot tests.
- **CI runs `--help` as an exit-code smoke test only** — `.github/workflows/release.yml:226` and `:356`. No content grep.
- **No `cargo doc` job in `.github/workflows/ci.yml`** and **no `[lints]` section in `Cargo.toml`** — so rustdoc lints (including `bare_urls`) are not a gate. Relevant to 6.4.

Conclusion: **zero pinning risk.** The help text is completely unguarded, which is also why it drifted.

---

## 4. CHANGELOG placement — correct, with better precedent than cited

`#### Changes` exists under `### Unreleased` at `CHANGELOG.md:143`. The plan cites the `--output-dir` entry as precedent; that exists (`CHANGELOG.md:168`, "Two collision messages advised `--output-dir`, which is not a valid flag"). Two stronger ones sit in the same section:

- `CHANGELOG.md:145` — "**The `--passthrough`-matches-R1/R2 message now says what it checks**" (#389, the immediately preceding commit `0f65048`). Message-only, same section. Closest structural match.
- `CHANGELOG.md:291-293` — "**Tidied the `--help` text**: removed developer-internal references … that had leaked into user-facing flag descriptions. No behaviour change." A pure `--help`-text entry under `#### Changes`. Direct precedent.

Placement approved. One content note: if the fix expands per 1.3, the entry should name the filename correction explicitly — a user who scripted against `_clumped_{1,2}.fq.gz` for uBAM output was misled, and that is the line most worth their attention.

---

## 5. Validation sufficiency

### 5.1 The three proposed checks cannot catch the failure mode that matters

Validation #1 greps for the removed phrases, #2 eyeballs the rendered help, #3 runs the suite. All three confirm the *mechanics* — that something changed and nothing broke. **None of them checks that the new sentences are true.** Given that 1.2 identifies two false statements the plan's own copy would introduce, the validation table has a hole exactly where the risk is.

The rebuild-before-render discipline (#2, called out in Self-Review against the session's prove-freshness rule) is genuinely good and worth keeping.

### 5.2 Two rows to add

- **#4 — Accuracy diff against the maintained truth.** Read the rewritten paragraph line-by-line against `docs/src/content/docs/modes/clump-only.md` §Compatibility (lines 69-105) + §Output filenames (39-45) + §Record fidelity (47-63). That page is accurate (verified in 1.1) and is the specification. Expected: every claim in the help text has a corresponding, non-contradicting line on the docs page.
- **#5 — No new claim contradicts a `bail!`.** For each flag named as composing, confirm no matching `anyhow::bail!` exists in either the `if self.clump_only` block (`cli.rs:742-883`) or the shared §3.4a `UBam` block (`cli.rs:571-616`). This is the check that catches 1.2(a) mechanically; it takes about two minutes and would have caught the defect at plan time.

### 5.3 Do not add a help-text regression test

Worth stating so it isn't proposed later: pinning help prose in a test is brittle and the repo has deliberately zero help-content assertions. The root cause here is doc-drift-after-feature, which no unit test catches. The durable mitigation is the single-source-of-truth question in 6.1, not an assertion.

---

## 6. Alternatives

### 6.1 Recommended: repair the paragraph, don't patch two sentences

Three options, in increasing scope:

- **(A) As planned — two statements.** Leaves 1.3(i) (wrong output filenames for uBAM) and 1.3(ii) (`@PG`) in place, and per 1.2 introduces two new inaccuracies. A second #400-shaped issue against the same paragraph becomes near-certain. Not recommended.
- **(B) Repair the whole block against the docs page. ← recommended.** Still small — roughly 8-10 lines of doc comment. Fixes the whole class in one PR, one CHANGELOG line, same review cost. Keeps the paragraph self-sufficient for offline `--help` readers.
- **(C) Shrink the help block and delegate the matrix to the docs page.** Structurally attractive: one maintained source, drift impossible by construction. But it removes information from `--help` for users without network access, is a larger editorial change than #400 asked for, and is Felix's call on register rather than a reviewer's. Worth raising as a follow-up question, not doing here.

### 6.2 On the docs-page cross-reference in step 2

Step 2's parenthetical — "(see the clump-only docs page for the paired input shapes)" — would be **the first docs-site cross-reference in any `cli.rs` help string** (zero hits for `docs page` / `documentation` / `trimgalore.com` across `cli.rs`). And the nearest precedent runs the other way: `CHANGELOG.md:291-293` records deliberately *removing* unresolvable pointers (internal `PLAN.md`/`§` references) from user-facing flag descriptions. "The clump-only docs page", with no URL, is closer to that pattern than to a resolvable reference.

There are only two paired shapes. Stating them inline costs one clause each and keeps `--help` self-contained:

> `--paired` takes two files (Shape A, multi-pair supported) or a single interleaved uBAM (Shape B, requires `--output-format ubam`).

If a pointer is preferred anyway, name the URL. No `cargo doc` job means a bare URL won't trip `rustdoc::bare_urls` (see 3), but backticks remain the local convention.

### 6.3 Efficiency

Nil, as the plan says. Doc-string only; no runtime, no allocation, no binary-size change worth measuring. Nothing to review here.

---

## 7. Action items

### Critical — fix before implementing

1. **Do not add `--output-format ubam` to the existing "Composes with" list while `--dont_gzip` sits in it** (`cli.rs:190-191`). `cli.rs:576-581` rejects that pair. Split the sentence, or qualify `--dont_gzip` as FASTQ-output-only. Model wording: `clump-only.md:78`, `:80`. The plan's Context already knows this fact — promote it into an implementation step.
2. **Fix the replacement copy's format matrix** (step 2). "FASTQ and uBAM are both supported on input and output" implies 4 legal combinations; uBAM-in → FASTQ-out is rejected (`main.rs:537-544`; `tests/integration_clump_only.rs:304`). Say that uBAM input requires `--output-format ubam`.

### Important

3. **Extend the fix to the output-filename sentence** (`cli.rs:184-188`). Under `--output-format ubam` the outputs are `*_clumped.bam` (SE) and **one** interleaved `*_clumped.bam` for PE — different count, not just different extension (`io.rs:417-458`). By #400's own stated criterion this is the same defect class, and it has the highest user impact of anything in the block. Fix here, or open a follow-up issue explicitly — do not skip silently.
4. **Add a uBAM clause to the contract-scope paragraph** (`cli.rs:200-204`). It closes on "**Both** normalizations", which reads as exhaustive, but the uBAM path adds a per-run `@PG` record so whole-file identity does not hold (`clump-only.md:61-63`). Step 2 plants a uBAM sentence directly beneath it, which makes the over-claim worse than it is today.
5. **Add validation rows #4 and #5** (5.2): diff the rewritten paragraph against `clump-only.md`, and confirm no newly-claimed flag has a matching `bail!` in `cli.rs:742-883` or `cli.rs:571-616`.

### Optional

6. **Qualify `--compression` as inert under uBAM output.** The `clump_only_*_to_bam` call sites take no compression argument (`main.rs:638-650`, `:695-707`, `:733+`), unlike the FASTQ arms which pass `cli.compression`. `clump-only.md:74` says "ignored for uBAM output". Nearly free while editing that sentence.
7. **Word `--preserve-tags` as an opt-in, not a bare list item** — "aux tags round-trip only with `--preserve-tags`" — which also covers 1.3(iii)'s tag-loss gap.
8. **Note the report count while touching that sentence**: PE FASTQ writes one report per mate, the uBAM shapes one per pair (`clump-only.md:122`, and #391's entry at `CHANGELOG.md:6-18`).
9. **Prefer inline shapes over "see the clump-only docs page"** (6.2), or name the URL.
10. **Correct the plan's line numbers**: block is 180-206, stale sentence at 206; 207 is the `#[clap]` attribute.
11. **Record in the plan that `--cores` in the composes-list is accepted-for-parity only** (single-threaded internally per `cli.rs:749-752` and `clump-only.md:76`). Pre-existing and *not* uBAM-caused, so genuinely out of #400's scope — flagged only so the decision to leave it is deliberate.
12. **State in step 3 that the two non-`src/` hits are historical changelog records** (`CHANGELOG.md:279`, `docs/…/reference/changelog.md:29`) and are intentionally not edited.

### Out of scope — noted, not for this PR

- `--phred64` + FASTQ-in + uBAM-out is accepted and threads `cli.phred_offset()` into the BAM writers (`main.rs:649`, `:706`, `:744`), whereas the FASTQ→FASTQ arm ignores phred entirely (byte copy). Quality *values* are preserved; the ASCII string necessarily is not, because BAM stores raw Phred. Semantically correct, asymmetric with the FASTQ arm's literal byte-copy, and immaterial to `--help` accuracy. `--phred64` + uBAM **input** is already rejected (#358, `CHANGELOG.md:235`). Mentioning only so it isn't rediscovered as a bug.

---

## 8. Verdict

The plan's premise is right and its homework is honest: both named defects are real, and every checkable claim it makes — A1, A2, the CHANGELOG heading, the docs page's accuracy — verified true, which is more than most plans this size manage. But it inherited the issue's "two-line fix" framing as a scope conclusion, and that framing does not survive reading the paragraph. Three more claims in the same doc block were falsified by the shipped uBAM arms, the worst of them telling users the wrong output filenames and the wrong output *count* for PE uBAM; and the two sentences the plan does propose would each introduce a fresh inaccuracy — one asserting `--dont_gzip` composes with `--output-format ubam` when `cli.rs:576` rejects exactly that, the other implying uBAM-in → FASTQ-out works when it is refused with a message telling users to add `--output-format ubam`. Fix the two Critical items, take the two Important ones (or defer item 3 to a named follow-up rather than silently), add the two accuracy validation rows, and this becomes the right change — still under a dozen lines of doc comment, with `docs/src/content/docs/modes/clump-only.md` serving as a specification that is already correct.
