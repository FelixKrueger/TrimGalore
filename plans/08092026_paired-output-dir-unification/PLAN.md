# Plan: one output directory per pair — reports and the passthrough carrier follow R1 (#398)

**Issue:** [#398](https://github.com/FelixKrueger/TrimGalore/issues/398) — follow-up from #391's plan review (B §1.3).
**Decisions (Felix):** reports follow the primaries into R1's parent; **`--passthrough`'s carrier output joins them**; the new report collisions are **accepted as loud pre-write refusals**, reconfirmed after being shown the mate-side-directory layout.
**Base:** `dev` @ `38deb1ee`.

## Revision history

- **r2 (2026-08-10):** dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — both REVISE, six Critical findings between them, no contradictions). Both independently verified the extraction is behaviour-preserving across all five namers, the arm table is complete, and A3 holds empirically. Six things changed. **The behavioural footprint is bidirectional**: a refusal also *disappears* (A §1.3), which r1 was silent about. **A second inverted guard exists** on the clump arm (A §1.2) — r1 built a whole warning around one of the pair and left the twin unannounced. **A new output-vs-input hazard** in #409's exact shape (A §1.4), which no r1 validation targeted. **Step 7 as written broke two single-end arms** r1 declared out of scope (B §1.8). The **refused population** is broader than "two mates of one pair" (both reviewers). And r1's **step-6 rationale for `--passthrough` was factually false** (both reviewers).
- r1 (2026-08-09): first draft. Written after two critical scope questions were resolved by Felix.

### Corrections to the reviews

- Both reviewers instruct changing the docs path prefix to `Docs/src/content/docs/`. **Both are wrong.** `git ls-files` tracks it as lowercase `docs/src/content/docs/guide/outputs.md`; the capital resolves only because APFS is case-insensitive, and `Docs/` would not resolve on a Linux checkout. r1's shorthand stands.
- A §1.7 asks that `pair_output_dir` "return the empty path rather than `.`". Correct in effect, but the safe instruction is narrower: copy the existing expression **verbatim**. See §Signature.

## Goal

When `--output_dir` is not given, every output of a paired run lands in one directory — R1's parent — instead of up to three. Today the primaries go to R1's parent, each trimming/clumping report goes beside its own mate, and `--passthrough`'s carrier output goes beside the passthrough input.

Single-end is untouched, paths included.

**Perl precedent, which r1 lacked.** v0.6.11 had no per-mate report directories at all: `trim{}` does `my $output_filename = (split (/\//,$filename))[-1];` — stripping the input's directory — then opens the report at `$output_dir.$report`, where `$output_dir` defaults to `''`, i.e. the CWD (verified at `git show 0.6.11:trim_galore`, the `sub trim` header and the `### OUTPUT DIR PATH` block). So Perl put every output of a run in **one** directory. The split layout is a v2 regression against v0.6.x behaviour, not a v0.6.x inheritance — which strengthens decision 1 and constrains what the changelog may claim (see §Integration).

## Context

### The asymmetry, precisely

Every **primary** paired namer in `io.rs` resolves its directory with a character-identical expression — `output_dir`, else `input_r1.parent()`, else `.`:

| Namer | Directory expression |
|---|---|
| `paired_bam_output_name` | `io.rs:296-298` |
| `paired_end_output_names` | `io.rs:337-339` |
| `unpaired_output_names` | `io.rs:368-370` |
| `clumped_paired_output_names` | `io.rs:473-475` |
| `clumped_paired_bam_output_name` | `io.rs:522-524` |

The **report** namers do not. `report_name` (`io.rs:551`), `json_report_name` (`io.rs:567`) and `clumping_report_name` (`io.rs:535`) each take one `input` and fall back to *that* input's parent. They are shared with the single-end paths, where that is correct — so the fix cannot live inside them. `passthrough_output_name` (`io.rs:404-409`) is a third rule, anchored on the carrier's own parent, while `io.rs:385` already calls it one of "the three pair outputs".

Reproduced on `38deb1ee`:

```
$ trim_galore --paired A/reads.fq B/reads.fq      # exit 0
A/  reads_val_1.fq  reads_val_2.fq  reads.fq_trimming_report.{txt,json}
B/                                  reads.fq_trimming_report.{txt,json}
```

### The behavioural footprint is bidirectional — r1 got this wrong

r1 described exactly one new refusal and was silent on refusals disappearing. Three classes of test change, and **each is a red or newly-green test**. The plan's own thesis is that unpredicted red tests are how this bug family recurs, so all three are enumerated here rather than left to be discovered.

**(a) Two inverted over-rejection guards, not one.** r1 named the trim guard (`tests/integration_output_collision.rs:1029-1041`) and missed its clump twin at **`:1206-1240`** (`clump_report_candidates_do_not_over_reject`), whose comment carries almost the same sentence:

> `Pins that the candidates honour output_dir = None; a builder that resolved reports into one directory would over-reject this.`

Both were written to catch precisely this change. Both must be **rewritten to assert the refusal, never narrowed** — narrowing restores the split layout on the candidate side only, while the writers move. The trim guard is three blocks and only block 3 changes; the clump twin is a single block, so its whole body is rewritten.

**(b) An existing deliberate refusal disappears.** `clump_paired_rejects_shared_mate_report_without_output_dir` (`:1247-1275`) covers `--clump_only --paired p/a.fq d/x.fq q/b.fq d/x.fq`, refused today because both pairs' reports for the shared mate `d/x.fq` resolve to one path. After the change pair 1 anchors on `p/` and pair 2 on `q/`, giving eight distinct paths, and **the run succeeds**. That is correct on the facts — two distinct reordered outputs in two distinct directories, nothing overwritten — so continuing to refuse would be over-rejection. The test must flip to asserting success. **The trim-mode twin of the same shape has no test at all** and would flip refuse→succeed silently; add one.

**(c) A new output-vs-input hazard, in #409's exact shape.** Moving a report into R1's parent can make it land *on an input*. Confirmed on `38deb1ee`:

```
$ trim_galore --paired A/reads.fq_trimming_report.txt B/reads.fq   # exit 0 today
A/  reads.fq_trimming_report.txt                      <- R1's input, intact
B/  reads.fq_trimming_report.{txt,json}               <- R2's report, harmless here
```

R2's report is `B/reads.fq_trimming_report.txt` today. After the change it becomes `A/reads.fq_trimming_report.txt` — **R1's own input path**. If the candidate builder moves with the writer, the pre-flight's output-vs-input branch (`io.rs:147-156`) refuses, which is correct. If it does not, the writer destroys a named input after the run has read it. That is #409 verbatim, on an arm this change creates. It needs its own test per paired report arm, in both the refusing and the `--no_report_file` accepting direction — the shapes at `tests/integration_output_collision.rs:1411` and `:1440` already model both.

### The refused population, stated correctly

r1 said "two mates share a filename", which is too narrow. The actual rule: **any two report-keyed inputs whose filenames fold-equal (`collision_key` is ASCII-case-folded) and whose pairs anchor on the same directory.** That spans pairs, not just mates within a pair.

The practically important instance is the **mate-side directory layout**, where a facility splits by mate rather than by sample. Confirmed running today:

```
$ trim_galore --paired R1/s1.fq R2/s1.fq R1/s2.fq R2/s2.fq        # exit 0
R1/  s1_val_{1,2}.fq  s2_val_{1,2}.fq  s1.fq_trimming_report.*  s2.fq_trimming_report.*
R2/  s1.fq_trimming_report.*  s2.fq_trimming_report.*
```

After the change the **entire invocation** is refused before any pair is processed — `main.rs:863-914` builds candidates for every chunk and calls the pre-flight once at `:915`, so this is not "pair 1 fails, pair 2 proceeds". Nothing runs.

The **per-sample** layout is unaffected, and saying so is worth as much as the warning: `sampleA/reads_1.fq sampleA/reads_2.fq sampleB/reads_1.fq sampleB/reads_2.fq` gives each pair its own R1 directory, so nothing collides.

### `OutputSource` attribution — a decision #397 explicitly deferred to #398

`main.rs:872-875`, written by #397:

> `Uniform Pair for paired primaries (#397 decision C2): _val_1 takes its stem from R1 but its directory from R1 too, and _val_2 mixes both — naming one mate would encode a claim about which supplies the directory, which is exactly what #398 may change.`

After this change the report and carrier candidates take their **filename** from one input and their **directory** from R1 — the condition `OutputSource::Pair` exists for (`io.rs:94-96`). **Decision: keep `OutputSource::Input`,** because switching degrades a message that shipped last week:

- `Input` → two distinct sources on one path → the **2a** branch (`io.rs:171-180`): *"the outputs named from A/reads.fq and B/reads.fq would be written to the same file"*. Actionable: rename one input.
- `Pair(A,B)` → both candidates share one source identity (`io.rs:103-111`) → the **2c** branch: *"List each input once — a file that appears in more than one pair is reported once per pair"*. That advice is **wrong here**: there are two inputs and renaming one does fix it.

### Writers and their pre-flight twins must move together

The pre-flight is only as complete as its candidate list, and this change alters *paths*, so every candidate builder moves in lockstep with its writer. Five recurrences of that family are fixed (#383, #388, #391, #409); #414 exists to make a sixth un-shippable.

| Arm | Writer | Candidates | Net effect |
|---|---|---|---|
| Trim paired → FASTQ | `main.rs:1729-1738` | `main.rs:906-912`; carrier `:893-903` | **changes** |
| Trim paired → uBAM (two FASTQ in) | `main.rs:2352-2358` | `main.rs:2026-2027` | **changes** |
| Clump paired → FASTQ | `clump_only.rs:535-536` | closure at `main.rs:609-622` | **changes** |
| Clump paired → BAM | `clump_only.rs:1080-1086` | `main.rs:698`, `:740` | **no-op** — `report_input = &inputs[0]`/`chunk[0]` is already R1, so both sides are behaviour-identical. Change it for uniformity, expect no path to move |

Three arms change behaviour; the fourth is a no-op recorded so nobody hunts for a diff that isn't there.

Out of scope, each by inspection:
- **Single-end arms** — one input, no second parent to disagree with.
- **Paired single-interleaved uBAM** — one input path, so the asymmetry is inexpressible.
- **`--clock` / `--implicon`** — write no reports (`grep -c report src/specialty.rs` → 0) **and**, load-bearingly, their primaries are CWD-anchored (`specialty.rs:488-491`, `:521-524` return a bare filename when `output_dir` is `None`), so R1-anchoring would *change* them. r1 gave only the weaker reason.
- **`--hardtrim5/3`** — the same CWD rule, with its own `CWD_OUTPUT_HINT`.
- **`--demux`** — single-end only, enforced at `cli.rs:1009-1012`.
- **`--clumpify`** (vs `--clump_only`) — writes no separate report; routes through the ordinary trim report path, so it is covered by the trim arms.
- **`--retain_unpaired`** — `unpaired_output_names` is already R1-anchored (`io.rs:368-370`).
- **`--fastqc`** — anchored on the *output* file (`fastqc.rs:37-58` passes `output_path` through), so its artifacts follow the primaries and the carrier for free. `--fastqc_args "-o DIR"` is already out of scope at `io.rs:1338-1340`.

## Behavior

1. One derivation decides the directory for every output of a paired run: `--output_dir` when given, else R1's parent.
2. The five primary namers consume it. **No primary path changes** — pure extraction of an expression they already share.
3. Both trimming reports of a pair (`.txt` and `.json`) resolve into that directory. Report **filenames** are unchanged — still each mate's own full input filename.
4. Both clumping reports of a `--clump_only --paired` run likewise.
5. `--passthrough`'s carrier output likewise. Filename unchanged.
6. Single-end behaviour is byte-identical, paths included.
7. **Refusals gained:** any two report-keyed inputs whose filenames fold-equal and whose pairs anchor on the same directory now collide, and the run is refused before anything is written, naming both inputs via the 2a branch. Includes the cross-pair mate-side-directory layout, which refuses the whole invocation.
8. **Refusals lost:** a mate shared between two pairs whose R1s live in different directories no longer collides, because the two reports now resolve into different directories. The run succeeds. Correct on the facts; nothing is overwritten.
9. **Output-vs-input:** a relocated report may now equal a named input. The pre-flight refuses it — which requires the candidate to have moved with the writer.
10. Every candidate builder for the arms above produces exactly the paths its writer will write, with the same `--no_report_file` gating.

### Edge cases

| Case | Expected |
|---|---|
| `--paired A/r1.fq A/r2.fq` (same dir) | unchanged |
| `--paired A/reads.fq B/reads.fq` | **refused** (was: exit 0, split reports) |
| `--paired R1/s1.fq R2/s1.fq R1/s2.fq R2/s2.fq` (mate-side layout) | **whole invocation refused**, no pair processed |
| `sampleA/r_1 sampleA/r_2 sampleB/r_1 sampleB/r_2` (per-sample layout) | unchanged, each pair anchors on its own R1 |
| `--paired A/reads.fq B/reads.fq --no_report_file` | succeeds; only reports collided |
| `--paired A/reads.fq B/reads.fq -o out` | unchanged — already refused today |
| `--paired A/r1.fq B/r2.fq --basename foo` | succeeds; reports keep per-input names, now both in `A/` |
| `--clump_only --paired p/a.fq d/x.fq q/b.fq d/x.fq` | **now succeeds** (was refused) |
| `--paired p/a.fq d/x.fq q/b.fq d/x.fq` (trim twin) | **now succeeds** (was refused, untested) |
| `--paired A/reads.fq_trimming_report.txt B/reads.fq` | **refused** — R2's report would land on R1's input (was: exit 0) |
| same, `--no_report_file` | succeeds; no report to collide |
| R1 as a bare filename (`reads.fq`) | `parent()` → `Some("")`, so paths stay bare; **stderr spelling must not gain a `./` prefix** |
| Multi-pair, distinct dirs | each pair anchors on **its own** R1 |

## Signature

```rust
/// The directory every output of a paired run lands in: `--output_dir` when
/// given, else R1's parent, else the current directory.
///
/// One derivation for primaries, reports and the `--passthrough` carrier alike.
/// #398 existed because the primaries used this rule while the reports used each
/// mate's own parent, so a pair could scatter across three directories.
pub fn pair_output_dir(input_r1: &Path, output_dir: Option<&Path>) -> PathBuf
```

**Copy the body verbatim from any of the five namers** — `output_dir.map(|d| d.to_path_buf()).unwrap_or_else(|| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf())`. For a bare filename `parent()` yields `Some("")`, so the `unwrap_or(".")` never fires and joined paths stay bare. "Simplifying" it to anything that yields `"."` would prefix every path and every error message with `./`, which is a user-visible diff in stderr. The `unwrap_or(".")` arm is reachable only for a root path.

Reports are then named by passing that directory to the existing namers — `report_name(input_r2, Some(&pair_dir))` — which keeps each report's per-input filename while relocating it. No new report namer and no signature change to the three shared report functions.

## Implementation outline

1. **`io.rs`** — add `pair_output_dir` beside `collision_key`, body copied verbatim.
2. **`io.rs`** — refactor the five primary namers to call it. Their existing unit tests must pass untouched.
3. **`main.rs` trim-paired writer** (`:1729-1738`) — compute `let pair_dir = naming::pair_output_dir(input_r1, output_dir);` once, then build both `PairedReportFile`s with `Some(&pair_dir)`.
4. **`main.rs` trim-paired uBAM writer** (`:2352-2358`) — same.
5. **`clump_only.rs` paired writers** (`:535-536`; `:1080-1086`'s `report_input` for uniformity, a no-op) — same.
6. **`--passthrough`** — fold `input_r1` into `passthrough_output_name` and derive via `pair_output_dir`. **r1's rationale for the alternative was false**: the namer already takes `output_dir`, and there are no single-end-shaped callers (`--passthrough` is paired-only, `cli.rs:906-916`). Folding `input_r1` in also keeps the four `output_dir = None` unit tests at `io.rs:854`, `:861`, `:875`, `:889` meaningful; leaving the namer alone strands them asserting a rule production can no longer reach.
7. **Pre-flight candidates, in the same commit as each writer:**
   - `main.rs:906-912` — paired trim reports, `Some(&pair_dir)`.
   - `main.rs:2026-2027` — paired uBAM trim reports.
   - `main.rs:893-903` — the carrier candidate.
   - **`clump_report_candidates`** (`main.rs:121-139`) — has **five** call sites, not two. `main.rs:616` (`&[r1, r2]`, one pair) and `:740` (`chunk[0]`) want the pair directory; `:653` and `:798` are handed **all** single-end inputs at once and need per-input parents when `output_dir` is `None`; `:698` is Shape B with one input. **Keep the parameter `Option<&Path>` with its existing "`None` → per-input parent" semantics and have only the paired callers pass `Some(&pair_dir)`.** A signature taking one required directory and applying it to every element would collapse every SE input's report into one directory, producing candidate paths the SE writers never write — the #383 family in mirror image, on an arm this plan declares out of scope.
   - Keep `OutputSource::Input` at every relocated candidate; see §Context.
8. **`PAIRED_REPORT_HINT`** (`main.rs:57-62`) — one hint serves three sites (`main.rs:608`, `:915`, `:2031`) **and is printed for primary collisions as well as report collisions**, so the reword must fit both. Its current text blames "`--output_dir` (or a shared input directory)"; after this change neither is required — two mates in different directories collide with no `-o` at all, so the stated cause is **absent in the failing case**. The reword must convey: paired outputs and reports all land in R1's parent unless `--output_dir` says otherwise, so inputs sharing a filename collide wherever they live; rename one, or pass `--output_dir`.
9. **Tests** — per Validation. Specifically: rewrite block 3 of the trim guard (`:1029-1041`) and the whole body of the clump twin (`:1206-1240`); flip `clump_paired_rejects_shared_mate_report_without_output_dir` (`:1247-1275`) to assert success; add the missing trim twin of that shape; add output-vs-input tests per paired report arm plus `--no_report_file` siblings. Use **`assert_dir_holds_only`** for "nothing was written" — **not** `assert_rejected_cleanly`, which asserts an *empty* directory and cannot express "unchanged beside the inputs"; reaching for it first produces a failure whose obvious fix damages eight other tests.
10. **Docs** (`docs/src/content/docs/`, lowercase) — `guide/outputs.md:102` is the one place the directory rule is stated and **it is already wrong today** ("writes outputs to `DIR/` instead of the current working directory" — trim outputs go to the *input's* parent, not the CWD); correct it and state the new rule there. `guide/outputs.md:6`'s "existing pipelines continue to work without changes" needs a qualifier: filenames still match v0.6.x, but the mate-side layout now requires `-o`. Add the rule to the Paired-end section (`:20-25`) and to `modes/passthrough.md`, since decision 2 moves the carrier. `guide/paired-end.md:16-19` and `modes/clump-only.md:45`/`:109`/`:122` were checked and carry no directory claim.
11. **CHANGELOG** — `#### Changes` does not yet exist under `### Unreleased` (only `#### Bug fixes`); the heading is precedented at `:205`, and `#### Behavioural notes (v2.x intentional widenings)` at `:998` is arguably the closer precedent for deliberately refusing a previously-working invocation. State: both layouts by name (mate-side refused, per-sample unaffected), `--output_dir` as the primary remediation — **not** `--no_report_file`, which throws the reports away — and that a shared-mate refusal disappears. The v0.6.11 fact bounds what may be claimed: this restores Perl's one-directory-per-run behaviour in spirit, but Perl's anchor was the CWD, not R1's parent, so "matches v0.6.x" is not available.
12. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`.

## Efficiency

Nil. One `PathBuf` per pair on the startup path, replacing one that was already being built per report.

## Integration

- **Reads:** `cli.input`, `cli.output_dir`, `cli.no_report_file`, `cli.passthrough`.
- **Writes:** relocates trimming reports (`.txt` + `.json`), clumping reports, and the carrier on paired arms only.
- **Report content knock-on:** the R2 text report's passthrough `Output:` line embeds the carrier's path (`report.rs:481`), so it moves with the carrier. JSON is unaffected.
- **Order:** all changes are at path-derivation time, before the pre-flight, which precedes every writer.
- **A free win from lockstep, worth recording:** moving a report candidate also moves its output-vs-input check, so a relocated report cannot silently land on a named input — provided step 7 is done with step 3.
- **Downstream:** two severities, not one. Scripts globbing `<mate>_trimming_report.txt` beside R2's input stop finding it; and runs on a mate-side layout **stop completing at all**. The changelog must state the second.
- **Perl parity:** no record content changes. Validation 7 is answered, not open — see below.

## Assumptions

- **A1:** The five primary namers' directory expressions are character-identical, so step 2 is a pure extraction. Verified by both reviewers independently.
- **A2:** The report namers derive filename and directory independently, so passing a directory relocates without renaming. Verified at `io.rs:535-580`.
- **A3:** Report paths are already pre-flight candidates on all paired report arms (#388, #391, #397), so new collisions are caught rather than silently overwritten. **Verified empirically by both reviewers** — this is the assumption the accepted-refusal decision rests on.
- **A4:** `--clock`/`--implicon` write no reports. Verified: `grep -c report src/specialty.rs` → 0.
- **A5:** Single-end paths never call the paired namers. Confirm via the single-end path assertions in validation 5. (r1 mis-cited `main.rs:909` here; that line is a paired candidate.)
- **A6:** `Path::parent()` on a bare filename yields `Some("")`, **not** `None` — so the `unwrap_or(Path::new("."))` arm never fires for that input and joined paths stay bare. r1 described this as the `unwrap_or` handling the case; it is the opposite, and the wrong reading invites a "simplification" that prefixes every path with `./`. The `unwrap_or` arm is reachable only for a root path.
- **A7 (new):** `collision_key` is ASCII-case-folded, so "fold-equal filenames" and not "identical filenames" is the collision condition.
- **A8 (new):** `--passthrough` has no single-end-shaped callers. Verified at `cli.rs:906-916`.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | The reported asymmetry is gone | `--paired A/r1.fq B/r2.fq`, no `-o` | both reports and both primaries in `A/`; `B/` holds only its input |
| 2 | Carrier joins the pair | `--paired A/r1 B/r2 --passthrough C/i1` | carrier output in `A/`; R2 text report's `Output:` line names the new path |
| 3 | **The accepted refusal fires, nothing written** | `--paired A/reads.fq B/reads.fq` | non-zero; 2a message naming both inputs; `assert_dir_holds_only` on **both** `A/` and `B/` |
| 4 | Mate-side layout refuses the whole run | `--paired R1/s1 R2/s1 R1/s2 R2/s2` | non-zero; **no pair processed** — neither `s1_val_1` nor `s2_val_1` exists |
| 5 | Per-sample layout unaffected | `sampleA/r_1 sampleA/r_2 sampleB/r_1 sampleB/r_2` | succeeds; each pair's outputs in its own directory |
| 6 | `--no_report_file` escape | validation 3 plus the flag | succeeds; primaries in `A/`; no reports |
| 7 | **Single-end untouched** | SE run, reports asserted at exact paths, plus a **bare-filename** input | byte-identical paths to `38deb1ee`; **stderr shows no `./` prefix** |
| 8 | **Output-vs-input refusal, per paired report arm** | `--paired A/reads.fq_trimming_report.txt B/reads.fq`, and the clump and uBAM equivalents | non-zero; nothing written; **R1's input intact and unmodified** (assert its content, not just existence) |
| 9 | …and its acceptance sibling | validation 8 plus `--no_report_file` | succeeds; the report-named input untouched |
| 10 | **Both inverted guards rewritten, not narrowed** | trim guard `:1029-1041` blocks 1–2 unchanged and passing, block 3 asserts refusal; clump twin `:1206-1240` body asserts refusal | as stated |
| 11 | **The lost refusals are recorded** | `clump_paired_rejects_shared_mate_report_without_output_dir` flipped to success; new trim twin asserts success | eight distinct paths; all four primaries and all four reports present |
| 12 | Writers and candidates agree | per arm, run it and compare every created file against that arm's candidate list. **Bound to runs without `--fastqc`** — FastQC paths are deliberately absent from candidate lists (`io.rs:1335-1340`); state that known-unplanned set rather than letting an exact-match assertion fail spuriously and get weakened | exact match on the three changing arms |
| 13 | Multi-pair anchors per pair | `--paired A/x1 A/x2 B/y1 B/y2` | pair 1 in `A/`, pair 2 in `B/` |
| 14 | Expected-fail control | validations 1, 2, 3, 4 and 8 against an unpatched build **rebuilt from HEAD** — `target/release/trim_galore` was found stale (stamped `3f1b0fe`) during review, so do not reuse whatever is in `target/` | 1, 2 show the split layout; 3, 4, 8 **succeed** there — the behaviour being removed |
| 15 | No collateral | full `cargo test` + fmt + clippy | green; count the delta against the **578** baseline rather than reading "ok" |

**Validation 7 of r1 is answered, not open.** Both reviewers checked `.github/workflows/ci.yml`: every output-producing invocation passes `-o`, and every fixture shares `test_files/`, so the Perl byte-identity matrix is unaffected twice over. Recorded here so nobody re-runs it.

## Questions or ambiguities

- **[Resolved — Felix]** Reports follow the primaries into R1's parent.
- **[Resolved — Felix]** `--passthrough`'s carrier joins them.
- **[Resolved — Felix]** New collisions are accepted as loud refusals, reconfirmed after seeing the mate-side layout.
- **[Resolved — this plan]** Keep `OutputSource::Input` for relocated candidates; `Pair` would route the message into the 2c branch whose advice is wrong for a two-input collision.
- **[Resolved — this plan]** The lost shared-mate refusal is correct behaviour and its test flips to asserting success. Flagged to Felix as a consequence he had not been shown.
- **[Answered by the Perl fact]** Should `--hardtrim5/3` share this helper? **No** — and now for a reason rather than a preference: v0.6.11 wrote *everything* to `$output_dir`/CWD, so the specialty modes' CWD anchor is the Perl-faithful one and the trim path's input-parent anchor is the v2 innovation. Folding them together would change the mode that currently matches Perl.
- **[Open, minor]** Whether to emit a one-time `NOTE` when a report's directory differs from its own mate's parent. Recommend no — it would fire on the common, now-correct case.
- **[Open, deferred]** Both reviewers raise a shape where writer/candidate drift becomes *unrepresentable* rather than tested-for: dedicated paired report namers taking `input_r1`, so no caller can reach the per-input rule by accident. After five recurrences and with #414 open precisely because tests are the weaker guarantee, this deserves a decision rather than a footnote — but it widens the diff beyond #398. **Recommend implementing #398 as specified, then evaluating the namer shape as part of #414.**

## Implementation notes (2026-08-10, branch `fix/398-pair-output-dir`, base `38deb1ee`)

Implemented as specified in r2. **The three predicted test changes were exactly the three that failed** — no unpredicted red test appeared, which was the point of enumerating them:

```
paired_report_candidates_do_not_over_reject                 (trim guard, block 3)
clump_report_candidates_do_not_over_reject                  (clump twin — A's Critical 1.2)
clump_paired_rejects_shared_mate_report_without_output_dir  (vanished refusal — A's Critical 1.3)
```

### What was done

| Step | Outcome |
|---|---|
| 1–2 | `pair_output_dir` added; all five primary namers refactored to it. The five directory expressions were character-identical, so a single `replace_all` was exact; their unit tests passed untouched, confirming the extraction is behaviour-preserving |
| 3–5 | Report writers on the trim-FASTQ, trim-uBAM and clump-FASTQ arms take `Some(&pair_dir)`. The clump-BAM arm routed through the same helper — a no-op as predicted, done for uniformity |
| 6 | `passthrough_output_name` gained `input_r1` and derives via `pair_output_dir`. Its six unit tests now put R1 in a **different directory** from the carrier, so a carrier-anchored implementation fails them |
| 7 | Candidates moved with their writers. **`clump_report_candidates` needed no signature change at all** — keeping `Option<&Path>` and passing `Some(&pair_dir)` from only the two paired sites was sufficient, so the two SE sites that hand it the whole input list keep per-input parents untouched |
| 8 | `PAIRED_REPORT_HINT` reworded to state the rule (all of a pair's outputs → `--output_dir` or R1's directory) and drop the false "shared input directory" precondition |
| 9 | Both guards rewritten to assert refusal; the vanished refusal flipped to assert success; 11 new tests |
| 10–11 | Docs and CHANGELOG per r2, including `outputs.md:102`'s pre-existing error and both layouts by name |

### Deviations

1. **No signature change to `clump_report_candidates`.** r2 said "keep the parameter `Option<&Path>`", and it turned out nothing else was needed — the paired callers pass `Some(&pair_dir)` and the SE callers are untouched. r1's version of this step would have changed the signature and broken `main.rs:653`/`:798`.
2. **The clump inverted guard was split into two tests, not just rewritten.** `clump_report_candidates_reject_shared_filename_without_output_dir` asserts the refusal, and a new `clump_shared_filename_without_output_dir_runs_without_reports` keeps the acceptance side alive — otherwise the arm would have had no test proving the primaries still work when only the reports collide.
3. **Both rewritten guards were renamed.** A test called `..._do_not_over_reject` that now asserts a rejection would read as a contradiction to the next person to open the file.

### Verification

| # | Validation | Result |
|---|---|---|
| 1 | Asymmetry gone | `paired_reports_follow_the_primaries_into_r1s_directory` — full listing of both dirs |
| 2 | Carrier joins the pair | `passthrough_output_follows_r1s_directory` — carrier in `A/`, asserted **absent** from `C/`, and the R2 report names the new path and not the old one |
| 3 | Accepted refusal, nothing written | trim guard block 3 + `clump_report_candidates_reject_shared_filename_without_output_dir`, both with `assert_dir_holds_only` on **both** directories |
| 4 | Mate-side layout refuses the whole run | `mate_side_directory_layout_refuses_the_whole_invocation` — asserts not even pair 1's primaries exist |
| 5 | Per-sample layout unaffected | `per_sample_directory_layout_is_unaffected` — each pair's outputs verified by read prefix, so a pair landing in the wrong directory fails |
| 6 | `--no_report_file` escape | `clump_shared_filename_without_output_dir_runs_without_reports` |
| 7 | SE untouched + bare filename | `bare_filename_r1_keeps_paths_unprefixed` asserts no `./` prefix in stderr **and** the exact directory contents; `test_pair_output_dir_precedence_and_edges` pins empty-parent, `--output_dir` precedence, and that `.` is reachable only for a root path |
| 8 | Output-vs-input per arm | `paired_report_that_would_land_on_r1s_input_is_refused` and `clump_paired_report_that_would_land_on_r1s_input_is_refused` — both assert the victim's **read count**, not just its existence |
| 9 | Acceptance sibling | `paired_report_input_alias_accepted_with_no_report_file` |
| 10 | Both guards rewritten, not narrowed | trim blocks 1–2 unchanged and passing; both bodies now assert refusal; each carries a comment naming its twin and forbidding the narrowing "fix" |
| 11 | Lost refusals recorded | `clump_paired_shared_mate_report_runs_when_pairs_anchor_apart` + the previously-untested `paired_shared_mate_report_runs_when_pairs_anchor_apart` |
| 12 | Writers/candidates agree | every refusal test asserts full directory listings; every acceptance test asserts the exact written set. No `--fastqc` run is asserted exactly, per r2 |
| 13 | Multi-pair anchors per pair | covered by validations 5 and 11 |
| 14 | Expected-fail control | Run against a build made **from `38deb1ee` in a fresh worktree** (not `target/`, per the stale-binary finding). All four behave as predicted: V1 shows the split layout (reports in both `A/` and `B/`); the **mate-side layout exits 0** pre-patch, writing all four pairs' outputs; the **report-on-input case exits 0** pre-patch and does *no* damage, because R2's report still goes to `B/` — confirming the hazard is created by this change and must be caught by the moved candidate, not that it pre-exists; and the **shared-mate case exits 1** pre-patch, confirming the refusal that this change removes |
| 15 | No collateral | `cargo test` **589 passed / 0 failed** (578 baseline + 11); `fmt --check` clean; `clippy --all-targets --release -D warnings` exit 0 with a genuine recompile |

### Iteration log

- **#1** First full run: 3 failures, exactly the three predicted. Fixed by rewriting the two guards and flipping the vanished refusal, per r2 §Context — not by touching candidate derivation.
- **#2** `cargo build --release` reported zero errors while the test targets were still broken by the `passthrough_output_name` signature change; `cargo test --no-run` surfaced 8. A release build does not compile test targets, so it cannot stand in for a compile check after a signature change.
- **#3** An `Edit` against the trim guard failed because r2's quoted comment wrapped differently from the file. Re-read the region and matched the real text rather than trusting the quote in the plan.

## Self-Review (r2)

- **What r1 got wrong, owned:** it described a one-directional behaviour change when the footprint is bidirectional; it built a prominent warning around one inverted guard and missed its near-identical twin one arm over; it never considered the output-vs-input direction, which is the shape of the data-loss bug fixed three commits earlier; and its step 7 would have broken two single-end arms it had itself declared out of scope. The pattern across all four is the same: r1 reasoned about the arm it was looking at and generalised without enumerating. The fix in each case came from enumerating call sites and test bodies rather than reasoning about them — which is the same lesson as #408's r1.
- **What the reviews added beyond corrections:** the v0.6.11 output-directory fact, which converts the `--hardtrim` open question from a preference into an answer and constrains the changelog; and the confirmation that `--fastqc` needs no change because it anchors on the output file.
- **Where both reviewers were wrong:** the docs path is lowercase `docs/`, not `Docs/`. Adopting their correction would have produced a path that fails on Linux.
- **Traps checked:** validation 3 asserts full directory listings on both sides, not one absent filename; validation 8 asserts R1's input *content*, since existence alone would pass if the file were truncated; validation 12 states its known-unplanned set so it cannot be quietly weakened when `--fastqc` trips it; validation 14 requires a control rebuilt from HEAD after a stale binary was found during review; validation 7 pins the bare-filename stderr spelling, which is the observable symptom of the A6 misreading.
- **Remaining risk:** the change refuses invocations that work today, including a real layout convention. The refusal is loud, pre-write, and names both inputs; there is no silent-loss path once step 7 lands with step 3. The three enumerated test flips are the live risk during implementation — each is a red or newly-green test whose obvious "fix" is the wrong one, which is why each is named with its line range and its correct resolution.
