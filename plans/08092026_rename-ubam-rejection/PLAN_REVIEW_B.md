# Plan Review B — reject `--rename` with `--output-format ubam` (#408)

**Reviewer:** B (independent; Reviewer A works in parallel with no shared state)
**Plan:** `plans/08092026_rename-ubam-rejection/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Verified against:** `cbb5ecd` (plan states base `d0bb76b`; HEAD has since advanced by #397)
**Method:** `git archive HEAD` into `$TMPDIR/revB408/src`, private `CARGO_TARGET_DIR`, release build,
plus the plan's own proposed patch applied to that copy and the full test suite run against it.

**Verdict: REVISE before implementing.** The deciding evidence holds up — but the plan
as written will fail `cargo test`, and its proposed error message says something that is
untrue for a whole class of users.

Scope note: the maintainer has decided to reject rather than repair, and to fold in the
`outputs.md` error. This review takes both as given.

---

## Summary of verification

| # | Claim | Result |
|---|---|---|
| 1 | Reproduction (FASTQ keeps annotation, uBAM drops it; data-dependent) | **TRUE**, reproduced byte-exactly |
| 2 | Perl v0.6.11 appends to the end of the whole ID line | **TRUE**, and more robustly than the plan claims |
| 3 | `append_to_id` is format-shared via `trim_read` | **TRUE** |
| 4 | §3.4a is the right home; shape matches; in-block ordering immaterial | **TRUE**, but see I2 — §3.4a pre-empts a *different* block |
| 5 | A2: nothing asserts the combination succeeds | **FALSE — one test does. Deferral was unsafe.** |
| 6 | `outputs.md:105` is wrong; compatibility list missing `--rename` | **TRUE**; nothing else in that list is stale |
| 7 | Over-rejection limited to whitespace-free headers | **Understated** — the whole uBAM-input class is lossless today |

---

## Logic review

### The reproduction is exact (claim 1 — verified)

Built the release binary from `cbb5ecd` and ran the plan's scenario on a two-record
fixture, one space-bearing header and one without:

```
FASTQ output  : @withspace 1:N:0:ACGTAC:clip5:ACG
                @nospace:clip5:ACG
uBAM QNAMEs   : withspace                  ← :clip5:ACG gone
                nospace:clip5:ACG          ← preserved
```

The plan's table and its data-dependence contrast are both correct, and the mechanism is
as described: `append_to_id` (`src/fastq.rs:173-181`) splices before the first `\t` only,
and `parse_name_and_data` (`src/bam.rs:707-726`) splits on the first `['\t', ' ']` and
discards the remainder on the space branch.

### The Perl-parity argument is sound (claim 2 — verified, and stronger than stated)

This is the plan's deciding evidence, so I checked whether `$l1` could already have been
stripped of its description before the append. It has not been. Every relevant variable is
a raw filehandle read, and the only mutation before the append is a line-ending strip:

| Perl site | Variable | Assigned at | Mutation before append |
|---|---|---|---|
| `:1098`, `:1115` | `$l1` | `:993` `my $l1 = <TRIM>;` | `s/\r|\n//g` only |
| `:1282`, `:1298` | `$l1` | `:1219` `my $l1 = <TRIM>;` | `s/\r|\n//g` only |
| `:1399`, `:1414` | `$l1` | `:1372` `my $l1 = <TRIM>;` | `s/\r|\n//g` only |
| `:2209`, `:2226` | `$id_1`/`$id_2` | `:2161`/`:2166` `= <IN1>`/`<IN2>` | `s/\r|\n//g` only |
| `:1910`, `:2000` | `$identifier` | `:1893`/`:1977` `= <$in>` | `s/\r|\n//g` only |

I grepped every `$l1` occurrence between the assignment at `:993` and the append at
`:1098` — there are four, and they are the assignment, a `print`, the line-ending strip,
and the append itself. No `split`, no whitespace handling. So Perl really does append
after any space-separated description, at all six sites, on both the trim path and the
hardtrim path. **The maintainer did not choose on a false premise.**

The plan's noted irony also checks out: `trim_galore:3456` warns the format will be
`readname:clip5:GATC` — the read *name* — which the implementation does not do.

### Claim 3 confirmed

`trim_read` (`src/trimmer.rs:89`) is reached from both format families, so no fix inside it
could be format-scoped:

- FASTQ: `run_single_end` (`:319`), `run_paired_end` (`:479`/`:480`), plus `parallel.rs:463/464/1027`
- uBAM: `run_single_end_to_bam` (`:645`), `run_paired_end_to_bam` (`:730`/`:731`)

### CRITICAL — A2 is false; the plan's own patch breaks a committed test

The plan defers A2 ("no test or CI grep asserts that `--rename` + uBAM **succeeds**") to
implementation. I resolved it instead: **there is such a test, and it is not a grep-level
coincidence — it is a deliberate regression guard.**

`tests/integration_ubam_out.rs:393-440`, `ubam_out_rename_with_preserve_tags_keeps_tags_intact`,
runs `--clip_R1 5 --rename --output-format ubam --preserve-tags CB,UB` and asserts
`status.success()` at `:420`. Its fixture `test_files/ubam_test_with_tags.bam` **is
committed and git-tracked**, so the `exists()` guard at `:400` does not skip it.

I applied the plan's exact proposed snippet to my copy of `src/cli.rs` §3.4a and ran the
full suite:

```
src/lib.rs unit tests ......................... 411 passed
integration_adapter2 / clump_only / clump_only_ubam / gzip_non_gz_extension /
no_args_help / non_restartable_input / output_collision / paired_format_guard /
passthrough / ubam ............................ all green
integration_ubam_out .......................... 27 passed; 1 FAILED
  ---- ubam_out_rename_with_preserve_tags_keeps_tags_intact ----
  panicked at tests/integration_ubam_out.rs:420:5: trim_galore exited non-zero
```

Blast radius is precisely one test, which is good news — but the plan's implementation
outline has no step for it, and step 6 (`cargo test`) would simply fail. The deferral was
**not** safe: this is exactly the #400 lesson the plan cites, and it landed on the wrong
side of it.

**Do not just delete the test.** Its C1 mechanism is still live on a legal path —
uBAM input → FASTQ output → `--rename` still produces a tab-tail ID where the splice
matters. Verified:

```
$ trim_galore --clip_R1 5 --rename --preserve-tags CB,UB --dont_gzip ubam_test_with_tags.bam
@SRR24827378.1:clip5:AATTA<TAB>CB:Z:ATCGATCG-1<TAB>UB:Z:GCTAGCTA
```

Without the C1 fix the suffix would land after `UB:Z:GCTAGCTA` and corrupt it. So the
correct remediation is to **port the test to FASTQ output** (swap `--output-format ubam`
for `--dont_gzip`, assert on the FASTQ header instead of BAM aux), keeping its regression
value, and add the new rejection test alongside.

### CRITICAL — "silent" is false, and the proposed message is false

Two separate problems, both in the wording the plan intends to ship.

**(a) The loss is already disclosed at runtime**, as of the plan's own stated base commit.
#406/#412 added `emit_description_dropped_once` (`src/bam.rs:1031-1045`), which fires on
exactly this path and echoes the dropped text:

```
NOTE: BAM read names cannot contain whitespace, so FASTQ header text after the first
space is not carried into uBAM output. First dropped: "1:N:0:ACGTAC:clip5:ACG"
```

That notice's own doc comment (`src/bam.rs:1027-1030`) even forward-references this issue:
"or — see #408 — `--rename`'s own `:clip5:` annotation." So the plan describes as "silent,
data-dependent loss" something the codebase already announces, in a message written
with #408 in mind. The Goal paragraph, the "silent" framing throughout, and especially
implementation step 5 ("CHANGELOG … it closes silent data loss") all need correcting —
otherwise the changelog ships a claim that is refuted by the program's own stderr.

The refusal is still an improvement (a one-time NOTE about "header text" does not tell a
`--rename` user their annotation specifically is gone, and it is easy to miss in a
pipeline log). But the justification is "an easily-missed notice becomes an actionable
refusal", not "silence becomes a refusal".

**(b) The message would be untrue for uBAM-input users.** The proposed text asserts the
annotation "would be silently dropped". For uBAM input it would not be dropped at all —
BAM QNAMEs cannot contain whitespace by spec, so a uBAM-derived ID never reaches the space
branch, and the tab branch handles the tag tail correctly. Verified both ways:

```
uBAM in → uBAM out, --rename --preserve-tags CB,UB
  QNAME SRR24827378.1:clip5:AATTA   CB:Z:ATCGATCG-1  UB:Z:GCTAGCTA   ← lossless
uBAM in → uBAM out, --rename (no --preserve-tags)
  QNAME SRR24827378.1:clip5:AATTA                                     ← lossless
  (and no description-dropped NOTE is emitted)
```

Shipping a rejection whose stated reason is counterfactual for the affected user is the
precise failure mode #406 was filed about — and the #406 entry sits immediately above the
new one in `CHANGELOG.md`. The message must be phrased in terms of what cannot be
*guaranteed* (a FASTQ description makes the annotation unrepresentable in QNAME), not as a
claim about what will happen.

### IMPORTANT — the over-rejection is materially larger than the plan concedes

The plan's Self-Review treats whitespace-free headers as the only over-reach and calls it
an accepted edge case. In fact the refusal removes **the entire uBAM-in → uBAM-out
`--rename` path**, which works losslessly today, has a committed fixture, and has a
dedicated regression test — and which is the flagship round-trip the `--output-format ubam`
feature exists to serve. That is not an edge case; it is the feature's core audience.

To be clear, **the decision can still stand**, and there is strong in-tree precedent for
it that the plan argues around without citing. `src/main.rs:265-269`, on the `--phred64`
rejection, records the house position almost verbatim:

> Rejected uniformly rather than per-mode. `--hardtrim5/3` and `--clump_only` currently
> accept the flag as an inert no-op, so this does remove working invocations — but a
> mode-dependent rule … is worse to document, worse to test, and one refactor away from
> being wrong.

The plan should cite that rather than re-deriving the argument, and should state the cost
accurately.

But the plan should also surface the alternative it never considers, because it has an
established home. A format-gated rule is not novel here: `main.rs:271-275` explains that
§3.4b lives in `main.rs` precisely because "validate() runs before input format detection
and cannot see `any_bam`". So `--rename` + uBAM output + any FASTQ input could be rejected
next to §3.4b, preserving the lossless uBAM-in case. That trades a wider blast radius for
a rule keyed on input format rather than mode — arguably a different and more defensible
axis than the "mode-dependent" shape `main.rs:265-269` warns against. **Felix chose
"reject" over "repair"; he was not offered "reject narrowly", and given the uBAM-in case
is lossless he should be.**

### IMPORTANT — §3.4a pre-empts the `--clump_only` block

Behavior item 4 says "Ordering within §3.4a is immaterial — every arm bails, so no arm
masks another." True as far as it goes, but it addresses the wrong ordering question.
§3.4a is the **first** thing in `validate()` (`src/cli.rs:587`); the `--clump_only` rename
rejection is at `:851`. So the new arm masks that one.

`--clump_only --output-format ubam` is a supported combination (`outputs.md:96`), so
`--clump_only --rename --output-format ubam` is reachable. Today:

```
Error: --clump_only preserves record contents byte-identically; --rename would mutate read IDs
```

After the change the user instead gets the uBAM message about `:clip5:`/`:clip3:` being
appended after a header description — a non-sequitur under `--clump_only`, which does no
clipping at all, so no annotation would ever exist. An accurate message regresses to an
inaccurate one.

The plan's Open question asks whether to *mention* the `--clump_only` rejection in the new
message. The real question is whether to *avoid pre-empting* it — e.g. `if self.rename && !self.clump_only`,
or siting the new arm after the clump-only block. Either is a one-line change; the plan
should decide it rather than leaving it to discovery.

### IMPORTANT — the hardtrim uBAM paths share the defect and go unmentioned

The plan's Context attributes the loss solely to `trim_read`. There are two more
independent call sites: `src/specialty.rs:166` (`hardtrim5_to_bam`) and `:219`
(`hardtrim3_to_bam`), each gated on `rename` and each writing through `BamWriter`.
`--hardtrim5`/`--hardtrim3` are **not** rejected in §3.4a, and they have real uBAM
drivers (`main.rs:450`, `:481`). Verified the same loss:

```
$ trim_galore --hardtrim5 10 --rename --output-format ubam in.fastq
NOTE: … First dropped: "1:N:0:ACGTAC:clip5:GTACGTACGTACGTACGTACGT"
QNAMEs: withspace                              ← annotation gone
        nospace:clip5:GTACGTACGTACGTACGTACGT
```

This cuts in the plan's favour — the blanket §3.4a arm closes the hardtrim variant too,
so the fix is broader than the plan realises. But three things follow: the Context section
is incomplete as written; the Validation table never exercises hardtrim, which is
plausibly where `--rename` sees most use; and the CHANGELOG bullet should say the refusal
covers hardtrim as well, or hardtrim users will not connect it to their runs.

### IMPORTANT — `bam.rs:1027-1030` goes stale on merge

That doc comment offers "`--rename`'s own `:clip5:` annotation" as an example of text the
notice may report as dropped, tagged "see #408". Once `--rename` + uBAM output is refused,
`:clip5:` can no longer reach `parse_name_and_data` by any route, so the example becomes
counterfactual — a comment citing the issue whose fix invalidated it. The plan has no step
for it. One-line edit, but it is the kind of thing #406 was about.

---

## Assumptions

- **A1 (rename absent from §3.4a) — TRUE.** §3.4a spans `src/cli.rs:584-632` and rejects
  `dont_gzip`, `clumpify`, `passthrough`, `clock`, `implicon`, `demux`, `retain_unpaired`.
  The only `self.rename` rejection in the file is the `--clump_only` one at `:851-855`.
  (Line numbers drifted slightly from the plan's "~590-631" — cosmetic.)
- **A2 — FALSE.** See above. This is the plan's one substantive error.
- **A3 (not in the byte-identity matrix) — TRUE.** `.github/workflows/ci.yml`'s only
  `rename` hits are the branch-rename comment at `:8` and a file rename in a diff step at
  `:1002`. Neither is a `--rename` invocation.
- **A4 (FASTQ path untouched, Perl parity untouched by construction) — TRUE.** The change
  cannot reach `append_to_id`.
- **Unstated assumption, now falsified:** that the loss is silent, and that the uBAM path
  loses the annotation generally rather than only for space-bearing (i.e. FASTQ-sourced)
  headers. Both are load-bearing for the Goal paragraph and the message text.
- **Unstated assumption:** that no *other* validation block depends on running before
  §3.4a. Only the `--clump_only` one does, and it now loses. See I2.

---

## Efficiency

Nil, as the plan says: one boolean test on the argument-parsing path. Nothing further to add.

---

## Validation sufficiency

The table is well-shaped — row 5's expected-fail control is genuinely informative here —
but it has gaps, and one row cannot fail.

- **Missing: the test that breaks.** No row anticipates
  `ubam_out_rename_with_preserve_tags_keeps_tags_intact`. Row 4 ("A2's grep before
  editing") would have caught it only if the grep were run; it was deferred instead.
- **Missing: the uBAM-input case.** Nothing checks what the refusal costs a
  uBAM-in → uBAM-out user, which is the combination that works today. A row asserting the
  refusal *and* acknowledging the loss of a working path would have surfaced I1 during
  planning.
- **Missing: hardtrim.** `--hardtrim5 N --rename --output-format ubam` should appear, both
  as a rejection case and because it exercises a different `append_to_id` call site.
- **Missing: `--clump_only` message stability.** Add a row asserting which message
  `--clump_only --rename --output-format ubam` produces, whichever way I2 is decided.
- **Row 3 cannot fail.** "Perl parity untouched — the byte-identity harness plus the
  validation matrix" is a no-op given A3 (correctly) establishes that `--rename` is not in
  the matrix. It is reassurance, not validation. Row 2 (FASTQ output still ends
  `:clip5:ACG`) is the row that actually protects parity, and it is good — I confirmed it
  passes on the patched build's FASTQ path.
- **Suggest matching sibling conventions** for the new unit test:
  `output_format_ubam_plus_rename_rejected`, beside `output_format_ubam_plus_demux_rejected`
  (`src/cli.rs:2146`). The plan's "assert `--rename` alone still validates" already has a
  model at `output_format_ubam_alone_validates` (`:2052`).

---

## Docs claims (claim 6 — verified)

- **`outputs.md:105` is wrong on every clause**, as the plan says. `--rename` is
  `pub rename: bool` (`src/cli.rs:303`); it takes no `PREFIX`, touches no filename, and
  the described behaviour is `--basename`'s.
- **The compatibility list (`:88` heading, `:90-96` items) is missing `--rename`** — and
  **nothing else in it is stale**. I checked all seven bullets against §3.4a: they match
  one-for-one, and the `--clump_only --output-format ubam` parenthetical at `:96` is
  accurate. The `--preserve-tags` hard error is covered in prose at `:87`. So the plan's
  "only one gap" reading is correct.
- **The plan's step-3 conditional resolves to "yes, a clause is needed."** `--basename`
  appears nowhere in the guide — only in `modes/passthrough.md:28`, `modes/clump-only.md:79`,
  and the changelog. The outputs page is exactly where it belongs, and it currently
  documents `--basename`'s behaviour under the wrong flag's name while never naming
  `--basename` itself.
- **One thing the plan misses: the heading is wrong too.** `## Renaming outputs` (`:103`)
  is `--basename`'s subject, not `--rename`'s. Replacing only the sentence leaves a
  correct sentence about read-ID annotation under a heading about output naming. Either
  keep the heading for `--basename` and give `--rename` its own, or retitle.

---

## Alternatives

1. **Format-gated rejection in `main.rs` beside §3.4b** (reject only when a FASTQ input is
   present). Preserves the lossless uBAM-in case; established precedent and an established
   home; costs a rule that reads off input format. **This is the alternative worth putting
   to Felix**, because the plan's cost estimate for blanket rejection was too low.
2. **Blanket rejection as planned**, with the message fixed and `main.rs:265-269` cited as
   precedent. Defensible, simplest to document and test, consistent with house style. My
   recommendation *if* Felix still prefers uniformity once he knows the uBAM-in case works.
3. **Warn instead of refuse**, upgrading the existing NOTE to name `--rename` specifically.
   Rejected — it keeps a data-dependent outcome, which the plan correctly argues against.
4. **Repair `append_to_id` with a format flag threaded through `trim_read`.** Correctly
   ruled out: it either breaks Perl parity or widens `trim_read`'s signature for one flag.

Options 1 and 2 are both reasonable; the plan should not proceed as though only 2 exists.

---

## Action items

### Critical

1. **Fix `tests/integration_ubam_out.rs:393-440`.** A2 is false; the plan's exact patch
   fails this test. Port it to FASTQ output (`--dont_gzip` in place of
   `--output-format ubam`, assert on the FASTQ header) so the C1 tab-splice guard survives
   — verified that path still exercises it. Add an explicit implementation step; do not
   leave it to `cargo test` to discover.
2. **Rewrite the error message.** "would be silently dropped" is false for uBAM input
   (lossless, verified) and overstated for FASTQ input (a NOTE already fires). Phrase it
   as unrepresentability, not prediction. Shipping the current wording regresses the exact
   principle #406 established one changelog entry above.
3. **Drop "silent" from the Goal and the CHANGELOG bullet.** `emit_description_dropped_once`
   (`src/bam.rs:1031-1045`) exists as of the plan's own base commit and echoes the lost
   annotation. The honest framing is "an easily-missed notice becomes an actionable refusal".

### Important

4. **Restate the cost of rejection.** The whole uBAM-in → uBAM-out `--rename` path is
   lossless today and would be refused — not just whitespace-free headers. Correct the
   Self-Review edge-case bullet, cite `main.rs:265-269` as the precedent for uniform
   rejection, and put alternative 1 (format-gated, beside §3.4b) to Felix explicitly.
5. **Decide the `--clump_only` pre-emption.** §3.4a (`:587`) runs before the clump-only
   block (`:851`), so the new arm replaces an accurate message with an inapplicable one
   for `--clump_only --rename --output-format ubam`. Behavior item 4 addresses in-block
   ordering only and misses this.
6. **Add the hardtrim paths to Context, Validation and the CHANGELOG.**
   `src/specialty.rs:166`/`:219` are independent `append_to_id` call sites with the same
   defect; the refusal covers them, which is a point in its favour, but nothing in the plan
   says so or tests it.
7. **Update `src/bam.rs:1027-1030`.** Its "see #408 — `--rename`'s own `:clip5:`
   annotation" example becomes unreachable once the combination is refused.
8. **Fix the `outputs.md:103` heading**, not just the `:105` sentence. And note that the
   step-3 conditional is already answered: `--basename` is undocumented in the guide, so
   the clause is required.

### Optional

9. Refresh the base commit and line numbers: base is now `cbb5ecd` (#397 landed after
   `d0bb76b`); §3.4a is `src/cli.rs:584-632`.
10. Replace Validation row 3 with something that can fail, or mark it explicitly as
    reassurance. Row 2 is the one doing the parity work.
11. Name the new test `output_format_ubam_plus_rename_rejected` to match the siblings at
    `src/cli.rs:2059-2160`.
