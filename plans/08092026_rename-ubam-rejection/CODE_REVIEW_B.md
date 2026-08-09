# Code Review B — `fix/408-rename-ubam-rejection` @ `3f1b0fe`

**Reviewer:** B (independent; no shared state with Reviewer A)
**Date:** 2026-08-09
**Target:** `git diff dev...HEAD` — `src/main.rs` +18, `src/bam.rs` +3/−4, `tests/integration_ubam_out.rs` +172, `docs/…/guide/outputs.md` +11/−5, `CHANGELOG.md` +18
**Verdict:** **APPROVE with recommendations.** Nothing Critical. Two Medium findings, both about
accuracy of user-facing text and one unenforced invariant; no defect in the executable change.
**No fixes applied** (per the review brief — two reviewers share this tree).

---

## Verification performed

Everything below was run, not read off the plan.

| Check | Method | Result |
|---|---|---|
| Full suite | `cargo test --release`, all 14 targets, per-target counts summed by hand | **575 passed / 0 failed** — matches the stated baseline exactly |
| The 5 new tests + the C1 guard | `cargo test --release --test integration_ubam_out` | **33 passed / 0 failed**; all five new names appear in the `^test` list; `ubam_out_rename_with_preserve_tags_keeps_tags_intact` passes |
| C1 guard is not taking its `SKIP` branch | `git ls-files --stage test_files/ubam_test_with_tags.bam` | tracked, blob `4cd4670`, 587 B → `.exists()` true, the guard really ran |
| `cargo fmt --all -- --check` | direct | exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | run in an **isolated** `CARGO_TARGET_DIR` (fresh 205 MB target, 45 `Compiling` units, log line 188 `Compiling trim-galore v2.3.0`) | exit 0, zero warnings — **and I proved the invocation can fail**: injected `let unused_probe = 1;` into a `git archive`-pinned copy of the tree and clippy errored `unused variable: unused_probe … could not compile trim-galore (test "integration_ubam_out")`. So the green result is real, not a cache hit. (The crate is `trim-galore` with a hyphen — grepping the log for `trim_galore` finds nothing and looks like a false negative.) |
| Siting precedes every writer | line numbers in `src/main.rs` | guard **:316–333**; `ensure_output_dir` :406; hardtrim5 :450; hardtrim3 :481; clock :512; implicon :534; clump_only :566; `run_ubam_output` :828. Nothing between :229 and :333 writes to disk (`sanity_check_any` :241 and `detect_input_format` :250 are read-only). **Confirmed.** |
| Plan's claim that `--hardtrim5 20 --rename` genuinely annotates on the FASTQ path | ran the built binary | `@withspace 1:N:0:ACGTAC:clip5:ACGTACGT` (space header) and `@nospacehere:clip5:ACGTACGT` (space-free). **Independently confirmed** — validation 5 refuses something that really worked. |
| Mechanism | read `fastq.rs::append_to_id` :173–181, `bam.rs::parse_name_and_data` :707–726 | append-to-end (splice before first TAB), then split on the **first** of `\t` or `' '`. Confirmed. |
| `--paired` + single interleaved uBAM + `--rename` + ubam out | ran on `test_files/ubam_paired_test.bam` | **accepted**; `ubam_paired_test_val.bam` written, 10 `clip5` occurrences in the QNAMEs. The gate correctly admits the N==1 interleaved case. |
| All eight compatibility bullets really are refused before any write | read `cli.rs` §3.4a (`:584–630`) for the seven; `main.rs:316` for the eighth | **preamble is accurate** — `--dont_gzip`, `--clumpify`, `--passthrough`, `--clock`, `--implicon`, `--demux`, `--retain_unpaired` all `bail!` inside `Cli::validate()`, which runs at `main.rs:229` before the banner and before `sanity_check_any`. Deviation 2 was necessary and is correct. |

---

## Issues by area

### 1. Logic

**Clean on the gate predicate itself.** `input_formats.iter().any(|f| !matches!(f, InputFormat::UnalignedBam))`
is the right shape for the three cases I probed:

- mixed FASTQ + uBAM (SE, multi-input) → refused. Correct: the FASTQ member would lose the annotation.
- `--paired` with one interleaved uBAM → accepted (`input_formats.len() == 1`, all uBAM). Confirmed by running it.
- a hypothetical future `InputFormat` variant → refused, i.e. conservative. `InputFormat` has exactly
  three variants today (`format.rs:31–40`); deviation 4's reasoning holds, and the negative form is
  the right default here.

`cli.rename` has exactly one definition (`cli.rs:302–303`) and is never set implicitly, so the
predicate's first conjunct can't fire spuriously.

**M-1 (Medium) — the gate is coarser than the loss, and nothing says so.**
`src/main.rs:318–333`, `docs/…/guide/outputs.md:98`.

Two classes of invocation are refused although the annotation would survive intact:

1. **FASTQ input with no header description.** CONFIRMED by running the binary:
   `--rename --clip_R1 3 --output-format ubam` on a header `@nospacehere` is refused, yet the
   annotation demonstrably round-trips — I fed the annotated ID `@nospacehere:clip5:ACGTACGT`
   through `--output-format ubam` and the BAM QNAME came out as
   `nospacehere:clip5:ACGTACGT`, intact. This is not a rare shape: simulated data, many
   `samtools fastq`-derived FASTQs, and Trim Galore's own output from a description-free input all
   look like this.
2. **`--rename` with no `--clip_*` / `--hardtrim*` flag at all.** CONFIRMED:
   `trim_galore --rename --output-format ubam ns.fastq` is refused. But `append_to_id` is only
   reached under `if let Some(n) = clip_5 / clip_3` (`trimmer.rs:257–275`) and the hardtrim
   equivalents, so with no clip flag `--rename` is a pure no-op and nothing could be lost. A
   pipeline that passes `--rename` unconditionally alongside *optional* `--clip_R1` now hard-fails
   on uBAM output even in its no-clip configuration.

Neither is a bug in the sense of wrong output — the design decision to gate on format rather than
per-record is right (a per-record decision would mean failing mid-write, which is worse). The defect
is in the **text**: both the error message and the docs bullet describe a *conditional,
data-dependent* mechanism ("whenever a FASTQ header carries text after the first space…") next to
*unconditional* behaviour, and neither says the refusal is deliberately coarser than the loss. The
user whose headers are clean reads the message, concludes it does not apply to their data, and has
no route forward — the message offers only "use FASTQ output, or drop `--rename`".

This matters more than a normal wording nit because the whole point of the #408 fold-in was to fix
a docs sentence that misdescribed a flag; and because the plan (Behavior 4) explicitly set out to
state "unrepresentability, not a prediction of loss" — which is exactly the distinction the current
wording blurs.

*Recommendation:* one clause in both places, e.g. after "…cannot be represented in the output":
"Trim Galore refuses the combination for any FASTQ input rather than deciding per record, because
whether the annotation survives would otherwise depend on each individual header." That makes the
behaviour and the text agree without changing a line of logic.

**M-2 (Medium, pre-existing — this PR inherits it rather than introducing it) — assumption A5 is a
spec property, not an enforced one, and the hole loses aux tags.**
`src/bam.rs:891–898` (`bam_record_to_fastq`).

The format gate's justification is that "BAM read names cannot carry a description by SAM spec".
`bam_record_to_fastq` performs no whitespace validation on QNAME, and samtools will happily build a
BAM that violates the rule. CONFIRMED end-to-end:

```
# samtools view -b from a SAM whose QNAME is literally "name with space"
trim_galore --rename --clip_R1 3 --output-format ubam --preserve-tags CB bad.bam
→ exit 0; output QNAME is "name"; :clip5:ACG gone; CB:Z:AAACCC gone
   NOTE: … First dropped: "with space:clip5:ACG	CB:Z:AAACCC"
```

So on the arm the gate deliberately *admits*, a spec-violating input silently drops the preserved
aux tag — the exact C1 corruption class the tab-splice in `append_to_id` exists to prevent. The
#406 NOTE does fire and echoes the dropped text, so it is disclosed rather than silent, and the
behaviour is identical with and without this commit — **this change regresses nothing.** But the
plan's "safe by construction, not by luck" is stronger than what the code guarantees, and a
one-line `bail!` on `name.contains(|c: char| c.is_ascii_whitespace())` in `bam_record_to_fastq`
would turn a silent aux-tag loss into an error for every mode, not just this one.

*Recommendation:* not a blocker for this PR. File a follow-up issue; if it lands, A5 becomes true
of the code and not only of the spec.

**L-1 (Low) — on a mixed `--paired` pair, the #408 message pre-empts the more fundamental one, and
its closing sentence points at a dead end.** `src/main.rs:316` vs `:377–392`.

CONFIRMED: `--paired --rename --clip_R1 3 --output-format ubam ns.fastq test_files/ubam_test.bam`
emits the #408 refusal, not `reject_bam_format_mismatch_in_pair`'s message. A mixed pair is
unsupported with or without `--rename`, so the user gets a two-step fix journey. Worse, the message
ends "uBAM input is unaffected", which invites converting the FASTQ mate to uBAM — and a two-BAM
pair is itself rejected (`ubam_out_two_bam_pair_rejected`), because paired uBAM wants one
interleaved file. Narrow, and arguably acceptable ordering, but worth knowing.

### 2. Efficiency

**Clean.** One boolean conjunction on the startup path, reusing `input_formats` that `:247–251`
already computed for `any_bam`. `.any()` short-circuits. Nothing measurable.

### 3. Errors

**Clean.** The only executable addition is an `anyhow::bail!`, sited ahead of `ensure_output_dir`
and every dispatch branch (verified by line number above, and asserted by
`rename_into_ubam_refused_for_fastq_input`'s `!dir.join("sp_trimmed.bam").exists()`). No new
unwraps, no new I/O, no ordering hazard. No security surface.

### 4. Structure and tests

**Test quality: good, none vacuous.** I checked each for the failure modes the brief names:

- All five assert on the *specific* message text or the *exact* ID, so none can pass because the
  run failed early for an unrelated reason. `clump_only_rename_keeps_its_own_message` is the
  strongest of the five: it asserts `byte-identically` **present** and `--rename cannot be honoured`
  **absent**, which pins the ordering between `Cli::validate()` and the new guard — something a
  library unit test genuinely could not do. Deviation 1 (dropping the planned `cli.rs` unit tests as
  unreachable) is correct: the guard needs `input_formats`, which only the `[[bin]]` has.
- No silent skips. None of the five has a `SKIP` branch. `rename_into_ubam_accepted_for_ubam_input`
  reads the tracked fixture directly and would fail loudly if it vanished — arguably better than the
  C1 guard's `SKIP`, which it sits next to.
- **Length filter cleared.** `write_one_record` (`:878–881`) writes `"I".repeat(len)` → Phred 40, so
  default `-q 20` trims nothing; 28 bp − `--clip_R1 3` = 25 bp ≥ default `--length 20`. The
  uBAM-accepted test uses the tracked fixture with `--clip_R1 5` and produced records. All confirmed
  by the tests passing with real assertions on output content.
- **The untouched C1 guard is genuinely still live.** `git diff -U0` on the test file is a single
  append hunk, the guard's name appears nowhere in the diff, the fixture is git-tracked (587 B, so
  the `SKIP` return is not what passed), and it passes by name in the 33-test run. Confirmed three
  ways. Its assertion (UB tag value free of `:clip5:`) still means what it claims.

**L-2 (Low) — `run_in` duplicates the existing `run_capturing_stderr`.**
`tests/integration_ubam_out.rs:983–995` vs `:883–890`. Same body, same `current_dir`, differing only
in whether the exit status is returned. Cleaner to give `run_capturing_stderr` a `(bool, String)`
return and update its three existing callers, or express one in terms of the other. Cosmetic.

**Comment policy: compliant.** The new `main.rs:316–317` comment is two lines, states the fact
(why it is sited there), and never references the state being fixed. It is the third comment in
`main()` making the same "validate() cannot see the format" point (`:271–275`, `:299–304`), so a
purist could shorten it to `// #408 — format-gated, so it cannot live in Cli::validate(). See §3.4b.`
— but as written it is within policy.

**`bam.rs` doc-comment rewording: sound.** I verified the removed example really is unreachable.
`emit_description_dropped_once` fires only on FASTQ→uBAM with a space in the header; with `--rename`
now refused for FASTQ input, and uBAM input carrying no description, the annotation can no longer be
part of the dropped text. The other two `append_to_id` families that could have reached it
(`specialty.rs:312/313` `--clock`, `:413/414` `--implicon` — note these are *not* `--rename`-gated)
are both rejected at CLI-validate with `--output-format ubam`, so they cannot reach it either. The
two surviving examples are real. Dropping the `#408` cross-reference rather than rewording it
matches the repo's "current state only" policy. (Aside: the plan's "six call sites in `specialty.rs`"
counts the clock/implicon pairs, which are unconditional rather than `--rename`-driven — a plan
inaccuracy only, no effect on the code.)

### 5. Docs

**L-3 (Low, but it is the same defect class the PR is correcting) — the new `--basename` paragraph
overclaims in two ways.** `docs/…/guide/outputs.md:106`.

"`--basename BASE` replaces the input filename stem in the output names" is unqualified. CONFIRMED
by running the binary:

- `--hardtrim5 20 --basename MYSAMPLE ns.fastq` → output is `ns.20bp_5prime.fq`. **`--basename` is
  silently ignored by the specialty modes** (no `basename` reference exists anywhere in
  `src/specialty.rs`, and `cli.rs` does not reject the combination).
- On the normal SE path the trimmed file is `MYSAMPLE_trimmed.fq` as documented, but the reports are
  `ns.fastq_trimming_report.{txt,json}` — **reports keep the input-derived name.** (That is correct,
  Perl-matching behaviour; it is the sentence that is loose.)

Also, the section sits directly under the uBAM output documentation yet omits the uBAM forms, which
`io.rs:282–291` implements: `BASE_trimmed.bam` and `BASE_val.bam`. `docs/…/modes/clump-only.md:79`
already documents its own `BASE_clumped.bam`, so the omission is visible.

The rest of the paragraph is accurate — I checked "only valid for one file (single-end) or one pair
(paired-end)" against `cli.rs:640` (`!paired && input.len() > 1`) and `:667`
(`paired && input.len() > 2`). ✓

**L-4 (Low) — the new "Annotating read IDs" section is accurate but incomplete.**
`docs/…/guide/outputs.md:110`. Every clause checks out against `cli.rs:301–303` and the call sites:
boolean ✓, does not affect filenames ✓, `:clip5:`/`:clip3:` ✓, the three clip-flag families ✓, "each
half only when that side was clipped" ✓ (`trimmer.rs:257–275`), anchor `#feature-compatibility`
resolves ✓. What is missing: `--rename` is also refused by `--clump_only` **regardless of output
format** (`cli.rs:851`, "preserves record contents byte-identically"). A section that names one
incompatibility and not the other reads as exhaustive.

**L-5 (Low) — the new compatibility bullet is the only multi-sentence entry in a seven-bullet
one-clause list**, and its final sentence ("checked after input-format detection rather than at
CLI-validate time, because the rule depends on the format") is implementation detail that the
generalised preamble already makes unnecessary. Style only.

### 6. CHANGELOG

Every mechanical claim checks out — I verified the append-to-end behaviour, the whitespace
constraint, the tail discard, the hardtrim coverage via `specialty.rs`, the uBAM-input aux-tag
round-trip, and the careful non-overclaiming framing of the #406 notice. The entry is
appropriately placed (newest-first, above #409) and correctly avoids "silently dropped".

**L-6 (Low) — one factual overclaim.** "`--basename` — a flag the guide did not document at all"
is not true. `grep -rn -- "--basename" docs/src/` finds it documented in
`docs/…/modes/clump-only.md:79` (including its uBAM naming) and `docs/…/modes/passthrough.md:28`.
What is true is that it had no dedicated section in `guide/outputs.md`. Also, the entry paraphrases
the old sentence as "replacing the **output** filename stem" where it said "replaces the **input**
filename stem in the output names" — trivial, but the entry is making a point about wrong
descriptions of flags.

**Not a finding:** `docs/…/reference/changelog.md` (the synced copy) was not updated. It is equally
stale for #397, #406 and #409 (`grep -c "issues/N"` → 0 for all four), so this PR follows existing
practice.

---

## Recommendations, by priority

**Critical:** none.

**High:** none.

**Medium**
1. **M-1** — add one clause to the error message (`src/main.rs:324–332`) and to the docs bullet
   (`outputs.md:98`) saying the refusal is deliberately format-level rather than per-record.
   Without it the text describes a conditional mechanism beside unconditional behaviour, and a user
   with description-free headers is refused a path I confirmed works, with no explanation and no
   alternative. Text-only; no logic change.
2. **M-2** — file a follow-up to validate whitespace in `bam_record_to_fastq` (`src/bam.rs:891–898`).
   Confirmed data loss (annotation **and** preserved aux tag) on a samtools-producible
   spec-violating uBAM, on the arm this gate admits. Pre-existing and unchanged by this commit —
   should not block the merge, but it is the one hole in the assumption the gate rests on.

**Low**
3. **L-3** — qualify the `--basename` paragraph (`outputs.md:106`): specialty modes ignore it
   (confirmed), reports keep input-derived names, and add the `BASE_trimmed.bam` / `BASE_val.bam`
   uBAM forms.
4. **L-6** — soften "a flag the guide did not document at all" in `CHANGELOG.md` to "had no
   dedicated section in the output guide"; `clump-only.md:79` and `passthrough.md:28` document it.
5. **L-4** — mention the `--clump_only` incompatibility in the new "Annotating read IDs" section.
6. **L-1** — consider whether the mixed-`--paired` case should reach
   `reject_bam_format_mismatch_in_pair` first; at minimum be aware the closing sentence of the
   message points at a combination that is itself rejected.
7. **L-2** — fold `run_in` and `run_capturing_stderr` into one helper.
8. **L-5** — trim the new compatibility bullet to the list's one-clause house style.

None of these blocks the merge. The executable change is correct, minimal, correctly sited, and
well covered.

---

## Notes for the caller

- I applied **no** fixes and modified **no** tracked file. `git status` shows only the pre-existing
  untracked entries plus this report.
- Scratch artefacts live under `$TMPDIR/rev408b*` and `$TMPDIR/tree408b` (a `git archive`-pinned
  copy used for the clippy negative control, with a deliberate warning appended to its *copy* of
  `tests/integration_ubam_out.rs` — the repo's file is untouched).
- Harness traps I hit and worked around, in case they bite the other reviewer: the crate is
  `trim-galore` (hyphen), so grepping a cargo log for `trim_galore` looks like a false negative;
  and `cargo clippy` on the shared target dir returns in well under a second as a cache hit, which
  is why I used an isolated `CARGO_TARGET_DIR` plus an injected-warning control.
