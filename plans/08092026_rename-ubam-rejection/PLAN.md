# Plan: refuse `--rename` into uBAM when the annotation cannot be represented (#408)

**Issue:** [#408](https://github.com/FelixKrueger/TrimGalore/issues/408) — filed from #406's plan review.
**Decisions (Felix):** reject rather than repair; fold in the `outputs.md` docs error; **r2: format-gated — refuse only when at least one input is FASTQ.**
**Base:** `dev` @ `d0bb76b` (r1 cited this; both reviewers audited at `cbb5ecd`, which additionally carries #397 — no finding depends on the difference).

## Revision history

- **r2 (2026-08-09):** dual plan review (PLAN_REVIEW_A.md, PLAN_REVIEW_B.md — both REVISE; B applied r1's exact patch to a copy and ran the suite). The Perl-parity evidence the decision rests on was **independently confirmed by both**. Four things changed: **scope is now format-gated** (Felix, on both reviewers' finding that uBAM→uBAM is lossless *by spec*), which **relocates the check from `cli.rs` §3.4a to `main.rs` §3.4b** and thereby resolves the `--clump_only` pre-emption for free; **A2 was false** and is corrected; the error message is reframed as unrepresentability rather than prediction; and "silent" is dropped throughout, because #406 — merged before r1 was written — already prints a `NOTE` naming the dropped text.
- r1 (2026-08-09): blanket rejection in `cli.rs` §3.4a. Superseded.

## Goal

`--rename` asks Trim Galore to record clipped bases in the read ID. On **FASTQ input → uBAM output** that record cannot survive: `append_to_id` appends the `:clip5:`/`:clip3:` suffix after any header description, and BAM read names cannot contain whitespace, so `parse_name_and_data` discards the tail — annotation included. Refuse that combination so an explicitly-requested flag is never quietly neutralised. **uBAM input is untouched**, because there the loss cannot arise.

## Context

### The loss, and its exact boundary

Reproduced (`--rename --clip_R1 3` on `@withspace 1:N:0:ACGTAC`), and confirmed by both reviewers:

| Path | Result |
|---|---|
| FASTQ in → FASTQ out | `@withspace 1:N:0:ACGTAC:clip5:ACG` — annotation present |
| FASTQ in → uBAM out | QNAME `withspace` — **`:clip5:ACG` gone** |
| FASTQ in (no space in header) → uBAM out | annotation survives |
| **uBAM in → uBAM out** | **lossless, by SAM spec** — see below |

**Why uBAM input is safe by construction, not by luck.** A BAM read name cannot contain whitespace (`bam.rs:703` notes the spec constraint), so a record read from uBAM never carries a space-separated description. `append_to_id` splices before the first TAB — deliberately, per code-review finding C1 (`fastq.rs:171`) — so the aux-tag tail survives and the annotation lands on the name. Both reviewers verified this end-to-end **including aux tags**. r1 mis-framed this as data-dependent good fortune; it is a guarantee, and refusing it would withdraw a working path for nothing.

### Where the check must live — and the pre-emption it fixes

Format-gating means the check needs `any_bam`, which is computed at `main.rs:240-242` *after* `detect_input_format`. `Cli::validate` cannot know the input format, so §3.4a in `cli.rs` is the wrong home; the right one is `main.rs` beside **§3.4b**, the `--preserve-tags` rule, which is format-dependent for exactly the same reason.

This relocation resolves a defect both reviewers found in r1 for free: `cli.validate()` runs at `main.rs:217`, **before** format detection, so `--clump_only --rename` keeps its accurate mode-specific message (`cli.rs:851`, "preserves record contents byte-identically") instead of being pre-empted by a mechanical one. r1's Behavior 4 reasoned only about ordering *within* §3.4a and missed that the block as a whole precedes the clump-only check.

### Why rejection rather than repair (confirmed evidence)

`append_to_id` is reached from `trim_read` (`trimmer.rs:89`, calls at `:264`/`:273`), shared by both output formats, **and independently from six call sites in `specialty.rs`** (A §1.3; `:166`/`:219` among them) — so no fix inside it could be format-scoped. And the FASTQ path currently matches Perl exactly. Perl v0.6.11 (read from the `0.6.11` tag by all three of us) appends to the end of the whole ID line:

```perl
$l1 .= ":clip5:$clipped_sequence\n";   # :1098, :1115, :1282, :1298, :1399, :1414, :2209, :2243, :2260
$identifier .= ":clip5:$clipped_sequence";   # :1910, :2000 — the hardtrim sites
```

So splicing before a space would diverge from shipped Perl on a Perl-era flag. Two honest footnotes: Perl's own warning at `:3456` states the intent as `readname:clip5:GATC` — the read *name* — which its implementation never does; and the poly-A branch at `:1255` folds whitespace to `_` before appending, so Perl is not even internally uniform. `--rename` is **not** in the byte-identity CI matrix, so this is a judgement about user-visible stability, not a gate.

### What is actually wrong today (r1 overstated this)

**Not "silent".** #406 — merged before r1 was written — added `emit_description_dropped_once` (`bam.rs:1031-1045`), which prints a `NOTE` echoing the dropped text, and `CHANGELOG.md:32-42` documents it in the very section r1 appended to. What is silent is narrower and worth stating precisely: **there is no flag-specific diagnostic**. A user who passed `--rename` sees a general note about header descriptions and has no reason to connect it to their annotation being gone. r1's framing would have regressed the principle #406 established one changelog entry above.

### Scope also covers the hardtrim paths

`--hardtrim5`/`--hardtrim3` with `--rename` reach `append_to_id` via `specialty.rs` and lose the annotation the same way. The refusal covers them — a point in its favour that r1 neither stated nor tested.

### The docs errors (folded in, larger than r1 found)

`docs/.../guide/outputs.md`:
- **`:103` — the heading itself.** "## Renaming outputs" sits over filename prose; it is a filename heading, so fixing only the sentence beneath leaves the section mis-titled.
- **`:105` — the sentence.** "`--rename PREFIX` replaces the input filename stem" is wrong in every clause: `--rename` is a **boolean** (`cli.rs:302-303`) that appends `:clip5:SEQ` to read IDs and never touches filenames.
- **`--basename` is undocumented in the guide entirely** (B §6.2) — which answers r1's step-3 conditional: the stem-replacement description has no home to be moved to and needs one written.
- **§Feature compatibility (`:90-97`, seven bullets)** enumerates §3.4a's rejections and is missing `--rename`.

## Behavior

1. When `--output-format ubam` is set, `--rename` is set, **and at least one input is FASTQ**, the run is refused after format detection and before any I/O.
2. All-uBAM input is **accepted** and behaves exactly as today (annotation lands on the QNAME; aux tags intact).
3. `--clump_only --rename` continues to be refused earlier, by `Cli::validate`, with its own mode-specific message.
4. The message states **unrepresentability**, not a prediction of loss: a BAM read name cannot hold the annotation when the FASTQ header carries a description, so the flag cannot be honoured. No "silently dropped" claim — that is both untrue for uBAM input and stale for FASTQ input.
5. The remediation must not name a route that reproduces the loss. **`samtools import` drops the description identically, with and without `-T '*'`** (A §1.8), so "convert afterwards" is out. What remains true: take FASTQ output, or drop `--rename`.
6. Everything else unchanged — FASTQ-path IDs keep their exact Perl-matching format; nothing reaches `append_to_id`.

## Implementation outline

1. `main.rs`, beside §3.4b (after `any_bam` at `:242`) — add the gated refusal, in that block's established shape. Message along the lines of: `--rename cannot be honoured with --output-format ubam when any input is FASTQ: the :clip5:/:clip3: annotation is appended after the header description, and BAM read names cannot contain whitespace, so it cannot be represented in the output. Use FASTQ output, or drop --rename. (uBAM input is unaffected — its read names carry no description.)`
2. **Do not touch `tests/integration_ubam_out.rs:393-440`.** Under format-gating the C1 regression guard keeps passing as-is, because its input is uBAM. This is the single largest practical win of the format-gated choice and must be verified, not assumed (validation 3).
3. `cli.rs` tests + integration test — refusal with a FASTQ input; **acceptance with a uBAM input**; and `--rename` with FASTQ output still working (no over-reach).
4. `bam.rs:1027-1030` — its comment cites "#408 — `--rename`'s own `:clip5:` annotation" as an example of what the description notice covers. Under this change that example becomes unreachable (FASTQ input refused; uBAM input has no description). Reword so it does not point at a case that can no longer occur.
5. Docs — fix the `:103` heading and the `:105` sentence, give `--basename` a home in the guide, and add `--rename` to the §Feature compatibility list.
6. `CHANGELOG.md` — one bullet under `#### Bug fixes`, framed per Behavior 4: an explicitly-requested flag was being neutralised without a flag-specific diagnostic; it is now refused where it cannot be honoured, and untouched where it can. Say plainly that uBAM input is unaffected.
7. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`.

## Efficiency / Integration

Nil. One conditional on the startup path, beside a check that already computes the same `any_bam`.

## Assumptions

- **A1:** `--rename` is absent from §3.4a and the only existing `self.rename` rejection is the `--clump_only` one — verified by both reviewers.
- **A2 (CORRECTED — r1 had this false, and deferred verifying it):** `tests/integration_ubam_out.rs:393` **does** assert this combination succeeds. It is the C1 tab-splice regression guard. Under format-gating it keeps passing untouched; under r1's blanket form it would have failed and its coverage would have needed porting. Deferring this grep was not safe.
- **A3:** `--rename` is not in the Perl byte-identity matrix — verified.
- **A4:** FASTQ-path behaviour is unchanged, so Perl parity is untouched by construction.
- **A5 (new):** uBAM input cannot carry a header description, per the SAM spec — the guarantee the format gate rests on.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Refused where unrepresentable | `--rename --output-format ubam` with a FASTQ input | non-zero exit, message per Behavior 4; nothing written |
| 2 | **Accepted where lossless** | same flags with a uBAM input | succeeds; QNAME carries `:clip5:…`; aux tags intact |
| 3 | **The C1 guard survives untouched** | `cargo test ubam_out_rename_with_preserve_tags_keeps_tags_intact` with no edit to that test | passes |
| 4 | No over-reach | `--rename --clip_R1 3`, FASTQ output, space-bearing header | succeeds; ID still ends `:clip5:ACG`, byte-identical to before |
| 5 | Hardtrim covered | `--hardtrim5 20 --rename --output-format ubam` with FASTQ input | refused by the same guard |
| 6 | `--clump_only` keeps its own message | `--clump_only --rename --output-format ubam` | the byte-identity message from `cli.rs:851`, not the new one |
| 7 | Perl parity untouched | byte-identity harness + validation matrix | unchanged |
| 8 | Expected-fail control | validations 1 and 5 against the unpatched build | they *succeed* there and lose the annotation — the behaviour being removed |
| 9 | Docs describe real flags | read `:103`, `:105`, the compatibility list and the new `--basename` text against `cli.rs` | accurate; `--rename` present in the rejected list |
| 10 | No collateral | full `cargo test` + fmt + clippy | green |

## Questions or ambiguities

- **[Resolved — Felix]** Reject not repair; docs folded in; **format-gated scope**.
- **[Open, minor]** Whether to file a follow-up for carrying the annotation as a BAM aux tag. Both reviewers expect the refusal to invite the request; A recommends filing it. Recommend yes, after this lands.

## Implementation notes (2026-08-09, base `aa9f764`)

Implemented as specified. Surface: `CHANGELOG.md` +18, `docs/…/guide/outputs.md` +11/−5, `src/bam.rs` +3/−4 (doc comment only), `src/main.rs` +18, `tests/integration_ubam_out.rs` +172. **No output-producing code path changed** — the only executable addition is one `anyhow::bail!`.

### Deviations

1. **Step 3's "cli.rs tests" was dropped as unreachable, not skipped.** The guard must see `input_formats`, so it lives in `main.rs` — the `[[bin]]`. Library unit tests in `cli.rs` can only reach `Cli::validate()`, which by design does not carry this rule. All five new tests are therefore binary-driven, in `tests/integration_ubam_out.rs`. Validation 6 (`--clump_only` keeps its own message) is the one that would have been a `validate()` unit test; as an integration test it additionally proves the *ordering* between `validate()` and the new guard, which a unit test could not.
2. **`outputs.md`'s compatibility-list preamble had to change too.** It read "Rejected at CLI-validate time", true of all seven existing bullets and false of this one. Appending a bullet under it would have reproduced the defect class the change is correcting. Now "Rejected at startup, before anything is written", which is true of all eight.
3. **The `bam.rs` doc comment lost the `#408` cross-reference entirely** rather than being reworded to mention the now-refused case. Two live examples remain (instrument identifier, Illumina `1:N:0:INDEX`); the history is in the commit message per the repo's comment policy.
4. **The gate tests `!matches!(f, InputFormat::UnalignedBam)`, not a positive FASTQ match.** Follows the existing precedent at `main.rs:332`, and means a future third input format is refused (conservative) rather than silently admitted. A positive `matches!` on the two FASTQ variants would have inverted that default; an exhaustive `match` would have tripped clippy's `match_like_matches_macro` under `-D warnings`.

### Verification

| # | Validation | Result |
|---|---|---|
| 1 | Refused where unrepresentable | `rename_into_ubam_refused_for_fastq_input` — non-zero exit, message asserted, `sp_trimmed.bam` asserted absent |
| 2 | Accepted where lossless | `rename_into_ubam_accepted_for_ubam_input` — exit 0, QNAME contains `:clip5:`, `UB` tag present |
| 3 | C1 guard survives untouched | proven twice: `git diff -U0` shows a **single** hunk (`@@ -980,0 +981,172 @@`, pure append) and the guard's name appears nowhere in the diff; and it passes by name. Fixture `test_files/ubam_test_with_tags.bam` confirmed **git-tracked** (587 B), so the test's `SKIP` branch is not what passed |
| 4 | No over-reach | `rename_with_fastq_output_still_annotates_the_id` — `assert_eq!` on the full ID, `@withspace 1:N:0:ACGTAC:clip5:ACG` |
| 5 | Hardtrim covered | `rename_into_ubam_refused_for_hardtrim_fastq_input` |
| 6 | `--clump_only` keeps its own message | `clump_only_rename_keeps_its_own_message` — asserts `byte-identically` present **and** the new text absent |
| 7 | Perl parity untouched | By construction (A4): no FASTQ-output code changed. Validation 4 pins the ID format byte-for-byte. The Perl md5 harness itself is CI-only and was not run locally |
| 8 | Expected-fail control | Ran both refusals against an **unpatched build** (detached worktree at `aa9f764`, guard string count 0). Trim arm: exit 0, wrote `sp_trimmed.bam`, first QNAME `withspace` — annotation gone. Hardtrim arm: exit 0, wrote `sp.20bp_5prime.bam`, first QNAME `withspace`. **Additional check the plan did not ask for:** confirmed `--hardtrim5 20 --rename` on the *FASTQ* path really does annotate (`:clip5:ACGTACGT`), so validation 5 refuses something that genuinely worked rather than a no-op |
| 9 | Docs describe real flags | `--rename` re-documented as the boolean it is under a new "Annotating read IDs" heading; `--basename` given its own section (it had none in the guide); compatibility bullet added |
| 10 | No collateral | `cargo fmt --all -- --check` clean (after one rustfmt rewrite of two `assert!` calls); `cargo clippy --all-targets --release -- -D warnings` exit 0; `cargo test` **575 passed / 0 failed** across 14 targets. The same suite on the unpatched control gave **570 passed / 0 failed**, so the delta is exactly the five tests added — counted rather than assumed, since equal counts are how a stale-tree verification hides |

### Iteration log

- **#1** First `fmt --check` failed on two `assert!` calls rustfmt wanted expanded. Ran `cargo fmt --all`; re-check clean. No semantic change.
- **#2** The by-name grep for the five new tests appeared to show `rename_into_ubam_accepted_for_ubam_input` missing. It had in fact run — the test binary's own stdout interleaved on the line, breaking the `^test` anchor. Re-ran with `-- --exact`: 1 passed, 32 filtered out. A grep artefact, not a missing test.
- **#3** `cargo clippy … | tail` reported success with `CLIPPY EXIT=` empty — `${PIPESTATUS[0]}` is bash, and zsh silently yielded nothing. Re-ran unpiped: exit 0, but in **0.18 s** with no mention of the edited target, i.e. a cache hit. Injected a deliberate `unused_variable` into the new test block; clippy flagged it at `:1001`, proving the target *is* linted. Removed the probe, re-ran with `-D warnings` (cache now invalid, 2.90 s, recompiled): exit 0, zero warnings.

## Self-Review (r2)

- **What r1 got wrong, owned:** it deferred A2's grep and A2 was false — the assumption I explicitly flagged as unverified was the one that broke the patch, for the second time this session. It mis-framed a spec guarantee as luck, which is what made the blanket scope look free. Its remediation advice named `samtools import`, which reproduces the loss. And it called the loss "silent" one merge after we shipped the notice that reports it — regressing #406's own principle in the changelog entry directly below it.
- **What the format gate bought beyond correctness:** the C1 regression guard survives with no edit, and the `--clump_only` pre-emption disappears because the check moves behind `validate()`. Neither was an argument for it; both are consequences worth recording.
- **Traps checked:** validation 3 pins the C1 guard explicitly rather than trusting that it still passes; validation 8's control is informative in the unusual direction (pre-patch the run *succeeds* while losing data); validation 2 exists so the gate cannot silently become blanket.
- **Remaining risk:** the message is prose and reviewable. The change cannot alter any output path.
