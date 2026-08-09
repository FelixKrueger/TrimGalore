# Code Review A — `fix/408-rename-ubam-rejection` @ `3f1b0fe`

**Reviewer:** A (independent; no shared state with Reviewer B)
**Target:** `git diff dev...HEAD` — one commit, `3f1b0fe`
**Base:** `aa9f764` (dev)
**Date:** 2026-08-09
**Constraint honoured:** no fixes applied, no tracked file modified, no branch switch. Control build made
via `git archive aa9f764 | tar -x -C $TMPDIR/tg408_ctrl` (no worktree, no `.git` mutation).

---

## Summary

The change is small, correctly sited, and does what it says: one `anyhow::bail!` on the startup path,
plus docs/CHANGELOG/test surface. No output-producing code path is touched. I independently confirmed
the core mechanism, the siting property ("nothing written"), the uBAM-input acceptance path (annotation
really does reach the QNAME, in SE *and* interleaved PE), the hardtrim arm, the `--clump_only`
ordering, and the plan's un-asked-for claim that `--hardtrim5 20 --rename` genuinely annotates on the
FASTQ path. The `bam.rs` doc-comment reasoning is sound. The docs `--basename` / `--rename` sections
match `cli.rs`. Every factual claim I could check in the CHANGELOG entry is true, including the
`#406`-notice claim.

Two findings are worth acting on, both about *diagnostics* rather than correctness:

1. **The refusal pre-empts two more-specific structural errors** (`--paired` with a single FASTQ;
   `--paired` with a mixed FASTQ+uBAM pair). The user gets a lecture about `:clip5:`
   representability when their real problem is a forgotten R2 or a mis-typed filename. This is the
   *same* defect class the plan congratulates itself for avoiding with `--clump_only` — and the fix is
   the same move: site the gate a few lines later, still ahead of `ensure_output_dir`. **CONFIRMED by
   running both invocations with and without `--rename`.**
2. **The message's headline clause is false for a describable share of real input.** A FASTQ whose
   headers carry no description loses nothing today — I ran it on the control build and the QNAME came
   out `plainhdr:clip5:ACG` — yet the message opens with "`--rename` cannot be honoured … when any
   input is FASTQ". The body hedges correctly ("whenever a FASTQ header carries text after the first
   space"), so the message contradicts itself within four lines. The *gate* is defensible (per-record
   knowledge isn't available up front); the *wording* asserts more than the code knows.

Everything else is nits.

---

## Verification performed

| # | What I ran | Result |
|---|---|---|
| V1 | `cargo build --release` on the branch | clean |
| V2 | control build at `aa9f764` from `git archive` | clean; guard-string count in `src/main.rs` = 0 |
| V3 | patched: `--rename --clip_R1 3 --output-format ubam sp.fastq` | exit 1, refusal message |
| V4 | patched, `-o $D/should_not_exist`, no other args | exit 1; **cwd unchanged, `-o` dir not created** |
| V5 | patched: `--clip_R1 5 --rename --output-format ubam test_files/ubam_test_with_tags.bam` | exit 0; QNAMEs `SRR24827378.1:clip5:AATTA` … ; `UB` present |
| V6 | patched: `--paired --clip_R1 5 --rename --output-format ubam test_files/ubam_paired_test.bam` | exit 0; interleaved out, R1 `…:clip5:AATTA` flag `0x4d`, R2 unannotated `0x8d` — correct (only R1 clipped) |
| V7 | patched: mixed FASTQ+uBAM, both argument orders | exit 1 both ways |
| V8 | patched: `--hardtrim3 20 --rename --output-format ubam` | exit 1, same message (hardtrim3 arm covered, not just hardtrim5) |
| V9 | control: description-**free** FASTQ + `--rename --clip_R1 3 --output-format ubam` | **exit 0, QNAME `plainhdr:clip5:ACG` — lossless** |
| V10 | control: `--rename` with **no clip flag** + ubam out | exit 0, QNAME `plainhdr` — annotation never attempted |
| V11 | control: `--hardtrim5 20 --rename` FASTQ out | `@withspace 1:N:0:ACGTAC:clip5:ACGTACGT` — plan's claim **confirmed** |
| V12 | control: FASTQ-with-description + `--rename` + ubam out, stderr | `NOTE: … First dropped: "1:N:0:ACGTAC:clip5:ACG"` — the #406 notice does fire and does echo the annotation |
| V13 | `test_files/ubam_test_with_tags.bam` | exists, git-tracked, 587 B → C1 guard's `SKIP` branch is not what passes |
| V14 | full `cargo test`, `cargo fmt --check`, `cargo clippy --all-targets --release -- -D warnings` | 575 passed / 0 failed across 14 targets; fmt exit 0; clippy exit 0 and genuinely recompiled — see "Test suite" below |

BAM QNAMEs were read with a standalone BGZF/BAM parser (python, `zlib` wbits=31 multi-member), not via
the crate under review, so the read side is independent of the code being reviewed.

---

## Logic

### L1 — CONFIRMED — The gate pre-empts two better error messages (Medium)

`src/main.rs:316-333` sits **before** the two structural format checks at `src/main.rs:350-356`
(`--paired` with a single non-BAM input) and `src/main.rs:377-392`
(`reject_bam_format_mismatch_in_pair`, issue #363).

Observed, `--output-format ubam` throughout:

```
$ trim_galore --paired --rename --clip_R1 3 --output-format ubam sp.fastq
Error: --rename cannot be honoured with --output-format ubam when any input is FASTQ. …

$ trim_galore --paired          --clip_R1 3 --output-format ubam sp.fastq
Error: --paired with a single input file is only legal if that file is a uBAM.
       Provide R1 and R2 as separate files for FASTQ paired-end mode.
```

```
$ trim_galore --paired --rename --clip_R1 3 --output-format ubam sp.fastq ubam_paired_test.bam
Error: --rename cannot be honoured with --output-format ubam when any input is FASTQ. …

$ trim_galore --paired          --clip_R1 3 --output-format ubam sp.fastq ubam_paired_test.bam
Error: --paired requires both inputs of a pair to be the same format. Got mixed: sp.fastq is
       FASTQ (plain) and …ubam_paired_test.bam is uBAM. Pass two FASTQ files, or a single
       interleaved uBAM. If you meant two FASTQ files, check for a mis-typed filename.
```

In both cases the run is broken *independently* of `--rename`, and the pre-empted message is the one
that names the actual defect — including a mis-typed-filename hint that the #363 author wrote
deliberately. A user has to drop a flag they wanted, re-run, and only then learn what is really wrong.

This is precisely the failure the plan celebrates avoiding for `--clump_only` (`PLAN.md:35`, "keeps its
accurate mode-specific message … instead of being pre-empted by a mechanical one"). The same argument
applies one level down and was not carried through.

**Recommendation:** move the block to immediately after `reject_bam_format_mismatch_in_pair` (i.e.
after `src/main.rs:392`). That is still ahead of `naming::ensure_output_dir` at `src/main.rs:406` and
ahead of every dispatch arm, so the "nothing written" property — which I verified at V4 — is
untouched. The `--clump_only` ordering test keeps passing (that guard is in `validate()`, far
earlier). Nothing else reads `cli.rename` in between.

### L2 — CONFIRMED — The refusal withdraws a path that is lossless today (Medium)

Two sub-cases, both verified against the control build:

**(a) FASTQ with no header description.** `append_to_id` (`src/fastq.rs:173-181`) appends to the end;
`parse_name_and_data` (`src/bam.rs:709`) splits on the first `['\t', ' ']`. With no space and no tab
there is nothing to split, so the whole annotated ID becomes the read name:

```
control @ aa9f764:  @plainhdr + --rename --clip_R1 3 --output-format ubam
                 →  QNAME  plainhdr:clip5:ACG      (exit 0, annotation intact)
patched @ 3f1b0fe:  same invocation → exit 1, refused
```

That is not an exotic shape — SRA/ENA `@SRR…​.N`, `samtools fastq`-derived FASTQ, and most synthetic
or simulator output have no description at all. For those files the flag was being honoured exactly
as asked.

**(b) `--rename` with no clipping flag.** `trimmer.rs:258-275` only calls `append_to_id` inside
`if let Some(n) = clip_5 / clip_3`, so without `--clip_R1/R2` or `--three_prime_clip_R1/R2` (and
outside hardtrim) `--rename` never touches the ID:

```
control: --rename --output-format ubam plain.fastq → exit 0, QNAME plainhdr (nothing appended)
patched: same → exit 1
```

A wrapper that threads `--rename` unconditionally into every invocation (the flag's stated use is UMI
bookkeeping, so pipeline-level defaults are plausible) now hard-fails on runs where it previously did
nothing at all.

`PLAN.md:26` records sub-case (a) in its boundary table and Behavior 1 chose format-level gating
anyway, so this is a *recorded* decision, not an oversight — I am flagging the consequence, not
claiming it was missed. Sub-case (b) is not in the plan at all.

I do **not** recommend making the gate per-record or per-header: reading the first record to decide
would make the refusal data-dependent and untestable, and letting the run proceed until the first
description appears would leave a partly-written BAM. The format-level gate is the right shape. What
needs fixing is what the message *asserts* — see M1.

### L3 — clean — Siting w.r.t. writers

Verified by reading `src/main.rs:229-480` end to end. Between `cli.validate()` (`:229`) and the new
block (`:316`) the only file access is `sanity_check_any` (`:241`) and `detect_input_format` (`:247`),
both read-only. The first write of any kind is `naming::ensure_output_dir` (`:406`), then the specialty
pre-flights and dispatch arms (`:450` onward). Empirically confirmed at V4: after a refusal the cwd
contains only the input, and an absent `-o` directory is **not** created.

Every arm that can reach a writer is downstream: `--hardtrim5` (`:450`), `--hardtrim3`, `--clock`,
`--implicon`, `--clump_only`, `run_ubam_output`, and the SE/PE trim dispatch. `--clock` / `--implicon`
never reach the gate in a way that matters — `cli.rs:606-616` rejects them with `--output-format ubam`
inside `validate()`, i.e. earlier. `--clump_only --rename` likewise (`cli.rs:851`).

### L4 — clean — The gate predicate

`input_formats.iter().any(|f| !matches!(f, InputFormat::UnalignedBam))` is the right test for the
three variants that exist (`FastqPlain`, `FastqGz`, `UnalignedBam`, `src/format.rs:31-40`):

- **Mixed FASTQ + uBAM** → refused, both argument orders (V7). Correct: SE mode processes each input
  independently, so the FASTQ member would silently lose the annotation while the BAM member kept it.
  Refusing the whole run matches how `any_bam` gates `--phred64` at `:276`.
- **Single interleaved uBAM under `--paired`** → accepted (V6), and the annotation really lands on the
  QNAME of the R1 half. `input_formats` has exactly one element and it is `UnalignedBam`, so `.any()`
  is false. Correct.
- **Empty input** → unreachable; `input` is `required = true`, so `.any()` on an empty slice can't be
  hit.
- **A future fourth variant** → refused, because the predicate is negative. Deviation note 4 in
  `PLAN.md:124` calls that conservative and cites `main.rs:332` as precedent; the live precedent is
  now `main.rs:350` (`!matches!(input_formats[0], InputFormat::UnalignedBam)`), same shape. I agree
  with the direction. Worth knowing that if the fourth variant were CRAM or BINSEQ — both of which
  also forbid whitespace in read names — the gate would refuse a lossless path. That is the safe
  failure and needs no change now.

### L5 — clean — `bam.rs` doc-comment reasoning (`src/bam.rs:1026-1029`)

The removed example ("`--rename`'s own `:clip5:` annotation") is genuinely unreachable now.
`emit_description_dropped_once` fires only from `parse_name_and_data`'s space branch
(`src/bam.rs:717-721`), which requires a space in `FastqRecord::id`. For that space to coexist with a
`--rename` annotation you need FASTQ input (refused) or a uBAM record carrying a description
(impossible per SAM spec, and `BamReader` builds the id from `name()` + tab-joined tags). The two
remaining `--rename`-adjacent routes are also closed: `--clump_only --rename` is rejected at
`cli.rs:851`, and `--clock`/`--implicon` — whose `append_to_id` calls at `specialty.rs:312/413` are
unconditional, not `--rename`-gated — are rejected with `--output-format ubam` at `cli.rs:606-616`.
Sound.

---

## Efficiency

Clean. One `Iterator::any` over a slice that is at most a handful of elements and was already
materialised for the `--phred64` and `--preserve-tags` rules. No allocation, no I/O, short-circuits on
`cli.rename` first. Nothing to say.

---

## Errors / message accuracy

### M1 — CONFIRMED — The message asserts unrepresentability where the code proves representability (Medium)

`src/main.rs:325-331`. Clause by clause:

| Clause | True? |
|---|---|
| "`--rename` cannot be honoured with `--output-format ubam` when any input is FASTQ" | **No** — V9/V10. Description-free FASTQ honours it exactly; `--rename` without a clip flag has nothing to honour |
| "The `:clip5:`/`:clip3:` annotation is appended to the end of the read ID" | Yes for FASTQ input (`fastq.rs:179`). Slightly loose in general — with a tab tail it splices *before* the tab (`fastq.rs:177`) — but the sentence is scoped to the FASTQ case it is describing |
| "so whenever a FASTQ header carries text after the first space the annotation lands in that text" | Yes — CONFIRMED V12, dropped text was `1:N:0:ACGTAC:clip5:ACG` |
| "and BAM read names cannot contain whitespace" | Yes, SAM spec; enforced by `parse_name_and_data`'s split |
| "so none of that tail can be represented in the output" | Yes |
| "Use FASTQ output, or drop `--rename`." | Correct, and correctly avoids `samtools import` |
| "uBAM input is unaffected: BAM read names carry no description, so there the annotation lands on the name itself." | Yes — CONFIRMED V5/V6 |

So the only defect is the headline sentence, and it is contradicted by the message's own next sentence:
"cannot be honoured … when any input is FASTQ" versus "**whenever** a FASTQ header carries text after
the first space". A user whose headers are bare will read this, check their headers, find no spaces,
and conclude the tool is wrong about their data.

The plan's Behavior 4 (`PLAN.md:69`) asks for exactly this property — "states **unrepresentability**,
not a prediction of loss" — and the headline breaks it in the other direction: it *over*-states,
asserting impossibility for inputs where the annotation demonstrably survives.

**Recommendation:** make the headline say what the check actually is — a format-level refusal, because
representability is per-header and cannot be known before the run. E.g. lead with "`--rename` is
refused with `--output-format ubam` when any input is FASTQ, because whether the annotation can be
represented depends on each header and cannot be established up front", then keep the existing
mechanism and remediation sentences unchanged. That is honest about the gate's granularity, keeps every
verified clause, and pre-empts the "but my headers have no spaces" bug report.

### M2 — no findings

No unhandled conditions, no panics, no security surface. `anyhow::bail!` on the startup path, exit
code 1 confirmed. No unwrap/expect added to library code.

---

## Structure / style

### S1 — Low — the two uBAM-acceptance tests are near-duplicates

`tests/integration_ubam_out.rs:1035` (`rename_into_ubam_accepted_for_ubam_input`, new) runs the
*identical* invocation to `:394` (`ubam_out_rename_with_preserve_tags_keeps_tags_intact`,
pre-existing) — same fixture, same `--clip_R1 5 --rename --output-format ubam --preserve-tags CB,UB`.
Only the assertions differ (annotation reaches the QNAME vs `UB` value not polluted). Both are worth
asserting; two full binary invocations to do it is duplication. Not worth churning the C1 guard for —
if anything, a doc-comment cross-reference between the two would stop a future reader deleting one as
redundant.

### S2 — clean — comment policy

The new comment (`src/main.rs:316-317`) is two lines, states the current-state fact ("sited here …
because the rule is format-gated"), narrates no prior state, and points at the sibling rule rather
than restating its rationale. Compliant, and notably leaner than the 20-line `--phred64` block
directly above it. The `bam.rs` edit removes text; nothing added.

### S3 — clean — idiom and naming

Fully-qualified `trim_galore::cli::OutputFormat::UBam` matches the surrounding usages at `:305`,
`:384`, `:412`. `matches!` over an exhaustive `match` matches the file's prevailing style (and, per
`PLAN.md:124`, avoids `clippy::match_like_matches_macro` under `-D warnings`). No new imports, no dead
code.

---

## Test quality

None of the five is vacuous. Each was checked for the specific failure modes named in the brief.

| Test | Verdict | Notes |
|---|---|---|
| `rename_into_ubam_refused_for_fastq_input` (`:1000`) | **Sound** | Asserts non-zero, message substring, **and** `sp_trimmed.bam` absent. Non-vacuous: the control build writes that file and exits 0 (V9-class). 28 bp seq, `--clip_R1 3` → 25 bp, clears default `--length 20`; qual is `I`×n so `-q 20` can't interfere. |
| `rename_into_ubam_accepted_for_ubam_input` (`:1035`) | **Sound, and load-bearing** | This is the test that stops the gate silently becoming blanket. I reproduced its assertions by hand (V5). No `SKIP` branch — a missing fixture would fail, not pass. Runs from crate root, so the relative `test_files/` path resolves. |
| `rename_with_fastq_output_still_annotates_the_id` (`:1074`) | **Sound** | `assert_eq!` on the whole ID (`@withspace 1:N:0:ACGTAC:clip5:ACG`) is the strongest form available and pins the Perl-matching format byte-for-byte. Reads `sp_trimmed.fq` (plain input → plain output), correct. |
| `rename_into_ubam_refused_for_hardtrim_fastq_input` (`:1098`) | **Sound, slightly weaker than its sibling** | Non-vacuous — control writes `sp.20bp_5prime.bam` and annotates (V11). It asserts non-zero + message but **not** that nothing was written, unlike the trim arm. See T1. |
| `clump_only_rename_keeps_its_own_message` (`:1126`) | **Sound** | Asserts the mode-specific text present *and* the new text absent — a genuine ordering assertion that a `validate()` unit test could not make. |

### T1 — Low — the hardtrim refusal test does not assert absence

`tests/integration_ubam_out.rs:1098-1124`. The trim arm asserts `!dir.join("sp_trimmed.bam").exists()`;
the hardtrim arm asserts nothing about the filesystem. The specialty arm is the one where the guard's
position matters most — `naming::preflight_output_collisions` and the writers sit at `main.rs:452+`,
downstream of the gate, and a future re-siting of the gate (including the L1 move I recommend) would
be caught by the trim arm but not here. Adding `assert!(!dir.join("sp.20bp_5prime.bam").exists())` —
the exact filename the control build produces — would close it for one line.

### T2 — Low — "nothing written" is asserted narrowly

The trim arm checks one filename. The run would also produce `sp.fastq_trimming_report.txt` /
`.json` if it got that far. Asserting the directory contains exactly `sp.fastq` would state the
property the comment claims ("the refusal must precede every writer") rather than a proxy for it. I
verified the strong form by hand (V4) — the cwd stays clean and an absent `-o` dir is not created — so
this is test-expressiveness, not a gap in behaviour.

### T3 — The untouched C1 guard is genuinely still exercising its claim

`tests/integration_ubam_out.rs:393-440`, unmodified by this commit (`git show --stat` lists
`tests/integration_ubam_out.rs | 172 +++`, all additions; `git show 3f1b0fe -- tests/` shows a single
`@@ -978,3 +978,175 @@` hunk, a pure append, and the guard's name appears nowhere in the diff).

Its `SKIP` branch (`:400-403`) requires `test_files/ubam_test_with_tags.bam` to be missing. The file
exists, is git-tracked, 587 B (V13), so `SKIP` is not what passes. Its input is uBAM, so the new gate
does not fire — and it should not: I confirmed the same invocation succeeds and the annotation reaches
the QNAME (V5). The guard still asserts what it claims (`UB` value free of `:clip5:`).

---

## Docs

`docs/src/content/docs/guide/outputs.md`

### D1 — clean — the `--basename` section (`:106`) is accurate

Checked against `cli.rs:156-159` (`Option<String>`, "replaces input filename stem"), `cli.rs:640-644`
and `cli.rs:667-671` (SE >1 file and PE >1 pair both refused, "ambiguous output naming"), and
`io.rs:326`/`:357`/`:462` (`BASE_val_1` / `BASE_val_2` naming, with the Perl v0.6.5+ note). Every
clause holds, including "longer input lists are refused, because the output naming would be
ambiguous".

### D2 — clean — the `--rename` section (`:110`) is accurate

Matches `cli.rs:300-301` almost word for word ("Appends `:clip5:SEQ` and/or `:clip3:SEQ` to the read
ID (each half only when that side was clipped)"), and "does not affect filenames" is confirmed by
there being no `cli.rename` reference anywhere in `io.rs`. Naming the three flag families that
actually drive it (`--clip_R1/R2`, `--three_prime_clip_R1/R2`, `--hardtrim5/3`) is right —
`--clock`/`--implicon`'s ID appends at `specialty.rs:312/413` are unconditional and not `--rename`'s
doing, so leaving them out is correct rather than an omission.

### D3 — clean — the compatibility-list preamble (`:90`) is true of all eight bullets

"Rejected at startup, before anything is written (v1 scope)". Seven bullets are `validate()`-time
(`cli.rs:587-632`). The eighth is `main.rs:316`, and I verified empirically (V4) that a refused run
writes nothing and does not even create an absent `--output_dir`. Deviation note 2 in `PLAN.md:122` is
the right call — appending under the old "Rejected at CLI-validate time" would have reproduced the
defect class being fixed.

### D4 — Low — the new bullet carries M1's overclaim, and breaks the list's register

`:98`. "**only when at least one input is FASTQ.** … so a FASTQ header carrying text after the first
space puts the annotation inside that text — and BAM read names cannot contain whitespace, so none of
that tail survives." Accurate as written, but like the error message it never says that a
description-free FASTQ is refused too, which is the case a reader will hit and be puzzled by. Should
track whatever wording M1 settles on.

Separately: every sibling bullet is one clause; this one is four sentences plus a meta-note about
*where* the check lives. The implementation detail ("checked after input-format detection rather than
at CLI-validate time") is now redundant with the preamble the same commit rewrote — the preamble no
longer claims CLI-validate time, so there is nothing left to except this bullet from. Trimming the
bullet to its user-visible rule and dropping the siting note would match the list.

### D5 — Low — `--rename`'s own `--help` text does not mention the restriction

`cli.rs:300-301` is unchanged. Compare `--preserve-tags` (`cli.rs:249-250`), which spells out "a hard
error when combined with `--output-format ubam`" right in its help. `--output-format`'s own help
(`cli.rs:259-260`) only says "some flag combinations are rejected (see `--help` and the startup
diagnostics for details)", so there is no enumerated list to add to — which leaves `--rename`'s help
as the only place a `--help` reader would look, and it is silent. One clause there would match the
established precedent.

---

## CHANGELOG

`CHANGELOG.md:8-24`. I checked every factual claim; all hold.

| Claim | Verdict |
|---|---|
| "`--rename` appends `:clip5:SEQ`/`:clip3:SEQ` to the end of the read ID" | True on the FASTQ path (`fastq.rs:179`). Same mild looseness as M1's second clause (tab tails splice earlier); acceptable in a changelog |
| "When a FASTQ header carries text after the first space the annotation lands inside that text … the whole tail was discarded" | **CONFIRMED** V12 — dropped text was `1:N:0:ACGTAC:clip5:ACG` |
| "now refused when at least one input is FASTQ, including on the `--hardtrim5`/`--hardtrim3` paths" | **CONFIRMED** for both — V8 covers `--hardtrim3`, which no test exercises (see T-note below) |
| "**uBAM input is unaffected and keeps working** … annotation lands on the name itself and preserved aux tags round-trip intact" | **CONFIRMED** V5 (SE, `UB` intact) and V6 (interleaved PE) |
| "The header-description notice added for #406 did report the dropped text" | **CONFIRMED** V12 — the notice fires and echoes the annotation verbatim. This is the claim most at risk of being wishful, and it is true |
| "a user who passed `--rename` had no reason to connect a general note about header descriptions to their own annotation being gone" | Fair characterisation of the V12 stderr |
| "it described `--rename` as `--rename PREFIX` replacing the output filename stem, which is `--basename` — a flag the guide did not document at all" | Substantively true. Minor misquote: the old line said "replaces the **input** filename stem in the output names", not "the output filename stem". Immaterial |

The entry does not overclaim "silent" anywhere, consistent with the plan's Behavior 4. Good.

One coverage note rather than an accuracy one: the changelog asserts `--hardtrim3` coverage and only
`--hardtrim5` has a test. The gate is mode-independent so the risk is nil, and V8 confirms the
behaviour — but the asserted arm is untested.

---

## Test suite, fmt, clippy — all verified on this branch

**`cargo test` from the crate root: 575 passed, 0 failed, 0 ignored, across 14 targets.** Matches the
stated baseline exactly. Per-target counts, summed rather than eyeballed:
`411 + 0 + 15 + 12 + 15 + 11 + 1 + 8 + 48 + 11 + 2 + 8 + 33 + 0 = 575`. The two zero-count targets are
empty harnesses, not filtered-out runs.

Each new test confirmed to actually execute, by exact name (the plan's iteration-log trap #2 —
interleaved stdout breaking a `^test` grep anchor — avoided by using `-- --exact` and reading the
`passed` count, not the "ok"):

```
rename_into_ubam_refused_for_fastq_input             1 passed; 32 filtered out
rename_into_ubam_accepted_for_ubam_input             1 passed; 32 filtered out
rename_with_fastq_output_still_annotates_the_id      1 passed; 32 filtered out
rename_into_ubam_refused_for_hardtrim_fastq_input    1 passed; 32 filtered out
clump_only_rename_keeps_its_own_message              1 passed; 32 filtered out
ubam_out_rename_with_preserve_tags_keeps_tags_intact 1 passed; 32 filtered out   ← the C1 guard
```

33 total in `integration_ubam_out`, consistent with the suite line above. The C1 guard run with
`--nocapture` prints **zero** `SKIP` lines, so it is passing on its real assertions (T3).

`cargo fmt --all -- --check` → exit 0.

`cargo clippy --all-targets --release -- -D warnings` → exit 0, zero `warning`/`error` lines. This was
a genuine compile, not a cache hit: the run emitted `Compiling trim-galore v2.3.0` and took 2.48 s
(the repo's known trap is a sub-second "success"; a second, cached invocation is what would be
instant).

---

## Recommendations, by priority

**Critical** — none.

**High** — none.

**Medium**

1. **L1 — move the gate behind the two structural format checks.** Relocate `src/main.rs:316-333` to
   just after `reject_bam_format_mismatch_in_pair` (`:392`). Restores the more useful message for
   `--paired`-with-one-FASTQ and mixed-pair invocations; keeps every property the current siting has
   (still ahead of `ensure_output_dir` at `:406` and every dispatch arm — V4 property preserved);
   costs nothing. This is the same move the plan already made once, for `--clump_only`.
2. **M1 — fix the message's headline clause.** It currently asserts impossibility for inputs where the
   annotation demonstrably survives (V9, V10), and contradicts its own following sentence. Reframe as
   a format-level refusal whose granularity is stated, keeping the mechanism and remediation text
   unchanged.
3. **L2 — decide consciously about the withdrawn-but-lossless cases and record it.** Description-free
   FASTQ (recorded in the plan) and `--rename` without any clip flag (not recorded) both worked before
   and are refused now. I think the format-level gate is right and neither should change behaviour —
   but sub-case (b) should be in the plan's boundary table, and if a bare-header user opens an issue
   the answer should already be written down.

**Low**

4. **T1** — assert `!sp.20bp_5prime.bam` in the hardtrim refusal test, matching the trim arm.
5. **D5** — add one clause to `--rename`'s `--help` text (`cli.rs:300-301`) naming the
   `--output-format ubam` restriction, matching `--preserve-tags`' precedent.
6. **D4** — trim the new docs bullet to its user-visible rule; the "checked after input-format
   detection rather than at CLI-validate time" note is now redundant with the preamble the same commit
   rewrote.
7. **T2** — strengthen "nothing written" from one filename to the whole directory listing.
8. **S1** — cross-reference the two near-identical uBAM-acceptance tests so neither looks redundant.
9. Optional: a one-line `--hardtrim3` case, since the CHANGELOG names that arm.

**Not recommended**

- Making the gate per-header or per-record. Data-dependent refusals are worse than a conservative
  format-level one, and a mid-stream refusal would leave a partial BAM.
- Touching `tests/integration_ubam_out.rs:393-440`. The C1 guard is correct as-is and still exercises
  its claim (T3).
