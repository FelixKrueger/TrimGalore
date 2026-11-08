# Code Review B — `--phred64` pre-flight widened to uBAM output

**Reviewer:** B (independent, fresh context)
**Target:** uncommitted changes on `dev`, base `49152c1`
**Files:** `src/main.rs`, `tests/integration_fastqc_phred64.rs`, `CHANGELOG.md`,
`docs/src/content/docs/guide/{quality,outputs}.md`
**Date:** 2026-11-08
**Mandate:** recommend only — no edits to `src/`, `tests/`, `docs/` or `CHANGELOG.md`.

---

## Verdict

**The code is correct.** The boolean is right, it is sited ahead of every dispatch branch and ahead
of `ensure_output_dir`, it does not break any currently-working invocation, and the four new tests
pin the pre-flight's own refusal rather than any non-zero exit. I found no functional defect at any
severity.

Every finding below is in **user-facing prose**, and one of them is substantive: the reworded
startup `NOTE` **replaced accurate wording with an overclaim**, in the same diff whose `PLAN.md` §13
records tightening the CHANGELOG away from that exact overclaim. Two secondary prose issues follow
from the same rewrite.

| # | Severity | Area | One line |
|---|---|---|---|
| B1 | **High** | `src/main.rs:567-569` | The `NOTE` now claims "a mislabelled input is refused"; the wording it replaced was byte-accurate |
| B2 | **Medium** | `src/main.rs:564-566` | The `NOTE`'s mode list now reads as exhaustive and omits `--clock` / `--implicon`, which are in the same class |
| B3 | Low | `docs/.../quality.md:37` | Same overclaim, mitigated by a forward reference to the caution |
| B4 | Low | `CHANGELOG.md:33` | `docs/guide/quality.md` is not a resolvable path (precedent exists at line 180) |
| B5 | Low | `src/main.rs:381-387` | The comment block is now 7 lines; the diff's own 2 are compliant, the pre-existing 5 are not |
| B6 | Low | tests | No acceptance case on the arm the diff adds; the #362 guards cover it from other files |
| B7 | Nit | `CHANGELOG.md:35-40` | The `#500` entry two bullets away is now stale about where the refusal happens |

---

## What I verified, and how

All runs used the existing `target/release/trim_galore` (mtime 10:09, newer than `src/main.rs` at
10:04, so built from this working tree). **I did not rebuild, and did not run `cargo test` or
`cargo clippy`** — both would have raced the two concurrent agents. I ran
`cargo fmt --all -- --check` (no `target/` writes): **exit 0**. Scratch:
`/tmp/claude-501/revB-phred64-76065`.

### Siting

The new gate is `src/main.rs:388-391`, under the comment at 381-387. Line numbers of everything downstream:

| Site | Line |
|---|---|
| new `--phred64` input pre-flight | **388** |
| `naming::ensure_output_dir` | 551 |
| the startup `NOTE` | 558 |
| `--hardtrim5` dispatch + `return` | 603 / 632 |
| `--hardtrim3` dispatch + `return` | 634 / 663 |
| `--clock` / `--implicon` | 665 / 687 |
| `--clump_only` dispatch + `return` | 719 / 970 |
| `run_ubam_output` | 988 |

Nothing is written before 388. Confirmed empirically — in every refusal below the output directory
was never created:

| Invocation (Phred+33 input, `--output-format ubam`) | exit | output dir |
|---|---|---|
| `--phred64 -q 0 --length 0` | 1 | absent |
| `--phred64 --hardtrim5 20` | 1 | absent |
| `--phred64 --hardtrim3 10` | 1 | absent |
| `--phred64 --clump_only` | 1 | absent |
| `--phred64 --paired R1 R2` (R1 clean, R2 low) | 1 | absent, blames `pe_R2.fastq` |
| gzipped input `p33.fastq.gz` | 1 | absent |

The paired and gzip cases matter because `lowest_quality_char` goes through
`open_sync_reader(path, &[])`, so compression and BGZF framing are handled by the same factory the
rest of the binary uses.

### Nothing currently working breaks

The four #362 guards all use `test_files/phred64_test{,_R1,_R2}.fastq`, whose lowest quality byte
is `B` = **ASCII 66** — above the floor, so the pre-flight accepts them. Re-ran all four shapes
against the current binary:

| Guard shape | exit | output |
|---|---|---|
| SE trim → uBAM | 0 | `phred64_test_trimmed.bam` + both reports |
| `--hardtrim5 20` → uBAM | 0 | `phred64_test.20bp_5prime.bam` |
| PE two-file trim → uBAM | 0 | `phred64_test_R1_val.bam` + 4 reports |
| `--clump_only` → uBAM | 0 | `phred64_test_clumped.bam` + clumping report |

Also confirmed the deliberately-unchecked path is untouched: `--phred64 -q 20 --length 20` on
Phred+33 input with FASTQ output still exits 0.

### The boolean

`cli.phred64 && (cli.fastqc_requested() || output_format == UBam)` is complete for its stated scope:

- `fastqc_requested()` is `self.fastqc || self.fastqc_args.is_some()` — both FastQC entry points.
- `OutputFormat::UBam` is the only way a BAM is written; every BAM writer call site is gated on it.
- uBAM **input** can never reach this check — `any_bam` bails at line 369, ahead of it.
- No false positives in the uBAM arm: with uBAM input excluded, every surviving input is FASTQ and
  `BamWriter::write_record` always does `b.saturating_sub(self.input_phred_offset)`
  (`src/bam.rs:627`), so the arithmetic always happens.
- `.chain(cli.passthrough.iter())` is dead in the uBAM arm (`--passthrough` + `--output-format
  ubam` bails in `Cli::validate()` §3.4a, `src/cli.rs:630`) but is load-bearing for the `--fastqc`
  arm. Correct as written.

The scope boundary itself holds up: on FASTQ output with `-q > 0`, a mislabelled input is still
accepted, but the damage there is *loud* (every score drops 31, so essentially the whole library is
trimmed away and the report says so), whereas on uBAM output it is a silent 31-point shift baked
into `QUAL`. That is a defensible line and the CHANGELOG states it.

### The claimed limit is real

`PLAN.md` §13's "the check cannot settle a uniformly high-quality file" reproduces exactly. A
Phred+33 file of uniform `I` (Q40, ASCII 73) carries no byte below 64:

```
--phred64 --output-format ubam  → exit 0, QUAL = "********************"   (ASCII 42 = Phred 9)
          --output-format ubam  → exit 0, QUAL = "IIIIIIIIIIIIIIIIIIII"   (Phred 40)
```

So the file is accepted and written 31 points low. `docs/guide/quality.md` states this; the
`NOTE` does not — see B1.

---

## Findings

### B1 — High — the reworded `NOTE` overclaims, replacing wording that did not

**`src/main.rs:567-569`.** The diff changes the closing sentence of the startup `NOTE` from

> `With --fastqc the FastQC step refuses a file holding any quality character below ASCII 64.`

to

> `A mislabelled input is refused before anything is written when the output is uBAM, or when --fastqc is given.`

The old sentence was wrong about the *actor* (post-#536 the pre-flight refuses, not FastQC) and the
rewrite correctly fixes that. But in doing so it swapped a **byte-level, provable** claim for a
**categorical** one that the check cannot deliver. A mislabelled input is refused only when it
*proves* itself mislabelled, within the first 10,000 reads. The uniform-Q40 file in the table above
is mislabelled, is not refused, and is written 31 points low — and the `NOTE` fires on that run and
tells the user the opposite.

Two things make this worth High rather than a style note:

1. **It contradicts the same diff's own correction.** `PLAN.md` §13 records: *"The CHANGELOG entry
   was tightened from 'a mislabelled Phred+33 input is refused' to 'an input holding any quality
   character below ASCII 64 is refused', because the first overclaims."* The CHANGELOG and
   `quality.md`'s caution were fixed; the `NOTE` — the one surface printed on **every** `--phred64`
   run, including the ones the check cannot settle — kept the overclaim.
2. **It undercuts the `NOTE`'s own opening.** The first sentence is `verify that it really is`.
   Telling the user four sentences later that a wrong claim is refused is the reason they would skip
   verifying. The `NOTE` exists to push verification; this sentence pushes the other way.

**Recommended wording** (two words longer, fully accurate, reuses the phrasing already settled on
for the CHANGELOG):

```
... writing uBAM they do not. An input holding any quality character below ASCII 64 is
refused before anything is written when the output is uBAM, or when --fastqc is given.
```

Note that `tests/integration_fastqc_phred64.rs:648` asserts on `"writing uBAM they do not"`, which
sits *before* the sentence to change, so this rewording does not touch the new test.

---

### B2 — Medium — the `NOTE`'s mode list now reads as exhaustive, and omits two modes in the same class

**`src/main.rs:564-566`.** Before: `Modes that do no quality arithmetic, such as --hardtrim5/3 and
--clump_only writing FASTQ, pass the scores through untouched.` After: `Writing FASTQ, the modes
that do no quality arithmetic — --hardtrim5/3 and --clump_only — pass the scores through untouched;
writing uBAM they do not.`

Fronting `Writing FASTQ,` genuinely fixes the old ambiguity (whether `writing FASTQ` attached to
`--clump_only` alone or to both), so the brief's question — *is it now unambiguous about which modes
and which output format* — is a yes on the output-format axis. On the **mode** axis it got worse:
dropping `such as` turns a non-exhaustive illustration into an em-dash apposition that reads as a
complete enumeration, and the enumeration is incomplete. `--clock` and `--implicon` slice quality
strings and never touch the values (`src/specialty.rs` makes no `quality::` call at all), so they
are in exactly the same pass-through class — and both are FASTQ-only, since
`--output-format ubam` rejects them in `Cli::validate()` (`src/cli.rs:633,639`).

Verified against the current binary, both with `--phred64` on Phred+33 input:

| Mode | input qual (first 20 B of the retained window) | output qual |
|---|---|---|
| `--clock` | `GGGGGGGGGGGGGFGGGGGG` | `GGGGGGGGGGGGGFGGGGGG` |
| `--implicon` | `CCCCCGGGGGGGGGGGGGGG` | `CCCCCGGGGGGGGGGGGGGG` |

Both runs exit 0 and both print this `NOTE`. A user running `--clock --phred64` reads a list that
excludes their mode and reasonably concludes their scores *are* being shifted, which is the exact
inverse of the truth. Either restore the non-exhaustive marker or complete the list:

```
Writing FASTQ, the modes that do no quality arithmetic — --hardtrim5/3, --clock,
--implicon and --clump_only — pass the scores through untouched; writing uBAM they do not.
```

Secondary, and smaller: `the modes that do no quality arithmetic` is used as a restrictive
identifier and then immediately contradicted by `writing uBAM they do not`, i.e. under uBAM they
*do* arithmetic. `the modes that do no quality trimming` would identify the same set without the
self-contradiction.

---

### B3 — Low — `quality.md` carries the same overclaim, with a mitigating pointer

**`docs/src/content/docs/guide/quality.md:37`**: *"Because that subtraction is arithmetic on the
declared offset, a mislabelled input is refused before anything is written — see the caution
below."*

Same categorical claim as B1. It is Low rather than High because `— see the caution below` points
straight at the paragraph that states both residuals (the 10,000-read bound and the
nothing-to-find class), and the caution itself is accurate and newly complete. Worth the same
`an input holding any quality character below ASCII 64` substitution for consistency with the
CHANGELOG, but it is not misleading as it stands.

While there: the caution is now a single ~200-word paragraph carrying five distinct points (verify
first; the 31-point shift; the `NOTE`; the byte-level check; two separate residuals). Splitting
after `...refused up front rather than part-processed.` would make the two residuals findable.
Prose-only, no accuracy issue.

---

### B4 — Low — the CHANGELOG points at a path that does not resolve

**`CHANGELOG.md:33`**: `see \`docs/guide/quality.md\` for what the check can and cannot settle.`
The file is `docs/src/content/docs/guide/quality.md`; `docs/guide/quality.md` resolves neither on
disk nor on GitHub nor on the rendered site.

Mitigating: **line 180 of the same file already does this** (`see \`docs/guide/outputs.md\``), so
the new entry follows an existing convention rather than inventing one. Against it: every other
doc reference in the recent entries is a live site link (`https://trimgalore.com/modes/poly/#stage-order`
at line 46, `.../performance/clumpy/#memory` at line 73), and line 667 of this very file is an
entry recording a previous bug of exactly this shape — a changelog pointing at a path that existed
in the repo but not in the renderings. `https://trimgalore.com/guide/quality/#phred-encoding` is the
consistent form, and the anchor is already verified (see below).

---

### B5 — Low — the comment block above the gate is now seven lines

**`src/main.rs:381-387`.** The project's `CLAUDE.md` sets one line as the default and two as the
maximum. The diff's own two lines are compliant in isolation and are the right *kind* of comment —
current state, no narrative, no apology:

```rust
// On uBAM output the declared offset feeds the QUAL arithmetic, and bases
// below Q31 floor at 0 rather than shifting.
```

Both halves check out: `BamWriter::write_record` uses `saturating_sub` (`src/bam.rs:627`), so
ASCII < 64 floors at 0 rather than wrapping, and `(33 + Q) - 64 = Q - 31` gives the stated shift.

The issue is the block they join. The pre-existing five lines from #536 include *"With `--fastqc`
the bundled engine contradicts a wrong one, but only after the trimmed output and the reports are on
disk"* — which is precisely the "justify the change in the source" shape `CLAUDE.md` rules out; that
argument belongs in the #536 commit message. Not this diff's defect, but the diff grows the block
from 5 to 7 lines, and the natural time to cut it is while touching it. A two-line block stating
only the invariant would carry everything a future reader needs:

```rust
// `--phred64` is a claim about the input; a Phred+64 file cannot hold a quality
// character below ASCII 64, and on uBAM output that offset feeds the QUAL arithmetic.
```

---

### B6 — Low — the new test block has no acceptance case for the arm it adds

The four new tests are **all refusals**. The complementary direction — genuine Phred+64 accepted on
the uBAM arm *without* `--fastqc` — is covered only from other files
(`integration_ubam_out.rs:632,666,799` and `integration_clump_only_ubam.rs:779`, all using the
ASCII-66 fixture). That is real coverage and I confirmed all four still pass by hand, so this is not
a gap in the suite, only in the block's self-containment: someone reading
`integration_fastqc_phred64.rs` alone sees a gate that only ever refuses.

`preflight_accepts_ascii_64_which_is_phred64_q0` pins the boundary on the `--fastqc` arm; the
cheapest closure here is the same fixture with `--output-format ubam` and no `--fastqc`, asserting
exit 0 and that the `.bam` exists.

---

### B7 — Nit — the `#500` entry two bullets down is now stale

**`CHANGELOG.md:35-40`**: *"A Phred+33 file mistakenly run with `--phred64` and `--fastqc` now fails
at the FastQC step rather than producing a report."* Since #536 it fails at the input pre-flight,
not at the FastQC step. Pre-existing (introduced by #536, not this diff) — flagged only because the
unreleased 2.4.0 section now carries three adjacent `--phred64` entries and this is the one that
contradicts the other two. The section is still unreleased (`v2.4.0` is not tagged;
`### Version 2.4.0` is accumulating, as the `49152c1` commit did), so it is cheap to fix now and
not after.

---

## Checks that came back clean

- **`cargo fmt --all -- --check`** — exit 0.
- **`OutputFormat::UBam` qualification style** — all 10 occurrences in `main.rs` use the full
  `trim_galore::cli::` path; `OutputFormat` is deliberately not imported. The new code matches.
- **Test traps of the `--length` shape.** `preflight_fires_on_clump_only_ubam_output` correctly
  omits both `--length` and `-q`, and says why. I checked the other three for the same shape:
  `--output-format ubam` with `-q 0 --length 0` (normal trim), with `--hardtrim5 20`, and the
  FASTQ-output `--hardtrim5 20` case all pass `Cli::validate()` cleanly, so none can exit non-zero
  for an unrelated reason. No second instance of the trap.
- **Refusals are attributed, not just non-zero.** Three of the four new tests go through
  `assert_preflight_refused`, which requires the literal `--phred64 declares Phred+64 input` *and*
  that the output directory was never created; the fourth adds `ASCII 35`. A `--clump_only`
  flag-conflict error or a FastQC refusal would fail all of them.
- **The new tests are genuine regressions.** Without the source change all four fail: the three
  refusal tests would see exit 0, and `note_qualifies_pass_through_by_output_format` asserts on
  `"writing uBAM they do not"`, a string the old `NOTE` does not contain. `PLAN.md` §13's mutation
  table reports the same three dying to the dropped `|| matches!(… UBam)` arm.
- **`docs/guide/outputs.md` bullet.** Added to a list headed *"Rejected at startup, before anything
  is written"* while `--phred64` is only conditionally rejected — but the `--rename` bullet
  directly above it is already the same qualified shape, so the pattern is established. The link
  `/guide/quality/#phred-encoding` resolves: `## Phred encoding` is at `quality.md:29`, and the
  root-relative-with-trailing-slash form matches `/modes/clump-only/` in the neighbouring bullet.
  The bullet is ~70 words against ~15 for every other entry in the list; trimming the mechanism
  explanation (the linked page has it) would fit the list better. Cosmetic.
- **`CHANGELOG` length.** The new entry is a bold headline plus one semicolon-joined sentence, ~55
  words — shorter than the #532 entry above it and shorter than most of the section. It conforms to
  the file's existing format. No finding.
- **Efficiency.** The gate adds one bounded scan (≤ 10,000 records) per input on `--phred64` uBAM
  runs, and re-calls `detect_input_format` per input although `input_formats` was computed at line
  340. Both are pre-existing from #536 and both are trivially cheap (24-byte peek plus at most one
  BGZF block); the re-detection is also the simpler option, since the loop chains `--passthrough`
  and so cannot index `input_formats`. Nothing to change.
- **The `continue` for `InputFormat::UnalignedBam` inside the loop** is dead in the uBAM arm (uBAM
  input bailed at 369, `--passthrough` is rejected with uBAM output) but live for
  `--phred64 --fastqc --paired R1.fq R2.fq --passthrough carrier.bam`, which reaches the loop
  before the passthrough format check at line ~455. Load-bearing as its comment claims, and
  unchanged by this diff.

---

## What I did not do

- No `cargo test` and no `cargo clippy` — both build into the shared `target/`, which would have
  raced the two concurrent agents. Local rustc is 1.95.0 against CI's 1.99 in any case, so a local
  clippy pass would not have been conclusive. `PLAN.md` §13 reports 842 tests passing and clippy
  clean on a cold run; I did not independently confirm either.
- No docs-site build (`npm install` rewrites `docs/package-lock.json`). I verified the one new link
  by hand against the heading and the sibling link form instead.
