# Plan: run the `--phred64` pre-flight whenever output is uBAM

**Date:** 2026-11-06 · **Revised:** 2026-11-08 after dual plan review (see §12)
**Branch:** `dev` @ `49152c1`
**Decision:** taken 2026-11-06 — extend the #532 pre-flight to uBAM output (option 1 of 4).
**Target release:** v2.4.0 (unreleased; `master` is 3 commits behind `dev`, no tag cut).

---

## 1. Goal

Make the `--phred64` input pre-flight added in #536 fire whenever `--output-format ubam`
is in effect, not only when `--fastqc` is requested.

Today, `--output-format ubam --phred64` on Phred+33 input exits 0 and writes a BAM whose every
`QUAL` is 31 points low, with every base below Q31 irreversibly flattened to Q0. After the change
the run is refused before anything is written.

**Scope boundary.** This closes the uBAM-output hole. It does **not** close the FASTQ-output hole
described in §2.4, which is real, measured, and deliberately left open. Do not read this plan as
making `--phred64` safe.

## 2. Context

### 2.1 Placement

One condition in `src/main.rs:386`, inside `main()`. Its position is already correct and must not
move: it sits after input-format detection (line 343, which produces `any_bam`) and before
`naming::ensure_output_dir` (line 546), so a refused run creates no directory, no report and no
output file. Both reviewers verified this independently — A empirically via the sibling guard at
line 369, B structurally across all six writing dispatch branches (`--hardtrim5` 597,
`--hardtrim3` 628, `--clock` 668, `--implicon` 690, `--clump_only` ~761, `run_ubam_output` 982),
every one of which begins after both line 386 and line 546.

### 2.2 Existing code this must align with

| Site | Role |
|---|---|
| `src/main.rs:369` | `--phred64` + uBAM **input** → hard error (#362). Runs first, so uBAM inputs never reach the scan. |
| `src/main.rs:386` | **The condition being widened.** Currently `cli.phred64 && cli.fastqc_requested()`. |
| `src/main.rs:26` | `PHRED64_PREFLIGHT_READS = 10_000` — the scan bound. Unchanged. |
| `src/main.rs:36` | `lowest_quality_char(path, max_reads)` — the scan. Unchanged. |
| `src/main.rs:456` | `--passthrough` uBAM-carrier refusal. **Seventy lines downstream of the scan** — see §11. |
| `src/main.rs:553` | The startup NOTE. Misleading as written (§5.3). |
| `src/main.rs:2436, 2595, 2715` | The three uBAM-path `fastqc::run` sites. All pass `phred64: false`. |
| `src/cli.rs:605` | `fastqc_requested()` — `self.fastqc \|\| self.fastqc_args.is_some()`. |
| `src/cli.rs:13` | `OutputFormat::{Fastq, UBam}`. |
| `src/cli.rs:630` | `--passthrough` + uBAM output → rejected, so that `.chain()` arm adds nothing in the uBAM case. Verified. |
| `src/cli.rs:826` | `--clump_only` + `--length` → rejected. **Trap for the test in §9 — see C-2.** |
| `src/bam.rs:534,550,615` | `BamWriter` subtracts `cli.phred_offset()` per byte (#362). This is the arithmetic a false claim feeds. |

### 2.3 Prior art, and why this is not a reversal

- **#358 → #362.** #358 reported the *inverse* bug: genuine Phred+64 input written 31 points too
  *high*, because the writer hardcoded `-33`. The fix made the writer honour the declared offset
  and added the uBAM-input guard. #358 never examined mislabelled Phred+33 FASTQ, so there is no
  earlier "NOTE over refusal" decision being overturned. The present exposure is a direct and
  intended consequence of #362: the writer now trusts the declaration, so a false declaration
  propagates into `QUAL`.
- **#532 → #536.** Added the verification, gated on `--fastqc`.

  **Correction (review I-1).** An earlier draft justified widening the gate by saying `--fastqc`
  was "the one configuration that already had a backstop". That is wrong for uBAM output: all
  three uBAM-path `fastqc::run` call sites pass `phred64: false`, so **uBAM output has no FastQC
  backstop with or without `--fastqc`**. The backstop at `src/fastqc.rs:79` exists only on the
  FASTQ-output path. #536's own commit body had already corrected this claim. The argument for
  widening does not depend on it and is strengthened without it.

### 2.4 Measurements

Fresh `--release` build at `49152c1`, Phred+33 input, `QUAL` read back with `samtools`. The whole
table was independently reproduced by both reviewers.

**uBAM output** (`I` = Q40 input):

| Invocation | Exit | Result |
|---|---|---|
| `--output-format ubam --phred64` (default `-q 20`) | 0 | reads discarded, empty BAM — loud |
| `--output-format ubam --phred64 -q 5` | 0 | Q40 → **Q9**, records written — silent |
| `--output-format ubam --phred64 -q 0` | 0 | Q40 → **Q9** — silent |
| `--hardtrim5 20 --output-format ubam --phred64` | 0 | Q40 → **Q9** — silent |
| `--clump_only --output-format ubam --phred64` | 0 | Q40 → **Q9** — silent |
| `--hardtrim5 20 --phred64` (FASTQ out) | 0 | quality bytes verbatim — correct |

Mixed-quality fixture spanning Q0–Q40 through `--output-format ubam --phred64 -q 0`:

```
in : 0  5 10 15 20 25 30 31 32 35 38 40 …
out: 0  0  0  0  0  0  0  0  1  4  7  9 …
7 of 12 distinct quality levels destroyed
```

`saturating_sub` floors at 0, so every base below Q31 collapses irreversibly. #358's
"recoverable via `samtools fastq -n`" does **not** carry over — that applied to scores written too
high, where nothing clamped.

> **Dropped figure (review O-1).** An earlier draft said "48% of bases flattened to Q0". That
> number is a property of the fixture's composition, not of the code — reviewers measuring with
> their own fixtures got different values. Cite the level-collapse, not a percentage.

**FASTQ output is NOT inert (review C-1).** An earlier draft measured only `-q 0` and `-q 20` and
interpolated across the gap, concluding that FASTQ output "passes scores through verbatim" and was
therefore safe. The quality *bytes* are indeed written verbatim — but the trimming *decision* is
made on shifted scores, and between those two points there is a band where good bases are silently
deleted. Measured on a 60 bp fixture (Q40×40 then Q35×20 at Phred+33), FASTQ output, no
`--fastqc`:

| `-q` | written | |
|---|---|---|
| 0, 1, 2 | 59 bp | nothing trimmed |
| **5, 9** | **39 bp** | **20 good bases silently deleted, exit 0** |
| 10, 11, 20 | discarded | loud |

The mechanism: the band is open whenever the read's *shifted* qualities straddle the cutoff —
`min_shifted < q ≤ max_shifted`. Below it nothing trims; above it everything is discarded, which
is loud. For typical Illumina data (max ASCII `I`=73 or `F`=70 → max shifted Q9 or Q6) the band
requires `-q` below about 10, so the default `-q 20` is self-defending — but a lowered `-q` is a
legitimate choice for some protocols, and there the loss is silent.

**Consequences for this plan:**

1. The scope is still correct — gating on output format protects the uBAM paths, and this plan
   does not claim to fix FASTQ.
2. **§10's "trim the NOTE" recommendation reverses.** The NOTE's `quality trimming` clause is the
   only warning a user in that band receives. Trimming it would remove the sole mitigation for a
   hole this plan leaves open.
3. A wider gate was considered and rejected on evidence: `cli.quality > 0` does not discriminate,
   because `-q` defaults to 20 and is silently accepted even under `--clump_only`, which does no
   quality trimming at all.

### 2.5 Which modes are exposed

The exposed uBAM modes are the ones doing **no quality arithmetic** — `--hardtrim5/3` and
`--clump_only` — because they have no threshold to trip. Reviewer A established this holds by
construction, not coincidence: every `cli.phred_offset()` call reaching a `BamWriter` sits behind
an `OutputFormat::UBam` match arm, so the output-format gate *is* the affected set rather than a
proxy for it.

## 3. Behavior

### Prerequisites

1. `cli.validate()` has run (`src/main.rs:322`).
2. `input_formats` / `any_bam` are computed (`src/main.rs:343`).
3. The uBAM-input guard at `src/main.rs:369` has already returned for any uBAM input.
4. `ensure_output_dir` has **not** yet run.

### Main logic

5. Enter the scan when `cli.phred64` **and** (`cli.fastqc_requested()` **or**
   `cli.output_format == OutputFormat::UBam`).
6. For each path in `cli.input` chained with `cli.passthrough`:
   a. If `detect_input_format(path)` is `UnalignedBam`, skip. **Load-bearing — see §11.**
   b. `lowest_quality_char(path, PHRED64_PREFLIGHT_READS)`; on `Some(c)` where `c < 64`, bail
      with the existing message, including the Solexa clause for ASCII 59–63.
7. Otherwise fall through unchanged.

### Edge cases

| Case | Handling |
|---|---|
| uBAM input + `--phred64` | Already refused at line 369; unreachable here. |
| Genuine Phred+64 FASTQ → uBAM | **Must keep working.** Lowest byte ≥ 64, so no bail. Four existing tests guard this — §9. |
| Empty file / no quality bytes | `lowest_quality_char` → `None` → no refusal. Unchanged; the sanity check owns that error. |
| `--hardtrim5/3` + uBAM, mislabelled | **Now refused** (was silent corruption). |
| `--clump_only` + uBAM, mislabelled | **Now refused** (was silent corruption). |
| Any mode, FASTQ output | **Unchanged** — still accepted, including the §2.4 band. NOTE only. |
| `--passthrough` + uBAM | Rejected at `cli.rs:630` before this point. |
| FIFO / process-substitution input | Already rejected outright, so the second read pass adds no new constraint (reviewer A). |
| Low character past read 10,000 | Still not caught. Near-moot on uBAM (§8). |

## 4. Signature

No signature changes. `lowest_quality_char` and `PHRED64_PREFLIGHT_READS` are reused as-is.

## 5. Implementation outline

### 5.1 Widen the condition — `src/main.rs:386`

```rust
if cli.phred64
    && (cli.fastqc_requested()
        || matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam))
{
```

### 5.2 Update the comment above it — `src/main.rs:381-385`

**At most two lines, stating the current rule, not justifying the change** (CLAUDE.md). An earlier
draft's example ran to three lines and argued the case in the source; both are breaches. Something
of this shape:

```rust
// On uBAM output the declared offset feeds the QUAL arithmetic, and bases
// below Q31 floor at 0 rather than shifting.
```

Do not mention the previous narrower gate, FASTQ behaviour, or the reasoning.

### 5.3 Disambiguate the NOTE — `src/main.rs:553-564`

The NOTE currently reads:

> Modes that do no quality arithmetic, such as `--hardtrim5/3` and `--clump_only` writing FASTQ,
> pass the scores through untouched.

Whether "writing FASTQ" distributes over both list items or only the second is undecidable from
the grammar. Read as applying to both, the sentence is true; read as applying only to
`--clump_only`, it is false, and reviewer A observed it printing on a run that then wrote Phred 9.

**The fix is to disambiguate**, so the qualifier unmistakably governs both modes. An earlier draft
called this a "factual error" and ranked it a standalone correctness defect; that overstates it
(review I-5) and the ranking is withdrawn. It remains worth fixing, and the fix is independent of
the rest of this plan.

Also update the trailing `--fastqc` sentence to state the widened gate. **Do not shorten the
`quality trimming` clause** — per §2.4 it is the only warning covering the FASTQ band.

### 5.4 Refresh the stale test comment — `tests/integration_fastqc_phred64.rs:304`

`preflight_does_not_fire_without_fastqc` passes FASTQ output and asserts
`p33_nofastqc_trimmed.fq` exists, so it **stays green**. Its comment ("The check is deliberately
scoped to `--fastqc`") becomes wrong. Update it, and rename the test to
`preflight_does_not_fire_without_fastqc_on_fastq_output`.

### 5.5 Add tests

Per §9, in `tests/integration_fastqc_phred64.rs`, reusing its existing helpers:
`write_phred33_fixture` (243), `run_raw` (255), `fixture_copy` (111), and — **required** —
`assert_preflight_refused` (536), which asserts both that stderr carries
`"--phred64 declares Phred+64 input"` *and* that the output directory was never created. Without
it a test cannot distinguish the pre-flight's refusal from any other non-zero exit. This is what
makes C-2 below a real trap rather than a style point.

### 5.6 CHANGELOG

Add one entry under `### Version 2.4.0`. An earlier draft said "one sentence, matching the file's
format", which is unsatisfiable as written (review O): the file's entries are bold-lede paragraphs
and the #532 entry runs seven lines. Resolve in favour of CLAUDE.md: **bold lede naming the change
plus one sentence of substance**, pointing at the docs rather than inlining migration guidance.

Then **amend the existing #532 entry**: its closing sentence, "`--phred64` without `--fastqc` is
unaffected", becomes false for uBAM output.

### 5.7 Docs

Two files, not one.

- **`docs/src/content/docs/guide/outputs.md`** — the canonical "Rejected at startup, before
  anything is written (v1 scope)" list for uBAM output (reviewer A). `--rename` is already there
  as a *conditional* entry ("only when at least one input is FASTQ, and only when…"), which is the
  precedent for how to phrase a refusal that fires only on mislabelled input. Add an entry in that
  shape.
- **`docs/src/content/docs/guide/quality.md`** — line 40's `:::caution` scopes the check to
  `--fastqc`; state that it also fires whenever the output is uBAM. Line 37 says `--phred64` is
  honoured for FASTQ input under `--output-format ubam` — still true for genuine Phred+64; add
  that a mislabelled input is refused. **Keep the caution's existing warning about quality
  trimming**, which per §2.4 still describes live behaviour on the FASTQ path.

Plain scientific register; no AI-slop phrasing.

## 6. Efficiency

One extra `detect_input_format` + bounded scan per input, only on runs that both pass `--phred64`
and write uBAM. The scan reads at most 10,000 records through `open_sync_reader` (gzip-aware,
64 KB buffered) and allocates nothing beyond one record at a time. Negligible against a trimming
run, and it happens before any output work, so a refused run does strictly less I/O than today.

The scan is a second pass over the first 10,000 reads on `--fastqc` runs that also write uBAM
(adapter auto-detection is a separate pass). Not worth merging; the bound makes it cheap.

## 7. Integration

**Reads:** `cli.phred64`, `cli.output_format`, `cli.fastqc_requested()`, `cli.input`,
`cli.passthrough`; the first 10,000 records of each non-BAM input.
**Writes:** nothing. On failure it returns `Err` before `ensure_output_dir`.

**Downstream:** the invariant block at `src/main.rs:415-424` concerns uBAM *input* and is
untouched. This change adds a parallel, weaker guarantee for uBAM *output*: when the declared
offset is 64, the input has been checked to contain no byte below ASCII 64 within the scan bound.

**Behaviour change:** new refusals in previously-succeeding runs (`--hardtrim5/3`, `--clump_only`,
and plain trim with low `-q`, all with uBAM output and mislabelled input). That fits the v2.4.0
minor bump, which already carries several new refusals.

## 8. Assumptions

**Fixed rules** (not configurable, derived from the encoding definition):
- Illumina 1.3+ Phred+64 starts at ASCII 64 (Q0 *is* ASCII 64), so genuine Phred+64 data cannot
  contain a byte below 64. This is what makes the check sound rather than heuristic.
- BAM stores raw Phred by spec, so BAM inputs are correctly skipped.

**Verified in-tree, not assumed:**
- `--passthrough` is rejected with uBAM output at `src/cli.rs:630`.
- `--clump_only` rejects `--length` at `src/cli.rs:826`.
- All three uBAM-path `fastqc::run` sites pass `phred64: false`.
- `phred64_test.fastq`'s lowest byte is `'B'` (66) — 24×Q40 (`'h'`) then 10×Q2 (`'B'`).

**Assumptions taken:**
- The 10,000-read bound stays. On uBAM the residual is near-moot: real Phred+33 data carries
  sub-ASCII-64 bytes in essentially every record. The bound remains a genuine weakness on the
  FASTQ/`--fastqc` path, where #536's own residual note applies unchanged.
- `--hardtrim5/3` and `--clump_only` remain permitted with uBAM output (not added to §3.4a). This
  plan verifies the claim rather than prohibiting the combination, because #362 deliberately made
  genuine Phred+64 → uBAM correct and a blanket refusal would take that back.

## 9. Validation

Run from the crate root.

| # | Verify | How | Expected |
|---|---|---|---|
| 1 | Mislabelled input refused before any write | `--output-format ubam --phred64 -q 0` on a Phred+33 fixture | `assert_preflight_refused` passes: stderr names the pre-flight, output dir never created |
| 2 | The silent hardtrim path is closed | `--hardtrim5 20 --output-format ubam --phred64` | same |
| 3 | The silent clump_only path is closed | `--clump_only --output-format ubam --phred64` — **no `--length`, no `-q`** | same |
| 4 | The NOTE no longer claims hardtrim is safe | assert on stderr of test 2 before the refusal, or on a FASTQ-output hardtrim run | stderr does not assert unqualified pass-through for `--hardtrim5/3` |
| 5 | Error names the offending file and character | assert on stderr | message contains the path, the character, its ASCII value |
| 6 | **#362 regression — existing guards stay green** | `phred64_ubam_out_se_stores_true_phred` (`integration_ubam_out.rs:632`), `…_specialty_…` (666), `…_pe_two_file_…` (799), `phred64_clump_only_ubam_out_stores_true_phred` (`integration_clump_only_ubam.rs:779`) | all pass unchanged |
| 7 | FASTQ output unaffected | `preflight_does_not_fire_without_fastqc_on_fastq_output` | stays green |

**C-2 — the trap in row 3.** Writing row 3 by copying the existing test pattern inherits
`--length 0`, which `--clump_only` rejects at `src/cli.rs:826`:

```
--clump_only --output-format ubam --phred64 -q 0 --length 0
  → exit 1, "--clump_only does not filter; --length is not compatible", 0 files written
```

Non-zero exit and nothing written is exactly what a naive assertion checks, so the test would pass
**green today, before any implementation**, and prove nothing. `assert_preflight_refused` is what
closes this, because it asserts on the refusal's *text*. Only row 3 carries the trap; rows 1 and 2
were checked and are clean.

**Row 6 replaces a redundant new test.** An earlier draft proposed writing a fresh genuine-Phred+64
→ uBAM guard. Four already exist (review I-3), in files this plan had not surveyed, and
`integration_ubam_out.rs:616` already provides a `quals()` helper for reading BAM `QUAL` — the
proposed new test would have duplicated `…:632` inside a file with no BAM reader at all. Note that
`…_specialty_…` at 666 exercises `--hardtrim5 20 --phred64 --output-format ubam` on the genuine
fixture, which is precisely the invocation this change affects. Baseline: 75 tests green across the
three files.

**Mutation check.** For each new test, confirm it *fails* under a deliberately broken condition:

| Mutation | Must fail |
|---|---|
| drop the `\|\| matches!(… UBam)` arm | rows 1–3 |
| change `\|\|` to `&&` | rows 1–3 |
| `lowest < 64` → `lowest <= 64` | **`preflight_accepts_ascii_64_which_is_phred64_q0`** (`tests/integration_fastqc_phred64.rs:508`) |

The third pairing is corrected from an earlier draft, which assigned it to the genuine-Phred+64
row. Both reviewers caught this independently and it was confirmed in the main session: that
fixture's minimum byte is 66, and `66 <= 64` is false, so the mutant does not bail and the test
passes either way. The mutation was right, the target was wrong — exactly the failure mode the
handoff warns about, committed while citing it. The ASCII-64 boundary test already covers it.

Dry-run the mutator and check its exit status before trusting the loop, and **rebuild the release
binary after restoring** — a leftover mutant binary has produced a false "does not reproduce" here.

**Suite:** `cargo test` (crate root), `cargo fmt --all -- --check`,
`cargo clippy --all-targets --release -- -D warnings` **cold**. Prefix cargo with
`DEVELOPER_DIR=/Library/Developer/CommandLineTools`. Local rustc is 1.95.0 against CI's 1.99, so
clippy-clean locally is not clippy-clean in CI.

**Not run by either reviewer:** no source mutation, no full suite, no `fmt`, no cold clippy. The
corrected mutation pairing above is arithmetic, not a measured result.

**Not covered:** the Perl byte-identity validation matrix is untouched — uBAM output has no Perl
counterpart, and no FASTQ path changes.

## 10. Questions or ambiguities

- **Resolved — NOTE length.** An earlier draft recommended trimming the NOTE. **Reversed** on
  C-1: the `quality trimming` clause is the only warning covering the FASTQ band this plan leaves
  open. Keep it; only disambiguate the mode list and update the gate sentence.
- **Open — the FASTQ band (§2.4).** Out of scope here, and no good gate is available
  (`cli.quality > 0` does not discriminate). Worth its own issue: a mislabelled Phred+33 input with
  `-q` in the band silently loses good bases on the FASTQ path, which this plan does not refuse.
  Recommend filing rather than widening scope — it meets the data-integrity bar.
- **Open — `-q 0` on FASTQ output.** Passes through with scores intact but a wrong
  `Quality encoding type selected: ASCII+64` line in the report. Harmless to the data.

## 11. Self-Review

**Efficiency.** One bounded scan on a narrow condition; no new allocation; refused runs do less I/O
than today.

**Logic.** Ordering walked: detection (343) → uBAM-input guard (369) → this scan (386) → invariant
block (415) → `ensure_output_dir` (546) → NOTEs (553) → dispatch (597+).

**The in-loop BAM skip is load-bearing.** An earlier draft called it "belt-and-braces for the
`--passthrough` arm rather than load-bearing" — inverted, and both reviewers caught it. The
`--passthrough` uBAM-carrier refusal lives at `src/main.rs:456`, seventy lines *downstream* of the
scan. The skip is the only thing that lets a uBAM carrier reach that correct message; without it
the scan bails first, naming a quality character the BAM does not contain. The plan retains the
skip, so no defect ships, but §11 is where a future refactor will look — hence this correction
rather than a silent fix.

**Remaining risks.**
- The 10,000-read bound is unchanged and still unsound in principle. Accepted; see §8.
- The FASTQ band stays open. Now stated rather than claimed away.
- `--hardtrim5/3` and `--clump_only` stay off §3.4a's exclusion list, so any future uBAM path
  reading `phred_offset()` without passing through `main()`'s scan would reopen this. The invariant
  comment at `src/main.rs:415` is the right place to record the output-side guarantee.
- CI cannot reproduce 1.99-only lints locally.

## 12. Revision History

**2026-11-08 — revised after dual plan review** (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`).
Both reviewers endorsed the code change and its siting; every finding below concerns the plan's
justification or validation. Changes folded in:

| From | Change |
|---|---|
| B C-1 | §2.4 rewritten: FASTQ output silently deletes good bases in a measured `-q` band. Premise corrected, band characterised by sweep, §10's NOTE-trimming recommendation **reversed**, new issue recommended. |
| B C-2 | §9 row 3 stripped of `--length`; §5.5 now requires `assert_preflight_refused` (536). |
| A 1 / B I-2 | §9 mutation pairing corrected to `preflight_accepts_ascii_64_which_is_phred64_q0` (508). |
| B I-1 | §2.3 backstop claim corrected — all three uBAM `fastqc::run` sites pass `phred64: false`. |
| B I-3 | §9 row 6 replaces a redundant new test with four existing guards; `quals()` helper noted. |
| A 3 / B I-4 | §11 BAM-skip inversion corrected. |
| B I-5 | §5.3 "factual error" softened to ambiguity; standalone-defect ranking withdrawn. |
| A 4 | §5.7 adds `docs/guide/outputs.md`, with `--rename` as the conditional-entry precedent. |
| A opt. | §5.2 comment example cut to two lines, justification removed. |
| A opt. | §5.6 "one sentence matching the file's format" resolved as unsatisfiable. |
| A / B | §2.4 "48%" figure dropped as fixture-dependent. |
| A | §3 edge-case table gains the FIFO row. |

Independently re-verified in the main session before folding: C-1 (39 bp vs 59 bp, plus the `-q`
sweep), C-2 (exit 1 on `--length`), I-1 (three `false` call sites), the ASCII-64 boundary test, the
four existing uBAM guards, `assert_preflight_refused`, and the `outputs.md` list.

## 13. Implementation notes — 2026-11-08

All seven steps of §5 implemented on `dev` (uncommitted). Verification: **842 tests pass**
(`cargo test --release`, exit 0), `cargo fmt --all -- --check` clean, `cargo clippy --all-targets
--release -- -D warnings` clean on a cold run (29 crates checked from scratch; a forced re-check of
the crate alone also exits 0).

### Mutation results — all three killed their predicted targets

Source backed up and restored by checksum (`md5 5d3acffbfc4a6a89f8bd0953063f85d7`) between each
mutation; every mutation was asserted to have applied before its run, and the binary was rebuilt
after the final restore.

| Mutation | Predicted | Observed |
|---|---|---|
| drop the `\|\| matches!(… UBam)` arm | rows 1–3 fail | exactly those 3 failed, 15 passed |
| `\|\|` → `&&` | rows 1–3 fail | those 3 **plus 5 pre-existing `--fastqc` tests** failed |
| `lowest < 64` → `<= 64` | `preflight_accepts_ascii_64_which_is_phred64_q0` fails | exactly that 1 failed, 17 passed |

The second mutation killing five existing tests confirms the suite also guards the `--fastqc` arm,
which the plan had not claimed. The third confirms the §9 correction: the original draft's target
would indeed have survived.

### Deviation 1 — §9 row 4 uses the FASTQ-output variant

The NOTE is printed at `src/main.rs:558`, *after* the pre-flight at 386. A refused run therefore
emits no NOTE at all, so "assert on stderr of test 2 before the refusal" is not implementable. §9
offered a FASTQ-output hardtrim run as the alternative and that is what
`note_qualifies_pass_through_by_output_format` does. Within plan intent; recorded because the
ordering is not obvious from the plan text.

### Finding — the check cannot settle a uniformly high-quality file

**This is a limitation of the approach, discovered during end-to-end verification, not a bug.**

A Phred+33 file whose every base is Q31 or better carries no byte below ASCII 64. Every such byte
is legal Phred+64, so no content-based check can refuse the file without also rejecting genuine
Phred+64 input — which is exactly what #362 made work. The check is therefore maximally sound: it
refuses precisely those inputs that *prove* they are not Phred+64.

The sharp consequence: **§2.4's headline fixture is in that class.** The uniform `I` (Q40) file
used to demonstrate Q40 → Q9 is *not* refused by this change, and on uBAM output it still writes
every `QUAL` 31 points low. Verified after implementation:

| Input | uBAM + `--phred64` |
|---|---|
| uniform `IIII…` (all Q40, no byte < 64) | exit 0, `QUAL` = Q9 — **unchanged** |
| `IIII…####IIII…` (contains ASCII 35) | exit 1, output dir never created — refused |

Real Phred+33 data essentially always falls in the second class: Q2 `B`-run tails and binned low
scores (NovaSeq Q12/Q23) are all below ASCII 64. Synthetic, heavily filtered, or uniformly
high-quality files can reach the first.

Two documentation changes followed from this rather than from the plan:
- The CHANGELOG entry was tightened from "a mislabelled Phred+33 input is refused" to "an input
  holding any quality character below ASCII 64 is refused", because the first overclaims.
- `guide/quality.md`'s caution gained a sentence stating what the check cannot settle. This is a
  *second*, distinct residual from the 10,000-read bound already recorded there: the bound is about
  where the scan stops looking, this is about files that contain nothing to find.

`guide/outputs.md`'s entry was already phrased conditionally ("a file that does is refused") and
needed no change.

### Not done

- **The docs site was not built.** `npm install` rewrites `docs/package-lock.json` on every branch
  switch, and the changes are prose plus one link in an existing list. The link target
  (`/guide/quality/#phred-encoding`) was verified against `## Phred encoding` at `quality.md:29`
  and against the `/guide/quality/` form used in `guide/overview.md:18` and `guide/flags.md:49`.
- **No test pins the uniform-high-quality residual.** `preflight_accepts_genuine_phred64_input`
  already pins the soundness direction (lowest ≥ 64 is accepted); a test asserting that a
  mislabelled uniform file is accepted and written shifted would encode the limitation as if it
  were desired behaviour.

## 14. Post-review fix round — 2026-11-08

Dual code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) and the coverage audit
(`COVERAGE.md`, **verdict COMPLETE**) all found the functional change correct: no defect at any
severity from either reviewer. Every finding was in user-facing prose. Fix scope chosen:
correctness-of-claim only.

### Correction to §13 — the "finding" was already documented

§13 presents the uniform-high-quality blind spot as discovered during implementation. It was not.
**`test_files/README.md:23-28` already records it**, in those words:

> The heuristic's actual weakness runs the other way, and `phred64_test.fastq` demonstrates it:
> uniformly high-quality data leaves no byte below 64, so the minimum tells you nothing and the
> guess falls back to Phred+64 even when the data is Phred+33.

The measurement in §13 stands; the framing does not. It was a known, written-down property of this
repository that the plan failed to look up — the same root cause as the §9 and C-2 errors, which
also came from not reading test-adjacent material.

### Fixes applied

| # | Finding | Change |
|---|---|---|
| 1 | A F-1 / B B1 — NOTE overclaimed | `src/main.rs` closing sentence restored to the ASCII-64 condition and extended to name the residual: "An input holding any quality character below ASCII 64 is refused … uniformly high-quality Phred+33 data holds no such character and is accepted." |
| 2 | B B2 — mode list read as exhaustive | "such as" restored, so `--clock` and `--implicon` are no longer implicitly excluded. Verified: `--clock --phred64` prints this NOTE, and `specialty.rs` contains no `quality::` or `saturating_sub` call. |
| 3 | B B3 — same overclaim at `quality.md:37` | Replaced with the ASCII-64 condition, pointing at the caution. |
| 4 | A F-2 — `outputs.md` misfiled the entry | Removed from the "Rejected at startup … (v1 scope)" list and rewritten as a following paragraph stating `--phred64` *is* supported and only invalid input is refused. A is right that `--rename` was not a precedent: it is refused for a broad class of legitimate runs, this is not. |

### A F-3 — refuted, reduced to one clause

A reported that the fixture's minimum-byte-66 invariant is recorded nowhere. It is:
`test_files/README.md:14` states it and names `tests/integration_fastqc_phred64.rs` as the
consumer. The genuine residual is narrower — that line attributed the dependency to one test file,
and this change adds two more (`integration_ubam_out.rs`, `integration_clump_only_ubam.rs` now hit
the same pre-flight via uBAM output). The README clause was extended to say so and to state the
failure mode a lower tail would produce. No new invariant record was needed.

### Not fixed (out of chosen scope)

B B4 (`docs/guide/quality.md` path in the CHANGELOG does not resolve — but `CHANGELOG.md:180`
already does this, so it follows convention), B B5 (the comment block is 7 lines; the 2 added are
compliant, the pre-existing 5 carry a justify-in-source shape), B B6 (no acceptance case for the
new arm inside this test file; covered from `integration_ubam_out.rs`), B B7 (`CHANGELOG.md:35-40`
still says the refusal happens "at the FastQC step" — stale since #536, predates this change), and
A's six Lows (stale module doc, a redundant `detect_input_format` that `input_formats` already
holds, `-q 0 --length 0` carried into one test, a brittle five-word prose assertion, two CHANGELOG
entries for one `if`).

### Re-verification after the fix round

`cargo test --release` **842 passed, exit 0**; `cargo fmt --all -- --check` exit 0; clippy
re-checked `trim-galore` clean (exit 0).

**Clippy caveat, stated because the first phrasing of it was wrong.** The run was *not* cold. A
scratch `CARGO_TARGET_DIR` was used, borrowing the auditor's method, but **subagents share the
parent session's scratchpad directory**, so that path already held 447 dependency artifacts from
the auditor's 10:11 run and only `trim-galore` was re-checked (1.46 s). That is adequate coverage
here — the only Rust change in this round is one string literal, and the dependency graph is
untouched since the auditor's genuinely cold 148-crate run — but it is a warm run and should not be
cited as a cold one. Local clippy is inconclusive for CI anyway: local rustc is 1.95.0 against CI's
1.99.
