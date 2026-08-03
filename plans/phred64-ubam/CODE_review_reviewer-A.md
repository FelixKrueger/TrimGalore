# CODE REVIEW — Reviewer A — issue #358 `--phred64` on the uBAM paths

**Target:** uncommitted working tree on `fix/phred64-ubam`, base `dev` @ `2e3da11`
**Scope:** 9 modified files + new `test_files/phred64_test.fastq`
**Plan:** `plans/phred64-ubam/PLAN.md` (v2)
**Reviewer:** A (independent; Reviewer B report not consulted)

## Verdict

**APPROVE WITH REVISIONS**

The two code fixes are correct, the wiring is complete, and the ordering property the
plan identified as load-bearing genuinely holds — I verified the guard fires ahead of
`ensure_output_dir` and ahead of the early-returning specialty/clump-only dispatch by
execution, not by reading. Every production `BamWriter::create` site threads a dynamic
offset; every test site passes `33` deliberately. Both false comment claims at the Bug 1
site are gone and their replacements check out against `fastq.rs`.

The revisions are concentrated in **user-facing prose**, and one of them is a factual
claim in the release notes that I empirically disproved (H1). A second wrong claim
appears in three places at once, including a `NOTE:` string users will read (M1). Both
were introduced by carrying the plan's §7 predictions into shipping text without picking
up §13's own corrections. Nothing here blocks the approach; it blocks the wording.

Two test defects were unambiguous and low-risk, so I fixed them in place (F1, F2). Full
suite re-verified after the edits: **418 passed, 0 failed** (`cargo test --release`),
`fmt --check` clean, `clippy --all-targets --release -- -D warnings` clean.

---

## What I verified by execution

The release binary at `target/release/trim_galore` already carries the fix (built from
`2e3da11` + working tree), so post-fix behaviour is directly observable. Pre-fix
behaviour I reproduced by exploiting a coincidence: **running the Phred+64 fixture
*without* `--phred64` produces a BAM byte-identical to the pre-fix artefact** (writer
subtracts 33 from `'h'`/`'B'` → raw 71/33, exactly what the bug produced). That gave me a
pre-fix specimen without rebuilding `dev`.

| Claim | Method | Result |
|---|---|---|
| Guard precedes `ensure_output_dir` | `--phred64 -q 20 ubam_test.bam -o <new dir>` | exit 1, **no directory created** ✅ |
| Guard precedes early-returning specialty dispatch | `--hardtrim5 20 --phred64 ubam_test.bam` | exit 1, `--phred64` message ✅ |
| Bug 1, trim SE | `-q 0 --length 0`, `samtools view` col 11 | `IIII…####` → raw `[40×24, 2×10]` ✅ |
| Bug 1, `--hardtrim5 20` | same | raw 40 × 20 ✅ |
| Bug 1, `--clump_only` | same | raw `[40×24, 2×10]` ✅ |
| O-3 `NOTE:` reaches the hardtrim path | stderr grep | fires before dispatch ✅ |
| No `BamReader` route bypasses `any_bam` | traced every `BamReader::*` / `open_*_reader` call | all descend from `cli.input`; `--passthrough` is `FastqReader`-only ✅ |
| Golden fixtures unmoved | `git status test_files/` | only `README.md` + new fixture ✅ |
| No CI/test passes `--phred64` | `grep -rn phred64 .github/workflows/ scripts/` | zero hits ✅ |
| `truncated.fq.gz` really is Phred+64 | `od` over quality lines: min 66, max 99, `HWUSI-EAS611` names | ✅ classification correct |

### The fixture deviation is justified — mechanism independently confirmed

I reproduced the `fastqc-rust` min-byte heuristic against the real dependency rather than
taking §13's word for it:

| BAM under test | `Encoding` reported | Base-1 mean |
|---|---|---|
| Pre-fix byte pattern (raw 71/33) | `Illumina 1.5` | **40.0** |
| Post-fix, mixed fixture (raw 40/2) | `Sanger / Illumina 1.9` | 40.0 (base 30: 2.0) |

The heuristic is the classic FastQC one (lowest observed quality char < 64 → Sanger, else
offset 64). A uniform-Q40 BAM reads back as ASCII 73 throughout, min 73 ≥ 64 → guessed
offset 64 → Q9, exactly as §13 recorded. The mixed fixture's `'B'`-run drops the minimum
to ASCII 35 and restores detection. **The deviation from the plan's uniform-`'h'`
fixture is correct and the new fixture is strictly better** — it pins a mapping rather
than a single constant.

The `-q 0` flag is load-bearing as claimed: without it the default `-q 20` trims the
`B`-run and only one quality value reaches the writer (verified — 24 bp out, all Q40).
`--length 0` is **not** load-bearing (see L4).

---

## Issues

### H1 — CHANGELOG asserts a FastQC symptom that does not occur

`CHANGELOG.md:71-72`:

> `--fastqc` on uBAM output was affected too — per-base-quality plots were centred on Q71.

**This is false, and the implementation's own notes say so.** PLAN §13 line 396 records:
"pre-fix it read ≈40 *by accident* (two errors cancelling, exactly as Reviewer A
described in #358)". The CHANGELOG carried forward the plan's §7 prediction, which §13
had already withdrawn.

Empirically, on a BAM holding the pre-fix byte pattern, FastQC reports **mean 40.0** with
`Encoding  Illumina 1.5`. The writer's 31-byte inflation and fastqc-rust's misdetected
offset cancel exactly: read adds 33, heuristic subtracts 64, net −31.

This matters beyond pedantry — it is the release-note sentence a user consults to decide
whether their archived reports are wrong. As written it tells them to distrust reports
that were numerically fine, and says nothing about the label that *was* wrong.

**Replace with something like:**

> `--fastqc` on uBAM output was mislabelled rather than misplotted: the inflated `QUAL`
> pushed every byte above ASCII 64, so fastqc-rust's encoding heuristic guessed
> `Illumina 1.5` and subtracted 64 — cancelling the writer's 31-byte inflation, so the
> plotted means happened to be right while the reported `Encoding` was wrong. Consumers
> that read BAM `QUAL` as raw Phred (samtools, aligners, Bismark) saw the inflated
> values directly.

### M1 — "all-zero-quality BAM" is wrong, in three places

The same overstatement appears in a `NOTE:` string, the docs, and the CHANGELOG:

- `src/main.rs:301-302` — "running Phred+33 data with `--phred64` yields an all-zero-quality BAM"
- `docs/src/content/docs/guide/quality.md:40` — "every quality score is zero, because the subtraction underflows and floors"
- `CHANGELOG.md:106-107` — "would silently write an all-zero-quality BAM"

The arithmetic is `q + 33 − 64 = q − 31`, floored at 0. So scores are **reduced by 31,
not zeroed**: everything at or below Q31 floors to 0, but Q32–Q41 survive as Q1–Q10.
Verified on `illumina_10K.fastq.gz` (`BBBBBFFFF…GGGG…HHHH`):

```
output QUAL: #####''''''''(((('((((())))))))…   →  raw 2, 6, 7, 8
```

Only a fixture whose maximum score is ≤ Q31 gives literal all-zero (which is why
`BS-seq_10K_R1.fastq.gz`, uniform `'?'` = Q30, looked all-zero when I first checked).

This weakens the diagnostic in the exact case it exists for. A user who reads "every
quality score is zero", then opens a BAM showing Q0–Q10, concludes the warning does not
describe their situation. It also *understates* the hazard: a Q0–Q10 BAM reads as
plausible bad-quality data, whereas an all-Q0 BAM at least looks degenerate.

**Suggested wording** (`NOTE:` and docs): "…collapses every score by 31 and floors at
zero — typical Illumina Q2–Q41 becomes Q0–Q10, which is indistinguishable from
genuinely poor data."

### M2 — the `saturating_add` residual is not benign; the comment says it is

`src/bam.rs:957-962` states the residual (a mixed-sentinel byte saturating to 255 instead
of mapping to `QUAL_MISSING_REPLACEMENT`) is out of remit and that "the point here is
only to remove the overflow". The overflow is removed, but what replaces it is not
inert. `b.saturating_add(33) as char` for `b = 0xFF` yields `U+00FF`, which is **two
UTF-8 bytes inside a `String` whose byte length is assumed equal to `seq.len()`**:

```
qual_raw = [40, 40, 255, 40]
→ chars=4 bytes=5, as_bytes: [73, 73, 195, 191, 73]
→ qual[3..]          panics: "byte index 3 is not a char boundary; it is inside 'ÿ'"
→ qual.truncate(3)   panics: "assertion failed: self.is_char_boundary(new_len)"
```

Those two operations are `FastqRecord::clip_5prime` (`src/fastq.rs:88`) and
`FastqRecord::truncate` (`src/fastq.rs:78`) — i.e. `--clip_R1`, `--hardtrim5`, and every
adapter/quality trim. So for raw bytes ≥ 223 the change trades a debug-profile
arithmetic panic for a **both-profile char-boundary panic further downstream**; in
release specifically, the old code wrapped to `' '` (1 byte, length-safe, and coincidentally
round-tripped to Q0). If no clipping occurs, the alternative outcome is a record whose
qual is longer than its seq — invalid FASTQ, and on the write side `saturating_sub` turns
`0xC3 0xBF` into two garbage scores.

Two honest framings, either acceptable:

1. **Close it.** One line, consistent with the module's existing bail-on-malformed style
   (it already bails on `'='`, on unrecognised bases, on qual/seq length mismatch):
   ```rust
   if let Some(&bad) = qual_raw.iter().find(|&&b| b > 93 && b != QUAL_MISSING_SENTINEL) {
       bail!("BAM record QUAL contains out-of-spec raw score {bad} (SAM permits 0–93)");
   }
   ```
   plus per-byte `0xFF → QUAL_MISSING_REPLACEMENT`, which is not a "semantic change" so
   much as making the mixed case agree with the all-`0xFF` case two lines above.
2. **Keep the code, fix the comment.** State that the residual can emit a non-ASCII
   quality byte and that `FastqRecord::truncate` / `clip_5prime` will panic on it, and
   file a follow-up.

What is not acceptable is the current comment, which reads as "handled, remainder
harmless". Note the multi-byte class already existed for raw 95–222 before this change —
that is a pre-existing hole, and worth saying so, but this change widens it to 223–255.

### M3 — `test_files/README.md` contradicts itself on the min-byte heuristic

`test_files/README.md:17-20`:

> `truncated.fq.gz` in particular is easy to misclassify because a naive
> minimum-quality-byte heuristic reads it as high-quality Phred+33 rather than
> moderate-quality Phred+64.

Its minimum quality byte is 66 (verified). Under the min-byte heuristic — the very one
the table invokes three lines earlier to explain the `phred64_test.fastq` deviation, and
which I confirmed fastqc-rust implements (it labelled an ASCII-min-66 BAM `Illumina 1.5`)
— 66 ≥ 64 classifies it as **Phred+64, correctly**. The heuristic gets this fixture
right.

What actually misclassifies it is *applying the default offset without checking*, which
is not a heuristic at all. The paragraph uses one mechanism to explain two fixtures and
is wrong about one of them — in a document whose entire purpose is to stop the next
person misclassifying these files.

**Suggested:** "…easy to misclassify because its scores are plausible under either
offset if you simply assume the default: min byte 66 reads as Q33 at offset 33 (looks
like good data) and Q2 at offset 64 (the truth). Note that a minimum-byte heuristic
*does* classify this one correctly; what fails is assuming Phred+33 without looking."

### M4 — `guide/quality.md` line 35 is self-contradictory and unscoped

> …before it was rejected it caused every read to be discarded as low-quality. Drop the
> flag for BAM input; results are unchanged.

Those two clauses cannot both hold. "Results are unchanged" is true only for the modes
that do no quality arithmetic (`--hardtrim5/3`, `--clump_only`) — which is exactly how
the CHANGELOG scopes the identical sentence at line 103, correctly. The docs version
drops the scope, so it simultaneously tells the reader their reads were all discarded and
that dropping the flag changes nothing.

Also "every read" holds only at `-q > 0`.

**Suggested:** "…with quality trimming it discarded every read as low-quality; in
`--hardtrim5/3` and `--clump_only` it was inert. Drop the flag for BAM input — those
modes' output is unchanged, and a trimming run now produces the result it should have
produced all along."

---

## Fixes applied

### F1 — `phred64_mixed_input_rejected` did not discriminate (`tests/integration_ubam_out.rs`)

As written it ran `--paired --phred64 phred64_test.fastq ubam_test.bam` and asserted only
`!status.success()`, with no stderr check. That input is **independently rejected** by the
pre-existing two-BAM check in the paired path (`src/main.rs:1256`):

```
$ trim_galore --paired test_files/phred64_test.fastq test_files/ubam_test.bam   # no --phred64
exit 1: "--paired with two BAM files is not supported. …"
```

So the assertion passed pre-fix, passes post-fix, and would pass with the `cli.phred64 &&
any_bam` clause deleted. Plan validation row 7 claimed "fails pre-fix"; it does not.

Rewritten to the **single-end multi-input** shape, which is otherwise a legal run
(verified: exits 0 without the flag) and is therefore attributable to the guard, plus the
stderr assertion its three sibling rejection tests already carry. Also documented why
`--paired` is the wrong shape here, so it does not get "helpfully" restored.

### F2 — `phred64_ubam_out_specialty_stores_true_phred` could pass vacuously

It iterated `for q in &quals(&out)` with no emptiness check, unlike its SE and
clump-only siblings which both assert non-empty. A regression that produced a
header-only BAM would pass it silently. Added `assert_eq!(qs.len(), 4, …)`, which also
pins the record count the loop implicitly assumes.

Both fixes: green, `fmt` clean, `clippy -D warnings` clean.

---

## Test-guard audit (as requested)

| Test | Fails pre-fix? | Notes |
|---|---|---|
`bam.rs::bam_writer_subtracts_input_phred_offset` (Phred+64 arm) | Cannot compile pre-fix | Pins the arithmetic. Honestly labelled. Strong forward regression guard. |
same test, Phred+33 arm | **No** | Correctly labelled in-source as "an over-correction guard, not a bug guard". ✅ honest |
`phred64_ubam_out_se_stores_true_phred` | **Yes** — `[71×24, 33×10]` vs expected `[40, 2]` | Genuine bug guard |
`phred64_ubam_out_specialty_stores_true_phred` | **Yes** — raw 71 | Genuine (after F2) |
`phred64_bam_input_rejected_fastq_out` | **Yes** — pre-fix exits 0, 10/10 discarded | Asserts stderr names `--phred64` **and** `BAM` ✅ |
`phred64_bam_input_rejected_ubam_out` | **Yes** — the silent all-Q0 path | ✅ |
`phred64_bam_input_rejected_early_returning_mode` | **Yes** — flag inert pre-fix | Earns its place: `main.rs:271` records the prior mis-siting this guards against |
`phred64_mixed_input_rejected` | **No → yes after F1** | See F1 |
`phred64_clump_only_ubam_out_stores_true_phred` | **Yes** — `[71, 33]` | ✅ |
`phred64_clump_only_bam_input_rejected` | **Yes** — pre-fix this path worked correctly | Sharpest case in the issue; asserts stderr ✅ |

Net: 8 of 10 are genuine pre-fix-failing guards (after F1), 1 is an explicitly-labelled
over-correction guard, 1 is compile-gated. That is an honest ledger — the in-source
labelling of the Phred+33 arm as "not a bug guard" is exactly the discipline that was
asked for.

---

## Wiring audit

All 7 production `BamWriter::create` sites thread a dynamic offset; all 7 test sites pass
`33` deliberately. No site left on a literal where it should be dynamic, and none
converted that should have stayed literal:

```
src/main.rs:1852, 1994, 2126        cli.phred_offset()   (trim SE / PE-2-file / PE-interleaved)
src/specialty.rs:133, 183           input_phred_offset   (hardtrim5/3_to_bam)
src/clump_only.rs:772, 989          input_phred_offset   (clump-only SE / PE)
src/main.rs:327, 352, 518, 585, 617 cli.phred_offset()   (5 wrapper-function callers)
src/bam.rs:1395,1445,1473,1528,1568,1627,1744  33        (test callers)
src/clump_only.rs:1546              33                   (test caller)
```

19 edit points, matching §13's corrected count (the plan's headline "7" undercounts;
§13 already owns that).

### Guard coverage

`any_bam` is computed over **all** of `cli.input`, and every `BamReader` construction in
the codebase descends from a `cli.input` path via `format::open_sync_reader` /
`open_threaded_reader` / `BamReader::open_paired_interleaved*`. `--passthrough` is the
only third input stream and is opened exclusively through `FastqReader`
(`src/main.rs:1273, 1308`), so it cannot reach `BamReader`. Paired-interleaved,
multi-pair, mixed FASTQ+BAM, specialty, and clump-only are all downstream of line 211.
**No bypass found.**

### Comment accuracy at the Bug 1 site

Both former falsehoods are gone, and the replacements hold:

- `FastqReader::sanity_check` (`src/fastq.rs:472-503`) checks **only** the `@` prefix and
  colorspace digits — no quality-range validation. ✅
- `next_record` (`src/fastq.rs:390-427`) detects only truncation. ✅

So `saturating_sub` is correctly described as load-bearing rather than defensive. Nothing
new and false introduced. `PHRED_OFFSET`'s new doc comment is sound modulo L10.

---

## Low-priority

**L1 — invariant enforced at a distance.** `clump_only_single_to_bam`'s doc comment
(`src/clump_only.rs:723-727`) declares the load-bearing byte-identity invariant on
`R.qual`. That invariant now depends on `input_phred_offset == 33` whenever the input is
uBAM — enforced 500 lines away in `main.rs`. Nothing local records the dependency. PLAN
§7 says "record this in the PR body"; it belongs in the doc comment too, since that is
where the next reader of the invariant will look. Same for
`clump_only_paired_to_bam_one_pair`.

**L2 — the four new wrapper parameters are undocumented.** `hardtrim5_to_bam`,
`hardtrim3_to_bam`, `clump_only_single_to_bam`, `clump_only_paired_to_bam_one_pair` each
gained a bare `input_phred_offset: u8` with no doc bullet, while `BamWriter::create` —
the one that needed it least, being the site the compiler already forces — got a
six-line explanation. One line each pointing at `BamWriter::create` would do.

**L3 — no value-domain check on the offset.** `create` accepts any `u8`. `cli.phred_offset()`
can only return 33 or 64, so this is future-proofing, but note a `debug_assert!` would be
useless here: `[profile.release]` sets no `debug-assertions` and CI runs `cargo test
--release` (PLAN §3.2 makes this point about the existing `debug_assert_eq!`). If added,
it must be a `match`/`bail!`.

**L4 — `--length 0` is redundant.** Verified: `-q 0` alone preserves the full 34 bp
(length filter default 20 < 34). The test comment
(`tests/integration_ubam_out.rs:636-638`) and `test_files/README.md:14` both present
`-q 0 --length 0` as jointly necessary. Harmless as belt-and-braces, but the stated
reason is wrong for half of it.

**L5 — SE test does not pin the record count.** `phred64_ubam_out_se_stores_true_phred`
asserts `!qs.is_empty()`; 4 is the expected count. (F2 pinned this on the specialty test.)

**L6 — new singular CHANGELOG heading.** `#### Behaviour changes` is the only instance in
1141 lines, and it sits *after* `#### Fixes`, so the one entry that removes working
behaviour is the last thing a reader reaches. Its nearest sibling — `--dont_gzip` +
`--output-format ubam` now rejected, an identically-shaped "previously accepted, now an
error" item — lives in `#### Changes` at the top. Either move the entry there or accept
the new heading deliberately; two adjacent entries of the same kind under different
headings will confuse the release-notes pass.

**L7 — per-byte field load.** `.map(|b| b.saturating_sub(self.input_phred_offset))`
captures `&self` inside the closure. LLVM should hoist it, but `let offset =
self.input_phred_offset;` above the map costs nothing and removes the question. Mentioned
only because this repo has form on exactly this class (#248, #287).

**L8 — `docs/src/content/docs/reference/changelog.md` not synced.** It is a manual copy
(no sync step in `docs.yml`, which merely triggers on `CHANGELOG.md`) and already lags by
two `Unreleased` entries, so skipping it is consistent with practice. Flagged because the
release-time sync will otherwise propagate H1 to the published site.

**L9 — read-side overflow framed too narrowly.** CHANGELOG:79-82 and the in-code comment
attribute the overflow to sentinel mixing. Any raw byte ≥ 223 overflows; `0xFF` is the
plausible instance, not the only one. One clause ("…or any other out-of-spec raw score
above 222").

**L10 — `PHRED_OFFSET` doc nuance.** "Fixed at 33 by the SAM spec" is loose: the SAM
spec's ASCII+33 rule governs SAM *text*, while BAM stores raw — 33 is fixed here because
TrimGalore's internal `FastqRecord` is Sanger. The distinction is small but this comment
exists specifically to prevent conflating two offsets, so precision earns its keep.

---

## Confirmed non-issues

- **Backward compatibility.** Independently verified, not taken on trust: zero `phred64`
  hits in `.github/workflows/` or `scripts/`; `test_files/ubam_out_{se,pe}_REFERENCE.bam`
  unmodified in `git status`; the `assert_ubam_eq` reference tests pass unmodified; the
  README regen recipe uses BAM input with no `--phred64`, so the fixtures cannot move.
- **Uniform rejection over mode-dependent.** The right call, and the O-2 reconciliation
  (reject on *input format*, not on *mode indifference*) is a principle that will survive
  refactors. The `--clump_only` case is the one that makes it non-optional.
- **Required parameter over defaulted setter.** Vindicated in practice — §13 records two
  call sites missed by a regex pass and caught by the compiler.
- **Register.** Clean throughout. No conversational headings, no "win"/"just
  works"/"sweet spot". `:::caution` matches the Starlight asides used in
  `paired-end.md` and `clumpy.md`.
- **rustdoc.** No new intra-doc warnings (16 pre-existing; the new
  `BamWriter::input_phred_offset` and `PHRED_OFFSET` links resolve).
- **Efficiency.** No complexity or allocation change, as §6 states.

---

## Recommended order

1. **H1** — CHANGELOG FastQC claim. Release-note text asserting a symptom that does not occur.
2. **M1** — "all-zero" in three places, one of them a runtime `NOTE:`.
3. **M2** — decide: close the residual, or correct the comment that calls it harmless.
4. **M4**, **M3** — docs contradiction, then the README heuristic paragraph.
5. **L1**, **L2** — the two maintainability items with real half-lives.
6. Remaining L items at discretion.

None of these require re-running the plan-review gate. H1/M1/M4/M3 are prose; M2 is a
scoping decision the author is entitled to make either way provided the comment matches
the code.

---

*Report: `/Users/fkrueger/Github/TrimGalore/plans/phred64-ubam/CODE_review_reviewer-A.md`*
*Files modified by this review: `tests/integration_ubam_out.rs` (F1, F2).*
