# PLAN Review — Reviewer B

**Target:** `plans/phred64-ubam/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore` @ `dev` / `472e6a2`
**Issue:** #358
**Method:** every structural claim in the plan re-derived from source; both bugs and all edge cases below reproduced with the prebuilt `target/release/trim_galore` + `samtools`.

## Verdict: **APPROVE WITH REVISIONS**

The design is right and I could not break it. The split (guard on read, offset on write) is the only correct shape, the guard's placement in the T17 block is forced and verified, the required-parameter choice is well argued, and Assumption 8 resolves cleanly at every one of the four sites it hedges on. What the plan is missing is not design — it is **the causal link between its two halves**. The plan presents the guard as an independent improvement bundled for convenience (§5 "Ordering note"). It is not: the Bug 1 fix *creates* a silent all-Q0 corruption path on BAM input, and the guard is the only thing that blocks it. Because the plan doesn't see that, its highest-risk combination is untested (C-2) and its rationale for "hard error, not warning" is weaker than the facts support (C-1).

Everything below is additive. No task in §5 needs to change.

---

## Critical

### C-1. The two fixes are coupled in the correctness direction — §5's "Ordering note" has it backwards

§5 states: *"Step 8 is independent of 1–7 and could land first … step 8 alone leaves Bug 1 live but is not a compile dependency."* True about compilation, wrong about correctness, and the dependency runs the opposite way.

Reproduced on `dev`:

```
$ trim_galore --phred64 -q 0 --length 0 --output-format ubam test_files/ubam_test.bam
Reads written (passing filters):         10 (100.0%)
$ samtools view out/ubam_test_trimmed.bam | head -1 | cut -f11
???????????????????????????????????????????????????????????????   # raw 30 — CORRECT

$ trim_galore --clump_only --phred64 --output-format ubam test_files/ubam_test.bam
$ samtools view out/ubam_test_clumped.bam | head -1 | cut -f11
?????????????????????????????????????????????????????????????????  # raw 30 — CORRECT
```

Both of these **work today**. Bug 2 only bites when quality trimming is active (`-q` > 0 discards the reads loudly). With `-q 0`, or on `--clump_only` where no quality arithmetic exists at all, records pass through and the writer's `saturating_sub(33)` on Phred+33 ASCII produces the right answer by accident.

Apply steps 1–7 without step 8 and the writer subtracts 64 from the reader's Phred+33 bytes (`'?'` = 63) → saturates to raw **0** → **an all-Q0 BAM, every read retained, no warning, no non-zero exit**. That is strictly worse than either bug being fixed: Bug 2's current failure is loud (100% discarded), and the new one is silent.

Three revisions follow:

1. **Rewrite the §5 ordering note.** The guard is a *prerequisite* for the writer fix on BAM-input paths, not an independent addition. State that it can never be relaxed to a warning-and-ignore later without re-opening a silent-corruption path.
2. **Strengthen Assumption 5.** Its rationale ("the current behaviour destroys 100 % of reads") is the pre-fix argument. The post-fix argument is stronger and should be the one recorded: warn-and-ignore would silently destroy *quality* rather than loudly destroy *reads*.
3. **Say it in the PR body.** A reviewer looking at the diff will see two changes that appear separable and may ask why they're in one PR. This is the answer.

### C-2. The highest-risk combination has no test

Validation checks 5 and 6 exercise `--phred64` + BAM input with **FASTQ output** only. Nothing in §9 covers:

- `--phred64 --output-format ubam test_files/ubam_test.bam`
- `--clump_only --phred64 --output-format ubam test_files/ubam_test.bam`

These are exactly the two invocations where a guard regression stops being a visible read-loss and becomes the silent all-Q0 BAM of C-1. Add both to the `tests/integration_ubam_out.rs` rejection set (the existing idiom at `tests/integration_ubam_out.rs:443` `ubam_out_two_bam_pair_rejected` is the right shape — `.output()`, `!status.success()`, stderr substring). For the clump-only one, `tests/integration_clump_only_ubam.rs` is the better home.

Cost: two ~12-line tests. Value: they are the only guards standing between a future refactor and silent quality destruction.

### C-3. The Bug 1 fix violates `clump_only`'s documented lossless invariant on BAM input; only the guard restores it

`src/clump_only.rs:721-726`, the doc comment on `clump_only_single_to_bam`:

> *"Load-bearing invariant: every input record R appears in the output BAM with `R.id` (name + preserved aux tags), `R.seq`, `R.qual` byte-identical through the FASTQ-record intermediate."*

Under `--phred64` + BAM input, the writer fix turns `R.qual` into all-zeros (C-1) — a direct breach of the mode's headline property, and `--clump_only`'s entire pitch is losslessness. §3.1's guard is what keeps the invariant intact.

The plan touches `clump_only.rs` in step 5 but never mentions the invariant. Add a sentence to §7 (Integration) recording that the guard is what preserves it. Whoever reviewed #356 will look for this, and it is the kind of claim that shouldn't be left implicit.

---

## Important

### I-1. §7's fixture-encoding claim is false — `test_files/truncated.fq.gz` is Phred+64

§7 asserts *"Every fixture in `test_files/` is Phred+33"* and offers it as one of two independent backward-compatibility checks. I scanned every FASTQ fixture's min/max quality byte independently. The claim does not hold:

```
min= 66 max= 99 recs= 10000  truncated.fq.gz
```

`truncated.fq.gz` is Illumina 1.5 Phred+64, and the evidence is unambiguous:

```
@HWUSI-EAS611_0001:5:1:1028:17510#0/1
NATGGGCTGTGTAGCTCTATGTGCAACCAACTGTCAAATG
+HWUSI-EAS611_0001:5:1:1028:17510#0/1
BIJJIGIKKI``````````VVYVW`````BBBBBBBBBB
```

- `HWUSI-EAS611_...#0/1` read names are Illumina GA/GAII vintage — precisely the era that shipped Phred+64.
- The trailing `B` runs are the **Illumina 1.5 read-segment quality control indicator**: `B` = 66 = Q2 at offset 64. The backtick `` ` `` = 96 = Q32 at offset 64.
- Under Phred+33 those bytes read as Q33 and Q63. Q63 does not exist in Illumina data.

§7 also miscounts the ambiguous set. It names six files as *"the four fixtures reporting min 73"*; the actual min-73 set is three — `PolyA`, `PolyT`, `BS-seq_10K_I1`. `polyAT_R1`, `polyAT_R2` and `illumina10K_with_polyA` all report min **45**, i.e. proven Phred+33 outright by the plan's own test.

**Functional impact: none.** `truncated.fq.gz` has *zero* references in `src/`, `tests/`, `.github/` — it is an orphan fixture (`CLAUDE.md` still lists it among the negative cases), so it never reaches `BamWriter`. And §7's *other* check carries the conclusion on its own: I verified `grep -rn "phred64" tests/ .github/workflows/` returns zero hits, so nothing in the test suite or CI can be affected by an offset-conditional change. The no-golden-churn conclusion stands.

But two things need fixing:

- **Correct §7's wording.** Presenting a false claim as one of two independent checks weakens the argument it's meant to support, on the one axis (quality encoding) a reviewer will read hardest.
- **§9 tells the implementer to record fixture encodings in `test_files/README.md`.** Writing "all fixtures are Phred+33" into the repo would commit a demonstrably false statement about in-tree data. If the README is going to carry an encoding note, `truncated.fq.gz` must be listed as the Phred+64 exception.

Bonus for the PR body: this file is in-tree proof that real Phred+64 data is what #358's reporter is working with. It also means a Phred+64 fixture already exists — just an unusable one (deliberately truncated), so §9's new fixture is still needed.

### I-2. §3.1's stated justification is wrong for the modes that never interpret quality — and it breaks currently-working invocations

§3.1 point 3 justifies rejecting regardless of output format with *"the defect is on the **read** side."* That is true for the trimming pipeline. It is **not** true for the modes that do no quality arithmetic. Both of these succeed today:

```
$ trim_galore --hardtrim5 20 --phred64 test_files/ubam_test.bam
Finished writing 10 sequences                       # --phred64 completely inert

$ trim_galore --clump_only --phred64 test_files/ubam_test.bam
                                                    # --phred64 completely inert
```

There is no read-side defect on these paths — the flag is a no-op. The plan will convert two working invocations into hard errors.

I think rejecting uniformly is still the right call, but for a different reason: a mode-dependent rule ("rejected unless you're in hardtrim or clump-only-to-FASTQ") is worse to document, worse to test, and one refactor away from being wrong. Restate the justification as *"the flag is meaningless for BAM input in every mode; a uniform rule is safer than a mode-dependent one"* and drop the read-side claim as the reason.

Two consequences to record:

- **CHANGELOG needs a behaviour-change line, not just a fix line.** A pipeline wrapper that passes `--phred64` unconditionally for legacy data will hard-fail the first time an input becomes uBAM, even in modes where the flag never did anything. That's a legitimate break and users should read it in the changelog rather than discover it.
- **This contradicts O-2's principle.** O-2 declines extra rejection for `--clump_only` + FASTQ on the grounds that *"the flag is inert where it is inert."* §3.1 rejects `--clump_only` + BAM where the flag is *equally* inert. Both decisions are defensible; holding both principles simultaneously is not. Pick one and name it — I'd suggest "reject where the flag is meaningless **and** the input format makes it so; leave inert-but-harmless FASTQ cases alone" — and note that O-2's case survives because FASTQ genuinely *has* an ASCII encoding to declare, whereas BAM does not.

### I-3. `bam.rs` already has a `PHRED_OFFSET` constant — the plan never mentions it

```
src/bam.rs:60   const PHRED_OFFSET: u8 = 33;
src/bam.rs:920      .map(|&b| (b + PHRED_OFFSET) as char)     # read side, its only use
src/bam.rs:580  ... .map(|b| b.saturating_sub(33))            # write side, bare literal
```

The module already has a name for the read-side offset and the write side simply didn't use it. Post-fix the file will carry two near-identical names with **different semantics**:

- `PHRED_OFFSET` — fixed 33, the offset the *reader adds* (correct, must stay fixed)
- `BamWriter.phred_offset` — 33 or 64, the offset the *writer subtracts*

That's a maintenance trap: the next person to touch quality handling in this file will reasonably assume they're the same thing, and "conflating two offsets" is the exact bug class being fixed. Recommend naming the field `input_phred_offset`, or at minimum a one-line comment at the const stating it is read-side-only and unrelated to the writer field.

Assumption 7's rationale ("rather than a named constant") should also acknowledge that a named constant already exists in the same module — the reasoning still holds, but as written it reads as though there wasn't one.

### I-4. The comment §5 step 3 rewrites carries a *second* falsehood the plan doesn't name

`src/bam.rs:576-579`:

```rust
// FastqRecord.qual is Sanger ASCII (+33); BAM stores raw Phred bytes.
// `saturating_sub` is defensive — under normal flow every byte is
// ≥33 (FastqReader rejects sub-33 input; BamReader synthesises '!'
// for missing qual, never sub-33).
```

The plan correctly identifies the first clause as the false premise behind the bug. The parenthetical is **also false**. `FastqReader` performs no such rejection — there is no quality-range validation in `src/fastq.rs` at all. Verified:

```
input qual: II IIIIIIIII...        (ASCII 32 space at position 3)
output QUAL: II!IIIIIIIII...       (raw 0 — reached the writer, saturated)
```

So `saturating_sub` is **load-bearing, not defensive**, even at offset 33. The plan's §3.2 edge case inherits the wrong framing ("Its defensive role is unchanged"). Name both halves in §5 step 3, so the rewritten comment doesn't carry the second falsehood forward — this is the one place in the codebase where quality-encoding assumptions get written down, and it has already been wrong once.

### I-5. `debug_assert_eq!` is compiled out in the profile CI uses

`Cargo.toml:58-62` sets no `debug-assertions` in `[profile.release]`, and `.github/workflows/ci.yml:74` runs `cargo test --release`. The length assertion at `src/bam.rs:581` therefore **never executes in CI**. §3.2 cites it as an edge-case safety net ("`debug_assert_eq!` on length still holds") and §11 lists empty-qual among "edge cases covered" partly on that basis.

Not a blocker — the fix cannot change lengths, and I confirmed the reachable degenerate case is benign. A read fully consumed by adapter trimming with `--length 0` does reach the writer:

```
$ trim_galore --phred64 --length 0 -a AGATCGGAAGAGC --output-format ubam <fastq>
name=allad seq=* qual=*        # zero-length record, noodles emits '*'/'*' — no panic
name=keep  seq=TTTT...         # 33bp survivor
```

Just stop citing an inert assertion as coverage. If the invariant is worth enforcing, `assert_eq!` or an explicit `bail!` would actually run.

---

## Minor / Optional

**M-1. Genuine Solexa / Illumina 1.0 data floors to Q0 post-fix.** Reproduced: a qual line of `;` (ASCII 59) currently stores raw 26; post-fix it saturates to raw 0. `--phred64` is precisely the flag old-data users reach for, and true Solexa files legitimately contain ASCII 59–63 (negative Solexa scores). Post-fix is still strictly better than today's silent Q26 fabrication and matches `quality.rs` exactly, so no change to the fix — but §3.2's edge-case list names only "a Phred+33 file mistakenly run with `--phred64`" and omits the more plausible occurrence. One line.

**M-2. Validation checks 3 and 4 don't state their input.** Both must use **FASTQ** input — with BAM input the §3.1 guard rejects them. An implementer who reaches for `test_files/ubam_test.bam` on check 4 will see a rejection and may conclude the guard is over-broad and "fix" it. Name the input explicitly.

**M-3. Which tests actually fail pre-fix.** Checks 2 and 9 pass before *and* after by construction. That's intentional — they're no-op proofs, and §9 says so for 9 but not for 2. Worth being explicit so nobody mistakes them for bug guards. The genuinely pre-fix-failing artefacts are: the offset-64 unit test (can't compile pre-fix — say so, it's the arithmetic pin), checks 1/3/4, and rejection checks 5/6/7. The offset-33 unit-test counterpart is an over-correction guard, not a bug guard; label it as such.

**M-4. No test asserts the guard fires for the early-returning modes.** I verified by reading that the T17 block (`main.rs:178-221`) precedes specialty dispatch (`main.rs:256`) and clump-only dispatch (`main.rs:349`), so the guard does cover them. But this exact class of bug has already bitten this codebase — `main.rs:237-240` records code-review round-2 finding B-NIT-2: the uBAM startup NOTEs *"never fired on the hardtrim BAM paths (which return earlier)."* Given that precedent, add `--hardtrim5 20 --phred64 <bam>` to the rejection tests. Cheap insurance against a future guard being placed one block too late.

**M-5. The report prints `ASCII+64` beside a BAM holding raw Phred.** Verified: `--phred64 --output-format ubam` writes *"Quality encoding type selected: ASCII+64"* into the trimming report (`main.rs:1855` → `report.rs:212-215`). That correctly describes the **input**, so it isn't a bug — but a user reading the report next to the BAM could conclude the BAM is ASCII+64. If O-3's NOTE lands, wording it as *"input is Phred+64; output BAM `QUAL` stores true Phred scores per the SAM spec"* resolves both at once for one `eprintln!`. Recommend adopting O-3 with that wording.

**M-6. FastQC on uBAM output is an unlisted beneficiary of the fix.** `--fastqc` with `--output-format ubam` runs fastqc-rust directly on the `.bam` (`main.rs:1887-1896`; `docs/…/outputs.md:107`), and fastqc-rust reads BAM `QUAL` as raw Phred. So pre-fix, `--phred64 --output-format ubam --fastqc` yields per-base-quality plots centred on **Q71**; post-fix they're correct. The plan doesn't mention this consumer at all. Add one validation row — a cheap end-to-end confirmation on a user-visible artefact that the unit test can't reach — and one CHANGELOG clause.

**M-7. Archived bad output cannot be recovered by re-running TrimGalore.** A pre-fix `--phred64` uBAM stores raw = true + 31. Reading it back adds 33 → Q71, and `--phred64` (which would subtract 64, not 31) is now rejected anyway. §11's informational note says users "know to regenerate" — make it specific: regenerate **from the original FASTQ**; re-processing the bad BAM cannot fix it, in either direction.

**M-8. §4's doc-comment block silently deletes existing history.** The current `create` doc comment (`src/bam.rs:524-527`) carries a *"**PLAN §4 deviation:**"* paragraph explaining why `command_line` is a 4th parameter. §4's replacement block omits it. An implementer copying §4 verbatim drops a deliberate historical note. Keep it, and consider adding the same treatment for `phred_offset` as the 5th.

**M-9. Accuracy nits.**
- §2/§5 claim `#[allow(clippy::too_many_arguments)]` precedent covers both files. `src/specialty.rs` has **none** — and needs none: `hardtrim5_to_bam` goes 6 → 7 params, at clippy's threshold rather than over it. `src/clump_only.rs:733` and `:881` do already carry it. Outcome fine, claim inaccurate.
- Validation 8 says "count ≥ current 407"; `dev` @ `472e6a2` carries **409** `#[test]` attributes. Use the figure `cargo test` reports, or drop it.
- §2 cites `src/cli.rs:968-970` for `phred_offset()`; the fn body is at 970-972.
- `CHANGELOG.md` has `### Unreleased` → `#### Changes` but no `#### Fixes` heading yet; step 9 creates it.

---

## Verified correct — don't re-derive

Everything below I checked independently and found accurate, so the implementer can take it as settled:

- **All 13 call sites exist at the exact claimed lines.** 7 production (`main.rs:1797, 1938, 2069`; `clump_only.rs:766, 981`; `specialty.rs:127, 175`) and 6 test (`bam.rs:1350, 1400, 1428, 1483, 1522, 1646`). The 8th `bam.rs` hit is indeed a doc comment (`bam.rs:809`), as §11 says.
- **Assumption 8 resolves cleanly — no restructuring anywhere.** `run_ubam_output_single` (`main.rs:1783`) already takes `cli: &Cli`, so `cli.phred_offset()` is immediate. `hardtrim5_to_bam` / `hardtrim3_to_bam` (`specialty.rs:106, 152`) and `clump_only_single_to_bam` (`clump_only.rs:734`) / `clump_only_paired_to_bam_one_pair` (`clump_only.rs:882`) take flat parameter lists with no `cli` — one added `phred_offset: u8` plus caller updates in `main.rs`. **No intermediate function to thread through** in any of the four. Assumption 8 can be marked resolved rather than left open.
- **Both bugs reproduce exactly as described.** Bug 1: trim SE path and `--hardtrim5` path both store raw 71 for input `'h'`; Phred+33 control stores raw 40. Bug 2: `--phred64 -q 20 test_files/ubam_test.bam` → 10/10 "too short", 0 written.
- **Guard placement is forced and correct.** `Cli::validate()` runs at `main.rs:164`; format detection at `main.rs:181`. The guard genuinely cannot live in `cli.rs` §3.4a. The T17 block precedes *every* dispatch branch. §3.4b at `cli.rs:544-546` is a real precedent for the pattern.
- **`grep -rn "phred64" tests/ .github/workflows/` → zero hits.** Verified. This single check is sufficient for the no-churn conclusion (see I-1).
- **Golden fixtures unaffected.** `ubam_out_se_REFERENCE.bam` / `ubam_out_pe_REFERENCE.bam` are generated without `--phred64` (regen recipe in `test_files/README.md`), so an offset-conditional change cannot move them. Note the README lists *"a bug fix in seq/qual handling"* as a regen trigger — this one qualifies textually but not materially; worth one line in the PR body pre-empting the question.
- **v0.6.11 byte-identity matrix genuinely untouched.** `bam.rs` holds the only ASCII↔raw conversion in the codebase — the whole file has exactly two qual-arithmetic sites (580 write, 920 read) — and `src/fastq.rs` performs no offset arithmetic at all, so the FASTQ pass-through contract in Assumption 3 is confirmed by absence.
- **Round-trip integrity is *improved*, not compromised.** Post-fix, `--phred64 <fastq> --output-format ubam` writes true Phred; reading that BAM back needs no flag (the reader normalises to +33). There is no legal "Phred+64 BAM" — the SAM spec forbids it — so refusing `--phred64` on BAM input creates **no** write/read asymmetry. The only BAMs the fix "refuses to re-read correctly" are pre-fix TrimGalore outputs, and `--phred64` couldn't fix those anyway (they need −31, not −64). The rejection is not over-reach on this axis.
- **`saturating_sub`'s failure direction is the right one.** A Phred+33 file run with `--phred64` floors to raw 0 — a visibly absurd all-`!` BAM that every downstream QC tool flags immediately. The alternative (wrapping, or the current silent inflation into the legal 0..93 window) is what made this bug survive to a user report in the first place. Keep it. Flooring is loud; 71-instead-of-40 was silent.
- **The proposed fixture works unaided.** 4 × 34 bp all-`'h'` runs fine with no `-a` — auto-detect falls back to Illumina, 0 adapters found, 4/4 written. `test_files/README.md` exists (44 lines, sectioned) and has a natural home for the encoding note.
- **O-1 confirmed:** `--phred33 --phred64` together → report says `ASCII+64`.
- **Empty input can't reach the writer:** `sanity_check_any` bails with *"Input file … is empty"* well before format detection.

---

## Action items

**Critical — resolve before implementing**
1. Rewrite §5's ordering note: the guard is a correctness prerequisite for the writer fix on BAM-input paths, not an independent addition. Strengthen Assumption 5 with the post-fix rationale. (C-1)
2. Add rejection tests for `--phred64 --output-format ubam <bam>` and `--clump_only --phred64 --output-format ubam <bam>`. (C-2)
3. Record in §7 that the guard is what preserves `clump_only`'s documented lossless-qual invariant. (C-3)

**Important**
4. Correct §7's fixture claim — `truncated.fq.gz` is Phred+64 — and fix the min-73 miscount. Ensure the §9 `test_files/README.md` note doesn't commit the false claim. (I-1)
5. Restate §3.1 point 3's justification; add a CHANGELOG behaviour-change line for the inert-flag modes; reconcile with O-2's principle. (I-2)
6. Rename the field `input_phred_offset` or document it against the existing `PHRED_OFFSET` const. (I-3)
7. In §5 step 3, name *both* falsehoods in the comment being rewritten — including "FastqReader rejects sub-33 input", which is also untrue. (I-4)
8. Stop citing `debug_assert_eq!` as edge-case coverage; it's compiled out under `cargo test --release`. (I-5)

**Optional**
9. M-1 Solexa note · M-2 name check 3/4's input · M-3 label which tests fail pre-fix · M-4 hardtrim rejection test · M-5 adopt O-3 with SAM-spec wording · M-6 FastQC-on-BAM validation row · M-7 "regenerate from FASTQ" · M-8 keep the §4 deviation paragraph · M-9 line-number and count nits.

**Mark resolved:** Assumption 8 — verified, one added parameter at each of the four sites, no intermediate threading.
