# PLAN Review — Reviewer A

**Target:** `plans/phred64-ubam/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev` @ `472e6a2`
**Baseline:** `cargo check --all-targets` clean before review.
**Verdict:** **APPROVE WITH REVISIONS**

The mechanism is right, the scope is right, both bugs reproduce exactly as described, and the plan's one acknowledged unknown (Assumption 8) resolves cleanly in its favour. Three revisions are load-bearing: **invert the step ordering** (§5), **promote three manual checks to integration tests** (§9), and **delete the fixture-encoding heuristic** (§7) — it is wrong on three counts and adds nothing the grep check doesn't already prove.

---

## 1. Verified correct

Everything below I checked directly; the plan states it accurately.

**Both bugs reproduce.**
```
Bug 1  --phred64 --output-format ubam on a 34bp 'h'-qual FASTQ
       samtools renders QUAL as 'hhhhhhhhhhhh'  → stored raw 71   (expected 40)
       control (Phred+33 'I' input, no flag)    → 'IIIIIIIIIIII'  → stored raw 40
Bug 2  --phred64 -q 20 test_files/ubam_test.bam
       Reads too short:  10 (100.0%)   Reads written (passing filters): 0 (0.0%)
```

**Mechanism.** Subtracting `phred_offset` at the BAM boundary is correct. SAM/BAM `QUAL` holds true Phred `0..93`; `bam.rs:580` is the sole ASCII→raw conversion, and it is the only place in the codebase that must know the input encoding. `fastq.rs` contains no offset arithmetic at all (confirmed by grep), so the FASTQ pass-through contract is genuinely untouched.

**The premise "`FastqRecord.qual` reaching `BamWriter` is in the input's declared encoding" holds.** I traced it:
- The mode set that can reach `BamWriter` is bounded to three families — normal trim, `--hardtrim5/3`, `--clump_only`. `--clock`, `--implicon`, `--demux`, `--passthrough`, `--retain_unpaired`, `--clumpify` are all rejected at `src/cli.rs:553-588`.
- None of the three re-encode. `hardtrim5_to_bam`/`hardtrim3_to_bam` only `truncate`/`clip_5prime`; `clump_only` copies verbatim; the trimmer slices (`fastq.rs:88`, `:135`).
- The only production qual *synthesis* anywhere is `bam.rs:909` (`'!'` for missing BAM qual) — read-side only, and §3.1 rejects that combination. The other `"I".repeat(...)` sites (`adapter.rs:749`, `clump_only.rs:1142`) are test helpers.

**Scope counts (production).** Exactly 7, exactly as enumerated. `bam.rs:809` is indeed a doc comment, not a call. 6 test sites in `bam.rs` at exactly 1350, 1400, 1428, 1483, 1522, 1646.

**Read side.** `bam.rs:904-922` is correct as asserted — raw + 33 → Sanger ASCII, with the all-`0xFF` sentinel handled.

**The §3.1 guard-placement ordering claim.** Confirmed: `cli.validate()` at `main.rs:164`; the `detect_input_format` loop at `main.rs:182-186`; `any_bam` at `:187-189`. `Cli::validate()` genuinely cannot see the format. The cited precedent is real — the §3.4b rule sits at `main.rs:191-208` for exactly this reason.

**The guard site catches every route.** It sits at ~`main.rs:190`, upstream of specialty dispatch (`:263`), clump-only dispatch (`:400-573`), and all trim dispatch. Paired-interleaved, multi-pair, and mixed input all flow through the same `any_bam`. `sanity_check_any` at `:176` runs earlier but performs no quality interpretation (the BAM branch reads one record purely to test for emptiness). `naming::ensure_output_dir` is at `:235`, *after* the guard — so §9 check 5's "no output written" holds, including no output directory created. Good.

**Assumption 8 (§8) holds — this is the plan's main unknown and it resolves cleanly.** All four sites create the writer *directly* inside a function called from `main()` with `cli` in scope. Zero intermediate functions to thread through:

| Site | Function | Params today | Called from |
|---|---|---|---|
`specialty.rs:127` | `hardtrim5_to_bam` | 6 | `main.rs:275` |
`specialty.rs:175` | `hardtrim3_to_bam` | 6 | `main.rs:299` |
`clump_only.rs:766` | `clump_only_single_to_bam` | 10, already `#[allow(clippy::too_many_arguments)]` | `main.rs:557` |
`clump_only.rs:981` | `clump_only_paired_to_bam_one_pair` | 10, same allow | `main.rs:460`, `:526` |

One extra `phred_offset: u8` per function. No restructuring. Assumption 8 can be promoted from "unverified" to verified.

**No `--phred64` in tests or CI.** Repo-wide, the only `phred64` hits are `cli.rs:144-145`, `cli.rs:971`, `quality.rs:231` (a unit test of `quality_trim_3prime` in isolation — no BAM involvement), and two lines in `docs/.../quality.md`. Confirmed.

---

## 2. Critical

### C-1. §5's "Ordering note" is factually wrong, and the gap it hides is a real regression risk

The plan says:

> Step 8 is independent of 1–7 and could land first […] Implement 1–7 then 8

Step 8 is **not** independent — it is a *correctness prerequisite* for steps 1–7 on the BAM-in → BAM-out path. I verified this empirically:

```
$ ./target/release/trim_galore --clump_only --phred64 --output-format ubam \
      --output_dir $D test_files/ubam_test.bam
  → succeeds, 10 records written

$ samtools view test_files/ubam_test.bam        | cut -f11 | sort | md5
  31878eb235407ac33699756350a951a7
$ samtools view $D/ubam_test_clumped.bam        | cut -f11 | sort | md5
  31878eb235407ac33699756350a951a7      ← identical, i.e. CORRECT today
```

`--clump_only` never interprets quality (records round-trip verbatim), and `BamReader` normalises to Phred+33, so the current hardcoded `saturating_sub(33)` is *right* for this path — even with `--phred64` set. **Apply the writer fix without the guard and this command silently floors every quality byte to ~0.** It is the one path where step 1–7 alone makes correct output incorrect.

Three consequences:

1. **Invert the order — implement step 8 first.** The plan's stated rationale ("so the tree compiles at each step") is satisfied either way; step 8 has no compile dependency in either direction. Doing 8 first means no intermediate commit in the branch history silently corrupts a working path.
2. **This is the real justification for the one-PR decision** — stronger than the plan gives. The plan frames one-PR as a shape preference. It is actually a correctness constraint: shipping the writer fix without the guard would introduce a *new* silent-corruption bug. Say so in §10 and in the PR body; it pre-empts any reviewer request to split.
3. **§7 and the CHANGELOG need a caveat.** §7 says "No consumer of a Phred+33 run sees any change." That is true for output bytes but not for exit status: `--clump_only --phred64 <bam> --output-format ubam` goes from **working** to **hard error**. It is the only case in the whole change where the fix removes behaviour that currently works correctly. The rejection is still the right call (accepting the flag would require `BamWriter` to know that clump-only doesn't interpret quality — a fragile special case), but it must be stated plainly in the CHANGELOG so anyone with that invocation in a pipeline knows to drop the flag rather than debug an error.

§10 **O-2** is the place this belongs and currently misses it: O-2 reasons only about `--clump_only` with *FASTQ* input, concluding "the flag is inert where it is inert." With *BAM* input under clump-only the flag is equally inert today — and §3.1 rejects it anyway. Extend O-2 to say so and to record why rejecting-anyway is preferred.

---

## 3. Important

### I-1. §9 leaves all 7 production call sites without automated coverage

The bug class being fixed is, in the plan's own words, "a site that never considered it." A required `create` parameter forces every site to pass *something* — it does not stop an implementer writing `33` at a production site. Yet §9 checks 1, 3 and 4 (the only checks that exercise a production site end-to-end) are **manual shell steps**, and §9 concludes "a CI step is not required."

That reasoning ("the unit test pins the arithmetic more precisely than a shell assertion would") conflates two different things. The unit test pins `BamWriter`'s arithmetic. Nothing pins the *wiring* — which is where the bug lives.

The infrastructure is already there and this is cheap. `tests/integration_ubam_out.rs` already has `bam_tuples()`, which extracts quality as **raw Phred bytes**:

```rust
let qual: Vec<u8> = rec.quality_scores().as_ref().to_vec();
```

So three tests, each `Command::new(binary())` per the file's existing idiom, asserting `qual == vec![40u8; 34]` on the new Phred+64 fixture, cover all three families:

- trim SE → `main.rs:1797`
- `--hardtrim5 20` → `specialty.rs:127`
- `--clump_only` → `clump_only.rs:766`

(PE `main.rs:1938`/`:2069` and `specialty.rs:175`/`clump_only.rs:981` share the same writer and offset-plumbing shape; two of the three above plus the unit test is a reasonable stopping point, but a PE case is nearly free if the fixture is generated as a pair.)

**Recommendation: promote §9 checks 1/3/4 from manual to `tests/integration_ubam_out.rs`.** This is the single highest-value revision to the plan after C-1.

Note the proposed tests *are* genuine regression guards, not tautologies — I confirmed the current binary emits raw 71 on this exact input, so they fail on the buggy code (the offset-64 unit test additionally cannot compile pre-fix, and would fail if `bam.rs:580` were ever reverted to a literal). The offset-33 counterpart is tautological w.r.t. this bug, but the plan says so and it is worth keeping as a no-regression pin.

### I-2. §7's fixture-encoding heuristic is wrong on three counts — delete it

I re-ran the min/max quality-byte scan across every FASTQ fixture in `test_files/`. The plan's second check does not survive contact:

| Plan claim | Actual |
|---|---|
"min values of 33/35/45/62/63" | The set also includes **66** — `truncated.fq.gz` (min 66, max 99), omitted entirely |
"the four fixtures reporting min 73" — then lists **six** filenames | Only **three** have min 73: `BS-seq_10K_I1`, `PolyA`, `PolyT`. `polyAT_R1`, `polyAT_R2`, `illumina10K_with_polyA` all have **min 45**, which proves Phred+33 outright — no synthetic-uniform-Q40 reasoning needed for them |

And the one genuinely ambiguous fixture is the one the plan missed. `truncated.fq.gz` has **min byte 66** (`'B'` — the classic Illumina 1.5 read-segment-quality-control indicator) and **max byte 99** (`'c'`). Under Phred+64 that reads Q2..Q35, a textbook Illumina 1.5 profile. Under Phred+33 it reads Q33..Q66 — above anything an Illumina instrument emits. Its byte profile is **more** consistent with Phred+64 than with Phred+33. (Harmless in practice: it is a negative fixture used only to test truncation rejection.)

Separately, the heuristic's own logic is weaker than "proves outright": byte 63 rules out Illumina 1.3+ Phred+64 (floor 64) but sits inside the Solexa/Illumina-1.0 range, which allowed down to byte 59 for Q=−5.

**The conclusion is nonetheless sound** — but it rests entirely on check (a), which is airtight and self-sufficient: `phred_offset()` returns 33 unless `--phred64` is passed, and nothing in `src/`, `tests/`, or `.github/workflows/` passes it. The writer change is therefore provably a no-op for every existing test and fixture *regardless* of what encoding any fixture is in.

**Recommendation:** delete the min-byte check from §7 and lean on the grep. Bonus: this also retires the plan's own §11 "remaining risk" about documenting the new Phred+64 fixture in `test_files/README.md` — that risk exists only because §7 introduced a heuristic a future auditor might trust. (Still document the fixture; `test_files/README.md` has an established per-fixture provenance convention and the new file should follow it. Just no longer a *risk*.)

Also worth noting on the fixture question: committing `test_files/phred64_test.fastq` is the right call — no integration test in this repo generates its inputs inline, so a committed, README-documented fixture matches convention.

### I-3. The comment rewrite has a second false premise the plan doesn't flag

§5 step 3 correctly identifies that `bam.rs:576` ("`FastqRecord.qual` is Sanger ASCII (+33)") is the false premise behind the bug. The *same comment* has a second false claim, at `:577-579`:

> `saturating_sub` is defensive — under normal flow every byte is ≥33 (**FastqReader rejects sub-33 input**; BamReader synthesises '!' for missing qual, never sub-33).

`FastqReader` does **no** per-byte quality validation. `FastqReader::sanity_check` (`src/fastq.rs:474-503`) checks exactly two things: the id starts with `@`, and the sequence contains no ASCII digits (colorspace). The record-read path (`:382`, `:416`) only detects truncation. There is no quality-byte validation anywhere in `fastq.rs`.

So `saturating_sub` is guarding *real* malformed input, not a theoretical case. Fix both sentences while you're in there — the plan is rewriting this comment anyway, and leaving a second false invariant in place invites the same class of mistake next time.

This slightly weakens §3.2's edge-case reasoning, too. The plan says a Phred+33 file mistakenly run with `--phred64` fails "in the correct direction (visible zero-quality output, not silent garbage)". A BAM full of raw 0 renders as `!!!!!` under `samtools view`, which many pipelines will consume as legitimate Q0 — "visible" is optimistic. §10 **O-3**'s decision to emit a `NOTE:` for `--phred64` + uBAM output mitigates this and is the right call for exactly this reason; make that the stated rationale for O-3 rather than "consistent with that block's style."

### I-4. The docs target is the wrong file

§5 step 10 targets `docs/src/content/docs/guide/outputs.md` §uBAM output. Two problems:

1. The new rule is about **input**, and §3.1 point 3 explicitly applies it with FASTQ output too. An output-format page is the wrong home.
2. `outputs.md:88` heads its exclusion list with "Rejected at **CLI-validate time** (v1 scope)". The new guard is deliberately *not* at CLI-validate time (that is the whole point of §3.1's placement discussion). Filing it under that heading states something false about where the check lives.

**The primary home is `docs/src/content/docs/guide/quality.md` §"Phred encoding" (lines 29-31)** — where `--phred64` is actually documented, and where a user reaching for the flag will look. A one-line cross-reference in `outputs.md` is fine as a secondary, but it should not join the CLI-validate list.

---

## 4. Minor

- **M-1. Call-site count is 14, not 13.** `src/clump_only.rs:1531` is a unit-test caller of `clump_only_single_to_bam` that must be updated when the function gains a param. Compiler-caught (baseline `cargo check --all-targets` is clean, so the compiler *will* catch it), so no correctness risk — but §4 and §2 are precise about counts and this one is short.
- **M-2. Line-reference drift.** `phred_offset()` is at `cli.rs:970-971` (plan says 968-970). The uBAM startup-`NOTE` block that O-3 targets is at `main.rs:237-260` (plan says 218-234). Cosmetic, but O-3's implementer will go to the wrong place.
- **M-3. Test-count baseline.** 409 `#[test]` attributes across `src/` + `tests/`, not 407. §9 check 8's "count ≥ current 407" would silently pass a regression that deleted two tests. Use the actual number, or drop the numeric gate in favour of "no test removed."
- **M-4. §9 check 4 is under-specified.** `--clump_only --phred64 --output-format ubam` needs "on a Phred+64 **FASTQ**" — with BAM input the same command is rejected by §3.1 (see C-1), so as written the check is ambiguous about which outcome it expects.
- **M-5. §9's unit test doesn't say at which layer to read back.** "read back, assert stored raw 40" — prefer raw-BAM-layer inspection via `bam::io::Reader` + `bam::Record::quality_scores()`. There is an in-file precedent at `src/bam.rs:1522` (`bam_writer_aux_typed_int_and_float_round_trip` does exactly this to avoid the reader's normalisation). A `BamReader::open` round-trip also catches the bug (`'h'` in → raw 40 → `'I'` out, vs `'h'` out today) but asserts the property indirectly.

---

## 5. Open questions in §10 — assessment

- **O-1** (`--phred33 --phred64` both → 64, no `conflicts_with`). Confirmed at `cli.rs:144-149` and `:970-971`. Deferring is defensible: it is pre-existing, not uBAM-specific, and fixing it changes FASTQ-path behaviour on a path no CI covers. Interaction with the new guard is benign (passing both with BAM input still rejects, since `self.phred64` is true). No hidden problem. Agree with the assumption taken.
- **O-2.** Hides a real problem — see **C-1**. The reasoning about FASTQ input is sound; the BAM-input case it omits is the one behaviour regression in the change. Must be extended.
- **O-3.** Assumption taken (include the `NOTE:`) is right, and for a better reason than the plan gives — see **I-3**. Also fix the line reference (**M-2**).

---

## 6. Efficiency

§6 is correct and needs no revision. A field load in place of an immediate, on a loop the optimiser will hoist it out of, on a path dominated by BGZF compression. No allocation change. The guard is two boolean tests once per run on already-computed values. Nothing to add.

---

## 7. Out of scope — one latent bug on the neighbouring line

Not this plan's problem, but it is two lines from a site the plan touches and it is the same defect family, so flag or file it:

`src/bam.rs:920` does `(b + PHRED_OFFSET)` with a plain `+`. The all-`0xFF` sentinel check at `:908` only catches records where **every** byte is `0xFF`. A BAM with a *mixed* QUAL — some real values, some `0xFF` — slips through and then overflows `u8`: panic in debug, wrap to 32 (`' '`) in release. Either fold a `saturating_add` in here (one character, zero risk, symmetric with the write side this plan is fixing) or open a follow-up issue.

---

## 8. Action items

**Critical**
1. Invert §5's ordering — implement step 8 (guard) **before** steps 1–7; delete the "Step 8 is independent of 1–7" sentence. **(C-1)**
2. Record in §10/O-2 and the PR body that the guard is a *correctness prerequisite* for the writer fix on BAM-in → BAM-out, with the verified `--clump_only --phred64 <bam> --output-format ubam` case as evidence. This is the real justification for the one-PR decision. **(C-1)**
3. Add a CHANGELOG note that `--clump_only --phred64 <bam> --output-format ubam` changes from working to a hard error — the only behaviour in the change that goes from correct to rejected. **(C-1)**

**Important**
4. Promote §9 checks 1/3/4 from manual shell steps to integration tests in `tests/integration_ubam_out.rs`, using the existing `bam_tuples()` raw-Phred extraction. Three tests, one per writer family. **(I-1)**
5. Delete §7's min-quality-byte fixture check; rely on the grep check alone. Downgrade the §11 fixture-documentation "remaining risk" accordingly (still document the fixture in `test_files/README.md`). **(I-2)**
6. In §5 step 3, fix **both** false premises in the `bam.rs:576-579` comment — including "FastqReader rejects sub-33 input", which is false. Restate O-3's rationale as mitigating the silent-floor-to-Q0 case. **(I-3)**
7. Retarget §5 step 10 to `docs/src/content/docs/guide/quality.md` §"Phred encoding"; keep `outputs.md` as an optional cross-reference, and not inside its "Rejected at CLI-validate time" list. **(I-4)**

**Optional**
8. Correct the call-site count to 14 (add `clump_only.rs:1531`). **(M-1)**
9. Fix line references: `cli.rs:970-971`, `main.rs:237-260`; test baseline 409. **(M-2, M-3)**
10. Specify "Phred+64 FASTQ input" in §9 check 4; specify raw-BAM-layer read-back in the unit test, citing the `bam.rs:1522` precedent. **(M-4, M-5)**
11. Promote Assumption 8 from "unverified, low risk" to verified — all four sites take one extra parameter, zero intermediates. **(§1)**
12. Fold in or file the `bam.rs:920` `+` → `saturating_add` overflow on mixed-`0xFF` QUAL. **(§7)**
