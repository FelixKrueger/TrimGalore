# Plan Review A — `plans/10202026_polya-behind-polyg/PLAN.md`

**Reviewer:** A (independent, fresh context)
**Date:** 2026-10-20
**Tree:** `/Users/fkrueger/Github/TrimGalore`, branch `wiki/init` @ `d471269` (`src/` at `dev` state)
**Method:** read `src/trimmer.rs`, `src/quality.rs`, `src/alignment.rs`, `src/cli.rs`, `src/library.rs`,
`.github/workflows/ci.yml`, `docs/.../migration.md`, `docs/.../poly.md`, `git show 0.6.10:trim_galore`;
plus four throwaway probe tests (written, run, deleted) that re-derived both of the plan's number tables
and swept the candidate over five real fixtures.

**Verdict: do not implement as written.** B1 is correct, verified, and sufficient for the reported bug.
**B2 as specified is a regression larger than the bug it fixes** — it newly trims 21–25% of reads on
ordinary libraries and cuts 142 bp off a read in the very fixture the plan uses to disqualify its own
rejected alternative. Recommend shipping B1, dropping B2, and replacing it (if anything) with the
terminal-homopolymer-run skip measured in §5 below.

---

## 1. What I verified and confirmed (plan is right)

| Claim | Status |
|---|---|
| `trimmer::trim_read` is the only implementation of the stage order | **VERIFIED.** `grep` for `trim_read(` gives only `src/parallel.rs:495,496,1133` plus in-module tests. `grep` for `homopolymer_trim_index\|poly_a_trim_index` across `src/` gives only `src/trimmer.rs:200` and `:229`. `specialty.rs`, `clump*.rs`, `bam.rs` hold no copy. One edit site, as stated. |
| uBAM output path is covered by the same edit | **VERIFIED.** `run_single_end_to_bam` / `run_paired_end_to_bam` call `trim_read` at `src/trimmer.rs:645` / `:730,731` — inside the same function, not a parallel implementation. |
| `--clump_only` is out of scope | **VERIFIED, and stronger than the plan says.** `src/cli.rs:863-871` *rejects* `--poly_a` and `--poly_g` under `--clump_only` at the CLI layer. The stage can never run there. |
| `bp_after_cutadapt` / the #232 invariant is untouched | **VERIFIED.** Captured at `src/trimmer.rs:161`, before RRBS (2.5) and both poly stages. Reordering 2.7/2.8 cannot move it. |
| No pinned CI output can move | **VERIFIED, no counterexample found.** `grep -rn "poly_a\|polyA\|poly-a" .github/workflows/` returns only two hits, both `--no_poly_g` on the demux step (`ci.yml:967-968`). `poly_a: true` appears nowhere in `src/`, `tests/`, or `examples/`. All six `TrimConfig` construction sites in test code (`trimmer.rs:842`, `trimmer.rs:1111 poly_g_config`, `report.rs:1353`, `parallel.rs:1425`, `parallel.rs:2168`) set `poly_a: false`; there is no `Default` impl to leak a `true` through `..Default::default()`. `src/library.rs` (the `--library` presets) does not mention poly at all — I checked specifically because a preset would have been the one invisible counterexample. The only `true` source is `src/main.rs:1294 poly_a: cli.poly_a`, default false. **The plan's riskiest-sounding claim is the one that holds up best.** |
| 509 lib tests is the right baseline for Validation 4 | **VERIFIED.** `cargo test --release --lib` → `509 passed`. |
| Perl `--polyA` defaulted to `-a "A{20}"` / `-a2 "T{150}"` | **VERIFIED** at `0.6.10:trim_galore:2754,2761,3357`. |
| `--2colour` workaround mechanism | **VERIFIED.** `quality::quality_trim_3prime_nextseq` (`src/quality.rs:61-90`) scores G as `cutoff-1`, step 1 runs before everything, and I hand-traced both G cases to 60 and 40. (Pedantic: it is `cutoff-1`, i.e. "+1 to the running sum", not literally "quality 0" — matching Cutadapt. The plan's prose simplification is harmless.) |
| B1's "R2 is unchanged" claim | **VERIFIED EMPIRICALLY — upgrade this from inference to fact.** The plan justifies it as "disjoint ends", which is *not* a proof: both stages mutate the same record and a short read can have both ends consumed. I probed it properly: 33,012 reads across `polyAT_R2`, `PolyT`, `BS-seq_10K_R2`, `illumina_10K` plus eight adversarial degenerate cases (`TTTGGG`, `T15G15`, `T20G3`, `T3G20`, `T5G3T5G5`, …) — **zero divergence** in output length *or* in either stat counter. Mechanism: the 5' `revcomp` gate is `errors*5 <= i+1`, which depends on `i` only, never on `n`, so a 3' truncation cannot perturb the 5' scan. Keep Validation 5; the claim is sound. |

---

## 2. CRITICAL — B2 is a data-loss regression on ordinary libraries

The plan rests B2's safety on a single synthetic row — `| internal 16-A run | 76 | 76 | no over-trim |`
— and a Behavior-section assertion that "B2 does not re-open the v0.6.x over-trim". I implemented the
candidate exactly as specified (`for offset in 0..=10`, scan `seq[..n-offset]`, first offset that cuts
wins) and swept it over real fixtures. The plan's synthetic table reproduces exactly; **the real-data
behaviour does not match the plan's conclusion.**

| Fixture | reads | current trims | candidate trims | **newly trimmed (0 → >0)** | extra bp | worst extra cut |
|---|---|---|---|---|---|---|
| `illumina_10K.fastq.gz` | 10,000 | 1,537 (15.37%) | 4,047 (40.47%) | **2,510 (25.10%)** | 30,519 | **142 bp** |
| `illumina10K_with_polyA.fastq.gz` | 10,900 | 2,437 (22.36%) | 4,947 (45.39%) | **2,510 (23.03%)** | 30,519 | **142 bp** |
| `BS-seq_10K_R1.fastq.gz` | 10,000 | 668 (6.68%) | 2,799 (27.99%) | **2,131 (21.31%)** | 21,652 | 61 bp |
| `SRR24766921_RRBS_R1.fastq.gz` | 10,000 | 258 (2.58%) | 1,325 (13.25%) | **1,067 (10.67%)** | 10,054 | 33 bp |

### The plan's own disqualifying metric, applied to the plan's own candidate

The plan rejects reverting to `-a "A{20}"` because, on `illumina10K_with_polyA.fastq.gz`, that path
"trimmed more than the current scan on 3,009 reads, 1,965 of them where the current scan trims nothing",
and because "the long-cut class is a real data-loss risk on bisulfite libraries". I measured both paths
side by side on that fixture:

| | **candidate B2** | **`A{20}` (rejected)** |
|---|---|---|
| trims more than current | 2,510 | 3,886 |
| of which 0 → >0 | **2,510** | 2,731 |
| of which >30 bp extra | **49** | 478 |
| extra bp total | 30,519 | 53,770 |

B2 is ~8% better than the alternative the plan calls disqualifying, on the plan's own headline metric.
It still has 49 reads losing >30 bp and a worst case of **142 bp — worse than the 117 bp the plan cites
as the reason `A{20}` is unacceptable.** The plan's §"Decision" argument therefore refutes B2 as well as
`A{20}`; the plan does not notice because it never ran the candidate over the fixture it used to build
the argument.

(Caveat on reproducibility: I could not reproduce the plan's exact `A{20}` figures. I get 3,886/2,731 on
`illumina10K_with_polyA` and 2,986/2,731 on `illumina_10K`; the plan reports 3,009/1,965 over "10,898
reads". The shapes agree, the digits do not. Whatever harness produced the plan's numbers no longer
exists, so the discrepancy cannot be resolved — flagging it because the plan's central decision rests on
those digits.)

### Mechanism — B2 is not "step over junk", it is a global error-gate relaxation

This is the part the plan misunderstands, and it is why the collateral is so large. Here is the
worst-case read, per-offset:

```
read = TACCTAAC AAAAAAAA…(A-rich 130 bp)…CGGAAGAGCGTTTACCCAGGAATCACGAAACCGCC   (150 bp)

  offset  0: sublen 150  idx 150        dropped=''
  offset  1: sublen 149  idx 149        dropped='C'
  offset  2: sublen 148  idx   8  CUT   dropped='CC'      <-- 140 bp removed
```

Dropping **two** terminal `C`s relocates the cut point 140 bases upstream. The stepped-over bases have
nothing to do with where the cut lands. The gate inside `homopolymer_trim_index` is
`errors * 5 <= n - i`. Removing one non-A base from the window shrinks the left side by **5** and the
right side by **1** — a net 4-point relaxation of the 20% error budget *for the entire scan*, not just
near the 3' terminus. Eleven offsets therefore means eleven progressively laxer scans, and
"first offset that cuts wins" means the loop fires on whichever scan the budget happens to admit first.

The firing-offset histogram is the fingerprint: `0:1537 1:366 2:315 3:286 4:255 5:263 6:225 7:220
8:204 9:195 10:181`. If B2 were genuinely tolerating 3' artifacts of a given length, offsets would
concentrate at that length. A near-flat tail across 1–10 is what a budget relaxation looks like.

**Action (Critical): drop B2 as specified.** If a 3'-tolerance margin is still wanted, use §5.

---

## 3. CRITICAL — the offset loop's "did it cut?" predicate is undefined, and the wrong choice silently truncates every read

B2 says "The first offset that yields a cut wins" but never defines *yields a cut*. There are two
readings and the plan's Signature section does not disambiguate:

```rust
if idx < sub.len()      { return idx; }   // correct
if idx < sequence.len() { return idx; }   // silently wrong
```

`homopolymer_trim_index` returns `n` (its own argument's length) on "no tail". At `offset > 0` that
return value is `sequence.len() - offset`, which is `< sequence.len()`. Under the second reading,
**every read with no poly-A tail at offset 0 gets truncated by `offset` bases**. I implemented both; the
`buggy` column of my probe shows the signature — every single case comes back exactly one base short:

```
case                        len  current  candidate  buggy
no tail (insert60)           60       60         60     59
insert60 + A20 + G20        100      100        100     99
internal A16 + insert16      92       92         92     91
```

Worse, at `offset == n` the slice is empty, `homopolymer_trim_index` returns `0`, and the buggy
predicate reads that as "cut at 0" → `record.truncate(0)` → **the read is emptied**. The plan's edge
table says "offset loop is bounded by `n`", which guards the slice but not this.

Two further implementation hazards the plan does not state:

- **Slice panic.** `&sequence[..n - offset]` underflows for `n < offset`. The edge table asserts the
  loop is bounded by `n`, but the Signature and Implementation-outline sections — which is what an
  implementer works from — do not carry the bound. Needs `if offset > n { break; }` or
  `0..=MAX_TAIL_OFFSET.min(n)` written into the plan.
- **No index remapping needed** — `sub` is a prefix, so `idx` is already an index into the original
  sequence, and `record.truncate(idx)` removes the stepped-over bases as intended. The plan is right
  here, but it is right by accident of slicing rather than by stated reasoning. Worth one line, because
  a reader who assumes remapping is needed will add a wrong `+ offset`.

**Action (Critical):** if any offset loop survives review, the plan must state the predicate verbatim,
the `offset > n` guard, and the no-remapping property. These are not implementation details — the wrong
predicate is a whole-library silent truncation that no proposed test would catch (see §6).

---

## 4. IMPORTANT — B1 alone fixes every reported case; the plan undersells it

The plan says `A20 + G20` "still needs B1", and frames B1+B2 as jointly necessary. I measured B1 alone
(poly-G first, then the *current* offset-0 poly-A scan, no B2):

| case | len | current (A→G) | **B1 only (G→A)** |
|---|---|---|---|
| `insert60 + A40` | 100 | 60 | 60 |
| `insert60 + A40 + G15` | 115 | 100 | **60** ✅ |
| `insert40 + A60 + G20` | 120 | 100 | **40** ✅ |
| `insert60 + A20 + G20` | 100 | 80 | **60** ✅ |
| `insert60 + A20 + G5` | 85 | 60 | 60 |
| `insert60 + A20 + junk10` | 90 | 90 | 90 |
| internal 16-A | 92 | 92 | 92 |

**B1 alone fixes all three poly-G cases**, with zero collateral on real data (it cannot change any read
where `poly_a` is off, and it is a provable no-op on R2 per §1). The only case it does not fix is
`insert60 + A20 + junk10` — non-homopolymer junk, which came from a synthetic probe, **not from the
reporter**. The reporter's symptom ("adapter and polyg are successfully removed, but the polya is only
partially trimmed") is the G case.

So the plan's own risk/benefit inverts: B1 carries essentially no risk and fixes the issue; B2 carries a
21–25%-of-reads regression and fixes only a case nobody reported. The plan calls B2 "the robustness
margin" — on this evidence it is the opposite.

**Action (Important):** restructure as "B1 is the fix, shipped alone", and move any 3'-tolerance work to
a separate, separately-validated change.

---

## 5. IMPORTANT — the third option the plan did not consider, and it is better than B2

The plan weighs exactly two algorithms (current scan vs `A{20}` DP). The review brief asks about
"making the poly-A scan skip a homopolymer run of any base" — I implemented and measured it:

> Before the 3' poly-A scan, if there is no tail at offset 0, measure the terminal run of whatever
> single non-A base ends the read. If that run is ≥ `min_run`, re-scan `seq[..n-run]`. Bounded by the
> actual run, not by a constant.

| case | len | current | **alt-C (min_run=3)** | B2 |
|---|---|---|---|---|
| `insert60 + A40 + G15` | 115 | 115 | **60** ✅ | 60 |
| `insert40 + A60 + G20` | 120 | 120 | **40** ✅ | 40 |
| `insert60 + A20 + G20` | 100 | 100 | **60** ✅ | 100 ❌ |
| `insert60 + A20 + junk10` | 90 | 90 | 90 | 60 |
| internal 16-A | 92 | 92 | **92** ✅ | 92 |
| **worst real read (150 bp)** | 150 | 150 | **150** ✅ | **8** ❌ |

Real-data collateral, same metric as §2:

| | reads trims-more | newly trimmed | >30 bp extra | worst | extra bp |
|---|---|---|---|---|---|
| **B2** on `illumina_10K` | 2,510 | 2,510 (25.10%) | 49 | 142 bp | 30,519 |
| **alt-C min_run=3** | 82 | 82 (0.82%) | 1 | 44 bp | 890 |
| **alt-C min_run=5** | 8 | 8 (0.08%) | 0 | 19 bp | 91 |
| **alt-C min_run=3** on `BS-seq_10K_R1` | 134 | 134 (1.34%) | 0 | 30 bp | 1,286 |
| **alt-C min_run=5** on `BS-seq_10K_R1` | 50 | 50 (0.50%) | 0 | 30 bp | 598 |

Alt-C is strictly better than B2 on every axis I measured: it fixes `A20+G20` (which B2 does not), it
leaves the 142 bp read alone, and it costs 0.08–1.34% of reads instead of 21–25%. It is also
*self-limiting* by construction — the skip length is the observed run, so it cannot relax the gate by
an arbitrary amount on a read with no artifact.

Two honest caveats: alt-C does not help the non-homopolymer `junk10` case (nor should it, in my view —
that is what quality trimming and `--2colour` are for), and because B1 already fixes every G case
end-to-end, **alt-C is probably unnecessary too**. I report it because it is the right shape for the
problem if a margin is wanted, and because it shows B2's design is not the only way to get tolerance.

Two further alternatives, dismissed with reasons:

- **Run the poly-A stage twice** (A → G → A). Equivalent to B1 for the G case, strictly more work, and
  the second pass re-opens the question of which pass owns `poly_a_trimmed`. No advantage over B1.
- **Handle it in quality trimming.** Already available as `--2colour`, already the reporter's
  workaround, and already verified sound (§1). Worth documenting, not worth re-implementing. Note it is
  not universal: it only works while the G-run is at the very 3' end — if adapter remnants sit behind
  the G-run, `--2colour` will not reach it. One sentence of hedging in `poly.md` would be honest.

---

## 6. IMPORTANT — every proposed assertion is one-sided; none would catch B2's regression

The proposed guard is *"assert no output read ends in ≥10 A"* (`tests/integration_poly.rs`), repeated
for the `--2colour` case. That is an **under-trim detector only**. Over-trimming makes reads shorter,
which satisfies "does not end in A" more easily. Specifically:

- It would **not** catch B2's 25%-of-reads regression.
- It would **not** catch the §3 predicate bug (a read truncated by 1–10 bases still does not end in A).
- It would **not** catch `truncate(0)` emptying reads — a zero-length read trivially passes, and may be
  dropped by the length filter before the assertion ever sees it.

Validation 3's expectation — *"no read loses >tail length"* — is not a measurable postcondition, because
B2's whole purpose is to redefine what counts as the tail. As written it can be declared satisfied by
inspection.

**Action (Important): add a two-sided guard.** Concretely: run `--poly_a` over
`test_files/illumina_10K.fastq.gz` (a library with no biological poly-A) before and after, and assert a
ceiling on `Reads with poly-A/T trimmed` and on `Poly-A/T bases removed` from the trimming report. On
today's code those are 1,537 reads / ~? bp; any change that moves the read count by thousands fails
loudly. This is the single assertion that would have caught B2, and it is cheap — the fixture and the
report fields already exist (`src/report.rs:401-411`).

Other validation gaps:

- **No test pins the stage order where it lives.** Commit 1 step 4 *renames* `test_poly_g_after_poly_a_removal`
  → `test_poly_g_tail`, which deletes the only in-tree record of the ordering. `trimmer.rs` has a
  ready-made `poly_g_config()` helper (`src/trimmer.rs:1111`) that does not set `poly_a`; the plan should
  add a `trim_read`-level test with **both** flags on asserting `[insert60][A40][G15] → 60`. Unit tests in
  `quality.rs` cannot express a stage order — only `trim_read` can.
- **Commit 3 item 1 bakes in a known deficiency**: "tail behind 20 G (expect *not* trimmed at this level
  — B1 covers it)". That test asserts a bug. If B2 is dropped it is wrong; if B2 is ever improved it
  blocks the improvement. Delete it, or state the postcondition at the `trim_read` level where it is
  actually true.
- **Validation 4 is right to be called the one that matters**, and §1 confirms its premise holds. Good.

Would the proposed tests have caught the original bug? **Yes** — `tests/integration_poly.rs` with
`[insert][polyA][polyG]` fails on v2.3.0. That part of the plan is well targeted, and the "no test under
`tests/` exercises `--poly_a` at all" diagnosis is correct and is the real root cause.

---

## 7. OPTIONAL — accuracy and documentation findings

- **Context table arithmetic.** Row 3, `insert40 + A60 + G20`, is listed as 140 bp. 40+60+20 = **120**,
  and the "emits 100" column is only consistent with 120. Row 4 (`insert60 + A20 + 10 bp junk` → "emits
  88") does not reproduce from the stated inputs either — the current scan returns 90 (no cut) on
  non-A junk, so the 2 bp must come from the junk's composition or from quality trimming. The probes
  that produced these numbers no longer exist, so neither row can be checked. Since the plan's whole
  Context section is the diagnosis, the numbers should be regenerated by a probe that is kept (e.g. as
  the `tests/integration_poly.rs` fixtures) rather than discarded.
- **Commit 4 misses a third wrong docs claim.** The plan catches `migration.md:28` and `:51-53`. It does
  not catch `docs/src/content/docs/modes/poly.md`, whose opening sentence reads *"v2 ships with two
  built-in poly-X trimmers **that the Perl version did not have**"*. That is the same factual error as
  the migration row, in the page commit 4 is already editing. Fix it in the same commit.
- **Perl `--polyA` was a separate run-and-exit mode, not a pipeline stage.** `0.6.10:trim_galore:1130`
  is an `elsif($polyA)` branch — a second Cutadapt pass over an already-trimmed file, documented at
  `:3375` as `trim_galore --polyA file_trimmed.fq.gz`, and it printed `POLY-A TRIMMING MODE;
  EXPERIMENTAL!!`. It also auto-detected poly-A vs poly-T orientation (`autodetect_polyA_type()`,
  `:2304`). The plan's "Perl `--polyA` was a Cutadapt invocation with `-a \"A{20}\"`" is true but flat;
  the migration-doc correction should say *experimental separate second pass*, or a migrating user will
  read "v0.6.x had --polyA" and expect the same integration semantics.
- **`MAX_TAIL_OFFSET` placement.** Per CLAUDE.md ("put an invariant on the thing it constrains"), the
  const and its one-line doc belong adjacent to the scan it bounds — the plan has this right.
- **Report line order** (`src/report.rs:401-425` prints poly-A before poly-G). The plan leaves it alone
  "so the reporter's existing reports stay comparable". Agreed, and worth keeping that reasoning out of
  the source comment.
- **Q2 (no Perl `--polyA` validation step): concur, and the reasoning is stronger than stated.** Perl's
  `--polyA` is a separate mode over a pre-trimmed file, so there is no single command whose md5 could be
  compared against a v2 single-pass run. Byte-identity is not merely the wrong gate — it is not
  expressible.
- **Q1 (`--no_poly_g --poly_a`): concur**, leave it, document it.
- **Q3 (`MAX_TAIL_OFFSET = 10`): moot if B2 is dropped.** If it survives, note that the measurements in
  §2 show the constant is not the lever that matters — even `offset ≤ 2` relocates a cut by 140 bases,
  so lowering it does not make the design safe.
- **Assumption 1 is correctly flagged as the gating unknown**, and the plan is right to say re-diagnose
  if `Reads with poly-A/T trimmed` comes back high. Good discipline; keep it.

---

## 8. Action items

### Critical

1. **Drop B2 as specified.** Measured regression: 2,131–2,510 reads newly trimmed per 10,000 on
   `illumina_10K` / `BS-seq_10K_R1` / `SRR24766921_RRBS_R1`, 49 reads losing >30 bp, worst single cut
   142 bp on `illumina10K_with_polyA.fastq.gz` — the same fixture and the same failure class the plan
   uses to reject `A{20}`. Rewrite the "Decision: keep the Cutadapt `--poly-a` semantics" section to
   account for this, since its argument disqualifies B2 too.
2. **If any offset loop survives, specify the cut predicate, the `offset > n` guard, and the
   no-remapping property verbatim.** `idx < sequence.len()` instead of `idx < sub.len()` truncates every
   tail-free read, and at `offset == n` empties it. No proposed test would catch either (§6).

### Important

3. **Ship B1 alone as the fix for #504.** Measured: fixes `A40+G15` → 60, `A60+G20` → 40,
   `A20+G20` → 60, all three with zero collateral. R2 no-op confirmed over 33,012 reads plus eight
   adversarial cases.
4. **Add a two-sided over-trim guard** on `--poly_a` over `test_files/illumina_10K.fastq.gz`, asserting
   a ceiling on `Reads with poly-A/T trimmed` and `Poly-A/T bases removed` from the trimming report.
   Every assertion the plan proposes detects under-trimming only.
5. **Add a `trim_read`-level ordering test** with both `poly_a` and `poly_g` set (extend the existing
   `poly_g_config()` at `src/trimmer.rs:1111`). Renaming `test_poly_g_after_poly_a_removal` removes the
   only in-tree record of the order without replacing it.
6. **Drop the "tail behind 20 G, expect not trimmed" unit test** — it asserts a known deficiency at a
   level where the postcondition is false after B1.
7. **Fix the third wrong docs claim**: `docs/src/content/docs/modes/poly.md`, "two built-in poly-X
   trimmers that the Perl version did not have".

### Optional

8. **If a 3'-tolerance margin is still wanted, use the terminal-homopolymer-run skip (§5)** rather than a
   blind offset loop: fixes one case more than B2, leaves the 142 bp read alone, and costs 0.08–1.34% of
   reads instead of 21–25%. Prefer `min_run = 5`.
9. Correct the Context table's row-3 arithmetic (140 → 120) and regenerate row 4, or keep the probe that
   produced them as the integration fixture.
10. In the migration doc, describe Perl `--polyA` as an experimental separate second pass with
    orientation auto-detection, not as an inline stage.
11. Hedge the `--2colour` recommendation in `poly.md`: it only reaches G-runs that are flush with the
    3' end.

---

## 9. Scope note

Temporary probes were written to `tests/zz_review_a_probe{,2,3,4}.rs` and **deleted** (`tests/` is clean
of them; `tests/zz_review_b_probe.rs` belongs to Reviewer B and was not touched). No file under `src/`,
`docs/` or `.github/` was modified. The only build artifacts are release test binaries under `target/`.
