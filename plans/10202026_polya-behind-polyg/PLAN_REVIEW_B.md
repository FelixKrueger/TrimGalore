# Plan Review B — poly-A tails survive behind a poly-G tail (#504)

**Plan:** `plans/10202026_polya-behind-polyg/PLAN.md`
**Reviewer:** B (independent, fresh context)
**Date:** 2026-10-20
**Tree state:** branch `wiki/init` @ `d471269`, `src/` at `dev` state. No source edits made.

**Method.** Every load-bearing claim was re-derived from the source rather than taken from the
plan. Where the plan's probes no longer exist I wrote my own (`tests/zz_review_b_probe.rs`, since
deleted) and measured against the real fixtures. Numbers below marked **[measured]** come from that
probe; **[read]** means verified by reading the code; **[inferred]** means reasoning, not measurement.

---

## Verdict

**B1 (stage reorder) is correct, well-scoped and low-risk — ship it.** Verified independently:
it fixes the reported case and changes 12 of 10,900 reads on `illumina10K_with_polyA.fastq.gz`
and 16 of 10,000 on `BS-seq_10K_R1.fastq.gz`. **[measured]**

**B2 (bounded 3' offset) is not safe as specified and should not ship in this form.** The plan's
central safety claim — "the internal-A case is unchanged, so B2 does not re-open the v0.6.x
over-trim" — is false. Measured end-to-end (quality Q20 → Illumina adapter → poly stages),
B2 alone changes **3,460 of 10,900 reads (31.7%)** on the poly-A fixture and **2,253 of 10,000
(22.5%)** on the BS-seq fixture. The entire blast radius of the plan lives in B2, not B1. **[measured]**

The plan reached the opposite conclusion from a single synthetic read. The mechanism is not subtle,
and I think it was simply missed: the stepped-over bases are **excluded from the slice**, so they
are not charged against the 20% error budget. That is what makes the offset loop work on
`A40+G15` — and it is the same thing that lets a bare `AAA` run sitting ≤10 bases from the 3' end
authorise the deletion of those 10 bases of real sequence.

---

## 1. Stage order and blast radius — verified, with one gap

| Plan claim | Status |
|---|---|
| `trimmer::trim_read` is the only implementation of the ordering | **Confirmed** [read] |
| `parallel.rs` holds no copy, one edit site | **Confirmed** [read] |
| `bp_after_cutadapt` captured ahead of both poly stages → #232 untouched | **Confirmed** [read] |
| `--rrbs` unaffected (step 2.5 precedes both) | **Confirmed** [read] |
| R2 output unchanged by the reorder | **Confirmed, and now measured** [measured] |
| No specialty mode / `--clump_only` path touched | **Confirmed** [read] |
| uBAM output unaffected as a separate path | **Confirmed** [read] |

Details worth recording:

- All six production call sites (`src/trimmer.rs:319`, `:479-480`, `:645`, `:730-731`;
  `src/parallel.rs:495-496`, `:1133`) funnel through the single `trim_read`. The uBAM-output
  entry points `run_single_end_to_bam` (`src/trimmer.rs:634`) and `run_paired_end_to_bam` (`:704`)
  are among them, so uBAM output inherits the change with no second edit — and `--poly_a` is *not*
  in `cli.rs::validate`'s uBAM rejection list (only `--clump_only` rejects it, `src/cli.rs:863-871`),
  so `--output-format ubam --poly_a` is a reachable combination that the change reaches. Correct
  as the plan states, but the plan should say so explicitly rather than leaving uBAM unmentioned.
- `specialty.rs` and `clump*.rs` contain no `poly_a` / `poly_g` / `homopolymer` references at all. [read]
- **R2:** the plan argues disjointness ("5' vs 3'"). The argument is incomplete — poly-G truncating
  the 3' end does shorten the sequence the 5' scan then sees, and the revcomp gate
  (`errors * 5 <= i + 1`, `src/quality.rs:135`) is prefix-local, so the two stages genuinely do
  commute. I verified empirically rather than relying on that: **0 of 13,000 reads** in
  `polyAT_R2.fastq.gz` change output under either B1 or B1+B2 with `is_r2 = true`. **[measured]**

### Gap: the reorder regresses the mirror arrangement `[insert][polyG][polyA]`

Not mentioned anywhere in the plan. On `insert60 + G15 + A40`: **[measured]**

```
poly-A first (today) -> 60 bp   (poly-A cuts the A40, poly-G then cuts the now-terminal G15)
poly-G first (B1)    -> 75 bp   (poly-G is blocked by the A40; poly-A cuts; the G15 survives)
```

This is the same blocking bug, mirrored, and B1 moves it from the A-stage to the G-stage. The plan's
Assumption 3 ("dark-cycle poly-G is strictly terminal") is the defence, and for sequencer-induced
G-runs it holds. It does not hold for template-derived G-runs (G-quadruplex, genomic poly-G ahead of
a genuine poly-A tail). The plan should either state this as an accepted consequence under
Assumption 3, or adopt the fixed-point alternative in §6, which removes the trade entirely.
In the fixtures this costs 5 reads where output gets *longer* after the change. **[measured]**

---

## 2. "No pinned CI output can move" — verified, claim holds

Checked exhaustively:

- `.github/workflows/*.yml`: a case-insensitive search for `poly` returns exactly two hits,
  `.github/workflows/ci.yml:967-968`, both `--no_poly_g` in the demux validation step. No
  `--polyA` / `--poly_a` / `--poly-a` anywhere in CI, including `validation` and
  `validation-ubam`. **[read]**
- Every `poly_a:` construction site: `src/main.rs:1294`, `:1528`, `:2318`, `:2697` all take
  `cli.poly_a`; `src/trimmer.rs:866`, `src/parallel.rs:1446`, `src/report.rs:1370` are test
  helpers hard-coded to `false`. The only place any poly flag is forced true in a test is
  `src/trimmer.rs:1113` (`c.poly_g = true`, the #321 R2 poly-G test) — poly-G, not poly-A. **[read]**
- `--library` presets resolve to clip flags only (`src/cli.rs:1233-1253`, `src/library.rs`); no
  preset can switch `poly_a` on. I checked this specifically because a preset would have been the
  one plausible counterexample. **[read]**
- No test under `tests/` passes `--poly_a` in any spelling. **[read]**

Both poly blocks are independently gated (`if config.poly_a && …` at `src/trimmer.rs:197`,
`if config.poly_g && …` at `:217`), so with `poly_a == false` the reorder is a no-op on both the
sequence and the stats. The claim is sound. **This is the one part of the plan that needed
verification most and survives it intact.**

---

## 3. Algorithm claims — re-derived; three numbers in the plan are wrong

Re-derived by reading `src/quality.rs::homopolymer_trim_index` (the 3' branch is `:149-169`) and by
probe. The scoring is +1/−2 with gate `errors * 5 <= n - i` and a strict `score > best_score`, so a
cut is only recorded where the running score is **positive** *and* the error fraction over the
window is ≤ 20%. With the A-run at index `p`, a run of `a` A's and `g` non-A bases behind it, a cut
at `p` requires `a - 2g > 0` **and** `4g <= a`. The second condition binds. **[read]**

Candidate rows, re-measured: **[measured]**

| Read | plan: current | mine | plan: candidate | mine |
|---|---|---|---|---|
| `insert60 + A40 + G15` | 115 | 115 | 60 | 60 ✓ |
| `insert40 + A60 + G20` | 120 | 120 | 40 | 40 ✓ |
| `insert60 + A20 + G20` | 100 | 100 | 100 | 100 ✓ |
| `insert60 + A20 + junk11` | 91 | 91 | 60 | 60 ✓ |
| `insert60 + A20 + junk10` | **90** | **60** ✗ | 60 | 60 |
| `junk5 + T40` (5') | 45 | 45 | — | 45 ✓ |
| `junk11 + T40` (5') | 51 | 51 | — | 51 ✓ |

The load-bearing rows reproduce. Three corrections:

1. **`insert60 + A20 + junk10`: the plan's "current = 90" is not a property of the code.** With
   random-ACGT junk the *current* code already cuts this read (60), because ~25% of the junk is A
   and the error budget covers the rest. Whether this row is a bug at all depends entirely on the
   junk's base composition. The row is presented as evidence that B2 is needed; it isn't.
2. **The first table's `insert40 + A60 + G20` is listed as 140 bp.** It is 120 bp. The same read is
   120 in the second table. One of the two is a transcription error.
3. **The first table mixes two measurement levels.** Column 3 is binary output length
   (quality+adapter+poly), the second table's columns are raw `poly_a_trim_index` return values.
   The `junk10` read appears as 88 in one and 90 in the other. Label the levels or a reader will
   take them as contradictory.

### 3a. The Perl comparison is misattributed (affects Commit 4, which writes it into the docs)

The plan says Perl `--polyA` was "a Cutadapt invocation with `-a "A{20}"`". Reading
`git show 0.6.10:trim_galore`: **[read]**

- `autodetect_polyA_type()` (`:2304`) counts reads containing `AAAAAAAAAA` vs `TTTTTTTTTT` in the
  first 1M reads and returns `$adapters{'PolyA'}->{seq}`, which is **`'AAAAAAAAAA'` — 10 A's**.
  `extend_adapter_sequence('A', 20)` is reached **only** in the inconclusive branch where *both*
  counts are zero. So on any library that actually has poly-A — the reporter's included — Perl ran
  `-a AAAAAAAAAA`, not `-a A{20}`. The `A{20}` figure appears only in the `--help` text
  (`:3358`) and in a validation error message (`:2750`); the help text and the code disagree.
- Measured with `find_3prime_adapter(.., 0.1, 1)` on the poly-A fixture after quality+adapter
  trimming: `A{10}` cuts **6,401 reads / 460,805 bp**; `A{20}` cuts **3,818 reads / 128,130 bp** —
  3.6× the bases. **[measured]** So real v0.6.x `--polyA` was considerably more destructive than the
  plan's stand-in. This *strengthens* the "Decision: keep the Cutadapt semantics" section; it also
  means the figures quoted there (117 bp cut, 3,009 reads) understate the case they are making.
- **Perl's `--polyA` was an exclusive mode, not a stage.** `elsif($polyA)` at `:1130` sits in the
  same `if ($rrbs) / elsif(...)` dispatch chain as `--hardtrim5` and `--clock`, and its Cutadapt
  invocation carries no `-q` and no `--length`. The `--help` text says so outright: "The poly-A
  trimming mode expects that sequences were both adapter and quality trimmed before looking for
  Poly-A tails". Issue #504 confirms the reporter did exactly that — two runs,
  `--quality 20 --phred33` then `--polyA`. The plan never says this, and it matters: the user's
  v0.6.10 and v2.3.0 command lines are not equivalent workflows, and `migration.md` should say
  "separate experimental second pass" rather than implying a renamed stage.
- Perl's `--polyA` also **rewrote read IDs** (`32:A:` prefix, `_PolyA:32` suffix, whitespace → `_`;
  `:3362-3368`). v2 has no equivalent — `grep -rn "PolyA:" src/` is empty. A 3'seq pipeline that
  greps read IDs for `PolyA:` silently gets nothing under v2. This belongs in the migration note
  Commit 4 is already rewriting. **[read]**

---

## 4. The candidate's correctness and edge cases

### 4a. CRITICAL — B2 deletes real 3' sequence at scale

Implementing B2 exactly as specified (`for k in 0..=MAX_TAIL_OFFSET.min(n)`, first offset whose
scan returns `idx < slice.len()` wins) and sweeping the fixtures: **[measured]**

Function level, `illumina10K_with_polyA.fastq.gz` (10,900 reads):

```
current trims                          : 2437
candidate trims                        : 4947
candidate trims MORE than current      : 2510  (23.03%)
  ...of which the read does NOT end in AAA : 2510   <-- all of them
worst extra: 142 bp | len 150, current keeps 150, candidate keeps 8
```

Pipeline level (Q20 quality trim → `-a AGATCGGAAGAGC` → poly stages), which is what a user sees:

| fixture | reads | B1 only | **B2 only** | B1+B2 | bases lost | reads losing >13 bp | dropped by `--length 20` |
|---|---|---|---|---|---|---|---|
| `illumina10K_with_polyA` | 10,900 | **12** | **3,460** | 3,475 | 40,341 | 952 | 90 → **131** |
| `BS-seq_10K_R1` | 10,000 | **16** | **2,253** | 2,277 | 22,721 | 348 | 1 → **4** |

Mechanism, in one line: the stepped-over bases are removed from the slice, so they are never
charged as errors. A read ending `…AAACCTTTAACAAAATTCAAC` has an `AAAA` run 10 bases from the 3'
end; offset 10 presents the scan with a clean A-tail, it cuts, and 10 bases of genuine sequence go
with it. Minimal constructed cases: **[measured]**

```
insert60 + A3  + junk10  ->  current 73, candidate 60   (3 A's authorise deleting 10 real bases)
insert60 + A16 + junk11  ->  current 87, candidate 60   (the plan's "internal A16 unchanged" case,
                                                         with the run 11 bases from the end)
```

The plan's `internal 16-A run | 76 | 76 | **no over-trim**` row is true only because that particular
read's A-run sits >10 bases from the 3' end. Move it closer and the row inverts. One read is not
evidence for a safety property.

For scale: the plan rejects reverting to `A{20}` because it "trimmed more than the current scan on
3,009 reads, 1,965 of them where the current scan trims nothing". B2 trims more on 2,510 reads,
**2,510** of them where the read does not even end in A. B2 re-opens the exact failure mode the
plan's own Decision section rules out, at the same order of magnitude, on a bisulfite fixture —
the library class the plan names as the data-loss risk.

**Action:** do not ship B2 as specified. Options in §6.

### 4b. CRITICAL — the implementation is under-specified in two ways that are silent-wrong, not loud-wrong

The Signature block gives only the doc comment; Commit 2 says "rewrite `poly_a_trim_index`'s body as
the offset loop". Two details decide whether the result is correct, and neither is stated:

1. **The no-cut sentinel must be compared against `slice.len()`, not `n`.**
   `homopolymer_trim_index` returns the length of what it was *given*. The natural-looking
   `let idx = homopolymer_trim_index(&seq[..n - k], b'A', false); if idx != n { return idx; }`
   returns `n - 1` at `k = 1` for **every read with no poly-A tail** — a silent universal 1 bp 3'
   trim for all `--poly_a` users, passing every test in the plan's Commit 3 list (all of which
   assert on reads that *do* have tails, plus `n < 3`). This is the single most likely way to
   implement the plan as written and get it wrong.
2. **The loop bound must be `MAX_TAIL_OFFSET.min(n)` (or `n.saturating_sub(k)`).**
   `&seq[..n - k]` with `n < k` underflows; in release that becomes a huge `usize` and the slice
   panics. Reads below 11 bp reach `trim_read` routinely — length filtering runs *after* it
   (`src/filters.rs`), so `--poly_a` on short-insert data would panic. The plan mentions the bound
   in the edge-case table prose but not in the Signature or Commit 2, and the test list jumps from
   `n < 3` straight to full-length reads, leaving `3 ≤ n ≤ 10` — the exact panic window — untested.

Both are cheap to fix in the plan: give the body, and add one test at `n = 5`.

### 4c. `MAX_TAIL_OFFSET = 10` and "first offset that cuts"

- **The constant is not the problem; the un-charged skip is.** At any offset the 20% gate still
  applies to the slice, so raising or lowering 10 only scales the over-trim — it does not change
  its character. Q3's framing ("no data behind the exact value") understates this: the value is not
  a tolerance knob, it is the number of real bases the stage may discard without justification.
- **"First offset that cuts = most conservative" is true but weaker than it sounds.** The cut index
  is determined by where the A-run starts, which does not move with `k`; what moves is whether a cut
  happens at all. So the selection rule is effectively "cut at the smallest `k` that makes a cut
  legal", and the alternative framings (earliest cut, largest `k`) differ only in pathological
  reads. The rule is fine. It is not what makes B2 unsafe.
- No panic, no `revcomp` mishandling (the revcomp branch delegates unchanged), and
  `poly_a_trimmed` is computed at the call site as `seq_len_before - idx`
  (`src/trimmer.rs:209-211`), so it cannot double-count — the plan's accounting note is correct. [read]

### 4d. IMPORTANT — B3 is not established, and B2 *creates* the asymmetry it claims to remove

The plan's two 5' probes reproduce exactly (`junk5 + T40 → 45`, `junk11 + T40 → 51`) **[measured]**,
but they do not support the conclusion. Those "junk" bases are random ACGT, so ~25% of them are T
and are scored +1, not as errors. With a *homopolymer* prefix the 5' branch fails identically to
the 3' branch: **[measured]**

```
G15 + T40 + insert60  ->  revcomp clip = 0     (5' mirror of #504; nothing is clipped)
```

Algebraically the two branches have **identical** tolerance before B2: 3' needs `4g <= a`, 5' needs
`4j <= t`. **[read]** So:

- "the 5' branch is already tolerant" is true only for mixed junk, not as a property of the branch;
- B2 gives the 3' branch a free 10-base un-charged skip that the 5' branch does not get, which
  **introduces** the asymmetry the plan says it removes.

Practical relevance: on R2 of a 3'seq library the poly-T head is the payload end, and anything 5' of
it — a UMI, a barcode, residual primer — blocks the scan. `--clip_r2` would remove it, but fixed
clipping is step 4, *after* the poly stages, so it cannot help. Worth a sentence in B3 and a
follow-up issue; not worth expanding this plan's scope.

### 4e. Edge cases the plan's table misses

| Case | Status |
|---|---|
| `3 ≤ n ≤ MAX_TAIL_OFFSET` | **missing** — the panic window (§4b) |
| Short tail *behind* an offset (`A3 + junk10`) | **missing** — the table's "tail < 3bp ignored at every offset" is wrong: a 3bp tail is not < 3, and at offset 10 it cuts 13 bases |
| Read with no poly-A at all under B2 | **missing** — the fixture sweep shows 23% of such reads get cut |
| `[insert][polyG][polyA]` under B1 | **missing** (§1) |
| Homopolymer prefix on the 5' branch | **missing** (§4d) |
| `n < 3`, empty, all-A, R2, `--rrbs` | present and correct [read] |

---

## 5. Validation sufficiency

**Would the proposed tests have caught #504?** Yes. `tests/integration_poly.rs` asserting "no output
read ends in ≥10 A" over `[insert][polyA][polyG]` would have failed on v2.3.0. Good test, right
level, and the plan is right that this is the gap that let the bug ship.

**Would they catch a regression?** For the ordering, yes. For over-trimming, no — and that is the
failure mode this change actually introduces. Every assertion in Commit 3 is **one-sided**: it checks
that A's are *gone*, never that non-A sequence is *kept*. A build that truncated every read to 1 bp
would pass all of them.

**Validation #3 is not a usable gate.** "no read loses >tail length" is undefined for a candidate
whose stated purpose is to remove bases beyond the tail, and `>tail length` is unmeasurable without
first agreeing what the tail is. As written it would have passed B2 while B2 shortens 31.7% of reads.

Concretely missing:

1. A **no-poly-A control read** that must come out byte-identical. (The plan's "clean
   `[insert][polyA]` control" still has a tail.)
2. A **fixture differential sweep** as a committed test or a recorded pre-merge measurement:
   run `poly_a_trim_index` (or the full binary) over `illumina10K_with_polyA.fastq.gz` before and
   after, and assert bounds on `reads where new < old` and on the maximum extra bases. This is the
   only check that would have caught §4a, it needs no Perl and no conda, and it runs in 0.2 s.
   My probe is a working template; it is deleted, but the shape is: for each read, compare old and
   new kept-length, histogram the difference, and assert the tail of that histogram is empty.
3. A test at `3 ≤ n ≤ 10` (§4b).
4. `--paired --poly_a` R2 byte-identity as a *test*, not a manual step — Validation #5 is a manual
   before/after diff and will not survive the next refactor. (My sweep says it passes today:
   0/13,000 on `polyAT_R2.fastq.gz`.)
5. The `[insert][polyG][polyA]` mirror, asserted in whichever direction the plan decides is correct.

**Q2 (add a Perl `--polyA` step to `validation`): agree, no.** §3a strengthens the reasoning —
Perl's real adapter was 10 A's, so byte-identity is not merely the wrong gate, the behaviour being
compared is a different workflow entirely. But declining byte-identity is not a reason to decline a
*behavioural* regression gate; item 2 above is that gate.

---

## 6. Alternatives

The plan considers exactly one alternative (revert to `A{20}`) and dismisses it correctly — more
strongly than it knows, per §3a. Three it did not consider, in my order of preference:

**(A) Ship B1 alone. Recommended.** B1 fixes the reported symptom, costs 12 reads in 10,900, and
carries no over-trim risk. The reporter's workaround (`--2colour 20`) already covers the residual
case, and B2 is explicitly described in the plan as "the robustness margin", not the fix. Shipping
B1 now resolves #504; B2 can be designed properly in its own plan with the fixture sweep as its
acceptance gate. This also keeps the fix small enough to backport.

**(B) Replace the blind offset with a homopolymer skip.** Instead of stepping over up to 10
*arbitrary* bases, skip a maximal 3'-terminal run of a single non-A base (G/C/T/N), of any length,
then scan. Properties: it cannot delete arbitrary sequence, because what it skips is provably a
homopolymer artefact; it fixes `A20 + G20`, which B2 cannot; it needs no magic constant; and it
makes the poly-A stage self-sufficient rather than order-dependent. The classes I measured in §4a
all arise from skipping *mixed* sequence and would disappear. Worth probing against the same
fixtures before adopting — I did not measure it.

**(C) Iterate the two poly stages to a fixed point** (alternate poly-G / poly-A until neither cuts,
bounded at two rounds). Order-independent, so it subsumes B1 *and* removes the `[polyG][polyA]`
regression in §1 instead of trading one blocking case for another. Cost is at most two extra O(n)
scans on a flag that is off by default. Termination is trivial — each round only truncates. This is
strictly better than B1 on correctness grounds; B1 is better on diff size. Either is defensible;
(C) + (B) together would be the complete fix.

**(D) Handle it in quality trimming** — i.e. what `--2colour` already does. Not general (template
G-runs carry real quality scores), and it conflates two orthogonal flags. Correctly left alone,
though the plan should say so since it is the workaround being recommended to the reporter.

---

## 7. Smaller notes

- **Stale line anchors.** The plan's implementation outline will send an implementer to the wrong
  place: `trim_read`'s doc comment is at **`src/trimmer.rs:79-88`** (plan says 95-101);
  `bp_after_cutadapt` is at **`:163`** (plan says 161); the poly-A block starts at **`:194`** and
  poly-G at **`:216`** (plan says ~192-212 / ~220-232, roughly right). `src/quality.rs:501`
  (`test_poly_g_after_poly_a_removal`) and `src/cli.rs:365-368` are correct. [read]
- **"No test exercises `--poly_a`" is true of `tests/` but not of the unit tests.** There are 11
  `poly_a_trim_index` tests at `src/quality.rs:347-427`. They all use tails flush with the 3' end,
  which is precisely why they missed #504 — a more useful framing for the "Why it was never caught"
  section than "no tests". I checked all 11 against the candidate by hand and by probe: none of them
  changes result under B2, so Commit 2 will not show up as a red test. That is worth knowing — it
  means the existing suite gives **no** signal on this change in either direction.
- **Report attribution shifts.** After B1, bases that are counted today under "Poly-G/C bases
  removed" move to "Poly-A/T bases removed" (`src/report.rs:401-425`). Harmless, but it is a
  user-visible report change for `--poly_a` users and the CHANGELOG sentence should not imply the
  reports stay comparable. The plan's note about report *line order* is about something else.
- **Assumption 1 is the right gate and should stay blocking.** Issue #504's body confirms the
  two-run v0.6.10 workflow and the symptom description, and Felix's reply already asks for the
  "Reads with poly-A/T trimmed" count. B1 is justified by the code regardless of what comes back,
  so the assumption blocks only the interpretation, not the fix.
- **Q4 / #245 item B.** Agree it should be re-opened, and §3a adds a reason: the fallback being
  recommended to users (`-a "A{20}"`) is both the thing #245 item B says diverges *and* not what
  v0.6.x actually ran.

---

## Action items

### Critical

1. **Do not ship B2 as specified.** It shortens 31.7% / 22.5% of reads on the two fixtures and
   removes 3' sequence from 2,510 reads that do not end in A at all, worst case 142 bp of 150.
   Replace the plan's single-read safety check with the fixture differential sweep, and adopt
   alternative (A) or (B) from §6. (§4a)
2. **If any offset loop ships, specify its body in the plan** — the no-cut sentinel compared against
   `slice.len()` and the `MAX_TAIL_OFFSET.min(n)` bound. As written, the two most natural
   implementations give a silent universal 1 bp trim and a panic on sub-11 bp reads respectively,
   and the proposed tests catch neither. (§4b)
3. **Make the Commit 3 assertions two-sided.** Add a no-poly-A control that must be byte-identical,
   and a test at `3 ≤ n ≤ 10`. Every assertion currently proposed would pass a build that truncates
   every read. (§5)

### Important

4. **Decide the `[insert][polyG][polyA]` mirror explicitly.** B1 moves the blocking bug from the
   A-stage to the G-stage; measured 60 bp → 75 bp on a constructed read, 5 reads per 10,900 in the
   fixtures. Either accept it under Assumption 3 and say so, or take alternative (C). (§1)
5. **Correct the Perl facts before Commit 4 writes them into `migration.md`:** the auto-detected
   adapter was `AAAAAAAAAA` (10 A's), `A{20}` was the zero-count fallback and the help text's claim;
   `--polyA` was an exclusive mode with no quality trimming or length filtering, run as a second
   pass; and it rewrote read IDs with `32:A:` / `_PolyA:32`, which v2 does not. (§3a)
6. **Fix B3.** The 5' branch is not inherently tolerant — `G15 + T40 + insert60` clips 0 — and B2
   gives the 3' branch an un-charged skip the 5' branch lacks, so it introduces an asymmetry rather
   than removing one. State the real boundary; file the 5'-prefix case separately. (§4d)
7. **Replace Validation #3 with a measurable criterion.** "no read loses >tail length" is undefined
   for this candidate and would have passed B2. (§5)
8. **Add `--paired --poly_a` R2 byte-identity as a committed test**, not a manual before/after. (§5)

### Optional

9. Fix the three numeric/labelling errors in the Context and Behavior tables: `insert40+A60+G20` is
   120 bp not 140; the `junk10` "current = 90" row depends on junk composition and is 60 for random
   ACGT; label which table is output length and which is `poly_a_trim_index`'s return. (§3)
10. Update the stale line anchors in the implementation outline. (§7)
11. Reframe "Why it was never caught": 11 unit tests do cover `poly_a_trim_index`, all with flush
    tails. Also note that none of them changes under B2, so the existing suite gives no signal
    either way. (§7)
12. Note the poly-A/poly-G report-attribution shift for `--poly_a` users. (§7)
13. Say in Integration that `--output-format ubam --poly_a` is reachable and inherits the change via
    the shared `trim_read`. (§1)
