# Plan: poly-A tails survive behind a poly-G tail (`--poly_a`)

**Issue:** [#504](https://github.com/FelixKrueger/TrimGalore/issues/504) — "Less effective PolyA
trimming with v2.3.0", reported by @marenmsk (3'seq libraries, v0.6.10 → v2.3.0).
**Date:** 2026-10-20 (rev 2, after dual review)
**Branch to cut:** `fix/polya-behind-polyg` off `dev`

## Revision history

**rev 2 (2026-10-20)** — after `PLAN_REVIEW_A.md` and `PLAN_REVIEW_B.md`. **B2 (bounded 3' offset)
is removed from the plan.** Both reviewers independently implemented it and swept the real fixtures:
it newly trims 21–25% of reads on ordinary libraries, removes 3' sequence from 2,510 reads that do
not end in `AAA` at all, and cuts 142 bp off a 150 bp read in `illumina10K_with_polyA.fastq.gz` —
the same fixture and failure class rev 1 used to disqualify the `A{20}` alternative. Rev 1's safety
argument for B2 rested on a single synthetic row and was wrong. Scope is now B1 alone, which both
reviewers verified is sufficient for every reported case. rev 1 is kept at `PLAN.md.v1`.

Also corrected in rev 2: the Perl `--polyA` facts (it used **10** A's, not 20, and was an exclusive
second pass), the claim that the 5' branch is inherently tolerant (it is not), two numeric errors,
and the one-sided validation criteria.

## Goal

`--poly_a` must remove poly-A tails that sit behind a 3' poly-G run. Today it does not, so 2-colour
3'seq libraries come out of v2.3.0 with full poly-A tails that v0.6.10 removed.

Non-goals, both deliberate:

- **Byte-identity with Perl `--polyA`.** See "Decision" below.
- **Tolerating arbitrary non-A bases at the 3' end.** Attempted in rev 1 and withdrawn; see
  "Withdrawn: the bounded 3' offset" and the successor option in "Follow-up".

## Context

### What changed between v0.6.x and v2

Perl `--polyA` ran Cutadapt with a poly-A homopolymer as the adapter. A 3' adapter match removes the
match **and everything downstream of it**, so `[insert][polyA][anything]` lost the tail and the
trailing bases together.

Three details of the Perl mode matter for the docs commit, and rev 1 got all three wrong
(`PLAN_REVIEW_B.md` §3a, re-read from `git show 0.6.10:trim_galore`):

- The adapter was **`A{20}`**, with `T{150}` as `-a2` on R2. `autodetect_polyA_type()` (`:2304`)
  returns a 10-mer, but the caller overwrites it two lines later with
  `extend_adapter_sequence("A", 20)` in **both** the paired and single-end branches (`:501-530`);
  the 10-mer only picks the direction via `$adapter =~ /A+/`, and its value never reaches Cutadapt.
  The v0.6.x `--help` text (`:3358`) was accurate. **Rev 1 had this right; rev 2 initially
  "corrected" it to 10 A's on a plan reviewer's reading and both code reviewers caught the
  inversion.**
- **It was an exclusive mode, not a stage.** `elsif($polyA)` (`:1130`) sits in the same dispatch
  chain as `--hardtrim5` and `--clock`, and its Cutadapt invocation carries no quality cutoff. The
  help text says so: "The poly-A trimming mode expects that sequences were both adapter and quality
  trimmed before looking for Poly-A tails". This is why the reporter ran two passes. Length
  filtering did still apply (`:1243`). Her v0.6.10 and v2.3.0 command lines are **not** equivalent
  workflows.
- **It rewrote read IDs** — `s/\s+/_/g` on every read, unconditionally, plus a `32:A:` prefix and
  `_PolyA:32` suffix on reads that lost bases (`:1206-1213`). v2 does neither
  (`grep -rn "PolyA:" src/` is empty).

v2's `--poly_a` is a different algorithm: `quality::homopolymer_trim_index` walks inward from the 3'
end (+1 per A, −2 per non-A) and records a cut only where the running score is positive **and** the
cumulative error gate `errors * 5 <= n - i` holds. With an A-run of length `a` and `g` non-A bases
behind it, a cut needs `a - 2g > 0` **and** `4g <= a`; the second binds.

`src/trimmer.rs::trim_read` runs poly-A at step 2.7 and poly-G at step 2.8 — **poly-A first**. On a
2-colour read `[insert][polyA][polyG]` the poly-A scan hits the G-run, fails the gate, and the
poly-G stage then strips the G-run, leaving the poly-A tail exposed in the output. That is the
reporter's symptom ("adapter and polyg are successfully removed, but the polya is only partially
trimmed").

### Reproduction

Measured three times independently — once building this plan, once by each reviewer. Figures below
are `poly_a_trim_index` return values (bases kept), not end-to-end output lengths:

| Read | input | current | after B1 |
|---|---|---|---|
| `insert60 + A40` | 100 | 60 | 60 |
| `insert60 + A40 + G15` | 115 | **115** (tail kept) | 60 |
| `insert40 + A60 + G20` | 120 | **120** (tail kept) | 40 |
| `insert60 + A20 + G20` | 100 | **100** (tail kept) | 60 |

`--2colour 20` in place of `--quality 20` works around it (G scores as `cutoff-1`, so the G-run goes
in step 1, ahead of everything) and is the workaround given to the reporter.

### Why it was never caught

- The `validation` job has **no `--polyA` step**. No Perl baseline was ever taken for this flag.
- No test under `tests/` passes `--poly_a`. `src/quality.rs` has 11 unit tests covering
  `poly_a_trim_index`, but **every one of them uses a tail flush with the 3' end** — the blocked
  case is untested at any level.
- `--poly_a` shipped as a *new* capability (CHANGELOG line 964), so nobody looked for a v0.6.x
  predecessor. `docs/.../migration.md:28` still records `| Poly-A trimming | None | Built-in
  --poly_a |`, and `docs/.../modes/poly.md:6` opens with "two built-in poly-X trimmers that the Perl
  version did not have". Both are wrong.
- Issue #245 item B was deferred partly on the reasoning that *"poly-A trimming has its own
  dedicated `--poly_a` flag that doesn't go through this path"*. That separation is the defect.

### Decision: keep the Cutadapt `--poly-a` semantics

v2's scan matches modern Cutadapt's `--poly-a`. It stays. Switching the stage to a homopolymer
adapter match would re-introduce v0.6.x's over-trimming: on
`test_files/illumina10K_with_polyA.fastq.gz` the `A{20}` match cut **117 bp** off one read because of
an internal 16-A stretch, where the current scan leaves it alone. The `A{10}` comparison that rev 2
briefly carried here is void — Perl never used a 10-mer as the adapter, so `A{20}` is the only
Perl-parity baseline.

### Withdrawn: the bounded 3' offset (rev 1's B2)

Recorded so it is not re-proposed. Rev 1 proposed letting the scan step over up to 10 bases at the
3' end. Both reviewers implemented it and swept the fixtures:

| Fixture | reads | currently trims | with the offset | **newly trimmed** | worst extra cut |
|---|---|---|---|---|---|
| `illumina_10K.fastq.gz` | 10,000 | 1,537 | 4,047 | **2,510 (25.1%)** | **142 bp** |
| `BS-seq_10K_R1.fastq.gz` | 10,000 | 668 | 2,799 | **2,131 (21.3%)** | 61 bp |
| `SRR24766921_RRBS_R1.fastq.gz` | 10,000 | 258 | 1,325 | **1,067 (10.7%)** | 33 bp |

In all 2,510 newly-trimmed reads the read does not end in `AAA` at all. The mechanism is not "step
over junk": excluding a terminal non-A base from the window drops `errors` by 1, so the gate's left
side falls by 5 while its right side falls by 1 — a **net relaxation of the global error budget**,
not a local skip. The firing-offset histogram is near-flat across 1–10, which is that signature.
Consequently a bare `AAA` run sitting ≤10 bases from the 3' end can authorise deleting 10 bases of
real sequence, and rev 1's reassuring `internal 16-A run | 76 | 76` row holds only because that
read's A-run sits >10 bases from the end; at 11 it inverts to 87 → 60.

## Behavior

### B1 — run poly-G before poly-A

`trim_read` runs poly-G then poly-A. All other stages keep their position: quality (1), adapter (2),
RRBS (2.5), N-trim (3), fixed clipping (4).

Rationale, one line in the source: poly-G is a terminal sequencer artifact, so nothing of value sits
behind it and a poly-A tail behind it is invisible to a 3'-anchored scan.

Scope, all verified by both reviewers:

- **R1 / single-end only.** Both stages compete for the 3' end on R1/SE. On R2, poly-T is removed
  from the 5' end and poly-G from the 3' end. Reviewer A measured 33,012 reads across four fixtures
  plus eight adversarial degenerate cases (`TTTGGG`, `T15G15`, `T20G3`, `T3G20`, …) and Reviewer B
  measured 13,000 on `polyAT_R2` — **zero divergence** in output length or either stat counter. The
  mechanism is exact: the 5' gate is `errors * 5 <= i + 1`, which depends on `i` only and never on
  `n`, so a 3' truncation cannot perturb it.
- **No existing output can move.** Both reviewers hunted for a counterexample and found none. Only
  two `poly` hits exist in all workflows (`ci.yml:967-968`, both `--no_poly_g`); `poly_a: true`
  appears nowhere in `src/`, `tests/` or `examples/`; all test-side `TrimConfig` sites hard-code
  `false`; there is no `Default` impl to leak a `true` through `..Default::default()`; and
  `src/library.rs`'s `--library` presets touch clip flags only. The sole `true` source is
  `cli.poly_a`, at four `TrimConfig` sites in `src/main.rs` (`:1294`, `:1528`, `:2318`, `:2697`).
  Both poly blocks are independently gated, so the reorder is a no-op whenever `poly_a` is false.
- **Collateral:** 12 reads in 10,900 on `illumina10K_with_polyA`, 16 in 10,000 on `BS-seq_10K_R1`
  (Reviewer B, end-to-end).

### Accepted consequence — the mirror arrangement `[insert][polyG][polyA]`

B1 moves the blocking case from the A-stage to the G-stage. Measured on `insert60 + G15 + A40`:
60 bp today, **75 bp** after (poly-G is blocked by the A-run, poly-A cuts, the G-run survives). Five
reads per 10,900 in the fixtures get *longer*.

This is accepted under Assumption 3: dark-cycle poly-G is terminal by construction, so
`[polyG][polyA]` is not a real 2-colour arrangement. It is not accepted as *correct* — it is a known
residual, it under-trims rather than deleting data, and the fixed-point option in "Follow-up"
removes it. A test pins the behaviour so the trade is visible rather than accidental.

### Edge cases

| Case | Handling |
|---|---|
| Read < 3 bp | existing `n < 3` guard in `homopolymer_trim_index` |
| Empty record | existing `!record.is_empty()` guards at both call sites |
| Read is all A, or all G | unchanged |
| Tail < 3 bp | unchanged — still ignored |
| `--no_poly_g --poly_a` on 2-colour data | tails behind G-runs survive. Explicit opt-out; documented, not worked around (Q1) |
| `[insert][polyG][polyA]` | accepted residual, above |
| `--poly_a` on R2 | unchanged (measured) |
| `--rrbs` | RRBS is step 2.5, ahead of both poly stages |
| `--output-format ubam --poly_a` | reachable, and inherits the change through the shared `trim_read` (`src/trimmer.rs:645`, `:730-731`) |
| `--clump_only --poly_a` | unreachable — rejected at `src/cli.rs:863-871` |

## Implementation outline

Three commits on `fix/polya-behind-polyg`. No change to `src/quality.rs`'s algorithm.

### Commit 1 — reorder the poly stages

1. `src/trimmer.rs`: move the poly-G block above the poly-A block; renumber the comment headers so
   poly-G is 2.7 and poly-A is 2.8. Keep the existing #321 explanation verbatim — it is still
   accurate and it is the only in-tree record of why poly-G is 3'-on-both-reads.
2. Add one line above the poly-G block giving the invariant (current state, no history, per
   CLAUDE.md "Code comments").
3. Update the `trim_read` doc comment, whose processing-order list omits the poly stages entirely.
4. `src/quality.rs`: the comment on `test_poly_g_after_poly_a_removal` states the old order in prose
   while its assertion is order-independent. Reword the comment; rename to `test_poly_g_tail`.
   **Do not leave the stage order recorded nowhere** — commit 2's `trim_read`-level test replaces it.

### Commit 2 — tests

The gap that let this ship. Note that none of the 11 existing `poly_a_trim_index` unit tests changes
under either B1 or the withdrawn B2, so the existing suite gives no signal in either direction.

1. `src/trimmer.rs` unit tests at `trim_read` level, using the existing `poly_g_config()` helper
   (`src/trimmer.rs:1111`) with `poly_a` also set: `[insert][polyA][polyG]` → tail removed;
   `[insert][polyG][polyA]` → pins the accepted residual; `--no_poly_g` → pins Q1.
2. `src/quality.rs` unit tests for the window the current list skips: `3 ≤ n ≤ 10`.
3. New `tests/integration_poly.rs` — synthetic FASTQ via `std::fs::write` into a directory from the
   house `fresh_tmpdir` helper (`tests/integration_partial_output.rs:26`), binary run with
   `--poly_a --poly_g -a AGATCGGAAGAGC --dont_gzip`.
4. **Every assertion two-sided.** "No output read ends in ≥10 A" is an under-trim detector that a
   build truncating every read to 1 bp would pass. Each case pairs it with an exact expected length,
   and the suite includes a no-`--poly_a` control that must be byte-identical to today's output.
5. `--paired --poly_a` R2 byte-identity as a committed test, not a manual before/after.

### Commit 3 — docs and CHANGELOG

1. `docs/.../reference/migration.md:28` — `| Poly-A trimming | None | Built-in --poly_a |` is false.
2. `docs/.../reference/migration.md:51-53` — the `-a 'A{N}'` note calls that form "a v0.6.x
   workaround that pre-dates `--poly_a`". It *was* `--polyA`. Rewrite with the three corrected Perl
   facts from Context: 10 A's, exclusive second pass with no `-q`/`--length`, and the `32:A:` /
   `_PolyA:32` read-ID rewriting that v2 does not do.
3. `docs/.../modes/poly.md:6` — "two built-in poly-X trimmers that the Perl version did not have" is
   wrong for poly-A. Add the stage order and the 2-colour interaction; add a 3'seq row to the "When
   to use which" table pointing at `--2colour`.
4. `src/cli.rs` `--poly_a` help text — extend "Runs after adapter trimming" to name poly-G.
5. `CHANGELOG.md` — one sentence, repo format, crediting @marenmsk.

## Efficiency

B1 is a statement reorder: zero cost, zero allocation, no complexity change.

## Integration

- **Reads:** `config.poly_a`, `config.poly_g`, `record.seq`. **Writes:** `record`,
  `TrimResult::poly_a_trimmed`, `TrimResult::poly_g_trimmed`.
- **Order:** quality → adapter → RRBS → **poly-G → poly-A** → N-trim → fixed clipping.
- `bp_after_cutadapt` is captured at `src/trimmer.rs:161`, ahead of both poly stages, so the #232
  invariant and the cutadapt-section report figures are untouched.
- `trimmer::trim_read` is the only implementation; `src/parallel.rs:495-496,1133` are the only
  external call sites and aggregate stats only. One edit site.
- **Report attribution shifts for `--poly_a` users on 2-colour data.** Affected reads move a few
  bases from the poly-G count to the poly-A count: measured on `polyAT_R1`, poly-A reads
  4,229 → 4,239 and poly-G bases 888 → 882. One line in the CHANGELOG entry.
- Report line order (`src/report.rs:401-425`) prints poly-A before poly-G. Left alone so existing
  reports stay comparable.

## Assumptions

1. The reporter's reads are `[insert][polyA][polyG]` — inferred from her description, **not
   confirmed on her data**. The reply asks for her `Reads with poly-A/T trimmed` count. If that is
   already high, the diagnosis is wrong and this plan does not apply.
2. Poly-G auto-detection fired on her run (she passed neither `--poly_g` nor `--no_poly_g`).
3. Dark-cycle poly-G is strictly terminal, so `[polyG][polyA]` is not a real 2-colour arrangement.
   This is what makes the mirror residual acceptable. It does **not** hold for template-derived
   G-runs (G-quadruplex, genomic poly-G ahead of a genuine tail).
4. Her v0.6.10 two-pass workflow and her v2.3.0 single command are not equivalent invocations — the
   Perl mode did no quality trimming or length filtering of its own.

## Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| 1 | The reported case is fixed end-to-end | `--quality 20 --phred33 --polyA --poly_g -a AGATCGGAAGAGC` on the synthetic set | `A40+G15` → 60, `A60+G20` → 40, `A20+G20` → 60 |
| 2 | Nothing else moves | `cargo test` (509 lib tests) + the full `validation` job | all green, no md5 step changes |
| 3 | **No over-trim, measured as a differential** | Sweep `illumina_10K`, `illumina10K_with_polyA`, `BS-seq_10K_R1`, `SRR24766921_RRBS_R1` before/after; count reads whose output shortens and the worst single delta | ≤ ~20 reads per fixture, no delta beyond the trailing G-run |
| 4 | R2 is byte-identical | `--paired --poly_a --poly_g` on `polyAT_R1/R2`, md5 both mates. **`--poly_g` must be forced** — auto-detection fires on 0 of 13,000 reads there, so without it the check is vacuous | R2 identical; R1 differs only where a G-run preceded a tail |
| 5 | The documented workaround still works | integration case with `--2colour 20` | same postcondition as 1 |
| 6 | The mirror residual is where we think it is | `[insert][polyG][polyA]` unit test | 75 bp, pinned |

Validation 3 replaces rev 1's "no read loses >tail length", which was not measurable when the tail
definition was itself in question and **would have passed the withdrawn B2**. The differential sweep
is the gate any future change to this stage must clear; it is the measurement rev 1 omitted.

## Questions or ambiguities

- **Q1 (Open).** `--no_poly_g --poly_a` on 2-colour data still leaves tails behind G-runs.
  Recommendation: leave it, honour the opt-out, document and pin with a test.
- **Q2 (Open).** Perl `--polyA` step in the `validation` job? Recommendation: **no**, and more
  firmly than rev 1 put it. Perl's mode is an exclusive second pass with no `-q`, so there is no
  single comparable v2 command — byte-identity is not merely the wrong gate, it is inexpressible.
- **Q3 (Open).** Take the fixed-point alternative instead of B1? See Follow-up. Recommendation: no,
  not in this fix — B1 is a reorder with a measured 12-read footprint and is backportable.
- **Q4 (Open).** #245 item B stays unresolved and is adjacent, since `-a 'A{N}'` is the fallback
  being offered to the reporter. Worth re-opening now that "the dedicated flag covers it" has failed.

## Follow-up (separate plans, not this fix)

- **Terminal-homopolymer skip.** Both reviewers independently proposed the same successor to the
  withdrawn B2: skip a maximal 3'-terminal run of a *single* non-A base of any length, then scan.
  It cannot delete arbitrary sequence because what it skips is provably a homopolymer; it needs no
  magic constant; and it fixes `A20 + G20` without depending on stage order. Reviewer A measured the
  collateral at **0.82% of reads at `min_run = 3`, 0.08% at `min_run = 5`**, against the withdrawn
  offset's 25%, and it leaves the 142 bp read untouched. Worth its own plan with the Validation 3
  sweep as the acceptance gate.
- **Fixed-point iteration.** Alternate poly-G / poly-A until neither cuts, bounded at two rounds.
  Order-independent, so it subsumes B1 and removes the mirror residual rather than trading one
  blocking case for another. Termination is trivial — each round only truncates.
- **5' prefix case.** `G15 + T40 + insert60` clips 0 on R2 today. Rev 1 claimed the 5' branch was
  inherently tolerant; it is not — both branches need `4·errors ≤ run`, and rev 1's probes passed
  only because random ACGT junk is ~25% T. File separately.

## Self-Review (rev 2)

**What changed and why.** B2 is gone. Rev 1 justified it on a single synthetic control
(`internal 16-A run`) and never swept the real fixtures; both reviewers did, and it removes 3'
sequence from ~23% of reads on ordinary libraries. Rev 1's own "Decision" argument against `A{20}`
applies to it with a worse worst case (142 bp vs 117 bp). The lesson is recorded in Validation 3 as
a standing gate rather than as prose.

**What I verified in the reviews rather than accepting.** The gate arithmetic — dropping a terminal
non-A base lowers `errors * 5` by 5 and `n - i` by 1, so the offset loop relaxes the global error
budget rather than skipping junk locally. That is arithmetic from the gate, and it explains both the
near-flat firing-offset histogram and the 142 bp cut. The two reviewers describe the same mechanism
in different words; they do not disagree.

**Known divergence between the reviews.** Their `A{20}` baseline figures differ from rev 1's
(3,886/2,731 and 2,986/2,731 vs rev 1's 3,009/1,965) because rev 1 measured after quality and adapter
trimming — i.e. at the pipeline position where the stage actually runs — and the reviewers measured
on raw reads. Different measurements, both internally consistent. It does not affect the decision,
which rests on the candidate sweep rev 1 never ran.

**Remaining risks.**

1. Assumption 1 is unconfirmed on the reporter's data. The reply asks for the one number that
   settles it; if it comes back high, re-diagnose before implementing.
2. The mirror residual is accepted on a biological argument (Assumption 3), not eliminated.
3. Her two workflows are not equivalent (Assumption 4), so some of the difference she sees may be
   quality/length filtering rather than poly-A trimming. The report figures will show this.
