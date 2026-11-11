# PLAN_REVIEW_A — `total_bp_after_trim` parity between the serial and parallel paths

**Reviewer:** A (independent; no coordination with Reviewer B)
**Plan under review:** `plans/11102026_bp-after-trim-cores-parity/PLAN.md`
**Repo state verified against:** `dev` @ `c0bbee4` (confirmed via `git rev-parse`/`git log`)
**Date:** 2026-11-10
**Method:** source reading only; no build, no test run (noted where a claim is therefore unverified).

---

## Verdict

**The fix is correct and the approach is sound.** The target semantics, the three sites, the
value-not-plumbing claim, the accumulation placement, the site inventory and the "nothing pins the
current value" claim all check out in the source. Steps 1 and 2 can be implemented as written.

What needs work is the *diagnosis* and the *guard*, not the fix:

- Two factual claims about the surrounding code are **wrong**, and both produce wrong follow-up
  actions: there *is* a paired serial-vs-parallel stats parity test (§2.3, §11), and #529 does *not*
  add its own accumulation (§1, §7, §11.3 — it calls the very functions this plan fixes).
- The new tests reproduce the structural weakness they exist to replace: nothing in them asserts
  that the chosen config still discriminates, so a later edit returns them to "cannot fail".
- No CHANGELOG step, and no end-to-end check of the paired figure — the two paired sites get unit
  coverage only.

Nine findings below: 0 Critical, 5 Important, 4 Optional.

---

## 1. Claims verified against the source

Everything in this table I read in the tree; none of it is taken from the plan.

| Plan claim | Verdict | Evidence |
|---|---|---|
| `bp_after_cutadapt` captured post-quality/post-adapter, pre-RRBS/poly/N-trim/clip | **True** | `src/trimmer.rs:165` sits after the adapter round loop (`:120-156`) and before RRBS (`:167`), poly-G (`:205`), poly-A (`:220`), `trim_ns` (`:246`) and fixed clipping (`:249-277`) |
| Contract documented at `trimmer.rs:64-71`, captured `:161-165` | **True** | exact lines |
| `TrimResult` in scope at `parallel.rs:542`, `:543`, `:1157` | **True** | `res_r1`/`res_r2` bound in `classify_paired` and used at `:498-531`; `result` bound at `:1133` and used at `:1134-1149`. Value change, not plumbing |
| Accumulation placement already identical between paths | **True** | serial SE `trimmer.rs:341-350` vs parallel `parallel.rs:1151-1157`; serial PE `trimmer.rs:521-537` vs parallel `parallel.rs:534-543`. Both sit after the `--discard_untrimmed` `continue`/early return and before `filters::filter_*` |
| `total_bp_written` is a separate field and cannot be perturbed | **True** | `src/report.rs:21` vs `:24`; written only in the `Pass` arms (`parallel.rs:552-553`, `:1163`) |
| JSON report does not carry the field | **True** | `report.rs:962-971` emits `total_bp_written` under `basepair_processing`; `total_bp_after_trim` appears only in the text writer at `:632-637` |
| Nothing in `src/` or `tests/` asserts a `total_bp_after_trim` value | **True** | `grep -rn total_bp_after_trim src/ tests/` → definition, merge, two format args, the 9 accumulation sites, one doc comment. No CI step diffs report text either: the `validation` job md5s only `gzip -dc`'d FASTQ, plus one adapter-line `grep` (`ci.yml:604`) and `test -f` existence checks (`:623-626`) |
| `parity_test_config()` sets no post-adapter stage | **True** | `parallel.rs:1424-1452`: `trim_n:false`, all four clips `None`, `rrbs:false`, `poly_a:false`, `poly_g:false` |
| `parity_test_config()` has 12 callers | **True** | 14 textual matches = definition `:1424` + doc mention `:2165` + the 12 call sites the plan lists. Step 3's reasoning holds, and understates it: `:2169` is inside `passthrough_test_config()`, which has 5 further callers, so mutating the shared config would reach ~16 tests |
| Site inventory is complete (6 serial + 3 parallel) | **True** | `grep -rn total_bp_after_trim src/` finds no other writer. `demux.rs`, `specialty.rs`, `clump_only.rs` never touch the field |
| Per-batch merge cannot double-count or drop | **True** (plan assumes it implicitly) | per-batch `TrimStats` built at `parallel.rs:348-349`/`:1094`, folded at `:267-268`/`:1037` through `TrimStats::merge` (`report.rs:66-71`), which is plain field-wise addition. The change is per-read inside the per-batch accumulator, so the invariant is order-independent and preserved |
| Clumpify and passthrough need no separate fix | **True** (plan does not state it) | both funnel through `classify_single` (`parallel.rs:1194`, sole caller) and `classify_paired` (`:628`, sole caller) |
| §11 risk 2: clipping may not apply to short reads | **True — now confirmed, not a suspicion** | `FastqRecord::clip_5prime`/`clip_3prime` (`src/fastq.rs:110-137`) are no-ops when `seq.len() <= n` (strict `>`). See F6 |
| §2.3 "no paired-end stats parity test at all" | **FALSE** | see F1 |
| §1/§7 "the uBAM output type does not exist on `dev` yet"; "#529 adds its own accumulation" | **FALSE** | see F2 |
| §2 "six serial sites, each carrying a comment that cites #232" | **False for 3 of 6** | see F7 |

Not verified (no build): validation 6's "842 tests pass", and the two measured bp figures in §2 —
though the §2 arithmetic is internally consistent once F6 is applied, which is indirect support.

---

## 2. Logic review

### F1 — Important. There *is* a paired serial-vs-parallel stats parity test; §2.3 and §11 are wrong about it

**Concerns:** §2 "Why it survived" item 3; §5 Step 5; §11 "Adjusted during review" bullet 2.
**Confidence:** high — read directly.

`test_passthrough_serial_parallel_parity` (`src/parallel.rs:2271`) drives
`trimmer::run_paired_end` (`:2288`) against `parallel::run_paired_end_parallel` (`:2311`) on the same
fixture and asserts, at `:2336-2338`:

```rust
assert_eq!(serial_stats.0, parallel_stats.0);
assert_eq!(serial_stats.1, parallel_stats.1);
```

Those are whole-struct `TrimStats` parity assertions for R1 and R2 across the two paired paths —
precisely the shape Step 5 proposes to create. So `:542`/`:543` are not uncovered; they are covered
by an assertion that **cannot fail**, for exactly the reason the single-end one cannot:
`passthrough_test_config()` (`:2168-2173`) is `parity_test_config()` plus `length_cutoff = 20` and
`is_paired = true` — no clipping, no RRBS, no poly, no N-trim — so `bp_after_cutadapt ==
record.seq.len()` by construction.

Three concrete consequences:

1. **§2's "Why it survived" is missing its most telling item.** The lesson is not "one test was
   built on a clean config and a second case had no test"; it is that **two** independent parity
   tests, written years apart, both inherited the same clean config and were both vacuous on this
   field. That is a fixture-design problem, and it is the argument for F3's in-test guard.
2. **Step 5 is cheaper than §11 bills it.** "New-test work, not an extension — a materially larger
   piece than the issue implies" is overstated: there is a working paired harness at `:2271-2325`
   (reader/writer wiring, `UnpairedLengths`, the 10-argument `run_paired_end` call, the 13-argument
   `run_paired_end_parallel` call) to copy verbatim. The estimate should come down.
3. **The implementer needs to be told this test stays green.** Under its own config both expressions
   are equal, so the fix does not move it — but someone fixing `:542`/`:543` and then seeing a
   second paired parity test in the same file will stop to work out why it was not already failing.
   One sentence in §7 saves that.

**Action:** correct §2.3 and §11; add `src/parallel.rs:2271` to the "why it survived" list; state in
§7 that it is unaffected; re-scope Step 5.

### F2 — Important. The #529 scope rationale is wrong in both directions, and produces a wrong follow-up

**Concerns:** §1 "Out of scope"; §7 "#529"; §11 remaining risk 3.
**Confidence:** high — read `dev` and the PR diff.

Two separate errors:

**(a) The uBAM output path exists on `dev` and already does the right thing.** §1 says "the uBAM
output path's own accumulation added by #529 (that type does not exist on `dev` yet)".
`trimmer::run_single_end_to_bam` (`src/trimmer.rs:638`) and `trimmer::run_paired_end_to_bam`
(`:708`) are on `dev` and accumulate `result.bp_after_cutadapt` correctly at `:672`, `:780`, `:781`.
Those three *are* half of the plan's own "six serial sites" — §2's table labels them "Serial
single-end" / "Serial paired" without noting they are the uBAM variants, so the plan simultaneously
cites them as evidence and declares them nonexistent.

**(b) #529 adds no accumulation; it inherits this one.** PR #529 is *"feat(ubam): multithreaded uBAM
input and output with `--cores`"*, +430 lines in `src/parallel.rs`, open against `dev`. Its diff
contains **no** `total_bp_after_trim` line. Its new batch functions delegate:

- `process_single_batch_bam` → `classify_single(record, config, &mut stats)` (diff line 1402)
- `process_paired_batch_bam` → `classify_paired(...)` (diff line 1430)

So the parallel uBAM path routes through `parallel.rs:1157` and `:542`/`:543` — the three lines this
plan fixes. Rebasing #529 onto this fix corrects its uBAM path for free.

**Why it matters concretely:** §7 and §11.3 direct a follow-up ("once it lands, the same fix will be
needed there", said to be noted on #541). Acting on that means hunting for a site that does not
exist, or adding a second accumulation that double-counts. The note on #541 should be corrected to
the opposite statement: *#529 inherits the fix through `classify_*` and needs no separate change.*

The scope **decision** (don't touch #529) is still right. The right residual risk to record is
merge ordering, not recurrence: #529 touches `src/parallel.rs` heavily, so whichever lands second
rebases — the three one-line changes are far from #529's insertions, but Steps 3–5's test additions
and #529's test additions both land in the same `#[cfg(test)]` module.

### F6 — Optional (but fix the prose). §2's explanation of the measured delta is wrong, and the data already settles §11 risk 2

**Concerns:** §2 "Measured, on `c0bbee4`" last sentence; §9 validation 1; §11 remaining risk 2.
**Confidence:** high for the mechanism; the measured numbers themselves I did not reproduce.

§2 says "The `--clip_R1` delta is 10,000 reads × 5 bp less the filtered reads." Filtering cannot
explain any shortfall: accumulation happens **before** `filters::filter_single_end`
(`parallel.rs:1157` then `:1159`), so filtered reads still contribute in full. The plan states this
correctly two paragraphs earlier and then contradicts it here.

The actual mechanism is the one §11 flags as unconfirmed: `FastqRecord::clip_5prime`
(`src/fastq.rs:110-113`) is a **no-op when `seq.len() <= n`**, strict `>`, so a read already ≤ 5 bp
after adapter trimming loses nothing to `--clip_R1 5`. And the measurement confirms it: 49,980 / 5 =
9,996 — exactly the written-read count, so precisely 4 reads were ≤ 5 bp, and those same 4 were the
ones the length filter then dropped. So §11 risk 2's "worth confirming on first run rather than
assuming" is already answered by the plan's own data: **clipping is skipped, not clamped, for short
reads**, and the delta is `clipped reads × 5`, not `all reads × 5 less filtered`.

Harmless for Steps 3–5 as written — `parity_test_config()`'s fixture is 30/34 bp post-adapter,
comfortably above 5+3 — but it matters for whoever picks a fixture for Step 5: a short-read fixture
weakens the discriminator read-by-read.

### F7 — Optional. Three of the six serial sites have no comment at all

**Concerns:** §2 ("each carrying a comment that cites #232"); §5 Step 2.
**Confidence:** high.

`trimmer.rs:350` and `:536-537` carry the three-line #232 comment. `:672`, `:780` and `:781` (the
`to_bam` functions) carry **none** — the accumulation sits on a bare line after the
`--discard_untrimmed` block. This slightly undercuts Step 2's story ("the serial comment with the
#232 clause dropped — plausibly how the wrong expression got written"), and it means Step 2 leaves
the codebase with three commented serial sites, three bare serial sites and three one-line parallel
comments. Not worth fixing in this plan; worth not asserting the opposite.

### F4, F5 — covered in the action items

Two findings are process/scope gaps rather than logic flaws and are stated in full, with evidence, as
action items 5 and 6: the missing CHANGELOG step (F4) and the unmentioned `--clumpify` impact plus
the `--cores 1` edge-case citation (F5).

### F9 — Optional. The Step 3 helper's name promises more than it delivers

**Concerns:** §5 Step 3 vs Step 5.
**Confidence:** high.

`parity_test_config_post_adapter()` as drafted sets only `clip_r1` and `three_prime_clip_r1`, yet
Step 5 needs `is_paired`, `clip_r2` and `three_prime_clip_r2` as well — stated in prose, absent from
the code block. Either widen the helper (set all four clips and leave `is_paired` to the caller) or
have Step 5 show the three extra assignments explicitly. Low risk, but it is the kind of gap that
produces a paired test where only R1 discriminates — which is exactly what validation 5 is designed
to catch, so the two steps are at cross purposes as written.

---

## 3. Assumptions

| # | Plan's assumption | Assessment |
|---|---|---|
| 1 | `bp_after_cutadapt` is the correct target, not merely the serial behaviour | **Sound, and well grounded.** The field exists for nothing else (`trimmer.rs:64-71`), `test_bp_after_cutadapt_skips_rrbs_truncation` (`:1076`) locks its capture point, and #232 is the stated reason. The plan does not invert |
| 2 | `trim_read` mutates in place, so `record.seq.len()` afterwards is the final length | **Verified.** `trim_read(&mut record, …)` (`trimmer.rs:1133` in the parallel caller); every stage mutates through `truncate`/`clip_*` |
| 3 | Accumulation point already correct on both paths | **Verified** (table §1) |
| 4 | `total_bp_written` stays as it is | **Verified** |
| 5 | `parity_test_config()` must stay byte-for-byte | **Verified, and stronger than stated** (≈16 tests via `passthrough_test_config`) |
| 6 | `clip_r1: Some(5)` makes the expressions differ for every surviving read | **False as stated, harmless as used.** Not "every" — only reads longer than 5 bp (F6). On the fixtures in play it differs for essentially all of them |
| 7 | `--cores 2` suffices to exercise the pool | **Sound**; the existing tests use 4 |

**Unstated assumptions I had to check myself:**

- *Batch-level merging preserves the invariant* — holds (§1 table). Worth one line in §8 since the
  whole fix lives inside a per-batch accumulator.
- *No other consumer reads the field* — holds; `report.rs:632-637` is the only one.
- *The per-read stats blocks are otherwise identical between paths under a clipping config, so the
  whole-struct assertion will pass after the fix* — holds. `rrbs_r2_clipped_5prime` is the only
  other clip-sensitive field and both paths gate it on `config.rrbs` (`parallel.rs:535`,
  `trimmer.rs:505`-region), which the proposed config leaves `false`. This is load-bearing for
  Steps 4–5 passing at all and deserves stating.
- *`--clumpify` reaches the parallel path* — see F5.

---

## 4. Efficiency

§6 is correct and complete. `bp_after_cutadapt` is a `usize` already materialised in `TrimResult`
for every read on both paths; a field read replaces an O(1) `String::len()`. No allocation, no extra
pass, no change to complexity, memory or the worker/batch structure. `TrimResult` is constructed
unconditionally either way, so nothing becomes newly reachable. Nothing further to add.

---

## 5. Validation sufficiency

The validation table is unusually good for a plan this small — in particular, splitting validation 5
into two separate reverts is exactly right, and I can confirm it is **both necessary and sound**:
`:542` writes `stats_r1` and `:543` writes `stats_r2`, which are distinct `TrimStats` instances
(`parallel.rs:348-349`), merged into distinct totals (`:267-268`), and surfaced as distinct tuple
elements. A test asserting both R1 and R2 therefore fails on either revert independently, and a
single-revert check would pass a half-applied fix.

Two gaps where a wrong result could still ship unnoticed:

### F3 — Important. The new tests have no guard against silently becoming vacuous again

**Concerns:** §5 Steps 3–5; §9 validations 4–5.
**Confidence:** high on the gap; medium on my exact remedy being the best form.

Steps 4–5 produce tests that discriminate **today**, verified once by hand in Step 6. Nothing in the
tests themselves records *why* they discriminate. Change `parity_test_config_post_adapter()` later —
drop a clip, lower a fixture's read length below the clip width (F6), or alter clip semantics — and
both tests return to the "cannot fail" state, with no failure anywhere. That is the defect this plan
exists to fix, reproduced one layer up. F1 shows it is not hypothetical: it has now happened twice in
this file.

The file already has the idiom. `test_parallel_serial_trim_stats_parity` ends with:

> `// Sanity-check that the test fixture actually exercised adapter trimming — otherwise we'd be
> passing trivially.` (`src/parallel.rs:1517-1524`)

Add the matching guard for *this* invariant to each new/extended test, asserting on the **serial**
(reference) stats:

```rust
assert!(
    stats_serial.total_bp_after_trim > stats_serial.total_bp_written,
    "fixture must exercise a post-adapter trim stage, or this parity assertion cannot fail"
);
```

I checked this is satisfiable and non-trivial on the described fixture: 30 reads at 30 bp / 34 bp
post-adapter, `clip_r1: Some(5)` + `three_prime_clip_r1: Some(3)` = 8 bp each, `length_cutoff: 0` so
nothing is filtered ⇒ `total_bp_after_trim` 960, `total_bp_written` 720. Under the *clean* config the
two are equal, so the assertion also encodes the distinction between the two configs rather than just
asserting a number. It must be on the serial stats: on the parallel ones it would be a second way of
asserting the fix rather than a statement about the fixture.

This converts Step 6's one-off manual mutation check into something that keeps holding.

### F8 — Important. No end-to-end check of the paired figure; validations 1–3 are all single-end

**Concerns:** §9 validations 1–3.
**Confidence:** high.

Validation 1 re-runs `--clip_R1 5` on `BS-seq_10K_R1.fastq.gz` and `--rrbs` on
`SRR24766921_RRBS_R1.fastq.gz` — both single-end, both exercising `:1157` only. `:542` and `:543`
get unit coverage and nothing else, so the plan never confirms that *the artifact a user reads* is
core-count-invariant for paired runs. Given the fix is to a report line and the report is the
MultiQC-facing surface, that is the wrong place to stop.

Add one row: `--paired --clip_R1 5 --clip_R2 5` on `test_files/BS-seq_10K_R{1,2}.fastq.gz` at
`--cores 1` and `--cores 4`, comparing `Total written (filtered)` in **both** `_trimming_report.txt`
files. Two runs, a couple of minutes, and it closes the unit-test-to-artifact gap for the two sites
the plan itself identifies as the riskiest.

### Smaller notes on §9

- Validation 2 ("identical to the pre-fix output") requires keeping the pre-fix outputs; say so, or
  it gets skipped and the record-identity claim rests on the post-fix run alone.
- Validation 6's "842 tests pass" I did not verify (no build). Treat the number as indicative.
- Validation 7 is right to call clippy/fmt necessary-not-sufficient given the rustc version gap.

---

## 6. Alternatives

**(A) Make the wrong expression unspellable — recommended as a follow-up, not here.** Add
`TrimStats::record_post_trim(&mut self, result: &TrimResult)` (or `accumulate_post_trim`) beside
`merge` in `src/report.rs` and route all nine sites through it. ~15 lines. It removes the *class* of
defect rather than this instance, and #529 would have been immune by construction. Trade-off: it
touches the six serial sites this plan deliberately leaves alone. Behaviour-identical, so no
Perl-parity exposure — but it does widen the diff on the byte-identity-sensitive file, and this plan
is right to keep the bug fix minimal. Worth filing as its own issue; the F2 correction is the
argument for it (the duplication is why this diverged, and `bases_quality_trimmed` is copy-pasted
across the same five blocks, `trimmer.rs:324/486-487/650/737-738`, `parallel.rs:498-499/1134`).

**(B) Have the serial entry points delegate to `classify_single`/`classify_paired`.** Those are
already documented "sink-agnostic" (`parallel.rs:1123-1124`) and are the single implementation the
parallel and uBAM-parallel paths share. Folding the serial paths into them would make serial/parallel
stats parity **structural** rather than test-enforced, retiring this entire defect class and the
parity tests' reason to exist. This is a real refactor across the byte-identity-critical surface and
must not ride along with a bug fix; noting it because the plan's §7 correctly observes that "the unit
tests are the only real guard", and (B) is the alternative to needing a guard at all.

**(C) Parameterise the parity test over a config matrix** (clip / rrbs / poly-A / poly-G / trim-n)
rather than one hand-picked config. Lower value than it looks: after the fix the expression is
`result.bp_after_cutadapt` unconditionally, so any single post-adapter stage discriminates all three
sites equally. A matrix buys breadth against *other* per-field drift, not against this one. Dongze's
standing proptest-harness offer is the natural home for it; out of scope here.

**(D) Mutating `parity_test_config()` instead of adding a sibling** — the plan considered and
rejected this, correctly, and the caller count makes the case stronger than stated (see the
`parity_test_config()` row in §1). No reason to revisit.

---

## 7. Action items

### Critical

None. The fix as specified in Steps 1–2 is correct and safe to implement as written.

### Important

1. **F2 — correct the #529 rationale, and fix the note on #541.** §1's "that type does not exist on
   `dev` yet" is false (`trimmer.rs:638`, `:708`, accumulating correctly at `:672`, `:780`, `:781`),
   and §7/§11.3's "#529 adds its own accumulation / the same fix will be needed there" is false —
   #529's `process_{single,paired}_batch_bam` call `classify_single`/`classify_paired`, so it
   inherits this fix on rebase. Replace the residual risk with the merge-ordering one. *(Prevents a
   wasted follow-up or a double-counting "fix" on the uBAM path.)*
2. **F1 — correct "there is no paired-end stats parity test".** `test_passthrough_serial_parallel_parity`
   (`src/parallel.rs:2271`, assertions at `:2336-2338`) is one, and is vacuous for the same reason as
   the single-end test. Add it to §2's "why it survived", state in §7 that it stays green, and
   re-scope Step 5 down — its harness is copyable.
3. **F3 — add an in-test guard that the fixture still discriminates**, on the serial stats, in both
   new/extended tests (`assert!(stats_serial.total_bp_after_trim > stats_serial.total_bp_written)`),
   following the existing idiom at `src/parallel.rs:1517-1524`. Without it, Step 6's manual mutation
   check is a one-time measurement of a property the tests do not hold.
4. **F8 — add a paired end-to-end validation row**: `--paired --clip_R1 5 --clip_R2 5` on
   `BS-seq_10K_R{1,2}.fastq.gz` at `--cores 1` vs `4`, comparing `Total written (filtered)` in both
   report files. `:542`/`:543` currently get no artifact-level check at all.
5. **F4 — add a CHANGELOG step.** `CHANGELOG.md` has an empty `### Unreleased` section and the repo
   gives every user-visible change a linked entry (the original #232 fix has one —
   `docs/src/content/docs/reference/changelog.md:1388`). This changes a number users read in the
   MultiQC-facing report on every `--cores ≥ 2` run. The plan has no step for it. *(Format: match
   the file's existing entries; note the tension with the global "one sentence" rule and let Felix
   settle it.)*

### Optional

6. **F5 — name `--clumpify` in §7's impact statement and fix §3's `--cores 1` citation.** Dispatch is
   `cli.cores > 1 || cli.clumpify` (`src/main.rs:1503`, `:1736`, `:2030`), and `--clumpify` requires
   `--cores >= 2` (`cli.rs:745-747`), so **every Clumpify run to date has reported the post-clip
   figure** — the largest affected population, and unmentioned. §3's edge-case row cites
   `src/cli.rs:358`, which is a doc comment, not the dispatch; cite `main.rs:1503` and note the
   `|| cli.clumpify` disjunct (moot in practice, but the plan should not rest a correctness claim on
   prose).
7. **F6 — fix §2's delta arithmetic and close §11 risk 2.** Filtering cannot explain the shortfall
   (accumulation precedes it); `clip_5prime` being a no-op for `len <= n` (`src/fastq.rs:110-113`)
   can, and 49,980/5 = 9,996 = the written-read count confirms it. Restate as "clipped reads × 5"
   and mark risk 2 resolved rather than "worth confirming".
8. **F9 — make the Step 3 helper and Step 5 consistent.** The helper sets only the R1 clips while
   Step 5 needs all four plus `is_paired`; show the extra assignments or widen the helper, so the
   paired test cannot end up discriminating on R1 only.
9. **F7 — drop the "each carrying a comment that cites #232" claim.** `trimmer.rs:672`, `:780`,
   `:781` carry no comment. Cosmetic, but it is part of Step 2's rationale.

### Also worth adding to §8 (one line each)

- The per-batch `TrimStats` merge (`report.rs:66-71`) is field-wise addition, so a per-read change
  inside a per-batch accumulator is order-independent and cannot double-count — the structural
  reason the fix is safe at this level.
- `rrbs_r2_clipped_5prime` is the only other clip-sensitive field and both paths gate it on
  `config.rrbs`, which the new config leaves `false` — the reason the whole-struct assertion can pass
  at all under a clipping config.

---

## 8. Summary

The plan correctly identifies the defect, the three sites, the target semantics and the test gap, and
its validation design — especially the two separate reverts — is better than most plans of this size.
The implementation steps need no change. What needs change is three factual claims about the
surroundings (F1, F2, F7), one piece of arithmetic (F6), and two additions that decide whether this
fix stays fixed: an in-test discrimination guard (F3) and a paired end-to-end check (F8), plus the
missing CHANGELOG step (F4).
