# PLAN_REVIEW_B — `total_bp_after_trim` parity between the serial and parallel paths

**Reviewer:** B (independent; no contact with Reviewer A)
**Plan under review:** `plans/11102026_bp-after-trim-cores-parity/PLAN.md`
**Repo state:** `dev` @ `c0bbee4`, clean with respect to `src/`
**Date:** 2026-11-10

---

## Verdict

**The fix itself is correct and the plan's core technical claims hold up.** I verified the three
parallel sites, the value-not-plumbing claim, the identical accumulation placement, the separateness
of `total_bp_written`, the `bp_after_cutadapt` contract, the absence of any pinned expectation, and
the `parity_test_config()` blast radius — all as stated. Steps 1 and 2 can be implemented exactly as
written.

**The test half of the plan — which the plan itself calls "the part that carries the real risk" — has
two defects.** One claim it leans on twice is false: a paired-end serial-vs-parallel `TrimStats`
parity test *does* exist (`src/parallel.rs:2271`), asserting both R1 and R2 on a 12,000-pair fixture.
And the Step 3 helper as written cannot discriminate the `:543` site that Step 5 is supposed to cover,
because it sets no R2 clipping — reproducing, inside this plan, exactly the "test that can only ever
pass" defect the plan exists to fix.

I also found the plan's root-cause story to be wrong in a verifiable way, which matters less for the
fix than for the comment Step 2 writes and for the CHANGELOG step the plan does not have.

Nine findings below; two Critical, five Important, six Optional.

---

## 1. Claim verification

Everything the brief asked me to check, plus what I checked around it. "Verified" means I read the
source, not the plan.

### Verified true

| Plan claim | Evidence | Status |
|---|---|---|
| `TrimResult` in scope at all three sites | `res_r1` / `res_r2` bound at `src/parallel.rs:495-496`; `result` at `:1133`. Both already feed a dozen stat fields immediately below. | ✅ Value change, not plumbing |
| Accumulation placement already identical | Serial SE `src/trimmer.rs:342-345` (exit) → `:350` (accum) → `:353` (filter); parallel SE `:1151-1154` → `:1157` → `:1159`. Serial PE `:523-531` → `:536-537` → `:540`; parallel PE `:534-540` → `:542-543` → `:545`. | ✅ Identical; no reordering needed |
| `parity_test_config()` at `:1424` with 12 callers at the listed lines | `grep -n` returns exactly `1483, 1553, 1590, 1620, 1663, 1664, 1758, 1826, 1909, 1980, 2025, 2169`. | ✅ (see Optional O3 for a nuance that *strengthens* the conclusion) |
| Nothing asserts a `total_bp_after_trim` value | `grep -rn total_bp_after_trim src/ tests/` → only the field decl (`report.rs:24`), the merge (`:70`), the format (`:635-636`), the nine accumulation sites and one doc comment. No `"Total written"` string in `tests/` at all. | ✅ No expectation to update |
| `total_bp_written` is a separate field, untouchable by this change | `src/report.rs:21`; written only inside the `Pass` arms (`trimmer.rs:357`, `parallel.rs:1162`, `:551-552`). | ✅ |
| `bp_after_cutadapt` is the correct target semantics | Contract `src/trimmer.rs:64-71` cites #232 explicitly; captured `:165`, before RRBS (`:167`), poly-G/A, N-trim (`:245`) and clipping (`:249-279`). | ✅ Serial path is the correct one; the plan does **not** invert |
| Batch merging preserves the invariant | `TrimStats::merge` (`report.rs:70`) is a plain `+=`. The change is per-read inside `classify_*`, merged once per batch. Summation is associative, so no double-count or drop. | ✅ Right level |
| No other accumulation site | Only three `trim_read` call sites exist in `parallel.rs` (`:495`, `:496`, `:1133`). Both batch loops funnel through one `classify_*` call each (`:628`, `:1194`); `process_pairs` has two callers (`:386`, `:415`) but both reach `classify_paired`. Clumpify and passthrough route through the same functions. | ✅ Clumpify/passthrough covered by the same three lines |
| `--clumpify` does not open a fourth site | `--clumpify` forces the parallel path (`main.rs:1503`) but reuses `classify_single`/`classify_paired`. | ✅ |
| `TrimStats` supports whole-struct comparison | `#[derive(Debug, Default, PartialEq)]` at `report.rs:8`. | ✅ |

### False or overstated

| # | Plan claim | What the source says |
|---|---|---|
| **F1** | §2.3 "There is **no paired-end stats parity test at all**, so `:542` and `:543` have no coverage of any kind"; §11 "no paired stats parity test in the file at all, so Step 5 is new-test work" | **False.** `test_passthrough_serial_parallel_parity` (`src/parallel.rs:2271`) drives `trimmer::run_paired_end` against `run_paired_end_parallel` on a 12,000-pair fixture and asserts `assert_eq!(serial_stats.0, parallel_stats.0)` and `assert_eq!(serial_stats.1, parallel_stats.1)` — both R1 and R2 `TrimStats`. See Critical C1. |
| **F2** | §8 assumption 6, "`clip_r1: Some(5)` makes the two expressions differ for every read that survives to the accumulation point" | **False as stated.** `clip_5prime` / `clip_3prime` (`src/fastq.rs:110`, `:122`) are guarded on `self.seq.len() > n` — a read of length ≤ `n` is **skipped entirely**, not clamped. See Important I2. |
| **F3** | §2 "The `--clip_R1` delta is 10,000 reads × 5 bp less the filtered reads" | **Wrong reasoning, right number.** The accumulation sits *before* the filter on both paths, so filtered reads contribute to `total_bp_after_trim` and cannot reduce the delta. See Important I3. |
| **F4** | §5 Step 2, the parallel comment is "the serial comment with the #232 clause dropped — plausibly how the wrong expression got written" | **Not what happened.** The #232 fix (`f4b544d`, 2026-04-27) touched **0 lines** of `src/parallel.rs`. See Important I4. |
| **F5** | §2 "The serial path honours the contract at six sites, each carrying a comment that cites #232" | Only three do (`:350`, `:536`, `:537`). The uBAM-output sites `:672`, `:780`, `:781` carry **no comment at all**. Optional O2. |
| **F6** | §3 edge case, "`--cores 1` ... Bypasses the worker pool entirely (`src/cli.rs:358`)" | `cli.rs:358` is a doc comment on the `--cores` field, not control flow. The dispatch is `main.rs:1503`: `if cli.cores > 1 || cli.clumpify`. Conclusion survives only because of `cli.rs:747`. Optional O1. |

---

## 2. Logic review

### Critical

#### C1 — §2.3 / §5 Step 5 / §11: a paired serial-vs-parallel parity test already exists

`src/parallel.rs:2271`, `test_passthrough_serial_parallel_parity`, lines ~2334-2339:

```rust
assert_eq!(serial_stats.2, parallel_stats.2, ...);   // PairValidationStats
// TrimStats parity (R1 + R2) — already a documented invariant for
// the non-passthrough path; this just confirms passthrough doesn't
// perturb it.
assert_eq!(serial_stats.0, parallel_stats.0);
assert_eq!(serial_stats.1, parallel_stats.1);
```

It calls `crate::trimmer::run_paired_end` (serial, `cores` irrelevant) and then
`run_paired_end_parallel(..., 4, true, ...)` on the same 12,000-pair multiome fixture, and compares
R1 and R2 `TrimStats` separately. That is structurally the test Step 5 proposes to write.

Why it does not catch the bug: `passthrough_test_config()` (`:2169`) is built as
`let mut c = parity_test_config();` — so it inherits `clip_r1: None`, `clip_r2: None`,
`three_prime_clip_*: None`, `rrbs: false`, `poly_a: false`, `poly_g: false`, `trim_n: false`. Same
root cause as the single-end test: the *coverage* exists, the *config* is non-discriminating.

Concrete consequences for the plan:

1. §2.3's "no coverage of any kind" is wrong. `:542` and `:543` are exercised by a 12,000-pair
   serial-vs-parallel comparison today; it is blind to this one field, not absent.
2. §11's "Adjusted during review" bullet 2 — presented as one of the plan's three corrections and as
   the reason Step 5 is "a materially larger piece than the issue implies" — is a false correction.
   Step 5 is a config variant of an existing 80-line test, not new fixture work.
3. The implementer should be pointed at `:2271`. I would **not** simply add clipping to that test:
   it carries a passthrough reader/writer and a `~⅓ dropped` fixture, which adds moving parts to a
   stats assertion. A sibling modelled on it, with `None` for both passthrough params and
   `UnpairedLengths { r1: 35, r2: 35 }`, is cleaner and is maybe 30 lines once
   `gen_multiome_rows` / `write_multiome_fixture` are reused or replaced with the 30-read inline
   fixture the single-end test uses.

**Confidence: certain** (read the assertions directly). **Impact:** misdirected effort and a false
premise in the plan's self-review, not a wrong fix.

#### C2 — §5 Step 3 vs Step 5: the proposed helper cannot discriminate `:543`

Step 3's code block is:

```rust
fn parity_test_config_post_adapter() -> TrimConfig {
    TrimConfig {
        clip_r1: Some(5),
        three_prime_clip_r1: Some(3),
        ..parity_test_config()
    }
}
```

`is_paired` stays `false` and **no R2 clipping is set**. Step 5's prose then says to run "under
`parity_test_config_post_adapter()` with `is_paired: true` and `clip_r2`/`three_prime_clip_r2` also
set" — a requirement that lives only in prose, contradicted by the code block beside it.

If the implementer uses the helper as written for the paired test, then for R2:
`trimmer::trim_read(r2, config, true)` selects `config.clip_r2` / `config.three_prime_clip_r2`
(`src/trimmer.rs:250-259`), both `None`, nothing trims after the adapter stage, and
`res_r2.bp_after_cutadapt == r2.seq.len()`. Reverting `:543` alone then changes nothing, the paired
test passes, and `:543` is covered in name only.

That is the same shape as the defect under repair: an assertion that cannot fail under its own
fixture. Validation 5 is the backstop and would expose it — which is why this is a plan defect rather
than a guaranteed code defect — but the plan should not be shipping an artifact whose happy path
reproduces its own failure mode.

**Fix:** put `clip_r2: Some(5)`, `three_prime_clip_r2: Some(3)` and `is_paired: true` in the helper
and have the single-end test override `is_paired: false` (or, cleaner, give the helper a
`paired: bool` parameter, or declare the paired config inline in Step 5 and drop the pretence of
sharing). Either way the plan must contain the R2 values in code, not prose.

**Confidence: certain.** **Impact:** high — it is the exact failure class the plan targets.

### Important

#### I1 — All proposed assertions are differential; nothing pins the absolute semantics

Every assertion in Steps 4 and 5 has the form `assert_eq!(serial_stats, parallel_stats)`. That pins
*agreement between paths*, not *which value is right*. Three ways the suite stays green while the
behaviour is wrong:

- A future change flips the serial line to `record.seq.len()` — all nine sites agree, every test
  passes, #232 is silently reverted on **every** path including the Perl-parity one.
- Someone meeting a parity failure from the new tests "fixes" it in the wrong direction (the
  failure message is symmetric; nothing in it says which side is canonical).
- The refactor in Alternative A2 below unifies the paths onto `seq.len()`.

`test_bp_after_cutadapt_skips_rrbs_truncation` (`src/trimmer.rs:1076`) pins `TrimResult`'s field, but
**nothing pins that `total_bp_after_trim` is fed from it**. Given §1 calls the test work "the part
that carries the real risk" and §7 says the unit tests "are the only real guard", a guard satisfiable
by breaking the serial path is a hole in the load-bearing part of the plan.

The remedy is one line. Under `parity_test_config_post_adapter()` the existing 30-read fixture has
`length_cutoff: 0` (nothing filtered) and post-adapter lengths of 30 and 34 bp, all > 5 and > 3, so
every read is clipped by exactly 8 bp:

```rust
assert_eq!(
    stats_serial.total_bp_after_trim,
    stats_serial.total_bp_written + 8 * 30,
    "total_bp_after_trim must be the pre-clip length (#232), not the written length"
);
```

A weaker but refactor-proof version is `assert!(total_bp_after_trim > total_bp_written)`. Add one of
them to Step 4 and the paired equivalent to Step 5.

**Confidence: high.** I did not run the fixture; the arithmetic follows from `length_cutoff: 0`, the
two fixture read lengths, and clip order (5' first, then 3' on the shortened read).

#### I2 — Assumption 6 is false: clipping is skipped, not clamped, for short reads

`src/fastq.rs:110-120`:

```rust
pub fn clip_5prime(&mut self, n: usize) -> Option<String> {
    if n > 0 && self.seq.len() > n { ... } else { None }
}
```

Strict `>`. A read of length exactly 5 under `clip_r1: Some(5)` is left untouched — not truncated to
zero. `clip_3prime` (`:122`) is identical.

§11 risk 2 flags this as an open uncertainty and reaches the right conclusion ("still non-zero on the
fixtures used, so the test still discriminates"): on the 30-read fixture every read is 30 or 34 bp
after adapter trim, far above both clip lengths, so the Step 4 test discriminates unconditionally.
The plan is safe. But the assumption should be recorded as **resolved false** rather than left open,
because it is the premise of §9 validation 1's arithmetic, and because a future fixture with short
reads would quietly weaken the test.

**Confidence: certain** (read the guard). **Impact:** low on the fix, real on §9's expectations.

#### I3 — §2's explanation of the 49,980 bp delta is wrong, and the plan's own data shows why

§2: "The `--clip_R1` delta is 10,000 reads × 5 bp less the filtered reads." Filtered reads cannot
reduce this delta — the accumulation sits *before* `filters::filter_single_end` on both paths
(`trimmer.rs:350` then `:353`; `parallel.rs:1157` then `:1159`), which §2 itself establishes two
paragraphs earlier.

The arithmetic that does work: 49,980 = 5 × 9,996, so 9,996 of 10,000 reads were clipped and **4 were
not**, because they were ≤ 5 bp after quality and adapter trimming (finding I2). Those same 4 reads
are then necessarily below the default `length_cutoff`, which is why 9,996 records were written — the
two 9,996s are the same four reads seen twice, not a coincidence.

This is inherited verbatim from the issue text, so the issue is wrong the same way. It matters
because §11 risk 2 asks for exactly this to be "confirmed on first run rather than assuming", and the
confirmation is already sitting in the plan's own measurement table.

**Confidence: high.** Derived from the plan's numbers plus the verified guard; I did not re-run the
fixture.

#### I4 — The root cause is not a dropped comment; `f4b544d` never touched `parallel.rs`

Step 2 speculates that the wrong expression came from writing the serial comment without its #232
clause. The history says otherwise:

- `82d1e34` (2026-04-17, *"Fix total_bp_after_trim missing in parallel path and DRY cleanup"*) added
  all three parallel sites. Its diff introduces exactly `stats_r1.total_bp_after_trim += r1.seq.len();`,
  `stats_r2... += r2.seq.len();`, `stats... += record.seq.len();`. At that point `bp_after_cutadapt`
  **did not exist**, and the serial path also used `seq.len()` — so the parallel sites were correct
  and consistent when written.
- `f4b544d` (2026-04-27, the #232 fix, PR #235) introduced `bp_after_cutadapt` and switched the
  serial sites from `seq.len()` to `result.bp_after_cutadapt`. `git show f4b544d --stat` lists
  **no `src/parallel.rs`**. The regression is an orphaned call site, not a mis-transcribed comment.

Two consequences:

1. Step 2's rationale is fiction. The output it prescribes (a one-line, current-state comment) is
   still the right output under the project's comment rules, so **the step stands as written** — only
   its justification needs deleting. Do not let the corrected history leak into the comment itself.
2. `f4b544d`'s commit message and the CHANGELOG entry it added both claim the fix applies "in both
   single-end and paired-end pipelines" — phrasing that reads as full coverage while silently
   meaning *serial* SE and *serial* PE. The user-facing record of #232 has therefore been inaccurate
   for `--cores ≥ 2` since April. Worth one clause in the new CHANGELOG entry (see I5).

**Confidence: certain** (read both diffs).

#### I5 — No CHANGELOG step

`CHANGELOG.md` has an empty `### Unreleased` section, and every comparable user-visible fix in the
2.4.0 block has an entry. This change moves a number in the MultiQC-facing text report for every
`--cores ≥ 2` user; it is exactly the kind of thing someone deciding whether to upgrade wants to see.
The plan has no CHANGELOG step and §7's integration list does not mention the file.

Follow the file's existing bold-headline-plus-short-paragraph format and keep it to the one sentence
the project's conventions call for.

**Confidence: certain.** **Impact:** an omitted deliverable, not a correctness risk.

#### I6 — Step 5 is the only step with no code and no signature check

Steps 1-3 give exact code; Step 5, which the plan itself calls the materially larger piece, gives
prose. The two entry points it names are not trivial:

- `trimmer::run_paired_end` (`src/trimmer.rs:385-396`) takes 10 arguments — two
  `&mut dyn RecordSource`, an `Option<&mut dyn RecordSource>` passthrough, three `&mut FastqWriter`,
  two `Option<&mut FastqWriter>` unpaired, `&TrimConfig`, `UnpairedLengths` — and returns
  `(TrimStats, TrimStats, PairValidationStats)`.
- `run_paired_end_parallel` (`src/parallel.rs:93-107`) takes 13, including `cores`, `gzip`,
  `UnpairedLengths` and `clump_layout: Option<ClumpLayout>`, returning the same triple.
- It also needs a second input file (the existing single-end fixture is one file) and the `open_fq`
  helper at `:1346`.

Finding C1 largely dissolves this: `:2271` already wires both calls correctly and can be copied. But
as the plan stands, the step with the most unknowns has the least specification, and §10's "Taken:
implementer's choice" covers Step 4, not Step 5.

**Confidence: certain.**

---

## 3. Assumptions

| # | Plan's assumption | Assessment |
|---|---|---|
| 1 | `bp_after_cutadapt` is the correct target semantics | **Verified.** Documented contract `trimmer.rs:64-71`, captured `:165` ahead of RRBS/poly/N-trim/clip, #232-backed. The plan does not invert. |
| 2 | `trim_read` mutates in place; `seq.len()` after the call is final | **Verified.** `trim_read(record: &mut FastqRecord, ...)`, `trimmer.rs:91`. |
| 3 | Accumulation point already correct on both paths | **Verified**, both SE and PE, both relative to `--discard_untrimmed` and to the filter. |
| 4 | `total_bp_written` is intentionally the final length | **Verified.** `report.rs:21`, written only in `Pass` arms. |
| 5 | `parity_test_config()` must stay byte-for-byte | **Verified and understated** — see O3. |
| 6 | `clip_r1: Some(5)` differs for every surviving read | **False** — finding I2. Harmless on the fixture in use. |
| 7 | `--cores 2` suffices to exercise the pool | Reasonable; the existing tests use 4 and there is no per-core branch on this field. The round-robin dispatch (`parallel.rs:106-113`) means 2 cores still exercises multi-worker merge. |

### Unstated assumptions the plan relies on

- **That the whole-struct `assert_eq!` will not fail on an unrelated field once clipping is on.**
  I checked the plausible candidates: `rrbs_r2_clipped_5prime` is gated on `config.rrbs`
  (`parallel.rs:514`, `trimmer.rs` equivalent), which stays false; `length_cutoff: 0` means no
  filter-count divergence; `clip_5prime_applied` is only consumed under `rrbs`. I expect clean
  parity post-fix, but this is an inference, not a test run.
- **That `--cores 1` is the only serial entry.** True today via `main.rs:1503`, but only because
  `cli.rs:747` rejects `--clumpify` below `--cores 2`. The plan cites a doc comment for this (O1).
- **That the report line is the only consumer.** Verified — `report.rs:635-636` is the sole read of
  the field, and the JSON path does not carry it.
- **That `843` is the baseline, not `842`.** `grep -c '#\[test\]'` over `src/*.rs` + `tests/*.rs`
  gives 843 attributes (some may be `cfg`-gated). §9's "842 tests pass" should be measured on the
  base commit before the change rather than asserted (O5).

---

## 4. Efficiency

§6 is correct and needs nothing added. `bp_after_cutadapt` is a `usize` already computed and stored
on `TrimResult` for every read on both paths (`trimmer.rs:165`, moved into the struct at `:285`);
reading a field costs less than the `String::len()` call it replaces. No allocation, no extra pass,
no complexity or memory change, nothing touching the hot alignment loop.

The new tests add cost: Step 4 re-runs a 30-read fixture (negligible). If Step 5 is modelled on
`:2271` it inherits a 12,000-pair fixture at `cores 4` with gzip — measurable but in line with what
that file already pays, and the plan's 30-read alternative would be cheaper. Either is fine; prefer
the small fixture, since stats parity does not need batch-boundary coverage for this field (the merge
is a plain sum).

---

## 5. Validation sufficiency

### What §9 gets right

- **Validation 5 is sound and necessary.** `:542` and `:543` write to different `TrimStats`
  instances — `stats_r1: &mut TrimStats, stats_r2: &mut TrimStats` at `parallel.rs:484-485` — so a
  half-applied fix is a real possibility and a single revert would not expose it. Reverting each line
  separately, with the paired test asserting R1 and R2 separately, is the right discipline. The plan
  is right to call this the step most likely to be skimped.
- **Validation 2 and 3** (record identity, `total_bp_written` unchanged) correctly bound the blast
  radius, and both are checkable without trusting the new tests.
- **§7's observation that the CI validation matrix "proves nothing about this change"** is correct
  and unusually honest: the matrix runs `--cores 1`.

### Where a wrong result still slips through

1. **The serial direction** — finding I1. The highest-value gap, and the cheapest to close.
2. **The `:543` site under the helper as specified** — finding C2. Validation 5 catches it, but only
   if run; the plan's own text predicts this is the step that gets skimped.
3. **The rendered report text is never asserted by the suite.** Every proposed check is on `TrimStats`
   or on a manual re-run. `report.rs:635` is a straight format of the field so the exposure is small,
   but there is no regression guard on the line that the issue is actually about. Optional O4.
4. **Validation 1's expected numbers depend on I2/I3.** The figures (644,525 / 490,527) are right;
   the arithmetic offered to justify them is not. Keep the figures as the expectation and drop the
   `reads × 5` reasoning, or correct it to "5 bp per read longer than 5 bp after adapter trim".
5. **Nothing checks the other direction on `--rrbs`.** Validation 1 covers it manually at the binary
   level, but no unit test uses `rrbs: true` for stats parity. The RRBS case is the one #232 was
   filed about. Not required — clipping is the stronger discriminator — but a second config
   (`rrbs: true, non_directional: false`) would cost one more `parity_test_config` spread and cover
   the historically-regressing path directly. Optional.

### Does the Step 4 test fail for the right reason?

Yes, as far as I can determine by reading. Pre-fix, with clip 5+3 on 30- and 34-bp post-adapter
reads and `length_cutoff: 0`, serial accumulates 30/34 per read and parallel accumulates 22/26 —
`total_bp_after_trim` differs by 240 bp over 30 reads and every other field is unaffected (see the
unstated assumption above). The `{:#?}` dump in the existing assertion names the field on failure.
Post-fix both sides read `bp_after_cutadapt`. I did not run it.

---

## 6. Alternatives

### A1 — Use the existing paired test as the Step 5 vehicle (recommended)

Covered under C1. Not a different approach so much as the cheaper route to the same place: copy
`:2271`'s call wiring, drop passthrough, swap the config. Reduces Step 5 from "the materially larger
piece" to the smallest step in the plan.

### A2 — Make the invariant structural instead of test-enforced

The plan fixes three expressions and adds tests, leaving the mechanism that produced the bug intact.
That mechanism is duplication: `trimmer.rs:320-352` and `parallel.rs:1125-1163` are near-identical
per-read stats-accumulation blocks, as are `trimmer.rs:483-537` and `parallel.rs:495-543` — **nine
accumulation sites for one semantic, maintained by hand in two files.** `82d1e34`'s own subject line
was "and DRY cleanup", and the field has now drifted between the paths twice (once as "missing
entirely", once as "wrong expression").

Two shapes, cheapest first:

- **Funnel the writes.** Give `TrimStats` a single `fn record_trimmed(&mut self, result: &TrimResult)`
  that owns `total_bp_after_trim += result.bp_after_cutadapt` (and whichever neighbouring fields are
  purely `TrimResult`-derived), and make it the only way to touch the field. Small diff, no control
  flow moves, and it kills the whole class — a future field added to `TrimResult` gets one home
  instead of nine call sites. This is a genuinely attractive follow-up.
- **Have the serial path call `classify_single` / `classify_paired`.** Those functions are already
  documented as "sink-agnostic" (`parallel.rs:1123-1124`) and return a write/drop verdict, so the shape
  fits. But the serial loops also carry passthrough lockstep, `--retain_unpaired` writer handling and
  `--rename`, and they write records inside the match arms — so this is a real refactor on the path
  the Perl byte-identity matrix guards. Not for a bug fix.

I agree with excluding both from *this* plan. But §2's "Why it survived" lists two reasons and
neither is the duplication, and §11's "Remaining risks" does not mention it. One line in §7 and a
follow-up issue for the funnel version would be proportionate.

### A3 — `--rrbs` as the discriminator instead of clipping

Weaker: the RRBS delta is data-dependent (only reads with the MspI signature), which is why the
issue's RRBS delta is 6,560 bp against 49,980 for clipping. The plan's choice of `clip_r1` is right,
and adding `three_prime_clip_r1` to widen it to both ends is a good touch. No change recommended —
though see validation gap 5 for adding RRBS as a *second* config rather than a replacement.

### A4 — Scope boundary on #529 and #542

I agree with both exclusions. `#529`'s uBAM accumulation is not on `dev` (the serial `*_to_bam`
sites at `trimmer.rs:672`, `:780-781` are already correct), and `#542`'s `MultithreadedWriter` is
orthogonal. The plan's note that Step 5 gives #529 a template is right and becomes more right under
A1. One addition: when #529 lands, the funnel from A2 would make its path correct by construction
rather than by remembering.

---

## 7. Action items

### Critical

1. **Correct §2.3 and §11: a paired serial-vs-parallel `TrimStats` parity test exists at
   `src/parallel.rs:2271`** (`test_passthrough_serial_parallel_parity`, asserting R1 and R2 on 12,000
   pairs; blind to this field only because `passthrough_test_config()` at `:2169` inherits
   `parity_test_config()`). Rewrite Step 5 to build on it, and strike §11's "no paired stats parity
   test in the file at all" correction. *(Finding C1 — confidence: certain)*

2. **Put `clip_r2` / `three_prime_clip_r2` / `is_paired` into the Step 3 helper in code, not into
   Step 5's prose.** As specified, the helper leaves R2 untrimmed after the adapter stage, so a
   paired test using it cannot fail when `:543` is reverted — the plan's own failure mode, inside the
   plan. *(Finding C2 — confidence: certain)*

### Important

3. **Add one absolute-value assertion** so the tests pin the semantics, not just path agreement —
   e.g. `assert_eq!(stats_serial.total_bp_after_trim, stats_serial.total_bp_written + 8 * 30)` in
   Step 4 and its paired equivalent in Step 5. Without it, a regression on the *serial* side keeps
   the whole suite green and silently reverts #232 everywhere. *(I1 — high)*

4. **Record assumption 6 as resolved false.** `clip_5prime`/`clip_3prime` (`src/fastq.rs:110`, `:122`)
   skip reads of length ≤ `n` rather than clamping them. Harmless on the 30-read fixture; it is the
   premise of §9 validation 1's arithmetic. *(I2 — certain)*

5. **Fix §2's delta explanation.** Filtered reads contribute to `total_bp_after_trim` (the
   accumulation precedes the filter), so they cannot reduce the delta. 49,980 = 5 × 9,996 means 4
   reads were ≤ 5 bp after adapter trim and so went unclipped. Worth correcting on #541 too, since
   the plan inherited it from there. *(I3 — high)*

6. **Replace Step 2's rationale with the actual history** (and keep the prescribed one-line comment
   unchanged): `82d1e34` wrote all three parallel sites correctly for their time, and `f4b544d` — the
   #232 fix — switched the serial sites while touching **0 lines** of `parallel.rs`. Consider one
   clause in the CHANGELOG noting that #232's fix had not applied at `--cores ≥ 2`. *(I4 — certain)*

7. **Add a CHANGELOG step.** `### Unreleased` is empty and every comparable user-visible fix in the
   2.4.0 block has an entry; this one moves a MultiQC-facing number for all `--cores ≥ 2` users.
   *(I5 — certain)*

### Optional

8. **§3 edge case citation.** `cli.rs:358` is a doc comment; the dispatch is `main.rs:1503`
   (`cli.cores > 1 || cli.clumpify`), safe only because `cli.rs:747` requires `--cores >= 2` with
   `--clumpify`. Cite the dispatch and the gate. *(O1)*

9. **§2's "six sites, each carrying a comment that cites #232"** — only `:350`, `:536`, `:537` do;
   `:672`, `:780`, `:781` carry no comment. Minor, and it further undercuts Step 2's causal story
   (three uncommented sites got the expression right). *(O2)*

10. **§5 Step 3's "12 tests"** is 12 call *sites* across 11 tests, one of which (`:2169`) is a second
    helper with its own callers — so the blast radius of mutating `parity_test_config()` is larger
    than stated, and the plan's conclusion is stronger than its evidence. Reword to "12 call sites"
    and note the transitive helper. *(O3)*

11. **Consider one suite-level assertion on the rendered report line.** Nothing in `src/` or `tests/`
    asserts the text `Total written (filtered)` — the artifact the issue is about. `report.rs:635` is
    a straight format, so exposure is small. *(O4)*

12. **Measure the test baseline rather than asserting it.** §9 validation 6 says 842; I count 843
    `#[test]` attributes on `c0bbee4`. Record the number from a run on the base commit. *(O5)*

13. **File a follow-up for funnelling the writes** — a single
    `TrimStats::record_trimmed(&TrimResult)` owning `total_bp_after_trim` would end a class of bug
    that has now bitten twice across nine hand-maintained accumulation sites, and would make #529
    correct by construction. Out of scope here; worth a line in §7 and an issue. *(A2)*

---

## 8. Separation of verified from suspected

**Verified by reading the source or git history:** every row of §1; the `:2271` assertions; the
`clip_5prime`/`clip_3prime` guards; the `f4b544d` and `82d1e34` diffs; the `main.rs:1503` dispatch and
`cli.rs:747` gate; the `TrimStats` derives and `merge`; the absence of any `total_bp_after_trim` or
`"Total written"` assertion in `src/` and `tests/`; the three `trim_read` call sites in
`parallel.rs`; the entry-point signatures; the empty `### Unreleased` block.

**Inferred, not executed:** that the Step 4 test fails pre-fix and passes post-fix on exactly
`total_bp_after_trim` with no collateral field divergence; the `+ 8 * 30` figure in action item 3;
the identification of the 4 unclipped reads in I3 (derived from the plan's own measurements plus the
verified guard). I built nothing and ran no tests, per the brief.

---

`/Users/fkrueger/Github/TrimGalore/plans/11102026_bp-after-trim-cores-parity/PLAN_REVIEW_B.md`
