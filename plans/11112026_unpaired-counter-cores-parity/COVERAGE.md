# Plan Coverage Report — #544 `r1_unpaired` / `r2_unpaired` cores parity

**Mode:** B — code vs. plan (post-implementation)
**Plan(s):** PLAN.md (revision 2)
**Date:** 2026-11-11
**Verdict:** COMPLETE

**Code audited:** working tree, branch `dev` @ `86aa301`, change uncommitted
(`git diff` hash `438c58d8…`). Repo left on `dev` with the change in place; no source file
was edited by this audit.

## Summary

- Total items: 21 (10 implementation steps + 11 validations)
- DONE: 20
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 1 (Step 10 — deviation declared in §12 and confirmed accurate)

Diff confined to three files, matching §12's claim byte for byte:

| File | +/− |
|---|---|
| `src/parallel.rs` | +203 / −32 |
| `CHANGELOG.md` | +7 / −0 |
| `docs/src/content/docs/reference/changelog.md` | +7 / −0 |

`src/trimmer.rs` is unmodified (`git diff --stat src/trimmer.rs` is empty), so §12's claim that it
was mutated during the checks and reverted exactly is verified.

## Coverage ledger

| # | Item | Source | Status | Evidence |
|---|------|--------|--------|----------|
| 1 | Two unconditional increments removed from `classify_paired`'s `TooShort` arm | §5 Step 1 | DONE | `src/parallel.rs:567-572` — arm now sets only `pairs_removed`, `stats_r1.too_short`, `stats_r2.too_short`, then returns `PairOutcome::Unpaired { r1_ok, r2_ok }` |
| 2 | Increments added inside both `if let Some(w) && r*_ok` guards in `process_pairs`, next to the writes | §5 Step 2 | DONE | `src/parallel.rs:643` (after `r1.write_to`), `:649` (after `r2.write_to`) |
| 3 | Arm's comment replaced with a one-line invariant | §5 Step 3 | DONE | `src/parallel.rs:637-638` — the six-line historical comment is gone |
| 4 | **Sibling** test, four fixture shapes, unequal R1/R2 totals | §5 Step 4 | DONE | `src/parallel.rs:1729` `test_parallel_serial_unpaired_counter_parity`; shapes `:1743-1749`; `length_cutoff: 30` `:1770`; `UnpairedLengths { r1: 20, r2: 20 }` `:1773` |
| 5 | Assertions under both rescue modes, in a non-passthrough host | §5 Step 5 | DONE | `src/parallel.rs:1775` `for retain in [false, true]`; passthrough reader and writer both `None` at `:1843` and `:1846` in the `run_paired_end_parallel` call, and `None` at the matching serial args |
| 6 | Hand-derived absolute anchors, not read back from a run | §5 Step 6 | DONE | `src/parallel.rs:1750-1751` consts `5 + 2` / `3 + 2`; asserted `:1871-1884` |
| 7 | Fixture discrimination guard | §5 Step 7 | DONE | `src/parallel.rs:1899-1905` — both non-zero **and** unequal |
| 8 | `gen_multiome_rows` edited directly; doc comment's false claim removed | §5 Step 8 | DONE | `src/parallel.rs:2452-2458` (new doc comment), `:2473-2475` (R2 now 48 bp unconditionally) |
| 9 | CHANGELOG entry naming both effects, mirrored into the docs copy | §5 Step 9 | DONE | `CHANGELOG.md:8-13`; `docs/src/content/docs/reference/changelog.md:14-19` — byte-identical |
| 10 | Mutation checks | §5 Step 10 | DEVIATED (declared) | §12 deviation verified accurate — see below |
| V1 | No-rescue run core-invariant | §9 | DONE | Re-derived: text block **absent** at `--cores 1` and `4`; JSON `0`/`0` both |
| V2 | Same with clipping | §9 | DONE | Re-derived: same, with `--clip_R1 5 --clip_R2 5` |
| V3 | Rescue run unchanged | §9 | DONE | Re-derived: `Unpaired read 2 kept: 4` at both core counts; 2 files written |
| V4 | Counts match the files | §9 | DONE | Re-derived: R2 file holds exactly 4 records, R1 file 0, at both core counts |
| V5 | JSON values move | §9 | DONE (post-fix half) | Re-derived: keys present, `0`/`0` no-rescue, `0`/`4` rescue, both core counts. Pre-fix `0`/`4` at `--cores 4` rests on §2's table, not on my own run |
| V6 | Records unchanged | §9 | DONE (post-fix half) | Re-derived: R1 `53d619cd94e561f2`, 9,996 records — **matches §12's recorded hash exactly**; identical across cores 1/4 and both rescue modes |
| V7 | R1 site covered | §9 | DONE | Independently reproduced — see mutation table |
| V8 | R2 site covered | §9 | DONE | Independently reproduced |
| V9 | "Always zero" is not parity | §9 | DONE | Independently reproduced; the **anchor** fires, not parity |
| V10 | #245 not reintroduced | §9 | DONE | Shape C (2 pairs, 26/26 bp) present and asserted; its contribution is visible in the measured totals (7 = 5 + 2, 5 = 3 + 2) |
| V11 | Suite and gates | §9 | DONE | `cargo test` exit 0, **845 passed / 0 failed**; `cargo fmt --all -- --check` exit 0; `cargo clippy --all-targets --release -- -D warnings` exit 0, 0 warnings |

## Independent verification performed

I re-derived every validation rather than reading §12's figures back.

### Mutation checks (V7–V9, Step 10)

Run in an isolated `git archive HEAD` + working-diff copy under the scratchpad, so no repo source
file was touched. Baseline in that copy passes all six parity tests, which rules out a dead harness.

| Mutation | Result | Matches §12? |
|---|---|---|
| Parallel R1 increment removed (`:643`) | **FAILED** — parity assertion, retain mode; serial `r1_unpaired: 7` vs parallel `0` | yes |
| Parallel R2 increment removed (`:649`) | **FAILED** — parity assertion; serial `r2_unpaired: 5` vs parallel `0` | yes, including the "serial 5 vs parallel 0" figure |
| Both parallel increments removed (Step 10 **as written**) | **FAILED** via the **parity** assertion (7/5 vs 0/0), *not* the retain assertion | yes — this is the deviation §12 declares |
| All four increments removed (serial *and* parallel) | **FAILED** via the **absolute anchor**: `serial r1_unpaired (retain) must equal the records written`, left `0` / right `7`. Parity held at 0/0 | yes |
| Original bug reintroduced in `classify_paired` | **TWO** tests failed: `test_parallel_serial_unpaired_counter_parity` **and** `test_passthrough_serial_parallel_parity` | yes — Step 8's un-blinding is real |

The fixture's discriminating power is confirmed empirically, not just by reading the shape table: the
failure output shows `pairs_analyzed: 20` (3+5+2+4+6), `pairs_removed: 14` (3+5+2+4),
`r1_unpaired: 7` (B+C), `r2_unpaired: 5` (A+C). All four shapes land in their intended branches and
the two totals are unequal, so Step 4's requirement and Step 7's guard both bite.

### Call-graph and parity surface

- Exactly four increment sites exist in the crate — `src/trimmer.rs:588`, `:594` (serial) and
  `src/parallel.rs:643`, `:649` (parallel) — plus the field-wise `merge` at `src/report.rs:148-149`.
  No fourth increment site, as §2 claims.
- The parallel predicate now matches the serial one exactly: `Some(writer) && r*_ok`, increment
  immediately after the write.
- Perl 0.6.11 grounding independently checked (`git show 0.6.11:trim_galore`): `++$read_1_printed`
  sits inside `if ($retain)` immediately after the four `print UNPAIRED1` statements, and the report
  line reads `Number of unpaired read 1 reads printed`. Assumption 1 is upstream-sanctioned.
- `retain_unpaired = unpaired_r1_path.is_some()` confirmed at `src/parallel.rs:108`.

### Scope boundaries held

- No code from #542, #529 or #541 was implemented. The only `#541` strings in the diff are unchanged
  CHANGELOG context lines beneath the new entry.
- `test_parallel_serial_bp_after_trim_parity_paired` (`:1614`) is untouched and still builds its
  config from `parity_test_config_post_adapter()`, which inherits `length_cutoff: 0` — so the plan's
  stated reason for a sibling rather than a mutation holds, and #541's anchors are undisturbed.
- Exactly one `#[test]` added, none removed, consistent with §12's "844 + 1 new = 845".

## §12's claims vs. the code

**§12 matches the code.** Every row of its outcome table, its five mutation results, its six
verification rows, its numstat, its "`src/trimmer.rs` untouched" claim and its recorded output hash
`53d619cd…` all reproduce. The declared Step 10 deviation is accurately described: reverting only the
two *parallel* increments leaves serial counting, so it fails through the parity assertion exactly as
the single-side reverts do and proves nothing extra; the informative mutation removes all four, which
makes `serial == parallel` hold at 0/0 so only Step 6's anchor can fire. I confirmed both halves of
that reasoning by running both mutations.

Four presentational notes, none of them coverage gaps:

1. **Step 3 comment length.** The comment is two physical lines wrapping one sentence
   (`src/parallel.rs:637-638`); §12 calls it "one line". The substance matches the plan and the
   repo's comment conventions — one short sentence, current state, no historical narrative.
2. **Step 10 count.** §12's table lists four mutations while the row says "Five performed"; the fifth
   (both parallel increments only) appears in the deviation paragraph below the table rather than in
   it. The count is correct; only the presentation is split. I ran all five.
3. **§1 vs. Step 4.** §1 says "Both [blind tests] have to change or the fix ships unverified", while
   Step 4 forbids mutating `test_parallel_serial_bp_after_trim_parity_paired`. The implementation
   followed §5: one blind test was un-blinded (the passthrough one) and the second's coverage is
   carried by the new sibling instead. The underlying goal — the fix does not ship unverified — is
   met, demonstrated by the five mutations above. The tension is internal to the plan document.
4. **§10's third open question** ("grep `docs/` for 'unpaired' before finishing") is not recorded as
   done in §12. I ran it: the only matches are the two lines of the changelog entry just added, so no
   docs page needed editing and the plan's taken assumption was correct.

### Work beyond the plan — all declared in §12

- A fifth fixture shape ("P", 6 pairs at 40/40 bp) that passes the length filter, giving the test
  non-filtered pairs alongside the four `TooShort` shapes. Declared as "plus 6 passing pairs".
- Step 7's guard strengthened from "both non-zero" to "both non-zero **and** unequal". Declared.
- Four extra assertions comparing each counter against the record count of its unpaired file, for
  both the serial and the parallel run (`src/parallel.rs:1892-1895`). Within Step 5's
  "equals the records actually written".
- Two `needless_borrows_for_generic_args` fixes in the new test. Declared, and structurally
  consistent: `FastqWriter::create` is generic over `AsRef<Path>` (so temporaries are passed
  unborrowed) while `run_paired_end_parallel` takes plain `&Path` (so `&dir.join(…)` is correct
  there).

## Gaps

None. No item is PARTIAL or MISSING.

Two validations are verifiable only in their post-fix half from the current tree, and I say so rather
than claiming more: V5's "was 0/4 at `--cores 4`" and V6's "identical to pre-fix" compare against a
pre-fix binary I do not have. V6's recorded hash reproduces exactly, the fix's mechanism (moving two
`usize` increments inside guards already evaluated for the writes) cannot alter records, and the
reintroduce-the-bug mutation reproduces the pre-fix behaviour in the test harness — so the residual
risk is nil, but the comparison itself rests on §2's and §12's measurements.

## Verdict

**COMPLETE.** All ten implementation steps and all eleven validations are satisfied. The single
deviation (Step 10) is declared in §12 and its description is accurate.
