# Plan Coverage Report

**Mode:** B (code vs. plan — the design plan *is* the implementation plan; there is no `IMPL.md`)
**Plan:** `plans/08172026_rename-ubam-inert-clip-refusal/PLAN.md` (v2, incl. §12)
**Code:** uncommitted working tree in `/Users/fkrueger/Github/TrimGalore/wt-450`, branch
`fix/450-inert-clip-refusal`, `HEAD = 4f7b907` (= `origin/dev` base, confirmed)
**Date:** 2026-08-17
**Verdict:** **COMPLETE** — 0 MISSING, 0 PARTIAL, 1 DEVIATED (documented as D1), 1 item deferred to CI
by construction (V9)

## Summary

- Total ledger items: **52**
  - Implementation outline: 18
  - B1 truth-table rows (each promised a test): 12
  - Validation V1–V9 (incl. V6b): 10
  - B3 edge-case rows: 6
  - Context six-site fold-in: 6
- DONE: **50**
- DEVIATED (documented): **1** (outline step 10 → §12 D1)
- DEFERRED (unverifiable locally): **1** (V9, the `Validate vs Perl TrimGalore` CI job)
- PARTIAL: **0**
- MISSING: **0**

### Independently established facts (not taken from §12)

| Check | Command | Result |
|---|---|---|
| Test count | `CARGO_TARGET_DIR=$TMPDIR/pm-target cargo test --no-fail-fast` | **788 passed, 0 failed, 0 ignored across 18 targets; exit 0** — matches §12's 788 |
| Net new test functions | `git grep -c -F -- "#[test]"` worktree vs `4f7b907` | 789 vs 782 → **delta +7**, matches §12's "7 new" |
| Formatting | `cargo fmt --all -- --check` | exit 0, no output |
| Lints | `CARGO_TARGET_DIR=… cargo clippy --all-targets --release -- -D warnings` | **exit 0, 0 warnings / 0 errors** |
| Diff size | `git diff --stat` | 8 files, **+438 / −88** — matches §12 |
| Tree is pristine (no mutation left applied) | read `applies`, `read_2`, `any_effective` bodies; `git ls-files --others` | `applies` = `paired \|\| !self.read_2()`; `read_2` = `matches!(self, ClipR2 \| ThreePrimeClipR2)`; `any_effective` = `applies && value.is_some()`; **no untracked files** (no leftover snapshot) |
| Nothing pre-existing weakened | `git diff -U0 \| grep '^-' \| grep -E 'fn \|assert'` | Only 3 removals, all plan-authorised: `any_set` (step 3), `every_override_knows_which_read_it_belongs_to` (step 11a, explicit "replace"), and that test's own `assert_eq!(r.overrides.len(), 4)` |

## Coverage ledger — implementation outline (18 steps)

| # | Item | Status | Evidence / notes |
|---|---|---|---|
| 1 | `ClipFlag` with `ALL`, `name()`, `read_2()`, `applies(paired)`; widen module doc | DONE | `library.rs:30` (enum), `:40` (`ALL` = R1, 3'R1, R2, 3'R2 — as the plan's snippet specified), `:49`, `:59` (`matches!`, no string inspection — Assumption 4 held), `:64`. Module doc gained the `ClipFlag` paragraph |
| 2 | `ResolvedClips::value` and `UserClips::value` | DONE | `library.rs:189`, `:210` — plain variant matches |
| 3 | Replace `any_set()` with `any_effective(paired)`; delete `any_set` | DONE | `library.rs:201`. `git grep any_set` over `src tests docs` → **no hits** |
| 4 | `ClipOverride.flag: ClipFlag`, drop `read_2`; `pick` reshaped | DONE | `library.rs:136`. `read_2:` and `o.read_2` → **no hits** anywhere. Four `pick` sites pass variants (§12 D4 records that dropping the 4th param makes #447's call-site-bool invariant unrepresentable rather than merely satisfied) |
| 5 | `flag_summary` uses `applies()`/`value()`, keeps its own order | DONE | `library.rs:168-186`, order `ClipR1, ClipR2, ThreePrimeClipR1, ThreePrimeClipR2` — unchanged from base. Pins green: `flag_summary_lists_all_four_for_paired`, `flag_summary_omits_read_2_for_single_end`, `integration_library.rs` |
| 6 | Extract `Cli::user_clips()`, `pub(crate)`, with the plan's doc | DONE | `cli.rs:1234`; `effective_clips()` now delegates, so still exactly one `UserClips` literal |
| 7 | Rewrite the warning ladder over its own order | DONE | `cli.rs:1156-1176`. Order `ClipR1, ThreePrimeClipR1, ClipR2, ThreePrimeClipR2` — identical to base. **Form differs, logic does not:** the plan wrote `value.is_some() && (all_four \|\| f.read_2())`; the code uses `if !all_four && !flag.read_2() { continue; }` then `if let Some(n)`, which is the same predicate. Assumption 3 held — it reads `user_clips()`, not `effective_clips()` |
| 8 | `main.rs:1128` + `report.rs:287` filters → `o.flag.applies(...)`; `.name()` at both print sites | DONE | `main.rs:1128`, `report.rs:287`; `main.rs:1131`/`report.rs:292` → `o.flag.name()` |
| 9 | `report.rs:898` → `json_escape_string(o.flag.name())` | DONE | Escape call kept as specified. JSON strings pinned by `integration_library.rs:313-323` (asserts `overrides.len() == 4` **and** each `flag` string) plus `report::tests::test_json_report_records_all_four_overrides` — both green, so B2's byte-identity claim is CI-held |
| 10 | `main.rs:398` → `any_effective(cli.paired)`; collapse the 7-line comment to one | **DEVIATED** | Predicate landed (`main.rs:393`). Comment is **two** lines (`main.rs:389-390`), not one — documented as **§12 D1** with its reason (criterion + siting will not fit on one line; two is CLAUDE.md's maximum). Deviation is documented, so it counts as DEVIATED, not PARTIAL |
| 11 | Drift guards (a) name-vs-`read_2` oracle replacing the old test, (b) `ALL` completeness, (c) exact names | DONE | (a) `library.rs:362` `every_clip_flag_knows_which_read_it_belongs_to` — old test removed as instructed; (b) `:377` `all_holds_every_clip_flag_once`, exhaustive match + name dedup, and the code comment states the honest limit the plan asked for; (c) `:395` `clip_flag_names_are_the_command_line_spellings`, all four with leading `--` |
| 12 | Rows 1–2 tests, separate output dirs | DONE | `integration_ubam_out.rs:1233`. Asserts exit 0, warning text present, **absence of** `--rename is refused`, `sp_trimmed.bam` written, no `:clip5:`/`:clip3:` in the read name. `fresh_tmpdir("tg_450_inert_{tag}")` per case — collision hazard handled |
| 13 | Rows 4, 5, 6, 7 tests | DONE | `integration_ubam_out.rs:1290`, four `Case`s: `three_prime_clip_r1_se`, `clip_r2_paired`, `tp_clip_r2_paired`, `r1_live_beside_inert_r2_se`. Each asserts non-zero exit, the refusal text, **and** that no BAM was written |
| 14 | Row 8 as a five-preset loop; `:1044` stays | DONE | `integration_ubam_out.rs:1360` loops `emseq, accel, zymo, scbs, pbat`. `rename_into_ubam_refused_for_fastq_input_via_library_preset` still at `:1044`, green |
| 15 | Pin the ladder order in `integration_inert_clip_flags.rs`; file otherwise unmodified | DONE | Added `find("--clip_R2 5") < find("--three_prime_clip_R2 7")` inside `each_inert_flag_gets_its_own_line`, plus new `all_four_inert_flags_warn_in_clipping_order` (`:250`) asserting 4 warning lines and 4 ascending positions. **No existing assertion in that file was changed or removed** (verified against the diff) |
| 16 | `outputs.md:100` and `:145` | DONE | Both hunks present; each adds the "applies to this run" qualifier, `:145` also spelling out the SE Read 2 case |
| 17 | `cli.rs:312-313` `--help` overstatement | DONE | Now "…when any input is FASTQ **and a clipping flag applies to the run**". No `--help` snapshot test exists, as the plan noted, so this is unpinned by design |
| 18 | CHANGELOG bullet citing #450 | DONE | Own bullet at the top of Unreleased → Changes, bold headline + **one** sentence, cites #450. Matches the file's existing bullet shape and the global one-sentence rule |

## Coverage ledger — B1 truth table (12 rows)

Every row confirmed **by test name in the run log**, not inferred from the total. Assertion direction
checked against the row's "After" column in each case.

| Row | Invocation | Expected after | Test (confirmed `... ok`) | Status |
|---|---|---|---|---|
| 1 | `--clip_R2 5` SE | runs + warning | `rename_into_ubam_accepted_for_inert_read_2_clip_flag` (case `clip_r2`) | DONE |
| 2 | `--three_prime_clip_R2 7` SE | runs + warning | same test, case `tp_clip_r2` | DONE |
| 3 | `--clip_R1 3` SE | refused | `rename_into_ubam_refused_for_fastq_input` (`:1006`) | DONE |
| 4 | `--three_prime_clip_R1 3` SE | refused | `rename_into_ubam_still_refused_where_a_clip_can_apply` (case `three_prime_clip_r1_se`) | DONE |
| 5 | `--clip_R2 5 --paired` | refused | same test, case `clip_r2_paired` | DONE |
| 6 | `--three_prime_clip_R2 7 --paired` | refused | same test, case `tp_clip_r2_paired` | DONE |
| 7 | `--clip_R1 3 --clip_R2 5` SE | **refused** | same test, case `r1_live_beside_inert_r2_se` | DONE |
| 8 | `--library <all five>` SE | refused | `rename_into_ubam_refused_for_every_library_preset` + `…_via_library_preset` (`:1044`) | DONE |
| 9 | `--hardtrim5` | refused | `rename_into_ubam_refused_for_hardtrim_fastq_input` (`:1172`) | DONE |
| 10 | `--hardtrim3` | refused | `rename_into_ubam_refused_for_hardtrim3_fastq_input` (`:1204`) | DONE |
| 11 | `--rename`, no clip flag | runs | `rename_without_clip_flag_is_accepted_into_ubam` (`:1078`) | DONE |
| 12 | uBAM input + clip flag | runs | `rename_into_ubam_accepted_for_ubam_input` (`:1109`, `--clip_R1 5`) | DONE |

**Row 4's pre-existing hole is closed and measured:** `git grep -c -F -- "three_prime_clip_R1" --
tests/integration_ubam_out.rs` returns **0 at `4f7b907`** and **1 now**. The plan's claim that the second
entry of the R1 pair was untested was true, and it no longer is.

## Coverage ledger — Validation V1–V9

| # | Status | Basis |
|---|---|---|
| V1 (defect reproduces at base) | DONE per §12 | Executed by the implementer (warning + refusal + exit 1, no BAM). **Not re-run by me**: it needs a build of `4f7b907`, and the audit constraint forbids modifying this worktree or touching the main one |
| V2 (narrowed case runs, identical to omitting the flag) | DONE per §12; half independently green | §12 records exit 0, warning kept, refusal absent, BAM byte-identical to a `--rename`-only run, plus a positive control (`--clip_R1 5` reports DIFFER, so the IDENTICAL verdict is not vacuous). The committed half — exit 0, warning, BAM written, no `:clip5:`/`:clip3:` — passed in my run |
| V3 (all twelve rows are tests) | **DONE, independently verified** | See the table above: 12/12 by name |
| V4 (narrowing is load-bearing) | DONE per §12; direction verified analytically | Mutation tightens (`any_effective` counts R2 always). Under it, `rename_into_ubam_accepted_for_inert_read_2_clip_flag` asserts `ok` → red, and `any_effective_discounts_read_2_on_single_end` asserts `!any_effective(false)` → red. §12's "2 failures" is consistent with exactly those two |
| V5 (paired direction is load-bearing) | DONE per §12; direction verified analytically | Mutation loosens (`applies` drops `paired \|\|`). Rows 5/6 assert `!ok` → red; `flag_summary_lists_all_four_for_paired`, `paired_still_announces_r2_overrides`, `test_json_report_records_all_four_overrides` also depend on it. §12's "6 failures across five sites" is consistent, and V4/V5 failing **different** tests is the empirical proof they are no longer the same mutation (v1's defect) |
| V6 (read ownership load-bearing) | DONE per §12 | Flipping `read_2()` for `ClipR2` breaks the `_R2` oracle plus every paired display pin — §12's 11 failures |
| V6b (R1 arm load-bearing) | DONE per §12 | Flipping `ThreePrimeClipR1` breaks the oracle **and** row 4 — the row v1 omitted |
| V7 (sites 1, 2, 4, 5, 6 are pure refactors) | DONE per §12 | 27 cells, 102 artifacts, 1 differing cell (the #450 shape), 0 vacuous. Required two normalisations (§12 **D2**: `--version` provenance line, and `command_line`'s argv[0]) and a `CONTROLS` set (§12 **D3**). **Not reproducible from the tree** — the harness was not committed, by design. The plan's own Remaining-risk 1 anticipated this and converted the top risk to CI via step 15's committed order assertions, which are green |
| V8 (nothing else moved) | **DONE, independently re-run** | 788 passed / 0 failed / exit 0; `fmt --check` clean; `clippy --all-targets --release -D warnings` exit 0, 0 warnings. See the arithmetic note below |
| V9 (Perl parity CI job) | **DEFERRED — CI only** | Not runnable in this sandbox. `Cli::validate()` was edited (step 7), so this is the one item to watch on the PR, exactly as §12 says |

**Arithmetic note on V8.** V8's expectation cell reads "≥ 782 + new tests passing". The observed total is
788 with **+7** net test functions, i.e. a base of 781 — one below the plan's estimate. The "782" was the
plan author's own rough figure (the same cell also estimates "~283 integration results" where the actual
integration count is 294); §12's "781 base + 7 new" is self-consistent with the measurement. V8's
substance — 0 failures, fmt and clippy clean — is satisfied. Not a gap.

## Coverage ledger — B3 edge cases (6 rows)

| Case | Plan's handling | Status |
|---|---|---|
| `--library <preset>` on SE | Every preset sets `clip_r1 ≠ 0`, so the Read 1 term keeps the refusal; step 11 tests all five | DONE — five-preset loop green |
| `--rrbs` directional auto-clip | Set in `setup_trimming`, after the guard, so it cannot reach it — **and therefore also escapes it, a separate under-refusal tracked as #459, deliberately not fixed here. "No change in this commit"** | **Row text verified verbatim; matches Q-5 and §12.** Code confirms no change: nothing in the diff touches the `--rrbs`/`setup_trimming` clip path |
| `--paired` + one interleaved uBAM | `paired == true` so R2 counts; the input-format term excludes the shape anyway | DONE — `any_effective(cli.paired)` receives `true`; input-format term untouched |
| `--clip_R2 0` on SE | Runs, as any inert Read 2 value now does | DONE — `applies(false)` is false for `ClipR2`, so the value is not consulted; the row holds whatever clap yields |
| `--clip_R1 0` (any run) | Still refused, still wrong; out of scope (Q-1) | DONE as specified — code unchanged on this path. **Out of scope by the owner's decision; not reported as a gap** |
| `--clock` / `--implicon` + clip flags | Refused for uBAM output before the guard | DONE — `cli.rs:631-641` intact and unmodified |

## Coverage ledger — the six-site fold-in

The plan's decisive argument for widening scope was that v1 rewrote the one site no test pins while
leaving the two best-pinned ones alone. All six are folded in:

| Site | Before | After | Status |
|---|---|---|---|
| 1 `flag_summary` | `if read_2 && !paired { continue }` over a string/bool tuple table | `if !flag.applies(paired)` over `ClipFlag` variants, `self.value(flag)` | DONE |
| 2 inert-flag warning ladder | raw tuple table with a nested `if all_four` per entry | `ClipFlag` list + `user_clips().value(flag)` + one `read_2()` guard | DONE |
| 3 `any_set()` → `any_effective(paired)` | no read awareness — **the defect** | `ALL.iter().any(\|f\| f.applies(paired) && value(f).is_some())` | DONE (the only behaviour change) |
| 4 `ClipOverride.flag` + `.read_2` | `&'static str` + four hand-written bools at `pick` sites | single `ClipFlag` field; wrong-bool now unrepresentable | DONE |
| 5 `main.rs:1128` filter | `\|o\| cli.paired \|\| !o.read_2` | `\|o\| o.flag.applies(cli.paired)` | DONE |
| 6 `report.rs:287` filter | `\|o\| config.paired \|\| !o.read_2` | `\|o\| o.flag.applies(config.paired)` | DONE |

**The predicate is now spelled exactly once.** `git grep -n "paired || !" -- src` returns a single hit,
`library.rs:65` inside `applies`. All four former spellings route through it.

## Observations (not gaps, no action implied)

1. **`ClipFlag::ALL`'s order coincides with the warning ladder's** (`R1, 3'R1, R2, 3'R2`). That is what
   the plan's own `ALL` snippet specified, and both display sites do keep their own explicit variant
   lists as required — but the coincidence means a future edit that switched the ladder to iterate `ALL`
   would look output-neutral. `all_four_inert_flags_warn_in_clipping_order` would catch a change to the
   ladder's order itself, and `flag_summary`'s order (which differs from `ALL`) is pinned 3–4 ways.
2. **Step 12's test inspects the first record's name**, where the plan says "any read name". The fixture
   is a one-record FASTQ, so first == any here.
3. **The removed drift test's incidental `overrides.len() == 4` assertion is covered elsewhere** —
   `integration_library.rs:313` asserts it, and `:323` additionally pins each JSON `flag` string, which
   is what makes step 9's byte-identity claim CI-held rather than measured-once.
4. **`ClipOverride.flag` no longer carries a per-field doc comment** (the old `&'static str` field had
   one). The plan's Signature snippet shows the field without one.
5. **V1, V2, V4–V7 are pre-merge measurements** recorded in §12 and not reproducible from the tree. I
   verified what is checkable without editing the worktree: the predicates are pristine, no snapshot or
   backup file was left behind, and each mutation's expected red test does in fact assert the direction
   §12 claims. The plan itself converted the highest-value part of V7 into committed assertions (step
   15), and those are green.

## Verdict

**COMPLETE.** All 18 outline steps landed; all 12 truth-table rows have a test that was confirmed to run
by name and whose assertions match the row's expected column; all 10 validation items are accounted for;
all 6 B3 edge cases are handled as written; all 6 read-ownership sites are folded into one `ClipFlag`.
Five deviations D1–D5 are documented in §12, so the one that touches an outline step (D1, the two-line
comment) counts as DEVIATED rather than PARTIAL. Independent re-runs agree with §12: 788 tests passing,
0 failures, fmt clean, clippy silent.

Nothing is outstanding for the implementer. The single unverified ledger item is **V9** (the
`Validate vs Perl TrimGalore` CI job), which cannot be run locally; `Cli::validate()` was edited, so it
should be checked on the PR.

Two known defects were deliberately excluded by the owner and are **not** counted as gaps: `--clip_R1 0`
over-refusing (Q-1, unfiled) and the `--rrbs` under-refusal (#459). B3's `--rrbs` row was read verbatim
and does state the under-refusal, the issue number, and "No change in this commit", as Q-5 and §12 claim.
