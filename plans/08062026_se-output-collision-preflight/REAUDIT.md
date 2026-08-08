# Plan Coverage Re-Audit — D7–D13 findings vs. current `dev`

**Mode:** B (code vs. plan + the D7–D13 findings ledger)
**Plan(s):** `PLAN.md` §13 (deviations D1–D13), `CODE_REVIEW_D13.md` (findings 1–4),
`PROGRESS.md` (Outstanding / Left undone)
**Date:** 2026-08-08
**Tree audited:** branch `dev` @ `ce13761`
**Verdict:** **INCOMPLETE — 3 items unresolved** (1 MISSING, 2 PARTIAL; none behavioural)

## Why this re-audit exists

`COVERAGE.md`, `COVERAGE_round2.md` and `COVERAGE_round3.md` all ran against the tree
*before* the post-code-review round; all three report 521/530 tests. `PROGRESS.md` records
the caveat verbatim: "no independent coverage audit has run against the final state". This
file is that audit. Its ledger is scoped to what the three historical audits could not have
seen — **D7–D13** in `PLAN.md` §13, the **four findings** in `CODE_REVIEW_D13.md`, and the
items `PROGRESS.md` listed as Outstanding / Left undone — verified against the code as it
now stands, seven merged PRs later.

The historical `COVERAGE*.md` files are left untouched.

**Supersession is in scope.** Five PRs merged after #385 touch this exact surface, and
several D7–D13-era requirements were satisfied by them rather than by #385 itself:

| PR | Commit | What it changed on this surface |
|---|---|---|
| #385 | `533582a` | the pre-flight itself (this plan) |
| #390 (#384) | `c4599f5` | case-insensitive FASTQ extension matching |
| #392 (#388) | `ad2639a` | paired report paths join the pre-flight |
| #393 (#386) | `03d793f` | `--hardtrim5` + `--hardtrim3` rejected |
| #394 (#387) | `7ea2741` | CWD-output naming documented on three mode pages |
| #395 (#389) | `0f65048` | `--passthrough` alias check: identity vs. collision split |
| #396 (#391) | `ce13761` | clump reports join the pre-flight; property tests re-keyed |

## Summary

- Total items: **20**
- DONE: **14**
- PARTIAL: **2**
- MISSING: **1**
- DEVIATED: **3** (2 documented, 1 undocumented but mutation-verified as covered)

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `collision_key` folds `..` lexically: absolutise, fold `ParentDir` against the stack, POSIX `/..` == `/`, preserve an unfoldable leading `..`; A8 narrowed to symlinks only | D7 | DONE | `io.rs:52-70`. Both directions pinned: `preflight_rejects_dotdot_alias` (`io.rs:1135`) and `preflight_keeps_unfoldable_paths_distinct` (`io.rs:1143`) |
| 2 | One definition of file identity — `collision_key` `pub`, and all four `cli.rs` checks use it instead of raw `==` / `norm_path` | D8 | DEVIATED | Superseded by D13 and #389. Three identity checks use `path_identity_key` (case-preserved); the `--passthrough` check uses **both** keys. Deviation documented in D13 itself. The behavioural goal survives: `--paired ./a_R1.fq a_R1.fq` is still refused, because both keys share `lexical_normalise` |
| 3 | Secondary outputs checked against inputs: reports + `--demux` per-barcode files join the candidate list; `guarded_inputs` gains the barcode file; `demux_output_paths` shares the writer's helpers | D9 (the unanimous High) | DONE | `planned_secondary_outputs` `main.rs:91-112`; `guarded_inputs` includes `cli.demux` `main.rs:77-79`; `demux::demux_output_paths` `demux.rs:143`; pinned by `demux_writes_exactly_the_planned_paths` (`tests/integration_output_collision.rs:806`) |
| 4 | CWD hint covers all four CWD-naming modes and **replaces** the generic advice; `HARDTRIM_HINT` → `CWD_OUTPUT_HINT`; `run_specialty_paired` takes a `hint` | D10 | DONE | `CWD_OUTPUT_HINT` `main.rs:65`, passed at `main.rs:425`, `:456` (hardtrim FASTQ + uBAM) and `:486`, `:501` (clock, implicon). Replacement, not append: `io.rs:122` `hint.unwrap_or(GENERIC_ADVICE)`, asserted in both directions by `preflight_appends_hint_only_when_given` (`io.rs:1149`) |
| 5 | `--clump_only --paired` passes `None` (its namer uses `input.parent()`) | D10 | DEVIATED | Superseded by #391: now `Some(PAIRED_REPORT_HINT)` at `main.rs:552`. Documented — the code comment names #391 and the reason (reports carry no `_clumped_N` discriminator), and #391 has its own plan directory |
| 6 | The V5 fixture can fail: add `d/same.fastq.gz` + `e/same.fastq.gz` so the A2 implication is not a tautology | D11 | DONE | `io.rs:1232-1241`, with the reason stated in-test. Strengthened by #391: a fold-equal twin `d/SAME.fastq.gz` plus comparison on `collision_key` rather than `PathBuf` equality, and the three clump primary namers added |
| 7 | `demux_base_name` lifted out of `demultiplex`'s doc-comment | D12 | DONE | `demux.rs:109-112` carries its own doc block; `demultiplex` has its own at `:167-169` |
| 8 | Stale key descriptions deleted at the paired and `run_specialty_paired` sites | D12 | DONE | `run_specialty_paired` now reads "see `io::collision_key` for the key" (`main.rs:2541`); the paired site carries no key description |
| 9 | `assert_rejected_cleanly`'s `unwrap_or_default()` replaced with a panic so a `read_dir` failure cannot pass vacuously | D12 | DONE | `tests/integration_output_collision.rs:100-110`, with the rationale in the comment |
| 10 | CHANGELOG: defect-2 overclaim corrected, plus entries for the CWD hint, the single identity definition, and A1's case-fold trade-off | D12 | PARTIAL | All four entries exist (`CHANGELOG.md:109-121`, `:200-204`, `:206-211`, `:213-217`). The "One definition of 'the same file'" entry was **not revised after D13** — see Gap 1 |
| 11 | The identity key is split: `path_identity_key` (case-preserved) → the three `cli.rs` input-identity checks; `collision_key` (case-folded) → the pre-flight and the `--passthrough` alias check | D13 | DONE | `io.rs:75`, `:81`; `cli.rs:537`, `:549-550`, `:635`/`:638`; pre-flight `io.rs:99`, `:104`. Pinned by `identity_key_is_case_sensitive_but_collision_key_is_not` (`io.rs:1178`), whose doc block records the regression |
| 12 | Finding 1 — the APFS false-negative: add the per-pair report paths to the paired pre-flight candidate list | `CODE_REVIEW_D13` §Findings 1 | DONE | Superseded by **#388 / PR #392** (`ad2639a`), which implemented exactly the recommended fix. Paired FASTQ `main.rs:824-831`, paired uBAM `main.rs:1938-1941`; skipped under `--no_report_file` to match the writers. Tests `paired_report_candidates_do_not_over_reject` and `paired_ubam_rejects_same_filename_r1_r2_into_shared_output_dir` (both PASS). Disclosed at `CHANGELOG.md:126-141`, including the two behaviour changes |
| 13 | Finding 2 — the `--passthrough` alias check asks an identity question on a collision key, and false-positives on a case-sensitive filesystem | `CODE_REVIEW_D13` §Findings 2 | DONE | Superseded by **#389 / PR #395** (`0f65048`), which took the second of the two options the finding offered and then some: `cli.rs:931-963` now tries `path_identity_key` first and `collision_key` second, with a distinct, accurate message per subcase, and the comment states the intent ("case-folded on purpose, not an identity check"). Test `test_passthrough_rejects_case_variant_of_input` (`cli.rs:1935`); `CHANGELOG.md:145-152` |
| 14 | Finding 3 — the D8 regression has no unit-level guard; a `Cli::parse_from(…).validate()` test over the four case-only names would have caught it | `CODE_REVIEW_D13` §Findings 3 | **MISSING** | See Gap 2. No such test exists anywhere in `src/` or `tests/`; only `.github/workflows/ci.yml:625-630` covers the wiring |
| 15 | Finding 4a — `both_keys_normalise_spelling` does not pin absolutisation | `CODE_REVIEW_D13` §Findings 4 | DEVIATED | The literal ask was not implemented — the test is unchanged and still all-relative (`io.rs:1194-1206`). **Mutation-verified** (below): the finding's claim is correct, *and* the same mutation is caught by `preflight_rejects_absolute_versus_relative_alias` (`io.rs:1091`), which drives the shared `lexical_normalise`. The regression is guarded; the deviation is undocumented |
| 16 | Finding 4b — `norm_path`'s doc block never mentions `path_identity_key` or the identity/collision distinction | `CODE_REVIEW_D13` §Findings 4 | PARTIAL | The block was rewritten and now names `collision_key` and #389 (`io.rs:38-44`), but still not `path_identity_key`. Residual tracked in open issue **#399** (demote or inline `norm_path`), which would remove the surface entirely. See Gap 3 |
| 17 | Symlink aliasing — the sole remaining A8 residual — disclosed in the CHANGELOG with its cost | `PROGRESS.md` Outstanding | DONE | `CHANGELOG.md:106-108`: "A path reached through a **symlink** still aliases undetected; that is the one remaining case, and it costs the same silent read loss this entry describes." Also stated on `lexical_normalise` (`io.rs:51`) |
| 18 | `--fastqc_args "-o DIR"` — A2's named exception, documented rather than covered | `PROGRESS.md` Outstanding | DONE | Recorded in the V5 test's own doc block as a known residual, explicitly not asserted (`io.rs:1221-1226`) |
| 19 | §11 Open 3 — `--hardtrim5 N --hardtrim3 M` silently ignores the 3′ trim; "wants its own issue" | `PROGRESS.md` Left undone | DONE | Superseded by **#386 / PR #393** (`03d793f`): rejected at `cli.rs:968`, with the §3.4a family-precedent message shape; `CHANGELOG.md:162-166` |
| 20 | §11 Open 4 — a `modes/hardtrim.md` sentence about CWD output; "wants its own issue" | `PROGRESS.md` Left undone | DONE | Superseded by **#387 / PR #394** (`7ea2741`), which covered all three CWD modes rather than only hardtrim: `docs/src/content/docs/modes/hardtrim.md:42`, `clock.md:49`, `implicon.md:32` |

## Gaps (detail)

### Gap 1 — Item 10: the CHANGELOG's identity entry contradicts D13

**Expected:** D12 added a CHANGELOG entry for "the single identity definition". D13 then
split that single definition in two, and #389 refined the `--passthrough` half again. The
entry should describe the shipped behaviour.

**Found:** `CHANGELOG.md:206-211` still reads:

> **One definition of "the same file".** `--paired`'s R1≠R2 check, its duplicate-pair
> check, the duplicate-input check and the `--passthrough` alias check all compared paths
> more weakly than the collision pre-flight did … **All four now use the same key as the
> pre-flight.**

After D13, that final sentence is not accurate for any of the four. The pre-flight uses
`collision_key` (case-folded); the three input-identity checks use `path_identity_key`
(case-preserved), and the `--passthrough` check uses both. Nothing in the CHANGELOG
discloses the split — `grep` for `path_identity_key` and "identity key" over `CHANGELOG.md`
returns nothing.

**Gap:** the entry needs its last sentence corrected to describe the two keys and why they
differ (case-folding is right for output paths, wrong for input identity — D13's own
reasoning). No code change; documentation accuracy only. The bug the entry claims to fix
(`--paired ./a_R1.fq a_R1.fq` exiting 0) *is* fixed, because both keys share
`lexical_normalise`, so the entry is stale rather than false about the outcome.

### Gap 2 — Item 14: the D8/D13 regression still has no unit-level guard

**Expected (verbatim from `CODE_REVIEW_D13.md` finding 3):** a test in the existing
`Cli::parse_from(…).validate()` style — no files on disk, so it runs identically on macOS
and Linux —

```rust
let cli = Cli::parse_from(["trim_galore", "--paired",
    "Sample_R1.fastq.gz", "Sample_R2.fastq.gz", "SAMPLE_R1.fastq.gz", "SAMPLE_R2.fastq.gz"]);
assert!(cli.validate().is_ok(), "case-only variants are distinct inputs");
```

**Found:** nothing. `grep -rn "SAMPLE_R1\|SAMPLE_R2" src/ tests/` returns only
`io.rs:1174-1180` (the *key*-level test, which is precisely the insufficiency the finding
named: "the two new tests pin the keys; what broke was the wiring in `cli.rs`"),
`specialty.rs:826`/`:840` (unrelated `.FASTQ.GZ` naming tests from #384), and
`.github/workflows/ci.yml:625-630` — the #216 validation guard, still the only thing
covering the wiring, exactly as the finding described.

**Gap:** the four-line test above, in `cli.rs`'s test module. Not tracked by any open issue
(#397–#400 cover unrelated follow-ups), so it is currently neither implemented nor deferred
on the record.

### Gap 3 — Item 16: `norm_path`'s doc block still omits the split

**Expected:** one line in `norm_path`'s doc block naming `path_identity_key` and the
identity-vs-collision distinction, so the split is discoverable from the function a reader
hits first.

**Found:** the block was rewritten during #389 and now reads "Case-folded (ASCII lowercase)
string view of a path — the folding component of `collision_key`, which absolutises first
and is what every collision check uses (issues #216, #383, #389)" (`io.rs:38-44`). It names
`collision_key` but not `path_identity_key`, so a reader landing on `norm_path` still cannot
see that a second, case-preserving key exists. `path_identity_key`'s own block (`io.rs:72-74`)
does explain the rationale.

**Gap:** either the one-line addition, or #399's resolution (making `norm_path` private or
inlining it), which removes the reader-facing surface altogether. Tracked in #399.

### Note on Item 15 (DEVIATED, undocumented)

Finding 4a is the one deviation with no artifact recording the decision, so its evidence is
set out here. I did not modify the repo to test this; I compiled a standalone copy of
`lexical_normalise` with and without the `absolute()` call
(`$TMPDIR/reaudit/norm.rs`) and evaluated both test predicates against each variant:

| variant | `preflight_rejects_absolute_versus_relative_alias` | `both_keys_normalise_spelling` |
|---|---|---|
| with `absolute()` (shipped) | PASS — keys collide, `unwrap_err` succeeds | PASS |
| without `absolute()` (mutant) | **FAIL** — keys differ, pre-flight returns `Ok`, `unwrap_err` panics | PASS |

So the reviewer was right that `both_keys_normalise_spelling` is blind to absolutisation
(all three of its paths reduce to `x.fq` with or without the call), and also that the nit
was never applied. But absolutisation is not unpinned in the tree: the existing #383
regression test catches its removal through the shared helper. The residual is that
`path_identity_key` has no *direct* absolutisation test of its own — it inherits the
coverage via `lexical_normalise`.

## Test verification

Run on `dev` @ `ce13761`. The full suite was verified at 561 green on this commit earlier
today; these are the targeted runs for this surface.

| Command | Result |
|---|---|
| `cargo test --test integration_output_collision` | **45 passed, 0 failed** |
| `cargo test --lib io::tests` | **54 passed, 0 failed** (356 filtered out) |

| Test | File | Status | Ledger item |
|---|---|---|---|
| `preflight_rejects_dotdot_alias` | `src/io.rs:1135` | PASS | 1 |
| `preflight_keeps_unfoldable_paths_distinct` | `src/io.rs:1143` | PASS | 1 |
| `preflight_rejects_dot_slash_alias` | `src/io.rs:1079` | PASS | 1 |
| `preflight_rejects_absolute_versus_relative_alias` | `src/io.rs:1091` | PASS | 1, 15 |
| `preflight_rejects_output_that_aliases_an_input` | `src/io.rs:1107` | PASS | 3 |
| `preflight_rejects_input_alias_across_spellings` | `src/io.rs:1125` | PASS | 3 |
| `demux_writes_exactly_the_planned_paths` | `tests/integration_output_collision.rs:806` | PASS | 3 |
| `preflight_appends_hint_only_when_given` | `src/io.rs:1149` | PASS | 4 |
| `distinct_primary_outputs_imply_distinct_secondary_outputs` | `src/io.rs:1228` | PASS | 6 |
| `primary_output_key_is_coarser_than_secondary_keys` | `src/io.rs:1320` | PASS | 6 |
| `identity_key_is_case_sensitive_but_collision_key_is_not` | `src/io.rs:1178` | PASS | 11 |
| `both_keys_normalise_spelling` | `src/io.rs:1195` | PASS (but see Item 15) | 15 |
| `paired_report_candidates_do_not_over_reject` | `tests/integration_output_collision.rs` | PASS | 12 |
| `paired_ubam_rejects_same_filename_r1_r2_into_shared_output_dir` | `tests/integration_output_collision.rs` | PASS | 12 |
| `test_passthrough_rejects_case_variant_of_input` | `src/cli.rs:1935` | PASS | 13 |
| *case-only four-file `validate()` guard* | — | **MISSING** | 14 |

One bookkeeping correction: `tests/integration_output_collision.rs` holds **45** tests, not
47 — `cargo test --test integration_output_collision` reports "45 passed" and
`grep -c '#\[test\]'` agrees.

## Verdict

**INCOMPLETE — 3 items unresolved.** No behavioural gap: every D7–D13 requirement that
changes what the binary does is implemented and tested, and the two findings the reviewer
flagged as reachable harm (findings 1 and 2) were both closed by later PRs — #392 with the
exact fix the review recommended, #395 with more than it asked for. All 20 ledger items
have code or an artifact behind them except the three below.

What remains, in the order I would take it:

1. **Gap 2 (MISSING, Item 14).** Add the four-line
   `Cli::parse_from([… Sample_R1, Sample_R2, SAMPLE_R1, SAMPLE_R2]).validate().is_ok()`
   test to `cli.rs`'s test module. This is the only unresolved item that guards a *silent*
   failure mode: the D8→D13 regression is currently caught by nothing but the CI validation
   job, so the same overreach could land again and pass a full local `cargo test`. Untracked
   by any issue.
2. **Gap 1 (PARTIAL, Item 10).** Correct the last sentence of `CHANGELOG.md:206-211` — the
   four checks do *not* all use the pre-flight's key any more, and the split is disclosed
   nowhere in the CHANGELOG. Documentation only; the described bug is genuinely fixed.
3. **Gap 3 (PARTIAL, Item 16).** One line in `norm_path`'s doc block naming
   `path_identity_key`, or close it via **#399**, which proposes removing the surface.
   Already tracked.

Item 15 is recorded as DEVIATED rather than a gap: the nit was not applied, but the
absolutisation it worried about is provably pinned by an existing test (mutation evidence
above). Items 2 and 5 are DEVIATED with their deviations documented in D13 and #391
respectively.
