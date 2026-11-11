# Plan Coverage Report

**Mode:** B, adapted — code audited against the design plan directly (no `IMPL.md` stage ran)
**Plan:** `plans/10202026_polya-behind-polyg/PLAN.md` (rev 2)
**Code:** `fix/polya-behind-polyg` @ `869912c`, three commits off `dev` (`a71d9ff`)
**Date:** 2026-10-20
**Verdict:** INCOMPLETE — 2 items unresolved (both one-line edits in `CHANGELOG.md`)

Everything functional in rev 2 landed. The two open items are documentation obligations the plan
stated explicitly and the entry does not meet. Four of the six validation checks were independently
reproduced from scratch against the pre-change binary; the figures match the commit messages exactly.

## Summary

- Total items: 34 (14 implementation-outline, 10 edge cases, 6 validation checks, Q1, rev-1 residue,
  B2-scope check, the Integration-section CHANGELOG obligation)
- DONE: 30
- PARTIAL: 1
- MISSING: 1
- DEVIATED: 2 (both trivial, one documented in its commit message)

## Coverage ledger

### Commit 1 — reorder the poly stages (`b3ad41a`)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1.1 | Move poly-G block above poly-A; renumber headers so poly-G is 2.7, poly-A 2.8; keep the #321 explanation verbatim | Outline C1.1 | DONE | `src/trimmer.rs:195-218` (poly-G), `:220-240` (poly-A). The #321 paragraph is character-identical across the move. |
| 1.2 | One line above the poly-G block giving the invariant, current-state only | Outline C1.2 | DONE | `src/trimmer.rs:197-199` — one sentence over two wrapped lines, no history. Within CLAUDE.md's two-line maximum. |
| 1.3 | Update the `trim_read` doc comment, whose processing-order list omits the poly stages | Outline C1.3 | DONE | `src/trimmer.rs:85-86` adds `2.7. Poly-G trim (if --poly_g)` / `2.8. Poly-A/poly-T trim (if --poly_a)`. |
| 1.4 | `src/quality.rs`: reword `test_poly_g_after_poly_a_removal`'s comment; rename to `test_poly_g_tail` | Outline C1.4 | DEVIATED | Renamed correctly (`src/quality.rs:517`). The comment was **deleted** rather than reworded. Intent is met — the misstated order is gone, and the order is now recorded in the `trim_read` doc comment, the inline invariant, and six `trim_read` tests, as the plan required ("do not leave the stage order recorded nowhere"). |

### Commit 2 — tests (`23cb2a6`)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 2.1 | `trim_read` unit tests via `poly_g_config()` with `poly_a` set: tail removed behind a G-run; mirror residual pinned; `--no_poly_g` pins Q1 | Outline C2.1 | DONE | Six tests at `src/trimmer.rs:1178-1274`, all three named cases present plus a flush-tail case, a no-tail two-sided control, and an R2 case. Helper `poly_ag_config()` wraps `poly_g_config()` as specified. |
| 2.2 | `src/quality.rs` unit tests for the window the existing list skips: `3 ≤ n ≤ 10` | Outline C2.2 | DONE | `test_poly_a_short_reads` (`src/quality.rs:501-514`) sweeps `n = 3..=10` on all-A, all-T, no-tail and no-head inputs — four assertions per `n`, both scan directions. |
| 2.3 | New `tests/integration_poly.rs`, synthetic FASTQ via `std::fs::write`, binary run with `--poly_a --poly_g -a AGATCGGAAGAGC --dont_gzip` | Outline C2.3 | DONE | 237 lines, 4 tests, all pass. The plan named `tempfile::TempDir`; the implementation uses `std::env::temp_dir()`. That is **not** a deviation from the repo — the cited reference (`tests/integration_partial_output.rs:525`) and all 15 other integration test files use `std::env::temp_dir()`; none uses `tempfile::TempDir`. The plan named the crate wrongly, the code followed the real pattern. |
| 2.4 | Every assertion two-sided — exact expected length, plus a no-`--poly_a` control | Outline C2.4 | DONE | Every integration case asserts `seq == INSERT` (exact 60 bp), never "does not end in A". `reads_without_poly_a_are_untouched` is the control and asserts the untrimmed read verbatim. `test_stage_order_does_not_touch_reads_without_tails` is the unit-level equivalent. |
| 2.5 | `--paired --poly_a` R2 byte-identity as a committed test, not a manual before/after | Outline C2.5 | DEVIATED (documented) | `paired_r2_5prime_clip_is_independent_of_3prime_trimming` pins a different property: every R2 record from the `--poly_g` run is a prefix of the same record from the `--no_poly_g` run, with a `shortened > 0` guard so it cannot pass vacuously. Byte-identity *across the change* cannot be expressed by a committed test (it needs the old binary). The substitution is stated plainly in the commit message. I reproduced the actual byte-identity out of band — see Validation 4. |

### Commit 3 — docs and CHANGELOG (`869912c`)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 3.1 | `migration.md` table row `\| Poly-A trimming \| None \| Built-in --poly_a \|` is false | Outline C3.1 | DONE | Now `` `--polyA`, a separate second pass `` → `` Built-in `--poly_a` stage ``. |
| 3.2 | Rewrite the `-a 'A{N}'` note with the three corrected Perl facts | Outline C3.2 | DONE | All three present and correct: `AAAAAAAAAA` (10 A's) auto-detected with `A{20}` only the no-match fallback and the `--help` text flagged as disagreeing with its own code; no `-q`/`--length`, documented as the second of two runs; the `32:A:` prefix / `_PolyA:32` suffix read-ID rewriting with the consequence for ID-grepping pipelines. The `-a 'A{N}'` material is retained under its own subheading. |
| 3.3 | `modes/poly.md:6` claim is wrong for poly-A; add stage order, the 2-colour interaction, and a 3'seq row | Outline C3.3 | DONE | Lede reworded; new `### Stage order` section with both consequences (`--no_poly_g`, the mirror arrangement) and the `--2colour 20` recommendation; `\| 3'seq on a 2-colour instrument \|` row added to "When to use which". |
| 3.4 | `src/cli.rs` `--poly_a` help text — extend "Runs after adapter trimming" to name poly-G | Outline C3.4 | DONE | `src/cli.rs:365-367`; confirmed in rendered `--help`. |
| 3.5 | `CHANGELOG.md` — one sentence, repo format, **crediting @marenmsk** | Outline C3.5 | PARTIAL | Entry present, repo format, one sentence after the issue link, with a working `#stage-order` deep link. **No credit to @marenmsk** — `@marenmsk` appears nowhere in the tree. The repo does credit reporters in exactly this position (`CHANGELOG.md:791` "Reported by @K81ta in [#321]", and @an-altosian throughout). |
| I1 | "Report attribution shifts for `--poly_a` users on 2-colour data … worth one line in the CHANGELOG entry" | Integration §, last bullet | MISSING | The entry says nothing about it. This is a real user-visible change, measured below: on `polyAT_R1` with both flags forced, R1 poly-A reads 4,229 → 4,239 and R1 poly-G reads 217 → 215 (888 → 882 bases). |

### Edge cases (plan's ten-row table)

| # | Case | Status | Test or guard |
|---|------|--------|---------------|
| E1 | Read < 3 bp | DONE | `n < 3` guard at `src/quality.rs:313-315`, unchanged. No test at `n ∈ {1,2}` — the plan asked only for the guard, and scoped the new tests to `3 ≤ n ≤ 10`. |
| E2 | Empty record | DONE | `!record.is_empty()` present at both call sites (`src/trimmer.rs:208`, `:222`); `test_poly_a_empty` covers `n = 0` at the function level. |
| E3 | Read is all A, or all G | DONE | Pre-existing `test_poly_a_all_a`, `test_poly_g_all_g`; `test_poly_a_short_reads` adds all-A and all-T for `n = 3..10`. |
| E4 | Tail < 3 bp | DONE | Pre-existing `test_poly_a_short_tail_ignored`, `test_poly_g_short_tail_ignored`, `test_poly_c_short_head_ignored`. |
| E5 | `--no_poly_g --poly_a` on 2-colour data (Q1) | DONE | `test_no_poly_g_leaves_tails_behind_a_g_run` asserts 100 bp kept and both counters at 0. |
| E6 | `[insert][polyG][polyA]` mirror residual | DONE | `test_poly_g_behind_poly_a_survives` asserts exactly 75 bp, `poly_a_trimmed == 40`, `poly_g_trimmed == 0`, tail ends in G. |
| E7 | `--poly_a` on R2 | DONE | `test_poly_a_r2_head_is_unaffected_by_the_stage_order` (unit) plus the paired integration case. |
| E8 | `--rrbs` | DONE | RRBS is still step 2.5, ahead of both poly blocks. No test — the plan asked for none. I confirmed `--rrbs --paired` output is md5-identical to `dev`. |
| E9 | `--output-format ubam --poly_a` | DONE | `run_single_end_to_bam` (`src/trimmer.rs:639`) and `run_paired_end_to_bam` (`:709`) both call the shared `trim_read` at `:650` and `:735-736`, so the reorder is inherited. No test — the plan asked for none. |
| E10 | `--clump_only --poly_a` | DONE | Still rejected at `src/cli.rs:863-865` (`--clump_only does not trim; --polyA is not compatible`). |

### Q1 — both halves

**DONE.** Rev 2 required `--no_poly_g --poly_a` to be "documented and pinned with a test".
Documented at `docs/src/content/docs/modes/poly.md:49` ("With `--no_poly_g`, a poly-A tail behind a
G-run is not removed. The opt-out is honoured, and the G-run still blocks the scan.") and pinned by
`test_no_poly_g_leaves_tails_behind_a_g_run`.

### Rev-1 / B2 scope check

**Clean, both directions.**

- **Rev 2 never names B2 as in-scope.** All seven mentions (`PLAN.md:10, 15, 106, 198, 268, 286, 301`)
  are revision history, the Non-goals list, the "Withdrawn" section, a note that the existing suite
  gives no signal either way, the Validation-3 rationale, the Follow-up successor, and the
  self-review. No outline item, behavior clause, edge case or validation row depends on it.
- **No half-carried rev-1 code.** `homopolymer_trim_index` at HEAD is byte-identical to `dev` — both
  `src/quality.rs` hunks on this branch are inside `mod tests`. There is no offset loop, no bound
  constant, and no change to either gate (`errors * 5 <= n - i` / `errors * 5 <= i + 1`). The plan's
  "no change to `src/quality.rs`'s algorithm" holds literally.

## Validation table — independently reproduced

I re-ran four of the six checks from scratch against the pre-change binary at
`/tmp/claude-501/tg_dev_b1/target/release/trim_galore` rather than copying `PROGRESS.md`'s results.
That binary reports its hash as `unknown` (built outside a work tree), so I confirmed it behaviourally
instead: on the Validation-1 inputs it returns the pre-change lengths, which only a pre-reorder build
does.

| # | Plan's check | Reproduced? | Result |
|---|---|---|---|
| 1 | The reported case is fixed end-to-end | **Yes, exactly** | Synthetic set through both binaries with `--quality 20 --phred33 --polyA --poly_g -a AGATCGGAAGAGC`. Old: `[100, 100, 80]`. New: `[60, 40, 60]` — the plan's expected values to the base, for the plan's own constructions (`insert60+A40+G15` → 60, `insert40+A60+G20` → 40, `insert60+A20+G20` → 60). Also confirms the `--polyA` alias reaches the new path. |
| 2 | Nothing else moves | **Partly** | 516 lib tests pass — exactly the plan's 509 + the 7 new ones (6 in `trimmer.rs`, 1 in `quality.rs`). `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release -- -D warnings` exits 0, but in 0.30 s, i.e. a replay of a cached run rather than a fresh compile; I did not force a rebuild because two reviewers share this tree. **Integration tests total 307 across 17 binaries, not the 219 `PROGRESS.md` claims** — all pass, so this is a tally error in the progress file, not a gap (the 88 difference is the size of three of the binaries). The Perl `validation` job is CI-only and was not run locally. |
| 3 | No over-trim, measured as a differential | **Yes, exactly** | Both binaries, `--quality 20 --phred33 --poly_a --poly_g --dont_gzip --length 0 -a AGATCGGAAGAGC`, compared per record by index: `illumina_10K` 10 shorter / 2 longer, worst cut 16 bp; `illumina10K_with_polyA` (10,900 reads) 10 / 2, worst 16 bp; `BS-seq_10K_R1` 7 / 9, worst 15 bp; `SRR24766921_RRBS_R1` 15 / 5, worst 8 bp. Matches commit 2's message figure for figure. Plan's gate was "≤ ~20 reads per fixture" — met with margin. |
| 4 | R2 is byte-identical | **Yes, and extended** | `--paired --poly_a` alone on `polyAT_R1/R2` is vacuous: poly-G auto-detection does not fire on that fixture (0 of 13,000), so both mates are trivially identical. Forcing `--poly_a --poly_g` gives the real answer: **R2 md5-identical** (zero differing lines), R1 differs in exactly 12 records. `PROGRESS.md`'s broader claim also reproduces — default-flag output is md5-identical to `dev` for SE, PE, `--rrbs --paired`, `--clock`, `--hardtrim5 30`, and for `--poly_g` alone (the `poly_a = false` no-op). |
| 5 | The documented workaround still works | **Yes** | `two_colour_quality_trimming_also_exposes_the_tail` passes, and manually `--2colour 20 --poly_a` gives `[60, 40, 60]` on **both** binaries — the workaround was unaffected, as the plan assumed. |
| 6 | The mirror residual is where we think it is | **Yes** | `test_poly_g_behind_poly_a_survives` asserts 75 bp and passes. |
| 7 | (`PROGRESS.md` extra) Docs build, 32 pages, `#stage-order` anchor | **Yes, from artifacts** | `docs/dist/` from the 22:23 build holds 31 `index.html` pages plus `404.html` = 32, and `docs/dist/modes/poly/index.html` contains `id="stage-order"`, so the CHANGELOG's deep link resolves. I did not rebuild — that writes into the shared tree. |

## Gaps (detail)

### Item 3.5: CHANGELOG does not credit the reporter

**Expected:** "`CHANGELOG.md` — one sentence, repo format, crediting @marenmsk."
**Found:** The entry at `CHANGELOG.md:8-13`. Correct format, one sentence, working deep link. No
mention of @marenmsk anywhere in the tree (`git grep marenmsk` is empty).
**Gap:** Add the credit. The house form in this file is `Reported by @handle in [#NNN](url)`
(`CHANGELOG.md:791`).

### Item I1: the report-attribution shift is not in the CHANGELOG

**Expected:** Integration §: "Report attribution shifts for `--poly_a` users on 2-colour data…
Expected, and worth one line in the CHANGELOG entry."
**Found:** Nothing. The entry describes the trimming change only.
**Gap:** One clause saying that on 2-colour data some reads that previously reported poly-G only will
now also report poly-A. Measured on `polyAT_R1` with both flags forced: R1 poly-A reads 4,229 → 4,239,
R1 poly-G reads 217 → 215, R1 poly-G bases 888 → 882. Note the plan predicted bases "counted as
poly-G will be counted as poly-G still"; 6 bases in fact move out of poly-G attribution, which is the
mirror residual (2 reads where poly-G is now the blocked stage) rather than a surprise.

## Test verification

| Test | File | Status |
|---|---|---|
| `test_poly_a_tail_behind_poly_g_is_trimmed` (3 cases) | `src/trimmer.rs:1187` | PASS |
| `test_poly_a_tail_flush_with_3prime_end_is_trimmed` | `src/trimmer.rs:1206` | PASS |
| `test_poly_g_behind_poly_a_survives` | `src/trimmer.rs:1217` | PASS |
| `test_no_poly_g_leaves_tails_behind_a_g_run` | `src/trimmer.rs:1232` | PASS |
| `test_stage_order_does_not_touch_reads_without_tails` | `src/trimmer.rs:1246` | PASS |
| `test_poly_a_r2_head_is_unaffected_by_the_stage_order` | `src/trimmer.rs:1257` | PASS |
| `test_poly_a_short_reads` (`3 ≤ n ≤ 10`) | `src/quality.rs:501` | PASS |
| `test_poly_g_tail` (renamed) | `src/quality.rs:517` | PASS |
| `poly_a_tail_behind_poly_g_is_trimmed` | `tests/integration_poly.rs:52` | PASS |
| `two_colour_quality_trimming_also_exposes_the_tail` | `tests/integration_poly.rs:95` | PASS |
| `reads_without_poly_a_are_untouched` | `tests/integration_poly.rs:131` | PASS |
| `paired_r2_5prime_clip_is_independent_of_3prime_trimming` | `tests/integration_poly.rs:178` | PASS |

Suite totals: 516 lib (was 509) + 307 integration + 0 doc-tests, all green. `cargo fmt` clean.
`cargo clippy -D warnings` exits 0 from cache.

## Notes on plan accuracy (not gaps)

Two statements in rev 2 are narrower or looser than the code, and both bear on arguments the plan
rested conclusions on. Neither changes a verdict, because I verified the conclusions empirically.

- **"The sole `true` source is `src/main.rs:1294 poly_a: cli.poly_a`"** (Behavior §, B1 scope). There
  are four: `src/main.rs:1294, 1528, 2318, 2697`. The safety claim the plan drew from this still holds
  — `--poly_g`-only and default output are md5-identical to `dev` — but it holds because the block is
  flag-gated, not because there is one construction site.
- **"`tempfile::TempDir` (pattern at `tests/integration_partial_output.rs:525`)"** (Outline C2.3). The
  cited file does not use `tempfile::TempDir`, and neither does any other integration test. The
  implementation followed the actual pattern.

One incidental observation, outside this plan's ledger: `--poly_g`'s own help text still says it trims
"poly-C heads from the 5' end of Read 2", which #321 established is not what the code does. Untouched
by this branch and not in scope here; noted in case the code reviewers want it.

## Verdict

**INCOMPLETE — 2 items unresolved**, both single-line edits to the `CHANGELOG.md` entry at lines 8-13:

1. **Credit @marenmsk** (outline item C3.5, stated verbatim in the plan). House form in this file:
   `Reported by @marenmsk in [#504](https://github.com/FelixKrueger/TrimGalore/issues/504)`.
2. **Name the report-attribution shift** (Integration §, "worth one line in the CHANGELOG entry"). One
   clause: on 2-colour data, reads that previously reported poly-G only will now also report poly-A.

Nothing else is outstanding. The reorder, all six `trim_read` tests, the `3 ≤ n ≤ 10` unit coverage,
the four-case integration file, all five docs edits, every one of the ten edge cases, both halves of
Q1, and all six validation checks are in place; four of the six validations reproduce to the figure
against a pre-change binary. The two DEVIATED items (a deleted rather than reworded test comment, and
a prefix-invariance test standing in for a cross-build byte-identity comparison a committed test
cannot make) are equivalent in effect and need no action. No rev-1 or B2 residue exists in either the
plan or the code.
