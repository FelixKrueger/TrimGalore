# Plan Coverage Report — second independent re-audit

**Mode:** B (code vs. implementation plan — the plan's own §5 is the implementation ledger; no separate `IMPL.md`, by design)
**Plan(s):** `plans/08062026_se-output-collision-preflight/PLAN.md` (v2, §13 through "Post-audit gap closures")
**Date:** 2026-08-07
**Verdict:** COMPLETE

**Provenance.** This is a *third* coverage pass, run concurrently with and independently of the
re-audit in `COVERAGE.md` (neither read the other while auditing; I read `COVERAGE.md` only after
finishing, to reconcile). `COVERAGE_round1.md` is the first pass. **`COVERAGE.md` remains the
canonical report** — `PROGRESS.md` cites it, it is correct, and this pass concurs with it. This
file exists to record an independent cross-check at finer granularity (67 items vs 25) and one
classification disagreement that I resolved in `COVERAGE.md`'s favour after checking §10 V7.

Audited state: branch `fix/383-output-collision-preflight` off `dev` @ `72624c7`, uncommitted —
**7 modified files** (`ci.yml`, `CHANGELOG.md`, `src/cli.rs`, `src/demux.rs`, `src/io.rs`,
`src/main.rs`, `src/specialty.rs`) plus untracked `tests/integration_output_collision.rs`;
`git diff --stat` +595/−111.

> **The tree moved mid-audit.** This pass began against the tree as the brief described it
> (6 modified files, +503/−99, 520 tests) and the `src/io.rs` + new `src/demux.rs` changes landed
> at 12:54–12:55, closing two of the three gaps this pass had found independently — the same two
> `COVERAGE_round1.md` raised. Every finding below is re-verified against the **current** tree
> with a release binary rebuilt at 13:18 (51 s compile — confirmed a real build, not a stale
> artifact).

## Summary

- Total items: 67
- DONE: 63
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 4 — D1 and D6 documented in §13 and acceptable; the V2.5/V2.6 test-shape change
  undocumented but behaviourally covered by V7 (see Reconciliation)

Gates re-run independently, from scratch, on the current tree:

| Gate | Result |
|---|---|
| `cargo test` (crate root) | **521 passed, 0 failed** — 399 lib + 122 integration across 12 binaries |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean |
| §2.2 reproductions (a)–(f), 8 invocations | all exit 1 with the correct message, nothing written |
| §10 V3.1 positive control | exit 0; `dirA` 40 own/0 foreign, `dirB` 40/0 — zero crosstalk |
| Three new `ci.yml` steps, executed locally | all pass every assertion |
| `ci.yml` YAML parse | valid; `validation` job has **41** steps (§13 still says 40) |
| Over-rejection probes: `--demux`, `--clock` ×2 pairs | both exit 0 — the widened input set does not over-reject |
| `--clock` R1==R2 | keeps its own precise message (D1 verified behaviourally) |
| Cross-pair duplicate `--paired a1 a2 b1 a2` | exit 1 via the pre-flight |

These figures were produced without reference to `COVERAGE.md` and match its independently
obtained figures exactly, including the 41-vs-40 CI step discrepancy.

## Coverage ledger

### §5 — Implementation outline (13 steps)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `io.rs`: `collision_key` + `preflight_output_collisions` after `norm_path`; extend `norm_path`'s doc-comment | Step 1 | DONE | `collision_key` `io.rs:48`, `preflight_output_collisions` `:57`. **Doc-comment sub-clause closed in the 12:55 change**: `:37-38` now names `collision_key`; the stale "`main::run` paired-end" locator is gone |
| 2 | Delete `preflight_collision_bam`; re-point 3 calls; doc-comment rejoins `resolve_clump_layout` | Step 2 | DONE | Absent from `src/`; calls at `:548`, `:585`, `:633`; clumpify-budget doc-comment reattached |
| 3 | Convert the 4 hand-rolled copies, preserving semantics | Step 3 | DONE | `:513` clump_only SE FASTQ; `:734` paired FASTQ (candidates accumulated across **all** `chunks(2)` into one `Vec`, so cross-pair collisions still fire — verified live; `--retain_unpaired` and `--passthrough` extras preserved); `:1825` paired uBAM; `:2430` `run_specialty_paired`. All pass `guarded_inputs(…)`. No `--output-dir` left in any source message (§2.5) |
| 4 | `HardtrimEnd` + `pub` on both namers; 4 internal sites + the bgz test | Step 4 | DONE | Enum `specialty.rs:19`; namers `pub` at `:442`, `:465`; sites `:45`, `:88`, `:139`, `:191`; `bgz_stem_reaches_specialty_output_names` updated |
| 5 | `planned_hardtrim_outputs` | Step 5 | DONE | `main.rs:71-89`, `match cli.output_format` |
| 6 | Guard `--hardtrim5`, first statement in the block | Step 6 | DONE | `main.rs:359-364`, ahead of the writer loop; passes `HARDTRIM_HINT` (`:56-59`) naming both real remedies |
| 7 | Guard `--hardtrim3` with `HardtrimEnd::Three` | Step 7 | DONE | `main.rs:390-395` |
| 8 | Guard SE trim FASTQ before the loop | Step 8 | DONE | `main.rs:791-799`; A9 verified — candidate call argument-identical to the writer's |
| 9 | Guard SE trim uBAM | Step 9 | DONE | `main.rs:1863-1869`; A9 likewise |
| 10 | `cli.rs`: reject duplicate SE inputs, own message | Step 10 | **DEVIATED** (documented, D1) | `cli.rs:625-642`. Predicate `!self.paired && !self.clock && self.implicon.is_none()`, not the planned `!self.paired`. Verified necessary: the `--clock`/`--implicon` `validate_paired_input` calls sit *after* this check, and `--clock s.fastq.gz s.fastq.gz` still yields "Read 1 and Read 2 appear to be the same file". Acceptable |
| 11 | Tests: `io.rs` (V1), `cli.rs` (V4), new integration file; `current_dir` cases absolutised | Step 11 | DONE | 12 + 3 + 20 tests. The absolutised-fixture clause is met by a stronger route: all fixtures written in-test, so no `current_dir` case references a crate-root-relative path |
| 12 | CI: three validation steps | Step 12 | DONE | Steps 23/24/25 of `validation`; idiom matches the existing guards |
| 13 | CHANGELOG: bug fixes + changes | Step 13 | DONE | `#### Bug fixes` covers #383 and both defects with reachability; `#### Changes` covers A6 (flagged as a behaviour change, with the clock/implicon carve-out) and the `--output-dir` → `--output_dir` correction |

### §4 — Specified signatures

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 14 | `collision_key` / `preflight_output_collisions` | §4.1 | DONE | Character-identical. Input keys in `HashMap<String, &PathBuf>` rather than §4.1's `HashSet<String>` — **required**, since §3.2 demands the alias message name the aliased input, which a set of keys cannot recover. Doc-comment enumerates call sites plus "a new dispatch path needs a call here too" |
| 15 | `pub enum HardtrimEnd`; both namers `pub`, taking it | §4.2 | DONE | Exact; `as_str` left private (unspecified by the plan) |
| 16 | `planned_hardtrim_outputs(...) -> Vec<std::path::PathBuf>` | §4.3 | DONE | Exact, `PathBuf` spelled in full per §4.3/§8 |
| 17 | A6 in `Cli::validate`, before `ensure_output_dir`, own message | §4.4 | DONE | `validate()` `main.rs:155` vs `ensure_output_dir` `:314`; verified empirically that a rejected duplicate leaves no `-o` directory (A7). Sited at `cli.rs:625` rather than "beside `:534-558`" — cosmetic |

### §1 / §2.4 — Call-site count

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 18 | 11 call sites: 3 re-pointed, 4 converted, 4 new; no path missed | §1, §5 | DONE | Exactly **11**. Re-pointed `:548`/`:585`/`:633`; converted `:513`/`:734`/`:1825`/`:2430`; new `:360`/`:391`/`:799`/`:1869`. Cross-checked against all 16 `cli.input` iterations: `:365`←`:360`, `:396`←`:391`, `:514`←`:513`, `:591`←`:585`, `:634`←`:633`, `:744`←`:734`, `:800`←`:799`, `:1828`←`:1825`, `:1870`←`:1869`, `:2433`←`:2430`, the rest being candidate builders. `--clump_only --paired` FASTQ delegates to `run_specialty_paired` (guarded `:2430`); `--demux` runs inside the SE loop (guarded `:799`). No unguarded writing loop |

### §10 V1 — unit, `src/io.rs` (9 sub-cases)

All 9 present and individually located; all PASS. 1→`preflight_accepts_distinct_outputs`;
2→`preflight_rejects_identical_outputs`; 3→`preflight_rejects_case_only_variants`;
4→`preflight_rejects_dot_slash_alias`; 5→`preflight_rejects_absolute_versus_relative_alias`;
6→`preflight_rejects_output_that_aliases_an_input` (asserts the alias wording **and** the absence
of the duplicate wording); 7→`preflight_accepts_dotdot_alias_known_limitation` (pins A8, with a
comment saying to retire A8 deliberately if it ever fails); 8→`preflight_appends_hint_only_when_given`
(both directions); 9→`preflight_accepts_empty_and_single`. Bonus:
`preflight_rejects_input_alias_across_spellings`. Items 19–27 — all DONE.

### §10 V2 — integration rejections (11 sub-cases)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 28 | SE trim FASTQ, shared stem | V2.1 | DONE | `se_trim_rejects_shared_stem` |
| 29 | SE trim FASTQ, `./sample` + `sample` | V2.2 | DONE | `se_trim_rejects_dot_slash_alias`; runs **without** `-o`, asserts exact directory contents (D3) |
| 30 | `dirA`+`dirB` same basename **with `-o`** | V2.3 | DONE | `se_trim_rejects_same_basename_across_dirs_with_output_dir`. Paired against V3.1's no-`-o` acceptance, so the `output_dir`-vs-`None` slip in Step 8 is pinned from both sides — the plan's stated highest-value case, correctly built |
| 31 | SE trim uBAM, colliding stems | V2.4 | DONE | `se_trim_ubam_rejects_shared_stem` |
| 32 | `--hardtrim5 20` FASTQ, same basename two dirs, `current_dir(tempdir)` | V2.5 | **DEVIATED** (undocumented) | `hardtrim5_rejects_same_basename_across_dirs` sets `current_dir` but **also** passes `-o`, so it exercises the `-o` arm rather than the CWD arm. Asserts the hint. Behaviour covered by V7 repro (c) — see Reconciliation |
| 33 | `--hardtrim3 20` FASTQ, same shape | V2.6 | **DEVIATED** (undocumented) | `hardtrim3_rejects_same_basename_across_dirs`, same shape change |
| 34 | `--hardtrim5 20 --output-format ubam` | V2.7 | DONE | `hardtrim5_ubam_rejects_shared_stem` |
| 35 | `--hardtrim3 20 --output-format ubam` | V2.8 | DONE | `hardtrim3_ubam_rejects_shared_stem` |
| 36 | SE trim, output aliases an input | V2.9 | DONE | `se_trim_rejects_output_that_aliases_an_input`; asserts alias wording, the negative, **and** that the named input keeps its own 40 reads and none of the source's |
| 37 | `--paired`, output aliases an input | V2.10 | DONE | `paired_rejects_output_that_aliases_an_input`; proves the converted paired site gained defect-2 protection |
| 38 | `.gz` + `.bgz` sharing a stem | V2.11 | DONE | `se_trim_rejects_gz_and_bgz_sharing_a_stem` |

Bonus: `se_trim_rejects_absolute_versus_relative_alias` (the mixed absolute/relative list §2.2(e)
names) and `duplicate_input_gets_a_precise_message`.

### §10 V3 — integration acceptances (7 sub-cases)

Items 39–45, all DONE, all asserting content or filename:
V3.1→`se_trim_accepts_same_basename_across_dirs` (four assertions: each output holds 40 of its own
and 0 of the other's); V3.2→`se_trim_accepts_single_input`;
V3.3→`hardtrim5_accepts_distinct_stems_and_names_output`;
V3.4→`hardtrim3_accepts_distinct_stems_and_names_output` (first output-producing coverage
`--hardtrim3` has ever had); V3.5→`hardtrim3_ubam_accepts_distinct_stems_and_names_output`;
V3.6→`se_trim_ubam_accepts_distinct_stems`; V3.7→`paired_accepts_two_distinct_pairs`.
V3.3/V3.4 assert `.fq` rather than the plan's `.fq.gz` because the in-test fixtures are plain
FASTQ (`gzip=false`); the `5prime`/`3prime` discriminator — the point of the test under A9 and
§4.2 — is asserted.

### §10 V4–V7

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 46 | Duplicate SE input rejected with the dedicated message; distinct inputs validate; existing tests pass | V4 | DONE | `test_validate_single_end_duplicate_input_rejected` (positive **and** `!contains("APFS/NTFS")`), `test_validate_single_end_distinct_inputs_accepted`, bonus `test_validate_hardtrim_duplicate_input_rejected`. All pre-existing paired/clock/implicon validation tests green |
| 47 | A2 as a table that can fail, across the flag matrix, incl. demux stem and FastQC zip | V5 | DONE | **Closed in the 12:54–12:55 change.** `distinct_primary_outputs_imply_distinct_secondary_outputs` asserts the plan's stated direction (primaries differ ⇒ secondaries differ) over all 8 `--basename`/`--dont_gzip`/`-o` combinations for `report_name`, `json_report_name`, `clumping_report_name` and the demux stem via the newly extracted `demux::demux_base_name` (D6). FastQC's zip is **documented as not separately assertable** — the bundled crate derives it from the primary path handed to it, so there is no second key of ours to compare — with A2's `--fastqc_args -o` exception named so the row cannot be deleted as dead weight. The former inverse-premise test is retained as an explicit complement |
| 48 | CI: SE trim collision | V6.1 | DONE | Executed locally: rc=1, grep OK, `out/` empty |
| 49 | CI: `--hardtrim5 30` collision with `-o` | V6.2 | DONE | Also greps `one invocation per input`. Executed locally: rc=1, both greps OK, `out/` empty |
| 50 | CI: output-aliases-input | V6.3 | DONE | Greps the alias-specific wording plus an md5 byte-identity check on the named input and a file-count check — stronger than V6's stated idiom. Executed locally: rc=1, grep OK, md5 unchanged, count=2 |
| 51 | Gates, rebuild, all §2.2 reproductions, V3.1 positive control | V7 | DONE | Re-verified independently on a freshly rebuilt binary; all 8 reproductions reject, positive control clean. **Repro (c) is the no-`-o` hardtrim case** and it rejects with the hint, writing nothing to the CWD |

### §3.5 — Edge-case table (15 rows)

Items 52–66. All 15 rows trace. Abbreviated, with the two that need comment given in full:

1→V1.9 + V3.2 · 2→V2.1, V6.1, repro (a) · 3→V2.11, repro (b) · 4→V1.4, V2.2, repro (e-dot) ·
5→V1.5, bonus test, repro (e-abs) · 6→V1.7 (A8 residual, pinned) · 7→V1.6, V2.9, V2.10, V6.3,
repros (f-se)/(f-pe) · 9→V1.3 (unit only, correctly — two case variants cannot coexist on APFS) ·
10→V3.1 + the V7 positive control · 12→V2.3 · 13→`cli.rs:620-624` untouched, precedence unchanged ·
14→clap `required` (A5) · 15→§11 Open 3, explicitly left undone in §13.

- **Row 8 (same file listed twice) — DONE.** V4 + `duplicate_input_gets_a_precise_message`. Under
  `--paired`/`--clock`/`--implicon` the precise message comes from `validate_paired_input` instead
  (D1); a cross-pair duplicate (`--paired a1 a2 b1 a2`) falls through to the pre-flight's collision
  message. All three shapes verified rejecting, so the row's claim holds on every path even though
  the message differs by mode.
- **Row 11 (same basename, different dirs, no `-o` — hardtrim → reject) — DONE.** Traces to
  **V7 repro (c)**, which the plan mandates re-running and which I ran: rc=1, hint present, CWD
  left holding only `dirA dirB`. Not pinned by `cargo test` — every hardtrim test adds `-o` — so
  the pin is a manual gate rather than CI. See Reconciliation.

### Beyond-plan change

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 67 | `src/demux.rs`: extract `pub fn demux_base_name(&Path) -> String` | not in plan | **DEVIATED** (documented, D6) | Pure extraction so V5 can assert the demux stem against the real derivation rather than a copy of the formula; called from both `demultiplex` (`demux.rs:153`) and the V5 test, so the test cannot drift from the code. Verified behaviour-preserving: `--demux` on `test_files/demux_test.fastq.gz` still emits `demux_test_trimmed_sample{1..6}.fq.gz`, `_NoCode.fq.gz` and `_demultiplexing_summary.txt`. §13's D6 also records that the `trimmed_name` binding feeding the `Processed sequences from file >…<` line was restored after the first attempt dropped it — material, since `--demux` is in the Perl byte-identity matrix. Acceptable |

## Gaps (detail)

**No unresolved coverage gaps.** Every §5 step, every §10 sub-case, all 11 call sites, all 15 §3.5
rows and all four §4 signature groups are implemented as specified; every deviation is either
documented in §13 (D1, D6) or behaviourally covered by another plan item (V2.5/V2.6 → V7).

### Reconciliation — where this pass first disagreed with `COVERAGE.md`, and why it no longer does

This pass initially recorded the hardtrim `-o` shape change as a gap: V2.5 and V2.6 DEVIATED and
§3.5 row 11 PARTIAL, on the grounds that row 11 had no trace to a test, an assumption, or an
out-of-scope entry. **That was wrong, and I withdrew it.** §2.2 reproduction (c) *is* the no-`-o`
hardtrim collision — `trim_galore --hardtrim5 20 dirA/same.fastq.gz dirB/same.fastq.gz`, with the
plan's own note "(outputs land in CWD)" — and §10 V7 mandates re-running (a)–(f), explicitly
adding "Run (c) from a scratch directory, not the repo root." Row 11 therefore traces to a
plan-specified validation item, which I executed and which passes. `COVERAGE.md` reached the
correct classification.

Two independent passes converged on the same underlying fact from different directions, which is
worth more than either verdict alone:

- **Agreed fact.** All seven hardtrim tests and the new CI step pass `-o`, so `cargo test` never
  exercises `hardtrim_output_name`'s bare-filename/CWD arm (`specialty.rs:447`). The widest
  collision class in the change is pinned only by a manual gate.
- **Agreed classification.** Not a coverage gap: V2.5/V2.6 as written require `current_dir(tempdir)`,
  which the tests do set; `-o` was added for the reason the plan itself gives at V6.2 (the
  "nothing written" residue check needs a dedicated directory); and the candidate builder and the
  writer share one namer through `HardtrimEnd`, so the untested arm cannot silently drift out of
  sync (A9 + §4.2).
- **Agreed residual risk.** Low, and behavioural correctness is verified today.

Independently of the classification, `PROGRESS.md` already carries this as an open code-review
finding ("Hardtrim without `-o` — the widest collision class — is untested"), so it is captured
for the decision the reviewers own. Two cheap, non-blocking follow-ups:

1. Add one hardtrim rejection case with **no** `-o` — two same-basename inputs in different
   directories under `current_dir(tempdir)` — asserting the CWD holds only its two input
   directories. `assert_dir_holds_only` already exists in the file for exactly this shape; D3
   introduced it when the SE defect-1 cases had to drop `-o` for the same reason. This moves the
   pin from a manual gate into CI.
2. Add one line to §13 recording the V2.5/V2.6 shape change, so the only two undocumented
   deviations in the change become documented ones.

### Observations (not coverage gaps)

1. **§13's headline numbers are stale.** It records "520 tests pass (398 lib + 122 integration)"
   and "6 files changed, +503/−99"; the tree is now 521 (399 lib) and 7 files, +595/−111. The
   negative-control table says "40 steps in `validation`" where a YAML parse gives 41.
2. **Three comments went stale with Step 1**, all comment-only, none a plan item: `main.rs:271`
   ("the three output-collision pre-flights" — there are eleven); `main.rs:692-698` (describes the
   paired key as "the full path, case-folded", omitting the absolutisation that is the whole of
   defect 1's fix); `tests/integration_clump_only_ubam.rs:625` (cites the deleted
   `preflight_collision_bam`). `COVERAGE.md` independently flags a fourth, `demux.rs`'s
   re-parented `demultiplex` doc-comment.
3. **The hardtrim hint reads oddly under `-o`.** With `-o out` the message prints
   `out/same.30bp_5prime.fq.gz` while the hint says output "is written to the current working
   directory". Each half is defensible; the sentence contradicts the paths shown above it.
   `COVERAGE.md` and two of three code reviewers reached the same place — the hint should replace
   the generic advice rather than append to it.
4. **The duplicate-output message names the same path twice** ("X and X would be written to the
   same file") when the two planned paths are byte-identical. Pre-existing #216 wording, unchanged
   here.
5. **The `HashSet` → `HashMap` change at §4.1 is required, not incidental** (§3.2 needs the message
   to name the input). Worth a §13 line so nobody "restores" the plan's data structure.
6. **§3.5 has 15 rows, not the 14 the audit brief cited.** All 15 audited.
7. **D5 confirmed.** All three new CI step names truncate at `(issue` in the parsed YAML, matching
   the pre-existing `(issue #216)` step. Cosmetic.

## Test verification (Mode B)

`cargo test`, crate root, single run, no filters: **521 passed, 0 failed, 0 ignored** (399 lib +
122 integration). New: 12 unit in `src/io.rs`, 3 unit in `src/cli.rs`, 20 integration in
`tests/integration_output_collision.rs`.

Per-test status is identical to the table in `COVERAGE.md` — all 35 named new tests PASS, and I
verified the V1/V2/V3 sub-case counts independently as 9/9, 11/11 and 7/7 by name rather than by
file existence. The one row worth restating:

| Test | File | Status |
|---|---|---|
| hardtrim rejection with **no** `-o` (§3.5 row 11) | — | not in `cargo test`; covered by V7 repro (c), executed and passing |

Driver scripts for this pass:
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d0ee0171-fc72-476b-bdba-c46c90a5bacd/scratchpad/coverage/{ci_steps,repros,extra}.sh`

## Verdict

**Verdict:** COMPLETE

All 67 audited items resolve: 63 DONE, 0 PARTIAL, 0 MISSING, 4 DEVIATED — D1 and D6 documented in
§13 and verified against the code, and the V2.5/V2.6 test-shape change behaviourally covered by
V7 repro (c), which I executed. The implementation matches the plan, the 11 call sites reconcile
against `72624c7`, no dispatch path is unguarded, and the two gaps raised by
`COVERAGE_round1.md` — Step 1's doc-comment sub-clause and V5's A2 table — are closed, which I
confirmed independently rather than taking on report.

This pass concurs with `COVERAGE.md`, including on the one point where it initially differed
(§3.5 row 11 / V2.5 / V2.6), having found that row 11 traces to §10 V7. Nothing needs implementing
to satisfy the plan.

Two non-blocking follow-ups, both already visible to the reviewers via `PROGRESS.md`: move the
hardtrim no-`-o` pin from V7's manual gate into `cargo test`, and refresh §13's counts
(520→521, 6→7 files, 40→41 CI steps) plus the V2.5/V2.6 shape note.
