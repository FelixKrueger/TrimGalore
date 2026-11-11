# Plan Coverage Report

**Mode:** B (code vs. plan), post-implementation
**Plan(s):** `plans/08122026_fastqc-silent-noop/PLAN.md`
**Date:** 2026-08-13
**Verdict:** COMPLETE — every plan item is implemented; 4 non-blocking documentation defects are
listed under Gaps (2 in the plan text, 2 in `integration_preflight_tripwire.rs` comments)

Audited in a detached worktree (`3ffcef4` + the uncommitted working diff applied via
`git diff | git apply`), so the two concurrent code reviewers' tree was never touched.

## Summary

- Total items: 48
- DONE: 45
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 3 (all documented in the plan's Implementation notes; all supersets or idiom
  changes — none reduces coverage)

**Claimed numbers, independently re-measured in the worktree — all confirmed:**

| Claim | Measured | Verdict |
|---|---|---|
| `cargo test` 680 tests green | 680 passed / 0 failed / 0 ignored across 15 suites, exit 0 | ✅ |
| `cargo test --release` green | 680 passed / 0 failed, exit 0 | ✅ |
| `cargo nextest run` green | `Summary [21.146s] 680 tests run: 680 passed, 0 skipped`, exit 0 | ✅ |
| +17 over 663 | 17 added `#[test]` fns (cli.rs 12, paired_format_guard 3, clump_only 1, tripwire 1) | ✅ |
| clippy `--all-targets --release -D warnings` clean | exit 0 after a real 29.6 s compile (not a cache hit) | ✅ |
| `cargo fmt --all -- --check` clean | exit 0 | ✅ |
| docs build 31 pages | `[build] 31 page(s) built in 12.91s`, exit 0 | ✅ |

The docs build needed two attempts: a symlinked `node_modules` makes astro mangle its own module
paths (`…/cov2-421/docs/Users/fkrueger/…`) and fails in `vite`. That is a harness artefact, not a
docs defect — with a real (APFS-cloned) `node_modules` the build is clean. Recording it so the
first failure is not mistaken for a finding.

## Coverage ledger

### §3 Behavior — the four numbered behaviours

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| B1 | Add `Cli::fastqc_requested()` | §3.1 | DONE | `src/cli.rs:585-593` |
| B2 | Five `bail!` guards in `validate`, sited late | §3.2 | DONE | `src/cli.rs:1059-1100`, one `bail!` per arm, inside a single `if self.fastqc_requested()`. Sited **after** the `--adapter2` parse (`:1054`), i.e. the corrected §2.4/§5.1 location — see deviation D4 for the plan text that still names the older site |
| B3 | Pass `fastqc_requested()` at the five `main.rs` sites | §3.3 | DEVIATED | Ten sites, not five — documented deviation #1. The five `--clump_only` pass-throughs (`:636 :671 :715 :778 :819`) **and** the five gate expressions (`:1458 :1777 :2241 :2398 :2521`). Verified: `grep -rn "fastqc \|\| " src/` now matches exactly once (`cli.rs:592`) and no bare `cli.fastqc` remains in `main.rs` |
| B4 | Nothing else changes; no `fastqc::run` site gained or lost | §3.4 | DONE | `fastqc::run` line counts identical to base: `main.rs` 7→7, `clump_only.rs` 6→6. `git diff --name-only -- src/` = `cli.rs`, `main.rs` only; `clump_only.rs`, `specialty.rs`, `fastqc.rs` untouched as §7 requires |

### §3 Behavior — every row of the edge-case table

| # | Input | Expected | Status | Evidence |
|---|---|---|---|---|
| E1 | `--hardtrim5 20 --fastqc a.fq` | error naming `--hardtrim5` | DONE | `fastqc_plus_hardtrim5_rejected` (unit) + tripwire case `hardtrim5` (binary) |
| E2 | `--hardtrim5 20 --fastqc_args "--quiet" a.fq` | error via the disjunction | DONE | `fastqc_args_alone_also_rejected`, which also asserts `!cli.fastqc && cli.fastqc_requested()` |
| E3 | `--clock --paired --fastqc r1 r2` | error naming `--clock` | DONE | `fastqc_plus_clock_rejected` + tripwire case `clock` |
| E4 | `--implicon --paired --fastqc r1 r2` | error naming `--implicon` | DONE | `fastqc_plus_implicon_rejected` + tripwire case `implicon` |
| E5 | `--paired --fastqc one.bam` | error naming the shape | DONE | `fastqc_plus_paired_single_file_rejected` + `fastqc_refused_on_paired_interleaved_fastq_output` + tripwire case `paired_interleaved` |
| E6 | `--paired --output-format ubam --fastqc one.bam` | **accepted** | DONE | `fastqc_plus_paired_ubam_output_single_file_accepted` (unit), new `fastqc_accepted_on_paired_interleaved_ubam_output` (binary, asserts exactly one `_fastqc.zip`), existing `ubam_out_pe_fastqc_produces_exactly_one_report` |
| E7 | `--paired --output-format ubam --fastqc r1 r2` | accepted | DONE | `fastqc_plus_paired_ubam_output_two_files_accepted` (unit, as §5.3 specified) |
| E8 | `--output-format ubam --fastqc one.bam` (SE) | accepted | DONE | existing `ubam_out_se_fastqc_produces_report` (unchanged, green) |
| E9 | `--clump_only --fastqc a.fq` | accepted, unchanged | DONE | `fastqc_plus_clump_only_accepted` (unit); existing `fastqc_produces_report_on_ubam_out` covers the BAM-out driver |
| E10 | `--clump_only --fastqc_args "--quiet" a.fq` | accepted **and now reports** | DONE | `clump_only_fastqc_args_alone_produces_a_report` asserts exactly 2 artefacts; I confirmed it is red pre-fix (see V1b) |
| E11 | `--clump_only --paired --fastqc one.bam` | `main.rs:597`'s message still wins | DONE | `clump_only_paired_single_file_keeps_its_own_message` asserts `"requires two FASTQ input files"`, the wording at `main.rs:599` |
| E12 | `--hardtrim5 0 --fastqc a.fq` | range error wins | DONE | `hardtrim5_range_error_precedes_the_fastqc_refusal` |
| E13 | `--clock --paired --fastqc r1 r2 r3` | odd-count error wins | DONE | `odd_input_count_precedes_the_fastqc_refusal` |
| E14 | `--paired --fastqc single.fastq` | refusal whose wording is *true* for this user | DONE | Wording pinned by two tests (`… && err.contains("two files")`, plus `!err.contains("your input is")`). No test passes a plain FASTQ here, so I probed the built release binary directly: exit 1, the shape-naming message, output dir not created |
| E15 | `--hardtrim5 20 a.fq` (no FastQC flag) | unchanged | DONE | tripwire `hardtrim5_fastq` / `hardtrim5_ubam` still green with the marker removed |

### §4 Signature

| # | Item | Status | Notes |
|---|---|---|---|
| S1 | `pub fn fastqc_requested(&self) -> bool { self.fastqc \|\| self.fastqc_args.is_some() }` on `Cli`, with the "doc text, not a clap `requires`" rationale | DONE | `cli.rs:585-593`; body and doc comment match §4 verbatim in substance. No message helper added, as §4 required |

### §5 Implementation outline — the nine numbered steps

| # | Step | Status | Notes |
|---|---|---|---|
| I1 | `cli.rs`: accessor + five guards after the `--adapter2` parse; arm 5's message names the *shape* | DONE | Guard block at `:1059`, directly below `parse_adapter_specs_quiet(&self.adapter2)` at `:1054`. Arm-5 message asserts no format claim and offers three remedies, per §5.1 |
| I2 | `main.rs`: substitutions at `:636 :671 :715 :778 :819`; no `clump_only.rs` signature change | DEVIATED | Superset — see B3. `clump_only.rs` signatures unchanged as required |
| I3 | `cli.rs` unit tests: 5 refusals + 3 acceptance twins + 1 `--fastqc_args`-only + 2 precedence | DONE | 12 delivered = the 11 enumerated + `malformed_adapter2_precedes_the_fastqc_refusal` (documented deviation #5) |
| I4 | Tripwire: delete the machinery, collapse the match, rewrite the module doc, add five refusal cases | DEVIATED | All deletions complete; five arms are one looping `#[test]` (documented deviation #2) — see T-rows |
| I5 | `integration_paired_format_guard.rs`: arm-5 refusal + `--clump_only` precedence | DONE | Both present, plus an unrequested acceptance twin (`fastqc_accepted_on_paired_interleaved_ubam_output`) — a superset |
| I6 | `--clump_only --fastqc_args` regression test | DONE | `integration_clump_only.rs:405-426`; asserts `_fastqc.html` + `_fastqc.zip` |
| I7 | Help text: one clause on `--fastqc`; `--fastqc_args` unchanged | DONE | `cli.rs:330-332` adds the restriction; `:336` keeps "Implies --fastqc", which B3 makes true everywhere (plan-sanctioned, listed as deviation #3 but §5.7 predicted it) |
| I8 | Docs: four lines (`outputs.md:143`, `:145`, `quickstart.md:53`, `flags.md` cross-flag) | DONE | All four addressed. `outputs.md:143`'s blanket promise is qualified by a new paragraph immediately below it rather than by editing the sentence; `:145` now reads "once per trimmed primary output" and adds the `--retain_unpaired` exception. The two `/guide/outputs/#fastqc` links resolve — `## FastQC` exists at `outputs.md:141` |
| I9 | `CHANGELOG.md`: one sentence, `**Breaking:**` prefix | DONE | Present, `**Breaking:**` prefix, issue-linked in the file's format. It is one sentence spanning four wrapped lines and also states the `--fastqc_args` half — factually one sentence, but longer than the house one-clause norm |

### §2.5 — tripwire deletion set and its replacement

| # | Item | Status | Notes |
|---|---|---|---|
| T1 | Delete the `fastqc_capable` struct field | DONE | Gone from `Tripwire` and from `Tripwire::new` |
| T2 | Delete `.fastqc_incapable()` and its eight call sites | DONE | Builder gone; exactly 8 removals (`hardtrim5_fastq`, `hardtrim5_ubam`, `hardtrim3_fastq`, `hardtrim3_ubam`, `clock_multipair`, `implicon_pair`, `orphan_paired_fastq_from_interleaved_bam`, `orphan_paired_fastq_retain_unpaired`) |
| T3 | Delete `compare()`'s 6th positional parameter | DONE | `compare(root, created, planned, args, arm)` — 5 params |
| T4 | Delete the argument at `Tripwire::run` | DONE | Call collapsed to one line |
| T5 | Delete the literal `true` arguments in the `comparison_*`/helper callers | DONE | 8 sites updated, not the 9 the plan predicted: the file has 9 `compare()` call sites total, one of which is T4's. Deletion is provably complete — it compiles under `-D warnings` and all 54 cases pass |
| T6 | Collapse the match to `if active && allowed.is_empty()` | DONE | `:347-357`; the `(true, false)` arm and its "#421 may have been fixed" message are gone |
| T7 | Rewrite the module-doc paragraph (`:25-36`) | DONE | Now states the arms are *refused* by `Cli::validate` and that the surviving rule is unconditional |
| T8 | Add five *active* refusal cases (non-zero exit, message, nothing written) | DEVIATED | One `#[test]` looping over five labelled cases (documented deviation #2). All three assertions are made per arm, so the substance is complete; the case count is 53 → 54, not → 58 |

### §9 Validation — all twelve V-rows

| V | Check exists? | Verifies what the row says? | Status | Evidence |
|---|---|---|---|---|
| V1 | yes | yes | DONE | I re-measured this rather than trusting the recorded matrix: built the **unpatched** `3ffcef4` binary in a second worktree (dep cache shared) and ran `--hardtrim5 20 --fastqc BS-seq_10K_R1.fastq.gz` → exit 0, zero `*_fastqc.*`, zero case-insensitive "fastqc" in the log |
| V1b | yes | yes | DONE | Same base binary: `--clump_only --fastqc_args "--quiet"` → exit 0, zero `*_fastqc.*`, zero "fastqc" mentions; wrote only `_clumped.fq.gz` + `_clumping_report.txt`. The sixth arm's defect is independently reproduced |
| V2 | yes | yes | DONE | Five unit tests, each asserting the message names its arm; arm 5 asserts `"single input file"` (it names a shape by design, per §5.1) |
| V3 | yes | yes | DONE | `fastqc_plus_paired_ubam_output_single_file_accepted` — would fail if the `output_format` term were dropped |
| V4 | yes | yes | DONE | `fastqc_args_alone_also_rejected` (no `--fastqc` given) |
| V5 | yes | yes | DONE | `clump_only_fastqc_args_alone_produces_a_report` green; red on the base binary per V1b |
| V6 | yes | yes | DONE | Unit acceptance twins ×3 green; existing end-to-end guards green: `integration_ubam_out` 16/16, `integration_clump_only_ubam` 36/36 |
| V7 | yes | yes | DONE | Both precedence tests green; a third (`--adapter2`) exists beyond the row's scope |
| V8 | yes | yes | DONE | `fastqc_refused_on_paired_interleaved_fastq_output` asserts non-zero + the message + `!out.exists()` |
| V9 | yes | yes | DONE | `clump_only_paired_single_file_keeps_its_own_message` |
| V10 | yes | substance yes, count no | DEVIATED | `cargo test --test integration_preflight_tripwire` → `54 passed`. The row predicts 53 + 5 = 58 cases; deviation #2 makes it 54 with all five arms asserted inside one test |
| V11 | yes | yes | DONE | All six gates re-run and green — see the Summary table |
| V12 | yes | yes | DONE | `cargo test --test integration_ubam_out --test integration_clump_only_ubam` → 16 + 36 passed, exit 0 |

## Gaps (detail)

No plan requirement is missing or partially implemented. What follows is four documentation
defects found while auditing — three are stale text, one is a misattached comment in committed
test source.

### D1: a doc comment in the tripwire now documents the wrong test

**Expected:** five refusal cases added near `a_refused_run_writes_nothing`, per §5.4.
**Found:** the new test was inserted *between* `a_refused_run_writes_nothing` and its doc comment.
`tests/integration_preflight_tripwire.rs:1312-1316` now reads:

```rust
/// The tripwire tests the pre-flight's *coverage*; this pins its *purpose*. Two inputs
/// collapsing onto one output under `-o` must be refused before anything is written —
/// the contract #409 violated.
#[test]
fn refuses_fastqc_on_every_arm_that_cannot_run_it() {
```

**Gap:** that paragraph describes the `-o` output-collision contract (#409), which is
`a_refused_run_writes_nothing` at `:1357` — that test is now undocumented, and the #421 test
carries a comment about #409. Compiles and passes; the defect is purely in what a reader is told,
in the file whose stated purpose is anti-decay.

### D2: the tripwire's `fastqc_active` doc comment cites a code shape this change removed

**Expected:** not a plan item.
**Found:** `:219-223` still says "Production's gate is `cli.fastqc || cli.fastqc_args.is_some()`
(`main.rs:1458` and four siblings)".
**Gap:** after B3, `main.rs` spells that disjunction nowhere — those five sites call
`cli.fastqc_requested()`, and deviation #1 relies on the disjunction appearing exactly once, in
`cli.rs`. The comment's cross-reference is now false.

### D3: §3 item 2 and the Implementation notes still name the pre-correction guard site

**Expected (per the task's deviation 4):** the plan text and the code agree on siting.
**Found:** they agree in two places and disagree in two. §2.4 and §5 step 1 carry the correction
("after the `--adapter2` parse … v2 said the loop and that was one check short"), and the code
matches (`cli.rs:1059`, below `:1054`). But **§3 item 2** still reads "after the
`check_restartable_input` loop (`cli.rs:~1020`)", and the Implementation notes' *"What was built"*
paragraph still reads "five guards sited after the `check_restartable_input` loop".
**Gap:** both statements are technically true-but-imprecise (the block *is* after that loop, just
not adjacent to it) and they point a future reader at `~1020`, the location iteration #5 proved
insufficient. This is the "corrective copy contains the defect it corrects" pattern: the
correction landed in the two sections review touched, not in the behaviour spec or the summary.
Nothing in the code needs to change.

### D4: §9 has no V-row for the `--adapter2` precedence case

**Expected:** V7 covers precedence.
**Found:** V7 names only the range error and the odd-count error — the two cases §2.4's table
lists. Iteration #5 added a third requirement (a malformed `-a2` value must report first) and a
test for it (`malformed_adapter2_precedes_the_fastqc_refusal`, green), but no V-row.
**Gap:** plan-side only; the code and test exist. Worth noting because the mid-flight correction
to §2.4 was not propagated to §9 either — the same omission as D3.

## Test verification

All runs in the worktree `$TMPDIR/cov2-421` (= `3ffcef4` + working diff).

| Test | File | Status |
|---|---|---|
| `fastqc_plus_hardtrim5_rejected` | `src/cli.rs` | PASS |
| `fastqc_plus_hardtrim3_rejected` | `src/cli.rs` | PASS |
| `fastqc_plus_clock_rejected` | `src/cli.rs` | PASS |
| `fastqc_plus_implicon_rejected` | `src/cli.rs` | PASS |
| `fastqc_plus_paired_single_file_rejected` | `src/cli.rs` | PASS |
| `fastqc_args_alone_also_rejected` | `src/cli.rs` | PASS |
| `fastqc_plus_paired_ubam_output_single_file_accepted` | `src/cli.rs` | PASS |
| `fastqc_plus_paired_ubam_output_two_files_accepted` | `src/cli.rs` | PASS |
| `fastqc_plus_clump_only_accepted` | `src/cli.rs` | PASS |
| `hardtrim5_range_error_precedes_the_fastqc_refusal` | `src/cli.rs` | PASS |
| `malformed_adapter2_precedes_the_fastqc_refusal` | `src/cli.rs` | PASS |
| `odd_input_count_precedes_the_fastqc_refusal` | `src/cli.rs` | PASS |
| `fastqc_refused_on_paired_interleaved_fastq_output` | `tests/integration_paired_format_guard.rs` | PASS |
| `fastqc_accepted_on_paired_interleaved_ubam_output` | `tests/integration_paired_format_guard.rs` | PASS |
| `clump_only_paired_single_file_keeps_its_own_message` | `tests/integration_paired_format_guard.rs` | PASS |
| `clump_only_fastqc_args_alone_produces_a_report` | `tests/integration_clump_only.rs` | PASS |
| `refuses_fastqc_on_every_arm_that_cannot_run_it` (5 arms) | `tests/integration_preflight_tripwire.rs` | PASS |
| Whole tripwire suite (V10) | `tests/integration_preflight_tripwire.rs` | PASS — 54/54 |
| `integration_ubam_out` (V12) | `tests/integration_ubam_out.rs` | PASS — 16/16 |
| `integration_clump_only_ubam` (V12) | `tests/integration_clump_only_ubam.rs` | PASS — 36/36 |
| Full suite (V11) | — | PASS — 680/680 debug, 680/680 release, 680/680 nextest |

Pre-fix probes (V1, V1b) were run against a binary built from unpatched `3ffcef4` in a separate
worktree; both defects reproduced exactly as the plan's matrix records them.

## Verdict

**COMPLETE.** Every §3 behaviour, every edge-case row, the §4 signature, all nine §5 steps, all
twelve §9 V-rows and the whole §2.5 deletion set are implemented and verified. The three
deviations (ten substitution sites instead of five; the five tripwire refusal cases as one looping
test; `--fastqc_args`' help text untouched) are each documented in the plan's Implementation notes
and each is a superset or an idiom choice — none removes coverage. Every claimed number
re-measured true, including the 680-test count and the 31-page docs build.

Four documentation defects remain, none blocking and none requiring a behaviour change:

1. **D1 (committed source)** — `tests/integration_preflight_tripwire.rs:1312`: the #409
   output-collision doc comment now sits on `refuses_fastqc_on_every_arm_that_cannot_run_it`,
   leaving `a_refused_run_writes_nothing` undocumented. This is the one worth fixing before the
   commit lands.
2. **D2 (committed source)** — `:219-223`: `fastqc_active`'s comment cites
   `cli.fastqc || cli.fastqc_args.is_some()` at `main.rs:1458`, a shape this change deleted.
3. **D3 (plan)** — §3 item 2 and the Implementation notes' "What was built" paragraph still site
   the guard block "after the `check_restartable_input` loop (`cli.rs:~1020`)"; §2.4 and §5 step 1
   carry iteration #5's correction, and the code matches those two.
4. **D4 (plan)** — §9 has no V-row for the `--adapter2` precedence requirement that iteration #5
   added; the test exists and passes.
