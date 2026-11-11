# PROGRESS — #414 pre-flight empirical exact-set tripwire

**Feature:** empirical exact-set tripwire over the output-collision pre-flight, all arms
**Issue:** [#414](https://github.com/FelixKrueger/TrimGalore/issues/414)
**Status:** 🚧 Implemented and reviewed — **awaiting Felix's call on commit + PR**
**Base:** `dev` @ `3008b37`

## Pipeline

| Step | Artifact | State |
|---|---|---|
| Spike | `spikes/SPIKE_observed_write_sets.md` | ✅ 40 combinations probed, baseline green |
| Plan v1 | `PLAN.md` | ✅ written, reviewed by Felix, approved for agent review |
| Plan review (dual) | `PLAN_REVIEW_A.md` (461 ln), `PLAN_REVIEW_B.md` (543 ln) | ✅ complete |
| Plan v2 | `PLAN.md` | ✅ all accepted findings folded in; see §13 |
| Implementation plan | `IMPL.md` | ➖ skipped — v2 §5 is already task-ordered |
| Implementation | `src/io.rs` + `tests/integration_preflight_tripwire.rs` | ✅ 643 tests pass in all three CI runs |
| Code review (dual) | `CODE_REVIEW_A.md` (488 ln), `CODE_REVIEW_B.md` (411 ln) | ✅ no Criticals; both ran independent mutation tests |
| Coverage audit | `COVERAGE.md` (106 items) | ✅ INCOMPLETE — 5 items; 4 resolved in round 2, 1 deferred by Felix |
| Round-2 fixes | `PLAN.md` §13 | ✅ both Highs + all doc/plan corrections |
| Commit + PR | — | ⏸️ **awaiting Felix** |

## Review outcome

Both reviewers verified all 11 site line numbers, re-traced every leaf of `main()`
and independently confirmed no third orphan arm exists, and confirmed the arm→site
mapping. Reviewer A also checked T1 against all five historical instances.

**Two of v1's three assertions were weaker than advertised** — the same shape both
times: asserting more than the mechanism delivered.

| Finding | Origin | Resolution in v2 |
|---|---|---|
| T2 unfalsifiable (`created = after \ before` excludes an overwritten input) | both | T2 deleted, T3 promoted to "no pre-existing file's bytes changed" |
| Dump carries no site identity → §3.4 unimplementable, case↔arm unbound | both | `#[track_caller]` + `Location::caller()`; new T5 + site-distinctness test |
| FastQC non-empty guard false on 5 arms | both | scoped to `FASTQC_CAPABLE_ARMS`; predicate derived once from args |
| T1 liveness rested on an accident | both (A found the clap exit-0 mechanism) | new T4 non-emptiness assertions |
| `env::set_var` unit test illegal (edition 2024 + `ci.yml:52-55` invariant) | both | pure `dump_planned_lines` + end-to-end gate test |
| V2's edit range would not compile (918-924, not 919-925) | B only | corrected, plus five §10 ranges |
| `--cores > 1` swaps in a second writer implementation | B only | new §2.5 + cases 15, 17; v1's claim retracted |
| fastqc-rust does **not** create its `-o` directory | A only — **contradicted B** | settled by re-running without the `mkdir`; A correct |

**Two late items, both declined after checking the tree.** B appended ~1.7 KB after
its report was first consumed: a **retraction** of its §10 range finding (correct —
an interim revision had wrongly rewritten all five citations into block spans; now
reverted to the namer-call spans), and a new O-g arguing the FastQC predicate should
key on `cli.fastqc` rather than the args. O-g reads `cli.rs:333`'s "Implies --fastqc"
as an implemented constraint; it is a doc comment, and `--fastqc_args "--quiet"`
alone runs FastQC with `cli.fastqc` false (verified). Following it would have made the
tripwire go red on a legitimate case. §3.3 records the evidence so nobody re-simplifies
the predicate.

## Decisions taken (Felix)

| Question | Answer |
|---|---|
| Two arms with no pre-flight | Weak coverage (T3) + follow-up issue; #414 stays test-only |
| How the test sees the candidate list | Env-var dump inside `io::preflight_output_collisions`, no `cfg` gate |
| #398's deferred paired report namers | Its own issue after #414; design in `PLAN.md` §10 |
| Converse assertion (stale candidates) | Leave out — a second allowance is a second decay surface |
| `--fastqc` silent no-op | File as its own issue → **#421 filed** |

## Shape of the change

- `src/io.rs` — 2 `pub const`s, a pure `dump_planned_lines`, `#[track_caller]`, ~8
  lines of gated dump, 1 unit test. No candidate-assembly code touched, so candidate
  order and message content are unchanged.
- `tests/integration_preflight_tripwire.rs` — new: 40 arm cases + 5 self-tests.

## Follow-ups this plan spawns

1. **#421** — `--fastqc` silently ignored on five arms. **Filed.**
2. Dedicated paired report namers (`paired_report_names`, `interleaved_report_names`)
   — design ready in `PLAN.md` §10. Not yet filed.
3. Add the pre-flight to the two interleaved paired arms. Currently safe only because
   neither reads `cli.basename` (A11) — that dependency must go in the issue, since
   threading basename through would make it a live data-loss defect. Not yet filed.
4. Possible: the converse/stale-candidate assertion (§11 O2).
