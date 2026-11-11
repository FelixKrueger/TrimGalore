# PROGRESS — #428: a mid-stream refusal must leave no output file

**Plan:** `PLAN.md` v2 (revised after dual review)
**Issue:** [#428](https://github.com/FelixKrueger/TrimGalore/issues/428)
**Branch base:** `dev` @ `285afba`
**Test baseline:** 684

| Step | State | Notes |
|---|---|---|
| Scope measured | **done** | 9 arms probed at `285afba`; residue on every output encoding, §2.4 |
| Plan written | **done** | `PLAN.md` v1 |
| Design decision (O1) | **awaiting Felix** | Option B (temp + rename) recommended over A (unlink) |
| Felix's manual review | **done** | 2026-08-14: "launch reviewers" |
| Dual plan review | **done** | `PLAN_REVIEW_A.md` + `PLAN_REVIEW_B_gen2.md`. 3 Critical (all A), 2 Critical (B2), one direct contradiction on O5 |
| Plan revised to v2 | **done** | Every Critical and Important folded in; §12 records what changed |
| Implementation | **done** | Branch `fix/428-no-partial-output-on-refusal`. 684 → **707** tests, fmt + clippy clean. 4 deviations in `PLAN.md` §12 |
| Mutation test | **done** | Pre-#428 behaviour planted → 15 tests red; `src/io.rs` restored md5-verified |
| Dual code review + coverage audit | **done** | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`, `COVERAGE.md`. No Critical; coverage INCOMPLETE 4/42, nothing MISSING |
| Review fixes applied | **done** | All High/Medium taken; §5 step 6's commit-naming promise withdrawn as unimplementable; two new tests proved load-bearing by re-running the reviewers' mutations |
| PR / merge / close | **done** | PR #433, CI 9/9 green by name, squash-merged as `d2ea293`; #428 closed manually |
| Follow-up filed | **done** | #434 — teardown errors discarded across all three sinks, plus the two parallel-path cases |

## Open items carried in the plan

- **O1** — both reviewers answered: Option B, conditional on §5's corrections. Felix's
  sign-off still wanted, since §3 now states A's flatter risk profile as a real argument.
- **O2** — **answered**: dot-prefixed sibling `.{file_name}.partial`, decided on glob
  visibility (`*.fq*` matches `x_trimmed.fq.gz.partial`; a dot-prefixed name matches nothing).
- **O3** — do *not* relax `main.rs:395-405`'s refusal here; edit only its stated reason.
- **O4** — **dissolved**: `test_files/ubam_ws_qname_paired.bam` already exists and is already
  driven by an existing test. No decision, no new fixture.
- **O5 (new)** — the sink enum that would make teardown errors propagable. **Deferred**; the
  one review finding that clears the filing bar, to be filed as a single issue after the PR.

## The thing most likely to go wrong

**The commit point is not a free choice.** Four consumers read the output by its *final* path
inside the same run — FastQC, `--demux`, `write_clump_only_report`, and
`std::fs::metadata(&output_path)?` at `clump_only.rs:827`/`:1071`. That last is
`?`-propagating on the **success** path, so committing any later than today's teardown point
fails a **good** run with `ENOENT`. v1 asserted no such consumers existed.
