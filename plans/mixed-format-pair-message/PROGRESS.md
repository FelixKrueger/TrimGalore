# Progress: Accurate per-pair format rejection (#363)

**Last updated:** 2026-07-26

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2 + implementation notes (§13) + verification round (§14) |
| Plan Review | ✅ Complete | `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` — 1 Critical, 11 Important |
| Impl Plan | ✅ Complete | Folded into `PLAN.md` §4 — no separate `IMPL.md` for a change this size |
| Implementation | ✅ Complete | Branch `fix/mixed-format-pair-message`; 441 tests (was 421); fmt + clippy clean |
| Code Review | ✅ Complete | `CODE_review_reviewer-A.md`, `CODE_review_reviewer-B.md` — 4 defects found and fixed, 5 recommendations adopted |
| Coverage | ✅ Complete | `COVERAGE.md` — 45 items: 42 DONE, 0 PARTIAL, 1 MISSING (user-gated Step 12), 2 DEVIATED (documented) |

## Where this stands

Code, tests and CHANGELOG are done and independently verified. **441 tests pass, `fmt` and `clippy` clean, and byte-identity holds** for all six baseline outputs (PE ×2, SE, hardtrim5, clock ×2) — re-checked after the final edit.

Reviewer A upgraded V11 from deductive to empirical: it built `dev` @ `a4ffd47` in a throwaway worktree and md5-compared 15 outputs across PE, SE, `--hardtrim5`, `--clock`, `--demux` and `--clump_only`. All identical. That closes the one gap the coverage audit had declared.

## Closed out 2026-07-26

1. ✅ **§4 Step 12** — filed as [#364](https://github.com/FelixKrueger/TrimGalore/issues/364), labelled `bug`. Covers `--hardtrim5/3` silently accepting mixed pairs and `--clock`/`--implicon` rejecting them for the wrong reason, with the per-mode decision framed separately (hardtrim is genuinely harmless; clock/implicon is a real gap whose current rejection is incidental on record counts). This clears the only MISSING coverage item.
2. ✅ **Correcting comment on #363** — [issuecomment-5082658849](https://github.com/FelixKrueger/TrimGalore/issues/363#issuecomment-5082658849). Records that `main.rs:542` was the reference implementation rather than a third defect, that the late firing was not purely cosmetic (partial multi-pair output), and the `samtools collate` measurement.
3. ✅ **`docs/quickstart.md:57`** corrected — the false "Paired reads may come as two BAM files" claim. Verified it was the only occurrence across `docs/src` and `README.md`; the README was already correct.
4. ✅ **Committed** as `e20e4d4` on `fix/mixed-format-pair-message` — 7 files, +1466/−87. `PLAN.md` committed; `PROGRESS.md` and the five review reports stay session-local per the `plans/` convention.

5. ✅ **Pushed, PR opened, CI green, merged.** [#365](https://github.com/FelixKrueger/TrimGalore/pull/365) squash-merged to `dev` as `4276810`; remote branch deleted. All 8 checks passed, including `Validate vs Perl TrimGalore` (the authoritative byte-identity check) and `Validate uBAM input`. Verified `conclusion: success` on the run rather than reading check buckets alone.
6. ✅ **#363 closed** with `--reason completed` and a comment pointing at the merge commit. Manual close was required — auto-close fires only on the default branch (`master`), and this merged into `dev`.

## Feature complete

Nothing outstanding on this work. `dev` is now 32 commits ahead of `master`.

Carried forward: **#364** (specialty modes, filed and unstarted).

## Environment correction

**SSH to `github.com` works.** The previous session's handoff records it as timing out, with a recommendation to push via an explicit HTTPS URL. Re-tested 2026-07-26: `ssh -T git@github.com` authenticates on both port 22 and port 443, and `git push` over the configured SSH remote succeeded first try. The timeout was transient, so the HTTPS workaround is no longer needed — worth correcting in the next handoff so it does not propagate.

## Notes worth carrying forward

- **The most serious finding was a false user-facing claim in the CHANGELOG** — it stated the old two-BAM message "suggested `samtools collate -O r1.bam r2.bam`". It never did; the shipped message carried no `samtools` command at all, and that form existed only in this plan's v1 draft. Both reviewers caught it independently via `git show`. This is the second occurrence of this class of error (the #362 round caught three), so CHANGELOG claims about *prior* behaviour warrant a `git show` check before writing, not after.
- **The v1 plan's own Critical** (guard preempting the `--clump_only` FASTQ aux-tag diagnosis) was invisible in both the code and the plan text — it only appeared when someone ran that shape at N=2. Reproduction, not reading, is what found it.
- **Concurrent reviewer edits collided.** Two of B's `Edit` calls hit A's changes, and one fix was applied twice. B correctly warned neither report is a complete account of the tree. Future rounds: review read-only, or give each reviewer its own worktree.

## History

- 2026-07-26: Code Review + Coverage → ✅ Complete (4 defects fixed, 5 recommendations adopted, 9 plan-text corrections; verdict INCOMPLETE only on the user-gated Step 12)
- 2026-07-26: Implementation → ✅ Complete (441 tests, 17/17 validations, V11 byte-identity held)
- 2026-07-26: Plan → ✅ Complete (v2: 1 Critical + 11 Important incorporated; 13 v1 claims corrected)
- 2026-07-26: Plan Review → ✅ Complete (both reviewer reports written)
- 2026-07-25: Plan → ✅ Complete (PLAN.md v1 created)
