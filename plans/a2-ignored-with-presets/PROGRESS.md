# Progress: `-a2` must always win over preset R2 adapters (#369)

**Last updated:** 2026-07-29

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2 + implementation notes (§13) |
| Plan Review | ✅ Complete | `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` (read-only) |
| Impl Plan | ✅ Complete | Folded into `PLAN.md` §4 |
| Implementation | ✅ Complete | Branch `fix/a2-honoured-with-presets`; 441 → **456** tests; fmt + clippy clean; all 23 validations pass |
| Code Review | ✅ Complete | `CODE_review_reviewer-{A,B}.md` — 4 defects found and fixed, read-only |
| Coverage | ✅ Complete | `COVERAGE.md` — 53 items: 48 DONE, 1 PARTIAL, 1 MISSING (both low-severity) |

## Where this stands

`-a2` now wins on all seven resolution paths. **456 tests, fmt and clippy clean**, and no accepted path moved — all 24 pre-change baseline files (FASTQ *and* reports) are byte-identical.

Two checks carried the weight:

- **V16, the differential oracle** — post-fix `--illumina -a2 SEQ` and auto-detect `-a2 SEQ` are byte-identical to the *pre-fix* `-a AGATCGGAAGAGC -a2 SEQ` output, on both reads. Needs no Perl and pins both reads at once.
- **V15** — the new CI gate was run against a `dev` build in a throwaway worktree and **fails** there, so it is not vacuous. That mattered: a gate that passes before the fix is exactly how this defect survived three releases.

## Scope decisions taken with the user

1. `-a` + preset stays **permissive** (Perl rejects it; it is the #369 workaround).
2. A displaced preset R2 default **announces itself** with a NOTE.
3. `--consider_already_trimmed` suppression **wins over** `-a2`, with a warning.

## Deviations from the plan (all in `PLAN.md` §13; a fourth was added after review)

1. **Suppression detected from the R1 adapter, not `DetectionResult.suppressed`.** The plan said `suppressed` was "already in scope" — true inside `resolve_adapter`, but not in the caller where the override ended up, and propagating it would have reinstated the signature change the plan had just removed. Keys on "R1 adapter is a single empty sequence" instead, which is the actual semantics and what `trimmer.rs` already skips on.
2. **No `md5` dev-dependency** — the test compares decompressed bytes directly.
3. **One extra docs edit** — `--bgiseq` needed the same hedge as `--small_rna`.
4. **No unit test of the override helper**, contrary to §4 Step 5. `main.rs` has no test module and the helper is private, so the plan's own escape hatch applies — but the omission was undisclosed until the coverage audit caught it.

## Shipped

- **Committed** `04c9463`, pushed, **[PR #370](https://github.com/FelixKrueger/TrimGalore/pull/370)** open against `dev`, MERGEABLE, CI running.
- **Commented on #369** ([issuecomment-5119791043](https://github.com/FelixKrueger/TrimGalore/issues/369#issuecomment-5119791043)) — the reporter can drop their workaround, and is warned that Read 1 output changes too, so the fixed run should be compared against their 0.6.11 results rather than their 2.3.0 ones.
- Verification round found **4 defects**, all fixed; 2 of my own claims corrected. Full record in `PLAN.md` §14. Tests 441 → **456**.

## Closed out

- **#370 squash-merged** to `dev` as `3ec0ff9`; remote branch deleted; local `dev` fast-forwarded. All 8 checks passed with `conclusion: success` on the run — verified at run level, not from check buckets.
- **`Validate vs Perl TrimGalore` passed** (4m30s). That was the **first real execution of the new `-a2` gate**, against conda-installed Perl 0.6.11 on a Linux runner. Neither the conda cutadapt build nor Linux-vs-macOS tie-breaking in adapter alignment diverged from the local result, which was the open risk.
- **#369 closed** `completed`, with the merge commit and the parity result.

## Feature complete

Nothing outstanding on this work. `dev` is now 33 commits ahead of `master`. Carried forward: **#364** (specialty modes, filed and unstarted).

## Notes worth carrying forward

- **`zsh` does not word-split unquoted variables** the way bash does. A test loop passing `m="--hardtrim5 30"` as `$m` handed clap one bogus argument, producing two false "NO WARNING" readings while the direct invocation warned correctly. Second time this session a broken harness produced a confident wrong conclusion (the first was `mktemp` under the sandbox). Verify the harness before believing a negative result.
- **The Read 1 guarantee here is structural, not tested.** No branch's Read 1 computation is in the diff — the only edit inside `resolve_adapter` is the `-a` branch's Read 2 slot — so Read 1 cannot move. That is why the caller-side override was the right call over the plan's original seven-branch restructure: the plan's own Medium risk disappeared rather than being mitigated. (Both code reviewers flagged an earlier phrasing of this, "zero branches were touched", as false — the `-a` branch *was* edited. Precision matters here because the claim is what carries the Read 1 argument.)

## History

- 2026-07-29: Implementation → ✅ Complete (454 tests, 23/23 validations, V13 byte-identity held, V15 proven non-vacuous)
- 2026-07-29: Plan → ✅ Complete (v2: 1 Critical + ~15 Important; 13 v1 claims corrected; design changed)
- 2026-07-29: Plan Review → ✅ Complete (both reports, read-only)
- 2026-07-29: Plan → ✅ Complete (PLAN.md v1 created)
