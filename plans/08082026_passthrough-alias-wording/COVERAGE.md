# Plan Coverage Report

**Mode:** B (code vs design plan — no separate IMPL.md; ledger built from PLAN.md r2 Implementation outline [5 steps], Behavior [5 items], and Validation [6 rows])
**Plan(s):** plans/08082026_passthrough-alias-wording/PLAN.md (r2, including Implementation notes recording deliberate deviations)
**Implementation commit:** `1f43363` — fix(cli): say what the --passthrough alias check actually checks (#389); base: dev @ `7ea2741`
**Date:** 2026-08-08
**Verdict:** COMPLETE

## Summary

- Total items: 16
- DONE: 14
- PARTIAL: 0
- MISSING: 0
- DEVIATED (documented, non-behavioural): 2

Test evidence: full `cargo test` = **549 passed, 0 failed** (matches the Implementation notes' "549 green"); `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean (finished, zero diagnostics). Both 1.ix messages were additionally rendered end-to-end via the real binary and byte-match the plan's Behavior 2 texts.

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | 1.ix block rewrite: keep `len == 2` guard, `.iter().find(...)` on collision-key equality, branch on `path_identity_key`, two `bail!`s, new comment, trailing-space-before-`\` continuation | Outline 1 | DONE | `src/cli.rs:931-956`. Guard kept; `find` yields `matched`; identity branch first, case-variant branch second; both arms bail. Every continuation line carries a trailing space before `\` (verified with `cat -A`) — B-O3 satisfied; rendered output has correct spacing. |
| 2 | Companion fixes: cli.rs preamble + io.rs `norm_path` rustdoc | Outline 2 / Behavior 4 | DONE | Preamble (`src/cli.rs:885-887`): "(issue #216 protection)" removed, replaced with "1.ix is the input dual-consume guard (#389) — see its own comment for the key choice". Rustdoc (`src/io.rs:38-40`): stale "Used by `Cli::validate()`" bullet dropped; `norm_path` now described as "the folding component of `collision_key`"; no "aliasing" verb anywhere in shipped text. |
| 3 | Test updates: `_r1`/`_r2` pin identity-branch prefix + matched path + negative pin (helper allowed); new portable case-variant test; stale shared comment rewritten | Outline 3 | DONE | Helper `assert_passthrough_identity_rejection` (`src/cli.rs:1878-1893`) asserts prefix `"--passthrough must be a third file"`, `contains(matched)` (R1 in `_r1`, R2 in `_r2` — the discriminating assertion), and the negative `!contains("aliases an input")` with the plan's exact panic message. New `test_passthrough_rejects_case_variant_of_input` (`src/cli.rs:1912-1940`): writes only `<tmp>/r1.fq`, passes unwritten `<tmp>/R1.fq` as input 1; pins `"case-insensitively"` + negative pin. Rewritten comment states byte-equal ⊂ case-folded and points at `io::tests::test_norm_path_case_folds`; the forbidden "`collision_key` is a single `to_ascii_lowercase()`" substitution was not made (sentence removed). Note: the comment wraps to three physical lines (plan said two) and adds a pointer to the case-variant test — content matches the spec; cosmetic only. |
| 4 | CHANGELOG entry under `#### Changes` (not Bug fixes/Fixes), modelled on the `--output-dir` message-only entry | Outline 4 | DONE | `CHANGELOG.md:131-139`, inside `### Unreleased` (line 4) → `#### Changes` (line 129). States the message is now subcase-specific and "Which runs are accepted or rejected is unchanged." |
| 5 | Hygiene: fmt, clippy `-D warnings`, cargo test | Outline 5 | DONE | All three run in this audit: fmt clean, clippy clean, 549/549 tests pass. |
| 6 | Rejection surface unchanged: 1.ix fires exactly when `collision_key(pt)` equals either input's key | Behavior 1 | DONE | Same key, same comparisons (`find` over the same two inputs), same position after 1.viii; both arms bail so the inner branch cannot change accept/reject. Acceptance path re-confirmed green (`test_passthrough_paired_pair_accepted`). |
| 7 | Two-tier message discriminated by `path_identity_key`, exact texts | Behavior 2 | DONE | Both messages rendered via `cargo run` and byte-match the plan's texts, including "{pt} is input {matched}" ordering in the identity arm and "matches input {matched} case-insensitively (for APFS/NTFS safety): {pt}" plus rename advice in the case-variant arm. |
| 8 | New two-line 1.ix comment keeping the 1.ii clause | Behavior 3 | DONE | `src/cli.rs:931-932` — verbatim match to the plan's specified two lines, including "len == 2 per 1.ii." |
| 9 | Companion sites: no identity claim, no #216 mislabel | Behavior 4 | DONE | Same evidence as item 2; grep confirms no "aliasing R1 or R2" anywhere in shipped text. |
| 10 | No new edge cases; pre-existing ordering artifact stands (not reordered) | Behavior 5 | DONE | 1.viii `check_restartable_input` still at `src/cli.rs:930` (before 1.ix); general input existence still checked later (`src/cli.rs:996`). No reordering in the diff. |
| 11 | V1: `cargo test` all pass; three passthrough-rejection tests pin prefix, matched path, negative assertion | Validation 1 | DONE | 549 passed / 0 failed. `_r1`/`_r2` pin all three via the helper; the case-variant test pins adverb + negative (matched-path pin applies to the identity-branch tests per plan step 3). |
| 12 | V2: case-only branch reachable and truthful on both filesystem families | Validation 2 | DONE | `test_passthrough_rejects_case_variant_of_input` passes (lexical check; verified here on APFS, ext4-identical by construction since neither key touches the filesystem). Runtime render of the branch confirmed. |
| 13 | V3: stale-string grep over `src/ tests/ docs/ .github/ README.md CHANGELOG.md` — plan said zero hits | Validation 3 | DEVIATED (documented) | Actual: exactly 3 hits — `src/cli.rs:1890` and `:1937` (the two negative assertions `!err.contains("aliases an input")`) and `CHANGELOG.md:133` (describing the old text). Matches the Implementation notes' deviation verbatim; no live message, comment, or doc carries the old framing. |
| 14 | V4: count greps — plan said `"must be a third file"` ×3 (message + 2 test pins) and `"case-insensitively (for APFS/NTFS"` ×2 (message + 1 test pin); pins must not straddle `\` breaks | Validation 4 | DEVIATED (documented) | Actual: ×2 (message + shared helper — pin centralized) and ×1 (message only; the test pins the bare adverb `"case-insensitively"`). Matches the Implementation notes' revised counts exactly. Both pin prefixes verified to sit on a single source line before any `\` continuation. |
| 15 | V5: prose sanity — common case gets the third-file instruction, rare case gets rename advice, neither claims identity unconditionally | Validation 5 | DONE | Both messages rendered end-to-end (real binary, real files). Identity arm: names the matched input, normative "must be a third file (e.g. the index read)". Case-variant arm: states the case-insensitive comparison and its APFS/NTFS rationale, gives the rename remediation, and conditions the identity claim on "On a case-insensitive filesystem". |
| 16 | V6: fmt + clippy `-D warnings` clean | Validation 6 | DONE | Both clean at `1f43363`. |

## Gaps (detail)

None. The two DEVIATED items (13, 14) are recorded in the plan's own Implementation notes with rationale ("the assertions ARE the enforcement"; "centralizing the pin means one test-side occurrence") and are non-behavioural; audit re-ran both greps and reproduced the documented results exactly.

## Test verification

| Test name | File | Status |
|-----------|------|--------|
| `cli::tests::test_passthrough_rejects_pointing_at_r1` | src/cli.rs | PASS |
| `cli::tests::test_passthrough_rejects_pointing_at_r2` | src/cli.rs | PASS |
| `cli::tests::test_passthrough_rejects_case_variant_of_input` (new) | src/cli.rs | PASS |
| `cli::tests::test_passthrough_paired_pair_accepted` (accept surface) | src/cli.rs | PASS |
| `io::tests::test_norm_path_case_folds` (fold pin, referenced by new comment) | src/io.rs | PASS |
| `io::tests::identity_key_is_case_sensitive_but_collision_key_is_not` (key-contrast pin) | src/io.rs | PASS |
| Full suite (`cargo test`, all targets + integration) | — | PASS — 549 passed, 0 failed, 0 ignored |
| `cargo fmt --all -- --check` | — | PASS |
| `cargo clippy --all-targets --release -- -D warnings` | — | PASS |

Runtime renders (beyond unit pins):

- Identity branch: `Error: --passthrough must be a third file (e.g. the index read), not one of the R1/R2 inputs: test_files/BS-seq_10K_R1.fastq.gz is input test_files/BS-seq_10K_R1.fastq.gz`
- Case-variant branch: `Error: --passthrough matches input <tmp>/R1.fq case-insensitively (for APFS/NTFS safety): <tmp>/r1.fq. On a case-insensitive filesystem these are the same file and the stream would be consumed twice; if they are genuinely two files, rename one so the paths differ by more than letter case.`

## Verdict

**COMPLETE.** All 5 Implementation outline steps, all 5 Behavior items, and all 6 Validation rows are satisfied at `1f43363`. The two grep-count deviations (Validation 3 and 4) exactly match the deviations documented in the plan's Implementation notes and change no behaviour, message text, or enforcement — nothing remains to be addressed.

---

## Post-audit note (2026-08-08, caller)

The audit's COMPLETE verdict covered commit `1f43363`. A follow-up review batch
(`478e1d0`) then applied the dual code review's Low items: identity-first
matched-input selection (rejection surface unchanged — identity-equal implies
collision-equal), positional "is input {matched}" test pins plus a new
./-spelled discrimination test, two comment trims, and two rustdoc touches.
Suite green (14/14 binaries), fmt + clippy clean after the batch.
