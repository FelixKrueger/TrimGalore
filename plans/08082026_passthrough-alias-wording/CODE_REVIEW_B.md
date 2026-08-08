# Code Review B — #389 `--passthrough` alias-check wording (commit `1f43363`)

**Reviewer:** B (independent, fresh context)
**Target:** branch `fix/389-passthrough-wording` @ `1f43363` vs base `dev` @ `7ea2741`
**Plan:** `plans/08082026_passthrough-alias-wording/PLAN.md` (r2, incl. Implementation notes)
**Files:** `src/cli.rs`, `src/io.rs`, `CHANGELOG.md` (confirmed via `git diff --stat` — no stray changes)

## Summary

The implementation matches the approved plan, including its two documented (non-behavioural) deviations. The accept/reject surface is provably unchanged: both arms of the new two-tier branch `bail!`, and the branch is reached only inside the pre-existing collision-key match under the retained `len == 2` guard — the reject set is identical to the old `||` check. Each message is truthful for its subcase, the `./`-spelling case lands in the identity branch (correctly), the new portable case-variant test is sound on both ext4 and APFS, and all four backslash-continuation joins were verified byte-level (`cat -A`) to carry the required space. Targeted tests (12 passthrough-filter tests, including the three at issue) pass; `cargo fmt --check` is clean; validation greps 3 and 4 reproduce exactly the counts the plan's Implementation notes document.

**Verdict: APPROVE.** Four Low-priority recommendations, none blocking; three are trivial fixes.

## Area 1 — Logic

### Verified sound

- **Accept/reject surface preserved.** Old: `pt_norm == key(input[0]) || pt_norm == key(input[1])` → bail. New: `find(|p| collision_key(p) == pt_key)` → `Some` → bail via either arm. Under the retained `self.input.len() == 2` guard these are extensionally identical; the inner `path_identity_key` comparison selects only the sentence (`src/cli.rs:933-957`).
- **Case-variant branch is exactly "differs only by ASCII case".** `path_identity_key` = `lexical_normalise(p).to_string_lossy()`; `collision_key` = the same string `to_ascii_lowercase()`d (`src/io.rs:75-83`). Same lossy conversion of the same normalised `PathBuf`, so collision-equal + identity-different ⟺ ASCII-case-only difference (non-UTF-8 bytes become U+FFFD in both, consistently). The case-variant message's claims are therefore truthful for every path that reaches that arm.
- **Identity branch truthful everywhere.** Identity-key equality means the two spellings lexically name one path; "{pt} is input {matched}" holds on every filesystem. `./`- and `..`-spellings fold in `lexical_normalise` (CurDir dropped, ParentDir popped, `src/io.rs:52-70`), so `./r1.fq` vs `r1.fq` takes the identity branch, where the message is truthful for it — as the plan requires.
- **Portable case-variant test is sound on both filesystem families** (`src/cli.rs:1914-1941`). 1.ix is lexical end to end; the only file opened before the bail is the passthrough itself (1.viii `check_restartable_input(pt, …)` at `cli.rs:930`); input existence is checked later (`cli.rs:995-997`); the earlier duplicate-pair check (`cli.rs:537`) is identity-keyed, so temp `R1.fq` vs the R2 fixture passes it. On APFS the unwritten `R1.fq` would alias the written `r1.fq` on disk, but the test never opens it. The pid-suffixed temp dir is unique per test process and the prefix is unique in the suite; cleanup runs before the asserts, so an assertion failure does not leak the dir (only an unexpected `Ok` from `validate()` would, via the `unwrap_err` panic — acceptable in a test).

### Findings

- **L1 (Low): subcase mis-prioritisation when both inputs collision-match the passthrough.** Inputs `["R1.fq", "r1.fq"]` are legal (duplicate detection is identity-keyed per #383; on ext4 they are genuinely two files). With `--passthrough r1.fq` — byte-equal to input 2 — `find` returns the *first* collision match, input 1 (`R1.fq`), so the **case-variant** message fires with rename advice, masking that the passthrough *is* input 2 verbatim. Every sentence emitted is still factually true, and the run is still rejected, but the remediation is wrong for the situation (renaming `R1.fq` just moves the user to the identity error on the next run). Recommend preferring an identity match among collision matches, e.g. find on `path_identity_key` equality first, falling back to the collision match. Exotic input set; message-selection only; does not block.
- **L2 (Low): the helper's matched-path assertion is not discriminating in the two identity tests.** In `test_passthrough_rejects_pointing_at_r1`/`_r2`, `matched` is textually equal to `pt`, so `err.contains(matched)` is satisfied by the `{pt}` placeholder alone — it cannot detect a swapped or omitted `{matched}` argument, which is the failure A§4.4's "discriminating assertion" was meant to catch. A `./`-spelled-pt test (`--passthrough ./test_files/BS-seq_10K_R1.fastq.gz R1 R2`) would (a) genuinely discriminate the two placeholders (`./…` ≠ `test_files/…` as displayed) and (b) pin the currently-untested plan/comment claim that `./`-spellings take the identity branch. Portable (fixture exists; unit tests run from crate root). Trivial fix — one added test.

## Area 2 — Efficiency

Nil, as the plan predicted. Error-path only: `find` short-circuits (marginally fewer `collision_key` evaluations than the old unconditional pair on an input-1 hit); the two `path_identity_key` evaluations occur only on the bail path, after which the process exits.

## Area 3 — Errors

- **Backslash-continuation whitespace verified byte-level** (`sed | cat -A` over `src/cli.rs:940-956`): all four joins carry a space before the `\` (`not·\`, `APFS/NTFS·\`, `same·\`, `genuinely·\`) → renders "not one", "APFS/NTFS safety", "same file", "genuinely two". No silently joined words. Both pinned prefixes (`"--passthrough must be a third file"`, `"case-insensitively"`) sit on single source lines, satisfying Validation 4's no-straddle requirement.
- **No stale pins or framing in shipped text.** Validation 3's grep returns exactly the three plan-documented deliberate hits (two negative assertions at `cli.rs:1890/1937` + the CHANGELOG description of the old text). Validation 4 counts reproduce: `"must be a third file"` ×2 (message + shared helper), `"case-insensitively (for APFS/NTFS"` ×1 (message).
- **No CI-grep interference.** `ci.yml:636-639` pins the *comma* form `"case-insensitive, for APFS/NTFS safety"` on an output-collision scenario; the new passthrough string uses the paren form `"case-insensitively (for APFS/NTFS safety)"` and is unreachable from that CI step. The plan's deliberate near-but-not-identical phrasing (B-O4) is honoured.
- Targeted runs: 12 passthrough-filtered unit tests green (including all three at issue); `cargo fmt --all -- --check` clean. Full suite (549) + clippy `-D warnings` verified green before commit per the caller; not re-run in full here.

## Area 4 — Structure

- 1.ix comment is exactly the plan's two-liner, keeping the 1.ii clause for the retained redundant guard. Preamble (`cli.rs:885-887`) corrected — the false "#216 protection" label is gone, replaced with input-dual-consume framing and a pointer.
- `norm_path` rustdoc (`io.rs:38-44`) now accurate: verified `norm_path` has no callers outside `collision_key` and its own tests, and both live collision checks (1.ix, `preflight_output_collisions`) hash `collision_key` — the "every collision check" claim holds.
- Helper `assert_passthrough_identity_rejection` is well-named and centralizes the pins; the resulting count deviation from Validation 4 is documented in the plan's Implementation notes and is an improvement, not a gap.
- CHANGELOG entry sits under Unreleased → `#### Changes` as directed, models the message-only entry style, states "Which runs are accepted or rejected is unchanged", and keeps the plain register.
- **S1 (Low, trivial fix): two new 3-line test comments** exceed the house ≤2-line comment rule and the plan's own "two lines" spec for the shared test comment: `cli.rs:1897-1899` (byte-equal ⊂ case-folded) and `cli.rs:1915-1917` (case-variant test rationale). Each compresses to two lines without losing the load-bearing facts, e.g. `// Byte-equal ⊂ case-folded (fold pinned by io::tests::test_norm_path_case_folds); / // the case-variant branch has its own lexical test below.` and `// Lexical check: only the passthrough must exist (1.viii); input existence is / // validated later, so the unwritten case-variant behaves identically on ext4/APFS.`
- **S2 (Low, trivial fix): `path_identity_key` rustdoc** (`io.rs:72-74`) cites only #383; check 1.ix (#389) is now its second caller. The sentence stays literally true ("input identity in `Cli::validate`"), so this is a completeness nit — append `#389`.

## Recommendations (priority order)

| # | Priority | What | Where | Trivial fix? |
|---|----------|------|-------|--------------|
| L1 | Low | Prefer identity match among collision matches so a passthrough byte-equal to input 2 isn't reported as a case-variant of input 1 | `src/cli.rs:935-938` | Near-trivial (identity-first `find`, collision fallback) |
| L2 | Low | Add a `./`-spelled-pt identity test to discriminate `{pt}` vs `{matched}` and pin the `./`-spelling→identity-branch claim | `src/cli.rs` tests | Trivial fix |
| S1 | Low | Compress the two 3-line test comments to 2 lines (house rule; plan spec) | `src/cli.rs:1897-1899`, `1915-1917` | Trivial fix |
| S2 | Low | Append #389 to `path_identity_key`'s issue list | `src/io.rs:74` | Trivial fix |

No Critical, High, or Medium findings. No fixes were applied (review-only worktree, per instructions).
