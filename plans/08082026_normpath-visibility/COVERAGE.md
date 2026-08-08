# Plan Coverage Report

**Mode:** B (code vs plan — no separate `IMPL.md`; ledger built from the plan's Implementation outline, Behavior, Assumptions A1–A5, Validation table + caveat, and the one code-affecting item in Questions/Resolved)
**Plan(s):** `plans/08082026_normpath-visibility/PLAN.md` (r2)
**Date:** 2026-08-08
**Verdict:** COMPLETE

Branch audited: `fix/399-normpath-private` @ `315557d` (base `dev` @ `ac08a97`).
Diff: `src/io.rs` only — 10 insertions, 12 deletions, single commit
`315557d refactor(io): make norm_path private and re-home the fold trade-off (#399)`.

## Summary

- Total items: 19
- DONE: 19 (one of them, outline step 5, was marked optional and was done anyway)
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 0

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `pub fn norm_path` → `fn norm_path` | Outline 1 | DONE | `src/io.rs:39` is `fn norm_path(p: &Path) -> String`. Line moved from `:45` to `:39` because the doc block above it shrank (step 3). Body byte-identical. |
| 2 | Move the trade-off paragraph from `norm_path` onto `collision_key` | Outline 2 | DONE | Removed from `norm_path`; present verbatim (3 lines) at `src/io.rs:77-79` inside `collision_key`'s block. |
| 3 | Reduce `norm_path`'s doc to one line; no "(private)" framing | Outline 3 | DONE | `src/io.rs:38` — ``/// Case-folded (ASCII lowercase) view of a path; the folding half of `collision_key`.`` — the plan's own suggested wording, one line, no visibility narration. |
| 4 | Add the one-line note that input identity uses the case-**preserving** `path_identity_key`; plain backticks, no `[…]` intra-doc links | Outline 4 | DONE | `src/io.rs:74-75`, appended to `collision_key`'s first paragraph: "Input identity asks a different question and uses the case-preserving `path_identity_key`." Plain backticks; no bracketed link anywhere in the file. Placement resolves the plan's ambiguous "while in the block" — it landed on `collision_key`, which is the only reading consistent with step 3's one-line cap on `norm_path`. |
| 5 | Optional drive-by: fix the test comment's imprecision | Outline 5 (**explicitly optional**) | DONE | `src/io.rs:765-766` rewritten to "The fold `collision_key` is built on (callers reach it through that key, never directly)." The false "pre-flight + `--passthrough` both use `norm_path`" claim is gone. The plan's second sub-point (that `--passthrough` also uses the case-sensitive key) is resolved by deleting the `--passthrough` reference rather than by expanding the comment; that contrast now lives on `collision_key`'s doc (item 4) and in the pre-existing test `identity_key_is_case_sensitive_but_collision_key_is_not` (`src/io.rs:1176`). No claim in the comment is now imprecise. |
| 6 | Run the four verification commands | Outline 6 | DONE | All four run, all green — see Test verification. |
| 7 | No CHANGELOG entry | Outline 7 (deliberate decision) | DONE | `git diff --name-only ac08a97..315557d` → `src/io.rs` only. The plan's corrected rationale checks out: `#### Infrastructure (contributor-facing)` does exist (`CHANGELOG.md:481` and `:1148`), and `grep -E 'collision_key\|norm_path\|path_identity_key' CHANGELOG.md` over 1586 lines returns **zero** hits, so the changelog has never named these helpers. |
| 8 | No behavioural change: same fold, same caller, same tests, same rendered output | Behavior | DONE | Diff is one visibility keyword plus doc/comment lines. `norm_path`'s body (`p.to_string_lossy().to_ascii_lowercase()`) and `collision_key`'s body (`norm_path(&lexical_normalise(p))`) are untouched; all five assertions in `test_norm_path_case_folds` unchanged; 561 tests green. |
| 9 | **A1** — zero callers outside `io.rs`, aliases included | Assumptions | DONE | `git grep norm_path` over the tracked tree: `io.rs:39` (def), `io.rs:81` (sole production caller, inside `collision_key`), `io.rs:764-777` (the test), and `cli.rs:1903` — a `//` comment naming the *test*, unaffected. Repeating the grep over **all** files including untracked ones (excluding `plans/`, `target/`, `.git/`) adds nothing. `naming::norm_path`: zero hits, while both alias sites still exist (`src/specialty.rs:12`, `src/clump_only.rs:52`). |
| 10 | **A2** — `mod tests` reaches the private item | Assumptions | DONE | Verified by execution, not by citing the rule: `io::tests::test_norm_path_case_folds ... ok` (1 passed). |
| 11 | **A3** — published library surface is uncurated and treated as an implementation detail; the removal is technically semver-breaking and accepted | Assumptions | DONE | `src/lib.rs` has 17 `pub mod` and **no** `pub use`; `Cargo.toml` has no `publish = false` and no `[lib]` section — so the plan's accurate framing (rather than r1's denial) holds. No code obligation in this plan. The durable follow-up was in fact filed separately: issue **#402** "lib.rs publishes 17 uncurated pub modules — declare the library surface unstable". |
| 12 | **A4** — no behavioural coupling to the fold's *name* | Assumptions | DONE | `--passthrough` validation (`src/cli.rs:934-952`) uses `path_identity_key` **and** `collision_key`; the pre-flight (`src/io.rs:98`, `:103`) uses `collision_key`; `src/main.rs:790` and `:2541` reference `collision_key` in comments. No caller names `norm_path`. |
| 13 | **A5** — no doc-tests reference the symbol | Assumptions | DONE | The only doc fence in `src/` is the ```` ```text ```` block at `src/specialty.rs:235-240`, and the full suite reports `Doc-tests trim_galore … running 0 tests` — the crate has no doc-tests at all, so none can reference the symbol. |
| 14 | **V1** — no missed caller, via `cargo clippy --all-targets --release -- -D warnings` | Validation | DONE | Exit 0, zero `warning:`/`error:` lines. Harness independently verified capable of seeing an external caller (below). |
| 15 | **V2** — fold still pinned, via `cargo test test_norm_path_case_folds` | Validation | DONE | Exit 0; `1 passed`. |
| 16 | **V3** — no collateral: full `cargo test` + `cargo fmt --all -- --check` | Validation | DONE | 561 passed / 0 failed across 14 result lines; `cargo fmt --all -- --check` exit 0 (so step 3's shortened doc line did not disturb rustfmt). |
| 17 | **V4** — trade-off note landed where readers arrive | Validation | DONE | `collision_key`'s block (`src/io.rs:73-79`) now states the *what* (output collision, case-folded, #216/#383), the contrast with input identity, and the false-positive cost; `norm_path` is one line. |
| 18 | **Caveat** — untracked `examples/fastqc_only.rs` means a local `--all-targets` run compiles a file CI lacks | Validation caveat | DONE | Confirmed the example is untracked and *is* compiled locally (it appears as `kind:["example"] name:"fastqc_only"` in the clippy target enumeration). It does not reference `norm_path` (the all-files grep in item 9 covers it). Because the extra target is purely additive and the lib/bin/`tests/` set is identical, the local green run is a **superset** of CI's, so the tracked-tree conclusion holds. |
| 19 | Commit message notes `pub(crate)` was considered and rejected | Questions/Resolved (r2) | DONE | `315557d` body ¶3: "pub(crate) was considered and rejected: specialty.rs and clump_only.rs both alias crate::io as naming, so pub(crate) would keep norm_path reachable as naming::norm_path from exactly the modules whose readers the change is meant to protect." — B's argument, as the plan settled it. |

## Gaps (detail)

None. No item is PARTIAL, MISSING, or DEVIATED.

Two items are worth reading before signing off, not as gaps but as judgement calls the auditor made explicit:

### Item 4 — where the `path_identity_key` note landed

**Expected:** outline step 4 says "While in the block, add the one-line note…". "The block" is not named.
**Found:** the note is on `collision_key` (`src/io.rs:74-75`), not on `norm_path`.
**Assessment:** DONE, not DEVIATED. Putting it on `norm_path` would have made that doc two lines and contradicted step 3's explicit one-line cap; and the note's substance (two keys, two questions) is a property of the key-comparison, which is `collision_key`. The reading chosen is the only self-consistent one.

### Item 5 — how the test comment's second imprecision was fixed

**Expected:** step 5 flags two faults in the old comment — (a) it credits `norm_path` with what the pre-flight and `--passthrough` validation do, true only indirectly; (b) it omits that the `--passthrough` check also uses the case-sensitive key for a different question.
**Found:** the new comment fixes (a) directly and dissolves (b) by dropping the `--passthrough` reference entirely, rather than by adding the case-sensitive-key sentence there.
**Assessment:** DONE. Step 5's stated purpose is removing imprecision, and no imprecise claim survives; the omitted contrast is now stated on `collision_key`'s doc (item 4), where the plan wanted it. This was also the optional step.

## Test verification

| Check | Command | Result |
|---|---|---|
| Clippy, all targets, deny warnings | `cargo clippy --all-targets --release -- -D warnings` | **PASS** — exit 0, no `warning:`/`error:` lines |
| Fold pinned (private-fn reachability) | `cargo test test_norm_path_case_folds` | **PASS** — `io::tests::test_norm_path_case_folds ... ok`, 1 passed, 409 filtered out |
| Full suite | `cargo test` | **PASS** — 561 passed, 0 failed, 0 ignored, across lib (410), bin (0), and 11 integration binaries; `Doc-tests trim_galore` 0 tests |
| Formatting | `cargo fmt --all -- --check` | **PASS** — exit 0 |

**Harness verified capable of failing (not assumed).** The plan's V1 rests on `--all-targets` compiling `tests/`, which is the only place an out-of-module caller of a `pub` item could hide. Rather than take that on faith, the clippy invocation was re-run with `--message-format=json` and its target list enumerated: it covers `lib trim_galore`, `bin trim_galore`, the `custom-build` script, the untracked `example fastqc_only`, and **all 11 test targets** (`integration_adapter2`, `integration_clump_only`, `integration_clump_only_ubam`, `integration_gzip_non_gz_extension`, `integration_no_args_help`, `integration_non_restartable_input`, `integration_output_collision`, `integration_paired_format_guard`, `integration_passthrough`, `integration_ubam`, `integration_ubam_out`). So an `E0603` from a planted external caller would be seen by exactly this command. The commit message additionally records that the implementer ran that positive control on this branch (planted caller → `E0603` under `--all-targets`, invisible to plain `cargo check`).

## Verdict

**COMPLETE — 19/19 items DONE, 0 unresolved.**

Every step of the Implementation outline is present in `315557d`, including the optional step 5; step 7's "no CHANGELOG entry" is a satisfied decision rather than an omission (the diff touches `src/io.rs` alone). The Behavior contract holds — the only executable change in the diff is the removal of the `pub` keyword. All five assumptions were re-verified against the code rather than carried over from the plan's prose, and all four validation rows were executed green in this working tree, with the validation caveat about the untracked example resolved in the change's favour (local target set is a superset of CI's; the example does not touch the symbol).

Nothing remains for the implementer.
