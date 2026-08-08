# Plan: `--passthrough` alias check — keep the fold, say what it checks (#389)

**Issue:** [#389](https://github.com/FelixKrueger/TrimGalore/issues/389) — filed from the D13 focused review.
**Decision (Felix, 2026-08-08):** keep `collision_key`; fix the words. Rationale: both keys are lexical and `path_identity_key` preserves case — on APFS/NTFS a case-variant spelling IS the input file, the case-preserving key can't see it, and `--passthrough` would silently dual-consume one stream at exit 0. Reviewer B additionally verified 1.ix is the **sole** guard for this (the output pre-flight checks output-vs-input and output-vs-output only, never input-vs-input).
**Base:** `dev` @ `7ea2741`.

## Revision history

- **r2 (2026-08-08):** incorporates dual plan-review (PLAN_REVIEW_A.md, PLAN_REVIEW_B.md; both REVISE). Fixes the false `norm_path()` premise (A§1.1/B-C3); widens scope to `io.rs:38-47` and `cli.rs:885-889` (A§1.6/B-C2); replaces the single message with a two-tier message that names the matched input and gives each subcase its own remediation (A§1.4-1.5/5.3, B-I2/I3); hardens the test pins with a negative assertion (B-C1); adds a portable case-variant test (supersedes the "cross-case is untestable" caveat); scopes and extends the grep validations (B-I4/I5, A§4.3); shrinks the comment to two lines keeping the 1.ii clause (A§1.2-1.3/B-I1); restates the behaviour claim precisely (A#7/B§1.5). The single-plan-review suggestion is withdrawn — B's closing note ("small diff, not small review": three criticals hid in a two-string plan) is accepted.
- r1 (2026-08-08): initial plan.

## Goal

Make check 1.ix's words match its semantics. `Cli::validate` §1.ix (`src/cli.rs:933-949`) compares `--passthrough` against R1/R2 with the case-folded `collision_key`, but its comment calls that catching "case-only INPUT aliases" and its message asserts "{} **aliases an input file**" — an identity claim that is false on a case-sensitive filesystem. **The accept/reject surface does not change**: same key, same comparisons, same position in the validate chain. What changes: the message text (now subcase-specific and naming the matched input), three comments, and the two rustdoc/preamble sites that carry the same mis-framing.

## Context

- `src/cli.rs:930-949` — check 1.ix. Runs after 1.viii (`check_restartable_input`), so the passthrough file always exists when 1.ix fires; R1/R2 existence is checked only later (`cli.rs:987-989`), so 1.ix is lexical end to end.
- `io::collision_key` (`io.rs:84-86`) = `norm_path(&lexical_normalise(p))`. **`norm_path` (`io.rs:48-50`) is a live component function, not a stale name** — r1 got this wrong; both reviewers caught it. The fold is already pinned platform-independently by `io::tests::test_norm_path_case_folds` (`io.rs:768`) and `io::tests::identity_key_is_case_sensitive_but_collision_key_is_not` (`io.rs:1181`).
- Message pins: only `src/cli.rs:1880/1887` (both `contains("cannot point at R1 or R2")`, byte-equal fixtures). No pins in `tests/integration_passthrough.rs` (its only stderr pin is 1.i), no CI grep (`ci.yml:638`'s grep targets the *output*-collision message, unreachable from a passthrough run), no docs/README copies. Verified independently by both reviewers.
- **In-scope companion sites carrying the same mis-framing** (review scope finding):
  - `src/io.rs:38-47` — `norm_path`'s rustdoc claims direct use by `Cli::validate()` (false post-#385: validate calls `collision_key`) and uses the "aliasing R1 or R2" identity verb.
  - `src/cli.rs:885-889` — the passthrough-envelope preamble labels 1.ix "(issue #216 protection)"; #216 is the *output*-collision issue, 1.ix guards *input* dual-consume.

## Behavior

1. Rejection surface unchanged: 1.ix still fires exactly when `collision_key(pt)` equals the collision key of either input.
2. **Two-tier message.** Inside the existing match, discriminate the subcase with `path_identity_key` (case-preserving; `./`/`..` spellings of one file compare equal, so they take the identity branch, which is truthful for them on every filesystem):
   - **Same file (identity keys equal)** — the common slip (R1 passed where the I1 index read was meant): name the matched input, state the requirement, give the correct remediation:
     `--passthrough must be a third file (e.g. the index read), not one of the R1/R2 inputs: {pt} is input {matched}`
   - **Case-variant only (collision keys equal, identity keys differ)**: state the case-insensitive comparison honestly and give the rename remediation:
     `--passthrough matches input {matched} case-insensitively (for APFS/NTFS safety): {pt}. On a case-insensitive filesystem these are the same file and the stream would be consumed twice; if they are genuinely two files, rename one so the paths differ by more than letter case.`
   Both arms `bail!`; the branch chooses only the sentence. This satisfies: no unconditional identity claim; the matched input is named (PR #219 both-paths precedent); the common case keeps a normative, actionable instruction; the rare case gets truthful advice. Qualifier phrasing deliberately close to the house idiom but not byte-identical to the CI-pinned output-collision string (B-O4 noted; A-#10 accommodated).
3. **New 1.ix comment**, two lines, fact-stated, keeping the invariant clause that explains the (deliberately retained) redundant `len == 2` guard:
   ```rust
   // 1.ix — case-folded on purpose, not an identity check (#389): a case-variant
   // passthrough IS R1/R2 on APFS/NTFS and would be consumed twice. len == 2 per 1.ii.
   ```
   Derivation goes to the commit message and #389.
4. **Companion fixes:** `io.rs:38-47` — drop the stale first "Used by" bullet; describe `norm_path` as the case-folding component of `collision_key` (no "aliasing" verb). `cli.rs:885-889` — replace "(issue #216 protection)" with input-dual-consume framing or a bare "see 1.ix" pointer.
5. **Edge cases:** none new; the pre-existing ordering artifact stands (a missing R1 plus a case-variant passthrough yields the 1.ix message rather than "input not found", because input existence is checked later) — acknowledged, not reordered (B-O5).

## Implementation outline

1. `src/cli.rs` 1.ix block: keep the `len == 2` guard; find the matched input (`.iter().find(...)` on collision-key equality); branch on `path_identity_key` equality; the two `bail!`s from Behavior 2. Replace the comment per Behavior 3. Watch the `\` string-continuation whitespace (a trailing space belongs *before* each backslash — B-O3).
2. `src/cli.rs:885-889` preamble and `src/io.rs:38-47` rustdoc per Behavior 4.
3. **Tests** (`src/cli.rs` unit tests):
   - `test_passthrough_rejects_pointing_at_r1` / `_r2`: pin the identity-branch prefix `"--passthrough must be a third file"`, pin the matched path (R1 in the first, R2 in the second — the discriminating assertion A§4.4 asked for), and add the negative pin `assert!(!err.contains("aliases an input"), "1.ix must not assert identity (#389); got: {err}")` (B-C1). Factor the shared asserts into a small helper if it reads better.
   - **New, portable:** `test_passthrough_rejects_case_variant_of_input` — write only the passthrough file (`<tmp>/r1.fq`) and pass the unwritten case-variant `<tmp>/R1.fq` as input 1. The check is lexical, so this hits the case-only branch identically on ext4 and APFS (1.viii needs only the *passthrough* file to exist; input existence is checked later). Pin `contains("case-insensitively")` and the negative pin. This supersedes r1's "cross-case is untestable" framing — that limitation belonged to the fixture-reuse approach, not the check.
   - Rewrite the stale shared test comment (do **not** substitute `collision_key` into the "single `to_ascii_lowercase()`" sentence — that would be false): two lines stating byte-equal ⊂ case-folded and pointing at `io::tests::test_norm_path_case_folds` for the fold.
4. `CHANGELOG.md` — Unreleased → `#### Changes` (named explicitly; not "Bug fixes"/"Fixes"): 1.ix message reworded and made subcase-specific (#389); accept/reject behaviour unchanged. Model: the `--output-dir` message-only entry at `CHANGELOG.md:145-146`.
5. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`.

## Efficiency

Error-path only: one extra `path_identity_key` evaluation and one branch, then the process exits. Nothing measurable.

## Integration

No keys, no ordering, no callers change. The diff is no longer string-only (a `find` + one `if` on the error path) — accepted by both reviewers as inside "wording" since both arms bail identically. Only the two named unit tests pin the message; both are updated in the same diff. `plans/` history files retain the old strings by design (see Validation).

## Assumptions

- **A1 (narrowed per B-I5):** no *test* pins of the old text outside `cli.rs:1880/1887` — verified by both reviewers across `src/ tests/ docs/ .github/ README.md CHANGELOG.md`.
- **A2:** D13 taxonomy unchanged — no third key, no platform probe (rejected in #385: syscall/symlink semantics; a probe would also make the message platform-dependent and untestable).
- **A3 (new):** nothing machine-parses TrimGalore stderr in-repo (verified) — what licenses a free-form rewrite; wrapper scripts grepping the old string are served by the CHANGELOG entry.
- **A4 (new):** `pt.display()` naming paths as spelled (not normalised) is the codebase-wide diagnostic convention (`CHANGELOG.md:109-110`) — kept.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Behaviour unchanged + pins updated | `cargo test` | all pass; the three passthrough-rejection tests pin prefix, matched path, and the negative "aliases an input" assertion |
| 2 | Case-only branch reachable and truthful | new portable unit test (impl step 3) | rejected via the case-variant message on both filesystem families |
| 3 | No stale pins/framing left in shipped text | `grep -rn "cannot point at R1\|aliases an input\|aliasing R1 or R2" src/ tests/ docs/ .github/ README.md CHANGELOG.md` | zero hits (`plans/` deliberately excluded — historical artifacts keep the old strings; do not edit them) |
| 4 | New text landed everywhere intended | grep `"must be a third file"` (expect: message + 2 test pins) and `"case-insensitively (for APFS/NTFS"` (expect: message + 1 test pin); chosen pin prefixes must not straddle a `\` line break | exact expected counts |
| 5 | Prose sanity (human) | read both rendered messages against their subcases | common case gets the third-file instruction; rare case gets rename advice; neither claims identity unconditionally |
| 6 | Hygiene | fmt + clippy `-D warnings` | clean |

## Questions or ambiguities

- **[Open] Message copy** is contract-shaped now (prefixes are pinned) but reviewers may still tighten words that aren't pinned.
- **[Resolved r2]** Review weight: dual review ran; its findings are folded in here.

## Follow-ups (out of scope, to file separately)

- Consider making `norm_path` private or inlining it — `pub` with no callers outside `collision_key` and its own tests; the public name is what enabled r1's false-premise error (A-#12).

## Implementation notes (2026-08-08, branch `fix/389-passthrough-wording` @ `1f43363`)

Implemented as planned: two-tier message with `path_identity_key` branch, 2-line 1.ix comment keeping the 1.ii clause, preamble + `norm_path` rustdoc companions, shared test-assert helper, new portable `test_passthrough_rejects_case_variant_of_input` (verified green on APFS; lexical so ext4-identical), CHANGELOG under `#### Changes`. Full suite 549 green; fmt + clippy `-D warnings` clean.

**Deviations (documented, none behavioural):**
- Validation 3's grep returns three deliberate hits: the two negative assertions (`!err.contains("aliases an input")`) and the CHANGELOG entry describing the old text. The grep's intent — no *live* message or comment carries the framing — is met; the assertions ARE the enforcement.
- Validation 4 counts: `"must be a third file"` ×2 (message + the shared assert helper — centralizing the pin means one test-side occurrence, not two); `"case-insensitively (for APFS/NTFS"` ×1 (message; the new test pins the bare adverb).

**Iteration log:**
- #1: initial implementation (`1f43363`); no failed iterations.
- #2 (`478e1d0`): review batch per dual code review — identity-first matched-input selection (both reviewers' L1; surface unchanged), positional `is input {matched}` pins + `./`-spelled discrimination test (B-L2), two comment trims (A#2/B-S1), rustdoc verb fix (A#3) and #389 on `path_identity_key` (B-S2). Not taken: A's two-tier-vs-single message alternatives beyond what shipped (already satisfied), norm_path privatization (filed as follow-up).

- **Logic:** the two-tier branch cannot change the rejection surface (both arms bail; the branch is reached only inside the existing match). The identity-branch discriminator (`path_identity_key`) deliberately groups `./`-spellings with byte-equal — truthful for them everywhere.
- **Edge cases:** re-checked the new portable test's reachability — 1.viii requires only the passthrough file; input existence is validated later (`cli.rs:987-989`), so the unwritten case-variant input is fine on both platforms.
- **Traps:** every reviewer-cited anchor re-verified against the tree before adoption (`io.rs:768/1181`, `CHANGELOG.md:129/145`); grep validations scoped so "zero hits" is achievable (B-I4).
- **Remaining risks:** none beyond copy quality; the deliverable (honest text) is now held by a negative assertion, two positive pins, and a count-checked grep.
