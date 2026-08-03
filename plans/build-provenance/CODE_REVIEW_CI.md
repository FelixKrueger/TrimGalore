# Code Review — CI / Dependabot / release.yml (single reviewer)

**Scope:**
- `/Users/fkrueger/Github/TrimGalore/.github/workflows/ci.yml` (modified) — full read
- `/Users/fkrueger/Github/TrimGalore/.github/workflows/release.yml` (modified) — new header block only
- `/Users/fkrueger/Github/TrimGalore/.github/dependabot.yml` (new) — full read

**Plan reference:** `plans/build-provenance/PLAN.md` §CI, §Dependabot, §release.yml.

**Verdict:** **PASS — every checklist item satisfied.** One Low-priority observation on the reproducibility job's cache `restore-keys` fallback and one Low-priority nit on the audit job's missing `schedule` trigger. No Critical, High, or Medium issues.

---

## 1. Checklist coverage

| # | Check | Location | Status |
|---|---|---|---|
| 1 | `lint` job runs `cargo fmt --all -- --check` + `cargo clippy --all-targets --release -- -D warnings`, own cache key, runs on push + PR | ci.yml:84-104 (triggers ci.yml:3-7) | PASS |
| 2 | `validation.needs` includes `lint` (downstream gates block on lint failure) | ci.yml:122 — `needs: [rust-tests, lint]` | PASS |
| 3 | `reproducibility` job: dedicated, separate cache key (not generic `cargo-<platform>-`); excludes `target/` from cache; uses `SOURCE_DATE_EPOCH=1700000000`; builds twice; diffs sha256; malformed-epoch hard-fail | ci.yml:47-82 | PASS (+Low note) |
| 4 | `audit` job: `rustsec/audit-check@v2` with `checks: write` + `issues: write`; runs on push + PR (minimum acceptable) | ci.yml:106-117 (triggers ci.yml:3-7) | PASS (+Low note) |
| 5 | `--version` verification: content regex `^[0-9a-f]{7,40} — `; negative check that `-V` lacks the provenance line | ci.yml:37-45 | PASS |
| 6 | `release.yml` header comment documents intentional non-gating with hotfix-tag caveat; NO actual lint gate added (documented deviation) | release.yml:15-22 | PASS |
| 7 | `dependabot.yml`: cargo + github-actions ecosystems, weekly Monday, limit 5, `FelixKrueger` assignee, no `labels:` field | dependabot.yml:1-33 | PASS |

---

## 2. Detailed findings

### Logic

**2.1 `--version` negative check is correctly anchored.**
`ci.yml:44-45` uses `! cmd | grep -qE '^[0-9a-f]{7,40} — '`. In bash pipelines, `!` inverts the exit status of the last pipeline element (grep). Under GitHub Actions' default `bash -e`, the negation correctly asserts grep found no match. The regex anchors to start-of-line and requires the literal em-dash separator, so false negatives are implausible. Correct.

**2.2 `--version` positive regex.**
`ci.yml:42` — `^[0-9a-f]{7,40} — [a-z0-9_]+/[a-z0-9_]+ — built [0-9T:Z-]+$`. The trailing `-` inside the final character class is a literal (not a range), so the class correctly matches ISO-8601 fragments like `2026-04-19T12:34:56Z`. Anchored at both ends — tight. Correct.

**2.3 Reproducibility job's `cargo clean` scope.**
`ci.yml:67-71` runs `cargo clean && cargo build --release` twice under identical `SOURCE_DATE_EPOCH`, then sha256s the binary. Because `target/` is excluded from this job's cache `path:` (ci.yml:59-62), the `cargo clean` calls don't pollute the saved cache. Cache save still writes `~/.cargo/registry` + `~/.cargo/git` only. Correct.

**2.4 Malformed-epoch assertion is sound.**
`ci.yml:72-82` captures the pipeline RC via `PIPESTATUS[0]` (the `cargo build`, not the `tee`), asserts non-zero, and greps the log for `SOURCE_DATE_EPOCH`. This correctly distinguishes a panic-with-message failure from an unrelated build failure. Correct.

**2.5 `validation` does not depend on `audit`.**
Plan §Step 4 is explicit that `audit` advisories are signals, not gates. `ci.yml:122` has `needs: [rust-tests, lint]` — audit correctly omitted. Matches plan.

**2.6 `release.yml` intentional non-gate.**
The header comment (release.yml:15-22) documents the deviation and its hotfix-tag caveat exactly as the plan prescribes (plan §Step 2, line 215: "Document this split in release.yml as a comment near the first job definition"). The comment also calls out the social-convention mitigation. No code-level lint dependency added — matches the "documented deviation, not fix" intent.

### Efficiency

**2.7 Lint cache key isolation.**
`cargo-lint-` prefix (ci.yml:99) is distinct from the generic `cargo-` used by `rust-tests` (ci.yml:27) and from `cargo-repro-` used by the reproducibility job (ci.yml:63). This prevents cross-job stomping on save, while `restore-keys: ${{ runner.os }}-cargo-` (ci.yml:100) still lets lint seed from any prior cargo cache. Sensible.

**2.8 Parallel job matrix.**
`rust-tests`, `lint`, `reproducibility`, and `audit` all run on `ubuntu-latest` in parallel. `validation` fans in on `[rust-tests, lint]`. No wall-clock regression vs. pre-change CI.

### Errors / risks

**2.9 [Low] Reproducibility `restore-keys` could seed target/ from rust-tests.**
`ci.yml:64` has `restore-keys: ${{ runner.os }}-cargo-` — this falls back to the generic `rust-tests` cache key prefix, which **does** include `target/` in its `path:` list (ci.yml:23-26). On save, the reproducibility job's `path:` excludes `target/`, so nothing leaks *out*. On restore, however, a matching generic-prefix cache *could* restore a stale `target/` into the workspace before `cargo clean` runs. In practice this is harmless because step 1 of the reproducibility run is `cargo clean`, which wipes any restored `target/`. But it's a minor waste of restore bandwidth, and it's also a subtle way for rust-tests cache state to affect this job's first build. **Recommendation:** either drop the generic `restore-keys` fallback entirely (force cold warm-up on this job) or narrow to `${{ runner.os }}-cargo-repro-`. Not blocking.

**2.10 [Low] `audit` job has no `schedule` trigger.**
The plan's verification checklist mentioned "Runs on push+PR+schedule or at minimum push+PR." As-implemented, `audit` inherits the workflow's push + PR triggers (ci.yml:3-7) — the minimum acceptable. This means a new CVE published against a dependency while no one is pushing won't surface until the next push. A weekly `schedule:` would close that gap. **Not required by plan; flagging for future consideration.**

**2.11 [Low] `audit` job + fork PRs.**
Plan §Step 4 already flagged that `checks: write` / `issues: write` can fail on fork PRs. As-implemented, nothing mitigates that yet (no `pull_request_target` variant, no annotations-only fallback). This is per-plan ("Note this to watch for"), not a gap. Flagging so it's visible if a first community PR arrives.

### Structure / style

**2.12 Dependabot config is clean.**
Both ecosystems use identical schedule/limit/assignee structure. Comments (dependabot.yml:1-13) clearly document the "no labels field" decision and the rationale for weekly Monday + limit-5. Matches plan §Step 3 verbatim.

**2.13 CI job ordering.**
`rust-tests` (with `--version` check) → `reproducibility` → `lint` → `audit` → `validation`. Reading order matches logical dependency order. Fine.

**2.14 Inline comments explain "why".**
Reproducibility job's `path:` comment (ci.yml:62) explains the `target/` exclusion. Lint job's cache key comment is in the plan but omitted in the file — minor, cache key name itself telegraphs intent. Release.yml's Lint-gating block is well-written and self-contained. Good.

---

## 3. Fixes applied

None. All findings are Low-priority observations. No edits made to CI / Dependabot / release.yml in this review.

---

## 4. Recommendations (Low priority)

1. **[Low]** Narrow `reproducibility` job's `restore-keys` from `${{ runner.os }}-cargo-` to `${{ runner.os }}-cargo-repro-` (ci.yml:64) to avoid pulling `target/`-heavy cache from `rust-tests` on cold starts. Cosmetic.
2. **[Low]** Consider adding `schedule: - cron: '0 6 * * 1'` (weekly Monday) to the `audit` job so new CVEs surface even during quiet weeks. Optional per plan.
3. **[Low]** If the `dependencies` label convention matters long-term, pre-create it via `gh label create dependencies --color 0366d6` before the first Dependabot PR lands — plan §Step 3 already flags this; re-surfacing for the merge checklist.

---

## 5. Summary

Implementation matches the plan on every checked item. The `lint` → `validation` gate, the isolated reproducibility job with `target/`-excluded cache, the content-regex `--version` assertion (plus negative `-V` check), and the `rustsec/audit-check@v2` job with the correct permissions are all present and correct. Dependabot config is minimal, labelled-as-planned (i.e. no `labels:` field), and matches the sample in plan §Step 3 line-for-line. `release.yml`'s intentional non-gating is documented in the header comment exactly as plan §Step 2 requires. No blocking issues.

**File:** `/Users/fkrueger/Github/TrimGalore/plans/build-provenance/CODE_REVIEW_CI.md`
