# Plan Coverage Report

**Mode:** B (code vs. plan, post-implementation — §5 of `PLAN.md` carried the implementation outline;
there is no separate `IMPL.md`)
**Plan:** `plans/11112026_release-notes-policy/PLAN.md` (revision 3)
**Implementation:** `8f984ac` (code) + `97288b6` (plan/progress), on `dev`
**Date:** 2026-11-18
**Verdict:** **INCOMPLETE — 2 items unresolved**

Both unresolved items are small and in the two hand-written documents, not in the mechanism. The
mechanism — Steps 2, 2a, 3 and 5a — is complete, matches the plan's snippets, and its edge-case
behaviour was executed rather than reasoned about.

## Summary

- Total items: 26
- DONE: 22
- PARTIAL: 2 (Step 1's worked example; Step 4's contributor acknowledgements)
- MISSING: 0
- DEVIATED: 1 (Step 2's condition — declared in §14, accurate)
- Correctly absent: 1 (Step 5b, dropped by §13 Open-5)

## Coverage ledger

### §5 — implementation steps

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `.github/release-notes/README.md` states the §3 policy, the four-part shape, **and one worked example** | Step 1 | **PARTIAL** | Policy, shape and mechanism all present and correct. No worked example, and no pointer to `v2.4.0.md`, which §4 designates as the worked example. See Gap 1. |
| 2 | README states the *mechanism*, not only the prose shape — job-0 failure, and that the error names the path | Step 1, §8.5 | DONE | "The mechanism, first" is the file's first section, before the prose shape, and names the v2.3.0 hazard. |
| 3 | `Require curated GA release notes` step in `check-release`, immediately after `Detect release`, gated on `is_release == 'true' && is_prerelease == 'false'` | Step 2 | DONE | `release.yml:148-160`, directly after `Detect release` (ends `:146`). Gate is the plan's string verbatim. |
| 4 | The check's condition rejects missing / empty / whitespace-only | Step 2 | **DEVIATED (declared)** | One `grep`, not `[ ! -s ] \|\| grep`. Executed both forms over six arms — missing, empty, whitespace-only, bare newline, directory, with content — and they agree on every one. See Deviation 1. |
| 5 | Failure message names the exact path and points at the README | Step 2 | DONE | `::error::GA releases require curated notes at $NOTES (see .github/release-notes/README.md)`. |
| 6 | `notes_path` written to `$GITHUB_OUTPUT` by that step | Step 2a | DONE | `release.yml:159`, step carries `id: notes`. |
| 7 | `notes_path` declared on `check-release`'s `outputs:` | Step 2a | DONE | `release.yml:75`. |
| 8 | Step 3 consumes the output rather than rebuilding the path | Step 2a, Step 3 | DONE | `release.yml:419` is the only other occurrence of the path; `.github/release-notes/` appears nowhere else in the workflow. |
| 9 | `Create release` branches on `is_prerelease`; GA arm uses `--notes-file`, prerelease arm keeps `--generate-notes --prerelease` | Step 3 | DONE | `release.yml:405-420`, byte-equivalent to the plan's snippet (only an added two-line comment). |
| 10 | The replacement defines no shell variable — interpolation only | Step 3 | DONE | No `$VAR` in the step; `env:` holds `GH_TOKEN` alone. The defect §11 calls the finding that matters is absent. |
| 11 | `PRERELEASE_FLAG` dropped | Step 3 | DONE | Removed from `env:`; no occurrence anywhere in the workflow. |
| 12 | `.github/release-notes/v2.4.0.md` written to the new policy | Step 4 | DONE | 2,070 bytes, matching §14's figure exactly. Conformance below. |
| 13 | Carry v2.4.0's three first-time-contributor acknowledgements into the curated body, **or** leave v2.4.0 alone | Step 4 | **PARTIAL** | The file carries none, and the plan offered only those two ways out. The `gh release edit` is deferred, so nothing is lost yet — but the artifact is not in a state the plan permits applying. Undeclared in §14. See Gap 2. |
| 14 | `gh release edit v2.4.0 --notes-file` **not** run without explicit go-ahead | Step 4, §9.1 | DONE (correctly withheld) | Not run. §14 says so plainly. |
| 15 | A static CI assertion that the notes check sits in `check-release` | Step 5a | DONE | `ga_notes_check_lives_in_check_release`. |
| 16 | …and that no status function appears in the workflow | Step 5a | DONE | `release_workflow_has_no_status_functions` over `always()`, `cancelled()`, `success()`, `failure()`; grep confirms zero occurrences in `release.yml`. |
| 17 | The assertion runs in CI | Step 5a, §4 | DONE | `tests/integration_release_workflow.rs`; `ci.yml:76` runs `cargo test --release`. |
| 18 | Step 5b — the absent arm folded into the next GA prep | Step 5b / §13 Open-5 | **Correctly absent** | §13 Open-5 drops 5b explicitly ("drop 5b. 5a carries the proof. §7's single-commit convention stands"), with the reasoning that 5b forces a real prep through a deliberately red run and proves the arm that fails *safely*. Confirmed the plan says this. |

### §3 — policy components and exclusions

Checked against `.github/release-notes/README.md` (the policy's statement) and
`.github/release-notes/v2.4.0.md` (the worked example).

| # | Item | Status | Notes |
|---|------|--------|-------|
| 19 | Component order: framing → 3–6 highlights → CHANGELOG link → install block, nothing else | DONE | README presents them in that order; `v2.4.0.md` follows it. |
| 20 | Highlight format `- **Headline** ([#NNN](url)) — one sentence.` | DONE, with one inconsistency | `v2.4.0.md` has exactly 6 highlights, each one sentence, each with the bold headline, the bracketed number and the em-dash. The README's template writes `/pull/NNN`; all six links in the backfill are `/issues/NNN`. See Observation 3. |
| 21 | The install block is **defined** here, not referenced | DONE | README defines the three lines; §3's "Revision 1 said *the existing* lines" correction is honoured. `cargo install trim-galore` matches `README.md:46`; the ghcr path matches what `IMAGE_NAME` lowercases to. |
| 22 | Exclusions: full PR list, migration instructions, rationale, per-caller mappings | DONE | All four named in README's "What a body excludes"; none present in the backfill. |
| 23 | 1,500–2,500 chars, stated as a sanity range and not a gate | DONE | README says exactly that; the backfill is 2,070. Nothing in the workflow enforces a length. |

### §3 — edge cases (six rows)

| # | Row | Status | Evidence |
|---|-----|--------|----------|
| 24a | GA, notes file absent → `check-release` fails, nothing published, message names the path | DONE | `grep` exits 2 on a missing file, `!` inverts, `exit 1`. Executed. Job 4 is gated on `is_release == 'true'` and every job names `check-release` in `needs`, so nothing downstream runs. |
| 24b | GA, present but empty or whitespace → same failure | DONE | Executed all four shapes: empty (exit 1), spaces+tabs (exit 1), bare newline (exit 1), directory (exit 2) all fail; content passes. This row was flagged as a likely miss; it is correct. |
| 24c | Prerelease → `--generate-notes`, unchanged, no notes file required or read | DONE | Step's `if:` excludes `is_prerelease == 'true'`, so it is skipped; Create release's prerelease arm reaches `--generate-notes --prerelease`. |
| 24d | Dry run → the check still runs | DONE | `Detect release` sets `is_release=true` for every `workflow_dispatch` irrespective of `inputs.dry_run` (`:103`), and the notes step has no dry-run gate. This row was flagged as a likely miss; it is correct — a GA dry run from `master` does require the file. |
| 24e | Re-dispatch after a partial ship → unchanged by this plan | DONE, premise narrowed | Nothing in the diff touches `tag_preexists` or the idempotency guard. But the row's reasoning ("the notes file being committed means a re-dispatch reads the same body") holds only for releases cut *under* the policy. See Observation 5. |
| 24f | Body edited on GitHub afterwards → permitted and untracked | DONE | Nothing reads the body back; README's "What these files are not" states it. |

### §4 — signature

| # | File | Status |
|---|------|--------|
| 25a | `.github/workflows/release.yml` — one new step, one new output, one rewritten step | DONE — exactly that, nothing else changed in the file by `8f984ac` |
| 25b | A CI assertion for 5a | DONE — `tests/integration_release_workflow.rs` (§4 left the venue open) |
| 25c | `.github/release-notes/README.md` | DONE (content PARTIAL, item 1) |
| 25d | `.github/release-notes/v2.4.0.md` | DONE (content PARTIAL, item 13) |
| 25e | No Rust changes to the product | DONE — the only Rust added is the test |

### §8 — assumptions, re-checked against the code as it stands

| # | Assumption | Plan's mark | Holds today? |
|---|---|---|---|
| 26a | `check-release` can see the repo | Verified | **Yes.** `release.yml:76-80` checks out before `Detect release`, same job; the notes step needs nothing of its own. |
| 26b | `gh release create --notes-file` rejects a missing file | **Unverified** | Still unverified, and the code does not depend on it. §14's "Not covered" says so. No contradiction. |
| 26c | Prereleases keep auto-notes | Policy choice | **Honoured.** |
| 26d | One file per GA, never rewritten | Policy choice | **Honoured** and stated in the README. |
| 26e | No contributor-facing release checklist to update | Verified | **Still true** — no `CONTRIBUTING.md` at the root or under `.github/`, and no page under `docs/src/content/docs/` describes the release process. But see Observation 4: `release.yml`'s own preamble is a second de facto home, and it was not updated. |
| 26f | The file proven present in job 0 is still present in job 4 | Verified | **Yes, re-verified.** All six `actions/checkout` steps (`:76`, `:186`, `:262`, `:389`, `:430`, `:459`) still omit `ref:` — zero `ref:` lines in the file. §8.6's `cancel-in-progress` caveat is itself now remedied at `:52-54` (`57317f3`). |

## Gaps (detail)

### Gap 1 — Step 1's worked example

**Expected:** Step 1 — "The policy from §3 ('what a GA body contains'), the four-part shape, **and one
worked example**." §4 adds that `v2.4.0.md` "doubles as the worked example".

**Found:** `.github/release-notes/README.md` has a highlight-line *template* and an install-block
*template*, but no example body, and never mentions `v2.4.0.md`. Someone opening the directory finds
the example by noticing the sibling file, not by being sent there.

**Gap:** one clause in the README pointing at `v2.4.0.md` as the worked example, or a short inline
example. The second reading — that §4 discharges this requirement by making the backfill itself the
example — is defensible, which is why this is PARTIAL and not MISSING; but on that reading the README
should still name it, since the README is by §8.5 the rule's only home.

### Gap 2 — Step 4's first-time-contributor acknowledgements

**Expected:** Step 4, "One cost revision 1 missed" — v2.4.0's auto-generated body carries three
first-time-contributor acknowledgements, which `--notes-file` deletes. "Carry them into the curated
body, or leave v2.4.0 alone — the policy does not need a backfilled example badly enough to drop
someone's first contribution off a release page."

**Found:** `.github/release-notes/v2.4.0.md` contains no acknowledgements. The file was written, so
the second option (leave v2.4.0 alone) was not taken either.

**Gap:** the file is not in a state the plan permits applying. Before the deferred `gh release edit`,
either add the three acknowledgements (they are recoverable from the live body — `gh release view
v2.4.0 --json body`) or decide to drop the backfill. §14 does not record this as outstanding, so it
is the item most likely to be lost: the `gh release edit` is otherwise described as needing nothing
but go-ahead.

## Deviations

### Declared in §14 — all three checked

1. **Step 2's condition is one `grep`, not two tests.** Real, and accurately described.
   `grep -q '[^[:space:]]'` returns non-zero for a missing file (exit 2), an empty one and a
   whitespace-only one. Executed both forms over missing / empty / whitespace / bare-newline /
   directory / content: they agree on all six. The claim that the `!` is not the errexit-exempt form
   §13 Open-6 found is correct — it is inside an `if` condition.
   *One part of this edit is undeclared:* the implemented form adds `2>/dev/null`, which the plan's
   snippet does not have. It suppresses grep's own "No such file or directory" line (harmless — the
   `::error::` names the path) but would also mask a genuine read error behind the
   "require curated notes" message. Inconsequential; noted for completeness.
2. **5a is a Rust test, not a step in `ci.yml`.** Accurately described, but **not strictly a
   deviation** — §4 says "Where it lives (a test, or a step in `ci.yml`) is an implementation detail".
   The implementation took one of the two options the plan offered. §14 says as much itself
   ("§4 left the venue open"), so the label is generous rather than wrong.
3. **5a asserts *direct* `needs`.** Real, accurate and genuinely stronger. All eight non-job-0 jobs
   name `check-release` directly in `needs:`, so the test asserts the direct property; the plan's
   §5 wording says "transitively", which `every_other_job_needs_check_release` would reject.

### Undeclared

| Where | What |
|---|---|
| `release.yml:155` | `2>/dev/null` on the `grep`, as above |
| `README.md` (release-notes) | No worked example — Gap 1 |
| `v2.4.0.md` | No contributor acknowledgements — Gap 2 |
| `v2.4.0.md` vs README template | `/issues/NNN` against the README's `/pull/NNN` — Observation 3 |
| `v2.4.0.md` | The CHANGELOG link has no section anchor — Observation 2 |

## Observations (not plan gaps)

1. **Step 2a's structural claim is airtight.** `notes_path` is empty whenever the notes step is
   skipped, so the question is whether the GA arm of `Create release` can run with it empty. It
   cannot: job 4 requires `is_release == 'true'`, `Detect release` sets `is_prerelease` to exactly
   `true` or `false`, and when it is `false` the notes step either wrote the output or failed job 0.
   The `--notes-file ""` case is unreachable.

2. **The backfill's CHANGELOG link is to the page, not the section.** §3 component 3 says "a link to
   the CHANGELOG **section**". `v2.4.0.md` links
   `https://www.trimgalore.com/reference/changelog/` with the text "CHANGELOG, Version 2.4.0".
   v2.4.0 starts around a fifth of the way down a ~1,000-line page, and Starlight generates an id for
   `### Version 2.4.0 (Release on 9 November 2026)`, so an anchor was available.

3. **`/issues/` vs `/pull/`.** All six backfill links use `/issues/NNN`, which is what `CHANGELOG.md`
   uses for the same six numbers — so the backfill is consistent with the material it was trimmed
   from, and the README's `/pull/NNN` template is the odd one out. Worth settling in one direction,
   since the backfill is meant to be the example the next release copies.

4. **The release procedure has a second home that was not updated.** §8.5's premise is intact — no
   `CONTRIBUTING.md`, no docs page. But `release.yml:3-33` opens with a "Release invariant" comment
   whose "canonical sequence is: 1. Bump version … 2. Push master 3. Dispatch" is the release
   checklist in practice, and it does not mention writing `.github/release-notes/<version>.md`.
   Anyone following that comment hits the new job-0 failure without warning. One line in that block
   closes it.

5. **A re-dispatch of a pre-policy release now fails in job 0.** Cargo.toml is at `2.5.0` and there is
   no `v2.5.0.md`, so a GA dispatch from `master` today would be refused for missing notes. That is
   the policy working as designed and it fails safely — job 0 publishes nothing — but edge-case row
   24e's reasoning ("the notes file being committed means a re-dispatch reads the same body") assumes
   a release cut under the policy. For v2.5.0 and earlier, recovery via re-dispatch now requires
   writing the notes file first. Worth knowing before it is needed in a hurry.

## Test verification

`cargo test` from the crate root: **912 passed, 0 failed** — matching §14's figure exactly.

| Test | File | Status |
|---|---|---|
| `ga_notes_check_lives_in_check_release` | `tests/integration_release_workflow.rs` | PASS |
| `release_workflow_has_no_status_functions` | same | PASS |
| `every_other_job_needs_check_release` | same | PASS |
| `job_parsing_sees_the_real_jobs` | same | PASS |

The job parser resolves 9 jobs (`check-release` plus the eight the plan counts) and excludes the
two-space `push:` key under `on:`, which is what `job_parsing_sees_the_real_jobs` exists to prove.

One note on `release_workflow_has_no_status_functions`: it matches the four strings against the whole
file, comments included, so a comment *mentioning* `always()` would fail it. That is a false-positive
risk, not a coverage hole.

## Verdict

**INCOMPLETE — 2 items.** Neither is in the mechanism.

1. **Step 1 — the worked example.** `.github/release-notes/README.md` states the policy and the
   mechanism correctly but contains no worked example and does not point at `v2.4.0.md`, which §4
   nominates for that role.
2. **Step 4 — the contributor acknowledgements.** `.github/release-notes/v2.4.0.md` drops the three
   first-time-contributor acknowledgements that Step 4 requires be carried over, and §14 does not
   record this as outstanding. It must be settled before the deferred `gh release edit`, not after.

Everything else the plan specifies is present: Steps 2, 2a, 3 and 5a as written, all six edge-case
rows behaving as the table says (verified by execution, not by reading), all four §4 artifacts, and
every §8 assumption marked "verified" still true against the file today. Step 5b's absence is correct
and the plan says so in §13. Step 4's unrun `gh release edit` is correctly withheld, not a gap.
