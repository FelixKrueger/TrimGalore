# PLAN_REVIEW_B — curated release notes + workflow guard

**Reviewer B** (independent, fresh context). **Target:** `plans/11112026_release-notes-policy/PLAN.md`
**Date:** 2026-11-18. **Repo state at review:** `dev` @ `e7ec1e6`, `origin/master` @ `b0a83b0` (= `v2.5.0`), both `Cargo.toml` at `2.5.0`.

**Verdict:** the core design decision — validate in `check-release` (job 0), not in `create-tag-and-release` (job 4) — is **correct and verified**, and for a stronger reason than the plan gives. But Step 5's proof is partly **vacuous**, Step 3's snippet as written **reintroduces the exact post-tag-push failure the plan exists to prevent**, and Step 5 **cannot be executed at all from the branch the work will land on**. Three Critical items, all in the mechanism rather than the policy.

---

## 0. Verification ledger

I checked every cited line number and every structural claim against `.github/workflows/release.yml` and the repo. Separated by how I established it.

### Verified by reading the source (plan correct)

| Plan claim | Result |
|---|---|
| `check-release` checks out before `Detect release` (`:73-76`, step at `:78`) | ✓ exact |
| `VERSION="v$RAW_VERSION"` at `:95`, so path is `v2.5.0.md` | ✓ exact |
| Tag push (`:374-376`) precedes `gh release create` (`:389`) **in the same job** | ✓ exact — §2's hazard is real |
| `Create tag` `:373`, `Create release` `:384`, `upload-binaries` `:399`, `publish-crate` `:425` all gated on `is_dry_run` | ✓ exact |
| Six `actions/checkout` at `:73 :169 :245 :370 :403 :432`, none with `ref:` | ✓ exact |
| `is_prerelease` computed in `Detect release` (plan cites `:84`) | ✓ — `:83` defaults false, `:84` sets true on a hyphen |
| No `.github/release-notes/` directory exists | ✓ |
| No `.github/release.yml` exists (§9.3 is genuinely untouched ground) | ✓ |
| §8.5: no `CONTRIBUTING.md` at root or `.github/`, no docs page covers the release process | ✓ — `grep -rln "workflow_dispatch\|release workflow\|dry_run" docs/ *.md` returns nothing |
| §9.2: Astro content collection is rooted at `docs/src/content/docs/` | ✓ — only `docs/` exists under `content/` |
| §3's register: "Edition" framing lines exist in `CHANGELOG.md` | ✓ at `:938` (v2.3.0), `:1054` (v2.2.0), `:1104`/`:1716` (v2.1.0) |
| Step 5: `master` is at `2.5.0` with no notes file, so the "absent" arm is unset-up | ✓ on content; ✗ on reachability — see **C2** |

### Verified by execution

The §10 "verified sound" claim about Step 2's two-part test — I ran the truth table rather than re-reasoning it:

```
missing.md -> FAILS   empty.md -> FAILS   ws.md -> FAILS
nl.md      -> FAILS   ok.md    -> PASSES  dir.md -> FAILS
```

All six correct, including the two the plan does not mention (bare-newline file; path occupied by a directory). **The check logic is sound.** This is the one part of the plan I can state is right without qualification.

### Verified, and the plan is wrong or incomplete

| Claim | Finding |
|---|---|
| "job 4 of 6" | There are **9 jobs**. The plan follows the file's comment numbering (`0.`…`6.`, which skips `1a`/`3a`), so it is internally consistent, but "of 6" misreads as a job count — and the §8.6 list of "all six checkouts" sitting next to it invites conflating the two. Three jobs (`smoke-test-binaries`, `docker-merge`, `smoke-test-docker`) check out nothing. |
| §5 "A dry run reaches none of them" | True, but **understated to the point of changing Step 5's expected outcome** — see **C1**. `create-tag-and-release` is not merely gated step-by-step on a dry run; the whole job is unreachable. |
| Step 5 "the only one that proves the check sits before the tag push" | **False.** See **C1**. |
| §3 item 4, "the existing cargo / conda / docker lines" | There is no install block in any workflow or template; the cargo line lives at `README.md:46`. If §2 is right that every body is pure `--generate-notes` output, then no shipped body contains an install block and "existing" has no referent in the release bodies. One of §2 and §3.4 is wrong. One `gh release view v2.4.0` settles it; `gh` is sandbox-blocked here (`~/.config/gh/hosts.yml` denied), so **unverified**. |
| §2's body-size table (16,084 / 2,043 / 2,285 / 6,416 chars) | **Unverified** — same `gh` block. Not load-bearing for the mechanism. |

### Inferred (documented behaviour, not executed)

- **§8.6 — `actions/checkout` with no `ref` pins `github.sha`.** I could not read the action's source (pinned by digest, not vendored). But the claim is right per documented behaviour, and the workflow already stakes three other things on the same invariant: `--target "${{ github.sha }}"` (`:390`), `git tag "$VERSION"` with no commit argument (`:375`, tags HEAD), and the idempotency guard comparing `git rev-parse "$VERSION^{}"` to `git rev-parse HEAD` *in job 0* (`:118-119`). If the checkouts did not all land on `github.sha`, the idempotency guard added after the v2.3.0 incident would already be comparing the wrong commits. **High confidence — and note it is a pre-existing load-bearing invariant, not one this plan introduces.**
- **Skip propagation through `needs`** (basis of **C1**): a job whose `if` contains no status function is skipped when any `needs` dependency is skipped. Confirmable in 30 seconds against any past dry-run log (`gh run list --workflow release.yml`), which I could not do.

---

## 1. Logic review

### C1 (Critical) — Step 5's dry run proves materially less than the plan claims, and one of its two assertions cannot fail

Two compounding problems.

**(a) The "no tag pushed" assertion is vacuous.** Step 5 bullet 2 asserts `git ls-remote --tags` unchanged, and the plan calls the absent-file dispatch "the only one that proves the check sits before the tag push". It cannot prove that. `Create tag` is gated on `is_dry_run != 'true'` (`:373`), so **no dry run ever pushes a tag regardless of where the notes check sits**. Run the same dispatch with the check deliberately mis-placed inside job 4 and `ls-remote` is still unchanged. The assertion cannot distinguish the correct design from the mistake it exists to catch — the exact failure pattern the brief asked me to hunt.

What the dry run *does* prove, and it is worth having: the check fires on a GA dispatch, names the path, and fails job 0. That is the useful half.

What actually proves placement is the job graph, and it is stronger than any dry run: **all eight other jobs transitively `needs: check-release`** (verified at `:150 :213 :234 :291 :347 :365 :400 :426`), so a job-0 failure skips every one of them — including `docker-build`/`docker-merge`, which are the only two jobs with no `is_release` gate. The plan's edge-case row "Nothing tagged, built, pushed or published" is true *because of this*, and the plan never says so. Replace the vacuous assertion with the static argument.

**(b) `create-tag-and-release` is skipped entirely on a dry run, so Step 5 bullet 4's expected outcome is wrong.** `smoke-test-docker` carries `if: needs.check-release.outputs.is_dry_run != 'true'` (`:348`) → skipped on a dry run. `create-tag-and-release` lists it in `needs` (`:365`) and its own `if` (`:364`) contains no status function → the skip propagates → **job 4 is skipped, not green**. Step 5 bullet 4 says "expect green through `create-tag-and-release`". An implementer following that will see a skipped job and start debugging a non-problem.

The consequence that matters more: **no dry run has ever executed, or ever can execute, the `Create release` step.** Step 3 — the one change this plan makes inside the hazardous job — has zero dry-run coverage, now and permanently. Its first execution will be a live GA release, immediately after the tag push. That makes **C3** much sharper than it looks.

### C2 (Critical) — Step 5 cannot run from where the work lands

Step 5 requires a `workflow_dispatch` from `master`, because `Detect release` hard-errors on any other combination:

- `master` + no hyphen → GA, check fires. The only path that works.
- `dev` + no hyphen → `::error::Release must be from master (GA) or dev (prerelease only)` at `:91`, exit 1, **before** the new step.
- any other ref → same error.
- `dev` + hyphen → prerelease, gate skips the check.

`dev` is currently 8 commits ahead of `master` and also at version `2.5.0`, so a dry run from `dev` hits the hard error. The proof is therefore only executable once Steps 1–3 are **on `master`** — which, under this repo's branch model (`dev` canonical, `master` fast-forwarded only at release time), is a decision the plan never names. Per the standing rule, an FF of `master` is not something to initiate without being asked for it, and Step 5 silently depends on one.

Three ways out; pick one explicitly in the plan:

1. **FF `master` for a CI-only change** — flag it as its own go-ahead point, alongside Step 4's `gh release edit`.
2. **Defer Step 5 to the next release prep.** But this collides with §7: if the version bump, CHANGELOG section and notes file land as "a single commit", the notes file exists the moment `master` is FF'd and the absent arm is gone. To keep it, the next prep must bump and FF *first*, dispatch, then commit the notes file — i.e. §7's single-commit convention and Step 5's ordering are in direct conflict, and §7 wins by default because it is the convention someone will follow without reading this plan.
3. **Drop the dispatch proof** in favour of the static job-graph argument (which is what actually proves the property, per C1) plus a local truth-table test of the check. Cheapest, and loses nothing real.

My recommendation: 3, with 1 available if a live dispatch is wanted anyway.

### C3 (Critical) — Step 3's snippet, taken literally, fails after the tag push

The snippet uses `$VERSION`, `$SHA` and `$TITLE`. The existing step defines none of them — it interpolates `${{ }}` directly (`:389-392`) and sets only `GH_TOKEN` and `PRERELEASE_FLAG` in `env:`. GitHub's default shell is `bash -e`, not `-u`. I tested it:

```
would run: gh release create [] --notes-file [.github/release-notes/.md]
step exit=0
```

Undefined expands to empty, the step does not abort on the expansion, and `gh release create ""` then fails non-zero — **in job 4, after `git push origin "$VERSION"` has already run**. Tag pushed, no release: the v2.3.0 partial ship, reintroduced by the plan's own remedy. There is no actionlint or zizmor in `.github/` to catch it, and per C1(b) no dry run can.

The plan reads as near-literal implementation code elsewhere (Step 2's YAML is copy-ready), so this will be copied. Either add the three `env:` entries, or — better — keep the existing `${{ }}` style and make Step 3 the smallest possible diff:

```yaml
      - name: Create release
        if: needs.check-release.outputs.is_dry_run != 'true'
        env:
          GH_TOKEN: ${{ github.token }}
        run: |
          if [ "${{ needs.check-release.outputs.is_prerelease }}" = "true" ]; then
            gh release create "${{ needs.check-release.outputs.version }}" \
              --target "${{ github.sha }}" \
              --title "Trim Galore ${{ needs.check-release.outputs.version }}" \
              --generate-notes --prerelease
          else
            gh release create "${{ needs.check-release.outputs.version }}" \
              --target "${{ github.sha }}" \
              --title "Trim Galore ${{ needs.check-release.outputs.version }}" \
              --notes-file "${{ needs.check-release.outputs.notes_path }}"
          fi
```

Note `PRERELEASE_FLAG` drops out: §3 says "keeping `PRERELEASE_FLAG` for the prerelease arm only", but the snippet hard-codes `--prerelease`, leaving the `env:` entry dead. Decide one; dead release-workflow config is how the next person gets confused.

### I1 (Important) — compute the notes path once, in job 0, and pass it down

§8.2 and §8.6 argue that jobs 0 and 4 *cannot* disagree about the file. That argument is sound today (same commit, same `version` output, same literal path). It is sound by reasoning, not by construction: the path string is written twice, in two jobs, 300 lines apart. Change one — move to `docs/release-notes/`, add a `.md` → `.mdx`, insert a subdirectory — and job 0 passes while job 4 fails after the tag push.

Make it structural. `check-release` already exposes seven outputs; add an eighth:

```yaml
    outputs:
      notes_path: ${{ steps.check.outputs.notes_path }}
```

emitted by the same step that builds it, consumed by Step 3. One source, no possible drift, and §8.2's "belt-and-braces" becomes an actual brace instead of an argument. This also kills the only residual case in which `gh release create --notes-file` can fail post-tag-push for a path reason.

### I2 (Important) — §8.2's belt-and-braces is the wrong direction of brace

§8.2 rests the fallback on `gh release create --notes-file` rejecting a missing file. It does — *after* the tag push. A fallback that only fires once the damage is done is not a fallback. If anything is wanted inside job 4, it must be a `test -s` as the job's **first** step, before `Create tag`. Cheap, cannot hurt, closes the window that C3 and I1 both live in. With I1 applied, even this is belt-and-braces on a brace — but it is the only guard that sits on the right side of line 376.

### I3 (Important) — `v2.5.0.md`'s status is undefined, and contradicts §8.4

Step 5 bullet 3 writes `v2.5.0.md` "since the dispatch keys off `master`'s current `Cargo.toml` version". v2.5.0 shipped on 13 November with `--generate-notes` output as its body. §8.4 frames the directory as "a record of what was announced". A `v2.5.0.md` written solely to green a dry run is **not** what was announced, and nothing in the plan says whether `gh release edit v2.5.0` follows. The directory then holds three different kinds of file — a real record (future GAs), a corrected record (v2.4.0 after Step 4), and a test fixture (v2.5.0) — with nothing distinguishing them. Either backfill v2.5.0's body too (consistent, one more public edit), or drop Step 5's dispatch per C2 so the file is never needed. Do not leave it unstated.

### I4 (Important) — the concurrency group can strand a tag by cancellation

Not introduced by this plan, but §8.6 invites the wrong conclusion. `concurrency: group: release-${{ github.ref }}` with `cancel-in-progress: true` (`:49-51`): a dispatch from `master` and a push to `master` share the group `release-refs/heads/master`, so **a push to `master` during a release run cancels it**. §8.6 reassures that such a push "cannot remove the notes file between the check and its use" — true of the file, false of the run. Cancellation between `:376` and `:389` produces the identical partial ship. #333's idempotency makes re-dispatch recoverable, so this is a documentation gap in §2's hazard analysis rather than a defect, but §8.6 should not read as "a mid-run push is harmless".

### O1 (Optional) — v2.4.0 has no framing line to draw on

§3 item 1 cites the "Edition" register. Verified: `CHANGELOG.md` carries one for v2.1.0, v2.2.0 and v2.3.0 — **not for v2.4.0 or v2.5.0**. So Step 4's backfill has to invent its framing from scratch, which is the single most subjective part of the policy, and the plan neither names what v2.4.0's framing is nor flags that the worked example will be the only one not derived from existing CHANGELOG prose.

### O2 (Optional) — Step 4 erases what was actually announced

`gh release edit v2.4.0 --notes-file` overwrites the body with no history kept by GitHub. Given §8.4's "record of what was announced" framing, capture the old body first (`gh release view v2.4.0 --json body > …`) before editing. Ten seconds, and it is the only copy.

### O3 (Optional) — §3 item 3's CHANGELOG link needs a defined target

"A link to the CHANGELOG section" is underspecified in a way that matters for a README template. The canonical source is `CHANGELOG.md` (the docs copy at `docs/src/content/docs/reference/changelog.md` says so in its own note, and I found no sync script — it is hand-maintained). Either target produces a **date-bearing anchor** from `### Version 2.5.0 (Release on 13 November 2026)`, so the link cannot be a template constant and must be rebuilt per release. Say which target and give the anchor rule, or each release improvises.

### O4 (Optional) — interpolating `version` into a `run:` body

Step 2 puts `${{ steps.check.outputs.version }}` directly inside the script. The value comes from `Cargo.toml` on a ref only writers can dispatch, so this is hygiene, not a vulnerability — but `env:` is the convention, and the repo has no workflow linter to nudge it later.

### Gaps and edge cases §3 does not cover

- **GA dispatched from `master` with a hyphenated version** (a prerelease cut from `master` — permitted; `:86-87`'s `: # always allowed (GA)` comment is misleading, `master` allows anything). The gate correctly skips, `--generate-notes` applies, and Docker's `latest`/major/minor tags are already gated the same way (`:312-314`). **Consistent — no action, worth recording as checked.**
- **A notes file for a version that is never dispatched** (bump abandoned, version re-bumped). Harmless clutter; accumulates silently. Nothing needed, but the README should say the directory is append-only and that a stale file is not an error.
- **Manual `gh release create` outside the workflow** bypasses the policy entirely. Out of scope, but the README is the only home for the rule (§8.5) and that is where it would be said in a clause.

---

## 2. Assumptions

| # | Plan's status | Mine |
|---|---|---|
| 1 | Verified (checkout before `Detect release`) | ✓ confirmed, `:73-78` |
| 2 | "`--notes-file` rejects a missing file — belt-and-braces" | True but **on the wrong side of the tag push** (I2). With I1 the premise stops mattering. |
| 3 | Prereleases keep auto-notes — a policy choice | ✓ sound, and the gate implements exactly it (verified against all five reachable dispatch combinations) |
| 4 | One file per GA, never rewritten | **Contradicted by the plan's own Steps 4 and 5** — one is a rewrite of a shipped body, one is a fixture (I3) |
| 5 | Verified: no release checklist exists | ✓ confirmed |
| 6 | Verified: all checkouts pin `github.sha` | ✓ high confidence (documented behaviour + three in-repo dependencies). Does **not** cover run cancellation (I4). |

**Unstated assumptions worth promoting into §8:**

- **A7 — every job transitively depends on `check-release`.** This, not the dry run, is what makes job-0 placement safe. Verified. It is the plan's single most important fact and it appears nowhere.
- **A8 — Steps 1–3 will be on `master` when Step 5 runs.** Currently false and not addressable without an FF (C2).
- **A9 — the implementer will use the two code snippets verbatim.** Step 2's is safe to copy; Step 3's is not (C3).
- **A10 — `--notes-file` resolves relative to the workspace root.** True: job 4 checks out at `:370` with default working-directory. Verified, but it is why the relative path works and the plan does not say it.

---

## 3. Validation sufficiency

For each validation, can it fail on the mistake it exists to catch?

| Validation | Catches | Can it fail? |
|---|---|---|
| Step 2 check, GA + missing/empty/whitespace file | Shipping a GA with no curated notes | **Yes** — truth table executed, all six cases correct |
| Step 2 check, GA + file for the wrong version | Bumped the version, forgot to rename the notes file | **Yes** — path is version-keyed, so a mismatch reads as absent. The realistic human error is covered. |
| Step 2 check, GA + stale or placeholder content | Last release's notes copy-pasted | **No** — inherent; a content check cannot judge curation. Acceptable, but a one-line guard against the README's own template markers leaking through would cost nothing. |
| Step 5 dispatch 1, "fails naming the path" | The check does not fire on a GA dispatch | **Yes** — this is the useful half |
| Step 5 dispatch 1, "`git ls-remote --tags` unchanged" | The check sitting after the tag push | **No — vacuous.** Dry run suppresses the push unconditionally (C1a) |
| Step 5 dispatch 2, "green through `create-tag-and-release`" | Step 3's GA arm being broken | **No — doubly so.** Job 4 is skipped on a dry run (C1b), and even if it ran, `Create release` is `is_dry_run`-gated |
| `gh release create --notes-file` rejecting a missing file (§8.2) | Disagreement between jobs 0 and 4 | Fires, but only after the tag is pushed (I2) |

**Net:** one validation is genuinely sound and executed; one is vacuous; one is impossible; the fallback is mistimed. The plan's own "only test that matters" is the vacuous one.

**What would actually cover Step 3**, given that no dry run can reach it: nothing dynamic. So eliminate the failure modes by construction instead — I1 (single path output) plus C3 (no bare shell variables) together reduce Step 3 to a branch on one boolean, with every other token interpolated exactly as it is today. That is the only form of coverage available for a step the dry-run mode cannot reach.

---

## 4. Alternatives

1. **Let `gh release create` create the tag; delete the `Create tag` step.** `gh release create <tag> --target <sha>` creates the tag itself. One atomic operation, and the partial-ship fault line that motivates this plan's whole job-0 argument disappears. Trade-offs: a lightweight tag created by `gh` rather than a pushed git tag; interacts with #333's idempotency guard; and a failed create then leaves *nothing* (strictly better than a stranded tag). Out of scope here, but it is the structural fix to the hazard §2 describes, and worth filing separately.
2. **`--draft` first, publish last.** Create the release as a draft in job 4, publish it after `upload-binaries`. A mid-pipeline failure leaves a draft rather than a release with no binaries. Independent of this plan.
3. **Extract the notes from the CHANGELOG instead of committing a second file.** §7 already describes the notes file as "a trimmed restatement of the CHANGELOG section" — i.e. the plan accepts hand-duplication of material that already exists, with the usual drift risk. An `awk` extract between `### Version X.Y.Z` and the next `### Version` needs no new artifact and cannot drift. It is rejected for a good reason — the CHANGELOG section is a full record (v2.5.0's runs ~100 lines), not 3–6 highlights — but the plan should *say* that, because "why not just read the CHANGELOG" is the first question any reader asks. A middle path exists (a `<!-- highlights -->` marker inside the CHANGELOG section) and trades a second file for a parsing dependency; the plan's choice is defensible, just unargued.
4. **Make the directory `v<version>/` with `body.md`** rather than `<version>.md`, so future per-release assets (an image, a migration snippet) have a home. Speculative; the flat form is fine.

---

## 5. Action items

### Critical

1. **C3 — fix Step 3's snippet.** `$VERSION`/`$SHA`/`$TITLE` are undefined; under `bash -e` they expand to empty and the step fails *after* the tag push, recreating the v2.3.0 partial ship. Use the minimal-diff form above, or add the `env:` entries. Also resolve `PRERELEASE_FLAG`: keep it and use it, or delete it.
2. **C1 — rewrite Step 5's claims.** Drop "`git ls-remote --tags` unchanged" (vacuous — the dry run suppresses the push regardless of where the check sits) and "proves the check sits before the tag push" (it cannot). Replace with the static argument: all eight other jobs transitively `needs: check-release`. Correct bullet 4's expected outcome — `create-tag-and-release` is **skipped** on a dry run, because `smoke-test-docker` is skipped and the skip propagates. State plainly that Step 3 has no dry-run coverage and never will.
3. **C2 — resolve how Steps 1–3 reach `master`.** The dispatch proof is only executable from `master`; `dev` hard-errors at `Detect release` before the new step. Name the choice: FF `master` for a CI-only change (its own go-ahead point), defer to the next prep (and then reconcile with §7's single-commit convention, which destroys the absent arm), or drop the dispatch in favour of the static argument. Recommend the last.

### Important

4. **I1 — emit `notes_path` from `check-release` and consume it in job 4.** Converts §8.2/§8.6's *argued* invariant into a structural one and removes the only way the two jobs can disagree.
5. **I2 — if anything stays in job 4, make it a `test -s` as the job's first step**, before `Create tag`. §8.2's current fallback only fires after the tag is pushed.
6. **I3 — decide `v2.5.0.md`'s status** and reconcile it with §8.4. Either backfill v2.5.0's body too, or remove the need for the file.
7. **I4 — add the cancellation hazard to §2** and soften §8.6: `concurrency` + `cancel-in-progress` means a push to `master` mid-run cancels the release, which can strand a tag exactly as a step failure would.

### Optional

8. **O1** — say what v2.4.0's framing line is; it is the only backfill with no CHANGELOG "Edition" line to draw on.
9. **O2** — archive v2.4.0's current body before `gh release edit` overwrites it.
10. **O3** — define the CHANGELOG link target and its date-bearing anchor rule in the README template.
11. **O4** — pass `version` into Step 2 via `env:` rather than interpolating into the script body.
12. **§3.4 install block** — reconcile with §2: one `gh release view v2.4.0` settles whether shipped bodies contain install lines at all. I could not run it (`gh` sandbox-blocked).
13. **Alternative 1** — file the "let `gh release create` create the tag" idea separately; it removes the fault line this plan routes around.

### Verified sound, no change needed

The job-0 placement decision; Step 2's check logic (executed, six cases); the GA/prerelease gate across all five reachable dispatch combinations; §8.1, §8.5, §8.6 and §9.2's Astro correction; the relative `--notes-file` path resolving in job 4; dry-run safety of the pre-existing `v2.5.0` tag (the idempotency guard at `:117` is itself dry-run-gated, which §5 does not mention but which is the actual reason that dispatch is safe).
