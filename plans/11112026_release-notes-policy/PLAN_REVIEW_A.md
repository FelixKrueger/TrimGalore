# PLAN_REVIEW_A — curated release notes + workflow guard

**Reviewer:** A (independent; no shared state with Reviewer B)
**Plan:** `plans/11112026_release-notes-policy/PLAN.md`
**Repo state at review:** `dev` @ `e7ec1e6`, `origin/master` 8 commits behind `dev`, `Cargo.toml` = `2.5.0` on both, `.github/release-notes/` does not exist.
**Date:** 2026-11-18

---

## 0. What I verified vs. what I inferred

### Verified against source (read, counted, or executed)

| Claim | Where | Verdict |
|---|---|---|
| `create-tag-and-release` is job 4, `release.yml:362` | `:362` | ✅ exact |
| `Create tag` pushes the tag at `:372-376` | push at `:376` | ✅ exact |
| `Create release` at `:383-392`, `gh release create` at `:389` | | ✅ exact |
| Tag push (`:376`) precedes `gh release create` (`:389`) **in the same job** | | ✅ exact — the §2 hazard is real |
| `check-release` checks out at `:73-76`, `fetch-depth: 1`, before `Detect release` (`:78`) | | ✅ exact — assumption 1 holds |
| `IS_PRERELEASE` computed at `:84` | `[[ "$RAW_VERSION" == *-* ]]` | ✅ exact |
| `VERSION="v$RAW_VERSION"` at `:95` → path is `v2.5.0.md` | | ✅ exact |
| Six `actions/checkout` at `:73, :169, :245, :370, :403, :432`, **none** with `ref:` | `grep -n "ref:" → (none)` | ✅ exact |
| Dry-run gates: `Create tag` `:373`, `Create release` `:384`, `upload-binaries` `:399`, `publish-crate` `:425` | | ✅ exact |
| No `always()` / `!cancelled()` anywhere in `release.yml` → a failed `check-release` **skips every downstream job** | | ✅ verified (the 7 hits are all in `ci.yml`) |
| No `CONTRIBUTING.md` at root, `.github/`, or `docs/`; no release-process docs page | | ✅ verified — §8.5 holds |
| Astro content root is `docs/src/content/docs/` | `docs/astro.config.mjs`, dir listing | ✅ verified — §9.2's correction is right |
| `.github/release.yml` does not exist | | ✅ verified — §9.3 consistent |
| Step 2's two-part test: `[ ! -s ]` catches missing+zero-byte, `grep` catches whitespace-only, `\|\|` short-circuits | reasoned + shell semantics | ✅ correct |
| v2.4.0 body = 16,084 chars, 110 PR lines, 29 dependabot lines | `gh release view` | ✅ exact on all three |
| Every release body is pure `--generate-notes` output | v2.4.0 / v2.3.0 / v2.2.0 inspected | ✅ verified — and see **L4** |

### Verified-with-a-correction

- **§2's dependabot char count.** I measure 4,375 chars across the 29 dependabot lines, leaving 11,709; the plan says 4,647 / ~11.4k. The plan is internally consistent (16,084 − 4,647 = 11,437) so this is a measurement-convention difference, not an error in the argument. Immaterial — both figures support "filtering alone does not get there". No action.
- **Branch base is stale.** §Header says `dev @ b0a83b0 (= origin/master = the v2.5.0 tag)`. `dev` is now `e7ec1e6`, 8 commits ahead of `origin/master`. This matters only because Step 5's premise depends on it — see **L1**.

### Could not verify (no safe local test)

- **Assumption 2** — that `gh release create --notes-file <missing>` exits non-zero *before* any API call. Reading `gh`'s `release create` source from memory, `NotesFile` is `os.ReadFile`'d inside `RunE` ahead of `createRun`, so a missing file errors before the network. I am ~85% confident but did **not** execute it (the only honest test creates a real release). It is belt-and-braces either way, so this is not load-bearing — but §8.2 should say "unverified" rather than asserting it.

---

## 1. Logic review

### L1 (Critical) — Step 5's "absent arm is free right now" is false; the proof cannot be run after landing on `dev`

The plan's branch base is `dev`. Step 5 says to dispatch `dry_run: true` **from `master`**. A `workflow_dispatch` runs the workflow file *as it exists on the dispatched ref*. `origin/master` has no notes-check step (`git show origin/master:.github/workflows/release.yml | grep -c release-notes` → `0`), so a dispatch from `master` after landing Steps 1–3 on `dev` exercises the **old** workflow and comes back green. The expected failure never appears.

And the GA arm cannot be reached from `dev` at all. `Detect release` (`:86-93`) allows `dev` **only** when `IS_PRERELEASE = true`; `dev`'s `Cargo.toml` is `2.5.0` (no hyphen) → `IS_PRERELEASE=false` → the job hard-errors at `:91` with "Release must be from master (GA) or dev (prerelease only)" **before** the new step runs. Bumping `dev` to a prerelease version to get past that gate sets `is_prerelease=true`, which *skips* the new step. There is no dispatch from `dev` that reaches the GA arm.

So the only ref that can exercise this check is literally `master`, which means the change must be on `master` first. Two consequences the plan does not price:

1. Fast-forwarding `master` outside a release contradicts the project's branch model (`master` = release pin, FF from `dev` only at release time) and would put 8 commits of unreleased work — including `--fastqc_raw` (#576) — on the release pin.
2. FF-ing `master` pushes commits carrying `Closes #N` onto the **default branch**, which is the only branch where GitHub's auto-close fires. Issues referenced in those 8 commits would close prematurely.

**Recommended fix:** drop the "free right now / must be run first" framing and relocate Step 5 into the next GA release prep, where FF-ing `master` is already the correct move:

1. Land Steps 1–3 on `dev` (with Step 4's `v2.4.0.md` as the worked example).
2. At the next GA: bump version + CHANGELOG, FF `master` as usual, but **do not yet commit `vX.Y.Z.md`**.
3. Dispatch `dry_run: true` from `master` → expect the refusal naming `.github/release-notes/vX.Y.Z.md`.
4. Commit the notes file, dispatch `dry_run: true` again → green.
5. Dispatch for real.

This is strictly better proof: it exercises the real release's version path, needs no branch-model deviation, and the "absent" window reopens at **every** future GA rather than closing forever. It also deletes the need for a throwaway `v2.5.0.md` (see **L3**).

### L2 (Critical) — the dry run gives *zero* coverage of Step 3, which is the step that sits past the tag push

This is the "validation that cannot fail on the mistake it exists to catch" pattern.

`Create release` is gated on `is_dry_run != 'true'` (`:384`). So in Step 5's second dispatch, "green through `create-tag-and-release`" means **every step in that job was skipped** — `Create tag`, `Reuse existing tag`, and `Create release` alike. The dry run proves Step 2 and proves nothing whatsoever about Step 3.

That is the wrong half to leave uncovered. Step 2's failure mode (missing file) is caught in job 0 and costs nothing. Step 3's failure mode — a typo'd path, an unexpanded variable, a malformed `if` — fires in job 4 **after** `Create tag` has already pushed (`:376` is the first step in the job), which is exactly the v2.3.0 partial-ship shape the plan exists to prevent. The plan moves the *file-existence* check out of the fault line and then introduces *new shell code* into it with no test.

This is not hypothetical, because the plan's own Step 3 snippet is broken — see **L5**.

Worse, #333's idempotency does **not** recover it. After a Step-3 bug strands a pushed tag, fixing the bug requires a commit on `master`, which moves `HEAD`; the next dispatch then hits `TAG_SHA != HEAD_SHA` at `:124-127` and hard-errors ("refusing to release from mismatched refs"). Recovery requires deleting the pushed tag by hand.

**Recommended fix** — make the dry run actually exercise the GA arm. Add a dry-run companion to `Create release` that resolves the same variables and reads the same file without calling the API:

```yaml
      - name: Preview release (dry run)
        if: needs.check-release.outputs.is_dry_run == 'true'
        env:
          VERSION: ${{ needs.check-release.outputs.version }}
          NOTES_FILE: ${{ needs.check-release.outputs.notes_file }}
        run: |
          if [ "${{ needs.check-release.outputs.is_prerelease }}" = "true" ]; then
            echo "Would run: gh release create $VERSION --generate-notes --prerelease"
          else
            echo "Would run: gh release create $VERSION --notes-file $NOTES_FILE"
            test -s "$NOTES_FILE" || { echo "::error::$NOTES_FILE unreadable in job 4"; exit 1; }
            head -5 "$NOTES_FILE"
          fi
```

This turns a wholly-skipped job into one that proves the GA arm's path resolves, the file is readable *from job 4's checkout*, and the branch selects correctly. It also makes §8.6 (same-commit checkout) an *asserted* property rather than an argued one. Mirrors the `Dry-run manifest preview` step already at `:336-342`, so it matches the workflow's existing idiom.

### L3 (Important) — Step 4/5 manufactures a `v2.5.0.md` that is a test artifact pretending to be a record

Step 5.3 writes `.github/release-notes/v2.5.0.md` solely so the second dispatch goes green. But v2.5.0 is already tagged and released with a `--generate-notes` body, and Step 4 backfills **v2.4.0 only**. The result is a file in a directory whose stated purpose (assumption 4) is "a record of what was announced" that records something never announced.

Three ways out, in order of preference:

1. **Adopt L1's relocation** — then no `v2.5.0.md` is needed at all and the problem evaporates.
2. Backfill **v2.5.0** instead of, or in addition to, v2.4.0. It is the release users see *now*, the file has to be written anyway under the current plan, and `gh release edit v2.5.0 --notes-file` is one extra command. Backfilling the older release while leaving the current one at 16k-style auto-notes is the odd choice.
3. Keep it and have the README state explicitly which files are post-hoc reconstructions.

The plan picks none of these; as written it silently produces outcome 3 without the disclosure.

### L4 (Important) — §3's "the install block — the existing cargo / conda / docker lines" does not exist

§2 states, and I confirmed, that every body in the project's history is raw `--generate-notes` output. I read v2.4.0, v2.3.0 and v2.2.0 in full: each is exactly `## What's Changed` → `## New Contributors` → `**Full Changelog**`. **There is no install block in any release body.** The word "existing" in §3.4 is wrong.

The cargo/conda/docker lines do exist — in `README.md:39-94` (`cargo install trim-galore`, `conda install -c bioconda trim-galore`, `docker run … ghcr.io/felixkrueger/trimgalore:latest`). So the component is authorable, but the plan never supplies its text, while §4 declares the README "the only place the rule will live" and §3 says "In order, nothing else". A mandatory, unspecified section in a policy whose sole home is one file is a gap that will be filled differently by whoever writes the next release.

**Fix:** §3.4 should name `README.md`'s Installation section as the canonical source and pin the exact three lines in Step 1's README as part of the worked example.

### L5 (Important) — Step 3's snippet is not runnable as written

```yaml
gh release create "$VERSION" --target "$SHA" --title "$TITLE" --generate-notes --prerelease
```

`$VERSION`, `$SHA` and `$TITLE` are **never assigned**. The current step (`:389-392`) inlines `${{ needs.check-release.outputs.version }}`, `${{ github.sha }}` and a literal title. Pasted as-is, this runs `gh release create "" --target "" --title ""`. Also missing from the snippet: the `env: GH_TOKEN: ${{ github.token }}` block (`:385-386`) without which `gh` is unauthenticated.

Secondary: §Step 3 says "Keeping `PRERELEASE_FLAG` for the prerelease arm only", but the snippet hardcodes `--prerelease` and never references `$PRERELEASE_FLAG`. Pick one — a retained-but-unused `env:` entry is dead config.

Given **L2**, a non-runnable snippet in the one step that executes past the tag push deserves to be corrected in the plan rather than left to implementation.

### L6 (Important) — nothing ties Step 2's path to Step 3's path

Step 2 builds `.github/release-notes/${{ steps.check.outputs.version }}.md`; Step 3 builds `.github/release-notes/$VERSION.md`. Today those resolve identically, so the plan is correct that the two jobs "cannot disagree about whether the file exists" — but they *can* disagree about **which file**. Any later edit to one string and not the other yields: job 0 passes, tag pushed, job 4 fails. The guard would then be actively misleading.

**Fix:** emit the path once. Add `notes_file` to `check-release`'s `outputs:` (`:63-70`), set it in `Detect release` alongside `version`, consume `${{ needs.check-release.outputs.notes_file }}` in job 4. One source of truth, and it removes a whole class of drift for about four lines.

### L7 (Optional) — the policy silently drops `## New Contributors`, and the backfill deletes it

v2.4.0's body credits three first-time contributors (@BenjaminDEMAILLE, @life4sucks, @sujeito-operator). §3's "In order, nothing else" removes that section from GA bodies, and Step 4's `gh release edit v2.4.0 --notes-file` would **delete** an acknowledgement that has been public since v2.4.0 shipped and that those three may have been pointed at. For a project that takes external contributions, that is a social cost, not a formatting one.

Suggest: either add "first-time contributors, if any" as a fifth policy element, or at minimum carry the existing `## New Contributors` lines through in the v2.4.0 backfill specifically. The same argument applies weakly to `**Full Changelog**: …/compare/A...B`, which is the only machine-navigable diff in the body and costs one line to keep.

### L8 (Optional) — `cancel-in-progress: true` remains a live partial-ship vector

`concurrency: group: release-${{ github.ref }}, cancel-in-progress: true` (`:51-53`) means a second dispatch from `master` cancels the first. If the first is between `:376` (tag pushed) and `:389` (release created), that is a partial ship with no failed step to point at. Entirely out of this plan's scope and correctly not claimed otherwise — but since the plan's framing is "a mechanism that cannot reintroduce a partial-ship hazard", one line in §7 noting this as a *separate, unaddressed* vector would keep a future reader from concluding the hazard class is closed.

### L9 (Optional) — the guard itself has no regression test

Once `v2.5.0.md` (or the next GA's file) exists, the "absent" arm is untestable until the next GA. A change that breaks the guard *permissively* — inverted `if:`, a path that always matches, a dropped step — produces no symptom until a GA ships with an empty body. Everything else in this repo that matters is covered by `cargo test` or a CI job.

Cheapest fix: extract the check into `.github/scripts/require-release-notes.sh` taking the path as `$1`, call it from the workflow, and add a two-case job in `ci.yml` (fixture present → exit 0, fixture absent → exit 1). That is also the only form of this guard that can be tested without dispatching a release at all, and it would have caught the class of bug in **L5** too.

---

## 2. Assumptions

| # | Plan's assumption | My verdict |
|---|---|---|
| 1 | `check-release` can see the repo | ✅ **Verified** — `:73-76`, `fetch-depth: 1`, before `Detect release` at `:78`. Full tree at depth 1; no sparse-checkout. |
| 2 | `gh release create --notes-file` rejects a missing file non-zero | ⚠️ **Unverified** (see §0). Reasoning supports it; I did not execute it. Label it unverified in §8.2 — the plan currently asserts it. Not load-bearing. |
| 3 | Prereleases keep auto-notes | ✅ Policy choice, internally consistent. Verified there is no reachable dispatch that wrongly catches a prerelease: the gate uses the same `is_prerelease` output that drives `PRERELEASE_FLAG`. |
| 4 | One file per GA, never rewritten | ⚠️ Contradicted in practice by Step 5.3's `v2.5.0.md` — see **L3**. |
| 5 | No contributor-facing release checklist exists | ✅ **Verified** — no `CONTRIBUTING.md` anywhere, no release docs page. One nuance: `release.yml:3-33` *is* a de-facto runbook ("Release invariant", the canonical 3-step sequence). A one-line pointer from that comment block to `.github/release-notes/README.md` would put the rule where someone cutting a release actually looks. |
| 6 | Job 0 and job 4 read the same commit | ✅ **Verified, and stronger than stated.** All six checkouts omit `ref:`; with `ref` empty on the workflow's own repo, `actions/checkout` sets `result.commit = github.context.sha` and checks out that commit, not the branch tip. Additionally the action is **pinned by commit SHA** (`3d3c42e5…`), so this behaviour cannot change without a deliberate bump — worth adding, since it is what makes the assumption durable rather than merely currently-true. |

### Unstated assumptions I found

- **A1.** That a dispatch can reach the GA arm from the branch the work lands on. It cannot — **L1**.
- **A2.** That "green through `create-tag-and-release`" means job 4's steps *ran*. They are all skipped in a dry run — **L2**.
- **A3.** That the install block already exists somewhere in release bodies — **L4**.
- **A4.** That a failed `check-release` skips everything downstream. True, but only because no job uses `always()` / `!cancelled()`. ✅ Verified today; worth a comment in the workflow, since adding `!cancelled()` to any downstream job later would silently reopen the fault line.
- **A5.** That the dry run is safe despite `v2.5.0` already being tagged. ✅ Verified, for a reason the plan does not name: besides the step gates it cites, the idempotency block at `:122` is itself guarded by `inputs.dry_run != "true"`, so a dry run never evaluates `TAG_SHA`/`HEAD_SHA` and cannot hard-error on a moved `master`.
- **A6.** That the check's content test is sufficient. `[ ! -s ] || ! grep -q '[^[:space:]]'` passes a file containing `TODO`. Since the filename is version-keyed, a stale copy from a prior release cannot masquerade; the realistic residual failure is a placeholder. The plan's "1,500–2,500 chars is a sanity range, not a gate" is a deliberate and, I think, correct call — over-gating prose creates its own refused-release failure mode. Noting it only so the choice is on the record.
- **A7.** `${{ steps.check.outputs.version }}` is interpolated into a shell string. The value comes from `Cargo.toml` via `grep`/`sed`, so injection requires landing a crafted `Cargo.toml` on `master` — which already implies workflow control. Consistent with the file's existing style (it inlines `${{ }}` everywhere). Flagging for completeness only; the `env:`-indirection form suggested in **L2** is the hardened shape if you want it.

---

## 3. Validation sufficiency

Taking each validation and asking whether it can fail on the mistake it exists to catch:

| Validation | Mistake it targets | Can it fail? |
|---|---|---|
| Step 5 dispatch #1 (file absent) | Guard does not fire on GA | ❌ **Not as sequenced** — runs `master`'s old workflow, comes green (**L1**). ✅ once the step is on `master`. |
| Step 5 dispatch #1, "`git ls-remote --tags` unchanged" | Check sits before the tag push | ❌ **Cannot fail.** `Create tag` is gated on `is_dry_run != 'true'` (`:373`), so no dry run ever pushes a tag — the assertion holds identically if the check lived in job 4 and failed there. The ordering property is static (job `needs` graph + no `always()`), proven by reading, not by this dispatch. The plan should say so instead of crediting the dry run with it. |
| Step 5 dispatch #2 (file present) | Step 3's `--notes-file` arm is wrong | ❌ **Cannot fail** — `Create release` is skipped under `is_dry_run` (`:384`). Zero coverage of the one step past the tag push (**L2**). |
| Step 2's `[ ! -s ] \|\| ! grep -q` | Missing / empty / whitespace-only file | ✅ Correct for all three. Short-circuit verified. |
| Assumption 2 (`gh` rejects missing file) | Belt-and-braces backstop | ⚠️ Untested, and under **L6** the two jobs could point at *different* paths, which this would catch only as a post-tag-push failure. |
| §3 edge case "nothing tagged, built, pushed or published" | Guard fires but something still ships | ✅ **Verified** — `needs: [check-release]` on every job, no `always()`/`!cancelled()` in `release.yml`. |

**Net:** two of the plan's three dispatch-based assertions cannot fail on the mistake they name, and the uncovered step is the dangerous one. The plan's validation strength is the opposite of its risk distribution.

---

## 4. Alternatives

1. **Generate the body instead of gating on it.** Have `create-tag-and-release` build the body from the CHANGELOG section for `$VERSION` (the CHANGELOG already carries the "Edition" framing line and one-sentence entries — `CHANGELOG.md:938`, `:1054`). No new artifact, no new required step, single source of truth. Rejected implicitly by the plan, and I agree with rejecting it: the CHANGELOG section is a complete record, the release body is a selection, and no parser makes that choice. Worth one line in §9 recording that it was considered — it is the first thing a reader will propose.
2. **`notes_file` as a `check-release` output** — see **L6**. Recommend adopting.
3. **A dry-run preview step in job 4** — see **L2**. Recommend adopting; it is the only way the plan's own validation reaches Step 3.
4. **Extract the check to a script + `ci.yml` fixtures** — see **L9**. Recommend; cheapest durable coverage, and it also makes the guard testable without a dispatch.
5. **Draft release instead of a hard gate.** `gh release create --draft` for GA, publish by hand after reviewing the body. Removes the missing-notes failure mode entirely (a draft with a bad body is edited, not re-released) but adds a manual publish step to every GA and changes when `upload-binaries` can attach assets. Mention and reject explicitly — the plan's automated-gate choice is defensible, but "why not a draft" is the obvious alternative and is currently unanswered.
6. **Require the file at PR time, not dispatch time.** A `ci.yml` check: if the PR changes `Cargo.toml`'s version to a non-prerelease, require `.github/release-notes/v<new>.md` in the same PR. Catches the omission days earlier, on a branch, where fixing it is free — and is complementary to, not a replacement for, the job-0 gate. Strong candidate for a follow-up.

---

## 5. Action items

### Critical

1. **Re-sequence Step 5.** The "absent" arm cannot be run from `master` until Steps 1–3 are *on* `master`, and cannot be run from `dev` at all (hard error at `release.yml:91` before the new step). Replace "free right now / must be run first" with the next-GA-prep sequence in **L1**. Remove the claim that the window closes once a notes file is committed — it reopens at every GA.
2. **Give Step 3 some coverage.** Add the dry-run preview step from **L2** so the GA arm's path, variable expansion and file readability are exercised in job 4 *without* publishing. Without it, Step 3's first real execution is a live GA release, past the tag push, with no test behind it.
3. **Fix the Step 3 snippet** (**L5**): assign `VERSION`/`SHA`/`TITLE` (or keep the existing `${{ }}` inlining), restore the `GH_TOKEN` `env:` block, and resolve the `PRERELEASE_FLAG`-retained-but-unused contradiction.
4. **Drop the "no tag pushed" claim as proof of ordering** (§Step 5, §2). It cannot fail. State instead that ordering is guaranteed statically by the `needs` graph plus the absence of `always()`/`!cancelled()`, both verifiable by reading.

### Important

5. **Emit `notes_file` from `check-release`** and consume it in job 4 (**L6**) — removes path drift between the guard and its consumer.
6. **Resolve the `v2.5.0.md` artifact** (**L3**) — adopting item 1 dissolves it; otherwise backfill v2.5.0 too, or disclose it in the README.
7. **Specify the install block** (**L4**). §3.4's "existing" is false — no release body has ever contained one. Point at `README.md:39-94` and pin the three lines verbatim in Step 1's README.
8. **Decide on `## New Contributors`** (**L7**) before Step 4 deletes three public acknowledgements from v2.4.0.
9. **Mark assumption 2 unverified** rather than asserted.

### Optional

10. Add `notes_file` / `require-release-notes.sh` + `ci.yml` fixture coverage for the guard itself (**L9**).
11. One line in `release.yml:3-33`'s header comment pointing at `.github/release-notes/README.md`, so the rule is reachable from the de-facto runbook (assumption 5 nuance).
12. A comment on the downstream jobs noting that adding `always()`/`!cancelled()` would reopen the fault line (**A4**).
13. One line in §7 acknowledging `cancel-in-progress` as a separate, unaddressed partial-ship vector (**L8**).
14. Record in §9 why CHANGELOG-derived generation was rejected (alternative 1) and why a `--draft` GA was rejected (alternative 5).
15. Strengthen §8.6 by noting `actions/checkout` is SHA-pinned, so the `github.sha` behaviour cannot drift silently.
16. Consider the PR-time check (alternative 6) as a follow-up — it moves the failure from release day to review day.

---

## 6. Overall

The policy half is sound and the central architectural claim — validate in job 0, which publishes nothing, because job 4 pushes the tag at `:376` before creating the release at `:389` — is **correct and verified**, as is the §8.6 reasoning that makes job-0 placement safe (all six checkouts pin `github.sha`, and the action is SHA-pinned so that cannot drift).

The weakness is not the placement; it is that the plan treats "the file exists" as the whole risk. It moves that one check to safety and then adds untested branching shell into the exact step it identified as the fault line, with a validation plan that (a) cannot be executed as sequenced and (b) exercises none of that new code. Items 1–4 close that gap; items 5–9 close the specification gaps a future release-cutter would otherwise hit.
