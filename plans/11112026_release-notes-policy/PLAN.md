# PLAN — curated release notes, and a workflow that cannot ship without them

**Issue:** none filed (policy decision taken in session, 2026-11-11)
**Branch base:** `dev` @ `b0a83b0` (= `origin/master` = the `v2.5.0` tag; was `ee675c1` when drafted)
**Date:** 2026-11-11, manual review 2026-11-13
**Decision taken:** curated highlights only — drop `--generate-notes` in favour of a committed
notes file. Chosen over "highlights + folded full list" and "filter the auto list".
**Status:** **Implemented, 2026-11-18** (§14), except Step 4's public `gh release edit`, which the
plan reserves for explicit go-ahead. Revision 3 settled §12's three open questions (§13) and
corrected one claim the plan made twice about `actionlint`. Dual agent review covered revision 2
(`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`, see §11); **revision 3's own changes have not been
re-reviewed.**

---

## 1. Goal

Make a Trim Galore GA release body something a reader can use to decide whether the upgrade affects
them, and make the release workflow refuse to run a GA release that has no such body — early, before
anything is published.

Two separate deliverables: a **policy** (what the body says) and a **mechanism** (how it reaches
GitHub without reintroducing a partial-ship hazard).

## 2. Context

### What the bodies look like today

Every release body in the project's history is `gh release create --generate-notes` output — a flat
`## What's Changed` list of every merged PR, one line each, newest tag to this one. There is **no
curated precedent to restore**; this is new ground.

| Tag | Body | PR lines |
|---|---|---|
| v2.4.0 | 16,084 chars | 110 |
| v2.3.0 | 2,043 chars | 14 |
| v2.2.0 | 2,285 chars | 15 |
| v2.1.0 | 6,416 chars | ~45 |

v2.3.0 and v2.2.0 are short because the release was small, not because anyone curated them. So the
length is a function of release size, and v2.4.0 is what the current mechanism produces for a
release of that size — not an anomaly to be corrected by filtering.

### Filtering alone does not get there

Of v2.4.0's 110 lines, 29 are `dependabot[bot]` and account for 4,647 chars. A
`.github/release.yml` excluding bots would leave ~11.4k chars and 81 lines. The remaining lines
*are* the real work; there is no mechanical rule that separates "the poly-G/poly-A ordering fix" from
"tidied the `--help` text". Hence curation rather than configuration.

### The hazard the mechanism must avoid

`create-tag-and-release` (`.github/workflows/release.yml:362`) is job 4 of 6. Within it:

1. `Create tag` — **pushes the tag to origin** (`:372-376`)
2. `Create release` — `gh release create … --generate-notes` (`:383-392`)

and `publish-crate` is job 6, after `upload-binaries`. So a step that fails *inside*
`create-tag-and-release` leaves the tag pushed and no release — the shape of the v2.3.0 partial-ship
incident. A naive "read the notes file here, fail if absent" would sit exactly on that fault line.

Therefore: **the notes file is validated in `check-release` (job 0), which publishes nothing.** By the
time job 4 runs, the file's existence is already proven.

### Prereleases are a different case

Prereleases are cut from `dev` and have been days apart. Requiring a hand-written body for every beta
would be a tax on exactly the releases nobody reads as a summary. The policy applies to **GA only**;
prereleases keep `--generate-notes`. `check-release` already computes `is_prerelease` (`:84`), so the
distinction is free.

## 3. Behavior

### Policy — what a GA body contains

In order, nothing else:

1. **One line of framing.** What this release is. v2.2.0's "The Clumpify Edition" / v2.3.0's "The
   Formats Edition" line from the CHANGELOG is the register.
2. **Up to six highlights**, one sentence each, in the CHANGELOG's own voice — `- **Headline**
   ([#NNN](url)) — one sentence.` Selected as "a user would change what they do because of this".
   Dependency bumps, CI, docs and internal refactors are never highlights, **except a bump that
   clears a security advisory**. Three to six is the usual shape; write fewer when the release
   genuinely has fewer. **Both clauses amended 2026-11-18** — see
   `plans/11182026_v2.6.0-release/RELEASE.md` §1 for the release that forced them.
3. **A link to the CHANGELOG section** for the full record.
4. **Contributor acknowledgements**, when the release has first-time contributors, under GitHub's
   own `## New Contributors` heading and wording. **Added in the 2026-11-18 code review** (§15) —
   revision 3 shipped a four-component policy that had no slot for them, and `--generate-notes`
   was supplying this section for free. See §15 for why the list cannot be regenerated later.
5. **An install block** — cargo, conda and docker lines. **Note these do not exist to copy.**
   Revision 1 said "the existing ... lines"; no prior release body has any, because every one is pure
   `--generate-notes` output. The cargo line nearest to canonical is `README.md:46`. Step 1's
   `README.md` must therefore *define* this block, not reference it.

Explicitly excluded: the full PR list, migration instructions (they live in `docs/`), rationale,
per-caller mappings. The one-sentence-per-entry rule that governs `CHANGELOG.md` governs this too —
same reader, same question.

A GA body should land near 1,500–2,500 chars. That is a sanity range, not a gate.

### Mechanism

- GA notes live at `.github/release-notes/<version>.md`, committed on the release-prep branch
  alongside the version bump and the CHANGELOG section. One file per GA; they accumulate as a record.
- `check-release` gains a step: when `is_release == true` and `is_prerelease == false`, assert the
  file exists and is non-empty, and fail with a message naming the exact path if not.
- `create-tag-and-release` passes `--notes-file` for GA and `--generate-notes` for prereleases.

### Edge cases

| Case | Behaviour |
|---|---|
| GA, notes file absent | `check-release` fails. Nothing tagged, built, pushed or published. Message names the path to create. |
| GA, notes file present but empty or whitespace | Same failure. `gh` would accept an empty body silently, so the check tests for content, not existence. |
| Prerelease | `--generate-notes`, unchanged. No notes file required or read. |
| Dry run | The check still runs — a dry run that passes while the real release would fail on a missing file is a dry run that proves nothing. |
| Re-dispatch after a partial ship (#333 idempotency) | `gh release create` on an existing release fails as it does now; unchanged by this plan. The notes file being committed means a re-dispatch reads the same body. |
| Someone edits the body on GitHub afterwards | Permitted and untracked. The file is the source for *creation*, not a sync target. |

## 4. Signature

No Rust changes. Three files, plus one per release thereafter:

- `.github/workflows/release.yml` — one new step plus one new output (`notes_path`) in
  `check-release`, one rewritten step in `create-tag-and-release`.
- A CI assertion for 5a — that the notes check sits in `check-release` and that no status function
  appears in the workflow. Where it lives (a test, or a step in `ci.yml`) is an implementation
  detail; that it exists is the point, since it replaces the dispatch that could not prove it.
- `.github/release-notes/README.md` — new, stating the policy in §3 so the rule lives next to the
  files it governs rather than only in this plan. **This is the only place the rule will live** —
  see §8.5.
- `.github/release-notes/v2.4.0.md` — the backfill (Step 4), which doubles as the worked example.
- Thereafter one `.github/release-notes/vX.Y.Z.md` per GA, accumulating as a record.

## 5. Implementation outline

### Step 1 — `.github/release-notes/README.md`

The policy from §3 ("what a GA body contains"), the five-part shape (four as revision 3 had it;
§15 added the fifth), and one worked example. Short: this is a rule, not an essay.

It must also state the **mechanism**, not only the prose shape: that `check-release` fails a GA
dispatch when the file is missing or blank, and that the error names the path to create. There is no
`CONTRIBUTING.md` and no docs page covering the release process (§8.5), so this file is where someone
facing a refused release finds out why.

### Step 2 — validation step in `check-release`

After the existing `Detect release` step (`:78`), gated on
`steps.check.outputs.is_release == 'true' && steps.check.outputs.is_prerelease == 'false'`:

```yaml
      - name: Require curated GA release notes
        if: steps.check.outputs.is_release == 'true' && steps.check.outputs.is_prerelease == 'false'
        run: |
          NOTES=".github/release-notes/${{ steps.check.outputs.version }}.md"
          if [ ! -s "$NOTES" ] || ! grep -q '[^[:space:]]' "$NOTES"; then
            echo "::error::GA releases require curated notes at $NOTES (see .github/release-notes/README.md)"
            exit 1
          fi
          echo "Using curated release notes from $NOTES"
```

`check-release` already checks out the repo (`release.yml:73-76`, `fetch-depth: 1`), so the step
needs nothing else and can sit immediately after `Detect release`. `version` is emitted as
`v$RAW_VERSION`, so the path interpolates to `.github/release-notes/v2.5.0.md`.

### Step 2a — emit the path as a job-0 output

The step writes `notes_path=$NOTES` to `$GITHUB_OUTPUT` and `check-release` adds it to its `outputs:`
block. Step 3 then consumes it rather than rebuilding the path.

This is not tidiness. With the path constructed twice, job 4 can disagree with the job that
validated it — a rename, a directory change or a version-prefix change updates one and not the other,
and the failure lands *after* the tag push. One output makes "job 4 uses the file job 0 checked" a
structural fact instead of a convention two places happen to share.

### Step 3 — branch the `Create release` step

**Revision 2. The revision-1 snippet would have reintroduced the v2.3.0 partial ship**, which is the
incident this plan exists to prevent — both reviewers found it independently. It used `$VERSION`,
`$SHA` and `$TITLE`, which the real step never defines: it interpolates `${{ }}` directly and its
`env:` block holds only `GH_TOKEN` and `PRERELEASE_FLAG` (`release.yml:384-392`). GitHub runs
`bash -e`, not `-u`, so the undefined names expand to empty, the expansion exits 0, and
`gh release create ""` fails — **after** `git push origin` at `:376`. A tag on the remote, no release.
#333's idempotency does not recover it: fixing the bug moves `HEAD`, so the re-dispatch refuses at
`:124` on `TAG_SHA != HEAD_SHA`.

So the replacement defines nothing new and keeps every value in the interpolation the existing step
already uses. The notes path comes from a **job-0 output** (`notes_path`, §5 step 2a) rather than
being rebuilt here, so the path exists in exactly one place and job 4 cannot disagree with the job
that validated it:

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

`PRERELEASE_FLAG` is dropped: with an explicit branch it was only ever a way to avoid one, and
keeping it would mean two sources for the same decision.

**This step has no dry-run coverage and cannot acquire any** (§5 step 5), so it is the one piece of
new shell that reaches production unexercised. Two mitigations, not three: the interpolation-only
form above has nothing to be undefined, and a GA is always dispatched by hand, so the first real use
is watched.

**The third — `actionlint` — does not exist.** The claim that it "would catch the class outright"
was tested in the 2026-11-18 session and is false. `actionlint` prepends a preamble declaring the
environment before handing a `run:` block to shellcheck, so SC2154 (*referenced but not assigned*),
the check that would flag `$VERSION`/`$SHA`/`$TITLE`, can never fire on workflow shell. Verified both
ways: shellcheck reports five SC2154 hits on the revision-1 snippet as a standalone script and none
on the same snippet inside a workflow, with `-o all`, while a sibling optional check (SC2250) fires
in both — so the suppression is specific to SC2154 and is not a configuration oversight. Nothing in
CI covers this step. The interpolation-only form carrying no variables is the entire defence, which
is a reason to keep it that way.

### Step 4 — backfill v2.4.0

**Decided: in** (§9.1). Write `.github/release-notes/v2.4.0.md` to the new policy and
`gh release edit v2.4.0 --notes-file` it. Turns the policy's first example into a real one and fixes
the 16k body that prompted this. It is a public edit to a shipped release, so it is the one step here
that needs explicit go-ahead at the time.

Ordering: it is written **after** 5b's first dispatch and **before** its second, so it doubles as the
"file present" arm rather than a separate exercise.

**One cost revision 1 missed:** v2.4.0's auto-generated body carries three first-time-contributor
acknowledgements, which `--notes-file` would delete. Carry them into the curated body, or leave
v2.4.0 alone — the policy does not need a backfilled example badly enough to drop someone's first
contribution off a release page.

### Step 5 — what actually proves the placement

**Revision 2. Revision 1's dry-run proof was wrong in three independent ways** — both reviewers
found all three, and each is worth keeping on the record because each is a different kind of mistake.

1. **Its central assertion could not fail.** "No tag pushed" holds on *every* dry run, because
   `Create tag` is gated on `is_dry_run != 'true'` (`:373`). It would hold identically with the check
   in job 4, failing there — which is the arrangement it was written to rule out.
2. **Its expected outcome was impossible.** "Green through `create-tag-and-release`" cannot happen on
   a dry run: `smoke-test-docker` is dry-run-gated (`:348`), job 4 lists it in `needs`, and the
   workflow contains **no** status functions (`always()`/`!cancelled()`: 0 occurrences). So job 4 is
   **skipped**, not green.
3. **It could not be run from where the work lands.** A `workflow_dispatch` runs the workflow as it
   exists *on the dispatched ref*. `Detect release` hard-errors at `:91` for a non-prerelease version
   on `dev`, so the GA arm is reachable only from `master` — and `master` has no notes step until the
   FF, so dispatching there beforehand runs the *old* workflow and returns a meaningless green. The
   only way to run the intended proof early is to FF `master` outside a release, which puts
   unreleased commits on the release pin and prematurely fires every `Closes #N` among them.

**The placement does not need a dispatch to prove, because it is static.** All eight other jobs
transitively `needs: check-release` (`:150`, `:213`, `:234`, `:291`, `:347`, `:365`, `:400`, `:426`),
and with no status functions anywhere a failed job 0 skips all of them. That is a property of the
`needs` graph, readable in the file, and it is strictly stronger evidence than any single run.

So the proof is replaced by:

**5a — a static assertion, in CI.** A test that reads `release.yml` and fails if the notes check is
not in `check-release`, or if any status function appears in the workflow (which would let a
dependent job run after job 0 failed). This is the invariant, asserted where it cannot rot.

**5b — the absent arm, folded into the next GA prep.** At that point FF-ing `master` is correct
anyway, so the proof costs nothing extra. Order within the prep:

1. FF `master`, bump the version, do **not** write `vX.Y.Z.md` yet.
2. Dispatch `dry_run: true` → expect `check-release` to **fail**, naming the path.
3. Write the notes file, commit.
4. Dispatch `dry_run: true` → expect `check-release` to **pass** and jobs 1–3 to run.

Note what step 4 does and does not show: jobs 1–3 green and job 4 **skipped**, which is the correct
dry-run shape. Expecting anything from job 4 is what revision 1 got wrong.

**§7's "single commit" convention has to bend for this.** Writing the version bump and the notes file
in one commit makes the absent arm unreachable forever. Either the prep is two commits, or 5b is
dropped in favour of 5a alone — a real choice, now §10 Open-5.

## 6. Efficiency

Nothing hot. One `grep` in a job that already runs.

## 7. Integration

- **Release process.** Adds one required artifact to release prep. The version bump, the CHANGELOG
  section and the notes file become a single commit; the notes file is a trimmed restatement of the
  CHANGELOG section, so it is written *from* material that already exists.
- **`docs/`.** No change needed — there is no contributor-facing release checklist to update
  (§8.5). The consequence is that Step 1's `README.md` carries the rule alone.
- **No interaction** with the tag-idempotency work in #333, the Docker or crates.io jobs, or the
  branch model (GA from `master`, prereleases from `dev`).

## 8. Assumptions

1. ~~**`check-release` can see the repo.**~~ **Verified** — `release.yml:73-76` checks out at
   `fetch-depth: 1` before `Detect release`. Step 2 needs no checkout of its own.
2. **`gh release create --notes-file` rejects a missing file** with a non-zero exit. **Unverified —
   the only honest test creates a real release.** Worth knowing it is a weak fallback regardless: it
   fires *after* the tag push at `:376`, so it converts a silent wrong body into the partial-ship
   failure. The job-0 check plus Step 2a's single path is the real guard; this is not a second line
   of defence so much as a louder way to fail.
3. **Prereleases keep auto-notes.** A policy choice, not a constraint — revisit if beta bodies start
   being read as summaries.
4. **One file per GA version, never rewritten.** They are a record of what was announced, not a
   living document.
5. ~~**A contributor-facing release checklist may need the new step** (§7).~~ **Verified — there is
   none.** No `CONTRIBUTING.md` at the repo root or under `.github/`, and no docs page describes the
   release process. Nothing to update, but it makes Step 1's `README.md` the *sole* home of the rule,
   which is why Step 1 must state the mechanism and not only the prose shape.
6. ~~**The file proven present in job 0 is still present in job 4.**~~ **Verified, and it is what
   makes job-0 placement sound.** Every `actions/checkout` in `release.yml` (`:73`, `:169`, `:245`,
   `:370`, `:403`, `:432`) omits `ref:`. With `ref` empty on the workflow's own repository,
   `actions/checkout` sets `result.commit = github.context.sha` — the commit that triggered the run,
   not the branch tip. So `check-release` and `create-tag-and-release` read the *same commit*, and a
   push to `master` during the ~30-minute build cannot remove the notes file between the check and
   its use. **One caveat revision 1 implied away:** `concurrency` with `cancel-in-progress: true`
   (`:48-51`) means such a push *cancels* the run instead — possibly mid-job-4, between the tag push
   and the release creation. The file cannot change under the run, but the run can die inside the
   window; see §12's open questions. Without this, the §2 claim that the file's existence is "already proven" by job 4 would
   not hold, and a future reader may be tempted to move the check back into job 4 to close a gap that
   does not exist.

## 9. Resolved questions

All three settled at manual review, 2026-11-13.

1. **Step 4 (backfilling v2.4.0) — in or out?** → **In**, and promoted from "optional" to the
   mechanism's own proof: the backfill is what exercises Step 5's "file present" arm. Editing a
   shipped body is user-visible, so the `gh release edit` itself still needs go-ahead at the time,
   but the file is written either way as the policy's worked example.
2. **`.github/release-notes/` or `docs/release-notes/`?** → **`.github/release-notes/`**, for
   proximity to the workflow that reads it. Note the Astro half of the original argument was wrong:
   the content collection is rooted at `docs/src/content/docs/`, so a sibling `docs/release-notes/`
   would never be globbed as content and would need no ignore rule. Proximity is the whole case, and
   it is sufficient on its own.
3. **Should the prerelease arm also drop dependabot lines** via a `.github/release.yml`? → **Out of
   this plan**, as originally proposed. Independent and cheap, and not filed as an issue either: it
   is neither a data-integrity nor a user-visible defect.

## 10. Manual review findings

Raised 2026-11-13 and folded into the sections named. Recorded here so the agent reviewers can see
what has already been checked rather than re-deriving it.

| # | Finding | Landed in |
|---|---|---|
| 1 | Job-0 placement is sound for an unstated reason — all six checkouts pin `github.sha`, so jobs 0 and 4 cannot disagree about the file. The plan called Step 3 "belt-and-braces" without saying what against. | §8.6, §8.2 |
| 2 | §7's "verify whether a release checklist exists" resolves to *there is none* — no `CONTRIBUTING.md`, no release docs page. Step 1's `README.md` is therefore the only home for the rule, so it must state the mechanism, not just the prose shape. | §8.5, Step 1, §4 |
| 3 | Step 5's sequencing was unstated, and the "absent" arm — the plan's own "only test that matters" — is free only until a notes file for the current version is committed. Now an explicit four-step order. | Step 5, Step 4 |
| 4 | §4's signature said "two files" but the change is three plus one per release. | §4 |

**Verified sound, no change needed:** the §2 hazard (tag push at `release.yml:376` genuinely precedes
`gh release create` at `:389`, same job); assumption 1 (`:73-76` checks out before `Detect release`);
the `v`-prefix (`VERSION="v$RAW_VERSION"` at `:95`, so the path is `v2.5.0.md`); and Step 2's
two-part test — `[ ! -s ]` catches missing and zero-byte, the `grep` catches whitespace-only which
`-s` passes, and `||` short-circuits so `grep` never runs on a missing file.

---

## 11. Agent review round 1 — findings and fixes

Dual independent review, 2026-11-18 (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). **No contradictions**,
and every line number revision 1 cited was verified exact. Three Critical findings were found by both
reviewers independently.

### The one that matters

**Step 3's snippet would have reintroduced the v2.3.0 partial ship — the incident this plan exists to
prevent.** It used `$VERSION`/`$SHA`/`$TITLE`; the real step defines none of them. GitHub runs
`bash -e`, not `-u`, so they expand to empty, the expansion exits 0, and `gh release create ""` fails
*after* the tag push at `:376`. #333's idempotency cannot recover it, because fixing the bug moves
`HEAD` and the re-dispatch refuses at `:124`.

That a plan whose entire subject is "do not strand a tag" proposed a remedy that strands a tag is the
finding worth remembering. The mechanism was never run — it was read, by an author who already knew
what it was supposed to do.

### Fixed

| Finding | Where | Fix |
| --- | --- | --- |
| Step 3 undefined variables → partial ship | both | Rewritten using only the interpolations the existing step already uses; `notes_path` comes from a job-0 output |
| The path was built twice, so job 4 could disagree with job 0 | B | New Step 2a emits `notes_path` from `check-release` |
| Step 5's "no tag pushed" assertion cannot fail | both | Removed; replaced by 5a, a static assertion on the `needs` graph |
| Step 5 bullet 4 expected an impossible outcome — job 4 is **skipped** on a dry run, not green | B | Corrected, and the correct dry-run shape stated |
| Step 5's "absent arm is free right now" is false — a dispatch runs the workflow on the *dispatched ref*, and `master` has no notes step until the FF | A | 5b folded into the next GA prep, where FF-ing is already correct |
| §3's "the **existing** cargo/conda/docker lines" describes something that never existed | A | §3 now says the block must be *defined*, and names `README.md:46` |
| Backfilling v2.4.0 would delete three first-time-contributor acknowledgements | A | Called out in Step 4 as a cost, with two ways out |
| §8.2's `gh` fallback fires after the tag push | B | Recorded as a louder failure rather than a second line of defence; also marked unverified |
| §8.6 read as if a mid-run push is harmless | B | `cancel-in-progress` caveat added |

### The stronger argument the plan was missing

Both reviewers noted the placement does not need a dispatch to prove. All eight other jobs
transitively `needs: check-release`, and the workflow contains no status functions, so a failed job 0
skips every one of them. That is readable in the file and is better evidence than any run. It is now
5a, asserted in CI.

### Verified, not taken on trust

Every cited line number. §8.1, §8.5, §8.6 (verified and stronger than stated — the checkouts are
SHA-pinned, so the `github.sha` behaviour cannot drift). §9.2's Astro correction. The GA/prerelease
gate, across all five reachable dispatch combinations — neither reviewer found a reachable GA the
gate skips or a prerelease it wrongly catches. B executed Step 2's truth table rather than re-reasoning
it: missing, empty, whitespace, bare-newline and directory all fail; content passes.

### Still unverified

Assumption 2 — that `gh release create --notes-file` rejects a missing file. The only honest test
creates a real release. It now matters less, since it was demoted from a guard to a failure mode.

## 12. Open questions after review

**All three resolved in the 2026-11-18 session — see §13.** Kept here as the record of what was asked.

**Open-5 — 5b versus §7's single-commit convention.** Writing the version bump and the notes file in
one commit makes the "absent" arm unreachable forever. Either the prep is two commits, or 5b is
dropped and 5a carries the proof alone. 5a is the stronger evidence either way.

**Open-6 — add `actionlint` to CI?** There is no workflow linter in `.github/`, and it would have
caught Step 3's defect outright. Out of scope here, but this is the second workflow defect to reach
review; worth its own issue.

**Open-7 — `cancel-in-progress` and the tag-push window.** A push to `master` can cancel a release run
between `:376` and `:389`, producing the v2.3.0 shape with nobody making a mistake. Older than this
plan and independent of it, but its hazard analysis is where one would notice. Probably wants
`concurrency` scoped to exclude release dispatches.

## 13. Resolutions, 2026-11-18

**Line numbers.** Every `release.yml:` citation above this section predates the `actions-lint`
commit, which added two lines at `:333`. Anything the plan cites beyond that point has moved by two:
the tag push is now `:378` (was `:376`) and `gh release create` is now `:391` (was `:389`).
Citations below `:333` are unchanged.

### Open-5 — **drop 5b. 5a carries the proof. §7's single-commit convention stands.**

The plan framed this as a choice between two commits and dropping 5b, and said "5a is the stronger
evidence either way". Taking that at its word settles it, and two things decide which way:

**5b asks a real GA prep to pass through a deliberately invalid state.** The absent arm is only
reachable by dispatching the release workflow *between* the version bump and the notes file, which
means a failed release run sits in `master`'s history for every future audit of that release to
recognise and discount. The project has already been bitten once by a release whose history had to
be reconstructed after the fact (v2.3.0); adding a run that is *supposed* to be red is the wrong
direction.

**The two arms fail in opposite directions.** The arm 5a proves — the check sits in `check-release`,
so nothing can run past it — fails *dangerously*: a check in job 4 ships a tag with no release. The
arm 5b proves — the check fires when the file is missing — fails *safely*: the worst case is a GA
whose body is blank, which is embarrassing and immediately visible, not a partial ship. Spending
release-prep choreography on the safe arm while the dangerous one is covered statically is backwards.

So: 5a as specified, 5b dropped, and the version bump, CHANGELOG section and notes file stay one
commit as §7 says.

**Follow-up, not a requirement.** If the absent arm's *message* is later judged worth guarding, the
way to do it is to lift Step 2's condition into a script with its own test, which covers both arms on
every PR and needs no release-prep choreography at all. That is strictly better than 5b and can be
added whenever, so it is not a reason to hold the plan.

### Open-6 — **`actionlint` added. Its premise was wrong, and the plan said so twice.**

Added in the 2026-11-18 session as a separate `actions-lint` job in `ci.yml` — separate rather than
folded into `lint` so that job's status-check name stays as branch protection has it.

**It does not catch Step 3's defect.** See §5 Step 3 for the correction and the evidence; the short
version is that `actionlint` declares the environment to shellcheck before linting a `run:` block, so
SC2154 can never fire there. The sentence in Step 3 and the sentence in Open-6 were the same claim,
and both were wrong. Nothing in CI covers Step 3.

**What it does catch, on this repository, today.** Twenty-four findings, of which three were real and
twenty-one were quoting and style noise. The three:

| Where | What |
| --- | --- |
| `ci.yml:88` | `! ./trim_galore -V \| grep -q …` — a bare `!` is exempt from errexit, so the guard that `-V` carries no provenance line could not fail |
| `ci.yml:765` | `! command -v fastqc` — same defect, in the guard CLAUDE.md names as the proof that the bundled FastQC library did the work |
| `ci.yml:396` | a `case` arm `''` that cannot match a word holding two literal colons; the arm was also missing the trailing field |

Both inert guards were mutation-tested before and after: under the condition each exists to catch,
the current form exits 0 and the fix exits 1. The gate itself was mutation-tested the same way.

Two caveats worth keeping. **`actionlint` reported `ci.yml:765` but not `ci.yml:88`** — SC2251 fires
on the simple-command form and not the pipeline form, so the second was found by grepping for the
pattern the first revealed; the tool narrowed the search, it did not complete it. And **bare
`actionlint` exits 0 when discovery finds nothing**, which is a vacuous pass one wrong working
directory away, so the job globs the files explicitly and asserts the list is non-empty.

The configuration excludes `SC2086`, `SC2012` and `SC2129` **by code and never by severity**:
`SC2251`, which found the first inert guard, is only info-level, so any severity floor would have
hidden the one finding that justified the job.

No issue filed: the defects are fixed in the same session, and §Open-6's "worth its own issue" was
premised on actionlint being out of scope here, which it no longer is.

### Open-7 — **real, verified, and the exposure is wider than the plan's two line numbers.**

`release.yml:49-51` is `group: release-${{ github.ref }}` with `cancel-in-progress: true`. A
`workflow_dispatch` and a `push` on the same branch produce the same `github.ref`, hence the same
group, so a push cancels an in-flight release run on that branch. Confirmed by reading the trigger
block: `on.push.branches` is `[dev, master]`.

Three corrections to the plan's framing:

1. **The window is not `:376`–`:389`** (now `:378`–`:391`). It opens at the tag push and stays open through
   `upload-binaries` and `publish-crate`, the last job. A cancellation anywhere in that span leaves
   the tag on the remote with some channels unshipped. (`docker-merge` pushes to GHCR *before* the
   tag, so a cancellation before the tag leaves images without a tag — untidy, not a partial ship.)
2. **#333's idempotency does not recover it.** The push that cancelled the run also moved the branch,
   so a re-dispatch sees `TAG_SHA != HEAD_SHA` and hits the hard error at `:120-122` by design. Recovery
   needs the tag deleted by hand. This is the same shape §5 Step 3 describes for its own defect.
3. **The sharp case is a prerelease, not a GA.** GA is dispatched from `master`, which receives one
   FF per release. Prereleases are dispatched from `dev`, which takes routine merges, and the run is
   long — three binary targets, two Docker platforms, two smoke-test jobs. A merge to `dev` during a
   prerelease is an ordinary day, and it silently kills the release.

What the push-triggered run does in the meantime is nothing: `check-release` sets `is_release=false`
for any event that is not a `workflow_dispatch`, so every later job is skipped. A no-op run cancels a
publishing one.

**Remedy — applied, 2026-11-18.** Two changes, because the group alone leaves the defect reachable
from the second direction:

```yaml
concurrency:
  group: release-${{ github.ref }}-${{ github.event_name }}
  cancel-in-progress: ${{ github.event_name != 'workflow_dispatch' }}
```

The group keeps a push out of a dispatch's group. The predicate keeps a *second dispatch* from
cancelling a publishing first one — the group fix does nothing about that, and a double-dispatch when
nothing seems to be happening is an ordinary mistake with the same partial-ship outcome. Together a
dispatched release becomes uncancellable, while push-triggered runs still collapse redundantly, which
is the only thing `cancel-in-progress` was ever doing here.

A second dispatch now queues instead of cancelling. It then finds the tag at the same commit, takes
#333's idempotency path, and fails at `gh release create` because the release exists — a safe
failure against an already-complete release, not a partial one.

`actionlint` type-checks `concurrency.cancel-in-progress`, so the expression is validated in CI;
confirmed by mutation (a bogus `github.*` property is reported at that line).

No issue filed: fixed in the same session as Open-6.

## 14. Implementation notes, 2026-11-18

Implemented on `dev` from `ccae784`. Steps 1, 2, 2a, 3 and 5a are done; Step 4's file is written but
its `gh release edit` is not run.

| Step | What landed |
| --- | --- |
| 1 | `.github/release-notes/README.md` — the §3 policy, and the mechanism first, since §8.5 established this is the rule's only home |
| 2 | `Require curated GA release notes` in `check-release`, after `Detect release` |
| 2a | `notes_path` declared on `check-release`'s `outputs:`, written by that step |
| 3 | `Create release` branched on `is_prerelease`; `PRERELEASE_FLAG` dropped |
| 4 | `.github/release-notes/v2.4.0.md` written — 2,070 chars against the 16,084 the auto-notes produced. **`gh release edit v2.4.0 --notes-file` not run**: a public edit to a shipped release, reserved by §9.1 for go-ahead at the time |
| 5a | `tests/integration_release_workflow.rs`, four tests |

### Deviations

**Step 2's condition is one `grep`, not two tests.** The plan had
`[ ! -s "$NOTES" ] || ! grep -q '[^[:space:]]' "$NOTES"`. `grep -q '[^[:space:]]'` already returns
non-zero for a missing file (exit 2), an empty one and a whitespace-only one, so `-s` was redundant.
Checked against all four arms — missing, empty, whitespace-only, with content — and the two forms
agree on every one. The `!` is inside an `if` condition, where it is consumed explicitly, so it is
not the errexit-exempt form §13 Open-6 describes.

**5a is a Rust test, not a step in `ci.yml`.** §4 left the venue open. It runs under `cargo test`,
which is the project's primary gate, and it can be mutation-tested precisely. The deciding argument
is that a shell assertion in `ci.yml` would have been a fourth place to write shell that cannot fail,
which is the defect Open-6 found twice on the same day.

**5a asserts *direct* `needs`, which is stronger than the plan's wording.** §5 says the other jobs
"transitively" need `check-release`. They do not — all eight name it directly, so the test asserts
that, and a job that stops naming it fails even if it keeps an indirect path.

### Mutation evidence

Every assertion was mutated before being trusted, each anchor's match count asserted at 1 first
(§5: a mutation whose anchor is not unique silently does not apply):

| Mutation | Caught by |
| --- | --- |
| the notes step **relocated into `create-tag-and-release`**, present exactly once | `ga_notes_check_lives_in_check_release` |
| the notes step deleted outright | same test's first assertion |
| `if: always()` added to `publish-crate` | `release_workflow_has_no_status_functions` |
| `needs: [check-release]` removed from `build-binaries` | `every_other_job_needs_check_release` |
| the `jobs:` anchor broken | `job_parsing_sees_the_real_jobs`, plus two others, which is the point of having it |

The relocation case is the one that matters, since deleting the step is not the hazard — a check in
job 4 that *looks* present is. `release.yml` was restored byte-identically after each.

### Verification

912 tests passing, 0 failed (908 before, plus 5a's four). `cargo fmt --check` clean, clippy clean at
`-D warnings`, `actionlint` clean. CI green on `ccae784` across all eleven jobs before these changes.

### Not covered

Step 3 still has no CI coverage and cannot acquire any — see §5 Step 3. §8.2 stays unverified: that
`gh release create --notes-file` rejects a missing file is only testable by making a real release,
and it is a weak fallback regardless, firing after the tag push.

## 15. Code review, 2026-11-18

Dual code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) plus a Mode B coverage audit
(`COVERAGE.md`), all three in fresh contexts against commit `8f984ac`.

**Both reviewers cleared the mechanism.** A traced five dispatch shapes, B traced nine; neither found
a path this commit introduces to "tag pushed, no release". Both confirmed job-0 placement, the eight
direct `needs: check-release` edges and the absence of status functions. All five mutations §14
claims do genuinely fail. B credits one benefit the commit did not claim: `--generate-notes`
recomputes "changes since the previous tag" at run time, so a committed file makes a re-dispatch
deterministic.

### Taken: the fifth policy component

All three reports converged on it, and A put it best — the missing acknowledgements in `v2.4.0.md`
are not a backfill slip but a **policy gap every future GA inherits**, because `--generate-notes`
appended `## New Contributors` automatically and `--notes-file` does not. §3 now has a fifth
component; `.github/release-notes/README.md` states it; `v2.4.0.md` carries the three lines
byte-identically to the published body, so no wording anyone was credited with changes.

**The list cannot be regenerated later, which is the part worth keeping.** Asking the
`releases/generate-notes` API for v2.4.0 in November 2026 returns **two** of the three contributors
its published body credits — with or without an explicit `previous_tag_name` — because someone who
was a first-time contributor then is not one now. A policy that said "regenerate it" would have
dropped a contributor, which is the failure this component exists to prevent. So the README says to
capture the list at release time and treats the published body as the only authoritative record for
past releases.

Also worth noting against §5 Step 4's framing: deferring `gh release edit` was about
*authorisation*, not content. The file was already in the state it would have been applied in, and
§14 did not record the acknowledgements as outstanding — which is exactly how this would have been
lost, since the remaining step was described as needing only go-ahead.

### Findings that changed the record

**§8.2 assumption 2 is now verified, and its framing corrected.** `gh` 2.92.0 rejects a *missing*
`--notes-file` path at the file read, before any API call — no real release needed to establish it.
But `--notes-file ""` passes argument parsing and reaches the API, producing a release with an empty
body and exit 0. For the empty-path case the assumption's "a louder way to fail" is wrong: it is
silent. What prevents it is the notes step's `if:` being the exact complement of `Create release`'s
branch, which `the_notes_path_reaches_create_release` now pins.

§8.5's premise (no `CONTRIBUTING.md`, no docs page) is intact.

### Process note

Both reviewers were given the same working tree, so A found `release.yml` mid-mutation from B's
harness. A worked from `git show 8f984ac:...` in `$TMPDIR` and touched nothing; B restored
byte-identically and its sha256 matched. Concurrent reviewers that mutate a shared tree want
worktrees or a pinned copy.

### Taken, 2026-11-18 — the preamble and two comments

`release.yml:3-33`'s "Release invariant" block is the release checklist in practice, and its
canonical sequence would have walked into the job-0 refusal unwarned. Step 1 now names the notes
file, and a short paragraph states the refusal and points at the policy.

Both comments the review flagged are replaced. `# Interpolation only: a shell variable here would
expand to empty under bash -e` was wrong as written — a *defined* variable expands fine, the hazard
was an *undefined* one — and it contradicted `Detect release` twenty lines above, which is built
entirely from shell variables. Both also restated commit-message paragraphs in source, against
CLAUDE.md, so each is now one line stating the constraint it guards.

Comment-only: `actionlint` clean, YAML parses, the four tests pass. **Test hardening is still
outstanding** — it touches a test file, so it needs an explicit implementation trigger.

### Taken, 2026-11-18 — test hardening and the docs items

**The suite now catches all five escapes**, and one more I found while fixing them. `is_job_header`
no longer filters job names at all — any two-space key under `jobs:` that is not a comment is a job,
because GitHub allows `-`, `_` and digits in a job id and a narrower filter silently merges such a
job into the block above it. `job_parsing_sees_every_job` asserts the **exact nine-job list** rather
than two memberships, so a new job fails the test until someone adds it and checks its `needs:`.

`release_workflow_cannot_run_past_a_failed_job_zero` replaces the status-function test and covers
`continue-on-error` at both job and step level. It matches calls by name with any spacing inside the
parentheses, so `always ()` — which GitHub honours and actionlint accepts — is caught.

`the_notes_path_reaches_create_release` is new, and exists because of §8.2's corrected framing: an
empty `--notes-file` is not refused, so the chain from the check to the flag is pinned at every link,
together with the gate keeping the check on exactly the dispatches that take the `--notes-file` arm.

`needs_of` reads block-form `needs:` as well as inline, so that form no longer produces a false
failure.

**Thirteen mutations, each restored byte-identically.** Twelve behaved correctly on the first run.
The thirteenth did not, and it is worth recording: renaming `id: notes` to `id: notes_check` passed,
because `id: notes` is a *prefix* of `id: notes_check` and the assertion used `contains`. That is the
same class of defect as everything else in §15 — an assertion that cannot fail — reintroduced inside
the fix for it. The fragment is now line-exact, and the rename is caught by the test as well as by
actionlint.

**Docs.** The README now points at `v2.4.0.md` as the worked example Step 1 asked for, states what
the CI test actually asserts (four properties, not two), and links highlights as `/issues/NNN` to
match the backfill and `CHANGELOG.md`. `v2.4.0.md` links the CHANGELOG section anchor
(`#version-240-release-on-9-november-2026`, confirmed present in the built site) rather than the page
root.

**Verification.** 913 tests passing, 0 failed. `cargo fmt --check` clean, clippy clean at
`-D warnings`, `actionlint` clean.

Everything §15 reported is now either taken or recorded as permanently uncovered (§5 Step 3).
