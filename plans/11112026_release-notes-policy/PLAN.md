# PLAN — curated release notes, and a workflow that cannot ship without them

**Issue:** none filed (policy decision taken in session, 2026-11-11)
**Branch base:** `dev` @ `b0a83b0` (= `origin/master` = the `v2.5.0` tag; was `ee675c1` when drafted)
**Date:** 2026-11-11, manual review 2026-11-13
**Decision taken:** curated highlights only — drop `--generate-notes` in favour of a committed
notes file. Chosen over "highlights + folded full list" and "filter the auto list".
**Status:** **Revision 2** — dual agent review complete (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`),
findings folded in; see §11. Awaiting re-review or an implementation trigger. No implementation.

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
2. **Three to six highlights**, one sentence each, in the CHANGELOG's own voice — `- **Headline**
   ([#NNN](url)) — one sentence.` Selected as "a user would change what they do because of this".
   Dependency bumps, CI, docs and internal refactors are never highlights.
3. **A link to the CHANGELOG section** for the full record.
4. **An install block** — cargo, conda and docker lines. **Note these do not exist to copy.**
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

The policy from §3 ("what a GA body contains"), the four-part shape, and one worked example. Short:
this is a rule, not an essay.

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
new shell that reaches production unexercised. Mitigations, all cheap: the interpolation-only form
above has nothing to be undefined; `actionlint` in CI would catch the class outright (there is no
linter in `.github/` today — see §10 Open-6); and a GA is always dispatched by hand, so the first
real use is watched.

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
