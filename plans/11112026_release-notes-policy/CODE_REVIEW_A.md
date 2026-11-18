# CODE REVIEW A — `8f984ac` *ci(release): require curated notes for a GA release*

**Reviewer:** A (independent)
**Date:** 2026-11-18
**Target:** commit `8f984ac` on `dev` — `.github/workflows/release.yml`, `.github/release-notes/README.md`, `.github/release-notes/v2.4.0.md`, `tests/integration_release_workflow.rs`
**Plan:** `plans/11112026_release-notes-policy/PLAN.md` (revision 3)

**Method note.** `release.yml` was being mutated in the working tree by another reviewer while this
review ran (`git diff HEAD` showed one deleted line). Nothing here touched the shared tree. All
mutation work ran against `git show 8f984ac:.github/workflows/release.yml` copied into
`$TMPDIR/relrev_A/`, driven by a standalone Rust harness (`harness.rs`) that embeds the test's
`jobs()` parser and the four assertion bodies verbatim. The harness reproduces the real suite exactly
on the pristine file: 4 pass, 9 jobs parsed.

---

## Summary

The central question — *can the workflow still reach a state where a tag is pushed but no release
exists?* — is **no, not through anything this commit introduces**. The job-0 placement is correct,
the reasoning behind it holds, and the `Create release` rewrite is sound. Job 0 and job 4 demonstrably
read the same commit (every `actions/checkout` omits `ref:`, so all six pin `github.sha`), so the file
proven present in job 0 is present in job 4.

The shell is correct on every arm it must handle. The `Create release` branch condition is the exact
complement of the notes step's `if:` across all five reachable dispatch shapes, so the GA branch
cannot currently receive an empty `notes_path`.

Two things are wrong, neither of them in the workflow:

1. **`v2.4.0.md` drops three first-time-contributor acknowledgements** that the plan itself said were
   not worth dropping, and the policy in `README.md` has no slot for them at all — so this is not a
   backfill slip, it is a gap in the policy that every future GA body inherits.
2. **The test suite is weaker than its own module doc claims.** Five mutations pass all four tests,
   including the exact relocation hazard the suite exists to catch, reached by giving a job a name
   with an underscore in it.

One correction to the plan: **§8.2's assumption 2 is now verified** (no real release needed), and its
framing is half wrong — `gh` fails fast on a *missing* file, but treats an *empty path* as if the flag
were absent and creates a bodiless release with exit 0.

---

## Verified correct — answers to the five questions

### Q1. The GA notes check's shell

```bash
NOTES=".github/release-notes/${{ steps.check.outputs.version }}.md"
if ! grep -q '[^[:space:]]' "$NOTES" 2>/dev/null; then
  echo "::error::GA releases require curated notes at $NOTES (…)"
  exit 1
fi
echo "notes_path=$NOTES" >> "$GITHUB_OUTPUT"
```

Simulated under `bash -e` with the real script body, four arms, `GITHUB_OUTPUT` captured:

| Arm | exit | `GITHUB_OUTPUT` |
|---|---|---|
| file missing | 1 | empty |
| file empty (0 bytes) | 1 | empty |
| whitespace-only (`  \n\t\n`) | 1 | empty |
| real content | 0 | `notes_path=.github/release-notes/v9.2.0.md` |

It refuses every case it must and accepts a valid file, and writes the output *only* on the pass arm.
Underlying `grep` exits: missing → 2, empty → 1, whitespace → 1, bare newline → 1, directory →
non-zero, content → 0.

**The `!` is not inert here.** SC2251's errexit exemption applies to `! cmd` as a *simple command*;
inside an `if` condition the `!` is consumed by the conditional. §14's deviation note is right, and the
one-`grep` form is genuinely equivalent to the plan's two-test form on all four arms.

**`2>/dev/null` hides nothing load-bearing**, but it does make the message misleading for the failure
modes that are *not* "file absent" — permission denied, an I/O error, or the path being a directory all
print "GA releases require curated notes at $NOTES", which reads as "write the file" when the file is
already there. Minor; see L5.

### Q2. Which dispatches run the step

The step's gate is `is_release == 'true' && is_prerelease == 'false'`. Enumerated:

| Dispatch | notes step | `notes_path` | `Create release` |
|---|---|---|---|
| `push` to dev/master | skipped (`is_release=false`) | empty | job skipped (`is_release != 'true'`) |
| dispatch, GA, `dry_run: false` | **runs** | set | GA branch, consumes it |
| dispatch, GA, `dry_run: true` | **runs** (validates) | set | step skipped (`is_dry_run`) |
| dispatch, prerelease (dev *or* master), either dry-run | skipped | empty | prerelease branch, does not read it |
| dispatch, dev with a GA version / any other branch | n/a | n/a | `Detect release` exits 1 at `:91`, job 0 fails, all skipped |

**The GA branch cannot receive an empty `notes_path`.** `is_prerelease` is written on every path
through `Detect release` and is always exactly `"true"` or `"false"`, never empty; `Create release`
only runs when job 0 succeeded *and* `is_release == 'true'`. So `is_prerelease == 'false'` and
`is_prerelease != 'true'` select the same set, and the step's gate is the branch condition's exact
complement. Correct — but load-bearing and untested; see **M3**.

### Q3. The `Create release` rewrite

No invocation can produce a wrong or empty body in the current file, and nothing in the rewrite can
fail after the tag push that could not before:

- Both branches parse (`bash -n` on the substituted form) and are pure interpolation — no shell
  variable to expand to empty, which was the revision-1 defect.
- `--notes-file` resolves against `$GITHUB_WORKSPACE` (no `working-directory`), and `notes_path` is
  repo-relative. Correct.
- `--prerelease` now appears only on the prerelease branch; `PRERELEASE_FLAG`'s removal is behaviour-
  preserving.
- Residual tag-push-without-release exposures — a `gh` API failure, and a re-dispatch against an
  already-existing release — are pre-existing and unchanged by this commit. Open-7's cancellation
  window is already remedied in the file (`concurrency` keyed on `github.event_name`).

`actionlint` is clean on the committed workflow (exit 0, no diagnostics, with
`SHELLCHECK_OPTS='-e SC2086 -e SC2012 -e SC2129'`).

### Q4. The tests — see M1, M2, M3 and the Low items. Five mutations escape.

### Q5. The markdown files — see H1, H2, L1, L2, L6.

---

## High

### H1 — `v2.4.0.md` drops three first-time-contributor acknowledgements

Plan Step 4 states the cost and names the two acceptable resolutions:

> v2.4.0's auto-generated body carries three first-time-contributor acknowledgements, which
> `--notes-file` would delete. Carry them into the curated body, **or leave v2.4.0 alone** — the policy
> does not need a backfilled example badly enough to drop someone's first contribution off a release
> page.

The file as written takes neither. Verified against the live release body
(`api.github.com/repos/FelixKrueger/TrimGalore/releases/tags/v2.4.0`), which carries:

```
## New Contributors
* @BenjaminDEMAILLE made their first contribution in …/pull/336
* @life4sucks made their first contribution in …/pull/347
* @sujeito-operator made their first contribution in …/pull/436
```

`.github/release-notes/v2.4.0.md` contains none of them. Nothing is lost *yet* — §14 correctly
reserves `gh release edit` for go-ahead — but deferring the publish was about authorisation, not about
content, so the file is in the state it will be used in. Running `--notes-file` with it as it stands
executes the outcome the plan ruled out.

**Fix:** append the three lines to `v2.4.0.md` before any `gh release edit`, and note in `README.md`
that a backfill carries forward whatever acknowledgements the auto-body held.

### H2 — the policy has no slot for contributor acknowledgements

`README.md`'s "What a body contains" is explicitly closed: *"In order, nothing else"*, four items
(framing, highlights, CHANGELOG link, install block). None of them is contributor credit.

H1 is therefore not a one-off. Every GA body produced under this policy silently drops the
`## New Contributors` section that `--generate-notes` supplied for free, on the one page where GitHub
surfaces a first contribution. The plan judged this cost serious enough to consider abandoning the
backfill over it; the policy should carry the same judgement.

**Fix:** add a fifth item to `README.md` — a `## New Contributors` section, carried over from the
auto-notes for that range, omitted when there are none. It is mechanical to produce
(`gh api …/releases/generate-notes` or the compare view), so it costs the author nothing and is
exactly the kind of content the one-sentence rule was never aimed at.

---

## Medium

### M1 — the job parser silently ignores job IDs with `_` or a digit, and two tests go vacuous for them

`is_job_header` requires every character of the name to satisfy
`c.is_ascii_lowercase() || c == '-'`. GitHub job IDs legally allow `_` and digits. A job whose ID uses
either is **not recognised as a header**, so its lines are absorbed into the preceding job's block and
it disappears from `jobs()` entirely.

Verified — each of these parses to exactly 9 jobs, the real nine, with the added job invisible:

| Mutant | What | Result |
|---|---|---|
| **M8** | notes step relocated into `publish_notes:` placed after `check-release`, present exactly once | **all 4 tests pass** |
| **M6** | extra job `verify_notes:` with no `needs:` at all | **all 4 tests pass** |
| **M7** | extra job `docker-build2:` with no `needs:` at all | **all 4 tests pass** |

M8 is the one that matters: it is §14's "relocation" mutation — *"the one that matters, since deleting
the step is not the hazard — a check in job 4 that looks present is"* — reached by nothing more than
naming the destination job with an underscore. The suite's headline guarantee has a one-character
bypass.

There is a second-order effect in `every_other_job_needs_check_release`: it takes
`.find(|l| … starts_with("needs:"))` over the *merged* block, so an absorbed job inherits the
preceding job's `needs:` line and passes on it. The job is never tested and nothing says so.

**Fix**, both parts:

```rust
.all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
```

and make `job_parsing_sees_the_real_jobs` assert the *whole* set rather than two memberships, so an
unparsed job fails loudly instead of vanishing:

```rust
assert_eq!(
    names,
    ["check-release", "build-binaries", "smoke-test-binaries", "docker-build",
     "docker-merge", "smoke-test-docker", "create-tag-and-release",
     "upload-binaries", "publish-crate"],
    "job parsing disagrees with release.yml"
);
```

That also turns adding a job into a deliberate test update, which is the right friction for this file.

### M2 — `continue-on-error` defeats the stated invariant and nothing tests for it

The module doc asserts sufficiency:

> Nothing can run after a failed job 0 as long as two things hold: every other job needs
> `check-release`, and the workflow has no status functions.

That is not sufficient. `continue-on-error: true` is a third lever, at two levels:

| Mutant | What | Result |
|---|---|---|
| **M9** | `continue-on-error: true` on the `check-release` **job** — dependents treat it as successful | all 4 tests pass |
| **M10** | `continue-on-error: true` on the notes **step** — the check cannot fail the job | all 4 tests pass |

**Fix:** add `"continue-on-error"` to the forbidden-token scan. It appears nowhere in `release.yml`
today, so the assertion starts green and can genuinely fail. Also soften the module doc from "as long
as two things hold" to naming the three it actually checks.

### M3 — an empty `notes_path` fails *silently*, not loudly, and the coupling that prevents it is untested

Probed `gh` 2.92.0 directly (scratch config dir, unresolvable host, nothing real contacted):

| Invocation | Behaviour |
|---|---|
| `--notes-file /definitely/missing/notes.md` | `open …: no such file or directory` — **fails before any API call** |
| `--notes-file ""` | goes straight to `POST …/releases` — **identical to passing no notes flag** |
| no notes flag at all | goes straight to `POST …/releases` |

Two consequences.

**Plan §8.2 can be resolved, and its framing corrected.** Assumption 2 — *"`gh release create
--notes-file` rejects a missing file with a non-zero exit"* — is **verified true**, and it did not take
a real release: the file read precedes the API call. But §8.2 calls it "a louder way to fail", and that
only holds for a missing file at a *non-empty* path. For an *empty* path there is no failure at all:
gh reached the POST without prompting even though stdout was not a TTY, so the result is a GA release
with an empty body and exit 0 — no signal anywhere.

**So the only thing standing between a mistake and a bodiless GA release is the structural coupling in
Q2** — the notes step's `if:` being the exact complement of the `Create release` branch. That coupling
is correct today and completely untested:

| Mutant | What | Result |
|---|---|---|
| **M12** | one token in the step's `if:` flipped (`is_prerelease == 'true'`), so the check never fires on a GA | all 4 tests pass |
| **M15** | step name kept, shell body gutted so it validates nothing | all 4 tests pass |

**Fix — in the test, not in job 4.** Assert that the notes step's `if:` and the `Create release` GA
branch agree, e.g. that the notes step's condition contains `is_prerelease == 'false'` *and* that the
`Create release` body branches on `is_prerelease`. A guard *inside* `Create release` is tempting and is
the wrong move: it would fire after the tag push, trading a silently bad body for a stranded tag, which
is strictly worse. Fixing it in job 0 is what the step already does; what is missing is the assertion
that it still does it.

M15 also marks the suite's honest scope limit: it proves the step's *placement*, never its *efficacy*.
Worth one line in the module doc, since §13 Open-5 deliberately dropped 5b and this is now the only
place a reader learns the absent arm is uncovered.

---

## Low

- **L1 — `README.md`'s highlight template contradicts its own worked example.** The template at
  `README.md:31` prescribes `…/pull/NNN`; `v2.4.0.md` uses `…/issues/NNN` for all six highlights, as
  does `CHANGELOG.md` (75 `/issues/` vs 22 `/pull/`). The template is the outlier — change it to
  `issues`.
- **L2 — `v2.4.0.md` links the changelog *page*, not the *section* `README.md` requires.** Per-version
  anchors exist and are confirmed in the built site
  (`docs/dist/reference/changelog/index.html` → `id="version-240-release-on-9-november-2026"`). Since
  this file is the worked example, the omission propagates to every future body.
- **L3 — the two new comments in `release.yml` restate the commit message and justify the change in
  source**, against CLAUDE.md's explicit rule (*"Never justify the change in the source"*, *"If the
  comment restates its own commit message, delete the comment"*). Both map near-verbatim onto commit-
  message paragraphs 2 and 4. Suggested: cut the `Create release` comment entirely (the commit message
  carries it), and reduce the notes-step comment to one line stating the invariant on the thing it
  constrains, e.g. `# Must stay in job 0: create-tag-and-release pushes the tag before creating the
  release.` The test file's 9-line module doc is the same instinct; defensible for a test, but it would
  read better at four lines.
- **L4 — `needs:` in YAML block form is a false positive.** Mutant **M11** rewrites
  `needs: [check-release]` as a block sequence — equivalent YAML — and
  `every_other_job_needs_check_release` fails with *"does not need `check-release`: needs:"*. It fails
  safely and loudly, so this is brittleness rather than a hole, but anyone reformatting the file will
  hit it. Worth either handling the block form or saying in the test that only the flow form is
  supported.
- **L5 — `2>/dev/null` makes the error message wrong for non-absent-file failures** (permission, I/O,
  path is a directory): it tells the reader to create a file that already exists. Dropping the
  redirect would let grep's own diagnostic through alongside the `::error::` line.
- **L6 — `README.md` undersells the CI test.** It says a test asserts two properties; there are three
  tests, and the third (`every_other_job_needs_check_release`) is half of what makes job-0 placement
  sufficient.
- **L7 — `NOTES_STEP`'s doc comment points at `plans/11112026_release-notes-policy/`**, which is in
  `Cargo.toml`'s `exclude` list and so is absent from the published crate. A dangling reference for any
  reader outside this checkout.
- **L8 — two theoretical evasions of the token scan**, both verified, neither likely to be typed by
  accident: `always ()` with a space is accepted by `actionlint` (exit 0, no diagnostics) and missed by
  `src.contains("always()")` (**M13b**); and `needs.contains("check-release")` is a substring match, so
  `needs: [check-release-lint]` passes (**M14**).

---

## Checked and dismissed — recorded so they are not re-derived

- **Job 0 and job 4 cannot disagree about the file.** All six `actions/checkout` invocations omit
  `ref:`, so each pins `github.sha`. §8.6 holds as written.
- **The one-`grep` deviation (§14) is sound**, and the `!` is not the errexit-exempt form — confirmed
  by execution, not reasoning.
- **The plan's five claimed mutations all genuinely fail.** M1 relocation → t1; M2 deletion → t1's
  first assertion; M3 `if: always()` → t2; M4 dropped `needs` → t3; M5 broken `jobs:` anchor → panic in
  all three dependent tests. §14's mutation table is honest.
- **`job_parsing_sees_the_real_jobs` is not itself vacuous** — `jobs:` is the last top-level key, so
  `on:`'s two-space `push:` cannot be reached, and the assertion can only fail if the anchoring breaks,
  which it is there to detect.
- **The suite is actually gated** — `ci.yml:76` runs `cargo test --release`.
- **`actionlint` clean** on the committed workflow.
- **The test depends on `.github/`, which `Cargo.toml` excludes from the published crate** — so
  `cargo test` from a crates.io tarball panics. Dismissed as pre-existing practice, not a regression:
  17 of the other integration tests already depend on `test_files/`, which is excluded the same way.
  `cargo publish --no-verify` is unaffected.
- **A BOM-only file would pass the content check under GNU grep** (the BOM bytes are non-space).
  Outside the specified arms and below the noise floor.
- **Version-string injection through `${{ … }}`** (a `"` in `Cargo.toml`'s version breaking out of the
  quoting) is pre-existing — the old step interpolated `version` too — and needs repo write plus a
  dispatch. Not introduced here.

---

## Recommended order

1. **H1** — append the three acknowledgements to `v2.4.0.md` before any `gh release edit`.
2. **H2** — give the policy a contributor slot, so H1 cannot recur.
3. **M1** — widen the job-name char class; assert the full job set.
4. **M2** — add `continue-on-error` to the token scan; correct the module doc's sufficiency claim.
5. **M3** — assert the `if:`/branch coupling; update §8.2 (assumption 2 verified, framing corrected).
6. **L1, L2, L6** — the three markdown contradictions, one line each.
7. **L3** — trim the two comments.
8. **L4, L5, L7, L8** — as convenient.

No fixes were applied. Nothing in the working tree was modified.
