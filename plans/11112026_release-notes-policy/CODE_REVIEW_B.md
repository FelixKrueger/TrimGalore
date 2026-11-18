# CODE REVIEW B — `8f984ac` (`ci(release): require curated notes for a GA release`)

**Reviewer:** B (independent; no shared state with Reviewer A)
**Date:** 2026-11-18
**Target:** commit `8f984ac` on `dev` — `.github/workflows/release.yml`,
`.github/release-notes/README.md`, `.github/release-notes/v2.4.0.md`,
`tests/integration_release_workflow.rs`
**Plan:** `plans/11112026_release-notes-policy/PLAN.md` (read in full; §3, §5, §8, §14)
**Tree state:** `release.yml` was mutated nine times during this review and restored
byte-identically after each; final `shasum -a 256` is
`c9fc18d7a4885f626a8c7a6d700401c0c7326cda21c509a387607a667a926533`, `git diff` empty.
No other file was touched. Local `actionlint` is **1.7.12**, the exact version
`ci.yml:172` pins, so every actionlint result below transfers to CI.

---

## Summary

**The workflow behaviour is correct.** I traced all seven reachable dispatch cases and
found none that publishes a wrong body, publishes nothing after tagging, or fails
undiagnosably. The job-0 placement argument holds, and `--notes-file` is a genuine
improvement on `--generate-notes` for re-dispatch determinism (noted below, the commit
does not claim it).

**The tests are weaker than the commit message says.** The commit's central claim is that
`tests/integration_release_workflow.rs` "asserts the properties that make job-0 placement
sufficient", and that the relocation mutation "is the hazard rather than deleting it". I
reproduced that mutation — it is real and it is caught. But the same relocation, into a job
whose name contains an underscore or a digit, **passes all four tests and `actionlint`**. I
proved this with a complete, coherent refactor, not a contrivance. Two of the three
substantive assertions rest on one 7-line name predicate, and `job_parsing_sees_the_real_jobs`
— the test written to stop exactly that — checks two hardcoded names and does not notice that
a job went missing from the parse.

One narrow path to "tag pushed, no release" remains, and it is on the dangerous side of a
change §13 itself recommends making next.

Verdict: **ship the workflow; fix the test before anyone relies on it.** Findings 1–3 are one
afternoon's work and two of them are a one-line predicate swap.

---

## Part 1 — Control-flow trace (the question asked)

`check-release` is job 0 and all eight other jobs name it **directly** in `needs:`
(verified by reading, not by the test). A GitHub `if:` that contains no status function has
`success()` applied implicitly, at job level as well as step level, so a failed job 0 skips
every dependent job. The workflow contains no status function (confirmed: `grep` finds
none), so nothing overrides that. The placement argument is sound.

| # | Dispatch | `check-release` outputs | What runs | `Create release` invokes | Verdict |
|---|---|---|---|---|---|
| 1 | GA from `master`, notes file present | `is_release=true`, `is_prerelease=false`, `notes_path=.github/release-notes/vX.Y.Z.md` | all jobs | `--notes-file "<path>"` | correct |
| 2 | GA from `master`, notes file **absent** | job 0 **fails** at `Require curated GA release notes` | nothing — all eight jobs skipped | — | correct; error names the path and the README |
| 3 | GA from `master`, notes file whitespace-only | job 0 **fails** (same step) | nothing | — | correct — `grep -q '[^[:space:]]'` catches it where `-s` would not |
| 4 | Prerelease from `dev` | `is_prerelease=true`; notes step skipped; `notes_path` empty | all jobs | `--generate-notes --prerelease` | correct — empty `notes_path` is never read |
| 5 | Dry-run GA | notes check **runs** and gates; `Create tag`/`Create release` skipped | jobs 1–3 run, job 4 **skipped** (`smoke-test-docker` is dry-run-gated and job 4 needs it) | — | correct, and matches §5's corrected expectation |
| 6 | Dry-run prerelease | notes step skipped | as above | — | correct |
| 7 | Re-dispatch after partial ship (tag exists at HEAD) | `tag_preexists=true`; notes check passes (file is committed) | `Create tag` skipped, `Reuse existing tag` logs | `--notes-file` — deterministic | **improvement**, see below |
| 8 | GA version dispatched from `dev` | `Detect release` hard-errors at the branch gate | nothing | — | correct (pre-existing) |
| 9 | `push` to `dev`/`master` | `is_release=false`; notes step skipped | docker jobs only (pre-existing) | — | correct |

### Can a tag still be pushed with no release?

Yes, by four routes. Three are pre-existing and untouched by this commit: `gh release create`
failing on its own (API error, or the release already existing from an earlier partial run),
the runner dying between the two steps, and a manual cancellation inside the window §13
Open-7 describes. The fourth is new and is **Finding 3**.

### An improvement the commit does not claim

`--generate-notes` computes "changes since the previous tag" **at the moment it runs**. On a
re-dispatch after a partial ship, or after any intervening tag, that computation can produce
a different body than the first attempt would have. `--notes-file` reads a committed file, so
case 7 now reproduces the same body byte-for-byte. The README's "the file is the source for
*creation* only" covers the policy; the determinism is worth knowing.

### Checked and clean

- `Create release`'s `else` branch is reachable only when `is_prerelease == 'false'`, which is
  exactly the condition under which the notes step ran — the two gates cannot disagree.
- Both `actions/checkout` calls (`:73`, `:403`) omit `ref:`, so job 0 and job 4 read the same
  commit. §8.6 verified this; I re-read it rather than taking it.
- `grep -q '[^[:space:]]'` returns non-zero for missing (2), empty (1) and whitespace-only (1),
  and the `!` is consumed inside an `if` condition, so it is not the errexit-exempt form
  §13 Open-6 found twice. §14's deviation note is accurate.
- The GA dispatch fails **before** any build starts, so a missing notes file costs seconds, not
  the ~30-minute build.

---

## Part 2 — Findings

### HIGH 1 — The relocation mutation the commit relies on can be made to pass

`tests/integration_release_workflow.rs:33` keys job headers on
`chars().all(|c| c.is_ascii_lowercase() || c == '-')`. Any job name containing an underscore
or a digit fails that predicate, and its lines are **absorbed into the preceding job's block**.

I moved the GA notes check out of `check-release` into a new job `gate_notes` placed directly
after it, gave `gate_notes` the `notes_path` output, and pointed job 4's `--notes-file` at
`needs.gate_notes.outputs.notes_path`. A complete refactor, nothing left dangling:

```
actionlint exit=0
test every_other_job_needs_check_release ... ok
test job_parsing_sees_the_real_jobs ... ok
test ga_notes_check_lives_in_check_release ... ok
test release_workflow_has_no_status_functions ... ok
test result: ok. 4 passed; 0 failed
```

Both CI gates green with the check outside `check-release`. For comparison, the same
relocation into `create-tag-and-release` (hyphen-named) **does** fail
`ga_notes_check_lives_in_check_release` — the commit's claimed mutation is real. The test
works for hyphen-named jobs and is blind to every other name, and this repository already
uses underscores and digits in adjacent identifiers (`linux-x86_64`, `ubuntu-24.04-arm`,
`binary-linux-x86_64`).

`gate_notes` as I built it publishes nothing, so that specific arrangement is harmless. The
finding is that the assertion does not hold the property it is written to hold.

### HIGH 2 — A job with a digit or underscore in its name is never audited for `needs:`

Same root cause, sharper consequence. I appended a job with **no `needs:` at all**:

```yaml
  publish-crate-2:
    name: Publish mirror
    runs-on: ubuntu-latest
    steps:
      - run: echo "no needs at all, publishes things"
```

`actionlint` exit 0 — a job without `needs` is legal — and all four tests green. It was
absorbed into `publish-crate`'s block, where `find(… starts_with("needs:"))` returned
`publish-crate`'s own `needs:` line, which contains `check-release`. So the test asserted the
wrong job's dependency and reported success.

`every_other_job_needs_check_release` exists to make "a failed job 0 stops everything" a
checked fact. For any non-conforming job name it is vacuous — the shape the project's standing
lesson is about.

**Fix for 1 and 2 (one predicate).** Drop the character filter and exclude comment lines
instead:

```rust
let is_job_header = |l: &str| {
    l.starts_with("  ")
        && !l.starts_with("   ")
        && l.trim_end().ends_with(':')
        && !l.trim_start().starts_with('#')
};
```

I simulated both predicates against the current file: they return the identical nine-job list,
so the swap is a no-op today. Anchoring at `jobs:` already excludes `on:`'s `push:`, which is
what the character filter was defending against — and `job_parsing_sees_the_real_jobs` keeps
asserting that. A quoted key (`"publish-crate":`) still slips through; `trim_matches('"')`
closes it if you want the belt.

### HIGH 3 — One uncovered mutation still reaches `--notes-file ""` after the tag push

Step 2a's single-path argument is the right design, and I tested all three ways to break the
link:

| Mutation | `actionlint` | tests |
|---|---|---|
| `notes_path:` removed from `check-release.outputs` | **catches it** — `property "notes_path" is not defined in object type {…}` | pass |
| step `id: notes` renamed, `outputs:` unchanged | **catches it** — `property "notes" is not defined` | pass |
| `echo "notes_path=$NOTES" >> "$GITHUB_OUTPUT"` deleted from the `run:` block | exit 0 | pass |

So the linkage is covered — by `actionlint`, not by the new test — except for the one case
where the step stops writing the output. `actionlint` cannot type what a `run:` block writes,
and nothing else looks. The consequence is a GA run reaching
`gh release create … --notes-file ""`, which fails when `gh` tries to open the empty path —
**after** `git push origin` in the same job. The v2.3.0 shape, from a one-line edit.

This is not hypothetical maintenance. §13 Open-5's own follow-up recommends lifting the check
"into a script with its own test". A script that validates the file and forgets to emit the
output lands exactly this.

**Fix.** Assert the three halves of the link in the test — that `check-release`'s `outputs:`
block declares `notes_path`, that the notes step's `run:` writes `notes_path=` to
`$GITHUB_OUTPUT`, and that `Create release` reads `needs.check-release.outputs.notes_path`.
Three `contains` assertions on text already in hand; the first two are what `actionlint`
cannot see.

### MEDIUM 4 — `continue-on-error` is not in the banned-token list

`STATUS_FUNCTIONS` bans `always()`, `cancelled()`, `success()` and `failure()`. It does not
ban `continue-on-error`, which defeats the same property more directly: `continue-on-error:
true` on the notes step makes the step's `exit 1` not fail the job, the `echo notes_path=`
line never runs, and the GA proceeds to `--notes-file ""` after the tag push — Finding 3's
outcome, reached by a flag whose whole purpose is "do not let this failure stop the run". At
job level it makes `check-release` report success to every dependent job.

There is no `continue-on-error` in any workflow today, so adding the token costs nothing. Rename
the const (`FORBIDDEN_TOKENS`) and add it.

### MEDIUM 5 — Block-form `needs:` produces a false failure with a false message

```yaml
    needs:
      - check-release
      - build-binaries
```

```
job `create-tag-and-release` does not need `check-release`: needs:
```

The job *does* need `check-release`. The parser takes the `needs:` line and looks for the
dependency on it; in block form the dependency is on the following lines. A YAML reformat — or
anyone who finds a five-element inline list long — gets a red test asserting something untrue,
and the obvious response to an assertion that is wrong is to weaken it. Fails safe, but the
message should not lie: collect the `needs:` line plus any following `      - ` lines, or
assert the block form is absent with a message that says so.

### MEDIUM 6 — The `Create release` comment is inaccurate and contradicts the step above it

```yaml
      # Interpolation only: a shell variable here would expand to empty under
      # `bash -e` and fail `gh release create` after the tag push.
```

A *defined* shell variable expands to its value, not to empty. The hazard §11 records was an
**undefined** variable — `$VERSION`/`$SHA`/`$TITLE`, which nothing assigned. As written the
comment states a general rule about shell variables that is false, and that rule is
contradicted twenty lines above by `Detect release`, which is built entirely out of
`RAW_VERSION`, `IS_PRERELEASE`, `VERSION`, `MAJOR`, `MINOR`, `TAG_PREEXISTS` and
`LOCK_VERSION`. A maintainer who believes the comment will not know why one step may use shell
variables and the other may not.

It is also change-justification in source — "if you find yourself pre-empting *but doesn't this
break X?*, that argument belongs in the commit message" — and the commit message already
carries it verbatim ("PRERELEASE_FLAG is gone: with an explicit branch it was only ever a way
to avoid one"). Per CLAUDE.md, "if the comment restates its own commit message, delete the
comment". My recommendation is to delete it.

The sibling comment on the notes step (`# In job 0, which publishes nothing. A check in
create-tag-and-release would fail after the tag push and leave a tag with no release.`) is
the defensible case — it is an invariant written on the thing it constrains, and it is why
someone must not move the step. It is two lines, which is the stated maximum. Keep it; consider
trimming the second sentence to "would strand the tag".

### MEDIUM 7 — The backfill drops the contributor acknowledgements §5 Step 4 said must be carried

Step 4: *"v2.4.0's auto-generated body carries three first-time-contributor acknowledgements,
which `--notes-file` would delete. Carry them into the curated body, or leave v2.4.0 alone —
the policy does not need a backfilled example badly enough to drop someone's first
contribution off a release page."*

`.github/release-notes/v2.4.0.md` carries none, and §14's Deviations do not record the choice.
Nothing is lost yet — §9.1 reserves the `gh release edit` for go-ahead, correctly — but the
file as written is the thing that would be applied.

The policy has the same hole going forward, which matters more. `--generate-notes` appends a
`## New Contributors` section automatically; GA bodies no longer get one, and
`.github/release-notes/README.md`'s "In order, nothing else" does not mention contributors at
all. Every future GA silently drops first-contribution credit unless the author remembers.
Add a fifth element to the README's shape, or state explicitly that credit lives in the
CHANGELOG and link it.

### LOW 8 — The status-function ban blocks a legitimate pattern, with a wrong explanation

A release-failure notification job (`if: failure()`) is an ordinary thing to add and would
publish nothing. The test fails it with "which can let a job run after `check-release`
failed" — true, and in that case intended. A maintainer whose correct addition is rejected by
an assertion that explains itself wrongly deletes the assertion. Either scope the ban to the
publishing jobs, or name the escape hatch in the test's message.

### LOW 9 — `job_parsing_sees_the_real_jobs` does not notice a job going missing

It is described as the guard that stops the other two being vacuous, and it checks two
hardcoded names plus the absence of `push`. Both are true in every mutation above. The
assertion that would have caught Findings 1 and 2 is that the parse accounts for every
two-space key under `jobs:` — which the permissive predicate makes automatic, or which
`assert_eq!(names.len(), 9)` makes explicit at the cost of a test edit per new job. The
permissive predicate is the better buy.

### LOW 10 — Smaller things

- **Test module doc.** Ten lines of narrative in committed source, including the v2.3.0
  incident and the reasoning chain. CLAUDE.md's one-line default covers tests explicitly, and
  the commit message already carries all of it. Two lines would do: what the file asserts, and
  that the check must stay in the job that publishes nothing.
- **`2>/dev/null` is an undocumented deviation.** The plan's snippet has no redirect; §14's
  Deviations covers dropping `[ ! -s ]` but not adding this. It also means a file that exists
  but is unreadable reports "GA releases require curated notes at <path>", which sends the
  maintainer to create a file that is already there. Niche, but one `[ -r "$NOTES" ]` ahead of
  the `grep` would separate the two.
- **The record already has a hole.** `.github/release-notes/` holds `v2.4.0.md` and nothing
  else, while `v2.5.0` is tagged and shipped (Cargo.toml is at `2.5.0`). The README calls these
  files a record that accumulates. Either backfill `v2.5.0.md` or say the record starts at the
  first release the policy gated.
- **Shell interpolation of `version`.** `NOTES=".github/release-notes/${{ steps.check.outputs.version }}.md"`
  interpolates a Cargo.toml-derived value into a shell script. The source is a protected branch
  and job 4 already does the same in `git tag "${{ … }}"`, so this adds no exposure class —
  noting it only because `env:`-passing is the form that would not need the argument.

---

## Part 3 — What I verified rather than assumed

- All four tests pass on the committed tree; `actionlint` (1.7.12, CI's pin) exits 0 on
  `release.yml`.
- `ci.yml:164-184` lints `release.yml` on every non-schedule event and asserts the glob is
  non-empty; `ci.yml:76` runs `cargo test --release`. Both gates are real, so the new test and
  the actionlint backstop both actually run.
- Nine mutations, each restored byte-identically; `release.yml`'s sha256 and `git diff` are
  clean at the end of this review.
- All eight non-`check-release` jobs name `check-release` directly in `needs:` — §14's
  "stronger than the plan's wording" claim is accurate.
- No `continue-on-error` and no status function exists in any workflow today.
- The strict and permissive header predicates return the identical nine-job list on the
  current file.

## Not verified

§8.2 — that `gh release create --notes-file` rejects a missing or empty path with a non-zero
exit. Finding 3 depends on it producing a *failure* rather than a silent empty body; the plan
is right that the only honest test creates a real release. Either way the step is past the tag
push, so neither outcome is good: a failure strands the tag, a silent pass ships a blank GA
body. That is the reason Finding 3 matters, not a reason to discount it.
