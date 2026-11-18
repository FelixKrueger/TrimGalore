# RELEASE — v2.6.0

**Date drafted:** 2026-11-18
**Base:** `origin/master` = `b0a83b0` (the `v2.5.0` tag). `dev` = `96e7cb1`, **5 commits unpushed**.
**Status:** not started. **No release has been initiated and none should be without Felix asking.**

---

## 0. What is in this release

19 commits on `dev` not in `master`. By type: 5 docs, 5 ci, 4 build(deps), 2 test, 1 fix, 1 feat,
1 chore.

**User-facing content is two items:** `--fastqc_raw` (#576) and the `rustix` 1.1.6 bump clearing
RUSTSEC-2026-0486 (#575). Everything else is dependency bumps, CI hardening, tests, plan artifacts,
and an examples fix. A new flag with no breaking change makes this a **minor** bump: v2.6.0.

No first-time contributors — confirmed by asking `releases/generate-notes` for the unmade tag against
`dev`; every PR in the range is `@FelixKrueger` or `@dependabot[bot]`, and the response carries no
`## New Contributors` section. The policy's fifth component therefore contributes nothing here.

## 1. Two policy questions this release forced — **both amended, 2026-11-18**

Both came from drafting `.github/release-notes/v2.6.0.md`. Neither was settled by bending the body:
`.github/release-notes/README.md` and §3 of the release-notes plan now carry both clauses, so the
drafted body is policy-compliant as it stands.

**The three-highlight floor.** `.github/release-notes/README.md` says "three to six highlights". This
release has **two** that pass the policy's own test — "a user would change what they do because of
this". Reaching three means promoting a dependency bump, CI work or a docs change, which the same
paragraph forbids by name, and padding is the thing this policy exists to prevent. Recommended
amendment: *fewer when the release genuinely has fewer; the floor is a guard against a thin body, not
a quota.*

**Dependency bumps versus security fixes.** The policy says dependency bumps "are never highlights".
RUSTSEC-2026-0486 is formally a dependency bump, and it is also the single thing in this release a
user scanning advisories would act on. The exclusion exists to stop routine bumps padding the list,
not to hide security fixes. Recommended amendment: *a dependency bump is not a highlight unless it
clears a security advisory.* The draft body currently includes it on that reading.

The draft is at 884 characters, below the 1,500–2,500 range. That range now says in as many words
that it describes a release of normal size and that a genuinely small release lands under it, so the
draft is not an exception to anything.

## 2. Content and housekeeping, before the version bump

- [x] **`.github/release-notes/v2.6.0.md`** — drafted, two highlights, policy-compliant after §1's
      amendments. Verified to pass `check-release`'s condition.
- [ ] **Its CHANGELOG link points at the page root.** Once the CHANGELOG heading's date is fixed, add
      the section anchor, as `v2.4.0.md` has
      (`#version-240-release-on-9-november-2026`). A wrong anchor degrades to the page top, so this is
      precision rather than correctness — but leave it right.
- [ ] **CHANGELOG:** `### Unreleased` → `### Version 2.6.0 (Release on <date>)`, in **both** copies.
      They must stay byte-identical: `tail -n +3 CHANGELOG.md` vs
      `tail -n +9 docs/src/content/docs/reference/changelog.md`. Verified identical on 2026-11-18.
- [ ] **Decide on a CHANGELOG entry for the rustix bump.** `Unreleased` currently holds only the
      `--fastqc_raw` entry. v2.5.0's section records no dependency or security bumps either, so
      omitting it matches practice — but see §1. One sentence if yes.
- [ ] **`README.md:35`** — the FastQC feature bullet says "optional **post-trimming** quality
      reports", which predates `--fastqc_raw`. The usage example at `:118` already shows the flag, so
      only the bullet needs a few words about before-and-after.
- [ ] **`wiki/usage/` does not mention `--fastqc_raw`** — only `wiki/reference/cli-flags.md` does.
      No page there is its home: `trim-reads.md` does not cover FastQC at all. Either add a short
      FastQC section or decide the reference page is sufficient and close the item. **Deferred three
      times; decide rather than defer again.**

## 3. The release-prep commit

One commit, per §7 of `plans/11112026_release-notes-policy/PLAN.md`:

- [ ] Version in **`Cargo.toml` and `Cargo.lock`** — both read 2.5.0. `check-release` hard-fails if
      they disagree.
- [ ] The CHANGELOG section heading (both copies).
- [ ] `.github/release-notes/v2.6.0.md`.

## 4. Dispatch

- [ ] Push the 5 outstanding `dev` commits first.
- [ ] Fast-forward `master` from `dev`. Verified clean: `master` is an ancestor of `dev`.
- [ ] Push `master`.
- [ ] **Dispatch `dry_run: true` with `--ref master` before the real dispatch.** This release is the
      first to exercise the notes gate and the new concurrency expression, and a dry run runs the
      notes check while suppressing every publishing side effect. **Expect jobs 1–3 green and job 4
      skipped** — that is the correct dry-run shape, not a failure, because `smoke-test-docker` is
      dry-run-gated and job 4 lists it in `needs` with no status functions anywhere.
- [ ] Dispatch for real with `--ref master`.

Do **not** `git tag` manually between the push and the dispatch. The concurrency fix means a push can
no longer cancel the dispatch, but nothing protects against a pre-pushed tag at a different commit —
that is still a hard error by design.

## 5. Verify — four channels, not one

`gh run watch` returning is not success, and a watcher can fail while the thing it watches is fine
(one did on 2026-11-18, timing out while CI was green). Query the run's `conclusion` by branch and
workflow, then check each channel:

- [ ] Tag `v2.6.0` on the remote, at the FF commit.
- [ ] GitHub Release exists, **carrying the curated body** rather than auto-notes — this is the first
      release where that distinction is visible.
- [ ] crates.io shows 2.6.0.
- [ ] Docker `:latest` and `:2.6.0` resolve, multi-arch.

## 6. Downstream

- [ ] **bioconda** — and `about.description` was deliberately queued for the v2.6.0 bump, to avoid a
      three-platform rebuild for prose. Same PR. Note bioconda's title convention is
      `Update <pkg> to <version>`, and an upstream author cannot self-approve.
- [ ] **nf-core module container bump** — the point of the release. `--fastqc_raw` cannot be
      referenced until a release carries it, and the module then needs *only* a container bump; the
      MultiQC behaviour was verified empirically, so no pipeline-side config.
- [ ] **nf-core/methylseq** after the module.
- [ ] **Dependabot alert #29** (fflate, medium) clears at the FF — it is computed from the default
      branch and the fix is on `dev`. Still open as of 2026-11-18. Note a plain
      `gh api .../dependabot/alerts` query misses it: #29 is older than the first page, so
      `--paginate`, or query `/alerts/29` directly.

## 7. Carried over, not blocking

- [ ] **Step 4 of the release-notes plan** — `gh release edit v2.4.0 --notes-file
      .github/release-notes/v2.4.0.md`. The file now carries the three first-time-contributor
      acknowledgements byte-identically to the published body, so applying it loses nobody's credit.
      It is a public edit to a shipped release and needs go-ahead at the time.
- [ ] `plans/11122026_pr-batch-547-566/` — untracked draft comments for PRs merged without them.
      Moot; keep or drop.
- [ ] `AGENTS.md` and `SESSION_HANDOFF.md` are untracked at the repo root.

## 8. State at the time of writing

Zero open issues, zero open PRs. 913 tests passing, 0 failed. `cargo fmt --check`, clippy at
`-D warnings` and `actionlint` all clean. CI green on `ccae784` across all eleven jobs, including the
new `Lint (workflows)` job — which also settles that `go` is on the `ubuntu-latest` runner.

`.github/workflows/release.yml` is at sha256
`2b8a30f1…` — recorded because two review rounds mutated it and each restoration was checked against
a hash.
