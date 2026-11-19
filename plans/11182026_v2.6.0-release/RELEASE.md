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
- [x] **CHANGELOG anchor set** — `#version-260-release-on-18-november-2026`, confirmed against the
      generated `id=` in the built site rather than inferred from the heading.
- [x] **CHANGELOG:** heading is `### Version 2.6.0 (Release on 19 November 2026)` in both copies,
      re-verified byte-identical. **Moved from the 18th to the 19th deliberately** — the prep ran up
      against midnight, so the date is the one it ships on rather than a race. The anchor in
      `v2.6.0.md` moved with it and was re-checked against the generated `id=` in the built site.
- [x] **rustix entry added.** §1's amendment makes a security-clearing bump worth telling users
      about, so leaving it out of the CHANGELOG while listing it in the release body would have been
      inconsistent. One sentence, matching the file's format.
- [x] **`README.md:35`** — the FastQC bullet now names both flags and says the two sets reach MultiQC
      as distinct samples without pipeline-side configuration. "post-trimming" is gone.
- [x] **Wiki** — the gap was not in `usage/` at all. `wiki/reference/outputs.md` has a "Report names"
      section that did not mention FastQC in any form, so both report sets are now in its table, with
      the `_raw` alias and the MultiQC sample-naming reason. `usage/` needed nothing: `cli-flags.md`
      already carries the flag and no `usage/` page covers FastQC by design.

## 3. The release-prep commit

One commit, per §7 of `plans/11112026_release-notes-policy/PLAN.md`:

- [x] Version in **`Cargo.toml` and `Cargo.lock`** — both at 2.6.0, the lock synced by a real
      `cargo build --release` rather than by hand, and the binary reports `trim_galore 2.6.0`.
- [x] The CHANGELOG section heading (both copies).
- [x] `.github/release-notes/v2.6.0.md` — committed earlier (`f145ba3`) rather than in the bump
      commit. §7's single-commit convention is about the version bump and the notes arriving
      together; the notes landing *first* is the safe order, and 5b being dropped (§13 Open-5) means
      nothing depends on the absent arm being reachable.

**`check-release`'s own logic was run locally against the bumped tree:** version parses to 2.6.0,
`is_prerelease` false, lock and toml agree, the notes gate passes, no pre-existing `v2.6.0` tag. A GA
dispatch would clear job 0.

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

**State as of 2026-11-19.** The chain is complete except for the description PR, which never gated
anything.

- [x] **bioconda package — 2.6.0 is on the channel.** Not via our PR: `BiocondaBot` autobump opened
      #71877 about five hours after the tag and `ewels` merged it at 05:43Z, with the same sha256 we
      had computed. The `Upload` workflow succeeded, so the package built for all three platforms.
- [ ] **bioconda #71879 — still open, and only carries prose.** Once autobump landed the version, our
      PR silently reduced to `about.description` plus `build.number` 0→1; retitled
      `Update trim-galore recipe description`. `thanhleviet` reviewed (COMMENTED, not APPROVED) and
      asked for Rust third-party licence bundling, which is now in: `cargo-bundle-licenses` as a build
      dep, the generation step in `build.sh`, and a two-entry `license_file`. CI green on all
      platforms; `REVIEW_REQUIRED` stands because an upstream author cannot self-approve.
      - **Open, non-blocking:** `thanhleviet` notes `host: zlib` is unused — `Cargo.lock` carries
        `zlib-rs`/`libz-rs-sys` 0.6.8 and no `libz-sys`, so the `libzlib` run_export is spurious. It
        was not always: `host: zlib` came from #65446 when feature unification via
        `fastqc-rust` → `zip` really did pull C zlib. Left to him to fold in or defer.
- [x] **nf-core module — MERGED 2026-11-19, by `mashehu`** (nf-core/modules#13831, merge commit
      `6b4353bb6`). `bioconda::trim-galore=2.6.0`, both container strings, all four `meta.yml`
      container entries, both conda lock files and the snapshot. Verified on their master afterwards,
      including that the stub entries survived.
      - **The earlier "BLOCKED on bioconda #71879" note was wrong**, and the correction it replaced
        was also wrong. It was never blocked on *our* PR — only on a 2.6.0 package existing, which
        autobump delivered. The containers are built from `environment.yml` by
        `nf-core modules containers create`, which **only exists in nf-core/tools 4.x**; hand-editing
        `meta.yml` is not viable, because the per-architecture `build_id`/`scan_id` fields and the
        `.conda-lock/` filenames are machine-generated. That command can also half-finish on a
        transient network error, leaving arm64 absent — re-run it and check all four entries.
      - **Snapshot gotcha, cost most of the effort.** The JSON report's checksum depends on how the
        binary was invoked, because `command_line` records `argv[0]` — now
        [#578](https://github.com/FelixKrueger/TrimGalore/issues/578). A local docker run of the amd64
        image gives an absolute path; CI through Nextflow gives a bare `trim_galore`. **Take the
        values from CI's own diff rather than regenerating locally.** And note
        `nf-test --update-snapshot` **writes a snapshot even when the tests fail to start** — a run
        whose containers died hollowed both stub entries from ~46 lines to 7, which only showed up as
        a diff stat too large for the edit made.
- [ ] **nf-core/methylseq — the last link, not started.** `nf-core modules update trimgalore` in the
      pipeline repo, which pulls the merged module and updates `modules.json`. Local checkouts:
      `~/Github/methylseq` and `~/Github/methylseq-pr638`.
- [x] **Dependabot alert #29** (fflate, medium) — **`fixed`**, cleared by the FF to `master` as
      predicted. Note a plain `gh api .../dependabot/alerts` query misses it: #29 is older than the
      first page, so use `--paginate` or query `/alerts/29` directly.

## 7. Carried over, not blocking

- [x] **Step 4 of the release-notes plan — applied 2026-11-19.** v2.4.0's body is now the curated
      2,439-character version; all three acknowledgements verified live afterwards. The 16,084-byte
      original is committed at `plans/11112026_release-notes-policy/v2.4.0_superseded_body.md`,
      because regenerating it returns only two of the three contributors.
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
