# PLAN v2 — Build the docs site in PR CI

**Branch:** `ci/docs-build-pr-check` off `dev`
**Motivation:** found while merging dependabot PR #368 (`astro` 7.1.1 → 7.1.3, `@astrojs/starlight` 0.41.3 → 0.41.4, **`satori` 0.28 → 0.29**). All eight PR checks were green and **not one of them built the docs site.**
**Reviews incorporated:** `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` (dual independent, read-only). See §11 for what v1 got wrong.

---

## 1. Goal

Make a PR that changes the docs — content or dependencies — prove the site still builds, before merge.

`docs.yml` triggers only on `push` to `master`/`dev`. Pushes to `dev` auto-deploy to www.trimgalore.com, so the **first** validation of any docs change is the live deploy on the trunk.

Secondary: assert the build produced its artefacts. A satori renderer failure can leave the build exiting 0 with zero-byte PNGs, so "build succeeded" is not sufficient evidence.

Non-goal: changing what `docs.yml` deploys, or when.

---

## 2. Context

### 2.1 The gap

`docs.yml:3-10` triggers on `push` to `[master, dev]` with a `paths` filter, plus `workflow_dispatch`. No `pull_request`. Confirmed by both reviewers: `grep -n pull_request .github/workflows/*.yml` returns exactly one hit, `ci.yml:6`, and neither `ci.yml` nor `release.yml` contains any `npm`/`node`/`astro`/`docs/` reference. **`ci.yml` is the only PR-triggered workflow, and nothing on a PR builds the site.**

PR #368's eight checks were `Rust Tests` ×2, `Lint`, `Coverage`, `Reproducibility`, `Security audit`, `Validate uBAM input`, `Validate vs Perl TrimGalore` — all Rust. The npm bump was entirely unvalidated.

`satori` 0.28 → 0.29 is a 0.x **minor**, which under semver-for-0.x may break. I validated #368 by hand before merging. This plan automates that step.

### 2.2 Why the job goes in `ci.yml`, not a `pull_request` trigger on `docs.yml`

**The decisive reason, which v1 missed:** `docs.yml:6-9` carries a `paths` filter. A path-filtered workflow triggered on `pull_request` reports **nothing at all** on a PR that touches no `docs/**` path — not a skipped check, not a neutral one. If that check were ever made required, those PRs would block forever, and the standard workaround is a duplicate no-op workflow with the same name. That is worse than the "skipped reads ambiguously" objection in §2.3, and it settles the question.

Two further grounds, both real but *soluble* — v1 overstated them as decisive:

- **Permissions.** `docs.yml:12-15` grants `pages: write` and `id-token: write` at workflow level, and for same-repo `pull_request` runs the token really does receive them (the read-only downgrade applies to forks and Dependabot). But job-level `permissions:` fully overrides workflow-level, so `permissions: {contents: read}` on `build` would close it.
- **Concurrency.** `docs.yml:17-19` (`group: pages`, `cancel-in-progress: false`) would serialise PR builds behind live deploys — but `concurrency.group` accepts expressions, so a per-event group would fix it.

The deploy job would also need an event guard, which is one more way to get the deploy path wrong. That path is the one thing worth not touching.

Accepted cost: the build steps exist in two workflows, so a Node bump needs both. Mitigated by a cross-reference comment in each. Dependabot's `github-actions` group already updates pinned SHAs across all three workflow files together (observed on #367; both files use `actions/checkout@3d3c42e…` today).

**Also apply the least-privilege point to the new job:** `ci.yml` has no top-level `permissions:` block (`permissions` appears only at `ci.yml:160`, job-level on `audit`), so `docs-build` would otherwise inherit the repo-default token scope. Add `permissions: {contents: read}`.

### 2.3 Always-run, not path-filtered

No job in `ci.yml` uses a `paths` filter — the five `grep` hits are all comment prose. Every job runs on every PR, so a path-filtered job would produce a *skipped* check on Rust-only PRs, which reads ambiguously.

Note the branch-protection half of v1's argument is currently **vacuous**: verified via `gh api` that neither `dev` nor `master` is protected and `rulesets` returns `[]`. There are no required checks to add. Always-run still stands as the shape that survives protection being introduced later, but it should not be justified on a constraint that does not exist today.

### 2.4 The docs build has no dependency outside `docs/`

Checked because `docs.yml:8` lists `CHANGELOG.md` as a trigger path. It does not: `docs/src/content/docs/reference/changelog.md` is a hand-synced copy, and nothing under `docs/` reads `../CHANGELOG.md`. So `working-directory: docs` is safe.

Side effect worth recording: a PR touching only the root `CHANGELOG.md` will run this check and validate nothing new, and `docs.yml`'s `CHANGELOG.md` trigger path is effectively vestigial — it redeploys byte-identical output. That is the drift this repo already tracks as a separate follow-up (the changelog mirror automation); not this plan's job.

---

## 3. Behavior

1. On every PR to `master`/`dev`, and on every push to those branches, run a `docs-build` job. Skip on `schedule`.
2. Steps mirror `docs.yml`'s `build` job **minus** `upload-pages-artifact`: checkout → setup Node 24 with npm cache → `npm ci` → `npm run build`.
3. Assert the build produced its artefacts (§3.2).
4. Fail the job — and the PR check — if any assertion fails.

### 3.1 Paths are relative to `docs/`

`defaults.run.working-directory: docs` is in force, so **every path in the assertion is `dist/…`, never `docs/dist/…`.** v1 wrote `docs/dist/…` throughout, which would have produced an assertion pointing at a directory that does not exist from the job's cwd — and that fails **silently**, because GitHub's default shell is `bash -e {0}` without `pipefail`:

```console
$ bash -e -c 'n=$(find docs/dist/og -name "*.png" | wc -l); echo $n'   # run from docs/
find: docs/dist/og: No such file or directory
0
$ echo $?
0          # inert, and green
```

Verified. This is the single most dangerous mis-implementation available here, because the job would pass while asserting nothing.

Note `docs.yml` itself mixes conventions — `docs.yml:40,43` are cwd-relative while `docs.yml:48`'s `path: docs/dist` is repo-root-relative, because `with:` inputs ignore `working-directory`. Easy to copy the wrong one.

### 3.2 The assertion

```bash
test -d dist/og
html=$(find dist -name '*.html' | wc -l)
png=$(find dist/og -name '*.png' | wc -l)
echo "$html HTML pages, $png OG images"
test "$html" -ge 1
test "$png"  -ge 1
runt=$(find dist/og -name '*.png' -size -10240c)
if [ -n "$runt" ]; then echo "Undersized OG images:"; echo "$runt"; exit 1; fi
```

Four properties, each earned from a verified failure mode:

- **`test -d dist/og` first.** A bare `-d` test cannot be pipeline-swallowed, so a wrong path is a loud failure rather than a zero count (§3.1).
- **Counts before the size check.** With `dist/og` absent, the size check alone exits 0 — `find`'s non-zero status is discarded once its output is consumed by `$(…)`, and `bash -e` does not abort. Verified.
- **Recursive `find`, not a glob.** Of 30 OG PNGs, only **3** are at top level (`index.png`, `install.png`, `quickstart.png`); 27 are nested under `og/guide/`, `og/modes/`, `og/performance/`, `og/reference/`, `og/rrbs/`, following `docs/src/pages/og/[...route].ts`'s `slug: ${toSlug(entry.id)}.png`. v1's A2 wrote the path as a top-level glob, which would have covered 3 of 30.
- **`-size -10240c`, not `-size 0` or `-size -1k`.** `-size 0` catches only truly empty files and misses truncation. `-size -1k` **also misses a 500-byte file** — both BSD and GNU `find` round `-1k` up to a block, so a sub-block file counts as one block and escapes. Verified on BSD `find`: `-size -1k` matched only the 0-byte file, `-size -10240c` matched both it and a 500-byte truncation. The smallest real OG image is 23,221 bytes, so a 10 kB floor is a 2× margin that never needs updating as pages are added.

**Do not** use the inverted-grep form `find … -size 0 | grep -q .` — it exits 1 when nothing matches, i.e. **it fails on a healthy build**. Verified.

### 3.3 Why the OG assertion is in scope

It is the reason this job is worth more than "the build didn't crash". Verified both halves by patching `docs/src/pages/og/[...route].ts` in an isolated copy:

| Simulated failure | `astro build` exit | On disk |
|---|---|---|
| Endpoint **throws** | **1** — build aborts | no `dist/og` |
| Renderer returns an **empty buffer** | **0** — `31 page(s) built`, `Complete!` | a **0-byte** PNG |

Hard failures are already caught by build exit status. The residual gap is the silent zero-byte write — the same silent-success shape that let the `-a2` defect (#369) survive three releases behind a green parity gate.

---

## 4. Implementation outline

1. **Add a `docs-build` job to `.github/workflows/ci.yml`**, placed after `lint` and before `audit`. Job name `Docs build (Astro)`. Include `if: github.event_name != 'schedule'` and `permissions: {contents: read}` (§2.2).

2. **Copy the four build steps from `docs.yml:28-43` verbatim** — including the pinned action SHAs and their `# v7.0.1`-style comments, so dependabot's group update keeps both files in step. Set `defaults.run.working-directory: docs`.

3. **Add the assertion from §3.2** as one `run:` step. Paths relative to `docs/` (§3.1).

4. **Add a cross-reference comment** in both files, one line each per `CLAUDE.md`'s comment rule: `ci.yml`'s job notes `docs.yml` builds the same site for deploy; `docs.yml`'s `build` job notes `ci.yml` builds it for PR validation.

5. **Add a CHANGELOG entry.** v1 said to skip this on the precedent of #247 and #346; that precedent is the **reverse** of what v1 claimed. #346 (`e81fc7f`, a one-line ignore-list tweak) got nothing, but **#247 got a dedicated `#### Infrastructure (contributor-facing…)` section** at `CHANGELOG.md:720`, added at PR time in #249, documenting among other things two new CI jobs. `#### Infrastructure` sections appear three times in this CHANGELOG, so it is an established pattern. The rule is "trivial tweak → nothing; new job or new capability → a contributor-facing line", and `docs-build` is a new job. The current `Unreleased` block has only `#### Changes` and `#### Fixes`, so this means adding the heading.

6. **Fix `justfile`'s `cd Docs`** (out of scope for the gap, in scope for this PR's coherence). `justfile:86`, `:90` and `:94` say `cd Docs`; the tracked directory is `docs/` — `git ls-files` has **0** entries under `Docs/` and 154 under `docs/`. It works on case-insensitive macOS and fails on Linux. The `justfile` exists specifically for local CI parity (`CHANGELOG.md:742-747`, "#247 item 7"), so shipping a docs-build CI gate while the local equivalent is broken on Linux would be incoherent. Three one-character fixes.

---

## 5. Efficiency

Timed locally: `npm ci` 3.8 s (warm cache, 402 packages), `npm run build` 13.1 s.

On a hosted `ubuntu-latest` runner expect **1.5–2.5 min** including `setup-node`'s toolchain download and cache restore/save — v1's "roughly one minute" would visibly miss on the first run. It does not change the conclusion: the job runs in parallel with `Rust Tests (macos-latest)` at ~9 min and `validation` at ~4½, so wall-clock PR time is unchanged. `ci.yml` already starts six jobs concurrently; seven is nowhere near the limit. Actions minutes are free on a public repo.

Two points in the plan's favour v1 did not claim:

- **The npm cache is shared, not duplicated.** `setup-node`'s key derives from the lockfile hash and platform, and `docs-build` uses the same `cache-dependency-path`, so it reuses `docs.yml`'s existing entry rather than adding pressure to the cache budget.
- **Accepted duplication:** `ci.yml` also triggers on `push` to `master`/`dev`, so a docs-touching push to `dev` builds the site twice — once here, once in `docs.yml` before deploy. ~13 s of duplicated runner time, free tier. Left deliberately: an `if: github.event_name == 'pull_request'` guard buys nothing and adds a condition.

---

## 6. Integration

**Reads:** `docs/package-lock.json`, `docs/package.json`, `docs/src/**`. **Writes:** nothing outside the runner.

**New check on every PR:** `Docs build (Astro)`, taking the count from 8 to 9. No branch protection exists (§2.3), so nothing needs adding there today.

**Does not touch** `docs.yml`'s deploy path, triggers, permissions or concurrency, beyond one comment line. Note that comment lands on a `docs.yml` trigger path, so merging to `dev` fires a real deploy of byte-identical content — harmless, but not a surprise.

---

## 7. Assumptions

- **A1.** `npm run build` is the whole validation `docs.yml` performs before deploying — the only step between build and artifact upload is the upload. Verified by both reviewers.
- **A2.** satori's output is **30 PNGs under `dist/og/**/*.png`, of which 3 are top-level and 27 nested**. Verified by building. (v1 wrote this as a flat `docs/dist/og/*.png` glob; see §3.2.)
- **A3.** Node 24 matches `docs.yml:35` and its ≥22.12.0 comment.
- **A4.** PR CI has network for `npm ci` — `docs.yml` already does this on push, and `validation` fetches conda packages.
- **A5.** **Corrected.** v1 claimed `ci.yml` jobs have no `needs:` between them; that is false — `validation` (`ci.yml:247`) and `validation-ubam` (`ci.yml:836`) both declare `needs: [rust-tests, lint]`. The conclusion survives on a narrower premise: *no existing job `needs:` the job being added, and `docs-build` needs nothing, so the dependency graph is unchanged.*
- **A6.** **Corrected.** v1 claimed every `ci.yml` job carries `if: github.event_name != 'schedule'`; `audit` (`ci.yml:157`) deliberately has none — that is why the weekly cron exists. `docs-build` follows the majority and takes the guard.
- **A7.** The lockfile carries Linux-specific optional deps for `sharp` and `@resvg/resvg-js`. Proven by `docs.yml` building on `ubuntu-latest` today; would break only if a lockfile were regenerated macOS-only.

---

## 8. Validation

**Rehearse V3/V4 from `docs/`, not the repo root** — the job's cwd. v1's rows used `docs/`-prefixed paths, so they would have passed from the repo root while the workflow step was inert (§3.1). This is the single most important instruction in this section.

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | The job builds the site on a PR | Push the branch, open a PR | `Docs build (Astro)` passes |
| V2 | **It fails on a broken docs build** | Scratch commit removing a required `title:` from one page, or a bad `slug:` in `astro.config.mjs`'s sidebar | Check **fails**. v1 proposed an unknown frontmatter key — **verified to build with exit 0**, so v1's V2 would have gone green and proved nothing. The sidebar-slug break is the more realistic real-world case (a page rename that misses the config) |
| V3 | Fails on a **nested** empty PNG | From `docs/`: build, `: > dist/og/guide/adapters.png`, run the §3.2 script | Non-zero, naming the file. v1 targeted `dist/og/index.png` — top-level, so it would pass against a non-recursive check |
| V3b | Fails on a **truncated** PNG | From `docs/`: `head -c 500 dist/og/index.png > dist/og/guide/adapters.png`, run the script | Non-zero. `-size 0` and `-size -1k` both miss this; `-size -10240c` catches it |
| V3c | Fails when `dist/og` is **absent** | From `docs/`: remove `dist/og`, run the script | Non-zero from `test -d`. Without it, the size check alone exits 0 |
| V4 | Passes on a clean build | From `docs/`: `npm ci && npm run build`, run the script | Pass; reports 32 HTML pages, 30 OG images |
| V5 | No existing check disturbed | Compare against #370's 8 checks | 9 checks, the other 8 unchanged |
| V6 | YAML valid before pushing | `ruby -ryaml -e 'p YAML.load_file(".github/workflows/ci.yml")["jobs"].keys'` | `docs-build` appears between `lint` and `audit`. v1 proposed `python3 -c "import yaml"` — **no `yaml` module on this machine**, and no `yamllint`/`actionlint` on `$PATH` |
| V7 | Skips on schedule | Read the `if:` guard | Matches the majority of sibling jobs (§7 A6) |
| V8 | `just docs-build` works | Run the `justfile` recipe after the §4 step 6 fix | Builds; previously `cd Docs` |

V2, V3 and V3c are the rows that would quietly not hold, and all three were wrong in v1. This plan exists because a suite of green checks did not exercise the thing at risk; a gate whose own tests cannot fail would repeat that in a new place.

---

## 9. Questions or ambiguities

**Resolved since v1:**

- **R1 — branch protection.** Neither trunk is protected; `rulesets` is empty. Nothing to add.
- **R2 — CHANGELOG.** An entry *is* the precedent for a new job (§4 step 5).
- **R3 — hosting.** Settled by the path-filter argument (§2.2), not by the permissions/concurrency grounds v1 leaned on.

**Open (assumption taken, no blocker):**

1. **Duplication with `docs.yml`.** Taken: accept it with cross-reference comments. A composite action shared by both would remove it but adds a file and an indirection for less than it saves.
2. **Whether to also assert pagefind and the sitemap.** The clean build emits `dist/pagefind/pagefind-entry.json` and `dist/sitemap-index.xml`; if Starlight's pagefind step degrades across a version bump, site search dies with a green build — the same failure shape as the OG images, at the cost of one `test -s` each. Taken: out of scope, recorded as a cheap follow-up.
3. **Broken internal links.** Needs `starlight-links-validator`, a new dependency. Out of scope.

---

## 10. Self-Review

**Logic.** Traced both hosting options; `ci.yml` wins, and the reason that settles it is the path-filter-reports-nothing behaviour, not the permissions and concurrency grounds v1 relied on — both of those are soluble with one line each. Confirmed the new job cannot perturb the existing graph on the corrected premise in A5.

**Adjusted after review.** The assertion changed shape entirely: recursive rather than a top-level glob, `test -d` first, byte-unit size floor rather than `-size 0`, and all paths relative to `docs/`. Three of the four validation rows meant to prove the gate works were wrong and are replaced.

**Edge cases.** Docs-unrelated PRs (job runs anyway — cheap, avoids skipped checks); a PR deleting docs pages (assertions are "at least one"); `schedule` runs (guarded); fork PRs (needs no elevated permissions, unlike the `docs.yml` option); npm registry outage (fails the check — the same exposure `docs.yml` already has); `dist/og` absent (V3c); truncated-but-non-empty PNG (V3b).

**Remaining risks.**

- *Low:* the duplicated build drifts from `docs.yml`. Mitigated by cross-reference comments and dependabot updating both files' pins together.
- *Low:* a valid PNG that renders blank or garbled — outside any cheap check, and accepted.
- *Low:* 1.5–2.5 min per PR, parallel with a 9-minute job, free tier.

---

## 11. Revision history

### v2 — 2026-07-30, after dual independent plan review

Both reviewers reproduced the core premise and confirmed the OG silent-failure mode empirically. Both were read-only; the tracked tree was untouched.

**Two Criticals, from Reviewer B, both on the checks meant to prove the gate works:**

1. **V2's deliberate break did not break the build.** v1 proposed a malformed frontmatter key; verified — Starlight is not strict about *unknown* keys, so `npm run build` exits **0**. V2 would have gone green, from which both available conclusions are wrong. Replaced with removing a required `title:` (verified: exit 1, "does not match collection schema") or a bad sidebar slug.
2. **The assertion could be a silent no-op that §8 could not detect.** `working-directory: docs` means paths must be `dist/…`, but v1 wrote `docs/dist/…` everywhere — and a wrong path under `bash -e` without `pipefail` prints to stderr and exits **0**. Worse, V3/V4 were to be rehearsed locally from the repo root, a different cwd from the job, so both could pass while the step was inert. Fixed by stating the cwd once, leading with `test -d dist/og`, and rehearsing from `docs/`.

**Reviewer A's Important findings:**

3. **The OG path form covered 3 of 30 images.** v1's A2 wrote a top-level glob; only `index.png`, `install.png` and `quickstart.png` are top-level, and 27 are nested. A verified that a non-recursive check *passes* with a nested 0-byte PNG present, and that v1's V3 — which truncated a top-level file — would pass against it.
4. **The emptiness check is vacuous when `dist/og` is absent.** Both candidate formulations exit 0. Fixed by ordering and `test -d`.
5. **`-size -1k` misses a 500-byte file.** Block rounding; verified. `-size -10240c` catches both truncation and emptiness at no maintenance cost.
6. **The decisive hosting argument was missing.** A path-filtered `pull_request` workflow reports *nothing* on non-docs PRs — worse than skipped. This converts §2.2 from a judgement call into a clear-cut one, and correspondingly downgrades v1's permissions and concurrency grounds, both of which are soluble.

**Corrections to v1's own claims:**

| v1 said | Actually |
|---|---|
| `ci.yml` jobs have no `needs:` (A5) | `validation` `:247` and `validation-ubam` `:836` both do. Conclusion survives on a narrower premise |
| Every job carries the `schedule` guard (§2.2) | `audit` `:157` deliberately does not — that is why the cron exists |
| #247 and #346 set a precedent for **no** CHANGELOG entry | #346 yes; **#247 got a `#### Infrastructure (contributor-facing)` section** at `CHANGELOG.md:720`. Three such sections exist. A new job takes an entry |
| Always-run avoids complicating branch protection | Neither trunk is protected and `rulesets` is empty, so that argument is currently vacuous. Always-run stands on the thinner "skipped reads ambiguously" plus future-proofing |
| ~1 minute per PR | 1.5–2.5 min on a hosted runner including `setup-node` |
| V6: `python3 -c "import yaml"` | No `yaml` module on this machine; use the `ruby -ryaml` form |

**Added:** `permissions: {contents: read}` on the new job (`ci.yml` has no top-level block); V3b, V3c, V8; the `justfile` `cd Docs` fix (three sites, 0 tracked files under `Docs/`, breaks on Linux — the same case mismatch already recorded as an environment gotcha, here a real committed bug); §2.4 recording that the docs build has no dependency outside `docs/`; the shared-npm-cache and double-build-on-push notes in §5.

**Recorded as out of scope:** asserting pagefind and the sitemap (§9 item 2) — same failure shape, one `test -s` each; internal-link validation, which needs a new dependency.

### v1 — 2026-07-30

Initial plan. Correct on the gap, the hosting conclusion, the OG-assertion rationale, A1/A3/A4, and §2.3's no-paths-filter claim — all confirmed by both reviewers.
