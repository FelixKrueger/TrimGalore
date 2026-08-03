# PLAN — Generate the published changelog page from `CHANGELOG.md`

**Branch:** `docs/generate-changelog-page` off `dev`
**Problem:** `docs/src/content/docs/reference/changelog.md` is a hand-maintained copy of the root `CHANGELOG.md`. Nothing enforces the copy, so it silently lags.

---

## 1. Goal

Make the published changelog page a **build artefact of `CHANGELOG.md`**, and delete the checked-in copy. Drift then becomes unrepresentable rather than merely detectable.

Non-goal: changing the root `CHANGELOG.md`'s format, or how releases are cut.

---

## 2. Context

### 2.1 The current drift, measured

| File | Role | Lines | Last touched |
|---|---|---|---|
| `CHANGELOG.md` | source of truth | 1331 | 2026-07-31, `7201f03` |
| `docs/src/content/docs/reference/changelog.md` | the page at www.trimgalore.com/reference/changelog/ | 1141 | **2026-07-25**, `4f86c0f` |

Diffing the bodies — page lines 1-9 (frontmatter, banner, blanks) and root lines 1-3 (`# Trim Galore Changelog` and its two blank lines) stripped from each side:

- **196 lines exist only in the root file**
- **0 lines exist only in the docs page**
- **no `###` release section is missing** from the page

The page is purely *behind*, never divergent. Nothing has ever been hand-authored into it, so a one-way generate-and-overwrite is **lossless** — there is nothing to merge. That is the fact that makes this plan safe.

**On the 196 vs 199.** An earlier count of this gap said 199, from a diff that stripped the page's frontmatter and banner but compared against the *unstripped* root file. The extra 3 are the root's H1 and its two blank lines — which the current page does not carry and the generator must therefore remove (§3.0). The H1 was hiding inside the drift figure. 196 is the body gap; 199 is the raw line difference between the two files' post-frontmatter content.

Currently unpublished: `--clump_only` uBAM in/out, the `--dont_gzip` + `--output-format ubam` rejection, both `--phred64` fixes, the mixed FASTQ+BAM pair diagnosis (#363), the `-a2` fix (#369), and the whole `#### Infrastructure (contributor-facing)` section. #369 is the notable one — a public issue whose fix is invisible on the public changelog.

This has happened before: #359's own description records the page as having been *"frozen at v2.1.0-beta.5"* until that PR resynced it.

### 2.2 Why the copy exists at all

Astro content collections require frontmatter. `docs/src/content.config.ts`:

```ts
docs: defineCollection({ loader: docsLoader(), schema: docsSchema() }),
```

`docsSchema()` requires a `title`. A root `CHANGELOG.md` carrying `---\ntitle: …\n---` would look wrong on GitHub and in release tooling, so the page cannot simply *be* the root file. The copy bridges that, and the frontmatter is the entire obstacle.

### 2.3 `docs.yml` is already wired for this

`docs.yml:8` lists `CHANGELOG.md` in its trigger paths, so editing the root file **already triggers a docs deploy** — which currently redeploys byte-identical output, because nothing propagates the change. The trigger was wired for an automation that was never built. This plan makes it mean what it looks like it means.

### 2.4 This is only safe because #371 landed

Deleting the tracked page means the build **must** run the generator or the page vanishes. Until today nothing built the docs on a PR, so a broken generator would have been discovered by the live deploy. #371's `docs-build` job now builds the site on every PR, so a generator that fails, or a page that goes missing, fails the check before merge.

---

## 3. Behavior

1. A generator script reads `CHANGELOG.md` from the repo root, **strips the leading `# Trim Galore Changelog` H1 and the blank lines following it** (§3.0), prepends frontmatter and a short banner, and writes `docs/src/content/docs/reference/changelog.md`.
2. It runs before `astro build` **and** before `astro dev` — `astro dev` reads the collection too, so a dev server without it would 404 the page.
3. The output path is exactly `reference/changelog`, matching `astro.config.mjs:105`'s `slug: 'reference/changelog'`.
4. The output is **gitignored and untracked**. The only committed copy of the changelog is the root file.
5. The generator fails loudly — non-zero, with a message naming the path — if `CHANGELOG.md` is missing or empty. A silent empty page would be worse than a stale one.
6. It is idempotent: running twice produces a byte-identical file. **Write only when the content differs** — read the existing file first and skip the write if unchanged. Both reviewers suggested this: it makes idempotency true at the filesystem level (`mtime` unchanged, not merely bytes equal), so a repeat `npm run dev` is a genuine no-op and does not churn `astro dev`'s watcher or invalidate the `.astro/` content-layer cache. Two lines.
7. It asserts its own output before exiting — see §4 step 2a. A generator that silently emits a truncated page is the failure mode with no signal.

### 3.0 The leading H1 must be stripped

`CHANGELOG.md` line 1 is `# Trim Galore Changelog`. Starlight renders frontmatter `title` as the page's `<h1>`, so a generated body that keeps that line produces **two H1s**. The current hand-maintained page strips it, which is why this has never surfaced.

The generator removes the first line if it is an H1, plus any blank lines immediately after it. Conditional rather than a fixed `sed 1,3d`, so a future reformat of the root file's header does not silently drop a real entry.

This is the one transformation the generator performs. Everything after it is verbatim, which keeps §3.2's "the generator does no transformation" reasoning intact for the body — the constraint being that the strip is anchored to the very top of the file and cannot reach content.

Consequence for validation: V1 cannot diff the generated body against the raw root file. See V1 as restated.

### 3.1 Frontmatter and banner

```yaml
---
title: Changelog
description: Release history for Trim Galore. The canonical source is CHANGELOG.md in the repo root.
---
```

The existing banner says *"The version below is a copy synced with the docs"*, which becomes misleading once generated — it implies a manual step that no longer exists. Reword to state that the page is generated from the root file at build time. Register per `feedback_docs_register`: plain statement, no filler.

**Retarget the link from `master` to `dev` at the same time.** The current banner links to `blob/master/CHANGELOG.md`. `docs.yml:5` deploys on push to **both** `master` and `dev`, and `dev` is the branch that auto-deploys to www.trimgalore.com — and is ahead of `master` (37 commits at time of writing; the argument needs only "nonzero", so the count is deliberately not load-bearing). So the generated page body will normally come from `dev` while the link points at `master`, which does not contain the `### Unreleased` entries the page is displaying.

Automation makes this worse rather than better: today the stale page happens to sit roughly at `master`, so the mismatch is masked. Once the page tracks the deploying branch, a reader following the link lands on *older* content than the page they came from — the exact failure the banner exists to prevent.

`dev` is the right target because it is what gets published. Recorded as a decision rather than a detail, because it means the banner link and the deploy branch are now coupled: if the deploying branch ever changes, this link has to change with it (A6).

### 3.2 Content compatibility

`[REVISED BY REVIEW]` — the original version of this section was wrong on **both** of its counts, and its governing rule was falsified. Restated below in a form that survives.

The root file must survive Astro's Markdown pipeline. What actually holds:

- **Braces and angle brackets are safe in `.md` regardless of code-span context.** Not, as originally claimed, "safe *because* they are all inside code spans" — they are not. `CHANGELOG.md:1145` (`… ending in .{INT}bp_3prime.fq(.gz).`) and `:1185` (`e.g. -a "A{10}" to trim poly-A tails`) are bare prose, and both already render literally on the published page. The empirical proof is stronger than the structural claim: these lines are in the *already-synced* portion and have been rendering correctly for months.
- **Counts, with provenance so they reproduce:** `grep -c '{' ` = 17; `grep -c '<[a-zA-Z]'` = 22; plain `grep -c '<'` = **24** (the two extra are bare comparisons such as `matches < 100` at `:722`, correctly escaped to `&#x3C;` in the output). Of the 196 new lines, 4 carry `<…>` — three in an indented fence at `:39-41`, one in a backtick span at `:50`.
- **One of the 24 is real HTML, not markup-in-prose.** `CHANGELOG.md:1164` is a raw `<img title="Multi-threading benchmark" … src="docs/Images/pigz_bench.png" >` outside any fence, which renders as a **live element**. `docs/Images/` exists nowhere in the repo; the real asset is `docs/src/assets/screenshots/pigz_bench.png`. The image is broken on the published site today, independent of this plan. **Fixed at source by the sibling plan — see §3.4.** The generator stays verbatim.

**The rule "anything valid on GitHub renders here" is false**, and A5 no longer states it. This pipeline is not GitHub-flavoured Markdown: `astro.config.mjs:16-20` adds `remarkMath` + `rehypeKatex`, so `$…$` is a math delimiter. That is a live bug in the current content and is the subject of §3.4. Treat the pipeline's divergences as an explicit list, not a general guarantee:

| Construct | GitHub | This pipeline |
|---|---|---|
| `{`, `<…>` in prose | literal | literal (safe) |
| raw `<img …>` in prose | literal element | literal element — relative `src` resolves against the page URL |
| `$…$` | literal | **consumed as inline math** |

**The page must stay `.md`, never `.mdx`.** The reason is stronger than originally stated: under MDX the failure is a **hard build error**, not a mis-render. `CHANGELOG.md:1164`'s `<img … >` is not self-closed, which is invalid JSX. Name *that* line in the generator comment — the originally cited `<input.bam>` is inside a fence and would be safe.

### 3.3 What the first run does

Generating closes the 196-line gap as a side effect, so **no separate manual sync is needed** — the automation is the sync. That is worth stating because the obvious instinct is to hand-sync first and automate second, which would just create a diff to reconcile.

### 3.4 Dependency: the render fixes must land first

Two defects in the current content are made worse by this plan, and are fixed separately in `plans/docs-render-fixes/PLAN.md`:

1. `$…$` KaTeX consumption (`CHANGELOG.md:451` is already mangled on the live site).
2. The broken `<img>` path at `CHANGELOG.md:1164`.

**Why the ordering matters.** Today a human hand-copying the changelog is a weak reader of the rendered result. This plan makes the page an *unattended* artefact, so every future `$` in a changelog entry mangles with nobody in the loop. Automating first would ship a known silent-corruption path and remove the last human check on it.

Neither fix touches the generator, `package.json`, or the collection — they are `astro.config.mjs` and `CHANGELOG.md` edits — so the two plans do not conflict. This plan should merge *after* that one.

---

## 4. Implementation outline

1. **Add `docs/scripts/sync-changelog.mjs`.** Follows the existing `docs/scripts/generate-png-logos.mjs` precedent — same directory, same plain-Node shape, no new dependency. Resolve paths with `node:path` and `import.meta.url` rather than assuming cwd, so it works from either the repo root or `docs/`.

2. **Wire it into `package.json` scripts**, before both `build` and `dev`:
   ```json
   "sync:changelog": "node scripts/sync-changelog.mjs",
   "build": "npm run sync:changelog && astro build",
   "dev":   "npm run sync:changelog && astro dev",
   "start": "npm run sync:changelog && astro dev",
   ```
   Prefer explicit `&&` over npm's `pre*` lifecycle hooks: `pre` hooks are silent when they fail to fire and are easy to miss when reading `package.json`.

2a. **Assert inside the generator, after writing.** `[ADDED BY REVIEW]` Both reviewers independently identified this as the highest-value change in the plan. V1/V1a are one-time manual checks: they protect the first commit and nothing after it. A later edit that truncates the output — an off-by-one in the strip, a partial read — produces a file that exists, resolves the sidebar slug, exits 0, satisfies `ci.yml`'s `test "$html" -ge 1`, and deploys a **half-length changelog**. Two assertions on data already in memory close that:

   - emitted body line count **==** input line count after the §3.0 strip → fails loud on any truncation, and subsumes behaviour item 5's "or empty" with a much stronger floor;
   - **zero** `/^# /` lines in the emitted body → fails loud on a duplicate H1, *including* the insert-above-the-H1 case that a line-1-anchored strip cannot see (§3.0, A7).

   These convert the plan's two acknowledged blind spots into build failures on every PR through the existing `docs-build` job, with no new CI wiring.

3. **Untrack the generated page and ignore it.** `git rm --cached docs/src/content/docs/reference/changelog.md`, then add `src/content/docs/reference/changelog.md` to `docs/.gitignore`. The comment must say two things: which script generates it, **and** that `CHANGELOG.md` in the repo root is where to grep. `ripgrep` honours `.gitignore`, so after this change `rg '<changelog phrase>' docs/` returns nothing and a contributor will reasonably conclude the page does not exist.

4. **Reword the banner** per §3.1.

4a. **Emit `editUrl` frontmatter.** `[ADDED BY REVIEW]` `astro.config.mjs:36-39` sets `editLink.baseUrl: '…/edit/master/docs/'`, so the built page currently links to `/edit/master/docs/src/content/docs/reference/changelog.md`. Once that file is untracked and deleted the path exists on no branch — and GitHub's `/edit/` on a missing path opens the **new-file editor**, so "Edit page" would invite a contributor to recreate the tracked copy this plan just deleted, resurrecting the exact drift being eliminated. `editUrl` is a supported field (`@astrojs/starlight/schema.ts:30`, `z.union([z.url(), z.boolean()]).optional().default(true)`), so emit:

   ```yaml
   editUrl: https://github.com/FelixKrueger/TrimGalore/edit/dev/CHANGELOG.md
   ```

   `editUrl: false` is the acceptable minimum. Same branch reasoning as §3.1.

5. **CHANGELOG entry** under the existing `#### Infrastructure (contributor-facing)` heading in `Unreleased` — added yesterday in #371, so the heading already exists. Contributor-facing, per the #247 precedent.

6. **Rewrite `docs/README.md` — mandatory, not conditional.** `[REVISED BY REVIEW]` The original step said "if that file documents the build; check first". It does, in a dedicated section, and three separate statements become false:

   - `:37-41` — the whole `## Sync with the repo CHANGELOG.md` section. It calls the page "a copy of the top-level `CHANGELOG.md` … and one inline image path rewritten" and instructs "Keep them in sync when releasing." Rewrite as: generated at build time by `scripts/sync-changelog.mjs`; do not edit, do not commit.
   - `:41` — *"and one inline image path rewritten"* names a **second transformation the plan did not know about**. Both reviewers verified it is **stale**: `CHANGELOG.md:1164` and the page's line 974 are byte-identical, and the diff has zero change hunks, which independently rules out any surviving rewrite. A hand-sync destroyed it at some point — a neat illustration of this plan's thesis. §3.0's "one transformation" therefore survives, but it survives by luck and the check is worth recording.
   - `:18` — *"The site auto-deploys from the `master` branch"*, already false per `docs.yml:5` `[master, dev]`. Same `master`-vs-`dev` confusion §3.1 diagnoses for the banner, sitting in the file being edited.

   Also note the deviation from house style here: `docs/public/logos/*.png` are generated by `npm run logos` and **committed**, so the repo's existing convention for generated artefacts is generate-and-commit. This plan is generate-and-ignore, deliberately — committing the artefact is the drift being eliminated. Say so, or the next contributor "fixes" the inconsistency in the wrong direction.

7. **Add one line to `ci.yml`'s existing `Assert build output` step.** `[ADDED BY REVIEW]` See §9 Q1, which this overturns:

   ```bash
   test -f dist/reference/changelog/index.html
   ```

---

## 5. Why a generator script, not a custom content loader

The more elegant option is a custom Astro content loader that injects the root file as a `docs` entry with no file on disk at all. Rejected for now:

- Starlight's `docsLoader()` is not documented as composable. Wrapping or merging it means depending on internals.
- Starlight 0.41.x is actively moving — 0.41.3 → 0.41.4 landed yesterday in #368. A loader built on internals is the kind of thing that breaks on a patch bump, and this repo takes dependabot bumps continuously.
- A 30-line Node script is boring and inspectable; a contributor debugging a missing page can read it in one sitting.

Revisit if Starlight ever exposes a documented composition point. Recorded here so the choice reads as deliberate.

Both reviewers independently endorsed this rejection and strengthened it: `docsLoader()` is imported from `@astrojs/starlight/loaders`, a subpath export whose return shape is not a documented contract; and the plain-file approach gets a **hard build failure for a missing page for free** from Starlight's own sidebar resolution, which a loader-based design would have to reproduce deliberately.

### 5.1 Two further alternatives, both rejected

`[ADDED BY REVIEW]` Both reviewers proposed the first of these independently, and noted §5 read as if the loader were the only option considered.

**Keep the page tracked, add a CI drift check.** Generate into the tracked path and have `docs-build` fail if the working tree is dirty afterwards — the `cargo fmt --check` pattern, already idiomatic here:

```bash
npm run sync:changelog && git diff --exit-code -- docs/src/content/docs/reference/changelog.md
```

This is the **nearest rival**, not the loader. In its favour: the changelog stays browsable and greppable under `docs/`, the Starlight edit link keeps working (so §4 step 4a becomes unnecessary), and nothing is deleted. Against: every `CHANGELOG.md` edit then *requires* a second committed file, which is precisely the manual step this plan exists to remove, and the failure mode is a red PR rather than an impossible state. **Rejected** — "drift unrepresentable" beats "drift detectable" — but it is the closer call, and it is the option that would have avoided §4 step 4a entirely.

**A plain `src/pages/` route outside the docs collection.** Sidesteps the frontmatter requirement altogether, but loses sidebar integration and Starlight page chrome, and would need the sidebar entry switched from `slug:` to `link:` — which per §9 Q1 also discards the build-time slug validation this plan leans on. Worse on both counts. **Rejected.**

**A symlink** from the content path to the root file fails on the frontmatter requirement (§2.2) exactly as the copy does, and adds a Windows-checkout hazard. Not viable.

---

## 6. Integration

**Reads:** `CHANGELOG.md`. **Writes:** `docs/src/content/docs/reference/changelog.md` (gitignored).

**Order:** before `astro build`/`astro dev`, every time.

**Visible effects:**

1. The published changelog page gains the 196 missing body lines on first deploy — including the `-a2` fix.
2. The banner's source link moves from `master` to `dev`, matching the branch the page is generated from (§3.1).
3. `docs/src/content/docs/reference/changelog.md` disappears from `git status` and from the repo.
4. `docs.yml`'s `CHANGELOG.md` trigger path becomes functional: a root-only changelog edit now changes the deployed site.
5. A contributor running `astro build` directly — or `npm run astro build`, which is the same passthrough — gets **no site at all**, not merely no changelog page. `[REVISED BY REVIEW]` Verified two ways: Reviewer A deleted the page and built, getting exit 1 and `AstroUserError: The slug "reference/changelog" … does not exist` thrown from `getSidebar` while rendering `/404.html`; Reviewer B found the throw site at `@astrojs/starlight/utils/navigation.ts:145-159`. Because `getSidebar` runs for every page, a missing page breaks the whole build rather than degrading one route. This is **fail-loud**, so §4 step 6's README note is a convenience rather than a mitigation — worth keeping only because the error names a *slug*, not a missing generator, which is obscure.
6. The Starlight "Edit page" link retargets to `CHANGELOG.md` on `dev` rather than a path that no longer exists (§4 step 4a).

**Interacts with the release flow:** at release time the root `CHANGELOG.md` gets its `### Version X.Y.Z` heading, and the page follows automatically on the next deploy. That removes one manual step from the release checklist — worth noting since `CITATION.cff` and five version strings are already on it.

---

## 7. Assumptions

- **A1.** The page is a pure lagging copy with nothing hand-authored. **Verified**: 196 body lines root-only, 0 lines page-only, no missing sections — with the root H1 stripped from both sides (§2.1). This is what makes overwrite-without-merge safe; re-verify before deleting the tracked file, since the whole plan rests on it. When re-verifying, **prove the diff can fail first** (append a sentinel line and confirm it is reported): a `diff` using process substitution fails under a sandbox with `Operation not permitted` and still reports zero differences.
- **A2.** `docsSchema()` requires only `title` (`@astrojs/starlight/schema.ts:15`, `z.string()`, no default); `description` is optional. Keeping `description` is **not cosmetic**: `docs/src/pages/og/[...route].ts:131-141` calls `getCollection('docs')` and derives `description: entry.data.description || …`, so it feeds the per-page satori OG image that `ci.yml` asserts on. Dropping it would regress a checked artefact.
- **A3.** `astro dev` reads the same collection as `astro build`, so both need the generator. **Verified**: with the page removed, `astro dev` returns **HTTP 500 on every route** (`/reference/changelog/` *and* `/`), both bodies carrying `AstroUserError`.
- **A4.** `[REVISED BY REVIEW — was wrong]` The page path is referenced in **two** places, not one: `docs/astro.config.mjs:105` and **`docs/README.md:39`**. Both reviewers found the second independently. `docs/README.md:37-41` is a whole section about this file and becomes false — see §4 step 6, which is consequently mandatory. Grep from the **repo root** before deleting; note that `git grep` and `ripgrep` are the right tools here and that any wrapper honouring `.gitignore` will stop finding the page after step 3.
- **A5.** `[REVISED BY REVIEW — the original rule was false]` The root file contains no construct that breaks `.md` processing **except the two catalogued in §3.4**, which are fixed at source by the sibling plan. The original formulation — "anything valid in the root file that renders on GitHub should render here" — is **demonstrably untrue**: this pipeline adds `remarkMath`, so `$…$` is consumed as math and `CHANGELOG.md:451` is already mangled on the live site. Do not reason from a general GitHub-equivalence rule; reason from §3.2's explicit divergence table. The generator still transforms only the leading H1 (§3.0).
- **A6.** `dev` is the branch whose push publishes the site. **Verified**: `docs.yml:5` lists `[master, dev]` and `dev` auto-deploys to www.trimgalore.com. Both the banner link (§3.1) and `editUrl` (§4 step 4a) are coupled to this; if the deploying branch changes, both must change. Nothing enforces the coupling — it is a comment in the generator, not a check. *(Reviewer A proposed removing the coupling instead by emitting `blob/<GITHUB_SHA>/CHANGELOG.md`, falling back to `dev` locally, so the link points at exactly the bytes the page was built from. Deferred — see §9 Q4.)*
- **A7.** `CHANGELOG.md` begins with exactly one H1 and it is the file title, not content. **Verified**: one ATX H1 at line 1, and no setext headings anywhere (`^=+$` and `^-{3,}$` both return nothing), so V1a's `grep -c '^# '` is adequate against the current file — though blind to a setext H1 in principle. The generator's strip is conditional on line 1 being an H1 (§3.0), which degrades safely for a *reformatted* header but **not** for a line inserted *above* it — see §9 Q5. §4 step 2a's zero-`^# ` assertion is what actually closes that.
- **A8.** `[ADDED BY REVIEW]` Line endings are stable and pose no hazard to V1/V2: both files end in a single `\n`, neither contains `\r`, and `.gitattributes` sets only `linguist-*` attributes (no `text`/`eol`).

---

## 8. Validation

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | Generated body matches the root file, H1 aside | Generate; strip frontmatter + banner from the output; `diff` against `CHANGELOG.md` **with its leading H1 and following blank lines removed**. Prove the diff can fail first (append a sentinel, confirm it is reported) | No differences |
| V1a | **Exactly one H1 on the rendered page** | `grep -c '^# ' ` the generated body; then count `<h1` in `dist/reference/changelog/index.html` | 0 in the body (Starlight supplies it); exactly 1 in the HTML |
| V1b | Banner link points at the deploying branch | `grep` the generated body for the source link | Resolves to `blob/dev/CHANGELOG.md`, not `master` (§3.1, A6) |
| V2 | Idempotent | Run the generator twice, compare | Byte-identical |
| V3 | `npm run build` produces the page | Remove the generated file, `npm run build`, then `test -f dist/reference/changelog/index.html` | The file exists. `[REVISED]` The old expectation was "32+ HTML pages" — today's count is *exactly* 32, so it has zero headroom, and 32 pages can be met while this page is absent. Target the artefact, not the count |
| V4 | **The page renders, not just builds** | Inspect `dist/reference/changelog/index.html` for content from the newly-added lines (e.g. the `-a2` entry); assert **no `application/x-tex`** anywhere in it; assert **no `src="docs/`** survives; confirm `A{N}` renders literally | Content present, correctly escaped, no math nodes, no relative image src. `[REVISED]` The old form checked only `{` and `<input.bam>` and would have caught **neither** live bug in §3.4 |
| V5 | `npm run dev` also generates | Remove the generated file, start `dev`, request the page and `/` | With the generator: 200. Without it: **every route 500s** with `AstroUserError`. `[REVISED]` The old expectation was "200, not 404" — the observed failure is a 500 on all routes, so a tester hunting a 404 would misfile it |
| V6 | **The build fails loudly on bad input** | Rename the root file (missing); `chmod 000` it (unreadable); truncate it (partial) | Missing and unreadable → non-zero, message names the path. Partial → caught by §4 step 2a's line-count assertion. `[REVISED]` The old form tested missing input only |
| V7 | The file is untracked and ignored | `git status --short docs/` from the **repo root** (pathspecs from inside `docs/` silently match nothing — the `Docs`/`docs` case trap) | Page absent from `git status`; `git ls-files` does not list it |
| V8 | Sidebar link still resolves | Build, follow `reference/changelog` from the sidebar | Page loads at the same URL as before |
| V9 | The 196-line gap is closed | Compare post-strip **body line counts**, and **`^#### ` counts** (40 root vs 39 page today), generated page vs root | Equal. `[REVISED — the old check could not fail]` `^### ` is already **49 in both files** on the stale page, so the original V9 was satisfied before any change and was blind to the one heading actually missing. `^####` is the count that discriminates; the line count is the one that catches truncation |
| V10 | `docs-build` CI passes | Push, open a PR | The #371 job goes green, now including the §4 step 7 `test -f` assertion |
| V11 | `[ADDED]` The edit link is not a dead path | `grep` `dist/reference/changelog/index.html` for the edit-link `href` | Points at `edit/dev/CHANGELOG.md`, **not** `edit/master/docs/src/content/…` (§4 step 4a) |

V6, V7 and V1a remain the three that would quietly not hold — but §4 step 2a now converts V1's and V1a's content into permanent generator assertions, which is the difference between validating the first commit and validating the mechanism. V6 because a generator that writes an empty file on missing input turns a stale page into a blank one; V7 because a still-tracked file would keep drifting while looking automated; V1a because a duplicate H1 builds, deploys, and renders with no failure signal at all.

**A harness note that generalises A1's.** Any *negative* result in this repo needs the harness proven capable of a positive first — not just `diff`. During review a working `astro dev` reported `HTTP=000` because it hit `nice(5) failed: operation not permitted` under the sandbox, which reads exactly like a page-missing failure. Prove the check can fire before believing it did not.

---

## 9. Questions or ambiguities

**Resolved by review:**

1. ~~**Whether to also assert the generated page in the `docs-build` job.** Taken: rely on the build failing, no extra assertion.~~ **`[OVERTURNED]`** The protection is real — verified empirically and in Starlight's source — but it is a **soft dependency on Starlight internals**. It fires only because the sidebar entry is `slug: 'reference/changelog'` (`astro.config.mjs:105`). Starlight does **not** validate `link:` targets, so a future `slug:`→`link:` switch, or removal of the sidebar entry, silently converts a hard build failure into a 404. One line in `ci.yml`'s existing `Assert build output` step removes the dependency for free — now §4 step 7. Reviewer A endorsed the original decision; Reviewer B identified the specific innocuous edit that breaks it. B's mechanism is concrete and the fix is one line, so it wins.

**Open (assumption taken, no blocker):**

2. **Whether contributors will trip over `astro build` directly.** Taken: document in `docs/README.md` (§4 step 6) rather than adding a guard. Strengthened by review — the failure is fail-loud (exit 1, no site), so the note is a legibility aid for an obscure error message, not a safety net.
3. **Custom loader instead of a script** (§5). Taken: script now, revisit if Starlight documents a composition point. Both reviewers endorsed.
4. **`[NEW]` Commit-pinned banner link instead of hardcoded `dev`.** Reviewer A proposed reading `GITHUB_SHA` — set in both workflows — and emitting `blob/<sha>/CHANGELOG.md`, falling back to `dev` when unset. The link would then point at exactly the bytes the page was generated from, and A6's unenforceable branch coupling would disappear rather than being documented. Same implementation cost. **Taken: hardcode `dev` for now**, because a SHA URL is less useful to a human who wants the current file and the coupling is one line in one place. Recorded because it is the better engineering answer if the coupling ever bites.
5. **`[NEW]` Whether to widen the H1 strip beyond line 1.** Reviewer A: the strip is anchored to line 1, so a line inserted *above* the H1 — a badge row, an HTML comment, a `<!-- markdownlint-disable -->` pragma — makes the strip no-op and ships two H1s with no signal. **Taken: do not widen the strip** (scanning for an H1 anywhere before the first `###` risks reaching content, which §3.0 exists to prevent); rely instead on §4 step 2a's zero-`^# ` assertion, which fails the build on exactly this case. Naming which mechanism covers it, per A's request.

**Critical:** none. Two live content bugs were found, but both are pre-existing and are split out to §3.4's sibling plan rather than blocking this one — subject to the ordering constraint stated there.

---

## 10. Self-Review

**Logic.** The plan hinges on one measured fact — the page is a pure lagging copy — and A1 says to re-verify it immediately before deleting the tracked file, because that is the step that cannot be undone by a rebuild.

**Adjusted while writing.** Two changes: dropped the idea of hand-syncing the page first (the generator's first run *is* the sync, and doing both creates a diff to reconcile); and moved from npm `pre*` hooks to explicit `&&` chaining, because a `prebuild` hook that silently does not fire is the same silent-success shape this session has been chasing all week.

**Adjusted after manual review.** Two substantive gaps, both found by checking the plan's own claims against the files:

- **The leading H1 (§3.0, V1, V1a, A7).** The plan described the generator as prepending frontmatter to a verbatim copy. That would have shipped two H1s, because Starlight renders `title` as the page `<h1>` and the root file opens with `# Trim Galore Changelog`. V1 as originally written — diff the stripped output against the raw root file — was *unsatisfiable alongside correct rendering*: passing it required keeping the duplicate. The 199-line drift figure had the same H1 folded into it; the real body gap is 196.
- **The banner's branch (§3.1, V1b, A6).** The link pointed at `blob/master/CHANGELOG.md` while the page will be generated from `dev`, which auto-deploys and is 34 commits ahead. Automation would have made this worse, not better — today's stale page sits near `master` by accident, masking it.

Both were in the plan's *verified* claims rather than its open questions, which is where the risk was: §2.1 and A1 read as settled and each carried an error. Worth remembering that "verified" in a plan means "the author checked", not "the check was falsifiable" — the A1 re-verification note now says to prove the diff can fail before trusting it.

**Adjusted after dual agent review.** Reports at `PLAN_review_reviewer-A.md` / `-B.md`. Both reviewers independently confirmed A1 with falsifiable harnesses and both went further than the plan: every diff hunk is a pure append, so the page body is a **strict subsequence** of the root body — "purely behind, never divergent" is exact, not an approximation. The design survived; the defects were in the supporting claims.

Both reviewers found, independently:

- **A4 is wrong** — `docs/README.md:39` also references the page, in a section that becomes false and that names a *second transformation* the plan had not checked (stale, verified). §4 step 6 promoted from conditional to mandatory.
- **§3.2 is wrong on both counts** — braces are not all in code spans (`:1145`, `:1185`), and `CHANGELOG.md:1164` is a raw `<img>` in prose. §3.2 restated around what actually holds, with a divergence table replacing the false GitHub-equivalence rule.
- **V9 could not fail** — `^### ` is already 49/49 on the stale page. Replaced.
- **§6 effect 5 understated its own failure** — a bare `astro build` exits 1 with no site, verified empirically by A and from `navigation.ts:145-159` by B.
- **The one-time checks should be generator assertions** — both called this the highest-value change. Now §4 step 2a.
- V3's page count, the 34→37 commit figure, "write only if content differs", and the missing CI-drift alternative (§5.1).

Unique to **Reviewer A**: the line-1-anchored strip hole (§9 Q5); V5's mis-specified expectation; **two leftover `199`s** that my own manual-review fold had missed; V6's narrow input coverage; the commit-pinned link (§9 Q4); the generalised harness note.

Unique to **Reviewer B**: **the `editUrl` breakage** (§4 step 4a) — the largest single miss in the plan, and one neither A nor I caught; **the KaTeX `$` bug** (§3.4), a live production defect on two published pages; the `slug:`-vs-`link:` soft dependency that overturned §9 Q1; `description`'s role in OG generation (A2); the generate-and-commit house-style deviation (§4 step 6); and that gitignoring the page breaks `ripgrep` discovery (§4 step 3).

They contradicted each other on one point of substance: **how to fix the broken `<img>`**. A left generator-side rewriting on the table as one of two options; B forecloses it, on the grounds that a broken-link repair is not a format change and so does not violate §1's non-goal, whereas generator rewriting destroys §3.0's "one transformation, cannot reach content" invariant that makes the generator auditable. **B's position adopted** — fixed at source in §3.4's sibling plan, generator stays verbatim.

**Not changed**, raised and consciously deferred: `npm run astro build` bypasses the generator (now folded into §6 effect 5); `astro dev` does not live-reload on `CHANGELOG.md` edits; the `logos` script is a precedent for style but is not itself wired into `build`; V9 overlaps V1 but both are kept since they now check different things.

**Open, not folded in.** `plans/docs-build-pr-check/PLAN.md:54` records that *"nothing under `docs/` reads `../CHANGELOG.md`, so `working-directory: docs` is safe"*. This plan retires that: the generator does read it, from a `working-directory: docs` context. Both reviewers confirmed it still holds in practice — both workflows use bare `actions/checkout` with no `sparse-checkout`, so `../CHANGELOG.md` is present, and `import.meta.url` resolution is cwd-independent — but the sibling plan's stated reasoning is stale, and no validation step here exercises the generator specifically under CI's cwd (V10 covers it only incidentally).

**Edge cases.** Missing/empty `CHANGELOG.md` (V6); a root file whose first line is *not* an H1 (§3.0 degrades to no-strip, A7); `astro build` invoked directly; `astro dev` (V5); a stale generated file left from a previous branch (the generator overwrites unconditionally); running the generator from the repo root vs `docs/` (resolve paths from `import.meta.url`); the `.md`-not-`.mdx` constraint (§3.2); pathspec checks run from inside `docs/` giving false answers (V7); a diff harness that cannot fail (A1).

**Remaining risks.**

- *Low:* a contributor confused by a page that exists locally but not in git. Mitigated by the `.gitignore` comment naming the generator, and the README note.
- *Low:* Starlight changes how collections resolve and the generated file stops being picked up. Caught by `docs-build` on the dependabot PR that bumps it — which is precisely the case #371 was built for.
- *Very low:* the root file gains a construct that renders on GitHub but not through remark. V4 checks the current content; a future entry could still break the docs build, which is now a PR-blocking check rather than a live-deploy surprise.
