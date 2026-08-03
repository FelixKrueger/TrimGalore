# PLAN review — Reviewer B

**Plan:** `plans/changelog-mirror-automation/PLAN.md`
**Reviewed:** 2026-08-01, repo root `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Method:** every load-bearing factual claim re-derived from the files. Comparisons run with real temp files under `$TMPDIR`, never process substitution; all `git` invocations from the repo root.

**Verdict:** the plan's central thesis is correct and I verified it independently — the page is a pure lagging copy and a one-way generate-and-overwrite is lossless. The design is sound. But three of the plan's *"Verified"* claims are wrong in ways that matter, one user-visible link breaks as a direct consequence of the change and is not mentioned anywhere in the plan, and one validation step (V9) cannot fail even on today's broken state.

---

## 0. Harness proof (so the numbers below can be trusted)

Before believing any comparison I proved the comparison could fail.

```
$ tail -n +10 docs/src/content/docs/reference/changelog.md > $TMPDIR/clrev/page_body.txt
$ tail -n +4  CHANGELOG.md                                > $TMPDIR/clrev/root_body.txt
    1132 page_body.txt
    1328 root_body.txt

# sentinel run — MUST report a difference
$ cat root_body.txt > root_sentinel.txt; printf 'ZZZ_SENTINEL_MUST_BE_REPORTED\n' >> root_sentinel.txt
$ diff page_body.txt root_sentinel.txt > sent.diff
diff exit=1
sentinel lines found in diff output: 1        # harness is falsifiable
```

Also noted for anyone re-running these checks: in this harness `grep` is a shell function wrapping **`ugrep --ignore-files`**, which honours `.gitignore`. Two consequences — (a) `grep -r` results here are not raw grep results, and (b) see finding O4 below.

---

## 1. Logic review

### 1.1 The core claim holds (verified)

| Claim (§2.1) | Verified? | Evidence |
|---|---|---|
| `CHANGELOG.md` 1331 lines, page 1141 | ✅ | `wc -l` |
| root last touched `7201f03` 2026-07-31; page `4f86c0f` 2026-07-25 | ✅ | `git log -1 --format='%h %ad' -- <path>` |
| strip page 1–9 / root 1–3 aligns the bodies | ✅ | page:10 and root:4 are both `### Unreleased` |
| **196 lines root-only, 0 page-only** | ✅ | 1328 − 1132 = 196; `diff` gives 196 `>` lines, 0 `<` lines, 201 total output lines = 196 + 5 hunk headers |
| no `###` release section missing | ✅ | 49 `^### ` in both files |

Stronger than the plan states: all five hunks are `a` (append) commands — `4a5,49`, `36a82,136`, `44a145,152`, `45a154,160`, `46a162,242`. There is **not a single `c` or `d` hunk**. "Purely behind, never divergent" is not an approximation, it is exact. A1 is safe and the §2.1 arithmetic reconciliation (196 body vs 199 raw) is right.

Unpublished-content spot checks all confirm §2.5: `#369` at `CHANGELOG.md:87`, absent from the page; `#363` 1 vs 0; `phred64` 9 vs 0.

### 1.2 §3.0's H1 problem is real, and I can show the mechanism

`docs/dist/reference/changelog/index.html` contains exactly one `<h1`, and it is `<h1 id="_top" class="…">Changelog`. Starlight supplies it from frontmatter `title`. `CHANGELOG.md:1` is `# Trim Galore Changelog`. So a verbatim body would ship two H1s exactly as §3.0 says, and V1a's expectation (0 in body, 1 in HTML) is correct.

A7 verified and slightly strengthened: there is exactly one ATX H1 (line 1), **and** no setext headings at all — `grep -nE '^=+[[:space:]]*$'` and `^-{3,}[[:space:]]*$'` both return nothing. Worth recording, because V1a's `grep -c '^# '` would be blind to a setext H1 (`Title` / `=====`); today none exists, so the check is adequate as written.

Heading census for context: 1 × H1, 0 × H2, 49 × H3, 40 × H4.

### 1.3 § 6 effect 5 contradicts §9 Q1 — and §9 Q1 is the correct one

§6 effect 5: *"A contributor running `astro build` directly (not via `npm run build`) gets no changelog page."*
§9 Q1: *"the build already fails if the page is missing (Starlight errors on the sidebar slug not resolving)."*

These cannot both be true. §9 Q1 is right. From `docs/node_modules/@astrojs/starlight/utils/navigation.ts:145-159`:

```ts
const route = routes.find((entry) => localizedSlug === entry.id);
if (!route) {
  …
  throw new AstroError(
    `The slug \`"${item.slug}"\` specified in the Starlight sidebar config does not exist.`,
    'Update the Starlight config to reference a valid entry slug in the docs content collection.\n' + …
  );
}
```

A bare `astro build` **hard-fails with a named, actionable error**, it does not silently omit a page. That is strictly better than §6 claims, and it means §4 step 6's README note is a convenience rather than a mitigation for a silent failure. §6 effect 5 needs rewording — leaving it as-is mis-calibrates the plan's own risk assessment.

Caveat that makes §9 Q1's decision less safe than it reads — see I2.

### 1.4 The CI/deploy wiring does what the plan needs (verified)

- `.github/workflows/docs.yml:5` → `branches: [master, dev]`; `:8` → `- 'CHANGELOG.md'`. Both §2.3 and A6 confirmed.
- `docs.yml` build job: `defaults.run.working-directory: docs`, step `Build site` → `run: npm run build`.
- `.github/workflows/ci.yml:158` `docs-build`, same `working-directory: docs`, same `npm run build`, and **no `paths:` filter** on `ci.yml`'s `pull_request` trigger (`ci.yml:6-10`) — so it genuinely runs on every PR to `master`/`dev`. §2.4's safety argument holds.

So wiring the generator into `npm run build` covers both workflows with one edit. Correct.

### 1.5 §4's mechanics check out

- `docs/package.json` scripts today: `dev`/`start` = `astro dev`, `build` = `astro build`, `astro` = `astro`, `logos` = `node scripts/generate-png-logos.mjs`. §4 step 2's proposed diff applies cleanly.
- `docs/scripts/generate-png-logos.mjs` uses `const __dirname = dirname(fileURLToPath(import.meta.url))` then `join(__dirname, '..', …)` — exactly the precedent §4 step 1 cites. Correct.
- `#### Infrastructure (contributor-facing)` exists at `CHANGELOG.md:226`, inside `### Unreleased` (lines 4–245). §4 step 5 correct.
- `docs/.gitignore` has no `src/content/…` entries yet; a path-with-slash pattern there is anchored to `docs/`, so `src/content/docs/reference/changelog.md` matches correctly. Correct.
- `docs/src/content.config.ts:6` matches §2.2's quote verbatim.
- `astro.config.mjs:105` → `{ label: 'Changelog', slug: 'reference/changelog' }`. §3 item 3 correct.
- `#368` = `c5c1834`, and its `docs/package.json` diff is `-"@astrojs/starlight": "^0.41.3"` / `+"^0.41.4"`. §5's dependency-churn argument is factually grounded.

---

## 2. Assumptions

### A1 — ✅ verified, and the plan's re-verification instruction is right

Confirmed above with a falsifiable harness. The instruction to prove the diff can fail before trusting it is the correct discipline and I followed it.

### A2 — ✅ verified

`docs/node_modules/@astrojs/starlight/schema.ts:15` → `title: z.string()` (no default, required). `description` is optional. Keeping `description` is not merely cosmetic: `docs/src/pages/og/[...route].ts:131-141` does `getCollection('docs')` and derives `description: entry.data.description || …`, so it feeds the per-page OG image. All 30 content entries produce 30 OG PNGs in `dist/og/`. A2's stated reason is correct.

### A3 — plausible, untested by me

`astro dev` reading the same collection is standard Astro behaviour and I did not start a dev server. The proposed `&&` chaining makes it moot either way. No objection.

### A4 — ❌ **the plan is wrong about a fact**

> *"Nothing else references the page's path except `astro.config.mjs:105` and inbound links. Grep before deleting."*

`git grep -n 'reference/changelog'` from the repo root returns **two** hits:

```
docs/README.md:39:`docs/src/content/docs/reference/changelog.md` is a copy of the
docs/astro.config.mjs:105:            { label: 'Changelog', slug: 'reference/changelog' },
```

`docs/README.md:37-41` is an entire section about this file:

```markdown
## Sync with the repo `CHANGELOG.md`

`docs/src/content/docs/reference/changelog.md` is a copy of the
top-level `CHANGELOG.md` with a Starlight frontmatter block prepended
and one inline image path rewritten. Keep them in sync when releasing.
```

The plan does gesture at this in §4 step 6, but as a conditional — *"if that file documents the build; check first"*. It does document it, in a dedicated section that becomes actively false. This is required work, not a check-and-maybe.

**And that section names a transformation the plan does not know about**: *"and one inline image path rewritten."* I verified this claim is **stale** — `CHANGELOG.md:1164` and page:974 are byte-identical (my diff found zero change hunks, which independently rules out any rewrite). So no second transformation exists today and §3.0's "this is the one transformation" survives. But the plan inherited that risk unknowingly; it should record that it checked and the README claim is obsolete.

Separately, `docs/README.md:18` says *"The site auto-deploys from the `master` branch"* — already false given `docs.yml:5` `[master, dev]` and A6. Same `master`-vs-`dev` confusion §3.1 diagnoses for the banner, sitting in the file the plan is about to edit.

### A5 — ❌ **falsified twice; this is my most significant finding**

> *"anything valid in the root file that renders on GitHub should render here"*

That rule is not safe, and I found two live counterexamples in the current file.

**(a) `$` is a math delimiter in this pipeline.**

`docs/astro.config.mjs:16-20`:

```js
markdown: {
  smartypants: false,
  remarkPlugins: [remarkMath],
  rehypePlugins: [rehypeKatex, [addClasses, { ".katex": "not-content" }]],
},
```

§3.2 checks `{` and `<`. It never checks `$`. `CHANGELOG.md:451`:

```
  default) — ~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour
```

This line is already published (page:261), and in the **already-built** `docs/dist/reference/changelog/index.html` it has been eaten by KaTeX:

```html
<annotation encoding="application/x-tex">7 vs ~</annotation>
… <span class="mord">7</span><span class="mord mathnormal">v</span><span class="mord mathnormal">s</span> …
```

`command grep -c '1000-sample cohort'` on that HTML returns **0** — the phrase is not contiguous in the output. `remarkMath` consumed `$7 vs ~$` as inline math, italicising "7vs" and swallowing both dollar signs. The page on www.trimgalore.com is wrong right now.

Same bug on a second page — `docs/dist/performance/benchmarks/index.html` contains
`<annotation encoding="application/x-tex">0.05/vCPU-hour, trimming 84M PE reads at the nf-core default costs roughly **`,
i.e. a whole sentence including its `**` bold markers swallowed into math. Notably `docs/src/content/docs/performance/benchmarks.md:129` already hand-escapes the tildes (`\~$0.05`, `\~$41`) — somebody has fought this before and only half-won. The changelog copy is the unescaped one.

KaTeX cannot simply be removed: `docs/src/content/docs/performance/clumpy.md:93-100` uses real display math, and `dist/performance/clumpy/index.html` shows `application/x-tex">\begin{aligned}`. But clumpy.md contains **only** `$$` fences (its only two `$`-bearing lines are the fences themselves). So the clean fix is to disable single-dollar inline math and keep `$$`:

```js
remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]],
```

Why this belongs in *this* plan rather than a separate one: today a human hand-copying the changelog is a (weak) reader of the rendered result. After automation, every future `$` in a changelog entry — a price, a `$var`, a `$cores` outside backticks — silently mangles with nobody in the loop. The plan's whole premise is that the page becomes an unattended artefact, which raises the cost of a silent transform bug. And no listed validation would catch it: V4 looks only for `{` and `<input.bam>`.

**(b) A repo-relative image path that works on GitHub and not on the site.**

`CHANGELOG.md:1164` (page:974, identical):

```html
<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="docs/Images/pigz_bench.png" >
```

- `docs/Images/` does not exist anywhere: `git ls-files | grep -i pigz_bench` → `docs/src/assets/screenshots/pigz_bench.png` only; `git ls-files | grep -i '/Images/'` → nothing.
- `dist/reference/changelog/index.html` emits `src="docs/Images/pigz_bench.png"` verbatim, so on `/reference/changelog/` it resolves to `/reference/changelog/docs/Images/pigz_bench.png`.
- `find docs/dist -iname '*pigz*'` → nothing. The asset is not in the build output.

So that image is broken on the live site (and, given no `docs/Images/`, on GitHub too). Fix it in `CHANGELOG.md` — a broken-link repair is not a "format change" and so does not violate §1's non-goal. Do **not** teach the generator to rewrite it; that reopens "the generator transforms content" and §3.0's one-transformation invariant is worth protecting.

### A6 — ✅ verified; one stale number

`docs.yml:5` = `branches: [master, dev]`. Confirmed. §3.1's argument that the banner must point at `dev` is correct and well-reasoned.

The count is stale: §3.1 says *"currently 34 commits ahead of `master`"*; `git rev-list --count master..dev` = **37** (and `dev..master` = 0). Direction and conclusion unaffected.

### A7 — ✅ verified

Exactly one H1, no setext headings, no `---`/`===` lines, no `:::` directives, no HTML comments in `CHANGELOG.md`. The conditional strip degrades safely.

### Implicit assumptions the plan does not state

1. **That the page carries no Starlight-supplied chrome tied to its being a tracked file.** It does — see C1 (`editUrl`).
2. **That `.md` raw-HTML passthrough stays enabled.** `CHANGELOG.md:1164` is raw HTML in prose; Astro allows it in `.md` by default. Worth folding into §3.2, because it is the real reason `.mdx` would be fatal: `<img … >` is not self-closed, so MDX would **fail the build**, not merely mis-render. That is a stronger argument than the one §3.2 makes and should replace it.
3. **That line endings are stable.** Verified benign: both files end in a single `0a`, `grep -c $'\r' CHANGELOG.md` = 0, and `.gitattributes` sets only `linguist-*` attributes (no `text`/`eol`). No CRLF hazard for V1/V2.

---

## 3. Efficiency analysis

Non-issue, and I can put numbers on it. `CHANGELOG.md` is **85,525 bytes / 1331 lines**. The generator is one `readFile`, one regex-ish head strip, one `writeFile` — single-digit milliseconds. The surrounding `astro build` renders 30 content entries into 32 HTML pages and rasterises 30 satori OG PNGs (`dist/og/`, all >10 KB per CI's runt check). The generator is three orders of magnitude below the noise floor. No streaming, no caching, no scalability concern. Nothing to optimise, and the plan is right not to discuss it.

One micro-improvement with a real (small) benefit: **write only if the content differs.** Unconditional overwrite bumps `mtime` on every `npm run dev`, and `astro dev`'s watcher keys on file events. A read-compare-write makes a no-op run a genuine no-op and hardens V2 by construction. Optional.

---

## 4. Validation sufficiency

The plan nominates V6, V7, V1a as the three that would quietly not hold. I agree with all three and I verified V1a's expected value against `dist` (exactly one `<h1 id="_top">Changelog`). But there are gaps.

### V9 cannot fail — it already passes on the broken state

> *V9: Count `###` sections in the generated page vs the root file … Expected: Equal*

`^### ` count is **49 in both files today**, before any change. V9 is satisfied by the stale page. Worse, it is insensitive to exactly the regression it is meant to catch: the count that actually differs is **`^#### ` — 40 in root vs 39 in the page** (the missing `#### Infrastructure (contributor-facing)`). A generator that dropped an entire `####` subsection would pass V9.

Fix: assert post-strip **body line-count equality** and **`^#### ` count equality** alongside `^### `. Or drop V9 — V1 already subsumes it (the plan notes the overlap in §10 but keeps the weaker check).

This also weakens §2.1's evidence: "no `###` release section is missing" is true but was never a discriminating test.

### Nothing validates the two A5 failure modes

- No check for `$`/KaTeX. Add: `command grep -c 'application/x-tex' dist/reference/changelog/index.html` → expect **0** (after the `singleDollarTextMath: false` fix). This is a one-line assertion for a defect class that builds, deploys, and renders wrong.
- No check that the `<img>` resolves. Add an assertion that no `src="docs/` survives in the built page, or that the referenced asset exists in `dist/`.

V4 as written (*"no raw `{` or stray `<input.bam>` markup leaked"*) would catch neither, and both are present in the current content.

### V3's threshold encodes today's page count

V3 expects "32+ HTML pages". `find docs/dist -name '*.html' | wc -l` = **exactly 32**. So the threshold has zero headroom and, more importantly, it is not specific — 32 pages can be met while the changelog page is absent. Replace with the thing you actually care about:

```bash
test -f dist/reference/changelog/index.html
```

§10 already notes V3 is stricter than CI's `-ge 1`; the answer is not "looser" or "stricter" but "targeted".

### V1b, V6, V7 are well specified

V1b's expectation (`blob/dev/CHANGELOG.md`) is right — the current banner is `docs/src/content/docs/reference/changelog.md:7` linking `blob/master/CHANGELOG.md`. V7's warning about pathspecs from inside `docs/` is real; I hit the `docs`/`Docs` resolution oddity twice during this review (a `Docs/.astro/data-store.json` hit, and a relative-path `grep` reporting "No such file or directory" for a file that exists) — absolute paths cured it both times.

Also worth noting for V1: the page's existing frontmatter (`docs/src/content/docs/reference/changelog.md:1-4`) already has **exactly** the `title` and `description` §3.1 proposes, character for character. So the frontmatter is a no-op change and V1's stripping logic can be written against the current file with confidence.

### Missing: no validation of the deleted file's replacement chrome

Nothing checks the rendered page's Starlight-supplied links after the file stops existing in git. See C1.

---

## 5. Alternatives

**§5's rejection of a custom content loader — I would make the same call.** The reasoning is sound and better-grounded than the plan claims: `docsLoader()` is imported from `@astrojs/starlight/loaders` (`content.config.ts:2`), a subpath export whose return shape is not a documented contract, and the churn argument is concrete (0.41.3 → 0.41.4 inside six days, `c5c1834`). Agree, no change.

**Alternative A — symlink the root file into the collection.** Dead on arrival, and §2.2 is right about why: frontmatter is mandatory (`schema.ts:15`), so the file in the collection cannot be the root file. Not worth mentioning in the plan.

**Alternative B — keep the page tracked, add a CI drift check** (`npm run sync:changelog && git diff --exit-code -- docs/src/content/docs/reference/changelog.md`). This deserves an explicit mention-and-reject in §5, because it is the one alternative that **dodges C1 and O4 entirely**: the file stays in git, so the Starlight edit link keeps working and the page stays greppable. Its cost is exactly what §1 argues against — drift becomes *detectable* rather than *unrepresentable*, and a human has to run the generator and commit. I agree with the plan's preference, but §5 currently reads as if the loader were the only alternative considered, and this one is the closer call.

**Alternative C — a plain `src/pages/` route outside the docs collection.** Sidesteps the schema entirely but loses sidebar integration and Starlight page chrome, and would need the sidebar entry switched from `slug:` to `link:` — which, per I2, also throws away the build-time slug validation the plan is relying on. Worse on both counts. Reject.

**Alternative D — let the generator normalise the known GitHub-vs-Astro divergences** (escape lone `$`, rewrite the image `src`). Tempting given A5, but reject: it destroys §3.0's "one transformation, anchored to the top of the file, cannot reach content" invariant, which is the property that makes the generator auditable. Fix both at source in `CHANGELOG.md` and in `astro.config.mjs` instead.

---

## 6. Action items

### Critical

**C1. The generated page will carry a broken "Edit page" link. Not mentioned anywhere in the plan.**
`docs/astro.config.mjs:36-39` sets `editLink.baseUrl: 'https://github.com/FelixKrueger/TrimGalore/edit/master/docs/'`, and the current built page contains:

```html
href="https://github.com/FelixKrueger/TrimGalore/edit/master/docs/src/content/docs/reference/changelog.md"
```

Once the file is untracked and deleted, that path exists on no branch. GitHub's `/edit/` route on a missing path does not 404 for an authenticated user — it opens the **new-file** editor at that path, so a contributor clicking "Edit page" would be invited to recreate the tracked copy the plan just deleted, resurrecting the exact drift this plan eliminates. (I could not confirm the post-auth behaviour anonymously: `/edit/…` redirects to `/login` for both existing and nonexistent paths, so HTTP status does not discriminate. The fix is correct under either behaviour.)

Fix — `editUrl` is a supported frontmatter field, `docs/node_modules/@astrojs/starlight/schema.ts:30`:

```ts
editUrl: z.union([z.url(), z.boolean()]).optional().default(true),
```

Emit it from the generator, pointing at the real source and at the deploying branch (same reasoning as §3.1's banner, which the plan already got right for the banner and missed here):

```yaml
editUrl: https://github.com/FelixKrueger/TrimGalore/edit/dev/CHANGELOG.md
```

`editUrl: false` is the acceptable minimum. Add a validation step asserting the built page's edit link is not the old collection path.

**C2. §3.2 / A5 miss `$`, and the current content already mis-renders because of it.** Detail and evidence in A5(a). Concretely: add `$` to §3.2's checked constructs; set `remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]]` in `docs/astro.config.mjs:18` (safe — `clumpy.md:93-100` uses only `$$`); and add a validation asserting `application/x-tex` does not appear in `dist/reference/changelog/index.html`. Also weaken A5's rule from "renders on GitHub ⇒ renders here" to an explicit list of pipeline-specific hazards (`$…$` math, raw-HTML relative paths, MDX-if-ever-converted), because the general rule is demonstrably false.

**C3. V9 is a check that passes on the broken state.** `^### ` is already 49 = 49; the discriminating count is `^#### ` (40 vs 39). Replace V9 with post-strip body-line-count equality plus `^#### ` equality, or fold it into V1. As written it provides no signal and would not notice a dropped `####` subsection. See §4.

### Important

**I1. `docs/README.md` needs a rewrite, and A4 needs correcting.** A4's "nothing else references the page's path" is false — `docs/README.md:39`. Three edits, not one note:
- `:37-41` — the whole `## Sync with the repo CHANGELOG.md` section becomes false. Rewrite as "generated at build time by `scripts/sync-changelog.mjs`; do not edit, do not commit".
- `:41` — *"and one inline image path rewritten"* is stale (verified: `CHANGELOG.md:1164` == page:974 byte-identical). Delete it, and record in the plan that this claimed second transformation was checked and does not exist.
- `:18` — *"auto-deploys from the `master` branch"* is already wrong per `docs.yml:5` and A6. Fix while you are in the file.

**I2. Fix §6 effect 5, and overturn §9 Q1.** §6 effect 5 says a bare `astro build` yields "no changelog page"; `navigation.ts:154-158` shows it throws an `AstroError` naming the slug. Correct the plan text.

Then reconsider §9 Q1's "rely on that, no extra assertion". That protection is a **soft dependency on Starlight internals**: it fires only because the sidebar entry is `slug: 'reference/changelog'` (`astro.config.mjs:105`). Starlight does not validate `link:` targets, so a future switch from `slug:` to `link:` — an innocuous-looking edit — silently converts a hard build failure into a 404. One line in `ci.yml`'s existing `Assert build output` step removes the dependency for free:

```bash
test -f dist/reference/changelog/index.html
```

**I3. Fix the broken image at source.** `CHANGELOG.md:1164` points at `docs/Images/pigz_bench.png`, which exists nowhere; the real asset is `docs/src/assets/screenshots/pigz_bench.png`, and no PNG lands in `dist`. Repoint it at a URL that works in both renderers (a `raw.githubusercontent.com` absolute URL is the simplest). Evidence in A5(b). Keep the generator verbatim.

**I4. Replace §3.2's evidence — it is wrong on both counts, and the truth is a stronger argument.**
- *"17 lines contain `{`, all inside backtick code spans"* — count is right, "all inside code spans" is false. `CHANGELOG.md:1145` (`… ending in .{INT}bp_3prime.fq(.gz).`) and `:1185` (`e.g. -a "A{10}" to trim poly-A tails`) are bare prose. Harmless in `.md` — verified in `dist`, both render literally at page:955/995 — but the stated reason does not hold.
- *"22 lines contain `<`-prefixed text … all inside fenced blocks or code spans"* — the count is **24**, and `:1164` is a raw HTML `<img … >` tag in prose, in neither a fence nor a code span.
- Restate as: braces and angle brackets are safe in `.md` **regardless** of code-span context, empirically proven by the already-published lines; and `.mdx` would be fatal rather than merely wrong, because `<img … >` is not self-closed and MDX would fail the parse. That argument survives A5's counterexamples; the current one does not.

**I5. Strengthen V3 and V4.** Replace V3's "32+ HTML pages" (today's count is exactly 32, and the number is not specific to the page under test) with `test -f dist/reference/changelog/index.html`. Extend V4 beyond `{`/`<input.bam>` to assert no `application/x-tex` and no surviving `src="docs/` in the built page.

### Optional

**O1.** §3.1's "34 commits ahead of `master`" is now 37. Prefer "ahead of `master`" without a count — a volatile number in a plan invites exactly the "verified means the author checked once" problem §10 diagnoses.

**O2.** Note the departure from the repo's existing generated-artefact convention. `docs/public/logos/*.png` are **generated by `npm run logos` and committed** (`git ls-files docs/public/logos/` lists the PNGs). So the house style for generated files is generate-and-commit; this plan is generate-and-ignore. The plan's choice is right here — committing a generated artefact is the drift being eliminated — but say so, or the next contributor "fixes" the inconsistency in the wrong direction. One sentence in §4 step 3 next to the `.gitignore` comment.

**O3.** Make the generator write only when the content differs, so re-runs are true no-ops (no `mtime` churn for `astro dev`'s watcher) and V2 holds by construction rather than by observation.

**O4.** Note the search-tooling consequence of gitignoring the page. `ripgrep` and this harness's `grep` (a wrapper around `ugrep --ignore-files`) both honour `.gitignore`, so after this change `rg 'some changelog phrase' docs/` finds nothing. Contributors will grep for changelog text and conclude the page does not exist. Worth a clause in the `.gitignore` comment: "generated by scripts/sync-changelog.mjs — grep CHANGELOG.md in the repo root instead".

**O5.** Add the CI-drift-check alternative (Alternative B above) to §5's mention-and-reject list. It is the nearest rival to this design and the one that would have avoided C1 — recording why it lost makes the choice read as deliberate, which is §5's stated purpose.

---

## 7. What I did not verify

- **A3** (`astro dev` reads the same collection). I did not start a dev server. Standard Astro behaviour; the `&&` chaining makes it moot.
- **Post-authentication GitHub `/edit/` behaviour on a missing path** (C1). Anonymous requests redirect to `/login` for existing and nonexistent paths alike, so status codes do not discriminate. The recommended fix is correct either way.
- **A live build with the generator in place.** My local `docs/node_modules` holds Starlight **0.41.3** while `docs/package.json` declares `^0.41.4`, so my `schema.ts` / `navigation.ts` reads are against 0.41.3. Both `editUrl` and the sidebar-slug `AstroError` are long-standing, so the conclusions carry; flagging it for completeness.
- **The sibling plan quote** at `plans/docs-build-pr-check/PLAN.md:54` (§10's "Open, not folded in"). Not read, per the standing instruction on `plans/`. The plan's own account of it — that the reasoning is now stale but the conclusion still holds via `actions/checkout` plus `import.meta.url` — is consistent with what I verified about `working-directory: docs` in both workflows.
