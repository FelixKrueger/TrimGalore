# PLAN — Two docs-render bugs: KaTeX eating prose, and a broken image path

**Branch:** `docs/fix-katex-and-image-path` off `dev`
**Problem:** two constructs in published content render wrongly on www.trimgalore.com today. Both were found by the dual agent review of `plans/changelog-mirror-automation/PLAN.md`; neither is caused by that plan, and both get worse once it lands.

**Revision status:** `[REVISED AFTER DUAL AGENT REVIEW]`. Both reviewers returned *do not implement as written* — not because the diagnosis was wrong, but because one instructed step introduced a new bug and three of eight validations could not fail. Both defects are corrected below; see §9 for the full account.

---

## 1. Goal

Fix two live rendering defects at source:

1. `$…$` in prose is consumed as inline math by `remarkMath`, silently corrupting prose on two published pages — one of which renders a full sentence in **red error text**.
2. `CHANGELOG.md:1164` carries a raw HTML `<img>` whose `src` points at a path that exists nowhere.

Non-goal: changing how the changelog page is produced (that is the sibling plan), removing KaTeX, or reformatting `CHANGELOG.md`.

---

## 2. Context

### 2.1 Why this is a separate plan

The sibling plan (`plans/changelog-mirror-automation/`) makes the published changelog a build artefact of `CHANGELOG.md`. These bugs live in `astro.config.mjs`, `CHANGELOG.md` and one asset path — none in the generator, `package.json`, or the content collection — so the diffs do not overlap.

Split out for two reasons. Scope: a docs-automation change should not carry a Markdown-pipeline configuration change. And **ordering**: today a human hand-copying the changelog is a weak reader of the rendered result. Automation removes that human, so every future `$` in a changelog entry mangles unattended. Fixing the silent-corruption path *before* automating is the correct sequence. See the sibling plan's §3.4.

### 2.2 Bug 1 — `$` is a math delimiter in this pipeline

`docs/astro.config.mjs:16-20`:

```js
markdown: {
  smartypants: false,
  remarkPlugins: [remarkMath],
  rehypePlugins: [rehypeKatex, [addClasses, { ".katex": "not-content" }]],
},
```

`remarkMath` defaults to enabling single-dollar inline math. There are **three** affected prose sites, and the damage differs at each. Verified by AST sweep and against built output by both reviewers.

**`CHANGELOG.md:451` / `changelog.md:261`** — `~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. `$7 vs ~$` is consumed as math. The rendered line becomes:

```
default) — ~7vs 7 vs ~7vs 41 per 1000-sample cohort at AWS $0.05/vCPU-hour
```

Two `$` deleted, one `~` deleted, `7vs` italicised as math variables, and content duplicated by KaTeX's MathML+HTML double emission.

`[CORRECTED]` An earlier version of this plan claimed `grep -c '1000-sample cohort'` returns **0** on that HTML, "the phrase is no longer contiguous". **It returns 1.** Both reviewers measured it independently on genuine pre-fix builds. The math node is `<annotation encoding="application/x-tex">7 vs ~</annotation>`, which ends *before* that phrase. The claim was inherited from an earlier review without verification, and it made V4 a check that could not fail. The pattern that *does* discriminate is the full dollar-bearing phrase — see V4.

**`benchmarks.md:129`** — **two** `inlineMath` nodes, not one, and the second is worse:

```html
<span class="katex-error"
      title="ParseError: KaTeX parse error: Expected group as argument to '\~' at end of input: … scales to **\~"
      style="color:#cc0000">0.007 with the Rust v2 build** — a 5.9× saving per sample.
      Across a 1000-sample cohort that scales to **\~</span>
```

A full sentence rendered in `#cc0000`, a literal `\~` leaking through, orphaned `**`, and `<strong>` spanning the wrong range. This is the more visible of the two defects and the original plan did not mention it.

**KaTeX cannot simply be removed.** `performance/clumpy.md:93-100` contains real display math (`dist/performance/clumpy/index.html` shows `application/x-tex">\begin{aligned}`), and its only `$`-bearing lines are its `$$` fences. Disabling single-dollar math keeps every legitimate use.

### 2.3 Bug 2 — an image path that resolves nowhere

`CHANGELOG.md:1164` (and, identically, `changelog.md:974`):

```html
<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="docs/Images/pigz_bench.png" >
```

- `docs/Images/` **does not exist**. `git ls-files | grep -i '/Images/'` returns nothing.
- The real asset is `docs/src/assets/screenshots/pigz_bench.png` (PNG, 1810×541, 80,712 bytes). Both reviewers opened it: a two-panel bar chart, *"SE trimming with Python3 + pigz (~17.5M reads, 2x75bp)"* and the PE equivalent, runtime versus `# of cores`. It is unambiguously the image the v0.6.0 `-j/--cores` entry at `CHANGELOG.md:1162` refers to.
- **The asset is an orphan.** `git grep pigz_bench` returns only the two broken lines, and it is absent from the build's 22 optimised images — Astro only emits `src/assets` files that its image pipeline references. So it can be moved freely.
- From `/reference/changelog/` under `base: '/'` the relative `src` resolves to `https://www.trimgalore.com/reference/changelog/docs/Images/pigz_bench.png` → 404. Broken on GitHub too, for the same reason.

### 2.4 Why fix at source rather than in the generator

The changelog reviewers disagreed here; B's position is adopted. A broken-link repair in `CHANGELOG.md` is not a format change, so it does not violate the sibling plan's non-goal, whereas teaching the generator to rewrite content destroys its "one transformation, anchored to the top of the file, cannot reach content" invariant.

Strengthened by measurement: the `<img>` at `:1164` is the **only** raw-HTML construct outside code in the entire file, and there are **zero** relative link or image targets. A generator rewrite rule would exist to serve exactly one line forever.

### 2.5 `[NEW]` A third hazard class: GFM single-tilde strikethrough

`[ADDED BY REVIEW — this is what nearly went wrong]`

`micromark-extension-gfm-strikethrough` defaults `singleTilde` to **true**, and Astro enables gfm by default (`@astrojs/markdown-remark/dist/index.js:41,48`), which `astro.config.mjs` does not disable. So a bare `~` can open strikethrough.

Bare `~` appears on roughly **50 prose lines** across the content and `CHANGELOG.md` (`~75% larger`, `~27 threads`, `~552 MiB`, …). Today **zero** produce `<del>` — verified across all 32 built pages — but that holds by luck, not design. It depends on attention-style flanking: in `**~$41` the `~` is preceded by `*`, non-whitespace, so it is a valid right-flanking closer and pairs with the opener in `at ~$0.05`. In `— ~$7 vs ~$41` the second `~` is space-preceded and cannot close.

This is why the `\~` escapes at `benchmarks.md:129` exist, and why removing them is a regression rather than a cleanup — see §3.3.

---

## 3. Behavior

### 3.1 KaTeX

Disable single-dollar inline math, keep `$$` display math:

```js
remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]],
```

Verified at four layers by the reviewers, so A1's hedge is dropped:

- `remark-math@6.0.0` is the only copy installed and forwards options verbatim (`lib/index.js` → `micromarkExtensions.push(math(settings))`).
- `micromark-extension-math/lib/math-text.js` reads `singleDollarTextMath` (default `true`) and enforces it: `if (sizeOpen < 2 && !single) return nok(code)`.
- Astro honours the tuple form (`@astrojs/markdown-remark/dist/load-plugins.js:12-14`, applied at `dist/index.js:56-58`).
- **`$$` display math is structurally immune**, not merely observed to survive: `micromark-extension-math/index.js` passes options only to `mathText`, never to `mathFlow`.

Real builds by both reviewers: changelog and benchmarks drop to zero math *and* zero error nodes; clumpy retains its one display node; 32 HTML pages; build green.

### 3.2 The image — `[REVISED]` self-host it

`[REVISED BY REVIEW]` The original plan proposed a `raw.githubusercontent.com` URL pinned to `dev`, and both reviewers rejected the branch choice: the two existing in-repo precedents (`README.md:3-4`, `astro.config.mjs:38`) both pin `master`, and the asset serves 200 on both branches so `dev` gained nothing.

Rather than answer the branch question, **dissolve it.** Move the orphaned asset into `docs/public/` and use a site-absolute URL:

```
https://www.trimgalore.com/images/pigz_bench.png
```

Why this over a `master`-pinned raw URL:

- No branch pin at all, so "the asset moves and the image silently breaks" cannot happen.
- Self-hosted: no third-party dependency, no `cache-control: max-age=300`, no raw.githubusercontent rate limiting, and it avoids introducing the site's *only* cross-origin asset.
- One line, serving both renderers.
- The pattern is already established — `docs/public/logos/` exists and is what the README's raw URLs point at.
- Free: the asset is an orphan (§2.3), so moving it breaks nothing.

Cost: a brief window between merge and deploy where the URL 404s. Accepted — the image is already broken, so this is not a regression.

**Also add an `alt` attribute** while the tag is being edited. Both reviewers raised it: the tag has only `title`, and the chart carries real information. The malformed `style="float:right;margin:20px 20 20 600px"` (three unitless values) and the space-bearing `id` are deliberately left alone — pre-existing, harmless, out of scope.

### 3.3 `[REVERSED]` The `\~` escapes are load-bearing and must stay

The original plan said the `\~$0.05` / `\~$41` escapes at `benchmarks.md:129` "become unnecessary once single-dollar math is off", that leaving them was "harmless — but they are now misleading", and instructed removing them.

**That premise was wrong, and acting on it would have shipped a new bug.** The escapes guard GFM single-tilde strikethrough (§2.5), not `$` math, and are unaffected by the `remarkMath` change. Both reviewers measured it independently through Astro's exact plugin chain:

```
WITH escapes (current source):
  <del> count = 0
  … <strong>~$41 with TG vs ~$7 with Rust v2</strong>, with proportional savings

WITHOUT escapes (what the original plan instructed):
  <del> count = 1        <<< strikethrough introduced
  On AWS at <del>$0.05/vCPU-hour, … scales to **</del>$41 with TG vs ~$7 with Rust v2**
```

So the change would have traded a KaTeX corruption for a strikethrough corruption of precisely the same character — a sentence swallowed into a wrapper element and a `**` bold pair collapsed to literal asterisks — in the commit that claims to fix silent prose corruption.

**And every listed validation passed on it.** `application/x-tex` 0, `katex-error` 0, V4 and V7 were no-ops, and V5 as written was satisfied by the two surviving `<strong>` elements. Hence V9 (§8), which is the only check that catches this class.

**Decision:** keep the escapes. Add a one-line source comment stating what they guard — `guards GFM single-tilde strikethrough, not $ math` — so the next reader does not repeat this mistake.

Provenance, for the record: `git log -L 129,129` shows they were added in `3af218b` (#282), alongside the "Oxidized" → "the Rust v2 build" rename.

### 3.4 `[NEW]` Why not fix strikethrough at the config layer too

§2.2 argues that hand-escaping is the wrong layer because it requires every future author to know. That argument applies verbatim to `\~`, so the symmetric fix would be `remark-gfm`'s `singleTilde: false`. **Checked and rejected:** Astro exposes `markdown.gfm` as a **boolean only**, and it is deprecated in this version (`astro/dist/types/public/config.d.ts:2267-2283`). Getting sub-options means `gfm: false` plus manually re-adding `remark-gfm`, which risks tables (used on `benchmarks.md` itself), autolinks and Starlight's expectations — to defend one line.

Recorded as a deliberate asymmetry rather than an oversight: `$` is fixed at the config layer, `~` is fixed by keeping an escape and making the hazard *visible* via V9.

---

## 4. Implementation outline

1. **`docs/astro.config.mjs:18`** — `remarkPlugins: [remarkMath]` → `remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]]`.
2. **Move the asset and repoint the `<img>`.** `git mv docs/src/assets/screenshots/pigz_bench.png docs/public/images/pigz_bench.png`, then update `CHANGELOG.md:1164`'s `src` to `https://www.trimgalore.com/images/pigz_bench.png` and add an `alt`.
3. **`docs/src/content/docs/performance/benchmarks.md:129`** — **leave the `\~` escapes in place**; add the one-line comment from §3.3. *(This step formerly instructed removing them. See §9.)*
4. **`docs/src/content/docs/reference/changelog.md`** — **one line, not two.** `[CORRECTED]` The math defect at line 261 is fixed by step 1 alone, since that page goes through the same pipeline as every other. Only line 974's `<img src>` needs a hand edit. Disappears entirely if the sibling automation lands first — see §7 Q1.
5. **Add an `Assert no silent render corruption` step to *both* workflows.** `[REVISED AFTER CODE REVIEW]` This is the only change that prevents *recurrence* rather than fixing one instance — but the first version of this step did not, for two reasons both code reviewers found independently:

   - **`katex-error` is blind to the defect this plan fixes.** `rehype-katex` emits an error span only when KaTeX *fails*; a successful accidental parse emits `class="katex"` and no error span. The pre-fix changelog page renders 3 `class="katex"` and **0** `katex-error`, so the original three-pattern step exited **0** on a tree whose only defect was the changelog corruption. §8 had already split V1 (`class="katex"`) from V2 (`katex-error`) for exactly this reason; only V2's pattern was carried into CI. **V1 needs a CI counterpart**, which has to be an allowlist because `clumpy.md`'s `$$` maths is legitimate.
   - **`ci.yml` gates pull requests; `docs.yml` publishes the site.** `docs.yml` deploys on every push touching `docs/**` or `CHANGELOG.md`, has no assertions, and no `needs:`/`workflow_run` gate on `ci.yml`. Both defects fixed here reached production by that path, so the step must exist in both — which is also what `docs.yml:22` ("keep the steps in step") already asked for.

   The step, identical in both files (`ci.yml` after `Assert build output`; `docs.yml` between `Build site` and `Upload Pages artifact`):

   ```bash
   test -d dist                       # grep exits 2 on a missing dir, which `if` swallows
   fail=0
   for pat in 'katex-error' '<img[^>]*src="[Dd]ocs/' '<del>'; do
     if grep -rlE --include='*.html' -- "$pat" dist; then
       echo "::error::'$pat' in built output (files above)"
       fail=1
     fi
   done
   stray=$(grep -rlF --include='*.html' -- 'class="katex' dist \
           | grep -v '^dist/performance/clumpy/index.html$' || true)
   if [ -n "$stray" ]; then fail=1; fi
   exit $fail
   ```

   Three hardening choices, each from a specific review finding: `fail=0` accumulates rather than fail-fast, so co-occurring defects surface in one run; the image pattern is anchored to `<img[^>]*src="[Dd]ocs/` rather than the bare string, because double quotes survive unescaped inside `<code>` so prose *documenting* the bug would otherwise trip the guard — the same own-goal as iteration #1, one layer out; and `[Dd]` covers this filesystem's case-insensitivity.

   The `grep -v` allowlist is a maintenance point: a second legitimate maths page must be added to it. That is the right trade against a sentinel that cannot see the defect.
6. **CHANGELOG entry** under `#### Bug fixes` in `Unreleased`. User-visible docs defects, so *not* the `#### Infrastructure (contributor-facing)` heading.

---

## 5. Integration

**Reads/writes:** `docs/astro.config.mjs`, `CHANGELOG.md`, `docs/public/images/pigz_bench.png` (moved), `docs/src/content/docs/performance/benchmarks.md` (comment only), `.github/workflows/ci.yml`, and conditionally `changelog.md` (§4 step 4).

**Visible effects:**

1. `~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour` renders in full on the changelog page.
2. **Both** swallowed passages on `/performance/benchmarks/` reappear — including the one currently rendered in red error text. `[CORRECTED — was singular]`
3. The multi-threading benchmark image loads on `/reference/changelog/`.
4. `$$` display math on `/performance/clumpy/` is unaffected.
5. Three new CI assertions make all three defect classes non-recurring.

**Interacts with:** `plans/changelog-mirror-automation/` — see §2.1 and §7 Q1. No file overlap.

---

## 6. Assumptions

- **A1.** `remark-math` supports `singleDollarTextMath`. **Verified** at four layers (§3.1), including that `$$` is structurally immune. Hedge dropped. *Note:* `// @ts-check` at `astro.config.mjs:1` gives **no** protection against a typo'd option name — Astro types plugin options as `any`, and there is no `astro check` or `tsc` anywhere in the repo or either workflow. That is the concrete reason V1 must assert against built output, not source.
- **A2.** Nothing relies on single-dollar inline math. **Verified exhaustively** — both reviewers parsed all 30 content entries plus `CHANGELOG.md` with the real `remark-parse` + `remark-math` and enumerated every math node under both settings. Only `benchmarks.md` (2 nodes), `changelog.md` (1) and `CHANGELOG.md` (1) change; the sole survivor is `clumpy.md`'s `$$` display node. Every other entry has zero `$`. *One latent site:* `index.mdx:21` has a lone `$` in `saves time, $, CO₂` — inert today because a single `$` cannot pair, but a live tripwire under the current default. MDX inherits the config (`@astrojs/mdx/dist/index.js:106` sets `extendMarkdownConfig: true`), so it is covered by the fix. This is an argument *for* the fix that the original plan missed.
- **A3.** `docs/src/assets/screenshots/pigz_bench.png` is the intended image. **Verified by inspection**, not inferred — see §2.3. Also an orphan, which is what makes §3.2's move free.
- **A4.** A site-absolute `https://www.trimgalore.com/images/…` URL works in both renderers. **No CSP exists anywhere** to interfere: no `_headers`, `netlify.toml`, `vercel.json` or `staticwebapp.config.json`; no CSP `<meta>` in `docs/src` (including `Head.astro`), `astro.config.mjs` or `docs/public`; GitHub Pages injects none. *Correction to the original:* the plan cited the README as precedent for `raw.githubusercontent` — true, but the README pins **`master`**, which undercut the `dev` choice rather than supporting it. §3.2 now avoids the question.
- **A5.** These are the only pipeline-divergence bugs in the current content. **Substantially narrowed, not closed.** Verified clean by both reviewers: the math class (fully enumerated), the raw-HTML class in `CHANGELOG.md` (exactly one node — the three other angle-bracket candidates are inside *multi-line* inline code spans that a per-line grep cannot see, and all render correctly), relative link/image targets (zero), and broken `src` attributes site-wide (the one known). **Not swept:** emoji shortcodes, footnotes, bare-URL autolinking. `{`-brace/JSX hazards cannot bite the changelog page — it is `.md`, not `.mdx`.
- **A6.** `[NEW]` GFM single-tilde strikethrough is a fourth divergence class, live in the source today, currently producing zero `<del>` by flanking luck rather than design (§2.5). It is the class most likely to bite the sibling automation, and the one the original audit never considered.
- **A7.** `[NEW]` `singleDollarTextMath: false` requires **two or more** dollars to open text math — it does not remove the construct. `$$` in prose still becomes inline math after the fix. No such case exists today, but `$$` is the shell PID idiom and the sibling plan pipes `CHANGELOG.md` in unattended, so §2.1's claim that this fix "closes the silent-corruption path" is narrower than it sounds.

---

## 7. Questions or ambiguities

**Critical — needs a decision before implementing:**

1. **Which plan lands first?** The sibling plan's §3.4 argues this one should, so automation does not ship a silent-corruption path unattended. **Recommendation: this plan first.** Its cost is one throwaway line in the tracked changelog page (§4 step 4) — `[CORRECTED]` one line, not two, since the math defect there is fixed by the config change alone. That strengthens the recommendation.

**Resolved by review:**

2. ~~Branch in the image URL (`dev` vs `master`).~~ **Dissolved.** Both reviewers rejected `dev`; rather than switch to `master`, §3.2 self-hosts from `docs/public/`, which removes branch pinning as a concept. Reviewer A recommended `master`, Reviewer B leaned toward this option; taken as B's, and reversible to a `master`-pinned raw URL if self-hosting proves awkward.
3. **Whether to move the image into the docs site properly.** The original rejection addressed only the Astro-import member of that family (`![](../../assets/…)`, which GitHub cannot resolve — correct, and `remarkCollectImages` running after user plugins confirms it *would* be optimised). `[REVISED]` The `docs/public/` member is now the chosen approach; §3.2 records why the whole family no longer loses.

---

## 8. Validation

`[HEAVILY REVISED]` Three of the original eight checks could not fail, and none caught the §3.3 regression. Method notes first, because two of them determine whether any number below is trustworthy.

**Capture baselines to a scratch directory, never `docs/dist/`.** It is gitignored shared mutable state and it moved three times during review — at one point holding a *post-fix* build against *pre-fix* source. Use `npx astro build --outDir "$TMPDIR/before"` (verified working, leaves `docs/dist` untouched) and diff two immutable trees.

**Count occurrences, not lines.** These pages put whole documents on very few lines — both math nodes on `benchmarks` are on line 439 — so `grep -c` under-reports. Use `grep -o … | wc -l`, and `grep -F` for any pattern containing `\`, `{` or `}`.

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | No accidental math survives | `grep -o 'class="katex' <page> \| wc -l` on changelog and benchmarks. `[REVISED]` sentinel changed from `application/x-tex`, which only appears when KaTeX **succeeds** and is blind to error spans | 0 on both. Pre-fix baseline to prove the check fires: changelog 3, benchmarks 4 |
| V2 | `katex-error` is gone site-wide | `grep -ro 'katex-error' dist --include='*.html' \| wc -l` | **0**. Pre-fix: 1 (on benchmarks). A standalone check, because a failed parse emits no `<annotation>` |
| V3 | Display math still works | `grep -cF 'begin{aligned}' dist/performance/clumpy/index.html`. `[REVISED]` the original used `\begin{aligned}` as a regex and returns 0 | ≥1, and `clumpy` retains exactly one `application/x-tex` node |
| V4 | The corrupted phrase is restored | `grep -cF '~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour'` on the built changelog. `[REVISED]` the original grepped `1000-sample cohort`, which returns **1 pre-fix** and so could not fail | 0 pre-fix → 1 post-fix. Both verified by review |
| V5 | The benchmarks passages are intact | Compare the **whole** sentence, not a fragment, against the post-fix-with-escapes render. `[REVISED]` the original said "byte-identical apart from the `$` reappearing", which is unachievable — pre-fix contains a KaTeX tree, a red error span, a literal `\~` and mis-scoped `<strong>` | Full sentence present, three correct `<strong>` spans, no `\~` |
| V6 | `$` renders literally where eaten | Inspect `dist/performance/benchmarks/index.html` for `$0.05/vCPU-hour`. `[REVISED]` the original pointed at the changelog page, where that `$` is the unpaired third one and already literal pre-fix — a no-op | Literal, no `<span class="mord">` |
| V7 | The image resolves in the browser | `curl -sI https://www.trimgalore.com/images/pigz_bench.png` after deploy; confirm no `src="docs/` in any built page | 200, `content-type: image/png`. Prove the check can fail: a bogus filename must 404 |
| V8 | The image resolves **on GitHub** | `[RE-SEQUENCED]` View `CHANGELOG.md` on github.com **after the deploy completes**, not before merging. The original wording was unachievable: §3.2 chose a site-absolute URL, which 404s on GitHub until the site serves `/images/…`, so the check would fail on a correct implementation. V7 and V8 are therefore one post-deploy check. Also: GitHub proxies external images through **camo, which caches** — if camo fetches mid-deploy it can cache the 404 past it, so a broken image on the first look is a cache, not a bug | Image renders in the GitHub view |
| V9 | **No strikethrough anywhere** | `grep -ro '<del>' dist --include='*.html' \| wc -l` | **0**. This is the only check that catches §3.3's regression. **Prove it fires**: it must return 1 on a build with the escapes removed, or it is decoration |
| V10 | Only the intended pages changed | `diff -rq "$TMPDIR/before" dist` | Exactly two HTML files differ (`benchmarks`, `changelog`). `[NOTE]` `dist/pagefind/*` also churns under new content hashes and `pagefind-entry.json` differs — expected, since the indexed text changed. `dist/index.html` must be unchanged |
| V11 | Source-side guard | `grep -c 'src="docs/' CHANGELOG.md docs/src/content/docs/reference/changelog.md` | 0. `[ADDED]` catches a future hand-edit reintroducing the path at source, not just in output |
| V12 | CI enforces all four classes, in both workflows | Push, open a PR | `docs-build` green with §4 step 5's step active. Locally verified: extracted from the parsed YAML and run under `bash -e` against four real directories — fixed build **0**, pre-fix build **1**, changelog-corruption-only build **1** (the case the first version passed), missing `dist` **1**; each pattern also proven to fire in isolation |
| V13 | `[ADDED]` The two workflows' step bodies are identical | Extract both `run:` blocks from the parsed YAML and compare | Byte-identical, per `docs.yml:22` |

**Do not add a build-log KaTeX-warning assertion.** Reviewer B measured the real build log at **zero** KaTeX warnings despite KaTeX demonstrably failing (`grep -ic 'LaTeX-incompatible'` → 0), concluding Astro swallows them; Reviewer A reported the warning as observable. The two disagree on the fact and agree on the remedy — grep the built output. Recorded so nobody spends time on the log path. Relatedly, `rehypeKatex` with `strict: 'error'` cannot work: `rehype-katex` catches every error and degrades to a `katex-error` span, so no KaTeX option can turn accidental math into a non-zero exit.

**`docs-build` asserts less than you might assume.** `ci.yml:184-194` *echoes* the HTML page count but only asserts `test "$html" -ge 1`; the OG >10 KB check is genuinely asserted. So V12 means "the build compiles and the three new greps pass", not a page-count guarantee.

V9 and V2 are the two that would quietly not hold — V9 because the regression it catches builds, deploys and renders with no error signal; V2 because failed math is invisible to the obvious sentinel.

---

## 9. Self-Review

**Logic.** Both fixes are small source edits with build-output assertions. The risk was never in the edits but in believing them: an unsupported plugin option and a still-broken image URL both fail *silently* behind a green build. That instinct was right. What the review exposed is that the instinct had not been applied to the plan's own validation table.

**Adjusted after dual agent review.** Reports at `PLAN_review_reviewer-A.md` / `-B.md`. Both returned *do not implement as written*, and both confirmed the diagnosis and the core fix were correct — `singleDollarTextMath: false` verified at four layers, `$$` structurally immune, A2 exhaustively closed across all 30 content entries.

Three critical defects, all found independently by both:

- **§4 step 3 would have introduced a new bug.** The `\~` escapes guard GFM single-tilde strikethrough, not math. Removing them renders the sentence inside `<del>` with a collapsed `**` pair — and every one of the eight original checks passed on it. §3.3 now says the opposite of what it said, and V9 exists solely to catch this.
- **§2.2's headline evidence was false.** `grep -c '1000-sample cohort'` returns 1, not 0. Inherited from an earlier review and stated as verified without checking, which made V4 unable to fail.
- **A second math node on benchmarks renders as a red `katex-error` span** — the more visible defect, absent from the original plan, and structurally invisible to an `application/x-tex` sentinel.

Also corrected: V6/V7 (was V4/V7) were no-ops; V5's criterion was unachievable; V3's pattern failed as a regex; V8's `docs-build` claim over-stated what CI asserts; `grep -c` counts lines; the `dev` branch pin was wrong on both reviewers' evidence; §4 step 4 is one line not two; and Pagefind churn makes "no other file changes" false as written.

**The pattern worth remembering.** Every defect was in a claim the plan presented as settled, and the single worst one came from *trusting another reviewer's measurement without re-deriving it*. That is the same failure the plan itself warns about one level down — and it is now the second time in this feature's history that "verified" has meant "someone checked once". A1's discipline (prove the check can fail) is now attached to V4, V9 and V7 explicitly, not just asserted in prose.

**Not adopted.** Reviewer A's suggestion to drop the `<img>` entirely or degrade it to a Markdown link — defensible, since the image decorates a v0.6.0 entry about Cutadapt/pigz threading that no longer describes the product, but repointing it is cheap and losing the chart loses real information. Recorded rather than dismissed.

**Edge cases.** An unsupported plugin option (V1, A1); other content relying on single-dollar math (A2, V3); the escapes removed by a future cleanup (V9, plus the source comment); `$$` reaching prose (A7); the image 404ing between merge and deploy (§3.2); a stale `docs/dist` poisoning a before/after (§8 method note); a future hand-edit reintroducing `src="docs/` (V11).

**Remaining risks.**

- *Medium:* A5 is narrowed but not closed — emoji shortcodes, footnotes and autolinking were never swept. V10's whole-tree diff catches unexpected page changes, which is the broadest net available without an exhaustive construct audit.
- *Low:* A7's `$$` hazard survives the fix, and nothing watches for it.
- *Low:* the brief post-merge 404 window for the moved asset.
- *Very low:* `smartypants: false` is already set, so no interaction with quote or dash transformation.

---

## 10. Implementation notes

**Implemented 2026-08-01.** Diff: 6 files, +31/−3. Branch not yet created; changes are in the working tree on `dev`.

| Step | File | Done |
|---|---|---|
| 1 | `docs/astro.config.mjs:18` | `remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]]` + one-line comment |
| 2 | `docs/public/images/pigz_bench.png` | `git mv` from `src/assets/screenshots/` (80,712 bytes preserved, recorded as `R` not add+delete) |
| 2 | `CHANGELOG.md:1164`, `changelog.md:974` | `src` → `https://www.trimgalore.com/images/pigz_bench.png`, `alt` added |
| 3 | `benchmarks.md:129` | Escapes **retained**; guard comment added above |
| 5 | `.github/workflows/ci.yml` | New `Assert no silent render corruption` step, three patterns |
| 6 | `CHANGELOG.md` | Two entries under `#### Fixes` in `Unreleased` |

### Iteration log

**#1 — first validation run: two FAILs, both self-inflicted.** `<del>` returned 1 and a leaked `\~` returned 1 post-fix. Cause: the guard comment added in step 3 contained the literal strings `\~` and `<del>`, and HTML comments pass through to built output in `.md` — so the comment tripped the very checks it documented, and the new CI step would have failed the build on this change. Reworded to convey the same fact without the trigger tokens ("Tilde escapes … guard GFM single-tilde strikethrough, not dollar maths"). Both checks clean on rebuild.

**#2 — V9 proven falsifiable before being trusted.** With the config fix applied and the escapes temporarily removed, `<del>` = **1** while `katex-error` = 0 and `class="katex"` on benchmarks = 0. That is a live reproduction of the exact scenario both reviewers described: every math-based check passes while the page ships a struck-through sentence. The escapes were restored and the restoration verified by diff (1 insertion, 0 deletions). The new CI assertion was separately proven to fail on a planted `<del>` and clean after its removal.

### Validation results

Baselines built to a scratch directory, never `docs/dist`, per §8.

| # | pre → post | Result |
|---|---|---|
| V1 | `class="katex"` changelog 3→0, benchmarks 4→0 | PASS |
| V2 | `katex-error` site-wide 1→0 | PASS |
| V3 | clumpy `begin{aligned}` 1→1, `application/x-tex` 1 | PASS |
| V4 | full corrupted phrase 0→1 | PASS |
| V5 | benchmarks sentence intact, 0 leaked `\~` | PASS |
| V6 | benchmarks `class="mord"` 46→0, `$0.05/vCPU-hour` literal | PASS |
| V9 | `<del>` site-wide 0→0, and proven to return 1 when escapes removed | PASS |
| V10 | exactly 2 HTML files differ; `images/` added; 13 pagefind files churn as predicted; `index.html` unchanged | PASS |
| V11 | `src="docs/` in source 2→0 | PASS |
| V7 | image resolves over HTTP | **Deferred** — needs the deploy |
| V8 | image renders on github.com | **Deferred** — needs the deploy, then a camo-cache-aware re-check (V8 as re-sequenced) |
| V12 | CI enforces all four classes in both workflows | **Deferred** — needs a PR. Verified locally against four real directories; see V12 row in §8 |
| V13 | Both workflows' step bodies identical | PASS — extracted from parsed YAML, byte-identical |

Build: 31 pages / 32 HTML files, 6.5–7.8 s, unchanged across the fix.

### Deviations from the plan

1. **`#### Fixes`, not `#### Bug fixes`.** §4 step 6 named a heading that does not exist; `Unreleased` has `#### Changes` / `#### Fixes` / `#### Infrastructure (contributor-facing)`. Used the existing one.
2. **Three `\~` escapes, not two.** Both reviewers cited `\~$0.05` and `\~$41`; there is also `\~$7`. All three retained.
3. **The two new `CHANGELOG.md` entries were *not* mirrored into `changelog.md`.** That page is **197** lines behind (re-derived during coverage: 211 missing lines, of which 14 are the new entries), so syncing only these two would be arbitrary, and the sibling automation closes the whole gap at once. Consistent with the precedent recorded in `plans/phred64-ubam/PLAN.md`. Only the `<img src>` was hand-edited there, exactly as §4 step 4 specifies. Both reviewers independently judged this defensible; the one thing worth tightening is the page's own `:::note` claim that it is "a copy synced with the docs", if the automation is more than a release away.
4. **The `<del>` CI assertion is a blanket site-wide ban.** Strikethrough is a legitimate Markdown feature and nothing in the docs uses it today, so this is a deliberate trade: it makes the §3.3 regression impossible at the cost of blocking any future intentional strikethrough. The failure message now names strikethrough and the `\~` remedy, because the original named only the *built* file — which is not where the fix goes. Revisit if anyone wants intentional strikethrough. Reviewer B's narrower alternative — assert at source that the escapes are still present — was considered and not taken: it guards the one known instance rather than the class.
5. **The CI assertions went into a new step, not the existing `Assert build output`.** `[ADDED AFTER CODE REVIEW]` §4 step 5 said to extend the existing step. A separate step gives clearer failure attribution and keeps the satori-specific comment from acquiring unrelated concerns — both reviewers called it the better choice — but it was a divergence recorded nowhere, which is what this list exists to prevent. All three stage-5 agents flagged the omission independently.

### Iteration log, continued

**#3 — stage-5 review found the guard did not guard.** Coverage returned COMPLETE (36 items, 0 MISSING) because the step conformed exactly to §4 step 5; both code reviewers independently returned *do not ship as written*, because §4 step 5 was itself wrong. The `katex-error` sentinel cannot see a *successful* accidental parse, which is precisely the changelog defect this plan fixes — Reviewer A constructed that tree and the original step exited **0** on it. §8 had already drawn the V1/V2 distinction; only V2 reached CI.

Fixed by adding the `class="katex"`-allowlist check, copying the whole step into `docs.yml` (the publishing path, previously unguarded), adding `test -d dist`, anchoring the image pattern to `<img[^>]*src="[Dd]ocs/`, switching to a `fail` accumulator, and correcting the step comment — `<del>` never shipped to the live site, so "each of these shipped" was false. Re-verified against four real directories: fixed **0**, pre-fix **1**, changelog-corruption-only **1**, missing `dist` **1**.

The §3.3 guard comment was also trimmed to one line; its second sentence narrated CI tooling, against `CLAUDE.md`'s convention, and any future edit naming a trigger token would self-fail the build.

**What this says about the process.** The plan's own §8 was right and its §4 was wrong, and a faithful implementation of the wrong instruction scores as DONE. Coverage cannot catch that by construction — which is the argument for running code review and coverage as independent passes rather than treating either as a substitute for the other.
