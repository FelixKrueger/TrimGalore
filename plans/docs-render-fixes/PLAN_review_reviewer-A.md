# PLAN review — Reviewer A

**Target:** `plans/docs-render-fixes/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Verdict:** the two headline fixes are correct and I verified them end to end. But **step 3 (removing the `\~` escapes) introduces a new silent-corruption bug**, and **three of the eight validations are no-ops that pass before the fix**. Both diagnoses in §2.2 are also factually overstated in one direction and understated in another.

Everything below was re-derived independently. Evidence is cited as file:line or as a command plus its observed output.

---

## Method note (and one environment hazard)

I did not rely on the checked-in `docs/dist`. Timeline observed during this review:

| Time | `docs/dist` state | `astro.config.mjs` state |
|---|---|---|
| 09:42 | pre-fix build (my baseline) | pre-fix |
| 13:05 | rebuilt **post-fix** | fix applied (by another session) |
| 13:10 | rebuilt pre-fix again | reverted to pre-fix |

`docs/dist` is a shared mutable artefact and it was rebuilt three times underneath me, at one point holding a *post-fix* build while the source was *pre-fix*. My "before" numbers were captured from the 09:42 build; my "after" numbers come from a build I ran myself to a separate directory (`npx astro build --outDir <scratchpad>` — this works and leaves `docs/dist` untouched). **Anyone doing a before/after on `docs/dist` in a concurrent session will get a false reading.** See Important #8.

I made no changes to the repo; `git status` is clean for `docs/` and `CHANGELOG.md`.

---

## 1. Logic review

### 1.1 The core fix is right, and I verified it at every layer

A1 is the plan's stated key risk. It holds. Three independent confirmations:

1. **The option is real in the installed tree**, not just in the docs. `docs/node_modules/micromark-extension-math@3.1.0`, `dev/lib/math-text.js:88`:
   ```js
   // Not enough markers in the sequence.
   if (sizeOpen < 2 && !single) {
     return nok(code)
   }
   ```
   with `single` defaulting to `true` at lines 21-25. So `singleDollarTextMath: false` makes a one-dollar opening sequence *reject*, which is exactly the required semantics.
2. **remark-math forwards it.** `node_modules/remark-math@6.0.0/lib/index.js` does `micromarkExtensions.push(math(settings))` — options pass straight through. Only **one** copy of `remark-math` is installed (`find node_modules -type d -name remark-math` → single hit, v6.0.0), so there is no version-skew trap between `package.json`'s `^6.0.0` and what actually runs.
3. **Astro honours the tuple form.** `@astrojs/markdown-remark@7.2.1`, `dist/load-plugins.js:12-14` destructures `[plugin, opts]` for array entries, and `dist/index.js:56-58` applies it as `parser.use(plugin, pluginOpts)`. The `[[remarkMath, {...}]]` syntax in §3.1 is correct.

I then ran a real `astro build` with the fix in place:

| Page | `application/x-tex` before → after | `katex-error` before → after |
|---|---|---|
| `reference/changelog` | 1 → **0** | 0 → 0 |
| `performance/benchmarks` | 1 → **0** | **1 → 0** |
| `performance/clumpy` | 1 → **1** (unchanged) | 0 → 0 |

32 HTML pages both times, build green, 6.67 s. Post-fix rendered text on the changelog page is exactly what §5 promises: `default) — ~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. Display math on `clumpy` survives (`x-tex">\begin{aligned}`, one node). So §3.1, §5.1, §5.2 and §5.4 are all confirmed.

### 1.2 Step 3 is actively harmful — the `\~` escapes are load-bearing

This is my most significant finding. §3.3 asserts:

> `docs/src/content/docs/performance/benchmarks.md:129`'s `\~$0.05` / `\~$41` become unnecessary once single-dollar math is off. `\~` still renders as `~`, so leaving them is harmless — but they are now misleading […] Remove them in the same change.

The premise is wrong. Those escapes are not defending against **math**, they are defending against **GFM single-tilde strikethrough**, which is a completely separate mechanism and is unaffected by the `remarkMath` change:

- `micromark-extension-gfm-strikethrough/readme.md:141` — `singleTilde` (`boolean`, **default: `true`**), confirmed in `dev/lib/syntax.js:23`.
- Astro enables gfm by default (`@astrojs/markdown-remark/dist/index.js:41,48` — `gfm = markdownConfigDefaults.gfm`, then `parser.use(remarkGfm)`), and `docs/astro.config.mjs` does not disable it.

I mirrored Astro's exact plugin order (`remarkParse → remarkGfm → remarkMath{singleDollarTextMath:false} → remarkRehype`) and ran line 129 both ways:

```
=== benchmarks WITH \~ escapes
  <del> count = 0
  <p>On AWS at ~$0.05/vCPU-hour, … roughly <strong>$0.041 with TG</strong> vs
  <strong>$0.007 with the Rust v2 build</strong> — … scales to
  <strong>~$41 with TG vs ~$7 with Rust v2</strong>, with proportional savings …</p>

=== benchmarks WITHOUT escapes (plan step 3)
  <del> count = 1   <<< STRIKETHROUGH INTRODUCED
  <p>On AWS at <del>$0.05/vCPU-hour, … roughly <strong>$0.041 with TG</strong> vs
  <strong>$0.007 with the Rust v2 build</strong> — a 5.9× saving per sample.
  Across a 1000-sample cohort that scales to **</del>$41 with TG vs ~$7 with
  Rust v2**, with proportional savings …</p>
```

So step 3 would trade a KaTeX corruption for a strikethrough corruption of *precisely the same character*: a sentence swallowed into a wrapper element, and a `**` bold pair collapsing to literal asterisks. It would land in the same commit that claims to fix silent prose corruption.

**Worse: none of the plan's checks catch it.** After step 3, `application/x-tex` is 0 (V1, V3 pass) and `katex-error` is 0. V5 is the only candidate, and as written — "grep the built page for the full sentence including `**`-derived `<strong>`" — it is too loose: the first two `<strong>` elements survive intact, so a partial grep passes. Only a grep for the *whole* sentence, or an assertion that `<del>` count is 0, fails.

Note the plan's own framing works against it here: §2.2 reads the escapes as evidence that "someone has hit this before and half-fixed it". `git log -L 129,129` shows they were added in `3af218b` (#282, the `--clumpify` feature commit), in the same edit that renamed "Oxidized" → "the Rust v2 build". Whatever the author's reason, the escapes are doing real work today and removing them regresses the page.

### 1.3 §2.2's headline evidence is wrong, and V4 inherits the error

§2.2 states:

> `grep -c '1000-sample cohort'` on that HTML returns **0** — the phrase is no longer contiguous.

It returns **1**. On the genuine pre-fix build:

```
$ grep -c '1000-sample cohort' dist/reference/changelog/index.html
1
$ grep -o '1000-sample cohort' dist/reference/changelog/index.html | wc -l
1
```

The reason is visible in the annotation: the math node is `<annotation encoding="application/x-tex">7 vs ~</annotation>` — it consumes `$7 vs ~$`, which sits *entirely before* "1000-sample cohort". The pre-fix rendered line (tags stripped, `dist/reference/changelog/index.html:317`) is:

```
default) — ~7vs 7 vs ~7vs 41 per 1000-sample cohort at AWS $0.05/vCPU-hour
```

The prose *is* mangled — two `$` deleted, `7vs` italicised as math variables, content duplicated by KaTeX's MathML+HTML double emission — but "1000-sample cohort" is contiguous throughout. **V4 therefore passes before and after the fix.** Its "Expected: ≥1 (currently **0**)" column is wrong on the parenthetical, which is the half that made it a meaningful check.

This is the exact failure mode §8 demands vigilance about for V1 ("Prove the check can fail first"). The discipline was applied to V1 and not to V4.

### 1.4 V7 is a no-op for the same reason

V7 inspects the built changelog for `$0.05/vCPU-hour` and expects "Literal text, no `<span class="mord">` wrapping". Look at the pre-fix rendered line above: `at AWS $0.05/vCPU-hour` is *already* literal, because it is the unpaired third `$` on that line. V7 passes pre-fix. It would have been a real check if pointed at `performance/benchmarks`, where `$0.05` *is* consumed.

### 1.5 §2.2 and §5 understate the benchmarks damage — there are two nodes, and one renders red

The plan describes one swallowed sentence on `/performance/benchmarks/`. An AST sweep finds **two** `inlineMath` nodes on line 129:

```
docs/src/content/docs/performance/benchmarks.md  dollars=5  math_default=2  math_off=0
   inlineMath @L129: "0.05/vCPU-hour, trimming 84M PE reads at the nf-core default costs roughly **"
   inlineMath @L129: "0.007 with the Rust v2 build** — a 5.9× saving per sample. Across a 1000-sample cohort that scales to **\\~"
```

The second one is not in the plan at all, and it manifests differently: KaTeX cannot parse it, so `rehype-katex` degrades it to an error span. From the pre-fix build:

```html
<span class="katex-error"
      title="ParseError: KaTeX parse error: Expected group as argument to '\~' at end of input: … scales to **\~"
      style="color:#cc0000">0.007 with the Rust v2 build** — a 5.9× saving per sample.
      Across a 1000-sample cohort that scales to **\~</span>
```

So a full sentence is currently rendered in **red error text** on www.trimgalore.com. That is a more serious and more visible defect than the one the plan leads with, and §5's "Visible effects" #2 ("The swallowed sentence … reappears", singular) should be plural. The fix does clear both — post-fix `katex-error` is 0.

### 1.6 §4 step 4 is one line, not two

Step 4 says the tracked page "carries the same two defects at lines 261 and 974" and asks for both to be fixed by hand. But line 261 is the *math* defect, and `docs/src/content/docs/reference/changelog.md` goes through the same Astro pipeline as every other page — so step 1 fixes it with no per-file edit. Confirmed on my post-fix build, where only the image defect survives:

```
$ grep -o 'src="docs/[^"]*"' dist_after/reference/changelog/index.html
src="docs/Images/pigz_bench.png"
$ grep -o 'cohort at AWS \$0.05'  dist_after/reference/changelog/index.html
cohort at AWS $0.05
```

Step 4 reduces to editing line 974's `<img src>`. That marginally strengthens §7 Q1's recommendation (a) — the throwaway edit is one line, not two.

### 1.7 A residual hazard the plan does not mention

`singleDollarTextMath: false` requires **two or more** dollars to open text math (`math-text.js:88`). It does not remove the construct. Confirmed:

```
=== synthetic: prose with two adjacent dollars — "Cost was $$5 and later $$9 per run."
  mathnodes default=1  singleDollarOff=1
```

So `$$` in prose still becomes inline math after the fix. There is no such case in the content today, but `$$` is exactly what a shell snippet uses for a PID, and the sibling plan's whole purpose is to pipe `CHANGELOG.md` in unattended. Worth one sentence in §9's remaining-risks list, since §2.1 claims this fix closes the silent-corruption path before automation.

---

## 2. Assumptions

| ID | Verdict | Evidence |
|---|---|---|
| **A1** | **Confirmed** — stronger than the plan claims | Three layers verified (§1.1). The plan can drop the hedge. |
| **A2** | **Confirmed, and I closed the gap** | Exhaustive sweep, below. |
| **A3** | **Confirmed** | I opened the asset. |
| **A4** | **Confirmed** | No CSP exists anywhere. |
| **A5** | **Mostly closed by me; one new hazard class found** | Below. |

### A2 — exhaustively verified (the plan marks this "partially verified")

I parsed all **30** content entries plus `CHANGELOG.md` with the real `remark-parse` + `remark-math` and enumerated every `math` / `inlineMath` node under both option settings. Complete result:

```
docs/src/content/docs/index.mdx                    dollars=1  math_default=0  math_off=0
docs/src/content/docs/install.md                   dollars=1  math_default=0  math_off=0
docs/src/content/docs/modes/passthrough.md         dollars=5  math_default=0  math_off=0
docs/src/content/docs/performance/benchmarks.md    dollars=5  math_default=2  math_off=0  <<< CHANGES
docs/src/content/docs/performance/clumpy.md        dollars=4  math_default=1  math_off=1
docs/src/content/docs/reference/changelog.md       dollars=6  math_default=1  math_off=0  <<< CHANGES
CHANGELOG.md                                       dollars=6  math_default=1  math_off=0  <<< CHANGES
```

Every other entry has zero `$`. The only node that survives the change is `clumpy.md:89`, and it is a **display** (`math`) node — `$$`-fenced, a different construct, provably unaffected. **Nothing anywhere relies on single-dollar inline math.** A2 can be promoted from "partially verified" to verified.

Two extra points the plan does not make:
- **MDX inherits the config.** `@astrojs/mdx/dist/index.js:106` sets `extendMarkdownConfig: true` by default and line 57 pulls `config.markdown`, so `index.mdx` is covered too. Its lone `$` (line 21, `saves time, $, CO₂`) is currently unpaired and safe — but it is a live tripwire under the current default: any future second `$` in that paragraph would swallow it. A genuine argument *for* the fix that the plan misses.
- The other four `$` in `CHANGELOG.md` (L591 `${basename}`, L931 `$trim_n`, L1132 `-j $cores`) are all inside inline code and immune.

### A3 — confirmed by inspection

I opened `docs/src/assets/screenshots/pigz_bench.png` (PNG, 1810×541, 80712 bytes). It is a two-panel bar chart, "SE trimming with Python3 + pigz (~17.5M reads, 2x75bp)" and "PE trimming with Python3 + pigz", runtime versus `# of cores` (1,2,3,4,6,8). `CHANGELOG.md:1162` is the v0.6.0 entry "Added multi-threading support with the new option `-j/--cores INT` […] if parallel gzip (`pigz`) is installed:", with the `<img>` immediately after. It is unambiguously the intended image. A3 is no longer an inference.

Also worth recording: `git grep pigz_bench` returns only the two broken lines. The asset is referenced nowhere else, so it is currently dead weight in the repo, and it is absent from `dist` (`find dist -iname '*pigz*'` → nothing) because `src/assets` files are only emitted when Astro's image pipeline references them.

### A4 — confirmed, no CSP exists

- No `_headers`, `netlify.toml`, `vercel.json`, or `staticwebapp.config.json`; `docs/public/` contains only `CNAME`, `favicon.svg`, `logos/`.
- No CSP `<meta>`: a grep for `content-security|Content-Security-Policy|csp` across `docs/src`, `docs/astro.config.mjs` and `docs/public` returns nothing, including in the custom `src/components/Head.astro`.
- GitHub Pages adds no CSP.
- The URL is reachable on **both** branches, and I proved the check can fail:
  ```
  dev    → http=200 type=image/png
  master → http=200 type=image/png
  bogus filename → http=404
  ```

So A4 holds, with one unstated consequence worth surfacing: every other image on the site is self-hosted and Astro-optimised. This introduces the site's only cross-origin asset dependency, which is a (minor) privacy and offline-rendering regression. See Optional #11.

### A5 — I closed the math and raw-HTML classes; found a third

**Math class: clean.** Fully enumerated above.

**Raw-HTML class in `CHANGELOG.md`: clean.** Exactly one raw-HTML node — the `<img>` at 1164. I checked the three other angle-bracket candidates and they are false positives, all inside *multi-line* inline code spans (which a per-line grep cannot see):

- `CHANGELOG.md:75` — `provided: <INPUT>...` closes a code span opened on line 74.
- `CHANGELOG.md:978-979` — `` `<git-hash> — <target> — built <ISO-8601 UTC timestamp>` `` spans two lines.

Both render correctly; the built page contains `&#x3C;INPUT` (escaped, 1 occurrence), `git-hash` and `ISO-8601 UTC timestamp`. No content is dropped.

**Broken-asset class: clean apart from the known one.** I enumerated every `src` in the built site. `src="docs/Images/pigz_bench.png"` is the only broken one. The six suspicious relative entries (`src="hero-dark.svg"` etc.) live in `dist/logos/preview.html`, which sits in the same directory as those SVGs, so they resolve — a false lead I chased and cleared.

**New hazard class: GFM single-tilde strikethrough.** Not in the plan and not in the sibling plan's `{` / `<…>` / `$` list. Bare `~` appears on roughly 50 prose lines across the content and `CHANGELOG.md` (`~75% larger`, `~27 threads`, `~552 MiB`, `~−23% wall`, …). Today **zero** produce `<del>` — I verified `<del>` count is 0 across all 32 built pages. But §1.2 shows how narrowly that holds: one line, with its escapes removed, produces a `<del>` immediately. This is a better-evidenced instance of A5's residual risk than the math-only framing, and it is the one most likely to bite the sibling automation.

---

## 3. Efficiency analysis

A non-issue, as expected, but I measured rather than assumed. Full build: **6.67 s**, 31 pages / 32 HTML files, unchanged before and after. The option change makes the tokenizer `nok()` earlier on single dollars and removes two `katex.renderToString` calls, so the fixed build is marginally *cheaper*. `CHANGELOG.md` is ~1200 lines and parses as part of a build that already completes in seconds. No memory or scalability concern at any plausible content size. Nothing to optimise; the plan is right not to discuss it.

---

## 4. Validation sufficiency

The plan's instinct — assert against built output, not source, because both failure modes are silent — is correct. The execution has four holes, three of which mean a check cannot fail.

| # | Status | Problem |
|---|---|---|
| V1 | **Weak sentinel** | `application/x-tex` cannot see failed math (below). Otherwise sound, and the "prove it fires pre-fix" instruction is right. |
| V2 | OK | Verified: `clumpy` keeps 1 node, `x-tex">\begin{aligned}` present on the same line, so a single-line grep works. |
| V3 | **Weak sentinel + method** | Same blind spot as V1; also see Important #8 on the baseline. |
| V4 | **No-op** | Passes pre-fix (returns 1, not 0). §1.3. |
| V5 | **Too loose to rely on** | It is the only check that could catch the step-3 regression, and only if the grep covers the *entire* sentence. "including `**`-derived `<strong>`" is satisfied by the two surviving `<strong>`s. §1.2. |
| V6 | OK | I confirmed 200 + `image/png`, and confirmed a bogus path 404s. Add the GitHub-side check (Important #9). |
| V7 | **No-op** | `$0.05/vCPU-hour` already renders literally on that page pre-fix. §1.4. |
| V8 | **Asserts less than claimed** | Below. |

### The sentinel is wrong

`application/x-tex` only appears when KaTeX *succeeds*. `rehype-katex@7.0.1/lib/index.js` never throws — it calls `katex.renderToString(…, {throwOnError: true})`, catches, retries with `{strict:'ignore', throwOnError:false}`, and if that also fails emits `<span class="katex-error" style="color:#cc0000">`. That span carries **no** `<annotation>`. Pre-fix ground truth makes the gap concrete:

| Page | `application/x-tex` | `katex-error` |
|---|---|---|
| `reference/changelog` | 1 | 0 |
| `performance/benchmarks` | 1 | **1** |
| `performance/clumpy` | 1 | 0 |

Three math nodes are visible to the plan's grep; the fourth — the red sentence on the benchmarks page — is invisible to it. So V1 could read 0 while prose is still being eaten. Use `class="katex` (covers `katex`, `katex-mathml`, `katex-display`, `katex-error`) for the "did any math happen" question, and assert `katex-error` is **0 site-wide** as a check in its own right. For reference, `class="katex` counts pre-fix are clumpy 4, benchmarks 4, changelog 3.

Also, a precision point: the plan uses `grep -c`, which counts *matching lines*, not occurrences. These pages put a whole document on very few lines — `dist/reference/changelog/index.html` is 1224 lines for a 1200-line changelog, and both math nodes on the benchmarks page are on line 439. It happens to agree today (one node per page), but V3's "count per file before and after" will silently under-report the "before" side the moment a page has two nodes on one line. Use `grep -o … | wc -l`.

### V8 asserts much less than the plan believes

V8 expects "Green; 32 HTML pages; OG images still >10 KB". The actual job, `.github/workflows/ci.yml:184-194`:

```yaml
test -d dist/og
html=$(find dist -name '*.html' | wc -l)
png=$(find dist/og -name '*.png' | wc -l)
echo "$html HTML pages, $png OG images"
test "$html" -ge 1
test "$png"  -ge 1
runt=$(find dist/og -name '*.png' -size -10240c)
if [ -n "$runt" ]; then …; exit 1; fi
```

The page count is **echoed, not asserted** — `test "$html" -ge 1` stays green if the site collapses to a single page. The OG >10 KB half is genuinely asserted. `.github/workflows/docs.yml` has no assertions at all beyond a successful build. So V8 is "the build compiles", not the content guarantee the plan implies.

### Where this could still silently pass every check

1. **Step 3 lands and regresses the benchmarks page.** `x-tex` 0, `katex-error` 0, V4/V7 no-ops, V5 loose. Nothing fails. This is the concrete realisation of the risk §9 is worried about, and it comes from the plan's own step 3.
2. **The image URL is fixed for Astro but not checked on GitHub.** V6 curls the URL and greps `dist`. The entire justification for choosing an absolute URL (§3.2, §7 Q3) is that one line must serve two renderers — but only one renderer is validated.
3. **A `$$` pair reaches prose later.** Still becomes math (§1.7); no check watches for it.

---

## 5. Alternatives

**§2.4 — fix at source, not in the generator.** Adopting Reviewer B's position is right, and for the reason given: keeping the generator's "one transformation, cannot reach content" invariant is worth more than saving one line in `CHANGELOG.md`. No disagreement.

**§7 Q2 — the branch in the image URL. I would choose `master`, not `dev`.** The plan justifies `dev` as "consistent with the sibling plan's branch decision". But the repo's two existing precedents for exactly this — an absolute GitHub URL embedded in a file that renders on github.com — both use `master`:

- `README.md:3-4` — `raw.githubusercontent.com/FelixKrueger/TrimGalore/master/docs/public/logos/hero-*.svg`. A4 cites the README as precedent; the README uses `master`.
- `docs/astro.config.mjs:38` — `editLink.baseUrl: 'https://github.com/FelixKrueger/TrimGalore/edit/master/docs/'`.

I verified the asset exists and serves 200 on **both** branches, so there is no availability argument for `dev`. `dev` is by construction the less stable ref, and the line being edited is in `CHANGELOG.md`, which people read on `master`. This is a judgement call rather than a correctness defect, but the plan's stated rationale ("consistent with…") is the weaker of the two consistency arguments available.

**§7 Q3 — Astro-optimised image.** Correctly rejected, and the mechanism check supports it: `remarkCollectImages` runs *after* user remark plugins (`@astrojs/markdown-remark/dist/index.js:59`), so `![](../../assets/screenshots/pigz_bench.png)` really would be optimised — but GitHub cannot resolve it, and the line has to serve both. The reasoning holds.

Two options the plan did not consider, both cheap:

- **Drop the `<img>`, or degrade it to a Markdown link.** The image has been broken in both renderers "for some time" (§2.3 — confirmed; it is referenced nowhere else in the repo), it decorates a v0.6.0 entry about Cutadapt/pigz threading that no longer describes the product, and a link renders correctly in both renderers with no external dependency and no branch pin. Worth one line in §7 Q3 as the rejected-but-considered baseline.
- **Fix the last-mile invariant in CI instead of trusting a one-off grep.** §9 defers "whether the docs pipeline should validate rendered output against known-hazardous constructs" as "a bigger idea than either bug". The narrow version is ~4 lines appended to the existing "Assert build output" step, and it makes both of today's bug classes permanently non-recurring:
  ```bash
  ! grep -rq 'katex-error' dist --include='*.html'
  ! grep -rq 'src="docs/'  dist --include='*.html'
  ```
  Optionally a `<del>`-count guard, given §1.2. I would take this over V1/V3/V4/V7 combined: it is the only proposal here that catches the *next* instance rather than this one.

One alternative I checked and would **not** recommend: making KaTeX itself fail the build (`rehypeKatex` with `strict: 'error'`). It cannot work — `rehype-katex` catches every error and degrades to `katex-error` (source quoted in §4), so no katex option can turn accidental math into a non-zero exit. It does call `file.message(…)`, so Astro emits a build warning; on the pre-fix content that warning is real and observable (`LaTeX-incompatible input and strict mode is set to 'warn': Unrecognized Unicode character "—"`), and it disappears after the fix. But a warning nobody fails on is what let this ship in the first place. The post-build grep is the robust form.

---

## 6. Action items

### Critical

1. **Drop step 3. Do not remove the `\~` escapes in `benchmarks.md:129`.** They guard against GFM single-tilde strikethrough (`singleTilde` default `true`), not against math, so the change is not a no-op cleanup — it renders the sentence as `<del>…</del>` and collapses a `**` bold pair to literal asterisks (§1.2). If the "misleading escape" concern in §3.3 is worth addressing, address it with a comment or an HTML entity, not by deleting the escape. Whatever is decided, §3.3's premise ("now unnecessary… verify the rendered text is unchanged") must be corrected — the rendered text is *not* unchanged.
2. **Replace V4 and V7; they pass before the fix.** `grep -c '1000-sample cohort'` returns 1 pre-fix, and `$0.05/vCPU-hour` is already literal on the changelog page pre-fix (§1.3, §1.4). Assertions that do flip: on `reference/changelog`, that the built page contains the literal `~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour` (verified present post-fix, absent pre-fix); for the `$0.05`-eaten case, point the check at `performance/benchmarks`. Apply §8's own "prove the check can fail" rule to every row, not just V1.
3. **Change the sentinel in V1 and V3** from `application/x-tex` to `class="katex`, and add a standalone assertion that `katex-error` is 0 across `dist`. As written, V1 cannot see a math node that KaTeX failed to parse — which is exactly what one of the two live defects is (§1.5, §4).

### Important

4. **Correct §2.2 and §5.** The benchmarks page has **two** `inlineMath` nodes, and the second renders as a red `katex-error` span containing a full sentence — arguably the worse of the two live defects and currently unmentioned. §5's "Visible effects" #2 should be plural.
5. **§4 step 4 is one line, not two.** The math defect on `changelog.md:261` is fixed by the config change alone; only the `<img src>` at line 974 needs a hand edit (§1.6). Worth updating §7 Q1's cost argument accordingly — it strengthens recommendation (a).
6. **Fix V8's expectation.** `ci.yml:184-194` asserts `html -ge 1`, not 32 pages. Either assert the page count there or stop citing it as one (§4).
7. **Reconsider the branch pin (§7 Q2) in favour of `master`.** Both in-repo precedents (`README.md:3-4`, `astro.config.mjs:38`) use `master`; the asset serves 200 on both branches, so nothing is gained by the less stable ref (§5).
8. **Specify how the V3 baseline is captured.** `docs/dist` is shared and mutable; it was rebuilt three times during this review, once holding a post-fix build against pre-fix source. Build the baseline to a scratch directory (`npx astro build --outDir <tmp>` — verified working, leaves `docs/dist` alone) and diff two immutable trees using real temp files. Also switch `grep -c` to `grep -o … | wc -l`, since these pages put whole documents on very few lines.
9. **Add a GitHub-side check for the image.** V6 validates only the Astro side. Since the whole point of the absolute URL is dual-renderer support, view `CHANGELOG.md` on github.com on the branch before merging.

### Optional

10. **Add the two-line CI guard** (`! grep -rq 'katex-error' dist`, `! grep -rq 'src="docs/' dist`) to the existing "Assert build output" step. This is the only item here that prevents recurrence rather than fixing one instance, and it is a much smaller change than §9's framing suggests (§5).
11. **Add `alt` to the `<img>`.** §3.2 deliberately leaves the malformed `style` and space-bearing `id` alone, which is a reasonable scope call, but it does not mention that the tag has no `alt` — only `title`. One word, in a tag already being edited. Related: note in A4 that this becomes the site's only cross-origin asset.
12. **Promote A2 to verified and record the new hazard class in A5.** A2 is exhaustively confirmed across all 30 entries plus `CHANGELOG.md` (§2). For A5, add GFM single-tilde strikethrough to the known-divergence list alongside `{` / `<…>` / `$` — bare `~` appears on ~50 prose lines, currently producing zero `<del>` by luck rather than design, and it is the class most likely to bite the sibling automation.
13. **Note the residual `$$` hazard in §9.** `singleDollarTextMath: false` requires two-or-more dollars, so `$$` in prose still becomes math (§1.7). Relevant because §2.1 sells this fix as closing the silent-corruption path before automation.
