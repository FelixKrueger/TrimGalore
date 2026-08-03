# PLAN review — Reviewer B

**Target:** `/Users/fkrueger/Github/TrimGalore/plans/docs-render-fixes/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Date:** 2026-08-01

## Method / what I actually ran

Everything below is measured, not read off documentation.

- Read the installed packages in `docs/node_modules`, not just `docs/package.json`.
- Reproduced both bugs in the existing `docs/dist/` build.
- Ran a faithful before/after render probe through Astro's own `createMarkdownProcessor` (`@astrojs/markdown-remark` from `docs/node_modules`), with the same plugin set as `astro.config.mjs`.
- Applied the §3.1 config change for real, ran `npm run build`, diffed the whole `dist/` tree against a `tar` snapshot of the pre-fix build, then reverted `astro.config.mjs` and rebuilt. `git status` on `docs/` and `CHANGELOG.md` is clean; `dist/` is back in its pre-fix state.
- Proved every comparison could fail before believing a pass: appended a `<!--SENTINEL-->` to `dist/index.html` in the snapshot and confirmed `diff -q` reported it (`rc=1`). No process substitution used; real temp files under `$TMPDIR`. `command grep` throughout (`grep` here is a shell function, and `dist/` is gitignored via `docs/.gitignore:2`).

**Verdict: the two fixes are correct and the diagnosis of Bug 1 and Bug 2 is real — but one instructed step (§4 step 3) introduces a new rendering bug, and the validation table would pass green on that bug. Do not implement as written.**

---

## 1. Logic review

### L1 — CRITICAL: §3.3 / §4 step 3 is actively harmful. The `\~` escapes guard GFM strikethrough, not math.

The plan says the `\~$0.05` / `\~$41` escapes in `docs/src/content/docs/performance/benchmarks.md:129` "become unnecessary once single-dollar math is off", that "`\~` still renders as `~`, so leaving them is harmless — but they are now misleading", and instructs (§4 step 3) to "drop the now-redundant `\~` escapes".

They are not redundant. They defend against `micromark-extension-gfm-strikethrough`'s single-tilde strikethrough. Removing them with the math fix in place turns the sentence into a struck-through one.

Measured, `singleDollarTextMath: false` applied, same processor Astro uses:

```
--- benchmarks line WITH \~ escapes (current source)
    del=false  katexError=false  texNodes=0
    <p>On AWS at ~$0.05/vCPU-hour, … costs roughly <strong>$0.041 with TG</strong> vs
       <strong>$0.007 with the Rust v2 build</strong> — … scales to
       <strong>~$41 with TG vs ~$7 with Rust v2</strong>, with proportional savings.</p>

--- benchmarks line WITHOUT \~ escapes (plan §3.3 proposes this)
    del=true   katexError=false  texNodes=0
    <p>On AWS at <del>$0.05/vCPU-hour, … costs roughly <strong>$0.041 with TG</strong> vs
       <strong>$0.007 with the Rust v2 build</strong> — … scales to **</del>$41 with TG
       vs ~$7 with Rust v2**, with proportional savings.</p>
```

Why it bites this line but not the changelog line, precisely: strikethrough uses attention-style flanking. In `**~$41`, the `~` is preceded by `*` — non-whitespace, so it is a valid right-flanking closer, and it pairs with the opener in `at ~$0.05`. In `— ~$7 vs ~$41` (`CHANGELOG.md:451`, `changelog.md:261`) the second `~` is preceded by a space, so it is not right-flanking and cannot close. Confirmed empirically: that line renders `del=false` after the fix, and the full string `~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour` appears verbatim in the rebuilt `dist/reference/changelog/index.html`.

**No check in V1–V8 catches this.** After removing the escapes: V1 → 0 `application/x-tex` ✓; V3 → benchmarks drops to 0 ✓; V5 ("grep the built page for the full sentence including `**`-derived `<strong>`") → `<strong>$0.041 with TG</strong>` is present ✓; V7 → `$0.05/vCPU-hour` literal with no `<span class="mord">` ✓. All eight go green while the published page ships a struck-through sentence with two orphaned `**`.

**Action:** delete §4 step 3 and rewrite §3.3 to say the opposite — the escapes are load-bearing and must stay. If a comment is wanted next to them, it should say "guards GFM single-tilde strikethrough", not "escaping for `$`". Add a `<del>` count to validation (see §4 below).

Note for the record: this is a *third* GitHub-vs-Astro-relevant construct class, which is exactly what A5 warned might exist. A5's caution was justified — for a construct the audit never considered.

### L2 — CRITICAL: §2.2's central evidence claim is false, and V4 is therefore a check that cannot fail.

§2.2 states: "`grep -c '1000-sample cohort'` on that HTML returns **0** — the phrase is no longer contiguous." V4 restates it as "≥1 (currently **0**)".

Measured on `docs/dist/reference/changelog/index.html`:

| build | `command grep -c '1000-sample cohort'` |
|---|---|
| pre-fix (as shipped) | **1** |
| post-fix (my real build) | **1** |

The reason is visible in the emitted HTML: the math node is `<annotation encoding="application/x-tex">7 vs ~</annotation>`, which ends *before* `41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. That text is a single contiguous run in the pre-fix output. The phrase was never split. (The plain `grep` shell function and `command grep` both return 1, so this is not a gitignore-honouring-wrapper artefact — I checked that specifically.)

What is actually lost pre-fix is two `$` characters and one `~`, plus italicised "7vs". The assertion that does discriminate, verified on both builds:

| pattern | pre-fix | post-fix |
|---|---|---|
| `~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour` | 0 | 1 |
| `$7 vs` | 0 | 1 |
| `$41 per` | 0 | 1 |

**Action:** correct §2.2's claim, and replace V4's pattern with the full dollar-bearing phrase (use `grep -F`).

### L3 — CRITICAL: Bug 1 is worse than §2.2 says, and the whole class is invisible to V1/V3.

§2.2 describes the benchmarks page as "an entire sentence … swallowed into a math node". There is a **second** math node on that page, and it is a KaTeX parse failure rendered in red on the live site. From `docs/dist/performance/benchmarks/index.html`:

```html
<strong><span class="katex-error"
  title="ParseError: KaTeX parse error: Expected group as argument to '\~' at end of input: … scales to **\~"
  style="color:#cc0000">0.007 with the Rust v2 build** — a 5.9× saving per sample. Across a
  1000-sample cohort that scales to **\~</span>41 with TG vs ~$7 with Rust v2</strong>
```

So the published page currently shows a sentence in `#cc0000`, a literal `\~` leaking through, orphaned `**`, and `<strong>` spanning the wrong range.

The load-bearing point: **`katex-error` spans emit no `application/x-tex` annotation.** V1, V3 and V7 all key on `application/x-tex`, so they are structurally blind to error nodes. Measured across all 32 pages:

| page | `application/x-tex` | `katex-error` |
|---|---|---|
| `dist/performance/benchmarks/index.html` | 1 | **1** |
| `dist/performance/clumpy/index.html` | 1 | 0 |
| `dist/reference/changelog/index.html` | 1 | 0 |

V3's "count `application/x-tex` per file before and after" would have reported benchmarks as having exactly one math node before and zero after — never noticing the error node existed. If a future `$` produced only an error node, V1 would report 0 and pass.

I also tested the obvious alternative signal and it does not work: the real `astro build` log contains **zero** KaTeX warnings (`grep -ic 'LaTeX-incompatible'` on the full pre-fix build log → 0), even though KaTeX demonstrably failed. My standalone probe *did* print `LaTeX-incompatible input … Unrecognized Unicode character "—"`, so Astro swallows them. Do not add a build-log assertion; grep the output for `katex-error`.

**Action:** V1 and V3 must count `katex-error` and `class="katex` alongside `application/x-tex`.

### L4 — IMPORTANT: §4 step 3 / V5's success criterion is unachievable as written.

§3.3 says "verify the rendered text is unchanged"; §4 step 3 says "confirm the rendered sentence is byte-identical apart from the `$` characters reappearing". Nothing about that sentence is byte-identical across the fix. Pre-fix it contains a KaTeX MathML/HTML tree, a red error span, a literal `\~`, and mis-scoped `<strong>`; post-fix it is clean prose with three correct `<strong>` spans. The comparison as specified will "fail" on a correct fix, which invites an implementer to conclude the fix is wrong.

The comparison the plan presumably means — and the only one that is meaningful — is **post-fix-with-escapes vs post-fix-without-escapes**. That is precisely the comparison that exposes L1. State the baseline explicitly.

### L5 — IMPORTANT: V3's "no other file changes" is false, benignly.

`diff -rq` of the pre-fix snapshot against the post-fix build. Exactly two HTML files differ:

```
Files …/dist/performance/benchmarks/index.html and dist/performance/benchmarks/index.html differ
Files …/dist/reference/changelog/index.html  and dist/reference/changelog/index.html  differ
```

but the Pagefind search index also churns — `dist/pagefind/fragment/*.pf_fragment` and `dist/pagefind/index/*.pf_index` are added/removed under new content hashes, `pagefind-entry.json` differs, and `pagefind.en_*.pf_meta` is renamed. That is expected (the indexed text changed) and not a regression. Say so in V3, or the implementer will chase it. Notably `dist/index.html` (the MDX homepage, which has a lone `$` — see A2) is **unchanged**, which is the right result.

### L6 — OPTIONAL: V2's expected string is regex-fragile.

V2 expects the clumpy annotation "containing `\begin{aligned}`". `{`, `}` and `\` are regex-significant. My literal attempt `command grep -c 'application/x-tex">\\begin{aligned}' dist/performance/clumpy/index.html` returned **0**, while `command grep -c 'begin{aligned}'` returned **1**. The annotation really is `\begin{aligned}\nn_{\text{bins}} …`. Specify `grep -F`.

### L7 — OPTIONAL: V1/V3 use `grep -c`, which counts matching *lines*, not occurrences.

Starlight emits a whole paragraph on one line, so two math nodes in one paragraph count as 1. Today it happens not to mislead — I verified `occurrences == lines == 1` for all three affected files — but the check is one added `$` away from under-counting. Use `command grep -o … | wc -l`.

---

## 2. Assumptions

### A1 — CONFIRMED, and stronger than the plan claims

The plan flags this as the key risk. It holds.

- `docs/node_modules/remark-math/package.json` → `"version": "6.0.0"`, consistent with the `^6.0.0` declaration in `docs/package.json:26`. I read the installed tree, not the manifest alone.
- `docs/node_modules/remark-math/lib/index.js:39` passes options straight through: `micromarkExtensions.push(math(settings))`.
- `docs/node_modules/micromark-extension-math/lib/math-text.js:18-21` reads `options_.singleDollarTextMath` (defaulting to `true`), and line 79 enforces it: `if (sizeOpen < 2 && !single) { return nok(code) }`.
- Declared in the type surface: `micromark-extension-math/index.d.ts:33` and `mdast-util-math/lib/index.d.ts:29-38` (`ToOptions` is exactly `{ singleDollarTextMath?: boolean | null | undefined }` — which is the type `remark-math` uses for its own options parameter).

Stronger than §3.1's empirical argument: `micromark-extension-math/index.js` hands options **only** to `mathText`, never to `mathFlow` (`{ flow: { 36: mathFlow }, text: { 36: mathText(options) } }`). `$$` display math is therefore *structurally* immune to this option, not merely observed to survive. Confirmed by real build: `dist/performance/clumpy/index.html` keeps its one math node containing `\begin{aligned}`; `$$a=b$$` inline also still renders.

End-to-end verification with the change applied: build green, `31 page(s) built`, `Found 32 HTML files`, changelog and benchmarks both drop to zero math and zero error nodes, clumpy retains its node.

**One thing to add to A1.** The `// @ts-check` pragma at `docs/astro.config.mjs:1` gives **no** protection against a typo'd option name: Astro types plugin options as `any` (`docs/node_modules/@astrojs/internal-helpers/dist/markdown.d.ts:24` — `RemarkPlugins = (string | [string, any] | RemarkPlugin | [RemarkPlugin, any])[]`), and there is no `astro check` or `tsc` anywhere — `docs/package.json` scripts are `dev/start/build/preview/astro/logos`, and both workflows run only `npm run build` (`.github/workflows/ci.yml:182`, `.github/workflows/docs.yml:44`). The plan's insistence that V1 assert against built output is correct, and this is the concrete reason.

### A2 — CONFIRMED exhaustively (this closes V2's gap at plan time)

I swept all 30 content entries plus `CHANGELOG.md`, blanking fenced code blocks, indented code blocks and inline code spans (line-wise, so it over-reports rather than hides). Complete inventory of prose-level `$`:

| file | line | count | nature |
|---|---|---|---|
| `docs/src/content/docs/performance/benchmarks.md` | 129 | 5 | hazard — currently broken |
| `docs/src/content/docs/performance/clumpy.md` | 93, 100 | 4 | `$$` fences — **the only legitimate math in the site** |
| `docs/src/content/docs/reference/changelog.md` | 261 | 3 | hazard — currently broken |
| `CHANGELOG.md` | 451 | 3 | hazard — currently broken |
| `docs/src/content/docs/index.mdx` | 21 | 1 | lone `$` in `saves time, $, CO₂` — inert, latent |

16 prose-level `$` total. **Nothing anywhere relies on single-dollar inline math.** A2 can be promoted from "partially verified" to verified, and V2 reduced to a regression guard.

The plan does not mention `index.mdx:21`. It is harmless today (a single `$` cannot pair; verified `dist/index.html` is byte-identical across the fix) but it belongs in the inventory: it is the one place where a *future* second `$` in the same paragraph would silently pair, and it sits in `.mdx`, which the plan never establishes is covered. It is — Astro's `markdown.remarkPlugins` extends to MDX by default and there is no `extendMarkdownConfig: false` in `astro.config.mjs`.

### A3 — CONFIRMED, not merely inferred

A3 says the asset is "inferred from the filename … not confirmed against the original intent". I opened `docs/src/assets/screenshots/pigz_bench.png`: two bar charts, *"SE trimming with Python3 + pigz (~17.5M reads, 2x75bp)"* and *"PE trimming with Python3 + pigz"*, y-axis `time [mm:ss]`, x-axis `# of cores` at 1/2/3/4/6/8. That is exactly `title="Multi-threading benchmark"` and exactly the subject of the `CHANGELOG.md:1164` entry (`-j/--cores`, Python 3, pigz). A3 is verified; drop the hedge.

Extra fact worth having: the asset is an **orphan**. Its only references in the repo are the two broken `<img src="docs/Images/…">` lines (`CHANGELOG.md:1164`, `docs/src/content/docs/reference/changelog.md:974`), and it is absent from the build's 22 optimised images (`ls dist/_astro | grep -ci pigz` → 0). That is *why* `find docs/dist -iname '*pigz*'` is empty — Astro only processes `src/assets` files that are imported. It also means the asset can be moved into `public/` without breaking anything (relevant to the alternative in §5).

### A4 — CONFIRMED

```
$ curl -sSI https://raw.githubusercontent.com/FelixKrueger/TrimGalore/dev/docs/src/assets/screenshots/pigz_bench.png
HTTP/2 200
content-type: image/png
cache-control: max-age=300
content-security-policy: default-src 'none'; style-src 'unsafe-inline'; sandbox
x-frame-options: deny
```

Reachable and correctly typed. The CSP in that response governs the fetched resource, not third-party `<img>` embedding; `x-frame-options: deny` blocks framing, not images. No site-side restriction exists: no `Content-Security-Policy` or `img-src` anywhere in `docs/src`, `docs/public`, `docs/astro.config.mjs` or `docs/scripts`; no `_headers`, `netlify.toml` or `vercel.json`; and the site deploys via GitHub Pages (`.github/workflows/docs.yml`), which injects none.

One correction to A4's supporting claim. It says raw.githubusercontent "is already how the repo's README serves images" — true, but the README pins **`master`**, not `dev` (`README.md:3-4`, the `hero-dark.svg`/`hero-light.svg` `<source>`/`<img>` pair). That matters for §7 Q2; see §5.

### A5 — substantially narrowed, but it missed the construct that actually bites

I swept `CHANGELOG.md` and `docs/src/content/docs/reference/changelog.md` for the other divergence classes, again with code blocks and inline code spans blanked:

- **Raw HTML outside code:** exactly one hit in each file — the `<img>` at `CHANGELOG.md:1164` / `changelog.md:974`. Nothing else. Every `<input.bam>`, `<INPUT>`, `<stem>`, `<git-hash>`, `<FILE>`-style placeholder is inside backticks or a fence and therefore immune.
- **Relative (non-`http`, non-anchor) link and image targets outside code:** **zero** in both files.

Combined with the A2 table, A5 holds for the `$`, raw-HTML and relative-link classes, and can be stated that way with evidence instead of "not verified". Classes I did **not** sweep: emoji shortcodes (`:warning:`), footnotes, bare-URL autolinking. `{`-brace/JSX hazards cannot bite the changelog page — it is `.md`, not `.mdx`.

But A5's own hedge is vindicated by L1: **GFM single-tilde strikethrough is a fourth divergence class**, it is live in the source today, and the plan's §4 step 3 would activate it. A5 should name it.

---

## 3. Efficiency analysis

Non-issue, and I have numbers rather than a shrug. Full `npm run build`: **7.3–8.7 s**, 31 pages / 32 HTML files, all 22 optimised images served from cache, Pagefind index in 279 ms. The config change *removes* work — three KaTeX renders and one KaTeX parse failure become one. Nothing here scales with content in a way that matters.

The only build-time consideration worth a line in the plan is L5's Pagefind churn, which is a correctness-of-expectation point, not a cost one.

---

## 4. Validation sufficiency

V1–V8 are the right *shape* — assert against built output, prove the check can fail — but as written the table has one vacuous check, one blind spot that hides a live defect, and one impossible criterion, and it passes green on the bug §4 step 3 introduces.

Summary of where it silently fails:

| # | Problem | Finding |
|---|---|---|
| V1 | Blind to `katex-error` nodes (they carry no `application/x-tex`) | L3 |
| V1/V3 | `grep -c` counts lines, not occurrences | L7 |
| V2 | `\begin{aligned}` used as a regex; returns 0 as written | L6 |
| V3 | Blind to `katex-error`; "no other file changes" contradicted by Pagefind churn | L3, L5 |
| V4 | Pattern returns 1 pre-fix and post-fix — cannot fail | L2 |
| V5 | Passes on the struck-through output; baseline unachievable | L1, L4 |
| — | **No check anywhere for GFM strikethrough** | L1 |

Checks to add:

1. **`<del>` guard (this is the one that matters).** `command grep -c '<del>' ` across all of `dist/**/*.html` → expect **0**. Measured: 0 in both my pre-fix and post-fix builds with the escapes retained; **1** if the escapes are removed. This is the only assertion that catches L1.
2. **`katex-error` in V1/V3.** Count `application/x-tex`, `katex-error` and `class="katex` per file, before and after. Pre-fix baseline to assert against: benchmarks 1/1, clumpy 1/0, changelog 1/0. Post-fix: clumpy 1/0 only, everything else zero.
3. **Occurrence counting, and `-F`.** `command grep -o -F … | wc -l` for V1/V2/V3/V4.
4. **V4's replacement pattern.** `grep -cF '~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour'` → 0 pre-fix, 1 post-fix. Both verified.
5. **Source-side guard for Bug 2.** V6 already greps built pages for `src="docs/`; add the same grep on `CHANGELOG.md` and `docs/src/content/docs/reference/changelog.md` so a future hand-edit reintroducing the path is caught at source. Pre-fix baseline: `command grep -c 'src="docs/' docs/dist/reference/changelog/index.html` → 1.
6. **Do NOT add a build-log KaTeX-warning assertion.** I tried it: the real build log has zero KaTeX warnings even though KaTeX demonstrably failed (`grep -ic 'LaTeX-incompatible'` → 0 on the full pre-fix log). Astro swallows them. Recording this so nobody spends time on it.

One process note the plan gets right and should keep: it demands the V1 grep be shown to fire on the pre-fix build. Extend that discipline to the new `<del>` check — it must be shown to return 1 on an escapes-removed build, or it is decoration.

---

## 5. Alternatives

**§2.4 — fix at source, not in the generator. Agree with Reviewer B's position, and I can strengthen it.** The plan argues from the generator's auditability invariant. Add a measurement: the `<img>` at `CHANGELOG.md:1164` is the *only* raw-HTML construct outside code in the entire file, and there are zero relative link/image targets. So a source fix is a single-site repair with nothing to generalise, whereas a generator rewrite rule would exist to serve exactly one line forever. The rejection is correct.

**§7 Q2 — branch pin. I would take `master`, not `dev`.** This is a judgement difference, but the plan's stated rationale is contradicted by the repo. Facts: both URLs return 200 and `git rev-parse dev:… master:…` gives the same blob `06baf8c2635e887e295ebb2350ba543a90b3b0ca`. Two existing precedents both pin `master` — `README.md:3-4` (the very precedent A4 cites, without noting its branch) and `astro.config.mjs:38` `editLink.baseUrl` → `/edit/master/docs/`. And this URL is embedded in a *historical* changelog entry that is also read from release tarballs and at tags, where `master` (release pin) is the more stable target than `dev` (work trunk). The plan says `dev` is "consistent with the sibling plan's branch decision"; that consistency argument applies to `editUrl`/banner behaviour, not to a static asset in a frozen historical entry.

**§7 Q3 — the rejection is under-argued: it refutes one option and dismisses a family.** The plan rejects "move the image into the docs site properly" on the grounds that `![](../../assets/screenshots/pigz_bench.png)` needs Astro's asset pipeline, which does not exist on GitHub. Correct — but that is only one member. It does not consider **`docs/public/` plus a site-absolute URL**:

```
https://www.trimgalore.com/images/pigz_bench.png
```

That needs no Astro pipeline, is one line serving both renderers, is self-hosted (no third-party host, no `cache-control: max-age=300`, no raw.githubusercontent rate limiting), and **eliminates §7 Q2's branch-pinning risk entirely** — the failure mode "the asset moves on `dev` and the image silently breaks" disappears. The pattern is already established: `docs/public/logos/` exists and is what the README's raw URLs point at. Cost: the PNG must be copied or moved into `docs/public/`, which breaks nothing because the asset is an orphan (A3), plus a brief post-merge/pre-deploy window where the URL 404s. I would take this; either is defensible. The actionable point is that Q3 should say why the whole family loses, not just the import-based member.

**Not-considered alternative: disable GFM single-tilde at the config layer.** §2.2 argues, rightly, that "hand-escaping is the wrong layer: it requires every future author to know". That argument applies verbatim to `\~`, so the symmetric fix would be `remark-gfm`'s `singleTilde: false`. **I checked, and it is not worth it here:** Astro exposes `markdown.gfm` as a boolean only, and it is deprecated in this version (`docs/node_modules/astro/dist/types/public/config.d.ts:2267-2283`, `@astrojs/internal-helpers/dist/markdown.d.ts:41-42`). Getting sub-options would mean `gfm: false` plus re-adding `remark-gfm` manually — which risks tables (used in `benchmarks.md`), autolinks and Starlight's expectations, to defend one line. Recommendation: keep the escapes, and rely on the `<del>` check from §4 to make the hazard visible instead of invisible. Worth one sentence in the plan so the asymmetry with the `$` fix is a recorded decision rather than an oversight.

**Minor, since the line is being edited anyway:** the `<img>` has no `alt`. §3.2 deliberately leaves the malformed `style` and the space-bearing `id` alone, which I agree with — but `alt` is an accessibility gap on a chart that carries real information, and the fix costs nothing while the `src` is already being rewritten. Optional, and a legitimate scope call to decline.

---

## 6. Action items

### Critical — fix before implementing

1. **Delete §4 step 3; rewrite §3.3.** Keep the `\~` escapes in `benchmarks.md:129`. They guard GFM single-tilde strikethrough, not `$` math. Removing them renders the sentence inside `<del>` with two orphaned `**`. (L1)
2. **Add a `<del>` count to validation** — 0 across all of `dist/**/*.html`, and prove it returns 1 on an escapes-removed build. This is the only check that catches item 1. (L1)
3. **Correct §2.2 and V4.** `grep -c '1000-sample cohort'` returns **1** pre-fix, not 0; the phrase is contiguous. Replace with `grep -cF '~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour'` (0 → 1, both verified). (L2)
4. **Make V1 and V3 count `katex-error` and `class="katex`, not just `application/x-tex`.** A live KaTeX ParseError on `/performance/benchmarks/` is invisible to the current checks, and any future one would be too. (L3)

### Important

5. **Restate §4 step 3 / V5's baseline** as post-fix-with-escapes vs post-fix-without-escapes. "Byte-identical apart from the `$`" is unachievable and would read as a failed fix. (L4)
6. **Document §2.2's full damage on the benchmarks page:** a second math node rendering as red `katex-error`, a literal `\~` leaking through, orphaned `**`, mis-scoped `<strong>`. (L3)
7. **Note the Pagefind index churn in V3** so "no other file changes" is not read literally. Exactly two HTML files change. (L5)
8. **Promote the assumptions with the evidence in §2 above:** A2 → verified with the five-site table (add `index.mdx:21`, and state that MDX inherits `markdown.remarkPlugins`); A3 → verified by inspecting the image, plus the orphan-asset fact; A4 → verified with the curl headers, and correct the README precedent to `master`; A5 → narrowed to named classes, with GFM strikethrough added as the class it missed.
9. **Reconsider §7 Q2 in favour of `master`,** or record why `dev` beats the two existing `master` precedents (`README.md:3-4`, `astro.config.mjs:38`). (§5)

### Optional

10. Use `grep -F` for V2 (`\begin{aligned}` fails as a regex) and `grep -o … | wc -l` instead of `grep -c` everywhere. (L6, L7)
11. Add the source-side `src="docs/` guard on `CHANGELOG.md` and `changelog.md`, not just built pages. (§4 item 5)
12. Assess `docs/public/` + `https://www.trimgalore.com/images/pigz_bench.png` for §7 Q3 — it dissolves Q2 entirely and the asset is an orphan, so moving it is free. (§5)
13. Record the `singleTilde` decision explicitly (Astro only exposes `gfm` as a deprecated boolean, so the config-layer fix costs tables/autolinks — not worth it). (§5)
14. Add an `alt` to the `<img>` while rewriting its `src`. (§5)
15. Note in A1 that `// @ts-check` provides no protection against a misspelled plugin option (Astro types them `any`; no `astro check` in CI) — this is the concrete justification for V1's built-output assertion. (A1)

---

## 7. What the plan gets right

Worth stating, because the defects above are all in the validation and one instructed step, not in the diagnosis:

- Both bugs are real and I reproduced both in built output.
- `singleDollarTextMath: false` is the correct fix, supported by the installed version, and structurally cannot touch `$$` display math. Verified by a real build.
- No content anywhere relies on single-dollar inline math, so the fix has no legitimate casualty.
- The image is genuinely broken in both renderers, the identified replacement asset is the right one, and the proposed URL works.
- The instinct that drives the whole plan — assert against built output because both failure modes leave a green build — is right, and §9's self-review names the correct top risk. The irony is that the plan's own validation table falls to exactly that failure mode in three places.
