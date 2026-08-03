# Plan Coverage Report — docs-render-fixes

**Mode:** B (code vs. plan — no `IMPL.md`; ledger derived from `PLAN.md` §3, §4, §6, §8)
**Plan(s):** `plans/docs-render-fixes/PLAN.md`
**Date:** 2026-08-02
**Code state audited:** uncommitted working tree on `dev` @ `7201f03` — 6 files, +31/−3
**Verdict:** **COMPLETE** — 0 MISSING, 0 PARTIAL. 3 validations legitimately DEFERRED, 2 DEVIATED (1 documented, 1 undocumented-but-equivalent). One defect found in the *plan's* V8 method, not in the code.

## Summary

- Total items: **36**
- DONE: **31**
- PARTIAL: **0**
- MISSING: **0**
- DEVIATED: **2** (items 12, 13)
- DEFERRED: **3** (items 20, 21, 25 — V7, V8, V12)

### Method

Every number below was re-derived from scratch; none was taken from §10 on trust.

- **`before`** — `git archive HEAD docs` extracted to a scratch tree (provably pristine: `remarkPlugins: [remarkMath]`, `pigz_bench.png` still in `src/assets/screenshots/`, no `public/images/`), `node_modules` copied in (a symlink breaks Astro's path resolution), built with `npx astro build --outDir "$TMPDIR/.../before"`.
- **`after`** — the real working tree, built with `--outDir` to scratch.
- **`noescape`** — `after` sources with the three `\~` escapes stripped from the benchmarks line, built to scratch. This is the V9 falsifiability build.
- `docs/dist` was never written: its mtime is still `Aug 1 14:06:46 2026` after all three builds.
- The repository is unmodified by this audit: `git diff HEAD --shortstat` still reports `6 files changed, 31 insertions(+), 3 deletions(-)`.
- All three builds: 31 pages / 32 HTML files, build green.

**One harness false-negative was caught and corrected.** The first run of the CI step used `dist` as a *symlink* to a build tree; BSD `grep -r` (and even `-R`) does not descend into it, so the step reported **exit 0 / PASS against the pre-fix build**, which contains both `katex-error` and `src="docs/`. Re-run against real directories (`cp -Rl`), it fails correctly. Every grep-based result below was obtained against a real directory and paired with a demonstration that the check can fail.

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Disable single-dollar inline math, keep `$$` | §3.1 | DONE | `docs/astro.config.mjs:19` — `remarkPlugins: [[remarkMath, { singleDollarTextMath: false }]]`, exactly as specified |
| 2 | Move orphaned asset to `docs/public/`, use site-absolute URL | §3.2 | DONE | `docs/public/images/pigz_bench.png`; `src="https://www.trimgalore.com/images/pigz_bench.png"` in both files |
| 3 | Add an `alt` attribute | §3.2 | DONE | `alt="Runtime versus core count for SE and PE trimming with Python 3 and pigz"` |
| 4 | Leave the malformed `style` and space-bearing `id` alone | §3.2 | DONE | Both byte-identical in the diff |
| 5 | Keep the `\~` escapes (reversal of the superseded step) | §3.3 | DONE | 3 escapes at HEAD, 3 in the working tree — none removed |
| 6 | One-line source comment stating what the escapes guard | §3.3 | DONE | `benchmarks.md:129`, one line, and free of the trigger tokens that broke iteration #1 |
| 7 | Do **not** fix strikethrough at the config layer | §3.4 | DONE | No `gfm` key in `astro.config.mjs` (0 matches); tables still render (8 `<table>` in both builds) |
| 8 | Step 1 — config tuple | §4.1 | DONE | See item 1 |
| 9 | Step 2 — `git mv` + repoint `CHANGELOG.md` `<img>` | §4.2 | DONE | Recorded as `R` (rename), not add+delete; md5 `28882bb9…` identical to `HEAD:docs/src/assets/screenshots/pigz_bench.png`, 80,712 bytes |
| 10 | Step 3 — escapes retained + comment (inverted step) | §4.3 | DONE | Code matches the **current** text, not the superseded instruction |
| 11 | Step 4 — `changelog.md`: the `<img src>` line only | §4.4 | DONE | Diff is exactly `2 +/-` on that one line |
| 12 | Step 5 — three CI assertions | §4.5 | **DEVIATED** | All three patterns present and proven to fire, but added as a **new** step rather than into the existing `Assert build output` step, and as a `for` loop rather than three `! grep -rq` lines. Undocumented in §10. See Gaps. |
| 13 | Step 6 — CHANGELOG entry under `#### Bug fixes` | §4.6 | **DEVIATED** | Placed under `#### Fixes`. Documented as §10 deviation 1 and factually justified: `#### Bug fixes` occurs **0** times in `CHANGELOG.md` |
| 14 | V1 — no accidental math survives | §8 | DONE | `class="katex` changelog **3 → 0**, benchmarks **4 → 0** |
| 15 | V2 — `katex-error` gone site-wide | §8 | DONE | **1 → 0**; pre-fix hit was `performance/benchmarks/index.html` |
| 16 | V3 — display math still works | §8 | DONE | clumpy `begin{aligned}` **1 → 1**; `application/x-tex` **1 → 1** |
| 17 | V4 — corrupted phrase restored | §8 | DONE | Full phrase **0 → 1**; a one-character-altered pattern returns 0 on the `after` build, so the check can fail |
| 18 | V5 — benchmarks passages intact | §8 | DONE | Whole sentence present, exactly 3 correct `<strong>` spans, leaked `\~` **3 → 0** |
| 19 | V6 — `$` renders literally where eaten | §8 | DONE | benchmarks `class="mord` **46 → 0**; `$0.05/vCPU-hour` literal **0 → 1** |
| 20 | V7 — image resolves over HTTP | §8 | **DEFERRED** | Genuine: needs the deploy. The verifiable half — no `src="docs/` in any built page — is **PASS** (**1 → 0**) |
| 21 | V8 — image renders on github.com | §8 | **DEFERRED** | Deferral genuine (nothing pushed), but V8's *stated method* cannot pass as written. See Gaps. |
| 22 | V9 — no strikethrough anywhere | §8 | DONE | `<del>` **0 → 0**, and **proven falsifiable**: the `noescape` build returns **1** |
| 23 | V10 — only the intended pages changed | §8 | DONE | Exactly **2** HTML files differ (benchmarks, changelog); `Only in after: images`; **13** pagefind entries churn; `index.html` byte-identical |
| 24 | V11 — source-side guard | §8 | DONE | `src="docs/` in source **2 → 0** (1 per file) |
| 25 | V12 — CI enforces all three classes | §8 | **DEFERRED** | Genuine: needs a PR. Every local claim independently re-derived — see Validation verification |
| 26 | A1 — `remark-math` supports the option (4 layers) | §6 | DONE | `remark-math` **6.0.0**, exactly **1** copy installed; `micromark-extension-math/lib/math-text.js:18` reads it and `:79` enforces `if (sizeOpen < 2 && !single)`; Astro tuple form at `load-plugins.js:12`; `syntax.js:21` `[36]: mathFlow` (no options) vs `:24` `[36]: mathText(options)` — `$$` structurally immune |
| 27 | A2 — nothing relies on single-dollar math | §6 | DONE | Latent site confirmed: `docs/src/content/docs/index.mdx:21` — `saves time, $, CO₂` |
| 28 | A3 — the asset is the intended image, and an orphan | §6 | DONE | md5 matches HEAD; **0** `pigz` hits among the **22** optimised `_astro` images |
| 29 | A4 — no CSP anywhere | §6 | DONE | No `_headers` / `netlify.toml` / `vercel.json` / `staticwebapp.config.json`; **0** files matching `content-security-policy` in `docs/src`, `astro.config.mjs`, `docs/public` |
| 30 | A5 — the raw-HTML and relative-target classes | §6 | DONE | Exactly **1** `<img` in `CHANGELOG.md`; **0** relative `src="`; **0** relative markdown images |
| 31 | A6 — GFM single-tilde is a live fourth class | §6 | DONE | Directly demonstrated: `before` **0** `<del>`, `noescape` **1** |
| 32 | A7 — `$$` in prose still becomes inline math after the fix | §6 | DONE | Re-derived through the real chain: `singleDollarTextMath:false` → single dollars produce no math node, but `$$ … $$` in prose still yields `inlineMath("and later")` |
| 33 | §10 deviation 1 — `#### Fixes` | §10 | DONE | Documented and accurate |
| 34 | §10 deviation 2 — three escapes, not two | §10 | DONE | Documented; confirmed `\~$0.05`, `\~$41`, `\~$7`, all retained |
| 35 | §10 deviation 3 — new entries not mirrored | §10 | DONE | Documented, and **sanctioned by §4 step 4**, which scopes the `changelog.md` edit to the `<img src>` line alone. Lag claim re-derived: the mirror is missing **211** lines, of which 14 are the new entries → **197** pre-existing (§10 said 196) |
| 36 | §10 deviation 4 — blanket `<del>` ban | §10 | DONE | Documented; this is **exactly** what §4 step 5 specified, so it is a trade-off note rather than a deviation |

## Gaps (detail)

### Item 12: Step 5 — the three CI assertions (DEVIATED, undocumented)

**Expected (§4 step 5):** "add three assertions **to the existing `Assert build output` step**", given as:

```bash
! grep -rq 'katex-error' dist --include='*.html'
! grep -rq 'src="docs/'  dist --include='*.html'
! grep -rq '<del>'       dist --include='*.html'
```

**Found:** a separate step, `.github/workflows/ci.yml:196-206`, named `Assert no silent render corruption`, iterating the same three patterns in a `for` loop and printing the offending files before `exit 1`.

**Gap:** form only, and it is not listed among §10's four deviations. Enforcement is not weakened — I extracted the step body from the parsed YAML and ran it verbatim under `bash -e` against three real build trees: PASS on `after`, FAIL on `before` (`katex-error`), FAIL on `noescape` (`<del>`). Each of the three patterns was additionally proven to fire in isolation, including `src="docs/`, which the loop's short-circuit skips on the `before` build. The step sits in the same `docs-build` job, after `Build site`, and inherits the job-level `working-directory: docs`, so it greps the same `dist` the plan intended. No action required beyond recording the deviation if §10 is meant to be exhaustive.

### Item 21: V8 — image renders on github.com (DEFERRED, and the method as written cannot pass)

**Expected (§8 V8):** "View `CHANGELOG.md` on github.com **on the branch before merging**. Expected: Image renders in the GitHub view."

**Found:** the deferral is genuine — no commit, no local branch `docs/fix-katex-and-image-path`, nothing matching `katex`/`image-path` on `origin`, `dev` level with `origin/dev`, so there is nothing to view on github.com.

**Gap:** V8's stated method is unachievable given §3.2's chosen approach. The `src` is now `https://www.trimgalore.com/images/pigz_bench.png` — a site-absolute URL, not a `raw.githubusercontent.com` one. GitHub will fetch that URL, and it 404s until the site actually serves `/images/pigz_bench.png`, which happens on deploy from the deploying branch. So viewing `CHANGELOG.md` on the branch *before* merging will show a **broken image**, and V8 run as specified will fail even though the implementation is correct. §3.2 already concedes the underlying fact ("a brief window between merge and deploy where the URL 404s... Accepted") but V8 was not re-sequenced to match. V7 and V8 have therefore collapsed into a single post-deploy check. This is a defect in the plan's validation method, not in the code — no code change is implied.

## Validation verification

Every row re-derived independently; §10's reported figures reproduced exactly except where noted.

| # | Check | §10 claimed | Re-derived | Falsifiability proven | Status |
|---|---|---|---|---|---|
| V1 | `class="katex` on changelog / benchmarks | 3→0 / 4→0 | **3→0 / 4→0** | pre-fix build returns non-zero | PASS |
| V2 | `katex-error` site-wide | 1→0 | **1→0** | pre-fix build returns 1 | PASS |
| V3 | clumpy `begin{aligned}` (`-F`) / `application/x-tex` | 1→1 / 1 | **1→1 / 1→1** | n/a (survival check) | PASS |
| V4 | full corrupted phrase on built changelog (`-F -f`) | 0→1 | **0→1** | altered phrase → 0 on `after` | PASS |
| V5 | sentence intact, 3 `<strong>`, leaked `\~` | intact, 0 leaked | **3 `<strong>`, `\~` 3→0**, whole sentence present | pre-fix render is a KaTeX tree | PASS |
| V6 | benchmarks `class="mord`; `$0.05/vCPU-hour` literal | 46→0, literal | **46→0**, literal **0→1** | pre-fix build returns 46 | PASS |
| V7 | image over HTTP / no `src="docs/` in built pages | deferred | HTTP **deferred**; built-page half **1→0** | pre-fix build returns 1 | DEFERRED (half PASS) |
| V8 | image renders on github.com | deferred | **deferred**; method defective as written | n/a | DEFERRED |
| V9 | `<del>` site-wide | 0→0, returns 1 without escapes | **0→0**, and `noescape` build returns **1** | **yes — reproduced live** | PASS |
| V10 | `diff -rq before after` | 2 HTML, images added, 13 pagefind, index unchanged | **2 HTML** (benchmarks, changelog), `Only in after: images`, **13** pagefind entries, `index.html` **identical** | `diff` shown to report a difference on a mismatched pair | PASS |
| V11 | `src="docs/` in source | 2→0 | **2→0** (1 per file) | HEAD returns 1+1 | PASS |
| V12 | `docs-build` green with the new step | deferred; verified locally | **deferred**; YAML parses, step present in `docs-build`, verbatim body PASSes on `after` and FAILs on `before` and `noescape` | all three patterns proven to fire | DEFERRED (local PASS) |

### The V9 reproduction, in full

With the config fix applied and only the escapes removed, the `noescape` build gives `<del>` = **1**, `katex-error` = **0**, `class="katex"` on benchmarks = **0**, and renders:

```
On AWS at <del>$0.05/vCPU-hour, trimming 84M PE reads at the nf-core default costs
roughly <strong>$0.041 with TG</strong> vs <strong>$0.007 with the Rust v2 build</strong>
— a 5.9× saving per sample. Across a 1000-sample cohort that scales to **</del>$41 with
TG vs ~$7 with Rust v2**, …
```

versus the `after` build, which renders the sentence in full with three correct `<strong>` spans and no `<del>`. This is the §3.3 regression exactly as both reviewers described it: every math-based check passes while the page ships a struck-through sentence and a collapsed `**` pair. Item 5 (keep the escapes) is therefore confirmed load-bearing, and the inverted step 3 is confirmed to match the current plan text rather than the superseded instruction.

Also confirmed from §10's iteration log #1: the guard comment **does** reach the built HTML (`<!-- Tilde escapes below guard GFM single-tilde strikethrough, not dollar maths. … -->` is present in `after/performance/benchmarks/index.html`), and the shipped wording contains **0** occurrences of `\~` and **0** of `<del>` — so the rewording that fixed iteration #1 is genuinely in place.

## Verdict

**COMPLETE.** Every behaviour in §3, every step in §4, and every assumption marked verified in §6 is present in the working tree and independently confirmed. Nine of the twelve validations were re-derived and pass with their falsifiability demonstrated; the numbers in §10 reproduce exactly.

Nothing is missing and nothing needs implementing. Three things are worth knowing before this is pushed:

1. **V7, V8 and V12 are genuinely deferred** — nothing is committed, no branch exists locally or on `origin`, and `dev` is level with `origin/dev`, so a deploy, a push and a PR are all real prerequisites. No verifiable check was folded into them to avoid running it: V7's built-output half was checked and passes, and all of V12's local claims were re-derived.

2. **V8 cannot pass as written** (see Gaps, item 21). Because §3.2 chose a site-absolute URL over a `raw.githubusercontent.com` one, the GitHub view will show a broken image until the site is deployed. V8's instruction to check "on the branch before merging" should be re-sequenced to post-deploy, at which point it and V7 are the same check. This is a plan-text correction, not a code change.

3. **Item 12 is an undocumented deviation in form.** The three CI assertions live in a new `Assert no silent render corruption` step using a loop, rather than being appended to `Assert build output` as three `! grep -rq` lines. Enforcement is equivalent — I proved all three patterns fail correctly against real defect builds — so this only matters if §10's deviation list is intended to be exhaustive.

One methodological caution for anyone repeating this audit on macOS: BSD `grep -r` will not descend into a symlinked `dist`, and returns exit 0 with zero hits on a tree that plainly contains the pattern. That produced a false PASS of the CI step against the pre-fix build here before being caught. Use real directories.
