# CODE review — `plans/docs-render-fixes` — Reviewer B

**Target:** uncommitted working tree at `/Users/fkrueger/Github/TrimGalore`, branch `dev` (`7201f03`), +31/−3 across six files.
**Plan:** `plans/docs-render-fixes/PLAN.md` (read in full, including §10).
**Method:** real Astro build to `$TMPDIR/revB-after` (never `docs/dist`); the production remark/rehype chain driven directly under both math settings; every assertion proved falsifiable on planted sentinels before being trusted.
**Modifications made:** none to any tracked file. One accident, disclosed and reverted — see "Experimental modifications" at the end.

---

## Summary

The two user-visible fixes are correct, and I verified them independently rather than taking §10 on trust. The config change does exactly what it claims at every layer I could measure; the asset move is byte-identical and lands at the URL the tag names; the `alt` text is accurate against the actual image; the retained `\~` escapes and the reworded guard comment are clean.

The CI step's shell is also correct — `exit 1` does propagate out of the `for` loop, all three quoted patterns survive parsing, and `--include` after the operand works on both grep flavours. I proved each of the three patterns fires independently.

Where the change falls short is not in the shell but in **what the shell is pointed at and when it runs**:

1. **It never runs on the path that publishes the site.** `docs.yml` builds and deploys the same site on every push to `dev` and has no assertions at all. The new step lives only in `ci.yml`, which gates PRs. Both defects being fixed reached production via that unguarded path.
2. **`katex-error` cannot see the failure mode this change is mostly about.** The changelog price corruption — the headline of the plan's bug 1 and of the CHANGELOG entry — emits `class="katex"` and **zero** `katex-error`. The check fires today only because one unrelated construct on `benchmarks.md` happens to break KaTeX's parser.

Both are small additive fixes. Neither is a regression — the change makes nothing worse — but together they mean the plan's claim that this step "prevents *recurrence*" (§4 step 5) is not yet true.

One factual error in a committed comment (MEDIUM-2) and one false-positive class that the sibling automation will make live (MEDIUM-1) round out the substantive findings.

---

## Verified correct

Recording these because several are claims I would otherwise have had to take from §10.

| Claim | Evidence |
|---|---|
| Asset move preserved bytes | `git rev-parse HEAD:docs/src/assets/screenshots/pigz_bench.png` and `git rev-parse :docs/public/images/pigz_bench.png` both `06baf8c2635e887e295ebb2350ba543a90b3b0ca`; sha256 `875e00d7ec97…` on both sides; 80,712 bytes; recorded as `R` in `git status --short` |
| Build emits it at the URL's path | `$TMPDIR/revB-after/images/pigz_bench.png`, sha256 `875e00d7ec97…` — identical. So `https://www.trimgalore.com/images/pigz_bench.png` will resolve once deployed |
| `alt` text is accurate | Opened the PNG: two-panel bar chart, *"SE trimming with Python3 + pigz (~17.5M reads, 2x75bp)"* and the PE equivalent, `time [mm:ss]` vs `# of cores` (1,2,3,4,6,8). "Runtime versus core count for SE and PE trimming with Python 3 and pigz" describes it correctly |
| `docs/Images/` never existed | `git ls-files \| grep 'docs/Images'` (case-sensitive) → empty. CHANGELOG's "exists in neither the repository nor the built site" holds |
| Config fix works, display maths survives | Real build: `class="katex` present in exactly one file, `./performance/clumpy/index.html`. Benchmarks and changelog: 0 |
| §10's V1 numbers are right | Production plugin chain, `singleDollarTextMath` on → off: `benchmarks.md` 4→0, `changelog.md` 3→0, `clumpy.md` 4→4 |
| CHANGELOG entry 1's rendering claim is exact | Visible-text extraction, pre-fix: `— ~7vs 7 vs ~7vs 41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. Post-fix: `— ~$7 vs ~$41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. Two `$` deleted, `7vs` as a maths variable — as written |
| V4 / V5 hold on the real build | Full phrase count 1; three `<strong>` spans `['$0.041 with TG', '$0.007 with the Rust v2 build', '~$41 with TG vs ~$7 with Rust v2']`; 0 leaked `\~` |
| Three `\~` escapes retained | `sed -n '129,130p' … \| grep -o '\\~' \| wc -l` → 3. §10 deviation 2 accurate |
| Guard comment is clean | Contains none of `katex-error`, `src="docs/`, `<del>`; full assertion passes on the real build. Iteration #1's own-goal is genuinely fixed |
| The comment does not swallow the paragraph | Built output: `<p>On AWS at ~$0.05/vCPU-hour, …` — the HTML comment closes on its own line, paragraph renders normally |
| `exit 1` escapes the `for` loop | Ran the block verbatim under `bash -e -c` with a trailing `echo REACHED-END-OF-LOOP`. On each planted sentinel: exit **1**, sentinel echo never printed |
| All three patterns fire | Planted `katex-error`, `src="docs/Images/x.png"`, `<del>struck</del>` into a copy of the build one at a time — each produced `Found … in built output:` + the file list + exit 1 |
| Patterns survive shell parsing | Single-quoted, so the `"` in `src="docs/` and the `<`/`>` in `<del>` are literal; `"$pat"` is double-quoted so no globbing. None contain BRE metacharacters |
| `--include` after the operand works | BSD grep 2.6.0-FreeBSD: found→0, absent→1, `--include='*.css'`→1 (i.e. honoured). GNU grep on `ubuntu-latest` permutes options by documented behaviour. No portability issue |
| Root `CHANGELOG.md` is safe for the sibling automation | Whole file through the post-fix chain: `katex-error=0  class="katex"=0  <del>=0  src="docs/=0` |
| `docs/public/` precedent is real | `docs/public/logos/` already holds 6 un-optimised PNGs, so §3.2's "the pattern is already established" is correct — the asset is not the site's only unoptimised image |
| The move left no orphan directory | `docs/src/assets/screenshots/` still holds 6 other PNGs |

---

## Issues by area

### Logic

#### HIGH-1 — The assertion never runs on the path that publishes the site

`docs.yml` is the deploy workflow. It builds the same site and has **no output assertions whatsoever**:

```
$ grep -n 'Assert\|grep' .github/workflows/docs.yml
(nothing)
$ grep -n 'name:' .github/workflows/docs.yml
29:  - name: Checkout        32: - name: Setup Node       40: - name: Install dependencies
43:  - name: Build site      46: - name: Upload Pages artifact    62: - name: Deploy to GitHub Pages
```

It triggers on `push: branches: [master, dev]` with `paths: ['docs/**', 'CHANGELOG.md', …]` (`docs.yml:3-9`), and its own comment at `docs.yml:53-56` states the consequence plainly:

> `dev` is the canonical live-content source — every push that touches docs/\*\* ships to www.trimgalore.com immediately.

`ci.yml` also runs on push to `dev`, but as a **separate workflow**. The deploy has no `needs:` on it and no `workflow_run` gate. So a direct push to `dev` deploys whatever the build produced; the new step can only put a red X on the commit afterwards, or block a PR if one is opened.

Both defects this plan fixes reached the live site. The plan cites this step as "the only change that prevents *recurrence* rather than fixing one instance" (§4 step 5). On the publishing path, it prevents nothing.

`docs.yml:22` also carries an explicit in-repo instruction that this change disregards:

> `# ci.yml's docs-build job builds the same site to gate PRs; keep the steps in step.`

In fairness, the drift is **pre-existing** — the older `Assert build output` step is `ci.yml`-only too — so this widens a gap rather than opening one. But it widens it exactly where it matters.

**Fix (exact):** insert into `docs.yml`'s `build` job, between `Build site` (`:43`) and `Upload Pages artifact` (`:46`):

```yaml
      - name: Assert no silent render corruption
        run: |
          for pat in 'katex-error' 'src="docs/' '<del>'; do
            if grep -rq "$pat" dist --include='*.html'; then
              echo "Found '$pat' in built output:"
              grep -rl "$pat" dist --include='*.html'
              exit 1
            fi
          done
```

`working-directory: docs` is already the job default, so `dist` resolves identically. This blocks the artifact upload, which is the only place a block has teeth.

#### HIGH-2 — `katex-error` is blind to the silent failure mode; the check fires today only by coincidence

Driving the production chain with `singleDollarTextMath` restored to its default:

```
singleDollarTextMath ON (pre-fix):
  benchmarks.md    katex-error=1  class="katex"=4
  changelog.md     katex-error=0  class="katex"=3      <<<
  CHANGELOG.md     katex-error=0  class="katex"=3      <<<
  clumpy.md        katex-error=0  class="katex"=4
```

The changelog price corruption — the plan's bug 1, the first thing the CHANGELOG entry describes, and the defect whose visible text I reproduced above — produces **zero** `katex-error`. It renders wrongly and silently. The only reason the CI pattern catches a config revert at all is that `benchmarks.md:129`'s `\~$` sequence happens to break KaTeX's parser.

That coupling is fragile in a specific, foreseeable way: `benchmarks.md:129` is a *pricing* sentence, exactly the kind of prose someone updates. Reword it so KaTeX parses cleanly, revert `singleDollarTextMath: false`, and you get corrupted prices on two pages with a green `docs-build`.

§8 of the plan is explicit about this distinction — V1 uses `class="katex`, V2 uses `katex-error`, and V2 is described as "a standalone check, because a failed parse emits no `<annotation>`". Only V2's pattern reached CI. **V1 has no CI counterpart.**

A blanket `class="katex` ban is impossible (clumpy legitimately has 4 nodes). A page-scoped assertion is feasible **today** — on the real build, `grep -rl 'class="katex' dist --include='*.html'` returns exactly one line, `./performance/clumpy/index.html`.

**Fix (exact):** append to the new step:

```bash
          # clumpy.md's $$ display maths is the only legitimate KaTeX on the site.
          katex=$(grep -rl 'class="katex' dist --include='*.html' | sort | tr '\n' ' ')
          if [ "$katex" != "dist/performance/clumpy/index.html " ]; then
            echo "Unexpected KaTeX pages (single-dollar maths re-enabled?): $katex"
            exit 1
          fi
```

This is the check that actually pins the config change. Without it, `singleDollarTextMath: false` is unprotected against the silent case.

#### MEDIUM-1 — `src="docs/` false-positives on any page that documents the attribute, and the sibling automation makes that live

Double quotes are **not** escaped inside code elements. Confirmed in the real Astro build (Expressive Code, not just my probe):

```
$ grep -rho '<code[^>]*>[^<]*"[^<]*"[^<]*</code>' dist --include='*.html' | head -3
<code dir="auto">-a " SEQ1 -a SEQ2"</code>
<code dir="auto">--fastqc_args "..."</code>
```

Direct test — a markdown file containing inline code `` `src="docs/Images/pigz_bench.png"` ``, a fenced `html` block with `<img src="docs/Images/pigz_bench.png">`, and `~~this was removed~~` produced:

```
no match -> katex-error
MATCH   -> src="docs/
MATCH   -> <del>
```

Note the fenced block escapes `<` to `&#x3C;` but leaves `"` alone — so `src="docs/` matches while `<img` does not.

The current CHANGELOG entry escapes this **only by phrasing**: it writes `` `docs/Images/pigz_bench.png` `` without the `src="` prefix. That is why the whole-file probe returns `src="docs/=0`. Safe by luck, not by construction — and the sibling `changelog-mirror-automation` plan pipes that exact file into the build, so a future entry that quotes the full attribute (the natural way to describe an `<img src>` bug) fails `docs-build` on a docs-only change, with a message pointing at built HTML rather than the source line.

This is the same own-goal as §10 iteration #1, one layer out: a guard whose pattern collides with prose *about* the guard.

**Fix (exact):** anchor the pattern to a live attribute —

```bash
for pat in 'katex-error' '<img[^>]*src="docs/' '<del>'; do
  if grep -rqE "$pat" dist --include='*.html'; then
```

(`grep -rq` → `grep -rqE` in both invocations.) The escaped-code form `&#x3C;img src="docs/…` cannot match, so documentation of the bug stays legal while a live tag still fails.

#### MEDIUM-3 — The blanket `<del>` ban is aimed at the wrong target, and its message is thin

The regression it defends against is narrow and specific: "someone deletes the three `\~` escapes at `benchmarks.md:129`". The check it installs bans strikethrough across all 32 pages, permanently.

§10 deviation 4 documents the trade honestly and says "Revisit if that ever happens", which is a reasonable posture — so most of this is *I would have done it differently* rather than a defect. Two things are objectively improvable:

- **The message does not name the problem.** A contributor who writes `~~foo~~` in `guide/adapters.md` sees `Found '<del>' in built output:` and `dist/guide/adapters/index.html`. No mention of strikethrough, no source path, no hint that the remedy is either to drop it or to escape a tilde. §10 claims "a clear, greppable CI failure naming the file" — it names the *built* file, which is not where the fix goes.
- **A targeted alternative exists** that costs nothing and bans nothing: assert the escapes are still present at source, e.g. `test "$(grep -c '\\~' docs/src/content/docs/performance/benchmarks.md)" -ge 1`. That catches the actual §3.3 regression and leaves strikethrough legal.

**Minimum fix:** make the message say what it means —

```bash
              echo "Found '$pat' in built output:"
              grep -rl "$pat" dist --include='*.html'
              echo "If this is '<del>': GFM turns a bare '~' into strikethrough."
              echo "Escape it as '\\~' in the source .md, or remove the strikethrough."
              exit 1
```

### Errors

#### LOW-1 — The assertion fails open when `dist` is absent

`grep` exits **2** on error (missing directory); `if` treats non-zero as false; the loop completes and the step succeeds. Proved:

```
$ bash -e -c 'for pat in "katex-error" "src=\"docs/" "<del>"; do
    if grep -rq "$pat" no-such-dist --include="*.html"; then echo match; exit 1; fi; done
    echo "STEP PASSED WITH NO dist/ AT ALL"'
STEP PASSED WITH NO dist/ AT ALL
--- exit: 0 ---
```

Mitigated in practice: the preceding `Assert build output` step's `test -d dist/og` fails first. But that is a *separate* step and the coupling is implicit — reorder them, or move this step into `docs.yml` (HIGH-1) where no such predecessor exists, and the guard silently becomes a no-op.

**Fix:** make the new step's first line `test -d dist`. One line, and it makes the step self-sufficient — which matters if HIGH-1 is actioned.

### Efficiency

Non-issue, as expected. Up to 6 `grep -r` passes over 32 HTML files (a few MB) against a ~9 s build; unmeasurable. `--include='*.html'` correctly excludes `dist/og/*.png`, `dist/pagefind/*` binaries and the sitemap XML. The failure path greps twice (`-q` then `-l`) but runs once, on failure.

### Structure

#### MEDIUM-2 — The ci.yml comment states something that is false for one of the three patterns

`.github/workflows/ci.yml:197-198`:

```
# Each of these shipped to the live site once. All three fail silently:
# the build stays green and the page renders wrongly.
```

`<del>` **never shipped.** The plan is unambiguous on this:

- §2.5: "Today **zero** produce `<del>` — verified across all 32 built pages"
- §10 V9: "`<del>` site-wide **0→0**"
- §10 iteration #2 reproduced it *only* by deliberately deleting the escapes

My own build confirms 0. So the comment tells the next reader that strikethrough corruption was observed in production, when in truth it is a pre-emptive guard against a regression the plan's reviewers caught *before* implementation. That inversion matters: a maintainer reading "shipped once" will not think to question whether the blanket ban (MEDIUM-3) is proportionate.

Secondarily, this is two lines where `CLAUDE.md` defaults to one, and "All three fail silently: the build stays green and the page renders wrongly" is the reasoning-chain the convention says belongs in the commit message.

**Fix:** one line, accurate — `# Two of these shipped to the live site; all three keep the build green.`
Or drop the archaeology entirely — `# Three render defects that leave the build green.`

#### LOW-2 — Undocumented deviation: a new step instead of extending the existing one

Plan §4 step 5 says "add three assertions **to the existing `Assert build output` step**". The implementation created a new step, `Assert no silent render corruption` (`ci.yml:196`). §10's deviation list has four entries and does not record this.

I think the new step is the better choice — clearer attribution in the Actions UI, and it keeps the satori-specific comment on the old step from acquiring unrelated concerns. But §10 exists to catch exactly this, and a plan-manager coverage pass reading step 5 literally will score it as partially unimplemented.

**Fix:** add it as deviation 5 in §10.

#### LOW-3 — The guard comment publishes to production HTML and exceeds what the plan specified

Verified in the built output:

```
$ grep -o 'Tilde escapes[^>]*' dist/performance/benchmarks/index.html
Tilde escapes below guard GFM single-tilde strikethrough, not dollar maths. CI asserts no strikethrough survives in the built output. --
```

Two points:

- Plan §3.3 specified "**a one-line** source comment stating what they guard — `guards GFM single-tilde strikethrough, not $ math`". The second sentence narrates the CI tooling, which `CLAUDE.md`'s "state the fact, don't narrate" convention argues against. It is also the sentence with no reader value at the point of use.
- Structurally, an HTML comment is the wrong vehicle: it is published in page source, and it remains permanently subject to the three CI patterns — which is precisely what bit iteration #1. Any future edit to it that mentions a trigger token self-fails the build.

**Fix:** trim to the first sentence. If it should not publish at all, the CommonMark link-reference form `[note]: # "Tilde escapes below guard GFM single-tilde strikethrough, not dollar maths."` emits nothing — worth a one-build confirmation against Astro before adopting.

#### LOW-5 — The check's blast radius includes hand-written static HTML

`docs/public/logos/preview.html` is copied verbatim into `dist` and is scanned by the new pass. (It is the source of the relative `src="hero-dark.svg"` values I found while auditing — they resolve correctly from `/logos/`, so no defect there; it also means the plan's A5 "relative image targets: zero" is a statement about *content*, and there is a static-file population outside that scope.) Not a problem today; worth knowing the check is not limited to rendered Markdown.

---

## Specific questions I was asked

**Is the absolute URL right?** For the root `CHANGELOG.md`, yes — GitHub needs it. For `docs/src/content/docs/reference/changelog.md`, it is unnecessary: that file is rendered only by Astro, so a root-relative `/images/pigz_bench.png` would be byte-identical in production *and* correct under `npm run dev` and any PR preview, where the absolute URL silently fetches from production instead. **But the counter-argument wins:** the sibling automation makes `changelog.md` a build artefact of `CHANGELOG.md`, at which point the two must be identical and only the absolute form works in both renderers. Knowing trade, not a defect.

One correction to the plan's *reasoning* (not the code): §3.2 justifies self-hosting as avoiding a "third-party dependency". For the **GitHub** renderer, `www.trimgalore.com` *is* the external dependency — and per `docs/public/CNAME` it is a custom domain, so the image on github.com now depends on DNS the repository does not control. Still the better call than a branch-pinned raw URL, but the argument is weaker than stated.

**The merge→deploy 404 window.** Short and self-closing: `docs.yml` triggers on `docs/**` *and* `CHANGELOG.md`, both of which this change touches, so the push that merges it also deploys the asset. The image is already broken, so §3.2's "not a regression" holds.

**§10 deviation 3 (entries not mirrored) — defensible?** Yes. Verified staleness: `CHANGELOG.md` 1,345 lines vs `changelog.md` 1,141 — the plan's "196 lines behind" is right to within frontmatter. Adding two entries to a file 200 lines behind is arbitrary, precedent exists (`plans/phred64-ubam/PLAN.md`), and the sibling plan closes the gap wholesale. There is a mild oddity — `changelog.md` receives the `<img src>` *fix* without the entry describing it — but the published page was already silent about ~200 lines of history, so this is not a new class of wrongness.

**Does `$$` display maths still work?** Yes, and structurally so, not just observably. `clumpy.md` holds 4 `class="katex"` nodes under both settings, and after the fix it is the **only** file in `dist` containing any. That invariant is what HIGH-2's proposed assertion would pin.

**Deferred validations.** V7/V8/V12 genuinely need a deploy, a pushed branch, and a PR; I could not close them. I did reduce V7 to "does GitHub Pages serve `docs/public/` verbatim" by confirming the build emits the asset at the exact path the URL names, byte-identical — and `docs/public/logos/*` already demonstrates that it does.

---

## Fixes applied

**None.** Two other agents are working this tree concurrently and it holds uncommitted work under review. Every finding above carries the exact change I would make.

### Experimental modifications — disclosed and reverted

One accident. An early build attempt used `astro build --root <docs>` from outside `docs/`, which created an untracked `/Users/fkrueger/Github/TrimGalore/.astro/` at the repository root — **not** gitignored there (only `docs/.gitignore` covers `docs/.astro/`), so it appeared as `?? .astro/` in `git status`. I relocated it out of the repository to `$TMPDIR/stray-astro-revB`. Verified afterwards: `git diff HEAD --stat` is exactly the 6 files / +31 −3 under review, and no stray untracked entries remain. No tracked file was touched at any point. All builds went to `$TMPDIR`; `docs/dist` was never written.

---

## Recommendations

### Critical
None. Nothing here is harmful, and nothing regresses existing behaviour.

### High
1. **HIGH-1 — Add the assertion to `docs.yml`'s `build` job** (exact YAML above), between `Build site` and `Upload Pages artifact`. Without it the guard never runs on the path that publishes the site, which is the path both original defects took. Also honours `docs.yml:22`'s "keep the steps in step".
2. **HIGH-2 — Add the clumpy-only KaTeX assertion** (exact bash above). `katex-error` is blind to the silent corruption that is the change's headline defect; today the check fires only because one unrelated construct breaks KaTeX's parser.

### Medium
3. **MEDIUM-2 — Correct the `ci.yml:197` comment.** `<del>` never shipped to the live site; the plan says so three times. Cheapest fix in the list and it prevents a maintainer mis-weighing MEDIUM-3.
4. **MEDIUM-1 — Anchor `src="docs/` to `<img[^>]*src="docs/` with `grep -rqE`.** Double quotes survive unescaped inside `<code>` in the real build, so a page documenting the bug trips the guard — and the sibling automation pipes `CHANGELOG.md` straight in.
5. **MEDIUM-3 — Extend the `<del>` failure message** to name strikethrough and point at the source remedy. Optionally replace the blanket ban with a source-side assertion that the three `\~` escapes are still present.

### Low
6. **LOW-1** — `test -d dist` as the step's first line; the assertion currently passes when there is no output to check.
7. **LOW-3** — Trim the `benchmarks.md:129` comment to its first sentence; consider a form that does not publish to page source.
8. **LOW-2** — Record the new-step-vs-extended-step choice as §10 deviation 5.
9. **LOW-4** — Optionally note in §10 why `changelog.md` uses the absolute rather than root-relative URL (forward compatibility with the sibling automation), since it costs local-preview fidelity.

---

## Verdict

**Don't ship as-is — ship after HIGH-1, HIGH-2 and MEDIUM-2.**

The content and config fixes are correct and I verified them independently at every layer I could reach; they could ship today on their own merits. The objection is narrow and specific: the CI step exists solely to prevent recurrence, and as written it does not run where the site is published (HIGH-1) and cannot detect the failure mode the change is mostly about (HIGH-2). Those are roughly fifteen lines of YAML and bash between the step as written and the step doing its stated job. MEDIUM-2 is a one-line comment correction that stops the next reader being actively misinformed.

Landing the whole thing now and following up immediately would also be defensible — nothing here regresses anything. But given that both original defects reached production through exactly the gap HIGH-1 describes, closing it in the same commit is the cheaper end of that trade.
