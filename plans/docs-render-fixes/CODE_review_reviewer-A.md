# CODE review — docs-render-fixes (Reviewer A)

**Target:** uncommitted working tree at `/Users/fkrueger/Github/TrimGalore`, branch `dev` (`7201f03`), +31/−3 across six files.
**Plan:** `plans/docs-render-fixes/PLAN.md` (including §10 Implementation notes).
**Date:** 2026-08-02.
**Nothing in the repository was modified.** All experiments ran in `$TMPDIR`. Verified at the end: `git diff HEAD --stat` still reports 6 files / +31/−3, no untracked files under `docs/` or `.github/`.

---

## Summary

**The two content fixes are correct and I independently reproduced every claim behind them.** The diagnosis, the `singleDollarTextMath: false` remedy, the load-bearing `\~` escapes, the asset move, and all eleven locally-checkable plan validations hold up under my own measurement. The CHANGELOG entries are factually accurate down to the specific rendered strings they quote.

**The CI step is the weak part of the change, and it fails in the exact way the change exists to prevent.** The three patterns are `katex-error`, `src="docs/` and `<del>`. `katex-error` only appears when KaTeX **fails**. On the real pre-fix build the changelog page carried **3 successful accidental-maths nodes and 0 `katex-error`** — so on a build whose only defect is the changelog corruption this commit is fixing, the new assertion **exits 0**. I built that tree and ran it (evidence E1). The plan itself identified this blindness and split V1 (`class="katex"`) from V2 (`katex-error`) precisely because neither sentinel sees the other's defect — but only V2's sentinel was carried into `ci.yml`.

Two further gaps: the assertion gates **pull requests but not deployment** (`docs.yml` deploys on every push touching `docs/**` and has no assertion step, despite its own comment saying "keep the steps in step"), and it **passes vacuously if `dist` is absent** (proven: exit 0).

Verdict: **ship the content fixes; amend the CI step before ship.** The amendment is additive YAML with no content risk, and I have a verified replacement below.

---

## Evidence base

Post-fix build: `npx astro build --outDir "$TMPDIR/rvA_after"` from `docs/` — 31 pages, 32 HTML files, green.
Pre-fix build: staged a copy of `docs/` in `$TMPDIR/rvA_pre` with `astro.config.mjs`, `benchmarks.md` and `changelog.md` restored from `git show HEAD:…`, real (not symlinked) `node_modules`, built to `$TMPDIR/rvA_before` — 31 pages, green.
Plugin-chain probes: standalone `unified` + `remark-parse` + `remark-gfm` + `remark-math` + `remark-rehype` + `rehype-katex` runs against `docs/node_modules` (remark-math 6.0.0, single installed copy).

| Measurement | pre-fix | post-fix |
|---|---|---|
| `class="katex` — `reference/changelog` | 3 | 0 |
| `class="katex` — `performance/benchmarks` | 4 | 0 |
| `class="katex` — `performance/clumpy` | 4 | **4** (display maths preserved) |
| `katex-error` site-wide | 1 (benchmarks only) | 0 |
| `<del>` site-wide | **0** | 0 |
| `src="docs/` site-wide | 1 (changelog only) | 0 |
| literal `\~` leaking into benchmarks HTML | 3 | 0 |
| `<strong>` on the benchmarks sentence | 1 (mis-scoped) | 3 (correct) |
| HTML files differing before→after | — | exactly 2 (`benchmarks`, `changelog`) |

Whole-tree diff done with `cmp` over real temp files (process substitution is unreliable here) and **proven falsifiable** — appending `<!-- SENTINEL -->` to a copy of `index.html` was reported as a difference; `index.html` itself is byte-identical before/after. No non-HTML file differs except the 6 predicted pagefind churn artefacts. (Caveat: my "before" tree was tarred from the current worktree, so it already contained `docs/public/images/pigz_bench.png`; that asset therefore shows as present in both trees rather than "added". It does not affect any HTML comparison.)

---

## Issues by area

### Logic

**A-1 (High) — the `katex-error` sentinel cannot detect the defect this commit fixes.**

`.github/workflows/ci.yml:196-205`. `rehype-katex` emits `<span class="katex-error" …>` **only on a parse failure**; a successful parse emits `class="katex"` and no error span. I confirmed the two are disjoint: a forced KaTeX failure produces `class="katex-error"` with `contains class="katex" (non-error): false`.

The pre-fix changelog page is the pure "success" case — 3 `class="katex"`, 0 `katex-error` — and it shipped `— ~7vs 7 vs ~7vs 41 per 1000-sample cohort at AWS $0.05/vCPU-hour` to the live site. On the current build it happens to be caught by the *unrelated* `src="docs/` pattern, because the image bug lives on the same page. Remove that coincidence and the check goes green.

**E1 — isolated demonstration.** Copied the pre-fix `dist`, rewrote the changelog page's `src="docs/…` to the fixed URL, and deleted the benchmarks page — leaving a tree whose only defect is the changelog maths corruption (3 accidental nodes, 0 `katex-error`). Ran the extracted `run:` block verbatim under `bash -e`:

```
=== CURRENT ci.yml assertion ===
  -> exit=0   <-- the live bug ships green
```

The same tree under my proposed assertion exits 1 and names the file.

This also leaves plan assumption **A7 completely unguarded**. I confirmed A7 is real: post-fix, `The shell PID idiom is $$ and costs $$5 per run.` yields `katex=3 err=0` and renders as `The shell PID idiom is andcosts5 per run.` — silent corruption, zero `katex-error`, and the sibling automation pipes `CHANGELOG.md` through this pipeline unattended. The plan grades A7 "Low… nothing watches for it"; adding one sentinel makes it watched for free.

**A-2 (Medium) — the guard gates PRs, not deployment.**

`docs.yml` triggers on `push` to `master`/`dev` filtered on `docs/**`, `CHANGELOG.md`, `.github/workflows/docs.yml`. Its `build` job steps are `Checkout → Setup Node → Install dependencies → Build site → Upload Pages artifact` — **no assertion step at all**, and no `needs:` on `ci.yml`. Its own inline comment reads:

> `# ci.yml's docs-build job builds the same site to gate PRs; keep the steps in step.`

and the `deploy` job comment states the intended flow explicitly:

> `` # `dev` is the canonical live-content source — every push that touches docs/** ships to www.trimgalore.com immediately. ``

So a direct push to `dev` — the documented normal path for docs content — deploys unguarded. `ci.yml`'s `docs-build` runs in parallel and cannot block the deploy. The plan (§4 step 5) only ever specified `ci.yml`, so this is a gap in the plan as much as the implementation.

**A-3 (Low) — assessment of the blanket `<del>` ban (§10 deviation 4): the trade is sound, the message is not.**

I agree with the call. The accidental class is genuinely silent, and I verified intentional use is zero today — `~~` appears nowhere in `CHANGELOG.md`, `docs/src`, `docs/public` or `README.md`, and `<del>` is 0 across all 32 pre-fix pages. Blocking loudly beats corrupting silently.

Discoverability is the weak half. A contributor who writes a legitimate `~~withdrawn~~` gets exactly:

```
Found '<del>' in built output:
dist/reference/changelog/index.html
```

Nothing says strikethrough is banned, why, or what to do. Fix: one line in the failure message (below). Same argument applies less urgently to the other two patterns, which are self-explanatory.

### Errors

**A-4 (Medium) — vacuous pass when `dist` is missing.** Proven: with `dist` removed, the extracted block exits **0**. `grep` returns 2 on a missing directory, and because the call sits in an `if` condition, `set -e` cannot see it, so "could not look" is indistinguishable from "found nothing".

Currently masked in `ci.yml` by `Assert build output`'s `test -d dist/og` on the preceding step. It becomes **live** the moment A-2 is fixed by copying the step into `docs.yml`, which has no such guard. One line (`test -d dist`) closes it.

**A-5 (Low) — the step's comment is factually wrong for one of its three patterns.** `ci.yml:197`:

> `# Each of these shipped to the live site once.`

`<del>` never shipped. Measured 0 `<del>` across all 32 pages of the genuine pre-fix build; plan §2.5 says the same ("Today **zero** produce `<del>`"), and §10 iteration #2 shows it appears only when the escapes are removed. Two of three shipped; the third is prophylactic. Per `CLAUDE.md` ("state the fact"), reword.

**A-6 (Low) — fail-fast hides co-occurring defects.** Proven: a file containing all three patterns reports only `katex-error`, then exits. A contributor fixes one, pushes, waits, discovers the next. Collect-all costs nothing here (whole step runs in 0.09 s).

**A-7 (Low) — pattern precision.** `src="docs/` is the exact string that broke once, not the class. It misses `src='docs/`, `src="./docs/`, `src="Docs/` (this filesystem is case-insensitive — `git` pathspecs from a subdirectory already produced a false clean tree in this repo) and `href="docs/`. Cheap broadening: `src="[Dd]ocs/`. Similarly `<del>` would miss `<del class=…>`; not producible by the current pipeline, so noted rather than urged.

**A-8 (Low, robustness only) — `grep` invocation relies on getopt permutation and regex-vs-literal luck.** `--include='*.html'` sits *after* the `dist` operand, which works only because GNU and BSD getopt permute; under `POSIXLY_CORRECT` it becomes a file operand. All three patterns are fixed strings, so `-F` removes any BRE/ERE question, and `--` guards a future pattern beginning with `-`. I verified the current form behaves correctly on two implementations locally (BSD grep 2.6.0-FreeBSD, which is what `bash` resolves here, and ugrep 7.5.0 behind the shell's `grep` function) and that `--include` genuinely restricts — a sentinel planted in a `.json` file inside `dist` was ignored. So this is hardening, not a defect.

### Efficiency

Non-issue, as expected, and measured rather than assumed: the whole step is **0.09 s** across 238 files / 6.6 MB of `dist`. `--include` means non-HTML files are skipped rather than read. Build time is unchanged by the config edit (7.29 s pre-fix, 7.99 s post-fix — noise). No further comment.

### Structure

**A-9 (Low) — `§10` "Deviations" is incomplete.** §4 step 5 instructed adding the three assertions *to the existing `Assert build output` step*; the implementation created a new `Assert no silent render corruption` step. That is the **better** choice — failure attribution is clearer and the two steps have unrelated failure modes — but it is a divergence from the plan recorded only implicitly in the §10 table, not in the four-item Deviations list that exists to catch exactly this.

**A-10 (Low) — the guard comment ships to the public page.** Confirmed in the built HTML:

```html
<!-- Tilde escapes below guard GFM single-tilde strikethrough, not dollar maths. CI asserts no strikethrough survives in the built output. -->
<p>On AWS at ~$0.05/vCPU-hour, …
```

Harmless, and I verified the current wording is clean of all three trigger tokens (§10 iteration #1's own-goal is genuinely fixed). Two residual notes: it is an internal CI note on a published page, and it must forever avoid naming the tokens it documents — a constraint no future editor will know about. `CLAUDE.md` would put this in the commit message. Not worth changing on its own; worth knowing.

**I tested and disproved one hazard I expected to find here.** I thought a future multi-line rewrite of that comment could swallow the following paragraph into the HTML block, since there is no blank line between them. It cannot: an HTML block ends at the line containing `-->`, so one-line, two-line, and blank-line-separated variants all render `<p>Costs <strong>$41</strong> here.</p>` identically. No change needed.

**A-11 (Low) — split asset convention, unrecorded in source.** `docs/src/assets/screenshots/` still holds six tracked PNGs, all consumed through Astro's image pipeline (`guide/adapters.md:62`, `guide/quality.md:16-17`). `pigz_bench.png` now lives alone in `docs/public/images/` on the opposite convention. Correct — GitHub needs a fetchable URL, which `src/assets` cannot provide — but nothing at either location says so. The plan records it; the tree does not.

**A-12 (Low) — the docs page's own "synced" claim.** `docs/src/content/docs/reference/changelog.md:6-8` tells readers "The version below is a copy synced with the docs." It is 1141 lines against `CHANGELOG.md`'s 1345, and its `Unreleased > Fixes` section contains the `--fastqc`/uBAM entry but **not** the `-a2` fix that is already in `CHANGELOG.md`. Pre-existing, and this change makes it two entries worse while declining to sync.

On the task's specific question — **is §10 deviation 3 defensible?** Yes. The gap is pre-existing and structural, not created here; syncing only the two newest entries while `-a2` stays missing would be arbitrary; the page carries a caveat pointing at the canonical file; and the sibling automation closes the whole gap at once. It leaves the published changelog stale, not *wrong in a new way*. The one thing I would tighten is the word "synced" in that note if the automation is not imminent.

---

## Verified correct — not defects

Recorded because the plan asserts them and I re-derived each independently rather than trusting it.

- **Asset move is byte-exact.** `git diff --raw -M` shows `R100`, blob `06baf8c2635e887e295ebb2350ba543a90b3b0ca` identical on both sides and in the worktree; 80,712 bytes; `PNG 1810 x 541`. The build emits it at `dist/images/pigz_bench.png` with matching md5 (`28882bb97b2d776aec334f037efc2260`), so `https://www.trimgalore.com/images/pigz_bench.png` will resolve once deployed.
- **The old path really was broken in both renderers.** `raw.githubusercontent.com/.../{master,dev}/docs/Images/pigz_bench.png` → **404** on both branches; the real pre-move location → **200 image/png**. `git ls-files | grep -i '/Images/'` returns only the new `docs/public/images/` entry.
- **`$$` display maths survives.** `clumpy` retains 4 `class="katex"`, 1 `application/x-tex`, and `begin{aligned}` — unchanged before→after. Structurally guaranteed, per plan §3.1.
- **The `\~` escapes are load-bearing.** Independently reproduced §3.3: with the config fix applied and the escapes removed, `<del>` = 1 and `<strong>` drops 3 → 2, rendering `**$41 with TG vs ~$7 with Rust v2**` with literal asterisks. Keeping them is right.
- **`exit 1` propagates out of the `for` loop.** Under `bash -e` (GitHub Actions' default `bash -e {0}`) the extracted block exits 1 on each planted pattern and 0 on a clean tree. Quoting survives parsing: `'src="docs/'` keeps its double quote, `'<del>'` produces no redirection.
- **`working-directory` is inherited.** `docs-build` sets `defaults.run.working-directory: docs`, which applies to the new step, so the bare `dist` resolves. YAML parses cleanly and the step is registered in `docs-build`.
- **CHANGELOG entry 1 is accurate in every specific.** "two `$` deleted and `7vs` italicised" — the pre-fix render is `— ~7vs 7 vs ~7vs 41 per 1000-sample cohort at AWS $0.05/vCPU-hour`. "a whole sentence on the benchmarks page rendered in red as a KaTeX parse error" — 1 `katex-error` span, `style="color:#cc0000"`. "`$$` display maths is unaffected" — confirmed.
- **CHANGELOG entry 2 is accurate.** Path absent from repo and from `dist`; broken in both renderers (404s above); asset now in `docs/public/images/`; one absolute tag serves both.
- **The new entries are themselves immune.** Rendered through the pipeline under *both* settings: `katex 0, katex-error 0, <del> 0, src="docs/ 0`. Every `$`, `$$` and `~` in them sits inside a code span — including the `` `$$` `` at `CHANGELOG.md:92`, which is otherwise a live A7 tripwire. Good discipline, and it closes the §10 iteration #1 class of own-goal.
- **List looseness unchanged.** I expected the blank-line-separated new entries to flip the `#### Fixes` list from tight to loose and change `<p>` wrapping for every sibling entry. They do not: the existing `-a2` item already contains multiple paragraphs, so the list was already loose. No side effect.
- **`git diff HEAD --check`** — clean, no whitespace errors.
- **`docs/public/logos/preview.html`** is copied verbatim into `dist` and is therefore scanned by the new assertion. Currently clean of all three patterns. Noted only so nobody is surprised that a hand-written asset is in scope.

---

## Fixes applied

**None.** Two other agents are working against this tree and it holds uncommitted work under review. Every experiment ran under `$TMPDIR`; the pre-fix tree was reconstructed with `git show HEAD:…` into scratch rather than by stashing or editing. Final `git status`/`git diff --stat` confirm the tree is byte-for-byte as I found it.

One harness note, since it produced a wrong answer mid-review and could mislead the next person: my first test of the proposed assertion used a **symlinked** `dist` and reported a false pass on the pre-fix tree. BSD `grep -r` does not follow symlinks (`-R` does); GNU `-r` follows command-line arguments. Re-run against real directories, which is what CI has. All numbers in this report come from real directories.

---

## Recommendations

### Critical
None.

### High

**R1 — add the missing sentinel, and make the step fail non-vacuously.** This is the one change I would block on: as written, the assertion added to prevent recurrence does not detect the headline instance it was written for (A-1, evidence E1), and it also leaves A7 unguarded (A-4, A-6 folded in). Replace `ci.yml:196-205` with:

```yaml
      - name: Assert no silent render corruption
        # katex-error catches maths that failed; class="katex" catches maths
        # that parsed but should never have been maths at all.
        run: |
          test -d dist
          fail=0
          for pat in 'katex-error' 'src="[Dd]ocs/' '<del>'; do
            if grep -rl --include='*.html' -- "$pat" dist; then
              echo "::error::'$pat' in built output (see files above)"
              fail=1
            fi
          done
          # $$ display maths is legitimate on this one page only.
          stray=$(grep -rlF --include='*.html' -- 'class="katex' dist \
                  | grep -v '^dist/performance/clumpy/index.html$' || true)
          if [ -n "$stray" ]; then
            echo "::error::unexpected rendered maths — a bare \$ or \$\$ in prose?"
            echo "$stray"
            fail=1
          fi
          exit $fail
```

Verified (with `-F` on the three literals; swap to `-F` if you keep plain `src="docs/`):

| tree | current step | proposed step |
|---|---|---|
| pre-fix build | exit 1 (`katex-error` only) | exit 1, names all three problems |
| pre-fix, changelog maths only (E1) | **exit 0** | exit 1 |
| post-fix build | exit 0 | exit 0 |
| `dist` absent | **exit 0** | exit 1 |

The `grep -v` allowlist is a maintenance point: a second legitimate maths page must be added to it. That is the correct trade — an explicit allowlist beats a sentinel that cannot see the defect.

### Medium

**R2 — put the same step in `docs.yml` before `Upload Pages artifact`.** (A-2.) `docs.yml` is the workflow that actually ships to www.trimgalore.com, it has no assertion, and its own comment says to keep the steps in step. Copy R1 verbatim; the `test -d dist` line matters more there, since `docs.yml` has no `Assert build output` step to mask a missing `dist`. This is the change that turns the guard from "cannot be merged" into "cannot be deployed".

**R3 — correct the step comment.** (A-5.) `<del>` never shipped live. Suggested:

```yaml
        # Two of these shipped to the live site; all three fail silently —
        # the build stays green and the page renders wrongly.
```

### Low

**R4 — name the rule in the `<del>` failure.** (A-3.) Add after the loop, or inline in the error line:

```bash
          # <del> means a bare ~ opened GFM strikethrough. Escape it as \~;
          # intentional strikethrough is not currently supported site-wide.
```

**R5 — add the two new entries to `docs/src/content/docs/reference/changelog.md`, or soften its `:::note`.** (A-12.) Deviation 3 is defensible as-is; if the sibling automation is more than a release away, "a copy synced with the docs" is the line that misleads, not the missing entries.

**R6 — record why `pigz_bench.png` lives in `public/` while its six siblings live in `src/assets/`.** (A-11.) One line in the commit message is enough; `CLAUDE.md` prefers that to a source comment.

**R7 — log the §10 deviation.** (A-9.) The new-step choice is an improvement over §4 step 5 and should be in the Deviations list, which exists to catch undocumented divergence.

**R8 — after merge, verify the image on github.com and re-check if it is still broken.** Plan §3.2 accepts a brief 404 window; the sharper edge is that **GitHub proxies external images through camo, which caches**. If camo fetches while the deploy is still running it can cache the 404 past the deploy. Merge, let `docs.yml` finish, *then* load `CHANGELOG.md` on github.com (V8). If it is broken, it is a cache, not a bug — but do not diagnose it as the latter.

---

## Ship / don't-ship

**Content fixes: ship.** `astro.config.mjs`, the asset move, both `<img>` repoints, the retained escapes, and both CHANGELOG entries are correct, and I verified every claim behind them against real pre- and post-fix builds. The change makes the site strictly better on the two pages it targets and provably does not touch the other 30.

**CI step: don't ship as written — apply R1.** The step's own comment says these defects "fail silently: the build stays green and the page renders wrongly," and that is precisely what it does to the changelog manifestation of bug 1 (exit 0 on evidence-E1's tree). Shipping it as-is buys a guard that would not have caught the bug on the page named in the commit message, and creates the impression that the class is closed. R1 is four added lines, verified against four trees, and no content risk.

**R2 is the highest-value follow-up** if you want one thing beyond R1: the guard currently cannot stop a deploy, only a merge, on a repo whose documented flow is "every push that touches `docs/**` ships immediately."
