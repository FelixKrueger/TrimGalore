# PLAN review — Reviewer A

**Target:** `plans/changelog-mirror-automation/PLAN.md`
**Repo state:** `dev` @ `7201f03` (`ci(docs): build the docs site on pull requests (#371)`)
**Method:** every load-bearing claim re-derived from the files. Diff harness proven falsifiable before use (sentinel line appended, confirmed reported). `git` run from repo root throughout. Two empirical builds and one dev-server run executed.

**Verdict:** the plan's central thesis is correct and its riskiest design decision (delete the tracked file, rely on the build failing) is *empirically sound* — I confirmed it rather than took it on trust. The core defect is that §3.2/A5's content-compatibility audit missed a raw HTML `<img>` tag, and `docs/README.md` documents a second transformation the plan never reconciled. Secondary theme: the plan's two self-identified silent-failure modes (V1a duplicate H1, V6 empty output) are mitigated only by one-time manual checks, when both could be permanent generator assertions.

---

## 1. Logic review

### 1.1 What I confirmed (the plan is right about these)

These were checked because the plan labels them verified, and they hold:

| Claim | Evidence |
|---|---|
| §2.1 table: 1331 / 1141 lines; `7201f03` 2026-07-31 / `4f86c0f` 2026-07-25 | `wc -l`, `git log -1 -- <path>` — all four values exact |
| **A1: pure lagging copy** | Real-file `diff` of stripped bodies: **196** `<` lines, **0** `>` lines, and all five hunks are `NdM` deletion form (`5,49d4`, `82,136d36`, `145,152d44`, `154,160d45`, `162,242d46`). No `c`/`a` hunks ⇒ the page body is a **strict subsequence** of the root body. This is the strongest-verified claim in the plan |
| §2.1 "no `###` release section is missing" | `grep -c '^### '` = 49 in both bodies |
| §2.1 the missing section is `#### Infrastructure (contributor-facing)` | `####` count 40 root vs 39 page; the only heading among the 196 root-only lines is exactly that one |
| §2.1 "#369 is invisible publicly" | `CHANGELOG.md:87` carries it; the page's only `369` hits are read-ID digits (lines 929, 1122) |
| A7 / §3.0 leading H1 | `CHANGELOG.md:1` = `# Trim Galore Changelog`; `grep -c '^# '` = 1 root, 0 page |
| §3.0's two-H1 reasoning | Built the site: `dist/reference/changelog/index.html` has exactly one `<h1`, and it is `<h1 id="_top" …>Changelog` — Starlight-supplied from frontmatter `title`. The premise is correct |
| §2.2 `content.config.ts` quote | `docs/src/content.config.ts:6` matches the plan's snippet character-for-character |
| §3 item 3 slug at `astro.config.mjs:105` | `105: { label: 'Changelog', slug: 'reference/changelog' },` — exact line |
| A6 / §2.3 `docs.yml:5` and `:8` | `branches: [master, dev]` at line 5; `- 'CHANGELOG.md'` at line 8 |
| **The generator's entry point is actually invoked by both workflows** | `docs.yml:44` `run: npm run build` and `ci.yml:182` `run: npm run build`, both under `working-directory: docs`. Wiring into the `build` npm script reaches deploy *and* PR gate |
| §2.4 `docs-build` runs on PRs | `ci.yml:3-9` triggers `pull_request: branches: [master, dev]`; `docs-build` job at `ci.yml:158` |
| V3's "32+ HTML pages" | Clean build produced exactly **32** |
| A2 / §3.1 frontmatter | The proposed block is **byte-identical** to the current page's lines 1-4, so no OG/meta regression |
| §4 step 1 precedent | `docs/scripts/generate-png-logos.mjs:14` uses `dirname(fileURLToPath(import.meta.url))` — exactly the shape the plan proposes; and it is wired only to `logos`, not `build` (plan's §10 note is right) |
| §4 step 5 heading already exists | `CHANGELOG.md:226` `#### Infrastructure (contributor-facing)` inside `Unreleased` |
| A4 inbound links | `grep -rn -i changelog docs/src` (proven falsifiable via a `Starlight` control that returned 3 files) finds **no** reference outside the page itself |
| No existing sync/drift check to break | `grep -rn CHANGELOG .github/ scripts/` returns exactly one hit: `docs.yml:8` |
| Text hygiene | Both files end in `\n`; no CRLF; the 39 tab-bearing lines are all *inline* tabs (no leading tabs ⇒ no accidental indented code blocks) and none fall in the new 196 |

### 1.2 §9 Q1 — the load-bearing safety claim, empirically verified

The plan bets everything on "the build already fails if the page is missing". I tested it: backed up the page, deleted it, ran `npm run build`.

```
BUILD EXIT CODE: 1
[ERROR] AstroUserError: The slug `"reference/changelog"` specified in the
        Starlight sidebar config does not exist.
  Hint: Update the Starlight config to reference a valid entry slug …
    at linkFromInternalSidebarLinkItem …
    at getSidebar …  at generateRouteData …
```

**Confirmed.** Deciding against an extra CI assertion is sound. Two refinements fall out of *where* it failed:

- It failed while rendering **`/404.html`** — a route unrelated to the changelog. The error originates in `getSidebar`, which every page renders. So a missing page does not degrade one route; it **breaks the entire build**. That is a stronger guarantee than the plan claims, and it changes the risk classification of §6 effect 5 (below).
- The net depends on `astro.config.mjs:105` keeping the sidebar entry. If that entry is ever removed, the safety property silently evaporates. Worth one line in the generator comment.

### 1.3 A3 — also verified empirically

`astro dev` with the page removed: `/reference/changelog/` → **HTTP 500**, `/` → **HTTP 500**, both bodies containing `AstroUserError`. Server start confirmed via `localhost:4399` in its log (first attempt returned `HTTP=000` because `astro dev` hit `nice(5) failed: operation not permitted` under the sandbox — a false negative I discarded and re-ran unsandboxed rather than reporting).

A3 holds: `dev` needs the generator. But see V5 in §4 — the *expected* result is mis-stated.

### 1.4 Logic gaps

**L1 — §3.0's strip is anchored to line 1 only, and that is a silent hole.** The plan says "removes the first line if it is an H1, plus any blank lines immediately after it", and argues the conditional degrades safely to "no strip". It degrades safely only for the failure it imagined (a *reformatted* header). It does **not** degrade safely for the likelier change: a line *inserted above* the H1 — an HTML comment, a badge row, a `<!-- markdownlint-disable -->` pragma. Then line 1 is not an H1, the strip no-ops, the H1 survives at line 2, and the page ships **two H1s**. That is precisely the defect §8 calls "a correctness defect with no failure signal at all", reached by a one-line edit to a file whose format the plan declares a non-goal to control. Mitigation in §2.

**L2 — §3.3 and §3.2 still carry the superseded 199.** §2.1 goes to some length to establish 196 as the body gap, and §6/V9 use 196. But `§3.2:102` says "Only the **199** new lines are unproven" and `§3.3:108` says "Generating closes the **199**-line gap". Both should be 196: the generator *strips* the H1, so the first run adds 196 body lines and never adds the other 3. Since §10 explicitly claims this correction was made, the leftovers read as the correction being incomplete.

**L3 — §6 effect 5 understates its own failure mode, in the safe direction.** "A contributor running `astro build` directly (not via `npm run build`) gets no changelog page." They get **no site**: exit 1, no `dist/`. The README note (§4 step 6) is still worth having — the error message names a *slug*, not a missing generator, so it is obscure — but this is fail-loud, not a silent-degradation risk.

---

## 2. Assumptions

### A5 / §3.2 is wrong, and it is the most consequential finding

§3.2 states: *"**22 lines contain `<`-prefixed text**, e.g. `trim_galore --phred64 <input.bam>`, all inside fenced blocks or code spans."*

The count is reproducible — `grep -c '<[a-zA-Z]' CHANGELOG.md` = **22** (a plain `grep -c '<'` gives 24; the two extra are line 722's `matches < 100` and one other bare `<`). So the plan's filter is recoverable. **But that 22-line set contains `CHANGELOG.md:1164`:**

```html
<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="docs/Images/pigz_bench.png" >
```

I classified all 24 `<`-lines with an indentation-tolerant fence detector (my first attempt anchored `/^```/` at column 0 and mis-reported the indented fences at lines 38-42; corrected before drawing any conclusion). Line 1164 sits after **8** fence toggles — **outside any fence**, and it contains no backtick. It is raw HTML in prose.

What it actually does, from the built output:

```html
<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="docs/Images/pigz_bench.png">
```

It renders as a **live element**, not escaped text. And:

- `docs/Images/pigz_bench.png` **does not exist** — no `docs/Images/` directory in the repo.
- `astro.config.mjs:9-11` sets `site: 'https://www.trimgalore.com'`, `base: '/'`. From the page at `/reference/changelog/`, the relative `src` resolves to `https://www.trimgalore.com/reference/changelog/docs/Images/pigz_bench.png` → 404 regardless.

So the published page carries a broken image **today** (pre-existing, not introduced here). The findings against the plan are:

1. A5's "Verified for braces and `<…>`" is **false** for the `<…>` half. One of the 22 is a real HTML element.
2. §3.0's *"This is the one transformation the generator performs"* is correct only by accident — see below.
3. §3.2's MDX warning is better-founded than the plan realises, and for the wrong reason. It cites `<input.bam>` as the MDX hazard; that line is inside a fence and would be safe. The line that would actually break an MDX conversion is 1164 — `<img …>` without a self-closing slash is invalid JSX and a hard build error. Worth naming *that* line in the generator comment, since it is the one the warning exists for.
4. **V4 would not catch it.** V4 checks "no raw `{` or stray `<input.bam>` markup leaked". The `<img>` is neither.

Contrast: the plan's assessment of the **new** content is correct. The `<`-lines among the 196 root-only lines are exactly 39, 40, 41, 50 — the "4 occurrences" the plan cites. Lines 39-41 are in an indented fence; line 50's `<bam>` is inside a backtick span. All safe. Its brace claim also holds: 17 `{` lines, and `A{N}` renders literally in the built HTML. Bare `< 100` at line 722 is correctly escaped to `&#x3C;`. The audit is sound for everything except the one raw-HTML line, which happens to live in the already-published portion — which is why the drift masked it.

### A4 is incomplete: `docs/README.md` references the page

A4: *"Nothing else references the page's path except `astro.config.mjs:105` and inbound links. Grep before deleting."* I grepped. `docs/README.md:37-41`:

```markdown
## Sync with the repo `CHANGELOG.md`

`docs/src/content/docs/reference/changelog.md` is a copy of the
top-level `CHANGELOG.md` with a Starlight frontmatter block prepended
and one inline image path rewritten. Keep them in sync when releasing.
```

Three consequences:

1. **§4 step 6 is mandatory, not conditional.** The plan hedges — "Note it in `docs/README.md` *if* that file documents the build; check first". It does, and the existing text becomes actively false in three ways: the file is no longer "a copy", "keep them in sync when releasing" is the retired manual step, and the path named no longer exists in git. This is a correction, not an addition.
2. **"and one inline image path rewritten" is a second transformation the plan never reconciled.** It is now stale — my diff proves the `<img>` line is byte-identical in both files (0 page-only lines, deletion-only hunks), so no rewrite survives. The rewrite existed once and a hand-sync destroyed it. That is a neat illustration of the plan's thesis, but it also means **a verbatim generator cements the broken path permanently**, whereas the hand-maintained copy at least once had it right. The plan should decide this deliberately rather than inherit it.
3. Separately, `docs/README.md:18` already says the site "auto-deploys from the `master` branch", which `docs.yml:5` (`[master, dev]`) contradicts. Since §3.1's whole argument is that `dev` is what publishes, and the plan is editing this file anyway, fix it in the same pass.

### Assumptions that hold

- **A1** — verified above; the strongest claim in the plan.
- **A2** — `docsSchema()` requires `title`; `description` present and identical to what the plan proposes.
- **A3** — verified empirically (§1.3).
- **A6** — `docs.yml:5` confirmed. The plan is candid that the banner-link/deploy-branch coupling is unenforced. §5 proposes removing the coupling instead of documenting it.
- **A7** — verified for the current file; the degradation argument is only *partly* sound (L1).

### Unverifiable here

§2.1's *"#359's own description records the page as having been frozen at v2.1.0-beta.5"* — `gh` fails in this sandbox with `tls: failed to verify certificate: x509: OSStatus -26276`. **Not refuted, just unchecked.** It is colour, not load-bearing; §2.4's substance I verified directly from `ci.yml` instead.

### Stale figure

§3.1 and §10 both say `dev` is "34 commits ahead of `master`". `git rev-list --count master..dev` = **37**. The argument is unaffected (it needs only "nonzero"), but it is presented as a measured fact.

---

## 3. Efficiency analysis

**No concern, and nothing to optimise.** `CHANGELOG.md` is 85,525 bytes. One read + one write of 85 KB is sub-millisecond against an observed `astro build` wall time of ~7 s (`Completed in 6.33s` + `6.90s` across the two vite passes, 32 pages). Memory: one 85 KB string. Even at 10× the changelog's size this stays invisible. The plan's "30-line Node script" framing is the right complexity budget.

Two cheap improvements, both correctness-adjacent rather than performance-driven:

- **Write only if the content differs.** Read the existing file first and skip the write when identical. This makes V2's idempotency true *at the filesystem level* (mtime unchanged, not merely byte-identical), and avoids needlessly invalidating Astro's `.astro/` content-layer cache on every `dev` restart. Two lines.
- **Assert rather than trust** — see §4. The assertions are O(lines) on data already in memory; free.

One thing *not* worth doing: watching `CHANGELOG.md` for live reload in `dev` (the plan defers this). Agreed — `astro dev` restarts are cheap and the flag surface isn't worth it.

---

## 4. Validation sufficiency

V1-V10 cover the right *implementation-time* ground. The gap is that the two checks the plan itself identifies as silent-failure-prone are **one-time manual checks**, so they protect the first commit and nothing after it.

### The silent-failure path that passes every listed check

Suppose a later edit to the generator truncates its output — an off-by-one in the strip, a partial read, a regex that eats more than line 1. Then: the file exists → the sidebar slug resolves → the build exits 0 → `ci.yml:191`'s `test "$html" -ge 1` passes (31 other pages) → `docs-build` is green → the deploy ships a **half-length changelog**. V1 would have caught it, but V1 is run once by hand at implementation time. Nothing in the repo re-runs it.

**Fix (highest-value item in this review): move V1 and V1a inside the generator as assertions.** After writing, the generator asserts

- emitted body line count **==** input line count after the §3.0 strip → fails loud on any truncation, and subsumes behaviour item 5's "or empty" with a far stronger floor;
- **zero** `/^# /` lines in the emitted body → fails loud on the duplicate-H1 hole, including L1's insert-above-the-H1 case that the line-1-anchored strip cannot see.

Both are ~3 lines, run on data already in memory, and convert the plan's two acknowledged blind spots into build failures on every PR via the existing `docs-build` job — no new CI wiring, which is consistent with §9 Q1's decision to lean on the build.

### Specific defects in the listed checks

- **V9's first half is vacuous.** "Count `###` sections in the generated page vs the root file → Equal." Already equal: 49 and 49, on the *stale* page, before any change. It cannot detect the gap it is meant to verify. Only the body-line-count half discriminates. Either drop the section count or replace it with the `####` count (40 vs 39 today — that one does discriminate).
- **V5's expected result is mis-specified.** "Remove the generated file, start `dev`, request the page → 200, not 404." The observed failure is **HTTP 500 on every route**, not 404 on one (`/` returned 500 too). A tester looking for 404 sees 500 and may mis-file it. Restate as: 200 with the generator; without it, every route 500s with `AstroUserError`.
- **V4 will not catch the `<img>`.** It looks for `{` and `<input.bam>`. Add: assert the built HTML contains no `<img` with a relative `src`, or at minimum inspect the rendered `pigz_bench` element and record the decision from §2 finding 2.
- **V6 tests missing input only.** Unreadable (permissions) and partial input are untested. The line-count assertion above covers partial; a plain `try/catch` covers unreadable.
- **V1's own harness.** The plan is right to demand a falsifiability proof, and right about the mechanism. I hit the sibling trap instead — `nice(5) failed: operation not permitted` under the sandbox made a working dev server report `HTTP=000`, which read exactly like a 404-adjacent failure. Worth generalising the A1 note: *any* negative result in this repo needs the harness proven capable of a positive first, not just `diff`.
- **`working-directory: docs`** — the plan flags that nothing exercises the generator under CI's cwd. I confirmed the reasoning holds: both workflows use bare `actions/checkout` with no `sparse-checkout`, so `../CHANGELOG.md` is present, and `import.meta.url` resolution is cwd-independent. Low risk, and `plans/docs-build-pr-check/PLAN.md:54` should get the one-line correction the plan already promises.

---

## 5. Alternatives

**§5's rejection of a custom content loader is well-reasoned, and my testing strengthens it.** The plain-file approach gets a hard build failure for a missing page *for free*, from Starlight's own sidebar resolution — a loader-based design would have to reproduce that guarantee deliberately. Combined with the `docsLoader()`-composability and dependabot-churn arguments, the choice stands. Keep it.

**Better than hardcoding `dev`: pin the banner link to the built commit.** §3.1 retargets `master` → `dev` and A6 concedes the resulting coupling is unenforced ("a comment in the generator, not a check"). Instead, read `GITHUB_SHA` — set in both `docs.yml` and `ci.yml` — and emit `blob/<sha>/CHANGELOG.md`, falling back to `dev` when unset (local builds). The link then points at *exactly the bytes the page was generated from*, the master/dev question stops mattering, and A6's unenforceable invariant disappears rather than being documented. Same cost as the `dev` hardcode.

**The middle ground the plan does not consider: keep the file tracked, add a CI drift check.** Generate into the tracked path and have `docs-build` fail if the working tree is dirty afterwards — the `cargo fmt --check` pattern, already idiomatic in this repo. Trade-offs: the changelog stays browsable on GitHub under `docs/`, no untracked-file confusion (§10's "Low" risk), and no need to delete anything. Against: every `CHANGELOG.md` edit now *requires* a second committed file, which is the manual step the plan exists to remove, and the failure mode is a red PR rather than an impossible state. **I would still choose the plan's approach** — "drift unrepresentable" beats "drift detectable" — but this option belongs next to §5's rejected loader so the decision space is on the record.

**Not viable, for completeness:** a symlink from the content path to `../../../../../CHANGELOG.md`. Fails on the frontmatter requirement (§2.2) exactly as the copy does, and adds a Windows-checkout hazard. Correctly not considered.

---

## 6. Action items

### Critical

1. **Fix §3.2 / A5 and decide the `<img>` question.** `CHANGELOG.md:1164` is a raw HTML `<img>` **outside** any fence (8 fence toggles precede it), inside the plan's own 22-line set. It renders as a live element with `src="docs/Images/pigz_bench.png"`; that file does not exist and the relative path would 404 from `/reference/changelog/` under `base: '/'`. Correct A5's "Verified", and choose explicitly: rewrite the `src` in the generator (making it a *second* transformation, with §3.0's "the one transformation" reworded), or record that the image stays broken and why. Extend V4 to assert no relative-`src` `<img>` in the built HTML.

2. **Move V1 and V1a into the generator as assertions.** (a) emitted body line count == input line count after the §3.0 strip; (b) zero `/^# /` in the emitted body. Without these, a later truncation or a duplicate H1 builds, passes `docs-build` (`test "$html" -ge 1`), and deploys silently. This is the difference between validating the first commit and validating the mechanism.

3. **Rewrite `docs/README.md:37-41`; amend A4.** The section names the page path, calls it "a copy … with one inline image path rewritten", and instructs "Keep them in sync when releasing" — all false after this change, and the image-rewrite clause is already stale (my diff shows no surviving rewrite). §4 step 6's "check first" resolves to *yes*: the edit is mandatory. Add `docs/README.md` to A4's list of referencing files. While there, fix `README.md:18`'s "auto-deploys from the `master` branch" against `docs.yml:5`.

### Important

4. **Close L1: the strip is anchored to line 1.** A badge row or HTML comment inserted above the H1 makes the strip no-op and ships two H1s with no signal. Either scan for a leading H1 anywhere before the first `###`, or rely on action item 2(b) — but say which.

5. **Fix V9 and V5.** V9's `###`-count half is already equal (49/49) on the stale page and cannot detect the gap; use the `####` count (40 vs 39) or drop it. V5's expectation should be "200 with the generator; without it, *every* route 500s", not "not 404" — verified: `/` returned 500 as well.

6. **Replace the two leftover `199`s** at `§3.2:102` and `§3.3:108` with 196. §10 claims this correction was already made.

7. **Prefer a commit-pinned banner link** (`GITHUB_SHA` → `blob/<sha>/CHANGELOG.md`, fallback `dev`) over hardcoding `dev`. Removes A6's unenforced coupling instead of documenting it.

### Optional

8. "34 commits ahead of `master`" is now **37**. Argument unaffected.

9. State §3.2's count provenance — `grep -c '<[a-zA-Z]'` = 22, plain `grep -c '<'` = 24 — so the number reproduces. And "their 4 `<…>` occurrences are all inside fences" is 3 in a fence + 1 in a code span (line 50); safe either way, but the phrasing is loose.

10. **Write only if content differs**, so repeat builds don't touch mtime or invalidate the `.astro/` cache.

11. **Soften §6 effect 5.** Direct `astro build` doesn't yield "no changelog page" — it exits 1 with no site at all. Fail-loud; the README note is still justified because the error names a slug, not a missing generator.

12. Name `CHANGELOG.md:1164` in the generator's `.md`-not-`.mdx` comment. It is the line that would actually break an MDX conversion (`<img …>` without a self-closing slash is invalid JSX); the cited `<input.bam>` is fenced and safe.

13. Generalise A1's harness note beyond `diff`: I got a false `HTTP=000` from a *working* dev server because `astro dev` hit `nice(5) failed: operation not permitted` under the sandbox. Any negative needs the harness proven capable of a positive.

14. `§2.1`'s "#359 … frozen at v2.1.0-beta.5" is **unverified** — `gh` fails here with `tls: failed to verify certificate: x509: OSStatus -26276`. Not refuted.

---

## 7. Note on repo state

All experiments restored their inputs; `git status --short docs/ CHANGELOG.md` is clean at the end of this review. The page was backed up and restored at 1141 lines after each of the three destructive tests. `docs/dist/` was rebuilt (gitignored). Confirmed the case-insensitivity trap is live in this working copy — `ls -d Docs docs` resolves both — and ran all `git` commands from the repo root accordingly (`git ls-files | grep -c '^docs/'` = 154, `'^Docs/'` = 0).
