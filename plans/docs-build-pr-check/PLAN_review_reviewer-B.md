# PLAN review — Reviewer B

**Plan:** `plans/docs-build-pr-check/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Reviewer:** B (independent; no shared state with Reviewer A)

Every finding below is labelled **[verified]** (I read the file or ran the command, output quoted) or
**[inferred]** (reasoning from documented behaviour, not executed here).

All builds were run in an isolated copy at `$TMPDIR/docsB` (rsync of `docs/` excluding `dist/` and
`.astro/`). `git status --porcelain docs .github` is empty — I did not touch the repo tree.

---

## Verdict

The change is worth making and the hosting decision is right. Two things in §8 do not hold as
written, and one of them means the plan's own validation table cannot tell a working assertion from
a silent no-op. Four factual claims are wrong; three of them do not change any conclusion, one does
(CHANGELOG). Fix the §8 items before implementing; everything else is bookkeeping.

---

## 1. Logic review

### 1.1 Core premise (§2.1) — confirmed correct

**[verified]** `.github/workflows/` contains exactly three files: `ci.yml`, `docs.yml`,
`release.yml`. `docs.yml:3-10` triggers on `push` to `[master, dev]` + `workflow_dispatch` only —
no `pull_request`.

```
$ grep -nE "npm|setup-node|astro|docs/" .github/workflows/ci.yml
(no output)
$ grep -nE "npm|node|docs/|astro" .github/workflows/release.yml
(no output)
```

Zero references. **No workflow other than `docs.yml` touches the docs site, and nothing builds it on
a PR.** The gap is exactly as described.

The eight checks listed in §2.1 also reconcile against the file: seven job definitions
(`rust-tests` `ci.yml:18`, `reproducibility` `:89`, `lint` `:134`, `audit` `:157`, `coverage` `:202`,
`validation` `:243`, `validation-ubam` `:832`) with a 2-entry OS matrix on `rust-tests`
(`ci.yml:27-28`) = 8 checks. V5's "9 checks" after the change is right.

### 1.2 CRITICAL — V2's deliberate break does not break the build

§8 V2 proposes "a malformed frontmatter key in one page". **[verified]** — it does not fail:

```
$ # install.md frontmatter: title, bogus_key: yes, description
$ npm run build; echo $?
[build] 31 page(s) built in 7.28s
[build] Complete!
0
```

Starlight's content schema is not strict about unknown keys. Running V2 as written yields a **green**
check, from which the only available conclusions are both wrong ("the gate can't fail" → panic, or
"nothing to worry about" → the gate is never actually exercised).

Two breaks I confirmed **do** fail:

```
$ # install.md with `title:` removed
$ npm run build; echo $?
[InvalidContentEntryDataError] docs → install data does not match collection schema.
  title**: **title: Required
1

$ # astro.config.mjs: { label: 'Installation', slug: "does-not-exist" }
$ npm run build; echo $?
[ERROR] AstroUserError: The slug `"does-not-exist"` specified in the Starlight
        sidebar config does not exist.
1
```

Replace V2's example with **removing a required `title:`** (tests the content-schema path) or a
**bad sidebar slug** (tests the config path, and is the more realistic real-world break — a page
rename that misses `astro.config.mjs`).

### 1.3 CRITICAL — the assertion can be a silent no-op, and §8 cannot detect it

Two facts collide.

**[verified]** §4 item 2 sets `defaults.run.working-directory: docs`, so every `run:` step — including
the assertion step — executes with cwd `docs/`. Paths in the assertion must therefore be
`dist/og`, **not** `docs/dist/og`. But §3 item 3, §7 A2, V3 and V4 all write `docs/dist/…`.
(`docs.yml` genuinely mixes both conventions: `docs.yml:40,43` are cwd-relative while
`docs.yml:48`'s `path: docs/dist` is repo-root-relative, because `with:` inputs ignore
`working-directory`. Easy to copy the wrong one.)

**[verified]** a wrong path fails *silently* in the natural pipeline form, because GitHub's default
shell is `bash -e {0}` without `pipefail`:

```
$ bash -e -c 'find dist/NOPE -name "*.png" | wc -l'; echo "exit=$?"
find: dist/NOPE: No such file or directory
       0
exit=0
```

So `find docs/dist/og -name '*.png' | wc -l` run from `docs/` prints an error to stderr, exits 0, and
**the job stays green while asserting nothing.**

Worse, §8 cannot catch this. V2 exercises the *build* step, not the assertion. V3 and V4 are
rehearsed locally with `docs/`-prefixed paths — i.e. from the repo root, a different cwd from the
job — so both can pass while the workflow step is inert.

Mitigations, all cheap:

1. State the paths in §3/§7/V3/V4 as `dist/og` and note the cwd once.
2. Make the first line of the assertion `test -d dist/og` so a wrong path is a loud failure rather
   than a zero count. (A bare `-d` test cannot be pipeline-swallowed.)
3. Rehearse V3/V4 **from `docs/`**, matching the job's cwd.

Also worth fixing in the same step: the classic inverted-grep form is a live trap here.
**[verified]**

```
$ bash -e -c 'find dist -name "*.png" -size 0 | grep -q .; echo reached end'; echo "exit=$?"
exit=1          # no empty PNGs → grep exits 1 → step FAILS on a healthy build
```

The correct shape (**[verified]** passes clean, fails on an empty file):

```bash
test -d dist/og
png=$(find dist/og -name '*.png' | wc -l)
html=$(find dist -name '*.html' | wc -l)
test "$png"  -ge 1
test "$html" -ge 1
empty=$(find dist/og -name '*.png' -size 0)
if [ -n "$empty" ]; then echo "empty OG PNG: $empty"; exit 1; fi
```

`-size 0` is portable for "empty" — **[verified]** on BSD `find` (macOS, where V3/V4 are rehearsed) a
1-byte and a 600-byte file are both *not* matched, only the 0-byte one; GNU `find` on the runner
behaves the same. `-empty` is equivalent and clearer. Avoid `stat`: `stat -f%z` (BSD) vs `stat -c%s`
(GNU) would make the local rehearsal and CI diverge.

### 1.4 §2.2's hosting argument — right answer, wrong lead reason

The three grounds for putting the job in `ci.yml` are not equally strong.

- **Permissions** (the plan's first and most emphasised reason) is the *weakest*. **[verified]**
  `docs.yml:12-15` does grant `pages: write` / `id-token: write` at workflow level. But
  **[inferred, documented GitHub behaviour]**: job-level `permissions:` overrides workflow-level, so
  `permissions: {contents: read}` on the build job removes the concern in one line; fork PRs get a
  read-only token regardless; and Dependabot `pull_request` runs get a read-only token with no
  secrets — meaning the motivating case (#368, a Dependabot PR) would never have had elevated
  permissions in the first place. Real only for the maintainer's own same-repo branch PRs, and
  trivially fixable.
- **Concurrency** is the decisive reason, and the plan understates it. **[verified]**
  `docs.yml:17-19` is `group: pages, cancel-in-progress: false`. **[inferred, documented behaviour]**
  with a static group and `cancel-in-progress: false`, GitHub keeps at most one *pending* run per
  group: a newly queued run **cancels the older pending one**. So PR docs builds would not merely
  queue behind deploys — a second PR pushed while one is queued would show a *cancelled* check.
  That is a worse failure mode than the plan's "serialising PR feedback".
- **Deploy guard** is real and correctly stated. `docs.yml:50-59` — `deploy` has `needs: build` and
  `environment: github-pages`, which carries a branch allowlist (cf.
  `feedback_branch_rename_hidden_settings`). A PR-triggered `deploy` would either be blocked by the
  environment or need a guard whose condition is a new thing to get wrong.

**Recommendation:** reorder §2.2 to lead with concurrency, and demote permissions to "avoidable but
another line to get right". The conclusion does not change; the reasoning should not rest on its
weakest leg.

### 1.5 §2.3 — first claim true, second moot

**[verified]** `ci.yml:3-15` has no `paths:` on either trigger; no job uses one. True.

**[verified]** — and this is the substantive part — **neither trunk has branch protection at all**:

```
$ gh api repos/FelixKrueger/TrimGalore/branches/dev/protection
{"message":"Branch not protected","status":"404"}
$ gh api repos/FelixKrueger/TrimGalore/branches/master/protection
{"message":"Branch not protected","status":"404"}
```

Three consequences:

1. §9 item 1 (open question, "flag it to the maintainer after merge") **resolves to a no-op** — there
   is no required-checks list to add to. Close the question with the command above.
2. §2.3's "complicates branch protection" justification for always-run is moot. Always-run is still
   the right call (12 s of build; a skipped check reads ambiguously), but that argument should be
   dropped rather than left as load-bearing.
3. **§1's goal is overstated.** "prove the site still builds, **before merge**" — with no protection,
   a red `Docs build (Astro)` does not block the merge button. The job *surfaces* a failure; it does
   not *gate* one. That is still a large improvement over #368, but §1 and §6 should say
   "surfaces a red check before merge" unless the maintainer also enables protection.

**[verified]** minor: §2.2 says "every job carries `if: github.event_name != 'schedule'`". `audit`
(`ci.yml:157-159`) deliberately does not — that is the whole point of the weekly cron
(`ci.yml:11-15`). Jobs do not *inherit* the guard; each declares it. Since §4 item 1 places
`docs-build` immediately **before `audit`**, an implementer copying the adjacent job's shape copies
the one job without a guard. V7 covers it, but the sentence should be corrected.

### 1.6 §4 item 5 / §9 item 2 — the CHANGELOG precedent is wrong

The plan asserts "prior CI-only changes such as #247 and #346 did not get entries". **[verified] —
half wrong, and it is the half that matters.**

```
$ grep -n "rustybuzz\|ttf-parser\|RUSTSEC" CHANGELOG.md
(no output)                                   # #346 → no entry. Correct.

$ sed -n '720p' CHANGELOG.md
#### Infrastructure (contributor-facing, since v2.1.0-beta.5) — CI hardening (#247)

$ sed -n '968,972p' CHANGELOG.md
#### Infrastructure (contributor-facing, no runtime effect)
- New CI gates on every PR: `cargo fmt --check` + `cargo clippy -D warnings`
  (lint), a dedicated reproducibility job that builds the release binary twice
  under a fixed `SOURCE_DATE_EPOCH` and asserts bit-identity, and a weekly
  `rustsec/audit-check` for dependency advisories.
```

#247 has a **dedicated multi-item section**. And `CHANGELOG.md:969` is the closest possible
precedent to this plan: "New CI gates on every PR", one entry per new job. The real pattern is not
"CI changes get no entry" — it is:

- audit-ignore-list tweaks (#346, #342, `d22b4b3`) → no entry;
- **new CI jobs / new PR gates → entry**, under the existing heading
  `#### Infrastructure (contributor-facing, no runtime effect)`.

`CHANGELOG.md:4` has an `### Unreleased` section, so there is a home for it. **Recommendation:** add
a one-line entry under that heading in Unreleased. Rewrite §4 item 5 and §9 item 2 accordingly —
as written they instruct the implementer to "verify" a claim that is already refuted.

### 1.7 A5 is false

§7 A5 and §10 both claim `ci.yml` jobs have no `needs:`, §10 calling it "Confirmed". **[verified]
false:**

```
$ grep -n "needs:" .github/workflows/ci.yml
247:    needs: [rust-tests, lint]      # validation
836:    needs: [rust-tests, lint]      # validation-ubam
```

The *conclusion* survives — a job with no `needs:` that nothing depends on cannot perturb the graph —
but the stated basis is wrong and §10 asserts it as verified. Restate A5 as: "no existing job will
`need` `docs-build`, and `docs-build` needs nothing, so the two `needs:` edges at `ci.yml:247` and
`:836` are untouched."

### 1.8 Assertion strength: the right half is strong, the weak half is the one that got the prose

§3.1 justifies the OG check via zero-byte PNGs. That is the *less* likely failure mode:
**[inferred]** if `satori()` or `Resvg` throws — the normal shape of a breaking 0.28→0.29 change, or
a font load failure — the static build fails outright and the build step already catches it. A
0-byte PNG requires `asPng()` to return an empty buffer without throwing: possible, narrow.

The genuinely valuable half is **"at least one"**. `docs/src/pages/og/[...route].ts:131`
(`await getCollection('docs')`) feeds `getStaticPaths()` at `:142`. If a future Astro major changes
the content-collection API such that `getCollection` returns `[]`, then `getStaticPaths` returns
`[]`, **zero OG routes are emitted, `dist/og` never exists, and every page still builds green.**
That is the real silent failure, and `≥1` catches it. §3.1 should say so — it is a stronger argument
than the zero-byte one.

Partial emission (5 of 30) still slips through. §10 flags this and dismisses it on maintenance
grounds ("a stricter check … would need updating whenever pages are added"). That is avoidable —
the count is derivable from the source, so a floor can be self-updating. **[verified]** it is exactly
1:1 today:

```
$ find src/content/docs -name '*.md' -o -name '*.mdx' | wc -l   →  30
$ find dist/og -name '*.png' | wc -l                            →  30
$ find dist -name '*.html' | wc -l                              →  32
        # 30 content pages + 404.html + dist/logos/preview.html (copied from public/)
```

So `test "$png" -ge "$(find src/content/docs -name '*.md' -o -name '*.mdx' | wc -l)"` needs no
maintenance when pages are added. Caveat: it would fire if a `draft:` or otherwise page-less entry
were ever added. If that feels brittle, a fixed floor (`-ge 25`, current 30) catches gross partial
loss and only needs touching if the site shrinks by five pages. Maintainer's call — but "≥1" and
"exact count" are not the only two options, and §10's dismissal implies they are.

Related: `test -d dist/og` (§1.3) is what makes the ≥1 check honest rather than vacuous. Note also
that "≥1 HTML page" is close to a tautology — `dist/404.html` and `dist/logos/preview.html` are both
emitted more or less unconditionally. `test -f dist/index.html` is strictly stronger for the same
cost, and is the one page whose absence unambiguously means the content collection broke.

---

## 2. Assumptions

| # | Status | Evidence |
|---|--------|----------|
| A1 | **Confirmed** | `docs.yml:39-48` — `npm ci`, `npm run build`, upload. No test/lint/link step. Reproducing the build is equivalent coverage. |
| A2 | **Confirmed with an important caveat** | See below. |
| A3 | **Confirmed** | `docs.yml:31-37` pins `node-version: 24` with the `>=22.12.0` comment. `docs/package.json` has no `engines` field, so 24 is a workflow choice, not a manifest constraint — copying it is right. |
| A4 | **Confirmed** | `docs.yml:40` already does `npm ci` on push; `ci.yml`'s `validation` job installs conda/apt packages (`ci.yml:258-261`). |
| A5 | **FALSE as stated** | §1.7. Conclusion survives, justification does not. |

**A2 caveat — the OG tree is nested, and the plan's glob only sees 3 of 30 files.** **[verified]**

```
$ find dist/og -name '*.png' | wc -l              → 30
$ find dist/og -maxdepth 1 -name '*.png' | wc -l  →  3    # index, install, quickstart
$ find dist/og -mindepth 2 | head -3
dist/og/guide/paired-end.png
dist/og/guide/reports.png
dist/og/guide/adapters.png
```

A2's literal `docs/dist/og/*.png` matches **three** files. §3 item 3's "no zero-byte PNG anywhere
under it" is correct and recursive; A2's glob and V3's `docs/dist/og/index.png` (top-level, so it
happens to work) both read as if the tree were flat. An implementer who takes A2 literally — `ls
dist/og/*.png`, or `find … -maxdepth 1` — leaves 27 of 30 PNGs, i.e. every `guide/`, `modes/`,
`rrbs/`, `performance/` and `reference/` card, unchecked. Fix A2's wording to `dist/og/**/*.png`
and say "nested, mirroring the page tree" explicitly.

**Unstated assumption, resolved in the plan's favour:** the docs build does **not** read the root
`CHANGELOG.md`. **[verified]** `docs/src/content/docs/reference/changelog.md` is a 74 KB **manual
copy** (`docs/README.md:37-40`: "Sync with the repo `CHANGELOG.md`" — copy with frontmatter
prepended); no build-time sync exists (`grep -rn CHANGELOG docs/src docs/scripts astro.config.mjs
package.json` finds only prose references). So:

- the new job needs nothing outside `docs/` — good, `actions/checkout` defaults suffice;
- a PR touching only root `CHANGELOG.md` produces an identical site, so "is it validated?" is moot;
- but **a PR that updates root `CHANGELOG.md` and forgets to sync the docs copy passes this job and
  ships a stale changelog page.** Out of the plan's stated scope, and I would not widen it — but §6
  should say so, because "the docs build is now checked on PRs" invites the assumption that it is
  covered. `docs.yml:8`'s `CHANGELOG.md` trigger path is itself vestigial for the same reason.

---

## 3. Efficiency

**§5's figures reproduce.** **[verified]** in the isolated copy:

```
$ npm ci        → added 402 packages … in 4s      (real 4.09 s)
$ npm run build → [build] 31 page(s) built in 9.63s   (real 11.77 s)
```

4 s / 12 s, exactly as claimed. **Caveat the plan should state:** these are Apple-Silicon local
figures with a warm `~/.npm`. A `ubuntu-latest` 4-core runner is typically 2–3× slower on this shape
of work — expect ~40–60 s for install+build plus ~20 s of runner setup. "Roughly one minute" holds as
an order of magnitude; the two precise numbers should be labelled as local.

**Cache behaviour — better than the plan claims, for a reason it does not mention.**
**[inferred, documented behaviour]** `setup-node`'s `cache: npm` caches `~/.npm` (the download
cache), not `node_modules` — so `npm ci` still does the full install work; only the network fetch is
saved. PR runs can read caches created on the base branch. Because §3 item 1 keeps the `push` trigger
(inherited from `ci.yml:4-5`), every push to `master`/`dev` refreshes the base-branch cache, so PRs
reliably restore. Worth one line in §5 — it is a real benefit of the always-run/push choice.

**Unmentioned cost:** on a docs push to `dev`, the site is now built **twice** — once by
`ci.yml`'s `docs-build`, once by `docs.yml`'s `build`. About one extra runner-minute per docs push.
Negligible, but §5 accounts only for PR cost.

**Wall clock:** unaffected. **[verified]** `docs-build` has no `needs:`, and nothing needs it, so it
runs concurrently; `validation`/`validation-ubam` (`ci.yml:247,836`) are the critical path.

**Free coverage worth naming in §3, since it costs nothing:** **[verified]** `npm ci` fails hard if
`package.json` and `package-lock.json` disagree — so a hand-edited dependency without a lock refresh
becomes a red check. And a bad sidebar slug fails the build (§1.2), which catches the common
page-rename-without-config-update bug. Both are real wins beyond "satori still renders".

---

## 4. Validation sufficiency

Where the job can be **green while protecting nothing**:

| Scenario | Caught? | Evidence |
|---|---|---|
| Astro/Starlight/satori bump breaks the build | **Yes** | The point of the job. |
| `getCollection` returns `[]` → zero OG cards, pages fine | **Yes**, by `≥1` — provided `test -d dist/og` is there (§1.3) | `[...route].ts:131,142` |
| satori renders 5 of 30 cards | **No** | §1.8 — self-updating floor available |
| Assertion written with the wrong cwd prefix | **No, and §8 cannot detect it** | §1.3 **[verified]** |
| Broken internal link | **No** | **[verified]**: added `[nonexistent page](/no/such/page/)` → `[build] Complete!`, exit 0. Starlight ships no link validator. |
| Root `CHANGELOG.md` edited, docs copy not synced | **No** | §2, manual copy |
| Page renamed, `astro.config.mjs` sidebar not updated | **Yes** | **[verified]** exit 1, `AstroUserError` |
| Required frontmatter field removed | **Yes** | **[verified]** exit 1, `InvalidContentEntryDataError` |
| Unknown/typo'd frontmatter key | **No** | **[verified]** exit 0 — this is V2's example |
| Empty pagefind search index / missing sitemap | **No** | Both are produced (`dist/pagefind/`, `dist/sitemap-index.xml`), neither asserted |

The first four rows are the ones that matter. Rows 1–2 are the plan's value and it delivers them.
Rows 3–4 are §8's gaps. Rows 5–6 are honest non-coverage that §6 should name so nobody assumes
otherwise — I would **not** widen scope to chase them.

**§8 table, item by item:** V1 fine. **V2 broken** (§1.2). **V3 works only if rehearsed from the
job's cwd** (§1.3) — and it is the only validation that can catch a no-op assertion, so its cwd must
be pinned. V4 fine but weak (passing on a healthy build proves nothing about the failure path; that
is V3's job). V5 fine — baseline 8 confirmed structurally. **V6 is good practice but the wrong
tool:** `yaml.safe_load` proves the file parses, not that GitHub accepts the schema — `actionlint`
(via `pipx`/`brew`) would catch an invalid `defaults` key or a bad `if:` expression that YAML happily
parses. Optional. V7 fine.

**Missing validation:** nothing proves the assertion step is not inert. Add
**V8: "temporarily `mv dist/og dist/og.bak` after the build and re-run the assertion → non-zero exit
naming the missing directory."** That is the single check that distinguishes a working gate from a
green no-op, and it is the failure mode §1.3 makes most likely.

---

## 5. Alternatives

1. **`pull_request` on `docs.yml` with `permissions: {contents: read}` on the build job** (§2.2's
   dismissed option, in its strongest form). Removes the permissions objection entirely
   (**[inferred]** job-level overrides workflow-level) and removes the duplication. Still loses on
   the `group: pages` concurrency semantics (§1.4) and the deploy guard. **Plan's rejection stands
   — but for the concurrency reason, not the permissions reason.**
2. **A third workflow file, `docs-ci.yml`, `on: pull_request` + `push`, no permissions block, no
   concurrency group.** Not considered by the plan. Fully avoids `docs.yml`'s baggage *and* avoids
   growing `ci.yml`, which is already **58 401 bytes** of entirely-Rust CI — a Node/npm job there is
   a mild cohesion cost. Costs ~6 lines of duplicated trigger config that `ci.yml` gives for free,
   and one more file to notice. **Roughly a wash; `ci.yml` is defensible on the "one PR-check
   surface" argument.** Worth one line in §2.2 so the reader sees it was weighed.
3. **Reusable workflow (`workflow_call`)** rather than the composite action §9 item 3 dismisses.
   Extract `docs.yml`'s build job into `docs-build.yml` with `on: workflow_call`, called by both.
   Eliminates the duplication that §2.2 accepts as its main cost, and is *smaller* than a composite
   action (no new action directory, no input plumbing). Trade-off: **[inferred]** a called workflow's
   jobs render as nested names in the checks list (`CI / docs-build / …`), which is uglier than a
   flat `Docs build (Astro)`, and `defaults` cannot be set at the `uses:` job level. Mention and
   reject on naming grounds, rather than only considering the heavier option.
4. **Skip the OG assertion, add `starlight-links-validator` instead.** Not recommended — link
   validation is a different (and noisier) concern, and the ≥1 OG check guards a failure mode
   nothing else can see (§1.8). Listed only because it is the obvious "what else could a docs gate
   do" alternative.

---

## 6. `CLAUDE.md` compliance

**Comment rule.** §4 item 4's "one line each" complies. Note the existing register is looser —
`docs.yml:53-56` is a four-line comment and `ci.yml:169-195` is a 27-line block — so one line is
consistent with the rule and conservative against the file. Fine as prescribed.

**Register.** No user-facing strings are added. The job name `Docs build (Astro)` is plain and
descriptive. If a CHANGELOG entry is added (§1.6), keep it to the one-line style of
`CHANGELOG.md:969`.

**"No external runtime deps"** does not apply — this is build tooling, not the shipped binary.

---

## 7. Action items

### Critical

1. **§8 V2 — replace the deliberate break.** An unknown frontmatter key builds green
   (**[verified]** exit 0). Use a removed `title:` or a bad sidebar slug (**[verified]** both exit 1).
   As written, V2 gives a false green on the plan's most important validation.
2. **§3/§7 A2/V3/V4 — resolve the cwd/path contradiction and make a wrong path loud.** With
   `working-directory: docs`, the assertion paths are `dist/…`; `docs/dist/…` inside a `find |
   wc -l` pipeline exits 0 under `bash -e` (**[verified]**), so the job would pass while checking
   nothing — and §8 as written cannot detect that. Start the step with `test -d dist/og`, use the
   `if [ -n "$empty" ]` form rather than `find … | grep -q .` (**[verified]** the grep form fails on a
   *healthy* build), and rehearse V3/V4 from `docs/`.
3. **Add V8:** `mv dist/og dist/og.bak` after a clean build, re-run the assertion, expect non-zero.
   This is the only check that proves the gate is not inert.

### Important

4. **§4 item 5 / §9 item 2 — the CHANGELOG precedent is refuted, not unverified.** #247 has a
   dedicated section (`CHANGELOG.md:720`) and `CHANGELOG.md:969` is literally "New CI gates on every
   PR" (**[verified]**). Add a one-line entry under
   `#### Infrastructure (contributor-facing, no runtime effect)` in `### Unreleased`. #346 having no
   entry is correct but is the wrong comparator — it was an ignore-list tweak, not a new gate.
5. **§7 A5 — false.** `ci.yml:247` and `:836` both carry `needs: [rust-tests, lint]`
   (**[verified]**). Restate the assumption; drop "Confirmed … no interdependencies" from §10.
6. **§1/§6/§9 item 1 — no branch protection exists on either trunk** (**[verified]** 404 "Branch not
   protected" for `dev` and `master`). Close §9 item 1 as a no-op, drop §2.3's branch-protection
   argument, and soften §1 from "prove … before merge" to "surface a failing check before merge" —
   unless the maintainer wants protection enabled, which is the separate and larger decision.
7. **§7 A2 — the OG tree is nested;** the literal `dist/og/*.png` glob matches 3 of 30
   (**[verified]**). Say `dist/og/**/*.png` and "nested, mirroring the page tree" so the assertion is
   written recursively.
8. **§2.2 — reorder the argument.** Concurrency is decisive (a newly queued run cancels the pending
   one, so PR builds could show as *cancelled*); permissions is one job-level line away from being a
   non-issue and would not even have applied to the motivating Dependabot PR. Right conclusion,
   wrong lead.

### Optional

9. **§3.1 — lead with the stronger argument:** `getCollection` returning `[]` yields zero OG cards
   with a green build (`[...route].ts:131,142`), which `≥1` catches; the zero-byte case is narrower
   because satori/resvg failures normally throw and fail the build outright.
10. **§10's residual risk is avoidable.** A floor tied to `find src/content/docs -name '*.md*' |
    wc -l` (**[verified]** 30 = 30 today) is self-updating and catches partial render, so "≥1 vs
    exact count" is a false dichotomy. Also swap "≥1 HTML page" for `test -f dist/index.html` —
    `404.html` and `logos/preview.html` are emitted near-unconditionally, so the current floor is
    nearly a tautology.
11. **§2.2 — `audit` (`ci.yml:157`) has no schedule guard** (**[verified]**); "every job carries" is
    wrong, and `audit` is precisely the neighbour §4 item 1 places the new job next to.
12. **§2.3's mechanics.** GitHub has no job-level `paths:` — path filtering would need a separate
    workflow (where required checks hang pending, the known hazard) or an in-job `git diff` guard.
    Right conclusion, imprecise premise.
13. **§6 — state the non-coverage:** broken internal links build green (**[verified]** exit 0), and
    the docs changelog page is a manual copy of root `CHANGELOG.md` (`docs/README.md:37-40`), so a
    forgotten sync also passes. Naming these prevents false confidence; do not widen scope for them.
14. **§5 — label the timings as local** (4 s / 12 s measured on Apple Silicon with a warm `~/.npm`;
    a `ubuntu-latest` runner will be 2–3× slower), note that `cache: npm` caches `~/.npm` and not
    `node_modules`, and note that the retained `push` trigger is what keeps the base-branch cache
    warm for PRs.
15. **§2.2 — mention the two alternatives not considered:** a standalone `docs-ci.yml`
    (`ci.yml` is already 58 KB of pure-Rust CI) and a `workflow_call` reusable workflow (smaller than
    the composite action §9 item 3 rejects; reject it on nested-check-naming instead).
16. **`docs/README.md:22`** describes the docs workflow triggers; a line about the new PR check would
    keep it accurate.
17. **V6** — `actionlint` would catch schema errors that `yaml.safe_load` cannot. Nice-to-have.

### Verified correct — no action

- §2.1's premise: no workflow other than `docs.yml` references npm/node/astro/docs (**[verified]**).
- §2.2's dependabot mitigation: #367 (`6df02ac`) touched `ci.yml`, `docs.yml` **and** `release.yml` in
  one grouped PR (**[verified]** `git show --stat`).
- §2.3's "no `paths` filter in `ci.yml`" (**[verified]** `ci.yml:3-15`).
- A1, A3, A4 (**[verified]**).
- §5's parallel-execution claim (**[verified]** no `needs:` on the new job).
- §4 item 1's placement between `lint` (ends `ci.yml:155`) and `audit` (starts `:157`) is clean;
  `defaults.run` already has precedent at `ci.yml:248-250`.
- The comment prescription in §4 item 4 complies with `CLAUDE.md`.

---

## 8. Reproduction notes

Isolated build environment (repo tree untouched — `git status --porcelain docs .github` empty):

```bash
mkdir -p "$TMPDIR/docsB"
rsync -a --exclude 'dist/' --exclude '.astro/' /Users/fkrueger/Github/TrimGalore/docs/ "$TMPDIR/docsB/"
cd "$TMPDIR/docsB" && npm ci && npm run build
```

Node 24 / npm from the local toolchain. Five builds run: clean baseline, unknown frontmatter key
(exit 0), missing `title:` (exit 1), bad sidebar slug (exit 1), broken internal link (exit 0).
`gh api …/branches/{dev,master}/protection` required `dangerouslyDisableSandbox` (keychain TLS
verification is blocked in the sandbox; `OSStatus -26276`).
