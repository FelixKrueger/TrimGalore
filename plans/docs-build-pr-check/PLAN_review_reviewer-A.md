# PLAN review — Reviewer A

**Plan:** `plans/docs-build-pr-check/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Reviewer:** A (independent; no shared state with Reviewer B)
**Method:** read every workflow file, queried the GitHub API for branch protection, and reproduced the docs build three times (clean, simulated empty-PNG, simulated renderer throw) in `$TMPDIR/revA-docs`. No repo file was modified (`git status --porcelain` clean for tracked files at the end).

---

## Verdict

**The plan is sound and its core premise is real.** I confirmed the gap it describes, confirmed the silent-failure mode the OG assertion targets (this was the claim most likely to be hand-waving, and it holds), and confirmed the hosting choice. Recommend proceeding.

Five findings need a plan edit before implementation. One is a genuine correctness bug in the assertion *as specified in A2/V3* — an implementer following the plan literally would ship a check that covers 3 of 30 OG images and would pass on the exact failure it exists to catch. I verified that miss end-to-end. The other four are factual corrections to stated premises, and one CHANGELOG precedent that is the reverse of what the plan claims.

---

## 1. Logic review

### 1.1 Core premise (§2.1) — confirmed

Read all three workflow files.

- `docs.yml:3-10` — `push` to `[master, dev]` with a `paths` filter, plus `workflow_dispatch`. **No `pull_request`.** Confirmed.
- `grep -n "pull_request" .github/workflows/*.yml` returns exactly one hit: `ci.yml:6`. `ci.yml` is the *only* workflow that runs on PRs.
- `release.yml` contains no `npm`, `node` or `docs/` reference at all (grep returned nothing). Nothing else builds the site.

So: today no PR check builds the docs. The premise is not overstated.

### 1.2 The OG-image premise (§3.1) — confirmed empirically, both halves

This was worth testing rather than believing, and it survives. In `$TMPDIR/revA-docs` I patched `src/pages/og/[...route].ts` twice:

| Simulated failure | `astro build` exit | On-disk result |
|---|---|---|
| Endpoint **throws** (`throw new Error('simulated satori failure')`) | **1** — `[ERROR] [build] Caught error rendering /og/guide/adapters.png` | build aborts, no `dist/og` |
| Renderer returns **empty buffer** (`png = Buffer.alloc(0)`) | **0** — `[build] 31 page(s) built`, `Complete!` | `dist/og/guide/adapters.png`, **0 bytes** |

That is precisely the shape §3.1 describes: hard failures are already caught by "the build didn't crash", and the residual gap is the zero-byte write, which is silent. The assertion earns its place. Note the failing file landed in a **nested** directory — see 1.3.

### 1.3 IMPORTANT — A2's path form covers 3 of 30 images, and V3 cannot detect that

Built the site clean: **32 HTML files, 30 OG PNGs, smallest 23,221 bytes** (`og/reference/logo.png`) — the plan's "30 PNGs, smallest 23 kB" reproduces exactly. But the layout is not flat:

```
find dist/og -name '*.png'            → 30
find dist/og -maxdepth 1 -name '*.png' → 3     # index.png, install.png, quickstart.png
```

The other 27 sit under `og/guide/`, `og/modes/`, `og/performance/`, `og/reference/`, `og/rrbs/`. This follows directly from `docs/src/pages/og/[...route].ts:135` — `slug: ${toSlug(entry.id)}.png`, where `entry.id` carries the collection subdirectory.

§3 item 3 says "no zero-byte PNG **anywhere under it**", which is correct. But **A2 states the path as `docs/dist/og/*.png`**, a top-level glob, and **V3 tests `: > docs/dist/og/index.png`**, a top-level file. Both are consistent with a non-recursive implementation, and V3 would pass against one. I built the non-recursive form and ran it against a tree with a nested empty PNG:

```
# A) top-level glob — what A2's wording implies
test $(ls dist/og/*.png | wc -l) -ge 1 && test $(find dist/og/*.png -size 0 | wc -l) -eq 0
→ exit 0        # PASSES with dist/og/guide/adapters.png at 0 bytes

# B) recursive — what §3 intends
find dist/og -name '*.png' -size 0 -print   → dist/og/guide/adapters.png
→ exit 1        # correct
```

**Fix:** restate A2 as `docs/dist/og/**/*.png` (30 files, 3 at top level, 27 nested), and change V3 to truncate a **nested** PNG (`dist/og/guide/adapters.png`). As written, V3 is a test that passes while the bug is present.

### 1.4 IMPORTANT — the "no empty PNG" check is vacuous if `dist/og` is absent

The caller asked what happens if `dist/og` does not exist at all. I hit this by accident (a build failed and left no `dist`) and then confirmed it deliberately: with `dist/og` missing, **both** candidate "no zero-byte PNG" formulations exit **0**. `find` writes `No such file or directory` to stderr and exits non-zero, but its status is discarded once the output is consumed by `$(...)`/a pipeline, and `bash -e` does not abort — I verified this under `bash -e -c`.

The plan does specify an "at least one OG PNG" check alongside, so the *composite* is sound. The risk is entirely in how the single `run:` step gets written. A form that works (verified against a tree with a nested 0-byte PNG → exit 1, and against a clean tree → exit 0):

```yaml
      - name: Assert build output
        run: |
          html=$(find dist -name '*.html' | wc -l)
          png=$(find dist/og -name '*.png' | wc -l)
          echo "$html HTML pages, $png OG images"
          test "$html" -ge 1
          test "$png" -ge 1
          runt=$(find dist/og -name '*.png' -size -10240c -print)
          if [ -n "$runt" ]; then echo "Undersized OG images:"; echo "$runt"; exit 1; fi
```

Two things to keep: the count checks come **first** (they are what fails when the directory is missing), and the emptiness check reports the offending filenames rather than just a count — V3 asks for "non-zero exit naming the empty file".

### 1.5 IMPORTANT — A5 and §2.2 both state something false about `ci.yml`

Both claims are wrong on the file as it stands:

- **A5 / §10:** "Jobs in `ci.yml` are independent with no `needs:` between them." False. `ci.yml:247` (`validation`) and `ci.yml:836` (`validation-ubam`) both declare `needs: [rust-tests, lint]`.
- **§2.2:** "every job carries `if: github.event_name != 'schedule'`". False. `audit` (`ci.yml:157-163`) deliberately has **no** guard — that is the entire reason the weekly cron exists, as the comment at `ci.yml:11-15` says.

Neither breaks the plan's conclusions. A new leaf job that nothing `needs:` and that needs nothing still cannot perturb the existing graph, and the plan correctly prescribes the `schedule` guard for `docs-build` (§3 item 4, V7). But these are load-bearing sentences in the justification, and leaving them wrong misleads the next person who reads the plan as a description of the file. Reword to: "no existing job `needs:` a job I am adding, and `docs-build` needs nothing, so the graph is unchanged"; and "every job except `audit` skips on `schedule`; `docs-build` follows the majority".

### 1.6 IMPORTANT — the #247 CHANGELOG precedent is the reverse of what §4 step 5 claims

§4 step 5 and §9 item 2 assert that "prior CI-only changes such as #247 and #346 did not get entries", flagged for verification. I verified. Result: **mixed, and the closer precedent argues for an entry.**

- **#346** — `e81fc7f`, `--stat` shows `.github/workflows/ci.yml | 10 +++++++++-`, one file, no CHANGELOG. Confirms the plan for that case. It was a one-line addition to an ignore list.
- **#247** — got a dedicated section: `CHANGELOG.md:720` `#### Infrastructure (contributor-facing, since v2.1.0-beta.5) — CI hardening (#247)`, running to ~line 767. It documents, among other items, **two new CI jobs** (`#247 item 4` → the `coverage` job; `#247 item 2` → the macOS matrix). `git log -S` traces that section to `76a58de` (PR #249) — a feature PR, not the release commit, so it was written at PR time.

The pattern is "trivial CI tweak → nothing; new job or new capability → a contributor-facing line", and `docs-build` is a new job. `CHANGELOG.md:4-84` shows the current `Unreleased` block has only `#### Changes` and `#### Fixes`, so an entry means adding an `#### Infrastructure (contributor-facing)` heading. That is the maintainer's call, but the plan should not skip on a precedent that does not exist. One line is enough, e.g.:

> **New `docs-build` CI job.** Every PR now builds the Astro docs site and asserts the site produced HTML pages and non-empty satori OG images. Closes the gap where docs dependency bumps merged unvalidated.

### 1.7 Path-form trap between §3 and §4

§3 item 3 names `docs/dist` and `docs/dist/og`; §4 step 2 sets `defaults.run.working-directory: docs`. The assertion step must therefore use `dist/...`, not `docs/dist/...`. One line in §4 step 3 saying so removes an easy mis-implementation.

### 1.8 §4 step 2 checks out

`docs.yml:28-43` is exactly the four steps to copy (Checkout, Setup Node, Install dependencies, Build site), with `upload-pages-artifact` starting at `docs.yml:45`. Line range is right.

---

## 2. Assumptions

| # | Claim | Verdict | How |
|---|---|---|---|
| A1 | `npm run build` is the whole validation `docs.yml` does before deploy | **Correct** | Read `docs.yml:39-48`; the only step between build and artifact upload is the upload itself |
| A2 | satori output lands in `docs/dist/og/*.png` | **Wrong as a path form** | Built: 30 PNGs, only 3 at top level, 27 nested. See 1.3 |
| A3 | Node 24 is the correct pin | **Correct** | `docs.yml:34-35` pins `node-version: 24` with the ≥22.12.0 comment. Local Node is v24.14.0 |
| A4 | PR CI has network for `npm ci` | **Correct** | `docs.yml` already does this on push; `validation` fetches conda packages. (My own sandbox blocks `registry.npmjs.org`, which is why I used `npm ci --offline`; irrelevant to CI) |
| A5 | `ci.yml` jobs have no `needs:` | **Wrong** | `ci.yml:247`, `ci.yml:836`. Conclusion survives; see 1.5 |
| §2.3 | No `ci.yml` job uses a `paths` filter | **Correct** | `grep -n "paths" .github/workflows/ci.yml` returns five hits, all inside comment prose (`ci.yml:53,568,828,878,1017`). No trigger-level or job-level filter exists |

### Unstated assumptions worth surfacing

- **The docs build has no dependency outside `docs/`.** I checked, because `docs.yml:8` lists `CHANGELOG.md` as a trigger path. **It does not.** `docs/src/content/docs/reference/changelog.md` is a hand-synced copy — its own note says "The version below is a copy synced with the docs" — and nothing in `docs/` reads `../CHANGELOG.md` (grep over `docs/` for `CHANGELOG` hits only `docs/README.md` and that page). The copy is currently 1141 lines against the root's 1312, i.e. already lagging. Consequence for the plan: none, and `working-directory: docs` is safe. Consequence worth recording: a PR that changes only root `CHANGELOG.md` will run the new check but validate nothing new, and `docs.yml`'s `CHANGELOG.md` trigger path is effectively vestigial (it redeploys byte-identical output). Follow-up territory, not this plan's job.
- **The lockfile carries linux-specific optional deps** for `sharp` and `@resvg/resvg-js`. `docs.yml` already builds on `ubuntu-latest` on every docs push, so this is proven, but it is the thing that would break if a lockfile were ever regenerated macOS-only.
- **Cache sharing.** `setup-node`'s npm cache key is derived from the lockfile hash and the runner platform, and `ci.yml`'s `docs-build` would use the same `cache-dependency-path`, so it shares `docs.yml`'s existing cache entry rather than creating a second one. This is a real point in the plan's favour that §5 does not claim — no new pressure on the repo's 10 GB cache budget alongside the Rust caches.

---

## 3. Efficiency

Timed on this machine (Apple Silicon, Node 24.14.0), source copied to `$TMPDIR/revA-docs`:

| Step | Measured |
|---|---|
| `npm ci --offline` (warm npm cache, no `node_modules`) | **3.8 s** — 402 packages |
| `npm run build` | **13.1 s** wall — Astro reports `31 page(s) built in 10.85s`, incl. 22 image transforms, pagefind index, sitemap |

§5's "4 s / 12 s" reproduces. Two caveats on extrapolating to CI:

1. A GitHub `ubuntu-latest` runner is materially slower than this laptop, and `setup-node` itself (toolchain download + cache restore/save) adds time `npm ci` does not include. Realistic job wall-clock is **1.5-2.5 min**, not "roughly one minute". This does not change the conclusion — it runs in parallel with a 9-minute `Rust Tests (macos-latest)` — but the plan should not promise a figure it will visibly miss on the first run.
2. **Free for public repos**, so the marginal cost is zero minutes, not just "negligible". Concurrency is also fine: `ci.yml` currently starts 6 jobs at once (`rust-tests` ×2, `reproducibility`, `lint`, `audit`, `coverage`) with two more gated behind `needs`. Seven concurrent is nowhere near the runner limit.

**Unmentioned duplication.** `ci.yml` triggers on `push` to `master`/`dev` as well as `pull_request` (`ci.yml:4-10`), and §3 item 1 embraces that. So every docs-touching push to `dev` will build the site **twice** — once in `ci.yml`'s `docs-build`, once in `docs.yml`'s `build` before deploy. ~13 s of duplicated runner time, on a free tier. I would leave it (an `if: github.event_name == 'pull_request'` guard buys nothing and adds a condition), but the plan should name it so it reads as a decision rather than an oversight.

---

## 4. Validation sufficiency

The validation table is unusually good for a CI change — V2 and V3 are exactly the right instinct, and the plan says so. Gaps:

1. **V3 cannot fail on the bug in 1.3.** Truncating `dist/og/index.png` exercises only the top-level path. Retarget it at `dist/og/guide/adapters.png`.
2. **No case covering "`dist/og` absent entirely".** This is the highest-value missing row, because I verified the naive emptiness check passes silently in that state (1.4). Add: *V3b — `rm -rf dist/og`, run the assertion step's command, expect non-zero.*
3. **V6's command does not run on this machine.** `python3 -c "import yaml"` → `ModuleNotFoundError: No module named 'yaml'` (verified; there is no `yamllint` or `actionlint` on `$PATH` either). A verified working substitute:
   ```
   ruby -ryaml -e 'p YAML.load_file(".github/workflows/ci.yml")["jobs"].keys'
   ```
   which prints `["rust-tests", "reproducibility", "lint", "audit", "coverage", "validation", "validation-ubam"]` today — so the "job count +1" check is directly usable, expecting `docs-build` to appear between `lint` and `audit`.
4. **V2's chosen breakage should be an OG-route breakage, not frontmatter.** A malformed frontmatter key proves Astro fails the build, which was never in doubt. The interesting question is whether the *satori path* can fail the job — and per 1.2, a throwing endpoint does (exit 1) while an empty write does not. If V2 must pick one, breaking the OG route exercises the path this job exists for.
5. **V5's baseline is right:** 8 checks today (`rust-tests` ×2, `reproducibility`, `lint`, `audit`, `coverage`, `validation`, `validation-ubam`), so 9 after.

### Where the job could pass while protecting nothing

Ranked by likelihood, all specific to this design:

1. Non-recursive glob → 27 of 30 images unchecked. **Verified reachable** (1.3).
2. `dist/og` missing and only the emptiness check runs → silent pass. **Verified reachable** (1.4).
3. A **truncated but non-empty** PNG (say 500 bytes of a 23 kB image). `-size 0` and `-empty` both miss it. See 5.1 for a zero-maintenance fix.
4. A valid PNG that renders **blank or garbled** — outside any cheap check. §10 acknowledges this and the trade-off is right.

---

## 5. Alternatives

### 5.1 Strengthen the size predicate (cheap, recommended)

`-size 0` / `-empty` only catch truly empty files. A byte-exact floor catches truncation too, at zero maintenance cost, because the floor is nowhere near any real value (smallest real OG image: 23,221 bytes).

Careful with units — I tested this, and the obvious form is wrong:

```
# tree contains a 0-byte PNG and a deliberately truncated 500-byte PNG
find dist/og -name '*.png' -size -1k      → adapters.png            # MISSES the 500-byte file
find dist/og -name '*.png' -size -10240c  → index.png, adapters.png # catches both
```

Both GNU and BSD `find` round `-size -1k` up to 512/1024-byte blocks, so a 500-byte file counts as one block and escapes. Use `c` (bytes): `-size -10240c`. A 10 kB floor is a 2× margin under the smallest real image and never needs touching as pages are added — strictly better than "non-empty" for the same one line.

### 5.2 `pull_request` on `docs.yml` — the plan's rejection is right, but for a reason it did not give

Assessing §2.2's three grounds on the merits:

- **Permissions — real but soluble.** For same-repo `pull_request` runs (which is how this repo's branches and human PRs work), the `GITHUB_TOKEN` *does* receive the workflow-level `pages: write` / `id-token: write` from `docs.yml:12-15`; the read-only downgrade applies to forks and Dependabot. So the concern is genuine. But **job-level `permissions:` fully overrides workflow-level**, so `permissions: {contents: read}` on `build` closes it completely. The plan dismisses this variant a little too fast.
- **Concurrency — real but soluble.** `docs.yml:17-19` (`group: pages`, `cancel-in-progress: false`) would indeed serialise PR builds behind live deploys. But `concurrency.group` accepts expressions, so a per-event group would fix it.
- **Deploy event guard — real, and the plan's framing is fair.** One more condition on the deploy path is a genuine correctness surface, and the deploy path is the one thing worth not touching.

**The argument §2.2 missed, which is decisive:** `docs.yml:6-9` carries a `paths` filter. On a `pull_request` that touches no `docs/**` path, a path-filtered workflow reports **nothing at all** — not a skipped check, not a neutral one. If that check were ever made required, such PRs would block forever, and the standard workaround is a duplicate no-op workflow with the same name. That is strictly worse than §2.3's "skipped reads ambiguously" complaint and it kills the option outright. Worth adding: it converts §2.2 from a judgement call into a clear-cut one.

**Also worth doing given this framing:** `ci.yml` has **no** top-level `permissions:` block (`grep` finds `permissions` only at `ci.yml:160`, job-level on `audit`), so `docs-build` would inherit the repo-default token scope. Adding `permissions: {contents: read}` to the new job makes §2.2's least-privilege argument actually true of the result rather than incidentally true. One line.

### 5.3 Shared composite action (§9 item 3) — agree with the plan

Four steps duplicated across two workflows, with Dependabot's `github-actions` group updating both pins together (verified: `.github/workflows/*.yml` all use the same pinned SHAs, e.g. `actions/checkout@3d3c42e…` at both `ci.yml:32` and `docs.yml:29`). A composite action would add a file and an indirection for less than it saves. The cross-reference comments are the right weight, and one line each complies with `CLAUDE.md`'s comment rule.

One factual note on that: editing `docs.yml` to add its cross-reference comment matches `docs.yml:9`'s own trigger path, so merging to `dev` will fire a real deploy of byte-identical content. Harmless; just not a surprise.

### 5.4 Additional silent-failure surfaces (optional, same cost as the PNG check)

The clean build log shows two more artefacts nothing would notice losing:

```
[starlight:pagefind] Found 32 HTML files. → dist/pagefind/pagefind-entry.json
[@astrojs/sitemap] `sitemap-index.xml` created at `dist` → dist/sitemap-index.xml
```

If Starlight's pagefind step ever degrades across a version bump, site search dies with a green build — the same failure shape as the OG images, and the same one-line check (`test -s dist/pagefind/pagefind-entry.json`). Not required; the plan's scope is defensible as drawn. Mentioning because the incremental cost is a single `test`.

Broken internal links are the other classic Astro silent failure, but that needs a plugin (`starlight-links-validator`) and a new dependency — correctly out of scope here.

### 5.5 Adjacent pre-existing bug (out of scope, flagging only)

`justfile:89-94` runs `cd Docs && npm run build` / `cd Docs && npm run dev`, capital `D`. The tracked directory is `docs/` — `git ls-files | grep -c '^Docs/'` returns **0**. It works on case-insensitive macOS and fails on Linux. Since the `justfile` exists specifically for local CI parity (documented at `CHANGELOG.md:742-747` as "#247 item 7 — `justfile` for local CI parity"), and this plan is about making the docs build a CI gate, a one-character fix in the same PR would be coherent. Entirely the maintainer's call; not a defect in the plan.

---

## 6. Action items

### Critical

None. Nothing in the plan can break the build, the deploy path, or an existing check.

### Important

1. **Fix A2's path form and V3's target** (1.3). A2 → `docs/dist/og/**/*.png`, 30 files of which 27 are nested; V3 → truncate `dist/og/guide/adapters.png`. Verified: the top-level form passes with a nested 0-byte PNG present.
2. **Specify the assertion so the missing-directory case fails** (1.4). Count checks first, emptiness check second, offending filenames echoed. Snippet in 1.4. Verified: the emptiness check alone exits 0 when `dist/og` is absent, even under `bash -e`.
3. **Add V3b: `rm -rf dist/og`, run the assertion, expect non-zero.**
4. **Correct A5 and §2.2's schedule-guard sentence** (1.5). `ci.yml:247` and `ci.yml:836` both have `needs:`; `audit` at `ci.yml:157` has no schedule guard. Conclusions unchanged, premises wrong.
5. **Revisit §4 step 5 / §9 item 2 on the CHANGELOG** (1.6). #247 *did* get a contributor-facing section at `CHANGELOG.md:720`, added in PR #249, covering two new CI jobs. Either add a one-line `#### Infrastructure (contributor-facing)` entry under `Unreleased`, or state the actual rule; do not cite #247 as precedent for skipping.
6. **Note in §4 step 3 that the assertion paths are relative to `docs/`** (1.7), since `working-directory: docs` is in force and §3 writes them as `docs/dist/...`.

### Optional

7. **Use `-size -10240c` instead of `-size 0`** (5.1) — catches truncation as well as emptiness, never needs updating. Verified that `-size -1k` misses a 500-byte file.
8. **Add `permissions: {contents: read}` to the `docs-build` job** (5.2) — `ci.yml` has no top-level `permissions:` block, so the new job otherwise inherits the repo default.
9. **Add the path-filter argument to §2.2** (5.2) — a path-filtered `pull_request` workflow reports *nothing* on non-docs PRs, which is worse than skipped and decisive against the `docs.yml` option.
10. **Adjust §5's cost figure** to 1.5-2.5 min on a hosted runner, and note the build is free on a public repo and shares `docs.yml`'s existing npm cache entry (§3 above).
11. **Name the push-path duplication** (§3) — docs-touching pushes to `dev` build the site twice, ~13 s, accepted deliberately.
12. **Replace V6's command** with the verified `ruby -ryaml` one-liner; `python3` has no `yaml` module on this machine.
13. **Point V2 at the OG route** rather than frontmatter (§4.4), so it exercises the path the job exists for.
14. **Close §6/§9 item 1:** verified via `gh api` that **neither `dev` nor `master` is protected** (both return `Branch not protected`, HTTP 404) and `repos/…/rulesets` returns `[]`. There are no required checks to add; the open question can be marked resolved. Side effect: §2.3's "complicates branch protection" support for always-run is currently vacuous — the remaining "a skipped check reads ambiguously" is thinner, though always-run is still the shape that survives protection being introduced later, so the decision stands.
15. **Optionally assert `dist/pagefind/pagefind-entry.json` and `dist/sitemap-index.xml`** (5.4) — same failure shape as the OG images, one `test` each.
16. **Pre-existing, out of scope:** `justfile:89,93` say `cd Docs`; the directory is `docs/`. Broken on Linux (5.5).

---

## 7. What I confirmed vs inferred

**Confirmed by running commands:** the trigger gap and that `ci.yml` is the only PR-triggered workflow; `release.yml` never builds docs; no `paths` filter in `ci.yml`; `needs:` at `ci.yml:247`/`ci.yml:836`; `audit` has no schedule guard; 32 HTML / 30 OG PNGs with 27 nested; smallest PNG 23,221 bytes; empty-buffer render → 0-byte PNG with exit 0; throwing endpoint → exit 1; top-level glob passes with a nested empty PNG; emptiness check passes with `dist/og` absent; `-size -1k` misses a 500-byte file while `-size -10240c` catches it; `npm ci --offline` 3.8 s, `npm run build` 13.1 s; #346 touched only `ci.yml`; #247's CHANGELOG section at line 720 traced to `76a58de`; `dev` and `master` unprotected with no rulesets; `python3` lacks `yaml`; `git ls-files` has no `Docs/` entries.

**Inferred, not confirmed:** hosted-runner timing (extrapolated from local timings plus typical `setup-node` overhead); GitHub's same-repo-versus-fork `GITHUB_TOKEN` behaviour under `pull_request` (documented behaviour, not something I could exercise here); that a path-filtered `pull_request` workflow reports no check rather than a skipped one (documented behaviour); that the ubuntu build produces the same 30 PNGs as macOS (strongly implied by `docs.yml` deploying successfully today).
