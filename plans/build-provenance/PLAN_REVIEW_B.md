# Plan Review B — build-provenance

Reviewer B, independent context. Target: `plans/build-provenance/PLAN.md`.

---

## Logic

- **Implementation order risk (Section 7)**: order is build.rs → clippy fixes → workflow. Safe for CI only because the `lint` job does not exist until step 3. BUT: reviewers running `cargo clippy -D warnings` locally on any mid-series commit will see red until step 2 lands. Not a blocker, just document that intermediate commits on the feature branch are not lint-clean.
- **`validation` gating on `lint` only, not `audit`**: the plan is right that `audit` is a signal not a gate (new CVEs should not block parity work). Confirmed.
- **`release.yml` not touched**: plan says nothing about it. `release.yml` has its own job DAG (`check-release → build-binaries → …`, lines 116/179/200/257/313/331/361/387). Should release builds also fail-fast on lint drift? Release builds ship to crates.io/GHCR — a lint-dirty release would be unfortunate. Either (a) add `lint` gate in `release.yml` too, or (b) document the intentional split. Plan is silent. **Important.**
- **Banner/stderr interaction**: the plan notes risk #10 (extra stderr line is additive) but CI step "Validate multi-pair paired-end" greps `grep -c 'Auto-detecting adapter type'` — the new banner line doesn't collide. Safe, but a scan of the full ci.yml for any startup-banner regex would be prudent. I didn't find any.
- **Step 1 banner replaces lines 25–26**: verified against actual file — lines 25–26 match exactly. Safe.
- **Clippy table line numbers**: I ran `cargo clippy --all-targets --release` at HEAD `aa92493`. All 13 warnings listed are present, exact line numbers match (`src/parallel.rs:58/223/300`, `src/trimmer.rs:69/163/332/639`, `src/main.rs:159/310`, `src/demux.rs:232`, `src/report.rs:735`, `src/alignment.rs:58`, `src/quality.rs:135`). Warning count total: 13 distinct, `--all-targets` surfaces no additional test-only warnings. Plan is accurate.

## Assumptions

- **`git rev-parse --short` length**: Git 2.11+ auto-scales `core.abbrev` based on repo object count. Local git reports 7-char hash at HEAD; a CI runner with a different Git version or a shallow clone may produce a different length (uncommon but possible). The plan relies on visual stability ("`abc1234 — …`" in the example output), not programmatic parsing, so not load-bearing. **Optional**: pin with `git rev-parse HEAD | cut -c1-7` for stability.
- **`SOURCE_DATE_EPOCH` malformed handling**: risk #6 covers overflow/parse-fail with `unwrap_or_else(fallback_to_now())`. This silently falls back to wall-clock on malformed input, which breaks reproducibility silently. Debian reproducible-builds guidance recommends **failing the build** on malformed `SOURCE_DATE_EPOCH` (see reproducible-builds.org/docs/source-date-epoch/). **Important**: make the fallback an error, or at least `println!("cargo:warning=SOURCE_DATE_EPOCH malformed, falling back to now")`.
- **`cargo:rerun-if-changed=.git/index`**: this triggers on any `git add`/`git rm`/staging change, not just commits. In an iterative dev loop where a dev stages unrelated files, `build.rs` reruns and `BUILD_TIMESTAMP` changes → the bin relinks. Low-cost (~100ms) but surprising. `.git/HEAD` alone covers commit + branch switch; consider dropping `.git/index`. **Optional**.
- **clap `long_version` with derive + `concat!`**: clap v4 derive accepts `long_version = <str>` where the RHS is evaluated as a `&'static str` expression. `concat!(env!("A"), "\n", env!("B"))` resolves at compile time to a `&'static str` — compiles fine. Clap-internally, `long_version` is shown ONLY for `--version`; `-V`/`--help`'s short line uses `version`. Verified per clap 4 docs. The plan's claim holds.
- **Dependabot `dependencies` label**: Dependabot DOES auto-create the `dependencies` label on first PR (documented behaviour per GitHub docs at docs.github.com/en/code-security/dependabot/working-with-dependabot/dependabot-options-reference). Plan is correct; risk #8 (preflight note) is well-placed.
- **`build.rs` autodetect at package root**: `cargo` has autodetected `build.rs` at crate root since pre-1.0. No `build = "build.rs"` key needed. `exclude` list in `Cargo.toml` does NOT list `build.rs`, so the published tarball will include it. **Verified correct.**
- **Implicit assumption**: the plan silently assumes Rust edition 2024 supports `env!` in derive-macro `concat!` — this is unchanged across editions, safe.

## Efficiency

- **`lint` cold cache cost**: plan estimates 4-5 min. Realistic on a fresh runner; `--release` forces full optimized compile of all deps. Running parallel to `rust-tests` means no wall-clock regression in the normal case. Acceptable.
- **Toolchain install duplication**: both `lint` and `rust-tests` run `dtolnay/rust-toolchain@stable` from scratch. ~20-30s overhead. Merging into a single job would save the download but costs separation-of-concerns. Plan's choice is defensible.
- **`rustsec/audit-check@v2` caching**: per the action's source (rustsec/audit-check on GitHub), it runs `cargo install cargo-audit` if not cached, then fetches the advisory DB each run. No built-in DB caching. First run ~90s, subsequent with toolchain cache ~30s. Acceptable but not free. **Optional**: explicit `actions/cache@v4` around `~/.cargo/bin/cargo-audit` and `~/.cargo/advisory-db`.
- **`--release` in clippy**: doubles `target/` cache size. Plan notes this. An alternative: run clippy in dev profile (fastest), and rely on `rust-tests` for release warnings. Debatable — shipped warnings are the ones that matter for correctness.

## Validation sufficiency

- **`wc -l | grep -q 2` assertion**: fragile. `--version` output ends with a trailing `\n` → `wc -l` reports 2 for a 2-line string; OK. But if clap ever emits a trailing blank line (seen in some versions), it becomes 3 and the grep fails. **Important**: use a targeted grep: `./target/release/trim_galore --version | grep -qE '^[0-9a-f]{7,40} — '` to assert presence of provenance line, not line count.
- **Reproducibility check via `strings | grep`**: `strings` will find the timestamp literal unless the linker has deduped it with another identical string (vanishingly unlikely for a `YYYY-MM-DDTHH:MM:SSZ` pattern). Fine, but direct: compare `sha256sum` of two successive `SOURCE_DATE_EPOCH=X` builds and assert equal. Stronger. **Important**.
- **`.git` missing simulation**: `mv .git /tmp/savegit && cargo clean && cargo build --release`. `cargo clean` wipes `target/`, which forces `build.rs` re-execution. Sufficient. Restoring `.git` after the test is NOT in the plan — add `mv /tmp/savegit .git` to the recipe or the dev corrupts their working tree.
- **No test for `SOURCE_DATE_EPOCH` malformed**: plan describes the fallback in risk #6 but no verification step exercises it. **Optional** CI step: `SOURCE_DATE_EPOCH=garbage cargo build --release` should not panic.
- **No test for fork-PR `audit-check` permission case**: plan flags the risk (#9) but no way to exercise it pre-merge. Acceptable; hard to test ahead of time.

## Alternatives

- **`cargo-deny` vs `rustsec/audit-check@v2`**: cargo-deny is broader (licence checks, banned crates, duplicate-version checks) and widely used in the Rust ecosystem (tokio, serde, bevy). `audit-check` is narrower (CVE-only). For v2.1.0 GA hygiene, `audit-check` is right-sized — adding cargo-deny is a policy decision that should be a separate PR. Keep `audit-check`.
- **Single combined lint job vs split**: three jobs (`lint`, `audit`, `rust-tests`) keep failure-surface clean and parallel-friendly, at the cost of ~1 min total CI compute across runners. A combined job would halve wall-clock for small PRs. Plan's choice is reasonable and matches ruSTAR. Accept.
- **`cargo install cargo-audit && cargo audit` vs `rustsec/audit-check@v2`**: the dedicated action gives PR annotations and GitHub check integration — worth the opacity. Keep plan's choice.
- **`build.rs` also emits `CARGO_CFG_TARGET_ENV`**: would distinguish `-gnu` vs `-musl` in `--version` output (visible to users debugging glibc/musl-specific issues). Cheap addition (one `env::var`, one `rustc-env=`). Version line could become `<hash> — <os>/<arch>/<env> — built <ts>`. **Optional** but recommended — matches what `rustc -Vv` shows.

## Action items

### Critical
_(none — plan is broadly sound)_

### Important
1. **`SOURCE_DATE_EPOCH` malformed handling** (Section 4 Step 1, risk #6): fail-build or emit `cargo:warning=` — do not silently fall back to `SystemTime::now()`. Breaks reproducibility contract.
2. **CI `--version` assertion** (Section 5 Step 1): replace `wc -l | grep -q 2` with `grep -qE '^[0-9a-f]{7,40} — '` to check provenance-line presence, not line count.
3. **Reproducibility verification** (Section 5 Step 1): compare `sha256sum` of two `SOURCE_DATE_EPOCH=X` builds instead of `strings | grep`.
4. **`release.yml` lint gating decision** (Section 4 Step 2): either extend the `lint` gate to `release.yml`'s `build-binaries` or document why releases don't need the same pre-flight.
5. **Clean up `.git` restore in verification recipe** (Section 5 Step 1, edge cases): add `mv /tmp/savegit .git` after the "`.git` missing" simulation step.

### Optional
6. **Drop `cargo:rerun-if-changed=.git/index`** (Section 4 Step 1): `.git/HEAD` alone is enough for commits + branch-switch; `.git/index` causes rebuilds on any staging change.
7. **Pin hash length** (Section 4 Step 1): `git rev-parse HEAD | cut -c1-7` for cross-git-version stability.
8. **Emit `CARGO_CFG_TARGET_ENV`** (Section 4 Step 1): enrich version line with `-gnu` vs `-musl` distinction.
9. **Cache `cargo-audit` / advisory-db** (Section 4 Step 4): `actions/cache@v4` around `~/.cargo/bin/cargo-audit` + advisory DB, shaves ~60s off cold audit runs.
10. **Add CI smoke for malformed `SOURCE_DATE_EPOCH`** (Section 5 Step 1): `SOURCE_DATE_EPOCH=garbage cargo build --release` should succeed (or fail cleanly depending on Important #1).

---

**Verified against repo at commit `aa92493`**: clippy table accurate, file locations present, test fixture valid, Dockerfile uses `rust:1.85-bookworm` (git available in Debian bookworm rust image), `validation` job dependency confirmed at `.github/workflows/ci.yml:40`.
