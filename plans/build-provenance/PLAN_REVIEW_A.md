# Plan Review A — build-provenance

Reviewer A, fresh context, no coordination with Reviewer B.

Plan file: `/Users/fkrueger/Github/TrimGalore/plans/build-provenance/PLAN.md`

## Repo-state verification (claims vs reality)

- **`Cargo.toml` auto-detects `build.rs`**: Confirmed. Cargo has auto-detected `build.rs` at package root since 2015; edition `2024` does not change this. The `exclude` list (`test_files/`, `plans/`, `docs/`, `optimus_prime/`, `.github/`, `.claude/`, `CLAUDE.md`, `CHANGELOG.md`) does **not** exclude `build.rs`. Plan section 3 is correct.
- **`src/cli.rs` current attr** (line 11): `#[clap(name = "trim_galore", version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)"), about)]` — matches the plan's pre-image.
- **`src/main.rs:25–26`**: Matches verbatim. Only one `FastqReader::sanity_check` consumes the banner (stderr, no downstream parser). Adding one more `eprintln!` is safe.
- **`.github/workflows/ci.yml`**: Two jobs (`rust-tests`, `validation: needs: [rust-tests]`). Structure matches plan.
- **Clippy warning table**: Ran `cargo clippy --all-targets --release` — all 13 plan entries match actual output (file, line, lint name). **No missing warnings; no extra warnings.** Clippy counts in output: lib = 10, bin = 2, tests add 1 unique (`len_zero` at `trimmer.rs:639`). Total 13 — matches plan.
- **`test_files/illumina_10K.fastq.gz`**: Exists. OK.
- **`Dockerfile`**: Confirmed `FROM rust:1.85-bookworm AS builder` — git is present in bookworm's default apt layer of the rust image. Plan risk #1 is correct.
- **`build.rs` and `.github/dependabot.yml`**: Do not yet exist — correctly treated as new files.

## Logic review

- **Banner insertion preserves ordering**: no downstream code grepping the stderr header found in the repo. Additive change is safe.
- **Job graph change** (`validation: needs: [rust-tests, lint]`): correct — `lint` runs in parallel with `rust-tests`, joins before `validation`. `audit` intentionally not a gate. Consistent.
- **`needless_range_loop` parity risk** is well-contained: parity via `validation` job. Good.
- **Gap — `cargo fmt --check` in `lint` job but no `rustfmt.toml` in repo**: if contributors have different default rustfmt configs this could thrash. Low risk but worth noting (`rustfmt` from `dtolnay/rust-toolchain@stable` pins to toolchain default — acceptable).
- **Gap — `lint` job has no `Cargo.lock` dependency path** for `restore-keys`; plan uses `key: ${{ runner.os }}-cargo-lint-${{ hashFiles('Cargo.lock') }}`. OK but `restore-keys: ${{ runner.os }}-cargo-` will pull the `rust-tests` cache on first run (different target artifacts) — initial lint run will be slow regardless, then diverge. Matches plan's own caveat at risk #2.
- **Gap — `--version` CI assertion uses `wc -l`**: `clap --version` writes to **stdout**, so `./target/release/trim_galore --version | wc -l` captures it. However, `clap` currently does **not** append a trailing newline after the last line (verify locally before pushing). `wc -l` counts newlines, so a two-line string `"A\nB\n"` → 2, but `"A\nB"` → 1. **This assertion could false-fail.** Use `grep -c '' file` against captured output, or use `./target/release/trim_galore --version | awk 'END{print NR}'` to get *line count* regardless of trailing newline.
- **Gap — reproducibility check**: `strings target/release/trim_galore | grep -o '2023-11-14T..:..:..Z' | head -1`. With `strip = "debuginfo"` in release profile, the `VERSION_BODY` string itself (used by `long_version`) is retained in `.rodata` — safe. But if Rust const-folds the `concat!()` into a single merged string, `grep -o` for just the timestamp substring still matches. OK.
- **Gap — `SOURCE_DATE_EPOCH` spec says decimal integer** (https://reproducible-builds.org/specs/source-date-epoch/). Signed vs unsigned: the spec explicitly says "non-negative decimal number", so `u64::from_str` is correct. Plan is right.
- **Gap — detached HEAD / worktrees / packed-refs**: `cargo:rerun-if-changed=.git/HEAD` + `.git/index` covers most cases, but:
  - On a worktree, `.git` is a file pointing to `worktrees/<name>`; `.git/HEAD` does not exist as a regular file. Build script will see a stale `GIT_SHORT_HASH` if a new commit arrives. Low impact (devs rarely use worktrees for this repo) but worth a sentence in the plan.
  - On detached HEAD with `git checkout <sha>`, `.git/HEAD` changes — covered.
  - `packed-refs`: only matters if a local branch advances via `git pack-refs` without index churn. Extremely rare in normal dev.
  - Mitigation: add `cargo:rerun-if-changed=.git/refs` as a belt-and-braces trigger. Optional.

## Assumptions (stated + implicit)

Stated and verified:
- `git` available in CI and Docker build stage (confirmed).
- `strip = "debuginfo"` retains string constants (correct — it strips debug info, not `.rodata`).
- Cargo auto-detects `build.rs` (correct for edition 2024).

Implicit / under-examined:
- **clap derive — `long_version`**: Assumption that `long_version` does NOT override `-V`. Verify against clap v4 derive docs: `version` maps to `Command::version()`, `long_version` to `Command::long_version()`. `-V` prints `version`, `--version` prints `long_version` if set, else falls back. Confirmed correct behaviour in clap v4.x. Plan is right.
- **`CARGO_CFG_TARGET_OS` / `CARGO_CFG_TARGET_ARCH` available in build script**: yes, cargo guarantees these.
- **`wc -l` behaviour** — see gap above. Implicit and load-bearing.
- **Dependabot auto-creates `dependencies` label** — correct as of GitHub's current behaviour (Dec 2024+), but fragile. Plan flags this.
- **`rustsec/audit-check@v2` network access in CI**: needs outbound HTTPS to advisory-db mirror. Ubuntu runners allow it; not mentioned in plan.

## Efficiency analysis

- **Lint cache footprint claim (2× cache)**: Plan's "doubles cache footprint" is accurate for `target/` because dev and release profiles build separate artifact trees (`target/debug/`, `target/release/`). `~/.cargo/registry` is shared — that dominates cache size (~400MB) and is not duplicated. Net footprint impact: +~150MB (release `target/` build). Plan's "2×" is a small overestimate but not misleading.
- **`--all-targets --release` duplicate work vs `rust-tests`**: `rust-tests` currently does `cargo test` (debug) + `cargo build --release` (no tests). `lint` job does `cargo clippy --all-targets --release` — compiles release + test binaries together. Distinct compilation units; no duplicate work on the same cache key since the lint job has its own cache key. If we *merged* lint into `rust-tests`, we'd save the release-test-compile cost (est. ~30–60 s on warm cache). See "Alternatives" for the trade-off.
- **`audit-check` runtime**: ~10–15 s on warm cache. Parallel to lint/tests. Negligible.

## Validation sufficiency

Gaps:
1. **`wc -l` CI assertion** (see Logic). Actionable.
2. **No test that `-V` output length < `--version` output length** — a regression where someone accidentally unifies the two wouldn't trip `wc -l = 2` if both become 2-line. Suggest: `[ "$(./target/release/trim_galore -V | awk 'END{print NR}')" -lt "$(./target/release/trim_galore --version | awk 'END{print NR}')" ]`.
3. **No reproducibility check in CI** — the plan only covers local verification. A one-shot CI matrix job (build twice with same `SOURCE_DATE_EPOCH`, diff binaries) would catch regressions. Optional but recommended for the "embedded-in-nf-core" use case cited in Section 1.
4. **No check that `GIT_SHORT_HASH=unknown` fallback actually works in `cargo install` flow** — risk #1 claims it; no verification step. Add: `cargo install --path . --locked` in a dir with `.git` temporarily hidden, assert `--version` succeeds and contains "unknown".
5. **Parity validation for `needless_range_loop` refactors**: the existing `validation` job md5s single-end + paired-end + clock + demux + hardtrim5. That covers the hot paths these refactors touch. Sufficient.

## Alternatives worth considering

- **`vergen` vs hand-rolled**: `vergen` ~20 transitive deps just for git/build info. Plan's rejection is correct; hand-rolled `Command::new("git")` is ~20 lines and no deps.
- **`cargo-deny` vs `rustsec/audit-check@v2`**: `cargo-deny` adds license-allowlist + duplicate-crate detection on top of advisory checks. More value per CI-minute, but heavier config surface. For v2.1.0 GA, `audit-check` is the right minimal choice. Revisit post-GA if Altos needs license governance.
- **Clippy in same job as `cargo test` vs separate**: merging saves ~30–60 s of toolchain download + one cache restore. Loses: (a) separate failure signal, (b) parallelism with `rust-tests`. Plan's "separate job" choice is right for a public-facing repo; the parallelism cancels the toolchain-download cost.
- **`chrono` for ISO-8601 formatting**: plan avoids it correctly. `format!("{}-{:02}-{:02}T{:02}:{:02}:{:02}Z", ...)` with manual epoch-to-datetime math (leap-year aware) is ~30 lines of stdlib code. Acceptable. Alternative: `humantime` crate (one small dep). Plan's choice is fine.

## Action items

### Critical
1. **Fix the `wc -l` CI assertion** (Plan §5 Step 1, CI assertion block around line 348–350). Replace with `awk 'END{print NR}'` or similar to be robust to trailing-newline presence. Without this fix, the provenance assertion can false-fail on any platform where clap omits the trailing newline.

### Important
2. **Add `-V` length < `--version` length check** to CI (Plan §5 Step 1). Protects against the two forms collapsing to the same output.
3. **Document worktree / detached-HEAD behaviour** in Plan §4 Step 1 build.rs design. Add `cargo:rerun-if-changed=.git/refs` as a belt-and-braces trigger, or document the known limitation.
4. **Add `cargo install --path .` + stripped-`.git` smoke test** to §5 Step 1 verification. Plan claims the fallback works (§6 risk 1) but never exercises it.
5. **Reproducibility CI assertion** (Plan §5 Step 1 Reproducibility block is local-only). Add a CI matrix job that builds twice with the same `SOURCE_DATE_EPOCH` and diffs the binaries. Strongest signal for the nf-core embedded-pipeline use case.

### Optional
6. Add a note in Plan §4 Step 2 that `rustfmt` uses toolchain-default config (no `rustfmt.toml`), so any future style preferences need a committed config file.
7. Consider `cargo-deny` as a post-GA follow-up (license allowlist + duplicate detection). Out of scope for this plan; worth a line in §2 Non-goals.
8. Plan §4 Step 3 (Dependabot): default to creating the `dependencies` label explicitly via `gh label create` before merge, instead of relying on Dependabot's auto-create. Removes a race condition with the first PR.
