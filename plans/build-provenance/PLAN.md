# Plan: Build Provenance, CI Lint Gates, Dependabot, and Cargo Audit

Feature branch: `build-provenance` (off `optimus_prime`)
Target merge: `optimus_prime` (v2.1.0 beta window; GA ≥ 2026-05-03)
Inspiration: https://github.com/Psy-Fer/ruSTAR/pull/2 (Ewels, merged 2026-04-17)

**Revision history**
- 2026-04-19: Dual plan-review pass (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`).
  Folded in 5 edits: (1) replace fragile `wc -l` CI assertion with a
  content-targeted grep; (2) strengthen reproducibility check to sha256sum
  diff; (3) add `.git`-restore + `cargo install --path .` fallback smoke in
  the verification recipe; (4) change `SOURCE_DATE_EPOCH` malformed-input
  handling from silent `now()` fallback to hard-fail; (5) document that
  `release.yml` intentionally does NOT gate on `lint`. Optional items
  (hash-length pinning, `CARGO_CFG_TARGET_ENV`, audit-db caching,
  `-V`/`--version` length check, worktree behaviour) deferred — see
  `PLAN_REVIEW_A.md` / `PLAN_REVIEW_B.md` for details.

---

## 1. Context

TrimGalore Optimus Prime is a single-binary Rust rewrite (Cargo crate at repo
root, `src/main.rs`). We are in the v2.1.0 beta window and need pre-GA hygiene
that (a) lets users embedding the binary in published pipelines (nf-core/rnaseq)
pin builds reproducibly, and (b) gives CI the means to catch the class of
regressions we've been chasing during beta.

Four related improvements, all lifted from ruSTAR's PR #2 and scoped down for
our repo:

1. `build.rs` that embeds git hash + build timestamp + target triple into a
   rich `--version` and startup banner.
2. CI lint gate (`cargo fmt --check`, `cargo clippy -D warnings`) that gates
   the expensive `validation` job.
3. Dependabot for `cargo` and `github-actions` ecosystems.
4. `rustsec/audit-check` for CVE advisories in the dependency tree.

---

## 2. Goals / Non-goals

### Goals

- `trim_galore --version` prints a 2-line form: version + `<git-hash> — <os>/<arch> — built <iso-8601-utc>`.
- `trim_galore -V` stays as a single line (PEP-8-esque short form).
- Startup banner gains the provenance line after the existing `"Trim Galore - Oxidized Edition v…"` line.
- Reproducible builds: honour `SOURCE_DATE_EPOCH` (Debian/Conda convention).
- CI fails on fmt drift or new clippy warnings (`-D warnings`, `--all-targets`, `--release`).
- CI flags new dependency CVEs via `rustsec/audit-check@v2`.
- Dependabot opens weekly PRs for cargo deps and GH action pins.

### Non-goals (explicit)

- **SIMD binary variants** (x86-64-v3/v4, apple-m1, neoverse-v1 SVE) — deferred
  post-GA. TrimGalore is I/O-bound at >16 cores; the 9× build matrix is not
  justified for a single-digit-percent speedup.
- **CPU runtime detection (`src/cpu.rs`)** — coupled to SIMD variants; deferred.
- **MSRV bump to 1.88** — we don't use let-chains; stay at 1.85 (current
  `rust-version` in `Cargo.toml`).
- **`cargo install` name gating** — we already own `trim-galore` on crates.io.
- **Moving release profile `lto = true` → `lto = "fat"`** — semantically identical.

---

## 3. Critical files

- `/Users/fkrueger/Github/TrimGalore/Cargo.toml` — may need `build = "build.rs"` (likely auto-detected; document which)
- `/Users/fkrueger/Github/TrimGalore/build.rs` — new file (repo root)
- `/Users/fkrueger/Github/TrimGalore/src/cli.rs` — `long_version` + `version` on `#[clap]`
- `/Users/fkrueger/Github/TrimGalore/src/main.rs` — banner update + each `type_complexity` warning
- `/Users/fkrueger/Github/TrimGalore/.github/workflows/ci.yml` — lint + audit jobs, job dependencies
- `/Users/fkrueger/Github/TrimGalore/.github/dependabot.yml` — new file

Supporting touch-ups (concrete clippy resolutions, all small):

- `src/parallel.rs` — three `too_many_arguments` allows
- `src/trimmer.rs` — one `too_many_arguments` allow + `collapsible_if` fix + `doc_lazy_continuation` fix + test `len_zero` fix
- `src/alignment.rs` — `needless_range_loop` fix (enumerate)
- `src/quality.rs` — `needless_range_loop` fix (enumerate)
- `src/report.rs` — `int_plus_one` fix (`range_start < len`)
- `src/demux.rs` — `write!` → `writeln!` (strip trailing `\n`)

---

## 4. Design

### Step 1 — `build.rs` + rich `--version` / banner

**New file** `/Users/fkrueger/Github/TrimGalore/build.rs`:

- Shell out to `git rev-parse --short HEAD` via `std::process::Command`. No
  `[build-dependencies]` needed. (Decision: avoid pulling in `vergen` or
  `git2` — a 6-line Command call is cheaper than a crate.)
- Emit three `cargo:rustc-env=` vars:
  * `GIT_SHORT_HASH` — git output, stripped; fallback to `"unknown"` when
    `.git` missing (cargo-install from crates.io, tarball builds), git not on
    PATH, or command exits non-zero.
  * `BUILD_TIMESTAMP` — ISO 8601 UTC. If `SOURCE_DATE_EPOCH` is set,
    parse with `u64::from_str(env)`; on parse failure, **hard-fail the
    build** via `panic!` with a descriptive message (matches Debian
    reproducible-builds guidance — silent fallback to `now()` would break
    the reproducibility contract without warning). If the env var is
    unset, use `SystemTime::now()`. Format as `YYYY-MM-DDTHH:MM:SSZ` via
    a tiny inline formatter (no `chrono` dependency).
  * `VERSION_BODY` — single pre-formatted line:
    ```
    <short-hash> — <os>/<arch> — built <timestamp>
    ```
    (em-dash for ruSTAR parity; no SIMD label — we are not shipping SIMD
    variants for v2.1.0).
- Read `CARGO_CFG_TARGET_OS` and `CARGO_CFG_TARGET_ARCH` from build-script env.
- Emit rebuild triggers:
  * `cargo:rerun-if-changed=.git/HEAD`
  * `cargo:rerun-if-changed=.git/index`
  * `cargo:rerun-if-env-changed=SOURCE_DATE_EPOCH`

  Without these, `BUILD_TIMESTAMP` would bake in at first build and never
  update. With only the first two, a fresh commit triggers a rebuild.
- **Do not** append `-dirty` — matches ruSTAR; keeps the hash stable across
  uncommitted edits during development.

**`Cargo.toml`**:

- No `[build-dependencies]` section (we only use `std`).
- No explicit `build = "build.rs"` needed — Cargo auto-detects `build.rs` at
  the package root. Call this out in the commit message so reviewers know it's
  intentional.
- `exclude` list: no change. `build.rs` is neither in `exclude` nor does it
  need to be — we want it shipped in the crate tarball so `cargo install
  trim-galore` gets a version string (with `GIT_SHORT_HASH=unknown`).

**`src/cli.rs`** — replace the current:

```rust
#[clap(name = "trim_galore", version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)"), about)]
```

with:

```rust
#[clap(
    name = "trim_galore",
    version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)"),
    long_version = concat!(
        env!("CARGO_PKG_VERSION"), " (Oxidized Edition)\n",
        env!("VERSION_BODY")
    ),
    about
)]
```

Effect: `-V` → one line (unchanged); `--version` → two lines.

**`src/main.rs`** — banner update at lines 25–26:

```rust
eprintln!("\nTrim Galore - Oxidized Edition v{}", env!("CARGO_PKG_VERSION"));
eprintln!("{}", env!("VERSION_BODY"));
eprintln!("==================================================\n");
```

### Step 2 — CI lint gates in `.github/workflows/ci.yml`

**Decision: split into a new `lint` job** (not steps in `rust-tests`).

Rationale:
- Separate failure surface — a lint failure shouldn't obscure a test failure,
  and vice versa. The GitHub checks panel reads cleaner.
- Parallelism — `rust-tests` and `lint` run simultaneously on the same
  ubuntu-latest pool; no wall-clock cost.
- Matches ruSTAR's pattern.

**New job**:

```yaml
lint:
  name: Lint (fmt + clippy)
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@stable
      with:
        components: rustfmt, clippy
    - uses: actions/cache@v4
      with:
        path: |
          ~/.cargo/registry
          ~/.cargo/git
          target
        key: ${{ runner.os }}-cargo-lint-${{ hashFiles('Cargo.lock') }}
        restore-keys: ${{ runner.os }}-cargo-
    - name: cargo fmt --check
      run: cargo fmt --all -- --check
    - name: cargo clippy (release, all targets, -D warnings)
      run: cargo clippy --all-targets --release -- -D warnings
```

- `--all-targets` catches test-code warnings too (e.g. the current
  `len_zero` in `src/trimmer.rs:639`).
- `--release` matches the profile we ship. Yes, this doubles cache footprint,
  but the lint cache key is separate so `rust-tests` isn't blown away.
- Cache key distinct from `rust-tests` to avoid cross-job stomping.

**Job dependency update**: change `validation: needs: [rust-tests]` →
`validation: needs: [rust-tests, lint]`. Rationale: validation is the
expensive job (12-min Conda+Perl install); no point running it if lint fails.

**`release.yml` is intentionally NOT gated on `lint`** — rationale: releases
are cut from tagged commits that have already passed `ci.yml` on the merge
commit; re-running lint at release time would duplicate work without adding
signal. If a hotfix release bypasses the merge-gate path (tag pushed without
a PR), the lint drift would slip — mitigation is to require PRs for all
release-bearing merges (social convention, not enforced here). Document this
split in `release.yml` as a comment near the first job definition.

**Concrete clippy resolutions** (every current warning, enumerated):

| File:line | Lint | Resolution |
|---|---|---|
| `src/parallel.rs:58` | `too_many_arguments` (11/7) | `#[allow(clippy::too_many_arguments)]` on `run_paired_end_parallel` |
| `src/parallel.rs:223` | `too_many_arguments` (8/7) | `#[allow(clippy::too_many_arguments)]` on `process_paired_batch` |
| `src/parallel.rs:300` | `too_many_arguments` (12/7) | `#[allow(clippy::too_many_arguments)]` on `process_pairs` |
| `src/trimmer.rs:332` | `too_many_arguments` (9/7) | `#[allow(clippy::too_many_arguments)]` on `run_paired_end` (matches existing pattern at `src/main.rs:497`) |
| `src/main.rs:159` | `type_complexity` | **Fix** — introduce `type SetupResult = Result<(String, Vec<(String, String)>, Vec<(String, String)>, trimmer::TrimConfig)>` alias above `setup_trimming` |
| `src/main.rs:310` | `type_complexity` | **Fix** — introduce `type ResolvedAdapter = Result<(String, Vec<(String, String)>, Vec<(String, String)>, Option<(usize, usize)>)>` alias above `resolve_adapter` |
| `src/trimmer.rs:163` | `collapsible_if` | **Fix** — collapse into single `if !(paired && is_r2) && len >= 2 && had_adapter { … }` |
| `src/trimmer.rs:69` | `doc_lazy_continuation` | **Fix** — add 3-space indent to doc list item 2.5 |
| `src/demux.rs:232` | `write_with_newline` | **Fix** — `write!(w, "      ¯\\_(ツ)_/¯\n")?` → `writeln!(w, "      ¯\\_(ツ)_/¯")?` |
| `src/report.rs:735` | `int_plus_one` | **Fix** — `range_start <= len - 1` → `range_start < len` |
| `src/alignment.rs:58` | `needless_range_loop` | **Fix** — `for (i, row) in dp.iter_mut().enumerate().take(m + 1).skip(1)`; adjust row indexing. Verify byte-identical parity via the `validation` job. |
| `src/quality.rs:135` | `needless_range_loop` | **Fix** — `for (i, &q) in quals.iter().enumerate().take(n)`; same parity check. |
| `src/trimmer.rs:639` (tests) | `len_zero` | **Fix** — `assert!(result.adapter_matches.len() >= 1)` → `assert!(!result.adapter_matches.is_empty())` |

**Policy**: prefer fix over `#[allow]` where the fix is mechanical and
< 5 lines. Prefer `#[allow]` when the function signature is load-bearing
(the `too_many_arguments` cases — breaking them up would hurt readability).

**Risk**: the two `needless_range_loop` fixes touch hot paths (alignment DP
fill, quality-trim inner loop). Parity is covered by the existing
`validation` job (byte-identical md5 checks vs Perl v0.6.11). If parity
breaks, revert to `#[allow(clippy::needless_range_loop)]` on the function.

### Step 3 — Dependabot config

**New file** `/Users/fkrueger/Github/TrimGalore/.github/dependabot.yml`:

```yaml
version: 2
updates:
  - package-ecosystem: "cargo"
    directory: "/"
    schedule:
      interval: "weekly"
      day: "monday"
    open-pull-requests-limit: 5
    assignees:
      - "FelixKrueger"
    # labels: omitted — see preflight note below

  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "weekly"
      day: "monday"
    open-pull-requests-limit: 5
    assignees:
      - "FelixKrueger"
```

**Label preflight (done during planning)**: `gh label list` confirmed
existing labels: `bug`, `duplicate`, `enhancement`, `help wanted`, `invalid`,
`question`, `wontfix`, `do-not-merge-until-GA`, `beta-feedback`. **`dependencies`
does not exist.**

Decision: ship without `labels:` (Dependabot's default is to auto-create the
`dependencies` label; we can add it ourselves post-merge if we prefer explicit
control). If we want labels, create the `dependencies` label first via `gh
label create dependencies --color 0366d6` before merging this plan — flag
this to the user.

**Day choice**: Monday. Common practice; aligns with a pre-week triage flow.

**Concurrency cap**: `open-pull-requests-limit: 5` per ecosystem. Prevents a
flood on first activation while still surfacing genuine updates.

### Step 4 — `rustsec/audit-check@v2`

**Decision: new job** (`audit`), parallel to `lint` and `rust-tests`.

Rationale:
- Transient advisory-db lookup failures shouldn't red the lint signal.
- A new CVE filed against a released version starts flagging PRs after merge
  — we want that signal visible but not tangled with lint output.

**New job in `.github/workflows/ci.yml`**:

```yaml
audit:
  name: Security audit
  runs-on: ubuntu-latest
  permissions:
    contents: read
    checks: write
    issues: write
  steps:
    - uses: actions/checkout@v4
    - uses: rustsec/audit-check@v2
      with:
        token: ${{ secrets.GITHUB_TOKEN }}
```

- Triggers: inherits `on: push` + `on: pull_request` from the workflow (no
  separate trigger needed).
- `validation` does **not** depend on `audit` — an advisory shouldn't block a
  parity run. (Advisories are signals, not gates.)
- **Caveat**: PRs from forks may hit permission issues with `checks: write` /
  `issues: write`. If we see that, switch the fork case to annotations-only
  or add a `pull_request_target` variant. Note this to watch for.

**Job matrix summary (post-change)**:

```
rust-tests ─┐
lint ───────┼─→ validation
audit ──────┘ (not needed by validation)
```

---

## 5. Verification

### Step 1 (`build.rs`)

**Local**:

```bash
cargo build --release
./target/release/trim_galore --version
# Expected two lines:
#   trim_galore 2.1.0-beta.1 (Oxidized Edition)
#   abc1234 — linux/x86_64 — built 2026-04-19T12:34:56Z

./target/release/trim_galore -V
# Expected one line:
#   trim_galore 2.1.0-beta.1 (Oxidized Edition)

./target/release/trim_galore test_files/illumina_10K.fastq.gz -o /tmp/op
# Startup banner (stderr) should now include the VERSION_BODY line.
```

**Reproducibility** (strong form — assert byte-identical binaries):

```bash
SOURCE_DATE_EPOCH=1700000000 cargo clean && cargo build --release
sha256sum target/release/trim_galore > /tmp/sha_a
SOURCE_DATE_EPOCH=1700000000 cargo clean && cargo build --release
sha256sum target/release/trim_galore > /tmp/sha_b
diff /tmp/sha_a /tmp/sha_b   # must match exactly
```

(The weaker `strings target/release/trim_galore | grep -o ...` check was
considered and rejected — it only proves the timestamp *string* is present,
not that the binary is bit-for-bit identical across runs.)

**Edge cases to exercise**:
- Shallow clone (`git clone --depth 1`): should still produce a short hash.
- `.git` missing + release build:
  ```bash
  mv .git "$TMPDIR/savegit"
  cargo clean && cargo build --release
  ./target/release/trim_galore --version | grep -q '^unknown — '
  mv "$TMPDIR/savegit" .git   # CRITICAL: restore before anything else
  ```
  Must emit `GIT_SHORT_HASH=unknown`, not fail the build.
- `cargo install --path .` fallback smoke (exercises the tarball path that
  crates.io users take):
  ```bash
  mv .git "$TMPDIR/savegit"
  cargo install --path . --locked --force --root "$TMPDIR/tg-install"
  "$TMPDIR/tg-install/bin/trim_galore" --version | grep -q '^unknown — '
  mv "$TMPDIR/savegit" .git
  ```
- Git not on PATH: same — `Command::new("git").output()` returns `Err`, fallback to `"unknown"`.
- Dirty tree: short hash unchanged (we don't add `-dirty`).
- **Malformed `SOURCE_DATE_EPOCH`** should fail cleanly (not panic the
  compiler, but exit with the build-script's descriptive error):
  `SOURCE_DATE_EPOCH=garbage cargo clean && cargo build --release 2>&1 | grep -q 'SOURCE_DATE_EPOCH'`.

**CI assertion**: add a step to the `rust-tests` job that asserts
provenance-line presence (not line count — robust to trailing-newline
quirks in clap's output):

```yaml
- name: Verify --version carries provenance
  run: |
    # Presence check: --version must emit the <hash> — <triple> — built <ts> line
    ./target/release/trim_galore --version \
      | grep -qE '^[0-9a-f]{7,40} — [a-z0-9_]+/[a-z0-9_]+ — built [0-9T:Z-]+$'
    # -V must stay the short form (must NOT carry the provenance line)
    ! ./target/release/trim_galore -V | grep -qE '^[0-9a-f]{7,40} — '
```

### Step 2 (lint gates)

**Local**:

```bash
cargo fmt --all -- --check     # should pass after touch-ups
cargo clippy --all-targets --release -- -D warnings   # zero warnings
```

**CI assertion**: the `lint` job itself — it fails on any warning.

### Step 3 (Dependabot)

**Local**: YAML validation only (`yamllint .github/dependabot.yml` or GitHub's
preview render).

**Post-merge verification**: within ~24h, Dependabot should post a first
status comment on the repo (Insights → Dependency graph → Dependabot). If the
first PR arrives labelled with an auto-created `dependencies` label, we can
either accept that or override via explicit `labels:` in the config.

### Step 4 (audit-check)

**Local**:

```bash
cargo install cargo-audit --locked
cargo audit   # should exit 0 against current Cargo.lock
```

**CI assertion**: the `audit` job itself. New advisories post as PR
annotations.

---

## 6. Risks / edge cases to double-check

1. **`build.rs` on minimal runners**: Ubuntu runners have `git` in PATH, so
   fine. Debian slim / scratch images (used in the Docker build stage) — our
   `Dockerfile` uses `rust:1.85-bookworm` which has git; confirmed safe.
   Users running `cargo install trim-galore` from crates.io get
   `GIT_SHORT_HASH=unknown` — this is intentional and documented in the
   fallback path.

2. **Clippy `--release` doubles CI time for the lint job**: independent
   cache key mitigates this, but first-run after cache eviction will be
   slow (~4-5 min instead of ~30 s). Acceptable — lint runs in parallel
   with `rust-tests`. If it becomes a pain, drop `--release` from clippy
   and rely on the `rust-tests` job's release build for release-profile
   coverage.

3. **`needless_range_loop` fixes in hot paths** (`alignment.rs`, `quality.rs`):
   the `validation` job's byte-identical md5 checks against Perl v0.6.11
   catch any behavioural drift. If parity breaks, revert to
   `#[allow(clippy::needless_range_loop)]` on the function and move on.

4. **Dependabot first-activation flood**: with 7 direct cargo deps + N
   transitive pinned ones + ~6 GitHub Actions pins, expect 3-8 PRs in the
   first week. Cap of 5 PRs/ecosystem prevents runaway.

5. **`audit-check` false-positives from unmaintained advisories**: the
   RustSec database flags `RUSTSEC-*` entries including `unmaintained`
   crates. If our dep tree pulls in an unmaintained crate transitively,
   the action will annotate on every push. Mitigation: if it becomes
   noise, add an `audit.toml` ignore-list at repo root (e.g. `ignore =
   ["RUSTSEC-2024-XXXX"]`) committed alongside the justification in a
   comment.

6. **`SOURCE_DATE_EPOCH` parsing**: **hard-fail on malformed input** via
   `panic!("SOURCE_DATE_EPOCH must be a u64 seconds-since-epoch, got {env:?}")`
   — this is a deliberate choice against silent fallback. Rationale: the
   Debian reproducible-builds spec requires the variable to hold a
   non-negative decimal integer; a caller who sets it but passes garbage
   has a scripting bug, and a silent fallback to `SystemTime::now()` would
   produce a non-reproducible binary while *appearing* to honour the
   request. A panic in `build.rs` surfaces cleanly as a compile error with
   our message, not an internal compiler crash. Unset env var → normal
   `now()` path (unchanged).

7. **`long_version` format mismatch with `-V`**: clap treats `version` as
   the short form and `long_version` as the expanded form, used by
   `--version` only. Verify both paths in local test before push.

8. **Dependabot label preflight** (flag for user): `dependencies` label
   doesn't exist. Decision needed: (a) ship without labels (default
   Dependabot behaviour auto-creates `dependencies`), or (b) create the
   label first via `gh label create dependencies --color 0366d6`. I
   default to (a) — simpler.

9. **`audit-check` + fork PRs**: `checks: write` / `issues: write` may
   be restricted for fork-sourced PRs. Monitor; switch to
   `pull_request_target` or drop to annotations-only if it breaks.

10. **Banner line added to stderr**: external parsers that scrape the
    stderr of `trim_galore` (unlikely but possible in older pipelines)
    will see an extra line. Non-breaking — the banner is strictly
    additive and the `==…==` separator stays where it was.

11. **`release.yml` not gated on `lint`** (intentional split). Risk:
    a hotfix tag pushed outside the normal PR merge flow could ship with
    lint drift. Mitigation: Felix is the sole release cutter, and
    release cuts happen from the same commit that just passed `ci.yml`
    in a PR. If this ever changes (team grows, release automation
    rewires), revisit and add `lint` as a `check-release` dependency in
    `release.yml`.

---

## 7. Implementation order

1. `build.rs` + `src/cli.rs` + `src/main.rs` banner — feature-complete,
   locally testable.
2. All clippy resolutions (one commit per category: allows, type aliases,
   mechanical fixes, loop refactors). Each commit keeps `cargo clippy
   --all-targets --release -- -D warnings` cleaner than the last.
3. `.github/workflows/ci.yml` — add `lint` job + `audit` job + update
   `validation: needs:`. Gate flips once all prior commits are green.
4. `.github/dependabot.yml` — last, because it triggers downstream PR
   activity that's easier to handle after the feature branch is merged.
