# Release Plan: Trim Galore - Oxidized Edition v2.1.0

## Goal

Release the Rust rewrite as **Trim Galore v2.1.0** — a drop-in replacement for the Perl original. Users installing via bioconda, crates.io, or downloading binaries should feel nothing but a speed boost. Same binary name (`trim_galore`), same CLI flags, same output filenames, same report format.

## Guiding Principle

**Seamless upgrade.** Existing nf-core pipelines, Snakemake workflows, and shell scripts that call `trim_galore` should work without modification. The only visible difference: it's faster and has no external dependencies (no Python, no Cutadapt, no pigz).

## Rollout Strategy

Staged rollout — build confidence before broad distribution:

1. **Repo restructure + crates.io** — make it installable via `cargo install trim-galore`
2. **Author validation** — Felix runs it on real-world datasets (WGBS, RRBS, RNA-seq, etc.), comparing outputs byte-for-byte against Perl TG on data he knows well
3. **Beta testing** — share with trusted testers (Phil, nf-core contributors), gather feedback across instruments, library types, and operating systems
4. **GitHub Release** — prebuilt binaries for Linux/macOS once confidence is established
5. **Bioconda** — update the `trim-galore` recipe only after Stages 2-3 confirm it's solid

Bioconda is the last step, not the first. A broken bioconda release affects thousands of users — the cost of a few weeks of personal validation is negligible by comparison.

---

## Phase 1: Repository Restructuring

### 1.1 Promote Rust code to repo root

Currently the Rust project lives in `optimus_prime/`. For a proper release, Cargo.toml should be at the repo root (standard Rust convention, required for clean `cargo install`).

**Changes:**
- Move `optimus_prime/src/` → `src/`
- Move `optimus_prime/Cargo.toml` → `Cargo.toml`
- Move `optimus_prime/Cargo.lock` → `Cargo.lock` (must remain committed — binary crate convention)
- Archive Perl script: `trim_galore` → `legacy/trim_galore` (keep for reference, not in PATH)
- Move `optimus_prime/SUMMARY.md` → `docs/SUMMARY.md` or similar
- Move benchmark PNGs to `docs/benchmarks/`
- Update `.gitignore` for Rust target directory at root (do NOT gitignore `Cargo.lock`)
- **Update both CI workflows** (`.github/workflows/ci.yml` and `release.yml`) to remove all `optimus_prime/` path prefixes — `working-directory`, cache paths, binary paths

### 1.2 Update Cargo.toml

```toml
[package]
name = "trim-galore"
version = "2.1.0"
edition = "2021"
rust-version = "1.74"
description = "Fast adapter and quality trimming for NGS data — Oxidized Edition"
license = "GPL-3.0-only"
authors = ["Felix Krueger"]
repository = "https://github.com/FelixKrueger/TrimGalore"
homepage = "https://github.com/FelixKrueger/TrimGalore"
keywords = ["bioinformatics", "ngs", "adapter-trimming", "fastq", "quality-trimming"]
categories = ["command-line-utilities", "science"]
readme = "README.md"
exclude = [
    "test_files/",
    "plans/",
    "docs/",
    "legacy/",
    ".github/",
    ".claude/",
]

[[bin]]
name = "trim_galore"
path = "src/main.rs"

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
strip = "debuginfo"
```

**Notes:**
- `exclude` is inside `[package]` (Cargo rejects it elsewhere)
- `license = "GPL-3.0-only"` — matches the LICENSE file (standard FSF GPL v3 text, not "or later")
- `strip = "debuginfo"` — strips DWARF debug info (saves most size) but preserves symbol names for backtraces in crash reports
- Crate name `trim-galore` (hyphenated, crates.io convention) produces binary `trim_galore` (underscored, as specified in `[[bin]]`)
- **Pre-requisite:** Verify `trim-galore` is available on crates.io before proceeding

### 1.3 Version string and branding

Apply consistent branding across all surfaces:

| Surface | Current | Target |
|---------|---------|--------|
| `--version` output | `trim_galore 0.1.0` | `Trim Galore - Oxidized Edition 2.1.0` |
| Startup banner (main.rs) | `Trim Galore v0.1.0` | `Trim Galore - Oxidized Edition v2.1.0` |
| Report header (report.rs) | `Trim Galore version: 0.1.0 (Optimus Prime)` | `Trim Galore version: 2.1.0 (Oxidized Edition)` |
| Clap `about` text | No edition mention | Include "Oxidized Edition" |

**Critical:** Verify MultiQC's version-line regex before changing the report header format. Check `multiqc/modules/trim_galore/trim_galore.py` for the parsing pattern. The version line must remain parseable.

### 1.4 Deprecated CLI flags (no-op with warning)

For drop-in compatibility, the Rust version must accept all Perl CLI flags — even those that no longer apply. Flags that have no Rust equivalent should be accepted silently with a deprecation warning printed to stderr:

```
WARNING: --path_to_cutadapt is deprecated in Trim Galore v2.0 (no external Cutadapt needed). Ignoring.
```

**Flags requiring no-op handling** (confirmed by full flag audit):

| Perl flag | Deprecation message |
|---|---|
| `--gzip` | "Output is gzipped by default in v2.0. Use --dont_gzip to disable." |
| `--path_to_cutadapt` | "No longer needed — Cutadapt is built in." |
| `--cutadapt_args` | "No longer needed — Cutadapt is built in." |
| `--suppress_warn` | "No longer needed — no Cutadapt subprocess." |
| `--report` | Silent no-op (reports are on by default; opposite of --no_report_file) |
| `--keep` | "Not yet supported in v2.0. RRBS reads below length cutoff will be removed." |
| `--hulu` | Silent no-op (Easter egg) |

**Aliases to add:** `--casio` and `--breitling` as aliases for `--clock` (Easter egg aliases from Perl version).

**Compatibility fix needed:** `--implicon` is a boolean flag in Perl but takes a `usize` argument (UMI length) in Rust. Make the value optional with a sensible default (e.g., 8, matching the original hardcoded value) so bare `--implicon` without a value still works.

Implementation: Add these as hidden clap arguments that print a warning to stderr on use.

---

## Phase 2: Pre-release Validation

### 2.1 Feature completeness audit

Complete diff of every CLI flag from Perl `trim_galore --help` vs Rust implementation:
- All flags accepted (matched, deprecated no-op, or documented as removed)
- Output filenames identical for all modes
- Report format compatible with MultiQC (verify regex parsing)
- Exit codes match

### 2.2 Byte-identical validation

Run the full test matrix on the final release binary (with LTO):
- All test files in `test_files/`
- Single-end, paired-end, RRBS, NextSeq, poly-A, poly-G, hardtrim, clock, implicon, demux
- Verify decompressed output is byte-identical to Perl TG (with Cutadapt >=4.4)
- Run with `--cores 1`, `--cores 4`, `--cores 8` to verify deterministic output

### 2.3 MultiQC compatibility

- Run MultiQC on trimming reports from both Perl TG and Oxidized Edition
- Verify MultiQC parses both identically
- Check report header format matches what MultiQC expects

### 2.4 Edge cases

- Empty input files
- Truncated gzip files (should error gracefully)
- Colorspace files (should reject with clear message)
- Very short reads (1-2bp)
- Single read in file
- Mismatched paired-end files

---

## Phase 3: CI/CD Infrastructure

### 3.1 CI workflow (`.github/workflows/ci.yml`)

Runs on every push and PR:
- `cargo test` — all unit tests
- `cargo clippy -- -D warnings` — lint
- `cargo fmt --check` — formatting
- MSRV check — separate job using `dtolnay/rust-toolchain@1.74.0` with `cargo check`
- Test on: Linux x86_64, macOS aarch64
- **Expanded validation matrix:** Add RRBS (directional + non-directional), NextSeq, poly-A, poly-G, `--retain_unpaired`, `--clip_R1/R2` to the byte-identical validation job

### 3.2 Release workflow (`.github/workflows/release.yml`)

Triggered on `workflow_dispatch` (dispatched from `master` for GA or `optimus_prime` for prereleases). Sequence: bump version in Cargo.toml → commit Cargo.toml + Cargo.lock → (for GA) merge PR to master → dispatch Release workflow → workflow creates and pushes the tag.

**Step 1: Build binaries**

| Platform | Target | Runner | Notes |
|----------|--------|--------|-------|
| Linux x86_64 | `x86_64-unknown-linux-musl` | `ubuntu-latest` | Static binary, no glibc dependency |
| Linux aarch64 | `aarch64-unknown-linux-musl` | `ubuntu-24.04-arm` | **Native ARM runner** (no cross-compilation) |
| macOS x86_64 | `x86_64-apple-darwin` | `macos-13` | Intel runner |
| macOS aarch64 | `aarch64-apple-darwin` | `macos-14` | Native Apple Silicon |

Each produces: `trim_galore-{version}-{target}.tar.gz` + `.sha256` (add `sha256sum` step to workflow)

**Step 2: Create GitHub Release**
- Create GitHub Release with auto-generated release notes
- Upload all tarballs + SHA256 checksums

**Step 3: Publish to crates.io** (add to workflow)
- Use `rust-lang/crates-io-auth-action` for OIDC trusted publishing
- `cargo publish`
- Alternative: manual `cargo publish` for first release, automate later

---

## Phase 4: Bioconda Recipe

### 4.1 Strategy: prebuilt binaries vs source build

Bioconda supports both approaches. Two options:

**Option A: Prebuilt binaries** (simpler, following RustQC pattern)
- Download platform-specific tarballs from GitHub Releases
- Trivial `build.sh` (just `cp trim_galore $PREFIX/bin/`)
- Skip lint: `should_be_noarch_generic`

**Option B: Build from source** (bioconda-preferred)
- Use `{{ compiler('rust') }}` in build requirements
- `cargo build --release` in `build.sh`
- No prebuilt binary dependency, bioconda handles linking
- Slower builds, potential cross-compilation issues in bioconda CI

**Decision: Start with Option A** (prebuilt binaries), as RustQC has validated this path with bioconda maintainers. Fall back to Option B if the PR receives pushback.

### 4.2 Recipe (prebuilt binary approach)

```yaml
{% set version = "2.1.0" %}

package:
  name: trim-galore
  version: '{{ version }}'

source:
  - url: https://github.com/FelixKrueger/TrimGalore/releases/download/v{{ version }}/trim_galore-{{ version }}-x86_64-unknown-linux-musl.tar.gz  # [linux and x86_64]
    sha256: PLACEHOLDER
  - url: https://github.com/FelixKrueger/TrimGalore/releases/download/v{{ version }}/trim_galore-{{ version }}-aarch64-unknown-linux-musl.tar.gz  # [linux and aarch64]
    sha256: PLACEHOLDER
  - url: https://github.com/FelixKrueger/TrimGalore/releases/download/v{{ version }}/trim_galore-{{ version }}-x86_64-apple-darwin.tar.gz  # [osx and x86_64]
    sha256: PLACEHOLDER
  - url: https://github.com/FelixKrueger/TrimGalore/releases/download/v{{ version }}/trim_galore-{{ version }}-aarch64-apple-darwin.tar.gz  # [osx and arm64]
    sha256: PLACEHOLDER

build:
  number: 0

requirements:
  build:
  run:
    # No runtime dependencies — single static binary
    # FastQC is optional (only needed with --fastqc flag)

test:
  commands:
    - trim_galore --version
    - trim_galore --help | grep "Oxidized Edition"

about:
  home: https://github.com/FelixKrueger/TrimGalore
  license: GPL-3.0-only
  license_family: GPL
  license_file: LICENSE
  summary: >
    Fast adapter and quality trimming for NGS data.
    Oxidized Edition — complete Rust rewrite, drop-in replacement.
  dev_url: https://github.com/FelixKrueger/TrimGalore

extra:
  additional-platforms:
    - linux-aarch64
    - osx-arm64
  identifiers:
    - usegalaxy-eu:trim_galore
  skip-lints:
    - should_be_noarch_generic  # prebuilt platform-specific binaries
```

**build.sh:**
```bash
#!/bin/bash
mkdir -p $PREFIX/bin
cp trim_galore $PREFIX/bin/
chmod +x $PREFIX/bin/trim_galore
```

### 4.3 Dependency impact

| | v0.6.11 (Perl) | v2.1.0 (Rust) |
|---|---|---|
| Runtime deps | perl, cutadapt ==5.2, python >=3, fastqc | (none required) |
| Optional deps | pigz | fastqc (for --fastqc) |
| Solver impact | Pins Python + Cutadapt versions | Zero Python dependency |

Removing the `cutadapt ==5.2` pin alone eliminates a major source of conda solver conflicts for users.

---

## Phase 5: crates.io Publication

### 5.1 Pre-publish checklist

- [x] Verify `trim-galore` crate name is available — **confirmed available** (2026-04-13)
- [ ] `cargo package --list` — review what gets included, verify no large binary files
- [ ] `cargo publish --dry-run` — verify it builds clean
- [ ] Set up OIDC trusted publishing on crates.io (link to GitHub repo)

### 5.2 Ongoing

After initial publish, each new version:
- Bump version in Cargo.toml
- `cargo publish` via release workflow (or manually)

---

## Phase 6: Documentation & Communication

### 6.1 README.md overhaul

- Lead with "Trim Galore - Oxidized Edition" branding
- Installation: `cargo install trim-galore` / bioconda / download binary
- Quick start (same as before — the CLI is identical)
- "Upgrading from v0.6.x" section: document it's a drop-in replacement; list deprecated flags that now print warnings; note `conda install trim-galore==0.6.11` for rollback
- Performance highlights with benchmark charts
- New features (poly-G auto-detection, poly-A trimming)
- Link to detailed SUMMARY.md for benchmarks

### 6.2 CHANGELOG.md

Create a changelog with the v2.1.0 entry. The release workflow extracts this for GitHub Release notes.

### 6.3 nf-core module update

The `nf-core/modules` trimgalore module currently:
- Installs via bioconda (will automatically pick up v2.1.0)
- Calculates `-j` as `task.cpus - 4` (paired) or `task.cpus - 3` (single)
- The `-4` subtraction exists to reserve cores for pigz/Cutadapt subprocess overhead

**For v2.1.0:** The Oxidized Edition can use `--cores` equal to `task.cpus` directly since there's no subprocess overhead (single process, N+4 threads from N cores). The module should be updated to pass the full CPU allocation instead of subtracting 4.

**Decided:** Update separately, after bioconda + biocontainer are available. Felix to coordinate with nf-core team. The current `-j (cpus-4)` still works — just underutilizes slightly.

---

## Phase 7: Staged Release Checklist

### Stage 1 — Repo + crates.io (make it installable)
- [x] Crate name `trim-galore` verified available on crates.io (confirmed 2026-04-13)
- [ ] Repo restructured (Cargo.toml at root, Perl script archived in `legacy/`)
- [ ] CI workflows updated for new paths (no more `optimus_prime/` prefix)
- [ ] Branding consistent: `--version`, startup banner, report header all say "Oxidized Edition"
- [ ] Deprecated Perl flags accepted as no-ops with warnings
- [ ] All tests pass with release binary (LTO)
- [ ] MultiQC version-line parsing verified
- [ ] Cargo.toml metadata complete (`exclude`, license, description)
- [ ] `cargo publish --dry-run` succeeds
- [ ] Publish to crates.io → `cargo install trim-galore` works
- [ ] CI workflow passing: tests, clippy, fmt, MSRV check

### Stage 2 — Author validation (weeks, not days)
- [ ] Felix runs on real WGBS datasets, compares output to Perl TG
- [ ] Felix runs on real RRBS datasets (directional + non-directional)
- [ ] Felix runs on RNA-seq / other library types
- [ ] Test with different instruments (NovaSeq, NextSeq, HiSeq, etc.)
- [ ] Verify poly-G auto-detection fires correctly on 2-colour data
- [ ] MultiQC compatibility confirmed on real reports
- [ ] Test at various `--cores` values on real-world file sizes
- [ ] Edge cases: very large files, many samples, unusual read lengths

### Stage 3 — Beta testing (trusted users)
- [ ] Share with Phil / nf-core contributors
- [ ] Gather feedback across different OS / platforms
- [ ] Address any issues found
- [ ] README.md updated
- [ ] CHANGELOG.md written

### Stage 4 — GitHub Release (prebuilt binaries)
- [ ] Update release workflow: native ARM runner for aarch64, add SHA256 generation
- [ ] Trigger release via tag push (`git tag v2.1.0 && git push --tags`)
- [ ] Verify GitHub Release created with binaries + checksums for all 4 platforms
- [ ] Download and test binaries on Linux + macOS
- [ ] Verify `cargo publish` succeeds (automated or manual)

### Stage 5 — Bioconda (only after confidence is established)
- [ ] Validate bioconda recipe locally with `conda-build` before submitting
- [ ] Submit bioconda recipe PR
- [ ] Test bioconda installation after merge
- [ ] Verify `conda install trim-galore` gives v2.1.0
- [ ] Verify `conda install trim-galore==0.6.11` still works (rollback path)
- [ ] Monitor for 2 weeks post-merge (GitHub issues, bioconda feedback)

### Post-release
- [ ] Announce (Twitter/Mastodon, Biostars, nf-core Slack)
- [ ] Monitor bioconda/GitHub issues for any compatibility reports
- [ ] Open nf-core module PR (update `--cores` allocation, remove Cutadapt/pigz deps)

---

## Rollback Plan

If v2.1.0 causes issues after bioconda release:

**Immediate (within hours):**
- Users can pin: `conda install trim-galore==0.6.11` — bioconda retains all historical versions
- Document this in a GitHub issue pinned to the repo

**Short-term (within days):**
- Submit a bioconda recipe PR reverting to v0.6.11 as the latest
- Post advisory on nf-core Slack and Biostars

**Rollback triggers:**
- Any byte-non-identical output on a common use case (SE/PE Illumina, RRBS)
- MultiQC parsing failures on v2.1.0 reports
- Crash or hang on valid input data
- Felix is the decision maker for rollback

**Prevention:**
- Stage 2 author validation is the primary defense (weeks of real-world testing)
- Stage 3 beta testing catches cross-platform and edge-case issues
- The staged rollout means issues are caught before they reach bioconda

---

## Decisions (resolved)

1. **Perl script:** Archive in `legacy/` directory. Preserved in git history and accessible via `v0.6.11` tag.
2. **SIMD tiers:** No — single binary per platform. `zlib-rs` handles SIMD dispatch internally.
3. **nf-core module:** Update separately, after bioconda + biocontainer are available. Felix to coordinate with nf-core team. The current `-j (cpus-4)` still works — just underutilizes slightly.
4. **Linux linking:** Static (musl) for maximum portability across HPC environments. No glibc version dependency.
5. **Docker images:** Not needed — BioContainers automatically builds Docker + Singularity images from every bioconda package. Seqera Wave can also resolve bioconda packages on-the-fly.
6. **Minimum Rust version:** Keep 1.74 (Dec 2023). Only affects `cargo install` users; 2.5 years old is accessible enough.
7. **License:** GPL-3.0-only (matches LICENSE file). Consistent across Cargo.toml, bioconda recipe, and README.
8. **Deprecated flags:** Accept as no-ops with deprecation warnings to stderr. Never hard-error on a flag that worked in v0.6.x.
9. **Bioconda approach:** Prebuilt binaries (Option A), following RustQC precedent. Fall back to source build if PR receives pushback.
10. **Release trigger:** Tag push (`v*`), matching existing workflow. Sequence: version bump → commit → tag → push.
11. **Strip mode:** `strip = "debuginfo"` — preserves symbol names for backtraces while removing bulk debug info.
12. **aarch64 Linux:** Native ARM runner (`ubuntu-24.04-arm`), not cross-compilation.

---

## Review Amendments

Changes incorporated from dual independent review (Reviewers A and B, 2026-04-13):

- Fixed bioconda recipe URLs: `linux-gnu` → `linux-musl` to match release workflow targets (H1/A, H1/B)
- Specified native ARM runner for aarch64-linux builds, removing fragile cross-compilation (H4/A, H2/B)
- Added explicit rollback plan with triggers and timeline (M2/A, H3/B)
- Moved `exclude` inside `[package]` in Cargo.toml template (M3/B)
- Changed `strip = true` to `strip = "debuginfo"` for backtrace preservation (L6/B)
- Added CI workflow path update to Phase 1.1 (M4/A)
- Added clippy, fmt, MSRV check to CI plan (M5/A, L4/A)
- Added SHA256 generation to release workflow (M4/B)
- Expanded CI validation matrix (M7/A, L4/B)
- Clarified license as GPL-3.0-only (M1/B)
- Added branding consistency table (M1/A, M2/B)
- Added deprecated CLI flags section with no-op strategy (M6/B)
- Removed `run_exports` from bioconda recipe (L5/A, L3/B)
- Noted Cargo.lock must be committed (L5/B)
- Documented both bioconda recipe strategies with decision rationale (H1/A)
- Clarified release trigger matches existing workflow (H2/A)
- Added `cargo publish` step to release workflow plan (M3/A)
- Added crate name check as blocking prerequisite (L3/A, M5/B)
