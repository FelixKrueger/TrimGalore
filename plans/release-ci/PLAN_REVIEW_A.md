# Plan Review A: Release CI -- Automated Binaries & Containers

**Reviewer**: Plan Reviewer A  
**Date**: 2026-04-13  
**Plan file**: `plans/release-ci/plan.md`

---

## 1. Logic Review

### 1.1 Trigger strategy has an internal contradiction

The plan text states two different trigger designs that conflict with each other:

- **Section "Trigger Strategy"** (line 73-74): `push` triggers only on `optimus_prime`.
- **Resolved Decision #2** (line 158-161): States push triggers on **both** `optimus_prime` and `master`, with `workflow_dispatch` restricted to `master`.

The YAML snippet only shows `optimus_prime`. The resolved decision adds `master`. The implementation needs to pick one and be consistent. If `workflow_dispatch` is intended to be restricted to `master` only (as Decision #2 implies), then the YAML trigger block needs a `branches` filter under `workflow_dispatch` -- but GitHub Actions `workflow_dispatch` does not support a `branches` restriction at the trigger level (it runs on whatever branch the user selects in the UI, defaulting to the repo default branch). The plan should clarify how to enforce that releases only happen from `master`.

### 1.2 macOS x86_64 cross-compilation is not "native"

The build matrix lists both macOS targets on `macos-latest`:

| Name | OS | Target |
|------|----|--------|
| macos-x86_64 | macos-latest | x86_64-apple-darwin |
| macos-aarch64 | macos-latest | aarch64-apple-darwin |

As of 2025, `macos-latest` on GitHub Actions resolves to an Apple Silicon (ARM64) runner. Building `x86_64-apple-darwin` on an ARM64 runner IS cross-compilation. The plan states "each arch builds on its own native runner, no cross-compilation needed" -- this is incorrect for the macOS x86_64 case.

This matters because:
- Cargo will cross-compile, but the resulting binary cannot be tested on the runner (no Rosetta in CI by default).
- The existing `release.yml` correctly uses `macos-13` (Intel) for x86_64 and `macos-14` (ARM) for aarch64.
- The plan should either (a) use `macos-13` for x86_64 (matching the current workflow), or (b) acknowledge this is cross-compilation and accept that the binary cannot be smoke-tested in CI.

Note: GitHub is deprecating Intel macOS runners. If `macos-13` becomes unavailable, cross-compiling on ARM is the fallback, but the plan should acknowledge this explicitly.

### 1.3 No smoke test for built binaries

The current `release.yml` has a `test-install` job that runs `cargo install --path .` and verifies `trim_galore --help`. The new plan drops this entirely. For dev pushes, untested binaries are uploaded as artifacts and baked into the `:dev` Docker image. For releases, untested binaries are attached to a GitHub Release.

At minimum, each native binary should be smoke-tested with `--help` or `--version` before upload. Cross-compiled binaries (macOS x86_64 if built on ARM) cannot be tested, but that should be documented as a known limitation.

### 1.4 Docker build uses `cargo build --release` without specifying target

The Dockerfile runs `cargo build --release` without a `--target` flag. This builds for the host architecture, which is correct for native runners. However, the plan should confirm that the Docker build stage inherits the platform from `docker/build-push-action`'s `--platform` flag (it does via BuildKit), so this is fine -- but worth calling out to avoid confusion since the binary build jobs DO specify targets explicitly.

### 1.5 Missing `--locked` flag on cargo build

Neither the Dockerfile nor the binary build jobs specify `cargo build --release --locked`. Without `--locked`, Cargo may update `Cargo.lock` during the build, potentially pulling in different dependency versions than what was tested in CI. This is a correctness risk.

### 1.6 Artifact retention policy is underspecified

The plan says "retention: 1 day for dev, until release for releases" but doesn't specify how this is implemented. GitHub Actions `upload-artifact` has a `retention-days` parameter. "Until release" is not a built-in concept -- artifacts are retained for a fixed number of days (default 90). The plan should specify a concrete retention period for release artifacts (e.g., 7 days, since they get attached to the Release anyway).

### 1.7 Job dependency ordering gap

The job graph shows `upload-binaries` depends on `create-tag-and-release` (needs the release to exist) and on `build-binaries` (needs the artifacts). But this dependency chain is not explicitly stated. If both are listed as `needs: [check-release]` without depending on each other, `upload-binaries` could run before the release exists.

### 1.8 No version consistency check between Cargo.toml and Cargo.lock

Job 0 (`check-release`) says it will "verify `Cargo.lock` version matches." But `Cargo.lock` does not have a top-level "version" field that mirrors `Cargo.toml`. The lock file has a `version = 4` (lock file format version) and individual package entries. The check should verify that the `trim-galore` entry in `Cargo.lock` matches the version in `Cargo.toml`, which requires parsing the lock file. This is a minor point but the plan should be precise about what is being checked.

---

## 2. Assumptions

### 2.1 Stated assumptions -- validated

| Assumption | Status |
|------------|--------|
| `edition = "2024"` requires Rust 1.85+ | Correct. Cargo.toml has `rust-version = "1.85"` and the Dockerfile uses `rust:1.85-bookworm`. |
| Binary name is `trim_galore` | Correct. `[[bin]] name = "trim_galore"` in Cargo.toml. |
| `procps` is needed for Nextflow | Correct. Nextflow uses `ps` for process monitoring. |
| No C build dependencies | Correct. Dependencies are pure Rust (flate2 with zlib-rs, gzp with deflate_rust). No system libraries needed. |
| GHCR uses `GITHUB_TOKEN` | Correct for public repos. |

### 2.2 Implicit assumptions -- surfaced

1. **`ubuntu-24.04-arm` runner availability**: The plan assumes GitHub Actions provides `ubuntu-24.04-arm` runners. This is a relatively new runner type. If the repo is under a free plan, ARM runners may not be available (they are currently only in public beta for public repos). The plan should note this dependency.

2. **`rust:1.85-bookworm` Docker image exists**: This image should exist on Docker Hub since Rust 1.85 is released, but the plan hardcodes it rather than using, e.g., `rust:1.85` (which defaults to bookworm). Minor, but pinning to `bookworm` is good practice.

3. **GHCR package naming**: `ghcr.io/felixkrueger/trimgalore` assumes the package name will be `trimgalore` (lowercase). GitHub Container Registry lowercases the owner but the package name comes from the image tag pushed. This is fine as long as all push commands consistently use this name.

4. **Cargo.lock is committed**: The plan does not mention `Cargo.lock`. For a binary crate, `Cargo.lock` should be committed (it is -- verified). But the Dockerfile should `COPY Cargo.lock .` or the full `COPY . .` must include it, which it does.

5. **No glibc version pinning**: The plan chooses `gnu` targets and states "glibc 2.31+, Ubuntu 20.04+". The `debian:bookworm-slim` runtime image has glibc 2.36. Binaries built on `ubuntu-latest` (currently Ubuntu 22.04, glibc 2.35) will link against glibc 2.35 and won't run on older distros. This is acceptable for the stated audience but should be documented.

6. **`workflow_dispatch` with no inputs**: The plan specifies `workflow_dispatch: {}` with no inputs. This means no way to override the version, skip steps, or do a dry run. Consider whether an `input` for version override or dry-run mode would be useful.

---

## 3. Efficiency Analysis

### 3.1 Docker build efficiency

The Dockerfile does `COPY . .` which copies the entire repo into the build stage. With `lto = true` and `codegen-units = 1`, release builds are slow (potentially 5-10 minutes). Every push to `optimus_prime` triggers this. The plan uses GHA cache for Docker layers, but `COPY . .` invalidates the cache on any file change, meaning the full `cargo build --release` runs every time.

**Recommendation**: Use a two-step COPY pattern to leverage Docker layer caching:
```dockerfile
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main(){}' > src/main.rs && cargo build --release && rm -rf src
COPY . .
RUN cargo build --release
```
This caches the dependency build layer and only rebuilds the application code. This can save 3-5 minutes per build.

### 3.2 Building binaries AND Docker on every push is expensive

Every push to `optimus_prime` triggers 4 binary builds + 2 Docker builds = 6 parallel jobs. For a development branch that may see multiple pushes per day, this is expensive in CI minutes. The Docker `:dev` image is the stated goal for nf-core testing -- the binary artifacts are less useful for dev pushes.

**Recommendation**: Consider building only Docker on push to `optimus_prime`, and building binaries only on `workflow_dispatch` (release). Or gate binary builds behind a path filter that only triggers when `src/` or `Cargo.*` change.

### 3.3 No concurrency control

Multiple rapid pushes to `optimus_prime` will trigger multiple workflow runs that all push to `ghcr.io/.../trimgalore:dev`. These can race and the final `:dev` tag may not correspond to the latest push.

**Recommendation**: Add a `concurrency` group to cancel in-progress runs:
```yaml
concurrency:
  group: release-${{ github.ref }}
  cancel-in-progress: true
```

### 3.4 Release LTO build time

With `lto = true` and `codegen-units = 1` in `[profile.release]`, each binary build will be slow. This is fine for releases but may be unnecessary for dev pushes. The plan does not differentiate build profiles between dev and release modes.

---

## 4. Validation Sufficiency

### 4.1 No binary validation in the workflow

The biggest validation gap: no built binary is ever executed in the workflow. The current release.yml has `test-install` + `trim_galore --help`. The new plan removes this. A corrupted or non-functional binary could be released.

### 4.2 No Docker image validation

The plan says to manually "Pull container, run `trim_galore --help`" in Implementation Order step 4, but the workflow itself doesn't include an automated smoke test of the Docker image. A `docker run ghcr.io/.../trimgalore:dev --help` step after the merge would catch runtime issues (missing libraries, wrong entrypoint, etc.).

### 4.3 No validation that the Docker image actually works with real data

The plan mentions testing with nf-core/rnaseq as step 5 of implementation, but this is manual. Consider adding a minimal smoke test job that runs the container against a tiny test FASTQ to catch issues like missing shared libraries or broken I/O.

### 4.4 SHA256 checksums mentioned but no verification

The plan mentions `.sha256` files alongside tarballs, but doesn't describe how they're generated or verified. The implementation should generate them with `sha256sum` and ideally verify them after download in the upload step.

### 4.5 Tag existence check may race

Job 0 checks that the tag doesn't exist, and Job 4 creates it. Between these two steps, another process could create the tag. This is a very unlikely race but the `create-tag-and-release` job should handle the "tag already exists" error gracefully.

---

## 5. Alternatives Considered

### 5.1 Use `cross` for cross-compilation instead of native ARM runners

ARM runners (`ubuntu-24.04-arm`) are new and may have availability issues. The `cross` tool (or `cargo-cross`) with QEMU can build ARM binaries on x86_64 runners. The trade-off is slower builds (~3-5x) but guaranteed runner availability. Given that native ARM runners are available for public repos, the plan's approach is fine, but `cross` should be documented as the fallback.

### 5.2 Single workflow vs. separate ci.yml and release.yml

The plan keeps ci.yml (tests) and release.yml (builds/releases) separate. This is clean separation of concerns. An alternative would be to have the release workflow depend on CI passing first (via `workflow_run` trigger or required status checks). Currently, a push to `optimus_prime` triggers both CI and release independently -- if CI fails but release succeeds, a broken `:dev` image is pushed.

**Recommendation**: Add `needs`-like gating, or at minimum run `cargo test` before `cargo build --release` in the binary build jobs. Or use a `workflow_run` trigger so release only runs after CI passes.

### 5.3 Use GitHub Attestations for supply chain security

GitHub now supports artifact attestations (`actions/attest-build-provenance@v2`). For a bioinformatics tool, supply chain provenance can be valuable. This is optional but forward-looking.

### 5.4 Use `cargo-dist` for release automation

[cargo-dist](https://opensource.axo.dev/cargo-dist/) automates binary releases for Rust projects, including cross-compilation, checksums, GitHub Releases, and installers. It would replace most of the custom workflow. Trade-off: less control but much less maintenance burden. Worth evaluating.

---

## 6. Action Items

### Critical

1. **Fix macOS x86_64 runner**: Either use `macos-13` for x86_64 builds (matching current workflow) or acknowledge cross-compilation and document that the binary cannot be smoke-tested. (`macos-latest` is ARM64.)

2. **Add `--locked` to all cargo build commands**: Both in the Dockerfile and in binary build jobs. Prevents dependency drift between CI test and release build.

3. **Add binary smoke test**: At minimum, run `trim_galore --help` (or `--version`) on each natively-built binary before uploading. This catches linking errors, missing dependencies, and startup crashes.

4. **Resolve trigger contradiction**: The trigger YAML and Resolved Decision #2 disagree on whether `master` push triggers the workflow. Pick one and update both sections.

5. **Gate release workflow on CI passing**: Either add `cargo test` to the release workflow, or use `workflow_run` to trigger release only after CI succeeds. Without this, a broken `:dev` image can be pushed to GHCR.

### Important

6. **Add Docker layer caching optimization**: Use the two-step COPY pattern in the Dockerfile to avoid rebuilding dependencies on every push. Saves significant CI time.

7. **Add concurrency control**: Use `concurrency` groups with `cancel-in-progress: true` to prevent race conditions on rapid pushes.

8. **Add Docker smoke test job**: After `docker-merge`, run `docker run ghcr.io/.../trimgalore:dev --help` to validate the image works end-to-end.

9. **Specify concrete artifact retention**: Replace "until release" with a concrete number (e.g., 7 days for releases, 1 day for dev).

10. **Clarify `upload-binaries` job dependencies**: Explicitly state it depends on both `build-binaries` and `create-tag-and-release`.

### Optional

11. **Consider building binaries only on release**: Dev pushes only need the Docker image for nf-core testing. Skipping 4 binary builds on every push saves CI minutes.

12. **Add `workflow_dispatch` inputs**: A `dry_run` boolean and/or version override input would make the release process more flexible and safer.

13. **Evaluate `cargo-dist`**: May eliminate most of the custom workflow with less maintenance.

14. **Document glibc version floor**: State explicitly that Linux binaries require glibc 2.35+ (Ubuntu 22.04+).

15. **Consider GitHub Attestations**: `actions/attest-build-provenance@v2` for supply chain security on release binaries.

---

## Summary

The plan is well-structured and makes sound architectural decisions (GHCR, gnu targets, no SIMD, multi-stage Docker). The major issues are: (1) macOS x86_64 builds on an ARM runner without acknowledging cross-compilation, (2) no binary or Docker smoke tests in the workflow, (3) missing `--locked` flag risking dependency drift, (4) contradictory trigger specifications, and (5) no gating on CI success before pushing dev containers. The efficiency recommendations around Docker layer caching and concurrency control would meaningfully improve the development experience.
