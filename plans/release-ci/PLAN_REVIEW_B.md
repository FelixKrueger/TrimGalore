# Plan Review B: Release CI -- Automated Binaries & Containers

**Reviewer:** Plan Reviewer B
**Plan:** `/plans/release-ci/plan.md`
**Date:** 2026-04-13

---

## 1. Logic Review

### 1.1 Trigger strategy: stated vs. designed conflict

The plan's Trigger Strategy section specifies:

```yaml
on:
  push:
    branches: [optimus_prime]
```

But Resolved Decision #2 says:

> Push to `optimus_prime` -> `:dev` container
> Push to `master` -> `:dev` container (in case of direct pushes/merges)
> `workflow_dispatch` on `master` -> full release

These two contradict. The trigger YAML only has `optimus_prime`, not `master`. The implementation must include both branches in the push trigger or the "push to master -> :dev container" behavior described in Decision #2 will not work. This needs to be reconciled before implementation.

### 1.2 workflow_dispatch has no branch constraint

The plan says `workflow_dispatch: {}` with no branch filter. This means anyone with write access can trigger a full release from *any* branch -- `optimus_prime`, `dev`, a feature branch, etc. Decision #2 says releases should only happen from `master`. There is no guard in the plan to abort `workflow_dispatch` unless the ref is `master`. The `check-release` job should verify `github.ref == 'refs/heads/master'` and fail fast otherwise.

### 1.3 macOS x86_64 cross-compilation is actually needed

The build matrix says:

| Name | OS | Target |
|------|----|--------|
| macos-x86_64 | macos-latest | x86_64-apple-darwin |
| macos-aarch64 | macos-latest | aarch64-apple-darwin |

Both use `macos-latest`, which since late 2024 resolves to Apple Silicon (M-series) runners. Building `x86_64-apple-darwin` on an ARM runner IS cross-compilation. The plan claims "no cross-compilation needed" -- this is incorrect for the macOS x86_64 target. Note the existing release.yml correctly uses `macos-13` (Intel) for x86_64 and `macos-14` (ARM) for aarch64. The plan should either:
- Use `macos-13` for x86_64 (Intel runner), or
- Acknowledge the cross-compilation requirement and add the `x86_64-apple-darwin` target via `rustup target add`.

Since `macos-13` runners are being deprecated (already limited availability), the safer long-term approach is to cross-compile on `macos-latest` with the target explicitly added.

### 1.4 Docker build on native ARM runners: runner availability

The plan specifies `ubuntu-24.04-arm` for both the linux-aarch64 binary build and the linux/arm64 Docker build. GitHub's ARM runners are available for public repos but have historically had lower capacity and occasional availability issues. The plan should note a fallback strategy (e.g., QEMU emulation on amd64 runners) even if native is the primary approach, since build failures from runner unavailability would block every push to `optimus_prime`.

### 1.5 Missing: `Cargo.lock` not committed / not mentioned

The plan mentions "verify `Cargo.lock` version matches" in the check-release job, but does not mention whether `Cargo.lock` is committed. For binary crates, the Rust convention is to commit `Cargo.lock`. If it is not in the repo, the version check will fail. This assumption should be made explicit.

### 1.6 Missing: concurrency control

No `concurrency` key is defined. Rapid pushes to `optimus_prime` will spawn parallel workflow runs, potentially pushing competing `:dev` tags to GHCR simultaneously. This can cause manifest corruption or wasted runner minutes. A concurrency group with `cancel-in-progress: true` should be added.

### 1.7 Docker tag on every push: overwriting semantics

Every push to `optimus_prime` pushes to `:dev`. This is the stated intent, but there is no mechanism to distinguish which commit a `:dev` image corresponds to. Adding a `:dev-<short-sha>` tag alongside `:dev` would cost nothing and enable debugging when nf-core tests fail against a specific image.

### 1.8 Job graph: upload-binaries depends on create-tag-and-release

`upload-binaries` needs the release ID to attach assets. The plan lists them as separate jobs but does not specify the dependency. `upload-binaries` must `needs: [build-binaries, create-tag-and-release]` and receive the release ID as an output. This dependency should be explicit in the plan.

### 1.9 Artifact retention: "until release" is vague

The plan says retention is "1 day for dev, until release for releases." GitHub Actions artifact retention is specified as a number of days (max 90). "Until release" is not a valid retention policy. For release builds, the artifacts are attached to the GitHub Release (which has no retention limit), so the build artifacts themselves can use the same short retention. This should be clarified.

---

## 2. Assumptions

### 2.1 Stated assumptions -- validation

| Assumption | Validated? | Notes |
|-----------|-----------|-------|
| Rust 1.85 is available in `rust:1.85-bookworm` | Partially | Rust edition 2024 requires 1.85+. The Docker image tag `rust:1.85-bookworm` pins to 1.85.x. However, this means the Dockerfile will NOT pick up newer Rust versions automatically. If a dependency starts requiring 1.86+, the Dockerfile breaks silently until updated. Consider `rust:1-bookworm` for auto-updates or document the pin-and-update policy. |
| `GITHUB_TOKEN` has GHCR write permissions | Yes | For public repos with `packages: write` permission. The plan does not show the `permissions` block. Must include `packages: write` and `contents: write`. |
| Trusted publishing for crates.io is configured | Assumed | Cannot verify from plan alone. If not yet configured on crates.io, the `publish-crate` job will fail. The plan should include a pre-implementation verification step. |
| `procps` is sufficient for Nextflow | Yes | Standard Nextflow requirement for resource monitoring. |
| No C dependencies in the Rust crate | Yes | Verified: flate2 uses `zlib-rs` (pure Rust), gzp uses `deflate_rust`. No system library linkage. |

### 2.2 Implicit assumptions surfaced

1. **`lto = true` in release profile**: The Cargo.toml has `lto = true` and `codegen-units = 1`. On GitHub Actions runners with limited memory, LTO for a full build can be slow (~5-10 min). This is acceptable but not acknowledged in the plan's time estimates (there are none).

2. **Binary name vs. crate name mismatch**: The crate is `trim-galore` (hyphen) but the binary is `trim_galore` (underscore). This is correctly handled by the `[[bin]]` section in Cargo.toml, but the plan's packaging step references `trim_galore` as the binary name. This is correct but worth calling out since it could confuse contributors.

3. **No smoke test of built binaries**: The plan builds binaries for 4 targets but never runs them (except implicitly via Docker). The linux binaries could be smoke-tested with `./trim_galore --help` after build. macOS binaries cannot be tested on Linux runners, but linux binaries can. The existing release.yml has a `test-install` job; the new plan drops this.

4. **No `--locked` flag on `cargo build`**: Without `--locked`, a build could use different dependency versions than what was tested in CI. The plan should use `cargo build --release --locked` everywhere.

5. **GHCR image path case sensitivity**: The plan uses `ghcr.io/felixkrueger/trimgalore`. GHCR requires lowercase paths. The GitHub username is `FelixKrueger` (mixed case). The plan correctly lowercases this, but it should be noted that `github.repository_owner` in Actions outputs `FelixKrueger` and must be lowercased explicitly (e.g., via shell `tr` or the `docker/metadata-action`).

6. **Edition 2024 requirement**: The `Cargo.toml` uses `edition = "2024"` which requires Rust 1.85+. The Dockerfile pins `rust:1.85-bookworm`, which is correct. However, the `dtolnay/rust-toolchain@stable` used in the binary build jobs will use whatever the current stable is. If a future Rust stable introduces breaking changes to edition 2024 semantics (unlikely but possible), builds could break. This is standard practice and acceptable, but the plan should be aware of the difference: Docker builds pin Rust 1.85, CI binary builds use latest stable.

---

## 3. Efficiency Analysis

### 3.1 Docker layer caching

The plan mentions "Uses GHA cache for layer caching" but the Dockerfile's `COPY . .` before `cargo build --release` means ANY file change (including README edits) invalidates the cargo build cache layer. This is a significant efficiency concern since the Rust compilation with LTO will take several minutes.

**Recommendation**: Use a two-stage COPY pattern:
```dockerfile
COPY Cargo.toml Cargo.lock ./
RUN cargo fetch
COPY src/ src/
RUN cargo build --release
```
This ensures dependency fetching is cached when only source files change.

### 3.2 Redundant builds: binary job + Docker job both compile

On every push, both `build-binaries` (linux-x86_64 + linux-aarch64) and `docker-build` (amd64 + arm64) compile the same binary independently. That is 4 Rust compilations for 2 unique binaries. The Docker build could instead COPY the pre-built binary from the `build-binaries` artifacts, skipping the in-container compilation entirely. This would:
- Halve the total build time for Linux targets
- Reduce runner minute consumption
- Ensure binary parity between the tarball and the container

Trade-off: slightly more complex workflow (Docker build depends on binary build). But the time savings are substantial given LTO builds.

### 3.3 Build on every push to optimus_prime is aggressive

Building 4 binary targets + 2 Docker images on every push to a development branch could consume significant Actions minutes. Consider:
- Only building binaries on release (workflow_dispatch) and PR merge
- Only building Docker `:dev` image on push (which is the stated need)
- Using path filters to skip builds when only docs/plans change

### 3.4 SHA256 checksums mentioned but not detailed

The plan says binaries are packaged as `.tar.gz + .sha256` but does not show how the SHA256 file is generated. This is a minor implementation detail but should be specified (e.g., `sha256sum trim_galore-*.tar.gz > trim_galore-*.tar.gz.sha256`).

---

## 4. Validation Sufficiency

### 4.1 No functional test of the Docker image

The plan says to "Pull container, run `trim_galore --help` to verify" as step 4 of Implementation Order, but this is a manual step, not part of the CI workflow. The workflow should include an automated smoke test job that:
- Pulls the just-pushed `:dev` image
- Runs `docker run ghcr.io/felixkrueger/trimgalore:dev --help`
- Optionally runs a minimal trim on a small test file

Without this, a Docker image could be pushed that fails at runtime (e.g., missing shared library, wrong ENTRYPOINT) and nf-core tests would break with an opaque error.

### 4.2 No validation that the tag does not already exist on GHCR

The `check-release` job verifies the git tag does not exist, but does not check whether the Docker tag (e.g., `v2.0.0`) already exists on GHCR. If a release is re-triggered after a partial failure, the Docker push might overwrite or conflict with an existing image.

### 4.3 No binary architecture validation

There is no step to verify that the built binary is actually the correct architecture. On macOS cross-compilation (see issue 1.3), it is easy to accidentally produce an ARM binary when targeting x86_64. A `file trim_galore` check on Linux or `lipo -info trim_galore` on macOS would catch this.

### 4.4 crates.io publish uses --no-verify

The plan specifies `cargo publish --no-verify`, which skips the pre-publish build check. This means a broken crate could be published if the build succeeded but the package metadata is wrong (e.g., missing files in the include list). Since the binary build already succeeded, the risk is low, but it should be a conscious decision documented in the plan.

---

## 5. Alternatives Considered

### 5.1 Use pre-built binaries in Docker (recommended)

As noted in 3.2, rather than compiling inside Docker, copy the binary from the build-binaries job. This is how many Rust projects handle it (e.g., ripgrep, bat). The Dockerfile becomes:

```dockerfile
FROM debian:bookworm-slim
COPY trim_galore /usr/local/bin/trim_galore
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
ENTRYPOINT ["trim_galore"]
```

The multi-stage build Dockerfile is still useful for local development (`docker build .`), so keep it but have CI use the simpler approach.

### 5.2 Use `cargo-zigbuild` for cross-compilation

Instead of requiring native runners for each architecture, `cargo-zigbuild` can cross-compile to any glibc target from a single runner. This would simplify the matrix to a single `ubuntu-latest` runner for all Linux targets and a single `macos-latest` for all macOS targets. Trade-off: adds a build dependency but reduces runner complexity.

### 5.3 Consider `workflow_call` for reusable release logic

If CI and release share compilation steps, a reusable workflow (`workflow_call`) could avoid duplication between `ci.yml` and `release.yml`. Not critical now but worth noting for maintainability.

### 5.4 Nightly `:dev` builds instead of every-push

Instead of building on every push, a nightly schedule (`cron: '0 4 * * *'`) building from `optimus_prime` HEAD would reduce runner usage while still catching issues before release. The trade-off is slower feedback, but for a `:dev` container used by nf-core tests (which are not blocking PRs), nightly is often sufficient.

---

## 6. Action Items

### Critical

| # | Issue | Section | Action |
|---|-------|---------|--------|
| C1 | macOS x86_64 builds on ARM runner | 1.3 | Fix the build matrix: either use `macos-13` for x86_64 or explicitly handle cross-compilation with `rustup target add`. Decide which approach before implementation. |
| C2 | Trigger YAML contradicts Decision #2 | 1.1 | Add `master` to the push branches trigger, or update Decision #2 to match the actual trigger. |
| C3 | workflow_dispatch lacks branch guard | 1.2 | Add a check in `check-release` that fails if `github.ref` is not `refs/heads/master`. |
| C4 | Missing `permissions` block | 2.1 | The workflow must declare `packages: write` for GHCR push and `contents: write` for creating releases/tags. Without this, GHCR pushes will fail with 403. |

### Important

| # | Issue | Section | Action |
|---|-------|---------|--------|
| I1 | No concurrency control | 1.6 | Add `concurrency: { group: release-${{ github.ref }}, cancel-in-progress: true }` at the workflow level. |
| I2 | Docker layer cache invalidated by any file change | 3.1 | Restructure Dockerfile to COPY Cargo.toml/Cargo.lock first, then `cargo fetch`, then COPY src. |
| I3 | No automated Docker smoke test | 4.1 | Add a job after `docker-merge` that pulls the image and runs `--help`. |
| I4 | Missing `--locked` flag | 2.2 | Use `cargo build --release --locked` in all build steps to ensure reproducible builds. |
| I5 | upload-binaries dependency on release not specified | 1.8 | Make explicit that `upload-binaries` needs `create-tag-and-release` and receives the release ID. |
| I6 | GHCR path case sensitivity | 2.2 | Ensure `github.repository_owner` is lowercased before use in image paths. |
| I7 | No smoke test of built linux binaries | 2.2 | Add `./trim_galore --help` after building linux targets (cheap sanity check). |
| I8 | Consider `:dev-<sha>` tag | 1.7 | Push both `:dev` and `:dev-<short-sha>` tags for debuggability. |

### Optional

| # | Issue | Section | Action |
|---|-------|---------|--------|
| O1 | Redundant compilation in Docker vs. binary jobs | 3.2, 5.1 | Consider using pre-built binaries in Docker image instead of recompiling. Saves ~50% of Linux build time. |
| O2 | Build on every push is aggressive | 3.3 | Consider path filters to skip builds for docs-only changes, or limit binary builds to releases only. |
| O3 | Rust version pin in Dockerfile | 2.1 | Document the policy for updating the Rust version pin in the Dockerfile. Consider `rust:1-bookworm` for auto-updates. |
| O4 | `cargo publish --no-verify` risk | 4.4 | Document why `--no-verify` is acceptable (binary already built) as a comment in the workflow. |
| O5 | Architecture validation | 4.3 | Add `file` command check on built binaries to verify correct architecture. |

---

## Summary

The plan is well-structured and covers the key requirements (dev containers for nf-core, release binaries, crates.io publishing). The core architecture -- dual-trigger workflow with native builds and multi-arch Docker -- is sound.

The four critical issues must be resolved before implementation: the macOS x86_64 cross-compilation gap (will produce wrong-arch binaries or fail), the trigger/decision contradiction (master branch won't get dev containers), the unguarded workflow_dispatch (releases from any branch), and the missing permissions block (GHCR pushes will 403).

The most impactful efficiency improvement is avoiding redundant Rust compilation by reusing binary artifacts in Docker builds (O1/5.1), though this adds workflow complexity and can be deferred to a follow-up iteration.
