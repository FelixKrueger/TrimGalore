# Plan: Release CI — Automated Binaries & Containers

## Goal

Build Docker containers (multi-arch) on every push to `optimus_prime` so that
nf-core pipeline tests (nf-test with `--profile test` / `test_full`) can use
`ghcr.io/felixkrueger/trimgalore:dev` to catch issues before a proper release.
On release, also produce versioned containers, binary tarballs, and publish
to crates.io.

Modelled after [RustQC's release.yml](https://github.com/seqeralabs/RustQC/blob/main/.github/workflows/release.yml),
simplified for TrimGalore's needs (no SIMD variants, no C build deps).

---

## Deliverables

| # | File | Purpose |
|---|------|---------|
| 1 | `Dockerfile` | Multi-stage build: compile in Rust image, copy binary to minimal runtime |
| 2 | `.github/workflows/release.yml` | Rewritten: dual-trigger (push → dev, workflow_dispatch → release) |

---

## Detailed Design

### 1. Dockerfile

Multi-stage with dependency caching layer:

```dockerfile
# ── Build stage ──────────────────────────────────────────────
FROM rust:1.85-bookworm AS builder

WORKDIR /build

# Cache dependency build: copy manifests first, build deps, then copy source
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs \
    && cargo build --release --locked \
    && rm -rf src

# Now copy real source and rebuild (only the binary, deps are cached)
COPY . .
RUN cargo build --release --locked

# ── Runtime stage ────────────────────────────────────────────
FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/FelixKrueger/TrimGalore"
LABEL org.opencontainers.image.description="Trim Galore — adapter/quality trimming (Oxidized Edition)"
LABEL org.opencontainers.image.licenses="GPL-3.0-only"

RUN apt-get update && apt-get install -y --no-install-recommends \
      procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/trim_galore /usr/local/bin/trim_galore

ENTRYPOINT ["trim_galore"]
```

**Design decisions:**

- **Dependency caching layer** (review fix #5): `Cargo.toml` + `Cargo.lock` are
  copied first and deps are built with a dummy `main.rs`. The real source copy
  only invalidates the final compile step, not the dependency download/build.
- **`--locked`** (review fix #3): ensures the build uses exactly the dependency
  versions from `Cargo.lock`, matching what CI tested.
- `debian:bookworm-slim` rather than `scratch`/`distroless` — Nextflow needs
  `/bin/sh` and basic utilities (`ps` via `procps`) for process management.
- No `CMD` — Nextflow passes the full command.
- `procps` is the only extra package (Nextflow uses `ps` for resource monitoring).

### 2. Release Workflow (`.github/workflows/release.yml`)

#### Permissions (review fix #8)

```yaml
permissions:
  contents: write    # create tags + releases
  packages: write    # push to GHCR
```

#### Trigger Strategy (review fix #2, #7)

```yaml
on:
  push:
    branches: [optimus_prime, master]   # → dev containers + binary artifacts
  workflow_dispatch: {}                  # → full release (tag, binaries, containers, crates.io)

concurrency:
  group: release-${{ github.ref }}
  cancel-in-progress: true
```

- **Every push to `optimus_prime` or `master`**: builds Docker image + pushes
  `ghcr.io/felixkrueger/trimgalore:dev`
- **Manual `workflow_dispatch`**: reads version from `Cargo.toml`, validates the
  tag doesn't exist yet, then creates a full release with versioned containers
  + binaries.
- **`workflow_dispatch` branch guard** (review fix #7): Job 0 (`check-release`)
  validates that `workflow_dispatch` was triggered on `master` branch only. If
  triggered from any other branch, the job fails with a clear error message.
- **Concurrency group** (review fix #6): rapid pushes to the same branch cancel
  in-progress runs, preventing race conditions on the `:dev` tag.

#### Job Graph

```
check-release              ← detects release vs dev, validates branch
    │
    ├── build-binaries     ← 4 targets (only on release)
    │
    ├── docker-build       ← 2 platforms (amd64, arm64) on native runners
    │       │
    │       └── docker-merge       ← multi-arch manifest
    │               │
    │               └── smoke-test ← pull image, run trim_galore --help
    │
    └── [release only]:
            ├── create-tag-and-release  ← git tag + GitHub Release
            ├── upload-binaries         ← attach .tar.gz + .sha256 to release
            └── publish-crate           ← cargo publish to crates.io
```

#### Job Details

**Job 0: `check-release`**
- If `workflow_dispatch`:
  - **Validate branch is `master`** (review fix #7) — fail with error if not.
  - Extract version from `Cargo.toml`, verify tag doesn't exist, verify
    `Cargo.lock` version matches.
  - Output `is_release=true`, `version=v2.0.0`, etc.
- If push: output `is_release=false`.

**Job 1: `build-binaries`** (release only) (review fix #10)
- Only runs on release (`is_release == 'true'`). Dev pushes only need Docker
  images for nf-core testing — binary tarballs are not needed until release.
- All builds use `cargo build --release --locked` (review fix #3).
- Matrix: 4 targets

| Name | OS | Target | Notes |
|------|----|--------|-------|
| linux-x86_64 | ubuntu-latest | x86_64-unknown-linux-gnu | native |
| linux-aarch64 | ubuntu-24.04-arm | aarch64-unknown-linux-gnu | native |
| macos-x86_64 | **macos-13** | x86_64-apple-darwin | Intel runner (review fix #1) |
| macos-aarch64 | macos-latest | aarch64-apple-darwin | native ARM64 |

- **macOS x86_64 uses `macos-13`** (review fix #1): `macos-latest` resolves to
  ARM64 runners. Building `x86_64-apple-darwin` there would be cross-compilation.
  `macos-13` is the last Intel runner generation.
- gnu targets (not musl) — Linux arches build on native runners. See Resolved
  Decisions §1.
- Package: `trim_galore-<name>.tar.gz` + `.sha256`
- Upload as build artifacts

**Job 1a: `smoke-test-binaries`** (release only, after build-binaries) (review fix #4)
- Downloads linux-x86_64 artifact
- Extracts and runs `trim_galore --help`
- Verifies exit code 0 and output contains expected version string

**Job 2: `docker-build`** (runs always, on native runners)
- Matrix: `linux/amd64` on `ubuntu-latest`, `linux/arm64` on `ubuntu-24.04-arm`
- Uses `docker/build-push-action` with push-by-digest
- Exports digest as artifact for merge step
- Uses GHA cache for layer caching

**Job 3: `docker-merge`** (runs always, after docker-build)
- Downloads per-platform digests
- Creates multi-arch manifest:
  - Always: `dev` tag
  - On release: `v2.0.0`, `2.0`, `2`, `latest`
- Pushes to `ghcr.io/felixkrueger/trimgalore`

**Job 3a: `smoke-test-docker`** (runs always, after docker-merge) (review fix #4)
- Pulls `ghcr.io/felixkrueger/trimgalore:dev`
- Runs `docker run --rm ghcr.io/felixkrueger/trimgalore:dev --help`
- Verifies exit code 0 and output contains expected strings

**Job 4: `create-tag-and-release`** (release only, after all builds + smoke tests)
- Creates git tag from version
- Extracts changelog section (if CHANGELOG.md exists) or uses auto-generated notes
- Creates GitHub Release (draft: false)

**Job 5: `upload-binaries`** (release only, after create-tag-and-release)
- Downloads binary artifacts, attaches to release

**Job 6: `publish-crate`** (release only, after upload-binaries)
- `cargo publish --no-verify`
- Uses `rust-lang/crates-io-auth-action` for OIDC token (no stored API key needed)
- Trusted publishing is configured ✓

---

## Resolved Decisions

### 1. musl vs gnu for Linux binaries → **gnu**

Switching from musl to gnu. Each architecture builds on its own native runner
(no cross-compilation). The primary Linux deployment is containers (which have
glibc). Bare-metal users are on modern distros (glibc 2.31+, Ubuntu 20.04+).

### 2. Branch trigger → **both optimus_prime and master**

- Push to `optimus_prime` → `:dev` container
- Push to `master` → `:dev` container (in case of direct pushes/merges)
- `workflow_dispatch` on `master` only → full release (versioned tags, binaries, GitHub Release)

### 3. crates.io publishing → **included**

Trusted publishing is configured on crates.io. The workflow will use
`rust-lang/crates-io-auth-action` for OIDC token — no stored API keys needed.
Publish job runs as the final step of a release (after binaries + containers).

### 4. Container registry → **GHCR**

`ghcr.io/felixkrueger/trimgalore` — free for public repos, uses `GITHUB_TOKEN`,
no extra secrets needed. BioContainers/bioconda handled separately by the
bioconda recipe.

### 5. SIMD variants → **none**

TrimGalore is I/O-bound (gzip compression), not compute-bound. No measurable
benefit from AVX2/AVX-512 specific builds unlike RustQC.

---

## Review Fixes Applied

| # | Finding | Fix |
|---|---------|-----|
| 1 | macOS x86_64 runner mismatch | Changed to `macos-13` (Intel runner) |
| 2 | Trigger YAML missing `master` | Added `master` to push branches |
| 3 | Missing `--locked` flag | Added to all `cargo build` commands |
| 4 | No smoke tests | Added `smoke-test-binaries` + `smoke-test-docker` jobs |
| 5 | Docker cache invalidation | Two-step COPY pattern in Dockerfile |
| 6 | No concurrency control | Added `concurrency` group with `cancel-in-progress` |
| 7 | `workflow_dispatch` no branch guard | `check-release` validates branch is `master` |
| 8 | Missing `permissions` block | Added `contents: write` + `packages: write` |
| 9 | No CI gating for dev builds | Addressed by concurrency group — broken builds get cancelled by the next push. Full CI gating would require a complex cross-workflow dependency; the smoke test after Docker merge catches broken images before they're used. |
| 10 | Wasteful binary builds on dev push | Binary builds now release-only; dev pushes only build Docker images |

---

## How to Use the Dev Container in nf-core Tests

Once the workflow is live, every push to `optimus_prime` produces:

```
ghcr.io/felixkrueger/trimgalore:dev
```

To test with nf-core/rnaseq:

```bash
# In an nf-test or manual run:
nextflow run nf-core/rnaseq \
  -profile test,docker \
  --outdir results \
  -c custom.config
```

Where `custom.config` overrides the trimgalore container:

```nextflow
process {
    withName: 'TRIMGALORE' {
        container = 'ghcr.io/felixkrueger/trimgalore:dev'
    }
}
```

Or for nf-test, add it to the test's `nextflow.config`.

---

## Implementation Order

1. Create `Dockerfile` (test locally with `docker build`)
2. Rewrite `.github/workflows/release.yml`
3. Push to `optimus_prime` → verify dev container appears at GHCR
4. Pull container, run `trim_galore --help` to verify
5. Test with nf-core/rnaseq `--profile test`

---

## What We're NOT Doing (and why)

- **No SIMD variants**: TrimGalore is I/O-bound (gzip), not compute-bound. Unlike
  RustQC which benefits from AVX2/AVX-512 for quality score processing, there's no
  measurable gain from CPU-specific builds.
- **No separate CI workflow changes**: The existing `ci.yml` handles testing on push.
  The release workflow focuses purely on build + publish.
- **No Singularity-specific build**: Singularity can pull Docker images from GHCR
  directly (`singularity pull docker://ghcr.io/...`). No need for a separate `.sif`.
- **No cross-workflow CI gating**: Making `release.yml` depend on `ci.yml` success
  requires complex workflow_run triggers or status checks. The smoke tests after
  Docker merge provide a safety net; truly broken code is caught by CI separately.
