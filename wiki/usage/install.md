---
type: usage
title: Getting Trim Galore
description: The distribution channels, which one to use in a pipeline, and what the binary needs at runtime.
tags: [install, bioconda, docker, crates]
page_id: install
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: Cargo.toml
    commit: 2ab0918
  - repo: TrimGalore
    path: Dockerfile
    commit: 8de2bb2
  - repo: TrimGalore
    path: src/fastqc.rs
    commit: 242546a
---

## Summary

Trim Galore ships as one binary called `trim_galore`, with no external runtime
dependencies. The developer recommends bioconda for pipeline users. Every
channel is public and needs no access request.

## Channels

| Channel | Command | Notes |
|---|---|---|
| bioconda | `conda install -c bioconda trim-galore` | **Recommended for pipeline users** (developer's recommendation). |
| crates.io | `cargo install trim-galore` | Builds locally; needs the Rust toolchain at or above the crate's `rust-version` floor. |
| Homebrew | `brew install trim-galore` | macOS and Linux, from homebrew-core. |
| Docker | `ghcr.io/felixkrueger/trimgalore:latest` | Multi-arch, amd64 and arm64. |
| Prebuilt binaries | GitHub Releases | Linux x86_64 and aarch64, macOS Apple Silicon. |
| Source | `cargo build --release` | Binary lands at `target/release/trim_galore`. |

The development branch can be installed directly with
`cargo install --git … --branch dev trim-galore --force`; `--force` is needed
to replace an existing `trim_galore` from a release install.

Intel Macs have no prebuilt binary — build locally through crates.io, or use
the amd64 Docker image.

## Nothing to install alongside it

Adapter detection, alignment, quality trimming, filtering, gzip compression
*and* FastQC reporting all run inside the process. There is no Java, no
Python, no external `cutadapt`, `pigz` or `fastqc` to provide. FastQC comes
from the bundled `fastqc-rust` library, wired up in `fastqc.rs`, so `--fastqc`
works in a container that has nothing else in it.

That is what makes the single-binary channels viable on HPC and in minimal
images: the binary plus its inputs is the whole requirement. CI proves the
FastQC half of it by asserting that no external `fastqc` is on `PATH` during
the run that produces a report.

## Pinned dependency worth knowing about

`fastqc-rust` is pinned to an exact version in `Cargo.toml` rather than a
range. The pin is deliberate, so treat any bump as a deliberate test event
and re-verify report output rather than a routine dependency update.

## Links

- `usage/trim-reads.md` — the first run once it is installed.
- `reference/cli-flags.md` — the flag surface.
- `overview.md` — whether this is the right tool for the job.

## Open questions

- None.
