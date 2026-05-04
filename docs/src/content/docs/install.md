---
title: Installation
description: Install Trim Galore via cargo, bioconda, Docker, or prebuilt binaries.
---

Trim Galore v2 ships as a single static binary. No Python, no Perl, no Cutadapt, no Java, no `igzip`, no `pigz`, no external FastQC. Pick whichever channel fits the rest of your stack.

## From crates.io

Requires the [Rust toolchain](https://rustup.rs/) (1.88+):

```bash
cargo install trim-galore
```

## From bioconda

```bash
conda install -c bioconda trim-galore
```

## Build from source

```bash
git clone https://github.com/FelixKrueger/TrimGalore.git
cd TrimGalore
cargo build --release
# Binary is at target/release/trim_galore
```

## Latest development version

To install the latest unreleased changes from the development branch:

```bash
cargo install --git https://github.com/FelixKrueger/TrimGalore --branch optimus_prime trim-galore --force
```

`--force` overwrites any existing `trim_galore` binary (e.g. a v2.0.0 install from crates.io).

## Docker

Multi-arch images (`amd64` and `arm64`) are published to the GitHub Container Registry:

```bash
docker run --rm -v "$PWD":/data -w /data \
    ghcr.io/felixkrueger/trimgalore:latest \
    trim_galore input.fastq.gz
```

FastQC is built in via the bundled [`fastqc-rust`](https://crates.io/crates/fastqc-rust) library, so `--fastqc` works without Java or an external `fastqc` install.

| Tag | Updates |
| --- | --- |
| `:latest` | latest stable release (currently v2.1.0) |
| `:v2.1.0` | pinned to a specific release |
| `:beta` | latest prerelease — only set during an active beta cycle |
| `:dev` | every push to the `optimus_prime` development branch |

## Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon) are on the [Releases page](https://github.com/FelixKrueger/TrimGalore/releases). On Intel Mac, install via `cargo install trim-galore` (local build) or use the Docker `amd64` image.

## Runtime dependencies

None. Trim Galore v2 is a single static binary. Adapter trimming, gzip, and FastQC reporting (`--fastqc`) all run in-process. Python, Perl, Cutadapt, Java, `igzip`, `pigz`, and the external FastQC tarball are no longer required.

## Verifying the install

```bash
trim_galore --version
```

Should print:

```
trim_galore 2.1.0 (Oxidized Edition)
<git-hash> — <os>/<arch> — built <ISO-8601-UTC>
```

The first line is also what `-V` prints. The second line is the build provenance (commit, target triple, deterministic build timestamp).
