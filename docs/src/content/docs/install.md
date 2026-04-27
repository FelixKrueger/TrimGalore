---
title: Installation
description: Install Trim Galore via cargo, bioconda, Docker, or prebuilt binaries.
---

Trim Galore v2 ships as a single static binary. No Python, no Cutadapt, no `pigz`, no Java, no external FastQC. Pick whichever channel fits the rest of your stack.

:::caution[Beta testing v2.1.0-beta.5]
The current stable release on crates.io is **v2.0.0**. v2.1.0 is in beta. `cargo install trim-galore` without `--version` installs v2.0.0. To install the beta:

```bash
cargo install trim-galore --version 2.1.0-beta.5
docker pull ghcr.io/felixkrueger/trimgalore:beta
```

To send feedback on the beta, open an issue with the `beta-feedback` label. This notice goes away at v2.1.0 GA.
:::

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
    ghcr.io/felixkrueger/trimgalore \
    trim_galore input.fastq.gz
```

FastQC is built in via the bundled [`fastqc-rust`](https://crates.io/crates/fastqc-rust) library, so `--fastqc` works without Java or an external `fastqc` install. The `dev` tag tracks the development branch. Versioned tags (e.g. `v2.1.0`) are published on release.

## Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon) are on the [Releases page](https://github.com/FelixKrueger/TrimGalore/releases). On Intel Mac, install via `cargo install trim-galore` (local build) or use the Docker `amd64` image.

## Runtime dependencies

None. Trim Galore v2 is a single static binary. Adapter trimming, gzip, and FastQC reporting (`--fastqc`) all run in-process. Python, Cutadapt, `pigz`, Java, and the external FastQC tarball are no longer required.

## Verifying the install

```bash
trim_galore --version
```

Should print something like:

```
                Trim Galore version: 2.1.0-beta.5 (Oxidized Edition)
                Cutadapt version:    (built-in adapter trimmer)
                Python version:      (not required)
                FastQC version:      0.12.1-compatible (bundled fastqc-rust)
```
