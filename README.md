# Trim Galore

Consistent quality and adapter trimming for next-generation sequencing data, with special handling for RRBS libraries.

[![CI](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml/badge.svg?branch=optimus_prime)](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/trim-galore)](https://crates.io/crates/trim-galore)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/trim-galore/README.html)

> [!NOTE]
> **Trim Galore v2.0** is a complete rewrite in Rust — a single binary with zero external dependencies, producing byte-identical output to v0.6.x. Same CLI, same output filenames, same report format. For details on what changed, benchmarks, and migration notes, see the [v2.0 writeup](docs/SUMMARY.md).

## Features

- **Adapter auto-detection** — automatically identifies Illumina, Nextera, and Small RNA adapters from the first 1M reads. Stranded Illumina and BGI/DNBSEQ adapters are selectable via explicit flags (`--stranded_illumina`, `--bgiseq`)
- **Multi-adapter support** — specify multiple adapters via `-a " SEQ1 -a SEQ2"` or `-a "file:adapters.fa"`, with optional multi-round trimming (`-n`)
- **Quality trimming** — Phred-based trimming from the 3' end (BWA algorithm)
- **Paired-end** — single-pass processing of both reads with automatic pair validation
- **RRBS** — MspI end-repair artifact removal, directional and non-directional libraries
- **Poly-G trimming** — sequence-based removal of no-signal G-runs at the 3' end of Read 1 (and poly-C at the 5' end of Read 2) from 2-colour instruments (NovaSeq, NextSeq, NovaSeq X). Auto-detected from the data; opt-out with `--no_poly_g`
- **NextSeq / 2-colour quality trim** — `--nextseq N` / `--2colour N` applies 2-colour-aware quality trimming (opt-in; replaces `-q`)
- **Poly-A trimming** — built-in removal of poly-A tails without external tools; recommended for mRNA-seq / poly-A-selected RNA-seq libraries
- **Parallel processing** — `--cores N` runs trimming and gzip compression in worker threads for near-linear speedup on multi-core systems
- **FastQC integration** — optional post-trimming quality reports (FastQC v0.12.1 bundled in the Docker image; for `cargo`/source installs, requires FastQC on `$PATH`)
- **MultiQC compatible** — trimming reports parse cleanly in MultiQC dashboards (text + JSON)
- **Demultiplexing** — 3' inline barcode demultiplexing

## Installation

> [!IMPORTANT]
> **Beta testing v2.1.0-beta.1.** The current stable release on crates.io is v2.0.0; v2.1.0 is in beta testing. Running `cargo install trim-galore` without `--version` installs v2.0.0 (the stable). To install the beta:
>
> ```bash
> cargo install trim-galore --version 2.1.0-beta.1
> docker pull ghcr.io/felixkrueger/trimgalore:beta
> ```
>
> Feedback on the beta is welcome — open an issue with the `beta-feedback` label. This section will be removed at v2.1.0 GA.

### From crates.io

Requires the [Rust toolchain](https://rustup.rs/) (1.85+):

```bash
cargo install trim-galore
```

### From bioconda

```bash
conda install -c bioconda trim-galore
```

### Build from source

```bash
git clone https://github.com/FelixKrueger/TrimGalore.git
cd TrimGalore
cargo build --release
# Binary is at target/release/trim_galore
```

### Latest development version

To install the latest unreleased changes directly from the development branch:

```bash
cargo install --git https://github.com/FelixKrueger/TrimGalore --branch optimus_prime trim-galore --force
```

The `--force` flag overwrites any existing `trim_galore` binary (e.g. a v2.0.0 install from crates.io).

### Docker

Multi-arch images (amd64 + arm64) are available from GitHub Container Registry:

```bash
docker run --rm -v "$PWD":/data -w /data ghcr.io/felixkrueger/trimgalore trim_galore input.fastq.gz
```

The image bundles FastQC v0.12.1 and a Java runtime, so `--fastqc` works out of the box. The `dev` tag tracks the latest development branch; versioned tags (e.g. `v2.1.0`) are published on release.

### Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon) are available on the [Releases](https://github.com/FelixKrueger/TrimGalore/releases) page. Intel Mac users: install via `cargo install trim-galore` (local build) or use the Docker `amd64` image.

## Usage

```bash
# Single-end
trim_galore input.fastq.gz

# Paired-end
trim_galore --paired file_R1.fastq.gz file_R2.fastq.gz

# Parallel processing (recommended for large files)
# Near-linear speedup up to ~16 cores; diminishing returns beyond ~20.
trim_galore --cores 8 --paired file_R1.fastq.gz file_R2.fastq.gz

# RRBS mode
trim_galore --rrbs --paired file_R1.fastq.gz file_R2.fastq.gz

# Run FastQC on trimmed output
trim_galore --fastqc input.fastq.gz
```

For the complete list of options:

```bash
trim_galore --help
```

### Output files

| Mode | Trimmed output | Reports |
|------|---------------|---------|
| Single-end | `*_trimmed.fq.gz` | `*_trimming_report.txt` + `*_trimming_report.json` |
| Paired-end | `*_val_1.fq.gz` / `*_val_2.fq.gz` | per-read text + JSON reports |
| Unpaired (with `--retain_unpaired`) | `*_unpaired_1.fq.gz` / `*_unpaired_2.fq.gz` | |

The JSON report contains the same statistics as the text report in a structured format (schema v1), designed for native parsing by MultiQC.

## Documentation

- [User Guide](docs/Trim_Galore_User_Guide.md) — full reference for all options and modes
- [RRBS Guide](docs/RRBS_Guide.pdf) — bisulfite sequencing and RRBS-specific guidance
- [v2.0 Writeup](docs/SUMMARY.md) — benchmarks, architecture, and what changed in v2.0

## Credits

_Trim Galore_ was developed at The Babraham Institute by [@FelixKrueger](https://github.com/FelixKrueger/), now part of [Altos Labs](https://altoslabs.com/).

## License

[GPL-3.0](LICENSE)
