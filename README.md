# Trim Galore

Consistent quality and adapter trimming for next-generation sequencing data, with special handling for RRBS and bisulfite-seq libraries.

[![CI](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml/badge.svg?branch=optimus_prime)](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/trim-galore)](https://crates.io/crates/trim-galore)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/trim-galore/README.html)

> [!NOTE]
> **Trim Galore v2.0** is a complete rewrite in Rust — a single binary with zero external dependencies, producing byte-identical output to v0.6.x. Same CLI, same output filenames, same report format. For details on what changed, benchmarks, and migration notes, see the [v2.0 writeup](docs/SUMMARY.md).

## Features

- **Adapter auto-detection** — automatically identifies Illumina, Nextera, Small RNA, Stranded, and BGI adapters
- **Quality trimming** — Phred-based trimming from the 3' end (BWA algorithm)
- **Paired-end** — single-pass processing of both reads with automatic pair validation
- **RRBS / bisulfite-seq** — MspI end-repair artifact removal, directional and non-directional libraries
- **NextSeq / 2-colour trimming** — poly-G removal for NovaSeq, NextSeq, and NovaSeq X
- **Poly-A / poly-G trimming** — built-in tail trimming without external tools
- **Parallel compression** — `--cores N` for faster gzip I/O on multi-core systems
- **FastQC integration** — optional post-trimming quality reports
- **MultiQC compatible** — trimming reports parse cleanly in MultiQC dashboards
- **Demultiplexing** — 3' inline barcode demultiplexing

## Installation

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

### Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon) will be available on the [Releases](https://github.com/FelixKrueger/TrimGalore/releases) page.

## Usage

```bash
# Single-end
trim_galore input.fastq.gz

# Paired-end
trim_galore --paired file_R1.fastq.gz file_R2.fastq.gz

# With parallel compression (recommended for large files)
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

| Mode | Trimmed output | Report |
|------|---------------|--------|
| Single-end | `*_trimmed.fq.gz` | `*_trimming_report.txt` |
| Paired-end | `*_val_1.fq.gz` / `*_val_2.fq.gz` | per-read report |
| Unpaired (with `--retain_unpaired`) | `*_unpaired_1.fq.gz` / `*_unpaired_2.fq.gz` | |

## Documentation

- [User Guide](docs/Trim_Galore_User_Guide.md) — full reference for all options and modes
- [RRBS Guide](docs/RRBS_Guide.pdf) — bisulfite sequencing and RRBS-specific guidance
- [v2.0 Writeup](docs/SUMMARY.md) — benchmarks, architecture, and what changed in v2.0

## Credits

_Trim Galore_ was developed at The Babraham Institute by [@FelixKrueger](https://github.com/FelixKrueger/), now part of [Altos Labs](https://altoslabs.com/).

## License

[GPL-3.0](LICENSE)
