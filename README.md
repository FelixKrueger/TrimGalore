# Trim Galore - Oxidized Edition

A complete Rust rewrite of [Trim Galore](https://github.com/FelixKrueger/TrimGalore/tree/master) — the widely used adapter and quality trimmer for next-generation sequencing data. This is a **drop-in replacement** that produces **byte-identical output** to the Perl original, with no external dependencies.

[![CI](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml/badge.svg?branch=optimus_prime)](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml)

> **Status: Pre-release (beta testing).** The Oxidized Edition passes all validation tests and produces byte-identical output to Trim Galore v0.6.11 across all modes. It is currently undergoing real-world testing before replacing the Perl version on bioconda. Feedback and bug reports are welcome.

## What's different?

- **Zero dependencies.** No Python, no Cutadapt, no pigz — a single static binary.
- **Faster.** 1.9x faster wall time at 8 cores, up to 4.4x at 24 cores. Uses 2.3-5x less CPU time.
- **Same CLI.** All existing flags work. Same output filenames. Same report format (MultiQC-compatible).
- **Single-pass paired-end.** Both reads processed together — no temp files, guaranteed synchronization.
- **Built-in poly-G trimming.** Auto-detected for 2-colour instruments (NovaSeq, NextSeq).

For detailed benchmarks, see [docs/SUMMARY.md](docs/SUMMARY.md).

## Installation

### From crates.io (recommended)

Requires the [Rust toolchain](https://rustup.rs/) (1.74+):

```bash
cargo install trim-galore
```

### Build from source

```bash
git clone https://github.com/FelixKrueger/TrimGalore.git
cd TrimGalore
git checkout optimus_prime
cargo build --release
# Binary is at target/release/trim_galore
```

### Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon) will be available on the [Releases](https://github.com/FelixKrueger/TrimGalore/releases) page once beta testing completes.

## Usage

The CLI is identical to the Perl version:

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

## Validating output

To verify byte-identical output against the Perl version:

```bash
# Run both versions on the same input
trim_galore -o /tmp/oxidized test_files/illumina_10K.fastq.gz
legacy/trim_galore -o /tmp/perl test_files/illumina_10K.fastq.gz

# Compare decompressed output (compressed bytes differ, content is identical)
diff <(gzip -dc /tmp/perl/illumina_10K_trimmed.fq.gz) \
     <(gzip -dc /tmp/oxidized/illumina_10K_trimmed.fq.gz)
```

The `legacy/` directory contains the original Perl script for comparison testing. It requires Perl, Cutadapt, and optionally FastQC.

## Documentation

For general Trim Galore usage, see the [User Guide](Docs/Trim_Galore_User_Guide.md).

For benchmarks and technical details of the Oxidized Edition, see [docs/SUMMARY.md](docs/SUMMARY.md).

## Credits

_Trim Galore_ was developed at The Babraham Institute by [@FelixKrueger](https://github.com/FelixKrueger/), now part of [Altos Labs](https://altoslabs.com/).

## License

[GPL-3.0](LICENSE)
