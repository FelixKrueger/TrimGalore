<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/FelixKrueger/TrimGalore/master/docs/public/logos/hero-dark.svg">
    <img alt="Trim Galore" src="https://raw.githubusercontent.com/FelixKrueger/TrimGalore/master/docs/public/logos/hero-light.svg" width="800">
  </picture>
</h1>

Consistent quality and adapter trimming for next-generation sequencing data, with special handling for RRBS libraries.


[![CI](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/FelixKrueger/TrimGalore/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/trim-galore)](https://crates.io/crates/trim-galore)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/trim-galore/README.html)

<h3 align="center"><a href="https://www.trimgalore.com/">https://www.trimgalore.com/</a></h3>

> [!NOTE]
> **Trim Galore v2.0** is a faithful Rust rewrite — a single binary with zero external dependencies, designed as a drop-in replacement for v0.6.x scripts and pipelines. Same CLI, same output filenames, same report format. Adds poly-G auto-detection and trimming for 2-colour instruments, a generic poly-A trimmer, per-pair adapter auto-detection, and cleaner multi-adapter invocation (repeatable `-a`/`-a2` instead of Perl's embedded-string syntax) — among other extensions. For details on what changed, benchmarks, and migration notes, see the [v2.0 migration notes](https://www.trimgalore.com/reference/migration/).

## Features

- **Adapter auto-detection** — automatically identifies Illumina, Nextera, Small RNA, and BGI/DNBSEQ adapters from the first 1M reads. Stranded Illumina remains explicit (`--stranded_illumina`) because its sequence is ambiguous with Nextera.
- **Multi-adapter support** — specify multiple adapters by repeating `-a`/`-a2` or via `-a "file:adapters.fa"`, with optional multi-round trimming (`-n`)
- **Quality trimming** — Phred-based trimming from the 3' end (BWA algorithm)
- **Paired-end** — single-pass processing of both reads with automatic pair validation
- **RRBS** — MspI end-repair artifact removal, directional and non-directional libraries
- **Poly-G trimming** — sequence-based removal of no-signal G-runs at the 3' end of Read 1 (and poly-C at the 5' end of Read 2) from 2-colour instruments (NovaSeq, NextSeq, NovaSeq X). Auto-detected from the data; opt-out with `--no_poly_g`
- **NextSeq / 2-colour quality trim** — `--nextseq N` / `--2colour N` applies 2-colour-aware quality trimming (opt-in; replaces `-q`)
- **Poly-A trimming** — built-in removal of poly-A tails without external tools; recommended for mRNA-seq / poly-A-selected RNA-seq libraries
- **Parallel processing** — `--cores N` runs trimming and gzip compression in worker threads for near-linear speedup on multi-core systems
- **FastQC integration** — optional post-trimming quality reports built in via the bundled [fastqc-rust](https://crates.io/crates/fastqc-rust) library; produces FastQC 0.12.1-compatible HTML + ZIP outputs without requiring Java or an external `fastqc` on `$PATH`
- **MultiQC compatible** — trimming reports parse cleanly in MultiQC dashboards (text + JSON)
- **Demultiplexing** — 3' inline barcode demultiplexing

## Installation

### From crates.io

Requires the [Rust toolchain](https://rustup.rs/) (1.88+):

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
cargo install --git https://github.com/FelixKrueger/TrimGalore --branch dev trim-galore --force
```

The `--force` flag overwrites any existing `trim_galore` binary (e.g. a v2.0.0 install from crates.io).

### Docker

Multi-arch images (amd64 + arm64) are available from GitHub Container Registry:

```bash
docker run --rm -v "$PWD":/data -w /data ghcr.io/felixkrueger/trimgalore:latest trim_galore input.fastq.gz
```

FastQC is built into the binary itself via the bundled fastqc-rust library — no external `fastqc` or Java runtime needed in the image. Tags published: `:latest` (latest stable, currently `v2.1.0`), `:v2.1.0` (pinned to a specific release), `:beta` (latest prerelease — only set during an active beta cycle), and `:dev` (every push to the `dev` branch). See the [docs site install page](https://www.trimgalore.com/install/) for the full table.

### Prebuilt binaries

Prebuilt binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon) are available on the [Releases](https://github.com/FelixKrueger/TrimGalore/releases) page. Intel Mac users: install via `cargo install trim-galore` (local build) or use the Docker `amd64` image.

## Usage

```bash
# Single-end
trim_galore input.fastq.gz

# Paired-end
trim_galore --paired file_R1.fastq.gz file_R2.fastq.gz

# Parallel processing (recommended for large files)
# Near-linear speedup up to ~8 cores on v2.1.0-beta.8; beyond that the
# gzip-output I/O on the storage layer typically becomes binding.
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

Output compression mirrors the input: gzipped input (`*.fastq.gz`) produces gzipped output (`*.fq.gz`); plain input (`*.fastq`) produces plain output (`*.fq`). Pass `--dont_gzip` to force plain output regardless. By default, gzip output is written at compression level 1 (fastest); pass `--compression <N>` (1–9) to override — decompressed content is byte-identical regardless of level, but level-1 `.fq.gz` files are roughly 75% larger than level-9 in exchange for substantially faster trimming.

Pass `--clumpify` to reorder reads inside each gzip member by canonical 16-mer minimizer so reads sharing similar sequence land adjacent on disk, letting gzip's 32 KB dictionary find longer back-references. The right configuration depends on what the trimmed FASTQ is used for:

```bash
# Pipeline intermediates (deleted after the run): reorder is essentially free
# and shrinks the file 15–35% on most data — net I/O win for the next step.
trim_galore --clumpify <input>

# Long-term storage / disk-constrained workdirs: add gzip L6 for 15–50% saving
# at 4–6× plain wall.
trim_galore --clumpify --compression 6 <input>

# Archival use, max compression
trim_galore --clumpify --compression 9 <input>

# Smaller output without the reorder cost (e.g. for 10x scRNA-seq, see docs)
trim_galore --compression 6 <input>
```

No information loss — only the on-disk order of records changes. Output records are byte-identical to the unsorted output and trimming reports are unaffected. `--clumpify` requires `--cores >= 2`.

Memory budget is controlled by the global `--memory` flag (default `1G`); bigger budgets give bigger per-gzip-member sort runs and better compression up to roughly the uncompressed input size, with sharply diminishing returns above ~2 GB. With enough memory, `--clumpify --compression 9` gets you within 1–2 percentage points of `bbmap clumpify` and `stevekm/squish` on the same data.

Intended for short reads (Illumina, AVITI). Long-read inputs (Oxford Nanopore, PacBio) and 10x scRNA-seq typically see no benefit or a small negative result; see [Clumpy compression](https://www.trimgalore.com/performance/clumpy/) for per-data-type recommendations.

The JSON report contains the same statistics as the text report in a structured format (schema v1), designed for native parsing by MultiQC.

## Documentation

Full documentation is published at [https://www.trimgalore.com/](https://www.trimgalore.com/)

- [User Guide](https://www.trimgalore.com/guide/overview/) — full reference for all options and modes
- [RRBS Guide](https://www.trimgalore.com/rrbs/guide/) — bisulfite sequencing and RRBS-specific guidance
- [v2.0 migration notes](https://www.trimgalore.com/reference/migration/) — what changed in the Rust rewrite
- [Benchmarks](https://www.trimgalore.com/performance/benchmarks/)

## Credits

_Trim Galore_ was developed at The Babraham Institute by [@FelixKrueger](https://github.com/FelixKrueger/), now part of [Altos Labs](https://altoslabs.com/).

## License

[GPL-3.0](LICENSE)
