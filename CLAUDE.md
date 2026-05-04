# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Trim Galore — Oxidized Edition — is a Rust rewrite of the original Perl Trim Galore (the v0.6.x script lives upstream and is no longer in this repo). It is a single-binary, single-pass adapter and quality trimmer for NGS FASTQ data. There are no external runtime dependencies: adapter detection, alignment, quality trimming, filtering, gzip compression, **and** FastQC reporting (via the bundled `fastqc-rust` library) all run in-process. No Java, no Python, no Cutadapt, no external `fastqc`.

`master` is the stable trunk (v2.1.0 GA and onwards); `dev` is the active development branch where prereleases are cut from. (Pre-v2.1.0 GA the layout was inverted — `optimus_prime` was the dev trunk and `master` held the legacy Perl release line; PR #214 collapsed that on 2026-05-04 by merging the Rust trunk into `master` and renaming `optimus_prime` to `dev`. The legacy Perl 0.6.11 source remains accessible via the `0.6.11` git tag.) The CI workflow `validation` job md5-compares Oxidized output against Perl Trim Galore 0.6.11 (installed from raw GitHub) for several core flag combinations — preserving byte-identity to v0.6.11 is a hard invariant for those flag paths.

## Build, lint, test

Requires Rust toolchain (the crate's `rust-version` floor is `1.88`; CI uses `dtolnay/rust-toolchain@stable`). Edition is 2024.

```bash
cargo build --release          # binary lands at target/release/trim_galore
cargo test                     # unit tests (cwd must be crate root — fixtures resolve relative to test_files/)
cargo test <name>              # run a single test by name substring
cargo fmt --all -- --check     # CI fails on unformatted code
cargo clippy --all-targets --release -- -D warnings   # CI fails on any clippy warning
```

To reproduce the `validation` CI job locally (md5-compare against Perl Trim Galore 0.6.11), see the steps in `.github/workflows/ci.yml` — it installs the Perl `trim_galore` from a pinned upstream URL via conda and diffs gzip-decompressed outputs.

### Reproducible builds

`build.rs` reads `SOURCE_DATE_EPOCH` to stamp the binary with a deterministic build timestamp; missing → wall-clock time, malformed → hard panic. The CI `reproducibility` job builds twice with the same `SOURCE_DATE_EPOCH` and asserts the resulting binaries are bit-identical. Don't introduce wall-clock time, hostnames, or absolute paths into the release binary.

`./target/release/trim_galore --version` (long form) must emit a provenance line `<git-hash> — <os>/<arch> — built <ISO-8601-UTC>` (literal em-dash); `-V` must NOT include that line. CI greps for both. The format is constructed in `build.rs` and printed via `env!("VERSION_BODY")`.

## Test fixtures

`test_files/` holds gzipped FASTQ fixtures used by both `cargo test` and the CI validation matrix: BS-seq paired-end (`BS-seq_10K_R{1,2}.fastq.gz`), RRBS (`SRR24766921_RRBS_R{1,2}.fastq.gz`), Clock-mode (`clock_10K_R{1,2}.fastq.gz`), demux (`demux_test.fastq.gz` + `demux_test_samplesheet.txt`), poly-A/T, plus negative cases (`colorspace_file.fastq`, `truncated.fq.gz`, `empty_file.fastq`). Tests in `src/cli.rs` reference these by relative path, so `cargo test` must be run from the crate root.

## Architecture

Single binary (`src/main.rs` is the only `[[bin]]`); the rest is a library crate (`src/lib.rs`) so unit tests can reach internals. `main()` parses CLI, runs a sanity check on the first input, then dispatches:

1. **Specialty modes** (run-and-exit, bypass the trimming pipeline): `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`. All four accept multi-pair input (an even number of files; per-pair "Pair N of M" headers + the same output-collision pre-flight that `--paired` runs).
2. **Paired mode** (`--paired`): consecutive files form R1/R2 pairs. Before any I/O, a pre-flight hashes prospective output paths case-folded (ASCII lowercase) so collisions on APFS/NTFS — and case-only aliases — fail loudly rather than silently overwriting. Adapter auto-detection runs **per pair** (intentional deviation from Perl v0.6.x, which detected once on `$ARGV[0]`).
3. **Single-end**: each input is processed independently in a loop.

### Module map (under `src/`)

- `cli.rs` — clap definitions, `rewrite_perl_short_flags()` pre-parser (fixes Perl-era `-r1`/`-r2`), `Cli::validate()` (paired-input sanity, duplicate-pair detection, mutually-exclusive flag conflicts beyond what clap expresses).
- `adapter.rs` — built-in adapter sequences (Illumina `AGATCGGAAGAGC`, Nextera `CTGTCTCTTATA`, smallRNA `TGGAATTCTCGG`, BGI/DNBSEQ, stranded Illumina) and auto-detection by scanning the first 1 M reads. BGI is probed; `--stranded_illumina` stays explicit-only because its sequence is ambiguous with Nextera.
- `alignment.rs` — semi-global unit-cost DP that re-implements Cutadapt's adapter matching algorithm. Byte-for-byte parity with Cutadapt is what makes the `validation` CI job pass.
- `quality.rs` — BWA-style 3' quality trimming (Li & Durbin 2009; matches Cutadapt's algorithm).
- `trimmer.rs` — orchestrator: wires quality trim → adapter trim → clip → filter → report for both single- and paired-end. Public `run_single_end` / `run_paired_end` entry points used by `main.rs`.
- `parallel.rs` — `--cores N` worker pool. Each worker trims **and** gzip-compresses its own chunk, and the chunks are concatenated in order — RFC 1952 permits gzip-member concatenation, so the output is a valid `.gz` file. This replaces the older readers→main→writers pipeline.
- `fastq.rs` — `FastqReader`/`FastqWriter` with gzip awareness and 64 KB buffered I/O; `FastqReader::sanity_check` is the entry-side guard against truncated/empty/colorspace input.
- `filters.rs` — length, max-length, max-N, and unpaired-rescue filters.
- `io.rs` — output naming: `*_trimmed.fq(.gz)` (SE), `*_val_{1,2}.fq(.gz)` (PE), `*_unpaired_{1,2}.fq(.gz)`, `*_trimming_report.txt` + `*_trimming_report.json`.
- `demux.rs` — 3' inline barcode demultiplexing (single-end only, matching the original).
- `specialty.rs` — `--hardtrim5/3`, `--clock` (Epigenetic Clock UMI), `--implicon` (UMI from R2). Each mode owns its output naming.
- `fastqc.rs` — bundled FastQC integration via the `fastqc-rust` crate (exact-pinned to `=1.0.1`; treat upstream bumps as deliberate test events). `--fastqc_args` accepts a curated subset of flags; unknown flags warn-and-ignore for forward compatibility.
- `report.rs` — generates the MultiQC-compatible text + JSON trimming reports.

### Adapter shorthand

`-a` / `-a2` accept `A{N}` shorthand (e.g. `-a A{10}` → `AAAAAAAAAA`), repeatable for multi-adapter (`-a SEQ1 -a SEQ2`), or `file:adapters.fa` to load from a FASTA file. The Perl-era embedded-string syntax (`-a " SEQ -a SEQ"`) still parses for back-compat.

## Bundled FastQC dependency

`fastqc-rust` is exact-pinned (`=1.0.1`). The pin is deliberate: the crate is brand-new (v1.0.0 published 2026-04-26) and any bump should be taken as a deliberate-test event with output-byte-identity re-verified against Java FastQC 0.12.1. The CI `validation` job has a smoke test that asserts (a) `fastqc` is **not** on `$PATH` (so a green run proves the bundled library did the work) and (b) the produced `.zip` contains `summary.txt`, `fastqc_data.txt`, `Images/`, `Icons/`.

## Conventions worth knowing

- **No external runtime deps.** Anything that would shell out to `cutadapt`, `pigz`, `fastqc`, or `java` is a regression — the v2.x story is "single static binary".
- **CI is `-D warnings`.** New clippy warnings will fail CI; fix them rather than `#[allow]`-ing without justification.
- **Validation matrix is load-bearing.** The CI `validation` job md5-checks Oxidized output against Perl 0.6.11 for SE, PE, hardtrim5, clock, and demux. If you change any of those code paths, expect to either preserve byte-identity or update the validation job with an explicit reason.
- **Output-collision pre-flight.** Don't bypass it — it catches issue #216 (case-only aliases on APFS/NTFS silently overwriting).
- **Em-dashes in user-facing strings.** `--version` provenance uses literal em-dashes; CI grep is content-targeted on that character.
