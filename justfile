# Trim Galore — local CI parity
#
# Run `just` (or `just default`) to run the same checks CI runs on every
# push, in the same order. `just <target>` to run individual checks.
#
# Install just via `cargo install just` or your package manager.
# See https://github.com/casey/just for more.
#
# Notes:
# - All targets are run from the crate root; `just` automatically `cd`s here.
# - The validate target requires `cutadapt`, `conda`, and a Perl Trim Galore
#   on PATH; without those, only `just lint test` is portable. Use `just ci`
#   to run the portable subset locally before pushing.

# Default: run the full portable check suite (lint + tests in both profiles).
default: ci

# Portable CI parity: everything `rust-tests` and `lint` jobs run on PR.
ci: fmt clippy test test-release
    @echo "✓ all portable CI checks passed"

# Format check (matches CI's `cargo fmt --all -- --check`).
fmt:
    cargo fmt --all -- --check

# Apply formatting fixes in-place.
fmt-fix:
    cargo fmt --all

# Clippy with -D warnings (matches CI's lint job).
clippy:
    cargo clippy --all-targets --release -- -D warnings

# Debug-profile test suite.
test:
    cargo test

# Release-profile test suite (matches the new `Run tests (release)` step).
test-release:
    cargo test --release

# Build the release binary (target/release/trim_galore).
build:
    cargo build --release

# Run the bundled binary's --version (provenance + git-hash check).
version: build
    ./target/release/trim_galore --version

# Reproducibility check: build twice with the same SOURCE_DATE_EPOCH and
# diff the binaries. Mirrors the `reproducibility` CI job.
reproduce:
    @echo "Building binary twice under SOURCE_DATE_EPOCH=1700000000..."
    cargo clean
    SOURCE_DATE_EPOCH=1700000000 cargo build --release
    cp target/release/trim_galore /tmp/trim_galore.repro1
    cargo clean
    SOURCE_DATE_EPOCH=1700000000 cargo build --release
    cp target/release/trim_galore /tmp/trim_galore.repro2
    diff /tmp/trim_galore.repro1 /tmp/trim_galore.repro2 && echo "✓ bit-identical"

# Validation matrix vs Perl 0.6.11 — requires Perl trim_galore + cutadapt
# on PATH. The CI `validation` job orchestrates this with conda; use this
# target to spot-check a single comparison locally if those tools are
# already installed.
validate-paired-end:
    rm -rf /tmp/tg_validate /tmp/rust_validate
    mkdir -p /tmp/tg_validate /tmp/rust_validate
    trim_galore_perl --paired -o /tmp/tg_validate \
        test_files/BS-seq_10K_R1.fastq.gz test_files/BS-seq_10K_R2.fastq.gz
    cargo build --release
    ./target/release/trim_galore --paired -o /tmp/rust_validate \
        test_files/BS-seq_10K_R1.fastq.gz test_files/BS-seq_10K_R2.fastq.gz
    @for f in BS-seq_10K_R1_val_1.fq.gz BS-seq_10K_R2_val_2.fq.gz; do \
        TG_MD5=$(gzip -dc /tmp/tg_validate/$f | md5sum | cut -d' ' -f1); \
        RUST_MD5=$(gzip -dc /tmp/rust_validate/$f | md5sum | cut -d' ' -f1); \
        if [ "$TG_MD5" = "$RUST_MD5" ]; then \
            echo "✓ $f matches"; \
        else \
            echo "✗ $f DIFFERS (Perl=$TG_MD5 Rust=$RUST_MD5)"; exit 1; \
        fi; \
    done

# Regenerate PNG logos from the SVG sources (in Docs/).
logos:
    cd Docs && npm run logos

# Build the docs site (in Docs/).
docs:
    cd Docs && npm run build

# Serve the docs site locally for preview (in Docs/).
docs-dev:
    cd Docs && npm run dev

# Show all available recipes.
help:
    @just --list
