# ── Build stage ──────────────────────────────────────────────
FROM rust:1.88-bookworm AS builder

WORKDIR /build

# Cache dependency build: copy manifests first, build deps with dummy source
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs \
    && cargo build --release --locked \
    && rm -rf src

# Copy real source and rebuild (only the binary recompiles, deps are cached)
COPY . .
RUN cargo build --release --locked

# ── Runtime stage ────────────────────────────────────────────
FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/FelixKrueger/TrimGalore"
LABEL org.opencontainers.image.description="Trim Galore — adapter/quality trimming (Oxidized Edition)"
LABEL org.opencontainers.image.licenses="GPL-3.0-only"

# FastQC is now built into the trim_galore binary itself via the bundled
# fastqc-rust library — no Java runtime, no FastQC tarball, no perl.
# Only minimal diagnostic packages remain.
RUN apt-get update && apt-get install -y --no-install-recommends \
      procps \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/trim_galore /usr/local/bin/trim_galore
