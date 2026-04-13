# ── Build stage ──────────────────────────────────────────────
FROM rust:1.85-bookworm AS builder

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

RUN apt-get update && apt-get install -y --no-install-recommends \
      procps \
      default-jre-headless \
      perl \
      unzip \
      wget \
      ca-certificates \
    && wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip -O /tmp/fastqc.zip \
    && unzip -q /tmp/fastqc.zip -d /opt \
    && chmod +x /opt/FastQC/fastqc \
    && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
    && rm /tmp/fastqc.zip \
    && apt-get purge -y unzip wget \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/trim_galore /usr/local/bin/trim_galore
