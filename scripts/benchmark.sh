#!/usr/bin/env bash
# Trim Galore — Buckberry-scale benchmark driver.
#
# Reproducible perf comparison across:
#   - Perl Trim Galore 0.6.11 + Cutadapt 5.2 (cores 1, 4, 8, 16; igzip + pigz)
#   - Rust v2.1.0-beta.5 (cores 1, 4, 8, 16, 24)  — pre-Buckberry-audit
#   - Rust v2.1.0-beta.7 (cores 1, 4, 8, 16, 24)  — post-Myers' prefilter
#
# Methodology mirrors @an-altosian's #248 audit:
#   - hyperfine --warmup 1 --runs 10 per condition
#   - --export-json per condition for downstream plotting / aggregation
#   - filesystem cache warmed once up front via cat ... > /dev/null
#   - md5 sanity-check at the end: beta.5 vs beta.7 decompressed output
#     must match (Myers' is byte-identity-preserving by construction)
#
# Designed to run on the Altos `oxy` EKS instance (`dcli ssh oxy`)
# where `/home/fkrueger/benchmark_TG_oxy/` is the canonical fixture
# directory. Override via env vars if running elsewhere — see below.
#
# Wall-time estimate at Buckberry scale (84M PE reads, 2.1 GiB
# gzipped): roughly 13 hours total. Run under tmux/nohup. Skip the
# Perl cores=1 row (PERL_CORES="4 8 16") to save ~5 hours if you've
# already got that datapoint from a previous run.

set -euo pipefail

# -- Config (override via env) --

# Fixture: Buckberry SRR24827373, the Dongze-#248-audit dataset
DATA_DIR="${DATA_DIR:-/home/fkrueger/benchmark_TG_oxy}"
R1="${R1:-$DATA_DIR/SRR24827373_GSM7445361_32F_NB3_p28_p2n2p_p10_rep1_Homo_sapiens_Bisulfite-Seq_R1.fastq.gz}"
R2="${R2:-$DATA_DIR/SRR24827373_GSM7445361_32F_NB3_p28_p2n2p_p10_rep1_Homo_sapiens_Bisulfite-Seq_R2.fastq.gz}"

# Binaries -- install via:
#   cargo install trim-galore --version 2.1.0-beta.5 --root ~/.cargo-tg-beta5 --locked
#   cargo install trim-galore --version 2.1.0-beta.7 --root ~/.cargo-tg-beta7 --locked
# (Override paths if you've installed elsewhere or built from source.)
RUST_BETA5="${RUST_BETA5:-$HOME/.cargo-tg-beta5/bin/trim_galore}"
RUST_BETA7="${RUST_BETA7:-$HOME/.cargo-tg-beta7/bin/trim_galore}"
PERL_TG="${PERL_TG:-$HOME/TrimGalore/trim_galore}"

# Cores ladder
RUST_CORES="${RUST_CORES:-1 4 8 16 24}"
PERL_CORES="${PERL_CORES:-1 4 8 16}"

# Hyperfine settings -- match Dongze's #248 methodology
WARMUP="${WARMUP:-1}"
RUNS="${RUNS:-10}"

# Output dirs
WORK_ROOT="${WORK_ROOT:-/tmp/tg_bench}"
RESULTS_DIR="${RESULTS_DIR:-$PWD/perf_data/$(date -u +%Y-%m-%d)}"

# -- Pre-flight --

for tool in hyperfine md5sum gzip; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not on PATH"; exit 1; }
done
[[ -r "$R1" && -r "$R2" ]] || { echo "ERROR: cannot read R1/R2 ($R1, $R2)"; exit 1; }
[[ -x "$RUST_BETA5" ]] || { echo "ERROR: Rust beta.5 binary not executable at $RUST_BETA5"; exit 1; }
[[ -x "$RUST_BETA7" ]] || { echo "ERROR: Rust beta.7 binary not executable at $RUST_BETA7"; exit 1; }
[[ -x "$PERL_TG" ]] || { echo "ERROR: Perl trim_galore not executable at $PERL_TG"; exit 1; }

mkdir -p "$WORK_ROOT" "$RESULTS_DIR"

echo "=========================================="
echo "Trim Galore benchmark driver"
echo "=========================================="
echo "Date:              $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Host:              $(hostname)"
echo "Cores (host):      $(nproc)"
echo "CPU model:         $(grep 'model name' /proc/cpuinfo | head -1 | sed 's/.*: //')"
echo "Memory:            $(free -h | awk '/^Mem:/{print $2}')"
echo "R1:                $R1"
echo "R2:                $R2"
echo "Rust beta.5:       $RUST_BETA5 ($($RUST_BETA5 --version | head -1))"
echo "Rust beta.7:       $RUST_BETA7 ($($RUST_BETA7 --version | head -1))"
echo "Perl trim_galore:  $PERL_TG ($($PERL_TG --version 2>&1 | grep -i 'version:' | head -1 | tr -s ' '))"
echo "Cutadapt:          $(cutadapt --version 2>&1 | head -1)"
echo "hyperfine:         $(hyperfine --version | head -1)"
echo "Cores ladder:      Rust=[$RUST_CORES]  Perl=[$PERL_CORES]"
echo "Hyperfine:         --warmup $WARMUP --runs $RUNS"
echo "Results dir:       $RESULTS_DIR"
echo "=========================================="
echo ""

# -- Cache warming --
echo "[1/5] Warming filesystem cache (one-time, ~30s on EFS)..."
time { cat "$R1" > /dev/null; cat "$R2" > /dev/null; }
echo ""

# -- Helpers --

# Run one hyperfine sweep for a single (binary, cores, label) tuple.
# Each iteration cleans the output dir so we measure end-to-end work,
# not no-op overwrites.
run_sweep() {
    local label="$1"; shift
    local cores="$1"; shift
    local cmd_template="$1"; shift  # e.g. "$BIN --paired --cores {cores} -o {outdir} $R1 $R2"

    local outdir="$WORK_ROOT/${label}_c${cores}"
    local json="$RESULTS_DIR/${label}_c${cores}.json"
    local md="$RESULTS_DIR/${label}_c${cores}.md"

    echo "----- ${label} cores=${cores} -----"
    rm -rf "$outdir"
    mkdir -p "$outdir"

    # Substitute {outdir} and {cores} in the template.
    local cmd="${cmd_template//\{outdir\}/$outdir}"
    cmd="${cmd//\{cores\}/$cores}"

    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --prepare "rm -rf $outdir && mkdir -p $outdir" \
        --command-name "${label}-c${cores}" \
        --export-json "$json" \
        --export-markdown "$md" \
        "$cmd"

    echo ""
}

# -- Sweeps --

echo "[2/5] Rust v2.1.0-beta.5 (pre-Buckberry audit)"
for c in $RUST_CORES; do
    run_sweep "rust-beta5" "$c" "$RUST_BETA5 --paired --cores {cores} -o {outdir} $R1 $R2"
done

echo "[3/5] Rust v2.1.0-beta.7 (post-Myers' prefilter)"
for c in $RUST_CORES; do
    run_sweep "rust-beta7" "$c" "$RUST_BETA7 --paired --cores {cores} -o {outdir} $R1 $R2"
done

echo "[4/5] Perl Trim Galore 0.6.11 + Cutadapt 5.2 (igzip + pigz)"
for c in $PERL_CORES; do
    run_sweep "perl-0.6.11" "$c" "$PERL_TG --paired -j {cores} -o {outdir} $R1 $R2"
done

# -- Cross-version byte-identity sanity --

echo "[5/5] Byte-identity cross-check: beta.5 vs beta.7 (decompressed)"

# Pick cores=8 outputs from each -- should be byte-identical when
# decompressed (Myers' is a transparent prefilter; gzip-level is
# the same in both releases).
B5_DIR="$WORK_ROOT/rust-beta5_c8"
B7_DIR="$WORK_ROOT/rust-beta7_c8"
SUMMARY="$RESULTS_DIR/byte_identity_summary.txt"
{
    echo "Cross-version byte-identity check (cores=8 outputs)"
    echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""
    for f in $(ls "$B5_DIR" | grep -E '\.fq\.gz$'); do
        if [[ -f "$B7_DIR/$f" ]]; then
            B5_MD5=$(gzip -dc "$B5_DIR/$f" | md5sum | awk '{print $1}')
            B7_MD5=$(gzip -dc "$B7_DIR/$f" | md5sum | awk '{print $1}')
            if [[ "$B5_MD5" == "$B7_MD5" ]]; then
                echo "MATCH    $f  $B5_MD5"
            else
                echo "MISMATCH $f  beta.5=$B5_MD5  beta.7=$B7_MD5"
            fi
        else
            echo "MISSING  $f  (only in beta.5)"
        fi
    done
} | tee "$SUMMARY"

echo ""
echo "=========================================="
echo "DONE. Results landed in: $RESULTS_DIR"
echo "  - per-condition: <label>_c<cores>.{json,md}"
echo "  - byte-identity: byte_identity_summary.txt"
echo "Run 'tar czf bench-$(date -u +%Y%m%d).tar.gz $RESULTS_DIR' to ship."
echo "=========================================="
