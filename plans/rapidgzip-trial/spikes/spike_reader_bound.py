#!/usr/bin/env python3
"""Spike: is TrimGalore's parallel path reader-bound on gzipped input?

QUESTION
    TrimGalore decompresses input on ONE background thread (fastq.rs:276
    `open_threaded` -> flate2 MultiGzDecoder, feeding a sync_channel(4)), while
    trimming AND output compression run on N workers. Does that single
    decompression thread cap throughput, and from what --cores?

DESIGN — attribution by differential input format
    Run identical reads twice: once as .fastq.gz, once as plain .fastq.
    If the gzipped arm plateaus while the plain arm keeps scaling, the serial
    decompression thread is the ceiling. No source instrumentation needed.

CONTROLS (each of these would otherwise confound the result)
    1. --dont_gzip in BOTH arms. TrimGalore's output compression follows the
       input's by default (CHANGELOG #245), so without this the plain arm would
       also skip OUTPUT compression -- two variables, not one.
    2. Explicit -a. Adapter auto-detection scans up to 1M reads FROM THE INPUT,
       so it is itself sensitive to decompression speed.
    3. --cores 1 is reported but excluded from the scaling fit: main.rs routes
       `cores > 1` to the worker pool and cores==1 to the sequential path, so
       they are different architectures.
    4. Single-end only. Paired-end spawns TWO reader threads, which doubles
       decompression throughput and muddies the attribution.

NOT IN SCOPE
    Integrating rapidgzip; the uBAM/BGZF path; output-compression tuning.
"""

import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TG = REPO / "target" / "release" / "trim_galore"
BASE_GZ = REPO / "10K_150bp.fastq.gz"          # 10k x 150bp, the realistic-length fixture
SCRATCH = Path(os.environ.get("SPIKE_SCRATCH", "/tmp")) / "spike_reader_bound"

COPIES = int(os.environ.get("SPIKE_COPIES", "80"))   # ~290 MB uncompressed
CORE_LEVELS = [1, 2, 4, 8, 10]
REPEATS = int(os.environ.get("SPIKE_REPEATS", "2"))  # min-of-N, wall clock
ADAPTER = "AGATCGGAAGAGC"                            # pin it; skip auto-detection


def build_inputs() -> tuple[Path, Path, int, int]:
    """Concatenate the fixture into a large .gz (multi-member) and a plain .fastq."""
    SCRATCH.mkdir(parents=True, exist_ok=True)
    # Cache key MUST include COPIES: keying on the bare name silently reused a
    # smaller input when COPIES changed, and the run reported the stale size.
    gz = SCRATCH / f"big_{COPIES}.fastq.gz"
    plain = SCRATCH / f"big_{COPIES}.fastq"

    if not gz.exists():
        # Concatenated gzip members are a valid gzip stream (RFC 1952) and are
        # what TrimGalore's own parallel writer produces, so this also exercises
        # the multi-member path rapidgzip has a fast path for.
        with gz.open("wb") as out:
            src = BASE_GZ.read_bytes()
            for _ in range(COPIES):
                out.write(src)

    if not plain.exists():
        with plain.open("wb") as out:
            subprocess.run(["gzip", "-dc", str(gz)], stdout=out, check=True)

    reads = sum(1 for _ in plain.open("rb")) // 4
    expected = COPIES * 10_000
    if reads != expected:
        sys.exit(f"input built wrong: {reads:,} reads, expected {expected:,} "
                 f"({COPIES} copies x 10k). Stale cache?")
    return gz, plain, reads, plain.stat().st_size


def run_once(inp: Path, cores: int, outdir: Path) -> float:
    """Wall-clock seconds for one trim run. Raises on non-zero exit."""
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    cmd = [
        str(TG),
        "-a", ADAPTER,        # control 2: no auto-detection scan
        "--dont_gzip",        # control 1: identical output format in both arms
        "--cores", str(cores),
        "-o", str(outdir),
        str(inp),
    ]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"FAILED (cores={cores}, {inp.name}):\n{r.stderr[-2000:]}")
    return dt


def main() -> None:
    if not TG.exists():
        sys.exit(f"binary missing: {TG}\nrun: cargo build --release")
    if not BASE_GZ.exists():
        sys.exit(f"fixture missing: {BASE_GZ}")

    print(f"repo        : {REPO}")
    print(f"cpus        : {os.cpu_count()}")
    print(f"copies      : {COPIES}   repeats: {REPEATS} (min wall-clock)")

    gz, plain, reads, unc_bytes = build_inputs()
    print(f"input       : {reads:,} reads, {unc_bytes/1e6:.0f} MB uncompressed, "
          f"{gz.stat().st_size/1e6:.0f} MB gzipped\n")

    results: dict[str, dict[int, float]] = {"gz": {}, "plain": {}}
    for cores in CORE_LEVELS:
        for label, inp in (("gz", gz), ("plain", plain)):
            times = [run_once(inp, cores, SCRATCH / f"out_{label}_{cores}")
                     for _ in range(REPEATS)]
            best = min(times)
            results[label][cores] = best
            mbps = unc_bytes / 1e6 / best
            spread = (max(times) - best) / best * 100 if best else 0
            print(f"  cores={cores:<3} {label:<5} {best:6.2f}s  "
                  f"{mbps:6.1f} MB/s  (spread {spread:.1f}%)")

    # ---- analysis ----
    print("\n" + "=" * 66)
    print("SCALING (worker-pool path only; cores=1 is the sequential path)")
    print("=" * 66)
    pool = [c for c in CORE_LEVELS if c > 1]
    base = {k: results[k][pool[0]] for k in results}
    print(f"{'cores':>6} {'gz s':>8} {'gz x':>6} {'plain s':>9} {'plain x':>8} "
          f"{'gz/plain':>9}")
    for c in pool:
        g, p = results["gz"][c], results["plain"][c]
        print(f"{c:>6} {g:>8.2f} {base['gz']/g:>6.2f} {p:>9.2f} "
              f"{base['plain']/p:>8.2f} {g/p:>9.2f}")

    hi = pool[-1]
    gz_speedup = base["gz"] / results["gz"][hi]
    pl_speedup = base["plain"] / results["plain"][hi]
    overhead = (results["gz"][hi] / results["plain"][hi] - 1) * 100
    ideal = hi / pool[0]

    print("\n" + "=" * 66)
    print("VERDICT")
    print("=" * 66)
    print(f"  ideal speedup {pool[0]}->{hi} cores : {ideal:.2f}x")
    print(f"  gzipped input achieved         : {gz_speedup:.2f}x")
    print(f"  plain input achieved           : {pl_speedup:.2f}x")
    print(f"  gz slower than plain at {hi} cores : {overhead:+.1f}%")

    # Decision rule stated up front, applied mechanically.
    if pl_speedup > gz_speedup * 1.15 and overhead > 15:
        verdict = ("READER-BOUND. The gzipped arm scales materially worse than "
                   "the plain arm and costs a large constant penalty at high "
                   "core counts -> parallel decompression (rapidgzip) has "
                   "headroom.")
    elif overhead < 5:
        verdict = ("NOT reader-bound. Gzipped and plain input perform within "
                   "5% at max cores -> decompression is not the ceiling and "
                   "rapidgzip would buy little or nothing.")
    else:
        verdict = ("INCONCLUSIVE / partially reader-bound. Some penalty from "
                   "decompression but it is not the dominant term; needs a "
                   "bigger input or per-stage instrumentation before "
                   "committing.")
    print(f"\n  => {verdict}")
    print(f"\nscratch: {SCRATCH}  (delete when done)")


if __name__ == "__main__":
    main()
