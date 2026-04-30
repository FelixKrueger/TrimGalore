#!/usr/bin/env python3
"""Generate benchmark charts from hyperfine JSON outputs.

Reads `<data>/{rust-beta5,rust-beta7,perl-0.6.11}_cN.json` (the per-condition
hyperfine outputs from `scripts/benchmark.sh`), produces three PNGs into
`<out>/`:

    benchmark_wall_time.png    — wall time vs cores, three engines
    benchmark_threads_cpu.png  — CPU time vs cores, three engines
    benchmark_scaling.png      — speedup factor vs Perl `-j 1` baseline

Conditions that haven't completed yet (zero-byte JSONs) are skipped silently,
so this can be re-run as the matrix fills in.

Usage:
    python3 docs/scripts/generate-benchmark-charts.py \\
        --data /path/to/perf_data/2026-04-29 \\
        --out  docs/src/assets/benchmarks
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# Use a writable matplotlib config dir (avoids the default ~/.matplotlib
# that may not be writable in some sandboxed environments).
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"),
)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ----- Engine display config ---------------------------------------------

ENGINES = {
    "perl-0.6.11": {
        "label": "Perl 0.6.11 + Cutadapt 5.2",
        "color": "#d62728",
        "marker": "o",
    },
    "rust-beta5": {
        "label": "Rust v2.1.0-beta.5 (pre-audit)",
        "color": "#7f7f7f",
        "marker": "D",
    },
    "rust-beta7": {
        "label": "Rust v2.1.0-beta.7 (current)",
        "color": "#1f77b4",
        "marker": "s",
    },
}
ENGINE_ORDER = ["perl-0.6.11", "rust-beta5", "rust-beta7"]


# ----- Data loading ------------------------------------------------------

JSON_NAME_RE = re.compile(r"^(.+)_c(\d+)\.json$")


def load_data(data_dir: Path) -> dict[str, dict[int, dict[str, float]]]:
    data: dict[str, dict[int, dict[str, float]]] = defaultdict(dict)
    for jp in sorted(data_dir.glob("*.json")):
        if jp.stat().st_size == 0:
            continue
        m = JSON_NAME_RE.match(jp.name)
        if not m or m.group(1) not in ENGINES:
            continue
        engine, cores = m.group(1), int(m.group(2))
        try:
            payload = json.loads(jp.read_text())
        except json.JSONDecodeError:
            continue
        results = payload.get("results")
        if not results:
            continue
        r = results[0]
        data[engine][cores] = {
            "wall_mean": r["mean"],
            "wall_stddev": r["stddev"],
            "cpu_total": r["user"] + r["system"],
        }
    return data


def fmt_mmss(seconds: float) -> str:
    m, s = divmod(round(seconds), 60)
    return f"{m}:{s:02d}"


# ----- Common figure styling ---------------------------------------------

NF_CORE_DEFAULT = 8


def style_xaxis(ax, all_cores: list[int]) -> None:
    """Linear-spaced cores axis with explicit tick labels at the matrix points."""
    ticks = sorted(set(all_cores))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel("Cores requested  (-j  /  --cores)")
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_axisbelow(True)


def add_nfcore_guide(ax, ymax: float, label_y_frac: float = 0.95) -> None:
    ax.axvline(NF_CORE_DEFAULT, color="gray", linestyle=":", alpha=0.5, linewidth=1)
    ax.text(
        NF_CORE_DEFAULT + 0.15,
        ymax * label_y_frac,
        "nf-core default",
        color="gray",
        fontsize=8,
        style="italic",
    )


# ----- Chart 1: wall time ------------------------------------------------


def plot_wall_time(data, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6))

    all_cores: list[int] = []
    ymax = 0.0
    for engine in ENGINE_ORDER:
        if engine not in data:
            continue
        cores = sorted(data[engine].keys())
        walls = [data[engine][c]["wall_mean"] for c in cores]
        cfg = ENGINES[engine]
        ax.plot(
            cores,
            walls,
            marker=cfg["marker"],
            markersize=8,
            linewidth=2,
            label=cfg["label"],
            color=cfg["color"],
        )
        all_cores.extend(cores)
        ymax = max(ymax, max(walls))

        # Inline mm:ss callouts at c1 and c8 if present.
        for c in (1, 8):
            if c in data[engine]:
                w = data[engine][c]["wall_mean"]
                ax.annotate(
                    f"{fmt_mmss(w)} min",
                    xy=(c, w),
                    xytext=(c + 1.3, w + ymax * 0.04),
                    fontsize=9,
                    fontweight="bold",
                    color=cfg["color"],
                    arrowprops=dict(arrowstyle="-", color=cfg["color"], alpha=0.5),
                )

    ax.set_ylabel("Wall time (seconds)")
    ax.set_title(
        "Wall Time — 84M PE Reads (Buckberry 2023, Xeon 6975P-C)",
        fontweight="bold",
    )
    style_xaxis(ax, all_cores)
    add_nfcore_guide(ax, ymax)
    ax.set_ylim(0, ymax * 1.12)
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


# ----- Chart 2: CPU time -------------------------------------------------


def plot_cpu_time(data, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6))

    all_cores: list[int] = []
    ymax = 0.0
    for engine in ENGINE_ORDER:
        if engine not in data:
            continue
        cores = sorted(data[engine].keys())
        cpus = [data[engine][c]["cpu_total"] for c in cores]
        cfg = ENGINES[engine]
        ax.plot(
            cores,
            cpus,
            marker=cfg["marker"],
            markersize=8,
            linewidth=2,
            label=cfg["label"],
            color=cfg["color"],
        )
        all_cores.extend(cores)
        ymax = max(ymax, max(cpus))

    # Annotate CPU savings vs Perl at cores=8 if all three are present.
    if (
        "perl-0.6.11" in data
        and 8 in data["perl-0.6.11"]
        and "rust-beta7" in data
        and 8 in data["rust-beta7"]
    ):
        perl_c8 = data["perl-0.6.11"][8]["cpu_total"]
        beta7_c8 = data["rust-beta7"][8]["cpu_total"]
        savings = perl_c8 / beta7_c8
        ax.annotate(
            f"{savings:.1f}× less\nCPU at c8",
            xy=(8, beta7_c8),
            xytext=(9.5, (perl_c8 + beta7_c8) / 2),
            fontsize=10,
            fontweight="bold",
            color="#2ca02c",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#2ca02c"),
            arrowprops=dict(arrowstyle="->", color="#2ca02c"),
        )

    ax.set_ylabel("CPU time, user + system  (seconds)")
    ax.set_title(
        "CPU Time — 84M PE Reads (Buckberry 2023, Xeon 6975P-C)",
        fontweight="bold",
    )
    style_xaxis(ax, all_cores)
    add_nfcore_guide(ax, ymax)
    ax.set_ylim(0, ymax * 1.12)
    # Perl c1 sits in the upper-left at ~4400s, so anchor the legend
    # to the right where the Rust lines stay flat (~700s max).
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


# ----- Chart 3: scaling --------------------------------------------------


def plot_scaling(data, outpath: Path) -> None:
    """Speedup vs Perl -j 1 baseline. Every engine is normalized to perl_c1."""
    if "perl-0.6.11" not in data or 1 not in data["perl-0.6.11"]:
        print("skipping scaling plot — perl c1 baseline missing")
        return
    baseline = data["perl-0.6.11"][1]["wall_mean"]

    fig, ax = plt.subplots(figsize=(11, 6))

    all_cores: list[int] = []
    ymax = 1.0
    for engine in ENGINE_ORDER:
        if engine not in data:
            continue
        cores = sorted(data[engine].keys())
        speedups = [baseline / data[engine][c]["wall_mean"] for c in cores]
        cfg = ENGINES[engine]
        ax.plot(
            cores,
            speedups,
            marker=cfg["marker"],
            markersize=8,
            linewidth=2,
            label=cfg["label"],
            color=cfg["color"],
        )
        for c, s in zip(cores, speedups):
            ax.annotate(
                f"{s:.1f}×",
                xy=(c, s),
                xytext=(0, 9),
                textcoords="offset points",
                ha="center",
                fontsize=8,
                fontweight="bold",
                color=cfg["color"],
            )
        all_cores.extend(cores)
        ymax = max(ymax, max(speedups))

    # Ideal linear scaling (dashed gray).
    ideal_max = max(all_cores) if all_cores else 24
    ax.plot(
        [1, ideal_max],
        [1, ideal_max],
        linestyle="--",
        color="lightgray",
        linewidth=1.2,
        label="Ideal linear scaling",
    )

    ax.set_ylabel("Speedup vs Perl `-j 1`")
    ax.set_title(
        "Scaling — Speedup Relative to Perl `-j 1` (84M PE Reads, Xeon 6975P-C)",
        fontweight="bold",
    )
    style_xaxis(ax, all_cores)
    add_nfcore_guide(ax, ymax)
    ax.set_ylim(0, ymax * 1.15)
    ax.legend(loc="upper left", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


# ----- Entrypoint --------------------------------------------------------


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, type=Path, help="Directory with hyperfine JSON outputs")
    p.add_argument("--out", required=True, type=Path, help="Directory to write PNG files into")
    args = p.parse_args()

    data = load_data(args.data)
    if not data:
        print(f"no data found in {args.data}", file=sys.stderr)
        return 1

    print(f"loaded data for engines: {sorted(data.keys())}")
    for engine in ENGINE_ORDER:
        if engine in data:
            print(f"  {engine}: cores {sorted(data[engine].keys())}")

    args.out.mkdir(parents=True, exist_ok=True)
    plot_wall_time(data, args.out / "benchmark_wall_time.png")
    plot_cpu_time(data, args.out / "benchmark_threads_cpu.png")
    plot_scaling(data, args.out / "benchmark_scaling.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
