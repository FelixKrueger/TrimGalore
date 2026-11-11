#!/usr/bin/env python3
"""Phase 2 calibration for #439: fit the clumpify memory model on a
distinct-read fixture, as an upper envelope rather than a mean.

Design notes (all of these are review findings, not preferences):

* `n_bins = max(16, 4*cores)`, so the bin-pool coefficient (sigma) and the
  per-core coefficient (k) are separately identifiable only where n_bins is
  pinned. `--clumpify` requires `--cores >= 2` (cli.rs), so the identifying
  core counts are {2, 3, 4} — all give n_bins = 16.
* cores {8, 16} are FALSIFICATION cells: they are predicted from the fit, not
  fitted, because there n_bins moves with cores.
* The promise is one-sided ("actual <= predicted"), so a least-squares mean is
  the wrong estimator. Each cell is summarised by its MAX over reps, the fit is
  taken over those maxima, and an explicit margin is reported.
* Ratios are computed in MiB on both sides. The banner prints MiB.
* Every cell asserts the banner appeared: a silent fallback to plain mode would
  report a beautiful ratio of ~0.1.
* The fixture's distinct-read count is recorded in the output header. The
  previous calibration was invalidated by a fixture of 30,000 distinct reads
  cycled 227x, and nothing logged the sample size.

usage: calibrate.py <binary> <r1> <r2> <outdir> [tag]
"""

import os
import re
import shutil
import subprocess
import sys
import time

MIB = 1024 * 1024
STATIC_MIB = 512  # shipped STATIC_OVERHEAD_BYTES
BANNER = re.compile(r"clumpify: (\d+) bins × (\d+) MiB; predicted peak ≈ (\d+) MiB")
REPS = 2

# (cores, [target bin sizes in MiB], role)
CELLS = [
    (2, [12.8, 19.0, 27.4, 37.9, 56.9], "fit"),
    (3, [12.8, 19.0, 27.4, 37.9, 56.9], "fit"),
    (4, [12.8, 19.0, 27.4, 37.9, 56.9], "fit"),
    (8, [19.0, 37.9], "falsify"),
    (16, [19.0], "falsify"),
]


def n_bins_for(cores):
    return max(16, 4 * cores)


def denom_for(cores):
    return 5 * n_bins_for(cores) + 7 * cores


def budget_mib_for(target_b, cores):
    return round(target_b * denom_for(cores) / 4) + STATIC_MIB


def exact_b_mib(budget_mib, cores):
    return 4 * (budget_mib - STATIC_MIB) / denom_for(cores)


def fixture_stats(path):
    """Reads, distinct sequences, mean length.

    Counted entirely in a shell pipeline: on Linux a child's `ru_maxrss` is
    inherited from the forking parent, so holding the sequences in this process
    would report the harness's own footprint as every measured peak.
    """
    out = subprocess.run(
        ["sh", "-c",
         f"gzip -dc {path!r} | awk 'NR%4==2{{n++; t+=length($0); print $0}} "
         f"END{{print \"__STATS__\", n, t > \"/dev/stderr\"}}' | sort -u | wc -l"],
        capture_output=True, text=True, check=True)
    distinct = int(out.stdout.strip())
    _, reads, total = out.stderr.strip().split()
    return int(reads), distinct, int(total) / int(reads)


def parent_footprint_mib():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
        MIB if sys.platform == "darwin" else 1024)


def run_cell(binary, r1, r2, cores, budget_mib, rep, workdir):
    outdir = os.path.join(workdir, f"c{cores}_{budget_mib}_{rep}")
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    log = os.path.join(workdir, f"log_c{cores}_{budget_mib}_{rep}.txt")
    argv = [binary, "--paired", "--cores", str(cores), "--clumpify",
            "--compression", "6", "--memory", f"{budget_mib}M", "-o", outdir, r1, r2]
    with open(log, "wb") as lf:
        t0 = time.monotonic()
        p = subprocess.Popen(argv, stdout=lf, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.monotonic() - t0
    code = os.waitstatus_to_exitcode(status)
    text = open(log, encoding="utf-8", errors="replace").read()
    m = BANNER.search(text)
    if code != 0 or m is None:
        raise SystemExit(f"cores={cores} budget={budget_mib}M rep={rep}: exit={code}, "
                         f"banner={'yes' if m else 'MISSING'} — see {log}")
    bins, bin_mib, pred_mib = (int(x) for x in m.groups())
    assert bins == n_bins_for(cores), (bins, n_bins_for(cores))
    rss_mib = ru.ru_maxrss / (MIB if sys.platform == "darwin" else 1024)
    out_bytes = sum(os.path.getsize(os.path.join(outdir, f))
                    for f in os.listdir(outdir) if f.endswith(".fq.gz"))
    shutil.rmtree(outdir, ignore_errors=True)
    return {"cores": cores, "budget": budget_mib, "bins": bins,
            "b_exact": exact_b_mib(budget_mib, cores), "b_banner": bin_mib,
            "pred": pred_mib, "rss": rss_mib, "ratio": rss_mib / pred_mib,
            "wall": wall, "cpu": ru.ru_utime + ru.ru_stime, "out": out_bytes}


def lstsq(points):
    n = len(points)
    sx = sum(p[0] for p in points); sy = sum(p[1] for p in points)
    sxx = sum(p[0] ** 2 for p in points); sxy = sum(p[0] * p[1] for p in points)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return slope, (sy - slope * sx) / n


def main():
    binary, r1, r2, workdir = sys.argv[1:5]
    tag = sys.argv[5] if len(sys.argv) > 5 else "cal"
    os.makedirs(workdir, exist_ok=True)

    reads, distinct, mean_len = fixture_stats(r1)
    footprint = parent_footprint_mib()
    print(f"# platform      {sys.platform}")
    print(f"# binary        {binary}")
    print(f"# fixture R1    {reads:,} reads, {distinct:,} distinct sequences "
          f"({distinct / reads * 100:.1f}%), mean length {mean_len:.0f}")
    print(f"# reps/cell     {REPS}, cell summary = MAX over reps (one-sided promise)")
    print(f"# harness rss   {footprint:.1f} MiB")
    # On Linux a child inherits the forking parent's peak, so a fat harness
    # reports its own footprint as every measured cell.
    if footprint > 200:
        sys.exit(f"harness footprint {footprint:.0f} MiB is too large to measure "
                 f"children reliably — keep the fixture scan out of this process")
    print()

    tsv = os.path.join(workdir, f"results_{tag}_{sys.platform}.tsv")
    rows = []
    with open(tsv, "w") as fh:
        fh.write("role\tcores\tn_bins\tbudget_mib\tb_exact_mib\tb_banner_mib\t"
                 "predicted_mib\trep\trss_mib\tratio\twall_s\tcpu_s\tout_bytes\n")
        print(f"{'role':8s} {'cores':>5s} {'bins':>4s} {'budget':>7s} {'B':>7s} "
              f"{'pred':>7s} {'maxRSS':>9s} {'ratio':>6s} {'wall':>7s}")
        for cores, targets, role in CELLS:
            for tb in targets:
                budget = budget_mib_for(tb, cores)
                reps = []
                for rep in range(1, REPS + 1):
                    r = run_cell(binary, r1, r2, cores, budget, rep, workdir)
                    reps.append(r)
                    fh.write("\t".join(str(x) for x in [
                        role, cores, r["bins"], budget, f"{r['b_exact']:.3f}",
                        r["b_banner"], r["pred"], rep, f"{r['rss']:.1f}",
                        f"{r['ratio']:.4f}", f"{r['wall']:.2f}", f"{r['cpu']:.2f}",
                        r["out"]]) + "\n")
                    fh.flush()
                worst = max(reps, key=lambda x: x["rss"])
                rows.append({**worst, "role": role,
                             "spread": max(x["rss"] for x in reps) - min(x["rss"] for x in reps)})
                print(f"{role:8s} {cores:>5d} {worst['bins']:>4d} {budget:>6d}M "
                      f"{worst['b_exact']:6.2f}M {worst['pred']:6d}M {worst['rss']:8.1f}M "
                      f"{worst['ratio']:6.3f} {worst['wall']:6.2f}s", flush=True)

    print("\n=== per-cores slope (n_bins pinned at 16 for cores 2,3,4) ===")
    slopes = {}
    for cores in (2, 3, 4):
        pts = [(r["b_exact"], r["rss"]) for r in rows if r["cores"] == cores]
        s, a = lstsq(pts)
        slopes[cores] = s
        print(f"  cores={cores}: RSS = {a:8.1f} MiB + {s:6.3f} * B   "
              f"(model charges {denom_for(cores) / 4:6.3f} * B)")

    # slope(c) = sigma*16 + k*c  ->  least squares over the three (c, slope) pairs
    k, intercept = lstsq([(c, s) for c, s in slopes.items()])
    sigma = intercept / 16
    print(f"\n=== fitted coefficients (upper-envelope cells) ===")
    print(f"  sigma (per bin)  = {sigma:5.2f}   (model assumes 1.25)")
    print(f"  k     (per core) = {k:5.2f}   (model assumes 1.75)")
    print(f"  honest denominator = 4*(sigma*n_bins + k*cores)")
    for cores in (2, 3, 4, 8, 16):
        d = 4 * (sigma * n_bins_for(cores) + k * cores)
        print(f"    cores={cores:2d}: {d:7.1f}  (shipped {denom_for(cores)})")

    print("\n=== falsification: predicted vs measured on unfitted cells ===")
    max_static = max(r["rss"] - r["b_exact"] * (sigma * r["bins"] + k * r["cores"])
                     for r in rows if r["role"] == "fit")
    print(f"  fixed term implied by the fitted cells (max): {max_static:7.1f} MiB")
    for r in rows:
        if r["role"] != "falsify":
            continue
        pred = max_static + r["b_exact"] * (sigma * r["bins"] + k * r["cores"])
        err = (r["rss"] - pred) / pred * 100
        print(f"  cores={r['cores']:2d} B={r['b_exact']:5.1f}M: model {pred:7.1f}M vs "
              f"measured {r['rss']:7.1f}M ({err:+5.1f}%)")

    print("\n=== margin required ===")
    print(f"  worst ratio over all cells: {max(r['ratio'] for r in rows):.3f} "
          f"(cores={max(rows, key=lambda r: r['ratio'])['cores']})")
    print(f"  largest within-cell spread: {max(r['spread'] for r in rows):.1f} MiB")
    print(f"\nwrote {tsv}")


if __name__ == "__main__":
    main()
