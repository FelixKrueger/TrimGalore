#!/usr/bin/env python3
"""Spike (#457) driver: is the per-allocation overhead resident in full, or does
it ride the text's 0..B oscillation?

Two experiments.

**reuse** — peak RSS cannot observe a decrease, so the discriminator is live-set
size versus cumulative allocation. The cycle run allocates CYCLES times the
pool's worth of records while never holding more than one pool live. If freed
pages are reused, peak tracks the pool; if every cycle took fresh pages from the
OS, peak would approach CYCLES x pool. That is a factor-of-8 signal and needs no
`ps` (which the sandbox blocks) and no mach FFI.

**cycle** — a pool of bins filled pseudo-randomly and flushed at `raw_bytes >=
budget`, exactly as `PairedBin::push` does. Peak RSS is taken at several bin
budgets, so the slope in B gives

    coef(r) = d(peak)/dB / n_bins = sigma_effective(r)

directly, with no static term and no worker term to confound it. Regressing
coef on 1/r then gives both coefficients at once:

    sigma(r) = a + s/r     intercept a, slope s

and `s` is the discriminator: if it matches the fully-resident 122-135 B/record
from SPIKE_slot_cost, the allocator overhead is resident in full; if it comes
back near 74 + a*alloc/2, it rides occupancy.

usage: spike_occupancy.py <repo_root> [reps]
"""

import os
import re
import shutil
import subprocess
import sys

READ_LENGTHS = (25, 51, 150)
BUDGETS_MIB = (12, 22, 40)
N_BINS = 16
CYCLES = 8
DARWIN = sys.platform == "darwin"
MIB = 1024 * 1024

# From SPIKE_slot_cost.md, fully resident: 148 B spine + allocator, per pair.
ALLOC_PER_PAIR = {25: 107.0, 36: 127.2, 51: 131.6, 100: 128.4, 150: 122.2}
S_FULLY_RESIDENT = {25: 122.49, 36: 132.60, 51: 134.81, 100: 133.20, 150: 130.07}
CYCLE_LINE = re.compile(
    r"MODE=cycle PAIRS=(\d+) READ_LEN=(\d+) N_BINS=(\d+) BUDGET=(\d+) "
    r"FLUSHES=(\d+) LIVE=(\d+) RESIDENT=(\d+)")


def analytic_r(pairs, read_len):
    total, width, lo = 0, 1, 0
    while lo < pairs:
        hi = min(pairs, 10 ** width)
        total += (hi - lo) * width
        lo, width = hi, width + 1
    return 34 + total / pairs + 2 * read_len + 5


def lstsq(xs, ys):
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return slope, (sy - slope * sx) / n


def run_cycle(binary, pairs, read_len, budget_bytes):
    p = subprocess.Popen([binary, "cycle", str(pairs), str(read_len),
                          str(N_BINS), str(budget_bytes)],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.stdout.read().decode()
    p.stdout.close()
    _, status, ru = os.wait4(p.pid, 0)
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit(f"cycle child failed at len={read_len} B={budget_bytes}")
    m = CYCLE_LINE.search(out)
    if m is None:
        raise SystemExit(f"could not parse cycle output: {out!r}")
    flushes = int(m.group(5))
    rss = ru.ru_maxrss if DARWIN else ru.ru_maxrss * 1024
    return rss, flushes


def main():
    repo = os.path.abspath(sys.argv[1])
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    here = os.path.dirname(os.path.abspath(__file__))
    staged = os.path.join(repo, "examples", "spike_occupancy.rs")
    binary = os.path.join(repo, "target", "release", "examples", "spike_occupancy")

    if os.path.exists(staged):
        raise SystemExit(f"{staged} already exists — refusing to overwrite")
    shutil.copyfile(os.path.join(here, "spike_occupancy.rs"), staged)
    try:
        build = subprocess.run(
            ["cargo", "build", "--release", "--example", "spike_occupancy"],
            cwd=repo, capture_output=True, text=True)
        if build.returncode != 0:
            raise SystemExit("build failed:\n" + build.stderr[-3000:])
        print(f"# platform {sys.platform}, mimalloc, n_bins {N_BINS}, "
              f"{CYCLES} flush cycles per bin, reps {reps}\n")

        print("=== cycle: sigma_effective(r) from the slope in B ===")
        print(f"{'len':>4s} {'r':>7s} {'flushes':>8s} {'coef':>7s} "
              f"{'shipped 1.4375':>15s} {'plan a+74/r':>12s}")
        coefs, reuse = {}, {}
        for read_len in READ_LENGTHS:
            bs, peaks, flushes_seen = [], [], []
            for bmib in BUDGETS_MIB:
                budget = bmib * MIB
                bt = 2 * analytic_r(1_000_000, read_len)
                pairs = int(N_BINS * budget * CYCLES / bt)
                per_rep = []
                for _ in range(reps):
                    rss, fl = run_cycle(binary, pairs, read_len, budget)
                    per_rep.append(rss)
                    flushes_seen.append(fl)
                bs.append(budget)
                peaks.append(max(per_rep))
            slope, _ = lstsq(bs, peaks)
            coef = slope / N_BINS
            r = analytic_r(int(N_BINS * BUDGETS_MIB[-1] * MIB * CYCLES
                               / (2 * analytic_r(1_000_000, read_len))), read_len)
            coefs[read_len] = (coef, r)
            # One pool at full occupancy, for the reuse check: text plus the
            # fully-resident spine and allocator cost from SPIKE_slot_cost.
            big = BUDGETS_MIB[-1] * MIB
            slots = big / (2 * r)
            pool_full = N_BINS * (big + (148 + ALLOC_PER_PAIR[read_len]) * slots)
            reuse[read_len] = (coef, r, peaks[-1], pool_full)
            plan_sigma = 0.937 + 74 / r
            print(f"{read_len:>4d} {r:>7.1f} {min(flushes_seen):>8d} {coef:>7.3f} "
                  f"{coef - 1.4375:>+15.3f} {coef - plan_sigma:>+12.3f}", flush=True)

        print("\n=== coefficients, from sigma(r) = a + s/r ===")
        xs = [1 / v[1] for v in coefs.values()]
        ys = [v[0] for v in coefs.values()]
        s_fit, a_fit = lstsq(xs, ys)
        print(f"  a (occupancy)          = {a_fit:6.3f}   plan 0.937, must be 0 < a <= 1")
        print(f"  s (B per record)       = {s_fit:6.2f}   plan 74")
        mean_full = sum(S_FULLY_RESIDENT[k] for k in coefs) / len(coefs)
        print(f"  fully-resident s (spike 1) = {mean_full:6.2f} over the same lengths")
        print(f"\n=== verdict on the occupancy question ===")
        if abs(s_fit - mean_full) < 8:
            print(f"  s matches the fully-resident value ({s_fit:.1f} vs {mean_full:.1f})")
            print("  => the allocator overhead is RESIDENT IN FULL; charge it at "
                  "full slots, not a*slots")
        else:
            implied = 74 + (mean_full - 74) * max(a_fit, 0.0)
            print(f"  s = {s_fit:.1f} vs fully-resident {mean_full:.1f}; "
                  f"a-scaled would predict {implied:.1f}")
            print("  => the overhead does NOT sit at full residency; it is "
                  "attenuated like the text")
        print("\n=== reuse check (peak vs pool, at the largest budget) ===")
        for read_len, (coef, r, peak, pool) in reuse.items():
            print(f"  {read_len:3d} bp: peak {peak / MIB:7.1f}M, one pool "
                  f"{pool / MIB:7.1f}M, {CYCLES}x pool "
                  f"{CYCLES * pool / MIB:8.1f}M -> peak/pool = {peak / pool:.2f}")
        print("  peak/pool near 1 means freed pages are reused; near "
              f"{CYCLES} would mean every cycle took fresh pages")
    finally:
        if os.path.exists(staged):
            os.remove(staged)
            print(f"\n# removed staged {staged}")


if __name__ == "__main__":
    main()
