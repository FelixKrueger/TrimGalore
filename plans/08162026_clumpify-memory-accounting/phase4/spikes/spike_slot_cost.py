#!/usr/bin/env python3
"""Spike (#457) driver: what is a paired-bin slot's non-text resident cost?

Method. `spike_slot_cost.rs` builds `PairedBin`'s exact shape and holds it past
the process peak. This runs it at four record counts per read length and takes
the SLOPE of peak RSS in record count. A slope cancels the process baseline, so
the answer does not depend on the static term — which is the whole reason the
calibration matrix could not settle this (static, `a` and `b` are degenerate at
one budget, and forcing `a`/`b` to structural values moves the excess into
`static` instead of `s`).

Two modes decompose the cost:

  spine  three empty Strings per record -> slope should be 72+72+4 = 148 B/pair
  full   real text                      -> slope is 148 + heap text + allocator

  heap text (exact) = 2*(r - 5)     `estimated_record_bytes` counts 5 B of '+'
                                    and newlines that exist on disk, not on heap
  allocator         = (full - spine) - heap text
  s (as clump.rs uses it) = (full slope - 2r) / 2

usage: spike_slot_cost.py <repo_root> [reps]
"""

import os
import re
import shutil
import subprocess
import sys

READ_LENGTHS = (25, 36, 51, 100, 150)
POINTS = (0.25, 0.5, 0.75, 1.0)
TARGET_BYTES = 1_200_000_000          # per run, so 150 bp stays ~1.2 GB
SPINE_SIZEOF = 148.0                  # 72 + 72 + 4
S_SIZEOF, S_FITTED = 74.0, 108.74
DARWIN = sys.platform == "darwin"
LINE = re.compile(r"PAIRS=(\d+) READ_LEN=(\d+) TEXT_BYTES=(\d+) MODE=(\w+) GUARD=(\d+)")


def analytic_r(pairs, read_len):
    """Mean bytes per record for this record count — the header carries the index."""
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


def run(binary, pairs, read_len, mode):
    p = subprocess.Popen([binary, str(pairs), str(read_len), mode],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.stdout.read().decode()
    p.stdout.close()
    _, status, ru = os.wait4(p.pid, 0)
    code = os.waitstatus_to_exitcode(status)
    if code != 0:
        raise SystemExit(f"spike exited {code} at pairs={pairs} len={read_len} {mode}")
    m = LINE.search(out)
    if m is None:
        raise SystemExit(f"could not parse spike output: {out!r}")
    if (int(m.group(1)) != pairs or int(m.group(2)) != read_len
            or m.group(4) != mode):
        raise SystemExit(f"spike echoed the wrong parameters: {out!r}")
    rss = ru.ru_maxrss if DARWIN else ru.ru_maxrss * 1024
    return rss, int(m.group(3))


def slope_for(binary, counts, read_len, mode):
    xs, ys, texts = [], [], []
    for n in counts:
        rss, text = run(binary, n, read_len, mode)
        if mode == "full":
            # The spike's records must carry the same text as the calibration
            # fixture, or the overhead is measured against the wrong baseline.
            expect = analytic_r(n, read_len)
            if abs(text / n / 2 - expect) > 0.05:
                raise SystemExit(f"spike text {text / n / 2:.2f} B/record != "
                                 f"fixture model {expect:.2f} at {read_len} bp")
        xs.append(n)
        ys.append(rss)
        texts.append(text)
    rss_slope, _ = lstsq(xs, ys)
    text_slope, _ = lstsq(xs, texts)
    return rss_slope, text_slope


def main():
    repo = os.path.abspath(sys.argv[1])
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    here = os.path.dirname(os.path.abspath(__file__))
    staged = os.path.join(repo, "examples", "spike_slot_cost.rs")
    binary = os.path.join(repo, "target", "release", "examples", "spike_slot_cost")

    if os.path.exists(staged):
        raise SystemExit(f"{staged} already exists — refusing to overwrite")
    shutil.copyfile(os.path.join(here, "spike_slot_cost.rs"), staged)
    try:
        build = subprocess.run(
            ["cargo", "build", "--release", "--example", "spike_slot_cost"],
            cwd=repo, capture_output=True, text=True)
        if build.returncode != 0:
            raise SystemExit("build failed:\n" + build.stderr[-3000:])
        print(f"# platform {sys.platform}, mimalloc, reps {reps}, "
              f"slope over {len(POINTS)} record counts per read length")
        print(f"# size_of: spine = {SPINE_SIZEOF:.0f} B/pair, s = {S_SIZEOF:.0f}; "
              f"matrix fit implies s = {S_FITTED:.2f}\n")
        print(f"{'len':>4s} {'spine':>7s} {'vs 148':>7s} {'full':>8s} {'2r':>7s} "
              f"{'heaptext':>9s} {'alloc':>7s} {'alloc/str':>9s} {'s':>8s}")

        table = {}
        for read_len in READ_LENGTHS:
            r_guess = 34 + 8 + 2 * read_len + 5
            n_max = int(TARGET_BYTES / (2 * r_guess + 220))
            counts = [max(1000, int(n_max * f)) for f in POINTS]
            rows = []
            for _ in range(reps):
                spine, _ = slope_for(binary, counts, read_len, "spine")
                full, text_slope = slope_for(binary, counts, read_len, "full")
                r = text_slope / 2
                heap_text = 2 * (r - 5)
                alloc = (full - spine) - heap_text
                rows.append((spine, full, text_slope, heap_text, alloc,
                             (full - text_slope) / 2))
            spine = sum(x[0] for x in rows) / reps
            full = sum(x[1] for x in rows) / reps
            text_slope = sum(x[2] for x in rows) / reps
            heap_text = sum(x[3] for x in rows) / reps
            alloc = sum(x[4] for x in rows) / reps
            s = sum(x[5] for x in rows) / reps
            table[read_len] = (spine, full, text_slope, heap_text, alloc, s,
                               max(x[5] for x in rows) - min(x[5] for x in rows))
            print(f"{read_len:>4d} {spine:>7.1f} {spine - SPINE_SIZEOF:>+7.1f} "
                  f"{full:>8.1f} {text_slope:>7.1f} {heap_text:>9.1f} "
                  f"{alloc:>7.1f} {alloc / 6:>9.1f} {s:>8.2f}", flush=True)

        print("\n=== verdict ===")
        spines = [v[0] for v in table.values()]
        print(f"  spine-only slope: {min(spines):.1f} to {max(spines):.1f} B/pair "
              f"against size_of's {SPINE_SIZEOF:.0f} -> "
              f"{'CONFIRMED' if max(abs(x - SPINE_SIZEOF) for x in spines) < 3 else 'NOT confirmed'}")
        allocs = [v[4] / 6 for v in table.values()]
        print(f"  allocator cost per String: {min(allocs):.1f} to {max(allocs):.1f} B "
              f"(6 Strings per pair)")
        ss = [v[5] for v in table.values()]
        rep_spread = max(v[6] for v in table.values())
        print(f"  s = {min(ss):.2f} to {max(ss):.2f} B/record across read lengths; "
              f"run-to-run spread {rep_spread:.2f}")
        print(f"  criterion (run-to-run <= 10 B/pair): "
              f"{'MET' if 2 * rep_spread <= 10 else 'NOT MET'}")
        print(f"  size_of s = {S_SIZEOF:.0f}: "
              f"{'consistent' if min(ss) <= S_SIZEOF <= max(ss) else 'EXCLUDED'}")
        print(f"  matrix  s = {S_FITTED:.2f}: "
              f"{'consistent' if min(ss) <= S_FITTED <= max(ss) else 'EXCLUDED'}")
        print(f"  single constant exact at every read length: "
              f"{'yes' if max(ss) - min(ss) < 5 else f'NO — s varies by {max(ss) - min(ss):.1f} B/record'}")
    finally:
        if os.path.exists(staged):
            os.remove(staged)
            print(f"\n# removed staged {staged}")


if __name__ == "__main__":
    main()
