#!/usr/bin/env python3
"""Spike: do the two one-line structural levers make #439's promise satisfiable?

Question
--------
At `--memory 1G --cores 4` the shipped binary peaks at ~1050 MiB against a
1024 MiB promise, so no coefficient change can satisfy "actual <= predicted <=
budget" with bin size held. Two candidate levers were proposed by review:

  A  work-channel depth 2 -> 1   (parallel.rs:119, :943) — fewer in-flight bins
  B  output reserve sized from actual record bytes instead of len * 300
     (parallel.rs:349, :1069) — the buffers hold *compressed* data, so
     300 B/record over-reserves ~2x on 64 bp reads, twice per pair

Success criteria
----------------
  1. At 1G/c4, ratio (peak RSS / banner prediction, both MiB) <= 1.0 with at
     least 25 MiB of headroom.
  2. At 2G/c4 (worst measured breach, 1.21x) ratio <= 1.0.
  3. Output bytes identical to baseline — these levers must not change results.
  4. No wall-clock regression beyond 3% (lever A blocks the reader sooner).

Out of scope
------------
Model fitting, coefficient derivation, Linux, FastQC, single-end, --clump_only,
real (non-cycled) libraries. Throwaway: source is patched, built, and reverted.

Run
---
  python3 plans/08162026_clumpify-memory-accounting/spikes/spike_memory-levers.py
"""

import os
import shutil
import subprocess
import sys
import time

REPO = os.environ.get("TG_REPO", "/Users/fkrueger/Github/TrimGalore")
SCRATCH = os.environ.get(
    "TG_SCRATCH",
    "/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/"
    "720d3ced-ec4c-473c-9e82-e2807221f1a7/scratchpad/rss",
)
WORK = os.path.join(SCRATCH, "spike")
R1 = os.path.join(SCRATCH, "data", "div_R1.fastq.gz")
R2 = os.path.join(SCRATCH, "data", "div_R2.fastq.gz")
TARGET = os.path.join(REPO, "src", "parallel.rs")
MIB = 1024 * 1024

# Exact source substitutions. Each variant is applied to a pristine tree.
LEVER_A = [
    ("mpsc::sync_channel::<PairedWork>(2)", "mpsc::sync_channel::<PairedWork>(1)"),
    ("mpsc::sync_channel::<SingleWork>(2)", "mpsc::sync_channel::<SingleWork>(1)"),
]
LEVER_B = [
    ("let cap = reads_r1.len() * 300;",
     "let cap = reads_r1.iter().map(estimated_record_bytes).sum::<usize>();"),
    ("let cap = reads.len() * 300;",
     "let cap = reads.iter().map(estimated_record_bytes).sum::<usize>();"),
]
VARIANTS = {
    "base": [],
    "A_chan1": LEVER_A,
    "B_reserve": LEVER_B,
    "AB_both": LEVER_A + LEVER_B,
}

CELLS = [("c4_1G", 4, "1025M"), ("c4_2G", 4, "2049M"), ("c8_1G", 8, "1025M")]
REPS = 2


def sh(cmd, **kw):
    return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, **kw)


def restore():
    sh(["git", "checkout", "--", "src/parallel.rs"])


def build_variant(name, subs):
    """Patch, build, stash the binary, always restore the source."""
    binary = os.path.join(WORK, f"tg_{name}")
    try:
        if subs:
            src = open(TARGET).read()
            for old, new in subs:
                if old not in src:
                    raise SystemExit(f"{name}: pattern absent, refusing to guess: {old!r}")
                src = src.replace(old, new)
            open(TARGET, "w").write(src)
        r = sh(["cargo", "build", "--release"])
        if r.returncode != 0:
            raise SystemExit(f"{name}: build failed\n{r.stderr[-3000:]}")
        shutil.copyfile(os.path.join(REPO, "target", "release", "trim_galore"), binary)
        os.chmod(binary, 0o755)
    finally:
        restore()
    return binary


def run_cell(binary, label, cores, mem, rep, keep=False):
    outdir = os.path.join(WORK, f"out_{label}_{os.path.basename(binary)}_{rep}")
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    log = os.path.join(WORK, f"log_{label}_{os.path.basename(binary)}_{rep}.txt")
    argv = [binary, "--paired", "--cores", str(cores), "--clumpify",
            "--compression", "6", "--memory", mem, "-o", outdir, R1, R2]
    with open(log, "wb") as lf:
        t0 = time.monotonic()
        p = subprocess.Popen(argv, stdout=lf, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.monotonic() - t0
    code = os.waitstatus_to_exitcode(status)
    text = open(log, encoding="utf-8", errors="replace").read()
    if code != 0 or "clumpify:" not in text:
        raise SystemExit(f"{label}/{binary}: exit={code}, banner "
                         f"{'present' if 'clumpify:' in text else 'MISSING'} — see {log}")
    pred = int(text.split("predicted peak ≈ ")[1].split(" MB")[0])
    out_files = sorted(f for f in os.listdir(outdir) if f.endswith(".fq.gz"))
    out_bytes = sum(os.path.getsize(os.path.join(outdir, f)) for f in out_files)
    res = {"rss": ru.ru_maxrss / MIB, "pred": pred, "wall": wall,
           "cpu": ru.ru_utime + ru.ru_stime, "out": out_bytes, "dir": outdir}
    if not keep:
        shutil.rmtree(outdir, ignore_errors=True)
    return res


def main():
    for p in (R1, R2):
        if not os.path.exists(p):
            sys.exit(f"missing fixture: {p} (regenerate the diverse pair first)")
    os.makedirs(WORK, exist_ok=True)
    if sh(["git", "status", "--porcelain", "--untracked-files=no"]).stdout.strip():
        sys.exit("tracked files are dirty — refusing to patch")

    print("building variants (source is reverted after each)...", flush=True)
    binaries = {n: build_variant(n, s) for n, s in VARIANTS.items()}
    for n, b in binaries.items():
        print(f"  {n:10s} {os.path.getsize(b):,} bytes")

    print(f"\n{'cell':8s} {'variant':10s} {'rep':>3s} {'pred':>7s} {'rss':>9s} "
          f"{'ratio':>6s} {'wall':>7s} {'cpu':>7s} {'out MB':>8s}")
    results = {}
    for label, cores, mem in CELLS:
        for name, binary in binaries.items():
            for rep in range(1, REPS + 1):
                r = run_cell(binary, label, cores, mem, rep)
                results.setdefault((label, name), []).append(r)
                print(f"{label:8s} {name:10s} {rep:>3d} {r['pred']:6d}M "
                      f"{r['rss']:8.1f}M {r['rss']/r['pred']:6.3f} {r['wall']:6.2f}s "
                      f"{r['cpu']:6.2f}s {r['out']/1e6:7.1f}", flush=True)

    print("\n--- worst-case (max over reps) ratio per cell ---")
    print(f"{'cell':8s} " + "".join(f"{n:>12s}" for n in VARIANTS))
    for label, _, _ in CELLS:
        row = f"{label:8s} "
        for name in VARIANTS:
            rs = results[(label, name)]
            row += f"{max(x['rss']/x['pred'] for x in rs):12.3f}"
        print(row)

    print("\n--- deltas vs base (max-over-reps RSS, MiB) ---")
    for label, _, _ in CELLS:
        b = max(x["rss"] for x in results[(label, "base")])
        for name in VARIANTS:
            if name == "base":
                continue
            v = max(x["rss"] for x in results[(label, name)])
            print(f"{label:8s} {name:10s} {v - b:+8.1f} MiB ({(v-b)/b*100:+5.1f}%)")

    print("\n--- criterion 4: wall vs base (mean) ---")
    for label, _, _ in CELLS:
        b = sum(x["wall"] for x in results[(label, "base")]) / REPS
        for name in VARIANTS:
            if name == "base":
                continue
            v = sum(x["wall"] for x in results[(label, name)]) / REPS
            print(f"{label:8s} {name:10s} {(v-b)/b*100:+6.2f}%")

    print("\n--- criterion 3: output byte-identity at c4_1G ---")
    keep = {n: run_cell(b, "id", 4, "1025M", 9, keep=True) for n, b in binaries.items()}
    base_dir = keep["base"]["dir"]
    for name, r in keep.items():
        if name == "base":
            continue
        same = True
        for f in sorted(os.listdir(base_dir)):
            if not f.endswith(".fq.gz"):
                continue
            a = os.path.join(base_dir, f)
            c = subprocess.run(["cmp", "-s", a, os.path.join(r["dir"], f)])
            same &= c.returncode == 0
        print(f"  {name:10s} identical={same}")
    for r in keep.values():
        shutil.rmtree(r["dir"], ignore_errors=True)

    dirty = sh(["git", "status", "--porcelain", "--untracked-files=no"]).stdout.strip()
    print(f"\nsource restored: {'YES' if not dirty else 'NO — ' + dirty}")


if __name__ == "__main__":
    main()
