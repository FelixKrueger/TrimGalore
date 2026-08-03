#!/usr/bin/env python3
"""Spike: WHICH serial stage caps TrimGalore's parallel path, and by how much?

ANTECEDENT
    plans/rapidgzip-trial/spikes/SPIKE_reader_bound.md established that parallel
    gzip decompression is worth only ~1.12x to us (vs 1.57x to fastqc-rust),
    because our decompression is already overlapped with parallel trimming. Its
    real finding was that BOTH gzipped and plain input plateau at 8 cores
    (~360 MB/s) against an ideal 5x from 2->10 cores (achieved 3.0x). It
    attributed that plateau to "the single reader thread as a whole".

    That attribution was made at 10 cores, where no serial stage can be
    separated from any other. This spike decomposes it.

QUESTION
    Single-end at --cores N there are THREE serial stages, not one:
      A  fastq.rs:286   MultiGzDecoder + read_next_direct + 4096-batching,
                        all fused on ONE thread.          -> per RECORD
      C  parallel.rs:947 round-robin dispatch into per-worker sync_channel(2).
                        1 send per 4096 records.          -> per BATCH
      B  parallel.rs:955 main thread: BTreeMap reorder + serial write_all.
                                                          -> per BYTE
    Which one is the ceiling, and what is actually available if it is fixed?

    Note our own --help already claims "beyond --cores 8, gzip-output I/O on the
    storage layer typically becomes binding" -- i.e. it asserts B, contradicting
    the antecedent's attribution to A. This spike settles that.

SUCCESS CRITERIA (fixed before running)
    1. Reproduce the plateau off the laptop and on real data, or show it absent.
    2. Attribute it: report each stage's share at max cores, with the residual
       stage named. A stage is "the ceiling" if removing it moves the plateau.
    3. Quantify: serial fraction S from an Amdahl fit, and the implied ceiling.
    Decision rule: if B >= A_parse at max cores, the antecedent's "profile the
    reader" recommendation is wrong and the work is writer-side instead.

CONTROLS CARRIED FORWARD (each would otherwise invalidate the result)
    1. --dont_gzip in the gz/plain comparison. Output compression follows input
       by default (#245), so without it the plain arm also skips OUTPUT
       compression -- two variables, not one.
    2. Explicit -a. Auto-detection scans up to 1M reads FROM THE INPUT, so it is
       itself sensitive to decompression speed.
    3. cores=1 reported but EXCLUDED from every fit: main.rs routes cores>1 to
       the worker pool and cores==1 to the sequential path -- different
       architectures, not two points on one curve.
    4. plain is built by decompressing the gz, never independently, so both arms
       hold byte-identical reads.

CONTROLS ADDED FOR THIS RUN
    5. CPU pinning to node0's physical cores. The host is 64 physical / 128
       logical over 2 NUMA nodes (node0 = 0-31,64-95). Unpinned, crossing 32
       cores adds remote-memory latency to the very channel handoffs we are
       measuring, and crossing 64 lands on SMT siblings worth ~1.2x not 2x.
       Either fabricates a plateau. Sweep stops at 32 = one node's physical
       cores. numactl is absent on this host and sudo needs a password, so this
       uses `taskset -c 0-31`: CPU affinity is pinned, memory is NOT explicitly
       bound. With every thread on node0 and Linux first-touch policy,
       allocations land on node0 anyway -- close to --membind=0 but not
       guaranteed under memory pressure. Documented as a deviation.
    6. Sentinel cell. (plain, cores=8) runs first and last; >5% drift means the
       machine changed under us and the sweep is void. `ps` cannot see
       co-tenants from inside this container, so loadavg is the only other
       signal and our own runs pollute it.
    7. Read-length discriminator. Real data here is 65 bp; the synthetic control
       is 150 bp. If the plateau is constant in READS/s the ceiling is
       per-record (stage A parse, or C); if constant in MB/s it is per-byte
       (stage B write, or memcpy). Those need different fixes, so both units are
       reported everywhere.
    8. Input read counts asserted after construction, keyed on the parameter
       that generates them (antecedent iteration #2 silently measured a stale
       68 MB input while reporting 300 copies).
    9. Output writeback forced inside the timed window (SPIKE_PERSIST, default
       on). Added after iteration #3: deleting each output before the kernel
       flushed it made the serial writer look free (-1% of runtime) when it is
       in fact 70% with --dont_gzip. This is the single control that decided the
       spike's answer, and the first sweep shipped without it.

NOT IN SCOPE
    Integrating rapidgzip; uBAM/BGZF; paired-end (2x stage A muddies
    attribution); --clumpify; --fastqc.
"""

import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

SCRATCH = Path(os.environ.get("SPIKE_SCRATCH", "/tmp/spike_reader_oxy"))
TG = Path(os.environ.get("SPIKE_TG", str(Path.home() / "TrimGalore/target/release/trim_galore")))
SRC_GZ = SCRATCH / "SRR24827373_1.fastq.gz"      # 84M x 65bp WGBS, ENA
SYNTH_SRC = SCRATCH / "10K_150bp.fastq.gz"       # 10k x 150bp, the antecedent's fixture

REAL_READS = int(os.environ.get("SPIKE_READS", "40000000"))
CORE_LEVELS = [1, 2, 4, 8, 16, 32]
# The synth arm needs the SAME levels as the real arms: its fitted S feeds the
# per-record/per-byte solve, and a 3-point fit is too noisy for that.
SYNTH_LEVELS = CORE_LEVELS
REPEATS_PRIMARY = int(os.environ.get("SPIKE_REPEATS", "3"))
REPEATS_SECONDARY = 2
# Default ON: see run_once. Off measures a memcpy, not a write.
PERSIST = os.environ.get("SPIKE_PERSIST", "1") != "0"
ADAPTER = "AGATCGGAAGAGC"
# node0 physical cores only: excludes SMT siblings (64-95) and all of node1.
PIN_CPUS = os.environ.get("SPIKE_PIN", "0-31")
PIN = ["taskset", "-c", PIN_CPUS]

# arm -> (input key, extra flags, what it isolates)
ARMS = {
    "gz":         ("real_gz",    ["--dont_gzip"],                     "A(+decomp) + full B"),
    "plain":      ("real_plain", ["--dont_gzip"],                     "A(-decomp) + full B"),
    "writer_off": ("real_plain", ["--dont_gzip", "--length", "10000"], "A(-decomp), B ~= 0"),
    "writer_cmp": ("real_gz",    [],                                  "A(+decomp) + B/3.5 + worker compress"),
    "synth":      ("synth",      ["--dont_gzip"],                     "150bp: read-length discriminator"),
}
PRIMARY = ("gz", "plain")


def sh(cmd: str) -> str:
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {cmd}\n{r.stderr[-2000:]}")
    return r.stdout.strip()


def loadavg() -> float:
    return float(Path("/proc/loadavg").read_text().split()[0])


def count_reads(p: Path) -> int:
    """Line count / 4. wc -l on 8 GB is ~2 s and is the only honest check."""
    return int(sh(f"wc -l < {p}")) // 4


def build_inputs() -> dict:
    """Build all inputs and assert their size. plain is derived FROM gz."""
    SCRATCH.mkdir(parents=True, exist_ok=True)
    inputs = {}

    # --- real arms: 40M reads, gz and plain holding identical records ---
    rg = SCRATCH / f"real_{REAL_READS}.fastq.gz"
    rp = SCRATCH / f"real_{REAL_READS}.fastq"
    if not rg.exists():
        # head closes the pipe early; pigz -dc's SIGPIPE exit is expected.
        sh(f"set +o pipefail; pigz -dc {SRC_GZ} | head -n {REAL_READS * 4} | pigz -p 16 > {rg}")
    if not rp.exists():
        sh(f"pigz -dc {rg} > {rp}")

    n = count_reads(rp)
    if n != REAL_READS:
        sys.exit(f"real input built wrong: {n:,} reads, expected {REAL_READS:,}. Stale cache?")
    inputs["real_gz"] = (rg, n, rp.stat().st_size)
    inputs["real_plain"] = (rp, n, rp.stat().st_size)

    # --- synthetic 150bp control, byte-matched to the real plain arm ---
    one = SCRATCH / "synth_one.fastq"
    if not one.exists():
        sh(f"pigz -dc {SYNTH_SRC} > {one}")
    per_copy = one.stat().st_size
    copies = max(1, round(rp.stat().st_size / per_copy))
    sp = SCRATCH / f"synth_{copies}.fastq"
    if not sp.exists():
        sh(f"for i in $(seq {copies}); do cat {one}; done > {sp}")
    sn = count_reads(sp)
    if sn != copies * 10_000:
        sys.exit(f"synth input built wrong: {sn:,} reads, expected {copies * 10_000:,}")
    inputs["synth"] = (sp, sn, sp.stat().st_size)

    for k, (p, reads, size) in inputs.items():
        bp = size / reads / 4  # crude: bytes per record / 4 lines
        print(f"  {k:<11} {reads:>12,} reads  {size/1e9:6.2f} GB  ~{bp:.0f} B/line  {p.name}")
    return inputs


def run_once(arm: str, cores: int, inputs: dict) -> tuple:
    """One trim run. Returns (seconds, loadavg_before).

    PERSIST=1 is mandatory for any conclusion about the serial writer. Deleting
    the output straight after the run discards its dirty pages before the kernel
    flushes them, so the write costs a memcpy instead of a disk write. The first
    full sweep ran with PERSIST=0 and measured the writer at -1% of runtime; with
    PERSIST=1 the same cell measured it at 70%. PERSIST=0 is retained only for
    isolating the compute side, and needs ~8 GB of free scratch either way.
    """
    key, extra, _ = ARMS[arm]
    inp, _, _ = inputs[key]
    outdir = SCRATCH / f"out_{arm}_{cores}"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    cmd = PIN + [str(TG), "-a", ADAPTER, *extra,
                 "--cores", str(cores), "-o", str(outdir), str(inp)]
    la = loadavg()
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if PERSIST:
        subprocess.run(["sync"], check=False)      # inside the timed window
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"FAILED ({arm}, cores={cores}):\n{' '.join(cmd)}\n{r.stderr[-2000:]}")
    shutil.rmtree(outdir, ignore_errors=True)
    return dt, la


def solve(A: list, b: list) -> list:
    """Gaussian elimination with partial pivoting. Avoids a numpy dependency."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-15:
            return []
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def fit(cores: list, times: list, basis) -> list:
    """Least squares over the given basis functions of N."""
    k = len(basis)
    A = [[sum(basis[i](n) * basis[j](n) for n in cores) for j in range(k)] for i in range(k)]
    b = [sum(basis[i](n) * t for n, t in zip(cores, times)) for i in range(k)]
    return solve(A, b)


def main() -> None:
    for p in (TG, SRC_GZ, SYNTH_SRC):
        if not p.exists():
            sys.exit(f"missing: {p}")
    if shutil.which("taskset") is None:
        sys.exit("taskset absent -- CPU pinning is control 5, do not run without it")
    # Prove the pin actually restricts before trusting any number it produces.
    seen = int(sh(f"taskset -c {PIN_CPUS} nproc"))
    lo, hi_cpu = (int(x) for x in PIN_CPUS.split("-"))
    want = hi_cpu - lo + 1
    if seen != want:
        sys.exit(f"pin not effective: taskset -c {PIN_CPUS} sees {seen} cpus, expected {want}")

    print(f"binary  : {TG}")
    print(f"version : {sh(f'{TG} --version | head -1')}")
    print(f"cpus    : {os.cpu_count()} logical; pinned to {PIN_CPUS} ({seen} cpus, verified)")
    print(f"loadavg : {Path('/proc/loadavg').read_text().strip()}")
    print(f"repeats : {REPEATS_PRIMARY} primary / {REPEATS_SECONDARY} secondary (min wall clock)\n")

    print("INPUTS")
    inputs = build_inputs()

    # Sentinel runs between every arm, not just at the ends: co-tenants roam all
    # 128 CPUs while we are pinned to 32 of them, so we are maximally exposed and
    # an 85-minute window can drift in the middle and recover by the end.
    print("\nSENTINEL (plain, cores=8) baseline")
    sent_base, _ = run_once("plain", 8, inputs)
    print(f"  {sent_base:.2f}s")
    sentinels = []

    results = {a: {} for a in ARMS}
    loads = []
    print("\nSWEEP")
    for arm in ARMS:
        levels = SYNTH_LEVELS if arm == "synth" else CORE_LEVELS
        reps = REPEATS_PRIMARY if arm in PRIMARY else REPEATS_SECONDARY
        _, _, isolates = ARMS[arm]
        print(f"\n  [{arm}] {isolates}")
        for cores in levels:
            r = 1 if cores == 1 else reps       # cores=1 is a reference point, not a fit point
            ts = []
            for _ in range(r):
                dt, la = run_once(arm, cores, inputs)
                ts.append(dt)
                loads.append((arm, cores, la))
            best = min(ts)
            results[arm][cores] = best
            _, reads, size = inputs[ARMS[arm][0]]
            spread = f"spread {(max(ts) - best) / best * 100:.1f}%" if r > 1 else "spread n/a"
            print(f"    cores={cores:<3} {best:7.2f}s  {size/1e6/best:7.1f} MB/s  "
                  f"{reads/1e6/best:6.2f} M reads/s  (n={r}, {spread})")
        st, _ = run_once("plain", 8, inputs)
        sentinels.append((arm, st))
        print(f"    sentinel after [{arm}]: {st:.2f}s")

    # ---------------- analysis ----------------
    bar = "=" * 74
    pool = [c for c in CORE_LEVELS if c > 1]

    print(f"\n{bar}\nSTAGE ATTRIBUTION at cores={pool[-1]}\n{bar}")
    hi = pool[-1]
    t_gz, t_pl = results["gz"][hi], results["plain"][hi]
    t_wo, t_wc = results["writer_off"][hi], results["writer_cmp"][hi]
    print(f"  gz          {t_gz:7.2f}s   baseline, decompress + parse + full write")
    print(f"  plain       {t_pl:7.2f}s   decompression removed")
    print(f"  writer_off  {t_wo:7.2f}s   serial write removed too")
    print(f"  writer_cmp  {t_wc:7.2f}s   write volume /3.5, worker compress added")
    print()
    print(f"  decompression (gz - plain)       {t_gz - t_pl:6.2f}s  {(t_gz-t_pl)/t_gz*100:5.1f}% of gz")
    print(f"  serial write  (plain - wr_off)   {t_pl - t_wo:6.2f}s  {(t_pl-t_wo)/t_pl*100:5.1f}% of plain")
    print(f"  residual      (writer_off)       {t_wo:6.2f}s  {t_wo/t_pl*100:5.1f}% of plain  <- A_parse + C + trim/N")

    print(f"\n{bar}\nSCALING (worker-pool path only)\n{bar}")
    print(f"{'arm':<11}" + "".join(f"{c:>9}" for c in pool) + f"{'2->' + str(hi):>10}")
    for arm in ARMS:
        levels = [c for c in (SYNTH_LEVELS if arm == "synth" else CORE_LEVELS) if c > 1]
        row = "".join(f"{results[arm].get(c, float('nan')):>9.2f}" for c in levels)
        sp = results[arm][levels[0]] / results[arm][levels[-1]]
        print(f"{arm:<11}{row}{sp:>9.2f}x")
    print(f"  ideal 2->{hi}: {hi/2:.2f}x")

    print(f"\n{bar}\nAMDAHL FIT  T(N) = S + P/N   (2-param, and +C*N for contention)\n{bar}")
    fitted = {}
    for arm in ARMS:
        levels = [c for c in (SYNTH_LEVELS if arm == "synth" else CORE_LEVELS) if c > 1]
        ts = [results[arm][c] for c in levels]
        two = fit(levels, ts, [lambda n: 1.0, lambda n: 1.0 / n])
        if not two:
            continue
        S, P = two
        fitted[arm] = S
        pred = [S + P / n for n in levels]
        resid = [(m - p) / m * 100 for m, p in zip(ts, pred)]
        _, reads, size = inputs[ARMS[arm][0]]
        ceil_r = reads / 1e6 / S if S > 0 else float("inf")
        print(f"  {arm:<11} S={S:6.2f}s  P={P:7.2f}s  serial@{levels[-1]}c="
              f"{S/ts[-1]*100:5.1f}%  ceiling={ceil_r:5.2f} M reads/s "
              f"({size/1e6/S if S > 0 else 0:.0f} MB/s)")
        print(f"  {'':<11} residuals: " + " ".join(f"{c}c:{r:+.1f}%" for c, r in zip(levels, resid)))
        if len(levels) >= 4:
            three = fit(levels, ts, [lambda n: 1.0, lambda n: 1.0 / n, lambda n: float(n)])
            if three:
                # A positive C is only meaningful if its term is a real share of
                # runtime at max cores; +0.0001 s/core is noise, not contention.
                share = three[2] * levels[-1] / ts[-1] * 100
                verdict = (f"contention {share:.1f}% of T at {levels[-1]}c"
                           if share > 2 else "no material contention term")
                print(f"  {'':<11} 3-param: S={three[0]:6.2f}s P={three[1]:7.2f}s "
                      f"C={three[2]:+.4f}s/core  ({verdict})")

    print(f"\n{bar}\nSERIAL-STAGE COST MODEL  S = a*records + b*bytes\n{bar}")
    print("  Solved from the two read lengths' fitted S. Comparing raw throughput")
    print("  would be confounded: the arms differ in per-read trim work too, and")
    print("  only S isolates the serial stage from that.")
    if "plain" in fitted and "synth" in fitted:
        _, rr, rs = inputs["real_plain"]
        _, sr, ss = inputs["synth"]
        sol = solve([[rr, rs], [sr, ss]], [fitted["plain"], fitted["synth"]])
        if sol and sol[0] > 0 and sol[1] > 0:
            a, b = sol
            print(f"\n  a = {a*1e9:7.1f} ns/record      b = {b*1e9:8.3f} ns/byte "
                  f"({1/b/1e6:.0f} MB/s per-byte ceiling)")
            for lbl, key in (("real 65bp", "real_plain"), ("synth 150bp", "synth")):
                _, n, sz = inputs[key]
                rec_t, byt_t = a * n, b * sz
                tot = rec_t + byt_t
                print(f"  {lbl:<12} per-record {rec_t:6.2f}s ({rec_t/tot*100:4.1f}%)   "
                      f"per-byte {byt_t:6.2f}s ({byt_t/tot*100:4.1f}%)   "
                      f"{sz/n:.0f} B/record")
            print("\n  => optimise the per-BYTE path" if b * (rs / rr) > a
                  else "\n  => optimise the per-RECORD path")
            print("  (2 read lengths = 2 equations, exactly determined: no residual to")
            print("   check, so treat the split as indicative. A third read length")
            print("   would make it overdetermined and testable.)")
        else:
            print("  solve failed or gave a negative coefficient -- S estimates too noisy")

    print(f"\n{bar}\nSENTINEL DRIFT (co-tenancy)\n{bar}")
    # Drift is measured against the MEDIAN of the in-sweep sentinels, not against
    # the baseline. The first run of the session competes with whatever writeback
    # the input build left pending, so it reads ~16% slow and previously voided a
    # sweep whose five in-sweep sentinels agreed to 2.2%.
    med = statistics.median(st for _, st in sentinels) if sentinels else 0.0
    print(f"  median in-sweep sentinel {med:.2f}s  (the drift reference)")
    print(f"  baseline {sent_base:.2f}s  ({(sent_base/med - 1)*100:+.1f}% vs median) "
          f"-- cold-start reference, NOT part of the drift test")
    for arm, st in sentinels:
        d = (st / med - 1) * 100
        print(f"  after {arm:<11} {st:6.2f}s  {d:+5.1f}%" + ("   !! VOID" if abs(d) > 5 else ""))
    worst = max((abs(st / med - 1) * 100 for _, st in sentinels), default=0.0)
    print(f"  worst drift {worst:.1f}%  -> "
          + ("sweep is internally consistent" if worst <= 5 else "TREAT THE SWEEP AS VOID"))
    print(f"\n  loadavg over sweep: min {min(l for _, _, l in loads):.1f}  "
          f"max {max(l for _, _, l in loads):.1f}  (our own runs inflate this)")
    print(f"\nscratch: {SCRATCH}  (delete when done)")


if __name__ == "__main__":
    main()
