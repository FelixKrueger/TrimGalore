#!/usr/bin/env python3
"""Phase 4 calibration for #457: confirm the code-derived per-record spine cost
and fit the two occupancy terms, with read length as a first-class axis.

Derived from `phase2/calibrate.py`, which is left untouched so the TSVs recorded
beside it stay reproducible by the script that produced them.

What changed from Phase 2, and why (each is a PLAN.md finding, not a preference):

* The banner regex is rebuilt with named groups against the shipped banner
  (§2.6). Phase 2's pattern predates #451/#468 and matches nothing, so every
  cell aborted. `self_check` now runs the real pattern *and* the rotted one
  against real captured output, so a future banner change fails once, up front,
  and the negative control proves the abort path still has teeth.

* `r` (bytes per record, as `clump::estimated_record_bytes` computes it) is
  measured per fixture and cross-checked against the analytic value for that
  pair count. Phase 4 v1 measured `r` on 2 000-pair fixtures and was 3.4 B low
  at every read length, because `mk_clump_fixture` writes the record index into
  the header (§2.4a).

* Read length is an axis, and cells are equalised on **total bytes**, not pair
  count. Holding pairs constant gives short reads far fewer flush cycles, which
  depresses their peak and biases the per-record term downward — the unsafe
  direction, on the axis being fitted (§4 Fault 1).

* The fit is `RSS = static + B·(a·n_bins + b·cores + s·(n_bins + 2·cores)/r)`.
  `s` is not being discovered: `size_of::<FastqRecord>()` = 72, `MinimizerKey` =
  4, so a paired slot is 148 B and `s = 74` B per record-byte of `r`. The fit
  confirms it and fits `a`, `b` (§2.3).

* `static` is **pinned**, not fitted. At a single budget the columns `1`,
  `B·n_bins` and `B·cores` all depend on `cores` alone, so they span a
  2-dimensional space and the three parameters are not separable. `--fit-static`
  plus more than one `--budgets` entry breaks the degeneracy; the rank report
  prints either way so the deficiency cannot hide.

usage:
  calibrate.py <binary> <workdir> [options]

options:
  --smoke                  tiny fixtures, 1 rep, 2 cells — proves the harness,
                           measures nothing
  --budgets 1024[,2048]    fit-cell budgets in MiB
  --reps N                 reps per cell (default 3)
  --static-mib X           pinned fixed term (default: platform measurement)
  --fit-static             fit the fixed term instead of pinning it
  --distinct exact|none|N  distinct-sequence scan: full, skipped, or first N
  --mk-fixture PATH        mk_clump_fixture binary (default: alongside <binary>)
"""

import os
import re
import shutil
import subprocess
import sys
import time

MIB = 1024 * 1024

# The shipped banner, verified against the binary before this pattern was
# written: "clumpify: 16 bins × 22 MiB; predicted peak ≈ 930 of 1024 MiB budget
# (gzip level 6)".
BANNER = re.compile(
    r"clumpify: (?P<bins>\d+) bins × (?P<bin_mib>\d+) MiB; "
    r"predicted peak ≈ (?P<pred_mib>\d+) of (?P<budget_mib>\d+) MiB budget "
    r"\((?P<encoding>[^)]*)\)"
)

# Phase 2's pattern, kept as the negative control: it is the one that rotted, so
# it is the one worth proving we would now notice.
ROTTED = re.compile(r"clumpify: (\d+) bins × (\d+) MiB; predicted peak ≈ (\d+) MiB")

# `prediction = budget × MARGIN_NUM/MARGIN_DEN` by construction (clump.rs).
MARGIN_NUM, MARGIN_DEN = 10, 11

# size_of::<FastqRecord>() + size_of::<MinimizerKey>() for a paired slot, halved
# to per-record-byte-of-r. Pinned in clump.rs, not fitted here.
S_EXPECTED = 74.0
A_EXPECTED, B_EXPECTED = 0.937, 2.999

# Measured fixed term on 10 M-pair paired runs (clump.rs:248-251).
STATIC_MIB = {"darwin": 205.0, "linux": 173.0}

# What the binary *charges*, not what it measures — needed to reproduce
# `resolve_layout`'s exact bin budget. clump.rs: STATIC_TRIM_BYTES, MARGIN.
STATIC_CHARGED_MIB = 224
SIGMA16, K16 = 33, 92

# ~3.0 GB per mate, so 25 bp gets ~31 M pairs and 150 bp ~8.7 M (§4).
TARGET_BYTES_PER_MATE = 3_000_000_000

READ_LENGTHS = (25, 36, 51, 100, 150)
FIT_CORES = (2, 4)
FALSIFY = [(rl, cores, 2048) for rl in (25, 100) for cores in (8, 16)]

# The fixture's header is "@SIM:1:FCX:1:15:6329:{i} {mate}:N:0:ATCCGA": 21
# leading chars, then the decimal index, then 13 trailing.
HEADER_FIXED = 34


class BannerError(RuntimeError):
    pass


def n_bins_for(cores):
    return max(16, 4 * cores)


def exact_bin_mib(budget_mib, cores):
    """`resolve_layout`'s bin budget, in MiB, with its integer arithmetic.

    The banner prints this floored to whole MiB. Regressing on the *banner*
    value is a ~5% error in the regressor at the smallest bin (c16), which is
    what phase 2 avoided by recording both and phase 4 v1 reintroduced by
    recording only the banner.
    """
    budgeted = budget_mib * MIB * 10 // 11
    usable = budgeted - STATIC_CHARGED_MIB * MIB
    denom = SIGMA16 * n_bins_for(cores) + K16 * cores
    return (usable * 16 // denom) / MIB


def mean_index_digits(pairs):
    """Mean decimal width of 0..pairs-1 — exact, not approximated."""
    total, width, lo = 0, 1, 0
    while lo < pairs:
        hi = min(pairs, 10**width)
        total += (hi - lo) * width
        lo, width = hi, width + 1
    return total / pairs


def analytic_r(pairs, read_len):
    """`estimated_record_bytes` for this fixture: id + seq + qual + 5."""
    header = HEADER_FIXED + mean_index_digits(pairs)
    return header + 2 * read_len + 5


def solve_pairs(target_bytes, read_len):
    """Pairs giving `target_bytes` per mate. `r` depends on the pair count
    through the header, so this is a fixed point rather than a division."""
    pairs = max(1, round(target_bytes / analytic_r(1, read_len)))
    for _ in range(8):
        nxt = max(1, round(target_bytes / analytic_r(pairs, read_len)))
        if nxt == pairs:
            break
        pairs = nxt
    return pairs


def parse_banner(text, cores, budget_mib, log_path, pattern=BANNER):
    """Parse and *validate*. A pattern that matches but yields the wrong numbers
    is the failure Phase 2 could not distinguish from a healthy parse."""
    m = pattern.search(text)
    if m is None:
        raise BannerError(f"banner did not match {pattern.pattern!r} — see {log_path}")
    try:
        g = m.groupdict()
    except Exception:  # a positional-group pattern, i.e. the rotted one
        raise BannerError(f"pattern has no named groups — see {log_path}")
    for key in ("bins", "bin_mib", "pred_mib", "budget_mib"):
        if key not in g or g[key] is None or g[key] == "":
            raise BannerError(f"group {key!r} did not capture — see {log_path}")
    bins, bin_mib = int(g["bins"]), int(g["bin_mib"])
    pred_mib, banner_budget = int(g["pred_mib"]), int(g["budget_mib"])
    if bins != n_bins_for(cores):
        raise BannerError(f"bins {bins} != max(16, 4*{cores}) — see {log_path}")
    if banner_budget != budget_mib:
        raise BannerError(
            f"banner budget {banner_budget} != requested {budget_mib} — see {log_path}"
        )
    if bin_mib <= 0 or pred_mib <= 0:
        raise BannerError(f"non-positive bin/pred: {bin_mib}/{pred_mib} — see {log_path}")
    expect_pred = budget_mib * MARGIN_NUM // MARGIN_DEN
    if abs(pred_mib - expect_pred) > 1:
        raise BannerError(
            f"predicted {pred_mib} MiB is not {MARGIN_NUM}/{MARGIN_DEN} of "
            f"{budget_mib} MiB (expected ~{expect_pred}) — see {log_path}"
        )
    if not g["encoding"]:
        raise BannerError(f"encoding group is empty — see {log_path}")
    return {"bins": bins, "bin_mib": bin_mib, "pred": pred_mib,
            "encoding": g["encoding"]}


def self_check(text, cores, budget_mib, log_path):
    """Prove the parser works and that its failure path fires, both against
    output the shipped binary actually produced."""
    good = parse_banner(text, cores, budget_mib, log_path)
    try:
        parse_banner(text, cores, budget_mib, log_path, pattern=ROTTED)
    except BannerError:
        pass
    else:
        raise SystemExit(
            "self-check failed: Phase 2's rotted pattern still parses the current "
            "banner, so this harness cannot detect a banner change"
        )
    mutated = re.compile(BANNER.pattern.replace("clumpify", "clumpifx"))
    try:
        parse_banner(text, cores, budget_mib, log_path, pattern=mutated)
    except BannerError:
        pass
    else:
        raise SystemExit("self-check failed: a deliberately broken pattern matched")
    print(f"# self-check    banner parses ({good['bins']} bins, {good['bin_mib']} MiB, "
          f"pred {good['pred']} MiB, {good['encoding']}); "
          f"rotted and mutated patterns both rejected")
    return good


def free_mib(path):
    return shutil.disk_usage(path).free / MIB


def require_disk(path, need_mib, what):
    have = free_mib(path)
    if have < need_mib:
        raise SystemExit(
            f"{what}: need {need_mib:.0f} MiB free under {path}, have {have:.0f}. "
            f"An out-of-space run reports inflated peak RSS, so this is fatal, "
            f"not a warning."
        )


def fixture_stats(path, distinct_mode):
    """Reads, distinct count, mean read length and measured `r`.

    Counted in a shell pipeline: on Linux a child's `ru_maxrss` is inherited
    from the forking parent, so holding sequences in this process would report
    the harness's own footprint as every measured peak.

    `r` is `len(id) + len(seq) + len(qual) + 5`, matching
    `clump::estimated_record_bytes` (the '+' line is 2 B there).
    """
    # The sequence is line 2 of the block, so it is held until the block closes;
    # printing $0 at the qual line would count distinct *qualities*, which are
    # uniform in this fixture and would read as "1 distinct sequence".
    if distinct_mode == "none":
        emit = ""
    elif distinct_mode.isdigit():
        # Limited inside awk: `head` closing the pipe early SIGPIPEs awk before
        # its END block, and the stats line is what END prints.
        emit = f"if (n <= {distinct_mode}) print sq;"
    else:
        emit = "print sq;"
    awk = (
        "awk 'NR%4==1{h=length($0)} NR%4==2{sq=$0; s=length($0)} "
        f"NR%4==0{{n++; ts+=s; tr+=h+s+length($0)+5; {emit}}} "
        "END{print \"__STATS__\", n, ts, tr > \"/dev/stderr\"}'"
    )
    if distinct_mode == "none":
        pipeline = f"gzip -dc {path!r} | {awk} >/dev/null; echo 0"
    else:
        pipeline = (
            f"gzip -dc {path!r} | {awk} | "
            f"sort -u -T {os.path.dirname(path)!r} | wc -l"
        )
    out = subprocess.run(["sh", "-c", pipeline], capture_output=True, text=True,
                         check=True)
    distinct = int(out.stdout.strip())
    stats = [ln for ln in out.stderr.strip().splitlines() if "__STATS__" in ln]
    if not stats:
        raise SystemExit(f"fixture scan produced no stats line for {path}")
    _, reads, total_seq, total_rec = stats[-1].split()
    reads = int(reads)
    if reads == 0:
        raise SystemExit(f"fixture scan counted 0 records in {path}")
    return reads, distinct, int(total_seq) / reads, int(total_rec) / reads


def parent_footprint_mib():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
        MIB if sys.platform == "darwin" else 1024)


def make_fixture(mk, pairs, read_len, prefix):
    p = subprocess.run([mk, str(pairs), str(read_len), prefix],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"mk_clump_fixture failed: {p.stderr.strip()}")
    return f"{prefix}_R1.fastq.gz", f"{prefix}_R2.fastq.gz"


def run_cell(binary, r1, r2, cores, budget_mib, rep, workdir, checked):
    outdir = os.path.join(workdir, f"out_c{cores}_{budget_mib}_{rep}")
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    log = os.path.join(workdir, f"log_c{cores}_{budget_mib}_{rep}.txt")
    argv = [binary, "--paired", "--cores", str(cores), "--clumpify",
            "--compression", "6", "--memory", f"{budget_mib}M", "-o", outdir,
            r1, r2]
    with open(log, "wb") as lf:
        t0 = time.monotonic()
        p = subprocess.Popen(argv, stdout=lf, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.monotonic() - t0
    code = os.waitstatus_to_exitcode(status)
    text = open(log, encoding="utf-8", errors="replace").read()
    if code != 0:
        raise SystemExit(f"cores={cores} budget={budget_mib}M rep={rep}: "
                         f"exit={code} — see {log}")
    if not checked:
        banner = self_check(text, cores, budget_mib, log)
    else:
        banner = parse_banner(text, cores, budget_mib, log)
    rss_mib = ru.ru_maxrss / (MIB if sys.platform == "darwin" else 1024)
    out_bytes = sum(os.path.getsize(os.path.join(outdir, f))
                    for f in os.listdir(outdir)
                    if f.endswith(".fq.gz") or f.endswith(".fq"))
    shutil.rmtree(outdir, ignore_errors=True)
    return {"cores": cores, "budget": budget_mib, "bins": banner["bins"],
            "b_banner": banner["bin_mib"],
            "b_exact": exact_bin_mib(budget_mib, cores), "pred": banner["pred"],
            "rss": rss_mib, "ratio": rss_mib / banner["pred"], "wall": wall,
            "cpu": ru.ru_utime + ru.ru_stime, "out": out_bytes}


def solve(a, b):
    """Gaussian elimination with partial pivoting. Returns (x, singular)."""
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            return None, True
        m[col], m[piv] = m[piv], m[col]
        for r in range(n):
            if r == col:
                continue
            f = m[r][col] / m[col][col]
            for c in range(col, n + 1):
                m[r][c] -= f * m[col][c]
    return [m[i][n] / m[i][i] for i in range(n)], False


def least_squares(rows, targets, labels):
    """Normal equations, with the rank check printed rather than assumed."""
    k = len(labels)
    ata = [[sum(r[i] * r[j] for r in rows) for j in range(k)] for i in range(k)]
    atb = [sum(r[i] * t for r, t in zip(rows, targets)) for i in range(k)]
    x, singular = solve(ata, atb)
    if singular:
        raise SystemExit(
            "design matrix is rank deficient for " + ", ".join(labels) +
            " — at one budget the columns 1, B*n_bins and B*cores depend on "
            "cores alone. Pass more than one --budgets value, or pin the fixed "
            "term (drop --fit-static)."
        )
    return x


def main():
    argv = sys.argv[1:]
    opts = {"smoke": False, "budgets": [1024], "reps": 3, "distinct": "exact",
            "static": None, "fit_static": False, "mk": None}
    pos = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--smoke":
            opts["smoke"] = True
        elif a == "--fit-static":
            opts["fit_static"] = True
        elif a == "--budgets":
            i += 1; opts["budgets"] = [int(x) for x in argv[i].split(",")]
        elif a == "--reps":
            i += 1; opts["reps"] = int(argv[i])
        elif a == "--distinct":
            i += 1; opts["distinct"] = argv[i]
        elif a == "--static-mib":
            i += 1; opts["static"] = float(argv[i])
        elif a == "--mk-fixture":
            i += 1; opts["mk"] = argv[i]
        elif a.startswith("--"):
            raise SystemExit(f"unknown option {a}")
        else:
            pos.append(a)
        i += 1
    if len(pos) != 2:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    binary, workdir = pos
    binary = os.path.abspath(binary)
    os.makedirs(workdir, exist_ok=True)

    mk = opts["mk"] or os.path.join(os.path.dirname(binary), "examples",
                                    "mk_clump_fixture")
    if not os.access(mk, os.X_OK):
        raise SystemExit(f"mk_clump_fixture not executable at {mk} — build it with "
                         f"`cargo build --release --example mk_clump_fixture` "
                         f"or pass --mk-fixture")

    read_lengths = (25, 51) if opts["smoke"] else READ_LENGTHS
    reps = 1 if opts["smoke"] else opts["reps"]
    target_bytes = 4_000_000 if opts["smoke"] else TARGET_BYTES_PER_MATE
    falsify = [] if opts["smoke"] else FALSIFY
    static_pin = opts["static"]
    if static_pin is None:
        static_pin = STATIC_MIB.get(sys.platform, STATIC_MIB["linux"])

    footprint = parent_footprint_mib()
    print(f"# platform      {sys.platform}")
    print(f"# binary        {binary}")
    print(f"# mode          {'SMOKE — proves the harness, measures nothing' if opts['smoke'] else 'measurement'}")
    print(f"# budgets       {opts['budgets']} MiB, reps/cell {reps}, "
          f"cell summary = MAX over reps (one-sided promise)")
    print(f"# bytes/mate    {target_bytes:,} — cells equalised on total bytes, "
          f"not pair count")
    print(f"# fixed term    {static_pin:.1f} MiB "
          f"({'FITTED, pin ignored' if opts['fit_static'] else 'pinned'})")
    print(f"# harness rss   {footprint:.1f} MiB")
    if footprint > 200:
        raise SystemExit(f"harness footprint {footprint:.0f} MiB is too large to "
                         f"measure children reliably")
    if opts["distinct"] == "none":
        print("# distinct      NOT MEASURED (--distinct none)")
    elif opts["distinct"].isdigit():
        print(f"# distinct      SAMPLE of first {opts['distinct']} reads only")
    print()

    tsv = os.path.join(workdir, f"results_phase4_{sys.platform}.tsv")
    rows, checked = [], False
    with open(tsv, "w") as fh:
        fh.write("role\tread_len\tpairs\tr_measured\tr_analytic\tdistinct\t"
                 "cores\tn_bins\tbudget_mib\tb_exact_mib\tb_banner_mib\t"
                 "predicted_mib\trep\trss_mib\tratio\twall_s\tcpu_s\tout_bytes\n")
        print(f"{'role':8s} {'len':>4s} {'pairs':>10s} {'r':>7s} {'cores':>5s} "
              f"{'bud':>6s} {'B':>5s} {'pred':>6s} {'maxRSS':>9s} {'ratio':>6s} "
              f"{'wall':>7s}")
        for read_len in read_lengths:
            pairs = solve_pairs(target_bytes, read_len)
            r_analytic = analytic_r(pairs, read_len)
            # Two mates of compressed fixture plus one run's worth of output.
            require_disk(workdir, 3 * pairs * r_analytic / MIB,
                         f"fixture at {read_len} bp")
            prefix = os.path.join(workdir, f"fx_{read_len}")
            r1, r2 = make_fixture(mk, pairs, read_len, prefix)
            reads, distinct, mean_len, r_measured = fixture_stats(r1, opts["distinct"])
            if reads != pairs:
                raise SystemExit(f"fixture at {read_len} bp holds {reads} records, "
                                 f"asked for {pairs}")
            if abs(r_measured - r_analytic) > 0.05:
                raise SystemExit(
                    f"measured r {r_measured:.2f} != analytic {r_analytic:.2f} at "
                    f"{read_len} bp — the header model is wrong, and r is the axis "
                    f"being fitted"
                )
            print(f"# fixture {read_len:3d} bp: {pairs:,} pairs, r = {r_measured:.2f} B "
                  f"(analytic {r_analytic:.2f}), mean len {mean_len:.0f}, "
                  f"distinct {distinct:,}", flush=True)

            cells = [(c, bud, "fit") for bud in opts["budgets"] for c in FIT_CORES]
            cells += [(c, bud, "falsify") for rl, c, bud in falsify if rl == read_len]
            for cores, budget, role in cells:
                reps_out = []
                for rep in range(1, reps + 1):
                    require_disk(workdir, 2 * pairs * r_analytic / MIB,
                                 f"cell {read_len}bp/c{cores}/{budget}M")
                    res = run_cell(binary, r1, r2, cores, budget, rep, workdir,
                                   checked)
                    checked = True
                    reps_out.append(res)
                    fh.write("\t".join(str(x) for x in [
                        role, read_len, pairs, f"{r_measured:.3f}",
                        f"{r_analytic:.3f}", distinct, cores, res["bins"], budget,
                        f"{res['b_exact']:.4f}", res["b_banner"], res["pred"], rep,
                        f"{res['rss']:.1f}",
                        f"{res['ratio']:.4f}", f"{res['wall']:.2f}",
                        f"{res['cpu']:.2f}", res["out"]]) + "\n")
                    fh.flush()
                worst = max(reps_out, key=lambda x: x["rss"])
                rows.append({**worst, "role": role, "read_len": read_len,
                             "r": r_measured, "pairs": pairs,
                             "spread": max(x["rss"] for x in reps_out)
                                       - min(x["rss"] for x in reps_out)})
                print(f"{role:8s} {read_len:>4d} {pairs:>10,d} {r_measured:7.2f} "
                      f"{cores:>5d} {budget:>5d}M {worst['b_banner']:>4d}M "
                      f"{worst['pred']:>5d}M {worst['rss']:8.1f}M "
                      f"{worst['ratio']:6.3f} {worst['wall']:6.2f}s", flush=True)
            for f in (r1, r2):
                os.remove(f)

    report(rows, static_pin, opts["fit_static"], tsv)


def design_row(r, fit_static):
    """[1,] B*n_bins, B*cores, B*(n_bins + 2*cores)/r — the §2.3 decomposition.

    `b_exact`, never `b_banner`: the banner floors to whole MiB.
    """
    b, n, c = r["b_exact"], r["bins"], r["cores"]
    cols = [b * n, b * c, b * (n + 2 * c) / r["r"]]
    return ([1.0] + cols) if fit_static else cols


def report(rows, static_pin, fit_static, tsv):
    fit_rows = [r for r in rows if r["role"] == "fit"]
    if len(fit_rows) < (4 if fit_static else 3):
        print(f"\nwrote {tsv}\n(too few fit cells to fit — smoke run?)")
        return

    labels = (["static", "a", "b", "s"] if fit_static else ["a", "b", "s"])
    design = [design_row(r, fit_static) for r in fit_rows]
    targets = [r["rss"] - (0.0 if fit_static else static_pin) for r in fit_rows]

    # Rank is printed, not assumed: at one budget the first three columns depend
    # on cores alone and cannot separate.
    distinct_cores_b = {(r["cores"], r["b_exact"]) for r in fit_rows}
    print(f"\n=== design ===")
    print(f"  {len(fit_rows)} fit cells, {len(distinct_cores_b)} distinct "
          f"(cores, B) combinations, {len(labels)} parameters")
    if fit_static and len(distinct_cores_b) < 3:
        print("  WARNING: fewer than 3 distinct (cores, B) combinations — "
              "'static', 'a' and 'b' are not separable")

    x = least_squares(design, targets, labels)
    got = dict(zip(labels, x))
    static_fit = got.get("static", static_pin)

    print("\n=== fitted coefficients (upper-envelope cells) ===")
    print(f"  a (per bin, occupancy)  = {got['a']:6.3f}   expected {A_EXPECTED:.3f}, "
          f"and < 1 excluded from below")
    print(f"  b (per worker)          = {got['b']:6.3f}   expected {B_EXPECTED:.3f}")
    print(f"  s (B per record-byte)   = {got['s']:6.2f}   size_of-derived "
          f"{S_EXPECTED:.0f} ({(got['s'] - S_EXPECTED) / S_EXPECTED * 100:+.1f}%)")
    print(f"  static                  = {static_fit:6.1f} MiB "
          f"({'fitted' if fit_static else 'pinned'})")
    if got["a"] >= 1.0:
        print("  NOTE: a >= 1 contradicts §2.4 — the spine is resident in full "
              "from the first push, so text occupancy should sit below one copy.")

    print("\n=== residuals by read length (curvature in 1/r shows as a U) ===")
    resid = {}
    for r, row in zip(fit_rows, design):
        pred = sum(c * v for c, v in zip(row, x)) + (0.0 if fit_static else static_pin)
        resid.setdefault(r["read_len"], []).append(r["rss"] - pred)
    for rl in sorted(resid):
        v = resid[rl]
        print(f"  {rl:3d} bp: mean {sum(v) / len(v):+7.1f} MiB over {len(v)} cells "
              f"(worst {max(v, key=abs):+7.1f})")

    print("\n=== zero-parameter check: s from size_of, a and b from the plan ===")
    for r in fit_rows:
        n, c, b = r["bins"], r["cores"], r["b_exact"]
        pred = static_pin + b * (A_EXPECTED * n + B_EXPECTED * c
                                 + S_EXPECTED * (n + 2 * c) / r["r"])
        print(f"  {r['read_len']:3d} bp c{r['cores']:<2d}: model {pred:7.1f}M vs "
              f"measured {r['rss']:7.1f}M ({(r['rss'] - pred) / pred * 100:+5.1f}%)")

    falsify_rows = [r for r in rows if r["role"] == "falsify"]
    if falsify_rows:
        print("\n=== falsification: predicted from the fit, not fitted to ===")
        for r in falsify_rows:
            n, c, b = r["bins"], r["cores"], r["b_exact"]
            pred = static_fit + b * (got["a"] * n + got["b"] * c
                                     + got["s"] * (n + 2 * c) / r["r"])
            print(f"  {r['read_len']:3d} bp c{r['cores']:<2d} B={b}M: model "
                  f"{pred:7.1f}M vs measured {r['rss']:7.1f}M "
                  f"({(r['rss'] - pred) / pred * 100:+5.1f}%)")

    print("\n=== margin required ===")
    worst = max(rows, key=lambda r: r["ratio"])
    print(f"  worst ratio over all cells: {worst['ratio']:.3f} "
          f"({worst['read_len']} bp, cores={worst['cores']})")
    print(f"  largest within-cell spread: {max(r['spread'] for r in rows):.1f} MiB")
    over = [r for r in rows if r["ratio"] > 1.0]
    print(f"  cells over the printed prediction: {len(over)}"
          + ("".join(f"\n    {r['read_len']} bp c{r['cores']} {r['budget']}M: "
                     f"{r['ratio']:.3f}" for r in over) if over else ""))
    print(f"\nwrote {tsv}")


if __name__ == "__main__":
    main()
