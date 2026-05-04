# Buckberry-scale benchmark artifacts (2026-04-29 / 2026-04-30)

Raw `hyperfine` outputs and run logs from the Trim Galore Buckberry-scale benchmark — the data behind the [Benchmarks page](../../src/content/docs/performance/benchmarks.md).

## Test conditions

84M PE reads, [`SRR24827373`](https://www.ncbi.nlm.nih.gov/sra/?term=SRR24827373) (Buckberry et al. 2023, whole-genome bisulfite sequencing). `hyperfine --warmup 1 --runs 10` per condition. Hardware: Intel Xeon 6975P-C (Granite Rapids), 128 cores, 991 GiB RAM, dockyard-oxy-0 (Altos Labs Dockyard EKS).

## Layout

| Pattern | Description |
|---|---|
| `<engine>_c<N>.json` | Hyperfine raw output: 10 sample wall-time runs, mean/min/max/stddev, user + system CPU time. Canonical reference data. |
| `<engine>_c<N>.md` | Hyperfine human-readable markdown summary table. Same data, easier to skim. |
| `byte_identity_summary.txt` | md5 cross-check between Rust beta.5 and beta.7 outputs at cores=8. Both R1 and R2 confirmed `MATCH` — Myers' prefilter is byte-identity-preserving by construction. |
| `memory_summary.txt` | Peak RSS measurements via `/usr/bin/time -v` for v2.1.0-beta.7 at cores 1/2/4/8 (single run per condition). Wall times here match the hyperfine numbers in the corresponding `*.json` files within ~1%. |
| `logs-and-scripts.tar.gz` | Verbose stdout logs from the main run (`bench_20260429_155450.log`) + extras run (`bench_extras_20260430_120230.log`) + the three wrapper scripts (`run_bench.sh`, `run_bench_extra.sh`, `watch_then_extras.sh`). Bundled to keep the directory listing readable. |

## Engines

- `perl-0.6.11`: Perl Trim Galore 0.6.11 + Cutadapt 5.2 + igzip + pigz (bioconda channel)
- `rust-beta5`: Trim Galore v2.1.0-beta.5 — pre-Buckberry-audit Rust baseline
- `rust-beta7`: Trim Galore v2.1.0-beta.7 — current, post-Myers' prefilter

## Cores ladder

- Perl: 1, 4, 8, 16
- Rust beta.5: 1, 4, 8, 10, 12, 14, 16, 24
- Rust beta.7: 1, 4, 8, 10, 12, 14, 16, 24

20 conditions total, 220 trim_galore invocations (10 runs + 1 warmup × 20).

## Reproducer

`scripts/benchmark.sh` on the `dev` branch drives the same matrix end-to-end. The fixture path defaults to `~/benchmark_TG_oxy/SRR24827373_*_R{1,2}.fastq.gz`; override via `R1=...` `R2=...` env vars.

## Headline results

| Comparison | Wall speedup | CPU savings |
|---|---:|---:|
| Rust v2.1.0-beta.7 vs Perl 0.6.11 — single-thread (`-j 1` / `--cores 1`) | 8.26× | **13.49×** |
| Rust v2.1.0-beta.7 vs Perl 0.6.11 — nf-core default (`-j 8` / `--cores 8`) | 4.54× | **5.93×** |
| Rust v2.1.0-beta.7 vs Rust v2.1.0-beta.5 (post-audit speedup, cores=1) | 2.74× | 2.73× |

## Charts

The three benchmark page charts (`benchmark_wall_time.png`, `benchmark_threads_cpu.png`, `benchmark_scaling.png`) are regenerated from this data via:

```bash
python3 docs/scripts/generate-benchmark-charts.py \
    --data docs/perf_data/buckberry-2026-04-29 \
    --out  docs/src/assets/benchmarks
```

(`docs/scripts/generate-benchmark-charts.py` is matplotlib-only, no extra dependencies.)
