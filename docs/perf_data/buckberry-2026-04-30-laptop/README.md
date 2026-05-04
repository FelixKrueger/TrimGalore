# Buckberry-scale benchmark — Apple M1 Pro laptop (2026-04-30)

Cross-platform companion run to the [server-side Buckberry benchmark](../buckberry-2026-04-29/) — same fixture (`SRR24827373`, 84M PE reads, 4.4 GiB gzipped), different hardware, lighter methodology since this is a laptop.

## Run config

- **Hardware**: Apple M1 Pro, 10 cores, 32 GiB RAM (macOS, Apple Silicon)
- **Binary**: Trim Galore v2.1.0-beta.7 (Apple Silicon native, `cargo install`)
- **Methodology**: `hyperfine --warmup 1 --runs 3` per condition (lighter than the server-side `--runs 10`; goal is a directional cross-platform datapoint, not paper-grade rigor)
- **Cores ladder**: 1, 2, 4, 6 (10 physical cores total; intentionally stops at 6 to leave headroom for the OS and other apps)
- **Engine**: Rust v2.1.0-beta.7 only (Perl baseline omitted — Perl c1 on Buckberry takes 45 min per iteration, prohibitive on laptop time budgets)

## Files

| Pattern | Description |
|---|---|
| `rust-beta7_c<N>.json` | Hyperfine raw output (3 sample runs, mean/min/max/stddev, user + system CPU time) |
| `rust-beta7_c<N>.md` | Hyperfine human-readable markdown summary |
| `laptop_bench.log` | Verbose stdout log from the run wrapper, including hardware probe + per-condition completion timestamps |

## Headline numbers

| `--cores` | Wall | CPU | Speedup vs c1 |
|---:|---:|---:|---:|
| 1 | 456 s | 430 s | 1.00× |
| 2 | 221 s | 482 s | 2.07× |
| 4 | 112 s | 487 s | 4.07× |
| 6 | 80 s | 501 s | 5.73× |

Apple M1 Pro is ~38% slower per-core than the Xeon 6975P-C in the [server run](../buckberry-2026-04-29/) (laptop c1: 456s vs server c1: 329s). The relative gap holds at higher core counts (c4: 112s laptop vs 81s server, same ~38% gap). At cores=6, the laptop finishes 84M PE reads in 80s — within ~40% of the server's cores=8 saturation point of 57s.
