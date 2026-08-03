# SPIKE — Is TrimGalore reader-bound on gzipped input?

**Date:** 2026-08-01
**Prompted by:** [ewels/FastQC-Rust#7](https://github.com/ewels/FastQC-Rust/pull/7), which adds optional parallel gzip decompression via [rapidgzip-rust](https://github.com/COMBINE-lab/rapidgzip-rust) and measures **1.57× end-to-end** on an 8 GB file.

---

## 1. Question and success criteria

**Question.** TrimGalore decompresses input on **one** background thread (`fastq.rs:276` `open_threaded` → flate2 `MultiGzDecoder`, feeding a `sync_channel(4)`), while trimming *and* output compression run on N workers. Does that single decompression thread cap throughput, and from what `--cores`?

**Success criteria.** A scaling curve across `--cores`, plus an attribution test isolating decompression. Decision rule, fixed before running:

- gzipped arm scales materially worse than plain **and** costs >15% at max cores → **reader-bound**, rapidgzip worth pursuing
- the two arms within 5% at max cores → **not reader-bound**, rapidgzip buys nothing
- in between → partially reader-bound, quantify before committing

**Out of scope.** Integrating rapidgzip; the uBAM/BGZF path; output-compression tuning.

---

## 2. Strategy — attribution by differential input format

Run identical reads twice: once as `.fastq.gz`, once as plain `.fastq`. If the gzipped arm plateaus while the plain arm keeps scaling, the serial decompression thread is the ceiling. No source instrumentation required.

**Three controls, each of which would otherwise have invalidated the result:**

1. **`--dont_gzip` in both arms.** TrimGalore's output compression follows the input's by default (#245), so a plain-input arm would also skip *output* compression — two variables changed, not one.
2. **Explicit `-a`.** Adapter auto-detection scans up to 1 M reads *from the input*, so it is itself sensitive to decompression speed.
3. **`--cores 1` excluded from the scaling fit.** `main.rs` routes `cores > 1` to the worker pool and `cores == 1` to the sequential path — different architectures, not two points on one curve.

Single-end is the primary experiment (one reader vs N workers, cleanest attribution); paired-end measured separately because it spawns two reader threads.

---

## 3. Script

`plans/rapidgzip-trial/spikes/spike_reader_bound.py`

```bash
export SPIKE_SCRATCH=/path/to/scratch SPIKE_COPIES=300 SPIKE_REPEATS=2
python3 plans/rapidgzip-trial/spikes/spike_reader_bound.py
```

Input is built by concatenating `10K_150bp.fastq.gz` 300× — concatenated gzip members are a valid stream (RFC 1952), which is also what TrimGalore's own parallel writer emits, so this exercises the multi-member path rapidgzip has a fast path for.

**Environment:** macOS, Apple Silicon, 10 logical CPUs. 3,000,000 reads × 150 bp = **1014 MB uncompressed / 295 MB gzipped**. Min-of-2 wall clock; observed spread ≤2.7%.

---

## 4. Results

### Iteration log

**#1 — smoke run (20 copies, 68 MB).** Harness worked; gz penalty visible and *widening* with core count (1.13× → 1.29×), the signature of a serial stage becoming dominant. Runs of 0.25–0.7 s were startup-dominated, so not trustworthy for a verdict.

**#2 — cache-key bug, caught by reading my own output.** A "300 copies" run reported `68 MB uncompressed` — `if not gz.exists()` ignored `SPIKE_COPIES`, so it silently re-measured the 20-copy file. The verdict had flipped to READER-BOUND purely from noise crossing a threshold on the *same* input. Fixed by keying the filename on `COPIES` and adding a hard read-count assertion that exits if the built input is not the requested size.

**#3 — real run (300 copies, 1014 MB), plus a paired-end probe.** Numbers below.

### Single-end

| `--cores` | gz | plain | gz penalty | gz scaling | plain scaling |
|---|---|---|---|---|---|
| 1 *(sequential path)* | 17.26 s — 58.7 MB/s | 15.57 s — 65.1 MB/s | +10.9% | — | — |
| 2 | 8.87 s — 114.2 MB/s | 8.41 s — 120.6 MB/s | +5.5% | 1.00× | 1.00× |
| 4 | 4.88 s — 207.8 MB/s | 4.45 s — 227.6 MB/s | +9.7% | 1.82× | 1.89× |
| 8 | 3.37 s — 300.5 MB/s | 2.77 s — 365.9 MB/s | +21.7% | 2.63× | 3.03× |
| 10 | 3.19 s — 317.8 MB/s | 2.80 s — 361.7 MB/s | **+13.8%** | 2.78× | 3.00× |

Ideal scaling 2 → 10 cores is 5.00×. Achieved: **2.78× gzipped, 3.00× plain.**

### Paired-end (two reader threads, twice the data)

| `--cores` | gz | plain | gz penalty |
|---|---|---|---|
| 8 | 6.03 s | 5.30 s | +13.6% |
| 10 | 5.77 s | 5.22 s | **+10.5%** |

Two reader threads do **not** reduce the relative penalty, because there is also twice as much to decompress — per-reader load is unchanged.

---

## 5. Findings

**1. Decompression costs 10–14% at high core counts, and that is the entire upper bound on what perfect parallel decompression could recover.** Not 1.5×. At the core counts most users actually run (2–4), it is 5–10%.

**2. The plateau is not decompression.** This is the finding that matters. The **plain** arm — with zero decompression — also stops scaling at 8 cores: 2.77 s → 2.80 s, marginally *worse*. So the dominant serial bottleneck is the rest of the single reader thread: file read, `memchr` record parsing, 4096-record batching, and the bounded `sync_channel(4)` handoff. Parallelising gzip leaves all of that untouched.

**3. Why the same crate buys Phil 1.57× and us ~1.12×.** The architectures differ in one decisive way:

| | FastQC-Rust | TrimGalore |
|---|---|---|
| Decompression | serial, **on the analysis thread** | serial, **on its own thread** |
| Main work | serial (12 modules, ~137 MB/s floor) | **parallel** across N workers |
| Consequence | `T = T_decomp + T_analysis`, so decompression is ~36% of total | `T = max(T_reader, T_trim/N)`, so decompression is **already overlapped** with trimming |

Phil's decompression was pure additive serial time — parallelising it removed ~36% of his runtime. Ours is already hidden behind parallel trimming, so only the residue shows up. **Being further along on parallelism is exactly why we have less to gain here.**

**4. This inverts my pre-spike reasoning.** I argued TrimGalore might have *more* headroom than FastQC-Rust because its analysis is already parallel. The opposite holds, for that very reason.

**5. The higher-value target is the reader, not the codec.** Both arms cap near 360 MB/s. Parallelising record parsing and batching — or running multiple readers over byte ranges — attacks the serial stage that actually dominates. rapidgzip becomes *more* attractive after that work, not before, because decompression's share of a shrunken serial stage grows.

---

## 6. Reference snippets worth carrying forward

Not the codec integration — the measurement design, which is the reusable part:

```python
# Attribution by differential input format. Both arms MUST pin output
# compression, or TrimGalore's "output matches input" default (#245) changes
# two variables at once and the result means nothing.
cmd = [TG, "-a", ADAPTER,      # no auto-detect scan (it reads the input)
       "--dont_gzip",          # identical output format in both arms
       "--cores", str(cores), "-o", str(outdir), str(inp)]
```

```python
# Any cached artefact must be keyed on the parameter that generates it, and
# asserted after construction. Iteration #2 measured a stale 68 MB input while
# reporting 300 copies.
gz = SCRATCH / f"big_{COPIES}.fastq.gz"
if reads != COPIES * 10_000:
    sys.exit(f"input built wrong: {reads:,} reads, expected {COPIES*10_000:,}")
```

Architectural facts established, worth not re-deriving:

- `fastq.rs:276` `open_threaded` — one thread, `MultiGzDecoder`, `READER_BATCH_SIZE = 4096`, `sync_channel(READER_CHANNEL_BATCHES = 4)`
- `parallel.rs:336` — each worker owns its own `GzEncoder`, so output compression is already parallel
- `cores == 1` is the sequential path, not a one-worker pool

---

## 7. Recommendation

**Do not adopt rapidgzip on throughput grounds now.** The measured ceiling is ~10–14% at 8–10 cores and 5–10% at 2–4 cores, against a real cost: a new dependency, a feature-flag matrix, a git-rev pin until the crate is published, and the byte-identity burden of 24 `validation` assertions plus Perl 0.6.11 md5 parity.

Three options, in the order I would take them:

1. **Park it, and profile the reader instead.** The serial reader is worth ~2× (both arms plateau at 3.0× against an ideal 5×), which is an order of magnitude more headroom than the codec. Revisit rapidgzip afterwards, when decompression is a larger share of a smaller serial stage.
2. **If we want it anyway** — and there are non-throughput reasons, below — copy Phil's shape exactly: off-by-default feature, git rev pinned *with* a `version` key so `cargo publish` keeps working, transparent flate2 fallback, and a runtime `A/B` env switch. Put the feature combination in the CI matrix on day one; Phil flags that his own PR does not build the feature in CI, and we just landed a job for precisely that class of gap.
3. **Reply to Phil's PR** either way, with these numbers. The asymmetry is genuinely interesting — the same crate is worth 1.57× to him and ~1.12× to us, for a reason that is about architecture rather than about the crate — and it is the kind of thing worth having on the record before someone else re-derives it.

**Non-throughput reasons that could still justify it later:** rapidgzip decodes BGZF, so a future unified reader could serve FASTQ and uBAM through one path; and `zlib-rs`-only builds keep the no-C-toolchain property intact.

---

## 8. Limitations

Stated plainly, because several would change the numbers:

- **Synthetic input.** 300 identical copies of one 10 k fixture. Per-read work should be representative, but adapter content, quality distribution and duplication all repeat, and cache behaviour is likely flattering to both arms. A real multi-GB library could move the ratio.
- **10 cores only.** Phil measured on 4 vCPUs and saw his largest win on 8 GB. On a 32- or 64-core node `T_trim/N` shrinks further, so `T → T_reader` and decompression's *share of the serial stage* is what matters — I did not measure that regime.
- **Apple Silicon, macOS.** Different memory bandwidth and a different `zlib-rs` SIMD path from a Xeon runner.
- **Decompression never measured in isolation.** The gz-vs-plain gap is an *attribution*, not a direct measurement. I deliberately avoided `gzip -dc` as a proxy: Phil measured `zlib-rs` at ~3.7× the `gzip` tool even single-threaded, so that proxy would have overstated decompression cost and biased toward a false READER-BOUND verdict.
- **Disk I/O favours the gz arm** (295 MB read vs 1014 MB), so the measured penalty is if anything a *lower* bound on decompression's true cost. The real ceiling may be slightly above 14%.
- **Untested:** `--clumpify`, uBAM input, `--fastqc`, and the interaction between a parallel decoder's window buffers and the `--memory` budget.
- **1 GB, not 8 GB.** Phil's largest effect appeared at 8 GB; I could not build an input that size in this environment without disproportionate time cost.
