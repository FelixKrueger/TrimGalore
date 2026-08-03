# Post-GA follow-ups

Lightweight running list of items deferred during the v2.1.0 beta/GA push.
Not a plan — just a parking lot so ideas don't get lost.

---

## 1. Empirically test SIMD binary variants (x86-64-v3 / v4 / apple-m1)

**Status**: Deferred in the `build-provenance` plan, but only on the
basis of reasoning about TrimGalore's workload (I/O-bound decode → pattern
match → encode). Phil pointed out this is **speculation** — better to test
with real data.

### What to do

1. Build three one-off variants **outside public CI and the release workflow**
   (private branch, no crates.io publish, no GHCR push):
   - `x86_64-unknown-linux-gnu` baseline (what we ship today)
   - `x86_64-unknown-linux-gnu` with `RUSTFLAGS='-C target-cpu=x86-64-v3'` (AVX2)
   - `x86_64-unknown-linux-gnu` with `RUSTFLAGS='-C target-cpu=x86-64-v4'` (AVX-512)
   - Optionally `aarch64-apple-darwin` with `-C target-cpu=apple-m1` for the
     macOS side.

2. Run against a realistic paired-end FASTQ load, e.g.:
   - 10-50M read pairs, 2×150bp, gzipped
   - Both `--cores 1` and `--cores 8` (single-thread is where SIMD gains are
     least dilutable by I/O saturation)
   - Cold and warm page cache
   - Both default Illumina auto-detect and explicit `--nextseq 20` (the 2-colour
     quality-trim path has arithmetic that may auto-vectorise)

3. Compare:
   - Wall-clock (hyperfine 5 runs, reject outliers)
   - `perf stat -e instructions,cycles,cache-misses,branches` for the pattern
     matcher if available
   - Flamegraph of the baseline (`cargo flamegraph --release`) — if the hot
     functions aren't already SIMD-saturated via `zlib-rs`, a tuned binary
     might move the needle

### Decision criterion

- **<2% delta**: confirm the current "defer SIMD variants" stance; close
  this item with a note.
- **2-5% delta**: borderline. Weigh against the infra cost (9× build matrix,
  bioconda CPU-detection wrapper, CPU runtime detection logic).
- **>5% delta**: promote to a real plan and ship for v2.2 or post-GA patch.

### Why this matters

Phil reports 5-10% for ruSTAR (STAR aligner, genuinely compute-bound). Our
reasoning says TrimGalore is dominated by zlib decode/encode which `zlib-rs`
already SIMD-accelerates internally. But:
- The Cutadapt-equivalent pattern matcher in `src/alignment.rs` is an Aho-
  Corasick + banded DP variant — this *could* auto-vectorise with AVX2/512.
- 2-colour quality trimming (`src/quality.rs`) has arithmetic hot loops.
- RRBS MspI end-repair scanning is tight byte-level work.

Reasoning is cheap; measurement is cheaper still at this scale. One afternoon
of benchmarking resolves the question definitively.

### Coupled item: CPU runtime detection

If the benchmark justifies shipping SIMD variants, we also need Phil's
runtime feature-detection pattern (`src/cpu.rs` equivalent in ruSTAR):
- Upgrade hint in the startup banner when a faster variant exists for the
  host CPU
- Pre-parse `check_cpu_compat()` guard so a SIMD-tuned binary on an
  incompatible CPU fails friendly instead of SIGILL

Scope this as part of the same decision — CPU runtime detection alone (no
SIMD variants) is dead code.

---

<!-- Future items: add new ## headings above this line -->
