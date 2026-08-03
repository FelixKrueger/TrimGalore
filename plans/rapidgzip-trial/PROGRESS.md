# Progress: rapidgzip-rust trial

**Last updated:** 2026-08-01

## Status

| Step | Status | Notes |
|------|--------|-------|
| Spike | ✅ Complete | `spikes/SPIKE_reader_bound.md` — verdict: **park it** |
| Plan | ❌ Excluded | Not written; the spike does not support an implementation on throughput grounds |

## Outcome

Prompted by [ewels/FastQC-Rust#7](https://github.com/ewels/FastQC-Rust/pull/7), where Phil measures **1.57×** from parallel gzip decompression. Spiked whether TrimGalore would see the same.

**It would not — ~1.12× at best.** Decompression costs 10–14% of runtime at 8–10 cores, 5–10% at the 2–4 cores most users run. The upper bound on perfect parallel decompression is that number.

**The plateau is not decompression.** With input as plain FASTQ (zero decompression) throughput *still* stops scaling at 8 cores. The dominant serial stage is the rest of the single reader thread — file read, `memchr` parsing, 4096-record batching, the `sync_channel(4)` handoff. Both arms cap near 360 MB/s against an ideal 5× from 2→10 cores.

**Why the same crate is worth 1.57× to Phil and ~1.12× to us:** FastQC-Rust decompresses *on its analysis thread*, so decompression was ~36% of its total. TrimGalore decompresses on its own thread, already overlapped with parallel trimming. Being further along on parallelism is exactly why there is less to gain.

This inverted my pre-spike reasoning, which had argued we'd have *more* headroom because our trimming is already parallel.

## Recommended next step

Profile the **reader**, not the codec — it is worth ~2× versus the codec's ~1.12×. Revisit rapidgzip afterwards, when decompression is a larger share of a smaller serial stage.

Non-throughput reasons that could still justify adoption later: rapidgzip decodes BGZF, so a future unified reader could serve FASTQ and uBAM through one path.

## Carried forward

- Findings recorded in memory (`reference_rapidgzip_rust`) so a future session does not re-derive them.
- If adopted anyway, copy Phil's shape exactly — see the spike report §7.
- **Not yet done:** replying to Phil's PR with these numbers. The asymmetry is architectural rather than about the crate, and worth putting on the record.

## History

- 2026-08-01: Spike → ✅ Complete (3 iterations; one cache-key bug caught by reading the reported input size)
