# Progress: Reader / serial-stage attribution

**Last updated:** 2026-08-01

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | 📋 Planned | No `PLAN.md`. Spike complete; its recommendation is the input to a future plan |
| Plan Review | 📋 Planned | — |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | — |
| Code Review | 📋 Planned | — |
| Coverage | 📋 Planned | — |

## Spike (outside the plan pipeline)

| Artifact | Status | Notes |
|---|---|---|
| `spikes/spike_serial_stages_oxy.py` | ✅ Complete | Run on `dockyard-oxy-0`, 40 M × 65 bp real ENA WGBS, pinned to node0 physical cores. Two full 5-arm × 6-core matrices (compute-only and `SPIKE_PERSIST=1`) |
| `spikes/SPIKE_serial_stages_oxy.md` | ✅ Complete | Verdict: the ceiling is the **serial output writer**, not the reader |

**Headline:** stage B (main-thread `write_all`) is **61.60 s of 88.38 s = 69.7%** of runtime under `--dont_gzip` and ~34% under default gzipped output, once output actually reaches storage. Confirmed on the full matrix at n=3, sentinel drift 2.3%. Amdahl serial fraction **99.7%** at 32 cores; 2→32 scaling **1.01×** against an ideal 16×.

The writer is already at the device limit (125 MB/s vs `dd`'s 113 MB/s) but is *serialised* against compute behind a 26 ms queue — so **~1.38×** is available from moving it to its own thread, **conditional on slow/shared storage** (1.41× at 113 MB/s, ~1.16× on a cloud volume, ~1.02× on fast local NVMe). rapidgzip drops to **~1.04×**.

**Second, independent defect:** `--cores 4` beats `--cores 8` in every plain-input arm, both regimes, both read lengths — `writer_off` +17.3%, `plain` +4.7%, `synth` +4.6%. Absent in both gz-input arms. Compute-path, not writer. Unexplained; diagnosing it needs `READER_BATCH_SIZE`/channel-depth source edits.

**Awaiting user decision** on whether to open a plan for the writer decoupling, and on the 4→8 regression.

## History

- 2026-08-01: `SPIKE_PERSIST=1` full matrix complete — confirmed 69.7% at n=3; corrected decompression share (3.7%) and established that the 4→8 regression survives writeback. Report §4 now carries both regimes; limitation #1 closed
- 2026-08-01: Spike complete — writer-side attribution, supersedes the antecedent's reader hypothesis
- 2026-08-01: Directory created; antecedent is `plans/rapidgzip-trial/spikes/SPIKE_reader_bound.md`
