# SPIKE — Which serial stage caps TrimGalore's parallel path?

**Date:** 2026-08-01
**Host:** `dockyard-oxy-0` — Xeon 6975P-C, 64 physical / 128 logical, 2 NUMA nodes, 991 GB RAM, AL2023 container on an overlay filesystem
**Binary:** `trim_galore 2.3.0` @ `7201f03` (oxy checkout detached at `origin/dev`)
**Antecedent:** [`plans/rapidgzip-trial/spikes/SPIKE_reader_bound.md`](../../rapidgzip-trial/spikes/SPIKE_reader_bound.md)

---

## 1. Question, success criteria, strategy

The antecedent measured that parallel gzip decompression is worth ~1.12× to us (vs 1.57× to `fastqc-rust`) and concluded the real bottleneck is "the single reader thread as a whole." That attribution was made at 10 cores on a 150 bp synthetic fixture, where no serial stage can be separated from any other.

**Question.** Single-end at `--cores N` there are three serial stages:

| | Thread | Work | Granularity |
|---|---|---|---|
| **A** | `fastq.rs:286` | `MultiGzDecoder` + `read_next_direct` + 4096-batching, **fused on one thread** | per **record** |
| **C** | `parallel.rs:947` | round-robin dispatch into per-worker `sync_channel(2)` | per **batch** |
| **B** | `parallel.rs:955` | main thread: `BTreeMap` reorder + serial `write_all` | per **byte** |

Which is the ceiling, and what is available if it is fixed?

Confirmed at runtime: **35 threads at `--cores 32`** = 32 workers + A + C + main. The `--help` text's "N+4" counts two decompressors, which is the paired-end case.

**Success criteria** (fixed before running): reproduce or refute the plateau on real data off the laptop; attribute it to a named stage by removing that stage and seeing the plateau move; quantify the serial fraction. Decision rule: *if B ≥ A_parse at max cores, the antecedent's "profile the reader" recommendation is wrong and the work is writer-side.*

**Strategy.** Five arms over `--cores [1,2,4,8,16,32]`, each removing or shrinking one stage, all pinned to node0's physical cores. Real ENA data (`SRR24827373_1`, Buckberry WGBS, 40 M × 65 bp) as primary; the antecedent's 150 bp fixture retained as a second read length so the serial cost can be split into per-record and per-byte terms.

**Out of scope.** Integrating rapidgzip; uBAM/BGZF; paired-end; `--clumpify`; `--fastqc`.

---

## 2. Script

`plans/reader-serial-stages/spikes/spike_serial_stages_oxy.py`

```bash
# on oxy, after building inputs into /tmp/spike_reader_oxy
SPIKE_READS=40000000 SPIKE_REPEATS=3 python3 -u spike_serial_stages_oxy.py
```

`SPIKE_PERSIST=1` (now the default) is load-bearing — see iteration #3.

---

## 3. Results per iteration

### #1 — smoke run, 2 M reads

Harness worked end-to-end. Two outcomes, one a finding and one a bug.

The real 65 bp arm was **flat from `--cores 2`** (1.89 s → 1.85 s, 1.02× to 32 cores) while the 150 bp synthetic still scaled 1.91×. Not an input-size artifact: scaling the input scales `S` and `P` together and leaves their ratio fixed. 65 bp reads carry too little per-read trim work to keep workers fed. The laptop only saw 3.0× scaling because its fixture was 150 bp.

The bug: my MB/s-vs-reads/s discriminator compared arms at raw wall clock, where the synth arm is nowhere near its ceiling. Replaced with a two-term solve on each arm's *fitted* `S`, which isolates the serial stage from per-read trim work. Also: synth arm promoted to the full core sweep (a 3-point fit is too noisy to feed the solve), sentinels interleaved between every arm, contention verdict gated on the term being >2% of runtime.

### #2 — full sweep, 40 M reads

Completed, ~45 min, 69 runs. Full output in §4 of this report's data tables below.

The sentinel logic reported **VOID** at −15.8% to −17.6% drift. That verdict was wrong: all five in-sweep sentinels cluster at 26.59–27.18 s (2.2% spread) and agree exactly with `plain cores=8` measured in the sweep proper (27.14 s). Only the **baseline** was the outlier at 32.28 s — it ran immediately after the input build had written ~17 GB and was competing with kernel writeback. The sweep is internally consistent; the reference point was cold.

### #3 — the confirmation loop that invalidated the sweep's headline number

Re-running `plain` at 4 and 8 cores interleaved, to check the surprising 4-core peak:

```
rep1 cores=4 22.90s    rep1 cores=8 26.81s     <- reproduces the sweep exactly
rep2 cores=4 89.37s    rep2 cores=8 88.94s     <- 3.3x slower
rep3 cores=4 85.35s
```

The difference between rep1 and rep2+ is my harness, not the machine. `run_once` calls `shutil.rmtree(outdir)` immediately after each run, so **each 7.9 GB output was deleted before the kernel ever flushed it** — the write cost a memcpy into page cache that was then discarded. My confirmation loop overwrote the same file instead, forcing real writeback.

The sweep therefore measured stage B in the most favourable possible regime, which is why `plain − writer_off` came out at **−1.0%**.

### #4 — stage B under real writeback

Alternating, `sync` inside the timed window:

```
writer_off (0 B output)      27.20s   27.07s
plain      (7.84 GB output)  89.58s   89.98s
```

**B = 62.7 s of 89.8 s — 70% of runtime.** Not −1%. This also vindicates our own `--help` text ("beyond `--cores 8`, gzip-output I/O on the storage layer typically becomes binding"), which the sweep appeared to contradict only because the sweep never let the storage layer see the bytes.

### #5 — is 125 MB/s our code or the device?

```
dd if=/dev/zero of=… bs=1M count=8000 conv=fdatasync   ->  113 MB/s
TrimGalore writer: 7.84 GB / 62.7 s                    ->  125 MB/s
```

The writer is **already faster than raw sequential `dd`**. There is no inefficiency inside the write. But it is *serialised against compute*:

| config | time |
|---|---|
| gz→gz, compute only (`--length 10000`, ~0 output) | 29.12 / 29.51 s |
| gz→gz, output persisted | 40.29 / 40.48 s |
| delta | **10.9 s** ≈ 1.367 GB / 113 MB/s = 12.1 s |

Additive, not overlapped. `result_rx` is `sync_channel(cores * 2)` = 16 slots at `cores=8`, ≈2.9 MB of compressed batches ≈ **0.026 s** of write time — so the instant the main thread blocks in `write_all`, the queue drains, workers stall, and thread A stalls behind them.

### #6 — full matrix under writeback

Re-ran all five arms × six core levels with `SPIKE_PERSIST=1` (§4a). Sentinel drift **2.3%** worst case against the in-sweep median, so internally consistent.

Two corrections to #4–#5 came out of it, recorded in findings 5 and 7: decompression is 3.7% (not the 2.4% I extrapolated), and the 4→8 core regression **survives** writeback rather than being masked by it — my initial reading came from two n=1 samples that straddled the noise.

It also validated the sentinel fix. The baseline again fell outside the in-sweep cluster (−4.9%, this time *faster*, since the inputs were already page-cached from the previous sweep), and referencing drift to the in-sweep median correctly reported consistency instead of the spurious `VOID` the baseline-referenced version produced in #2. The baseline run is structurally unlike the rest of the sweep in whichever direction the cache happens to favour — which is the argument for never using it as a reference.

---

## 4. Data tables (40 M reads, cores pinned 0-31)

Both regimes were run as full 5-arm × 6-core matrices. Wall clock, min of n (n=3 primary, n=2 secondary, n=1 at `cores=1`). Ideal 2→32 is **16.00×** in both.

### 4a. Writeback forced (`SPIKE_PERSIST=1`) — the realistic regime

`sync` inside the timed window, so every number is end-to-end.

| arm | 1 | 2 | 4 | 8 | 16 | 32 | 2→32 |
|---|---|---|---|---|---|---|---|
| gz | — | 92.87 | 92.09 | 93.27 | 91.64 | 91.82 | 1.01× |
| plain | 96.82 | 88.99 | **84.68** | 88.87 | 88.98 | 88.38 | 1.01× |
| writer_off | 53.90 | 26.52 | **22.86** | 26.81 | 26.95 | 26.78 | 0.99× |
| writer_cmp | 104.29 | 47.43 | **40.43** | 40.54 | 40.40 | 40.46 | 1.17× |
| synth 150 bp | 92.13 | 88.75 | **80.71** | 84.40 | 82.92 | 81.38 | 1.09× |

Attribution at 32 cores: decompression **3.43 s (3.7%)**, serial write **61.60 s (69.7%)**, residual `A_parse + C + trim/N` **26.78 s (30.3%)**.

Amdahl: `plain` S=88.16 s, serial **99.7%** of runtime at 32 cores. `writer_off` S=26.42 s, 98.6%. Under writeback the parallel term is not merely small, it is invisible.

### 4b. Output discarded before flush — the compute-side regime

Retained because it isolates the compute path. Every absolute number understates real runtime.

| arm | 1 | 2 | 4 | 8 | 16 | 32 | 2→32 |
|---|---|---|---|---|---|---|---|
| gz | 70.08 | 30.48 | 28.99 | 29.18 | 29.11 | 29.59 | 1.03× |
| plain | 59.34 | 28.79 | **22.88** | 27.14 | 26.91 | 27.41 | 1.05× |
| writer_off | 56.01 | 27.77 | **22.81** | 27.13 | 27.10 | 27.69 | 1.00× |
| writer_cmp | 106.73 | 48.86 | 30.42 | 30.18 | 30.41 | 30.34 | 1.61× |
| synth 150 bp | 94.90 | 52.46 | **26.67** | 31.31 | 29.33 | 26.80 | 1.96× |

The two regimes agree on the compute path to under 1% where they should: `writer_off` at 4 cores is 22.86 s vs 22.81 s, and at 8 cores 26.81 s vs 27.13 s — despite host load average differing by roughly 5× between the runs (6.5–53.5 vs 28–60). That cross-regime agreement is the strongest evidence the `taskset 0-31` pin plus min-of-n controls work.

### 4c. Serial-stage cost model, both regimes

Solved from the two read lengths' fitted `S`. The inputs are **byte-matched by construction** (7.887 GB vs 7.886 GB), which is what makes `a` well-conditioned — it is determined almost entirely by the 16.66 M difference in record count at equal bytes.

| regime | a (ns/record) | b (ns/byte) | implied per-byte ceiling |
|---|---|---|---|
| compute-only | 155.4 | 2.546 | 393 MB/s |
| writeback | 415.8 | 9.069 | **110 MB/s** |

**`b` under writeback is the device, not our code.** 110 MB/s against `dd conv=fdatasync`'s measured 113 MB/s. The byte term has stopped describing the reader and started describing the disk, which is the correct behaviour of the model and a useful confirmation that it is measuring what it claims.

**`a` under writeback is NOT interpretable, and should not be quoted as a per-record cost.** It tripled (155 → 416 ns) when a purely per-*byte* stage was added, which is a contradiction, not a finding. The cause is a missing term: the model has no *output*-bytes variable, and the two arms' output volumes differ (real WGBS carries real adapter content and trims more than the repeated synthetic fixture). That difference is collinear with nothing else in the model, so it lands on `a`. Use the compute-only `a = 155 ns/record` for reasoning about the reader; treat the writeback `a` as an error bar on the fit's specification.

Per-byte share rises from 76–85% (compute-only) to 81–88% (writeback), in both cases pointing at the byte path.

---

## 5. Findings

**1. The plateau is real and far worse on real data than the laptop implied.** 2→32 cores buys 1.03–1.05× on real 65 bp WGBS against an ideal 16×. The antecedent's 3.0× came from a 150 bp fixture whose longer alignment DP keeps workers fed; short reads are serial-bound from two cores.

**2. Which stage depends entirely on whether output reaches storage — and the first sweep got this wrong.** With output discarded before flush, stage A's per-byte path is the ceiling (~27 s, ~290 MB/s) and B is free. With output persisted, **B is 61.60 s of 88.38 s = 69.7%** under `--dont_gzip` (full matrix, n=3, sentinel drift 2.3%) and ~34% under the default gzipped output. The decision rule set in §1 fires: **B ≥ A_parse, so the antecedent's "profile the reader" recommendation is wrong.**

Under writeback the Amdahl serial fraction is **99.7%** at 32 cores and 2→32 scaling is **1.01×** against an ideal 16×. Cores buy essentially nothing on this storage.

**3. The writer is at the device limit, so the fix is not a faster write — it is an overlapped one.** 125 MB/s vs `dd`'s 113 MB/s leaves nothing on the table inside `write_all`. The waste is that 10.9 s of write time lands *additively* on 29.3 s of compute because the main thread both reorders and writes, behind a 26 ms buffer.

**4. The available win, and its condition.** For the default gz→gz path, compute (29.3 s) and write (10.9 s) are both under the device rate, so a dedicated writer thread with a modest queue would fully overlap them: T → 29.3 s, **~1.38×**. No new dependency. That is materially larger than rapidgzip's ~1.12×.

   The size of this win is **storage-dependent**, because the hidden term is `bytes / device_bandwidth`:

   | device | write time | speedup from overlapping |
   |---|---|---|
   | 113 MB/s (this container overlay) | 12.1 s | **1.41×** |
   | ~300 MB/s (NFS / cloud volume) | 4.6 s | 1.16× |
   | ~2 GB/s (fast local NVMe) | 0.7 s | 1.02× |

   So this is worth doing to the extent that users write to shared or networked storage — which in bioinformatics is the common case, not the exception. On a fast local NVMe it is nearly worthless.

   Under `--dont_gzip` the output (7.84 GB) exceeds what the device can absorb within compute time, so that path stays device-bound: floor 69.4 s against a measured 89.8 s.

**5. rapidgzip is now weaker than the antecedent found.** Decompression = `gz − plain` = 2.17 s = 7.3% of `gz` in the compute-only regime, and **3.43 s = 3.7%** measured directly under writeback. The ceiling on perfect parallel decompression is **~1.04×**, against the antecedent's ~1.12×. Park it harder.

**6. Within stage A, optimise bytes not records.** 155 ns/record against 2.55 ns/byte, with per-byte at 76–85% of the serial cost. 393 MB/s is slow for a byte path, which points at record construction copying into owned `String`s rather than at allocation count — but this was fitted in the compute-only regime and is not itself a profile.

**7. A reproducible 4→8 core regression that survives writeback, unexplained.** `--cores 4` is the fastest setting in **every plain-input arm, in both regimes**:

| arm (writeback) | c4 | c8 | penalty at 8 |
|---|---|---|---|
| `writer_off` | 22.86 | 26.81 | **+17.3%** |
| `plain` | 84.68 | 88.87 | +4.7% |
| `synth` 150 bp | 80.71 | 84.40 | +4.6% |

Absent in both gz-input arms (`gz` 92.09 → 93.27 = +1.3%; `writer_cmp` 40.43 → 40.54 = +0.3%), in both regimes. That split — three of three plain-input arms, zero of two gz-input arms, across two independent sweeps and two read lengths — makes it a property rather than an oddity.

The penalty is largest where the writer is absent and shrinks in proportion to how much runtime the writer adds, which places the defect in the **compute path** (A / C / workers) with the writer diluting rather than causing it. Corroborating: `writer_off`'s 2-parameter Amdahl residual at 4 cores is −13.1%, the largest of any arm, and its 3-parameter contention term is 9.5% of runtime at 32 cores.

Candidates, unresolved: round-robin head-of-line blocking on the depth-2 per-worker channels (stage C), or L3 pressure from more concurrent 4096-record batches. Distinguishing them needs `READER_BATCH_SIZE` / channel-depth changes, i.e. source edits, which were out of scope here.

**This is a user-visible defect independent of any optimisation: `--cores 4` beats `--cores 8` end-to-end on plain short-read input, by ~5%.**

**8. Both harnesses that flagged a problem were wrong before the data was.** The sentinel's VOID verdict came from a cold baseline; the sweep's "B is free" came from `rmtree` beating writeback. In both cases the instrument, not the machine, produced the anomaly.

---

## 6. Reference snippets worth carrying forward

Not the measurement plumbing — the two controls that decided the answer.

```python
# The control the first sweep shipped without. Deleting the output before the
# kernel flushes it turns a disk write into a discarded memcpy: the same cell
# measured the serial writer at -1% of runtime without this, and 70% with it.
r = subprocess.run(cmd, capture_output=True, text=True)
if PERSIST:
    subprocess.run(["sync"], check=False)   # inside the timed window
dt = time.perf_counter() - t0
```

```python
# Compare sentinels to EACH OTHER, never to a baseline taken right after the
# input build -- that one competes with ~17 GB of pending writeback and reads as
# a 16% machine drift that did not happen.
sentinels = [(arm, t, (t / sent_base - 1) * 100) for ...]
```

```bash
# Before concluding anything about a write path, ask whether the number is the
# code or the device. TrimGalore writes at 125 MB/s; dd manages 113 MB/s.
dd if=/dev/zero of=probe bs=1M count=8000 conv=fdatasync
```

Architectural facts established (worth not re-deriving):

- Single-end is **N + 3** threads: N workers + A (decompress+parse+batch, fused) + C (dispatch) + main-as-writer. Verified by thread count at runtime.
- `result_rx` is `sync_channel(cores * 2)` — 16 slots ≈ 2.9 MB ≈ 26 ms of write time at `cores=8`. Far too shallow to decouple the writer.
- `--length 10000` discards 100% of reads post-trim and writes a **0-byte** output, which zeroes stage B while retaining all of A and all trim work. Zero-source-change writer isolation.
- `taskset -c 0-31` pins to node0's physical cores, excluding SMT siblings (64-95). Affinity is inherited by every spawned thread.
- There is **no `--gzip` flag**; forcing compressed output means giving gz input and omitting `--dont_gzip`.

---

## 7. Recommendation

**Decouple the output writer from the main thread, and drop the reader-parsing work down the list.**

1. **Move `write_all` onto its own thread behind a deeper queue.** The main thread keeps the `BTreeMap` reordering; the writer thread drains ordered buffers. Worth ~1.4× on slow/shared storage, ~1.15× on cloud volumes, ~nothing on fast local NVMe — so scope it as an I/O-bound-user win and measure on the target storage before committing. This is the only change here with a favourable cost/benefit: no new dependency, no feature-flag matrix, and the byte-identity invariants are untouched because the bytes and their order do not change.
2. **Investigate the 4→8 core regression before any reader work.** `--cores 4` beating `--cores 8` on plain short-read input is a user-visible defect independent of any optimisation, and if it is head-of-line blocking on the depth-2 dispatch channels it may be a small fix. It also has to be understood before a scaling curve can be trusted.
3. **Reader/parse work third, framed as a per-byte problem.** 393 MB/s through the byte path with per-byte at 76–85% of serial cost. Worth attacking only after (1), because until the writer is overlapped its cost masks any gain here.
4. **rapidgzip: park it harder.** 2.4% of realistic runtime. The antecedent's ~1.12× ceiling was itself measured without writeback; the honest figure is ~1.02–1.08×.

The antecedent's option 1 — "park the codec, profile the reader" — was right about the codec and wrong about the target. The serial stage that matters on real data is the one writing bytes out, not the one parsing them in.

---

## 8. Limitations

Stated plainly, because several would change the numbers.

- ~~The main sweep's absolute times are compute-only.~~ **Closed.** The full 5-arm × 6-core matrix was re-run with `SPIKE_PERSIST=1` (§4a). It confirmed the hand-measured 70% at n=3, and corrected two things: decompression is 3.7% not 2.4%, and the 4→8 regression *survives* writeback rather than being masked by it (finding 7) — my initial hand-measurement suggested the opposite from two n=1 samples that straddled the noise.
- **The writeback cost model's `a` term is mis-specified**, and quoting it as a per-record cost would be wrong — see §4c. The model lacks an output-bytes variable, so the two arms' differing output volumes land on `a`. Adding a third arm at a different output volume would identify it.
- **113 MB/s is a slow device**, and the headline 1.38× is proportional to write time. This is a container overlay filesystem on a shared host, not representative of a dedicated NVMe. The finding is conditional on storage and I have not measured any other storage class.
- **Shared host, co-tenant load 28–60 throughout.** `ps` cannot see other containers' processes, so load average was the only signal. Interleaved sentinels agree to 2.2%, which bounds drift but does not eliminate it.
- **`taskset`, not `numactl`.** CPU affinity pinned; memory not explicitly bound. Relies on first-touch locality, which is not guaranteed under pressure. numactl is absent and `sudo` needs a password.
- **The per-record/per-byte split is exactly determined**, so it has no residual and cannot be validated. A third read length would make it testable. The two read lengths also differ in content, not just length.
- **Single-end only.** Paired-end doubles stage A and writes two outputs; the writer decoupling argument should be stronger there, but it is untested.
- **The 4→8 regression has a mechanism-shaped hypothesis and no evidence.** I did not vary `READER_BATCH_SIZE` or channel depth, which would need source changes.
- **Untested:** `--clumpify`, uBAM in or out, `--fastqc`, `--cores` above 32, and whether the writer decoupling interacts with the `--memory` budget.
- The synthetic 150 bp arm is 2,334 concatenated copies of one 10 k fixture, so its cache behaviour flatters it and its trim work repeats.
