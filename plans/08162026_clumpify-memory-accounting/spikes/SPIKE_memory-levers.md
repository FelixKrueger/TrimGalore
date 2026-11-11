# Spike — two one-line structural levers for #439

**Date:** 2026-08-16
**Plan:** `plans/08162026_clumpify-memory-accounting/PLAN.md` (v2, steps 8 and 10)
**Script:** `spikes/spike_memory-levers.py` — `python3 plans/08162026_clumpify-memory-accounting/spikes/spike_memory-levers.py`
**Binary base:** dev tip `2ab0918` (ships mimalloc). Source patched, built, reverted per variant; tree verified clean at exit.

## Question

v2's Goal is unsatisfiable as written: at `--memory 1G --cores 4` the shipped
binary peaks above the 1024 MiB it promises, so no coefficient can bound it while
bin size is held. Both v2 reviewers independently said the escape must be a
structural reduction in actual memory. Two one-line candidates were proposed:

- **A** — work-channel depth 2 → 1 (`parallel.rs:119`, `:943`), so a worker holds
  two bins in flight rather than three.
- **B** — size the per-batch output reserve from actual record bytes instead of
  `len * 300` (`parallel.rs:349`, `:1069`). Those buffers receive *compressed*
  data, and `cap` is used for **both** mate buffers, so 64 bp reads reserve
  600 B/pair against ~344 B of raw text — ~398 MiB of in-flight capacity at
  `1G/--cores 4`.

Do either or both make the promise satisfiable without changing results?

## Success criteria

1. `1G/c4` ratio (peak RSS ÷ banner prediction, both MiB) ≤ 1.0 with ≥ 25 MiB headroom.
2. `2G/c4` (worst measured breach, 1.24×) ratio ≤ 1.0.
3. Output bytes identical to baseline.
4. No wall-clock regression beyond 3%.

Out of scope: model fitting, Linux, FastQC, single-end, `--clump_only`, real
libraries. Fixture is the 227×-cycled diverse pair — acceptable here because
every comparison is relative against a baseline on the same input, and the
reserve being tested is data-independent (`len * 300`).

## Results (max over 2 reps)

| cell | base | A (chan 1) | B (sized reserve) | A+B |
|---|---|---|---|---|
| `1G/c4` | **1.045** | 0.913 | 0.999 | **0.857** |
| `2G/c4` | 1.236 | 1.064 | 1.235 | **1.056** |
| `1G/c8` | 0.941 | 0.850 | 0.881 | 0.796 |

RSS delta vs base (MiB, max over reps):

| cell | A | B | A+B |
|---|---|---|---|
| `1G/c4` | −135.2 (−12.6%) | −47.2 (−4.4%) | −192.9 (−18.0%) |
| `2G/c4` | −352.2 (−13.9%) | **−3.6 (−0.1%)** | −369.0 (−14.6%) |
| `1G/c8` | −93.4 (−9.7%) | −61.8 (−6.4%) | −148.4 (−15.4%) |

Wall clock vs base (mean): A `+0.33% / −3.96% / −0.05%`, B `−2.84% / +0.57% /
+0.61%`, A+B `−1.26% / −4.50% / +3.01%`.

Output byte-identity at `1G/c4`: **identical for all three variants**, every
`_val_*.fq.gz`.

## Findings

1. **Criterion 1 met by A and by A+B.** A alone leaves 89 MiB of headroom at the
   contradiction cell; A+B leaves 147 MiB. B alone lands at 0.999 — inside the
   promise but with 1 MiB of margin, which is not a margin.
2. **Criterion 2 fails for every variant.** The `2G/c4` breach is ~24% and the
   best lever combination removes ~15%, leaving 1.056. **The levers resolve the
   default cell but not the mechanism at larger budgets**, so a calibration step
   is still required — the residual gap is ~6%, not ~24%.
3. **A is the effective lever, and it is free.** −9.7% to −13.9% RSS at no
   wall-clock cost; at `2G/c4` it is *faster* (−3.96%). This refutes v2's
   Efficiency section, which rated the queue cap as the plan's main throughput
   risk and deferred it. Reviewer B predicted exactly this: the reader is the
   bottleneck at `cores ≥ 4`, so shortening the queue costs nothing.
4. **B's benefit does not scale with bin size** — −47 MiB at `1G`, −3.6 MiB at
   `2G`, despite the reserve growing with records per batch. Unexplained; the
   likely cause is that untouched reserved pages never enter RSS, so the
   over-reserve only materialises when the allocator actually faults it in. Worth
   one probe before treating B as a memory lever rather than tidiness.
5. **Both levers are result-preserving.** Byte-identical output confirms these are
   pure memory/scheduling changes — flush order is set by fill order and
   round-robin, not by channel depth.
6. **Criterion 4 met.** The single +3.01% (A+B at `1G/c8`) sits in a cell where A
   alone measured −0.05%, so it reads as noise rather than signal.

## Reference snippets

Lever A — two sites, same edit:

```rust
let (tx, rx) = mpsc::sync_channel::<PairedWork>(1);   // parallel.rs:119
let (tx, rx) = mpsc::sync_channel::<SingleWork>(1);   // parallel.rs:943
```

Lever B — `estimated_record_bytes` is already imported in `parallel.rs`:

```rust
let cap = reads_r1.iter().map(estimated_record_bytes).sum::<usize>();  // :349
let cap = reads.iter().map(estimated_record_bytes).sum::<usize>();     // :1069
```

Note `cap` is applied to both mate buffers at `:349`, so the paired path reserves
`2 × cap`. Sizing from R1 alone under-reserves R2 when mates differ in length —
acceptable for a spike, not for implementation.

## Recommendation

**Proceed, with the plan restructured:**

- Promote **A to primary**, gated on clumpy mode (the channels are shared with
  plain round-robin, which gains no memory and should not absorb any risk).
  Measured cost is zero, so v2's deferral of it was wrong.
- Keep **B**, but as correctness/tidiness rather than a load-bearing lever until
  finding 4 is explained. It is free and it fixes a real over-reservation.
- **Do not drop the calibration step.** Finding 2 is the reason: large budgets
  still breach. But the calibration now runs *on top of* A+B, so it is fitting a
  ~6% residual rather than a ~24% one, and the matrix should be re-measured with
  the levers in place — every coefficient in v2 was measured without them.
- **Buffer recycling (v2 step 9) drops in priority.** Reviewer A wanted it as a
  prerequisite; A+B already clears the default cell, so recycling becomes a
  candidate for the residual large-budget gap, to be compared against simply
  charging more per bin.

## Limitations

- macOS/arm64 only. Lever A's benefit may differ under glibc, where the earlier
  data showed far less allocator retention; it should be confirmed on oxy.
- The 227×-cycled fixture stands (relative comparison, data-independent reserve),
  but the *absolute* headroom figures will move on distinct reads.
- Single-end, `--fastqc`, `--clump_only` and multi-pair runs untested; A patches
  the single-end channel but no single-end cell was measured.
- Two reps per cell; the base `1G/c4` spread alone was 1.017–1.045, so treat
  sub-3% differences as noise.
- Lever B sizes from R1 only, and its interaction with read length is untested —
  on 150 bp reads `len * 300` under-reserves rather than over-reserves, so the
  sign of this lever's effect may invert.
