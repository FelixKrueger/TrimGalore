# Spike — is the per-allocation overhead resident in full, or does it ride occupancy? (#457, Phase 4)

**Date:** 2026-08-19 · **Platform:** darwin/arm64, 10 CPUs · **Binary:** `dev` @ `3fa4de3`
**Follows:** `SPIKE_slot_cost.md`, which measured the overhead but held everything resident by
construction and so could not speak to occupancy.

## 1. Question, criteria, strategy

**Question.** `SPIKE_slot_cost` found ~20 B per `String` allocation, 107–132 B/pair. In the real flush
cycle, is that charged at full `slots` per bin (allocator retains its pages) or does it ride the text's
`0..B` oscillation at occupancy `a < 1`? It decides whether the model charges `120·slots` or
`a·120·slots`.

**Success criteria.** (1) Establish whether freed pages are reused, with a clear signal. (2) Measure
`σ_effective(r)` under production-like cycling and say which candidate charge it supports.

**Strategy.** A pool of `n_bins` bins, filled **pseudo-randomly** and flushed at
`raw_bytes >= budget`, mirroring `PairedBin::push` including the one `reserve_exact` per bin per cycle.
Round-robin filling was rejected: it makes every bin reach its budget on the same record and flush in
lockstep, which would report occupancy 1.0 by construction. Peak RSS is taken at three bin budgets, so
the slope in `B` gives `σ_effective(r) = d(peak)/dB / n_bins` with no static and no worker term to
confound it. Regressing `σ` on `1/r` then yields `a` (intercept) and the per-record charge (slope).

**Out of scope.** Threading, gzip, channels, the real dispatcher's batching, other platforms.

## 2. Script and how to run

- `plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_occupancy.rs`
- `plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_occupancy.py`

```bash
python3 plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_occupancy.py . 2
```

Stages the `.rs` into `examples/`, builds, runs, and removes it in a `finally`.

## 3. Results per iteration

**#1 — blocked on `ps`.** The first design sampled *current* RSS from outside at checkpoints, because
`ru_maxrss` is a high-water mark and cannot observe a decrease. The sandbox blocks `ps`
(`PermissionError: Operation not permitted`). Rather than bypass the sandbox or add mach FFI, the
experiment was redesigned.

**#2 — redesigned so peak RSS is sufficient.** The discriminator became *live-set size versus cumulative
allocation*: the run allocates `CYCLES ×` the pool's worth of records while never holding more than one
pool live. Reuse shows as peak tracking the pool; no reuse would show as peak approaching `CYCLES ×`
pool. That is a factor-of-8 signal needing no external observer. The `retain` mode remains in the `.rs`
but is not driven.

```
 len       r  flushes    coef  vs shipped 1.4375  vs plan a+74/r
  25    96.6      128   3.342             +1.905          +1.639
  51   148.4      128   2.183             +0.745          +0.747
 150   345.9      123   1.347             -0.090          +0.196

  a (occupancy)     =  0.513     plan 0.937
  per-record charge = 266.90 B   plan 74; fully resident (spike 1) 129.12

 reuse check, largest budget:
   25 bp: peak 1973.0M, one pool 1484.7M, 8x pool 11877.3M -> peak/pool = 1.33
   51 bp: peak 1389.9M, one pool 1242.9M, 8x pool  9943.5M -> peak/pool = 1.12
  150 bp: peak  924.7M, one pool  890.0M, 8x pool  7120.0M -> peak/pool = 1.04
```

**#3 — checked that 8 cycles is steady state, not a transient.** The arc has been bitten by ramp effects
before (§4 Fault 1), so 25 bp at 40 MiB was rerun at 8, 16 and 32 cycles:

```
 cycles        pairs  flushes   peak MiB  peak/8cyc
      8   27,962,026      128     1973.0      1.000
     16   55,924,053      256     1973.0      1.000
     32  111,848,106      512     1974.0      1.001
```

Peak is flat to 0.05% while allocation volume rises 4×. Steady state at 8 cycles, and page reuse
confirmed a second way.

## 4. Findings

1. **Freed pages are reused — that candidate is dead.** 4× the allocation volume produces an identical
   peak (1973.0 vs 1974.0 MiB), and peak/pool is 1.04–1.33 rather than 8. The overhead does not
   accumulate across flush cycles.

2. **Text occupancy measures `a = 0.513`.** That is the mean of a uniform distribution on `(0, B)`,
   which is exactly what random bin assignment plus flush-at-`B` produces. Note this contradicts the
   plan's own §2.4, which argues "text occupancy is ~uniform on `(0, B)`" and then concludes
   `a = 0.937` — uniform on `(0, B)` has mean 0.5, so the stated mechanism and the stated value do not
   agree. The mechanism is right; 0.937 must be absorbing something else.

3. **The per-record charge behaves as fully resident, and then some.** 266.9 B/record in the cycling
   model against 129.1 B/record measured fully resident in `SPIKE_slot_cost`. Cycling adds a
   fragmentation penalty *on top of* full residency, and it grows with allocation density: peak/pool is
   1.33 at 25 bp, 1.12 at 51 bp, 1.04 at 150 bp.

   Mechanism, consistent with all three numbers: allocations from different bins are interleaved in
   time, so a mimalloc page holds blocks belonging to many bins. When one bin flushes, its blocks are
   scattered across pages that other bins still occupy, so those pages cannot be released or cleanly
   reused. Shorter reads mean more allocations per byte, hence more scattering.

4. **Neither candidate charge is correct.** It is not `a·slots` (finding 3) and it is not plain `slots`
   (finding 3 again, from the other side — plain `slots` still under-charges). The true form is
   `slots + fragmentation(r)`, with the fragmentation term rising as `r` falls.

5. **The isolated model overstates production, so its magnitude is not transferable.** Against the
   production matrix's 10-cell fit (`σ(r) = 0.786 + 108.74/r`):

   | read length | spike `σ` | production-model `σ` | ratio |
   |---:|---:|---:|---:|
   | 25 bp | 3.342 | 1.912 | 1.75× |
   | 51 bp | 2.183 | 1.519 | 1.44× |
   | 150 bp | 1.347 | 1.100 | 1.22× |

   The **shape** (falling in `r`) and the **sign** (fragmentation adds) match; the **level** does not.
   The isolated model is single-threaded, drops each batch immediately rather than handing it to a
   worker, and does not exercise mimalloc's thread-local heaps. Attribute the mechanism here; take the
   magnitude from the shipping configuration.

## 5. Reference snippets

Pseudo-random bin choice is load-bearing, not incidental — with round-robin every bin crosses its budget
on the same record:

```rust
let b = (rng.next_u64() >> 33) as usize % n_bins;
```

The slope-in-`B` trick, which is what lets a spike answer what the pipeline fit could not:

```python
slope, _ = lstsq(budgets_bytes, peaks_bytes)
coef = slope / N_BINS          # = sigma_effective(r), no static, no worker term
s_fit, a_fit = lstsq([1 / r for r in rs], coefs)   # intercept a, slope s
```

## 6. Recommendation

**This closes the open question from `SPIKE_slot_cost` finding 6, and it closes it against both
candidates.** The charge is `slots` plus a read-length-dependent fragmentation term.

The consequence for §6 is the same one that killed Option D, arriving a second time by a different
route: **there is no per-record constant that is correct at more than one read length.** D died because
the spine is only ~55% of the per-slot cost; the proposed "third option" (charge the allocation count at
a measured constant) dies because the fragmentation term on top of that count is itself
read-length-dependent.

**What survives is a conservative bound calibrated at the shortest supported read length — Option C.**
Two measurements now point the same way: C's logic is untouched by any of this, and only its numbers
move. Recommend re-opening §6 with C as the presumptive answer and the floor recomputed against
measured `σ` at the shortest supported `r`, on the shipping configuration rather than either spike.

## 7. Limitations

- **Magnitude is not transferable to production** (finding 5). This spike attributes a mechanism; it
  does not set a constant. Any number that reaches `clump.rs` must come from the pipeline.
- **Single-threaded.** mimalloc has thread-local heaps and cross-thread frees are the classic
  fragmentation case; production allocates on the reader and frees on workers. That is very likely part
  of the 1.22–1.75× gap, and it is untested here.
- **Batches are dropped at flush, not handed to a worker,** so live sets are shorter than production's.
- **One platform**, 16 KiB pages. Linux/x86_64 has 4 KiB pages and a different fragmentation profile.
- **Fixed-length records**, as with everything else in Phase 4 (PLAN §10 Finding 2).
- **`a = 0.513` is measured in the isolated model**, so it inherits the same non-transferability. It is
  evidence about the *mechanism* of the occupancy factor, not a replacement value for `clump.rs`.
