# Spike — what does a paired-bin slot actually cost? (#457, Phase 4)

**Date:** 2026-08-19 · **Platform:** darwin/arm64, 10 CPUs · **Binary:** `dev` @ `3fa4de3`
**Branch:** `perf/clumpify-spine-accounting`

## 1. Question, criteria, strategy

**Question.** Is a paired-bin slot's non-text resident cost 148 B — `size_of::<FastqRecord>()`×2 +
`size_of::<MinimizerKey>()` = 72+72+4, as PLAN §2.3 asserts — or ~217 B, as the calibration matrix's
fitted `s = 108.74` implies? Equivalently: is `s ≈ 74` or `s ≈ 109`?

**Why a spike and not another fit.** The matrix cannot answer it. At one budget the columns `1`,
`B·n_bins` and `B·cores` depend on `cores` alone, so `static`, `a` and `b` are degenerate (PLAN §10
Finding 1). Concretely: with `a` and `b` free, `s = 108.74` at every static pin from 173 to 300 MiB; with
`a` and `b` forced to their structural values, `static = 260 MiB` returns `s = 73.9`. The fit cannot
separate a per-record cost from a constant.

**Success criteria.**
1. Non-text resident B/pair at read lengths {25, 36, 51, 100, 150}.
2. Run-to-run spread ≤ 10 B/pair, against a 69.5 B/pair gap between the two predictions.
3. A verdict that names one of: 148 confirmed / ~217 confirmed / neither.

**Strategy — measure a slope, not a level.** Build `PairedBin`'s exact shape (`r1`, `r2`:
`Vec<FastqRecord>`; `keys`: `Vec<MinimizerKey>`; all `reserve_exact`ed once) at four record counts per
read length, hold it resident past the process peak, and regress peak RSS on record count. **A slope
cancels the process baseline**, which is precisely the static term the matrix could not separate.

**Out of scope.** The static term, the occupancy factor `a`, the worker term `b`, channel depths, gzip,
threading, and the flush threshold itself.

## 2. Script and how to run

- Rust measurement: `plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_slot_cost.rs`
- Driver: `plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_slot_cost.py`

```bash
python3 plans/08162026_clumpify-memory-accounting/phase4/spikes/spike_slot_cost.py . 3
```

The driver stages the `.rs` into `examples/`, builds it, runs the matrix, and removes the staged file in
a `finally` — the spike leaves no source file behind. It refuses to start if
`examples/spike_slot_cost.rs` already exists rather than overwriting it.

## 3. Results per iteration

**#1 — first run, and the result was measured against the wrong allocator.** Slope method worked
immediately and reproducibly (0.0 B/pair spread within each read length), giving overhead 246–271 B/pair
and `s ≈ 132`. But `#[global_allocator]` is declared in `src/main.rs`, which is the *binary* — a cargo
example is its own binary and does not inherit it. So iteration 1 measured macOS libmalloc, not the
shipping allocator. Criteria not met: wrong configuration.

**#2 — declared mimalloc in the example.** `s` moved from 131.57 to 130.64; per-read-length overhead
from 246.0–270.7 to 245.0–269.6. **The allocator choice accounts for 1–4 B/pair, so it is not the
driver** — worth recording, because it was the leading hypothesis going in. Also split the reported
spread into run-to-run (within a read length) versus across-read-length, which the first version
conflated into one number and would have hidden the key finding.

**#3 — decomposed with a `spine` mode.** Rather than keep theorising about where ~100 B/pair came from,
added a mode that pushes records whose three `String`s are empty, so there is no heap text and the slope
is the `size_of` claim alone. This is the run that answers the question:

```
 len   spine  vs 148     full      2r  heaptext   alloc alloc/str        s
  25   148.0    +0.0    436.7   191.8     181.8   107.0      17.8   122.49
  36   148.0    +0.0    500.9   235.7     225.7   127.2      21.2   132.60
  51   148.0    +0.0    565.2   295.6     285.6   131.6      21.9   134.81
 100   148.0    +0.0    757.5   491.1     481.1   128.4      21.4   133.20
 150   148.0    -0.0    950.7   690.6     680.6   122.2      20.4   130.07

  spine-only slope: 148.0 to 148.0 B/pair against size_of's 148 -> CONFIRMED
  allocator cost per String: 17.8 to 21.9 B (6 Strings per pair)
  s = 122.49 to 134.81 B/record; run-to-run spread 0.02
  criterion (run-to-run <= 10 B/pair): MET
  size_of s = 74:    EXCLUDED
  matrix  s = 108.74: EXCLUDED
  single constant exact at every read length: NO — s varies by 12.3 B/record
```

## 4. Findings

1. **The spine is exactly 148.0 B/pair, at every read length, to the resolution of the measurement
   (0.0).** §2.3's `size_of` reasoning is exactly right about the structure it describes. It was never
   wrong — it was **incomplete**.

2. **Each of the six `String` allocations per pair costs 17.8–21.9 B beyond its text.** That is
   107–132 B/pair, and it is what both fits were groping at. Neither `size_of` (inline struct bytes) nor
   `estimated_record_bytes` (text length) can see it.

3. **`s` is 122.5–134.8 B/record. Both predictions are excluded** — 74 by 49–61 B, and the matrix's
   108.74 by 14–26 B.

4. **The excess is allocator-independent.** libmalloc and mimalloc agree within 1–4 B/pair. So this is
   not a mimalloc size-class artefact, and switching allocators will not recover it.

5. **`s` is not constant in read length** (12.3 B/record spread, ~10%). The per-allocation overhead is
   roughly flat (~20 B) while text grows, so the *text inflation factor* falls with read length —
   1.59× at 25 bp, 1.46× at 51 bp, 1.18× at 150 bp. The right shape is a **constant per record**, not a
   multiplier on text.

6. **This reconciles the matrix's 108.74 with the spike's 122–135.** The spike holds everything resident;
   in production the spine is resident in full from the first push after each flush, but the text — *and
   the allocator overhead attached to it* — oscillates over `0..B`. So the matrix's single `s` blends a
   fully-resident 148 B with a partially-resident ~120 B, and lands between the two. The two numbers are
   not in conflict; they are measuring different mixtures.

7. **The root cause of the breach, stated plainly.** The bin flushes when `raw_bytes` — a sum of
   `estimated_record_bytes`, i.e. *text* — reaches `B`. Actual resident is that text plus 148 B/pair of
   spine plus ~120 B/pair the allocator takes. At 25 bp the text term alone is understated by 59%.
   Short reads breach because they buy more allocations per byte.

## 5. Reference snippets

The spine measurement, which is the part worth carrying — an example needs its own allocator
declaration or it silently measures the platform's:

```rust
// An example is its own binary and does not inherit main.rs's allocator, so the
// shipping allocator has to be declared here or the spike measures libmalloc.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

Heap text is *not* `estimated_record_bytes` — that counts 5 B of `+` and newlines which exist on disk,
never on the heap. Any accounting that equates the two is 5 B/record out before it starts:

```python
heap_text = 2 * (r - 5)          # r = estimated_record_bytes, per record
allocator = (full_slope - spine_slope) - heap_text
```

## 6. Recommendation

**Re-plan §6, do not pick C or D on the current tables.** Both options' arithmetic is computed from
`s = 74`, so every cost figure in §6b and §6c is stale.

- **Option D as written is dead.** Its claim is `text + spine = A` *exactly, at every read length*, with
  a `size_of`-derived 148 B/slot. The spine is indeed exactly 148 — but it is not the whole per-slot
  cost, the remainder is ~120 B/pair, and that remainder varies by 10% across read lengths. Charging 148
  in the threshold leaves ~45% of the real per-slot cost uncharged. D can charge ~268 B instead, but then
  the constant is measured rather than derived, which removes D's headline advantage, and it is not exact
  at any read length but one.
- **Option C survives intact and gets more expensive.** C is a conservative bound, not a model, so a
  larger true per-record cost changes its numbers without touching its logic. Recomputing C's floor
  against `s ≈ 135` (the worst read length) rather than 74 is mechanical.
- **A third option now exists that neither reviewer proposed:** charge the *allocation count* rather than
  bytes. Resident is `text + 148·slots + ~20·(3·records)`, and the dispatcher already knows the record
  count it is about to reserve. This is D's shape — account rather than fit — but against the quantity
  that is actually measured to be constant (per-allocation overhead), not the one `size_of` happens to
  expose. It should be costed before §6 is re-opened.

**Next measurement, if wanted:** the occupancy question this spike deliberately left out — does the
allocator overhead ride the text's `0..B` oscillation, or does the allocator retain freed pages so it
behaves as fully resident? Finding 6 assumes the former, and it decides whether the charge is
`120·slots` or `a·120·slots`.

## 7. Limitations

- **The mechanism behind ~20 B per allocation is not attributed.** It is measured, reproducible to
  0.02 B, and allocator-independent, but it is not size-class rounding (both allocators agree, and the
  overhead is flat rather than tracking distance to a class boundary). Candidates not tested:
  16-byte quantum plus minimum block size, page-level accounting in `ru_maxrss`, or `mi_usable_size`
  differing from the request. This does not affect the verdict, which rests on the totals.
- **One platform.** darwin/arm64 only, with 16 KiB pages. Linux/x86_64 has 4 KiB pages and a different
  libmalloc; PLAN Assumption 6 already flags platform dependence.
- **Occupancy is out of scope**, as above — the spike measures fully-resident cost by construction.
- **Fixed-length records only.** Same blindness as the calibration matrix (PLAN §10 Finding 2): nothing
  here speaks to variable-length input, which is every real library after trimming.
- **`keys` is charged at 4 B/pair** as the code does. If a future minimizer key widens, this changes.
