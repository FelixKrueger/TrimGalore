# Plan review — Phase 4 (σ refit, #457) — Reviewer A

**Target:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v1
**Repo state:** `dev` at `c7dd106` (checked; no branch switch, no tracked file modified)
**Measurements run:** none. Every number below is either **code-verified** (read from `src/`,
`.github/`, `examples/`) or **arithmetic-verified** (integer arithmetic reproduced by hand/`awk`
against values the repo already pins). Claims I did not check are marked **unverified**.

**Arithmetic self-check.** Before trusting my own numbers I reproduced three values the repo owns
independently: the `1G/c4` pin `19_006_356` B, the `4 GiB/c32` pin `11_761_649` B, and the `--cores 2`
floor `281` MiB. All three match to the byte, so the integer pipeline I use below (`floor(budget×10/11)
− static`, `×16/denom`, `div_ceil` on the floor) is the same one the binary uses.

---

## Verdict in one paragraph

The diagnosis is right and the mechanism is real. But **§2.4's two constants are wrong by ~2.6× on
`s`, and the sign of the error is against the plan**: the true per-record cost is not a quantity that
needs a 72-run matrix to discover — it is `size_of::<FastqRecord>()` and it is already written into
`parallel.rs`, where each bin calls `reserve_exact(budget / bytes_of_first_record)` on every flush
cycle. From that alone, with zero free parameters, the observed 25-bp-vs-51-bp penalty is predicted to
within 4% on Linux and 7% on macOS. Substituting the code-derived `s` moves every cost figure in the
plan: the headline "~5–6% bin shrink" that §6 says "makes this option viable" becomes ~16%, `σ×16`
goes from 25 to 28–31, and the floor arithmetic in §6/§7a is wrong in a second, independent way. §7b's
load-bearing "misses by 0.001" does not survive being computed with the static term held fixed, and it
omits the second banded cell entirely. Option C is defensible, but it is being chosen against a cost
estimate that is 3× too favourable, and there is a fourth option — charge the spine the dispatcher
already reserved, in the flush threshold — that removes `r` from the model instead of fitting it, and
so deletes `r_min`, the warning, the bounded promise and the docs boundary at a stroke. That deserves
an explicit compare before Step 4 is written.

---

## 1. Logic

### 1.1 §2.3's form is right, and `s` is not an unknown — **Critical** (code-verified)

`src/parallel.rs:897-909` (`PairedBin::push`), and identically `:1286-1295` (`SingleBin::push`):

```rust
let bytes_this = estimated_record_bytes(&r1) + estimated_record_bytes(&r2);
if self.keys.capacity() == 0 {
    let predicted = self.budget.div_ceil(bytes_this).max(1);
    self.r1.reserve_exact(predicted);
    self.r2.reserve_exact(predicted);
    self.keys.reserve_exact(predicted);
}
```

`take()` (`:920`) does `mem::take` on all three vectors, so `keys.capacity()` is 0 again after every
flush and the reservation is re-made **once per bin per flush cycle**. Consequences:

- The spine is `reserve_exact`, so it is *exactly* `B/r` slots — not a doubling cascade, not an
  estimate. `FastqRecord` is three `String`s = **72 B**; `MinimizerKey` is `u32` = **4 B**. Paired:
  `72 + 72 + 4 = 148 B` per slot, and slots = `B/(2r)` because `bytes_this` is both mates. So
  **`s = 74 B` per record-byte-of-`r`** (paired), `76` (single-end). This is `σ(r) = a + 74/r`, with
  `a ≤ 1` rather than `≥ 1` — see §1.2.
- The spine is resident **in full from the first push after each flush**, while a bin's *text*
  oscillates over `0..B`. That is what makes `a < 1`: text is on average well below `B`, the spine
  never is.
- The same per-record cost rides the **worker** term: a batch in flight is a whole bin, text *and*
  spine. So `1/r` multiplies `n_bins` **and** `workers`, which no table in the plan reflects (§1.5).

**Zero-parameter check against the arc's own measurements.** Both `1G_c4` cells have the *same* `B`
(`B` does not depend on `r`), so their ratio difference isolates the per-record term:

```
predicted Δpeak = s × (n_bins + 2·workers) × B × (1/r25 − 1/r51)
                = 74 × (16 + 8) × 19,006,356 × (1/95.9 − 1/147.9)
                = 118.0 MB = 112.5 MiB
measured Δpeak  = (0.914 − 0.788) × 930.9 MiB = 117.3 MiB   (Linux, 27 obs/cell)
                = (1.07  − 0.939) × 930.9 MiB = 121.9 MiB   (macOS, phase3 §"First runner results")
```

Within 4% and 7% respectively, with **no fitted constant** — and the `+2·workers` factor is not
fudged: `k = 4` is documented as "queued batch + batch in hand + 2 result slots", and exactly the
first two carry *records*; the result slots carry compressed bytes. The mechanism is therefore
established analytically, and §4's matrix should be re-scoped from *discovering* `s` to *confirming*
it (§4 below).

**What this does to §2.4.** Fitting the two points the plan uses gives `s = 28.87` (the plan says
28.7 — recompute: `0.1125 / (1/92.4 − 1/144.4) = 0.1125/0.00389730 = 28.87`), i.e. **2.6× below the
structural value**. The reason is the choice of points, not the algebra:

- `1.4375` is not "σ fitted at 51 bp". It is `23/16`, the **rounded-up envelope** of Phase 2's
  measured `1.29` (Linux) / `1.41` (macOS) (`src/clump.rs:281-287`). Mixing a safety-rounded model
  constant with a raw measurement at the other end of the axis biases the slope toward zero.
- `1.55` at 25 bp cannot be reconciled with phase 3's own measured ratios. Back-solving σ from
  macOS `25bp_1G_c4 = 1.07` of prediction (static 205 MiB, `B = 18.126` MiB, `k = 4`) gives
  **σ ≈ 1.73**; from the c2 breach (1.02–1.06× *budget* = 1.12–1.17× prediction) gives **≈1.86**. Using
  1.73 with 1.4375 gives `s ≈ 80, a ≈ 0.88` — which is the structural answer (74, ≲1) to within 8%.

**Action:** replace §2.4 with the structural derivation, `s = 148 B/pair` from `reserve_exact` +
`size_of`, and pin it with `assert_eq!(size_of::<FastqRecord>(), 72)` next to `estimated_record_bytes`
so the model's constant and the struct cannot drift apart silently.

### 1.2 `a` is not "bytes held more than once" — **Important** (code-verified + arithmetic)

§2.4 reads `a ≈ 1.238` as "the pool appears to hold ~24% more than one copy of each bin's bytes
*before* any per-record spine". That interpretation is not available:

- It is arithmetically inconsistent with its own `s`. A 72-byte spine per record is a hard floor on the
  per-record term; `s = 28.7` is *below* it. Any (a, s) pair with `s < 72` is telling you the fit is
  absorbing something else, not measuring copies.
- The mechanism points the other way. Bins flush at `raw_bytes ≥ budget`, so instantaneous text
  occupancy per bin is ~uniform on `(0, B)` and the peak of the sum over 16 near-independent bins is
  ≈`0.6·16·B`, not `1.0·16·B`. Solving the structural form against the shipped constant at the fitting
  fixture's `r` (147.9 B at 10 M pairs) gives `a = 1.4375 − 74/147.9 = 0.937` — i.e. **~94% of one
  copy**, the rest of `k`'s and σ's magnitude being in-flight duplication already charged separately.
- σ and `k` are entangled in the Phase 2 fit (`n_bins` pinned at 16 while `cores` varied), so the
  split between them is a fit artefact. Reading either as a physical count over-claims.

Keep the plan's own caveat ("order-of-magnitude sanity check only") and delete the interpretive
sentence; report `a` as an occupancy factor, and expect the fit to **exclude 1.0 from below**, which
is the opposite of what §2.4 primes the implementer to look for. Validation 2 as written ("report `a`
with its interval and whether 1.0 is excluded") is fine — it is the framing above it that will make a
correct `a < 1` look like a failed fit.

### 1.3 §7a's `r` values are measured on the wrong fixture size — **Important** (arithmetic-verified)

`examples/mk_clump_fixture.rs:57` writes `@SIM:1:FCX:1:15:6329:{i} {mate}:N:0:ATCCGA`, so
`id.len() = 34 + digits(i)` and `estimated_record_bytes = id + 2L + 5`.

| pairs | mean `digits(i)` | `r` at 25 bp | `r` at 51 bp |
|---|---|---|---|
| 2,000 (the plan's table) | 3.445 | **92.4** | **144.4** |
| 10,000,000 (§4's and CI's fixture) | 6.889 | **95.9** | **147.9** |

The plan's table reproduces exactly at 2000 pairs (92.445, 144.445 — good, its provenance is sound),
but §4 mandates 10 M-pair fixtures and the CI job already uses them, where `r` is 3.5 B higher.
§2.3's own rule ("any fixture used for fitting must have its `r` measured rather than assumed") is
violated by its own table. The effect on `σ×16` is small (24.6 → 25 at `r_min = 92`), but it flows
into `r_min`, which is *the* remaining decision, and the header is 34–41 B where a real Illumina
header is ~50–60 B — so the "short headers are conservative" claim holds, just with less margin than
the table implies.

### 1.4 §7b does not survive being computed — **Critical** (arithmetic-verified)

Two independent problems.

**(a) It scales the whole ratio, when only the pool scales.** `ratio = peak/prediction`, and
`peak = static_actual + (terms proportional to B)`. `static_actual` on Linux is 173 MiB
(`src/clump.rs:250-252`) against a 930.9 MiB prediction, so **18.6% of every ratio does not move when
`B` moves**. The plan multiplies the whole ratio by the shrink. Recomputing with the static term held
fixed (and, in the third column, also scaling lever A itself — a work-channel-depth revert adds
`k·B` per worker, so its magnitude is proportional to `B` too):

| `σ×16` | `B` shrink (c4) | plan's drop | drop, static held | new mean | plan's verdict | lever A fixed at +0.126 | lever A scaled |
|---:|---:|---:|---:|---:|---|---|---|
| 25 | 4.878% | 0.037 | **0.029** | 0.759 | 0.877 — misses by 0.001 | 0.885 — **catches by 0.007** | 0.879 — catches by 0.001 |
| 26 | 7.143% | 0.054 | **0.043** | 0.745 | 0.860 — misses by 0.018 | 0.871 — misses by 0.007 | 0.866 — misses by 0.012 |
| 28 | 11.364% | 0.086 | **0.068** | 0.720 | 0.828 — misses by 0.050 | 0.846 — misses by 0.032 | 0.838 — misses by 0.040 |

At `σ×16 = 25` — the plan's leading candidate — the conclusion **inverts**. At 26 and 28 it survives.
So the claim "at *every* candidate the band would go green on the very regression it exists for" is
not established; what is established is that the outcome at the leading candidate sits inside ±0.01,
i.e. inside the modelling error, which is a *different* and weaker statement.

**(b) It omits half the band.** `BANDED="1G_c4 25bp_1G_c4"` (`ci.yml:351`), and the job fails if
**either** cell fires. §7b analyses only `1G_c4`. For `25bp_1G_c4` (baseline 0.914, lever A measured
+0.128 per the TSV's "fires by 0.038", trigger 1.004):

| `σ×16` | new mean (static held) | + lever A (0.128) | verdict |
|---:|---:|---:|---|
| 25 | 0.878 | **1.006** | catches by 0.003 (0.999 → misses by 0.005 if lever A scales) |
| 26 | 0.862 | 0.990 | misses by 0.014 |
| 28 | 0.831 | 0.959 | misses by 0.045 |

Same shape: marginal at 25, clear at 26+. Since the two cells are near-independent draws, "the band
misses" needs *both* to miss, which roughly squares the risk the single-cell analysis reports.

**What this changes.** Not the direction — a band whose detection margin is 0.001–0.007 under one
modelling assumption and −0.007 under another is not a band anyone should rely on, so the plan's
*decision* (do not carry the old centres forward) is right. What changes is the justification: it must
be "the new centres cannot be predicted to better than ±0.01 and the margin is that size", not "the
band misses by 0.001". §9 explicitly praises this table for turning a caution into a hard constraint
by computing it — that is exactly the move this arc has been burned by three times, and it happened
again here.

### 1.5 A single σ cannot cover c2 and c4 at short `r` — **Important** (arithmetic-verified)

Because the `1/r` term multiplies `n_bins + 2·workers` (§1.1) and `n_bins = max(16, 4·cores)`, the
`1/r` weight per unit of denominator is *not* constant across core counts: at c4 it is 24 units
against 16 bins, at c2 it is 20 units against 16 bins. Folding all of it into σ (which multiplies
`n_bins` only) therefore needs a **different** σ at c2 than at c4. Computing the required
denominator at `r_min = 92` with `a = 0.937, s = 74, k_text = 3`:

| cell | current denom | required denom | implied `σ×16` | `B` shrink |
|---|---:|---:|---:|---:|
| c4, 16 bins, 4 workers | 624 | 740.6 | 30.3 → **31** | **15.7%** |
| c2, 16 bins, 2 workers | 496 | 593.2 | 29.1 → **30** | **16.4%** |
| `--clump_only` (16 bins, 1 worker) | 432 | 519.4 | 28.5 → **29** | **16.8%** |

Three things follow. **(i)** The honest σ at `r_min = 92` is 28–31, not 25 — so §7a's table is
optimistic by 3–6 integer steps and the shrink is ~16%, not 4.9–6.1%. **(ii)** Taking the max (31)
over-charges c2 by ~5%, which is the price of keeping one constant; that is acceptable but must be
stated. **(iii)** `--clump_only` (workers = 1) is the **worst-affected** mode at every candidate
(6.9% / 10.0% / 15.6% at σ×16 = 25/26/28 respectively, versus c4's 4.9/7.1/11.4), and it appears in
no cost table in the plan even though §9 names it as an edge case and Step 5 edits its docs.

### 1.6 The floor arithmetic in §6 and §7a is wrong — **Important** (arithmetic-verified)

§6 Option B: "Floors rise by the same ~5% … the `--cores 2` floor moves 281 → ~290 MiB." Both halves
are wrong, because `clumpify_min_memory_bytes` (`clump.rs:243-246`) is
`ceil((static + ceil(MIN_BIN_BYTES·denom/16)) · 11/10)` and **`static` (224 MiB) dominates the 31 MiB
dynamic part** — only the dynamic part scales.

| `σ×16` | `dyn_min` | c2 floor | vs today |
|---:|---:|---:|---:|
| 23 (today) | 31 MiB | **281 MiB** (pinned in `clumpify_floor_at_the_default_cell_is_pinned_absolutely`) | — |
| 25 | 33 MiB | **283 MiB** | +0.7% |
| 26 | 34 MiB | **284 MiB** | +1.1% |
| 28 | 36 MiB | **286 MiB** | +1.8% |

`290 MiB` is today's **c4** floor, which is where the number appears to have come from. Note the
direction: this one overstates the cost of the option the plan chose, so it did not bias the decision
— but Step 5's "the floor range … move[s] ~5–14%" is wrong in the docs-facing direction and would put
a wrong figure into `clumpy.md` if the implementer scales rather than re-derives (which is precisely
the risk §9(f) names).

### 1.7 The chosen candidate probably does not close #457 — **Critical** (arithmetic-verified from committed data)

The plan already has the data to pre-screen `r_min` and does not use it. `.github/clumpify-memory-baseline.tsv`
records `25bp_1G_c2` at **n=25, mean 1.030, min 1.007, max 1.050** — the #457 cell. The promise is
gated per-run, so what must clear 1.0 is the **upper tail**, not the mean:

| `σ×16` | c2 `B` shrink | new mean | new max | margin at max |
|---:|---:|---:|---:|---:|
| 25 | 6.06% | 0.979 | **0.999** | 0.001 |
| 26 | 8.82% | 0.956 | 0.976 | 0.024 |
| 28 | 13.89% | 0.913 | 0.933 | 0.067 |

(static held at 173 MiB; the plan's whole-ratio scaling gives 0.986/0.960/0.905 at max — still under
5 points at `σ×16 = 25`.) And this is the *Linux* cell; macOS runs ~17 points hotter at 51 bp (0.939 vs
0.770, per `SESSION_HANDOFF.md` §5) and phase 3 measured 25 bp/c2 at **1.02–1.06× budget** = 1.12–1.17×
prediction, which needs ≥18% of pool cut just to reach 1.000 — more than **any** candidate in §7a
delivers. §6's Option B caveat gestures at this ("may not bring 25 bp back under the printed
prediction"); §7a then presents `r_min = 92 B` as the natural first row anyway, and §6's decisive
sentence — "**A ~5–6% bin shrink is a much smaller price than §2.5 led me to expect**, and it is the
number that makes this option viable rather than merely simple" — is anchored on the one candidate
the committed data says is insufficient. §2.5's instinct was right; the revised estimate was not.

**Action:** add a pre-screen step *before* Step 2 that solves the required shrink from the recorded
upper tails on both platforms. It costs nothing, it narrows `r_min` to one or two candidates, and it
prevents the 72-run matrix from being spent around the wrong point.

### 1.8 §3.1 contradicts the chosen option; §3.2's premise is moot — **Important**

§3 item 1 requires that "`dyn_denominator` becomes a function of `r`". Option C does **not** do that —
it changes one integer. Either §3 is pre-decision text that should be re-scoped to "the model charges
the spine correctly at the documented floor", or Option C does not satisfy the stated behaviour. As
written, a plan-manager coverage pass will report §3.1 as unimplemented. Same for §3.2's "a rational
`s/r` term must keep the same overflow-safety and exact-inverse properties" — under C there is no
rational term, so the whole worry evaporates (which is worth saying explicitly, since it is one of
Option C's genuine advantages).

### 1.9 Validation 6's threshold is guaranteed to fail — **Important** (arithmetic-verified)

"agreement within the current ~312 B across cores ≤ 32". The exact-inverse error is bounded by
`denom/16 + 1`, and `312` is precisely `dyn_denominator(128, 32)/16 = 4992/16`. Raising σ raises the
denominator, so the bound becomes **328 B** (σ×16 = 25), **336** (26), **352** (28). Stated as a fixed
312 B, Validation 6 fails at every candidate for a reason that has nothing to do with correctness.
Restate as `≤ dyn_denominator(n_bins, cores)/16 + 1`, which is provable rather than magic, and note
that `floor_and_sizing_round_consistently_without_overflow` already sweeps `cores` to 65536 (bound
≈624 KiB there) while asserting only `predicted ≤ budget`. **The property itself is safe**: σ×16 stays
an integer, so `saturating_mul` and the u64 headroom are unchanged (largest product at
`budget = u64::MAX` is unchanged in magnitude).

### 1.10 The CI blast radius is incompletely enumerated — **Important** (code-verified)

Step 4.6 says "recompute all six cells of the 4a pin table and the absolute floor pin. Every one
moves." Verified — and these also move, unnamed:

| site | what breaks at σ×16 = 25 |
|---|---|
| `clump.rs:556-561` `resolve_layout_default` | asserts `112 MiB ≤ B < 114`; B goes 112.89 → **106.05** MiB |
| `clump.rs:564-570` `resolve_layout_scales_with_cores` | asserts `44 ≤ B < 46`; B goes 44.87 → **42.68** MiB |
| `clump.rs:573-580` `resolve_layout_pins_the_default_cell` | pins `19_006_356` → 18_079_217 (a 7th pinned value, outside the "six cells") |
| `clump.rs:657-675` `resolve_layout_fastqc_reserves_more_and_shrinks_bins` | hard-codes `/ 624` in `expected` |
| `clump.rs:171-172`, `:267-289` doc comments | "σ = 23/16", the written-out `16 × usable / (23 × n_bins + 64 × workers)`, and the whole Phase-2 provenance paragraph |
| `ci.yml:441-444` | the comment "25 bp is above the prediction today … Recorded, not gated on the promise, **until that refit lands**" plus `measure 25bp_1G_c4 **0**` |
| `.github/clumpify-memory-baseline.tsv` comment block | 27-observation stats for all **six** cells, plus the `25bp_1G_c2` note "its value will move when sigma is refitted" |

The `ci.yml` row is the sharpest: **the plan never flips `gating` from 0 to 1 on the 25 bp cells.**
That flag is the only CI-visible proof that #457 is closed, the in-tree comment promises it will
happen when the refit lands, and Validation 8 only asks that the cell "falls by the predicted amount".
If `25bp_1G_c4` is not promise-gated after the refit, the fix ships unenforced. (Whether `25bp_1G_c2`
can also be gated depends on §1.7's margin — which is an argument for the larger σ, not for leaving it
ungated.)

### 1.11 The empty-`BANDED` mechanism already exists — **Important** (code-verified)

§7b's diagnosis is correct: with `BANDED=""` the armed-check loop at `ci.yml:363-366` iterates zero
times and reports success, and with the rows deleted the per-cell `[ -n "$want" ]` guard skips the
band silently. Verified.

But the remedy already exists and is tested: `ci.yml:353-357` greps `^# BOOTSTRAP` in the baseline
file and, when present, prints `::warning::baseline is in BOOTSTRAP mode - drift band not enforced`,
skips the armed-check, skips the band, and tells the reader how to re-arm — while the promise gate at
`:420-423` keeps running (so §7b's "the promise assertion keeps gating throughout" is code-verified).
Phase 3 shipped in exactly this state on purpose. So Step 4.8 should be **"add the `# BOOTSTRAP`
line"**, not "empty `BANDED` and add a new empty-`BANDED` warning": no new CI code, no deleted rows,
the ratios still print, and the re-arm instruction is already in the warning text. Writing a second,
parallel disarm path is the more expensive and less tested of the two.

### 1.12 Option C's warning site does not do what §6 says — **Important** (code-verified)

§6 sells C as "a warning from the record already read" — reusing `sanity_check_with`'s discarded
record. Step 4.7 instead samples afresh inside `resolve_clump_layout`. Those are different designs,
and the second has costs §6's "For" list does not carry:

- `resolve_clump_layout(cli: &Cli)` (`main.rs:178`) has no reader; it would have to **open and
  decompress `cli.input[0]` a second time**, adding an I/O + error path to a function whose only
  current failure is `parse_memory_size`. Ordering is fine (`sanity_check_any` at `:270` precedes all
  three call sites at `:1386`, `:1623`, `:1902`), so a failure there is a *second* failure of an input
  already validated — decide whether it warns-and-continues or propagates, and say so.
- Which input? It is called **per pair** (the CI job's `head -1` comment exists because the banner
  prints once per pair), so a multi-pair run should sample each pair; and in paired mode the pool
  holds R1 **and** R2, whose lengths can differ (already-trimmed input, 2×different-length runs). `r`
  for the memory model is the mean over both mates, not R1's.
- "Sample a fixed number of records (deterministic, so the message is stable)" — deterministic per
  input, yes; but the honest statistic is not the mean. The dispatcher's reservation keys off the
  **first record after each flush** (§1.1), so the quantity that actually drives memory is an order
  statistic over ~80 draws per bin, not the mean. A mean-based warning will stay silent on exactly the
  input that hurts (see §1.13).

None of this sinks Option C. But §6's "the sample influences only a message" is true only of the
sizing; the *reading* is new behaviour on the startup path and should be listed as such.

### 1.13 An unmodelled hazard the refit does not touch — **Important** (code-verified, magnitude unverified)

`predicted = self.budget.div_ceil(bytes_this)` is uncapped and keys off **one** record. If that record
is shorter than the bin's mean — an adapter-dimer read at the head of a flush cycle in a 150 bp
library, say — the bin over-reserves the spine by `r_mean/r_first`. There are ~80 flush cycles × 16
bins per 10 M-pair run, so an unusual first record is not a rare event; the over-reservation is
stochastic, input-dependent, and charged nowhere. It is a plausible contributor to the 0.043–0.081
within-cell ranges the TSV records, and it is invisible to §4's matrix because **every**
`mk_clump_fixture` output is constant-length. A running mean of the previous cycle's record size
(free, deterministic per input) or a clamp would remove it. Worth one line in the plan even if it is
deferred — it is the same defect class as #457 (a per-record cost sized from the wrong quantity) and
fitting σ on constant-length fixtures will not see it.

---

## 2. Assumptions

- **Assumption 3 (integer arithmetic / exact inverse) — safe**, see §1.9. Verified against
  `predicted_peak_bytes` and the u64 sweep.
- **Assumption 2 — verified**, though the reference is `cli.rs:743`, not `:747`.
- **Assumption 5 (one `(a, s)` for both platforms) — open, and the structural result makes it
  *likelier* than the plan expects.** `s` is `size_of` — it is platform-invariant by construction on
  any 64-bit target. What differs between platforms is `a` (occupancy/allocator retention) and
  `static`. So the fit should be structured as "shared `s`, per-platform `a`", which is a stronger and
  cheaper hypothesis than fitting both per platform.
- **Assumption 6 (linearity in record count) — the plan tests for curvature; it should also expect
  non-monotone scatter.** `reserve_exact` removes the `Vec` doubling sawtooth I would otherwise have
  flagged here (good — the code is better than the plan's Assumption 6 fears), but three heap
  allocations per record still round to mimalloc size classes, so `s` has a small step component in
  read length. Residuals may scatter rather than curve; Validation 2's pass/fail needs a tolerance, not
  a shape verdict.
- **§2.3's "header length moves σ exactly as read length does" — true, and the control cell is
  right.** But it is not a *pure* `r`-mover: lengthening the id changes the id allocation's size class
  too. So Validation 3 / Step 3's binary stop condition ("lands off the fit ⇒ the hypothesis is wrong
  ⇒ stop and re-plan") can fire on a second-order effect. State the tolerance up front — with
  `s = 74` known a priori, "on the fit" means "within ~5% of the predicted `74·(n+2w)·B/r` delta",
  which is a real test rather than a judgement call.
- **Unstated assumption, now the load-bearing one: that measured peak is a deterministic function of
  `B`.** It is not — the promise is a fit to a distribution with an observed range of 0.034–0.081 per
  cell, plus the stochastic term in §1.13. Validation 5's "≥3 reps" is thin against that; the arc's
  own baseline file argues a handful of draws estimates dispersion low by construction. Size `r_min`
  against the recorded upper tail (§1.7), not against 3 fresh reps.

---

## 3. Efficiency and cost of the plan itself

- **§4's run count is wrong: 6 lengths × 2 core counts × 3 reps = 36 per platform**, 72 across both.
  The plan says "72 runs per platform". Add the falsification (12/platform) and the header control
  (3) and it is ~51 per platform, each a 10 M-pair paired run. **Unverified** but worth stating: at
  a few minutes per run that is multiple hours per platform, and the `clumpify-memory` job carries
  `timeout-minutes: 30` (`ci.yml:277`) — so Step 2's "throwaway draft-PR branch" route cannot host
  this matrix without splitting it. Phase 2 got away with it on a smaller grid.
- **Six 10 M-pair fixtures at 25–150 bp is several GB of gzip** on top of per-run outputs. §4's
  disk-assertion rule covers correctness; the plan should still say the number so nobody discovers it
  mid-matrix (the arc has one false finding from disk exhaustion already).
- **The matrix can be much smaller.** With `s` known (§1.1), the free parameters are `a` and
  `k_text` — 2 unknowns, both of which Phase 2 already fitted at one length. Three read lengths
  (25/51/150) × 2 core counts × 3 reps = 18/platform *verifies* the form and refits `a`; six lengths
  buys a linearity check that the structural derivation makes near-tautological. Spend the saved runs
  on the axes that are genuinely unmeasured: **variable-length input** (§1.13) and **c2 on macOS at
  25 bp with enough reps to see the tail** (§1.7).
- **Step 1's harness repair list is incomplete.** `calibrate.py` is stale well beyond the banner
  regex (which the plan correctly verified — `:35` cannot match the current banner):
  `STATIC_MIB = 512` (`:34`) against the shipped 224 MiB, and `denom_for(cores) = 5*n_bins + 7*cores`
  with `/4` scaling (`:52-61`) — that is σ = 1.25, k = 1.75, i.e. the **pre-Phase-2** model. Those
  helpers pick each cell's `--memory` (`budget_mib_for`) and report `b_exact`, so repairing only the
  regex leaves the grid placed by a two-generations-old model and the printed "shipped" comparisons
  meaningless. Same failure class §2.6 is about, one function down.

---

## 4. Alternatives the plan does not consider

### D. Charge the spine in the flush threshold — removes `r` from the model instead of fitting it

The dispatcher already computes the exact spine it is about to reserve
(`148 × budget.div_ceil(bytes_this)`). Size the reservation from an **accounted** threshold instead of
a text threshold: `slots = A / (bytes_this + 148)`, flush at `raw_bytes ≥ A − 148·slots`. Then
`text + spine = A` per bin **at every read length**, and:

- σ becomes genuinely `r`-invariant (actual/A varies 0.95→0.97 across 25–150 bp — a 2% residual,
  versus the 13–17 point spread that is #457). **No `r_min`, no warning, no bounded promise, no docs
  boundary, no Assumption 7, and §7a's remaining decision disappears.**
- Re-deriving both coefficients in accounted units (σ ≈ 1.0, and `k` drops because an in-flight batch
  is now `A` inclusive of its spine — roughly 2.4A rather than 4B) puts the 51 bp text budget at
  ~19 MiB against today's 18.1: **long-read members do not shrink at all**, and 25 bp pays ~11%. Under
  Option C at the honest σ (§1.5) *everyone* pays ~16%. D is cheaper at every read length.
- The band story improves rather than degrading: if `A` is chosen so the 51 bp text budget is
  unchanged, `1G_c4`'s peak and prediction are both unchanged, so **its baseline row survives** and
  only the 25 bp row needs re-basing — and its new value is *predictable* (it must equal the 51 bp
  cell's, by construction of `r`-invariance), which is exactly what §7b says cannot be done for
  Option C.

Costs, honestly: ~6 lines across `PairedBin::push`/`SingleBin::push` plus one field, in a hot path;
output bytes move at every read length (as they do under C); σ's *meaning* changes, so the docs
rewrite is larger than Step 5's; and the coefficient re-derivation needs the same matrix. My
arithmetic here is a first pass — the 19 MiB/11% figures are **unverified** and would need the fit to
confirm. But D dominates C on the two axes the plan says it cares about (compression cost, and how
much of the promise stays bounded), so it deserves a row in §6 rather than silence.

### E. Cap the reservation (§1.13) — cheap, independent, and fixes a defect no option addresses

Derive `predicted` from a running mean of the previous cycle rather than one record, or clamp it at
`A/r_min`. Deterministic, no CLI surface, no model change. Worth doing under any option.

### F. Keep the band armed by re-centring, rather than disarming it

The baseline file's double-digit-sample rule governs the **width** ("the half-range must come from a
double-digit sample"), because width is what noise fires. The **centre** is a mean, and its error
enters the trigger linearly: post-refit dispersion should be near-unchanged (same runner, same
fixture, `B` a few percent smaller), so keeping `width = 0.09` and taking the centre from the refit
PR's own 2–4 runs costs ~0.01–0.02 of centre error against a lever-A margin of 0.036–0.038. That
keeps both cells armed continuously, deletes Step 6 (PR C) and the entire inert-but-looks-armed
hazard, and the follow-up becomes "re-mean at ≥10 observations" rather than "re-arm". The plan should
say why this is not preferred; as written it reads as if the double-digit rule forbids a provisional
centre, and it does not.

---

## 5. Validation sufficiency

| # | Verdict |
|---|---|
| 1 | Good, and cheap. Extend it: break `denom_for`/`STATIC_MIB` too, not just the regex (§3). |
| 2 | Needs a tolerance and a stated expectation that `a < 1` (§1.2, §2). As written a correct fit reads as a failure. |
| 3 | Right cell, wrong stop rule — make it quantitative against `74·(n+2w)·B/r` (§2). |
| 4 | Fine. But `cores ∈ {8,16}` with `n_bins = 4·cores` holds the `1/r` weight **constant**; it cannot falsify the σ-vs-`k` mis-attribution. Only **c2** (where `n_bins` floors at 16) can, and c2 is a *fitted* cell, so nothing holds it out. Move one c2 cell to the falsification set (§1.5). |
| 5 | Under-powered: 3 reps against a cell whose recorded range is 0.043 (§1.7). Use the recorded tails to set `r_min`, then confirm. |
| 6 | Threshold wrong (§1.9). |
| 7 | Good, and the one item with no substitute. |
| 8 | Incomplete — the two banded cells, not one, and the `gating` flags must flip (§1.10). |
| 9 | Fine. `--clumpify` is on no Perl-parity path; asserting rather than assuming is right. |
| **missing** | Nothing tests that the **printed prediction is unchanged**. It is the plan's own load-bearing invariant (§2.5, §7 A1) and three integration tests already assert `predicted peak ≈ 930 of 1024 MiB budget` (`integration_preflight_tripwire.rs:605`, `integration_clump_only.rs:126`, `integration_clump_only_ubam.rs:169`). Name them as the check — they turn "assert rather than assume" into a green/red. |
| **missing** | Nothing measures a **variable-length** input (§1.13), so the refit's constants are verified only on the one input shape that cannot expose the reservation's actual behaviour. |

---

## 6. Action items

### Critical

1. **Replace §2.4's constants with the structural derivation.** `s = 148 B/pair` (`72·2 + 4`) from
   `reserve_exact` in `parallel.rs:897-909`; `a ≈ 0.94` from `1.4375 − 74/147.9`. The two-point fit
   gives `s = 28.87` (not 28.7) and it is ~2.6× low. Verify by the zero-parameter prediction of the
   25-vs-51 bp penalty: 112.5 MiB predicted vs 117.3 (Linux) / 121.9 (macOS) measured. Add a
   `size_of::<FastqRecord>()` pin so the constant and the struct cannot drift.
2. **Recompute every cost figure that depends on `s`** — §6's Option B table, §7a's σ and shrink
   columns, §7b's drops, §6's headline "~5–6%". At `r_min = 92 B` the honest σ×16 is 28–31 and the
   shrink is ~16% at c4/c2 and ~17% on `--clump_only`. Then re-read §6's decision: it was made against
   a 3×-too-favourable number.
3. **Pre-screen `r_min` against the recorded upper tails before spending the matrix.** `25bp_1G_c2`
   max 1.050 over 25 observations, plus macOS 1.12–1.17× prediction, says `σ×16 = 25` leaves 0.001 of
   margin on Linux and does not clear macOS at all. Free, and it decides the only open sub-decision.
4. **Rebuild §7b with the static term held fixed and both banded cells.** At `σ×16 = 25` the band
   *catches* by 0.007 rather than missing by 0.001; at 26/28 the conclusion survives. Re-argue the
   sequencing from "the new centres cannot be predicted to better than ±0.01, and the margin is that
   size" — which is both true and sufficient.

### Important

5. Use the existing **`# BOOTSTRAP`** marker instead of emptying `BANDED` plus new warning code
   (§1.11). Consider alternative F (provisional centre, band stays armed) and say why not.
6. **Flip `gating` 0 → 1 on the 25 bp CI cells** and delete `ci.yml:441-443`'s "until that refit
   lands" comment. Without this the fix has no CI enforcement.
7. Add the unnamed blast-radius sites from §1.10 to Step 4 — four more tests, two doc comment blocks,
   the TSV's six-cell statistics block.
8. Fix the floor arithmetic (§1.6): 281 → 283/284/286 MiB, ~1% not ~5%; `290` is today's c4 figure.
   Correct Step 5's "~5–14%" for floors.
9. Restate Validation 6 as `≤ dyn_denominator/16 + 1` (§1.9); note that Option C makes §3.2's
   rational-term worry moot, and re-scope §3.1 so Option C can satisfy it (§1.8).
10. Give the `1/r` term its own coefficient on the **worker** side, and hold out a **c2** cell
    (§1.5, Validation 4). One σ cannot cover c2 and c4 at short `r`; the plan must choose the max and
    say what it over-charges.
11. Extend Step 1's harness repair to `STATIC_MIB` and `denom_for` (§3), and reconcile §2.3's `r`
    table with the 10 M-pair fixture (95.9 / 147.9, not 92.4 / 144.4) (§1.3).
12. Reconcile §6's "warning from the record already read" with Step 4.7's fresh sample; decide the
    R2 / multi-pair / failure-to-open behaviour (§1.12).
13. Add **alternative D** to §6 as a compared option (§4).

### Optional

14. State §4's real run count (36/platform, ~51 with falsification) and its disk and wall-clock cost;
    check it against the job's 30-minute timeout if the draft-PR route is used.
15. Shrink the fitted grid to 3 lengths and spend the difference on variable-length input and macOS
    c2 reps (§3).
16. Note the uncapped `predicted` hazard (§1.13, alternative E) even if deferred.
17. Fix the small reference drifts: `cli.rs:743` (not `:747`), `calibrate.py:64-79` (not `:65-73`).
    §9's claim that no user-facing string cites an issue number is **verified** — the `#NNN` hits in
    `src/` are all test assertion messages.

---

## 7. What I did not check

- `#446`'s `+0.27% / +0.24%` output-size figures, and the CHANGELOG/docs claims about them.
- Phase 2's `results_cal2_linux.tsv` and `CALIBRATION.md` beyond `calibrate.py`'s constants; I did
  not re-derive Phase 2's σ/k fit.
- Anything requiring a `trim_galore` invocation. I ran none, so all peak-RSS values here are read from
  committed artefacts (`.github/clumpify-memory-baseline.tsv`, `phase3/PLAN.md`, `src/clump.rs`
  comments), and my structural predictions are arithmetic on top of them.
- The 19 MiB / 11% cost figures for alternative D — first-pass arithmetic, explicitly unverified.
- Whether mimalloc's size-class rounding materially inflates `s` above 148 B/pair; I assert only that
  148 is a floor and that it already predicts the measured penalty to within 4–7%.
