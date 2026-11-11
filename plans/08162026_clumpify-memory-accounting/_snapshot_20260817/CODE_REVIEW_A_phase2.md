# Code Review A — Phase 2, `--clumpify` memory accounting (#439)

**Reviewer:** A (independent)
**Date:** 2026-08-17
**Target:** uncommitted diff on `fix/439-clumpify-peak-memory`, HEAD `70e75e3`
**Touched:** `src/clump.rs`, `src/clump_only.rs`, `src/main.rs`, `src/parallel.rs`, `CHANGELOG.md`
**Inputs read:** `phase2/CALIBRATION.md`, `PLAN.md` (Implementation notes — Phase 2), `phase2/results_cal_darwin.tsv`, `phase2/results_cal2_linux.tsv`

## Verdict

The core arithmetic is **correct and better than the previous model**. Both properties the
change stakes itself on hold up under brute force: the predictor never exceeds the budget, and
the advertised floor is byte-exactly the smallest resolvable budget. Those are real results, and
the floor exactness is load-bearing rather than cosmetic (see V2).

The problems are all on the **FastQC term**, which is the one part of the model that
`CALIBRATION.md` explicitly said was unmeasured and explicitly recommended building a
*different way*. It became a per-core addend to the static reservation; the calibration spec
said `max`, not addend. That single choice produces one silent feature regression at the default
budget, one hard-abort regression on `--clump_only`, and one path where the model
under-reserves by ~670 MiB. None of it is caught by the tests, because the per-core scaling has
no test that can see it.

Counts: 3 High, 5 Medium, 7 Low. No Critical — no path corrupts output; the failures are a
silently disabled feature, an abort, and a breachable promise.

## Verification performed

| check | result |
|---|---|
| `cargo test --release` | **732 passed, 0 failed**, exit 0 (matches recorded state) |
| `cargo clippy --all-targets --release -- -D warnings` | clean, exit 0 |
| `cargo fmt --all -- --check` | clean, exit 0 |

The first clippy run finished in 0.46 s — a cache hit, not a compile. I `touch`ed the four
changed sources and re-ran; the second run emitted `Checking trim-galore v2.3.0` and took
2.18 s, so exit 0 reflects a real check of the current tree.

### V1 — `predicted_peak_bytes()` is a sound inverse (no finding)

I modelled both functions exactly (u64, `saturating_*`, floor division) and swept
`budget × cores × fastqc`: every per-cell floor for `cores` 1–64 both fastqc settings, every MiB
from 1–4096 MiB, GiB powers up to 64 TiB, and `2^60 … 2^63`, `U64`, `U64/2`, `U64/10`,
`U64/11`, `U64/16`, crossed with `cores ∈ 1..64, 128, 256, 1024, 65536, 2^20`.

**Zero violations of `predicted ≤ budget`. Zero overflows of `B × denom`.**

It is provable, not just empirical. `B = ⌊usable·16/denom⌋ ⟹ B·denom ≤ usable·16 ⟹
dyn = ⌊B·denom/16⌋ ≤ usable`, so `peak = static + dyn ≤ static + usable = budgeted ≤ budget`.
This holds for *any* σ, k, static and any margin ≤ 1 — which is exactly why the test that
asserts it is vacuous (M1).

Overflow is safe in the right direction. `saturating_mul` caps *below* the true product, so a
saturated `usable·16` yields a **smaller** `B` and therefore a smaller prediction. And
`B × denom` can never overflow when `B` came from a division by `denom`: `B·denom ≤ X ≤ U64`.
The only way to overflow it is via the public `ClumpLayout::new` (L11).

### V2 — the floor is byte-exactly the smallest resolvable budget (no finding)

For `cores` 1–256 × `fastqc ∈ {false, true}`: `resolve_layout(floor)` succeeds and
`resolve_layout(floor − 1 byte)` fails, in **every** cell. The repo test only probes
`floor − 1 MiB`; the byte-exact assertion also passes (see L13).

The earlier `budget/11*10` vs `×11/10` bug really is fixed: multiplying first makes the two
exact inverses. It is exact because `MIN_BIN_BYTES/16 = 65536` divides evenly, so
`dyn_min = 65536 × denom` carries no rounding, and `ceil(11S/10)` is the minimal `F` with
`⌊10F/11⌋ ≥ S`.

**This exactness is load-bearing, not decorative.** `main.rs:184-199` gates on the floor and then
calls `resolve_layout` with a bare `?` (`main.rs:201`). If the floor were even one byte loose,
a budget that passed the gate could bail and turn the intended warn-and-fall-back into a hard
abort. Worth a comment saying so, since a future edit to either function could break it silently.

The **bail message's** suggested minimum is achievable: `div_ceil(1024*1024)` rounds up, so
`printed_MiB × 1 MiB ≥ floor`. Correct direction. `main.rs`'s sibling message is not — that is H1.

---

## High

### H1 — the fallback warning advertises a minimum that does not work

`src/main.rs:189-193`

```rust
"WARNING: --memory {} is too small for --clumpify at --cores {} (need ≥ {} MiB).",
cli.memory, cli.cores, min_required / (1024 * 1024),
```

Plain division **truncates**. The floor is almost never a whole number of MiB, so the printed
value is below the real floor and a user who follows the advice verbatim gets the identical
warning again.

Measured across the sampled cells — `resolves?` is whether passing exactly the advertised MiB
clears the floor:

| cores | fastqc | true floor | `main.rs` says | resolves? | bail says | resolves? |
|---|---|---|---|---|---|---|
| 2 | no | 280.500 MiB | `≥ 280` | **no** | `≥ 281` | yes |
| 2 | yes | 333.300 MiB | `≥ 333` | **no** | `≥ 334` | yes |
| 4 | no | 289.300 MiB | `≥ 289` | **no** | `≥ 290` | yes |
| 4 | yes | 394.900 MiB | `≥ 394` | **no** | `≥ 395` | yes |
| 8 | no | 332.200 MiB | `≥ 332` | **no** | `≥ 333` | yes |
| 8 | yes | 543.400 MiB | `≥ 543` | **no** | `≥ 544` | yes |
| 16 | no | 418.000 MiB | `≥ 418` | yes | `≥ 418` | yes |
| 16 | yes | 840.400 MiB | `≥ 840` | **no** | `≥ 841` | yes |

7 of 8. **Failure scenario:** `trim_galore --clumpify --cores 8 --memory 200M` prints
"need ≥ 332 MiB"; the user reruns with `--memory 332M` and gets the same warning, with no hint
that 333 would have worked.

The truncation predates this change, but this diff edited that exact call and fixed the
*sibling* message in `resolve_layout` to `div_ceil` — leaving the two inconsistent with each
other is new. One-line fix:

```rust
min_required.div_ceil(1024 * 1024),
```

### H2 — the FastQC term is an addend where the calibration said `max`, and that silently disables clumping at the default budget

`src/clump.rs:174-180`, `src/clump.rs:196-199`

`CALIBRATION.md:114-118` is unambiguous:

> `f0` becomes ~205 MiB with FastQC handled as a **`max`, not an addend** (FastQC runs after
> trimming, so the peak is `max(trim, fastqc)`)

The implementation makes it an addend to the *static* reservation, and per-core. On the main
`--clumpify` path the spec is mechanically right: the pool is a local `Vec<PairedBin>` /
`Vec<SingleBin>` inside `parallel.rs::read_pairs_clumpy` (`parallel.rs:813`) and
`::read_single_clumpy` (`parallel.rs:1218`), so it is dropped when those reader functions return
— well before the trim call at `main.rs:1351` returns, and long before `fastqc::run` at
`main.rs:1482`. Charging both concurrently shrinks the pool to make room for memory that is
never co-resident.

The cost is a user-visible regression at the **default** `--memory 1G`:

| | old floor | new floor | default 1 GiB? |
|---|---|---|---|
| `--fastqc --cores 16` | 620 MiB | 840 MiB | still fits |
| `--fastqc --cores 21` | 654 MiB | **1026 MiB** | **falls back to plain mode** |
| `--fastqc --cores 32` | 728 MiB | **1434 MiB** | falls back |
| `--fastqc --cores 64` | 944 MiB | **2622 MiB** | falls back |

**Failure scenario:** `trim_galore --clumpify --fastqc --cores 24 <fastqs>` on a 32-core node
with default memory. Before: clumped, ~24% smaller output. Now: a warning scrolls past and the
run silently produces un-clumped output. The user's only signal is one line of stderr in the
middle of a pipeline log. Nothing is wrong with the data — the *feature they asked for* is off.

The term is also **unbounded in `cores`** while being measured only to 16 (76/139/283 MiB at
4/8/16). Extrapolating 24 MiB/core to 64 cores predicts 1536 MiB of FastQC overhead; FastQC's
per-thread state is bounded by the work available, and each `fastqc::run` call processes a
single file, so linearity past 16 is an assumption, not a measurement.

There is a real counter-argument I want to state fairly: under mimalloc, freeing the pool does
not necessarily return RSS, so peak may be closer to additive than `max` suggests. But that
argues for a **retention** term, not a *per-core FastQC* term — retention scales with the pool,
which is already in the dynamic term. Recommend one of:

1. adopt the spec's `max(trim_static + dyn, fastqc_static)` form on the trim path; or
2. cap the term (`24 MiB × min(cores, 16)`) so an unmeasured extrapolation cannot switch the
   feature off; or
3. keep the addend and record it in the plan's **Deviation** section with the retention
   rationale — right now the plan documents only the `--clump_only` deferral, so this deviation
   from an explicit spec item is invisible to the next reader.

Note it is defensible *as written* on the `--clump_only` path — see M3, where the pool genuinely
is still live.

### H3 — `--fastqc_args "-t N"` decouples FastQC's real thread count from the `cores` the model charges

`src/fastqc.rs:85-88` lets `--fastqc_args` override `config.threads`:

```rust
"-t" | "--threads" => { ... config.threads = v; ... }
```

The reservation is `24 MiB × cli.cores`, and `fastqc_requested()` is `self.fastqc ||
self.fastqc_args.is_some()` — so `--fastqc_args` alone turns the term on while independently
setting the thread count the term is supposed to proxy.

**Failure scenario (under-reserve, breaches the promise):**
`trim_galore --clumpify --cores 4 --memory 2G --fastqc_args "-t 32"`. The model reserves
`224 + 4×24 = 320 MiB` and sizes the pool against the remaining ~1.5 GiB. FastQC then runs 32
threads, which by the model's *own* coefficient costs ~768 MiB. Under-reserved by ~670 MiB —
against the exact guarantee the CHANGELOG now advertises.

**Failure scenario (over-reserve):** `--cores 16 --fastqc_args "-t 1"` charges 384 MiB for a
single-threaded FastQC and lifts the floor to 840 MiB for nothing.

Cheapest fix: derive the FastQC thread count the same way `fastqc.rs` does (parse
`--fastqc_args` for `-t`/`--threads`, else `cli.cores`) and charge *that*. If the term is capped
per H2, the under-reserve shrinks but does not vanish.

---

## Medium

### M1 — `resolve_layout_predicted_peak_fits_budget` is vacuous, confirmed by mutation

`src/clump.rs:534-558`

I proved analytically in V1 that the assertion holds for any constants with margin ≤ 1, then
confirmed it empirically. I re-implemented the eight layout tests against a parameterised model
and perturbed each constant:

| mutation | default | scales | pins_cell | fastqc_shrinks | **predicted_fits** | bails_small | bails_static | floor_smallest |
|---|---|---|---|---|---|---|---|---|
| σ 23→20 | FAIL | FAIL | FAIL | FAIL | **pass** | pass | pass | pass |
| σ 23→26 | FAIL | FAIL | FAIL | FAIL | **pass** | pass | pass | pass |
| k 64→56 | FAIL | FAIL | FAIL | FAIL | **pass** | pass | pass | pass |
| k 64→72 | FAIL | FAIL | FAIL | FAIL | **pass** | pass | pass | pass |
| k 64→28 (revert to old 7/4) | FAIL | FAIL | FAIL | FAIL | **pass** | pass | pass | pass |
| static 224→205 | pass | pass | FAIL | pass | **pass** | pass | pass | pass |
| static 224→256 | FAIL | pass | FAIL | pass | **pass** | pass | pass | pass |
| static 224→512 (revert) | FAIL | FAIL | FAIL | pass | **pass** | pass | pass | pass |
| fastqc 24→0 (term removed) | pass | pass | pass | FAIL | **pass** | pass | pass | pass |
| fastqc 24→16 | pass | pass | pass | FAIL | **pass** | pass | pass | pass |
| fastqc 24→32 | pass | pass | pass | FAIL | **pass** | pass | pass | pass |
| **fastqc per-core → flat 96 MiB** | pass | pass | pass | **pass** | **pass** | pass | pass | pass |
| fastqc per-core → flat 288 MiB | pass | pass | pass | FAIL | **pass** | pass | pass | pass |
| margin 10/11 → 1/1 (removed) | FAIL | FAIL | FAIL | pass | **pass** | pass | pass | pass |
| margin 10/11 → 9/10 | FAIL | pass | FAIL | pass | **pass** | pass | pass | pass |

**It survives all 15.** So: still vacuous as a calibration test — it cannot fail for any value
of any constant this change introduced.

It is not *entirely* worthless: it catches an inverted margin (`11/10`) and a predictor that
double-counts `static`. But its guard on the inverse is only good to within the margin — a
predictor using `k = 4.5` while the sizer uses `k = 4.0` **passes**, because the 5% denominator
mismatch hides inside the 9% margin. So it cannot detect banner/sizer drift below ~10%, which
is the drift most likely to appear in a future edit.

Recommend replacing `predicted ≤ budget` with the exact relationship, which would catch all
three:

```rust
// predicted peak is the budgeted total, less only the two floor roundings
let budgeted = budget * 10 / 11;
assert!(budgeted - predicted_peak <= layout.n_bins as u64 * 32); // slack ≤ rounding
```

`resolve_layout_bails_when_too_small`, `resolve_layout_bails_below_static_overhead` and
`clumpify_floor_matches_the_smallest_resolvable_budget` also survive all 15 mutations. That is
acceptable for the two bail tests (they guard saturation paths, not constants), and the floor
test guards a *relationship* — the one that caught the real `×11/10` bug — rather than a value.
Worth knowing that `resolve_layout_pins_the_default_cell` is doing nearly all the constant-pinning
work: it is the only test that fails under every non-FastQC perturbation.

### M2 — the per-core-ness of the FastQC term has no test

`src/clump.rs:512-532`

Row 12 above is the finding: **`fastqc per-core → flat 96 MiB` passes every test.** The test
probes a single core count (4), where `24 × 4 = 96`, so a flat 96 MiB addend is
indistinguishable from the per-core form. The design decision that drives H2's regression is
therefore unverified.

The test also hardcodes `24u64 * 4 * 1024 * 1024 * 16 / 624`, restating the constant and the
denominator rather than deriving them, so it pins magnitude at one point and nothing else.

Fix is small — assert the delta at two core counts and that it scales:

```rust
let d4 = plain4.bin_byte_budget - qc4.bin_byte_budget;
let d8 = plain8.bin_byte_budget - qc8.bin_byte_budget;
// per-core: the static delta doubles from 4 to 8 cores
```

### M3 — `--clump_only` charges 64 units per worker for a single-threaded mode, and the FastQC addend now pushes a working invocation into a hard abort

`src/clump_only.rs:275, 413, 748, 905`

There is no `spawn` or `thread::` anywhere in `clump_only.rs` — the mode is single-threaded, as
its own module doc says (`clump_only.rs:33`). It nonetheless passes `cores.max(1)` into
`resolve_layout`, so the denominator charges `64 × cores` for workers that do not exist. At
`--cores 16` that is 1024 of a 2496 denominator, sizing bins ~5.8× smaller than the one real
worker needs.

That over-charge predates this change, but the change **amplifies it twice**: the per-core
coefficient went from `7/4 = 1.75` to `64/16 = 4.0`, and the FastQC addend now adds another
`24 MiB × cores` on the same conflated parameter. Result:

**Failure scenario:** `trim_galore --clump_only --fastqc --cores 16 --memory 800M in.fq.gz`.
Old floor 620.0 MiB → ran fine. New floor 840.4 MiB → `resolve_layout` bails and
`clump_only.rs:275` propagates it with a bare `?`. Hard abort, no output, on a command that
worked before.

The plan's deferral is **half defensible**. Deferring the *fallback* is right — falling back to
plain mode is meaningless for a mode whose entire purpose is reordering, so an error is the
correct behaviour and the only question is its wording. But the **mis-sizing** is a separate,
fixable bug that the deferral note folds in with it, and it is what moved the floor.

The obvious fix — pass `1` — is **wrong**, and worth recording so nobody tries it:
`layout.n_bins` feeds `clump::bin_for` (`clump_only.rs:312, 461, 794`), so lowering `cores`
lowers `n_bins`, which changes bin assignment and therefore output order. Byte-identity is this
mode's load-bearing invariant. The change needs three separable inputs, not one `cores`:

- `n_bins` — must stay `max(16, 4 × cli.cores)` to preserve output order;
- pool/worker term — `1` for `clump_only`;
- FastQC threads — the user's value (H3).

With the worker term at 1 and `n_bins` pinned at 64, the `c16 --fastqc` floor is 774.4 MiB and
the `--memory 800M` invocation above keeps working.

Note the FastQC **addend** is correct on this path, unlike the trim path in H2: `bins` is
declared at `clump_only.rs:293` and never explicitly dropped, so the whole
`n_bins × bin_byte_budget` of pre-sized capacity is still live when `fastqc::run` is called at
`clump_only.rs:375`. If H2 adopts a `max` form, adding `drop(bins);` before the `fastqc::run`
calls would make `max` sound here too — and would itself cut this path's real peak.

### M4 — CHANGELOG scope, and it omits both user-visible costs

`CHANGELOG.md:8-11`

> `--clumpify` now keeps peak memory inside the `--memory` budget it prints, **at every core
> count and budget measured**

Checked against the data. The post-change cells are the six in `PLAN.md:380-387`: `1G/c4`,
`2G/c4`, `1G/c8`, `1G/c4 --fastqc`, `1G/c8 --fastqc`, `2G/c16 --fastqc` — ratios 0.862, 0.887,
0.828, 0.739, 0.617, 0.682. All inside budget, so the sentence is **literally true**. But the
scope is: three core counts, two budgets, **macOS/arm64 only**. `CALIBRATION.md`'s Linux data is
all *pre*-change, and its own Linux `c16` cell missed by 22% (`CALIBRATION.md:79`). A changelog
reader will not decode "measured" as "six cells on one laptop".

The previous entry was rated High for over-claiming; this is milder because the qualifier is
present and honest, hence Medium. But two facts a reader "deciding whether an upgrade affects
them" needs are missing:

- **bin sizes shrink at large budgets**: `4G/c8` goes 66.37 → 44.87 MiB (−32.4%), `2G/c4`
  56.89 → 42.00 MiB (−26.2%), `4G/c2` 152.51 → 112.89 MiB (−26.0%). That is the price of the
  fix and it affects compression.
- **`--fastqc` now needs materially more memory** (H2), to the point of switching clumping off
  at the default budget above 20 cores.

House style is one sentence per entry; this is one sentence carrying three claims across four
lines, which is borderline. Given the standing feedback on changelog verbosity, I would cut it
to the claim plus the cost:

```
- **`--clumpify` peak memory now stays inside the printed `--memory` budget in every cell
  measured** ([#439](…)) — bins shrink 25–30% at budgets ≥ 2 GiB, and `--fastqc` runs
  reserve more; see `docs/performance/clumpy.md`.
```

### M5 — the plan's deviation note compares floors from different core counts

`PLAN.md:417-419`

> the FastQC floor rises with cores (~841 MiB at `--cores 16 --fastqc` against 536 MiB before)

536 MiB is the **`cores 2`** floor under the old model. At `cores 16` the old floor was
**620.0 MiB** (`512 MiB + ceil(1 MiB × 432 / 4)`). The real c16 delta is 620 → 841 MiB (+36%),
not 536 → 841 (+57%). The note reads as a like-for-like comparison and is not one. Worth fixing
because this paragraph is the record of a decision being handed to Felix.

---

## Low

**L1 — stale comment naming a deleted constant.** `src/clump_only.rs:1188`:
"`resolve_layout` requires STATIC_OVERHEAD (512 MiB)". `STATIC_OVERHEAD_BYTES` no longer exists;
`grep -rn STATIC_OVERHEAD src/` returns this line only. The `1 GiB` value it justifies is still
comfortably above the floor, so only the comment is wrong.

**L2 — `clumpify_min_memory_bytes` uses unchecked `*` where `resolve_layout` saturates.**
`src/clump.rs:188`: `(static_bytes_for(..) + dyn_min) * MARGIN_DEN`. `--cores` has no clap upper
bound (`cli.rs:351`), and the product exceeds `u64` for `cores > ~4.7e10`: panics in a debug
build, wraps in release. It degrades to an error rather than bad output either way, so this is
consistency, not risk — but the asymmetry with `resolve_layout`'s all-saturating arithmetic is
the kind of thing that reads as an oversight later.

**L3 — `ClumpLayout::new` is public and unvalidated.** `src/clump.rs:150`. `predicted_peak_bytes`
uses plain `*`, so a caller-supplied `(bin_byte_budget, cores)` pair can overflow and wrap in
release. Values from `resolve_layout` provably cannot (V1). Only test code calls `new` today; if
it stays public, `debug_assert` the product or make it `pub(crate)`.

**L4 — `parallel.rs:1615` hardcodes `224 * 1024 * 1024`.** `STATIC_TRIM_BYTES` is private to
`clump.rs`, so the test helper duplicates it and will drift silently when the constant is
retuned. `pub(crate)` on the constant, or a `ClumpLayout::for_test(n_bins, budget, cores)` that
fills in the real value.

**L5 — `clumpify_floor_matches_the_smallest_resolvable_budget` under-claims its own name.**
`src/clump.rs:579-595`. The name says "smallest resolvable budget"; the body probes
`floor − 1 MiB`. I verified the byte-exact form (`floor − 1`) also fails for `cores` 1–256 in
both fastqc settings, so tightening it costs nothing and earns the name. Also note
`floor - 1024*1024` would underflow if a floor ever dropped below 1 MiB — not reachable today
(minimum floor is 280.5 MiB) but it is an unguarded subtraction in a test that sweeps `cores`.

**L6 — the banner's 9% margin is invisible.** `main.rs:202-208` prints
`predicted peak ≈ 930 MiB` for `--memory 1G` — and the ratio is exactly 0.909 at every budget
and core count I sampled (1G/c4, 2G/c4, 4G/c8, 16G/c16). Not a contradiction; the prediction is
honest. But a user reading "930 of 1024" may conclude 94 MiB is slack they can reclaim by
lowering `--memory`, when it is deliberate headroom sized by the 96 MiB macOS spread. One
word — "predicted peak ≈ 930 MiB (9% margin held back)" — would close it.

**L7 — constants' doc comments carry per-platform measurements, against house style.**
`src/clump.rs:191-204, 221-226`. `CLAUDE.md` asks for one line stating the current fact, with
"no measurements, no per-platform numbers", and its flagged example is precisely a per-platform
figure list ("370 under docker, 374 under conda"). `STATIC_TRIM_BYTES`'s comment is three lines
ending "205 MiB on macOS/arm64 and 173 MiB on Linux/x86_64"; the `resolve_layout` docstring runs
~20 lines including which cells were held out and what they predicted.

My judgment: **trim the numbers, keep the pointer.** The pointer to `CALIBRATION.md` is already
there (`clump.rs:226`) and is the durable mechanism — the per-platform figures are a snapshot of
one fixture at one read length, which `CALIBRATION.md:129-139` itself flags as the model's main
limitation, so they will rot faster than anything else in the file. Two lines would carry the
load: what σ and k are, that they are the worse of two platforms, and where the fit lives.

Rated Low, not Medium, for two reasons: the constant this replaced had the same style, so it is
not a regression; and the letter of the rule is arguably in tension with its own carve-out
("put an invariant on the thing it constrains") for a calibrated constant, where provenance *is*
the invariant. Felix's call.

---

## Phase 3 scope — what `docs/performance/clumpy.md` now states incorrectly

Out of scope for this diff, listed so nothing is missed. Every line below is currently wrong.

| line | published | actual now |
|---|---|---|
| 96 | `usable = memory − 512 MiB` | `memory × 10/11 − 224 MiB` (plus `24 MiB × cores` if FastQC) |
| 96 | "(FastQC + allocator + runtime overhead)" | FastQC is now a separate conditional per-core term, not part of the fixed reservation |
| 99 | `B = 4 × usable / (5·n_bins + 7·cores)` | `B = 16 × usable / (23·n_bins + 64·cores)` |
| 102 | `--cores 8 --memory 1G` → "32 bins × 12 MB" | 32 bins × **9 MiB** |
| 102 | `--cores 8 --memory 4G` → "32 bins × 66 MB" | 32 bins × **45 MiB** |
| 108 | "minimum of ~535 MiB" | **280.5 MiB** at `--cores 2` |
| 108 | "mostly from a fixed 512 MiB reservation for FastQC, allocator, and runtime" | 224 MiB, and FastQC is conditional |
| 108 | "stays in the 535–730 MiB range for any sensible core count" | **280.5 – 1434.4 MiB** (c2 no-FastQC → c32 with FastQC) — wrong at both ends, and the upper end is now ~2× the published ceiling |
| 113 | sample warning "at --cores 6 (need ≥ 552 MiB)" | 310 MiB (311 once H1 is fixed) |

Reference floors for the rewrite, MiB:

| cores | old | new, no FastQC | new, FastQC |
|---|---|---|---|
| 2 | 535.5 | 280.5 | 333.3 |
| 4 | 539.0 | 289.3 | 394.9 |
| 6 | 552.5 | 310.8 | 469.2 |
| 8 | 566.0 | 332.2 | 543.4 |
| 16 | 620.0 | 418.0 | 840.4 |
| 32 | 728.0 | 589.6 | 1434.4 |

Line 108's claim is the one to prioritise: it is the only published number that a user sizing a
cluster job would rely on, and the FastQC column now exceeds the default `--memory 1G`.

Also stale, one line: `docs/src/content/docs/reference/changelog.md:199` describes the 512 MiB
reservation. That is a synced historical entry, so per `CLAUDE.md` leave it alone.

---

## Recommendations, in order

1. **H1** — `min_required.div_ceil(1024 * 1024)` in `main.rs:192`. One line, no design question.
2. **H2** — decide the FastQC form: `max`, capped addend, or documented deviation. This is the
   only finding with a behaviour question in it, and it gates the c≥21 default-budget fallback.
3. **H3** — charge FastQC threads, not `cli.cores`, or state the limitation.
4. **M3** — split `cores` into `n_bins` / worker-count / FastQC-threads so `--clump_only` stops
   paying for 15 workers it does not have; keep `n_bins` pinned to preserve byte-identity.
5. **M1, M2** — make the two weak tests bite: exact inverse, and per-core scaling at two counts.
6. **M4, M5** — trim the CHANGELOG to the claim plus the cost; fix the c16-vs-c2 floor comparison
   in the plan.
7. **L1–L7** — housekeeping; L1 is a one-word edit.

Nothing here argues against the change. The refit is the right move and the two properties it
claims are genuinely true. The FastQC term is the part that was never measured the way the rest
was, and it shows.
