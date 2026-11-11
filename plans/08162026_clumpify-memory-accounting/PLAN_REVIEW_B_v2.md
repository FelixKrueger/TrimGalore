# Plan Review B — v2 of the `--clumpify` memory-accounting plan (#439)

**Reviewer:** B (independent; did not read `PLAN_REVIEW_A_v2.md`)
**Target:** `plans/08162026_clumpify-memory-accounting/PLAN.md` (v2)
**Repo state:** `dev` @ `2ab0918`, tree clean. All line numbers re-derived at this tip.
**Data re-analysed:** `scratchpad/rss/results_sweep.tsv` (19 cells), `results_investigate.tsv`,
`results_div_diverse_darwin.tsv`, `oxy_results.tsv`, `results_fixedbin.tsv`, `sweep.py`.

## Verdict

**Do not implement Phase 0 as scoped.** v2 is a large improvement on v1 — it stops asserting
coefficients it cannot support, and staging the measurement first is right. But three things are
wrong at the load-bearing level:

1. The plan's **central diagnostic claim is a fitting artefact.** "The intercept is core-dependent,
   so there is a per-core cost the model has no term for" does not survive a joint fit. A
   3-parameter model with **no** per-core fixed term fits all 19 macOS cells *including* the c8
   cell the plan says contradicts the extrapolation. The per-series intercepts of 151/302 MiB are
   slope↔intercept trade-off inside each series, not evidence of a missing term.
2. The **Goal is unsatisfiable as written.** "Prediction ≥ actual in every cell" and "bin size at
   `1G --cores 4` unchanged" cannot both hold: that cell measures 1.025–1.081× today with B fixed
   at 18.96 MiB. Validation rows 1 and 3 contradict each other. Something has to give, and the
   plan does not say which.
3. **The mechanism is already known and the plan un-knows it.** `FINDINGS.md:61-78` derives the
   per-core term as *(queue depth + 1) × spine factor* and confirms it arithmetically on Linux
   (695 predicted vs 672 measured). v2 demotes that to Step 8 "Phase 0 decides whether it is
   needed at all" and replaces a derivable coefficient with an ~80-run curve-fit. The 80 runs are
   mostly redundant, and the one cell that would actually resolve the open question (`cores 3`) is
   not in the matrix.

The fix I would ship instead: **cut the actual memory (queue depth + the 300 B/pair output
reserve), derive σ and k from the mechanism, drop the fixed reservation from 512 → ~250 MiB, and
make FastQC a `max` rather than an addend.** That is bounded work with ~15 confirming runs, and it
lets bin size *grow* at the default instead of shrinking.

## What I probed and found correct

Every arithmetic claim in "What the measurements establish" reproduces:

| Plan claim | My recomputation | Verdict |
|---|---|---|
| `c2: 151.1 + 32.71·B` | 151.14 ± 29.55 + 32.709 ± 0.809·B, R²=0.9970 | ✅ |
| `c4: 301.6 + 38.42·B` | 301.63 ± 13.30 + 38.416 ± 0.364·B, R²=0.9996 | ✅ |
| under-charge 1.39× / 1.42× | 32.709/23.5 = 1.3919; 38.416/27.0 = 1.4228 | ✅ |
| `301.6 + 38.42·B ≤ 1024 → B ≤ 18.8` vs shipped 18.96 | 18.804 vs 4·512/108 = 18.963 | ✅ |
| `B = 28 MiB` at 1 GiB ≈ 1.37× | interpolating 27.41→1378.2 and 37.89→1752.2 gives 1399 MiB → 1.366 | ✅ |
| breach grows with budget: 1.03× @1G, 1.21× @2G (c4) | 1.025; 1.215/1.211 | ✅ |
| `--compression` does not move RSS | z1 1033.1/1079.7 vs z6 1050.2 at the same cell | ✅ |
| bin size pays compression: 44.5 → 20.2 MB | out_bytes 44,515,492 → 20,159,662 | ✅ |
| ~75 MiB/core extrapolation contradicts measured 0.94× at 1G/c8 | 1009.5 MB/2²⁰ = 962.7 MiB ÷ 1023 = 0.941 | ✅ |
| v1-era harness inflated ratios 4.86% | `bench_div.py` divides by 1e6, banner is MiB | ✅ |

Code-level claims, all verified at this tip:

- `ClumpLayout` is built with a struct literal in test code — `src/parallel.rs:1611-1616`. A
  private field would be E0451 across the module boundary. **Step 2 is right.**
- `Cli::fastqc_requested()` exists at `src/cli.rs:592-594` and is `self.fastqc ||
  self.fastqc_args.is_some()`. **Step 4 is right.**
- Four `resolve_layout` call sites in `src/clump_only.rs` — lines **275, 413, 748, 905** — and each
  enclosing function does carry `fastqc_requested: bool` (params at 258, 396, 743, 893).
  **Step 6 is right.**
- `resolve_layout_predicted_peak_fits_budget` (`src/clump.rs:465-480`) **is** a tautology:
  `predicted = STATIC + floor(4·usable/denom)·denom/4 ≤ STATIC + usable = budget` for any
  coefficients, because the same `denom` appears in both directions and the floor only helps.
  **Right, and worth deleting rather than keeping.**
- Docs `--cores 8 --memory 1G` → "32 bins × 12 MB" is wrong: `4·(1024−512)/216 = 9.481 MiB`, and
  the banner would print `9`. (`docs/.../clumpy.md:102`.) The neighbouring 4G value, 66 MB, is
  correct. **Right.**
- Banner prints MiB labelled MB — `src/main.rs:203`. Floor warning at `:189` correctly says MiB.
- The `--fastqc --memory 570M --cores 8` demotion example is well chosen: today's floor at c8 is
  `512 + 1 MiB·216/4 = 566 MiB`, so 570M resolves today by 4 MiB. **Right.**
- `--clumpify` requires `--cores >= 2` (`src/cli.rs:734`), so there is no c1 cell to worry about on
  the pipeline path; `--clump_only` is explicitly exempt and does reach `resolve_layout(mem, 1)`.

Also correct and worth keeping as-is: the `BTreeMap` needs no cap; `--compression` is irrelevant to
RSS; Phase 0 having no code changes.

---

# Critical

## C1. The "missing per-core term" is a fitting artefact; a 3-parameter model already fits everything

The plan fits each core count separately, then reads the two intercepts (151, 302) as a per-core
fixed cost of ~75 MiB/core, and concludes "no coefficient can substitute for a missing term"
(PLAN.md:50-51). That inference does not survive pooling the data.

Fitting **one** shared structure to all 19 macOS cells (7 × c2, 11 × c4 incl. the z1 reps, 1 × c8):

```
RSS_MiB = f0 + B·(σ·n_bins + k·cores)
        = 244.2 + B·(1.328·n_bins + 4.614·cores)      rmse 44.8 MiB (≈4%)
```

Residual on the c8 cell the plan calls a contradiction: **−34.3 MiB (3.5%, conservative side).**
No per-core intercept is needed. Adding one (4 params) changes rmse from 44.8 → 44.1 — i.e. buys
nothing, and pushes the c8 residual the *wrong* way (−65.5).

Why the per-series fits mislead: within one series B spans only 12.8–56.9 MiB while rep-to-rep
noise is up to 4.9%, so intercept and slope trade off almost freely. Two different contrasts on
the same data give wildly different per-core coefficients:

| estimator | per-core B-coefficient |
|---|---|
| slope difference `(s4−s2)/2` (what `sweep.py` prints, what the plan uses) | **2.85** |
| fixed-B pair `(1050.2 − 813.2)/2 ÷ 19 MiB` | **6.24** |
| joint LS over all 19 cells | **4.61** |

A 2.2× spread in the parameter the plan wants to calibrate. That is the real identification
problem, and it is *not* "three parameters from two core counts" (three are identified). It is that
`cores ∈ {2,4}` at `n_bins = 16` yields exactly **one** contrast to split between a per-core
*intercept* and a per-core *slope*. Four unknowns, four fitted numbers, zero degrees of freedom, no
test of form.

**Consequence for the model shipped by Phase 1:** the pooled model says the shipped decomposition is
wrong in *two* directions that cancel at the default and diverge with budget — the fixed term is
roughly **double** what it should be (512 vs ~244 MiB) and the B-coefficients are ~18% too small.
That is exactly why the breach grows with budget (PLAN.md:55): at 1G, `usable = 512` and the two
errors cancel; at 2G, `usable = 1536` and the slope error dominates.

Honest bin size under the pooled mean fit, `--memory 1G`:

| cores | pooled-fit honest B | shipped B | ratio |
|---|---|---|---|
| 2 | 25.6 MiB | 21.79 | 1.17× |
| 4 | 19.6 MiB | 18.96 | **1.04×** |
| 8 | 9.8 MiB | 9.48 | 1.04× |
| 16 | 4.9 MiB | 4.74 | 1.04× |

So the mean-fit model *raises* B slightly everywhere at 1G and only bites at large budgets — which
is a far better story than the plan's framing ("today's B is right, the decomposition is wrong").
The plan gets to the same place for the wrong reason and then over-constrains itself.

**Action:** replace measurement 2 with the pooled fit (or the mechanism derivation in C3), and drop
the premise that Phase 0 must discover a term structure. What Phase 0 must do is *test* the form,
which needs one extra core count inside the clamp region, not a full factorial.

## C2. Least squares cannot produce an upper bound — the promise needs an envelope plus a margin

The Goal is one-sided: printed prediction ≥ actual peak. Every fit in the plan (and in `sweep.py`)
is least squares **through the mean**. By construction ~half the runs land above the line, so
Validation row 1's "ratio ≤ 1.0 in **every** cell" will fail on the very data used to fit it. This
is not a nuance; it is the difference between a promise and an average.

From the pooled fit, per-cell ratios of measured ÷ predicted run 0.865 … **1.081**, with the worst
under-predictions clustered at c4 near the default (max positive residual +83.2 MiB). To cover
every measured macOS cell the B-term must be inflated **1.108×** (σ 1.328 → 1.471, k 4.614 →
5.110). And the observed rep-to-rep spread — 4.9% at c2/B56.9, 4.5% at c4/z1_1G, 2.3% between two
independent measurements of the *same* cell on two binaries (1026.5 MiB in
`results_div_diverse_darwin.tsv` vs 1050.2 MiB in `results_sweep.tsv`, same fixture, same cell,
both mimalloc) — means an envelope over 1–2 reps will itself be exceeded. A ~10% margin on top is
the minimum defensible choice.

The price of the promise, at `--memory 1G --cores 4`:

| model | B | vs shipped | est. output size on the sweep fixture |
|---|---|---|---|
| shipped | 18.96 | — | baseline |
| pooled mean fit | 19.6 | +3% | −0.5% |
| + envelope (cover all cells) | 17.7 | −7% | +1.0% (z1) … +2.6% (z6) |
| + envelope + 10% margin | 15.6 | −18% | ~+2.7% (z1) … +8% (z6) |

(sensitivity `dOut/dB` measured locally from `out_bytes`: −0.55 MB/MiB at z1, −0.80 MB/MiB at z6
near B = 19.)

**Action:** state explicitly which quantile the model targets, fit the **upper envelope** (max over
reps) not the mean, and name the safety margin as a number in the plan. Then decide the
compression trade with eyes open — see C4, which can pay for the margin and more.

## C3. Goal and Validation rows 1 & 3 are mutually unsatisfiable; the structural fix is not optional

At `1G --cores 4` with B = 18.96 MiB the measured ratios are 1.025 (z6), 1.008 and 1.053 (z1),
1.003 (the other harness). The prediction there *is* the budget (1025 MiB). So:

- Validation row 1 demands ratio ≤ 1.0 in every cell.
- Validation row 3 demands B stays "within a few percent of 18.96 MiB".
- With actual memory unchanged, no arithmetic in `resolve_layout` can satisfy both. Shrinking B by
  a few percent buys ~40 MiB of headroom against a 4.9% noise band ≈ 50 MiB.

The plan half-acknowledges this (Remaining risk (c)) but never resolves it, and its Goal sentence
("measurement says today's `B` is already within 1% of the honest maximum") asserts a 0.85% margin
inside a 2.3% same-cell spread. That precision is not available in this data.

Only three ways out, and the plan must pick one:

- **(a) Reduce actual memory** so the promise holds at today's B — C4. Preferred.
- **(b) Accept a bin-size reduction** at the default (7–18%, C2) and take the compression hit.
- **(c) Re-word the promise** — make the banner say "target" not "predicted peak", and state a
  documented tolerance. Cheapest, but it is a docs change, not a fix for #439.

Read literally, the current Goal picks (a) but relegates the (a)-work to Step 9 as a "candidate".
That is the wrong ranking.

## C4. The primary fix should be the queue depth and the 300 B/pair output reserve — not buffer recycling

The per-core cost is not a mystery to be fitted. `FINDINGS.md:61-78` already derives it, and the
code confirms it:

- `src/parallel.rs:119` — `mpsc::sync_channel::<PairedWork>(2)` per worker. A worker can hold **2
  queued bins + 1 in hand = 3 full bins** of text + spine outside the bin pool.
- `src/parallel.rs:895-907` — `PairedBin::push` `reserve_exact`s the whole spine up front, and
  `take()` (`:913-919`) hands those *fully touched* Vecs to the worker. So the in-flight bytes are
  resident, not merely reserved.
- Hence `k ≈ (queue_depth + 1) × spine_factor = 3 × 1.39 = 4.17` at 65 bp, against a joint-fit
  4.61 and a charged **1.75**. And `σ ≈ spine_factor = 1.39` against a joint-fit 1.328 and a
  charged 1.25. Both fitted coefficients are within noise of the mechanism.

That is a *derivation*, not a calibration: the model's per-core term is wrong because it charges one
in-flight bin per worker while the channel permits three. Two consequences the plan misses:

1. **Reducing the depth reduces the coefficient.** `sync_channel(1)` (or a clumpy-only depth of 1)
   takes `k` from ~4.2 to ~2.8; depth 0 (rendezvous) to ~1.4. Carrying that through the envelope
   (×1.108) and a 10% margin on the whole prediction, honest B at 1G/c4 becomes
   `(1024/1.1 − 244)/(1.471·16 + 3.08·4) = 19.2 MiB` at depth 1 and **23.1 MiB** at depth 0 — i.e.
   at or above today's 18.96 with the margin already paid, versus 15.6 MiB if the depth stays at 2.
   This is the single change that dissolves the Goal/Validation conflict in C3.
2. **The output reserve is the other half.** `src/parallel.rs:349-351` — `cap = reads_r1.len() *
   300` for *each* of `buf_r1` and `buf_r2`, i.e. 600 B reserved per pair while the pair's own text
   is 378 B at 65 bp (294 B at 51 bp) and the gzip output is ~2% of that. In bin units that is a
   **1.59·B** reservation per in-flight batch at 65 bp, and **2.04·B** at 51 bp, against a charged
   0.5·B. Whether the untouched tail becomes resident is allocator-dependent — that is precisely
   the cheap experiment to run (below) — but the plan calls this a "cheap adjacent fix" (Step 10)
   when it is a candidate for the largest single term in the per-core coefficient.

By contrast, **Step 9 (buffer recycling) is the weakest of the three structural candidates**, and
the plan nominates it as the one that "may subsume most of the recalibration". Against the code:

- The churn it removes is small: at B = 19 MiB over a 2.57 GB fixture there are ~135 flushes × 3
  Vecs ≈ 400 large cross-thread free/alloc pairs in a 10-second run. That is not a retention
  stress test, and mimalloc already recovered the macOS overshoot (1.48 → 1.05) without it.
- It does not touch any of the three actual peak components: resident bins (`n_bins·B`, by design),
  in-flight copies (bounded by channel depth — C4.1), or the output reserve (C4.2). The buffers it
  recycles are the very ones already counted.
- To bound anything, the reader must **block** on a returned buffer, which converts the pipeline to
  bin-granularity backpressure — a throughput risk of exactly the kind the plan (rightly) worries
  about for Step 8. With a fallback pool it bounds nothing.

**Action:** promote Steps 8 and 10 to Phase 1 primary work; keep Step 9 as an optional allocator
experiment. Reword Step 8: the question is not "is a cap needed", it is "the model must charge
`(depth+1)·spine` per core, or the depth must come down and be charged honestly". Gating on clumpy
mode is the right call and is cheap (`clump_layout.is_some()` is already threaded to the reader at
`src/parallel.rs:188`).

## C5. FastQC is sequential, so it belongs in a `max`, not in the fixed term

`src/main.rs:1798` runs `fastqc::run` **after** the trimming pipeline has finished and its outputs
are closed (same shape at `:1481`, `:2269`, `:2427`, `:2546`; `src/clump_only.rs:374, 554, 845,
1095`). The FastQC working set therefore does not coexist with the bin pool. The honest model is

```
peak ≈ max( trim_peak(B, cores), fastqc_peak + retained_after_trim )
```

Step 4 keeps it additive ("it belongs in the fixed term"), which is wrong in both directions:

- **Over-conservative** whenever `trim_peak` dominates — it subtracts ~405 MiB from the bin budget
  for memory that is spent after the bins are gone. At 1G that is the difference between a workable
  and an unworkable FastQC configuration, and it is what produces the demotion the plan flags in
  Step 5.
- **Insufficient** if `fastqc_peak` alone exceeds the budget. The one measurement available
  (`results_investigate.tsv`, `pe_1g_c4_fastqc`, **stock allocator**) is 1936.9 MB — i.e. the
  process high-water mark with `--fastqc` at a 1 GiB budget was ~1847 MiB, *above* the 1512 MB
  no-FastQC peak. Whatever the split between FastQC's own working set and retained trim memory, no
  choice of B fixes it, because neither addend scales with B.

Also note the ~425 MB FastQC figure the plan quotes is a **single cell on the pre-mimalloc
binary**. mimalloc moved the no-FastQC baseline 1512 → 1061 MB at the same cell; the FastQC delta
under the shipped allocator is unmeasured.

**Action:** re-specify Step 4 as a `max`, and add a FastQC-phase cell (mimalloc, with and without
`--clumpify`) to the measurement work. Validation row 5 ("ratio ≤ 1.0" on the FastQC path) may be
infeasible at 1G for reasons unrelated to bin sizing — decide in advance whether the promise covers
the FastQC phase at all, or whether the banner's number is a trim-phase promise with FastQC's cost
stated separately.

---

# Important

## I1. Phase 0's matrix is a full factorial that omits the two factors that matter

`cores ∈ {2,4,8,16} × B ∈ {12.8, 19, 27.4, 37.9, 56.9}` is a superset of what identification needs,
so it is not *insufficient* for σ/k — but the redundancy is in the wrong place and two declared
factors are missing.

What each region can and cannot distinguish, given `n_bins = max(16, 4·cores)`:

| cells | slope contributed | identifying power |
|---|---|---|
| `cores ∈ {2,3,4}` (n_bins clamped at 16) | `16σ + c·k` | **the only region where σ and k separate.** Three core counts here give 6 fitted numbers for 4 parameters → the form becomes testable |
| `cores ≥ 4` | `4c·σ + c·k = c·(4σ + k)` | σ and k are **perfectly collinear** — c8 and c16 both measure the same combination `(4σ+k)`. Useful as extrapolation checks, worth ~2 cells each, not 10 |

So: `cores = 3` is the **single cheapest cell in the whole design** and it is not in the matrix. Add
a B-sweep at c3 (5 runs) and the "is there a per-core fixed term?" question — the one the plan says
Phase 0 exists to answer (PLAN.md:80-83) — is answered with 2 spare degrees of freedom. Trim c8 to
3 bin sizes and c16 to 2, and the matrix shrinks from ~40 to ~20 core runs while answering *more*.

Missing factors, both of which the plan itself declares necessary elsewhere:

- **Read length.** Assumption 4 (PLAN.md:180-181) says the calibration "must cover more than one
  read length or the model will be tuned to one fixture" — and then the cell list has exactly one
  fixture. It is 65 bp (`data/div_R1.fastq.gz`, verified), spine factor 1.39. TrimGalore's core
  audience runs 51 bp RRBS (spine ≈ 1.50 by the same arithmetic) and 150 bp WGBS (≈ 1.21). Since
  both σ and k are ≈ the spine factor (C4), a model fitted at 65 bp under-charges 51 bp input by
  ~8% and over-charges 150 bp by ~13%. Concretely: honest B at 1G/c4 moves from 19.6 MiB (65 bp) to
  18.2 MiB (51 bp) on the mean fit, and to 16.4 MiB once the envelope is applied — *another*
  independent reason the Goal's "hold B at 18.96" cannot hold. Two extra fixtures × 2 cells each =
  4 runs; make it a factor, not a footnote.
- **`--clump_only`.** All four of its `resolve_layout` calls (`src/clump_only.rs:275, 413, 748, 905`)
  use the same model, and Step 7 asks for an explicit `cores` decision — but no Phase 0 cell and no
  Validation row measures that path at all. Its peak shape is structurally different: bins resident
  + **one** serialised member built in a `Vec` (`src/clump_only.rs:213-227`), no worker batches, no
  queues. Today it is over-charged (safe). Note the direction of Step 7's proposed fix: passing
  `cores = 1` *raises* B by 24% (denominator 87 vs 108) while removing the per-core reservation, so
  it makes the actual peak **higher** on the path nobody has measured. Measure before deciding.
- **uBAM output.** `clump_only_single_to_bam` / `..._paired_to_bam_one_pair` share the layout and
  buffer BGZF differently. One cell each, or state explicitly that they are out of scope.

## I2. Several matrix corners are not measurable on a 16 GiB host

Budgets in the matrix are derived from the target B, so the high-core/high-B corners get large:

| cell | required `--memory` | predicted actual peak (pooled model) |
|---|---|---|
| c8 × 56.9 MiB | 3585 MiB | ~4.8 GiB |
| c16 × 37.9 MiB | 4605 MiB | ~6.3 GiB |
| c16 × 56.9 MiB | 6657 MiB | **~9.3 GiB** |

On the 16 GiB Mac the code itself references (`src/clump.rs:894`), a ~9.3 GiB peak alongside an IDE
and agent processes will hit the memory compressor — and compressed pages leave the resident set, so
`ru_maxrss` **understates** the peak. Those cells would silently return flattering ratios: exactly
the false-negative class that has bitten this workstream before. Either cap the matrix at ~4 GiB
actual, or record `vm_stat` compressor pages alongside each cell and discard any run where they
moved.

## I3. The promise is per-process but every measurement is per-pair

`resolve_clump_layout` is called once per input file inside each loop (`src/main.rs:1344`, `:1580`,
`:1859`), and `--memory` reads as a promise about the process. All measured cells are single-pair
runs. A user trimming 10 pairs in one invocation carries each pair's allocator retention into the
next pair's high-water mark, so the process peak exceeds any single-pair measurement. Add one
multi-pair cell (e.g. 3 pairs at the constraint cell) and, if it drifts upward, say so in the docs
rather than discovering it in an issue.

## I4. The Tests section is closer but still miscounted, and one listed test does not break

Checked each against the pooled-envelope coefficients (f0 = 244 MiB, σ = 1.471, k = 5.110):

| test | line | outcome |
|---|---|---|
| `resolve_layout_default` (4G/c2, expects 150–156 MiB) | `src/clump.rs:445` | **fails** — new B = 114.1 MiB |
| `resolve_layout_scales_with_cores` (4G/c8, expects 65–70 MiB) | `:455` | **fails** — new B = 43.8 MiB |
| `resolve_layout_bails_when_too_small` (560 MiB/c8) | `:483` | **fails** — new floor 332 MiB, so 560 MiB now resolves; the `assert!(res.is_err())` at `:487` fires before the `unwrap_err()` at `:488`, so it fails cleanly |
| `"static overhead"` substring | `:490` | **fails** on any bail-message rewrite (an assertion inside the test above, not a sixth test) |
| `resolve_layout_bails_below_static_overhead` (256 MiB/c4) | `:494` | **passes, by 32 MiB.** New floor 288 MiB > 256. It calls `unwrap_err()` with no preceding `is_err()` guard, so it panics if the fixed term lands below ~212 MiB |
| `resolve_layout_predicted_peak_fits_budget` | `:465` | **passes unchanged** — the plan lists it as affected and then correctly says it is a tautology. Contradiction within the same section |

So: **4 tests break plus 1 assertion**, not six tests. Two further groups are floor-sensitive and
unlisted:

- `tests/integration_preflight_tripwire.rs:581-593` runs `--clumpify --cores 2` at the default
  `--memory 1G` and asserts `clumpify:` appears in stderr — it will fail if the floor ever rises
  above 1 GiB at c2. (This is also the existing guard the plan correctly cites for its new banner
  cell.) `tests/integration_partial_output.rs:138-140` likewise runs `--clumpify --cores 2`.
- 14 call sites in `src/clump_only.rs` use `small_memory_budget()` = 1 GiB
  (`src/clump_only.rs:1187`); none passes `fastqc_requested = true`, so a FastQC-raised floor does
  not reach them. Worth stating as verified rather than left implicit.

Also: the tautology only *stays* a tautology if Phase 1 inverts exactly the same expression it
predicts with. If a per-core additive term is rounded differently in the two directions, that test
becomes a real check — cheap and worth keeping in that case.

## I5. `1G --cores 4` is not "the default cell", and pinning the constraint there protects the wrong thing

`--cores` defaults to **1** (`src/cli.rs:351`) and `--clumpify` requires `>= 2`
(`src/cli.rs:734`), so there is no default core count for clumpify at all; the docs' own worked
examples use `--cores 6` and `--cores 8` (`docs/.../clumpy.md:102, 146`). Measured ratios vary
sharply across that space — 0.85 at c2/1G, 1.03 at c4/1G, 0.94 at c8/1G — so a constraint anchored
to one cell (Goal, Validation row 3) locks the one point where the two model errors happen to cancel
and leaves the rest free to move. Under the pooled fit, B at 1G moves +17% at c2 and +4% at c8/c16
— i.e. the *shape* of the correction is not "hold everything, fix nothing".

Better constraint: bound the **compression outcome** across the grid ("total output bytes no more
than X% above baseline at `{c2,c4,c8} × {1G,2G}`"), which is what the Goal actually cares about, and
let B move where the data says it should.

## I6. The docs section under-scopes: the published savings table *will* go stale, and there are five unit labels

- The plan argues the savings table is "only stale if Phase 1 moves bin size, which the Goal forbids
  at the default". But the Goal only fixes the default; the table's headline claim is at the other
  end — `docs/.../clumpy.md:42` promises RRBS paired "~24% at default `--memory 1G`; up to ~31% at
  `--memory 4G+`", and `:61` advertises `--memory 4G` for archival. Under the envelope model, B at
  4G/c6 drops from 88.5 → 58.4 MiB (−34%), which is precisely where that "up to ~31%" was earned.
  The high-memory rows go stale under *any* honest recalibration. Say so, and plan the re-run (or
  soften the claim to a range).
- "Fix the label in the same change" is one site in the plan and **five** in the tree:
  `src/main.rs:203` plus `src/clump_only.rs:281, 420, 760, 999`, all printing
  `bin_byte_budget / (1024*1024)` under a `MB` label. The floor warning at `src/main.rs:189` is
  already correctly `MiB`.
- Minor: the docs' floor lower bound "~535 MiB" (`:108`) is the c2 value; c1 (reachable via
  `--clump_only`) is 533 MiB. Harmless, but it is in the same paragraph you are rewriting.

---

# Optional

## O1. Derive the spine factor at runtime instead of baking a constant

Both σ and k are ≈ the spine factor (C4), which is `1 + 148/text_bytes_per_pair` — a function of
read length, i.e. of the input. `FastqReader::sanity_check` already reads the first records before
sizing happens, so `resolve_layout` could take an observed `bytes_per_record` and compute the factor
instead of hard-coding one fixture's value. That removes Assumption 4's whole risk class (a model
tuned to one read length) rather than mitigating it with more fixtures, and it makes the banner's
number honest for ONT/amplicon input too, where the current constants are wildly conservative.

## O2. Retire `resolve_layout_predicted_peak_fits_budget` and put the assertion where it bites

It cannot fail (C-verified above), and it reads as coverage. Replace it with the regression test the
plan already proposes (pin resolved B at named `(budget, cores)` cells) plus one property test that
`clumpify_min_memory_bytes(c)` is the smallest budget for which `resolve_layout(_, c)` succeeds —
that one *can* fail, and it is exactly the sizing/floor disagreement Step 5 worries about.

## O3. `--clump_only`'s `cores` question has a third answer

Step 7 frames it as "pass `cores = 1` or the user's `cores`". There is a third: give the clump-only
path its own layout function, since its peak shape genuinely differs (no queues, one member
buffered, `src/clump_only.rs:213-227`) and `--cores` there only feeds `fastqc::run`. Sharing
`resolve_layout` between two different memory topologies is what makes Step 7 a dilemma at all.

---

# Answers to the five questions posed

**1. Does the sweep support v2's conclusions?** The arithmetic is right (every number reproduced);
the *inference* is not. R² of 0.997/0.9996 makes the linear form look solid, but that is 5 B-values
with 1–2 reps and residual patterns that wander (c4: −18, −12, +19, +24, −5, 0, −8), so the fit is
precise about the mean and silent about the maximum, which is the quantity being promised (C2).
Intercept extrapolation to B = 0 is only ~0.3 range-widths beyond the smallest point, so that is not
the problem — the problem is that the per-series intercept and slope are not separately identified
(C1), and the plan's headline claim ("today's B is within 1% of the honest maximum") asserts a
0.85% margin inside a 2.3% same-cell spread and a 4.9% rep spread. The claims that survive: the
model's decomposition is wrong; the breach grows with budget; compression pays for bin size;
`--compression` is irrelevant to RSS.

**2. Can Phase 0's matrix identify the model?** As a superset, yes for σ/k — but it spends ~10 runs
in the region where σ and k are perfectly collinear (`cores ≥ 4`), omits the one cell that makes the
form *testable* (`cores = 3`), omits read length despite Assumption 4 requiring it, and omits
`--clump_only` and the FastQC phase despite Steps 4/5/7 depending on them (I1). Runtime is not the
objection — the 18 sweep cells took ~4 minutes, so 80 runs is ~20 minutes — the objection is that
the cells cannot answer the questions the plan poses. And four corners are not measurable on a
16 GiB host without the memory compressor corrupting the reading (I2).

**3. Is the central constraint right?** No. Fixing B at `1G --cores 4` while promising ratio ≤ 1.0
is self-contradictory (C3), it is anchored to a cell that is neither the default nor representative
(I5), it is stated to a precision the data cannot support, and it becomes false as soon as a 51 bp
fixture enters the calibration (I1). Constrain the compression *outcome* across the grid instead,
and let B move.

**4. Should the structural fix be primary?** Yes — but not the one nominated. The queue depth
(`sync_channel(2)`, three live bins per worker against 1.25 charged) and the 300 B/pair output
reserve are the per-core coefficient, mechanistically and arithmetically (C4, and `FINDINGS.md:61-78`
already had it). Cutting the depth to 1 takes `k` from ~4.2 to ~2.8, which holds B at 19.2 MiB
(depth 1) or raises it to 23.1 MiB (depth 0) at 1G/c4 *with* the envelope and a 10% margin already
paid, against 15.6 MiB if the depth stays as-is — dissolving the Goal/Validation conflict instead of
negotiating it. Buffer recycling (Step 9) is the weakest candidate: ~400 large cross-thread
free/alloc pairs per run is not an allocator stress test, mimalloc already fixed the retention
symptom, and recycling touches none of the three real peak components. Keep it as an experiment,
not the headline.

**5. Code-level claims.** All verified correct: the struct-literal constraint, `fastqc_requested()`,
the four `clump_only.rs` call sites with `fastqc_requested` in scope, the tautological test, the
docs' wrong `12 MB` worked value, the floor/warning/LaTeX inventory, and the `--fastqc --memory 570M
--cores 8` demotion example (today's c8 floor is 566 MiB, so it resolves by 4 MiB). Corrections:
the Tests section is 4 tests + 1 assertion, not six (and it lists the tautology as affected while
also explaining that it is not); `resolve_layout_bails_below_static_overhead` survives by only
32 MiB; two integration tests and 14 `small_memory_budget()` sites are floor-sensitive and unlisted;
the unit-label fix is five sites, not one (I4, I6).

**6. Validation, adversarially.** Failure modes that pass all 9 rows and still ship a broken promise:

- **Fitted-to-the-mean.** Half the cells breach on the fitting data itself; if the fit is instead
  eyeballed to pass, rows 1–2 are checked against the same runs that produced the coefficients (row
  2's held-out cells are the only real defence, and 3 cells cannot cover a 2-D space).
- **Peak is a maximum of a stochastic process.** Queue occupancy varies run to run (4.9% observed);
  with 2 reps the observed max understates the true max, so a promise "verified" today breaks for a
  user who runs the same cell twenty times. Either raise reps on constraint cells to ≥ 5 **and**
  verify against the max, or remove the randomness by bounding the queues (C4).
- **Memory-compressor false negative** on the large corners (I2) — RSS drops out of `ru_maxrss`
  under pressure and the cell reads as a pass.
- **Fixture too small.** Rows 1/2/5/6 do not specify an input size. On a 10K-record fixture the bins
  never fill and every ratio is ~0.1. Specify a minimum (the 6.81 M-pair fixture, or ≥ 30 bin
  flushes per bin).
- **Read length.** Every row can pass at 65 bp and breach at 51 bp by ~8% (I1).
- **`--clump_only` / uBAM-out / multi-pair** are untested by any row, and all three use the same
  model (I1, I3).
- **Compression regression at large budgets** (row 4 does not name budgets): the honest model cuts B
  ~28% at 2G and ~34% at 4G/c6, so output grows ~15% at 2G on the sweep fixture while the default
  cell shows ~1%. A row-4 check run only at 1G passes and the docs' "up to ~31% at 4G+" claim
  silently becomes false (I6).

---

# Action items

## Critical

1. **Re-derive the model from the mechanism, not the curve.** Adopt `f0 + B·(σ·n_bins + k·cores)`
   with `σ ≈ spine_factor`, `k ≈ (work-channel depth + 1)·spine_factor`, and `f0 ≈ 250 MiB`
   (FastQC excluded). State that the pooled fit over all 19 existing cells already gives
   f0 = 244, σ = 1.33, k = 4.61 with rmse 4% *including* the c8 cell, so no per-core fixed term is
   warranted. Delete measurement 2 and the "two core counts cannot identify three parameters"
   premise. (C1)
2. **Make the promise one-sided on purpose.** Fit the upper envelope over reps, add an explicit
   named safety margin (~10%), and write the resulting ratio target into the plan. Say which
   quantile the banner's number means. (C2)
3. **Resolve the Goal/Validation contradiction before Phase 0 starts.** Pick: reduce actual memory
   (preferred), accept a bin-size cut at the default, or re-word the promise. As written, rows 1 and
   3 cannot both pass. (C3)
4. **Promote Steps 8 and 10 to Phase 1 primary work**, gated on clumpy mode; demote Step 9 to an
   experiment. Run the two discriminating measurements *first* — (a) work channel depth 2 → 1 at
   `1G/c4`, (b) output reserve sized from `raw_bytes` instead of `len·300` — each is a one-line
   patch and two runs, and together they decide whether any recalibration is needed. (C4)
5. **Re-specify FastQC as a `max`, not an addend**, and re-measure its cost under mimalloc. Decide
   whether the promise covers the FastQC phase before writing Validation row 5. (C5)

## Important

6. Redesign Phase 0's cells: add a **`cores = 3`** B-sweep (the only cell that makes the form
   testable), trim c8/c16 to 2–3 bin sizes each, add **51 bp and 150 bp** fixtures as a factor, add
   `--clump_only` and uBAM-out cells, add one multi-pair cell. Net: ~20–25 runs per platform doing
   strictly more than the proposed ~80. (I1, I3)
7. Cap the matrix at ~4 GiB predicted actual peak, or record compressor activity per cell and
   discard runs where it moved. (I2)
8. Fix the Tests section: 4 tests + 1 assertion; drop the tautology from the affected list; note
   `resolve_layout_bails_below_static_overhead`'s 32 MiB margin and add an `is_err()` guard before
   its `unwrap_err()`; list `tests/integration_preflight_tripwire.rs:581` and
   `tests/integration_partial_output.rs:138` as floor-sensitive. (I4)
9. Re-anchor the constraint from "B at one cell" to "output bytes across `{c2,c4,c8} × {1G,2G}`",
   and stop calling `1G --cores 4` the default. (I5)
10. Expand the docs scope: five `MB`→`MiB` label sites; plan for the 4G savings rows going stale;
    correct the `32 bins × 12 MB` worked value to 9. (I6)

## Optional

11. Derive the spine factor from the first records rather than a constant. (O1)
12. Replace `resolve_layout_predicted_peak_fits_budget` with a floor/sizing-agreement property
    test. (O2)
13. Give `--clump_only` its own layout function rather than answering Step 7's dilemma. (O3)

