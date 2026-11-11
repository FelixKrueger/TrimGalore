# Plan review A — honest `--clumpify` memory accounting (#439)

**Plan:** `plans/08162026_clumpify-memory-accounting/PLAN.md`
**Reviewer:** A (independent; no coordination with reviewer B)
**Tree state:** branch `dev`, tip `7ffe86f` — **not** `115d84f` as the plan's context states.
`7ffe86f` ("publish an output only if its teardown succeeded", #436) landed after the
plan's line references were taken and shifted `src/parallel.rs` below ~line 200.
**Raw data re-read:** `scratchpad/rss/*.tsv` + `bench_host.py` from session `720d3ced`
(still present on this host), so every measurement below is checked against the
primary data, not against `FINDINGS.md`'s summary of it.

---

## Verdict

The diagnosis is right and well evidenced. **The prescription does not follow from
the evidence.** The plan quotes the measurement that refutes its own central
constant (126–140 MiB per core at B ≈ 19 MiB) and then picks a coefficient worth
54 MiB at that bin size, without ever reconciling the two. It simultaneously spends
the 384 MiB reservation that `FINDINGS.md` identifies as the only thing currently
holding the promise on Linux. The net of "spend the cushion + under-correct the
per-core term" is a model that is *more* wrong than today's at the default
`--memory 1G --cores 4`, in the direction of breaking the promise rather than
keeping it.

Separately, the plan will not compile as specified (private field vs. an
out-of-module struct literal) and will fail two existing unit tests it does not
list.

Recommended action: do not implement as written. One cheap missing measurement
(a mimalloc bin-size sweep at fixed cores, 3 cells, ~2 min each) determines the
coefficient. Everything else in the plan is sound and can proceed.

---

## 1. Arithmetic — recomputed independently

All three headline numbers check out. `n_bins = max(16, 4·cores)`,
`usable = budget − static`, `B = 4·usable/(5·n_bins + 12·cores)`.

| cell | old B | new B (no FastQC) | plan claims | verdict |
|---|---|---|---|---|
| 1 GiB, 4 cores | `4·512/108` = **18.963 MiB** | `4·896/128` = **28.000 MiB** | 18.96 → 28.00 | ✅ exact |
| 4 GiB, 2 cores | `4·3584/94` = 152.51 MiB | `4·3968/104` = **152.615 MiB** | ~152.6 MiB | ✅ exact |
| 4 GiB, 8 cores | `4·3584/216` = 66.37 MiB | `4·3968/256` = **62.000 MiB** | 66 → ~62 MiB | ✅ exact |

- `resolve_layout_default` (clump.rs:445) window is `>= 150 MiB && < 156 MiB` →
  152.615 passes. ✅ The plan is right that this is a coincidence and right to
  assert the new arithmetic deliberately.
- `resolve_layout_scales_with_cores` (clump.rs:455) window `>= 65 && < 70` → 62.0
  fails. ✅ correctly flagged.
- `resolve_layout_predicted_peak_fits_budget` (clump.rs:465) passes for *any*
  coefficient: `B = floor(4u/d)` ⟹ `floor(B·d/4) ≤ u` ⟹ `predicted ≤ static + usable
  = budget`. ✅ passes — see §5 for why calling it "the load-bearing invariant"
  is the wrong reading of it.
- Fallback coefficient 3.75: `4·896/140` = 25.6 MiB. ✅ exact.
- Bin growth is 1.48× for cores ≥ 4 (`(1.75)·(27/32)`) and 1.58× at cores 2
  (`1.75·(94/104)`) — the plan's "~48%" is right and holds at every core count.
- Fallback floor moves from **539 MiB → 160 MiB** at 4 cores
  (`128 + ceil(1·128/4)`), and the whole range moves from 535–730 MiB to
  **154–256 MiB** (cores 2 → 16). The plan mentions the direction but never the
  numbers, which matters for docs (§4) and for §3.4.

**Not in the plan: with `--fastqc`, bins get *smaller*.** The static reservation
stays 512 MiB while the denominator grows, so
`B_new/B_old = 27/32 = 0.844` for cores ≥ 4 — a **15.6% bin shrink** on the
`--clumpify --fastqc` path (1 GiB/4 cores: 18.96 → 16.00 MiB). The Goal section's
claim that the change "moves the per-bin budget **up** … so this is expected to
*improve* clumpify's output size" is true only without FastQC. By the plan's own
elasticity figure (33% smaller bins → +28% output on duplicate-rich data), a 16%
shrink is a real compression regression on that path, and **no validation row
would catch it** (row 3 measures compression without `--fastqc`; row 7 measures
FastQC's memory but not its output size).

---

## 2. Facts verified against the code (what the plan got right)

Worth stating so the criticism below is not read as blanket doubt:

- `resolve_layout` clump.rs:196, `predicted_peak_bytes` :147, `clumpify_min_memory_bytes`
  :157, `STATIC_OVERHEAD_BYTES` :168, `MIN_BIN_BYTES` :38, `cores == 0` guard :197 —
  all correct.
- `Cli::fastqc_requested()` is cli.rs:592-594 and is exactly `self.fastqc ||
  self.fastqc_args.is_some()`. The `--fastqc_args`-only edge case is real and the
  predicate is the right one. ✅
- The four `clump_only.rs` `resolve_layout` sites are :275, :413, :748, :905 ✅, and
  **Assumption 6 can be closed now rather than hedged**: each enclosing function
  already takes a `fastqc_requested: bool` parameter (clump_only.rs:258, :396, :743,
  :893). No site needs to guess or pass `true` defensively.
- `--clumpify` has no Perl counterpart, so the byte-identity matrix is untouched. ✅
- Banner shape is preserved, which matters more than the plan says: the measurement
  harness parses it with a literal regex
  (`bench_host.py:20`: `clumpify: (\d+) bins × (\d+) MB; predicted peak ≈ (\d+) MB`).
  Any reflow silently zeroes out every validation row. Worth stating as a constraint.
- `PairedBin::push` counts both mates (clump.rs:250 via parallel.rs:896), so the
  budget is per-pair — `FINDINGS.md` checked and rejected this; confirmed. ✅

### Line references that are wrong and will misdirect an implementer

| plan says | actual |
|---|---|
| main.rs :178 and :195 (pass `fastqc_requested()`) | the two call sites are **:184** (`clumpify_min_memory_bytes`) and **:201** (`resolve_layout`); :178 is the `fn` signature |
| banner at main.rs:196 | :202-208 |
| single-end work channel parallel.rs:928 | **:943** (`sync_channel::<SingleWork>(2)`) |
| single-end result channel :933 | **:955** |
| reassembly maps :227, :972 | **:228**, **:998** |
| "`tests/` already covers this (`test_clumpy_single_end_is_a_permutation`, and the paired lockstep test at src/parallel.rs:1624)" | both are unit tests **inside src/parallel.rs** (:1636 and :1690); nothing in `tests/` covers clumpy permutation. :1624 is inside a helper |

---

## 3. Critical findings

### C1 — The plan will not compile: `ClumpLayout` gains a private field, but there is an out-of-module struct literal

`src/parallel.rs:1611-1616` builds the layout by literal:

```rust
fn small_clump_layout() -> ClumpLayout {
    ClumpLayout { n_bins: 16, bin_byte_budget: 1024 * 1024 }
}
```

The Signature section declares `static_bytes: u64` **without `pub`**. That is E0063
(missing field) *and* E0451 (private field in a struct literal outside the defining
module) — `cargo test` and `cargo clippy --all-targets` both fail. This is the only
out-of-module construction (verified: `grep -rn "ClumpLayout" src/ tests/`), and all
five clumpy tests in `parallel.rs` (:1636, :1690, :1793, :1839, :1913) route through
it, so the entire clumpy unit-test surface goes dark.

Fix: either make the field `pub`, or add a `ClumpLayout::for_test(n_bins, bin_bytes)`
constructor, or give the struct a `Default`-ish builder. Name the choice in the plan;
the "no way for the two to drift" rationale for privacy is still satisfied by a
constructor that derives `static_bytes` from a `fastqc: bool`.

### C2 — Two existing unit tests flip from `Err` to `Ok` and are not in the plan's list

The plan lists three tests to update and asserts the rest are unaffected. Two more
break, both because the floor moves (§1):

- `resolve_layout_bails_when_too_small` (clump.rs:483): 560 MiB, 8 cores. New model,
  `fastqc = false`: `usable = 432 MiB`, `denom = 256`, `B = 4·432/256 = 6.75 MiB` —
  well clear of the 1 MiB floor. `assert!(res.is_err())` **fails**. (It still bails
  if the test passes `fastqc = true`: `4·48/256 = 0.75 MiB`. So the test's outcome now
  silently depends on an argument the test author has to choose — that is worth a
  comment in the test itself.)
- `resolve_layout_bails_below_static_overhead` (clump.rs:493): 256 MiB, 4 cores. New
  model, `fastqc = false`: `usable = 128 MiB`, `denom = 128`, `B = 4 MiB`. The test
  does `res.unwrap_err()` → **panics**. Its premise ("below 512 MiB STATIC, usable
  saturates to 0") no longer exists as a concept.

Note the irony: the plan's proposed *new* test (c) — "a budget in the old dead zone
now resolves instead of falling back" — is the exact inverse of
`resolve_layout_bails_when_too_small`. The plan wrote the new test without noticing
it was contradicting an existing one. Also check the bail-message assertions
(`msg.contains("static overhead")`) survive step 3's rewording.

### C3 — The per-core coefficient contradicts the measurement the plan cites, by ~2×

This is the substantive objection.

Re-derived from `results_fixedbin.tsv` and `oxy_results.tsv` (bin size **and**
`n_bins` held fixed at 18 MiB / 16, only `cores` varying 2 → 4, so this is a clean
marginal):

| host / allocator | c2 actual | c4 actual | marginal per core | in units of B (18.87 MB) |
|---|---|---|---|---|
| macOS, mimalloc | 821.4 MB | 1073.8 MB | **126.2 MB** | **6.69 × B** |
| Linux, mimalloc | 794.6 MB | 1075.5 MB | **140.5 MB** | **7.44 × B** |
| macOS, system malloc | 1152.6 MB | 1452.2 MB | 149.8 MB | 7.94 × B |
| Linux, glibc | 705.9 MB | 672.3 MB | ≈ 0 (returns memory) | — |

Both platforms **agree** under mimalloc — the prerequisite allocator — at
**6.7–7.4 × B per core**, with work channels at capacity 2 (3 bins reachable per
worker). Removing one queued bin (step 8) subtracts at most `1.25 × B`, leaving
**≈ 5.4–6.2 × B**.

The plan's coefficient is **3.0 × B**. Its stated fallback is 3.75. Both are below
the measurement by ~2×. `FINDINGS.md`'s own option 1 recommended ≈ 5 per core
*before* the queue cap (`5·n_bins + 20·cores`); the plan silently landed on 12, not
20, and never says why it departed from the investigation's recommendation.

Where the 3.0 decomposition loses the bytes: the worker holds the input batch **and**
`buf_r1`/`buf_r2`, each pre-allocated at `reads_r1.len() * 300` bytes
(parallel.rs:349-351). For the 65 bp fixture (record 189 B, pair 378 B) that is
`300/378 = 0.794 × B` **each**, i.e. **1.59 × B of reserved output buffer per in-hand
batch** — not the 0.25 × B the plan budgets. Under glibc the untouched tail is never
faulted in, so only the compressed prefix is resident; under mimalloc (eager commit)
much more of it is. That single detail explains the platform split and predicts that
the output term gets *worse*, not better, under #437.

Adding it up structurally: `1.25` (queued) + `1.25` (in hand) + `1.59` (reserved out
buffers) ≈ **4.1 × B**, still short of the measured 6.7–7.4, with the residual
attributable to mimalloc per-thread retention and the `reserve_exact` re-reservation
in `PairedBin::push` (parallel.rs:897-902) that refills a bin's full spine the
instant it is flushed, while the flushed copy is still in flight.

**Consequence for the promise.** Correcting the harness unit bug (§F1), today's true
ratios under mimalloc on the diverse fixture are 1.003 (macOS c4), 0.918 (Linux c4),
0.941 (macOS c8), 0.955 (Linux c8) — i.e. the promise is already at its limit at
`c4`. The plan then adds 9.04 MiB to every one of 16 bins. **The resident bin pool
alone** grows by `16 × 9.04 × 1.25 = 181 MiB` (and by 201 MiB if the spine factor is
the fixture-correct 1.39, see I2). Ignoring every other term — no per-core growth,
no output growth, no change in static — that puts macOS at
`1026.5 + 181 = 1207 MiB` against a 1023 MiB promise (**1.18×**) and Linux at
`939.4 + 181 = 1120 MiB` (**1.09×**). Folding in the measured per-core response to B
puts both in the 1.2–1.35 range.

The promise is broken *by the fix*, at the default `--memory`, on both platforms,
with #437 landed. And the reason is exactly what `FINDINGS.md` warned about at
line 139: spending the static reservation "removes the very cushion currently masking
error 1, so it only makes sense together with option 1 or 2" — the plan does a
weakened option 1 *and* a partial option 2, which is not enough.

**What I cannot prove, and it matters.** The `n_bins` term and the per-core term
cannot be separated from the existing data, because `n_bins = max(16, 4·cores)` locks
them together for cores ≥ 4, and the only independent variation in the whole dataset
is the c2/c4 fixed-bin pair. There is a real possibility that the `1.25 × n_bins`
term is ~2× conservative (bins are independent random walks between flushes, so their
*simultaneous* fill averages ≈ B/2, not B) and that this conservatism is what has
been paying for the per-core error all along. If so, the honest model is
"smaller `n_bins` coefficient, much larger per-core coefficient", which could land
near the plan's 28 MiB — but for entirely different reasons, and with a different
core-count profile.

**Therefore: the missing measurement is a mimalloc bin-size sweep at fixed cores.**
There is not one cell of it in the dataset. The only B-sweep is system malloc
(`results_investigate.tsv`, c4, B = 18/56/132 → 1482/2916.5/4269 MB; the 18 MiB cell
is the mean of two measurements of the same config, 1512.0 and 1452.2, which also
tells you the run-to-run noise is ~4%), whose local slope is 37.8 MB per MiB-of-B
between 18 and 56 and 17.8 between 56 and 132 — versus the plan's modelled 32 and
today's 27. Three `tg_mi` cells at c4 with
`--memory 1G / 1.5G / 2G` settle the coefficient empirically in about five minutes,
and the harness already supports it (`bench_host.py` CELLS list). **Do this before
choosing 12.**

### C4 — Step 8 changes the queue depth for *every* `--cores N` run, not just clumpify

`parallel.rs:119` and `:943` are the per-worker channels for **both** routing modes —
`read_pairs_round_robin` (parallel.rs:665) and `read_pairs_clumpy` (:805) share them.
Capping them at 1 therefore changes jitter absorption for ordinary parallel trimming,
which is the overwhelmingly common path, has no memory promise to keep, and is the
subject of active performance work. In plain mode a batch is `BATCH_SIZE = 4096`
records (~1.5 MB), so the cap buys essentially no memory there — it is pure
throughput risk for zero benefit.

Validation row 4 measures throughput on "the same input" in a clumpify context only.
Either gate the capacity on `clump_layout.is_some()` (pass the capacity into the
channel construction), or add a plain-mode `--cores 8` throughput cell to row 4. As
written, a regression on the default path ships unmeasured.

(Minor: the Files-touched table calls this "Step 3"; the outline calls it Step 8.)

---

## 4. Important findings

### I1 — The docs page publishes the formula, the reservation, and the floor — all of which become *wrong*, not stale

The plan treats `docs/.../performance/clumpy.md` as a savings-table problem that
cannot be resolved without real libraries. That is only one line of that page.
The same page also publishes, in LaTeX, the thing this plan changes:

- `clumpy.md:93-100` — the full formula, including `usable = memory − 512 MiB` and
  the `5·n_bins + 7·cores` denominator.
- `:102` — worked examples: "with `--cores 8 --memory 1G` you get **32 bins × 12 MB**;
  with `--cores 8 --memory 4G` you get 32 bins × 66 MB". (The 12 MB figure is
  *already* wrong today — the banner prints 9 MB — so this page has rotted once
  before. New values: 14 MiB and 62 MiB.)
- `:108` — "needs a minimum of ~535 MiB to run (mostly from a fixed 512 MiB
  reservation) … stays in the 535–730 MiB range". New range: 154–256 MiB without
  FastQC, 538–640 MiB with.
- `:113` — a verbatim warning example, "need ≥ 552 MiB" at `--cores 6`. New value:
  176 MiB.

None of this needs a WGBS library. It is arithmetic, it is deterministic, and it must
land in the same PR or the docs actively mislead. The plan's open question should
shrink to the *savings table only*.

### I2 — `docs:104` contradicts the plan's central benefit claim

`clumpy.md:104`: "In theory, bigger budget → bigger per-gzip-member sort runs →
better compression. **In practice, increasing memory doesn't seem to make very much
difference in our tests.**"

The plan's justification for handing 384 MiB to the bin pool is that bigger bins
improve output size. The project's own published finding is that they mostly do not
(with RRBS at 4G+ the documented exception). Interpolating the published curve
(24% at B ≈ 19 MiB, 31% at B ≈ 132 MiB) puts B = 28 MiB at roughly 25–26% — i.e.
**~1–2 percentage points, ~2% smaller output**, not the "meaningfully smaller" that
validation row 3 expects. Row 3's stated expectation is also derived from a 9 → 6 MiB
datapoint, which is the steep end of a saturating curve and cannot be extrapolated to
19 → 28 MiB.

So the trade being proposed is: certain +181 MiB resident, for ~2% output. Row 3
should gate on "output no larger" only, and the plan should state the expected
magnitude honestly — because if the gain really is ~2%, the correct design is to
*keep bins where they are* and only fix the banner (see §6).

### I3 — Assumption 1's 1.25 spine factor is read-length dependent, and the validation fixture is at 1.39

Spine per pair is fixed at 148 B (`2 × 3 × 24` for the `String` headers + 4 B key);
text per pair is not. From the actual fixtures:

| fixture | header | seq | pair text | spine fraction | true factor |
|---|---|---|---|---|---|
| `SRR24766921_RRBS` (51 bp) | 54 | 51 | 322 B | 0.46 | **1.46** |
| `BS-seq_10K` (65 bp) — the huge fixture | 54 | 65 | 378 B | 0.39 | **1.39** |
| `clock_10K` (76 bp) | 58 | 76 | 430 B | 0.34 | 1.34 |
| typical 150 bp | ~54 | 150 | 718 B | 0.21 | 1.21 |

1.25 is right at ~150 bp and understates by 11–17% for the short-read libraries this
tool is most used on (RRBS, BS-seq, sRNA). The plan freezes it as "Fixed … the Linux
measurement fits it" — that measurement was taken on the 65 bp fixture, where the
true factor is 1.39, so the fit is absorbing the error elsewhere. Since the plan
*increases* the bin pool's share of the budget, this error gets multiplied.

Related and unmodelled: `PairedBin::push` (parallel.rs:895-907) sizes
`reserve_exact` from the **first pair pushed after each flush**
(`predicted = budget / bytes_this`). For heterogeneous read lengths the spine is
over-reserved by `median_size / first_size`, unbounded in principle. Every fixture in
the dataset is fixed-length, so this is entirely untested — and `--clump_only` is
precisely the mode aimed at already-trimmed (variable-length) FASTQ.

### I4 — Lowering the floor turns a safe, loud fallback into a silent over-budget run

The plan frames the dead zone opening as "the intended direction". Consider
`--memory 300M --cores 4` on a SLURM allocation with a hard `--mem=300M`:

- today: warning + fall back to plain mode, peak ~100 MB, job succeeds.
- after: `B = 4·172/128 = 5.375 MiB`, banner promises 300 MiB. Projected actual using
  the measured per-core term: `97 (static) + 16·5.375·1.25 (=107) + 4·6·5.375 (=129)`
  ≈ **333 MiB** → OOM-killed.

That is a strictly worse outcome for the user than the warning they get today, and it
lands in a range the plan is deliberately opening up. At minimum, row 6 of the
validation table should not be a unit test only — a newly-resolving budget needs one
real RSS measurement.

### I5 — The `--clump_only` path sizes bins for workers that do not exist

`cli.rs:779-782`: "`--clump_only` v1 is single-threaded internally … `--cores` is
accepted at any value >= 1". But all four sites pass the user's `--cores` into
`resolve_layout` (`cores.max(1)`), so the per-core term is fiction there. The plan
raises that term from 7 to 12 per core, which at `--clump_only --cores 8` means
**37.5% of the denominator** (`96/256`) is reserved for eight in-flight workers when
there is one. Passing `cores = 1` on those four sites (with `n_bins` chosen
separately) is a bigger and more certain win for that mode than the 128-vs-512 MiB
question the plan raises as its open item 4.

### I6 — The self-consistency of the coefficient is off even on the plan's own terms

Behavior step 3 derives `12 = 4 × (2 × 1.25 + 2 × 0.25)`, exactly. Assumption 2 then
says the plan "keeps that spirit by rounding the coefficient up rather than fitting
L6 exactly". It does not round up at all — 12 *is* the L6 fit. Meanwhile:

- **`--compression` defaults to 1**, not 6 (cli.rs:221-231; `test_compression_defaults_to_one`
  at :1744). The help text says L1 files are "75% larger".
- The existing code's doc comment (clump.rs:180-181, :187-188) deliberately uses the
  **L1** figure, `0.5 × B`, "as the upper bound".
- Every measurement in the dataset was taken at `--compression 6` (`bench_host.py:41`).

So the plan replaces a documented worst case with a best case, on the axis where the
shipped default is the worst case, and calls it conservative. Self-consistently it
should be `2 × 1.25 + 2 × 0.5 = 3.5` → `5·n_bins + 14·cores`. (The measurement says
higher still; see C3.)

---

## 5. Validation sufficiency

### F1 — The gate metric has a unit bug; "ratio ≤ 1.0" does not mean what it says

`bench_host.py:57-58`:

```python
rss_mb = ru.ru_maxrss * RSS_DIVISOR / 1e6     # true MB (10^6)
ratio  = rss_mb / float(pred)                  # pred is the banner's number
```

The banner divides by `1024*1024` (main.rs:205-206) and labels the result "MB", so
`pred` is **MiB**. The ratio is therefore `MB / MiB` and overstates the true ratio by
**4.86%** in every published cell. Corrected:

| cell (mi, diverse) | reported | true |
|---|---|---|
| macOS c4 | 1.05 | **1.003** |
| macOS c8 | 0.99 | 0.941 |
| Linux c4 | 0.96 | 0.918 |
| Linux c8 | 1.00 | 0.955 |

Direction is safe (the gate is 5% tighter than advertised), but a gate that is
miscalibrated in an unknown direction to its author is not a gate. Fix the harness
(divide by `1<<20`), and fix the banner itself — printing MiB while writing "MB" is a
user-facing inaccuracy in the very line this issue is about. Note this also revises
the plan's cited baselines: row 2's "today with #437: 0.98–1.05×" is really
0.92–1.00×, and Linux mimalloc c4 is 0.96 reported / 0.918 true, not 0.98.

### F2 — Failure modes that pass all eight rows

1. **Compression regression on `--clumpify --fastqc`** (§1): bins shrink 15.6% there.
   Row 3 has no FastQC; row 7 has no size check. Ships silently.
2. **Memory breach at the default `--compression 1`**: every row is at L6 or
   unspecified, and the harness hard-codes `--compression 6`. The output term is the
   one the plan shrank (0.5 → 0.25 × B), and L1 is where it is largest. Add a
   `--compression 1` cell at 1 GiB / 4 cores.
3. **Breach at `--cores 2`**: rows cover 4, 8, 12, 16. The per-core cost *in units of
   B* is worst at low core counts. Solving each mimalloc cell for the per-core term
   (assuming static = 97 MiB and spine 1.25, so these are model-dependent unlike the
   c2→c4 marginal in C3): **9.1 × B at c2, 7.9 at c4, 3.3 at c8, 3.5 at c12, 4.0 at
   c16** — the single reader thread cannot keep 8+ workers' queues full, so high-core
   cells self-limit and the plan's growth there is comparatively safe.
   `--clumpify` enforces `--cores >= 2`
   (`test_clumpify_requires_cores_at_least_two`, cli.rs:1626), so **c2 is the real
   worst case and it is untested**. `--clump_only` reaches `cores = 1`, worse again.
4. **Heterogeneous read lengths** (I3): both `reserve_exact` over-reservation and the
   spine factor are untested off fixed-length fixtures.
5. **Throughput regression on the plain `--cores N` path** (C4).
6. **Docs going factually wrong** (I1): nothing in the table checks the documented
   formula/floor/warning against the shipped values. A `--memory 100M --cores 6` run
   that emits a warning quoting a number different from `clumpy.md:113` would catch it.

### F3 — `resolve_layout_predicted_peak_fits_budget` is not "the load-bearing invariant"

It is a tautology of the shared formula (proof in §1): it passes for coefficient 7,
12, 15, or 1000. It passed continuously while macOS users were getting 2.8× their
budget. Treating it as the load-bearing check is how this defect survived. The
load-bearing check is rows 1/2/7/8 — measured RSS — and those are gated on hardware
and on #437, i.e. they cannot run in CI. That asymmetry deserves a sentence in the
plan: after this change, nothing in CI defends the promise, so the *measurement
must be recorded in the PR* and the coefficient's derivation must be in the commit
message where the next person can re-derive it.

Also: rows 1–4 and 7–8 all depend on regenerating a 2.6 GB fixture on two hosts, and
row 3's baseline must be captured *before* the change on the same fixture. Say so
explicitly — a "compare before vs after" row with no preserved baseline is not
runnable after the fact.

---

## 6. Alternative worth putting on the table

The plan considers exactly two shapes (its option 1 with a fitted coefficient, and a
macOS safety factor) and treats bigger bins as the natural consequence. There is a
third, and the data points at it:

**Fix the model, keep the bins.** Split the reservation as planned, *and* set the
per-core coefficient to the measured value (≈ 6 × B after the queue cap →
`5·n_bins + 24·cores`). Then at 1 GiB / 4 cores:

`B = 4 · 896 / (80 + 96) = 20.4 MiB`

— i.e. within 8% of today's 18.96 MiB. At 1 GiB / 8 cores: `4·896/352 = 10.2 MiB` vs
today's 9.48. Bins barely move, the banner becomes true rather than true-by-
construction, no compression regression is possible on any path (including
`--fastqc`, where the reservation is unchanged and the coefficient rise shrinks bins
— so cap that path's B at its current value or accept the shrink knowingly), and the
docs' savings table stays valid, which dissolves the plan's headline open question.

This is what the measurement actually says: **today's bin size is roughly right; only
the model is wrong.** The two errors have been cancelling almost exactly at the
default, which is why Linux has been fine. The plan's shape treats "the errors
cancel" as a coincidence to be broken in the direction of bigger bins; the safer read
is that the cancellation is load-bearing and should be replaced by a correct model
that lands in the same place.

Caveat, stated plainly: coefficient 24 still does not hold at `--cores 2` — `B` would
be 28 MiB there (`4·896/128`), while the fitted per-core at c2 is ≈9 × B, or ≈7.8
after the queue cap. No single constant fits all core counts, because the bottleneck
moves from the workers to the single reader thread as cores rise. Either fit the
coefficient at the worst core count (c2 → ≈31) and accept smaller bins everywhere, or
make the in-flight allowance `min(cores, k)`-shaped. That choice needs the mimalloc
B-sweep from C3.

---

## 7. Action items

### Critical
1. **Take the missing measurement before fixing the coefficient**: `tg_mi`, `--cores 4`,
   `--memory 1G / 1.5G / 2G` (B ≈ 19 / 29 / 39 MiB), plus one cell at `--cores 2`.
   Three-plus cells, existing harness, minutes. Without it the coefficient is a guess
   that the existing data contradicts (C3).
2. **Raise the per-core coefficient to the measured value** (≈ 6 × B post-cap →
   `5·n_bins + 24·cores`), or adopt §6 wholesale. 12 is refuted by the plan's own
   citation; 15 is too (C3, I6).
3. **Resolve the private-field/struct-literal break** at parallel.rs:1611 — name the
   mechanism in the plan (C1).
4. **Add `resolve_layout_bails_when_too_small` and `resolve_layout_bails_below_static_overhead`
   to the tests-to-update list**; both currently flip to `Ok` (C2).
5. **Do not cap the plain-mode work channel**, or benchmark it: gate capacity on
   `clump_layout.is_some()` (C4).
6. **Fix the harness unit bug** before any row can gate anything, and re-baseline the
   plan's quoted "today" ratios (F1).

### Important
7. State that `--clumpify --fastqc` bins *shrink* 15.6%, and add a compression check
   on that path (§1, F2.1).
8. Update `clumpy.md`'s formula, worked examples, floor prose and warning example in
   the same PR — these are arithmetic, not measurements, and are currently in the
   plan's "cannot re-measure" bucket by association (I1).
9. Reconcile row 3's expectation with `clumpy.md:104`; gate on "no larger" and state
   the expected magnitude as ~2% (I2).
10. Add validation cells at `--cores 2` and at `--compression 1` (F2.2, F2.3).
11. Re-derive the spine factor for the actual fixture (1.39 at 65 bp) or state that
    1.25 is a 150 bp assumption and quantify the shortfall for short-read libraries (I3).
12. Add one measured RSS cell in the newly-opened dead zone, and say in the plan that
    the change converts a safe fallback into a possible OOM under a hard cap (I4).
13. Close Assumption 6 — every `clump_only.rs` site already has `fastqc_requested`
    in scope (§2).
14. Fix the line references in §2 before handing this to an implementer.

### Optional
15. Pass `cores = 1` from `clump_only.rs`'s four sites, since that mode is
    single-threaded (I5).
16. Make the banner print MiB or divide by 10⁶ — the unit label is wrong in the exact
    line #439 is about (F1).
17. Bound the reassembly `BTreeMap` explicitly. It is *provably* bounded today
    (round-robin + capacity-1 channels means the reader can run at most ~2 batches per
    worker ahead of `expected`, so channel + map together hold ≤ 2·cores compressed
    batches — which is what the plan's 2 result slots per core already covers). Say
    that out loud in the plan rather than leaving "unbounded" as an open worry; it is
    the one risk area that turns out to be fine.
18. Note the `reserve_exact`-from-first-record hazard for variable-length input as a
    follow-up issue (I3).
19. Reconcile the "Step 3" / "Step 8" numbering for the parallel.rs change.
