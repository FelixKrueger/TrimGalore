# Plan review B — honest `--clumpify` memory accounting (#439)

**Plan:** `plans/08162026_clumpify-memory-accounting/PLAN.md`
**Evidence read:** `plans/08152026_clumpify-memory-overshoot/FINDINGS.md`
**Tree:** `dev` @ `115d84f`, `src/` and `docs/` clean (verified with `git status --porcelain`)
**Reviewer:** B (independent; no coordination with Reviewer A)

Units note: the code and the findings both print MiB but label it "MB". I follow
that convention below and mark it where the distinction changes a conclusion.

---

## Verdict up front

The diagnosis is right and the two-errors-in-opposite-directions framing is
correct. But the plan's **headline claim — that correcting both errors moves the
per-bin budget up 48% while making the banner honest — does not survive the
plan's own evidence.** The 48% gain comes almost entirely from spending the
384 MiB FastQC cushion, and the findings say in as many words that that cushion
"is precisely what keeps Linux under budget despite error 1"
(FINDINGS.md:97-98). The replacement per-core coefficient (3.0) is ~1.8× smaller
than the marginal per-core cost the same findings measured (§Coefficient below).
Spending a real 384 MiB cushion to buy a modelled 1.25 × B/core correction is a
net loss of margin, and my arithmetic says the post-change ratio lands at
**~1.2× on macOS+mimalloc and ~1.3× on Linux+mimalloc** — i.e. the promise is
*more* broken after this change than it is today with #437 alone.

An honest model calibrated on the measured marginal still yields bigger bins than
today (≈ +14%, not +48%). That is the change I would ship.

Separately, there are two things that will not compile or will fail on the first
`cargo test`, and a docs page that publishes the formula, the constant and the
floor range — none of which the plan accounts for.

---

## 1. Factual claims checked against the code

| Plan's claim | Verified? | Actual |
|---|---|---|
| `resolve_layout` at src/clump.rs:196 | ✅ | exact |
| `predicted_peak_bytes` at src/clump.rs:147 | ✅ | exact |
| `STATIC_OVERHEAD_BYTES` = 512 MiB at src/clump.rs:168 | ✅ | exact |
| `MIN_BIN_BYTES` at src/clump.rs:38 | ✅ | exact |
| `clumpify_min_memory_bytes` at src/clump.rs:157 | ✅ | exact |
| `cores == 0` guarded at src/clump.rs:197 | ✅ | exact (but dead code — see 4.6) |
| `Cli::fastqc_requested()` at src/cli.rs:592 | ✅ | exact; `self.fastqc \|\| self.fastqc_args.is_some()` |
| cli.rs test asserting `!cli.fastqc` on args-only | ✅ | assert is at src/cli.rs:2272 (plan says 2273) |
| paired work channel `sync_channel(2)` at parallel.rs:119 | ✅ | exact |
| paired result channel `cores * 2` at parallel.rs:128 | ✅ | exact |
| paired reassembly `BTreeMap` at parallel.rs:227 | ⚠️ | it is at **:228** (:227 is `let mut expected`) |
| SE work channel at parallel.rs:928 | ❌ | it is at **:943** |
| SE result channel at parallel.rs:933 | ❌ | it is at **:955** |
| SE `BTreeMap` at parallel.rs:972 | ❌ | it is at **:998** |
| main.rs call sites at :178 and :195 | ❌ | `clumpify_min_memory_bytes` is at **:184**, `resolve_layout` at **:201** (:178 is the `fn` line, :195 is inside an `eprintln!`) |
| four `clump_only.rs` `resolve_layout` sites at :275/:413/:748/:905 | ✅ | exact |
| stale 512 MiB comment at clump_only.rs:1188 | ✅ | exact (:1188-1189) |
| clump_only FastQC flows from "that mode's own `fastqc_args` parameter (:259)" | ⚠️ | :259 is `fastqc_args: Option<&str>`; the right thing to pass is the **`fastqc_requested: bool` param at :258** (and :396, :743, :893). Deriving it from `fastqc_args.is_some()` would miss plain `--fastqc`. |

The single-end `parallel.rs` references are all pre-`115d84f` — commit `115d84f`'s
#434 fix inserted the `Result<_>` comment block at :947-954, which is exactly the
offset. **Do not blind-edit by line number.**

Assumption 6 ("`--clump_only`'s FastQC state is knowable at each of its four call
sites") is **verified, not assumed**: all four functions already take
`fastqc_requested: bool` (clump_only.rs:258, :396, :743, :893) and it is in scope
at each `resolve_layout` call. Downgrade it from "Assumed" to a fact; the
`if any site cannot determine it, pass true` hedge is dead weight.

---

## 2. Arithmetic — recomputed independently

`n_bins = max(16, 4·cores)`; `B = floor(4·(mem − static) / (5·n_bins + k·cores))`.

| mem | cores | FastQC | n_bins | usable | denom (old→new) | B old | B new |
|---|---|---|---|---|---|---|---|
| 1 GiB | 4 | no | 16 | 512→896 MiB | 108 → 128 | 18.963 | **28.000** |
| 4 GiB | 2 | no | 16 | 3584→3968 | 94 → 104 | 152.51 | **152.615** |
| 4 GiB | 8 | no | 32 | 3584→3968 | 216 → 256 | 66.37 | **62.000** |
| 1 GiB | 8 | no | 32 | 512→896 | 216 → 256 | 9.481 | 14.000 |
| 1 GiB | 4 | **yes** | 16 | 512→512 | 108 → 128 | 18.963 | **16.000** |
| 1 GiB | 8 | **yes** | 32 | 512→512 | 216 → 256 | 9.481 | **8.000** |

Every number the plan asserts **holds**:

- 1G/4 cores: 18.96 → **28.00 MiB** ✅ (`4·896/128 = 28` exactly)
- 4 GiB/2 cores: **152.615 MiB** ✅, inside the existing 150–156 MiB window
- 4 GiB/8 cores: **62.000 MiB** ✅ exactly, so the 65–70 MiB window must change
- bins grow 28/18.963 = **+47.7%** ✅ ("~48%")
- step-8-reverted fallback (`5n + 15c`): 1G/4c → `4·896/140 = 25.6 MiB` ✅
- `predicted_peak ≤ budget` survives integer division: `B = floor(4u/d)` ⇒
  `B·d ≤ 4u` ⇒ `floor(B·d/4) ≤ u`. The load-bearing invariant test is safe. ✅

**But the 152.615 figure is only true for `fastqc = false`.** With `fastqc = true`
the same cell gives `4·3584/104 = 137.85 MiB`, which falls outside the existing
150–156 MiB window. The plan never says which value the updated unit tests pass.
Same problem in reverse for the two bail tests (§3.2).

### The FastQC path is a strict regression the plan does not mention

On the FastQC path the static reservation is unchanged, so only the denominator
moves: `108/128 = 0.844`. **Every `--fastqc --clumpify` configuration loses
15.6% of its bin budget.** By the plan's own row-3 calibration (33% smaller bins
cost +28% output on duplicate-rich data), a 15.6% shrink is worth roughly +13%
output size. Validation row 3 measures compression *without* `--fastqc`, so it
cannot see this.

### The fallback floor moves *up* on the FastQC path

`clumpify_min_memory_bytes = static + ceil(1 MiB · denom / 4)`:

| cores | today | new, no FastQC | new, FastQC |
|---|---|---|---|
| 2 | 536 MiB | 154 | **538** |
| 4 | 539 | 160 | **544** |
| 6 | 553 | 176 | **560** |
| 8 | 566 | 192 | **576** |
| 16 | 620 | 256 | **640** |
| 32 | 728 | 384 | **768** |

The plan's edge-case bullet says "At 128 MiB static the floor is easier to clear
than today". That is **false for every FastQC run**. Concrete regression:
`--clumpify --fastqc --memory 570M --cores 8` resolves today
(`4·58/216 = 1.07 MiB ≥ 1 MiB`) and **silently falls back to plain mode** after
the change (`4·58/256 = 0.906 MiB`). The plan asked whether any configuration
gets pushed somewhere worse — this is one, and it is silent (warning + exit 0).

---

## 3. Things that will not compile / will fail on first `cargo test`

### 3.1 CRITICAL — a private `static_bytes` field breaks `parallel.rs`

src/parallel.rs:1611-1616 constructs `ClumpLayout` with a struct literal:

```rust
fn small_clump_layout() -> ClumpLayout {
    ClumpLayout { n_bins: 16, bin_byte_budget: 1024 * 1024 }
}
```

`parallel::tests` is a different module from `clump`, so adding the **private**
field the plan's Signature block specifies (step 2: "Add the private
`static_bytes` field") is a hard compile error (E0451, plus a missing-field
error). This helper feeds `test_clumpy_single_end_is_a_permutation` and the
paired lockstep test — i.e. exactly the two tests the plan nominates as its
integration check in step 11, so the plan's own safety net is what breaks.

Fix: make it `pub(crate)`, or add a `ClumpLayout::for_tests(n_bins, budget)` /
`with_static(...)` constructor. Either is one line, but it must be in the plan.

### 3.2 CRITICAL — two existing tests the plan does not list will fail

`resolve_layout_bails_when_too_small` (src/clump.rs:483) and
`resolve_layout_bails_below_static_overhead` (src/clump.rs:494) both depend on
the 512 MiB reservation:

- 560 MiB / 8 cores with `fastqc = false`: usable 432 MiB, denom 256 →
  **6.75 MiB, no bail** → `assert!(res.is_err())` fails.
- 256 MiB / 4 cores with `fastqc = false`: usable 128 MiB, denom 128 →
  **4 MiB, no bail** → `res.unwrap_err()` **panics**.

They only survive if passed `fastqc = true`, and even then their comments'
arithmetic (`560 − 512 = 48 MiB`, `denom = 216`) is wrong. Additionally
`resolve_layout_bails_when_too_small` asserts `msg.contains("static overhead")`
— step 3 rewrites that bail message, so the reword must keep the phrase or the
test must change with it. Step 9 lists only two tests to update and declares one
other passes unchanged; it needs to enumerate all six touched tests.

### 3.3 IMPORTANT — new unit test (a) is not implementable as written

Step 10(a) wants "FastQC-on vs FastQC-off at the same budget yields the expected
384 MiB difference in usable bytes". `usable` is a local; `static_bytes` is
private; and `predicted_peak_bytes` is *equal* in both cases (both ≈ budget), so
the difference there is **0**, not 384 MiB. The only observable is
`bin_byte_budget`, whose difference is `4 · 384 MiB / denom` — at 1G/4 cores
`28 − 16 = 12 MiB`. State the assertion that way.

---

## 4. The per-core coefficient does not reconcile with the measurement

This is the heart of the review.

### 4.1 What the findings measured

FINDINGS.md:66-74 pins `B ≈ 18–19 MB` and `n_bins = 16` and varies only cores
(2 → 4), giving a **marginal** per-core cost:

| | per core | as multiple of B (19 MB) |
|---|---|---|
| old model (`1.75·B`) | 33–35 MB | 1.75 |
| **new model (`3.0·B`)** | **57 MB** | **3.00** |
| measured, macOS mimalloc | 126 MB | **6.63** |
| measured, macOS system malloc | 150 MB | 7.89 |

Dropping the work channel from 2 → 1 removes at most one in-flight bin,
`1.25·B ≈ 24 MB`, from the measured figure: **~102 MB ≈ 5.4·B**. The plan's 3.0
covers 56% of that. Its stated fallback of 3.75 covers 69%. **Neither closes the
gap**, and the plan presents 3.0 as derived-from-first-principles rather than
as a value it knows to be ~1.8× under the only marginal measurement taken.

Independent cross-check on a second cell (macOS mimalloc, 1G/8 cores, today):
predicted 1024 MB, ratio 0.99 ⇒ actual ≈ 1014. Subtract measured static
(101 MB) and the pool (`32 · 1.25 · 9.481 = 379`): 534 MB over 8 cores =
**66.8 MB = 7.0·B**. The 6.6–7.0 range is stable across cells.

### 4.2 The projected post-change ratios breach the promise

Using `actual ≈ static_real + 1.25·n_bins·B + 5.4·cores·B` (5.4 = measured
marginal minus one in-flight bin), at `--memory 1G`, no FastQC, macOS+mimalloc:

| cores | B new | predicted | projected actual | ratio |
|---|---|---|---|---|
| 4 | 28 MiB | 1024 | 101 + 560 + 605 = **1266** | **1.24×** |
| 8 | 14 MiB | 1024 | 101 + 560 + 605 = **1266** | **1.24×** |

Today, with #437 and nothing else, those cells read 1.05× and 0.99×. **This plan
takes macOS from ~1.0× to ~1.24×.**

Linux is worse, because #437 has already spent the glibc margin. FINDINGS.md:108-111
records Linux+mimalloc at **1.00×** at `--cores 8` today. Decomposing that cell
the same way (static 57, pool 379) leaves 73.5 MB/core = 7.75·B, so post-change:
`57 + 32·1.25·14 + 8·6.5·14 = 57 + 560 + 728 = 1345` against 1024 → **~1.31×**.
Linux/glibc alone stays inside (≈0.94×), but #437 makes mimalloc the shipped
allocator on both platforms, and Linux is where SLURM `--mem` actually kills the
job (the findings say exactly this at :116-117).

The plan quotes the findings' warning that the Linux margin "goes to zero" under
#437 (Context, line 47-48 of the plan) and then proceeds to spend 384 MiB of
reservation anyway. That is the logical gap: **the two errors were not merely
"partly cancelling" — the cushion is larger than the modelled correction, so
removing the cushion while under-correcting the coefficient is strictly worse
than doing neither.**

### 4.3 An honest model still improves bins — by ~14%, not ~48%

Solve for the bin size a *calibrated* model gives at 1G/4 cores, static 128 MiB,
coefficient 5.4: `denom = 80 + 21.6 = 101.6` → `B = 4·896/101.6 = 35.3`… no:
the coefficient multiplies `cores` in the pre-divided form, so
`denom = 5·16 + 4·5.4·4 = 80 + 86.4 = 166.4` → **B = 21.5 MiB**, i.e. **+13.5%**
over today's 18.96 MiB, with a banner that is actually true. If you also keep a
safety-honest 192 MiB static (measured 101 MB macOS is uncomfortably close to
128), you land at `4·832/166.4 = 20.0 MiB`, still +5%.

So the honest outcome is a modest bin increase, not 48%. The plan's Goal
sentence ("Correcting both moves the per-bin budget **up** … 28 MiB vs 18.96 MiB
… this is expected to *improve* clumpify's output size while making the banner
honest") is only reachable by leaving the model dishonest. **Rewrite that claim
before implementation, or the implementer will treat a validation failure at
row 1/2 as a surprise rather than the predicted outcome.**

### 4.4 Where the plan's decomposition double-counts, and where it under-counts

Round-robin dispatch is strict and blocking (parallel.rs:820-827:
`idx = *next_worker`, blocking `send`, then `(idx+1) % len`). That has two
consequences the plan misses:

- **Double-count (plan is conservative here):** a batch is either an input bin
  (`1.25·B`) or a compressed result (`~0.25·B`), never both. With strict
  round-robin the reader cannot get more than `capacity + 1` batches ahead of the
  slowest worker, so the total in-flight count is bounded by
  `cores · (capacity + 1)` = `2·cores` at capacity 1 — bounding *everything*
  downstream, including the "unbounded" `pending` BTreeMap. So `2 · 1.25 = 2.5`
  per core is a genuine upper bound on the *logical* live set, and the separate
  `+0.5` for result slots is slack, not additional. Good news for the derivation
  — and it also means the unbounded BTreeMap needs no explicit cap (worth saying
  so in the plan rather than listing it as an unbounded risk).
- **Under-count:** the worker *also* holds the output buffers for the batch in
  hand while building them (parallel.rs:349-351 allocates **two** `Vec`s of
  `reads.len() * 300` capacity each), and the single reader can be blocked inside
  `send` holding one further bin. `+0.25/core` and `+1.25` global.

Net: the *logical* bound is ≈2.75–3.0 per core, so **3.0 is a correct bound on
live bytes and still ~1.8× below measured RSS.** That is the tell: the missing
memory is not live bytes, it is allocator behaviour (see 4.5). No coefficient
derived from channel capacities can capture it.

### 4.5 The likely mechanism — and why it points at a different fix

`PairedBin::take()` (parallel.rs:913-919) `mem::take`s all three Vecs, so the
buffers travel to the worker and are **freed on the worker thread**, while the
next `push` (`keys.capacity() == 0` at :897) **re-allocates full-budget
capacity on the reader thread**. Every bin flush is therefore a cross-thread
allocate-here/free-there of three large buffers.

That is the worst case for both allocators under test: libmalloc retains
(the findings' stated mechanism, :100-106), and mimalloc's per-thread heaps
defer reclamation of cross-thread frees. It explains why the marginal cost
scales with **both** `B` and thread count at ~5-7× `B` with no corresponding
live bytes, and why it is ~2× larger on macOS than on glibc.

**Alternative worth putting in the plan: recycle the bin buffers.** Give each
worker a return channel that ships the drained `Vec`s back to the reader for
reuse (or hand the reader a small free-list). Allocation then happens
`n_bins + 2·cores` times for the whole run instead of once per flush, the
cross-thread free disappears, and the memory model becomes true *by
construction* rather than fitted per allocator. This plausibly fixes #439 on
macOS **without** #437 and removes the "3.0 vs 5.4" argument entirely. It is a
bigger change than the plan's step 8, but it is the only option on the table
that makes the promise structural. At minimum it belongs in the alternatives
the plan considered and rejected.

A cheaper partial: `cap = reads.len() * 300` (parallel.rs:349) is a fixed
300 bytes per pair regardless of read length — for 51 bp reads (RRBS, the very
library type the docs recommend clumpify for) that is **2.2 × B** of reserved
capacity per buffer pair per batch. Sizing it from `raw_bytes` would cut the
churn substantially at zero risk.

### 4.6 The compressed-output fraction is *not* rounded up — it is quietly relaxed

Plan assumption 2 says it "keeps that spirit by rounding the coefficient up
rather than fitting L6 exactly". Check the books:

- old: `cores × 0.5·B` for **one** output buffer at the **L1** worst case
- new: `2 slots × 0.25·B` = `cores × 0.5·B` for **two** buffers at the **L6** typical

The output allowance is **numerically identical** while covering twice as many
buffers, achieved by halving the per-buffer fraction and switching the
calibration point from worst case to typical. At `--compression 1` (which is the
**default** — `cli.rs` `--compression` defaults to 1, and the docs' own wall-time
table calls `--clumpify` compression 1 "essentially free", so it is the common
path) gzip L1 on FASTQ runs ~0.30–0.40 of raw, so the output term is understated
by ~1.5×. The plan should either state the fraction it is calibrating to and at
which level, or keep 0.5·B **per buffer** as the L1 bound.

### 4.7 The 1.25× spine factor is read-length dependent — and the plan calls it "Fixed"

`FastqRecord` is three `String`s = 72 bytes of spine per record (src/fastq.rs:58-65),
and `PairedBin::push` reserves `predicted` records in **both** `r1` and `r2`
(parallel.rs:898-901), i.e. `144 bytes per pair`:

| read length | pair text bytes | spine / text | implied factor |
|---|---|---|---|
| 2 × 150 bp (~350 B/record) | ~700 | 0.206 | 1.21 ✅ matches 1.25 |
| 2 × 51 bp (~137 B/record) | ~274 | **0.526** | **1.53** |

Assumption 1 lists 1.25 as "Fixed" and cites the Linux one-cell fit (695 vs 672)
as support. That fit cannot support it: it is one equation containing the pool
term, the in-flight term and the static term, none of which was measured
separately (FINDINGS.md:151-154 admits the per-core residual is attributed "by
elimination, not by instrumentation"). For 51 bp libraries the pool term — the
*dominant* term, 80 of 128 in the 1G/4-core denominator — is understated by ~22%.
Either derive the spine factor from the observed record size (it is available:
`bytes_this` at parallel.rs:896) or state 1.25 as an assumption that holds only
for ≥100 bp reads and size the reservation to absorb the short-read case.

---

## 5. Risk areas the plan flags, re-attacked

**Work channel 2 → 1 / head-of-line blocking.** The mechanism is real
(blocking send to a fixed worker index), but two mitigations the plan doesn't
mention make the throughput risk small: the reader thread does the minimizer
computation *and* `sort_paired_by_key` on the flush path (parallel.rs:819), so
the reader is very likely the bottleneck and worker queues are usually empty.
Also, blocking in `send` stalls reading entirely, so bins cannot overgrow while
blocked — no hidden memory risk. **The under-flagged risk is scope:** step 8
changes `sync_channel(2)` on the *shared* code path (parallel.rs:119, :943),
which also serves plain non-clumpify parallel trimming with `BATCH_SIZE = 4096`
batches (parallel.rs:33). Small, fast batches dispatched at a much higher rate
are where a capacity-1 channel is most likely to cost throughput, and validation
row 4 measures only clumpify cells. Either gate the capacity on
`clump_layout.is_some()`, or add plain-mode throughput to row 4.

**Unbounded reassembly BTreeMap.** Bounded implicitly by strict round-robin at
`cores · (capacity + 1)` batches in flight — see 4.4. Say so; don't leave it as
an open risk.

**Lowering the static reservation.** Two configurations get worse, both found
above: the FastQC path loses 15.6% of bin budget unconditionally (§2), and the
FastQC floor rises by up to 40 MiB, silently demoting a band of budgets to plain
mode (§2). Also note the 128 MiB choice has only 27 MiB of headroom over the
measured 101 MB macOS static, and the measured static was taken at `--cores 4`
only; `fastqc::run` receives `threads: cores.max(1)` (src/fastqc.rs:45) so the
FastQC component plausibly scales with cores while both reservations are flat.

**FastQC's 425 MB was measured on the allocator being replaced.** FastQC runs
*strictly after* trimming completes (src/main.rs:1481 for SE, :1797 for PE), so
the +425 MB delta in FINDINGS.md:88-94 is not concurrent demand — it is
retention of the trimming phase's memory across the phase boundary, and that
measurement is macOS **system malloc** (1512/1023 = 1.48× is exactly the
system-malloc cores-4 cell). Under mimalloc the two phases may not sum at all,
in which case 512 MiB on the FastQC path is over-reserved and is shrinking bins
for no reason. **The plan hard-codes the 512 MiB without re-measuring it under
the allocator it declares a prerequisite.** Validation row 7 checks only
`ratio ≤ 1.0`, which cannot detect over-reservation.

**`--clump_only`.** Verified: `clump_only.rs` contains **no threads at all**
(`grep -n "sync_channel\|spawn\|thread::"` → no matches), yet it sizes bins with
`resolve_layout(mem, cores.max(1))`, so the entire per-core term is fictional
there and raising it 7 → 12 charges that mode for queues it doesn't have. The
static drop still leaves it net-better at every cell I checked (1G/8 cores:
9.48 → 14.0 MiB), so this is not a regression — but it is an argument for a
`cores`-aware or mode-aware coefficient, and it makes the plan's open question
("is 128 MiB right for `--clump_only`?") the less interesting half of the
question. Note also that `--clump_only` propagates the bail with `?` (no
fallback), so its floor is fatal, not a warning — the plan's Behavior §5 says
"main.rs still warns and falls back", which is true only for `--clumpify`.

---

## 6. Validation sufficiency — what passes all 8 rows and still ships broken

1. **Rows 1, 2, 7, 8 pass most convincingly when clumpify silently doesn't
   run.** A configuration that falls back to plain mode peaks at ~100 MB and
   reports a gorgeous ratio of 0.1×. Given that this change *moves the floor*
   (up on the FastQC path — §2), that is not a theoretical risk. The repo
   already has the right guard: `tests/integration_preflight_tripwire.rs:591`
   asserts `stderr.contains("clumpify:")` for exactly this reason. **Every
   measured cell must assert the banner appeared**, and should record the
   resolved `bins × MB` so a silently-shrunk bin size is visible.
2. **Row 3 uses one duplication-rich fixture.** The project's own published
   guidance is that clumpify *hurts* WGBS and helps RRBS. Changing bin size
   moves that crossover. A row-3 pass on a 680×-duplicated or 3-library-cycled
   fixture tells you nothing about whether bigger bins make the already-negative
   low-duplication case worse. Add a low-duplication cell.
3. **No FastQC compression cell.** The one path that provably shrinks bins
   (−15.6%) is measured for RSS (row 7) and never for output size.
4. **All 8 rows are `--memory 1G` (or "same input") and paired-end.** No large
   budget, no single-end cell — despite step 8 changing the single-end channel at
   parallel.rs:943 and the single-end path having a different B semantics (one
   mate per record, not a pair).
5. **Row 4's ~3% threshold with 3 reps on a 128-core shared box** has no
   variance treatment; 3% is plausibly inside the noise on `dockyard-oxy-0`.
   Report a spread, and pick the threshold from the observed rep variance
   (`results_reps.tsv` exists from the investigation).
6. **Nothing validates the short-read regime**, where the spine factor is 1.53
   rather than 1.25 (§4.7). A 51 bp fixture at 1G/8 cores is the cheapest way to
   probe the pool term, which is the dominant one.
7. **No row confirms the docs' claims still hold** — see §7; row 3 is about
   bytes, not about the published formula.

Also: rows 1–4 and 7–8 "need the fixtures regenerated on both hosts". The
findings note (:147-149) that the oxy real-data fixtures were cleaned off the
box, so the whole matrix rests on synthetic/cycled input. That is honest but it
means **the plan cannot demonstrate the promise on any real library**, and it
should say so as a limitation of validation, not only as a docs problem.

---

## 7. Documentation impact is much larger than the plan states

The plan's Files-touched table and Integration section describe
`docs/src/content/docs/performance/clumpy.md` as having stale *savings figures*.
It publishes far more than that:

- **The formula itself, in LaTeX** (clumpy.md ~:95-105), including
  `usable = memory − 512 MiB` and `4 × usable / (5·n_bins + 7·cores)`.
- **Worked resolved values**: "with `--cores 8 --memory 1G` you get **32 bins ×
  12 MB**; with `--cores 8 --memory 4G` you get 32 bins × 66 MB". The 4G figure
  matches (66.37 MiB); **the 1G figure is already wrong today** — the shipped
  formula yields 9.48 MiB and the banner prints "9 MB", not 12. Worth fixing as
  part of this change.
- **The floor**: "needs a minimum of ~535 MiB … stays in the 535–730 MiB range
  for any sensible core count". After the change that becomes 154–384 MiB
  (no FastQC) / 538–768 MiB (FastQC) — one range becomes two.
- **A worked warning message** quoting "need ≥ 552 MiB" at `--cores 6`; the new
  value is 176 MiB (no FastQC) / 560 MiB (FastQC).

Leave `CHANGELOG.md:831-843` and its docs mirror alone — they describe released
v2.2.0 — but note that they publish "measured peak matches predicted peak within
1% on 31 M-record paired-end runs", which the findings have now falsified. The
new one-sentence entry should be phrased so a reader comparing the two is not
confused.

---

## 8. Smaller points

- **Second drift source left open.** The plan stores `static_bytes` so the banner
  and the sizing cannot disagree about the reservation, but
  `predicted_peak_bytes(&self, cores)` still takes `cores` as a parameter — so a
  caller can still evaluate the prediction against a different core count than
  the layout was sized for. Since you are touching the struct anyway, store the
  sized denominator (or `cores`) and drop the parameter; that closes the drift
  class rather than one instance of it.
- **`cores == 0` is dead code.** `--clumpify` requires `--cores >= 2`
  (src/cli.rs:732-734) and `--clump_only` passes `cores.max(1)`, so the guard at
  clump.rs:197 is unreachable from either path. Harmless; don't spend words on it
  in the plan's edge-case list.
- **`resolve_clump_layout` is called inside the per-file/per-pair loops**
  (main.rs:1344, :1580, :1859), so the banner prints once per file. Adding a
  `fastqc_requested()` argument is uniform across all three; the plan's "both
  :178 and :195" undercounts the sites (there are three `resolve_clump_layout`
  callers, one function to change — fine in practice, but the line references
  are wrong).
- **Plan step 11's integration check is the thing 3.1 breaks.** Reorder: fix the
  constructor first, then the tests are available as a net.

---

## Action items

### Critical

1. **Re-derive the per-core coefficient from the measured marginal, or change
   the plan's goal.** 3.0 is a valid bound on *live* bytes and ~1.8× below the
   measured marginal RSS (5.4·B with capacity 1; 6.6–7.0·B as shipped).
   Projected post-change ratios: ~1.24× macOS+mimalloc, ~1.31× Linux+mimalloc —
   both worse than today with #437 alone. Either use ≈5.5 (giving 21.5 MiB at
   1G/4 cores, +13.5% over today, honest banner) or state explicitly that the
   plan trades banner accuracy for bin size and drop the "makes the promise
   true" framing.
2. **Do not spend the 384 MiB cushion and under-correct the coefficient in the
   same change.** The findings say that cushion is the only thing keeping Linux
   in budget under error 1, and #437 has already taken Linux to 1.00×. If you
   want the cushion, land the calibrated coefficient with it.
3. **Fix the `ClumpLayout` constructibility break** — src/parallel.rs:1611
   builds the struct literally; a private field is a compile error. Add
   `pub(crate)` or a constructor, and do it before step 11's tests are relied on.
4. **Enumerate all six affected unit tests**, not two:
   `resolve_layout_default`, `resolve_layout_scales_with_cores`,
   `resolve_layout_bails_when_too_small` (fails outright with `fastqc = false`),
   `resolve_layout_bails_below_static_overhead` (**panics** with
   `fastqc = false`), `resolve_layout_predicted_peak_fits_budget` (needs the new
   argument), plus the `"static overhead"` substring assertion that step 3's
   message rewrite can break. State the `fastqc` value each test passes — the
   plan's 152.6 MiB claim is true only for `false`.
5. **Assert the banner appeared in every measured validation cell.** Rows 1, 2,
   7 and 8 all pass trivially if clumpify silently falls back, and this change
   moves the floor. Copy the guard from
   `tests/integration_preflight_tripwire.rs:591`.

### Important

6. **Own the FastQC-path regression.** −15.6% bin budget at every core count and
   a floor that rises up to 40 MiB (silently demoting e.g.
   `--fastqc --memory 570M --cores 8` to plain mode). Add a compression cell
   with `--fastqc`, and either lower the FastQC reservation in step with the
   denominator or say in the plan that the FastQC path pays for the correction.
7. **Re-measure the 512 MiB FastQC reservation under mimalloc.** FastQC runs
   strictly after trimming (main.rs:1481, :1797) and the +425 MB delta was
   measured on system malloc. Under mimalloc the phases may not sum, in which
   case 512 MiB is over-reserved. Also note it is flat in `cores` while
   `fastqc::run` gets `threads: cores` (fastqc.rs:45), and was only measured at
   `--cores 4`.
8. **Scope step 8 to the clumpify path or measure plain throughput.**
   `sync_channel(2)` at parallel.rs:119/:943 is shared with plain parallel
   trimming (`BATCH_SIZE = 4096`), which row 4 never measures.
9. **State the compression level the output term is calibrated to.** The new
   coefficient keeps the same `0.5·B` total while claiming to cover two buffers —
   that is an L1→L6 relaxation, not a round-up, and `--compression 1` is the
   default.
10. **Downgrade or qualify assumption 1 (1.25× spine).** It is 1.21 at 150 bp and
    1.53 at 51 bp (`144 bytes/pair`, parallel.rs:898-901, fastq.rs:58-65). The
    pool term is the dominant one; a 22% understatement on short-read libraries
    matters. Add a 51 bp validation cell.
11. **Fix the docs scope.** clumpy.md publishes the formula, the 512 MiB
    constant, worked bin sizes and the floor range — all wrong after this change,
    and its "32 bins × 12 MB" at 1G/8 cores is **already** wrong (shipped
    formula: 9 MB). This is not just a stale savings table.
12. **Refresh the stale line references** before implementation: parallel.rs
    SE sites are :943/:955/:998 (not 928/933/972), the paired BTreeMap is :228,
    main.rs sites are :184/:201, and the clump_only param to pass is
    `fastqc_requested` at :258/:396/:743/:893 — not `fastqc_args` at :259.

### Optional

13. **Consider bin-buffer recycling as the real fix** (return channel or reader
    free-list). Bins are allocated on the reader and freed on workers
    (parallel.rs:913-919 + :897-902), which is the pathological pattern for both
    allocators and the most likely source of the unexplained 5–7×`B` per core.
    Recycling makes the model true by construction and may fix #439 on macOS
    without #437. At minimum, record it as an alternative considered.
14. **Size `cap = reads.len() * 300` (parallel.rs:349) from `raw_bytes`.** For
    51 bp reads that is ~2.2×B of reserved capacity per batch, per worker, per
    buffer pair — pure churn, which is exactly what the findings say libmalloc
    retains.
15. **Store the sized `cores`/denominator in `ClumpLayout`** and drop the `cores`
    parameter from `predicted_peak_bytes`, closing the drift class rather than
    one instance.
16. **Say that no real library was measured** as a validation limitation, not
    only as a docs blocker (FINDINGS.md:147-149).
17. Assumption 6 is verified — all four `clump_only` sites already have
    `fastqc_requested` in scope. Drop the "pass `true` if unknowable" hedge.
18. Note in the plan that `--clump_only` is thread-free, so the per-core term is
    fictional there (net still better after the static drop), and that its bail
    is fatal rather than a fallback — Behavior §5 currently implies otherwise.

---

## What the plan gets right (probed and confirmed)

- Every arithmetic figure it asserts (28.00, 152.6, 62.0, 18.96, 25.6, +48%) is
  correct, and `predicted_peak ≤ budget` survives the integer division.
- `fastqc_requested()` is the right predicate, and the `--fastqc_args`-alone
  reasoning is exactly right (cli.rs:592, test at :2262-2279).
- Four independent `clump_only` call sites, each with the state it needs.
- The "3 in-flight bins per worker" reading of `sync_channel(2)` is right, and
  capacity 1 → 2 live bins is right.
- Rejecting the shrink-bins option on measured grounds (1.77× → 2.55×) was the
  correct call, and recording it in Self-Review is the right kind of artifact.
- Keeping step 8 as an independently revertible commit is good engineering.
- `--clumpify` has no Perl counterpart, so the byte-identity matrix is genuinely
  unaffected.
