# Plan review A — paired-end memory figures and cluster-allocation guidance

**Plan:** `plans/08182026_threading-docs-cluster-allocation/PLAN.md`
**Reviewer:** A (independent)
**Date:** 2026-08-18
**Verdict:** Do not implement as written. The plan's central factual premise is wrong, and the
structural change it proposes (a paired-end table beside "the single-end table") would publish two
mutually contradictory paired-end figures. The *goal* is right and the page-cache gap is real; the
narrative and the table shape need rebuilding on verified provenance first.

---

## 1. Logic review

### C1 — The existing memory table is already paired-end. The plan's cause chain is unsupported.

The plan's Goal says the published ~100 MB figure was "our figure, measured single-end on short
reads", and everything downstream (Behavior 2, step 3, the "add beside it" assumption, Open question 3)
follows from that. It is not single-end. Evidence, all in-repo:

- `docs/perf_data/buckberry-2026-04-29/memory_summary.txt`: *"Fixture: SRR24827373 (Buckberry 2023,
  84M **PE** reads) … Host: dockyard-oxy-0 … Method: /usr/bin/time -v, peak RSS, single run per
  condition"*, with exactly the four rows in `threading.md:53-58` (5.1 / 60.1 / 73.0 / 91.6).
- Every invocation in that matrix is `--paired`: `scripts/benchmark.sh:130` and the archived
  `scripts/run_bench_extra.sh:58,63` inside `logs-and-scripts.tar.gz` both read
  `--paired --cores {cores} -o {outdir} $R1 $R2`.
- The memory run *is* that run: its wall times match the hyperfine means within ~1% — c1 328.85 s vs
  `5:27.36`, c4 80.88 s vs `1:20.75`, c8 56.61 s vs `0:56.93` (`rust-beta7_c{1,4,8}.json`), which the
  perf_data README states as the tie-in.

Consequences:

- The Goal's arithmetic — *"Paired-end at `--cores 8` is roughly twice that, so the real headroom was
  ~5×"* — has no support. The published figure was paired.
- Behavior 2 / step 3 / the second Assumption would place a new paired table next to an existing
  paired table showing ~92 MB vs ~200 MiB for a nominally identical condition, with no explanation.
  That is worse for the reader than the current page.
- Open question 3 ("second table vs paired column") is moot as posed: there is no single-end/paired
  axis to split on. The real axis is fixture + platform + version.
- One thing this *does* settle cheaply: the plan's "To confirm before publishing" assumption about
  whether the existing figures are Linux. They are — dockyard-oxy-0, Xeon Granite Rapids, GNU time.

### C2 — Four existing plain-path paired short-read figures already disagree by ~2×, and the plan reconciles none of them.

| Figure | Condition | Source |
|---|---|---|
| 91.6 MB @ c8, 73.0 @ c4 | 84M pairs, 150 bp, real WGBS, gz in/out, Linux, beta.7 | `threading.md:53-58` + `memory_summary.txt` |
| 71–78 MiB @ c4 | 150 bp, 302 MB random fixture, gz in/out, `fdade7e` | issue #454, published at `threading.md:65` |
| 104 MiB @ c4 | 10M × 51 bp paired, random fixture | `clumpy.md:126` |
| 173 MiB (Linux/x86_64), 205 MiB (macOS/arm64) | 10M-pair paired runs | `src/clump.rs:250-252` (`STATIC_TRIM_BYTES` doc comment) |

Note the last row: the repo already holds a **Linux paired** peak-RSS figure of 173 MiB and a macOS
one of 205 MiB, the latter matching the plan's own macOS diagnostic (202 MiB @ c8) almost exactly. So
the plan's unpublishable macOS numbers are corroborated — and the gap to the published 91.6 MB is
*not* explained by platform (only ~1.2× there), nor by read length (51 bp scoring higher than 150 bp
inverts the record-count model the docs teach).

The most likely dominant factor is one the plan never mentions: **fixture compressibility.**
`examples/mk_clump_fixture.rs` emits fixed-seed random ACGT, which deflate barely compresses, so each
worker's output buffer and every slot in the `cores * 2` result channel holds ~4× the bytes a real
sample's would. The published 91.6 MB came from real WGBS reads; #454's and the clump model's figures
came from random fixtures. Publishing a fifth number without controlling for this gives the docs a
third contradiction, and if the new figures come from random fixtures they are an upper bound for real
data, not a measurement of it.

**Missing step:** a reconciliation run. Re-measure the exact published condition — the Buckberry pair,
whose path `scripts/benchmark.sh` still defaults to (`~/benchmark_TG_oxy/SRR24827373_*_R{1,2}.fastq.gz`)
— with today's binary. Two outcomes, and the plan has a branch for neither:

- It reproduces ~92 MB → the docs table is correct, the heap really is ~92 MiB paired at c8, and the
  #1903 story is *entirely* cgroup accounting. Behaviors 2 and 5 then need re-deriving, and a second
  table may not be warranted at all.
- It does not → there is a **memory regression** between beta.7 and master on the plain paired path,
  which is a code finding, not a docs finding, and should stop this plan rather than be papered over
  with a new table. (For what it's worth, the plain path's queue geometry is unchanged: beta.7's
  `sync_channel::<PairedWork>(2)` / `(cores * 2)` equals HEAD's `channel_depths(false)` = work 2,
  result_per_core 2, and `BATCH_SIZE` is still 4096. So a regression would be elsewhere — 2 250 lines
  of `parallel.rs` changed over that span.)

### C3 — Behavior 3 ("FastQC is not additive in memory") contradicts the repo's own calibrated model and a sibling docs page.

- `src/clump.rs:255-257`, verbatim: *"FastQC's additional resident cost per thread. Its phase follows
  trimming but **the pool's pages stay resident, so it adds rather than maxes**."* The constant is
  `STATIC_FASTQC_PER_THREAD_BYTES = 24 MiB`, charged up to `FASTQC_CHARGED_THREAD_CAP = 16` — i.e. up
  to **384 MiB** of additive FastQC charge in the calibrated model.
- `docs/src/content/docs/performance/clumpy.md:114` already publishes that additivity: floors of
  281→590 MiB without `--fastqc` and 334→1012 MiB with it, over `--cores` 2–32.
- `-t N` is not file-level parallelism. `src/fastqc.rs:35` documents `cores` as *"threading budget for
  `fastqc-rust`'s internal **rayon**"*, and `resolved_threads` lets `-t` override it, so the FastQC
  phase's own peak scales with `-t` regardless of how many files it reports on.

"Runs after trimming, therefore does not add" is the wrong mechanism to publish: allocator retention
is precisely why the project's own model charges it additively. Even if the plain path measures
`max()` rather than sum, the sentence as drafted will read as "`--fastqc_args -t 24` is free", and
nf-core runs FastQC in this process. The cluster recommendation must be sized *with* FastQC, or state
that it is not.

Validation 4 compounds this: its expected result ("difference small") points the same way as the
hypothesis, so it can only confirm. And the datapoint that motivated it is not a controlled comparison
— per the plan's own Assumptions, the 190 MiB `-t 24` figure was taken on 12M unique reads while the
202 MiB baseline came from the plateau runs. Compare at identical input, at both a low and a high `-t`,
and pre-commit to rewriting Behavior 3 if the delta tracks `-t`.

### C4 — The Goal states, as established fact, the conclusion step 1 is supposed to establish.

*"Paired-end at `--cores 8` is roughly twice that, so the real headroom was ~5×"* is written as the
premise of the whole plan, but rests on macOS figures the plan itself marks unpublishable plus the
mis-attribution in C1. Add an explicit decision rule up front: **if the Linux paired figures land near
the published table, the diagnosis is cgroup accounting and the deliverable shrinks to Behavior 6 + 7
plus a budgeting subsection.** As written, an implementer who measures ~92 MiB has no instruction other
than to publish a table that says the same thing the existing one does.

### I1 — Step 4 removes one invitation to the misreading and leaves four.

The plan names only `threading.md:63`'s *"well below the noise floor of cluster scheduler memory
allocation"*. Still standing after the planned edit:

- `:61` — *"Memory growth is bounded and predictable, **well-suited to cluster scheduling**: on these
  150 bp reads, even worst-case at the saturation point (`--cores 8`), the process stays under
  ~100 MB."* This is the sentence closest in content to #11541's rationale.
- `:63` (rest of the sentence) — *"stays under 100 MB across all reasonable core counts"*.
- `:80` — *"**Lower memory:** 5 MB single-threaded, ~10 MB per additional worker."* A bullet in
  "Beyond speed": the most quotable line on the page, and the one a module author lifts.
- `:84-88`, `## What this means for nf-core / Nextflow` — gives CPU guidance (*"If a Nextflow process
  has 12 CPUs, just pass `--cores 12`"*) and **no memory guidance at all**. This is the section the
  audience that produced #11541 actually reads.

Validation 5 (`grep -n "noise floor"`) passes with every one of those intact. So, answering the
question directly: as written, the plan mostly adds data alongside the sentences that invited the
misreading. Step 4 needs to enumerate all five sites, and Validation 5 needs to grep for `100 MB`,
`noise floor`, and `well-suited to cluster scheduling`.

### I2 — Cross-page contradiction with `clumpy.md`, already live since `0e86e2a`.

`clumpy.md:126` says *"Memory usage without `--clumpify` **depends on read length rather than on file
size**. Short-read runs sit around the 100 MB mark (104 MiB for 10 M × 51 bp paired at `--cores 4`)"*.
The first clause is exactly what #454 falsified this morning ("'not input file size' is wrong"), and
the second re-quotes the ~100 MB framing. The plan touches `clumpy.md` only to add a cross-link
(step 8), and Validation 7 checks consistency only *within* `threading.md`. Extend both.

---

## 2. Assumptions

**Validated:**

- The measurement harness exists as described. `.github/workflows/ci.yml:295-317` is the GNU-time
  capability probe (with the "presence is not capability" rationale in comments),
  `:318-325` the >= 8 GiB free-disk guard with the stalled-writer reasoning, `:375` the `measure`
  helper. Line numbers in the plan are accurate to a few lines.
- Plain-path memory model unchanged: `src/parallel.rs:33` `BATCH_SIZE = 4096`,
  `src/clump.rs:187-199` `channel_depths(false)` = work 2 / result_per_core 2.
- FastQC does run after trimming and report writes (`src/main.rs:1523`, `:1840-1850`; paired calls are
  sequential, one file at a time).
- `--output-format ubam` is single-threaded (`src/main.rs:452-453`).
- `--memory` is clumpify-only (`src/cli.rs:772`, `:924-925`).
- Section boundaries: `## Memory profile` is at `:50`, `## I/O and cache behaviour` at `:69`.
- No page of the docs currently mentions cgroups, `sacct`, or OOM — Behaviors 6 and 7 are genuinely new.

**Falsified:** the existing table is single-end (see C1). **Resolved for free:** the existing table is
Linux (`memory_summary.txt`).

**Unstated assumption that needs stating — the fixture.** Step 1 says "high-diversity reads, not a
replicated fixture" without naming what to use. `examples/mk_clump_fixture.rs` already exists, already
satisfies that requirement (fixed-seed SplitMix64, arbitrary read length), and already prints its
distinct-sequence count *because* a repeated fixture silently invalidated the #439 calibration once.
Name it, record the printed count — and then decide C2's compressibility question deliberately, because
its output is close to incompressible.

**Overstated: the cgroup mechanism (Behavior 6).** Assumptions hedges this correctly ("inference, not
measurement"), but Behavior 6's own wording asserts two things as fact, and one of them is
site-configuration-dependent:

1. "the kernel's page cache … is charged to the job's memory cgroup" — true for cgroup v1 and v2.
2. "SLURM's `jobacct_gather/cgroup` reports peak cgroup usage as `MaxRSS`" — true *of that plugin*, but
   which plugin a site runs is local config; under `jobacct_gather/linux`, `MaxRSS` is process RSS and
   would not include cache, while the kill still comes from the cgroup limit. Say "depending on the
   site's `jobacct_gather` plugin".

A third point is missing and matters: **clean page cache is reclaimed under pressure**, so cache
inflates *reported* usage far more readily than it causes a kill. That means `MaxRSS 1046420K` against
a 1 GiB cap is entirely consistent with cache accounting, but the exit-137 kill needs a further
ingredient — dirty writeback that cannot be reclaimed in time (which is Behavior 7's mechanism, stated
in the I/O section but not carried into the memory section), or additional heap such as the FastQC
phase. Behavior 6 should carry that caveat, otherwise the docs assert a cause the plan's own
Assumptions disclaim.

---

## 3. Validation sufficiency

Answering the measurement question directly: **no, step 1 as specified would not yield publishable
figures**, for five reasons beyond C2's fixture problem.

1. **Single run per cell.** Every existing figure in this repo is "single run per condition", and the
   clumpify calibration found within-cell spread up to 96 MiB (`src/clump.rs:263-265`, the `MARGIN`
   comment). A table a cluster budget hangs on needs >= 3 runs per cell, reporting the max.
2. **Exit status.** Validation 1 is right and the fix already exists — reuse it rather than re-derive
   it. `ci.yml`'s `measure` captures rc, checks it, then parses, with the rationale in a comment:
   *"a parseable log from a failed run looks like data"*.
3. **Nothing asserts work happened.** The harness also checks that output exists and is gzipped. A run
   that reads zero records exits 0 and emits a perfectly parseable time log with a small peak. Assert
   the output record count (or the trimming report's read count) alongside exit 0.
4. **The disk guard threshold is not transferable.** 8 GiB is calibrated for 10M × 51 bp. At 150 bp
   with "two input sizes far enough apart to show the plateau", input plus trimmed output can exceed
   that. Recompute for the chosen sizes; inheriting the constant is how a guard becomes a no-op — which
   is the exact failure Validation 2 is trying to prevent, one level up.
5. **The measurement does not exercise the condition being advised about.** `--cores 16` on a 128-core
   host with 991 GiB RAM tells you nothing about an 8-cpu, 1 GiB cgroup. Run at least one cell under a
   real cap (`systemd-run --scope -p MemoryMax=… -p CPUQuota=…`, or `docker run -m 1g --cpus 8`) and
   record both `/usr/bin/time -v` peak RSS and the cgroup's `memory.peak`. That single addition
   converts Behavior 6 from inference into measurement and retires the plan's headline Remaining Risk.

Validation 3 (plateau) and 6 (recommendation defensible both ways) are well-posed. Validation 4 is
mis-posed (C3). Validation 5 is too narrow (I1). Validation 7 should extend to `clumpy.md` (I2).
Validation 8 should name the command — `npm run build` in `docs/` (Astro 7 + Starlight, per
`docs/package.json`) — and check the wide tables render rather than overflow.

---

## 4. Edge cases

Routed correctly: `--cores 1` (worker-pool bypass), long reads (defer to #454), `--clumpify` (separate
model), uBAM (single-threaded). Three corrections and two additions:

- **`--cores 1`**: the *published* value is 5.1 MB (`threading.md:53`), so the plan's "~9 MiB paired"
  would be a third value for that one cell too, and `threading.md:35`, `:80` and `src/cli.rs:356` all
  say "5 MB". Whatever is measured, reconcile with those three sites or the page contradicts itself.
- **uBAM**: "single-threaded, so it inherits the `--cores 1` profile" is inference, not measurement —
  BAM writing carries its own `noodles` buffers, and the same notice at `src/main.rs:452-453` says
  FastQC still honours `--cores`. State that the path is single-threaded; do not attach a figure.
- **`--clumpify` cross-link**: point at `clumpy.md`'s below-floor table specifically (the 281→590 /
  334→1012 MiB floors), not the page generally. That table is the actionable one for a budget, and it
  is also the page whose `:126` needs fixing (I2).
- **Missing: `--cores` scales with the allocation.** The page's own advice is `--cores task.cpus`
  (`threading.md:86`), so memory scales with the cpus a site grants. A single absolute recommendation
  will be wrong on a 32-cpu band. Phrase it as a base plus per-core increment, anchored at `--cores 8`.
- **Missing: multi-pair / multi-sample invocations.** Pairs are processed sequentially, so peak is
  per-pair rather than summed. One clause; Snakemake users pass many pairs in one invocation.

---

## 5. Efficiency

The plan's own note ("a handful of minutes") is wrong once the matrix is real: the published beta.7
Buckberry `--cores 1` run alone was 5.5 minutes, and adding C2's reconciliation run against the 84M
fixture plus three repeats per cell puts this near an hour. Immaterial to the decision, but do not size
the session on "a handful of minutes". Measuring serially is the right call and for the right reason.

---

## 6. Alternatives worth weighing

1. **Publish the macOS figures, labelled, and ship today.** The plan bars them as platform-mixing, but
   the project already publishes cross-platform figures in a labelled block —
   `benchmarks.md:107` is an explicit Apple-Silicon section, with methodology caveats stated inline and
   raw data linked. And `src/clump.rs:250-252` already commits both platforms' numbers side by side.
   This is a live option if the Linux host is a bottleneck; the page-cache paragraph and the budgeting
   subsection are the parts that actually prevent a recurrence, and neither needs a new table.
2. **Lead with the mechanism, not the numbers.** Given C1 and C2, the most defensible change is: fix
   the five invitations (I1), add the cgroup/page-cache paragraph, add the budgeting subsection with a
   per-core formula, and cite the *existing* table rather than adding one. A new table is only worth its
   risk once the ~2× spread is explained.
3. **Put the number where it gets cited.** #11541 cited a sentence, not a section. The durable fix is a
   memory line in `## What this means for nf-core / Nextflow` next to the `--cores task.cpus` advice —
   more valuable than the choice between a `###` subsection and a top-level heading (Open question 4).
4. **Take Open question 1 (CI drift guard) more seriously than "separate follow-up".** Deferring is
   defensible, but note the base rate: this page has been materially wrong twice in one week (#454, and
   the figure #11541 cited). If the next number is going to be cited by config authors, the guard is
   what keeps it true. At minimum, commit the perf_data with host, commit SHA, fixture distinct-count
   and run count, so the next correction is cheap.

On Open questions generally: the four asked are all real, and their assumed answers are reasonable. But
two Critical items are absent or mis-filed — the single-end attribution is asserted as settled fact in
the Goal when it is wrong (C1), and the reconciliation of the existing 2× spread, including the
possibility that it is a code regression, is not raised at all (C2). Open question 2 (cite the incident)
I agree with as assumed: describe the mechanism, do not name the incident.

---

## 7. Action items

### Critical

1. **Correct the premise.** The existing table is paired-end, Linux, real WGBS data
   (`docs/perf_data/buckberry-2026-04-29/memory_summary.txt`; `--paired` in `scripts/benchmark.sh:130`;
   wall times match the hyperfine means within ~1%). Rewrite the Goal's cause chain, drop the "~2×,
   so ~5× headroom" arithmetic, and re-derive Behavior 2, step 3, Assumption 2 and Open question 3.
2. **Add a reconciliation step and a stop rule.** Re-measure the published condition (Buckberry pair,
   `--paired --cores 8`) with today's binary before publishing anything. If it does not reproduce
   ~92 MB, stop: that is a memory regression and a code plan, not a docs plan.
3. **Explain the ~2× spread across the four existing figures** (C2 table) before adding a fifth, and
   name the fixture explicitly. If random `mk_clump_fixture` data is used, label the figures an upper
   bound for real samples; if the recommendation is meant for real data, measure real data.
4. **Rewrite Behavior 3.** `src/clump.rs:255-257` says FastQC "adds rather than maxes" and charges up
   to 384 MiB; `clumpy.md:114` publishes that additivity; `-t` drives fastqc-rust's rayon pool
   (`src/fastqc.rs:35`). Do not publish "not additive". Re-pose Validation 4 as a controlled
   low-`t` vs high-`t` comparison at identical input, and size the cluster recommendation with FastQC.
5. **Add the decision rule for step 1's outcome** — what the deliverable becomes if the Linux figures
   land near the published table.

### Important

6. Extend step 4 to all five sites that invite the underestimate: `threading.md:61`, `:63` (both
   clauses), `:80`, and the `## What this means for nf-core / Nextflow` section at `:84-88`. Widen
   Validation 5's grep accordingly.
7. Fix `clumpy.md:126` — it still asserts "depends on read length rather than on file size", falsified
   by `0e86e2a` this morning, and re-quotes ~100 MB. Extend Validation 7 to cover it.
8. Strengthen step 1's methodology: >= 3 runs per cell reporting max; reuse `ci.yml`'s `measure`
   rc-then-parse discipline verbatim; assert output record count as well as exit 0; recompute the
   free-disk threshold for 150 bp at the chosen sizes.
9. Add one cell measured under an actual cgroup cap, recording both `time -v` peak RSS and
   `memory.peak`. This is what turns Behavior 6 from inference into measurement.
10. Soften Behavior 6's two asserted mechanisms: `MaxRSS`-includes-cache is `jobacct_gather` plugin
    dependent, and clean cache is reclaimed under pressure — so cache explains an inflated *report*
    more readily than a kill. Carry Behavior 7's dirty-writeback mechanism into the memory section.
11. Express the recommendation as a base plus per-core increment anchored at `--cores 8`, not one
    absolute number, because the page tells users to pass `--cores task.cpus`.
12. Note that #454's paragraph has no perf_data audit trail (`0e86e2a` touched only `CHANGELOG.md` and
    `threading.md`) and no platform label. If the new material is platform-labelled, that paragraph
    needs the same label — an addition, not a contradiction.
13. Add the multi-pair clause (peak is per pair, not summed) and drop the uBAM figure, keeping only
    "single-threaded".

### Optional

14. Unit hygiene: the section mixes MB (existing table) and MiB (#454 paragraph); GNU time reports
    KiB. State the conversion once, use MiB in all new material.
15. perf_data links in the docs are `tree/master/…` (`benchmarks.md:107,145`); a directory added on
    `dev` will 404 until the fast-forward. Use a commit permalink or the `dev` path.
16. Name Validation 8's command (`npm run build` in `docs/`) and check table rendering, not just build
    success.
17. Correct the Efficiency estimate (~an hour with repeats and the reconciliation run, not "a handful
    of minutes").
18. Consider Alternative 1 (labelled macOS figures, shipping now) if the Linux host is the bottleneck —
    it matches `benchmarks.md`'s existing precedent and `src/clump.rs`'s two-platform comment.
