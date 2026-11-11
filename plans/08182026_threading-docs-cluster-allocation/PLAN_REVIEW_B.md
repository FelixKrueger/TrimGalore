# Plan Review B — paired-end memory figures and cluster-allocation guidance

**Plan:** `plans/08182026_threading-docs-cluster-allocation/PLAN.md`
**Reviewer:** B (independent)
**Date:** 2026-08-18
**Verdict:** the *goal* is right and worth doing; the *premise it is argued from is factually wrong*,
and two of the eight validations cannot fail. Do not implement as written.

---

## 0. What I verified, and what checked out

Everything the plan asserts about file locations and the CI harness is accurate:

| Claim in plan | Verified |
|---|---|
| `## Memory profile` at ~line 50, `## I/O and cache behaviour` at ~line 71 | ✓ lines 50 and 69 |
| `clumpify-memory` job from ~line 271 | ✓ `.github/workflows/ci.yml:271` |
| GNU-time capability probe ~302-316; free-disk guard ~319-325; `measure` ~380 | ✓ 300-318, 319-325, 370-434 |
| FastQC runs *after* trimming and report writes | ✓ `src/main.rs:1523-1531` (SE), `:1840-1861` (PE); R1 then R2 sequentially |
| `-t N` in `--fastqc_args` is honoured, not ignored | ✓ `src/fastqc.rs:101-106` sets `config.threads` |
| Plain-path memory geometry unchanged since the published table | ✓ `git show v2.1.0-beta.7:src/parallel.rs` — `BATCH_SIZE 4096`, work depth 2, result `cores*2`; same at HEAD via `clump::channel_depths(false)` |
| uBAM output is synchronous | ✓ `src/main.rs:2218` uses `open_sync_reader` |

Two things I could **not** verify and that the plan asserts anyway are in §1.1 and §1.2 below. Both
change what gets published.

---

## 1. Logic review

### 1.1 CRITICAL — the existing table is paired-end. The plan's central premise is false.

The plan's Goal says the cited ~100 MB figure is *"our figure, measured single-end on short reads"*
and that *"paired-end at `--cores 8` is roughly twice that, so the real headroom was ~5×"*. The audit
trail says otherwise.

`docs/perf_data/buckberry-2026-04-29/memory_summary.txt`:

```
Fixture: SRR24827373 (Buckberry 2023, 84M PE reads)
| --cores | Peak RSS (MB) | Wall (s) |
| 1 | 5.1 | 5:27.36 |   | 2 | 60.1 | 2:39.38 |   | 4 | 73.0 | 1:20.75 |   | 8 | 91.6 | 0:56.93 |
```

Those wall times are the tell. The paired hyperfine means in the same directory are 328.85 s (c1),
80.88 s (c4), 56.61 s (c8) — against 327.36 / 80.75 / 56.93 in the memory run, i.e. within ~0.5%. And
every Rust sweep in `scripts/benchmark.sh:130-137` is `--paired`. A single-end run over one of the two
files would be roughly half those wall times. **The published 91.6 MB at `--cores 8` is a paired-end
figure on 84M pairs, measured on Linux (`dockyard-oxy-0`).**

I also pulled the consumer PR to check the wording. nf-core/modules#11541 says: *"Memory drops to 1 GB,
~10× the observed ~100 MB peak on 30M PE."* That was an **accurate** reading of a paired-end number
from our docs. The mis-citation the plan is built on did not happen.

Consequences, all of which need action:

1. The Goal's causal chain has to be rewritten. The docs defect is *not* SE-vs-PE. It is (a) a
   `getrusage` heap figure being read as a scheduler budget, and (b) the "noise floor of cluster
   scheduler memory allocation" framing that made 1 GB look generous. Behaviors 5, 6 and 7 survive
   intact; Behavior 2's *rationale* does not.
2. Behavior 2 as written would **introduce** an error: it instructs the implementer to place the new
   table "beside the existing single-end table". Labelling a paired table single-end publishes a fresh
   false claim on the page whose job is to stop false claims. (This is the "corrective copy contains
   the defect it corrects" pattern, and the fifth or sixth instance of it in this repo.)
3. The same false claim is already public. nf-core/rnaseq#1905's `<details>` block says *"That figure
   is from Trim Galore's own docs and was measured single-end on short reads; paired-end at `--cores 8`
   is roughly twice it."* #1905 is **still open** (`mergedAt: null` as of today), so the body is
   editable — correct it there rather than leaving our own repo cited for something it doesn't say.
   The plan's Goal also states the fix "shipped as nf-core/rnaseq#1905"; it has not merged yet.
4. The real puzzle the measurement must answer changes shape. Every *other* paired plain-path figure
   in this repo is 40-60% higher than the published table at the same core count: `docs/.../clumpy.md:126`
   gives **104 MiB at `--cores 4`** (10M × 51 bp paired), commit `2ab0918`'s table gives **103 MB at
   `--cores 4`** (4M pairs, macOS/arm64), the plan's own macOS diagnostic gives **116 MiB at `--cores 4`** —
   against threading.md's **73.0 MB**. The interesting question is why the published table is the low
   outlier, not why paired differs from single-end.

### 1.2 CRITICAL — version and allocator confound: the plan would publish figures for a build nobody runs

`mimalloc` became the global allocator in `2ab0918` (#437, 2026-08-16) — **after** the `v2.3.0` tag
(`git show v2.3.0:Cargo.toml | grep -c mimalloc` → 0). The published table is `v2.1.0-beta.7`, i.e.
the platform allocator. So a Linux run on current `master` would put two allocators and three
releases' worth of code in one docs section, and the numbers would describe a binary that no user can
install — the nf-core module pins the `trim-galore:2.3.0` container
(`modules/nf-core/trimgalore/main.nf`), which is stock-allocator.

This is not hypothetical for memory: #437's own commit message reports the allocator swap moving
`--clumpify` peak RSS by -24%, and mimalloc's page-retention policy is precisely what decides whether
RSS falls between the trim phase and the FastQC phase — the mechanism Behavior 3 depends on.

Pick one and say so in the table lead: measure the **released** `v2.3.0` binary (what readers can run,
and what #1903 hit), or measure `master` and label version *and* allocator on both tables with a note
that the figures move at the next release.

### 1.3 CRITICAL — Behavior 3 contradicts the shipped memory model, on a page the plan links to

The plan wants to publish *"FastQC threads are not additive in memory"*. `src/clump.rs:253-256` says
the opposite, in the source:

```rust
/// FastQC's additional resident cost per thread. Its phase follows trimming but
/// the pool's pages stay resident, so it adds rather than maxes.
const STATIC_FASTQC_PER_THREAD_BYTES: u64 = 24 * 1024 * 1024;
```

charged linearly to `FASTQC_CHARGED_THREAD_CAP = 16` (up to 384 MiB), pinned by
`fastqc_reservation_scales_with_threads_and_is_capped`, gated in CI by the `1G_c4_fastqc` cell
(`ci.yml:440`), and **published as a formula** at `docs/.../clumpy.md:96` and explained at `:102`
("the FastQC term is charged only when a report is requested, capped at 16 threads"). Step 8 of this
plan cross-links to that very page. A reader who follows the link finds a formula charging 24 MiB per
FastQC thread next to a sentence saying FastQC threads cost nothing.

The plan's structural argument is also only half the story, and the stronger half is missing:
`fastqc-rust`'s `-t` parallelises **across file groups** (`runner.rs:71-100`: one rayon pool,
`file_groups.par_iter()`, and `build_file_groups` makes each file its own group), while Trim Galore
hands `fastqc::run` exactly **one** file per invocation (`src/main.rs:1525`, `:1842`, `:1848`). So
`-t 24` has nothing to parallelise *by construction of how we call it* — a much better reason than
phase ordering, and one that also explains why the clumpify model's per-thread charge is a
conservative reservation rather than a measurement of additivity (its gate is `peak <= prediction`,
which passes whenever the charge is too generous).

Publishable reconciliation, in decreasing order of safety:

- State the mechanism, not a rule: FastQC's phase follows trimming, and TG hands it one file per
  invocation, so `-t` beyond 1 has no additional file to work on. Report the **measured whole-run peak
  with `--fastqc`** as a figure, and let the number carry the claim.
- Say explicitly that `--clumpify`'s budget still charges per FastQC thread, and why (the bin pool's
  pages stay resident, and the budget is a floor that must not under-promise) — otherwise the two
  pages read as a contradiction.
- Do **not** publish a general additivity rule across input types and flag combinations. Nothing in
  the plan's matrix tests uBAM input, `--nogroup` (which widens FastQC's per-base arrays), long reads,
  or multi-file invocations.

### 1.4 CRITICAL — Behaviors 4 and 6 contradict each other

Behavior 4: the peak plateaus with input size, *"so a bigger sample does not need a bigger allocation
at a given read length"*. Behavior 6: the accounted quantity includes page cache for the input read
and the `.fq.gz` write. Page cache grows with file size. Both cannot be true of the same number, and
the reader will apply whichever they read last to their config.

The fix is a units discipline the plan currently lacks: say which quantity each statement is about
(`ru_maxrss` heap vs cgroup `memory.peak`), and state the plateau as a property of the heap only.

### 1.5 Step 4 under-scopes the rewrite (Important)

`grep -n` on the target file finds four sites carrying the framing, not one:

- `:35` — *"each added core is exactly 1 worker thread plus ~10 MB of memory"*
- `:61` — *"Memory growth is bounded and predictable, **well-suited to cluster scheduling**: … the process stays under ~100 MB"*
- `:63` — *"…well below the noise floor of cluster scheduler memory allocation"*  ← the only one step 4 addresses
- `:80` — *"**Lower memory:** 5 MB single-threaded, ~10 MB per additional worker"*

#11541's wording ("~10× the observed ~100 MB peak") is closer to `:61` than to `:63`. And `:35`/`:80`
become *wrong* the moment a paired table lands: the plan's own diagnostics imply ~13-36 MB per added
core (9→116→202→304 MiB across c1/c4/c8/c16), not ~10 MB. Validation 5's grep must cover all four.

### 1.6 clumpy.md already publishes a conflicting plain-path figure (Important)

`docs/src/content/docs/performance/clumpy.md:126`:

> Memory usage without `--clumpify` **depends on read length rather than on file size**. Short-read
> runs sit around the 100 MB mark (**104 MiB for 10 M × 51 bp paired at `--cores 4`**) …

Two problems. The figure disagrees with threading.md's 73.0 MB at the same core count — and disagrees
in the awkward direction, with *shorter* reads and a *smaller* input using *more* memory. And the first
clause contradicts the #454 paragraph the plan is told not to contradict (`threading.md:67`: "At a
fixed read length it also grows with input size, until the pipeline saturates"). The plan edits this
file (step 8, the cross-link) but never reconciles it, and Validation 7 only re-reads threading.md.

---

## 2. Validation sufficiency

| # | Verdict |
|---|---|
| 1 | Good, and the reasoning is right (a failed run still emits a parseable log). Keep. |
| 2 | Good — a positive control on the guard. Note it now tests a *local* script, not the CI job (see §3). |
| 3 | Real. Fine. |
| 4 | **Cannot fail.** See §1.3: TG passes one file per `fastqc::run`, and `-t` parallelises across files, so "no difference" is structurally guaranteed regardless of whether the claim generalises. Same shape as the mutation-direction false negative in this repo's history: an experiment whose only possible outcome is the one asserted. |
| 5 | Under-scoped — one of four sites (§1.5). |
| 6 | **Circular.** Both legs are heap-side ("comfortably above measured peak", "not so large"), neither has a threshold, and the quantity being recommended (cgroup-accounted memory) is never measured. Any number in a very wide band passes. |
| 7 | Necessary but too narrow: docs-vs-docs within one file. Extend to `clumpy.md` (§1.6) and to docs-vs-source (§1.3). |
| 8 | Insufficient as written — see §2.2. |

### 2.1 CRITICAL — make Validation 6 a real measurement

The plan's own Assumptions concede the page-cache mechanism is *"inference, not measurement — no
cgroup reproduction was run"*, and then Validation 6 checks the recommendation against heap figures.
So the plan derives a cgroup budget from a quantity that excludes cache, using a mechanism it hasn't
observed. That is exactly the gap that produced the incident, one level up.

It is cheap to close. On `dockyard-oxy` the process already runs inside a cgroup, so after a run:

```bash
cat /sys/fs/cgroup/memory.peak          # cgroup v2: the quantity SLURM's jobacct_gather/cgroup reports
# or, for an isolated cell:
systemd-run --scope -p MemoryMax=1G -p MemoryAccounting=yes ./target/release/trim_galore …
```

Two results worth having, both falsifiable:

- **Both numbers side by side** for the same run — `/usr/bin/time -v` peak RSS vs cgroup `memory.peak`.
  That single pair of numbers is the entire argument of the new subsection, and it makes Behavior 6 a
  measurement instead of an inference.
- **The 1 GiB reproduction**: does a 30M-pair paired run at `--cores 8` with `--fastqc_args '-t 24'`
  die under a 1 GiB cap? If yes, the recommendation is anchored in an observed failure. If no, the
  page-cache story is not sufficient to explain #1903 and the docs must not lean on it — an equally
  valuable answer, and one that would stop us publishing a mechanism that isn't the one that bites.

Note the mechanics while writing Behavior 6, or the page will overshoot: clean page cache is
*reclaimed* under cgroup pressure, so cache inflates what the scheduler **reports** far more readily
than it causes a **kill**; a kill needs unreclaimable memory (anon plus dirty pages awaiting
writeback). Behavior 7 has the dirty-page half — Behavior 6 should not imply that reported MaxRSS and
kill risk are the same quantity, or authors will size `heap + file size` and reintroduce the
over-allocation #11541 was reducing.

### 2.2 IMPORTANT — Validation 8 must run the docs-build *assertions*, not `npm run build`

`.github/workflows/ci.yml:160-225` gates every docs change, and its two "silent render corruption"
steps exist because the build **exits 0** on both failure modes:

- a bare `~` opening GFM strikethrough → `<del>` in the built HTML (*"Two of these shipped to the live
  site; all three keep the build green"*);
- a bare `$` in prose parsed as maths → stray `class="katex"` on any page but `performance/clumpy`.

This plan adds a paragraph dense with `~NNN MiB` figures and a budgeting subsection likely to contain
currency- or dollar-adjacent prose. Run the two assertion steps, not just the build. Also: the plan's
Integration section says "no CI files touched … the `validation` CI job … untouched" — true, and the
byte-identity reasoning in §4 below holds, but "no CI exposure" is wrong in kind: `docs-build` gates
this change.

### 2.3 CRITICAL — one harness hazard the plan will hit, in the safe-looking direction

Implementation step 1 specifies *"Plain paired mode, 150 bp, `--gzip`"*. `--gzip` does not force gzip:
output compression follows `input[0]` (`src/main.rs:208`, `:438` → `naming::output_is_gzipped`), and
`--gzip` warns and is ignored (`src/cli.rs:1187`, the #453 fix). A plain-text fixture plus `--gzip`
therefore produces **plain-text output, no deflate state, an understated peak, and exit 0** — with only
a warning buried in the log. Every figure in the new table would be quietly low, and low is the
direction that gets someone OOM-killed.

Carry over the clumpify job's positive-engagement assertions, which exist for exactly this class of
failure:

- `ls "$out"/*.fq.gz` (`ci.yml:382`) — output is actually gzipped;
- the plain-path analogue of `grep -q '^clumpify:'` — assert the startup banner says
  `Using N worker threads` for the N you passed;
- record `trim_galore --version` output into the perf_data file per run, given §1.2.

### 2.4 IMPORTANT — single-run peaks are not a basis for a budget

Behavior 1 says *"one Linux measurement run"* per condition, matching the existing table ("single run
per condition"). But `src/clump.rs:263-265` records *"Largest observed within-cell spread is 96 MiB"*
as the reason the clumpify budget holds back 9%. Even allowing that that spread is a clumpify-path
observation, a project that has measured ~96 MiB of within-cell variation should not derive a cluster
recommendation from n=1. Run ≥3 per cell; publish the max (a budget is an upper bound, not a mean) and
state the spread.

---

## 3. Harness reuse: which clumpify idioms transfer, and which leak

Three of the four transfer cleanly, and the plan is right to reuse rather than hand-roll:

- **GNU-time capability probe** — transfers, and is worth keeping even on a known host. "Presence is
  not capability" is host-agnostic.
- **Capture-then-check-status-then-parse** in `measure` (`ci.yml:373-386`) — transfers, and is the
  single most important idiom for this task (it is Validation 1).
- **Positive-engagement assertions** — transfer *in shape but not in content*; the plain path needs its
  own (see §2.3). This is the one place where copying literally would leave the harness proving nothing.
- **The `>= 8 GiB` free-disk guard** — the *reasoning* transfers ("disk is not independent of memory: a
  stalled writer backs up the result channel … so an out-of-space run reports inflated peak RSS"), the
  *threshold does not*. 8 GiB was sized for two 10M-pair 51/25 bp fixtures. A 150 bp paired matrix at
  two input sizes, one of them large enough to show the plateau, plus per-run outputs, is several times
  that. Recompute the number from the actual fixture plan; a guard that passes because the threshold is
  too low is the same no-op Validation 2 exists to catch.

Two assumptions of that job that do **not** apply, and should not be carried in by reflex: the drift
band and `peak <= prediction` gate are meaningful only because `--clumpify` *prints* a prediction; the
plain path has no `--memory` model and no prediction, so there is nothing to ratio against. The plan
is right to defer a plain-path CI guard to a follow-up — but when that follow-up is written, the
baseline has to be an absolute band, not a `peak/pred` ratio.

---

## 4. Scope boundary: does "docs-only" hold?

Mostly yes, with one correction and one caveat.

- **Byte-identity / `validation` matrix: correct.** No source change, so no exposure. Independently:
  `2ab0918` re-verified `_val_1`/`_val_2` byte-identity across an allocator swap, so even the
  master-vs-2.3.0 question in §1.2 is a *figures* problem, not a byte-identity one.
- **Correction:** `docs-build` gates the change (§2.2). "No CI exposure" is not accurate.
- **Caveat:** the plan writes a new `docs/perf_data/` subdirectory. `Cargo.toml:14-22` excludes
  `docs/` from the published crate, so no crates.io package-size effect — but keep the raw data small
  (the existing directory bundles logs into a tarball for exactly this reason).
- **One scope item the plan should adopt, not defer:** `clumpy.md:126` (§1.6). It is one sentence and
  one figure in a file the plan already opens, and leaving it means the site publishes two
  contradictory plain-path paired numbers and two contradictory statements about input-size scaling.

---

## 5. How the published result could still mislead a workflow author

Even executed exactly as written:

1. **A number without its axes.** A single "Budgeting on a cluster" figure will be copied into configs
   regardless of `--cores`, read length, `--fastqc`, `--clumpify`, or output compression. Give it as a
   short row-wise table or a formula keyed on `--cores`, with the read-length caveat *in the same
   sentence* as the number.
2. **`--cores` is not `task.cpus`.** The nf-core module computes `cores = task.cpus - 4` (PE, capped at
   8, floored at 1) — so at `process_medium`'s 6 cpus the tool runs at `--cores 2`, and #1903's own log
   shows `--cores 8` only because rnaseq overrides cpus on that selector. Meanwhile `threading.md:85-89`
   tells authors *"If a Nextflow process has 12 CPUs, just pass `--cores 12`"*. A memory recommendation
   keyed on cpus is wrong under the module and wrong under our own advice at the same time. Key it on
   `--cores`.
3. **Wrong home.** `## What this means for nf-core / Nextflow` (`:85`) is the section a config author
   reads, and it currently gives cpu advice with **no** memory advice — that pairing is what invites a
   low memory label. Put the recommendation there, or put a one-line pointer there and make the
   subsection a linkable heading (the plan's fourth open question, which I'd resolve in favour of
   linkable-and-cross-referenced).
4. **Disagreeing with the config that shipped.** rnaseq#1905 proposes 4 GB. If our docs land on 1-2 GB,
   someone will "correct" the config down to match us and reproduce #1903 — and they will be able to
   cite us for it. Either land on a number consistent with 4 GB or state explicitly why the docs number
   is smaller (e.g. it is a floor for the tool, not a retry-scaled allocation).
5. **The plateau, read as a cgroup claim.** §1.4. "A bigger sample needs no bigger allocation" is true
   of the heap and false of what `sacct` prints.
6. **Provenance the audit trail cannot support.** `memory_summary.txt` records neither the read length
   nor whether "84M" is pairs or reads, yet threading.md publishes both ("84M-read", "150 bp"). The
   plan's confirm-before-publishing list has platform (confirmed: `dockyard-oxy-0`, Linux) but not
   these. State *pairs* explicitly in the new table's lead, and either confirm the 150 bp attribution
   or drop it from the old table's caption.
7. **uBAM.** `open_sync_reader` makes it genuinely single-threaded, but the path adds noodles record and
   bgzf writer state, so it does not literally inherit the `--cores 1` row. Say "closer to the
   `--cores 1` row than to the multi-core rows" without attaching a number, or measure it.

---

## 6. Efficiency

Nothing to add on the published change. On the measurement step the plan is right that serial
measurement is mandatory (concurrent runs contend for page cache — which, given Behavior 6, would
corrupt the one quantity that matters most). Two additions: state the cache policy per run (the
existing `benchmark.sh` warms the cache once up front; a cold-cache run and a warm-cache run will not
report the same cgroup peak, and the difference is the point), and budget for §2.4's repeats.

---

## 7. Alternatives worth considering

1. **Publish both quantities, not one.** The durable fix for a downstream consumer is a table with two
   columns — heap (`/usr/bin/time -v`) and cgroup peak (`memory.peak`) — because the mis-citation
   happened when a reader took a heap number for a scheduler number. Two columns make that mistake
   structurally impossible. This is a stronger deliverable than a recommendation paragraph, and §2.1's
   measurements produce it for free.
2. **A formula rather than a number.** `static + per-core × cores`, with the cache term named
   separately, generalises to core counts we did not measure and does not rot when the allocator
   changes. It also matches how `clumpy.md:96` already presents the clumpify budget, so the two pages
   would read as one system.
3. **Print it at runtime (follow-up, source change).** `--clumpify` and `--clump_only` both print a
   predicted peak and the budget it was sized against (#456). A plain-path banner would put the number
   in front of every user without anyone having to cite a docs page — the only fix that survives the
   next relabelling PR. Out of scope here; worth filing.
4. **Do the CI drift guard now.** The plan defers it and gives a good reason. I agree with deferring,
   with one caveat: a published figure with no guard is exactly how the current table became a four-month-old
   `beta.7` number that every other measurement in the repo disagrees with (§1.1 item 4).

---

## 8. Action items

### Critical

1. **Correct the premise.** The existing table is paired-end, Linux, 84M pairs
   (`docs/perf_data/buckberry-2026-04-29/memory_summary.txt`, wall times matching the `--paired`
   hyperfine runs, `scripts/benchmark.sh:130-137`). Rewrite the Goal, drop "single-end" from Behavior 2,
   and re-aim the plan at heap-vs-cgroup accounting and the `:61`/`:63` framing. **Do not label the
   existing table single-end.**
2. **Correct nf-core/rnaseq#1905's `<details>` block**, which carries the same false claim and is still
   open. Also fix the plan's "the fix shipped as #1905" — it has not merged.
3. **Resolve the version/allocator confound** (§1.2): measure `v2.3.0` (what users run), or measure
   `master` and label version + allocator on both tables. mimalloc post-dates `v2.3.0` (`2ab0918`, #437).
4. **Reframe Behavior 3** so it does not contradict `src/clump.rs:253-256` and `clumpy.md:96`. Publish
   the structural reason (one file per `fastqc::run`; `-t` parallelises across files) plus the measured
   with-`--fastqc` peak. Drop the general "not additive" rule. **Replace Validation 4**, which cannot fail.
5. **Make Validation 6 a measurement** (§2.1): cgroup `memory.peak` beside `/usr/bin/time -v` peak RSS
   for the same runs, plus the 1 GiB-cap reproduction. Without it the recommendation is a cgroup number
   inferred from heap data, checked against heap data.
6. **Reconcile Behaviors 4 and 6** by naming the quantity in every claim (§1.4).
7. **Harden the measurement against the `--gzip` no-op** (§2.3): gzipped input fixtures, assert `.gz`
   output and the worker-thread banner, record `--version` per run.

### Important

8. Extend step 4 and Validation 5 to all four framing sites: `threading.md:35`, `:61`, `:63`, `:80` (§1.5).
9. Reconcile `clumpy.md:126` — its 104 MiB at `--cores 4` and its "depends on read length rather than on
   file size" both conflict with threading.md (§1.6). Extend Validation 7 to that file and to docs-vs-source.
10. Validation 8 must run `docs-build`'s two render-corruption assertions, not just `npm run build`
    (§2.2). Correct the Integration section's "no CI exposure".
11. Key the recommendation on `--cores`, not cpus; put it where the cpu advice lives (`:85`) or
    cross-link from there; and reconcile the number with rnaseq#1905's 4 GB (§5 items 2-4).
12. Split Behavior 6 into "what the scheduler reports" and "what gets you killed" — reclaimable cache
    vs anon + dirty (§2.1, last paragraph).
13. ≥3 runs per cell, publish the max plus the spread (§2.4).
14. Recompute the free-disk threshold for this fixture matrix rather than copying 8 GiB (§3).
15. Add read length and pairs-vs-reads to the confirm-before-publishing list; platform is already
    confirmed as Linux/`dockyard-oxy-0` (§5 item 6).

### Optional

16. The CHANGELOG entry needs an issue/PR number — every `Unreleased` entry has one.
17. Soften the uBAM clause or measure it; `open_sync_reader` is single-threaded but the bgzf writer is
    not free (§5 item 7).
18. Consider publishing both heap and cgroup columns, and/or a formula (§7 items 1-2).
19. File the plain-path "predicted peak" banner as a follow-up (§7 item 3) — the only fix that does not
    depend on a downstream author reading our docs correctly.

---

## 9. What I would keep unchanged

The plan's instincts are good where it matters: refusing to publish the macOS figures beside Linux
ones; splitting Validation 1 into value-and-exit-status; adding the positive control on the disk guard;
routing `--cores 1`, long reads, `--clumpify` and uBAM in Behavior rather than leaving them to the
implementer; deferring the CI guard; and flagging the page-cache mechanism as inference in Assumptions
rather than asserting it. The self-review's three "adjusted during review" items are all real
improvements. The problem is not carelessness — it is that the plan inherited a wrong fact from
#1905's analysis and never checked it against `docs/perf_data/`, which is the one place that could
settle it.
