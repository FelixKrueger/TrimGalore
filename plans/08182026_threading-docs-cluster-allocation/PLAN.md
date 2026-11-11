# Plan: paired-end memory figures and cluster-allocation guidance in the threading docs

**Date:** 2026-08-18
**Type:** Documentation (no source, CLI, or behaviour change)

## Goal

Extend the **Memory profile** section of `docs/src/content/docs/performance/threading.md` so it carries
paired-end peak-RSS figures and an explicit cluster memory-allocation recommendation, and so it states
that cgroup-based schedulers charge page cache to the job.

Those three gaps are load-bearing rather than cosmetic. nf-core/modules#11541 sized `TRIMGALORE` at
1 GB on the stated rationale *"~10× the observed ~100 MB peak on 30M PE"* — our figure, measured
single-end on short reads. Paired-end at `--cores 8` is roughly twice that, so the real headroom was
~5×, against a number that counts only the heap. Every sample in nf-core/rnaseq#1903 was then
OOM-killed at exit 137. The fix shipped as nf-core/rnaseq#1905; this plan closes the documentation
gap that invited it.

## Context

**Target file:** `docs/src/content/docs/performance/threading.md`

- `## Memory profile` (currently ~line 50 through the paragraph ending `…than the short-read figures
  suggest.`) — the section to extend.
- `## I/O and cache behaviour` (currently ~line 71) — one caveat to add; its claim *"disk is not the
  bottleneck for any realistic deployment"* is true about bandwidth but reads as absolute, and slow
  writeback under a hard memory cap is exactly the mechanism implicated in rnaseq#1903.

**What already exists and must not be redone:**

- The single-end table (84M-read Buckberry fixture, v2.1.0-beta.7, `/usr/bin/time -v`): `--cores`
  1 / 2 / 4 / 8 at 5.1 / 60.1 / 73.0 / 91.6 MB, with the N+4 thread model explained beneath it.
- Commit `0e86e2a` (#454 / #465, 2026-08-18) already re-qualified the long-read ceiling: read-length
  progression at `--cores 4` (~75 MiB at 150 bp, ~245 MiB at 1 kb, ~800 MiB at 10 kb) and the
  input-size growth to a ~2.3 GiB plateau at 10 kb. **This plan builds on that paragraph and must not
  restate or contradict it.**

**Measurement harness to reuse — do not hand-roll one:** the `clumpify-memory` job in
`.github/workflows/ci.yml` (from ~line 271) already solves the hard parts on Linux:

- a capability probe for a `/usr/bin/time` that actually reports peak RSS, installing GNU time when
  the preinstalled one cannot (~line 302-316) — *"presence is not capability"*;
- a `>= 8 GiB` free-disk guard (~line 319-325), on the explicit reasoning that *"disk is not
  independent of memory: a stalled writer backs up the"* pipeline and inflates peak RSS;
- a `measure` helper (~line 380) wrapping `/usr/bin/time -v ./target/release/trim_galore …`.

**Related but separate:** `--clumpify` has its own predictive `--memory` budget model
(`src/clump.rs`, `channel_depths`, `LayoutInputs`) with its own CI baseline in
`.github/clumpify-memory-baseline.tsv`. The plain path has no `--memory` accounting. The docs must
point at the clumpify model rather than blur the two.

## Behavior

What the section must say once changed:

1. **Prerequisite.** Every published figure comes from one Linux measurement run on the plain
   (non-clumpify) paired path, produced with the ci.yml harness idioms, with each run's exit status
   asserted zero.
2. **Paired-end figures.** A paired-end peak-RSS table at 150 bp for `--cores` 1 / 2 / 4 / 8 / 16,
   stated as paired-end and platform-labelled, sitting beside the existing single-end table rather
   than replacing it.
3. **FastQC is not additive in memory.** `--fastqc` / `--fastqc_args -t N` runs *after* trimming
   completes (`src/main.rs` — the `fastqc::run` calls follow the trimming and report writes), so a
   large `-t` does not add to the trimming peak. State this, because `-t 24` on an 8-cpu allocation
   looks alarming in a log and is not the cause of a kill.
4. **Plateau at fixed read length.** Note that the paired figure stops growing with input size once
   the pipeline saturates, consistent with the existing #454 paragraph, so a bigger sample does not
   need a bigger allocation at a given read length.
5. **An explicit cluster recommendation.** A short "Budgeting on a cluster" subsection giving a
   number a workflow author can put in a config, with the reasoning for why it exceeds measured peak.
6. **Page cache is charged to the job.** Under cgroups, the kernel's page cache for the FASTQ read and
   the `.fq.gz` write is charged to the job's memory cgroup, and SLURM's `jobacct_gather/cgroup`
   reports peak cgroup usage as `MaxRSS`. So the number a scheduler shows is not the tool's heap, and a
   budget sized only from heap figures is unsafe.
7. **Writeback caveat in the I/O section.** Clean page cache is reclaimable; dirty pages are not until
   writeback flushes them. On slow shared cluster storage under a hard memory cap, the dirty pool can
   collide with the cap. Bandwidth is still not the bottleneck — the existing claim stands, scoped.

**Edge cases to keep correct:**

- `--cores 1` bypasses the worker pool entirely, so it is the one setting with no batch residency
  (~9 MiB paired). Do not let a "budget N cores × M MB" rule imply otherwise.
- Long reads: defer to the existing #454 paragraph. The cluster recommendation must be stated as
  short-read and must point at that paragraph for ONT/PacBio.
- `--clumpify`: cross-link to its `--memory` budget; do not imply the plain-path recommendation covers it.
- uBAM output (`--output-format ubam`) is always single-threaded, so it inherits the `--cores 1`
  profile. Worth one clause if it fits without bloating the section.

## Implementation outline

1. **Measure on Linux.** On a Linux host (dockyard-oxy per the usual `dcli` workflow; source
   `isoquant_env` first, and `~/.cargo/env` for the toolchain). Plain paired mode, 150 bp, `--gzip`:
   - `--cores` 1, 2, 4, 8, 16, at two input sizes far enough apart to show the plateau;
   - one `--cores 8` run adding `--fastqc_args '-t 24'`;
   - use high-diversity reads, not a replicated fixture — a fixture that repeats the same reads
     understates any per-sequence state in FastQC's duplicate tracking;
   - carry over the harness's GNU-time capability probe and free-disk guard, and assert exit 0 on
     every run before recording a figure.
2. **Record the raw data** under `docs/perf_data/threading-paired-2026-08-<dd>/`, matching the
   existing `docs/perf_data/` convention, so the published table has an audit trail like #287's.
3. **Add the paired-end table** immediately after the single-end table, with a one-line lead saying
   which fixture, which version, and which platform.
4. **Rewrite the sentence that invited the 1 GB cap.** Currently: *"…well below the noise floor of
   cluster scheduler memory allocation"*. That framing is what makes a 1 GB budget look generous.
   Replace with a statement that the heap is small but is not the quantity a cgroup scheduler
   accounts, and point at the new subsection.
5. **Add `### Budgeting on a cluster`** with the recommendation, the reasoning, the `--cores 1`
   exception, and the pointer to the long-read paragraph for ONT/PacBio.
6. **Add the page-cache paragraph** (Behavior 6), including that a `MaxRSS` far above the figures in
   these tables points at cache accounting rather than at the tool, and naming the
   `/usr/bin/time -v` vs `sacct` comparison as the way to tell them apart.
7. **Scope the I/O claim** (Behavior 7) with one sentence in `## I/O and cache behaviour`. Keep the
   bandwidth finding intact.
8. **Cross-link `--clumpify`'s `--memory`** so the two memory models are not conflated.
9. **CHANGELOG entry** — one sentence, following the existing docs-entry style in `CHANGELOG.md`.

## Efficiency

Not applicable to the published change (documentation only). One note on the measurement step: the
`--cores 16` and large-input runs dominate wall-clock; the whole matrix is a handful of minutes on a
release binary and needs no parallelisation. Measure serially — concurrent runs on one host contend
for page cache and would corrupt the very figures being collected.

## Integration

- **Reads:** nothing at runtime. The change is prose plus one new table.
- **Writes:** `docs/src/content/docs/performance/threading.md`, a new `docs/perf_data/` subdirectory,
  and `CHANGELOG.md`.
- **Downstream:** `docs/` publishes to www.trimgalore.com from `dev`, so this lands via the normal
  dev flow. The audience that matters is workflow authors sizing a scheduler budget — nf-core,
  Snakemake, CWL — which is the citation path that produced #11541.
- **No source, CLI, flag, or output change.** The `validation` CI job and all byte-identity
  invariants are untouched; `--memory` remains clumpify-only.

## Assumptions

- **Fixed:** docs-only. No code change, no new flag, no change to `--memory`'s scope.
- **Fixed:** the existing single-end table stays as it is; the paired figures are added beside it, not
  merged into it, so neither set is silently re-attributed to the other's fixture.
- **Fixed:** register is plain scientific prose per the project's docs-register guidance — no
  marketing adjectives, and no rewriting of neighbouring legacy sections while in the file.
- **To confirm before publishing:** that the existing table's figures are Linux. They are attributed
  to `/usr/bin/time -v`, which is GNU time and therefore almost certainly Linux, but the audit trail
  in `docs/perf_data/` should confirm it. If they are not Linux, both tables need explicit platform
  labels rather than one shared note.
- **Not for publication:** the macOS/aarch64 figures gathered while diagnosing rnaseq#1903 (paired
  150 bp: 9 MiB at `--cores 1`, 116 at 4, 202 at 8, 304 at 16; 190 MiB at `--cores 8` with
  `--fastqc_args '-t 24'` on 12M unique reads; plateau confirmed between 4M and 12M pairs). They are
  indicative and set the expected order of magnitude for step 1, but publishing them beside Linux
  figures in one section would mix platforms.
- The page-cache mechanism is the leading explanation for rnaseq#1903's ~1022 MiB `MaxRSS` against a
  measured ~200 MiB heap, but it is **inference, not measurement** — no cgroup reproduction was run.
  The docs must therefore describe how cgroup accounting works (which is documented behaviour) without
  asserting it as the confirmed cause of any particular incident.

## Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| 1 | Every published figure came from a successful run | Assert exit 0 per run before recording; keep the exit status beside each figure in the perf_data file | No figure recorded from a non-zero exit. A failed run still emits a parseable time log, so the status must be checked separately from the value |
| 2 | The disk guard actually guards | Run the free-disk check with a deliberately impossible threshold | Fails loudly. Without a positive control the guard may be a no-op, and a stalled writer inflates peak RSS |
| 3 | The plateau claim holds | Compare `--cores 8` at both input sizes | Within ~10%; if it still grows, the larger input was not yet saturating and the claim must be re-scoped |
| 4 | FastQC threads are not additive | Compare `--cores 8` with and without `--fastqc_args '-t 24'` | Difference small (order tens of MiB). If large, Behavior 3 is wrong and must be rewritten |
| 5 | The invited-misreading phrasing is gone | `grep -n "noise floor" docs/src/content/docs/performance/threading.md` | No match |
| 6 | The recommendation is defensible in both directions | Check it against measured peak and against a typical 30M-pair output volume | Comfortably above measured peak at `--cores 8`, and not so large it reintroduces the over-allocation #11541 was reducing |
| 7 | No contradiction with #454 | Re-read the read-length and input-size paragraphs after editing | Long-read ceiling and plateau statements remain consistent |
| 8 | Docs build | Build the Astro docs site | Builds clean; both tables render |

## Questions or ambiguities

- **Open — CI drift guard.** Should the `clumpify-memory` job be extended to cover the plain paired
  path, so these published figures are guarded against drift the way the clumpify baselines are?
  Recommend yes, as a separate follow-up: it is a CI change, not a docs change, and should not gate
  this one. Assumption taken: out of scope here.
- **Open — cite the incident?** Reference nf-core/rnaseq#1903 and modules#11541 in the docs as a
  worked example, or keep the prose impersonal? Recommend one neutral sentence naming the failure mode
  (a 1 GB cgroup cap with an 8-cpu band) without attributing fault, since the mechanism is the
  instructive part. Assumption taken: mechanism described, incident not named.
- **Open — table shape.** A second table (recommended, and assumed here) versus adding a paired
  column to the existing one. A second table keeps each fixture's provenance attached to its own
  figures; a merged table is more compact but re-attributes rows.
- **Open — where the recommendation lives.** `### Budgeting on a cluster` under Memory profile
  (assumed), versus a standalone top-level section that workflow authors can link to directly. If
  nf-core is expected to cite it, a linkable top-level heading may serve better.

## Self-Review

**Logic.** The chain from goal to change holds: the published figure was single-end and heap-only, a
consumer sized a budget from it, the budget was breached. Each of the three additions closes one link.
Step 4 is the one that matters most — adding a paired table while leaving *"well below the noise floor
of cluster scheduler memory allocation"* in place would leave the sentence that produced the problem
untouched, so the edit is explicitly called out rather than left implicit in "extend the section".

**Adjusted during review.** Three changes. (a) Step 1 originally reused the replicated 150 bp fixture;
that repeats identical reads, which understates per-sequence state in FastQC's duplicate tracking, so
it now requires high-diversity reads. (b) Validation 1 originally checked only the recorded value;
split into value-and-exit-status, because a failed run still writes a parseable time log and a
plausible number read out of one is the specific way this measurement goes wrong. (c) Validation 2
(positive control on the disk guard) was added after noticing ci.yml treats disk and memory as
coupled — a guard that silently never fires gives false confidence.

**Edge cases.** `--cores 1` (no batch residency), long reads (defer to #454), `--clumpify` (separate
`--memory` model), uBAM output (single-threaded, so `--cores 1` profile). All routed in Behavior rather
than left to the implementer.

**Integration.** No source or CI files touched, so no byte-identity or validation-matrix exposure. The
one coupling worth watching is prose consistency with #454, which landed the same day — Validation 7
covers it.

**Remaining risks.** Two. First, the page-cache explanation is inference; the mitigation is to document
how cgroup accounting works rather than to assert a cause, and Assumptions says so explicitly. Second,
the recommended number is a judgement call balancing headroom against the over-allocation #11541 was
reducing — Validation 6 checks it from both directions, but it stays a judgement call and is worth a
maintainer's eye rather than an implementer's.
