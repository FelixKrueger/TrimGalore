# Plan: a peak-RSS baseline for the plain trimming path in CI

**Date:** 2026-08-19
**Type:** CI (workflow + a new baseline data file). No library or binary source change.

## Goal

Give the plain (non-`--clumpify`) trimming path the same kind of peak-RSS guard that `--clumpify`
already has, so a change that moves plain-path memory is visible on the PR that makes it.

Two real events motivate this, both found this week and neither by CI:

- `2ab0918` (#437, mimalloc) raised plain paired peak RSS ~1.9x at `--cores 8` — 74 → 146 MiB on the
  65 bp Buckberry pair. It sat on `dev` for three days.
- `1595d59` (#442) changed a per-batch output-buffer reservation on the **shared** path while scoped to
  cutting `--clumpify` memory. Under the platform allocator that cost ~18 MiB at `--cores 4`; under
  mimalloc it saves ~15 MiB at `--cores 8`.

Both were found only because an unrelated nf-core OOM (nf-core/rnaseq#1903) sent someone looking. The
`clumpify-memory` job would not have caught either: it measures only cells that engage `--clumpify`.

## Context

**Files:**

- `.github/workflows/ci.yml` — the `clumpify-memory` job (from ~line 271) is the model to mirror.
- `.github/plain-memory-baseline.tsv` — new, alongside `.github/clumpify-memory-baseline.tsv`.
- `examples/mk_clump_fixture.rs` — the existing fixture generator (`cargo build --example`).

**What the model job already solves, and must be reused rather than reinvented:**

1. **GNU-time capability probe** with an install fallback, on the recorded reasoning that *"presence is
   not capability"* — the runner ships a `/usr/bin/time` whose `-v` is unrecognised.
2. **A `>= 8 GiB` free-disk assertion**, because *"disk is not independent of memory: a stalled writer
   backs up the result channel and the pending map, so an out-of-space run reports inflated peak RSS"*.
3. **Capture, check status, then parse** — never pipe into `grep`, because *"a parseable log from a
   failed run looks like data"*.
4. **Positive-engagement assertions** per cell — the clumpify job greps `^clumpify:` so a silent
   fallback cannot pass as a measurement.
5. **A `BOOTSTRAP` marker** in the baseline file that warns and suspends enforcement. This is the
   mechanism for the first pass; it does not need inventing.
6. **Empty-field parse guards** and KiB-side comparison to avoid a second truncation.

**The structural difference that shapes this plan.** The clumpify job gates on two things: the
*promise* (measured peak <= the peak the binary printed) and a *drift band* on peak/prediction. The
plain path **prints no predicted peak** — `--memory` governs `--clumpify` only. So:

- there is no absolute, platform-independent gate available here;
- the only gate is drift against a recorded **absolute** peak, in KiB;
- therefore the baseline is tied to `ubuntu-latest` *and* to the allocator in the tree, and both must be
  recorded in the file.

**The banding rule this must obey**, quoted from the existing baseline file:

> Band a cell only when the width satisfies both sides: `width >= 3x the observed half-range` (so noise
> cannot fire it), and `width <= the measured regression delta` minus a working margin (so the
> regression still fires it). The half-range must come from a double-digit sample. Two or four draws
> estimate dispersion LOW by construction.

The 3-rep figures gathered on `dockyard-oxy-0` this week therefore **cannot** set a CI band: wrong host,
and far too few draws. This is why the first pass records and does not enforce, matching the precedent
in `033cf1d` and `8e64346` (measure a cell rather than record it as a blind spot, then band it once
enough readings exist).

## Behavior

1. A job runs on `ubuntu-latest` (4 vCPU), in parallel with the rest, on no one's critical path.
2. It builds `trim_galore` and `mk_clump_fixture`, then asserts the harness can measure at all
   (GNU time reporting peak RSS; `>= 8 GiB` free) before generating anything.
3. It generates two fixtures with the existing generator: 10M pairs at 51 bp, and 10M pairs at 150 bp.
4. For each cell it runs `trim_galore` under `/usr/bin/time -v`, and records peak RSS in KiB.
5. Per-cell assertions, each of which must fail the cell rather than silently pass:
   - exit status zero, captured separately from the parse;
   - the startup banner says `Using N worker threads` for the N the cell passed — the plain-path
     analogue of the `^clumpify:` check, so a silent single-threaded fallback cannot pass;
   - output is actually gzipped (`ls "$out"/*.fq.gz`), because `--gzip` does not force gzip: compression
     follows `input[0]`, and a plain-text output would understate the peak at exit 0;
   - the parsed peak field is non-empty.
6. It prints every cell's peak in KiB and MiB in a copy-pasteable block for the baseline file.
7. Gating, once out of BOOTSTRAP: fail if a banded cell's peak rises above `baseline + width`, or falls
   below `baseline - width`. **The lower bound matters as much as the upper**: an improvement is also a
   change that invalidates the baseline, and `1595d59` was an improvement under mimalloc.
8. While `# BOOTSTRAP` is present, it warns and does not fail on drift — but every assertion in step 5
   still fails the job, since those are correctness checks, not drift checks.

**Edge cases:**

- `--cores 1` bypasses the worker pool entirely (~9 MiB), a genuinely different code path. Cheap to
  measure and worth a cell.
- The generator emits fixed-seed random ACGT, which deflate barely compresses, so every figure is an
  upper bound relative to real reads. Fine for a baseline, which needs consistency rather than realism,
  but it must be stated in the file or someone will compare these numbers to the docs' and find a
  contradiction.
- A missing baseline file must fail loudly; a label with no row must report that its band is not armed.
- An allocator change is a mandatory re-baseline event, not a red to be re-run.

## Implementation outline

1. **Add `.github/plain-memory-baseline.tsv`** with a `# BOOTSTRAP` line, the cell list, and a header
   recording: host (`ubuntu-latest`), allocator, fixture generator and its parameters, and the banding
   rule this file inherits. No baseline values yet.
2. **Add a `plain-memory` job to `ci.yml`**, mirroring `clumpify-memory`'s structure: its own cargo
   cache key, the capability probe, the disk assertion, fixture generation, then measure-and-gate.
   Reuse the model job's shell idioms verbatim where they apply rather than paraphrasing them.
3. **Cells** (label — invocation — why):
   - `plain_c4` — `--paired --cores 4`, 51 bp — the headline cell.
   - `plain_c2` — `--paired --cores 2`, 51 bp — guards the per-worker growth term.
   - `plain_c1` — `--paired --cores 1`, 51 bp — the worker-pool-bypass path.
   - `plain_c4_150bp` — `--paired --cores 4`, 150 bp — the read-length axis, which is the axis that
     produced #1903.
   - `plain_c4_fastqc` — `--paired --cores 4 --fastqc`, 51 bp — **record only, do not band.** The
     clumpify file already found its FastQC cell too wide to band (0.080 range, 1.2x at +0.09) because
     the allocator holds pool pages across the trim→FastQC boundary. Expect the same here.
4. **Bootstrap pass:** merge with BOOTSTRAP set, let it run on `dev` pushes and PRs, and collect
   readings. Do not set a band from the first few.
5. **Banding pass (separate PR):** once a cell has a double-digit sample, compute mean and half-range,
   choose a width satisfying both sides of the rule, record the arithmetic in the file the way the
   clumpify baseline does, and delete `# BOOTSTRAP`.
6. **Document the width's regression side** using the two known deltas: mimalloc at ~+72 MiB (c8) and
   the reservation change at ~+18 MiB (c4, stock allocator). The narrower of these is what the band must
   still catch, so the width has an evidenced upper bound rather than a guessed one.

## Efficiency

Runtime is dominated by fixture generation (two 10M-pair fixtures) and five measured cells. The model
job allows 30 minutes for four cells plus a cold cargo build; this should be budgeted similarly.
Measure serially — concurrent runs on one runner contend for page cache and would corrupt the figures.

If the wall-clock cost is unwelcome, the fallback is to drop `plain_c1` and `plain_c2` and keep the
`c4` pair, since `c4` is where a per-worker regression shows most clearly at 4 vCPU.

## Integration

- **Reads:** `.github/plain-memory-baseline.tsv`.
- **Writes:** nothing in the repo at runtime; logs and a copy-pasteable block only.
- **No source change**, so no exposure to the Perl 0.6.11 byte-identity `validation` matrix.
- Runs in parallel with existing jobs, with no `needs:`, so it does not extend the critical path.
- Deliberately **not** folded into `clumpify-memory`: that job's baseline file, header commentary and
  gating semantics are all built around peak/prediction ratios, which do not exist here. See the open
  question — sharing the job would save one build and one fixture generation.

## Assumptions

- **Fixed:** CI-only. No library, binary, flag, or output change; `--memory` stays `--clumpify`-scoped.
- **Fixed:** the first pass does not enforce drift. Anything that claims a band from fewer than ten
  readings per cell contradicts the rule the existing baseline file states explicitly.
- `ubuntu-latest` is 4 vCPU, so `--cores 8` is not a meaningful cell here; the interesting core-count
  axis on this runner is 1/2/4.
- 10M pairs saturates the plain path's resident capacity at short read lengths. Supported by this
  week's measurements — peak was flat between 4M and 12M pairs at 150 bp — but at 10M and 51 bp it is
  an inference and step 4 should confirm the cells are not input-limited before any band is set.
- Baselines recorded now are **mimalloc** baselines, because that is what is in `dev`. An allocator or
  `Cargo.lock` change invalidates them.
- The peak the job records is `ru_maxrss`, which is anon RSS. It is not the quantity a cgroup scheduler
  charges — that includes page cache. This job guards our own allocations, not a cluster budget.

## Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| 1 | The job fails a genuinely broken run | Point a cell at a nonexistent input | Fails on exit status, not on an empty parse |
| 2 | The worker-banner assertion discriminates | Run a `--cores 4` cell with the banner check set to `Using 8` | Cell fails. Without this control the check may be a no-op |
| 3 | The gzip assertion discriminates | Feed one cell a plain-text fixture | Fails rather than recording an understated peak at exit 0 |
| 4 | The disk guard fires | Run the free-space check against an impossible threshold | Fails loudly |
| 5 | BOOTSTRAP genuinely suspends only drift | Commit an absurd baseline value with BOOTSTRAP set | Job warns, still green; with BOOTSTRAP removed, red |
| 6 | Cells are not input-limited | Measure `plain_c4` at 10M and 20M pairs | Within noise. If the larger is higher, 10M is not saturating and the fixture must grow before banding |
| 7 | The band would catch the known events | Replay the two measured deltas against the chosen width | Both exceed it; if the reservation delta (~18 MiB) does not, the width is too wide to be useful |
| 8 | The lower bound is armed | Check an improvement fires it | Red, forcing a re-baseline |

## Questions or ambiguities

- **Open — separate job or extend `clumpify-memory`?** Assumed separate, for clean gating semantics and
  an independent baseline file. Extending it would save one cargo build and one fixture generation
  (~several minutes per run) at the cost of mixing ratio-gated and absolute-gated cells in one file.
  Worth overriding if CI minutes matter more than file clarity.
- **Open — which read lengths?** 51 bp and 150 bp assumed: 51 matches the clumpify job's fixture so the
  two are comparable, and 150 is what real short-read pipelines send. A 25 bp cell would match the
  other clumpify fixture but adds a cell for an axis the plain path has no known sensitivity to.
- **Open — should this gate PRs, or only report?** Assumed gating once banded, matching
  `clumpify-memory`. A report-only job is easier to live with but is the state we are already in.
- **Open — how to handle the allocator dependency mechanically?** Simplest is a header note plus
  reviewer discipline. A stricter option is to hash the allocator-relevant `Cargo.toml`/`main.rs` lines
  into the baseline file and fail when it changes, which is more machinery than the clumpify job uses
  for the same class of risk.

## Self-Review

**Logic.** The chain holds: two real plain-path memory events shipped unseen, the existing job cannot
see them by construction, so the gap is a plain-path measurement. The design's one genuine weakness is
that no absolute gate is available without a printed prediction, which is why the assertions in
Behavior 5 carry more weight here than in the model job — they are the only checks that work on the
first run, before any band exists.

**Adjusted during review.** Three changes. (a) The first draft banded cells using this week's oxy
figures; removed once I re-read the baseline file's own rule that a half-range from a handful of draws
is set too tight — and those draws were from the wrong host besides. (b) Added the lower-bound gate and
Validation 8 after noticing that `1595d59` was an *improvement* under mimalloc: a one-sided band would
have stayed green through a change that invalidated the baseline. (c) Added Validation 6 because the
plan asserts 10M pairs saturates, and an input-limited cell produces a baseline that moves when the
fixture is regenerated rather than when the code changes.

**Edge cases.** `--cores 1` bypass, random-ACGT incompressibility, missing baseline file, unarmed
labels, and allocator changes are all routed in Behavior or Assumptions rather than left to the
implementer.

**Integration.** No `src/` change, so no byte-identity exposure. The one coupling to watch is that
this job and `clumpify-memory` both build the same binary and generate overlapping fixtures; the open
question covers whether to share.

**Remaining risks.** Two. GitHub-hosted runner variance may prove wide enough that no width satisfies
both sides of the banding rule for some cells — in which case those cells are recorded, not banded,
exactly as the clumpify file does for its two widest. And the plan cannot say in advance how many runs
the bootstrap needs; "double-digit" is the floor the existing file sets, and the honest answer is that
step 5 waits until the data supports a width.
