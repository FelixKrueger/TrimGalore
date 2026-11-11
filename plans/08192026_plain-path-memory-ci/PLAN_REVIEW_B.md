# Plan Review B — `plans/08192026_plain-path-memory-ci/PLAN.md`

**Reviewer:** B (independent)
**Date:** 2026-08-19
**Target:** a peak-RSS baseline for the plain trimming path in CI
**Verdict:** sound premise, correctly read model job, **two defects that would break the job as
written** (one hard red on every run, one unspecified schema), and **one large missed opportunity**: a
platform-independent gate *is* available, so the job need not be inert for the whole bootstrap window.

Everything below was checked against the tree, the real `clumpify-memory` job, its last successful run,
and local measurement. Sources of each number are named.

---

## 0. Evidence gathered for this review

**Local measurement** (macOS/arm64, mimalloc, `target/release/trim_galore`, 1M pairs, fixtures from
`examples/mk_clump_fixture`):

| cell | wall | peak RSS | banner present |
|---|---|---|---|
| `--paired --cores 4`, 51 bp | 1.20 s | **9,420,800 B → 9.0 MiB** at c1 | — |
| `--paired --cores 1`, 51 bp | 4.02 s | c1 = **9.0 MiB** | **no** |
| `--paired --cores 2`, 51 bp | 2.00 s | c2 = **49.3 MiB** | yes |
| `--paired --cores 4`, 51 bp | 1.20 s | c4 = **79.9 MiB** | yes |
| `--paired --cores 4 --fastqc`, 51 bp | 3.07 s | — | yes |
| `--paired --cores 4`, 150 bp | 3.04 s | — | yes |

Fixture sizes per 1M pairs: 51 bp = 24.44 MB/mate; 150 bp = 60.62 MB/mate. Generation 1.49 s / 3.76 s.

**Real runner data**, last successful `dev` run (32179108913), job *Clumpify peak RSS (#439)*, total
**5 minutes**:

| step | duration |
|---|---|
| Build binary and fixture generator | **73 s** |
| Assert the harness can measure at all | 0 s |
| Generate fixtures (10M×51 + 10M×25) | **30 s** |
| Measure and gate (7 cells, one of them multi-pair) | **214 s** |

That step's log also prints `free in RUNNER_TEMP: 88190 MiB` and `peak-RSS probe OK: Maximum resident
set size (kbytes): 1260` — i.e. **86 GiB free**, and the preinstalled `/usr/bin/time` on the current
image *does* support `-v` (the apt fallback did not fire).

Repo facts: `FelixKrueger/TrimGalore` is **PUBLIC** (GitHub-hosted minutes are free), and **neither
`dev` nor `master` has branch protection** — `gh api .../branches/{dev,master}/protection` → 404
`Branch not protected`. There are no required status checks in this repo today.

---

## 1. Logic review

### 1.1 The central premise is TRUE — verified in source

> "The plain path prints no predicted peak — `--memory` governs `--clumpify` only."

Confirmed. `src/main.rs:178` `resolve_clump_layout()` returns `Ok(None)` immediately when
`!cli.clumpify`, before any print. The only `predicted peak ≈` sites in the tree are
`src/main.rs:215` (inside that function, after the clumpify guard) and four in `src/clump_only.rs`
(:301, :450, :800, :1050) — all `--clump_only`. `grep -niE 'memory|rss|peak' src/report.rs` returns
nothing, so neither the text nor the JSON report carries a memory field either. There is no printed
absolute promise to gate on, and the plan is right that this is the structural difference from the
model job.

**But "no absolute gate" ≠ "no platform-independent gate."** See §4.1 — the between-cell ordering is
platform-independent, armable on run #1, and my measurements show it has 62% headroom at its tightest
point. The plan's Self-Review calls the missing absolute gate "the design's one genuine weakness" and
then accepts it; it is partly recoverable, and recovering it is the difference between a job that
guards something in week 1 and a job that guards nothing until roughly week 3.

An *analytic* prediction (the other way to get an absolute gate) genuinely is not viable, and it is
worth writing that down so nobody re-proposes it: `src/parallel.rs` has exactly one size constant,
`BATCH_SIZE = 4096` (line 33). The analytic term for c4/51 bp is ≈ 4 workers × depth 2 × 2 mates ×
4096 × ~150 B ≈ 14 MiB against a measured 80 MiB. Five-sixths of the peak is allocator arenas and
deflate state, so a formula-derived ceiling would sit ~5x below the truth and gate nothing.

### 1.2 **CRITICAL** — the worker-banner assertion makes `plain_c1` fail on every run

Behavior 5 applies the banner assertion to every cell:

> the startup banner says `Using N worker threads` for the N the cell passed … so a silent
> single-threaded fallback cannot pass

`src/main.rs:858` guards it:

```rust
if cli.cores > 1 {
    eprintln!("Using {} worker threads (parallel trim + compress)", cli.cores);
}
```

At `--cores 1` **nothing is printed**. I confirmed empirically:
`grep -c "Using .* worker threads"` → `1` at c4, `1` at c2, `0` at c1. So `plain_c1` — a cell the plan
argues for explicitly, correctly, as the worker-pool-bypass path — fails its own engagement assertion
on every run, in BOOTSTRAP mode too, because Behavior 8 keeps step-5 assertions armed. The job is red
from the first merge.

The fix is a per-cell expectation, not one rule: assert *presence* of `Using $cores worker threads` for
cells with `cores >= 2`, and assert *absence* of `Using ` for `plain_c1`. The absence form is the
better check anyway — it is what proves the c1 cell took the serial branch rather than being silently
promoted.

### 1.3 **CRITICAL** — the banner is a flag echo, so it does not prove what the plan claims

Even where the banner exists, it is printed from `cli.cores`, i.e. from the parsed flag, before any
pool is constructed. `src/parallel.rs` prints nothing at all (`grep -n 'eprintln!\|println!'
src/parallel.rs` → no matches). So `Using 4 worker threads` proves clap parsed `--cores 4`. It does not
prove four workers ran.

This is materially weaker than the model's `^clumpify:` grep, which is emitted from the *resolved*
`ClumpLayout` after `resolve_layout()` succeeded and therefore cannot print if clumpify silently fell
back (`resolve_clump_layout` returns `None` before the print on the too-small-budget path). The plan
calls the banner "the plain-path analogue of the `^clumpify:` check". It is not an analogue; it is a
weaker instrument, and the plan's confidence in it matters because Self-Review leans on step-5
assertions as "the only checks that work on the first run".

The real analogue is the ordering invariant in §4.1: it is derived from the measured quantity itself,
so nothing about it can be satisfied by a flag that was parsed and then ignored.

### 1.4 **IMPORTANT** — the gzip assertion cannot fail as the job is constructed

`src/io.rs:41` `output_is_gzipped(dont_gzip, first_input)` = `!dont_gzip && is_gzipped(first_input)`,
and `is_gzipped` is **extension**-based (`src/io.rs` tests at :1636-1657 pin exactly that, including
`sample.gz.fastq` → false). The plan's reasoning is right on the mechanism, and the `--gzip` claim
checks out too (`src/cli.rs:1187`/`:1191` — warn and ignore, both branches).

But the job generates its own inputs, and `mk_clump_fixture` can only ever write
`<prefix>_R1.fastq.gz`. With `--dont_gzip` never passed, `ls "$out"/*.fq.gz` is true by construction.
It can only fire if a future edit changes the fixture path or adds `--dont_gzip` — which is a real if
small guard-rail, so keep it (it is one line and free). Just do not carry it in the plan as a
"positive-engagement assertion" beside the banner: of the four step-5 assertions, exactly one (exit
status) is unconditionally load-bearing, one is broken (§1.2), one is a flag echo (§1.3), and this one
is presently unfalsifiable.

I verified the glob itself is right: the plain paired path writes `f51_R1_val_1.fq.gz` /
`f51_R2_val_2.fq.gz`, which `*.fq.gz` matches.

### 1.5 **IMPORTANT** — the baseline file schema for widths is unspecified, and the model's shape does not carry over

The model hardcodes one global width in the shell (`+ 0.09` / `- 0.10`) and keeps the TSV at two
columns (`cell<TAB>baseline_peak_over_prediction`). That works because every clumpify cell's gated
quantity is a *ratio* around ~0.8-1.0, so one width fits all cells.

The plain job's quantity is absolute KiB, and my measurements span **9.0 MiB (c1) to ~80 MiB (c4) to
several hundred MiB (c4_150bp)** — nearly two orders of magnitude. A single global width therefore
cannot work in either form:

- **absolute** (e.g. ±20 MiB): more than 2x the entire c1 cell, so c1 is unbandable, while being far
  inside noise-plus-drift for c4_150bp;
- **fractional** (e.g. ±15%): the plan's own regression side breaks. The two evidenced deltas are
  absolute (+72 MiB, +18 MiB). +18 MiB is +24% of a 74 MiB c4 base but a much smaller fraction of the
  150 bp cell, so a fraction tuned to catch it at c4 will not catch it at c4_150bp — and the plan
  names read length as "the axis that produced #1903".

The plan says "fail if a banded cell's peak rises above `baseline + width`" without ever saying where
`width` lives or what unit it is in. Decide now and write it into the outline: a **three-column TSV**
(`cell<TAB>baseline_kb<TAB>width_kb`), per-cell absolute widths, arithmetic recorded in the header the
way the clumpify file does. Otherwise the banding PR inherits the ambiguity at the point where it is
most expensive to discover.

### 1.6 **IMPORTANT** — during BOOTSTRAP the job guards nothing, and BOOTSTRAP is weeks long

Behavior 8 is explicit: while `# BOOTSTRAP` is present, drift does not fail. Combined with §1.2-1.4,
the first-merge job asserts exit status and (broken) engagement, and *nothing at all* about memory.

How long? The clumpify file's own history is the calibration: banded "from 27 observations per cell
(31 job logs)" after being introduced in #439 — roughly two weeks of pushes and PRs. So the plan ships
a job that, for about two weeks, would not have caught either motivating event. `2ab0918` sat on `dev`
for three days.

This is fixable without violating the banding rule the plan quotes. The rule's noise side is
`width >= 3x the observed half-range`; the clumpify cells' *observed* half-ranges are ±2-4% (ranges
0.034-0.081 on means 0.71-1.03). A **coarse tripwire at ±50%** of the first reading is therefore
~12-25x half-range — it clears the noise side by an order of magnitude and needs no double-digit
sample to justify. Its regression side: it catches the mimalloc 1.9x event (the exact class of event
that motivated the plan) and misses the +24% reservation change, which the fine band will catch later.
Half a loaf, armed on day one, with arithmetic that satisfies the rule rather than dodging it.

### 1.7 **IMPORTANT** — the plain path is touched by most PRs, so the observation pool is much harder to accumulate than the clumpify one

The clumpify baseline file closes with:

> When `src/` moves, judge whether the change can reach the clump path before pooling new observations
> with these; the proxy `git diff --quiet <old> <new> -- src Cargo.lock` is deliberately conservative
> and exists to force that judgement, not to make it.

For the clump path that judgement usually answers "no". For the **plain** path — the shared trim +
compress path — it usually answers "yes": `1595d59` is the plan's own example of a clumpify-scoped PR
that moved the shared reservation. So the plain job's readings will be invalidated far more often than
the clumpify job's, and Implementation step 4 ("let it run on `dev` pushes and PRs, and collect
readings") quietly assumes poolability that the plain path does not have.

Say how the pool is defined. My suggestion: pool only `push`-on-`dev` readings, record the commit sha
beside each reading, and prefer the within-run ratios of §4.1 for early signal, since a ratio between
two cells measured in the same binary is far more robust to code drift than either absolute number.

### 1.8 **IMPORTANT** — a thin adapter-detection margin is an unrecorded baseline mover

Not in the plan at all, and it is a live risk on random-ACGT fixtures. Measured:

- 51 bp fixture: `Using Nextera adapter for trimming (count: 5). Second best hit was smallRNA (count: 2).`
- 150 bp fixture: `Using Nextera adapter for trimming (count: 10). Second best hit was smallRNA (count: 6).`

Both cells trim with Nextera, and it is deterministic for a fixed fixture (fixed-seed SplitMix64), and
the detection window is the first 1M reads regardless of the fixture's pair count — so there is no
run-to-run flake. But 5-vs-2 is a two-hit margin on noise. Any change to the detection heuristic, the
scan window, the built-in adapter set, or the fixture parameters can flip the detected adapter, which
changes trimming work and hence the peak — and it would present as unexplained baseline drift, not as
a detection change.

One line fixes it: assert `Adapter: Nextera` per cell alongside the other step-5 checks, so a flip is
a red with an obvious cause. I would keep auto-detection rather than pinning `--nextera`, because the
1M-read detection scan is part of the plain path's memory profile and pinning would silently remove it
from what the job measures.

### 1.9 Two model-job idioms that must NOT be copied

The plan instructs "reuse the model job's shell idioms verbatim where they apply". Two traps:

- **The prediction-margin assertion** (`pred_kb == budget_kb * 10 / 11 ± 2048`) is entirely
  clumpify-specific. There is no prediction here; drop it, and say so, because "verbatim where they
  apply" invites an implementer to keep it and then invent a plain-path meaning for it.
- **The empty-field guard's shape.** `case "${peak_kb}:${pred_mib}:${budget_mib}" in ''|*::*|:*)` is a
  three-field pattern. With one field it degenerates to `[ -z "$peak_kb" ]`. Copying the `case`
  verbatim still happens to work (`''` matches; the other two arms are unreachable), but the
  unreachable arms then read as if they check something. Write the one-field form.

Also, `measure` mutates `fail` in the caller's shell; it only works because the model calls it
unpipelined at top level. If a plain-path cell is ever written as `measure ... | tee`, `fail=1` is lost
in the subshell and every red becomes green. One caution line in the outline is worth it.

### 1.10 Correctly read — no action needed

I checked each reuse claim in Context and each cited fact:

- GNU-time capability probe with apt fallback, and the "presence is not capability" reasoning — exists
  verbatim, `ci.yml` "Assert the harness can measure at all".
- `>= 8 GiB` free-disk assertion and the stalled-writer reasoning — exists (`test "$free_kb" -ge
  8388608`), quoted accurately.
- Capture, check status, then parse; "a parseable log from a failed run looks like data" — exists,
  quoted accurately.
- `^clumpify:` positive-engagement grep — exists.
- `BOOTSTRAP` marker that warns and suspends the band — exists, and the `want`-lookup gate is indeed
  skipped when `bootstrap=1`, so Validation 5's expectation is right.
- KiB-side comparison — exists, with the unit-slip rationale.
- The banding-rule block quoted from the baseline file — accurate, including the double-digit-sample
  clause.
- `1G_c4_fastqc 27 0.7069 0.688 0.768 0.080 1.2x at +0.09` and the trim→FastQC pool-page reasoning —
  accurate; `plain_c4_fastqc` as record-only is the right call.
- `--gzip` warns and is ignored; compression follows `input[0]` — accurate (`src/cli.rs:1187`/`:1191`,
  `src/io.rs:41`).
- `--cores 1` is a genuinely different path — accurate: `src/main.rs:1382`/`:1614`/`:1900` all branch
  on `cli.cores > 1 || cli.clumpify`.
- The four cited commits all exist with matching subjects: `2ab0918` "perf(alloc): use mimalloc as the
  global allocator (#437)", `1595d59` "perf(clump): cut --clumpify peak RSS 13-17%…(#442)", `033cf1d`,
  `8e64346`. `1595d59`'s message confirms the plan's characterisation: "Each per-batch output buffer is
  reserved from its own mate's record bytes instead of a flat 300 bytes per record" — shared path,
  clumpify-scoped PR.
- The `~9 MiB` figure for `--cores 1` — measured 9,420,800 B = 8.98 MiB. Exact.

One omission in this class: the new job needs `if: github.event_name != 'schedule'`. The workflow has a
weekly cron (`0 6 * * 1`) whose stated purpose is the `audit` job, and every non-audit job carries that
guard. Without it the new job burns ~7 minutes every Monday and can post a Monday-morning red.

---

## 2. Assumptions

| Plan assumption | Verdict |
|---|---|
| CI-only, no source change | Holds. Nothing proposed touches `src/`, so no Perl-0.6.11 byte-identity exposure. |
| First pass does not enforce drift | Holds as stated, but see §1.6 — it need not mean *no* memory gate. |
| `ubuntu-latest` is 4 vCPU, so `--cores 8` is not meaningful; 1/2/4 is the axis | Correct today. Note the consequence the plan does not: `--cores 4` on 4 vCPU **is** oversubscribed once the OS and the reader/writer threads are counted, so queue occupancy — the thing being measured — is CPU-contention-sensitive. Empirically tolerable: the model's `1G_c4` shows range 0.055 on mean 0.788 (±3.5%) across 27 observations on the same runner class. Worth one line in the baseline header rather than a design change. |
| 10M pairs saturates at short read length | Correctly flagged as an inference, and Validation 6 is the right control. Note it is a weaker assumption here than for clumpify: there is no bin pool to fill on the plain path, and resident capacity is `cores × depth × BATCH_SIZE` records — which is ~32K records at c4, i.e. saturated within the first second at any fixture size. My c1 figure (9.0 MiB at 1M pairs) is already at the plan's quoted 9 MiB steady state. **Prediction:** Validation 6 will show 10M and 20M within noise, and 10M is likely far more than needed. If so, say it in the header — the next person will otherwise assume the fixture size is load-bearing and be afraid to shrink it. |
| Baselines are mimalloc baselines; allocator/`Cargo.lock` change invalidates them | Correct. `Cargo.toml:50 mimalloc = "=0.1.52"`, `src/main.rs:29-30` `#[global_allocator]`. |
| `ru_maxrss` is anon RSS, not what a cgroup charges | Correct and worth keeping, especially as #1903 was a container OOM. |

Unstated assumptions worth surfacing:

- **That a red will be acted on.** Neither `dev` nor `master` has branch protection (verified, 404), so
  no job in this repo mechanically blocks a merge. This cuts both ways and is the factual answer to
  open question 3 — see §3.2.
- **That the readings are poolable across commits** — §1.7.
- **That the detected adapter is stable** — §1.8.
- **That the baseline file's width lives somewhere** — §1.5.

---

## 3. Efficiency analysis

### 3.1 The cost concern is not real, and the plan's fallback is a bad trade

The plan writes: "The model job allows 30 minutes for four cells plus a cold cargo build". The model job
measures **seven** cells (`1G_c4`, `1G_c4_fastqc`, `25bp_1G_c4`, `multipair_c4` — itself two pair-runs —
`1G_c2`, `25bp_1G_c2`, `2G_c4`). "Four" is a stale comment inside `ci.yml` itself, which the plan has
inherited. The error is in the plan's favour, so no harm, but the budget should be argued from measured
numbers, and they are available:

Model job, last successful `dev` run: build 73 s, fixtures 30 s, measure 214 s, **5 min total**.
That is ~24 s per 10M-pair *clumpify* cell. My local timings put the plain cells at 12/20/40/30/31 s
per 10M pairs (c4/c2/c1/c4_150/c4_fastqc); the runner ran the clumpify fixture generation at ~2x my
local speed-adjusted rate, so scaling by ~2x:

| plain job step | estimate |
|---|---|
| Build (warm cache) | ~75 s |
| Fixtures (10M×51 + 10M×150) | ~70 s |
| Measure, 5 cells | ~270 s |
| **Total, warm** | **~7 min** |
| Total, cold cargo cache | ~12 min |

Against a 30-minute timeout that is comfortable, matching the model. **Recommendation: keep all five
cells and delete the Efficiency paragraph's fallback.** Dropping `plain_c1` to save ~80 s would remove
the only serial-path cell *and* the anchor of the ordering gate in §4.1; dropping `plain_c2` removes
the per-worker growth term the plan itself says the cell exists for. Both are false economies at these
numbers.

### 3.2 Separate job vs extending `clumpify-memory` — the open question has a decisive answer

Separate. The arithmetic:

- Extending saves one build (73 s) plus job setup (~15 s) ≈ **1.5 min of runner time**, and costs a
  single ~11-minute job instead of two ~7-minute jobs running in parallel.
- The repo is **PUBLIC**, so GitHub-hosted minutes are free. There is no minutes budget to protect.
- Neither job is on anyone's critical path (`rust-tests` alone allows 45 min and took 8-11 min).

So the entire "at the cost of mixing ratio-gated and absolute-gated cells in one file" trade buys
nothing. Close the question in favour of separate, and record *why* (public repo → free minutes) so it
does not get reopened.

### 3.3 Disk is a non-issue and can be stated as such

Measured fixture sizes: 10M×51 bp = 489 MB total, 10M×150 bp = 1.21 GB total → **1.70 GB** of fixtures.
Largest single-cell output (150 bp) ≈ 1.24 GB, and the model's `rm -f "$out"/*.fq.gz` at the end of
`measure` keeps only one cell's output resident. Peak ≈ **3 GB**. The real runner reports **88,190 MiB
free** in `RUNNER_TEMP`.

So the `>= 8 GiB` assertion passes with ~29x headroom and will never fire — it is decorative here.
Keep it anyway (it is two lines and it documents the stalled-writer reasoning), but the plan should not
carry disk as a risk, and reviewer emphasis on "fixtures filling the runner disk" can be retired with
this number.

### 3.4 Cross-job contention does not exist

GitHub-hosted jobs each get a dedicated VM, so `plain-memory` and `clumpify-memory` cannot contend for
page cache or CPU even when they run in the same workflow. The plan's "Measure serially — concurrent
runs on one runner contend for page cache" is correct *within* the job, and the model already satisfies
it with sequential `measure` calls. No change needed; just do not let the concern grow into a `needs:`
edge, which would put the job on the critical path for no benefit.

---

## 4. Alternatives worth considering

### 4.1 The between-cell ordering gate — the one substantive design addition I would make

The plan concludes no platform-independent gate exists. A *relative* one does, and it is free: the
cells are measured in the same job, on the same runner, with the same binary, so the ratios between
them are independent of host and allocator in a way each absolute number is not.

Measured locally (macOS/arm64, mimalloc, 1M pairs, 51 bp):

```
c1  9,420,800 B   (9.0 MiB)
c2 51,740,672 B  (49.3 MiB)   c2/c1 = 5.5x
c4 83,820,544 B  (79.9 MiB)   c4/c2 = 1.62x   c4/c1 = 8.9x
```

Against ±3-5% observed runner noise, the tightest gap has **62% headroom**. So

```
peak(c1) < peak(c2) < peak(c4)
```

is armable on run #1, needs no baseline file, needs no double-digit sample, and — unlike the banner —
is computed from the measured quantity, so a pool that silently collapsed to serial while `--cores 4`
was parsed cannot satisfy it. That is precisely the failure §1.3 shows the banner cannot see.

Two caveats to state honestly:
- It is a *structural* check, not a drift check. A regression that scales every cell by the same factor
  keeps the ratios flat, and only the absolute band catches it. Complement, not replacement.
- I would arm `c1 < c2` and `c1 < c4` immediately (order-of-magnitude gaps) and merely *record*
  `c2 < c4` until the bootstrap data shows the 1.62x gap survives on Linux under CPU contention.

Recording `c4/c2` and `c4/c1` per run during bootstrap costs one `awk` line and gives the banding PR a
second, drift-robust candidate metric for free.

### 4.2 A coarse day-one tripwire

Per §1.6: `±50%` of the first recorded reading, armed immediately, replaced by the fine band once the
sample supports it. Satisfies the noise side of the banding rule by 12-25x; catches the 1.9x class of
event; explicitly does not catch the +24% class. Cheap and honest.

### 4.3 A greppable machine line per cell

The model prints human-formatted `%-14s peak %5d MiB …`. The plan already fixes the precision half of
this by printing KiB as well as MiB — good, and necessary, since the baseline is absolute. Add a
single-token prefix, e.g.

```
PLAINMEM<TAB>plain_c4<TAB>81234
```

so accumulating N readings is `gh run list … | xargs -n1 gh run view --log | grep ^PLAINMEM`. The
clumpify bootstrap needed 31 job logs scraped by hand; the plain bootstrap will need at least as many,
and §1.7 means some will have to be discarded by sha. Making the scrape a one-liner is the cheapest
thing in this review.

### 4.4 Considered and rejected

- **Print a plain-path prediction from the binary** (giving a true absolute promise gate). Rejected:
  the plan fixes "no source change", and §1.1 shows the analytic term is ~1/6 of the real peak, so the
  printed number would not be a promise anyone could gate on.
- **Hash the allocator-relevant lines into the baseline file** (the plan's open question 4). Agreed
  with the plan's instinct: header note plus reviewer discipline, matching the clumpify precedent. But
  note §1.7 — the *judgement* the clumpify file's `git diff -- src Cargo.lock` proxy forces will fire
  on nearly every PR here, so the header should say what to do when it fires rather than only that it
  should be considered.
- **A single-end cell.** The plan measures only `--paired`. SE is a distinct entry point
  (`run_single_end`) and is what a lot of users run. I would not add a cell — a shared-pool regression
  shows in PE too, and both motivating events were paired — but the baseline file should name SE as a
  known blind spot, the way the clumpify file names its blind spots rather than leaving them implicit.
- **25 bp cell** (open question 2). Agree with the plan: skip it. There is no known plain-path read-length
  sensitivity below 51 bp, and 150 bp is the axis tied to #1903.

---

## 5. Validation sufficiency

The eight-row table is good, and better than the plan needed to be: Validations 6 and 8 are exactly
the two controls a reviewer would otherwise have to ask for. Specific notes:

- **V1, V3, V4** — sound, and each mutation points *against* its assertion (the failure mode from the
  harness-false-negative log where a mutation and its expectation pointed the same way and could only
  pass).
- **V2** — sound as a mutation test. Extend it: also assert that `plain_c1` *fails* if the presence-form
  check is applied to it (i.e. the §1.2 fix is itself tested), because that is the one defect that
  would otherwise be discovered by a red on `dev`.
- **V5** — matches the model's behaviour (the `want` gate is skipped when `bootstrap=1`).
- **V6** — right control, and see §2: I expect it to pass easily and to show 10M is far more than
  needed. Record the numbers either way.
- **V8** — the mechanism is unstated. The clumpify precedent used a deliberate throwaway commit; its
  baseline file notes the pool excluded "the two runs carrying Validation 8's deliberate regression".
  Say that here too, including the exclusion, or the verification runs will silently pollute the
  bootstrap pool the plan depends on.

Missing controls:

- **Nothing verifies the cells differ from one another** — §4.1 covers this and doubles as V2's
  replacement.
- **Nothing verifies a memory regression is caught during BOOTSTRAP**, because by construction none is
  — §1.6. If the tripwire in §4.2 is adopted, add: "commit a deliberate +2x allocation with BOOTSTRAP
  set → red on the tripwire, not on drift."
- **Nothing pins the adapter** — §1.8.

---

## 6. Action items

### Critical

1. **Fix the banner assertion per cell.** `--cores 1` prints no banner (`src/main.rs:858`,
   `if cli.cores > 1`), so `plain_c1` fails on every run including under BOOTSTRAP. Assert presence of
   `Using $cores worker threads` for `cores >= 2`; assert *absence* of `Using ` for `plain_c1`. §1.2
2. **Stop treating the banner as the `^clumpify:` analogue.** It is printed from the parsed flag, before
   any pool exists, and `src/parallel.rs` prints nothing — so it cannot detect a pool that silently
   collapsed. Adopt the ordering gate of §4.1 as the real engagement check, and reword Behavior 5 and
   Self-Review accordingly. §1.3, §4.1
3. **Specify where the drift width lives and in what unit.** The model's single global width works only
   because its quantity is a ratio; plain cells span 9 MiB → several hundred MiB, and the two evidenced
   regression deltas are absolute. Commit to a three-column TSV
   (`cell<TAB>baseline_kb<TAB>width_kb`), per-cell absolute widths. §1.5

### Important

4. **Add a coarse ±50% tripwire armed on run #1**, so the job is not inert for the ~2 weeks the
   bootstrap takes. It clears the banding rule's noise side by 12-25x on the observed ±3-5% dispersion,
   catches the mimalloc-class 1.9x event, and does not pretend to catch the +24% class. §1.6, §4.2
5. **Define the observation pool.** The plain path is the shared path, so most PRs invalidate readings —
   a much sharper problem than the clumpify job's. Pool `push`-on-`dev` only, record the sha per
   reading, exclude V8's deliberate-regression runs. §1.7, §5
6. **Assert the detected adapter per cell** (`Adapter: Nextera`). The detection margins are 5-vs-2
   (51 bp) and 10-vs-6 (150 bp) on random ACGT; a heuristic or fixture change can flip it and would
   present as unexplained baseline drift. Keep auto-detection rather than pinning `--nextera`, so the
   1M-read scan stays inside what is measured. §1.8
7. **Add `if: github.event_name != 'schedule'`.** Every non-audit job in this workflow carries it; the
   weekly cron exists for `audit` only. §1.10
8. **Name the two idioms that must not be copied** — the clumpify prediction-margin assertion, and the
   three-field `case` empty-guard — plus a caution that `measure` must never be called through a pipe
   (`fail=1` would be lost in the subshell). §1.9

### Optional

9. **Delete the "drop `plain_c1` and `plain_c2`" fallback** and replace the Efficiency section's cost
   argument with measured numbers: the model job runs seven cells in 5 min (build 73 s / fixtures 30 s /
   measure 214 s); the plain job projects to ~7 min warm, ~12 min cold, against a 30-min timeout. Also
   correct "four cells" — it is a stale `ci.yml` comment. §3.1
10. **Close open question 1 in favour of a separate job**, recording that the repo is public so the
    ~1.5 min saved by sharing has no budget to come out of. §3.2
11. **Retire disk as a risk** with the measured figures: 1.70 GB fixtures + ≤1.24 GB live output against
    88,190 MiB free. Keep the assertion, drop the worry. §3.3
12. **Reframe open question 3 with the fact that neither trunk is protected** — no check in this repo
    mechanically blocks a merge, so "gating" here means a visible red plus reviewer discipline. That
    lowers the flaky-red liability and raises the value of items 4-6. §2
13. **Print a greppable `PLAINMEM<TAB>label<TAB>peak_kb` line** to make the bootstrap scrape a
    one-liner. §4.3
14. **Record `c4/c2` and `c4/c1` during bootstrap** so the banding PR has a drift-robust candidate
    metric beside the absolute KiB. §4.1
15. **Name single-end as a known blind spot** in the baseline header, and state the random-ACGT
    incompressibility caveat there as the plan already intends. §4.4
16. **Note in the header that `--cores 4` is oversubscribed on a 4-vCPU runner** and that the model's
    `1G_c4` shows this is tolerable (±3.5% over 27 observations). §2

---

## 7. Summary

The plan is well-grounded: I verified every commit it cites, every idiom it claims to reuse, the
banding-rule quote, the FastQC-cell precedent, the `--gzip`/compression mechanics, the `~9 MiB` c1
figure, and its central premise that the plain path prints no prediction. The motivation is real and
the gap it names is real.

Three things stand between it and a job that works. `plain_c1` cannot satisfy the banner assertion and
would be red from the first merge. The width has no home and no unit, and the model's schema does not
carry over because absolute KiB spans two orders of magnitude where ratios spanned none. And the job as
designed guards nothing about memory for the ~2 weeks the bootstrap runs — during which the very class
of event that motivated it would ship unseen, exactly as `2ab0918` did.

The last one is the most recoverable, and the plan is a little too quick to accept it. "No printed
prediction" removes the *absolute* gate, not every platform-independent gate: the cells are measured
side by side in one job, and `c1 < c2 < c4` holds with 62% headroom at its tightest against ±3-5%
noise. That ordering, plus a coarse ±50% tripwire, arms the job on run #1 under arithmetic that
satisfies the plan's own banding rule rather than waiving it.

Cost, disk and cross-job contention are all non-issues on measured numbers, and the plan is more
conservative about them than it needs to be — the five-cell design should ship whole.
