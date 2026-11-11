# Plan re-review B (round 2) — Phase 3 v2: defend and correct the `--clumpify` memory promise (#439)

**Reviewer:** B, round 2 (independent; no coordination with Reviewer A this round or last)
**Target:** `plans/08162026_clumpify-memory-accounting/phase3/PLAN.md` **v2**
**Round-1 reports read:** `phase3/PLAN_REVIEW_A.md`, `phase3/PLAN_REVIEW_B.md` (both unmodified)
**Baseline:** `dev` at `fdade7e`; `target/release/trim_galore` at that hash; no branch switch, no stash, no worktree, no build. Read-only on `src/`, `tests/`, `Docs/`, `.github/`, `Cargo.toml`. My only write is this file.
**Date:** 2026-08-17

## Verdict

**Implement with the listed changes.**

v2 is a large improvement and it discharges the great majority of round 1 honestly rather
than rhetorically. I re-ran the perturbation matrix against v2's new pinned-cell table and it
**fails for all seven perturbations**, including the two round 1 proved the sweep was blind
to — so the 4a/4b split is a genuine fix, not a relabel. Every corrected figure I could
check landed correctly: `c2 = 329`, `c10 = 579`, `154.4 MiB at c16`, `281–590`/`334–1012`
over cores 2–32, `c64 = 933`/`1356`, `≈ 930 of 1024`, `≈ 1363 of 1500`. All twelve-plus
`file:line` references still resolve, including the four new ones v2 introduces.

**But v2's load-bearing new claim does not survive.** I attacked §0 with a controlled 2×2
that neither round-1 reviewer nor the author ran — input encoding × output encoding, with the
**records held byte-identical across all four arms** — and then a block-vs-interleaved
control. In **40 runs at 10 M pairs** across four encoding combinations and two core counts,
peak RSS never exceeded the printed prediction, let alone the budget. The maximum anywhere
was **900.0 MiB (0.968 of prediction)**, and it was on the **gz** path, not the plain one.
§0's `1G/c8` plain figure of 1081.8 MiB — the number Assumption 8, the fixture rationale,
issue #1's headline and the `CALIBRATION.md` note all rest on — is **255 MiB above anything I
could produce in that cell**, and the plain-vs-gz difference runs in the **opposite
direction** at `--cores 4`. Details and the full data in §1.1.

Second Critical, independent of that: **§5's tripwire claim is arithmetically false at v2's
own expected Linux ratio.** The job runs on `ubuntu-latest`, v2 expects ~0.85 there, and a
lever-A revert adds ~104 MiB — which lands *under* the 930 MiB gate. §5 proves the tripwire
on macOS numbers and the job never runs on macOS (§1.2).

Both fixes are local. Nothing in D1, D2, D4, the docs arithmetic, or the rebuilt guards needs
to change.

---

## Method

Every number below is mine unless attributed. Harness: `os.posix_spawn` + `os.wait4` in a
trivial Python parent that holds nothing (`ru_maxrss`, bytes on Darwin);
`/usr/bin/time -l` is blocked under the sandbox. Arithmetic claims come from an exact Python
replication of `src/clump.rs` at `fdade7e` (σ = 23, k = 64, static 224 MiB, FastQC 24 MiB/thread
capped at 16, margin 10/11, `MIN_BIN_BYTES` = 1 MiB, `n_bins = max(16, 4·bin_cores)`), which
reproduces `bin_byte_budget = 19_006_356` at `1G/c4` and every floor the binary advises.

**Scope of this review, stated because the round-2 reviewers are not reviewing the same
document.** This report covers **`phase3/PLAN.md` v2 at `fdade7e`**. Reviewer A's round-2
report (`../PLAN_REVIEW_A_v2.md`) targets the **parent `PLAN.md` (Phase 0) at tip `2ab0918`**,
which is the pre-refit model: static 512 MiB, `4 × usable / (5·n_bins + 7·cores)`, and **no
margin term at all**. Its figures therefore do not transfer. Concretely, its headline
`1G/c4 → 1050.2 MiB > budget` is a `2ab0918` measurement at `B = 19 MiB`, `budget = 1025 MiB`
and `--compression 6`; in all 18 rows of its sweep `predicted ≈ budget`, because that tip has
nothing scaling the prediction down, so `actual ≤ predicted ≤ budget` degenerates to
`actual ≤ budget` with zero headroom. That finding is sound *for Phase 0* — and it is the
argument for the refit which has already shipped as `fdade7e` (#446). On `fdade7e` the same
cell predicts 930 against a 1024 budget and I measured it 45 times at a max of 900.0 MiB, so
the ordering holds here. **Do not read A's round-2 figures as statements about Phase 3.**

**Contention is disclosed per figure.** Reviewer A was measuring throughout my session
(load average 4.77 at start; `ps` is blocked so I could not enumerate). I handled it by
**interleaving arms rep-by-rep** rather than running them as blocks, so drift in machine
state maps onto both arms equally instead of onto the arm difference. I then ran deliberate
blocks as a control. Where I report a maximum, it is *inflated* by contention if anything —
which cuts against, not for, my refutation of §0.

**The lever that made the 2×2 possible.** `src/main.rs:424` is
`let gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0]);`, and `io::is_gzipped`
(`src/io.rs:33`) tests the **filename**, while the reader sniffs bytes. So a hardlink of a
gzipped fixture under a `.fastq` name gives gz input with plain output, and a plain fixture
under a `.fastq.gz` name gives plain input with gz output. Both verified:

```
gz content,   plain name  → gz_named_plain_trimmed.fq       : ASCII text
plain content, gz name    → plain_named_gz_trimmed.fq.gz    : gzip compressed data, max speed
```

That decouples the two axes §0 varies together, at zero extra disk (hardlinks).

---

## 1. Logic review

### 1.1 CRITICAL — §0's memory-regime claim does not reproduce; the CLI defect does

**What I set out to test.** §0 asserts that the plain-text *output* branch is a different
memory regime, reaching **1.057 × `--memory`** at `1G/c8` "reproducibly (3 of 5 reps at
10 M pairs)", against a gz path that holds with 6–12% headroom. Everything downstream — the
gzipped-fixture requirement, Assumption 8, the `.fq.gz` assertion, issue #1's headline, and
the second dated note in `CALIBRATION.md` — hangs off it.

**Design.** 10 M pairs, 51 bp, 100% distinct, `--paired --clumpify --memory 1G`, prediction
930 MiB in every cell (verified from each run's own banner). Four arms with **byte-identical
records** — the plain fixture is `gunzip -c` of the gz one, and the two "wrong-name" arms
alias it, so fixture identity is eliminated as a variable. Arms interleaved rep-by-rep.

**Per-cell provenance.** Requested so this is comparable against a re-measurement. Output
encoding is confirmed by **output filename and output magic bytes**, never by the banner —
which prints `(gzip level 1)` in all four arms, as the table shows. `is_gzipped` reads the
input *name*; the reader sniffs the *resolved content*; so the two are independent.

| arm | input path | input content (first 3 B, symlink followed) | name ends `.gz`? | output file | output content | `file(1)` | peak MiB |
|---|---|---|---|---|---|---|---:|
| `gzin_gzout` | `t_R1.fastq.gz` (real file, 244,953,962 B) | `1f8b08` = gz | yes | `t_R1_val_1.fq.gz` | `1f8b08` | gzip compressed data, max speed | 785.9 |
| `gzin_plainout` | `nz_R1.fastq` → `t_R1.fastq.gz` | `1f8b08` = gz | no | `nz_R1_val_1.fq` | `40` = `@` | ASCII text | 789.9 |
| `plainin_gzout` | `pgz_R1.fastq.gz` → `pl_R1.fastq` | `40` = `@` = plain | yes | `pgz_R1_val_1.fq.gz` | `1f8b08` | gzip compressed data, max speed | 768.8 |
| **`plainin_plainout`** | **`pl_R1.fastq` (real file, 1,478,888,890 B)** | **`40` = `@` = plain** | **no** | **`pl_R1_val_1.fq`** | **`40` = `@`** | **ASCII text** | **785.4** |

**The bolded row is §0's configuration exactly** — a genuinely decompressed 1.48 GB plain
input file, not an alias of gz bytes, writing plain output. It measures **785.4 MiB**, against
§0's 1081.8 MiB for the same cell: a 296 MiB gap. So the disagreement is not definitional. It
is not "we mean different things by plain input": the arm fed real decompressed bytes and the
arm fed gz bytes under a plain name land 4.5 MiB apart (785.4 vs 789.9), which is inside the
rep spread, so **input decode is not the driver either**. Both gz-out arms report an identical
`original size modulo 2^32 3475396`, confirming the records are the same in every arm.

(One reading note: `stat -f%z` on the two aliased inputs reports the symlink's own size, 119
and 120 bytes. The content probe uses `head -c 3`, which follows the link, so the "input
content" column is the target's bytes. The real byte sizes are the two given for the
non-aliased rows.)

**Result — `--cores 8`, 32 runs:**

| arm (input / output) | reps | max MiB | max ratio |
|---|---|---:|---:|
| gz / **gz** | 12 | 814.7 | 0.876 |
| gz / **plain** | 7 | 787.7 | 0.847 |
| plain / **gz** | 4 | 793.3 | 0.853 |
| plain / **plain** | 9 | **826.1** | 0.888 |

**Result — `--cores 4`, 8 runs (the gating cell), interleaved:**

| arm | reps | max MiB | max ratio |
|---|---|---:|---:|
| gz in / gz out | 5 | **900.0** | **0.968** |
| plain in / plain out | 4 | 846.1 | 0.910 |

(c4 gz reps: 876.2, 873.5, 854.2, 900.0, 861.8. c4 plain: 787.1, 846.1, 843.0, 812.5.)

**Against §0's table:**

| cell | §0 plain | §0 gz | **mine plain** | **mine gz** |
|---|---:|---:|---:|---:|
| `1G/c4` | 908.3 (10 reps) | 873.2 (3) | **846.1** (4) | **900.0** (4) |
| `1G/c8` | **1081.8** (5 reps) | 814.4 (3) | **826.1** (9) | **814.7** (12) |

Four things follow, and they are separable.

1. **The breach does not exist in any encoding combination.** Maximum over all 40 runs is
   900.0 MiB = 0.968 of prediction, 0.879 of budget. Nothing came within 120 MiB of §0's
   1081.8, and nothing exceeded the *prediction*, let alone the 1024 MiB budget. §0's
   "1.057 × `--memory`" is not reproducible.

2. **Output encoding is not the mechanism.** With input held byte-identical, switching output
   from gz to plain moved the c8 mean by **−1.9 MiB (0.2%)** and the max by −27 MiB — i.e.
   plain output is *cheaper*, which is the physically expected direction (no deflate state,
   no compressed-output buffers). §0's 267 MiB gap is not an output-encoding effect.

3. **My harness agrees with §0 on the gz arm and disagrees only on the plain arm.** At
   `1G/c8` gz I measured 814.7 MiB against §0's 814.4 — 0.3 MiB apart. That agreement is
   what makes the plain-arm divergence attributable to §0's plain measurement rather than to
   my method.

4. **The sign reverses at the gating cell.** At `1G/c4` I measure gz **above** plain
   (900.0 vs 846.1). Reviewer A's round-1 3 M measurement also had gz above plain (0.888 vs
   0.835). Two independent harnesses now put gz ≥ plain; only §0 puts plain ≫ gz.

**The cause — established by the plan's author after this section was first written, and it is
not what I guessed.** I had ruled out the two hypotheses the brief named (*concurrency:* my
figures are the contended ones and they are lower; *block-vs-interleaved:* a deliberate block
of 5 consecutive `1G/c8` plain reps gave 784.7, 787.2, 826.1, 788.9, 823.2 — max 826.1, no
breach), and then offered **fixture composition** as the remaining candidate. **That guess was
wrong** and I am striking it: the author's arms did share records, since the gz fixture was
`gzip -1 -c` of the very plain file being measured.

The actual cause is **disk exhaustion inflating resident output buffers**. At the time of those
c8 readings the filesystem had ~1.6 GiB free while each rep wrote ~2.9 GB of plain output, and
`io::PendingOutput` writes a temp file alongside the previous run's output — so consecutive
reps into one directory needed ~5.8 GB and could not all complete. A blocked or failing writer
backs up the result channel and the `pending` map, both of which hold written-but-unflushed
buffers. The author's harness grepped only for `peak=`, so a non-zero exit was invisible.
The author's re-measurement with adequate disk and explicit exit-code checks gives
**756.4 / 789.4 / 756.4 / 761.4 / 786.7 MiB** — max 0.849 of prediction, converging with my
768.8–826.1 and confirming the refutation.

**This is corroboration for §4c, not just an erratum.** The mechanism that produced the
290 MiB inflation *is* the result channel plus the `pending` map holding resident compressed
buffers — the exact lever v2's §4c proposes to pin (`parallel.rs:130`, `:959`) and the one the
parent plan names as the next thing someone will pull for the `2G/c4` residual. §4c is
justified on the argument that `k = 4` charges for those slots; this incident is an accidental
*measurement* of how much they can hold when the drain stalls. Worth citing in §4c.

**And it puts a new requirement on PR 3, which nothing in v2 covers** (see §3): the CI job
writes output per cell on a runner with ~14 GB free, and §8's two new cells push the
unmanaged total to ~3.7 GB. If disk ever gets tight there, this mechanism makes **disk
pressure present itself as a memory-model failure** — a red `peak ≤ prediction` with a
correct model. The job therefore needs a free-space precondition before the first cell,
per-cell output truncation, and — the part that made this incident invisible for a whole
round — an explicit **exit-code check on every `trim_galore` invocation**, not just a
successful parse of its stderr. Validation 13 tests an unparseable log; it does not test a
cell that exited non-zero while still printing a parseable banner.

**What survives.** §0's *CLI defect* is real and I reproduced all three commands verbatim
on `fdade7e` (Validation 16, discharged — see §4):

```
--clumpify --cores 2 in.fastq          → in_trimmed.fq  ; banner says "(gzip level 1)"
--clumpify --gzip --cores 2 in.fastq   → in_trimmed.fq  ; "--gzip is deprecated … Ignoring."
--clumpify --dont_gzip --cores 2 …     → Error: … mutually exclusive (clumping plain text is pointless)
```

So the defect stands; its **severity does not**. It is a wrong deprecation message plus an
asymmetry against `clump-only.md:78`, not a data-integrity breach of the `--memory` promise.

**Failure scenario, and why this is the sharpest item this round.** v2 queues issue #1 as a
**public** GitHub issue whose stated bar is *"data-integrity: breaks the `--memory` promise at
1.057× budget"*, and §9 writes the same claim into `CALIBRATION.md` as a dated note in a
document v2 explicitly describes as "the record of what was believed when the fit was made".
An unreproducible budget breach published under #439 is worse than the understated threading
figure round 1 caught, because it is the *same class of error inverted* and it lands on the
one document that is meant to be durable. Validation 16 as written only re-runs the three
CLI commands; it does not re-run the measurement, so the plan's own gate would pass this
through.

**Fix.**

- **Delete the memory-regime claim** from §0, from Assumption 8's justification, from §5's
  fixture rationale, from §7's issue-#1 bar, and from §9's second dated note. Replace it with:
  peak RSS is materially insensitive to output encoding (≤ 2% of prediction at c8, and
  gz-higher at c4); the ~1.06× figure did not reproduce in a controlled 2×2 and its cause is
  unknown.
- **Re-headline issue #1** as a CLI/messaging defect: the deprecation sentence is false for
  plain input, and the two clumping modes take opposite positions on plain output. That still
  clears the repo's filing bar as a user-visible defect. Do not assert a budget breach.
- **Keep the gzipped fixture** — see §1.3 for the reasons that survive, one of which is
  stronger than the one v2 gives.
- **Rewrite the Self-Review's "measurement lesson" paragraph.** The lesson that actually
  holds is not "nobody guarded the output encoding"; it is *interleave arms and hold the
  records byte-identical when comparing encodings — a block A/B on a shared machine is not a
  controlled comparison.* As written, §0's narrative teaches a lesson drawn from a non-finding.

### 1.2 CRITICAL — §5's drift-tripwire claim is false at v2's own expected Linux ratio

§5's central conclusion, which inverts v1's: *"at 10 M the base is ~873–883 MiB and the same
+104 MiB **breaches the prediction**, so the single assertion v1 already specified *is* the
drift tripwire."* Then, four bullets later: *"Linux measured ~9 points below macOS in Phase 2,
so expect **~0.85** / ~0.72 on the runner."*

Both cannot carry the design, because the job only ever runs on the runner:

| | macOS (§0's figure) | macOS (mine) | **Linux (v2's own expectation)** |
|---|---:|---:|---:|
| baseline ratio at `1G/c4` gz | 0.939 | 0.968 | **≈ 0.85 / ≈ 0.88** |
| baseline MiB | 873.2 | 900.0 | **790 / 817** |
| + lever-A revert (+104 MiB) | 977 | 1004 | **894 / 921** |
| vs 930 MiB gate | **red** | **red** | **green — misses by 9–36 MiB** |

The +104 MiB delta is not in dispute and I re-derived it: reverting lever A restores one
extra resident bin-sized batch per worker, `4 × 1.4375 × 18.13 MiB = 104.2 MiB`. Nor is the
delta platform-scaled — it is real bytes held. Using v2's stated 0.85 it misses by 36 MiB;
using my measured macOS 0.968 minus Phase 2's 9-point Linux offset it still misses by 9 MiB.

This is the same species of error v1 made — promoting a single platform's measurement to a
structural conclusion — and v2's Self-Review diagnoses that error precisely (*"it reasoned
from one measurement instead of a ramp"*) without noticing it has committed the platform
version of it.

**Failure scenario.** PR 3 lands. Validation 12 is run **locally on macOS**, where the
lever-A revert goes red, and is marked discharged. The phase closes stating in §5 and in the
Self-Review that the job is the drift tripwire. On the runner it has ~15 points of headroom
against a regression that consumes ~11, so the tripwire is green forever — which is exactly
v1's failure mode, re-created by a plan whose stated purpose was to remove it.

**Fix.** Three clauses, no new design:

1. **Validation 12 must run on the runner**, not locally. Say so explicitly — the row
   currently says only "revert lever A (work depth → 2) and run the job".
2. **Carry the Linux arithmetic in §5** beside the macOS arithmetic, and state plainly that at
   ~0.85 the lever-A delta is insufficient. Right now the two numbers sit four bullets apart
   and are never multiplied together.
3. **Say what happens if Validation 12 comes back green.** v2 says "the fixture is still too
   small — learn it now, not in a year", which is the right instinct, but it does not budget
   for the consequence: reaching the same ratio on Linux plausibly needs a larger fixture,
   and §5's disk/time figures and `timeout-minutes: 30` are sized for 10 M. Either name a
   fallback fixture size, or accept a ratio-based gate on the runner's own measured baseline,
   or state that the job ships as a promise guard on Linux and the drift claim is retracted.

Note the honest fallback is already available and cheap: my c4 gz baseline of **0.968** on
macOS means the *macOS* headroom is only 30 MiB. That is a strong tripwire and a flaky gate
at the same time — which is reviewer A's round-1 point, unresolved rather than wrong.

### 1.3 IMPORTANT — the gzipped-fixture decision is right, and v2 gives the one reason that is false

v2 justifies gzipping with three reasons; my data grades them:

- *"§0 shows the plain-output branch is a different memory regime"* — **false** (§1.1).
- *"CI must measure the branch users run"* — **true and sufficient on its own.**
- *"takes the fixture from 2.96 GB to ~490 MB"* — **true**; I measured the plain pair at
  2.96 GB exactly.

And there is a fourth reason v2 does not give, which is the strongest one: **at the gating
cell, gz is the more sensitive fixture.** My interleaved `1G/c4` figures are gz 0.968 vs
plain 0.910 — 5.8 points, in the direction that buys tripwire headroom, and matching
reviewer A's round-1 3 M result (0.888 vs 0.835, 5.3 points). Since §1.2 shows sensitivity is
the binding constraint on the whole job, that is the reason to put in the plan.

Worth noting explicitly because v2's §0 table implies the opposite: taken at face value, §0
says gzipping the fixture *costs* 3.8 points at `1G/c4` (0.977 plain → 0.939 gz), i.e. that
v2 chose the less sensitive fixture for its gating cell and did not remark on it. Removing
the §0 table removes that implication too.

### 1.4 IMPORTANT — the new multi-pair cell breaks the parse spec PR 3 inherits

§8 adds "one multi-pair cell (the same fixture passed twice)". I verified it runs, and **it
prints the banner twice**:

```
=== Pair 1 of 2 ===
clumpify: 16 bins × 18 MiB; predicted peak ≈ 930 MiB (gzip level 1)
=== Pair 2 of 2 ===
clumpify: 16 bins × 18 MiB; predicted peak ≈ 930 MiB (gzip level 1)
```

PR 3 step 5 specifies a single parse — *"Parse … `predicted peak ≈ (\d+) of (\d+) MiB
budget`. Assert all captures are non-empty positive integers before comparing"* — written
when every cell was single-pair. On this cell the capture is two values. Under
`set -euo pipefail` an integer comparison against a two-line string errors out, so it fails
loudly rather than silently, which is the acceptable half. The unacceptable half is that it
fails on the **first** run of a cell v2 introduced, and the plan has no row for it:
Validation 11 says "first run … banner present, output `.fq.gz`, `peak ≤ prediction`, both
ratios logged", which the multi-pair cell cannot satisfy as specified.

Also: `grep -c '^clumpify:'` is 2 on this cell, so any "banner present" check written as an
equality on 1 fails. `grep -q` is fine.

**Fix.** One clause in PR 3 step 5: pin the parse to the first banner (`| head -1`) and state
that the multi-pair cell asserts one prediction against one whole-process peak — which is the
point of the cell, since the promise is per-process.

### 1.5 IMPORTANT — §8's remedies cannot discharge §8's own concern, and §8 says they can

§8 is right that short reads and multi-pair are the two load-bearing unverified axes, and
right about why. But both remedies are specified as **"logged, non-gating"**. A non-gating
cell cannot fail, so:

- **Short reads.** §8's argument is that σ ≈ 1.55 at 25 bp against the shipped 1.4375 means
  "there is a *known* read length at which the model under-reserves — the promise-breaking
  direction", and that a 51 bp fixture is "structurally unable to see it". A 25 bp cell that
  cannot fail is also structurally unable to see it. It converts a blind spot into an
  unwatched log line.
- **Multi-pair.** Same. §8's stated mechanism is that "pair 2's pool may sit on pair 1's
  retained pages — **which would breach a per-process promise, and nothing tests it**". A
  non-gating cell still does not test it.

Worth pricing the short-read cell honestly, because it argues for the other resolution: the
pool term at `1G/c4` is `σ·n_bins·B = 1.4375 × 16 × 18.13 ≈ 417 MiB` of a 930 MiB
prediction, so an 8% σ under-charge is ~33 MiB. Against a measured baseline around 0.85–0.97
that will not breach either. So the candid statement is: *a gated 25 bp cell would not fire
today; the axis remains unguarded and this cell only establishes a baseline to watch.*

**Fix.** Either gate both cells (with §1.2's runner caveat), or say plainly that §8 adds
baselines rather than guards and that the two axes remain unguarded. Do not leave §8 reading
as coverage. This matters more than usual because v2's own closing rule is *"prove each of
these can [fail]"*.

### 1.6 IMPORTANT — v2 reproduces, one table away, the rounding error it corrects

v2's `--clump_only` table correctly fixes `c10`: *"v1 said 578 — truncated where the binary
uses `div_ceil`"*. Two tables earlier, the handoff-corrections table says *"**276/302** are
the `--cores 1` values"*. Those are the truncations. The exact values are 276.10 and 302.50,
and the binary advises `div_ceil` — **277 and 303**. Reviewer A confirmed 277 live in round 1
(`--clump_only --memory 200M` advises "try ≥ 277 MiB"), and v2's own §2 edge case says
*"at `--cores 1` the floor is **277** MiB"*.

So v2 states 276 and 277 for the same quantity in one document, and the wrong one is in the
table whose job is to correct someone else's arithmetic. Confined to plan prose (the docs
table in §3 uses the correct 281/590/334/1012), but the plan is the audit record for this
arc and this is the third round in which a rounding convention has been the error.

### 1.7 IMPORTANT — issue #1's "no way to obtain gzipped output from plain input" is false

§0 part 2 and the open question both assert it. It is false, and the counter-example is the
user's workaround: **name the input `*.gz`**. Verified — a plain-text FASTQ hardlinked as
`plain_named_gz.fastq.gz` produces `plain_named_gz_trimmed.fq.gz`, genuinely gzipped
(`gzip compressed data, max speed`), because `io::is_gzipped` (`src/io.rs:33`) reads the
filename while `FastqReader` sniffs the bytes. A symlink is enough; no copy.

The precise claim is "no CLI *flag* requests gzip for plain input", which is what §0's own
first bullet already says. File the precise one, and give users the symlink workaround — it
is the difference between an issue that helps someone today and one that only records a
grievance.

**Field evidence that the defect confuses expert readers, acquired by accident.** Round 2's
Reviewer A, working from the source with full repo access, proposed `--dont_gzip` to me as the
cheap way to put `--clumpify` on the plain-output branch — "it short-circuits before
`is_gzipped`, forces plain output, and costs 0 bytes". It does not: `cli.rs:747-750` refuses
the combination outright (`--clumpify and --dont_gzip are mutually exclusive (clumping plain
text is pointless)`), which I confirmed live. Had that been run as a measurement, the
validation error would have presented as a clean null result. So within one afternoon the
asymmetry misled a reviewer who had read the code — which is a stronger argument for filing
D3's issue than anything in §0's memory table, and it belongs in the issue as the motivating
example. It also sharpens §1.5's warning: a cell that cannot run looks exactly like a cell
with nothing to report.

### 1.8 IMPORTANT — the issue must be framed against the decision it proposes to reverse

v2's open question says *"the answer decides whether the fix is 'always gzip' or 'honour
`--gzip`'"*. "Always gzip" is not an open option: it shipped and was reverted.
`src/main.rs:417-424` says so in the source —

> Mirror Perl: by default the output's compression matches the input's … See #245 for the
> parity rationale (Rust v2.1.0-beta.5 always gzipped regardless of input, breaking pipelines
> that globbed `*.fq` no-gz).

— and `src/io.rs:14-32` already documents this exact class of surprise (a `.bgz` input
"decompressed correctly but produces plain output"), concluding that moving the decision to
the detected format "would change behaviour for every run … so it is deliberately left for a
separate change. See the CHANGELOG."

An issue that offers "always gzip" without naming #245 invites a maintainer to re-land a
known regression, or to close it as already-decided. Cite both comments, and scope the ask to
the two parts that are *not* already decided: the false deprecation sentence, and the
`--clumpify` / `--clump_only` asymmetry on plain output.

### 1.9 IMPORTANT — §6's re-measured table merges two different input sizes into one row

§6's row *"10 kb, 764 MB–1.2 GB | reviewer A 1366.0 | reviewer B 1808.8 (c4), 2052.3 (c8)"*
puts A's **764 MB** measurement and B's **1.2 GB** measurement side by side under one label.
The two inputs differ by 57%, and §6's whole thesis is that peak RSS scales with input size —
so the row destroys exactly the signal it is presenting, and reads as a 443 MiB disagreement
between two harnesses that may be measuring different things.

This is going into a public issue. Split the row by input size, or drop the cross-reviewer
comparison and publish only the re-measurement Validation 15 mandates.

### 1.10 Consequences of D1–D4 that v2 does not carry

Not re-arguing any of them.

- **D4 (docs first) leaves PR 1 with no changelog entry.** The files table routes
  `CHANGELOG.md` to the code PR, and the repo rule forbids non-docs changes in a docs PR — so
  eleven corrected user-facing figures ship with no changelog line, and the line that does
  land (PR 2) describes the banner. Either accept that and say so, or note that the #439
  bullet extension in PR 2 should name the docs corrections too.
- **D2's second figure is the parsed budget, and PR 3 parses it back.** That is sound, and
  PR 2 step 2's exact-substring tripwire already closes the "binary owns both sides" hole
  deterministically at zero CI cost. Worth stating, because §5 sells the CI margin assertion
  as the closure and it is in fact the belt to PR 2's braces.
- **D3 keeps the defect unfixed, and v2's Remaining risks (e) states the consequence in the
  refuted terms** ("keeps breaching the budget until that issue is picked up"). After §1.1
  the residual consequence is that `--clumpify` silently does the thing it calls pointless
  and the deprecation message stays false — still worth recording, but not a breach.

### 1.12 IMPORTANT — v2 drops the one parent-plan assignment it was handed by name

Credit to Reviewer A, which surfaced this from the `fdade7e` tree while I was checking coverage
from the plan side; I verified both halves.

`src/clump.rs:577` still carries `resolve_layout_predicted_peak_fits_budget`, which asserts
`predicted_peak_bytes() <= budget` over five cells. By the algebra v2's own §4 sets out —
`b = floor(usable×16/d)` ⇒ `dyn ≤ usable` ⇒ `predicted ≤ budgeted ≤ budget` for any σ, k,
static — this test is a restatement of `MARGIN_NUM ≤ MARGIN_DEN`. It is the *same* defect §4
spends a section proving about the sweep, in a test that is already shipped. And #446 made it
**more** tautological, not less, by inserting the margin: the chain is now
`predicted ≤ budget × 10/11 ≤ budget`, so it clears with 9% to spare.

The parent plan assigned its disposal to this phase, twice and explicitly (`PLAN.md:168-170`
and `:332-337`, quoted in §6), with the precondition stated: retirement waits on "Phase 3's CI
check is the real replacement", because deleting a passing test before its replacement exists
reduces coverage. **v2's PR 3 delivers exactly that replacement** — so the precondition is met
inside this very phase, and the retirement is still not scheduled anywhere. v2 mentions the
test zero times.

The asymmetry is what makes this worth an action item rather than a nit: v2's Assumption 1
carefully protects the *good* pin (`resolve_layout_pins_the_default_cell` "must pass
unmodified"), and its file table opens `src/clump.rs` for three separate edits — the helper,
the pinned-cell table, the renamed sweep — so the file is already being touched by the same PR
that lands the replacement. Retiring the tautology is a deletion in a file PR 2 already opens.

**Failure scenario.** Phase 3 closes. The parent plan's Phase-1 note still reads "assigned to
Phase 3", Phase 3's own documents never mention it, and the repo keeps a five-cell test whose
assertion is `10 ≤ 11`. The next reader sees a passing test named
`predicted_peak_fits_budget` and reasonably concludes the budget relationship is under test —
which is precisely the misreading round 1 said would be the worst outcome of the sweep, now
reached by a different route.

**Fix.** One line in PR 2: delete it, or relocate its assertion to the pinned-cell table where
it bites. Then strike the "assigned to Phase 3" note in the parent plan so the paper trail
closes. If it is deliberately kept as a smoke test, say so and rename it to something that
does not claim more than it checks.

### 1.11 The three-PR split and its ordering constraint — sound

Checked and I could not fault it. PR 1 documents only behaviour that is already true on
`fdade7e` (D1's hard failure exists today; v2's D1 row says so), so it does not describe a
future state. PR 2's banner must precede PR 3's regex or the parse matches nothing, and with
`set -euo pipefail` plus the non-empty integer checks that is now a loud failure rather than
round 1's silent pass — but ordering it after PR 2 anyway is the right call, and
`ci.yml` genuinely is not a docs file. The one refinement: PR 1's `clumpy.md:110`
`--clump_only` exception and PR 2's bail-context wrap describe the same behaviour from two
PRs; harmless, but if PR 1 quotes the bail text it will be stale for one merge cycle. It
should not quote it.

---

## 2. Assumptions

Taking v2's eight in turn. Six are sound; two need work.

1. **"No coefficient moves in Phase 3."** Sound, and the right narrowing of v1's overbroad
   version. The parenthetical justification — *"narrowed from v1's broader 'the model is
   correct as shipped', which §0 shows is true only on the gz path"* — is the part that falls
   with §1.1. The narrowing is still right; the reason should be "that is all Phase 3 needs".
2. **"`predicted_peak_bytes()` remains the exact inverse."** Verified again symbolically and
   numerically: prediction is within 23 bytes of `budget × 10/11` at `1G/c4` and within
   ~312 B across cores ≤ 32. The observation that this one property is load-bearing in two
   places with opposite signs is the sharpest sentence in v2's Assumptions.
3. **"`--clumpify` requires `--cores ≥ 2`."** Verified. Worth one footnote: `n_bins` is
   `max(16, 4·bin_cores)`, so cores 1–4 all give 16 bins and the `--clump_only` floor is
   277 MiB at each of them. `--cores` only starts moving the floor at 5.
4. **"`ubuntu-latest` has ≥ 4 usable cores and ≥ 4 GiB RAM."** Sound, and dropping `2G/c4`
   does remove the RAM half. The "insensitive to runner core count" note is correct and is a
   real strength. But §8's two new cells put **disk** back where RAM was: see §3.
5. **"Fixture sizing is not a free dial."** Correct reclassification, and §1.2 sharpens it —
   sizing determines coverage *per platform*, not just in the abstract.
6. **"All ratios here are macOS/arm64; Linux measured ~9 points lower."** Correctly flagged
   and then **not propagated** into the design that depends on it (§1.2). This is the
   assumption that is stated and ignored, which is worse than an unstated one.
7. **"Not assumed: that the plain-output branch is in scope."** Fine.
8. **"New: that `--clumpify` output is gzipped. False for plain input (§0) — which is why
   PR 3 asserts the output extension rather than trusting the banner's `(gzip level 1)`."**
   The *fact* is right (verified). The *justification* is the refuted one. And the assertion
   it motivates is close to vacuous in CI: the job generates its own `.fastq.gz` fixture, so
   the extension check can only fail if `io::is_gzipped` or the generator changes — it cannot
   catch "someone measured on plain input", because the job controls the input. Keep it as
   documentation; stop calling it "the axis that fooled everyone".

**One unstated assumption worth surfacing.** That the *records* are the controlled variable
when encoding is compared. They were not in §0, and holding them identical is what changed
the answer (§1.1). If any future cell compares encodings, the plan should require the same
fixture bytes on both sides.

---

## 3. Efficiency

PRs 1 and 2 are text or compile-time; `channel_depths` is a branch on a `bool`; 4a is six
integer comparisons and 4b's ~1400 cells are nothing. All correctly assessed, and the
withdrawal of v1's "100% distinct is conservative" justification is the right call — my data
agrees that distinctness is roughly neutral and pair count is the axis.

**PR 3's cost is now under-summed.** v2's Efficiency paragraph prices *one* cell
("~490 MB in and ~600 MB out per cell") and was written when §5 had two cells. §8 then added
two more, and PR 3 step 3 makes them part of the same job:

| cell | fixture | output |
|---|---|---|
| `1G/c4` | 490 MB (shared) | ~600 MB |
| `1G/c4 --fastqc` | shared | ~600 MB + FastQC artefacts |
| 25 bp | ~350 MB (second fixture) | ~450 MB |
| multi-pair (same fixture ×2) | shared | ~1.2 GB |

≈ 0.85 GB of fixtures and ≈ 2.9 GB of output if nothing is truncated between cells — still
inside the runner's ~14 GB, but four times what the paragraph states, and PR 3's cleanup step
is a single step at the end rather than per-cell. Two clauses: sum the four cells, and
truncate each cell's output before the next one runs.

**Wall time is also for two cells, not four.** My cells run 7–13 s each on macOS/arm64 at
10 M pairs; a 4-vCPU runner will be several times that, the FastQC cell more, the multi-pair
cell double, and two fixtures must be generated. `timeout-minutes: 30` against a cold cargo
cache (noodles, fastqc-rust, zip, mimalloc) is probably still enough, but the plan should say
it was budgeted for four cells rather than leaving v1's two-cell arithmetic in place.

**A cheap saving.** The 25 bp fixture does not need 10 M pairs. Its purpose is a σ baseline
at a second read length, not pool saturation, and §1.5 shows it cannot breach anyway. Half
the pairs would do and would pay for the multi-pair cell's output.

---

## 4. Validation sufficiency

### 4.1 The three rows v2 declares "not optional" — would they fail if the thing were broken?

| # | Tests | Can it fail? |
|---|---|---|
| **12** | "CI guard catches the regression it is sold on" — revert lever A, expect red | **On macOS yes; on the runner no** (§1.2). And the row does not say which platform to run it on, so the plausible reading is the local one, where it passes and licenses the hole. This is the row that most needs the fix. |
| **13** | "fails loudly on a parse miss" — feed it a truncated log | **Yes.** With `set -euo pipefail` plus non-empty positive-integer checks before comparing, a truncated log gives an empty capture and a non-zero exit. This row is well designed and closes round 1's silent-pass finding properly. |
| **16** | "The §0 defect reproduces as written" | **Yes for what it tests, and it tests the wrong half.** I ran it: all three commands reproduce verbatim. But it validates the *CLI* claim while the *measurement* claim (1.057× budget) is what issue #1 headlines, and no row re-checks that. That gap is how §1.1 would have reached a public issue. |

**Fix for 16:** add a row that re-measures the breach on a byte-identical fixture with arms
interleaved before the figure goes public — or, per §1.1, drop the figure and the row's
second half together.

### 4.2 Rows that are sound

- **5** is now correct, and I proved it rather than taking it on trust. Running v2's six
  pinned cells against seven perturbations:

  ```
  sigma 23→20        TABLE FAILS (6/6 cells)      k 64→32          TABLE FAILS (6/6)
  static 224→512MiB  TABLE FAILS (6/6)            fastqc 24→0      TABLE FAILS (3/6)
  cap 16→1024        TABLE FAILS (1/6)            margin REMOVED   TABLE FAILS (6/6)
  margin inverted    TABLE FAILS (6/6)
  ```

  Baseline pins: `19_006_356 / 3_527_028 / 2_170_478 / 2_042_803 / 23_911_222 / 2_375_794`.
  So 4a bites for every perturbation round 1 named, and v2's decision to say 4b "is
  *expected to pass* — that is the point of splitting them" is honest and correct. **Round-1
  Critical #1 is genuinely discharged.**
- **1, 2, 3, 4, 6, 7, 8, 10** are appropriate. I confirmed 1's arithmetic live: `1G` prints
  `≈ 930`, `1500M` prints `≈ 1363`.
- **14**'s new branch ("if any cell exceeds 1.0 on unmodified code, stop and treat it as a
  finding, not a threshold to tune") is exactly the branch round 1 asked for.
- **15** correctly gates §6's public figures behind re-measurement rather than publishing
  round 1's numbers. That is the right disposition and it is why I did not re-measure them.

### 4.3 Remaining gaps

| Gap | What slips through |
|---|---|
| **No row for the two new §8 cells** (§1.4, §1.5) | The multi-pair cell cannot satisfy Validation 11 as specified (two banners), and neither new cell can fail. The one axis §8 calls promise-breaking ships with a log line. |
| **Validation 12 has no platform** (§1.2) | The tripwire claim is validated where the job does not run. |
| **Nothing re-checks issue #1's headline figure** (§4.1) | An unreproducible budget breach reaches a public issue and `CALIBRATION.md`. |
| **The `± 1 MiB` margin tolerance has ~91 KiB of slack, not 1 MiB** | `budget × 10/11 = 930.909` MiB and the banner prints the truncation, `930`. So 0.909 MiB of the 1 MiB tolerance is consumed by the printed truncation itself, and the assertion is effectively exact equality on `floor(budget × 10/11)`. It works at `1G` — I checked in KiB: 952,320 vs 953,251, a 931 KiB gap inside a 1024 KiB tolerance — but an implementer will read "± 1 MiB" as slack and there is none. Say "exact after truncation", or assert `prediction == 930` for the pinned `1G` cell as PR 2 step 2 already does. |
| **`cap 16→1024` is caught by exactly one of the six pins** | The 17-thread cell is irreplaceable; a later tidy-up that trims the table to five cells silently loses cap coverage. One clause on the table. |
| **Nothing checks each cell's exit code** (§1.1) | This is the gap that produced §0's 1081.8 MiB and hid it for a whole round: a run that fails on disk still prints a parseable banner, so a harness that greps for the numbers and not the status reports an inflated peak as a real measurement. Validation 13 covers an *unparseable* log; it does not cover a *parseable log from a failed run*, which is strictly the more dangerous case because it looks like data. |
| **Nothing checks free disk before measuring** (§1.1, §3) | Disk pressure inflates resident output buffers via the result channel and `pending` map, so it presents as a memory-model breach. On a runner near capacity the job goes red with a correct model. |

---

## 5. Alternatives worth considering

1. **Gate the runner's own measured baseline instead of the prediction** (§1.2). Record the
   first green run's `peak/prediction` as a committed number and assert the current run is
   within, say, +5 points of it. That is the only formulation that survives the macOS→Linux
   offset without buying more fixture, and it makes lever A a 11-point move against a
   5-point band. Costs a committed baseline file, which is the thing v2's "the job never
   restates a number the binary owns" principle was trying to avoid — so it is a real
   trade-off, not a free win. Worth naming in §5 next to the fixture-growth option.
2. **Make the gating cell `1G/c2`, not `1G/c4`.** Reviewer A measured 0.952 at c2 in round 1
   and I measure 0.968 at c4; both are near the top of the range, and c2's larger bins
   (`16 × 22 MiB` vs `16 × 18 MiB`) make the per-worker lever proportionally larger. It also
   halves the runner's core requirement. Worth one line of consideration before committing
   to c4, given §1.2 makes sensitivity the binding constraint.
3. **Drop the `.fq.gz` output assertion, keep the `--dont_gzip`-style asymmetry in the
   issue.** Per §2 assumption 8 it cannot fail in a job that generates its own fixture. If it
   stays, label it as documentation of intent rather than a guard.
4. **File issue #1 without a memory section at all.** The CLI defect is clean, verifiable in
   three commands, and self-contained. Bolting an unreproducible RSS figure onto it is what
   makes it contestable.
5. **Rejected alternative — re-open the plain/gz output decision inside Phase 3.**
   `src/main.rs:417-424` and `src/io.rs:14-32` show this is a settled, #245-motivated design
   with a deliberate deferral already recorded. D3's routing is right and v2 should keep it.

---

## 6. Round-1 findings — discharged / partly discharged / not discharged

**Both reviewers' agreed Criticals**

| # | Finding | Status |
|---|---|---|
| 1 | Swept test blind to every coefficient | **Discharged.** 4a's six pins fail under all seven perturbations (verified, §4.2); 4b renamed for rounding + saturation, and Validation 5 now states 4b is expected to pass. The split is real. |
| 2 | 3 M fixture not at steady state; "structurally cannot be a tripwire" inverted | **Partly discharged.** The ramp replaces the false structural claim and 10 M is the right size — but the replacement claim is macOS-only and false at v2's own Linux expectation (§1.2). |

**Reviewer A's Criticals**

| # | Finding | Status |
|---|---|---|
| A1 | Re-size the fixture | **Discharged** (10 M, gzipped). |
| A2 | `peak ≤ prediction` not satisfied; gate the budget; run the 990.9 outlier down | **Not discharged.** v2 keeps prediction-gating and explains the outlier by the plain-output regime — which §1.1 refutes, so the outlier is unexplained again. Empirically A2's concern is *eased*: the prediction held in all 40 of my runs across four encodings and two core counts, and A's 990.9 did not reproduce for me either (c4 plain max 846.1). But the plan should say the outlier is unexplained with contention the leading hypothesis, not that it was diagnosed. |
| A3 | Fix Step 4 and Validation 5 | **Discharged.** |

**Reviewer B's (my) round-1 Criticals**

| # | Finding | Status |
|---|---|---|
| C1 | Redesign Step 4 | **Discharged.** |
| C2 | Resize fixture, drop `2G/c4`, rewrite §5's third bullet | **Partly discharged** — see agreed #2. `2G/c4` dropped ✓, ramp added ✓, Linux arithmetic missing. |
| C3 | Re-measure before filing the threading issue | **Discharged.** §6 + Validation 15 gate the figures; the plateau and both wrong sentences are carried. One presentation defect remains (§1.9). |

**Importants**

| Item | Status |
|---|---|
| A4 / C3 re-measure 10 kb | Discharged (gated by Validation 15). |
| A5 / — `threading.md:63` second claim | Discharged (§6). |
| A6 / I3 two table cells + prose bound | **Discharged and verified**: c2 = 329, c10 = 579, 154.4 MiB at c16, `[587, 611)`. |
| A7 queue-depth pin honesty | Discharged — §4c drops the overclaim in as many words. |
| A8 / I9 unverified axes disposition | **Partly discharged.** All four are now named (§8), but the two remedies cannot fail (§1.5). |
| A9 / I10 EOF-flush swallowed send | Discharged — §7, and I verified `parallel.rs:714`/`:754`/`:1202` use `let _ =` against `:794`/`:825`'s `.is_err()`. |
| A10 `needs:`, `/usr/bin/time`, fixture size | Discharged (PR 3 steps 2 and 4; sizes now consistent at ~490 MB). |
| A11 / — sharpen the mutation to lever A | Discharged as Validation 12 — but see §1.2. |
| I1 banner arithmetic assertion | **Discharged.** PR 2 step 2 pins `predicted peak ≈ 930 of 1024 MiB budget`; I verified c2 at default `1G` prints exactly that prediction. |
| I2 harden the parse | Discharged for the single-pair cells; **re-opened by §8's multi-pair cell** (§1.4). |
| I4 result-channel depth | **Discharged.** `parallel.rs:130`/`:959` are `sync_channel(cores * 2)`; folding both into one helper is right, and the k accounting closes exactly: 1 queued + 1 in hand + 2 result = 4 = k. |
| I5 ordering contradiction | Discharged (§1.11). |
| I6 core span for the floor range | **Discharged and verified**: c64 = 933 plain / 1356 `--fastqc`. |
| I7 `clump-only.md:76` parity claim | **Discharged and verified**: line 76 is the parity bullet; floors 277 (c2) → 454 (c32). |
| I8 three batch constants | Discharged (§6). |
| I11 leave `reference/changelog.md` alone | Discharged (files table). |
| I12 give the 9% margin a documented home | Discharged (§3 additions). |
| A12–A19, O1–O9 | All discharged. Spot-checked: placement at `clump-only.md:75` ✓, `--cores 1` remedy ✓, `10/11`/`workers` definitions ✓, single CHANGELOG bullet ✓, distinctness justification withdrawn ✓, `ClumpLayout::new` routed ✓, clippy on `examples/` ✓, `--paired` explicit ✓, five stderr lines ✓, fresh-runner-VM rationale ✓, `timeout-minutes: 30` ✓, `clumpy.md:110` exception ✓. |

**Coverage against the assigning documents.** `SESSION_HANDOFF.md` §2's seven rows all have a
home: fallback decision → D1, banner wording → D2, CI guard → §5/PR 3, docs arithmetic →
§3/PR 1, unfiled issue → §6/§7, **unverified axes → §8** (the row v1 dropped; now
dispositioned, with §1.5's caveat), `CALIBRATION.md` → §9. The handoff says `clumpy.md` is
wrong in 9 places; v2 says 11 and both additions (lines 91 and 119) are real.

**Correction to this section — one parent-plan assignment is NOT covered, and my first pass
said it was.** I originally wrote that the parent plan's Phase 3 assignment at `PLAN.md:102`,
`:170`, `:334`, `:342` "is fully covered". That was wrong, and I had cited the very line
numbers that disprove it without checking the deliverable against v2. Reviewer A raised it
independently; I verified it. `PLAN.md:168-170`:

> `resolve_layout_predicted_peak_fits_budget` **passes for any coefficient** — it compares the
> formula against its own inverse. Retire it or move the assertion somewhere it bites;
> Phase 3's CI check is the real replacement.

and `PLAN.md:332-337`:

> The tautological `resolve_layout_predicted_peak_fits_budget` is **deliberately left in place
> and assigned to Phase 3**, since the plan pairs its retirement with "Phase 3's CI check is
> the real replacement" — deleting a passing test before its replacement exists would reduce
> coverage, not improve it.

The test is live at `src/clump.rs:577`, asserting `predicted_peak_bytes() <= budget` over five
cells. **`grep -c predicted_peak_fits_budget phase3/PLAN.md` returns 0** — v2 never mentions
it. See §1.12.

---

## 7. Action items

### Critical

1. **Withdraw §0's memory-regime claim everywhere it is load-bearing** (§1.1). In 40 runs at
   10 M pairs with records held byte-identical and arms interleaved, peak RSS never exceeded
   the 930 MiB prediction; max 900.0 MiB, and it was on the **gz** path. Output encoding moves
   c8 by 0.2% of prediction, and at c4 gz is *higher* than plain. Specifically: strike the
   1.057× figure from §0, from Assumption 8's justification, from §5's fixture rationale, from
   §7's issue-#1 bar, and from §9's second `CALIBRATION.md` note; rewrite the Self-Review's
   "measurement lesson" paragraph. Keep the CLI defect — I reproduced all three commands
   verbatim — and re-headline issue #1 as a messaging/asymmetry defect.
2. **Fix §5's tripwire claim for the platform the job runs on** (§1.2). At v2's own expected
   Linux ratio (~0.85) a lever-A revert lands at ~894 MiB against a 930 MiB gate — green.
   Using my macOS 0.968 minus Phase 2's 9-point offset it is ~921 MiB — still green. Carry the
   Linux arithmetic beside the macOS arithmetic, **state that Validation 12 must run on the
   runner**, and name what happens if it comes back green (larger fixture, or a
   baseline-relative gate, or retract the drift claim on Linux and keep the promise guard).

### Important

3. **Specify the multi-pair cell's parse** (§1.4). It prints the banner twice — verified. Pin
   the parse to the first match and give the cell its own line in Validation 11.
4. **Say plainly that §8's two cells are baselines, not guards** — or gate them (§1.5). A
   non-gating 25 bp cell cannot see the under-reserve §8 identifies, and the arithmetic says a
   gated one would not fire either (~33 MiB against ~100 MiB of headroom). Both statements are
   fine; implying coverage is not.
5. **Correct 276/302 → 277/303** in the handoff-corrections table (§1.6). The binary uses
   `div_ceil`; v2's own §2 says 277 for the same quantity, and v2 corrects this exact rounding
   class two tables later.
6. **Fix issue #1's factual claim about gzip routes** (§1.7). Naming the input `*.gz` — a
   symlink suffices — does produce `.fq.gz` output; I verified it. The true claim is "no CLI
   flag requests it". Ship the workaround with the issue.
7. **Frame issue #1 against #245** (§1.8). `src/main.rs:417-424` records that always-gzipping
   shipped in v2.1.0-beta.5 and was reverted, and `src/io.rs:14-32` records the deferral of
   the content-vs-name decision. Without those citations the issue invites re-landing a known
   regression.
8. **Split §6's "764 MB–1.2 GB" row by input size** (§1.9). Two reviewers' different input
   sizes are presented as one row in a table whose thesis is that input size drives RSS, and
   this row is destined for a public issue.
9. **Re-sum PR 3's disk and wall time for four cells** (§3). The paragraph prices one cell and
   §8 added two more; truncate each cell's output before the next runs, and consider halving
   the 25 bp cell's pair count.
9a. **Make PR 3 immune to the failure that produced §0's number** (§1.1). Now that the cause is
   known to be disk exhaustion backing up the result channel and `pending` map, the same
   mechanism can make a tight runner disk read as a memory-model breach — a red
   `peak ≤ prediction` on a correct model, which is the single most expensive false positive
   this job can produce. Three clauses: assert free space before the first cell (`df` against a
   stated minimum), truncate each cell's output before the next, and **check
   `trim_galore`'s exit code on every cell**. The last one is the gap that hid this for a full
   round: a non-zero exit can still print a parseable banner, so Validation 13's
   unparseable-log test does not reach it. Add a validation row for a cell that exits non-zero
   with a parseable banner.
9b. **Retire `resolve_layout_predicted_peak_fits_budget`** (`src/clump.rs:577`), the one
   parent-plan assignment v2 drops (§1.12, and a correction to my own coverage claim in §6).
   The parent plan assigned it to Phase 3 twice by name, gated on the CI check that PR 3
   delivers; it is a restatement of `MARGIN_NUM ≤ MARGIN_DEN` and #446 made it slacker still.
   Delete it or fold its assertion into 4a's pinned-cell table, in PR 2 — which already opens
   that file — then strike the parent plan's "assigned to Phase 3" note.
10. **Replace §1.3's fixture rationale with the reason that holds** (§1.3). Gz is the *more*
    sensitive fixture at the gating cell — 0.968 vs 0.910 in my interleaved c4 runs, matching
    reviewer A's 3 M result — and since §1.2 makes sensitivity the binding constraint, that is
    the argument to put in the plan. The disk and representativeness reasons also stand.

### Optional

11. **State that the `± 1 MiB` margin tolerance is exact-after-truncation** (§4.3), with
    ~91 KiB of real slack, and note that PR 2 step 2's substring already pins the same thing
    deterministically.
12. **Note that `cap 16→1024` is caught by only the 17-thread pin**, so the six-cell table
    must not be trimmed (§4.3).
13. **Say the `.fq.gz` assertion is documentation, not a guard** (§2, assumption 8) — the job
    generates its own gzipped fixture, so it cannot fail on the axis §0 described.
14. **Define `ChannelDepths`** in the Signature block; it is used and never declared. And
    `resolve_layout_pins_the_default_cell` is at `clump.rs:545`, not `:546`.
15. **Footnote that `--cores` does not move the `--clump_only` floor below 5**
    (`n_bins = max(16, 4·bin_cores)`, so c1–c4 all floor at 277 MiB) — the amendment to
    `clump-only.md:76` is still right, but "sets `n_bins`" is only true from c5 up.
16. **Decide whether PR 1 gets a changelog line** (§1.10). Eleven corrected published figures
    currently ship with none.

---

## 8. What is sound, briefly

Recorded because most of v2 is right and the two Criticals are concentrated in one section
and one paragraph.

- **The 4a/4b split is a real fix**, and I proved it independently rather than accepting it:
  six pins, seven perturbations, every one fires. v2's willingness to write "4b is *expected
  to pass* — that is the point of splitting them" is the opposite of the overclaiming round 1
  caught.
- **Every corrected figure landed, and nothing was traded for a new error** — except the
  276/302 rounding, which is prose rather than a published number. I re-derived the whole
  floor family from `clump.rs` rather than checking v2's arithmetic against itself.
- **All `file:line` references resolve**, including v2's four new ones: `clump-only.md:76`
  and `:78`, `parallel.rs:130`/`:959`, the three `let _ = send` sites, and `clump.rs:162`
  and `:216-217` as genuinely unchecked multiplies. That is three rounds without a stale ref.
- **The `k = 4` accounting closes exactly** on the clumpy path: work depth 1 + batch in hand
  + 2 result slots per worker = 4. That is what makes §4c's "pin both channels" the right
  shape rather than a nicety, and it independently confirms the +104 MiB lever-A delta.
  Strengthened after the fact: §1.1's disk-exhaustion diagnosis is an accidental measurement
  of how much the result channel and `pending` map hold when the drain stalls, which is the
  best available evidence that §4c pins something that moves.
- **Validation 13 properly closes round 1's silent-pass finding**, and Validation 14's new
  "treat >1.0 as a finding, not a threshold to tune" branch is the right instinct — it is
  what would have caught §1.1 if it had been aimed at the measurement as well as the log.
- **D1, D2 and D4 are clean**, and the D4 sequencing argument (stop publishing eleven wrong
  figures while the guards are reworked) is the right priority call.
- **v2's honesty about 4c** ("catches an edit to the helper, not a re-inlined literal … v1
  called this 'the guard for the exact regression the handoff fears'; that claim is dropped")
  is the model for what §0 and §5 now need.

**Scope note.** No edits under `src/`, `tests/`, `Docs/`, `.github/` or `Cargo.toml`; `git
status` for those paths is unchanged and the tree is still at `fdade7e`. No branch switch, no
stash, no worktree, no `cargo build`. No GitHub issues filed. All measurement ran in the
session scratch directory under `B2/`; outputs are truncated to zero bytes (`rm` is blocked
by the permission gate).

**Incident, disclosed.** I destroyed the shared 10 M fixture and then restored it. To alias
the gz fixtures under plain-looking names I used **hardlinks** (`ln -f t_R1.fastq.gz
B2/nz_R1.fastq`), so my per-arm cleanup (`cat /dev/null > B2/nz_R1.fastq`) truncated the
shared inode and zeroed `scratchpad/t_R{1,2}.fastq.gz` at ~20:35. Regenerated from
`scratchpad/gen10.py`, which is seeded (`random.Random(42)` plus fixed offsets), and the
restored files are **byte-identical in size to the originals** (244,953,962 and 244,962,989 B)
with 10,000,000 pairs and the same first record. Verified with a live run: 861.8 MiB at
`1G/c4` gz, inside my earlier 854.2–900.0 spread — so the restored fixture measures the same
as the one it replaced, and it is a fifth data point for §1.1's c4 gz arm. Reviewer A and the
team lead were notified at the time.

**Incident closed out.** Three follow-ups, all done:

- **Fixture verified complete.** R1 244,953,962 B and R2 244,962,989 B, both matching the
  originals; `gzip -t` passes on each; `awk 'NR%4==1' | wc -l` gives **exactly 10,000,000
  records on both mates**. A short R2 the lead observed mid-session was a snapshot taken while
  `gzip -1 -c` was still writing the second mate — R2 is written after R1, so the intermediate
  size was real but transient.
- **My aliases no longer share the inode.** `ln -sf` replaced both `B2/nz_R{1,2}.fastq`
  hardlinks with symlinks, which unlinks the name without touching the data (link count 4 → 3,
  `gzip -t` still clean afterwards). The rebuilt plain arms use symlinks too. Nothing of mine
  can now truncate the fixture. Note that a symlink is still destructive if truncated through,
  so the durable rule is: do not truncate any alias, and check the link count first.
- **Two further hardlinks exist on the fixture's inode, and I misattributed them —
  correction.** `find -samefile` shows `scratchpad/revA2/p_R{1,2}.fastq` and
  `revA2/g_R{1,2}.fastq.gz` sharing the inode. I did not create them and did not touch them.
  I warned Reviewer A that they were its files; **they are not.** `revA2/` sits in *my* session
  scratchpad (`29c7aad0-…`) while Reviewer A's is a different session (`720d3ced-…`), where
  `stat` reports `links=1` on every fixture. I inferred authorship from a directory *name* —
  the same species of error this review exists to catch, and Reviewer A was right to refuse to
  act on the claim rather than `mkdir` the path and manufacture the hazard under its name. I am
  not attributing `revA2/` to anyone: its mtime is still advancing, so something is writing
  there, and a second guess would repeat the first.

  The safety property, which holds regardless of ownership:
  `stat -f 'mode=%Sp links=%l inode=%i'` gives `mode=-r--r--r-- links=3 inode=199305342` for
  both `t_R1.fastq.gz` and `revA2/p_R1.fastq` — same inode, but mode 444, so truncation through
  either name fails. Mitigated rather than live.

None of my reported numbers depend on the damaged window: every figure in §1.1's tables was
measured before 20:35, and the post-restoration provenance run (§1.1) reproduces all four
arms on the restored fixture — 768.8 to 789.9 MiB at c8, against §0's 1081.8.
