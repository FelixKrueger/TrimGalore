# Plan review A (round 2, measured) — Phase 3, `--clumpify` memory promise (#439)

**Reviewer:** A (the lineage briefed on `phase3/PLAN.md` **v2** with a measurement budget)
**Target as briefed:** `phase3/PLAN.md` v2 · **checked against** `phase3/PLAN.md` **v3** (§A below)
**Baseline:** `dev` at `fdade7e`; `target/release/trim_galore` = `fdade7e — macos/aarch64 — built 2026-08-17T14:38:36Z`
**Round-1 reports read:** `phase3/PLAN_REVIEW_A.md`, `phase3/PLAN_REVIEW_B.md` (both unmodified by me)
**Date:** 2026-08-17

> ### Filename, and a correction to the record — read first
>
> This document was first written to `phase3/PLAN_REVIEW_A_v2.md` at 21:02 (43,496 B; extended to
> 50,153 B at 21:06). At 21:06 that path was **overwritten** by a different Reviewer-A lineage —
> the one briefed on the parent `../PLAN.md` at `2ab0918`, which states in its own header that it
> "ran zero `trim_galore` invocations". Its work is legitimate and I have not touched it; I have
> re-issued mine under a distinct filename so the two cannot collide again. This is the
> report-path lease hazard: one path, two generations.
>
> Consequently two statements now in the plan and progress file are wrong, and both are checkable:
>
> - `phase3/PLAN.md` §0.3 and `PROGRESS.md`: *"`PLAN_REVIEW_A_v2.md` for this phase was never
>   written. The agent assigned to it reviewed the **parent** plan at the **pre-refit tip
>   `2ab0918`** instead."* **It was written**, against `phase3/PLAN.md` v2, at HEAD `fdade7e`,
>   with **62 verified `trim_galore` runs** behind it. Two lineages were both called "Reviewer A";
>   only one of them reviewed the parent at `2ab0918`, and it says so itself.
> - §0.3's headline figure — *"`1G/c4` at 1050.2 MiB … `actual ≤ predicted ≤ budget`
>   unsatisfiable"* — is not mine and appears nowhere in my work. My `1G/c4` figures are
>   855.8–882.3 MiB (gz) and 813.7–843.8 MiB (plain). The strings `1050` and `unsatisfiable` do
>   not occur in this document.
>
> I make no claim about whether the parent directory's 16 August review was overwritten, by whom,
> or whether it is recoverable — `plans/` is untracked and I never wrote outside `phase3/`.
> **v3's practical conclusion still stands with one word changed:** round 2 had two reviewers of
> `phase3/PLAN.md` v2 (Reviewer B and me), and one of the three round-2 documents reviewed a
> different plan. v3's *new* material — §0.1–§0.3, §5.3, §5.4, §8's rewrite — has been seen by
> Reviewer B and, as of §A below, by me.

---

## Verdict

**Implement v3, with the four items in §A.2 fixed first.**

On v2 as briefed: *implement with the listed changes* — PR 1 (docs) was ready as written, most of
PR 2 was ready, and §0 did not survive. v3 has since retracted §0, fixed the platform hole in
§5, and independently reproduced the finding I had ranked Critical (§5.3's Linux table is my
Critical 4, arrived at separately). Credit where due: v3 is a better document than either round-2
review asked for.

**What is still wrong is now concentrated in one place, and it is new in v3.** §8 computes that a
25 bp under-charge is "~33 MiB" and concludes *"a gated cell would not fire today … the axis
remains unguarded."* Measured, on the gz path with 10 M pairs: **25 bp reads reach 1.054–1.064 ×
prediction at `1G/c4` and 1.044 × the 1024 MiB budget at `1G/c2`.** A gated cell fires
immediately, and the promise the CHANGELOG publishes is already broken on that axis. v3's
arithmetic understates the effect by roughly 3× and does so in the reassuring direction.

---

## A. Reconciliation with v3

### A.1 Discharged by v3 without my input

| My finding (against v2) | Where v3 fixes it |
|---|---|
| §0's regime claim is false; encoding is not a regime | §0.1 — retracted, with the disk-stall cause identified |
| "No other way to request gzip" is false; the rename/symlink route works | §0.2, second bullet — and it ships the workaround in the issue |
| The tripwire is blind on the runner (+104 MiB vs a ~0.85 Linux baseline) | §5.3's table — same delta, same conclusion, independently derived |
| Disk pressure can present as a memory breach | §5.4 — including "check every invocation's exit code", which is exactly the hole |
| `ChannelDepths` was never declared | Signature, line 410 |
| `--cores 1` floor quoted as 276/303 under mixed rounding | Context table — now "**277 / 303** — `div_ceil`, as the binary advises" |
| Encoding arms must share fixture bytes | Assumption at line 550 |

### A.2 Still live against v3 — four items, all measured

**1. §8's short-read arithmetic is wrong by ~3×, in the direction that produces inaction. [Critical]**

v3, §8: *"The pool term at `1G/c4` is `1.4375 × 16 × 18.13 ≈ 417 MiB` of a 930 MiB prediction, so
an 8% under-charge is ~33 MiB — **a gated cell would not fire today** against 0.85–0.97
baselines."* Measured instead, same binary, 10 M pairs, gz in / gz out:

| cell | 51 bp | **25 bp** | 25 bp vs prediction | 25 bp vs **budget** |
|---|---:|---:|---:|---:|
| `1G/c2` | max 938.0 | **1040.9, 1068.9** | **1.149** | **1.044** |
| `1G/c4` | max 882.3 | **979.8, 989.8** | **1.064** | 0.966 |
| `1G/c8` | max 787.1 | 904.6, 905.6 | 0.974 | 0.884 |

The measured 25 bp penalty at c4 is **~100–107 MiB**, not 33: scaling only the pool term by 8%
misses that at 25 bp the per-record spine and header bytes are a far larger fraction of record
bytes, so the effective σ is worse than the 1.55 the handoff recorded. Consequences for v3:

- The claim "a gated cell would not fire today" is false — it fires by 6 points at the gating
  cell's core count and by 15 at c2.
- "**Half the pair count suffices** (its purpose is not pool saturation)" now cuts the wrong way:
  at half the pairs the ramp puts the recorded baseline below steady state, so D5's committed
  band would be calibrated against a number that is not the cell's real ratio.
- This axis, not the CLI defect, is the data-integrity item in this arc: `CHANGELOG.md:30`
  publishes *"`--clumpify` keeps peak memory inside the `--memory` budget it prints"*, and at
  25 bp / `--cores 2` it does not. It deserves its own issue with these numbers, and the phase
  should say plainly that the promise is length-scoped until σ is refitted.

**2. §5.3's "Promise (absolute): `peak ≤ prediction`. Platform-independent, cannot flake." [Important]**

"Cannot flake" is the premise that lets the absolute assertion gate, and it does not hold:

- `1G/c2`, gz path, 51 bp, **12 reps**: 878.5 889.5 892.4 897.7 898.0 910.4 911.1 914.2 917.2
  926.5 927.4 **938.0** — one rep **above** the 930 MiB prediction (1.009), three within 0.5 %.
- The **multi-pair** cell v3 keeps: **929.5 MiB against 930** (2 pairs, `1G/c4`) — 0.9995.
- At the gating cell itself the margin is thin: my max 882.3 and B's 900.0 sit 48 and 30 MiB
  below the gate, against a **96 MiB** documented macOS within-cell spread (`CALIBRATION.md`) and
  a 93 MiB margin. The parent plan already recorded this at `PLAN.md:479-483`.

`--cores 2` is the *minimum legal* `--clumpify` cell, so this is not a corner. It does not argue
against v3's two-assertion design; it argues that the absolute half should be stated as "has not
flaked at c4 on 10 M pairs" rather than "cannot flake", and that a red on it is a finding about
the margin, not about CI.

**3. The multi-pair cell has no stated input, and the obvious one is rejected. [Important]**

v3 fixes the *parse* (banner printed twice, pin the first match — a real catch) but never says
what to pass. v2 said "the same fixture passed twice", and that is refused before any I/O:

```
$ trim_galore --paired --clumpify --cores 4 --memory 1G -o dup x_R1.fastq x_R2.fastq x_R1.fastq x_R2.fastq
Error: Pair 2 (x_R1.fastq, x_R2.fastq) is a duplicate of pair 1. Did you mean to pass different files?
```

Under `set -euo pipefail` that fails the whole job, gating cells included. Distinct *paths* to the
same bytes are accepted — I used hardlinks — so it is one `ln` in PR 3 step 4, but it has to be
written down. Measured for the baseline v3 wants to record: 2 pairs → 929.5 MiB (0.9995), 3 pairs
→ 908.9 MiB, against a single-pair max of 882.3 at the same cell. Retention across the pair
boundary is real, ~+47 MiB, and not cumulative beyond the second pair.

**4. Two smaller carries. [Optional]**

- **`clumpy.md:119`.** v3 line 198 still says to "strip line 119's '~100MB mark' figure". That
  figure is *correct* for the reads the page is about: the non-clumpify path at 10 M pairs ×
  51 bp, `--cores 4`, measures **85.1 MiB**. It is wrong only at long read lengths — which is §6's
  subject. Qualify by read length; stripping a true statement mis-describes the defect.
- **`phase2/calibrate.py:35`** is `re.compile(r"clumpify: (\d+) bins × (\d+) MiB; predicted peak ≈ (\d+) MiB")`.
  D2 breaks it and `calibrate.py:102` then raises `SystemExit` on every cell. Loud, not silent,
  but that harness is what Q1's RRBS re-measurement and Validation 15 depend on, and no section
  of v3 mentions it. One line, one Files-table row.
- **§6's "Lead with the gzip-output figure"** (line 339) survived §0's withdrawal without its
  justification. It happens to be right, for a reason v3 no longer states: the result channel and
  `pending` map hold *written* output, so the plain/gz gap scales with per-batch output bytes —
  tens of MiB at 51 bp under `--clumpify` (measured: 20–90 MiB, sign varying), hundreds at 10 kb
  (the parent plan's 802 vs 508). Restore that sentence as the reason, or drop the instruction.

---

## B. Method, and what it is worth

| Method | Used for |
|---|---|
| `os.posix_spawn` + `os.wait4`, `ru_maxrss` (bytes on Darwin), measuring parent holds nothing | every peak-RSS figure |
| Ran `target/release/trim_galore` at `fdade7e` | every banner/warning figure, the output-extension behaviour, the duplicate-pair rejection |
| Python re-implementation of `clump.rs` integer arithmetic (`div_ceil` where the binary uses it) + `git show 2ab0918:src/clump.rs` for the old formula | all floors, the regression table, the crossover sweep, window widths |
| Reading `src/`, `Docs/`, `tests/`, `phase2/calibrate.py`, both result TSVs | every `file:line` claim, and the §0 provenance question |

**Harness proved able to fail** before any negative was trusted: a child touching a 1500 MiB
allocation reports **1509.0 MiB**; a trivial child reports 8.2 MiB. No clamping, no inheritance.

**Contention is declared.** Reviewer B measured concurrently throughout; load averages ran 6–34
and at one point ~59 MB of pages were free. That matters in one direction: macOS reclaims under
pressure, so pressure **suppresses** peak RSS. Every breach reported here is therefore
conservative, and the one negative that mattered (§0's cell) was re-run at the quietest point
available (load 6.3–8.8, 3.9 GB free pages) with the same result.

### B.1 Provenance of every number

**Every reported cell is verified complete, not assumed complete.** `run.sh` captured each child's
exit code and full stderr; outputs were truncated afterwards, logs were not. Across all **62**
reported runs: `exit=0` in every one, and the trimming report says `Total reads processed:
10000000` for **both** mates in every one (4 of 4 mates for the 2-pair cell, 6 of 6 for the
3-pair cell). The only `exit=1` runs in my logs are six `g*` cells whose inputs had just been
zeroed by the shared-fixture truncation described in B.2 — they reported `rss_mib=7.8` with an
empty banner, were recognised as failures on sight, and appear nowhere in this report. This is the
precondition v3 §5.4 now demands of the CI job, applied to my own harness.

**Fixture lineage** — all from the shared fixture's *original, complete* content, captured before
the 20:35 truncation:

| my name | what it is | inode / bytes | used for |
|---|---|---|---|
| `../t_R{1,2}.fastq.gz` | the shared fixture; verified **40,000,000 lines = 10 M records** by `gzip -dc \| wc -l` at 20:26 | 199305342 / 244,953,962 | `g4a`, `g8a–c`, `p8a–c` |
| `g_R{1,2}.fastq.gz` | hardlink of the above | **same inode** | gz-in / gz-out |
| `p_R{1,2}.fastq` | hardlink of the above under a plain-looking name | **same inode** | gz-in / plain-out |
| `pp_R{1,2}.fastq` | `gzip -dc` of the shared fixture, 20:28 | 199529189 / 1,478,888,890 | plain-in / plain-out |
| `q_R{1,2}.fastq.gz` | hardlink of `pp_R*` under a `.gz` name | **same inode as `pp_R*`** | plain-in / gz-out |
| `gg_R{1,2}.fastq.gz` | `gzip -c pp_R*`, 20:41, level 6, 216 MB/mate | own inode | all `gg*` cells (gz-in / gz-out) |
| `s25_R{1,2}.fastq.gz` | `pp_R*` with sequence and quality truncated to 25 bp, then gzipped | own inode | all 25 bp cells |
| `h_*`, `i_*` | hardlinks of `gg_*` under distinct names | same inode as `gg_*` | multi-pair cells |

**Byte-identity across compared arms** — by identity where possible:

- plain-in/plain-out vs plain-in/gz-out (`pp` vs `q`): the two arms read the **same inode**
  (199529189, confirmed by `ls -i`). Not equal bytes — the same bytes.
- gz-in/gz-out vs gz-in/plain-out (`g` vs `p`): **same inode** (199305342).
- gz vs plain *input* arms (`gg`/`g` vs `pp`): one lineage, decompressed once from the shared
  fixture and recompressed from that decompression, so the record stream is identical by
  construction. **Verification status: construction plus per-run 10 M record counts, not a hash** —
  I truncated the fixtures on finishing, before the provenance request arrived, so post-hoc md5 is
  no longer available. Flagged rather than dropped.
- 25 bp vs 51 bp: same headers, same order, same count; sequence and quality truncated to 25.

**Reps were run as blocks, not interleaved.** This is the weakest part of my design and I do not
dress it up. Two things limit the damage: the `1G/c8` plain-output cell — §0's cell — was measured
in **two blocks 20 minutes apart** under very different machine states (load 17–34 with ~59 MB
free pages; then load 6.3–8.8 with 3.9 GB free pages) and returned 768.3–791.1 then 759.8–793.0,
the same distribution; and the cross-core sweep alternated the gz and plain arms within each core
count rather than running one arm to completion first. B's interleaved 2×2 is the stronger design
for the encoding question and I defer to it there.

**Free disk was never tight.** Plain-output runs wrote ~2.8 GB, gz-output runs ~560 MB. Free space,
logged per run in the first block and spot-checked after: **7.5–16 GiB**, minimum 7.5 GiB — at
least 2.6× the largest write. The run producing my highest plain-output figure (793.0 MiB) had
8.0 GiB free. This independently corroborates v3 §5.4 and §0.1: the withdrawn 1081.8 MiB was
measured at 1.6 GiB free while writing 2.9 GB, and a writer that cannot drain leaves the result
channel (`sync_channel(cores * 2)` — 16 slots at c8) and the `pending` map resident.

### B.2 Coordination hazards encountered

Both are worth carrying into how parallel reviews are run, not just this one:

1. The shared 10 M gz fixture was truncated by the other reviewer's cleanup at 20:35 while my
   hardlinks pointed at its inode, zeroing my inputs mid-matrix (six `exit=1` cells, discarded).
   Fixtures should be per-agent; hardlinks across agents are the specific trap.
2. `phase3/PLAN_REVIEW_A_v2.md` was overwritten by another lineage at 21:06 (see the header note).
   One path per generation, and a distinct filename per lineage.

`revA2/p_R{1,2}.fastq` and `revA2/g_R{1,2}.fastq.gz` are hardlinks to the shared fixture's inode
and are therefore now `-r--r--r--` with it. **Left alone, nothing chmod-ed:** read-only is what
makes them harmless — my one truncation attempt against `p_R1.fastq` failed with "permission
denied", as intended — and swapping in symlinks would re-point the same names at the same target.
My own fixtures (`pp_*`, `q_*`, `gg_*`, `s25_*`, each on its own inode) are truncated, returning
~3.7 GB.

### B.3 Master measurement table

`--paired --clumpify --memory 1G`, 10,000,000 pairs, 100 % distinct, macOS/arm64, `fdade7e`.
Printed prediction **930 MiB** in every cell; budget **1024 MiB**.

**51 bp** (the CI fixture's geometry):

| `--cores` | input → output | reps | peak RSS (MiB) | max / pred | max / budget |
|---:|---|---:|---|---:|---:|
| 2 | gz → gz | 12 | 878.5 889.5 892.4 897.7 898.0 910.4 911.1 914.2 917.2 926.5 927.4 **938.0** | **1.009** | 0.916 |
| 2 | plain → plain | 2 | 870.5 873.2 | 0.939 | 0.853 |
| 2 | plain → gz | 2 | 901.6 910.3 | 0.979 | 0.889 |
| 3 | gz → gz | 2 | 848.8 852.8 | 0.917 | 0.833 |
| 4 | gz → gz | 4 | 855.8 876.5 878.7 882.3 | 0.949 | 0.862 |
| 4 | plain → plain | 4 | 813.7 814.7 817.2 843.8 | 0.907 | 0.824 |
| 4 | plain → gz | 2 | 897.7 900.6 | 0.968 | 0.880 |
| 8 | gz → gz | 6 | 754.5 782.9 783.7 784.6 785.4 787.1 | 0.846 | 0.769 |
| 8 | **plain → plain** | **10** | 759.8 760.7 768.3 771.7 778.9 780.5 789.7 791.1 792.5 793.0 | **0.853** | **0.774** |
| 8 | gz → plain | 3 | 756.5 759.3 761.3 | 0.818 | 0.743 |
| 8 | plain → gz | 2 | 769.2 805.6 | 0.866 | 0.787 |
| 16 | plain → plain | 2 | 729.5 729.6 | 0.785 | 0.713 |
| 16 | plain → gz | 2 | 785.5 806.4 | 0.867 | 0.787 |

**25 bp** (same records, sequence and quality truncated; 10 M pairs processed, 0.03 % removed):

| `--cores` | input → output | reps | peak RSS (MiB) | max / pred | max / budget |
|---:|---|---:|---|---:|---:|
| 2 | gz → gz | 2 | 1040.9 **1068.9** | **1.149** | **1.044** |
| 4 | gz → gz | 2 | 979.8 **989.8** | **1.064** | 0.966 |
| 8 | gz → gz | 2 | 904.6 905.6 | 0.974 | 0.884 |
| 8 | gz → plain | 2 | 876.1 **939.3** | **1.010** | 0.917 |

**Multi-pair**, 51 bp gz → gz, `--cores 4`: 2 pairs → **929.5** MiB (0.9995 of prediction);
3 pairs → 908.9 MiB. Single-pair max at that cell: 882.3.

**Non-clumpify baseline**, 51 bp plain → plain, `--cores 4`: **85.1 MiB**.

All four input/output encoding combinations are reachable because output encoding comes from the
input's *extension* (`main.rs:424`, `naming::is_gzipped`, `io.rs:33`) while input decoding is
content-based (`format::detect_input_format`). Verified both directions: gz content named
`.fastq` → `x_R1_val_1.fq`, ASCII text; plain content named `.fastq.gz` → `y_R1_val_1.fq.gz`,
real gzip.

---

## C. Round-1 findings — discharged / partly / not

Judged against **v2**, which is what I was briefed on. v3 does not disturb any of these verdicts
except A-C2, which it improves.

**Reviewer A (round 1).**

| # | Finding | Status |
|---|---|---|
| A-C1 | CI fixture undersized, inverts §5 | **Discharged** — 10 M, gzipped, `2G/c4` dropped, structural claim replaced by the ramp. |
| A-C2 | `peak ≤ prediction` not satisfied; gate the budget instead | **Not discharged in v2**; **partly in v3.** v2 correctly showed the specific c4 rep does not reproduce, then generalised from c4/c8 to "the promise holds". It fails at c2 (938.0) and at 25 bp (1.149). v3's two-assertion design is better but still calls the absolute half "cannot flake" (§A.2 item 2). |
| A-C3 | Step 4's sweep cannot fail; Validation 5 self-contradictory | **Discharged** — 4a pinned cells + 4b renamed, Validation 5 now expects 4b to pass. The pinned set (default cell, `workers ≠ bin_cores`, FastQC cap at 16 **and** 17, c2, c32) bites on every constant in `clump.rs`. |
| A-I4 | Queue-depth pin covers the helper, not the call sites | **Discharged** — stated as a limitation. |
| A-I5 | 10 kb figure input-limited (~2283 plateau) | **Discharged** — §6 table, "~20–23×". |
| A-I6 | `threading.md:63`'s second wrong claim | **Discharged**. |
| A-I7 | Two wrong rows + "≤ 150 MiB" | **Discharged and re-verified**: c2 = 328.90 → 329, c10 = 578.05 → 579, widest = **154.4 MiB at c16** (global max over `--cores` 1–256). |
| A-I8 | Unverified axes need a disposition | **Discharged** in form; v3's §8 rewrite is better but carries the arithmetic error in §A.2 item 1. |
| A-I9 | EOF-flush swallowed send | **Discharged** — `parallel.rs:714`, `:754`, `:1202` vs `:794`, `:825`, all exact. |
| A-I10 | `needs:`, `/usr/bin/time`, fixture size quoted twice | **Discharged**. |
| A-I11 | Sharpen Validation 10 to a lever-A revert | **Discharged** — Validation 12. |
| A-O12…O19 | banner numbers, `clump-only.md` placement, `--cores 1` remedy, core span, define `10/11`+`workers`, CHANGELOG bullet, drop the distinctness justification, `ClumpLayout::new` | **All discharged**; every figure re-checked holds. |

**Reviewer B (round 1).** C1, C2, C3 discharged as above (C2 with the same caveat). I1–I12 all
discharged — including I11's grep claim, which is right: `clumpy.md:96`, `:108` and
`changelog.md:199` are the only 512 MiB references in `Docs/`. O1–O9 discharged, with the `:436`
uniformity nit still open (v3 line 470 keeps it).

**No traded errors.** Every figure corrected after round 1 is right: c2 = 329, c10 = 579, 154.4 MiB
at c16, `[587, 611)`, 281–590 / 334–1012 **for `--cores` 2–32**, c64 = 933 / 1356, crossover
exactly `--cores ≥ 11`, `620 → 352` at c16 without `--fastqc`, `need ≥ 311 MiB` at c6, `1500M →
1363`, and `--cores 32 --fastqc --memory 1011M` → "need ≥ 1012 MiB" (so Validation 9's method
works). All recomputed from `clump.rs` and, where printable, run.

---

## D. Validation sufficiency

**The three rows v2 declared not optional:**

- **Validation 12 (revert lever A → red).** Would fire on macOS at c4 (882 + 104 = 986 > 930); would
  **not** fire at c8 (787 + 104 = 891); marginal-to-green on the runner. v3 §5.3 now agrees and
  adds the relative band, which is the right fix. Run it *on the runner*, not locally.
- **Validation 13 (truncated log → non-zero).** Can fail only if the parse is factored out of
  inline YAML so a log can be fed to it. v3 §5.4 adds the sharper sibling — a *parseable* log from
  a *failed* run — which is the case that actually bit this arc, and it needs the exit-code check
  to be per-invocation, not per-step.
- **Validation 16 (the §0 defect reproduces).** The three commands do reproduce; the breach half is
  gone with §0. v3 splits it correctly.

**Gaps that remain:**

| Gap | What slips through |
|---|---|
| Nothing proves the lower bound bites | It does: a fallback to plain mode measures **85.1 MiB** = 0.09 of prediction, so D5's lower band catches it with 10× to spare. One line puts the number on record. |
| No row states what a `> 1.0` ratio on the 25 bp cell means | It will be 1.05–1.06 at c4 on day one (§A.2 item 1). Without a stated expectation this reads as a broken guard. |
| Nothing asserts the fixture is what it claims | The generator prints pairs / read length / distinct count to stderr and no assertion consumes it. One `grep` makes a degenerate fixture loud — the failure `CALIBRATION.md` records as having invalidated an entire earlier matrix. |
| The multi-pair cell's baseline is at the promise boundary | 929.5 of 930. D5 records a baseline from the first green run; that cell's baseline is 1.00, which leaves no room above it. |

---

## E. What is sound

- **Every docs figure in §3 is right.** All eleven corrections and the six additions recomputed
  from `clump.rs` and run where printable: 281–590 / 334–1012 for `--cores` 2–32, 933 / 1356 at
  c64, 311 at c6, 9 MiB at `1G/c8`, 44 MiB at `4G/c8`, and lines 91 / 96 / 96–97 / 98 / 102 / 104 /
  108 / 110 / 113 / 119 all exact.
- **The 4a/4b/4c split** is the right answer to round 1's agreed Critical and is honest about what
  each guard proves. 4b's saturation value is real: `clump.rs:162`, `:216`, `:217` all use unchecked
  `*`, verified.
- **D2's arithmetic assertion is the strongest available guard, and it is available.** `1G/c2` prints
  `≈ 930` deterministically, `1500M` prints `≈ 1363`, and the tripwire case at
  `integration_preflight_tripwire.rs:592/599` runs exactly `["--clumpify", "--cores", "2", "a.fastq"]`
  at the default `--memory 1G`. The `/(1000*1000)` mutation does break the substring.
- **D1's analysis holds.** The bail does not name the mode; `small_memory_budget()` is 1 GiB, above
  every `workers: 1` floor; nothing asserts on the bail text; the `--cores 1` remedy caveat is right;
  and `anyhow` printing the cause chain from `fn main() -> Result<()>` means `.with_context()` keeps
  the "try ≥ N MiB" advice.
- **Handoff §2 is fully dispositioned** — all eight rows. The two axes v3 promoted, short reads and
  multi-pair, are exactly the two where I measured a breach; the instinct was right even where the
  arithmetic under-estimated it.
- **v3's §5.4 is the best paragraph in the arc.** It generalises a harness failure into a property of
  the system under test — disk pressure presents as a memory breach — and then names the false
  positive it would produce and the wrong response it would invite. My independent free-disk record
  (7.5–16 GiB against ≤2.9 GB written) is the control that makes it checkable.

**Scope note:** read-only on source. `git status` on `src/`, `tests/`, `Docs/`, `.github/`,
`Cargo.toml` and `Cargo.lock` is clean; tree still at `fdade7e`; no branch switch, stash, worktree,
or `cargo build`; no GitHub issues filed. All measurement ran under
`…/scratchpad/revA2/`; the 62 run logs (`out_*/stderr.txt`, `out_*/measure.txt`, `mp*.txt`,
`nc*.txt`) are intact there if any figure needs re-deriving. My only writes in the repo are this
file and, earlier, `phase3/PLAN_REVIEW_A_v2.md` — which another lineage has since replaced.
