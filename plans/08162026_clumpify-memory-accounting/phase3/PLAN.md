# Plan — Phase 3: defend and correct the `--clumpify` memory promise (#439), v3

**Parent plan:** `../PLAN.md` v3 (Phase 3 = its "defend and correct" stage)
**Issue:** [#439](https://github.com/FelixKrueger/TrimGalore/issues/439) — deliberately still open
**Reviews:** round 1 → `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`; round 2 → `PLAN_REVIEW_B_v2.md`, and `PLAN_REVIEW_A_v2.md` (targets v2, delivered after v3 was written — see §0.3)
**Baseline:** `dev` at `fdade7e`, confirmed equal to `origin/dev`
**Prerequisites:** none outstanding.

## Revision history

- **v1 — reviewed, three Criticals, two agreed.** Its docs arithmetic and `file:line`
  references were right; both its *guards* were useless. The swept property test was blind to
  every coefficient (algebraically), and the CI fixture was undersized by 13–17 points of
  ratio, from which v1 drew the structural conclusion that the job "cannot be a drift
  tripwire".
- **v2 — guards rebuilt, but it introduced a false measurement of its own.** v2's §0 claimed
  the plain-text output branch was a distinct memory regime reaching 1.057 × `--memory`. It is
  not. That claim propagated into an assumption, the fixture rationale, a queued public issue
  and a `CALIBRATION.md` note before round 2 refuted it.
- **v3 — this document.** Retracts §0, fixes the platform hole in §5 that v2 inherited from
  v1's reasoning style, and records that round 2 had one reviewer, not two.

## 0. What round 2 established

### 0.1 v2's §0 is withdrawn

Reviewer B ran a controlled 2×2 — input encoding × output encoding, records held
byte-identical across all four arms, reps interleaved rather than blocked — and in 40 runs at
10 M pairs peak RSS never approached the printed prediction. I then reproduced it
independently. Both sets, macOS/arm64, `--memory 1G`, prediction 930 MiB:

| cell | v2's §0 claimed | reviewer B (40 runs) | mine, re-measured |
|---|---|---|---|
| `1G/c4` plain out | 908.3 MiB | 846.1 | — |
| `1G/c4` gz out | 873.2 MiB | **900.0** (max of all arms) | — |
| `1G/c8` plain out | **1081.8 MiB** | 826.1 | **756–789 MiB** (5 reps) |
| `1G/c8` gz out | 814.4 MiB | 814.7 | — |

Two things make this conclusive rather than a disagreement. B's gz figure at `1G/c8` is
**814.7 against my earlier 814.4** — the harnesses agree to 0.3 MiB wherever the measurement
was sound. And B's decisive arm used a genuinely decompressed 1.48 GB plain file in exactly
the cell v2 claimed 1081.8, measuring **785.4**.

**Output encoding is not a memory regime.** With records held identical, switching output from
gz to plain moves the c8 mean by 0.2% of prediction, and at c4 gz is *higher* than plain.

**Cause of v2's figure: my own harness.** Those readings were taken with ~1.6 GiB of free disk
while each rep wrote ~2.9 GB of plain output, and `io::PendingOutput` writes a temp file beside
the previous output, so consecutive reps into one directory needed ~5.8 GB. A stalled writer
backs up the result channel (`sync_channel(cores * 2)` — 16 slots at c8) and the `pending`
map, both of which hold written-but-unflushed buffers. **Disk exhaustion inflates peak RSS.**
The harness grepped for `peak=` and never for the exit status, so a failed run's parseable
banner read as data.

That has a consequence for this plan beyond the retraction, and it is §5.4.

### 0.2 What survives: a smaller, real CLI defect

All three commands reproduce verbatim on `fdade7e` (B verified independently):

```
--clumpify --cores 2 in.fastq          → in_trimmed.fq  ; banner still says "(gzip level 1)"
--clumpify --gzip --cores 2 in.fastq   → in_trimmed.fq  ; "--gzip is deprecated … Ignoring."
--clumpify --dont_gzip --cores 2 …     → Error: … mutually exclusive (clumping plain text is pointless)
```

Two of v2's three claims about this were wrong and are dropped:

- **"No way to obtain gzipped output from plain input" is false.** Naming the input `*.gz` — a
  symlink suffices — yields `.fq.gz` output, because `io::is_gzipped` (`src/io.rs:33`) reads
  the *filename* while `FastqReader` sniffs bytes. The true claim is narrower: no CLI *flag*
  requests it.
- **"Always gzip" is not an open question.** `src/main.rs:417-424` records that v2.1.0-beta.5
  always gzipped, broke pipelines globbing `*.fq`, and was reverted for Perl parity under
  **#245**; `src/io.rs:14-32` records the content-vs-filename decision as deliberately
  deferred. Output-follows-input is settled design.

So the defect is: a deprecation message that is false for plain input, plus `--clumpify`
rejecting `--dont_gzip` as "pointless" while `modes/clump-only.md:78` documents it as
supported. **No memory claim.** Filed per D6.

### 0.3 Round 2 had two reviewers of this phase, and a third document about another plan

**Correcting an earlier version of this section.** Three round-2 documents exist and only two
review this phase:

| document | target | tip | invocations |
|---|---|---|---|
| `PLAN_REVIEW_B_v2.md` | phase3 v2 | `fdade7e` | 40 runs, controlled 2×2 |
| `PLAN_REVIEW_A_v2_measured.md` | phase3 v2 | `fdade7e` | **62 verified runs** |
| `PLAN_REVIEW_A_v2.md` | **`../PLAN.md`** (parent) | **`2ab0918`** | **zero** |

Two Reviewer-A lineages were in flight. The measuring one wrote `phase3/PLAN_REVIEW_A_v2.md` at
21:02; that path was overwritten at 21:06 by the parent-plan lineage, so it re-issued under
`_measured`. **The `1050.2 MiB` / "unsatisfiable independent of coefficient choice" figure
belongs to the parent-plan lineage, not the measuring reviewer** — an earlier version of this
plan and my reporting attributed it wrongly. The measuring reviewer's c4 figures are 855.8–882.3
(gz) and 813.7–843.8 (plain).

Its one methodological weakness, self-disclosed: reps were **blocks, not interleaved**. Mitigated
by measuring the `1G/c8` plain cell in two blocks 20 minutes apart under very different load and
free-disk conditions, which gave the same distribution (768–791, then 760–793). Six `exit=1` runs
— cells whose inputs the 20:35 truncation zeroed — reported `rss 7.8` with an empty banner and
are excluded. Its harness was proved able to fail first: a deliberate 1500 MiB allocation
reported 1509.0.

### 0.3a The original claim that round 2 was single-reviewed

`PLAN_REVIEW_A_v2.md` for this phase was never written. The agent assigned to it reviewed the
**parent** plan at the **pre-refit tip `2ab0918`** instead, and wrote into the parent
directory's existing `PLAN_REVIEW_A_v2.md`, overwriting the 16 August review of that document
(unrecoverable — `plans/` is untracked; all artifacts are now snapshotted under
`../_snapshot_20260817/` and locked read-only).

Its headline figure — `1G/c4` at 1050.2 MiB against a 1024 MiB budget, making
`actual ≤ predicted ≤ budget` "unsatisfiable independent of any coefficient choice" — **does
not transfer**, and I verified why at source. At `2ab0918`:

```
usable    = memory_budget − STATIC_OVERHEAD_BYTES   (512 MiB, no margin term)
predicted = STATIC_OVERHEAD_BYTES + dyn_bytes   ⇒   predicted ≈ budget, zero headroom
```

Two corrections from reviewer A, since this section is the paper trail for a retraction: the
budget at that cell is **1025 MiB, not 1024** (the harness computes `round(19 × 108/4) + 512` and
the TSV records 1025), so the gap was 25.2 MiB rather than 26. And **the margin is not what
resolved it** — a margin only *shrinks* `predicted` and therefore *tightens* the bind. What
resolved it was **`B` moving 4.40%** (18.96 → 18.126 MiB). So the parent Goal's constraint
"don't move the bin size at the default" is the thing that had to give, which is the actionable
half of that finding.

**Reviewer A subsequently reviewed this phase** (`PLAN_REVIEW_A_v2.md`, targeting v2, delivered
after v3 was written). It reached v3's two Criticals independently — sever §0, and do not gate
on an unmeasured Linux baseline — and withdrew its own staging proposal in favour of D5. Four
findings from it survive against v3 and are folded in: pin the sum rather than the two literals
(§4c), run Validation 12 in both mutation forms, give the 25 bp cell a band and its own issue
(§8), and record `2G/c4` as a deliberate blind spot with a tracking issue (§5.2, §7).

It also re-derived the arithmetic layer independently — all seven trim-path floors, all six
`clump_only` cells, the `c1`/`c64` endpoints, the plain-path 352/454, the crossover, the
154.4 MiB window, both bin sizes and both banner values: **twenty-plus claims, zero errors**,
with `div_ceil` confirmed as the reproducing convention.

**Provenance caveat that matters.** A ran **zero** `trim_galore` invocations, so no RSS figure
in either of its reports is its own; everything is marked code-verified, arithmetic-verified, or
unverified, and every *measurement* in this plan is marked unverified by it. So the arithmetic
and code-reference layers are now genuinely dual-reviewed; the **measurements** rest on
reviewer B plus my own re-measurement.

### 0.4 A real breach, found where v2 dismissed it — short reads

**This supersedes §8's arithmetic and is the data-integrity finding of the arc.** v3 estimated
an 8% σ under-charge at 25 bp as "~33 MiB — a gated cell would not fire today". That estimate
was ~3× too small, in the reassuring direction. Reviewer A measured it; I reproduced it
independently on a 25 bp fixture with headers identical to the 51 bp one, so read length is the
only variable. 10 M pairs, gz in → gz out, ample disk, every run `exit=0` with
10,000,000 records on both mates and `.fq.gz` output:

| cell | mine (reps) | reviewer A | vs 930 prediction | vs **1024 MiB budget** |
|---|---|---|---|---|
| `1G/c2` | 1048.1 / 1049.6 / 1075.1 / 1082.6 | 1040.9 / 1068.9 | 1.13–1.16× | **1.024–1.057× — BREACHES** |
| `1G/c4` | 994.3 / 996.3 / 1000.0 | 979.8 / 989.8 | 1.07× | 0.97× |
| `1G/c8` | 903.0 | — | 0.97× | 0.88× |
| `1G/c16` | 904.3 | — | 0.97× | 0.88× |

Seven reps at `1G/c2` across two independent harnesses; **all seven exceed the budget.** The
gradient is mechanistically right: at c2, `n_bins = 16` and `B` is largest (22 MiB), so the bin
pool dominates and σ's under-estimate bites hardest; it narrows as `B` shrinks with cores.

**It is user-reachable.** `--small_rna` is a shipped preset and miRNA libraries trim to ~22 nt,
so short trimmed records are a supported case, not a synthetic one.

**It falsifies `CHANGELOG.md:30` as written** — "**`--clumpify` keeps peak memory inside the
`--memory` budget it prints**". That entry is still in `### Unreleased`, so nothing has shipped
with it; it must not ship unqualified. Note this is the same line that was rewritten twice
during Phase 1/2 review for over-claiming, and strengthened after Phase 2.

**Consequences for this plan.** §8's short-read cell is not a baseline that "cannot fire" — it
fires immediately, and at `--cores 2` it fires on the budget rather than the prediction. So
§8's remedy is wrong in both directions: the cell should be **gating**, and "half the pair count
suffices" now cuts the wrong way, because the ramp puts a half-size fixture below steady state
and D5's committed band would be calibrated against the wrong number.

**Scope decision required (see Questions).** Phase 3 cannot both claim to defend the promise and
ship knowing the promise fails at 25 bp / `--cores 2`.

## Decisions

| # | Decision | Chosen |
|---|---|---|
| D1 | `--clump_only` below the floor | **Documented hard failure.** Reordering is the whole job; a fallback would write a byte-copy of the input at full I/O cost and exit 0. Behaviour unchanged, so zero regression risk. |
| D2 | Banner margin visibility | **Name the budget**: `predicted peak ≈ 930 of 1024 MiB budget`. |
| D3 | The plain-output defect | **Its own issue**, not Phase 3. Premise narrowed by §0.2 — no breach. |
| D4 | Docs corrections | **Ship first**, as a docs-only PR. |
| D5 | The CI gate | **Baseline-relative band** on a committed ratio, plus the absolute promise assertion. See §5.3. |
| D6 | Issue scope | **File, narrowed** to the false message and the mode asymmetry; cite #245 and `io.rs:14-32`; ship the `*.gz` workaround. |
| D7 | The short-read breach (§0.4) | **Scope the claim and ship Phase 3.** Qualify `CHANGELOG.md:30` to the measured read lengths, make the 25 bp cell **gating** at a worsening-band, and file the σ refit as its own issue. Phase 3 still delivers the docs, the guards and the CI job; it stops claiming a promise that does not hold everywhere. The refit needs a calibration matrix across read lengths that Phase 3 was never scoped to produce. |

## Context

### Corrections to the handoff (all three verified by both round-1 reviewers)

| Handoff | Actual |
|---|---|
| banner prints `≈ 931 MiB` | **930 MiB** (integer truncation) |
| floor range `276–590` / `302–1012` MiB | **281–590** / **334–1012** MiB for `--cores` 2–32. The `--cores 1` values are **277 / 303** — `div_ceil`, as the binary advises, not the 276/302 truncations v2 carried |
| `--clump_only --fastqc --cores 16 --memory 800M` aborts | **resolves**; that floor is 775 MiB |

v2 corrected v1's `c10` truncation and then reproduced the same rounding error two tables
earlier. Fixed above.

### The `--clump_only` regression surface

`workers: 1`, `fastqc_threads = cores`, `div_ceil` to MiB:

| `--cores` | old floor | new floor | previously-working budgets that now abort |
|---:|---:|---:|---|
| 2 | 535 MiB | 329 MiB | none — floor dropped |
| 8 | 566 MiB | 513 MiB | none |
| 10 | 580 MiB | 579 MiB | none |
| 11 | 586 MiB | 611 MiB | [587, 611) MiB |
| 16 | 620 MiB | 775 MiB | [620, 775) MiB |
| 32 | 728 MiB | 876 MiB | [728, 876) MiB |

Widest window **154.4 MiB at `--cores 16`**. Crossover is exactly `--cores ≥ 11` (both
reviewers swept independently). Without `--fastqc` the floor dropped everywhere. This
narrowness is what makes D1 cheap.

## Behaviour

### 1. Banner (D2)

```
clumpify: 32 bins × 9 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
```

Second figure is the parsed `--memory` (`1500M` → 1500, verified). Sub-MiB budgets truncate to
0 but are unreachable — the smallest resolvable budget is 281 MiB.

**Guard consequence.** After D2 the predicted peak has no unit token of its own; the `MiB`
belongs to the budget. So a clause keyed on `"MiB budget (gzip level"` pins the *budget's*
divisor, and changing `predicted_peak_bytes() / (1024*1024)` to `/ (1000*1000)` would print
`≈ 975 of 1024 MiB budget` and pass. That re-opens the divisor bug Phase 2's post-review fix
#4 closed. Fixed by asserting the arithmetic (PR 2 step 2).

### 2. `--clump_only` below the floor (D1)

Behaviour unchanged — non-zero exit, no output file. The bail gains the mode's name, and the
behaviour is documented. B verified the context wrap cannot break a test: nothing asserts on
the bail text and `small_memory_budget()` is 1 GiB, above every `workers: 1` floor.

At `--cores 1` the floor is 277 MiB and lowering `--cores` is unavailable, so the docs must say
raising `--memory` is the only remedy there. Note `n_bins = max(16, 4·bin_cores)`, so `--cores`
does not move the floor at all below 5.

### 3. Docs arithmetic (D4)

Eleven wrong statements in `clumpy.md`; every line number and replacement value verified by
both round-1 reviewers.

| Line | Wrong | Correct |
|---|---|---|
| 91 | "in an attempt to get the predicted peak RSS ≤ `--memory`" | drop the hedge — coefficients are fitted to measured peak RSS on both platforms |
| 96 | `usable = memory − 512 MiB` | `usable = memory × 10/11 − static` |
| 96–97 | one flat term "(FastQC + allocator + runtime overhead)" | static = 224 MiB, **plus** 24 MiB per FastQC thread only when requested, capped at 16 charged threads |
| 98 | `4 × usable / (5·n_bins + 7·cores)` | `16 × usable / (23·n_bins + 64·workers)` |
| 98 | `cores` | `workers` — `--cores` on the trim path, 1 under `--clump_only` |
| 102 | "32 bins × 12 MB" at `1G/c8` | **32 bins × 9 MiB** |
| 102 | "32 bins × 66 MB" at `4G/c8` | **32 bins × 44 MiB** |
| 108 | "~535 MiB" | **281 MiB** at `--cores 2` |
| 108 | "535–730 MiB … any sensible core count" | **281–590** without `--fastqc`, **334–1012** with, **for `--cores` 2–32** |
| 108 | "fixed 512 MiB reservation for FastQC, …" | 224 MiB trim reservation; FastQC conditional and capped |
| 113 | "need ≥ 552 MiB" at `--cores 6` | **need ≥ 311 MiB** |

Additions: name the core span (at `--cores 64` the floor is **933 MiB** plain / **1356 MiB**
with `--fastqc`, above the 1 GiB default); define `× 10/11` and `workers`; give the ~9% margin
a documented home as deliberate run-to-run headroom (Phase 2's reviewer B predicted D2's
support question verbatim); **date** the stale savings figures rather than fabricating
replacements, and resolve the pre-existing contradiction between line 104 and the RRBS row;
add a `--clump_only` exception to line 110.

**Line 119 is not stripped** (reviewer A, correcting v3): "~100MB mark" is *correct* for short
reads — A measured 85.1 MiB non-clumpify at 10 M × 51 bp, `--cores 4`. **Qualify it by read
length** instead, pointing at the threading page for the long-read case. v3's instruction to
strip it would have removed a true statement.

In `clump-only.md`: amend line 76's "`--cores` accepted for interface parity" claim, which the
new below-floor paragraph contradicts — with the c5 caveat above. Leave
`reference/changelog.md:198-199` alone (synced page); it and `clumpy.md:96`/`:108` are the only
512 references in `Docs/`.

### 4. Deterministic drift guards

Round 1 proved v1's sweep blind to every constant. B verified v2's replacement fires — six
pinned cells against seven perturbations:

```
sigma 23→20  6/6 fail    k 64→32  6/6 fail    static 224→512MiB  6/6 fail
fastqc 24→0  3/6 fail    cap 16→1024  1/6 fail    margin removed  6/6 fail    margin inverted  6/6 fail
```

Baseline pins: `19_006_356 / 3_527_028 / 2_170_478 / 2_042_803 / 23_911_222 / 2_375_794`.

- **4a — pinned-cell table.** Six hand-computed `(budget, LayoutInputs) → exact
  bin_byte_budget` cells: the default cell, `workers ≠ bin_cores`, the FastQC cap at **16 and
  17** threads, `--cores 2`, `--cores 32`. The 17-thread cell is the *only* one that catches a
  cap change, so the table must not be trimmed — state that on it. Reviewer A also checked the
  `n_bins` mapping, which the perturbation matrix omits: changing the `16` breaks the c2 cell and
  changing the `4` breaks the default cell, so it is covered.

  **Two gaps sit one layer *above* the pins, and neither is caught by anything** (reviewer A):

  - **`LayoutInputs` construction is entirely unguarded.** `fastqc::resolved_threads` has five
    call sites and **zero tests** — `fastqc.rs`'s test module has eight `apply_fastqc_args_*`
    tests and nothing for it. So reverting to `fastqc_threads: Some(cores)` under
    `--fastqc_args "-t 2"`, or to `workers: cores` in `clump_only.rs`, passes all six pins, 4b,
    **and the CI job** — the last because §5.2's cells are `--cores 4`, where all three roles are
    numerically equal and the split is invisible. This undoes Phase 2's post-review item 2,
    which is the fix that made the prediction a bound under `--fastqc_args -t N`. Add direct
    tests for `resolved_threads`, and a construction-level assertion that the three roles differ
    where they should.
  - **Lowering `MIN_BIN_BYTES` is undetectable.** All three floor tests are relative to the
    function under test — the same vacuity §4 proves about v1's sweep. Raising it to 4 MiB is
    caught, because four of the six pins sit below 4 MiB; **lowering it to 256 KiB is caught by
    nothing**, yet it moves every floor figure PR 1 publishes. Fix with one absolute pin:
    `clumpify_min_memory_bytes(&ins(2, false)) == 281 MiB`, which pins the constant *and* the
    headline docs number in the same assertion.
- **4b — the sweep, renamed** for the two properties it does prove: rounding direction (it
  catches the real Phase 2 iteration-#2 bug at 28 violations) and saturation. The latter is not
  idle — `clump.rs:162` and `:216-217` use unchecked `*`, safe only via `resolve_layout`'s
  `saturating_mul(16)`. Validation 5 states 4b is *expected to pass* under coefficient
  perturbation; that is the point of splitting them.

- **4d — retire `resolve_layout_predicted_peak_fits_budget`** (`clump.rs:577`). This is the one
  deliverable the parent plan assigns to Phase 3 **by name, twice** — `../PLAN.md:168-170`
  ("Retire it or move the assertion somewhere it bites; Phase 3's CI check is the real
  replacement") and `:332-337` ("deliberately left in place and **assigned to Phase 3**"). v1
  and v2 both proposed *replacing* it with the sweep; v3 dropped the reference entirely, which
  is how a named assignment went missing (caught by both round-2 reviewers, from opposite
  directions — one from the tree, one from the plan).

  It asserts `predicted_peak_bytes() <= budget` over five cells, which by §4's own algebra is a
  restatement of `MARGIN_NUM ≤ MARGIN_DEN` — the identical defect §4 spends a section proving
  about the proposed sweep, already shipping. #446 made it **more** tautological: the chain is
  now `predicted ≤ budget × 10/11 ≤ budget`, clearing with 9% to spare.

  The parent plan's stated precondition is not merely "Phase 3 does it" — it gives a reason, at
  `../PLAN.md:332-336`: *"deleting a passing test before its replacement exists would reduce
  coverage, not improve it."* Reviewer A's Critical: v2 and an earlier v3 scheduled the deletion
  in **PR 2** while the replacement lands in **PR 3**, which is exactly the window the parent
  ruled out — and Validation 6a would report "discharged" inside it. Substantively harmless
  (the test is vacuous, so no real coverage is lost), but it should not be implicit.

  **So the deletion moves to PR 3**, after the CI check exists. Fold its intent into 4a's pinned
  table, and strike the parent's "assigned to Phase 3" note in the same pass — **which is also
  where `../PLAN.md:329`'s stale figure gets fixed**: it describes the pin as asserting
  19,884,107 bytes, where the shipped test asserts **19,006,356** (`clump.rs:551`). Parent-doc
  rot, confirmed by recomputation.
- **4c — channel depths, both channels.** `clump.rs:247-248` says `k = 4` pays for "the queued
  batch, the batch in hand, and the compressed output waiting on the result channel". The
  accounting closes exactly: work depth 1 + batch in hand + 2 result slots = 4. v1 pinned only
  the work channel; `sync_channel(cores * 2)` at `parallel.rs:130` and `:959` was unguarded,
  and the parent plan names it as the next lever someone will pull.

  **§0.1 is accidental evidence that 4c pins something that moves**: the mechanism that
  inflated my erroneous reading was exactly those slots filling when the writer stalled. Cite
  it there.

  **Pin the sum, not the two literals.** Reviewer A verified the decomposition at
  `clump.rs:246-248` exactly: work channel 1 + 1 in hand + `result_per_core` 2 = 4 = `k`. So
  assert `work + 1 + result_per_core == 4` with a message naming `k`, *alongside* the
  value pins. Two independent literals do not communicate that their **sum** is what `k`
  encodes, and the sum is the invariant a future editor would break. A's check also confirms
  4c is a pure refactor, so Integration's "no output bytes change" claim holds.

Stated honestly: 4c catches an edit to the helper, not a re-inlined literal at the call sites.
The extraction still removes duplication and puts the `k`-coupling note beside
`dyn_denominator`, where the repo's comment rule wants an invariant.

**And that concession is PR 3's strongest justification, which §5 currently omits** (reviewer
A). The inline conditional is what stands at `parallel.rs:117`/`:943` today, so re-inlining is
the likely shape of a future regression — and **the CI job is the only guard that catches the
re-inlined form**. §5 argues for the job from fixture saturation; it should lead with this.
Hence Validation 12 runs **both** mutations: edit the helper, and re-inline the literal.

### 5. CI guard

#### 5.1 Fixture: 10 M pairs, 51 bp, gzipped

B measured the ramp at `1G/c4`, monotone throughout: 500 K → 0.378, 1.5 M → 0.754,
3 M → 0.804, 6 M → 0.875, 10 M → ~0.95. Nothing is at steady state below ~10 M, so v1's 3 M
proved nothing.

Gzip the fixture for three reasons — and note the strongest one is not the one v2 gave:

1. **Gz is the *more* sensitive fixture at the gating cell** — B's interleaved c4 runs give
   gz 0.968 vs plain 0.910, matching reviewer A's round-1 3 M result (0.888 vs 0.835). Since
   §5.3 makes sensitivity the binding constraint, this is the reason that matters.
2. CI should measure the branch users run.
3. 490 MB rather than 2.96 GB on disk.

v2's justification ("plain output is a different memory regime") is withdrawn with §0.

#### 5.2 Cells

`1G/c4` gating, `1G/c4 --fastqc` logged. `2G/c4` dropped — its pool is 672 MiB, needing ~2.3×
the data to saturate, so it is the least informative per byte and the sole reason a runner-RAM
assumption was load-bearing. Plus §8's two cells.

**`2G/c4` therefore becomes a deliberate, recorded blind spot** (reviewer A). It is the cell
that carried the original breach, §4c names it as "the next lever someone will pull", and after
this phase it has neither a CI cell nor a tracking issue. Dropping it from CI stays the right
call on saturation grounds — but the residual gets an issue (§7), so it is rediscovered as a
known item rather than as novel.

#### 5.3 The gate (D5): a baseline-relative band

v2 claimed `peak ≤ prediction` is itself the drift tripwire. **That is false on the platform
the job runs on**, and it is v1's error in platform form — v2 stated the Linux expectation four
bullets from the arithmetic and never multiplied them together:

| | macOS | **Linux (v2's own expectation)** |
|---|---:|---:|
| baseline at `1G/c4` | 0.939–0.968 → 873–900 MiB | **~0.85 → 790 MiB** |
| + lever-A revert (+104.2 MiB) | 977–1004 | **894** |
| vs the 930 MiB gate | red | **green — misses by 35 MiB** |

**The delta is a range, not a point** (reviewer A). `+104.2 MiB` is `cores × σ × B`, but the
model's own per-worker unit is `1·B` rather than `σ·B` — `dyn = B × (1.4375·n_bins + 4·workers)`
— which gives `4 × 18.13 = 72.5 MiB`. A 44% spread. **State the range, because the conclusion
survives it and is therefore stronger:** 862 and 894 are both under the 930 gate, and 7.8 and
11.2 points both breach `+0.05`. Validation 12 records the *measured* delta so the band's real
headroom is known rather than assumed. Either way the delta is real bytes held and does not
scale with platform.

So the job asserts **two** things:

- **Promise (absolute):** `peak ≤ prediction`. Platform-independent, and **it has not flaked at
  `--cores 4` on 10 M pairs at 51 bp** — which is a narrower claim than v2's "cannot flake".
  Reviewer A measured `1G/c2` gz at **938.0 against a printed 930** (1.009, 1 of 12 reps, 3
  more within 0.5%), my own single c2 rep was 927.5, and multi-pair `1G/c4` came to 929.5 of
  930. Headroom at the gating c4 cell is 30–48 MiB against `CALIBRATION.md`'s 96 MiB macOS
  spread. So a red here is a **finding about the margin**, not automatically a regression — and
  `--cores 2` is too tight to gate on.
- **Drift (relative):** `peak/prediction` within a committed band of the recorded baseline —
  fail above `baseline + 0.05`, fail below `baseline − 0.10`. A lever-A revert is an ~11-point
  move and breaches the upper bound on any platform; a silent fallback or a fixture that
  stopped filling the pool is a multi-point *drop* and breaches the lower bound. Neither needs
  a guessed fixture size.

The baseline lives in a committed file, populated from the first green run and updated
deliberately in the same PR as any legitimate improvement. This is a real trade-off: it breaks
v2's "the job never restates a number the binary owns" principle. It is accepted because that
principle was protecting against the binary moving its own goalpost, which the margin
assertion below already closes, whereas no self-describing formulation survives the
macOS→Linux offset without sizing a fixture blind.

Also asserted: `prediction == budget × 10/11`, so deleting `MARGIN_NUM` cannot inflate the
prediction to the budget and pass. Note the `± 1 MiB` tolerance is **exact after truncation** —
`budget × 10/11 = 930.909` MiB prints as `930`, consuming 0.909 of the 1 MiB — so an
implementer must not read it as slack. PR 2 step 2's substring already pins the same thing
deterministically at zero CI cost; this is the belt to that braces.

#### 5.4 Disk is not independent of memory

**This is §0.1's consequence and nothing in v2 covers it.** The job writes output per cell; with
§8's cells the unmanaged total is ~3.7 GB against a runner's ~14 GB. If disk gets tight, the
same mechanism makes **disk pressure present as a memory-model breach** — a red gate on a
correct model, which is the most expensive false positive this job can produce, because the
obvious response is to loosen the gate.

Three clauses: assert free space before the first cell, truncate each cell's output before the
next runs, and **check every invocation's exit code**. Validation 13 tests an *unparseable*
log; a **parseable log from a failed run** is strictly more dangerous because it looks like
data. That is the hole my own harness fell through and PR 3 inherits it as specified.

### 6. The `threading.md` issue

v1's 812 MiB headline is input-limited — 400 MB of 10 kb reads is ~4.7 batches of 4096 records
against ~20 slots of resident capacity. Both round-1 reviewers measured past it. **Report by
input size, not merged** (v2's single "764 MB–1.2 GB" row combined two reviewers' different
inputs in a table whose thesis is that input size drives RSS):

| read length | input | c4 | c8 | source |
|---|---|---:|---:|---|
| 10 kb | 400 MB | 768.5 / 813.2 | 817.4 | A / B |
| 10 kb | 764 MB | 1366.0 | — | A |
| 10 kb | 1.2 GB | 1808.8 | 2052.3 | B |
| 10 kb | 1.5 GB | 2121.2 | — | A |
| 10 kb | 3.0 GB | **2282.6** ← plateau | — | A |

Headline: **~20–23× the documented ~100 MB**. Drop "flat in cores past 4" — at 1.2 GB, c8
exceeds c4 by 244 MiB. `--cores 1` at ~10 MiB reproduces across three harnesses and is the
sharpest fact: the pool is bypassed.

Also: `threading.md:63` carries a second wrong claim v1 missed — *"Peak RSS scales with read
length, **not input file size**"* — refuted at fixed read length by the table above. **Three**
batch constants, not two: `parallel.rs:33 BATCH_SIZE = 4096` sizes the round-robin batches in
the measured path, and its docstring makes the same read-length-blind assumption. Lead with the
gzip-output figure, quoting plain text as the worst case.

Scope: **file only.** The fix touches the path every mode shares.

### 7. Issues to file

| Issue | Content | Bar |
|---|---|---|
| **`--clumpify` plain output** (D6) | the false deprecation sentence; the `--clumpify` / `clump-only.md:78` asymmetry; cite **#245** and `io.rs:14-32` so nobody re-lands the reverted always-gzip; ship the `*.gz` symlink workaround. **No memory claim** | user-visible wrong message |
| **`threading.md` memory claims** | §6's table, the plateau, both wrong sentences, three batch constants | user-visible wrong documentation |
| **EOF bin flush swallows a failed send** | `parallel.rs:714`, `:754`, `:1202` use `let _ = …send(…)` where `:794`/`:825`/`:1193`/`:1228` check `.is_err()`; a dead worker silently drops the final partial batch and can exit 0. **`:1193` and `:1202` are nine lines apart in the same function and the neighbouring comment cites #434 — so this is a gap in #434's own fix**, which is worth saying in the issue | both round-1 reviewers: clears the filing bar more easily than the docs sentence |
| **`--clump_only` predicted peak** | prints bins × MiB but no prediction, so #439's visibility never reaches the mode D1 just made a hard failure | a recommendation living only in a plan rots |
| **Short-read σ refit** (D7) | **The `--memory` promise is breached at 25 bp / `--cores 2`** — 1.02–1.06× budget, 7 reps, two harnesses (§0.4). σ is treated as a constant but is a *spine* factor: 25 bp packs 1.57× more records per bin byte than 51 bp, so the under-charge scales with record count, not byte count. The fix is to charge the spine per record; needs a calibration matrix across read lengths | **data integrity** — the sharpest item in the arc |
| **`2G/c4` residual** | the cell that carried the original breach, dropped from CI on saturation grounds; §4c's "next lever" (result channel + `pending` map) is the named candidate | otherwise rediscovered as novel |

Also routed: `ClumpLayout::new` (`clump.rs:150`) is public, has no production caller, and
validates neither invariant `resolve_layout` enforces — any future caller can build a layout
whose `predicted_peak_bytes()` is arbitrary, and the CI guard would faithfully report the lie
and pass. Make it `pub(crate)` or assert `bin_byte_budget >= MIN_BIN_BYTES`.

### 8. Unverified axes — baselines, not guards

`SESSION_HANDOFF.md` §2's last row, dropped by v1. Two cells are added — and B is right that
as *non-gating* cells they cannot discharge the concern, so v3 says what they are:

- **Short reads (25 bp).** σ ≈ 1.55 measured against the shipped 1.4375, so there is a known
  length where the model under-reserves. The pool term at `1G/c4` is
  `1.4375 × 16 × 18.13 ≈ 417 MiB` of a 930 MiB prediction, so an 8% under-charge is ~33 MiB —
  **a gated cell would not fire today** against 0.85–0.97 baselines. Half the pair count
  suffices (its purpose is not pool saturation).

  **Superseded by §0.4 — the arithmetic above is ~3× too small and the conclusion is wrong.**
  Measured, the 25 bp cell fires immediately: 1.07× prediction at `--cores 4` and **1.02–1.06×
  the budget at `--cores 2`**, seven reps across two harnesses. So this cell must be
  **gating**, not logged, and it must be **full-size** — "half the pair count suffices" would
  calibrate D5's band below steady state. What it gates depends on the scope decision in
  Questions; if Phase 3 ships before the refit, the cell gates at a recorded band that fails on
  *worsening* and the plan states plainly that the axis is already breached.
- **Multi-pair.** By Phase 2's own mechanism — mimalloc retaining the pool's pages across a
  phase boundary, which is why FastQC is an addend — pair 2's pool may sit on pair 1's retained
  pages. **It prints the banner twice**, so PR 3's parse must pin the first match and the
  "banner present" check must be `grep -q`, not an equality on 1. One prediction against one
  whole-process peak, which is the point: the promise is per-process.
- **Single-end and low-duplication/WGBS:** explicitly deferred, named so they stay in the paper
  trail.

### 9. `CALIBRATION.md` amendment

Item 1 recommends FastQC as a `max`; the shipped addend is correct because mimalloc keeps the
pool's pages resident across the trim→FastQC boundary (at `--cores 16` with a 64 MiB pool the
process peaks at 438.8 MiB — the trim reservation plus ~215 MiB). Amend in place with a dated
note; do not rewrite. **v2's second note is not written** — the plain-output regime it recorded
does not exist.

**And the replacement note must not tell that document its numbers were plain-path.** Reviewer A
proved Phase 2's calibration was already gz-in/gz-out, from `phase2/calibrate.py:73`/`:109` plus
non-zero `out_bytes` in both TSVs. So the note records only that output encoding measures
neutral (0.2% at c8) and that the **sign is inverted at c2/c4/c16** — gz runs 20–90 MiB *above*
plain, which is the opposite of what v2's §0 implied.

**Separately: `phase2/calibrate.py:35`'s banner regex breaks under D2** and exits on every cell.
No other section mentions it, and **Q1's RRBS re-measurement depends on that harness** — so the
regex must be updated in the same breath as the banner, or the re-measurement cannot run.

## Signature

```rust
// src/clump.rs, beside dyn_denominator so the coupling is visible where k lives.

/// Work-channel and result-channel depths. A clumpy batch is a whole bin, so shallower
/// queues cap resident bins; `k` in `dyn_denominator` was fitted at these depths.
pub struct ChannelDepths {
    pub work: usize,
    pub result_per_core: usize,
}

pub fn channel_depths(clumpy: bool) -> ChannelDepths {
    if clumpy {
        ChannelDepths { work: 1, result_per_core: 2 }
    } else {
        ChannelDepths { work: 2, result_per_core: 2 }
    }
}
```

```rust
// src/main.rs — D2
eprintln!(
    "clumpify: {} bins × {} MiB; predicted peak ≈ {} of {} MiB budget (gzip level {})",
    layout.n_bins,
    layout.bin_byte_budget / (1024 * 1024),
    layout.predicted_peak_bytes() / (1024 * 1024),
    memory_bytes / (1024 * 1024),
    cli.compression,
);
```

```
examples/mk_clump_fixture.rs
    cargo run --release --example mk_clump_fixture -- <pairs> <read_len> <out_prefix>

Writes <prefix>_R{1,2}.fastq.gz. Fixed-seed inline LCG, no `rand` dependency, byte-identical
across platforms. Prints pairs, read length and distinct-sequence count to stderr before
writing. Must clear `cargo clippy --all-targets --release -- -D warnings`, which covers
`examples/` — a hand-rolled LCG is `needless_range_loop` / `cast_possible_truncation` country.
```

## Implementation outline

### PR 1 — docs only (D4)

1. `clumpy.md`: the eleven corrections plus the additions in §3.
2. `clump-only.md`: D1's below-floor subsection — **after** the "Composes with:" list, not
   inside it (line 75 is a bullet) or under its own `## Memory` heading — and amend line 76.
   Do **not** quote the bail text; PR 2 changes it and PR 1 would be stale for a cycle.
3. Leave `reference/changelog.md` alone.
4. `npm run build` in `Docs/`, including the no-silent-render-corruption check.

PR 1 carries no changelog line: `CHANGELOG.md` is not a docs file and the repo rule forbids
mixing. PR 2's #439 bullet extension should therefore name the docs corrections too.

### PR 2 — code and guards

1. D2 banner (`main.rs:207-213`).
2. Tripwire, arithmetic not labels (`integration_preflight_tripwire.rs:599`): the case runs at
   the default `--memory 1G` with `--cores 2`, printing `predicted peak ≈ 930 of 1024 MiB
   budget` deterministically. Assert that exact substring — it pins both divisors, the margin
   and the wording. Keep the `"MiB; predicted peak ≈"` clause for the bins label.
3. D1 bail context at the four `clump_only.rs` sites (275/421/764/929) via `anyhow::Context`.
   No fallback, no exit-path change.
4. Q4 banner tidy in the **same** PR, so no release ships mixed terminology:
   `clump_only.rs:289`, `:784`, `:1031` → `N bins × M MiB each` (`:436` already omits
   "budget"). Five stderr lines change in total.
5. 4c `channel_depths`; four call sites (`parallel.rs:117`, `:943`, `:130`, `:959`); pin both
   arms, message naming the `k` coupling.
6. 4a pinned-cell table (six cells, with the do-not-trim note).
7. 4b sweep, renamed.
8. Direct tests for `fastqc::resolved_threads` (five call sites, currently zero tests) and the
   absolute floor pin `clumpify_min_memory_bytes(&ins(2, false)) == 281 MiB` — §4a's two gaps.
9. `ClumpLayout::new` → `pub(crate)` or invariant-asserting.
10. **CHANGELOG (D7):** extend #439's `### Unreleased` bullet, and **qualify the "keeps peak
    memory inside the `--memory` budget it prints" claim to the measured read lengths** — it is
    false at 25 bp / `--cores 2` (§0.4). This must not ship unqualified.

### PR 3 — CI guard (after PR 2 merges; never in PR 1)

1. `examples/mk_clump_fixture.rs`.
2. Job `clumpify-memory`, `ubuntu-latest`, `if: github.event_name != 'schedule'`,
   `timeout-minutes: 30`, `shell: bash` with `set -euo pipefail`, **no `needs:`**. Verify
   `/usr/bin/time` exists (`command -v`, or install `time`).
3. **Assert free disk before the first cell** (§5.4).
4. Generate the 10 M × 51 bp gzipped fixture into `$RUNNER_TEMP`, plus the **full-size** 25 bp
   fixture (§0.4 — half-size would sit below steady state). **The multi-pair cell needs
   distinct paths to the same bytes**: `--paired A_R1 A_R2 A_R1 A_R2` is refused with
   "Pair 2 … is a duplicate of pair 1", which under `set -euo pipefail` kills the whole job.
   One `ln -s` per mate (reviewer A verified this route works).
5. Run each cell under `/usr/bin/time -v`, `--paired` explicit, stderr to its own log.
   **Check each invocation's exit code** — `/usr/bin/time -v` exits with the child's status, so
   one check covers it, **but piping into `grep`/`tee` replaces that status with the filter's**,
   and `grep` succeeding on a failed run's parseable banner is precisely §5.4's hole. Capture to
   a file, check the status, *then* parse. Truncate each cell's output before the next runs.
6. Parse `Maximum resident set size (kbytes)` and `predicted peak ≈ (\d+) of (\d+) MiB budget`,
   **pinning the first match** (multi-pair prints two). Assert captures are non-empty positive
   integers before comparing; compare in **KiB on both sides**.
7. Assert per gating cell: banner present (`grep -q`); `peak ≤ prediction`; `peak/prediction`
   inside `[baseline − 0.10, baseline + 0.05]`; `prediction == budget × 10/11`. Echo both ratios
   always. The `.fq.gz` output check stays as **documentation of intent** — the job generates
   its own gzipped fixture, so it cannot fail on that axis.
8. Comment that the measuring parent must stay trivial (Linux inherits a child's `ru_maxrss`
   from the forking parent; `/usr/bin/time` is ~2 MB and safe).
9. Cleanup — rationale is this job's own disk, not later jobs; each GitHub job gets a fresh VM.

### Issues and the `CALIBRATION.md` amendment — any time.

## Efficiency

PRs 1 and 2 are text or compile-time. `channel_depths` is a branch on a `bool`; 4a is six
integer comparisons; 4b's ~1400 cells are trivial.

PR 3, summed for **four** cells rather than v2's two: ~0.85 GB of fixtures and ~2.9 GB of
output if nothing is truncated between cells, which §5.4 now requires anyway. Wall time is
7–13 s per cell on macOS/arm64 at 10 M pairs; a 4-vCPU runner will be several times that, the
FastQC cell more, multi-pair double, and two fixtures must be generated — `timeout-minutes: 30`
against a cold cargo cache is budgeted for that.

The 100%-distinct fixture is deterministic and cheap. v1's claim that it is *conservative* is
withdrawn: B's 10 M/100%-distinct ratios sit on Phase 2's 10 M/52%-distinct ones, so
distinctness is neutral and pair count is the axis. The warning that it must never be reused
for a compression claim stays.

## Integration

- **Writes:** five changed stderr lines. No output-file bytes; Perl byte-identity untouched.
- **Downstream:** exactly one in-repo assertion depends on the banner (both round-1 reviewers
  verified by grep).
- **Ordering:** PR 3 after PR 2. If PR 3 landed first its regex would match nothing — now a
  loud failure under `set -euo pipefail`, but still wrong. `ci.yml` is not a docs file.

## Assumptions

1. **Fixed:** no coefficient moves in Phase 3 — `resolve_layout_pins_the_default_cell`
   (`clump.rs:545`) must pass unmodified. That is all Phase 3 needs; v2's broader "the model is
   correct as shipped" rested on the withdrawn §0.
2. **Fixed:** `predicted_peak_bytes()` is the exact inverse of `resolve_layout` — within 23 B
   of `budget × 10/11` at `1G/c4`, ~312 B across cores ≤ 32. This one property is load-bearing
   in two places with opposite signs: it makes the CI parse single-source and it makes 4b
   vacuous under coefficient change.
3. **Fixed:** `--clumpify` requires `--cores ≥ 2`, so the docs floor range starts there.
4. **Fixed:** `ubuntu-latest` has ≥ 4 usable cores. Dropping `2G/c4` removes the RAM half. The
   job is **insensitive to runner core count** — `--cores` is explicit and FastQC's threads
   derive from it via `resolved_threads`; only allocator page behaviour and I/O timing are
   runner-dependent.
5. **Reclassified:** fixture sizing determines coverage *per platform*, not a cost preference.
6. **Open:** every ratio here is macOS/arm64. Linux is unmeasured; the ~0.85 expectation is
   Phase 2's offset applied to a macOS figure, which is why D5 gates on a *recorded* baseline
   rather than a predicted one.
7. **Not assumed:** that the plain-output branch is in scope. D3.
8. **New:** that records are the controlled variable whenever encodings are compared. They were
   not in v2's §0, and holding them identical is what changed the answer. Any future cell
   comparing encodings must use the same fixture bytes on both sides.

## Validation

| # | Verify | How | Expected |
|---|---|---|---|
| 1 | Banner arithmetic | `--memory 1G`, `--memory 1500M` | `≈ 930 of 1024`; `≈ 1363 of 1500` |
| 2 | No coefficient moved | `resolve_layout_pins_the_default_cell` unmodified | passes |
| 3 | Tripwire pins the divisors | change the peak divisor to `1000*1000` | fails |
| 4 | Both channel depths pinned | set work → 2, then result → 4 | fails in both cases |
| 5 | 4a bites where 4b cannot | seven perturbations, one at a time | 4a fails each; **4b expected to pass** |
| 6 | 4b bites on its own property | invert the margin; make floor/resolve rounding disagree | 4b fails both |
| 6a | **The named assignment is discharged** | `grep -c resolve_layout_predicted_peak_fits_budget src/clump.rs`; and the parent plan's `:332-337` note | **0**; note struck |
| 7 | `--clump_only` below the floor | a below-floor cell | non-zero exit, no output, message names the mode |
| 8 | `--clump_only` above the floor, incl. both uBAM sites | existing suites + `--fastqc --cores 16 --memory 800M` | all resolve |
| 9 | Docs figures | run each cell; for the floor-range endpoints, which nothing prints, reimplement `clumpify_min_memory_bytes` or read the top end from `--cores 32 --fastqc --memory 1011M` | identical |
| 10 | Docs build | `npm run build` in `Docs/` | clean; LaTeX renders |
| 11 | CI guard reports | first run, per cell | banner present, exit 0, `peak ≤ prediction`, ratios logged; multi-pair cell has its own row |
| 12 | **The drift band catches the regression it is sold on, in both its forms** | revert lever A **two ways** — (a) edit the helper, (b) re-inline the literal at `parallel.rs:117`/`:943` — and run the job **on the runner** | **red both times**. Not locally: macOS goes red at any fixture size and would license the hole. Form (b) is the one only this job can catch |
| 13 | Fails loudly on an unparseable log | truncated log | non-zero exit |
| 14 | **Fails loudly on a parseable log from a failed run** | a cell run under artificially low free disk | non-zero exit, and the message names disk, not memory (§5.4) |
| 15 | Linux baseline | first green runs | record it; **if any cell exceeds 1.0 on unmodified code, treat it as a finding, not a threshold to tune** |
| 16 | Issue figures before they are public | re-measure the 10 kb cell at ≥ 8× the records, gz output, exit codes checked | plateau ≥ 2 GiB |
| 17 | The §0.2 defect reproduces as written | the three commands | as shown |

Validations 12, 13, 14 and 16 are not optional. Every false measurement in this arc came from
an instrument that could not report its own failure — prove each of these can.

## Questions or ambiguities

- **Resolved (D1–D6).**
- **Resolved — the ratio ceiling.** Superseded by D5's band. A fixed 0.92 ceiling was
  unworkable: at a stressful fixture the ratio spans 0.849–1.065, so it would be red most runs.
- **Resolved — stale savings figures:** date them, resolve the line-104 contradiction,
  re-measure the RRBS row before release.
- **Open — does the model under-reserve on short reads?** σ ≈ 1.55 at 25 bp against 1.4375.
  §8's cell establishes a baseline but cannot fail, and the arithmetic says a gated one would
  not fire either. Genuinely unguarded; needs its own measurement to resolve.
- **Open — is the Linux baseline close enough to be useful?** If the first runs come back near
  0.85, the band works but the promise assertion has 15 points of slack. Nothing to decide until
  it is measured.
- **CRITICAL, needs Felix — what does §0.4 do to this phase?** The `--memory` promise is
  measurably broken at 25 bp / `--cores 2` (7 reps, two harnesses, 1.02–1.06× budget), and
  `CHANGELOG.md:30` claims the opposite. Three ways forward:
  1. **Scope the claim and ship Phase 3.** Qualify the CHANGELOG to the measured read lengths,
     make the 25 bp cell gating at a worsening-band, file the refit as its own issue. Phase 3
     still delivers the docs, the guards and the CI job; it just stops claiming a promise that
     does not hold everywhere.
  2. **Hold Phase 3 and refit σ first.** σ is read-length dependent and the model treats it as a
     constant; the honest fix is to charge the record spine per record rather than per byte.
     That is a Phase 4-sized change with its own calibration matrix, and it would reopen the
     coefficients Phase 3 promises not to touch.
  3. **Ship Phase 3, then refit.** Same as 1 but with the refit scheduled rather than filed.
  My recommendation is **1**, because every other deliverable in the phase is independent of the
  refit and the CHANGELOG line has not shipped yet — but the choice is a release-content
  decision, not mine.

## Implementation notes — PR 1, docs only (2026-08-17)

Branch `docs/clumpify-memory-arithmetic-439` off `dev` at `fdade7e`. `git diff --stat`:
`clumpy.md` and `clump-only.md` only — 34 insertions, 10 deletions, no formatter noise on
neighbouring lines. Uncommitted pending review.

**Built:** all eleven corrections in §3, plus the six additions — the core span (2–32, with the
`--cores 64` figures given explicitly since they exceed the 1 GiB default), definitions for
`static` / `10/11` / `workers`, the margin's documented home ("which is why the predicted peak
the binary prints is always a little under the `--memory` you gave it"), the dated-savings
caution, the `--clump_only` exception at line 110, and line 119 **qualified by read length rather
than stripped**. In `clump-only.md`: line 76's parity claim amended with the c5 threshold, and a
new `## Memory` section placed **after** the "Composes with:" list so it does not split it.

**Two published figures were wrong and verification caught them.** This is the whole reason
Validation 9 exists:

1. **The quoted error text invented a byte count.** I wrote "the derived bin budget is 917504
   bytes"; the real output at `--clump_only --fastqc --cores 16 --memory 700M` is **309806**.
   The 917504 came from an illustrative preview written during the D1 decision and was never a
   real invocation. Corrected to the verbatim output.
2. **The plain-mode figure was a reviewer's, not mine.** I published "~85 MiB for 10 M × 51 bp
   paired at `--cores 4`" from reviewer A's measurement; my own run of that cell gives
   **104.2 MiB**. Corrected to 104 MiB. The surrounding claim ("around the 100 MB mark") holds
   either way, which is why the error was easy to miss.

**Verified against the binary at `fdade7e`** — every figure now published: clumpify floors 281
(c2) / 311 (c6) / 590 (c32) / 933 (c64) plain and 334 (c2) / 1012 (c32) / 1356 (c64) with
`--fastqc`; `clump_only` floors 277 (c1–c4) / 454 (c32) and the 775 MiB c16-`--fastqc` retry
advice; bin sizes 9 MiB at `1G/c8` and 44 MiB at `4G/c8`; plain-mode 104.2 MiB. All match.

**Docs build:** `npm run build` clean, 32 pages. The CI `docs-build` job's assertions replicated
locally and all pass — 33 HTML / 31 OG images, no undersized PNGs, no `katex-error`, no `<del>`,
no stray maths outside `clumpy.md`, and the **positive control matched** (clumpy's display maths
found), so the null results are meaningful rather than vacuous. Rendered output spot-checked for
every figure and both new sections.

## Implementation notes — PR 2, code and guards (2026-08-17)

Branch `feat/clumpify-guards-439` off `dev` at `4f7b907`. **PR [#451](https://github.com/FelixKrueger/TrimGalore/pull/451)**,
commits `1b2579e` + `44dc732`. 7 files, 211 insertions.

**Built:** all of steps 1–10 as planned. 790 tests (was 734), fmt and clippy clean.

**The perturbation matrix was run, not assumed** — each applied alone, the guards' verdicts
recorded:

| perturbation | 4a pin table | absolute floor pin | 4b sweep |
|---|---|---|---|
| σ 23→20 | **FAIL** | **FAIL** | pass |
| k 64→32 | **FAIL** | **FAIL** | pass |
| static 224→512 MiB | **FAIL** | **FAIL** | pass |
| FastQC 24→0 MiB/thread | **FAIL** | ok | pass |
| cap 16→1024 | **FAIL** | ok | pass |
| `MIN_BIN_BYTES` **lowered** to 256 KiB | ok | **FAIL** | pass |

Every perturbation is caught by at least one guard; 4b passes throughout, which is the point of
splitting them. The last row is reviewer A's §4a gap (b) — previously caught by nothing.

Also proved able to fail: the channel-depth pin (work 1→2 → FAIL) and the banner tripwire
(divisor `1024²`→`1000²` → FAIL, which the old label-only assertion passed).

**Iteration log:**

- `#1` `pub(crate) fn new` alone left `ClumpLayout::new` dead in the non-test build, which
  fails CI's `-D warnings`. Its only caller is inside `parallel.rs`'s `#[cfg(test)]` block, so
  `#[cfg(test)] pub(crate)` is the honest form rather than an `#[allow]`.
- `#2` **Self-inflicted loss.** Reverting a perturbation with `git checkout src/clump.rs` while
  the baseline was still *uncommitted* destroyed every `clump.rs` edit. Re-applied from context;
  the other five files were untouched. The rule that follows: never revert a perturbation via
  git unless the baseline is committed — and commit before starting a perturbation sweep.
- `#3` All six pin values in the first draft were placeholders. Computed properly, five of six
  were wrong (e.g. the `--clump_only` shape 31,744,594 → 27,453,626). Never let an
  illustrative number reach a committed assertion.

**Recorded deviation:** 4d (deleting `resolve_layout_predicted_peak_fits_budget`) is **not** in
this PR. See PR 3.

## Implementation notes — PR 3, CI guard (2026-08-17)

Branch `ci/clumpify-memory-guard-439` off **PR 2's tip**, not off `dev`, so 4d's replacement
exists in the same history. **PR [#452](https://github.com/FelixKrueger/TrimGalore/pull/452)**,
commits `4e47b53` + `992235e`. 4 files, 273 insertions.

**Built:** `examples/mk_clump_fixture.rs` (SplitMix64, no new dependency, streaming so generator
memory is flat in `pairs`), the `clumpify-memory` job, `.github/clumpify-memory-baseline.tsv`,
and 4d's deletion plus the parent plan's note struck and its stale 19,884,107 figure corrected.

**Iteration log:**

- `#1` Clippy caught `write_with_newline` in the generator — exactly the category reviewer A
  predicted for a hand-rolled example, since `--all-targets` covers `examples/`.
- `#2` **The parse failed on a *good* log.** `[0-9]\+` is GNU-only; BSD `sed` does not support
  it, so every field came back empty and the guard fired on everything. It would have worked on
  the Linux runner, but a pattern that cannot be tested locally cannot be trusted — rewritten
  as `[0-9][0-9]*`, which is portable, then verified against three synthetic logs: a good one
  parses, the multi-pair one correctly takes the *first* of two banners, and a truncated one
  fires the empty-field guard.
- `#3` The margin assertion's `± 1 MiB` had only 94 KiB of real slack, because the printed
  prediction is truncated to whole MiB (reviewer A's Optional 11). Widened to 2 MiB, which
  still catches margin *removal* — a ~93 MiB move — by ~45×.
- `#4` Committed `plans/…/PLAN.md` by accident; `plans/` is deliberately untracked in this repo.
  Amended it out.

**Honest limits of what shipped.** The baseline file is in `BOOTSTRAP` mode: it records ratios
and warns rather than enforcing the band. The marker prints on every run so it cannot sit
unnoticed, and Validation 15 removes it once two or three runs have landed. Until then the
drift half of D5 is inert and only the absolute promise gates. The 25 bp cell is recorded
rather than gated, per D7 and #457.

### First runner results (run 32065004600) — Validation 11 and 15

The job passes, and its numbers change two things the plan assumed:

| cell | runner peak/pred | peak/budget | macOS/arm64 |
|---|---|---|---|
| `1G_c4` | **0.770** | 0.700 | 0.939 |
| `1G_c4_fastqc` | 0.697 | 0.633 | 0.724 |
| `25bp_1G_c4` | **0.941** | 0.854 | 1.07 |
| `multipair_c4` | 0.792 | 0.719 | — |

1. **The platform offset is ~17 points, not the ~9 §5.3 assumed** (0.770 vs 0.939 at the gating
   cell). So the Linux arithmetic in §5.3 was optimistic in the plan's own favour: at 0.77 the
   lever-A delta of 72–104 MiB lands at 788–820 MiB against a 930 MiB gate, further under than
   the 894 MiB the plan predicted. That strengthens the case for D5's band over the absolute
   assertion, which was the reason for D5 in the first place.
2. **#457 does not reproduce on the runner at `--cores 4`** — 0.941, inside both bounds. The
   25 bp penalty is present and similar in size (~17 points against 51 bp on the same platform,
   matching macOS's ~13), but Linux starts low enough to absorb it. `--cores 2`, the worst macOS
   cell, is unmeasured on Linux. #457 has been commented with this rather than left implying a
   general defect.
3. **No cross-pair retention.** `multipair_c4` at 0.792 sits beside single-pair 0.770, so the
   worry that pair 2's pool would sit on pair 1's retained pages — the mechanism that forced
   FastQC to be an addend — does not show at two pairs.
4. Free space in `RUNNER_TEMP` is 88 GiB, so §5.4's 8 GiB assertion has ample margin and is
   cheap insurance rather than a constraint.

**Iteration log, continued:**

- `#5` **The harness-capability probe failed on its first real run, silently, in 42 ms** — the
  exact failure mode it exists to prevent. `/usr/bin/time` is present on `ubuntu-latest` but its
  `-v` is unrecognised, and `2>&1 | grep -q` sent the usage message *into the pipe* where
  `grep -q` discarded it. Now: capture stderr to a file, test it, print what the binary actually
  said, install GNU time, retry once, then fail with the captured output. Each of the three step
  scripts is also extracted from the YAML and `bash -n`-checked locally now.
- `#6` PR #448's `Validate vs Perl TrimGalore` failed on a transient conda error. Confirmed
  pre-existing-unrelated by checking that the same job passed on `dev` 30 minutes earlier, then
  re-ran the failed job; green.

## Issues filed (2026-08-17)

| # | Title |
|---|---|
| [#453](https://github.com/FelixKrueger/TrimGalore/issues/453) | `--clumpify` writes plain output on plain input, having just refused `--dont_gzip` as pointless |
| [#454](https://github.com/FelixKrueger/TrimGalore/issues/454) | Documented memory ceiling (~100 MB) is ~20× out on long reads, and "not input file size" is wrong |
| [#455](https://github.com/FelixKrueger/TrimGalore/issues/455) | EOF bin flush discards a failed send, so a dead worker silently truncates output |
| [#456](https://github.com/FelixKrueger/TrimGalore/issues/456) | `--clump_only` should print its predicted peak, like `--clumpify` does |
| [#457](https://github.com/FelixKrueger/TrimGalore/issues/457) | **`--clumpify` exceeds the `--memory` budget on short reads**: σ is fitted per byte, not per record |
| [#458](https://github.com/FelixKrueger/TrimGalore/issues/458) | Record `--memory 2G`/`--cores 4` as a deliberate CI blind spot, with its named next lever |

#453 ships the `*.gz` symlink workaround and cites #245 so nobody re-lands the reverted
always-gzip behaviour. #454 reports by input size rather than merging two reviewers' different
fixtures, and names all three batch constants. #457 is cross-referenced from the CHANGELOG and
from the CI job's 25 bp cell.

## Self-Review

**What changed, and the lesson that actually holds.** v2's Self-Review taught that "nobody
guarded the output encoding". That lesson was drawn from a non-finding and is withdrawn. The
lesson the evidence supports is narrower and more useful: **every false measurement in this arc
came from an instrument that could not report its own failure**, and the fix is procedural —
hold the records byte-identical when comparing conditions, interleave reps rather than blocking
them, check exit status alongside the value, and assert free disk before measuring memory,
because disk pressure presents as a memory result. Five incidents now support it: the fixture
cycled 227×, the inherited `ru_maxrss`, the banner asserting `(gzip level 1)` over ASCII text,
my own disk-starved reps, and a review agent that fluently reviewed the wrong document at the
wrong commit.

**Logic.** The one ordering hazard is named and handled. D5 resolves the platform hole by
measuring rather than predicting, which is the same move that fixed v1's fixture problem.

**Efficiency.** PR 3's costs are now summed for the cells it actually has, and per-cell
truncation is required rather than optional.

**Edge cases.** Named: the banner's lost unit token, `--cores 1` having no remedy under D1,
`--cores` not moving the floor below 5, sub-MiB budgets, the FastQC cap boundary, the
multi-pair double banner, a parseable log from a failed run, and short reads as a known
under-charging axis that this phase does not close.

**Remaining risks.** (a) **Round 2 was single-reviewed** (§0.3), so v3's new material —
D5's band, §5.4, the §8 reclassification — has had one reviewer and me. (b) All ratios are
macOS; the Linux baseline is unmeasured, which D5 accommodates but does not remove.
(c) 4a pins six cells, so an edit preserving all six passes — a tripwire, not a proof.
(d) 4c catches a helper edit, not a re-inlined literal. (e) The short-read axis is unguarded
and the plan says so rather than implying coverage. (f) The stale savings figures stay stale
until re-measured; dating them is honest, not a fix.
