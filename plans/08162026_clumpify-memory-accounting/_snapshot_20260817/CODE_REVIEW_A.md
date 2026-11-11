# Code Review A — Phase 1 of #439 (clumpify memory accounting)

**Reviewer:** A (independent; did not read `CODE_REVIEW_B.md`)
**Target:** uncommitted working-tree diff on `dev`, base `2ab0918`
**Touched:** `src/clump.rs`, `src/clump_only.rs`, `src/main.rs`, `src/parallel.rs`,
`tests/integration_preflight_tripwire.rs`, `CHANGELOG.md` (49 insertions, 24 deletions)
**Plan:** `PLAN.md` v3, "Implementation notes — Phase 1"
**Date:** 2026-08-16

## Verdict

**The code is correct.** All four parts do what the plan says, the gate is sound, and I found
no path that mixes the wrong queue depth with the wrong batching mode. Every finding below is
about *claims* (changelog, docstrings, comments) and about one test that pins less than it says
it does. Nothing here blocks the diff on correctness grounds; H1 does block it as written,
because the changelog states a promise that the plan's own measurements say is not yet kept.

## Verification (run, not assumed)

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 |
| `cargo test --release` | 729 passed, 0 failed (455+15+13+16+11+1+8+57+14+21+3+57+19+39) |
| `se_trim_clumpify` present in the run | `test se_trim_clumpify ... ok` |

Agrees with the caller's numbers exactly. Both cargo invocations were cache hits (no compile),
which is expected — the artefacts were already current for this source.

Live banner from the built binary (`--clumpify --cores 2`, default `--memory 1G`):

```
clumpify: 16 bins × 21 MiB; predicted peak ≈ 1023 MiB (gzip level 1)
```

## Correctness confirmations (the caller's specific probes)

**The gate is right, and it compiles for a duller reason than suspected.**
`parallel.rs:117` reads `clump_layout.is_some()`; `Option::is_some(&self)` *borrows*, so no move
happens and `Copy` is not what makes it compile. `ClumpLayout` *is* `Copy` (`clump.rs:136`), which
is what lets the later `s.spawn(move || … if let Some(layout) = clump_layout …)` at
`parallel.rs:190` / `:992` take it by copy. Either way, both sites read the **same** `Option`
that selects the reader (`:117`/`:190` paired, `:943`/`:992` single-end), so depth and batching
cannot disagree — there is no third path.

The two ways `clump_layout` can be `None` while `--clumpify` was requested both keep depth 2 and
round-robin, correctly: the below-floor fallback (`main.rs:185-199`) returns `None`, and
`--passthrough` never coexists with `--clumpify` at all — `cli.rs:939-941` is a hard bail, not a
downgrade. So the passthrough path always keeps its jitter absorption, as intended. Related
guards worth recording: `--clumpify` requires `--cores >= 2` (`cli.rs:732`) and rejects
`--dont_gzip` (`cli.rs:737`), so **every clumpy run is a gzip run** — which matters for M2/L5
below.

**Depth 1 cannot deadlock.** The wait graph is acyclic: reader → work channels (depth 1) →
workers → result channel (`cores * 2`) → main. Main never sends on a work channel, so shrinking
work-channel capacity can only shorten the reader's lead, never introduce a cycle. Head-of-line
blocking on a busy worker is pre-existing (`next_worker` advances only on a successful send);
depth 1 makes it bite sooner. Measured free at 4 and 8 cores in the spike; untested above 8, but
the reader is the bottleneck there, so the risk is low.

**The reserve can never be a correctness problem.** `Vec<u8>` grows on demand and `write_all`
into a `Vec` cannot fail, so exceeding the reserve costs one realloc. Bounds, per branch:

- **Plain (`gzip == false`)**: the reserve is *exact*. `FastqRecord::write_to` emits
  `id + 1 + seq + 1 + 2 + qual + 1` (`fastq.rs:87-94`), byte-for-byte the definition of
  `estimated_record_bytes` (`clump.rs:262-265`), and trimming only shrinks `seq`/`qual` while
  filtered records go to `buf_up_*` / are dropped. One exception: `--rename` appends
  `:clip5:<seq>` / `:clip3:<seq>` to the id (`trimmer.rs:264`, `:273`), which nets about +14 B
  per record against the shrunk seq/qual — so the bound is "exact except under `--rename`".
  Don't describe it as exact anywhere.
- **Gzip**: always a large over-reserve (the buffer receives compressed bytes). Deflate framing
  cannot realistically exceed it — worst case is stored blocks at +5 B/64 KB plus an 18 B
  header/trailer, ~110 B on a 1.2 MB batch, against a reserve of the full raw text.
- Unpaired and passthrough writes go to their own `Vec::new()` buffers (`parallel.rs:354-359`),
  so nothing else lands in `buf_r1`/`buf_r2`. `reads_r1.len() == reads_r2.len()` always (batches
  are pushed in lockstep), so the "unequal mates" the plan mentions is a *length* difference, not
  a count one.
- Batches are never empty (both readers flush only non-empty bins/batches), so no
  `with_capacity(0)` pathology.

**Part 3 changes no behaviour.** There are exactly two `predicted_peak_bytes` callers —
`main.rs:206`, whose layout came from `clump::resolve_layout(memory_bytes, cli.cores)` two lines
above at `:201`, and `clump.rs:482`, whose layout came from `resolve_layout(budget, cores)` in
the same loop iteration — and both previously passed the same `cores` the layout was sized with.
Construction sites are also exactly two: `clump.rs:227` and the test helper at
`parallel.rs:1616`. Nothing else builds a `ClumpLayout` by literal, in `src`, `tests` or
`examples`.

**All five relabelled figures really are MiB.** Divisors verified individually:
`main.rs:205`, `main.rs:206`, `clump_only.rs:285`, `:426`, `:764`, `:1003` — every one is
`/ (1024 * 1024)`. The already-correct neighbours (`main.rs:192` fallback warning,
`clump.rs:217-224` bail message) also say MiB. Sweeping `src` for every remaining
`/(1024*1024)` in user-facing output found no other mislabelled computed figure; the only
remaining `MB` in user-facing text is prose (`cli.rs:345`, "~5 MB RAM" in `--cores` help — an
approximation where the distinction is noise) and the docs (L2).

## Findings

### H1 — High: the CHANGELOG claims the bound that Phase 2 is still needed for

`CHANGELOG.md:8-11` opens with **"`--clumpify` now keeps its peak memory inside the `--memory`
budget it prints"**. The plan's own Phase-1 measurement table (`PLAN.md:250-255`) records
`2G/c4` at **ratio 1.014** — still *above* the printed budget — and `PLAN.md:284-287` says in
as many words that "`2G/c4` still sits at 1.014 … which is not a margin. Phase 2 remains
necessary." The model in `predicted_peak_bytes` also still charges only `1.75 · B` per core
(`(5·n_bins + 7·cores)/4` = `1.25·n_bins + 1.75·cores`), i.e. roughly one bin's worth of
spine per worker, while lever A leaves *two* bins per worker resident — so the formula continues
to under-charge by design, deliberately, until Phase 2 calibrates it.

Concrete failure: a user reads this entry, submits `--clumpify --memory 2G --cores 4` under a
2 GiB cgroup or scheduler allocation on the strength of it, and is OOM-killed. That is #439
happening again, with the changelog as the source of the wrong expectation.

This is also the specific failure mode the repo has been bitten by before (memory:
"corrective copy contains the defect it corrects") — a fix for dishonest memory figures should
not ship with a dishonest claim about itself.

**Recommend** rewording to what Phase 1 measurably did, and leaving the bound to Phase 2. For
example: "**`--clumpify` now holds fewer bins in flight, cutting its peak RSS by 13–17%**
([#439](…)); the banner's figures are also now labelled MiB, which is what they always were."
One sentence, matches the file's bold-headline-plus-issue-link format, and stays true when the
issue is still open.

### M1 — Medium: the model docstring now describes a quantity the diff changed

`clump.rs:190-193`, item 2 of the peak-RSS derivation, still reads "Worker input batches in
flight: **up to `cores × bin_byte_budget`** of text + matching spine." With depth 2 that was
actually up to `3 × cores × B` (two queued plus one in the worker); after lever A it is
`2 × cores × B`. The diff moved this number without touching the sentence that states it.

Why it matters more than a typo: `PLAN.md:73-90` makes this docstring the basis of Phase 2's
`k ≈ (work-channel depth + 1) · σ` calibration. Leaving the prose saying "depth 0" is how a
calibration gets fitted against the wrong mechanism.

**Recommend** the one-line fix: "Worker input batches in flight: up to
`2 × cores × bin_byte_budget` (one queued, one in the worker)."

### M2 — Medium: lever B is a data-dependent trade, and is described as a reduction

Measured mean `estimated_record_bytes` on real fixtures in this repo:

| fixture | mean B/record | vs the old 300 |
|---|---|---|
| `10K_150bp.fastq.gz` (150 bp) | **337.9** | **+13% reserve** |
| `test_files/BS-seq_10K_R1.fastq.gz` (~75 bp) | 192.6 | −36% reserve |

Break-even is around 125 bp for a 45-char header (`id + 2L + 5 = 300`). So on the single most
common Illumina configuration the new reserve is *larger* than the constant it replaced. The
spike anticipated exactly this (`SPIKE_memory-levers.md:132-134`: "on 150 bp reads `len * 300`
under-reserves rather than over-reserves, so the sign of this lever's effect may invert"), but
`PLAN.md:234-237` and `CHANGELOG.md:10` both present it as part of the memory fix.

It is still the right change — it is a genuine correctness fix, exact in the plain branch, and
it removes the doubling cascade that long reads used to trigger — but it should be described as
tidiness, not as a memory lever (which is what the spike recommended, `SPIKE:111-113`).

**Recommend** (a) drop it from the memory claim in the changelog and the plan's "Built" list,
and (b) optionally reserve `sum / 2` in the *gzip* branch only. Half the raw text is above every
compressed size we have measured (0.17× on my long-read fixture, ~0.3× typical), is smaller than
both the old and the new value at every read length, and removes the read-length sensitivity
entirely. The plain branch should keep the full `sum`, where it is an exact bound.

### M3 — Medium: two of the added comments break the repo's comment rules

`clump.rs:140-141`:

```rust
/// Worker count this layout was sized for. The prediction is only valid for
/// this count, so it is stored rather than passed in again.
```

"rather than passed in again" is a reference to the parameter this diff removed — written for a
reader who saw the old signature. That is the exact pattern `CLAUDE.md` calls out ("rather than
replace it", flagged twice by @pinin4fjords). Suggest: `/// Worker count this layout was sized
for; the prediction is only valid for that count.`

`parallel.rs:1073`, `// The buffer holds this batch, so reserve the batch's own bytes.` — this
narrates the line beneath it and says nothing the code doesn't. The paired twin at `:351`
("Each buffer holds one mate of this batch, so reserve that mate's own bytes") earns its keep
because *per-mate* is the non-obvious part; the single-end one has no such content. Suggest
deleting it, or reducing both to the fact that the reserve is the batch's own byte count.

The test comment at `integration_preflight_tripwire.rs:597` ("The banner divides by 1024², so
both figures are MiB") states the evidence rather than the fact; "Both banner figures are MiB"
would do. Minor.

### L1 — Low: the new assertion pins one label, and its comment claims two

`tests/integration_preflight_tripwire.rs:598-602`. Given the live banner
`clumpify: 16 bins × 21 MiB; predicted peak ≈ 1023 MiB (gzip level 1)`:

- `contains("MiB; predicted peak ≈")` pins the **first** figure. It would have failed before this
  change (the old banner read `… × 21 MB; predicted peak …`), so the test does bite — confirmed
  by inspection of the old format string, not just by assumption.
- `!contains(" MB;")` adds nothing the first clause doesn't already catch, and it **cannot** catch
  a revert of the second label: `predicted peak ≈ 1023 MB (gzip level 1)` is followed by `" ("`,
  never `" ;"`. So the comment's "both figures" is not what the assertion checks.

Spurious-failure risk: low. No other stderr line in this path prints `" MB;"`; the fallback
warning uses MiB and FastQC is not involved.

**Recommend** `contains("MiB; predicted peak ≈") && contains("MiB (gzip level")`, which pins
both and drops the redundant negative.

### L2 — Low: the docs now contradict the banner they quote

`docs/src/content/docs/performance/clumpy.md:102`: "So with `--cores 8 --memory 1G` you get
**32 bins × 12 MB**; with `--cores 8 --memory 4G` you get 32 bins × 66 MB. **The binary prints
those resolved values at startup.**" As of this diff it prints MiB. The figures are already
MiB-derived, so this is a two-word edit with no arithmetic attached — distinct from the docs
*arithmetic* that `PLAN.md:214-216` defers. Since the sentence explicitly claims to quote the
banner, it should move in the same change.

### L3 — Low: `ClumpLayout::new` is public API with no production caller and no invariants

`clump.rs:148-154`. Its doc says "for callers that bypass `resolve_layout`"; the only such caller
is the test helper at `parallel.rs:1616`. It also skips both invariants `resolve_layout` enforces
— `cores >= 1` and `bin_byte_budget >= MIN_BIN_BYTES` — so `ClumpLayout::new(16, 0, 0)` is
constructible and `predicted_peak_bytes()` will happily return `STATIC_OVERHEAD` for it. Nothing
divides by `cores`, so there is no panic today; the cost is a public constructor that can mint a
layout `resolve_layout` would have refused.

**Recommend** `pub(crate)`, or keep it public and say in the doc that it validates nothing. The
constructor itself was the right remedy for the E0063 hazard — this is only about its reach.

### L4 — Low: storing `cores` removes one drift surface and creates a smaller one

`clump_only.rs:275`, `:413`, `:748`, `:905` all size layouts with `resolve_layout(budget,
cores.max(1))` from the user's `--cores`, while that path is single-threaded in v1 (its own
comment at `:277` says so). The stored `cores` therefore does not describe the run it belongs to.
Harmless today — no `--clump_only` banner prints a prediction — but the field is now the obvious
thing for a future caller to trust. Concrete future failure: Phase 3 adds a predicted-peak line
to the `clump-only` banner, and it prints a figure sized for N workers that path never spawns.

**Recommend** either naming the field `sizing_cores`, or having `clump_only` pass `1` and keep
`n_bins` sizing separate — a decision that pairs naturally with the deferred item 5.

### L5 — Low, informational: measurement that answers the spike's open question 4

The spike left "why lever B's benefit does not scale with bin size" unexplained
(`SPIKE:73-77`), and `PLAN.md:213-215` still lists it open. Two runs of the built binary on a
20,000 × 10 kb-read fixture at `--cores 4` (`/usr/bin/time -l`, macOS/arm64) settle it:

| branch | reserve per buffer | bytes actually written | peak RSS |
|---|---|---|---|
| plain text (`.fastq` in → `.fq` out) | 82 MB | 82 MB (touched) | **802 MiB** |
| gzip (`.fastq.gz` in → `.fq.gz` out) | 82 MB | ~14 MB (0.17× ratio) | **508 MiB** |

Identical record shape and identical uncompressed size (400,608,890 B, 20,000 records, same
header form); the two fixtures differ only in sequence content, one being motif-repetitive and
the other diverse. Both runs therefore reserve the same 82 MB per buffer, and the gzip run — the
one with an *extra* background decompressor thread — peaks ~294 MiB **lower**. Had the untouched
capacity counted, the two would have peaked alike; instead the difference is precisely the
bytes-written difference, ~68 MB × ~4 in-flight buffers that were reserved and never faulted in.
So the over-reserve is lazily committed at large sizes (mimalloc mmaps it;
the OS faults pages in on write), and it *is* resident at small sizes, where the allocation comes
out of an already-committed segment. That is the whole of finding 4: −47 MiB at `1G` (buffers
~1.2 MB, resident) versus −3.6 MiB at `2G` (bigger buffers, lazier).

Two consequences for Phase 2, both free:

1. **Lever B is close to memory-neutral in clumpy mode**, because every clumpy run is a gzip run
   (`cli.rs:737`) and clumpy buffers are bin-sized, not batch-sized. Reinforces M2.
2. **The next term of the same kind is the result channel**, `mpsc::sync_channel(cores * 2)` at
   `parallel.rs:130` / `:959`, plus the `pending` `BTreeMap`. Those hold *written* (therefore
   resident) compressed buffers, up to ~3 × cores of them, and are ungated and untouched by this
   diff. If the `2G/c4` residual needs another structural cut, that is where I would look before
   reaching for buffer recycling.

Also worth knowing, and pre-existing rather than caused by this diff: 802 MiB peak in plain mode
at `--cores 4` is a sharp counterexample to `docs/…/performance/threading.md:61`, "even
worst-case at the saturation point (`--cores 8`), the process never exceeds ~100 MB". `BATCH_SIZE`
is 4096 *records* regardless of length, so plain-mode memory scales linearly with read length
with no bound at all. Not this diff's problem; relevant to #439's theme, and to the ONT/PacBio
fixture the clumpy docs already cite.

## Summary table

| # | Sev | Where | One line |
|---|---|---|---|
| H1 | High | `CHANGELOG.md:8-11` | Claims the `--memory` bound the plan measures as still breached at `2G/c4` (1.014) |
| M1 | Medium | `clump.rs:190-193` | Model docstring still says in-flight is `cores × B`; lever A made it `2 × cores × B` |
| M2 | Medium | `parallel.rs:352-353`, `:1074` | Reserve is a data-dependent trade (+13% at 150 bp), described as a reduction |
| M3 | Medium | `clump.rs:140-141`, `parallel.rs:1073` | Comment references the removed parameter; another narrates its own code |
| L1 | Low | `integration_preflight_tripwire.rs:598` | Pins the first unit label only; the negative clause is redundant |
| L2 | Low | `docs/…/performance/clumpy.md:102` | Docs quote the banner as `MB` after the banner became `MiB` |
| L3 | Low | `clump.rs:148-154` | Public constructor, no production caller, skips both `resolve_layout` invariants |
| L4 | Low | `clump_only.rs:275` et al. | Stored `cores` does not describe that path's single-threaded run |
| L5 | info | — | Measured: large over-reserves are lazily committed — answers the spike's open question 4 |

No Critical findings. No fixes applied (two reviewers share this tree; recommendations only).
