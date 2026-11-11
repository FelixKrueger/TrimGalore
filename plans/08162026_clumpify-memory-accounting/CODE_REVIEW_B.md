# Code review B — Phase 1 of #439 (`--clumpify` memory accounting)

**Reviewer:** B (independent; did not read `CODE_REVIEW_A.md`)
**Target:** uncommitted change on `dev`, base `2ab0918`
**Files:** `src/clump.rs`, `src/clump_only.rs`, `src/main.rs`, `src/parallel.rs`,
`tests/integration_preflight_tripwire.rs`, `CHANGELOG.md` (+49 / −24)
**Verdict:** the code is correct and the gating is right. Two findings should be
resolved before this ships: the CHANGELOG claims a bound the change does not
deliver, and lever B's memory effect inverts sign at 150 bp — the most common
Illumina read length, and one that was never measured.

## Verification (run, not assumed)

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | clean, exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 |
| `cargo test --release` | **729 passed, 0 failed**, exit 0, across 16 targets |

Agrees with the caller's record. The clippy and test runs were cache hits
(sub-second `Finished`), so I confirmed the cache covers the current source:
`cargo build --release -v` reports `Fresh trim-galore v2.3.0` against the current
source mtimes, and the built binary emits the new banner —
`clumpify: 32 bins × 9 MiB; predicted peak ≈ 1023 MiB (gzip level 1)` — so the
artifacts under test contain the change.

---

## High

### H1 — The CHANGELOG asserts a bound Phase 1 does not deliver

`CHANGELOG.md:8` — "**`--clumpify` now keeps its peak memory inside the
`--memory` budget it prints**".

It does not. The plan's own Phase-1 measurement (`PLAN.md:253`) records clumpy
`2G/c4` at peak 2077.2 MiB, **ratio 1.014**. The ratio is defined by the spike as
peak RSS ÷ banner prediction with ≤ 1.0 as the promise
(`SPIKE_memory-levers.md:27`), so 1.014 means the peak is ~29 MiB **over** the
figure the banner printed. `PLAN.md:284-287` says as much: "Phase 2 remains
necessary."

Failure scenario: a user reads the entry as a fix for #439, sets
`--memory 2G --cores 4` under a 2 GiB cgroup or an `sbatch --mem=2G`, and is
OOM-killed by the same overshoot the entry says is gone.

Recommend describing what Phase 1 did rather than what Phase 2 will do, e.g.:

> - **`--clumpify` peak memory reduced by ~14%** ([#439](…)) by holding fewer
>   bins in flight and sizing output reserves from each batch's own bytes; the
>   banner's figures are now labelled MiB.

Also drop "which is what they always were" — self-justifying historical register,
the class the repo's changelog rule exists to stop. Otherwise the entry is
correctly one sentence in the file's bullet format, matching the #429/#428
neighbours.

**Secondary, in the plan rather than the code:** `PLAN.md:286` reads 1.014 as
"inside the budget but with ~1% of margin". Under the ratio's own definition
1.014 is 1.4% *outside*, with no margin at all. Worth correcting so Phase 2 is
not planned against an inverted reading of its starting point.

### H2 — Lever B's sign inverts with read length; the inverting case was never measured

`src/parallel.rs:352-353` (paired), `:1074` (single-end).

The reserve changed from a flat 300 B/record to `Σ estimated_record_bytes` =
`id + seq + qual + 5` (`clump.rs:262-265`). In clumpy mode a bin flushes when
`raw_bytes >= bin_byte_budget`, and `raw_bytes` accumulates *the same*
`estimated_record_bytes` (`parallel.rs:841`, `:898`, `:905`). So per flushed bin:

```
new reserve = 1.0 · B          old reserve = B · 300 / rec_bytes
break-even at rec_bytes = 300
```

Measured against this repo's own fixtures:

| fixture | id | read | rec_bytes | new ÷ old |
|---|---|---|---|---|
| `test_files/SRR24766921_RRBS_R1.fastq.gz` | 54 | 51 bp | 161 | 0.54× (−46%) |
| `test_files/BS-seq_10K_R1.fastq.gz` | 54 | 65 bp | 189 | 0.63× (−37%) |
| `10K_150bp.fastq.gz` | 33 | 150 bp | 338 | **1.13× (+13%)** |
| 2×250 bp MiSeq | 33 | 250 bp | 538 | **1.79× (+79%)** |

Every measurement behind this change is on the short side of the break-even. The
spike's fixture is the 64 bp pair (`SPIKE:19-21`) and the spike itself flagged
this exact risk in its Limitations: *"on 150 bp reads `len * 300` under-reserves
rather than over-reserves, so the sign of this lever's effect may invert"*
(`SPIKE:132-134`). `PLAN.md:235-238` justifies shipping B **ungated** on
correctness grounds and does not record the inversion; the read-length axis is
validation row 11, deferred to Phase 2.

The reserve is also multiplied by the in-flight count. Capacity travels with a
moved `Vec` and nothing calls `shrink_to_fit`, so each reserved buffer is held
through the result channel (`cores * 2` slots, `parallel.rs:130`, `:959`) plus up
to `cores` in workers plus whatever sits in the main thread's `pending` map. At
`--clumpify --cores 8 --memory 4G` on 2×150 bp (B = 66 MiB) that is ~66 MiB of
reserved output capacity per result × up to 24 in flight ≈ **1.6 GiB**, against
the sizing formula's charge of `0.5 · B · cores` = 264 MiB. Note lever A shrank
the *work* channel and left this larger, ungated pool at depth `cores * 2`.

Failure scenario: 2×150 bp at `--memory 2G --cores 4` — the cell Phase 1 leaves
1.4% over (H1). Raising the output-reserve term ~13% moves that cell further
over, at the most common Illumina read length, with no measurement at that
length in the record.

Recommend either:

1. Measure one 150 bp cell before this lands (validation row 11, pulled forward —
   it is one run and it decides whether B ships ungated), **or**
2. Reserve a fraction on the gzip path only: the buffer receives *compressed*
   bytes there and deflate L1 on FASTQ lands at ~25–30%, so the full text size is
   a ~3–4× over-reserve in the default mode. Keep the full text size for the
   `--dont_gzip` path, where it is exact and where the old 300 B constant actually
   under-reserved (see L1).

Not a data-integrity issue either way — see L1 for why exceeding the reserve is
only ever a realloc.

---

## Medium

### M1 — Docs now contradict the binary's own banner

`docs/src/content/docs/performance/clumpy.md:102`:

> So with `--cores 8 --memory 1G` you get **32 bins × 12 MB**; with
> `--cores 8 --memory 4G` you get 32 bins × 66 MB. **The binary prints those
> resolved values at startup.**

The binary now prints `32 bins × 9 MiB`. Both the unit and the first number
disagree: `4 × 512 MiB / (5·32 + 7·8) = 9.48 MiB`, so "12" was already wrong
arithmetic; "66" is correct read as MiB. The sentence pointing readers at the
banner is what makes the mismatch user-visible, and this change is what created
it. `PLAN.md:213-216` already says the docs' formula, floor range, worked values
and unit labels "are pure arithmetic and can be fixed now" — only the published
savings figures need Phase 2. Two-number, two-label fix.

### M2 — The `resolve_layout` derivation comment is the accounting's only written record, and both levers move the code away from it

`src/clump.rs:190-199` still states:

- "Worker input batches in flight: up to `cores × bin_byte_budget`" — with depth
  1 the real count is **2 per worker** (one queued, one in process).
- "Worker output Vecs growing during gzip-encode: roughly `0.5 ×` input size at
  L1 … the 0.5 × output factor is the L1 worst case" — the *reserve* is now
  deterministically `1.0 ×` the input text regardless of what L1 achieves.

Neither is newly wrong (pre-change it was 3 per worker and 1.74·B of reserve at
64 bp, so both levers move *toward* the comment). But this docstring is the only
place the model is written down and it is what Phase 2 will calibrate against, so
leaving it describing the pre-lever code is how the next revision inherits a
wrong starting point. One edit stating the current in-flight counts.

### M3 — The `cores` remedy reintroduces the drift it removed, inside the test helper

`src/parallel.rs:1615-1617` — `ClumpLayout::new(16, 1024 * 1024, 2)`. Four of the
five call sites pass `cores = 4` to `run_*_parallel` (`:1671`, `:1764`, `:1823`,
`:1893`); only `:1935` passes 2. Harmless today because nothing on those paths
reads `layout.cores`, but the whole point of part 3 was that a layout must not
carry a core count the run does not use. A future assertion on
`predicted_peak_bytes()` in these tests would compute for 2 cores against a
4-core run — silently, and with the parameter gone there is no longer a call site
to notice it at. Use `4`, or have the helper take `cores`.

---

## Low

### L1 — The reserve is a heuristic, not a bound (two ways to exceed it; consequence is a realloc)

Answering the question directly: **nothing worse than a realloc.** Both sinks
grow on demand — plain mode writes into `&mut Vec<u8>`, gzip mode into
`GzEncoder<&mut Vec<u8>>` — so an under-reserve cannot truncate or corrupt
output. Two real ways to exceed it:

1. **Plain path with `--rename`.** `trimmer.rs:264` and `:273` append
   `:clip5:<seq>` / `:clip3:<seq>` (7 + n chars) to the id while seq and qual each
   lose n, so the net per record is `14 − (n1 + n2)`. `--dont_gzip --rename
   --clip_R1 2 --three_prime_clip_R1 2` writes ~10 B/record **more** than
   reserved. (`estimated_record_bytes` is evaluated pre-trim, so ordinary
   trimming can only shrink.)
2. **Gzip path on incompressible, untrimmed input:** 10-byte header + 8-byte
   trailer + ~5 B per stored block exceeds the text size. Not realistic for FASTQ.

Worth one line of acknowledgement rather than a fix, with one memory caveat: a
realloc of a nearly-full N-byte buffer allocates 2N and copies, i.e. a transient
3N for that buffer — in clumpy mode N ≈ B/2 per mate, the shape of spike this
plan exists to bound. Conversely, credit where due: at 150 bp in **plain** mode
the old 300 B constant under-reserved (338 needed) and reallocated every batch;
the new reserve removes that. So B is a genuine improvement on the `--dont_gzip`
path and a trade on the gzip path (H2) — not a strict improvement, as the plan's
implementation note implies.

Passthrough and unpaired writes go to their own `Vec::new()` buffers
(`parallel.rs:354-359`), so neither can overflow `buf_r1`/`buf_r2`.

### L2 — The gate is right, but duplicated and unobservable

`parallel.rs:117` and `:943`. It compiles because `Option::is_some` takes
`&self`; the `Copy` on `ClumpLayout` (`clump.rs:136`) is incidental — the later
`if let Some(layout) = clump_layout` (`:190`, `:992`) copies rather than moves,
but a move would also have been fine since nothing reads it afterwards.

No path can mismatch. The predicate is the *same* one that selects
`read_pairs_clumpy` / `read_single_clumpy` in the same function, so depth 1 ⟺
clumpy reader and depth 2 ⟺ round-robin, with no third arm. `--passthrough` is
covered twice over: `--clumpify` + `--passthrough` is rejected at `cli.rs:940`,
and `resolve_clump_layout` (`main.rs:178-209`) returns `None` unless
`cli.clumpify` — including the below-floor fallback, which therefore also gets
depth 2 alongside its round-robin reader. `--cores 1` cannot reach either path
(`--clumpify` requires `--cores >= 2`). Both work channels were patched; both
result channels correctly stayed at `cores * 2` (but see H2 on their reserve
footprint).

Two nits: the magic 1/2 and the identical comment are duplicated at both sites — a
named helper or two consts would make them provably identical; and queue depth is
not observable from any test, so the gate's only guard is the measurement harness.

### L3 — MiB relabelling is complete for what it touched; two adjacent surfaces remain

All six figures across the five strings divide by `1024 * 1024`, so all six are
genuinely MiB: `main.rs:205`, `:206`; `clump_only.rs:285`, `:426`, `:764`,
`:1003`. Verified end-to-end by running the binary.

Remaining user-facing, both one-liners in `--help`:

- `cli.rs:233-236` documents `--memory` as "e.g. `1G`, `512M`, `8G`" without
  saying the suffixes are 1024-based, though `parse_memory_size`
  (`clump.rs:244-247`) is KiB/MiB/GiB. The flag that produces the now-MiB banner
  does not state its own base.
- `cli.rs:345` still says "~5 MB RAM" in `--cores`' help.

Nothing else in `src/` carries a mislabelled memory figure: the other `MB` hits
(`fastq.rs:219`, `:925`, `bam.rs:47`, `parallel.rs:30`) are internal comments.

### L4 — The new assertion bites, but guards only one of the five strings

`tests/integration_preflight_tripwire.rs:598-602`.

**Would it have failed before?** Yes, on both clauses: the old banner read
`… × 9 MB; predicted peak ≈ …`, so `"MiB; predicted peak ≈"` was absent *and*
`" MB;"` was present.

**Spurious failure?** Low risk today — no `eprintln!` in `src/` emits `" MB;"`,
and `cli.rs`' help text is not on this path. But `!contains(" MB;")` scans the
whole of stderr, so it is a repo-wide tripwire that an unrelated future message
could trip; and folding both clauses into one `assert!` means the panic message
cannot say which one failed.

**Coverage gap:** the four `clump_only.rs` labels are asserted by nothing —
`--clump_only` never runs in this test, and the negative clause cannot reach
strings that were never printed.

### L5 — Part 3 is behaviour-neutral; confirmed

Only two call sites ever existed, and both passed the sizing `cores`:
`main.rs:201` + `:206` (`cli.cores` in both) and `clump.rs:481-482` (loop `cores`
in both). The sole struct literal outside `clump.rs` was `small_clump_layout()`,
now converted. So `predicted_peak_bytes()` returns the same value for every
existing caller. Two structural notes:

- `n_bins` and `bin_byte_budget` remain `pub` while `cores` is private, so the
  prediction can still be desynchronised by mutating a `mut` layout. The
  asymmetry weakens the invariant the change was added to protect; private fields
  with accessors, or `#[non_exhaustive]`, would close it.
- `ClumpLayout::new`'s doc says "for callers that bypass `resolve_layout`" —
  there are none in production. It exists for the test helper, so `pub(crate)`
  would say that honestly. (Also: adding a private field to a `pub` struct whose
  fields were all public is a semver break for any downstream literal
  construction — irrelevant unless someone depends on the library API.)

### L6 — Pre-existing, out of scope: EOF bin flush swallows a failed send

`parallel.rs:848-852` — the `break` on a failed `flush_bin` exits only the `for
bin in bins.iter_mut()` loop, then falls through to sending poison pills, so the
remaining non-empty bins are dropped rather than reported. A send only fails when
a worker's receiver is gone, in which case the main thread has an error of its
own, so this is not currently reachable as silent data loss. Unrelated to this
change (depth does not affect receiver liveness) — noted only because the repo
treats output completeness as an invariant.

---

## Comment style

Judged against "one line, current fact, never reference the state being fixed".

| Comment | Verdict |
|---|---|
| `parallel.rs:116`, `:942` — "A clumpy batch is a whole bin, so a **shallower** queue caps resident bins." | **Flag.** "Shallower" is comparative — it only means something to a reader who knows the old depth. Suggest "…so the queue depth caps resident bins." |
| `parallel.rs:351` — "Each buffer holds one mate of this batch, so reserve that mate's own bytes." | Good: one line, current fact. |
| `parallel.rs:1073` — "The buffer holds this batch, so reserve the batch's own bytes." | Good. |
| `clump.rs:140-141` — "Worker count this layout was sized for. The prediction is only valid for this count, so it is **stored rather than passed in again**." | **Flag** (in the diff, not in the caller's list of three). The trailing clause narrates the removed parameter — the state being fixed — and it takes two lines where one does. Suggest "Worker count this layout was sized for; the prediction is only valid for it." |
| `tests:597` — "The banner divides by 1024², so both figures are MiB." | Good: states the invariant on the thing it constrains. |

The `ClumpLayout::new` doc ("for callers that bypass `resolve_layout`") and the
`predicted_peak_bytes` doc ("so the startup banner cannot disagree with the
sizing") are both current-state and fine in register; see L5 on the former's
accuracy.

## CHANGELOG house style

One sentence, bold lead, issue link, matching the file's bullet format — compliant
(the neighbouring #434 entry is far longer, but #429 and #428 set the pattern this
one follows). The problem is factual accuracy, not length: see H1.

---

## Summary of recommendations, in order

1. **H1** — reword the CHANGELOG so it does not claim the bound; correct
   `PLAN.md:286`'s reading of ratio 1.014.
2. **H2** — measure one 150 bp cell, or reserve a fraction on the gzip path.
3. **M1** — fix `clumpy.md:102` (9 MiB, 66 MiB).
4. **M2** — update the `resolve_layout` derivation comment to the current
   in-flight counts.
5. **M3** — `small_clump_layout()` should say `cores = 4`.
6. **L2/L3/L4** — dedupe the queue-depth constant; state the `--memory` suffix
   base in `--help`; consider splitting the two assertion clauses.
7. **Comment style** — two edits: "shallower" → "the", and trim the `cores` field
   doc to one line.

Nothing here blocks on correctness of the shipped code paths: the gate is right,
part 3 is behaviour-neutral, output byte-identity holds by construction (flush
order comes from fill order and round-robin, not channel depth; capacity changes
cannot alter bytes), and 729 tests pass clean.
