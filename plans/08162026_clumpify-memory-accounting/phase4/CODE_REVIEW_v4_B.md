# Code review — #457 memory refit (Reviewer B)

Target: working-tree diff on `perf/clumpify-spine-accounting`, based on `origin/dev` @ `3fa4de3`.
Spec: `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v4 (§5, §6, §8, §13).

**The tree changed under me mid-review.** I read the diff at 226 insertions / 79 deletions
(`git diff | shasum -a 256` → `733141f5cefa2759`, 2026-08-20T06:14Z) and had already read an earlier
state at 224/77. Three of my findings were fixed by whoever made those edits while I was working; they
are kept below, marked **RESOLVED**, with what changed, so this report can be compared against
Reviewer A's line by line. Everything else was re-verified against the later state.

**No source edits were made by me.** M2/L1/L2 were the three I would otherwise have applied directly,
and the concurrent-edit collision above is exactly why I did not: the change is uncommitted and a second
reviewer is working the same working tree, so an edit from me would silently rewrite the diff both
reviews describe.

## Summary

The arithmetic is right. I rebuilt the model from first principles in an independent script (`n_bins`,
`static`, the `10/11` margin, `33·n_bins + 92·workers`, the `MIN_BIN_BYTES` bail and
`clumpify_min_memory_bytes`'s rounding) and every one of the six pinned budgets, both range tests, the
FastQC-delta expectation, all nine floor figures across the two docs pages, both quoted error messages
and both `--memory` restore budgets match — and each of the six pins fails under `(23, 64)`, so none is
vacuous. The renamed `channel_depths` test is not weaker than the one it replaced, and both new
`parallel.rs` tests are non-vacuous by construction.

What I found is mostly *around* the change rather than in it: the release notes still tell users short
reads are uncovered, and the new floor warning fires on the very fixture the floor was measured from —
so it will fire on both CI cells this PR exists to gate.

### Verified independently

| check | result |
|---|---|
| `cargo fmt --all -- --check` | clean (both tree states) |
| `cargo clippy --all-targets --release -- -D warnings` | clean; forced a real re-check by touching the five edited files rather than trusting a cache hit |
| `cargo test --release` | **805 passed, 0 failed**, 18 test binaries |
| `docs && npm run build` | clean, 32 pages, no MDX/KaTeX breakage from the `$$` edit |
| `ci.yml` via `docs/node_modules/js-yaml` | parses; `clumpify-memory` job present |
| `bash -n` on the extracted `Measure and gate` `run:` | clean |
| enforced band vs TSV prose | `+0.065` / `-0.072` in both — match |
| six pinned budgets | all six recomputed and matching; all six move under `(23, 64)` |
| docs floor table, layout examples, both error messages | reproduced from `target/release/trim_galore` |

Positive confirmations, since these were the plausible failure modes:

- **Every bin drain goes through `take()`** — `parallel.rs:825`, `:1246`, and the four `clump_only.rs`
  flush helpers (`:207`, `:222`, `:747`, `:759`). The EOF drains use the same helpers, so no path
  empties a bin behind the mean's back.
- **No stale mean crosses a boundary.** All six `*Bin::with_budget` call sites construct their `Vec`
  inside the per-run dispatcher, so a multi-pair run and a multi-input `--clump_only` run both start
  every file at `prev_mean_bytes == 0`.
- **No per-record cost.** The new branch sits inside the `self.keys.capacity() == 0` cycle-start guard;
  `take()` adds one integer division per flush. `bytes_this` was already computed per record.
- **The uBAM branch returns the same quantity.** `bam.rs:979` builds the id as `format!("@{id_body}")`,
  so `estimated_record_bytes` counts the same bytes on both arms of `sanity_check_any`. All eleven
  `sanity_check_any` call sites propagate with `?`; no error is dropped.
- **The single-end path is covered by a paired-end calibration.** Spine per text byte is `76/r` single
  against `74/r` paired — 2.7% worse, far inside the 1.435 the refit bought.
- **`estimated_record_bytes` and the plan's `r` are the same statistic** (`id + seq + qual + 5`), so the
  floor's *units* are right. See M1 for the statistic that is not.
- **Fixed-length input reserves exactly, with zero doublings.** `mean` from the previous cycle equals `r`
  exactly and the cycle ends at `ceil(B/r)` records, which is what `div_ceil` reserves.

---

## High

### H1 — the release notes still say this defect is unfixed, and nothing records the fix

`CHANGELOG.md:69-79`, in the **Unreleased** section, closes the #439 entry with:

> Very short reads are not yet covered: the per-bin coefficient is fitted per byte rather than per
> record, so a 25 bp library packs more records into the same bin bytes and can exceed the budget at
> `--cores 2` ([#457](...)).

That is now false, and it sits in the same unreleased cycle as the change that falsifies it — it would
ship as release notes describing behaviour the binary no longer has. Separately the diff adds no
changelog line at all, though Step 5.4 asked for one and §13's "Not done, and why" does not mention
dropping it. Three user-visible facts have no entry:

- bin sizes fall ~30.3% at **every** budget (verified: 1G/c4 18.13 → 12.62 MiB; 4G/c32 11.22 → 7.81 MiB);
- `--clump_only`'s hard floor rose (c16 352 → 398 MiB, c32 454 → 544 MiB), so budgets that ran now error;
- `.fq.gz` bytes move again — the third time this cycle, which is what 5.4 wanted said once.

Replace the stale clause with one sentence in the file's existing voice, e.g. *"`--clumpify` now holds
peak memory inside the printed budget on 25 bp reads too, at the cost of ~30% smaller bins at every
budget and higher `--memory` floors ([#457](...))."* The synced page
(`docs/.../reference/changelog.md`) does not yet carry the Unreleased block, so only `CHANGELOG.md`
needs editing.

---

## Medium

### M1 — the floor warning fires on the fixture the floor was measured from, including both newly-gated CI cells

The threshold and the estimator are different statistics:

- `CLUMPIFY_CALIBRATED_FLOOR_BYTES = 97` is `ceil(96.642)`, and 96.642 B is the **mean** record size of
  the byte-equalised 25 bp fixture — mean header 41.6 chars, `41.64 + 50 + 5 = 96.64`.
- `warn_below_clumpify_floor` compares the **first** record. `examples/mk_clump_fixture` writes the
  record index into the header (`@SIM:1:FCX:1:15:6329:{i} 1:N:0:ATCCGA`), so record 0 carries the
  shortest header in the file — 35 chars, giving **90 B** at 25 bp, 6.6 B under the fixture's own mean.

Measured, not inferred:

```
$ ./target/release/trim_galore --clumpify --cores 2 --memory 1G -o out ci25.fq.gz
WARNING: the first record is 90 bytes, below the 97 bytes the --memory reservation was
calibrated at, so peak memory may exceed the printed prediction. ...
```

Consequences:

1. **The CI job's log will contradict its own gate.** `25bp_1G_c4` and `25bp_1G_c2` are the two cells
   promoted to `gating 1` by this PR — the deliverable (4c.8). Both will print "peak memory may exceed
   the printed prediction" in the same log where the job asserts `peak ≤ prediction`. Nothing fails, but
   the next person reading that log has to work out which of the two to believe.
2. **The docs understate who warns.** `clumpy.md` says "25 bp with a terse header, or already-trimmed
   input". Solving `id + 2L + 5 < 97` gives the real boundary: at 25 bp *any* header under 42 chars
   warns, which includes ordinary `@SIM:…`- and SRA-style headers, not only terse ones. At 36 bp only
   headers under 20 chars warn; by 40 bp effectively none.
3. `main.rs:173-175` says the floor is "25 bp with a 41-byte header", which computes to **96**, not 97.
   The 41 is the fixture's *mean* header length (41.64), so the comment's own arithmetic cannot reach
   the constant it explains.

Recommendation (a judgement call, so not applied): measure both sides the same way — either set the
threshold from the first-record size of the calibration fixture (90 B) rather than its mean, or average
the first N records inside `sanity_check_with` so the runtime statistic matches the calibrated one. A
floor warning that fires at the calibration point is the one thing it must not do.

### M2 — RESOLVED during the review: two displaced doc comments

In the state I first read, the inserted items sat between two docstrings and the functions they
describe, leaving `resolve_clump_layout` (`main.rs`) and `estimated_record_bytes` (`clump.rs`)
undocumented while their text documented `CLUMPIFY_CALIBRATED_FLOOR_BYTES` and the `const _` assert.
Both are now correctly attached: the const + `warn_below_clumpify_floor` were moved above
`resolve_clump_layout`'s docstring, and the `const _` assert above `estimated_record_bytes`'s. Nothing
further needed — recorded because Reviewer A may have seen the earlier state.

### M3 — the docs' copy of the new warning is not what the binary prints

`docs/.../performance/clumpy.md:120-123` quotes:

```
WARNING: the first record is 68 bytes, below the 97 bytes the --memory reservation
         was calibrated at, so peak memory may exceed the printed prediction.
```

The binary emits one unindented line, and a second sentence this block drops entirely: *"Short reads
with terse headers and already-trimmed input reach here; raise --memory if the run is near its limit."*
The 9-space continuation indent is invented — the `--memory too small` block a few lines below
reproduces its own two `eprintln!`s faithfully, indent included, so a reader will reasonably take this
one as verbatim too. Quote both sentences, or mark the block as abridged.

### M4 — `1G_c2` keeps `gating 0` and lost the reason it had

`ci.yml:444-446`:

```
# c2 gives the largest bin, so it runs closest to its prediction.
measure 1G_c2        0 ...
measure 25bp_1G_c2   1 ...
```

The old comment stated a decision ("--cores 2 is too tight to gate"); the new one states a mechanism and
then sits above an ungated *loose* cell followed by a gated *tight* one — the inverse of what the
sentence implies. Post-refit, `1G_c2` lands near 0.72 of prediction (pre-refit 0.879 through the same
`1/f` scaling this PR uses elsewhere), so gating it is free coverage on the largest-bin configuration.
Either promote it to `gating 1` or record why it stays out; `.github/clumpify-memory-baseline.tsv` has a
whole "Measured but NOT banded, with the reason" section establishing that convention.

Same class, smaller: the `2G_c4` comment ("A 2 GiB pool needs ~2.3x this fixture to saturate") is a
pre-refit figure. The pool shrank ~30%, so the multiple is now ~1.6x — still above 1, so `gating 0`
still holds, but the number is stale in a file this PR edits, and it is the number that justifies the
cell staying ungated.

### M5 — the bootstrap instructions don't ask for the measurement the band's own rule requires

The TSV states the admission rule as `width >= 3x the observed half-range`, from *a double-digit
sample*. The new header block justifies `0.065` with "the half-ranges scale down by the same 1.435, so
0.065 is still 3.4x and 4.3x them" — a prediction, since no post-refit sample exists. That scaling
assumes the whole spread rides on `B`. The static/allocator component does not shrink, and because the
ratio's denominator (`budget × 10/11`) is unchanged while the pool term shrinks by 1.435, the static
component's *share* of the variance grows. So the post-refit half-range is plausibly wider than
`range/1.435`, and 3.4x could come in under 3x.

The `# BOOTSTRAP` note tells the next person to "Record the ratios from this PR's own Linux run, replace
both rows, and delete this line" — centres only. Add the range to that instruction and re-check the 3x
side before the marker comes out; otherwise the band arms at a width whose noise margin was never
observed, on a job that already has a flakiness history. (The file's other claim — that deleting the
marker without re-recording reds the band rather than passing silently — is correct: ~0.62 against a
0.788 baseline trips the `-0.072` side.)

### M6 — the residual exposure the change does not remove: one slot short costs a full `Vec` doubling

For fixed-length input the new reservation is exact (see the positive confirmations). For
variable-length input — the case the change exists for, and which §6.2 notes no fixture here can see —
`predicted` now lands within ±1 slot of the true count, so **any** downward drift in mean record size
overflows the `reserve_exact` and `Vec::push` grows the allocation to `2 × predicted`
(`RawVec::grow_amortized` → `max(2·cap, required)`). For symmetric variability that is roughly half of
all cycles, and the transient during the copy holds old + new capacity at once.

Scale, at 51 bp pairs and 1G/c4 (`B` = 12.62 MiB, `r` = 284 B, count ≈ 46,600, spine 148 B/pair):
6.9 MiB of spine per bin becomes 13.8 MiB, ~20.7 MiB transient. Across 16 bins that is 110 MiB of extra
spine against a 930 MiB prediction — **+0.118 of prediction, wider than the 0.065 drift band this PR
narrows to.** It cannot show up in CI, because every fixture is fixed-length.

This is not a regression: the old code under-reserved whenever the first record was longer than the
cycle mean, also about half the time, and could double repeatedly from a much lower base. But the change
gets nearly all the remaining benefit for a small headroom that removes the doubling outright:

```rust
let predicted = self.budget.div_ceil(est.max(1)).max(1);
let predicted = predicted + predicted / 16 + 1;   // absorb mean drift; a doubling costs 2x
```

6.25% of a bin's spine against a 100% overshoot on half the cycles. Worth doing in a PR about peak
memory; flagged rather than applied because the constant is a judgement call.

### M7 — three of the four patched sites have no test

Step 6.1 patched four `reserve_exact` sites; the two new tests cover exactly one
(`parallel.rs::PairedBin`). `parallel.rs::SingleBin`, `clump_only.rs::SingleBin` and
`clump_only.rs::PairedBin` are hand-duplicated copies of the same six lines, so the failure this guards
against — a copy that reads `bytes_this` where it should read `est` — is invisible to the suite and, per
§6.2, invisible to every fixture. I read all four and they are correct today; the gap is that nothing
keeps them that way. `clump_only` is the mode where bin size *is* the product. Two more tests in the
shape of the existing ones close it; `clump_only`'s structs are private, so they need an in-file
`mod tests`.

### M8 — the docs don't name the configuration that newly falls back at the default budget

Verified against the binary (`--cores 50` engages; `--cores 51` warns and falls back):

| at the default `--memory 1G` | before | after |
|---|---:|---:|
| max `--cores` keeping `--clumpify` engaged | 72 | **50** |
| same, with `--fastqc` | 33 | **23** |
| max `--cores` for `--clump_only` at 512M | 41 | **28** |

`clumpy.md` gives the floor table and says c64 is above the 1 GiB default, but the newly-crossed
threshold is `--fastqc` at `--cores 24–33`, which used to clump at the default budget and now silently
falls back to plain mode. One clause in the below-floor section covers it.

For the plan record: §5 row 12 predicts 52 and 29 where the binary gives 50 and 28. Nothing shipped
quotes those numbers, so no docs defect — but Validation 12's expectations need correcting before
someone "verifies" them against the binary and concludes the code is wrong.

---

## Low

- **L1 — RESOLVED during the review.** `prev_mean_bytes`'s doc comment said "Mean bytes per record" on
  both `PairedBin`s, where the value is bytes per pair. Now reads "Mean bytes per reserved slot over the
  last flushed cycle", which is accurate for all four structs (a paired slot is one pair).
- **L2 — RESOLVED during the review.** `clump_only.rs`'s "See `PairedBin::push`" pointed at a local
  `PairedBin::push` that carries no such comment; it now reads `parallel::PairedBin::push`.
- **L3** `first_record_bytes == 0` (`main.rs:181`) is unreachable: `estimated_record_bytes` is
  `id + seq + qual + 5` and `id` has already been checked to start with `@`, so the minimum is 6. It
  reads as a "no sample" sentinel for a state the function cannot produce — the residue of the
  unreachable `Ok(0)` §13 records fighting.
- **L4** `FastqReader::sanity_check` / `sanity_check_with` are `pub` on a crate published to crates.io
  (`trim-galore`), so `Result<()>` → `Result<usize>` is source-breaking for a library consumer. Layering
  too: the entry-side I/O guard now calls `crate::clump::estimated_record_bytes`, i.e. `fastq` depends
  on `clump` for a clumpify-specific quantity. Returning the record — as the plan phrased it — would
  keep the coupling out and serve the next caller that wants a different metric.
- **L5** `const _: () = assert!(2 * size_of::<FastqRecord>() + size_of::<MinimizerKey>() == 148)` fails
  *compilation* on any 32-bit target (`String` is 12 B → 76). Nothing builds 32-bit today and the repo
  claims no such support, so this is a note rather than a defect: `#[cfg(target_pointer_width = "64")]`
  keeps the pin without making the crate unbuildable elsewhere.
- **L6** The TSV's stats table still has `at +0.09` as its column header over two rows that now report
  at `+0.065`.
- **L7** `clumpy.md`'s new section states two model extrapolations as observations: "at 150 bp a run
  peaks around 0.61 of its printed prediction, against 0.94 at 25 bp". 0.61 is the pre-refit 0.775 put
  through the `1/f` scaling §13 explicitly calls "a conservative approximation, not a verified
  identity", and the measured post-change 25 bp values were 0.821–0.923, not 0.94. Both err
  conservatively (they overstate peak), so this is wording — "should peak around", not "peaks".
- **L8** The floor warning prints before the clumpify banner and before the below-floor fallback
  warning, so `--clumpify --memory 100M` on short-read input says peak "may exceed the printed
  prediction" and then falls back to plain mode, printing no prediction and reserving nothing (verified).
  Gating the warning on a resolved layout means moving it past `resolve_clump_layout`, which is called
  from three places — hence the current position; worth one line saying so.
- **L9** `assert_eq!(bin.keys.capacity(), from_mean)` over-specifies `reserve_exact`, which is documented
  to allocate *at least* the request. It passes today, and the `capacity() < from_short` assertion above
  it is the load-bearing one — this is only a note about which line can rot.
- **L10** `modes/clump-only.md` gets the new floors but not the "raise `--memory` to recover the bin
  size" remedy `clumpy.md` gained, and `clumpy.md`'s restore budgets are computed at `workers = cores`.
  At `workers = 1` they nearly coincide (1G→1363M at c4 is identical; 2G/c8 is 2832M against the printed
  2834M), so a sentence pointing at the new section is enough.
- **L11** `parallel.rs:900` is 92 chars in a docstring wrapped at 82 — an artefact of the concurrent
  edit to the two lines above it. `cargo fmt` does not reflow comments, so it stays green.

---

## Fixes applied

None by me — see the note at the top. M2, L1 and L2 were fixed in the tree by a concurrent editor while
this review was in progress; L6 is the remaining one-liner I would otherwise have applied.
