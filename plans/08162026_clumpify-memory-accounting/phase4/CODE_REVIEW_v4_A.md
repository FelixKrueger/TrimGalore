# Code Review — Phase 4 (#457) refit, Reviewer A

**Target:** working-tree diff on `perf/clumpify-spine-accounting`, based on `origin/dev` @ `3fa4de3`.
**Spec:** `plans/08162026_clumpify-memory-accounting/phase4/PLAN.md` v4, §5 / §6 / §8 / §13.
**Reviewed:** 2026-08-20. Independent of Reviewer B (their report was not read).

**Verdict: the arithmetic is right and the hot path is correct. Two blocking gaps are outside `src/`** —
Step 5.4 (CHANGELOG) was not implemented and not declared as skipped, and the new floor warning is
absent for every pair after the first in a multi-pair run. Neither is caught by any gate in the repo.

---

## Summary

| area | result |
|---|---|
| `prev_mean_bytes` hot path, all four sites | **correct**; denominator is bytes-per-pair against a bytes-per-pair estimator; `take()` is the only emptying path |
| `sanity_check_with` signature change, 11 call sites | **correct**; no error path lost, no warning, uBAM branch comparable |
| `warn_below_clumpify_floor` reachability | **gap** — first input only; silent for pairs 2..N |
| compile-time `148` pin | **correct on 64-bit**; hard compile error on 32-bit targets |
| integer arithmetic / overflow / inverse | **correct**; proven by an existing test at `cores` 65536 × `budget` u64::MAX under dev-profile overflow checks |
| six pinned budgets, floor, ranges | **all reproduce exactly** from the formula, recomputed independently |
| CI workflow + baseline file | YAML and bash valid; band genuinely disarmed; **5 stale `+0.09` references** and one flipped verdict remain in the baseline prose |
| docs figures | **every re-derivable figure verified**, and the docs *corrected* the plan on one (`2834M`, not §6d's `2833M`); three unverified/inconsistent items below |
| comment discipline | four defects, all fixed directly (below) |

Verification after my fixes: `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release
-- -D warnings` clean (forced a real recompile with `touch`, not a cache hit — 2.53 s of `Checking`);
**805 tests pass, 0 failures**.

---

## Critical

None.

---

## High

### H1 — Step 5.4 (CHANGELOG) was not implemented, and not listed under "Not done, and why"

`CHANGELOG.md:76-79`, in the **Unreleased** section, still reads:

> Very short reads are not yet covered: the per-bin coefficient is fitted per byte rather than per
> record, so a 25 bp library packs more records into the same bin bytes and can exceed the budget at
> `--cores 2` ([#457](https://github.com/FelixKrueger/TrimGalore/issues/457)).

This diff is the fix for #457, so that sentence is false in the tree it would ship from. Two things are
missing:

1. The stale sentence has to go (or become "now covered down to 97 B/record, warned below that").
2. A ~30% bin-size reduction at **every** budget including the default `--memory 1G` is user-visible —
   the docs' own caution box now says exactly that — and has no entry. The existing #439 entry's
   "`.fq.gz` output differs from earlier releases on identical input" line covers the byte-boundary half
   of Step 5.4's "cumulative about output bytes", but not the sizing change itself.

Step 5.4 is the only Step-5 sub-item absent from §13, and §13's "Not done, and why" names Validations
7/9/10/14 and the §10/§12 refit — not this. A skipped step that is neither done nor declared is the
failure mode the plan's own §13 structure exists to prevent.

`docs/src/content/docs/reference/changelog.md` does **not** yet carry the #439/#457 entry, so the fix is
confined to `CHANGELOG.md`. Recommend, not applied: the entry's wording is the maintainer's call, and
CLAUDE.md caps it at one sentence.

### H2 — `warn_below_clumpify_floor` is silent for pairs 2..N of a multi-pair run

`src/main.rs:294` is the only call, and it consumes only `sanity_check_any(&cli.input[0])`. Every
subsequent pair resolves its own clumpify layout and prints its own prediction, but is never checked.

Confirmed live — long first pair (BS-seq 10K), short second pair (68 B/record):

```
=== Pair 1 of 2 ===
clumpify: 16 bins × 15 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
=== Pair 2 of 2 ===
clumpify: 16 bins × 15 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
```

No warning at all, for a pair that is 29 B/record below the calibrated floor. The reverse ordering warns
*once, at the top*, for a run whose second pair is fine — so the signal is wrong in both directions.

The multi-pair loops already have the value in hand and discard it: `src/main.rs:809`, `:811`
(specialty/uBAM-out paired), `:1008`, `:1010` (paired trim), `:2150`, `:2152` (the shared multi-pair
driver), plus the single-end loops at `:1069` and `:2203`. Threading `warn_below_clumpify_floor(&cli, …)`
onto those closes it.

Recommend rather than applied: it changes user-visible stderr on a path that integration tests read, and
"warn per pair" vs "warn once for the worst pair" is a product choice.

---

## Medium

### M1 — `clumpy.md` quotes a warning the binary does not print

`docs/src/content/docs/performance/clumpy.md:120-123`:

```
WARNING: the first record is 68 bytes, below the 97 bytes the --memory reservation
         was calibrated at, so peak memory may exceed the printed prediction.
```

What the binary actually prints (captured, 25 bp with an SRA-style header, `--clumpify --paired
--cores 2 --memory 1G`) is one unwrapped ~290-character line, and it does not stop where the docs stop:

```
WARNING: the first record is 68 bytes, below the 97 bytes the --memory reservation was calibrated at, so peak memory may exceed the printed prediction. Short reads with terse headers and already-trimmed input reach here; raise --memory if the run is near its limit.
```

The `\` continuations in `src/main.rs:191-194` strip the newline *and* the following indentation, so the
9-space continuation the docs show never exists. §13 claims "the two error messages **captured from the
binary**" — the two pre-existing ones were; this new third one was not.

Note the sibling warning 10 lines further down in the same doc page *does* match its source, because
`resolve_clump_layout` emits it as two `eprintln!`s with a literal 9-space indent. So the house style for
a multi-line warning is two calls, and this one is the outlier — which makes "wrap the binary to match
the docs and the sibling" at least as defensible as "rewrite the docs to one long line". That choice is
why I did not apply either.

### M2 — the docs publish predicted ratios as observed behaviour

`clumpy.md:112`: "at 150 bp a run peaks around **0.61** of its printed prediction, against **0.94** at
25 bp".

Both numbers are §6c's `ratio new` column, which §5 row 8 and §6c itself label a *prediction* until
re-measured. §13 re-measured only the two 25 bp cells, and both came in **below** prediction:

| cell | §6c predicted | §13 measured | miss |
|---|---:|---|---:|
| 25 bp c8 2G | 0.939 | 0.923 / 0.909 / 0.907 | −1.6 to −3.2 pts |
| 25 bp c2 1G | 0.903 | 0.821 | −8.2 pts |

So the model runs optimistic-high by up to 8 points, and no 150 bp cell was re-measured at all. Quoting
0.61 as what a run does is stating an extrapolation from a relationship §13 explicitly records as "**not**
exactly linear … a conservative approximation, not a verified identity". Either attribute them ("the
model predicts ~0.6 at 150 bp") or use the measured 25 bp figure (~0.91) for the one cell that has one.

### M3 — the baseline file still quotes the retired `+0.09` width in five places, and one verdict flips

`.github/clumpify-memory-baseline.tsv`. The diff updated the two *banded* rows' annotations to `+0.065`
but left the rest of the file on the old width, including the column header immediately above them:

| line | text | problem |
|---|---|---|
| 45 | `#   cell … range   at +0.09` | header contradicts the two rows under it, which now say `at +0.065` |
| 51 | `multipair_c4 … 1.2x at +0.09; needs +0.121 for 3x` | at 0.065 with the half-range scaled by 1.435 it is ~2.3x |
| 58 | `1G_c4_fastqc … 1.2x at +0.09` | ~2.3x |
| 73 | `Lever A moves it ~+0.09 to +0.13` (`2G_c4`) | ~+0.063 to +0.091 post-refit, which now *straddles* the 0.065 width |
| 77 | `1G_c2 … 2.9x at +0.09 — just under the bar` | **verdict flips** — half-range 0.017 → 0.0118, so 0.065 is ~5.5x, well over the file's own 3x bar |
| 78 | `25bp_1G_c2 … 4.2x at +0.09 — eligible next pass` | ~4.3x; verdict unchanged, number stale |

Line 77 is the one that matters: this file is the operational reference a maintainer reads to decide what
to band next, and it now advises "just under the bar" about a cell that clears the bar comfortably. The
same 1.435 scaling the diff applies to the two banded cells applies here.

### M4 — `ci.yml`'s replacement comment no longer explains the `0` it sits above

```yaml
# c2 gives the largest bin, so it runs closest to its prediction.
measure 1G_c2        0 --clumpify --paired --cores 2 --memory 1G "$F51R1" "$F51R2"
measure 25bp_1G_c2   1 --clumpify --paired --cores 2 --memory 1G "$F25R1" "$F25R2"
```

The deleted comment ("`--cores 2` is too tight to gate; both cells record their ratio only") explained
both zeros. The replacement is a reason to *gate*, and the two cells now differ, so nothing in the file
says why `1G_c2` stays un-gated. Pre-refit it read 0.879 on Linux and post-refit ~0.70, so it looks
gate-able; if it is deliberately held back pending the Linux measurement, say that.

### M5 — nothing tests `warn_below_clumpify_floor` or the 97-byte threshold

`grep` finds no reference to `warn_below_clumpify_floor` or `CLUMPIFY_CALIBRATED_FLOOR_BYTES` outside
`src/main.rs`. Because it lives in the binary crate it is unreachable from `cargo test` unit tests, so
the only route is a `tests/integration_*.rs` stderr assertion — a pattern this repo already uses heavily
(`integration_partial_output.rs`, `integration_ubam.rs`).

The asymmetry is what makes this worth raising: the coefficient half of this change is pinned by six
exact assertions plus a `const` assert, and the user-facing half — the boundary at which the plan says
"the bound is extrapolation" — is pinned by nothing. A threshold that drifts to 96 or 98, or a guard
inverted to `>`, is invisible to CI.

---

## Low

### L1 — the `148` pin is 64-bit-only

`src/clump.rs`:

```rust
2 * std::mem::size_of::<FastqRecord>() + std::mem::size_of::<MinimizerKey>() == 148
```

`FastqRecord` is three `String`s; on a 32-bit target that is 36 B, giving `2·36 + 4 = 76`, and the
`const _: () = assert!(…)` becomes a **hard compile error** rather than a failed assumption. Every
shipped target is 64-bit (`release.yml` builds `x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`,
`aarch64-apple-darwin`), so this only bites `cargo install trim-galore` on i686/armv7 — hence Low. Fix
is one attribute:

```rust
#[cfg(target_pointer_width = "64")]
const _: () = assert!(…);
```

Not applied: gating the pin narrows the invariant's coverage, which is a deliberate call.

### L2 — one of four changed sites is tested

The two new tests cover `parallel::PairedBin` only. `parallel::SingleBin` and both bins in
`clump_only.rs` are hand-duplicated copies of the same eight lines — precisely the shape where a
copy-paste slip survives. I read all four and they are logically identical, so this is coverage, not a
defect. `bin_reservation_follows_the_previous_cycle_not_one_short_record` generalises to the other three
with a two-line helper each.

### L3 — the capacity assertion leans on an undocumented std guarantee

```rust
assert_eq!(bin.keys.capacity(), from_mean);
```

`Vec::reserve_exact` documents that "the allocator may give the collection more space than it requests".
It does not today (`RawVec` stores the requested `cap`), so the test passes and is a sharper check than
the `<` assertion above it. Worth knowing it is an implementation detail, not a contract.

### L4 — the floor warning fires before the prediction it refers to, and on a path with no prediction

Observed order:

```
Trim Galore v2.3.0
==================================================

WARNING: the first record is 68 bytes … may exceed the printed prediction.

Using 2 worker threads …
… 7 lines …
clumpify: 16 bins × 15 MiB; predicted peak ≈ 930 of 1024 MiB budget (gzip level 1)
```

Also: on the below-floor fallback path `resolve_clump_layout` returns `None` and no prediction is ever
printed, yet the floor warning has already promised the reader one. Moving the call to just after the
`clumpify:` line (and skipping it when the layout is `None`) fixes both. This interacts with H2 — a
per-pair call site would land it in the right place naturally.

### L5 — `first_record_bytes == 0` is unreachable

`estimated_record_bytes` is `id.len() + seq.len() + qual.len() + 5`, so it is ≥ 5 for any record, and
both arms of `sanity_check_any` return it for a record that exists (the empty-input case `bail!`s). The
guard reads as though `0` were a "size unknown" sentinel, which it never is.

### L6 — the uBAM estimate omits `--preserve-tags` bytes

`sanity_check_any` opens the BAM with `BamReader::open(path)`, which carries no preserved tags, so the
sampled size is the bare `@QNAME` form. A `--preserve-tags CB,UB` run's real records are larger. The
error is in the safe direction (a spurious warning, never a missing one), and the `@`-prefix convention
in `bam_record_to_fastq` means the quantity is otherwise directly comparable to the FASTQ branch's.

### L7 — `fastq.rs` now reaches into `clump.rs`

`FastqReader::sanity_check_with` calls `crate::clump::estimated_record_bytes`. The helper is a generic
"FASTQ text bytes for this record" computation with nothing clumpify-specific in it; `fastq.rs` is
arguably where it belongs, and moving it would let `clump.rs` re-export rather than have the reader
depend on the dispatcher.

### L8 — test doc comment carries provenance CLAUDE.md pushes to the commit message

```rust
/// A short record arriving first after a flush must not size the whole next
/// cycle. Every fixture in the memory calibration is fixed-length, so this
/// case is invisible to it — real libraries are variable-length after
/// trimming.
```

Sentence one is the assertion; sentences two and three are the §10 Finding 2 argument for why the test
exists. Committed tests are inside CLAUDE.md's brevity rule ("This covers code, tests and config that
gets committed"), and this is the "hard-won reasoning inflates the comment" pattern the file names.

---

## Fixes applied

All four are comment-only; zero behaviour change. `cargo fmt` clean, `cargo clippy --all-targets
--release -- -D warnings` clean on a forced recompile, 805 tests pass afterwards.

**1. `src/main.rs` — `resolve_clump_layout` had lost its doc comment to the new const.**
The new `CLUMPIFY_CALIBRATED_FLOOR_BYTES` was inserted *between* `resolve_clump_layout`'s five-line
`///` block and its signature, so the const carried a doc that opens "Resolve `--memory` into a
`(n_bins, bin_byte_budget)` layout when `--clumpify` is set, and emit a one-line startup notice…" and
the function had none. Moved the const and `warn_below_clumpify_floor` above the block and restored it
to the function.

**2. `src/clump.rs` — the same mistake, on a `pub` function.**
The `const _: () = assert!(…)` was inserted between `estimated_record_bytes`' two-line doc and its
signature, so a public function lost its rustdoc and `const _` gained a four-line doc that begins by
describing a different item. Reordered: the assert's own doc now sits on the assert, and
`estimated_record_bytes` has its doc back.

Both are invisible to `rustc`, `clippy` and `rustfmt` — doc comments on `const _` items are legal and
unwarned — which is why the same slip happened twice.

**3. `src/parallel.rs:895-897` — stale `PairedBin::push` docstring.**
It still said the reservation is extrapolated "from this record's text size and the bin's budget", which
is now the fallback rather than the rule. Now: "expected pair count from the previous cycle's mean pair
size and the bin's budget".

**4. `prev_mean_bytes` field doc was wrong on the two `PairedBin`s, and the cross-reference dangled.**
"Mean bytes per **record** over the last completed cycle" is bytes per *pair* on `PairedBin`
(`raw_bytes` accumulates `r1 + r2`, `keys.len()` counts pairs) — a factor-of-two misdescription of the
quantity the review brief specifically asked about. Replaced at all four sites with one wording that is
correct everywhere and matches the const assert's vocabulary: "Mean bytes per reserved slot over the
last flushed cycle; 0 initially." Also retargeted `clump_only.rs`'s `See PairedBin::push` to
`parallel::PairedBin::push`, since the local `PairedBin::push` it resolved to carries no explanation.

---

## Verified correct — recorded so it is not re-checked

**The `prev_mean_bytes` mechanism (review area 1).** Correct at all four sites, and they are genuinely
consistent now.

- The denominator is right. `PairedBin::raw_bytes` accumulates `estimated_record_bytes(r1) +
  estimated_record_bytes(r2)` per push and `keys.len()` counts pairs, so `raw_bytes / keys.len()` is
  bytes-per-*pair* — the same unit as `bytes_this`, and the same unit `budget.div_ceil(est)` needs to
  yield a slot count. `SingleBin`'s is bytes-per-record against a bytes-per-record `bytes_this`.
- No division by zero: guarded by `!self.keys.is_empty()`, and `estimated_record_bytes ≥ 5` means the
  quotient can never round to 0 either, so the `> 0` sentinel test in `push` cannot misfire.
- `take()` is the only path that empties a bin. All six flush sites go through it —
  `parallel.rs:825` and `:1246`, `clump_only.rs:207`, `:222`, `:747`, `:759` — and there is no `.clear()`
  or direct `mem::take` on any bin field outside `take()` itself. `std::mem::take` leaves a
  capacity-0 `Vec`, so `keys.capacity() == 0` correctly identifies the first push of the next cycle.
- First-cycle fallback preserves the previous behaviour exactly.
- `est.max(1)` also closes a latent `div_ceil(0)` panic that `parallel.rs` had (`div_ceil(bytes_this)`
  unguarded) and `clump_only.rs` did not. Unreachable in practice, but the four sites now agree.
- The residual asymmetry is inherent and not a defect: a mean that *over*-estimates (long cycle followed
  by a short one) under-reserves and admits one `Vec` doubling. The single-record estimator had the same
  exposure in both directions; the mean is strictly the more robust of the two.

**`sanity_check_with` (review area 2).** All three `bail!` paths survive; the match simply became the
tail expression. `usize` is not `#[must_use]`, so the ten call sites that discard the value with `?;`
are warning-free under `-D warnings` (confirmed by a forced recompile). The uBAM arm returns
`estimated_record_bytes` of the converted `FastqRecord`, and `bam_record_to_fastq` builds
`id = format!("@{}", name)`, matching the FASTQ reader's `@`-inclusive convention — so the two branches
return the same quantity, not one off by the sigil.

**Threshold applied to the right quantity (review area 3).** §6g's floor `r = 96.642 B` is bytes per
record (25 bp + a 41-byte header ⇒ `41 + 1 + 25 + 1 + 2 + 25 + 1 = 96`), and `sanity_check_any` returns
per-record bytes for R1. Strict `< 97` covers the whole calibrated range and excludes nothing above it.
Reachability: `--clumpify`, `--clump_only`, paired, single-end and uBAM input all pass through
`main.rs:293` before dispatch — confirmed live for `--clumpify --paired` and `--clump_only --paired`.
The only gap is multi-pair (H2).

**Integer arithmetic and the inverse (review area 5).** No new overflow exposure.
`floor_and_sizing_round_consistently_without_overflow` already sweeps `cores ∈ {1,2,3,4,8,16,17,32,1024,
65536}` × `fastqc ∈ {false,true}` × `budget ∈ {floor, floor+1, u64::MAX/2, u64::MAX}`, `Cargo.toml`
overrides nothing on the dev profile so `cargo test` runs with `overflow-checks = true`, and it passes
with the new constants. `predicted_peak_bytes` remains a sound inverse for a structural reason that does
not depend on the coefficients: `resolve_layout` floors `B = ⌊usable·16/d⌋` and the prediction computes
`⌊B·d/16⌋`, and `⌊x/d⌋·d ≤ x` for all `x`, so `predicted ≤ static + usable = budget·10/11` always. The
`saturating_mul`s upstream mean even `--memory 18446744073709551615` cannot push `B·d` past `u64::MAX`.

**All six pinned budgets and the floor.** Recomputed independently from the formula, not read off the
diff — every one matches:

| case | asserted | recomputed |
|---|---:|---:|
| `uniform(4)`, 1 GiB | 13_236_569 | 13_236_569 |
| `bin_cores 4, workers 1`, 1 GiB | 19_128_978 | 19_128_978 |
| `+ fastqc 16`, 1 GiB | 6_046_334 | 6_046_334 |
| `+ fastqc 17` (cap) | 6_046_334 | 6_046_334 |
| `uniform(2)`, 1 GiB | 16_657_256 | 16_657_256 |
| `uniform(32)`, 4 GiB | 8_191_148 | 8_191_148 |
| floor at `cores 2` | 296 MiB | 296 MiB |

Both range windows hold (4 GiB/c2 = 78.64 MiB ∈ [78, 80); 4 GiB/c8 = 31.25 MiB ∈ [31, 33)), the
in-comment denominators 712 / 1792 / 896 are right, and the FastQC-delta expectation
`24·4·MiB·16/896 = 1_797_558` is 1 byte off the true 1_797_559 — inside the test's `≤ 1024` tolerance.
§13's note that the plan's table was wrong by 1 on the floor (296, not 294/295) is confirmed.

**CI workflow and baseline file (review area 6).** `ci.yml` parses under `docs/node_modules/js-yaml`;
all four extracted `run` bodies pass `bash -n`. `# BOOTSTRAP - the two rows below are PRE-REFIT…`
matches the script's `grep -q '^# BOOTSTRAP'`, so the drift band really is disarmed and the `BANDED`
single-row assertion is skipped with it — the "safe failure" §13 describes is real. The `awk '$1==c
{print $2}'` lookup cannot be fooled by any comment line in the file (`$1` is `#` on all of them).
The enforced numbers and the two banded rows' prose agree at 0.065 / 0.072; the arithmetic behind the
narrowing checks out (`0.126/1.435 = 0.0878`, margin `0.0878 − 0.065 = 0.023`; `0.065/0.0192 = 3.4×`,
`0.065/0.0150 = 4.3×`). The remaining prose is M3.

**Docs figures (review area 7).** Every figure I could re-derive matches the binary or the formula:
1G/c8 → `32 bins × 6 MiB`; 4G/c8 → 31 MiB; floors 296 / 740 (no FastQC) and 349 / 1162, plus 1232 /
1655 at `--cores 64`; `clump-only` 290 MiB at c1–c4 and 544 at c32; the `--cores 16 --fastqc` error
message's `215908 bytes` and `try ≥ 821 MiB` (it is the `--memory 700M` cell, and 820 M genuinely
fails while 821 M passes); the `--memory 100M --cores 6` fallback's `need ≥ 339 MiB`; and
`1G → 1363M` restoring the old `B` at both c2 and c4.

Two things worth crediting. **The docs corrected the plan**: §6d says `2G/c8 → 2833M`, and 2833 M is
short by 2662 bytes — 2834 M is the first value that restores the old bin. That is Step 4c.3's
"re-derive against the binary, not the table" instruction doing its job a second time. And the caution
box's "~30% at **every** budget" is *exactly* right rather than approximately: `B_new/B_old` is
`denom_old/denom_new`, which cancels `usable` entirely, and it is 0.696 for every `(n_bins, workers)`
shape (496/712, 1248/1792, 432/620, 4992/7168 all agree to three places).

**No surviving stale constants.** Nothing in `src/`, `tests/`, `docs/src/`, `.github/`, `README.md` or
`CHANGELOG.md` still names `23`/`64` as the coefficients, `496`/`624`/`1248` as denominators,
`19_006_356`, `281 MiB`, `277 MiB`, `need ≥ 311`, `309806`, `σ = 23` or `k = 64`. `--help`'s `--memory`
text carries no figures, so it needed no change.

---

## Observations, not findings

- **Arming the promise gate on Linux before Linux has been measured post-refit is by design.** §5 row 9
  makes the Linux run the merge gate, and this PR's own CI run is that measurement. Pre-refit Linux read
  0.914 (`25bp_1G_c4`) and 1.030 (`25bp_1G_c2`); both predicted well under 1.0 after the refit, so a red
  here would be information rather than noise. Worth stating only because "gate armed on unmeasured
  platform" looks like a defect until you read §5.
- `channel_depths_are_pinned_to_the_calibrated_values` correctly stops asserting the false
  `work + 1 + result_per_core == 4` arithmetic and pins the depths themselves, which is what §6e's
  "lever A's second tripwire" needs. The docstring on `ChannelDepths` was updated to match.
- `resolve_layout`'s rewritten provenance paragraph is seven lines. CLAUDE.md's one-line rule is aimed
  at inline comments; a doc comment for the module's core sizing decision is the right place for it, and
  every clause is now true.
