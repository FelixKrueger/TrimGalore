# CODE REVIEW A — `r1_unpaired` / `r2_unpaired` parity fix (#544)

**Reviewer:** A (independent; no contact with Reviewer B or the coverage audit)
**Target:** uncommitted working tree at `/Users/fkrueger/Github/TrimGalore`, branch `dev`, base `86aa301`
**Diff:** 3 files, +217/−32 — `src/parallel.rs` (+203/−32), `CHANGELOG.md` (+7/−0),
`docs/src/content/docs/reference/changelog.md` (+7/−0)
**Date:** 2026-11-11

---

## Summary

**The fix is sound.** The two increments move from `classify_paired`'s `TooShort` arm into the two
existing `if let Some(ref mut w) = writer_up_r* && r*_ok` guards in `process_pairs`, which is
*exactly* the serial path's predicate and exactly Perl 0.6.11's placement. I verified this both by
reading the call graph and by running the shipped binary end to end on a purpose-built fixture, at
`--cores 1` and `--cores 4`, with and without `--retain_unpaired`, on **both** the plain and the gzip
output branches. Every figure now matches the records on disk.

I was asked to spend attention where the implementer's evidence was weakest. All four of the
substantive claims I was asked to challenge **hold**, and I confirmed the two most load-bearing ones
by direct observation rather than by re-deriving them on paper:

- **The four fixture shapes do land in their four intended filter branches**, shape C included. I
  proved this by read ID, not by arithmetic — see §1.
- **Adapter trimming is genuinely inert on the fixture.** The trimming report says
  `Reads with adapters: 0 (0.0%)` on both mates and no output record changed length — see §2.
- **The `gen_multiome_rows` change disturbs none of its five callers.** I read all five; none asserts
  a record length, an absolute bp figure, or a drop count that the 16 → 48 bp R2 change moves — see §3.
- **Nothing else is perturbed.** `pairs_removed`, `too_short`, the passthrough drop counter, the uBAM
  route and the clumpify route are all untouched or provably unreachable — see §4.

Findings are all **Low**. There is no Critical or High issue, and I am not manufacturing a Medium to
balance the ledger. Four are wording/coverage nits; one (§5.2) is a latent invariant the fix now
leans on that is asserted nowhere, which I would note in passing rather than fix in this change.

**Confidence that the fix is correct: high.** It is a two-statement move into a guard that already
existed and already governed the writes, and I observed the corrected output directly.

---

## What I verified, and how

I worked on a pinned copy of the working tree under the scratchpad for anything that compiles, and
ran the repo's existing `./target/release/trim_galore` (current with the fix) for the end-to-end
checks. **I made no edit to any file in `/Users/fkrueger/Github/TrimGalore` other than this report.**
`git diff --numstat` is unchanged at 3 files, +203/−32, +7/−0, +7/−0 (re-confirmed at the end).

| | Verified by running | Verified by reading | Taken on the implementer's word |
|---|---|---|---|
| Branch assignment of all five shapes | ✅ by read ID | ✅ against `src/filters.rs:88-107` | — |
| Adapter inertness on the fixture | ✅ report + lengths | — | — |
| No-rescue counters 0 at `--cores 1` / `4` | ✅ plain **and** gzip | ✅ | — |
| Rescue counters = records on disk | ✅ plain **and** gzip | ✅ | — |
| JSON keys present, values move 4 → 0 | ✅ | ✅ `src/report.rs:1073-1074` | — |
| `gen_multiome_rows` blast radius | ✅ suite green | ✅ all five callers | — |
| Call graph / no fourth increment site | — | ✅ `grep` over `src/` | — |
| `cargo test --lib` | ✅ 519 passed, 0 failed, 0 warnings | — | — |
| Full `cargo test` count (845) | ✅ see §6 | — | — |

---

## 1. Logic — do the four shapes land in the four intended branches?

**Verified by direct observation. They do, shape C included.**

The arithmetic self-check the brief mentions (7 and 5 rather than 5 and 3) is real but weak — it
cannot distinguish "C contributed to both" from "A and B each contributed one extra". So I built the
fixture as FASTQ with the shape encoded in each read ID and ran the real binary with flags matching
`parity_test_config()` plus `length_cutoff: 30` / `UnpairedLengths { r1: 20, r2: 20 }`:

```
--paired --retain_unpaired -a AGATCGGAAGAGC -a2 AGATCGGAAGAGC
-q 0 --stringency 1 -e 0.1 --length 30 -r 20 --length_2 20
```

The unpaired files then name which shape reached which side:

| Shape | (count, r1, r2) | In `_unpaired_1` | In `_unpaired_2` | Intended | Verdict |
|---|---|---|---|---|---|
| A | (3, 10, 40) | — | `A_read_0,1,2` | `r2_ok` only | ✅ |
| B | (5, 40, 10) | `B_read_3..7` | — | `r1_ok` only | ✅ |
| **C** | (2, 26, 26) | **`C_read_8,9`** | **`C_read_8,9`** | **both (#245)** | ✅ |
| D | (4, 10, 10) | — | — | neither | ✅ |
| P | (6, 40, 40) | — | — | `Pass` → `_val_{1,2}` | ✅ |

Shape C appears in **both** unpaired files by ID, so the `r1_ok == r2_ok == true` branch is genuinely
exercised. Totals: R1 = 5 + 2 = **7**, R2 = 3 + 2 = **5**, matching `EXPECTED_R1_UNPAIRED` /
`EXPECTED_R2_UNPAIRED` at `src/parallel.rs:1750-1751`. Identical at `--cores 1` and `--cores 4`.

This is also the strongest thing in the change: shape C is what stops a "fix" that reintroduces the
beta.5 #245 gate (`rescue only if the other mate failed`) from passing as parity, and it is now
pinned by an absolute anchor rather than only by cross-path equality.

**Confidence: high** (observed, not inferred).

---

## 2. Logic — is adapter trimming really inert on the fixture?

**Verified, and by a stronger signal than the implementer's reasoning.**

The implementer's argument is that `seq(n)` cycles `ACGT` and the chosen lengths avoid `n % 4 == 1`,
so no sequence ends in `A` and nothing can match the adapter's first base at `min_overlap: 1`. I
confirmed the premise — `seq(10)` and `seq(26)` end in `C`, `seq(40)` ends in `T` — but that argument
only covers the 1-base overlap, not longer suffixes at `error_rate: 0.1`. Rather than hand-compute
13 overlap lengths against Cutadapt's semi-global DP, I read it off the run:

```
=== Summary (Read 1) ===   Reads with adapters:   0 (0.0%)
=== Summary (Read 2) ===   Reads with adapters:   0 (0.0%)
```

and every output record retained its input length (40 → 40, 26 → 26, in `_val_*` and `_unpaired_*`
alike). So no suffix of any fixture sequence matches the Illumina adapter at `min_overlap: 1`,
`error_rate: 0.1` — not at length 1 and not at any longer overlap. Branch assignment cannot be
perturbed by trimming.

Two supporting points worth recording, because they are what make the lengths safe rather than lucky:

- `quality_cutoff: 0` means quality trimming is **skipped entirely**, not merely a no-op —
  `src/trimmer.rs:94` guards the whole stage on `> 0`.
- `parity_test_config()` leaves `clip_r*`, `three_prime_clip_r*`, `trim_n`, `rrbs`, `poly_a`,
  `poly_g`, `max_n` and `max_length` all off, so adapter trimming is the *only* length-modifying
  stage in play. That is why a single inertness check settles the whole fixture.

The 10 bp records never reach any output (D is dropped; A's R1 and B's R2 are not rescuable), but
their exact length is immaterial anyway — a trim could only shorten them, and they are already below
both the 30 bp discard cutoff and the 20 bp rescue cutoff.

**Confidence: high.**

---

## 3. Structure — does the `gen_multiome_rows` change disturb its other callers?

**Verified by reading all five, plus a green suite. It does not.**

The change makes R2 48 bp on every row (`src/parallel.rs:2473-2475`), where short-R1 rows previously
carried a 16 bp R2. With `passthrough_test_config()` (`length_cutoff: 20`) and
`UnpairedLengths { r1: 35, r2: 35 }`:

- Pair verdict is unchanged. `r1_short = 16 < 20` already forces `TooShort` on those rows regardless
  of R2, so `pairs_removed`, `stats_r1.too_short` and `stats_r2.too_short` are all unmoved.
- Drop rate is unchanged — still every 3rd row, governed entirely by R1.
- Only `r2_ok` moves: `16 >= 35` → false becomes `48 >= 35` → true. That is the whole point.

The five callers and what each asserts:

| Caller | Line | Asserts | Affected? |
|---|---|---|---|
| `test_passthrough_serial_smoke` | `:2567` | pt checked/kept/dropped counts, id lockstep, pt sequence verbatim | No |
| `test_passthrough_serial_parallel_parity` | `:2648` | serial-vs-parallel `TrimStats` **and** `PairValidationStats`, decoded-record identity, lockstep, pt sanity | No — see below |
| `test_passthrough_truncation_detected_mid_stream` | `:2760` | error message text | No |
| `test_passthrough_id_mismatch_detected` | `:2804` | error message text | No |
| `test_passthrough_dont_gzip_plain_output` | `:2851` | lockstep, `r1_recs.len() == passthrough_records_kept` | No |

`:2648` is the one that deserved a second look, because it *does* assert `TrimStats` for both mates
(`serial_stats.0 == parallel_stats.0`, `.1 == .1`) and `stats_r2.total_bp_after_trim` genuinely
shifts — short rows now contribute 48 bp where they contributed 16. But every assertion there is
**serial against parallel**, not against a golden figure, so both sides move together and the
comparison is undisturbed. There is no absolute bp anchor anywhere in the passthrough tests. Same
for `assert_eq!(sr2, pr2)`: decoded-record identity between paths, not against a fixture.

I also confirmed the un-blinding actually works in the direction claimed: with the bug reintroduced,
the parallel path would increment `r2_unpaired` on every short-R1 row (4,000 of 12,000) while the
serial path stays at 0, so `serial_stats.2 != parallel_stats.2` and `:2648` fails. That matches the
implementer's mutation #4 — see §6.

**Confidence: high.**

---

## 4. Logic — could the fix perturb anything else?

**No.** Each route checked:

- **`pairs_removed` / `too_short`.** Still set in `classify_paired`'s `TooShort` arm
  (`src/parallel.rs:567-572`); only the two increments were removed. Confirmed by reading the arm.
- **`passthrough_records_dropped`.** Sits *outside* both rescue guards in the `Unpaired` arm
  (`src/parallel.rs:651-653`) and is unchanged. The old comment's reasoning about passthrough
  interaction was about the *writes*, which still behave identically.
- **uBAM output.** `trimmer::run_paired_end_to_bam`'s `TooShort { .. }` arm (`src/trimmer.rs:805-810`)
  neither rescues nor counts — `--retain_unpaired` is rejected with `--output-format ubam`. Counters
  stay 0 there, which is now consistent with the FASTQ paths rather than divergent from them.
- **`--clumpify`.** `read_pairs_clumpy` (`src/parallel.rs:804`) is a different *reader*; it dispatches
  the same `PairedWork` to the same worker channels, so `process_paired_batch` → `process_pairs` is
  shared. Clumpify cannot bypass the fix. (Record *order* in the unpaired files may differ under
  clumpify by design; counts do not.)
- **No fourth increment site.** `grep -rn "r[12]_unpaired\s*+=" src/` returns exactly four — two
  serial (`src/trimmer.rs:588`, `:594`) and two parallel (`src/parallel.rs:643`, `:649`) — plus the
  field-wise `merge` at `src/report.rs:148-149`, which aggregates per-worker counts correctly.
- **`classify_paired` has exactly one caller** (`src/parallel.rs:622`), and `process_pairs` exactly
  two call sites (`:386` gzip, `:415` plain), both inside `process_paired_batch`, which itself has one
  caller (`:142`, the worker loop). The call graph in the plan is accurate.

At both `run_paired_end_parallel` call sites (`src/main.rs:1746`, `:2033`) the unpaired path pair is
built as `(Some, Some)` or `(None, None)` from a single `cli.retain_unpaired` / `outputs.unpaired`
branch, so `writer_up_r*.is_some()` really is "`--retain_unpaired` is on". See §5.2 for the one
caveat.

**Confidence: high.**

---

## 5. Issues

All findings are **Low**. None blocks the change.

### 5.1 (Low) The rewritten `gen_multiome_rows` doc comment names the wrong flag for R2

`src/parallel.rs:2454`

```rust
/// drop rate — against a 48 bp R2 that clears the `length_1 = 35` rescue
/// cutoff.
```

R2's rescue threshold is `--length_2` / `UnpairedLengths::r2`, not `length_1`
(`src/filters.rs:101-104` routes each mate against its own side). The passthrough tests pass
`UnpairedLengths { r1: 35, r2: 35 }`, so both happen to be 35 and the **behaviour described is
correct** — but the flag named is the wrong one, in a comment whose entire job is to explain which
threshold makes `r2_ok` true. Note the sibling comment four lines down (`:2467`, about R1's 16 bp)
names `length_1` correctly, so the two now disagree about what `length_1` governs.

Suggested: `… that clears the `length_2 = 35` rescue cutoff`.

**Failure scenario:** a future reader tunes `length_1` expecting it to move `r2_ok`, finds the test
unchanged, and concludes the fixture is inert when it is not. **Confidence: high** (read against
`src/filters.rs`).

### 5.2 (Low) `retain_unpaired` is derived from the R1 path alone, and nothing asserts the pair

`src/parallel.rs:108`

```rust
let retain_unpaired = unpaired_r1_path.is_some();
```

That single bool gates **both** rescue writers (`:370-379` gzip, `:430-438` plain), while the main
thread opens the two output files independently (`:227-228`) and writes each only if its own file is
open (`:257-262`). So a caller passing `unpaired_r1_path: Some(p), unpaired_r2_path: None` would get
`gz_up_r2 == Some(...)`, R2 rescues written into `buf_up_r2`, those bytes silently discarded for want
of an output file — and, **after this change, `r2_unpaired` incremented for them**. That is the #544
defect class in a different corner.

It is **not reachable today**: both call sites build the pair as `(Some, Some)` or `(None, None)` in
one branch (`src/main.rs:1706-1713` from `cli.retain_unpaired`; `src/main.rs:2021-2029` from
`outputs.unpaired`). I verified both. So this is a latent invariant, not a live bug, and it was
equally latent before the fix — the pre-fix code mis-counted that corner too, just from the other
direction.

What makes it worth a line: the function **already carries a `debug_assert!` for exactly this shape**
on the passthrough pair twelve lines below (`:109-113`, "must be both-Some or both-None"), and the
unpaired pair — which now governs a reported figure, not just a write — has none. I would not expand
this change's scope to fix it; a one-line `debug_assert!` mirroring the passthrough one, or deriving
`retain_unpaired` from `unpaired_r1_path.is_some() && unpaired_r2_path.is_some()`, would close it
whenever this file is next touched.

**Confidence: high** that the asymmetry exists; **high** that it is unreachable from the CLI today.

### 5.3 (Low) The new test compares only `PairValidationStats`, not `TrimStats`

`src/parallel.rs:1863-1868`

```rust
assert_eq!(serial.2, parallel.2, "PairValidationStats must match …");
```

`too_short` — which the brief asked me to check for perturbation — lives in `TrimStats`
(`stats_r1.too_short`, set at `src/parallel.rs:569-570`), not in `PairValidationStats`. The new test
therefore asserts nothing about it, although the fixture drives 14 of 20 pairs through the
`TooShort` arm and both `TrimStats` structs are already in hand as `serial.0` / `serial.1` and
`parallel.0` / `parallel.1`.

The fix does not touch `too_short`, and `test_passthrough_serial_parallel_parity` does compare both
`TrimStats` on a `length_cutoff: 20` fixture, so the invariant is covered elsewhere. This is a cheap
win left on the table rather than a hole:

```rust
assert_eq!(serial.0, parallel.0, "R1 TrimStats ({tag})");
assert_eq!(serial.1, parallel.1, "R2 TrimStats ({tag})");
```

Two lines, no new fixture, and it would make this test self-sufficient on the filter accounting it
already exercises. **Confidence: high.**

### 5.4 (Low) The new test exercises only the plain-output branch

`src/parallel.rs:1859` passes `gzip: false`, so `process_paired_batch`'s `GzEncoder` branch
(`:366-412`) — the real-world default, since output compression mirrors the input — is not covered
for these counters by any test.

The increments live in `process_pairs`, which both branches share, so the risk is low by
construction. I closed it by hand rather than leaving it asserted-by-reasoning: with gzipped input at
`--cores 1` and `--cores 4`, `--retain_unpaired` reports 7 / 5 and the two `.fq.gz` files hold
exactly 7 and 5 records; without it, no unpaired file is produced, no report block appears, and both
JSON values are 0. (Worth recording for anyone repeating this: passing no `--dont_gzip` is *not*
enough to get gzip output — Trim Galore mirrors the input's compression, so the input has to be
`.gz`. My first attempt silently measured the plain branch twice.)

**Confidence: high** (verified by running).

### 5.5 (Low) The new comment's first clause restates the guard it sits above

`src/parallel.rs:637-638`

```rust
// Counted only when the rescue writer is open, so the figure
// matches the records written (#544, and Perl 0.6.11).
```

Against the project rule "don't narrate what the next line does": `if let Some(ref mut w) =
writer_up_r1` is visible on the very next line, so "counted only when the rescue writer is open" is
the half a reader can already see. The load-bearing half — the invariant that the figure must equal
the records written, and that this is a Perl-parity surface — is the second clause. The comment does
satisfy the harder rules (it states the current fact, carries no historical narrative, and the
`(#NNN)` citation is house style, cf. `src/parallel.rs:541`), and two lines is within the limit.

A tighter one-liner, if it is being touched anyway:

```rust
// The figure must equal the records written (#544, Perl 0.6.11).
```

**Confidence: medium** — this is a judgement call on a comment that is already compliant with the
letter of the rule.

### 5.6 (Low) The CHANGELOG entry carries one clause of audit trail

`CHANGELOG.md:8-13` (mirrored identically — I diffed both additions, they match byte for byte).

It is **one sentence**, so it complies with the project rule literally, and it is no longer than the
adjacent #541 entry; the 2.4.0 section above it runs to three and four sentences per entry, so this
is on the short end of the file's own house style, not an outlier. Two observations:

- The semicolon joins two independent clauses, which is two sentences wearing one coat. Splitting at
  the semicolon would read as two sentences and breach the rule more visibly than the current form
  does, so the semicolon is arguably the right call — worth stating explicitly so nobody "fixes" it
  into a breach.
- `, as the serial path and Perl v0.6.x already did` is rationale, and the rule puts rationale in the
  commit message. Against that, naming Perl is also the clearest signal to the reader deciding
  whether the figure they parse was wrong before — which is the audience the rule names. I lean
  toward keeping it.

What I would **not** cut: naming both effects (vanished text block *and* changed JSON values). The
JSON keys are emitted unconditionally (`src/report.rs:1073-1074` — I confirmed this, the plan's
revision 1 had it wrong), so a machine consumer sees a value move from 4 to 0 with no structural
hint, and that is exactly what an upgrade-decision reader needs to know.

**Confidence: medium** (style judgement).

---

## 6. Gates and mutation checks — reproduced independently

Run on a pinned copy of the working tree under the scratchpad (source files byte-identical to the
repo's, SHA-256 verified before and after), with my own `CARGO_TARGET_DIR`. Local `rustc 1.95.0`
against CI's newer toolchain, so these are necessary, not sufficient.

| Gate | Result |
|---|---|
| `cargo test` (full) | **845 passed, 0 failed, 0 ignored**, exit 0 — 21 result lines summed; matches the implementer's figure exactly |
| `cargo test --lib` | 519 passed, 0 failed, **0 warnings** |
| `cargo fmt --all -- --check` | exit 0, no output |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0, **0 warnings**, crate genuinely checked (see note) |

**Harness note, recorded because it nearly fooled me.** My first clippy log showed 45 dependency
lines, `Finished`, and no crate line, which reads exactly like a false green. It was my grep: the
package is `trim-galore` with a **hyphen** (`Cargo.toml:2`) while the binary target is `trim_galore`
with an underscore (`:25`). The crate was checked, at line 180. Anyone re-running this should grep
the hyphenated form. Separately, the crate recompiles on every clippy invocation rather than caching —
expected, since `build.rs` stamps a build timestamp — so a ~7 s "Compiling trim-galore" is a real
check, not a cache hit. Corroborating evidence that `--all-targets` does lint the new test code: the
implementer's notes record clippy initially failing on two `needless_borrows_for_generic_args` inside
it.

### Mutation checks I ran myself

Both on the pinned copy; both reverted and the copy's two source files restored to their pre-mutation
SHA-256 before I continued.

| # | Mutation | Expected | Observed |
|---|---|---|---|
| M1 | Reintroduce the original bug (increments back in `classify_paired`, out of `process_pairs`) | 2 tests fail | ✅ **exactly 2**: `test_parallel_serial_unpaired_counter_parity` **and** `test_passthrough_serial_parallel_parity`. The latter reports `r2_unpaired: 4000` parallel against `0` serial on 12,000 pairs with 4,000 dropped — the predicted 1/3, and proof the Step 8 un-blinding works |
| M2 | Remove **all four** increments, serial *and* parallel | parity passes at 0/0; only the absolute anchor fires | ✅ **exactly 1** test fails, and the firing assertion is `serial r1_unpaired (retain) must equal the records written` — the hand-derived anchor, not the cross-path comparison. `test_passthrough_serial_parallel_parity` passes, as it must at 0/0 |

M2 is the one that matters, and it independently confirms the implementer's most important claim:
**`serial == parallel` alone would accept a "fix" that simply zeroes the counters.** The absolute
anchors at `src/parallel.rs:1750-1751` are the only thing standing between this change and that
degenerate pass. They should not be softened into a cross-path comparison by any later cleanup.

### End-to-end, on the shipped binary

Purpose-built fixture (20 pairs, shape in each read ID), `./target/release/trim_galore`:

| Scenario | `--cores 1` | `--cores 4` |
|---|---|---|
| No rescue, plain out | no block, no file, JSON 0/0 | **no block, no file, JSON 0/0** |
| No rescue, gzip out | no block, no file, JSON 0/0 | **no block, no file, JSON 0/0** |
| `--retain_unpaired`, plain out | 7 / 5, files hold 7 / 5 | 7 / 5, files hold 7 / 5 |
| `--retain_unpaired`, gzip out | 7 / 5, files hold 7 / 5 | 7 / 5, files hold 7 / 5 |
| `--retain_unpaired --clumpify` | n/a (needs ≥ 2) | **7 / 5, files hold 7 / 5** |

One reporting detail worth passing on: the `Unpaired read N kept` block lands in the **R2** trimming
report, not R1. Grepping only `*_R1*_trimming_report.txt` returns nothing even on a correct rescue
run, which could easily be misread as the fix having removed the block everywhere.

---

## 7. What the tests still cannot catch

Named inputs and states, not generic gaps.

1. **`--clumpify --retain_unpaired --paired`.** `--retain_unpaired` is rejected with
   `--output-format ubam` (`src/cli.rs:652`), `--clump_only` (`:921`) and `--passthrough` (`:948`) —
   but **not** with `--clumpify`. Since `--clumpify` forces `--cores >= 2`, every such run takes the
   path this fix touches, and no test covers the combination: the clumpify test
   (`src/parallel.rs:2233`) passes `clump_layout: Some(...)` with unpaired paths `None, None`, while
   the new test passes unpaired paths `Some(...)` with `clump_layout: None`. The two never meet. I
   closed this by hand — the run above reports 7 / 5 with 7 and 5 records on disk, and the record
   *multiset* matches the non-clumpify run while the *order* differs (C before B, as the clusterer
   intends). A test would need to compare multisets, not record sequences. This is the gap I would
   most want closed, because clumpify is both the mode that always takes the parallel path and the
   one that reorders batches.
2. **`unpaired_r1_path: Some(_), unpaired_r2_path: None`** — §5.2. No test and no `debug_assert`.
   Unreachable from the CLI today; a new caller of `run_paired_end_parallel` could reach it and would
   get a count for bytes that go nowhere.
3. **uBAM input → FASTQ output with `--retain_unpaired` at `--cores >= 2`**
   (`run_paired_ubam_single_file`, `src/main.rs:2033`). The de-interleaver feeds the same parallel
   path so the counters are correct by construction, but the combination is untested.
4. **A later fixture edit that moves a shape between branches while the totals still sum to 7 / 5.**
   The guard at `src/parallel.rs:1899-1906` requires both sides non-zero and unequal, which catches a
   dropped shape A or B, but not a rebalancing that preserves both sums. Remote, and I would not add
   machinery for it.
5. **Order and content of the unpaired files under any `--cores`.** The new test counts records
   (`:1886-1895`); it does not check *which* records. A hypothetical bug that rescued the wrong mate
   while keeping the count right would pass. The end-to-end runs above cover this by ID for the
   non-clumpify case; nothing automated does.

Worth stating plainly: **CI cannot catch a regression here at all.** The validation matrix runs
`--cores 1`, i.e. the already-correct path. The new test plus the un-blinded passthrough test are the
entire guard, which is why M2's finding about the absolute anchors matters more than it looks.

---

## 8. Efficiency

No change, and nothing to raise. Two `usize` increments move inside guards that were already being
evaluated for the writes; the arithmetic per rescued record is identical and the branch count is
unchanged. `PairValidationStats::merge` (`src/report.rs:148-149`) is a field-wise sum, so per-worker
counting aggregates correctly with no added synchronisation.

One small correctness benefit, probably incidental: the increment now sits **after** `write_to(...)?`,
so a mid-batch I/O failure cannot leave the counter ahead of the file. That matches both the serial
path and Perl's `++$read_1_printed` placement.

---

## 9. Error handling

Nothing to raise. The change introduces no new fallible operation, no new `unwrap`/`expect`, and no
new error path. The `?` on `r1.write_to(*w)` is unchanged and still propagates first.

---

## 10. Recommendations

### Critical
None.

### High
None.

### Medium
None.

### Low
1. **`src/parallel.rs:2454`** — change `length_1 = 35` to `length_2 = 35` in the `gen_multiome_rows`
   doc comment. It describes R2's threshold and names R1's flag (§5.1). One word; the only finding I
   would bother fixing before merge.
2. **`src/parallel.rs:1863-1868`** — add `assert_eq!(serial.0, parallel.0)` / `(serial.1, parallel.1)` to
   the new test, bringing `too_short` and the rest of `TrimStats` under the same parity assertion the
   fixture already exercises (§5.3). Two lines, no new fixture.
3. **`src/parallel.rs:108`** — consider a `debug_assert!` that the two unpaired paths are both-`Some`
   or both-`None`, mirroring the passthrough one immediately below it (§5.2). Latent only; fine to
   defer to the next change that touches this file.
4. **Coverage** — a `--clumpify` + `--retain_unpaired` case, asserting the unpaired record *multiset*
   rather than the sequence (§7.1). The only genuinely untested legal combination on the path this
   fix changes. I verified it by hand and it is correct today.
5. **`src/parallel.rs:637-638`** — optionally tighten the comment to one line that states only the
   invariant, dropping the clause that restates the visible guard (§5.5).
6. **`CHANGELOG.md:8-13`** — compliant as written and consistent with the file; if anything is
   trimmed, it is the trailing `, as the serial path and Perl v0.6.x already did` (§5.6). I would
   leave it.

### Explicitly not recommended
- **Do not** replace the absolute anchors (`src/parallel.rs:1750-1751`, `:1876-1885`) with a
  cross-path comparison. M2 shows they are the only check that fails against a counters-always-zero
  implementation.
- **No docs change is needed.** I grepped `docs/`, `wiki/` and `README.md` for `Unpaired read`,
  `r1_unpaired` and `r2_unpaired`: the only hits are the changelog entry this change adds. The guide
  pages (`docs/src/content/docs/guide/{outputs,paired-end,length,flags}.md`) document the unpaired
  *files* and the `--length_1` / `--length_2` cutoffs, never the report lines or the JSON keys. The
  plan's §10 open question is answered: no.

---

## 11. Verdict

**The change is sound and I recommend it as-is.** The fix is a two-statement move into a guard that
already existed and already governed the writes; it is the serial path's predicate exactly, and Perl
0.6.11's placement exactly. I verified the corrected behaviour end to end at both core counts, on
both output branches, with and without rescue, and under `--clumpify`; I reproduced the two
load-bearing mutation checks independently and both behaved as claimed; and all four of the
implementer's substantive claims that I was asked to challenge hold, two of them now by direct
observation rather than by argument.

The six findings are all Low, and only the first (a wrong flag name in one comment) is worth acting
on before merge. The most valuable thing in the change is not the fix — which is nearly mechanical —
but the absolute anchors and the un-blinded passthrough fixture, because CI runs `--cores 1` and
therefore cannot see this class of defect at all.

**Working tree confirmed unchanged:** `git diff --numstat` returns exactly `CHANGELOG.md` 7/0,
`docs/src/content/docs/reference/changelog.md` 7/0, `src/parallel.rs` 203/32, on `dev` @ `86aa301`.
All my mutation work was done on a scratchpad copy and reverted; no source file in the repository was
edited by me, and the only file I created is this report.
