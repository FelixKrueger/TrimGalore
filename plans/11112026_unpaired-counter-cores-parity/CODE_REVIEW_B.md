# CODE REVIEW B — #544 `r1_unpaired` / `r2_unpaired` cores-parity fix

**Reviewer:** B (independent; no contact with Reviewer A or the plan-manager)
**Target:** uncommitted working-tree change at `/Users/fkrueger/Github/TrimGalore`, branch `dev`, base `86aa301`
**Scope:** `src/parallel.rs` (+203/−32), `CHANGELOG.md` (+7), `docs/src/content/docs/reference/changelog.md` (+7)
**Date:** 2026-11-11

---

## Summary

**The fix is correct and I recommend it.** The two increments now sit inside exactly the predicate
the serial path uses (`Some(writer) && r*_ok`), immediately after the write, which is also where Perl
0.6.11 puts `++$read_1_printed`. I re-derived the Perl grounding, re-traced the call graph, and
confirmed the behaviour end-to-end against the shipped binary at four core/rescue combinations plus
the gzip, multi-batch and `--clumpify` routes the plan did not exercise. No `PairValidationStats` or
`TrimStats` field other than the two counters moves, and no output record changes.

I found **no Critical and no High issues.** Two Medium findings are about the *test*, not the fix: the
new fixture is 20 pairs and plain-output only, so it covers neither the gzip rescue-writer branch nor
cross-batch accumulation of these two fields — both of which I verified correct by hand but which the
suite still cannot catch. The rest are Low: one assertion in the new test is unreachable-as-a-failure,
and the `gen_multiome_rows` doc comment names the wrong flag and a length the filter never sees.

The plan's implementation notes hold up on every claim I checked, including the four branch
assignments, the adapter-inertness argument, and the nil blast radius of the `gen_multiome_rows` edit.

---

## What I verified by running

All runs used `./target/release/trim_galore` (current with the fix) on a fixture built outside the
repo reproducing the new test's shapes exactly — `(3,10,40) (5,40,10) (2,26,26) (4,10,10) (6,40,40)`
over `seq(n) = "ACGT"` cycled — with `-q 0 -a AGATCGGAAGAGC --stringency 1 -e 0.1 --length 30
--length_1 20 --length_2 20`, i.e. the unit test's `parity_test_config()` plus `length_cutoff: 30`
and `UnpairedLengths { r1: 20, r2: 20 }`. `src/main.rs:1736` routes on `cli.cores > 1` with no size
threshold, so `--cores 4` genuinely takes the parallel path at 20 records.

| Check | Result |
|---|---|
| **Branch assignment (brief §1)** | **Confirmed directly**, by read ID. `unpaired_1` holds reads 3–7 (shape B) **and 8–9 (shape C)**; `unpaired_2` holds reads 0–2 (shape A) **and 8–9 (shape C)**; reads 10–13 (D) appear in neither; 14–19 (P) in `val_1`/`val_2`. So C contributes to **both** sides — the #245 shape is exercised, not merely implied by the arithmetic. |
| **Anchors 7 / 5** | `Unpaired R1 kept: 7`, `Unpaired R2 kept: 5` at `--cores 1` **and** `4`, files holding 7 and 5 records. Matches the hand-derived constants. |
| **No-rescue case** | No `Unpaired` line in stderr and **no `Unpaired read N kept` block in the trimming report** at `--cores 1` or `4`, with plain and gzip output, and no unpaired files on disk. |
| **Adapter inertness (brief §2)** | **Confirmed empirically**, not just by arithmetic: every output record retains its input length (40, 26). Separately, a control record ending in `A` (`"TGCA"×12`, 48 bp) **was** trimmed to 47 bp under the same config — so the hazard the fixture avoids is real, and the `n % 4 == 1` reasoning at `src/parallel.rs:1730` is sound. Last bases are `C`, `T`, `C` for n = 10, 40, 26. |
| **JSON** | Keys present throughout. `r1_unpaired`/`r2_unpaired` = 0/0 no-rescue at both core counts; 7/5 with `--retain_unpaired` at both. |
| **Report text block** | Lives only in the **R2** trimming report (`input_R2.fq_trimming_report.txt:56-57`), not R1. |
| **gzip rescue branch** (untested in the suite — see M-1) | gzipped inputs → `--retain_unpaired --cores 1` and `4`: 7/5 reported, `.fq.gz` files decode to 7 and 5 records. |
| **Multi-batch / multi-worker** (untested in the suite — see M-1) | Same shapes ×300 = 6,000 pairs (≥ 2 batches at `BATCH_SIZE = 4096`): 2,100 / 1,500 at `--cores 1` and `4`, files matching. `PairValidationStats::merge` sums these fields correctly in practice. |
| **`--clumpify` route** | `--clumpify --retain_unpaired --cores 4` on the 6,000-pair fixture: 2,100 / 1,500, files matching. Without rescue: zero `Unpaired` lines. |
| **CHANGELOG mirror** | The six added lines are byte-identical between `CHANGELOG.md` and the docs copy. |
| **Working tree** | `git diff --numstat` is still exactly `CHANGELOG.md 7/0`, `docs/.../changelog.md 7/0`, `src/parallel.rs 203/32`, on `dev` @ `86aa301`. I ran no mutations and edited no source file. |

## What I verified by reading

- **Perl 0.6.11** (`git show 0.6.11:trim_galore`): `++$read_1_printed` / `++$read_2_printed` sit
  inside `if ($retain){ … }`, after the four `print UNPAIRED1` lines, under the `length($seq_1) >=
  $length_read_1` test (around :2325-2341); the report at :2396 calls the figure
  `Number of unpaired read 1 reads printed`. The post-fix Rust predicate and placement match.
- **Single increment site per side, crate-wide.** `grep` for `r1_unpaired|r2_unpaired` across `src/`
  returns exactly two `+= 1` sites in `src/parallel.rs` (`:643`, `:649`) and two in `src/trimmer.rs`
  (`:588`, `:594`). No fourth.
- **Call graph.** `classify_paired` → one caller `process_pairs`; `PairOutcome::Unpaired` → one
  consumer; `process_pairs` → two call sites (`:386` gzip, `:415` plain), both inside
  `process_paired_batch`, which has one caller (`:142`). `retain_unpaired` at `:108` is
  `unpaired_r1_path.is_some()`, from `src/main.rs:1706-1713` where both paths are set or neither.
  `--clumpify` changes only the *reader* (`read_pairs_clumpy` at `:194`); workers still go through
  `process_paired_batch` → `process_pairs`, so it cannot bypass the fix. Confirmed by the clumpify
  run above.
- **uBAM route untouched.** `run_paired_end_to_bam`'s `TooShort` arm (`src/trimmer.rs:805-810`)
  neither rescues nor counts.
- **No collateral counters.** The removed block in `classify_paired` (`:567-571`) took only the two
  unpaired increments with it; `pairs_removed`, `stats_r*.too_short` and the returned
  `PairOutcome::Unpaired { r1_ok, r2_ok }` are unchanged. In `process_pairs` the
  `passthrough_records_dropped` increment stays outside both writer guards, so passthrough drop
  accounting is unaffected.
- **Error ordering.** `r1.write_to(*w)?` precedes the increment, so a failed write cannot inflate the
  counter — same ordering as `src/trimmer.rs:587-588`.
- **Step 8 blast radius is nil.** All five `gen_multiome_rows` callers (`:2567`, `:2648`, `:2760`,
  `:2804`, `:2851`) use `length_cutoff = 20` and `UnpairedLengths { r1: 35, r2: 35 }`, and none
  asserts a record length, `pairs_removed`, a bp total, or R2 sequence content — only record counts,
  ID lockstep, passthrough `kept`/`dropped`/`checked`, and serial-vs-parallel equality. Short-R1 rows
  are dropped as pairs either way, so **R1/R2/passthrough output bytes are unchanged** by the 16→48 bp
  R2 edit; only `r2_ok`, `stats_r2.total_bp_after_trim` and `total_reads`-adjacent sums move, none of
  which is anchored absolutely. The R2 adapter-event count is also unchanged: both the old 16 bp and
  the new 48 bp R2 end in `A` and take the same 1 bp trim.
- **`--passthrough` + `--retain_unpaired` is rejected** at `src/cli.rs:948-953`, so the plan's
  Assumption 5 holds and the passthrough tests legitimately pin the no-rescue case.
- **Unpaired output bytes are core-count invariant.** On the 6,000-pair fixture, the decompressed
  `*_unpaired_{1,2}.fq.gz` are SHA-256-identical between `--cores 1` and `--cores 4`
  (`e4144fc7…` / `7390c018…`, both non-empty). They differ under `--clumpify`, which is expected —
  the clumpy reader reorders bins — and is not a regression from this change.

## What I did **not** re-run

`cargo clippy --all-targets --release -- -D warnings` and `cargo fmt --all -- --check`. I take the
implementer's result on those; CI gates them anyway, and the local toolchain is behind CI's.

---

## Issues by area

### Logic

Nothing wrong found. Three things I specifically tried to break and could not:

- **Asymmetric writers.** `process_pairs` would count R1 without R2 if handed
  `Some(writer_up_r1), None`. Unreachable: `process_paired_batch` gates both on the single
  `retain_unpaired` bool (`:370-379`, `:430-438`), itself `unpaired_r1_path.is_some()`, and
  `src/main.rs:1706-1713` sets both paths or neither. The only caller passing an asymmetric pair
  would be a future one.
- **Double counting.** There is no remaining increment upstream, and `PairOutcome::Unpaired` has one
  consumer.
- **Rescue regression (#245).** Shape C (both mates 26 bp under `--length 30`, both ≥ `--length_{1,2}`
  20) is rescued on **both** sides — reads 8 and 9 appear in both unpaired files. A fix that
  reintroduced the beta.5 "other mate must have failed" gate would drop them and fail the 7/5 anchors.

### Efficiency

No impact. Two `usize` increments move inside guards that were already evaluated for the writes; the
`TooShort` arm loses two branches. The new test adds a 20-record fixture driven four times (two paths
× two rescue modes) — negligible.

### Errors

- Write-then-count ordering is correct (`?` precedes `+= 1`), matching the serial path. In the
  parallel path the sink is an in-memory `Vec<u8>` or `GzEncoder`, so a failure there aborts the whole
  batch regardless.
- `fs::read_to_string(p)?.lines().count() / 4` in the new test floors a truncated trailing record
  rather than failing. Acceptable in a test whose writer it also controls.

### Structure

- The destination guard is literally the serial predicate, so the two paths are now
  line-for-line comparable — that is the right shape for a defect class that has now recurred twice
  (#541, #544).
- `classify_paired` is now purely classification plus pair-level accounting; all writer-dependent
  accounting lives with the writes. Good separation, and it is the separation whose absence caused
  the bug.

---

## Recommendations

### Critical

None.

### High

None.

### Medium

**M-1. The new test covers one batch and plain output only, so two routes the fix touches remain
untested.** `src/parallel.rs:1840-1864` passes `cores = 4, gzip = false` over 20 pairs.

- `BATCH_SIZE = 4096` (`src/parallel.rs:33`), so 20 pairs is **one batch dispatched to one worker**.
  Cross-batch accumulation of `r1_unpaired`/`r2_unpaired` through `PairValidationStats::merge`
  (`src/report.rs:148-149`) is never exercised with the rescue writer open — and `merge` has no direct
  unit test either (the only `merge` test, `src/report.rs:1937-1938`, is for `TrimStats`).
- `gzip = false` means the `GzEncoder` rescue sinks at `src/parallel.rs:370-379` — `gz_up_r1` /
  `gz_up_r2` creation, `finish()`, and the main-thread concatenation of those members — are
  **never driven with content by any test**. I checked every `run_paired_end_parallel` call in the
  module: the new one at `:1840` is the *only* test passing `Some` unpaired paths, and it is the
  plain branch.

  *Failure scenario:* a future change that accumulated into a per-batch field instead of merging, or
  that wrote the rescue members in the wrong order / dropped the final one, would pass the suite. Any
  real `--paired --retain_unpaired --cores 4` run on gzipped input over 4,096 pairs — i.e. the normal
  case — exercises both.

  *This is coverage, not a defect:* I verified both routes correct in the shipped binary (6,000 pairs
  → 2,100/1,500 at `--cores 1` and `4` with matching files; gzipped input → 7/5 with `.fq.gz` files
  decoding to 7 and 5).

  *Suggested change, cheap:* wrap the existing body in `for gzip in [false, true]` and replicate
  `shapes` by a factor `R` chosen so `20 * R > 4096` (e.g. `R = 300` → 6,000 pairs), scaling the two
  constants by `R`. I ran exactly that shape through the binary and it holds at 2,100/1,500.
  **Confidence: high** that the gap is real; high that the code is nonetheless correct.

**M-2. The rewritten `gen_multiome_rows` doc comment names the wrong flag and a length the filter
never sees** (`src/parallel.rs:2452-2458`).

- "*a 48 bp R2 that clears the `length_1 = 35` rescue cutoff*" — R2 is gated by `--length_2`
  (`UnpairedLengths.r2`), not `--length_1`. All five callers happen to pass `{ r1: 35, r2: 35 }`, so
  the number is right and the conclusion holds, but the flag name is wrong in a fixture comment whose
  whole job is to tell the next editor which thresholds the shape depends on.
- "*48 bp throughout*" (`:2473`) is the **input** length. The sequence
  `TGCATGCA…TGCA` ends in `A`, and under this module's `min_overlap: 1` + Illumina adapter the
  trailing base is trimmed — I measured it: a 48 bp `"TGCA"×12` record comes out at **47 bp**. The
  filter therefore sees 47, not 48. Still ≥ 35, so `r2_ok` is unaffected and nothing downstream
  changes.

  What makes this worth fixing rather than ignoring: the *sibling* comment added in the same diff
  (`:1730-1731`) goes out of its way to explain this exact hazard — "*lengths avoid `n % 4 == 1` so no
  sequence ends in `A`, which at `min_overlap 1` would match the adapter's first base*" — and then the
  fixture 700 lines later states a pre-trim length as if it were the filtered one. A reader who trusts
  both will conclude the hazard does not apply here.

  *Suggested wording:* "… against a 48 bp R2 (47 bp after the `min_overlap 1` trim of its trailing
  `A`) that clears the `--length_2 = 35` rescue cutoff." **Confidence: high** on both points
  (flag name read from all five callers; the 1 bp trim measured).

### Low

**L-1. The fixture-discrimination guard cannot fail** (`src/parallel.rs:1899-1905`). It is reached
only after `assert_eq!(got.r1_unpaired, 7)` and `assert_eq!(got.r2_unpaired, 5)` have passed for both
paths (`:1877-1885`) and after the four file-count assertions, so `r1_unpaired > 0 && r2_unpaired > 0
&& r1_unpaired != r2_unpaired` is unconditionally true at that point. Its stated purpose — catching a
future fixture edit that silently drops shape A or B — is already served by the absolute anchors,
which such an edit would break first.

This is the same "assertion that cannot fail" pattern the plan itself indicts in §2 for the two blind
parity tests, so it is worth not leaving one behind. Either delete it, or move it to where it bites:
asserting it against the **shape table** rather than the result, e.g. that the shapes yielding `r1_ok`
and those yielding `r2_ok` have non-zero and unequal totals, which is a statement about the fixture
and would fail on a fixture edit before any run happens. **Confidence: certain** (control flow read
directly).

**L-2. The CHANGELOG entry is one sentence, but a 56-word one spliced with a semicolon.**
`CHANGELOG.md:8-14`. It matches this file's house style — the neighbouring #541 entry is comparably
long, and the 2.4.0 entries run to several sentences, so this is not the nf-core one-line format and
should not be forced into it. But the clause after the semicolon is independent and carries the second
effect, so the entry is two sentences wearing one punctuation mark. Trimming
"*whether or not `--retain_unpaired` had opened a rescue file*" to "*without `--retain_unpaired`*"
would cut ~10 words with no loss. Content is accurate on every point I checked, including that the
JSON keys stay present while their values change. Mirror into the docs copy is byte-identical.

One wording note: the headline names `Unpaired read N kept`, which is the **trimming report** line
(`src/report.rs:469-479`, and I confirmed it appears only in the R2 report). The `--cores`-dependent
stderr summary line is `Unpaired R1 kept` / `Unpaired R2 kept` (`src/main.rs:1902-1910`) and also
disappears. The report line is the archived, MultiQC-visible one, so leading with it is the right
choice; no change needed.

**L-3. Source comment** (`src/parallel.rs:637-638`) is one wrapped sentence stating the current fact
and does not narrate the change — compliant. The issue-number citation matches house style
throughout this module. If tightening to a single line is wanted:
`// Counted at the write, so the figure matches the records actually written (Perl 0.6.11).`
No action required.

**L-4. Answer to brief §6: no docs change is needed.** I grepped `docs/` for the counters. Nothing
outside the mirrored changelog documents `r1_unpaired`/`r2_unpaired` or the `Unpaired read N kept`
block: `docs/src/content/docs/guide/reports.md` reproduces the paired report *without* that block
(`:103-108`), so removing the lines creates no doc drift.

Adjacent but pre-existing and **out of scope for this change**: `docs/src/content/docs/guide/outputs.md:35-36`
("Read 1 reads whose mate dropped below `--length_2`") and `docs/src/content/docs/guide/paired-end.md:65-66`
("R1 reads whose mate was too short") describe rescue as conditioned on the *other* mate failing. That
is true under the defaults (`--length 20` with `--length_1/2 35`, where a read failing `--length`
cannot clear 35) but false in general — shape C in the new fixture has both mates short and both
rescued. It is the same misconception as the beta.5 regression #245. Worth a separate issue if
anything, not a rider on this fix.

**L-5. Perl parity is restored for the counter, not for the block's gating.** Perl gates the report
block on `if ($retain)` (`0.6.11:trim_galore:2395-2399`) and prints two zero lines when nothing was
rescued; Rust gates on `> 0` (`src/report.rs:469`, `src/main.rs:1902`, `:2123`), so a
`--retain_unpaired` run that rescued nothing prints no block at all. That divergence predates this
change and the serial path has always had it, so the CHANGELOG's "*as the serial path and Perl v0.6.x
already did*" is accurate about **the counting**, which is what it claims. Flagging only so it is not
mistaken for a new gap if someone diffs against Perl output next. **Confidence: high.**

**L-6. The new test compares only `PairValidationStats`,** not `serial.0 == parallel.0` /
`serial.1 == parallel.1` (`src/parallel.rs:1866-1872`). It is the first parity test in the module with
a **non-zero `length_cutoff` and no passthrough**, so it is the natural place to pin `stats_r*.too_short`
and `total_bp_after_trim` under filtering — two fields the #541 fix touched. Today `too_short` parity
is covered only incidentally, by `test_passthrough_serial_parallel_parity` at `length_cutoff = 20`
(`:2710-2711`), and not at all by `test_parallel_serial_bp_after_trim_parity_paired`, whose
`length_cutoff: 0` makes the arm unreachable. Two added lines would close that.

**L-7. `fresh_tmpdir` uses a fixed slug and `remove_dir_all`s it first** (`src/parallel.rs:1406-1411`),
so two concurrent `cargo test` runs of this crate sharing a `TMPDIR` delete each other's fixtures
mid-run. Pre-existing module-wide pattern that the new test inherits, not introduced here; mentioned
because two reviewers running the suite at once on this machine is exactly the collision, and a
failure from it looks like a real regression.

**L-8. The rescue files are checked by record count, not by which records.** `src/parallel.rs:1889-1895`
counts lines / 4. A bug that wrote the correct number of *wrong* records (e.g. the mate's sequence)
would pass. The writes are untouched by this change, so this is not a gap the fix opens; I verified
externally by read ID that reads 3–9 land in `unpaired_1` and 0–2, 8–9 in `unpaired_2`. Asserting the
ID set rather than the count would make the test say what it means, at the cost of a few lines.

---

## What the tests still cannot catch — named states

1. **A paired `--retain_unpaired` run at `--cores >= 2` over more than 4,096 pairs.** One batch in the
   new test ⇒ `PairValidationStats::merge` of these two fields is never exercised with the writer open.
   (Verified correct by hand: 6,000 pairs → 2,100/1,500.) — M-1
2. **A paired `--retain_unpaired` run at `--cores >= 2` with gzip output.** No test anywhere passes
   `Some` unpaired paths with `gzip = true`, so `gz_up_r1`/`gz_up_r2` are never driven with content.
   (Verified correct by hand.) — M-1
3. **`--clumpify --retain_unpaired`.** No test passes `Some` unpaired paths together with a
   `clump_layout`. Counts are layout-independent and I confirmed 2,100/1,500 through that route, but
   the unpaired records come out in bin-flush order, so unpaired output is **not** byte-identical to
   `--cores 1` there (it is under round-robin — I checked both). If that ever needs to hold, nothing
   would tell you it broke.
4. **Which records reach the rescue files**, as opposed to how many. — L-8
5. **`--retain_unpaired` with zero rescued reads**, where the report block vanishes entirely rather
   than printing zeros. No test pins the empty-rescue shape on either path. — L-5

---

## Independent confirmation of the implementer's evidence

I re-ran the suite myself from the crate root with my own `CARGO_TARGET_DIR`:

```
cargo test  →  exit 0, 845 passed / 0 failed / 0 ignored across 21 test binaries, 0 compiler warnings
  parallel::tests::test_parallel_serial_unpaired_counter_parity ... ok
  parallel::tests::test_passthrough_serial_parallel_parity ... ok
```

| Implementer claim | My finding |
|---|---|
| 845 passed / 0 failed (844 + 1 new) | **Confirmed**, exit 0, independently run |
| No-rescue runs show no `Unpaired read` block at `--cores 1` or `4`, with and without clipping | **Confirmed** (plain and gzip output; also confirmed no unpaired files on disk) |
| Rescue runs report the same figure at both core counts, files matching | **Confirmed** at 7/5 on the unit fixture and 2,100/1,500 at 6,000 pairs; also through `--clumpify` |
| JSON values 0 throughout with keys present | **Confirmed** |
| Trimmed records byte-identical to pre-fix | Not re-derived — I did not build the pre-fix binary. The fix cannot reach any writer other than the two rescue streams, and the suite's decoded-record identity assertions pass, so I have no reason to doubt it. |
| Mutation: all four increments removed ⇒ parity passes at 0/0, only the absolute anchor fires | Not re-run (I ran no mutations, per the brief). Consistent with the control flow I read: `serial == parallel` holds trivially at 0/0, and `assert_eq!(got.r1_unpaired, 7)` is the only assertion that can then fire. |
| Mutation: reintroducing the bug fails two tests | Not re-run. Consistent: with the shared fixture now producing `r2_ok == true` on ~1/3 of rows and no rescue writer, the pre-fix code would give the parallel path `r2_unpaired ≈ 4000` against the serial path's 0, failing `test_passthrough_serial_parallel_parity`'s whole-struct comparison at `src/parallel.rs:2700-2705`. |
| Four shapes land in the four intended branches (brief §1) | **Confirmed directly by read ID**, including shape C on both sides |
| Adapter trimming inert on the fixture (brief §2) | **Confirmed by measurement**, and the stated mechanism verified by a positive control |
| `gen_multiome_rows` edit disturbs no other caller (brief §3) | **Confirmed** — all five callers read, no length/`pairs_removed`/bp/drop-rate assertion shifts, and output bytes are unchanged because those rows are dropped as pairs either way |
| Nothing else perturbed (brief §4) | **Confirmed** for `pairs_removed`, `too_short`, `passthrough_records_dropped`, and the uBAM and `--clumpify` routes |

One nuance on the plan's §12 "deviation" note: the guard idiom the new test borrowed from
`test_parallel_serial_bp_after_trim_parity_paired` (`:1714-1719`) *can* fail there, because
`total_bp_after_trim > total_bp_written` is not implied by that test's preceding equalities. In the
new test the analogous guard is implied — see L-1.

---

## Verdict

**Ship it.** The fix is minimal, correct, grounded in Perl 0.6.11 rather than in our own convention,
and closes the only route that could reach the counters without passing the write. The accompanying
test genuinely discriminates — it fails if either increment is dropped, and the absolute anchors
7/5 defeat a degenerate zero-everything fix — and Step 8 turns a fixture that had been shaped around
this divergence back into one that exposes it.

The Medium findings are both about reach, not correctness: I would take M-1 (loop `gzip`, scale the
fixture past `BATCH_SIZE`) before merge, since it is a few lines and covers the two routes a real
`--paired --retain_unpaired --cores 4` run actually takes, and M-2 (two words in a fixture comment) in
the same pass. L-1 is worth doing in the same breath — leaving an assertion that cannot fail inside
the fix for a class of bug caused by assertions that could not fail is the one thing here that would
read badly later. Everything else is optional.

**Report path:** `/Users/fkrueger/Github/TrimGalore/plans/11112026_unpaired-counter-cores-parity/CODE_REVIEW_B.md`
