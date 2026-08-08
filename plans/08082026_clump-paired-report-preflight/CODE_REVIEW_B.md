# Code Review B — #391 clumping reports join the `--clump_only` collision pre-flights

**Branch:** `fix/391-clump-report-preflight` @ `2704b52` (base `dev` @ `7ea2741`)
**Plan:** `plans/08082026_clump-paired-report-preflight/PLAN.md` (r2 + Implementation notes)
**Reviewer:** B (independent; Reviewer A reviews the same diff separately)
**Files reviewed:** `src/main.rs`, `src/io.rs`, `tests/integration_output_collision.rs`, `CHANGELOG.md`
**Files edited by this review:** none (recommend-only, per instructions)

## Summary

**Verdict: approve.** No correctness defect found. Every one of the five clump arms passes exactly
the report-keyed inputs its writer uses; the `no_report_file` gate is identical in polarity to all
four writers; the widened namer provably cannot change a written path; `--clock`/`--implicon` plan
byte-identical candidate lists to before. The reworded hint constant assembles correctly (verified
by extracting the literal from the built binary — no lost or doubled spaces at the five
continuations, em-dash intact).

Findings are all non-blocking: one **pre-existing docs statement that this PR proves false**
(Medium), four small test/property-test coverage gaps (Low), and four comment/style nits (Low).

### Independent verification performed

| Check | Result |
|---|---|
| `cargo test` (full) | 557 passed, 0 failed (matches the implementation note) |
| `cargo test --test integration_output_collision` | 43 passed (9 new + 34 pre-existing) |
| `cargo test --test integration_clump_only --test integration_clump_only_ubam` | 12 + 15 passed |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets -- -D warnings` | clean |
| Hint literal, extracted from `target/debug/trim_galore` | correct spacing, em-dash preserved |
| Assembled rejection message, run against live fixtures | correct (see Logic §5) |
| **Over-rejection probe** (not in the test suite) | legitimate no-`-o` same-filename paired run still exits 0 and writes both reports (see Logic §6) |

Not verified by me: the expected-fail control against the unpatched tree. The tree is shared with
other agents and I was told not to switch branches. I confirmed the *mechanism* by inspection
instead — pre-patch the paired closure returned only the two primaries, so the two identical report
paths were never hashed, and on the SE/Shape A alias shapes the report path was never compared
against the input list. The claim is sound; I am relying on the caller's run for the empirical half.

---

## Logic

### 1. Candidate inputs match the writers' keys — all five arms

Re-derived from the writers, not from the plan:

| Arm | Writer + key | Candidate call | Match |
|---|---|---|---|
| SE FASTQ | `clump_only.rs:365-366`, `clumping_report_name(input, …)` per input | `main.rs:590`, `&cli.input` | ✅ |
| Paired FASTQ | `clump_only.rs:534-536`, both `input_r1` and `input_r2` in one gate | `main.rs:558`, `&[r1, r2]` | ✅ |
| SE BAM | `clump_only.rs:836-837`, per input | `main.rs:727`, `&cli.input` | ✅ |
| Paired BAM Shape A | `clump_only.rs:1082-1084`, ONE report keyed `inputs[0]` where `inputs` = `chunk` (`main.rs:696`) | `main.rs:672`, `std::slice::from_ref(&chunk[0])` | ✅ |
| Paired BAM Shape B | same writer, `inputs[0]` where `inputs` = `&cli.input` len 1 (`main.rs:640`) | `main.rs:633`, `&cli.input` | ✅ |

`grep clumping_report_name src/` confirms no fifth writer (the only other hits are two `clump_only.rs`
test bodies and the helper itself). `--clumpify` on the trim path writes no clumping report — the
only mention outside `clump_only.rs`/`io.rs` is a doc string in `cli.rs:186`.

### 2. Gate parity — identical polarity and position

All four writers gate on `if !no_report_file` and on nothing else: no zero-record short-circuit, no
"only when stats non-empty" condition, no early `return` between the primary write and the report
write that would skip it. The helper's `if no_report_file { return Vec::new(); }` is the exact
complement. Position differs harmlessly: the writers write reports *after* the primary and only on
success, so the pre-flight is conservative (it can plan a report a failing run never writes) — the
same conservatism #388 already shipped on the trim path.

### 3. The widened namer cannot change a written path

`run_specialty_paired` consumes `pair_outputs`' return value at exactly one place — `planned.extend(…)`
at `main.rs:2545` — and `planned` feeds only `preflight_output_collisions` at `:2547`. `run_pair` is a
separate closure that re-derives its own paths inside `clump_only::*`. Verified by reading the whole
function (`main.rs:2531-2571`); there is no second consumer.

### 4. `--clock` / `--implicon` unchanged

`(a, b)` → `vec![a, b]` preserves both element order and content, and the old driver pushed `o1` then
`o2` in that same order. The planned list is therefore byte-identical, and their hint argument
(`CWD_OUTPUT_HINT`) is untouched. `clock_collision_gets_the_cwd_hint_not_the_false_advice` still
passes, including its negative assertion that the generic advice does not appear.

### 5. Hint: correct string, and a free improvement

Extracted from the binary:

> Outputs and reports are named from the input filename alone, so inputs sharing a filename can
> collide when --output_dir (or a shared input directory) sends them to one place — rename one input,
> or pass --no_report_file if only the reports collide.

Each of the five continuation lines ends with a space before the `\`, so no words are joined. The
statement is true at all three use sites (`main.rs:552`, `:838`, `:1942`); the clause it replaced
("collide **regardless of source directory**") was false at all three, since reports follow each
mate's own parent on the FASTQ arms.

Unremarked side effect worth knowing: on the clump-paired arm a *primary*-vs-primary collision
previously rendered `GENERIC_ADVICE` ("different source directories or `--output_dir`"), which is
false advice for a user who already passed `-o` and already has different source directories. Those
messages improve too. Nothing pins `GENERIC_ADVICE` (`grep` finds only the two negative assertions in
the hardtrim/clock tests and the `io.rs` unit test, all green), so this is safe.

Pre-existing wart, unchanged and correctly deferred to the plan's follow-up list: a report-vs-report
collision prints the *same path twice*, e.g.

```
Error: Output path collision (…): /…/out/reads.fq_clumping_report.txt and
/…/out/reads.fq_clumping_report.txt would be written to the same file. Outputs and reports are …
```

The reworded hint ("inputs sharing a filename") does soften it — the user can at least infer what to
rename, which the old text could not express.

### 6. Over-rejection: probed, and clean

The direction reviews usually miss. Two checks:

- **Cross-class impossible by construction.** A primary always ends `_clumped.fq`, `_clumped.fq.gz`
  or `_clumped.bam`; a report always ends `_clumping_report.txt`, and `.txt` is not stripped by
  `strip_fastq_extensions`. So a report candidate can never fold-equal a primary candidate, and the
  new entries cannot manufacture a collision between two things that are both written.
- **The legitimate same-filename shape still runs.** I ran the debug binary on
  `--clump_only --paired a/reads.fq b/reads.fq` with **no** `-o`: exit 0, primaries
  `a/reads_clumped_1.fq` + `a/reads_clumped_2.fq`, and **two** reports, one in each mate's own
  directory. That is the shape the directory asymmetry makes legal, and the fix leaves it legal.
  It is also the one behaviour this PR does not pin — see Recommendation R2.

### 7. Edge cases re-checked upstream of the change

- Odd input counts: rejected in `validate_paired_input` (`cli.rs:527`), so `chunk[1]` is always safe.
- `--clump_only --paired` N=1 FASTQ: rejected at `main.rs:539` before the driver.
- Behaviour 2c really does reach the pre-flight: `validate_paired_input` rejects only within-pair
  R1==R2 (`cli.rs:537`) and *exact* duplicate pairs (`cli.rs:549`), and the general duplicate scan is
  `!self.paired`-gated (`cli.rs:633`). A file reused as the mate of two different pairs passes — which
  is exactly what `clump_paired_rejects_shared_mate_report_without_output_dir` exercises.
- `--basename` with 2 paired inputs is allowed (`cli.rs:651` only rejects `> 2`), so the `--basename`
  variant of the repro is a valid invocation.
- `--no_report_file` is unconditionally accepted (`cli.rs:232`, no validate interaction), so the
  hint's remedy is valid at all three hint sites including the uBAM-out one.

---

## Efficiency

Matches the plan's O(paths) claim; nothing to flag. The paired closure's `vec![o1, o2]` + `extend`
may realloc once from capacity 2 to 4 per pair — on a code path about to read gigabytes. The
pre-flight remains two hash lookups per candidate. `impl AsRef<Path>` monomorphises to two instances
(`&Path`, `PathBuf`).

One cosmetic redundancy in the property test: `clumped_bam_output_name` ignores `gzip`, so its primary
set is rebuilt identically for `gzip=true` and `gzip=false` — 8 iterations where 4 would do. Unit-test
noise only, not worth changing.

---

## Errors

No bugs found. Specifics checked and clear: no `unwrap`/index added on a fallible path; the helper
cannot panic (empty slice yields an empty Vec); `from_ref(&chunk[0])` is bounds-safe because
`chunks(2)` guarantees a non-empty chunk; no new I/O ordering (the pre-flight still precedes every
reader open on all five arms).

---

## Structure and tests

### Tests match the file's conventions

All nine follow the established shape: unique `tempdir` tags (verified — all 45 tags in the file are
distinct, and the near-prefix pairs `clump_rep`/`clump_rep_norep` and
`clump_se_alias`/`clump_se_alias_ok` are sibling directories, so the `remove_dir_all`-at-start hazard
B-V5b warned about does not apply), attributable read prefixes (`CR1`/`CN1`/`PX`/…), rejection cases
paired with acceptance cases, acceptance cases asserting content rather than existence, and
`assert_dir_holds_only` for the no-`-o` cases where `assert_rejected_cleanly`'s empty-directory form
cannot apply. `clump_paired_rejects_shared_mate_report_without_output_dir` open-codes the
`!ok` + `DUP_MSG` asserts instead of using the helper, matching
`hardtrim5_rejects_same_basename_across_dirs_without_output_dir`.

The Shape A test does discriminate what it claims: keyed on `chunk[1]`, pair 1 would plan
`y.fq_clumping_report.txt`, which collides with nothing, and the run would proceed.

### `assert_dir_holds_only` message change

Correct and necessary — the helper now serves acceptance-side exact-set assertions, so
"a rejected run must write nothing" would have been a false diagnostic. `assert_rejected_cleanly`
keeps its own version of that wording for the empty-directory case, so nothing is lost.

### Items to consider (none blocking)

**I1 [Medium] — `docs/…/modes/clump-only.md:122` states something this PR proves false.**
"Multi-pair PE runs produce one report per pair, matching v1 SE's one-report-per-input convention."
True for the uBAM shapes (`clump_only.rs:1082`, one report keyed `inputs[0]`), **false for paired
FASTQ**, which writes one report *per mate* (`clump_only.rs:534-548` — and my live run produced both
`a/reads.fq_clumping_report.txt` and `b/reads.fq_clumping_report.txt`). The sentence sits in the
generic "Reorder report" section with no shape qualifier. Pre-existing, so out of the diff's strict
scope — but the per-mate paired-FASTQ report topology is the entire subject of #391, and a reader
reaching for the docs to understand the new rejection lands on the one sentence that contradicts it.
Cheapest fix, one line: *"Paired FASTQ runs write one report per mate; the uBAM shapes write one per
pair, matching v1 SE's one-report-per-input convention."* Precedent exists (no prose docs were added
for #383/#385/#388's rejections, so **not** adding rejection docs here is consistent — this is about
the false sentence only).

**I2 [Low] — property-test coverage: `clumped_paired_bam_output_name` missing from `primary_sets`.**
`io.rs:1251-1264` covers `single_end_output_name`, `clumped_output_name`, `clumped_bam_output_name`,
but not the Shape A primary namer — and that test carries the direction that licenses "hash primaries
alone" (distinct primaries ⇒ distinct secondaries). The sibling test at `:1368` does cover it, but
only in the coarser-than direction. The gap is latent rather than live, because
`clumped_paired_bam_output_name` currently ignores its `_input_r2` parameter and is therefore
identical in shape to `clumped_bam_output_name`; if anyone ever makes it *use* R2 (a combined stem,
say), Shape A's A3 licence would silently lapse. One line closes it.

**I3 [Low] — stale doc comment on the test that changed most.**
`io.rs:1211-1227` still reads "Single-end naming only: paired `_val_N` primaries break this (#388)"
and says nothing about now covering the clump namers or about the switch from `PathBuf` equality to
`collision_key`. The *sibling* test's doc got its #391 note (`:1310-1315`); this one did not. Strictly
speaking "single-end naming only" is still true of the namers tested, so this is an omission rather
than an error.

**I4 [Low] — `clump_paired_rejects_cross_pair_report_collision` is not attributable.**
It asserts `DUP_MSG` only. The plan's reason for skipping a filename pin was candidate order deciding
*which* report gets flagged — but the class suffix is order-independent: whichever report collides,
`_clumping_report.txt` appears. Adding `assert!(stderr.contains("_clumping_report.txt"))` buys
attribution without reintroducing the fragility. Trivial fix.

**I5 [Low] — no acceptance sibling on the Shape A dispatch path.**
The module doc (`tests/integration_output_collision.rs:11-14`) states the file's invariant: "Every
rejection case here is paired with an acceptance case on the same dispatch path and output format."
`clump_paired_bam_rejects_report_that_aliases_an_input` has none. The path *is* covered elsewhere —
`integration_clump_only_ubam.rs::multi_pair_pe_bam_produces_one_output_per_pair` runs two pairs with
reports on — so this is a consistency gap in the file's own contract, not a coverage hole. Either add
the `--no_report_file` sibling or note the cross-file coverage.

**I6 [Low] — comments exceeding the two-line maximum.**
`main.rs:630-632` (Shape B defensive-symmetry note) and `io.rs:1248-1250` (the `#391` note above
`primary_sets`) are three lines each; both condense to two without losing the fact. Trivial fix.

**I7 [Low, optional] — `clump_report_candidates`' doc third sentence is CHANGELOG material.**
"The pre-flight is only as complete as its candidate list: #391 was the paired arm planning primaries
alone while the writer also wrote per-mate reports." — history, recoverable from `git log -L`. That
said, `planned_secondary_outputs` directly above it carries an eight-line doc of exactly this kind, so
keeping it is the *locally consistent* choice. Flagging for the record only; my recommendation is to
leave it.

**I8 [Low] — stale line reference in a comment this diff re-indented.**
`io.rs:1270` says `cli.rs:620`; the multi-input `--basename` guard is at `cli.rs:624`. Pre-existing
text, but the diff moved it, so it is fair drive-by scope. Trivial fix.

**I9 [Low, informational] — the plan's "unit-testable" rationale for the helper is not realised.**
`src/main.rs` has no `#[cfg(test)]` module (verified: `cargo test` reports 0 tests for the `main.rs`
target), so `clump_report_candidates` is only exercised through the integration tests. That is
consistent with its neighbours `planned_secondary_outputs` and `planned_hardtrim_outputs`, so I would
**not** move it into the library just for testability — recording it so the plan's phrasing isn't
mistaken for coverage that exists.

**I10 [Low, informational] — a justification was dropped at the clump-paired call site.**
The replaced comment explained why the hint was `None` ("`clumped_paired_output_names` uses
`input.parent()`, not the CWD" — i.e. `CWD_OUTPUT_HINT` would be false here). A future reader of
`main.rs:550-552` can no longer see why the *other* hint constant is wrong for this arm. The two
constants' doc comments carry part of it (`:55-56`, `:64`). Optional half-line.

### Duplication note (no action recommended)

There are now three spellings of "collect the report paths a run plans":
`planned_secondary_outputs` (SE trim, returns `Result`, also handles demux), the two inline
`if !cli.no_report_file { for input in [&chunk[0], &chunk[1]] { … } }` blocks at `main.rs:827` and
`:1931` (paired trim FASTQ and uBAM, from #388), and now `clump_report_candidates`. Unifying them
would touch #388's shipped paths for no behavioural gain; the right time is whenever those two inline
blocks next need editing.

---

## Out-of-scope observations (surface, don't fix here)

- **`--clump_only`'s `--help` text is stale.** `cli.rs:195` lists `--output-format ubam` among the
  flags "rejected" under `--clump_only`, and `cli.rs:206` says "v1 is FASTQ in / FASTQ out only. uBAM
  in/out is a natural follow-up." Both predate the shipped uBAM arms this PR just extended. User-facing
  text; worth its own issue.
- **FastQC side-outputs remain outside every candidate list.** Under `--clump_only --fastqc`, an input
  named `s_clumped_fastqc.zip` can still be overwritten. Pre-existing and already acknowledged as a
  known residual in `io.rs:1222-1227`; unchanged by this diff.

---

## Recommendations, by priority

**Critical / High:** none.

**Medium**

- **R1 (I1)** — fix `docs/…/modes/clump-only.md:122`'s "one report per pair" sentence, which this
  change proves false for paired FASTQ. One line, in this PR or as a filed follow-up. *Trivial fix.*

**Low**

- **R2 (Logic §6)** — add the over-rejection acceptance test the family's own model provides
  (`paired_report_candidates_do_not_over_reject` at `tests/integration_output_collision.rs:1006`):
  same-filename mates, **no** `-o`, expect exit 0 with both primaries in R1's directory and one report
  in *each* mate's directory. I verified this behaviour by hand; nothing in CI pins it, and the
  per-mate report-directory rule is currently pinned only from the rejection side. Highest-value item
  in this list.
- **R3 (I2)** — add `clumped_paired_bam_output_name` to `primary_sets` in
  `distinct_primary_outputs_imply_distinct_secondary_outputs`. One line.
- **R4 (I4)** — add the order-independent `_clumping_report.txt` substring assertion to
  `clump_paired_rejects_cross_pair_report_collision`. *Trivial fix.*
- **R5 (I3)** — refresh `io.rs:1211-1227`'s doc comment: clump namers now covered, comparison now via
  `collision_key`. *Trivial fix.*
- **R6 (I6, I8)** — condense `main.rs:630-632` and `io.rs:1248-1250` to two lines; correct
  `cli.rs:620` → `cli.rs:624` at `io.rs:1270`. *Trivial fixes.*
- **R7 (I5)** — either add the Shape A `--no_report_file` acceptance sibling or note the cross-file
  coverage, so the module doc's stated invariant stays literally true.

**Informational (no action):** I7, I9, I10, and both out-of-scope observations.
