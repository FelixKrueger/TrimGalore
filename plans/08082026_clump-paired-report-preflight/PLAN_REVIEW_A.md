# Plan Review A — `--clump_only --paired` clumping reports join the collision pre-flight (#391)

**Reviewer:** A (independent)
**Plan:** `plans/08082026_clump-paired-report-preflight/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev` @ `03d793f` (verified: `git rev-parse HEAD` = `03d793f76f9993b43abf506f1dc90d870f1db08d`)
**Date:** 2026-08-08
**Scope note:** Task 2 (reports join the four sibling clump arms' pre-flights) is a maintainer-confirmed IN-scope decision and is NOT relitigated here.

**Verdict: approve with changes.** The diagnosis is correct, the mechanism is the right one, and every load-bearing claim I could check against the code holds — including the A3 derivation, which I re-derived independently rather than accepting. The gaps are all in the *validation* half of the plan, not the fix: one of Task 2's two claimed representative tests is never actually specified, the A3 derivation that licenses skipping three whole arms is never machine-checked despite the repo already owning the exact test scaffolding for it, and the plan's `-o`-only test set misses a filesystem-independent collision shape that its own analysis does not predict. Two factual over-claims (Shape B exposure; "duplicate inputs are rejected earlier") need correcting before they reach the PR body, and step 5's prescribed comment fix is literally wrong as written.

---

## 1. Logic review

### 1.1 Claims verified against code (all correct)

| Plan claim | Verified at | Result |
|---|---|---|
| Dispatch arm | `src/main.rs:526-550` | ✅ exact |
| Driver + pre-flight loop | `src/main.rs:2476-2518`, namer called at `2490` only | ✅ exact |
| Writer + gate | `src/clump_only.rs:532-549`, gate at `534` | ✅ exact |
| Namer keys on full filename, no `--basename` | `src/io.rs:469-482` | ✅ exact |
| Pre-flight semantics | `src/io.rs:96-130` | ✅ exact |
| `--clock` / `--implicon` sites | `src/main.rs:462-476`, `477-501` | ✅ exact |
| Four sibling clump arms | `main.rs:552-558`, `674-678`, `621-630`, `587-597` | ✅ exact |
| BAM report keyed on `inputs[0]` | `src/clump_only.rs:1080-1087` | ✅ exact |
| Hint constants | `src/main.rs:55-66` | ✅ (see §1.4) |
| #388 pattern | `src/main.rs:770-778` (plan says 771; comment opens at 770) | ✅ trivial off-by-one |

### 1.2 The namer really is pre-flight-only — no behavioural leak

Confirmed three ways, so widening `NameFn`'s return to `Vec<PathBuf>` cannot change a single written byte:

1. `output_names` is invoked exactly once in `run_specialty_paired`, at `main.rs:2490`, and its result is pushed into `planned` which flows only into `preflight_output_collisions` (`main.rs:2494`). It is never passed to `run_pair`.
2. `clump_only_paired` re-derives its own paths at `clump_only.rs:415-416` (`naming::clumped_paired_output_names(...)`) and its report paths at `535-536`. Same for `clock`/`implicon`, which take `(output_dir, gzip)` and name internally.
3. `run_specialty_paired` has exactly three call sites, all in `main.rs` (`463`, `478`, `527`), and no test references the signature (`tests/integration_clump_only_ubam.rs:409` mentions it only in a doc comment). So the change is genuinely `main.rs`-internal, as §Integration claims.

The proposed closure compiles: `&cli`, `basename` (`cli.basename.as_deref()`, `main.rs:507`) and the namer's `cli.no_report_file` are all *immutable* borrows of `cli`, and the run closure already borrows `cli.cores`/`cli.compression`/`cli.fastqc`/`cli.no_report_file` today (`main.rs:541-546`). No borrow conflict.

### 1.3 The `!no_report_file` gate matches the writers exactly

All four writer gates are `if !no_report_file`, with no inner condition:

- `clump_only.rs:365` — SE FASTQ, one report keyed on `input`
- `clump_only.rs:534` — paired FASTQ, **both** mates inside one gate, keyed on `input_r1` / `input_r2`
- `clump_only.rs:836` — SE BAM, one report keyed on `input`
- `clump_only.rs:1082` — paired BAM (Shapes A and B), **one** report keyed on `inputs[0]`

The plan's candidate additions mirror each of these one-for-one, including Shape B's `&cli.input[0]` ≡ `inputs[0]` when `inputs.len() == 1`. Gate parity is exact. `grep -rn clumping_report_name src/` confirms these are the *only* four write sites — in particular `--clumpify` on the trim path writes no clumping report, so the plan's enumeration of arms is complete and there is no fifth site hiding on the trim paths.

### 1.4 A3 (uBAM arms cannot collide report-vs-report) — airtight, and I re-derived it

The plan asked me not to inherit its argument, so I rebuilt it. The invariant needed is:

> report-path collision ⟹ primary-path collision (which the existing pre-flight already rejects).

It holds because for every one of the four arms, report and primary are functions of **the same input path**, with **the same directory rule**:

| Arm | Report dir | Primary dir | Report key | Primary key |
|---|---|---|---|---|
| SE FASTQ | `-o` else `input.parent()` | `-o` else `input.parent()` (`io.rs:430-433`… `clumped_output_name`) | `file_name(input)` | `strip_fastq_extensions(input)` |
| SE BAM | same | same (`io.rs:430-433`) | `file_name(input)` | `strip_fastq_extensions(input)` |
| Paired BAM Shape A | `-o` else `chunk[0].parent()` | `-o` else `input_r1.parent()` (`io.rs:456-458`) | `file_name(chunk[0])` | `strip_fastq_extensions(chunk[0])` |
| Paired BAM Shape B | N=1 — see §1.5 | | | |

And `strip_fastq_extensions` is **fold-equivariant**: its gzip-family strip uses `strip_suffix_ignore_ascii_case`, its `.fastq`/`.fq` strip likewise, and its fallback `Path::file_stem()` splits at the last `.` — a position, not a character class — so ASCII case-folding commutes with all three. Therefore fold-equal filenames ⟹ fold-equal stems ⟹ (same dir, same suffix) fold-equal primaries. `gzip` is global (`main.rs:351`: `let gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0]);`) so the SE-FASTQ suffix is uniform across inputs. Cross-shape collision (report vs primary) is impossible: `_clumping_report.txt` cannot equal `_clumped.fq{,.gz}` or `_clumped.bam`.

The `--basename` and mixed-extension shapes the lead flagged are covered: with `--basename` every primary in a multi-input SE run collapses to one path (and `cli.rs:624-628` rejects multi-input SE + `--basename` *before* that anyway), and every multi-pair Shape A primary collapses to `foo_clumped.bam`. Mixed extensions (`x.fq` vs `x.fastq`) move the *primary* key coarser, not the report key — the wrong direction to create an uncaught report collision. **A3 stands.** But see §4.2: it stands as prose, and nothing in CI holds it there.

### 1.5 ⚠️ Over-claim: paired-BAM **Shape B cannot be exposed at all**

The plan states, in §Context and again in §Behavior 5, that "all four sibling clump arms **ARE** exposed to report-vs-INPUT overwrites". For Shape B that is false, and the code proves it:

- Shape B is reached only when `cli.input.len() == 1` (`main.rs:581`).
- `--passthrough` and `--demux` — the only other contributors to `guarded_inputs` (`main.rs:70-79`) — are hard-rejected under `--output-format ubam` (`cli.rs:588`, `cli.rs:602`).
- So `guarded_inputs` is exactly `[cli.input[0]]`, and the report path is that same input's filename **plus a suffix** — strictly longer, hence never fold-equal to it.

Shape B's *existing* one-candidate pre-flight (`main.rs:587-597`) is therefore already unreachable, and Task 2's addition there is purely defensive symmetry. Keeping the edit is fine (five sites, one shape, easier to review than four-plus-an-exception), but the plan, the CHANGELOG and the PR body must not claim it closes a hole. This matters because it is currently the stated *justification* for the edit and for skipping a Shape B test.

### 1.6 ⚠️ Step 5's prescribed comment fix is wrong as written

The plan says the two doc comments at `main.rs:55-62` "are swapped … Swap them." They are not symmetrically swapped. Only the first is wrong:

```rust
/// Replaces the generic advice for the modes that name output into the CWD
/// (`--hardtrim5/3`, `--clock`, `--implicon`), where "use `--output_dir`" is false.
const PAIRED_REPORT_HINT: &str = …          // ← wrong: describes the OTHER constant

/// CWD-output modes' remediation; see `PAIRED_REPORT_HINT` for the paired sites.
const CWD_OUTPUT_HINT: &str = …             // ← already correct for its constant
```

A literal swap yields `PAIRED_REPORT_HINT` documented as "*see `PAIRED_REPORT_HINT` for the paired sites*" — self-referential nonsense. **Rewrite the first comment; leave the second alone.** And note this change gives `PAIRED_REPORT_HINT` a *third* user, so the rewritten comment should name them: paired trim FASTQ (`main.rs:784`), paired trim → uBAM (`main.rs:1888`), and now `--clump_only --paired`.

### 1.7 ⚠️ Behavior §6 over-claims what `validate` rejects

The plan says "duplicate inputs are rejected earlier by existing checks (`cli.rs::validate`, `main.rs:519-525`)". `main.rs:519-525` is the N=1 guard, not a duplicate check, and `validate_paired_input` rejects only two things: R1 == R2 *within* a pair, and an *exact* duplicate pair (same R1 **and** same R2). The general duplicate-input scan at `cli.rs:631` is gated on `!self.paired`, so it never runs here. **A file reused as a mate of two different pairs passes validate.** That is not a nit — it is exactly the shape that produces the uncovered collision class in §4.3.

### 1.8 Minor logic notes

- **A5 is half-true.** "The pre-flight error names both colliding paths" holds for the output-vs-**output** branch (`io.rs:120-126`) but not the output-vs-**input** branch (`io.rs:109-115`), which names only `input.display()`. This constrains Task 2's test 5: it can assert the alias filename (because that filename *is* the input) but cannot assert two paths. Fix the assumption's wording so the implementer doesn't write an assertion that can't pass.
- **Task 1 also closes report-vs-input for the paired FASTQ arm**, for free, because `preflight_output_collisions` runs the input check on every candidate. The plan attributes report-vs-input closure exclusively to Task 2, which undersells Task 1 in the CHANGELOG.
- **A4 is benign on this arm, not merely "pre-existing".** `--fastqc` side-outputs derive from `out_r1_path`/`out_r2_path` (`clump_only.rs:551-554`), which carry the `_clumped_1`/`_clumped_2` discriminators — so they inherit distinctness and cannot collide here. Worth saying, since "out of scope" reads like an unquantified risk.
- **`PAIRED_REPORT_HINT`'s "regardless of source directory" clause is loose on this arm.** Without `-o`, clump-paired reports go to *each mate's own* parent (`io.rs:478-481`) while both primaries go to `input_r1.parent()` (`io.rs:407-409`). So two same-named inputs in different directories do **not** collide on reports without `-o`. The actionable half of the hint ("rename one input, or pass `--no_report_file`") is still correct in every case where it prints, so this is cosmetic — but it is the same directory-rule asymmetry that produces §4.3's missing test, so it is worth understanding rather than waving through.
- **Parameter naming:** renaming `output_names` → `planned_outputs` puts it one character from the local `planned` (`main.rs:2488`) in a six-line body. `pair_outputs` reads better.
- **CHANGELOG:** the live section is `### Unreleased` → `#### Bug fixes` (`CHANGELOG.md:4,6`), not "fixed".

## 2. Assumptions

| # | Status | Note |
|---|---|---|
| A1 | ✅ verified | `clumping_report_name(input, output_dir)` — `io.rs:469`; no `basename` parameter exists to honour. |
| A2 | ✅ verified | `grep -c report src/specialty.rs` = 0. `--clock`/`--implicon` plan nothing but their two primaries. |
| A3 | ✅ verified, ❗untested | Re-derived independently (§1.4). Sound today; nothing in CI pins it. See §4.2. |
| A4 | ✅ verified, stronger than stated | Benign on this arm specifically (§1.8). |
| A5 | ⚠️ half-true | Output-vs-input branch names one path only (§1.8). |
| Gate is load-bearing | ✅ agreed | Exactly right, and the writers confirm it (§1.3). |

**Implicit assumptions the plan does not surface:**

1. **`gzip` is global and derived from `cli.input[0]`.** A3 leans on this. True (`main.rs:351`), but it is an assumption about a variable defined 150 lines from the arm it protects. Worth naming.
2. **`--output_dir` is the only thing that can put two mates' reports in one directory on this arm.** Not stated, and as §4.3 shows, not true — two pairs can share a mate's directory without `-o`.
3. **Shape A's report keys on `chunk[0]`, never `chunk[1]`.** True (`clump_only.rs:1083`: `let report_input = &inputs[0];`), and the plan's step 4 gets it right — but nothing in the proposed test set would catch getting it wrong. See §4.1.
4. **A `.txt`-named file with FASTQ content is accepted as FASTQ input.** Task 2's test 5 depends on this. It holds — `format::detect_input_format` peeks the first byte (`@` → plain FASTQ) and ignores the name — but it is the load-bearing premise of the test and should be written down.

## 3. Efficiency analysis

Nothing to flag. The candidate list grows from 2P to at most 4P paths for paired FASTQ and by ≤1 per input/pair elsewhere; `preflight_output_collisions` is two `HashMap` passes, so total cost stays O(paths) with the map pre-sized (`io.rs:104`). No new I/O — the pre-flight is pure path arithmetic and still runs before any reader opens (`main.rs:2494` precedes the `run_pair` loop).

The one new cost is one `Vec<PathBuf>` allocation per pair where a tuple used to be returned, i.e. P small allocations for P pairs of a run that is about to read gigabytes. Irrelevant, and the alternative (an out-parameter `&mut Vec<PathBuf>`) would trade a measurable nothing for a worse signature. Keep the `Vec`.

## 4. Validation sufficiency

The four Task-1 tests are well designed — in particular tests 1 and 2 are a genuine matched pair (same fixtures, differing only in `--no_report_file`), and test 1's assertion on the *report* filename in stderr is what makes it prove the report path rather than a primary. I checked test 1's shape end to end: `a/reads.fq` + `b/reads.fq` with `-o out` gives primaries `out/reads_clumped_1.fq` / `out/reads_clumped_2.fq` (distinct), so the reports are the *only* collision. That is the right isolation. `assert_rejected_cleanly` will work for it because the inputs live in `a/` and `b/`, not in `out`.

Four gaps, in descending severity.

### 4.1 ❗ Task 2's paired-BAM Shape A change ships untested

Step 6 says "One representative per shape family: SE FASTQ + paired-BAM Shape A", but the enumerated test list contains only `clump_se_rejects_report_that_aliases_an_input` and its `--no_report_file` sibling, and Validation row 5 lists only that test. **No Shape A test is actually specified.** The prose and the deliverables disagree.

This is the gap that matters most, because Shape A is the one arm whose edit sits inside a `for chunk in cli.input.chunks(2)` loop and must key on `chunk[0]` rather than `chunk[1]`. Push `chunk[1]`'s report by mistake and you get a pre-flight that guards a path the run never writes while leaving the path it *does* write unguarded — a silent, permanent hole, and not one test in the proposed set fires. Specify it:

```
clump_paired_bam_rejects_report_that_aliases_an_input
  inputs: a/x.fq  a/y.fq  a/x.fq_clumping_report.txt  a/z.fq   (all valid FASTQ)
  args:   --clump_only --paired --output-format ubam            (no -o)
  expect: rejected, ALIAS_MSG, all four inputs intact
```
Pair 1's report is `a/x.fq_clumping_report.txt` = input 3. Keying on `chunk[1]` instead would plan `a/y.fq_clumping_report.txt`, collide with nothing, and the run would proceed — so this test discriminates the exact mistake. Shape A accepts FASTQ input only (BAM pairs are rejected upstream by `reject_bam_format_mismatch_in_pair`), so plain FASTQ fixtures are correct here.

### 4.2 ❗ A3 licenses skipping three arms and is never machine-checked

A3 is the load-bearing argument that SE FASTQ, SE BAM and Shape B need no report-vs-report coverage. It is sound today (§1.4). But it is a naming invariant, and naming invariants in this repo have broken before — #382 widened `strip_fastq_extensions` and turned a naming inconsistency into data loss; #388 is *literally* the case where a primary key stopped being coarser than a report key. If a future change makes `clumped_bam_output_name` honour something `clumping_report_name` doesn't, A3 dies silently and no test in the plan notices.

The repo already owns the scaffolding. `src/io.rs:1288-1320`, `primary_output_key_is_coarser_than_secondary_keys`, does exactly this job for SE trim naming — and its doc comment already says *"paired `_val_N` inverts this, #388"*, i.e. the author of that test understood precisely the invariant #391 is about. Extend it (or add a sibling) over the clump namers:

- for each pair of inputs drawn from a table spanning `{.fq, .fastq, .fq.gz, .bgz, mixed case, same name/different dir}` × `{basename None/Some}` × `{gzip true/false}` × `{output_dir None/Some}`,
- assert: `collision_key(clumping_report_name(a, o)) == collision_key(clumping_report_name(b, o))` **implies** `collision_key(primary(a)) == collision_key(primary(b))`, for `primary ∈ {clumped_output_name, clumped_bam_output_name, clumped_paired_bam_output_name}`.

Note the existing test uses case-**sensitive** `assert_ne!` on `PathBuf`, which is weaker than the pre-flight's `collision_key`. Use `collision_key` so the fold dimension is actually exercised. This turns A3 from a paragraph a reviewer has to re-derive into a line CI holds — cheap, and it is the single highest-value addition to this plan.

### 4.3 ❗ A collision class the plan's analysis does not predict, and no test covers

Every proposed test passes `-o out`. That hides the directory-rule asymmetry from §1.8: reports use **each mate's own** `input.parent()`, primaries use **R1's** `input_r1.parent()` for both mates. Combined with §1.7 (a file may be the mate of two different pairs), there is a **case-free, filesystem-independent** collision that the fix catches and the plan never mentions:

```
pair 1: p/a.fq  d/x.fq
pair 2: q/b.fq  d/x.fq        # same R2 in both pairs — validate permits this
args:   --clump_only --paired  (no -o)

primaries: p/a_clumped_1.fq  p/x_clumped_2.fq  q/b_clumped_1.fq  q/x_clumped_2.fq   → all distinct
reports:   p/a.fq_…  d/x.fq_…  q/b.fq_…  d/x.fq_…                                   → COLLIDE
```

Today: exit 0, one report silently lost. With the fix: rejected. This is strictly better than the plan's test 4 as a "candidate list spans pairs" proof, because it needs no `-o` (so it also pins the report directory rule) and no case-folding (so it behaves identically on APFS and ext4 — test 4's cousin `paired_rejects_fold_equal_filenames_into_shared_output_dir` needed a filesystem caveat for exactly this reason). Add it alongside test 4, not instead of it.

While you are here: correct Behavior §6's claim about duplicate inputs (§1.7), since this test is a live counterexample to it.

### 4.4 Smaller validation items

- **Extend validation 6 (the expected-fail control) to all three rejection tests**, not just test 1. I checked each pre-patch: test 1 exits 0 (reports collide, no primary does), test 4 exits 0 (all four primaries distinct), test 5 exits 0 *and destroys an input*. All three genuinely fail pre-patch, so the control is cheap and it is the plan's own stated defence against fixture-selection blindness. Applying it to one of three tests is where that defence historically leaks.
- **Add a `--basename` variant of test 1.** A1 asserts `--basename` runs are "exposed identically" and fixed by the same candidates; nothing tests it. Three extra lines: `--basename foo` makes the primaries `out/foo_clumped_1.fq` / `out/foo_clumped_2.fq` — *guaranteed* distinct by construction — while the reports still collide. That is a harder isolation of the report path than test 1 achieves, since a reader no longer has to reason about stems at all.
- **Point test 5 at the existing model.** `assert_rejected_cleanly` cannot be used (it demands an empty directory, and test 5's inputs live in the directory being checked). The pattern to copy is `se_trim_rejects_output_that_aliases_a_report_input` (`tests/integration_output_collision.rs:674-691`): assert `ALIAS_MSG`, then `count_reads_from(&alias_file, "REPORTY") == 40` for content-verified intactness. Keeping the alias file valid FASTQ is essential and the plan already says so — without it, a future regression that removed the pre-flight would still fail the run (on a parse error) and the test would pass for the wrong reason.
- **Test 4 must not assert a filename.** Candidate order is `[p1_primary, p1_primary, p1_rep_r1, p1_rep_r2, p2_primary, …]`, so the first duplicate found is pair 2's R1 report — `y.fq_clumping_report.txt`, not `x.fq`. The plan correctly asserts only "rejected"; flagging so an implementer doesn't "improve" it into a flaky assertion.
- **No existing test regresses.** I checked the ones that could: `pe_byte_identity` (`tests/integration_clump_only.rs:141-171`) uses distinct `_R1`/`_R2` filenames; `multi_pair_pe_bam_produces_one_output_per_pair` (`integration_clump_only_ubam.rs:477-518`) copies `B_R{1,2}.fastq.gz` *into* its `-o` dir, but neither new report candidate equals an input, so no over-rejection; `pe_bam_collision_preflight_case_folded` (`:630`) still trips on pair 2's *primary* (which precedes its report in the candidate order), so its `contains("collision")` assertion holds. **And the hint change is safe:** the only tests pinning hint text are `integration_output_collision.rs:625` and `:661`, both on `--hardtrim5`/`--clock` (`CWD_OUTPUT_HINT`), and the two negative assertions on `"different source directories"` are on those same tests. Nothing pins `GENERIC_ADVICE` on any clump path. Validation 7's grep is a good instinct but it will come back clean.

## 5. Alternatives

**Keep the chosen shape.** Widening `NameFn` to `-> Vec<PathBuf>` is the right call and the plan's rationale is sound. It also matches how the trim arm already does this job — `main.rs:747` builds a local `candidates: Vec<PathBuf>` and extends `planned` — so after the change the two paired pre-flights read the same way. The rejected alternative (a second pre-flight inside the clump arm) is correctly rejected: it would run the check twice with two different hints and duplicate the chunking.

Two alternatives worth a decision rather than silence:

**(a) A shared `clump_report_candidates` helper instead of five inline blocks.** Task 1 + Task 2 plant the same `if !cli.no_report_file { push(clumping_report_name(…)) }` shape at five sites. `main.rs` already established the opposite convention for this exact job: `planned_secondary_outputs` (`main.rs:89-110`) and `planned_hardtrim_outputs` (`main.rs:113-130`) are named helpers precisely so the candidate-list↔writer parity is auditable in one place — and `planned_secondary_outputs`' doc comment spells out why. A `fn clump_report_candidates(cli: &Cli, inputs: &[PathBuf], output_dir: Option<&Path>) -> Vec<PathBuf>` would collapse all five (SE FASTQ passes `&cli.input`, SE BAM the same, Shape A passes `&chunk[..1]`, Shape B `&cli.input`, paired FASTQ `&[r1, r2]`). Trade-off: five two-line blocks are individually trivial to review, and the helper adds an indirection whose "which inputs" parameter is itself a place to get `chunk[0]` vs `chunk[1]` wrong. My recommendation: **take the helper**, because the writers it must stay in step with are already scattered across four functions in `clump_only.rs`, and one call site per arm is easier to diff against them than five hand-rolled copies — but it is a genuine judgement call and either answer is defensible.

**(b) Derive `clumping_report_name` from the primary output path.** `demux::demux_base_name` already takes this approach (deriving from the primary, not the input), which is why demux's secondaries inherit primary distinctness structurally. Applied here it would make report collisions *impossible* wherever primaries differ, retiring #391, #388's clump twin, and A3 in one stroke. **Reject it**, and the reason is worth recording: it changes report filenames, which are a published contract — `io.rs:463-468` documents that `*_clumping_report.txt` is deliberately shaped so nf-core/rnaseq's MultiQC `*_trimming_report.*` glob does not pick it up, and `--clumpify`'s per-input report convention depends on the current spelling. Not worth it for a pre-flight fix.

## 6. Action items

### Critical

1. **Specify the paired-BAM Shape A test.** Step 6's prose claims it as a representative; the deliverable list omits it, so Task 2's most error-prone edit would ship with zero coverage. Concrete shape in §4.1 — it discriminates the `chunk[0]` vs `chunk[1]` mistake, which nothing else in the plan does.

### Important

2. **Add the `collision_key`-based property test for A3** (§4.2), extending or siring a sibling to `src/io.rs:1288`'s `primary_output_key_is_coarser_than_secondary_keys`. A3 is what licenses skipping three arms; right now it lives only in the plan. Use `collision_key`, not `PathBuf` equality, so the fold dimension is real. Highest value-per-line item in this review.
3. **Correct the Shape B exposure claim** in §Context and §Behavior 5 (§1.5): Shape B is `N=1` with `--passthrough`/`--demux` rejected, so its report can never alias an input and its pre-flight is already unreachable. Keep the edit as defensive symmetry if you prefer, but say that — the current wording will otherwise land in the CHANGELOG and PR body as a fixed hole that never existed.
4. **Fix step 5's instruction** (§1.6): rewrite `PAIRED_REPORT_HINT`'s doc comment; do not swap the two. A swap makes it self-referential. The rewrite should list all three users, including the new `--clump_only --paired` site.
5. **Add the no-`-o` cross-pair test from §4.3.** Case-free, filesystem-independent, and it covers a collision class the plan's `-o`-only analysis does not predict — namely that reports use per-mate `input.parent()` while primaries use R1's.
6. **Extend validation 6's expected-fail control to tests 4 and 5** (§4.4). All three demonstrably exit 0 pre-patch, so the control costs one extra run each and is the plan's own named defence against fixture-selection blindness.
7. **Correct Behavior §6** on duplicate inputs (§1.7): only within-pair `R1 == R2` and *exact* duplicate pairs are rejected; a file reused across two different pairs passes `validate`, and item 5's test relies on that.
8. **Fix A5** (§1.8): the output-vs-input branch names only the input path. State it, so test 5's assertions are written against what the code actually prints.

### Optional

9. Add the `--basename` variant of test 1 (§4.4) — three lines, and it isolates the report path more sharply than test 1 does.
10. Decide alternative (a): one `clump_report_candidates` helper vs five inline blocks (§5). I lean helper; either is defensible, but the plan should choose deliberately rather than default into five copies.
11. Note in §Behavior/CHANGELOG that Task 1 also closes report-vs-input for the paired FASTQ arm for free (§1.8) — a real part of the fix that currently goes unadvertised.
12. Sharpen A4: FastQC side-outputs on this arm inherit the `_clumped_{1,2}` discriminators and so cannot collide (§1.8). "Out of scope" reads like an open risk; it isn't one here.
13. Record the two unstated premises: `gzip` is global from `cli.input[0]` (`main.rs:351`), and content-based format detection accepts FASTQ content in a `.txt`-named file (which is what makes test 5 possible).
14. Rename the widened parameter `pair_outputs` rather than `planned_outputs` — one character from the local `planned` (§1.8).
15. CHANGELOG heading is `#### Bug fixes` under `### Unreleased`, not "fixed" (`CHANGELOG.md:4,6`).
16. Plan says `main.rs:771-778` for the #388 pattern; the comment opens at 770.
