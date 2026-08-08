# Plan Review B — `--clump_only --paired` clumping reports join the collision pre-flight (#391)

**Reviewer:** B (independent; fresh context)
**Plan:** `plans/08082026_clump-paired-report-preflight/PLAN.md`
**Base verified at:** `dev` @ `03d793f` (`fix(cli): reject --hardtrim5 with --hardtrim3 (#386) (#393)`) — matches the plan's stated base.
**Date:** 2026-08-08

**Verdict up front:** the diagnosis is correct, the fix shape is correct, and every code citation in the plan checks out against the tree. I found **no Critical defects** — nothing that would make the patch wrong or the issue stay open. What I did find is a cluster of **Important** issues concentrated in two places: (a) the user-facing hint text the plan adopts is partly false for this arm, and (b) the validation section over-promises relative to what it enumerates, and one of its tests cannot be written with the helper the plan names. Plus a high-value, low-cost validation upgrade the plan misses.

---

## 0. Claim-by-claim verification

Every line/behaviour claim in the plan was checked against the tree. Result: **all verified, one sub-argument factually wrong (harmless conclusion).**

| Plan claim | Verdict |
|---|---|
| Dispatch at `main.rs:526-550`, namer returns only the two primaries | ✅ exact |
| `run_specialty_paired` at `main.rs:2476-2518`; `NameFn: FnMut(&Path,&Path) -> (PathBuf,PathBuf)` | ✅ exact (bound at 2484) |
| Namer used **only** in the pre-flight (2489-2493); real paths re-derived in `clump_only.rs:415-416` | ✅ verified — see §2.1 |
| Writer `clump_only.rs:532-549`, gated `!no_report_file`, keys `clumping_report_name(input_rN, output_dir)` | ✅ exact (gate at 534, names at 535-536) |
| `clumping_report_name` (`io.rs:469-482`) ignores `--basename` | ✅ signature is `(input, output_dir)`; body uses `file_name()` verbatim |
| `preflight_output_collisions` (`io.rs:96-130`): input check first, then dup check; names both paths; hint replaces `GENERIC_ADVICE` | ✅ exact |
| `--clock`/`--implicon` write nothing but two primaries | ✅ `grep -n "report\|fastqc" src/specialty.rs` → **zero matches** |
| uBAM paired report keyed on `inputs[0]` (`clump_only.rs:1080-1087`) | ✅ exact |
| Task 2 sites: SE FASTQ 552-558, SE BAM 674-678, Shape A 621-630, Shape B 587-597 | ✅ all four exact |
| `PAIRED_REPORT_HINT` / `CWD_OUTPUT_HINT` doc comments are mismatched (`main.rs:55-66`) | ✅ confirmed — the "modes that name output into the CWD (`--hardtrim5/3`, `--clock`, `--implicon`)" prose sits on `PAIRED_REPORT_HINT` |
| A3's premise "extension stripping is case-insensitive since #390" | ✅ `strip_fastq_extensions` (`io.rs:540`) uses `strip_suffix_ignore_ascii_case` for both suffix groups; the `.bam` fallback is `Path::file_stem`, which splits positionally — so the function is **fold-equivariant**, which is what A3 actually needs |
| A3's sub-argument "with `--basename`, multi-pair primaries collide immediately" | ❌ **wrong mechanism** — see A-4 |

Two additional facts I established that the plan does not state, both favourable:

- **No existing test breaks.** `grep` for `GENERIC_ADVICE`'s text ("distinct output paths") across `tests/` → zero hits. The only hint-text assertions in the suite are `stderr.contains("current working directory")` at `tests/integration_output_collision.rs:625` and `:661`, which pin `CWD_OUTPUT_HINT` on the hardtrim/clock paths and are untouched. Validation 7's grep is prudent but will come up empty — the plan can say so.
- **The only binary-driven `--clump_only --paired` test today is `tests/integration_clump_only.rs:141 pe_byte_identity`**, and it uses `BS-seq_10K_R1/R2.fastq.gz` — distinct filenames, so distinct reports. Not newly rejected. `tests/integration_output_collision.rs` contains **zero** `clump` matches today, so the plan is adding the first collision coverage for this mode, not extending existing coverage.

---

## 1. Logic review

### 1.1 The diagnosis is right, and the fix is complete for the issue as filed

Reports collide iff two candidates share a `collision_key`, i.e. same directory + fold-equal filename. `clumping_report_name` is injective in `file_name()`, and `clumped_paired_output_names` appends `_clumped_1` / `_clumped_2` — so the primary key is *finer* than the report key on this arm, which is exactly the inversion #388 documented for `_val_N`. Adding the two gated report paths to the candidate list closes it.

I checked the two reachable collision shapes and both are covered by the plan's candidate list:

1. **With `--output_dir`** (the issue's repro): `-o out a/reads.fq b/reads.fq` → both reports land on `out/reads.fq_clumping_report.txt`. Rejected. ✅
2. **Without `--output_dir`, via a case-only alias**: `a/x.fq a/X.fq`. `path_identity_key` is case-**sensitive** (`io.rs:78`), so `Cli::validate`'s same-file and duplicate-pair checks let this through; reports become `a/x.fq_clumping_report.txt` and `a/X.fq_clumping_report.txt` → fold-equal → rejected. Primaries (`x_clumped_1.fq` / `X_clumped_2.fq`) stay distinct, so nothing else catches it. ✅ The plan does not name this shape; it is covered by construction, but see V-4.

### 1.2 The `!no_report_file` gate placement matches the writer exactly

`main.rs:546` passes `cli.no_report_file` into `clump_only_paired`, which gates at `clump_only.rs:534`. The plan's closure gates on `!cli.no_report_file` — the same value, the same polarity, guarding the same two `clumping_report_name` calls with the same arguments. **Exact parity, including the directory asymmetry** noted in §1.3. No drift.

### 1.3 An asymmetry worth recording (not a defect in the plan)

On the paired FASTQ arm the report and primary do **not** share a directory rule when `-o` is absent: `clumped_paired_output_names` puts *both* primaries in `input_r1.parent()` (`io.rs:407-409`), while `clumping_report_name(input_r2, None)` puts R2's report in `input_r2.parent()`. So `trim_galore --clump_only --paired A/r1.fq B/r2.fq` writes both outputs into `A/` but R2's report into `B/`.

The plan's fix mirrors the writer exactly, so this causes no bug — but it matters twice:

- it is why A3's "same directory rule" argument is scoped to the SE/uBAM arms only (correct as written);
- it is why the hint text the plan adopts is inaccurate (see A-1).

The identical asymmetry exists on the #388 trim-paired path (`paired_end_output_names` also resolves to `input_r1.parent()`), so this is a shared, pre-existing quirk and a reasonable follow-up issue, not a blocker.

### 1.4 A3 is sound — and Task 2 is provably rejection-neutral because of it

A3's chain holds: report collision ⟹ fold-equal `file_name()` ⟹ (fold-equivariance of `strip_fastq_extensions`, verified) fold-equal stem ⟹ fold-equal primary in the same directory ⟹ already rejected. I probed the shapes the brief asked about:

- **Mixed extensions** (`x.fq` vs `x.FQ`, `x.fq` vs `x.fq.gz`): fold-equal filenames force identical extensions up to case, so this reduces to the case-alias shape and is covered. The *reverse* direction (`x.fq` + `x.fastq` → same stem, different reports) produces a primary collision without a report collision — irrelevant, already rejected.
- **Shape B**: `cli.input.len() == 1`, therefore exactly one report and one primary. Report-vs-report is structurally impossible. ✅
- **Shape A**: one report per pair keyed on `chunk[0]`, primary keyed on the same `chunk[0]` into the same directory. Multi-pair collisions coincide. ✅
- **Report-vs-primary**: a report path always ends `_clumping_report.txt`; every clump primary ends `_clumped.fq`, `_clumped.fq.gz`, `_clumped_N.fq(.gz)` or `_clumped.bam`. Disjoint suffixes, so this cross-kind collision cannot occur on any arm. ✅
- **`--basename`**: unreachable as a multi-input shape — see A-4.

A consequence the plan should state explicitly, because it is the strongest argument for Task 2's safety: **by A3, Task 2 adds zero new output-vs-output rejections.** Every report path it appends can only collide with another report path when a primary already collides (rejected today), and can never collide with a primary. So Task 2's *entire* behavioural delta is the report-vs-input direction it is designed to close. That makes it about as low-risk as a pre-flight widening can be — and it also means A3 is load-bearing for Task 2's blast-radius claim, not just for the decision to skip report-vs-report on the uBAM arms. Which is why I want A3 machine-checked (V-1).

### 1.5 Widening `NameFn` to `Vec<PathBuf>` really is pre-flight-only

Confirmed on all three call sites:

- `run_specialty_paired` calls `output_names` at `main.rs:2490` and nowhere else; the value is pushed into `planned`, consumed by `preflight_output_collisions`, and dropped.
- `clump_only_paired` re-derives via `naming::clumped_paired_output_names` at `clump_only.rs:415-416`.
- `specialty::clock` / `specialty::implicon` receive `(r1, r2, gzip, output_dir, cores, compression)` — no paths passed in, so they re-derive too.

So there is **no behavioural leak**: the widened return type cannot change a single written path. The plan's claim is airtight.

On compilation: the proposed closure captures `cli` immutably while `&cli` is also passed as argument 1 of the same call. That already happens today (the run closure reads `cli.cores`, `cli.fastqc`, `cli.no_report_file`), so it compiles. `naming` is already imported at these sites.

### 1.6 Coverage of the dispatch surface is complete

I enumerated every `--clump_only` dispatch branch (`main.rs:502-697`): paired FASTQ, SE FASTQ, Shape B, Shape A, SE BAM — five arms, matching the plan's five (one Task 1 + four Task 2). They map onto exactly four report-writing sites (`clump_only.rs:366`, `535/536`, `837`, `1084`); Shapes A and B share `1084`. Nothing is missed.

I also checked for a **sixth** exposure the plan might have inherited: `--clumpify` on the trim path. `grep -rn clumping_report_name src/` returns hits only in `io.rs` and `clump_only.rs`, so `--clumpify` writes no clumping report and the trim path's pre-flight (`planned_secondary_outputs`, `main.rs:89-110`) is already complete for its own secondaries after #388. **No fifth site.** The comment at `clump_only.rs:532` ("mirrors `--clumpify`'s per-input report convention") refers to the convention, not to a second writer.

### 1.7 The cross-pair test's arithmetic checks out

Pairs `(a/x.fq, a/y.fq)` and `(b/y.fq, b/x.fq)` with `-o out` give primaries `x_clumped_1`, `y_clumped_2`, `y_clumped_1`, `x_clumped_2` — four distinct names — while reports collide on both `x.fq_…` and `y.fq_…`, and collide **only across pairs** (within each pair the two filenames differ). So the test isolates exactly the property it claims. `validate`'s duplicate-pair check compares `(r1,r2)` key pairs and does not fire. ✅ Good test design.

---

## 2. Assumptions

### A-1 (Important) — `PAIRED_REPORT_HINT` contains a clause that is false on this arm

The plan changes the hint from `None` to `Some(PAIRED_REPORT_HINT)`. I agree the change is *needed* — `GENERIC_ADVICE` is actively unhelpful for the issue's repro, since it suggests "different source directories or `--output_dir`" when the inputs are already in different directories and `--output_dir` is what *caused* the collision. But the constant's text says:

> "Paired outputs and reports are named from the input filename alone, so two inputs sharing a filename collide **regardless of source directory** — rename one input, or pass `--no_report_file` if only the reports collide."

Two problems on the clump-paired arm:

1. **"regardless of source directory" is false without `-o`.** Per §1.3, dropping `--output_dir` genuinely *does* separate the two reports (they follow their own inputs' parents). The hint tells the user the one remedy that works is unavailable.
2. **"Paired outputs … are named from the input filename alone"** is the opposite of the defect: the *outputs* carry `_clumped_1`/`_clumped_2` and are fine; only the *reports* lack a discriminator. A user reading this will look for a primary collision.

The remedy clauses ("rename one input", "`--no_report_file`") are both correct here, including under `--basename` (which affects primaries only). So the damage is confined to the explanatory half.

Recommendation: introduce a third constant rather than reuse `PAIRED_REPORT_HINT`, e.g.

> "Clumping reports are named from the input filename alone (no `_clumped_N` discriminator), so two inputs sharing a filename collide on reports once `--output_dir` collects them in one place — rename one input, drop `--output_dir`, or pass `--no_report_file`."

If the maintainer prefers one constant for both paired sites (defensible — it keeps the message surface small), then fix the wording at the constant, since the "regardless of source directory" clause is equally false at #388's trim-paired site for the same reason. Either way this should be a deliberate decision recorded in the plan, not an inherited string.

### A-2 (Important) — the hint change is untested in both directions

No proposed validation asserts the hint appears. Given A-1, that cuts both ways: without an assertion a future edit can silently revert the hint; with an assertion on the wrong text you lock in the inaccuracy. Recommendation: settle A-1 first, then assert only the durable, actionable fragment (`--no_report_file`) in test 1 — not the whole sentence.

### A-3 (verified, but state the consequence) — A5 is true, and the message is unhelpfully symmetric

`preflight_output_collisions` names both paths — but for a report-vs-report collision `existing` and `p` are *the same string*, so stderr reads:

> `out/reads.fq_clumping_report.txt and out/reads.fq_clumping_report.txt would be written to the same file.`

Combined with a hint that says "rename one input", the user is told to rename an input the message never identifies. This is pre-existing (identical at the #388 site) and fixing it means threading input provenance into the pre-flight (`&[(PathBuf, &Path)]`), which is out of scope for a bug-fix PR. It does not break the plan's test assertion — `reads.fq_clumping_report.txt` is present and is attributable only to a report namer. Worth a follow-up issue and a one-line acknowledgement in the plan's Assumptions, so the next reader does not mistake it for a defect introduced here.

### A-4 (Important) — A3's `--basename` branch cites a mechanism that does not exist

The plan writes: "With `--basename`, multi-pair primaries collide immediately (all `foo_clumped.bam`)." That shape never reaches the pre-flight: `cli.rs:651` rejects `--basename` with more than one paired pair, and `cli.rs:624` rejects it with more than one SE input. So `--basename` is only ever seen with a single input (SE) or a single pair, in which case there is exactly one report and one primary and no collision of any kind is possible.

The conclusion (no silent report loss under `--basename`) survives intact; only the reason is wrong. But A3 is the plan's justification for a *deliberate non-fix* on four arms, so a reviewer inheriting a wrong mechanism is precisely how this rots. Replace the sentence with the actual guard citation.

### A-5 (Optional) — A4 is broader than it needs to be

A4 states FastQC side-outputs stay outside every pre-flight. True, and correctly out of scope — but on the clump paths the residual is *narrower* than the general case: `clump_only.rs:551-554` derives FastQC's inputs from `out_r1_path` / `out_r2_path`, which carry the `_clumped_N` discriminators, so a FastQC output cannot collide where the primaries do not. Sharpening A4 to say so removes a phantom risk from the next reader's list.

### A-6 (Optional) — the `no_report_file` gate rationale is right and worth keeping

The plan's "Fixed vs configurable" note ("a candidate list that includes reports the run won't write would reject legal `--no_report_file` runs") is exactly the correct framing and is the reason test 2 is a genuine control rather than a formality. No change needed; flagging it as a strength.

---

## 3. Efficiency analysis

Nothing to challenge. Concretely:

- Candidate list grows by ≤ 2 per pair (Task 1) and ≤ 1 per input/pair (Task 2). `preflight_output_collisions` is two `HashMap` builds over that list, so total work stays `O(planned + inputs)` with `String` keys of path length. For realistic invocations (tens of files) this is microseconds against a run that reads gigabytes.
- The `Vec<PathBuf>` return replaces a 2-tuple with one heap allocation per pair. Allocations equal to the pair count, on a path that already does per-pair `File::create` and gigabyte-scale I/O. Immeasurable.
- `guarded_inputs(cli)` clones `cli.input` on every call — pre-existing, unchanged, and once per dispatch.
- No new I/O, no new syscalls, no change on any accepted run.

One micro-note if the implementer cares: `planned.extend(planned_outputs(...))` will reallocate as it grows; `Vec::with_capacity(cli.input.len() * 3)` would avoid a couple of reallocs. Not worth the line.

---

## 4. Validation sufficiency

The proposed suite is well-designed in structure — rejection/acceptance pairs on the same dispatch path, a shared-fixture negative control, and an expected-fail-pre-patch step. Four specific gaps, one of them a plan-versus-plan inconsistency.

### V-1 (Important) — A3 is prose-only, and there is an existing test that should carry it

`src/io.rs:1229 distinct_primary_outputs_imply_distinct_secondary_outputs` already asserts, over a matrix of `basename × gzip × output_dir`, that distinct primaries imply distinct secondaries — and it **already includes `clumping_report_name`** in its namer list. But its primaries come from `single_end_output_name` only. Adding `clumped_output_name` and `clumped_bam_output_name` to that primary set (a few lines, no new fixtures, no new test) converts A3 from an argument in a plan file into a CI-enforced invariant. Given §1.4 — A3 is what makes both the uBAM non-fix *and* Task 2's rejection-neutrality true — this is the highest-value item in my review.

While there: `src/io.rs:1287-1291`'s doc comment says "the primary key is strictly coarser than the report key (in single-end naming — paired `_val_N` inverts this, #388)". `_clumped_N` inverts it the same way; add #391 so the next reader finds the second instance.

### V-2 (Important) — test 5 cannot be written with `assert_rejected_cleanly`, and the plan does not say what shape it takes

`assert_rejected_cleanly` (`tests/integration_output_collision.rs:95`) requires the passed directory to be **empty** after the run. For a report-vs-input collision the aliasing input must sit in the directory the report would be written to — i.e. inside `-o out`, or (cleaner) both inputs in one directory with no `-o` at all. Either way that directory contains the inputs, so `assert_rejected_cleanly` fails on its `leftovers.is_empty()` assertion.

The workable shape is the no-`-o` one, which the file already has a helper for:

- fixtures `dir/s.fq` (valid FASTQ) and `dir/s.fq_clumping_report.txt` (also valid FASTQ content, so the rejection cannot come from a format error);
- `s.fq`'s report resolves to `dir/s.fq_clumping_report.txt` == input 2 → `ALIAS_MSG`;
- primaries are `dir/s_clumped.fq` and `dir/s.fq_clumping_report_clumped.fq` — distinct, so nothing else fires;
- assert with `assert_dir_holds_only(&dir, &["s.fq", "s.fq_clumping_report.txt"])`, which is strictly stronger than "both inputs intact" because it also proves no primary was created.

Note `input[0]` must be valid FASTQ regardless: `sanity_check_any(&cli.input[0])` runs at main entry before dispatch. Also confirmed there is **no input-extension whitelist** in `Cli::validate`, so a `.txt`-named input is accepted and the fixture is viable.

### V-3 (Important) — the plan promises two Task 2 representatives and enumerates one

Step 6's last bullet says "One representative per shape family: SE FASTQ + paired-BAM Shape A", but only `clump_se_rejects_report_that_aliases_an_input` is named, and validation row 5 lists only that one. Either write the Shape A test or drop the claim. If written, it is cheap: Shape A takes two FASTQ inputs plus `--output-format ubam` (the existing `se_trim_ubam_rejects_shared_stem` at `:284` establishes that in-test plain FASTQ works fine on uBAM-output arms), and the aliasing input is a file named `<r1-filename>_clumping_report.txt` in the report's directory. My preference is to write it: Shape A is the only Task 2 arm whose report keys on `inputs[0]` rather than "each input", so it is the one whose candidate-list line differs in shape from the others.

### V-4 (Important) — strengthen the acceptance test to bind namer to writer

Test 3 asserts "primaries + both reports present". Upgrade it to `assert_dir_holds_only(&out, &[<the four exact filenames>])`. The helper already exists at `:77`. This costs one line and converts a weak existence check into the invariant whose *absence* is the root cause of this entire bug family: **the pre-flight's candidate list equals the set of files the run writes.** #388 and #391 are both instances of the namer and the writer disagreeing; an exact-set assertion is the only test shape that catches a future third instance. I would rate this the second-highest-value item after V-1.

(It also happens to be the cheapest available answer to Alternative 5.2 below, which is the structurally correct but expensive fix.)

### V-5 (Optional) — two smaller items

- **Test 1 needs two assertions, not one.** `assert_rejected_cleanly` takes a single `expect`. Pass `DUP_MSG` (so the collision *kind* is pinned, via the existing constant) and add a separate `assert!(stderr.contains("reads.fq_clumping_report.txt"))` for the report attribution. The plan currently implies one assertion covering both.
- **Unique `tempdir` tags.** `tempdir(tag)` (`:31`) keys on `std::process::id()`, which is identical for every test in the binary, and it `remove_dir_all`s first. Two tests sharing a tag will destroy each other's fixtures under `cargo test`'s thread-parallel execution. Five new tests need five unused tags. The nastier failure mode is not the obvious one: test 2's "no report files exist" assertion could pass **vacuously** if a same-tag test wiped the directory mid-run — and test 2 is the plan's whole defence against fixture-selection blindness. Worth a line in step 6.

### V-6 (Optional) — the paired arm's report-vs-input direction has no test

Task 1's widened list closes report-vs-input on the paired FASTQ arm too, and that is the more destructive direction (it overwrites a file the run was told to read). Reachable: `--paired dir/x.fq dir/x.fq_clumping_report.txt` with no `-o` → R1's report lands on input 2 while primaries stay distinct. Contrived, but it is one more test of the same shape as V-2's and it covers a different candidate-list construction. At minimum, note in the validation table that it is covered by construction so a coverage audit does not read it as a gap.

### V-7 — could any proposed test pass for the wrong reason?

I checked each against the repo's fixture-selection-blindness history and they hold up well:

- **Test 1**: `assert_rejected_cleanly` pins the `PREFIX` ("Output path collision…"), so a rejection from `validate`, from format detection, or from a sanity check fails the test. `path_identity_key` is case-sensitive and the two inputs differ in their parent, so neither the same-file nor the duplicate-pair check fires. The string `reads.fq_clumping_report.txt` is producible only by `clumping_report_name`. **Sound.**
- **Test 2**: the true negative control — identical fixtures, one flag, and the *only* thing that flag removes from the candidate list is the two report paths. **Sound, and the plan is right to lean on it.**
- **Test 3**: passes today, pre-patch (nothing collides), so it is not a fix-detector — it is an over-rejection detector, which is its stated job. Weak only in its assertions; V-4 fixes that.
- **Test 4**: verified in §1.7 — collisions occur only across pairs, so the cross-pair property is genuinely isolated.
- **Test 5**: sound *if* built per V-2; as sketched with `assert_rejected_cleanly` it will not run at all.
- **Validation 6** (run test 1 against the unpatched build) is the right guard and should be reported in the PR body, not just performed.

One thing no proposed test covers: `gzip` is derived from `is_gzipped(&cli.input[0])` (`main.rs:351`), so all-plain-FASTQ fixtures give `_clumped_1.fq` (no `.gz`) — which matches the filenames the plan asserts. Consistent; just confirming the plan's expected filenames are right.

---

## 5. Alternatives

### 5.1 `Vec<PathBuf>` namer vs. an optional secondary closure — plan's choice is right

Agreed, and for a reason the plan does not give: the tuple return type *encodes the bug*. "Return the two primaries" is a shape that cannot express a mode with three planned outputs, so every future secondary would need another parameter. `Vec<PathBuf>` with the doc contract "every path this pair writes" makes the pre-flight's completeness a property of the type rather than of the caller's diligence. Keep it, and keep the parameter rename (`output_names` → `planned_outputs`) — the name is doing real work. Also extend the function's doc comment's mode list: it currently says "(`--clock`, `--implicon`)" at `main.rs:2470-2471` and has silently gained `--clump_only --paired`.

### 5.2 Share one derivation between pre-flight and writer (structurally correct, expensive)

The residual weakness after this patch: `clumping_report_name` is still called independently in two places under two independently written `!no_report_file` gates. That is the same duplication that produced #391, just with the two copies now agreeing. The durable fix is one derivation — either pass the report paths into `clump_only_paired`, or have both sides call a shared `clump_paired_planned_outputs(...)`.

I do **not** recommend doing this in a bug-fix PR: it touches the writer, and byte-identity of accepted runs is the plan's strongest safety claim. V-4's exact-set assertion buys most of the protection for one line, and is the right trade here. Worth recording as the reason the duplication was accepted.

### 5.3 Task 2: five inline pushes vs. one `planned_clump_secondary_outputs` helper (worth considering)

Task 2 as written adds the same `if !cli.no_report_file { push(clumping_report_name(…)) }` at four sites, and Task 1 adds a fifth copy in a closure. The trim path already solved this shape once: `planned_secondary_outputs` (`main.rs:89-110`) centralises "the secondaries an SE trim run writes", and its doc comment even explains the A2/A3-style reasoning. A `fn planned_clump_report_outputs(cli, inputs, output_dir) -> Vec<PathBuf>` would give five call sites one line each, one gate, and — usefully — a unit-testable seam for A3-style properties without driving the binary.

Trade-off: the four arms key their reports differently (each input for SE; `inputs[0]` per chunk for the BAM paired shapes), so the helper needs either a small enum or two functions, which erodes some of the win. Given the maintainer confirmed Task 2 is in scope and the plan calls the additions "four trivially-reviewable lines", inline is defensible. But the repo's own history (#383 → #385 → #388 → #391, four rounds of "another dispatch path had no pre-flight") argues that the reviewable unit should be "one list every clump arm shares" rather than "four lists that currently agree". I'd raise it with the maintainer as a genuine fork rather than assume inline.

### 5.4 Rejected alternative worth naming explicitly: discriminate the report name instead of rejecting

Give clumping reports a positional discriminator (`reads_clumped_1.fq_clumping_report.txt`, or key the report off the *output* path rather than the input) and the collision disappears — both reports get written, no run fails. This is arguably what a user hitting the issue wants. Reasons to reject, which the plan should state since it is the first question a reviewer will ask:

- it changes output filenames, which is the one thing this repo treats as near-invariant (MultiQC/nf-core glob on `*_clumping_report.txt`; the naming is deliberately distinct from `*_trimming_report.*` per `io.rs:463-468`);
- it diverges from #388, which chose rejection for the structurally identical trim-report case three commits ago — two different answers to one defect class is worse than either answer;
- rejection is strictly safer: it never silently produces a file whose name the user did not predict.

The plan's "Rejected alternative" section only discusses a duplicate pre-flight inside the clump arm. Adding this one makes the design record complete.

### 5.5 Test placement (the plan's open question) — keep it in `integration_output_collision.rs`

Agreeing with the plan, with a reason: that file's module doc (`:1-18`) is written as the narrative of the pre-flight defect family, and its three shared constants (`DUP_MSG`, `ALIAS_MSG`, `PREFIX`) plus `assert_rejected_cleanly` / `assert_dir_holds_only` are exactly what these tests need. Splitting collision cases by mode would mean re-deriving those helpers in `integration_clump_only.rs`, whose stated job is byte-identity and determinism. Close the question in the plan's favour.

---

## 6. Action items

### Critical

None. The plan is implementable as written and closes the issue; every code citation verifies; the `Vec<PathBuf>` widening is provably pre-flight-only; the gate matches the writer exactly.

### Important

1. **A-1 — Fix the hint text before adopting it.** `PAIRED_REPORT_HINT`'s "regardless of source directory" is false on this arm (dropping `--output_dir` does separate the reports), and "Paired outputs … named from the input filename alone" points the user at the primaries, which are fine. Either add a clump-specific constant or fix the wording at `PAIRED_REPORT_HINT` — noting the same clause is equally false at #388's site. Record the choice in the plan.
2. **V-1 — Make A3 machine-checked.** Extend `src/io.rs:1229`'s primary set with `clumped_output_name` and `clumped_bam_output_name`. A3 is what justifies *not* fixing report-vs-report on four arms **and** what makes Task 2 rejection-neutral (§1.4); it should not live only in a plan file. Also add `_clumped_N` / #391 to the doc comment at `io.rs:1287-1291`.
3. **V-2 — Respecify test 5.** It cannot use `assert_rejected_cleanly` (the directory holds the inputs). Use the no-`-o` two-inputs-in-one-directory shape with `assert_dir_holds_only`, and keep the aliasing input valid FASTQ so the rejection is attributable to the pre-flight.
4. **V-4 — Turn test 3 into an exact-set assertion.** `assert_dir_holds_only(&out, &[<4 filenames>])`. One line; it is the only test shape that catches the next namer-vs-writer divergence, which is the root cause of this whole family.
5. **V-3 — Resolve the Task 2 representative count.** Write the paired-BAM Shape A test (it is the one arm keying on `inputs[0]` rather than per-input) or strike the claim from step 6.
6. **A-4 — Correct A3's `--basename` sentence.** Multi-input `--basename` never reaches the pre-flight; it is rejected at `cli.rs:624` (SE) and `cli.rs:651` (paired). Conclusion unchanged, mechanism wrong.
7. **A-2 — Decide the hint assertion deliberately.** After A-1, assert only the durable fragment (`--no_report_file`) in test 1, so the hint cannot silently revert and the inaccurate half is not pinned.
8. **5.3 — Put the "four inline pushes vs. one shared helper" fork to the maintainer** rather than assuming inline. The trim path's `planned_secondary_outputs` is the established precedent, and this repo has now litigated missing-pre-flight four times.

### Optional

9. **V-5a** — test 1 needs two assertions (`DUP_MSG` via the helper, plus the report filename separately).
10. **V-5b** — five unique `tempdir` tags; note the vacuous-pass risk for test 2 if tags collide.
11. **V-6** — add (or explicitly mark as covered-by-construction) the paired-arm report-vs-input case.
12. **A-3** — record that the report-vs-report message names the same path twice and cannot identify which input to rename; file a follow-up for input provenance in `preflight_output_collisions`.
13. **A-5** — sharpen A4: on clump paths FastQC's outputs derive from the discriminated primaries, so the residual is narrower than the general case.
14. **§1.3** — record the report/primary directory asymmetry on the paired arms (both primaries → `input_r1.parent()`, R2's report → `input_r2.parent()`) as a follow-up issue; it is pre-existing and shared with #388.
15. **5.4** — add "discriminate the report name instead of rejecting" to the rejected-alternatives list with the MultiQC-glob and #388-consistency reasons.
16. **Step 1** — the `run_specialty_paired` doc comment's mode list (`main.rs:2470-2471`) still says "(`--clock`, `--implicon`)"; add `--clump_only --paired`.
17. **Step 4** — the SE FASTQ arm's candidate list is an iterator chain (`main.rs:553-557`); adding reports needs a `for` loop or `flat_map`, not a `push`. Spell out which, since the plan's own guidance warns about post-`cargo fmt` anchor drift during scripted patching.
18. **Validation 7 is already answered:** no test in the suite pins `GENERIC_ADVICE`; the only hint assertions are the two `"current working directory"` checks at `tests/integration_output_collision.rs:625,661`, which cover `CWD_OUTPUT_HINT` on hardtrim/clock and are unaffected. The plan can state this as verified rather than as a step to perform.
