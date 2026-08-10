# Plan Review B — #398 one output directory per pair

**Reviewer:** B (independent; no shared state with Reviewer A)
**Plan:** `plans/08092026_paired-output-dir-unification/PLAN.md` r1
**Tree state:** `dev` @ `38deb1ee` (verified `git log --oneline -1`)
**Method:** every line-number claim in the plan was opened and read; the current
behaviour was reproduced with `target/release/trim_galore` in fresh directories
under `$TMPDIR`; `.github/workflows/ci.yml` and the docs tree were grepped
directly. Findings are marked **CONFIRMED** (I ran or read the specific thing) or
**SUSPICION** (reasoned, not executed).

**Verdict: REVISE** — the plan is sound in its core mechanism and its arm
enumeration is complete, but it (a) mis-states the real-world blast radius of the
accepted refusal, (b) leaves `passthrough_output_name` as a divergent sixth rule
with four unit tests pinning behaviour production will no longer reach, (c) gives
a step-7 instruction for `clump_report_candidates` that, taken literally, breaks
the multi-input single-end callers, and (d) points its docs step at the one table
that does *not* state the directory rule while missing the two sentences that do
— one of which is already wrong today. None of these are re-litigations of the
three fixed decisions.

---

## 1. Logic review

### 1.1 Arm enumeration is complete — I enumerated independently

I did not take the plan's five-row table on trust. I enumerated every
`preflight_output_collisions` call site (`grep -n preflight_output_collisions
src/*.rs`) and mapped each to its writer. There are exactly eleven production
sites, which matches the "all 11 arms" figure in the open #414 task:

| Pre-flight site | Arm | Report candidates | #398 impact |
|---|---|---|---|
| `main.rs:464` | `--hardtrim5` | none | none (CWD rule) |
| `main.rs:495` | `--hardtrim3` | none | none (CWD rule) |
| `main.rs:658` | clump SE FASTQ | per input | none (SE) |
| `main.rs:703` | clump paired BAM, Shape B (1 interleaved input) | 1, `inputs[0]` | **no-op** |
| `main.rs:746` | clump paired BAM, Shape A (two-file, multi-pair) | 1 per pair, `chunk[0]` | **no-op** |
| `main.rs:803` | clump SE BAM | per input | none (SE) |
| `main.rs:915` | trim paired FASTQ | 2 per pair (`:906-912`) + carrier (`:893-903`) | **changes** |
| `main.rs:994` | trim SE FASTQ | via `planned_secondary_outputs` | none (SE) |
| `main.rs:2031` | trim paired uBAM-out (two FASTQ in) | 2 per pair (`:2026-2027`) | **changes** |
| `main.rs:2093` | trim SE uBAM-out | `:2089-2090` | none (SE) |
| `main.rs:2653` | `run_specialty_paired` — `--clock`, `--implicon`, clump paired FASTQ | clump: 2 per pair (`:616-620`); clock/implicon: none | **changes** (clump only) |

CONFIRMED. Nothing that writes a report or a per-pair output in paired mode is
missing from the plan's table. Specifically:

- **`--retain_unpaired`** — candidates `main.rs:878-888`, writer
  `main.rs:1506-1511` via `unpaired_output_names`, which is already R1-anchored
  (`io.rs:368-370`). Correctly out of scope. CONFIRMED.
- **`--demux`** — genuinely single-end-only, enforced at `cli.rs:1009-1012`
  ("Demultiplexing is only allowed for single-end files"). CONFIRMED; the plan's
  implicit assumption holds.
- **`--clock` / `--implicon`** — `grep -c report src/specialty.rs` → 0, and
  `clock_output_name` (`specialty.rs:488-491`) / `implicon_output_name`
  (`specialty.rs:521-524`) return a **bare filename** when `output_dir` is
  `None`, i.e. CWD. CONFIRMED. Note the plan's stated reason for excluding them
  ("write no reports") is true but incomplete — the load-bearing reason is that
  their primaries are CWD-anchored, so R1-anchoring would change them too. The
  conclusion is right; the justification is weaker than it needs to be.
- **`--clumpify`** (as opposed to `--clump_only`) writes no separate report; it
  routes through the ordinary trim report path (`main.rs:1536`, `:1832`), so it
  is covered by the trim arms. CONFIRMED by grep — `clumping_report_name` appears
  only in `clump_only.rs`.
- **`--clump_only` in all output-format combinations** — all four are in the
  table above; two are no-ops (see 1.2).

### 1.2 Two of the plan's five rows are no-ops, and the plan presents them as changes

`clump_only.rs:1080-1086`:

```rust
let report_input = &inputs[0];
let report_path = naming::clumping_report_name(report_input, output_dir);
```

`inputs[0]` **is R1** on Shape A and is the single interleaved input on Shape B.
So `clumping_report_name(inputs[0], None)` already resolves to R1's parent, which
is exactly what `pair_output_dir(input_r1, None)` returns. The candidate side
matches: `main.rs:740-744` passes `std::slice::from_ref(&chunk[0])`, and
`main.rs:698-702` passes the single input. CONFIRMED by reading both.

Consequence the plan does not state: for the "Clump paired → BAM" row, both the
writer change and the candidate change are **behaviour-identical no-ops**. That
matters two ways:

- An implementer who edits it and observes no test move may conclude the change
  "didn't take" and go looking for something to break.
- More usefully, it is a free assertion: the clump-paired-BAM report **already**
  obeys the #398 rule, and a test pinning that is a genuine regression guard for
  the R1-anchoring invariant on that arm.

Recommendation: mark the row `no-op by construction (report already keys on R1)`
and keep the edit for uniformity, or drop the edit and add the assertion.

**Related internal inconsistency:** the table has **five** rows, Behavior §2
refactors **five** `io.rs` namers (a different five), and Validation 8 says "all
**four** arms". Three different fives-and-fours in one plan is exactly the kind of
ambiguity that lets an arm be skipped. Validation 8 must name its arms.

### 1.3 The inverted guard: the plan's instruction is correct, but the obvious helper will fail

I read `tests/integration_output_collision.rs:1004-1042`. CONFIRMED:

- **Block 1** (`:1005-1014`): `--paired --basename foo r1.fq r2.fq`, both inputs
  in **one** directory. R1-anchoring is a no-op there — R1's parent already *is*
  both parents. Unaffected, passes untouched.
- **Block 2** (`:1016-1027`): `sample.fq` + `sample.fastq`, also **one**
  directory. Unaffected.
- **Block 3** (`:1029-1041`): `A/reads.fq` + `B/reads.fq`, no `-o`. Both reports
  resolve to `A/reads.fq_trimming_report.txt` after the change → refusal. Must be
  rewritten.

So the plan's "rewrite block 3, keep blocks 1–2" is **correct**. CONFIRMED by
reading, and by running the suite on `38deb1ee`:
`cargo test --release --test integration_output_collision` → **48 passed, 0
failed** (0.43 s), so the guard including block 3 is green today and block 3 is
the only part that must flip. That 48 is the number to re-check after the change:
the file should stay at 48 (block 3 rewritten in place) or rise to 49 (block 3
split into its own test, my preference — see below).

**Trap the plan does not flag.** The natural way to write the rewrite is
`assert_rejected_cleanly(&dir3.join("A"), ok, &stderr, DUP_MSG)`. That helper
(`:100-121`) asserts the directory is **empty** — and `A/` still holds the input
`reads.fq`, so it fails. An implementer hitting that failure is one step from
"fixing" the helper, which would weaken the emptiness assertion for the eight
other tests that rely on it.

The right helper already exists and is documented for exactly this case:
`assert_dir_holds_only` (`:74-92`) — *"For rejected runs (nothing written beside
the inputs — the no-`--output_dir` twin of 'directory is empty')"*. The plan's
Validation 3 asks for the correct property ("assert the whole directory listing
of `A/` and `B/` is unchanged") but does not name the helper. Say
`assert_dir_holds_only(&dir3.join("A"), &["reads.fq"])` explicitly.

**Second trap:** after the rewrite, the test function is named
`paired_report_candidates_do_not_over_reject` and doc-commented "T6 —
over-rejection guards at the boundary A1 describes", while containing a rejection
assertion. Either split block 3 into its own `#[test]` with a #398 name (my
preference — the two properties are independent and a single failure then tells
you which one broke), or rename and re-comment the function. The plan says
"rewrite … keeping its first two blocks", which as written leaves a misnamed
test.

### 1.4 The refusal class is broader than "two mates of one pair" — and the layout it breaks is a real one

This is my most important finding, and it is not a re-litigation: the decision to
accept a loud refusal is fixed, but I do not believe the maintainer has been shown
*which* invocations stop working.

The plan frames the refusal through `--paired A/reads.fq B/reads.fq` — which reads
as an artificial edge case. The actual rule after the change is: **any two inputs
in the run that share a filename and whose pairs resolve to the same directory
now collide.** The practically important instance is the *mate-side directory*
layout, where a facility splits reads by mate rather than by sample:

```
R1/s1.fq  R1/s2.fq
R2/s1.fq  R2/s2.fq
```

I ran this on `38deb1ee`:

```
$ trim_galore --paired R1/s1.fq R2/s1.fq R1/s2.fq R2/s2.fq     # exit 0
R1/  s1_val_1.fq s1_val_2.fq s1.fq_trimming_report.{txt,json}
     s2_val_1.fq s2_val_2.fq s2.fq_trimming_report.{txt,json}
R2/  s1.fq_trimming_report.{txt,json}  s2.fq_trimming_report.{txt,json}
```

CONFIRMED (exit 0, layout as shown). After the change the **entire invocation** is
refused before any pair is processed: `main.rs:863-914` builds candidates for
every chunk and then calls the pre-flight once (`:915`), so this is not "pair 1
fails and pair 2 proceeds" — nothing runs at all. That is a
layout convention, not a corner case, and it is qualitatively different from the
plan's §Integration framing:

> "pipelines that glob `<mate>_trimming_report.txt` beside R2's input will stop
> finding it … the realistic blast radius is hand-rolled scripts"

The blast radius is not only "scripts that look in the wrong place afterwards" —
it is "runs that previously completed now exit non-zero before writing
anything". Those are different severities and the changelog needs the second one
stated plainly, with `-o` (not `--no_report_file`) as the primary remediation,
since `--no_report_file` throws the reports away.

Conversely, the **per-sample directory** layout is safe, and saying so is worth as
much as the warning:

```
sampleA/reads_1.fq sampleA/reads_2.fq  sampleB/reads_1.fq sampleB/reads_2.fq
```

Each pair's R1 is in its own sample directory, so nothing collides. CONFIRMED by
running the analogous `--paired A/x1.fq A/x2.fq B/y1.fq B/y2.fq` — pair 1's
outputs landed in `A/`, pair 2's in `B/`, no cross-talk.

**Action:** add the mate-side-directory layout to the plan's edge-case table as
its own row, and put both layouts in the changelog and docs. Felix accepted "a
loud refusal for two same-named mates"; he should see that this reads, in the
field, as "the R1/-and-R2/ directory layout now requires `-o`".

### 1.5 The refusal message will name a path in R1's directory and will not explain why

Reproduced on `38deb1ee` with `-o out` (the case that already refuses):

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety): the outputs
named from A/reads.fq and B/reads.fq would be written to the same file,
out/reads.fq_trimming_report.txt. Outputs and reports are named from the input
filename alone, so inputs sharing a filename can collide when --output_dir (or a
shared input directory) sends them to one place — rename one input, or pass
--no_report_file if only the reports collide.
```

CONFIRMED (exit 1, `out/` empty, both inputs named — #397's provenance works, so
A3 holds).

After the change the same message fires with **no `-o` given** and the path
`A/reads.fq_trimming_report.txt`. The user sees two inputs in two different
directories, and a colliding path inside the *first* one. The current hint's
stated cause — "`--output_dir` (or a shared input directory) sends them to one
place" — is then flatly false: neither is present. The plan's step 8 correctly
identifies that the hint must be reworded but supplies no target text and does
not say what the new text has to *explain*. It must state the R1-anchoring rule,
or the message actively misleads:

> …every output of a paired run is written to Read 1's directory (or
> `--output_dir`), and reports are named from the input filename alone, so two
> mates sharing a filename resolve to one report path — pass `--output_dir`,
> rename one input, or `--no_report_file` if only the reports collide.

Note `PAIRED_REPORT_HINT` is shared by three sites (`main.rs:915`, `:2031`,
`:608`), so one rewording covers all of them. CONFIRMED by grep.

### 1.6 Step 6's stated justification rests on a caller that does not exist

Step 6: *"Prefer adding an explicit `output_dir` argument at the call site over
changing the namer, **so the single-pair SE-shaped uses stay put**."*

Two problems, both CONFIRMED:

1. `passthrough_output_name` **already** takes `output_dir: Option<&Path>` as its
   second parameter (`io.rs:387-392`). No argument needs adding; the change is
   passing `Some(&pair_dir)` instead of `output_dir` at `main.rs:1500` and
   `main.rs:895-900`. The step reads as if a signature change were required.
2. There are **no SE-shaped uses**. `--passthrough` requires `--paired` and
   exactly one pair (`cli.rs:906-916`: "`--passthrough` requires `--paired`",
   "requires exactly one R1/R2 pair"). `grep -n passthrough_output_name src/*.rs`
   returns exactly two production call sites, both on the paired path. The stated
   reason for choosing the call-site shape is false on the facts.

The decision may still be the right one, but it needs a real reason — and it has
a cost the plan does not account for (1.7).

### 1.7 The plan leaves a sixth, divergent directory rule in `io.rs` — with four tests pinning it

After the change, `io.rs` contains five paired namers calling `pair_output_dir`
and one — `passthrough_output_name` (`io.rs:404-409`) — still deriving the
directory from the **carrier input's** parent, in a branch no production path can
reach. Its docstring (`io.rs:375-386`) even calls the carrier one of "the three
pair outputs" while naming it as if it were not, which the plan itself quotes as
evidence of the bug.

Worse, four unit tests pin the unreachable branch: `test_passthrough_output_name_bare`
(`io.rs:854`), `_plain` (`:861`), `_with_basename` (`:875`), `_all_extensions`
(`:889`) — all pass `output_dir = None` and assert `/data/I1_passthrough.fq.gz`,
i.e. *beside the carrier*. They will keep passing after the change while
asserting a rule the binary no longer produces. CONFIRMED by reading all four.

That is the same shape of defect #398 exists to fix: a namer whose documented
rule differs from what its callers do. It also means the passthrough arm has the
**weakest safety net of the five** — I checked `tests/integration_passthrough.rs`
and its only positive test (`passthrough_cli_cores_1_smoke`) passes
`--output_dir`, so **no existing test pins the carrier's no-`-o` location in
either direction**. A half-implementation (candidate moved, writer not, or vice
versa) is invisible to the current suite.

Recommendation, in order of preference:

1. Change `passthrough_output_name` to take `input_r1` and derive the directory
   via `pair_output_dir`, updating the four unit tests to the new rule. It is
   paired-only, so this is safe, and it makes the drift structurally impossible.
2. If keeping the call-site shape: rewrite the docstring to say the directory is
   the caller's to supply, and repoint the four `None` tests at the pair rule.

Either way, add an integration test for `--passthrough` **without** `-o`
(Validation 2 asks for this — make sure it lands as a committed test, not just a
manual check).

### 1.8 Step 7's `clump_report_candidates` instruction breaks the SE callers if taken literally

Step 7: *"give it an explicit directory argument so the paired callers pass the
pair directory and the SE caller passes `output_dir`. **Do not** let it keep
deriving per input."*

`clump_report_candidates` (`main.rs:121-139`) is called from **five** sites, not
two:

| Site | `report_inputs` | Needs |
|---|---|---|
| `main.rs:616` | `&[r1, r2]` (one pair) | the pair directory |
| `main.rs:653` | `&cli.input` — **all** SE FASTQ inputs at once | per-input parent when `output_dir` is `None` |
| `main.rs:698` | `&cli.input` (len 1, Shape B) | that input's parent |
| `main.rs:740` | `from_ref(&chunk[0])` | the pair directory (= R1's parent, no-op) |
| `main.rs:798` | `&cli.input` — **all** SE BAM inputs at once | per-input parent |

CONFIRMED by reading all five. The two SE sites hand it the whole input list, and
those inputs may live in **different** directories. A signature that takes one
required directory and applies it to every element is wrong for them: it would
collapse every SE input's report into one directory, silently producing candidate
paths the SE writers never write — the #383 family in the mirror direction, and
this time on a *single-end* arm the plan declares out of scope.

The workable shape is: keep the parameter `Option<&Path>` with the existing
"`None` → per-input parent" semantics, and have the paired callers pass
`Some(&pair_dir)`. Then "do not let it keep deriving per input" is true of the
paired callers only, which is what the plan means but not what it says. Reword
step 7 so the SE contract is explicit, or the sentence is an invitation to break
two arms.

### 1.9 The lockstep argument cuts in the plan's favour on one point worth recording

Because `guarded_inputs` (`main.rs:72-81`) feeds the same pre-flight, moving a
report candidate automatically moves its **output-vs-input** check too. So the new
risk of a relocated R2 report landing on top of a *named input* in R1's directory
is covered for free — provided the candidate moves with the writer. That is a
concrete correctness reason for the plan's insistence on same-commit lockstep,
stronger than the "instance #6" framing, and worth a sentence in the plan.

### 1.10 The report text embeds the carrier's output path

`report.rs:466-481` writes `Input:` and `Output:` lines for the passthrough block
into the R2 **text** report. Moving the carrier therefore changes report content,
not just file location. The JSON report carries counts only
(`report.rs:1010-1035`), so no JSON value changes. CONFIRMED by reading both.
Harmless, but §Integration says "Writes: relocates … the passthrough carrier" and
should add "and the `Output:` line of the R2 text report's passthrough block", so
a reviewer diffing report fixtures is not surprised.

---

## 2. Assumptions

| # | Claim | Verdict |
|---|---|---|
| A1 | All five primary paired namers share an identical directory expression | **CONFIRMED.** Read all five: `io.rs:296-298`, `337-339`, `368-370`, `473-475`, `522-524`. Byte-identical `output_dir.map(\|d\| d.to_path_buf()).unwrap_or_else(\|\| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf())`. Extraction is pure. |
| A2 | The three report namers derive filename from `file_name()` and directory independently | **CONFIRMED.** `io.rs:535-548`, `551-564`, `567-580`. Passing a directory relocates without renaming. |
| A3 | Report paths are already in the candidate list, so the new collision is caught not overwritten | **CONFIRMED empirically**, as the plan demands. `--paired -o out A/reads.fq B/reads.fq` on `38deb1ee` → exit 1, message names both inputs and the colliding report path, `out/` left empty. |
| A4 | `--clock`/`--implicon` write no reports | **CONFIRMED.** `grep -c report src/specialty.rs` → 0. |
| A5 | Single-end paths never call the paired namers | **CONFIRMED by grep** — no `*_paired_*` namer is reachable from an SE dispatch branch. But see 2.1: one of the line references the plan offers as evidence is wrong. |
| A6 | `Path::parent()` on a bare filename yields `Some("")` and the existing `unwrap_or(".")` handles it | **Conclusion CONFIRMED, mechanism wrong.** See 2.2. |

### 2.1 A5's supporting citation is wrong

The out-of-scope list cites *"Single-end arms (`main.rs:100-101`, **`:909`-adjacent
SE paths**, `:1299`, `:2089-2090`, `:2118`)"*. `main.rs:909` is
`candidates.push((naming::report_name(input, output_dir), src.clone()));` inside
the **paired** candidate loop (`:906-912`) — it is the very line the change edits,
not an SE path. The SE candidate builder is `planned_secondary_outputs`
(`:91-115`, called at `:993`). CONFIRMED by reading. The other four citations are
correct. A5's conclusion holds, but a list presented as "out of scope by
inspection, not assumption" containing a mis-citation is worth correcting — it is
the kind of slip that makes a reviewer wonder what else was skimmed.

### 2.2 A6's mechanism is not the one that fires

`Path::new("reads.fq").parent()` returns `Some("")`, so `unwrap_or(Path::new("."))`
is **not** reached; the empty path is joined, yielding the bare relative filename.
`unwrap_or(".")` fires only for a path with no parent at all (`/`). A6 says "so
the existing `unwrap_or(".")` path is unchanged", which describes the wrong
branch. The **conclusion** — no panic, no path change — is CONFIRMED empirically:
`cd bare && trim_galore --paired r1.fq r2.fq` on `38deb1ee` wrote
`r1_val_1.fq`, `r2_val_2.fq`, `r1.fq_trimming_report.{txt,json}`,
`r2.fq_trimming_report.{txt,json}` all into the CWD, exit 0. Fix the wording so a
future reader does not go looking for a `.`-fallback that never executes.

### 2.3 Implicit assumptions the plan does not state

- **`pair_output_dir` gets no unit test of its own.** The plan asserts A6 in prose
  and relies on the five namers' existing tests for A1, but adds no test for the
  new function's own contract (`output_dir` precedence; bare-filename parent;
  root path). Since it becomes the single source of truth for every paired
  directory decision, it should carry three direct assertions. Not listed in the
  validation table.
- **The pre-flight does not plan FastQC output paths, deliberately.** `io.rs:1335-1340`
  states this explicitly: *"FastQC's `<stem>_fastqc.zip` is not asserted
  separately … Its one real exception — `--fastqc_args "-o DIR"` … is recorded in
  the plan as a known residual."* CONFIRMED. This directly bounds Validation 8
  (see 4.1).
- **`--basename` neither creates nor masks the new collision.** With
  `--basename foo` the primaries become `foo_val_{1,2}` (distinct by suffix) while
  reports keep per-input names, so two same-named mates in different directories
  still collide **on the reports** — `--basename` is not an escape hatch. And
  `--basename` is rejected for multi-pair at `cli.rs:669-671` ("--basename cannot
  be used with multiple paired-end pairs"), so it cannot create a cross-pair
  collision. The plan's edge-case row is correct. CONFIRMED by reading
  `paired_end_output_names` (`io.rs:326-341`), `cli.rs:669-671`, and the guard
  test's block 1, which exercises exactly this and passes today.
- **Anchoring is per pair, and the primaries already behave that way.**
  CONFIRMED by running `--paired A/x1.fq A/x2.fq B/y1.fq B/y2.fq` (1.4): pair 1's
  `x1_val_1.fq` + `x2_val_2.fq` both in `A/`, pair 2's both in `B/`. The reports
  will match. `run_specialty_paired` (`main.rs:2648-2653`) builds candidates
  per-chunk via the closure, so the clump arm inherits the same per-pair
  behaviour for free.

---

## 3. Efficiency

Nothing to report. The plan's "nil" is correct: all five namers already allocate a
`PathBuf` for the directory (`d.to_path_buf()` / `.to_path_buf()`), so returning
one from `pair_output_dir` is the same allocation moved, not an extra one. One
additional `PathBuf` per pair at the writer sites, on a path that then does file
I/O. No new syscalls, no change in complexity or memory. CONFIRMED by reading the
five expressions.

Minor: the writer sites could pass `pair_dir.as_path()` down rather than
recomputing, but `pair_output_dir` is two branches and a clone — not worth
plumbing.

---

## 4. Validation sufficiency

Validations 1–11 cover the main risk (writer/candidate drift) better than most
plans in this repo, and Validation 9's inverted control is the right instinct.
Three gaps.

### 4.1 Validation 8 will fail spuriously, or get weakened, on any `--fastqc` run

Validation 8: *"for each arm, run it and compare every file created against that
arm's candidate list → exact match; no writer produces an unplanned path."*

That property is **false by design** for two flags: FastQC's `*_fastqc.html` /
`*_fastqc.zip` are never planned (`io.rs:1335-1340`, CONFIRMED), and `--demux`'s
per-barcode files are planned only via `planned_secondary_outputs` on the SE arm.
An implementer running Validation 8 with `--fastqc` on will see unplanned files
and either report a false failure or — the dangerous outcome — "fix" it by
relaxing the exactness to a subset check, which is precisely the assertion
strength #391/#408 fought to get.

**Action:** scope Validation 8 to runs without `--fastqc`, and state the
known-unplanned set inline so nobody re-derives it under time pressure.

### 4.2 Validation 7 — I checked it, so the plan can stop asking

The plan asks whether CI's SE/PE invocations rely on per-mate report paths. Read
`.github/workflows/ci.yml` directly:

- Every Perl-comparison invocation passes `-o`: SE at `:383-384`, PE at
  `:394-395`, `-a2` matrix at `:411-414`, hardtrim5 at `:764-765`, clock at
  `:773-774`, demux at `:784-786`.
- Every PE invocation uses `test_files/BS-seq_10K_R{1,2}.fastq.gz` — one shared
  input directory, distinct filenames. So the change is a no-op there twice over.
- The md5 oracles compare **FASTQ payloads only** (`gzip -dc … | md5sum`, e.g.
  `:396-400`); no report is md5'd.
- The only report-path assertions are `:422` (`grep -q 'a AAATCAAAAAAAC'
  /tmp/rust_a2/BS-seq_10K_R2.fastq.gz_trimming_report.txt`), `:441-444`
  (multi-pair `test -f` on four report paths in `/tmp/rust_multi`) and `:843-845`
  (clumping-report filename discipline in `/tmp/co-det1`). **All three are inside
  `-o` directories.**
- The invocations that omit `-o` are all **negative** tests that exit before any
  write: `:466` (odd count), `:481` (paired + single FASTQ), `:490` (R1==R2),
  `:504` (basename + multi-pair), `:551`/`:562` (clock dup / odd count).

**CONFIRMED: the Perl validation matrix is entirely unaffected.** Validation 7 can
be restated as a one-line pin rather than an open question. Worth keeping as a
*test*, not a question: if someone later drops `-o` from the PE step to save a
`mkdir`, the matrix starts refusing.

### 4.3 Gaps to add

- **A `pair_output_dir` unit test** (2.3). Three assertions, five lines.
- **A committed `--passthrough` no-`-o` test** (1.7) — nothing pins the carrier's
  location today, so this behaviour change is currently unobservable to the suite
  in either direction.
- **The mate-side-directory layout as a validation row** (1.4) — the multi-pair
  refusal is the shape a user will actually hit, and Validation 6 only covers the
  *accepting* multi-pair case.
- **An assertion that the clump-paired-BAM report is already R1-anchored** (1.2) —
  converts a no-op edit into a guard.
- Validation 11's "count the delta against the 578 baseline" is the right
  instruction. Add: also check the *passed* count, not just `ok` — a filter that
  matches nothing prints `test result: ok. 0 passed`.

---

## 5. Alternatives for the implementation shape

The three behaviour decisions are fixed; these are about how to land them.

### 5.1 Paired report namers instead of "everyone calls `pair_output_dir`" (recommended consideration)

The plan's shape leaves two sites per arm each free to pass either `output_dir` or
`Some(&pair_dir)`. Drift needs only one forgotten call site, and the sites are
600+ lines apart (`main.rs:895` vs `:1500`; `:906-912` vs `:1729-1738`).

A strictly more drift-proof shape: add one paired-report namer to `io.rs`, e.g.

```rust
pub fn paired_report_names(input_r1: &Path, input_r2: &Path, output_dir: Option<&Path>)
    -> [(PathBuf, PathBuf); 2]   // (txt, json) for R1 and R2
```

and make **both** the candidate builder and the writer call it. Then it is
impossible for the two sides to disagree, because there is only one expression.
Same for the carrier: `paired_passthrough_name(input_r1, input_passthrough, …)`.
Cost: two more functions in `io.rs` and a slightly wider diff. Benefit: the
failure mode the plan calls "the single most likely way to get this wrong" becomes
unrepresentable rather than merely tested for.

Given the #383/#388/#391/#409 history — five recurrences of exactly this family —
I think the structural fix earns its diff. If the plan keeps the
`pair_output_dir`-everywhere shape, at minimum add a comment at each of the four
candidate sites pointing at its writer and vice versa.

### 5.2 Make the no-op explicit rather than editing it

For the clump-paired-BAM arm (1.2), prefer an assertion over an edit. Editing
identical behaviour into identical behaviour adds diff surface with no test
movement to prove it landed.

### 5.3 Do not fold the CWD rule in

The plan's open question recommends *not* backing `--hardtrim5/3`'s CWD rule with
`pair_output_dir`. I agree, and I would go further: the docs already state the CWD
rule as a deliberate, separately-documented behaviour in three places
(`modes/hardtrim.md:42`, `modes/clock.md:49`, `modes/implicon.md:32`, all
CONFIRMED, all identical wording from #387). Folding them together would make one
function answer two questions. Keep them apart.

### 5.4 The other open question (transition NOTE) — agree, no note

Agreed with the plan's "recommend no". A NOTE that fires on the now-correct common
case is noise. The changelog and the refusal message carry the transition.

---

## 6. Documentation findings (step 10 is pointed at the wrong lines)

The plan's docs path shorthand `guide/outputs.md` resolves to
`Docs/src/content/docs/guide/outputs.md` (note capital `Docs`). Its claim about
the Paired-end table at `:20-25` is CONFIRMED — and the table indeed says nothing
about directories. But the two sentences that *do* speak to this are not in the
plan:

- **`guide/outputs.md:102`** — *"`--output_dir DIR` writes outputs to `DIR/`
  instead of **the current working directory**. The trimmed FASTQ filename stem is
  unchanged; only the parent directory differs."* This is the one place the
  directory rule is stated, and it is **already wrong today**: trim outputs go to
  the *input's* parent, not the CWD (only the specialty modes use CWD). This is
  where the new rule belongs, and it needs correcting either way. CONFIRMED by
  reading, and by the `bare/` run in 2.2 plus the `A/`-anchored run in 1.4.
- **`guide/outputs.md:6`** — *"File naming matches v0.6.x exactly, so existing
  pipelines continue to work without changes."* Filenames still match after
  #398, so the first clause survives; the second becomes false for the
  mate-side-directory layout. Needs a qualifier.

Checked the two files the plan asks about, both CONFIRMED clean of any
directory claim, so nothing to correct there — only the new rule to add if
desired:

- `guide/paired-end.md:16-19` lists the four paired output filenames, no directory.
- `modes/clump-only.md:45`, `:109`, `:122` describe report *filenames* and
  per-mate vs per-pair counts ("Paired FASTQ runs write one report per mate; the
  uBAM shapes write one per pair (keyed on the pair's first input)" — which
  independently confirms 1.2), no directory.
- `guide/outputs.md:125` says FastQC artifacts land "alongside the trimmed
  output" — correct and unaffected, since FastQC is anchored on the output file
  (`fastqc.rs:37-58` passes `output_path` and the trim-galore `output_dir`
  straight through). **Answering the reviewer question directly: `--fastqc` needs
  no plan change.** Its artifacts already follow the primaries into R1's parent.
  The one residual is `--fastqc_args "-o DIR"`, already recorded as out of scope
  at `io.rs:1338-1340`.

**CHANGELOG:** step 11's `#### Changes` heading does not yet exist in the
`### Unreleased` section (currently only `#### Bug fixes`, `CHANGELOG.md:6`), but
the heading style is precedented (`:205`), as is `#### Behavioural notes (v2.x
intentional widenings)` (`:998`) — arguably the closer precedent for a deliberate
refusal of a previously-working invocation. Either works; just note the heading is
new to this section.

---

## 7. Action items

### Critical

1. **Reword step 7's `clump_report_candidates` instruction** so the SE contract
   survives (1.8). Taken literally, "do not let it keep deriving per input" breaks
   `main.rs:653` and `:798`, which pass **all** SE inputs at once and need
   per-input parents when `output_dir` is `None`. Specify: parameter stays
   `Option<&Path>`; paired callers pass `Some(&pair_dir)`; SE callers keep passing
   `output_dir`. CONFIRMED by reading all five call sites.

2. **Name the layout the refusal breaks, in the plan, the changelog and the docs**
   (1.4). `R1/s1.fq R2/s1.fq R1/s2.fq R2/s2.fq` runs today (CONFIRMED, exit 0) and
   will be refused for every pair. Add it as an edge-case row and a validation
   row; lead the changelog remediation with `--output_dir`, not
   `--no_report_file`. This is the consequence I judge the maintainer most likely
   not to have seen.

3. **Specify what the reworded `PAIRED_REPORT_HINT` must explain** (1.5). Without
   the R1-anchoring rule in the text, the message names a path inside R1's
   directory while the user passed no `--output_dir` and has inputs in two
   directories — the stated cause is absent in the failing case. One hint serves
   all three sites (`main.rs:608`, `:915`, `:2031`).

### Important

4. **Resolve `passthrough_output_name`'s divergence** (1.6, 1.7). Step 6's stated
   justification is false — there are no SE-shaped callers (`cli.rs:906-916`), and
   the `output_dir` parameter already exists. Either move its directory derivation
   to `pair_output_dir(input_r1, …)` or fix its docstring and repoint the four
   `output_dir = None` unit tests (`io.rs:854`, `:861`, `:875`, `:889`), which
   otherwise keep asserting a rule production can no longer produce.

5. **Bound Validation 8 to no-`--fastqc` runs and state the known-unplanned set**
   (4.1). `io.rs:1335-1340` documents that FastQC paths are deliberately absent
   from candidate lists; an exact-match assertion will otherwise fail spuriously
   and invite weakening.

6. **Name `assert_dir_holds_only` in Validation 3 / step 9, and split block 3 into
   its own test** (1.3). `assert_rejected_cleanly` asserts an *empty* directory and
   cannot express "unchanged beside the inputs"; reaching for it first produces a
   failure whose obvious fix damages eight other tests.

7. **Fix step 10's docs targets** (§6). Add `guide/outputs.md:102` (states the
   directory rule, and states it wrongly today) and `:6` ("existing pipelines
   continue to work without changes"). Correct the path prefix to
   `Docs/src/content/docs/`.

8. **Mark the clump-paired-BAM row a no-op and reconcile five-vs-four** (1.2).
   `clump_only.rs:1083`'s `report_input = &inputs[0]` is already R1, so both its
   writer and its candidate change are behaviour-identical. Validation 8 says
   "four arms" against a five-row table; name them.

### Optional

9. **Consider the paired-report-namer shape** (5.1) — makes writer/candidate drift
   unrepresentable rather than tested-for, which is worth weighing after five
   recurrences of this family.

10. **Add a `pair_output_dir` unit test** (2.3): `output_dir` precedence, bare
    filename, root path.

11. **Fix A5's `:909` mis-citation and A6's `unwrap_or(".")` mechanism** (2.1,
    2.2). Both conclusions are correct; both citations point at the wrong thing.
    `io.rs:286-305` for `paired_bam_output_name` is also off — the function is
    `285-301` and its directory expression is `296-298`.

12. **Add the report-content knock-on to §Integration** (1.10): the R2 text
    report's passthrough `Output:` line moves with the carrier
    (`report.rs:481`). JSON is unaffected.

13. **Record the free win from lockstep** (1.9): moving a report candidate also
    moves its output-vs-input check, so a relocated report cannot silently land on
    a named input.

---

## 8. Verdict

**REVISE.**

The mechanism is right, the arm enumeration is complete (I checked all eleven
pre-flight sites independently), A1–A4 and A6's conclusion all verify, and the
inverted-guard instruction is correct. Validation 7's open question is now
answered: the Perl matrix is unaffected, twice over.

What needs to change before implementation: item 1 is a latent break of two
single-end arms the plan declares out of scope; item 2 is a user-facing
consequence materially larger than the plan's framing; item 3 leaves the refusal
message stating a cause that is absent in the failing case; item 4 leaves `io.rs`
with the same documented-rule-vs-caller-rule divergence that #398 exists to
remove. Items 5–8 are the difference between a validation suite that catches drift
and one that gets relaxed when it fires.
