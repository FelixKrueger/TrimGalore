# PLAN_REVIEW_A — #423 → #424 → #429, one PR

**Reviewer:** A (independent; Reviewer B reviewing the same plan with no shared state)
**Plan:** `plans/08142026_paired-report-namers-batch/PLAN.md`
**Tree:** `dev` @ `d2ea293`, read-only. Binary at `target/release/trim_galore` (no build run).
**Verdict:** **Not safe to implement exactly as written** — three gaps, all cheap to close in the
plan. With C1–C3 addressed the batch is low-risk and the order is sound.

---

## 0. What I verified empirically (so the rest can be read as findings, not opinion)

| Claim under test | Result |
|---|---|
| Candidate order at both candidate sites is `[R1.txt, R1.json, R2.txt, R2.json]` (§2.3) | **True.** Dumped via `TRIM_GALORE_DUMP_PLANNED_OUTPUTS=1`; exact order reproduced below |
| Neither interleaved arm reads `cli.basename` (§2.4) | **True.** `run_paired_ubam_single_file` 1807-1988 and `run_ubam_output_paired_single_file` 2411-2529 contain no `basename` |
| §2.6's write sets | **True and complete.** Verified for FASTQ / FASTQ+`--retain_unpaired` / uBAM / `--cores 4` / `--no_report_file` / uBAM+`--fastqc` |
| §3's adversarial names are harmless today (V4's premise) | **True.** `x_val_1.fq` (BAM content) → `x_val_1_val_1.fq`; `y_R1_trimming_report.txt` → `y_R1_trimming_report_R1_trimming_report.txt`; both inputs survived, exit 0 |
| #429's premise at `d2ea293` | **True.** `@` in a read name → `Error: failed to write BAM record / Caused by: invalid input parameter` |
| §2.1 — #428 removed the residue half | **True.** After the refused run the directory holds only the input; no `.bam`, no `.partial` |
| Line numbers in §2.2 for the four in-scope sites and the two hand-rolled blocks | **All correct at `d2ea293`** |
| `integration_output_collision.rs` pins candidate order | **Partly** — see F2 |
| The tripwire compares sets, not sequences | **True.** `planned_keys: BTreeSet<String>`, `tests/integration_preflight_tripwire.rs:417` |

Dump output, for the record (this is the V1 baseline for the two candidate sites):

```
src/main.rs:928   a/r1_val_1.fq                      src/main.rs:2045  b/r1_val.bam
src/main.rs:928   a/r2_val_2.fq                      src/main.rs:2045  b/r1.fq_trimming_report.txt
src/main.rs:928   a/r1_unpaired_1.fq                 src/main.rs:2045  b/r1.fq_trimming_report.json
src/main.rs:928   a/r2_unpaired_2.fq                 src/main.rs:2045  b/r2.fq_trimming_report.txt
src/main.rs:928   a/r1.fq_trimming_report.txt        src/main.rs:2045  b/r2.fq_trimming_report.json
src/main.rs:928   a/r1.fq_trimming_report.json
src/main.rs:928   a/r2.fq_trimming_report.txt
src/main.rs:928   a/r2.fq_trimming_report.json
```

Observed write sets for the two interleaved arms (input `ip.bam`, CWD output):

```
--paired ip.bam                    → ip_val_1.fq ip_val_2.fq
                                     ip_R{1,2}_trimming_report.{txt,json}
  + --retain_unpaired              → ip_unpaired_1.fq ip_unpaired_2.fq
  + --cores 4                      → identical set
  + --no_report_file               → primaries only
--paired --output-format ubam      → ip_val.bam  ip_R{1,2}_trimming_report.{txt,json}
  + --fastqc                       → ip_val_fastqc.{html,zip}   (tripwire allowance covers these)
```

`--gzip` does not change the extensions here: `gzip` is `!dont_gzip && is_gzipped(input[0])`
(`main.rs:416`), and a `.bam` input is never `is_gzipped`. So the arms' `ext` is `.fq` in practice,
but the candidate list must still take it from the same variable the writer uses (see F4).

---

## 1. Logic review

### C1 (Critical) — V1 is blind to the half of #423 that can actually change a path

V1 says the whole safety case for the refactor is a before/after diff of `DUMP_PLANNED_ENV` output.
That machinery exists only where a pre-flight runs. The two arms whose report paths **step 3 and
step 4 rewrite** have no pre-flight today — that is precisely why they are `Arm::Unplanned`, and
the tripwire asserts *nothing is dumped* for them. So V1 covers the two candidate sites (`:922`,
`:2040`) and cannot see `interleaved_report_names` at all.

Worse, the usual net does not close afterwards either:

- `grep -rn "_R1_trimming_report" tests/ docs/` → **zero hits.** Nothing in the repo pins those
  four filenames. `integration_ubam.rs:66` pins the interleaved *primaries*
  (`ubam_paired_test_val_{1,2}.fq` against committed references) and says nothing about reports.
- The promoted tripwire cases will not catch it. T1 is `created ⊆ planned`, and after #424 both
  sides come from the *same* namer. A joint rename — namer and writer both emitting
  `ip.bam_R1_trimming_report.txt` — keeps the subset relation true and the case green while every
  user's report filename has changed.

And the joint rename is the *likely* mistake, not an exotic one: `io.rs`'s existing report namers
key on `file_name()` (`report_name`, `io.rs:586-591` → `sample.fq.gz_trimming_report.txt`), while
these arms key on `file_stem()` (`main.rs:1817`, `:2422` → `ip_R1_trimming_report.txt`). An
implementer putting `interleaved_report_names` in `io.rs` next to `report_name` and composing it
from `report_name` produces `ip.bam_R1_trimming_report.txt` — a user-visible rename that breaks
nf-core/MultiQC glob expectations, with every proposed validation still green.

**Fix (cheap):** (a) capture the created-file baseline for both arms as step 0 — the two sets above;
(b) add `io.rs` unit tests for `interleaved_report_names` asserting literal paths, including the
`file_stem` convention, `--output_dir` override, a bare filename, and a two-dot name (`ip.bam.gz`);
(c) add one filename assertion at integration level (the promoted tripwire cases return
`Outcome.created` — assert the four report names appear, or use `assert_dir_holds_only` in
`integration_ubam.rs`). Restate V1 as covering the candidate sites only, and name (b)+(c) as the
safety case for steps 3-4.

### C2 (Critical) — V3 cannot be constructed; the plan contradicts itself

§2.4 (and #424's body) prove that on these arms `stem + suffix` can never equal
`stem + "." + ext`, that with one input there is no output-vs-output pair, and that
`guarded_inputs` is a singleton (`--passthrough` needs two inputs, `--demux` is SE-only —
`main.rs:74-83`). I re-derived it: a collision needs `file_name == file_stem + suffix`, i.e. the
stripped extension `".ext"` would have to equal a suffix beginning with `_`. Impossible.

V3 nonetheless asks for "an input named so an output would overwrite it … refused before any I/O".
**No such input exists.** As written V3 will be attempted, fail to reproduce, and either burn time
or get quietly downgraded — the worst outcome for a validation table. This also answers **O3**: the
refusal is unreachable, so the message question is academic for these arms. For the record, the
output-vs-input text does read acceptably for one input (`io.rs:193-202`): it names the source and
the colliding path and never names a second input. "one of its inputs" and "drop it from the input
list" are mildly odd with a singleton, but coherent, and `--output_dir` remains valid advice. **No
message variant needed — O3 closes as "no scope addition".**

**Fix:** replace V3 with (a) an assertion on the *content* of the new candidate lists — the promoted
tripwire cases plus a `DUMP_PLANNED_ENV` check that all 6 (uBAM) / 6-8 (FASTQ) expected paths are
dumped, and (b) an `io.rs` unit test driving `preflight_output_collisions` with a synthetic
one-input candidate list to pin the one-input message. Keep V4 as-is; it is real and it passed
today.

### C3 (Critical) — Step 6 is three edits, not one, and one of them is a design decision

Promoting the three cases *is* a flag flip (`Tripwire::new` defaults to `Arm::Planned`; drop
`.unplanned()` at `:1021`, `:1035`, `:1046`). But those are the **only** three call sites of
`unplanned()`, and no `comparison_*` self-test passes `Arm::Unplanned` to `compare`. After the
promotion:

1. `fn unplanned` is never called and `Arm::Unplanned` is never constructed → `dead_code` →
   **CI is `-D warnings` on `--all-targets`, so this fails.** Decision needed: delete `Arm`, the
   `arm` field, `unplanned()` and the `arm == Arm::Unplanned` branch (`:383-393`), threading the
   parameter out of `compare` and its 8 `Arm::Planned` call sites — or keep it with an `#[allow]`.
   I'd delete it: after #424 there are no unplanned arms, and a variant kept alive by an allow is
   how this file rots.
2. `the_site_inventory_is_complete` (`:1150-1163`) asserts **exactly 11** `preflight_output_collisions(`
   occurrences and its own failure message demands you add a case and a shape *before* touching the
   count. #424 makes it 13.
3. `every_site_is_reached_by_some_case` (`:1252-1334`) is `[(&str, &[&str]); 11]` plus
   `sites.len() == 11`. Two shapes must be added — `["--paired", "ip.bam"]` and
   `["--paired", "--output-format", "ubam", "ip.bam"]`; `ubam_paired_test.bam` is already staged
   into every case dir, so this is mechanical. Both new sites are distinct (each arm returns before
   the next pre-flight), so T5 is satisfied.

Also inside the same file: the module doc's Limits paragraph (`:63-65`) reasons *by analogy* to
"the two `Arm::Unplanned` cases", the `Arm::Unplanned` variant doc (`:214-221`) states the
`file_stem`/`basename` argument, the section header at `:1012` reads "the two arms with no
pre-flight", and `comparison_detects_a_modified_input`'s doc (`:1116-1118`) calls T3 "the
substantive assertion on the three `Arm::Unplanned` cases". Step 6 says "update the module doc";
it is four places plus a deleted-machinery decision.

### F1 — §2.2's inventory is wrong in two details (conclusion unaffected)

- "At `d2ea293` there are **six** `pair_output_dir` call sites" — there are **eight** outside
  `io.rs`: the six in `main.rs` plus `clump_only.rs:537` and `:1089`. Both new ones are
  `clumping_report_name` writers, i.e. the same "different scheme" bucket as `main.rs:619`/`:746`,
  so nothing enters scope. But the plan's inventory reads as exhaustive and is not, and the same
  sentence is what a reviewer would check against later.
- The two clump rows describe each other's behaviour. `:619` passes `&[r1, r2]` → **two** per-mate
  clumping reports (matching `clump_only_paired`, `clump_only.rs:535-551`). `:746` passes
  `from_ref(&chunk[0])` → **one** report keyed on the pair's first input (matching
  `clump_only_paired_to_bam_one_pair`, `clump_only.rs:1085-1092`). The plan attributes the
  one-report shape to `:619` and calls `:746` "same shape". They are not the same shape — which
  matters if anyone later tries to fold them into one helper.

No seventh in-scope site: `write_paired_reports` (`main.rs:2557`) constructs no paths, it takes
`PairedReportFile` descriptors; SE report paths call `report_name(input, output_dir)` directly and
are out of scope.

### F2 — the order guard is narrower than §2.3 implies

`integration_output_collision.rs:1040-1043` asserts `stderr.contains("reads.fq_trimming_report.txt")`
for a co-located-report refusal. That pins **txt-before-json within a mate**: swap them and the
collision fires on `.json`, and this test goes red. Good — that is the shape §2.3 warns about.

What it does **not** pin: mate order (R1 before R2), or report candidates coming after the
primaries. The 2a message renders `prev_source` then `source`, and no test asserts that order
(`:902` checks a substring, not a sequence). So V2 catches the specific reorder the plan fears and
nothing else — which is fine given the per-mate signature, but it means **V1's order-sensitive diff
is the real guard, not V2.** Say so, or the safety case looks broader than it is.

### F3 (Important) — the report paths must be hoisted, which the plan does not mention

In both interleaved arms the report paths are computed *after* the whole trimming run
(`main.rs:1951-1968`, `:2481-2498`). A pre-flight "before any I/O" therefore cannot reuse them
where they stand: step 5 requires **hoisting** `interleaved_report_names(...)` above the reader/
writer setup (above `:1844` and `:2436` respectively) and threading the four paths down into the
`if !cli.no_report_file` block. That is a real structural edit, it interacts with
`--no_report_file` (compute unconditionally, plan conditionally — or compute inside an `Option`),
and it is not in the step list.

### F4 (Important) — build #424's candidates from the writers' own variables

The plan says "build from `interleaved_report_names` plus the primaries in §2.6". If the candidate
list re-derives `stem`/`ext`/`dir` instead of reusing `output_r1`, `output_r2`,
`unpaired_r1_path`, `unpaired_r2_path` and the hoisted report paths, #424 reintroduces exactly the
drift #423 exists to delete — one arm, two derivations. Recommend stating explicitly: the
pre-flight goes **inside** each arm function, immediately after the paths are computed and before
any writer is opened, and the candidate list is assembled *from those bindings*. That also settles
where it goes (main's dispatch at `:854`/`:2004` would have to recompute everything).

Side effect worth stating: inside-the-arm means the refusal lands after `setup_trimming`'s adapter
auto-detection has already printed. Harmless (detection only reads) and consistent with the arms'
existing shape.

### F5 — #429: the wrap site has the name; the plan's "record ordinal" does not exist

`BamWriter::write_record` (`bam.rs:593-644`) already holds `name: Vec<u8>` from
`parse_name_and_data`, and the failing call is
`write_alignment_record(...).context("failed to write BAM record")` at `:641-643`. So **O2's
preferred option exists** — no predicate duplication needed.

Three corrections to Step 7 / V6:

1. **There is no record ordinal available and none is asked for.** `BamWriter` has no counter
   (`bam.rs:525-538`); adding one is new state. The cited precedent, `reject_whitespace_in_qname`
   (`bam.rs:883-895`), names *the offending name* and no ordinal; #429's body asks for "the record
   and the offending byte". Drop "ordinal" from the plan, or state that a counter is being added.
2. **The >254-byte case is a different error.** Verified: a 255-byte name fails with
   `Caused by: out of range integral type conversion attempted`, not `invalid input parameter`. A
   pre-check on the character class alone would leave that case opaque — a second, stronger reason
   to wrap rather than duplicate the predicate. It also means the wrap message must not claim
   `[!-?A-~]` as *the* cause; phrase the constraint to cover both (1-254 printable ASCII excluding
   `@`).
3. **The wrap must not misdiagnose I/O failures.** `write_alignment_record` also fails on disk-full,
   EPIPE and friends; unconditionally appending "read names must be…" would relabel a full disk as
   a read-name problem. Since the call returns `io::Result`, gate the explanation on
   `e.kind() == std::io::ErrorKind::InvalidInput` before wrapping — the codebase already matches on
   `io::ErrorKind` in `cli.rs:513` and `format.rs:284`, so no downcast and no new idiom.

V7 confirmed already: all three name shapes (`@`, non-ASCII, 255 bytes) run clean to
`*_trimmed.fq` on the FASTQ-output path, so #415's accepted set is intact and the fix must stay
write-side.

### F6 — dead bindings after step 2

Removing the `report_name`/`json_report_name` pairs orphans the `pair_dir` bindings at
`main.rs:883`, `:1740`, `:2037` and `:2365` — at `:883` in particular the binding sits 39 lines
above its only two uses, so it is easy to leave behind. `unused_variables` is a hard error under
`-D warnings`. Trivial, but naming it stops it being "discovered".

---

## 2. Assumptions

| # | Plan's status | My finding |
|---|---|---|
| 1 | Fixed: clump sites out of scope | **Holds**, with F1's corrections to the count and the descriptions |
| 2 | Fixed: no `basename` in either arm | **Verified true** |
| 3 | Fixed: order load-bearing, tripwire compares sets | **True**, but the pin is narrower than implied (F2) |
| 4 | Fixed: #429 is a message fix | **True**; scope is one wrap site, two distinct underlying error kinds (F5) |
| 5 | Assumed: #423 is a pure refactor, V1 makes it checkable | **Half-true.** V1 checks the candidate half only (C1) |
| 6 | Assumed: one-input message reads sensibly, V3 confirms | **V3 cannot run** (C2). The message does read acceptably; verified by reading `io.rs:193-202` |

Two implicit assumptions the plan should surface:

- **A7 — `interleaved_report_names` keys on `file_stem()`, not `file_name()`.** Unstated, and it is
  the difference between a pure refactor and a silent rename (C1).
- **A8 — the promoted cases' write sets are complete for every flag combination they will ever
  run under.** True for the flags in §2.6 plus `--cores` and `--no_report_file` (verified). Note
  `--fastqc` is refused on the FASTQ interleaved arm (#421,
  `refuses_fastqc_on_every_arm_that_cannot_run_it:1377`) but accepted on the uBAM one, where it
  writes `ip_val_fastqc.{html,zip}` — allowed by the tripwire's FastQC rule, so a future
  `--fastqc` variant of that case stays green.

---

## 3. Efficiency

§6 is right and needs no change: two `parent()` lookups instead of one, twice per run. The
suggestion to pass a pre-computed `pair_dir` should be declined for the reason the plan gives —
taking `input_r1` is the entire point. #424 adds one hash of ≤8 paths per run.

One non-obvious cost saving: F4's "reuse the bindings" formulation makes #424 *smaller* than the
plan's version (no second derivation of stem/ext/dir), not larger.

---

## 4. Alternatives

- **#424 without #423 (answers "is the order forced?").** No — the dependency is stylistic. Hoist
  the existing `format!` block (F3) and build the candidate list from those same bindings (F4), and
  #424 hand-rolls nothing new: it reuses the one construction that already exists. So the order is
  *preferable* (namer first keeps `main.rs` shorter and gives the arms one place to look) but not
  *forced*, and the PR could be split if CI or review load ever argued for it. Keeping one PR with
  three commits (O1) is fine either way. Worth correcting the plan's "the order is not a preference"
  to "the order is a preference, and here is why it is the right one".
- **A single `PairedReportPaths` return instead of a per-mate tuple.** Correctly rejected in §2.3;
  I verified the order it would break and the test that would catch it (F2).
- **#429 by pre-check instead of wrap.** Available but strictly worse: it duplicates
  `[!-?A-~]{1,254}`, and the character-class half would still miss the length case (F5.2). Wrap
  with an `ErrorKind::InvalidInput` gate (F5.3).
- **Deleting `Arm` vs `#[allow(dead_code)]`.** Recommend deleting (C3.1) — but say which in the
  plan so it is not decided under CI pressure.

---

## 5. Action items

### Critical

1. **C1 — Rewrite V1 and add real coverage for `interleaved_report_names`.** State that V1 covers
   the two candidate sites only; add `io.rs` unit tests pinning the literal interleaved report
   paths (`file_stem` convention, `--output_dir`, bare filename, two-dot name) plus one integration
   filename assertion; capture the §0 created-file baselines as step 0. Without this the riskiest
   edit in the batch has no guard, and the most likely wrong turn — composing the namer from
   `report_name` — is invisible to every check in §7.
2. **C2 — Replace V3.** It cannot be constructed; §2.4 and #424's body both prove the refusal is
   unreachable. Substitute a candidate-list-content assertion (`DUMP_PLANNED_ENV` on both arms) and
   an `io.rs`-level one-input message test. **Close O3 as "no scope addition"** — the existing
   output-vs-input text reads acceptably for a singleton input.
3. **C3 — Expand step 6 to four items** and decide the `Arm` question: drop `.unplanned()` ×3;
   delete (recommended) or `#[allow]` the now-dead `Arm::Unplanned` / `unplanned()` /
   `compare`'s `arm` branch; bump `the_site_inventory_is_complete` 11 → 13; add two shapes and
   bump the count in `every_site_is_reached_by_some_case`. Plus the four doc spots at
   `:63-65`, `:214-221`, `:1012`, `:1116-1118`.

### Important

4. **F3 — Add the hoist to step 5.** Report paths currently materialise after the trimming run in
   both arms; the pre-flight cannot precede I/O without moving them above the writer setup and
   threading them down. Say how `--no_report_file` interacts.
5. **F4 — Specify that the candidate list is assembled from the arm's own path bindings**, with the
   pre-flight call inside each arm function. Re-deriving stem/ext/dir would reintroduce the drift
   #423 deletes.
6. **F5 — Fix Step 7 and V6.** Drop "record ordinal" (no counter exists; the precedent and the
   issue both ask for the name); note that the 255-byte trigger surfaces
   `out of range integral type conversion attempted`, so wrap rather than pre-check; gate the added
   explanation on `io::ErrorKind::InvalidInput` so a disk-full failure is not relabelled a
   read-name problem.
7. **F1 — Correct §2.2.** Eight `pair_output_dir` sites outside `io.rs` (add `clump_only.rs:537`,
   `:1089`), and swap the `:619`/`:746` descriptions: `:619` plans two per-mate reports, `:746`
   plans one keyed on R1.
8. **F2 — Downgrade V2's claim.** It pins txt-before-json within a mate and nothing else; V1's
   ordered diff is the actual order guard.
9. **F6 — Name the dead-binding cleanup** at `main.rs:883`, `:1740`, `:2037`, `:2365`.

### Optional

10. **V1's diff must ignore the site column.** The dump line is
    `PREFIX\t<file:line>\t<path>`, and `#[track_caller]` makes that line number shift with every
    edit above the call — a byte-for-byte diff is guaranteed to differ for a benign reason. Use
    `… 2>&1 | grep TRIM_GALORE_PLANNED | cut -f3` and diff that, order preserved.
11. **V6 fixture strategy.** Generate the offending FASTQ in-test (a `@bad@name` record) rather than
    committing another BAM; `tests/integration_ubam.rs:261+` (the #415 block) is the home, and its
    `stderr.contains(...)` style is the model.
12. **`write_paired_reports`' doc comment names two callers; there are four** (`main.rs:2533-2536`).
    One line, current state only, while step 2 is touching all four.
13. **Tripwire section headers name site line numbers** (`── site 1006 ──`,
    `── sites 2047 / 2109 ──`); both shift under this batch.
14. **Reword "the order is not a preference"** (§1, §10) — it is one, for the reasons in §4.

No new GitHub issues recommended: every finding above is inside this PR's scope, and the two
inventory inaccuracies (F1) have no user-visible symptom.

---

## 6. Ambiguities recorded, not blocking

- Whether `interleaved_report_names` should live in `io.rs` beside `report_name` (consistent
  location, invites the `file_name()` composition mistake) or be visually separated with the
  `file_stem` invariant stated on it. I'd keep it in `io.rs` — that is the file that owns naming —
  with the invariant on the function, per step 3's own instruction.
- Whether the promoted `orphan_paired_fastq_retain_unpaired` case should also assert the two
  `_unpaired_N` paths are planned-but-possibly-uncreated. T1 permits over-planning, so it is not
  required; asserting it would pin the conditional branch of the new candidate list.
- Whether O1's three commits survive to `master`. Squash-merge collapses them; the plan already
  says so.
