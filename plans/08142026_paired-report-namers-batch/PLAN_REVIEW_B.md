# PLAN_REVIEW_B — #423 → #424 → #429, one PR

**Reviewer:** B (independent; Reviewer A ran in parallel, no shared state)
**Plan:** `plans/08142026_paired-report-namers-batch/PLAN.md` v1
**Tree:** `dev` @ `d2ea293`, read-only. Binary at `target/release/trim_galore` (not rebuilt).
**Baselines measured here:** prebuilt `integration_preflight_tripwire` 55/55 pass,
`integration_output_collision` 57/57 pass.

**Verdict: not safe to implement as written.** Four defects, all small and local: one
that would ship a nonsense refusal message (C1 — the plan actively directs it), two that
make CI red for reasons unrelated to the intended signal (C2, C3), and one that leaves
four of the six refactored constructions pinned by nothing at all (C4). Everything else
in the plan checked out — including the parts I expected to find wrong.

---

## What I verified, and what it cost the plan

| Plan claim | Status |
|---|---|
| §2.1 #428 landed; #429 is message-only, no residue | **Confirmed empirically.** Failed uBAM-out run: exit 1, `find` shows only the input — no `.partial`, no reports |
| §2.2 line numbers `922-923`, `1742-1748`, `2040-2041`, `2367-2373`, `1959-1966`, `2489-2496` | **All six exact.** `pair_output_dir` for the two candidate sites is at `:883` / `:2037`; the pushes are at the quoted lines |
| §2.2 four in scope, two clump sites out | **Right, descriptions swapped** — see N1 |
| §2.2 no seventh site, nothing bypasses the namers | **Confirmed.** `grep -rn "trimming_report" src` gives io.rs's two namers plus exactly the two hand-rolled blocks |
| §2.3 candidate order is `[R1.txt, R1.json, R2.txt, R2.json]` | **Confirmed by running the dump** (output below) |
| §2.3 the tripwire compares sets | **Confirmed** — `planned_keys: BTreeSet<String>`, containment only |
| §2.3 `integration_output_collision.rs` pins order | **True but narrower** — see N2 |
| §2.4 neither interleaved arm reads `cli.basename`; `file_stem` + `_`-prefixed suffixes | **Confirmed**, and it survives the `.gz` case: `one.bam.gz` → `one.bam_val_1.fq.gz` |
| §2.5 the three cases assert exit 0 / created non-empty / inputs unchanged / nothing dumped | **Confirmed** in `Tripwire::run` + `compare`'s `Unplanned` branch |
| §2.6 write sets | **Confirmed exactly, empirically** — see N6 |
| §3 `@` refused on uBAM out, accepted on FASTQ out; 254 ok, 255 refused | **Confirmed**, but 255 fails with a *different* cause — see I4 |
| §3 #415 admits `@`, non-ASCII, VT, >254 read-side | **Confirmed** (`is_ascii_whitespace` excludes VT; bam.rs:881 says so deliberately) |
| §7 promotion is a flag flip | **Yes for `Arm::Planned`'s requirements, no for the file** — see C2, C3 |
| §1 "hard order because each depends on the last" | **Overstated** — see N3 |

---

## 1. Logic review

### C1 — Step 5's provenance choice would ship a nonsense message (Critical)

Step 5 says: *"`OutputSource::Pair` for paired primaries (#397 decision C2)"*. These arms
have **one** input, so `Pair` can only be `Pair(input, input)`. `OutputSource`'s `Display`
(io.rs:150-157) joins with `" + "`, so the refusal would read:

```
the output named from one.bam + one.bam would be written to one.bam, which is also one of its inputs.
```

I did not infer that rendering — I produced it. Forcing a `Pair`-provenance alias refusal
on the existing two-file arm (`--paired --basename foo foo_val_1.fq r2.fq`) gives:

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety): the output named from
foo_val_1.fq + r2.fq would be written to foo_val_1.fq, which is also one of its inputs. …
```

Substitute one input for both mates and the `X + X` is what a user gets.

`Pair`'s stated reason (io.rs:130-132, and the comment at main.rs:879-881) is that *neither
mate alone explains the path*. With one input, it explains the path completely. Use
`OutputSource::Input(input)` for every candidate on both arms — primaries and reports alike.

**This also closes O3, which the plan flagged Critical-if-surprising.** With `Input`
provenance a one-input refusal reads correctly. Measured, by making the SE uBAM-out arm
alias its own input (`--output-format ubam --basename x x_trimmed.bam`):

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety): the output named from
x_trimmed.bam would be written to x_trimmed.bam, which is also one of its inputs. If that file
is an earlier run's output, drop it from the input list; otherwise use `--output_dir` to write
elsewhere.
```

No second input is named, nothing reads oddly, and the `--output_dir` remedy is actionable
for a single input. **No message redesign is needed and O3 is not a scope addition** — but
only if the provenance is `Input`. Assumption 6 can be marked resolved before implementation
rather than confirmed by V3.

### C2 — Step 6 as written turns CI red on `-D warnings` (Critical)

`fn unplanned` (tripwire `:267-270`) is the **only** constructor of `Arm::Unplanned`
(`grep`: line 268 is the single occurrence outside `compare` and the docs). Delete the three
`.unplanned()` call sites and the method is dead → `dead_code` → `clippy --all-targets
-- -D warnings` fails.

Verified on a minimal reproduction of the exact shape (`rustc --edition 2024`):

```
warning: method `unplanned` is never used
```

and a second probe confirms the useful half: a variant that is only ever *compared* does
**not** warn. So the minimal correct edit is to delete `fn unplanned` too, keep the
`Arm::Unplanned` variant and `compare`'s early-return branch (they are the trap for a future
no-pre-flight arm), and add one self-test that drives `compare(…, Arm::Unplanned)` with a
non-empty `planned` and expects the "promote the case" error — otherwise that branch is
reachable by no test at all.

### C3 — two required test edits the plan does not mention (Critical)

`the_site_inventory_is_complete` (`:1149-1164`) asserts `preflight_output_collisions(`
appears exactly **11** times in `main.rs`. I counted: 11 today, 13 after #424. Its own
failure message says *"If one was added, give it a case and add its shape to
`every_site_is_reached_by_some_case`, then update this count — do not just update the count."*

`every_site_is_reached_by_some_case` (`:1253-1334`) is a `[(&str, &[&str]); 11]` array plus
`assert_eq!(sites.len(), 11)`. It needs two new shapes — `["--paired", "ip.bam"]` and
`["--paired", "--output-format", "ubam", "ip.bam"]` — and both counts moved to 13. The
harness already stages `ip.bam` for every shape (`:1311`), so no new staging.

Consequence for V5: the checkpoint will show **five** red tests, not three. Say so in the
plan, or the two extra reds get misdiagnosed as the promote signal misfiring.

### C4 — V1 is blind to four of the six replaced constructions, and nothing else covers them (Critical)

The dump only exists inside `preflight_output_collisions`, so it covers the two **candidate**
sites. It says nothing about:

- the two **writer** sites (`:1742-1748`, `:2367-2373`) — no dump there;
- both **interleaved** hand-rolled blocks (`:1959-1966`, `:2489-2496`) — those arms are
  `Arm::Unplanned` today and dump *nothing*, so a before/after dump diff is empty-vs-something
  for exactly the two sites Step 3+4 rewrite.

The writers are partly covered in practice: the tripwire's T1 (`created ⊆ planned`) goes red
if a writer's path moves while its candidate does not. The interleaved blocks have **no
coverage whatsoever**:

```
$ grep -rln "_R1_trimming_report" tests/    →  (no matches)
```

No test anywhere pins `<stem>_R{1,2}_trimming_report.{txt,json}`. And after Step 5 the
candidate list and the writer both derive from `interleaved_report_names`, so a typo
(`_trimming_reports.txt`, `_R1_` on both mates) is *self-consistent*: T1 stays green,
`unchanged` stays green, and the run silently renames a user-visible file. That is the
plan's own §10 failure mode in the place it did not look.

Two additions close it:

1. **io.rs unit tests pinning both namers' exact strings** — the file already has that
   convention (`test_pair_output_dir_precedence_and_edges` at `:1083`, report-name assertions
   at `:942-957`). The plan's §7 has no row for the namers' own tests.
2. **Extend V1 to the created-file set** for the two interleaved arms, not just the dumps.
   Baseline captured here at `d2ea293`, so it can be diffed directly after the refactor:

| Arm | Files written (verified) |
|---|---|
| `--paired one.bam` | `one_val_{1,2}.fq`, `one_R{1,2}_trimming_report.{txt,json}` |
| `+ --retain_unpaired` | above `+ one_unpaired_{1,2}.fq` |
| `--paired --output-format ubam one.bam` | `one_val.bam`, `one_R{1,2}_trimming_report.{txt,json}` |
| `--paired --no_report_file one.bam` | `one_val_{1,2}.fq` only |
| `--paired one.bam.gz` | `one.bam_val_{1,2}.fq.gz`, `one.bam_R{1,2}_trimming_report.{txt,json}` |

### I1 — V1's byte-for-byte diff will show every line changed (Important)

The dump line's second field is the `#[track_caller]` site. Real output from this tree:

```
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R1_val_1.fq
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R2_val_2.fq
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R1_unpaired_1.fq
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R2_unpaired_2.fq
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R1.fastq_trimming_report.txt
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R1.fastq_trimming_report.json
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R2.fastq_trimming_report.txt
TRIM_GALORE_PLANNED	src/main.rs:928	r1/s_R2.fastq_trimming_report.json
```

(Also the direct confirmation of §2.3's order claim.) `src/main.rs:928` shifts as soon as
main.rs changes length, so a byte diff is all-noise — and noise is where a real one-path
change hides behind "only line numbers moved". V1 must compare **field 3 only, in order**
(`cut -f3`). Env var is `TRIM_GALORE_DUMP_PLANNED_OUTPUTS`, not the shorter name.

### I2 — #424 re-creates the duplication #423 exists to delete, in the same PR (Important)

Step 5 builds the candidate list from `interleaved_report_names` **plus the primaries in
§2.6**. Those primaries are hand-rolled inside the arms — `format!("{}_val_1{}", stem, ext)`
at `:1826-1827`, `_unpaired_{1,2}` at `:1835-1836`, `_val.bam` at `:2430`. Deriving them a
second time at the dispatch site leaves two copies of every primary string in one PR whose
first commit deletes two copies of the report strings. That is #383/#388/#391/#409/#398's
shape exactly: a writer and a candidate list that agree only by hand.

Give the primaries namers too — `interleaved_paired_output_names(input, output_dir, gzip)`,
the unpaired pair, and the `_val.bam` name — or (better, see A1) one struct built once and
consumed by both sides.

### I3 — pre-flight placement is unspecified, and the obvious choice is the slow one (Important)

Every other arm pre-flights in `main()` before any per-input work; the two interleaved arms
are entered *after* `setup_trimming` has already run at the dispatch site (`:855`, `:2005`),
and `setup_trimming` is adapter auto-detection — a scan of up to 1 M reads. Put the
pre-flight at the dispatch site **before** `setup_trimming`, not at the top of the arm
function, or a refusal arrives only after the user has waited through a full detection pass.
`cli`, `output_dir` and `gzip` are all in scope there.

### I4 — V6's expected outcome contradicts Step 7's chosen approach (Important)

V6 expects *"no `invalid input parameter`"*. anyhow prints the whole cause chain, so a
`with_context` wrap yields:

```
Error: failed to write BAM record <name>
Caused by:
    invalid input parameter
```

The string V6 forbids is still there. Either restate V6 as "the message names the offending
read name" (correct, and what #429 asks for) or accept that only a pre-check satisfies it as
written — and the plan rightly prefers not to pre-check.

Relevant discovery while reproducing: the 255-byte name fails with a **different** cause —

```
Caused by:
    out of range integral type conversion attempted
```

— so this one wrap site emits two distinct opaque errors, and a pre-check duplicating
`[!-?A-~]{1,254}` would have to reproduce the length rule as well as the charset. That
strengthens the wrap decision; O2 can be closed in the plan's favour.

### I5 — #429's ordinal does not exist, and the name is moved before the wrap (Important)

§3 and V6 say the error names *"the record ordinal and the offending name, in the shape
`reject_whitespace_in_qname` already uses"*. Two corrections:

- That shape (bam.rs:883-897) has **no ordinal** — it names the escaped name, the reason and
  a remedy.
- `BamWriter` (bam.rs:524-537) has **no record counter**. An ordinal means a new field
  incremented per write. Decide it in the plan: add the counter, or drop the ordinal.

The name *is* reachable at the wrap site — the plan's O2 worry about an opaque error is
unfounded, `name` is a local from `parse_name_and_data`. But it is `Vec<u8>` and is **moved**
into `set_name(name)` on the line above, so the closure must capture something else. Use
`record.id`, matching the precedent five lines up: `format!("while writing BAM record {:?}",
record.id)` (bam.rs:613). One sentence of *why* in the message (this name is legal in a
FASTQ header but not as a SAM QNAME) costs nothing and does not duplicate the predicate in
code, so it cannot drift.

### I6 — add a gzip-input tripwire case (Important)

Verified: an input named `one.bam.gz` (content still BGZF/BAM, so format detection is
unaffected) produces `one.bam_val_{1,2}.fq.gz` + `one.bam_R{1,2}_trimming_report.*`. `gzip`
comes from `!cli.dont_gzip && is_gzipped(&cli.input[0])` (main.rs:416) and `is_gzipped` is
extension-only. So a candidate list that hardcodes `.fq` passes all three promoted cases and
silently mis-plans every `.gz` run — the pre-flight guards a path nobody writes while the
real output is unguarded. One extra case (`--paired ip.bam.gz`) pins it.

---

## 2. Assumptions

Assumptions 1-4 are correctly marked Fixed and all four hold, with the wording corrections
in N1 and N2. Assumption 5 (#423 is a pure refactor) is *not* made checkable by V1 alone —
that is C4; with the two additions it becomes checkable. Assumption 6 is **resolved now**,
not pending V3 — see C1.

Unstated assumptions I had to test:

- **That a collision is reachable on these arms at all.** §2.4 argues no output can equal the
  input (`file_stem` strips one `.`-suffix; every appended suffix starts with `_`), and I
  confirmed it holds for the `.gz` case too. Output-vs-output is unreachable (`_val_1` vs
  `_val_2`, `_R1` vs `_R2`, all distinct under case-folding). So **V3 as written — "an input
  named so an output would overwrite it, on both arms" — cannot be constructed.** That is not
  a flaw in #424 (defence-in-depth is the stated point) but V3 needs restating: either drive
  the arms' candidate list directly from a unit test, or verify the machinery on a sibling
  arm where a collision *is* reachable (the `--basename` trick above), and note in the plan
  that the interleaved refusal is unreachable by construction until something honours
  `--basename` there.
- **That `Arm::Planned` needs no machinery the promoted cases lack.** It does not: planned
  non-empty ✓, single site ✓, all candidates inside the case root ✓ (a bare-filename input
  gives `parent() == ""`, so paths stay bare and `root.join(p)` resolves), `created ⊆ planned`
  ✓ given §2.6. The blockers are C2/C3, which are about the file, not the enum.
- **That `--fastqc` does not enlarge the write set into unplanned territory.** On the uBAM
  interleaved arm it adds `one_val_fastqc.{html,zip}`, allowed by `is_fastqc_artifact`; on the
  FASTQ interleaved arm `Cli::validate` refuses the flag (#421). Neither needs candidates.
- **That `--clumpify` is not a variant to plan for**: refused earlier with *"--clumpify
  requires --cores >= 2"*. `--cores 2` gives the same write set as serial.

---

## 3. Efficiency

§6 is right and the arithmetic is trivial: two `parent()` lookups per pair instead of one, on
a path taken twice per run. Do not take the "pass the already-computed `pair_dir`" escape
hatch — it re-opens the derivation the namer exists to own.

The one real latency point is I3: pre-flighting *after* `setup_trimming` makes a refused run
pay for a full adapter-detection scan first. That is user-visible, unlike the `parent()` calls.

---

## 4. Alternatives

**A1 — one struct per arm instead of one namer per report pair (recommended if the durable
fix is wanted).** Build `InterleavedOutputs { val_1, val_2, unpaired_1, unpaired_2, r1_txt,
r1_json, r2_txt, r2_json }` once in `io.rs`; the candidate list is `struct → Vec<PlannedOutput>`
and the writers read the same fields. This subsumes I2, removes C4's self-consistency hole
(one derivation, unit-tested once, cannot disagree with itself because there is only one),
and is what "drift becomes unrepresentable" would actually mean. Cost: one type. This is the
shape I would pick.

**A2 — the per-mate `paired_report_names` signature is right, for a better reason than the
plan gives.** Order is a weak argument (N2). The strong one: the two *writer* sites build a
`PairedReportFile` **per mate**, so a per-mate namer is the only shape that serves all four
sites with one function. Keep the signature; upgrade the justification.

**A3 — splitting the PR.** #429 touches only `bam.rs` and depends on nothing here; it could
ship in an hour rather than waiting out two namer/pre-flight items in review. Felix chose one
PR, so this is a note. O1's three-commit structure is the right call and matches the real
dependency graph (N3).

**A4 — a pre-check for #429.** Correctly rejected, and I4 makes the case stronger: two
distinct noodles failure modes, so a pre-check duplicates a charset *and* a length rule.

---

## 5. Notes (real, low or no user impact — recorded, not for filing)

Per the issue-filing bar, none of these is worth an issue; N5 and N7 are worth one line each
in the plan.

- **N1 — §2.2's clump descriptions are swapped.** `main.rs:619` passes `&[r1, r2]` to
  `clump_report_candidates` → **two** clumping reports per pair, one per mate. `main.rs:746`
  passes `std::slice::from_ref(&chunk[0])` → **one**, keyed on R1 (its comment says so).
  Classification (both out of scope, different scheme) is right; the shapes are as the table
  states only for `:746`. Also "six sites" is main.rs-only — `clump_only.rs:537` and `:1089`
  are two more `pair_output_dir` callers (the clump writers), and io.rs has seven internal uses.
- **N2 — what `integration_output_collision.rs` actually pins.** `:1040-1043` asserts the
  refusal names `reads.fq_trimming_report.txt`; with json first, the first duplicate insert
  would be the `.json` and that assertion goes red. So **txt-before-json is pinned**; mate
  order is pinned nowhere (both mates' txt paths render identically, and no test asserts the
  order of the two source names in the message). The shape that would actually break a test
  is grouping by *kind* (`[all txt, all json]`), not returning four paths per pair. The
  per-mate signature is still right — see A2 — but §10's "the thing most likely to go wrong"
  is narrower than stated.
- **N3 — the stated order is not fully forced.** #429 is independent of both (bam.rs only).
  #424 depends on #423 for tidiness, not correctness — it could hand-roll the strings. Only
  423 → 424 is real, and as a quality dependency. Keep the order, drop "each depends on the last".
- **N4 — "drift becomes unrepresentable" is an overstatement.** After the refactor both sides
  still call the namer independently; what is gone is the `pair_output_dir` + suffix
  derivation, not the possibility of passing the wrong `input`. V1 and the tripwire remain
  the guard. A1 is the version where the claim would be true.
- **N5 — both interleaved arms inline a copy of `pair_output_dir`'s body** (`:1823-1825`,
  `:2427-2429`). Identical semantics today. If `interleaved_report_names` calls
  `pair_output_dir` internally, the arm's primaries and its reports then derive their
  directory from two copies of one rule. One-line fix while you are in there: have the arms
  call `naming::pair_output_dir(input, output_dir)`.
- **N6 — §2.6 is exactly right**, verified against `test_files/ubam_paired_test.bam` across
  `--retain_unpaired`, `--no_report_file`, `--cores 2`, `--gzip`, a `.gz`-named input, and
  `--fastqc` (table in C4). §2.1 likewise: exit 1, zero residue, no `.partial`.
- **N7 — V5 will show five red tests, not three** (C3's two inventory tests join). Worth
  naming so the extras are not misread as the promote trap misfiring.
- **N8 — the `hint` argument for the two new pre-flight calls is unobservable**: only the
  duplicate-path messages use it, and those are unreachable with one input. Either choice is
  defensible; `Some(PAIRED_REPORT_HINT)` matches the siblings, though its "rename one input /
  give the pairs distinct `--output_dir` runs" text would read oddly if it ever fired.
- **N9 — `--gzip` (cli.rs:424, `hide = true`) is not consulted** by main.rs:416's gzip
  resolution, which uses input extension only. Pre-existing, out of scope, plausibly
  deliberate Perl-era compatibility. Mentioned only so I6's case is not mistaken for it.

---

## 6. Action items

### Critical — fix before implementing

1. **C1** — use `OutputSource::Input(input)`, not `Pair`, for every candidate on both
   interleaved arms. `Pair(x, x)` renders `x + x`. Mark assumption 6 and O3 resolved: the
   one-input message is correct as-is under `Input`.
2. **C2** — Step 6 must delete `fn unplanned` (tripwire `:267-270`) along with the three
   `.unplanned()` calls, or CI fails on `dead_code` under `-D warnings`. Keep the variant and
   `compare`'s branch; add a self-test driving `compare(…, Arm::Unplanned)` so the promote
   trap stays exercised.
3. **C3** — add to Step 6: `the_site_inventory_is_complete` 11 → 13, and
   `every_site_is_reached_by_some_case` gains the two interleaved shapes and both counts → 13.
   Also update the three stale doc spots: module Limits (`:63-65`), the `Arm::Unplanned`
   variant doc (`:214-221`), T3's doc (`:317-318`).
4. **C4** — add io.rs unit tests pinning both namers' exact strings, and extend V1 to diff the
   **created file set** for the interleaved arms (baseline table above). Nothing in the suite
   currently pins `_R{1,2}_trimming_report.*`, and post-Step-5 a typo is self-consistent
   across writer and candidate list, so the tripwire cannot see it.

### Important

5. **I1** — V1 compares dump field 3 only, in order; the site column carries a line number
   that the refactor shifts. Env var is `TRIM_GALORE_DUMP_PLANNED_OUTPUTS`.
6. **I2** — give the interleaved *primaries* namers as well (or adopt A1), so #424 does not
   re-introduce hand-rolled duplication in the PR that deletes it.
7. **I3** — pre-flight at the dispatch site (`main.rs` ~`:853`, ~`:2004`), before
   `setup_trimming`, so a refusal does not wait on adapter detection.
8. **I4** — restate V6 as "names the offending read name"; a wrap keeps `invalid input
   parameter` as the cause. Record that the 255-byte case fails with `out of range integral
   type conversion attempted`, and close O2 in favour of wrapping.
9. **I5** — decide the ordinal explicitly (`BamWriter` has no counter; the reference message
   shape has no ordinal). Capture `record.id` in the context closure — `name` is moved into
   `set_name`.
10. **I6** — add a `--paired ip.bam.gz` tripwire case so the candidate list cannot hardcode `.fq`.
11. **V3 is unconstructible as written** — no collision is reachable on either arm today
    (that is #424's premise). Restate it as a unit-level test of the candidate list, or verify
    the machinery on a sibling arm, and say in the plan that the refusal is unreachable by
    construction until something honours `--basename` there.

### Optional

12. **A1** — one `InterleavedOutputs` struct per arm instead of report-only namers. Subsumes
    I2 and removes C4's self-consistency hole. My preference if Felix wants the durable fix.
13. **N1/N2/N3/N4** — correct §2.2's clump descriptions, narrow §2.3's order claim to
    txt-before-json, drop "each depends on the last" from §1, and soften "unrepresentable".
14. **N5** — route both interleaved arms' `dir` through `naming::pair_output_dir`.
