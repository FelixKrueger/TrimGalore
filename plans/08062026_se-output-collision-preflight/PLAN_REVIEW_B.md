# Plan Review B — Output-collision pre-flight for the four uncovered dispatch paths (#383)

**Plan:** `plans/08062026_se-output-collision-preflight/PLAN.md` (v1, 2026-08-06)
**Reviewer:** B (independent; no shared state with Reviewer A)
**Repo state:** `dev` @ `72624c7`; verified against `target/release/trim_galore` as built at that commit
**Verdict:** the plan is accurate and implementable. One **Critical** gap: after the fix, a two-argument
single-end command line can still lose 100 % of one input's reads and exit 0. Assumption A2 is **false**
in one verified corner, and the plan's chosen error message gives remediation advice that is
**demonstrably wrong** on the `--hardtrim` path it newly guards.

---

## 0. Claim verification

Everything below was checked against source or by running the binary. Reproductions live under
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d0ee0171-fc72-476b-bdba-c46c90a5bacd/scratchpad/review-b/`.

### Confirmed as stated

| Plan claim | Verified |
|---|---|
| `preflight_collision_bam` is path-generic, `main.rs:64-80` | ✅ exactly as quoted, body identical to §4.1 |
| Its three call sites are `main.rs:533`, `:566`, `:614` | ✅ `grep` finds exactly three |
| Doc-comment at `:55-59` documents `resolve_clump_layout` but is attached to `preflight_collision_bam` | ✅ — `:55-59` is the clump-budget prose, `:60-63` the BAM prose, one `///` run, `resolve_clump_layout` at `:82` |
| Four unguarded loops: `main.rs:345`, `:370`, `:784`, `:1858` | ✅ I enumerated all ten `cli.input` loops in `main()`; the other six each sit behind a pre-flight (`:484`→`:499`, `:558`→`:572`, `:611`→`:615`, `:683`→`:737`, `:1803`→`:1823`, `:2415`→`:2432`). **The hole map is complete.** |
| `--clock`/`--implicon` already covered via `run_specialty_paired` | ✅ pre-flight at `:2412-2429` |
| `--clump_only` SE FASTQ already rejects | ✅ ran it: exit 1 + message |
| #383 repro (a): `sample.fastq.gz` + `sample.fq.gz` | ✅ exit 0, one `sample_trimmed.fq.gz`, ALPHA=0 / BETA=40, two reports each claiming 100 % written |
| Repro (c): `--hardtrim5 20` on same basename in two dirs | ✅ exit 0, one `same.20bp_5prime.fq.gz` in CWD, DIRA=0 / DIRB=40 |
| SE trim uBAM hole | ✅ exit 0, one `sample_trimmed.bam`, two reports |
| SE trim on same basename in different dirs is **accepted** today and correct | ✅ two outputs, DIRA→dirA only, DIRB→dirB only |
| A3: `hardtrim5` block `return`s at `main.rs:367` | ✅ and no `Cli::validate` conflict between `--hardtrim5`/`--hardtrim3` (Open 4 real) |
| A4: `gzip` decided once at `main.rs:295` from `cli.input[0]`, threaded to every writer | ✅ `run_single_file`/`specialty::hardtrim*` all take the one flag |
| A5: `cli.rs:620` rejects `--basename` + multi SE input | ✅ |
| A6: same file twice exits 0 with correct output today | ✅ ran it: exit 0, 40 reads |
| `cli.rs:549` cross-pair duplicate rejection | ✅ (and `:535` R1≠R2) |
| Namer signatures in §4.2 | ✅ byte-exact, `specialty.rs:423` and `:446` |
| Step 7/8 candidate expressions match what the writers use | ✅ `run_single_file` (`main.rs:1086-1087`) and `run_ubam_output_single` (`:1880-1881`) build their paths with the identical calls |
| Specialty namers resolve to CWD (`specialty.rs:435/456/471/504`) | ✅ |
| `report_name`/`json_report_name` key on the full input filename (`io.rs:401-430`) | ✅ |
| No CI validation step passes colliding stems | ✅ I read every `trim_galore` invocation in `ci.yml`; all multi-input ones are `--paired` with `A_`/`B_`/`C_` or the two deliberate collision guards |
| CI guards assert *empty*, not *absent* | ✅ `test -z "$(ls -A …)"` at `:615` and `:640` |
| `CHANGELOG.md` has `### Unreleased` → `#### Bug fixes` | ✅ at `:4`/`:6` |
| No crate-level `missing_docs`, no `[lints]` table | ✅ so `fn`→`pub fn` on an undocumented namer will not trip `-D warnings` |

### Wrong or imprecise

1. **"five hand-rolled copies" is four** (§2.4, and Out of scope 3 depends on the count).
   The list `(main.rs:482, :680, :1801, :2414, and the --clump_only SE FASTQ one at :486)` double-counts:
   `:482` **is** the `--clump_only` SE FASTQ copy and `:486` is the `insert` line inside that same block.
   The four distinct sites are `:482-498`, `:680-727`, `:1800-1820`, `:2412-2429`.

2. **Message-drift description is off by one.** Two sites say `--output-dir` (`:721` paired FASTQ,
   `:1815` paired uBAM), not one; three say `--output_dir` (`:493`, `:2423`, and the shared helper `:73`).
   More usefully: **`--output-dir` is not a valid flag.** `cli.rs:152` declares `long = "output_dir"` only,
   and `trim_galore --output-dir …` fails with clap's usage error. So the deferred consolidation commit
   (Out of scope 3) is *not* behaviour-neutral in the user-visible sense — it fixes a wrong flag name in
   two error messages. Worth saying so, because it raises that commit's priority. The plan's choice of
   `--output_dir` as the canonical wording is the correct one.

3. **§3.1 point 1 overstates the ordering.** "Runs before any reader is opened" is not literally true:
   `sanity_check_any(&cli.input[0])` (`main.rs:153`) opens a reader and `detect_input_format` runs over
   *every* input at `:159`, both before dispatch. §7's phrasing ("before every reader open on the guarded
   path") is the accurate one. Consequence worth one line in §3.3: a corrupt or missing input **N** is
   diagnosed before the collision, so error precedence is format/sanity → collision.

4. **§3.4 anchors.** `:604`/`:628` are the `trim_galore` invocation lines; the `test -z "$(ls -A …)"`
   assertions the sentence describes are at `:615`/`:640`.

---

## 1. Logic review

### 1.1 Critical — the fix does not close the reported failure class

§1 states the goal as making it impossible for one invocation to silently overwrite one input's output.
An output-vs-output check does not achieve that, because an output can also collide with an **input of
the same invocation**. Verified:

```
e4/s.fastq.gz          (40 reads, ALPHA)
e4/s_trimmed.fq.gz     (40 reads, BETAPRIOR — a previous run's output sitting in the directory)

$ trim_galore e4/s.fastq.gz e4/s_trimmed.fq.gz
exit=0
s_trimmed.fq.gz        : ALPHA=40 BETAPRIOR=0
s_trimmed_trimmed.fq.gz: ALPHA=40 BETAPRIOR=0     ← wrong sample's reads, exit 0
```

Input 1's output path *is* input 2. Input 1 is processed first, clobbers input 2 on disk, and input 2 is
then read back as ALPHA. `BETAPRIOR`'s reads are gone and **both** output files claim success. The
planned pre-flight accepts this invocation, because `s_trimmed.fq.gz ≠ s_trimmed_trimmed.fq.gz`.

This is not a contrived ordering: `trim_galore *.gz` in a directory that already contains a previous
run's output reproduces it, and the glob order (`s.fastq.gz` < `s_trimmed.fq.gz`, `.` = 0x2E < `_` = 0x5F)
puts the clobberer first. Severity is at least #383's: #383 loses one input's reads; this produces a
*second* output file populated with the wrong sample's reads, which is harder to notice than a missing file.

§10 Out of scope 2 does not cover this, or at least will not be read as covering it: it is written about
"files already on disk **from an earlier run**" and frames the answer as "a `--force` / no-clobber policy".
No policy is needed here. The user has explicitly named the file as an input, so writing it is
unambiguously wrong, and the check is a set intersection over data the pre-flight already has:

```rust
// alongside the duplicate-output loop, in the same helper or a sibling
let inputs: HashSet<String> = cli.input.iter().map(|p| naming::norm_path(p)).collect();
// then, per planned path: if inputs.contains(&norm_path(p)) → bail with a distinct message
```

Cost is ~6 lines and one `HashSet`, and it needs its **own** message (the "would be written to the same
file" wording is wrong for this case). It applies to all four newly guarded paths — and to `--paired`,
which has the same hole today, so the helper should take the input list and the whole family benefits.

**Recommendation:** add it to this plan (it is the same commit, the same helper, the same test file), or
if scope must stay tight, replace Out of scope 2 with two entries — one for the `--force`/no-clobber
policy, one naming *this* case explicitly with its own issue number — so that a future reader does not
conclude #383's class was closed.

### 1.2 The insertion points are correct and the candidate expressions match the writers

I traced each of the four:

- `main.rs:344` / `:369` — `output_dir` and `gzip` are both in scope (bound at `:295-296`), `cli` is owned
  so `&cli` is fine, and the check precedes `specialty::hardtrim*`, hence every reader and writer on that
  path. A3 holds, so a per-block check is sufficient.
- `main.rs:781` — precedes `setup_trimming` (which is where adapter auto-detection reads up to 1 M
  records) and `run_single_file`. The candidate expression in Step 7 is *character-identical* to
  `main.rs:1086-1087`, so the checked path cannot diverge from the written path.
- `main.rs:1857` — same property against `main.rs:1880-1881`.

This last property is the one that matters most and the plan gets it right for SE. It is **not** true of
the hardtrim insertions, see 1.3.

### 1.3 The `"5prime"` / `"3prime"` literal is now duplicated, and nothing tests the duplicate

After Steps 5–6 the discriminator string exists twice per mode: once in the collision candidate
(`planned_hardtrim_outputs(&cli, n, "5prime", …)`) and once in the writer (`specialty::hardtrim5` →
`hardtrim_output_name(input, keep, "5prime", …)`). A copy-paste slip in Step 6 that passes `"5prime"` for
`--hardtrim3` produces a pre-flight that hashes paths the run never writes. Every rejection test in §9 V2
still passes (both candidates carry the same wrong discriminator, so they still collide), and there is no
hardtrim3 acceptance test to catch it — see §4.2. Two ways out: assert the expected output filename in a
hardtrim3 acceptance test, or have `planned_hardtrim_outputs` take an enum/`cli.hardtrim5.is_some()`
rather than a `&str` so the literal is written once.

### 1.4 The generic error message's remediation advice is wrong on the hardtrim path

The message the plan adopts ends:

> Check that inputs produce distinct output paths (e.g., different source directories or `--output_dir`).

On `--hardtrim5/3`, **both** suggested remedies fail, because the namer ignores the input's parent and
`-o` joins the *stem-only* filename:

```
$ trim_galore --hardtrim5 20 -o e9out e3/dirA/same.fastq.gz e3/dirB/same.fastq.gz
exit=0 → e9out/same.20bp_5prime.fq.gz, DIRA=0 BETA/DIRB=40
```

"Different source directories" is precisely the input state that fails, and `--output_dir` reproduces the
collision inside the `-o` directory. After the fix the user gets a loud error whose two suggestions both
lead back to the same error, with no hint that the real remedies are one invocation per input, distinct
input basenames, or `cd`-ing between runs.

This matters more than a wording nit because the docs advertise multi-file hardtrim:
`docs/src/content/docs/modes/hardtrim.md` and the v0.6.x changelog text carried in
`docs/src/content/docs/reference/changelog.md:955`/`:997` both say the mode "processes **one or more
files**", and neither mentions that output lands in the CWD rather than beside the input. So
`--hardtrim5 20 */*.fastq.gz` over a per-sample-directory layout is a documented pattern that becomes a
hard failure. Rejecting it is right (it silently loses data today). Rejecting it with unusable advice is
not. §7 (Downstream) checks `ci.yml` but not `docs/` — it should.

**Recommendation:** give the hardtrim insertion its own message, or pass a remediation hint into the
helper. One line, e.g. `"…; hard-trimmed output is written to the current working directory, so identical
input basenames collide regardless of --output_dir — run one invocation per input."`

### 1.5 A6's rejection lands in the wrong layer and inherits the wrong message

For `trim_galore x.fastq.gz x.fastq.gz` the planned helper emits the same path twice —
"`x_trimmed.fq.gz` and `x_trimmed.fq.gz` would be written to the same file" — plus the
different-source-directories advice, which is inapplicable. The codebase has already decided against
exactly this: `cli.rs:504-516`'s doc-comment says the cross-pair duplicate check exists to "emit a
precise error rather than the case-insensitive output-collision pre-flight's APFS/NTFS message".

The consistent placement is `Cli::validate`, beside `cli.rs:534-558`, as a single-end analogue of the
duplicate-pair check. Three benefits over doing it in the pre-flight: a precise message; it runs before
`ensure_output_dir`, so A7's empty-directory residue does not apply to this case; and it covers
`--hardtrim`, `--clock`, `--implicon` and `--clump_only` in one place instead of per dispatch path. It
also makes A6 a *CLI-validation* decision rather than a side effect of the collision helper, which is
easier to revisit if a user objects.

I agree with the **decision** to reject (a duplicated positional is a mistyped glob in every realistic
case, and de-duplicating would need `fs::canonicalize`). Only the placement and the message need moving.

### 1.6 Smaller logic notes

- **A7 / §3.4.** Consistent with the paired contract and correct for the SE `-o` case. But the hardtrim
  V2 cases run with `current_dir(tempdir)` and no `-o`, so `ensure_output_dir` is a no-op there and the
  strongest available assertion is "the temp dir is empty" — stronger than the plan's "the prospective
  output file does not exist", and it also pins the absence of reports, which was #383's most confusing
  artifact (two reports, both claiming 100 % written). Worth upgrading V2 to that.
- **Step 7 will not compile as written.** `main.rs` imports `std::path::Path` but **not** `PathBuf`
  (`main.rs:5`); the existing helper spells it `std::path::PathBuf` in full. `let planned: Vec<PathBuf>`
  in Step 7 and `-> Vec<PathBuf>` in §4.3 need either the full path or a new import. Trivial, but the plan
  presents these as paste-ready.
- **`--basename` + `--hardtrim`** is accepted for a 2-file `--paired --hardtrim5` invocation
  (`cli.rs:620` only fires for `!paired`, `:625` only for `>2`), and the hardtrim namers ignore
  `basename` entirely. Pre-existing silent no-op, unrelated to collisions, and consistent between the
  check and the writer — so it does not affect this plan. Not worth fixing here; worth not being surprised by.

---

## 2. Assumptions

- **A1 (case-fold).** Sound, inherited, and CI-pinned (`ci.yml:634`). No comment.
- **A2 (primary output is a sufficient collision key) — FALSE as worded.** See §3 below. The
  *implementation* is unaffected; the *claim* and V6 need narrowing.
- **A3.** Verified. `main.rs:367` returns; no CLI conflict rule exists, so Open 4 is real.
- **A4.** Verified, and stronger than the plan claims: because `gzip` is global, argument *order* cannot
  produce a false accept. `x.fastq.gz x.fastq.bgz` → both `x_trimmed.fq.gz`; reversed → both
  `x_trimmed.fq`. Both collide, both get caught.
- **A5.** Verified.
- **A6.** Decision endorsed, placement and message disputed (§1.5).
- **A7.** Verified against `ci.yml:615`/`:640`.
- **Unstated assumption worth adding:** that the candidate expression and the writer's expression stay in
  sync. It holds by construction for SE (identical calls) and by convention for hardtrim (duplicated
  literal, §1.3). This is the assumption whose violation is *invisible to every rejection test*, so it
  deserves to be named.
- **Unstated assumption worth adding:** that no planned output path aliases a declared input (§1.1). It is
  false today.

---

## 3. Assumption A2 under test — one verified counterexample

The plan's reasoning about `report_name` is correct: reports key on the full input filename and the
trimmed output on the stripped stem, both resolving the directory the same way
(`output_dir` → else `input.parent()`), so `report_name(a) == report_name(b)` implies
`single_end_output_name(a) == single_end_output_name(b)`. I checked the same property for
`json_report_name` (identical shape, `io.rs:417-430`), for `single_end_bam_output_name` (`:109-123`), and
for `--demux` (`demux.rs:142-147` derives its directory as `output_dir` → else `trimmed_file.parent()`,
and its base name from the trimmed file's own name — so distinct primaries give distinct demux outputs).
All hold.

**FastQC does not.** `fastqc.rs:91-96` lets `--fastqc_args` override the output directory:

```rust
"-o" | "--outdir" => { config.output_dir = Some(PathBuf::from(v)); }
```

so the FastQC artifacts' directory is decoupled from the analysed file's directory, while the *filename*
is still derived from the trimmed output's file name. Two inputs whose primaries do not collide therefore
can, and do, produce colliding FastQC output:

```
$ trim_galore --fastqc --fastqc_args "-o e5/shared" e5/dirA/same.fastq.gz e5/dirB/same.fastq.gz
exit=0
primaries: e5/dirA/same_trimmed.fq.gz AND e5/dirB/same_trimmed.fq.gz   ← distinct, accepted
e5/shared: same_trimmed_fastqc.zip, same_trimmed_fastqc.html           ← ONE pair; dirA's was overwritten
```

Negative control (the check was capable of failing): same command with distinct stems `x`/`y` into the
same `-o` yields **four** files — `x_trimmed_fastqc.{zip,html}` and `y_trimmed_fastqc.{zip,html}`.

Assessment. The lost artifact is a QC report — derived and regenerable, not reads — and the hole is
pre-existing on every path including `--paired` (the #216 pre-flight does not hash FastQC paths either).
So this is **not** a reason to widen the candidate list, and I would not add `--fastqc_args`-aware entries:
that would mean parsing `--fastqc_args` twice and would still miss whatever the next `FastQCConfig` field
does. It *is* a reason to fix the plan's text, because A2 is stated as **(fixed)** and V6 is offered as its
test.

**Recommendation:** reword A2 to the property that is actually true and load-bearing — *for every
secondary output whose directory is resolved by the same `output_dir`-else-parent rule as the primary,
primary-path distinctness implies secondary-path distinctness* — and add the `--fastqc_args -o/--outdir`
override as a named, accepted exception (with the note that it costs a QC report, not reads, and pre-dates
this change on all paths). Then §10 gains one out-of-scope entry.

### V6 as specified is close to vacuous

V6 proposes asserting that where `report_name(a) == report_name(b)`, `single_end_output_name(a) ==
single_end_output_name(b)`. The only way to satisfy the premise with `a != b` is equal `file_name()` plus
a directory rule that erases the parent — i.e. `output_dir = Some(_)` — in which case both sides are equal
by construction. The test can essentially only pass, and it does not test A2's actual claim (that no
*secondary* collides while the primary does not).

A test that can fail, and that would have caught the FastQC case, is the contrapositive over a table:
for path pairs whose primaries differ, assert that every secondary name also differs — report, JSON
report, demux base, and `format!("{}_fastqc.zip", primary_stem)` in the *default* (no-`--fastqc_args`)
configuration — across `--basename` / `--dont_gzip` / `-o` on and off. Then add one comment naming the
`--fastqc_args --outdir` exception so a future reader does not "fix" the table by deleting the row.

---

## 4. Validation sufficiency

### 4.1 The claimed pairing property does not hold for two of the six guarded combinations

§9 asserts that every rejection test is paired with an acceptance test *on the same dispatch path*, so an
unconditionally-rejecting helper would fail the suite. Mapping the proposed tests:

| Guarded combination | Rejection (V2/V4) | Acceptance | Verdict |
|---|---|---|---|
| SE trim FASTQ | ✅ | ✅ V3 ×2 (one with read attribution) | pairing holds |
| SE trim uBAM | ✅ | ❌ none in the plan | relies on pre-existing single-input tests |
| `--hardtrim5` FASTQ | ✅ | ✅ V3 bullet 3 | holds |
| `--hardtrim3` FASTQ | ✅ | ❌ **none anywhere** | gap |
| `--hardtrim5` uBAM | ✅ | ✅ existing `ubam_out_hardtrim5_writes_bam` (single input) | holds |
| `--hardtrim3` uBAM | ❌ none | ❌ none | untested both ways |

`--hardtrim3` has **zero** output-producing coverage in the repository today: the only hits in `tests/`
are `integration_adapter2.rs:377-382`, which asserts a `-a2` rejection *message*. Nothing asserts that
`--hardtrim3 N` writes `*.{N}bp_3prime.*` at all. That is precisely the blind spot §1.3 describes, and
Step 6 is a fresh copy-paste of Step 5 into it.

For SE trim uBAM the practical risk is lower (existing single-input tests in `integration_ubam_out.rs`
would fail against a truly unconditional rejecter), but the plan's own stated property is not met, and
neither the plan nor the existing suite has a **multi-input** SE uBAM acceptance case — which is what
catches an over-rejecting candidate list, e.g. hashing the input path instead of the output path.

**Recommendation:** add three acceptance cases, each asserting the *expected output filename*, not just
"two files exist":
1. `--hardtrim3 20` with two distinct stems, `-o <tmp>` → exit 0, both `*.20bp_3prime.fq.gz` present.
2. `--hardtrim3 20 --output-format ubam` (rejection + acceptance) → both `*.20bp_3prime.bam`.
3. SE `--output-format ubam` with two distinct stems → exit 0, both `<stem>_trimmed.bam`.

### 4.2 V2's hardtrim cases need absolute input paths

The plan correctly pins `Command::current_dir(tempdir)` (and §11 records that the first draft littered the
repo root). It does not say that the fixtures must then be addressed absolutely: `test_files/…` is
relative to the crate root, so with `current_dir(tempdir)` the input is not found and the run fails with
"Input file not found". The rejection assertion (`!ok`) would pass for the wrong reason; only the
`stderr.contains("Output path collision")` grep saves the test. Worth one sentence in Step 10 —
`fixture()` must return an absolutised path (or the fixture must be copied into the tempdir) for any case
that sets `current_dir`.

### 4.3 Other validation notes

- **V2** should assert the output directory is *empty* rather than only that the primary is absent
  (§1.6); that also pins the two-reports artifact that made #383 hard to spot.
- **V3 bullet 1** is the right assertion and the right reason. Keep the read-attribution form.
- **V4** is well chosen: `.gz` + `.bgz` is the newest way to reach the bug (#382 shipped yesterday's
  stem-stripping widening) and pins `strip_fastq_extensions` against a future regression.
- **V5** mirrors the existing guards' shape correctly, including `rc=${PIPESTATUS[0]}` and the
  `test -z "$(ls -A …)"` residue check. The SE step should use the *same* file copied under two names
  (`sample.fastq.gz` + `sample.fq.gz`) as stated — note that the hardtrim CI step must `cd` or pass `-o`
  for the same CWD reason as V2, and since `-o` does not rescue the collision (§1.4), `-o` is the right
  choice there.
- **Missing from V7:** re-running the four §2.2 reproductions is right, but add the *positive* control
  from E2 above (SE, same basename, two directories, no `-o` → exit 0, two outputs, correct attribution).
  That is the one existing behaviour this change is most likely to break by over-rejecting, and it is not
  covered by any CI step.
- **No test covers the A6 rejection.** §9 lists no case for "same file listed twice", even though it is
  the plan's only intentional behaviour change. Add one (and if A6 moves to `Cli::validate` per §1.5, it
  belongs in `cli.rs`'s `mod tests` with the other validation cases).

---

## 5. Efficiency

Nothing to dispute. O(n) over a command line's worth of paths, one `String` + one `PathBuf` clone each, no
syscalls; `planned_hardtrim_outputs` allocates one short-lived `Vec`. The observation that the check
precedes the 1 M-record adapter scan is correct and is the only performance statement that matters — a
colliding invocation now fails in microseconds instead of after a full auto-detection pass.

Two micro-notes, neither worth acting on: the map could be a `HashSet<String>` if the error only named the
duplicate (it names both, so the `PathBuf` value is needed); and if the §1.1 input-alias check is added,
build the input `HashSet` once and share it across all four insertion points rather than per block.

---

## 6. Alternatives

1. **One candidate-list builder in `main()` instead of four insertion points.** Compute the full
   prospective-output list immediately after `ensure_output_dir`, dispatching on mode in one `match`, and
   check once. Cost: a chunky `match` in an already-long `main()`. Benefit: a *fifth* dispatch path cannot
   be added without the compiler forcing a decision about its outputs — which is exactly how #383 and the
   three sibling holes came to exist (each new mode grew its own loop, and three of eight remembered the
   check). Given that the plan is fixing the fourth instance of one omission, the structural fix deserves
   at least a paragraph in §10 rather than silence. I would still ship the plan's version now — four
   local checks are reviewable and match the existing idiom — but note the direction.

2. **A writer-level backstop.** A process-global registry of paths already opened for writing, consulted
   in `FastqWriter::create` / `BamWriter::create`, would have caught #383, the hardtrim CWD case, the
   FastQC `--outdir` case in §3, *and* the input-alias case in §1.1, at one choke point that no future
   mode can bypass. Its weakness is fatal for the primary role: it fires after the first file has been
   written, so it cannot deliver the "nothing written" contract the CI guards assert. Correct framing is
   complement, not replacement: keep the fail-fast pre-flight as the user-facing contract, and consider
   the registry later as a defence-in-depth assertion. Explicitly *not* a reason to change this plan.

3. **De-duplicating rather than rejecting for A6.** The plan's rejection of this (canonicalisation cost,
   symlink/bind-mount ambiguity) is well argued and I agree. Worth recording that a `HashSet` fold over
   byte-equal paths — no `canonicalize`, so no syscall — would silently accept `x.fq.gz x.fq.gz` while
   still rejecting `x.fq.gz ./x.fq.gz`, i.e. an inconsistency that is worse than either clean choice.
   That strengthens the plan's position; it is the reason not to take the "three-line change" escape hatch
   offered at the end of A6.

---

## 7. Action items

### Critical

1. **The plan does not close #383's failure class: an output path can alias an *input* of the same
   invocation.** Verified: `trim_galore s.fastq.gz s_trimmed.fq.gz` (reachable via `trim_galore *.gz`
   after a previous run) exits 0, loses 100 % of the second input's reads, and writes a
   `s_trimmed_trimmed.fq.gz` populated with the *first* input's reads. The planned pre-flight accepts it,
   because the two output paths differ. §10 Out of scope 2 is written about a `--force`/no-clobber policy
   and will not be read as covering this. Either intersect the planned outputs with `naming::norm_path`
   over `cli.input` and bail with a distinct message (~6 lines, same helper, applies to `--paired` too), or
   split Out of scope 2 into two entries and name this case explicitly with its own issue number. (§1.1)

### Important

2. **Assumption A2 is false for `--fastqc_args -o/--outdir`.** Verified: two inputs with distinct primary
   outputs produce a single `same_trimmed_fastqc.{zip,html}` in the overridden directory; the first is
   overwritten. Narrow A2 to secondaries that share the primary's `output_dir`-else-parent directory rule,
   record the FastQC override as an accepted exception (a QC artifact, not reads; pre-existing on all
   paths including `--paired`), and add the corresponding §10 entry. Do **not** widen the candidate list.
   (§3)

3. **V6 as specified cannot fail.** Its premise (`report_name(a) == report_name(b)` with `a != b`) forces
   the conclusion by construction. Replace with the contrapositive over a path table — primaries differ ⇒
   report, JSON, demux base and default-configuration `_fastqc.zip` names all differ — across
   `--basename` / `--dont_gzip` / `-o`, with the `--fastqc_args --outdir` exception noted in a comment.
   (§3)

4. **The adopted error message's remediation advice is wrong on the hardtrim path.** Verified: with
   `--hardtrim5 20 -o out dirA/same.fastq.gz dirB/same.fastq.gz` the collision persists inside `-o`, so
   both suggestions ("different source directories", "`--output_dir`") lead back to the same error. Give
   the hardtrim insertion a mode-specific hint. Also extend §7 (Downstream) to `docs/` — the hardtrim page
   and the carried v0.6.x changelog text both advertise "one or more files", so
   `--hardtrim5 N */*.fastq.gz` is a documented pattern that becomes a hard failure. (§1.4)

5. **Move A6's duplicate-input rejection to `Cli::validate`** beside the duplicate-pair check
   (`cli.rs:534-558`), with its own message. As planned, a duplicated positional produces
   "`x_trimmed.fq.gz` and `x_trimmed.fq.gz` would be written to the same file" plus inapplicable advice —
   and `cli.rs:504-516`'s own doc-comment states that duplicate inputs should get a precise error
   "rather than the case-insensitive output-collision pre-flight's APFS/NTFS message". Validating there
   also puts it ahead of `ensure_output_dir` (no empty-dir residue) and covers every mode at once. Add a
   test — §9 currently has none for the plan's only intentional behaviour change. (§1.5)

6. **Close the §9 pairing gaps.** `--hardtrim3` has zero output-producing coverage in the repo today
   (`integration_adapter2.rs:377-382` only checks a rejection message), and `--hardtrim3 --output-format
   ubam` gets neither a rejection nor an acceptance test. Because the `"5prime"`/`"3prime"` literal is now
   duplicated between the collision candidate and the writer, a Step 6 copy-paste slip would pass every
   proposed rejection test undetected. Add filename-asserting acceptance cases for hardtrim3 FASTQ,
   hardtrim3 uBAM (plus its rejection), and multi-input SE uBAM — or eliminate the duplicated literal.
   (§1.3, §4.1)

### Optional

7. §2.4: four hand-rolled copies, not five (`:482` and `:486` are the same block). Out of scope 3's
   framing should also note that consolidating them fixes a real user-facing bug: two of them advise
   `--output-dir`, which `cli.rs:152` does not define — the flag does not parse.
8. §3.1 point 1: "before any reader is opened" is not literally true; `sanity_check_any` (`main.rs:153`)
   and `detect_input_format` over all inputs (`:159`) precede dispatch. Use §7's phrasing, and note in
   §3.3 that sanity/format errors take precedence over collision errors.
9. Step 7 / §4.3 will not compile: `main.rs` does not import `PathBuf` (only `Path`, `main.rs:5`). Spell
   it `std::path::PathBuf` as the existing helper does.
10. Step 10: for any test that sets `Command::current_dir`, the fixture path must be absolutised —
    otherwise the run fails with "Input file not found" and the `!ok` half of the assertion passes for the
    wrong reason.
11. V2: assert the output directory is empty rather than only that the primary is absent — that pins the
    absence of the two misleading reports, which is #383's most confusing artifact.
12. V7: add the E2 positive control (SE, same basename in two directories, no `-o` → exit 0, two outputs,
    correct read attribution). It is the behaviour most at risk from an over-broad candidate list and no
    CI step covers it.
13. §3.4 anchors: the `test -z "$(ls -A …)"` assertions are at `ci.yml:615`/`:640`; `:604`/`:628` are the
    invocation lines.
14. Consider recording alternative 1 (single candidate-list builder in `main()`) in §10 as the structural
    direction, given that this plan is fixing the fourth instance of the same omission.

---

## 8. What the plan gets notably right

Worth saying, because these are the parts a reviewer would otherwise re-derive: the hole map is complete
(I enumerated all ten `cli.input` loops independently and found no fifth gap); the two naming asymmetries
in §2.3 are real and correctly characterised; the SE candidate expressions are character-identical to what
the writers use, which removes the largest class of pre-flight bug; V3's insistence on read attribution
over file existence is exactly the right lesson from #383; and the decision to move the helper into
`io.rs` is correct for the stated reason — the case-fold assertion is unreachable from any portable
integration test, since two case-variant paths cannot coexist on APFS.
