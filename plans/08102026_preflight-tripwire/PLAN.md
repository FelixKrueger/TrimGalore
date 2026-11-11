# PLAN — #414: empirical exact-set tripwire over the output-collision pre-flight

**Issue:** [#414](https://github.com/FelixKrueger/TrimGalore/issues/414)
**Branch base:** `dev` @ `3008b37` (all line numbers were read from that tree and
re-verified by both plan reviewers)
**Type:** test-only, plus one observability hook in `src/io.rs`
**Baseline established empirically:** green — see §2.2. This lands as a guard, not a fix.
**Revision:** v2 — incorporates `PLAN_REVIEW_A.md` and `PLAN_REVIEW_B.md`. See
§13 for what changed and why.

---

## 1. Goal

Make it impossible to ship a writer whose output path never entered the
output-collision pre-flight's candidate list.

The mechanism is empirical: snapshot the tree, run one dispatch arm, diff, and
assert every file the run created was in the candidate list *that arm's
pre-flight actually checked*. The expectation is derived by **observation**, never
by hand-enumeration — hand-enumerating relocates the original knowledge problem
into a test file, where it rots the same way.

The defect has recurred five times (#383, #388, #391, #409, and the report hole
#398 closed). Every instance was a candidate **nobody wrote**, so `OutputSource`
typing (#397) cannot catch it: a type describes candidates that exist. Reviewer A
checked T1 as specified against all five and confirms each would have failed.

**Non-goal:** the structural single-derivation fix (one derivation feeding both
writers and pre-flight, so a writer without a candidate stops compiling). It is
blocked on normalising candidate order without changing which path today's
messages name. This tripwire is the cheaper first step and produces the observed
write sets that a later restructure needs as its specification.

---

## 2. Context

### 2.1 Why the assertion needs a production hook

The candidate list is assembled in `src/main.rs` — a `[[bin]]`. Nothing in
`tests/` can reach it. Every past fix therefore compared writer output against a
*hand-written* expectation, which is why each of the five instances needed a
human to notice the omission first. This is an **observability gap, not a
diligence gap**.

`io::preflight_output_collisions` (`src/io.rs:148`) is the one function all 11
call sites funnel through. A dump emitted there, gated on an environment
variable, exposes the exact list the pre-flight checked — not a reconstruction of
it — from a single place, without touching any candidate-assembly code. Candidate
**order** is therefore untouched, which matters: order decides which path a
collision message names, and existing tests depend on it. That is why #409
deliberately did not reuse `planned_secondary_outputs`.

Alternatives rejected (decided with Felix, 2026-08-10):
- A hidden `clap` flag: `preflight_output_collisions` has no access to `Cli`, so
  the flag must be threaded as a 4th argument through all 11 call sites — editing
  exactly the code we want to leave alone — or stashed in a global, which is an
  env var with extra steps.
- Moving candidate assembly into the library: that *is* the deferred structural fix.
- **Dumping to a file named by the env var** rather than to stderr (Reviewer B's
  alternative 2). It would retire A3's newline caveat and remove the "parse a
  needle out of a noisy stream" surface. Rejected deliberately, not by inheritance:
  stderr reuses the existing `run_in`-style harness with no new plumbing, Reviewer A
  confirmed no other stderr writer can produce the prefix, and a gated `eprintln!`
  is a smaller surface in a released binary than a gated file write.

### 2.2 Observed write sets — the empirical baseline

Probed with a release binary rebuilt from `3008b37` (the one in `target/` was
stamped `38deb1e`; per the session handoff, a stale binary is a live trap). 40
arm/flag combinations, recursive before/after tree diff. Full transcript:
`spikes/` (probe scripts + raw output + spike report).

Two results shape this plan:

**(a) The baseline is green.** Across all 11 pre-flighted arms, every created
file was already in that arm's candidate list:

| Arm | Created | Planned? |
|---|---|---|
| SE trim | `a_trimmed.fq`, `a.fastq_trimming_report.{txt,json}` | ✔ site `1006` + `planned_secondary_outputs` |
| SE trim `--demux` | + `a_trimmed_{sampleA,sampleB,NoCode}.fq`, `a_trimmed_demultiplexing_summary.txt` | ✔ — the summary **is** planned (`demux.rs:156`) |
| paired trim (R1/R2 in different dirs) | `r1/s_R1_val_1.fq`, `r1/s_R2_val_2.fq`, both mates' reports in `r1/` | ✔ site `927`; confirms #398's anchoring |
| paired `--retain_unpaired` | + `s_R{1,2}_unpaired_{1,2}.fq` — created **even when empty** | ✔ |
| paired `--passthrough` | + `s_I1_passthrough.fq` | ✔ |
| hardtrim5/3, clock, implicon | CWD-anchored primaries, **no reports** | ✔ sites `466`/`497`/`2670` |
| clump_only SE / paired / uBAM ×3 | primaries + clumping reports | ✔ sites `661`/`2670`/`706`/`752`/`809` |
| uBAM-out SE / paired-two-file | `.bam` primary + reports | ✔ sites `2109` (#409's reports) / `2047` |

`--clumpify` leaves **no** temp files in the tree. `--basename` with multiple SE
inputs is refused at the CLI. `--clumpify` requires `--cores >= 2`.

**(b) `--fastqc` is the entire known-unplanned set,** and it is a *rule*, not a
list: for each primary output `P`, FastQC writes `<P-stem>_fastqc.html` and
`<P-stem>_fastqc.zip` beside `P` — including for the `--passthrough` carrier and
for a `.bam` primary. `--fastqc_args "-o DIR"` relocates both wholesale into
`DIR`, **which must already exist** — fastqc-rust does not create it. Verified
directly: without the directory the run trims, writes both reports, then fails at
FastQC with exit 1. (The v1 plan inherited a false negative here from a probe
script that created the directory itself; Reviewer A caught it, Reviewer B drew the
opposite conclusion from the same probe, and the disagreement was settled by
re-running without the `mkdir`.)

This matches the residual documented at `src/io.rs:1371-1376`, where the A2
assumption test declines to assert the FastQC formula because "writing the formula
out here would only test the formula." Reviewer A additionally confirmed
`--extract` is **not** in the curated `--fastqc_args` subset (`fastqc.rs:71-72`),
so no `_fastqc/` **directory** can appear and the two-suffix rule is complete.

**(c) `--fastqc` is a silent no-op on five arms.** `fastqc::run` is called only
from `main.rs:1459`, `:1778`, `:1784`, `:1796`, `:2242`, `:2399`, `:2522` and the
`clump_only` drivers. `specialty.rs` contains no reference to it, and
`run_paired_ubam_single_file` has none. No CLI conflict rejects the combination, so
`--hardtrim5 20 --fastqc a.fastq` is accepted, exits 0, and writes no QC output.
This is a product defect, filed separately (§11); here it only scopes the FastQC
guard (§3.3).

### 2.3 The two arms with no pre-flight at all

`grep -n preflight_output_collisions src/main.rs` returns 11 sites. Two
write-capable dispatch paths are not among them:

| Arm | Reached from | Writes |
|---|---|---|
| paired FASTQ-out from one interleaved BAM | `main.rs:854` → `run_paired_ubam_single_file` | `ip_val_{1,2}.fq`, `ip_R{1,2}_trimming_report.{txt,json}`, **and `ip_unpaired_{1,2}.fq` under `--retain_unpaired`** (`main.rs:1837-1839`) |
| paired uBAM-out from one interleaved BAM | `run_ubam_output:2006` → `run_ubam_output_paired_single_file` | `ip_val.bam`, `ip_R{1,2}_trimming_report.{txt,json}` |

`--retain_unpaired` is rejected only for uBAM *output* (`cli.rs:629`), so
`--paired --retain_unpaired ip.bam` is a legal FASTQ-out invocation — two extra
writers on an arm with no pre-flight. The v1 plan's write set missed them; the
follow-up pre-flight issue would have been planned from an incomplete list.

Both hand-roll `{stem}_R{N}_trimming_report.{txt,json}` (`main.rs:1962`, `:2492`)
— a naming scheme produced by no namer in `io.rs` and used nowhere else.

**Their missing pre-flight is safe by construction.** The proof, in the form both
reviewers converged on (crisper than v1's length argument): both arms derive the
stem with `Path::file_stem()` (`main.rs:1820`, `:2424`), which only ever removes a
`.`-prefixed extension, while every appended suffix begins with `_` — `_val_1`,
`_val_2`, `_val.bam`, `_unpaired_N`, `_R{1,2}_trimming_report.*`. So
`stem + suffix` can never equal `file_name = stem + "." + ext`: a `.` would have to
equal a `_`. This holds for any input, with or without `-o`, and case-folding
cannot change it. Confirmed empirically — a BAM named `x_val_1.fq` yields
`x_val_1_val_1.fq`, one named `y_R1_trimming_report.txt` yields
`y_R1_trimming_report_R1_trimming_report.txt`, and the input survived both runs.
With one input, `guarded_inputs` is a singleton (`--passthrough` needs two inputs,
`--demux` is single-end only), so no output-vs-output collision is expressible either.

**The proof rests on an unstated invariant: neither arm reads `cli.basename`** —
verified absent from `main.rs:1810-2004` and `:2413-2560`; both use `file_stem()`
unconditionally. `--basename` decouples the output stem from the input filename,
which is the premise. With basename honoured, `--basename ip` plus an input named
`ip_val_1.fq` aliases and destroys it. This dependency goes in the `Arm::Unplanned`
doc comment and in the follow-up issue, because it converts that issue from
defence-in-depth into a live data-loss defect the day someone threads basename through.

**Decision (Felix, 2026-08-10):** cover these arms with the weaker assertions —
T3 only, not T1 — so their write set is exercised and the absence is named in the
test rather than silently uncovered. Follow-up issue for adding a real pre-flight.

### 2.4 The 11 sites are 13 pre-flighted arms

Site `2670` is `run_specialty_paired`, which serves **three** modes with three
different closures — `--clock` (`:527`), `--implicon` (`:549`), `--clump_only
--paired` (`:605`). The unit of value is "a dispatch path that writes files", so the
inventory is 13 pre-flighted arms + the 2 orphans = **15 arms** over 11 sites:

| # | Arm | Site | Trigger |
|---|---|---|---|
| 1 | hardtrim5 | `466` | `--hardtrim5 N` (both output formats) |
| 2 | hardtrim3 | `497` | `--hardtrim3 N` (both output formats) |
| 3 | clump_only SE FASTQ | `661` | `--clump_only` |
| 4 | clump_only paired FASTQ | `2670` | `--clump_only --paired` |
| 5 | clock | `2670` | `--clock --paired` |
| 6 | implicon | `2670` | `--implicon --paired` |
| 7 | clump_only uBAM paired, shape B | `706` | `--clump_only --paired --output-format ubam`, 1 BAM |
| 8 | clump_only uBAM paired, shape A | `752` | same, 2 files |
| 9 | clump_only uBAM SE | `809` | `--clump_only --output-format ubam` |
| 10 | trim paired FASTQ | `927` | `--paired` |
| 11 | trim SE FASTQ | `1006` | (default) |
| 12 | trim uBAM-out paired two-file | `2047` | `--paired --output-format ubam` |
| 13 | trim uBAM-out SE | `2109` | `--output-format ubam` |
| 14 | **orphan** paired FASTQ-out, 1 interleaved BAM | — (`854`) | `--paired one.bam` |
| 15 | **orphan** paired uBAM-out, 1 interleaved BAM | — (`2006`) | `--paired --output-format ubam one.bam` |

Both reviewers traced every leaf of `main()` independently and confirmed **there is
no third orphan**: every non-test `File::create` / `FastqWriter::create` /
`BamWriter::create` in `src/` sits downstream of one of those 15. The two shapes
that could have been a third are closed — `--clump_only --paired` with one input is
rejected at `main.rs:597`, and `--clock`/`--implicon` go through
`validate_paired_input` (`cli.rs:1004`, `:1007`), which requires an even count.

### 2.5 Two writer implementations behind one candidate list

`cli.cores > 1 || cli.clumpify` (`main.rs:1316`, `:1548`, `:1847`) switches to
`parallel.rs`, whose paired writer opens its own five handles at `:212-216` —
primaries, passthrough, and both unpaired files — entirely separate code from the
serial block. `--cores` does not *add* a writer; it swaps in a **second writer
implementation behind the same hand-assembled candidate list**. That is exactly the
drift shape this tripwire exists to police, and v1's matrix ran every paired case
at the default `--cores 1`, so it never looked at it. The v1 assumption that
`--cores` "adds no writer" was wrong (Reviewer B, I4).

### 2.6 Existing precedents to align with

- `tests/integration_output_collision.rs` (1698 lines, 57 tests) already has
  `assert_dir_holds_only` — a *hand-enumerated* exact-set check. The tripwire is
  its generalisation: same shape, expectation taken from the binary.
- That file's header records why fixtures are built in-test: the specialty arms
  set `current_dir`, under which `test_files/…` does not resolve. The tripwire
  inherits this and **copies** the committed BAM fixtures into the case root.
- `io::pair_output_dir` (`src/io.rs:91`) is the working "one derivation, many
  callers" precedent, and the model for §10's namer appendix.
- `io::collision_key` (`:80`) and `io::PlannedOutput` (`:137`) are `pub`, so the
  test compares using the production equivalence relation.
- **`collision_key` → `lexical_normalise` → `std::path::absolute` (`io.rs:47`)
  absolutises against the *calling process's* CWD.** In the test that is the test
  process, not the child. Keying a raw relative dumped path would therefore
  compare against the wrong root and T1 would pass vacuously — which is why §3.2
  step 1 canonicalises the root and §3.3 keys on `root.join(p)`. `Path::join`
  discards the base when the argument is absolute, so absolute `--output_dir`
  candidates compare correctly too. (Reviewer B verified this is not merely
  correct but *necessary*.)

---

## 3. Behavior

### 3.1 The hook

1. `preflight_output_collisions` gains `#[track_caller]`. At entry — **before**
   building the input-key map or evaluating any collision — it reads
   `TRIM_GALORE_DUMP_PLANNED_OUTPUTS`.
2. If the variable is present with any value, it prints one line per candidate to
   **stderr**:
   `TRIM_GALORE_PLANNED\t<caller-file>:<caller-line>\t<path as constructed>`
3. If absent — every normal run — it does nothing.
4. No `cfg` gate (§8, A4).
5. Paths are dumped as constructed, not canonicalised. `OutputSource` is
   deliberately omitted: the tripwire compares sets of paths.

The call-site field is what makes §3.4's message implementable and what binds a
case to the arm it is meant to exercise (§4.4). `Location::caller()` yields the
source-relative `src/main.rs` spelling that every `unwrap()` already embeds, so the
reproducible-build job's no-absolute-paths rule is untouched. Stable since Rust
1.46; the crate floor is 1.88.

### 3.2 One tripwire case

1. Create a fresh case root under the system temp dir, canonicalised (macOS
   `temp_dir()` yields `/tmp/…` while the child's `getcwd` yields `/private/tmp/…`;
   without canonicalising the comparison passes vacuously — the reason
   `integration_output_collision.rs::tempdir` already does this). **Tags must be
   unique**: that helper keys on `std::process::id()`, identical for all cases under
   `cargo test --release`, and it opens with `remove_dir_all` — a duplicate tag
   deletes a sibling's tree mid-run, and a deletion that removes an unplanned file
   converts a real T1 miss into a pass. Derive the tag from the test function name.
2. Stage that arm's fixtures inside the root. Paired mates carry **identical read
   IDs**. Cases needing gzipped input write real gzip bytes (`write_fastq_gz`), not
   a plain file named `.gz`.
3. Record `before` = files under the root (recursive), and a content hash of every
   member of **`guarded_inputs`** — `cli.input` plus the `--passthrough` carrier
   plus the `--demux` barcode file (`main.rs:74-83`). A clobbered carrier or barcode
   file is exactly the harm class of #389 and #409.
4. Run the binary with `current_dir` = the case root and
   `TRIM_GALORE_DUMP_PLANNED_OUTPUTS=1`, set on the child `Command` only.
5. Require exit 0, printing stderr on failure.
6. `created` = files under the root after the run, minus `before`.
7. `planned` = the dumped `(site, path)` pairs parsed off stderr.
8. Assert T1, T3, T4, T5 below.

### 3.3 The assertions

**T1 — planned coverage** (arms 1–13): `created \ allowed ⊆ planned`, comparing on
`io::collision_key(root.join(p))`.

- The converse (`planned ⊆ created`) is **not** asserted: over-planning is
  legitimate — `--retain_unpaired` plans `_unpaired_N` whether or not anything
  strands, and `--demux` plans a file per barcode. Asserting the converse against a
  closed-form allowance was considered and declined (§11, O2).
- `allowed` is the FastQC rule from §2.2(b), derived **once**, inside `compare`,
  from the args slice as `--fastqc || --fastqc_args` — mirroring production's gate
  at `main.rs:1458`, `:1777`, `:2241`, `:2398`, `:2521`. v1 had the predicate in two
  places (a `fastqc: bool` parameter and args-based prose) with a disjunction on one
  side only, which is two sources of truth for the one predicate whose drift this
  plan calls fatal.
  **The disjunction is load-bearing, and must not be simplified to `cli.fastqc`.**
  `--fastqc_args`' "Implies --fastqc" (`cli.rs:333`) is **doc text only** — there is
  no clap `requires` and no `validate()` mutation, so `cli.fastqc` is `false` on an
  args-only run and the implication exists solely in the five call-site disjunctions.
  Verified: `--fastqc_args "--quiet" a.fastq` with no `--fastqc` runs FastQC and
  writes both `a_trimmed_fastqc.html` and `.zip`. A predicate mirroring `cli.fastqc`
  alone would call that run "FastQC inactive", leave its two artefacts unallowed, and
  go red on a legitimate case — the exact weakening pressure §3.4 exists to resist.
  (Reviewer B's late O-g argued the opposite from the doc comment; declined on this
  evidence.)
- The allowance is self-policing, **scoped to `FASTQC_CAPABLE_ARMS`**: within that
  set the allowed set must be non-empty when the predicate is true (otherwise
  FastQC never ran and the allowance is masking a broken case), and empty when
  false. Outside it — `--hardtrim5/3`, `--clock`, `--implicon`, orphan arm 14 —
  FastQC never runs (§2.2c), so an unscoped rule would go red on a legitimate case,
  and a provably-wrong red is the strongest argument for deleting the guard. The
  test file states the no-op explicitly so nobody adds a specialty `--fastqc` case
  and then "fixes" the invariant.

**T3 — no pre-existing file's bytes changed** (all 15 arms). This replaces v1's
T2 *and* T3. v1's T2 ("no created path aliases an input") was unfalsifiable: an
output that overwrites an input writes to a path already in `before`, so it is
excluded from `created` by construction, and on case-insensitive APFS a case-variant
alias resolves to the same inode and is likewise absent. Both reviewers caught this
independently. The replacement is strictly stronger than the v1 pair, one assertion
shorter, and it catches #409's harm directly. **T3 is what gives the two orphan arms
their coverage** — v1's claim that T2 did was false.

**T4 — liveness** (all arms): `created` is non-empty; and on `Arm::Planned`,
`planned` is non-empty. Not decoration: `main.rs:215-230` routes clap's
`DisplayHelpOnMissingArgumentOrSubcommand` to **stdout with exit 0**, so a case
whose positional list is swallowed by a preceding flag prints help, creates
nothing, and satisfies every `⊆` vacuously *including* step 5's exit-0 gate. Three
of the probe's own cases died on fixture bugs; a fourth landing on the help path
would have been invisible. These two assertions also make every subset relation in
the suite non-vacuous by construction, and turn "the hook still fires" from an
accident of the arithmetic into an asserted property.

**T5 — one arm per run** (all arms): every dumped line carries the **same** site.
This machine-checks the invariant `parse_planned` silently relies on — all 11 sites
sit outside their per-pair loops today, so a run dumps exactly one arm's list. The
day someone moves a pre-flight inside a loop, the parse would start unioning two
arms' lists and T1 would get quietly weaker; T5 makes that a visible failure.

`Arm::Unplanned` skips T1 and T5 and additionally asserts `planned` is **empty**, so
the day either orphan gains a pre-flight the test goes red and tells the implementer
to promote the case.

### 3.4 Failure messages

A T1 failure names the arm, the dumped **site** (from §3.1, no hand-maintained
table), each unplanned path, the full candidate list, and the correct fix:

> `<path>` was created but never planned. Add it to the candidate list at
> `<site>`. Do **not** add it to the FastQC allowance — that is how this tripwire
> dies.

The handoff records the mechanism by which this bug family recurs: when a
candidate-list test goes red, the tempting fix is to weaken the test. #398 hit the
same shape from the other side, where narrowing the candidate list was the wrong fix
for a red `..._do_not_over_reject`.

---

## 4. Signature

```rust
// src/io.rs — beside collision_key

/// Environment variable that makes the pre-flight print the candidate list it is
/// about to check, one path per line, prefixed with [`DUMP_PLANNED_PREFIX`] and
/// tagged with the call site.
///
/// Test hook for `tests/integration_preflight_tripwire.rs`, which compares a run's
/// observed writes against these lines. Any value activates it; unset in normal
/// runs, where the pre-flight prints nothing.
pub const DUMP_PLANNED_ENV: &str = "TRIM_GALORE_DUMP_PLANNED_OUTPUTS";

/// Line prefix for one dumped candidate.
pub const DUMP_PLANNED_PREFIX: &str = "TRIM_GALORE_PLANNED\t";

/// Render the dump lines for one candidate list. Pure — no environment access, no
/// I/O — so it is unit-testable without mutating process state.
pub fn dump_planned_lines(planned: &[PlannedOutput], site: &str) -> Vec<String>;
```

```rust
// tests/integration_preflight_tripwire.rs

/// One arm's observed-vs-planned comparison. Pure, so it can be unit-tested in
/// both directions.
fn compare(
    root: &Path,
    created: &[PathBuf],
    planned: &[(String, PathBuf)],
    guarded_inputs: &[PathBuf],
    args: &[&str],
    arm: Arm,
) -> Result<(), String>;
```

---

## 5. Implementation outline

1. **`src/io.rs`** — add the two `pub const`s and `dump_planned_lines` beside
   `collision_key` (`:80`); add `#[track_caller]` to `preflight_output_collisions`
   (`:148`) and, before `input_keys`:

   ```rust
   if std::env::var_os(DUMP_PLANNED_ENV).is_some() {
       let at = std::panic::Location::caller();
       let site = format!("{}:{}", at.file(), at.line());
       for line in dump_planned_lines(planned, &site) {
           eprintln!("{line}");
       }
   }
   ```

2. **`src/io.rs` unit test** — `dump_planned_lines_renders_one_line_per_candidate`,
   calling the pure function directly. **No `std::env::set_var`**: it is `unsafe`
   under edition 2024, it leaks to sibling threads under `cargo test --release`, and
   `ci.yml:52-55` records as a maintained invariant that *"no `static mut` /
   `OnceCell` / `env::set_var` globals exist — the codebase is nextest-clean."* v1's
   step 2 would have made that comment false. The env gate is tested end-to-end
   instead (step 5).

3. **`tests/integration_preflight_tripwire.rs`** — new file. Helpers:
   `case_root(tag)` (unique tag, canonicalised); `write_fastq(path, id_prefix)` and
   `write_fastq_gz`; `copy_fixture(name, dest)` resolved from
   `env!("CARGO_MANIFEST_DIR")`; `snapshot(root)` recursive, **files only**
   (`ensure_output_dir` at `main.rs:420` creates directories before dispatch, and
   directories are not writes); `hash_guarded_inputs`; `run_arm(root, args)` setting
   `current_dir` + the env var on the child; `parse_planned(stderr) -> Vec<(String,
   PathBuf)>`; `compare(...)` per §4; and `assert_tripwire(...)` driving T1/T3/T4/T5.

4. **Arm cases** — one `#[test]` per row, matching the repo's many-small-named-tests
   convention (nextest parallelises them; a failure names the arm):

   | # | Test | Args | Notes |
   |---|---|---|---|
   | 1 | `se_trim_two_inputs` | `a.fastq b.fastq` | site 1006 |
   | 2 | `se_trim_output_dir` | `-o out a.fastq b.fastq` | SE `-o` branch |
   | 3 | `se_trim_fastqc` | `--fastqc a.fastq` | exercises the allowance |
   | 4 | `se_trim_fastqc_args_outdir` | `--fastqc --fastqc_args "-o qc" a.fastq` | documented residual; **must pre-create `qc/`** or the run exits 1 (§2.2b) |
   | 5 | `se_trim_demux` | `--demux bc.txt a.fastq` | 4 extra writers; real tabs, 8 bp barcodes; carrier/barcode file hashed |
   | 6 | `se_trim_demux_fastqc` | + `--fastqc` | |
   | 7 | `se_trim_no_report_file` | `--no_report_file a.fastq` | SE report gate |
   | 8 | `se_trim_gzipped_input` | `a.fastq.gz` (real gzip) | `.fq.gz` naming + gz writer |
   | 9 | `se_trim_clumpify` | `--clumpify --cores 2 a.fastq` | confirms no temp writers; **assert the `clumpify:` notice on stderr** — `resolve_clump_layout` (`main.rs:172-204`) falls back to plain mode with a warning and exit 0 below the memory floor, which would silently turn this case into a no-op |
   | 10 | `se_trim_from_ubam_input` | `s.bam` | FASTQ out, BAM in |
   | 11 | `paired_trim_split_dirs` | `--paired r1/s_R1.fastq r2/s_R2.fastq` | #398 anchoring; `r2/` must gain nothing |
   | 12 | `paired_trim_multipair` | 2 pairs | per-pair loop |
   | 13 | `paired_trim_output_dir` | `--paired -o out r1/… r2/…` | **new** — `pair_output_dir`'s `output_dir.map` arm and the `Some(&pair_dir)` threading can diverge only here, and that is #398's exact shape |
   | 14 | `paired_trim_retain_unpaired` | `--retain_unpaired`, short R2 so reads really strand | |
   | 15 | `paired_trim_retain_unpaired_cores` | + `--cores 2` | **new** — §2.5's second writer implementation |
   | 16 | `paired_trim_passthrough` | `--passthrough s_I1.fastq …` | carrier candidate + carrier hashed |
   | 17 | `paired_trim_passthrough_cores` | + `--cores 2` | **new** — `parallel.rs:212-216` opens the carrier handle |
   | 18 | `paired_trim_fastqc` | `--paired --fastqc` | 4 allowed files |
   | 19 | `paired_trim_passthrough_fastqc` | both | carrier gets FastQC'd too |
   | 20 | `paired_trim_basename` | `--basename fixed` | primaries collapse, reports do not |
   | 21 | `paired_trim_no_report_file` | `--paired --no_report_file` | **new** — the paired half of the report gate |
   | 22–23 | `hardtrim5_fastq`, `hardtrim5_ubam` | `--hardtrim5 20` ±`--output-format ubam` | site 466, CWD anchor |
   | 24–25 | `hardtrim3_fastq`, `hardtrim3_ubam` | `--hardtrim3 20`, two inputs | site 497 |
   | 26 | `clock_multipair` | `--clock --paired` ×2 pairs | site 2670, no reports |
   | 27 | `implicon` | `--implicon --paired` | site 2670, `.fastq` extension |
   | 28 | `clump_only_se` | `--clump_only a.fastq` | site 661 |
   | 29 | `clump_only_paired_split_dirs` | `--clump_only --paired r1/… r2/…` | site 2670 |
   | 30 | `clump_only_se_fastqc` | + `--fastqc` | |
   | 31 | `clump_only_paired_no_report_file` | + `--no_report_file` | `clump_report_candidates`' early return (`main.rs:128-130`) |
   | 32 | `clump_only_ubam_se` | `--clump_only --output-format ubam s.bam` | site 809 |
   | 33 | `clump_only_ubam_paired_interleaved` | 1 BAM | site 706 |
   | 34 | `clump_only_ubam_paired_two_files` | 2 FASTQ | site 752 — one report, keyed on chunk[0] |
   | 35 | `ubam_out_se` | `--output-format ubam a.fastq` | site 2109 |
   | 36 | `ubam_out_se_fastqc` | + `--fastqc` | FastQC on a `.bam` primary |
   | 37 | `ubam_out_paired_two_files` | `--paired --output-format ubam r1/… r2/…` | site 2047 |
   | 38 | `orphan_paired_fastq_from_interleaved_bam` | `--paired ip.bam` | `Arm::Unplanned` |
   | 39 | `orphan_paired_fastq_retain_unpaired` | + `--retain_unpaired` | **new** — legal here (`cli.rs:629` refuses it only for uBAM out); records the arm's full write set for the follow-up issue |
   | 40 | `orphan_paired_ubam_from_interleaved_bam` | `--paired --output-format ubam ip.bam` | `Arm::Unplanned` |

5. **Self-tests** — five, all permanent:
   - `comparison_reports_an_unplanned_path` / `comparison_accepts_overplanning` —
     `compare` on synthetic sets, both directions.
   - `hook_is_inert_when_env_unset` — run one arm **twice**, var set then unset;
     assert prefixed-line count `> 0` then `== 0`. This is V1 done without any
     in-process env mutation, and it proves the hook fires *through the real binary*.
   - `end_to_end_goes_red_without_the_dump` — same unset run, then assert
     `compare(…, planned = &[], …)` returns `Err`. Drives hook → stderr → parse →
     compare and proves the **whole chain** can go red, permanently, in both CI
     profiles. v1 relied on V2, a manual one-off that is reverted immediately.
   - `every_site_is_reached_by_some_case` — collect the dumped site from one
     representative case per site and assert **11 distinct values**. No hard-coded
     line numbers, so it survives unrelated `main.rs` edits, and it converts "all 11
     sites covered" from prose into an assertion. Dispatch conditions in `main()`
     change often — the `--clump_only --paired` N=1 bail at `main.rs:597` exists
     precisely because a shape silently routed to the wrong arm.
   - A **refused-run case**: args provoking a genuine collision; assert non-zero
     exit **and** `created.is_empty()`. The tripwire tests the pre-flight's coverage;
     this pins its *purpose* — fail before any write, the contract #409 wanted.

6. **Module doc** — what the file guards, the five prior instances, that
   expectations come from the binary rather than the author, that weakening the
   allowance is the wrong fix, the five arms where `--fastqc` is a no-op, and
   `collision_key`'s case-folding limit (§8, A8).

---

## 6. Efficiency

- The hook costs one `var_os` per run. Nothing when unset.
- ~46 binary invocations × 2 profiles on 40-read fixtures. `integration_output_collision.rs`
  already runs 57 such tests, so the precedent is measured. `.config/nextest.toml`
  sets `slow-timeout = 60s, terminate-after = 4`, so a hung case is killed rather
  than hanging the job.
- Six `--fastqc` cases dominate; each was instant in the probe. Default is
  `--cores 1` (`cli.rs:347`), so concurrent children will not oversubscribe a
  2-core runner; the two `--cores 2` cases add negligible load.
- `snapshot` walks ≤ ~20 files; sets are `BTreeSet<String>` keyed on
  `collision_key`.

---

## 7. Integration

- **Reads:** `src/io.rs` (hook), `test_files/ubam_test.bam` and
  `ubam_paired_test.bam` (copied, never read in place).
- **Writes:** `src/io.rs` (2 consts, 1 pure helper, `#[track_caller]`, ~8 lines of
  gated dump, 1 unit test), one new test file. **No candidate-assembly code is
  touched**, so no candidate reordering and no message-content change.
- **CI:** new tests run in both `cargo nextest run` and `cargo test --release`. The
  `validation` md5 job is untouched (the hook is inert without the env var).
- **Reproducible builds:** unaffected — a runtime env read, and
  `Location::file()` is the source-relative spelling rustc already embeds.
- **Downstream:** the observed write sets in §2.2 and `spikes/` are the record the
  deferred structural fix needs as its specification. The tripwire asserts a subset
  relation; it does not itself publish a per-arm write set.

---

## 8. Assumptions

| | Assumption | Fixed or configurable |
|---|---|---|
| A1 | Every write-capable dispatch path is one of the 11 pre-flight sites or one of the 2 orphans. Independently re-traced by both reviewers; no third orphan exists. | Fixed; T3 still applies to a third, but it would not be *listed* |
| A2 | `--fastqc` is the only source of unplanned writes — across the probed set, plus `--extract` absent from the curated subset so no directory case. Residual: `--fastqc_args "--svg"` (`fastqc.rs:81`) was not probed; probe it once during implementation. | Fixed by the scoped self-policing allowance |
| A3 | Paths contain no newline; `path.display()` is lossy for non-UTF-8. No fixture uses either. | Fixed; documented limitation |
| A4 | The hook must **not** be `#[cfg(debug_assertions)]`-gated. CI runs `cargo nextest run` (`ci.yml:69`) **and** `cargo test --release` (`:76`); a debug gate would silently no-op the tripwire in the release run — a false negative that reads as a pass. Both reviewers name this the most important assumption in the plan. | Fixed |
| A5 | Over-planning is legitimate; only `created ⊆ planned` is asserted. | Fixed (§11 O2 declined) |
| A6 | A failed run may leave partial output; the tripwire asserts only on exit 0. A dump is therefore not evidence the run was accepted. | Fixed |
| A7 | Directories are not writes; the snapshot counts files only. | Fixed |
| A8 | Comparison uses `io::collision_key` — the pre-flight's own case-folded, lexically-normalised relation. **Limit:** because it folds case, the tripwire cannot see a candidate differing from the written path only in case, which on Linux is a real defect (writer writes `X`, pre-flight guards `x`). Using the production relation is still right; the limit is recorded in the module doc so it is not rediscovered. | Fixed |
| A9 | BAM fixtures are copied into the case root because specialty arms set `current_dir`. Both are committed. | Fixed |
| A10 | All 11 sites sit outside their per-pair loops, so a run reaches one pre-flight. `parse_planned` unions every dumped line and relies on this. **T5 asserts only the weaker half** — that no run reaches two *distinct* sites. One site invoked repeatedly is undetectable: the dump has no invocation boundary, so one call of 12 candidates and two of 6 are byte-identical. Harmless, because the union of per-pair lists equals the single list. | Assumed for the repeat case; asserted for the two-site case |
| A11 | Neither orphan arm reads `cli.basename` (`main.rs:1820`, `:2424` use `file_stem()` unconditionally). §2.3's safety proof depends on it. | Fixed; recorded in the `Arm::Unplanned` doc and the follow-up issue |

---

## 9. Validation

| | Verify | How | Expected |
|---|---|---|---|
| V1 | Hook fires through the real binary, and is inert unset | `hook_is_inert_when_env_unset` | line count `> 0` set, `== 0` unset |
| V2 | **The whole chain can go red** — one-off control | Comment out `main.rs:918-924` (the `if !cli.no_report_file { … }` block — **not** 919-925, which leaves the `if` unclosed and deletes `planned.extend(candidates)` at `:925`), rebuild **from HEAD in a fresh worktree**, run `paired_trim_split_dirs` | Fails naming all four report paths; revert |
| V3 | Same, permanently | `end_to_end_goes_red_without_the_dump` + the two `compare` self-tests | Pass with no rebuild |
| V4 | The FastQC allowance cannot silently widen | A created `_fastqc.zip` with no FastQC flag in args | Fails |
| V5 | Site coverage is real | `every_site_is_reached_by_some_case` | 11 distinct sites |
| V6 | The pre-flight's purpose, not just its coverage | the refused-run case | non-zero exit **and** nothing created |
| V7 | All cases green on unmodified `dev` | `cargo nextest run --test integration_preflight_tripwire`, then `cargo test --release --test …` | all pass, both profiles |
| V8 | No regression | full `cargo test`; `cargo clippy --all-targets --release -- -D warnings`; `cargo fmt --all -- --check` | pass, zero warnings. Note: the pre-change `#[test]` count is ~589–590 depending on how `cfg`-gated blocks are counted — do **not** treat 589 as an exact gate |
| V9 | Compile check after any signature change | `cargo test --no-run` — `cargo build --release` does **not** compile test targets | clean |
| V10 | Orphan arms stay honest | `orphan_*` cases | `planned` empty; T3 passes |

---

## 10. Appendix — the deferred namer shape (follow-up issue, not this PR)

Both #398 reviewers asked for dedicated paired report namers taking R1, so
writer/candidate drift becomes *unrepresentable* rather than tested-for.

**Decision: its own issue, landing after #414.** The refactor edits the
candidate-assembly code whose order is load-bearing — precisely what the tripwire
polices. Landing both in one PR means the tripwire has never been shown to pass on
unmodified code, so a green suite would prove nothing about the refactor. Different
risk classes too: #414 is test-only plus an inert hook.

```rust
/// Text + JSON trimming-report paths for one mate of a two-file pair.
///
/// Both land in the pair's output directory — R1's parent unless `--output_dir`
/// overrides — while the filename stays keyed on the mate's own input. Taking R1
/// explicitly means no site can forget the directory rule.
pub fn paired_report_names(
    input: &Path, input_r1: &Path, output_dir: Option<&Path>,
) -> (PathBuf, PathBuf)
```
Consumers — 2 candidate sites, 2 writer sites, each currently threading
`pair_output_dir` by hand: `main.rs:919-924` (the candidate loop), `1745-1751`,
`2040-2045`, `2369-2375`.

```rust
/// Text + JSON report paths for a single-file interleaved pair, whose two reports
/// share one input and are discriminated by an `_R1`/`_R2` infix.
pub fn interleaved_report_names(
    input: &Path, output_dir: Option<&Path>,
) -> ((PathBuf, PathBuf), (PathBuf, PathBuf))
```
Consumers — the two hand-rolled `format!` blocks at `main.rs:1962-1971` and
`2492-2501`. This namer is also what a future pre-flight on the orphan arms would
plan from.

These citations deliberately span the **namer-call lines**, not the enclosing blocks
— that is what a reader following them into the follow-up issue needs to see. Only
V2's range in §9 is a block span, because V2 is an instruction you execute.

---

## 11. Questions or ambiguities

Five questions were put to Felix and answered on 2026-08-10. Three before v1:
orphan-arm treatment (→ weak coverage + follow-up issue), the observability hook
(→ env-var dump in the pre-flight), the namer shape (→ its own issue). Two after
review:

- **O1 — the `--fastqc` silent no-op → filed as [#421](https://github.com/FelixKrueger/TrimGalore/issues/421).**
  Both reviewers found it independently; §2.2c documents it and §3.3 scopes the guard
  around it. The scoping is a workaround for the behaviour, not a fix.
- **O2 — the converse assertion (`planned \ created ⊆ {known over-planning rules}`)
  → leave out.** It would catch a *stale* candidate whose writer was removed —
  #398's shape from the other side, which no current test covers — but it needs a
  second allowance, and this review established that the existing allowance was
  already the weak point in the anti-decay design. Adding a second decay surface in
  the same PR is the wrong trade. Candidate for its own issue later.

**Open (non-blocking):**

- **O3 — case count.** 40 arm cases + 5 self-tests, from crossing the
  writer-relevant flag set with the 15 arms, not a full cross-product. `--rename`
  and `--length_1/2` are excluded as adding no writer; `--cores` is **not** excluded
  any more (§2.5).
- **O4 — CI wall-clock.** Six cases invoke bundled FastQC. If it becomes a concern,
  drop 6, 19, 30 and 36; keep 3 (exercises the allowance) and 4 (pins the residual).
- **O5 — partial output from a failed run.** A6 confines the tripwire to exit-0
  runs. Whether a failed run should clean up is a separate question; the refused-run
  case (§5 step 5) anchors the part that is decided.

---

## 12. Self-Review

**Logic.** The assertion direction was the main thing to get right: `created ⊆
planned`, never the converse, because over-planning is legitimate on at least two
arms. Comparing on `io::collision_key` rather than `PathBuf` equality was a second
correction — and Reviewer B established it is not merely tidier but necessary, since
`lexical_normalise` absolutises against the *test* process's CWD, so keying a raw
relative dumped path would compare against the wrong root and pass vacuously.

**Sets, not sequences.** The comparison is over `BTreeSet` keys throughout, and the
hook is placed where it cannot influence assembly. Candidate order is untouched —
the constraint #414 names, and the reason #409 did not reuse
`planned_secondary_outputs`.

**What review changed, and the pattern in it.** Two of the three v1 assertions were
weaker than advertised: T2 could never fail, and T1's liveness rested on an accident
of the arithmetic. Both defects had the same shape — *asserting more than the
mechanism delivered* — and both were invisible from inside the plan because the
prose was checked against intent rather than against set arithmetic. The other class
was factual: a wrong edit range that would have failed to compile in the one step
proving the tripwire works, and a claim that `--cores` adds no writer when it
substitutes an entire writer implementation.

**How this tripwire could die, and what stops it.** The failure mode is not a false
negative, it is a maintainer weakening the allowance when a legitimate new writer
appears. Four things guard it now: the allowance is a rule derived once from the
args; it is scoped to the arms where FastQC can actually run, so it cannot produce a
provably-wrong red; it must be non-empty when active and empty when not; and the
failure message names the correct fix and forbids the wrong one. v1's unscoped rule
would itself have been a decay vector inside the anti-decay mechanism.

**Verified, not assumed.** The baseline, the FastQC unplanned set, the absence of
clumpify temp files, the orphan arms' non-aliasing, and — after review disagreed on
it — that fastqc-rust does *not* create its `-o` directory, were all established by
running the binary, rebuilt from HEAD after finding `target/` stamped `38deb1e`. My
own first attempt to settle the `-o` question returned "exit 1, nothing created" and
looked like evidence against the binary; the `&&` chain had short-circuited on a
preceding `ls`, so the run never happened. Third harness-induced false negative of
the session, and the reason the rule is: prove the check can produce a positive
before believing its negative.

**Remaining risks.**
1. A4 remains the assumption that could silently defeat the exercise: any
   conditional compilation on the hook turns `cargo test --release` into a false pass.
2. A1 rests on an enumeration of write-capable entry points — now independently
   confirmed twice, but a third orphan added later would be unexercised rather than
   caught. The structural fix is what closes this properly.
3. A8's case-folding limit means a case-only candidate/writer divergence is
   invisible on Linux. Recorded, not fixed.
4. The hook is the first `std::env::var` read in `src/`. Inert without the variable,
   documented as a test hook, but a hidden surface in a released binary.

---

## 13. Implementation notes (2026-08-10)

Landed on `feature/414-preflight-tripwire` off `dev` @ `3008b37`.

**What was built.** `src/io.rs`: `DUMP_PLANNED_ENV`, `DUMP_PLANNED_PREFIX`, the pure
`dump_planned_lines`, `#[track_caller]` on `preflight_output_collisions`, the 7-line
gated dump, and one unit test — 68 lines, no candidate-assembly code touched.
`tests/integration_preflight_tripwire.rs`: 50 tests (41 arm cases + 9 self-tests).
Test total 589 → **640**.

**Verified.**

| Check | Result |
|---|---|
| V2 negative control | Disabling the report block at `main.rs:918-924` fails `paired_trim_split_dirs` naming all four report paths, the site (`src/main.rs:928`), and the two candidates the arm did plan. Reverted; `git diff src/main.rs` empty |
| V5 site coverage | `every_site_is_reached_by_some_case` → 11 distinct sites |
| V7 both profiles | `cargo test` 640 ✔ · `cargo test --release` 640 ✔ · `cargo nextest run` 640/640 ✔ (19.9 s) |
| V8 lint | `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean |
| Hook inertness | Release binary, no env var → 0 dump lines. `--version` still emits the em-dash provenance line; `-V` does not |

Tripwire runtime: 9.0 s debug, 0.68 s release. The six FastQC cases dominate.

**Deviations from §5, all additive.**
1. **Six self-tests, not two, for `compare`.** Every assertion inside it now has its own
   proof it can fire: unplanned path, over-planning accepted, empty run (T4), missing
   dump (T4), two sites in one run (T5), stray FastQC artefact (V4). The plan folded V4
   into a binary run; as a `compare` call it is permanent and needs no rebuild.
2. **`hook_is_inert_when_env_unset` and `end_to_end_goes_red_without_the_dump` merged**
   into one test — they share the same two runs, and splitting them would have run the
   binary twice more for no extra coverage.
3. **New case `se_trim_fastqc_args_only`.** §3.3 argues from evidence that
   `--fastqc_args` alone activates FastQC while `cli.fastqc` stays false; this puts that
   evidence in the suite rather than only in the plan.
4. **`Outcome.root` dropped** — nothing read it and CI is `-D warnings`.

**Iteration log.**
- *#1* Dead `Outcome.root` field failed the `-D warnings` gate. Removed rather than
  `#[allow]`-ed; the two tests that need extras read `stderr` and `created`.
- *#2* `cargo fmt` reflowed two expressions in the test file. Applied, rechecked clean.
- *#3* A doc comment on `parse_planned` referenced `one_arm_per_run`, a test name from
  an earlier draft. Repointed at `compare`.
- *#4* `cargo clippy` finished in 0.24 s and reported nothing, which per
  `feedback_shared_tree_branch_switch` is a cache hit rather than a pass, and grepping
  the log for the test-target name proved nothing (cargo prints `Checking <crate>` once,
  not per target). Planted a `len_zero` violation in the test file: clippy failed with
  exit 101 naming `test "integration_preflight_tripwire"`, so the clean pass was real.
  Probe removed, fmt and clippy rechecked.

**Deviations from §4's signatures, recorded after the coverage audit flagged them as
undocumented (G1, G5).**
5. **`compare` dropped `guarded_inputs` and gained `fastqc_capable`.** The gain follows
   §3.3's scoping. The drop moved T3 into the driver, which put it outside the self-test
   net — the substantive consequence, since T3 is the only assertion on the three
   `Arm::Unplanned` cases. Corrected in round 2 by extracting `unchanged()`.
6. **`hash_guarded_inputs` became `hash_file` plus an inline loop**, **`assert_tripwire`
   became the `Tripwire` builder**, and case 27 is named `implicon_pair` rather than
   `implicon`. Cosmetic, but a reader auditing §4 against the code would not find what §4
   describes.

**Not done, deliberately.** Nothing in the plan was skipped. The two unfiled follow-ups
(§10's namers, and a pre-flight for the orphan arms) remain unfiled pending Felix's
call; #421 is filed.

### Round 2 — after dual code review + coverage audit

Coverage verdict was **INCOMPLETE — 5 items** (106 items: 98 DONE, 4 PARTIAL, 1 MISSING,
3 DEVIATED); no item was a false-pass risk in the suite as it stood. Neither code reviewer
raised a Critical. Test total 640 → **643**.

**Both reviewers mutated `src/main.rs` independently** to add a writer with no candidate
and confirmed the suite goes red naming the unplanned paths, the site, and the candidate
list. One mutation was live in the shared tree while the other reviewer and the auditor were
running, so **two agents independently observed 11 cases go red on a writer neither had
planted** — an unplanned V2-class control nobody arranged. All probes were reverted;
`src/main.rs` verified byte-identical to HEAD.

Applied by Reviewer A during review:
- T3 hashes **every pre-existing file**, not only the declared `guarded` list (closes G2).
  The narrow version was complete only by a per-case convention nothing enforced.
- `compare` rejects an allowance that swallows the **entire** created set, which would
  otherwise leave T1 comparing two empty sets. Plus a self-test for it.
- Root containment is a path-boundary test, not a raw string prefix.

Applied in round 2:
- **`unchanged()` extracted as a pure checker with `comparison_detects_a_modified_input`**
  (G1, Reviewer A H1). T3 was the one assertion whose ability to fail rested on a probe
  that gets reverted; it now has a standing proof, covering rewrite **and** deletion.
- **`the_site_inventory_is_complete`** (Reviewer B H1) — asserts `main.rs` has exactly 11
  pre-flight call sites, so a 12th arm cannot appear without touching this suite. The
  oracle is sound: all 11 mentions are calls, zero prose mentions. This is the bug family
  one level up, from "a candidate is missing from a list" to "an arm is missing from the
  suite".
- **T5's comment and A10 corrected** (G3, Reviewer A M1). T5 detects two *distinct* sites;
  it cannot detect one site invoked repeatedly, because the dump carries no invocation
  boundary — one call of 12 candidates and two of 6 are byte-identical. Softened rather
  than adding a boundary to the hook: the union of per-pair lists equals the single list,
  so T1's strength is unaffected, and the hook stays minimal.
- **Module doc gained three limits**: writes escaping the case root are invisible on the
  created side (Reviewer B M2); a new arm with *no* pre-flight is still undetected; and the
  `--svg` residual is closed — the auditor probed it and it writes no `.svg`, so A2's
  two-suffix rule is complete (G4).

Left open, with reasons:
- **Per-case arm identity is unasserted** (A M2 / B M1 — the one finding both raised).
  Section headers name sites in prose; nothing verifies an individual case reached the site
  its header claims, and the three arms sharing site 2670 cannot be separated by the site
  test. Reviewer A hand-verified all 21 shapes route correctly today, so this is decay risk,
  not a present defect. Felix's call: not in this PR.
- Reviewer lows: allowance breadth (stem-matched artefacts), `planned_keys` computed twice,
  case roots not cleaned up, the trailing bare `bool`, `#[doc(hidden)]` on the hook's public
  items, and the hook activating on an *empty* env value.

**Round 2 verification.** `cargo test` 643 · `cargo test --release` 643 ·
`cargo nextest run` 643/643 · clippy `-D warnings` clean (with a real `Checking` line, not
a cache hit) · `cargo fmt --all -- --check` clean.

---

## 14. Revision History

**v2** (2026-08-10) — after dual independent plan review. Changes, by origin:

*Both reviewers, independently:*
- **T2 deleted, T3 promoted** to "no pre-existing file's bytes changed". v1's T2 was
  unfalsifiable (`created = after \ before` excludes an overwritten input), and v1's
  claim that T2 covered the orphan arms was false — T3 does.
- **Site identity in the dump** via `#[track_caller]` + `Location::caller()`, making
  §3.4's message implementable without a hand-maintained arm→site table (which §1
  forbids) and binding each case to its arm. New T5 + `every_site_is_reached_by_some_case`.
- **FastQC guard scoped** to `FASTQC_CAPABLE_ARMS`; predicate derived once from the
  args as `--fastqc || --fastqc_args`, matching production's disjunction.
- **T4 liveness assertions** (`created` non-empty; `planned` non-empty on planned arms).
- **`guarded_inputs`** named explicitly as the hashed set.
- **Unique case-root tags** required, with the shared-PID amplification stated.
- **Crisper non-aliasing proof** (`file_stem` strips a `.`-extension; every suffix
  begins with `_`), plus A11 recording the `--basename` dependency it rests on.

*Reviewer A, unique:*
- The **clap exit-0 mechanism** (`main.rs:215-230`) that makes T4 Critical rather
  than theoretical.
- **Case 4 must pre-create `qc/`** — verified against Reviewer B's contrary claim.
- Case 9 asserts the **`clumpify:` notice**, so a memory-floor fallback cannot
  silently void it.
- The **refused-run case** (pre-flight's purpose), and asserting planned paths
  resolve inside the case root.
- `--no_report_file` on paired and `clump_only` arms (cases 21, 31).
- Paired `--output_dir` case (13) — the only place `pair_output_dir`'s `-o` arm and
  the writer can diverge, i.e. #398's exact shape.

*Reviewer B, unique:*
- **V2's edit range corrected to `918-924`**; as written it left the `if` unclosed
  *and* deleted `planned.extend(candidates)` at `:925`, which drops every candidate
  for the arm — so T1 would have failed naming six paths, not the four V2 predicts,
  and the difference is the signal that the control did what you meant.
- **`--cores > 1` swaps in a second writer implementation** (`parallel.rs:212-216`)
  → new §2.5 and cases 15, 17; v1's "`--cores` adds no writer" retracted.
- The in-process env unit test **replaced** — `env::set_var` is `unsafe` under
  edition 2024 and contradicts `ci.yml:52-55`'s nextest-clean invariant. Split into
  a pure `dump_planned_lines` plus an end-to-end gate test.
- **A10** (one dump per run) surfaced. Round 2 narrowed the claim: T5 asserts the
  two-distinct-sites half only.
- **A8**'s case-folding limit, real gzip bytes for case 8, `--fastqc_args "--svg"`
  as an A2 residual, and the `#[test]`-count caveat in V8.
- Orphan arm 14's **`--retain_unpaired` writers** added to §2.3 and case 39.

*Felix, on review findings:* the converse assertion left out (§11 O2); the
`--fastqc` silent no-op filed as its own issue (§11 O1).

*Declined, with evidence:* Reviewer B appended two items after its report was first
consumed, and both were checked against the tree rather than applied.
- **B's own retraction of its §10 range finding was correct**, and an earlier
  revision of this section had wrongly "corrected" all five citations into block
  spans. Reverted: v1's ranges point exactly at the namer-call lines (`1745`/`1746`
  and `1750`/`1751`, etc.), which is what a reader following them needs. Only
  `919-925` was genuinely wrong, and only because §9's V2 executes it.
- **B's O-g** (key the allowance on `cli.fastqc`, not the args) rests on reading
  `cli.rs:333`'s "Implies --fastqc" as an implemented constraint. It is a doc
  comment; an args-only run leaves `cli.fastqc` false yet still produces both FastQC
  artefacts (verified). Following O-g would have produced a spurious red on a
  legitimate case, so §3.3 keeps the disjunction and now records why.

**v1** (2026-08-10) — first draft. Baseline, FastQC allowance, and orphan-arm safety
established by probe rather than by reading code.
