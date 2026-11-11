# Plan: the collision pre-flight names the inputs at fault (#397)

**Issue:** [#397](https://github.com/FelixKrueger/TrimGalore/issues/397) — "Collision
pre-flight cannot name which input produced a colliding report"
**Branch base:** `dev` @ `90dda11`
**Written:** 2026-08-09

---

## Goal

`io::preflight_output_collisions` refuses a run whose planned outputs collide, and its
message tells the user to "rename one input" — but it is given only a list of paths, so it
cannot say *which* inputs are at fault. When the colliding pair are two **reports**, the two
paths are the same string, so the message prints that string twice and names nothing:

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety):
out/reads.fq_clumping_report.txt and out/reads.fq_clumping_report.txt would be written to
the same file. Outputs and reports are named from the input filename alone, so inputs
sharing a filename can collide when --output_dir (or a shared input directory) sends them
to one place — rename one input, or pass --no_report_file if only the reports collide.
```

Thread the originating input(s) through the candidate list so every message can name them.
Behaviour is otherwise unchanged: the same runs are accepted and the same runs are refused,
with the same exit code and the same "nothing is written" guarantee.

Verified live at `90dda11` (debug binary, `--version` self-reports `90dda11`), five probes —
see [Context §Reproductions](#reproductions). One of them (Probe 3) found **live data loss on
an arm this change touches**; see Critical question **C4**.

---

## Context

### The function

`src/io.rs:92-126`. Two bail branches, one shared prefix, one optional caller hint:

| Branch | Trigger | Pinned substring |
|---|---|---|
| output-vs-input (`:104-112`) | a planned output equals a guarded input | `which is also one of its inputs` |
| output-vs-output (`:113-123`) | two planned outputs share a `collision_key` | `would be written to the same file` |

Both open with `Output path collision (case-insensitive, for APFS/NTFS safety): `. The
output-vs-output branch ends with `hint.unwrap_or(GENERIC_ADVICE)` — the hint **replaces**
the generic advice rather than appending, because on the CWD-naming modes the generic advice
("different source directories or `--output_dir`") is false. `GENERIC_ADVICE` is
`src/io.rs:85`; the two hints are `PAIRED_REPORT_HINT` (`src/main.rs:57`) and
`CWD_OUTPUT_HINT` (`src/main.rs:65`).

`collision_key` (`src/io.rs:80`) is case-folded and lexically normalised, so two candidates
can collide while their *rendered* spellings differ (`./x` vs `x`, absolute vs relative,
`Reads.fq` vs `reads.fq`). That distinction drives the message shapes below.

### Call sites — 11 physical, 13 logical

The doc comment's "~11 call sites" claim is exact. All are in `src/main.rs`; nothing outside
`src/main.rs` and `src/io.rs`'s own unit tests calls the function (tree-wide grep).

| # | Line | Dispatch path | Candidate list built by | Hint |
|---|---|---|---|---|
| 1 | `:422` | `--hardtrim5` | `planned_hardtrim_outputs` (`:133`) | `CWD_OUTPUT_HINT` |
| 2 | `:453` | `--hardtrim3` | `planned_hardtrim_outputs` | `CWD_OUTPUT_HINT` |
| 3 | `:595` | `--clump_only` SE FASTQ | inline + `clump_report_candidates` (`:118`) | none |
| 4 | `:637` | `--clump_only --paired` uBAM, Shape B (1 interleaved file) | inline + `clump_report_candidates` | none |
| 5 | `:677` | `--clump_only --paired` uBAM, Shape A (2-file pairs) | inline + `clump_report_candidates` | none |
| 6 | `:731` | `--clump_only` SE uBAM | inline + `clump_report_candidates` | none |
| 7 | `:834` | `--paired` trim, FASTQ out | inline (`:791-833`) | `PAIRED_REPORT_HINT` |
| 8 | `:905` | SE trim, FASTQ out | inline + `planned_secondary_outputs` (`:91`) | none |
| 9 | `:1938` | `--paired` trim, uBAM out | inline (`:1922-1937`) | `PAIRED_REPORT_HINT` |
| 10 | `:1986` | SE trim, uBAM out | inline — **reports missing, see C4** | none |
| 11 | `:2546` | `run_specialty_paired` (`:2530`) | the `pair_outputs` closure | caller-supplied |

Site 11 serves three logical callers, each passing its own closure and hint:
`--clock` (`:483-494`, `CWD_OUTPUT_HINT`), `--implicon` (`:498-519`, `CWD_OUTPUT_HINT`), and
`--clump_only --paired` FASTQ (`:547-581`, `PAIRED_REPORT_HINT`). Its `NameFn` bound
(`:2538`) is therefore part of the signature change.

One dispatch path deliberately has **no** pre-flight: `run_ubam_output`'s
paired-single-interleaved-file arm (`:1902-1913`). With one input there is no
output-vs-output pair, and neither `<stem>_val.bam` nor `<name>_trimming_report.*` can equal
the input. Correct as-is; noted so a reviewer does not read it as a twelfth site.

### Where provenance already exists

At every push site the originating input is in scope — this change adds no plumbing, only
an argument. Four cases are worth stating because one input maps to several candidates, or
one candidate to several inputs:

- **One input → many candidates.** SE FASTQ (site 8): primary + `*_trimming_report.txt` +
  `*_trimming_report.json` + one file per `--demux` barcode, all from the same `input`.
  `planned_secondary_outputs` already loops `for input in &cli.input`.
- **Per-mate reports.** Sites 7 and 9 loop `for input in [&chunk[0], &chunk[1]]`, so each
  report's mate is in hand. `clump_report_candidates` already *takes* the inputs it keys on
  (`report_inputs: &[impl AsRef<Path>]`), so it can return provenance for free.
- **Pair → one output.** `paired_bam_output_name` and `clumped_paired_bam_output_name` take
  both mates but name the file from R1 only (`src/io.rs:215-231`, `:441-457`); their second
  argument is `_input_r2`. `paired_end_output_names` / `clumped_paired_output_names` take the
  **filename** of `_val_2` from R2's stem but the **directory** from R1's parent
  (`src/io.rs:256-271`, `:392-407`). So no single mate explains a paired primary path.
- **Pair → one report.** Site 5 keys its per-pair report on `chunk[0]` alone, because
  `clump_only.rs` writes `clumping_report_name(inputs[0], …)`. R1 is the truthful source
  there; renaming R2 would not resolve a collision. See **C3**.

### Reproductions

All five run the `90dda11` debug binary; scripts under the session scratchpad.

**Probe 1 — the issue's own example.** `--clump_only --paired -o out a/reads.fq b/reads.fq`
→ refused, `out/reads.fq_clumping_report.txt` printed twice, neither `a/reads.fq` nor
`b/reads.fq` named.

**Probe 2 — the worst case in the family.** `--paired --dont_gzip A.fq B.fq B.fq C.fq`
(`B.fq` is R2 of pair 1 and R1 of pair 2 — legal, since `validate_paired_input`
(`src/cli.rs:528-575`) rejects duplicate *pairs* but not a file reused across pairs, and the
SE duplicate-positional check at `src/cli.rs:644-661` is gated on `!self.paired`):

```
Error: Output path collision (…): B.fq_trimming_report.txt and B.fq_trimming_report.txt
would be written to the same file. … rename one input, or pass --no_report_file if only the
reports collide.
```

Both colliding candidates come from **the same input**. Renaming `B.fq` cannot help — the
only followable remedies are listing it once or `--no_report_file`, and the message advises
neither. Refusing is right (before #388 this run silently wrote `B.fq`'s report twice and
kept only the second), but the advice is unfollowable. Reachable with and without `-o`.

**Probe 4 — fold-equal, different spellings.** `--clump_only --paired -o out a/Reads.fq
b/reads.fq` → prints `out/Reads.fq_clumping_report.txt` **and**
`out/reads.fq_clumping_report.txt`. Here printing both paths is the whole message: it is how
the user learns two different-looking names are one file. Pinned by
`clump_paired_rejects_fold_equal_report_names`
(`tests/integration_output_collision.rs:1145-1168`), which asserts
`stderr.contains("_clumping_report.txt")`.

**Probe 5 — output-vs-input.** `--dont_gzip sample.fastq sample.fastq_trimming_report.txt`
→ refused cleanly, the clobbered input named, the 40 reads in it intact. What the message
does **not** say is which input's report would have landed there.

**Probe 3 — out of scope, live data loss (C4).** The same shape as Probe 5 with
`--output-format ubam`: `--output-format ubam sample.fastq sample.fastq_trimming_report.txt`
is **accepted**, writes `sample_trimmed.bam`, then **overwrites**
`sample.fastq_trimming_report.txt` with input 1's trimming report — its 40 reads verified
gone — and only then fails reading it as input 2. Site 10's candidate list omits the reports
that `run_ubam_output_single` writes (`src/main.rs:2011`, `:2099`), while the SE FASTQ twin
covers them via `planned_secondary_outputs`. Fifth instance of the #383/#385/#388/#391
family.

---

## Behavior

### 1. Prerequisites

Unchanged: the pre-flight runs before any output file is opened, after `Cli::validate` and
after `reject_bam_format_mismatch_in_pair`, on every dispatch path that writes more than one
file.

### 2. Main logic

For each planned candidate, in list order (unchanged):

1. Compute `collision_key`.
2. **Output-vs-input.** If the key matches a guarded input, bail naming both the input that
   would be clobbered *and* the input the offending output is named from.
3. **Output-vs-output.** If the key was already seen, bail with one of three shapes, chosen
   by comparing the two candidates' rendered paths and their sources:

   | Shape | Condition | Example |
   |---|---|---|
   | **2a** | paths render equal, sources differ | Probes 1, 2's sibling case, the #388 report collision |
   | **2b** | paths render differently (fold-equal or spelling aliases) | Probe 4, `./x` vs `x`, case-only variants |
   | **2c** | paths render equal **and** sources are the same input | Probe 2; already pinned by `clump_paired_rejects_shared_mate_report_without_output_dir` |

4. Otherwise record the key and continue. `Ok(())` when the list is exhausted.

### 3. Message shapes

Every shape keeps the prefix `Output path collision (case-insensitive, for APFS/NTFS
safety): ` and the substring its tests pin. The pinned substrings are contiguous string
literals in the format template, not assembled from fragments.

**Output-vs-input** — keeps `which is also one of its inputs`:

> Output path collision (case-insensitive, for APFS/NTFS safety): the output named from
> `{source}` would be written to `{input}`, which is also one of its inputs. If that file is
> an earlier run's output, drop it from the input list; otherwise use `--output_dir` to write
> elsewhere.

**2a** — keeps `would be written to the same file`, names both sources, names the path once:

> Output path collision (case-insensitive, for APFS/NTFS safety): the outputs named from
> `{source_a}` and `{source_b}` would be written to the same file, `{path}`. `{hint or
> GENERIC_ADVICE}`

**2b** — keeps both spellings *and* the substring:

> Output path collision (case-insensitive, for APFS/NTFS safety): `{path_a}` (from
> `{source_a}`) and `{path_b}` (from `{source_b}`) would be written to the same file. `{hint
> or GENERIC_ADVICE}`

**2c** — keeps the substring; supplies its own advice and ignores `hint`, on the same
reasoning that already exempts the output-vs-input branch (the hint's "rename one input" is
false here):

> Output path collision (case-insensitive, for APFS/NTFS safety): two outputs named from
> `{source}` would be written to the same file, `{path}`. List each input once — a file that
> appears in more than one pair is reported once per pair — or pass `--no_report_file` if
> only the reports collide.

("in more than one pair", not "as R2 of one pair and R1 of the next": the committed test
`clump_paired_rejects_shared_mate_report_without_output_dir` reaches this shape with one file
passed as R2 of *both* pairs.)

Wording is **O1** (open); the shapes and their pinned substrings are not.

### 4. Source rendering

- one input → `dirA/reads.fq`
- a pair → `dirA/reads_R1.fq + dirA/reads_R2.fq`

### 5. Edge cases

| Case | Handling |
|---|---|
| empty `planned` | `Ok(())`, unchanged |
| one candidate | `Ok(())`, unchanged |
| output equals the input it is *named from* | cannot arise — every namer appends a suffix or replaces the stem (A3). Would render "the output named from X would be written to X"; odd, not wrong |
| same source, both candidates, non-repeated positional | cannot arise — the several candidates one input contributes (primary, `.txt` report, `.json` report, demux files) are pairwise distinct paths (A4). Site 4 is the closest: primary `<stem>_clumped.bam` vs report `<name>_clumping_report.txt`, never equal |
| `Pair(a, b)` vs `Input(a)` sources | sources differ → 2a or 2b; both are named, so the message stays actionable |
| a candidate list whose paths differ only in case | 2b (rendered spellings differ), so both are printed — unchanged from today |

### 6. Non-goals

No change to: which runs are accepted or refused, exit codes, the "a refused run writes
nothing" guarantee, `collision_key`, `path_identity_key`, `guarded_inputs`, the three hint
constants, `GENERIC_ADVICE`, or any output naming function.

---

## Signature

```rust
/// Which input(s) a planned output is named from. The pre-flight cannot tell a user to
/// "rename one input" without it (#397).
#[derive(Clone, Debug)]
pub enum OutputSource {
    /// Named from one input: reports, SE primaries, `--demux` per-barcode files.
    Input(PathBuf),
    /// Named from a pair. Paired primaries take the filename from one mate and the
    /// directory from R1, so neither mate alone explains the path.
    Pair(PathBuf, PathBuf),
}

impl OutputSource {
    /// Case-sensitive identity, for "did these two candidates come from the same input?".
    /// Folding here would call `X` and `x` one source on a case-sensitive filesystem.
    fn identity(&self) -> String { … }
}

impl std::fmt::Display for OutputSource { … }   // `a` | `a + b`

/// A prospective output path and the input(s) it is named from.
pub type PlannedOutput = (PathBuf, OutputSource);

/// Reject duplicate output paths, or an output that would overwrite an input, before any
/// file is opened. Called by every `main.rs` dispatch path that writes more than one file:
/// SE trim, `--paired` trim, `--hardtrim5/3`, `--clock`/`--implicon` and `--clump_only` —
/// all four output formats. A new dispatch path needs a call here too.
///
/// Each candidate carries the input it is named from, so the refusal can name the files the
/// user has to change (#397). The tuple is the point: a candidate cannot be added without
/// stating where it came from.
pub fn preflight_output_collisions(
    planned: &[PlannedOutput],
    inputs: &[PathBuf],
    hint: Option<&str>,
) -> Result<()>
```

Helper signature changes in `src/main.rs`:

```rust
fn planned_secondary_outputs(…) -> Result<Vec<PlannedOutput>>      // :91
fn clump_report_candidates(…)   -> Vec<PlannedOutput>              // :118
fn planned_hardtrim_outputs(…)  -> Vec<PlannedOutput>              // :133
    NameFn: FnMut(&Path, &Path) -> Vec<PlannedOutput>              // :2538
```

---

## Implementation outline

### Step 1 — `src/io.rs`: the type

Add `OutputSource`, its `Display`, its `identity`, and the `PlannedOutput` alias immediately
above `GENERIC_ADVICE` (`:85`). `identity` joins the members' `path_identity_key`s with a
separator that cannot appear in a key rendering, so `Pair(a,b)` and `Pair(b,a)` stay
distinct.

### Step 2 — `src/io.rs`: the function

```rust
pub fn preflight_output_collisions(
    planned: &[PlannedOutput],
    inputs: &[PathBuf],
    hint: Option<&str>,
) -> Result<()> {
    let input_keys: HashMap<String, &PathBuf> =
        inputs.iter().map(|p| (collision_key(p), p)).collect();
    let mut seen: HashMap<String, &PlannedOutput> = HashMap::with_capacity(planned.len());

    for candidate in planned {
        let (path, source) = candidate;
        let key = collision_key(path);
        if let Some(input) = input_keys.get(&key) {
            anyhow::bail!(
                "Output path collision (case-insensitive, for APFS/NTFS safety): \
                 the output named from {} would be written to {}, which is also one of \
                 its inputs. If that file is an earlier run's output, drop it from the \
                 input list; otherwise use `--output_dir` to write elsewhere.",
                source,
                input.display()
            );
        }
        if let Some((prev_path, prev_source)) = seen.insert(key, candidate) {
            // Three diagnoses: one input listed twice, two inputs onto one name, or two
            // spellings of one path.
            if prev_source.identity() == source.identity() {
                anyhow::bail!(/* 2c — own advice, `hint` deliberately unused */);
            }
            if prev_path == path {
                anyhow::bail!(/* 2a */ …, source_a, source_b, path.display(),
                              hint.unwrap_or(GENERIC_ADVICE));
            }
            anyhow::bail!(/* 2b */ …, prev_path.display(), prev_source,
                          path.display(), source, hint.unwrap_or(GENERIC_ADVICE));
        }
    }
    Ok(())
}
```

Order matters: **2c before 2a**, because 2c's condition is a subset of 2a's (same source
implies the paths render equal, given A4). The `seen` map holds `&PlannedOutput` rather than
a clone — the current `p.clone()` is no longer needed.

Update the doc comment. Keep the existing "A new dispatch path needs a call here too"
sentence; add one line on the tuple, per the one-line comment convention.

### Step 3 — `src/main.rs`: the four helpers

- `planned_hardtrim_outputs` (`:133`): the `.map` already has `input`; wrap each path as
  `(path, OutputSource::Input(input.clone()))`.
- `clump_report_candidates` (`:118`): same — the `.map` already has `input`.
- `planned_secondary_outputs` (`:91`): three push points inside `for input in &cli.input`,
  all `Input(input.clone())`; the `--demux` `v.extend(demux::demux_output_paths(…))` becomes
  a `.map` attributing every barcode file to the same `input`.
- `run_specialty_paired` (`:2530`): widen the `NameFn` bound; the body's
  `planned.extend(pair_outputs(&chunk[0], &chunk[1]))` is unchanged.

### Step 4 — `src/main.rs`: the 11 call sites

Per the assignment table below. Mechanical; no control flow changes.

| Site | Candidate | Source |
|---|---|---|
| 1, 2 | hardtrim primary | `Input(input)` |
| 3 | `clumped_output_name` primary; report | `Input(input)` both |
| 4 | `clumped_paired_bam_output_name(input[0], None, …)`; report | `Input(cli.input[0])` both — Shape B has one input, so there is no pair |
| 5 | `clumped_paired_bam_output_name(chunk[0], Some(chunk[1]), …)` | `Pair(chunk[0], chunk[1])` |
| 5 | per-pair report keyed on `chunk[0]` | `Input(chunk[0])` — see **C3** |
| 6 | `clumped_bam_output_name` primary; report | `Input(input)` both |
| 7 | `_val_1` / `_val_2` | `Pair(chunk[0], chunk[1])` — see **C2** |
| 7 | `_unpaired_1` / `_unpaired_2` | `Pair(chunk[0], chunk[1])` (same naming shape) |
| 7 | `passthrough_output_name(pt_input, …)` | `Input(pt_input)` |
| 7 | 4 reports (`.txt` + `.json` × 2 mates) | `Input(mate)` each |
| 8 | SE primary + everything from `planned_secondary_outputs` | `Input(input)` |
| 9 | `paired_bam_output_name` | `Pair(chunk[0], chunk[1])` |
| 9 | 4 reports | `Input(mate)` each |
| 10 | `single_end_bam_output_name` | `Input(input)` — reports still absent, **C4** |
| 11 | `--clock` R1/R2 outputs | `Input(r1)`, `Input(r2)` (each named from its own mate, into the CWD) |
| 11 | `--implicon` R1/R2 outputs | `Input(r1)`, `Input(r2)` |
| 11 | `--clump_only --paired` `_clumped_1` / `_clumped_2` | `Pair(r1, r2)` — see **C2** |
| 11 | `--clump_only --paired` per-mate reports | `Input(r1)`, `Input(r2)` |

### Step 5 — `src/io.rs` unit tests: re-type, don't re-assert

11 `#[test]` fns, 13 call expressions: `:1035`, `:1040`, `:1042`, `:1050`, `:1065`, `:1079`,
`:1093`, `:1108`, `:1126`, `:1135`, `:1143`, `:1149`, `:1153`. Add one local helper beside
`pb` so the arrays stay readable, and **leave every assertion untouched**:

```rust
fn from(path: &str, src: &str) -> PlannedOutput {
    (pb(path), OutputSource::Input(pb(src)))
}
```

Give each existing case a plausible distinct source, so none of them accidentally lands in
2c. Two need care:

- `preflight_rejects_identical_outputs` (`:1047`) asserts
  `err.contains("x_trimmed.fq.gz")` — 2a prints the path once, so it still holds. Sources
  must differ (e.g. `x.fastq` / `x.fq`) or the test moves to 2c.
- `preflight_rejects_output_that_aliases_an_input` (`:1104`) asserts
  `!err.contains("would be written to the same file")`. The new alias wording contains
  "would be written to `{input}`" — the pinned substring is contiguous and is **not**
  present. Do not reword it to "…would be written to the same file as one of its inputs";
  that flips this assertion. Flag in the PR description.

### Step 6 — new unit tests in `src/io.rs`

One per shape, plus rendering: `preflight_names_both_sources_when_one_path_collides` (2a),
`preflight_names_both_spellings_and_both_sources` (2b),
`preflight_diagnoses_one_input_listed_twice` (2c: asserts the new advice, asserts
`!contains("different source directories")` and that a supplied hint is absent),
`preflight_alias_message_names_the_source`, `output_source_pair_renders_both_mates`.

### Step 7 — new integration tests in `tests/integration_output_collision.rs`

Reuse `tempdir` / `write_fastq` / `run_in` / `assert_rejected_cleanly` / `DUP_MSG`. Three
cases, each asserting the *sources* appear in stderr, which is the property #397 buys:

- `clump_paired_report_collision_names_both_mates` — Probe 1: assert `a/reads.fq` and
  `b/reads.fq` both in stderr, plus `assert_rejected_cleanly(&out, …, DUP_MSG)`.
- `paired_report_collision_names_the_offending_pairs` — the #388 twin: two pairs, same
  filenames, different dirs, `-o out`; assert all four mates (or both pair renderings)
  appear.
- `input_reused_across_pairs_says_list_it_once` — Probe 2 (`--paired --dont_gzip A.fq B.fq
  B.fq C.fq`): assert the 2c advice, assert `B.fq` appears, assert the string `rename one
  input` does **not**. The clump analogue already exists at `:1215` and needs no edit — this
  one covers the trim arm, where the shape is currently untested.

### Step 8 — CHANGELOG

One entry under `### Unreleased` → `#### Bug fixes`, plain register per the docs convention:
what the message said, what it says now, and one sentence that `preflight_output_collisions`
took a new parameter shape (a library-API change for a `pub fn` in a published crate; no
external consumer is plausible, so no compatibility shim — **O2**).

### Step 9 — no CI change

All seven existing greps still match (see Integration). No new CI step: the three new
integration tests cover the new shapes, and the CI collision steps exist to catch
release-build-only behaviour, which this is not.

---

## Efficiency

O(n) in the number of planned candidates, unchanged — one `collision_key` and two hash
lookups per candidate. `n` is small: 4 candidates per pair plus 4 reports, so a 100-pair run
plans ~800 entries.

New cost: one or two `PathBuf` clones per candidate (~50-100 bytes each), i.e. tens of KB in
the worst realistic run, freed at the end of the pre-flight. Borrowing `&'a Path` instead
would avoid the clones but requires a lifetime on `OutputSource` threaded through
`run_specialty_paired`'s `FnMut` closures — which already borrow `cli` — and through three
helper return types. Not worth it against a cost that is invisible beside opening one FASTQ.

`seen` now holds `&PlannedOutput` instead of an owned `PathBuf`, removing today's `p.clone()`
at `src/io.rs:113`. No filesystem access is added; the pre-flight stays purely lexical.

---

## Integration

**Reads:** the same `cli.input`, `cli.passthrough`, `cli.demux`, `output_dir`, `gzip`,
`basename`, `no_report_file` already in scope at each site. No new inputs.

**Writes:** nothing. Error text only.

**Order:** unchanged — after `Cli::validate` and `reject_bam_format_mismatch_in_pair`, before
the first writer is opened.

### Test-pinning exposure

This is the highest-risk part of the change. Two kinds of pin.

**Compile-time (signature).** Breaks exactly: `src/io.rs`'s 13 call expressions in 11 test
fns (Step 5), and `src/main.rs`'s 11 call sites + 4 helper signatures (Steps 3-4). Nothing
else — verified by a tree-wide grep for `preflight_output_collisions`, which returns only
those two files.

**Runtime (message substrings).** Every one is preserved:

| Pin | Where | Uses | Branch reached | Safe? |
|---|---|---|---|---|
| `Output path collision (case-insensitive, for APFS/NTFS safety)` (`PREFIX`) | `tests/integration_output_collision.rs:96`, asserted inside `assert_rejected_cleanly` | all 45 tests in the file | every branch | yes — prefix unchanged |
| `would be written to the same file` (`DUP_MSG`) | same file `:94` | 23 | 2a, 2b, 2c | yes — contiguous literal in all three templates |
| `which is also one of its inputs` (`ALIAS_MSG`) | same file `:95` | 8 | output-vs-input | yes |
| `one invocation per input` | same file `:360` and the `--hardtrim3` twin | 2 | 2a via `CWD_OUTPUT_HINT` | yes — hint constants untouched |
| `!contains("different source directories")` | same file `:634`, `:670` | 2 | 2a with a hint | yes — hint still replaces `GENERIC_ADVICE`; 2c emits neither |
| the colliding **path** must appear in stderr | same file `:1037`, `:1137`, `:1165`, `:1237`, `:1270`, `:1291` | 6 | see next table | yes — all four templates render the colliding path at least once |
| `err.contains("x_trimmed.fq.gz")` | `src/io.rs:1057` | 1 | 2a | yes — 2a prints the path once |
| `!err.contains("would be written to the same file")` | `src/io.rs:1116` | 1 | output-vs-input | yes, but **fragile** — see Step 5 |
| `contains("collision") \|\| contains("would be written to the same file")` | `tests/integration_clump_only_ubam.rs:693` | 1 | 2a/2b | yes |
| `contains("Output path collision")` | `tests/integration_gzip_non_gz_extension.rs:424` | 1 | 2a | yes |
| `grep -q "Output path collision"` | `.github/workflows/ci.yml:610`, `:634`, `:655`, `:673` | 4 | 2a (paired ×2, SE, hardtrim) | yes |
| `grep -q "case-insensitive, for APFS/NTFS safety"` | `ci.yml:638` | 1 | 2a | yes |
| `grep -q "one invocation per input"` | `ci.yml:674` | 1 | 2a via `CWD_OUTPUT_HINT` | yes |
| `grep -q "which is also one of its inputs"` | `ci.yml:692` | 1 | output-vs-input | yes |

No pin exists in `docs/` (grep). Strings in `plans/` are review artifacts, not pins.

Those six path assertions classify as follows — worth stating because **all three duplicate
shapes are already reachable from the committed suite**, so the templates are regression-
constrained, not free:

| Test (`tests/integration_output_collision.rs`) | Shape reached | Why |
|---|---|---|
| `clump_paired_rejects_shared_report_name` `:1018` | 2a | two pairs, same filenames, different dirs, `-o out` |
| `clump_paired_rejects_cross_pair_report_collision` `:1114` | 2a | `a/y.fq` and `b/y.fq` are R2 of pair 1 and R1 of pair 2 — one path, two different sources |
| `clump_paired_rejects_fold_equal_report_names` `:1146` | **2b** | `Reads.fq` vs `reads.fq` render differently, so both spellings print |
| `clump_paired_rejects_shared_mate_report_without_output_dir` `:1215` | **2c** | `d/x.fq` is passed as R2 of *both* pairs, so one path arrives twice from **one source** |
| `clump_paired_rejects_shared_report_name_with_basename` `:1249` | 2a | `--basename` makes the primaries distinct, isolating the report collision |
| `clump_se_rejects_report_that_aliases_an_input` `:1280` | output-vs-input | a previous run's clumping report fed back in |

The 2c row is the one to note: shape 2c is **not** new behaviour and **not** untested. That
test asserts `!ok`, `DUP_MSG`, the colliding path, and that three directories are untouched —
all four of which the proposed 2c template satisfies unchanged. What #397 adds there is the
source name and followable advice. Probe 2 is the same class on the `--paired` trim arm.

(Contrast `multi_pair_names_the_offending_pair_and_writes_nothing`,
`tests/integration_paired_format_guard.rs:339`, which also reuses a file across pairs but
asserts the *format guard's* message — that guard fires first, so it never reaches the
pre-flight. Its comment predates #388.)

### Interaction with #398 — undecided, and this design does not presuppose it

#398 covers the paired-arm directory asymmetry: primaries land in R1's parent, reports beside
each mate's own. Four contact points, none of which forces a choice here:

1. **Provenance stays truthful under either resolution.** A report's *filename* always comes
   from its own mate (`report_name` uses `input.file_name()`, `src/io.rs:481-494`); only its
   *directory* is in question. `Input(mate)` is therefore correct whether reports follow the
   mate or R1.
2. **#398 changes which runs collide, not what the message says.** If reports move to R1's
   parent, `--paired dirA/reads.fq dirB/reads.fq` with no `-o` starts colliding on reports;
   shape 2a then fires and names both mates. Nothing here changes. Conversely
   `clump_report_candidates_do_not_over_reject`
   (`tests/integration_output_collision.rs:1174`) pins today's non-collision and is #398's to
   revisit, not this change's.
3. **Hint wording is #398's, not #397's.** `PAIRED_REPORT_HINT`'s parenthetical "(or a shared
   input directory)" describes today's directory rule. This plan leaves all three hint
   constants byte-identical, so #398 owns any rewording — which also keeps the
   `one invocation per input` and `different source directories` pins untouched.
4. **The one place a #398 answer would change this design** is **C2**: whether paired
   primaries carry `Pair(r1, r2)` or the more pointed `Input(r1)` for `_val_1`. `Input(r1)`
   encodes today's "directory comes from R1" rule; `Pair` does not, and survives either
   resolution. Recommendation is `Pair`, precisely so #397 can land first.

If #398 were to consolidate paired reports into one per pair (the shape site 5 already has),
the motivating example thins out but no shape, type, or call site changes — some `Input(mate)`
sites would become `Input(r1)`.

---

## Assumptions

**A1.** Only `src/io.rs` and `src/main.rs` call `preflight_output_collisions` (tree-wide
grep). Changing a `pub fn`'s parameter type in a published crate is a library-API break; no
plausible external consumer exists for this function, so no shim (**O2**).

**A2.** Every output namer derives its filename from exactly one input — its stem
(`strip_fastq_extensions`) or its full `file_name()`. Verified for all namers reachable from
the 11 sites. Paired primaries take the *directory* from R1 while the filename may come from
R2; that asymmetry is why `OutputSource::Pair` exists rather than a bare `PathBuf`.

**A3.** No planned output ever equals the input it is named from — every namer appends a
suffix or replaces the stem, including under `--basename`. Not enforced; if violated the
alias message reads oddly but stays true.

**A4.** Equal-source duplicates (2c) arise only from a positional listed more than once in a
`--paired`-family run. The several candidates one input contributes are pairwise distinct
paths, so 2c cannot fire on a well-formed single-input list. 2c's wording states the fact
("two outputs named from X") and offers the remedy without asserting the cause, so it stays
accurate even if A4 is ever violated.

**A5.** The three hint constants and `GENERIC_ADVICE` stay byte-identical (#398's territory,
and the CI/test pins depend on them).

**A6.** `cli.input.chunks(2)` pairing, `guarded_inputs`, `collision_key` and
`path_identity_key` are unchanged. Configurable behaviour touched: none. Fixed rules relied
on: hint-replaces-generic; the case-folded output key vs the case-sensitive input key.

---

## Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| V1 | 2a names both sources and the path once | new unit test, sources `x.fastq` / `x.fq` onto `x_trimmed.fq.gz` | both sources in the message; `x_trimmed.fq.gz` appears once; contains `PREFIX` + `DUP_MSG` |
| V2 | 2b keeps both spellings and adds both sources | new unit test, `./x_trimmed.fq` + `x_trimmed.fq` | both spellings and both sources present; contains `DUP_MSG` |
| V3 | 2c diagnoses a repeated input and ignores the hint | new unit test, identical source both sides, `Some("EXTRA-HINT.")` | contains the "List each input once" advice and `DUP_MSG`; contains neither `EXTRA-HINT.` nor `different source directories`; does not contain `rename one input` |
| V4 | the alias branch names the source and keeps its exemption | new unit test + the untouched `preflight_rejects_output_that_aliases_an_input` | contains `ALIAS_MSG` and the source; does **not** contain `DUP_MSG` |
| V5 | `Pair` renders both mates | new unit test | `a + b` |
| V6 | the issue's own repro is now actionable | new integration test = Probe 1 | non-zero exit, `out/` empty, stderr names `a/reads.fq` **and** `b/reads.fq` |
| V7 | the #388 twin names the offending pairs | new integration test, two same-named pairs + `-o out` | non-zero exit, `out/` empty, all four mates named |
| V8 | Probe 2 gets followable advice | new integration test = `--paired --dont_gzip A.fq B.fq B.fq C.fq` | non-zero exit, nothing written, `B.fq` named, "List each input once" present, `rename one input` absent |
| V9 | no assertion regressions | `cargo test` from the crate root | all 45 `tests/integration_output_collision.rs` tests and the 11 re-typed `src/io.rs` pre-flight tests pass with assertions unedited. The anchors for each new shape are `:1215` (2c), `:1146` (2b) and `:1018`/`:1114`/`:1249` (2a) — if a template drops `DUP_MSG` or stops rendering the colliding path, these fail first |
| V10 | the CI greps still match | `cargo build --release`, then run the five collision steps of `.github/workflows/ci.yml:595-696` by hand | all seven greps match; every `ls -A` residue check still empty |
| V11 | lint clean | `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D warnings` | clean; watch `clippy::type_complexity` on the `Vec<(PathBuf, OutputSource)>` returns — the `PlannedOutput` alias exists to keep it quiet |
| V12 | no accidental behaviour change | before/after the change, run the five probe scripts | identical accept/refuse verdicts and identical exit codes; only the message bodies differ (and Probe 3 stays broken until C4 is answered) |

V12 is the one that matters most: this change must be message-only. A candidate list that
gained or lost an entry during the mechanical Step 4 edit would show up here as a verdict
flip, and nowhere else.

---

## Decisions (2026-08-09)

All four Critical questions are answered; the plan is unblocked for implementation.

**C1 — Include shape 2c.** Accepted, as recommended. The deciding argument is not cost but
that the alternative ships a branch whose advice we know cannot be followed: today
`--paired A.fq B.fq B.fq C.fq` is refused with "rename one input", and no rename resolves a
file colliding with itself. A message that cannot be acted on is worse than a terse one. The
case is already reachable and already pinned
(`clump_paired_rejects_shared_mate_report_without_output_dir`), so this adds a message shape
for behaviour that exists rather than new behaviour. Sub-question also accepted: **2c ignores
the call-site `hint`**, on the precedent that the output-vs-input branch already does — the
hint's "rename one input" is precisely the false advice 2c exists to replace.

**C2 — Uniform `Pair` for paired primaries.** Accepted. `Input(r1)` is more precise for
`_val_1` today, but precision bought by encoding a claim about *which mate supplies the
directory* is precision we may have to retract: that claim is exactly what #398 is about. If
#398 later moves reports to follow the primaries, `Pair` remains true while `Input(r1)` may
not. Landing #397 before #398 is decided is worth a slightly wordier `_val_1` message, and
it removes the coupling that would otherwise make two open issues block each other.

**C3 — Keep `Input(chunk[0])` for the clump paired-BAM per-pair report.** Accepted, and
blessed explicitly so a reviewer does not read it as an oversight. The asymmetry is
truthful: `clump_only.rs` genuinely keys that report on `inputs[0]`, so renaming R2 would not
resolve the collision, and the message's only job is to name the file that must change.
Symmetry with the sibling candidates would cost accuracy for appearance.

**C4 — Resolved out of band; #397 stays message-only.** The live data loss this raised was
filed as **#409** and fixed in **PR #410** (merging separately), exactly as the plan
recommended. Two consequences for implementation: the candidate list at that site now
already carries both report paths, so #397 must **rebase onto the merged `dev`** rather than
the `90dda11` it was planned against; and #397's CHANGELOG entry stays a message-quality
entry with no data-loss claim.

**O2 is now easier than the plan assumed.** The `pub fn` signature change is technically
breaking for a library consumer, and #402 — landing in PR #411 — adds a `lib.rs` module doc
declaring the library surface an implementation detail with no semver guarantee. That
pre-authorises this change and its successors, so the CHANGELOG needs one sentence rather
than a justification. O1, O3 and O4 stand as written.

## Questions or ambiguities

### Critical

**C1 — Is shape 2c in scope?** The *case* is already reachable and already pinned
(`clump_paired_rejects_shared_mate_report_without_output_dir`,
`tests/integration_output_collision.rs:1215`, where `d/x.fq` is passed as R2 of both pairs).
What this change makes newly *visible* is that both candidates come from one input — and
today that case gets the family's worst message: `--paired A.fq B.fq B.fq C.fq` is refused
with advice that cannot be followed (Probe 2), since one file cannot be renamed apart from
itself. The fix is ~10 lines plus one test on the trim arm, and the existing clump test
passes the proposed 2c template unchanged.
*Recommendation: include it* — cutting it ships #397 with one branch's advice knowingly
unfollowable. Sub-question, same decision: 2c ignores the call-site `hint`
(as the output-vs-input branch already does) rather than appending to it —
*recommendation: ignore*, because "rename one input" is false in exactly this case.
Note it cannot instead be moved to `Cli::validate` alongside the SE duplicate-positional
check (`src/cli.rs:644-661`): reuse across pairs is harmless under `--no_report_file`, and
only the pre-flight is gated on that flag.

**C2 — `Pair` for paired primaries, uniformly, or the more pointed `Input(r1)` for `_val_1`?**
`_val_1` takes both its stem and its directory from R1, so `Input(r1)` is exact and names one
file instead of two. `_val_2` takes its filename from R2 and its directory from R1, so it
needs `Pair` either way. Uniform `Pair` costs a slightly wordier message for `_val_1`
collisions and buys #398-independence: it encodes no claim about which mate supplies the
directory. **This is the one design point a #398 answer would change.**
*Recommendation: uniform `Pair`, so #397 can land before #398 is decided.*

**C3 — Does the `--clump_only --paired` uBAM per-pair report keep `Input(chunk[0])`?**
`clump_only.rs` writes `clumping_report_name(inputs[0], …)`, so R1 genuinely is the source
and renaming R2 would not resolve a collision — but the user typed a pair, and every other
paired candidate at that site says `Pair`. *Recommendation: `Input(chunk[0])` — truthful
beats symmetric, because the message's only job is to name the file that must change.*
Worth an explicit blessing since a reviewer will read the asymmetry as an oversight.

**C4 — Live data loss on site 10, which this change edits.** Probe 3: `--output-format ubam
sample.fastq sample.fastq_trimming_report.txt` is accepted today, overwrites
`sample.fastq_trimming_report.txt` with input 1's trimming report (40 reads verified
destroyed), and only then fails reading it as input 2. Site 10's candidate list
(`src/main.rs:1981-1986`) omits the reports `run_ubam_output_single` writes
(`src/main.rs:2011`, `:2099`); the SE FASTQ twin covers them via `planned_secondary_outputs`
and refuses cleanly (Probe 5). Fifth instance of #383/#385/#388/#391.
Output-vs-**output** is not independently reachable on that arm (for SE, a report collision
implies a primary collision), so the exposure is exactly the output-vs-input case above.
*Recommendation: file it as its own issue today and let it land first, so #397 stays a
message-only change with a clean CHANGELOG entry — but the implementer will be rewriting
those very four lines, so folding it in is defensible if you prefer one PR.* Your call; it
changes #397's scope and its CHANGELOG framing either way.

### Open

**O1 — Wording of the four message templates.** Shapes and pinned substrings are fixed by
the test/CI inventory above; the prose between them is not. Register per the docs convention
(plain, no slop).

**O2 — CHANGELOG framing for the `pub fn` signature change.** `io` is a `pub mod` of a
published crate, so this is technically breaking for a library consumer. Proposal: one
sentence in the entry, no compatibility shim (A1).

**O3 — `--passthrough`'s separate alias check** (`src/cli.rs:966`, its own
"case-insensitively (for APFS/NTFS safety)" wording) already names both files and is not
touched. Recorded so a reviewer does not flag the two messages' different phrasing as an
inconsistency this change introduced.

**O4 — Does the `run_specialty_paired` doc comment need more than a type update?** It
currently says `pair_outputs` "must return EVERY path the pair will write". With the tuple it
must return every path *and its source*. One clause, not a rewrite.

---

## Self-Review

**Efficiency.** O(n) unchanged, n ≤ ~800 in a large run. The clone-vs-borrow trade is
explicit above and resolved toward clones because a lifetime on `OutputSource` would have to
thread through `run_specialty_paired`'s `FnMut` closures. One clone is *removed*
(`src/io.rs:113`) by holding `&PlannedOutput` in `seen`.

**Logic.** Checked branch ordering: 2c must precede 2a because same-source implies
same-rendered-path under A4, so testing 2a first would swallow every 2c case. Checked that
`identity()` uses the case-**sensitive** `path_identity_key`, not `collision_key` — folding
would call `X` and `x` one source and mis-diagnose a legitimate two-input collision as a
repeated positional on a case-sensitive filesystem. Adjusted the alias-branch wording after
noticing it sat one word away from tripping the negative assertion at `src/io.rs:1116`.

**Edge cases.** Empty and single-candidate lists unchanged. Enumerated the four cases where
one input yields several candidates and confirmed none can produce an equal-path pair
(A4) — Shape B at site 4 is the closest and still cannot. Confirmed the paired-single-file
arm's absence of a pre-flight is correct rather than a missed site. Confirmed no namer can
produce a path equal to its own source (A3).

**Integration.** Built the pin inventory by grep rather than by memory, and mapped each pin
to the branch that reaches it; the two that surprised me were `_clumping_report.txt` (reaches
2b, not 2a, because `Reads.fq` and `reads.fq` render differently) and the negative
`!contains(DUP_MSG)` on the alias branch. Traced all four #398 contact points and isolated
the single design point (**C2**) where its answer matters.

**Remaining risks.** (i) The Step 4 edit is mechanical across ~30 push sites; a mis-attributed
source produces a *wrong but plausible* message that no compile error catches — V6-V8 assert
specific source names, and V12 catches a candidate accidentally added or dropped. (ii) The
*advice* text of each shape is asserted only by the new tests; the existing suite pins the
prefix, `DUP_MSG` and the path, so a template could keep every pin and still give bad advice.
(iii) Probe 3 (**C4**) is a live data-loss bug in a line this change edits; leaving it is a
deliberate scope decision, not an oversight, and needs recording as such wherever it lands.

**Corrected during self-review.** A first pass through the pin inventory found only one
path-content assertion and concluded shape 2c was untested. A second grep found six, one of
which (`clump_paired_rejects_shared_mate_report_without_output_dir`, `:1215`) reaches 2c
today. All three duplicate shapes are therefore regression-constrained by the committed
suite, and 2c's wording was corrected from "R2 of one pair and R1 of the next" to "in more
than one pair", which is what that test actually exercises.
