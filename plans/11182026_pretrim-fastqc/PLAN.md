# PLAN — `--fastqc_raw`: FastQC on the input, before trimming

**Issue:** none filed. Originates from @ewels for nf-core/methylseq, which wants pre-trim QC without a
separate FastQC process — and without publishing a `fastqc-rust` nf-core module.
**Branch base:** `dev` @ `a197f06`
**Date:** 2026-11-18. **Revision 2** — see §12.
**Status:** revised after dual plan review. Awaiting re-review or implementation trigger.
**Scope:** Trim Galore only. With this revision the downstream work reduces to a **container bump**
in the nf-core module; see §7.

---

## 1. Goal

Add an opt-in `--fastqc_raw` that runs the bundled FastQC engine on each **input** file before
trimming starts, so one Trim Galore invocation produces both pre- and post-trim QC — and so that the
two report sets are distinguishable by every downstream consumer **without pipeline-side
configuration**.

That last clause is the whole design problem. The QC itself is nearly free: the engine is already
in-process (`fastqc-rust`, exact-pinned `=1.1.1`) and `fastqc::run` already takes an arbitrary path.

## 2. Context

### What exists today

`fastqc::run` (`src/fastqc.rs:55`) has this shape, and is called from **twelve** sites in `src/`:

```rust
pub fn run(
    output_path: &Path,
    fastqc_args: Option<&str>,
    output_dir: Option<&Path>,
    cores: usize,
    phred64: bool,
) -> Result<()>
```

Seven in `main.rs` (`:1679`, `:1992`, `:1999`, `:2012`, `:2488`, `:2664`, `:2806`) and five in
`clump_only.rs` (`:306`, `:497`, `:498`, `:779`, `:1034`). **Every one passes an output path** — the
parameter is even named `output_path`. This plan adds the first calls that pass an *input* path.

### The three mechanisms this design turns on

All verified at source during review; they are the reason revision 2 looks different from revision 1.

**M1 — MultiQC takes its sample name from inside the report, not from the filename.**
`multiqc/modules/fastqc/fastqc.py` seeds `s_name` from the zip name, then overrides it:

```python
fn_search = re.search(r"Filename\s+(.+)", file_contents)
s_name = SampleName(self.clean_s_name(fn_search.group(1), f))
```

So renaming the report file changes nothing. On a duplicate name MultiQC **skips the zip entirely**
(`log.debug("Skipping … as already parsed")`) or **overwrites** — in either case one report set is
silently dropped, not merged.

**M2 — that `Filename` line is `path.file_name()` of what fastqc-rust was handed.**
`runner.rs:283-286` sets `FileGroup.name` from `path.file_name()`; `:387` copies it to
`file_display_name`; `basic_stats.rs:403` writes `Filename\t{name}`. **Therefore the only lever on
the MultiQC sample name is the name of the path passed in.**

**M3 — `fastqc-rust` can be told the format, overriding extension sniffing.**
`FastQCConfig.sequence_format` (`config.rs:41`) is checked *first* at `bam.rs:400-431` and accepts
`bam`/`sam`/`bam_mapped`/`sam_mapped`/`fastq`. Only when it is `None` does the extension fallback run
(`.bam`/`.ubam` → BAM, `.sam` → SAM, `.fast5` → Fast5, **everything else assumed FASTQ, with no
error**).

TG already knows each input's format: `input_formats` is computed at `main.rs:343-347`.

### Dispatch order in `main()`

1. `cli.validate()?` — `:330`
2. `sanity_check_any(&cli.input[0])?` — `:337` — **first input only** (§3.7)
3. `input_formats` — `:343-347`
4. `--phred64` pre-flight — `:391`, gated on `fastqc_requested() || output_format == UBam`
5. `output_dir` bound — `:557`
6. Specialty dispatch — `:606`–`:722`, each branch `return`s
7. SE trim `:781`; paired path `:1048`

**Every `naming::preflight_output_collisions` runs at `:608` or later** (13 sites). That fact drives
§3.6, and it is why revision 1's single insertion point was wrong.

## 3. Behavior

### 3.1 The flag

`--fastqc_raw` (bool), independent of `--fastqc`.

| Invocation | Raw reports | Trimmed reports |
| --- | --- | --- |
| `--fastqc_raw` | yes | no |
| `--fastqc` | no | yes (today's behaviour) |
| `--fastqc_raw --fastqc` | yes | yes |
| `--fastqc_args X` alone | no | yes — `--fastqc_args` implies `--fastqc` |
| `--fastqc_raw --fastqc_args X` | yes | no — the args configure the raw run |

**Revised at code review.** Revision 2 documented the last row as "both sets", accepting that
`--fastqc_args` implies `--fastqc` unconditionally. That was wrong in a way the plan did not see: it
made `--fastqc_raw --fastqc_args X` refuse on all four specialty modes, with an error naming
`--fastqc` — a flag the user never passed — falsifying this flag's own help text. `fastqc_requested()`
now reads `self.fastqc || (self.fastqc_args.is_some() && !self.fastqc_raw)`: the args configure
whichever run was asked for, and conjure no report the mode cannot produce.

### 3.2 Naming — a `_raw`-suffixed symlink, reports written flat

For each input, TG stages a symlink whose filename carries the marker, runs FastQC on **the
symlink**, then removes it:

```
sample_R1.fastq.gz                     # the input
<output_dir>/sample_R1_raw.fastq.gz    # symlink -> the input, transient
<output_dir>/sample_R1_raw_fastqc.html # the report, flat
<output_dir>/sample_R1_raw_fastqc.zip
<output_dir>/sample_R1_val_1_fastqc.*  # trimmed, unchanged
```

By **M2**, the report's internal `Filename` becomes `sample_R1_raw.fastq.gz`, so by **M1** MultiQC
reads the sample as `sample_R1_raw` — distinct from both `sample_R1` and `sample_R1_val_1`, in
single-end *and* paired. **No MultiQC configuration is required downstream.** This is the property
revision 1 did not have and could not have obtained by renaming reports.

Suffix placement: inserted before the full extension, so `sample.fastq.gz` → `sample_raw.fastq.gz`.
`fastqc-rust` strips known extensions when naming the report (`runner.rs`, tested at `:1070`), giving
`sample_raw_fastqc.*`.

**Reports are written flat — no subdirectory — into `--output_dir`, or beside their own input when
none is given**, exactly as `single_end_output_name` and the report namers do (`io.rs`, #398). With
distinct filenames the separation a subdirectory bought is unnecessary, and flat output is what the
existing nf-core globs already collect (§7).

**The per-input directory is load-bearing, not cosmetic.** Pinning one directory for the whole run —
which revision 2 did — pools every input's reports together, so two inputs sharing a basename across
directories write the same report and one silently wins while the run exits 0. The trimmed outputs do
not collide in that shape, so the pre-flight has nothing to object to. See §13.

### 3.3 Format is declared, not inferred

Every raw run passes `sequence_format` derived from the `input_formats` entry TG already computed
(**M3**): `UnalignedBam → "bam"`, FASTQ family → `"fastq"`.

This removes the entire class of filename/content mismatch. A content-detected uBAM named `reads.dat`
is analysed correctly rather than being read as FASTQ, and a FASTQ misnamed `*.bam` is likewise read
correctly — the mirror case a suffix guard would have missed. **No refusal is needed and none is
added.**

*(Revision 1 refused non-`.bam`-suffixed uBAM. That guard, its error message, its test and its
one-sidedness are all deleted; see §12.)*

### 3.4 `-o` / `--outdir` inside `--fastqc_args`

`fastqc.rs:124-129` applies `--fastqc_args` **after** TG sets `config.output_dir`, so a user-supplied
`-o` silently overrides TG's choice. For the raw path this would move reports somewhere TG never
promised.

**Behaviour:** when `--fastqc_raw` is given together with `-o`/`--outdir` inside `--fastqc_args`,
refuse at the CLI layer, naming the conflict and pointing at `--output_dir`. Cheaper and clearer than
making the precedence conditional.

### 3.5 Threading

Sequential: the raw pass completes before trimming begins. It adds wall time but no concurrent thread
or memory demand, so `fastqc::resolved_threads()`, the `--clumpify` memory layout (`main.rs:249`) and
the nf-core module's `cores = task.cpus - 4` (`- 3` single-end) arithmetic all stay correct unchanged.
`cli.cores` is passed, as at the existing sites.

One call per input forfeits fastqc-rust's cross-file parallelism. Accepted for v1 and recorded as
§10 Open-5 rather than silently inherited.

### 3.6 Placement — after the output-collision pre-flight

`main.rs:368-371` states the convention: pre-flights sit *"ahead of every dispatch branch … so a
rejected run creates nothing"*, and CLAUDE.md calls the collision pre-flight load-bearing. A raw pass
at revision 1's site would have written reports and read every input **before** any collision refusal
— a far longer wait than the 1M-read adapter scan that convention already protects.

**The raw pass therefore runs inside each dispatch branch, immediately after that branch's
`preflight_output_collisions`.** This costs revision 1's "one insertion point serves every mode"
property; that property was the reason for the breach, so it is the thing to give up.

Sites: the four specialty branches (`:606`–`:722`), the SE trim path (`:781`), the paired path
(`:1048`), and `--clump_only`. Each becomes `if cli.fastqc_raw { run_pretrim_fastqc(...)?; }` after
its existing pre-flight call.

The raw reports' own paths are added to the planned-output set passed to
`preflight_output_collisions`, so a raw report colliding with another output is refused by the same
mechanism rather than by a parallel one.

### 3.7 Input validation

`sanity_check_any` runs pre-dispatch on `cli.input[0]` **only** (`main.rs:337`); every other input is
checked inside the drivers (`:982`, `:984`, `:1131`, `:1133`, `:1197`, `:2292`, `:2294`, `:2345`).

Without care, a bad **R2** would reach `fastqc-rust` before TG looked at it, and the error would come
from the QC engine rather than TG's own guard. `run_pretrim_fastqc` therefore calls
`sanity_check_any` over **every** input it is about to analyse, before analysing any of them.

**What that does and does not buy** (measured during implementation, §13): `sanity_check_any`
inspects the *first record*, so it catches colorspace and empty input with no report written at all.
Truncation deep in the stream is invisible to it — the failure surfaces from the QC pass, after
earlier inputs have been reported on. Trim Galore has the same blind spot while trimming, so this is
the existing standard rather than a regression, but the guard must not be described as covering
truncation.

### 3.8 Edge cases

| Case | Behaviour |
| --- | --- |
| Paired input | One report per mate. Distinct names, so no clash. |
| Multiple pairs / multi-file SE | One report per input file. Two inputs sharing a basename collide — now caught by §3.6's pre-flight rather than silently overwriting (revision 1's Open-3). |
| `--output_dir` unset | Reports and the transient symlink in the working directory. |
| `--basename` | Does not apply; raw reports are named after the input + `_raw`. |
| Specialty modes | Allowed. `--fastqc` stays refused there (§3.9). |
| `--clump_only` | Allowed; the five `clump_only.rs` sites are untouched. |
| `--demux` | Allowed; dispatches at `main.rs:1689-1695`. One pre-demux report per input, since demultiplexing happens after. |
| uBAM input, any filename | Analysed correctly (§3.3). |
| `--phred64` | Pre-flight gate widened (§5 step 4). |
| Input already named `*_raw.*` | The symlink would be `sample_raw_raw.fastq.gz`. Ugly but unambiguous and non-colliding. |
| Symlink creation fails (filesystem without symlinks) | Hard error naming the path. See §10 Open-2 and V9. |
| Empty or colorspace input, any position | Refused by §3.7's loop before any analysis; no report written. |
| Truncation deep in the stream | Not detectable without reading the file. The run fails from the QC pass, having reported on earlier inputs (§3.7). |
| FIFO / pipe / stdin input | Already refused by #379 (`cli.rs:536`, `format.rs:342`), which is what makes a second full read of each input safe at all. |

### 3.9 Interaction with the `#421` refusals

`fastqc_requested()` (`cli.rs:609-610`) and the `#421` refusals (`cli.rs:1113-1151`) are **not
modified**, so all existing behaviour and tests are untouched. `--fastqc_raw` is deliberately not
subject to them: they exist because those modes write no post-trim report, and the input is present
and analysable in every mode.

## 4. Signature

```rust
/// Run the bundled FastQC on each input, before trimming.
///
/// Stages a `_raw`-suffixed symlink per input and analyses that, so the
/// report's own `Filename` field carries the marker and MultiQC reads the raw
/// and trimmed sets as distinct samples with no pipeline-side configuration.
/// The symlink is removed afterwards; `*_raw_fastqc.{html,zip}` remain.
///
/// Declares each input's format from `input_formats` rather than letting
/// fastqc-rust infer it from the filename.
///
/// Sanity-checks every input before analysing any of them: `main.rs:337`
/// covers only the first.
fn run_pretrim_fastqc(
    cli: &Cli,
    inputs: &[PathBuf],
    input_formats: &[InputFormat],
    output_dir: Option<&Path>,
) -> Result<()>
```

`fastqc::run` gains one parameter:

```rust
pub fn run(
    output_path: &Path,
    fastqc_args: Option<&str>,
    output_dir: Option<&Path>,
    cores: usize,
    phred64: bool,
    sequence_format: Option<&str>,   // new — sets FastQCConfig::sequence_format
) -> Result<()>
```

All twelve existing call sites pass `None`, preserving today's extension-based dispatch — including
the three uBAM-output sites whose comments rely on it.

## 5. Implementation outline

1. **`cli.rs`** — add `#[clap(long = "fastqc_raw")] pub fastqc_raw: bool` beside `--fastqc` (`:345`).
   Doc comment: runs on the input before trimming; reports are `*_raw_fastqc.*`; allowed in the
   specialty modes that refuse `--fastqc`.
2. **`cli.rs`** — leave `fastqc_requested()` untouched (§3.9). Add
   `pub fn fastqc_raw_requested(&self) -> bool { self.fastqc_raw }` if it reads better at call sites.
3. **`cli.rs::validate()`** — refuse `--fastqc_raw` together with `-o`/`--outdir` inside
   `--fastqc_args` (§3.4). Pure string inspection, so it belongs here rather than in `main.rs`.
4. **`main.rs:391`** — widen the `--phred64` pre-flight gate from
   `cli.fastqc_requested() || output_format == UBam` to include `|| cli.fastqc_raw`.
5. **`fastqc.rs`** — add the `sequence_format` parameter and set `FastQCConfig::sequence_format`
   **before** `apply_fastqc_args` so an explicit format is not clobbered. Update all twelve callers
   to pass `None`.
6. **`main.rs`** — add `run_pretrim_fastqc`:
   a. `sanity_check_any` over every input (§3.7).
   b. Per input: build `<dir>/<stem>_raw<ext>`, `std::os::unix::fs::symlink` to the **absolute**
      input path, call `fastqc::run(&link, args, output_dir, cli.cores, cli.phred64, Some(fmt))`,
      then remove the link — including on the error path, so a failure does not strand it.
   c. `create_dir_all` the output directory first: fastqc-rust does **not** create it
      (`runner.rs:477`, `archive.rs:66` write without creating). Pattern to follow: `fastq.rs:739`.
7. **`main.rs`** — call it in each dispatch branch immediately after that branch's
   `preflight_output_collisions` (§3.6), and add the raw-report paths to the planned-output set.
8. **Tests**
   - `--fastqc_raw` parses; does not imply `--fastqc`; does not set `fastqc_requested()`.
   - Accepted with each of `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon` — the acceptance
     twins beside `cli.rs:2773`; `--fastqc` still refused there.
   - `--fastqc_raw` + `--fastqc_args "-o X"` refused (§3.4).
   - Integration, SE and PE: `*_raw_fastqc.zip` exists **and the in-zip `Filename` line carries
     `_raw`** — the property the whole design rests on (V2).
   - uBAM fixture copied to a non-`.bam` name: report is BAM-derived, not FASTQ-misread (V3).
   - A collision-refused run leaves **no** raw report (V4).
   - Truncated/colorspace R2 is refused by TG, not by fastqc-rust (V5).
9. **Docs** — `--help`, `docs/` flag reference, `wiki/reference/cli-flags.md`, and a one-sentence
   CHANGELOG entry.

## 6. Efficiency

One additional full read of each input, sequential with respect to trimming (§3.5). Safe because
#379 already refuses non-restartable inputs (`cli.rs:536`, `format.rs:342`) — without that, a second
read of a FIFO would consume a stream the trimmer still needs. No change to the trimming hot path.
The symlink is metadata only; no payload copy.

## 7. Integration — downstream reduces to a container bump

### 7.1 nf-core/modules `trimgalore`

Verified against current upstream. The module pins `trim-galore:2.5.0` and declares:

```groovy
tuple val(meta), path("*.html"), emit: html, optional: true
tuple val(meta), path("*.zip") , emit: zip , optional: true
```

These are **task-root globs and do not recurse**. Because §3.2 writes reports flat, they already
collect the raw set, and because §3.2's names are distinct, downstream consumers can tell the sets
apart. **So the only change the module strictly needs is the container bump** to a version that has
the flag.

Optional, not required: dedicated `html_raw` / `zip_raw` channels keyed on `*_raw_fastqc.*`, if
methylseq wants to route the sets separately rather than merely distinguish them.

*(Revision 1 put reports in a subdirectory, which these non-recursing globs would have collected
**not at all** — a silent drop, and strictly worse than the mixing revision 1 claimed to avoid.)*

### 7.2 MultiQC — no configuration needed

By **M1** and **M2** (§2), `sample_R1_raw`, `sample_R1` and `sample_R1_val_1` are three distinct
sample names, single-end and paired alike.

This is what the subdirectory could not deliver: nf-core's MULTIQC module stages inputs with
`path(multiqc_files, stageAs: "?/*")`, flattening them into numbered directories, so `fastqc_pre/`
would have been gone before MultiQC ran — defeating both `-d/--dirs` and `path_filters`. Filenames
survive staging; directories do not.

Worth stating for methylseq: had the names collided, the failure would have been a **silent drop** at
`log.debug`, not a visible merge.

### 7.3 nf-core/methylseq

Needs the module container bump and nothing else for correctness. `--skip_trimming` runs bypass TG
entirely and so get no raw QC from this route — unchanged by this revision.

### 7.4 Not affected

The `validation` CI job md5-compares against Perl 0.6.11 for SE, PE, hardtrim5, clock and demux.
`--fastqc_raw` is new and off by default, so no compared flag path changes (V8).

## 8. Assumptions

1. **`--fastqc_args` applies to both runs**, and still implies `--fastqc` (§3.1, last row).
2. **Sequential, not concurrent** (§3.5). Reversible; reversing means `resolved_threads` and the
   clumpify memory layout must account for the overlap.
3. **`_raw` is the marker**, inserted before the full extension. Fixed, not configurable, in v1.
4. **Symlinks are available.** Unix-only API; see §10 Open-2.
5. ~~fastqc-rust may create its report directory~~ — **verified it does not**
   (`runner.rs:477`, `archive.rs:66`). `create_dir_all` is mandatory, and V7 tests a real failure.
6. **`sequence_format` overrides extension dispatch** — verified at `config.rs:41` + `bam.rs:400-431`.
7. **MultiQC reads `Filename` from inside the report** — verified in `modules/fastqc/fastqc.py`.
8. **The `=1.1.1` pin does not move** for this work.
9. **A symlink's own name reaches the report; its target supplies the data** — **verified.**
   `FileGroup.name` is taken from `path.file_name()` of the path as given (`runner.rs:283-286`),
   *before* anything is opened; `open_group` passes the path through unchanged (`:356-360`); and
   nothing in the crate calls `canonicalize` or `read_link` on an input path — the one `canonicalize`
   (`runner.rs:321`) is on the report output directory, for a per-location write lock.
   `std::fs::metadata(path)` (`bam.rs:79`, `fastq.rs:416`) does follow the link, but only to read
   `.len()`, which should be the target's size.

## 9. Validation

Each states what would make it fail, because revision 1 shipped three that could not.

| # | What to verify | How | Expected / **how it fails** |
| --- | --- | --- | --- |
| V1 | Reports are produced, named and located | SE and PE with `--fastqc_raw` | `<stem>_raw_fastqc.{html,zip}` flat in `--output_dir`; no transient symlink left behind. **Fails if** the symlink survives or the name lacks `_raw` |
| V2 | **The in-zip `Filename` carries `_raw`** — the property the design rests on | Read `fastqc_data.txt` out of the zip, as `tests/integration_fastqc_phred64.rs` already does ("a wrong number in the report is invisible from outside") | `Filename	sample_R1_raw.fastq.gz`. **Fails if** FastQC is pointed at the input instead of the link — which existence checks cannot see |
| V3 | Format is declared, not sniffed | uBAM fixture **copied** (not hardlinked) to `reads.dat`, `--fastqc_raw` | Report shows BAM-derived stats. **Fails if** `sequence_format` is not plumbed: the file is read as FASTQ. Prove it can fail by passing `None` |
| V4 | A refused run creates nothing (§3.6) | Force an output collision with `--fastqc_raw` | Refusal, and **no `*_raw_fastqc.*` and no symlink**. **Fails if** the call is placed before the pre-flight |
| V5 | TG validates every input, not just the first (§3.7) | **Colorspace** R2 + valid R1, `--fastqc_raw` | TG's own colorspace message and **no R1 report**. **Fails if** the error comes from fastqc-rust, or R1 was reported on first. A *truncated* R2 cannot serve here — see §3.7 — and has its own test pinning the narrower scope |
| V6 | Specialty modes work **and write reports** | `--hardtrim5`, `--clock`, `--implicon` each with `--fastqc_raw` | The specific `*_raw_fastqc.zip` exists per mode. **Fails if** the call sits after the branch's `return` — "each completes" would not |
| V7 | The output directory is created | `--output_dir` at a fresh path | Created, reports inside. **Fails if** step 6c is skipped |
| V8 | Perl byte-identity intact | `cargo test`, then the `validation` CI job | **Exactly 883 tests plus the new ones** — a named number, since "883+" cannot fail on lost tests — and unchanged md5s |
| V9 | Symlink failure is loud | Output dir on a filesystem without symlink support, or a pre-existing path at the link name | Hard error naming the path; no partial report |
| V10 | `-o` in `--fastqc_args` is refused (§3.4) | `--fastqc_raw --fastqc_args "-o elsewhere"` | CLI refusal. **Fails if** reports land in `elsewhere/` |
| V11 | `--phred64` pre-flight widened (§5 step 4) | Phred+64 fixture, `--fastqc_raw` and **no** `--fastqc` | The pre-flight's own wording ("contains quality character", "within the first") **and no report written**. **Fails if** step 4 is skipped — note `fastqc.rs:85-95` has its own `--phred64` failure path, so asserting "the run fails" would pass either way |

## 10. Questions or ambiguities

**Open-1 — `--fastqc_args` implying `--fastqc` (§3.1).** Documented rather than fixed. Alternatives:
add `--fastqc_raw_args`, or narrow the implication when `--fastqc_raw` is present. Not blocking.

**Open-2 — symlink portability** (mechanism half now closed; see Assumption 9).
`std::os::unix::fs::symlink` is Unix-only and TG's Windows story is unclear. Fallbacks, worst-first:
hard error (v1), copy the input (doubles I/O), hardlink (see `feedback_hardlink_truncate_hazard` —
do not). Only portability is open; no spike is needed on whether the marker reaches the report.

**Open-3 — `--passthrough` carriers.** Assumption taken: **inputs only**, since a carrier is copied,
not sequencing data being QC'd. The `--phred64` pre-flight does iterate both.

**Open-4 — flag name.** Settled: `--fastqc_raw`, naming the data rather than the ordering, matching
how "raw reads QC" is said downstream.

**Open-5 — cross-file parallelism.** One `fastqc::run` per input forfeits fastqc-rust's own
multi-file threading. Deferred on the record; a `run_many` entry point would recover it.

## 11. Self-Review

**Checked.** Every source location re-verified after review corrected nine citations in revision 1.
The three mechanisms in §2 (M1/M2/M3) were each read in their own source rather than inferred, since
the entire design now rests on them. `main()`'s order re-traced *past* the insertion point, which is
where revision 1's gaps were concentrated.

**Adjusted.** Two design reversals and one structural fix; see §12.

**Remaining risks.**

1. **The design's load-bearing mechanism is now verified, not assumed** (Assumption 9). What remains
   of Open-2 is portability only: on a filesystem without symlinks the feature hard-errors rather
   than degrading. V9 covers that.
2. **Symlink lifecycle on the error path.** Step 6b must remove the link even when `fastqc::run`
   fails, or a failed run strands an artifact that the next run's collision pre-flight may trip over.
3. **§3.6 costs the single insertion point**, so the raw call now exists in ~7 places. Each is a site
   where a future refactor could drop it; V6 covers the specialty branches but not all seven.
4. **Three assumptions remain policy rather than fact**: 1, 2 and 3. All reversible.

## 12. Revision history

**Revision 2 (2026-11-18)** — after dual independent plan review (`PLAN_REVIEW_A.md`,
`PLAN_REVIEW_B.md`). Both reviewers agreed on five findings; no contradictions.

Two design reversals:

- **The `fastqc_pre/` subdirectory is gone**, replaced by a `_raw`-suffixed symlink and flat reports.
  The subdirectory was chosen to fix the MultiQC collision; it could not, because nf-core stages with
  `stageAs: "?/*"` and MultiQC keys on the report's internal `Filename`, not the path. An intermediate
  proposal — renaming the report file — fails for the same reason. Only changing the *analysed path's*
  filename reaches the `Filename` field. This also reduces the downstream requirement from
  "module bump + MultiQC config" to a container bump.
- **The uBAM refusal is gone**, replaced by declaring the format. `FastQCConfig.sequence_format`
  overrides extension dispatch, and TG already holds `input_formats`. This accepts `reads.dat` instead
  of refusing it and subsumes the guard's one-sidedness (`.sam`, `.fast5`, and FASTQ-misnamed-`.bam`
  were all unguarded).

One structural fix: **the single insertion point was an invariant breach.** Every
`preflight_output_collisions` runs at `:608` or later, so a raw pass at `:557-606` wrote reports and
read every input before any collision refusal, against the convention stated at `main.rs:368-371`.
The call moves into each dispatch branch, after its pre-flight.

Also corrected: `sanity_check_any` covers only `cli.input[0]` (§3.7, was claimed to cover all);
`specialty.rs:556` is test code, not production precedent (`#[cfg(test)]` at `:520`) — `fastq.rs:739`
cited instead; the `--phred64` gate also includes `|| output_format == UBam`; specialty dispatch is
`:606-722`; the paired path is `:781`/`:1048`; Assumption 5 upgraded from hedge to verified fact;
`-o` inside `--fastqc_args` added (§3.4); `--demux` added to §3.8; #379's FIFO refusal named in §6 as
the reason a second read is safe; and V1/V4/V11 rewritten because they could not fail on the mistakes
they existed to catch.

---

## 13. Implementation notes

Implemented 2026-11-18 on `feat/fastqc-raw`. **896 tests pass** (883 before, +13 new), `cargo fmt
--all -- --check` clean, `cargo clippy --all-targets --release -- -D warnings` clean.

### Deviations from §5, and why

1. **`SeqFormat` enum rather than `Option<&str>`** (§4). A sixth `Option<&str>` sitting next to
   `fastqc_args` would be silently transposable at twelve call sites. `fastqc::SeqFormat::{Fastq,Bam}`
   maps to the config string internally, so the two cannot be swapped.
2. **One `preflight_then_raw_qc` wrapper instead of thirteen `if cli.fastqc_raw` blocks** (§3.6).
   §11's risk 3 warned that seven-plus sites is somewhere a refactor drops one; routing every
   pre-flight through a wrapper makes "QC only after the refusal has passed" structural rather than
   conventional. **This broke a tripwire and the break was the point** — see below.
3. **`run_pretrim_fastqc(cli)` derives its own output directory** rather than taking `output_dir`.
   `run_specialty_paired` has no such parameter, and `output_dir` is always `cli.output_dir.as_deref()`.
4. **The format is re-detected rather than threading `input_formats` through every driver**
   (reviewer A's I-5 not taken). Threading it would change six driver signatures; re-detection is one
   block decompress against a full-file QC pass.
5. **A run-once guard on the raw pass.** Some drivers reach a pre-flight in more than one branch;
   an `AtomicBool` makes "at most once per process" true without auditing control flow in six
   functions.

### The tripwire, and what it caught

`tests/integration_preflight_tripwire.rs` asserts that `main.rs` has exactly thirteen pre-flight
sites **and that each is individually exercised**, the latter by parsing the `DUMP_PLANNED` dump,
which names its site via `#[track_caller]` on `preflight_output_collisions` ("so the dump can name
which dispatch arm built the list", `io.rs:178`).

Wrapping the pre-flight made all thirteen arms report the wrapper's line, collapsing the coverage
claim — `every_site_is_reached_by_some_case` failed with "pe and se both dispatched to
src/main.rs:179". Its sibling's message says plainly: *"do not just update the count."*

Fixed by marking the wrapper `#[track_caller]` too, so the location chains through to the real arm.
`every_site_is_reached_by_some_case` then passed unchanged, which is the proof that attribution
survived. The inventory test was updated to count the new call name, plus a **new assertion** that
`naming::preflight_output_collisions` appears exactly once in `main.rs` — a direct call would either
skip the raw pass or run it before the refusal, which is the whole invariant.

### Iteration log

- **#1** — `--implicon 8` in a CLI test parsed `8` as a positional input (the flag is
  `require_equals`), so validation failed on "got 3 input files". Changed to `--implicon=8`. Test
  error, not a product one.
- **#2** — `single_end_raw_report_names_the_sample_raw` asserted a 50 bp fixture; `illumina_10K.fastq.gz`
  is 150 bp. The `Filename` assertion — the one the design rests on — passed first time.
- **#3** — `a_bad_second_input_...` demanded that a truncated R2 be refused before R1 was reported
  on. Reproduced by hand: `sanity_check_any` passes on `truncated.fq.gz` (it reads the first record;
  the break is at bit 47,185,920), so R1's report is written and fastqc-rust then fails with
  "invalid DEFLATE data". **The test encoded an unachievable requirement**, and §3.7 had overstated
  the guard. Split into two tests — colorspace for "refused before any report", truncation for the
  narrower scope — and §3.7/§3.8/V5 corrected. Verified the three fixtures directly:
  colorspace → 0 reports, empty → 0 reports, truncated → 2 reports then failure.
- **#4** — the tripwire failure above.

### Not done

- `--fastqc_raw` is absent from `docs/src/content/docs/guide/outputs.md`'s FastQC section and from
  `wiki/usage/*`; only `guide/flags.md` and `wiki/reference/cli-flags.md` were updated.
- The CHANGELOG entry carries no PR link yet, since no PR exists. Add it when one is opened.
- §10's Open-1, Open-3 and Open-5 are unchanged and unaddressed, as planned.
- Portability (§10 Open-2) is untested: `std::os::unix::fs::symlink` is used directly, so this will
  not compile for a Windows target. No CI job covers one.

---

## 14. Code review round 1 — findings and fixes

Dual independent code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) plus a coverage audit
(`COVERAGE.md`, verdict **INCOMPLETE**). Both reviewers and the audit independently found the same
Critical. **904 tests pass** (896 before, +8), fmt and clippy clean.

### Fixed

| # | Finding | Fix |
| --- | --- | --- |
| **C1** | **Raw reports silently overwrite each other.** `run_pretrim_fastqc` pinned the directory to `output_dir.unwrap_or(".")` for the whole run. Two inputs sharing a basename across directories → distinct trimmed outputs (pre-flight passes), **one** raw report, exit 0. Reproduced SE and PE. | `io::raw_alias_name` now places the alias in `--output_dir` or beside its own input, like every other namer. Reports follow the alias. |
| **C2** | The raw reports were never pre-flight candidates — §3.6/§5 step 7 specified this and it was not built. | `io::raw_report_names` builds both artifacts as `PlannedOutput`s; `preflight_then_raw_qc` extends the candidate list with them. |
| **H1** | `--fastqc_raw` was outside the pre-flight tripwire entirely — two new writers per input on all 13 arms, none guarded. | `fastqc_active` knows the flag; three cases added, including `pe_trim_fastqc_raw_split_dirs`. |
| **H2** | `--fastqc_raw --fastqc_args X` refused on all four specialty modes, naming `--fastqc`. | `fastqc_requested()` narrowed (§3.1). |
| **M1** | The `_raw` alias can *be* another input — `s.fastq.gz` + `s_raw.fastq.gz` destroys input #2. §3.8 had dismissed this row as "non-colliding", which was wrong. | Checked against `guarded_inputs` before staging; refusal names the alias, the input it is for, and the clashing input. |
| **M2** | No unwind-safe cleanup: a panic inside fastqc-rust strands the alias and blocks every later run. | `AliasGuard`, a `Drop` guard. |
| **M3** | TOCTOU between the occupancy check and `symlink()`. | The check is gone; `symlink()`'s own `EEXIST` is mapped to the curated message. Atomic. |
| **M5** | `raw_alias_name` untested. | Four unit tests in `io.rs`, pinning extensions, the non-FASTQ fallback, the per-input directory, and both report names. |
| — | `exists()` missed a *dangling* symlink — the exact state a killed run leaves. | `symlink_metadata`, then subsumed by M3's atomic approach. |
| — | Two comments asserted pre-flight behaviour that does not exist, and `pub fn run` shipped undocumented because my insertion orphaned its rustdoc onto `resolved_threads`. | Corrected. |

### Deliberately not changed

**M4 — `--passthrough` carriers get no raw report.** §10 Open-3's decision stands: a carrier is copied,
not trimmed. Noted as an inconsistency with `guarded_inputs` and the `--phred64` pre-flight, which do
treat it as an input. Revisit if anyone asks for QC on a Multiome I1 read.

### The recurring failure, now four deep

Three validations in the plan could not fail on the mistake they named (caught at plan review), V7
could not either (caught by the audit), and `raw_and_trimmed_reports_coexist_with_distinct_sample_names`
asserted `assert_ne!` between a raw and a trimmed `Filename` — which differ under *any* scheme,
alias or not. Now pins both literal values.

**The fifth was nearly added while fixing the fourth.** The three new tripwire cases do not catch C1:
`raw_report_names` predicts candidates *through* `raw_alias_name`, so mutating that function moves
the prediction and the behaviour together and the tripwire stays green. Verified by mutation — the
tripwire passed with the bug restored. The guard that does discriminate is
`raw_reports_land_beside_each_input_not_pooled_in_the_cwd`, which spells the expected paths out
literally instead of deriving them; it fails (exit 101) under the same mutation. **A prediction
derived from the code under test cannot test that code** — the tripwire cases are still worth having
for the rest of the output set, but they are not C1's regression guard.

### V5 — MultiQC verified empirically (2026-11-18)

The last unproven link, and the only one no code review could close. Run on MultiQC **1.14.dev0**
(older than the 1.35 the reviewers read; the behaviour matched the source reading, so this point is
stable across both).

**Single-end — the collision case — passes.** Two reports found, two distinct samples:

| Sample | `Filename` |
| --- | --- |
| `illumina_10K` | `illumina_10K_trimmed.fq.gz` |
| `illumina_10K_raw` | `illumina_10K_raw.fastq.gz` |

For exactly the predicted reason: `_trimmed` **is** stripped — it is a terminal match in
`fn_clean_trim`, hence the bare `illumina_10K` — while `_raw` is not in the strip list at all. That
settles the half of §7.2 that was read from config rather than observed.

**Paired passes too.** Four reports, four distinct samples: `BS-seq_10K_R{1,2}_raw` and
`BS-seq_10K_R{1,2}_val_{1,2}`.

**The counterfactual is what makes this evidence.** Rewriting one line — the `Filename` field inside
the raw report's `fastqc_data.txt`, from `illumina_10K_raw.fastq.gz` back to `illumina_10K.fastq.gz`,
which is the alias's only effect — collapses it: MultiQC reports *"Found 1 reports"* and a single
sample survives. **The raw set is dropped silently**, with no warning on stderr. So the alias is
load-bearing rather than decorative, and the failure it averts is invisible rather than loud.

### Still open

Nothing blocking. The remaining items are §10's Open-1/3/5, M4 (§14), Windows portability, and the
two doc targets in §13's "Not done".
