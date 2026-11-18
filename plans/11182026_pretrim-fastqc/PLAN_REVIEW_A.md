# PLAN_REVIEW_A — `--fastqc_raw`

**Reviewer:** A (of two independent reviewers, no shared state)
**Plan:** `/Users/fkrueger/Github/TrimGalore/plans/11182026_pretrim-fastqc/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore` @ `dev` / `7055089`
**Date:** 2026-11-18

Verdict up front: the plan is unusually well sourced — I re-read every cited location and the
substance holds almost everywhere. It has two structural problems that matter more than anything in
its own Open-1..4 list: **the refusal in §3.3 is avoidable** (fastqc-rust has a format override the
plan did not find), and **the insertion point writes files ahead of the output-collision
pre-flights**, against a convention this repo states explicitly in a comment eight lines above the
proposed site. Separately, §7.1's premise about the nf-core globs is factually inverted, and §7.2 is
both understated and — I think — fixable from the Trim Galore side for the cost of one symlink,
which would retire the plan's own "necessary but not sufficient" caveat.

---

## 0. What I verified vs. inferred

### Verified by reading the source (not taken on trust)

| Claim | Status |
| --- | --- |
| `fastqc::run` signature, `src/fastqc.rs:55` | **Correct**, exact, parameter is named `output_path` |
| Twelve call sites; 7 in `main.rs` at `:1679 :1992 :1999 :2012 :2488 :2664 :2806` | **Correct**, all seven exact |
| 5 in `clump_only.rs` at `:306 :497 :498 :779 :1034` | **Correct**, all five exact |
| Every existing site passes an *output* path | **Correct** |
| `main()` at `main.rs:297`; `cli.validate()?` at `:330`; `sanity_check_any(&cli.input[0])?` at `:337` | **Correct** |
| `--phred64` pre-flight gated on `fastqc_requested()` at `main.rs:391` | **Correct**, exact |
| `let output_dir = cli.output_dir.as_deref();` at `main.rs:557` | **Correct**, exact |
| One insertion point between `output_dir` and specialty dispatch covers every mode | **Correct** — dispatch order is specialty `:606/:657/:668/:690` → `clump_only :722` → `run_ubam_output :741` → paired `:781` → SE `:783` |
| `#421` refusal block content, incl. the interleaved-uBAM arm | **Correct** (lines 1112–1151, see line-number note below) |
| `--fastqc` clap attribute at `cli.rs:345` | **Correct**, exact |
| `format.rs:329` — "classify by content (NOT by filename)" | **Correct**, verbatim |
| `main.rs:~2487` — "fastqc-rust dispatches on file EXTENSION (not content)" | **Correct** (comment at `:2484`, gate at `:2487`) |
| fastqc-rust `sequence/bam.rs:440` extension dispatch; `else` = "Everything else is assumed to be FASTQ", no error | **Correct**, verbatim, and the `else` genuinely does not error |
| fastqc-rust's only `create_dir_all` are `report/archive.rs:188,194` (extraction) | **Correct** (+ two `#[cfg(test)]` sites in `sequence/fastq.rs`) |
| nf-core `trimgalore` globs are `path("*.html")` / `path("*.zip")`, untargeted | **Correct** |
| Module pins `trim-galore:2.5.0` | **Correct** (`community.wave.seqera.io/library/trim-galore:2.5.0--e26c3be09dac93c9`) |
| MultiQC `fn_clean_exts` contains no `_val` / `_trimmed` | **Correct** |
| MultiQC `fn_clean_trim` contains both `_val` and `_trimmed` (and `.trimmed`) | **Correct** |
| **Paired survives, single-end collides** | **Correct** — I traced the whole chain, see §4 |
| `test_files/` has uBAM fixtures for V3 | **Correct** (`ubam_test.bam`, `ubam_paired_test.bam`, 16 others) |
| `wiki/reference/cli-flags.md` exists | **Correct** |
| Test count ≈ 883 | 884 `#[test]` attributes across `src/` + `tests/` today |

### Assumption 5 — upgraded from "flagged" to **verified**

The plan hedges ("confirm in V6 rather than assuming"). It can stop hedging. fastqc-rust writes the
HTML with `std::fs::write(&html_path, …)` at `runner.rs:477` and the zip with
`fs::File::create(zip_path)` at `report/archive.rs:66`, and `report_location` (`runner.rs:299-311`)
only *joins* — nothing between `runner::run`'s entry and those writes creates a directory. TG calls
`runner::run` directly (`fastqc.rs:78`), bypassing fastqc-rust's own `main.rs`, so any directory
setup there is irrelevant. **`create_dir_all` in step 6 is mandatory, not defensive**, and V6 is
testing a real failure rather than a hypothetical.

### Mechanism I verified that the plan asserts without showing

Report naming is locked to the input's basename in two places at once, and this is load-bearing for
both Open-3 and §7.2:

- `build_file_groups` (`runner.rs:283-286`) sets `group.name = path.file_name()` — the basename, not
  the path.
- `report_location` (`runner.rs:300`) → `strip_extensions(group.name)` (`runner.rs:739-750`), which
  strips `.gz .bz2 .txt .fastq .fq .csfastq .sam .bam .ubam .fast5`, **each once, in that order**.
- `file_display_name = group.name` (`runner.rs:387`) → `module.set_filename(…)` (`runner.rs:402`) →
  the `Filename` line inside `fastqc_data.txt`.

So §3.2's naming claim is exactly right, and Open-3's collision is exactly real.

### Could not verify

- Nothing. Every external fact was checkable locally: the vendored crate is at
  `~/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/fastqc-rust-1.1.1/`, MultiQC **1.35** is
  installed on this machine, and the nf-core module fetched cleanly.

### Wrong

1. **§7.1's premise is inverted** — see C-2 below. This is the only substantive factual error.
2. **"paired path from ~`:982`"** (§2) — `run_paired_trim` is defined at `main.rs:1048` and called at
   `:781`. Line `:982` is inside `run_clump_only`.
3. **§7.2's stated rule for `fn_clean_trim` is not the rule** — it strips a leading match as well as
   a terminal one. The paired/SE *conclusion* is unaffected; the rule as stated is incomplete.
4. **Line-number drift of 3–4 lines** in four citations: `fastqc_requested()` is at `cli.rs:613` (the
   doc comment opens at 605); the `#421` block opens at `cli.rs:1112`; specialty dispatch opens at
   `main.rs:605/606`; `main.rs:2487` is the gate, the quoted comment is at `:2484`. Immaterial
   individually, but the plan presents these as read-from-source, so they should be right.

---

## 1. Logic review

### C-1 (Critical) — the pre-trim run writes files ahead of the output-collision pre-flights

§5 step 7 puts the call "after `output_dir` is bound (`:557`) and before the specialty dispatch
(`:609`)". Every output-collision pre-flight runs *after* that window:
`naming::preflight_output_collisions` at `:608` (hardtrim5), `:641` (hardtrim3), inside
`run_specialty_paired` for clock/implicon, and the three trim-path pre-flights around `:492`.

The repo states the opposite convention in a comment eight lines above the proposed site
(`main.rs:368-371`):

> Placed ahead of every dispatch branch — including the early-returning specialty and clump-only
> paths — and ahead of `ensure_output_dir`, so a rejected run creates nothing.

Under the plan as written, a run that `preflight_output_collisions` is about to refuse first writes
a populated `fastqc_pre/` and burns a full decompress+scan of every input. CLAUDE.md calls the
collision pre-flight load-bearing ("Don't bypass it"); writing files before it is a soft bypass.

**Fix:** either move the pre-trim call to after the pre-flight for each branch (costs the "one site
serves every mode" elegance), or — better — hoist the pre-flights, which already only need
`output_dir` and `gzip`, both bound by `:557`. The second keeps the plan's single insertion point.
Note the plan already leans on exactly this ordering principle in §3.6 ("Already refused by
`sanity_check_any` at `:337`, before the new call site"), so it is applying the rule in one row of
the table and breaking it in step 7.

### C-2 (Critical, factual) — §7.1 describes a counterfactual as the current state

> Both globs are **untargeted**, so today they would sweep up pre- and post-trim reports into the
> same two channels with nothing to tell them apart.

`path("*.html")` is a glob at the **task root**. It does not recurse. With the plan's design, raw
reports in `fastqc_pre/` are **not collected at all** by the current module — they are silently
dropped, not mixed. The sweep-into-one-channel hazard is what would happen if TG wrote raw reports
*flat*, which is the design the plan rejects.

This matters because it changes the failure mode a methylseq user hits with an un-bumped module from
"two kinds of report in one channel" (confusing, visible) to "pre-trim QC requested and nothing
appears" (silent). Both argue for the subdirectory; the paragraph needs restating so the next reader
does not inherit a wrong model of the glob.

### I-1 (Important) — the §3.3 refusal is avoidable; fastqc-rust has a format override

`FastQCConfig` has `sequence_format: Option<String>` (`config.rs:41`), and
`open_sequence_file` checks it **before** extension dispatch (`sequence/bam.rs:401-431`), accepting
`"bam"`, `"sam"`, `"bam_mapped"`, `"sam_mapped"`, `"fastq"`.

Trim Galore already knows the answer: `input_formats: Vec<InputFormat>` is computed for every
`cli.input` at `main.rs:343-347`. Passing `sequence_format: Some("bam")` for a content-detected uBAM
makes **any** filename work, including `reads.dat`. The cost is one extra parameter on `fastqc::run`
(or a sibling entry point so the existing twelve sites stay untouched).

The plan's framing — "refuse where a sound report cannot be produced" — is sound policy, but here a
sound report *can* be produced. Accepting beats refusing, and §3.3's "the reason the guard is not
optional" is only true given the current `fastqc::run` signature, which this plan is free to change.
The refusal also costs users the one input shape that content-detection exists to support.

I'd still keep a guard for the residual case (detected format that fastqc-rust has no
`sequence_format` word for), but it stops being the headline of the feature.

### I-2 (Important) — the guard as specified is one-sided

fastqc-rust dispatches four suffixes, not two: `.bam`/`.ubam` → BAM, `.sam` → SAM, `.fast5` →
Fast5 (`bam.rs:440-450`). §3.3 guards only *uBAM content without a `.bam`/`.ubam` suffix*. The mirror
cases are unguarded:

- A content-detected **FASTQ** named `reads.bam` or `reads.sam` → read as BAM/SAM.
- A content-detected FASTQ named `reads.fast5` → routed to the HDF5 reader.

This is not purely hypothetical: CLAUDE.md states that `bgzip x.fq` and `samtools import` output
both look BGZF-framed, which is precisely why TG detects by content. A user who has a BGZF FASTQ and
named it badly is in scope.

**Fix:** state the invariant as *the detected format must agree with the suffix fastqc-rust will
dispatch on*, and check it in one direction only if I-1 is declined. With I-1 adopted, the whole
class disappears — set `sequence_format` from `input_formats[i]` and extension never matters.

### I-3 (Important) — §3.1's behaviour table is violated by `--fastqc_args`

§3.1 promises "`--fastqc_raw` alone → pre-trim reports only". But `--fastqc_args` makes
`fastqc_requested()` true (`cli.rs:613`), so `--fastqc_raw --fastqc_args "--nogroup"` in a normal
(non-specialty) mode produces **both** sets. §8 Assumption 1 discusses only the specialty-mode
*refusal* consequence; it does not mention that in the common case the user silently gets post-trim
reports they did not ask for.

Either the table needs a third row, or `--fastqc_args` needs to stop implying `--fastqc` when
`--fastqc_raw` is the only QC flag given (Open-1's second option — which this finding makes less
"magic" and more "necessary for the documented behaviour to be true").

### I-4 (Important) — `-o` inside `--fastqc_args` silently relocates the raw reports

`fastqc::run` sets `config.output_dir` from its parameter, **then** calls `apply_fastqc_args`
(`fastqc.rs:73-75`), and `apply_fastqc_args` handles `-o`/`--outdir` by overwriting
`config.output_dir` (`fastqc.rs:124-129`).

So `--fastqc_raw --fastqc_args "-o qc"` writes the raw reports to `qc/`, not `fastqc_pre/`, with no
warning. §3.2 asserts the location unconditionally and §7.1 builds a glob contract on it. The
contract is not enforceable as specified.

This is pre-existing behaviour on the post-trim path, where it is harmless (the user asked for it).
It becomes a correctness issue the moment a downstream module globs a fixed directory name.

**Fix:** refuse `-o`/`--outdir` in `--fastqc_args` when `--fastqc_raw` is given, or apply the args
first and let the TG-chosen directory win for the raw run. Refusing is more honest and one `bail!`.

### I-5 (Important) — step 5 re-detects a format already in hand

`main.rs:343-347` computes `input_formats` for every `cli.input` under the comment "Detect format of
all inputs once; the call is cheap". Step 5 proposes a fresh `detect_input_format(input)?` per input.
Besides the redundant I/O, two detections of the same path can in principle disagree (the file can
change between them), and the rest of `main()` acts on `input_formats`. Reuse it. `--passthrough`
genuinely needs its own call (as `main.rs:465` already does) if Open-2 resolves to "yes".

### O-1 (Optional) — `--demux` is missing from §3.6

Demultiplexing dispatches inside the SE trim path (`main.rs:1689-1695`), so `--fastqc_raw` yields one
raw report of the pre-demux input. That is the right behaviour, but `--demux` is one of the five
flag paths in the load-bearing validation matrix, and §3.6 lists `--clump_only` and the four
specialty modes while omitting it. One row.

### O-2 (Optional) — `fastqc::run`'s parameter name becomes actively wrong

The plan notes the parameter is "even *named* `output_path`" as evidence. Once a call site passes an
input, the name misleads at five of the thirteen reading positions. Rename to `target` (or add the
sibling entry point from I-1 and leave `run` alone). Pure readability; CI will catch nothing.

### O-3 (Optional) — opaque error if the directory is missing

`fastqc::run` maps `runner::run`'s `Err(code)` to `"FastQC analysis failed (exit code {code})"` with
no path. If `create_dir_all` is skipped or fails, that is all the user sees. Give step 6's
`create_dir_all` a `.with_context(|| …)` naming the directory.

### Not a problem (verified, worth recording)

**The extra full read of each input is safe, and the plan does not say why.** §6 adds a second (and
with adapter auto-detection, a third) read of every input. That is only safe because issue #379
already rejects non-restartable input — FIFOs, pipes, `/dev/stdin`, process substitution — at
`Cli::validate` (`cli.rs:536` via `check_restartable_input`) and again inside `detect_input_format`
(`format.rs:342`), with the message "cannot be re-read from the start". I went looking for this as a
likely silent hole and it is already closed. The plan should cite #379 in §6 so the next person
widening the input surface knows which guard the pre-trim read is leaning on.

---

## 2. Assumptions

| # | Plan's assumption | Assessment |
| --- | --- | --- |
| 1 | `--fastqc_args` applies to both runs | Holds, but the stated consequence is incomplete — see **I-3** (both report sets in normal modes) and **I-4** (`-o` moves the directory). Both are more user-visible than the specialty-refusal wart the plan flags. |
| 2 | Sequential, not concurrent | Sound, and the plan is right that it retires the thread-budget work. Verified: `resolve_clump_layout` gates `resolved_threads` on `fastqc_requested()` (`main.rs:248-250`), untouched by a sequential pre-trim run. |
| 3 | `fastqc_pre/` fixed, not configurable | Right call for a glob contract — but **I-4** means it is not actually fixed. |
| 4 | FastQC's own naming inside that directory | Verified, and structurally unavoidable (both the report filename and the in-report `Filename` derive from `group.name` = the input basename). This is the hinge of §7.2 — see §4. |
| 5 | fastqc-rust does not create its report directory | **Verified, not assumed.** Drop the hedge. |
| 6 | `.bam`/`.ubam` only | Verified — but incomplete: `.sam` and `.fast5` dispatch too (**I-2**). |
| 7 | `=1.1.1` pin does not move | Fine. |

### Unstated assumptions I surfaced

1. **Inputs are re-readable.** True today via #379; see "Not a problem" above. Make it explicit.
2. **The guard must be a refusal.** False — fastqc-rust accepts a format override (**I-1**).
3. **The fixed directory is actually fixed.** False (**I-4**).
4. **`fastqc_pre/` is relative to `--output_dir`, which `ensure_output_dir` has already created.**
   True only if the call lands after `main.rs:561`. §5 step 7 says "after `output_dir` is bound
   (`:557`)", which admits a placement at 558 where `--output_dir` itself does not yet exist. A
   `create_dir_all` on the nested path papers over it, but the plan should name `:561`, not `:557`,
   as the floor — and C-1 pushes it later still.
5. **Nothing downstream globs `*.html` recursively.** True for the nf-core module, and it is why
   C-2's correction matters.

---

## 3. Efficiency

§6 is directionally right: the trimming hot path is untouched, no new allocation in it, and §3.5's
sequencing genuinely removes the concurrency questions. Three refinements:

- **Quantify it.** "Adds wall time" understates. The pre-trim run is ≈ the cost of one post-trim
  FastQC per input — full decompress plus every module's per-read work. On the `test_files/` fixtures
  FastQC dominates the run; on real libraries it is a second full pass. For methylseq the plan's
  "cost-neutral against the status quo" framing is right (it replaces a separate FASTQC process), but
  for a bare `trim_galore --fastqc_raw` user it roughly doubles a `--fastqc` run's wall time. Say
  both.
- **Peak memory is unchanged only because of the sequencing.** `resolve_clump_layout` reserves for
  the clump buffer and `fastqc_requested()` gates its FastQC term; a sequential pre-trim run holds
  neither at the same time. One clause, not a paragraph.
- **(I-6, Important-for-scale) One `fastqc::run` per input forfeits fastqc-rust's file-level
  parallelism.** `runner::run` calls `plan_threads(config.threads, file_groups.len(),
  available_parallelism())` and builds an outer rayon pool of `outer_slots` — one slot per file group
  analysed in parallel, each driving `threads_per_file` inner threads. Called once per file, that
  machinery can never engage across inputs. Passing all inputs to one `runner::run` parallelises
  across files inside the same `--cores` budget, and also engages `report_locks` (`runner.rs:313-320`),
  which exists precisely for "`run1/S1.fastq.gz` and `run2/S1.fastq.gz` into one `-o`" — i.e. Open-3's
  case. For a 2-file PE run the gain is modest; for a multi-pair run (which §3.6 explicitly supports)
  it is the difference between N sequential FastQC runs and `min(N, cores)` concurrent ones.
  Requires a `fastqc::run_many`-shaped entry point. Worth doing; at minimum, worth recording as a
  deliberate deferral rather than an oversight.

Note that even with one call, `report_locks` only prevents *corruption*, not overwrite — Open-3's
last-writer-wins stands either way.

---

## 4. Validation sufficiency

I went hunting specifically for validations that cannot fail on the mistake they exist to catch.
Found two, plus three uncovered findings.

### V3 is the best-constructed validation in the plan

It has an explicit negative control — "Then prove it can fail: the same file as `reads.bam` must
succeed". That is the right instinct and the rest of the table should copy it. Two gaps:
- It does not cover the mirror case (**I-2**): a FASTQ named `*.bam`.
- "no report written and no trimming done" is prose. Make it an assertion that `fastqc_pre/` does
  **not exist** after the refusal, otherwise the test passes on a refusal that happens *after* the
  reports are written.

### V4 is largely satisfied by a no-op — **this is the dangerous one**

> V4 | Specialty modes accept it | `--hardtrim5`, `--clock`, `--implicon`, each with `--fastqc_raw` |
> Each completes and writes raw reports

The most likely implementation mistake here is placing the call *after* the specialty dispatch —
natural, because all twelve existing FastQC call sites are post-output, and the specialty branches
all `return Ok(())`. Under that mistake every specialty mode still **completes** cleanly; only the
reports are missing. Half the stated pass criterion ("each completes") is satisfied by the bug.

**Fix:** the assertion must be *the report files exist*, per mode, named explicitly —
`fastqc_pre/<input-stem>_fastqc.html` and `.zip` — not "completes and writes raw reports". And add a
negative control in V3's style: with the call deliberately moved after dispatch, V4 must fail.

### V6 can fail on its mistake — confirmed

I checked this one specifically because of Assumption 5's hedge. `ensure_output_dir`
(`main.rs:561` → `io.rs:273`) creates `--output_dir` but not `<output_dir>/fastqc_pre`, and
fastqc-rust creates nothing, so a missing `create_dir_all` produces a real failure. V6 is sound.

Gap: it only covers the `--output_dir`-**set** case. §3.6 promises `fastqc_pre/` in the working
directory when `--output_dir` is unset, and nothing tests that. Add the variant — it is the path
where `ensure_output_dir` does nothing at all, so it is strictly the weaker-covered of the two.

### V7's pass criterion cannot detect a regression

"883+ tests green" passes if tests are *lost*. The count today is 884 `#[test]` attributes, and
step 8 adds roughly six. Name the expected number, or drop the count and keep "all green" — the
`+` makes the assertion unfalsifiable in the direction that matters.

### V8 is thin

"Pre-flight fires as it does for `--fastqc`" can pass on any error. State the expected message
fragment (`contains("--phred64 declares Phred+64 input")`), so it cannot pass on an unrelated bail.
Also note the uBAM arm is unreachable: `--phred64` + uBAM is already refused at `main.rs:372`, and the
pre-flight skips uBAM at `:400`, so V8 is a FASTQ-only test by construction. Fine, but say so.

### Missing validations

| Gap | Finding |
| --- | --- |
| No run that trips `preflight_output_collisions` with `--fastqc_raw` and asserts `fastqc_pre/` was **not** created | **C-1** |
| No test of `--fastqc_raw --fastqc_args "-o qc"` | **I-4** — the documented directory is silently overridable |
| No test that `--fastqc_raw --fastqc_args X` in a *normal* mode does/doesn't also produce post-trim reports | **I-3** — §3.1's table is currently untested in its interesting case |
| No `--demux` case | **O-1** |
| No `--output_dir`-unset case | V6 gap above |
| No mirror-suffix case (FASTQ named `*.bam`) | **I-2** |

### V5 — I resolved most of it statically; it is also understated

The plan calls V5 its weakest evidence and asks for an empirical run. I traced MultiQC **1.35**'s
actual code path, which answers it without running anything:

- `base_module._clean_s_name` applies `fn_clean_exts` **first** as a truncate/remove/regex loop
  (`base_module.py:1030-1065`), **then** `fn_clean_trim` (`:1066-1070`) — and that loop strips a
  match at **either** end, each list entry applied once, in list order.
- Paired: `sample_R1_val_1.fq.gz` → truncate at `.fq` → `sample_R1_val_1`. `_val` is not terminal
  (it ends `_1`), nothing else matches → `sample_R1_val_1`. Raw → `sample_R1`. **Distinct. Plan
  correct.**
- SE: `sample_trimmed.fq.gz` → `sample_trimmed` → `_trimmed` *is* terminal → `sample`. Raw →
  `sample`. **Collide. Plan correct.**

Two things the plan does not have:

1. **The SE failure is a silent *drop*, not a merge, and it is logged at `debug`.** MultiQC's FastQC
   module computes an s_name from the zip filename and then `continue`s if that name is already
   present — `log.debug(f"Skipping '{f['fn']}' as already parsed '{s_name}'")`
   (`modules/fastqc/fastqc.py:225-229`) — and `parse_fastqc_report` logs its duplicate case at
   `log.debug` too (`:349-350`). So for single-end, one of the two reports is **never parsed**, which
   one depends on file-discovery order, and nothing appears at default verbosity. §7.2's "collides"
   is too kind. This strengthens §3.2's refusal to call the flag MultiQC-ready — keep that sentence.
2. **`-d`/`--dirs` will not work the way §7.2 hopes, under Nextflow.** `prepend_dirs`
   (`base_module.py:1017-1026`) prefixes `Path(root).parts`. But `path("fastqc_pre/*.zip")` emits the
   *files*; Nextflow stages them **flat** into the MultiQC task directory, so there is no directory
   left to prepend. The workable routes are emitting the directory itself (`path("fastqc_pre")`),
   `--replace-names`, or a module-scoped `fn_clean_exts` entry in the pipeline's `multiqc_config`.
   Fix §7.2 before methylseq plans around `-d`.

The module's own docs confirm the plan's §3.2 caveat verbatim: *"Sample names are discovered by
parsing the line beginning `Filename` in `fastqc_data.txt`, **not** based on the FastQC report
names."* (`modules/fastqc/fastqc.py:86-88`).

Keep V5 as an empirical confirmation, but it is now confirming a traced mechanism rather than
supplying the only evidence. And it does not gate Trim Galore.

---

## 5. Alternatives

### A-1 — Feed FastQC a renamed symlink, and §7.2 stops being a dependency (**recommended, high value**)

This is the finding I would most like the plan to consider, because it retires its own closing
caveat ("necessary but not sufficient").

The MultiQC sample name comes from `Filename` inside `fastqc_data.txt`, which is
`group.name = path.file_name()` of whatever file FastQC opened. The report *filename* comes from
`strip_extensions` of the same string. So renaming the report cannot help — but renaming the **input
FastQC sees** changes both, together.

TG can create `fastqc_pre/sample_raw.fastq.gz -> ../sample.fastq.gz` and run FastQC on the symlink.
`File::open` follows it; `file_name()` is the link name. Result:
`fastqc_pre/sample_raw_fastqc.{html,zip}`, with `Filename  sample_raw.fastq.gz` inside. MultiQC then
cleans that to `sample_raw` — distinct from `sample` *and* from `sample_trimmed`, in both SE and PE,
with **no** `-d`, no `--replace-names`, no second FastQC section, no methylseq MultiQC change at all.
§7.1's glob bump remains the only downstream work.

Cost: one symlink per input, removed after the run. Use `ln -s` semantics
(`std::os::unix::fs::symlink`), never a hard link — a hard link aliased to a fixture is the hazard
recorded in `feedback_hardlink_truncate_hazard`. Windows needs a fallback (copy, or skip the rename
and accept §7.2); TG's `#![cfg(unix)]` precedent in `tests/integration_non_restartable_input.rs`
shows the project is comfortable with unix-gated behaviour.

Trade-off: a symlink in the output tree is a small surprise, and `_raw` becomes a user-visible naming
convention rather than an internal detail. Against that, it converts the feature from "produces
reports methylseq cannot consume without extra config" to "works". Worth a spike before committing to
the plan's current shape.

### A-2 — `sequence_format` instead of a refusal (see I-1) — **recommended**

Accept the content-detected uBAM whatever it is called. One parameter, and §3.3 shrinks from the
feature's headline risk to a footnote. Also subsumes I-2 entirely.

### A-3 — One `runner::run` call with all inputs (see I-6)

Parallelism across inputs inside the existing `--cores` budget, plus `report_locks`. Needs a new
entry point. Defer explicitly rather than by omission.

### A-4 — Dead end worth closing off in the plan

"Give raw reports a distinguishing suffix (`sample_raw_fastqc.html`) to fix the MultiQC clash"
**does not work**, and someone will propose it. The MultiQC sample name is read from `Filename`
inside the report, which TG cannot set independently of the file FastQC opened (`runner.rs:387,402`).
Renaming the output has no effect on the sample name. One sentence in §7.2 saves a future round.
A-1 is the version of this idea that does work, because it changes the input name, not the output's.

### A-5 — Flag name (Open-4)

One data point: nf-core and MultiQC both already use "raw" for pre-trim FastQC (methylseq's FASTQC
process on raw reads vs. TRIMGALORE). `--fastqc_raw` matches the downstream vocabulary; `--fastqc_pre`
would introduce a second word for the same thing. Weak, but it is the only non-aesthetic argument
available and it points one way.

---

## 6. Action items

### Critical

1. **Move the pre-trim call after the output-collision pre-flights** (or hoist the pre-flights above
   it). As written, a refused run leaves `fastqc_pre/` populated and has read every input — against
   the convention stated at `main.rs:368-371`. §5 step 7 and §3.6. **(C-1)**
2. **Rewrite §7.1's first paragraph.** `path("*.html")` does not recurse; with `fastqc_pre/` the
   un-bumped module collects *nothing*, it does not mix the two sets. The conclusion (bump the
   module) is unchanged; the stated mechanism is wrong. **(C-2)**
3. **Add a validation that a pre-flight-refused run creates no `fastqc_pre/`.** Nothing in §9 covers
   C-1.
4. **Fix V4 so it cannot pass on the likeliest bug.** Assert the specific report files exist per
   specialty mode; "each completes" is satisfied by a call placed after the dispatch. **(§4)**

### Important

5. **Reconsider §3.3 against `FastQCConfig.sequence_format`** (`config.rs:41`, `bam.rs:401-431`). TG
   already has `input_formats`; passing the format explicitly accepts what the plan refuses, for one
   parameter. **(I-1)**
6. **State the suffix guard as format/suffix *agreement*,** covering `.sam` and `.fast5` and the
   FASTQ-named-`.bam` mirror — or adopt item 5, which removes the class. **(I-2)**
7. **Handle `-o`/`--outdir` in `--fastqc_args`.** It overwrites `config.output_dir` after TG sets it
   (`fastqc.rs:124-129`), so `fastqc_pre/` is not actually fixed and §7.1's glob contract is not
   enforceable. Refuse it with `--fastqc_raw`, or apply args before the TG-chosen directory. Add a
   test. **(I-4)**
8. **Close §3.1's gap for `--fastqc_raw --fastqc_args X`** — in normal modes it silently produces both
   report sets. Either document the third row or narrow the implication (Open-1, option 2). Add a
   test. **(I-3)**
9. **Reuse `input_formats` (`main.rs:343-347`) in step 5** instead of re-detecting. **(I-5)**
10. **Correct §7.2 on two points:** the SE case is a silent *drop* logged at `debug`
    (`modules/fastqc/fastqc.py:225-229, 349-350`), not a merge; and `-d` will not help once Nextflow
    stages the files flat. **(§4)**
11. **Spike A-1 (renamed symlink) before freezing the design.** If it works it removes the plan's
    "necessary but not sufficient" caveat entirely and leaves §7.1 as the only downstream work.
12. **Fix V7's pass criterion** — "883+" cannot fail on lost tests. Name the number.
13. **Drop Assumption 5's hedge; it is verified.** `runner.rs:477` + `archive.rs:66` write without
    creating. `create_dir_all` is mandatory.
14. **Decide I-6 explicitly** — one `fastqc::run` per input forfeits fastqc-rust's cross-file
    parallelism. Defer on the record, or add a `run_many` entry point.

### Optional

15. Add a `--demux` row to §3.6 (dispatches at `main.rs:1689-1695`; one pre-demux report). **(O-1)**
16. Cite #379 in §6 as the reason a second full read of each input is safe
    (`cli.rs:536`, `format.rs:342`). It is currently an unstated load-bearing assumption. **(§1)**
17. Add the `--output_dir`-unset variant to V6, and the expected message fragment to V8.
18. Rename `fastqc::run`'s `output_path` parameter, or add the sibling entry point from item 5.
    **(O-2)**
19. `.with_context()` on step 6's `create_dir_all` — fastqc-rust's failure is a bare exit code.
    **(O-3)**
20. Add the A-4 note to §7.2 so nobody proposes report-renaming as a MultiQC fix.
21. Fix the drifted citations: `fastqc_requested()` → `cli.rs:613`; `#421` block → `cli.rs:1112`;
    specialty dispatch → `main.rs:605`; paired path → `main.rs:781` (def at `:1048`), not `:982`;
    `main.rs:2487` is the gate, `:2484` the comment.
22. §5 step 9: the docs tree is an Astro site — flag reference lives under
    `docs/src/content/docs/guide/flags.md`, not a bare `docs/`. Also update `--fastqc`'s own clap doc
    (`cli.rs:338-344`), which currently tells specialty-mode users to "run FastQC separately".

---

## 7. Note on this plan's quality

The self-review in §11 is accurate about what it checked, and the two inversions it records (uBAM,
MultiQC) are both real and both correctly resolved. The citation accuracy is high — twelve call
sites, all exact. The failure mode this repo has seen before, *a validation that cannot fail on the
mistake it exists to catch*, appears once (V4) and is avoided deliberately in V3, which carries an
explicit negative control. That is the right habit; apply it to V4, V6 and V8.

The two structural items (C-1, I-1) are both cases of the plan reasoning correctly from a constraint
it did not question: that the call must sit at the single point every mode passes through, and that
`fastqc::run`'s signature is fixed. Neither constraint is real.
