# Plan Review B — `--fastqc_raw`

**Plan:** `/Users/fkrueger/Github/TrimGalore/plans/11182026_pretrim-fastqc/PLAN.md`
**Reviewer:** B (independent, fresh context; no shared state with Reviewer A)
**Repo state:** `dev` @ `7055089`
**Date:** 2026-11-18

---

## 0. Verification ledger

The plan asks to be checked rather than trusted, so this comes first. Everything below is read from
the actual file at the actual line, not inferred.

### Verified correct

| Plan claim | Verdict |
| --- | --- |
| `fastqc::run` signature, `src/fastqc.rs:55` | ✅ exact — `pub fn run(` is line 55, params and order as quoted |
| **Twelve** `fastqc::run` call sites: `main.rs` 1679/1992/1999/2012/2488/2664/2806, `clump_only.rs` 306/497/498/779/1034 | ✅ all twelve, all at the stated lines |
| Every existing call site passes an output path | ✅ |
| `Cli::fastqc_requested()` at `cli.rs:609`, body `self.fastqc \|\| self.fastqc_args.is_some()` at `:610` | ✅ |
| `--fastqc` clap attribute at `cli.rs:345` | ✅ |
| `#421` refusal block `cli.rs:1113-1151` | ✅ comment 1112-1113, `if self.fastqc_requested()` 1114, last bail ends 1150. Includes the interleaved-uBAM arm as described |
| `#421` acceptance-twin test block at `cli.rs:2773` | ✅ exact comment line |
| `main()` at `main.rs:297`; `cli.validate()?` at `:330`; `sanity_check_any(&cli.input[0])?` at `:337` | ✅ |
| `let output_dir = cli.output_dir.as_deref();` at `main.rs:557` | ✅ exact |
| `run_single_end_trim` call at `main.rs:783` | ✅ |
| **No early return between `:330` and `:606`** — so one insertion point really does serve every mode | ✅ verified by scanning the whole span; the plan's central structural claim holds |
| `format.rs:329` says *"Peek the first bytes of `path` and classify by content (NOT by filename)"* | ✅ verbatim, line 329 |
| `fastqc-rust` `sequence/bam.rs:440`: `if name_lower.ends_with(".bam") \|\| name_lower.ends_with(".ubam")`, final `else` = *"Everything else is assumed to be FASTQ"*, no error | ✅ exact line, exact comment |
| Assumption 5 — `fastqc-rust` does **not** create its configured report dir | ✅ **confirmed, not just plausible.** `create_dir_all` appears only in `report/archive.rs:188,194` (zip *extraction*) and in `sequence/fastq.rs` test code. `runner.rs:452-458` does `output_dir.join(...)` and writes straight to it. TG must create it. |
| Assumption 4 — raw reports named after the file read, `--basename` irrelevant | ✅ `runner.rs:299` `report_location` derives `base_name` from `strip_extensions(group.name)`, and `group.name` is `path.file_name()` (`runner.rs:284-290`) |
| nf-core/modules `trimgalore` emits untargeted `path("*.html")` / `path("*.zip")` | ✅ verified against **`origin/master` @ `cb523f66bc`**, lines 21-22 |
| The module pins `trim-galore:2.5.0` | ✅ verified upstream: `community.wave.seqera.io/library/trim-galore:2.5.0--e26c3be09dac93c9`. (My local `~/Github/modules` checkout was stale at 2.3.0 — the plan is right, my working copy was wrong.) |
| MultiQC `fn_clean_exts` has no `_val` / `_trimmed`; `fn_clean_trim` has both | ✅ checked in **MultiQC 1.35** (`config_defaults.yaml:413,415,416` are the `fn_clean_trim` entries) and in 1.14 — identical |
| **§7.2's whole paired-safe / SE-collides conclusion** | ✅ **verified from MultiQC source, not just config** — see §2.1 below. The plan calls this its weakest evidence; it is stronger than it thinks. |
| methylseq has `params.skip_trimming` | ✅ `nextflow.config:90`, used at `workflows/methylseq/main.nf:102,371` |
| `--phred64` pre-flight iterates `cli.input.iter().chain(cli.passthrough.iter())` (§10 Open-2) | ✅ `main.rs:395` |
| `--fastqc_raw` / `fastqc_pre` are free names | ✅ zero hits in `src/`, `docs/`, `wiki/` |
| "883+ tests" (V7) | ⚠️ consistent — 884 `#[test]` attributes across `src/` + `tests/`. Not run. |

### Verified wrong

| Plan claim | Reality |
| --- | --- |
| §3.2: *"several TG sites already do this, e.g. `specialty.rs:556`"* (create their own report dir) | ❌ **`specialty.rs:556` is inside a unit test** (`test_clock_happy_path_extracts_umi_and_clips`, `std::fs::create_dir_all(&dir)` on a temp dir). `specialty.rs` has **no** production `create_dir_all` at all. Across all of `src/`, the only non-test ones are `io.rs:273` (`ensure_output_dir`) and `fastq.rs:739` (writer creating its parent). Cite `fastq.rs:739` instead. |
| §2 table: *"`--phred64` input pre-flight — `main.rs:391`, gated on `fastqc_requested()`"*, and §5 step 4: *"Change the gate from `cli.fastqc_requested()` to `cli.fastqc_requested() \|\| cli.fastqc_raw`"* | ⚠️ **incomplete.** The actual gate (`main.rs:391-393`) is `cli.phred64 && (cli.fastqc_requested() \|\| matches!(cli.output_format, OutputFormat::UBam))`. A literal textual substitution still produces correct code, so this is a description error rather than a code error — but anyone reasoning from the plan's version will mis-model which runs already fire it. |
| §2 step 4: *"Specialty dispatch — `:609`–`:709`"* | ❌ `hardtrim5` at `606`, `hardtrim3` `637`, `clock` `668`, `implicon` `690` returning at `720`, `clump_only` `722`. Span is 606-722, off at both ends. Harmless drift, but step 7 ("before the specialty dispatch (`:609`)") points three lines into the `hardtrim5` block. |
| §2 step 5: *"paired path from ~`:982`"* | ❌ `run_paired_trim` is **called at `:781`** and defined at `:1048`. Line 982 is inside `run_clump_only`'s paired-to-BAM arm. |
| §3.5 / §7.1: *"the nf-core module's `cores = task.cpus - 4` arithmetic"* | ⚠️ `- 4` is the paired branch; single-end is `- 3` (`main.nf:31-33`). Conclusion unaffected. |

### Could not verify

- **No MultiQC was actually run.** §2.1 below traces `_clean_s_name` line by line, which is a strictly
  better class of evidence than reading `config_defaults.yaml`, but V5 (an end-to-end run) is still
  the only thing that proves the whole chain including report discovery and duplicate handling.
- The `fastqc-rust` `=1.1.1` pin's behaviour on a truncated or colorspace FASTQ (relevant to §1.2) —
  read the dispatch, not the parser.

---

## 1. Logic review

### 1.1 CRITICAL — the insertion point sits ahead of every collision pre-flight and every sanity check but one

The plan picks `:557`-`:606` because no mode escapes it. That is true (verified). What the plan does
not notice is **what else lives after `:606`**:

- every `naming::preflight_output_collisions` call — specialty modes at `:607+`, interleaved-uBAM at
  `:761`, and the ones inside `run_paired_trim` / `run_single_end_trim`;
- `run_clump_only`'s own guards, e.g. *"--clump_only --paired requires two FASTQ input files"*
  (`:814`);
- **`sanity_check_any` on every input other than `cli.input[0]`.** Verified: the only pre-dispatch
  call is `main.rs:337` on `cli.input[0]`. The rest are at `:982`, `:984`, `:1131`, `:1133`, `:1141`
  (passthrough), `:1197`, `:2292`, `:2294`, `:2345` — all inside the drivers.

Two consequences, one of which is a factually wrong claim in the plan:

**(a) §3.6's edge-case row is wrong.** The plan states:

> | Empty / truncated / colorspace input | Already refused by `sanity_check_any` at `:337`, before the new call site. |

That holds for the **first** input only. A paired run with a clean R1 and a colorspace or truncated R2
will now hand R2 to `fastqc-rust` before TG ever looks at it. The user gets whatever `fastqc-rust`
says about a colorspace file (which, per `bam.rs:440`'s final `else`, it accepts as FASTQ without
erroring) instead of TG's own `colorspace_file.fastq` diagnosis — and `test_files/colorspace_file.fastq`
exists precisely so that diagnosis is tested.

**(b) It breaks an invariant the codebase states in its own words.** `main.rs:366-370`:

> *"Placed ahead of every dispatch branch — including the early-returning specialty and clump-only
> paths — and ahead of `ensure_output_dir`, so a rejected run creates nothing."*

and `main.rs:756-757`:

> *"Ahead of setup_trimming, which scans up to 1 M reads for adapters — a refusal must not arrive
> after that wait."*

A full FastQC pass over every input is a far longer wait than the 1 M-read adapter scan the codebase
already guards against, and it creates `fastqc_pre/` plus two files per input before the refusal
lands. The plan does not break this silently out of malice — it simply never looks past `:606` — but
it must name the breach and decide on it rather than inherit it.

**Fix (cheap, keeps the single insertion point):**

1. Inside `run_pretrim_fastqc`, loop `sanity_check_any(input)?` over **all** of `cli.input` before
   the first `fastqc::run`. It is a peek plus one record, so the cost is nil, and it restores §3.6's
   claim to true. (`detect_input_format` has already been called for all inputs at `:346`, so the
   format verdict is in hand; the plan's step 5 re-detects — reuse `input_formats` instead.)
2. For the collision pre-flights, either accept-and-document ("a run refused for an output collision
   may already have written `fastqc_pre/`") or hoist. Accept-and-document is defensible — the reports
   are correct raw reports and regenerate identically — but it must be a stated decision, with a
   validation that pins it (see V11 below), not an accident.

### 1.2 IMPORTANT — the §3.3 guard is one-directional

§3.3 is correct and well-argued, but it states the condition too narrowly. The real correctness
requirement is **`fastqc-rust`'s extension dispatch must agree with TG's content verdict**, and
`bam.rs:433-452` dispatches on four suffix families, not two:

```rust
if   name_lower.ends_with(".bam") || name_lower.ends_with(".ubam") { BAMFile::open(path, true, false) }
else if name_lower.ends_with(".sam")   { BAMFile::open(path, false, false) }
else if name_lower.ends_with(".fast5") { Fast5File::open(path) }
else                                   { FastQFile::open(config, path) }   // "assumed to be FASTQ"
```

The plan guards only *uBAM content + non-`.bam`/`.ubam` name*. It leaves unguarded:

- **uBAM content named `*.sam`** → opened as SAM, not BAM. Not covered by the plan's suffix test
  (`.sam` ends in neither `.bam` nor `.ubam`, so it *is* refused — fine, accidentally correct).
- **FASTQ content named `*.bam` / `*.sam` / `*.fast5`** → misdispatched the other way. TG accepts it
  (content detection says FASTQ); `fastqc-rust` tries BAM/SAM/HDF5. This almost certainly errors
  loudly rather than lying, but *"almost certainly"* is the plan's own standard for §3.3 and it did
  not apply it here. §3.3's headline — *"This is the only way the feature can produce a silently
  wrong report"* — is asserted, not demonstrated, for this direction.

The cleanest resolution is not a wider guard. See Alternative A in §4 — `fastqc-rust` has a
`sequence_format` config field that overrides extension dispatch entirely, which makes the whole
class disappear.

Also: the plan says the suffix test should be **"ASCII-lowercased, matching `bam.rs:440`"**.
`bam.rs:437` uses `to_string_lossy().to_lowercase()`, which is **Unicode** lowercasing, not ASCII.
Practically identical for `.BAM`; the plan should not claim a match it does not have.

### 1.3 IMPORTANT — `--fastqc_args -o <dir>` silently defeats `fastqc_pre/`

`fastqc::run` builds the config with `output_dir` from its parameter (`fastqc.rs:60`), **then**
applies `--fastqc_args` on top (`fastqc.rs:72-74`). And `apply_fastqc_args` honours `-o` / `--outdir`
(`fastqc.rs:125-130`):

```rust
"-o" | "--outdir" => { config.output_dir = Some(PathBuf::from(v)); }
```

So `--fastqc_raw --fastqc_args "-o qc"` sends **both** report sets to `qc/`. On disk nothing is
overwritten (the names differ), but:

- `fastqc_pre/` does not exist, so the nf-core glob `path("fastqc_pre/*.zip")` matches nothing and,
  being `optional: true`, **emits an empty channel with no error** — the exact silent-downstream
  failure §7.1 exists to prevent;
- Assumption 1 ("`--fastqc_args` applies to both runs") makes this the *designed* behaviour, so it
  will not be caught by surprise during implementation.

Decide explicitly: refuse `-o` under `--fastqc_raw`, strip it for the raw run only, or document it as
the user's problem. Whichever — V2 must cover it (it currently does not).

### 1.4 IMPORTANT — §7.2's `-d` remedy does not work in the nf-core path

§7.2 tells methylseq to *"run MultiQC with `-d`/`--dirs` so directory names prefix sample names"*.
Verified against the nf-core MULTIQC module (`modules/nf-core/multiqc/main.nf:12`):

```groovy
tuple val(meta), path(multiqc_files, stageAs: "?/*"), ...
```

Every input file is staged into its own **numbered** directory — `1/`, `2/`, `3/`. The `fastqc_pre/`
component is gone by the time MultiQC sees the files. So:

- `-d` yields `1 | sample`, `2 | sample` — distinct, but numbered by staging order, i.e. unstable
  across runs and meaningless in the report;
- the other remedy ("declare two FastQC sections") in practice means `module_order` with
  `path_filters`, which matches against the staged path — also no `fastqc_pre` in it. methylseq has
  no `module_order`/`path_filters` today (verified: zero hits in `assets/`, `conf/`).

So **neither** of §7.2's two proposed remedies survives contact with nf-core's staging. §7.2 is right
that there is a real dependency and right about which shapes collide; it is wrong about the fix. This
matters because §7.2 is the section methylseq will act on. See Alternative B.

### 1.5 Minor logic gaps

- **`--demux` is unlisted.** It reaches `fastqc::run` at `main.rs:1689-1695` inside
  `run_single_end_trim`, i.e. after the insertion point, so `--fastqc_raw --demux` works and reports
  on the un-demultiplexed input. Correct, but §3.6 should say so — `--demux` is in the Perl-parity
  validation matrix.
- **`--clumpify` is unlisted.** `--clumpify` reorders inside the pipeline and does not rewrite
  `cli.input` (verified: `main.rs:219-291`, `1536`, `1763`, `2059` — no input substitution), so
  `--fastqc_raw` reports the as-delivered order. Correct and worth one row.
- **`--output-format ubam` needs no §3.4a entry**, because raw QC is a statement about the input and
  is independent of output format. The plan never says this; given how long the §3.4a exclusion list
  is, a reviewer will ask.
- **FIFO / process-substitution inputs are already refused** — `format.rs:258-268`
  (`non_restartable_kind`) bails on fifos and sockets inside `detect_input_format`, which already
  runs for all inputs at `:346`. So reading each input twice is safe. The plan never raises the
  double-read question; it happens to be fine, and §6 would be stronger for saying why.

---

## 2. Assumptions

### 2.1 §7.2 — stronger than the plan claims (verified from MultiQC 1.35 source)

§11 risk 1 calls this *"the weak evidence"*. It is not weak; it is right, and here is the trace
rather than the config file.

`multiqc/base_module.py:1046-1071` — the FastQC module takes `s_name` from the `Filename` header in
`fastqc_data.txt` (`modules/fastqc/fastqc.py:345-347`, and the module docstring at `:87` says so
explicitly: *"`Filename` in `fastqc_data.txt`, not based on the FastQC report names"*), then:

1. for each `fn_clean_exts` entry, **truncate** at its first occurrence anywhere:
   `trimmed_name.split(pattern, 1)[0]`;
2. then **one** pass over `fn_clean_trim`, stripping a match only as a **suffix or prefix**:
   `if trimmed_name.endswith(characters): trimmed_name = trimmed_name[:-len(characters)]`.

Tracing the nf-core shapes (the module symlinks inputs to `${prefix}.fastq.gz` /
`${prefix}_{1,2}.fastq.gz`, `trimgalore/main.nf:50,61-62`, so these are the real strings):

| | `Filename` | after `fn_clean_exts` | after `fn_clean_trim` |
| --- | --- | --- | --- |
| PE raw | `S1_1.fastq.gz` | `S1_1` | `S1_1` |
| PE trimmed | `S1_1_val_1.fq.gz` | `S1_1_val_1` | `S1_1_val_1` (ends `_1`, so `_val` is not terminal) |
| SE raw | `S1.fastq.gz` | `S1` | `S1` |
| SE trimmed | `S1_trimmed.fq.gz` | `S1_trimmed` | **`S1`** ← `_trimmed` is terminal |

**Paired is safe, single-end collides — confirmed at source level.** One thing to add: a MultiQC
sample-name collision is not a cosmetic merge, it is an **overwrite** (MultiQC logs *"Duplicate sample
name found! Overwriting"*), so in SE one of the two report sets is silently dropped from the report,
and which one depends on file-discovery order. §7.2 should say that — "merge" undersells it.

V5 is still worth running, but it is now confirming a source-level derivation rather than a guess.

### 2.2 Assumption 1 — the `--fastqc_args` wart is bigger than Open-1 admits

Open-1 frames it as *"a confusing refusal"* in specialty modes. §1.3 above is the second half: in
normal modes it silently relocates the raw reports out of the directory the whole integration story
depends on. Both come from the same root (`--fastqc_args` is one knob for two runs). Resolving Open-1
without resolving §1.3 fixes the smaller half.

### 2.3 Assumption 2 — sequential is asserted, never costed

§3.5 says sequential *"adds wall time but no concurrent thread or memory demand"*, then §6 repeats it,
then §3.5 concludes *"cost-neutral against the status quo it replaces"*. The last step is doing more
work than it can carry: it is cost-neutral **only after §7.1 and §7.2 land and methylseq actually
drops its separate FASTQC process**. Until then it is pure added wall time for everyone who sets the
flag, and §11 risk 4 already knows the downstream may not land. Say "cost-neutral *once methylseq
drops its FASTQC process*", not "cost-neutral".

Also unexamined: in the nf-core module `cores = task.cpus - 4` (PE), so on an 8-cpu task FastQC gets
4 threads during a phase where nothing else runs — half the allocation idles through the raw pass.
Out of TG's control, but it is the concrete shape of "adds wall time".

### 2.4 Assumption 5 — now confirmed, upgrade the prose

§3.2 hedges (*"almost certainly must create the report directory itself. Confirm in V6 rather than
assuming"*) and §8 lists it as *"Flagged for confirmation, not asserted"*. It is confirmed: see the
ledger. Keep the `create_dir_all` and keep V6, but drop the hedge — and fix the `specialty.rs:556`
citation, which is test code.

### 2.5 Unstated assumption — the raw report is of the *input*, and nothing checks that

The function being called is `fastqc::run(output_path, …)`. The parameter is named `output_path`;
twelve existing call sites pass an output path; three of them carry comments reasoning about
`output_path` ending in `.bam`. This plan passes an input to a parameter called `output_path` and
does not rename it. That is the single likeliest place for a later maintainer to "fix" the new call
site back into the old shape. Either rename the parameter to `path` (twelve mechanical call-site
edits, no behaviour change) or add a one-line comment at the new site. The plan's §4 docstring
explains the uBAM guard but not this.

---

## 3. Efficiency

Broadly fine; the plan's §6 is accurate as far as it goes.

- **One extra full decompress+parse per input, serialised ahead of trimming.** Correct. For PE that is
  two passes. No change to the trimming hot path. ✓
- **Memory:** one `fastqc-rust` analysis at a time, so peak is unchanged from the existing post-trim
  path. ✓ Verified that `resolved_threads` is only consulted for the `--clumpify` layout
  (`main.rs:249`), so §3.5's "leaves it untouched" holds. ✓
- **Redundant work the plan introduces:** step 5 calls `detect_input_format(input)?` again, when
  `main.rs:346-350` has already built `input_formats: Vec<InputFormat>` for exactly these paths. Reuse
  it. (Negligible cost — a 24-byte peek plus at most one BGZF block — but re-deriving a verdict the
  function already holds is how the two verdicts drift apart later.)
- **Not an issue, worth stating:** reading each input twice is safe because fifos and sockets are
  already refused in `detect_input_format` (`format.rs:258-268`), which runs for every input before
  this point. Add one line to §6 — otherwise the first reviewer to think of `<(zcat …)` has to
  re-derive it.

---

## 4. Validation sufficiency

The brief asks, for each validation, whether it could actually fail on the mistake it exists to
catch. Three cannot, or only half-can.

### CRITICAL — V8 cannot fail on the mistake it exists to catch

> | V8 | Phred+64 pre-flight extends | Phred+64 fixture with `--fastqc_raw` and no `--fastqc` | Pre-flight fires as it does for `--fastqc` |

The mistake V8 exists to catch is *forgetting step 4* (not extending the gate at `main.rs:391`). If
step 4 is forgotten, what happens? `fastqc::run` is still called with `phred64: true`, and
`fastqc.rs:85-95` already has a bespoke failure path for exactly this:

```rust
anyhow::anyhow!(
    "FastQC analysis failed (exit code {code}). --phred64 declares the \
     input as Phred+64, and FastQC rejects a file holding any quality \
     character below ASCII 64 …"
)
```

**The run fails either way.** A V8 that asserts "the run is refused", or asserts on exit status, or
asserts on the substring `--phred64`, passes identically with and without step 4. It is structurally
incapable of detecting the omission.

To make it able to fail, assert on what *only* the pre-flight produces:

- the distinctive pre-flight text — `contains quality character` **and** `within the first` (from
  `main.rs:371-379`), neither of which appears in the `fastqc.rs` message; and
- **that no report was written** — `fastqc_pre/` must be absent or empty. Without step 4, FastQC runs
  first and fails during analysis; with step 4, nothing is read at all. That filesystem assertion is
  the cleanest discriminator and it is one line.

`tests/integration_fastqc_phred64.rs` and `test_files/phred64_test.fastq` already exist, so this has a
home and a fixture.

### CRITICAL — V1 cannot catch a report computed on the wrong file

> | V1 | Reports are produced and land in the right place | … | `fastqc_pre/<input>_fastqc.{html,zip}` per input; trimmed reports absent |

V1 asserts existence and location. The feature's entire value proposition is that the report
*describes pre-trim data*. V1 cannot distinguish a correct raw report from a report of the wrong file
in any scenario where both files exist (`--fastqc_raw --fastqc`, V2's scenario, which is also the
scenario methylseq will run).

This repo has already learned this lesson and written it down —
`tests/integration_fastqc_phred64.rs:10-11`:

> *"Every other FastQC test asserts only that `*_fastqc.zip` exists. These read inside it, because a
> wrong number in the report is invisible from outside."*

Use the mechanism that file already provides (`basic_stat(zip, stem, measure)` reads
`fastqc_data.txt` out of the archive) and assert something that differs between raw and trimmed:

- `Sequence length` — a fixed value for untrimmed input (e.g. `100`), a range (`<min>-<max>`) after
  quality/adapter trimming. This is the sharpest single discriminator and it is one `basic_stat` call.
- optionally `Total Sequences` (differs whenever `--length` filtering drops reads).

Without this, V1 is an existence check wearing a correctness check's label.

### IMPORTANT — V2 does not cover the `-o` override

Per §1.3, `--fastqc_args "-o X"` silently collapses the two sets into one directory and empties the
nf-core glob. V2 is the validation that owns "the two sets coexist" and it does not exercise the one
input that breaks coexistence. Add it, asserting whatever behaviour §1.3 is resolved to.

### Good as written

- **V3 is the best-built validation in the plan.** It names the failure mode, and it explicitly
  includes the inverse — *"Then prove it can fail: the same file as `reads.bam` must succeed"*. That
  second clause is precisely the discipline V1 and V8 are missing. Keep it, and keep the
  "no report written and no trimming done" clause, which is a real filesystem assertion.
- **V5** can fail. Now confirming a source-level derivation (§2.1) rather than a guess — but see §1.4:
  it should also test the **nf-core staging shape** (`stageAs: "?/*"`), not just a flat directory,
  because the flat-directory result will look fine and mislead methylseq.
- **V6** can fail: with `--output_dir` fresh, `ensure_output_dir` (`io.rs:268-276`) creates the output
  dir at `main.rs:561` but `fastqc_pre` beneath it still will not exist. Add the `--output_dir`
  **unset** case (`fastqc_pre/` in cwd) — the plan lists it in §3.6 and does not validate it.
- **V7** can fail on an accidental default-path change. Reword the assertion: *"no previously-passing
  test fails"*, not *"883+ tests green"* — a count grows when you add tests and tells you nothing.
- **V4** can fail. Note `--clock` and `--implicon` are paired modes and need an even input count.

### Missing validations

| # | What it would catch |
| --- | --- |
| V9 | `--fastqc_args "-o X"` + `--fastqc_raw` — §1.3, the empty-glob silent failure |
| V10 | **Colorspace or truncated R2** + `--fastqc_raw` → TG's own refusal, no report written. §1.1(a); `test_files/colorspace_file.fastq` and `truncated.fq.gz` exist. This is the validation that would have caught the wrong §3.6 row. |
| V11 | A run refused for an output collision + `--fastqc_raw` → pins whether `fastqc_pre/` was written. §1.1(b). `tests/integration_output_collision.rs` and `integration_preflight_tripwire.rs` are the precedent. |
| V12 | Raw report *content* reflects untrimmed data — folded into V1 above |
| V13 | Two inputs sharing a basename (Open-3) → pins the behaviour instead of leaving it silent |
| V14 | `--fastqc_raw --demux` and `--fastqc_raw --clumpify` — unlisted modes that reach the new site |

---

## 5. Alternatives

### A. Set `fastqc-rust`'s `sequence_format` from TG's content verdict, instead of refusing (§3.3)

**This is the strongest alternative and it removes a whole section of the plan.**

`FastQCConfig` has `pub sequence_format: Option<String>` (`fastqc-rust-1.1.1/src/config.rs:41`), and
`bam.rs:400-431` checks it **before** extension dispatch:

```rust
if let Some(ref format) = config.sequence_format {
    match format_lower.as_str() {
        "bam" => return Ok(Box::new(BAMFile::open(path, true, false)?)),   // all records
        "fastq" => return Ok(Box::new(FastQFile::open(config, path)?)),
        …
    }
}
```

TG builds that config itself (`fastqc.rs:60-70`) and already knows the content verdict (`input_formats`
at `main.rs:346`). Passing `sequence_format: Some("bam")` for a content-detected uBAM makes
`reads.dat` **work** instead of being refused, and it fixes the inverse direction from §1.2 (FASTQ
content named `*.bam`/`*.sam`/`*.fast5`) at the same time. It is reachable: `apply_fastqc_args` does
not translate `-f`/`--format` (`fastqc.rs:105-137`), so no user input can contradict it.

Cost: one extra parameter on `fastqc::run` (or a sibling entry point), twelve call sites passing
`None`. Benefit: §3.3, its refusal message, its test, and the §1.2 asymmetry all disappear, replaced
by five lines that make the engine read what TG already knows the file is.

Counter-argument for keeping the refusal: it is more conservative, and `#421` establishes
"refuse where a sound report cannot be produced" as house style. But `#421` refuses where **no report
exists**; here a perfectly sound report is one config field away. Worth putting to Felix explicitly
rather than shipping the refusal by default.

### B. Separate the two sets by **name**, not only by directory

§1.4 shows the directory cannot reach MultiQC through nf-core's staging. MultiQC keys on the
`Filename` header, which `fastqc-rust` takes from `path.file_name()` (`runner.rs:284-290`). So if the
raw pass read the input through a staged alias — `fastqc_pre/S1_raw.fastq.gz` → `S1.fastq.gz` — then
both the report filenames **and** the in-report `Filename` carry `_raw`, and `_raw` is in neither
`fn_clean_exts` nor `fn_clean_trim` (verified in 1.35), so MultiQC separates `S1_raw` from `S1` with
**zero pipeline-side configuration** and no `-d`, no `module_order`, no `path_filters`.

Trade-offs: needs a symlink (note `feedback_hardlink_truncate_hazard` — symlink, never hardlink) and
the report says `S1_raw` rather than `S1`, which some will read as a wart. But it is the only option
canvassed that actually delivers what methylseq needs, and it composes with `fastqc_pre/` rather than
replacing it — keep the directory for the module glob, add the name suffix for MultiQC. Worth
costing before the directory-only design is locked, because the directory name ships permanently.

### C. nf-core-side: emit the directory, not the files

`path("fastqc_pre")` (the directory) rather than `path("fastqc_pre/*.zip")` would survive
`stageAs: "?/*"` as `1/fastqc_pre/…`, letting `path_filters: ['*fastqc_pre*']` work. Pipeline-side
only, nothing changes in TG. Cheaper than B but pushes config onto every consumer; B pushes it onto
none. Mention it in §7.1 as the fallback if B is rejected.

### D. Open-1 without a new flag

The plan offers `--fastqc_raw_args` or "narrow the implication (magic, but one line)". There is a
third, narrower option that touches neither `fastqc_requested()` nor the flag surface: change the
`#421` gate from `if self.fastqc_requested()` to
`if self.fastqc || (self.fastqc_args.is_some() && !self.fastqc_raw)`. One condition, in the block that
already owns this decision, and `fastqc_requested()` stays exactly as §3.4 requires.

### E. Concurrent raw pass

§3.5 presents sequential as settled policy; Assumption 2 concedes it is reversible but never costs the
alternative. Running the raw pass on a background thread while trimming proceeds would hide most of
the added wall time (inputs are regular files — fifos are already refused, §1.5 — so a second
concurrent read is safe). It costs exactly what §3.5 says: `resolved_threads` and the `--clumpify`
memory layout would have to account for the overlap. Not recommended for v1; worth one sentence in
§3.5 so the decision reads as a trade-off taken rather than a constraint discovered.

### F. Flag name (Open-4)

Low stakes, one observation: `--fastqc_pre` is coherent with the `fastqc_pre/` directory the plan
already fixes, so a user who sees one can guess the other. `--fastqc_raw` reads better in isolation
and matches how nf-core labels these sections ("FastQC (raw)"). If `fastqc_pre/` stays the directory
name, `--fastqc_pre` is the more self-documenting pair; if Alternative B lands and reports gain a
`_raw` suffix, `--fastqc_raw` becomes the coherent one. Decide the two together.

---

## 6. Action items

### Critical

1. **Make V8 capable of failing.** Assert on the pre-flight's own wording (`contains quality
   character` + `within the first`) **and** that no report was written. As specified, V8 passes
   identically whether or not step 4 is implemented, because `fastqc.rs:85-95` fails the run anyway.
2. **Make V1 capable of failing.** Assert raw report *content* — `Sequence length` fixed for raw vs a
   range for trimmed — using the in-zip reader already in
   `tests/integration_fastqc_phred64.rs`. Existence-and-location cannot detect a report computed on
   the wrong file.
3. **Fix §3.6's sanity-check row and close the gap.** `sanity_check_any` runs pre-dispatch on
   `cli.input[0]` only (`main.rs:337`); every other input is checked inside the drivers (`:1131`,
   `:1133`, `:1141`, `:1197`, `:2292`, `:2294`, `:2345`, `:982`, `:984`). Loop `sanity_check_any` over
   all inputs at the top of `run_pretrim_fastqc`, and add V10.
4. **Name the ordering decision.** The insertion point precedes every `preflight_output_collisions`
   call. `main.rs:366-370` states *"so a rejected run creates nothing"* and `main.rs:756-757` states
   *"a refusal must not arrive after that wait"*. Decide explicitly whether `--fastqc_raw` opts out,
   say so in §3.5/§3.6, and pin it with V11.

### Important

5. **Resolve `--fastqc_args -o`** (§1.3): it silently sends both sets to one directory and empties the
   nf-core `fastqc_pre/*` glob with no error. Add V9.
6. **Rewrite §7.2's remedy** (§1.4): `-d/--dirs` and `path_filters` both fail under nf-core's
   `stageAs: "?/*"` staging. V5 must test the nf-core staging shape, not a flat directory. Note also
   that a MultiQC name collision **overwrites**, it does not merge.
7. **Decide §3.3 vs Alternative A**: `FastQCConfig.sequence_format` (`config.rs:41`, honoured at
   `bam.rs:400-431`) makes `reads.dat` work rather than be refused, and closes the inverse direction
   in §1.2. Put it to Felix rather than defaulting to the refusal.
8. **Cost Alternative B before fixing `fastqc_pre/`** as the separation mechanism. Directory names do
   not reach MultiQC through nf-core; filenames do. The directory name ships permanently.
9. **Fix the `specialty.rs:556` citation** — it is test code, and `specialty.rs` has no production
   `create_dir_all`. Cite `fastq.rs:739`. Drop Assumption 5's hedge; it is confirmed.
10. **Reuse `input_formats`** (`main.rs:346`) in step 5 rather than re-calling `detect_input_format`.
11. **Address Open-3 rather than deferring it.** `planned_secondary_outputs`' own doc comment
    (`main.rs:120-124`) argues that *"the pre-flight's output-vs-input check is only as complete as
    the path list it is given"* — this feature adds outputs that are not in that list, which is the
    exact scenario the comment was written for. Add V13 even if the resolution is last-writer-wins.

### Optional

12. Fix the line-number drift: specialty dispatch is `606-722` not `609-709`; the paired path is
    `:781` / defined `:1048`, not `~:982`; the `--phred64` gate also includes
    `|| matches!(cli.output_format, UBam)`.
13. Note that the nf-core module's SE branch uses `cores = task.cpus - 3`, not `- 4`.
14. Add `--demux` and `--clumpify` rows to §3.6, and one line saying why no §3.4a `--output-format
    ubam` exclusion is needed.
15. Add a line to §6 recording that fifos/sockets are already refused in `detect_input_format`, so
    the second read of each input is safe.
16. Soften §3.5's "cost-neutral" to "cost-neutral once methylseq drops its FASTQC process" — §11 risk
    4 already concedes the downstream may not land.
17. Either rename `fastqc::run`'s `output_path` parameter to `path`, or comment the new call site:
    twelve call sites and the parameter name all say "output", and this is the first that does not.
18. Add the `--output_dir`-unset case to V6, and reword V7's assertion as "no previously-passing test
    fails".
19. Decide Open-4 together with Alternative B (§5F).

---

## 7. Overall

This is a careful plan. The structural claim it rests on — one insertion point between `:557` and the
specialty dispatch serves every mode — is **verified true**, and the three self-corrections in §11 are
all real and all read from source. §3.3 is a genuine catch that the originating sketch would have
missed, V3 is built the right way round, and §7.2's conclusion is not merely plausible but correct at
source level, which the plan undersells.

The gaps cluster in one place: the plan traced `main()` **down to** `:606` and did not trace what
comes after it. Everything in §1.1 — the collision pre-flights, the per-input sanity checks, the
"rejected run creates nothing" invariant — lives past that line. And the two validations that cannot
fail (V8, V1) share a shape: both assert an outcome that the broken implementation also produces.

Fix the four Critical items and this is ready to implement.
