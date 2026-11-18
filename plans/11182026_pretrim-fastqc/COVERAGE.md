# Plan Coverage Report — `--fastqc_raw`

**Mode:** B (code vs. plan, post-implementation; the plan carries its own implementation outline in §5)
**Plan:** `plans/11182026_pretrim-fastqc/PLAN.md` (revision 2, with §13 implementation notes)
**Code:** `feat/fastqc-raw` @ `cf3bc84`, diffed against `origin/dev`
**Date:** 2026-11-18
**Verdict:** **INCOMPLETE — 6 items unresolved** (1 behaviour the plan specified and the code does not
implement, 4 validations with no test at all, 1 validation that cannot fail as written)

## Summary

| | Count |
|---|---|
| Total ledger items | 74 |
| Implemented | 56 |
| Implemented with deviation — **disclosed** in §13 | 5 |
| Implemented with deviation — **not disclosed** | 2 |
| Partially implemented | 4 |
| Not implemented | 5 |
| Cannot determine | 2 |

Build state verified independently, not taken from §13:

- `cargo test` — **896 passed, 0 failed** across all suites (exit 0). §13's "883 before, +13 new"
  reconciles exactly: 6 new unit tests in `src/cli.rs`, 7 new integration tests in
  `tests/integration_fastqc_raw.rs`.
- `cargo fmt --all -- --check` — clean.
- `cargo clippy --all-targets --release -- -D warnings` — exit 0, re-run after touching
  `main.rs`/`cli.rs`/`fastqc.rs`/`integration_fastqc_raw.rs` so it was a real 1m check of all targets
  rather than a warm-cache no-op.

---

## Coverage ledger

### §3.1 — the flag

| # | Item | Status | Evidence |
|---|---|---|---|
| B1 | `--fastqc_raw` bool, independent of `--fastqc` | Implemented | `src/cli.rs:353-354`; test `fastqc_raw_does_not_imply_fastqc` (`src/cli.rs:2798`) |
| B2 | `--fastqc_raw` alone → raw yes, trimmed no | Implemented | `fastqc_requested()` left untouched (`src/cli.rs:617-619`); asserted at `src/cli.rs:2806-2810` |
| B3 | `--fastqc` alone → today's behaviour unchanged | Implemented | No change to the `fastqc_requested()` gate or the 12 post-trim call sites beyond the new `None` argument |
| B4 | `--fastqc_raw --fastqc` → both sets | Implemented | `raw_and_trimmed_reports_coexist_with_distinct_sample_names` (`tests/integration_fastqc_raw.rs:125`) |
| B5 | `--fastqc_args` still implies `--fastqc`, so `--fastqc_raw --fastqc_args X` produces both; documented in `--help` | Implemented (indirectly documented) | Behaviour: `src/cli.rs:617-619` unchanged. `--help` states the implication on `--fastqc_args` ("Implies --fastqc"), but the `--fastqc_raw` entry does not restate the consequence the plan called "a deliberate, documented consequence" |

### §3.2 — naming

| # | Item | Status | Evidence |
|---|---|---|---|
| B6 | A `_raw`-suffixed symlink is staged per input, FastQC is run on the symlink, the symlink is removed | Implemented | `src/main.rs:128-157`; `the_staged_alias_does_not_survive_the_run` (`tests/integration_fastqc_raw.rs:81`) |
| B7 | Suffix inserted before the *full* extension (`sample.fastq.gz` → `sample_raw.fastq.gz`) | Implemented | `raw_alias_path` (`src/main.rs:86-100`) via `naming::strip_fastq_extensions` (`src/io.rs:816`), with an `rfind('.')` fallback for non-FASTQ-family names |
| B8 | Reports written flat into `--output_dir`, no subdirectory | Implemented | `output_dir` passed straight to `fastqc::run` (`src/main.rs:145-152`); every integration test asserts the zip at `dir.path().join(...)` |
| B9 | The report's internal `Filename` carries `_raw`, so MultiQC reads a distinct sample — SE and PE | Implemented | `single_end_raw_report_names_the_sample_raw:68` and `paired_raw_reports_are_named_per_mate:112` both read `fastqc_data.txt` out of the zip |

### §3.3 — format declared, not inferred

| # | Item | Status | Evidence |
|---|---|---|---|
| B10 | Every raw run declares `sequence_format` from TG's own format knowledge | Implemented, **deviation disclosed** (§13 #4: re-detected per input rather than threading `input_formats`) | `src/main.rs:123-126` (`detect_input_format` → `fastqc::SeqFormat`), `src/fastqc.rs:94,100` |
| B11 | A content-detected uBAM named anything is analysed as BAM; a FASTQ misnamed `*.bam` likewise | Implemented, **untested** | Mechanism present as above; no test exercises either direction — see Gap 2 |
| B12 | No refusal is added for non-`.bam`-suffixed uBAM | Implemented | No such guard exists anywhere in `src/` |

### §3.4 / §3.5 / §3.6 / §3.7

| # | Item | Status | Evidence |
|---|---|---|---|
| B13 | `--fastqc_raw` + `-o`/`--outdir` inside `--fastqc_args` is refused at the CLI layer, naming the conflict and pointing at `--output_dir` | Implemented | `src/cli.rs:1161-1172` + `fastqc::args_set_outdir` (`src/fastqc.rs:65-71`), which mirrors `apply_fastqc_args`' tokenisation including its valueless-`-o` no-op; tests at `src/cli.rs:2860`, `:2880`, `:2892` |
| B14 | Sequential: the raw pass completes before trimming; `cli.cores` passed; `resolved_threads` and the clumpify memory layout unchanged | Implemented | `src/main.rs:149` passes `cli.cores`; `resolved_threads` (`src/fastqc.rs:73`) untouched; `src/main.rs:249` untouched |
| B15 | The raw pass runs **after** each dispatch branch's `preflight_output_collisions`, in every branch | Implemented, **deviation disclosed** (§13 #2 wrapper, §13 #5 run-once guard) | `preflight_then_raw_qc` (`src/main.rs:177-190`) called at all 13 sites (`:719, :751, :874, :989, :1035, :1081, :1138, :1227, :1312, :2356, :2403, :2466, :3083`); `naming::preflight_output_collisions` now appears exactly once in `main.rs`, enforced by `tests/integration_preflight_tripwire.rs:1225-1234` |
| B16 | **The raw reports' own paths are added to the planned-output set** passed to `preflight_output_collisions`, so a raw-report collision is refused by the same mechanism rather than a parallel one | **Not implemented** — and **not disclosed** in §13 | No planned-output builder references the raw report names; `grep -rn fastqc_raw src/` returns only `main.rs:184`, `main.rs:503`, and the `cli.rs` flag/guard. What exists instead is a narrower parallel mechanism: an `alias.exists()` bail on the *transient symlink* (`src/main.rs:129-135`). See Gap 1 |
| B17 | `run_pretrim_fastqc` sanity-checks **every** input before analysing any of them | Implemented | `src/main.rs:112-116`, ahead of the per-input loop at `:122` |
| B18 | The guard is not described as covering mid-stream truncation | Implemented | Plan text corrected (§3.7, §3.8, V5) and the narrower scope pinned by `mid_file_truncation_surfaces_from_the_qc_pass` (`tests/integration_fastqc_raw.rs:224`) |

### §3.8 — edge-case table (13 rows)

| # | Row | Status | Evidence |
|---|---|---|---|
| E1 | Paired input → one report per mate, distinct names | Implemented, tested | `paired_raw_reports_are_named_per_mate:96` |
| E2 | Multiple pairs / multi-file SE → one report per input; two inputs sharing a basename are refused by the pre-flight | Implemented, tested **for the SE-trim case only** | `a_collision_refused_run_writes_no_raw_report:156` passes because the two inputs' *trimmed* outputs collide, not because the raw report paths are in the planned set (B16). Two inputs whose trimmed names differ but whose raw reports coincide (e.g. `a/s.fastq` + `b/s.fastq.gz`) are not refused |
| E3 | `--output_dir` unset → reports and symlink in the working directory | Implemented, untested | `src/main.rs:120` (`output_dir.unwrap_or_else(\|\| Path::new("."))`); every test passes `--output_dir` |
| E4 | `--basename` does not apply to raw reports | Implemented, untested | `raw_alias_path` names from the input only; `cli.basename` is never read in `run_pretrim_fastqc` |
| E5 | Specialty modes allowed; `--fastqc` stays refused there | Implemented | `fastqc_raw_accepted_with_every_specialty_mode:2813` and `fastqc_still_refused_in_specialty_modes_alongside_fastqc_raw:2841`. That reports are actually *written* in those modes is untested — see Gap 3 (V6) |
| E6 | `--clump_only` allowed; the five `clump_only.rs` FastQC sites otherwise untouched | Implemented, untested | Raw pass reached via `src/main.rs:989, 1035, 1081, 1138`; the `clump_only.rs` sites changed only by the new trailing `None` |
| E7 | `--demux` allowed; one pre-demux report per input | Implemented, untested | Demux dispatches at `src/main.rs:1809`, downstream of `run_single_end_trim`'s pre-flight at `:1312` |
| E8 | uBAM input, any filename → analysed correctly | Implemented, **untested** | Same as B11 / Gap 2 |
| E9 | `--phred64` pre-flight gate widened | Implemented, **untested** | `src/main.rs:500-504`; no test — see Gap 5 (V11) |
| E10 | Input already named `*_raw.*` → `sample_raw_raw.fastq.gz`, ugly but non-colliding | Implemented by construction, untested | `raw_alias_path` appends unconditionally |
| E11 | Symlink creation failure → hard error naming the path | Implemented, **untested** | `src/main.rs:129-143` (both the pre-existing-path bail and the `symlink` `with_context`); no test — see Gap 4 (V9) |
| E12 | Empty or colorspace input in any position → refused before any analysis, no report written | Implemented; colorspace tested, empty **not** tested | `a_bad_second_input_is_refused_before_any_report_is_written:193` covers colorspace. The empty case is claimed as hand-verified in §13's iteration log, not as a test |
| E13 | Mid-stream truncation → run fails from the QC pass, having reported on earlier inputs | Implemented, tested | `mid_file_truncation_surfaces_from_the_qc_pass:224` |
| E14 | FIFO / pipe / stdin already refused by #379 | Implemented (pre-existing, untouched) | `src/cli.rs:536`, `src/format.rs:342` unchanged |

### §3.9 — interaction with the #421 refusals

| # | Item | Status | Evidence |
|---|---|---|---|
| B19 | `fastqc_requested()` and the #421 refusals are **not** modified | Implemented | `src/cli.rs:617-619` and the refusal block byte-identical to `origin/dev`; `--fastqc_raw` is not added to either |

### §5 — implementation outline, steps 1–9

| # | Step | Status | Evidence |
|---|---|---|---|
| S1 | `cli.rs`: add the flag beside `--fastqc`, with a doc comment naming (a) runs on the input before trimming, (b) `*_raw_fastqc.*`, (c) allowed in the specialty modes that refuse `--fastqc` | Implemented (all three clauses present) | `src/cli.rs:348-354`; verified in `--help` output. Undisclosed extra: `alias = "fastqc-raw"` — additive, no behaviour removed |
| S2a | Leave `fastqc_requested()` untouched | Implemented | As B19 |
| S2b | Optional `fastqc_raw_requested()` helper ("if it reads better at call sites") | Absent — **not a gap**, the plan made it conditional | `grep -rn fastqc_raw_requested src/ tests/` → no match |
| S3 | `validate()` refuses `-o`/`--outdir` inside `--fastqc_args` | Implemented | As B13 |
| S4 | Widen the `--phred64` pre-flight gate with `\|\| cli.fastqc_raw` | Implemented | `src/main.rs:503` |
| S5a | `fastqc.rs`: add the format parameter, set it **before** `apply_fastqc_args` | Implemented, **deviation disclosed** (§13 #1: `SeqFormat` enum rather than `Option<&str>`) | `src/fastqc.rs:88-110` — set in the struct literal at `:100`, `apply_fastqc_args` called at `:108`, and `apply_fastqc_args` never writes `sequence_format`, so an explicit format cannot be clobbered |
| S5b | Update all twelve callers to pass `None` | Implemented — exactly 12 | `src/main.rs:1801, 2115, 2123, 2137, 2616, 2793, 2936` (7) and `src/clump_only.rs:~309, ~505, ~513, ~803, ~1065` (5) |
| S6a | `run_pretrim_fastqc`: `sanity_check_any` over every input | Implemented | `src/main.rs:112-116` |
| S6b | Per input: build `<dir>/<stem>_raw<ext>`, symlink to the **absolute** input path, run FastQC, remove the link including on the error path | Implemented | `src/main.rs:136-157`; `std::fs::canonicalize` for the absolute target at `:136`; `remove_file` before `result?` at `:155-156` |
| S6c | `create_dir_all` the output directory first | Implemented | `naming::ensure_output_dir` (`src/io.rs:268-276`) called at `src/main.rs:119` |
| S6d | Signature `run_pretrim_fastqc(cli, inputs, input_formats, output_dir)` (§4) | Implemented, **deviation disclosed** (§13 #3) | `fn run_pretrim_fastqc(cli: &Cli)` at `src/main.rs:111` |
| S7a | Call it in each dispatch branch immediately after that branch's pre-flight | Implemented, **deviation disclosed** (§13 #2, #5) | As B15 |
| S7b | Add the raw-report paths to the planned-output set | **Not implemented**, **not disclosed** | As B16 / Gap 1 |
| S8.1 | Test: `--fastqc_raw` parses; does not imply `--fastqc`; does not set `fastqc_requested()` | Implemented | `fastqc_raw_does_not_imply_fastqc:2798` |
| S8.2 | Test: accepted with each of `--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`; `--fastqc` still refused there | Implemented | `fastqc_raw_accepted_with_every_specialty_mode:2813` (all four), `fastqc_still_refused_in_specialty_modes_alongside_fastqc_raw:2841` |
| S8.3 | Test: `--fastqc_raw` + `--fastqc_args "-o X"` refused | Implemented | `fastqc_raw_refuses_outdir_inside_fastqc_args:2860`, plus two negative-control tests at `:2880` and `:2892` |
| S8.4 | Test: integration SE **and** PE — zip exists and the in-zip `Filename` carries `_raw` | Implemented | `single_end_raw_report_names_the_sample_raw:49`, `paired_raw_reports_are_named_per_mate:96` |
| S8.5 | Test: uBAM fixture copied to a non-`.bam` name → report is BAM-derived | **Not implemented** | No test anywhere references a renamed uBAM under `--fastqc_raw`. See Gap 2 |
| S8.6 | Test: a collision-refused run leaves no raw report | Implemented (partial on the symlink half) | `a_collision_refused_run_writes_no_raw_report:156` asserts no `*_raw_fastqc.zip`; it does not assert the absence of the transient symlink, which V1/V4 both call for |
| S8.7 | Test: truncated/colorspace R2 refused by TG, not by fastqc-rust | Implemented, **deviation disclosed** (§13 iteration #3: the truncation half was found unachievable, the plan's §3.7/§3.8/V5 were corrected, and the test split in two) | `a_bad_second_input_is_refused_before_any_report_is_written:193` + `mid_file_truncation_surfaces_from_the_qc_pass:224` |
| S9a | Docs: `--help` | Implemented | Verified against `./target/debug/trim_galore --help` |
| S9b | Docs: `docs/` flag reference | Implemented | `docs/src/content/docs/guide/flags.md:21` |
| S9c | Docs: `wiki/reference/cli-flags.md` | Implemented | `wiki/reference/cli-flags.md:129` |
| S9d | Docs: a **one-sentence** CHANGELOG entry | Partially implemented | `CHANGELOG.md:4-10` — two sentences across three lines, and no PR link. §13 discloses the missing link only, not the length. It does match the surrounding entries' house style, which are multi-sentence with a link |

### §9 — validations V1–V11

| # | Validation | Status | Test that exercises it |
|---|---|---|---|
| V1 | Reports produced, named, located; no transient symlink left behind | Implemented | `single_end_raw_report_names_the_sample_raw:49` (name + location) + `the_staged_alias_does_not_survive_the_run:81` (symlink removed) |
| V2 | The in-zip `Filename` carries `_raw` — read out of the archive | Implemented | `single_end_raw_report_names_the_sample_raw:66-70`, `paired_raw_reports_are_named_per_mate:112`, `raw_and_trimmed_reports_coexist_with_distinct_sample_names:146-152` |
| V3 | Format is declared, not sniffed — uBAM copied to `reads.dat` shows BAM-derived stats | **Not implemented — no test exists** | — |
| V4 | A refused run creates nothing: no `*_raw_fastqc.*` **and no symlink** | Partially implemented | `a_collision_refused_run_writes_no_raw_report:156` covers the report half only |
| V5 | TG validates every input, not just the first — colorspace R2, TG's own message, no R1 report | Implemented | `a_bad_second_input_is_refused_before_any_report_is_written:193`, asserting the message *and* the absence of R1's report |
| V6 | Specialty modes work **and write reports** — `--hardtrim5`, `--clock`, `--implicon` each produce a `*_raw_fastqc.zip` | **Not implemented — no test exists** | The CLI tests prove acceptance at parse time; none runs a specialty mode end to end with `--fastqc_raw` |
| V7 | The output directory is created | Implemented, **not independently testable as written** | `naming::ensure_output_dir` also runs unconditionally at `src/main.rs:672`, before every dispatch branch, so skipping step 6c could not have produced the failure V7 was written to catch. No test passes a non-existent `--output_dir` with `--fastqc_raw` |
| V8a | Exactly 883 tests plus the new ones | Implemented — verified | 896 passed / 0 failed; 883 + 13 new (6 unit + 7 integration) reconciles exactly |
| V8b | Unchanged Perl 0.6.11 md5s (`validation` CI job) | Cannot determine | Needs the conda-installed Perl `trim_galore` and network; not runnable here. Supporting evidence only: `--fastqc_raw` defaults to false and `preflight_then_raw_qc` adds nothing when it is false, and the 12 existing `fastqc::run` sites pass `None`, preserving extension dispatch |
| V9 | Symlink failure is loud — unsupported filesystem, or a pre-existing path at the link name | **Not implemented — no test exists** | Both error paths are present in the code (`src/main.rs:129-143`); neither is exercised |
| V10 | `-o` in `--fastqc_args` refused | Implemented | `fastqc_raw_refuses_outdir_inside_fastqc_args:2860` |
| V11 | `--phred64` pre-flight widened — Phred+64 fixture with `--fastqc_raw` and **no** `--fastqc` fails with the pre-flight's own wording and writes no report | **Not implemented — no test exists** | `tests/integration_fastqc_phred64.rs` has the `--fastqc` twin (`preflight_refuses_mislabelled_input_before_writing_anything`) and the no-FastQC negative control, but no `--fastqc_raw` case |

### §8 — assumptions the code must honour

| # | Assumption | Status | Evidence |
|---|---|---|---|
| A1 | `--fastqc_args` applies to both runs and still implies `--fastqc` | Honoured | `cli.fastqc_args` passed to the raw run (`src/main.rs:147`); `fastqc_requested()` unchanged |
| A2 | Sequential, not concurrent | Honoured | One blocking `fastqc::run` per input inside the pre-flight wrapper, before any driver runs |
| A3 | `_raw` is the marker, inserted before the full extension, fixed not configurable | Honoured | `raw_alias_path` (`src/main.rs:86-100`); no CLI surface for it |
| A4 | Symlinks are available | Honoured as stated, **portability untested** | `std::os::unix::fs::symlink` (`src/main.rs:138`) — Unix-only, will not compile for a Windows target. §13 "Not done" discloses this |
| A5 | fastqc-rust does not create its report directory, so `create_dir_all` is mandatory | Honoured | `src/main.rs:119` |
| A6 | `sequence_format` overrides extension dispatch | Honoured as plumbed | `src/fastqc.rs:100`; never clobbered by `apply_fastqc_args`. The end-to-end consequence is the untested V3 |
| A7 | MultiQC reads `Filename` from inside the report | Honoured and tested | V2's assertions |
| A8 | The `=1.1.1` fastqc-rust pin does not move | Honoured | `Cargo.toml` and `Cargo.lock` are not in the diff |
| A9 | A symlink's own name reaches the report; its target supplies the data | Honoured and tested | `single_end_raw_report_names_the_sample_raw` asserts both halves — `Filename` = the link name, `Sequence length` = 150, the untrimmed input's length |

### §7 — downstream work described as required

| # | Item | Status | Evidence |
|---|---|---|---|
| D1 | nf-core/modules `trimgalore`: container bump is the only strictly required change | Out of scope for this repo — nothing to verify here | Enabled by B8 (flat reports) + B9 (distinct names), both implemented |
| D2 | MultiQC needs no configuration | Satisfied in-repo | B9 / V2 |
| D3 | nf-core/methylseq: needs the module container bump and nothing else | Out of scope for this repo | — |
| D4 | The Perl validation matrix is not affected | Cannot determine locally | V8b |

---

## Gaps (detail)

### Gap 1 — the raw-report paths are not in the planned-output set (§3.6, §5 step 7; ledger B16 / S7b)

**Expected.** §3.6: *"The raw reports' own paths are added to the planned-output set passed to
`preflight_output_collisions`, so a raw report colliding with another output is refused by the same
mechanism rather than by a parallel one."* §5 step 7 repeats it.

**Found.** Nothing. No planned-output builder in `src/main.rs` or `src/io.rs` constructs a
`*_raw_fastqc.*` path. What was built instead is exactly the parallel mechanism §3.6 set out to
avoid, and a narrower one: an `alias.exists()` bail (`src/main.rs:129-135`) that guards only the
transient symlink's name, not the two report files.

**Consequence, reachable today.** Two single-end inputs whose trimmed outputs differ but whose raw
reports coincide — `a/s.fastq` and `b/s.fastq.gz` give `s_trimmed.fq` vs `s_trimmed.fq.gz`, no
collision, but `s_raw.fastq` and `s_raw.fastq.gz` both reduce to `s_raw_fastqc.zip` — overwrite
silently. §3.8 row E2 asserts this class is "now caught by §3.6's pre-flight"; it is caught only when
the trimmed names happen to collide too, which is what
`a_collision_refused_run_writes_no_raw_report` actually exercises.

**Disclosure.** §13 lists five deviations; this is not among them, and "Not done" does not mention it.

### Gap 2 — §3.3's reversal is entirely untested (ledger B11 / E8 / S8.5 / V3)

**Expected.** §5 step 8: *"uBAM fixture copied to a non-`.bam` name: report is BAM-derived, not
FASTQ-misread (V3)."* V3 adds: *"Prove it can fail by passing `None`."*

**Found.** No test under `tests/` or in `src/` references a renamed uBAM with `--fastqc_raw`. The
plumbing is present and correct by inspection (`src/main.rs:123-126` → `src/fastqc.rs:100`), but the
single behaviour that justified deleting revision 1's uBAM refusal has no regression test. A future
caller passing `None` at `src/main.rs:151` would silently re-break it and the suite would stay green.

**Disclosure.** Not mentioned in §13.

### Gap 3 — V6: no specialty mode is run end to end with `--fastqc_raw` (ledger E5 / V6)

**Expected.** V6: *"`--hardtrim5`, `--clock`, `--implicon` each with `--fastqc_raw` → the specific
`*_raw_fastqc.zip` exists per mode. **Fails if** the call sits after the branch's `return`."*

**Found.** Only parse-time acceptance (`fastqc_raw_accepted_with_every_specialty_mode`). That test
calls `Cli::validate()` and never runs the binary, so it cannot detect a raw call placed after a
branch's `return` — the exact failure V6 names. The hardtrim arms (`src/main.rs:719, 751`) and the
`run_specialty_paired` arm that serves `--clock` and `--implicon` (`:3083`) are correct by
inspection; none is covered.

**Disclosure.** Not mentioned in §13.

### Gap 4 — V9: the symlink failure paths are untested (ledger E11 / V9)

**Expected.** V9: *"Output dir on a filesystem without symlink support, or a pre-existing path at the
link name → hard error naming the path; no partial report."*

**Found.** Both branches exist (`src/main.rs:129-135` pre-existing path, `:138-143` symlink failure)
and both name the path. Neither is exercised. The pre-existing-path half is cheap to test — create
`<dir>/<stem>_raw.fastq.gz` and run.

**Disclosure.** Not mentioned in §13.

### Gap 5 — V11: the widened `--phred64` gate is untested (ledger E9 / S4 / V11)

**Expected.** V11: *"Phred+64 fixture, `--fastqc_raw` and **no** `--fastqc` → the pre-flight's own
wording ('contains quality character', 'within the first') **and no report written**."*

**Found.** The gate is widened (`src/main.rs:503`). `tests/integration_fastqc_phred64.rs` holds the
`--fastqc` twin and a no-FastQC negative control, but nothing passes `--fastqc_raw`. V11 explicitly
warned that `fastqc.rs:85-95` has its own `--phred64` failure path, so the test has to assert the
pre-flight's wording — that distinction is currently unasserted in either direction.

**Disclosure.** Not mentioned in §13.

### Gap 6 — V7 cannot fail as written (ledger S6c / V7)

**Expected.** V7: *"`--output_dir` at a fresh path → created, reports inside. **Fails if** step 6c is
skipped."*

**Found.** Step 6c is implemented (`src/main.rs:119`), but `naming::ensure_output_dir` also runs
unconditionally at `src/main.rs:672`, ahead of every dispatch branch, so the directory exists however
step 6c went. V7 is therefore a validation that could not have failed on the mistake it was written
to catch — the same category of defect §9's preamble says revision 1 shipped three of. No test covers
it either way.

---

## Items noted but **not** counted as gaps

- **§5 step 2's `fastqc_raw_requested()` helper is absent.** The plan made it conditional ("if it
  reads better at call sites"). `cli.fastqc_raw` is read at two sites; the helper would add nothing.
- **V5's truncation half is inverted relative to revision 2's original text.** §13 iteration #3
  documents the discovery that the requirement was unachievable, the plan's §3.7, §3.8 and V5 were
  corrected to match, and the test was split in two. Disclosed and reconciled.
- **Portability (§10 Open-2) is untested.** §13 "Not done" discloses it; the plan itself chose "hard
  error (v1)" and left portability open.
- **`docs/guide/outputs.md` and `wiki/usage/*` are not updated.** Neither was named in §5 step 9;
  §13 "Not done" discloses both. Worth knowing that `guide/flags.md:21` links to
  `/guide/outputs/#fastqc`, a page that does not mention the new flag.
- **The CHANGELOG entry has no PR link.** Disclosed in §13; no PR exists yet.
- **`alias = "fastqc-raw"` was added beyond §5 step 1.** Additive, removes nothing.
- **Not a coverage item, noted once:** in `src/fastqc.rs:25-45` the new `SeqFormat` doc paragraph was
  appended to the end of the pre-existing doc block that documents `run`'s parameters, so that whole
  block (including `run`'s parameter list and the threading note) is now attached to the `SeqFormat`
  enum. Compiles and lints clean; it is a rustdoc attachment matter, not a behaviour one.

---

## Test verification

| Test | File | Result |
|---|---|---|
| `fastqc_raw_does_not_imply_fastqc` | `src/cli.rs:2798` | PASS |
| `fastqc_raw_accepted_with_every_specialty_mode` | `src/cli.rs:2813` | PASS |
| `fastqc_still_refused_in_specialty_modes_alongside_fastqc_raw` | `src/cli.rs:2841` | PASS |
| `fastqc_raw_refuses_outdir_inside_fastqc_args` | `src/cli.rs:2860` | PASS |
| `fastqc_raw_allows_fastqc_args_without_outdir` | `src/cli.rs:2880` | PASS |
| `bare_outdir_with_no_value_is_not_refused` | `src/cli.rs:2892` | PASS |
| `single_end_raw_report_names_the_sample_raw` | `tests/integration_fastqc_raw.rs:49` | PASS |
| `the_staged_alias_does_not_survive_the_run` | `tests/integration_fastqc_raw.rs:81` | PASS |
| `paired_raw_reports_are_named_per_mate` | `tests/integration_fastqc_raw.rs:96` | PASS |
| `raw_and_trimmed_reports_coexist_with_distinct_sample_names` | `tests/integration_fastqc_raw.rs:125` | PASS |
| `a_collision_refused_run_writes_no_raw_report` | `tests/integration_fastqc_raw.rs:156` | PASS |
| `a_bad_second_input_is_refused_before_any_report_is_written` | `tests/integration_fastqc_raw.rs:193` | PASS |
| `mid_file_truncation_surfaces_from_the_qc_pass` | `tests/integration_fastqc_raw.rs:224` | PASS |
| `the_site_inventory_is_complete` (updated: 13 `preflight_then_raw_qc` sites + exactly one `naming::preflight_output_collisions`) | `tests/integration_preflight_tripwire.rs:1215` | PASS |
| `every_site_is_reached_by_some_case` (unchanged — proves `#[track_caller]` preserved per-arm attribution) | `tests/integration_preflight_tripwire.rs` | PASS |
| V3 — uBAM under a non-`.bam` name | — | **MISSING** |
| V6 — specialty mode writes a raw report | — | **MISSING** |
| V9 — symlink failure / pre-existing alias | — | **MISSING** |
| V11 — `--phred64` pre-flight with `--fastqc_raw` alone | — | **MISSING** |

Whole suite: **896 passed, 0 failed** (`cargo test`, exit 0).

---

## Verdict

**INCOMPLETE — 6 items unresolved.**

One is a behaviour the plan specified and the code does not implement:

1. **Add the raw-report paths to the planned-output set** (§3.6, §5 step 7). Currently absent and
   replaced by a narrower `alias.exists()` guard on the transient symlink only. §3.8 row E2's claim
   that same-basename inputs are "caught by §3.6's pre-flight" does not hold in general. **Undisclosed
   in §13** — either implement it, or record the deviation and correct §3.6/§3.8/E2.

Four are validations with no test at all:

2. **V3** — uBAM fixture copied to a non-`.bam` name, asserting BAM-derived stats. This is the sole
   regression guard for the §3.3 reversal that deleted revision 1's uBAM refusal.
3. **V6** — at least one specialty mode run end to end with `--fastqc_raw`, asserting its
   `*_raw_fastqc.zip` exists. Parse-time acceptance cannot catch a call placed after a branch's
   `return`, which is the failure V6 was written for.
4. **V9** — a pre-existing file at the alias path produces a hard error naming it.
5. **V11** — Phred+64 fixture with `--fastqc_raw` and no `--fastqc`, asserting the pre-flight's own
   wording and that no report was written.

One is a validation that cannot fail as written:

6. **V7** — `ensure_output_dir` already runs unconditionally at `src/main.rs:672` before every
   dispatch branch, so V7 could not detect step 6c being skipped. Either retire V7 or rewrite it.

Also outstanding but either disclosed or minor: V4's symlink half is unasserted; the CHANGELOG entry
is two sentences and carries no PR link; `docs/guide/outputs.md` still omits the flag while
`guide/flags.md` links to it. V8's Perl-md5 half cannot be checked outside CI.

Everything else in §3, §5, §8 and §9 is implemented, and the five deviations §13 records
(`SeqFormat` enum, the `preflight_then_raw_qc` wrapper, the `run_pretrim_fastqc(cli)` signature,
per-input format re-detection, and the run-once `AtomicBool`) each check out against the source.
