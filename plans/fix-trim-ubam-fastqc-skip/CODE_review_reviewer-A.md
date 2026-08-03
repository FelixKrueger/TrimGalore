# Code Review — `fix(ubam): run FastQC on --output-format ubam output`

**Reviewer:** A (independent; dual-review protocol per `~/.claude/CLAUDE.md`)
**Branch:** `fix/trim-ubam-fastqc-skip`
**Commit:** `72540f0a1e1e35651fe6694de7b2f6638c0c460a`
**Base:** `edc8434` (`dev`)
**Diff size:** 3 files, +71 / −0 (`src/main.rs` +32, `.github/workflows/ci.yml` +29, `CHANGELOG.md` +10)
**Date:** 2026-07-25

---

## Summary

**Verdict: APPROVE** — merge-ready as-is. Recommendations below are additive hardening, none blocking.

The fix is correct, minimal, and faithfully mirrors the FASTQ-output path. All six invariants named in the review brief hold. I verified the change empirically rather than by inspection alone: I built the branch and exercised every reachable uBAM-output code path, and confirmed the strongest available correctness property — **the FastQC report generated from the trimmed BAM is byte-identical to the report generated from the equivalent trimmed FASTQ, apart from the `Filename` line.**

The single substantive gap is test *placement*, not test *correctness*: the regression guard exists only in CI, so `cargo test` gives a developer no signal. That matters more than usual here because the bug being fixed was precisely "a call was silently absent" — the failure mode is invisible without an explicit assertion.

No regression risk found. `cargo fmt` clean, `cargo clippy --all-targets --release -- -D warnings` clean, 331 unit tests + all four integration suites pass.

---

## Verification performed

### Invariants from the brief — all confirmed

| # | Invariant | Result |
|---|---|---|
| 1 | `fastqc::run` called *after* report writing | ✅ Code order (`main.rs:1637` > `1635`; `1773` > `1771`; `1894` > `1892`) and runtime stderr order (`Report:` / `JSON report:` precede `Running FastQC on`) both confirm. Report survives a FastQC failure. |
| 2 | Guard is `cli.fastqc \|\| cli.fastqc_args.is_some()` | ✅ Character-for-character identical to `main.rs:883` (SE FASTQ) and `main.rs:1189` (PE FASTQ). |
| 3 | PE calls `fastqc::run` exactly ONCE per pair | ✅ Verified at runtime: 1 zip emitted, and its `fastqc_data.txt` reports `Total Sequences 20` against `samtools view -c` = 20 — one report covering both mates, which is the correct count for a single interleaved BAM. |
| 4 | `&output_path` passed | ✅ Same binding used by the collision pre-flight and `writer.finish()` in each function. No re-derivation, so no drift possible. |
| 5 | `cli.cores` passed, not hardcoded | ✅ And semantically right — see note in Low-6 below. |
| 6 | Multi-pair invokes FastQC N times | ✅ **Intended, not a concern.** Verified with 2 pairs → 2 output BAMs → 2 reports, one per output file. This is exactly the FASTQ path's semantics (one FastQC per output file; the FASTQ path does 2N for N pairs because it emits 2N files). Reports interleave with pair processing, matching the FASTQ path. |

### Byte-identity of the FastQC report (strongest evidence)

```
diff <(grep -v ^Filename  BAM-derived/fastqc_data.txt) \
     <(grep -v ^Filename FASTQ-derived/fastqc_data.txt)
→ IDENTICAL
```

Both report `Total Sequences 10`, `Total Bases 644 bp`, `Sequence length 63-65`,
`Encoding Sanger / Illumina 1.9`. Every QC module's data matches. This proves the
BAM path feeds FastQC exactly the same sequence/quality stream as the FASTQ path.

### All code paths exercised

| Path | Driver | Result |
|---|---|---|
| SE, uBAM input | `run_ubam_output_single` | ✅ 1 report; matches FASTQ path byte-for-byte |
| SE, multi-file (2 inputs) | same, in loop | ✅ 2 reports, one per output |
| PE, single interleaved BAM in | `run_ubam_output_paired_single_file` | ✅ 1 report, 20 seqs (both mates) |
| PE, two FASTQ files in | `run_ubam_output_paired_two_files` | ✅ 1 report for `polyAT_R1_val.bam` |
| PE, multi-pair (2 pairs) | same, `chunks(2)` loop | ✅ 2 reports |
| **PE, two BAM files in** | — | **Does not exist** — rejected up-front at `main.rs:1437-1449`. The brief listed this as a fourth path to check; it is unreachable by construction. |

### Edge cases

| Case | Result |
|---|---|
| `--fastqc_args` without `--fastqc` | ✅ FastQC runs (guard's `\|\|` branch). Matches FASTQ path. |
| `--no_report_file` + `--fastqc` | ✅ FastQC still runs; trimming reports correctly absent. The guard sits *outside* the `if !cli.no_report_file` block, which is right. |
| Empty output (all reads filtered, `--length 500` → 0 records) | ✅ Graceful: html + zip produced, **exit 0**, no panic. |
| Specialty + uBAM (`--hardtrim5 20 --fastqc --output-format ubam`) | ✅ No FastQC — **and none on the FASTQ path either.** `src/specialty.rs` contains zero `fastqc` references, so run-and-exit modes have never invoked FastQC on either output format. Consistent pre-existing behaviour; **there is no fourth gap, the 3-driver scope is complete.** |

### Upstream dependency claim audited

The commit message asserts "fastqc-rust reads .bam natively (verified against the crate's runner.rs)". I independently confirmed this in the vendored crate (`fastqc-rust-1.0.1`):

- `src/sequence/bam.rs:430` — `open_sequence_file` routes `.bam` / `.ubam` to `BAMFile::open(path, true, false)` (`only_mapped = false`, so all records are read).
- `src/runner.rs:283` — `strip_extensions` strips `.bam`, so `ubam_test_trimmed.bam` → basename `ubam_test_trimmed` → `ubam_test_trimmed_fastqc.{zip,html}`. **The filenames the new CI steps assert are exactly correct.**
- **Reverse-complement hazard is a non-issue.** `bam.rs:297-300` revcomps sequence and reverses qualities whenever flag `0x10` is set. TrimGalore writes only `0x04` / `0x4D` / `0x8D` (`src/bam.rs:558-565`, asserted by lib tests at `src/bam.rs:1461-1465`) and *rejects* input records carrying `0x10` (`src/bam.rs:833`). The transform can never fire. Good — this was the most plausible silent-corruption vector and it is closed.

---

## Issues by area

### Logic
No defects found. Guard placement, call count, and path derivation are all correct.

### Errors / error handling
No defects. `?` propagation is consistent with the FASTQ path; report-before-FastQC ordering means a FastQC failure cannot cost the user their trimming report.

### Efficiency
No concerns. One FastQC pass per output file is minimal. `cli.cores` is forwarded so FastQC's rayon pool scales.

### Structure
Three near-identical 8-line blocks. **I recommend against factoring these into a helper** — the guard-plus-call idiom is already duplicated across the FASTQ drivers at `main.rs:883` and `1189`, so the copies are locally consistent with an established house pattern, and a 5-line helper taking `(cli, path, output_dir)` would add indirection without removing meaningful complexity. Leaving it is the right call.

---

## Recommendations

### MEDIUM-HIGH — 1. No `cargo test` coverage; the regression guard is CI-only

`grep -rln fastqc tests/` returns **zero hits**. The entire guard for this fix lives in `.github/workflows/ci.yml`.

Why this is worth closing:
- The bug fixed here is "a required call was silently missing." Its recurrence mode is a future refactor of the three drivers dropping the call again — caught only after push, in a job that needs `apt-get install samtools` and a release build.
- `tests/integration_ubam_out.rs` **already has the exact harness needed** — `binary()` (`:23`) and `fresh_tmpdir()` (`:27`). No new dependency, no new fixture.
- **I measured the cost: 209 ms** for a full `--output-format ubam --fastqc` run on `ubam_test.bam`. Negligible against the existing suite.

Suggested additions to `tests/integration_ubam_out.rs`:

```rust
#[test]
fn ubam_out_se_fastqc_produces_report() {
    let dir = fresh_tmpdir("tg_int_ubam_out_se_fastqc");
    let status = Command::new(binary())
        .args(["--output-format", "ubam", "--fastqc"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    assert!(dir.join("ubam_test_trimmed.bam").exists(), "output BAM missing");
    assert!(
        dir.join("ubam_test_trimmed_fastqc.zip").exists(),
        "FastQC zip missing — --fastqc silently skipped on the uBAM-output path"
    );
    assert!(dir.join("ubam_test_trimmed_fastqc.html").exists(), "FastQC html missing");
}

#[test]
fn ubam_out_pe_fastqc_produces_exactly_one_report() {
    let dir = fresh_tmpdir("tg_int_ubam_out_pe_fastqc");
    let status = Command::new(binary())
        .args(["--paired", "--output-format", "ubam", "--fastqc"])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    assert!(dir.join("ubam_paired_test_val_fastqc.zip").exists(), "FastQC zip missing");

    // PE uBAM output is ONE interleaved BAM, so exactly ONE report.
    let zips = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().ends_with("_fastqc.zip"))
        .count();
    assert_eq!(zips, 1, "expected exactly 1 FastQC zip for interleaved PE output, got {zips}");
}
```

### MEDIUM — 2. The new CI steps may not be fail-fast; some assertions are non-gating

The `validation-ubam` job overrides the shell at `.github/workflows/ci.yml:679`:

```yaml
defaults:
  run:
    shell: bash -l {0}
```

GitHub injects `-eo pipefail` only for the bare `bash` shell *keyword*; a custom template string is used verbatim. So `-e` and `-o pipefail` are very likely **absent**, meaning only the **last** command in each `run:` block determines the step's outcome.

Consequences for the two new steps:

- **SE step** (`:855-862`): the last command is `unzip -l <zip> | grep -q fastqc_data.txt`, which transitively requires the zip to exist — so the load-bearing assertion does gate, by luck of ordering. But `test -f ..._fastqc.html` (`:861`) is **non-gating**: the HTML could stop being produced and CI would stay green.
- **PE step** (`:870-877`): the last command *is* the `ZIP_COUNT` assertion, so this one gates correctly.

The job itself hints at the ambiguity — line 774 issues a bare `set -e` to "restore" after `set +e`, which only makes sense if the author assumed `-e` was already on.

Fix (2 lines, removes all doubt and makes every assertion gating):

```yaml
        run: |
          set -euo pipefail
          mkdir -p /tmp/ubam_out_se_fqc
          ...
```

Add to both new blocks. Worth considering for the whole job separately, but that is out of scope here.

### LOW — 3. Stale module doc in `src/fastqc.rs`

Lines 11-14 and 23 are now factually wrong. The module doc says:

> Public entry point: [`run`]. **The two callers** in `src/main.rs` (one per single-end output, two per paired-end output) invoke it with the **trimmed FASTQ path** …

There are now **seven** call sites, and three of them pass a BAM. Line 23 likewise: `` `output_path` — the FASTQ to analyse ``.

Suggested replacement for lines 11-14:

```rust
//! Public entry point: [`run`]. Callers in `src/main.rs` invoke it with a
//! trimmed output path — FASTQ on the default path, BAM on the
//! `--output-format ubam` path (`fastqc-rust` dispatches on file
//! extension and reads `.bam` natively). Output `*_fastqc.html` and
//! `*_fastqc.zip` land in the trim-galore `--output_dir` (or current dir
//! if unset).
```

and for line 23:

```rust
/// `output_path` — the trimmed output to analyse (FASTQ, or BAM when
///                 `--output-format ubam` is in effect).
```

*Not fixed directly:* Reviewer B is editing concurrently with no shared state, so a same-file edit risks a clobber; and per the project's phase separation this is a review, not an implementation, turn. Exact text supplied so it can be applied once.

### LOW — 4. The fix silently depends on the output filename ending in `.bam`

The whole mechanism rests on extension-based dispatch inside the upstream crate (`fastqc-rust-1.0.1/src/sequence/bam.rs:430`). This is worth an explicit note because it runs *against* this project's stated convention — `CLAUDE.md` is emphatic that uBAM detection is content-based, not filename-based. If the output naming ever changed to, say, `*_trimmed.ubam.gz` or an extensionless stem, FastQC would silently fall through to `FastQFile` and emit a garbage report rather than erroring.

Consider tightening the existing comment on the three new blocks:

```rust
// Run FastQC if requested. fastqc-rust dispatches on file EXTENSION
// (not content), so this relies on output_path ending in `.bam`.
```

### LOW — 5. Failure-artifact upload doesn't cover the new output directories

`Upload uBAM validation outputs on failure` (`:879-886`) uploads only `/tmp/via_samtools/` and `/tmp/via_trimgalore/`. The two new `_fqc` directories are not included — nor are the pre-existing `/tmp/ubam_out_se/`, `/tmp/ubam_out_pe/`, `/tmp/ubam_out_aux/`. A pre-existing gap that this PR widens slightly. Adding `/tmp/ubam_out_se_fqc/` and `/tmp/ubam_out_pe_fqc/` would make a red run diagnosable from the artifact alone.

### LOW — 6. The `--cores` startup NOTE is now marginally overbroad

`src/main.rs:220-225` prints:

> `NOTE: --output-format ubam uses single-threaded compression in v1; --cores {} is ignored.`

The qualifier "compression" is accurate, but the trailing clause reads absolute, and after this change `cli.cores` *is* honoured — by FastQC's rayon pool. Passing `cli.cores` is the correct choice (it matches the FASTQ path and the user did ask for N cores); only the message could mislead. Optional tweak: "… `--cores {}` is ignored for BAM compression (FastQC, if requested, still uses it)."

---

## Out of scope — pre-existing bug discovered during review

**`--phred64` + `--output-format ubam` writes spec-non-conformant quality scores.** Not introduced by this PR; flagging because I found it while auditing the quality path that FastQC now consumes. Worth a separate issue.

`BamWriter::write_record` (`src/bam.rs:~580`) unconditionally does:

```rust
let raw_qual: Vec<u8> = record.qual.bytes().map(|b| b.saturating_sub(33)).collect();
```

But `--phred64` input keeps Phred+64 ASCII all the way through the pipeline — only `quality.rs` applies the offset, for trimming decisions. So the subtraction is wrong by 31 on that path. Verified: a `'h'` (ASCII 104, true **Q40** under Phred+64) is stored in the BAM as raw **71** instead of 40. `--phred64` is **not** in `cli.rs`'s §3.4a uBAM exclusion list, so the combination is reachable.

**This does not corrupt the new FastQC report.** Two errors cancel: fastqc-rust adds 33 back (→ `'h'`), then its encoding auto-detector classifies the file as Phred+64 and subtracts 64, recovering Q40. I confirmed the report reads `Encoding: Illumina 1.5` with per-base mean `40.0` — correct. The BAM itself remains wrong for every other consumer, though: `samtools fastq` would emit `'h'` for downstream tools to read as Q71.

Suggested resolution (either is defensible): have the writer subtract `cli.phred_offset()` instead of the literal `33`, or add `--phred64` to the §3.4a rejection list.

---

## Quality gates

| Gate | Result |
|---|---|
| `cargo fmt --all -- --check` | ✅ clean |
| `cargo clippy --all-targets --release -- -D warnings` | ✅ clean |
| `cargo test --release` | ✅ 331 unit + 1 + 2 + 8 + 13 integration; 0 failed |
| CHANGELOG | ✅ correctly placed under `### Unreleased` → `#### Fixes`; accurate, user-facing, and the "covering both mates in the PE case" claim is empirically true (20 sequences in one report) |
| Commit message | ✅ Accurate. Explains the what, the why, the one-vs-two report reasoning, and the upstream-capability basis. The `fastqc-rust reads .bam natively` claim checks out. |

---

## Conclusion

Ship it. The fix does exactly what it claims, on every reachable path, and I have byte-level evidence that the resulting report matches the FASTQ path's. The three duplicated blocks are the right level of abstraction for this codebase.

Highest-value follow-up is **Recommendation 1** — move the regression guard into `cargo test`, where it costs 209 ms and ~30 lines and protects against exactly the failure mode this PR just repaired. **Recommendation 2** (`set -euo pipefail`) is a two-line change that converts "gates by luck of command ordering" into "gates by construction". Neither blocks merge.
