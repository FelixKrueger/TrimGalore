# Code Review — Reviewer B

**Target:** commit `72540f0a` on `fix/trim-ubam-fastqc-skip` — "fix(ubam): run FastQC on `--output-format ubam` output"
**Reviewer:** B (independent; did not read Reviewer A's report)
**Date:** 2026-07-25
**Scope:** `src/main.rs` (+32), `.github/workflows/ci.yml` (+29), `CHANGELOG.md` (+10). Insertions only, no deletions.

---

## Verdict

**APPROVE.** The fix is correct, minimal, and consistent with the FASTQ-output path it mirrors. I verified all five stated invariants **empirically** (built the branch and ran every code path), plus every edge case raised in the review brief. I found **no correctness defects**.

All findings below are about the *strength of the new regression guard* and *documentation accuracy* — none block merge. The most substantive is **M-1**: the new CI `test -f` assertions cannot fail their step under this job's declared shell, so the steps are stricter-looking than they actually are. The primary regression guard (FastQC ran at all) does still work in both steps.

---

## Verification performed

I did not review this on reading alone — I built the branch (`cargo build --release`, clean) and exercised each path.

### The five stated invariants — all confirmed

| # | Invariant | Result |
|---|---|---|
| 1 | `fastqc::run` called **after** report writing | ✅ All three sites sit after the `if !cli.no_report_file` block (`main.rs:1637` after 1635; `1773` after 1771; `1894` after 1892). Confirmed by observed stdout order: `Report:` → `JSON report:` → `Running FastQC on …`. |
| 2 | Guard is `cli.fastqc \|\| cli.fastqc_args.is_some()` | ✅ Textually identical to the FASTQ path (`main.rs:883`, `1189`). Also matches documented behaviour — `--fastqc_args` help says "Implies `--fastqc`" (`cli.rs:285`). |
| 3 | PE path calls `fastqc::run` exactly **once** per pair | ✅ One zip produced. Its `fastqc_data.txt` reports `Total Sequences 20` for the 10-pair fixture — i.e. one report covering **both mates pooled**, which is the correct count for a single interleaved BAM. |
| 4 | `&output_path` — same path as the collision preflight and `writer.finish()` | ✅ Same binding in all three functions; no re-derivation, so no drift possible. |
| 5 | `cli.cores` passed through, not hardcoded | ✅ Passed as `cli.cores` at all three sites. See **L-3** for a nuance. |

### Empirical results

All commands run against the branch build. Every case behaved correctly.

| Scenario | Command shape | Result |
|---|---|---|
| SE uBAM + `--fastqc` | `--output-format ubam --fastqc -o D ubam_test.bam` | ✅ `ubam_test_trimmed_fastqc.{zip,html}` produced; zip contains `fastqc_data.txt`; `Total Sequences 10` |
| PE interleaved + `--fastqc` | `--paired --output-format ubam --fastqc -o D ubam_paired_test.bam` | ✅ Exactly **one** zip `ubam_paired_test_val_fastqc.zip`; `Total Sequences 20` (both mates) |
| **`--fastqc_args` without `--fastqc`** | `--fastqc_args "--nogroup"` | ✅ FastQC runs. Guard's `\|\|` is intentional and matches help text |
| **`output_dir` is `None`** | no `-o` | ✅ Artifacts land **beside the output BAM** (input's parent dir). Nothing leaked into cwd |
| **All reads filtered → empty BAM** | `--fastqc --length 10000` | ✅ **Exit 0**, valid report, no panic. `BAMFile::read_next` maps EOF → `None`, not an error |
| **Multi-input SE** (2 BAMs) | two positional BAMs | ✅ 2 FastQC invocations; distinct `sampleA_trimmed_fastqc.zip` / `sampleB_trimmed_fastqc.zip` |
| **Multi-pair two-file** (4 FASTQs) | `--paired` × 2 pairs | ✅ 2 invocations; `p1_R1_val_fastqc.zip` / `p2_R1_val_fastqc.zip`; no collision |
| Toolchain gates | `fmt` / `clippy -D warnings` / `cargo test --release` | ✅ fmt clean, clippy clean, **331 unit + all integration tests pass** |

### Upstream claims independently checked

The commit message asserts "fastqc-rust reads .bam natively (verified against the crate's runner.rs)". I verified this against the pinned `fastqc-rust 1.0.1` source rather than taking it on trust:

- **BAM dispatch is real.** `src/sequence/bam.rs::open_sequence_file` checks `config.sequence_format` first; when it is `None` it falls through to extension detection, where `.bam`/`.ubam` → `BAMFile::open(path, true, false)` — BAM reader, **all** records (mapped + unmapped). `src/fastqc.rs` leaves `sequence_format` at its default `None`, so extension detection applies and picks BAM. Correct by construction, not by luck.
- **Output naming is right.** `runner.rs::strip_extensions` strips `.bam` (and `.ubam`), so `x_trimmed.bam` → base `x_trimmed` → `x_trimmed_fastqc.zip`. The CI steps assert exactly these names, and my runs confirm them. Had the crate *not* stripped `.bam`, every asserted filename in the new CI steps would have been wrong — this was the single highest-risk assumption in the change and it holds.
- **`output_dir: None` fallback.** `runner.rs` uses the **parent directory of the input file**, not the process cwd. Because the BAM output already lives in that directory, the FastQC artifacts land beside it — the intended placement. (This does mean a doc comment is inaccurate; see **L-1**.)

### Questions from the brief, answered directly

- **Half-written state if FastQC errors after the reports are written?** No. `writer.finish()?` fully finalizes the BAM (including the BGZF EOF block) *before* the report block, and the reports are written and flushed before `fastqc::run`. If FastQC fails, `?` propagates and the process exits non-zero, but the BAM, `.txt`, and `.json` are all complete and valid on disk. Only partial FastQC artifacts could remain, and those are written last inside `process_group`. The ordering invariant delivers exactly what it promises.
- **Race between `writer.finish()` and reading the BAM back?** None. The sequence is strictly sequential on one thread, with the whole report block interposed. `finish()` flushes and closes before any read. uBAM output is single-threaded by design, so there is no concurrent writer to race with.
- **Multi-pair outer-loop misbehaviour?** No. The `fastqc::run` call lives *inside* each per-pair/per-file driver, and the drivers are invoked once per pair (`main.rs:1474-1503`) or per input (`1509-1516`). Each pair therefore gets its own report, named from its own `output_path`. Verified with a real 2-pair run. The existing case-folded collision preflight (`main.rs:1451-1471`) already guarantees distinct output BAM paths, and since zip names derive deterministically from those paths, distinct BAMs imply distinct zips — the preflight transitively protects the FastQC artifacts too.
- **`ls | wc -l | tr -d ' '` portability / newline handling?** Non-issue. The job is `runs-on: ubuntu-latest` only (no macOS leg), GNU `wc -l` does not pad, and `$( )` strips trailing newlines anyway. `tr -d ' '` is harmless defensive habit. Zero-match globbing is handled correctly: `2>/dev/null` swallows `ls`'s error, stdout is empty, count is `0`, and the `||` branch reports it.
- **Could a test pass locally but fail in CI?** No meaningful divergence found. I ran on darwin, CI runs ubuntu-latest, but `fastqc-rust` is pure Rust so artifact names and contents don't vary by OS. `unzip` is present on ubuntu runners and already used at `ci.yml:487`. The one real asymmetry runs the *other* way — see **M-1**: the shell semantics make CI *laxer* than a naive local `bash -e` reading suggests.

### CI wiring hygiene — checked and correct

- Both tmp dirs (`/tmp/ubam_out_se_fqc`, `/tmp/ubam_out_pe_fqc`) are written by **only** these two new steps — no cross-step contamination that could skew the PE count assertion.
- `if: ${{ !cancelled() }}` matches the three immediately preceding uBAM-out steps exactly — established pattern in this job, not an ad-hoc choice.
- Step placement is correct: appended after the existing uBAM-out steps and *before* the `if: failure()` artifact-upload step.
- The upload glob `/tmp/ubam_out_*` (`ci.yml:886`) already captures both new dirs, so failures are debuggable with no extra wiring. Nice reuse.
- `CHANGELOG.md` is well-formed: a single `#### Fixes` heading correctly placed under `### Unreleased` after `#### Changes`. Its technical claims match observed behaviour, including the subtle "covering both mates in the PE case".

---

## Issues by area

### Logic
No issues. The three insertions are identical in shape to the reference FASTQ-path call, the guard matches exactly, ordering is correct in all three, and the paired-path call count (one, not two) is right for the single-interleaved-BAM output shape rather than a copy-paste of the FASTQ path's two calls. The two comments added are the good kind — they explain *why* one call suffices, which is precisely the non-obvious part a future reader would otherwise "fix" into a bug.

### Efficiency
No issues. FastQC dominates runtime here; one pooled pass over the interleaved BAM is strictly cheaper than two passes would be. Threading is delegated to `fastqc-rust`'s rayon pool via `cli.cores`.

### Errors
No issues. `?` propagation matches the FASTQ path. The empty-BAM edge case — the most plausible new failure mode, since the fix adds a read-back of a file that may contain zero records — exits 0 and produces a valid report; I verified this rather than assuming it.

### Structure
The three blocks are near-verbatim duplicates (10 lines × 3). I deliberately do **not** recommend extracting a helper: the FASTQ path already repeats this same block inline, the duplication is 3 lines of payload inside a 2-line guard, and a helper would add indirection for no readability gain. Matching the established local idiom is the right call here, and per the repo's "three similar lines is better than a premature abstraction" guidance this is correct as written.

---

## Fixes applied

**None — deliberately.**

Per this project's `CLAUDE.md`, editing source or CI without an explicit implementation trigger is out of bounds, and a reviewer must not mutate the branch under review. The one change that would otherwise qualify as "unambiguous and low-risk" (**M-1**'s `set -euo pipefail`) is a CI-semantics change on a shared workflow, so I am recommending it with an exact patch rather than applying it.

---

## Recommendations

### M-1 (Medium) — CI: the `test -f` assertions cannot fail their step

The `validation-ubam` job declares `shell: bash -l {0}` (`ci.yml:677-679`). GitHub Actions uses a **custom shell string verbatim**, so neither `-e` nor `-o pipefail` is active — unlike the bare `shell: bash` keyword, which expands to `bash --noprofile --norc -eo pipefail {0}`. I confirmed the consequence locally: a failing mid-script `test -f` under `bash -l` leaves the script's exit status at `0`, whereas the same script under `bash -e` exits `1`.

Per-step consequence:

- **SE step** — `test -f …_fastqc.html` cannot fail the step and nothing downstream depends on the HTML, so **a missing HTML file passes**. The `…_trimmed.bam` check is likewise unenforced. The *zip* check is effectively enforced only because the final `unzip -l <zip> | grep -q fastqc_data.txt` is the last command, and a pipeline's status is its last element's — so an absent zip makes `grep` return 1 and the step fails.
- **PE step** — both `test -f` lines are unenforced. The `ZIP_COUNT` check *is* enforced, because of its explicit `exit 1`. But the count globs `*_fastqc.zip`, so a **misnamed** zip (e.g. `…_val.bam_fastqc.zip`, the exact failure that would occur if `strip_extensions` ever stopped handling `.bam`) still counts as 1 and passes — the asserted *name* is unverified.

**The primary regression guard is intact**: if FastQC were skipped again (zero zips), the SE step fails at `unzip`/`grep` and the PE step fails the count check. So the fix is genuinely protected; the steps are simply stricter-looking than they are, which is a trap for whoever next edits them.

Suggested one-line hardening, as the first line of each new `run:` block:

```yaml
run: |
  set -euo pipefail
  mkdir -p /tmp/ubam_out_se_fqc
  …
```

Note that under `pipefail`, `ZIP_COUNT=$(ls … | wc -l | tr -d ' ')` would abort on the zero-match case *before* printing the friendly "Expected 1 fastqc zip, got 0" message. If that diagnostic is worth keeping, use `ZIP_COUNT=$(ls … 2>/dev/null | wc -l | tr -d ' ' || true)`.

This lax-shell shape is **pre-existing repo-wide** (`ci.yml:487-490` has the same form), so a broader hardening pass is a separate piece of work — but these two new steps are cheap to get right now, and they are specifically regression guards whose strength is the whole point.

### M-2 (Medium) — one of the three changed functions has no new CI coverage

The two new steps cover `run_ubam_output_single` and `run_ubam_output_paired_single_file`. **`run_ubam_output_paired_two_files` receives the same change but is not covered** — that is the FASTQ+FASTQ-in → interleaved-BAM-out path. I exercised it manually and it works correctly (2 pairs → 2 correctly-named zips), so this is a coverage gap, not a defect.

A third step closes it using fixtures already in the repo:

```yaml
- name: Validate PE uBAM-out from FASTQ pair — --fastqc produces one report
  if: ${{ !cancelled() }}
  run: |
    set -euo pipefail
    mkdir -p /tmp/ubam_out_pe2_fqc
    ./target/release/trim_galore --paired --output-format ubam --fastqc \
      -o /tmp/ubam_out_pe2_fqc \
      test_files/BS-seq_10K_R1.fastq.gz test_files/BS-seq_10K_R2.fastq.gz
    test -f /tmp/ubam_out_pe2_fqc/BS-seq_10K_R1_val_fastqc.zip
```

### M-3 (Medium) — no Rust-level test; the guard is CI-shell-only

`cargo test` will not catch a regression of this bug — a contributor could re-break it and see a fully green local run. `tests/integration_ubam_out.rs` is the natural home for an assertion that `--fastqc` yields a `*_fastqc.zip`. Cost is negligible: FastQC on the 10-read fixture completed in well under a second in my runs. Optional, but it would make the guard portable to local development instead of living only in CI.

### L-1 (Low) — stale module doc in `src/fastqc.rs:11-14`

The doc says: *"The two callers in `src/main.rs` (one per single-end output, two per paired-end output) invoke it with the trimmed **FASTQ** path"*. There are now **7** call sites and **three of them pass a BAM**. This PR is what makes the "FASTQ path" claim wrong, so correcting it belongs in this change.

While nearby: `src/fastqc.rs:29-30` states `None` means "current directory (FastQC's default)". `fastqc-rust` actually falls back to the **input file's parent directory** (`runner.rs`). Pre-existing inaccuracy, worth fixing in the same pass since the uBAM path now relies on that fallback landing artifacts beside the output BAM.

### L-2 (Low) — user-facing docs don't mention PE mate pooling

For PE uBAM output, one report covers **both mates pooled** (confirmed: 20 sequences from a 10-pair fixture). A user moving from FASTQ mode to uBAM mode goes from two reports with separate R1/R2 quality curves to one blended report — a different sample count in MultiQC, and the loss of the per-mate view. Since R2 quality is typically worse than R1, pooling can mask real mate-specific degradation. That is an unavoidable consequence of the single-interleaved-BAM output shape and the right behaviour, but it is currently explained **only** in the CHANGELOG. Worth one clause in `--fastqc`'s or `--output-format`'s help text (`cli.rs:279-281`) and/or the README uBAM section, so it is discoverable at the point of use.

### L-3 (Low) — `--cores` interaction with the "cores ignored" NOTE

`CLAUDE.md` and `main()`'s up-front NOTE state that uBAM output is always single-threaded and `--cores N` is silently ignored on this path. Passing `cli.cores` to FastQC is nonetheless **correct** — it matches the FASTQ path, and FastQC's rayon pool is independent of the BAM writer's single-threadedness. But a user who set `--cores 8` and read "ignored" now gets an 8-thread FastQC pass. Harmless and arguably desirable; the imprecision is in the NOTE's wording, not the code. Consider narrowing that NOTE to "ignored for trimming/compression" if it is ever touched. No action required.

### L-4 (Nit) — `tr -d ' '` is inert on this job

`runs-on: ubuntu-latest` only, GNU `wc -l` does not pad, and command substitution already strips the trailing newline. Harmless, and correct-by-default if the job ever gains a macOS leg. Keep or drop; no impact either way.

---

## Summary table

| ID | Priority | Area | Summary | Blocks merge |
|---|---|---|---|---|
| M-1 | Medium | CI | `test -f` assertions non-enforcing under `bash -l {0}`; add `set -euo pipefail`. Primary guard still works | No |
| M-2 | Medium | CI | `run_ubam_output_paired_two_files` has no new coverage (verified working manually) | No |
| M-3 | Medium | Tests | No Rust-level test; `cargo test` can't catch a regression | No |
| L-1 | Low | Docs | `src/fastqc.rs:11-14` "two callers … FASTQ path" now wrong (7 sites, 3 pass BAM) | No |
| L-2 | Low | Docs | PE mate pooling documented only in CHANGELOG, not in help/README | No |
| L-3 | Low | Docs | "`--cores` ignored" NOTE now imprecise for FastQC | No |
| L-4 | Nit | CI | `tr -d ' '` inert on ubuntu-only job | No |

**Bottom line:** the fix itself is clean and I could not break it. Every invariant holds, every edge case in the brief behaves correctly, and the two subtle judgment calls — call FastQC *after* the reports, and call it *once* for PE — are both right and correctly commented. The actionable follow-ups are all in the guard and the docs, with **M-1** the one I would most want addressed before merge, since a regression guard that reads stricter than it behaves invites future erosion.
