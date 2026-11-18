# Progress: `--fastqc_raw` — FastQC on the input, before trimming

**Last updated:** 2026-11-18 (code review round 1 applied)

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` **revision 2** — two design reversals after review (§12) |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — agreed findings folded in |
| Impl Plan | 📋 Planned | — |
| Implementation | 📋 Planned | Nothing blocking — the symlink mechanism is verified (PLAN §8.9) |
| Code Review | ✅ Complete | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` — 2 Critical + 2 High + 5 Medium; all fixed bar M4 (PLAN §14) |
| Coverage | ✅ Complete | `COVERAGE.md` — was INCOMPLETE; the code gap and 4 of 5 missing validations now closed |

## Notes

- Originates from @ewels for nf-core/methylseq: pre-trim QC in one command, without publishing a
  `fastqc-rust` nf-core module.
- **Revision 2 reduces the downstream requirement to a container bump.** Revision 1 needed a module
  change *and* MultiQC configuration; revision 2 needs neither for correctness.
- Both reviewers agreed on five findings and contradicted each other on none. The two that reversed
  design decisions:
  - **The subdirectory could not do its job.** `fastqc_pre/` was chosen to fix the MultiQC collision,
    but nf-core's MULTIQC stages with `stageAs: "?/*"` (flattening it away) and MultiQC takes its
    sample name from the report's internal `Filename` line, not the path. Renaming the *report* fails
    for the same reason. Only changing the *analysed path's* filename reaches that field — hence the
    `_raw` symlink.
  - **The uBAM refusal was unnecessary.** `FastQCConfig.sequence_format` (`config.rs:41`) overrides
    extension dispatch at `bam.rs:400-431`, and TG already holds `input_formats` at `main.rs:343`.
    Declaring the format accepts `reads.dat` instead of refusing it, and subsumes the guard's
    one-sidedness — `.sam`, `.fast5` and FASTQ-misnamed-`.bam` were all unguarded.
- **The structural fix:** revision 1's "one insertion point serves every mode" was an invariant
  breach. Every `preflight_output_collisions` runs at `:608` or later, so a raw pass at `:557-606`
  would write reports and read every input before any collision refusal — against the convention
  stated at `main.rs:368-371` and the one CLAUDE.md calls load-bearing. The call now sits in each
  dispatch branch after its pre-flight, at the cost of the single-site elegance.
- **Three validations could not fail on the mistake they existed to catch** (V1, V4, V11 in the new
  numbering) — the failure mode this repo has hit before. Each now states what makes it fail.
  V2 is new and is the one that matters: it reads the `Filename` line out of the zip, because the
  whole design rests on that field carrying `_raw`, and an existence check cannot see it.
- Nine source citations in revision 1 were wrong or drifted; all corrected. Notably
  `specialty.rs:556`, cited as production precedent for `create_dir_all`, is test code
  (`#[cfg(test)]` at `:520`) — `fastq.rs:739` is the real precedent.
- **The symlink mechanism is verified, so nothing gates implementation.** fastqc-rust takes
  `file_name()` from the path as given (`runner.rs:283-286`) before opening anything, `open_group`
  passes it through unchanged, and no input path is ever `canonicalize`d or `read_link`ed — the sole
  `canonicalize` (`runner.rs:321`) is on the report output directory, for a write lock. So the link's
  name reaches the report while the target supplies the data. What remains of Open-2 is portability
  only: no symlinks means a hard error, covered by V9.

## History

- 2026-11-18: **V5 verified empirically** — MultiQC 1.14.dev0 reads raw and trimmed as distinct
  samples, single-end and paired. Counterfactual confirms the alias is load-bearing: rewriting the
  report's `Filename` field back collapses it to one sample, dropped silently. The nf-core bump's
  last evidence gap is closed.

- 2026-11-18: Code Review + Coverage → ✅ (both reports + audit; C1/C2/H1/H2/M1/M2/M3/M5 fixed,
  M4 deliberately unchanged, 904 tests). Both reviewers and the audit independently found the same
  Critical: raw reports pooled into one directory and silently overwrote each other at exit 0.

- 2026-11-18: Plan → ✅ **Revision 2** (two design reversals, one invariant fix, 9 citations corrected,
  3 unfailable validations rewritten + V2/V9/V10/V11 added)
- 2026-11-18: Plan Review → ✅ Complete (both reports; 5 agreed findings, 0 contradictions)
- 2026-11-18: Plan → ✅ Complete (`PLAN.md` revision 1 written; 2 design questions resolved with Felix)
