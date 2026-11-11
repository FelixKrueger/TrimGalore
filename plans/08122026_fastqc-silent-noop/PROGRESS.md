# PROGRESS — #421: refuse `--fastqc` on the arms that cannot run it

**Plan:** `PLAN.md` v1
**Issue:** [#421](https://github.com/FelixKrueger/TrimGalore/issues/421)
**Branch base:** `dev` @ `3ffcef4`

| Step | State | Notes |
|---|---|---|
| Design decision | **done** | Felix, 2026-08-12: refuse at the CLI, not wire or warn |
| Plan written | **done** | `PLAN.md` v1 |
| Felix's manual review | **awaiting** | — |
| Dual plan review | not started | Launch after Felix approves |
| Implementation | not started | Requires the exact trigger |
| Dual code review + coverage audit | not started | Automatic after implementation |
| PR / merge / close | not started | Squash-merge into `dev`; manual close (master is default) |

## Open items carried in the plan

- **O1** — `--paired --fastqc single.fastq` reports the FastQC refusal before the
  "only legal if that file is a uBAM" error. Assumption: acceptable.
- **O2** — each message names a remedy. Assumption: yes, one clause.
- **O3** — `--help` gets one clause, not the full enumeration.

## The thing most likely to go wrong

Arm 5 must be keyed on `output_format` as well as `paired && input.len() == 1`.
`--paired --output-format ubam <one.bam>` is a **capable** arm (`main.rs:2522`) and
`.github/workflows/ci.yml:1173` asserts it produces a report. Dropping the `output_format`
term refuses a working feature and turns the `Validate uBAM input` job red. V3 exists to
catch exactly that.
