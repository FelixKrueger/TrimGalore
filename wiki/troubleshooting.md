---
type: troubleshooting
title: Refusals, warnings and recovery
description: What Trim Galore refuses, what it only warns about, and the remediation each message points at.
tags: [errors, refusals, recovery]
page_id: troubleshooting
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/fastq.rs
    commit: 7ffe86f
  - repo: TrimGalore
    path: src/format.rs
    commit: a71d9ff
  - repo: TrimGalore
    path: src/cli.rs
    commit: bf28cbd
  - repo: TrimGalore
    path: src/bam.rs
    commit: 115d84f
  - repo: TrimGalore
    path: src/io.rs
    commit: c7e147d
  - repo: TrimGalore
    path: src/main.rs
    commit: 26ec1cf
---

## Summary

Most failures are refusals raised before any output is written, and each one
names what to change. Two cases are warnings rather than refusals, and one
case changes behaviour silently by design. The refusals below are the ones the
code raises; the flag conflicts are enforced in `Cli::validate`.

## Input refused

| Message names | Cause | What to do |
|---|---|---|
| A completely empty file | The first input has no records | Check the file is the one you meant. |
| Not in FASTQ format | First line does not start with `@` | Confirm the file is FASTQ, not FASTA or a compressed archive of something else. |
| SOLiD colorspace | Digits in the first sequence | Colorspace needs separate processing; Trim Galore does not support it. |
| Truncated FASTQ | A record is missing its sequence, `+`, or quality line | Re-fetch or re-generate the file; a partial download is the usual cause. |
| An aligned BAM record | A record does not have the unmapped flag set | Filter first, as the message says, with `samtools view -h -f 4`. |
| Grouped rather than mate-adjacent BAM | Paired uBAM input whose mates are not adjacent | Re-order with `samtools collate`, which the message names. |

The aligned-BAM check runs per record, not only on the first, so a file that
is unmapped at the start and aligned later is still caught.

## Input cannot be re-read

`check_restartable_input` refuses pipes, FIFOs, sockets, and anything that
does not return to the start when reopened. This is the one refusal that
surprises people, because the shell idiom it rejects looks reasonable:

- Process substitution, such as feeding `<(zcat reads.fq.gz)` as the input.
- `/dev/stdin`, and named pipes created with `mkfifo`.

Trim Galore opens each input more than once — format detection, the sanity
check, and adapter auto-detection each open it independently — so it needs a
regular file. `not_restartable_message` says this and gives the fix: write the
stream to a file first, then trim that file. The same check guards the
`--passthrough` and `--demux` files.

## Outputs would collide

`io::preflight_output_collisions` runs before any file is created, and refuses
a run where two planned outputs resolve to one path. The message says the
outputs would be written to the same file, and names the inputs to change.
Paths are compared through `io::collision_key`, which folds case, so a
case-only difference between two inputs still collides on APFS and NTFS where
the filesystem itself would silently overwrite.

## Flags refused together

`--output-format ubam` rejects each of `--dont_gzip`, `--clumpify`,
`--passthrough`, `--clock`, `--implicon`, `--demux` and `--retain_unpaired`.
These are v1 exclusions enforced in `Cli::validate` with a message per flag,
not silent no-ops.

`--passthrough` requires `--paired` and exactly one pair, and is refused
alongside `--retain_unpaired` and `--clumpify`.

`--phred64` is refused with uBAM input: BAM stores raw quality values and the
reader always yields Phred+33 internally, so the flag could only corrupt the
output. A B (array) or H (hex) aux tag is refused by `bam::parse_tag_value`,
which supports only A, Z, i and f scalars — the FASTQ intermediate is textual
and cannot carry the other two.

`--fastqc` is refused on the four specialty modes, and on `--paired` given a
single interleaved uBAM with FASTQ output, because those paths write no
report. Adding `--output-format ubam` to the last one does produce a report.

## Warnings, not refusals

| Condition | Behaviour |
|---|---|
| Shortest input read below the calibrated floor | `main::warn_below_clumpify_floor` warns that peak memory may exceed the predicted peak, and suggests raising `--memory`. The run continues. |
| `--memory` below the usable floor | The run falls back to plain trimming with a loud warning rather than refusing. Clumping is skipped. |

The second case is worth knowing before you tune `--cores`: raising the core
count raises the floor, so a budget that worked at a lower core count can
silently drop clumping.

## A failed run leaves no half-written output

Output files are written through `io::PendingOutput`, which writes a temporary
and renames it into place only once every byte the file owes has been written.
A run refused part-way through therefore publishes nothing, and a previous
file is left as it was. A failed rename deliberately leaves the temporary on
disk rather than deleting the data it was protecting.

## Links

- `reference/inputs.md` — what is accepted, which is the other half of this page.
- `reference/cli-flags.md` — the flags named here, with their defaults.
- `usage/clumpify.md` — choosing `--memory` so the floor warnings do not fire.

## Open questions

- None.
