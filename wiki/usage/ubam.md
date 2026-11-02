---
type: usage
title: uBAM input and output
description: Running unaligned BAM input, writing uBAM output, carrying aux tags, and the v1 limits of the BAM path.
tags: [ubam, bam, aux-tags, interleaved]
page_id: ubam
subject_status: current
updated: 2026-11-02
sources:
  - repo: FelixKrueger/TrimGalore
    path: src/bam.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/main.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/format.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/cli.rs
    commit: 6303152
---

## Summary

uBAM input needs no flag: it is detected from the file's content. Output stays
FASTQ unless you ask for `--output-format ubam`. Paired uBAM is one
interleaved file in both directions, matching samtools, Picard and fgbio.

## Reading uBAM

```bash
trim_galore sample.bam                                  # single-end
trim_galore --paired interleaved.bam                     # paired, one file
```

Detection is by content, not by filename, so a `.bam` extension is neither
required nor sufficient — see `reference/inputs.md`. Every record must be
unmapped, and the check runs per record rather than only on the first.

Paired uBAM must be **one** interleaved file with mates adjacent, which is
what `samtools sort -n`, `samtools collate`, Picard and fgbio all produce. Two
separate uBAM files are not the paired-uBAM shape. Grouped-but-not-adjacent
input is refused with a `samtools collate` pointer.

`--phred64` is refused with uBAM input: BAM stores raw quality values and the
reader always yields Phred+33 internally, so the flag could only corrupt the
output.

## Carrying aux tags into FASTQ

```bash
trim_galore --preserve-tags CB,UB,RX sample.bam
```

Opt-in and comma-separated, compatible with `samtools -T`. Only A, Z, i and f
scalar types round-trip, because the FASTQ intermediate is textual;
`bam::parse_tag_value` rejects B (array) and H (hex) tags explicitly rather
than dropping them silently.

## Writing uBAM

```bash
trim_galore --output-format ubam sample.fastq.gz
trim_galore --output-format ubam --paired sample_R1.fq.gz sample_R2.fq.gz
```

Paired output is one interleaved `<stem>_val.bam`, not two files.

| Behaviour | Detail |
|---|---|
| Threading | Single-threaded compression in v1. `--cores` is ignored for BAM writing and an up-front note says so; FastQC still uses it. |
| Provenance | A `@PG` line naming the program, version and command line is appended on every run. |
| uBAM in, uBAM out | The input's header chain is preserved and our line is added, so output is **not** byte-identical to input. Provenance is kept by addition, not by rewriting. |
| Version bumps | A routine version bump changes the `@PG` version tag, which is why the project's golden-fixture comparisons ignore `@PG`. |

## v1 limits

`--output-format ubam` is refused together with `--dont_gzip`, `--clumpify`,
`--passthrough`, `--clock`, `--implicon`, `--demux` and `--retain_unpaired`.
The multi-output and header-encoding features are out of scope for v1, and the
refusals are explicit messages rather than silent no-ops.

`--clump_only` does support uBAM in and out, writing `<stem>_clumped.bam`.

## Links

- `reference/inputs.md` — how uBAM is detected and what is refused.
- `reference/outputs.md` — the BAM output names.
- `troubleshooting.md` — the refusals in full.
- `usage/clumpify.md` — `--clump_only` with BAM.

## Open questions

- None.
