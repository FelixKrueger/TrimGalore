---
type: usage
title: RRBS and non-directional libraries
description: Running MspI-digested RRBS data, what --rrbs changes, and the non-directional variant.
tags: [rrbs, msp1, bisulfite, methylation]
page_id: rrbs
subject_status: current
updated: 2026-11-02
sources:
  - repo: FelixKrueger/TrimGalore
    path: src/trimmer.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/main.rs
    commit: 6303152
  - repo: FelixKrueger/TrimGalore
    path: src/cli.rs
    commit: 6303152
---

## Summary

`--rrbs` removes the 2 bp artifact that the MspI end-repair reaction leaves
behind, so that methylation calls are not made on filled-in bases. It changes
three things: a 2 bp trim, an auto-set clip on R2, and the suppression of
N-trimming.

## Directional RRBS

```bash
trim_galore --rrbs sample.fastq.gz
trim_galore --rrbs --paired sample_R1.fq.gz sample_R2.fq.gz
```

| What happens | Where |
|---|---|
| R1 and single-end reads lose 2 bp from the 3' end — but only when an adapter was found | `trimmer::trim_read` |
| R2 instead gets `--clip_r2 2` set automatically, unless you set it yourself or a `--library` preset supplied it | `main.rs`, which prints a line saying it did this |
| `--trim-n` is ignored | `trimmer::trim_read` |

The 3' trim is conditional on adapter removal because that is the case where
the read ran into the fill-in bases. A read with no adapter did not reach
them, so trimming it would discard real sequence.

N-trimming is suppressed to keep output byte-identical to Perl v0.6.x, whose
RRBS path never passed `--trim-n` through to Cutadapt. Combining `--rrbs` with
`--trim-n` is therefore accepted but has no effect on the Ns.

## Non-directional RRBS

```bash
trim_galore --rrbs --non_directional --paired sample_R1.fq.gz sample_R2.fq.gz
```

Non-directional libraries contain reads from both strands, so the artifact can
sit at either end. `trimmer::trim_read` checks every read, R1 and R2 alike: a
read starting with `CAA` or `CGA` loses 2 bp from the **5'** end, and
otherwise the directional 3' rule applies when an adapter was found.

Trim Galore prints a line at startup saying which of the two rules is in
force, so the run's own log records what was done.

## Why the clip on R2 rather than a trim

R2 in a directional library starts in the filled-in region, so the bases to
remove are at its 5' end and are there whether or not an adapter was found.
A fixed clip is the right instrument for that; a conditional 3' trim is not.
This is why the two mates are treated differently rather than symmetrically.

## Links

- `internals/trimming-pipeline.md` — where the RRBS step sits relative to adapter trimming.
- `reference/cli-flags.md` — `--rrbs`, `--non_directional`, the clip flags and `--library`.
- `usage/trim-reads.md` — the ordinary run these options modify.

## Open questions

- None.
