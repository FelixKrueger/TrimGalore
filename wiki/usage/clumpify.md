---
type: usage
title: Compression and archival with clumpify
description: Using --clumpify while trimming, --clump_only for reorder-only archival, and how to pick --memory.
tags: [clumpify, compression, archival, memory]
page_id: clumpify
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/clump.rs
    commit: 26ec1cf
  - repo: TrimGalore
    path: src/clump_only.rs
    commit: 26ec1cf
  - repo: TrimGalore
    path: src/main.rs
    commit: 26ec1cf
---

## Summary

Clumping reorders reads so that similar sequences sit next to each other,
which lets gzip find much more redundancy. `--clumpify` does it while
trimming; `--clump_only` does it with no trimming at all, which is the
archival case. Both are bounded by `--memory`.

## While trimming

```bash
trim_galore --clumpify --cores 8 --paired sample_R1.fq.gz sample_R2.fq.gz
```

The run prints the resolved bin layout and its predicted peak memory at
startup, so the budget is visible before the work starts rather than after it
runs out.

## Reorder only, for archival

```bash
trim_galore --clump_only sample.fastq.gz
trim_galore --clump_only --paired sample_R1.fq.gz sample_R2.fq.gz
```

No trimming, no filtering, no adapter detection, no record modification. The
invariant is that every input record appears in the output with its header,
sequence and quality unchanged — only the order of records in the file
differs. That is what makes it safe for archiving data you may need to return
to untouched.

Two normalisations still apply, and they come from the reader and writer
rather than from this mode: the third line of each record is written as a bare
`+`, and CRLF line endings become LF.

`--clump_only` also supports uBAM in and out; see `usage/ubam.md`.

## Choosing --memory

`--memory` is the whole budget, not just the bin pool: the resolver first
subtracts a fixed reservation for the runtime, gzip state, I/O buffers and
allocator retention, charges more when `--fastqc` was requested, and divides
what is left across bins and in-flight workers.

Two consequences in practice:

- **Raising `--cores` raises the memory the run needs.** Bin count scales with
  cores, so a budget that was comfortable at 2 cores can be below the floor at
  16. If clumping silently stops happening after a core-count change, this is
  why.
- **Below the floor, the run does not fail.** It warns loudly and falls back
  to plain trimming, so you get valid output that is simply larger than you
  expected. Check the startup lines if output size surprises you.

Short reads are the other edge: the reservation was calibrated at and above a
floor read length, and anything shorter is outside that fit. The run warns
that peak memory may exceed the prediction and suggests raising `--memory`.

## What you get for it

The win comes from sorting within each bin, not from which bin a read lands
in. Sorted reads put similar sequences inside one gzip window, and the
reported effect on typical Illumina FASTQ is roughly a 20–35% reduction in
output size. How much you actually see depends on library composition rather
than on settings: see the project's own performance pages for the RRBS and
WGBS contrast, where the same mechanism pays off very differently.

## Links

- `internals/parallel-and-clumping.md` — the bin mechanism and the memory model.
- `troubleshooting.md` — the floor warning and the plain-mode fallback.
- `reference/cli-flags.md` — `--clumpify`, `--clump_only`, `--memory`, `--compression`.

## Open questions

- The 20–35% figure and the RRBS-versus-WGBS contrast come from the module
  documentation and the project's performance pages, not from a measurement
  made while writing this page.
