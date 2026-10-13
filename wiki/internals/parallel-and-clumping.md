---
type: internals
title: Worker pool and read clumping
description: How --cores distributes trimming and compression, how --clumpify reorders reads into bins, and how the memory budget is resolved.
tags: [parallel, clumpify, memory, gzip]
page_id: parallel-and-clumping
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/parallel.rs
    commit: 26ec1cf
  - repo: TrimGalore
    path: src/clump.rs
    commit: 26ec1cf
  - repo: TrimGalore
    path: src/clump_only.rs
    commit: 26ec1cf
---

## Summary

Under `--cores N` each worker trims *and* compresses its own chunk, and the
compressed chunks are concatenated in order. `--clumpify` adds a reordering
stage in front of that, grouping similar reads so gzip finds more redundancy.
Both are bounded by a memory budget that is resolved once at startup.

## The worker pool

`parallel::run_single_end_parallel` and `parallel::run_paired_end_parallel`
distribute reads across N workers. Each worker owns both halves of the cost:
it trims its chunk and gzip-compresses it, emitting an independently
compressed gzip member. Members are written in sequence order, which is valid
because RFC 1952 permits concatenated gzip members — so the output is an
ordinary `.gz` file that any reader accepts.

This replaced an earlier readers-to-one-writer pipeline. The reason is where
the time goes: compression dominates runtime, so funnelling it through a
single thread capped the speedup no matter how many readers were added.

One consequence worth knowing when comparing outputs: member boundaries depend
on how reads were chunked, so the same records can produce different output
bytes at different core counts. The records and their order are unchanged.

## Clumping

`--clumpify` routes each read to one of N in-memory bins. The key is a
canonical k-mer minimizer from `clump::canonical_minimizer`, and the bin is
chosen by `clump::bin_for`, an FNV-1a hash of that key modulo N.

When a bin's accumulated raw bytes exceed its budget, the bin is sorted by
minimizer key and shipped to a worker as one batch, becoming one gzip member.
The compression win comes from that in-bin sort rather than from which bin a
read landed in: sorting puts similar sequences adjacent, and gzip's 32 KB
window then finds the redundancy. Because bin assignment only has to balance
load, no locality-preserving bucketing is needed.

`--clump_only` is the same reordering without any trimming. Its invariant is
that every input record appears in the output with its header, sequence and
quality unchanged — only file-level order differs. Two normalisations still
apply, and they are reader and writer behaviour rather than anything specific
to this mode: the third line is written as a bare `+`, and CRLF line endings
become LF.

## Resolving the memory budget

`clump::resolve_layout` turns `--memory` into a bin count and a per-bin byte
budget. The shape of it:

- `LayoutInputs::n_bins` scales bins with `--cores`, with a floor, and is
  deliberately keyed to the user's `--cores` so that output layout does not
  shift when the other thread counts change.
- `clump::static_bytes_for` reserves the resident cost that sits outside the
  bin pool — runtime, gzip state, I/O buffers, allocator retention — and
  charges extra per FastQC thread when a report was requested, up to a cap.
- What remains is divided across bins and in-flight workers. The divisor
  carries a per-bin and a per-worker coefficient, which *bound* the real
  per-record cost rather than modelling it: part of that cost belongs to the
  allocator, and no single per-record constant is exact at more than one read
  length.

The coefficients were fitted on a byte-equalised paired-end matrix across
several read lengths and core counts, summarised by the worst repetition per
cell. Reads shorter than the calibrated floor are outside that fit, which is
why `main::warn_below_clumpify_floor` warns rather than silently trusting the
prediction.

If the derived per-bin budget falls below the minimum, `resolve_layout` fails
rather than producing a degraded layout. `clump::clumpify_min_memory_bytes`
exists so the caller can detect that case first and fall back to plain
trimming with a warning instead.

## Links

- `usage/clumpify.md` — choosing `--memory` and `--cores` in practice.
- `internals/trimming-pipeline.md` — the per-read work a worker performs.
- `troubleshooting.md` — the floor warning and the plain-mode fallback.

## Open questions

- None.
