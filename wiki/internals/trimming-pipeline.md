---
type: internals
title: Order of the trimming steps
description: The per-read step order in trim_read, which module owns each step, and where filtering and auto-detection attach.
tags: [trimming, pipeline, order]
page_id: trimming-pipeline
subject_status: current
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/trimmer.rs
    commit: bf0be75
  - repo: TrimGalore
    path: src/quality.rs
    commit: 95c8154
  - repo: TrimGalore
    path: src/alignment.rs
    commit: 2e00f32
  - repo: TrimGalore
    path: src/adapter.rs
    commit: d2ea293
  - repo: TrimGalore
    path: src/filters.rs
    commit: 32fb357
---

## Summary

One read passes through `trimmer::trim_read`, which applies every trimming
step in a fixed order and returns what it did as a `TrimResult`. Filtering is
not part of that function: it runs afterwards, on the trimmed read. Adapter
auto-detection runs earlier still, before any read is trimmed.

## Step order in trim_read

The order matters, because each step sees what the previous one left.

| Step | Owner | What it does |
|---|---|---|
| Quality trim | `quality::quality_trim_3prime`, or `quality::quality_trim_3prime_nextseq` under `--nextseq` | 3' quality trimming. The NextSeq variant treats G as quality zero, so a 2-colour poly-G tail is removed here, before adapter matching. |
| Adapter trim | `alignment::find_3prime_adapter`, looped in `trim_read` | Repeats up to `--times` rounds. Each round scores every configured adapter and takes the match starting earliest in the read. R2 uses the `-a2` set when one is configured. |
| RRBS | `trim_read` | Non-directional input starting `CAA` or `CGA` loses 2 bp from the 5' end; otherwise, and on directional R1, 2 bp come off the 3' end when an adapter was found. Directional R2 is handled by the auto-set `--clip_r2 2` instead. |
| Poly-A / poly-T | `quality::poly_a_trim_index` | R1 trims a 3' poly-A; R2 trims a 5' poly-T, because the artifact is template-side and flips with the strand. |
| Poly-G | `quality::homopolymer_trim_index` | 3' on both R1 and R2. The 2-colour artifact is sequencer-induced, so it does not flip with the strand the way poly-A does. |
| N trim | `FastqRecord::trim_ns` | Both ends, under `--trim-n`. Suppressed under `--rrbs`. |
| Fixed clipping | `trim_read` | `--clip_r1` / `--clip_r2` at the 5' end, `--three_prime_clip_r1` / `--three_prime_clip_r2` at the 3' end. Under `--rename` the removed bases are appended to the read ID. |

Between the adapter step and RRBS, `trim_read` records the read length as
`bp_after_cutadapt`. That value feeds the "Total written (filtered)" figure in
the report, and it is captured there so the figure matches Perl v0.6.x on RRBS
samples, where later steps would otherwise change it.

## Where filtering attaches

`filters::filter_single_end` and `filters::filter_paired_end` run after
`trim_read`, from the orchestrating functions in `trimmer.rs`. They decide
whether the trimmed read survives — length, maximum length, maximum N count,
and for paired input whether one mate can be rescued as unpaired. A read is
therefore always trimmed before it is judged, never filtered out first.

## Where auto-detection attaches

Adapter auto-detection is a separate pass in `adapter.rs` over the start of
the file, capped at one million reads, and it finishes before trimming begins.
In paired mode it runs per pair, which is a deliberate difference from Perl
v0.6.x, where detection ran once on the first file only.

The same pass counts reads carrying a long trailing poly-G run, and that count
is what decides whether poly-G trimming switches itself on. So the decision is
made on untrimmed reads while the trimming itself happens at the poly-G step
above, on reads that adapter trimming has already shortened.

## Links

- `reference/cli-flags.md` — the flags that switch each step on, with defaults.
- `internals/parallel-and-clumping.md` — how this per-read path is driven under `--cores`.
- `reference/outputs.md` — the report fields the steps feed.

## Open questions

- None.
