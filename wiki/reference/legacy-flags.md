---
type: reference
title: Deprecated and no-op flags
description: The v0.6.x flags Trim Galore still parses, and why each one no longer does anything.
tags: [cli, legacy, deprecated]
page_id: legacy-flags
subject_status: legacy
updated: 2026-10-13
sources:
  - repo: TrimGalore
    path: src/cli.rs
    commit: bf28cbd
---

## Summary

These flags are accepted so that command lines written for the Perl v0.6.x
releases keep parsing. Each one prints a deprecation note and changes nothing
about the run. The Perl version itself is deprecated (developer's statement).

| Flag | Why it does nothing |
|---|---|
| `--gzip` | Output is gzipped by default; `--dont_gzip` is the lever now. |
| `--path_to_cutadapt`, `--cutadapt_args` | Cutadapt is built in, not called out to. |
| `--report` | Reports are written by default. |
| `--keep` | The single-pass architecture has no quality-trim intermediate to keep. |
| `--suppress_warn` | Only essential progress output is emitted; redirect the stream if you want less. |
| `--hulu` | Easter egg. |

`--phred33` is *not* in this group: it names the default encoding and is still
meaningful as an explicit opt-in, with `--phred64` as the alternative.

## Links

- `reference/cli-flags.md` — the flags that do have an effect.
- `overview.md` — the relationship to the Perl v0.6.x releases.

## Open questions

- None.
