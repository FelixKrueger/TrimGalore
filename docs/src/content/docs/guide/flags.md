---
title: Flag reference notes
description: Cross-flag interactions and scenarios that need more context than the terse `--help` output.
---

Run `trim_galore --help` for the complete list of options with one-line descriptions and default values. This section focuses on flag interactions and scenarios that benefit from more context than the terse `--help` output.

## Cross-flag interactions

A few combinations are worth knowing about:

- **`--small_rna`** lowers the `--length` default to 18 bp (from 20) and auto-sets `--adapter2` to the Illumina small RNA 5' adapter (`GATCGTCGGACT`) for paired-end data.
- **`--rrbs`** in paired-end directional mode auto-sets `--clip_R2 2` to mask the 2 bp end-repair bias at the start of Read 2, unless the user provides their own `--clip_R2` value. `--non_directional` intentionally skips this auto-clip.
- **`--paired` + `--length`** discards the whole read pair unless *both* reads pass the length cutoff. To rescue the surviving read when only one becomes too short, add `--retain_unpaired`; the per-side cutoff is governed by `--length_1`/`--length_2` (default 35 bp each).
- **`--trim-n`** is suppressed under `--rrbs` (matches Perl v0.6.x; N-trimming interacts poorly with RRBS end-repair masking).
- **`--discard_untrimmed`** keeps only reads where at least one adapter match was found. For paired-end, the pair is discarded only if *neither* read had an adapter.
- **`--poly_g`** is auto-enabled when the data looks like it came from a 2-colour instrument (sequence-based detection on trailing G-runs). Use `--no_poly_g` to force-disable, or `--poly_g` to force-enable. It is independent of `--nextseq` (which is quality-score-based).

## RRBS-specific guidance

For library kits where `--rrbs` should **not** be used (Tecan Ovation RRBS Methyl-Seq, MseI-digested libraries), see [When NOT to use `--rrbs`](/TrimGalore/modes/rrbs/#when-not-to-use---rrbs).

## Adapter specification recap

See [Multiple adapter sequences](/TrimGalore/guide/adapters/#multiple-adapter-sequences) for the three equivalent syntaxes (repeatable `-a`/`-a2`, embedded-string form, and `file:adapters.fa`). Supplementary notes:

- Single-base expansion `-a A{N}` / `-a2 A{N}` repeats the base `N` times, matching Perl v0.6.x syntax.
- Adapter auto-detection runs **per pair** in v2.x (Perl ran it once per invocation), so a shell glob `trim_galore --paired *fastq.gz` correctly handles multiple samples with different library types or 2-colour/4-colour chemistries.
- `--consider_already_trimmed <INT>` suppresses adapter trimming entirely when no auto-detect probe exceeds that count (quality trimming still runs). Useful for feeding already-trimmed data through Trim Galore for QC reporting without over-trimming.

## Perl-era flag forms (still supported)

To keep v0.6.x scripts working unchanged, a small pre-parse hook recognises these Perl-era spellings and rewrites them to the Clap-friendly long-alias forms:

- `-r1` / `-r2` rewrites to `--r1` / `--r2`
- `-a2` rewrites to `--a2`
- `-a " SEQ1 -a SEQ2"` (embedded-string form)

Plus `A{N}` brace expansion in `-a` / `-a2`. See the [migration guide](/TrimGalore/reference/migration/) for the full picture.

## Where to find what

| You want to… | Look here |
|---|---|
| Tune Phred quality trimming | [Quality trimming](/TrimGalore/guide/quality/) |
| Specify or override adapters | [Adapter trimming](/TrimGalore/guide/adapters/) |
| Rescue under-length pairs | [Paired-end data](/TrimGalore/guide/paired-end/) |
| Set length and max-N cutoffs | [Length filtering](/TrimGalore/guide/length/) |
| Read a trimming report | [Trimming reports](/TrimGalore/guide/reports/) |
| Trim RRBS / bisulfite libraries | [RRBS mode](/TrimGalore/modes/rrbs/) and [Bisulfite & RRBS](/TrimGalore/rrbs/guide/) |
| Use a specialty UMI mode | [Mouse Epigenetic Clock](/TrimGalore/modes/clock/) and [IMPLICON](/TrimGalore/modes/implicon/) |
