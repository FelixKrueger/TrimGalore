---
title: Library presets
description: One kit name instead of four clip flags, with --library.
---

Most modern methylation library preps need a fixed number of bases clipped off each read end to mask chemistry artefacts. Getting that wrong biases methylation calls in a way that is easy to miss, and the numbers historically lived in vendor PDFs rather than in the tool.

`--library <preset>` replaces the four clip flags with the kit name:

```bash
trim_galore --paired --library emseq sample_R1.fq.gz sample_R2.fq.gz
```

is exactly:

```bash
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
  sample_R1.fq.gz sample_R2.fq.gz
```

## Presets

| `--library` | `--clip_R1` | `--clip_R2` | `--three_prime_clip_R1` | `--three_prime_clip_R2` |
|---|---|---|---|---|
| `emseq` | 10 | 10 | 10 | 10 |
| `accel` (aliases `swift`, `xgen`) | 10 | 15 | 10 | 10 |
| `zymo` | 10 | 10 | 10 | 10 |
| `scbs` (alias `single_cell`) | 6 | 6 | 6 | 6 |
| `pbat` | 8 | 8 | 8 | 8 |

`accel` carries aliases because the kit was resold twice: Accel-NGS Methyl-seq, then Swift, now IDT xGen Methyl-seq. The values match [nf-core/methylseq](https://nf-co.re/methylseq)'s presets of the same names, so trimming is identical whether you drive Trim Galore through that pipeline or directly.

## What a preset does *not* set

A preset is 5' and 3' end clipping, and nothing else:

- **No adapter sequence.** Auto-detection runs per pair and would generally pick better than a preset pinning one sequence.
- **No `--length`.** That is a filtering choice, not a property of the chemistry.
- **No mode changes.** `--rrbs`, `--clock` and `--implicon` are not presets and stay separate flags, because they change trimming behaviour rather than expanding into four numbers. Tecan/NuGEN Ovation RRBS is absent for the same reason: it needs diversity trimming, a different algorithm, and it must run *without* `--rrbs`. See [When NOT to use `--rrbs`](/modes/rrbs/#when-not-to-use---rrbs).
- **No effect in the specialty modes.** `--hardtrim5`, `--hardtrim3`, `--clock` and `--implicon` do not honour the clip flags, so `--library` is refused alongside them rather than accepted and silently ignored. `--clump_only`, which does no trimming at all, refuses it too.

## Overrides

An explicit clip flag wins over the preset value it displaces, and both numbers appear in the log:

```bash
trim_galore --library emseq --clip_R1 12 sample.fq.gz
```

```
Library preset 'emseq' selected: --clip_R1 12 --three_prime_clip_R1 10
--clip_R1 12 was given on the command line and overrides the emseq preset value 10
```

The example above is single-end, so only the Read 1 values are listed: `--clip_R2` and
`--three_prime_clip_R2` are reported for paired-end runs, where they apply.

Refusing the combination would send you back to writing all four values by hand, which is what the flag exists to remove.

## Reporting and versioning

Both trimming reports record the preset name and the clipping it resolved to. The text report lists
the values in force for the run, which for this single-end example is Read 1 only:

```
Library preset: emseq
Clipping in force: --clip_R1 12 --three_prime_clip_R1 10
--clip_R1 12 was given on the command line and overrides the emseq preset value 10
```

and the JSON report, under `parameters.library` (`null` when no preset was used):

```json
"library": {
  "preset": "emseq",
  "clip_r1": 12,
  "clip_r2": 10,
  "three_prime_clip_r1": 10,
  "three_prime_clip_r2": 10,
  "overrides": [{"flag": "--clip_R1", "preset_value": 10, "user_value": 12}]
}
```

The JSON block carries all four resolved values whatever the run type, so the whole configuration is
readable in one place; the top-level `mode` key (`single-end` or `paired-end`) says which of them the
reads received. That is why the same run's text and JSON reports differ above.

Preset values may change between releases when vendor guidance moves; every change gets a `CHANGELOG.md` entry. Because the expanded numbers are in the report, an old run stays reproducible from its own record.
