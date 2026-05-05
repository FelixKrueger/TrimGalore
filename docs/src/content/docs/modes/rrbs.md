---
title: RRBS mode
description: MspI-digested RRBS libraries. Directional and non-directional trimming.
---

Trim galore! also has an `--rrbs` option for DNA material that was digested with the restriction enzyme MspI. In this mode, Trim Galore identifies sequences that were adapter-trimmed and removes another 2 bp from the 3' end of Read 1, and for paired-end libraries also the first 2 bp of Read 2 (which is equally affected by the fill-in procedure). This is to avoid that the filled-in cytosine position close to the second MspI site in a sequence is used for methylation calls. Sequences which were merely trimmed because of poor quality will not be shortened any further.

```bash
# Single-end RRBS
trim_galore --rrbs sample.fq.gz

# Paired-end RRBS (auto-sets --clip_R2 2)
trim_galore --rrbs --paired sample_R1.fq.gz sample_R2.fq.gz
```

For the biology behind RRBS (strand types, fill-in artifacts, and why trimming is needed), see the [Bisulfite & RRBS guide](/rrbs/guide/).

## Non-directional mode

Trim Galore also has a `--non_directional` option, which will screen adapter-trimmed sequences for the presence of either CAA or CGA at the start of sequences and clip off the first 2 bases if found. If CAA or CGA are found at the start, no bases will be trimmed off from the 3' end even if the sequence had some contaminating adapter sequence removed (in this case the sequence read likely originated from either the CTOT or CTOB strand; refer to [the RRBS guide](/rrbs/guide/) for the meaning of CTOT and CTOB strands).

```bash
trim_galore --rrbs --non_directional --paired sample_R1.fq.gz sample_R2.fq.gz
```

`--non_directional` intentionally skips the auto `--clip_R2 2` that directional `--rrbs` applies, since the 5' bias on Read 2 is handled by the start-clip logic instead.

## What `--rrbs` changes vs default mode

| | Default | `--rrbs` |
|---|---------|---------|
| Quality trim | Yes | Yes |
| Adapter trim | Yes | Yes |
| Extra 2 bp clip on adapter-trimmed Read 1 (3' end) | No | Yes |
| Auto `--clip_R2 2` (paired-end) | No | Yes (directional only) |
| `--trim-n` | If passed | Suppressed (matches Perl v0.6.x) |

`--trim-n` is suppressed under `--rrbs` because N-trimming interacts poorly with RRBS end-repair masking. This restores byte-identical output with Perl v0.6.x for users who combine both flags.

## When NOT to use `--rrbs`

Two cases come up regularly.

**Tecan Ovation RRBS Methyl-Seq kit (with or without TrueMethyl oxBS 1-16):** the Tecan kit attaches a varying number of nucleotides (0-3) after each MspI site, so standard `--rrbs` 2 bp 3'-trimming over-trims. Run Trim Galore *without* `--rrbs` and handle the fill-in via Tecan's subsequent diversity-trimming step (see their manual).

**MseI-digested RRBS libraries:** if your DNA was digested with MseI (recognition motif `TTAA`) instead of MspI, you do **not** need `--rrbs` or `--non_directional`. Virtually all reads should start with `TAA`, and the end-repair of TAA-restricted sites does not involve cytosines, so no special treatment is required. Just run Trim Galore in the default (non-RRBS) mode.

## Related flags

| Flag | Purpose |
|------|---------|
| `--rrbs` | Enable RRBS mode (MspI fill-in handling). |
| `--non_directional` | Strip CAA/CGA at start of CTOT/CTOB-derived reads. |
| `--clip_R2 N` | Override the auto 2 bp clip on Read 2 (RRBS paired-end). |
| `--three_prime_clip_R1 N` / `--three_prime_clip_R2 N` | Manual 3' clip on either side. |
