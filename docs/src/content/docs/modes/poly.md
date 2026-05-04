---
title: Poly-G and Poly-A trimming
description: 2-colour instrument poly-G handling (auto-detected) and built-in poly-A trimming for mRNA-seq.
---

v2 ships with two built-in poly-X trimmers that the Perl version did not have: **`--poly_g`** for 2-colour-instrument no-signal G-runs, and **`--poly_a`** for poly-A-tailed mRNA-seq libraries.

## Poly-G on 2-colour instruments

On 2-colour Illumina instruments (NextSeq, NovaSeq, NovaSeq X), no signal is encoded as a high-quality **G**. This produces spurious G-runs at the 3' end of reads with little signal. They look like high-quality bases to the basecaller but represent dark cycles, not actual G calls.

`--poly_g` scans the 3' end of Read 1 (and the 5' end of Read 2, where the reverse complement reads as poly-C) for these runs and trims them. It is auto-enabled when the data looks 2-colour, based on sequence detection of trailing G-runs.

```bash
# Auto-detection runs by default. No flag needed.
trim_galore input.fq.gz

# Force-enable
trim_galore --poly_g input.fq.gz

# Force-disable
trim_galore --no_poly_g input.fq.gz
```

`--poly_g` is independent of `--nextseq` and `--2colour`:

- `--2colour N` (alias `--nextseq N`) is quality-score based and replaces `-q`. It treats trailing Gs as low-quality during quality trimming.
- `--poly_g` is sequence-based and runs after quality trimming. It looks for actual G-runs, regardless of their quality scores.

You can use them together. They target overlapping artifacts.

## Poly-A for mRNA-seq libraries

`--poly_a` trims poly-A tails from the 3' end of reads. Use this for **mRNA-seq and poly-A-selected RNA-seq libraries**, where the poly-A tail can survive into the read.

```bash
trim_galore --poly_a input.fq.gz
```

In paired-end mode, the corresponding **poly-T** is removed from the 5' end of the mate read.

`--poly_a` is opt-in. It is not auto-detected. Pass it when you know your library is poly-A selected.

## When to use which

| Library type | Recommendation |
|--------------|---------------|
| 2-colour instrument (NextSeq, NovaSeq, NovaSeq X) | Leave `--poly_g` auto-detection on. Add `--2colour 20` if you also want quality-score-based handling. |
| 4-colour instrument (HiSeq, MiSeq) | No poly-G needed. `--poly_g` will not auto-trigger. |
| Poly-A-selected mRNA-seq | Add `--poly_a`. |
| Total RNA / ribodepleted | Skip `--poly_a`. There are no poly-A tails to trim. |
