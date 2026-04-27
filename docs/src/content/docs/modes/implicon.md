---
title: IMPLICON UMI transfer
description: UMI extraction for IMPLICON libraries. The UMI lives at the 5' end of Read 2.
---

The option `--implicon` is a similar run-and-exit specialty mode for paired-end libraries where the UMI lives at the 5' end of Read 2 only. It transfers the first `N` bases of R2 into the read ID of both reads (as `:BARCODE`), then clips R2 by `N` bases. Read 1 is otherwise untouched. In Perl v0.6.x the UMI length was hardcoded at 8 bp; v2.x takes an optional value (`--implicon=12`, for example) with a default of 8.

Output filenames follow the pattern `{stem}_{N}bp_UMI_R1.fastq(.gz)` / `{stem}_{N}bp_UMI_R2.fastq(.gz)`. Like `--clock`, Trim Galore exits after the transfer; the resulting files are then run through a second Trim Galore invocation for adapter and quality trimming.

## Workflow

```bash
# Step 1: IMPLICON UMI transfer
trim_galore --implicon --paired sample_R1.fq.gz sample_R2.fq.gz

# Step 2: adapter / quality trim
trim_galore --paired sample_8bp_UMI_R1.fq.gz sample_8bp_UMI_R2.fq.gz
```

The downstream alignment and deduplication tooling parses the UMI back out of the read ID.

## Multi-pair input

Since v2.1.0-beta.4, `--implicon` accepts any even number of inputs as consecutive R1/R2 pairs:

```bash
trim_galore --implicon A_R1.fq.gz A_R2.fq.gz B_R1.fq.gz B_R2.fq.gz
```

Each pair gets a header (`=== IMPLICON pair N of M ===`) and the same output-collision pre-flight (case-insensitive, full-path) that `--paired` runs. Pairwise (R1, R2, R1, R2, …) order is required; don't pass `*R1.fq.gz *R2.fq.gz` globs (those produce all-R1s-then-all-R2s).
