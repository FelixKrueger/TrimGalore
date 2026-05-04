---
title: Hard-trimming
description: Fixed-length hard-clip from either end (`--hardtrim5` / `--hardtrim3`).
---

`--hardtrim5` and `--hardtrim3` are run-and-exit modes that hard-clip every read to a fixed length, then exit. They run before quality and adapter trimming and don't touch report output.

## Hard-trimming to leave bases at the 5' end

`--hardtrim5 INT` clips reads down to their first `INT` bases. Output goes to `{stem}.{INT}bp_5prime.fq(.gz)`. Useful for shortening reads to a uniform length:

```
before:         CCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAAT
--hardtrim5 20: CCTAAGGAAACAAGTACACT
```

## Hard-trimming to leave bases at the 3' end

`--hardtrim3 INT` keeps only the **last** `INT` bases. Output goes to `{stem}.{INT}bp_3prime.fq(.gz)`. Useful for removing a fixed bias at the 5' end (e.g. UMIs or inline barcodes):

```
before:         CCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAAT
--hardtrim3 20:                                                   TTTTTAAGAAAATGGAAAAT
```

## When to use this

Hard-trimming is a sledgehammer: quality is ignored, adapters are ignored, every read is clipped to the same length. It's the right tool for:

- Standardising read lengths for tools that require fixed-length input.
- Removing a known fixed bias at one end (e.g. UMI or inline barcode) before passing data through a second Trim Galore run for adapter and quality trimming.
- Producing equally-sized inputs for benchmarking.

For everything else, prefer the fine-grained `--clip_R1`, `--clip_R2`, `--three_prime_clip_R1`, `--three_prime_clip_R2` flags, which apply on top of normal quality and adapter trimming.

## Compatibility

Both modes accept multiple input files in a single invocation and process them sequentially. Paired-end mode is not required. `--hardtrim*` operates per-file independently.
