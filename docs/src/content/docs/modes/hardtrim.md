---
title: Hard-trimming
description: Fixed-length hard-clip from either end (`--hardtrim5` / `--hardtrim3`).
---

`--hardtrim5` and `--hardtrim3` are run-and-exit modes that hard-clip every read to a fixed length, then exit. They run before quality and adapter trimming and don't touch report output.

## Hard-trimming to leave bases at the 5'-end

The option `--hardtrim5 INT` allows you to hard-clip sequences from their 3' end. This option processes one or more files (plain FastQ or gzip compressed files) and produces hard-trimmed FastQ files ending in `.{INT}bp_5prime.fq(.gz)`. This is useful when you want to shorten reads to a certain read length. Here is an example:

```
before:         CCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAAT
--hardtrim5 20: CCTAAGGAAACAAGTACACT
```

## Hard-trimming to leave bases at the 3'-end

The option `--hardtrim3 INT` allows you to hard-clip sequences from their 5' end. This option processes one or more files (plain FastQ or gzip compressed files) and produces hard-trimmed FastQ files ending in `.{INT}bp_3prime.fq(.gz)`. We found this quite useful in a number of scenarios where we wanted to remove biased residues from the start of sequences. Here is an example:

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
