---
title: Mouse Epigenetic Clock trimming
description: UMI extraction for dual-UMI RRBS libraries used in the Mouse Epigenetic Clock.
---

The option `--clock` trims reads in a specific way that is currently used for the Mouse Epigenetic Clock (see here: [Multi-tissue DNA methylation age predictor in mouse, Stubbs et al., Genome Biology, 2017 18:68](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1203-5)). Following the trimming, Trim Galore exits.

In its current implementation, the dual-UMI RRBS reads come in the following format:

```
Read 1  5' UUUUUUUU CAGTA FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF TACTG UUUUUUUU 3'
Read 2  3' UUUUUUUU GTCAT FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF ATGAC UUUUUUUU 5'
```

Where UUUUUUUU is a random 8-mer unique molecular identifier (UMI), CAGTA is a constant region, and FFFFFFF... is the actual RRBS-Fragment to be sequenced. The UMIs for Read 1 (R1) and Read 2 (R2), as well as the fixed sequences (F1 or F2), are written into the read ID and removed from the actual sequence. Here is an example:

```
R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT
    ATCTAGTTCAGTACGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT
    CAATTTTGCAGTACAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA

R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA
```

Following clock trimming, the resulting files (`.clock_UMI.R1.fq(.gz)` and `.clock_UMI.R2.fq(.gz)`) should be adapter- and quality-trimmed with a second Trim Galore run. Even though the data is technically RRBS, it doesn't require the `--rrbs` option. Instead the reads need to be trimmed by 15 bp from their 3' end to get rid of potential UMI and fixed sequences. For a single sample:

```bash
trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 sample.clock_UMI.R1.fq.gz sample.clock_UMI.R2.fq.gz
```

For multiple samples, `--paired` requires input files in pairwise (R1, R2, R1, R2, …) order. Don't use `*R1.fq.gz *R2.fq.gz` globs which produce all-R1s-then-all-R2s. Either invoke once per sample, or interleave explicitly.

Following this, reads should be aligned with Bismark and deduplicated with UmiBam in `--dual_index` mode (see here: <https://github.com/FelixKrueger/Umi-Grinder>). UmiBam recognises the UMIs within this pattern: `R1:(`**ATCTAGTT**`):R2:(`**CAATTTTG**`):` as UMI R1 = **ATCTAGTT** and UMI R2 = **CAATTTTG**.
