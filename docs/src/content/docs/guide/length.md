---
title: Length filtering
description: Removing short reads after quality and adapter trimming.
---

Lastly, since quality and/or adapter trimming may result in very short sequences (sometimes as short as 0 bp), Trim Galore can filter trimmed reads based on their sequence length (default: 20 bp). This is to reduce the size of the output file and to avoid crashes of alignment programs which require sequences with a certain minimum length.

## Single-end

Reads below the `--length` cutoff are removed before being written to the trimmed output file. The number of removed sequences is logged in the trimming report:

```
RUN STATISTICS FOR INPUT FILE: SE.fastq.gz
=============================================
1000000 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp: 552 (0.1%)
```

## Paired-end

Note that it is not recommended to remove too-short sequences if the analysed FastQ file is one of a pair of paired-end files, since this confuses the sequence-by-sequence order of paired-end reads which is again required by many aligners. For paired-end files, Trim Galore has an option `--paired` which runs a paired-end validation on both trimmed `_1` and `_2` FastQ files once the trimming has completed. This step removes entire read pairs if at least one of the two sequences became shorter than a certain threshold. If only one of the two reads is longer than the set threshold, e.g. when one read has very poor qualities throughout, this singleton read can be written out to unpaired files (see option `--retain_unpaired`) which may be aligned in a single-end manner.

See [Paired-end data](/TrimGalore/guide/paired-end/) for the full validation flow.

## Other filters

- `--max_length INT`. Drop reads above this length.
- `--max_n INT`. Drop reads containing more than this many `N` bases.
- `--trim-n`. Trim trailing `N`s from both ends before length filtering.

These filters run after quality and adapter trimming, so they apply to the final trimmed read length.

## Small RNA libraries

`--small_rna` automatically lowers `--length` to **18 bp** to match the typical small-RNA library design. Pass `--length` explicitly to override.

## Related flags

| Flag | Purpose |
|------|---------|
| `--length INT` | Minimum read length (default 20 bp; 18 bp under `--small_rna`). |
| `--max_length INT` | Drop reads longer than this. |
| `--max_n INT` | Drop reads with more than `N` `N`-bases. |
| `--trim-n` | Trim trailing `N`s from both ends. Suppressed under `--rrbs`. |
| `--retain_unpaired` | Write singletons to `*_unpaired_*.fq.gz` instead of discarding the pair. |
| `--length_1 INT` / `--length_2 INT` | Per-side length cutoff for `--retain_unpaired` (default 35 bp each). |
