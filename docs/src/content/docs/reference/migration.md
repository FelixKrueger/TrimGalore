---
title: Migrating from v0.6.x
description: What changed between the Perl Trim Galore (v0.6.x) and the Rust Oxidized Edition (v2).
---

Trim Galore v2 is a Rust rewrite of the Perl wrapper, built as a drop-in replacement for v0.6.x. Outputs match the Perl original for the core flag set, verified end-to-end via the nf-core/rnaseq integration matrix. The rewrite adds capabilities the Perl version did not have.

## What stays the same

- Command name. Still `trim_galore`.
- Output filenames. `*_trimmed.fq.gz`, `*_val_1.fq.gz` / `*_val_2.fq.gz`, `*_unpaired_*.fq.gz`, `*_trimming_report.txt`. Drop-in for any pipeline that consumes these names.
- Trimming-report text format. Cutadapt-compatible header is retained for MultiQC backwards compatibility.
- CLI flags. `--paired`, `--rrbs`, `--non_directional`, `--length`, `--clip_R1`, `--three_prime_clip_R2`, and friends all behave the same way.
- Per-flag output bytes. Verified byte-identical to Perl for the core feature set across the nf-core/rnaseq matrix.

## What's new in v2

| Feature | Perl v0.6.x | Rust v2 |
|---|---|---|
| Architecture | Wrapper around Cutadapt and pigz | Single process, single binary |
| Dependencies | Perl, Python, Cutadapt, Java, igzip, pigz, external FastQC | None (FastQC bundled via fastqc-rust) |
| Paired-end | Two sequential Cutadapt runs + pair validation | Single pass |
| Threads | Up to ~3N+3 (Cutadapt + pigz + pigz/igzip) | Exactly N+4 |
| Adapter auto-detection | Once per invocation | Per pair (handles mixed-library globs) |
| BGI/DNBSEQ adapter | `--bgiseq` only | Auto-detected too |
| Multi-adapter syntax | Embedded-string only | Repeatable `-a` / `-a2`, plus embedded and FASTA |
| Poly-G handling | None | Auto-detected for 2-colour instruments |
| Poly-A trimming | None | Built-in `--poly_a` |
| `--implicon` UMI length | Hardcoded 8 bp | Configurable (default 8) |
| JSON trimming report | None | Included alongside text report |
| `--version` provenance | Just version | Version + git hash + target + build timestamp |
| Reproducible builds | not available | `SOURCE_DATE_EPOCH` produces a bit-identical binary |

## Backwards-compatibility shims

A small pre-parse hook recognises Perl-era flag spellings and rewrites them to the Clap-friendly long-alias forms, so existing scripts keep working without changes:

| Perl form | v2 alias |
|-----------|----------|
| `-r1 N` | `--r1 N` |
| `-r2 N` | `--r2 N` |
| `-a2 SEQ` | `--a2 SEQ` |
| `-a " SEQ -a SEQ"` (embedded) | Recognised directly |
| `-a A{N}` brace expansion | Recognised directly |

## What's not in v2 (intentional)

- Colorspace input. Rejected with a clear error message, matching v0.6.x behaviour. Colorspace data has not been generated in years.
- The per-flag `Cutadapt version` and `Python version` lines in the parameter summary. Cutadapt and Python are not subprocesses in v2, so those lines are dropped. The MultiQC parser already handles their absence.

## Note on `-a 'A{N}'` for poly-A trimming

For poly-A trimming, prefer the dedicated `--poly_a` flag — it has proper paired-end semantics (3' end on R1, 5' end on R2) and is the supported v2 path. The `-a 'A{N}'` form (a v0.6.x workaround that pre-dates `--poly_a`) still works and the brace expansion still produces `AAAA…` byte-for-byte, but the alignment DP may produce slightly different matches than Perl v0.6.x on highly-repetitive adapter patterns: every position in a poly-A read has potential matches, so the choice of "where the adapter starts" reduces to a tie-break that Cutadapt and the Rust port resolve differently. Both are valid global-best alignments under their respective tie-break rules — neither is wrong by spec — but the byte-identity guarantee in the table above does not extend to single-base repeating adapters.

## Migration checklist

To move an existing pipeline from v0.6.x to v2:

1. **Replace the install.** `cargo install trim-galore` or pull the Docker image. Drop the conda environment with Cutadapt and pigz.
2. **Threading.** Drop the `task.cpus - 4` Nextflow shim if you used it. Pass `--cores task.cpus` directly.
3. **Auto-detection scripts.** If you parsed the auto-detect line, note that v2 may now report BGI/DNBSEQ for BGI data instead of falling back to Illumina.
4. **MultiQC.** No changes needed. The report format is unchanged.
5. **JSON report.** Optional new artifact. Consume it for cleaner programmatic parsing if your dashboard wants structured data.

If you hit a regression vs v0.6.x output, open an issue with the `--version` output and a minimal reproducer. Byte-identity for the core flag set is a hard guarantee in v2.
