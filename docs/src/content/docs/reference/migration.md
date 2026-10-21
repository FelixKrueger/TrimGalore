---
title: Migrating from v0.6.x
description: What changed between the Perl Trim Galore (v0.6.x) and the Rust rewrite (v2).
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
| Poly-A trimming | `--polyA`, a separate second pass | Built-in `--poly_a` stage |
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

## Poly-A trimming: `--polyA` in v0.6.x vs `--poly_a` in v2

`--poly_a` is not a renamed `--polyA`. The two differ in what they do and in how they are invoked, so expect different output rather than identical output.

**In v0.6.x**, `--polyA` was an exclusive mode that ran a second Cutadapt pass with a poly-A homopolymer as the adapter. Three consequences:

- Its Cutadapt call carried no quality cutoff, so it expected input that had already been adapter- and quality-trimmed. The documented workflow was two runs: trim adapters and qualities first, then `--polyA`. Length filtering still applied.
- The adapter was `A{20}`, with `T{150}` as `-a2` on Read 2. Auto-detection decided only the orientation — whether a sample was poly-A or poly-T — by counting reads containing a 10-mer of each; the chosen direction was then expanded to the fixed 20 and 150 bp sequences.
- It rewrote read IDs. Whitespace became underscores on every read, and reads that lost bases also gained a `32:A:` prefix and a `_PolyA:32` suffix recording how many. v2 does neither, so a pipeline that greps read IDs for `PolyA:`, or that relies on the underscores, finds nothing.

**In v2**, `--poly_a` is a stage inside the single pass, applied after adapter removal and poly-G trimming, and it uses a 3'-anchored homopolymer scan rather than an adapter alignment. A 3' adapter match removes everything downstream of the match, which is why the v0.6.x mode also cleared any bases sitting behind the tail; the v2 scan does not, which is why poly-G is removed first. See [Poly-G and Poly-A trimming](/modes/poly/) for the stage order.

### `-a 'A{N}'` still works

The `-a 'A{N}'` form parses and the brace expansion produces `AAAA…` byte-for-byte, so it is the way to get the v0.6.x behaviour back — see [removing A-runs and everything after them](/modes/poly/#removing-a-runs-and-everything-after-them), which also covers why `A{20}` is the right length to pass. The alignment DP may produce slightly different matches than Perl v0.6.x on highly-repetitive adapter patterns: every position in a poly-A read has potential matches, so the choice of "where the adapter starts" reduces to a tie-break that Cutadapt and the Rust port resolve differently. Both are valid global-best alignments under their respective tie-break rules — neither is wrong by spec — but the byte-identity guarantee in the table above does not extend to single-base repeating adapters.

## Migration checklist

To move an existing pipeline from v0.6.x to v2:

1. **Replace the install.** `cargo install trim-galore` or pull the Docker image. Drop the conda environment with Cutadapt and pigz.
2. **Threading.** Drop the `task.cpus - 4` Nextflow shim if you used it. Pass `--cores task.cpus` directly.
3. **Auto-detection scripts.** If you parsed the auto-detect line, note that v2 may now report BGI/DNBSEQ for BGI data instead of falling back to Illumina.
4. **MultiQC.** No changes needed. The report format is unchanged.
5. **JSON report.** Optional new artifact. Consume it for cleaner programmatic parsing if your dashboard wants structured data.

If you hit a regression vs v0.6.x output, open an issue with the `--version` output and a minimal reproducer. Byte-identity for the core flag set is a hard guarantee in v2.
