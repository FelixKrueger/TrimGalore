---
title: Multiome passthrough
description: Carry a third FASTQ file (cell-barcode / I1) through trimming in lockstep with R1 and R2 via `--passthrough`.
---

`--passthrough <FILE>` adds a third FASTQ input to paired-end trimming that is **not trimmed** — instead, its records are carried through 1:1 to a parallel output file, kept in lockstep with R1/R2 such that any pair dropped by quality, adapter, length, or N-content filtering also drops the matching record from the passthrough output.

The headline use case is **10X Multiome ATAC-seq**, where the modified ATAC adapter is opaque to Cellranger Arc's auto-detection. Trim R1/R2 against your own `-a`/`-a2` while the cell-barcode read (I1 / I2 / R3 — whichever your library design names it) stays aligned to the surviving fragments. One invocation, no awk pipelines on the side.

## Usage

```bash
trim_galore --paired \
            --passthrough sample_I1.fastq.gz \
            sample_R1.fastq.gz \
            sample_R2.fastq.gz
```

Output files:

| File | Description |
|---|---|
| `sample_R1_val_1.fq.gz` | Trimmed R1 (standard `--paired` output) |
| `sample_R2_val_2.fq.gz` | Trimmed R2 (standard `--paired` output) |
| `sample_I1_passthrough.fq.gz` | Passthrough output — same record count and order as R1/R2 outputs |
| `sample_R{1,2}.fastq.gz_trimming_report.{txt,json}` | Trimming reports — R2's text report contains a new `=== Passthrough file ===` block, both JSON reports gain a `"passthrough"` object |

With `--basename foo` set, all three outputs use that basename uniformly: `foo_val_1.fq.gz`, `foo_val_2.fq.gz`, `foo_passthrough.fq.gz`.

## How the sync check works

Every record, Trim Galore extracts an ID prefix from R1, R2, and the passthrough read and compares them three-way. The prefix is *everything after the leading `@`, before the first whitespace, with a trailing `/1` / `/2` / `/3` stripped*. That covers both modern Illumina (`@HEADER 1:N:0:CGATCG` → prefix `HEADER`) and legacy SRA/ENA (`@read/1` `@read/2` `@read/3` → prefix `read`).

If any pair-wise comparison disagrees, the run fails loudly with a row-numbered error:

```
Read ID mismatch at record 4137: R1='SRR123.4137', R2='SRR123.4137',
passthrough='SRR123.4138'. Files are out of sync.
```

Mid-stream truncation in any of the three streams is caught the same way:

```
--passthrough file is truncated — R1/R2 have more reads than passthrough
at record 5001
```

The opposite case (passthrough longer than R1/R2) fails with a matching message. Files of unequal length never produce silently misaligned barcodes — that was the load-bearing concern that drove the per-record sync check.

## What's not supported in v1

`--passthrough` is rejected at startup in combination with:

- **`--retain_unpaired`** — passthrough requires strict pair semantics in v1; a half-rescued mate has no clean interpretation for the index file.
- **`--clumpify`** — the clumpy reordering shuffles records by minimizer, which breaks the lockstep contract.
- **Multi-pair input** — exactly one R1/R2 pair per invocation; for multiple samples, run Trim Galore once per pair.
- **Specialty modes** — `--clock`, `--implicon`, `--hardtrim5/3`, `--demux` each own their own input arity and output naming.

`--passthrough` is also rejected without `--paired`. Each rejection comes with a precise error message.

## FastQC on the passthrough output

When `--fastqc` is set, FastQC runs on all three outputs, including the passthrough. **The passthrough report will look noisier than R1/R2:** cell-barcode reads are 16–28 bp of uniformly-structured sequence by design (positional base bias from the barcode whitelist + natural duplication from multi-fragment cells), so FastQC's *Per base sequence content* and *Sequence Duplication Levels* modules typically FAIL even on perfectly good data. That's the data type, not a defect. The remaining modules (per-base quality, adapter content, N-content) carry real signal.

## Implementation notes

- **Untouched bytes.** The passthrough record's sequence and quality are written verbatim. The third line (`+...` description) and line endings are canonicalised to bare `+` and `\n` — the standard canonicalisation R1/R2 already go through. The id / seq / qual fields are byte-identical to the input.
- **Single-pass, in-memory lockstep.** The third stream is read in parallel with R1/R2 in the same reader thread (parallel path) or the same loop iteration (serial path); the writes happen in the same flush iteration as R1/R2 within each batch. No temp files.
- **Mid-stream errors leave partial output on disk.** If a sync error fires at record 5,000 of 10,000, files already-written to disk are not rolled back — re-run the input after fixing the source. Trim Galore prints a clear error and exits with non-zero status.
- **Output gzip is uniform** across all three outputs and follows R1's input compression (plain R1 in → plain `_passthrough.fq` out; gzipped R1 in → gzipped `_passthrough.fq.gz` out), regardless of the passthrough input's own extension.

## Worked example

Generate a synthetic I1 file from one of the bundled BS-seq fixtures and trim everything in one pass:

```bash
# Mint a synthetic 16 bp barcode read aligned to BS-seq_10K_R1's headers.
gunzip -c test_files/BS-seq_10K_R1.fastq.gz | awk 'NR % 4 == 1 {
  print $1
  print "AAAACCCCGGGGTTTT"
  print "+"
  print "IIIIIIIIIIIIIIII"
}' | gzip > /tmp/BS-seq_10K_I1.fastq.gz

trim_galore --paired \
            --passthrough /tmp/BS-seq_10K_I1.fastq.gz \
            --output_dir /tmp \
            test_files/BS-seq_10K_R1.fastq.gz \
            test_files/BS-seq_10K_R2.fastq.gz

# All three outputs should report the same record count.
for f in /tmp/BS-seq_10K_R1_val_1.fq.gz \
         /tmp/BS-seq_10K_R2_val_2.fq.gz \
         /tmp/BS-seq_10K_I1_passthrough.fq.gz; do
  printf "  %-50s %s records\n" "$(basename $f)" \
    "$(gunzip -c "$f" | awk 'NR%4==1' | wc -l | tr -d ' ')"
done

# The new block in the R2 text report:
grep -A6 "Passthrough file" /tmp/BS-seq_10K_R2.fastq.gz_trimming_report.txt
```

Out of 10,000 input pairs you should see 9,996 surviving (4 dropped at default settings), and the same count in all three outputs.
