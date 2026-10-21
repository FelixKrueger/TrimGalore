---
title: Poly-G and Poly-A trimming
description: 2-colour instrument poly-G handling (auto-detected) and built-in poly-A trimming for mRNA-seq.
---

v2 ships with two built-in poly-X trimmers: **`--poly_g`** for 2-colour-instrument no-signal G-runs, which the Perl version did not have, and **`--poly_a`** for poly-A-tailed mRNA-seq libraries, which replaces the Perl `--polyA` mode with a different algorithm — see [the migration notes](/reference/migration/).

## Poly-G on 2-colour instruments

On 2-colour Illumina instruments (NextSeq, NovaSeq, NovaSeq X), no signal is encoded as a high-quality **G**. This produces spurious G-runs at the 3' end of reads with little signal. They look like high-quality bases to the basecaller but represent dark cycles, not actual G calls.

`--poly_g` scans the 3' end of **both** Read 1 and Read 2 for these runs and trims them. The artifact is sequencer-induced (overrun → high-quality "G" called for no signal), so it appears at the 3' end of every read regardless of strand orientation — unlike poly-A, which is a template feature and presents at opposite ends of R1 and R2 in the reverse-complement orientation. `--poly_g` is auto-enabled when the data looks 2-colour, based on sequence detection of trailing G-runs.

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

### Stage order

Poly-G is removed before poly-A. The poly-A scan is anchored at the 3' end, so it only reaches a tail that sits at the end of the read; a poly-G run left in front of it would block the scan. On a 2-colour read of the form `[insert][poly-A][poly-G]` — common in 3'seq — this order is what allows both to go.

Two consequences worth knowing:

- With `--no_poly_g`, a poly-A tail behind a G-run is not removed. The opt-out is honoured, and the G-run still blocks the scan.
- In the reverse arrangement, `[insert][poly-G][poly-A]`, the G-run survives: poly-G is now the blocked stage. Dark cycles are terminal by construction, so this is not a real 2-colour read; it leaves sequence in rather than removing too much.

If poly-A tails still survive — with `--no_poly_g`, or behind a G-run the scan does not take — `--2colour 20` in place of `--quality 20` removes the G-run during quality trimming, ahead of both poly stages.

### What `--poly_a` does not remove

The scan is anchored at the 3' end. It removes an A-run that reaches the end of the read and stops at the first base that is not an A, so an A-run with sequence behind it is left in place:

```
insert (40 bp)                           poly-A (25 bp)            + 20 bp
GTCCAACGTTCAGGACTTCAGCGTAACCTGGTTCAGCGTC AAAAAAAAAAAAAAAAAAAAAAAAA CTCGTTGCCTGGTTCAGGCT

--poly_a     →  85 bp, the run is left in place
-a "A{20}"   →  38 bp
```

This is deliberate. Trim Galore does not treat sequence lying beyond an A-run as disposable, because the genome contains many genuine stretches of ten or more consecutive A's. Removing A-runs together with whatever follows them therefore depletes reads originating from those loci, which shows up as missing coverage rather than as anything visible in a QC report. The effect is largest in bisulfite and EM-seq data, where G-to-A conversion creates A-tracts that were never poly-A: on a converted human genome, around 11.6% of 150 bp windows contain a 12-A run ([#180](https://github.com/FelixKrueger/TrimGalore/issues/180)).

### Removing A-runs and everything after them

For 3' counting assays (3'seq, QuantSeq and similar), sequence past the tail carries no information and the more aggressive behaviour is the one you want. Pass a poly-A homopolymer as an adapter:

```bash
trim_galore --quality 20 -a "A{20}" input.fq.gz
```

An adapter match removes the match **and everything 3' of it**, which is what the Perl v0.6.x `--polyA` mode did — it was implemented as a Cutadapt adapter match rather than as a scan. Add your real adapter alongside it (`-a AGATCGGAAGAGC -a "A{20}"`) if you are not relying on auto-detection.

Two things to expect:

- **Prefer `A{20}` over `A{10}`.** The shorter form matches far more genomic A-runs, which is exactly the depletion described above. Measured on `illumina_10K.fastq.gz`, `A{10}` removes 213,308 bp against `A{20}`'s 28,064, and drops 1,358 reads below `--length 20` against 138 — roughly ten times the data loss for the same artifact.
- **A few bases of real insert go with the tail.** Matching is error-tolerant (`-e 0.1`), so a 20 bp A-adapter permits 2 mismatches and the alignment is placed 2 bases early to extend the match. That is why the read above returns 38 bp rather than 40.

## When to use which

| Library type | Recommendation |
|--------------|---------------|
| 2-colour instrument (NextSeq, NovaSeq, NovaSeq X) | Leave `--poly_g` auto-detection on. Add `--2colour 20` if you also want quality-score-based handling. |
| 4-colour instrument (HiSeq, MiSeq) | No poly-G needed. `--poly_g` will not auto-trigger. |
| Poly-A-selected mRNA-seq | Add `--poly_a`. |
| 3'seq / QuantSeq (3' counting) | Use `-a "A{20}"` rather than `--poly_a`, to remove the run and everything behind it. |
| Total RNA / ribodepleted | Skip `--poly_a`. There are no poly-A tails to trim. |
