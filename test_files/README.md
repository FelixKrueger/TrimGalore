# Test fixtures

Binary fixtures live in `test_files/` because `cargo test` runs from the crate
root and resolves these paths relative to it.

## Quality encoding

Most fixtures here are Phred+33 (Sanger / Illumina 1.8+), which is the default
Trim Galore assumes. Two are **not**, and neither should be used as a generic
input without passing the matching flag:

| Fixture | Encoding | Notes |
|---|---|---|
| `phred64_test.fastq` | **Phred+64** | 4 reads × 34 bp. Quality is 24 × `'h'` (ASCII 104 = Q40) followed by 10 × `'B'` (ASCII 66 = Q2) — a deliberate imitation of the read-segment-quality-control `B`-run tails that real Illumina 1.5 data carries. Added for issue #358: the `--phred64` + `--output-format ubam` tests in `tests/integration_ubam_out.rs` and `tests/integration_clump_only_ubam.rs` assert the writer stores raw Phred `[40 × 24, 2 × 10]`, where the pre-fix writer produced `[71, 33]`. The two distinct values matter — a uniform-quality fixture would only exercise one offset, and would additionally trip FastQC's encoding heuristic (upstream FastQC 0.12.1 behaviour, reproducible on either output format: it keys off the lowest observed quality character, and an all-Q40 record leaves nothing below ASCII 64, so the guess falls to Phred+64 and every mean comes out 31 too low). Tests pass `-q 0` so the low-quality tail survives to the writer; `--length 0` is passed alongside as belt-and-braces but is not required here, since the default length filter of 20 is below the fixture's 34 bp. |
| `truncated.fq.gz` | **Phred+64** | Illumina GA/GAII-era data (`HWUSI-EAS611…#0/1` read names; trailing `B` runs are the Illumina 1.5 read-segment quality-control indicator, `B` = 66 = Q2 at offset 64). Deliberately truncated mid-record — used only as a negative fixture for truncation rejection, so its encoding never matters in practice. Do not cite it as a Phred+33 example. |

Do not assert "all fixtures are Phred+33" anywhere — the two above are
counter-examples.

`truncated.fq.gz` is easy to misclassify because its bytes are plausible under
either offset if you simply assume the default: min byte 66 reads as Q33 at
offset 33 (looks like good data) and as Q2 at offset 64 (the truth). Note that a
minimum-quality-byte heuristic classifies this one **correctly** — 66 ≥ 64, so
it lands on Phred+64. What fails is assuming Phred+33 without looking. The
heuristic's actual weakness runs the other way, and `phred64_test.fastq`
demonstrates it: uniformly high-quality data leaves no byte below 64, so the
minimum tells you nothing and the guess falls back to Phred+64 even when the
data is Phred+33.

## uBAM fixtures

`ubam_test.bam` (SE, 10 reads) and `ubam_paired_test.bam` (PE, 10 pairs,
mate-adjacent interleaved) are committed for the uBAM integration tests in
`tests/integration_ubam.rs` and the format-detection unit tests in
`src/format.rs::tests`.

Recreate from source (requires `samtools`; not a runtime dep):

```sh
gunzip -c test_files/BS-seq_10K_R1.fastq.gz | head -40 > "$TMPDIR/ubam_se.fastq"
gunzip -c test_files/BS-seq_10K_R1.fastq.gz | head -40 > "$TMPDIR/ubam_pe_r1.fastq"
gunzip -c test_files/BS-seq_10K_R2.fastq.gz | head -40 > "$TMPDIR/ubam_pe_r2.fastq"
samtools import -0 "$TMPDIR/ubam_se.fastq" -o test_files/ubam_test.bam
samtools import -1 "$TMPDIR/ubam_pe_r1.fastq" -2 "$TMPDIR/ubam_pe_r2.fastq" -o test_files/ubam_paired_test.bam
```

Verify:
```sh
samtools flagstat test_files/ubam_test.bam        # 10 total, 10 unmapped, 0 secondary
samtools flagstat test_files/ubam_paired_test.bam # 20 total, 20 unmapped, 10 read1 + 10 read2
```

Provenance: SRR24827378 (RRBS). See `BS-seq_10K_R{1,2}.fastq.gz` for the source.

## uBAM output reference fixtures

`ubam_out_se_REFERENCE.bam` (SE) and `ubam_out_pe_REFERENCE.bam` (PE,
interleaved) are committed for the uBAM-output integration tests in
`tests/integration_ubam_out.rs::ubam_out_*matches_reference`.

These reference BAMs are compared via [`assert_ubam_eq`](../tests/integration_ubam_out.rs)
which IGNORES `@PG` lines (which carry `VN:<package-version>` and would
otherwise break the fixture on every routine release bump). Header
sections compared: `@HD` / `@SQ` / `@RG` / `@CO`. Record-stream
comparison: `(name, flags, seq, qual, sorted aux)` tuples per record.

Regenerate when the trim_galore record-level output legitimately changes
(e.g., a default-adapter change, a new aux-tag round-trip rule, a bug fix
in seq/qual handling). Routine `@PG VN:` bumps do NOT require regen.

```sh
cargo build --release
mkdir -p "$TMPDIR/tg_ubam_regen"
./target/release/trim_galore --output-format ubam \
    --output_dir "$TMPDIR/tg_ubam_regen" test_files/ubam_test.bam
cp "$TMPDIR/tg_ubam_regen/ubam_test_trimmed.bam" test_files/ubam_out_se_REFERENCE.bam

./target/release/trim_galore --paired --output-format ubam \
    --output_dir "$TMPDIR/tg_ubam_regen" test_files/ubam_paired_test.bam
cp "$TMPDIR/tg_ubam_regen/ubam_paired_test_val.bam" test_files/ubam_out_pe_REFERENCE.bam
```
