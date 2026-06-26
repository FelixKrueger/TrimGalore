# Test fixtures

Binary fixtures live in `test_files/` because `cargo test` runs from the crate
root and resolves these paths relative to it.

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
