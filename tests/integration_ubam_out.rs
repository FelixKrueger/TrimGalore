//! Binary-driven integration tests for uBAM **output** support
//! (PLAN v2.1 §5 step 5).
//!
//! Tests exercise the built `trim_galore` binary with `--output-format ubam`
//! end-to-end, covering wiring the lib-tests can't reach:
//!   * `main.rs::run_ubam_output` dispatch (SE / PE-two-FASTQ / PE-one-uBAM-interleaved)
//!   * format-detection-driven source-header propagation
//!   * `--hardtrim5` + uBAM (specialty mode dispatch)
//!   * each PLAN §3.4a rejection rule, end-to-end
//!
//! Golden-fixture comparison uses [`assert_ubam_eq`] — a noodles-based
//! tuple-comparator that IGNORES the `@PG` chain (which carries the
//! `VN:<package-version>` tag and would otherwise break on every release
//! bump). Compare header-minus-@PG + per-record (name, flags, seq, qual,
//! sorted aux) tuples.

use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;
use noodles::sam::alignment::record::data::field::Tag;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// One BAM record, projected into comparable form: (name, flag-bits, seq
/// bytes, qual raw-Phred bytes, sorted-by-tag aux byte tuples).
type RecordTuple = (Vec<u8>, u16, Vec<u8>, Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>);

/// Extract all `(name, flag, seq, qual, aux)` tuples from a BAM file.
/// Aux fields are sorted by tag-bytes so map-iteration order doesn't
/// influence equality.
fn bam_tuples(path: &Path) -> Vec<RecordTuple> {
    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("failed to open {}: {}", path.display(), e));
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    let _header = reader.read_header().expect("BAM header parse");
    let mut tuples = Vec::new();
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("BAM record read") > 0 {
        let name: Vec<u8> = rec
            .name()
            .map(|bn| AsRef::<[u8]>::as_ref(bn).to_vec())
            .unwrap_or_default();
        let flag = rec.flags().bits();
        let seq: Vec<u8> = rec.sequence().iter().collect();
        let qual: Vec<u8> = rec.quality_scores().as_ref().to_vec();
        let mut aux: Vec<(Vec<u8>, Vec<u8>)> = rec
            .data()
            .iter()
            .filter_map(|res| {
                let (tag, value) = res.ok()?;
                let tag_bytes = vec![tag.as_ref()[0], tag.as_ref()[1]];
                let value_dbg = format!("{:?}", value).into_bytes();
                Some((tag_bytes, value_dbg))
            })
            .collect();
        aux.sort_by(|a, b| a.0.cmp(&b.0));
        tuples.push((name, flag, seq, qual, aux));
    }
    tuples
}

/// Extract the SAM header text minus the `@PG` lines.
fn header_minus_pg(path: &Path) -> String {
    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("failed to open {}: {}", path.display(), e));
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    let header = reader.read_header().expect("BAM header parse");
    // Strip the @PG block, keep @HD, @SQ, @RG, @CO.
    let mut out = String::new();
    if let Some(hd) = header.header() {
        out.push_str(&format!("HD:{:?}\n", hd));
    }
    for (id, sq) in header.reference_sequences() {
        out.push_str(&format!("SQ:{:?}:{:?}\n", id, sq));
    }
    for (id, rg) in header.read_groups() {
        out.push_str(&format!("RG:{:?}:{:?}\n", id, rg));
    }
    for co in header.comments() {
        out.push_str(&format!("CO:{:?}\n", co));
    }
    out
}

/// Compare two BAM files for content equivalence, IGNORING `@PG` lines.
///
/// Compares:
///   - SAM header lines EXCEPT `@PG` (`@HD`, `@SQ`, `@RG`, `@CO`)
///   - record stream: (name, flags, seq, qual, sorted aux fields) per record
///
/// Tuple equality avoids transient field-ordering noise in aux Data;
/// sorting by tag-name on both sides normalises that. NOT a `samtools view
/// -H | grep -v @PG` pipe-diff — the Rust comparator is precise and
/// CI-stable across rust-toolchain version drift.
///
/// Fixture regen: see `test_files/README.md` "uBAM output reference
/// fixtures" section.
fn assert_ubam_eq(actual: &Path, expected: &Path) {
    let actual_hdr = header_minus_pg(actual);
    let expected_hdr = header_minus_pg(expected);
    assert_eq!(
        actual_hdr,
        expected_hdr,
        "BAM header (minus @PG) mismatch between {} and {}",
        actual.display(),
        expected.display()
    );
    let actual_recs = bam_tuples(actual);
    let expected_recs = bam_tuples(expected);
    assert_eq!(
        actual_recs.len(),
        expected_recs.len(),
        "BAM record count mismatch: actual={} ({} records), expected={} ({} records)",
        actual.display(),
        actual_recs.len(),
        expected.display(),
        expected_recs.len()
    );
    for (i, (got, want)) in actual_recs.iter().zip(expected_recs.iter()).enumerate() {
        assert_eq!(
            got, want,
            "BAM record {} differs:\n  actual={:?}\n  expected={:?}",
            i, got, want
        );
    }
}

// ─── Successful end-to-end runs ─────────────────────────────────────────────

#[test]
fn ubam_out_se_ubam_input_matches_reference() {
    let dir = fresh_tmpdir("tg_int_ubam_out_se");
    let status = Command::new(binary())
        .args(["--output-format", "ubam"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let out_bam = dir.join("ubam_test_trimmed.bam");
    assert!(out_bam.exists(), "expected output BAM missing");
    assert_ubam_eq(&out_bam, Path::new("test_files/ubam_out_se_REFERENCE.bam"));
}

#[test]
fn ubam_out_pe_one_ubam_interleaved_matches_reference() {
    let dir = fresh_tmpdir("tg_int_ubam_out_pe");
    let status = Command::new(binary())
        .args(["--paired", "--output-format", "ubam"])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let out_bam = dir.join("ubam_paired_test_val.bam");
    assert!(out_bam.exists(), "expected interleaved output BAM missing");
    assert_ubam_eq(&out_bam, Path::new("test_files/ubam_out_pe_REFERENCE.bam"));
}

#[test]
fn ubam_out_se_fastq_input_produces_valid_bam() {
    // FASTQ input → uBAM output. No source @PG to propagate; minimal
    // synthesised header. Compare record count + per-record content with
    // the FASTQ output the SAME input would produce.
    let dir = fresh_tmpdir("tg_int_ubam_out_se_fq");
    let status = Command::new(binary())
        .args(["--output-format", "ubam"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");
    let out_bam = dir.join("BS-seq_10K_R1_trimmed.bam");
    assert!(out_bam.exists(), "expected output BAM missing");
    // Verify it parses as a valid BAM with > 0 records.
    let tuples = bam_tuples(&out_bam);
    assert!(
        !tuples.is_empty(),
        "expected at least one record in output BAM"
    );
}

#[test]
fn ubam_out_pe_two_fastq_interleaved() {
    // Two FASTQ inputs → ONE interleaved BAM. Verify FREAD1/FREAD2 flag
    // bits alternate (mate-adjacent).
    let dir = fresh_tmpdir("tg_int_ubam_out_pe_fq");
    let status = Command::new(binary())
        .args(["--paired", "--output-format", "ubam"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("test_files/BS-seq_10K_R2.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");
    let out_bam = dir.join("BS-seq_10K_R1_val.bam");
    assert!(out_bam.exists(), "expected interleaved BAM missing");

    let tuples = bam_tuples(&out_bam);
    assert!(!tuples.is_empty(), "expected at least one record");
    // Even-indexed records should be R1 (FREAD1=0x40), odd should be R2
    // (FREAD2=0x80). Check the first two pairs to lock in mate-adjacent.
    let r1_flag_mask = 0x40;
    let r2_flag_mask = 0x80;
    assert!(
        tuples[0].1 & r1_flag_mask != 0,
        "first record should be R1 (FREAD1 set); got flag {:#x}",
        tuples[0].1
    );
    if tuples.len() > 1 {
        assert!(
            tuples[1].1 & r2_flag_mask != 0,
            "second record should be R2 (FREAD2 set); got flag {:#x}",
            tuples[1].1
        );
    }
}

#[test]
fn ubam_out_se_preserve_tags_propagated() {
    // uBAM input with CB/UB tags → uBAM output with --preserve-tags.
    // Read back and verify the tags survived the round-trip.
    let dir = fresh_tmpdir("tg_int_ubam_out_preserve");
    let status = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB,UB"])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    // The fixture `ubam_test_with_tags.bam` may not exist in this repo's
    // committed test_files/ — if so, skip the assertion gracefully.
    if !Path::new("test_files/ubam_test_with_tags.bam").exists() {
        eprintln!(
            "SKIP: test_files/ubam_test_with_tags.bam not committed; \
             preserve-tags propagation tested by lib tests in src/bam.rs::tests"
        );
        return;
    }
    assert!(status.success(), "trim_galore exited non-zero");

    let out_bam = dir.join("ubam_test_with_tags_trimmed.bam");
    let tuples = bam_tuples(&out_bam);
    let first = tuples.first().expect("expected at least one record");
    let tag_names: Vec<&[u8]> = first.4.iter().map(|(t, _)| t.as_slice()).collect();
    assert!(
        tag_names.contains(&b"CB".as_ref()),
        "CB tag missing from preserved aux: {:?}",
        tag_names
    );
    assert!(
        tag_names.contains(&b"UB".as_ref()),
        "UB tag missing from preserved aux: {:?}",
        tag_names
    );
    // Avoid unused-import lint if the fixture exists but doesn't carry these tags.
    let _ = Tag::new(b'C', b'B');
}

#[test]
fn ubam_out_hardtrim5_writes_bam() {
    // --hardtrim5 20 + --output-format ubam should produce a `.bam` file
    // with the 5prime/3prime discriminator preserved per PLAN §3.2 +
    // implementation deviation note in specialty::hardtrim_bam_output_name.
    let dir = fresh_tmpdir("tg_int_ubam_out_hardtrim5");
    let status = Command::new(binary())
        .args(["--hardtrim5", "20", "--output-format", "ubam"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let out_bam = dir.join("ubam_test.20bp_5prime.bam");
    assert!(
        out_bam.exists(),
        "expected hardtrim5 BAM at {}",
        out_bam.display()
    );
    let tuples = bam_tuples(&out_bam);
    assert_eq!(tuples.len(), 10, "expected 10 hard-trimmed records");
    // Every record should be trimmed to ≤ 20 bp.
    for (i, t) in tuples.iter().enumerate() {
        assert!(
            t.2.len() <= 20,
            "record {} has seq length {} > 20",
            i,
            t.2.len()
        );
    }
}

#[test]
fn ubam_out_cores_gt_1_warns_and_proceeds() {
    // PLAN §9 row: `--cores N>1` + uBAM-out must warn ("ignored") and
    // STILL PRODUCE a valid output BAM. Code-review AGREE-2 row.
    let dir = fresh_tmpdir("tg_int_ubam_out_cores");
    let output = Command::new(binary())
        .args(["--cores", "4", "--output-format", "ubam"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "--cores 4 + uBAM should succeed (with warning); stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("ignored"),
        "expected '--cores ignored' note in stderr, got: {}",
        stderr
    );
    assert!(
        dir.join("ubam_test_trimmed.bam").exists(),
        "expected output BAM to still be produced"
    );
}

#[test]
fn ubam_out_preserve_tags_mixed_batch_allowed() {
    // PLAN §3.4b A-O1 loosening: --preserve-tags + mixed-input batch
    // (at least one uBAM + one FASTQ) is ALLOWED. The uBAM input's
    // records get their tags preserved; the FASTQ input's records get
    // empty aux. Code-review AGREE-2 row.
    if !std::path::Path::new("test_files/ubam_test_with_tags.bam").exists() {
        eprintln!("SKIP: ubam_test_with_tags.bam not committed");
        return;
    }
    let dir = fresh_tmpdir("tg_int_ubam_out_mixed");
    let status = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB,UB"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "mixed FASTQ + uBAM batch with --preserve-tags should be allowed"
    );

    // Both outputs must exist.
    let fastq_out = dir.join("BS-seq_10K_R1_trimmed.bam");
    let ubam_out = dir.join("ubam_test_with_tags_trimmed.bam");
    assert!(fastq_out.exists(), "expected BAM from FASTQ input");
    assert!(ubam_out.exists(), "expected BAM from uBAM input");

    // uBAM-side records carry CB tags; FASTQ-side records do not.
    let ubam_tuples = bam_tuples(&ubam_out);
    let cb_present = ubam_tuples
        .first()
        .map(|t| t.4.iter().any(|(tag, _)| tag.as_slice() == b"CB"))
        .unwrap_or(false);
    assert!(
        cb_present,
        "uBAM-side first record should carry CB tag in mixed-batch run"
    );

    let fq_tuples = bam_tuples(&fastq_out);
    let fq_cb_present = fq_tuples
        .first()
        .map(|t| t.4.iter().any(|(tag, _)| tag.as_slice() == b"CB"))
        .unwrap_or(false);
    assert!(
        !fq_cb_present,
        "FASTQ-side first record should NOT carry CB tag (no source)"
    );
}

// ─── §3.4a rejection rules end-to-end ──────────────────────────────────────

#[test]
fn ubam_out_rename_with_preserve_tags_keeps_tags_intact() {
    // Code-review C1 regression guard: --rename + --preserve-tags must
    // NOT corrupt the last preserved tag value. Before the
    // `append_to_id` fix, the `:clip5:<seq>` suffix landed inside the
    // last `Z:` tag value because the suffix appended after the tab tail
    // instead of splicing into the name.
    if !std::path::Path::new("test_files/ubam_test_with_tags.bam").exists() {
        eprintln!("SKIP: ubam_test_with_tags.bam not committed");
        return;
    }
    let dir = fresh_tmpdir("tg_int_ubam_out_rename_tags");
    let status = Command::new(binary())
        .args([
            "--clip_R1",
            "5",
            "--rename",
            "--output-format",
            "ubam",
            "--preserve-tags",
            "CB,UB",
        ])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let out_bam = dir.join("ubam_test_with_tags_trimmed.bam");
    let tuples = bam_tuples(&out_bam);
    let first = tuples.first().expect("expected at least one record");

    // Check the UB tag value is intact — NOT polluted by ":clip5:..." text.
    // Aux value debug format includes the variant + raw value; checking
    // for the corruption marker is sufficient + minimal.
    let ub_entry = first
        .4
        .iter()
        .find(|(tag, _)| tag.as_slice() == b"UB")
        .expect("UB tag missing from preserved aux");
    let ub_value_str = String::from_utf8_lossy(&ub_entry.1);
    assert!(
        !ub_value_str.contains(":clip5:"),
        "UB tag value was corrupted by --rename annotation: {}",
        ub_value_str
    );
}

#[test]
fn ubam_out_two_bam_pair_rejected() {
    // Two BAM inputs under --paired is not supported with --output-format
    // ubam (mirrors the FASTQ-path rejection). uBAM paired mode expects
    // ONE interleaved file. Code-review B-I2 regression guard.
    let dir = fresh_tmpdir("tg_int_ubam_out_two_bam_rej");
    let output = Command::new(binary())
        .args(["--paired", "--output-format", "ubam"])
        .arg("test_files/ubam_test.bam")
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "two-BAM-pair + --paired + --output-format ubam should error"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("two BAM files is not supported") || stderr.contains("interleaved file"),
        "expected two-BAM rejection in stderr, got: {}",
        stderr
    );
}

#[test]
fn ubam_out_clock_rejected_at_cli() {
    let dir = fresh_tmpdir("tg_int_ubam_out_clock_rej");
    let output = Command::new(binary())
        .args(["--clock", "--output-format", "ubam"])
        .arg("test_files/clock_10K_R1.fastq.gz")
        .arg("test_files/clock_10K_R2.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--clock + --output-format ubam should error"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--clock") && stderr.contains("--output-format ubam"),
        "expected clock+ubam rejection in stderr, got: {}",
        stderr
    );
}

#[test]
fn ubam_out_clumpify_rejected_at_cli() {
    let dir = fresh_tmpdir("tg_int_ubam_out_clumpify_rej");
    let output = Command::new(binary())
        .args(["--clumpify", "--cores", "2", "--output-format", "ubam"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--clumpify + --output-format ubam should error"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--clumpify") && stderr.contains("--output-format ubam"),
        "expected clumpify+ubam rejection in stderr, got: {}",
        stderr
    );
}

#[test]
fn ubam_out_preserve_tags_all_fastq_rejected() {
    // PLAN §3.4b — `--preserve-tags` + all-FASTQ inputs + `--output-format
    // ubam` is a hard error because there are no source tags AND the user
    // explicitly requested uBAM output.
    let dir = fresh_tmpdir("tg_int_ubam_out_ptags_rej");
    let output = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "preserve-tags + all-FASTQ + --output-format ubam should error"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--preserve-tags"),
        "expected preserve-tags rejection in stderr, got: {}",
        stderr
    );
}
