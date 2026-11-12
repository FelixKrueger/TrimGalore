//! A `SEQ = *` uBAM record decodes to a zero-length read, which then has to
//! survive the trimming pipeline rather than only the reader.

use std::process::Command;

use trim_galore::bam::BamWriter;
use trim_galore::fastq::FastqRecord;

/// One empty record between two ordinary ones, so a pipeline that chokes on the
/// empty read loses the record after it too.
fn write_bam(path: &std::path::Path) {
    let mut w = BamWriter::create(path, None, &[], "trim_galore empty-seq test", 33)
        .expect("failed to open BAM writer");
    for (id, seq) in [
        ("@first", "ACGTACGTACGTACGTACGTACGTACGTACGT"),
        ("@empty", ""),
        ("@last", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
    ] {
        w.write_record(
            &FastqRecord {
                id: id.to_string(),
                seq: seq.to_string(),
                qual: "I".repeat(seq.len()),
            },
            None,
        )
        .expect("failed to write record");
    }
    w.finish().expect("failed to finish BAM");
}

#[test]
fn empty_sequence_record_survives_the_trimming_pipeline() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("empty_seq.bam");
    write_bam(&bam);

    let out = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(["--no_poly_g", "-o"])
        .arg(dir.path())
        .arg(&bam)
        .output()
        .expect("failed to run trim_galore");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "stderr: {stderr}");

    // uBAM input is not itself gzip, so the FASTQ output is written plain.
    let trimmed = dir.path().join("empty_seq_trimmed.fq");
    assert!(
        trimmed.exists(),
        "no output at {}; stderr: {stderr}",
        trimmed.display()
    );
    let text = std::fs::read_to_string(&trimmed).unwrap();

    // The zero-length read is below the default length cutoff, so it is
    // filtered; the records on either side of it must both be written.
    assert!(text.contains("@first\n"), "first read missing:\n{text}");
    assert!(text.contains("@last\n"), "last read missing:\n{text}");
    assert!(!text.contains("@empty\n"), "empty read written:\n{text}");
    assert_eq!(
        text.lines().count(),
        8,
        "expected exactly two records:\n{text}"
    );
}
