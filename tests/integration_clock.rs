//! `--clock` summary line: the count is of pairs WITHOUT the fixed sequence
//! CAGT in both reads, so the label must say so (Perl 0.6.11 said the opposite).

use std::process::Command;

#[test]
fn clock_summary_labels_the_count_of_pairs_without_cagt() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(["--clock", "-o"])
        .arg(dir.path())
        .args([
            "test_files/clock_10K_R1.fastq.gz",
            "test_files/clock_10K_R2.fastq.gz",
        ])
        .output()
        .expect("failed to run trim_galore");
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(out.status.success(), "stderr: {stderr}");
    assert!(
        stderr.contains("Sequences processed in total: 10000\n"),
        "stderr: {stderr}"
    );
    // 1073 of the 10000 pairs lack CAGT at bases 9-12 of R1 or R2.
    assert!(
        stderr.contains(
            "thereof did not have the fixed sequence CAGT in both R1 and R2: 1073 (10.73%)\n"
        ),
        "stderr: {stderr}"
    );
}
