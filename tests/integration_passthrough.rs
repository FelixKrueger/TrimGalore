//! Binary-driven integration tests for `--passthrough`.
//!
//! These tests close the §7 coverage gap (`--cores 1` dispatcher integration)
//! from the plan v2 validation matrix — exercising
//! `main.rs::run_paired` end-to-end via the built binary rather than the
//! library entry points. That covers wiring the lib tests can't reach:
//!   * output-path computation (`io::passthrough_output_name`)
//!   * the output-collision pre-flight extension (`main.rs:177–214`)
//!   * sanity-check ordering for the passthrough file
//!   * FastQC sweep extension
//!   * report block emission (text + JSON)
//!
//! `CARGO_BIN_EXE_trim_galore` is auto-set by cargo for integration tests;
//! it resolves to `target/{debug,release}/trim_galore` depending on the
//! cargo invocation. `cargo test` builds the binary before running these.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Locate the `trim_galore` binary built by cargo for this test target.
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// Create a fresh temp dir for the test, removing any leftover from a prior run.
fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Write a tiny three-file Multiome fixture (R1, R2, I1) to `dir`. Returns
/// the three paths. Uses plain `.fq` (no gzip) for test speed — the binary's
/// gzip mode is decided per the R1 input extension, so plain in → plain out.
fn write_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    let i1 = dir.join("i1.fq");

    // 30 records. Every 3rd has a short R1 (under length cutoff 20) so it
    // gets dropped — exercises both Pass and Discard branches and
    // demonstrates lockstep is preserved across drops. R2 also short on
    // dropped rows (avoids the pre-existing parallel/serial r2_unpaired
    // accounting divergence noted in the lib parity test).
    let mut r1_s = String::new();
    let mut r2_s = String::new();
    let mut i1_s = String::new();
    for i in 0..30 {
        let short = i % 3 == 0;
        let r1_seq = if short {
            "ACGTACGTACGTACGT" // 16 bp — dropped
        } else {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC" // 42 bp
        };
        let r2_seq = if short {
            "TGCATGCATGCATGCA" // 16 bp
        } else {
            "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
        };
        let pt_seq = "AAAACCCCGGGGTTTT"; // 16 bp cell barcode

        r1_s.push_str(&format!(
            "@read_{i}\n{r1_seq}\n+\n{}\n",
            "I".repeat(r1_seq.len())
        ));
        r2_s.push_str(&format!(
            "@read_{i}\n{r2_seq}\n+\n{}\n",
            "I".repeat(r2_seq.len())
        ));
        i1_s.push_str(&format!(
            "@read_{i}\n{pt_seq}\n+\n{}\n",
            "I".repeat(pt_seq.len())
        ));
    }
    std::fs::write(&r1, r1_s).unwrap();
    std::fs::write(&r2, r2_s).unwrap();
    std::fs::write(&i1, i1_s).unwrap();
    (r1, r2, i1)
}

/// Count records in a plain FASTQ file (line count / 4).
fn count_records(path: &Path) -> usize {
    let s = std::fs::read_to_string(path).expect("read FASTQ");
    s.lines().filter(|l| l.starts_with('@')).count()
}

/// §7 — end-to-end smoke via the built binary at `--cores 1`. Exercises:
///   * Cli parsing of `--passthrough`
///   * Cli::validate's 9-item compatibility envelope
///   * main.rs::run_paired dispatch (which routes to `trimmer::run_paired_end`
///     in the serial path at `cli.cores == 1`)
///   * Three output files actually land at the expected paths
///   * Each output has the same record count (= 20 survivors out of 30)
///   * Trimming reports get written (text + JSON)
#[test]
fn passthrough_cli_cores_1_smoke() {
    let dir = fresh_tmpdir("tg_pt_cli_cores1");
    let (r1, r2, i1) = write_fixture(&dir);

    let out = Command::new(binary())
        .arg("--paired")
        .arg("--length")
        .arg("20")
        .arg("--passthrough")
        .arg(&i1)
        .arg("--output_dir")
        .arg(&dir)
        .arg("--dont_gzip")
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("spawn trim_galore");

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.success(),
        "trim_galore exited with status {:?}\nstderr:\n{stderr}",
        out.status.code()
    );

    let out_r1 = dir.join("r1_val_1.fq");
    let out_r2 = dir.join("r2_val_2.fq");
    let out_pt = dir.join("i1_passthrough.fq");
    assert!(
        out_r1.exists(),
        "missing {}\nstderr:\n{stderr}",
        out_r1.display()
    );
    assert!(
        out_r2.exists(),
        "missing {}\nstderr:\n{stderr}",
        out_r2.display()
    );
    assert!(
        out_pt.exists(),
        "missing {}\nstderr:\n{stderr}",
        out_pt.display()
    );

    // Record counts: 30 input, 10 short (every 3rd: i = 0, 3, ..., 27 → 10),
    // 20 survivors. All three outputs must agree.
    let c1 = count_records(&out_r1);
    let c2 = count_records(&out_r2);
    let cp = count_records(&out_pt);
    assert_eq!(c1, 20, "R1 out: expected 20 records, got {c1}");
    assert_eq!(c2, 20, "R2 out: expected 20 records, got {c2}");
    assert_eq!(cp, 20, "passthrough out: expected 20 records, got {cp}");

    // Trimming reports written (one per pair-mate, in --output_dir).
    let r1_report = dir.join("r1.fq_trimming_report.txt");
    let r1_json = dir.join("r1.fq_trimming_report.json");
    assert!(
        r1_report.exists(),
        "missing text report at {}\nstderr:\n{stderr}",
        r1_report.display()
    );
    assert!(
        r1_json.exists(),
        "missing JSON report at {}\nstderr:\n{stderr}",
        r1_json.display()
    );

    // The R2 text report (where pair-validation + passthrough blocks land)
    // should contain the `=== Passthrough file ===` block.
    let r2_report_text =
        std::fs::read_to_string(dir.join("r2.fq_trimming_report.txt")).unwrap_or_default();
    assert!(
        r2_report_text.contains("=== Passthrough file ==="),
        "R2 text report missing passthrough block:\n{r2_report_text}"
    );
    assert!(
        r2_report_text.contains("Records carried through:"),
        "R2 text report missing 'Records carried through' line"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// Negative-path: `--passthrough` without `--paired` is rejected by
/// validation at startup, before any output is touched.
#[test]
fn passthrough_without_paired_rejected() {
    let dir = fresh_tmpdir("tg_pt_cli_no_paired");
    let (r1, _r2, i1) = write_fixture(&dir);

    let out = Command::new(binary())
        .arg("--passthrough")
        .arg(&i1)
        .arg(&r1)
        .output()
        .expect("spawn trim_galore");

    assert!(
        !out.status.success(),
        "trim_galore should have exited with error"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--passthrough requires --paired"),
        "unexpected stderr: {stderr}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
