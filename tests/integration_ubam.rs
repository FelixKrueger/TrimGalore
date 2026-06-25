//! Binary-driven integration tests for uBAM input support (#316).
//!
//! These tests exercise `main.rs::run_single_file` / `run_paired_ubam_single_file`
//! end-to-end via the built `trim_galore` binary, covering wiring the lib
//! tests can't reach:
//!   * format detection → reader factory dispatch
//!   * output-path computation for `.bam` inputs (single-file paired-uBAM)
//!   * `--clumpify` + BAM (worker-pool dispatch through `RecordSource` trait)
//!   * specialty mode (`--hardtrim5`) + BAM (alternate entry points)
//!   * aligned-BAM rejection (per-record check in `BamReader::next_record`)
//!
//! Output-content parity is checked via `(id, seq, qual)` tuples — resilient
//! to gzip framing drift between runs. See PLAN §6.3 tier-1 assertion.

use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Parse a plain FASTQ file into `(id, seq, qual)` tuples. Used for
/// content-tuple parity (gzip-framing-resilient).
fn read_fastq_tuples(path: &Path) -> Vec<(String, String, String)> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path.display(), e));
    let mut tuples = Vec::new();
    let mut iter = text.lines();
    while let Some(id) = iter.next() {
        let seq = iter.next().unwrap_or_default();
        let _plus = iter.next().unwrap_or_default();
        let qual = iter.next().unwrap_or_default();
        tuples.push((id.to_string(), seq.to_string(), qual.to_string()));
    }
    tuples
}

#[test]
fn single_end_ubam_matches_reference() {
    let dir = fresh_tmpdir("tg_int_se_ubam");
    let status = Command::new(binary())
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got = read_fastq_tuples(&dir.join("ubam_test_trimmed.fq"));
    let want = read_fastq_tuples(Path::new("test_files/ubam_test_trimmed_REFERENCE.fq"));
    assert_eq!(
        got, want,
        "SE uBAM output did not match committed reference (record-by-record content-tuple comparison)"
    );
    assert_eq!(got.len(), 10, "expected 10 trimmed records");
}

#[test]
fn paired_end_interleaved_ubam_matches_reference() {
    let dir = fresh_tmpdir("tg_int_pe_ubam");
    let status = Command::new(binary())
        .args(["--paired"])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got_r1 = read_fastq_tuples(&dir.join("ubam_paired_test_val_1.fq"));
    let want_r1 = read_fastq_tuples(Path::new("test_files/ubam_paired_test_val_1_REFERENCE.fq"));
    let got_r2 = read_fastq_tuples(&dir.join("ubam_paired_test_val_2.fq"));
    let want_r2 = read_fastq_tuples(Path::new("test_files/ubam_paired_test_val_2_REFERENCE.fq"));

    assert_eq!(got_r1, want_r1, "PE R1 output diverged from reference");
    assert_eq!(got_r2, want_r2, "PE R2 output diverged from reference");
    assert_eq!(got_r1.len(), 10);
    assert_eq!(got_r2.len(), 10);

    // Same-template invariant: matched-pair IDs must align between R1/R2.
    for (a, b) in got_r1.iter().zip(got_r2.iter()) {
        assert_eq!(a.0, b.0, "R1/R2 records out of sync at id {}", a.0);
    }
}

#[test]
fn clumpify_plus_ubam_runs_clean() {
    // PLAN §3.4: --clumpify + BAM is allowed (works transparently via the
    // shared FastqRecord shape inside the worker pool). Smoke test: just
    // verify it exits non-zero-free and produces output.
    let dir = fresh_tmpdir("tg_int_clumpify_ubam");
    let status = Command::new(binary())
        .args(["--clumpify", "--cores", "2"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "--clumpify + uBAM smoke test failed; PLAN §3.4 says it should work"
    );
    assert!(
        dir.join("ubam_test_trimmed.fq.gz").exists() || dir.join("ubam_test_trimmed.fq").exists(),
        "expected trimmed output not found under {}",
        dir.display()
    );
}

#[test]
fn hardtrim5_plus_ubam_runs_clean() {
    // PLAN §3.4: specialty modes (--hardtrim5/3 / --clock / --implicon)
    // work with BAM input because they loop over `next_record()` via the
    // same RecordSource path. Smoke test the simplest specialty mode.
    let dir = fresh_tmpdir("tg_int_hardtrim5_ubam");
    let status = Command::new(binary())
        .args(["--hardtrim5", "20"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "--hardtrim5 + uBAM smoke test failed; PLAN §3.4 says it should work"
    );
}

#[test]
fn passthrough_plus_ubam_rejected() {
    // PLAN §3.4: --passthrough + BAM is rejected in v1 (three-way ID sync
    // would need re-thinking for BAM record offsets). Confirm the error
    // surfaces with the expected message text.
    let dir = fresh_tmpdir("tg_int_passthrough_ubam_reject");
    let output = Command::new(binary())
        .args([
            "--paired",
            "--passthrough",
            "test_files/BS-seq_10K_I1.fastq.gz",
        ])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--passthrough + uBAM must be rejected"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    // Cli::validate's existing passthrough rule ("--passthrough requires
    // exactly one R1/R2 pair") fires FIRST when --paired+uBAM with N=1 is
    // combined with --passthrough, because validate() runs before
    // sanity_check_any. Either rejection path is acceptable — the contract
    // is "--passthrough + uBAM in v1 is rejected, with a clear pointer at
    // either fix". Accept either rejection message.
    assert!(
        stderr.contains("--passthrough")
            && (stderr.contains("not supported") || stderr.contains("requires exactly")),
        "expected some passthrough rejection message, got stderr: {}",
        stderr
    );
}

#[test]
fn preserve_tags_roundtrip_matches_golden() {
    // PLAN §3.2.5 + T23. The committed `ubam_test_with_tags.bam` carries
    // CB:Z:ATCGATCG-1 and UB:Z:GCTAGCTA aux tags on every record. With
    // `--preserve-tags CB,UB`, the FASTQ headers should be
    // `@<name>\tCB:Z:ATCGATCG-1\tUB:Z:GCTAGCTA` (samtools -T-compatible).
    // Verify against a committed golden reference.
    let dir = fresh_tmpdir("tg_int_preserve_tags");
    let status = Command::new(binary())
        .args(["--preserve-tags", "CB,UB"])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got = read_fastq_tuples(&dir.join("ubam_test_with_tags_trimmed.fq"));
    let want = read_fastq_tuples(Path::new(
        "test_files/ubam_test_with_tags_trimmed_REFERENCE.fq",
    ));
    assert_eq!(
        got, want,
        "--preserve-tags output diverged from committed golden reference"
    );
    // Spot-check the tag-format invariant on the first record.
    assert!(
        got[0].0.contains("\tCB:Z:ATCGATCG-1"),
        "first record id must carry CB tag in user-specified order: {:?}",
        got[0].0
    );
    assert!(
        got[0].0.contains("\tUB:Z:GCTAGCTA"),
        "first record id must carry UB tag: {:?}",
        got[0].0
    );
}

#[test]
fn preserve_tags_user_specified_order_honoured() {
    // Same fixture, but user requests UB then CB — the resulting header must
    // be `...\tUB:...\tCB:...`. Tag order is policy-defined by the flag, NOT
    // by BAM file aux-field order.
    let dir = fresh_tmpdir("tg_int_preserve_tags_order");
    let status = Command::new(binary())
        .args(["--preserve-tags", "UB,CB"])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success());
    let got = read_fastq_tuples(&dir.join("ubam_test_with_tags_trimmed.fq"));
    // UB must come BEFORE CB in the header now.
    let id = &got[0].0;
    let ub_pos = id.find("UB:Z:").expect("UB tag missing");
    let cb_pos = id.find("CB:Z:").expect("CB tag missing");
    assert!(
        ub_pos < cb_pos,
        "tag order must match --preserve-tags argument (UB before CB), got id: {:?}",
        id
    );
}

#[test]
fn paired_with_single_fastq_rejected() {
    // PLAN §3.3: --paired with a single input file is only legal if that
    // file is a uBAM. A single FASTQ in --paired mode must be rejected.
    let dir = fresh_tmpdir("tg_int_paired_single_fastq_reject");
    let output = Command::new(binary())
        .args(["--paired"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--paired with single FASTQ must be rejected"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("only legal if that file is a uBAM"),
        "expected single-FASTQ rejection message, got stderr: {}",
        stderr
    );
}
