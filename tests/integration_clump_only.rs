//! Binary-driven integration tests for `--clump_only`.
//!
//! Covers the load-bearing byte-identity invariant end-to-end by spawning
//! the built binary rather than exercising library code directly. Uses the
//! same "reconstitute 4-line FASTQ records then multiset-diff" pattern the
//! CI's byte-identity step will use (see `.github/workflows/ci.yml`), so
//! the on-disk assertion here is the same shape as the assertion in CI.
//!
//! Fixtures: reuses `test_files/BS-seq_10K_R{1,2}.fastq.gz` (real BS-seq
//! data, bare `+` on line 3, LF endings).

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_clump_only_it_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Slurp a plain-or-gzip FASTQ file and reconstitute the multiset of
/// 4-line records: each record → one string (tab-joined lines), then
/// collect into a sorted vec for permutation comparison.
///
/// Sorting *whole records* (not individual lines) is what makes this
/// resilient to the class of bug where lines from different records get
/// mixed up (e.g. R1's quality attached to R2's sequence) — a naïve
/// `zcat | sort | diff` would miss that but `paste - - - -` before sort
/// catches it. See PLAN.md §Implementation outline step 9.
fn record_multiset(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let boxed: Box<dyn Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(boxed);
    let lines: Vec<String> = reader.lines().collect::<std::io::Result<Vec<_>>>()?;
    assert_eq!(
        lines.len() % 4,
        0,
        "FASTQ line count not divisible by 4: {}",
        path.display()
    );
    let mut records: Vec<String> = lines
        .chunks(4)
        .map(|c| format!("{}\t{}\t{}\t{}", c[0], c[1], c[2], c[3]))
        .collect();
    records.sort();
    Ok(records)
}

fn fixture(name: &str) -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest_dir).join("test_files").join(name)
}

// ── Byte-identity: SE (gzip in / gzip out) ─────────────────────────

#[test]
fn se_byte_identity_gzip() -> Result<()> {
    let dir = tempdir("se_gz");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    assert!(input.exists(), "fixture missing: {}", input.display());

    let status = Command::new(binary())
        .args(["--clump_only", "--cores", "2", "-o", dir.to_str().unwrap()])
        .arg(&input)
        .status()?;
    assert!(status.success(), "trim_galore --clump_only failed");

    // Output filename: <stem>_clumped.fq.gz — matches io::clumped_output_name.
    let out = dir.join("BS-seq_10K_R1_clumped.fq.gz");
    assert!(out.exists(), "expected output missing: {}", out.display());

    let in_records = record_multiset(&input)?;
    let out_records = record_multiset(&out)?;
    assert_eq!(
        in_records.len(),
        out_records.len(),
        "record count preserved (SE gzip)"
    );
    assert_eq!(
        in_records, out_records,
        "record multiset preserved byte-identically (SE gzip)"
    );
    Ok(())
}

// ── Byte-identity: SE (gzip in / plain out via --dont_gzip) ────────

#[test]
fn se_byte_identity_dont_gzip() -> Result<()> {
    let dir = tempdir("se_plain");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    assert!(input.exists());

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--dont_gzip",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()?;
    assert!(
        status.success(),
        "trim_galore --clump_only --dont_gzip failed"
    );

    let out = dir.join("BS-seq_10K_R1_clumped.fq");
    assert!(
        out.exists(),
        "expected plain output missing: {}",
        out.display()
    );

    let in_records = record_multiset(&input)?;
    let out_records = record_multiset(&out)?;
    assert_eq!(
        in_records, out_records,
        "record multiset preserved (--dont_gzip)"
    );
    Ok(())
}

// ── Byte-identity: PE ──────────────────────────────────────────────

#[test]
fn pe_byte_identity() -> Result<()> {
    let dir = tempdir("pe");
    let r1 = fixture("BS-seq_10K_R1.fastq.gz");
    let r2 = fixture("BS-seq_10K_R2.fastq.gz");
    assert!(r1.exists() && r2.exists());

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&r1)
        .arg(&r2)
        .status()?;
    assert!(status.success(), "trim_galore --clump_only --paired failed");

    let out_r1 = dir.join("BS-seq_10K_R1_clumped_1.fq.gz");
    let out_r2 = dir.join("BS-seq_10K_R2_clumped_2.fq.gz");
    assert!(out_r1.exists(), "R1 output missing: {}", out_r1.display());
    assert!(out_r2.exists(), "R2 output missing: {}", out_r2.display());

    // Per-mate multiset preserved.
    assert_eq!(record_multiset(&r1)?, record_multiset(&out_r1)?);
    assert_eq!(record_multiset(&r2)?, record_multiset(&out_r2)?);
    Ok(())
}

// ── Cross-run determinism ──────────────────────────────────────────

#[test]
fn se_deterministic_across_runs() -> Result<()> {
    let dir1 = tempdir("det1");
    let dir2 = tempdir("det2");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    assert!(input.exists());

    let run = |out_dir: &Path| -> Result<()> {
        let status = Command::new(binary())
            .args([
                "--clump_only",
                "--cores",
                "2",
                "-o",
                out_dir.to_str().unwrap(),
            ])
            .arg(&input)
            .status()?;
        assert!(status.success());
        Ok(())
    };
    run(&dir1)?;
    run(&dir2)?;

    let out1 = dir1.join("BS-seq_10K_R1_clumped.fq.gz");
    let out2 = dir2.join("BS-seq_10K_R1_clumped.fq.gz");
    let mut b1 = Vec::new();
    let mut b2 = Vec::new();
    File::open(&out1)?.read_to_end(&mut b1)?;
    File::open(&out2)?.read_to_end(&mut b2)?;
    assert_eq!(b1, b2, "cross-run output byte-identity failed");
    Ok(())
}

// ── Report filename discipline ─────────────────────────────────────

#[test]
fn report_uses_clumping_not_trimming_filename() -> Result<()> {
    // Guard the nf-core scan-glob protection: the mode emits
    // `*_clumping_report.txt` and NOT `*_trimming_report.txt`.
    let dir = tempdir("rep");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let status = Command::new(binary())
        .args(["--clump_only", "--cores", "2", "-o", dir.to_str().unwrap()])
        .arg(&input)
        .status()?;
    assert!(status.success());

    // clumping_report exists.
    let clumping_report = dir.join("BS-seq_10K_R1.fastq.gz_clumping_report.txt");
    assert!(
        clumping_report.exists(),
        "clumping report missing: {}",
        clumping_report.display()
    );

    // trimming_report does NOT exist (glob-scan protection for nf-core).
    let trimming_report = dir.join("BS-seq_10K_R1.fastq.gz_trimming_report.txt");
    assert!(
        !trimming_report.exists(),
        "trimming report must NOT be created under --clump_only: {}",
        trimming_report.display()
    );

    // Sentinel "Mode: --clump_only" appears in the clumping report.
    let contents = std::fs::read_to_string(&clumping_report)?;
    assert!(
        contents.contains("Mode: --clump_only"),
        "clumping report missing sentinel line; contents:\n{contents}"
    );
    Ok(())
}

// ── Rejection matrix: representative flags ─────────────────────────

#[test]
fn rejects_length_flag() -> Result<()> {
    let dir = tempdir("rej_len");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--length",
            "20",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()?;
    assert!(
        !out.status.success(),
        "--clump_only --length must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--length") || stderr.contains("does not filter"),
        "stderr should explain the --length rejection; got: {stderr}"
    );
    Ok(())
}

#[test]
fn rejects_adapter_flag() -> Result<()> {
    let dir = tempdir("rej_adapt");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "-a",
            "AGATCGGAAGAGC",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()?;
    assert!(!out.status.success(), "--clump_only -a must be rejected");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--adapter") || stderr.contains("does not trim"),
        "stderr should explain the adapter rejection; got: {stderr}"
    );
    Ok(())
}

#[test]
fn rejects_ubam_input_without_ubam_output() -> Result<()> {
    // v2: uBAM input IS supported under `--clump_only`, but only when
    // paired with `--output-format ubam`. On the plain FASTQ-output path
    // (the default), uBAM input would drop aux tags — so it's rejected
    // with a message pointing the user at the uBAM-output path. See
    // tests/integration_clump_only_ubam.rs for the positive-path tests.
    let dir = tempdir("rej_ubam_no_output_fmt");
    let input = fixture("ubam_test.bam");
    assert!(input.exists(), "uBAM fixture missing: {}", input.display());
    let out = Command::new(binary())
        .args(["--clump_only", "--cores", "2", "-o", dir.to_str().unwrap()])
        .arg(&input)
        .output()?;
    assert!(
        !out.status.success(),
        "--clump_only on uBAM input (without --output-format ubam) must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--output-format ubam"),
        "stderr should route the user to --output-format ubam; got: {stderr}"
    );
    Ok(())
}

#[test]
fn rejects_rename_flag() -> Result<()> {
    let dir = tempdir("rej_rename");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--rename",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()?;
    assert!(
        !out.status.success(),
        "--clump_only --rename must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--rename") || stderr.contains("byte-identically"),
        "stderr should explain the --rename rejection; got: {stderr}"
    );
    Ok(())
}

/// Helper: run --clump_only with an extra flag/value pair and assert
/// (a) exit 0 and (b) byte-identity preserved. Used by the silent-accept
/// regression guards for `-q`, `--stringency`, `-e`.
fn assert_silently_accepts(tag: &str, extra_args: &[&str]) -> Result<()> {
    let dir = tempdir(tag);
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let mut cmd = Command::new(binary());
    cmd.args(["--clump_only", "--cores", "2", "-o", dir.to_str().unwrap()]);
    for a in extra_args {
        cmd.arg(a);
    }
    cmd.arg(&input);
    let status = cmd.status()?;
    assert!(
        status.success(),
        "--clump_only {extra_args:?} should succeed (silent-accept per Blocker 4)"
    );
    let out = dir.join("BS-seq_10K_R1_clumped.fq.gz");
    assert!(out.exists());
    assert_eq!(record_multiset(&input)?, record_multiset(&out)?);
    Ok(())
}

#[test]
fn silently_accepts_quality_flag() -> Result<()> {
    // -q has a clap default so `Cli::validate()` can't distinguish
    // user-set from default. Per Blocker 4 (ii), we silent-accept and
    // document ignore-behavior. Regression guard: `--clump_only -q 30`
    // must still succeed (not reject) and produce byte-identical output.
    assert_silently_accepts("silent_q", &["-q", "30"])
}

#[test]
fn silently_accepts_stringency_flag() -> Result<()> {
    // --stringency has a clap default (adapter-match min-overlap); same
    // silent-accept semantics as -q under --clump_only.
    assert_silently_accepts("silent_stringency", &["--stringency", "5"])
}

#[test]
fn silently_accepts_error_rate_flag() -> Result<()> {
    // -e / --error has a clap default (adapter-match max error rate);
    // same silent-accept semantics as -q under --clump_only.
    assert_silently_accepts("silent_e", &["-e", "0.05"])
}

/// #421 — `--fastqc_args` alone must activate FastQC on the `--clump_only` drivers.
#[test]
fn clump_only_fastqc_args_alone_produces_a_report() -> Result<()> {
    let dir = tempdir("fastqc_args_only");
    let status = Command::new(binary())
        .args(["--clump_only", "--fastqc_args", "--quiet"])
        .arg(fixture("BS-seq_10K_R1.fastq.gz"))
        .arg("-o")
        .arg(&dir)
        .status()?;
    assert!(status.success(), "--clump_only --fastqc_args must succeed");

    let reports: Vec<String> = std::fs::read_dir(&dir)?
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .filter(|n| n.contains("_fastqc."))
        .collect();
    assert_eq!(
        reports.len(),
        2,
        "expected _fastqc.html + _fastqc.zip, got {reports:?}"
    );
    Ok(())
}
