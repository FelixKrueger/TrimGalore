//! Binary-driven tests for the `--clumpify` calibrated-floor warning.
//!
//! The sizing coefficients are bounded at 25 bp, so shorter reads are an
//! extrapolation and get a warning. Two things need holding: the warning fires
//! below the floor, and it does *not* fire at it — the CI memory cells run a
//! 25 bp fixture whose early records carry short headers, and an earlier
//! record-size-based form of this check warned on every one of them.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_floor_warn_{tag}_{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// One record with a deliberately terse header, so the record's *byte* size is
/// far below the floor even when its read length is not.
fn write_fastq(dir: &Path, name: &str, read_len: usize) -> PathBuf {
    let p = dir.join(name);
    let seq = "ACGT".repeat(read_len.div_ceil(4))[..read_len].to_string();
    let qual = "I".repeat(read_len);
    fs::write(&p, format!("@SRR1234567.1\n{seq}\n+\n{qual}\n")).unwrap();
    p
}

fn stderr_of(args: &[&str]) -> String {
    let out = Command::new(binary())
        .args(args)
        .output()
        .expect("binary must run");
    String::from_utf8_lossy(&out.stderr).into_owned()
}

const NEEDLE: &str = "below the 25 bp the --memory reservation was calibrated at";

#[test]
fn warns_below_the_calibrated_floor() {
    let d = tempdir("below");
    let fq = write_fastq(&d, "short.fq", 19);
    let err = stderr_of(&[
        "--clumpify",
        "--cores",
        "2",
        "--memory",
        "1G",
        "-o",
        d.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        err.contains(NEEDLE) && err.contains("shortest input read is 19 bp"),
        "expected the floor warning, got:\n{err}"
    );
}

#[test]
fn silent_at_the_calibrated_floor() {
    let d = tempdir("at");
    // 25 bp with a 13-byte header is 68 record bytes, well under the 97 that an
    // earlier form of this check compared against.
    let fq = write_fastq(&d, "at_floor.fq", 25);
    let err = stderr_of(&[
        "--clumpify",
        "--cores",
        "2",
        "--memory",
        "1G",
        "-o",
        d.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        !err.contains(NEEDLE),
        "warned at the floor, not below it:\n{err}"
    );
}

#[test]
fn checks_every_input_not_just_the_first() {
    let d = tempdir("multi");
    let long = write_fastq(&d, "p1_R1.fq", 51);
    let long2 = write_fastq(&d, "p1_R2.fq", 51);
    let short = write_fastq(&d, "p2_R1.fq", 15);
    let short2 = write_fastq(&d, "p2_R2.fq", 15);
    let err = stderr_of(&[
        "--clumpify",
        "--paired",
        "--cores",
        "2",
        "--memory",
        "1G",
        "-o",
        d.to_str().unwrap(),
        long.to_str().unwrap(),
        long2.to_str().unwrap(),
        short.to_str().unwrap(),
        short2.to_str().unwrap(),
    ]);
    assert!(
        err.contains("shortest input read is 15 bp"),
        "a short second pair was hidden behind a long first pair:\n{err}"
    );
}

#[test]
fn clump_only_gets_the_warning_too() {
    let d = tempdir("only");
    let fq = write_fastq(&d, "short.fq", 19);
    let err = stderr_of(&[
        "--clump_only",
        "--cores",
        "2",
        "--memory",
        "1G",
        "-o",
        d.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        err.contains(NEEDLE),
        "--clump_only never reaches resolve_clump_layout, so it needs its own call:\n{err}"
    );
}
