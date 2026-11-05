//! `--phred64` must reach the bundled FastQC config.
//!
//! FastQC 0.13.0 reads Phred+33 unless the encoding is declared; 0.12.1 inferred
//! it from the lowest quality character and happened to get Phred+64 right. Trim
//! Galore writes FASTQ quality verbatim, so `--phred64` output is Phred+64 on disk
//! and the engine has to be told, or every mean comes out 31 too high. BAM declares
//! Phred+33 itself and the format's own encoding wins upstream, so BAM paths pass
//! `false`.
//!
//! Every other FastQC test asserts only that `*_fastqc.zip` exists. These read
//! inside it, because a wrong number in the report is invisible from outside.

use std::io::Read;
use std::path::Path;
use std::process::Command;

use tempfile::tempdir;

/// One `#Measure`/value row from the report's Basic Statistics section.
fn basic_stat(zip_path: &Path, stem: &str, measure: &str) -> String {
    let file = std::fs::File::open(zip_path)
        .unwrap_or_else(|e| panic!("opening {}: {e}", zip_path.display()));
    let mut archive = zip::ZipArchive::new(file).expect("valid zip");
    let mut data = String::new();
    archive
        .by_name(&format!("{stem}_fastqc/fastqc_data.txt"))
        .expect("fastqc_data.txt in archive")
        .read_to_string(&mut data)
        .expect("readable");
    data.lines()
        .find(|l| l.starts_with(measure))
        .unwrap_or_else(|| panic!("no '{measure}' row in:\n{data}"))
        .split('\t')
        .nth(1)
        .expect("tab-separated value")
        .to_string()
}

/// Mean quality at `base` from the Per base sequence quality section.
fn per_base_mean(zip_path: &Path, stem: &str, base: &str) -> String {
    let file = std::fs::File::open(zip_path).expect("report archive");
    let mut archive = zip::ZipArchive::new(file).expect("valid zip");
    let mut data = String::new();
    archive
        .by_name(&format!("{stem}_fastqc/fastqc_data.txt"))
        .expect("fastqc_data.txt in archive")
        .read_to_string(&mut data)
        .expect("readable");
    let section = data
        .split(">>Per base sequence quality")
        .nth(1)
        .expect("per-base-quality section");
    section
        .lines()
        .skip_while(|l| !l.starts_with('#'))
        .find(|l| l.split('\t').next() == Some(base))
        .unwrap_or_else(|| panic!("no row for base {base}"))
        .split('\t')
        .nth(1)
        .expect("mean column")
        .to_string()
}

fn run_trim_galore(dir: &Path, extra: &[&str]) {
    let mut args: Vec<&str> = vec!["--phred64", "--fastqc", "-q", "0", "--length", "0"];
    args.extend_from_slice(extra);
    let out = dir.to_str().unwrap();
    args.extend_from_slice(&["--output_dir", out, "test_files/phred64_test.fastq"]);

    let status = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(&args)
        .status()
        .expect("trim_galore runs");
    assert!(status.success(), "trim_galore failed with {args:?}");
}

#[test]
fn phred64_reaches_the_fastqc_config_on_fastq_output() {
    let dir = tempdir().unwrap();
    run_trim_galore(dir.path(), &[]);

    let zip = dir.path().join("phred64_test_trimmed_fastqc.zip");
    let stem = "phred64_test_trimmed";

    assert_eq!(
        basic_stat(&zip, stem, "Encoding"),
        "Illumina 1.5",
        "a Phred+64 file read as Phred+33 inflates every mean by 31"
    );
    // The fixture is 24 x Q40 then 10 x Q2 under Phred+64.
    assert_eq!(per_base_mean(&zip, stem, "1"), "40.0");
    assert_eq!(per_base_mean(&zip, stem, "25"), "2.0");
}

#[test]
fn bam_format_encoding_wins_over_the_phred64_flag() {
    let dir = tempdir().unwrap();
    run_trim_galore(dir.path(), &["--output-format", "ubam"]);

    let zip = dir.path().join("phred64_test_trimmed_fastqc.zip");
    let stem = "phred64_test_trimmed";

    // BAM declares Phred+33 by construction and that wins over the flag
    // upstream; this pins the contract the BAM call sites rely on.
    assert_eq!(basic_stat(&zip, stem, "Encoding"), "Sanger / Illumina 1.9");
    assert_eq!(per_base_mean(&zip, stem, "1"), "40.0");
    assert_eq!(per_base_mean(&zip, stem, "25"), "2.0");
}

/// Copy the Phred+64 fixture under `name`, so lockstep mates can be built.
fn fixture_copy(dir: &Path, name: &str) -> std::path::PathBuf {
    let dst = dir.join(name);
    std::fs::copy("test_files/phred64_test.fastq", &dst).expect("fixture copy");
    dst
}

fn trim_galore(dir: &Path, args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(args)
        .args(["--output_dir", dir.to_str().unwrap()])
        .status()
        .expect("trim_galore runs");
    assert!(status.success(), "trim_galore failed with {args:?}");
}

#[test]
fn phred64_reaches_both_mates_on_paired_output() {
    let dir = tempdir().unwrap();
    let r1 = fixture_copy(dir.path(), "p64_R1.fastq");
    let r2 = fixture_copy(dir.path(), "p64_R2.fastq");

    trim_galore(
        dir.path(),
        &[
            "--phred64",
            "--fastqc",
            "-q",
            "0",
            "--length",
            "0",
            "--paired",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
        ],
    );

    // Both mates must be told, not just R1.
    for (zip_name, stem) in [
        ("p64_R1_val_1_fastqc.zip", "p64_R1_val_1"),
        ("p64_R2_val_2_fastqc.zip", "p64_R2_val_2"),
    ] {
        let zip = dir.path().join(zip_name);
        assert_eq!(basic_stat(&zip, stem, "Encoding"), "Illumina 1.5", "{stem}");
        assert_eq!(per_base_mean(&zip, stem, "1"), "40.0", "{stem}");
    }
}

#[test]
fn phred64_reaches_the_passthrough_carrier() {
    let dir = tempdir().unwrap();
    let r1 = fixture_copy(dir.path(), "pt_R1.fastq");
    let r2 = fixture_copy(dir.path(), "pt_R2.fastq");
    let i1 = fixture_copy(dir.path(), "pt_I1.fastq");

    trim_galore(
        dir.path(),
        &[
            "--phred64",
            "--fastqc",
            "-q",
            "0",
            "--length",
            "0",
            "--paired",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--passthrough",
            i1.to_str().unwrap(),
        ],
    );

    // The carrier is copied untrimmed, so it is Phred+64 like the mates.
    let zip = dir.path().join("pt_I1_passthrough_fastqc.zip");
    assert_eq!(
        basic_stat(&zip, "pt_I1_passthrough", "Encoding"),
        "Illumina 1.5"
    );
}

#[test]
fn phred64_reaches_fastqc_under_clump_only() {
    let dir = tempdir().unwrap();
    let input = fixture_copy(dir.path(), "cl_p64.fastq");

    // --clump_only does no quality arithmetic, so the output is the input's
    // Phred+64 bytes verbatim.
    trim_galore(
        dir.path(),
        &[
            "--clump_only",
            "--phred64",
            "--fastqc",
            input.to_str().unwrap(),
        ],
    );

    let zip = dir.path().join("cl_p64_clumped_fastqc.zip");
    assert_eq!(
        basic_stat(&zip, "cl_p64_clumped", "Encoding"),
        "Illumina 1.5"
    );
    assert_eq!(per_base_mean(&zip, "cl_p64_clumped", "1"), "40.0");
}

#[test]
fn phred64_reaches_both_mates_under_clump_only() {
    let dir = tempdir().unwrap();
    let r1 = fixture_copy(dir.path(), "cp_R1.fastq");
    let r2 = fixture_copy(dir.path(), "cp_R2.fastq");

    trim_galore(
        dir.path(),
        &[
            "--clump_only",
            "--phred64",
            "--fastqc",
            "--paired",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
        ],
    );

    for (zip_name, stem) in [
        ("cp_R1_clumped_1_fastqc.zip", "cp_R1_clumped_1"),
        ("cp_R2_clumped_2_fastqc.zip", "cp_R2_clumped_2"),
    ] {
        let zip = dir.path().join(zip_name);
        assert_eq!(basic_stat(&zip, stem, "Encoding"), "Illumina 1.5", "{stem}");
    }
}
