//! `--fastqc_raw` must analyse the input and name the sample distinctly.
//!
//! The whole design rests on one field: MultiQC takes its sample name from the
//! `Filename` row inside `fastqc_data.txt`, not from the report's own name, and
//! fastqc-rust writes that row from `path.file_name()` of whatever it was handed.
//! Pointing FastQC at the input directly would therefore name the raw sample
//! identically to the trimmed one in single-end runs, where MultiQC silently
//! drops the duplicate rather than merging it.
//!
//! So these tests read inside the archive. Asserting that `*_raw_fastqc.zip`
//! exists would pass just as well with FastQC pointed at the wrong file.

use std::io::Read;
use std::path::Path;
use std::process::Command;

use tempfile::tempdir;

/// One Basic Statistics row from inside the report archive.
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

fn run(dir: &Path, extra: &[&str]) -> std::process::Output {
    let mut args: Vec<&str> = vec!["--output_dir", dir.to_str().unwrap()];
    args.extend_from_slice(extra);
    Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(&args)
        .output()
        .expect("trim_galore runs")
}

#[test]
fn single_end_raw_report_names_the_sample_raw() {
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &["--fastqc_raw", "test_files/illumina_10K.fastq.gz"],
    );
    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let zip = dir.path().join("illumina_10K_raw_fastqc.zip");
    assert!(zip.exists(), "no raw report at {}", zip.display());

    // The field MultiQC keys on. Without the staged alias this reads
    // `illumina_10K.fastq.gz`, which collides with the trimmed sample name.
    assert_eq!(
        basic_stat(&zip, "illumina_10K_raw", "Filename"),
        "illumina_10K_raw.fastq.gz"
    );

    // And it really is the untrimmed input: the fixture is uniform 150 bp, so a
    // report computed on the trimmed output would show a range instead.
    assert_eq!(
        basic_stat(&zip, "illumina_10K_raw", "Sequence length"),
        "150",
        "raw report must describe the input, not the trimmed output"
    );
}

#[test]
fn the_staged_alias_does_not_survive_the_run() {
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &["--fastqc_raw", "test_files/illumina_10K.fastq.gz"],
    );
    assert!(out.status.success());
    assert!(
        !dir.path().join("illumina_10K_raw.fastq.gz").exists(),
        "the staged alias must be removed; a stranded symlink would trip the \
         next run's output-collision pre-flight"
    );
}

#[test]
fn paired_raw_reports_are_named_per_mate() {
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &[
            "--fastqc_raw",
            "--paired",
            "test_files/BS-seq_10K_R1.fastq.gz",
            "test_files/BS-seq_10K_R2.fastq.gz",
        ],
    );
    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    for mate in ["R1", "R2"] {
        let stem = format!("BS-seq_10K_{mate}_raw");
        let zip = dir.path().join(format!("{stem}_fastqc.zip"));
        assert!(zip.exists(), "no raw report for {mate}");
        assert_eq!(
            basic_stat(&zip, &stem, "Filename"),
            format!("BS-seq_10K_{mate}_raw.fastq.gz")
        );
    }
}

#[test]
fn raw_and_trimmed_reports_coexist_with_distinct_sample_names() {
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &[
            "--fastqc_raw",
            "--fastqc",
            "test_files/illumina_10K.fastq.gz",
        ],
    );
    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let raw = dir.path().join("illumina_10K_raw_fastqc.zip");
    let trimmed = dir.path().join("illumina_10K_trimmed_fastqc.zip");
    assert!(
        raw.exists() && trimmed.exists(),
        "both sets must be written"
    );

    // Distinct in the field MultiQC reads — the point of the whole feature.
    assert_ne!(
        basic_stat(&raw, "illumina_10K_raw", "Filename"),
        basic_stat(&trimmed, "illumina_10K_trimmed", "Filename"),
    );
}

#[test]
fn a_collision_refused_run_writes_no_raw_report() {
    // The pre-trim pass runs after the output-collision pre-flight, so a refused
    // run must create nothing. Two inputs with the same basename in different
    // directories collide on their trimmed output.
    let dir = tempdir().unwrap();
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    std::fs::create_dir_all(&a).unwrap();
    std::fs::create_dir_all(&b).unwrap();
    for sub in [&a, &b] {
        std::fs::copy(
            "test_files/illumina_10K.fastq.gz",
            sub.join("same_name.fastq.gz"),
        )
        .unwrap();
    }

    let out = run(
        dir.path(),
        &[
            "--fastqc_raw",
            a.join("same_name.fastq.gz").to_str().unwrap(),
            b.join("same_name.fastq.gz").to_str().unwrap(),
        ],
    );
    assert!(
        !out.status.success(),
        "colliding outputs must be refused; stdout: {}",
        String::from_utf8_lossy(&out.stdout)
    );
    assert!(
        !dir.path().join("same_name_raw_fastqc.zip").exists(),
        "a refused run must not leave a raw report behind"
    );
}

#[test]
fn a_bad_second_input_is_refused_before_any_report_is_written() {
    // `main` sanity-checks only `cli.input[0]` before dispatch, so the pre-trim
    // pass validates every input itself — otherwise R1's report is written and
    // only then does R2 turn out to be unusable.
    //
    // Colorspace is the detectable case: `sanity_check_any` reads the first
    // record, so it catches this (and empty input) without reading the file.
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &[
            "--fastqc_raw",
            "--paired",
            "test_files/BS-seq_10K_R1.fastq.gz",
            "test_files/colorspace_file.fastq",
        ],
    );
    assert!(!out.status.success(), "a colorspace mate must be refused");

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.to_lowercase().contains("colorspace"),
        "expected trim-galore's own colorspace message, got: {stderr}"
    );
    assert!(
        !dir.path().join("BS-seq_10K_R1_raw_fastqc.zip").exists(),
        "R1's report must not be written before R2 has been validated"
    );
}

#[test]
fn mid_file_truncation_surfaces_from_the_qc_pass() {
    // Scope of the guard above, pinned so nobody reads it as broader than it is.
    // `sanity_check_any` inspects the first record, so truncation deep in the
    // stream is invisible to it — Trim Galore has the same blind spot while
    // trimming. The run still fails; it just fails during the QC pass, after
    // earlier inputs have been reported on.
    let dir = tempdir().unwrap();
    let out = run(
        dir.path(),
        &[
            "--fastqc_raw",
            "--paired",
            "test_files/BS-seq_10K_R1.fastq.gz",
            "test_files/truncated.fq.gz",
        ],
    );
    assert!(!out.status.success(), "a truncated mate must still fail");
    assert!(
        dir.path().join("BS-seq_10K_R1_raw_fastqc.zip").exists(),
        "R1 is reported on before R2's truncation is discoverable; if this ever \
         stops being true the guard has genuinely widened and the comment above \
         is stale"
    );
}
