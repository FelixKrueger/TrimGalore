//! Binary-driven tests for the output-collision pre-flight
//! (issues [#216](https://github.com/FelixKrueger/TrimGalore/issues/216),
//! [#383](https://github.com/FelixKrueger/TrimGalore/issues/383)).
//!
//! Before #383, four dispatch paths looped over `cli.input` with no pre-flight:
//! single-end trim (FASTQ and uBAM) and `--hardtrim5/3` (FASTQ and uBAM). Two
//! inputs whose output stems collided produced one output file at exit 0, with
//! the first input's reads gone and *two* trimming reports each claiming 100 %
//! written — which is why it went unnoticed.
//!
//! Every rejection case here is paired with an acceptance case on the same
//! dispatch path and output format, and the acceptance cases assert content or
//! filename rather than mere existence: a pre-flight that rejected everything,
//! or a candidate list built from the wrong namer, would otherwise pass.
//!
//! Fixtures are built in-test as plain FASTQ so read IDs identify their source
//! file, and so no case needs a fixture path relative to the crate root — the
//! hardtrim cases set `current_dir`, under which `test_files/…` would not resolve.

use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// Canonicalised on purpose: the collision key is lexical (assumption A8), and on
/// macOS `std::env::temp_dir()` yields `/tmp/…` while the child's `getcwd` yields
/// `/private/tmp/…`. Without this the absolute-vs-relative case compares two
/// spellings that differ by a symlink the key cannot see, and passes vacuously.
fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_coll_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    std::fs::canonicalize(&d).unwrap_or(d)
}

/// 40 reads whose IDs all carry `prefix`, so survivors are attributable.
fn write_fastq(path: &Path, prefix: &str) {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let mut s = String::new();
    for i in 0..40 {
        s.push_str(&format!(
            "@{prefix}_read{i}\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    std::fs::write(path, s).unwrap();
}

/// Run from `cwd` (matters for the specialty modes, which name output into the
/// current directory) and return (success, stderr).
fn run_in(cwd: &Path, args: &[&str]) -> (bool, String) {
    let out = Command::new(binary())
        .current_dir(cwd)
        .args(args)
        .output()
        .expect("failed to run trim_galore");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

fn count_reads_from(path: &Path, prefix: &str) -> usize {
    let body =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    body.lines()
        .filter(|l| l.starts_with(&format!("@{prefix}_")))
        .count()
}

/// `dir` must contain exactly `expected`. For rejected runs (nothing written
/// beside the inputs — the no-`--output_dir` twin of "directory is empty") and
/// for accepted runs where the written set must equal the planned set (#391).
fn assert_dir_holds_only(dir: &Path, expected: &[&str]) {
    let mut found: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    found.sort();
    let mut want: Vec<String> = expected.iter().map(|s| s.to_string()).collect();
    want.sort();
    assert_eq!(
        found,
        want,
        "{} must hold exactly the expected files",
        dir.display()
    );
}

const DUP_MSG: &str = "would be written to the same file";
const ALIAS_MSG: &str = "which is also one of its inputs";
const PREFIX: &str = "Output path collision (case-insensitive, for APFS/NTFS safety)";

/// A rejected run must write nothing at all — not just no primary output. Two
/// reports beside one data file was #383's most misleading artifact.
fn assert_rejected_cleanly(dir: &Path, ok: bool, stderr: &str, expect: &str) {
    assert!(
        !ok,
        "expected a non-zero exit, got success. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains(PREFIX),
        "missing collision prefix:\n{stderr}"
    );
    assert!(stderr.contains(expect), "expected {expect:?} in:\n{stderr}");
    // `.expect`, not `unwrap_or_default`: a read_dir failure would otherwise yield an
    // empty list and pass the assertion vacuously.
    let leftovers: Vec<String> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()))
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    assert!(
        leftovers.is_empty(),
        "a rejected run must write nothing, found {leftovers:?}"
    );
}

// ── single-end trim, FASTQ ────────────────────────────────────────────────

/// The reported bug: one stem, two FASTQ extensions.
#[test]
fn se_trim_rejects_shared_stem() {
    let dir = tempdir("se_stem");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    let (ok, stderr) = run_in(
        &dir,
        &["-o", out.to_str().unwrap(), "sample.fastq", "sample.fq"],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// REGRESSION (#383). A `./` prefix on one argument used to defeat the guard,
/// because the key was a raw string while both paths named one file.
#[test]
fn se_trim_rejects_dot_slash_alias() {
    let dir = tempdir("se_dot");
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    // No `-o`: the output path inherits the input's spelling, which is what makes
    // the two keys differ textually while naming one file.
    let (ok, stderr) = run_in(&dir, &["./sample.fastq", "sample.fq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(DUP_MSG),
        "expected duplicate wording:\n{stderr}"
    );
    assert_dir_holds_only(&dir, &["sample.fastq", "sample.fq"]);
}

/// REGRESSION (#383). Same defect through a mixed absolute/relative list, the
/// shape `find`/`xargs` and pipeline staging produce.
#[test]
fn se_trim_rejects_absolute_versus_relative_alias() {
    let dir = tempdir("se_abs");
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    let abs = dir.join("sample.fastq");
    // No `-o`, for the same reason as the `./` case above.
    let (ok, stderr) = run_in(&dir, &[abs.to_str().unwrap(), "sample.fq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(DUP_MSG),
        "expected duplicate wording:\n{stderr}"
    );
    assert_dir_holds_only(&dir, &["sample.fastq", "sample.fq"]);
}

/// The one slip nothing else would catch: passing `None` where `output_dir`
/// belongs. Without `-o` these two inputs do NOT collide (they land beside their
/// own inputs); with `-o` they must.
#[test]
fn se_trim_rejects_same_basename_across_dirs_with_output_dir() {
    let dir = tempdir("se_odir");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("dirA/same.fastq"), "DIRA");
    write_fastq(&dir.join("dirB/same.fastq"), "DIRB");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "-o",
            out.to_str().unwrap(),
            "dirA/same.fastq",
            "dirB/same.fastq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// REGRESSION (#382 widening). `.bgz` began sharing a stem with `.gz` when the
/// suffix stripping was fixed, turning a naming inconsistency into data loss.
#[test]
fn se_trim_rejects_gz_and_bgz_sharing_a_stem() {
    let dir = tempdir("se_bgz");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    // Content is irrelevant: the pre-flight runs before any read.
    write_fastq(&dir.join("wide.fastq.gz"), "GZA");
    write_fastq(&dir.join("wide.fastq.bgz"), "BGZB");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "-o",
            out.to_str().unwrap(),
            "wide.fastq.gz",
            "wide.fastq.bgz",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// REGRESSION (#383). A planned output that IS an input: the run would destroy
/// that input's reads and emit a second file holding the wrong sample.
#[test]
fn se_trim_rejects_output_that_aliases_an_input() {
    let dir = tempdir("se_alias");
    write_fastq(&dir.join("s.fastq"), "SRC");
    write_fastq(&dir.join("s_trimmed.fq"), "PRIOR");
    let (ok, stderr) = run_in(&dir, &["s.fastq", "s_trimmed.fq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert!(
        !stderr.contains(DUP_MSG),
        "must not use the duplicate-output wording:\n{stderr}"
    );
    // The named input must be untouched.
    assert_eq!(count_reads_from(&dir.join("s_trimmed.fq"), "PRIOR"), 40);
    assert_eq!(count_reads_from(&dir.join("s_trimmed.fq"), "SRC"), 0);
}

/// ACCEPTANCE. Same basename in different directories, no `-o`: outputs land
/// beside their own inputs and are genuinely distinct. Asserts attribution, not
/// existence — two files existed while #383 was live.
#[test]
fn se_trim_accepts_same_basename_across_dirs() {
    let dir = tempdir("se_ok_dirs");
    write_fastq(&dir.join("dirA/same.fastq"), "DIRA");
    write_fastq(&dir.join("dirB/same.fastq"), "DIRB");
    let (ok, stderr) = run_in(&dir, &["dirA/same.fastq", "dirB/same.fastq"]);
    assert!(ok, "expected success, got failure:\n{stderr}");
    let a = dir.join("dirA/same_trimmed.fq");
    let b = dir.join("dirB/same_trimmed.fq");
    assert_eq!(
        count_reads_from(&a, "DIRA"),
        40,
        "dirA must hold its own reads"
    );
    assert_eq!(
        count_reads_from(&a, "DIRB"),
        0,
        "dirA must not hold dirB's reads"
    );
    assert_eq!(
        count_reads_from(&b, "DIRB"),
        40,
        "dirB must hold its own reads"
    );
    assert_eq!(
        count_reads_from(&b, "DIRA"),
        0,
        "dirB must not hold dirA's reads"
    );
}

/// ACCEPTANCE. A single input must never trip the check.
#[test]
fn se_trim_accepts_single_input() {
    let dir = tempdir("se_ok_one");
    write_fastq(&dir.join("only.fastq"), "ONLY");
    let (ok, stderr) = run_in(&dir, &["only.fastq"]);
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(count_reads_from(&dir.join("only_trimmed.fq"), "ONLY"), 40);
}

// ── single-end trim, uBAM output ──────────────────────────────────────────

#[test]
fn se_trim_ubam_rejects_shared_stem() {
    let dir = tempdir("ubam_stem");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--output-format",
            "ubam",
            "-o",
            out.to_str().unwrap(),
            "sample.fastq",
            "sample.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// ACCEPTANCE. Catches an over-rejecting candidate list — e.g. one built from
/// the input paths, or from the FASTQ namer on the BAM arm.
#[test]
fn se_trim_ubam_accepts_distinct_stems() {
    let dir = tempdir("ubam_ok");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("x.fastq"), "XX");
    write_fastq(&dir.join("y.fastq"), "YY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--output-format",
            "ubam",
            "-o",
            out.to_str().unwrap(),
            "x.fastq",
            "y.fastq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert!(out.join("x_trimmed.bam").is_file(), "missing x_trimmed.bam");
    assert!(out.join("y_trimmed.bam").is_file(), "missing y_trimmed.bam");
}

// ── --hardtrim5 / --hardtrim3 ─────────────────────────────────────────────
//
// These name output into the CWD, so their collision class is wider than SE
// trim's: identical basenames collide even from different directories, and
// `--output_dir` does not rescue it. Hence the mode-specific hint.

#[test]
fn hardtrim5_rejects_same_basename_across_dirs() {
    let dir = tempdir("ht5");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("dirA/same.fastq"), "DIRA");
    write_fastq(&dir.join("dirB/same.fastq"), "DIRB");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim5",
            "20",
            "-o",
            out.to_str().unwrap(),
            "dirA/same.fastq",
            "dirB/same.fastq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
    assert!(
        stderr.contains("one invocation per input"),
        "hardtrim needs its own remediation — the generic advice does not work here:\n{stderr}"
    );
}

#[test]
fn hardtrim3_rejects_same_basename_across_dirs() {
    let dir = tempdir("ht3");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("dirA/same.fastq"), "DIRA");
    write_fastq(&dir.join("dirB/same.fastq"), "DIRB");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim3",
            "20",
            "-o",
            out.to_str().unwrap(),
            "dirA/same.fastq",
            "dirB/same.fastq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

#[test]
fn hardtrim5_ubam_rejects_shared_stem() {
    let dir = tempdir("ht5_ubam");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim5",
            "20",
            "--output-format",
            "ubam",
            "-o",
            out.to_str().unwrap(),
            "sample.fastq",
            "sample.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

#[test]
fn hardtrim3_ubam_rejects_shared_stem() {
    let dir = tempdir("ht3_ubam");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("sample.fastq"), "ALPHA");
    write_fastq(&dir.join("sample.fq"), "BETA");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim3",
            "20",
            "--output-format",
            "ubam",
            "-o",
            out.to_str().unwrap(),
            "sample.fastq",
            "sample.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// ACCEPTANCE, and the only coverage in the repository that `--hardtrim5`
/// produces `*.{N}bp_5prime.*` at all. Asserting the filename is what catches a
/// candidate list built with the wrong end discriminator.
#[test]
fn hardtrim5_accepts_distinct_stems_and_names_output() {
    let dir = tempdir("ht5_ok");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("x.fastq"), "XX");
    write_fastq(&dir.join("y.fastq"), "YY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim5",
            "20",
            "-o",
            out.to_str().unwrap(),
            "x.fastq",
            "y.fastq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert!(
        out.join("x.20bp_5prime.fq").is_file(),
        "missing x.20bp_5prime.fq"
    );
    assert!(
        out.join("y.20bp_5prime.fq").is_file(),
        "missing y.20bp_5prime.fq"
    );
}

/// ACCEPTANCE. `--hardtrim3` had no output-producing coverage before #383, and
/// its guard is a copy of `--hardtrim5`'s — so the `3prime` filename is the
/// assertion that matters.
#[test]
fn hardtrim3_accepts_distinct_stems_and_names_output() {
    let dir = tempdir("ht3_ok");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("x.fastq"), "XX");
    write_fastq(&dir.join("y.fastq"), "YY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim3",
            "20",
            "-o",
            out.to_str().unwrap(),
            "x.fastq",
            "y.fastq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert!(
        out.join("x.20bp_3prime.fq").is_file(),
        "missing x.20bp_3prime.fq"
    );
    assert!(
        out.join("y.20bp_3prime.fq").is_file(),
        "missing y.20bp_3prime.fq"
    );
}

#[test]
fn hardtrim3_ubam_accepts_distinct_stems_and_names_output() {
    let dir = tempdir("ht3_ubam_ok");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("x.fastq"), "XX");
    write_fastq(&dir.join("y.fastq"), "YY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--hardtrim3",
            "20",
            "--output-format",
            "ubam",
            "-o",
            out.to_str().unwrap(),
            "x.fastq",
            "y.fastq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert!(
        out.join("x.20bp_3prime.bam").is_file(),
        "missing x.20bp_3prime.bam"
    );
    assert!(
        out.join("y.20bp_3prime.bam").is_file(),
        "missing y.20bp_3prime.bam"
    );
}

// ── --paired: gained the input-alias check by sharing the helper ──────────

/// REGRESSION (#383). `--paired` had the #216 duplicate-output pre-flight since
/// 2.0 but no input-alias check, so a second pair naming the first pair's
/// prospective outputs destroyed them.
#[test]
fn paired_rejects_output_that_aliases_an_input() {
    let dir = tempdir("pe_alias");
    write_fastq(&dir.join("a_R1.fq"), "NEW1");
    write_fastq(&dir.join("a_R2.fq"), "NEW2");
    write_fastq(&dir.join("a_R1_val_1.fq"), "OLD1");
    write_fastq(&dir.join("a_R2_val_2.fq"), "OLD2");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "a_R1.fq",
            "a_R2.fq",
            "a_R1_val_1.fq",
            "a_R2_val_2.fq",
        ],
    );
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert_eq!(count_reads_from(&dir.join("a_R1_val_1.fq"), "OLD1"), 40);
}

/// ACCEPTANCE. Regression guard on the converted `--paired` pre-flight: two
/// ordinary pairs must still run to completion.
#[test]
fn paired_accepts_two_distinct_pairs() {
    let dir = tempdir("pe_ok");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    for (stem, tag) in [("p", "P"), ("q", "Q")] {
        write_fastq(&dir.join(format!("{stem}_R1.fq")), &format!("{tag}1"));
        write_fastq(&dir.join(format!("{stem}_R2.fq")), &format!("{tag}2"));
    }
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "-o",
            out.to_str().unwrap(),
            "p_R1.fq",
            "p_R2.fq",
            "q_R1.fq",
            "q_R2.fq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    for f in [
        "p_R1_val_1.fq",
        "p_R2_val_2.fq",
        "q_R1_val_1.fq",
        "q_R2_val_2.fq",
    ] {
        assert!(out.join(f).is_file(), "missing {f}");
    }
}

// ── duplicate positional (rejected in Cli::validate, not the pre-flight) ──

#[test]
fn duplicate_input_gets_a_precise_message() {
    let dir = tempdir("dup");
    write_fastq(&dir.join("dup.fastq"), "DUP");
    let (ok, stderr) = run_in(&dir, &["dup.fastq", "dup.fastq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains("was given more than once"),
        "expected the precise duplicate-input message:\n{stderr}"
    );
    assert!(
        !stderr.contains("APFS/NTFS"),
        "must not defer to the collision message:\n{stderr}"
    );
    // Rejected in validate(), before ensure_output_dir — nothing written.
    assert!(!dir.join("dup_trimmed.fq").exists());
}

// ── the widest hardtrim class: no `--output_dir`, so output lands in the CWD ──

/// The shape `--hardtrim5 30 */*.fastq.gz` takes. Every other hardtrim test passes
/// `-o`, which is the narrower class; this is the one that loses data today.
#[test]
fn hardtrim5_rejects_same_basename_across_dirs_without_output_dir() {
    let dir = tempdir("ht5_nocwd");
    write_fastq(&dir.join("dirA/same.fastq"), "DIRA");
    write_fastq(&dir.join("dirB/same.fastq"), "DIRB");
    let (ok, stderr) = run_in(
        &dir,
        &["--hardtrim5", "20", "dirA/same.fastq", "dirB/same.fastq"],
    );
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(DUP_MSG),
        "expected duplicate wording:\n{stderr}"
    );
    assert!(
        stderr.contains("current working directory"),
        "the CWD hint must replace the generic advice here:\n{stderr}"
    );
    assert!(
        !stderr.contains("different source directories"),
        "generic advice is false for this mode and must not appear:\n{stderr}"
    );
    // Nothing written into the CWD the run was launched from.
    assert!(!dir.join("same.20bp_5prime.fq").exists());
}

/// `--clock` names output into the CWD exactly as hardtrim does, so it needs the
/// same hint — `-o` cannot rescue the collision.
#[test]
fn clock_collision_gets_the_cwd_hint_not_the_false_advice() {
    let dir = tempdir("clock_hint");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    for d in ["dirA", "dirB"] {
        write_fastq(&dir.join(format!("{d}/same_R1.fq")), "R1");
        write_fastq(&dir.join(format!("{d}/same_R2.fq")), "R2");
    }
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clock",
            "-o",
            out.to_str().unwrap(),
            "dirA/same_R1.fq",
            "dirA/same_R2.fq",
            "dirB/same_R1.fq",
            "dirB/same_R2.fq",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains("current working directory"),
        "missing hint:\n{stderr}"
    );
    assert!(
        !stderr.contains("different source directories"),
        "the user already supplied -o; that advice is false here:\n{stderr}"
    );
}

// ── secondary outputs may not overwrite an input (all three reproduced routes) ──

/// A previous run's trimming report fed back in as an input.
#[test]
fn se_trim_rejects_output_that_aliases_a_report_input() {
    let dir = tempdir("alias_report");
    write_fastq(&dir.join("sample.fastq"), "SRC");
    write_fastq(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--dont_gzip",
            "sample.fastq",
            "sample.fastq_trimming_report.txt",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert_eq!(
        count_reads_from(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY"),
        40
    );
}

/// `--demux`'s barcode file, named like the prospective trim output.
#[test]
fn demux_rejects_output_that_aliases_the_barcode_file() {
    let dir = tempdir("alias_bc");
    write_fastq(&dir.join("sample.fastq"), "SRC");
    std::fs::write(
        dir.join("sample_trimmed.fq"),
        "sample1\tACGT\nsample2\tTGCA\n",
    )
    .unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--dont_gzip",
            "--demux",
            "sample_trimmed.fq",
            "sample.fastq",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    let bc = std::fs::read_to_string(dir.join("sample_trimmed.fq")).unwrap();
    assert!(
        bc.starts_with("sample1\t"),
        "barcode file was modified: {bc:?}"
    );
}

/// A previous demux run's per-barcode output fed back in as an input. This is the
/// route that needs the *secondary* output paths in the candidate list — every
/// primary stays distinct here.
#[test]
fn demux_rejects_output_that_aliases_a_per_barcode_input() {
    let dir = tempdir("alias_bcout");
    write_fastq(&dir.join("sample.fastq"), "SRC");
    write_fastq(&dir.join("sample_trimmed_sample1.fq"), "PRECIOUS");
    std::fs::write(dir.join("bc.txt"), "sample1\tACGT\nsample2\tTGCA\n").unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--dont_gzip",
            "--demux",
            "bc.txt",
            "sample.fastq",
            "sample_trimmed_sample1.fq",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert_eq!(
        count_reads_from(&dir.join("sample_trimmed_sample1.fq"), "PRECIOUS"),
        40
    );
}

/// REGRESSION (#383). `..` used to defeat the key: the reported bug reproduced
/// through its own fix whenever one argument carried a `..` component.
#[test]
fn se_trim_rejects_dotdot_alias() {
    let dir = tempdir("dotdot");
    write_fastq(&dir.join("data/s.fastq"), "ALPHA");
    write_fastq(&dir.join("data/s.fq"), "BETA");
    std::fs::create_dir_all(dir.join("work")).unwrap();
    let abs = dir.join("data/s.fq");
    let (ok, stderr) = run_in(
        &dir.join("work"),
        &["../data/s.fastq", abs.to_str().unwrap()],
    );
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(DUP_MSG),
        "expected duplicate wording:\n{stderr}"
    );
    assert!(
        !dir.join("data/s_trimmed.fq").exists(),
        "nothing may be written"
    );
}

/// `--paired` accepted one file as its own mate when the two spellings differed,
/// emitting two byte-identical files labelled a validated pair.
#[test]
fn paired_rejects_one_file_as_its_own_mate_across_spellings() {
    let dir = tempdir("selfpair");
    write_fastq(&dir.join("a_R1.fq"), "PAIR");
    let (ok, stderr) = run_in(&dir, &["--paired", "./a_R1.fq", "a_R1.fq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains("appear to be the same file"),
        "expected the precise R1==R2 message:\n{stderr}"
    );
    assert!(!dir.join("a_R1_val_1.fq").exists());
}

/// Pins `demux::demux_output_paths` against what `demultiplex` actually writes.
/// The pre-flight guards the planned set, so if the two ever drift the guard would
/// silently protect paths the run never touches — invisible to every rejection test.
#[test]
fn demux_writes_exactly_the_planned_paths() {
    let dir = tempdir("demux_pin");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("m.fastq"), "SRC");
    std::fs::write(dir.join("bc.txt"), "sample1\tACGT\nsample2\tTGCA\n").unwrap();

    let (ok, stderr) = run_in(
        &dir,
        &[
            "--dont_gzip",
            "--no_poly_g",
            "--demux",
            "bc.txt",
            "-o",
            out.to_str().unwrap(),
            "m.fastq",
        ],
    );
    assert!(ok, "demux run should succeed:\n{stderr}");

    let barcodes = trim_galore::demux::read_barcode_file(&dir.join("bc.txt")).unwrap();
    let trimmed = out.join("m_trimmed.fq");
    let mut planned: Vec<String> =
        trim_galore::demux::demux_output_paths(&trimmed, &barcodes, false, Some(out.as_path()))
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
    planned.sort();

    // Everything demux itself produced: the per-barcode files, NoCode, the summary.
    let mut actual: Vec<String> = std::fs::read_dir(&out)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("m_trimmed_"))
        .collect();
    actual.sort();

    assert_eq!(
        planned, actual,
        "demux_output_paths must predict exactly what demultiplex writes"
    );
    assert!(
        !planned.is_empty(),
        "precondition: the planned set must be non-empty"
    );
}

// ── #388: paired report paths join the pre-flight ─────────────────────────
//
// Paired primaries carry a positional discriminator (`_val_1`/`_val_2`) that
// report names do not, so two inputs with distinct primaries can collide on
// reports. Before #388 all three rejection shapes below exited 0 having
// silently overwritten one report pair.

/// T1 — the case-free shape: same filename as R1 AND R2 of one pair, shared
/// output dir. Runs identically on every filesystem.
#[test]
fn paired_rejects_same_filename_r1_r2_into_shared_output_dir() {
    let dir = tempdir("388_samename");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("A/reads.fq"), "R1SIDE");
    write_fastq(&dir.join("B/reads.fq"), "R2SIDE");
    let (ok, stderr) = run_in(&dir, &["--paired", "-o", "out", "A/reads.fq", "B/reads.fq"]);
    assert!(!ok, "expected rejection, got success:\n{stderr}");
    assert!(
        stderr.contains(DUP_MSG),
        "expected collision wording:\n{stderr}"
    );
    // #397 — the message must name BOTH inputs the user has to choose between.
    // Before provenance it printed the colliding path twice and named neither.
    assert!(
        stderr.contains("A/reads.fq") && stderr.contains("B/reads.fq"),
        "the refusal must name both source inputs:\n{stderr}"
    );
    assert!(
        std::fs::read_dir(&out).unwrap().next().is_none(),
        "a refused run must write nothing"
    );
}

/// #397 — one file passed as a mate of two different pairs. `validate` permits it
/// (only exact duplicate pairs are rejected), and the old message told the user to
/// "rename one input", which cannot be done for a file colliding with itself.
#[test]
fn paired_one_file_in_two_pairs_gets_followable_advice() {
    let dir = tempdir("397_same_source");
    write_fastq(&dir.join("A.fq"), "A");
    write_fastq(&dir.join("B.fq"), "B");
    write_fastq(&dir.join("C.fq"), "C");
    let (ok, stderr) = run_in(&dir, &["--paired", "A.fq", "B.fq", "B.fq", "C.fq"]);
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(stderr.contains(DUP_MSG), "got:\n{stderr}");
    assert!(
        stderr.contains("two outputs named from") && stderr.contains("B.fq"),
        "must attribute both candidates to the one input:\n{stderr}"
    );
    assert!(
        stderr.contains("List each input once"),
        "must give advice that can actually be followed:\n{stderr}"
    );
    assert!(
        !stderr.contains("rename one input"),
        "the unfollowable advice must be gone:\n{stderr}"
    );
}

/// T2 — fold-equal filenames as R1/R2. On a case-insensitive filesystem this is
/// a true positive; on a case-sensitive one it is the documented #216-style
/// false positive (loud error over silent loss). The REJECTION is asserted, so
/// the test is filesystem-independent.
#[test]
fn paired_rejects_fold_equal_filenames_into_shared_output_dir() {
    let dir = tempdir("388_foldeq");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("a/r1.fq"), "LOWER");
    write_fastq(&dir.join("b/R1.fq"), "UPPER");
    let (ok, stderr) = run_in(&dir, &["--paired", "-o", "out", "a/r1.fq", "b/R1.fq"]);
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// T3 — the cross-pair route: pair 2 reuses a filename from pair 1.
#[test]
fn paired_rejects_cross_pair_report_collision() {
    let dir = tempdir("388_crosspair");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("a/r1.fq"), "P1R1");
    write_fastq(&dir.join("a/r2.fq"), "P1R2");
    write_fastq(&dir.join("b/R2.fq"), "P2R1");
    write_fastq(&dir.join("b/x.fq"), "P2R2");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired", "-o", "out", "a/r1.fq", "a/r2.fq", "b/R2.fq", "b/x.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// T4 — the uBAM-output twin had the identical hole.
#[test]
fn paired_ubam_rejects_same_filename_r1_r2_into_shared_output_dir() {
    let dir = tempdir("388_ubam");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("A/reads.fq"), "R1SIDE");
    write_fastq(&dir.join("B/reads.fq"), "R2SIDE");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "--output-format",
            "ubam",
            "-o",
            "out",
            "A/reads.fq",
            "B/reads.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
}

/// T5 — the gate: with --no_report_file no reports will be written, so the
/// same inputs must run to completion. Pins candidates == writers.
#[test]
fn paired_same_filename_accepted_when_reports_disabled() {
    let dir = tempdir("388_gate");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    write_fastq(&dir.join("A/reads.fq"), "R1SIDE");
    write_fastq(&dir.join("B/reads.fq"), "R2SIDE");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "--no_report_file",
            "-o",
            "out",
            "A/reads.fq",
            "B/reads.fq",
        ],
    );
    assert!(ok, "expected success with reports disabled:\n{stderr}");
    assert!(out.join("reads_val_1.fq").is_file());
    assert!(out.join("reads_val_2.fq").is_file());
    assert_eq!(
        std::fs::read_dir(&out).unwrap().count(),
        2,
        "no reports may be written under --no_report_file"
    );
}

/// T6 — over-rejection guards at the boundary A1 describes.
#[test]
fn paired_report_candidates_do_not_over_reject() {
    // Single pair + --basename, no -o: reports keep per-input names beside
    // their own inputs; primaries take the basename. Nothing collides.
    let dir = tempdir("388_ok_basename");
    write_fastq(&dir.join("r1.fq"), "B1");
    write_fastq(&dir.join("r2.fq"), "B2");
    let (ok, stderr) = run_in(&dir, &["--paired", "--basename", "foo", "r1.fq", "r2.fq"]);
    assert!(ok, "basename pair must still run:\n{stderr}");
    assert!(dir.join("foo_val_1.fq").is_file());
    assert!(dir.join("r1.fq_trimming_report.txt").is_file());
    assert!(dir.join("r2.fq_trimming_report.txt").is_file());

    // Same stem, different extension as R1/R2: report keys differ (full
    // filename), primaries differ (_val_1/_val_2). Nothing collides.
    let dir2 = tempdir("388_ok_stem");
    write_fastq(&dir2.join("sample.fq"), "S1");
    write_fastq(&dir2.join("sample.fastq"), "S2");
    let (ok, stderr) = run_in(&dir2, &["--paired", "sample.fq", "sample.fastq"]);
    assert!(
        ok,
        "same-stem different-extension pair must still run:\n{stderr}"
    );
    assert!(dir2.join("sample.fq_trimming_report.txt").is_file());
    assert!(dir2.join("sample.fastq_trimming_report.txt").is_file());

    // #398 — T1 minus -o is now REFUSED: both reports resolve into R1's directory,
    // where their identical filenames collide. This block previously asserted the
    // split layout; co-locating reports is now the specified behaviour. Do NOT
    // narrow the candidate list to make it green again — that restores the old
    // paths on the pre-flight side while the writers move. Its clump twin is
    // `clump_report_candidates_reject_shared_filename_without_output_dir`.
    let dir3 = tempdir("388_reject_no_odir");
    write_fastq(&dir3.join("A/reads.fq"), "N1");
    write_fastq(&dir3.join("B/reads.fq"), "N2");
    let (ok, stderr) = run_in(&dir3, &["--paired", "A/reads.fq", "B/reads.fq"]);
    assert!(!ok, "co-located reports must be refused:\n{stderr}");
    assert!(
        stderr.contains("reads.fq_trimming_report.txt"),
        "the colliding path must be named:\n{stderr}"
    );
    assert_dir_holds_only(&dir3.join("A"), &["reads.fq"]);
    assert_dir_holds_only(&dir3.join("B"), &["reads.fq"]);
}

// ── --clump_only: clumping reports join the pre-flight (#391) ─────────────

/// The #391 repro: same filename mates, `-o` collecting both reports onto one
/// path. Primaries carry `_clumped_1`/`_clumped_2` and never collide; the
/// reports are the only collision, and the message must say so.
#[test]
fn clump_paired_rejects_shared_report_name() {
    let dir = tempdir("clump_rep");
    write_fastq(&dir.join("a/reads.fq"), "CR1");
    write_fastq(&dir.join("b/reads.fq"), "CR2");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "-o",
            "out",
            "a/reads.fq",
            "b/reads.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
    assert!(
        stderr.contains("reads.fq_clumping_report.txt"),
        "the colliding path must be the clumping report:\n{stderr}"
    );
    assert!(
        stderr.contains("--no_report_file"),
        "hint must offer the report-only remedy:\n{stderr}"
    );
}

/// Same fixtures, `--no_report_file`: the only paths the flag removes from the
/// candidate list are the two reports, so success here proves the rejection
/// above came from the reports — and that the gate matches the writer.
#[test]
fn clump_paired_accepts_shared_report_name_with_no_report_file() {
    let dir = tempdir("clump_rep_norep");
    write_fastq(&dir.join("a/reads.fq"), "CN1");
    write_fastq(&dir.join("b/reads.fq"), "CN2");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "--no_report_file",
            "-o",
            "out",
            "a/reads.fq",
            "b/reads.fq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(count_reads_from(&out.join("reads_clumped_1.fq"), "CN1"), 40);
    assert_eq!(count_reads_from(&out.join("reads_clumped_2.fq"), "CN2"), 40);
    assert_dir_holds_only(&out, &["reads_clumped_1.fq", "reads_clumped_2.fq"]);
}

/// Acceptance guard: distinct filenames → primaries + reports all written, and
/// nothing else. The exact-set assertion pins that the pre-flight's candidate
/// list and the writer's file set agree — the invariant whose absence is this
/// bug family's root cause (#388, #391).
#[test]
fn clump_paired_accepts_distinct_filenames() {
    let dir = tempdir("clump_ok");
    write_fastq(&dir.join("a/r1.fq"), "D1");
    write_fastq(&dir.join("b/r2.fq"), "D2");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "-o",
            "out",
            "a/r1.fq",
            "b/r2.fq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(count_reads_from(&out.join("r1_clumped_1.fq"), "D1"), 40);
    assert_eq!(count_reads_from(&out.join("r2_clumped_2.fq"), "D2"), 40);
    assert_dir_holds_only(
        &out,
        &[
            "r1_clumped_1.fq",
            "r2_clumped_2.fq",
            "r1.fq_clumping_report.txt",
            "r2.fq_clumping_report.txt",
        ],
    );
}

/// Cross-pair: all four primaries distinct (the positional suffixes cross),
/// reports collide across pairs. Pins that the candidate list spans pairs.
/// No filename assertion: candidate order decides which report is flagged.
#[test]
fn clump_paired_rejects_cross_pair_report_collision() {
    let dir = tempdir("clump_xpair");
    write_fastq(&dir.join("a/x.fq"), "X1");
    write_fastq(&dir.join("a/y.fq"), "Y1");
    write_fastq(&dir.join("b/y.fq"), "Y2");
    write_fastq(&dir.join("b/x.fq"), "X2");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "-o",
            "out",
            "a/x.fq",
            "a/y.fq",
            "b/y.fq",
            "b/x.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
    assert!(
        stderr.contains("_clumping_report.txt"),
        "whichever report is flagged first, the class is order-independent:\n{stderr}"
    );
}

/// Fold-equal report names: the primaries differ even case-folded (positional
/// suffixes), so the reports are the only collision, and only under the fold.
/// The REJECTION is asserted, so the test is filesystem-independent.
#[test]
fn clump_paired_rejects_fold_equal_report_names() {
    let dir = tempdir("clump_foldeq");
    write_fastq(&dir.join("a/Reads.fq"), "FU");
    write_fastq(&dir.join("b/reads.fq"), "FL");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "-o",
            "out",
            "a/Reads.fq",
            "b/reads.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
    assert!(
        stderr.contains("_clumping_report.txt"),
        "the fold-equal pair must be the reports:\n{stderr}"
    );
}

/// #398 — same filename in two dirs WITHOUT `-o` is refused: both reports now
/// resolve into R1's directory, where their identical filenames collide. The
/// primaries were never the problem (`_clumped_1`/`_clumped_2` keep them apart).
///
/// This is the clump twin of `paired_report_candidates_do_not_over_reject`'s third
/// block. Both were written to catch a builder that co-located reports; #398 makes
/// co-location the specified behaviour, so both assert the refusal instead. If this
/// goes red again, fix the writers or the candidates — do NOT narrow the candidate
/// list, which would restore the split layout on one side only.
#[test]
fn clump_report_candidates_reject_shared_filename_without_output_dir() {
    let dir = tempdir("clump_reject_no_odir");
    write_fastq(&dir.join("A/reads.fq"), "Q1");
    write_fastq(&dir.join("B/reads.fq"), "Q2");
    let (ok, stderr) = run_in(
        &dir,
        &["--clump_only", "--paired", "A/reads.fq", "B/reads.fq"],
    );
    assert!(!ok, "co-located reports must be refused:\n{stderr}");
    assert!(
        stderr.contains("reads.fq_clumping_report.txt"),
        "the colliding path must be named:\n{stderr}"
    );
    // Refused before any writer: both directories hold only their inputs.
    assert_dir_holds_only(&dir.join("A"), &["reads.fq"]);
    assert_dir_holds_only(&dir.join("B"), &["reads.fq"]);
}

/// #398 — the same pair with `--no_report_file` runs, because only the reports
/// collided. Keeps the acceptance side of the guard alive.
#[test]
fn clump_shared_filename_without_output_dir_runs_without_reports() {
    let dir = tempdir("clump_no_report_no_odir");
    write_fastq(&dir.join("A/reads.fq"), "Q1");
    write_fastq(&dir.join("B/reads.fq"), "Q2");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "--no_report_file",
            "A/reads.fq",
            "B/reads.fq",
        ],
    );
    assert!(ok, "only the reports collided:\n{stderr}");
    assert_eq!(
        count_reads_from(&dir.join("A/reads_clumped_1.fq"), "Q1"),
        40
    );
    assert_eq!(
        count_reads_from(&dir.join("A/reads_clumped_2.fq"), "Q2"),
        40
    );
    assert_dir_holds_only(
        &dir.join("A"),
        &["reads.fq", "reads_clumped_1.fq", "reads_clumped_2.fq"],
    );
    assert_dir_holds_only(&dir.join("B"), &["reads.fq"]);
}

/// #398 — one file as the mate of two pairs, no `-o`. This was refused before
/// #398, because both pairs' reports for the shared mate resolved to
/// `d/x.fq_clumping_report.txt`. Each pair now anchors on its own R1, so the two
/// reports land in `p/` and `q/` and all eight paths are distinct. Nothing is
/// overwritten, so continuing to refuse would be over-rejection.
///
/// The refusal this replaces was added by #391; the flip is deliberate, not a
/// weakened assertion.
#[test]
fn clump_paired_shared_mate_report_runs_when_pairs_anchor_apart() {
    let dir = tempdir("clump_mate");
    write_fastq(&dir.join("p/a.fq"), "M1");
    write_fastq(&dir.join("d/x.fq"), "MX");
    write_fastq(&dir.join("q/b.fq"), "M2");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "p/a.fq",
            "d/x.fq",
            "q/b.fq",
            "d/x.fq",
        ],
    );
    assert!(ok, "eight distinct paths must not be refused:\n{stderr}");
    // Pair 1 → p/, pair 2 → q/. The shared mate's own directory gains nothing.
    assert_dir_holds_only(
        &dir.join("p"),
        &[
            "a.fq",
            "a_clumped_1.fq",
            "x_clumped_2.fq",
            "a.fq_clumping_report.txt",
            "x.fq_clumping_report.txt",
        ],
    );
    assert_dir_holds_only(
        &dir.join("q"),
        &[
            "b.fq",
            "b_clumped_1.fq",
            "x_clumped_2.fq",
            "b.fq_clumping_report.txt",
            "x.fq_clumping_report.txt",
        ],
    );
    assert_dir_holds_only(&dir.join("d"), &["x.fq"]);
}

/// `--basename foo` forces `foo_clumped_1`/`foo_clumped_2` primaries — distinct
/// by construction — so this isolates the report path without any reasoning
/// about stems. Reports ignore `--basename` and still collide.
#[test]
fn clump_paired_rejects_shared_report_name_with_basename() {
    let dir = tempdir("clump_base");
    write_fastq(&dir.join("a/reads.fq"), "B1");
    write_fastq(&dir.join("b/reads.fq"), "B2");
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "--basename",
            "foo",
            "-o",
            "out",
            "a/reads.fq",
            "b/reads.fq",
        ],
    );
    assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG);
    assert!(
        stderr.contains("reads.fq_clumping_report.txt"),
        "with --basename the reports are the only possible collision:\n{stderr}"
    );
}

/// #385 defect-2 class on the SE clump arm: an input named like another input's
/// clumping report was silently overwritten (this run exited 0 and destroyed
/// input 2 before #391). Both files are valid FASTQ so the rejection is
/// attributable to the pre-flight, not a parse error.
#[test]
fn clump_se_rejects_report_that_aliases_an_input() {
    let dir = tempdir("clump_se_alias");
    write_fastq(&dir.join("s.fq"), "SRC");
    write_fastq(&dir.join("s.fq_clumping_report.txt"), "REPORTY");
    let (ok, stderr) = run_in(&dir, &["--clump_only", "s.fq", "s.fq_clumping_report.txt"]);
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert!(
        stderr.contains("s.fq_clumping_report.txt"),
        "the aliased path must be the clumping report:\n{stderr}"
    );
    assert_eq!(
        count_reads_from(&dir.join("s.fq_clumping_report.txt"), "REPORTY"),
        40
    );
    assert_dir_holds_only(&dir, &["s.fq", "s.fq_clumping_report.txt"]);
}

/// The same shape runs when reports are off — and writes exactly the primaries.
#[test]
fn clump_se_accepts_report_alias_with_no_report_file() {
    let dir = tempdir("clump_se_alias_ok");
    write_fastq(&dir.join("s.fq"), "SRC");
    write_fastq(&dir.join("s.fq_clumping_report.txt"), "REPORTY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--no_report_file",
            "s.fq",
            "s.fq_clumping_report.txt",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(count_reads_from(&dir.join("s_clumped.fq"), "SRC"), 40);
    assert_eq!(
        count_reads_from(&dir.join("s.fq_clumping_report_clumped.fq"), "REPORTY"),
        40
    );
    assert_dir_holds_only(
        &dir,
        &[
            "s.fq",
            "s.fq_clumping_report.txt",
            "s_clumped.fq",
            "s.fq_clumping_report_clumped.fq",
        ],
    );
}

/// Task 2's paired-BAM Shape A arm plans ONE report per pair, keyed on the
/// pair's FIRST input (clump_only.rs keys on inputs[0]). Keyed on chunk[1] by
/// mistake, pair 1 would plan y.fq_clumping_report.txt, collide with nothing,
/// and the run would proceed — this test discriminates exactly that.
/// Acceptance sibling for this dispatch path lives cross-file:
/// integration_clump_only_ubam.rs::multi_pair_pe_bam_produces_one_output_per_pair.
#[test]
fn clump_paired_bam_rejects_report_that_aliases_an_input() {
    let dir = tempdir("clump_bam_alias");
    write_fastq(&dir.join("x.fq"), "PX");
    write_fastq(&dir.join("y.fq"), "PY");
    write_fastq(&dir.join("x.fq_clumping_report.txt"), "PR");
    write_fastq(&dir.join("z.fq"), "PZ");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "x.fq",
            "y.fq",
            "x.fq_clumping_report.txt",
            "z.fq",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    assert_eq!(count_reads_from(&dir.join("x.fq"), "PX"), 40);
    assert_eq!(
        count_reads_from(&dir.join("x.fq_clumping_report.txt"), "PR"),
        40
    );
    assert_dir_holds_only(&dir, &["x.fq", "y.fq", "x.fq_clumping_report.txt", "z.fq"]);
}

// ── SE uBAM-output trim: reports join the pre-flight (#409) ───────────────

/// #409 — the uBAM twin of `se_trim_rejects_output_that_aliases_a_report_input`.
/// Before the fix this run overwrote input 2 with input 1's trimming report and
/// only then failed reading it, so the error blamed the input for not being
/// FASTQ when the run had just made that true.
#[test]
fn se_trim_ubam_rejects_output_that_aliases_a_report_input() {
    let dir = tempdir("409_ubam_alias_report");
    write_fastq(&dir.join("sample.fastq"), "SRC");
    write_fastq(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--output-format",
            "ubam",
            "sample.fastq",
            "sample.fastq_trimming_report.txt",
        ],
    );
    assert!(!ok, "expected rejection:\n{stderr}");
    assert!(
        stderr.contains(ALIAS_MSG),
        "expected alias wording:\n{stderr}"
    );
    // The data-loss assertion: the victim must still hold its own reads.
    assert_eq!(
        count_reads_from(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY"),
        40,
        "input 2 was overwritten — this is the #409 data loss"
    );
    assert_dir_holds_only(&dir, &["sample.fastq", "sample.fastq_trimming_report.txt"]);
}

/// Acceptance sibling: with reports off, nothing collides and the run proceeds.
#[test]
fn se_trim_ubam_accepts_report_alias_with_no_report_file() {
    let dir = tempdir("409_ubam_alias_ok");
    write_fastq(&dir.join("sample.fastq"), "SRC");
    write_fastq(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--output-format",
            "ubam",
            "--no_report_file",
            "sample.fastq",
            "sample.fastq_trimming_report.txt",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(
        count_reads_from(&dir.join("sample.fastq_trimming_report.txt"), "REPORTY"),
        40
    );
    assert!(dir.join("sample_trimmed.bam").is_file());
    assert!(
        dir.join("sample.fastq_trimming_report_trimmed.bam")
            .is_file()
    );
}

// ── #398: one output directory per pair ───────────────────────────────────

/// The layout the change is for: mates in different directories, no `-o`. Every
/// output — both primaries and both reports — lands in R1's directory.
#[test]
fn paired_reports_follow_the_primaries_into_r1s_directory() {
    let dir = tempdir("398_unified");
    write_fastq(&dir.join("A/r1.fq"), "U1");
    write_fastq(&dir.join("B/r2.fq"), "U2");
    let (ok, stderr) = run_in(&dir, &["--paired", "A/r1.fq", "B/r2.fq"]);
    assert!(ok, "distinct filenames must still run:\n{stderr}");
    assert_dir_holds_only(
        &dir.join("A"),
        &[
            "r1.fq",
            "r1_val_1.fq",
            "r2_val_2.fq",
            "r1.fq_trimming_report.txt",
            "r1.fq_trimming_report.json",
            "r2.fq_trimming_report.txt",
            "r2.fq_trimming_report.json",
        ],
    );
    assert_dir_holds_only(&dir.join("B"), &["r2.fq"]);
}

/// The mate-side directory layout (`R1/` and `R2/` rather than per-sample): the
/// whole invocation is refused, not just the offending pair, because candidates
/// for every chunk are collected before one pre-flight call.
#[test]
fn mate_side_directory_layout_refuses_the_whole_invocation() {
    let dir = tempdir("398_mate_side");
    write_fastq(&dir.join("R1/s1.fq"), "A1");
    write_fastq(&dir.join("R2/s1.fq"), "A2");
    write_fastq(&dir.join("R1/s2.fq"), "B1");
    write_fastq(&dir.join("R2/s2.fq"), "B2");
    let (ok, stderr) = run_in(
        &dir,
        &["--paired", "R1/s1.fq", "R2/s1.fq", "R1/s2.fq", "R2/s2.fq"],
    );
    assert!(!ok, "co-located reports must be refused:\n{stderr}");
    // No pair ran: not even the first pair's primaries exist.
    assert_dir_holds_only(&dir.join("R1"), &["s1.fq", "s2.fq"]);
    assert_dir_holds_only(&dir.join("R2"), &["s1.fq", "s2.fq"]);
}

/// The per-sample layout is unaffected — each pair's R1 sits in its own directory,
/// so nothing collides and each pair's outputs stay with it.
#[test]
fn per_sample_directory_layout_is_unaffected() {
    let dir = tempdir("398_per_sample");
    write_fastq(&dir.join("sA/reads_1.fq"), "S1");
    write_fastq(&dir.join("sA/reads_2.fq"), "S2");
    write_fastq(&dir.join("sB/reads_1.fq"), "T1");
    write_fastq(&dir.join("sB/reads_2.fq"), "T2");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "--no_report_file",
            "sA/reads_1.fq",
            "sA/reads_2.fq",
            "sB/reads_1.fq",
            "sB/reads_2.fq",
        ],
    );
    assert!(ok, "per-sample layout must run:\n{stderr}");
    // Each pair anchored on its OWN R1 — pair 2 did not land in sA/.
    assert_eq!(count_reads_from(&dir.join("sA/reads_1_val_1.fq"), "S1"), 40);
    assert_eq!(count_reads_from(&dir.join("sB/reads_1_val_1.fq"), "T1"), 40);
    assert!(!dir.join("sA/reads_1_val_1.fq.1").exists());
}

/// #398 creates a new output-vs-**input** path: R2's report, relocated into R1's
/// directory, can equal R1's own input. This is #409's shape on a new arm — if the
/// candidate list had not moved with the writer, the run would destroy an input.
#[test]
fn paired_report_that_would_land_on_r1s_input_is_refused() {
    let dir = tempdir("398_alias_input");
    write_fastq(&dir.join("A/reads.fq_trimming_report.txt"), "VICTIM");
    write_fastq(&dir.join("B/reads.fq"), "MATE");
    let (ok, stderr) = run_in(
        &dir,
        &["--paired", "A/reads.fq_trimming_report.txt", "B/reads.fq"],
    );
    assert!(
        !ok,
        "a report landing on an input must be refused:\n{stderr}"
    );
    // The data-loss assertion: content, not just existence — a truncated file exists.
    assert_eq!(
        count_reads_from(&dir.join("A/reads.fq_trimming_report.txt"), "VICTIM"),
        40,
        "R1's input was overwritten — this is #409's shape on the paired arm"
    );
    assert_dir_holds_only(&dir.join("A"), &["reads.fq_trimming_report.txt"]);
    assert_dir_holds_only(&dir.join("B"), &["reads.fq"]);
}

/// Acceptance sibling: with reports off there is no report to alias the input.
#[test]
fn paired_report_input_alias_accepted_with_no_report_file() {
    let dir = tempdir("398_alias_ok");
    write_fastq(&dir.join("A/reads.fq_trimming_report.txt"), "VICTIM");
    write_fastq(&dir.join("B/reads.fq"), "MATE");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--paired",
            "--no_report_file",
            "A/reads.fq_trimming_report.txt",
            "B/reads.fq",
        ],
    );
    assert!(ok, "expected success:\n{stderr}");
    assert_eq!(
        count_reads_from(&dir.join("A/reads.fq_trimming_report.txt"), "VICTIM"),
        40
    );
}

/// The same output-vs-input hazard on the `--clump_only --paired` arm.
#[test]
fn clump_paired_report_that_would_land_on_r1s_input_is_refused() {
    let dir = tempdir("398_clump_alias");
    write_fastq(&dir.join("A/reads.fq_clumping_report.txt"), "VICTIM");
    write_fastq(&dir.join("B/reads.fq"), "MATE");
    let (ok, stderr) = run_in(
        &dir,
        &[
            "--clump_only",
            "--paired",
            "A/reads.fq_clumping_report.txt",
            "B/reads.fq",
        ],
    );
    assert!(
        !ok,
        "a report landing on an input must be refused:\n{stderr}"
    );
    assert_eq!(
        count_reads_from(&dir.join("A/reads.fq_clumping_report.txt"), "VICTIM"),
        40,
        "R1's input was overwritten"
    );
}

/// The trim-mode twin of `clump_paired_shared_mate_report_runs_when_pairs_anchor_apart`,
/// which had no test before #398 and would otherwise have flipped from refuse to
/// succeed unobserved.
#[test]
fn paired_shared_mate_report_runs_when_pairs_anchor_apart() {
    let dir = tempdir("398_trim_mate");
    write_fastq(&dir.join("p/a.fq"), "M1");
    write_fastq(&dir.join("d/x.fq"), "MX");
    write_fastq(&dir.join("q/b.fq"), "M2");
    let (ok, stderr) = run_in(&dir, &["--paired", "p/a.fq", "d/x.fq", "q/b.fq", "d/x.fq"]);
    assert!(ok, "eight distinct paths must not be refused:\n{stderr}");
    assert!(dir.join("p/a.fq_trimming_report.txt").is_file());
    assert!(dir.join("p/x.fq_trimming_report.txt").is_file());
    assert!(dir.join("q/b.fq_trimming_report.txt").is_file());
    assert!(dir.join("q/x.fq_trimming_report.txt").is_file());
    // The shared mate's own directory gains nothing.
    assert_dir_holds_only(&dir.join("d"), &["x.fq"]);
}

/// A bare-filename R1 has an EMPTY parent, so paths stay bare. A `pair_output_dir`
/// that returned `.` would prefix every path here with `./`.
#[test]
fn bare_filename_r1_keeps_paths_unprefixed() {
    let dir = tempdir("398_bare");
    write_fastq(&dir.join("r1.fq"), "Z1");
    write_fastq(&dir.join("r2.fq"), "Z2");
    let (ok, stderr) = run_in(&dir, &["--paired", "r1.fq", "r2.fq"]);
    assert!(ok, "expected success:\n{stderr}");
    assert!(
        !stderr.contains("./r1.fq_trimming_report.txt"),
        "paths must not gain a ./ prefix:\n{stderr}"
    );
    assert_dir_holds_only(
        &dir,
        &[
            "r1.fq",
            "r2.fq",
            "r1_val_1.fq",
            "r2_val_2.fq",
            "r1.fq_trimming_report.txt",
            "r1.fq_trimming_report.json",
            "r2.fq_trimming_report.txt",
            "r2.fq_trimming_report.json",
        ],
    );
}
