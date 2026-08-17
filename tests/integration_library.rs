//! Binary-driven integration tests for `--library` presets
//! (issue [#440](https://github.com/FelixKrueger/TrimGalore/issues/440)).
//!
//! A preset is only worth having if `--library <kit>` and the four clip flags it
//! stands for produce the same bytes, so the load-bearing test here is a
//! differential one against hand-written flags rather than an assertion about
//! internal state. The reporting tests are end-to-end for the same reason the
//! `-a2` suite is: the defect class this feature can hit is "one dispatch path
//! never reads the preset", which no unit test of the resolver can catch.

use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_library_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const R1: &str = "test_files/BS-seq_10K_R1.fastq.gz";
const R2: &str = "test_files/BS-seq_10K_R2.fastq.gz";

fn run(args: &[&str]) -> (bool, String) {
    let out = Command::new(binary())
        .args(args)
        .output()
        .expect("failed to run trim_galore");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

fn read(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("{} unreadable: {e}", path.display()))
}

fn report_text(dir: &Path) -> String {
    String::from_utf8_lossy(&read(
        &dir.join("BS-seq_10K_R1.fastq.gz_trimming_report.txt"),
    ))
    .to_string()
}

fn report_json(dir: &Path) -> serde_json::Value {
    serde_json::from_slice(&read(
        &dir.join("BS-seq_10K_R1.fastq.gz_trimming_report.json"),
    ))
    .expect("JSON report is not valid JSON")
}

/// The contract in one test: a preset is exactly its four clip flags.
#[test]
fn preset_output_is_identical_to_the_flags_it_stands_for() {
    let preset_dir = tempdir("accel_preset");
    let flags_dir = tempdir("accel_flags");

    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--library",
        "accel",
        "-o",
        preset_dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "preset run failed:\n{err}");

    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--clip_R1",
        "10",
        "--clip_R2",
        "15",
        "--three_prime_clip_R1",
        "10",
        "--three_prime_clip_R2",
        "10",
        "-o",
        flags_dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "explicit-flag run failed:\n{err}");

    for name in ["BS-seq_10K_R1_val_1.fq", "BS-seq_10K_R2_val_2.fq"] {
        assert_eq!(
            read(&preset_dir.join(name)),
            read(&flags_dir.join(name)),
            "{name} differs between --library accel and its four clip flags"
        );
    }
}

/// Non-vacuity: the preset must actually change the reads. Without this, the
/// equivalence test above would still pass if presets did nothing at all.
#[test]
fn preset_changes_the_output_relative_to_no_clipping() {
    let with_preset = tempdir("scbs_on");
    let without = tempdir("scbs_off");

    for (dir, extra) in [
        (&with_preset, vec!["--library", "scbs"]),
        (&without, vec![]),
    ] {
        let mut args = vec!["--dont_gzip", "-o", dir.to_str().unwrap()];
        args.extend(extra);
        args.push(R1);
        let (ok, err) = run(&args);
        assert!(ok, "run failed:\n{err}");
    }

    assert_ne!(
        read(&with_preset.join("BS-seq_10K_R1_trimmed.fq")),
        read(&without.join("BS-seq_10K_R1_trimmed.fq")),
        "--library scbs left the reads untouched"
    );
}

/// Explicit wins, and the log names both numbers so the user can see which one
/// was dropped.
#[test]
fn explicit_clip_flag_overrides_the_preset_and_is_logged() {
    let dir = tempdir("override");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--library",
        "emseq",
        "--clip_R1",
        "12",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(
        err.contains("--clip_R1 12") && err.contains("preset value 10"),
        "override log line missing both values:\n{err}"
    );
    // The banner is the summary's other consumer. Only --clip_R1 is overridden here,
    // so no Read 2 flag reaches stderr from either the banner or an override line.
    assert!(
        !err.contains("--clip_R2") && !err.contains("--three_prime_clip_R2"),
        "stderr announces Read 2 clipping on a single-end run:\n{err}"
    );

    // And the bytes follow the override, not the preset.
    let expected = tempdir("override_oracle");
    let (ok, err2) = run(&[
        "--dont_gzip",
        "--clip_R1",
        "12",
        "--clip_R2",
        "10",
        "--three_prime_clip_R1",
        "10",
        "--three_prime_clip_R2",
        "10",
        "-o",
        expected.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "oracle run failed:\n{err2}");
    assert_eq!(
        read(&dir.join("BS-seq_10K_R1_trimmed.fq")),
        read(&expected.join("BS-seq_10K_R1_trimmed.fq"))
    );
}

/// Presets may change between releases, so a report has to be self-describing:
/// preset name plus the four values that were actually used.
#[test]
fn both_reports_record_the_preset_and_its_expanded_values() {
    let dir = tempdir("reporting");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--library",
        "pbat",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");

    let text = report_text(&dir);
    assert!(text.contains("Library preset: pbat"), "got:\n{text}");
    // Single-end: Read 2 clipping is not in force and is not listed.
    assert!(
        text.contains("Clipping in force: --clip_R1 8 --three_prime_clip_R1 8"),
        "got:\n{text}"
    );

    let json = report_json(&dir);
    assert_eq!(json["parameters"]["library"]["preset"], "pbat");
    assert!(
        json["parameters"]["library"]["overrides"]
            .as_array()
            .unwrap()
            .is_empty()
    );

    // The JSON records the resolved configuration, so the block and the top-level
    // keys must agree on every one of the four values.
    for key in [
        "clip_r1",
        "clip_r2",
        "three_prime_clip_r1",
        "three_prime_clip_r2",
    ] {
        assert_eq!(
            json["parameters"]["library"][key], json["parameters"][key],
            "{key} disagrees between the library block and the top level"
        );
        assert_eq!(json["parameters"]["library"][key], 8, "{key}");
    }
}

/// A run without a preset must look exactly as it did before this feature.
#[test]
fn a_run_without_a_preset_reports_no_library_block() {
    let dir = tempdir("no_preset");
    let (ok, err) = run(&["--dont_gzip", "-o", dir.to_str().unwrap(), R1]);
    assert!(ok, "run failed:\n{err}");

    assert!(!report_text(&dir).contains("Library preset"));
    let json = report_json(&dir);
    assert!(json["parameters"]["library"].is_null());
    // `Value::Null` is also what a missing key yields, so pin the key's presence.
    assert!(
        json["parameters"]
            .as_object()
            .expect("parameters is an object")
            .contains_key("library"),
        "the library key must be present and null, not absent"
    );
}

/// An unknown kit name has to stop the run. Trimming with no clipping at all
/// because of a typo is the failure mode this feature must not introduce.
#[test]
fn an_unknown_preset_name_is_refused() {
    let dir = tempdir("unknown");
    let (ok, err) = run(&["--library", "nugen", "-o", dir.to_str().unwrap(), R1]);
    assert!(!ok, "unknown preset was accepted");
    assert!(err.contains("nugen"), "unhelpful error:\n{err}");
}

/// The paired half: all four values apply, so the report and the banner both list
/// all four. Guards against a gate that drops Read 2 unconditionally.
#[test]
fn paired_report_lists_all_four_values() {
    let dir = tempdir("paired_all_four");
    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--library",
        "accel",
        "-o",
        dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "run failed:\n{err}");

    let expected = "--clip_R1 10 --clip_R2 15 --three_prime_clip_R1 10 --three_prime_clip_R2 10";
    let text = report_text(&dir);
    assert!(
        text.contains(&format!("Clipping in force: {expected}")),
        "got:\n{text}"
    );
    assert!(
        err.contains(expected),
        "the banner must list all four on a paired run:\n{err}"
    );
}

/// Four overrides at once — the only shape that drives the JSON `overrides` array
/// past a single element.
#[test]
fn preset_with_multiple_overrides_round_trips_through_json() {
    let dir = tempdir("multi_override");
    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--library",
        "accel",
        "--clip_R1",
        "1",
        "--clip_R2",
        "2",
        "--three_prime_clip_R1",
        "3",
        "--three_prime_clip_R2",
        "4",
        "-o",
        dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "run failed:\n{err}");

    let json = report_json(&dir);
    let overrides = json["parameters"]["library"]["overrides"]
        .as_array()
        .expect("overrides is an array");
    assert_eq!(overrides.len(), 4);
    for (i, (flag, preset_value, user_value)) in [
        ("--clip_R1", 10, 1),
        ("--clip_R2", 15, 2),
        ("--three_prime_clip_R1", 10, 3),
        ("--three_prime_clip_R2", 10, 4),
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(overrides[i]["flag"], *flag, "override {i} flag");
        assert_eq!(
            overrides[i]["preset_value"], *preset_value,
            "override {i} preset_value"
        );
        assert_eq!(
            overrides[i]["user_value"], *user_value,
            "override {i} user_value"
        );
    }
}

/// No surface announces a Read 2 clip that a single-end run did not perform. The
/// warning names the flag on purpose — it says the value was *not* used — so the
/// assertion targets the override sentence rather than the flag name.
#[test]
fn single_end_announces_no_read_2_clip() {
    let dir = tempdir("se_no_r2_announcement");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--library",
        "emseq",
        "--clip_R2",
        "5",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");

    assert!(
        !err.contains("--clip_R2 5 was given on the command line"),
        "stderr announces a Read 2 override on a single-end run:\n{err}"
    );
    assert!(
        err.contains("is not used in this mode"),
        "it should say the value was ignored:\n{err}"
    );

    let text = report_text(&dir);
    assert!(
        !text.contains("--clip_R2"),
        "the report must not mention Read 2 clipping at all:\n{text}"
    );
}

/// The paired half: Read 2 overrides do apply, so both surfaces keep naming them.
#[test]
fn paired_still_announces_r2_overrides() {
    let dir = tempdir("pe_keeps_r2_announcement");
    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--library",
        "emseq",
        "--clip_R2",
        "5",
        "-o",
        dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(
        err.contains("--clip_R2 5 was given on the command line"),
        "the banner must still name a Read 2 override on a paired run:\n{err}"
    );
    let text = report_text(&dir);
    assert!(
        text.contains("--clip_R2 5 was given on the command line"),
        "the report must still name it:\n{text}"
    );
}
