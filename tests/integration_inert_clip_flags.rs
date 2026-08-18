//! Binary-driven integration tests for clip flags a mode cannot honour
//! (issue [#440](https://github.com/FelixKrueger/TrimGalore/issues/440)).
//!
//! Two inert sets, and they differ: the specialty modes read no clip flag at all,
//! so all four are inert there, while a single-end run honours the Read 1 pair and
//! only the Read 2 pair is inert. A single `!paired` test cannot express that — it
//! reports "this is a single-end run" for a two-file `--clock` run, which sends the
//! user to add `--paired`, silencing the warning without making the flag work.
//!
//! Modelled on `integration_adapter2.rs`, which tests the same warn-and-ignore
//! contract for `-a2` and whose mode table these rows mirror. The specialty rows
//! deliberately pass two files and no `--paired`, since that is the shape a bare
//! `!paired` predicate gets wrong.

use std::path::PathBuf;
use std::process::Command;

const R1: &str = "test_files/BS-seq_10K_R1.fastq.gz";
const R2: &str = "test_files/BS-seq_10K_R2.fastq.gz";
const CLOCK_R1: &str = "test_files/clock_10K_R1.fastq.gz";
const CLOCK_R2: &str = "test_files/clock_10K_R2.fastq.gz";

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_inert_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

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

const INERT: &str = "is not used in this mode";

/// One row of the mode table: a mode, the clip flags it cannot honour, and the
/// reason its warning must give.
struct Case {
    tag: &'static str,
    mode_flags: &'static [&'static str],
    clip_flags: &'static [&'static str],
    reason: &'static str,
    clock_fixture: bool,
}

/// The reason has to be true for the mode, not merely true of some mode. Under the
/// specialty modes it is `--clip_R1` that proves the point: it applies perfectly
/// well on a single-end run, so a warning about it can only come from the mode.
#[test]
fn inert_clip_flag_warns_with_a_mode_accurate_reason() {
    const R1_PAIR: &[&str] = &["--clip_R1", "--three_prime_clip_R1"];
    const R2_PAIR: &[&str] = &["--clip_R2", "--three_prime_clip_R2"];
    const SPECIALTY: &str = "--hardtrim5/--hardtrim3 do not honour the clip flags";

    // Each row carries both flags of its inert set, so all four entries of the
    // warning's array are exercised rather than one representative.
    let cases = [
        Case {
            tag: "hardtrim5",
            mode_flags: &["--hardtrim5", "30"],
            clip_flags: R1_PAIR,
            reason: SPECIALTY,
            clock_fixture: false,
        },
        Case {
            tag: "hardtrim3",
            mode_flags: &["--hardtrim3", "30"],
            clip_flags: R1_PAIR,
            reason: SPECIALTY,
            clock_fixture: false,
        },
        Case {
            tag: "clock",
            mode_flags: &["--clock"],
            clip_flags: R1_PAIR,
            reason: "--clock does not honour the clip flags",
            clock_fixture: true,
        },
        Case {
            tag: "implicon",
            mode_flags: &["--implicon"],
            clip_flags: R1_PAIR,
            reason: "--implicon does not honour the clip flags",
            clock_fixture: true,
        },
        Case {
            tag: "single_end",
            mode_flags: &[],
            clip_flags: R2_PAIR,
            reason: "it applies to Read 2 of a pair, and this is a single-end run",
            clock_fixture: false,
        },
    ];

    for Case {
        tag,
        mode_flags,
        clip_flags,
        reason,
        clock_fixture,
    } in cases
    {
        let dir = tempdir(tag);
        let dir_s = dir.to_str().unwrap().to_string();
        let mut args: Vec<&str> = vec!["--dont_gzip"];
        args.extend(mode_flags.iter().copied());
        for f in clip_flags {
            args.extend([*f, "3"]);
        }
        args.extend(["-o", dir_s.as_str()]);
        if tag == "single_end" {
            args.push(R1);
        } else if clock_fixture {
            args.push(CLOCK_R1);
            args.push(CLOCK_R2);
        } else {
            args.push(R1);
            args.push(R2);
        }

        let (ok, err) = run(&args);
        assert!(ok, "[{tag}] should warn, not fail:\n{err}");
        assert!(
            err.contains(INERT),
            "[{tag}] expected a warning, got:\n{err}"
        );
        assert!(
            err.contains(reason),
            "[{tag}] warning must give a reason true for this mode; got:\n{err}"
        );
        for f in clip_flags {
            assert!(
                err.contains(*f),
                "[{tag}] the warning must name {f}; got:\n{err}"
            );
        }
        assert_eq!(
            err.lines().filter(|l| l.contains(INERT)).count(),
            clip_flags.len(),
            "[{tag}] expected one line per inert flag; got:\n{err}"
        );
    }
}

/// Read 1 clipping applies on a single-end run, so naming it would be false in both
/// halves — it is not a Read 2 flag and it is not being ignored.
#[test]
fn read_1_clip_flags_do_not_warn_on_single_end() {
    let dir = tempdir("se_r1_quiet");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--clip_R1",
        "3",
        "--three_prime_clip_R1",
        "4",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(
        !err.contains(INERT),
        "Read 1 clipping applies on a single-end run and must not be reported inert:\n{err}"
    );
}

/// Read 2 clipping applies on a paired run, so nothing is inert and nothing warns.
#[test]
fn paired_clip_flags_do_not_warn() {
    let dir = tempdir("paired_quiet");
    let (ok, err) = run(&[
        "--paired",
        "--dont_gzip",
        "--clip_R2",
        "5",
        "--three_prime_clip_R2",
        "7",
        "-o",
        dir.to_str().unwrap(),
        R1,
        R2,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(!err.contains(INERT), "paired must not warn:\n{err}");
}

/// A preset supplies all four values, so reading the resolved view instead of the
/// user's own flags would warn on every single-end preset run.
#[test]
fn preset_alone_on_single_end_does_not_warn() {
    let dir = tempdir("preset_quiet");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--library",
        "emseq",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(
        !err.contains(INERT),
        "a preset the user typed no clip flag for must not warn:\n{err}"
    );
}

/// One line per flag, so the sentence stays singular and each flag is named.
#[test]
fn each_inert_flag_gets_its_own_line() {
    let dir = tempdir("two_lines");
    let (ok, err) = run(&[
        "--dont_gzip",
        "--clip_R2",
        "5",
        "--three_prime_clip_R2",
        "7",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert_eq!(
        err.lines().filter(|l| l.contains(INERT)).count(),
        2,
        "expected one warning line per inert flag, got:\n{err}"
    );
    assert!(err.contains("--clip_R2 5"), "got:\n{err}");
    assert!(err.contains("--three_prime_clip_R2 7"), "got:\n{err}");
    // 5' before 3', so the pair reads in clipping order.
    assert!(
        err.find("--clip_R2 5") < err.find("--three_prime_clip_R2 7"),
        "expected --clip_R2 before --three_prime_clip_R2:\n{err}"
    );
}

/// All four at once, under a mode where all four are inert. The four-flag order
/// is otherwise unexercised — every other case sets at most two.
#[test]
fn all_four_inert_flags_warn_in_clipping_order() {
    let dir = tempdir("four_lines");
    let (ok, err) = run(&[
        "--hardtrim5",
        "20",
        "--clip_R1",
        "3",
        "--three_prime_clip_R1",
        "4",
        "--clip_R2",
        "5",
        "--three_prime_clip_R2",
        "7",
        "-o",
        dir.to_str().unwrap(),
        R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert_eq!(
        err.lines().filter(|l| l.contains(INERT)).count(),
        4,
        "expected one warning line per inert flag, got:\n{err}"
    );
    let positions: Vec<Option<usize>> = [
        "--clip_R1 3",
        "--three_prime_clip_R1 4",
        "--clip_R2 5",
        "--three_prime_clip_R2 7",
    ]
    .iter()
    .map(|f| err.find(f))
    .collect();
    assert!(
        positions.iter().all(|p| p.is_some()),
        "every flag must be named:\n{err}"
    );
    assert!(
        positions.windows(2).all(|w| w[0] < w[1]),
        "expected Read 1 5'/3' then Read 2 5'/3':\n{err}"
    );
}
