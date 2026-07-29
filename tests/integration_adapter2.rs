//! Binary-driven integration tests for `-a2`/`--adapter2` precedence
//! (issue [#369](https://github.com/FelixKrueger/TrimGalore/issues/369)).
//!
//! Before #369, adapter resolution was an early-return chain in which only the
//! explicit-`-a` branch read `cli.adapter2`. Every preset and the default
//! auto-detect path returned a hardcoded Read 2 adapter, so `-a2` was silently
//! discarded unless `-a` was also given. The contract restored here is Perl
//! 0.6.11's `unless (defined $a2)`: a user `-a2` always wins, and presets or
//! auto-detection supply a Read 2 default only in its absence.
//!
//! These are end-to-end and per-branch on purpose. A unit test of the override
//! helper cannot catch this bug class — the defect was never "the helper computes
//! the wrong answer", it was "six branches never called it". Only running each
//! branch proves the call happens.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The reporter's Read 2 adapter. Chosen because it discriminates on the
/// fixture: it moves R2 `Reads with adapters` from 48.8% to 54.4%, so a test
/// asserting output changed cannot pass vacuously.
const A2: &str = "AAATCAAAAAAAC";

/// The sequence `--illumina` selects, and what auto-detection picks on this
/// fixture. Used to build the differential oracle in [`oracle_equivalence`].
const ILLUMINA: &str = "AGATCGGAAGAGC";

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_a2_{tag}_{}", std::process::id()));
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

/// A named line from the R2 trimming report.
fn r2_report_line(dir: &Path, needle: &str) -> String {
    let report = std::fs::read_dir(dir)
        .expect("output dir missing")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| {
            let n = p.file_name().unwrap_or_default().to_string_lossy();
            n.contains("R2") && n.ends_with("trimming_report.txt")
        })
        .expect("no R2 trimming report");
    std::fs::read_to_string(&report)
        .expect("unreadable report")
        .lines()
        .find(|l| l.contains(needle))
        .unwrap_or_else(|| panic!("no '{needle}' line in R2 report"))
        .to_string()
}

/// The adapter the R2 trimming report says was handed to the trimmer.
fn r2_adapter(dir: &Path) -> String {
    let report = std::fs::read_dir(dir)
        .expect("output dir missing")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| {
            let n = p.file_name().unwrap_or_default().to_string_lossy();
            n.contains("R2") && n.ends_with("trimming_report.txt")
        })
        .expect("no R2 trimming report");
    let body = std::fs::read_to_string(&report).expect("unreadable report");
    let line = body
        .lines()
        .find(|l| l.contains("Command line parameters"))
        .expect("no command line in report");
    line.split_whitespace()
        .skip_while(|t| *t != "-a")
        .nth(1)
        .unwrap_or("")
        .to_string()
}

/// Decompressed contents, compared directly rather than hashed — no hashing
/// dependency needed, and a mismatch can be inspected.
fn gunzip(path: &Path) -> Vec<u8> {
    use std::io::Read;
    let f = std::fs::File::open(path).expect("missing output file");
    let mut d = flate2::read::MultiGzDecoder::new(f);
    let mut buf = Vec::new();
    d.read_to_end(&mut buf).expect("bad gzip");
    buf
}

fn val_outputs(dir: &Path) -> (Vec<u8>, Vec<u8>) {
    let mut v1 = None;
    let mut v2 = None;
    for e in std::fs::read_dir(dir)
        .expect("output dir missing")
        .flatten()
    {
        let p = e.path();
        let n = p
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        if n.ends_with("_val_1.fq.gz") {
            v1 = Some(gunzip(&p));
        } else if n.ends_with("_val_2.fq.gz") {
            v2 = Some(gunzip(&p));
        }
    }
    (v1.expect("no _val_1 output"), v2.expect("no _val_2 output"))
}

// ── -a2 wins on every resolution path ──────────────────────────────────

/// The core contract, one case per branch of adapter resolution. `--illumina`
/// is the reported case; the auto-detect row (no adapter flag at all) is the
/// one most likely to affect other users, since it is the default path.
#[test]
fn a2_wins_on_every_resolution_path() {
    let cases: [(&str, Vec<&str>); 7] = [
        ("explicit_a", vec!["-a", ILLUMINA]),
        ("illumina", vec!["--illumina"]),
        ("nextera", vec!["--nextera"]),
        ("stranded", vec!["--stranded_illumina"]),
        ("small_rna", vec!["--small_rna"]),
        ("bgiseq", vec!["--bgiseq"]),
        ("autodetect", vec![]),
    ];
    for (tag, flags) in cases {
        let dir = tempdir(&format!("win_{tag}"));
        let mut args = vec!["--paired"];
        args.extend(flags.iter().copied());
        args.extend(["-a2", A2, R1, R2, "-o"].iter().copied());
        let dir_s = dir.to_str().unwrap().to_string();
        let args: Vec<&str> = args
            .into_iter()
            .chain(std::iter::once(dir_s.as_str()))
            .collect();
        let (ok, stderr) = run(&args);
        assert!(ok, "[{tag}] run failed: {stderr}");
        assert_eq!(
            r2_adapter(&dir),
            A2,
            "[{tag}] Read 2 must be trimmed with the -a2 sequence, not a preset default"
        );
    }
}

/// The differential oracle. Once `-a2` is honoured, `--illumina -a2 SEQ` and
/// `-a AGATCGGAAGAGC -a2 SEQ` are the same invocation in everything that reaches
/// the FASTQ output — same R1 sequence, same R2 sequence, same length cutoff
/// (keyed on the R1 sequence), same poly-G decision. So their output must match
/// byte for byte, on both reads. Stronger than reading the report line, and it
/// needs no Perl installation.
#[test]
fn oracle_equivalence() {
    let explicit = tempdir("oracle_explicit");
    let (ok, err) = run(&[
        "--paired",
        "-a",
        ILLUMINA,
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        explicit.to_str().unwrap(),
    ]);
    assert!(ok, "explicit run failed: {err}");

    for (tag, flag) in [("illumina", Some("--illumina")), ("autodetect", None)] {
        let dir = tempdir(&format!("oracle_{tag}"));
        let mut args = vec!["--paired"];
        if let Some(f) = flag {
            args.push(f);
        }
        let dir_s = dir.to_str().unwrap().to_string();
        args.extend(["-a2", A2, R1, R2, "-o"].iter().copied());
        let args: Vec<&str> = args
            .into_iter()
            .chain(std::iter::once(dir_s.as_str()))
            .collect();
        let (ok, err) = run(&args);
        assert!(ok, "[{tag}] run failed: {err}");
        let (got_r1, got_r2) = val_outputs(&dir);
        let (want_r1, want_r2) = val_outputs(&explicit);
        assert!(
            got_r1 == want_r1 && got_r2 == want_r2,
            "[{tag}] output must be identical to the explicit -a/-a2 equivalent, both reads"
        );
    }
}

/// Read 1's output file changes for these invocations — by design, via the joint
/// pair-length filter — so this pins the property that actually holds: Read 1's
/// own per-read trimming is untouched.
#[test]
fn read1_trimming_is_unaffected_by_the_r2_adapter() {
    let with = tempdir("r1_with_a2");
    let without = tempdir("r1_without_a2");
    run(&[
        "--paired",
        "--illumina",
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        with.to_str().unwrap(),
    ]);
    run(&[
        "--paired",
        "--illumina",
        R1,
        R2,
        "-o",
        without.to_str().unwrap(),
    ]);

    let stat = |dir: &Path| -> String {
        let report = std::fs::read_dir(dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| {
                let n = p.file_name().unwrap_or_default().to_string_lossy();
                n.contains("R1") && n.ends_with("trimming_report.txt")
            })
            .expect("no R1 report");
        std::fs::read_to_string(report)
            .unwrap()
            .lines()
            .find(|l| l.contains("Reads written (passing filters)"))
            .unwrap()
            .to_string()
    };
    assert_eq!(
        stat(&with),
        stat(&without),
        "Read 1's own trimming must not depend on the Read 2 adapter"
    );
}

// ── absence of -a2 leaves defaults intact ──────────────────────────────

#[test]
fn preset_r2_defaults_survive_when_no_a2_given() {
    for (tag, flag, expected) in [
        ("small_rna", "--small_rna", "GATCGTCGGACT"),
        (
            "bgiseq",
            "--bgiseq",
            "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
        ),
    ] {
        let dir = tempdir(&format!("default_{tag}"));
        let (ok, err) = run(&["--paired", flag, R1, R2, "-o", dir.to_str().unwrap()]);
        assert!(ok, "[{tag}] run failed: {err}");
        assert_eq!(
            r2_adapter(&dir),
            expected,
            "[{tag}] the preset R2 default must still apply when -a2 is absent"
        );
    }
}

/// `--illumina` has no Read 2 default, so Read 2 falls back to Read 1's adapter.
#[test]
fn r2_falls_back_to_r1_when_no_default_and_no_a2() {
    let dir = tempdir("fallback");
    run(&[
        "--paired",
        "--illumina",
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert_eq!(r2_adapter(&dir), ILLUMINA);
}

// ── announcement, and the cases where -a2 cannot apply ─────────────────

/// The NOTE fires only when a real default was displaced. `--illumina` has none,
/// so displacing nothing must not be announced.
#[test]
fn note_fires_only_on_a_real_displacement() {
    let dir = tempdir("note_small_rna");
    let (_, stderr) = run(&[
        "--paired",
        "--small_rna",
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(
        stderr.contains("NOTE: Read 2 adapter taken from -a2") && stderr.contains("GATCGTCGGACT"),
        "displacing the small RNA default should name it; got: {stderr}"
    );

    let dir = tempdir("note_illumina");
    let (_, stderr) = run(&[
        "--paired",
        "--illumina",
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(
        !stderr.contains("NOTE: Read 2 adapter taken from -a2"),
        "--illumina has no R2 default, so nothing is displaced; got: {stderr}"
    );
}

/// `--consider_already_trimmed` suppression wins over `-a2`: trimming R2 while
/// R1 is left alone would be asymmetric under a banner promising quality
/// trimming only.
#[test]
fn suppression_wins_over_a2() {
    let dir = tempdir("suppressed");
    let (ok, stderr) = run(&[
        "--paired",
        "--consider_already_trimmed",
        "10000",
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(ok, "run should succeed: {stderr}");
    assert!(
        stderr.contains("adapter trimming is suppressed"),
        "must say -a2 was not applied; got: {stderr}"
    );
    // Asserted on behaviour, not on the echoed command line: the suppressed
    // adapter is empty, which `split_whitespace` collapses away.
    assert!(
        r2_report_line(&dir, "Reads with adapters").contains(" 0 (0.0%)"),
        "Read 2 must stay untrimmed, symmetric with Read 1"
    );
}

/// Modes that perform no adapter trimming must say `-a2` is unused, with a
/// reason true for that mode — and must still succeed.
#[test]
fn unusable_a2_warns_with_a_mode_accurate_reason() {
    // `--implicon` uses require_equals, hence `--implicon=8`.
    let cases: [(&str, Vec<&str>, &str); 5] = [
        ("clock", vec!["--paired", "--clock"], "--clock"),
        ("implicon", vec!["--paired", "--implicon=8"], "--implicon"),
        (
            "hardtrim5",
            vec!["--paired", "--hardtrim5", "30"],
            "--hardtrim5/--hardtrim3",
        ),
        (
            "hardtrim3",
            vec!["--paired", "--hardtrim3", "5"],
            "--hardtrim5/--hardtrim3",
        ),
        ("single_end", vec!["--illumina"], "single-end"),
    ];
    for (tag, flags, reason) in cases {
        let dir = tempdir(&format!("unusable_{tag}"));
        let mut args = flags.clone();
        args.extend(["-a2", A2].iter().copied());
        if tag == "single_end" {
            args.push(R1);
        } else if tag == "clock" || tag == "implicon" {
            args.push("test_files/clock_10K_R1.fastq.gz");
            args.push("test_files/clock_10K_R2.fastq.gz");
        } else {
            args.push(R1);
            args.push(R2);
        }
        let dir_s = dir.to_str().unwrap().to_string();
        args.push("-o");
        let args: Vec<&str> = args
            .into_iter()
            .chain(std::iter::once(dir_s.as_str()))
            .collect();
        let (ok, stderr) = run(&args);
        assert!(ok, "[{tag}] should warn, not fail: {stderr}");
        assert!(
            stderr.contains("not used in this mode") && stderr.contains(reason),
            "[{tag}] warning must give a reason true for this mode; got: {stderr}"
        );
    }
}

/// Single-end must not announce a Read 2 adapter it is about to ignore, and must
/// not contradict its own warning. The FASTQ is byte-identical either way, so only
/// stderr can catch this.
///
/// `--small_rna` is the load-bearing case: it has a Read 2 default, so it is the
/// only shape that can produce the displacement NOTE. An earlier version of this
/// fix gated the `Adapter 2 (Read 2):` line but not the NOTE, so a single-end run
/// warned that `-a2` was unused and then announced that `-a2` had displaced the
/// preset default. `--illumina` cannot catch that — it has no default to displace.
#[test]
fn single_end_does_not_announce_a_read2_adapter() {
    for (tag, preset) in [("illumina", "--illumina"), ("small_rna", "--small_rna")] {
        let dir = tempdir(&format!("se_no_r2_line_{tag}"));
        let (ok, stderr) = run(&[preset, "-a2", A2, R1, "-o", dir.to_str().unwrap()]);
        assert!(ok, "[{tag}] run failed: {stderr}");
        assert!(
            !stderr.contains("Adapter 2 (Read 2)"),
            "[{tag}] single-end must not advertise a Read 2 adapter; got: {stderr}"
        );
        assert!(
            !stderr.contains("NOTE: Read 2 adapter taken from -a2"),
            "[{tag}] single-end must not claim -a2 displaced anything; got: {stderr}"
        );
        assert!(
            stderr.contains("not used in this mode"),
            "[{tag}] single-end should still warn that -a2 is unused; got: {stderr}"
        );
    }
}

// ── newly-validated input ──────────────────────────────────────────────

/// `-a2` is now parsed on every path, so malformed values are rejected instead
/// of silently discarded. This is a deliberate behaviour change.
#[test]
fn malformed_a2_is_rejected_on_preset_paths() {
    // Asserting the message, not just the exit code: a bare `!ok` would also pass
    // if the run failed for an unrelated reason.
    for (bad, expected) in [
        ("ZZZQQQ", "must contain only DNA characters"),
        ("", "Empty adapter sequence"),
        ("file:definitely_missing.fa", "adapter FASTA"),
    ] {
        let dir = tempdir("malformed");
        let (ok, stderr) = run(&[
            "--paired",
            "--illumina",
            "-a2",
            bad,
            R1,
            R2,
            "-o",
            dir.to_str().unwrap(),
        ]);
        assert!(
            !ok,
            "-a2 '{bad}' must be rejected, not silently ignored; stderr: {stderr}"
        );
        assert!(
            stderr.contains(expected),
            "-a2 '{bad}' should fail with '{expected}'; got: {stderr}"
        );
    }
}

/// A value the run has just announced it will ignore must not fail the run. These
/// invocations completed before #369 and must keep completing.
#[test]
fn malformed_a2_is_tolerated_where_the_mode_ignores_it() {
    let cases: [(&str, Vec<&str>); 2] = [
        ("clock", vec!["--paired", "--clock"]),
        ("single_end", vec!["--illumina"]),
    ];
    for (tag, flags) in cases {
        let dir = tempdir(&format!("tolerated_{tag}"));
        let mut args = flags.clone();
        args.extend(["-a2", "ZZZQQQ"].iter().copied());
        if tag == "single_end" {
            args.push(R1);
        } else {
            args.push("test_files/clock_10K_R1.fastq.gz");
            args.push("test_files/clock_10K_R2.fastq.gz");
        }
        let dir_s = dir.to_str().unwrap().to_string();
        args.push("-o");
        let args: Vec<&str> = args
            .into_iter()
            .chain(std::iter::once(dir_s.as_str()))
            .collect();
        let (ok, stderr) = run(&args);
        assert!(
            ok,
            "[{tag}] a malformed -a2 that the mode ignores must not fail the run; got: {stderr}"
        );
        assert!(
            stderr.contains("not used in this mode"),
            "[{tag}] should still warn; got: {stderr}"
        );
    }
}

/// Repeated `-a2` reaches Read 2 on a preset path, and sizes the per-adapter
/// stats correctly — a shape previously reachable only via `-a`.
#[test]
fn repeated_a2_applies_on_a_preset_path() {
    let dir = tempdir("repeated");
    let (ok, stderr) = run(&[
        "--paired",
        "--illumina",
        "-a2",
        A2,
        "-a2",
        "GGGGCCCCTTTT",
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(ok, "run failed: {stderr}");
    assert!(
        stderr.contains("Adapters R2 (2 sequences)"),
        "both -a2 sequences should reach Read 2; got: {stderr}"
    );
    // The stats vector is sized from adapters_r2, so the report must break down two.
    assert!(
        r2_report_line(&dir, "GGGGCCCCTTTT").contains("GGGGCCCCTTTT"),
        "the R2 report should carry a per-adapter entry for the second sequence"
    );
}

/// `A{N}` expansion works for `-a2` on a preset path. Documented as supported
/// since v2, but true only on the `-a` branch before this fix.
#[test]
fn brace_expansion_works_for_a2_on_a_preset_path() {
    let dir = tempdir("brace");
    let (ok, stderr) = run(&[
        "--paired",
        "--illumina",
        "-a2",
        "A{10}",
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(ok, "run failed: {stderr}");
    assert_eq!(r2_adapter(&dir), "AAAAAAAAAA");
}

/// The poly-G scan piggybacks on auto-detection and must not be duplicated or
/// lost. Invisible to output md5 on this fixture (2/10000 is below threshold),
/// so the scan-line count is the only signal.
#[test]
fn polyg_scan_is_not_duplicated_or_lost() {
    let dir = tempdir("polyg_preset");
    let (_, preset) = run(&[
        "--paired",
        "--illumina",
        "-a2",
        A2,
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    let dir = tempdir("polyg_auto");
    let (_, auto) = run(&["--paired", "-a2", A2, R1, R2, "-o", dir.to_str().unwrap()]);
    let count = |s: &str| s.matches("Scanning for poly-G content").count();
    assert_eq!(
        count(&preset),
        1,
        "a preset run scans for poly-G separately"
    );
    assert_eq!(
        count(&auto),
        0,
        "auto-detection piggybacks the poly-G scan; a second pass would be wasted I/O"
    );
}

/// uBAM output resolves adapters through the same `setup_trimming`, so it inherits
/// the fix. Paired uBAM output is one interleaved BAM, so the check is on the report.
#[test]
fn ubam_output_inherits_the_a2_fix() {
    let dir = tempdir("ubam_out");
    let (ok, stderr) = run(&[
        "--paired",
        "--illumina",
        "-a2",
        A2,
        "--output-format",
        "ubam",
        R1,
        R2,
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(ok, "run failed: {stderr}");
    assert_eq!(r2_adapter(&dir), A2);
}
