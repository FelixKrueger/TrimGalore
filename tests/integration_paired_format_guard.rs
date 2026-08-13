//! Binary-driven integration tests for the `--paired` per-pair format guard
//! (issue [#363](https://github.com/FelixKrueger/TrimGalore/issues/363)).
//!
//! Before #363, a pair containing exactly *one* BAM was reported as "two BAM
//! files" and offered the single-interleaved-file remediation — advice that does
//! not apply when the real cause is a mis-typed filename. The guard now
//! distinguishes the shapes, and runs in the cheap-validation window in `main()`
//! rather than inside the per-pair worker.
//!
//! Four properties are pinned here that unit tests cannot reach:
//!
//! 1. The rejection has **no side effects** — no output directory, no adapter
//!    auto-detection, and on multi-pair input no partial output from earlier
//!    pairs. All three happened before #363.
//! 2. `--clump_only` on the FASTQ-output path keeps its own aux-tag diagnosis.
//!    This path had zero coverage and is where the first draft of the fix
//!    regressed: it would have printed `--clump_only --paired interleaved.bam`,
//!    a command the binary itself rejects.
//! 3. A plain+gzip FASTQ pair is still **accepted** — the guard keys on BAM
//!    count, not format equality.
//! 4. `--fastqc` is refused on this same `--paired`-with-one-input family, whose
//!    FASTQ-output path writes no QC report (#421) — and accepted one flag away,
//!    under `--output-format ubam`, which does.
//!
//! Fixtures: `phred64_test.fastq` (plain FASTQ), `BS-seq_10K_R{1,2}.fastq.gz`
//! (gzipped FASTQ), `ubam_test.bam` / `ubam_paired_test.bam` (two distinct
//! uBAMs — distinct matters, since `Cli::validate` rejects R1 == R2 first).

use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// A pre-created scratch directory. Note that tests asserting the *absence* of
/// an output directory must pass a `nested` child of this (see
/// [`nonexistent_out`]) — handing this path to `-o` would make the check
/// vacuous, because it already exists.
fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_pfg_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A path under `dir` that does not exist yet, so `!exists()` after a rejected
/// run is a real assertion about `ensure_output_dir` never having run.
fn nonexistent_out(dir: &Path) -> PathBuf {
    let p = dir.join("nested_out");
    assert!(!p.exists(), "precondition: {} must not exist", p.display());
    p
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from("test_files").join(name)
}

/// Run the binary and return (success, stderr).
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

// ── mixed FASTQ + BAM within one pair ──────────────────────────────────

/// The reported bug: one BAM in the pair must not be described as two, and must
/// not be offered the interleaving remediation.
#[test]
fn mixed_pair_reports_format_mismatch_not_two_bams() {
    for (r1, r2) in [
        (fixture("phred64_test.fastq"), fixture("ubam_test.bam")),
        (fixture("ubam_test.bam"), fixture("phred64_test.fastq")),
    ] {
        let dir = tempdir("mixed");
        let (ok, stderr) = run(&[
            "--paired",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "-o",
            nonexistent_out(&dir).to_str().unwrap(),
        ]);
        assert!(!ok, "a mixed FASTQ+BAM pair must be rejected");
        assert!(
            stderr.contains("same format") && stderr.contains("mixed"),
            "expected a format-mismatch message; got: {stderr}"
        );
        // #363's substance. Without this assertion the two messages can
        // silently re-converge in a later refactor.
        assert!(
            !stderr.contains("uBAM paired mode expects"),
            "mixed pair must not be told to pass a single interleaved file; got: {stderr}"
        );
        assert!(
            !stderr.contains("two BAM files"),
            "mixed pair must not be described as two BAM files; got: {stderr}"
        );
        // Both inputs named, with human-readable format labels.
        assert!(
            stderr.contains("phred64_test.fastq") && stderr.contains("ubam_test.bam"),
            "both inputs should be named; got: {stderr}"
        );
        assert!(stderr.contains("uBAM"), "got: {stderr}");
        assert!(
            !stderr.contains("UnalignedBam") && !stderr.contains("FastqPlain"),
            "internal enum names must not leak into user-facing text; got: {stderr}"
        );
    }
}

/// The rejection must happen before anything is created or scanned. Before
/// #363 this run auto-detected adapters over the inputs, printed the trimming
/// banner, and created the output directory before failing.
#[test]
fn mixed_pair_rejection_has_no_side_effects() {
    let dir = tempdir("noside");
    let out = nonexistent_out(&dir);
    let (ok, stderr) = run(&[
        "--paired",
        fixture("phred64_test.fastq").to_str().unwrap(),
        fixture("ubam_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(
        !out.exists(),
        "a rejected run must not create the output directory: {}",
        out.display()
    );
    assert!(
        !stderr.contains("Auto-detecting adapter"),
        "adapter auto-detection must not run before a guaranteed rejection; got: {stderr}"
    );
    assert!(
        !stderr.contains("Trimming (paired-end)"),
        "the trimming banner must not print before a guaranteed rejection; got: {stderr}"
    );
}

// ── two BAM files within one pair ──────────────────────────────────────

/// The genuine two-BAM case keeps its message and its remediation. This also
/// closes a coverage gap: the plain trim path (no `--output-format ubam`) had
/// no test at all.
#[test]
fn two_bam_pair_keeps_interleave_remediation() {
    let dir = tempdir("twobam");
    let (ok, stderr) = run(&[
        "--paired",
        fixture("ubam_test.bam").to_str().unwrap(),
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        nonexistent_out(&dir).to_str().unwrap(),
    ]);
    assert!(
        !ok,
        "two distinct BAM files under --paired must be rejected"
    );
    assert!(
        stderr.contains("two BAM files is not supported"),
        "got: {stderr}"
    );
    assert!(
        stderr.contains("uBAM paired mode expects") && stderr.contains("single interleaved"),
        "got: {stderr}"
    );
    // `samtools collate` cannot combine two files — its second positional is
    // the temp-file prefix, so it exits 0 having emitted only the first file's
    // records. The suggested command must be one that actually works.
    assert!(
        stderr.contains("samtools merge -n"),
        "expected a verified combine command; got: {stderr}"
    );
    // Asserts the property, not one literal spelling of its violation: an
    // earlier draft pinned the full `collate -O <fixture> <fixture>` string,
    // which a fixture rename would have quietly turned into a tautology.
    assert!(
        !stderr.contains("samtools collate"),
        "collate takes one input — it must never be suggested for combining two \
         files on this path; got: {stderr}"
    );
    // Both files named, not just one.
    assert!(
        stderr.contains("ubam_test.bam") && stderr.contains("ubam_paired_test.bam"),
        "got: {stderr}"
    );
}

/// Same rejection on the uBAM-output path, with the remediation echoing the
/// user's own mode so it can be pasted.
#[test]
fn two_bam_pair_rejected_with_ubam_output_and_echoes_mode() {
    let dir = tempdir("twobam_ubamout");
    let (ok, stderr) = run(&[
        "--paired",
        "--output-format",
        "ubam",
        fixture("ubam_test.bam").to_str().unwrap(),
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        nonexistent_out(&dir).to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(
        stderr.contains("trim_galore --paired --output-format ubam interleaved.bam"),
        "remediation should reproduce the user's mode; got: {stderr}"
    );
}

#[test]
fn mixed_pair_rejected_with_ubam_output() {
    let dir = tempdir("mixed_ubamout");
    let (ok, stderr) = run(&[
        "--paired",
        "--output-format",
        "ubam",
        fixture("phred64_test.fastq").to_str().unwrap(),
        fixture("ubam_test.bam").to_str().unwrap(),
        "-o",
        nonexistent_out(&dir).to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(
        stderr.contains("same format") && !stderr.contains("uBAM paired mode expects"),
        "got: {stderr}"
    );
}

// ── --clump_only, FASTQ output ─────────────────────────────────────────

/// LOAD-BEARING. This path had zero coverage, and the first draft of the #363
/// fix regressed it: the hoisted guard preempted the aux-tag diagnosis and
/// printed `trim_galore --clump_only --paired interleaved.bam` instead — a
/// command `main()` rejects outright, so the user would paste it and get a
/// second error, having lost the one piece of information they needed.
///
/// In this mode the predicate is "any BAM in the pair", not "the pair
/// disagrees": both a mixed pair and two BAMs are the same user error, namely
/// BAM input on a FASTQ-output path.
#[test]
fn clump_only_fastq_output_keeps_aux_tag_diagnosis() {
    let cases: [(PathBuf, PathBuf); 3] = [
        (fixture("BS-seq_10K_R1.fastq.gz"), fixture("ubam_test.bam")),
        (fixture("ubam_test.bam"), fixture("BS-seq_10K_R1.fastq.gz")),
        (fixture("ubam_test.bam"), fixture("ubam_paired_test.bam")),
    ];
    for (r1, r2) in cases {
        let dir = tempdir("clump_fq");
        let (ok, stderr) = run(&[
            "--clump_only",
            "--paired",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "-o",
            nonexistent_out(&dir).to_str().unwrap(),
        ]);
        assert!(
            !ok,
            "uBAM input on the clump-only FASTQ path must be rejected"
        );
        assert!(
            stderr.contains("--output-format ubam"),
            "must name the flag that fixes it; got: {stderr}"
        );
        assert!(
            stderr.contains("drop aux tags"),
            "must give the reason, not just the flag; got: {stderr}"
        );
        assert!(
            !stderr.contains("interleaved.bam"),
            "must not suggest an interleaved uBAM — `--clump_only --paired \
             interleaved.bam` is itself rejected; got: {stderr}"
        );
    }
}

/// The clump-only uBAM-output path had the only correct implementation of this
/// distinction before #363; confirm it survived being replaced by the shared
/// helper, including the mode-specific remediation.
#[test]
fn clump_only_ubam_output_distinguishes_mixed_from_two_bam() {
    let dir = tempdir("clump_ubam_mixed");
    let (ok, stderr) = run(&[
        "--clump_only",
        "--paired",
        "--output-format",
        "ubam",
        fixture("BS-seq_10K_R1.fastq.gz").to_str().unwrap(),
        fixture("ubam_test.bam").to_str().unwrap(),
        "-o",
        nonexistent_out(&dir).to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(stderr.contains("same format"), "got: {stderr}");
    assert!(
        !stderr.contains("uBAM paired mode expects"),
        "got: {stderr}"
    );

    let dir = tempdir("clump_ubam_two");
    let (ok, stderr) = run(&[
        "--clump_only",
        "--paired",
        "--output-format",
        "ubam",
        fixture("ubam_test.bam").to_str().unwrap(),
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        nonexistent_out(&dir).to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(
        stderr.contains("two BAM files is not supported"),
        "got: {stderr}"
    );
    assert!(
        stderr.contains("trim_galore --clump_only --paired --output-format ubam interleaved.bam"),
        "remediation should reproduce the user's mode; got: {stderr}"
    );
}

// ── multi-pair ─────────────────────────────────────────────────────────

/// The offending pair is named by index, and — the larger behaviour change —
/// no earlier pair is processed. Before #363 the guard lived inside the per-pair
/// loop, so pair 1's `_val_*` files and trimming reports were written to disk
/// before pair 2 failed.
///
/// The four arguments must yield four *distinct output* paths. `BS-seq_10K_R2`
/// appears twice, which is fine — as pair 1's R2 and pair 2's R1 it produces
/// `..._val_2` and `..._val_1` respectively. What must not happen is an output
/// collision, because the pre-flight that detects one now runs *after* the
/// guard: the test would still fail, but on precedence rather than on pair
/// indexing.
#[test]
fn multi_pair_names_the_offending_pair_and_writes_nothing() {
    let dir = tempdir("multipair");
    let out = nonexistent_out(&dir);
    let (ok, stderr) = run(&[
        "--paired",
        fixture("BS-seq_10K_R1.fastq.gz").to_str().unwrap(),
        fixture("BS-seq_10K_R2.fastq.gz").to_str().unwrap(),
        fixture("BS-seq_10K_R2.fastq.gz").to_str().unwrap(),
        fixture("ubam_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!ok);
    // "Pair 2 of 2 is mixed", not bare "Pair 2 of 2" — the per-pair progress
    // banner prints `=== Pair 2 of 2 ===`, so the looser substring would also
    // be satisfied by a run that got as far as processing pair 2.
    assert!(
        stderr.contains("Pair 2 of 2 is mixed"),
        "the offending pair should be identified by index; got: {stderr}"
    );
    assert!(
        !out.exists(),
        "no pair may be processed before the rejection — pair 1's outputs and \
         trimming reports were written here before #363"
    );
}

// ── invocations that must keep working ─────────────────────────────────

/// The guard keys on BAM count, not format equality. `InputFormat` has three
/// variants, so a naive `formats[0] != formats[1]` would reject this — a pair
/// that works today.
#[test]
fn plain_plus_gzip_fastq_pair_is_still_accepted() {
    let dir = tempdir("plainplusgz");
    let plain = dir.join("plain_R1.fastq");
    {
        use std::io::Write;
        let f = std::fs::File::open(fixture("BS-seq_10K_R1.fastq.gz")).unwrap();
        let mut dec = flate2::read::MultiGzDecoder::new(f);
        let mut buf = Vec::new();
        std::io::Read::read_to_end(&mut dec, &mut buf).unwrap();
        std::fs::File::create(&plain)
            .unwrap()
            .write_all(&buf)
            .unwrap();
    }
    let out = dir.join("out");
    let (ok, stderr) = run(&[
        "--paired",
        plain.to_str().unwrap(),
        fixture("BS-seq_10K_R2.fastq.gz").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        ok,
        "a plain + gzipped FASTQ pair is legal and must not be rejected; stderr: {stderr}"
    );
    assert!(out.join("plain_R1_val_1.fq").exists(), "R1 output missing");
}

/// Single-file interleaved uBAM under `--paired` is the supported shape and
/// must not be caught by a guard aimed at two-file pairs.
#[test]
fn single_interleaved_ubam_still_accepted() {
    let dir = tempdir("interleaved");
    let out = dir.join("out");
    let (ok, stderr) = run(&[
        "--paired",
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        ok,
        "single interleaved uBAM under --paired must be accepted; stderr: {stderr}"
    );
}

/// `--hardtrim5` processes each input independently and never consults
/// `--paired`, so a mixed pair is harmless there. Left working deliberately;
/// see the #363 follow-up for whether the specialty modes should change.
#[test]
fn hardtrim_still_accepts_a_mixed_pair() {
    let dir = tempdir("hardtrim");
    let out = dir.join("out");
    let (ok, stderr) = run(&[
        "--hardtrim5",
        "10",
        "--paired",
        fixture("phred64_test.fastq").to_str().unwrap(),
        fixture("ubam_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        ok,
        "the paired format guard must not fire for --hardtrim5; stderr: {stderr}"
    );
}

/// #421 — `--paired` with one input file is the interleaved-uBAM path, which writes
/// FASTQ output through `run_paired_ubam_single_file` and reaches no `fastqc::run`.
#[test]
fn fastqc_refused_on_paired_interleaved_fastq_output() {
    let dir = tempdir("fastqc_arm5");
    let out = nonexistent_out(&dir);
    let (ok, stderr) = run(&[
        "--paired",
        "--fastqc",
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!ok, "--fastqc must be refused on this arm");
    assert!(
        stderr.contains("--fastqc is not supported") && stderr.contains("single input file"),
        "expected the FastQC refusal naming the shape; got: {stderr}"
    );
    // The remedy has to work for a user who simply forgot Read 2.
    assert!(
        stderr.contains("two files"),
        "message must offer the two-file remedy; got: {stderr}"
    );
    assert!(
        !out.exists(),
        "a refused run must not create its output dir"
    );
}

/// One flag apart from the case above, and capable: uBAM output routes above the
/// interleaved-FASTQ arm and runs FastQC on the single `_val.bam`.
#[test]
fn fastqc_accepted_on_paired_interleaved_ubam_output() {
    let dir = tempdir("fastqc_arm5_ubam");
    let (ok, stderr) = run(&[
        "--paired",
        "--output-format",
        "ubam",
        "--fastqc",
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        dir.to_str().unwrap(),
    ]);
    assert!(ok, "uBAM output is FastQC-capable; got: {stderr}");
    let zips: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .filter(|n| n.ends_with("_fastqc.zip"))
        .collect();
    assert_eq!(zips.len(), 1, "expected exactly one report, got {zips:?}");
}

/// `--clump_only --paired` with one input keeps its own diagnosis, which names the
/// actual remedy; the FastQC guard excludes it so that message still wins.
#[test]
fn clump_only_paired_single_file_keeps_its_own_message() {
    let dir = tempdir("fastqc_clump_prec");
    let out = nonexistent_out(&dir);
    let (ok, stderr) = run(&[
        "--clump_only",
        "--paired",
        "--fastqc",
        fixture("ubam_paired_test.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(!ok, "the combination is still refused");
    assert!(
        stderr.contains("requires two FASTQ input files"),
        "expected --clump_only's own message, not the FastQC one; got: {stderr}"
    );
    // Contrast with the sibling above: `Cli::validate` refuses before `ensure_output_dir`
    // (`main.rs:231` vs `:420`), this refusal after it. Which layer answers is the point.
    assert!(
        out.exists(),
        "main.rs's refusal runs after ensure_output_dir; if this flips, the guard changed layer"
    );
}
