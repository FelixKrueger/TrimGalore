//! #428 — a refused run must publish no output file.
//!
//! The five per-record bails in `bam_record_to_fastq` sit on the read side, so
//! they fire on every output path. Each case here drives one output arm to a
//! mid-stream refusal and asserts that neither the final name nor the
//! temporary survives.
//!
//! Two properties are easy to conflate and both matter: that the final name is
//! absent (not merely empty — the `--cores 2` arm used to leave a valid 0-byte
//! FASTQ), and that no `.partial` is left in the output directory.
//!
//! The success-path cases guard the other direction: publishing happens by
//! renaming *after* the writer is closed, so a slip in that order would ship a
//! trailerless gzip or an EOF-less BAM under the final name.

use std::path::{Path, PathBuf};
use std::process::Command;

use trim_galore::fastq::FastqReader;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Skip adapter auto-detection and the poly-G scan so the trimming loop is
/// reached; both pre-scans read up to 1M records and would catch the offender
/// before any writer exists.
const SKIP_PRESCANS: [&str; 3] = ["-a", "AGATCGGAAGAGC", "--no_poly_g"];

/// Every file left in `dir`, sorted.
fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .collect();
    names.sort();
    names
}

/// Assert the refusal published nothing: no final name, and no temporary.
fn assert_nothing_published(dir: &Path, expected_output: &str) {
    let out = dir.join(expected_output);
    assert!(
        !out.exists(),
        "a refused run must publish no output file; {} exists at {} bytes",
        expected_output,
        std::fs::metadata(&out).map(|m| m.len()).unwrap_or(0)
    );
    let leftovers: Vec<String> = listing(dir)
        .into_iter()
        .filter(|n| n.ends_with(".partial"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "the temporary must be removed on the error path, found: {:?}",
        leftovers
    );
}

/// Run the binary on the offending fixture with `extra` flags, asserting it is
/// refused at `record` — i.e. inside the trimming loop, past writer creation.
///
/// The ordinal is what keeps every absence assertion in this file non-vacuous:
/// if the whitespace check ever moved to `sanity_check_any`, they would all pass
/// with no writer ever created.
fn refuse(dir: &Path, extra: &[&str], fixture: &str, record: u32) {
    let output = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(extra)
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(dir)
        .arg(fixture)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "the whitespace read name must be refused; stdout: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace"),
        "expected the whitespace refusal, got: {stderr}"
    );
    assert!(
        stderr.contains(&format!("BAM record {record}")),
        "must be refused inside the trimming loop at record {record}, got: {stderr}"
    );
}

/// Trailing bytes of a complete BGZF stream.
const BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const OFFENDER: &str = "test_files/ubam_ws_qname_late.bam";
/// 5000 clean records ahead of the offender, so the unpublished residue would
/// have been large and plausible rather than a single record.
const OFFENDER_BULK: &str = "test_files/ubam_ws_qname_bulk.bam";

// ─── T3: every measured arm publishes nothing ───────────────────────────────

#[test]
fn refused_se_fastq_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_se_fastq");
    refuse(&dir, &[], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.fq");
}

#[test]
fn refused_se_ubam_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_se_ubam");
    refuse(&dir, &["--output-format", "ubam"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.bam");
}

/// The parallel path created the output before any worker produced a chunk, so
/// the residue was a valid *empty* FASTQ — "no reads survived", not "failed".
#[test]
fn refused_parallel_publishes_nothing_not_an_empty_file() {
    let dir = fresh_tmpdir("tg_428_cores2");
    refuse(&dir, &["--cores", "2"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.fq");
}

/// `--clumpify` is refused only under `--output-format ubam`, so with FASTQ
/// output on a uBAM input it runs and is a live arm.
#[test]
fn refused_clumpify_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_clumpify");
    refuse(&dir, &["--clumpify", "--cores", "2"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.fq");
}

#[test]
fn refused_hardtrim5_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_ht5");
    refuse(&dir, &["--hardtrim5", "10"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late.10bp_5prime.fq");
}

#[test]
fn refused_hardtrim3_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_ht3");
    refuse(&dir, &["--hardtrim3", "10"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late.10bp_3prime.fq");
}

#[test]
fn refused_hardtrim5_ubam_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_ht5_bam");
    refuse(
        &dir,
        &["--hardtrim5", "10", "--output-format", "ubam"],
        OFFENDER,
        2,
    );
    assert_nothing_published(&dir, "ubam_ws_qname_late.10bp_5prime.bam");
}

#[test]
fn refused_clump_only_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_clump_only");
    let output = Command::new(binary())
        .args([
            "--preserve-tags",
            "CB",
            "--clump_only",
            "--output-format",
            "ubam",
        ])
        .arg("-o")
        .arg(&dir)
        .arg(OFFENDER)
        .output()
        .expect("trim_galore failed to run");
    assert!(!output.status.success(), "must be refused");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("BAM record 2") && stderr.contains("read name contains whitespace"),
        "expected the mid-loop whitespace refusal, got: {stderr}"
    );
    assert_nothing_published(&dir, "ubam_ws_qname_late_clumped.bam");
}

/// `--demux` reads the trimmed output back by its final path and fans out more
/// files from it, so it is the arm where a late commit is most likely to bite.
/// The refusal happens at the trimming stage, before any barcode file exists.
#[test]
fn refused_demux_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_demux");
    refuse(
        &dir,
        &["--demux", "test_files/demux_test_samplesheet.txt"],
        OFFENDER,
        2,
    );
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.fq");
}

// ─── T6: the paired arm, and a large plausible residue ──────────────────────

/// For paired uBAM the bail site is the de-interleaver's producer thread, which
/// converts records before the trimming loop ever sees them.
#[test]
fn refused_paired_fastq_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_paired_fq");
    refuse(
        &dir,
        &["--paired"],
        "test_files/ubam_ws_qname_paired.bam",
        3,
    );
    assert_nothing_published(&dir, "ubam_ws_qname_paired_val_1.fq");
    assert_nothing_published(&dir, "ubam_ws_qname_paired_val_2.fq");
}

#[test]
fn refused_paired_ubam_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_paired_bam");
    refuse(
        &dir,
        &["--paired", "--output-format", "ubam"],
        "test_files/ubam_ws_qname_paired.bam",
        3,
    );
    assert_nothing_published(&dir, "ubam_ws_qname_paired_val.bam");
}

/// The shape that makes #428 a data-integrity bug rather than a cosmetic one:
/// thousands of good records are written before the bail, so the residue was
/// large, well-formed and indistinguishable from a short successful run.
#[test]
fn refused_after_thousands_of_good_records_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_bulk");
    refuse(&dir, &[], OFFENDER_BULK, 5001);
    assert_nothing_published(&dir, "ubam_ws_qname_bulk_trimmed.fq");
    // Nothing at all beyond the temporary's absence — no reports either, since
    // those are written after the loop.
    assert!(
        listing(&dir).is_empty(),
        "a refused run leaves the output directory empty, found: {:?}",
        listing(&dir)
    );
}

// ─── T4: a previous good output survives a failed re-run ────────────────────

/// The property Option A could not have delivered. `File::create` used to
/// truncate the previous output before the first record was read, so a re-run
/// that then failed destroyed a complete result and left a plausible short one
/// in its place. The re-run must target the *same* output path, which means the
/// same input filename — otherwise the two runs write different files and the
/// test passes vacuously.
#[test]
fn failed_rerun_leaves_the_previous_good_output_intact() {
    let dir = fresh_tmpdir("tg_428_rerun");
    let subject = dir.join("subject.bam");

    std::fs::copy("test_files/ubam_test_with_tags.bam", &subject).unwrap();
    let good = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg(&subject)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        good.status.success(),
        "the baseline run must succeed; stderr: {}",
        String::from_utf8_lossy(&good.stderr)
    );

    let out = dir.join("subject_trimmed.fq");
    let report = dir.join("subject.bam_trimming_report.txt");
    let before_out = std::fs::read(&out).expect("baseline run wrote an output");
    let before_report = std::fs::read(&report).expect("baseline run wrote a report");
    assert!(!before_out.is_empty(), "baseline output must be non-empty");

    // Same input name, so the same output path; different content, so it fails.
    std::fs::copy(OFFENDER, &subject).unwrap();
    let failed = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg(&subject)
        .output()
        .expect("trim_galore failed to run");
    assert!(!failed.status.success(), "the re-run must be refused");

    assert_eq!(
        before_out,
        std::fs::read(&out).expect("the previous output must still exist"),
        "a failed re-run must not touch the previous run's output"
    );
    assert_eq!(
        before_report,
        std::fs::read(&report).expect("the previous report must still exist"),
        "a failed re-run must not touch the previous run's report"
    );
}

// ─── T2: publishing must not lose the trailer ───────────────────────────────

/// Count records by reading the output through the same decoder the tool uses.
/// A gzip stream missing its CRC/ISIZE trailer fails here, which record-counting
/// against a permissive decoder would not catch.
fn read_back(path: &Path) -> usize {
    let mut reader = FastqReader::open(path).expect("output must be readable");
    let mut n = 0;
    while reader
        .next_record()
        .expect("the output must decode cleanly to EOF")
        .is_some()
    {
        n += 1;
    }
    n
}

/// Serial gzip: the trailer comes from `GzEncoder`'s own teardown, so a rename
/// ordered before that teardown would publish a truncated stream.
#[test]
fn published_serial_gzip_output_decodes_to_eof() {
    let dir = fresh_tmpdir("tg_428_gz_serial");
    let output = Command::new(binary())
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .output()
        .expect("trim_galore failed to run");
    assert!(output.status.success(), "the run must succeed");
    let out = dir.join("BS-seq_10K_R1_trimmed.fq.gz");
    assert!(out.exists(), "gzipped output must be published");
    assert!(read_back(&out) > 9000, "expected ~10K records");
}

/// `--cores 2` with gzip output builds a `gzp::ParCompress`, whose teardown is
/// a different mechanism from `GzEncoder`'s — it joins compressor threads and
/// flushes their blocks. This is the only arm that exercises it.
#[test]
fn published_parallel_gzip_output_decodes_to_eof() {
    let dir = fresh_tmpdir("tg_428_gz_par");
    let output = Command::new(binary())
        .args(["--hardtrim5", "30", "--cores", "2"])
        .arg("-o")
        .arg(&dir)
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "the run must succeed; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let out = dir.join("BS-seq_10K_R1.30bp_5prime.fq.gz");
    assert!(out.exists(), "gzipped output must be published");
    assert!(read_back(&out) > 9000, "expected ~10K records");
}

/// uBAM: the BGZF EOF marker is written by `try_finish` before the rename.
#[test]
fn published_ubam_output_carries_the_bgzf_eof_marker() {
    let dir = fresh_tmpdir("tg_428_bgzf_eof");
    let output = Command::new(binary())
        .args(["--output-format", "ubam"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_test_with_tags.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "the run must succeed; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let out = dir.join("ubam_test_with_tags_trimmed.bam");
    let bytes = std::fs::read(&out).expect("uBAM output must be published");
    assert!(
        bytes.ends_with(BGZF_EOF),
        "published BAM must end with the BGZF EOF marker ({} bytes)",
        bytes.len()
    );
}

// ─── T10: the published name is what downstream stages consume ──────────────

/// FastQC derives its report name from the path it is handed, so a commit that
/// happened after `fastqc::run` would name the report from the temporary.
#[test]
fn fastqc_report_is_named_from_the_published_output() {
    let dir = fresh_tmpdir("tg_428_fastqc_name");
    let output = Command::new(binary())
        .args(["--fastqc"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_test_with_tags.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "the run must succeed; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let names = listing(&dir);
    assert!(
        names
            .iter()
            .any(|n| n == "ubam_test_with_tags_trimmed_fastqc.html"),
        "FastQC must be named from the published output, got: {names:?}"
    );
    assert!(
        !names
            .iter()
            .any(|n| n.starts_with('.') || n.contains(".partial")),
        "no temporary may survive a successful run, got: {names:?}"
    );
}

// ─── §2.4's remaining arm ───────────────────────────────────────────────────

/// Same `gzip == false` branch as the default arm for a `.bam` input, but §2.4
/// lists it as its own row.
#[test]
fn refused_dont_gzip_publishes_nothing() {
    let dir = fresh_tmpdir("tg_428_dont_gzip");
    refuse(&dir, &["--dont_gzip"], OFFENDER, 2);
    assert_nothing_published(&dir, "ubam_ws_qname_late_trimmed.fq");
}

// ─── T7: a legitimately empty output is still published ─────────────────────

/// The one way this change could break a green pipeline: "publish nothing on a
/// refusal" must not become "publish nothing when every read was filtered".
#[test]
fn all_reads_filtered_still_publishes_an_empty_output() {
    let dir = fresh_tmpdir("tg_428_all_filtered");
    let output = Command::new(binary())
        .args(["--length", "500"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "filtering every read is success, not failure; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("100.0%"),
        "the fixture must filter every read for this test to mean anything"
    );

    let out = dir.join("BS-seq_10K_R1_trimmed.fq.gz");
    assert!(out.exists(), "an empty output must still be published");
    assert_eq!(read_back(&out), 0, "the published output holds no records");
    assert!(
        !listing(&dir).iter().any(|n| n.ends_with(".partial")),
        "no temporary may survive"
    );
}

/// Same property on the uBAM path, where "empty" is header-only.
#[test]
fn all_reads_filtered_still_publishes_an_empty_ubam() {
    let dir = fresh_tmpdir("tg_428_all_filtered_bam");
    let output = Command::new(binary())
        .args(["--output-format", "ubam", "--length", "500"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_test_with_tags.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "filtering every read is success, not failure; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = std::fs::read(dir.join("ubam_test_with_tags_trimmed.bam"))
        .expect("an empty uBAM must still be published");
    assert!(
        bytes.ends_with(BGZF_EOF),
        "the empty uBAM is still a complete BAM ({} bytes)",
        bytes.len()
    );
}

// ─── The commit that no test used to assert ─────────────────────────────────

/// `--retain_unpaired`'s two writers are the one set whose missed commit the
/// suite could not see: the pre-flight tripwire asserts `created ⊆ planned`, and
/// over-planning is legal there by design, so an absent output passes.
#[test]
fn retain_unpaired_publishes_both_rescued_outputs() {
    let dir = fresh_tmpdir("tg_428_retain_unpaired");
    let r1 = dir.join("s_R1.fastq");
    let r2 = dir.join("s_R2.fastq");
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..40 {
        a.push_str(&format!(
            "@p{i}\n{}\n+\n{}\n",
            "A".repeat(80),
            "I".repeat(80)
        ));
        b.push_str(&format!("@p{i}\n{}\n+\n{}\n", "A".repeat(5), "I".repeat(5)));
    }
    std::fs::write(&r1, a).unwrap();
    std::fs::write(&r2, b).unwrap();

    let output = Command::new(binary())
        .args(["--paired", "--retain_unpaired"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "the run must succeed; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Every R2 is 5 bp against the default --length 20, so every pair strands
    // and every R1 is rescued.
    let unpaired_r1 = dir.join("s_R1_unpaired_1.fq");
    assert!(
        unpaired_r1.exists(),
        "the rescued R1 output must be published, got: {:?}",
        listing(&dir)
    );
    assert_eq!(read_back(&unpaired_r1), 40, "every R1 is rescued");
    assert!(
        dir.join("s_R1_val_1.fq").exists() && dir.join("s_R2_val_2.fq").exists(),
        "the validated outputs are published even though they are empty, got: {:?}",
        listing(&dir)
    );
}
