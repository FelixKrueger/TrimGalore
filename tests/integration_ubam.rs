//! Binary-driven integration tests for uBAM input support (#316).
//!
//! These tests exercise `main.rs::run_single_file` / `run_paired_ubam_single_file`
//! end-to-end via the built `trim_galore` binary, covering wiring the lib
//! tests can't reach:
//!   * format detection → reader factory dispatch
//!   * output-path computation for `.bam` inputs (single-file paired-uBAM)
//!   * `--clumpify` + BAM (worker-pool dispatch through `RecordSource` trait)
//!   * specialty mode (`--hardtrim5`) + BAM (alternate entry points)
//!   * aligned-BAM rejection (per-record check in `BamReader::next_record`)
//!
//! Output-content parity is checked via `(id, seq, qual)` tuples — resilient
//! to gzip framing drift between runs. See PLAN §6.3 tier-1 assertion.

use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Parse a plain FASTQ file into `(id, seq, qual)` tuples. Used for
/// content-tuple parity (gzip-framing-resilient).
fn read_fastq_tuples(path: &Path) -> Vec<(String, String, String)> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path.display(), e));
    let mut tuples = Vec::new();
    let mut iter = text.lines();
    while let Some(id) = iter.next() {
        let seq = iter.next().unwrap_or_default();
        let _plus = iter.next().unwrap_or_default();
        let qual = iter.next().unwrap_or_default();
        tuples.push((id.to_string(), seq.to_string(), qual.to_string()));
    }
    tuples
}

#[test]
fn single_end_ubam_matches_reference() {
    let dir = fresh_tmpdir("tg_int_se_ubam");
    let status = Command::new(binary())
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got = read_fastq_tuples(&dir.join("ubam_test_trimmed.fq"));
    let want = read_fastq_tuples(Path::new("test_files/ubam_test_trimmed_REFERENCE.fq"));
    assert_eq!(
        got, want,
        "SE uBAM output did not match committed reference (record-by-record content-tuple comparison)"
    );
    assert_eq!(got.len(), 10, "expected 10 trimmed records");
}

#[test]
fn paired_end_interleaved_ubam_matches_reference() {
    let dir = fresh_tmpdir("tg_int_pe_ubam");
    let status = Command::new(binary())
        .args(["--paired"])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got_r1 = read_fastq_tuples(&dir.join("ubam_paired_test_val_1.fq"));
    let want_r1 = read_fastq_tuples(Path::new("test_files/ubam_paired_test_val_1_REFERENCE.fq"));
    let got_r2 = read_fastq_tuples(&dir.join("ubam_paired_test_val_2.fq"));
    let want_r2 = read_fastq_tuples(Path::new("test_files/ubam_paired_test_val_2_REFERENCE.fq"));

    assert_eq!(got_r1, want_r1, "PE R1 output diverged from reference");
    assert_eq!(got_r2, want_r2, "PE R2 output diverged from reference");
    assert_eq!(got_r1.len(), 10);
    assert_eq!(got_r2.len(), 10);

    // Same-template invariant: matched-pair IDs must align between R1/R2.
    for (a, b) in got_r1.iter().zip(got_r2.iter()) {
        assert_eq!(a.0, b.0, "R1/R2 records out of sync at id {}", a.0);
    }
}

#[test]
fn clumpify_plus_ubam_runs_clean() {
    // PLAN §3.4: --clumpify + BAM is allowed (works transparently via the
    // shared FastqRecord shape inside the worker pool). Smoke test: just
    // verify it exits non-zero-free and produces output.
    let dir = fresh_tmpdir("tg_int_clumpify_ubam");
    let status = Command::new(binary())
        .args(["--clumpify", "--cores", "2"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "--clumpify + uBAM smoke test failed; PLAN §3.4 says it should work"
    );
    assert!(
        dir.join("ubam_test_trimmed.fq.gz").exists() || dir.join("ubam_test_trimmed.fq").exists(),
        "expected trimmed output not found under {}",
        dir.display()
    );
}

#[test]
fn hardtrim5_plus_ubam_runs_clean() {
    // PLAN §3.4: specialty modes (--hardtrim5/3 / --clock / --implicon)
    // work with BAM input because they loop over `next_record()` via the
    // same RecordSource path. Smoke test the simplest specialty mode.
    let dir = fresh_tmpdir("tg_int_hardtrim5_ubam");
    let status = Command::new(binary())
        .args(["--hardtrim5", "20"])
        .arg("test_files/ubam_test.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "--hardtrim5 + uBAM smoke test failed; PLAN §3.4 says it should work"
    );
}

#[test]
fn passthrough_plus_ubam_rejected() {
    // PLAN §3.4: --passthrough + BAM is rejected in v1 (three-way ID sync
    // would need re-thinking for BAM record offsets). Confirm the error
    // surfaces with the expected message text.
    let dir = fresh_tmpdir("tg_int_passthrough_ubam_reject");
    let output = Command::new(binary())
        .args([
            "--paired",
            "--passthrough",
            "test_files/BS-seq_10K_I1.fastq.gz",
        ])
        .arg("test_files/ubam_paired_test.bam")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--passthrough + uBAM must be rejected"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    // Cli::validate's existing passthrough rule ("--passthrough requires
    // exactly one R1/R2 pair") fires FIRST when --paired+uBAM with N=1 is
    // combined with --passthrough, because validate() runs before
    // sanity_check_any. Either rejection path is acceptable — the contract
    // is "--passthrough + uBAM in v1 is rejected, with a clear pointer at
    // either fix". Accept either rejection message.
    assert!(
        stderr.contains("--passthrough")
            && (stderr.contains("not supported") || stderr.contains("requires exactly")),
        "expected some passthrough rejection message, got stderr: {}",
        stderr
    );
}

#[test]
fn preserve_tags_roundtrip_matches_golden() {
    // PLAN §3.2.5 + T23. The committed `ubam_test_with_tags.bam` carries
    // CB:Z:ATCGATCG-1 and UB:Z:GCTAGCTA aux tags on every record. With
    // `--preserve-tags CB,UB`, the FASTQ headers should be
    // `@<name>\tCB:Z:ATCGATCG-1\tUB:Z:GCTAGCTA` (samtools -T-compatible).
    // Verify against a committed golden reference.
    let dir = fresh_tmpdir("tg_int_preserve_tags");
    let status = Command::new(binary())
        .args(["--preserve-tags", "CB,UB"])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success(), "trim_galore exited non-zero");

    let got = read_fastq_tuples(&dir.join("ubam_test_with_tags_trimmed.fq"));
    let want = read_fastq_tuples(Path::new(
        "test_files/ubam_test_with_tags_trimmed_REFERENCE.fq",
    ));
    assert_eq!(
        got, want,
        "--preserve-tags output diverged from committed golden reference"
    );
    // Spot-check the tag-format invariant on the first record.
    assert!(
        got[0].0.contains("\tCB:Z:ATCGATCG-1"),
        "first record id must carry CB tag in user-specified order: {:?}",
        got[0].0
    );
    assert!(
        got[0].0.contains("\tUB:Z:GCTAGCTA"),
        "first record id must carry UB tag: {:?}",
        got[0].0
    );
}

#[test]
fn preserve_tags_user_specified_order_honoured() {
    // Same fixture, but user requests UB then CB — the resulting header must
    // be `...\tUB:...\tCB:...`. Tag order is policy-defined by the flag, NOT
    // by BAM file aux-field order.
    let dir = fresh_tmpdir("tg_int_preserve_tags_order");
    let status = Command::new(binary())
        .args(["--preserve-tags", "UB,CB"])
        .arg("test_files/ubam_test_with_tags.bam")
        .arg("-o")
        .arg(&dir)
        .status()
        .expect("trim_galore failed to run");
    assert!(status.success());
    let got = read_fastq_tuples(&dir.join("ubam_test_with_tags_trimmed.fq"));
    // UB must come BEFORE CB in the header now.
    let id = &got[0].0;
    let ub_pos = id.find("UB:Z:").expect("UB tag missing");
    let cb_pos = id.find("CB:Z:").expect("CB tag missing");
    assert!(
        ub_pos < cb_pos,
        "tag order must match --preserve-tags argument (UB before CB), got id: {:?}",
        id
    );
}

#[test]
fn paired_with_single_fastq_rejected() {
    // PLAN §3.3: --paired with a single input file is only legal if that
    // file is a uBAM. A single FASTQ in --paired mode must be rejected.
    let dir = fresh_tmpdir("tg_int_paired_single_fastq_reject");
    let output = Command::new(binary())
        .args(["--paired"])
        .arg("test_files/BS-seq_10K_R1.fastq.gz")
        .arg("-o")
        .arg(&dir)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "--paired with single FASTQ must be rejected"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("only legal if that file is a uBAM"),
        "expected single-FASTQ rejection message, got stderr: {}",
        stderr
    );
}

// ---------------------------------------------------------------------------
// #415 — whitespace in a BAM read name, and framing bytes in a Z tag value.
//
// Each case names the entry point it pins. Three entry points run before the
// trimming loop (sanity check, adapter detection, poly-G scan), so a fixture
// whose first record is the offender never reaches the writer; `-a` plus
// `--no_poly_g` skips both scans and is the arm a >1M-read file always takes.
// ---------------------------------------------------------------------------

/// Skip adapter auto-detection and the poly-G scan so the trimming loop is reached.
const SKIP_PRESCANS: [&str; 3] = ["-a", "AGATCGGAAGAGC", "--no_poly_g"];

#[test]
fn ws_qname_refused_at_sanity_check_entry_point() {
    let dir = fresh_tmpdir("tg_415_sanity");
    let output = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB"])
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_ws_qname.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(!output.status.success(), "whitespace QNAME must be refused");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace") && stderr.contains("name with space"),
        "message must name the defect and the offending name, got: {}",
        stderr
    );
    // Record 1 is the offender, so the refusal precedes any writer.
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .filter(|n| n.contains("_trimmed") || n.contains("_val"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "refusal at the sanity check must write nothing, found: {:?}",
        leftovers
    );
}

#[test]
fn ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial() {
    let dir = fresh_tmpdir("tg_415_loop_ubam");
    let output = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_ws_qname_late.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(!output.status.success(), "whitespace QNAME must be refused");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("BAM record 2"),
        "must be caught in the trimming loop at record 2, not at the sanity check, got: {}",
        stderr
    );
    // The writer opens before the read loop, so record 1 is already on disk.
    // Asserted, not desired — the cure spans all per-record bails and is tracked
    // separately.
    let partial = dir.join("ubam_ws_qname_late_trimmed.bam");
    let bytes = std::fs::read(&partial).expect("a mid-stream refusal leaves a partial uBAM behind");
    // What makes the residue a data-integrity problem is that it is indistinguishable
    // from a complete BAM: bgzf's Drop finalises it, EOF marker included.
    const BGZF_EOF: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    assert!(
        bytes.ends_with(BGZF_EOF),
        "the partial carries a valid BGZF EOF marker, so nothing downstream flags it ({} bytes)",
        bytes.len()
    );
}

#[test]
fn ws_qname_refused_in_trimming_loop_fastq_out() {
    let dir = fresh_tmpdir("tg_415_loop_fastq");
    let output = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_ws_qname_late.bam")
        .output()
        .expect("trim_galore failed to run");
    // The widened scope: FASTQ output loses no tag text today, and is refused anyway.
    assert!(
        !output.status.success(),
        "the refusal is read-side, so FASTQ output is refused too"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace"),
        "expected the whitespace message on the FASTQ-output path, got: {}",
        stderr
    );
}

#[test]
fn ws_qname_refused_via_threaded_reader() {
    let dir = fresh_tmpdir("tg_415_threaded");
    let output = Command::new(binary())
        .args(["--cores", "2"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_ws_qname_late.bam")
        .output()
        .expect("trim_galore failed to run");
    // `--cores 2` is the only route to the threaded single-stream reader.
    assert!(!output.status.success(), "whitespace QNAME must be refused");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace"),
        "expected the whitespace message from the threaded reader, got: {}",
        stderr
    );
}

#[test]
fn ws_qname_refused_via_interleaved_deinterleaver() {
    let dir = fresh_tmpdir("tg_415_deinterleaved");
    let output = Command::new(binary())
        .arg("--paired")
        .args(SKIP_PRESCANS)
        .args(["-a2", "AGATCGGAAGAGC"])
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_ws_qname_paired.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(!output.status.success(), "whitespace QNAME must be refused");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace"),
        "expected the whitespace message, got: {}",
        stderr
    );
    // Record 3 is the first read of the second pair — proves the de-interleaver
    // reported it, not the sanity check (which only ever sees record 1).
    assert!(
        stderr.contains("BAM record 3"),
        "must come from the de-interleaver at record 3, got: {}",
        stderr
    );
}

#[test]
fn lf_in_qname_refused_end_to_end() {
    let dir = fresh_tmpdir("tg_415_lf_qname");
    let output = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_lf_qname.bam")
        .output()
        .expect("trim_galore failed to run");
    // A newline in a read name would split one record across five lines,
    // desynchronising every 4-line-block reader after it.
    assert!(!output.status.success(), "LF in a QNAME must be refused");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace") && stderr.contains("readB\\nEVIL"),
        "message must show the newline escaped, got: {}",
        stderr
    );
}

#[test]
fn lf_in_tag_value_refused_end_to_end() {
    let dir = fresh_tmpdir("tg_415_lf_tag");
    let output = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_lf_tagvalue.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "a newline in a preserved Z value corrupts framing the same way a QNAME newline does"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    // `aux tag 'CB'` comes from the call site, not the predicate — asserting the bare
    // prefix would also match the predicate's own wording and pin nothing.
    assert!(
        stderr.contains("aux tag 'CB'") && stderr.contains("AAA\\nCCC"),
        "message must name the tag and show the escaped value, got: {}",
        stderr
    );
}

#[test]
fn lf_in_character_tag_value_refused_end_to_end() {
    // `A` values reach the same id as `Z` values through a different match arm.
    let dir = fresh_tmpdir("tg_415_lf_atag");
    let output = Command::new(binary())
        .args(["--preserve-tags", "XA"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_lf_atag.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        !output.status.success(),
        "a newline in an A tag value corrupts framing the same way a Z value does"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("aux tag 'XA'") && stderr.contains("cannot pass through"),
        "message must name the tag and the framing failure, got: {}",
        stderr
    );
}

#[test]
fn printable_character_tags_round_trip() {
    // The acceptance twin for the A arm — green whether or not the guard is present,
    // so it pins carriage rather than refusal.
    let dir = fresh_tmpdir("tg_415_atag_ok");
    let output = Command::new(binary())
        .args(["--preserve-tags", "XA"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_atag_ok.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "printable A values must not be refused, stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let tuples = read_fastq_tuples(&dir.join("ubam_atag_ok_trimmed.fq"));
    assert_eq!(tuples.len(), 3);
    let tags: Vec<&str> = tuples.iter().map(|t| t.0.as_str()).collect();
    assert!(
        tags[0].ends_with("XA:A:+") && tags[1].ends_with("XA:A:-"),
        "A values must reach the header intact, got: {:?}",
        tags
    );
}

#[test]
fn clean_tag_value_with_space_still_accepted() {
    // A space is legal in a `Z` value, and the tab introducing the tag precedes
    // it — so the tag tail still parses. Guards against reusing the QNAME set.
    let dir = fresh_tmpdir("tg_415_tag_space");
    let output = Command::new(binary())
        .args(["--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg("test_files/ubam_tagvalue_space.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "a space in a Z tag value must stay legal, stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let tuples = read_fastq_tuples(&dir.join("ubam_tagvalue_space_trimmed.fq"));
    assert_eq!(tuples.len(), 1);
    assert!(
        tuples[0].0.contains("CB:Z:has a space"),
        "the space-bearing tag must survive, got id: {:?}",
        tuples[0].0
    );

    // uBAM out is the harder direction: `parse_name_and_data` has to re-split the tail
    // around the space rather than treat it as a description boundary.
    let ubam_dir = fresh_tmpdir("tg_415_tag_space_ubam");
    let output = Command::new(binary())
        .args(["--output-format", "ubam", "--preserve-tags", "CB"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&ubam_dir)
        .arg("test_files/ubam_tagvalue_space.bam")
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "uBAM round-trip of a space-bearing Z value must succeed, stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut decoded = Vec::new();
    std::io::Read::read_to_end(
        &mut noodles::bgzf::Reader::new(
            std::fs::File::open(ubam_dir.join("ubam_tagvalue_space_trimmed.bam")).unwrap(),
        ),
        &mut decoded,
    )
    .unwrap();
    assert!(
        decoded.windows(15).any(|w| w == b"CBZhas a space\0"),
        "the space-bearing Z value must round-trip into the output BAM's aux data"
    );
}

#[test]
fn description_notice_wording_is_fastq_input_only() {
    // No space reaches `parse_name_and_data` from BAM input, so #406's "FASTQ
    // header text" wording is unconditionally accurate; this pins the one
    // direction that still reaches it.
    let dir = fresh_tmpdir("tg_415_notice");
    let fq = dir.join("desc.fastq");
    std::fs::write(
        &fq,
        "@readA 1:N:0:CGATCG\n\
         ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
         +\n\
         IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args(["--output-format", "ubam"])
        .args(SKIP_PRESCANS)
        .arg("-o")
        .arg(&dir)
        .arg(&fq)
        .output()
        .expect("trim_galore failed to run");
    assert!(
        output.status.success(),
        "a FASTQ description is legitimate and must still be dropped with a notice, stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("FASTQ header text"),
        "the #406 notice must still fire on FASTQ input, got: {}",
        stderr
    );
}
