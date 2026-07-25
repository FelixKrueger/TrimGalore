//! Binary-driven integration tests for `--clump_only --output-format ubam` (v2).
//!
//! Extends the v1 FASTQ-only integration tests with the uBAM in/out paths.
//! Uses the same BAM-decode helpers as `integration_ubam_out.rs` (noodles
//! `bam::Reader` + tuple projection + `@PG`-ignoring header comparison).
//!
//! Fixture: `test_files/ubam_test.bam` (SE uBAM),
//! `test_files/ubam_paired_test.bam` (PE interleaved uBAM),
//! `test_files/BS-seq_10K_R{1,2}.fastq.gz` (paired FASTQ).

use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("tg_clump_only_ubam_{slug}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("test_files")
        .join(name)
}

/// One BAM record projected into comparable form: (name, flag-bits, seq bytes,
/// qual raw-Phred bytes, sorted-by-tag aux byte tuples). Sorting aux tags
/// normalises map-iteration order in equality checks.
type RecordTuple = (Vec<u8>, u16, Vec<u8>, Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>);

fn bam_tuples(path: &Path) -> Vec<RecordTuple> {
    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("failed to open {}: {}", path.display(), e));
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    let _header = reader.read_header().expect("BAM header parse");
    let mut tuples = Vec::new();
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("BAM record read") > 0 {
        let name: Vec<u8> = rec
            .name()
            .map(|bn| AsRef::<[u8]>::as_ref(bn).to_vec())
            .unwrap_or_default();
        let flag = rec.flags().bits();
        let seq: Vec<u8> = rec.sequence().iter().collect();
        let qual: Vec<u8> = rec.quality_scores().as_ref().to_vec();
        let mut aux: Vec<(Vec<u8>, Vec<u8>)> = rec
            .data()
            .iter()
            .filter_map(|res| {
                let (tag, value) = res.ok()?;
                let tag_bytes = vec![tag.as_ref()[0], tag.as_ref()[1]];
                let value_dbg = format!("{:?}", value).into_bytes();
                Some((tag_bytes, value_dbg))
            })
            .collect();
        aux.sort_by(|a, b| a.0.cmp(&b.0));
        tuples.push((name, flag, seq, qual, aux));
    }
    tuples
}

/// Sort BAM record tuples canonically for multiset comparison after a reorder.
/// The mode is a *permutation* — sort-both-sides-then-compare is the load-
/// bearing test.
fn sorted_bam_tuples(path: &Path) -> Vec<RecordTuple> {
    let mut t = bam_tuples(path);
    t.sort();
    t
}

/// Count `@PG` records in a BAM header. Used to verify the input `@PG` chain
/// is preserved (existing @PG lines) AND our TrimGalore `@PG` is appended.
fn count_pg_lines(path: &Path) -> usize {
    let file = std::fs::File::open(path).expect("open BAM");
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    let header = reader.read_header().expect("BAM header parse");
    header.programs().as_ref().len()
}

/// Assert the output BAM header carries a `@PG ID:trim_galore` record.
fn has_trim_galore_pg(path: &Path) -> bool {
    let file = std::fs::File::open(path).expect("open BAM");
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    let header = reader.read_header().expect("BAM header parse");
    header
        .programs()
        .as_ref()
        .keys()
        .any(|id| AsRef::<[u8]>::as_ref(id) == b"trim_galore")
}

// ── SE positive path ──────────────────────────────────────────────

#[test]
fn se_ubam_out_from_fastq_in() {
    let dir = fresh_tmpdir("se_fq_to_bam");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    assert!(input.exists(), "fixture missing: {}", input.display());

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()
        .expect("trim_galore failed to run");
    assert!(
        status.success(),
        "--clump_only --output-format ubam SE failed"
    );

    let out = dir.join("BS-seq_10K_R1_clumped.bam");
    assert!(out.exists(), "expected BAM missing: {}", out.display());

    // Record multiset preserved: the output is a permutation of the input.
    // Reference: decoded input via a quick FASTQ→BAM conversion isn't
    // available cheaply here, so we just verify the record count matches
    // (10K records) and that reads decode cleanly.
    let output_tuples = bam_tuples(&out);
    assert_eq!(
        output_tuples.len(),
        10_000,
        "expected 10000 records in BAM output, got {}",
        output_tuples.len()
    );
    // Sanity: every record should have SE flag bits set (paired=0, unmapped=1
    // → 0x4). The BAM SE flag is exactly 0x0004.
    for rec in &output_tuples {
        assert_eq!(
            rec.1, 4,
            "expected SE unmapped flag (0x4), got {:#x}",
            rec.1
        );
    }
}

#[test]
fn se_ubam_in_ubam_out_pg_chain() {
    let dir = fresh_tmpdir("se_pg_chain");
    let input = fixture("ubam_test.bam");
    assert!(input.exists());

    // Count input @PG lines so we can assert our appended line is additive.
    let input_pg_count = count_pg_lines(&input);

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()
        .expect("trim_galore failed");
    assert!(status.success());

    let out = dir.join("ubam_test_clumped.bam");
    assert!(out.exists());

    let output_pg_count = count_pg_lines(&out);
    assert_eq!(
        output_pg_count,
        input_pg_count + 1,
        "expected output @PG count = input + 1 (TrimGalore's appended line); got in={} out={}",
        input_pg_count,
        output_pg_count
    );
    assert!(
        has_trim_galore_pg(&out),
        "output BAM header should contain @PG ID:trim_galore"
    );
    // Record multiset preserved.
    assert_eq!(sorted_bam_tuples(&input), sorted_bam_tuples(&out));
}

// ── PE positive path ──────────────────────────────────────────────

#[test]
fn pe_ubam_out_interleaved_from_fastq_pair() {
    let dir = fresh_tmpdir("pe_fq_to_bam");
    let r1 = fixture("BS-seq_10K_R1.fastq.gz");
    let r2 = fixture("BS-seq_10K_R2.fastq.gz");
    assert!(r1.exists() && r2.exists());

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&r1)
        .arg(&r2)
        .status()
        .expect("trim_galore failed");
    assert!(status.success());

    let out = dir.join("BS-seq_10K_R1_clumped.bam");
    assert!(out.exists(), "expected interleaved PE BAM missing");

    let tuples = bam_tuples(&out);
    assert_eq!(tuples.len(), 20_000, "expected 10K pairs = 20K records");

    // Mate-adjacency: record[2i] has PAIRED|READ1 flag = 0x4D (77);
    // record[2i+1] has PAIRED|READ2 flag = 0x8D (141). Names match at
    // pair boundaries.
    for i in 0..10_000 {
        let r1_rec = &tuples[i * 2];
        let r2_rec = &tuples[i * 2 + 1];
        assert_eq!(
            r1_rec.1,
            0x4D,
            "record {} (R1 slot) flag 0x{:x} != expected 0x4D",
            i * 2,
            r1_rec.1
        );
        assert_eq!(
            r2_rec.1,
            0x8D,
            "record {} (R2 slot) flag 0x{:x} != expected 0x8D",
            i * 2 + 1,
            r2_rec.1
        );
        // Mate names identical (bcl2fastq etc. strip the /1 /2 suffixes when
        // writing to BAM; whatever the convention, the name portion should
        // match for a pair).
        assert_eq!(
            r1_rec.0, r2_rec.0,
            "pair at position {} has mismatched names: {:?} vs {:?}",
            i, r1_rec.0, r2_rec.0
        );
    }
}

#[test]
fn pe_ubam_in_interleaved_output() {
    let dir = fresh_tmpdir("pe_bam_in_out");
    let input = fixture("ubam_paired_test.bam");
    assert!(input.exists());

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()
        .expect("trim_galore failed");
    assert!(
        status.success(),
        "--clump_only --paired --output-format ubam on interleaved uBAM failed"
    );

    let out = dir.join("ubam_paired_test_clumped.bam");
    assert!(out.exists());

    // Record multiset preserved (record count matches; sorted tuples match).
    assert_eq!(sorted_bam_tuples(&input), sorted_bam_tuples(&out));
}

// ── Rejection matrix ──────────────────────────────────────────────

#[test]
fn rejects_dont_gzip_with_ubam_output() {
    // Hoisted into §3.4a — applies to all uBAM output paths (not just
    // --clump_only). This is the specific --clump_only invocation.
    let dir = fresh_tmpdir("rej_dont_gzip");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--output-format",
            "ubam",
            "--dont_gzip",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()
        .expect("trim_galore failed");
    assert!(
        !out.status.success(),
        "--dont_gzip + --output-format ubam must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--dont_gzip") && stderr.contains("BGZF"),
        "stderr should mention --dont_gzip and BGZF; got: {stderr}"
    );
}

#[test]
fn rejects_two_bam_paired() {
    // Shape A + two DISTINCT BAM files is ambiguous — samtools convention is
    // single-file interleaved. Match trim uBAM path's rejection.
    //
    // Regression guard for Reviewer A finding C-1: the previous version of
    // this test passed the same BAM twice, tripping Cli::validate's R1==R2
    // dup-check FIRST — leaving the actual two-BAM Shape A guard at
    // main.rs untested. Copy the fixture to a second distinct path so R1
    // and R2 are different files, exercising the real guard.
    let dir = fresh_tmpdir("rej_two_bam");
    let bam = fixture("ubam_test.bam");
    let bam_copy = dir.join("ubam_test_copy.bam");
    std::fs::copy(&bam, &bam_copy).expect("failed to copy fixture");

    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&bam)
        .arg(&bam_copy)
        .output()
        .expect("trim_galore failed");
    assert!(
        !out.status.success(),
        "two distinct BAM files under --paired must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    // The actual two-BAM guard's message mentions the "single interleaved"
    // convention — this pins the rejection to the guard we added, not to
    // the R1==R2 dup check (which mentions "Read 1 and Read 2 appear to be").
    assert!(
        stderr.contains("single interleaved") || stderr.contains("uBAM paired mode expects"),
        "stderr should mention single-interleaved-file convention; got: {stderr}"
    );
}

/// Regression guard for Reviewer B finding C-1: `--clump_only --paired
/// <bam>` (WITHOUT `--output-format ubam`) previously panicked with
/// index-out-of-bounds — the general N=1 carve-out under --paired
/// accepted a single BAM, dispatch entered `run_specialty_paired`, and
/// `chunk[1]` on a length-1 chunk crashed. Fix: reject at dispatch time
/// in the FASTQ output arm.
#[test]
fn rejects_paired_single_bam_without_output_format_ubam() {
    let dir = fresh_tmpdir("rej_paired_single_bam_no_fmt");
    let input = fixture("ubam_test.bam");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()
        .expect("trim_galore failed");
    // Must exit cleanly with a clear error, not panic.
    assert!(
        !out.status.success(),
        "--clump_only --paired single-BAM (no --output-format ubam) must be rejected cleanly"
    );
    // Reject-message must route the user to --output-format ubam.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--output-format ubam"),
        "stderr should route the user to --output-format ubam; got: {stderr}"
    );
    // Panic marker check: no "panicked at" in stderr means we caught it
    // cleanly via bail! rather than falling through to the index-out-of-bounds.
    assert!(
        !stderr.contains("panicked at") && !stderr.contains("index out of bounds"),
        "regression: dispatch should not panic; got: {stderr}"
    );
}

#[test]
fn rejects_non_bam_n1_paired() {
    let dir = fresh_tmpdir("rej_non_bam_n1");
    let input = fixture("BS-seq_10K_R1.fastq.gz");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .output()
        .expect("trim_galore failed");
    assert!(
        !out.status.success(),
        "--paired --output-format ubam with N=1 non-BAM must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("interleaved uBAM") || stderr.contains("FASTQ"),
        "stderr should explain the format mismatch; got: {stderr}"
    );
}

#[test]
fn multi_pair_pe_bam_produces_one_output_per_pair() {
    // Regression guard for v1's N=4, 6, … multi-pair FASTQ input support.
    // Two pairs (four FASTQ inputs) → two output BAMs, one per pair.
    let dir = fresh_tmpdir("multi_pair");
    let r1 = fixture("BS-seq_10K_R1.fastq.gz");
    let r2 = fixture("BS-seq_10K_R2.fastq.gz");
    assert!(r1.exists() && r2.exists());

    // Copy the fixture pair to a `B_R{1,2}.fq.gz` alongside the first pair
    // so the two pairs have distinct output stems.
    let b_r1 = dir.join("B_R1.fastq.gz");
    let b_r2 = dir.join("B_R2.fastq.gz");
    std::fs::copy(&r1, &b_r1).unwrap();
    std::fs::copy(&r2, &b_r2).unwrap();

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&r1)
        .arg(&r2)
        .arg(&b_r1)
        .arg(&b_r2)
        .status()
        .expect("trim_galore failed");
    assert!(status.success(), "multi-pair PE BAM failed");

    let a_out = dir.join("BS-seq_10K_R1_clumped.bam");
    let b_out = dir.join("B_R1_clumped.bam");
    assert!(a_out.exists(), "pair-A output missing");
    assert!(b_out.exists(), "pair-B output missing");
    // Each output has 20K records (10K pairs × 2).
    assert_eq!(bam_tuples(&a_out).len(), 20_000);
    assert_eq!(bam_tuples(&b_out).len(), 20_000);
}

// ── Aux-tag round-trip (load-bearing v2 invariant) ────────────────

/// Regression guard for Reviewer B finding H-1: without this test, the
/// aux-tag round-trip claim in the plan / CHANGELOG / docs would be
/// completely untested. Fixture `ubam_test_with_tags.bam` carries `CB:Z`
/// and `UB:Z` tags on every record; running with `--preserve-tags CB,UB`
/// must preserve both tag values on the same records in the output.
#[test]
fn ubam_aux_tag_roundtrip_via_preserve_tags() {
    let dir = fresh_tmpdir("aux_tag_roundtrip");
    let input = fixture("ubam_test_with_tags.bam");
    assert!(input.exists(), "tagged uBAM fixture missing");

    // Collect input `(name, CB, UB)` tuples via samtools view — cheap
    // baseline that doesn't depend on noodles' internal aux representation.
    let baseline = std::process::Command::new("samtools")
        .args(["view", input.to_str().unwrap()])
        .output()
        .expect("samtools view failed on input fixture");
    assert!(baseline.status.success(), "samtools view failed");
    let baseline_stdout = String::from_utf8_lossy(&baseline.stdout);
    // Extract (name, CB:Z:..., UB:Z:...) per line.
    let extract_tags = |sam_line: &str| -> (String, Option<String>, Option<String>) {
        let mut fields = sam_line.split('\t');
        let name = fields.next().unwrap_or("").to_string();
        let mut cb = None;
        let mut ub = None;
        for f in fields {
            if let Some(rest) = f.strip_prefix("CB:Z:") {
                cb = Some(rest.to_string());
            } else if let Some(rest) = f.strip_prefix("UB:Z:") {
                ub = Some(rest.to_string());
            }
        }
        (name, cb, ub)
    };
    let mut input_tags: Vec<_> = baseline_stdout
        .lines()
        .filter(|l| !l.is_empty())
        .map(extract_tags)
        .collect();
    input_tags.sort();
    // Sanity: the fixture should have some CB tags. If this assertion
    // fails, the fixture changed and this test needs updating.
    assert!(
        input_tags.iter().any(|(_, cb, _)| cb.is_some()),
        "fixture doesn't carry CB tags — test scaffolding is wrong"
    );

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--output-format",
            "ubam",
            "--preserve-tags",
            "CB,UB",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()
        .expect("trim_galore failed");
    assert!(
        status.success(),
        "--clump_only --output-format ubam --preserve-tags CB,UB failed"
    );

    let out = dir.join("ubam_test_with_tags_clumped.bam");
    assert!(out.exists(), "output BAM missing: {}", out.display());

    // Extract output tags via samtools view.
    let out_view = std::process::Command::new("samtools")
        .args(["view", out.to_str().unwrap()])
        .output()
        .expect("samtools view failed on output");
    assert!(out_view.status.success());
    let out_stdout = String::from_utf8_lossy(&out_view.stdout);
    let mut output_tags: Vec<_> = out_stdout
        .lines()
        .filter(|l| !l.is_empty())
        .map(extract_tags)
        .collect();
    output_tags.sort();

    // Load-bearing invariant: after reorder + tag round-trip, the multiset
    // of (name, CB, UB) tuples must be identical.
    assert_eq!(
        input_tags.len(),
        output_tags.len(),
        "record count mismatch after aux-tag round-trip"
    );
    assert_eq!(
        input_tags, output_tags,
        "aux-tag multiset diverged — round-trip broken"
    );
    // Belt-and-braces: every output record MUST have CB and UB present.
    for (name, cb, ub) in &output_tags {
        assert!(
            cb.is_some() && ub.is_some(),
            "record {name} lost aux tags (CB={cb:?} UB={ub:?})"
        );
    }
}

// ── PE input-shape rejections (missing plan tests) ────────────────

/// Regression guard for plan-manager gap #2: the mixed-format Shape A
/// guard at `main.rs::dispatch` (rejecting FASTQ + BAM in the same pair)
/// had zero test coverage.
#[test]
fn rejects_mixed_format_paired() {
    let dir = fresh_tmpdir("rej_mixed_fmt");
    let r1_fastq = fixture("BS-seq_10K_R1.fastq.gz");
    let r2_bam = fixture("ubam_test.bam");
    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&r1_fastq)
        .arg(&r2_bam)
        .output()
        .expect("trim_galore failed");
    assert!(
        !out.status.success(),
        "--paired with mixed FASTQ+BAM formats must be rejected"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("same format") || stderr.contains("mixed"),
        "stderr should explain the mixed-format rejection; got: {stderr}"
    );
}

/// Regression guard for plan-manager gap #2 + Reviewer A finding C-2:
/// the `preflight_collision_bam` helper had zero test coverage. Two
/// input pairs whose case-folded output paths collide (via APFS/NTFS
/// case-insensitivity) should bail at pre-flight, before opening any
/// reader.
#[test]
fn pe_bam_collision_preflight_case_folded() {
    // Two pairs where the case-folded output stems collide: A/a differ
    // only in case → on a case-insensitive filesystem (or under case-
    // folded comparison), they map to the same output path.
    let dir = fresh_tmpdir("collision");
    let a_r1 = dir.join("A_R1.fq.gz");
    let a_r2 = dir.join("A_R2.fq.gz");
    let b_r1 = dir.join("a_R1.fq.gz"); // case-folded collision with A_R1
    let b_r2 = dir.join("a_R2.fq.gz");
    let src_r1 = fixture("BS-seq_10K_R1.fastq.gz");
    let src_r2 = fixture("BS-seq_10K_R2.fastq.gz");
    // On a case-insensitive filesystem (APFS default), the second copy
    // may resolve to the same inode as the first. Detect and skip in
    // that case — the collision-preflight logic still needs to work,
    // but this specific test needs a case-sensitive workspace.
    std::fs::copy(&src_r1, &a_r1).expect("copy A_R1");
    std::fs::copy(&src_r2, &a_r2).expect("copy A_R2");
    if a_r1.canonicalize().ok() == b_r1.canonicalize().ok().or(None) {
        // Try creating b_r1 — if it succeeds and is distinct, filesystem
        // is case-sensitive.
    }
    let write_result = std::fs::copy(&src_r1, &b_r1);
    if write_result.is_err() {
        // Case-insensitive FS: A_R1 and a_R1 are the same file, the
        // copy will succeed but they'll be the same inode. Skip this
        // test with a soft acknowledgement.
        eprintln!("Skipping pe_bam_collision_preflight_case_folded on case-insensitive FS");
        return;
    }
    std::fs::copy(&src_r2, &b_r2).expect("copy b_R2");
    // Distinct paths on disk? If not, skip.
    if std::fs::canonicalize(&a_r1).unwrap_or_default()
        == std::fs::canonicalize(&b_r1).unwrap_or_default()
    {
        eprintln!("Skipping pe_bam_collision_preflight_case_folded on case-insensitive FS");
        return;
    }

    let out = Command::new(binary())
        .args([
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&a_r1)
        .arg(&a_r2)
        .arg(&b_r1)
        .arg(&b_r2)
        .output()
        .expect("trim_galore failed");
    // Under case-folded collision, both pairs would want to write
    // `A_r1_clumped.bam` — should bail at preflight.
    assert!(
        !out.status.success(),
        "case-folded output-path collision must be caught at preflight"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("collision") || stderr.contains("would be written to the same file"),
        "stderr should mention collision; got: {stderr}"
    );
}

// ── --fastqc on BAM output ────────────────────────────────────────

#[test]
fn fastqc_produces_report_on_ubam_out() {
    // Q1 lock: fastqc-rust reads BAM natively; --fastqc should produce
    // an HTML + ZIP report on the reordered BAM.
    let dir = fresh_tmpdir("fastqc_bam");
    let input = fixture("BS-seq_10K_R1.fastq.gz");

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "--output-format",
            "ubam",
            "--fastqc",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
        ])
        .arg(&input)
        .status()
        .expect("trim_galore failed");
    assert!(
        status.success(),
        "--clump_only --output-format ubam --fastqc failed"
    );

    // FastQC report file naming: `<stem>_fastqc.html` / `.zip` where the
    // stem strips the .bam extension per fastqc-rust's extension stripper.
    let expected_html = dir.join("BS-seq_10K_R1_clumped_fastqc.html");
    let expected_zip = dir.join("BS-seq_10K_R1_clumped_fastqc.zip");
    assert!(
        expected_html.exists(),
        "FastQC HTML missing: {} (fastqc-rust reads BAM natively; if this fails, \
         the BAM input path may need investigation)",
        expected_html.display()
    );
    assert!(
        expected_zip.exists(),
        "FastQC ZIP missing: {}",
        expected_zip.display()
    );
}
