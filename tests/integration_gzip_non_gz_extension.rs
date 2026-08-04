//! Binary-driven integration tests for gzipped input whose filename does not
//! end in `.gz`.
//!
//! These exist because the unit tests did not catch the real gap. The library
//! tests exercise `FastqReader` directly, so they only cover the paths that
//! were already fixed; the callers that still derived gzip-ness from the
//! filename live in `main.rs` and `clump_only.rs` and are only reachable
//! through the binary. Nine green checks on the first version of this change
//! missed a `--paired --cores 1` regression for exactly that reason.
//!
//! The matrix below is therefore chosen by *dispatch path*, not by flag
//! aesthetics:
//!
//! | Test | Path it covers |
//! |---|---|
//! | single-end | `main::run_single`, sync reader |
//! | paired, `--cores 1` | `main::run_paired` serial, `FastqReader::open` |
//! | paired, `--cores 2` | `main::run_paired` parallel, threaded reader |
//! | `--clump_only` | `clump_only.rs`, its own reader construction |
//!
//! `CARGO_BIN_EXE_trim_galore` is auto-set by cargo for integration tests.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Locate the `trim_galore` binary built by cargo for this test target.
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// Create a fresh temp dir for the test, removing any leftover from a prior run.
fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Write `body` to `path`, gzip-compressed, under whatever name is asked for.
///
/// The point of every test here is that the name and the content disagree, so
/// the caller chooses them independently.
fn write_gz(path: &Path, body: &str) {
    let mut enc = flate2::write::GzEncoder::new(
        std::fs::File::create(path).unwrap(),
        flate2::Compression::default(),
    );
    enc.write_all(body.as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// A few reads: 40 bases of ordinary sequence followed by the start of the
/// Illumina adapter.
///
/// The sequence half has to be long enough to survive `--length` (20 by
/// default) once the adapter is cut, or every read is discarded and the output
/// is a valid but empty file, which would make these tests pass without
/// checking anything.
fn sample_reads(tag: &str) -> String {
    const BODY: &str = "ACGTTGCAACCGGTTAACGTACGTTGCAACCGGTTAACGT"; // 40 bp
    const ADAPTER: &str = "AGATCGGAAGAGC"; // Illumina, trimmed off
    let mut s = String::new();
    for i in 0..8 {
        let seq = format!("{BODY}{ADAPTER}");
        let qual = "I".repeat(seq.len());
        s.push_str(&format!("@{tag}_{i}\n{seq}\n+\n{qual}\n"));
    }
    s
}

/// Read a trimmed output back, whether or not it ended up compressed.
///
/// Which of the two it is depends on `io::is_gzipped`, which still looks at
/// the input *filename*, so a `.bgz` input currently yields plain output. That
/// asymmetry is deliberate and documented; these tests care that the reads are
/// correct, not how they are framed.
///
/// Note the stems the callers pass: a `.bgz` input is not recognised by
/// `io::strip_fastq_extensions`, so `pair_R1.fq.bgz` yields
/// `pair_R1.fq_val_1.fq`, with the inner `.fq` still in the name. That is a
/// known cosmetic wart (the report file disagrees with the output file, which
/// matters to MultiQC) and is filed for a separate change; these tests encode
/// current behaviour rather than the behaviour we would prefer.
fn read_output(dir: &Path, stem: &str) -> String {
    let text = read_output_raw(dir, stem);
    assert!(
        !text.trim().is_empty(),
        "output {stem} is empty: every read was filtered, so this test would \
         pass without checking anything"
    );
    text
}

/// The bytes, with no emptiness check. Split out so the assertion above has
/// something to call.
fn read_output_raw(dir: &Path, stem: &str) -> String {
    let plain = dir.join(format!("{stem}.fq"));
    let gz = dir.join(format!("{stem}.fq.gz"));
    if plain.exists() {
        std::fs::read_to_string(&plain).unwrap()
    } else {
        let f = std::fs::File::open(&gz).unwrap_or_else(|e| {
            panic!("no output at {} or {}: {e}", plain.display(), gz.display())
        });
        let mut out = String::new();
        std::io::Read::read_to_string(&mut flate2::read::MultiGzDecoder::new(f), &mut out).unwrap();
        out
    }
}

/// Single-end `.bgz`. The path that already worked before this change; kept so
/// a future regression in the sniff is caught here too.
#[test]
fn single_end_bgz_input_is_decompressed() {
    let dir = fresh_tmpdir("tg_bgz_se");
    let input = dir.join("sample.fq.bgz");
    write_gz(&input, &sample_reads("se"));

    let out = Command::new(binary())
        .args(["-o", dir.to_str().unwrap(), input.to_str().unwrap()])
        .output()
        .expect("binary must run");

    assert!(
        out.status.success(),
        "single-end .bgz must succeed. stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let trimmed = read_output(&dir, "sample.fq_trimmed");
    assert!(
        trimmed.contains("@se_0"),
        "trimmed output must contain the reads, got:\n{trimmed}"
    );
}

/// Paired-end at `--cores 1`, which is the default and was the regression.
///
/// Before the content sniff this failed, and worse than on base `dev`: the
/// entry guard passed, so the output files were created and truncated before
/// the first read failed. Asserting on file contents rather than only on the
/// exit status is what would catch that returning.
#[test]
fn paired_bgz_input_serial_is_decompressed() {
    let dir = fresh_tmpdir("tg_bgz_pe_serial");
    let r1 = dir.join("pair_R1.fq.bgz");
    let r2 = dir.join("pair_R2.fq.bgz");
    write_gz(&r1, &sample_reads("pe1"));
    write_gz(&r2, &sample_reads("pe2"));

    let out = Command::new(binary())
        .args([
            "--paired",
            "--cores",
            "1",
            "-o",
            dir.to_str().unwrap(),
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");

    assert!(
        out.status.success(),
        "paired .bgz at --cores 1 must succeed. stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(read_output(&dir, "pair_R1.fq_val_1").contains("@pe1_0"));
    assert!(read_output(&dir, "pair_R2.fq_val_2").contains("@pe2_0"));
}

/// Paired-end at `--cores 2`, which takes the worker-pool path and a different
/// reader construction. This one passed before the sniff, which is precisely
/// why it is worth pinning next to the serial case: the two must not diverge.
#[test]
fn paired_bgz_input_parallel_is_decompressed() {
    let dir = fresh_tmpdir("tg_bgz_pe_parallel");
    let r1 = dir.join("pair_R1.fq.bgz");
    let r2 = dir.join("pair_R2.fq.bgz");
    write_gz(&r1, &sample_reads("pp1"));
    write_gz(&r2, &sample_reads("pp2"));

    let out = Command::new(binary())
        .args([
            "--paired",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");

    assert!(
        out.status.success(),
        "paired .bgz at --cores 2 must succeed. stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(read_output(&dir, "pair_R1.fq_val_1").contains("@pp1_0"));
    assert!(read_output(&dir, "pair_R2.fq_val_2").contains("@pp2_0"));
}

/// Serial and parallel must agree on the reads, not merely both succeed.
#[test]
fn paired_bgz_serial_and_parallel_agree() {
    let dir = fresh_tmpdir("tg_bgz_pe_agree");
    let body1 = sample_reads("ag1");
    let body2 = sample_reads("ag2");

    let mut outputs = Vec::new();
    for cores in ["1", "2"] {
        let sub = dir.join(format!("cores{cores}"));
        std::fs::create_dir_all(&sub).unwrap();
        let r1 = sub.join("pair_R1.fq.bgz");
        let r2 = sub.join("pair_R2.fq.bgz");
        write_gz(&r1, &body1);
        write_gz(&r2, &body2);

        let out = Command::new(binary())
            .args([
                "--paired",
                "--cores",
                cores,
                "-o",
                sub.to_str().unwrap(),
                r1.to_str().unwrap(),
                r2.to_str().unwrap(),
            ])
            .output()
            .expect("binary must run");
        assert!(out.status.success(), "--cores {cores} must succeed");
        outputs.push(read_output(&sub, "pair_R1.fq_val_1"));
    }

    assert_eq!(
        outputs[0], outputs[1],
        "--cores 1 and --cores 2 must produce identical reads from .bgz input"
    );
}

/// `--clump_only` builds its own readers and was the third gap.
#[test]
fn clump_only_bgz_input_is_decompressed() {
    let dir = fresh_tmpdir("tg_bgz_clump");
    let input = dir.join("sample.fq.bgz");
    write_gz(&input, &sample_reads("cl"));

    let out = Command::new(binary())
        .args([
            "--clump_only",
            "-o",
            dir.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");

    assert!(
        out.status.success(),
        "--clump_only on .bgz must succeed. stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// The other direction, which this change also alters: a plain FASTQ misnamed
/// `.fastq.gz` previously failed with `invalid gzip header` and now reads.
///
/// Pinned as a test because it is a user-visible behaviour change, and because
/// nothing else in the suite would notice if it silently reverted.
#[test]
fn plain_fastq_misnamed_gz_is_read_as_plain() {
    let dir = fresh_tmpdir("tg_bgz_misnamed");
    let input = dir.join("sample.fastq.gz");
    std::fs::write(&input, sample_reads("mn")).unwrap();

    let out = Command::new(binary())
        .args(["-o", dir.to_str().unwrap(), input.to_str().unwrap()])
        .output()
        .expect("binary must run");

    assert!(
        out.status.success(),
        "plain FASTQ misnamed .gz must now be read. stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    // `sample.fastq.gz` IS in the strip list, so this one keeps the tidy stem.
    assert!(read_output(&dir, "sample_trimmed").contains("@mn_0"));
}
