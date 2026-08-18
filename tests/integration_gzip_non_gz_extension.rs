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
/// The stems the callers pass are the tidy ones: since #381
/// `io::strip_fastq_extensions` recognises the `.bgz` forms, so
/// `pair_R1.fq.bgz` yields `pair_R1_val_1.fq` with no inner `.fq` left in the
/// name. `bgz_output_and_report_agree_on_the_stem` below pins that directly.
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
/// #453 — `--clumpify` accepts `--dont_gzip`, as `--clump_only` already did. The
/// refusal it replaced called clumping plain text pointless while the tool did
/// exactly that on plain input, so the reason was never true.
#[test]
fn clumpify_accepts_dont_gzip_and_writes_plain_output() {
    let dir = fresh_tmpdir("tg_453_clumpify_plain");
    let input = dir.join("sample.fq.gz");
    write_gz(&input, &sample_reads("cy"));

    let out = Command::new(binary())
        .args([
            "--clumpify",
            "--dont_gzip",
            "--cores",
            "2",
            "-o",
            dir.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(
        out.status.success(),
        "--clumpify --dont_gzip must be accepted:\n{err}"
    );
    assert!(
        !err.contains("mutually exclusive"),
        "the refusal must be gone:\n{err}"
    );

    // Plain output, and every record still present: --dont_gzip changes the
    // container, not the contents.
    let plain = dir.join("sample_trimmed.fq");
    assert!(
        plain.exists(),
        "expected plain output, got: {:?}",
        std::fs::read_dir(&dir).map(|d| d
            .filter_map(|e| e.ok())
            .map(|e| e.file_name())
            .collect::<Vec<_>>())
    );
    let body = std::fs::read_to_string(&plain).expect("read plain output");
    assert_eq!(
        body.lines().count() / 4,
        8,
        "all eight records must survive:\n{body}"
    );
}

/// #453 — output compression follows the input, so `--gzip`'s deprecation notice
/// must not promise gzipped output for a plain-text input, where no flag delivers it.
#[test]
fn gzip_deprecation_notice_matches_what_the_run_will_do() {
    let plain_claim = "produces plain-text output and no flag changes that";
    let gz_claim = "Output is gzipped by default";

    // Plain input: the notice must describe follow-the-input.
    let dir = fresh_tmpdir("tg_453_plain");
    let input = dir.join("sample.fastq");
    std::fs::write(&input, sample_reads("se")).expect("write plain input");
    let out = Command::new(binary())
        .args([
            "--gzip",
            "-o",
            dir.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(out.status.success(), "plain run must succeed:\n{err}");
    assert!(
        err.contains(plain_claim),
        "plain input must not be told its output is gzipped:\n{err}"
    );
    assert!(
        !err.contains(gz_claim),
        "the gzipped-by-default claim is false here:\n{err}"
    );
    assert!(
        dir.join("sample_trimmed.fq").exists(),
        "plain input gives plain output, which is what the notice now says: {:?}",
        std::fs::read_dir(&dir).map(|d| d.count())
    );

    // Gzipped input: the original wording is correct and must be kept.
    let dir = fresh_tmpdir("tg_453_gz");
    let input = dir.join("sample.fq.gz");
    write_gz(&input, &sample_reads("se"));
    let out = Command::new(binary())
        .args([
            "--gzip",
            "-o",
            dir.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .expect("binary must run");
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(out.status.success(), "gzipped run must succeed:\n{err}");
    assert!(
        err.contains(gz_claim),
        "gzipped input keeps the original wording:\n{err}"
    );
    assert!(
        !err.contains(plain_claim),
        "the follow-the-input wording belongs to plain input only:\n{err}"
    );
}

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
    let trimmed = read_output(&dir, "sample_trimmed");
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
    assert!(read_output(&dir, "pair_R1_val_1").contains("@pe1_0"));
    assert!(read_output(&dir, "pair_R2_val_2").contains("@pe2_0"));
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
    assert!(read_output(&dir, "pair_R1_val_1").contains("@pp1_0"));
    assert!(read_output(&dir, "pair_R2_val_2").contains("@pp2_0"));
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
        outputs.push(read_output(&sub, "pair_R1_val_1"));
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

/// The reproduction from
/// [#381](https://github.com/FelixKrueger/TrimGalore/issues/381): trimmed
/// output and trimming report must agree on the sample name.
///
/// The report filename is the input filename plus a suffix, while the output
/// filename is the *stripped* stem, so the two only line up when the stripper
/// knows the input's extension. It did not know `.bgz`, so `sample.fastq.bgz`
/// gave `sample.fastq_trimmed.fq` next to
/// `sample.fastq.bgz_trimming_report.txt`. MultiQC takes the sample name from
/// the report filename, so the pair can land under different samples.
///
/// The same bytes are run twice, once under each name, and the assertion is
/// the same for both: this is a contrast test, not a hard-coded expectation.
#[test]
fn bgz_output_and_report_agree_on_the_stem() {
    for (name, stem) in [
        ("sample.fastq.bgz", "sample"),
        ("sample.fq.bgz", "sample"),
        ("sample.fastq.gz", "sample"), // the control: always worked
    ] {
        let dir = fresh_tmpdir(&format!("tg_bgz_stem_{stem}_{}", name.replace('.', "_")));
        let input = dir.join(name);
        write_gz(&input, &sample_reads("st"));

        let out = Command::new(binary())
            .args([
                "--dont_gzip",
                "-o",
                dir.to_str().unwrap(),
                input.to_str().unwrap(),
            ])
            .output()
            .expect("binary must run");
        assert!(
            out.status.success(),
            "{name} must trim. stderr:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );

        let trimmed = dir.join(format!("{stem}_trimmed.fq"));
        assert!(
            trimmed.exists(),
            "{name}: expected output {}, found {:?}",
            trimmed.display(),
            std::fs::read_dir(&dir)
                .unwrap()
                .map(|e| e.unwrap().file_name())
                .collect::<Vec<_>>()
        );
        for report in [
            dir.join(format!("{name}_trimming_report.txt")),
            dir.join(format!("{name}_trimming_report.json")),
        ] {
            assert!(
                report.exists(),
                "{name}: expected report {}",
                report.display()
            );
        }
    }
}

// ── #384: uppercase extensions behave as their lowercase equivalents ──────

/// V5. Literal filename assertions, computed without calling any function under
/// test; no `--dont_gzip`, so both halves of the change are observable. The two
/// halves fail independently: a wrong stem fails at path resolution, a wrong
/// compression decision fails on the magic bytes.
#[test]
fn uppercase_extension_names_and_compresses_like_lowercase() {
    let dir = fresh_tmpdir("tg_384_upper");
    let input = dir.join("SAMPLE.FASTQ.GZ");
    write_gz(&input, &sample_reads("U"));
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();

    let status = Command::new(binary())
        .args(["-o", out.to_str().unwrap(), input.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());

    // Naming half: exact literals. Before #384 the stem kept `.FASTQ` and the
    // output was plain (`SAMPLE.FASTQ_trimmed.fq`).
    let trimmed = out.join("SAMPLE_trimmed.fq.gz");
    assert!(
        trimmed.is_file(),
        "expected SAMPLE_trimmed.fq.gz, found: {:?}",
        std::fs::read_dir(&out)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.file_name())
            .collect::<Vec<_>>()
    );
    assert!(out.join("SAMPLE.FASTQ.GZ_trimming_report.txt").is_file());

    // Compression half: gzip magic, checked on bytes — the name is what is under test.
    let bytes = std::fs::read(&trimmed).unwrap();
    assert!(
        bytes.starts_with(&[0x1f, 0x8b]),
        "output must be gzip-compressed, got {:02x?}",
        &bytes[..bytes.len().min(4)]
    );
}

/// V4. The refusals the fold newly creates: mixed extension *spellings* whose
/// stems now agree. (Pure case-variants are already refused today via the
/// case-folded report keys, so they cannot discriminate — see the #384 plan.)
#[test]
fn uppercase_and_lowercase_spellings_of_one_stem_now_collide() {
    let dir = fresh_tmpdir("tg_384_collide");
    // On case-insensitive filesystems (default APFS) the last two names alias ONE
    // file, so the second write_gz overwrites the first and the loop's second
    // iteration feeds the same inode under a different spelling. Both iterations
    // still refuse on both filesystem kinds; do not add content assertions here.
    write_gz(&dir.join("SAMPLE.FASTQ.GZ"), &sample_reads("A"));
    write_gz(&dir.join("sample.fq.gz"), &sample_reads("B"));
    write_gz(&dir.join("SAMPLE.FQ.GZ"), &sample_reads("C"));
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();

    for second in ["sample.fq.gz", "SAMPLE.FQ.GZ"] {
        let output = Command::new(binary())
            .current_dir(&dir)
            .args(["-o", "out", "SAMPLE.FASTQ.GZ", second])
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "{second}: expected collision refusal"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("Output path collision"),
            "{second}: {stderr}"
        );
        assert!(
            std::fs::read_dir(&out).unwrap().next().is_none(),
            "{second}: a refused run must write nothing"
        );
    }
}

/// V8. The least-guarded consumer: `--clump_only`'s report derives its
/// `Input: … (gzip|plain)` label and its `Compression ratio:` gate from
/// `is_gzipped`. Before #384 a `.GZ`-named input read `plain` with no ratio line.
#[test]
fn clump_only_uppercase_gz_reports_gzip_and_ratio() {
    let dir = fresh_tmpdir("tg_384_clump");
    let input = dir.join("SAMPLE.FASTQ.GZ");
    write_gz(&input, &sample_reads("K"));
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();

    let status = Command::new(binary())
        .args([
            "--clump_only",
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    assert!(
        out.join("SAMPLE_clumped.fq.gz").is_file(),
        "clump stem must fold too"
    );
    let report = std::fs::read_to_string(out.join("SAMPLE.FASTQ.GZ_clumping_report.txt")).unwrap();
    assert!(
        report.contains("gzip"),
        "Input label must say gzip:\n{report}"
    );
    assert!(!report.contains("plain"), "must not say plain:\n{report}");
    assert!(
        report.contains("Compression ratio:"),
        "ratio line must appear:\n{report}"
    );
}

/// #384 on the paired path, restoring this file's dispatch-path matrix: paired
/// mode adds the coupling that `input[0]` drives the run-wide `gzip` flag.
#[test]
fn paired_uppercase_extensions_name_and_compress_like_lowercase() {
    let dir = fresh_tmpdir("tg_384_paired");
    write_gz(&dir.join("P_R1.FASTQ.GZ"), &sample_reads("R1"));
    write_gz(&dir.join("P_R2.FASTQ.GZ"), &sample_reads("R2"));
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();

    let status = Command::new(binary())
        .current_dir(&dir)
        .args(["--paired", "-o", "out", "P_R1.FASTQ.GZ", "P_R2.FASTQ.GZ"])
        .status()
        .unwrap();
    assert!(status.success());

    for stem in ["P_R1_val_1", "P_R2_val_2"] {
        let f = out.join(format!("{stem}.fq.gz"));
        assert!(f.is_file(), "missing {stem}.fq.gz");
        let bytes = std::fs::read(&f).unwrap();
        assert!(bytes.starts_with(&[0x1f, 0x8b]), "{stem} must be gzip");
    }
}
