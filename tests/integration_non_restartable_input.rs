//! Binary-driven integration tests for issue #379 — input that cannot be
//! re-read from the start.
//!
//! Trim Galore opens each input several times (format detection, sanity check,
//! adapter auto-detection, the trim pass), so a stream cannot be processed.
//! Before the fix these inputs either blamed the user's data for the format
//! (`doesn't seem to be in FastQ format`), blamed the compression (`Failed to
//! decompress first block`), or — for a FIFO — hung with no message at all.
//!
//! Every child here is bounded and killed rather than waited on:
//! `Command::output()` reads both pipes to EOF, so a child blocked in `open()`
//! would never return and the test would look like a slow job, not a failure.

#![cfg(unix)]

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const REJECTION: &str = "cannot be re-read from the start";

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// One directory per test: a shared slug would let one nextest process unlink
/// another's FIFO.
fn fresh_dir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_fastq(path: &Path) {
    let mut body = String::new();
    for i in 0..8 {
        body.push_str(&format!(
            "@read{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    std::fs::write(path, body).unwrap();
}

fn make_fifo(path: &Path) {
    let st = Command::new("mkfifo")
        .arg(path)
        .status()
        .expect("mkfifo must run");
    assert!(st.success(), "mkfifo {} failed", path.display());
    use std::os::unix::fs::FileTypeExt;
    assert!(
        std::fs::metadata(path).unwrap().file_type().is_fifo(),
        "fixture must be a FIFO; a regular file would test nothing"
    );
}

/// Run `cmd` with a hard time bound. Returns `(succeeded, stderr)`.
fn run_bounded(mut cmd: Command, what: &str) -> (bool, String) {
    let mut child = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("failed to spawn {what}: {e}"));

    // Drain stderr concurrently so a child writing more than the pipe buffer
    // cannot deadlock against our wait loop.
    let mut pipe = child.stderr.take().expect("stderr is piped");
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let mut s = String::new();
        let _ = pipe.read_to_string(&mut s);
        let _ = tx.send(s);
    });

    let deadline = Instant::now() + Duration::from_secs(20);
    let status = loop {
        match child.try_wait().expect("try_wait failed") {
            Some(st) => break st,
            None if Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                // Drain before panicking: on the failure that matters most, the
                // binary's own output is the only useful diagnostic.
                let stderr = rx.recv_timeout(Duration::from_secs(5)).unwrap_or_default();
                panic!("{what} never exited; it must fail, not hang. stderr was: {stderr}");
            }
            None => std::thread::sleep(Duration::from_millis(25)),
        }
    };
    let stderr = rx.recv_timeout(Duration::from_secs(5)).unwrap_or_default();
    (status.success(), stderr)
}

/// A FIFO with no writer at all. `File::open` on this blocks forever, so before
/// the fix this hung; the path `stat` in `Cli::validate` is what makes it fail.
#[test]
fn fifo_input_is_rejected_and_creates_no_output() {
    let dir = fresh_dir("tg_379_fifo_in");
    let fifo = dir.join("stream.fq");
    make_fifo(&fifo);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("-o").arg(&out).arg(&fifo);
    let (ok, stderr) = run_bounded(cmd, "trim_galore on a FIFO");

    assert!(!ok, "a FIFO must be rejected; stderr was: {stderr}");
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
    assert!(
        stderr.contains("is a pipe or FIFO"),
        "unexpected message: {stderr}"
    );
    assert!(
        !out.exists(),
        "rejection must happen before any output directory is created"
    );
}

/// Plain FASTQ down a pipe. The old failure blamed the data's format.
#[test]
fn dev_stdin_from_a_plain_pipe_is_rejected() {
    let dir = fresh_dir("tg_379_stdin_plain");
    let fq = dir.join("reads.fastq");
    write_fastq(&fq);
    let out = dir.join("out");

    let mut cmd = Command::new("sh");
    cmd.arg("-c").arg(format!(
        "cat '{}' | '{}' -o '{}' /dev/stdin",
        fq.display(),
        binary().display(),
        out.display()
    ));
    let (ok, stderr) = run_bounded(cmd, "trim_galore on /dev/stdin from a pipe");

    assert!(!ok, "a pipe must be rejected; stderr was: {stderr}");
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
    // Negative control for this row specifically.
    assert!(
        !stderr.contains("doesn't seem to be in FastQ format"),
        "must not blame the data's format: {stderr}"
    );
}

/// Gzipped FASTQ down a pipe. This row's old failure was a *different* wrong
/// message — the gzip branch re-opens and `MultiGzDecoder` chokes on the
/// consumed prefix — so it needs its own negative control.
#[test]
fn dev_stdin_from_a_gzipped_pipe_is_rejected() {
    let dir = fresh_dir("tg_379_stdin_gz");
    let fq = dir.join("reads.fastq");
    write_fastq(&fq);
    let out = dir.join("out");

    let mut cmd = Command::new("sh");
    cmd.arg("-c").arg(format!(
        "gzip -c '{}' | '{}' -o '{}' /dev/stdin",
        fq.display(),
        binary().display(),
        out.display()
    ));
    let (ok, stderr) = run_bounded(cmd, "trim_galore on a gzipped pipe");

    assert!(!ok, "a gzipped pipe must be rejected; stderr was: {stderr}");
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
    assert!(
        !stderr.contains("Failed to decompress first block"),
        "must not blame the compression: {stderr}"
    );
}

/// `--passthrough` is a separate CLI surface from `cli.input`, and it reached
/// none of the format-detection guards. Both plan reviewers found this
/// independently; before the fix a FIFO here hung *after* the R1/R2 writers had
/// been created, leaving partial output on disk.
#[test]
fn passthrough_fifo_is_rejected_and_creates_no_output() {
    let dir = fresh_dir("tg_379_passthrough");
    let r1 = dir.join("r1.fastq");
    let r2 = dir.join("r2.fastq");
    write_fastq(&r1);
    write_fastq(&r2);
    let fifo = dir.join("index.fq");
    make_fifo(&fifo);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("--paired")
        .arg("-o")
        .arg(&out)
        .arg("--passthrough")
        .arg(&fifo)
        .arg(&r1)
        .arg(&r2);
    let (ok, stderr) = run_bounded(cmd, "trim_galore --passthrough on a FIFO");

    assert!(
        !ok,
        "a FIFO passthrough must be rejected; stderr was: {stderr}"
    );
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
    assert!(
        !out.exists(),
        "rejection must precede output creation, or a hang leaves partial output"
    );
}

/// A uBAM passthrough must be rejected before any output exists. The existing
/// uBAM rejection computes `any_bam` from `cli.input` only, so routing the
/// passthrough file through `sanity_check_any` made it *pass* — the BAM arm
/// accepts it — and the run then died in the FASTQ reader with three empty
/// output files on disk.
#[test]
fn passthrough_ubam_is_rejected_before_output() {
    let dir = fresh_dir("tg_379_passthrough_ubam");
    let r1 = dir.join("r1.fastq");
    let r2 = dir.join("r2.fastq");
    write_fastq(&r1);
    write_fastq(&r2);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("--paired")
        .arg("-o")
        .arg(&out)
        .arg("--passthrough")
        .arg("test_files/ubam_test.bam")
        .arg(&r1)
        .arg(&r2);
    let (ok, stderr) = run_bounded(cmd, "trim_galore --passthrough on a uBAM");

    assert!(!ok, "a uBAM passthrough must be rejected; stderr: {stderr}");
    assert!(
        stderr.contains("not supported with uBAM input"),
        "unexpected message: {stderr}"
    );
    assert!(
        !stderr.contains("valid UTF-8"),
        "must not fail deep in the FASTQ reader: {stderr}"
    );
    assert!(
        !out.join("r1_val_1.fq").exists(),
        "must reject before writing output"
    );
}

/// The class layer 1 cannot see: `/dev/fd/N` over a regular file `stat`s as a
/// regular file, so only the second-open check rejects it. macOS-only because
/// Linux resolves `/dev/stdin` through `/proc` and genuinely restarts there.
#[cfg(target_os = "macos")]
#[test]
fn passthrough_dev_stdin_is_rejected() {
    let dir = fresh_dir("tg_379_passthrough_devfd");
    let r1 = dir.join("r1.fastq");
    let r2 = dir.join("r2.fastq");
    let idx = dir.join("idx.fastq");
    write_fastq(&r1);
    write_fastq(&r2);
    write_fastq(&idx);
    let out = dir.join("out");

    let mut cmd = Command::new("sh");
    cmd.arg("-c").arg(format!(
        "'{}' --paired -o '{}' --passthrough /dev/stdin '{}' '{}' < '{}'",
        binary().display(),
        out.display(),
        r1.display(),
        r2.display(),
        idx.display()
    ));
    let (ok, stderr) = run_bounded(cmd, "trim_galore --passthrough /dev/stdin");

    assert!(!ok, "stderr was: {stderr}");
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
}

/// The barcode file is read once, so restartability does not apply — but a
/// writer-less FIFO blocked in `File::open` after trimming and reporting had
/// already finished.
#[test]
fn demux_fifo_barcode_file_is_rejected() {
    let dir = fresh_dir("tg_379_demux");
    let fq = dir.join("reads.fastq");
    write_fastq(&fq);
    let fifo = dir.join("barcodes.txt");
    make_fifo(&fifo);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("-o").arg(&out).arg("--demux").arg(&fifo).arg(&fq);
    let (ok, stderr) = run_bounded(cmd, "trim_galore --demux on a FIFO");

    assert!(
        !ok,
        "a FIFO barcode file must be rejected; stderr: {stderr}"
    );
    assert!(stderr.contains(REJECTION), "unexpected message: {stderr}");
    assert!(
        !out.exists(),
        "rejection must precede trimming, not follow it"
    );
}

/// Positive control: the guard must not reject an ordinary regular file. Without
/// this, every assertion above would pass on a build that rejected everything.
#[test]
fn regular_file_still_runs() {
    let dir = fresh_dir("tg_379_control");
    let fq = dir.join("reads.fastq");
    write_fastq(&fq);
    let out = dir.join("out");

    let mut cmd = Command::new(binary());
    cmd.arg("-o").arg(&out).arg("--dont_gzip").arg(&fq);
    let (ok, stderr) = run_bounded(cmd, "trim_galore on a regular file");

    assert!(ok, "a regular file must still trim; stderr was: {stderr}");
    assert!(
        !stderr.contains(REJECTION),
        "regular file wrongly rejected: {stderr}"
    );
}
