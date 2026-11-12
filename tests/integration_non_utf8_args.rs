//! A non-UTF-8 argument (legal in a Unix file name) must produce a normal
//! error, not a panic: `std::env::args()` panics on one, `args_os()` does not.

#![cfg(unix)]

use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;
use std::process::Command;

#[test]
fn non_utf8_file_name_gets_a_normal_error_not_a_panic() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .current_dir(dir.path())
        .arg(OsStr::from_bytes(b"\xff.fq"))
        .output()
        .expect("failed to run trim_galore");
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert_eq!(out.status.code(), Some(1), "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "stderr: {stderr}");
    assert!(stderr.contains("\u{FFFD}.fq"), "stderr: {stderr}");
}
