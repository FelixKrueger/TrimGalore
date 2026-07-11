//! Binary-driven integration test for the bare-invocation help path.
//!
//! Running `trim_galore` with no arguments should behave like most modern
//! CLIs: print the full help to **stdout** and exit **0**, rather than emit a
//! terse "the following required arguments were not provided" usage error to
//! stderr with a non-zero status. See `main.rs` (the `try_parse_from` handler
//! that routes clap's help/version outcomes to stdout + success) and
//! `cli.rs` (`arg_required_else_help = true`).
//!
//! `CARGO_BIN_EXE_trim_galore` is auto-set by cargo for integration tests.

use std::path::PathBuf;
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

#[test]
fn bare_invocation_prints_help_to_stdout_and_exits_zero() {
    let out = Command::new(binary())
        .output()
        .expect("failed to run trim_galore with no args");

    assert!(
        out.status.success(),
        "bare `trim_galore` should exit 0, got {:?}",
        out.status.code()
    );

    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);

    // Full help lands on stdout: usage line plus the options list.
    assert!(
        stdout.contains("Usage: trim_galore") && stdout.contains("--adapter"),
        "expected help text on stdout, got: {stdout}"
    );

    // No error message, and nothing on stderr for a bare invocation.
    assert!(
        !stdout.contains("error:") && stderr.is_empty(),
        "bare invocation should not emit an error; stderr was: {stderr:?}"
    );
}
