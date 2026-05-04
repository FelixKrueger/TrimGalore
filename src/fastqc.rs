//! FastQC integration via the bundled `fastqc-rust` library.
//!
//! Replaces the v0.6.x-style shell-out to an external `fastqc` binary
//! with an in-process call. Targets byte-equivalent output to Java
//! FastQC 0.12.1 — the same version we previously bundled in the Docker
//! image — see <https://github.com/ewels/FastQC-Rust> for the upstream
//! crate. The bundled approach removes Java + the FastQC tarball as
//! runtime dependencies, completing the "single static binary, zero
//! deps" story for the v2.x rewrite.
//!
//! Public entry point: [`run`]. The two callers in `src/main.rs` (one
//! per single-end output, two per paired-end output) invoke it with the
//! trimmed FASTQ path; output `*_fastqc.html` and `*_fastqc.zip` land
//! in the trim-galore `--output_dir` (or current dir if unset).

use std::path::Path;

use anyhow::Result;
use fastqc_rust::{config::FastQCConfig, runner};

/// Run FastQC on a single trimmed output file.
///
/// `output_path` — the FASTQ to analyse (one of the trimmed outputs).
/// `fastqc_args` — raw user-supplied `--fastqc_args` string. A subset
///                 of common FastQC CLI flags is translated to
///                 [`FastQCConfig`] mutations; unknown flags emit a
///                 warning and are ignored.
/// `output_dir`  — where the `*_fastqc.html` / `*_fastqc.zip` artifacts
///                 land. `None` means current directory (FastQC's
///                 default).
/// `cores`       — threading budget for `fastqc-rust`'s internal rayon
///                 pool. Always at least 1.
pub fn run(
    output_path: &Path,
    fastqc_args: Option<&str>,
    output_dir: Option<&Path>,
    cores: usize,
) -> Result<()> {
    let mut config = FastQCConfig {
        output_dir: output_dir.map(|p| p.to_path_buf()),
        threads: cores.max(1),
        // trim-galore prints "Running FastQC on …" itself; suppress
        // fastqc-rust's per-file progress chatter to avoid
        // double-logging.
        quiet: true,
        ..FastQCConfig::default()
    };

    if let Some(args) = fastqc_args {
        apply_fastqc_args(&mut config, args);
    }

    eprintln!("\nRunning FastQC on {}", output_path.display());
    runner::run(&config, &[output_path.to_path_buf()])
        .map_err(|code| anyhow::anyhow!("FastQC analysis failed (exit code {})", code))
}

/// Translate a subset of Java FastQC CLI flags into [`FastQCConfig`]
/// mutations.
///
/// Supports the common flags users pass via `--fastqc_args`. Unknown
/// flags emit a warning and are ignored — deliberate forward-compat:
/// `fastqc-rust` may add flags we don't list yet, and we don't want
/// every new upstream flag to require a trim-galore release.
///
/// Flags handled:
///   `--nogroup`, `--expgroup`, `--quiet`, `--svg`, `--nano`,
///   `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`.
fn apply_fastqc_args(config: &mut FastQCConfig, raw: &str) {
    let tokens: Vec<&str> = raw.split_whitespace().collect();
    let mut i = 0;
    while i < tokens.len() {
        match tokens[i] {
            "--nogroup" => config.nogroup = true,
            "--expgroup" => config.expgroup = true,
            "--quiet" => config.quiet = true,
            "--svg" => config.svg_output = true,
            "--nano" => config.nano = true,
            "--nofilter" => config.nofilter = true,
            "--casava" => config.casava = true,
            "-t" | "--threads" => {
                if let Some(v) = tokens.get(i + 1).and_then(|t| t.parse().ok()) {
                    config.threads = v;
                    i += 1;
                }
            }
            "-o" | "--outdir" => {
                if let Some(v) = tokens.get(i + 1) {
                    config.output_dir = Some(std::path::PathBuf::from(v));
                    i += 1;
                }
            }
            unknown => {
                eprintln!(
                    "Warning: --fastqc_args flag '{}' is not translated to the bundled FastQC; ignored",
                    unknown
                );
            }
        }
        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &str) -> FastQCConfig {
        let mut config = FastQCConfig::default();
        apply_fastqc_args(&mut config, args);
        config
    }

    #[test]
    fn test_apply_fastqc_args_nogroup() {
        let config = parse("--nogroup");
        assert!(config.nogroup);
    }

    #[test]
    fn test_apply_fastqc_args_threads() {
        let config = parse("-t 4");
        assert_eq!(config.threads, 4);
    }

    #[test]
    fn test_apply_fastqc_args_threads_long_form() {
        let config = parse("--threads 8");
        assert_eq!(config.threads, 8);
    }

    #[test]
    fn test_apply_fastqc_args_outdir() {
        let config = parse("-o /tmp/somewhere");
        assert_eq!(
            config.output_dir,
            Some(std::path::PathBuf::from("/tmp/somewhere"))
        );
    }

    #[test]
    fn test_apply_fastqc_args_combined() {
        let config = parse("--nogroup --svg -t 2");
        assert!(config.nogroup);
        assert!(config.svg_output);
        assert_eq!(config.threads, 2);
    }

    #[test]
    fn test_apply_fastqc_args_unknown_flag_ignored() {
        // --quiet + an unknown flag: the unknown one warns and is
        // skipped; the known one still applies.
        let config = parse("--unknown-flag --quiet");
        assert!(config.quiet);
    }

    #[test]
    fn test_apply_fastqc_args_threads_without_value() {
        // `-t` with no following token: leave threads unchanged.
        let config = parse("-t");
        assert_eq!(config.threads, FastQCConfig::default().threads);
    }

    #[test]
    fn test_apply_fastqc_args_threads_non_numeric() {
        // `-t abc` — non-numeric value; leave threads unchanged.
        // (The token "abc" then becomes a positional in the next loop
        // iteration and falls through to the unknown-flag branch.)
        let config = parse("-t abc");
        assert_eq!(config.threads, FastQCConfig::default().threads);
    }
}
