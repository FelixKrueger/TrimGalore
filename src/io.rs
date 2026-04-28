//! Output file naming conventions and I/O utilities.
//!
//! Implements TrimGalore-compatible output naming:
//! - Single-end: *_trimmed.fq(.gz)
//! - Paired-end: *_val_1.fq(.gz) / *_val_2.fq(.gz)
//! - Unpaired: *_unpaired_1.fq(.gz) / *_unpaired_2.fq(.gz)
//! - Reports: *_trimming_report.txt

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

/// True iff `path` ends with a `.gz` extension. The same heuristic that
/// `FastqReader` uses to decide whether to wrap the input in a gzip
/// decoder, so output naming + reader behaviour stay consistent.
///
/// Used to mirror input compression in the output filename / writer:
/// `plain.fastq` → `plain_trimmed.fq` (plain), `plain.fastq.gz` →
/// `plain_trimmed.fq.gz`. Matches Perl v0.6.x behaviour.
pub fn is_gzipped(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext == "gz")
}

/// Ensure the user-supplied `--output_dir` exists, creating it (and any
/// missing ancestors) if not. No-op when no `--output_dir` was passed or
/// the directory already exists.
///
/// Why this lives at the top of `main()` rather than at the per-file
/// writer: the parallel paired-end path opens its outputs via raw
/// `File::create` (see `parallel.rs`), which fails immediately on a
/// missing parent. By that point reader and worker threads have already
/// been spawned, so the early-return drops the receiver channel and the
/// process deadlocks with reader+workers stuck producing into a queue
/// nobody is draining. Hoisting the directory-create here covers every
/// downstream code path (parallel, single-threaded, paired, single-end,
/// every specialty mode) in one shot.
///
/// Matches Perl v0.6.x behaviour: "If an output directory which was
/// specified with -o output_directory did not exist, it will be created
/// for you" (v0.6.0 changelog).
pub fn ensure_output_dir(dir: Option<&Path>) -> Result<()> {
    let Some(dir) = dir else { return Ok(()) };
    if dir.exists() {
        return Ok(());
    }
    std::fs::create_dir_all(dir)
        .with_context(|| format!("Failed to create output directory: {}", dir.display()))?;
    Ok(())
}

/// Generate the trimmed output filename for single-end mode.
///
/// Follows TrimGalore convention:
/// - `.fastq.gz` → `_trimmed.fq.gz`
/// - `.fastq` → `_trimmed.fq`
/// - `.fq.gz` → `_trimmed.fq.gz`
/// - `.fq` → `_trimmed.fq`
pub fn single_end_output_name(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input));

    let ext = if gzip {
        "_trimmed.fq.gz"
    } else {
        "_trimmed.fq"
    };
    let filename = format!("{}{}", stem, ext);

    match output_dir {
        Some(dir) => dir.join(&filename),
        None => input.parent().unwrap_or(Path::new(".")).join(&filename),
    }
}

/// Generate the validated output filenames for paired-end mode.
///
/// Returns (val_1_path, val_2_path).
pub fn paired_end_output_names(
    input_r1: &Path,
    input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> (PathBuf, PathBuf) {
    let ext = if gzip { ".fq.gz" } else { ".fq" };

    let (stem1, stem2) = match basename {
        Some(b) => (format!("{}_R1", b), format!("{}_R2", b)),
        None => (
            strip_fastq_extensions(input_r1),
            strip_fastq_extensions(input_r2),
        ),
    };

    let f1 = format!("{}_val_1{}", stem1, ext);
    let f2 = format!("{}_val_2{}", stem2, ext);

    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf());

    (dir.join(&f1), dir.join(&f2))
}

/// Generate the unpaired output filenames for paired-end mode with --retain_unpaired.
pub fn unpaired_output_names(
    input_r1: &Path,
    input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> (PathBuf, PathBuf) {
    let ext = if gzip { ".fq.gz" } else { ".fq" };

    let (stem1, stem2) = match basename {
        Some(b) => (format!("{}_R1", b), format!("{}_R2", b)),
        None => (
            strip_fastq_extensions(input_r1),
            strip_fastq_extensions(input_r2),
        ),
    };

    let f1 = format!("{}_unpaired_1{}", stem1, ext);
    let f2 = format!("{}_unpaired_2{}", stem2, ext);

    let dir = output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf());

    (dir.join(&f1), dir.join(&f2))
}

/// Generate the trimming report filename.
pub fn report_name(input: &Path, output_dir: Option<&Path>) -> PathBuf {
    let input_name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    let report = format!("{}_trimming_report.txt", input_name);

    match output_dir {
        Some(dir) => dir.join(&report),
        None => input.parent().unwrap_or(Path::new(".")).join(&report),
    }
}

/// Generate the JSON trimming report filename.
pub fn json_report_name(input: &Path, output_dir: Option<&Path>) -> PathBuf {
    let input_name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    let report = format!("{}_trimming_report.json", input_name);

    match output_dir {
        Some(dir) => dir.join(&report),
        None => input.parent().unwrap_or(Path::new(".")).join(&report),
    }
}

/// Strip FASTQ extensions from a filename, returning just the base stem.
///
/// Handles: .fastq.gz, .fastq, .fq.gz, .fq
pub fn strip_fastq_extensions(path: &Path) -> String {
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    // Strip extensions in order of specificity
    for ext in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if name.ends_with(ext) {
            return name[..name.len() - ext.len()].to_string();
        }
    }

    // Fallback: strip .gz then whatever extension remains
    let name = if name.ends_with(".gz") {
        &name[..name.len() - 3]
    } else {
        &name
    };

    Path::new(name)
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_fastq_extensions() {
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fastq.gz")),
            "sample"
        );
        assert_eq!(strip_fastq_extensions(Path::new("sample.fq.gz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.fastq")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.fq")), "sample");
        assert_eq!(
            strip_fastq_extensions(Path::new("sample_R1.fq.gz")),
            "sample_R1"
        );
    }

    #[test]
    fn test_single_end_output_name() {
        let input = Path::new("/data/sample.fq.gz");
        let out = single_end_output_name(input, None, None, true);
        assert_eq!(out, PathBuf::from("/data/sample_trimmed.fq.gz"));

        let out = single_end_output_name(input, None, None, false);
        assert_eq!(out, PathBuf::from("/data/sample_trimmed.fq"));
    }

    #[test]
    fn test_single_end_with_basename() {
        let input = Path::new("/data/sample.fq.gz");
        let out = single_end_output_name(input, None, Some("custom"), true);
        assert_eq!(out, PathBuf::from("/data/custom_trimmed.fq.gz"));
    }

    #[test]
    fn test_paired_end_output_names() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let (o1, o2) = paired_end_output_names(r1, r2, None, None, true);
        assert_eq!(o1, PathBuf::from("/data/sample_R1_val_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/data/sample_R2_val_2.fq.gz"));
    }

    #[test]
    fn test_report_name() {
        let input = Path::new("/data/sample.fq.gz");
        let out = report_name(input, None);
        assert_eq!(out, PathBuf::from("/data/sample.fq.gz_trimming_report.txt"));
    }

    #[test]
    fn test_json_report_name() {
        let input = Path::new("/data/sample.fq.gz");
        let out = json_report_name(input, None);
        assert_eq!(
            out,
            PathBuf::from("/data/sample.fq.gz_trimming_report.json")
        );

        let out = json_report_name(input, Some(Path::new("/output")));
        assert_eq!(
            out,
            PathBuf::from("/output/sample.fq.gz_trimming_report.json")
        );
    }

    #[test]
    fn test_ensure_output_dir_creates_missing() {
        // Pick a path that does not exist; ensure_output_dir should create
        // it (along with any missing intermediate components). Regression
        // test for the parallel-path deadlock on missing --output_dir
        // (reported via beta.5 user feedback).
        let base = std::env::temp_dir().join("tg_ensure_dir_creates_missing");
        let _ = std::fs::remove_dir_all(&base);
        let nested = base.join("a").join("b").join("c");
        assert!(!nested.exists());

        ensure_output_dir(Some(&nested)).unwrap();
        assert!(nested.exists());
        assert!(nested.is_dir());

        std::fs::remove_dir_all(&base).unwrap();
    }

    #[test]
    fn test_ensure_output_dir_idempotent_on_existing() {
        // Calling on an already-existing directory must succeed silently
        // without altering the directory's mtime in a surprising way.
        let dir = std::env::temp_dir().join("tg_ensure_dir_idempotent");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        ensure_output_dir(Some(&dir)).unwrap();
        assert!(dir.exists());
        // Second call also succeeds.
        ensure_output_dir(Some(&dir)).unwrap();
        assert!(dir.exists());

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn test_ensure_output_dir_none_is_noop() {
        // No --output_dir was passed; helper must be a clean no-op.
        ensure_output_dir(None).unwrap();
    }

    #[test]
    fn test_is_gzipped() {
        // Common extensions that indicate gzip.
        assert!(is_gzipped(Path::new("sample.fastq.gz")));
        assert!(is_gzipped(Path::new("sample.fq.gz")));
        assert!(is_gzipped(Path::new("/some/dir/x.gz")));
        // Plain FASTQ — no gzip.
        assert!(!is_gzipped(Path::new("sample.fastq")));
        assert!(!is_gzipped(Path::new("sample.fq")));
        assert!(!is_gzipped(Path::new("/some/dir/x")));
        // Heuristic is extension-based (matches FastqReader's gzip detection),
        // so a misnamed file doesn't trigger.
        assert!(!is_gzipped(Path::new("sample.gz.fastq")));
    }
}
