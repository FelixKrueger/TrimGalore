//! Output file naming conventions and I/O utilities.
//!
//! Implements TrimGalore-compatible output naming:
//! - Single-end: *_trimmed.fq(.gz)
//! - Paired-end: *_val_1.fq(.gz) / *_val_2.fq(.gz)
//! - Unpaired: *_unpaired_1.fq(.gz) / *_unpaired_2.fq(.gz)
//! - Reports: *_trimming_report.txt

use std::path::{Path, PathBuf};

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
    let stem = basename.map(|b| b.to_string()).unwrap_or_else(|| {
        strip_fastq_extensions(input)
    });

    let ext = if gzip { "_trimmed.fq.gz" } else { "_trimmed.fq" };
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
        assert_eq!(strip_fastq_extensions(Path::new("sample.fastq.gz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.fq.gz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.fastq")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.fq")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample_R1.fq.gz")), "sample_R1");
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
        assert_eq!(out, PathBuf::from("/data/sample.fq.gz_trimming_report.json"));

        let out = json_report_name(input, Some(Path::new("/output")));
        assert_eq!(out, PathBuf::from("/output/sample.fq.gz_trimming_report.json"));
    }
}
