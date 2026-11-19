//! JSON reports use a canonical executable name; BAM provenance keeps argv[0].

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::Command;

use noodles::bam;
use noodles::sam::header::record::value::map::program::tag as program_tag;

const BINARY: &str = env!("CARGO_BIN_EXE_trim_galore");
const FASTQ: &str = "@read1\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

fn run(command: &mut Command) {
    let output = command.output().expect("run trim_galore");
    assert!(
        output.status.success(),
        "{command:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn snapshot(
    dir: &Path,
    expected_report: &str,
    expected_raw: &str,
    expected_reports: usize,
) -> BTreeMap<String, Vec<u8>> {
    let mut files = BTreeMap::new();
    let mut reports = 0;
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_owned();
        let bytes = match path.extension().and_then(|ext| ext.to_str()) {
            Some("json") => {
                let bytes = fs::read(&path).unwrap();
                let json: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
                assert_eq!(json["command_line"], expected_report, "{name}");
                reports += 1;
                bytes
            }
            Some("bam") => {
                let mut reader =
                    bam::io::Reader::new(std::io::BufReader::new(fs::File::open(&path).unwrap()));
                let header = reader.read_header().unwrap();
                let program = header
                    .programs()
                    .as_ref()
                    .values()
                    .last()
                    .expect("appended program");
                let command_line: &[u8] = program
                    .other_fields()
                    .get(&program_tag::COMMAND_LINE)
                    .unwrap()
                    .as_ref();
                assert_eq!(command_line, expected_raw.as_bytes(), "{name}");
                let program_name: &[u8] = program
                    .other_fields()
                    .get(&program_tag::NAME)
                    .unwrap()
                    .as_ref();
                assert_eq!(program_name, b"trim_galore");
                let mut records = Vec::new();
                for result in reader.records() {
                    let record = result.unwrap();
                    records.push((
                        record
                            .name()
                            .map(|name| AsRef::<[u8]>::as_ref(name).to_vec()),
                        record.flags().bits(),
                        record.sequence().iter().collect::<Vec<_>>(),
                        record.quality_scores().as_ref().to_vec(),
                        format!("{:?}", record.data()),
                    ));
                }
                assert!(!records.is_empty(), "{name}");
                format!("{records:?}").into_bytes()
            }
            Some("gz") => {
                let mut bytes = Vec::new();
                flate2::read::MultiGzDecoder::new(fs::File::open(&path).unwrap())
                    .read_to_end(&mut bytes)
                    .unwrap();
                assert!(!bytes.is_empty(), "{name}");
                bytes
            }
            _ => fs::read(&path).unwrap(),
        };
        files.insert(name, bytes);
    }
    assert_eq!(reports, expected_reports);
    assert!(!files.is_empty());
    files
}

#[cfg(unix)]
fn compare_invocations(cwd: &Path, out: &Path, args: &[&str], reports: usize) {
    use std::os::unix::process::CommandExt;
    let expected_report = format!("trim_galore {}", args.join(" "));
    let mut baseline = None;
    for executable in [
        "trim_galore",
        "./trim_galore",
        "/opt/wave/bin/trim_galore",
        "/path with spaces/renamed-trimmer",
    ] {
        run(Command::new(BINARY)
            .current_dir(cwd)
            .arg0(executable)
            .args(args));
        let files = snapshot(
            out,
            &expected_report,
            &format!("{executable} {}", args.join(" ")),
            reports,
        );
        if let Some(expected) = &baseline {
            assert_eq!(&files, expected, "argv[0]={executable}");
        } else {
            baseline = Some(files);
        }
    }
}

#[cfg(unix)]
#[test]
fn reports_and_reads_are_independent_of_argv0_across_dispatches() {
    for paired in [false, true] {
        for bam_input in [false, true] {
            if bam_input && !paired {
                continue;
            }
            for bam_output in [false, true] {
                for cores in ["1", "2"] {
                    let dir = tempfile::tempdir().unwrap();
                    fs::write(dir.path().join("read_R1.fq"), FASTQ).unwrap();
                    fs::write(dir.path().join("read_R2.fq"), FASTQ).unwrap();
                    fs::copy(
                        Path::new(env!("CARGO_MANIFEST_DIR"))
                            .join("test_files/ubam_paired_test.bam"),
                        dir.path().join("paired.bam"),
                    )
                    .unwrap();
                    let mut args = vec![
                        "-a",
                        "AGATCGGAAGAGC",
                        "--no_poly_g",
                        "--cores",
                        cores,
                        "-o",
                        "out",
                    ];
                    if bam_output {
                        args.extend(["--output-format", "ubam"]);
                    } else {
                        args.push("--gzip");
                    }
                    if paired {
                        args.extend(["--paired", "-r1", "20", "-r2", "20"]);
                    }
                    if bam_input {
                        args.push("paired.bam");
                    } else {
                        args.push("read_R1.fq");
                        if paired {
                            args.push("read_R2.fq");
                        }
                    }
                    compare_invocations(
                        dir.path(),
                        &dir.path().join("out"),
                        &args,
                        if paired { 2 } else { 1 },
                    );
                }
            }
        }
    }
}

#[cfg(unix)]
#[test]
fn specialty_bam_headers_preserve_original_invocation() {
    for mode in ["--hardtrim5", "--hardtrim3", "--clump_only"] {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("reads.fq"), FASTQ).unwrap();
        let mut args = vec![mode];
        if mode == "--clump_only" {
            args.extend(["--cores", "2", "--memory", "320M", "--no_report_file"]);
        } else {
            args.push("20");
        }
        args.extend(["--output-format", "ubam", "-o", "out", "reads.fq"]);
        compare_invocations(dir.path(), &dir.path().join("out"), &args, 0);
    }
}

#[cfg(unix)]
#[test]
fn appended_program_keeps_raw_invocation_with_existing_history() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("reads.fq"), FASTQ).unwrap();
    run(Command::new(BINARY).current_dir(dir.path()).args([
        "--output-format",
        "ubam",
        "-a",
        "AGATCGGAAGAGC",
        "--no_poly_g",
        "-o",
        "source",
        "reads.fq",
    ]));
    compare_invocations(
        dir.path(),
        &dir.path().join("out"),
        &[
            "--output-format",
            "ubam",
            "-a",
            "AGATCGGAAGAGC",
            "--no_poly_g",
            "-o",
            "out",
            "source/reads_trimmed.bam",
        ],
        1,
    );
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(
        fs::File::open(dir.path().join("out/reads_trimmed_trimmed.bam")).unwrap(),
    ));
    let header = reader.read_header().unwrap();
    assert_eq!(header.programs().as_ref().len(), 2);
}

#[test]
fn relative_and_absolute_executable_paths_produce_identical_reports() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fq");
    let out = dir.path().join("out");
    fs::write(&input, FASTQ).unwrap();
    let args = [
        "-a",
        "AGATCGGAAGAGC",
        "--no_poly_g",
        "--gzip",
        "-o",
        out.to_str().unwrap(),
        input.to_str().unwrap(),
    ];
    let absolute = Path::new(BINARY);
    let relative = Path::new(".").join(absolute.file_name().unwrap());
    let expected_report = format!("trim_galore {}", args.join(" "));
    let mut baseline = None;
    for executable in [absolute, relative.as_path()] {
        run(Command::new(executable)
            .current_dir(absolute.parent().unwrap())
            .args(args));
        let files = snapshot(&out, &expected_report, "unused for FASTQ output", 1);
        if let Some(expected) = &baseline {
            assert_eq!(&files, expected);
        } else {
            baseline = Some(files);
        }
    }
}
