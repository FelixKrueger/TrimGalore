//! The paired-end text report's "shorter than the length cutoff" line must
//! count only length-cutoff removals. Pairs dropped for N content,
//! `--max_length` or `--discard_untrimmed` are reported on their own lines.

use std::path::Path;
use std::process::Command;

const ADAPTER: &str = "AGATCGGAAGAGC";
const SEQ: &str =
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACG";

/// Four pairs, each mate identical: one passes (30 bp after trimming), one is
/// too short (10 bp), one has no adapter (discarded as untrimmed), one is too
/// long (70 bp > `--max_length 60`).
fn write_fastq(path: &Path) {
    let reads = [
        format!("{}{}TT", &SEQ[..30], ADAPTER),
        format!("{}{}{}", &SEQ[..10], ADAPTER, &SEQ[30..57]),
        SEQ[..52].to_string(),
        format!("{}{}TT", &SEQ[..70], ADAPTER),
    ];
    let mut text = String::new();
    for (i, seq) in reads.iter().enumerate() {
        text.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(path, text).unwrap();
}

fn run(cores: &str) -> String {
    let dir = tempfile::tempdir().unwrap();
    let r1 = dir.path().join("s_R1.fq");
    let r2 = dir.path().join("s_R2.fq");
    write_fastq(&r1);
    write_fastq(&r2);
    let out = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args([
            "--paired",
            "-a",
            ADAPTER,
            "--no_poly_g",
            "--max_length",
            "60",
        ])
        .args(["--discard_untrimmed", "--cores", cores, "-o"])
        .arg(dir.path())
        .args([&r1, &r2])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::read_to_string(dir.path().join("s_R2.fq_trimming_report.txt")).unwrap()
}

#[test]
fn pe_report_shorter_line_counts_only_length_removals() {
    for cores in ["1", "2"] {
        let report = run(cores);
        assert!(
            report.contains(
                "Number of sequence pairs removed because at least one read was shorter than the length cutoff (25.00%):\n                                      1\n"
            ),
            "--cores {cores}:\n{report}"
        );
    }
}
