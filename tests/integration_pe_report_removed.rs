//! Each removal reason in the paired-end text report is counted under its own
//! label: the length cutoff, N content and `--max_length` have a line each, and
//! `--discard_untrimmed` removals appear under none of them.

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

fn run(cores: &str) -> (String, String) {
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
    (
        std::fs::read_to_string(dir.path().join("s_R2.fq_trimming_report.txt")).unwrap(),
        std::fs::read_to_string(dir.path().join("s_R2.fq_trimming_report.json")).unwrap(),
    )
}

#[test]
fn pe_report_shorter_line_counts_only_length_removals() {
    for cores in ["1", "2"] {
        let (report, json) = run(cores);
        assert!(
            report.contains(
                "Number of sequence pairs removed because at least one read was shorter than the length cutoff (25.00%):\n                                      1\n"
            ),
            "--cores {cores}:\n{report}"
        );
        // The over-length pair is reported on its own line, not under the
        // length cutoff, and the line is emitted whatever the count.
        assert!(
            report.contains(
                "Number of pairs removed for exceeding the maximum length (25.00%):\n                                      1\n"
            ),
            "--cores {cores}:\n{report}"
        );
        assert!(
            report.contains("Number of pairs removed for N content (0.00%):"),
            "--cores {cores}:\n{report}"
        );

        // `pairs_removed` stays the total that MultiQC reads, with the
        // per-reason figures as its sub-reasons: of the four pairs, one is
        // short, one untrimmed and one over --max_length.
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let pv = &v["pair_validation"];
        assert_eq!(pv["pairs_analyzed"], 4, "--cores {cores}: {pv}");
        assert_eq!(pv["pairs_removed"], 3, "--cores {cores}: {pv}");
        assert_eq!(pv["pairs_removed_too_short"], 1, "--cores {cores}: {pv}");
        assert_eq!(pv["pairs_removed_too_long"], 1, "--cores {cores}: {pv}");
        assert_eq!(pv["pairs_removed_n"], 0, "--cores {cores}: {pv}");
    }
}
