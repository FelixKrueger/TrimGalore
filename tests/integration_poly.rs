//! #504 — a poly-A tail behind a 3' poly-G run must be trimmed.
//!
//! Each case asserts exact output sequences and a record count. "No read ends in
//! poly-A" alone is an under-trim detector that an empty file also satisfies.

use std::path::PathBuf;
use std::process::Command;
use trim_galore::fastq::FastqReader;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

fn fresh_tmpdir(slug: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(slug);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// 60 bp of insert with no poly-A or poly-G content of its own.
const INSERT: &str = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCT";

/// Pin the adapter so auto-detection does not read the fixture and pick its own.
const ADAPTER: [&str; 2] = ["-a", "AGATCGGAAGAGC"];

fn write_fastq(path: &PathBuf, reads: &[(&str, String)]) {
    let mut body = String::new();
    for (name, seq) in reads {
        body.push_str(&format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(path, body).unwrap();
}

/// Output sequences, in input order.
fn sequences(path: &PathBuf) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .step_by(4)
        .map(str::to_string)
        .collect()
}

fn run(args: &[&str]) {
    let status = Command::new(binary()).args(args).status().unwrap();
    assert!(status.success(), "trim_galore failed: {args:?}");
}

#[test]
fn poly_a_tail_behind_poly_g_is_trimmed() {
    let dir = fresh_tmpdir("tg_poly_behind_g");
    let input = dir.join("reads.fastq");
    write_fastq(
        &input,
        &[
            ("flush_A40", format!("{INSERT}{}", "A".repeat(40))),
            (
                "A40_G15",
                format!("{INSERT}{}{}", "A".repeat(40), "G".repeat(15)),
            ),
            (
                "A20_G20",
                format!("{INSERT}{}{}", "A".repeat(20), "G".repeat(20)),
            ),
            (
                "A60_G20",
                format!("{INSERT}{}{}", "A".repeat(60), "G".repeat(20)),
            ),
            ("no_tail", INSERT.to_string()),
        ],
    );

    run(&[
        "--quality",
        "20",
        "--phred33",
        "--poly_a",
        "--poly_g",
        "--dont_gzip",
        "-o",
        dir.to_str().unwrap(),
        ADAPTER[0],
        ADAPTER[1],
        input.to_str().unwrap(),
    ]);

    let out = sequences(&dir.join("reads_trimmed.fq"));
    assert_eq!(out.len(), 5, "every read must be written");
    for (i, seq) in out.iter().enumerate() {
        assert_eq!(seq, INSERT, "read {i}: expected the insert alone");
    }
}

#[test]
fn two_colour_quality_trimming_also_exposes_the_tail() {
    // `--2colour` scores G as no-signal, so the G-run goes during quality
    // trimming, ahead of both poly stages.
    let dir = fresh_tmpdir("tg_poly_2colour");
    let input = dir.join("reads.fastq");
    write_fastq(
        &input,
        &[
            (
                "A40_G15",
                format!("{INSERT}{}{}", "A".repeat(40), "G".repeat(15)),
            ),
            (
                "A60_G20",
                format!("{INSERT}{}{}", "A".repeat(60), "G".repeat(20)),
            ),
        ],
    );

    run(&[
        "--2colour",
        "20",
        "--phred33",
        "--poly_a",
        "--dont_gzip",
        "-o",
        dir.to_str().unwrap(),
        ADAPTER[0],
        ADAPTER[1],
        input.to_str().unwrap(),
    ]);

    let out = sequences(&dir.join("reads_trimmed.fq"));
    assert_eq!(out.len(), 2, "every read must be written");
    for seq in &out {
        assert_eq!(seq, INSERT);
    }
}

#[test]
fn poly_g_without_poly_a_still_trims_the_g_run() {
    // The configuration most 2-colour users reach through auto-detection: the
    // G-run goes, the tail stays.
    let dir = fresh_tmpdir("tg_poly_g_only");
    let input = dir.join("reads.fastq");
    let tail = "A".repeat(40);
    write_fastq(
        &input,
        &[("A40_G15", format!("{INSERT}{tail}{}", "G".repeat(15)))],
    );

    run(&[
        "--quality",
        "20",
        "--phred33",
        "--poly_g",
        "--dont_gzip",
        "-o",
        dir.to_str().unwrap(),
        ADAPTER[0],
        ADAPTER[1],
        input.to_str().unwrap(),
    ]);

    let out = sequences(&dir.join("reads_trimmed.fq"));
    assert_eq!(out.len(), 1);
    assert_eq!(out[0], format!("{INSERT}{tail}"));
}

#[test]
fn reads_without_poly_a_are_untouched() {
    // The other direction: without `--poly_a` the same input must come through
    // with its tails intact, so a regression that over-trims is visible here
    // rather than hiding behind a tail-free assertion.
    let dir = fresh_tmpdir("tg_poly_control");
    let input = dir.join("reads.fastq");
    let reads = [
        (
            "A40_G15",
            format!("{INSERT}{}{}", "A".repeat(40), "G".repeat(15)),
        ),
        ("no_tail", INSERT.to_string()),
    ];
    write_fastq(&input, &reads);

    run(&[
        "--quality",
        "20",
        "--phred33",
        "--no_poly_g",
        "--dont_gzip",
        "-o",
        dir.to_str().unwrap(),
        ADAPTER[0],
        ADAPTER[1],
        input.to_str().unwrap(),
    ]);

    let out = sequences(&dir.join("reads_trimmed.fq"));
    assert_eq!(
        out[0], reads[0].1,
        "poly-A tail must survive without --poly_a"
    );
    assert_eq!(out[1], INSERT);
}

/// Sequences in file order. `polyAT_R2.fastq.gz` repeats read IDs (250 copies of
/// `@seq7` among others), so records are compared positionally, not by ID.
fn records(path: &PathBuf) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| c[1].to_string())
        .collect()
}

#[test]
fn paired_r2_5prime_clip_is_independent_of_3prime_trimming() {
    // R2 takes poly-T off the 5' end and poly-G off the 3'. The 5' gate
    // `errors * 5 <= i + 1` is indexed from the 5' end, so a 3' truncation cannot
    // move it: each R2 record with poly-G is a prefix of the same record without.
    let dir = fresh_tmpdir("tg_poly_paired_r2");
    let with_g = dir.join("with_g");
    let without_g = dir.join("without_g");
    std::fs::create_dir_all(&with_g).unwrap();
    std::fs::create_dir_all(&without_g).unwrap();

    for (out, g_flag) in [(&with_g, "--poly_g"), (&without_g, "--no_poly_g")] {
        run(&[
            "--paired",
            "--poly_a",
            g_flag,
            "--dont_gzip",
            "--length",
            "0",
            "-o",
            out.to_str().unwrap(),
            ADAPTER[0],
            ADAPTER[1],
            "test_files/polyAT_R1.fastq.gz",
            "test_files/polyAT_R2.fastq.gz",
        ]);
    }

    let trimmed = records(&with_g.join("polyAT_R2_val_2.fq"));
    let untrimmed = records(&without_g.join("polyAT_R2_val_2.fq"));
    assert_eq!(trimmed.len(), 13_000, "every R2 record must be compared");
    assert_eq!(trimmed.len(), untrimmed.len(), "same records in both runs");

    // Adapter and quality trimming only touch the 3' end, so an output that is
    // not a prefix of its input lost bases at the 5' end — a poly-T clip.
    let mut input = Vec::new();
    let mut reader = FastqReader::open("test_files/polyAT_R2.fastq.gz").unwrap();
    while let Some(rec) = reader.next_record().unwrap() {
        input.push(rec.seq);
    }
    assert_eq!(input.len(), trimmed.len());

    // The invariant only has content where both ends moved. Count those.
    let mut both_ends_moved = 0;
    for (i, seq) in untrimmed.iter().enumerate() {
        let with = &trimmed[i];
        assert!(
            seq.starts_with(with.as_str()),
            "record {i}: R2 5' clip moved — {with} is not a prefix of {seq}"
        );
        if with.len() < seq.len() && !input[i].starts_with(seq.as_str()) {
            both_ends_moved += 1;
        }
    }
    assert!(
        both_ends_moved > 0,
        "no record had both a 5' clip and a 3' poly-G trim, so this proves nothing"
    );
}
