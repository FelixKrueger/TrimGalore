//! Input-format detection for FASTQ vs uBAM.
//!
//! Two-stage check (PLAN §3.1 / §5 step 2):
//!
//! 1. **Cheap byte peek.** First byte `@` → plain FASTQ. First three bytes
//!    `1F 8B 08` → gzip family; fall through.
//! 2. **Decompress first block + payload check.** For the gzip family,
//!    decompress and read 4 bytes. If those equal `BAM\1` → unaligned BAM,
//!    otherwise gzipped FASTQ.
//!
//! The decompress step is **load-bearing**. Heuristic-only detection (matching
//! on the BGZF `BC` extra-field subfield) would misclassify `bgzip x.fq` —
//! BGZF-framed FASTQ — as BAM, because the framing is identical. The only
//! safe discriminator is the decompressed payload. Both plan reviewers
//! caught this (A-C2 + B-Crit-2).

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Classification of an input file's container format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    /// Plain (uncompressed) FASTQ.
    FastqPlain,
    /// Gzip-compressed FASTQ — includes both `gzip x.fq` (plain multi-member
    /// gzip) and `bgzip x.fq` (BGZF-framed gzip). Both decompress identically
    /// via `MultiGzDecoder` downstream.
    FastqGz,
    /// Unaligned BAM — BGZF-framed, first decompressed bytes are `BAM\1`.
    UnalignedBam,
}

/// Human-readable label for an input format, for use in user-facing messages.
///
/// `InputFormat` deliberately has no `Display` impl: `{:?}` would put internal
/// variant names (`FastqGz`, `UnalignedBam`) into error text. Promoted here
/// from `clump_only.rs` so every module renders formats the same way.
pub fn input_format_label(fmt: InputFormat) -> &'static str {
    match fmt {
        InputFormat::FastqPlain => "FASTQ (plain)",
        InputFormat::FastqGz => "FASTQ (gzip)",
        InputFormat::UnalignedBam => "uBAM",
    }
}

/// Which paired-input shape is being validated. Determines both the predicate
/// and the remediation, which differ by mode: `ClumpOnlyFastqOut` rejects any
/// BAM in the pair (its problem is BAM-on-a-FASTQ-path, not heterogeneity),
/// while the others reject only a format mismatch or a two-BAM pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairedShape {
    /// `--paired`, FASTQ output.
    Trim,
    /// `--paired --output-format ubam`.
    TrimUbamOut,
    /// `--clump_only --paired`, FASTQ output. Keeps the aux-tag diagnosis that
    /// `clump_only.rs` uses for this shape: the fix here is a *flag*, not a
    /// re-shaped input, and suggesting `--paired interleaved.bam` would print a
    /// command that `main.rs`'s clump-only N=1 guard rejects.
    ClumpOnlyFastqOut,
    /// `--clump_only --paired --output-format ubam`.
    ClumpOnlyUbamOut,
}

impl PairedShape {
    /// The mode flags to echo back in a remediation command.
    fn mode(self) -> &'static str {
        match self {
            PairedShape::Trim | PairedShape::TrimUbamOut => "--paired",
            PairedShape::ClumpOnlyFastqOut | PairedShape::ClumpOnlyUbamOut => {
                "--clump_only --paired"
            }
        }
    }

    /// The output-format flag to echo back, so a pasted remediation reproduces
    /// the user's own mode.
    fn fmt_flag(self) -> &'static str {
        match self {
            PairedShape::TrimUbamOut | PairedShape::ClumpOnlyUbamOut => " --output-format ubam",
            PairedShape::Trim | PairedShape::ClumpOnlyFastqOut => "",
        }
    }
}

/// Reject paired inputs whose two members disagree on being BAM.
///
/// The predicate is a **BAM count per pair**, not format equality:
/// `InputFormat` has three variants, and a `FastqPlain` + `FastqGz` pair is
/// legal. Only BAM-vs-FASTQ within one pair is an error.
///
/// Operates on the formats detected once in `main()`, so it performs no I/O and
/// runs before any output directory is created — which is the point. The sites
/// this replaces fired only after adapter auto-detection had scanned the
/// inputs, the output directory existed, and (on multi-pair input) earlier
/// pairs had already been written to disk.
///
/// Callers must skip the specialty modes (`--hardtrim5/3`, `--clock`,
/// `--implicon`), which accept mixed pairs today, and must not call this for
/// `N == 1` (single-file interleaved uBAM, validated separately in `main()`).
///
/// Assumes `Cli::validate()` has already enforced an even input count and
/// R1 != R2 per pair. `inputs.len() == formats.len()` is `debug_assert`ed,
/// since the `chunks(2)` indexing depends on it.
pub fn reject_bam_format_mismatch_in_pair(
    inputs: &[std::path::PathBuf],
    formats: &[InputFormat],
    shape: PairedShape,
) -> Result<()> {
    // A hard check, not a `debug_assert`: release builds set no
    // `debug-assertions`, so a debug-only assert would compile out of every
    // shipped binary — and both violations fail badly. A short `formats` makes
    // `zip` truncate, so later pairs go unexamined and this returns `Ok(())`
    // *without having checked*: a validation guard that silently passes, which
    // is the failure class this whole function exists to prevent. An odd length
    // would panic on `paths[1]` instead of erroring. Neither is reachable today
    // (`input_formats` is built by mapping over `cli.input`, and
    // `Cli::validate` rejects odd N for every paired mode), so this is
    // defence-in-depth for future callers of a `pub` function.
    anyhow::ensure!(
        inputs.len() == formats.len() && inputs.len().is_multiple_of(2),
        "internal error: reject_bam_format_mismatch_in_pair got {} inputs and {} formats \
         (must be equal and even). Please report this at \
         https://github.com/FelixKrueger/TrimGalore/issues",
        inputs.len(),
        formats.len()
    );

    let total_pairs = inputs.len() / 2;
    let is_bam = |f: &InputFormat| matches!(f, InputFormat::UnalignedBam);

    for (pair_idx, (paths, fmts)) in inputs
        .chunks_exact(2)
        .zip(formats.chunks_exact(2))
        .enumerate()
    {
        let n_bam = fmts.iter().filter(|f| is_bam(f)).count();
        if n_bam == 0 {
            continue;
        }

        // "Pair 3 of 5 is " when there are several pairs, plain "Got " when
        // there is one — "Pair 1 of 1" is noise on the common case, but the
        // clause still needs a subject to read as a sentence.
        let pair_prefix = if total_pairs > 1 {
            format!("Pair {} of {} is ", pair_idx + 1, total_pairs)
        } else {
            "Got ".to_string()
        };
        let (r1, r2) = (paths[0].display(), paths[1].display());

        // --clump_only with FASTQ output: any BAM in the pair is the error,
        // regardless of whether the other side is a BAM too. The user's
        // problem is BAM input on a FASTQ-output path, so the actionable
        // information is the flag to add — not the shape of the input.
        if matches!(shape, PairedShape::ClumpOnlyFastqOut) {
            let first_bam = if is_bam(&fmts[0]) { r1 } else { r2 };
            bail!(
                "uBAM input under --clump_only requires --output-format ubam \
                 (using the FASTQ output path with uBAM input would drop aux tags). \
                 Input: {first_bam}"
            );
        }

        let (mode, fmt_flag) = (shape.mode(), shape.fmt_flag());

        // Pinned-substring note: the three existing rejection tests assert on
        // "two BAM files is not supported", "uBAM paired mode expects",
        // "single interleaved" and "same format". Rust's `\`-continuation eats
        // the newline AND the next line's leading whitespace, so each of those
        // phrases must stay on one source line and every break must carry its
        // space BEFORE the backslash.
        if n_bam == 2 {
            bail!(
                "{mode} with two BAM files is not supported. \
                 uBAM paired mode expects a single interleaved file: \
                 `trim_galore {mode}{fmt_flag} interleaved.bam`. \
                 {pair_prefix}two BAM files: {r1} and {r2}. \
                 Combine them into one mate-adjacent file first: \
                 `samtools merge -n -o interleaved.bam {r1} {r2}`."
            );
        }

        // Exactly one BAM. Deliberately NOT offered the single-interleaved-file
        // remediation: this is overwhelmingly a mis-typed filename, and issue
        // #363 is precisely about that advice pointing the wrong way.
        bail!(
            "{mode} requires both inputs of a pair to be the same format. \
             {pair_prefix}mixed: {r1} is {} and {r2} is {}. \
             Pass two FASTQ files, or a single interleaved uBAM. \
             If you meant two FASTQ files, check for a mis-typed filename.",
            input_format_label(fmts[0]),
            input_format_label(fmts[1])
        );
    }

    Ok(())
}

/// Peek the first bytes of `path` and classify by content (NOT by filename).
///
/// See module-level docs for the algorithm and the rationale for the
/// decompress-and-check step.
pub fn detect_input_format(path: &Path) -> Result<InputFormat> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open input file: {}", path.display()))?;
    let mut peek = [0u8; 4];
    let n = file
        .read(&mut peek)
        .with_context(|| format!("Failed to read from {}", path.display()))?;

    if n == 0 {
        bail!("Input file '{}' is empty", path.display());
    }

    if peek[0] == b'@' {
        return Ok(InputFormat::FastqPlain);
    }

    // Gzip family — magic bytes `1F 8B 08`. The third byte is the
    // compression method; flate is what every modern gzip implementation
    // uses, and the BGZF spec also fixes it at 08.
    if n >= 3 && peek[0] == 0x1F && peek[1] == 0x8B && peek[2] == 0x08 {
        // Re-open from byte 0 (cheaper than seeking and sidesteps any
        // interaction between `MultiGzDecoder` and the already-consumed
        // prefix; the alternative `PeekReader` wrapper in PLAN §5 step 2.3
        // is a documented future optimization, not required for correctness).
        let file = File::open(path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut payload = [0u8; 4];
        let nread = decoder.read(&mut payload).with_context(|| {
            format!(
                "Failed to decompress first block of '{}' for format detection",
                path.display()
            )
        })?;
        if nread == 4 && &payload == b"BAM\x01" {
            return Ok(InputFormat::UnalignedBam);
        }
        return Ok(InputFormat::FastqGz);
    }

    bail!(
        "Input '{}' is not recognised as FASTQ (plain or gzipped) or unaligned BAM",
        path.display()
    );
}

/// Open a sync (single-threaded) reader for `path`, dispatching by detected
/// format. Returns the reader as a boxed trait object so callers can stay
/// reader-source-agnostic. `preserve_tags` is honoured for BAM input;
/// silently ignored for FASTQ.
pub fn open_sync_reader(
    path: &Path,
    preserve_tags: &[String],
) -> Result<Box<dyn crate::fastq::RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => {
            Ok(Box::new(crate::fastq::FastqReader::open(path)?))
        }
        InputFormat::UnalignedBam => Ok(Box::new(
            crate::bam::BamReader::open(path)?.with_preserved_tags(preserve_tags),
        )),
    }
}

/// Open a threaded (background-decompression) reader for `path`, dispatching
/// by detected format. `preserve_tags` is honoured for BAM input; silently
/// ignored for FASTQ.
pub fn open_threaded_reader(
    path: &Path,
    preserve_tags: &[String],
) -> Result<Box<dyn crate::fastq::RecordSource>> {
    match detect_input_format(path)? {
        InputFormat::FastqPlain | InputFormat::FastqGz => {
            Ok(Box::new(crate::fastq::FastqReader::open_threaded(path)?))
        }
        InputFormat::UnalignedBam => Ok(Box::new(crate::bam::BamReader::open_threaded_with_tags(
            path,
            preserve_tags,
        )?)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn fresh_tmpdir(slug: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(slug);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn detect_fastq_plain_from_at_sign() -> Result<()> {
        let dir = fresh_tmpdir("tg_format_plain");
        let p = dir.join("a.fq");
        std::fs::write(&p, b"@read1\nACGT\n+\nIIII\n")?;
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqPlain);
        Ok(())
    }

    #[test]
    fn detect_fastq_gz_from_plain_gzip() -> Result<()> {
        let dir = fresh_tmpdir("tg_format_gz");
        let p = dir.join("a.fq.gz");
        let mut e = GzEncoder::new(std::fs::File::create(&p)?, Compression::default());
        e.write_all(b"@read1\nACGT\n+\nIIII\n")?;
        e.finish()?;
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqGz);
        Ok(())
    }

    #[test]
    fn detect_unaligned_bam_via_decompressed_magic() -> Result<()> {
        let p = std::path::PathBuf::from("test_files/ubam_test.bam");
        assert_eq!(detect_input_format(&p)?, InputFormat::UnalignedBam);
        Ok(())
    }

    #[test]
    fn detect_paired_ubam_via_decompressed_magic() -> Result<()> {
        let p = std::path::PathBuf::from("test_files/ubam_paired_test.bam");
        assert_eq!(detect_input_format(&p)?, InputFormat::UnalignedBam);
        Ok(())
    }

    /// LOAD-BEARING (PLAN_REVIEW_A C2 + PLAN_REVIEW_B Crit-2). `bgzip x.fq`
    /// produces BGZF-framed FASTQ — the same framing as BAM. The only safe
    /// discriminator is the decompressed payload (`@` for FASTQ vs `BAM\1`
    /// for BAM). This test uses noodles' BGZF writer to produce a genuine
    /// BGZF-framed file with FASTQ contents and confirms it classifies
    /// as `FastqGz`, NOT `UnalignedBam`.
    #[test]
    fn detect_bgzipped_fastq_is_fastq_not_bam() -> Result<()> {
        use noodles::bgzf;
        let dir = fresh_tmpdir("tg_format_bgzf_fq");
        let p = dir.join("a.fq.bgz");
        {
            let mut w = bgzf::Writer::new(std::fs::File::create(&p)?);
            w.write_all(b"@read1\nACGT\n+\nIIII\n")?;
            w.finish()?;
        }
        assert_eq!(detect_input_format(&p)?, InputFormat::FastqGz);
        Ok(())
    }

    #[test]
    fn detect_rejects_random_binary() {
        let dir = fresh_tmpdir("tg_format_neg");
        let p = dir.join("noise.bin");
        std::fs::write(&p, [0u8, 1, 2, 3, 4, 5, 6, 7]).unwrap();
        let r = detect_input_format(&p);
        assert!(r.is_err(), "non-FASTQ non-BAM input must error");
    }

    #[test]
    fn detect_empty_file_errors() {
        let dir = fresh_tmpdir("tg_format_empty");
        let p = dir.join("empty.fq");
        std::fs::write(&p, b"").unwrap();
        let r = detect_input_format(&p);
        assert!(r.is_err());
    }

    // ── reject_bam_format_mismatch_in_pair (#363) ──────────────────────
    //
    // Pure decision logic over pre-detected formats, so these need no
    // fixtures. The end-to-end messages are covered by
    // tests/integration_paired_format_guard.rs.

    fn paths(n: usize) -> Vec<std::path::PathBuf> {
        (0..n)
            .map(|i| std::path::PathBuf::from(format!("in{i}.fq")))
            .collect()
    }

    fn err(fmts: &[InputFormat], shape: PairedShape) -> String {
        let p = paths(fmts.len());
        reject_bam_format_mismatch_in_pair(&p, fmts, shape)
            .expect_err("expected rejection")
            .to_string()
    }

    #[test]
    fn pair_guard_accepts_all_fastq() {
        for shape in [
            PairedShape::Trim,
            PairedShape::TrimUbamOut,
            PairedShape::ClumpOnlyFastqOut,
            PairedShape::ClumpOnlyUbamOut,
        ] {
            let fmts = [InputFormat::FastqGz, InputFormat::FastqGz];
            assert!(
                reject_bam_format_mismatch_in_pair(&paths(2), &fmts, shape).is_ok(),
                "all-FASTQ pair must be accepted for {shape:?}"
            );
        }
    }

    /// The predicate is BAM *count*, not format equality. `InputFormat` has
    /// three variants, and a plain+gzip pair is legal — naming the helper for
    /// "format mismatch" invites `formats[0] != formats[1]`, which would
    /// reject this working invocation.
    #[test]
    fn pair_guard_accepts_plain_plus_gzip_fastq() {
        let fmts = [InputFormat::FastqPlain, InputFormat::FastqGz];
        assert!(
            reject_bam_format_mismatch_in_pair(&paths(2), &fmts, PairedShape::Trim).is_ok(),
            "plain + gzipped FASTQ is a legal pair and must not be rejected"
        );
    }

    #[test]
    fn pair_guard_rejects_mixed_in_both_orders() {
        for fmts in [
            [InputFormat::FastqPlain, InputFormat::UnalignedBam],
            [InputFormat::UnalignedBam, InputFormat::FastqPlain],
        ] {
            let msg = err(&fmts, PairedShape::Trim);
            assert!(msg.contains("same format"), "got: {msg}");
            assert!(msg.contains("mixed"), "got: {msg}");
            // The substance of #363: the mixed case must NOT be offered the
            // single-interleaved-file remediation.
            assert!(
                !msg.contains("uBAM paired mode expects"),
                "mixed-pair message must not carry the two-BAM remediation; got: {msg}"
            );
            // Labels, not enum variant names.
            assert!(msg.contains("uBAM"), "got: {msg}");
            assert!(!msg.contains("UnalignedBam"), "got: {msg}");
            assert!(!msg.contains("FastqPlain"), "got: {msg}");
        }
    }

    #[test]
    fn pair_guard_rejects_two_bam_with_interleave_remediation() {
        let fmts = [InputFormat::UnalignedBam, InputFormat::UnalignedBam];
        let msg = err(&fmts, PairedShape::Trim);
        assert!(msg.contains("two BAM files is not supported"), "got: {msg}");
        assert!(msg.contains("uBAM paired mode expects"), "got: {msg}");
        assert!(msg.contains("single interleaved"), "got: {msg}");
        // Verified remediation: `samtools collate` cannot combine two files —
        // its second positional is the temp prefix, so it silently emits only
        // the first file's records.
        assert!(msg.contains("samtools merge -n"), "got: {msg}");
        // Both files named, not just one.
        assert!(
            msg.contains("in0.fq") && msg.contains("in1.fq"),
            "got: {msg}"
        );
    }

    /// `--clump_only` with FASTQ output rejects on ANY BAM in the pair, not on
    /// mismatch: the problem is BAM-on-a-FASTQ-path. Suggesting a single
    /// interleaved uBAM there would print a command the binary rejects.
    #[test]
    fn pair_guard_clump_only_fastq_out_keeps_aux_tag_diagnosis() {
        for fmts in [
            [InputFormat::FastqGz, InputFormat::UnalignedBam],
            [InputFormat::UnalignedBam, InputFormat::FastqGz],
            [InputFormat::UnalignedBam, InputFormat::UnalignedBam],
        ] {
            let msg = err(&fmts, PairedShape::ClumpOnlyFastqOut);
            assert!(msg.contains("--output-format ubam"), "got: {msg}");
            assert!(msg.contains("drop aux tags"), "got: {msg}");
            assert!(
                !msg.contains("interleaved.bam"),
                "clump-only FASTQ output must not suggest an interleaved uBAM \
                 (main() rejects that invocation); got: {msg}"
            );
        }
    }

    #[test]
    fn pair_guard_clump_only_fastq_out_names_the_bam_side() {
        let fmts = [InputFormat::FastqGz, InputFormat::UnalignedBam];
        assert!(err(&fmts, PairedShape::ClumpOnlyFastqOut).contains("in1.fq"));
        let fmts = [InputFormat::UnalignedBam, InputFormat::FastqGz];
        assert!(err(&fmts, PairedShape::ClumpOnlyFastqOut).contains("in0.fq"));
    }

    #[test]
    fn pair_guard_reports_the_offending_pair_index() {
        // N=6: pairs 1 and 3 are clean, pair 2 is mixed.
        let fmts = [
            InputFormat::FastqGz,
            InputFormat::FastqGz,
            InputFormat::FastqGz,
            InputFormat::UnalignedBam,
            InputFormat::FastqGz,
            InputFormat::FastqGz,
        ];
        let msg = err(&fmts, PairedShape::Trim);
        assert!(msg.contains("Pair 2 of 3"), "got: {msg}");
        assert!(
            msg.contains("in2.fq") && msg.contains("in3.fq"),
            "got: {msg}"
        );
    }

    #[test]
    fn pair_guard_omits_index_for_a_single_pair() {
        let fmts = [InputFormat::FastqGz, InputFormat::UnalignedBam];
        let msg = err(&fmts, PairedShape::Trim);
        assert!(msg.contains("Got mixed"), "got: {msg}");
        assert!(!msg.contains("Pair 1 of 1"), "got: {msg}");
    }

    #[test]
    fn pair_guard_remediation_echoes_the_users_mode() {
        let two_bam = [InputFormat::UnalignedBam, InputFormat::UnalignedBam];
        assert!(
            err(&two_bam, PairedShape::TrimUbamOut)
                .contains("trim_galore --paired --output-format ubam interleaved.bam")
        );
        assert!(
            err(&two_bam, PairedShape::ClumpOnlyUbamOut)
                .contains("trim_galore --clump_only --paired --output-format ubam interleaved.bam")
        );
        assert!(err(&two_bam, PairedShape::Trim).contains("trim_galore --paired interleaved.bam"));
    }
}
