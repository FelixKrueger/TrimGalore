//! Output file naming conventions and I/O utilities.
//!
//! Implements TrimGalore-compatible output naming:
//! - Single-end: *_trimmed.fq(.gz)
//! - Paired-end: *_val_1.fq(.gz) / *_val_2.fq(.gz)
//! - Unpaired: *_unpaired_1.fq(.gz) / *_unpaired_2.fq(.gz)
//! - Reports: *_trimming_report.txt

use anyhow::{Context, Result};
use std::path::{Component, Path, PathBuf};

/// True iff `path` ends with a `.gz` extension.
///
/// **Filename-based on purpose, and no longer the same test the reader uses.**
/// `FastqReader` decides whether to decompress by sniffing the file's first
/// three bytes, because the name is unreliable in both directions. This
/// function keeps looking at the name, because it answers a different
/// question: not "is this input compressed" but "should the output be".
///
/// Used to mirror input compression in the output filename / writer:
/// `plain.fastq` → `plain_trimmed.fq` (plain), `plain.fastq.gz` →
/// `plain_trimmed.fq.gz`. Matches Perl v0.6.x behaviour.
///
/// The two can therefore disagree, and today they do: a `.bgz` input is
/// decompressed correctly but produces plain output, because this returns
/// false for it. Moving the output decision to the detected format would
/// change behaviour for every run, including a plain file misnamed `.gz`, so
/// it is deliberately left for a separate change. See the CHANGELOG.
///
/// Matching is ASCII-case-insensitive (#384) and must stay in step with
/// `strip_fastq_extensions`: folding one but not the other names the output
/// like the gzipped convention while writing it plain, or vice versa.
pub fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

/// Case-folded (ASCII lowercase) view of a path; the folding half of `collision_key`.
fn norm_path(p: &Path) -> String {
    p.to_string_lossy().to_ascii_lowercase()
}

/// Lexically normalised path: absolutised, with `..` folded against the component
/// stack, so `./x`, `<cwd>/x` and `a/../x` name one path. No filesystem access, so a
/// symlinked path is not resolved.
fn lexical_normalise(p: &Path) -> PathBuf {
    let abs = std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf());
    let mut stack: Vec<Component> = Vec::new();
    for c in abs.components() {
        match c {
            // `/..` is `/` on POSIX; `..` above a relative root has to be kept.
            Component::ParentDir => match stack.last() {
                Some(Component::Normal(_)) => {
                    stack.pop();
                }
                Some(Component::RootDir) | Some(Component::Prefix(_)) => {}
                _ => stack.push(c),
            },
            Component::CurDir => {}
            other => stack.push(other),
        }
    }
    stack.iter().collect()
}

/// Is this the same file? Case-**sensitive**, because on a case-sensitive filesystem
/// `X` and `x` are two files and calling them one would reject valid input. Used for
/// input identity in `Cli::validate` (issues #383, #389).
pub fn path_identity_key(p: &Path) -> String {
    lexical_normalise(p).to_string_lossy().into_owned()
}

/// Would these two paths be the same *output* file? Case-**folded**, because outputs
/// differing only in case alias each other on APFS/NTFS (issues #216, #383). Input
/// identity asks a different question and uses the case-preserving `path_identity_key`.
///
/// Pragmatic trade-off: on a case-sensitive filesystem (Linux, or opt-in
/// case-sensitive APFS) two outputs differing only in case are distinct, and this
/// rejects them — a loud early error in preference to silent data loss.
pub fn collision_key(p: &Path) -> String {
    norm_path(&lexical_normalise(p))
}

/// Set this to any value to make [`preflight_output_collisions`] print the candidate
/// list it is about to check. Test hook for `tests/integration_preflight_tripwire.rs`,
/// which asserts that every file a run creates was in that list (#414).
pub const DUMP_PLANNED_ENV: &str = "TRIM_GALORE_DUMP_PLANNED_OUTPUTS";

/// Line prefix for one dumped candidate.
pub const DUMP_PLANNED_PREFIX: &str = "TRIM_GALORE_PLANNED\t";

/// Render the dump lines for one candidate list, tagged with the call site.
///
/// `site` identifies which of the dispatch arms built this list, so a tripwire failure
/// can name the candidate list to add to. Pure, so it is testable without touching the
/// environment.
///
/// A path containing a newline would split across lines and defeat the parse; none can
/// occur in the fixtures, and `display()` is lossy for non-UTF-8 paths either way.
pub fn dump_planned_lines(planned: &[PlannedOutput], site: &str) -> Vec<String> {
    planned
        .iter()
        .map(|(path, _)| format!("{}{}\t{}", DUMP_PLANNED_PREFIX, site, path.display()))
        .collect()
}

/// The directory every output of a paired run lands in: `--output_dir` when given,
/// else R1's parent.
///
/// One derivation for primaries, reports and the `--passthrough` carrier alike (#398).
/// `parent()` yields `Some("")` for a bare filename, so joined paths stay bare and the
/// `"."` arm is reachable only for a root path — returning `"."` here would prefix every
/// path, and every message naming one, with `./`.
pub fn pair_output_dir(input_r1: &Path, output_dir: Option<&Path>) -> PathBuf {
    output_dir
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf())
}

/// Remediation offered when nothing mode-specific applies.
const GENERIC_ADVICE: &str = "Check that inputs produce distinct output paths \
                              (e.g., different source directories or `--output_dir`).";

/// Which input(s) a planned output is named from. The pre-flight cannot tell a user to
/// "rename one input" without it (#397).
#[derive(Clone, Debug)]
pub enum OutputSource {
    /// Named from one input: reports, SE primaries, `--demux` per-barcode files.
    Input(PathBuf),
    /// Named from a pair. Paired primaries take the filename from one mate and the
    /// directory from R1, so neither mate alone explains the path.
    Pair(PathBuf, PathBuf),
}

impl OutputSource {
    /// Case-**sensitive** identity, for "did these two candidates come from the same
    /// input?". Folding here would call `X` and `x` one source on a case-sensitive
    /// filesystem — the D8/D13 mistake, in a new place.
    fn identity(&self) -> String {
        match self {
            OutputSource::Input(p) => path_identity_key(p),
            // `\0` cannot occur in a path, so it cannot make two distinct pairs collide.
            OutputSource::Pair(a, b) => {
                format!("{}\0{}", path_identity_key(a), path_identity_key(b))
            }
        }
    }
}

impl std::fmt::Display for OutputSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OutputSource::Input(p) => write!(f, "{}", p.display()),
            OutputSource::Pair(a, b) => write!(f, "{} + {}", a.display(), b.display()),
        }
    }
}

/// A prospective output path and the input(s) it is named from.
pub type PlannedOutput = (PathBuf, OutputSource);

/// Reject duplicate output paths, or an output that would overwrite an input, before
/// any file is opened. Called by every `main.rs` dispatch path that writes more than
/// one file: SE trim, `--paired` trim, `--hardtrim5/3`, `--clock`/`--implicon` and
/// `--clump_only` — all four output formats. A new dispatch path needs a call here too.
///
/// Each candidate carries the input it is named from, so a refusal can name the files the
/// user has to change (#397). Note what the tuple does and does not buy: it makes an
/// *unattributed* candidate impossible, but it cannot detect a *missing* one — which is
/// what #383, #385, #388, #391 and #409 all were.
/// `#[track_caller]` so the dump can name which dispatch arm built the list.
#[track_caller]
pub fn preflight_output_collisions(
    planned: &[PlannedOutput],
    inputs: &[PathBuf],
    hint: Option<&str>,
) -> Result<()> {
    if std::env::var_os(DUMP_PLANNED_ENV).is_some() {
        let at = std::panic::Location::caller();
        let site = format!("{}:{}", at.file(), at.line());
        for line in dump_planned_lines(planned, &site) {
            eprintln!("{line}");
        }
    }

    let input_keys: std::collections::HashMap<String, &PathBuf> =
        inputs.iter().map(|p| (collision_key(p), p)).collect();
    let mut seen: std::collections::HashMap<String, PlannedOutput> =
        std::collections::HashMap::with_capacity(planned.len());

    for (path, source) in planned {
        let key = collision_key(path);
        if let Some(input) = input_keys.get(&key) {
            anyhow::bail!(
                "Output path collision (case-insensitive, for APFS/NTFS safety): the output \
                 named from {} would be written to {}, which is also one of its inputs. \
                 If that file is an earlier run's output, drop it from the input list; \
                 otherwise use `--output_dir` to write elsewhere.",
                source,
                input.display()
            );
        }
        if let Some((prev_path, prev_source)) = seen.insert(key, (path.clone(), source.clone())) {
            let same_path = prev_path == *path;
            let same_source = prev_source.identity() == source.identity();
            if same_path && same_source {
                // 2c — one input feeding two candidates. The hint is deliberately ignored:
                // "rename one input" cannot be followed when there is only one to rename.
                anyhow::bail!(
                    "Output path collision (case-insensitive, for APFS/NTFS safety): two \
                     outputs named from {} would be written to the same file, {}. List each \
                     input once — a file that appears in more than one pair is reported once \
                     per pair — or pass `--no_report_file` if only the reports collide.",
                    source,
                    path.display()
                );
            } else if same_path {
                // 2a — two sources, one rendered path: print the path once.
                anyhow::bail!(
                    "Output path collision (case-insensitive, for APFS/NTFS safety): the \
                     outputs named from {} and {} would be written to the same file, {}. {}",
                    prev_source,
                    source,
                    path.display(),
                    hint.unwrap_or(GENERIC_ADVICE)
                );
            } else {
                // 2b — the spellings differ (fold-equal or `./x` vs `x`), so show both.
                anyhow::bail!(
                    "Output path collision (case-insensitive, for APFS/NTFS safety): {} \
                     (from {}) and {} (from {}) would be written to the same file. {}",
                    prev_path.display(),
                    prev_source,
                    path.display(),
                    source,
                    hint.unwrap_or(GENERIC_ADVICE)
                );
            }
        }
    }
    Ok(())
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

/// Generate the trimmed output filename for single-end mode, uBAM variant.
///
/// `<stem>_trimmed.bam` — compression is implicit (BGZF, always on for BAM).
/// Per PLAN v2.1 §3.2.
pub fn single_end_bam_output_name(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input));
    let filename = format!("{}_trimmed.bam", stem);

    match output_dir {
        Some(dir) => dir.join(&filename),
        None => input.parent().unwrap_or(Path::new(".")).join(&filename),
    }
}

/// Generate the validated output filename for paired-end mode, uBAM variant.
///
/// Returns a SINGLE interleaved BAM path (`<stem>_val.bam`) — no `_1`/`_2`
/// suffix. Matches samtools/Picard/fgbio convention and #317's
/// `BamReader::open_paired_interleaved` reader-side expectation. Per PLAN
/// v2.1 §3.2.
///
/// When `--basename foo` is set, the path is `foo_val.bam`; otherwise the
/// stem derives from R1's filename (R2's is ignored — same convention as
/// the FASTQ path).
pub fn paired_bam_output_name(
    input_r1: &Path,
    _input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input_r1));
    let filename = format!("{}_val.bam", stem);

    let dir = pair_output_dir(input_r1, output_dir);

    dir.join(&filename)
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

    // When --basename is supplied, both stems are exactly that basename:
    // `foo` → `foo_val_1.fq.gz` / `foo_val_2.fq.gz`. The `_val_{1,2}`
    // suffix carries the pair-side distinction; there's no `_R1`/`_R2`
    // segment between the basename and the suffix. Matches Perl v0.6.5+
    // behaviour (CHANGELOG: "In a `--paired --basename BASE` scenario,
    // the output files will now be called `BASE_val_1.fq.gz
    // BASE_val_2.fq.gz` as described in the documentation"). Beta.5
    // had a `_R1`/`_R2` interpolation here that broke the documented
    // contract — silently miss-pathed outputs for any nf-core /
    // Snakemake pipeline globbing the documented `${basename}_val_*`.
    // See #244.
    let (stem1, stem2) = match basename {
        Some(b) => (b.to_string(), b.to_string()),
        None => (
            strip_fastq_extensions(input_r1),
            strip_fastq_extensions(input_r2),
        ),
    };

    let f1 = format!("{}_val_1{}", stem1, ext);
    let f2 = format!("{}_val_2{}", stem2, ext);

    let dir = pair_output_dir(input_r1, output_dir);

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

    // Same `--basename` semantic as paired_end_output_names — the
    // `_unpaired_{1,2}` suffix carries the pair-side, no `_R{1,2}`
    // segment is interpolated. Matches Perl v0.6.5+ behaviour. See #244.
    let (stem1, stem2) = match basename {
        Some(b) => (b.to_string(), b.to_string()),
        None => (
            strip_fastq_extensions(input_r1),
            strip_fastq_extensions(input_r2),
        ),
    };

    let f1 = format!("{}_unpaired_1{}", stem1, ext);
    let f2 = format!("{}_unpaired_2{}", stem2, ext);

    let dir = pair_output_dir(input_r1, output_dir);

    (dir.join(&f1), dir.join(&f2))
}

/// Generate the passthrough output filename for paired-end mode with
/// `--passthrough` (10X Multiome / scATAC cell-barcode carrier).
///
/// The naming follows the same `--basename` semantic as
/// `paired_end_output_names`: when `--basename foo` is set the output is
/// `foo_passthrough.{fq,fq.gz}`; otherwise the stem derives from the
/// passthrough input filename via `strip_fastq_extensions`.
///
/// Note: `gzip` is the resolved global flag (driven by R1's input
/// compression in `main.rs`), not the passthrough input's own extension —
/// output compression is uniform across the three pair outputs, see plan
/// v2 §Assumptions §11.
///
/// The directory comes from `pair_output_dir`, not from the carrier's own parent:
/// the carrier is one of the pair's three outputs and lands with the other two (#398).
pub fn passthrough_output_name(
    input_r1: &Path,
    input_passthrough: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input_passthrough));

    let ext = if gzip {
        "_passthrough.fq.gz"
    } else {
        "_passthrough.fq"
    };
    let filename = format!("{}{}", stem, ext);

    let dir = pair_output_dir(input_r1, output_dir);
    dir.join(&filename)
}

/// Generate the reorder-only (`--clump_only`) output filename for single-end mode.
///
/// Follows the same input-extension-stripping convention as
/// `single_end_output_name`, but emits `_clumped` in place of `_trimmed` so
/// downstream globbing pipelines don't mistake the output for a trim result:
///
/// - `.fastq.gz` → `_clumped.fq.gz` (or `_clumped.fq` if `--dont_gzip`)
/// - `.fastq`    → `_clumped.fq`(.gz)
/// - `.fq.gz`    → `_clumped.fq`(.gz)
/// - `.fq`       → `_clumped.fq`(.gz)
///
/// When `--basename BASE` is supplied the output is `BASE_clumped.fq(.gz)`.
pub fn clumped_output_name(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input));

    let ext = if gzip {
        "_clumped.fq.gz"
    } else {
        "_clumped.fq"
    };
    let filename = format!("{}{}", stem, ext);

    match output_dir {
        Some(dir) => dir.join(&filename),
        None => input.parent().unwrap_or(Path::new(".")).join(&filename),
    }
}

/// Generate the reorder-only output filenames for paired-end mode.
///
/// Returns `(clumped_1_path, clumped_2_path)`. Follows the same
/// `--basename` semantic as `paired_end_output_names`: with basename "foo"
/// both mates use `foo_clumped_{1,2}` (no `_R1`/`_R2` interpolation).
pub fn clumped_paired_output_names(
    input_r1: &Path,
    input_r2: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
    gzip: bool,
) -> (PathBuf, PathBuf) {
    let ext = if gzip { ".fq.gz" } else { ".fq" };

    let (stem1, stem2) = match basename {
        Some(b) => (b.to_string(), b.to_string()),
        None => (
            strip_fastq_extensions(input_r1),
            strip_fastq_extensions(input_r2),
        ),
    };

    let f1 = format!("{}_clumped_1{}", stem1, ext);
    let f2 = format!("{}_clumped_2{}", stem2, ext);

    let dir = pair_output_dir(input_r1, output_dir);

    (dir.join(&f1), dir.join(&f2))
}

/// Generate the reorder-only (`--clump_only --output-format ubam`) output
/// filename for single-end mode.
///
/// Output: `<stem>_clumped.bam`. Mirrors `hardtrim_bam_output_name`'s
/// convention — `strip_fastq_extensions` handles both FASTQ and BAM
/// input extensions (BAM via `file_stem()` fallback). Honors `--basename`.
pub fn clumped_bam_output_name(
    input: &Path,
    output_dir: Option<&Path>,
    basename: Option<&str>,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input));
    let filename = format!("{}_clumped.bam", stem);

    match output_dir {
        Some(dir) => dir.join(&filename),
        None => input.parent().unwrap_or(Path::new(".")).join(&filename),
    }
}

/// Generate the reorder-only output filename for paired-end BAM output.
///
/// Returns a SINGLE interleaved BAM path (`<stem>_clumped.bam`) — no
/// `_1`/`_2` suffix. Matches samtools/Picard/fgbio mate-adjacent
/// convention and mirrors `paired_bam_output_name` in shape.
///
/// The stem derives from `input_r1`'s filename; `_input_r2` is accepted
/// for API-parity with `paired_bam_output_name` but not used. When
/// `--basename foo` is set, the path is `foo_clumped.bam`.
pub fn clumped_paired_bam_output_name(
    input_r1: &Path,
    _input_r2: Option<&Path>,
    output_dir: Option<&Path>,
    basename: Option<&str>,
) -> PathBuf {
    let stem = basename
        .map(|b| b.to_string())
        .unwrap_or_else(|| strip_fastq_extensions(input_r1));
    let filename = format!("{}_clumped.bam", stem);

    let dir = pair_output_dir(input_r1, output_dir);

    dir.join(&filename)
}

/// Generate the `--clump_only` reorder report filename.
///
/// Deliberately distinct from `report_name`'s `*_trimming_report.txt` so
/// downstream tools (nf-core/rnaseq's MultiQC integration) that scan
/// `*_trimming_report.*` don't mis-classify an empty-of-trim-stats file.
/// Text-only (no JSON): the reorder report is short enough to grep.
pub fn clumping_report_name(input: &Path, output_dir: Option<&Path>) -> PathBuf {
    let input_name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    let report = format!("{}_clumping_report.txt", input_name);

    match output_dir {
        Some(dir) => dir.join(&report),
        None => input.parent().unwrap_or(Path::new(".")).join(&report),
    }
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

/// Text + JSON trimming-report paths for one mate of a **two-file** pair.
///
/// Both land in the pair's output directory — R1's parent unless `--output_dir`
/// overrides — while the filename stays keyed on the mate's own input. Taking
/// `input_r1` rather than an already-computed directory is the point: a caller
/// cannot get the directory rule wrong without passing the wrong mate.
///
/// Per-mate, not all four paths at once, because the two writer sites build a
/// `PairedReportFile` per mate and the two candidate sites push txt-then-json
/// per mate — an order `tests/integration_output_collision.rs` pins.
pub fn paired_report_names(
    input: &Path,
    input_r1: &Path,
    output_dir: Option<&Path>,
) -> (PathBuf, PathBuf) {
    let dir = pair_output_dir(input_r1, output_dir);
    (
        report_name(input, Some(&dir)),
        json_report_name(input, Some(&dir)),
    )
}

// ─── Single-file interleaved paired output (#423) ───────────────────────────
//
// One derivation per arm, consumed by both the writers and the collision
// pre-flight. Keyed on `file_stem()`, not `file_name()`: these
// arms have always named their outputs `ip_val_1.fq` from `ip.bam`, where
// `report_name` would give `ip.bam_trimming_report.txt`. Composing these from
// `report_name` renames every interleaved report.
//
// Every field added here must appear in that type's `planned()`; the exhaustive
// destructure there makes a field that does not a build error.

/// The two mates' trimming-report paths for a single interleaved pair.
///
/// `_R1`/`_R2` infix rather than per-input names, because both mates come from
/// one input and would otherwise share a path.
#[derive(Clone, Debug)]
pub struct InterleavedReports {
    pub r1_txt: PathBuf,
    pub r1_json: PathBuf,
    pub r2_txt: PathBuf,
    pub r2_json: PathBuf,
}

impl InterleavedReports {
    fn new(stem: &str, dir: &Path) -> Self {
        Self {
            r1_txt: dir.join(format!("{stem}_R1_trimming_report.txt")),
            r1_json: dir.join(format!("{stem}_R1_trimming_report.json")),
            r2_txt: dir.join(format!("{stem}_R2_trimming_report.txt")),
            r2_json: dir.join(format!("{stem}_R2_trimming_report.json")),
        }
    }

    /// txt before json within a mate, the order the two-file paired sites use and
    /// `tests/integration_output_collision.rs` pins.
    fn planned(&self, src: &OutputSource) -> Vec<PlannedOutput> {
        let Self {
            r1_txt,
            r1_json,
            r2_txt,
            r2_json,
        } = self;
        vec![
            (r1_txt.clone(), src.clone()),
            (r1_json.clone(), src.clone()),
            (r2_txt.clone(), src.clone()),
            (r2_json.clone(), src.clone()),
        ]
    }
}

/// Every path a FASTQ-output single-file interleaved paired run writes.
#[derive(Clone, Debug)]
pub struct InterleavedFastqOutputs {
    pub val_1: PathBuf,
    pub val_2: PathBuf,
    /// `Some` under `--retain_unpaired`.
    pub unpaired: Option<(PathBuf, PathBuf)>,
    /// `None` under `--no_report_file`.
    pub reports: Option<InterleavedReports>,
}

impl InterleavedFastqOutputs {
    pub fn new(
        input: &Path,
        output_dir: Option<&Path>,
        gzip: bool,
        retain_unpaired: bool,
        reports: bool,
    ) -> Self {
        let (stem, dir) = interleaved_stem_and_dir(input, output_dir);
        let ext = if gzip { ".fq.gz" } else { ".fq" };
        Self {
            val_1: dir.join(format!("{stem}_val_1{ext}")),
            val_2: dir.join(format!("{stem}_val_2{ext}")),
            unpaired: retain_unpaired.then(|| {
                (
                    dir.join(format!("{stem}_unpaired_1{ext}")),
                    dir.join(format!("{stem}_unpaired_2{ext}")),
                )
            }),
            reports: reports.then(|| InterleavedReports::new(&stem, &dir)),
        }
    }

    /// The pre-flight candidate list. `src` is the lone input, so every path
    /// carries the same one.
    pub fn planned(&self, src: &OutputSource) -> Vec<PlannedOutput> {
        let Self {
            val_1,
            val_2,
            unpaired,
            reports,
        } = self;
        let mut v = vec![(val_1.clone(), src.clone()), (val_2.clone(), src.clone())];
        if let Some((up_1, up_2)) = unpaired {
            v.push((up_1.clone(), src.clone()));
            v.push((up_2.clone(), src.clone()));
        }
        if let Some(reports) = reports {
            v.extend(reports.planned(src));
        }
        v
    }
}

/// Every path a uBAM-output single-file interleaved paired run writes.
#[derive(Clone, Debug)]
pub struct InterleavedBamOutputs {
    pub val: PathBuf,
    /// `None` under `--no_report_file`.
    pub reports: Option<InterleavedReports>,
}

impl InterleavedBamOutputs {
    pub fn new(input: &Path, output_dir: Option<&Path>, reports: bool) -> Self {
        let (stem, dir) = interleaved_stem_and_dir(input, output_dir);
        Self {
            val: dir.join(format!("{stem}_val.bam")),
            reports: reports.then(|| InterleavedReports::new(&stem, &dir)),
        }
    }

    /// The pre-flight candidate list. `src` is the lone input, so every path
    /// carries the same one.
    pub fn planned(&self, src: &OutputSource) -> Vec<PlannedOutput> {
        let Self { val, reports } = self;
        let mut v = vec![(val.clone(), src.clone())];
        if let Some(reports) = reports {
            v.extend(reports.planned(src));
        }
        v
    }
}

/// Stem and directory shared by both interleaved arms. `pair_output_dir` with
/// the lone input standing in for R1.
fn interleaved_stem_and_dir(input: &Path, output_dir: Option<&Path>) -> (String, PathBuf) {
    let stem = input
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    (stem, pair_output_dir(input, output_dir))
}

/// ASCII-case-insensitive suffix strip, preserving the case of what remains (#384).
///
/// Byte-wise on purpose: slicing `&str` at `len - suffix.len()` panics when the
/// index lands inside a multi-byte character, and extensionless non-ASCII names
/// are valid input. A matched all-ASCII tail guarantees the index is a boundary.
fn strip_suffix_ignore_ascii_case<'a>(name: &'a str, suffix: &str) -> Option<&'a str> {
    let idx = name.len().checked_sub(suffix.len())?;
    name.as_bytes()[idx..]
        .eq_ignore_ascii_case(suffix.as_bytes())
        .then(|| &name[..idx])
}

/// Strip FASTQ extensions from a filename, returning just the base stem.
///
/// Handles `.fastq` / `.fq`, each optionally followed by a gzip-family
/// suffix: `.gz`, or the `.bgz` / `.bgzf` names `bgzip` output is often
/// given. The two suffix groups are stripped in sequence rather than
/// enumerated as a combined list, so adding one compression name does not
/// multiply the cases. The `.bgz` forms are here because such a file is read
/// correctly since #374, and its output name has to follow
/// ([#381](https://github.com/FelixKrueger/TrimGalore/issues/381)).
///
/// Anything else (a `.bam` input, an unrelated extension) falls through to
/// `Path::file_stem`, which drops the final component only.
pub fn strip_fastq_extensions(path: &Path) -> String {
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    // Strip one gzip-family suffix, if present. These are disjoint: `.bgz`
    // does not end in `.gz`, so the order within the list does not matter.
    let stem = [".gz", ".bgz", ".bgzf"]
        .iter()
        .find_map(|ext| strip_suffix_ignore_ascii_case(&name, ext))
        .unwrap_or(&name);

    // Then the FASTQ extension itself, longest first.
    for ext in [".fastq", ".fq"] {
        if let Some(base) = strip_suffix_ignore_ascii_case(stem, ext) {
            return base.to_string();
        }
    }

    Path::new(stem)
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string()
}

// ─── Deferred publication of record output (#428) ───────────────────────────

/// Temporary path a record output is written to before it is published.
///
/// The marker wraps the *whole* file name, dot-prefixed: `x_trimmed.fq.gz` →
/// `.x_trimmed.fq.gz.partial`. No suffix may be appended to a name that has
/// been through this — an output namer must derive from the input path, never
/// from a temporary.
///
/// The mapping is byte-exact and therefore injective, which is what
/// `tests/integration_preflight_tripwire.rs` rests on when it excludes these
/// from the collision pre-flight.
pub fn partial_output_name(final_path: &Path) -> PathBuf {
    let mut tmp = std::ffi::OsString::from(".");
    tmp.push(final_path.file_name().unwrap_or_default());
    tmp.push(".partial");
    match final_path.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir.join(tmp),
        _ => PathBuf::from(tmp),
    }
}

/// A record output whose final name appears only once its writer is closed.
///
/// Writes land on [`partial_output_name`]; [`PendingOutput::commit`] renames
/// onto the final path. Dropped without committing, the temporary is removed
/// and the final path is never created — so a refused run leaves no output,
/// and a previous run's output at that path survives untouched.
///
/// **Commit where the writer is closed, never later.** FastQC, `--demux`,
/// `write_clump_only_report` and `clump_only`'s `output_bytes` all read the
/// output by its *final* path within the same run, and the last of those is a
/// `?`-propagating `metadata` call on the success path — so a late commit
/// fails a good run with `ENOENT`.
pub struct PendingOutput {
    tmp: PathBuf,
    final_path: PathBuf,
    /// Whether `Drop` still owns `tmp`. Set *before* the rename is attempted,
    /// so a failed rename leaves the temporary on disk rather than deleting
    /// the data it was protecting.
    released: bool,
}

impl PendingOutput {
    /// Create the temporary and hand back its handle alongside the guard.
    pub fn create(final_path: &Path) -> Result<(Self, std::fs::File)> {
        let tmp = partial_output_name(final_path);
        let file = std::fs::File::create(&tmp).with_context(|| {
            format!(
                "Failed to create output {} (writing to {})",
                final_path.display(),
                tmp.display()
            )
        })?;
        Ok((
            Self {
                tmp,
                final_path: final_path.to_path_buf(),
                released: false,
            },
            file,
        ))
    }

    /// Publish the temporary under its final name. Consumes `self`.
    pub fn commit(mut self) -> Result<()> {
        self.released = true;
        std::fs::rename(&self.tmp, &self.final_path).with_context(|| {
            format!(
                "Failed to publish {} (the trimmed data is intact at {})",
                self.final_path.display(),
                self.tmp.display()
            )
        })
    }

    /// The temporary path, for error messages.
    pub fn tmp_path(&self) -> &Path {
        &self.tmp
    }
}

impl Drop for PendingOutput {
    fn drop(&mut self) {
        if !self.released {
            let _ = std::fs::remove_file(&self.tmp);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ─── #423: the paired report namers, pinned to literal paths ───────────
    //
    // Nothing else in the repo pins the interleaved `_R{1,2}_trimming_report`
    // names, and the writer and the pre-flight both derive them from here — so a
    // rename would be self-consistent and invisible everywhere else.

    #[test]
    fn paired_report_names_share_r1s_directory_and_keep_each_mates_filename() {
        let (txt, json) = paired_report_names(
            Path::new("R2/s_R2.fastq.gz"),
            Path::new("R1/s_R1.fastq.gz"),
            None,
        );
        assert_eq!(txt, PathBuf::from("R1/s_R2.fastq.gz_trimming_report.txt"));
        assert_eq!(json, PathBuf::from("R1/s_R2.fastq.gz_trimming_report.json"));
    }

    #[test]
    fn paired_report_names_honours_output_dir() {
        let (txt, json) = paired_report_names(
            Path::new("R2/s_R2.fastq.gz"),
            Path::new("R1/s_R1.fastq.gz"),
            Some(Path::new("out")),
        );
        assert_eq!(txt, PathBuf::from("out/s_R2.fastq.gz_trimming_report.txt"));
        assert_eq!(
            json,
            PathBuf::from("out/s_R2.fastq.gz_trimming_report.json")
        );
    }

    /// `file_stem`, NOT `file_name`: `report_name` would give
    /// `ip.bam_trimming_report.txt` and rename every interleaved report.
    #[test]
    fn interleaved_outputs_key_on_the_file_stem() {
        let o = InterleavedFastqOutputs::new(Path::new("d/ip.bam"), None, false, false, true);
        assert_eq!(o.val_1, PathBuf::from("d/ip_val_1.fq"));
        assert_eq!(o.val_2, PathBuf::from("d/ip_val_2.fq"));
        let r = o.reports.unwrap();
        assert_eq!(r.r1_txt, PathBuf::from("d/ip_R1_trimming_report.txt"));
        assert_eq!(r.r1_json, PathBuf::from("d/ip_R1_trimming_report.json"));
        assert_eq!(r.r2_txt, PathBuf::from("d/ip_R2_trimming_report.txt"));
        assert_eq!(r.r2_json, PathBuf::from("d/ip_R2_trimming_report.json"));
    }

    /// A `.gz`-*named* BAM: `gzip` is decided from the extension, so the
    /// primaries gain `.fq.gz` while the stem keeps the inner `.bam`.
    #[test]
    fn interleaved_outputs_two_dot_name_with_gzip() {
        let o = InterleavedFastqOutputs::new(Path::new("ip.bam.gz"), None, true, true, true);
        assert_eq!(o.val_1, PathBuf::from("ip.bam_val_1.fq.gz"));
        assert_eq!(
            o.unpaired.unwrap().1,
            PathBuf::from("ip.bam_unpaired_2.fq.gz")
        );
        assert_eq!(
            o.reports.unwrap().r1_txt,
            PathBuf::from("ip.bam_R1_trimming_report.txt")
        );
    }

    /// A bare input name yields bare output names: `Path::new("ip.bam").parent()`
    /// is `Some("")`, so the `unwrap_or(".")` in `pair_output_dir` never fires.
    #[test]
    fn interleaved_outputs_bare_filename_stay_bare() {
        let o = InterleavedFastqOutputs::new(Path::new("ip.bam"), None, false, false, false);
        assert_eq!(o.val_1, PathBuf::from("ip_val_1.fq"));
        assert!(o.unpaired.is_none(), "--retain_unpaired was not requested");
        assert!(o.reports.is_none(), "--no_report_file was in effect");
    }

    #[test]
    fn interleaved_outputs_honours_output_dir() {
        let o = InterleavedFastqOutputs::new(
            Path::new("d/ip.bam"),
            Some(Path::new("out")),
            false,
            true,
            true,
        );
        assert_eq!(o.val_1, PathBuf::from("out/ip_val_1.fq"));
        assert_eq!(o.unpaired.unwrap().0, PathBuf::from("out/ip_unpaired_1.fq"));
        assert_eq!(
            o.reports.unwrap().r2_json,
            PathBuf::from("out/ip_R2_trimming_report.json")
        );
    }

    #[test]
    fn interleaved_bam_outputs_name_one_interleaved_primary() {
        let o = InterleavedBamOutputs::new(Path::new("d/ip.bam"), None, true);
        assert_eq!(o.val, PathBuf::from("d/ip_val.bam"));
        assert_eq!(
            o.reports.unwrap().r1_txt,
            PathBuf::from("d/ip_R1_trimming_report.txt")
        );
        assert!(
            InterleavedBamOutputs::new(Path::new("d/ip.bam"), None, false)
                .reports
                .is_none()
        );
    }

    /// Renders one arm's candidate list, checking every entry names the lone input.
    fn planned_paths(planned: &[PlannedOutput], input: &str) -> Vec<String> {
        planned
            .iter()
            .map(|(path, src)| {
                assert_eq!(src.to_string(), input, "every candidate names one input");
                path.display().to_string()
            })
            .collect()
    }

    #[test]
    fn interleaved_fastq_planned_lists_every_field_in_every_option_state() {
        let src = OutputSource::Input(PathBuf::from("ip.bam"));
        let paths = |retain, reports| {
            let o = InterleavedFastqOutputs::new(Path::new("ip.bam"), None, false, retain, reports);
            planned_paths(&o.planned(&src), "ip.bam")
        };

        assert_eq!(paths(false, false), ["ip_val_1.fq", "ip_val_2.fq"]);
        assert_eq!(
            paths(true, false),
            [
                "ip_val_1.fq",
                "ip_val_2.fq",
                "ip_unpaired_1.fq",
                "ip_unpaired_2.fq"
            ]
        );
        assert_eq!(
            paths(false, true),
            [
                "ip_val_1.fq",
                "ip_val_2.fq",
                "ip_R1_trimming_report.txt",
                "ip_R1_trimming_report.json",
                "ip_R2_trimming_report.txt",
                "ip_R2_trimming_report.json"
            ]
        );
        assert_eq!(
            paths(true, true),
            [
                "ip_val_1.fq",
                "ip_val_2.fq",
                "ip_unpaired_1.fq",
                "ip_unpaired_2.fq",
                "ip_R1_trimming_report.txt",
                "ip_R1_trimming_report.json",
                "ip_R2_trimming_report.txt",
                "ip_R2_trimming_report.json"
            ]
        );
    }

    #[test]
    fn interleaved_bam_planned_lists_every_field_in_every_option_state() {
        let src = OutputSource::Input(PathBuf::from("ip.bam"));
        let paths = |reports| {
            let o = InterleavedBamOutputs::new(Path::new("ip.bam"), None, reports);
            planned_paths(&o.planned(&src), "ip.bam")
        };

        assert_eq!(paths(false), ["ip_val.bam"]);
        assert_eq!(
            paths(true),
            [
                "ip_val.bam",
                "ip_R1_trimming_report.txt",
                "ip_R1_trimming_report.json",
                "ip_R2_trimming_report.txt",
                "ip_R2_trimming_report.json"
            ]
        );
    }

    /// A refusal on these arms names one input. `OutputSource::Pair` would have
    /// only `Pair(x, x)` to offer, which renders `ip.bam + ip.bam`. The duplicate
    /// is synthetic: `planned()` cannot produce one (#424).
    #[test]
    fn interleaved_refusal_names_one_input() {
        let src = OutputSource::Input(PathBuf::from("ip.bam"));
        let one = InterleavedBamOutputs::new(Path::new("ip.bam"), None, false).planned(&src);
        let err = preflight_output_collisions(&[one[0].clone(), one[0].clone()], &[], None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("named from ip.bam"), "got: {err}");
        assert!(!err.contains("ip.bam + ip.bam"), "got: {err}");
    }

    #[test]
    fn partial_output_name_wraps_the_whole_file_name() {
        assert_eq!(
            partial_output_name(Path::new("/d/x_trimmed.fq.gz")),
            PathBuf::from("/d/.x_trimmed.fq.gz.partial")
        );
        assert_eq!(
            partial_output_name(Path::new("x_trimmed.fq")),
            PathBuf::from(".x_trimmed.fq.partial")
        );
    }

    /// `commit` sets `released` before attempting the rename, so a rename that
    /// fails leaves the trimmed data at the temporary instead of `Drop`
    /// removing it. Provoked with a non-empty directory at the final path.
    #[test]
    fn a_failed_commit_leaves_the_temporary_on_disk() {
        let dir = std::env::temp_dir().join("tg_428_commit_fail");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let final_path = dir.join("out.fq");
        let (pending, mut file) = PendingOutput::create(&final_path).unwrap();
        let tmp = pending.tmp_path().to_path_buf();
        std::io::Write::write_all(&mut file, b"@r\nA\n+\nI\n").unwrap();
        drop(file);

        std::fs::create_dir(&final_path).unwrap();
        std::fs::write(final_path.join("blocker"), b"x").unwrap();

        let err = pending
            .commit()
            .expect_err("renaming onto a directory must fail");
        assert!(
            tmp.exists(),
            "a failed commit must leave the trimmed data at the temporary"
        );
        let msg = format!("{err:#}");
        assert!(
            msg.contains("out.fq") && msg.contains(".partial"),
            "the error must name both paths, got: {msg}"
        );
    }

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

    /// REGRESSION ([#381](https://github.com/FelixKrueger/TrimGalore/issues/381)).
    /// `bgzip` output is commonly named `.bgz`/`.bgzf`, and since #374 such a
    /// file is read correctly. The stem has to follow, or `sample.fq.bgz`
    /// produces `sample.fq_trimmed.fq` next to
    /// `sample.fq.bgz_trimming_report.txt` and the two disagree about the
    /// sample name (which is what MultiQC groups on).
    #[test]
    fn test_strip_fastq_extensions_bgz() {
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fastq.bgz")),
            "sample"
        );
        assert_eq!(strip_fastq_extensions(Path::new("sample.fq.bgz")), "sample");
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fastq.bgzf")),
            "sample"
        );
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fq.bgzf")),
            "sample"
        );
        assert_eq!(
            strip_fastq_extensions(Path::new("sample_R1.fq.bgz")),
            "sample_R1"
        );
        // Bare compression suffix, no inner .fastq/.fq. The old fallback
        // already handled this one via `file_stem`; pinned so the rewrite
        // below does not lose it.
        assert_eq!(strip_fastq_extensions(Path::new("sample.bgz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.bgzf")), "sample");
        // Multi-component name with no inner .fastq/.fq: the fallback now runs
        // on the already-stripped name, so one more component goes than before
        // (this was `sample.txt`). That is what `.txt.gz` has always done.
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.txt.bgz")),
            "sample"
        );
    }

    /// The non-FASTQ names that reach this function must keep their existing
    /// stems: `.bam` inputs (uBAM), bare `.gz`, and unrelated extensions all
    /// fall through to `file_stem`.
    #[test]
    fn test_strip_fastq_extensions_leaves_other_names_alone() {
        assert_eq!(strip_fastq_extensions(Path::new("sample.bam")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.gz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.txt.gz")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample")), "sample");
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

    /// #244 regression: `--basename foo --paired` must produce
    /// `foo_val_1.fq.gz` / `foo_val_2.fq.gz`, NOT `foo_R1_val_1.fq.gz`
    /// / `foo_R2_val_2.fq.gz`. Beta.5 silently interpolated `_R1`/`_R2`
    /// between the basename and the `_val_{1,2}` suffix, which broke
    /// every nf-core / Snakemake pipeline globbing the documented
    /// Perl path. Matches Perl v0.6.5+ semantics.
    #[test]
    fn test_paired_end_output_names_with_basename() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        // No --output_dir: falls back to input parent dir.
        let (o1, o2) = paired_end_output_names(r1, r2, None, Some("foo"), true);
        assert_eq!(o1, PathBuf::from("/data/foo_val_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/data/foo_val_2.fq.gz"));

        // With --output_dir
        let out = Path::new("/tmp/out");
        let (o1, o2) = paired_end_output_names(r1, r2, Some(out), Some("foo"), true);
        assert_eq!(o1, PathBuf::from("/tmp/out/foo_val_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/tmp/out/foo_val_2.fq.gz"));

        // gzip=false suffix
        let (o1, o2) = paired_end_output_names(r1, r2, None, Some("foo"), false);
        assert_eq!(o1, PathBuf::from("/data/foo_val_1.fq"));
        assert_eq!(o2, PathBuf::from("/data/foo_val_2.fq"));
    }

    /// #244 regression for the `--retain_unpaired` companion path —
    /// same `_R1`/`_R2` interpolation bug, same fix. With basename "foo"
    /// the unpaired files must be `foo_unpaired_1.fq.gz` /
    /// `foo_unpaired_2.fq.gz`.
    #[test]
    fn test_unpaired_output_names_with_basename() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let (o1, o2) = unpaired_output_names(r1, r2, None, Some("foo"), true);
        assert_eq!(o1, PathBuf::from("/data/foo_unpaired_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/data/foo_unpaired_2.fq.gz"));
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
    fn test_norm_path_case_folds() {
        // The fold that `collision_key` is built on; production callers only reach it
        // through that key, which absolutises first — hence the direct test here.
        assert_eq!(norm_path(Path::new("foo.fq.gz")), "foo.fq.gz");
        assert_eq!(norm_path(Path::new("FOO.FQ.GZ")), "foo.fq.gz");
        assert_eq!(norm_path(Path::new("Foo.Fq.Gz")), "foo.fq.gz");
        assert_eq!(
            norm_path(Path::new("/Some/Dir/Sample_R1.fq.gz")),
            "/some/dir/sample_r1.fq.gz"
        );
        // Case-folded aliases compare equal under norm_path.
        assert_eq!(
            norm_path(Path::new("R1.fq.gz")),
            norm_path(Path::new("r1.fq.gz"))
        );
    }

    // ── passthrough_output_name (plan v2 Step 2) ──────────────────────────

    /// #398 — stem from the carrier, directory from R1. R1 is deliberately in a
    /// different directory here; a carrier-anchored implementation passes only when
    /// the two happen to coincide.
    #[test]
    fn test_passthrough_output_name_takes_r1s_directory() {
        let r1 = Path::new("/reads/R1.fq.gz");
        let pt = Path::new("/data/I1.fq.gz");
        let out = passthrough_output_name(r1, pt, None, None, true);
        assert_eq!(out, PathBuf::from("/reads/I1_passthrough.fq.gz"));
    }

    #[test]
    fn test_passthrough_output_name_plain() {
        let r1 = Path::new("/reads/R1.fq.gz");
        let pt = Path::new("/data/I1.fq.gz");
        let out = passthrough_output_name(r1, pt, None, None, false);
        assert_eq!(out, PathBuf::from("/reads/I1_passthrough.fq"));
    }

    #[test]
    fn test_passthrough_output_name_with_output_dir() {
        let r1 = Path::new("/reads/R1.fq.gz");
        let pt = Path::new("/in/I1.fq.gz");
        let out = passthrough_output_name(r1, pt, Some(Path::new("/out")), None, true);
        assert_eq!(out, PathBuf::from("/out/I1_passthrough.fq.gz"));
    }

    #[test]
    fn test_passthrough_output_name_with_basename() {
        let r1 = Path::new("/reads/R1.fq.gz");
        let pt = Path::new("/data/I1.fq.gz");
        let out = passthrough_output_name(r1, pt, None, Some("foo"), true);
        assert_eq!(out, PathBuf::from("/reads/foo_passthrough.fq.gz"));
    }

    #[test]
    fn test_passthrough_output_name_basename_with_output_dir() {
        let r1 = Path::new("/reads/R1.fq.gz");
        let pt = Path::new("/in/I1.fq.gz");
        let out = passthrough_output_name(r1, pt, Some(Path::new("/out")), Some("foo"), true);
        assert_eq!(out, PathBuf::from("/out/foo_passthrough.fq.gz"));
    }

    #[test]
    fn test_passthrough_output_name_all_extensions() {
        // Every FASTQ extension Trim Galore recognises strips cleanly.
        let r1 = Path::new("/d/R1.fq.gz");
        for (input, expected) in [
            ("/d/s.fastq.gz", "/d/s_passthrough.fq.gz"),
            ("/d/s.fq.gz", "/d/s_passthrough.fq.gz"),
            ("/d/s.fastq", "/d/s_passthrough.fq.gz"),
            ("/d/s.fq", "/d/s_passthrough.fq.gz"),
        ] {
            let out = passthrough_output_name(r1, Path::new(input), None, None, true);
            assert_eq!(out, PathBuf::from(expected), "input={input}");
        }
    }

    // ── pair_output_dir (#398) ────────────────────────────────────────────

    #[test]
    fn test_pair_output_dir_precedence_and_edges() {
        // --output_dir wins over R1's parent.
        assert_eq!(
            pair_output_dir(Path::new("/reads/R1.fq"), Some(Path::new("/out"))),
            PathBuf::from("/out")
        );
        // Otherwise R1's parent.
        assert_eq!(
            pair_output_dir(Path::new("/reads/R1.fq"), None),
            PathBuf::from("/reads")
        );
        // A bare filename has an EMPTY parent, not `.`, so joined paths stay bare —
        // returning `.` here would prefix every path and message with `./`.
        assert_eq!(pair_output_dir(Path::new("R1.fq"), None), PathBuf::from(""));
        assert_eq!(
            pair_output_dir(Path::new("R1.fq"), None).join("R1.fq_trimming_report.txt"),
            PathBuf::from("R1.fq_trimming_report.txt")
        );
        // The `.` fallback is reachable only for a root path, which has no parent.
        assert_eq!(pair_output_dir(Path::new("/"), None), PathBuf::from("."));
    }

    // ── dump_planned_lines (#414) ─────────────────────────────────────────

    /// The pure renderer, tested without touching the environment: `env::set_var` is
    /// `unsafe` under edition 2024 and would leak into sibling threads under
    /// `cargo test --release`. The env gate is driven end-to-end from
    /// `tests/integration_preflight_tripwire.rs` instead.
    #[test]
    fn dump_planned_lines_renders_one_line_per_candidate() {
        let planned: Vec<PlannedOutput> = vec![
            (
                PathBuf::from("out/a_trimmed.fq"),
                OutputSource::Input(PathBuf::from("a.fastq")),
            ),
            (
                PathBuf::from("out/a.fastq_trimming_report.txt"),
                OutputSource::Input(PathBuf::from("a.fastq")),
            ),
        ];
        let lines = dump_planned_lines(&planned, "src/main.rs:1006");

        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            "TRIM_GALORE_PLANNED\tsrc/main.rs:1006\tout/a_trimmed.fq"
        );
        // The site is repeated on every line so the parse needs no state.
        for line in &lines {
            assert!(line.starts_with(DUMP_PLANNED_PREFIX));
            assert_eq!(line.matches('\t').count(), 2);
        }
        // An empty list dumps nothing, so a run that plans nothing is distinguishable
        // from one that never reached a pre-flight.
        assert!(dump_planned_lines(&[], "src/main.rs:1").is_empty());
    }

    // ── clumped_output_name / clumped_paired_output_names / clumping_report_name ──

    #[test]
    fn test_clumped_output_name_gzip() {
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_output_name(input, None, None, true);
        assert_eq!(out, PathBuf::from("/data/sample_clumped.fq.gz"));
    }

    #[test]
    fn test_clumped_output_name_plain() {
        // --dont_gzip is allowed under --clump_only (diverges from --clumpify).
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_output_name(input, None, None, false);
        assert_eq!(out, PathBuf::from("/data/sample_clumped.fq"));
    }

    #[test]
    fn test_clumped_output_name_with_basename() {
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_output_name(input, None, Some("archive"), true);
        assert_eq!(out, PathBuf::from("/data/archive_clumped.fq.gz"));
    }

    #[test]
    fn test_clumped_output_name_with_output_dir() {
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_output_name(input, Some(Path::new("/tmp/out")), None, true);
        assert_eq!(out, PathBuf::from("/tmp/out/sample_clumped.fq.gz"));
    }

    #[test]
    fn test_clumped_paired_output_names_bare() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let (o1, o2) = clumped_paired_output_names(r1, r2, None, None, true);
        assert_eq!(o1, PathBuf::from("/data/sample_R1_clumped_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/data/sample_R2_clumped_2.fq.gz"));
    }

    #[test]
    fn test_clumped_paired_output_names_with_basename() {
        // Same #244-shaped semantic: `--basename foo` yields `foo_clumped_{1,2}`
        // with no `_R1`/`_R2` interpolation.
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let (o1, o2) = clumped_paired_output_names(r1, r2, None, Some("foo"), true);
        assert_eq!(o1, PathBuf::from("/data/foo_clumped_1.fq.gz"));
        assert_eq!(o2, PathBuf::from("/data/foo_clumped_2.fq.gz"));
    }

    #[test]
    fn test_clumped_paired_output_names_plain() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let (o1, o2) = clumped_paired_output_names(r1, r2, None, None, false);
        assert_eq!(o1, PathBuf::from("/data/sample_R1_clumped_1.fq"));
        assert_eq!(o2, PathBuf::from("/data/sample_R2_clumped_2.fq"));
    }

    // ── BAM-output filename helpers (v2) ──────────────────────────────

    #[test]
    fn test_clumped_bam_output_name_from_fastq_input() {
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_bam_output_name(input, None, None);
        assert_eq!(out, PathBuf::from("/data/sample_clumped.bam"));
    }

    #[test]
    fn test_clumped_bam_output_name_from_bam_input() {
        let input = Path::new("/data/sample.bam");
        let out = clumped_bam_output_name(input, None, None);
        // strip_fastq_extensions falls through to file_stem() for .bam.
        assert_eq!(out, PathBuf::from("/data/sample_clumped.bam"));
    }

    #[test]
    fn test_clumped_bam_output_name_with_basename() {
        let input = Path::new("/data/sample.bam");
        let out = clumped_bam_output_name(input, None, Some("archive"));
        assert_eq!(out, PathBuf::from("/data/archive_clumped.bam"));
    }

    #[test]
    fn test_clumped_bam_output_name_with_output_dir() {
        let input = Path::new("/data/sample.fq.gz");
        let out = clumped_bam_output_name(input, Some(Path::new("/tmp/out")), None);
        assert_eq!(out, PathBuf::from("/tmp/out/sample_clumped.bam"));
    }

    #[test]
    fn test_clumped_paired_bam_output_name_single_file_per_pair() {
        // PE-BAM produces ONE interleaved file per pair (samtools convention).
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let out = clumped_paired_bam_output_name(r1, Some(r2), None, None);
        assert_eq!(out, PathBuf::from("/data/sample_R1_clumped.bam"));
    }

    #[test]
    fn test_clumped_paired_bam_output_name_from_interleaved_bam() {
        // Shape B: single interleaved uBAM input.
        let bam = Path::new("/data/interleaved.bam");
        let out = clumped_paired_bam_output_name(bam, None, None, None);
        assert_eq!(out, PathBuf::from("/data/interleaved_clumped.bam"));
    }

    #[test]
    fn test_clumped_paired_bam_output_name_with_basename() {
        let r1 = Path::new("/data/sample_R1.fq.gz");
        let r2 = Path::new("/data/sample_R2.fq.gz");
        let out = clumped_paired_bam_output_name(r1, Some(r2), None, Some("foo"));
        assert_eq!(out, PathBuf::from("/data/foo_clumped.bam"));
    }

    #[test]
    fn test_clumping_report_name_distinct_from_trimming_report() {
        // Load-bearing invariant: the clump-only report filename must NOT
        // match the `*_trimming_report.*` glob that nf-core / MultiQC scan.
        let input = Path::new("/data/sample.fq.gz");
        let out = clumping_report_name(input, None);
        assert_eq!(out, PathBuf::from("/data/sample.fq.gz_clumping_report.txt"));
        // Confirm it's NOT the trimming-report shape.
        assert_ne!(out, report_name(input, None));
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
        // Heuristic is extension-based on purpose (the reader sniffs content;
        // this names the output), so a misnamed file doesn't trigger.
        assert!(!is_gzipped(Path::new("sample.gz.fastq")));
    }

    /// #384. Case variants of `.gz` must fold; nothing else may start matching.
    /// This is the committed guard for the compression half of the change.
    #[test]
    fn test_is_gzipped_folds_case() {
        assert!(is_gzipped(Path::new("SAMPLE.FASTQ.GZ")));
        assert!(is_gzipped(Path::new("sample.fastq.Gz")));
        assert!(is_gzipped(Path::new("SAMPLE.FQ.gz")));
        // .bgz is not .gz in ANY case — the pre-existing bgz asymmetry stays,
        // and the fold must not widen is_gzipped to the stem's gzip-family list.
        assert!(!is_gzipped(Path::new("sample.fastq.bgz")));
        assert!(!is_gzipped(Path::new("SAMPLE.FASTQ.BGZ")));
        assert!(!is_gzipped(Path::new("SAMPLE.FASTQ")));
        assert!(!is_gzipped(Path::new("sample.GZ.fastq")));
    }

    /// #384. Extension matching folds; the retained part of the name must not.
    #[test]
    fn test_strip_fastq_extensions_folds_case() {
        for (input, stem) in [
            ("SAMPLE.FASTQ.GZ", "SAMPLE"),
            ("Sample.FastQ.Gz", "Sample"),
            ("SAMPLE.FQ", "SAMPLE"),
            ("SAMPLE.FASTQ", "SAMPLE"),
            ("SAMPLE.FQ.GZ", "SAMPLE"),
            ("SAMPLE.FASTQ.BGZ", "SAMPLE"),
            ("SAMPLE.FASTQ.BGZF", "SAMPLE"),
            ("sample_R1.FQ.GZ", "sample_R1"),
        ] {
            assert_eq!(
                strip_fastq_extensions(Path::new(input)),
                stem,
                "input {input}"
            );
        }
    }

    /// #384 regression. The naive slice form of the case-insensitive strip
    /// panics mid-character on multi-byte names; these inputs work end-to-end
    /// today and must keep doing so.
    #[test]
    fn test_strip_fastq_extensions_non_ascii_and_short_names() {
        assert_eq!(strip_fastq_extensions(Path::new("😀")), "😀");
        assert_eq!(strip_fastq_extensions(Path::new("é.fq")), "é");
        assert_eq!(strip_fastq_extensions(Path::new("a")), "a");
        // Non-FASTQ extensions keep the file_stem fallback, un-folded.
        assert_eq!(strip_fastq_extensions(Path::new("sample.bam")), "sample");
        assert_eq!(strip_fastq_extensions(Path::new("sample.txt")), "sample");
    }

    // --- preflight_output_collisions (issues #216, #383) ---
    //
    // Lexical throughout, so none of these touch the filesystem.

    fn pb(s: &str) -> PathBuf {
        PathBuf::from(s)
    }

    /// A planned output attributed to its own synthetic input (#397). Distinct
    /// sources on purpose: these cases are about two *different* inputs racing for
    /// one path, which is shape 2a/2b — `po_same` covers the one-input shape 2c.
    fn po(s: &str) -> PlannedOutput {
        (pb(s), OutputSource::Input(pb(&format!("src_of_{s}"))))
    }

    /// Two planned outputs from the *same* input — shape 2c.
    fn po_from(src: &str, s: &str) -> PlannedOutput {
        (pb(s), OutputSource::Input(pb(src)))
    }

    #[test]
    fn preflight_accepts_distinct_outputs() {
        let planned = [po("a_trimmed.fq.gz"), po("b_trimmed.fq.gz")];
        let inputs = [pb("a.fastq.gz"), pb("b.fastq.gz")];
        assert!(preflight_output_collisions(&planned, &inputs, None).is_ok());
    }

    #[test]
    fn preflight_accepts_empty_and_single() {
        assert!(preflight_output_collisions(&[], &[], None).is_ok());
        assert!(
            preflight_output_collisions(&[po("a_trimmed.fq.gz")], &[pb("a.fastq.gz")], None)
                .is_ok()
        );
    }

    #[test]
    fn preflight_rejects_identical_outputs() {
        let planned = [po("x_trimmed.fq.gz"), po("x_trimmed.fq.gz")];
        let err = preflight_output_collisions(&planned, &[], None)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("would be written to the same file"),
            "got: {err}"
        );
        assert!(err.contains("x_trimmed.fq.gz"), "must name the path: {err}");
    }

    /// Issue #216. Unreachable from an integration test: two paths differing only
    /// in case cannot coexist on APFS.
    #[test]
    fn preflight_rejects_case_only_variants() {
        let planned = [po("Sample_trimmed.fq.gz"), po("SAMPLE_trimmed.fq.gz")];
        let err = preflight_output_collisions(&planned, &[], None)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("would be written to the same file"),
            "got: {err}"
        );
    }

    /// REGRESSION (#383). A raw-string key let `./x` and `x` through while naming
    /// one file, so the reported bug survived its own fix.
    #[test]
    fn preflight_rejects_dot_slash_alias() {
        let planned = [po("./x_trimmed.fq.gz"), po("x_trimmed.fq.gz")];
        let err = preflight_output_collisions(&planned, &[], None)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("would be written to the same file"),
            "got: {err}"
        );
    }

    /// REGRESSION (#383). Same defect via a mixed absolute/relative argument list.
    #[test]
    fn preflight_rejects_absolute_versus_relative_alias() {
        let abs = std::env::current_dir().unwrap().join("x_trimmed.fq.gz");
        let planned = [
            (abs, OutputSource::Input(pb("src_abs"))),
            po("x_trimmed.fq.gz"),
        ];
        let err = preflight_output_collisions(&planned, &[], None)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("would be written to the same file"),
            "got: {err}"
        );
    }

    /// #397 shape 2c — two candidates from ONE input. "Rename one input" cannot be
    /// followed here, so this branch supplies its own advice and ignores the
    /// call-site hint, exactly as the output-vs-input branch already does.
    #[test]
    fn preflight_same_source_collision_replaces_the_advice() {
        let planned = [
            po_from("d/x.fq", "d/x.fq_report.txt"),
            po_from("d/x.fq", "d/x.fq_report.txt"),
        ];
        let err = preflight_output_collisions(&planned, &[], Some("EXTRA-HINT."))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("would be written to the same file"),
            "got: {err}"
        );
        assert!(err.contains("List each input once"), "got: {err}");
        assert!(
            !err.contains("EXTRA-HINT.") && !err.contains("different source directories"),
            "2c supplies its own advice and must not carry the hint: {err}"
        );
        // The one source is named, and the path appears once, not twice.
        assert!(err.contains("d/x.fq"), "got: {err}");
    }

    /// REGRESSION (#383). An output that would overwrite an input gets its own
    /// message: "same file" is untrue here and its remedies do not apply.
    #[test]
    fn preflight_rejects_output_that_aliases_an_input() {
        let planned = [po("s_trimmed.fq.gz")];
        let inputs = [pb("s.fastq.gz"), pb("s_trimmed.fq.gz")];
        let err = preflight_output_collisions(&planned, &inputs, None)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("which is also one of its inputs"),
            "got: {err}"
        );
        assert!(
            !err.contains("would be written to the same file"),
            "must not fall back to the duplicate-output wording: {err}"
        );
    }

    /// The alias check is keyed the same way, so a `./` spelling still catches it.
    #[test]
    fn preflight_rejects_input_alias_across_spellings() {
        let planned = [po("./s_trimmed.fq.gz")];
        let inputs = [pb("s_trimmed.fq.gz")];
        assert!(preflight_output_collisions(&planned, &inputs, None).is_err());
    }

    /// REGRESSION (#383). `..` used to defeat the key, so `trim_galore
    /// ../data/s.fastq /abs/data/s.fq` lost one input's reads at exit 0 — the
    /// reported bug reproducing through its own fix.
    #[test]
    fn preflight_rejects_dotdot_alias() {
        let planned = [po("a/../x_trimmed.fq.gz"), po("x_trimmed.fq.gz")];
        assert!(preflight_output_collisions(&planned, &[], None).is_err());
    }

    /// `..` that cannot be folded away must not be silently dropped: two paths
    /// differing only below a kept `..` stay distinct.
    #[test]
    fn preflight_keeps_unfoldable_paths_distinct() {
        let planned = [po("../x_trimmed.fq.gz"), po("../y_trimmed.fq.gz")];
        assert!(preflight_output_collisions(&planned, &[], None).is_ok());
    }

    #[test]
    fn preflight_appends_hint_only_when_given() {
        // Two *different* inputs racing for one path — shape 2a, which is where the
        // hint applies. Shape 2c (one input, two candidates) deliberately carries
        // neither hint nor generic advice; pinned separately below.
        let planned = [
            po_from("srcA.fq", "x_trimmed.fq.gz"),
            po_from("srcB.fq", "x_trimmed.fq.gz"),
        ];
        let with = preflight_output_collisions(&planned, &[], Some("EXTRA-HINT."))
            .unwrap_err()
            .to_string();
        assert!(with.contains("EXTRA-HINT."), "got: {with}");
        let without = preflight_output_collisions(&planned, &[], None)
            .unwrap_err()
            .to_string();
        assert!(!without.contains("EXTRA-HINT."), "got: {without}");
        // ...and the generic advice appears only when no hint was given.
        assert!(
            without.contains("different source directories"),
            "got: {without}"
        );
        assert!(
            !with.contains("different source directories"),
            "a hint must REPLACE the generic advice, not append to it: {with}"
        );
    }

    /// The two keys differ in exactly one respect, and that difference is load-bearing.
    ///
    /// REGRESSION: re-keying `Cli::validate`'s input-identity checks on the case-folded
    /// key made the #216 CI guard fail — it feeds four genuinely distinct files
    /// (`Sample_R1` / `SAMPLE_R1`) on a case-sensitive filesystem and asserts the
    /// *output* pre-flight refuses them. Case-folding input identity called them
    /// duplicates instead, which on Linux is false.
    #[test]
    fn identity_key_is_case_sensitive_but_collision_key_is_not() {
        let lower = Path::new("Sample_R1.fastq.gz");
        let upper = Path::new("SAMPLE_R1.fastq.gz");
        assert_ne!(
            path_identity_key(lower),
            path_identity_key(upper),
            "two files differing only in case are distinct files"
        );
        assert_eq!(
            collision_key(lower),
            collision_key(upper),
            "but their outputs alias each other on APFS/NTFS"
        );
    }

    /// Both keys normalise spelling — that is the part #383 needed.
    #[test]
    fn both_keys_normalise_spelling() {
        for a in ["./x.fq", "a/../x.fq"] {
            assert_eq!(
                path_identity_key(Path::new(a)),
                path_identity_key(Path::new("x.fq"))
            );
            assert_eq!(
                collision_key(Path::new(a)),
                collision_key(Path::new("x.fq"))
            );
        }
    }

    /// Assumption A2, in the direction that makes the pre-flight sufficient:
    /// where two inputs' PRIMARY output paths differ, every secondary output
    /// path must differ too — so a secondary can never collide unless a primary
    /// already has, and hashing primaries alone is enough. Covers the SE trim
    /// and clump primary namers, compared on `collision_key` (the pre-flight's
    /// own metric). Paired `_val_N` / `_clumped_N` primaries break the
    /// implication (#388, #391), which is why the paired pre-flights carry
    /// report candidates explicitly.
    ///
    /// Checked across the flag matrix (`--basename` / `--dont_gzip` / `-o`, each
    /// on and off) because `--basename` and `-o` are exactly the flags that
    /// collapse distinct inputs onto one path.
    ///
    /// FastQC's `<stem>_fastqc.zip` is not asserted separately: the bundled
    /// crate derives it from the primary path we hand it, so there is no second
    /// key of ours to compare, and writing the formula out here would only test
    /// the formula. Its one real exception — `--fastqc_args "-o DIR"`, which
    /// overrides the output directory where the pre-flight cannot see it — is
    /// recorded in the plan as a known residual, not covered here.
    #[test]
    fn distinct_primary_outputs_imply_distinct_secondary_outputs() {
        // The last two share a basename across directories. Without that pair the
        // assertion is a tautology: every namer embeds `file_name()`, so distinct
        // basenames make it true regardless of whether a namer honours the directory.
        let inputs = [
            Path::new("d/alpha.fastq.gz"),
            Path::new("d/beta.fastq.gz"),
            Path::new("e/gamma.fq.gz"),
            Path::new("d/same.fastq.gz"),
            Path::new("e/same.fastq.gz"),
            // Fold-equal twin of d/same: the skip must use collision_key, not
            // PathBuf equality, or the assertion below fails on this pair (#391).
            Path::new("d/SAME.fastq.gz"),
        ];
        let out = PathBuf::from("shared_out");

        for basename in [None, Some("fixed")] {
            for gzip in [true, false] {
                for output_dir in [None, Some(out.as_path())] {
                    // #391 — same implication for every clump primary namer, on the
                    // pre-flight's own metric (collision_key) so the fold is real.
                    let primary_sets: Vec<Vec<PathBuf>> = vec![
                        inputs
                            .iter()
                            .map(|p| single_end_output_name(p, output_dir, basename, gzip))
                            .collect(),
                        inputs
                            .iter()
                            .map(|p| clumped_output_name(p, output_dir, basename, gzip))
                            .collect(),
                        inputs
                            .iter()
                            .map(|p| clumped_bam_output_name(p, output_dir, basename))
                            .collect(),
                        inputs
                            .iter()
                            .map(|p| clumped_paired_bam_output_name(p, None, output_dir, basename))
                            .collect(),
                    ];
                    for primaries in &primary_sets {
                        for (i, a) in primaries.iter().enumerate() {
                            for (j, b) in primaries.iter().enumerate().take(i) {
                                if collision_key(a) == collision_key(b) {
                                    // --basename collapses every input onto one primary; the
                                    // pre-flight rejects that, and cli.rs:624 rejects it earlier
                                    // still for multi-input SE. Nothing to prove here.
                                    continue;
                                }
                                // Primaries differ, so every secondary must differ too.
                                for namer in [report_name, json_report_name, clumping_report_name] {
                                    assert_ne!(
                                        collision_key(&namer(inputs[i], output_dir)),
                                        collision_key(&namer(inputs[j], output_dir)),
                                        "secondary collided while primaries {a:?} / {b:?} differ \
                                         (basename={basename:?} gzip={gzip} out={output_dir:?})"
                                    );
                                }
                            }
                        }
                    }

                    // The demux stem, via the same function `demultiplex` uses — a
                    // property of the SE trim primary only. Only meaningful when the
                    // two primaries share a directory: demux resolves its output dir
                    // to `-o` else the primary's parent, so differing parents already
                    // separate the paths and the stem is allowed to repeat.
                    let se_primaries = &primary_sets[0];
                    for (i, a) in se_primaries.iter().enumerate() {
                        for b in se_primaries.iter().take(i) {
                            if collision_key(a) != collision_key(b) && a.parent() == b.parent() {
                                assert_ne!(
                                    crate::demux::demux_base_name(a),
                                    crate::demux::demux_base_name(b),
                                    "demux stems collided in one directory while primaries \
                                     differ: {a:?} / {b:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /// Complement to the above: the primary key is strictly *coarser* than the
    /// report key (in single-end naming — paired `_val_N` inverts this, #388,
    /// as does `clumped_paired_output_names`' `_clumped_N`, #391; the clump
    /// primaries asserted below all collapse), which is why checking SE
    /// primaries covers reports rather than merely coinciding with them.
    /// Three spellings of one sample share a primary while keeping three
    /// distinct report names.
    #[test]
    fn primary_output_key_is_coarser_than_secondary_keys() {
        let variants = [
            Path::new("d/sample.fastq.gz"),
            Path::new("d/sample.fq.gz"),
            Path::new("d/sample.fastq.bgz"),
        ];

        // All three collapse to one primary output...
        let primaries: Vec<PathBuf> = variants
            .iter()
            .map(|p| single_end_output_name(p, None, None, true))
            .collect();
        assert!(
            primaries.windows(2).all(|w| w[0] == w[1]),
            "expected one shared primary, got {primaries:?}"
        );

        // ...while every secondary name stays distinct, so a secondary can never
        // collide unless the primary already has.
        for namer in [report_name, json_report_name, clumping_report_name] {
            let secondaries: Vec<PathBuf> = variants.iter().map(|p| namer(p, None)).collect();
            for (i, a) in secondaries.iter().enumerate() {
                for b in &secondaries[..i] {
                    assert_ne!(a, b, "secondary names must be distinct: {a:?} vs {b:?}");
                }
            }
        }

        // Same property for the uBAM primary.
        let bam: Vec<PathBuf> = variants
            .iter()
            .map(|p| single_end_bam_output_name(p, None, None))
            .collect();
        assert!(bam.windows(2).all(|w| w[0] == w[1]), "got {bam:?}");

        // #391 — and for every clump primary namer: the three spellings share
        // one primary each, so a clump report can only collide where the
        // primary already has.
        let clumped: Vec<PathBuf> = variants
            .iter()
            .map(|p| clumped_output_name(p, None, None, true))
            .collect();
        assert!(clumped.windows(2).all(|w| w[0] == w[1]), "got {clumped:?}");
        let clumped_bam: Vec<PathBuf> = variants
            .iter()
            .map(|p| clumped_bam_output_name(p, None, None))
            .collect();
        assert!(
            clumped_bam.windows(2).all(|w| w[0] == w[1]),
            "got {clumped_bam:?}"
        );
        let clumped_pe_bam: Vec<PathBuf> = variants
            .iter()
            .map(|p| clumped_paired_bam_output_name(p, None, None, None))
            .collect();
        assert!(
            clumped_pe_bam.windows(2).all(|w| w[0] == w[1]),
            "got {clumped_pe_bam:?}"
        );
    }
}
